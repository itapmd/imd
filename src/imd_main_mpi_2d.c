
/******************************************************************************
*
* imd_main_mpi_2d.c -- main loop, MPI specific part, two dimensions
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"


/******************************************************************************
*
* calc_forces 
*
* The forces of the atoms are calulated here. To achive this, atoms on
* the surface of a cpu are exchanged with the neigbours.
*
* The force calculation is split into those steps:
*
* i)   send atoms positions of cells on surface neighbours, 
*      receive atom positions from neigbours
* ii)  zero forces on all cells (local and buffer)
* iii) calculate forces in local cells, use lower half of neigbours 
*      for each cell and use actio==reactio
* iv)  calculate forces also for upper half of neighbours for all cells
*      that are on the upper surface
*
******************************************************************************/

void calc_forces(void)
{
  int  n, k;
  real tmpvec1[4], tmpvec2[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

  /* clear global accumulation variables */
  tot_pot_energy = 0.0;
  virial         = 0.0;
  vir_x          = 0.0;
  vir_y          = 0.0;

  /* clear per atom accumulation variables */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (k=0; k<nallcells; ++k) {
    int i;
    cell *p;
    p = cell_array + k;
    for (i=0; i<p->n; ++i) {
      p->kraft X(i) = 0.0;
      p->kraft Y(i) = 0.0;
      p->pot_eng[i] = 0.0;
#ifdef ORDPAR
      p->nbanz[i] = 0;
#endif
#ifdef NVX
      p->heatcond[i] = 0.0;
#endif      
#ifdef STRESS_TENS
      p->presstens X(i) = 0.0;
      p->presstens Y(i) = 0.0;
      p->presstens_offdia[i] = 0.0;
#endif      
    }
  }

  /* What follows is the standard one-cpu force 
     loop acting on our local data cells */

  /* compute forces for all pairs of cells */
  for (n=0; n<nlists; ++n) {
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tot_pot_energy,virial,vir_x,vir_y)
#endif
    for (k=0; k<npairs[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n] + k;
      pbc.x = P->ipbc[0] * box_x.x + P->ipbc[1] * box_y.x;
      pbc.y = P->ipbc[0] * box_x.y + P->ipbc[1] * box_y.y;
      do_forces(cell_array + P->np, cell_array + P->nq, pbc,
                &tot_pot_energy, &virial, &vir_x, &vir_y, &vir_z);
    }
  }

#ifndef AR

  /* If we don't use actio=reactio accross the cpus, we have do do
     the force loop also on the other half of the neighbours for the 
     cells on the surface of the CPU */

  /* compute forces for remaining pairs of cells */
  for (n=0; n<nlists; ++n) {
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tot_pot_energy,virial,vir_x,vir_y)
#endif
    for (k=npairs[n]; k<npairs2[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n] + k;
      pbc.x = P->ipbc[0] * box_x.x + P->ipbc[1] * box_y.x;
      pbc.y = P->ipbc[0] * box_x.y + P->ipbc[1] * box_y.y;
      /* potential energy and virial are already complete;          */
      /* to avoid double counting, we update only the dummy tmpvec2 */
      do_forces(cell_array + P->np, cell_array + P->nq, pbc,
                tmpvec2, tmpvec2+1, tmpvec2+2, tmpvec2+3, tmpvec2+4);
    }
  }

#endif /* AR */

  /* sum up results of different CPUs */
  tmpvec1[0] = tot_pot_energy;
  tmpvec1[1] = virial;
  tmpvec1[2] = vir_x;
  tmpvec1[3] = vir_y;

  MPI_Allreduce( tmpvec1, tmpvec2, 4, REAL, MPI_SUM, cpugrid); 

  tot_pot_energy = tmpvec2[0];
  virial         = tmpvec2[1];
  vir_x          = tmpvec2[2];
  vir_y          = tmpvec2[3];

}

#ifdef SR

/******************************************************************************
*
* send_atoms  -  only used for fix_cells
*
******************************************************************************/

void send_atoms()
{
  MPI_Status  stat;

  if (cpu_dim.x > 1) {
    /* send east, receive west, move atoms from west to cells */
    sendrecv_buf( &send_buf_east, nbeast, &recv_buf_west, nbwest, &stat);
    MPI_Get_count( &stat, REAL, &recv_buf_west.n );
    process_buffer( &recv_buf_west, (cell *) NULL );

    /* send west, receive east, move atoms from east to cells */
    sendrecv_buf( &send_buf_west, nbwest, &recv_buf_east, nbeast, &stat );
    MPI_Get_count( &stat, REAL, &recv_buf_east.n );
    process_buffer( &recv_buf_east, (cell *) NULL );

    /* append atoms from east & west to north send buffer */
    copy_atoms_buf( &send_buf_north, &recv_buf_west );
    copy_atoms_buf( &send_buf_north, &recv_buf_east );
    /* check special case cpu_dim.y==2 */ 
    if (nbsouth!=nbnorth) {
      /* append atoms from east & west to south send buffer */
      copy_atoms_buf( &send_buf_south, &recv_buf_east );
      copy_atoms_buf( &send_buf_south, &recv_buf_west );
    }
  }

  if (cpu_dim.y > 1) {
    /* send north, receive south, move atoms from south to cells */
    sendrecv_buf( &send_buf_north, nbnorth, &recv_buf_south, nbsouth, &stat);
    MPI_Get_count( &stat, REAL, &recv_buf_south.n );
    process_buffer( &recv_buf_south, (cell *) NULL );   

    /* send south, receive north, move atoms from north to cells */
    sendrecv_buf( &send_buf_south, nbsouth, &recv_buf_north, nbnorth, &stat);
    MPI_Get_count( &stat, REAL, &recv_buf_north.n );
    process_buffer( &recv_buf_north, (cell *) NULL );   
  }

}

#else /* not SR */

/******************************************************************************
*
* send_atoms  -  only used for fix_cells
*
******************************************************************************/

void send_atoms()
{
  MPI_Status  stateast[2], statwest[2], statnorth[2], statsouth[2];
  MPI_Request  reqeast[2],  reqwest[2],  reqnorth[2],  reqsouth[2];

  if (cpu_dim.x > 1) {
    /* send east */
    irecv_buf( &recv_buf_west, nbwest, &reqwest[1] );
    isend_buf( &send_buf_east, nbeast, &reqwest[0] );

    /* send west */
    irecv_buf( &recv_buf_east, nbeast, &reqeast[1] );
    isend_buf( &send_buf_west, nbwest, &reqeast[0] );

    /* wait for atoms from west, move them to cells */
    MPI_Waitall(2, reqwest, statwest);
    MPI_Get_count( &statwest[1], REAL, &recv_buf_west.n );
    process_buffer( &recv_buf_west, (cell *) NULL );

    /* wait for atoms from east, move them to cells */
    MPI_Waitall(2, reqeast, stateast);
    MPI_Get_count( &stateast[1], REAL, &recv_buf_east.n );
    process_buffer( &recv_buf_east, (cell *) NULL );

    /* append atoms from east & west to north send buffer */
    copy_atoms_buf( &send_buf_north, &recv_buf_west );
    copy_atoms_buf( &send_buf_north, &recv_buf_east );
    /* check special case cpu_dim.y==2 */ 
    if (nbsouth!=nbnorth) {
      /* append atoms from east & west to south send buffer */
      copy_atoms_buf( &send_buf_south, &recv_buf_east );
      copy_atoms_buf( &send_buf_south, &recv_buf_west );
    }
  }

  if (cpu_dim.y > 1) {
    /* send atoms north */
    irecv_buf( &recv_buf_south, nbsouth, &reqsouth[1] );
    isend_buf( &send_buf_north, nbnorth, &reqsouth[0] );

    /* send atoms south */
    irecv_buf( &recv_buf_north, nbnorth, &reqnorth[1] );
    isend_buf( &send_buf_south, nbsouth, &reqnorth[0] );

    /* wait for atoms from south, move them to cells */
    MPI_Waitall(2, reqsouth, statsouth);
    MPI_Get_count( &statsouth[1], REAL, &recv_buf_south.n );
    process_buffer( &recv_buf_south, (cell *) NULL );   

    /* Wait for atoms from north, move them to cells */
    MPI_Waitall(2, reqnorth, statnorth);
    MPI_Get_count( &statnorth[1], REAL, &recv_buf_north.n );
    process_buffer( &recv_buf_north, (cell *) NULL );   
  }

}

#endif



