
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
  int    n, k;
  real   tmp, tmpvir;
  vektor tmpvec;

  tot_pot_energy = 0.0;
  virial         = 0.0;
  vir_vect.x     = 0.0;
  vir_vect.y     = 0.0;

  /* Zero Forces */
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
  for (n=0; n<4; ++n) {
    for (k=0; k<npairs[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n] + k;
      pbc.x = P->ipbc[0] * box_x.x + P->ipbc[1] * box_y.x;
      pbc.y = P->ipbc[0] * box_x.y + P->ipbc[1] * box_y.y;
      do_forces(cell_array + P->np, cell_array + P->nq, pbc);
    }
  }

#ifndef AR

  /* If we don't use actio=reactio accross the cpus, we have do do
     the force loop also on the other half of the neighbours for the 
     cells on the surface of the CPU */

  /* potential energy and virial are already complete; to avoid double
     counting, we keep a copy of the current value, which we use later */

  tmp      = tot_pot_energy;
  tmpvir   = virial;
  tmpvec.x = vir_vect.x;
  tmpvec.y = vir_vect.y;

  /* compute forces for remaining pairs of cells */
  for (n=0; n<4; ++n) {
    for (k=npairs[n]; k<npairs2[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n] + k;
      pbc.x = P->ipbc[0] * box_x.x + P->ipbc[1] * box_y.x;
      pbc.y = P->ipbc[0] * box_x.y + P->ipbc[1] * box_y.y;
      do_forces(cell_array + P->np, cell_array + P->nq, pbc);
    }
  }

  /* use the previously saved values of potential energy and virial */
  tot_pot_energy = tmp;
  virial     = tmpvir;
  vir_vect.x = tmpvec.x;
  vir_vect.y = tmpvec.y;

#endif /* AR */

  /* sum up results of different CPUs */
  MPI_Allreduce( &tot_pot_energy, &tmp, 1, REAL, MPI_SUM, cpugrid); 
  tot_pot_energy = tmp; 

  MPI_Allreduce( &virial,   &tmpvir,    1, REAL, MPI_SUM, cpugrid);
  virial = tmpvir;

  MPI_Allreduce( &vir_vect, &tmpvec,  DIM, REAL, MPI_SUM, cpugrid);
  vir_vect.x = tmpvec.x;
  vir_vect.y = tmpvec.y;

}


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
    isend_buf( &send_buf_east, nbeast, &reqwest[0] );
    irecv_buf( &recv_buf_west, nbwest, &reqwest[1] );

    /* send west */
    isend_buf( &send_buf_west, nbwest, &reqeast[0] );
    irecv_buf( &recv_buf_east, nbeast, &reqeast[1] );

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
    isend_buf( &send_buf_north, nbnorth, &reqsouth[0] );
    irecv_buf( &recv_buf_south, nbsouth, &reqsouth[1] );

    /* send atoms south */
    isend_buf( &send_buf_south, nbsouth, &reqnorth[0] );
    irecv_buf( &recv_buf_north, nbnorth, &reqnorth[1] );

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




