
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2011 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

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

void calc_forces(int steps)
{
  int  n, k;
  real tmpvec1[5], tmpvec2[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  /* fill the buffer cells */
  if ((steps == steps_min) || (0 == steps % BUFSTEP)) setup_buffers();
  send_cells(copy_cell,pack_cell,unpack_cell);

  /* clear global accumulation variables */
  tot_pot_energy = 0.0;
  virial = 0.0;
  vir_xx = 0.0;
  vir_yy = 0.0;
  vir_xy = 0.0;
  nfc++;

  /* clear per atom accumulation variables */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (k=0; k<nallcells; ++k) {
    int i;
    cell *p;
    p = cell_array + k;
    for (i=0; i<p->n; ++i) {
      KRAFT(p,i,X) = 0.0;
      KRAFT(p,i,Y) = 0.0;
      POTENG(p,i)  = 0.0;
#ifdef NNBR
      NBANZ(p,i) = 0;
#endif
#if defined(STRESS_TENS)
      PRESSTENS(p,i,xx) = 0.0;
      PRESSTENS(p,i,yy) = 0.0;
      PRESSTENS(p,i,xy) = 0.0;
#endif      
    }
  }

#ifdef RIGID
  /* clear total forces */
  if ( nsuperatoms>0 ) 
    for(k=0; k<nsuperatoms; k++) {
      superforce[k].x = 0.0;
      superforce[k].y = 0.0;
    }
#endif

  /* What follows is the standard one-cpu force 
     loop acting on our local data cells */

  /* compute forces for all pairs of cells */
  for (n=0; n<nlists; ++n) {
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime) \
  reduction(+:tot_pot_energy,virial,vir_xx,vir_yy,vir_xy)
#endif
    for (k=0; k<npairs[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n] + k;
      pbc.x = P->ipbc[0] * box_x.x + P->ipbc[1] * box_y.x;
      pbc.y = P->ipbc[0] * box_x.y + P->ipbc[1] * box_y.y;
      do_forces(cell_array + P->np, cell_array + P->nq, pbc,
                &tot_pot_energy, &virial, &vir_xx, &vir_yy, &vir_zz,
                                          &vir_yz, &vir_zx, &vir_xy);
    }
  }

#ifndef AR

  /* If we don't use actio=reactio accross the cpus, we have do do
     the force loop also on the other half of the neighbours for the 
     cells on the surface of the CPU */

  /* compute forces for remaining pairs of cells */
  for (n=0; n<nlists; ++n) {
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
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
                tmpvec2, tmpvec2+1, tmpvec2+2, tmpvec2+3, tmpvec2+4,
                                    tmpvec2+5, tmpvec2+6, tmpvec2+7);
    }
  }

#endif /* AR */

  /* sum up results of different CPUs */
  tmpvec1[0] = tot_pot_energy;
  tmpvec1[1] = virial;
  tmpvec1[2] = vir_xx;
  tmpvec1[3] = vir_yy;
  tmpvec1[4] = vir_xy;

  MPI_Allreduce( tmpvec1, tmpvec2, 5, REAL, MPI_SUM, cpugrid); 

  tot_pot_energy = tmpvec2[0];
  virial         = tmpvec2[1];
  vir_xx         = tmpvec2[2];
  vir_yy         = tmpvec2[3];
  vir_xy         = tmpvec2[4];

#ifdef AR
  send_forces(add_forces,pack_forces,unpack_forces);
#endif

}

/******************************************************************************
*
*  fix_cells
*
*  check if each atom is in the correct cell and on the correct CPU; 
*  move atoms that have left their cell or CPU
*
******************************************************************************/

void fix_cells(void)
{
  int i,j,l,clone;
  cell *p, *q;
  ivektor coord, lcoord, dcpu, to_coord;
  msgbuf *buf;

  empty_mpi_buffers();

  /* apply periodic boundary conditions */
  do_boundaries();

  /* for each cell in bulk */
  for (i=cellmin.x; i < cellmax.x; ++i)
    for (j=cellmin.y; j < cellmax.y; ++j) {

      p = PTR_2D_V(cell_array, i, j, cell_dim);

      /* loop over atoms in cell */
      l=0;
      while( l < p->n ) {

        coord = cell_coord( ORT(p,l,X), ORT(p,l,Y) );
        lcoord = local_cell_coord( coord );
	/* see if atom is in wrong cell */
        if ((lcoord.x == i) && (lcoord.y == j)) {
          l++;
        } else {

          /* Calculate distance on CPU grid */
          to_coord = cpu_coord_v( coord );
          dcpu.x = to_coord.x - my_coord.x;
          dcpu.y = to_coord.y - my_coord.y;

          /* Consider PBC */
          if (pbc_dirs.x == 1) {
            if (cpu_dim.x == 1) dcpu.x = 0; 
            else dcpu.x -= ((int) (dcpu.x / (cpu_dim.x/2)) * cpu_dim.x);
          }
          if (pbc_dirs.y == 1) {
            if (cpu_dim.y == 1) dcpu.y = 0;
            else dcpu.y -= ((int) (dcpu.y / (cpu_dim.y/2)) * cpu_dim.y);
          }

          /* Check, if atom is on my cpu */
          /* If not, copy into send buffer else move to correct cell */
          buf   = NULL;
          if      ((0<dcpu.x) && (cpu_dim.x>1)) {
            buf = &send_buf_west; 
          }
          else if ((0>dcpu.x) && (cpu_dim.x>1)) { 
            buf = &send_buf_east;
          }
          else if (0<dcpu.y) { 
            buf = &send_buf_south;
          }
          else if (0>dcpu.y) { 
            buf = &send_buf_north;
          }
          else { /* atom is on my cpu */
            q = PTR_VV(cell_array,lcoord,cell_dim);
            MOVE_ATOM(q, p, l);
#ifdef CLONE
            if (l < p->n-nclones)
              for (clone=1; clone<nclones; clone++) 
                MOVE_ATOM(q, p, l+clone);
            else /* we are dealing with the last in the stack */
              for (clone=1; clone<nclones; clone++) 
                MOVE_ATOM(q, p, l); 
#endif
          }
          if (buf != NULL) {
            int to_cpu = cpu_coord( coord );
            copy_one_atom( buf, to_cpu, p, l, 1); 
#ifdef CLONE
            if (l < p->n-nclones)
              for (clone=1; clone<nclones; clone++)
                copy_one_atom( buf, to_cpu, p, l+clone, 1);
            else /* we are dealing with the last in the stack */
              for (clone=1; clone<nclones; clone++)
                copy_one_atom( buf, to_cpu, p, l, 1);
#endif
	  }
        }
      }
    }

  /* send atoms to neighbbour CPUs */
  send_atoms();

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
    process_buffer( &recv_buf_west );

    /* send west, receive east, move atoms from east to cells */
    sendrecv_buf( &send_buf_west, nbwest, &recv_buf_east, nbeast, &stat );
    MPI_Get_count( &stat, REAL, &recv_buf_east.n );
    process_buffer( &recv_buf_east );

    if (cpu_dim.y > 1) {
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
  }

  if (cpu_dim.y > 1) {
    /* send north, receive south, move atoms from south to cells */
    sendrecv_buf( &send_buf_north, nbnorth, &recv_buf_south, nbsouth, &stat);
    MPI_Get_count( &stat, REAL, &recv_buf_south.n );
    process_buffer( &recv_buf_south );   

    /* send south, receive north, move atoms from north to cells */
    sendrecv_buf( &send_buf_south, nbsouth, &recv_buf_north, nbnorth, &stat);
    MPI_Get_count( &stat, REAL, &recv_buf_north.n );
    process_buffer( &recv_buf_north );   
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
    process_buffer( &recv_buf_west );

    /* wait for atoms from east, move them to cells */
    MPI_Waitall(2, reqeast, stateast);
    MPI_Get_count( &stateast[1], REAL, &recv_buf_east.n );
    process_buffer( &recv_buf_east );

    if (cpu_dim.y > 1) {
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
    process_buffer( &recv_buf_south );   

    /* Wait for atoms from north, move them to cells */
    MPI_Waitall(2, reqnorth, statnorth);
    MPI_Get_count( &statnorth[1], REAL, &recv_buf_north.n );
    process_buffer( &recv_buf_north );   
  }

}

#endif



