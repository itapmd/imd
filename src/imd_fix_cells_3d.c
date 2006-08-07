
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2006 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_fix_cells_3d.c -- code for fixing cell distribution, three dimensions
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#define INDEXED_ACCESS
#include "imd.h"

/******************************************************************************
*
*  fix_cells
*
*  check if each atom is in the correct cell and on the correct CPU;
*  move atoms that have left their cell or CPU
*
*  this also uses Plimpton's communication scheme
*
******************************************************************************/

void fix_cells(void)
{
  int i,j,k,l,clone,to_cpu;
  minicell *p, *q;
  ivektor coord, lcoord;
  msgbuf *buf;

#ifdef MPI
  empty_mpi_buffers();
#endif

  /* apply periodic boundary conditions */
  do_boundaries();

  /* for each cell in bulk */
  for (i=cellmin.x; i < cellmax.x; ++i)
    for (j=cellmin.y; j < cellmax.y; ++j)
      for (k=cellmin.z; k < cellmax.z; ++k) {

	p = PTR_3D_V(cell_array, i, j, k, cell_dim);

	/* loop over atoms in cell */
	l=0;
	while( l<p->n ) {

          coord  = cell_coord( ORT(p,l,X), ORT(p,l,Y), ORT(p,l,Z) );
	  lcoord = local_cell_coord( coord );

 	  /* see if atom is in wrong cell */
	  if ((lcoord.x == i) && (lcoord.y == j) && (lcoord.z == k)) {
            l++;
          } 
          else {

            to_cpu = cpu_coord(coord);
            buf    = NULL;

            /* atom is on my cpu */
	    if (to_cpu==myid) {
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
#ifdef MPI
            /* west */
            else if ((cpu_dim.x>1) && 
               ((to_cpu==nbwest) || (to_cpu==nbnw)  || (to_cpu==nbws) ||
                (to_cpu==nbuw  ) || (to_cpu==nbunw) || (to_cpu==nbuws)||
                (to_cpu==nbdw  ) || (to_cpu==nbdwn) || (to_cpu==nbdsw))) {
              buf = &send_buf_west;
            }
            
            /* east */
            else if ((cpu_dim.x>1) &&
                ((to_cpu==nbeast) || (to_cpu==nbse)  || (to_cpu==nben) ||
                 (to_cpu==nbue  ) || (to_cpu==nbuse) || (to_cpu==nbuen)||
                 (to_cpu==nbde  ) || (to_cpu==nbdes) || (to_cpu==nbdne))) {
              buf = &send_buf_east;
            }
                   
            /* south  */
            else if ((cpu_dim.y>1) &&
                ((to_cpu==nbsouth) || (to_cpu==nbus)  || (to_cpu==nbds))) {
              buf = &send_buf_south;
            }
                   
            /* north  */
            else if ((cpu_dim.y>1) &&
                ((to_cpu==nbnorth) || (to_cpu==nbun)  || (to_cpu==nbdn))) {
              buf = &send_buf_north;
            }
            
            /* down  */
            else if ((cpu_dim.z>1) && (to_cpu==nbdown)) {
              buf = &send_buf_down;
            }
            
            /* up  */
            else if ((cpu_dim.z>1) && (to_cpu==nbup)) {
              buf = &send_buf_up;
            }
            
            else {
#ifdef SHOCK
              /* remove atom from simulation */
              buf = &dump_buf;
              dump_buf.n = 0;
              natoms  -= nclones;
              nactive -= nclones * DIM;
              num_sort [ SORTE(p,l)] -= nclones;
              num_vsort[VSORTE(p,l)] -= nclones;
              warning("Atom jumped multiple CPUs");
#else
              error("Atom jumped multiple CPUs");
#endif
	    }

            if (buf != NULL) {
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
#endif /* MPI */
	  }
	}
      }

#ifdef MPI
  /* send atoms to neighbbour CPUs */
  send_atoms();
#endif

#ifdef NBLIST
  /* tag neighbor list as outdated */
  have_valid_nbl = 0;
#endif

}

#ifdef MPI

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
    sendrecv_buf( &send_buf_west, nbwest, &recv_buf_east, nbeast, &stat);
    MPI_Get_count( &stat, REAL, &recv_buf_east.n );
    process_buffer( &recv_buf_east );

    if (cpu_dim.y > 1) {
      /* append atoms from west and east to north send buffer */
      copy_atoms_buf( &send_buf_north, &recv_buf_west );
      copy_atoms_buf( &send_buf_north, &recv_buf_east );
      /* check special case cpu_dim.y==2 */ 
      if (nbsouth!=nbnorth) {
        /* append atoms from west and east to south send buffer */
        copy_atoms_buf( &send_buf_south, &recv_buf_west );
        copy_atoms_buf( &send_buf_south, &recv_buf_east );
      }
    }

    if (cpu_dim.z > 1) {
      /* append atoms from west and east to up send buffer */
      copy_atoms_buf( &send_buf_up, &recv_buf_east );
      copy_atoms_buf( &send_buf_up, &recv_buf_west );
      /* check special case cpu_dim.z==2 */ 
      if (nbdown!=nbup) {
        /* append atoms from west and east to down send buffer */
        copy_atoms_buf( &send_buf_down, &recv_buf_east );
        copy_atoms_buf( &send_buf_down, &recv_buf_west );
      }
    }
  }

  if (cpu_dim.y > 1) {
    /* send north, receive south, move atoms from south to cells */
    sendrecv_buf(  &send_buf_north, nbnorth, &recv_buf_south, nbsouth, &stat);
    MPI_Get_count( &stat, REAL, &recv_buf_south.n );
    process_buffer( &recv_buf_south );

    /* send south, receive north, move atoms from north to cells */
    sendrecv_buf( &send_buf_south, nbsouth, &recv_buf_north, nbnorth, &stat);
    MPI_Get_count( &stat, REAL, &recv_buf_north.n );
    process_buffer( &recv_buf_north );

    if (cpu_dim.z > 1) {
      /* append atoms from north and south to up send buffer */
      copy_atoms_buf( &send_buf_up, &recv_buf_north );
      copy_atoms_buf( &send_buf_up, &recv_buf_south );
      /* check special case cpu_dim.z==2 */ 
      if (nbdown!=nbup) {
        /* append atoms from north and south to down send buffer */
        copy_atoms_buf( &send_buf_down, &recv_buf_north );
        copy_atoms_buf( &send_buf_down, &recv_buf_south );
      }
    }
  }

  if (cpu_dim.z > 1) {
    /* send up, receive down, move atoms from down to cells */
    sendrecv_buf( &send_buf_up, nbup, &recv_buf_down, nbdown, &stat);
    MPI_Get_count( &stat, REAL, &recv_buf_down.n );
    process_buffer( &recv_buf_down );
  
    /* send down, receive up, move atoms from up to cells */
    sendrecv_buf( &send_buf_down, nbdown, &recv_buf_up, nbup, &stat );
    MPI_Get_count( &stat, REAL, &recv_buf_up.n );
    process_buffer( &recv_buf_up );
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
  MPI_Status  stateast[2],  statwest[2];
  MPI_Status statnorth[2], statsouth[2];
  MPI_Status    statup[2],  statdown[2];

  MPI_Request  reqeast[2],   reqwest[2];
  MPI_Request reqnorth[2],  reqsouth[2];
  MPI_Request    requp[2],   reqdown[2];

  if (cpu_dim.x > 1) {
    /* send east */
    irecv_buf( &recv_buf_west, nbwest, &reqwest[1] );
    isend_buf( &send_buf_east, nbeast, &reqwest[0] );

    /* send west */
    irecv_buf( &recv_buf_east, nbeast, &reqeast[1] );
    isend_buf( &send_buf_west, nbwest, &reqeast[0] );

    /* Wait for atoms from west, move them to cells */
    MPI_Waitall(2, reqwest, statwest);
    MPI_Get_count( &statwest[1], REAL, &recv_buf_west.n );
    process_buffer( &recv_buf_west );

    /* Wait for atoms from east, move them to cells */
    MPI_Waitall(2, reqeast, stateast);
    MPI_Get_count( &stateast[1], REAL, &recv_buf_east.n );
    process_buffer( &recv_buf_east );

    if (cpu_dim.y > 1) {
      /* append atoms from west and east to north send buffer */
      copy_atoms_buf( &send_buf_north, &recv_buf_west );
      copy_atoms_buf( &send_buf_north, &recv_buf_east );
      /* check special case cpu_dim.y==2 */ 
      if (nbsouth!=nbnorth) {
        /* append atoms from west and east to south send buffer */
        copy_atoms_buf( &send_buf_south, &recv_buf_east );
        copy_atoms_buf( &send_buf_south, &recv_buf_west );
      }
    }

    if (cpu_dim.z > 1) {
      /* append atoms from west and east to up send buffer */
      copy_atoms_buf( &send_buf_up, &recv_buf_east );
      copy_atoms_buf( &send_buf_up, &recv_buf_west );
      /* check special case cpu_dim.z==2 */ 
      if (nbdown!=nbup) {
        /* append atoms from west and east to down send buffer */
        copy_atoms_buf( &send_buf_down, &recv_buf_east );
        copy_atoms_buf( &send_buf_down, &recv_buf_west );
      }
    }
  }

  if (cpu_dim.y > 1) {
    /* Send atoms north */
    irecv_buf( &recv_buf_south, nbsouth, &reqsouth[1] );
    isend_buf( &send_buf_north, nbnorth, &reqsouth[0] );
  
    /* Send atoms south */
    irecv_buf( &recv_buf_north, nbnorth, &reqnorth[1] );
    isend_buf( &send_buf_south, nbsouth, &reqnorth[0] );

    /* Wait for atoms from south, move them to cells */
    MPI_Waitall(2, reqsouth, statsouth);
    MPI_Get_count( &statsouth[1], REAL, &recv_buf_south.n );
    process_buffer( &recv_buf_south );

    /* Wait for atoms from north, move them to cells */
    MPI_Waitall(2, reqnorth, statnorth);
    MPI_Get_count( &statnorth[1], REAL, &recv_buf_north.n );
    process_buffer( &recv_buf_north );

    if (cpu_dim.z > 1) {
      /* append atoms from north and south to up send buffer */
      copy_atoms_buf( &send_buf_up, &recv_buf_north );
      copy_atoms_buf( &send_buf_up, &recv_buf_south );
      /* check special case cpu_dim.z==2 */ 
      if (nbdown!=nbup) {
        /* append atoms from north and south to down send buffer */
        copy_atoms_buf( &send_buf_down, &recv_buf_north );
        copy_atoms_buf( &send_buf_down, &recv_buf_south );
      }
    }
  }

  if (cpu_dim.z > 1) {
    /* Send atoms up */
    irecv_buf( &recv_buf_down , nbdown, &reqdown[1]);
    isend_buf( &send_buf_up   , nbup  , &reqdown[0]);
  
    /* Send atoms down */
    irecv_buf( &recv_buf_up  , nbup  , &requp[1] );
    isend_buf( &send_buf_down, nbdown, &requp[0] );

    /* Wait for atoms from down, move them to cells */
    MPI_Waitall(2, reqdown, statdown);
    MPI_Get_count( &statdown[1], REAL, &recv_buf_down.n );
    process_buffer( &recv_buf_down );

    /* Wait for atoms from up, move them to cells */
    MPI_Waitall(2, requp, statup);
    MPI_Get_count( &statup[1], REAL, &recv_buf_up.n );
    process_buffer( &recv_buf_up );
  }

}

#endif /* not SR */

#endif /* MPI */
