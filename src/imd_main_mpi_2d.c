
/******************************************************************************
*
* imd_main_mpi_2d.c -- main loop, MPI specific part, two dimensions
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
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
* In 2D, we don't use actio=reactio across CPUs yet
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
#ifdef TRANSPORT
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

  /* Since we don't use actio=reactio accross the cpus, we have do do
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

  /* sum up results of different CPUs */
  MPI_Allreduce( &tot_pot_energy, &tmp, 1, MPI_REAL, MPI_SUM, cpugrid); 
  tot_pot_energy = tmp; 

  MPI_Allreduce( &virial,   &tmpvir,    1, MPI_REAL, MPI_SUM, cpugrid);
  virial = tmpvir;

  MPI_Allreduce( &vir_vect, &tmpvec,  DIM, MPI_REAL, MPI_SUM, cpugrid);
  vir_vect.x = tmpvec.x;
  vir_vect.y = tmpvec.y;

}


/******************************************************************************
*
*  global_pbc tells if a local buffer cell is across the boundaries
*
*  no longer used
*
******************************************************************************/

vektor global_pbc(int i, int j)
{
  ivektor global_coord;
  ivektor local_coord;
  vektor pbc = {0.0, 0.0};

  local_coord.x = i;
  local_coord.y = j;

  /* Cannot use global_cell_coord function, this includes already pbc */
  global_coord.x = local_coord.x - 1 + my_coord.x * (cell_dim.x - 2);
  global_coord.y = local_coord.y - 1 + my_coord.y * (cell_dim.y - 2);

  if (global_coord.x < 0) {
    pbc.x -= box_x.x;      
    pbc.y -= box_x.y;
  };

  if (global_coord.x >= global_cell_dim.x) {
    pbc.x += box_x.x;      
    pbc.y += box_x.y;
  };

  if (global_coord.y < 0) {
    pbc.x -= box_y.x;      
    pbc.y -= box_y.y;
  };

  if (global_coord.y >= global_cell_dim.y) {
    pbc.x += box_y.x;      
    pbc.y += box_y.y;
  };

  return pbc;

}


/******************************************************************************
*
* set up mpi
*
******************************************************************************/

void init_mpi(int argc,char *argv[])
{
  /* Initialize MPI */
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&num_cpus);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  if (0 == myid) { 
    printf("%s\n", argv[0]);
    printf("Starting up MPI on %d nodes.\n",num_cpus);
  }

  /* Setup send/receive buffers */
  send_buf_north.n_max = 0;
  send_buf_south.n_max = 0;
  send_buf_east.n_max  = 0;
  send_buf_west.n_max  = 0;
  recv_buf_north.n_max = 0;
  recv_buf_south.n_max = 0;
  recv_buf_east.n_max  = 0;
  recv_buf_west.n_max  = 0;

}


/******************************************************************************
*
* shut down mpi
*
******************************************************************************/

void shutdown_mpi(void)
{
  MPI_Barrier(MPI_COMM_WORLD);   /* Wait for all processes to arrive */
  MPI_Finalize();                /* Shutdown */
}


/******************************************************************************
*
* send_atoms  -  only used for fix_cells; we always have mode == ATOMS
*
******************************************************************************/

void send_atoms(int mode)
{
  MPI_Status  stateast[2], statwest[2], statnorth[2], statsouth[2];
  MPI_Request  reqeast[2],  reqwest[2],  reqnorth[2],  reqsouth[2];

  /* send east */
  isend_buf( &send_buf_east, nbeast, &reqwest[0] );
  irecv_buf( &recv_buf_west, nbwest, &reqwest[1] );

  /* send west */
  isend_buf( &send_buf_west, nbwest, &reqeast[0] );
  irecv_buf( &recv_buf_east, nbeast, &reqeast[1] );

  /* wait for atoms from west, move them to cells */
  MPI_Waitall(2, reqwest, statwest);
  MPI_Get_count( &statwest[1], MPI_REAL, &recv_buf_west.n );
  process_buffer( &recv_buf_west, mode );

  /* wait for atoms from east, move them to cells */
  MPI_Waitall(2, reqeast, stateast);
  MPI_Get_count( &stateast[1], MPI_REAL, &recv_buf_east.n );
  process_buffer( &recv_buf_east, mode );

  /* append atoms from east & west to north send buffer */
  copy_atoms_buf( &send_buf_north, &recv_buf_west );
  copy_atoms_buf( &send_buf_north, &recv_buf_east );
  /* check special case cpu_dim.y==2 */ 
  if (nbsouth!=nbnorth) {
    /* append atoms from east & west to south send buffer */
    copy_atoms_buf( &send_buf_south, &recv_buf_east );
    copy_atoms_buf( &send_buf_south, &recv_buf_west );
  }

  /* send atoms north */
  isend_buf( &send_buf_north, nbnorth, &reqsouth[0] );
  irecv_buf( &recv_buf_south, nbsouth, &reqsouth[1] );

  /* send atoms south */
  isend_buf( &send_buf_south, nbsouth, &reqnorth[0] );
  irecv_buf( &recv_buf_north, nbnorth, &reqnorth[1] );

  /* wait for atoms from south, move them to cells */
  MPI_Waitall(2, reqsouth, statsouth);
  MPI_Get_count( &statsouth[1], MPI_REAL, &recv_buf_south.n );
  process_buffer( &recv_buf_south, mode );   

  /* Wait for atoms from north, move them to cells */
  MPI_Waitall(2, reqnorth, statnorth);
  MPI_Get_count( &statnorth[1], MPI_REAL, &recv_buf_north.n );
  process_buffer( &recv_buf_north, mode );   

}


/******************************************************************************
*
*  copy_cell_force  -  copy contents of one cell to another cell
*
******************************************************************************/

void copy_cell_force( int j, int k, int l, int m)
{
  int i, tmp_n;
  cell *from, *to;

  from = PTR_2D_V(cell_array, j, k, cell_dim);
  to   = PTR_2D_V(cell_array, l, m, cell_dim);

  tmp_n = from->n;
  if (tmp_n >= to->n_max) {
    to->n = 0;
    alloc_cell(to, tmp_n);
  }
  
  to->n = tmp_n;
  for (i=0; i<to->n; ++i) {
    to->ort X(i) = from->ort X(i);
    to->ort Y(i) = from->ort Y(i);
    to->sorte[i] = from->sorte[i];
  }
}


/******************************************************************************
*
*  move_atoms_force  -  move atoms from MPI buffer to cells
*
******************************************************************************/

void move_atoms_force( msgbuf *b, int j, int k)
{
  int i;
  int tmp_n;
  cell *to;

  to = PTR_2D_V(cell_array, j, k, cell_dim);

  tmp_n = (int) b->data[ b->n++ ];
  
  if (tmp_n >= to->n_max) {
    to->n = 0;
    alloc_cell(to, tmp_n);
  }
  
  to->n = tmp_n;
  for (i=0; i<to->n; ++i) {
    to->ort X(i) = b->data[ b->n++ ];
    to->ort Y(i) = b->data[ b->n++ ];
    to->sorte[i] = (shortint) b->data[ b->n++ ];
  }
  if (b->n_max <= b->n) error("Buffer overflow in move_atoms_force");
}


/******************************************************************************
*
*  copy_atoms_force  -  copy atoms to MPI buffer
*
******************************************************************************/

void copy_atoms_force( msgbuf *b, int j, int k)
{
  int i;
  cell *from;
    
  from = PTR_2D_V(cell_array, j, k, cell_dim);

  b->data[ b->n++ ] = (real) from->n;
    
  for (i=0; i<from->n; ++i) {
    b->data[ b->n++ ] = from->ort X(i);
    b->data[ b->n++ ] = from->ort Y(i);
    b->data[ b->n++ ] = (real) from->sorte[i];
  }
  if (b->n_max <= b->n)  error("Buffer overflow in copy_atoms_force");
}


/******************************************************************************
*
* send_atoms_force  -  send atoms for force computation (no AR)
*
******************************************************************************/

void send_atoms_force()
{
  int i;

  MPI_Status  stateast[2], statwest[2], statnorth[2], statsouth[2];
  MPI_Request  reqeast[2],  reqwest[2],  reqnorth[2],  reqsouth[2];

  empty_mpi_buffers();
  empty_buffer_cells();

  /* exchange east/west */
  if (cpu_dim.x==1) {
    /* simply copy east/west atoms to buffer cells */
    for (i=cellmin.y; i<cellmax.y; ++i) {
      copy_cell_force( 1, i, cell_dim.x-1, i );
      copy_cell_force( cell_dim.x-2, i, 0, i );
    }
  } else {
    /* copy east atoms into send buffer, send east */
    for (i=cellmin.y; i<cellmax.y; ++i)
      copy_atoms_force( &send_buf_east, 1, i );
    isend_buf( &send_buf_east, nbeast, &reqwest[0] );
    irecv_buf( &recv_buf_west, nbwest, &reqwest[1] );

    /* copy west atoms into send buffer, send west */
    for (i=cellmin.y; i<cellmax.y; ++i) 
      copy_atoms_force( &send_buf_west, cell_dim.x-2, i);
    /* send west */
    isend_buf( &send_buf_west, nbwest, &reqeast[0] );
    irecv_buf( &recv_buf_east, nbeast, &reqeast[1] );

    /* Wait for atoms from west, move them to buffer cells */
    MPI_Waitall(2, reqwest, statwest);
    recv_buf_west.n = 0;
    for (i=cellmin.y; i<cellmax.y; ++i)
      move_atoms_force( &recv_buf_west, cell_dim.x-1, i );

    /* Wait for atoms from east, move them to buffer cells */
    MPI_Waitall(2, reqeast, stateast);
    recv_buf_east.n = 0;
    for (i=cellmin.y; i<cellmax.y; ++i) 
      move_atoms_force( &recv_buf_east, 0, i);
  }

  /* exchange north/south */
  if (cpu_dim.y==1) {
    /* simply copy north/south atoms to buffer cells */
    for (i=0; i<cell_dim.x; ++i) {
      copy_cell_force( i, 1, i, cell_dim.y-1 );
      copy_cell_force( i, cell_dim.y-2, i, 0 );
    }
  } else {
    /* copy north atoms into send buffer, send north */
    for (i=0; i<cell_dim.x; ++i)
      copy_atoms_force( &send_buf_north, i, 1);
    isend_buf( &send_buf_north, nbnorth, &reqsouth[0] );
    irecv_buf( &recv_buf_south, nbsouth, &reqsouth[1] );

    /* copy south atoms into send buffer, send south */
    for (i=0; i<cell_dim.x; ++i)
      copy_atoms_force( &send_buf_south, i, cell_dim.y-2);
    isend_buf( &send_buf_south, nbsouth, &reqnorth[0] );
    irecv_buf( &recv_buf_north, nbnorth, &reqnorth[1] );

    /* wait for atoms from south, move them to buffer cells */
    MPI_Waitall(2, reqsouth, statsouth);
    recv_buf_south.n = 0;
    for (i=0; i<cell_dim.x; ++i)
      move_atoms_force( &recv_buf_south, i, cell_dim.y-1);

    /* Wait for atoms from north, move them to buffer cells */
    MPI_Waitall(2, reqnorth, statnorth);
    recv_buf_north.n = 0;
    for (i=0; i<cell_dim.x; ++i)
      move_atoms_force( &recv_buf_north, i, 0);
  }
}





