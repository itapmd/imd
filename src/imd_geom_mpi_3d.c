
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
* imd_geom_mpi_3d.c -- Geometry routines for MPI, three dimensional version
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/******************************************************************************
*
*  set up the MPI communication topology
*
******************************************************************************/

#ifdef MPI

void setup_mpi_topology( void )
{
  int period[3] = { 1, 1, 1 };
  ivektor cpuc, nbcoord;

  /* Set up process topology */
  if (0==myid) {
    printf("MPI process array dimensions: %d %d %d\n",
	   cpu_dim.x,cpu_dim.y,cpu_dim.z);
  }
  
  MPI_Cart_create(MPI_COMM_WORLD, 3, (int *) &cpu_dim, period, 1, &cpugrid);
  MPI_Comm_rank(cpugrid, &myid);
  MPI_Comm_size(cpugrid, &num_cpus);
  MPI_Cart_coords(cpugrid, myid, 3, (int *) &my_coord);
  
  cpu_ranks = (int *) malloc( cpu_dim.x * cpu_dim.y * cpu_dim.z * sizeof(int));

  for (cpuc.x=0; cpuc.x < cpu_dim.x; ++cpuc.x)
    for (cpuc.y=0; cpuc.y < cpu_dim.y; ++cpuc.y)
      for (cpuc.z=0; cpuc.z < cpu_dim.z; ++cpuc.z) 
        MPI_Cart_rank(cpugrid, (int *) &cpuc, PTR_VV(cpu_ranks,cpuc,cpu_dim));

  /* set up neighbour cpus ids */
  /* faces */
  nbcoord = my_coord; --nbcoord.x; nbeast = cpu_grid_coord( nbcoord );
  nbcoord = my_coord; ++nbcoord.x; nbwest = cpu_grid_coord( nbcoord );
  nbcoord = my_coord; --nbcoord.y; nbnorth= cpu_grid_coord( nbcoord );
  nbcoord = my_coord; ++nbcoord.y; nbsouth= cpu_grid_coord( nbcoord );
  nbcoord = my_coord; --nbcoord.z; nbup   = cpu_grid_coord( nbcoord );
  nbcoord = my_coord; ++nbcoord.z; nbdown = cpu_grid_coord( nbcoord );

  /* edges */
  nbcoord = my_coord; ++nbcoord.x; --nbcoord.y; nbnw = cpu_grid_coord(nbcoord);
  nbcoord = my_coord; ++nbcoord.x; ++nbcoord.y; nbws = cpu_grid_coord(nbcoord);
  nbcoord = my_coord; --nbcoord.x; ++nbcoord.y; nbse = cpu_grid_coord(nbcoord);
  nbcoord = my_coord; --nbcoord.x; --nbcoord.y; nben = cpu_grid_coord(nbcoord);
  nbcoord = my_coord; --nbcoord.y; --nbcoord.z; nbun = cpu_grid_coord(nbcoord);
  nbcoord = my_coord; ++nbcoord.x; --nbcoord.z; nbuw = cpu_grid_coord(nbcoord);
  nbcoord = my_coord; ++nbcoord.y; --nbcoord.z; nbus = cpu_grid_coord(nbcoord);
  nbcoord = my_coord; --nbcoord.x; --nbcoord.z; nbue = cpu_grid_coord(nbcoord);
  nbcoord = my_coord; --nbcoord.y; ++nbcoord.z; nbdn = cpu_grid_coord(nbcoord);
  nbcoord = my_coord; --nbcoord.x; ++nbcoord.z; nbde = cpu_grid_coord(nbcoord);
  nbcoord = my_coord; ++nbcoord.y; ++nbcoord.z; nbds = cpu_grid_coord(nbcoord);
  nbcoord = my_coord; ++nbcoord.x; ++nbcoord.z; nbdw = cpu_grid_coord(nbcoord);

  /* corners */
  nbcoord = my_coord; ++nbcoord.x; --nbcoord.y; --nbcoord.z; nbunw = cpu_grid_coord( nbcoord );
  nbcoord = my_coord; ++nbcoord.x; ++nbcoord.y; --nbcoord.z; nbuws = cpu_grid_coord( nbcoord );
  nbcoord = my_coord; --nbcoord.x; ++nbcoord.y; --nbcoord.z; nbuse = cpu_grid_coord( nbcoord );
  nbcoord = my_coord; --nbcoord.x; --nbcoord.y; --nbcoord.z; nbuen = cpu_grid_coord( nbcoord );
  nbcoord = my_coord; --nbcoord.x; --nbcoord.y; ++nbcoord.z; nbdne = cpu_grid_coord( nbcoord );
  nbcoord = my_coord; --nbcoord.x; ++nbcoord.y; ++nbcoord.z; nbdes = cpu_grid_coord( nbcoord );
  nbcoord = my_coord; ++nbcoord.x; ++nbcoord.y; ++nbcoord.z; nbdsw = cpu_grid_coord( nbcoord );
  nbcoord = my_coord; ++nbcoord.x; --nbcoord.y; ++nbcoord.z; nbdwn = cpu_grid_coord( nbcoord );

  init_io();

}

#else

void setup_mpi_topology(void)
{
  cpu_ranks = (int *) calloc(1, sizeof(int));
  if (NULL==cpu_ranks) error("cannot allocate cpu_ranks");
}

#endif

/******************************************************************************
*
*  local_cell_coord computes the local cell coordinates of a position
*
******************************************************************************/
#ifdef LOADBALANCE
ivektor local_cell_coord(ivektor global_coord)
{
  ivektor cellc;

  cellc.x = global_coord.x - lb_cell_offset.x;
  cellc.y = global_coord.y - lb_cell_offset.y;
  cellc.z = global_coord.z - lb_cell_offset.z;

  return cellc;
}
#else
ivektor local_cell_coord(ivektor global_coord)
{
  ivektor cellc;

  cellc.x = global_coord.x - my_coord.x * (cell_dim.x - 2) + 1;
  cellc.y = global_coord.y - my_coord.y * (cell_dim.y - 2) + 1;
  cellc.z = global_coord.z - my_coord.z * (cell_dim.z - 2) + 1;

  return cellc;
}
#endif


/******************************************************************************
*
*  cpu_coord computes the rank of a CPU from the coordinates of a cell
*
******************************************************************************/

int cpu_coord(ivektor cellc)
{
  ivektor coord;
  
  /* map cell to CPU grid */
  coord.x = (int) (cellc.x * cpu_dim.x / global_cell_dim.x);
  coord.y = (int) (cellc.y * cpu_dim.y / global_cell_dim.y);
  coord.z = (int) (cellc.z * cpu_dim.z / global_cell_dim.z);

  /* return CPU rank */
  return *PTR_3D_VV(cpu_ranks, coord, cpu_dim); 
}


/******************************************************************************
*
*  cpu_grid_coord computes the CPU rank from the CPU coordinates.
*  Used to determine the rank of neighbor CPUs, so periodic
*  boundary conditions must always be applied.
*
******************************************************************************/

int cpu_grid_coord(ivektor cellc)
{
  /* Apply PBC */
  if (cellc.x < 0)          cellc.x += cpu_dim.x;
  if (cellc.x >= cpu_dim.x) cellc.x -= cpu_dim.x;
  if (cellc.y < 0)          cellc.y += cpu_dim.y;
  if (cellc.y >= cpu_dim.y) cellc.y -= cpu_dim.y;
  if (cellc.z < 0)          cellc.z += cpu_dim.z;
  if (cellc.z >= cpu_dim.z) cellc.z -= cpu_dim.z;

  /* return CPU rank */
  return *PTR_3D_VV(cpu_ranks, cellc, cpu_dim); 
}


/******************************************************************************
*
*  cpu_coord_v computes the coordinates of a CPU from the cell coordinates.
*
******************************************************************************/

ivektor cpu_coord_v(ivektor cellc)
{
  ivektor coord;

  /* Map cell to cpugrid */
  coord.x = (int) (cellc.x * cpu_dim.x / global_cell_dim.x);
  coord.y = (int) (cellc.y * cpu_dim.y / global_cell_dim.y);
  coord.z = (int) (cellc.z * cpu_dim.z / global_cell_dim.z);

  return coord;
}

/******************************************************************************
*
*  calc_cpu_dim chooses the CPU dimensions if no array dimensions are given
*  or if cpu_dim does not correspond to the available number of CPUs.
*  The number of CPUs is factorized evenly.
*
******************************************************************************/

void calc_cpu_dim(void)
{
  int trial, n;
  int *cpu_dim_ptr[3];
  int *tmpptr, tmp;
  ivektor fctr;

  /* sort cpu_dim by size. cpu_dim_ptr[0] points to the largest component */
  cpu_dim_ptr[0] = &cpu_dim.x;
  cpu_dim_ptr[1] = &cpu_dim.y;
  cpu_dim_ptr[2] = &cpu_dim.z;
  if ( *cpu_dim_ptr[2] > *cpu_dim_ptr[1] ) {
    tmpptr         = cpu_dim_ptr[1];
    cpu_dim_ptr[1] = cpu_dim_ptr[2];
    cpu_dim_ptr[2] = tmpptr;
  }
  if ( *cpu_dim_ptr[1] > *cpu_dim_ptr[0] ) {
    tmpptr         = cpu_dim_ptr[0];
    cpu_dim_ptr[0] = cpu_dim_ptr[1];    
    cpu_dim_ptr[1] = tmpptr;
  }
  if ( *cpu_dim_ptr[2] > *cpu_dim_ptr[1] ) {
    tmpptr         = cpu_dim_ptr[1];
    cpu_dim_ptr[1] = cpu_dim_ptr[2];
    cpu_dim_ptr[2] = tmpptr;
  }  

   /* factorize num_cpus evenly */
  n = num_cpus;

  trial = (int) ceil(pow(n, 1.0/3.0));

  for( fctr.x=trial; fctr.x>0; fctr.x--)
    if ( n%fctr.x == 0 )
      break;
  n /= fctr.x;

  trial = (int) ceil( sqrt( (double) n ) );

  for( fctr.y=trial; fctr.y>0; fctr.y--)
    if ( n%fctr.y == 0 )
      break;

  fctr.z = n/fctr.y;
  
  /* sort factorization by size */
  if ( fctr.z > fctr.y ) {
    tmp = fctr.y;
    fctr.y = fctr.z;
    fctr.z = tmp;
  }
  if ( fctr.y > fctr.x ) {
    tmp = fctr.x;
    fctr.x = fctr.y;
    fctr.y = tmp;
  }
  if ( fctr.z > fctr.y ) {
    tmp = fctr.y;
    fctr.y = fctr.z;
    fctr.z = tmp;
  }

  *cpu_dim_ptr[0] = fctr.x;
  *cpu_dim_ptr[1] = fctr.y;
  *cpu_dim_ptr[2] = fctr.z;
}
