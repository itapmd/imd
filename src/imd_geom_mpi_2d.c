
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_geom_mpi_2d.c -- Geometry routines for MPI, two dimensional version
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

void setup_mpi_topology( void )
{
  int tmp;
  int period[2] = { 1, 1 };
  ivektor cpuc, nbcoord;

  /* set up process topology */
  if (0==myid) {
    printf("Global cpu array dimensions: %d %d\n", cpu_dim.x, cpu_dim.y);
    tmp = cpu_dim.x * cpu_dim.y;
    printf("Want %d cpus, have %d cpus.\n", (int) tmp, num_cpus);
    if (tmp > num_cpus) error("Not enough cpus."); 
  }

  MPI_Cart_create( MPI_COMM_WORLD, 2, (int *) &cpu_dim, period, 1, &cpugrid );
  MPI_Comm_rank(cpugrid,&myid);
  MPI_Comm_size(cpugrid,&num_cpus);
  MPI_Cart_coords(cpugrid, myid, 2, (int *) &my_coord);

  cpu_ranks = (int *) malloc( cpu_dim.x * cpu_dim.y * sizeof(int));
  for (cpuc.x=0; cpuc.x < cpu_dim.x; ++cpuc.x)
    for (cpuc.y=0; cpuc.y < cpu_dim.y; ++cpuc.y)
      MPI_Cart_rank( cpugrid, (int *) &cpuc, PTR_VV(cpu_ranks, cpuc, cpu_dim));

  /* set up neighbour cpus ids */
  nbcoord.x = my_coord.x - 1;
  nbcoord.y = my_coord.y;
  nbeast = cpu_grid_coord( nbcoord );

  nbcoord.x = my_coord.x + 1;
  nbcoord.y = my_coord.y;
  nbwest = cpu_grid_coord( nbcoord );

  nbcoord.x = my_coord.x;
  nbcoord.y = my_coord.y - 1;
  nbnorth = cpu_grid_coord( nbcoord );

  nbcoord.x = my_coord.x;
  nbcoord.y = my_coord.y + 1;
  nbsouth = cpu_grid_coord( nbcoord );

}


/******************************************************************************
*
*  local_cell_coord computes the local cell coordinates of a position
*
******************************************************************************/

ivektor local_cell_coord(real x, real y)
{
  ivektor global_coord, local_coord;

  global_coord  = cell_coord(x,y);
  local_coord.x = global_coord.x - my_coord.x * (cell_dim.x - 2) + 1;
  local_coord.y = global_coord.y - my_coord.y * (cell_dim.y - 2) + 1;

  return local_coord;
}


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

  /* return CPU rank */
  return *PTR_2D_VV(cpu_ranks, coord, cpu_dim); 
}


/******************************************************************************
*
*  cpu_grid_coord computes the CPU rank from the CPU coordinates.
*  Used to determine the rank of neighbor CPUs, so periodic
*  boundary conditions must alway be applied.
*
******************************************************************************/

int cpu_grid_coord(ivektor cellc)
{
  /* Apply PBC */
  if (cellc.x < 0)          cellc.x += cpu_dim.x;
  if (cellc.x >= cpu_dim.x) cellc.x -= cpu_dim.x;
  if (cellc.y < 0)          cellc.y += cpu_dim.y;
  if (cellc.y >= cpu_dim.y) cellc.y -= cpu_dim.y;

  /* return CPU rank */
  return *PTR_2D_VV(cpu_ranks, cellc, cpu_dim); 
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
 
  return coord;
}
