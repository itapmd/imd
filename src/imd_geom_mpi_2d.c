/******************************************************************************
*
* imd_geom_mpi.c -- Geometry routines for mpi two dimensional version
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/
#include "imd.h"


/* set up the MPI communication topology */

void setup_mpi_topology( void )

{
  int tmp;
  int period[2] = { 1, 1 };
  ivektor cpuc, nbcoord;

/* set up process topology */

  if (0==myid) {
    printf("Global cpu array dimensions: %d %d\n",
	   cpu_dim.x,cpu_dim.y);
  };

  tmp = cpu_dim.x*cpu_dim.y;
  if ( 0 == myid ) {
    printf("Want %d cpus, have %d cpus.\n",(int) tmp,num_cpus);
    if (tmp > num_cpus) error("Not enough cpus."); 
  };

  MPI_Cart_create( MPI_COMM_WORLD, 2, (int *) &cpu_dim, period, 1, &cpugrid );

  MPI_Comm_rank(cpugrid,&myid);
  MPI_Comm_size(cpugrid,&num_cpus);
  MPI_Cart_coords(cpugrid, myid, 2, (int *) &my_coord);

  cpu_ranks = (int *) malloc( cpu_dim.x * cpu_dim.y * sizeof(int));
  if ( 0 == myid )
  if (NULL == cpu_ranks) error("Can't allocate memory for cpu_ranks");

  for (cpuc.x=0; cpuc.x < cpu_dim.x; ++cpuc.x)
    for (cpuc.y=0; cpuc.y < cpu_dim.y; ++cpuc.y)
	MPI_Cart_rank( cpugrid, (int *) &cpuc, 
                       PTR_2D_VV(cpu_ranks, cpuc, cpu_dim));

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


/*
*
* cpu_coord gives the (integral) cpu_rank of a cell
*
*/

int cpu_coord(ivektor cellc)

{
  ivektor coord;
  int rank;
  
  /* Map cell to cpugrid */
  coord.x = (int) (cellc.x * cpu_dim.x / global_cell_dim.x);
  coord.y = (int) (cellc.y * cpu_dim.y / global_cell_dim.y);

  /* Apply PBC */

  if (coord.x < 0)          coord.x += cpu_dim.x;
  if (coord.x >= cpu_dim.x) coord.x -= cpu_dim.x;
  if (coord.y < 0)          coord.y += cpu_dim.y;
  if (coord.y >= cpu_dim.y) coord.y -= cpu_dim.y;

  /* Get Cpu Rank */

  rank = *PTR_2D_VV(cpu_ranks, coord, cpu_dim); 

  return(rank);
}

/*
*
* cpu_grid_coord gives the (integral) cpu_rank from cpu_coords
*
*/

int cpu_grid_coord(ivektor cellc)

{
  int rank;

  /* Apply PBC */

  if (cellc.x < 0)          cellc.x += cpu_dim.x;
  if (cellc.x >= cpu_dim.x) cellc.x -= cpu_dim.x;
  if (cellc.y < 0)          cellc.y += cpu_dim.y;
  if (cellc.y >= cpu_dim.y) cellc.y -= cpu_dim.y;

  /* Get Cpu Rank */

  rank = *PTR_2D_VV(cpu_ranks, cellc, cpu_dim); 
  
  return(rank);
}


/*
*
* local_cell_coord gives local coordinates of a position
*
*/

ivektor local_cell_coord(real x, real y)

{
  ivektor global_coord;
  ivektor local_coord;

  global_coord = cell_coord(x,y);

  local_coord.x = global_coord.x - my_coord.x * (cell_dim.x - 2) + 1;
  local_coord.y = global_coord.y - my_coord.y * (cell_dim.y - 2) + 1;

  return(local_coord);
}

/*
*
* global_cell_coord gives global coordinates of a cell
* PBC are applied.
*
* local cells are [1..max-2] 
* buffer/border cells are [0] and [max-1]
*
*/

ivektor global_cell_coord(ivektor local_coord)
     
{
  
  ivektor global_coord;

  global_coord.x = local_coord.x - 1 + my_coord.x * (cell_dim.x - 2);
  global_coord.y = local_coord.y - 1 + my_coord.y * (cell_dim.y - 2);

  if (global_coord.x < 0)                  global_coord.x += global_cell_dim.x;
  if (global_coord.x >= global_cell_dim.x) global_coord.x -= global_cell_dim.x;
  if (global_coord.y < 0)                  global_coord.y += global_cell_dim.y;
  if (global_coord.y >= global_cell_dim.y) global_coord.y -= global_cell_dim.y;

  return(global_coord);

}


/*
*
* cpu_coord gives the coordinates of a cell 
* does not apply pbc!!!
*
*/

ivektor cpu_coord_v(ivektor cellc)

{
  ivektor coord;

  /* Map cell to cpugrid */
  coord.x = (int) (cellc.x * cpu_dim.x / global_cell_dim.x);
  coord.y = (int) (cellc.y * cpu_dim.y / global_cell_dim.y);
 
  return(coord);
}

