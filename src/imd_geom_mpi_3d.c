/******************************************************************************
*
* imd_geom_mpi.c -- Geometry routines for mpi
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
  int period[3] = { 1, 1, 1 };
  ivektor cpuc, nbcoord;

/* Set Up process topology */

  if (0==myid) {
    printf("Global cpu array dimensions: %d %d %d\n",
	   cpu_dim.x,cpu_dim.y,cpu_dim.z);
  };
  
  tmp = cpu_dim.x*cpu_dim.y*cpu_dim.z;
  if ( 0 == myid ) {
    printf("Want %d cpus, have %d cpus.\n",(int) tmp,num_cpus);
    if (tmp > num_cpus) error("Not enough cpus."); 
  };

#ifndef PACX
  MPI_Cart_create( MPI_COMM_WORLD, 3, (int *) &cpu_dim, period, 1, &cpugrid );
  MPI_Comm_rank(cpugrid,&myid);
  MPI_Comm_size(cpugrid,&num_cpus);
  MPI_Cart_coords(cpugrid, myid, 3, (int *) &my_coord);
#else
  /* kopiere Communicator */
  cpugrid=MPI_COMM_WORLD;
  /* ermittle int myid, unnoetig, da MPI_COMM_WOLRD = cpu_grid */
  /* MPI_Comm_rank(cpugrid,&myid); */
  /* MPI_Comm_size(cpugrid,&num_cpus); */
  /* berechne aus int myid eigene CPU Koordinaten ivektor my_coord*/
  my_coord = my_cart_coords( myid );
#endif
  
  cpu_ranks = (int *) malloc( cpu_dim.x * cpu_dim.y * cpu_dim.z * sizeof(int));
  if ( 0 == myid )
    if (NULL == cpu_ranks) error("Can't allocate memory for cpu_ranks");

      for (cpuc.x=0; cpuc.x < cpu_dim.x; ++cpuc.x)
         for (cpuc.y=0; cpuc.y < cpu_dim.y; ++cpuc.y)
            for (cpuc.z=0; cpuc.z < cpu_dim.z; ++cpuc.z) 

#ifndef PACX
              MPI_Cart_rank( cpugrid, (int *) &cpuc, 
                            PTR_3D_VV(cpu_ranks, cpuc, cpu_dim));
#else
  /* durchlaufe alle CPUs (ivektor cpuc) und berechne pointer cpu_ranks */
              my_cart_rank( cpuc );
#endif

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
  
  /* Map cell to cpugrid*/
  coord.x = (int) (cellc.x * cpu_dim.x / global_cell_dim.x);
  coord.y = (int) (cellc.y * cpu_dim.y / global_cell_dim.y);
  coord.z = (int) (cellc.z * cpu_dim.z / global_cell_dim.z);

  /* Apply PBC */

  if (coord.x < 0)          coord.x += cpu_dim.x;
  if (coord.x >= cpu_dim.x) coord.x -= cpu_dim.x;
  if (coord.y < 0)          coord.y += cpu_dim.y;
  if (coord.y >= cpu_dim.y) coord.y -= cpu_dim.y;
  if (coord.z < 0)          coord.z += cpu_dim.z;
  if (coord.z >= cpu_dim.z) coord.z -= cpu_dim.z;

  /* Get Cpu Rank */

  /*  MPI_Cart_rank( cpugrid, (int *) &coord, &rank);  */

  rank = *PTR_3D_VV(cpu_ranks, coord, cpu_dim); 

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
  if (cellc.z < 0)          cellc.z += cpu_dim.z;
  if (cellc.z >= cpu_dim.z) cellc.z -= cpu_dim.z;

  /* Get Cpu Rank */

  rank = *PTR_3D_VV(cpu_ranks, cellc, cpu_dim); 
  
  return(rank);
}


/*
*
* local_cell_coord gives local coordinates of a position
*
*/

ivektor local_cell_coord(real x, real y, real z)

{
  ivektor global_coord;
  ivektor local_coord;
  
  global_coord = cell_coord(x,y,z);

  local_coord.x = global_coord.x - my_coord.x * (cell_dim.x - 2) + 1;
  local_coord.y = global_coord.y - my_coord.y * (cell_dim.y - 2) + 1;
  local_coord.z = global_coord.z - my_coord.z * (cell_dim.z - 2) + 1;

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
  global_coord.z = local_coord.z - 1 + my_coord.z * (cell_dim.z - 2);


  if (global_coord.x < 0)                  global_coord.x += global_cell_dim.x;
  if (global_coord.x >= global_cell_dim.x) global_coord.x -= global_cell_dim.x;
  if (global_coord.y < 0)                  global_coord.y += global_cell_dim.y;
  if (global_coord.y >= global_cell_dim.y) global_coord.y -= global_cell_dim.y;
  if (global_coord.z < 0)                  global_coord.z += global_cell_dim.z;
  if (global_coord.z >= global_cell_dim.z) global_coord.z -= global_cell_dim.z;

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
  coord.z = (int) (cellc.z * cpu_dim.z / global_cell_dim.z);

  return(coord);
}






