
/******************************************************************************
*
* imd_comm_direct_3d.c -- communication for force computation, 3D
* 
* This version sends also in diagonal directions. It is (at least on
* the T3E) slower than the Plimpton scheme, and thus is obsolete.
* It works only in AR mode (not all cells are communicated).
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"


/******************************************************************************
*
* This sends the atoms in surface cells
* to the neighbouring cpus and copies them to
* the buffer cells
*
* The data sent is:
*
* Faces:   up, south, down, west
* Edges:   up-south, down-south, up-west, north-west, down-west, west-south
* Corners: down-north-west, down-south-west, up-north-west, up-west-south
*
* The data is received on the canonical 'opposite side' of the data sent.
*
* A processor's area:
*    
*  from top    from botton
*
*      N           N
*    |---|       |---|
*  W | U | E    E| D |W    the coordinates origin is in 
*    |---|       |---|     the upper, north, west corner (unw)
*      S           S
*
*  Labeling:
*
*  Faces:    north, west, south, east, up, down
*  Edges:    north-west, west-south, south-east, east-north
*            up-north, up-west, up-south, up-east
*            down-north, down-east, down-south, down-west
*  Corners:  up-north-west, up-west-south, up-south-east, up-east-north
*            down-north-east, down-east-south, down-south-west, down-west-north
*   
*  Abbreviatations are e.g. dne for down-north-east etc.
*
******************************************************************************/

void send_cells_direct(void)
{
  int i,j,k,l;

  MPI_Status  stateast[2],  statwest[2];
  MPI_Status statnorth[2], statsouth[2];
  MPI_Status    statup[2],  statdown[2];

  MPI_Request  reqeast[2],   reqwest[2];
  MPI_Request reqnorth[2],  reqsouth[2];
  MPI_Request    requp[2],   reqdown[2];

  /* Clear buffers */
  empty_mpi_buffers();
  empty_buffer_cells();
  
  /* Exchange faces */
  /* Fill send buffers */
  /* up, down */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.y-1; ++j) {
      (*pack_cell)( &send_buf_up  , i, j,            1 );
      (*pack_cell)( &send_buf_down, i ,j, cell_dim.z-2 );
  }
  /* south, north */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.z-1; ++j) {
      (*pack_cell)( &send_buf_south, i, cell_dim.y-2, j );
      (*pack_cell)( &send_buf_north, i,            1, j );
  }
  /* west, east */
  for (i=1; i < cell_dim.y-1; ++i)
    for (j=1; j < cell_dim.z-1; ++j) {
      (*pack_cell)( &send_buf_west, cell_dim.x-2, i, j );
      (*pack_cell)( &send_buf_east,            1, i, j );  
  }

  /* Exchange faces */
  irecv_buf( &recv_buf_west , nbeast , &reqwest[1] );
  irecv_buf( &recv_buf_south, nbnorth, &reqsouth[1]);
  irecv_buf( &recv_buf_down , nbup   , &reqdown[1] );
  irecv_buf( &recv_buf_up   , nbdown , &requp[1]   );
  irecv_buf( &recv_buf_north, nbsouth, &reqnorth[1]);
  irecv_buf( &recv_buf_east , nbwest , &reqeast[1]);

  isend_buf( &send_buf_west , nbwest  , &reqwest[0] );
  isend_buf( &send_buf_south, nbsouth , &reqsouth[0]);
  isend_buf( &send_buf_down , nbdown  , &reqdown[0] );
  isend_buf( &send_buf_up   , nbup    , &requp[0]   );
  isend_buf( &send_buf_north, nbnorth , &reqnorth[0]);
  isend_buf( &send_buf_east , nbeast  , &reqeast[0] );

  MPI_Waitall(2, reqwest, statwest);
  recv_buf_west.n = 0;
  MPI_Waitall(2, reqsouth, statsouth);
  recv_buf_south.n = 0;
  MPI_Waitall(2, reqdown, statdown);
  recv_buf_down.n = 0;
  MPI_Waitall(2, requp, statup);
  recv_buf_up.n = 0;
  MPI_Waitall(2, reqnorth, statnorth);
  recv_buf_north.n = 0;
  MPI_Waitall(2, reqeast, stateast); 
  recv_buf_east.n = 0;

  /* Add forces from faces */
  /* up, down */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.y-1; ++j) {
      (*unpack_cell)( &recv_buf_up  , i, j, cell_dim.z-1 );
      (*unpack_cell)( &recv_buf_down, i, j,            0 );
  }
  /* south, north */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.z-1; ++j) {
      (*unpack_cell)( &recv_buf_south, i,            0, j);
      (*unpack_cell)( &recv_buf_north, i, cell_dim.y-1, j);
  }
  /* west, east */
  for (i=1; i < cell_dim.y-1; ++i)
    for (j=1; j < cell_dim.z-1; ++j) {
      (*unpack_cell)( &recv_buf_west,            0, i, j );      
      (*unpack_cell)( &recv_buf_east, cell_dim.x-1, i, j ); 
  }

  /* send edges Part 1 */
  /* mapping buffers to edges 

     buffer     send        recv
     west    -  up-south    down-north
     east    -  down-south  up-north
     north   -  up-west     down-east
     south   -  north-west  south-east
     up      -  down-west   up-east
     down    -  west-south  east-north
  */

  empty_mpi_buffers();
  /* north-west, west-south */
  for (i=1; i < cell_dim.z-1; ++i) {
    (*pack_cell)( &send_buf_south, cell_dim.x-2,            1, i);
    (*pack_cell)( &send_buf_down,  cell_dim.x-2, cell_dim.y-2, i);
  }
  /* up-west, down-west */
  for (i=1; i < cell_dim.y-1; ++i) {
    (*pack_cell)( &send_buf_north, cell_dim.x-2, i,            1);
    (*pack_cell)( &send_buf_up   , cell_dim.x-2, i, cell_dim.z-2);
  }
  /* up-south, down-south */
  for (i=1; i < cell_dim.x-1; ++i) {
    (*pack_cell)( &send_buf_west , i, cell_dim.y-2,            1);
    (*pack_cell)( &send_buf_east , i, cell_dim.y-2, cell_dim.z-2);
  }

  /* Exchange edges */
  irecv_buf( &recv_buf_west , nbdn , &reqwest[1] );
  irecv_buf( &recv_buf_east , nbun , &reqeast[1] );
  irecv_buf( &recv_buf_north, nbde , &reqnorth[1]);
  irecv_buf( &recv_buf_south, nbse , &reqsouth[1]);
  irecv_buf( &recv_buf_up   , nbue , &requp[1]   );
  irecv_buf( &recv_buf_down , nben , &reqdown[1] );

  isend_buf( &send_buf_west , nbus , &reqwest[0] );
  isend_buf( &send_buf_east , nbds , &reqeast[0] );
  isend_buf( &send_buf_north, nbuw , &reqnorth[0]);
  isend_buf( &send_buf_south, nbnw , &reqsouth[0]);
  isend_buf( &send_buf_up   , nbdw , &requp[0]   );
  isend_buf( &send_buf_down , nbws , &reqdown[0] );

  MPI_Waitall(2, reqwest , statwest);
  MPI_Waitall(2, reqeast , stateast);
  MPI_Waitall(2, reqnorth, statnorth);
  MPI_Waitall(2, reqsouth, statsouth);
  MPI_Waitall(2, requp   , statup);
  MPI_Waitall(2, reqdown , statdown);

  recv_buf_west.n  = 0;
  recv_buf_east.n  = 0;
  recv_buf_north.n = 0;
  recv_buf_south.n = 0;
  recv_buf_up.n    = 0;
  recv_buf_down.n  = 0;

  /* Add edges */
  /* south-east, east-north  */
  for (i=1; i < cell_dim.z-1; ++i) {
    (*unpack_cell)( &recv_buf_south, 0, cell_dim.y-1, i);
    (*unpack_cell)( &recv_buf_down , 0,            0, i);
  }
  /* down-east, up-east */
  for (i=1; i < cell_dim.y-1; ++i) {
    (*unpack_cell)( &recv_buf_north, 0, i, cell_dim.z-1);
    (*unpack_cell)( &recv_buf_up   , 0, i,            0);
  }
  /* down-north, up-north */
  for (i=1; i < cell_dim.x-1; ++i) {
    (*unpack_cell)( &recv_buf_west, i, 0, cell_dim.z-1);
    (*unpack_cell)( &recv_buf_east, i, 0,            0);
  }

  /* send edges Part 2  */
  /* mapping buffers to edges 

     buffer     send        recv
     west    -  down-north  up-south
     east    -  up-north    down-south
     north   -  south-east  north-west  
     south   -  east-north  west-south 
     up      -  up-east     down-west
     down    -  down-east   up-west
  */

  empty_mpi_buffers();
  /* down-north, up-north */
  for (i=1; i < cell_dim.x-1; ++i) {
    (*pack_cell)( &send_buf_west, i, 1, cell_dim.z-2);
    (*pack_cell)( &send_buf_east, i, 1,            1);
  }
  /* south-east, east-north */
  for (i=1; i < cell_dim.z-1; ++i) { 
    (*pack_cell)( &send_buf_north, 1, cell_dim.y-2, i); 
    (*pack_cell)( &send_buf_south, 1,            1, i); 
  }
  /* up-east, down-east */ 
  for (i=1; i < cell_dim.y-1; ++i) { 
    (*pack_cell)( &send_buf_down, 1,  i,  cell_dim.z-2); 
    (*pack_cell)( &send_buf_up  , 1,  i,             1); 
  }

  /* Exchange edges */
  irecv_buf( &recv_buf_west , nbus , &reqwest[1] );
  irecv_buf( &recv_buf_east , nbds , &reqeast[1] );
  irecv_buf( &recv_buf_north, nbnw , &reqnorth[1]); 
  irecv_buf( &recv_buf_south, nbws , &reqsouth[1]); 
  irecv_buf( &recv_buf_up   , nbdw , &requp[1]   ); 
  irecv_buf( &recv_buf_down , nbuw , &reqdown[1] ); 

  isend_buf( &send_buf_west , nbdn , &reqwest[0] );
  isend_buf( &send_buf_east , nbun , &reqeast[0] );
  isend_buf( &send_buf_north, nbse , &reqnorth[0]); 
  isend_buf( &send_buf_south, nben , &reqsouth[0]); 
  isend_buf( &send_buf_up   , nbue , &requp[0]   ); 
  isend_buf( &send_buf_down , nbde , &reqdown[0] ); 

  MPI_Waitall(2, reqwest , statwest);
  MPI_Waitall(2, reqeast , stateast);
  MPI_Waitall(2, reqnorth, statnorth);
  MPI_Waitall(2, reqsouth, statsouth);
  MPI_Waitall(2, requp   , statup);
  MPI_Waitall(2, reqdown , statdown);

  recv_buf_west.n  = 0;
  recv_buf_east.n  = 0;
  recv_buf_north.n = 0;
  recv_buf_south.n = 0;
  recv_buf_up.n    = 0;
  recv_buf_down.n  = 0;

  /* Add edges */
  /* down-north, up-north */
  for (i=1; i < cell_dim.x-1; ++i) {
    (*unpack_cell)( &recv_buf_west, i, cell_dim.y-1,            0);
    (*unpack_cell)( &recv_buf_east, i, cell_dim.y-1, cell_dim.z-1);
  }
  /* south-east, east-north */
  for (i=1; i < cell_dim.z-1; ++i) {
    (*unpack_cell)( &recv_buf_north, cell_dim.x-1,            0, i);
    (*unpack_cell)( &recv_buf_south, cell_dim.x-1, cell_dim.y-1, i);
  }
  /* up-east, down-east */
  for (i=1; i < cell_dim.y-1; ++i) {
    (*unpack_cell)( &recv_buf_down, cell_dim.x-1, i,            0); 
    (*unpack_cell)( &recv_buf_up  , cell_dim.x-1, i, cell_dim.z-1); 
  }

  /* send corners */
  /* mapping buffers - to corners

     buffer    send             recv 
     west   -  up-north-west    down-east-south
     east   -  up-west-south    down-north-east
     south  -  down-south-west  up-east-north
     north  -  down-west-north  up-south-east

  */

  empty_mpi_buffers();
  /* unw, uws, dsw, dwn */

  (*pack_cell)( &send_buf_west , 1,            1,            1 );
  (*pack_cell)( &send_buf_east , 1, cell_dim.y-2,            1 );
  (*pack_cell)( &send_buf_south, 1, cell_dim.y-2, cell_dim.z-2 );
  (*pack_cell)( &send_buf_north, 1,            1, cell_dim.z-2 );

  /* Exchange corners */
  irecv_buf( &recv_buf_west ,  nbdsw , &reqwest[1] );
  irecv_buf( &recv_buf_east ,  nbdwn , &reqeast[1] );
  irecv_buf( &recv_buf_south,  nbunw , &reqsouth[1]);
  irecv_buf( &recv_buf_north,  nbuws , &reqnorth[1]);

  isend_buf( &send_buf_west ,  nbuen , &reqwest[0] );
  isend_buf( &send_buf_east ,  nbuse , &reqeast[0] );
  isend_buf( &send_buf_south,  nbdes , &reqsouth[0]);
  isend_buf( &send_buf_north,  nbdne , &reqnorth[0]);

  MPI_Waitall(2, reqwest, statwest);
  recv_buf_west.n = 0;
  MPI_Waitall(2, reqeast, stateast);
  recv_buf_east.n = 0;
  MPI_Waitall(2, reqnorth, statnorth);
  recv_buf_north.n = 0;
  MPI_Waitall(2, reqsouth, statsouth);
  recv_buf_south.n = 0;

  /* add corners */
  (*unpack_cell)( &recv_buf_west ,cell_dim.x-1, cell_dim.y-1, cell_dim.z-1);
  (*unpack_cell)( &recv_buf_east ,cell_dim.x-1,            0, cell_dim.z-1);
  (*unpack_cell)( &recv_buf_south,cell_dim.x-1,            0, 0);
  (*unpack_cell)( &recv_buf_north,cell_dim.x-1, cell_dim.y-1, 0);

}


/******************************************************************************
*
* This sends the forces accumulated by atoms in buffer cells
* back to the original cpus of these atoms and adds this 
* forces to each atoms total.
*
* The data sent is:
*
* Faces:   up, south, down, west
* Edges:   up-south, down-south, up-west, north-west, down-west, west-south
* Corners: down-north-west, down-south-west, up-north-west, up-west-south
*
* The data is received on the canonical 'opposite side' of the data sent.
*
* A processor's area:
*    
*  from top    from botton
*
*      N           N
*    |---|       |---|
*  W | U | E    E| D |W    the coordinates origin is in 
*    |---|       |---|     the upper, north, west corner (unw)
*      S           S
*
*  Labeling:
*
*  Faces:    north, west, south, east, up, down
*  Edges:    north-west, west-south, south-east, east-north
*            up-north, up-west, up-south, up-east
*            down-north, down-east, down-south, down-west
*  Corners:  up-north-west, up-west-south, up-south-east, up-east-north
*            down-north-east, down-east-south, down-south-west, down-west-north
*   
*  Abbreviatations are e.g. dne for down-north-east etc.
*
******************************************************************************/

void send_forces_direct(void)
{
  int i,j,k,l;

  MPI_Status  stateast[2],  statwest[2];
  MPI_Status statnorth[2], statsouth[2];
  MPI_Status    statup[2],  statdown[2];

  MPI_Request  reqeast[2],   reqwest[2];
  MPI_Request reqnorth[2],  reqsouth[2];
  MPI_Request    requp[2],   reqdown[2];

  /* Clear buffers */
  empty_mpi_buffers();

  /* Fill send buffers */
  /* up, down */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.y-1; ++j) {
      (*pack_forces)( &send_buf_up  , i, j,            0 );
      (*pack_forces)( &send_buf_down, i ,j, cell_dim.z-1 );
  }
  /* south, north */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.z-1; ++j) {
      (*pack_forces)( &send_buf_south, i, cell_dim.y-1, j );
      (*pack_forces)( &send_buf_north, i,            0, j );
  }
  /* west */
  for (i=1; i < cell_dim.y-1; ++i)
    for (j=1; j < cell_dim.z-1; ++j) {
      (*pack_forces)( &send_buf_west, cell_dim.x-1, i, j );
      /* (*pack_forces)( &send_buf_east, 0, i, j );  */
  }

  /* Exchange faces */
  irecv_buf( &recv_buf_west , nbeast , &reqwest[1] );
  irecv_buf( &recv_buf_south, nbnorth, &reqsouth[1]);
  irecv_buf( &recv_buf_down , nbup   , &reqdown[1] );
  irecv_buf( &recv_buf_up   , nbdown , &requp[1]   );
  irecv_buf( &recv_buf_north, nbsouth, &reqnorth[1]);
  /* irecv_buf( &recv_buf_east, nbwest, &reqeast[1] ); */

  isend_buf( &send_buf_west , nbwest  , &reqwest[0] );
  isend_buf( &send_buf_south, nbsouth , &reqsouth[0]);
  isend_buf( &send_buf_down , nbdown  , &reqdown[0] );
  isend_buf( &send_buf_up   , nbup    , &requp[0]   );
  isend_buf( &send_buf_north, nbnorth , &reqnorth[0]);
  /* isend_buf( &send_buf_east, nbeast, &reqeast[0]); */

  MPI_Waitall(2, reqwest, statwest);
  recv_buf_west.n = 0;
  MPI_Waitall(2, reqsouth, statsouth);
  recv_buf_south.n = 0;
  MPI_Waitall(2, reqdown, statdown);
  recv_buf_down.n = 0;
  MPI_Waitall(2, requp, statup);
  recv_buf_up.n = 0;
  MPI_Waitall(2, reqnorth, statnorth);
  recv_buf_north.n = 0;
  /* MPI_Waitall(2, reqeast, stateast); */
  recv_buf_east.n = 0;

  /* Add forces from faces */
  /* up, down */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.y-1; ++j) {
      (*unpack_forces)( &recv_buf_up  , i, j, cell_dim.z-2 );
      (*unpack_forces)( &recv_buf_down, i, j,            1 );
  }
  /* south, north */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.z-1; ++j) {
      (*unpack_forces)( &recv_buf_south, i,            1, j);
      (*unpack_forces)( &recv_buf_north, i, cell_dim.y-2, j);
  }
  /* west, east */
  for (i=1; i < cell_dim.y-1; ++i)
    for (j=1; j < cell_dim.z-1; ++j) {
      (*unpack_forces)( &recv_buf_west,            1, i, j );      
     /* (*unpack_forces)( &recv_buf_east, cell_dim.x-2, i, j ); */
  }

  /* send edges Part 1 */
  /* mapping buffers to edges 

     buffer     send        recv
     west    -  up-south    down-north
     east    -  down-south  up-north
     north   -  up-west     down-east
     south   -  north-west  south-east
     up      -  down-west   up-east
     down    -  west-south  east-north
  */

  /* fill buffers */
  empty_mpi_buffers();
  /* north-west, west-south */
  for (i=1; i < cell_dim.z-1; ++i) {
    (*pack_forces)( &send_buf_south, cell_dim.x-1,            0, i);
    (*pack_forces)( &send_buf_down,  cell_dim.x-1, cell_dim.y-1, i);
  }
  /* up-west, down-west */
  for (i=1; i < cell_dim.y-1; ++i) {
    (*pack_forces)( &send_buf_north, cell_dim.x-1, i,            0);
    (*pack_forces)( &send_buf_up   , cell_dim.x-1, i, cell_dim.z-1);
  }
  /* up-south, down-south */
  for (i=1; i < cell_dim.x-1; ++i) {
    (*pack_forces)( &send_buf_west , i, cell_dim.y-1,            0);
    (*pack_forces)( &send_buf_east , i, cell_dim.y-1, cell_dim.z-1);
  }

  /* Exchange edges */
  irecv_buf( &recv_buf_west , nbdn , &reqwest[1] );
  irecv_buf( &recv_buf_east , nbun , &reqeast[1] );
  irecv_buf( &recv_buf_north, nbde , &reqnorth[1]);
  irecv_buf( &recv_buf_south, nbse , &reqsouth[1]);
  irecv_buf( &recv_buf_up   , nbue , &requp[1]   );
  irecv_buf( &recv_buf_down , nben , &reqdown[1] );

  isend_buf( &send_buf_west , nbus , &reqwest[0] );
  isend_buf( &send_buf_east , nbds , &reqeast[0] );
  isend_buf( &send_buf_north, nbuw , &reqnorth[0]);
  isend_buf( &send_buf_south, nbnw , &reqsouth[0]);
  isend_buf( &send_buf_up   , nbdw , &requp[0]   );
  isend_buf( &send_buf_down , nbws , &reqdown[0] );

  MPI_Waitall(2, reqwest , statwest );
  MPI_Waitall(2, reqeast , stateast );
  MPI_Waitall(2, reqnorth, statnorth);
  MPI_Waitall(2, reqsouth, statsouth);
  MPI_Waitall(2, requp   , statup   );
  MPI_Waitall(2, reqdown , statdown );

  recv_buf_west.n  = 0;
  recv_buf_east.n  = 0;
  recv_buf_north.n = 0;
  recv_buf_south.n = 0;
  recv_buf_up.n    = 0;
  recv_buf_down.n  = 0;

  /* Add edges */
  /* south-east, east-north  */
  for (i=1; i < cell_dim.z-1; ++i) {
    (*unpack_forces)( &recv_buf_south, 1, cell_dim.y-2, i);
    (*unpack_forces)( &recv_buf_down , 1,            1, i);
  }
  /* down-east, up-east */
  for (i=1; i < cell_dim.y-1; ++i) {
    (*unpack_forces)( &recv_buf_north, 1, i, cell_dim.z-2);
    (*unpack_forces)( &recv_buf_up   , 1, i,            1);
  }
  /* down-north, up-north */
  for (i=1; i < cell_dim.x-1; ++i) {
    (*unpack_forces)( &recv_buf_west, i, 1, cell_dim.z-2);
    (*unpack_forces)( &recv_buf_east, i, 1,            1);
  }

  /* send edges Part 2  */
  /* mapping buffers to edges 

     buffer     send        recv
     west    -  down-north  up-south
     east    -  up-north    down-south
     north   -  south-east  north-west  
     south   -  east-north  west-south 
     up      -  up-east     down-west
     down    -  down-east   up-west
  */

  empty_mpi_buffers();
  /* down-north, up-north */
  for (i=1; i < cell_dim.x-1; ++i) {
    (*pack_forces)( &send_buf_west, i, 0, cell_dim.z-1);
    (*pack_forces)( &send_buf_east, i, 0,            0);
  }

  /* Exchange edges */
  irecv_buf( &recv_buf_west , nbus , &reqwest[1]);
  irecv_buf( &recv_buf_east , nbds , &reqeast[1]);
  isend_buf( &send_buf_west , nbdn , &reqwest[0]);
  isend_buf( &send_buf_east , nbun , &reqeast[0]);

  MPI_Waitall(2, reqwest , statwest);
  MPI_Waitall(2, reqeast , stateast);

  recv_buf_west.n  = 0;
  recv_buf_east.n  = 0;

  /* Add edges */
  /* down-north, up-north */
  for (i=1; i < cell_dim.x-1; ++i) {
    (*unpack_forces)( &recv_buf_west, i, cell_dim.y-2,            1);
    (*unpack_forces)( &recv_buf_east, i, cell_dim.y-2, cell_dim.z-2);
  }

  /* send corners */
  /* mapping buffers - to corners

     buffer    send             recv 
     west   -  up-north-west    down-east-south
     east   -  up-west-south    down-north-east
     south  -  down-south-west  up-east-north
     north  -  down-west-north  up-south-east

  */

  empty_mpi_buffers();

  /* unw, uws, dsw, dwn */
  (*pack_forces)( &send_buf_west , cell_dim.x-1,            0,            0 );
  (*pack_forces)( &send_buf_east , cell_dim.x-1, cell_dim.y-1,            0 );
  (*pack_forces)( &send_buf_south, cell_dim.x-1, cell_dim.y-1, cell_dim.z-1 );
  (*pack_forces)( &send_buf_north, cell_dim.x-1,            0, cell_dim.z-1 );

  /* Exchange corners */
  irecv_buf( &recv_buf_west , nbdes, &reqwest[1] );
  irecv_buf( &recv_buf_east , nbdne, &reqeast[1] );
  irecv_buf( &recv_buf_south, nbuen, &reqsouth[1]);
  irecv_buf( &recv_buf_north, nbuse, &reqnorth[1]);

  isend_buf( &send_buf_west , nbunw, &reqwest[0] );
  isend_buf( &send_buf_east , nbuws, &reqeast[0] );
  isend_buf( &send_buf_south, nbdsw, &reqsouth[0]);
  isend_buf( &send_buf_north, nbdwn, &reqnorth[0]);  

  MPI_Waitall(2, reqwest, statwest );
  recv_buf_west.n = 0;
  MPI_Waitall(2, reqeast, stateast );
  recv_buf_east.n = 0;
  MPI_Waitall(2, reqnorth, statnorth );
  recv_buf_north.n = 0;
  MPI_Waitall(2, reqsouth, statsouth );
  recv_buf_south.n = 0;

  /* add corners */
  (*unpack_forces)( &recv_buf_west , 1, cell_dim.y-2, cell_dim.z-2);
  (*unpack_forces)( &recv_buf_east , 1,            1, cell_dim.z-2);
  (*unpack_forces)( &recv_buf_south, 1,            1, 1);
  (*unpack_forces)( &recv_buf_north, 1, cell_dim.y-2, 1);

}




