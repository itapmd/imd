
/******************************************************************************
*
* imd_comm_force_3d.c -- communication for force computation, three dimensions
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

#ifdef SR
/******************************************************************************
*
*  Send cells into buffer cells on neighboring CPUs for the force computation.
*  What exactly is sent is determined by the parameter functions.
*  We use Steve Plimptons communication scheme: we send only along
*  the main axis of the system, so that edge cells travel twice,
*  and corner cells three times. In AR mode, one cell wall (including
*  adjacent edge and corner cells) is not needed in the buffer cells.
*
******************************************************************************/

void send_cells(void (*copy_func)  (int, int, int, int, int, int),
                void (*pack_func)  (msgbuf*, int, int, int),
                void (*unpack_func)(msgbuf*, int, int, int))
{
  int i,j;

  MPI_Status  stateast[1],  statwest[1];
  MPI_Status statnorth[1], statsouth[1];
  MPI_Status    statup[1],  statdown[1];

  empty_mpi_buffers();

  /* exchange up/down */
  if (cpu_dim.z==1) {
    /* simply copy up/down atoms to buffer cells */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j) {
        (*copy_func)( i, j, 1, i, j, cell_dim.z-1 );
        (*copy_func)( i, j, cell_dim.z-2, i, j, 0 );
      }
  } else {
    /* copy up atoms into send buffer, send up */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_up, i, j, 1 );
    isendrecv_buf(&send_buf_up, nbup, &recv_buf_down, nbdown, statup);

   
    /* copy down atoms into send buffer, send down */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_down, i, j, cell_dim.z-2 );
    isendrecv_buf(&send_buf_down, nbdown, &recv_buf_up, nbup, statdown);
 
    /* atoms from down, move them to buffer cells */
    recv_buf_down.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*unpack_func)( &recv_buf_down, i, j, cell_dim.z-1 );

    /* atoms from up, move them to buffer cells */
    recv_buf_up.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*unpack_func)( &recv_buf_up, i, j, 0 );
  }

  /* exchange north/south */
  if (cpu_dim.y==1) {
    /* simply copy north/south atoms to buffer cells */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j) {
        (*copy_func)( i, 1, j, i, cell_dim.y-1, j );
        (*copy_func)( i, cell_dim.y-2, j, i, 0, j );
      }
  } else {
    /* copy north atoms into send buffer, send north */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_north, i, 1, j );
    isendrecv_buf(&send_buf_north, nbnorth, &recv_buf_south, nbsouth, statnorth);

    /* copy south atoms into send buffer, send south */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_south, i, cell_dim.y-2, j );
    isendrecv_buf(&send_buf_south, nbsouth, &recv_buf_north, nbnorth, statsouth);

    /* atoms from south, move them to buffer cells */
    recv_buf_south.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_south, i, cell_dim.y-1, j );

    /* atoms from north, move them to buffer cells */
    recv_buf_north.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_north, i, 0, j );
  }

   
  /* exchange east/west */
  if (cpu_dim.x==1) {
    /* simply copy east/west atoms to buffer cells */
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) { 
        (*copy_func)( 1, i, j, cell_dim.x-1, i, j );
#if !defined(AR) || defined(COVALENT) 
        (*copy_func)( cell_dim.x-2, i, j, 0, i, j );
#endif
      }
  } else {
    /* copy east atoms into send buffer, send east */
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) 
        (*pack_func)( &send_buf_east, 1, i, j );
    isendrecv_buf(&send_buf_east, nbeast, &recv_buf_west, nbwest, stateast);



#if !defined(AR) || defined(COVALENT) 
    /* copy west atoms into send buffer, send west */
    for (i=0; i < cell_dim.y; ++i) 
      for (j=0; j < cell_dim.z; ++j) 
        (*pack_func)( &send_buf_west, cell_dim.x-2, i, j );
    isendrecv_buf(&send_buf_west, nbwest, &recv_buf_east, nbeast, statwest);
#endif

    /* wait for atoms from west, move them to buffer cells */
    recv_buf_west.n = 0;
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) 
        (*unpack_func)( &recv_buf_west, cell_dim.x-1, i, j );

#if !defined(AR) || defined(COVALENT) 
    /* atoms from east, move them to buffer cells */
    recv_buf_east.n = 0;
    for (i=0; i < cell_dim.y; ++i) 
      for (j=0; j < cell_dim.z; ++j) 
        (*unpack_func)( &recv_buf_east, 0, i, j );
#endif
  }  
}
#endif

#ifndef SR
/******************************************************************************
*
*  Send cells into buffer cells on neighboring CPUs for the force computation.
*  What exactly is sent is determined by the parameter functions.
*  We use Steve Plimptons communication scheme: we send only along
*  the main axis of the system, so that edge cells travel twice,
*  and corner cells three times. In AR mode, one cell wall (including
*  adjacent edge and corner cells) is not needed in the buffer cells.
*
******************************************************************************/

void send_cells(void (*copy_func)  (int, int, int, int, int, int),
                void (*pack_func)  (msgbuf*, int, int, int),
                void (*unpack_func)(msgbuf*, int, int, int))
{
  int i,j;

  MPI_Status  stateast[2],  statwest[2];
  MPI_Status statnorth[2], statsouth[2];
  MPI_Status    statup[2],  statdown[2];

  MPI_Request  reqeast[2],   reqwest[2];
  MPI_Request reqnorth[2],  reqsouth[2];
  MPI_Request    requp[2],   reqdown[2];

  empty_mpi_buffers();

  /* exchange up/down */
  if (cpu_dim.z==1) {
    /* simply copy up/down atoms to buffer cells*/
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j) {
        (*copy_func)( i, j, 1, i, j, cell_dim.z-1 );
        (*copy_func)( i, j, cell_dim.z-2, i, j, 0 );
      }
  } else {
    /* copy up atoms into send buffer, send up */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_up, i, j, 1 );
    irecv_buf( &recv_buf_down , nbdown, &reqdown[1]);
    isend_buf( &send_buf_up   , nbup  , &reqdown[0]);

    /* copy down atoms into send buffer, send down */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_down, i, j, cell_dim.z-2 );
    irecv_buf( &recv_buf_up  , nbup  , &requp[1] );
    isend_buf( &send_buf_down, nbdown, &requp[0] );

    /* wait for atoms from down, move them to buffer cells */
    MPI_Waitall(2, reqdown, statdown);
    recv_buf_down.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*unpack_func)( &recv_buf_down, i, j, cell_dim.z-1 );

    /* wait for atoms from up, move them to buffer cells*/
    MPI_Waitall(2, requp, statup);
    recv_buf_up.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*unpack_func)( &recv_buf_up, i, j, 0 );
  }

  /* exchange north/south */
  if (cpu_dim.y==1) {
    /* simply copy north/south atoms to buffer cells */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j) {
        (*copy_func)( i, 1, j, i, cell_dim.y-1, j );
        (*copy_func)( i, cell_dim.y-2, j, i, 0, j );
      }
  } else {
    /* copy north atoms into send buffer, send north*/
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_north, i, 1, j );
    irecv_buf( &recv_buf_south, nbsouth, &reqsouth[1] );
    isend_buf( &send_buf_north, nbnorth, &reqsouth[0] );

    /* copy south atoms into send buffer, send south*/
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_south, i, cell_dim.y-2, j );
    irecv_buf( &recv_buf_north, nbnorth, &reqnorth[1] );
    isend_buf( &send_buf_south, nbsouth, &reqnorth[0] );

    /* wait for atoms from south, move them to buffer cells */
    MPI_Waitall(2, reqsouth, statsouth);
    recv_buf_south.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_south, i, cell_dim.y-1, j );

    /* wait for atoms from north, move them to buffer cells*/
    MPI_Waitall(2, reqnorth, statnorth);
    recv_buf_north.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_north, i, 0, j );
  }

  /* exchange east/west*/
  if (cpu_dim.x==1) {
  /* simply copy east/west atoms to buffer cells*/
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) { 
        (*copy_func)( 1, i, j, cell_dim.x-1, i, j );
#if !defined(AR) || defined(COVALENT) 
        (*copy_func)( cell_dim.x-2, i, j, 0, i, j );
#endif
      }
  } else {
  /* copy east atoms into send buffer, send east*/
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) 
        (*pack_func)( &send_buf_east, 1, i, j );
    irecv_buf( &recv_buf_west, nbwest, &reqwest[1] );
    isend_buf( &send_buf_east, nbeast, &reqwest[0] );

#if !defined(AR) || defined(COVALENT) 
    /* copy west atoms into send buffer, send west*/
    for (i=0; i < cell_dim.y; ++i) 
      for (j=0; j < cell_dim.z; ++j) 
        (*pack_func)( &send_buf_west, cell_dim.x-2, i, j );
    irecv_buf( &recv_buf_east, nbeast, &reqeast[1] );
    isend_buf( &send_buf_west, nbwest, &reqeast[0] );
#endif

    /* wait for atoms from west, move them to buffer cells*/
    MPI_Waitall(2, reqwest, statwest);
    recv_buf_west.n = 0;
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) 
        (*unpack_func)( &recv_buf_west, cell_dim.x-1, i, j );

#if !defined(AR) || defined(COVALENT) 
    /* wait for atoms from east, move them to buffer cells*/
    MPI_Waitall(2, reqeast, stateast);
    recv_buf_east.n = 0;
    for (i=0; i < cell_dim.y; ++i) 
      for (j=0; j < cell_dim.z; ++j) 
        (*unpack_func)( &recv_buf_east, 0, i, j );
#endif
  }  
}
#endif

#ifdef SR
/******************************************************************************
*
*  Add forces in buffer cells on neighboring CPUs back to the original cells.
*  What exactly is sent is determined by the parameter functions.
*  We use Steve Plimptons communication scheme: we send only along
*  the main axis of the system, so that edge cells travel twice,
*  and corner cells three times. In AR mode, one buffer cell wall 
*  (including adjacent edge and corner cells) contains no forces.
*
******************************************************************************/

void send_forces(void (*add_func)   (int, int, int, int, int, int),
                 void (*pack_func)  (msgbuf*, int, int, int),
                 void (*unpack_func)(msgbuf*, int, int, int))
{
  int i,j;

  MPI_Status  stateast[2],  statwest[2];
  MPI_Status statnorth[2], statsouth[2];
  MPI_Status    statup[2],  statdown[2];

  empty_mpi_buffers();

  /* send forces east/west */
  if (cpu_dim.x==1) {
    /* simply add east/west forces to original cells */
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) { 
#if !defined(AR) || defined(COVALENT) 
        (*add_func)( 0, i, j, cell_dim.x-2, i, j );
#endif
        (*add_func)( cell_dim.x-1, i, j, 1, i, j );
      }
  } else {
#if !defined(AR) || defined(COVALENT) 
    /* copy east forces into send buffer, send east */
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_east, 0, i, j );
    isendrecv_buf(&send_buf_east, nbeast, &recv_buf_west, nbwest, stateast);
#endif

    /* copy west forces into send buffer, send west */
    for (i=0; i < cell_dim.y; ++i) 
      for (j=0; j < cell_dim.z; ++j) 
        (*pack_func)( &send_buf_west, cell_dim.x-1, i, j );
    isendrecv_buf(&send_buf_west, nbwest, &recv_buf_east, nbeast, statwest);

#if !defined(AR) || defined(COVALENT) 
    /* wait for forces from west, add them to original cells */
    recv_buf_west.n = 0;
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) 
        (*unpack_func)( &recv_buf_west, cell_dim.x-2, i, j );
#endif

    /* wait for forces from east, add them to original cells */
    recv_buf_east.n = 0;
    for (i=0; i < cell_dim.y; ++i) 
      for (j=0; j < cell_dim.z; ++j) 
        (*unpack_func)( &recv_buf_east, 1, i, j );
  }

  /* send forces north/south */
  if (cpu_dim.y==1) {
    /* simply add north/south forces to original cells */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j) {
        (*add_func)( i, 0, j, i, cell_dim.y-2, j );
        (*add_func)( i, cell_dim.y-1, j, i, 1, j );
      }
  } else {
    /* copy north forces into send buffer, send north */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_north, i, 0, j );
    isendrecv_buf(&send_buf_north, nbnorth, &recv_buf_south, nbsouth, statnorth);

    /* copy south forces into send buffer, send south */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_south, i, cell_dim.y-1, j );
    isendrecv_buf(&send_buf_south, nbsouth, &recv_buf_north, nbnorth, statsouth);

    /* wait for forces from south, add them to original cells */
    recv_buf_south.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_south, i, cell_dim.y-2, j );

    /* wait for forces from north, add them to original cells */
    recv_buf_north.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_north, i, 1, j );
  }

  /* send forces up/down */
  if (cpu_dim.z==1) {
    /* simply add up/down forces to original cells */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j) {
        (*add_func)( i, j, 0, i, j, cell_dim.z-2 );
        (*add_func)( i, j, cell_dim.z-1, i, j, 1 );
      }
  } else {
    /* copy up forces into send buffer, send up */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_up, i, j, 0 );
    isendrecv_buf(&send_buf_up, nbup, &recv_buf_down, nbdown, statup);
  
    /* copy down forces into send buffer, send down */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_down, i, j, cell_dim.z-1 );
    isendrecv_buf(&send_buf_down, nbdown, &recv_buf_up, nbup, statdown);

    /* wait for forces from down, add them to original cells */
    recv_buf_down.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*unpack_func)( &recv_buf_down, i, j, cell_dim.z-2 );

    /* wait for forces from up, add them to original cells */
    recv_buf_up.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*unpack_func)( &recv_buf_up, i, j, 1 );
  }
}
#endif 

#ifndef SR
/******************************************************************************
*
*  Add forces in buffer cells on neighboring CPUs back to the original cells.
*  What exactly is sent is determined by the parameter functions.
*  We use Steve Plimptons communication scheme: we send only along
*  the main axis of the system, so that edge cells travel twice,
*  and corner cells three times. In AR mode, one buffer cell wall 
*  (including adjacent edge and corner cells) contains no forces.
*
******************************************************************************/

void send_forces(void (*add_func)   (int, int, int, int, int, int),
                 void (*pack_func)  (msgbuf*, int, int, int),
                 void (*unpack_func)(msgbuf*, int, int, int))
{
  int i,j;

  MPI_Status  stateast[2],  statwest[2];
  MPI_Status statnorth[2], statsouth[2];
  MPI_Status    statup[2],  statdown[2];

  MPI_Request  reqeast[2],   reqwest[2];
  MPI_Request reqnorth[2],  reqsouth[2];
  MPI_Request    requp[2],   reqdown[2];

  empty_mpi_buffers();

  /* send forces east/west */
  if (cpu_dim.x==1) {
    /* simply add east/west forces to original cells */
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) { 
#if !defined(AR) || defined(COVALENT) 
        (*add_func)( 0, i, j, cell_dim.x-2, i, j );
#endif
        (*add_func)( cell_dim.x-1, i, j, 1, i, j );
      }
  } else {
#if !defined(AR) || defined(COVALENT) 
    /* copy east forces into send buffer, send east */
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_east, 0, i, j );
    irecv_buf( &recv_buf_west, nbwest, &reqwest[1] );
    isend_buf( &send_buf_east, nbeast, &reqwest[0] );
#endif

    /* copy west forces into send buffer, send west */
    for (i=0; i < cell_dim.y; ++i) 
      for (j=0; j < cell_dim.z; ++j) 
        (*pack_func)( &send_buf_west, cell_dim.x-1, i, j );
    irecv_buf( &recv_buf_east, nbeast, &reqeast[1] );
    isend_buf( &send_buf_west, nbwest, &reqeast[0] );

#if !defined(AR) || defined(COVALENT) 
    /* wait for forces from west, add them to original cells */
    MPI_Waitall(2, reqwest, statwest);
    recv_buf_west.n = 0;
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) 
        (*unpack_func)( &recv_buf_west, cell_dim.x-2, i, j );
#endif

    /* wait for forces from east, add them to original cells */
    MPI_Waitall(2, reqeast, stateast);
    recv_buf_east.n = 0;
    for (i=0; i < cell_dim.y; ++i) 
      for (j=0; j < cell_dim.z; ++j) 
        (*unpack_func)( &recv_buf_east, 1, i, j );
  }

  /* send forces north/south */
  if (cpu_dim.y==1) {
    /* simply add north/south forces to original cells */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j) {
        (*add_func)( i, 0, j, i, cell_dim.y-2, j );
        (*add_func)( i, cell_dim.y-1, j, i, 1, j );
      }
  } else {
    /* copy north forces into send buffer, send north */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_north, i, 0, j );
    irecv_buf( &recv_buf_south, nbsouth, &reqsouth[1] );
    isend_buf( &send_buf_north, nbnorth, &reqsouth[0] );

    /* copy south forces into send buffer, send south */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_south, i, cell_dim.y-1, j );
    irecv_buf( &recv_buf_north, nbnorth, &reqnorth[1] );
    isend_buf( &send_buf_south, nbsouth, &reqnorth[0] );

    /* wait for forces from south, add them to original cells */
    MPI_Waitall(2, reqsouth, statsouth);
    recv_buf_south.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_south, i, cell_dim.y-2, j );

    /* wait for forces from north, add them to original cells */
    MPI_Waitall(2, reqnorth, statnorth);
    recv_buf_north.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_north, i, 1, j );
  }

  /* send forces up/down */
  if (cpu_dim.z==1) {
    /* simply add up/down forces to original cells */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j) {
        (*add_func)( i, j, 0, i, j, cell_dim.z-2 );
        (*add_func)( i, j, cell_dim.z-1, i, j, 1 );
      }
  } else {
    /* copy up forces into send buffer, send up */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_up, i, j, 0 );
    irecv_buf( &recv_buf_down , nbdown, &reqdown[1]);
    isend_buf( &send_buf_up   , nbup  , &reqdown[0]);

    /* copy down forces into send buffer, send down */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_down, i, j, cell_dim.z-1 );
    irecv_buf( &recv_buf_up  , nbup  , &requp[1] );
    isend_buf( &send_buf_down, nbdown, &requp[0] );

    /* wait for forces from down, add them to original cells */
    MPI_Waitall(2, reqdown, statdown);
    recv_buf_down.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*unpack_func)( &recv_buf_down, i, j, cell_dim.z-2 );

    /* wait for forces from up, add them to original cells */
    MPI_Waitall(2, requp, statup);
    recv_buf_up.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*unpack_func)( &recv_buf_up, i, j, 1 );
  }
}
#endif

/******************************************************************************
*
*  copy contents of one cell to another (buffer) cell (for force comp.)
*
******************************************************************************/

void copy_cell( int k, int l, int m, int r, int s, int t )
{
  int i, tmp_n;
  cell *from, *to;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  to   = PTR_3D_V(cell_array, r, s, t, cell_dim);

  tmp_n = from->n;
  if (tmp_n > to->n_max) {
    to->n = 0;
    alloc_cell(to, tmp_n);
  }
  
  to->n = tmp_n;
  for (i=0; i<to->n; ++i) {
    to->ort X(i) = from->ort X(i);
    to->ort Y(i) = from->ort Y(i);
    to->ort Z(i) = from->ort Z(i);
#ifndef MONOLJ
    to->sorte[i] = from->sorte[i];
#ifdef UNIAX
    to->achse X(i) = from->achse X(i);
    to->achse Y(i) = from->achse Y(i);
    to->achse Z(i) = from->achse Z(i);
    to->shape X(i) = from->shape X(i);
    to->shape Y(i) = from->shape Y(i);
    to->shape Z(i) = from->shape Z(i);
    to->pot_well X(i) = from->pot_well X(i);
    to->pot_well Y(i) = from->pot_well Y(i);
    to->pot_well Z(i) = from->pot_well Z(i);
#endif
#endif 
#ifdef EAM
    to->nummer[i] = from->nummer[i];
#endif
  }
}


/******************************************************************************
*
*  pack cell into MPI send buffer (for force comp.)
*
******************************************************************************/

void pack_cell( msgbuf *b, int k, int l, int m )
{
  int i;
  cell *from;
    
  from = PTR_3D_V(cell_array, k, l, m, cell_dim);

  b->data[ b->n++ ] = (real) from->n;
    
  for (i=0; i<from->n; ++i) {
    b->data[ b->n++ ] = from->ort X(i);
    b->data[ b->n++ ] = from->ort Y(i);
    b->data[ b->n++ ] = from->ort Z(i);
#ifndef MONOLJ
    b->data[ b->n++ ] = (real) from->sorte[i];
#ifdef UNIAX
    b->data[ b->n++ ] = from->achse X(i);
    b->data[ b->n++ ] = from->achse Y(i);
    b->data[ b->n++ ] = from->achse Z(i);
    b->data[ b->n++ ] = from->shape X(i);
    b->data[ b->n++ ] = from->shape Y(i);
    b->data[ b->n++ ] = from->shape Z(i);
    b->data[ b->n++ ] = from->pot_well X(i);
    b->data[ b->n++ ] = from->pot_well Y(i);
    b->data[ b->n++ ] = from->pot_well Z(i);
#endif
#endif
#ifdef EAM
    b->data[ b->n++ ] = (real) from->nummer[i];
#endif
  }
  if (b->n_max < b->n)  error("Buffer overflow in pack_cell");
}


/******************************************************************************
*
*  unpack cell from MPI buffer to buffer cell (for force comp.)
*
******************************************************************************/

void unpack_cell( msgbuf *b, int k, int l, int m )
{
  int i;
  int tmp_n;
  cell *to;

  to = PTR_3D_V(cell_array, k, l, m, cell_dim);

  tmp_n = (int) b->data[ b->n++ ];
  
  if (tmp_n > to->n_max) {
    to->n = 0;
    alloc_cell(to, tmp_n);
  }
  
  to->n = tmp_n;
  for (i=0; i<to->n; ++i) {
    to->ort X(i) = b->data[ b->n++ ];
    to->ort Y(i) = b->data[ b->n++ ];
    to->ort Z(i) = b->data[ b->n++ ];
#ifndef MONOLJ
    to->sorte[i] = (shortint) b->data[ b->n++ ];
#ifdef UNIAX
    to->achse X(i) = b->data[ b->n++ ];
    to->achse Y(i) = b->data[ b->n++ ];
    to->achse Z(i) = b->data[ b->n++ ];
    to->shape X(i) = b->data[ b->n++ ];
    to->shape Y(i) = b->data[ b->n++ ];
    to->shape Z(i) = b->data[ b->n++ ];
    to->pot_well X(i) = b->data[ b->n++ ];
    to->pot_well Y(i) = b->data[ b->n++ ];
    to->pot_well Z(i) = b->data[ b->n++ ];
#endif
#endif 
#ifdef EAM
    to->nummer[i] = (integer) b->data[ b->n++ ];
#endif
  }
  if (b->n_max < b->n) error("Buffer overflow in unpack_cell");
}


/******************************************************************************
*
*  add forces of one cell to those of another cell
*
******************************************************************************/

void add_forces( int k, int l, int m, int r, int s, int t )
{
  int i;
  cell *from, *to;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  to   = PTR_3D_V(cell_array, r, s, t, cell_dim);

  for (i=0; i<to->n; ++i) {
    to->kraft X(i) += from->kraft X(i);
    to->kraft Y(i) += from->kraft Y(i);
    to->kraft Z(i) += from->kraft Z(i);
#ifndef MONOLJ
    to->pot_eng[i] += from->pot_eng[i];
#endif
#ifdef NVX
    to->heatcond[i] += from->heatcond[i];
#endif
#ifdef STRESS_TENS
    to->presstens X(i) += from->presstens X(i);
    to->presstens Y(i) += from->presstens Y(i);
    to->presstens Z(i) += from->presstens Z(i);
    to->presstens_offdia X(i) += from->presstens_offdia X(i);
    to->presstens_offdia Y(i) += from->presstens_offdia Y(i);
    to->presstens_offdia Z(i) += from->presstens_offdia Z(i);
#endif
#ifdef ORDPAR
    to->nbanz[i] += from->nbanz[i];
#endif
  }
}


/******************************************************************************
*
*  pack forces from buffer cell into MPI buffer
*
******************************************************************************/

void pack_forces( msgbuf *b, int k, int l, int m)
 {
  int i;
  cell *from;
    
  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  for (i=0; i<from->n; ++i) {
    b->data[ b->n++ ] = from->kraft X(i);
    b->data[ b->n++ ] = from->kraft Y(i);
    b->data[ b->n++ ] = from->kraft Z(i);
#ifndef MONOLJ
    b->data[ b->n++ ] = from->pot_eng[i];
#endif
#ifdef NVX
    b->data[ b->n++ ] = from->heatcond[i];
#endif
#ifdef STRESS_TENS
    b->data[ b->n++ ] = from->presstens X(i);
    b->data[ b->n++ ] = from->presstens Y(i);
    b->data[ b->n++ ] = from->presstens Z(i);
    b->data[ b->n++ ] = from->presstens_offdia X(i);
    b->data[ b->n++ ] = from->presstens_offdia Y(i);
    b->data[ b->n++ ] = from->presstens_offdia Z(i);
#endif
#ifdef ORDPAR
    b->data[ b->n++ ] = (real) from->nbanz[i];
#endif
  }
  if (b->n_max < b->n) error("Buffer overflow in pack_forces.");
}


/******************************************************************************
*
*  unpack forces from MPI buffer, and add them to those of the original cell
*
******************************************************************************/

void unpack_forces( msgbuf *b, int k, int l, int m )
 {
  int i;
  cell *to;

  to = PTR_3D_V(cell_array, k, l, m, cell_dim);
  for (i=0; i<to->n; ++i) {
    to->kraft X(i) += b->data[ b->n++ ];
    to->kraft Y(i) += b->data[ b->n++ ];
    to->kraft Z(i) += b->data[ b->n++ ];
#ifndef MONOLJ
    to->pot_eng[i] += b->data[ b->n++ ];
#endif
#ifdef NVX
    to->heatcond[i] += b->data[ b->n++ ];
#endif
#ifdef STRESS_TENS
    to->presstens X(i)        += b->data[ b->n++ ];
    to->presstens Y(i)        += b->data[ b->n++ ];
    to->presstens Z(i)        += b->data[ b->n++ ];
    to->presstens_offdia X(i) += b->data[ b->n++ ];
    to->presstens_offdia Y(i) += b->data[ b->n++ ];
    to->presstens_offdia Z(i) += b->data[ b->n++ ];
#endif
#ifdef ORDPAR
    to->nbanz[i] += (shortint) b->data[ b->n++ ];
#endif
  }
  if (b->n_max < b->n) error("Buffer overflow in unpack_forces.");
}


#ifdef EAM2

/******************************************************************************
*
*  copy eam2_rho_h of one cell to another cell
*
******************************************************************************/

void copy_rho_h( int k, int l, int m, int r, int s, int t )
{
  int i;
  cell *from, *to;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  to   = PTR_3D_V(cell_array, r, s, t, cell_dim);

  for (i=0; i<to->n; ++i) {
    to->eam2_rho_h[i] = from->eam2_rho_h[i];
  }
}


/******************************************************************************
*
*  pack eam2_rho_h into MPI buffer
*
******************************************************************************/

void pack_rho_h( msgbuf *b, int k, int l, int m )
{
  int i;
  cell *from;
    
  from = PTR_3D_V(cell_array, k, l, m, cell_dim);

  for (i=0; i<from->n; ++i) {
    b->data[ b->n++ ] = from->eam2_rho_h[i];
  }
}


/******************************************************************************
*
*  unpack eam2_rho_h from MPI buffer into cell
*
******************************************************************************/

void unpack_rho_h( msgbuf *b, int k, int l, int m )
{
  int i;
  cell *to;

  to = PTR_3D_V(cell_array, k, l, m, cell_dim);

  for (i=0; i<to->n; ++i) {
    to->eam2_rho_h[i] = b->data[ b->n++ ];
  }
}

#endif /* EAM2 */
