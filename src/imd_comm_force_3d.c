
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
* imd_comm_force_3d.c -- communication for force computation, three dimensions
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#define INDEXED_ACCESS
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

void send_cells(void (*copy_func)  (int, int, int, int, int, int, vektor),
                void (*pack_func)  (msgbuf*, int, int, int, vektor),
                void (*unpack_func)(msgbuf*, int, int, int))
{
  int i,j;

  vektor uvec={0,0,0}, dvec={0,0,0};
  vektor nvec={0,0,0}, svec={0,0,0};
  vektor evec={0,0,0}, wvec={0,0,0};

#ifdef MPI
  MPI_Status  stat;
  empty_mpi_buffers();
#endif

#ifdef VEC
  atoms.n_buf = atoms.n;
#endif
#ifdef NBLIST
  if (pbc_dirs.x==1) {
    if (my_coord.x==0) evec = box_x;
    if (my_coord.x==cpu_dim.x-1) {
      wvec.x = -box_x.x; wvec.y = -box_x.y; wvec.z = -box_x.z;
    }    
  }
  if (pbc_dirs.y==1) {
    if (my_coord.y==0) nvec = box_y;
    if (my_coord.y==cpu_dim.y-1) {
      svec.x = -box_y.x; svec.y = -box_y.y; svec.z = -box_y.z;
    }    
  }
  if (pbc_dirs.z==1) {
    if (my_coord.z==0) uvec = box_z;
    if (my_coord.z==cpu_dim.z-1) {
      dvec.x = -box_z.x; dvec.y = -box_z.y; dvec.z = -box_z.z;
    }    
  }
#endif

  /* exchange up/down */
  if (cpu_dim.z==1) {
    /* simply copy up/down atoms to buffer cells */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j) {
        (*copy_func)( i, j, 1, i, j, cell_dim.z-1, uvec );
        (*copy_func)( i, j, cell_dim.z-2, i, j, 0, dvec );
      }
  }
#ifdef MPI
  else {
    /* copy up atoms into send buffer */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_up, i, j, 1, uvec );

    /* send up, receive down */
    sendrecv_buf(&send_buf_up, nbup, &recv_buf_down, nbdown, &stat);

    /* upack atoms from down */
    recv_buf_down.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*unpack_func)( &recv_buf_down, i, j, cell_dim.z-1 );

    /* copy down atoms into send buffer */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_down, i, j, cell_dim.z-2, dvec );

    /* send down, receive up */
    sendrecv_buf(&send_buf_down, nbdown, &recv_buf_up, nbup, &stat);
 
    /* unpack atoms from up */
    recv_buf_up.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*unpack_func)( &recv_buf_up, i, j, 0 );
  }
#endif

  /* exchange north/south */
  if (cpu_dim.y==1) {
    /* simply copy north/south atoms to buffer cells */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j) {
        (*copy_func)( i, 1, j, i, cell_dim.y-1, j, nvec );
        (*copy_func)( i, cell_dim.y-2, j, i, 0, j, svec );
      }
  }
#ifdef MPI
  else {
    /* copy north atoms into send buffer */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_north, i, 1, j, nvec );

    /* send north, receive south */
    sendrecv_buf(&send_buf_north, nbnorth, &recv_buf_south, nbsouth, &stat);

    /* unpack atoms from south */
    recv_buf_south.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_south, i, cell_dim.y-1, j );

    /* copy south atoms into send buffer */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_south, i, cell_dim.y-2, j, svec );

    /* send south, receive north */
    sendrecv_buf(&send_buf_south, nbsouth, &recv_buf_north, nbnorth, &stat);

    /* unpack atoms from north */
    recv_buf_north.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_north, i, 0, j );
  }
#endif

  /* exchange east/west */
  if (cpu_dim.x==1) {
    /* simply copy east/west atoms to buffer cells */
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) { 
        (*copy_func)( 1, i, j, cell_dim.x-1, i, j, evec );
#if !defined(AR) || defined(COVALENT) 
        (*copy_func)( cell_dim.x-2, i, j, 0, i, j, wvec );
#endif
      }
  }
#ifdef MPI
  else {
    /* copy east atoms into send buffer */
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) 
        (*pack_func)( &send_buf_east, 1, i, j, evec );

    /* send east, receive west */
    sendrecv_buf(&send_buf_east, nbeast, &recv_buf_west, nbwest, &stat);

    /* unpack atoms from west */
    recv_buf_west.n = 0;
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) 
        (*unpack_func)( &recv_buf_west, cell_dim.x-1, i, j );

#if !defined(AR) || defined(COVALENT) 
    /* copy west atoms into send buffer */
    for (i=0; i < cell_dim.y; ++i) 
      for (j=0; j < cell_dim.z; ++j) 
        (*pack_func)( &send_buf_west, cell_dim.x-2, i, j, wvec );

    /* send west, receive east */
    sendrecv_buf(&send_buf_west, nbwest, &recv_buf_east, nbeast, &stat);

    /* unpack atoms from east */
    recv_buf_east.n = 0;
    for (i=0; i < cell_dim.y; ++i) 
      for (j=0; j < cell_dim.z; ++j) 
        (*unpack_func)( &recv_buf_east, 0, i, j );
#endif
  }  
#endif
}

#else /* not SR */

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

void send_cells(void (*copy_func)  (int, int, int, int, int, int, vektor),
                void (*pack_func)  (msgbuf*, int, int, int, vektor),
                void (*unpack_func)(msgbuf*, int, int, int))
{
  int i,j;

  vektor uvec={0,0,0}, dvec={0,0,0};
  vektor nvec={0,0,0}, svec={0,0,0};
  vektor evec={0,0,0}, wvec={0,0,0};

#ifdef MPI
  MPI_Status  stateast[2],  statwest[2];
  MPI_Status statnorth[2], statsouth[2];
  MPI_Status    statup[2],  statdown[2];

  MPI_Request  reqeast[2],   reqwest[2];
  MPI_Request reqnorth[2],  reqsouth[2];
  MPI_Request    requp[2],   reqdown[2];

  empty_mpi_buffers();
#endif

#ifdef VEC
  atoms.n_buf = atoms.n;
#endif
#ifdef NBLIST
  if (pbc_dirs.x==1) {
    if (my_coord.x==0) evec = box_x;
    if (my_coord.x==cpu_dim.x-1) {
      wvec.x = -box_x.x; wvec.y = -box_x.y; wvec.z = -box_x.z;
    }    
  }
  if (pbc_dirs.y==1) {
    if (my_coord.y==0) nvec = box_y;
    if (my_coord.y==cpu_dim.y-1) {
      svec.x = -box_y.x; svec.y = -box_y.y; svec.z = -box_y.z;
    }    
  }
  if (pbc_dirs.z==1) {
    if (my_coord.z==0) uvec = box_z;
    if (my_coord.z==cpu_dim.z-1) {
      dvec.x = -box_z.x; dvec.y = -box_z.y; dvec.z = -box_z.z;
    }    
  }
#endif

  /* exchange up/down */
  if (cpu_dim.z==1) {
    /* simply copy up/down atoms to buffer cells*/
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j) {
        (*copy_func)( i, j, 1, i, j, cell_dim.z-1, uvec );
        (*copy_func)( i, j, cell_dim.z-2, i, j, 0, dvec );
      }
  }
#ifdef MPI
  else {
    /* copy up atoms into send buffer, send up */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_up, i, j, 1, uvec );
    irecv_buf( &recv_buf_down , nbdown, &reqdown[1]);
    isend_buf( &send_buf_up   , nbup  , &reqdown[0]);

    /* copy down atoms into send buffer, send down */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_down, i, j, cell_dim.z-2, dvec );
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
#endif

  /* exchange north/south */
  if (cpu_dim.y==1) {
    /* simply copy north/south atoms to buffer cells */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j) {
        (*copy_func)( i, 1, j, i, cell_dim.y-1, j, nvec );
        (*copy_func)( i, cell_dim.y-2, j, i, 0, j, svec );
      }
  }
#ifdef MPI
  else {
    /* copy north atoms into send buffer, send north */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_north, i, 1, j, nvec );
    irecv_buf( &recv_buf_south, nbsouth, &reqsouth[1] );
    isend_buf( &send_buf_north, nbnorth, &reqsouth[0] );

    /* copy south atoms into send buffer, send south*/
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_south, i, cell_dim.y-2, j, svec );
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
#endif

  /* exchange east/west*/
  if (cpu_dim.x==1) {
  /* simply copy east/west atoms to buffer cells*/
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) { 
        (*copy_func)( 1, i, j, cell_dim.x-1, i, j, evec );
#if !defined(AR) || defined(COVALENT) 
        (*copy_func)( cell_dim.x-2, i, j, 0, i, j, wvec );
#endif
      }
  }
#ifdef MPI
  else {
  /* copy east atoms into send buffer, send east*/
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) 
        (*pack_func)( &send_buf_east, 1, i, j, evec );
    irecv_buf( &recv_buf_west, nbwest, &reqwest[1] );
    isend_buf( &send_buf_east, nbeast, &reqwest[0] );

#if !defined(AR) || defined(COVALENT) 
    /* copy west atoms into send buffer, send west*/
    for (i=0; i < cell_dim.y; ++i) 
      for (j=0; j < cell_dim.z; ++j) 
        (*pack_func)( &send_buf_west, cell_dim.x-2, i, j, wvec );
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
#endif
}

#endif /* not SR */

#ifdef SR

/******************************************************************************
*
*  Add forces in buffer cells on neighboring CPUs back to the original cells.
*  What exactly is sent is determined by the parameter functions.
*  We use Steve Plimptons communication scheme: we send only along
*  the main axis of the system, so that edge cells travel twice,
*  and corner cells three times. If not COVALENT, one buffer cell wall 
*  (including adjacent edge and corner cells) contains no forces.
*
******************************************************************************/

void send_forces(void (*add_func)   (int, int, int, int, int, int),
                 void (*pack_func)  (msgbuf*, int, int, int),
                 void (*unpack_func)(msgbuf*, int, int, int))
{
  int i,j;

#ifdef MPI
  MPI_Status  stat;
  empty_mpi_buffers();
#endif

  /* send forces east/west */
  if (cpu_dim.x==1) {
    /* simply add east/west forces to original cells */
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) { 
#ifdef COVALENT 
        (*add_func)( 0, i, j, cell_dim.x-2, i, j );
#endif
        (*add_func)( cell_dim.x-1, i, j, 1, i, j );
      }
  }
#ifdef MPI
  else {
#ifdef COVALENT 
    /* copy east forces into send buffer */
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_east, 0, i, j );

    /* send east, receive west */
    sendrecv_buf(&send_buf_east, nbeast, &recv_buf_west, nbwest, &stat);

    /* add forces from west to original cells */
    recv_buf_west.n = 0;
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) 
        (*unpack_func)( &recv_buf_west, cell_dim.x-2, i, j );
#endif

    /* copy west forces into send buffer */
    for (i=0; i < cell_dim.y; ++i) 
      for (j=0; j < cell_dim.z; ++j) 
        (*pack_func)( &send_buf_west, cell_dim.x-1, i, j );

    /* send west, receive east */
    sendrecv_buf(&send_buf_west, nbwest, &recv_buf_east, nbeast, &stat);

    /* add forces from east to original cells */
    recv_buf_east.n = 0;
    for (i=0; i < cell_dim.y; ++i) 
      for (j=0; j < cell_dim.z; ++j) 
        (*unpack_func)( &recv_buf_east, 1, i, j );
  }
#endif

  /* send forces north/south */
  if (cpu_dim.y==1) {
    /* simply add north/south forces to original cells */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j) {
        (*add_func)( i, 0, j, i, cell_dim.y-2, j );
        (*add_func)( i, cell_dim.y-1, j, i, 1, j );
      }
  }
#ifdef MPI
  else {
    /* copy north forces into send buffer */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_north, i, 0, j );

    /* send north, receive south */
    sendrecv_buf(&send_buf_north, nbnorth, &recv_buf_south, nbsouth, &stat);

    /* add forces from south to original cells */
    recv_buf_south.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_south, i, cell_dim.y-2, j );

    /* copy south forces into send buffer */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*pack_func)( &send_buf_south, i, cell_dim.y-1, j );

    /* send south, receive north */
    sendrecv_buf(&send_buf_south, nbsouth, &recv_buf_north, nbnorth, &stat);

    /* add forces from north to original cells */
    recv_buf_north.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=0; j < cell_dim.z; ++j)
        (*unpack_func)( &recv_buf_north, i, 1, j );
  }
#endif

  /* send forces up/down */
  if (cpu_dim.z==1) {
    /* simply add up/down forces to original cells */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j) {
        (*add_func)( i, j, 0, i, j, cell_dim.z-2 );
        (*add_func)( i, j, cell_dim.z-1, i, j, 1 );
      }
  }
#ifdef MPI
  else {
    /* copy up forces into send buffer */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_up, i, j, 0 );

    /* send up, receive down */
    sendrecv_buf(&send_buf_up, nbup, &recv_buf_down, nbdown, &stat);
  
    /* add forces from down to original cells */
    recv_buf_down.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*unpack_func)( &recv_buf_down, i, j, cell_dim.z-2 );

    /* copy down forces into send buffer */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*pack_func)( &send_buf_down, i, j, cell_dim.z-1 );

    /* send down, receive up */
    sendrecv_buf(&send_buf_down, nbdown, &recv_buf_up, nbup, &stat);

    /* add forces from up to original cells */
    recv_buf_up.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j)
        (*unpack_func)( &recv_buf_up, i, j, 1 );
  }
#endif
}

#else /* not SR */

/******************************************************************************
*
*  Add forces in buffer cells on neighboring CPUs back to the original cells.
*  What exactly is sent is determined by the parameter functions.
*  We use Steve Plimptons communication scheme: we send only along
*  the main axis of the system, so that edge cells travel twice,
*  and corner cells three times. If not COVALENT, one buffer cell wall 
*  (including adjacent edge and corner cells) contains no forces.
*
******************************************************************************/

void send_forces(void (*add_func)   (int, int, int, int, int, int),
                 void (*pack_func)  (msgbuf*, int, int, int),
                 void (*unpack_func)(msgbuf*, int, int, int))
{
  int i,j;

#ifdef MPI
  MPI_Status  stateast[2],  statwest[2];
  MPI_Status statnorth[2], statsouth[2];
  MPI_Status    statup[2],  statdown[2];

  MPI_Request  reqeast[2],   reqwest[2];
  MPI_Request reqnorth[2],  reqsouth[2];
  MPI_Request    requp[2],   reqdown[2];

  empty_mpi_buffers();
#endif

  /* send forces east/west */
  if (cpu_dim.x==1) {
    /* simply add east/west forces to original cells */
    for (i=0; i < cell_dim.y; ++i)
      for (j=0; j < cell_dim.z; ++j) { 
#ifdef COVALENT 
        (*add_func)( 0, i, j, cell_dim.x-2, i, j );
#endif
        (*add_func)( cell_dim.x-1, i, j, 1, i, j );
      }
  }
#ifdef MPI
  else {
#ifdef COVALENT 
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

#ifdef COVALENT 
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
#endif

  /* send forces north/south */
  if (cpu_dim.y==1) {
    /* simply add north/south forces to original cells */
    for (i=1; i < cell_dim.x-1; ++i)
      for (j=0; j < cell_dim.z; ++j) {
        (*add_func)( i, 0, j, i, cell_dim.y-2, j );
        (*add_func)( i, cell_dim.y-1, j, i, 1, j );
      }
  }
#ifdef MPI
  else {
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
#endif

  /* send forces up/down */
  if (cpu_dim.z==1) {
    /* simply add up/down forces to original cells */
    for (i=1; i < cell_dim.x-1; ++i) 
      for (j=1; j < cell_dim.y-1; ++j) {
        (*add_func)( i, j, 0, i, j, cell_dim.z-2 );
        (*add_func)( i, j, cell_dim.z-1, i, j, 1 );
      }
  }
#ifdef MPI
  else {
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
#endif
}

#endif /* not SR */

/******************************************************************************
*
*  copy contents of one cell to another (buffer) cell (for force comp.)
*
******************************************************************************/

void copy_cell( int k, int l, int m, int r, int s, int t, vektor v )
{
  int i, count;
  minicell *from, *to;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  to   = PTR_3D_V(cell_array, r, s, t, cell_dim);

  /* increase minicell size if necessary */
  if (from->n > to->n_max) {
    to->n = 0;
    ALLOC_MINICELL(to, from->n);
  }

#ifdef VEC
  /* increase cell size if necessary */
  if (atoms.n_buf + from->n > atoms.n_max) 
    alloc_cell( &atoms, atoms.n_buf + from->n); 
#endif

  to->n = from->n;
#ifdef VEC
  count = atoms.n_buf;
#pragma cdir nodep
#endif
  for (i=0; i<to->n; ++i) {
#ifdef VEC
    to->ind[i]    = count++;
#endif
    ORT(to,i,X)   = ORT(from,i,X) + v.x;
    ORT(to,i,Y)   = ORT(from,i,Y) + v.y;
    ORT(to,i,Z)   = ORT(from,i,Z) + v.z;
#ifndef MONO
    SORTE(to,i)   = SORTE(from,i);
#endif
#ifdef UNIAX
    ACHSE(to,i,X) = ACHSE(from,i,X);
    ACHSE(to,i,Y) = ACHSE(from,i,Y);
    ACHSE(to,i,Z) = ACHSE(from,i,Z);
#endif
  }
#ifdef VEC
  atoms.n_buf = count;
#endif
}

/******************************************************************************
*
*  pack cell into MPI send buffer (for force comp.)
*
******************************************************************************/

void pack_cell( msgbuf *b, int k, int l, int m, vektor v )
{
  int i, j = b->n;
  minicell *from;
    
  from = PTR_3D_V(cell_array, k, l, m, cell_dim);

  b->data[ j++ ] = (real) from->n;

#ifdef VEC
#pragma cdir nodep
#endif    
  for (i=0; i<from->n; ++i) {
    b->data[ j++ ] = ORT(from,i,X) + v.x;
    b->data[ j++ ] = ORT(from,i,Y) + v.y;
    b->data[ j++ ] = ORT(from,i,Z) + v.z;
#ifndef MONO
    b->data[ j++ ] = (real) SORTE(from,i);
#endif
#ifdef UNIAX
    b->data[ j++ ] = ACHSE(from,i,X);
    b->data[ j++ ] = ACHSE(from,i,Y);
    b->data[ j++ ] = ACHSE(from,i,Z);
#endif
  }
  if (b->n_max < b->n)  error("Buffer overflow in pack_cell");
  b->n = j;
}

/******************************************************************************
*
*  unpack cell from MPI buffer to buffer cell (for force comp.)
*
******************************************************************************/

void unpack_cell( msgbuf *b, int k, int l, int m )
{
  int i, j = b->n, count, tmp_n;
  minicell *to;

  to = PTR_3D_V(cell_array, k, l, m, cell_dim);

  tmp_n = (int) b->data[ j++ ];
  
  /* increase minicell size if necessary */
  if (tmp_n > to->n_max) {
    to->n = 0;
    ALLOC_MINICELL(to, tmp_n);
  }

#ifdef VEC
  /* increase cell size if necessary */
  if (atoms.n_buf + tmp_n > atoms.n_max) 
    alloc_cell( &atoms, atoms.n_buf + tmp_n); 
#endif

  /* copy indices and atoms */
  to->n = tmp_n;
#ifdef VEC
  count = atoms.n_buf;
#pragma cdir nodep
#endif
  for (i=0; i<to->n; ++i) {
#ifdef VEC
    to->ind[i]  = count++;
#endif
    ORT(to,i,X) = b->data[ j++ ];
    ORT(to,i,Y) = b->data[ j++ ];
    ORT(to,i,Z) = b->data[ j++ ];
#ifndef MONO
    SORTE(to,i) = (shortint) b->data[ j++ ];
#endif
#ifdef UNIAX
    ACHSE(to,i,X) = b->data[ j++ ];
    ACHSE(to,i,Y) = b->data[ j++ ];
    ACHSE(to,i,Z) = b->data[ j++ ];
#endif
  }
  if (b->n_max < b->n) error("Buffer overflow in unpack_cell");
#ifdef VEC
  atoms.n_buf = count;
#endif
  b->n = j;
}

/******************************************************************************
*
*  add forces of one cell to those of another cell
*
******************************************************************************/

void add_forces( int k, int l, int m, int r, int s, int t )
{
  int i;
  minicell *from, *to;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  to   = PTR_3D_V(cell_array, r, s, t, cell_dim);

#ifdef VEC
#pragma cdir nodep
#endif
  for (i=0; i<to->n; ++i) {
    KRAFT(to,i,X)  += KRAFT(from,i,X);
    KRAFT(to,i,Y)  += KRAFT(from,i,Y);
    KRAFT(to,i,Z)  += KRAFT(from,i,Z);
#ifndef MONOLJ
    POTENG(to,i)   += POTENG(from,i);
#endif
#ifdef NVX
    HEATCOND(to,i) += HEATCOND(from,i);
#endif
#ifdef STRESS_TENS
    PRESSTENS(to,i,xx) += PRESSTENS(from,i,xx);
    PRESSTENS(to,i,yy) += PRESSTENS(from,i,yy);
    PRESSTENS(to,i,zz) += PRESSTENS(from,i,zz);
    PRESSTENS(to,i,yz) += PRESSTENS(from,i,yz);
    PRESSTENS(to,i,zx) += PRESSTENS(from,i,zx);
    PRESSTENS(to,i,xy) += PRESSTENS(from,i,xy);
#endif
#ifdef NNBR
    NBANZ(to,i) += NBANZ(from,i);
#endif
#ifdef UNIAX
    DREH_MOMENT(to,i,X) += DREH_MOMENT(from,i,X);
    DREH_MOMENT(to,i,Y) += DREH_MOMENT(from,i,Y);
    DREH_MOMENT(to,i,Z) += DREH_MOMENT(from,i,Z);
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
  int i, j = b->n;
  minicell *from;
    
  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
#ifdef VEC
#pragma cdir nodep
#endif
  for (i=0; i<from->n; ++i) {
    b->data[ j++ ] = KRAFT(from,i,X);
    b->data[ j++ ] = KRAFT(from,i,Y);
    b->data[ j++ ] = KRAFT(from,i,Z);
#ifndef MONOLJ
    b->data[ j++ ] = POTENG(from,i);
#endif
#ifdef NVX
    b->data[ j++ ] = HEATCOND(from,i);
#endif
#ifdef STRESS_TENS
    b->data[ j++ ] = PRESSTENS(from,i,xx);
    b->data[ j++ ] = PRESSTENS(from,i,yy);
    b->data[ j++ ] = PRESSTENS(from,i,zz);
    b->data[ j++ ] = PRESSTENS(from,i,yz);
    b->data[ j++ ] = PRESSTENS(from,i,zx);
    b->data[ j++ ] = PRESSTENS(from,i,xy);
#endif
#ifdef NNBR
    b->data[ j++ ] = (real) NBANZ(from,i);
#endif
#ifdef UNIAX
    b->data[ j++ ] = DREH_MOMENT(from,i,X);
    b->data[ j++ ] = DREH_MOMENT(from,i,Y);
    b->data[ j++ ] = DREH_MOMENT(from,i,Z);
#endif
  }
  if (b->n_max < b->n) error("Buffer overflow in pack_forces.");
  b->n = j;
}

/******************************************************************************
*
*  unpack forces from MPI buffer, and add them to those of the original cell
*
******************************************************************************/

void unpack_forces( msgbuf *b, int k, int l, int m )
 {
  int i, j = b->n;
  minicell *to;

  to = PTR_3D_V(cell_array, k, l, m, cell_dim);
#ifdef VEC
#pragma cdir nodep
#endif
  for (i=0; i<to->n; ++i) {
    KRAFT(to,i,X)  += b->data[ j++ ];
    KRAFT(to,i,Y)  += b->data[ j++ ];
    KRAFT(to,i,Z)  += b->data[ j++ ];
#ifndef MONOLJ
    POTENG(to,i)   += b->data[ j++ ];
#endif
#ifdef NVX
    HEATCOND(to,i) += b->data[ j++ ];
#endif
#ifdef STRESS_TENS
    PRESSTENS(to,i,xx) += b->data[ j++ ];
    PRESSTENS(to,i,yy) += b->data[ j++ ];
    PRESSTENS(to,i,zz) += b->data[ j++ ];
    PRESSTENS(to,i,yz) += b->data[ j++ ];
    PRESSTENS(to,i,zx) += b->data[ j++ ];
    PRESSTENS(to,i,xy) += b->data[ j++ ];
#endif
#ifdef NNBR
    NBANZ(to,i) += (shortint) b->data[ j++ ];
#endif
#ifdef UNIAX
    DREH_MOMENT(to,i,X) += b->data[ j++ ];
    DREH_MOMENT(to,i,Y) += b->data[ j++ ];
    DREH_MOMENT(to,i,Z) += b->data[ j++ ];
#endif
  }
  if (b->n_max < b->n) error("Buffer overflow in unpack_forces.");
  b->n = j;
}


#ifdef EAM2

/******************************************************************************
*
*  copy eam_dF and eeam_dM of one cell to another cell
*
******************************************************************************/

void copy_dF( int k, int l, int m, int r, int s, int t, vektor v )
{
  int i;
  minicell *from, *to;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  to   = PTR_3D_V(cell_array, r, s, t, cell_dim);

#ifdef VEC
#pragma cdir nodep
#endif
  for (i=0; i<to->n; ++i) {
    EAM_DF    (to,i)    = EAM_DF    (from,i);
#ifdef EEAM
    EAM_DM    (to,i)    = EAM_DM    (from,i);
#endif
#ifdef ADP
    ADP_MU    (to,i,X)  = ADP_MU    (from,i,X);
    ADP_MU    (to,i,Y)  = ADP_MU    (from,i,Y);
    ADP_MU    (to,i,Z)  = ADP_MU    (from,i,Z);
    ADP_LAMBDA(to,i,xx) = ADP_LAMBDA(from,i,xx);
    ADP_LAMBDA(to,i,yy) = ADP_LAMBDA(from,i,yy);
    ADP_LAMBDA(to,i,zz) = ADP_LAMBDA(from,i,zz);
    ADP_LAMBDA(to,i,yz) = ADP_LAMBDA(from,i,yz);
    ADP_LAMBDA(to,i,zx) = ADP_LAMBDA(from,i,zx);
    ADP_LAMBDA(to,i,xy) = ADP_LAMBDA(from,i,xy);
#endif
  }
}

/******************************************************************************
*
*  add eam_rho of one cell to another cell
*  add eeam_p_h of one cell to another cell
*
******************************************************************************/

void add_rho( int k, int l, int m, int r, int s, int t )
{
  int i;
  minicell *from, *to;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  to   = PTR_3D_V(cell_array, r, s, t, cell_dim);

#ifdef VEC
#pragma cdir nodep
#endif
  for (i=0; i<to->n; ++i) {
    EAM_RHO(to,i) += EAM_RHO(from,i);
#ifdef EEAM
    EAM_P(to,i)   += EAM_P(from,i);
#endif
#ifdef ADP
    ADP_MU    (to,i,X)  += ADP_MU    (from,i,X);
    ADP_MU    (to,i,Y)  += ADP_MU    (from,i,Y);
    ADP_MU    (to,i,Z)  += ADP_MU    (from,i,Z);
    ADP_LAMBDA(to,i,xx) += ADP_LAMBDA(from,i,xx);
    ADP_LAMBDA(to,i,yy) += ADP_LAMBDA(from,i,yy);
    ADP_LAMBDA(to,i,zz) += ADP_LAMBDA(from,i,zz);
    ADP_LAMBDA(to,i,yz) += ADP_LAMBDA(from,i,yz);
    ADP_LAMBDA(to,i,zx) += ADP_LAMBDA(from,i,zx);
    ADP_LAMBDA(to,i,xy) += ADP_LAMBDA(from,i,xy);
#endif
  }
}

/******************************************************************************
*
*  pack eam_dF and eeam_dM into MPI buffer
*
******************************************************************************/

void pack_dF( msgbuf *b, int k, int l, int m, vektor v )
{
  int i, j = b->n;
  minicell *from;
    
  from = PTR_3D_V(cell_array, k, l, m, cell_dim);

#ifdef VEC
#pragma cdir nodep
#endif
  for (i=0; i<from->n; ++i) {
    b->data[ j++ ] = EAM_DF    (from,i);
#ifdef EEAM
    b->data[ j++ ] = EAM_DM    (from,i);
#endif
#ifdef ADP
    b->data[ j++ ] = ADP_MU    (from,i,X);
    b->data[ j++ ] = ADP_MU    (from,i,X);
    b->data[ j++ ] = ADP_MU    (from,i,X);
    b->data[ j++ ] = ADP_LAMBDA(from,i,xx);
    b->data[ j++ ] = ADP_LAMBDA(from,i,yy);
    b->data[ j++ ] = ADP_LAMBDA(from,i,zz);
    b->data[ j++ ] = ADP_LAMBDA(from,i,yz);
    b->data[ j++ ] = ADP_LAMBDA(from,i,zx);
    b->data[ j++ ] = ADP_LAMBDA(from,i,xy);
#endif
  }
  b->n = j;
}

/******************************************************************************
*
*  pack eam_rho into MPI buffer
*  pack eeam_p_h into MPI buffer
*
******************************************************************************/

void pack_rho( msgbuf *b, int k, int l, int m )
{
  int i, j = b->n;
  minicell *from;
    
  from = PTR_3D_V(cell_array, k, l, m, cell_dim);

#ifdef VEC
#pragma cdir nodep
#endif
  for (i=0; i<from->n; ++i) {
    b->data[ j++ ] = EAM_RHO(from,i);
#ifdef EEAM
    b->data[ j++ ] = EAM_P(from,i);
#endif
#ifdef ADP
    b->data[ j++ ] = ADP_MU    (from,i,X);
    b->data[ j++ ] = ADP_MU    (from,i,X);
    b->data[ j++ ] = ADP_MU    (from,i,X);
    b->data[ j++ ] = ADP_LAMBDA(from,i,xx);
    b->data[ j++ ] = ADP_LAMBDA(from,i,yy);
    b->data[ j++ ] = ADP_LAMBDA(from,i,zz);
    b->data[ j++ ] = ADP_LAMBDA(from,i,yz);
    b->data[ j++ ] = ADP_LAMBDA(from,i,zx);
    b->data[ j++ ] = ADP_LAMBDA(from,i,xy);
#endif
  }
  b->n = j;
}

/******************************************************************************
*
*  unpack eam_dF and eeam_dM from MPI buffer into cell
*
******************************************************************************/

void unpack_dF( msgbuf *b, int k, int l, int m )
{
  int i, j = b->n;
  minicell *to;

  to = PTR_3D_V(cell_array, k, l, m, cell_dim);

#ifdef VEC
#pragma cdir nodep
#endif
  for (i=0; i<to->n; ++i) {
    EAM_DF    (to,i)    = b->data[ j++ ];
#ifdef EEAM
    EAM_DM    (to,i)    = b->data[ j++ ];
#endif
#ifdef ADP
    ADP_MU    (to,i,X)  = b->data[ j++ ];
    ADP_MU    (to,i,Y)  = b->data[ j++ ];
    ADP_MU    (to,i,Z)  = b->data[ j++ ];
    ADP_LAMBDA(to,i,xx) = b->data[ j++ ];
    ADP_LAMBDA(to,i,yy) = b->data[ j++ ];
    ADP_LAMBDA(to,i,zz) = b->data[ j++ ];
    ADP_LAMBDA(to,i,yz) = b->data[ j++ ];
    ADP_LAMBDA(to,i,zx) = b->data[ j++ ];
    ADP_LAMBDA(to,i,xy) = b->data[ j++ ];
#endif
  }
  b->n = j;
}

/******************************************************************************
*
*  unpack and add eam_rho from MPI buffer into cell
*  unpack and add eeam_p_h from MPI buffer into cell
*
******************************************************************************/

void unpack_add_rho( msgbuf *b, int k, int l, int m )
{
  int i, j = b->n;
  minicell *to;

  to = PTR_3D_V(cell_array, k, l, m, cell_dim);

#ifdef VEC
#pragma cdir nodep
#endif
  for (i=0; i<to->n; ++i) {
    EAM_RHO(to,i) += b->data[ j++ ];
#ifdef EEAM
    EAM_P(to,i)   += b->data[ j++ ];
#endif
#ifdef ADP
    ADP_MU    (to,i,X)  += b->data[ j++ ];
    ADP_MU    (to,i,Y)  += b->data[ j++ ];
    ADP_MU    (to,i,Z)  += b->data[ j++ ];
    ADP_LAMBDA(to,i,xx) += b->data[ j++ ];
    ADP_LAMBDA(to,i,yy) += b->data[ j++ ];
    ADP_LAMBDA(to,i,zz) += b->data[ j++ ];
    ADP_LAMBDA(to,i,yz) += b->data[ j++ ];
    ADP_LAMBDA(to,i,zx) += b->data[ j++ ];
    ADP_LAMBDA(to,i,xy) += b->data[ j++ ];
#endif
  }
  b->n = j;
}

#endif /* EAM2 */
