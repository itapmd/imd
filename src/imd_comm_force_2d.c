/******************************************************************************
*
* imd_comm_force_2d.c -- communication for force computation, two dimensions
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"


/******************************************************************************
*
*  Send cells into buffer cells on neighboring CPUs for the force computation.
*  What exactly is sent is determined by the parameter functions.
*  We use Steve Plimptons communication scheme: we send only along
*  the main axis of the system, so that corner cells travel twice.
*  In AR mode, one cell wall (including adjacent corner cells) is not 
*  needed in the buffer cells.
*
******************************************************************************/

void send_cells(void (*copy_func)  (int, int, int, int),
                void (*pack_func)  (msgbuf*, int, int),
                void (*unpack_func)(msgbuf*, int, int))
{
  int i;

  MPI_Status  stateast[2],  statwest[2];
  MPI_Status statnorth[2], statsouth[2];

  MPI_Request  reqeast[2],   reqwest[2];
  MPI_Request reqnorth[2],  reqsouth[2];

  empty_mpi_buffers();

  /* exchange north/south */
  if (cpu_dim.y==1) {
    /* simply copy north/south atoms to buffer cells */
    for (i=1; i<cell_dim.x-1; ++i) {
      (*copy_func)( i, 1, i, cell_dim.y-1 );
      (*copy_func)( i, cell_dim.y-2, i, 0 );
    }
  } else {
    /* copy north atoms into send buffer, send north */
    for (i=1; i<cell_dim.x-1; ++i)
      (*pack_func)( &send_buf_north, i, 1 );
    isend_buf( &send_buf_north, nbnorth, &reqsouth[0] );
    irecv_buf( &recv_buf_south, nbsouth, &reqsouth[1] );

    /* copy south atoms into send buffer, send south */
    for (i=1; i<cell_dim.x-1; ++i)
      (*pack_func)( &send_buf_south, i, cell_dim.y-2 );
    isend_buf( &send_buf_south, nbsouth, &reqnorth[0] );
    irecv_buf( &recv_buf_north, nbnorth, &reqnorth[1] );

    /* wait for atoms from south, move them to buffer cells */
    MPI_Waitall(2, reqsouth, statsouth);
    recv_buf_south.n = 0;
    for (i=1; i<cell_dim.x-1; ++i)
      (*unpack_func)( &recv_buf_south, i, cell_dim.y-1 );

    /* Wait for atoms from north, move them to buffer cells */
    MPI_Waitall(2, reqnorth, statnorth);
    recv_buf_north.n = 0;
    for (i=1; i<cell_dim.x-1; ++i)
      (*unpack_func)( &recv_buf_north, i, 0 );
  }

  /* exchange east/west */
  if (cpu_dim.x==1) {
    /* simply copy east/west atoms to buffer cells */
    for (i=0; i<cell_dim.y; ++i) {
      (*copy_func)( 1, i, cell_dim.x-1, i );
#ifndef AR
      (*copy_func)( cell_dim.x-2, i, 0, i );
#endif
    }
  } else {
    /* copy east atoms into send buffer, send east */
    for (i=0; i<cell_dim.y; ++i)
      (*pack_func)( &send_buf_east, 1, i );
    isend_buf( &send_buf_east, nbeast, &reqwest[0] );
    irecv_buf( &recv_buf_west, nbwest, &reqwest[1] );

#ifndef AR
    /* copy west atoms into send buffer, send west */
    for (i=0; i<cell_dim.y; ++i) 
      (*pack_func)( &send_buf_west, cell_dim.x-2, i );
    /* send west */
    isend_buf( &send_buf_west, nbwest, &reqeast[0] );
    irecv_buf( &recv_buf_east, nbeast, &reqeast[1] );
#endif

    /* Wait for atoms from west, move them to buffer cells */
    MPI_Waitall(2, reqwest, statwest);
    recv_buf_west.n = 0;
    for (i=0; i<cell_dim.y; ++i)
      (*unpack_func)( &recv_buf_west, cell_dim.x-1, i );

#ifndef AR
    /* Wait for atoms from east, move them to buffer cells */
    MPI_Waitall(2, reqeast, stateast);
    recv_buf_east.n = 0;
    for (i=0; i<cell_dim.y; ++i) 
      (*unpack_func)( &recv_buf_east, 0, i );
#endif
  }
}


/******************************************************************************
*
*  Add forces in buffer cells on neighboring CPUs back to the original cells.
*  What exactly is sent is determined by the parameter functions.
*  We use Steve Plimptons communication scheme: we send only along
*  the main axis of the system, so that corner cells travel twice.
*  In AR mode, one buffer cell wall (including adjacent corner cells) 
*  contains no forces.
*
******************************************************************************/

void send_forces(void (*add_func)   (int, int, int, int),
                 void (*pack_func)  (msgbuf*, int, int),
                 void (*unpack_func)(msgbuf*, int, int))
{
  int i;

  MPI_Status  stateast[2],  statwest[2];
  MPI_Status statnorth[2], statsouth[2];

  MPI_Request  reqeast[2],   reqwest[2];
  MPI_Request reqnorth[2],  reqsouth[2];

  empty_mpi_buffers();

  /* send forces east/west */
  if (cpu_dim.x==1) {
    /* simply add east/west forces to original cells */
    for (i=0; i < cell_dim.y; ++i) {
#ifndef AR 
      (*add_func)( 0, i, cell_dim.x-2, i );
#endif
      (*add_func)( cell_dim.x-1, i, 1, i );
    }
  } else {
#ifndef AR 
    /* copy east forces into send buffer, send east */
    for (i=0; i < cell_dim.y; ++i)
      (*pack_func)( &send_buf_east, 0, i );
    irecv_buf( &recv_buf_west, nbwest, &reqwest[1] );
    isend_buf( &send_buf_east, nbeast, &reqwest[0] );
#endif

    /* copy west forces into send buffer, send west */
    for (i=0; i < cell_dim.y; ++i) 
      (*pack_func)( &send_buf_west, cell_dim.x-1, i );
    irecv_buf( &recv_buf_east, nbeast, &reqeast[1] );
    isend_buf( &send_buf_west, nbwest, &reqeast[0] );

#ifndef AR
    /* wait for forces from west, add them to original cells */
    MPI_Waitall(2, reqwest, statwest);
    recv_buf_west.n = 0;
    for (i=0; i < cell_dim.y; ++i)
      (*unpack_func)( &recv_buf_west, cell_dim.x-2, i );
#endif

    /* wait for forces from east, add them to original cells */
    MPI_Waitall(2, reqeast, stateast);
    recv_buf_east.n = 0;
    for (i=0; i < cell_dim.y; ++i) 
      (*unpack_func)( &recv_buf_east, 1, i );
  }

  /* send forces north/south */
  if (cpu_dim.y==1) {
    /* simply add north/south forces to original cells */
    for (i=1; i < cell_dim.x-1; ++i) {
      (*add_func)( i, 0, i, cell_dim.y-2 );
      (*add_func)( i, cell_dim.y-1, i, 1 );
    }
  } else {
    /* copy north forces into send buffer, send north */
    for (i=1; i < cell_dim.x-1; ++i) 
      (*pack_func)( &send_buf_north, i, 0 );
    irecv_buf( &recv_buf_south, nbsouth, &reqsouth[1] );
    isend_buf( &send_buf_north, nbnorth, &reqsouth[0] );

    /* copy south forces into send buffer, send south */
    for (i=1; i < cell_dim.x-1; ++i) 
      (*pack_func)( &send_buf_south, i, cell_dim.y-1 );
    irecv_buf( &recv_buf_north, nbnorth, &reqnorth[1] );
    isend_buf( &send_buf_south, nbsouth, &reqnorth[0] );

    /* wait for forces from south, add them to original cells */
    MPI_Waitall(2, reqsouth, statsouth);
    recv_buf_south.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      (*unpack_func)( &recv_buf_south, i, cell_dim.y-2 );

    /* wait for forces from north, add them to original cells */
    MPI_Waitall(2, reqnorth, statnorth);
    recv_buf_north.n = 0;
    for (i=1; i < cell_dim.x-1; ++i) 
      (*unpack_func)( &recv_buf_north, i, 1 );
  }
}


/******************************************************************************
*
*  copy contents of one cell to another (buffer) cell (for force comp.)
*
******************************************************************************/

void copy_cell( int j, int k, int l, int m )
{
  int i, tmp_n;
  cell *from, *to;

  from = PTR_2D_V(cell_array, j, k, cell_dim);
  to   = PTR_2D_V(cell_array, l, m, cell_dim);

  tmp_n = from->n;
  if (tmp_n > to->n_max) {
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
*  pack cell into MPI send buffer (for force comp.)
*
******************************************************************************/

void pack_cell( msgbuf *b, int j, int k )
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
  if (b->n_max < b->n)  error("Buffer overflow in pack_cell");
}


/******************************************************************************
*
*  unpack cell from MPI buffer to buffer cell (for force comp.)
*
******************************************************************************/

void unpack_cell( msgbuf *b, int j, int k )
{
  int i;
  int tmp_n;
  cell *to;

  to = PTR_2D_V(cell_array, j, k, cell_dim);

  tmp_n = (int) b->data[ b->n++ ];
  
  if (tmp_n > to->n_max) {
    to->n = 0;
    alloc_cell(to, tmp_n);
  }
  
  to->n = tmp_n;
  for (i=0; i<to->n; ++i) {
    to->ort X(i) = b->data[ b->n++ ];
    to->ort Y(i) = b->data[ b->n++ ];
    to->sorte[i] = (shortint) b->data[ b->n++ ];
  }
  if (b->n_max < b->n) error("Buffer overflow in unpack_cell");
}


/******************************************************************************
*
*  add forces of one cell to those of another cell
*
******************************************************************************/

void add_forces( int j, int k, int l, int m )
{
  int i;
  cell *from, *to;

  from = PTR_2D_V(cell_array, j, k, cell_dim);
  to   = PTR_2D_V(cell_array, l, m, cell_dim);

  for (i=0; i<to->n; ++i) {
    to->kraft X(i) += from->kraft X(i);
    to->kraft Y(i) += from->kraft Y(i);
    to->pot_eng[i] += from->pot_eng[i];
#ifdef NVX
    to->heatcond[i] += from->heatcond[i];
#endif
#ifdef STRESS_TENS
    to->presstens X(i)      += from->presstens X(i);
    to->presstens Y(i)      += from->presstens Y(i);
    to->presstens_offdia[i] += from->presstens_offdia[i];
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

void pack_forces( msgbuf *b, int j, int k )
 {
  int i;
  cell *from;
    
  from = PTR_2D_V(cell_array, j, k, cell_dim);
  for (i=0; i<from->n; ++i) {
    b->data[ b->n++ ] = from->kraft X(i);
    b->data[ b->n++ ] = from->kraft Y(i);
    b->data[ b->n++ ] = from->pot_eng[i];
#ifdef NVX
    b->data[ b->n++ ] = from->heatcond[i];
#endif
#ifdef STRESS_TENS
    b->data[ b->n++ ] = from->presstens X(i);
    b->data[ b->n++ ] = from->presstens Y(i);
    b->data[ b->n++ ] = from->presstens_offdia[i];
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

void unpack_forces( msgbuf *b, int j, int k )
 {
  int i;
  cell *to;

  to = PTR_2D_V(cell_array, j, k, cell_dim);
  for (i=0; i<to->n; ++i) {
    to->kraft X(i) += b->data[ b->n++ ];
    to->kraft Y(i) += b->data[ b->n++ ];
    to->pot_eng[i] += b->data[ b->n++ ];
#ifdef NVX
    to->heatcond[i] += b->data[ b->n++ ];
#endif
#ifdef STRESS_TENS
    to->presstens X(i)      += b->data[ b->n++ ];
    to->presstens Y(i)      += b->data[ b->n++ ];
    to->presstens_offdia[i] += b->data[ b->n++ ];
#endif
#ifdef ORDPAR
    to->nbanz[i] += (shortint) b->data[ b->n++ ];
#endif
  }
  if (b->n_max < b->n) error("Buffer overflow in unpack_forces.");
}
