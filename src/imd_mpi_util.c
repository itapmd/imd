
/******************************************************************************
*
* imd_mpi_util.c -- MPI utility routines
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/******************************************************************************
*
* isend_buf lean version of send_cell for force_loop
*
******************************************************************************/

int isend_buf(msgbuf *b, int to_cpu, MPI_Request *req)
{
  return MPI_Isend(b->data, b->n, REAL, to_cpu, BUFFER_TAG, cpugrid, req);
}

/******************************************************************************
*
* irecv_buf lean version of recv_cell for force_loop
*
******************************************************************************/

int irecv_buf(msgbuf *b, int from, MPI_Request *req)
{
  return MPI_Irecv(b->data, b->n_max, REAL, from, BUFFER_TAG, cpugrid, req);
}

/******************************************************************************
*
* copy an atom from a cell into a buffer, and possibly delete it in the cell
*
******************************************************************************/

void copy_one_atom(msgbuf *to, cell *from, int index, int delete)
{
  /* Check the parameters */
  if ((0 > index) || (index >= from->n)) {
    printf("%d: i %d n %d\n",myid,index,from->n);
    error("copy_one_atom: index argument out of range.");
  }

  /* See if we need some space */
  if (to->n >= to->n_max) error("copy_one_atom: buffer overflow.");

  /* copy atom */
  /* data is packed in the same order as in the cell data structure */
  to->data[ to->n++ ] = from->ort X(index); 
  to->data[ to->n++ ] = from->ort Y(index); 
#ifndef TWOD
  to->data[ to->n++ ] = from->ort Z(index); 
#endif
#ifndef MONOLJ
  to->data[ to->n++ ] = from->nummer[index];
  to->data[ to->n++ ] = from->sorte[index];
  to->data[ to->n++ ] = from->masse[index];
  to->data[ to->n++ ] = from->pot_eng[index];
#endif
#ifdef EAM2
  /* eam2_rho_h is not sent */
#endif
#ifdef DISLOC
  to->data[ to->n++ ] = from->Epot_ref[index];
  to->data[ to->n++ ] = from->ort_ref X(index); 
  to->data[ to->n++ ] = from->ort_ref Y(index); 
#ifndef TWOD
  to->data[ to->n++ ] = from->ort_ref Z(index); 
#endif
#endif 
#if defined(ORDPAR) && !defined(TWOD)
  to->data[ to->n++ ] = from->nbanz[index]; 
#endif
#ifdef REFPOS
  to->data[ to->n++ ] = from->refpos X(index);
  to->data[ to->n++ ] = from->refpos Y(index);
#ifndef TWOD
  to->data[ to->n++ ] = from->refpos Z(index);
#endif
#endif
#ifdef NVX
  to->data[ to->n++ ] = from->heat_cond[index];
#endif
#ifdef STRESS_TENS
  /* stress tensor is not sent */
#endif
  to->data[ to->n++ ] = from->impuls X(index); 
  to->data[ to->n++ ] = from->impuls Y(index); 
#ifndef TWOD
  to->data[ to->n++ ] = from->impuls Z(index); 
#endif
  /* force is not sent */
#ifdef COVALENT
  /* neighbor table is not sent */
#endif
#ifdef UNIAX
  to->data[ to->n++ ] = from->traeg_moment[index];
  to->data[ to->n++ ] = from->achse X(index); 
  to->data[ to->n++ ] = from->achse Y(index); 
  to->data[ to->n++ ] = from->achse Z(index); 
  to->data[ to->n++ ] = from->shape X(index); 
  to->data[ to->n++ ] = from->shape Y(index); 
  to->data[ to->n++ ] = from->shape Z(index); 
  to->data[ to->n++ ] = from->pot_well X(index); 
  to->data[ to->n++ ] = from->pot_well Y(index); 
  to->data[ to->n++ ] = from->pot_well Z(index); 
  to->data[ to->n++ ] = from->dreh_impuls X(index); 
  to->data[ to->n++ ] = from->dreh_impuls Y(index); 
  to->data[ to->n++ ] = from->dreh_impuls Z(index); 
  /* dreh_moment is not sent */
#endif

  /* Delete atom in original cell? */
  if (delete==1) {

    --from->n;

    /* we treat the data in the same order as in the cell data structure */
    if (0 < from->n) {

      from->ort X(index) = from->ort X(from->n); 
      from->ort Y(index) = from->ort Y(from->n); 
#ifndef TWOD
      from->ort Z(index) = from->ort Z(from->n); 
#endif
#ifndef MONOLJ
      from->nummer[index]  = from->nummer[from->n]; 
      from->sorte[index]   = from->sorte[from->n];
      from->masse[index]   = from->masse[from->n]; 
      from->pot_eng[index] = from->pot_eng[from->n]; 
#endif
#ifdef EAM2
      /* eam2_rho_h need not be copied */
#endif
#ifdef DISLOC
      from->Epot_ref[index]    = from->Epot_ref[from->n]; 
      from->ort_ref X(index)   = from->ort_ref X(from->n); 
      from->ort_ref Y(index)   = from->ort_ref Y(from->n); 
#ifndef TWOD
      from->ort_ref Z(index)   = from->ort_ref Z(from->n); 
#endif
#endif
#if defined(ORDPAR) && !defined(TWOD)
      from->nbanz[index]    = from->nbanz[from->n]; 
#endif
#ifdef REFPOS
      from->refpos X(index) = from->refpos X(from->n);
      from->refpos Y(index) = from->refpos Y(from->n);
#ifndef TWOD
      from->refpos Z(index) = from->refpos Z(from->n);
#endif
#endif
#ifdef NVX
      from->heatcond[index] = from->heatcond[from->n]; 
#endif
#ifdef STRESS_TENS
      /* stress tensor need not be copied */
#endif
      from->impuls X(index) = from->impuls X(from->n); 
      from->impuls Y(index) = from->impuls Y(from->n); 
#ifndef TWOD
      from->impuls Z(index) = from->impuls Z(from->n); 
#endif
      /* force need not be copied */
#ifdef COVALENT
      /* neighbor table need not be copied */
#endif
#ifdef UNIAX
      from->traeg_moment[index] = from->traeg_moment[from->n];
      from->achse X(index) = from->achse X(from->n); 
      from->achse Y(index) = from->achse Y(from->n); 
      from->achse Z(index) = from->achse Z(from->n); 
      from->shape X(index) = from->shape X(from->n); 
      from->shape Y(index) = from->shape Y(from->n); 
      from->shape Z(index) = from->shape Z(from->n); 
      from->pot_well X(index) = from->pot_well X(from->n); 
      from->pot_well Y(index) = from->pot_well Y(from->n); 
      from->pot_well Z(index) = from->pot_well Z(from->n); 
      from->dreh_impuls X(index) = from->dreh_impuls X(from->n); 
      from->dreh_impuls Y(index) = from->dreh_impuls Y(from->n); 
      from->dreh_impuls Z(index) = from->dreh_impuls Z(from->n); 
      /* dreh_moment need not be copied */
#endif
    }
  } /* delete or not delete */
}


/******************************************************************************
*
*  process buffer unpacks particles from mpi buffer.
*  if p == NULL, they are sorted into the cell array,
*  otherwise they are put into the cell pointed to by p.
*  mpi buffers have been packed with copy_one_atom.
*
******************************************************************************/

void process_buffer(msgbuf *b, cell *p)
{
  ivektor coord, coord2;
  int to_cpu, j;
  static cell *input = NULL;
  cell *to;

  if (NULL == input) {
    input = (cell *) malloc(sizeof(cell));
    if (0==input) error("Cannot allocate buffer cell.");
    input->n_max=0;
    alloc_cell(input, 1);
  }

  j=0;
  /* we treat the data in the same order as in the cell data structure */
  if (b->n > 0) do {     
    input->n        = 1;
    input->ort X(0) = b->data[j++];
    input->ort Y(0) = b->data[j++];
#ifndef TWOD
    input->ort Z(0) = b->data[j++];
#endif
#ifndef MONOLJ
    input->nummer[0]  = b->data[j++];
    input->sorte[0]   = b->data[j++];
    input->masse[0]   = b->data[j++];
    input->pot_eng[0] = b->data[j++];
#endif
#ifdef EAM2
    /* don't send eam2_rho_h */
#endif
#ifdef DISLOC
    input->Epot_ref[0]     = b->data[j++];
    input->ort_ref X(0)    = b->data[j++];
    input->ort_ref Y(0)    = b->data[j++];
#ifndef TWOD
    input->ort_ref Z(0)    = b->data[j++];
#endif
#endif
#if defined(ORDPAR) && !defined(TWOD)
    input->nbanz[0]        = b->data[j++];
#endif
#ifdef REFPOS
    input->refpos X(0)     = b->data[j++];
    input->refpos Y(0)     = b->data[j++];
#ifndef TWOD
    input->refpos Z(0)     = b->data[j++];
#endif
#endif
#ifdef NVX
    input->heatcond[0]     = b->data[j++];
#endif
#ifdef STRESS_TENS
    /* don't send stress tensor */
#endif
    input->impuls X(0)     = b->data[j++];
    input->impuls Y(0)     = b->data[j++];
#ifndef TWOD
    input->impuls Z(0)     = b->data[j++];
#endif
    /* don't send force */
#ifdef COVALENT
    /* don't send neighbor table */
#endif
#ifdef UNIAX
    input->traeg_moment[0] = b->data[j++];
    input->achse X(0) = b->data[j++];
    input->achse Y(0) = b->data[j++];
    input->achse Z(0) = b->data[j++];
    input->shape X(0) = b->data[j++];
    input->shape Y(0) = b->data[j++];
    input->shape Z(0) = b->data[j++];
    input->pot_well X(0) = b->data[j++];
    input->pot_well Y(0) = b->data[j++];
    input->pot_well Z(0) = b->data[j++];
    input->dreh_impuls X(0) = b->data[j++];
    input->dreh_impuls Y(0) = b->data[j++];
    input->dreh_impuls Z(0) = b->data[j++];
    /* don't send dreh_moment */
#endif

    if (p==NULL) {  /* distribute atom into cell array */
#ifdef TWOD
      coord =local_cell_coord( input->ort X(0), input->ort Y(0) );
      coord2=cell_coord(       input->ort X(0), input->ort Y(0) );
#else
      coord =local_cell_coord(input->ort X(0),input->ort Y(0),input->ort Z(0));
      coord2=cell_coord(      input->ort X(0),input->ort Y(0),input->ort Z(0));
#endif
      to_cpu = cpu_coord( coord2 );
      to = PTR_VV(cell_array,coord,cell_dim);
      if (to_cpu == myid) move_atom(to, input, 0);
    } else {  /* put atom into cell pointed by p */
      move_atom(p, input, 0);
    }
  } while (j<b->n);

}


/******************************************************************************
*
* copy_atoms_buf copies atoms from one message buffer to another
*
******************************************************************************/

void copy_atoms_buf(msgbuf *to, msgbuf *from)
{
  int i;

  if ((to->n_max - binc * CSTEP) < (to->n + from->n)) 
      error("Buffer overflow in copy_buf.");

  for (i=0; i<from->n; ++i) 
    to->data[ to->n++ ] = from->data[i];
}


/******************************************************************************
*
* setup_buffers sets up the send/receive buffers
*
* This is called periodically to check the buffers are large enough
*
* The buffers never shrink. 
* 
******************************************************************************/

void setup_buffers(void)
{
  int largest_cell, largest_local_cell=0;
  int k;
  int size_east;
  int size_north;
  int size_up;
  int binc1, binc2;

  /* determine buffer size per atom */
  if (binc==0) {

    /* for communication to buffer cells */
#ifdef TWOD
    binc1=2;     /* position */
#else
    binc1=3;     /* position */
#endif
#ifndef MONOLJ
    binc1++;     /* sorte */
#endif
#ifdef UNIAX
    binc1+=9;    /* direction, etc. */
#endif
#ifdef EAM
    binc1++;     /* nummer */
#endif

    /* for communication from buffer cells */
#ifdef TWOD
    binc2=2;     /* force */
#else
    binc2=3;     /* force */
#endif
#ifndef MONOLJ
    binc2++;     /* pot_eng */
#endif
#ifdef NVX
    binc2++;     /* heatcond */
#endif
#ifdef STRESS_TENS
#ifdef TWOD
    binc2+=3;    /* presstens */
#else
    binc2+=6;    /* presstens */
#endif
#endif
#ifdef ORDPAR
#ifndef TWOD
    binc2++;     /* nbanz */
#endif
#endif

    /* one way or two ways */
#ifdef AR
    binc=MAX(binc1,binc2);
#else
    binc=binc1;
#endif
  }

  /* Find largest cell */
  for (k=0; k<ncells; ++k) {
    int n;
    n = (cell_array + CELLS(k))->n;
    if (largest_local_cell < n) largest_local_cell = n;
  }

  MPI_Allreduce( &largest_local_cell, &largest_cell, 1, 
                 MPI_INT, MPI_MAX, cpugrid);

  /* Add security */
  largest_cell += CSTEP;

#ifndef TWOD
  size_east  = largest_cell * cell_dim.y * cell_dim.z * binc;
  size_north = largest_cell * cell_dim.x * cell_dim.z * binc;
  size_up    = largest_cell * cell_dim.x * cell_dim.y * binc;
#ifdef DEBUG
  if (0==myid) 
     printf("Max. cell is %d size east %d size north %d size up %d.\n", 
		      largest_cell,size_east,size_north,size_up);
#endif
#else
  size_east  = largest_cell * cell_dim.y * binc;
  size_north = largest_cell * cell_dim.x * binc;
#ifdef DEBUG
  if (0==myid) printf("Max. cell is %d size east %d size north %d.\n", 
		      largest_cell,size_east,size_north);
#endif
#endif

  /* Allocate east/west buffers */
  if (size_east > send_buf_east.n_max) {

    /* Buffer size is zero only at startup */
    if (0 != send_buf_east.n_max) {
      free(send_buf_east.data);
      free(send_buf_west.data);
      free(recv_buf_east.data);
      free(recv_buf_west.data);
    }

    /* Allocate buffers */
    send_buf_east.data = (real *) malloc( size_east * sizeof(real) );
    send_buf_west.data = (real *) malloc( size_east * sizeof(real) );
    recv_buf_east.data = (real *) malloc( size_east * sizeof(real) );
    recv_buf_west.data = (real *) malloc( size_east * sizeof(real) );

    /* Always check malloc results */
    if ((NULL == send_buf_east.data ) ||
	(NULL == send_buf_west.data ) ||
	(NULL == recv_buf_east.data ) ||
	(NULL == recv_buf_west.data )) 
        error("Cannot allocate east/west buffer.");

    send_buf_east.n_max = size_east;
    send_buf_west.n_max = size_east;
    recv_buf_east.n_max = size_east;
    recv_buf_west.n_max = size_east;
  }

  /* Allocate north/south buffers */
  if (size_north > send_buf_north.n_max) {

    /* Buffer size is zero only at startup */
    if (0 != send_buf_north.n_max) {
      free(send_buf_north.data);
      free(send_buf_south.data);
      free(recv_buf_north.data);
      free(recv_buf_south.data);
    }

    /* Allocate buffers */
    send_buf_north.data = (real *) malloc( size_north * sizeof(real) );
    send_buf_south.data = (real *) malloc( size_north * sizeof(real) );
    recv_buf_north.data = (real *) malloc( size_north * sizeof(real) );
    recv_buf_south.data = (real *) malloc( size_north * sizeof(real) );

    /* Always check malloc results */
    if ((NULL == send_buf_north.data ) ||
	(NULL == send_buf_south.data ) ||
	(NULL == recv_buf_north.data ) ||
	(NULL == recv_buf_south.data )) 
        error("Cannot allocate north/south buffer.");

    send_buf_north.n_max = size_north;
    send_buf_south.n_max = size_north;
    recv_buf_north.n_max = size_north;
    recv_buf_south.n_max = size_north;

  }

#ifndef TWOD
  /* Allocate up/down buffers */
  if (size_up > send_buf_up.n_max) {

    /* Buffer size is zero only at startup */
    if (0 != send_buf_up.n_max) {
      free(send_buf_up.data);
      free(send_buf_down.data);
      free(recv_buf_up.data);
      free(recv_buf_down.data);
    }

    /* Allocate buffers */
    send_buf_up.data   = (real *) malloc( size_up * sizeof(real) );
    send_buf_down.data = (real *) malloc( size_up * sizeof(real) );
    recv_buf_up.data   = (real *) malloc( size_up * sizeof(real) );
    recv_buf_down.data = (real *) malloc( size_up * sizeof(real) );

    /* Always check malloc results */
    if ((NULL == send_buf_up.data ) ||
	(NULL == send_buf_down.data ) ||
	(NULL == recv_buf_up.data ) ||
	(NULL == recv_buf_down.data )) 
        error("Cannot allocate up/down buffer.");

    send_buf_up.n_max   = size_up;
    send_buf_down.n_max = size_up;
    recv_buf_up.n_max   = size_up;
    recv_buf_down.n_max = size_up;
  }

#endif

}


/******************************************************************************
*
* empty mpi buffers -- set number of atoms in mpi buffers to zero
*
******************************************************************************/

void empty_mpi_buffers(void)
{
 /* Empty MPI buffers */
  send_buf_north.n = 0;
  send_buf_south.n = 0;
  send_buf_east.n  = 0;
  send_buf_west.n  = 0;

  recv_buf_north.n = 0;
  recv_buf_south.n = 0;
  recv_buf_east.n  = 0;
  recv_buf_west.n  = 0;

#ifndef TWOD
  recv_buf_down.n  = 0;
  recv_buf_up.n    = 0;
  send_buf_down.n  = 0;
  send_buf_up.n    = 0;
#endif

}


/******************************************************************************
*
*  send_cell packs the data of all particles in the cell
*  into a buffer, and sends the buffer to the destination cpu.
*
******************************************************************************/

void send_cell(cell *p, int to_cpu, int tag)
{
  int i;
  msgbuf *b;

  b = &send_buf_east;
  b->n = 0;

  for (i=0; i<p->n; i++) copy_one_atom(b, p, i, 0);
  if  (b->n_max < b->n)  error("Buffer overflow in send_cell");

  MPI_Send(b->data, b->n, REAL, to_cpu, tag, cpugrid);
}


/******************************************************************************
*
*  recv_cell receives in a buffer the data of the particles in one cell,
*  and upacks the buffer into the cell pointed to by p.
*
******************************************************************************/

void recv_cell(cell *p, int from_cpu, int tag)
{
  MPI_Status status;
  msgbuf *b;

  b = &recv_buf_east;

  /* check message size */
  MPI_Probe( from_cpu, tag, cpugrid, &status );
  MPI_Get_count( &status, REAL, &(b->n) );
  if (b->n_max < b->n) error("Buffer overflow in recv_cell");

  /* upack data */
  p->n = 0;
  MPI_Recv(b->data, b->n, REAL, from_cpu, tag, cpugrid, &status);
  process_buffer(b,p);
}


/******************************************************************************
*
*  send_cell_old send cell data to another process  (old version)
*
******************************************************************************/

void send_cell_old(cell *p, int to_cpu, int tag)
{
  MPI_Ssend( p->ort, DIM*p->n, REAL,    to_cpu, tag + ORT_TAG,    cpugrid);
#ifndef MONOLJ
  MPI_Ssend( p->sorte,  p->n,  SHORT,   to_cpu, tag + SORTE_TAG,  cpugrid);
  MPI_Ssend( p->masse,  p->n,  REAL,    to_cpu, tag + MASSE_TAG,  cpugrid);
  MPI_Ssend( p->nummer, p->n,  INTEGER, to_cpu, tag + NUMMER_TAG, cpugrid);
#ifdef UNIAX
  MPI_Ssend( p->traeg_moment, p->n, REAL,to_cpu,tag+TRAEG_MOMENT_TAG,cpugrid);
  MPI_Ssend( p->achse, DIM*p->n,    REAL, to_cpu, tag + ACHSE_TAG,    cpugrid);
  MPI_Ssend( p->shape, DIM*p->n,    REAL,to_cpu, tag + SHAPE_TAG,  cpugrid);
  MPI_Ssend( p->pot_well, DIM*p->n, REAL,to_cpu, tag + POT_WELL_TAG,  cpugrid);
  MPI_Ssend( p->dreh_impuls, DIM*p->n, REAL, to_cpu, tag + DREH_IMPULS_TAG, 
             cpugrid);
  MPI_Ssend( p->dreh_moment, DIM*p->n, REAL, to_cpu, tag + DREH_MOMENT_TAG,  
             cpugrid);
#endif
  MPI_Ssend( p->pot_eng,p->n,     REAL, to_cpu, tag + POT_TAG,    cpugrid);
#ifdef ORDPAR
#ifndef TWOD
  MPI_Ssend( p->nbanz,p->n,      SHORT, to_cpu, tag + NBA_TAG,    cpugrid);
#endif
#endif
#ifdef REFPOS
  MPI_Ssend( p->refpos, DIM*p->n, REAL, to_cpu, tag + REFPOS_TAG, cpugrid);
#endif
#ifdef DISLOC
  MPI_Ssend( p->Epot_ref,p->n,    REAL, to_cpu, tag + POT_REF_TAG,cpugrid);
  MPI_Ssend( p->ort_ref,DIM*p->n, REAL, to_cpu, tag + ORT_REF_TAG,cpugrid);
#endif
#endif /* not MONOLJ */
  MPI_Ssend( p->impuls, DIM*p->n, REAL, to_cpu, tag + IMPULS_TAG, cpugrid);
  MPI_Ssend( p->kraft,  DIM*p->n, REAL, to_cpu, tag + KRAFT_TAG,  cpugrid);
}


/******************************************************************************
*
* recv_cell_old receive cell data for cell  (old version)
*
******************************************************************************/

void recv_cell_old(cell *p, int from_cpu, int tag)
{
  int size;
  int newsize;
  MPI_Status status;

  MPI_Probe( from_cpu, tag + ORT_TAG, cpugrid, &status );
  MPI_Get_count( &status, REAL, &size );
  size /= DIM;

  /* realloc cell if necessary */
  newsize = p->n_max; 
  while ( newsize < size ) newsize += incrsz;
  if (newsize > p->n_max) {
    /* Inihibit superfluous copy operation */
    p->n = 0;
    alloc_cell(p,newsize);
  }
  p->n = size;

  MPI_Recv(p->ort,     DIM * size, REAL,     from_cpu, tag + ORT_TAG,
                                             cpugrid, &status );
#ifndef MONOLJ
  MPI_Recv(p->sorte,         size, SHORT,    from_cpu, tag + SORTE_TAG , 
                                             cpugrid, &status);
  MPI_Recv(p->masse,         size, REAL,     from_cpu, tag + MASSE_TAG , 
                                             cpugrid, &status);
  MPI_Recv(p->nummer,        size, INTEGER,  from_cpu, tag + NUMMER_TAG, 
                                             cpugrid, &status);
#ifdef UNIAX
  MPI_Recv(p->traeg_moment,      size, REAL, from_cpu, tag + TRAEG_MOMENT_TAG, 
                                             cpugrid, &status);
  MPI_Recv(p->achse,       DIM * size, REAL, from_cpu, tag + ACHSE_TAG,
                                             cpugrid, &status );
  MPI_Recv(p->shape,       DIM * size, REAL, from_cpu, tag + SHAPE_TAG,
                                             cpugrid, &status );
  MPI_Recv(p->pot_well,    DIM * size, REAL, from_cpu, tag + POT_WELL_TAG,
                                             cpugrid, &status );
  MPI_Recv(p->dreh_impuls, DIM * size, REAL, from_cpu, tag + DREH_IMPULS_TAG, 
                                             cpugrid, &status);
  MPI_Recv(p->dreh_moment, DIM * size, REAL, from_cpu, tag + DREH_MOMENT_TAG,  
                                             cpugrid, &status);
#endif
  MPI_Recv(p->pot_eng, size, REAL,  from_cpu, tag + POT_TAG, cpugrid, &status);
#ifdef ORDPAR
#ifndef TWOD
  MPI_Recv(p->nbanz,   size, SHORT, from_cpu, tag + NBA_TAG, cpugrid, &status);
#endif
#endif
                                            
#ifdef REFPOS
  MPI_Recv(p->refpos,  DIM * size, REAL, from_cpu, tag + REFPOS_TAG, 
                                         cpugrid, &status);
#endif
#ifdef DISLOC
  MPI_Recv(p->Epot_ref,      size, REAL, from_cpu, tag + POT_REF_TAG , 
                                         cpugrid, &status);
  MPI_Recv(p->ort_ref, DIM * size, REAL, from_cpu, tag + ORT_REF_TAG, 
                                         cpugrid, &status);
#endif
#endif /* not MONOLJ */
  MPI_Recv(p->impuls,  DIM * size, REAL, from_cpu, tag + IMPULS_TAG, 
                                         cpugrid, &status);
  MPI_Recv(p->kraft,   DIM * size, REAL, from_cpu, tag + KRAFT_TAG,  
                                         cpugrid, &status);
}
