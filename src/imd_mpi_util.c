
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
* set up mpi
*
******************************************************************************/

void init_mpi(int *argc_pointer, char **argv)
{
  /* Initialize MPI */
  MPI_Init(argc_pointer, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&num_cpus);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  if (0 == myid) { 
    printf("%s\n", argv[0]);
    printf("Starting up MPI on %d nodes.\n", num_cpus);
  }
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

#ifdef SR
/******************************************************************************
*
* sendrecv_buf lean version of sendrecv for force_loop
*
******************************************************************************/

int sendrecv_buf(msgbuf *send, int to_cpu, 
                 msgbuf *recv, int from_cpu, MPI_Status *status)
{
  return MPI_Sendrecv(send->data, send->n, REAL, to_cpu, BUFFER_TAG,
		      recv->data, recv->n_max, REAL, from_cpu, BUFFER_TAG,
		      cpugrid, status );
}
#endif

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
  if (to->n + MAX_ATOM_SIZE > to->n_max) {
    to->n_max += BUFFER_SIZE_INC;
    to->data = (real *) realloc( to->data, to->n_max * sizeof(real) );
    if (NULL == to->data) error("Cannot allocate buffer in copy_one_atom");
  }

  /* copy atom */
  /* data is packed in the same order as in the cell data structure */
  to->data[ to->n++ ] = ORT(from,index,X); 
  to->data[ to->n++ ] = ORT(from,index,Y); 
#ifndef TWOD
  to->data[ to->n++ ] = ORT(from,index,Z); 
#endif

#ifndef MONOLJ
  to->data[ to->n++ ] = NUMMER(from,index);
  to->data[ to->n++ ] = VSORTE(from,index);
  to->data[ to->n++ ] = MASSE (from,index);
  to->data[ to->n++ ] = POTENG(from,index);
#endif
#ifdef EAM2
  /* eam2_rho_h is not sent */
#endif

#ifdef CG
  to->data[ to->n++ ] = CG_H(from,index,X); 
  to->data[ to->n++ ] = CG_H(from,index,Y); 
#ifndef TWOD
  to->data[ to->n++ ] = CG_H(from,index,Z); 
#endif
  to->data[ to->n++ ] = CG_G(from,index,X); 
  to->data[ to->n++ ] = CG_G(from,index,Y); 
#ifndef TWOD
  to->data[ to->n++ ] = CG_G(from,index,Z); 
#endif
  to->data[ to->n++ ] = OLD_ORT(from,index,X); 
  to->data[ to->n++ ] = OLD_ORT(from,index,Y); 
#ifndef TWOD
  to->data[ to->n++ ] = OLD_ORT(from,index,Z); 
#endif
#endif /* CG */

#ifdef DISLOC
  to->data[ to->n++ ] = EPOT_REF(from,index);
  to->data[ to->n++ ] = ORT_REF (from,index,X); 
  to->data[ to->n++ ] = ORT_REF (from,index,Y); 
#ifndef TWOD
  to->data[ to->n++ ] = ORT_REF (from,index,Z); 
#endif
#endif 
#ifdef AVPOS
  to->data[ to->n++ ] = AV_EPOT(from,index);
  to->data[ to->n++ ] = AV_POS (from,index,X); 
  to->data[ to->n++ ] = AV_POS (from,index,Y);
  to->data[ to->n++ ] = SHEET  (from,index,X);
  to->data[ to->n++ ] = SHEET  (from,index,Y);
#ifndef TWOD
  to->data[ to->n++ ] = AV_POS (from,index,Z); 
  to->data[ to->n++ ] = SHEET  (from,index,Z);
#endif
#endif
#ifdef ORDPAR  
  to->data[ to->n++ ] = NBANZ(from,index); 
#endif
#ifdef REFPOS
  to->data[ to->n++ ] = REF_POS(from,index,X);
  to->data[ to->n++ ] = REF_POS(from,index,Y);
#ifndef TWOD
  to->data[ to->n++ ] = REF_POS(from,index,Z);
#endif
#endif
#ifdef NVX
  to->data[ to->n++ ] = HEATCOND(from,index);
#endif
#ifdef STRESS_TENS
  to->data[ to->n++ ] = PRESSTENS(from,index,xx);   
  to->data[ to->n++ ] = PRESSTENS(from,index,yy);   
  to->data[ to->n++ ] = PRESSTENS(from,index,xy);   
#ifndef TWOD
  to->data[ to->n++ ] = PRESSTENS(from,index,zz);   
  to->data[ to->n++ ] = PRESSTENS(from,index,yz);   
  to->data[ to->n++ ] = PRESSTENS(from,index,zx);   
#endif
#endif /* STRESS_TENS */
  to->data[ to->n++ ] = IMPULS(from,index,X); 
  to->data[ to->n++ ] = IMPULS(from,index,Y); 
#ifndef TWOD
  to->data[ to->n++ ] = IMPULS(from,index,Z); 
#endif

  /* force is not sent */
#ifdef COVALENT
  /* neighbor table is not sent */
#endif
#ifdef UNIAX
  to->data[ to->n++ ] = TRAEG_MOMENT(from,index);
  to->data[ to->n++ ] = ACHSE(from,index,X); 
  to->data[ to->n++ ] = ACHSE(from,index,Y); 
  to->data[ to->n++ ] = ACHSE(from,index,Z); 
  to->data[ to->n++ ] = SHAPE(from,index,X); 
  to->data[ to->n++ ] = SHAPE(from,index,Y); 
  to->data[ to->n++ ] = SHAPE(from,index,Z); 
  to->data[ to->n++ ] = POT_WELL(from,index,X); 
  to->data[ to->n++ ] = POT_WELL(from,index,Y); 
  to->data[ to->n++ ] = POT_WELL(from,index,Z); 
  to->data[ to->n++ ] = DREH_IMPULS(from,index,X); 
  to->data[ to->n++ ] = DREH_IMPULS(from,index,Y); 
  to->data[ to->n++ ] = DREH_IMPULS(from,index,Z); 
  /* dreh_moment is not sent */
#endif

  /* Delete atom in original cell? */
  if (delete==1) {

    --from->n;

    /* we treat the data in the same order as in the cell data structure */
    if (0 < from->n) {

      ORT(from,index,X)  = ORT(from,from->n,X); 
      ORT(from,index,Y)  = ORT(from,from->n,Y); 
#ifndef TWOD
      ORT(from,index,Z)  = ORT(from,from->n,Z); 
#endif
#ifndef MONOLJ
      NUMMER(from,index) = NUMMER(from,from->n); 
      VSORTE(from,index) = VSORTE(from,from->n);
      MASSE (from,index) = MASSE (from,from->n); 
      POTENG(from,index) = POTENG(from,from->n); 
#endif
#ifdef EAM2
      /* eam2_rho_h need not be copied */
#endif

#ifdef CG
      CG_H(from,index,X) = CG_H(from,from->n,X); 
      CG_H(from,index,Y) = CG_H(from,from->n,Y); 
#ifndef TWOD
      CG_H(from,index,Z) = CG_H(from,from->n,Z); 
#endif
      CG_G(from,index,X) = CG_G(from,from->n,X); 
      CG_G(from,index,Y) = CG_G(from,from->n,Y); 
#ifndef TWOD
      CG_G(from,index,Z) = CG_G(from,from->n,Z); 
#endif
      OLD_ORT(from,index,X) = OLD_ORT(from,from->n,X); 
      OLD_ORT(from,index,Y) = OLD_ORT(from,from->n,Y); 
#ifndef TWOD
      OLD_ORT(from,index,Z) = OLD_ORT(from,from->n,Z); 
#endif
#endif /* CG */

#ifdef DISLOC
      EPOT_REF(from,index)   = EPOT_REF(from,from->n); 
      ORT_REF (from,index,X) = ORT_REF (from,from->n,X); 
      ORT_REF (from,index,Y) = ORT_REF (from,from->n,Y); 
#ifndef TWOD
      ORT_REF (from,index,Z) = ORT_REF (from,from->n,Z); 
#endif
#endif
#ifdef AVPOS
      AV_EPOT(from,index)   = AV_EPOT(from,from->n); 
      AV_POS (from,index,X) = AV_POS (from,from->n,X); 
      AV_POS (from,index,Y) = AV_POS (from,from->n,Y); 
      SHEET  (from,index,X) = SHEET  (from,from->n,X); 
      SHEET  (from,index,Y) = SHEET  (from,from->n,Y); 
#ifndef TWOD
      AV_POS (from,index,Z) = AV_POS (from,from->n,Z); 
      SHEET  (from,index,Z) = SHEET  (from,from->n,Z); 
#endif
#endif
#ifdef ORDPAR 
      NBANZ(from,index)    = NBANZ(from,from->n); 
#endif
#ifdef REFPOS
      REF_POS(from,index,X) = REF_POS(from,from->n,X);
      REF_POS(from,index,Y) = REF_POS(from,from->n,Y);
#ifndef TWOD
      REF_POS(from,index,Z) = REF_POS(from,from->n,Z);
#endif
#endif
#ifdef NVX
      HEATCOND(from,index)  = HEATCOND(from,from->n); 
#endif
#ifdef STRESS_TENS
      PRESSTENS(from,index,xx) = PRESSTENS(from,from->n,xx);   
      PRESSTENS(from,index,yy) = PRESSTENS(from,from->n,yy);   
      PRESSTENS(from,index,xy) = PRESSTENS(from,from->n,xy);   
#ifndef TWOD
      PRESSTENS(from,index,zz) = PRESSTENS(from,from->n,zz);   
      PRESSTENS(from,index,yz) = PRESSTENS(from,from->n,yz);   
      PRESSTENS(from,index,zx) = PRESSTENS(from,from->n,zx);   
#endif
#endif /* STRESS_TENS */
      IMPULS(from,index,X) = IMPULS(from,from->n,X); 
      IMPULS(from,index,Y) = IMPULS(from,from->n,Y); 
#ifndef TWOD
      IMPULS(from,index,Z) = IMPULS(from,from->n,Z); 
#endif
      /* force need not be copied */
#ifdef COVALENT
      /* neighbor table need not be copied */
#endif
#ifdef UNIAX
      TRAEG_MOMENT(from,index) = TRAEG_MOMENT(from,from->n);
      ACHSE(from,index,X) = ACHSE(from,from->n,X); 
      ACHSE(from,index,Y) = ACHSE(from,from->n,Y); 
      ACHSE(from,index,Z) = ACHSE(from,from->n,Z); 
      SHAPE(from,index,X) = SHAPE(from,from->n,X); 
      SHAPE(from,index,Y) = SHAPE(from,from->n,Y); 
      SHAPE(from,index,Z) = SHAPE(from,from->n,Z); 
      POT_WELL(from,index,X) = POT_WELL(from,from->n,X); 
      POT_WELL(from,index,Y) = POT_WELL(from,from->n,Y); 
      POT_WELL(from,index,Z) = POT_WELL(from,from->n,Z); 
      DREH_IMPULS(from,index,X) = DREH_IMPULS(from,from->n,X); 
      DREH_IMPULS(from,index,Y) = DREH_IMPULS(from,from->n,Y); 
      DREH_IMPULS(from,index,Z) = DREH_IMPULS(from,from->n,Z); 
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

  if (b->n > b->n_max) error("buffer overflow in process_buffer");

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
    ORT(input,0,X)  = b->data[j++];
    ORT(input,0,Y)  = b->data[j++];
#ifndef TWOD
    ORT(input,0,Z)  = b->data[j++];
#endif
#ifndef MONOLJ
    NUMMER(input,0) = b->data[j++];
    VSORTE(input,0) = b->data[j++];
    MASSE (input,0) = b->data[j++];
    POTENG(input,0) = b->data[j++];
#endif
#ifdef EAM2
    /* don't send eam2_rho_h */
#endif

#ifdef CG
    CG_H(input,0,X) = b->data[j++];
    CG_H(input,0,Y) = b->data[j++];
#ifndef TWOD
    CG_H(input,0,Z) = b->data[j++];
#endif
    CG_G(input,0,X) = b->data[j++];
    CG_G(input,0,Y) = b->data[j++];
#ifndef TWOD
    CG_G(input,0,Z) = b->data[j++];
#endif
    OLD_ORT(input,0,X) = b->data[j++];
    OLD_ORT(input,0,Y) = b->data[j++];
#ifndef TWOD
    OLD_ORT(input,0,Z) = b->data[j++];
#endif
#endif /* CG */

#ifdef DISLOC
    EPOT_REF(input,0)   = b->data[j++];
    ORT_REF (input,0,X) = b->data[j++];
    ORT_REF (input,0,Y) = b->data[j++];
#ifndef TWOD
    ORT_REF (input,0,Z) = b->data[j++];
#endif
#endif
#ifdef AVPOS
    AV_EPOT(input,0)    = b->data[j++];
    AV_POS (input,0,X)  = b->data[j++];
    AV_POS (input,0,Y)  = b->data[j++];
    SHEET  (input,0,X)  = b->data[j++];
    SHEET  (input,0,Y)  = b->data[j++];
#ifndef TWOD
    AV_POS (input,0,Z)  = b->data[j++];
    SHEET  (input,0,Z)  = b->data[j++];
#endif
#endif
#ifdef ORDPAR
    NBANZ(input,0)      = b->data[j++];
#endif
#ifdef REFPOS
    REF_POS(input,0,X)  = b->data[j++];
    REF_POS(input,0,Y)  = b->data[j++];
#ifndef TWOD
    REF_POS(input,0,Z)  = b->data[j++];
#endif
#endif
#ifdef NVX
    HEATCOND(input,0)   = b->data[j++];
#endif
#ifdef STRESS_TENS
    PRESSTENS(input,0,xx) = b->data[j++];   
    PRESSTENS(input,0,yy) = b->data[j++];   
    PRESSTENS(input,0,xy) = b->data[j++];   
#ifndef TWOD
    PRESSTENS(input,0,zz) = b->data[j++];   
    PRESSTENS(input,0,yz) = b->data[j++];   
    PRESSTENS(input,0,zx) = b->data[j++];   
#endif
#endif /* STRESS_TENS */
    IMPULS(input,0,X)     = b->data[j++];
    IMPULS(input,0,Y)     = b->data[j++];
#ifndef TWOD
    IMPULS(input,0,Z)     = b->data[j++];
#endif
    /* don't send force */
#ifdef COVALENT
    /* don't send neighbor table */
#endif
#ifdef UNIAX
    TRAEG_MOMENT(input,0) = b->data[j++];
    ACHSE(input,0,X) = b->data[j++];
    ACHSE(input,0,Y) = b->data[j++];
    ACHSE(input,0,Z) = b->data[j++];
    SHAPE(input,0,X) = b->data[j++];
    SHAPE(input,0,Y) = b->data[j++];
    SHAPE(input,0,Z) = b->data[j++];
    POT_WELL(input,0,X) = b->data[j++];
    POT_WELL(input,0,Y) = b->data[j++];
    POT_WELL(input,0,Z) = b->data[j++];
    DREH_IMPULS(input,0,X) = b->data[j++];
    DREH_IMPULS(input,0,Y) = b->data[j++];
    DREH_IMPULS(input,0,Z) = b->data[j++];
    /* don't send dreh_moment */
#endif

    if (p==NULL) {  /* distribute atom into cell array */
#ifdef TWOD
      coord =local_cell_coord( ORT(input,0,X), ORT(input,0,Y) );
      coord2=cell_coord(       ORT(input,0,X), ORT(input,0,Y) );
#else
      coord =local_cell_coord( ORT(input,0,X), ORT(input,0,Y), ORT(input,0,Z));
      coord2=cell_coord(       ORT(input,0,X), ORT(input,0,Y), ORT(input,0,Z));
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

  if (to->n_max < to->n + from->n) {
    while (to->n_max < to->n + from->n) to->n_max += BUFFER_SIZE_INC;
    to->data = (real *) realloc( to->data, to->n_max * sizeof(real) );
    if (NULL == to->data) error("Cannot allocate buffer in copy_atoms_buf");
  }

  for (i=0; i<from->n; ++i) 
    to->data[ to->n++ ] = from->data[i];
}


/******************************************************************************
*
* setup_buffers sets up the send/receive buffers
* This is called periodically to check the buffers are large enough
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
    binc1 = DIM;     /* position */
#ifndef MONOLJ
    binc1++;         /* sorte */
#endif
#ifdef UNIAX
    binc1 += 9;      /* direction, etc. */
#endif

    /* for communication from buffer cells */
    binc2 = DIM;     /* force */
#ifndef MONOLJ
    binc2++;         /* pot_eng */
#endif
#ifdef NVX
    binc2++;         /* heatcond */
#endif
#ifdef STRESS_TENS
#ifdef TWOD
    binc2 += 3;      /* presstens */
#else
    binc2 += 6;      /* presstens */
#endif
#endif
#ifdef ORDPAR
    binc2++;         /* nbanz */
#endif
#ifdef UNIAX
    binc2 += 3;      /* dreh_moment */
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
  if (b->n_max < b->n) {
    if (0 != b->n_max) free(b->data);
    while (b->n_max < b->n) b->n_max += BUFFER_SIZE_INC;
    b->data = (real *) malloc( b->n_max * sizeof(real) );
    if (NULL == b->data) error("Cannot allocate buffer in recv_cell");
  }

  /* upack data */
  p->n = 0;
  MPI_Recv(b->data, b->n, REAL, from_cpu, tag, cpugrid, &status);
  process_buffer(b,p);
}


















