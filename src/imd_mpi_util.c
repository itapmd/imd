
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
* imd_mpi_util.c -- MPI utility routines
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"
#ifdef BGL
#include <rts.h>
#include <bglpersonality.h>
#endif

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
    printf("Starting up MPI with %d processes.\n", num_cpus);
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

/******************************************************************************
*
* initialize I/O parameters
*
******************************************************************************/

void init_io(void)
{

#ifdef BGL
  BGLPersonality personality;
  int *tmp, i;

  tmp     = (int *) malloc( num_cpus * sizeof(int) );
  io_grps = (int *) malloc( num_cpus * sizeof(int) );
  if ((NULL==tmp) || (NULL==io_grps)) 
    error("cannot allocate io_grps array");

  rts_get_personality(&personality, sizeof(personality));

  /* I/O group as a function of the rank */
  for (i=0; i<num_cpus; i++) tmp[i] = 0;
  tmp[myid] = BGLPersonality_psetNum(&personality);
  MPI_Allreduce(tmp, io_grps, num_cpus, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  /* input parameters */
  if (parallel_input==1) {
    n_inp_grps   = BGLPersonality_numPsets(&personality);
    my_inp_grp   = BGLPersonality_psetNum (&personality);
    inp_grp_size = BGLPersonality_numNodesInPset(&personality);
    if (BGLPersonality_virtualNodeMode(&personality)) inp_grp_size *= 2;
    my_inp_id = 0;  while (my_inp_grp != io_grps[my_inp_id]) my_inp_id++;
  }
  else {
    n_inp_grps   = 1;
    my_inp_grp   = 0;
    my_inp_id    = 0;
    inp_grp_size = num_cpus;
  }

  /* output parameters */
  if (parallel_output==1) {
    n_out_grps   = BGLPersonality_numPsets(&personality);
    my_out_grp   = BGLPersonality_psetNum (&personality);
    out_grp_size = BGLPersonality_numNodesInPset(&personality);
    if (BGLPersonality_virtualNodeMode(&personality)) out_grp_size *= 2;
    my_out_id = 0;  while (my_out_grp != io_grps[my_out_id]) my_out_id++;
  }
  else {
    n_out_grps   = 1;
    my_out_grp   = 0;
    my_out_id    = 0;
    out_grp_size = num_cpus;
  }

#ifdef DEBUG
  /* some debug output */
  if (myid==my_inp_id)
    printf("myid=%d inpgrp=%d inpgrp_size=%d n_inpgrps=%d n_cpus=%d\n", 
      myid, my_inp_grp, inp_grp_size, n_inp_grps, num_cpus);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myid==my_out_id)
    printf("myid=%d outgrp=%d outgrp_size=%d n_outgrps=%d n_cpus=%d\n", 
      myid, my_out_grp, out_grp_size, n_out_grps, num_cpus);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myid==0) {
    for (i=0; i<num_cpus; i++) printf("%d ", io_grps[i]);
    printf("\n");
  }
#else
  if (0==myid) 
    printf("%d input group(s), %d output group(s)\n", n_inp_grps, n_out_grps);
#endif

  free(tmp);

#else /* not BGL */

  /* input parameters */
  if (parallel_input==1) {
    n_inp_grps   = num_cpus;
    my_inp_grp   = myid;
    my_inp_id    = myid;
    inp_grp_size = 1;
  }
  else {
    n_inp_grps   = 1;
    my_inp_grp   = 0;
    my_inp_id    = 0;
    inp_grp_size = num_cpus;
  }

  /* output parameters */
  if (parallel_output==1) {
    n_out_grps = num_cpus;
    my_out_grp = myid;
    my_out_id  = myid;
    out_grp_size = 1;
  }
  else {
    n_out_grps   = 1;
    my_out_grp   = 0;
    my_out_id    = 0;
    out_grp_size = num_cpus;
  }

#endif

}

/******************************************************************************
*
* allocate/deallocate a message buffer
*
******************************************************************************/

void alloc_msgbuf(msgbuf *b, int size)
{
#ifdef MPI2
  if (b->data) MPI_Free_mem(b->data);
  MPI_Alloc_mem( size * sizeof(real), MPI_INFO_NULL, &(b->data) );
#else
  free(b->data);
  b->data = (real *) malloc( size * sizeof(real) );
#endif
  if (NULL == b->data) error("cannot allocate message buffer");
  b->n = 0;
  b->n_max = size;
}

void realloc_msgbuf(msgbuf *b, int size)
{
  int  i;
  real *new;

#ifdef MPI2
  MPI_Alloc_mem( size * sizeof(real), MPI_INFO_NULL, &new );
#else
  new = (real *) malloc( size * sizeof(real) );
#endif
  if (NULL == new) error("cannot allocate message buffer");
  for (i=0; i<b->n; i++) new[i] = b->data[i];
#ifdef MPI2
  if (b->data) MPI_Free_mem(b->data);
#else
  free(b->data);
#endif
  b->data = new;
  b->n_max = size;
}

void free_msgbuf(msgbuf *b)
{
#ifdef MPI2
  if (b->data) MPI_Free_mem(b->data);
#else
  free(b->data);
#endif
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
* copy an atom from a cell into a buffer
*
******************************************************************************/

void copy_atom(msgbuf *to, int to_cpu, cell *p, int ind )
{
  /* Check the parameters */
  if ((0 > ind) || (ind >= p->n)) {
    printf("%d: i %d n %d\n",myid,ind,p->n);
    error("copy_atom: index argument out of range.");
  }

  /* See if we need some space */
  if (to->n + MAX_ATOM_SIZE > to->n_max) {
    realloc_msgbuf(to, to->n_max + BUFFER_SIZE_INC);
  }

  /* copy atom */
  /* data is packed in the same order as in the cell data structure */
  to->data[ to->n++ ] = to_cpu; 
  to->data[ to->n++ ] = ORT(p,ind,X); 
  to->data[ to->n++ ] = ORT(p,ind,Y); 
#ifndef TWOD
  to->data[ to->n++ ] = ORT(p,ind,Z); 
#endif

#ifndef MONOLJ
  to->data[ to->n++ ] = NUMMER(p,ind);
#ifndef MONO
  to->data[ to->n++ ] = SORTE (p,ind);
#endif
  to->data[ to->n++ ] = VSORTE(p,ind);
  to->data[ to->n++ ] = MASSE (p,ind);
  to->data[ to->n++ ] = POTENG(p,ind);
#endif
#ifdef EAM2
  to->data[ to->n++ ] = EAM_RHO(p,ind);
  /* eam_dF  is not sent */
#ifdef EEAM
  to->data[ to->n++ ] = EAM_P  (p,ind);
  /* eam_dM  is not sent */
#endif
#endif
#ifdef ADP
  /* adp_mu and adp_lambda are not sent */
#endif
#ifdef CG
  to->data[ to->n++ ] = CG_H(p,ind,X); 
  to->data[ to->n++ ] = CG_H(p,ind,Y); 
#ifndef TWOD
  to->data[ to->n++ ] = CG_H(p,ind,Z); 
#endif
  to->data[ to->n++ ] = CG_G(p,ind,X); 
  to->data[ to->n++ ] = CG_G(p,ind,Y); 
#ifndef TWOD
  to->data[ to->n++ ] = CG_G(p,ind,Z); 
#endif
  to->data[ to->n++ ] = OLD_ORT(p,ind,X); 
  to->data[ to->n++ ] = OLD_ORT(p,ind,Y); 
#ifndef TWOD
  to->data[ to->n++ ] = OLD_ORT(p,ind,Z); 
#endif
#endif /* CG */
#ifdef DAMP
  to->data[ to->n++ ] = DAMPF(p,ind);
#endif
#ifdef DISLOC
  to->data[ to->n++ ] = EPOT_REF(p,ind);
  to->data[ to->n++ ] = ORT_REF (p,ind,X); 
  to->data[ to->n++ ] = ORT_REF (p,ind,Y); 
#ifndef TWOD
  to->data[ to->n++ ] = ORT_REF (p,ind,Z); 
#endif
#endif 
#ifdef AVPOS
  to->data[ to->n++ ] = AV_EPOT(p,ind);
  to->data[ to->n++ ] = AV_POS (p,ind,X); 
  to->data[ to->n++ ] = AV_POS (p,ind,Y);
  to->data[ to->n++ ] = SHEET  (p,ind,X);
  to->data[ to->n++ ] = SHEET  (p,ind,Y);
#ifndef TWOD
  to->data[ to->n++ ] = AV_POS (p,ind,Z); 
  to->data[ to->n++ ] = SHEET  (p,ind,Z);
#endif
#endif
#ifdef NNBR  
  to->data[ to->n++ ] = NBANZ(p,ind); 
#endif
#ifdef REFPOS
  to->data[ to->n++ ] = REF_POS(p,ind,X);
  to->data[ to->n++ ] = REF_POS(p,ind,Y);
#ifndef TWOD
  to->data[ to->n++ ] = REF_POS(p,ind,Z);
#endif
#endif
#ifdef NVX
  to->data[ to->n++ ] = HEATCOND(p,ind);
#endif
#ifdef STRESS_TENS
  to->data[ to->n++ ] = PRESSTENS(p,ind,xx);   
  to->data[ to->n++ ] = PRESSTENS(p,ind,yy);   
  to->data[ to->n++ ] = PRESSTENS(p,ind,xy);   
#ifndef TWOD
  to->data[ to->n++ ] = PRESSTENS(p,ind,zz);   
  to->data[ to->n++ ] = PRESSTENS(p,ind,yz);   
  to->data[ to->n++ ] = PRESSTENS(p,ind,zx);   
#endif
#endif /* STRESS_TENS */
#ifdef SHOCK
  to->data[ to->n++ ] = PXAVG(p,ind);
#endif
  to->data[ to->n++ ] = IMPULS(p,ind,X); 
  to->data[ to->n++ ] = IMPULS(p,ind,Y); 
#ifndef TWOD
  to->data[ to->n++ ] = IMPULS(p,ind,Z); 
#endif

  /* force is not sent */
#ifdef COVALENT
  /* neighbor table is not sent */
#endif
#ifdef NBLIST
  /* neighbor list reference positions are not sent */
#endif
#ifdef UNIAX
  to->data[ to->n++ ] = ACHSE(p,ind,X); 
  to->data[ to->n++ ] = ACHSE(p,ind,Y); 
  to->data[ to->n++ ] = ACHSE(p,ind,Z); 
  to->data[ to->n++ ] = DREH_IMPULS(p,ind,X); 
  to->data[ to->n++ ] = DREH_IMPULS(p,ind,Y); 
  to->data[ to->n++ ] = DREH_IMPULS(p,ind,Z); 
  /* dreh_moment is not sent */
#endif
#ifdef VEC
  /* ind is not sent */
#endif

}

/******************************************************************************
*
* copy an atom from a minicell into a buffer, and possibly delete it
*
******************************************************************************/

void copy_one_atom(msgbuf *to, int to_cpu, minicell *from, int index, int del)
{
  cell *p;
  int  ind;

#ifdef VEC
  p   = &atoms;
  ind = from->ind[index];
#else
  p   = from;
  ind = index;
#endif

  /* copy the atom to the message buffer */
  copy_atom(to, to_cpu, p, ind);

  /* Delete atom in original cell? */
  if (del==1) {

    p->n--;

#ifdef VEC
    /* we move the last atom to slot ind, so we have to correct */
    /* the minicell entry pointing to the atom to be moved      */
    if (0 < p->n) {
      ivektor  coord;
      minicell *last;
      coord = cell_coord( ORT(p,p->n,X), ORT(p,p->n,Y), ORT(p,p->n,Z) );
      coord = local_cell_coord( coord );
      last  = PTR_VV(cell_array,coord,cell_dim);
      last->ind[ p->ind[p->n] ] = ind;
      p->ind[ind] = p->ind[p->n];
    }
    from->n--;
    if (0 < from->n) from->ind[index] = from->ind[from->n];
#endif

    /* we treat the data in the same order as in the cell data structure */
    if (0 < p->n) {

      ORT(p,ind,X)  = ORT(p,p->n,X); 
      ORT(p,ind,Y)  = ORT(p,p->n,Y); 
#ifndef TWOD
      ORT(p,ind,Z)  = ORT(p,p->n,Z); 
#endif
#ifndef MONOLJ
      NUMMER(p,ind) = NUMMER(p,p->n); 
#ifndef MONO
      SORTE (p,ind) = SORTE (p,p->n);
#endif
      VSORTE(p,ind) = VSORTE(p,p->n);
      MASSE (p,ind) = MASSE (p,p->n); 
      POTENG(p,ind) = POTENG(p,p->n); 
#endif
#ifdef EAM2
      EAM_RHO(p,ind) = EAM_RHO(p,p->n); 
      /* eam_dF  need not be copied */
#ifdef EEAM
      EAM_P  (p,ind) = EAM_P  (p,p->n); 
      /* eam_dM  need not be copied */
#endif
#endif
#ifdef ADP
  /* adp_mu and adp_lambda need not be copied */
#endif
#ifdef CG
      CG_H(p,ind,X) = CG_H(p,p->n,X); 
      CG_H(p,ind,Y) = CG_H(p,p->n,Y); 
#ifndef TWOD
      CG_H(p,ind,Z) = CG_H(p,p->n,Z); 
#endif
      CG_G(p,ind,X) = CG_G(p,p->n,X); 
      CG_G(p,ind,Y) = CG_G(p,p->n,Y); 
#ifndef TWOD
      CG_G(p,ind,Z) = CG_G(p,p->n,Z); 
#endif
      OLD_ORT(p,ind,X) = OLD_ORT(p,p->n,X); 
      OLD_ORT(p,ind,Y) = OLD_ORT(p,p->n,Y); 
#ifndef TWOD
      OLD_ORT(p,ind,Z) = OLD_ORT(p,p->n,Z); 
#endif
#endif /* CG */
#ifdef DAMP
      DAMPF(p,ind) = DAMPF(p,p->n); 
#endif
#ifdef DISLOC
      EPOT_REF(p,ind)   = EPOT_REF(p,p->n); 
      ORT_REF (p,ind,X) = ORT_REF (p,p->n,X); 
      ORT_REF (p,ind,Y) = ORT_REF (p,p->n,Y); 
#ifndef TWOD
      ORT_REF (p,ind,Z) = ORT_REF (p,p->n,Z); 
#endif
#endif
#ifdef AVPOS
      AV_EPOT(p,ind)   = AV_EPOT(p,p->n); 
      AV_POS (p,ind,X) = AV_POS (p,p->n,X); 
      AV_POS (p,ind,Y) = AV_POS (p,p->n,Y); 
      SHEET  (p,ind,X) = SHEET  (p,p->n,X); 
      SHEET  (p,ind,Y) = SHEET  (p,p->n,Y); 
#ifndef TWOD
      AV_POS (p,ind,Z) = AV_POS (p,p->n,Z); 
      SHEET  (p,ind,Z) = SHEET  (p,p->n,Z); 
#endif
#endif
#ifdef NNBR 
      NBANZ(p,ind)     = NBANZ(p,p->n); 
#endif
#ifdef REFPOS
      REF_POS(p,ind,X) = REF_POS(p,p->n,X);
      REF_POS(p,ind,Y) = REF_POS(p,p->n,Y);
#ifndef TWOD
      REF_POS(p,ind,Z) = REF_POS(p,p->n,Z);
#endif
#endif
#ifdef NVX
      HEATCOND(p,ind)  = HEATCOND(p,p->n); 
#endif
#ifdef STRESS_TENS
      PRESSTENS(p,ind,xx) = PRESSTENS(p,p->n,xx);   
      PRESSTENS(p,ind,yy) = PRESSTENS(p,p->n,yy);   
      PRESSTENS(p,ind,xy) = PRESSTENS(p,p->n,xy);   
#ifndef TWOD
      PRESSTENS(p,ind,zz) = PRESSTENS(p,p->n,zz);   
      PRESSTENS(p,ind,yz) = PRESSTENS(p,p->n,yz);   
      PRESSTENS(p,ind,zx) = PRESSTENS(p,p->n,zx);   
#endif
#endif /* STRESS_TENS */
#ifdef SHOCK
      PXAVG(p,ind)    = PXAVG(p,p->n); 
#endif
      IMPULS(p,ind,X) = IMPULS(p,p->n,X); 
      IMPULS(p,ind,Y) = IMPULS(p,p->n,Y); 
#ifndef TWOD
      IMPULS(p,ind,Z) = IMPULS(p,p->n,Z); 
#endif
      /* force need not be copied */
#ifdef COVALENT
      /* neighbor table need not be copied */
#endif
#ifdef NBLIST
      /* neighbor list reference positions are not copied */
#endif
#ifdef UNIAX
      ACHSE(p,ind,X) = ACHSE(p,p->n,X); 
      ACHSE(p,ind,Y) = ACHSE(p,p->n,Y); 
      ACHSE(p,ind,Z) = ACHSE(p,p->n,Z); 
      DREH_IMPULS(p,ind,X) = DREH_IMPULS(p,p->n,X); 
      DREH_IMPULS(p,ind,Y) = DREH_IMPULS(p,p->n,Y); 
      DREH_IMPULS(p,ind,Z) = DREH_IMPULS(p,p->n,Z); 
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
    to_cpu          = (int) (b->data[j++] + 0.1);
    ORT(input,0,X)  = b->data[j++];
    ORT(input,0,Y)  = b->data[j++];
#ifndef TWOD
    ORT(input,0,Z)  = b->data[j++];
#endif
#ifndef MONOLJ
    NUMMER(input,0) = b->data[j++];
#ifndef MONO
    SORTE (input,0) = b->data[j++];
#endif
    VSORTE(input,0) = b->data[j++];
    MASSE (input,0) = b->data[j++];
    POTENG(input,0) = b->data[j++];
#endif
#ifdef EAM2
    EAM_RHO(input,0) = b->data[j++];
    /* don't send eam_dF  */
#ifdef EEAM
    EAM_P  (input,0) = b->data[j++];
    /* don't send eam_dM  */
#endif
#endif
#ifdef ADP
  /* don't send adp_mu and adp_lambda */
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
#ifdef DAMP
    DAMPF(input,0) = b->data[j++];
#endif
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
#ifdef NNBR
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
#ifdef SHOCK
    PXAVG(input,0)        = b->data[j++];
#endif
    IMPULS(input,0,X)     = b->data[j++];
    IMPULS(input,0,Y)     = b->data[j++];
#ifndef TWOD
    IMPULS(input,0,Z)     = b->data[j++];
#endif
    /* don't send force */
#ifdef COVALENT
    /* don't send neighbor table */
#endif
#ifdef NBLIST
    /* neighbor list reference positions are not sent */
#endif
#ifdef UNIAX
    ACHSE(input,0,X) = b->data[j++];
    ACHSE(input,0,Y) = b->data[j++];
    ACHSE(input,0,Z) = b->data[j++];
    DREH_IMPULS(input,0,X) = b->data[j++];
    DREH_IMPULS(input,0,Y) = b->data[j++];
    DREH_IMPULS(input,0,Z) = b->data[j++];
    /* don't send dreh_moment */
#endif

    if (p==NULL) {  /* distribute atom into cell array */
#ifdef TWOD
      coord2=cell_coord( ORT(input,0,X), ORT(input,0,Y) );
#else
      coord2=cell_coord( ORT(input,0,X), ORT(input,0,Y), ORT(input,0,Z) );
#endif
      coord =local_cell_coord( coord2 );
      /* to_cpu = cpu_coord( coord2 ); determined by sendig CPU */
      if (to_cpu == myid) {
        minicell *to;
        to = PTR_VV(cell_array,coord,cell_dim);
        INSERT_ATOM(to, input, 0);
      }
    } else {  /* put atom into cell pointed by p */
#ifdef VEC
      /* we arrive here only in outdated socket communication code */
      error("not supported in vector mode");
#else
      move_atom(p, input, 0);
#endif
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
    realloc_msgbuf(to, to->n + from->n + BUFFER_SIZE_INC);
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
  int binc1, binc2, binc3;

  /* determine buffer size per atom */
  if (binc==0) {

    /* for communication to buffer cells */
    binc1 = DIM;     /* position */
#ifndef MONOLJ
    binc1++;         /* sorte */
#endif
#ifdef UNIAX
    binc1 += 3;      /* achse */
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
#ifdef NNBR
    binc2++;         /* nbanz */
#endif
#ifdef UNIAX
    binc2 += 3;      /* dreh_moment */
#endif

    /* communication of host electron density, adp_mu, adp_lambda */
#ifdef EAM2
#ifdef EEAM
    binc3 = 2;
#else 
    binc3 = 1;
#endif
#ifdef ADP
    binc3 += 9;
#endif
#endif /* EAM2 */

    /* one way or two ways */
#ifdef AR
    binc=MAX(binc1,binc2);
#else
    binc=binc1;
#endif
#ifdef EAM2
    binc=MAX(binc,binc3);
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
    alloc_msgbuf(&send_buf_east, size_east);
    alloc_msgbuf(&send_buf_west, size_east);
    alloc_msgbuf(&recv_buf_east, size_east);
    alloc_msgbuf(&recv_buf_west, size_east);
  }

  /* Allocate north/south buffers */
  if (size_north > send_buf_north.n_max) {
    alloc_msgbuf(&send_buf_north, size_north);
    alloc_msgbuf(&send_buf_south, size_north);
    alloc_msgbuf(&recv_buf_north, size_north);
    alloc_msgbuf(&recv_buf_south, size_north);
  }

#ifndef TWOD
  /* Allocate up/down buffers */
  if (size_up > send_buf_up.n_max) {
    alloc_msgbuf(&send_buf_up,   size_up);
    alloc_msgbuf(&send_buf_down, size_up);
    alloc_msgbuf(&recv_buf_up,   size_up);
    alloc_msgbuf(&recv_buf_down, size_up);
  }

#endif

#ifdef SHOCK
  if (0==dump_buf.n_max) {
    alloc_msgbuf(&dump_buf, BUFFER_SIZE_INC);
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

#ifdef USE_SOCKETS

/******************************************************************************
*
*  send_cell packs the data of all particles in the cell
*  into a buffer, and sends the buffer to the destination cpu.
*
******************************************************************************/

void send_cell(minicell *p, int to_cpu, int tag)
{
  int i;
  msgbuf *b;

  b = &send_buf_east;
  b->n = 0;

  for (i=0; i<p->n; i++) copy_one_atom(b, to_cpu, p, i, 0);
  MPI_Send(b->data, b->n, REAL, to_cpu, tag, cpugrid);
}


/******************************************************************************
*
*  recv_cell receives in a buffer the data of the particles in one cell,
*  and upacks the buffer into the cell pointed to by p.
*
******************************************************************************/

void recv_cell(minicell *p, int from_cpu, int tag)
{
  MPI_Status status;
  msgbuf *b;

  b = &recv_buf_east;

  /* check message size */
  MPI_Probe( from_cpu, tag, cpugrid, &status );
  MPI_Get_count( &status, REAL, &(b->n) );
  if (b->n_max < b->n) {
    alloc_msgbuf(b, b->n + BUFFER_SIZE_INC); 
  }

  /* upack data */
  p->n = 0;
  MPI_Recv(b->data, b->n, REAL, from_cpu, tag, cpugrid, &status);
  process_buffer(b,p);
}

#endif
