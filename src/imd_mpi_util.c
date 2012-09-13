
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2012 Institute for Theoretical and Applied Physics,
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
#ifdef BG
#undef INIT
#endif
#ifdef BGL
#include <rts.h>
#include <bglpersonality.h>
#endif
#ifdef BGP
#include <spi/kernel_interface.h>
#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>
#endif
#ifdef AFF
#define __USE_GNU
#include <sched.h>
void set_affinity_mask();
#endif

/******************************************************************************
*
* set up mpi
*
******************************************************************************/

#ifndef NEB
void init_mpi(void)
{
  /* Initialize MPI */
  MPI_Comm_size(MPI_COMM_WORLD,&num_cpus);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  if (0 == myid) { 
    printf("Starting up MPI with %d processes.\n", num_cpus);
#ifdef AFF
    set_affinity_mask();
#endif
#ifdef MPI2
    printf("Using MPI2\n");
#endif
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
#ifdef MPELOG
  MPE_Log_sync_clocks();
#ifdef NO_LLMPE
  MPE_Finish_log( progname );
#endif
#endif
  MPI_Finalize();                /* Shutdown */
}
#endif /* NEB */
/******************************************************************************
*
* initialize I/O parameters
*
******************************************************************************/

void init_io(void)
{

#ifdef BGL
#define   get_personality                rts_get_personality
#define   get_processor_id               rts_get_processor_id
#define   Personality                    BGLPersonality
#define   Personality_getLocationString  BGLPersonality_getLocationString
#define   Personality_numIONodes         BGLPersonality_numIONodes
#define   Personality_numPsets           BGLPersonality_numPsets
#define   Personality_numNodesInPset     BGLPersonality_numNodesInPset
#define   Personality_rankInPset         BGLPersonality_rankInPset
#define   Personality_psetNum            BGLPersonality_psetNum
#endif

#ifdef BGP
#define   get_personality                Kernel_GetPersonality
#define   get_processor_id               Kernel_PhysicalProcessorID
#define   Personality                    _BGP_Personality_t
#define   Personality_getLocationString  BGP_Personality_getLocationString
#define   Personality_numPsets           BGP_Personality_numIONodes
#define   Personality_numNodesInPset     BGP_Personality_psetSize
#define   Personality_rankInPset         BGP_Personality_rankInPset
#define   Personality_psetNum            BGP_Personality_psetNum
#endif

#ifdef BG
  Personality personality;
  int *tmp, i, node_config, mult=1;

  tmp     = (int *) malloc( num_cpus * sizeof(int) );
  io_grps = (int *) malloc( num_cpus * sizeof(int) );
  if ((NULL==tmp) || (NULL==io_grps)) 
    error("cannot allocate io_grps array");

  get_personality(&personality, sizeof(personality));
#ifdef BGL
  if (BGLPersonality_virtualNodeMode(&personality))    mult = 2;
  else                                                 mult = 1;
#endif
#ifdef BGP
  node_config = personality.Kernel_Config.ProcessConfig;
  if      (node_config == _BGP_PERS_PROCESSCONFIG_VNM) mult = 4;
  else if (node_config == _BGP_PERS_PROCESSCONFIG_2x2) mult = 2;
  else                                                 mult = 1;
#endif

  /* I/O group as a function of the rank */
  for (i=0; i<num_cpus; i++) tmp[i] = 0;
  tmp[myid] = Personality_psetNum(&personality);
  MPI_Allreduce(tmp, io_grps, num_cpus, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  /* input parameters */
  if (parallel_input==1) {
    n_inp_grps   = Personality_numPsets(&personality);
    my_inp_grp   = Personality_psetNum (&personality);
    inp_grp_size = Personality_numNodesInPset(&personality);
    inp_grp_size*= mult;
    my_inp_id    = 0;  while (my_inp_grp != io_grps[my_inp_id]) my_inp_id++;
  }
  else {
    n_inp_grps   = 1;
    my_inp_grp   = 0;
    my_inp_id    = 0;
    inp_grp_size = num_cpus;
  }

  /* output parameters */
  if (parallel_output==1) {
    n_out_grps   = Personality_numPsets(&personality);
    my_out_grp   = Personality_psetNum (&personality);
    out_grp_size = Personality_numNodesInPset(&personality);
    out_grp_size*= mult;
    my_out_id    = 0;  while (my_out_grp != io_grps[my_out_id]) my_out_id++;
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

#else /* not BG */

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

  /* find size of space an atom occupies in message buffer */
  /* we do this before having read or generated any atoms  */
  { 
    msgbuf b = {NULL, 0, 0};
    minicell c;
    alloc_msgbuf( &b, 256 );
    c.n_max = 0;
    ALLOC_MINICELL( &c, 1 );
    c.n = 1;
#ifdef VEC
    c.ind[0] = 0;
    atoms.n_max = 0;
    alloc_cell( &atoms, 1 );
    atoms.n = 1;
#endif
    copy_one_atom( &b, 0, &c, 0, 0);
    atom_size = b.n;
    free_msgbuf( &b );
    ALLOC_MINICELL( &c, 0 );
#ifdef VEC
    alloc_cell( &atoms, 0 );
#endif
  }

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
* copy an atom from a (mini)cell into a buffer
*
******************************************************************************/

void copy_atom_cell_buf(msgbuf *to, int to_cpu, cell *p, int ind )
{
  /* Check the parameters */
  if ((0 > ind) || (ind >= p->n)) {
    printf("%d: i %d n %d\n", myid, ind, p->n);
    error("copy_atom_cell_buf: index argument out of range.");
  }

  /* check if buffer is large enough; we can't just increase it, */ 
  /* as destination buffer wouldn't know about it                */
  if (to->n + atom_size > to->n_max) {
    error("buffer overflow in copy_atom_cell_buf");
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
#ifdef VARCHG
  to->data[ to->n++ ] = CHARGE(p,ind);
#endif
#ifdef DIPOLE 
  /* dp_E_stat, dp_E_ind and dp_p_stat are not sent */
/*   to->data[ to->n++ ] = DP_P_STAT(p,ind,X); */
/*   to->data[ to->n++ ] = DP_P_STAT(p,ind,Y); */
/*   to->data[ to->n++ ] = DP_P_STAT(p,ind,Z); */
  to->data[ to->n++ ] = DP_P_IND(p,ind,X);
  to->data[ to->n++ ] = DP_P_IND(p,ind,Y);
#ifndef TWOD
  to->data[ to->n++ ] = DP_P_IND(p,ind,Z);
#endif
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
#ifdef HC
  to->data[ to->n++ ] = HCAVENG(p,ind);
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
#ifdef AVPOS
  to->data[ to->n++ ] = AVPRESSTENS(p,ind,xx);   
  to->data[ to->n++ ] = AVPRESSTENS(p,ind,yy);   
  to->data[ to->n++ ] = AVPRESSTENS(p,ind,xy);   
#ifndef TWOD
  to->data[ to->n++ ] = AVPRESSTENS(p,ind,zz);   
  to->data[ to->n++ ] = AVPRESSTENS(p,ind,yz);   
  to->data[ to->n++ ] = AVPRESSTENS(p,ind,zx);   
#endif
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
  copy_atom_cell_buf(to, to_cpu, p, ind);

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
    if (index < from->n) from->ind[index] = from->ind[from->n];
#endif

    /* copy data of last atom into empty slot */
    if (ind < p->n) copy_atom_cell_cell(p, ind, p, p->n);

  } /* delete or not delete */
}

/******************************************************************************
*
*  Append atom from message buffer to cell
*
******************************************************************************/

void copy_atom_buf_cell(minicell *p, msgbuf *b, int start)
{
  int  ind, j = start + 1;  /* the first entry is the CPU number */
  cell *to;

#ifdef VEC
  if (p->n >= p->n_max) alloc_minicell(p,p->n_max+incrsz);
  p->ind[p->n] = atoms.n;
  if (atoms.n >= atoms.n_max) alloc_cell(&atoms, 2*atoms.n_max);
#ifdef MPI
  atoms.ind[atoms.n] = p->n; /* pointer from atom to index in its minicell */
#endif
  ind = atoms.n++;
  p->n++;
  to = &atoms;
#else
  to = p;
  if (to->n >= to->n_max) alloc_cell(to,to->n_max+incrsz);
  ind = to->n++;
#endif

  ORT(to,ind,X)  = b->data[j++];
  ORT(to,ind,Y)  = b->data[j++];
#ifndef TWOD
  ORT(to,ind,Z)  = b->data[j++];
#endif
#ifndef MONOLJ
  NUMMER(to,ind) = b->data[j++];
#ifndef MONO
  SORTE (to,ind) = b->data[j++];
#endif
  VSORTE(to,ind) = b->data[j++];
  MASSE (to,ind) = b->data[j++];
  POTENG(to,ind) = b->data[j++];
#endif
#ifdef EAM2
  EAM_RHO(to,ind) = b->data[j++];
  /* don't send eam_dF  */
#ifdef EEAM
  EAM_P  (to,ind) = b->data[j++];
  /* don't send eam_dM  */
#endif
#endif
#ifdef ADP
  /* don't send adp_mu and adp_lambda */
#endif
#ifdef VARCHG
  CHARGE(to,ind)     = b->data[j++];
#endif
#ifdef DIPOLE
  /* don't send p_stat, E_stat, E_ind */
  DP_P_IND(to,ind,X) = b->data[j++];
  DP_P_IND(to,ind,Y) = b->data[j++];
#ifndef TWOD
  DP_P_IND(to,ind,Z) = b->data[j++];
#endif
#endif /* DIPOLE */
#ifdef CG
  CG_H(to,ind,X) = b->data[j++];
  CG_H(to,ind,Y) = b->data[j++];
#ifndef TWOD
  CG_H(to,ind,Z) = b->data[j++];
#endif
  CG_G(to,ind,X) = b->data[j++];
  CG_G(to,ind,Y) = b->data[j++];
#ifndef TWOD
  CG_G(to,ind,Z) = b->data[j++];
#endif
  OLD_ORT(to,ind,X) = b->data[j++];
  OLD_ORT(to,ind,Y) = b->data[j++];
#ifndef TWOD
  OLD_ORT(to,ind,Z) = b->data[j++];
#endif
#endif /* CG */
#ifdef DAMP
  DAMPF(to,ind) = b->data[j++];
#endif
#ifdef DISLOC
  EPOT_REF(to,ind)   = b->data[j++];
  ORT_REF (to,ind,X) = b->data[j++];
  ORT_REF (to,ind,Y) = b->data[j++];
#ifndef TWOD
  ORT_REF (to,ind,Z) = b->data[j++];
#endif
#endif
#ifdef AVPOS
  AV_EPOT(to,ind)    = b->data[j++];
  AV_POS (to,ind,X)  = b->data[j++];
  AV_POS (to,ind,Y)  = b->data[j++];
  SHEET  (to,ind,X)  = b->data[j++];
  SHEET  (to,ind,Y)  = b->data[j++];
#ifndef TWOD
  AV_POS (to,ind,Z)  = b->data[j++];
  SHEET  (to,ind,Z)  = b->data[j++];
#endif
#endif
#ifdef NNBR
  NBANZ(to,ind)      = b->data[j++];
#endif
#ifdef REFPOS
  REF_POS(to,ind,X)  = b->data[j++];
  REF_POS(to,ind,Y)  = b->data[j++];
#ifndef TWOD
  REF_POS(to,ind,Z)  = b->data[j++];
#endif
#endif
#ifdef HC
  HCAVENG(to,ind)    = b->data[j++];
#endif
#ifdef STRESS_TENS
  PRESSTENS(to,ind,xx) = b->data[j++];   
  PRESSTENS(to,ind,yy) = b->data[j++];   
  PRESSTENS(to,ind,xy) = b->data[j++];   
#ifndef TWOD
  PRESSTENS(to,ind,zz) = b->data[j++];   
  PRESSTENS(to,ind,yz) = b->data[j++];   
  PRESSTENS(to,ind,zx) = b->data[j++];   
#endif
#ifdef AVPOS
  AVPRESSTENS(to,ind,xx) = b->data[j++];   
  AVPRESSTENS(to,ind,yy) = b->data[j++];   
  AVPRESSTENS(to,ind,xy) = b->data[j++];   
#ifndef TWOD
  AVPRESSTENS(to,ind,zz) = b->data[j++];   
  AVPRESSTENS(to,ind,yz) = b->data[j++];   
  AVPRESSTENS(to,ind,zx) = b->data[j++];   
#endif
#endif
#endif /* STRESS_TENS */
#ifdef SHOCK
  PXAVG(to,ind)        = b->data[j++];
#endif
  IMPULS(to,ind,X)     = b->data[j++];
  IMPULS(to,ind,Y)     = b->data[j++];
#ifndef TWOD
  IMPULS(to,ind,Z)     = b->data[j++];
#endif
  /* don't send force */
#ifdef COVALENT
  /* don't send neighbor table */
#endif
#ifdef NBLIST
  /* neighbor list reference positions are not sent */
#endif
#ifdef UNIAX
  ACHSE(to,ind,X) = b->data[j++];
  ACHSE(to,ind,Y) = b->data[j++];
  ACHSE(to,ind,Z) = b->data[j++];
  DREH_IMPULS(to,ind,X) = b->data[j++];
  DREH_IMPULS(to,ind,Y) = b->data[j++];
  DREH_IMPULS(to,ind,Z) = b->data[j++];
  /* don't send dreh_moment */
#endif

}

/******************************************************************************
*
*  unpack atoms from message buffer
*
******************************************************************************/

void process_buffer(msgbuf *b)
{
  int i;
  minicell *to;
  ivektor coord;

  for (i=0; i<b->n; i+=atom_size) {
    if (myid != (int) (b->data[i] + 0.1)) continue;
#ifdef TWOD
    coord = cell_coord( b->data[i+1], b->data[i+2] );
#else
    coord = cell_coord( b->data[i+1], b->data[i+2], b->data[i+3] );
#endif
    coord = local_cell_coord( coord );
    to = PTR_VV(cell_array, coord, cell_dim);
    copy_atom_buf_cell(to, b, i);
  }
}

/******************************************************************************
*
* copy atoms from one message buffer to another
*
******************************************************************************/

void copy_atoms_buf(msgbuf *to, msgbuf *from)
{
  int i, j;

  /* copy only atoms which have to go elsewhere */
  for (i=0; i<from->n; i+=atom_size) {
    if (myid != from->data[i]) {
      if (to->n_max < to->n + atom_size) 
        error("buffer overflow in copy_atoms_buf");
      for (j=0; j<atom_size; j++) to->data[to->n++] = from->data[i+j];
    }
  }
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
  int k, tmp=0;
  int size_east;
  int size_north;
  int size_up;
  int binc1, binc2, binc3, binc4;

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
#ifdef VARCHG
    binc1++;         /* charge */   
#endif

    /* for communication from buffer cells */
    binc2 = DIM;     /* force */
#ifndef MONOLJ
    binc2++;         /* pot_eng */
#endif
#ifdef HC
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
    /* communication of induced dipoles etc */
#ifdef DIPOLE
    binc4 = 4*DIM; 		/* 4 vector fields */
#endif

    /* one way or two ways */
#ifdef AR
    binc=MAX(binc1,binc2);
#else
    binc=binc1;
#endif
#ifdef EAM2
    binc=MAX(binc,binc3);
#endif
#ifdef DIPOLE
    binc=MAX(binc,binc4);
#endif
#ifdef NYETENSOR
/* Nye tensors are not included in standard communication, but are send independently in their own subroutine */
/* Buffer must be large enough for this */
    binc=MAX(binc,19); 
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
  largest_cell = (int) largest_cell * msgbuf_size;

#ifndef TWOD
  size_east  = largest_cell * cell_dim.y * cell_dim.z * binc;
  size_north = largest_cell * cell_dim.x * cell_dim.z * binc;
  size_up    = largest_cell * cell_dim.x * cell_dim.y * binc;
#ifdef DEBUG
  if (1==myid) 
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
    alloc_msgbuf(&dump_buf, 1024);
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

#ifdef AFF

/******************************************************************************
*
* set processor affinity
*
******************************************************************************/

void set_affinity_mask()
{
  int  n=0, i, ht, ncores=1, siblings=1;
  char *res, *str, line[1024];
  FILE *inp;

  /* count number of CPUs */
  inp = fopen("/proc/cpuinfo", "r");
  if (NULL==inp) return;
  while (!feof(inp)) {
    res=fgets(line,1024,inp);
    if (NULL==res) break;
    if (strncmp(line, "processor", 9) == 0) n++;
    if (strncmp(line, "cpu cores", 8) == 0) {
      str = strstr(line,":") + 1;
      sscanf(str, "%d", &ncores);
    }
    if (strncmp(line, "siblings", 8) == 0) {
      str = strstr(line,":") + 1;
      sscanf(str, "%d", &siblings);
    }
  }
  fclose(inp);

  /* if hyperthreading, confine processes to even CPUs */
  ht = siblings / ncores; 
  if (ht > 1) {
    cpu_set_t my_cpu_mask;
    CPU_ZERO(&my_cpu_mask);
    for (i=0; i<n; i+=ht) CPU_SET(i,&my_cpu_mask);
    sched_setaffinity(0,sizeof(cpu_set_t),&my_cpu_mask);
  }
}

#endif
