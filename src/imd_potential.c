
/*****************************************************************************
*
*  Routines for reading and accessing tabulated potentials and functions
*  Potentials are stored in a data structure of type pot_table_t.
*  When reading a new potential, cellsz must be set to the maximum
*  of cellsz and the maximum of the new squared potential cutoffs.
*  If necessary, cellsz must then be broadcasted again.
*
*  $Revision$
*  $Date$
*
******************************************************************************/

#include "imd.h"


/*****************************************************************************
*
*  read potential in first format: each line contains
*
*  r**2 V00 V01 V02 ... V10 V11 V12 ... VNN
*
*  N is the number of different atom types
*
*  Note that it is assumed that Vij == Vji and that the r**2 are aequidistant.
*
******************************************************************************/

void read_pot_table1( pot_table_t *pt, char *filename )
{
  FILE *infile;
  int i, k;
  int size, tablesize, npot=0;
  real val, numstep, delta;
  real r2, r2_start, r2_step;
  str255 msg;

  /* allocate info block of function table */
  size = ntypes*ntypes;
  pt->maxsteps = 0;
  pt->begin    = (real *) malloc(size*sizeof(real));
  pt->end      = (real *) malloc(size*sizeof(real));
  pt->step     = (real *) malloc(size*sizeof(real));
  pt->invstep  = (real *) malloc(size*sizeof(real));
  if ((pt->begin   == NULL) || (pt->end == NULL) || (pt->step == NULL) || 
      (pt->invstep == NULL)) {
    sprintf(msg,"Cannot allocate info block for function table %s.",filename);
    error(msg);
  }

#ifdef MPI
  /* read table only on master processor? */
  if ((0==myid) || (1==parallel_input)) {
#endif

    infile = fopen(filename,"r");
    if (NULL==infile) {
      sprintf(msg,"Cannot open file %s.",filename);
      error(msg);
    }

    /* allocate the function table */
    pt->maxsteps = PSTEP;
    tablesize = size * pt->maxsteps;
    pt->table = (real *) malloc(tablesize*sizeof(real));
    if (NULL==pt->table) {
      sprintf(msg,"Cannot allocate memory for function table %s.",filename);
      error(msg);
    }

    /* input loop */
    while (!feof(infile)) {

      /* still some space left? */ 
      if (((npot%PSTEP) == 0) && (npot>0)) {
        pt->maxsteps += PSTEP;
        tablesize = size * pt->maxsteps;
        pt->table = (real *) realloc(pt->table, tablesize*sizeof(real));
        if (NULL==pt->table) {
          sprintf(msg,"Cannot extend memory for function table %s.",filename);
          error(msg);
        }
      }

      /*  read in potential */
      if ( 1 != fscanf(infile,"%lf",&r2) ) break;
      if (npot==0) r2_start = r2;  /* catch first value */
      for (i=0; i<size; ++i) {
	if (( 1 != fscanf(infile,"%lf", &val)) && (myid==0)) 
           error("Line incomplete in potential file.");
	*PTR_2D(pt->table,npot,i,pt->maxsteps,size) = val;
        if (val!=0.0) pt->end[i] = r2; /* catch last non-zero value */
      }
      ++npot;
    }

    fclose(infile);

    r2_step = (r2 - r2_start) / npot;

    if (0==myid) {
      printf("Read potential %s with %d lines.\n",filename,npot);
      printf("Starts at r2_start: %f, r_start: %f\n",r2_start,sqrt(r2_start));
      printf("Ends at r2_end:     %f, r_end:   %f\n",r2,      sqrt(r2));
      printf("Step is r2_step:    %f\n",r2_step);
    }

    /* fill info block, and shift potential to zero */
    for (i=0; i<size; ++i) {
      pt->begin[i] = r2_start;
      pt->step[i] = r2_step;
      pt->invstep[i] = 1.0 / r2_step;
      delta = *PTR_2D(pt->table,(npot-1),i,pt->maxsteps,size);
      if (delta!=0.0) {
        if (0==myid)
          printf("Potential %1d%1d shifted by %f\n",
                 (i/ntypes),(i%ntypes),delta);
        for (k=0; k<npot; ++k) *PTR_2D(pt->table,k,i,pt->table,size) -= delta;
      } else {
        pt->end[i] += r2_step;
      }
      cellsz = MAX(cellsz,pt->end[i]);
    }
    if (0==myid) printf("\n");

    /* The interpolation uses k+1 and k+2, so we add zeros at end of table */
    for(k=1; k<=5; ++k) {
      /* still some space left? */ 
      if (((npot%PSTEP) == 0) && (npot>0)) {
        pt->maxsteps += PSTEP;
        tablesize = size * pt->maxsteps;
        pt->table = (real *) realloc(pt->table, tablesize*sizeof(real));
        if (NULL==pt->table) {
          sprintf(msg,"Cannot extend memory for function table %s.",filename);
          error(msg);
        }
      }
      for (i=0; i<size; ++i)
	*PTR_2D(pt->table,npot,i,pt->table,size) = 0.0;
      ++npot;
    }

#ifdef MPI
  }
  if (0==parallel_input) {
    /* Broadcast table to other CPUs */
    MPI_Bcast( &(pt->maxsteps),  1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast( pt->begin,    size, MPI_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast( pt->end,      size, MPI_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast( pt->step,     size, MPI_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast( pt->invstep,  size, MPI_REAL, 0, MPI_COMM_WORLD);
    tablesize = pt->maxsteps * size;
    if (0 != myid) {
      pt->table = (real *) malloc(tablesize*sizeof(real));
      if (NULL==pt->table)
        error("Cannot allocate memory for function table");
    }
    MPI_Bcast( pt->table, tablesize, MPI_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast( &cellsz,           1, MPI_REAL, 0, MPI_COMM_WORLD);
  }
#endif

}


/*****************************************************************************
*
*  read potential in second format: at the beginning ntypes*ntypes times 
*  a line of the form
*
*  r_begin r_end r_step,
*  
*  then the values of the potential (one per line), first those 
*  for atom pair  00, then an empty line (for gnuplot), then 01 and so on.
*
*  Note that it is assumed that Vij == Vji and that the r**2 are aequidistant.
*
******************************************************************************/

void read_pot_table2( pot_table_t *pt, char *filename, int mode )
{
  FILE *infile;
  int i, k, *len;
  int size, tablesize;
  real val, numstep;
  str255 msg;

  /* allocate info block of function table */
  size = ntypes;
  if (mode==POTTAB) size *= ntypes;
  pt->maxsteps = 0;
  pt->begin    = (real *) malloc(size*sizeof(real));
  pt->end      = (real *) malloc(size*sizeof(real));
  pt->step     = (real *) malloc(size*sizeof(real));
  pt->invstep  = (real *) malloc(size*sizeof(real));
  len          = (int  *) malloc(size*sizeof(real));
  if ((pt->begin   == NULL) || (pt->end == NULL) || (pt->step == NULL) || 
      (pt->invstep == NULL) || (len == NULL)) {
    sprintf(msg,"Cannot allocate info block for function table %s.",filename);
    error(msg);
  }

#ifdef MPI
  /* read table only on master processor? */
  if ((0==myid) || (1==parallel_input)) {
#endif

    infile = fopen(filename,"r");
    if (NULL==infile) {
      sprintf(msg,"Cannot open file %s.",filename);
      error(msg);
    }

    /* read the info block of the function table */
    for(i=0; i<size; i++) {
      if ((0==myid) && (3!=fscanf(infile, "%lf %lf %lf",
                              &pt->begin[i], &pt->end[i], &pt->step[i]))) { 
        sprintf(msg, "Info line in %s corrupt.", filename);
        error(msg);
      }
      cellsz = MAX(cellsz,pt->end[i]);
      pt->invstep[i] = 1.0 / pt->step[i];
      numstep        = 1 + (pt->end[i] - pt->begin[i]) / pt->step[i];
      len[i]         = (int) (numstep+0.5);  
      pt->maxsteps   = MAX(pt->maxsteps, len[i]);

      /* some security against rounding errors */
      if ((fabs(len[i] - numstep) >= 0.1) && (0==myid)) {
        fprintf(stderr,"WARNING: numstep = %f rounded to %d in file %s.\n",
                numstep, len[i], filename);
      }
    }

    /* allocate the function table */
    tablesize = size * pt->maxsteps;
    pt->table = (real *) malloc(tablesize*sizeof(real));
    if (NULL==pt->table) {
      sprintf(msg,"Cannot allocate memory for function table %s.",filename);
      error(msg);
    }

    /* input loop */
    for (i=0; i<size; i++)
      for (k=0; k<len[i]; k++) {
        if ((1 != fscanf(infile,"%lf", &val)) && (0==myid)) {
          sprintf(msg,"wrong format in file %s.",filename);
          error(msg);
        }
        *PTR_2D(pt->table,k,i,pt->maxsteps,size) = val;
      }

    fclose(infile);

    if (0==myid) {
      if (mode==FUNTAB) {
        printf("Read tabulated function %s for %d atoms types.\n",
               filename,size);
      } else {
        printf("Read tabulated function %s for %d pairs of atoms types.\n",
               filename,size);
      }
      printf("Maximal length of table is %d.\n",pt->maxsteps);
    }

#ifdef MPI
  }
  if (0==parallel_input) {
    /* Broadcast table to other CPUs */
    MPI_Bcast( &(pt->maxsteps), 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast( pt->begin,    size, MPI_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast( pt->end,      size, MPI_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast( pt->step,     size, MPI_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast( pt->invstep,  size, MPI_REAL, 0, MPI_COMM_WORLD);
    tablesize = pt->maxsteps * size;
    if (0 != myid) {
      pt->table = (real *) malloc(tablesize*sizeof(real));
      if (NULL==pt->table)
        error("Cannot allocate memory for function table");
    }
    MPI_Bcast( pt->table, tablesize, MPI_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast( &cellsz,           1, MPI_REAL, 0, MPI_COMM_WORLD);
  }
#endif

}


/*****************************************************************************
*
*  Evaluate monolj potential 
*
******************************************************************************/

int pair_int_monolj(real *pot, real *grad, real r2)
{
  real sig_d_rad2, sig_d_rad6, sig_d_rad12;

  sig_d_rad2  = 2.0 / r2;
  sig_d_rad6  = sig_d_rad2 * sig_d_rad2 * sig_d_rad2;
  sig_d_rad12 = sig_d_rad6 * sig_d_rad6;

  *grad = - 6 * sig_d_rad2 * ( sig_d_rad12 - sig_d_rad6 );
  *pot  = sig_d_rad12 - 2.0 * sig_d_rad6 - monolj_shift;
  return 0;
}


/*****************************************************************************
*
*  Evaluate potential table with quadratic interpolation. 
*  Returned are potential value V and gradient (dV/dr)/r.
*  col is p_typ * ntypes + q_typ
*
******************************************************************************/

int pair_int2(real *pot, real *grad, pot_table_t *pt, int col, real r2)
{
  real r2a, istep, chi, p0, p1, p2, dv, d2v;
  real *ptr;
  int  k, inc, is_short=0;

  /* check for distances shorter than minimal distance in table */
  r2a = r2 - pt->begin[col];
  if (r2a < 0) {
    r2a   = 0;
    is_short = 1;
  }

  /* Indices into potential table */
  istep = pt->invstep[col];
  k     = (int) (r2a * istep);
  chi   = (r2a - k * pt->step[col]) * istep;

  /* A single access to the potential table involves two multiplications 
     We use a intermediate pointer to avoid this as much as possible.
     Note: This relies on layout of the pot-table in memory!!! */

  inc = ntypes * ntypes;
  ptr = PTR_2D(pt->table, k, col, pt->maxsteps, inc);
  p0  = *ptr; ptr += inc;
  p1  = *ptr; ptr += inc;
  p2  = *ptr;
  dv  = p1 - p0;
  d2v = p2 - 2 * p1 + p0;

  /* norm of gradient and potential energy */
  *grad = 2 * istep * (dv + (chi - 0.5) * d2v);
  *pot  = p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;
  return is_short;

}


/*****************************************************************************
*
*  Evaluate potential table with quartic interpolation. 
*  Returned are potential value and gradient.
*  col is p_typ * ntypes + q_typ
*
******************************************************************************/

int pair_int4(real *pot, real *grad, pot_table_t *pt, int col, real r2)
{
  return 0;
}


/*****************************************************************************
*
*  Evaluate function table with quadratic interpolation. 
*
******************************************************************************/

int table_func2(real *val, pot_table_t *pt, int col, real x)
{
  return 0;
}


/*****************************************************************************
*
*  Evaluate function table with quartic interpolation. 
*
******************************************************************************/

int table_func4(real *val, pot_table_t *pt, int col, real x)
{
  return 0;
}
