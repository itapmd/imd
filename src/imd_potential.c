
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/*****************************************************************************
*
*  Routines for reading tabulated potentials and functions.
*  The tables are accessed with the macros in potaccess.h.
*  The acces functions in this file are currently not used.
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
*  read potential table; choose format according to header
*
*****************************************************************************/

void read_pot_table( pot_table_t *pt, char *filename, int ncols )
{
  FILE *infile;
  char buffer[1024], msg[255];
  char *token, *res;
  int  have_header=0, have_format=0, end_header;
  int  size=ncols, tablesize, npot=0;
  int  format=DEFAULT_POTFILE_TYPE;  /* 2 for EAM2, 1 otherwise */
  int  i, k;

  /* read header only on master processor? */
  if ((0==myid) || (1==parallel_input)) {

    /* open file */
    infile = fopen(filename,"r");
    if (NULL == infile) {
      sprintf(msg,"Could not open file %s\n",filename);
      error(msg);
    }

    /* read the header */
    do {
      /* read one line */
      res=fgets(buffer,1024,infile);
      if (NULL == res) {
        sprintf(msg,"Unexpected end of file in %s",filename);
        error(msg);
      }

      /* see if it is a header line */
      if (buffer[0]=='#') {
        have_header = 1;
        /* stop after last header line */
        end_header = (buffer[1]=='E');
        /* see if it is the format line */
        if (buffer[1]=='F') {
          /* format complete? */
          if (2!=sscanf( (const char*)(buffer+2), "%d%d", &format, &size )) {
            sprintf(msg,"Corrupted format header line in file %s",filename);
            error(msg);
          }
          /* right number of columns? */
          if (size!=ncols) {
            sprintf(msg,"Wrong number of data columns in file %s",filename);
            error(msg);
          }
          /* recognized format? */
          if ((format!=1) && (format!=2)) {
            sprintf(msg,"Unrecognized format specified for file %s",filename);
            error(msg);
          }
          have_format=1;
	}
      } else if (have_header) { 
        /* header does not end properly */
        sprintf(msg,"Corrupted header in file %s",filename);
        error(msg);
      } else {
        /* we have no header, stop reading further */
	end_header=1;
      }
    } while (!end_header);

    /* did we have a format in the header */
    if ((have_header) && (!have_format)) {
      sprintf(msg,"Format not specified in header of file %s",filename);
      error(msg);
    }

    /* rewind if there was no header */
    if (!have_header) rewind(infile);

    /* warn if we have no header */
    if ((0==myid) && (!have_header)) {
      sprintf(msg,"File %s has no header",filename);
      warning(msg);
    }

  } /* have read header */

  /* allocate info block of function table */
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

  /* catch the case where potential is identically zero */
  for (i=0; i<size; ++i) {
    pt->end[i] = 0.0;
  }

  /* read table only on master processor? */
  if ((0==myid) || (1==parallel_input)) {
    if (format==1) read_pot_table1(pt, size, filename, infile);
    if (format==2) read_pot_table2(pt, size, filename, infile);
  }
  fclose(infile);

#ifdef MPI
  if (0==parallel_input) {
    /* Broadcast table to other CPUs */
    MPI_Bcast( &(pt->maxsteps), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( pt->begin,    size, REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast( pt->end,      size, REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast( pt->step,     size, REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast( pt->invstep,  size, REAL, 0, MPI_COMM_WORLD);
    tablesize = pt->maxsteps * size;
    if (0 != myid) {
      pt->table = (real *) malloc(tablesize*sizeof(real));
      if (NULL==pt->table)
        error("Cannot allocate memory for function table");
    }
    MPI_Bcast( pt->table, tablesize, REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast( &cellsz,           1, REAL, 0, MPI_COMM_WORLD);
  }
#endif

}


/*****************************************************************************
*
*  read potential in first format: each line contains
*
*  r**2 V00 V01 V02 ... V10 V11 V12 ... VNN
*
*  N is the number of different atom types
*
*  Note that it is assumed that the r**2 are aequidistant.
*
******************************************************************************/

void read_pot_table1(pot_table_t *pt, int ncols, char *filename, FILE *infile)
{
  int i, k;
  int tablesize, npot=0;
  real val, numstep, delta;
  real r2, r2_start, r2_step;
  str255 msg;

  /* allocate the function table */
  pt->maxsteps = PSTEP;
  tablesize = ncols * pt->maxsteps;
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
      tablesize = ncols * pt->maxsteps;
      pt->table = (real *) realloc(pt->table, tablesize*sizeof(real));
      if (NULL==pt->table) {
        sprintf(msg,"Cannot extend memory for function table %s.",filename);
        error(msg);
      }
    }

    /*  read in potential */
    if ( 1 != fscanf(infile,"%lf",&r2) ) break;
    if (npot==0) r2_start = r2;  /* catch first value */
    for (i=0; i<ncols; ++i) {
      if (( 1 != fscanf(infile,"%lf", &val)) && (myid==0)) 
        error("Line incomplete in potential file.");
      *PTR_2D(pt->table,npot,i,pt->maxsteps,ncols) = val;
      if (val!=0.0) pt->end[i] = r2; /* catch last non-zero value */
    }
    ++npot;
  }

  r2_step = (r2 - r2_start) / (npot-1);

  if (0==myid) {
    printf("Read potential %s with %d lines.\n",filename,npot);
    printf("Starts at r2_start: %f, r_start: %f\n",r2_start,sqrt(r2_start));
    printf("Ends at r2_end:     %f, r_end:   %f\n",r2,      sqrt(r2));
    printf("Step is r2_step:    %f\n",r2_step);
  }

  /* fill info block, and shift potential to zero */
  for (i=0; i<ncols; ++i) {
    pt->begin[i] = r2_start;
    pt->step[i] = r2_step;
    pt->invstep[i] = 1.0 / r2_step;
    delta = *PTR_2D(pt->table,(npot-1),i,pt->maxsteps,ncols);
    if (delta!=0.0) {
      if (0==myid)
        printf("Potential %1d%1d shifted by %f\n",
               (i/ntypes),(i%ntypes),delta);
      for (k=0; k<npot; ++k) *PTR_2D(pt->table,k,i,pt->table,ncols) -= delta;
    } else {
      pt->end[i] += r2_step;
    }
    if (ncols==ntypes*ntypes) cellsz = MAX(cellsz,pt->end[i]);
  }
  if (0==myid) printf("\n");

  /* The interpolation uses k+1 and k+2, so we add zeros at end of table */
  for (k=1; k<=5; ++k) {
    /* still some space left? */ 
    if (((npot%PSTEP) == 0) && (npot>0)) {
      pt->maxsteps += PSTEP;
      tablesize = ncols * pt->maxsteps;
      pt->table = (real *) realloc(pt->table, tablesize*sizeof(real));
      if (NULL==pt->table) {
        sprintf(msg,"Cannot extend memory for function table %s.",filename);
        error(msg);
      }
    }
    for (i=0; i<ncols; ++i)
      *PTR_2D(pt->table,npot,i,pt->table,ncols) = 0.0;
    ++npot;
  }
}


/*****************************************************************************
*
*  read potential in second format: at the beginning <ncols> times 
*  a line of the form
*
*  r_begin r_end r_step,
*  
*  then the values of the potential (one per line), first those 
*  for atom pair  00, then an empty line (for gnuplot), then 01 and so on.
*  Analogously, if there is only one column per atom type.
*
*  Note that it is assumed that the r**2 are aequidistant.
*
******************************************************************************/

void read_pot_table2(pot_table_t *pt, int ncols, char *filename, FILE *infile)
{
  int i, k, *len;
  int tablesize;
  real val, numstep;
  str255 msg;

  len = (int  *) malloc(ncols * sizeof(real));
  if (len==NULL) error("allocation failed in read_pot_table");

  /* read the info block of the function table */
  for(i=0; i<ncols; i++) {
    if (3 != fscanf(infile, "%lf %lf %lf",
                  &pt->begin[i], &pt->end[i], &pt->step[i])) {
      if (0==myid) { 
        sprintf(msg, "Info line in %s corrupt.", filename);
        error(msg);
      }
    }
    if (ncols==ntypes*ntypes) cellsz = MAX(cellsz,pt->end[i]);
    pt->invstep[i] = 1.0 / pt->step[i];
    numstep        = 1 + (pt->end[i] - pt->begin[i]) / pt->step[i];
    len[i]         = (int) (numstep+0.5);  
    pt->maxsteps   = MAX(pt->maxsteps, len[i]);

    /* some security against rounding errors */
    if ((fabs(len[i] - numstep) >= 0.1) && (0==myid)) {
      char msg[255];
      sprintf(msg,"numstep = %f rounded to %d in file %s.",
              numstep, len[i], filename);
      warning(msg);
    }
  }

  /* allocate the function table */
  /* allow some extra values at the end for interpolation */
  tablesize = ncols * (pt->maxsteps+3);
  pt->table = (real *) malloc(tablesize * sizeof(real));
  if (NULL==pt->table) {
    sprintf(msg,"Cannot allocate memory for function table %s.",filename);
    error(msg);
  }

  /* input loop */
  for (i=0; i<ncols; i++) {
    for (k=0; k<len[i]; k++) {
      if (1 != fscanf(infile,"%lf", &val)) {
        if (0==myid) {
          sprintf(msg, "wrong format in file %s.", filename);
          error(msg);
        }
      }
      *PTR_2D(pt->table,k,i,pt->maxsteps,ncols) = val;
    }
    /* make some copies of the last value for interpolation */
    for (k=len[i]; k<len[i]+3; k++)
      *PTR_2D(pt->table,k,i,pt->maxsteps,ncols) = val;
  }

  if (0==myid) {
    if (ncols==ntypes) {
      printf("Read tabulated function %s for %d atoms types.\n",
             filename,ncols);
    } else {
      printf("Read tabulated function %s for %d pairs of atoms types.\n",
             filename,ncols);
    }
    printf("Maximal length of table is %d.\n",pt->maxsteps);
  }
}

/*****************************************************************************
*
*  Evaluate monolj potential 
*
******************************************************************************/

void pair_int_monolj(real *pot, real *grad, real r2)
{
  real sig_d_rad2, sig_d_rad6, sig_d_rad12;

  sig_d_rad2  = 2.0 / r2;
  sig_d_rad6  = sig_d_rad2 * sig_d_rad2 * sig_d_rad2;
  sig_d_rad12 = sig_d_rad6 * sig_d_rad6;

  *grad = - 6 * sig_d_rad2 * ( sig_d_rad12 - sig_d_rad6 );
  *pot  = sig_d_rad12 - 2.0 * sig_d_rad6 - monolj_shift;
}

#ifdef STIWEB

/*****************************************************************************
*
*  Evaluate pair potential for Stillinger-Weber potential 
*
******************************************************************************/

void pair_int_stiweb(real *pot, real *grad, int p_typ, int q_typ, real r2)
{
  real radius, phi_r, phi_a, inv_c, inv_r, f_cut;

  radius = sqrt(r2);
  phi_r  =   sw_a[p_typ][q_typ] * pow( radius, - sw_p[p_typ][q_typ] );
  phi_a  = - sw_b[p_typ][q_typ] * pow( radius, - sw_q[p_typ][q_typ] );
  inv_c  = 1.0 / ( radius - sw_a1[p_typ][q_typ] );
  inv_r  = 1.0 / radius;
  f_cut  = exp( sw_de[p_typ][q_typ] * inv_c );

  *pot  = ( phi_r + phi_a ) * f_cut;
  *grad = ( - *pot * sw_de[p_typ][q_typ] * inv_c * inv_c 
	    - f_cut * inv_r * ( sw_p[p_typ][q_typ] * phi_r 
			      + sw_q[p_typ][q_typ] * phi_a ) ) * inv_r; 

}

#endif

/*****************************************************************************
*
*  Evaluate potential table with quadratic interpolation. 
*  Returns the potential value and twice the derivative.
*  Note: we need (1/r)(dV/dr) = 2 * dV/dr^2 --> use with equidistant r^2 
*  col is p_typ * ntypes + q_typ
*
******************************************************************************/

void pair_int2(real *pot, real *grad, int *is_short, pot_table_t *pt, 
               int col, int inc, real r2)
{
  real r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;
  int  k;

  /* check for distances shorter than minimal distance in table */
  r2a = MIN(r2,pt->end[col]);
  r2a = r2a - pt->begin[col];
  if (r2a < 0) {
    r2a   = 0;
    *is_short = 1;
  }

  /* indices into potential table */
  istep = pt->invstep[col];
  k     = (int) (r2a * istep);
  chi   = (r2a - k * pt->step[col]) * istep;

  /* intermediate values */
  ptr = PTR_2D(pt->table, k, col, pt->maxsteps, inc);
  p0  = *ptr; ptr += inc;
  p1  = *ptr; ptr += inc;
  p2  = *ptr;
  dv  = p1 - p0;
  d2v = p2 - 2 * p1 + p0;

  /* potential and twice the derivative */
  *pot  = p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;
  *grad = 2 * istep * (dv + (chi - 0.5) * d2v);
}

/*****************************************************************************
*
*  Evaluate potential table with cubic interpolation. 
*  Returns the potential value and twice the derivative.
*  Note: we need (1/r)(dV/dr) = 2 * dV/dr^2 --> use with equidistant r^2 
*  col is p_typ * ntypes + q_typ
*
******************************************************************************/

void pair_int3(real *pot, real *grad, int *is_short, pot_table_t *pt, 
               int col, int inc, real r2)
{
  real r2a, istep, chi, p0, p1, p2, p3, *ptr;
  real fac0, fac1, fac2, fac3, dfac0, dfac1, dfac2, dfac3;
  int  k;

  /* check for distances shorter than minimal distance in table */
  /* we need one extra value at the lower end for interpolation */
  r2a = MIN(r2,pt->end[col]);
  r2a = r2a - pt->begin[col] + pt->step[col];
  if (r2a < 0) {
    r2a = 0;
    *is_short = 1;
  }

  /* indices into potential table */
  istep = pt->invstep[col];
  k     = (int) (r2a * istep);
  chi   = (r2a - k * pt->step[col]) * istep;

  /* factors for the interpolation */
  fac0 = -(1.0/6.0) * chi * (chi-1.0) * (chi-2.0);
  fac1 =        0.5 * (chi*chi-1.0) * (chi-2.0);
  fac2 =       -0.5 * chi * (chi+1.0) * (chi-2.0);
  fac3 =  (1.0/6.0) * chi * (chi*chi-1.0);

  /* factors for the interpolation of the derivative */
  dfac0 = -(1.0/6.0) * ((3.0*chi-6.0)*chi+2.0);
  dfac1 =        0.5 * ((3.0*chi-4.0)*chi-1.0);
  dfac2 =       -0.5 * ((3.0*chi-2.0)*chi-2.0);
  dfac3 =    1.0/6.0 * (3.0*chi*chi-1.0);

  /* intermediate values */
  ptr = PTR_2D(pt->table, k, col, pt->maxsteps, inc);
  p0  = *ptr; ptr += inc;
  p1  = *ptr; ptr += inc;
  p2  = *ptr; ptr += inc;
  p3  = *ptr;

  /* potential energy */ 
  *pot = fac0 * p0 + fac1 * p1 + fac2 * p2 + fac3 * p3;

  /* twice the derivative */
  *grad = 2 * istep * (dfac0 * p0 + dfac1 * p1 + dfac2 * p2 + dfac3 * p3);
}

/*****************************************************************************
*
*  Evaluate tabulated function with quadratic interpolation. 
*
******************************************************************************/

void val_func2(real *val, int *is_short, pot_table_t *pt, 
               int col, int inc, real r2)
{
  real r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;
  int  k;

  /* check for distances shorter than minimal distance in table */
  r2a = MIN(r2,pt->end[col]);
  r2a = r2a - pt->begin[col];
  if (r2a < 0) {
    r2a = 0;
    *is_short = 1;
  }

  /* indices into potential table */
  istep = pt->invstep[col];
  k     = (int) (r2a * istep);
  chi   = (r2a - k * pt->step[col]) * istep;

  /* intermediate values */
  ptr = PTR_2D(pt->table, k, col, pt->maxsteps, inc);
  p0  = *ptr; ptr += inc;
  p1  = *ptr; ptr += inc;
  p2  = *ptr;
  dv  = p1 - p0;
  d2v = p2 - 2 * p1 + p0;

  /* the function value */
  *val = p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;
}

/*****************************************************************************
*
*  Evaluate tabulated function with cubic interpolation. 
*
******************************************************************************/

void val_func3(real *val, int *is_short, pot_table_t *pt, 
               int col, int inc, real r2)
{
  real r2a, istep, chi, p0, p1, p2, p3;
  real fac0, fac1, fac2, fac3, *ptr;
  int  k;

  /* check for distances shorter than minimal distance in table */
  /* we need one extra value at the lower end for interpolation */
  r2a = MIN(r2,pt->end[col]);
  r2a = r2a - pt->begin[col] + pt->step[col];
  if (r2a < 0) {
    r2a = 0;
    *is_short = 1;
  }

  /* indices into potential table */
  istep = pt->invstep[col];
  k     = (int) (r2a * istep);
  chi   = (r2a - k * pt->step[col]) * istep;

  /* factors for the interpolation */
  fac0 = -(1.0/6.0) * chi * (chi-1.0) * (chi-2.0);
  fac1 =        0.5 * (chi*chi-1.0) * (chi-2.0);
  fac2 =       -0.5 * chi * (chi+1.0) * (chi-2.0);
  fac3 =  (1.0/6.0) * chi * (chi*chi-1.0);

  /* intermediate values */
  ptr = PTR_2D(pt->table, k, col, pt->maxsteps, inc);
  p0  = *ptr; ptr += inc;
  p1  = *ptr; ptr += inc;
  p2  = *ptr; ptr += inc;
  p3  = *ptr;

  /* the function value */ 
  *val = fac0 * p0 + fac1 * p1 + fac2 * p2 + fac3 * p3;
}

/*****************************************************************************
*
*  Evaluate the derivative of a function with quadratic interpolation. 
*  Returns *twice* the derivative.
*
******************************************************************************/

void deriv_func2(real *grad, int *is_short, pot_table_t *pt, 
                 int col, int inc, real r2)
{
  real r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;
  int  k;

  /* check for distances shorter than minimal distance in table */
  r2a = MIN(r2,pt->end[col]);
  r2a = r2a - pt->begin[col];
  if (r2a < 0) {
    r2a   = 0;
    *is_short = 1;
  }

  /* indices into potential table */
  istep = pt->invstep[col];
  k     = (int) (r2a * istep);
  chi   = (r2a - k * pt->step[col]) * istep;

  /* intermediate values */
  ptr = PTR_2D(pt->table, k, col, pt->maxsteps, inc);
  p0  = *ptr; ptr += inc;
  p1  = *ptr; ptr += inc;
  p2  = *ptr;
  dv  = p1 - p0;
  d2v = p2 - 2 * p1 + p0;

  /* twice the derivative */
  *grad = 2 * istep * (dv + (chi - 0.5) * d2v);
}

/*****************************************************************************
*
*  Evaluate the derivative of a function with cubic interpolation. 
*  Returns *twice* the derivative.
*
******************************************************************************/

void deriv_func3(real *grad, int *is_short, pot_table_t *pt, 
                 int col, int inc, real r2)
{
  real r2a, istep, chi, p0, p1, p2, p3, *ptr;
  real dfac0, dfac1, dfac2, dfac3;
  int  k;

  /* check for distances shorter than minimal distance in table */
  /* we need one extra value at the lower end for interpolation */
  r2a = MIN(r2,pt->end[col]);
  r2a = r2a - pt->begin[col] + pt->step[col];
  if (r2a < 0) {
    r2a = 0;
    *is_short = 1;
  }

  /* indices into potential table */
  istep = pt->invstep[col];
  k     = (int) (r2a * istep);
  chi   = (r2a - k * pt->step[col]) * istep;

  /* factors for the interpolation of the 1. derivative */
  dfac0 = -(1.0/6.0) * ((3.0*chi-6.0)*chi+2.0);
  dfac1 =        0.5 * ((3.0*chi-4.0)*chi-1.0);
  dfac2 =       -0.5 * ((3.0*chi-2.0)*chi-2.0);
  dfac3 =    1.0/6.0 * (3.0*chi*chi-1.0);

  /* intermediate values */
  ptr = PTR_2D(pt->table, k, col, pt->maxsteps, inc);
  p0  = *ptr; ptr += inc;
  p1  = *ptr; ptr += inc;
  p2  = *ptr; ptr += inc;
  p3  = *ptr;

  /* twice the derivative */
  *grad = 2 * istep * (dfac0 * p0 + dfac1 * p1 + dfac2 * p2 + dfac3 * p3);
}
