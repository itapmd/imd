
/******************************************************************************
*
*  IMD -- The ITAP Molecular Dynamics Program
*
*  Copyright 1996-2001 Institute for Theoretical and Applied Physics,
*  University of Stuttgart, D-70550 Stuttgart
*
*  $Revision$
*  $Date$
*
******************************************************************************/

/******************************************************************************
*
*  This utility program can be used to shorten pair potentials. For each 
*  data column, a cutoff radius (or radius squared) must be specified. 
*  Up to this radius, the returned potential is equal to the input 
*  potential, but beyond this radius it is replaced by a function which 
*  tends smoothly to zero. The resulting potential is continuous and has 
*  a continuous first derivative.  
*
*  Each line of the input potential starts either with a radius or a squared
*  radius, which is followed by the corresponding potential values. The
*  radii or squared radii, whichever are used, must be equally spaced.
*  For the output potential, only equally spaced, squared radii are used.
*  This format is required by IMD. The cutoff values must be given in the
*  same scale as the input potential, i.e., as radii for linear spacing,
*  and as squared radii for quadratic spacing.
*
*  There are two tail functions possible, a quadratic polynomial b(a-x)^2,
*  and an exponential function b*exp(-c/(a-x)). In the latter case, for
*  each column a further parameter, the steepness, can be specified, which
*  controls how fast the potential goes to zero. The default steepness
*  value is 2.0.
*
*  cutpot requires the following parameters (here with example values), 
*  which must be given in the parameter file:
*
*  ncols      4                  # number of data (potential) columns
*  nsteps     2000               # number of lines in output potential
*  infile     input.pot          # input potential file
*  outfile    output.pot         # output potential file
*  spacing    lin                # spacing of input potential; lin or sqr
*  tailtype   sqr                # type of tail function; sqr or exp
*  rcut       1.4 1.6 1.6 1.6    # cutoff values, where tail function starts
*  steepness  1.5 1.5 1.5 1.5    # steepness of exponential tail function
*
*  The optional keyword 'symmetric' in the parameter file indicates that
*  the potential is symmetric in the atom types, and the input file 
*  contains only the ntypes*(ntypes+1)/2 potential values V_ij with i<=j.
*  cutpot then writes all ntypes*ntypes values to the output file, as
*  required for IMD.
*
*  Compilation:  gcc -o cutpot -O cutpot.c -lm
*
*  Usage:        cutpot <paramfile>
*
******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "param.h"

#define SQRTAIL 0
#define EXPTAIL 1
#define LINSPACING 0
#define SQRSPACING 1

#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define PTR_2D(var,i,j,dim_i,dim_j) (((var) + ((i)*(dim_j)) + (j)))

/* memory allocation increment for potential */
#define PSTEP 50

typedef char str255[255];

typedef double real;

/* data structure to store a potential table or a function table */ 
typedef struct {
  real *begin;      /* first value in the table */
  real *end;        /* last value in the table (followed by extra zeros) */
  real *step;       /* table increment */
  real *invstep;    /* inverse of increment */
  int  maxsteps;    /* physical length of the table */
  real *table;      /* the actual data */
} pot_table_t;

/* global variables */
str255 infilename="\0", outfilename="\0";
int    ncols=0, ntypes=0, sym=0, nsteps=0, spacing=-1, tailtype=-1;
real   *rcut=NULL, *a=NULL, *b=NULL, *c=NULL, *steepness=NULL;

/******************************************************************************
*
* error -- Complain and abort
*
******************************************************************************/

void error(char *msg)
{
  fprintf(stderr,"Error: %s\n",msg);
  fflush(stderr);
  exit(2);
}


/*****************************************************************************
*
*  read tag-based parameter file
*  lines beginning with comment characters '#' or blank lines are skipped
*
*****************************************************************************/

void getparamfile(char *paramfname)
{
  FILE *pf;
  char buffer[1024];
  char *token;
  char *res;
  str255 tmp;

  pf = fopen(paramfname,"r");
  if (NULL == pf) {
    sprintf(tmp,"Cannot open parameter file %s",paramfname);
    error(tmp);
  }

  do {
    res=fgets(buffer,1024,pf);
    if (NULL == res) break; /* probably EOF reached */
    token = strtok(buffer," =\t\r\n");
    if (NULL == token) continue; /* skip blank lines */
    if (token[0]=='#') continue; /* skip comments */

    if (strcasecmp(token,"infile")==0) {
      /* file name for input potential */
      getparam(token,infilename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"outfile")==0) {
      /* file name for output potential */
      getparam(token,outfilename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"ncols")==0) {
      /* number of columns */
      getparam(token,&ncols,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"symmetric")==0) {
      /* symmetric potential, ntypes*(ntypes+1)/2 columns */
      sym=1;
    }
    else if (strcasecmp(token,"nsteps")==0) {
      /* number of output potential lines */
      getparam(token,&nsteps,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"spacing")==0) {
      /* linear or square spacing of input table */
      getparam(token,tmp,PARAM_STR,1,255);
      if      (strcmp((const char *)tmp,"lin")==0) spacing = LINSPACING;
      else if (strcmp((const char *)tmp,"sqr")==0) spacing = SQRSPACING;
      else error("parameter 'spacing' must be 'lin' or 'sqr'");
    }
    else if (strcasecmp(token,"tailtype")==0) {
      /* quadratic or exponetial potential tail */
      getparam(token,tmp,PARAM_STR,1,255);
      if      (strcmp((const char *)tmp,"sqr")==0) tailtype = SQRTAIL;
      else if (strcmp((const char *)tmp,"exp")==0) tailtype = EXPTAIL;
      else error("parameter 'tailtype' must be 'sqr' or 'exp'");
    }
    else if (strcasecmp(token,"rcut")==0) {
      /* cutoff radii */
      if (ncols==0) error("specify ncols before rcut");
      rcut = (real *) malloc( ncols * sizeof(real) );
      if (rcut==NULL) error("allocation of rcut failed");
      getparam(token,rcut,PARAM_REAL,ncols,ncols);
    }
    else if (strcasecmp(token,"steepness")==0) {
      /* c */
      if (ncols==0) error("specify ncols before steepness");
      steepness = (real *) malloc( ncols * sizeof(real) );
      if (rcut==NULL) error("allocation of steepness failed");
      getparam(token,steepness,PARAM_REAL,ncols,ncols);
    }
  } while (!feof(pf));
  fclose(pf);

  /* check whether we have everything */
  printf("\nParameters:\n\n");
  if (ncols==0) 
    error("ncols must be specified in parameter file");
  else
    printf("Using %d data column(s)\n",ncols);

  if (sym==1) {
    ntypes = (int) sqrt(2.0*ncols);
    if (ncols!=ntypes*(ntypes+1)/2)
      error("keyword \"symmetric\" implies ntypes*(ntypes+1)/2 input columns");
  }

  if (nsteps==0) 
    error("nsteps must be specified in parameter file");
  else
    printf("Writing %d potential lines\n",nsteps);

  if (infilename[0]=='\0') 
    error("infilename must be specified in parameter file");
  else
    printf("Reading input from %s\n",infilename);

  if (outfilename[0]=='\0') 
    error("outfilename must be specified in parameter file");
  else
    printf("Writing output to %s\n",outfilename);

  if (spacing==-1) 
    error("input spacing type must be specified in parameter file");
  else
    if (spacing==LINSPACING)
      printf("Assuming linear input spacing\n");
    else
      printf("Assuming square input spacing\n");

  if (tailtype==-1) 
    error("type of potential tail must be specified in parameter file");
  else {
    if (tailtype==SQRTAIL)
      printf("Writing quadratic potential tail\n");
    else
      printf("Writing exponential potential tail\n");
  }

  if (rcut==NULL)
    error("rcut must be specified in parameter file");

  if (steepness==NULL) {
    int i;
    steepness = (real *) malloc( ncols * sizeof(real) );
    if (steepness==NULL) error("data allocation failed");
    for (i=0; i<ncols; i++) steepness[i] = 2.0;
  }

} /* getparamfile */


/******************************************************************************
*
*   read the command line parameters, and the parameter file given 
*   on the command line
*
******************************************************************************/

void read_parameters(int argc,char **argv)
{
  str255 progname, paramfilename, msg;

  strcpy(progname,argv[0]);
  if (argc<2) {
    sprintf(msg,"Parameter file missing!\nUsage: %s <paramfile>",progname);
    error(msg);
  }
  strcpy(paramfilename,argv[1]);
  getparamfile(paramfilename);
} 


/*****************************************************************************
*
*  read potential table (routine from imd_param.c)
*
*****************************************************************************/

void read_pot_table1( pot_table_t *pt, char *filename )
{
  FILE *infile;
  int i, k;
  int size, tablesize, npot=0;
  real val, numstep, delta;

  real r2, r2_start, r2_step;
  str255 msg;

  /* allocate info block of function table */
  pt->maxsteps = 0;
  pt->begin    = (real *) malloc(ncols*sizeof(real));
  pt->end      = (real *) malloc(ncols*sizeof(real));
  pt->step     = (real *) malloc(ncols*sizeof(real));
  pt->invstep  = (real *) malloc(ncols*sizeof(real));
  if ((pt->begin   == NULL) || (pt->end == NULL) || (pt->step == NULL) || 
      (pt->invstep == NULL)) {
    sprintf(msg,"Cannot allocate info block for function table %s.",filename);
    error(msg);
  }

  /* catch the case where potential is identically zero */
  for (i=0; i<ncols; ++i) {
    pt->end[i] = 0.0;

  }

    infile = fopen(filename,"r");
    if (NULL==infile) {
      sprintf(msg,"Cannot open file %s.",filename);
      error(msg);
    }

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
	if ( 1 != fscanf(infile,"%lf", &val)) 
           error("Line incomplete in potential file.");
	*PTR_2D(pt->table,npot,i,pt->maxsteps,ncols) = val;
        if (val!=0.0) pt->end[i] = r2; /* catch last non-zero value */
      }
      ++npot;
    }

    fclose(infile);

    r2_step = (r2 - r2_start) / (npot-1);

    printf("\nRead input potential %s with %d lines\n",filename,npot);
    printf("Input potential starts at %f, ends at %f\n",r2_start,r2);

    /* fill info block, and shift potential to zero */
    for (i=0; i<ncols; ++i) {
      pt->begin[i] = r2_start;
      pt->step[i] = r2_step;
      pt->invstep[i] = 1.0 / r2_step;
      pt->end[i] += r2_step;
    }

    /* The interpolation uses k+1 and k+2, so we add zeros at end of table */
    for (k=1; k<=3; ++k) {
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
*  Evaluate potential table with quadratic interpolation. 
*  Returns the potential value and the derivative.
*  col is p_typ * ntypes + q_typ
*
******************************************************************************/

void pair_int2(real *pot, real *grad, pot_table_t *pt, 
               int col, int inc, real r2)
{
  real r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;
  int  k;

  /* check for distances shorter than minimal distance in table */
  r2a = MIN(r2,pt->end[col]);
  r2a = r2a - pt->begin[col];
  if (r2a < 0) {
    r2a   = 0;
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
  *grad = istep * (dv + (chi - 0.5) * d2v);
}


/*****************************************************************************
*
*  main - cut potential 
*
******************************************************************************/

int main(int argc, char **argv)
{
  int  i, j, k, col, *off;
  real x, r2, rmax=0, r2step;
  pot_table_t pt; 
  real grad, *pot;    
  FILE *out;

  read_parameters(argc, argv);

  a   = (real *) malloc( ncols * sizeof(real) );
  b   = (real *) malloc( ncols * sizeof(real) );
  c   = (real *) malloc( ncols * sizeof(real) );
  pot = (real *) malloc( ncols * sizeof(real) );
  off = (int  *) malloc( ncols * sizeof(int ) );
  if ((a==NULL) || (b==NULL) || (c==NULL) || (pot==NULL) || (off==NULL)) 
    error("data allocation failed");

  /* read input potential */
  read_pot_table1(&pt,infilename);
  
  /* prepare data for tail computation */
  printf("\nCutoff Data:\n\n");
  printf("col  potential   gradient       rcut       rmax\n");
  for (col=0; col<ncols; col++) {

    pair_int2(pot+col,&grad,&pt,col,ncols,rcut[col]);

    if (tailtype==SQRTAIL) {
      a[col] = rcut[col]-2*pot[col]/grad; 
      b[col] = grad*grad*0.25/pot[col];
    } else {
      a[col] = rcut[col]-steepness[col]*pot[col]/grad;
      b[col] = pot[col] * exp(steepness[col]);
      c[col] = (a[col]-rcut[col]) * steepness[col];
    }

    if (spacing==SQRSPACING) {
      if (rmax<sqrt(a[col])) rmax=sqrt(a[col]);
    } else {
      if (rmax<a[col]) rmax=a[col];
    }
    printf("%3d %10.6f %10.6f %10.6f %10.6f\n",
           col,pot[col],grad,rcut[col],a[col]);
    if (pot[col]*grad>=0) {
      fflush(stdout);
      error("potential and gradient must have different sign!");
    }
  }
  printf("\n");

  /* compute step width */
  if (spacing==SQRSPACING) {
    r2step = (rmax*rmax - pt.begin[0])/(nsteps-10);
  } else {
    r2step = (rmax*rmax - pt.begin[0]*pt.begin[0])/(nsteps-10);
  }

  /* compute offsets */
  if (sym==1) {
    off[0]=0;
    for (i=1;i<ntypes;i++) off[i]=off[i-1]+i;
  }

  /* write output potential */
  out=fopen(outfilename,"w");
  if (sym==1) fprintf(out,"#F 1 %d\n#E\n",ntypes*ntypes);
  else        fprintf(out,"#F 1 %d\n#E\n",ncols);

  for (i=0; i<nsteps; i++) {

    /* radius */
    if (spacing==SQRSPACING) {
      r2 = i*r2step+pt.begin[0];
      x  = r2;
    } else {
      r2 = i*r2step+pt.begin[0]*pt.begin[0];
      x  = sqrt(r2);
    }
    fprintf(out,"%e", r2);

    /* compute potential values */
    for (col=0; col<ncols; col++) {
      if (x<rcut[col])
        pair_int2(pot+col,&grad,&pt,col,ncols,x);
      else if (x<a[col])
        if (tailtype==SQRTAIL) {
          pot[col] = b[col]*(x-a[col])*(x-a[col]);
	} else {
          pot[col] = b[col]*exp(-(c[col]/(a[col]-x)));
	}
      else pot[col] = 0;
    }

    /* write potential values */
    if (sym==0) {
      for (col=0; col<ncols; col++) fprintf(out," %e", pot[col]);
    } else {
      for (j=0;j<ntypes;j++)
        for (k=0;k<ntypes;k++)
	  if (j<=k) fprintf(out," %e", pot[j*ntypes-off[j]+k]);
          else      fprintf(out," %e", pot[k*ntypes-off[k]+j]);
    }
    fprintf(out,"\n");
  }

  fclose(out);
  return 1;
}

