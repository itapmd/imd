
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2007 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/*****************************************************************************
* $Revision$
* $Date$
******************************************************************************/

/*****************************************************************************
*
* Autocorrelation of heat current
*
******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <fftw3.h>

#define SQR(x) ((x)*(x))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

typedef char str255[255];

/* error message  */
void error(char *msg)
{
  fprintf(stderr,"Error: %s\n",msg);
  exit(2);
}

/* error message built from two strings */
void error_str(char *msg, char *str)
{
  char buf[255];
  sprintf(buf, msg, str);
  error(buf);
}

void usage(char *progname)
{
  printf("\n");
  printf("   Usage: %s [options] <file.hc>\n\n", progname);
  printf("   Computes heat current autocorrelation.\n\n");
  printf("             -c <t_cut> Time cutoff [ps]\n");
  printf("             -s <start> starting point\n");
  printf("             -e <end>   end point\n");
  printf("             -h         This help\n\n");
}

/* compute autocorrelation of real time series */
void auto_corr(int n, double *x, double *y)
{
  fftw_plan p1, p2;
  double *r, *c;
  int i, m;

  /* allocate data */
  m=1; while (m<n) m*=2;
  r = (double *) fftw_malloc( 2 *  m    * sizeof(double) );
  c = (double *) fftw_malloc( 2 * (m+1) * sizeof(double) );
  if ((r==NULL) || (c==NULL)) error("cannot allocate data");

  /* make plans */
  p1 = fftw_plan_dft_r2c_1d( 2*m, r, (fftw_complex *) c, FFTW_ESTIMATE );
  p2 = fftw_plan_dft_c2r_1d( 2*m, (fftw_complex *) c, r, FFTW_ESTIMATE );

  /* forward transform */
  for (i=0; i<n;   i++) r[i] = x[i];
  for (i=n; i<2*m; i++) r[i] = 0.0;
  fftw_execute(p1);

  /* backward transform */
  for (i=0; i<2*(m+1); i+=2) {
    c[i  ] = SQR(c[i]) + SQR(c[i+1]);
    c[i+1] = 0.0;
  }
  c[0] = 0.0; /* eliminate dc part */
  fftw_execute(p2);
  for (i=0; i<n; i++) y[i] = r[i] / (n - i) / (2.0*m);

  /* cleanup */
  fftw_free(r);
  fftw_free(c);
  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);

}

int main(int argc, char **argv)
{
  double *t, *x, *y, *z, v[4], t_cut = 10.0, t_unit = 0.01018;
  int    i, n, dim, k, max, s=0, e=0;
  FILE   *out, *inf;
  str255 fname, line;
  char   *progname, *filename, *tmpc;

  progname = strdup(argv[0]);
  while ((argc > 1) && (argv[1][0] =='-')) {
    if (strcasecmp(argv[1],"-c" )==0) {
      if (argc<3) {
        printf("\n   Not enough arguments!\n");
        usage(progname);
        exit(-1);
      }
      t_cut = atof(argv[2]);
      argc -= 2;
      argv += 2;
    }
    else if (strcasecmp(argv[1],"-s" )==0) {
      if (argc<3) {
        printf("\n   Not enough arguments!\n");
        usage(progname);
        exit(-1);
      }
      s     = atoi(argv[2]);
      argc -= 2;
      argv += 2;
    }
    else if (strcasecmp(argv[1],"-e" )==0) {
      if (argc<3) {
        printf("\n   Not enough arguments!\n");
        usage(progname);
        exit(-1);
      }
      e     = atoi(argv[2]);
      argc -= 2;
      argv += 2;
    }
    else if (strcasecmp(argv[1],"-h" )==0) {
      usage(progname);
      exit(-1);
    }
    else {
      printf("\n   Unknown option!\n");
      usage(progname);
      exit(-1);
    }
  }
  if (argc<2) {
    printf("\n   Not enough arguments!\n");
    usage(progname);
    exit(-1);
  }

  /* open input file */
  filename = strdup(argv[1]);
  inf = fopen(filename, "r");
  if (NULL == inf) error_str("Cannot open input file %s",filename);

  /* count number of lines */
  fgets(line, 255, inf);
  while ('#'==line[0]) fgets(line, 255, inf);
  dim = sscanf( line, "%lf %lf %lf %lf", v, v+1, v+2, v+3);
  dim--;
  n = 1;
  while (!feof(inf)) {
    fgets(line, 255, inf);
    n++;
  }

  /* allocate data */
  t = (double *) malloc( 2 * n * sizeof(double) );
  x = (double *) malloc( 2 * n * sizeof(double) );
  y = (double *) malloc( 2 * n * sizeof(double) );
  z = (double *) malloc( 2 * n * sizeof(double) );
  if ((t==NULL) || (x==NULL) || (y==NULL) || (z==NULL))
    error("cannot allocate data");

  /* get input */
  rewind(inf);
  fgets(line, 255, inf);
  while ('#'==line[0]) fgets(line, 255, inf);
  i = 0;
  while (!feof(inf)) {
    k = sscanf( line, "%lf %lf %lf %lf", t+i, x+i, y+i, z+i);
    if (k==dim+1) i++;
    fgets(line, 255, inf);
  }
  fclose(inf);
  if (i>n) error("array overflow!");
  n = i;
  printf("dim = %d, %d samples\n", dim, n);
  if ((e==0) || (e>n)) e = n;

  /* transform data */
  auto_corr(e-s, x+s, x+s);
  auto_corr(e-s, y+s, y+s);
  if (dim==3)
    auto_corr(e-s, z+s, z+s);

  /* write output */
  sprintf(filename, "%s.ac", filename);
  out = fopen(filename, "w");
  if (NULL == out) error_str("Cannot open output file %s", filename);
  t_cut = t_cut / t_unit;
  i = 0;
  if (dim==2) {
    while (t[i] < t_cut) {
      fprintf(out, "%e %e %e %e %e %e %e\n", t[i]*t_unit, x[i+s], y[i+s], (x[i+s]+y[i+s])/2., x[i+s]/x[s], y[i+s]/y[s], (x[i+s]/x[s]+y[i+s]/y[s])/2.);
      i++;
    }
  }
  else {
    while (t[i] < t_cut) {
      fprintf(out, "%e %e %e %e %e %e %e %e %e\n", t[i]*t_unit, x[i+s], y[i+s], z[i+s], (x[i+s]+y[i+s]+z[i+s])/3., x[i+s]/x[s], y[i+s]/y[s], z[i+s]/z[s], (x[i+s]/x[s]+y[i+s]/y[s]+z[i+s]/z[s])/3.);
      i++;
    }
  }
  fclose(out);

  return 1;

}
