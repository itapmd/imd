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
* averaging of temperature profiles
* 
******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

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
  printf("   Usage: %s <n1> <n2> <margin> <file.hcprof>\n\n", progname);
  printf("   Writes average of profiles n1 to n2 to stdout\n\n");
}

int main(int argc, char **argv)
{
  double *x, *y, *s, a, b, kappa, fact;
  double heatcurr, Sx, Sy, Sxy, Sx2;
  int    i, j, n, n1, n2, m, cnt;
  FILE   *out, *inf;
  str255 fname;
  char   *progname, *filename, c;

  progname = strdup(argv[0]);
  if (argc<5) {
    printf("\n   Not enough arguments!\n");
    usage(progname);
    exit(-1);
  }
  n1 = atoi(argv[1]);
  n2 = atoi(argv[2]);
  m  = atoi(argv[3]);
  filename = strdup(argv[4]);

  /* open input file */
  inf = fopen(filename, "r");
  if (NULL == inf) error_str("Cannot open input file %s",filename);

  /* read header */
  fscanf(inf, "%c %d %lf\n", &c, &n, &heatcurr);

  /* allocate arrays */
  x = (double *) malloc( n * sizeof(double) );
  y = (double *) malloc( n * sizeof(double) );
  s = (double *) malloc( n * sizeof(double) );
  if ((NULL==x) || (NULL==y) || (NULL==s)) error("Cannot allocate arrays.");
  for (j=0; j<n; j++) s[j] = 0.0;

  /* read away histograms 0 to n1-1 */
  for (i=0; i<n1; i++)
    for (j=0; j<n; j++)
      fscanf(inf, "%lf %lf\n", x+j, y+j);

  /* add up histograms n1 to n2 */
  for (i=n1; i<=n2; i++)
    for (j=0; j<n; j++) {
      fscanf(inf, "%lf %lf\n", x+j, y+j);
      s[j] += y[j];
    }
  for (j=0; j<n; j++) s[j] /= (n2-n1+1);
  fclose(inf);

  /* fit gradient of averaged profile */
  cnt = 0;
  for (j=m; j<n-m; j++) {
    Sx  += x[j];
    Sy  += s[j];
    Sxy += x[j] * s[j];
    Sx2 += x[j] * x[j];
    cnt++;
  }
  Sx  /= cnt;
  Sy  /= cnt;
  Sxy /= cnt;
  Sx2 /= cnt;
  a = (Sxy - Sx * Sy) / (Sx2 - Sx * Sx);
  b = Sy - a * Sx;

  kappa = heatcurr / a;
  /* conversion factor of kappa to SI units */
  fact  = 1.6022e-19 / (1.0179e-14 * 1e-10 * 11605);

  printf("# gradT deltaT kappa kappa[W/mK]\n");
  printf("# %10.4e %10.4e %10.4e %10.4e\n", a, b, kappa, fact*kappa);
  for (j=0; j<n; j++)
    printf("%10.4e %10.4e %10.4e\n", x[j], s[j], a*x[j]+b);

  return 0;
}
