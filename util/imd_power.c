
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
*  imd_power -- calculate power spectrum of IMD eng file
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

/*****************************************************************************
*
*  include files
*
*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "fft.c"

#define SQR(x) (x)*(x)

/*****************************************************************************
*
*  prototypes
*
*****************************************************************************/

int  main (int argc, char **argv);
void read_input(int argc, char **argv);
void read_engfile(void);
void usage(void);
void error(char *);
void write_output(void);

/*****************************************************************************
*
*  global variables
*
*****************************************************************************/

int n_min, n_max;      /* numbers of first and last line to be treated */   
int n_lines, n_cols;   /* number of lines and columns */
unsigned long n_alloc; /* length of FFT (power of 2) */
char progname[255];    /* name of executable (argv[0]) */
char fbasename[255];   /* base name of input/output files */
double *cols[16];      /* columns of engfile */

/*****************************************************************************
*
*  main
*
*****************************************************************************/

int main(int argc, char **argv)
{
  int i;
  read_input(argc,argv);
  for (i=1; i<n_cols; i++) 
    gsl_fft_real_radix2_transform(cols[i], 1, n_alloc);
  write_output();
  return 0;
}

/******************************************************************************
*
*  Usage -- educate users
*
******************************************************************************/

void usage(void)
{
  fprintf(stderr, "%s <n_min> <n_max> <engfilename>[.eng]\n", progname); 
  exit(1); 
}

/******************************************************************************
*
*  error -- complain and abort
*
******************************************************************************/

void error(char *mesg)
{
  fprintf(stderr,"Error: %s\n",mesg);
  exit(2);
}

/******************************************************************************
*
*  reads the command line parameters, read input 
*
******************************************************************************/

void read_input(int argc,char **argv)
{
  int p;

  /* get command line parameters */
  strcpy(progname,argv[0]);
  if (argc<4) {
    fprintf(stderr,"wrong number of arguments\n");
    usage();
  }
  n_min = atoi(argv[1]);
  n_max = atoi(argv[2]);
  n_lines = n_max - n_min + 1;
  n_alloc = 2;
  while (n_alloc < n_lines) n_alloc *= 2;

  /* get base file name, strip trailing .eng */
  strcpy(fbasename, argv[3]);
  p = strlen(fbasename) - 4;
  if (NULL!=strstr(fbasename + p, ".eng")) fbasename[p] = '\0';

  /* read .eng file */
  read_engfile();
} 

/******************************************************************************
*
*  read engfile
*
******************************************************************************/

void read_engfile()
{
  char engfilename[255], msg[255];
  FILE *infile;
  char buf[512];
  double r[16];
  double sum[16];
  int i,j,p;

  /* open .eng file */
  sprintf(engfilename, "%s.eng", fbasename);
  infile = fopen(engfilename,"r");
  if (NULL==infile) {
    sprintf(msg, "Cannot open file %s", engfilename);
    error(msg);
  }

  /* read the first line */
  buf[0] = '\0';
  fgets(buf,sizeof(buf),infile);

  /* eat comments */
  while ('#'==buf[0]) fgets(buf,sizeof(buf),infile);

  /* advance to line n_min */
  for (i=0; i<n_min; i++) fgets(buf,sizeof(buf),infile);

  /* read first line */
  n_cols = sscanf( buf,
    "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
    r,r+1,r+2,r+3,r+4,r+5,r+6,r+7,r+8,r+9,r+10,r+11,r+12,r+13,r+14,r+15);

  /* allocate data */
  for (j=0; j<n_cols; j++) {
    cols[j] = (double *) calloc(n_alloc, sizeof(double));
    cols[j][0] = r[j];
    sum[j] = r[j];
  }

  /* read further lines */
  for (i=1; i<n_lines; i++) {
    fgets(buf,sizeof(buf),infile); 
    p = sscanf( buf,
      "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
      r,r+1,r+2,r+3,r+4,r+5,r+6,r+7,r+8,r+9,r+10,r+11,r+12,r+13,r+14,r+15);
    if ((p<n_cols) || (feof(infile))) error("not enough data");
    for (j=0; j<n_cols; j++) {
      cols[j][i] = r[j];
      sum[j] += r[j];
    }
  }
  fclose(infile);  

  /* subtract average */
  for (j=1; j<n_cols; j++) {
    sum[j] /= n_lines;
    for (i=0; i<n_lines; i++) cols[j][i] -= sum[j];
  }
}

/******************************************************************************
*
*  write output
*
******************************************************************************/

void write_output() 
{
  char outfilename[255], msg[255];
  FILE *outfile;
  int i,j;
  double omega;

  /* get frequency scale */
  n_alloc /= 2;
  omega = 4 * atan(1.0) / (n_alloc*(cols[0][1]-cols[0][0]));

  /* open output file */
  sprintf(outfilename, "%s.fft", fbasename);
  outfile = fopen(outfilename,"w");
  if (NULL==outfile) {
    sprintf(msg, "Cannot open output file %s", outfilename);
    error(msg);
  }

  /* write lowest frequency */
  fprintf(outfile, "%f ", 0.0);
  for (j=1; j<n_cols; j++) fprintf(outfile, "%f ", SQR(cols[j][0]));
  fprintf(outfile, "\n");

  /* write intermediate frequencies */
  for (i=1; i<n_alloc; i++) {
    fprintf(outfile, "%f ", i * omega);
    for (j=1; j<n_cols; j++) 
      fprintf(outfile, "%f ", SQR(cols[j][i])+SQR(cols[j][2*n_alloc-i]));
    fprintf(outfile, "\n");
  }

  /* write highest frequency */
  fprintf(outfile, "%f ", n_alloc * omega);
  for (j=1; j<n_cols; j++) fprintf(outfile, "%f ", SQR(cols[j][n_alloc]));
  fprintf(outfile, "\n");

  /* close output file */
  fclose(outfile);
}
