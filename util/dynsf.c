
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
* Postprocessing of dynamical structure factor files
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

double  t_unit = 0.01018;  /* IMD time unit in ps */
typedef char str255[255];
typedef struct {double x, y, z; } vektor;

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
  printf("   Usage: %s [options] <file.dsf>\n\n", progname);
  printf("   Compute dynamical structure factor.\n\n");
  printf("   Options:  -w <width> Frequency smoothing width (FWHM)\n");
  printf("             -i <int>   Frequency writing interval\n");
  printf("             -c <cut>   Frequency cutoff\n");
  printf("             -e         Frequency and width in meV (default THz)\n");
  printf("             -s         Separate reords by two emtpy lines\n");
  printf("             -h         This help\n\n");
}

int main(int argc, char **argv)
{
  double *data, *sum, *klen, *cdat, *rdat, *dat;
  double width=-1.0, smooth, ts, f_unit=1.0, freq, omega_cut;
  int    *k0, *kdir, *kmax, *koff;
  int    i, j, jmax, k, kk, nk, nks, dim, n1, n2;
  int    count, cont, input_endian, omega_int = 8;
  FILE   *out, *inf;
  vektor tbox_x, tbox_y, tbox_z;
  str255 fname;
  char   *progname, *filename, *separator="\n", *tmpc;
  fftw_plan p1, p2, p3;

  progname = strdup(argv[0]);
  while ((argc > 1) && (argv[1][0] =='-')) {
    if (strcasecmp(argv[1],"-w" )==0) { 
      if (argc<3) {
        printf("\n   Not enough arguments!\n");
        usage(progname);
        exit(-1);
      }
      width = atof(argv[2]);
      argc -= 2;
      argv += 2;
    }
    else if (strcasecmp(argv[1],"-s" )==0) { 
      if (argc<2) {
        printf("\n   Not enough arguments!\n");
        usage(progname);
        exit(-1);
      }
      separator = "\n\n";
      argc -= 1;
      argv += 1;
    }
    else if (strcasecmp(argv[1],"-e" )==0) { 
      if (argc<2) {
        printf("\n   Not enough arguments!\n");
        usage(progname);
        exit(-1);
      }
      f_unit = 4.1375;
      argc -= 1;
      argv += 1;
    }
    else if (strcasecmp(argv[1],"-i" )==0) { 
      if (argc<3) {
        printf("\n   Not enough arguments!\n");
        usage(progname);
        exit(-1);
      }
      omega_int = atoi(argv[2]);
      argc -= 2;
      argv += 2;
    }
    else if (strcasecmp(argv[1],"-c" )==0) { 
      if (argc<3) {
        printf("\n   Not enough arguments!\n");
        usage(progname);
        exit(-1);
      }
      omega_cut = atof(argv[2]);
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
  if (width<0.0) error("smoothing width must be specified");

  filename = strdup(argv[1]);
  tmpc = strstr( filename, ".dsf" );
  if (NULL!=tmpc) tmpc[0] = '\0';

  /* open input file */
  sprintf(fname, "%s.dsf", filename);
  inf = fopen(fname, "r");
  if (NULL == inf) error_str("Cannot open input file %s",fname);

  /* read header */
  cont = 1; nk = 0;
  do {
    str255 str, line;
    fgets(line, 255, inf);
    if (line[0]!='#') {
      error("file header corrupt!");
    } else {
      if (line[1]=='F') {
        sscanf(line+2, "%s %d %d", str, &dim, &nks );
        k0   = (int *) malloc( dim * nks * sizeof(int) );
        kdir = (int *) malloc( dim * nks * sizeof(int) );
        kmax = (int *) malloc(       nks * sizeof(int) );
        koff = (int *) malloc(       nks * sizeof(int) );
        if ((NULL==k0) || (NULL==kdir) || (NULL==kmax) || (NULL==koff))
          error("cannot allocate data");
        if      (str[0]=='B') input_endian = 1;
        else if (str[0]=='L') input_endian = 0;
        else error("file header corrupt (F)!"); 
      } else  if (line[1]=='T') {
        sscanf(line+2, "%lf", &ts);
      } else  if (line[1]=='X') {
        if (dim!=sscanf(line+2,"%lf %lf %lf",&tbox_x.x,&tbox_x.y,&tbox_x.z))
          error("file header corrupt (X)!");
      } else  if (line[1]=='Y') {
        if (dim!=sscanf(line+2,"%lf %lf %lf",&tbox_y.x,&tbox_y.y,&tbox_y.z))
          error("file header corrupt (Y)!");
      } else  if (line[1]=='Z') {
        if (dim!=sscanf(line+2,"%lf %lf %lf",&tbox_z.x,&tbox_z.y,&tbox_z.z))
          error("file header corrupt (Z)!");
      } else  if (line[1]=='K') {
        if (nk>=nks) error("too many k-points!");
        if (dim==2) {
          if (5 != sscanf(line+2, "%d %d %d %d %d",
            k0+2*nk, k0+2*nk+1, kdir+2*nk, kdir+2*nk+1, kmax+nk))
            error("file header corrupt (K)!");
          nk++;
        } else if (dim==3) {
          if (7 != sscanf(line+2, "%d %d %d %d %d %d %d",
            k0  +3*nk, k0  +3*nk+1, k0  +3*nk+2,
            kdir+3*nk, kdir+3*nk+1, kdir+3*nk+2, kmax+nk))
            error("file header corrupt (K)!");
          nk++;
        }
      } else  if (line[1]=='E') {
	cont=0;
      }
    }
  } while (cont);

  /* initialize k-vectors */
  if (nk<nks) error("too few k-points!");
  klen = (double *) malloc( nks * sizeof(double) );
  if (NULL==klen) error("Cannot allocate data!");
  nk = 0;
  for (i=0; i<nks; i++) {
    koff[i] = nk;
    nk += kmax[i] + 1;
    if (dim==2) {
      klen[i]  = SQR(kdir[2*i] * tbox_x.x + kdir[2*i+1] * tbox_y.x);
      klen[i] += SQR(kdir[2*i] * tbox_x.y + kdir[2*i+1] * tbox_y.y);
      klen[i]  = sqrt(klen[i]);
    } else if (dim==3) {
      klen[i]  = SQR(kdir[3*i] * tbox_x.x + kdir[3*i+1] * tbox_y.x 
                                          + kdir[3*i+2] * tbox_z.x);
      klen[i] += SQR(kdir[3*i] * tbox_x.y + kdir[3*i+1] * tbox_y.y 
                                          + kdir[3*i+2] * tbox_z.y);
      klen[i] += SQR(kdir[3*i] * tbox_x.z + kdir[3*i+1] * tbox_y.z 
                                          + kdir[3*i+2] * tbox_z.z);
      klen[i]  = sqrt(klen[i]);
    }
  }
  dat = (double *) malloc( 2 * nk * sizeof(double) );
  if (NULL==dat) error("cannot allocate data");

  /* get number of samples */
  count = 0;
  while (!feof(inf)) {
    if (2 * nk == fread(dat, sizeof(double), 2 * nk, inf)) count++;
  }
  fclose(inf);
  printf("%d samples\n", count);
  n1 = count;
  n2 = 2*n1;
  freq   = f_unit / (n2 * ts * t_unit);
  width *= (1.34 / freq);
  jmax   = omega_cut > 0 ? omega_cut / freq : n1 / 2;

  /* get input */
  data = (double *) malloc( 2 * nk * n1 * sizeof(double) );
  sum  = (double *) malloc( 2 * nk      * sizeof(double) );
  if ((NULL==data) || (sum==NULL)) error("cannot allocate data");
  for (j=0; j<2*nk; j++) sum[j] = 0.0;

  inf = fopen(fname, "r");
  if (NULL == inf) error_str("Cannot open input file %s",fname);

  /* discard header */
  cont=1;
  do {
    str255 line;
    fgets(line, 255, inf);
    if ((line[0]=='#') && (line[1]=='E')) cont=0;
  } while (cont);

  /* read data */
  count = 0;
  while (!feof(inf)) {
    if (2*nk == fread(dat, sizeof(double), 2*nk, inf)) {
      for (j=0; j<nk; j++) {
        data[2*j*n1   + 2*count] = dat[2*j  ]; 
        data[2*j*n1+1 + 2*count] = dat[2*j+1]; 
        sum [2*j  ] += dat[2*j  ];
        sum [2*j+1] += dat[2*j+1];
      }
      count++;
    }
  }
  fclose(inf);
  free(dat);
  for (j=0; j<2*nk; j++) sum[j] /= n1;

  /* allocate arrays */
  cdat = (double *) fftw_malloc( 2 * n2   * sizeof(double) );
  rdat = (double *) fftw_malloc( (n2 + 2) * sizeof(double) );
  if ((NULL== cdat) || (NULL==rdat)) error("cannot allocate arrays");

  /* create fftw plans */
  p1 = fftw_plan_dft_1d( n2, (fftw_complex *) cdat, (fftw_complex *) cdat, 
                         FFTW_BACKWARD, FFTW_ESTIMATE );
  p2 = fftw_plan_dft_r2c_1d(n2, rdat, (fftw_complex *) rdat, FFTW_ESTIMATE);
  p3 = fftw_plan_dft_c2r_1d(n2, (fftw_complex *) rdat, rdat, FFTW_ESTIMATE);

  /* transform and write data */
  for (i=0; i<nks; i++) {
    sprintf(fname, "%s.%d.plot", filename, i);
    out = fopen(fname, "w");
    if (NULL == out) error_str("Cannot open output file %s", fname);
    for (k=0; k<=kmax[i]; k++) {
      kk = 2 * (k + koff[i]);
      for (j=0; j<n1; j++) {
        cdat[2*j  ] = data[kk*n1 + 2*j  ] - sum[kk  ];
        cdat[2*j+1] = data[kk*n1 + 2*j+1] - sum[kk+1];
      }
      for (j=0; j<n2; j++) cdat[j+n2] = 0.0;
      fftw_execute(p1);
      for (j=0; j<n2; j++) rdat[j] = (SQR(cdat[2*j]) + SQR(cdat[2*j+1])) / n2;
      fftw_execute(p2);
      for (j=0; j<n1+1; j++) {
        smooth = exp( -0.5 * SQR( width * j / n1 ) ) / (n1 + 1 - j);
        rdat[2*j  ] *= smooth;
        rdat[2*j+1] *= smooth;
      }
      fftw_execute(p3);
      for (j=0; j<jmax; j+=omega_int) 
        fprintf(out, "%e %e %e\n", k * klen[i], j * freq, 
                     MAX(rdat[j] * 1e5, 1e-8) );
      fprintf(out, separator);
    }
    fclose(out);
  }
}
