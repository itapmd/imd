
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
*  The utility program sqd2hist converts .sqd files (containing the square 
*  displacements for each atom) to displacement histograms. On each invocation,
*  histograms for one given atom type are produced.
*
*  Compilation:  gcc -o sqd2hist -O sqd2hist.c -lm
*
*  Usage:        sqd2hist <sqd_file> <atom_type> [<num_bins>]
*
*  The default value for <num_bins> is 200 histogram bins. 
*  The format of each output line is:
*
*  displacement num_tot num_x num_y [num_z]
*
*  num_tot is the number of atoms with this displacement, num_x the number
*  of atoms with this displacement in x-direction, etc.
* 
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))

void error(char *msg)
{
  fprintf(stderr,"Error: %s\n",msg);
  exit(2);
}

int main(int argc, char **argv) 
{
  FILE   *infile, *outfile;
  char   line[255];
  double dr, x2, y2, z2=0.0, max, xmax=0.0, ymax=0.0, zmax=0.0;
  int    size=200, ty, t, i, *hist, *x_hist, *y_hist, *z_hist;

  /* scan and check arguments */
  if (argc<3) {
    printf("Usage:  %s <sqd_file> <atom_type> [<n_bins>]\n",argv[0]);
    exit(1);
  }
  sscanf(argv[2], "%d", &ty);
  if (argc==4) sscanf(argv[3], "%d", &size);

  /* get histogram size */
  if (NULL==(infile=fopen(argv[1],"r"))) error("Cannot open input file");
  while (!feof(infile)) {
    fgets(line,255,infile);
    if (line[0]!='#') {
      sscanf(line,"%d %lf %lf %lf",&t,&x2,&y2,&z2);
      if (t==ty) {
        xmax = MAX(xmax,x2);
	ymax = MAX(ymax,y2);
	zmax = MAX(zmax,z2);
      }
    }
  }
  fclose(infile);

  /* allocate histograms */
  max    = sqrt(xmax+ymax+zmax);
  dr     = max/(size-2);
    hist = (int *) calloc(size,sizeof(int));
  x_hist = (int *) calloc(size,sizeof(int));
  y_hist = (int *) calloc(size,sizeof(int));
  z_hist = (int *) calloc(size,sizeof(int));
  if ((NULL==hist) || (NULL==x_hist) || (NULL==y_hist) || (NULL==z_hist))
    error("Cannot allocate histograms");

  /* process sqd-file */
  if (NULL==(infile=fopen(argv[1],"r"))) error("Cannot open input file");
  while (!feof(infile)) {
    fgets(line,255,infile);
    if (line[0]!='#') {
      sscanf(line,"%d %lf %lf %lf",&t,&x2,&y2,&z2);
      if (t==ty) {
        hist[ (int) (sqrt(x2+y2+z2)/dr) ]++;
        x_hist[ (int) (sqrt(x2)/dr) ]++;
        y_hist[ (int) (sqrt(y2)/dr) ]++;
        z_hist[ (int) (sqrt(z2)/dr) ]++;
      }
    }
  }
  fclose(infile);

  /* write histograms */
  sprintf(line,"%s.%d.hist",argv[1],ty);
  if (NULL==(outfile=fopen(line,"w"))) error("Cannot open output file");
  if (zmax==0.0) { /* 2d */
    for (i=0; i<size; i++)
      fprintf(outfile,"%f %d %d %d\n",
              i*dr,hist[i],x_hist[i],y_hist[i]);
  } else { /* 3d */
    for (i=0; i<size; i++)
      fprintf(outfile,"%f %d %d %d %d\n",
              i*dr,hist[i],x_hist[i],y_hist[i],z_hist[i]);
  }
  fclose(outfile);  

  return 0;  
}

