/******************************************************************************
*
*  sample2dist is a utility program to convert checkpoint files to
*  atoms distribution files.
*
*  Compilation:
*
*    gmake sample2dist
*
*  Usage:
*
*     sample2dist <dimx> <dimy> <dimz> <chkpt-file> [--normalize]
*
*  dimx, etc. are the distribution sizes in the three directions
*  if --normalize is used the bins are normalize by the number of atoms
*  they contain
*
*  Two files <chkpt-file>.pot.dist and ...kin.dist are created
*
******************************************************************************/              
#include <stdio.h>
#include <string.h>
#include <math.h>

int main(int argc, char **argv) {
  FILE *fp, *pot, *kin;
  char line[255], str[255], potstr[255], kinstr[255];
  float m, x, y, z, vx, vy, vz, e, boxx, boxy, boxz;
  int n, t, i, j, k, pos, sx, sy, sz, size, norm = 0;
  int *anz;
  float *smpp, *smpk;

  /* Usage */
  if ((argc!=5)&&(argc!=6)) {
    fprintf(stderr, "Usage: ./sampledist dimx dimy dimz chkpt-file [--normalize]!");
    exit(-1);
  }

  /* Dimension of dist */
  sscanf(argv[1],"%d",&sx);
  sscanf(argv[2],"%d",&sy);
  sscanf(argv[3],"%d",&sz);
  size = sx*sy*sz;

  /* allocation */
  anz = (int *)calloc(sx*sy*sz, sizeof(int));
  smpp = (float *)calloc(sx*sy*sz, sizeof(float));
  smpk = (float *)calloc(sx*sy*sz, sizeof(float));

  /* normalize */
  if ((argc==6)&&(!strcmp(argv[5], "--normalize")))
      norm = 1;

  /* init */
  for (i=0;i<sx*sy*sz;i++) {
      anz[i] = 0;
      smpp[i] = 0.0;
      smpk[i] = 0.0;
  }

  /* read chkpt file */
  fp=fopen(argv[4], "r");
  /* header */
  fgets(line, 255, fp);
  fgets(line, 255, fp);
  fgets(line, 255, fp);
  sscanf(line, "%s %f %f %f", str, &boxx, &boxy, &boxz);
  fgets(line, 255, fp);
  sscanf(line, "%s %f %f %f", str, &boxy, &boxy, &boxz);
  fgets(line, 255, fp);
  sscanf(line, "%s %f %f %f", str, &boxz, &boxz, &boxz);
  fgets(line, 255, fp);
  fgets(line, 255, fp);

  /* atom lines */
  while(fgets(line, 255, fp)) {
    sscanf(line, "%d %d %f %f %f %f %f %f %f %f\n",
	   &n, &t, &m, &x, &y, &z, &vx, &vy, &vz, &e);

    /* Coordinates of bin */
    i = (int)floor(sx*x/boxx);
    j = (int)floor(sy*y/boxy);
    k = (int)floor(sz*z/boxz);
    i = (i>0)?i:0;
    j = (j>0)?j:0;
    k = (k>0)?k:0;
    i = (i<sx)?i:sx-1;
    j = (j<sy)?j:sy-1;
    k = (k<sz)?k:sz-1;
    pos = i*sy*sz + j*sz + k;

    /* values into bin */
    smpp[pos] += e;
    smpk[pos] += vx*vx+vy*vy+vz*vz;
    anz[pos]++;
  }
  fclose(fp);
  
  /* output, 1st header */
  sprintf(potstr, "%s.pot.dist", argv[4]);
  sprintf(kinstr, "%s.kin.dist", argv[4]);
  pot = fopen(potstr, "w");
  kin = fopen(kinstr, "w");

  fprintf(pot, "#F A 3 0 1\n");
  fprintf(pot, "#C Epot\n");
  fprintf(pot, "#D %d %d %d\n", sx, sy, sz);
  fprintf(pot, "#S %f %f %f\n", boxx/(float)sx, boxy/(float)sy, boxz/(float)sz);
  fprintf(pot, "#E\n");
  fprintf(kin, "#F A 3 0 1\n");
  fprintf(kin, "#C Ekin\n");
  fprintf(kin, "#D %d %d %d\n", sx, sy, sz);
  fprintf(kin, "#S %f %f %f\n", boxx/(float)sx, boxy/(float)sy, boxz/(float)sz);
  fprintf(kin, "#E\n");

  /* dist */
  for (i=0;i<sx*sy*sz;i++)
      if (norm) {
	  fprintf(pot, "%f\n", smpk[i]/anz[i]);
	  fprintf(kin, "%f\n", smpp[i]/anz[i]);
      } else {
	  fprintf(pot, "%f\n", smpp[i]);
	  fprintf(kin, "%f\n", smpk[i]);
      }

  fclose(pot);
  fclose(kin);
}


