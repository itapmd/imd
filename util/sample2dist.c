#include <stdio.h>
#include <math.h>

int main(int argc, char **argv) {
  FILE *fp;
  char line[255], str[100];
  float m, x, y, z, vx, vy, vz, e, corr, boxx, boxy, boxz;
  int n, t, i, j, k, sx, sy, sz;
  int anz[50][50][12];
  float smp[50][50][12];

  if (argc!=5) {
    fprintf(stderr, "Usage: ./sampledist dimx dimy dimz chkpt-file!");
    exit(-1);
  }

  sscanf(argv[1],"%d",&sx);
  sscanf(argv[2],"%d",&sy);
  sscanf(argv[3],"%d",&sz);

  for (i=0;i<sx;i++)
    for (j=0;j<sy;j++)
      for (k=0;k<sz;k++) {
	anz[i][j][k] = 0;
	smp[i][j][k] = 0.0;
      }

  fp=fopen(argv[4], "r");
  /* header */
  fgets(line, 255, fp);
  fgets(line, 255, fp);
  fgets(line, 255, fp);
  sscanf(line, "%s %s %f %f %f", str, str, &boxx, &boxy, &boxz);
  fgets(line, 255, fp);
  sscanf(line, "%s %s %f %f %f", str, str, &corr, &boxy, &boxz);
  fgets(line, 255, fp);
  sscanf(line, "%s %s %f %f %f", str, str, &boxz, &boxz, &boxz);
  fgets(line, 255, fp);
  fgets(line, 255, fp);

  while(fgets(line, 255, fp)) {
    sscanf(line, "%d %d %f %f %f %f %f %f %f %f\n",
	   &n, &t, &m, &x, &y, &z, &vx, &vy, &vz, &e);
    x -= corr/boxx*y;
    i = (int)floor(sx*x/boxx);
    j = (int)floor(sy*y/boxy);
    k = (int)floor(sz*z/boxz);
    i = (i>0)?i:0;
    j = (j>0)?j:0;
    k = (k>0)?k:0;
    i = (i<sx)?i:sx-1;
    j = (j<sy)?j:sy-1;
    k = (k<sz)?k:sz-1;
    smp[i][j][k] += vx*vx+vy*vy+vz*vz;
    anz[i][j][k]++;
  }
  fclose(fp);
  
  for (i=0;i<sx;i++)
    for (j=0;j<sy;j++)
      for (k=0;k<sz;k++)
	if (anz[i][j][k]>0)
	  printf("%d %d %d %f\n", i, j, k, smp[i][j][k]/anz[i][j][k]);
	else
	  printf("%d %d %d 0\n", i, j, k);
}


