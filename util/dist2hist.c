
#define NBIN 20

#include "dist_tools.c"

void usage(char *progname)
{
  printf("\n");
  printf("   Usage: %s <infile> [<n>]\n\n", progname);
  printf("   Writes histogram of entry <n> (default 0) of distribution to stdout\n\n");
}

int main(int argc, char **argv) 
{
  float_dist_t fl;
  float d, min;
  int hist[NBIN];
  int n=0, i, j, k, ix, iy, iz;

  if ((argc<2) || (argv[1][0]=='-')) {
    usage(argv[0]);
    exit(-1);
  }

  init_float_dist(&fl);
  if (-1==read_float_dist(&fl,argv[1])) error("cannot read input file");

  for (i=0; i<NBIN; i++) hist[i]=0;

  if (argc==3) n = atoi(argv[2]);
  if (n>=fl.n) error("Not enough components in distribution!");

  printf("# Histogram of %s\n", fl.cont[n]);
  d = (NBIN-1) * 1.0001 / (fl.max[n]-fl.min[n]);  
  min = fl.min[n];

  for (ix=0; ix<fl.dimx; ix++)
    for (iy=0; iy<fl.dimy; iy++)
      for (iz=0; iz<fl.dimz; iz++) {
        j = (int) (d * (ELM(&fl,ix,iy,iz,n) - min));
        hist[j] += 1;
      }

  for (i=0; i<NBIN; i++) {
    printf("%e %d\n", min + i / d, hist[i]);
  }
  return 0;
}

