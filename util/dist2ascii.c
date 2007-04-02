
#include "dist_tools.c"

void usage(char *progname)
{
  printf("\n");
  printf("   Usage: %s <infile>\n\n", progname);
  printf("   Writes ASCII version of distribution in <infile> to stdout\n\n");
}

int main( int argc, char **argv ) 
{
  int ix, iy, iz, i, j;
  float_dist_t fl;

  if ((argc<2) || (argv[1][0]=='-')) {
    usage(argv[0]);
    exit(-1);
  }

  /* read distribution */
  init_float_dist( &fl );
  if (-1==read_float_dist( &fl, argv[1] )) error("Cannot read distribution");

  /* write as ascii to stdout */

  /* first the header */
  printf("#F A %d %d %d\n", fl.dim, fl.dim, fl.n);
  printf("#C");
  for (i=0; i<fl.n; i++) printf(" %s", fl.cont[i]);
  printf("\n");
  if (fl.dim==2)
    printf("#D %d %d\n", fl.dimx, fl.dimy);
  else 
    printf("#D %d %d %d\n", fl.dimx, fl.dimy, fl.dimz);
  if (fl.dim==2)
    printf("#S %f %f\n", fl.sx, fl.sy);
  else 
    printf("#S %f %f %f\n", fl.sx, fl.sy, fl.sz);
  printf("#E\n");

  /* then the data */
  for (ix=0; ix<fl.dimx; ix++)
    for (iy=0; iy<fl.dimy; iy++)
      for (iz=0; iz<fl.dimz; iz++) {
        if (fl.dim==2) printf("%d %d",    ix, iy);
        else           printf("%d %d %d", ix, iy, iz);
        for (j=0; j<fl.n; j++) printf(" %e", ELM(&fl,ix,iy,iz,j) );
        printf("\n");
      }
  return 0;
}
