
#include "dist_tools.c"

void usage(char *progname)
{
  printf("\n");
  printf("   Usage: %s <infile> <n>\n\n", progname);
  printf("   Extracts data entry <n> from distribution\n\n");
}

int main( int argc, char **argv ) 
{
  int i, j;
  float_dist_t fl1, fl2;
  char *tmp, *tmp2, *infile, outfile[255];

  if ((argc<3) || (argv[1][0]=='-')) {
    usage(argv[0]);
    exit(-1);
  }

  infile = strdup(argv[1]);
  j      = atoi(argv[2]);

  /* read distribution */
  init_float_dist( &fl1 );
  if (-1==read_float_dist( &fl1, argv[1] )) error("Cannot read distribution");

  j = atoi(argv[2]);
  if (j>=fl1.n) error("Not enough components in distribution");

  init_float_dist( &fl2 );
  fl2.dim  = fl1.dim;
  fl2.dimx = fl1.dimx; 
  fl2.sx   = fl1.sx;
  fl2.dimy = fl1.dimy; 
  fl2.sy   = fl1.sy;
  if (fl1.dim==3) {
    fl2.dimz = fl1.dimz; 
    fl2.sz   = fl1.sz;
  }
  fl2.nbin = fl1.nbin;
  fl2.len  = fl1.nbin;
  fl2.n    = 1;

  /* (re)allocate distribution */
  alloc_float_dist( &fl2 );  
  fl2.cont[0] = fl1.cont[j];

  /* copy the data */
  for (i=0; i<fl1.nbin; i++) fl2.dat[i] = fl1.dat[fl1.n*i+j];

  /* construct outfile name, replacing last suffix by name of data entry */
  tmp = infile;
  do {
    tmp2 = strstr( tmp, ".");
    if (tmp2) tmp=tmp2+1;
  } while (tmp2);
  if (tmp!=infile) *tmp='\0';
  sprintf(outfile, "%s%s", infile, fl1.cont[j]);

  /* write the distribution */
  write_float_dist(&fl2, outfile);

  return 0;
}
