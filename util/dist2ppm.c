
#include "dist_tools.c"

void usage(char *progname)
{
  printf("\n");
  printf("   Usage: %s [options] <infile>\n\n", progname);
  printf("   Options:  -l <llx> <lly>   lower left  corner\n");
  printf("             -u <urx> <ury>   upper right corner\n");
  printf("             -i <min> <max>   data cutuff interval\n");
  printf("             -n <n>           index of data slot\n");
  printf("             -h               this help\n\n");
}

int main( int argc, char **argv ) 
{
  char *progname, *infile, outfile[255];
  char *tmp, *tmp2;
  int have_ur=0, have_minmax=0;
  int llx=0, lly=0, urx, ury, n=0;
  float min, max;
  float_dist_t fl;
  byte_dist_t bt;

  /* parse command line options */
  progname = strdup(argv[0]);
  while ((argc > 1) && (argv[1][0] =='-')) {
    if (argv[1][1]=='l') {
      llx = atoi(argv[2]);
      lly = atoi(argv[3]);
      argc -= 3;
      argv += 3;
    }
    else if (argv[1][1]=='u') {
      have_ur=1;
      urx = atoi(argv[2]);
      ury = atoi(argv[3]);
      argc -= 3;
      argv += 3;
    }
    else if (argv[1][1]=='i') {
      have_minmax=1;
      min = atof(argv[2]);
      max = atof(argv[3]);
      argc -= 3;
      argv += 3;
    }
    else if (argv[1][1]=='n') {
      n = atoi(argv[2]);
      argc -= 2;
      argv += 2;
    }
    else if (argv[1][1]=='h') {
      usage(progname);
      exit(-1);
    }
    else {
      printf("\n   Illegal option %s \n", argv[1]);
      usage(progname);
      exit(-1);
    }
  }
  if (argc<2) {
      usage(progname);
      exit(-1);
  }

  infile = strdup(argv[1]);

  /* read distribution */
  init_float_dist( &fl );
  if (-1==read_float_dist( &fl, infile )) error("Cannot read distribution");
  if (fl.dim!=2) error("Need 2D distribution!");
  if (n>=fl.n) error("Not enough data entries!");

  /* defaults for options not specified */
  if (!have_ur) {
    urx = fl.dimx;
    ury = fl.dimy;
  }
  if (!have_minmax) {
    min = fl.min[n];
    max = fl.max[n];
  }

  /* convert to rgb distribution */
  init_byte_dist( &bt );
  float2ppm_dist( &fl, &bt, min, max, n, llx, lly, urx, ury, rgbcolor);

  /* construct outfile name, replacing last suffix by name of data entry */
  tmp = infile;
  do {
    tmp2 = strstr( tmp, ".");
    if (tmp2) tmp=tmp2+1;
  } while (tmp2);
  if (tmp!=infile) *tmp='\0';
  sprintf(outfile, "%s%s.ppm", infile, fl.cont[n]);

  /* write ppm file */
  write_ppm( &bt, outfile );

  return 0;
}

