
#include "dist_tools.c"

/* supported output formats */
#define UV8   1    /* scalar  8-bit - UltraVis, VolRend */
#define UV12  2    /* scalar 12-bit - UltraVis          */
#define RVF   3    /* scalar  8-bit - Virvo (no movies) */
#define XVF   4    /* scalar  8-bit - Virvo             */ 

void usage(char *progname)
{
  printf("\n");
  printf("   Usage: %s [options] <infilebase> <suffix> \n\n", progname);
  printf("   Options:  -l <llx> <lly> <llz>  lower left  corner\n");
  printf("             -u <urx> <ury> <urz>  upper right corner\n");
  printf("             -i <min> <max>        data cutuff interval\n");
  printf("             -f <min> <max>        range of time steps\n");
  printf("             -v <speed>            speed (seconds/frame)\n");
  printf("             -s <int>              time step interval\n");
  printf("             -t uv8|uv12|rvf|xvf   output format\n");
  printf("             -n <n>                index of data slot\n");
  printf("             -h                    this help\n\n");
}

void  write_header_uv( int dimx, int dimy, int dimz, int nframe, 
  char *file, char *voxeltype) 
{
  FILE *fp;
  char fname[255];

  sprintf(fname, "%s.uhd", file );
  if (NULL==(fp=fopen(fname, "w"))) error("Cannot open parameter file");

  fprintf(fp, "HP Labs UltraVis Volume Data Set Header File, Hewlett-Packard Company, 2000\r\n");
  fprintf(fp, "VERSION: 1\r\n");
  fprintf(fp, "DATASETNAME: %s.uvd\r\n", file );
  fprintf(fp, "MANUFACTURER: IMD\r\n");
  fprintf(fp, "COPYRIGHT: ITAP, Stuttgart University\r\n");
  fprintf(fp, "DESC: %s\r\n", file);
  fprintf(fp, "MENUNAME: %s\r\n", file);
  fprintf(fp, "XSIZE: %d\r\n", dimx);
  fprintf(fp, "YSIZE: %d\r\n", dimy);
  fprintf(fp, "ZSIZE: %d\r\n", dimz);
  fprintf(fp, "VOXELTYPE: %s\r\n", voxeltype);
  fprintf(fp, "TIMESTEPS: %d\r\n", nframe);
  fprintf(fp, "VOLUMEPATH: %s.uvd\r\n", file);
  fprintf(fp, "COMPRESSION: NONE\r\n");
  fprintf(fp, "GEOPATH:\r\n");
  fprintf(fp, "PARAMPATH: %s.ups\r\n", file);
  fprintf(fp, "USERLIBRARY:\r\n");
  fprintf(fp, "COMMENTS:\r\n");

  fclose(fp);
}

void write_header_rvf( FILE *fp, int dx, int dy, int dz)
{
  char header[6], *str;
  unsigned short us;

  str = header;
  us = (unsigned short) dx; copy_2_bytes(str,&us); str +=2; /* dim_x */
  us = (unsigned short) dy; copy_2_bytes(str,&us); str +=2; /* dim_y */
  us = (unsigned short) dz; copy_2_bytes(str,&us); str +=2; /* dim_y */
  if (6!=fwrite(header, sizeof(char), 6, fp)) error("Cannot write header!");
}

void write_header_xvf( FILE *fp, float_dist_t *fl, 
  int dx, int dy, int dz, int nf, float min, float max, float v)
{
  char header[69], *str;
  unsigned short us;
  unsigned int ui;
  unsigned char uc;
  float f;

  /* make volume file header */
  str = header;
  sprintf(str,"%s","VIRVO-XVF");      str +=9;   /* file type */
  us = 69; copy_2_bytes(str,&us);     str +=2;   /* header size */
  ui = dx; copy_4_bytes(str,&ui);     str +=4;   /* dim_x */
  ui = dy; copy_4_bytes(str,&ui);     str +=4;   /* dim_y */
  ui = dz; copy_4_bytes(str,&ui);     str +=4;   /* dim_z */
  ui = nf; copy_4_bytes(str,&ui);     str +=4;   /* number of frames */
  uc =  8; (unsigned char) *str = uc;   str++;   /* bits per voxel */
  f  = fl->sx; copy_4_bytes(str,&f);  str +=4;   /* x-length of voxel */
  f  = fl->sy; copy_4_bytes(str,&f);  str +=4;   /* y-length of voxel */
  f  = fl->sz; copy_4_bytes(str,&f);  str +=4;   /* z-length of voxel */
               copy_4_bytes(str,&v);  str +=4;   /* secs per frame */
  f  = min; copy_4_bytes(str,&f);     str +=4;   /* minimum data value */
  f  = max; copy_4_bytes(str,&f);     str +=4;   /* maximum data value */
  f  = 0.0; copy_4_bytes(str,&f);     str +=4;   /* x-pos of volume center */
  f  = 0.0; copy_4_bytes(str,&f);     str +=4;   /* y-pos of volume center */
  f  = 0.0; copy_4_bytes(str,&f);     str +=4;   /* z-pos of volume center */
  uc =   0; (unsigned char) *str = uc;  str++;   /* no run length encoding */
  us =   0; copy_2_bytes(str,&us);    str +=2;   /* number of transf. func. */
  us =   0; copy_2_bytes(str,&us);    str +=2;   /* type of transf. func. */
  if (69!=fwrite(header, sizeof(char), 69, fp)) error("Cannot write header!");
}

void write_header_xvf_sep( char *basename, float_dist_t *fl, 
  int dx, int dy, int dz, int nf, float min, float max, float v)
{
  FILE *fp;
  char fname[255];

  sprintf(fname, "%s.xvf-header", basename );
  if (NULL==(fp=fopen(fname, "w"))) error("Cannot open parameter file");
  write_header_xvf( fp, fl, dx, dy, dz, nf, min, max, v);
  fclose(fp);
}

int main( int argc, char **argv ) 
{
  char *progname, *infilebase, *suffix, *ext;
  char tmpstr[255], infile[255], outfile[255];
  FILE *fp;
  int have_ur=0, have_minmax=0, have_type=0, type, step=1, nsteps;
  int llx=0, lly=0, llz=0, urx, ury, urz, fmin=0, fmax=0, n=0, i;
  float min, max, v=0.04;
  float_dist_t fl;
  byte_dist_t bt;

  /* parse command line options */
  progname = strdup(argv[0]);
  while ((argc > 1) && (argv[1][0] =='-')) {
    if (argv[1][1]=='l') {
      llx = atoi(argv[2]);
      lly = atoi(argv[3]);
      llz = atoi(argv[4]);
      argc -= 4;
      argv += 4;
    }
    else if (argv[1][1]=='u') {
      have_ur=1;
      urx = atoi(argv[2]);
      ury = atoi(argv[3]);
      urz = atoi(argv[4]);
      argc -= 4;
      argv += 4;
    }
    else if (argv[1][1]=='i') {
      have_minmax=1;
      min = atof(argv[2]);
      max = atof(argv[3]);
      argc -= 3;
      argv += 3;
    }
    else if (argv[1][1]=='f') {
      fmin = atoi(argv[2]);
      fmax = atoi(argv[3]);
      argc -= 3;
      argv += 3;
    }
    else if (argv[1][1]=='s') {
      step = atoi(argv[2]);
      argc -= 2;
      argv += 2;
    }
    else if (argv[1][1]=='v') {
      v     = atof(argv[2]);
      argc -= 2;
      argv += 2;
    }
    else if (argv[1][1]=='t') {
      have_type=1;
      if      (strcasecmp(argv[2],"uv8" )==0) { type = UV8;  ext="uvd"; }
      else if (strcasecmp(argv[2],"uv12")==0) { type = UV12; ext="uvd"; }
      else if (strcasecmp(argv[2],"rvf" )==0) { type = RVF;  ext="rvf"; }
      else if (strcasecmp(argv[2],"xvf" )==0) { type = XVF;  ext="xvf"; }
      else {
        printf("Unsopported output type %s", argv[2]);
        usage(progname);
        exit(-1);
      }
      argc -= 2;
      argv += 2;
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
  if ((fmin!=fmax) && (type==RVF)) 
    error("Virvo format rvf supports only single images"); 
  if (argc<3) {
      usage(progname);
      exit(-1);
  }

  infilebase = strdup(argv[1]);
  suffix     = strdup(argv[2]);

  /* read the first distribution */
  init_float_dist( &fl );
  init_byte_dist(  &bt );
  sprintf(infile, "%s.%d.%s", infilebase, fmin, suffix);
  if (-1==read_float_dist( &fl, infile )) error("Cannot read distribution!");
  if (fl.dim!=3) error("Need 3D distribution!");
  if (n>=fl.n) error("Not enough components in distribution!");

  /* set some defaults */
  if (!have_ur) {
    urx = fl.dimx;
    ury = fl.dimy;
    urz = fl.dimz;
  }
  if (!have_type) {
    type = UV8;
    ext = "uvd";
  }
  nsteps = (fmax-fmin+step) / step;

  /* open output file, write header */
  sprintf(tmpstr, "%s.%s", infilebase, fl.cont[n] );
  sprintf(outfile, "%s.%s", tmpstr, ext);
  fp = fp=fopen(outfile,"w");
  if (NULL==fp) error("Cannot open output file");
  if (type==UV8) {
    write_header_uv(urx-llx,ury-lly,urz-llz,nsteps,tmpstr,"SCALAR8" );
    write_header_xvf_sep(tmpstr,&fl,urx-llx,ury-lly,urz-llz,nsteps,min,max,v);
  } else if (type==UV12) { 
    write_header_uv(urx-llx,ury-lly,urz-llz,nsteps,tmpstr,"SCALAR12");
  } else if (type==RVF) {
    write_header_rvf(fp, urx-llx, ury-lly, urz-llz);
  } else if (type==XVF ) {
    write_header_xvf(fp, &fl, urx-llx, ury-lly, urz-llz, nsteps, min, max, v);
  }

  /* loop over frames */
  for (i=fmin; i<=fmax; i+=step) {

    /* defaults if min and max are not specified */
    if (!have_minmax) {
      min = fl.min[n];
      max = fl.max[n];
    }

    /* convert the distribution */
    if ((type==UV8) || (type==RVF) || (type==XVF)) {
      float2scalar8_dist( &fl, &bt, min, max, n, llx, lly, llz, urx, ury, urz);
    } else if (type==UV12) {
      float2scalar12_dist(&fl, &bt, min, max, n, llx, lly, llz, urx, ury, urz);
    }

    /* write the distribution */
    fwrite( bt.dat, sizeof(unsigned char), bt.len , fp ); 

    /* read next distribution */
    if (i<=fmax-step) {
      sprintf(infile, "%s.%d.%s", infilebase, i+step, suffix);
      if (-1==read_float_dist( &fl,infile)) error("Cannot read distribution!");
    }

  }
  fclose(fp);

  return 0;
}
