
#define MAX(x,y) (((x)<(y)) ? (y) : (x))
#define PIXEL(pic,i,j,ioff,joff) \
  (3 * ( (pic)->dimx * ((j) + (joff)) + (i) +(ioff) ) )

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef struct {
  int dimx, dimy, size;
  unsigned char *data;
} ppm_t;

void usage(char *progname)
{
  printf("\n");
  printf("   Usage: %s [options] <pic1.ppm> <pic2.ppm> <outpic.ppm>\n\n", 
         progname);
  printf("   Copies pic2 on top of pic1 on a canvas sufficiently big to accomodate both\n\n");
  printf("   Options:  -offset <xoff> <yoff>  Offsets of pic2 in pic1 (can be negative)\n");
  printf("             -bg     <r> <g> <b>    Background color of canvas (default white)\n");
  printf("                                    Color values must be in 0..255\n");
  printf("             -h                     This help\n\n");
}

void error(char *msg)
{
  printf("Error: %s\n", msg);
  exit(2);
}

void error_str(char *msg, char *str)
{
  char buf[255];
  sprintf(buf, msg, str);
  error(buf);
}

void error_str_str(char *msg, char *str1, char *str2)
{
  char buf[255];
  sprintf(buf, msg, str1, str2);
  error(buf);
}

void read_ppm( char *filename, ppm_t *pic )
{
  FILE *infile;
  char str[255];
  int format, depth, nread=0;

  /* open input file */
  if (NULL == (infile = fopen( filename, "r" )))
    error_str( "Cannot open input file %s", filename );

  /* read ppm header */
  do fgets(str, 255, infile); while (str[0]=='#');
  if      ((str[0]=='P') && (str[1]=='3')) { format = 3; }
  else if ((str[0]=='P') && (str[1]=='6')) { format = 6; }
  else error_str_str( "Unknown file format %s of file %s", str, filename );
  do fgets(str, 255, infile); while (str[0]=='#');
  sscanf(str, "%d %d", &(pic->dimx), &(pic->dimy) );
  do fgets(str, 255, infile); while (str[0]=='#');
  sscanf(str, "%d", &depth);
  if (255 != depth) {
    error("Can deal only with 8 bit data");
  }

  /* alloc data */
  pic->size = 3 * pic->dimx * pic->dimy;
  pic->data = (unsigned char *) malloc( pic->size );
  if (NULL == pic->data) error("Cannot allocate data");

  /* read data */
  if (format == 6) {
    nread = fread( pic->data, 1, pic->size, infile ); 
  } else if (format == 3) {
    do {
      fscanf( infile, "%d", pic->data + nread );
      nread++;
    } while ((!feof(infile)) || (nread < pic->size));
  }
  if (nread < pic->size) error_str("Premature end of file %s", filename);

  /* close input file */
  fclose(infile);
}

void write_ppm( char *filename, ppm_t *pic )
{
  FILE *out;

  /* open output file */
  if (NULL == (out = fopen( filename, "w" )))
    error_str( "Cannot open output file %s", filename );

  /* write ppm header */
  fprintf( out, "P6\n%d %d\n255\n", pic->dimx, pic->dimy );

  /* write data */
  fwrite( pic->data, 1, 3 * pic->dimx * pic->dimy, out ); 

  /* close output file */
  fclose(out);
}

void combine_ppm( ppm_t *pic, ppm_t *pic1, ppm_t *pic2, int xoff, int yoff,
                  unsigned char *bg)
{
  int i, j, ind, ind1, ind2, xoff1=0, yoff1=0, xoff2=0, yoff2=0;

  /* prepare canvas */
  if (xoff < 0) xoff1 = -xoff; else xoff2 = xoff;
  if (yoff < 0) yoff1 = -yoff; else yoff2 = yoff;
  pic->dimx = MAX( pic1->dimx + xoff1, pic2->dimx + xoff2 );
  pic->dimy = MAX( pic1->dimy + yoff1, pic2->dimy + yoff2 );
  pic->size = 3 * pic->dimx * pic->dimy;
  if (NULL == (pic->data = (unsigned char *) malloc( pic->size )))
    error("Cannot allocate data");

  /* put background color */
  for (i=0; i < pic->dimx * pic->dimy; i++) {
    pic->data[3*i  ] = bg[0];
    pic->data[3*i+1] = bg[1];
    pic->data[3*i+2] = bg[2];
  }

  /* copy pic1 on canvas */
  for (j=0; j<pic1->dimy; j++)
    for (i=0; i<pic1->dimx; i++) {
      ind  = PIXEL(pic, i,j,xoff1,yoff1);
      ind1 = PIXEL(pic1,i,j,0,    0    );
      pic->data[ind  ] = pic1->data[ind1  ];
      pic->data[ind+1] = pic1->data[ind1+1];
      pic->data[ind+2] = pic1->data[ind1+2];
    }

  /* copy pic2 on canvas */
  for (j=0; j<pic2->dimy; j++)
    for (i=0; i<pic2->dimx; i++) {
      ind  = PIXEL(pic, i,j,xoff2,yoff2);
      ind2 = PIXEL(pic2,i,j,0,    0    );
      pic->data[ind  ] = pic2->data[ind2  ];
      pic->data[ind+1] = pic2->data[ind2+1];
      pic->data[ind+2] = pic2->data[ind2+2];
    }
}

int main( int argc, char **argv ) 
{
  ppm_t pic, pic1, pic2;
  int xoff, yoff;
  unsigned char bg[3] = {255, 255, 255};
  char *progname;

  /* parse command line options */
  progname = strdup(argv[0]);
  while ((argc > 1) && (argv[1][0] =='-')) {
    if (strcasecmp(argv[1],"-offset" )==0) { 
      if (argc<4) {
        printf("\n   Not enough arguments!\n");
        usage(progname);
        exit(-1);
      }
      xoff = atoi(argv[2]);
      yoff = atoi(argv[3]);
      argc -= 3;
      argv += 3;
    }
    else if (strcasecmp(argv[1],"-bg" )==0) { 
      if (argc<5) {
        printf("\n   Not enough arguments!\n");
        usage(progname);
        exit(-1);
      }
      bg[0] = atoi(argv[2]);
      bg[1] = atoi(argv[3]);
      bg[2] = atoi(argv[4]);
      argc -= 4;
      argv += 4;
    }
  }
  if (argc<4) {
    printf("\n   Not enough arguments!\n");
    usage(progname);
    exit(-1);
  }

  /* read, process and write pictures */
  read_ppm(argv[1], &pic1);
  read_ppm(argv[2], &pic2);
  combine_ppm(&pic, &pic1, &pic2, xoff, yoff, bg);
  write_ppm(argv[3], &pic );

  return 0;
}
