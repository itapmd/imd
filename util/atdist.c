#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))

float *atoms_dist;
int dim, dimx, dimy, dimz, ntypes, size;

/******************************************************************************
*
*  return error message and stop
*
******************************************************************************/

void error(char *msg)
{
  fprintf(stderr,"Error: %s\n",msg);
  exit(2);
}

/******************************************************************************
*
*  endian returns 1 if system is big endian, 0 if little endian
*
******************************************************************************/

int endian(void)
{
  unsigned short int word = 0x0001;
  unsigned char  *byte    = (unsigned char *) &word;
  return (byte[0] ? 0 : 1);
}

/******************************************************************************
*
*  do endian swap for a four byte quantity
*
******************************************************************************/

void do_endian_swap_4(char *str) {
  char *c;
  *c = * str;    * str    = *(str+3); *(str+3) = *c;
  *c = *(str+1); *(str+1) = *(str+2); *(str+2) = *c;
}

/******************************************************************************
*
*  read atoms distribution
*
******************************************************************************/

int read_atoms_dist(char *fname)
{
  FILE *infile;
  char line[255];
  int  cont, input_endian, i;

  /* open file */
  if (NULL==(infile=fopen(fname,"r"))) error("Cannot open input file");

  /* read file header */
  cont=1;
  do {
    fgets(line,255,infile);
    if (line[0]!='#') {
      error("file header corrupt!");
    } else {
      if (line[1]=='F') {
        sscanf(line+2,"%d %d %d",&input_endian,&dim,&ntypes);
      } else  if (line[1]=='D') {
        if (dim!=sscanf(line+2,"%d %d %d",&dimx,&dimy,&dimz))
          error("file header corrupt (dimension)!");
      } else  if (line[1]=='E') {
	cont=0;
      }
    }
  } while (cont);
  
  /* compute array size */
  size = dimx * dimy;
  if (dim==3) size *= dimz;
  size *= ntypes;

  /* allocate atoms distribution array */
  atoms_dist = (float *) malloc( size * sizeof(float) );
  if (NULL==atoms_dist) error("Cannot allocate atoms distribution array");

  /* read atoms distribution data */
  if (size!=fread(atoms_dist,sizeof(float),size,infile))
    error("Cannot read histogram data");

  /* close input file */
  fclose(infile);

  /* do endian swap if necessary */
  if (input_endian!=endian()) {
    for (i=0; i<size; i++) do_endian_swap_4((char *)(atoms_dist+i));
  }

  return dim;
}

/******************************************************************************
*
*  project histogram on the three main axis
*
******************************************************************************/

void axis_projections_3d()
{
  float *histx, *histy, *histz;
  int dimxy, dimxyz, t, i, j, k, l;

  histx = (float *) calloc( ntypes * dimx, sizeof(float) );
  histy = (float *) calloc( ntypes * dimy, sizeof(float) );
  histz = (float *) calloc( ntypes * dimz, sizeof(float) );
  if ((histx==NULL) || (histy==NULL) || (histx==NULL))
    error("out of memory");

  dimxy  = dimx  * dimy;
  dimxyz = dimxy * dimz;

  for (t=0; t<ntypes; t++)
    for (i=0; i<dimx; i++)
      for (j=0; j<dimy; j++)
        for (k=0; k<dimz; k++) {
          l = t * dimxyz + i * dimxy + j * dimx + k;
          histx[ntypes * dimx + i] += atoms_dist[l];
          histy[ntypes * dimy + j] += atoms_dist[l];
          histz[ntypes * dimz + k] += atoms_dist[l];
	}

  printf("# Projection on x-Axis:\n");
  for (i=0; i<dimx; i++) {
    printf("%d",i);
    for (t=0; t<ntypes; t++) printf(" %e", histx[ntypes * dimx + i]);
    printf("\n");
  }

  printf("\n\n# Projection on y-Axis:\n");
  for (i=0; i<dimy; i++) {
    printf("%d",i);
    for (t=0; t<ntypes; t++) printf(" %e", histy[ntypes * dimy + i]);
    printf("\n");
  }

  printf("\n\n# Projection on z-Axis:\n");
  for (i=0; i<dimz; i++) {
    printf("%d",i);
    for (t=0; t<ntypes; t++) printf(" %e", histz[ntypes * dimz + i]);
    printf("\n");
  }
}

/******************************************************************************
*
*  write pictures of xy-slices
*
******************************************************************************/

void xy_pictures_3d(char *infile, int min, int max)
{
  float *hist, fmax;
  char  *pix, *pix2;
  int dimxy, dimxyz, t, i, j, k, l;
  FILE *out;
  char fname[255]; 

  hist = (float *) calloc( ntypes * dimx * dimy, sizeof(float) );
  pix  = (char  *) calloc( ntypes * dimx * dimy, sizeof(char ) );
  pix2 = (char  *) calloc( 3      * dimx * dimy, sizeof(char ) );
  if ((hist==NULL) || (pix==NULL) || (pix2==NULL)) error("out of memory");

  dimxy  = dimx  * dimy;
  dimxyz = dimxy * dimz;

  for (t=0; t<ntypes; t++)
    for (i=0; i<dimx; i++)
      for (j=0; j<dimy; j++)
        for (k=min; k<max+1; k++) {
          l = t * dimxyz + i * dimxy + j * dimx + k;
          hist[t * dimxy + dimx * i + j] += atoms_dist[l];
	}

  /* renormalize atoms distribution for pgm files*/
  fmax = 0.0;
  for (i=0; i<ntypes*dimxy; i++) fmax = MAX(fmax,hist[i]);
  fmax = 255/fmax;
  for (i=0; i<ntypes*dimxy; i++) pix[i] = (char) (255-hist[i]*fmax);

  /* write pgm files */
  for (t=0; t<ntypes; t++) {
    sprintf(fname,"%s.xy.%d.pgm",infile,t);
    if (NULL==(out=fopen(fname,"w"))) 
      error("Cannot open pgm-file.");
    fprintf(out,"P5\n%d %d\n255\n", dimx, dimy);
    if (dimxy!=fwrite(pix+t*dimxy,sizeof(char),dimxy,out))
      error("Cannot write pgm-file.");
    fclose(out);
  }

  /* write ppm-file */
  if (ntypes<4) {
    /* set everything to white */
    for (i=0; i<3*dimxy; i++) pix2[i] = 255;
    /* get renormalization factor for color pictures */
    fmax = 0.0;
    for (i=0; i<dimxy; i++) {
      float tmp=0.0;
      for (t=0; t<ntypes; t++) tmp += hist[t*dimxy+i];
      fmax = MAX(fmax,tmp);
    }
    fmax = 255/fmax;
    /* make color picture */
    for (t=0; t<ntypes; t++) 
      for (i=0; i<dimxy; i++) 
        for (k=0; k<3; k++) 
          if (k!=t) pix2[3*i+k] -= (char) (hist[t*dimxy+i]*fmax);
    /* write color picture */
    sprintf(fname,"%s.xy.ppm",infile);
    if (NULL==(out=fopen(fname,"w"))) error("Cannot open pgm-file.");
    fprintf(out,"P6\n%d %d\n255\n", dimx, dimy);
    if (3*dimxy!=fwrite(pix2,sizeof(char),3*dimxy,out))
      error("Cannot write pgm-file.");
    fclose(out);
  }
}

/******************************************************************************
*
*  write pictures of xz-slices
*
******************************************************************************/

void xz_pictures_3d(char *infile, int min, int max)
{
  float *hist, fmax;
  char  *pix, *pix2;
  int dimxy, dimxz, dimxyz, t, i, j, k, l;
  FILE *out;
  char fname[255]; 

  hist = (float *) calloc( ntypes * dimx * dimz, sizeof(float) );
  pix  = (char  *) calloc( ntypes * dimx * dimz, sizeof(char ) );
  pix2 = (char  *) calloc( 3      * dimx * dimz, sizeof(char ) );
  if ((hist==NULL) || (pix==NULL) || (pix2==NULL)) error("out of memory");

  dimxy  = dimx  * dimy;
  dimxz  = dimx  * dimz;
  dimxyz = dimxy * dimz;

  for (t=0; t<ntypes; t++)
    for (i=0; i<dimx; i++)
      for (j=min; j<max+1; j++)
        for (k=min; k<max+1; k++) {
          l = t * dimxyz + i * dimxy + j * dimx + k;
          hist[t * dimxz + dimx * i + k] += atoms_dist[l];
	}

  /* renormalize atoms distribution for pgm files*/
  fmax = 0.0;
  for (i=0; i<ntypes*dimxz; i++) fmax = MAX(fmax,hist[i]);
  fmax = 255/fmax;
  for (i=0; i<ntypes*dimxz; i++) pix[i] = (char) (255-hist[i]*fmax);

  /* write pgm files */
  for (t=0; t<ntypes; t++) {
    sprintf(fname,"%s.xz.%d.pgm",infile,t);
    if (NULL==(out=fopen(fname,"w"))) 
      error("Cannot open pgm-file.");
    fprintf(out,"P5\n%d %d\n255\n", dimx, dimz);
    if (dimxz!=fwrite(pix+t*dimxz,sizeof(char),dimxz,out))
      error("Cannot write pgm-file.");
    fclose(out);
  }

  /* write ppm-file */
  if (ntypes<4) {
    /* set everything to white */
    for (i=0; i<3*dimxz; i++) pix2[i] = 255;
    /* get renormalization factor for color pictures */
    fmax = 0.0;
    for (i=0; i<dimxz; i++) {
      float tmp=0.0;
      for (t=0; t<ntypes; t++) tmp += hist[t*dimxz+i];
      fmax = MAX(fmax,tmp);
    }
    fmax = 255/fmax;
    /* make color picture */
    for (t=0; t<ntypes; t++) 
      for (i=0; i<dimxz; i++) 
        for (k=0; k<3; k++) 
          if (k!=t) pix2[3*i+k] -= (char) (hist[t*dimxz+i]*fmax);
    /* write color picture */
    sprintf(fname,"%s.xz.ppm",infile,t);
    if (NULL==(out=fopen(fname,"w"))) error("Cannot open pgm-file.");
    fprintf(out,"P6\n%d %d\n255\n", dimx, dimz);
    if (3*dimxz!=fwrite(pix2,sizeof(char),3*dimxz,out))
      error("Cannot write pgm-file.");
    fclose(out);
  }
}

/******************************************************************************
*
*  write pictures of yz-slices
*
******************************************************************************/

void yz_pictures_3d(char *infile, int min, int max)
{
  float *hist, fmax;
  char  *pix, *pix2;
  int dimxy, dimyz, dimxyz, t, i, j, k, l;
  FILE *out;
  char fname[255]; 

  hist = (float *) calloc( ntypes * dimy * dimz, sizeof(float) );
  pix  = (char  *) calloc( ntypes * dimy * dimz, sizeof(char ) );
  pix2 = (char  *) calloc( 3      * dimy * dimz, sizeof(char ) );
  if ((hist==NULL) || (pix==NULL) || (pix2==NULL)) error("out of memory");

  dimxy  = dimx  * dimy;
  dimyz  = dimy  * dimz;
  dimxyz = dimxy * dimz;

  for (t=0; t<ntypes; t++)
    for (i=min; i<max+1; i++)
      for (j=0; j<dimy; j++)
        for (k=0; k<dimz; k++) {
          l = t * dimxyz + i * dimxy + j * dimx + k;
          hist[t * dimyz + dimy * j + k] += atoms_dist[l];
	}

  /* renormalize atoms distribution for pgm files*/
  fmax = 0.0;
  for (i=0; i<ntypes*dimyz; i++) fmax = MAX(fmax,hist[i]);
  fmax = 255/fmax;
  for (i=0; i<ntypes*dimyz; i++) pix[i] = (char) (255-hist[i]*fmax);

  /* write pgm files */
  for (t=0; t<ntypes; t++) {
    sprintf(fname,"%s.yz.%d.pgm",infile,t);
    if (NULL==(out=fopen(fname,"w"))) 
      error("Cannot open pgm-file.");
    fprintf(out,"P5\n%d %d\n255\n", dimy, dimz);
    if (dimyz!=fwrite(pix+t*dimyz,sizeof(char),dimyz,out))
      error("Cannot write pgm-file.");
    fclose(out);
  }

  /* write ppm-file */
  if (ntypes<4) {
    /* set everything to white */
    for (i=0; i<3*dimyz; i++) pix2[i] = 255;
    /* get renormalization factor for color pictures */
    fmax = 0.0;
    for (i=0; i<dimyz; i++) {
      float tmp=0.0;
      for (t=0; t<ntypes; t++) tmp += hist[t*dimyz+i];
      fmax = MAX(fmax,tmp);
    }
    fmax = 255/fmax;
    /* make color picture */
    for (t=0; t<ntypes; t++) 
      for (i=0; i<dimyz; i++) 
        for (k=0; k<3; k++) 
          if (k!=t) pix2[3*i+k] -= (char) (hist[t*dimyz+i]*fmax);
    /* write color picture */
    sprintf(fname,"%s.yz.ppm",infile);
    if (NULL==(out=fopen(fname,"w"))) error("Cannot open pgm-file.");
    fprintf(out,"P6\n%d %d\n255\n", dimy, dimz);
    if (3*dimyz!=fwrite(pix2,sizeof(char),3*dimyz,out))
      error("Cannot write pgm-file.");
    fclose(out);
  }
}

/******************************************************************************
*
*  write pictures of 2d histogram
*
******************************************************************************/

void pictures_2d(char *infile)
{
  float *hist, fmax;
  char  *pix, *pix2;
  int dimxy, t, i, j, k, l;
  FILE *out;
  char fname[255]; 

  pix  = (char  *) calloc( ntypes * dimx * dimy, sizeof(char) );
  pix2 = (char  *) calloc( 3      * dimx * dimy, sizeof(char) );
  if ((hist==NULL) || (pix==NULL) || (pix2==NULL)) error("out of memory");
  dimxy  = dimx  * dimy;

  /* renormalize atoms distribution for pgm files*/
  fmax = 0.0;
  for (i=0; i<ntypes*dimxy; i++) fmax = MAX(fmax,atoms_dist[i]);
  fmax = 255/fmax;
  for (i=0; i<ntypes*dimxy; i++) pix[i] = (char) (255-atoms_dist[i]*fmax);

  /* write pgm files */
  for (t=0; t<ntypes; t++) {
    sprintf(fname,"%s.%d.pgm",infile,t);
    if (NULL==(out=fopen(fname,"w"))) error("Cannot open pgm-file.");
    fprintf(out,"P5\n%d %d\n255\n", dimx, dimy);
    if (dimxy!=fwrite(pix+t*dimxy,sizeof(char),dimxy,out))
      error("Cannot write pgm-file.");
    fclose(out);
  }

  /* write ppm-file */
  if (ntypes<4) {
    /* set everything to white */
    for (i=0; i<3*dimxy; i++) pix2[i] = 255;
    /* get renormalization factor for color pictures */
    fmax = 0.0;
    for (i=0; i<dimxy; i++) {
      float tmp=0.0;
      for (t=0; t<ntypes; t++) tmp += atoms_dist[t*dimxy+i];
      fmax = MAX(fmax,tmp);
    }
    fmax = 255/fmax;
    /* make color picture */
    for (t=0; t<ntypes; t++) 
      for (i=0; i<dimxy; i++) 
        for (k=0; k<3; k++) 
          if (k!=t) pix2[3*i+k] -= (char) (atoms_dist[t*dimxy+i]*fmax);
    /* write color picture */
    sprintf(fname,"%s.ppm",infile);
    if (NULL==(out=fopen(fname,"w"))) error("Cannot open pgm-file.");
    fprintf(out,"P6\n%d %d\n255\n", dimx, dimy);
    if (3*dimxy!=fwrite(pix2,sizeof(char),3*dimxy,out))
      error("Cannot write pgm-file.");
    fclose(out);
  }
}

/******************************************************************************
*
*  main
*
******************************************************************************/

int main(int argc, char **argv) 
{
  int dim, dir, min, max;

  /* check number of arguments */
  if ((argc!=2) && (argc!=5)) {
     printf("Usage:  %s infile [dir min max]\n",argv[0]);
     exit(1);
  }

  /* scan further arguments */
  if (argc==5) {
    sscanf(argv[2], "%d", &dir);
    sscanf(argv[3], "%d", &min); 
    sscanf(argv[4], "%d", &max);
  }

  /* read atom distribution */
  dim = read_atoms_dist(argv[1]);

  /* write summary or pictures */
  if ((argc==2) && (dim==3)) {
    axis_projections_3d();
  } else  if ((argc==2) && (dim==2)) {
      pictures_2d(argv[1]);
  } else {
    if (dir==3) {
      xy_pictures_3d(argv[1],min,max);
    } else if (dir==2) { 
      xz_pictures_3d(argv[1],min,max);
    } else if (dir==1) {
      yz_pictures_3d(argv[1],min,max);
    }
  }
  return 0;
}
