
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
*  diffpat is a utility program to analyse/convert 3D diffraction patterns.
*
*  Compilation: 
*
*    make diffpat
*
*  Usage:
*
*    diffpat -i file 
*
*  projects the diffraction pattern onto the three axis, and writes these 
*  1D histograms to the terminal. This helps to find out where the dense
*  planes are.
*
*    diffpat [-v] [-c fmax] [-l fmin] file dir min max
*
*  adds the slices min..max-1 perpendicular to dir, and writes 
*  a pgm file, or a virvo file if option -v is given. If fmax is given, 
*  the intensity is cut at fmax. By default, the intensities are mapped 
*  linearly to the gray values. If the option -l is given, the intensities 
*  are mapped logarithmically to the gray values, and intensities smaller 
*  than fmin are set to zero. 
*
*    diffpat [-c fmax] [-l fmin] file xmin ymin zmin xmax ymax zmax
*
*  cuts a rectangular block from the volume and writes it in virvo xvf 
*  format. If fmax is given, the intensity is cut at fmax. By default, 
*  the intensities are mapped linearly to the gray values. If the option 
*  -l is given, the intensities are mapped logarithmically to the gray 
*  values, and intensities smaller than fmin are set to zero. 
*
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

float *diffpat, ddx, ddy, ddz;
int   dim, dimx, dimy, dimz, nx, ny, nz, size, my_endian;

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

void swap_4_bytes(void *data)
{
  char c, *str;
  str = (char *) data;
  c = str[0]; str[0] = str[3]; str[3] = c;
  c = str[1]; str[1] = str[2]; str[3] = c;
}

/******************************************************************************
*
*  copy 2 bytes to big endian
*
******************************************************************************/

void copy_2_bytes(char *str, void *source)
{
  char *c;
  c = source;
  if (my_endian==1) {
    str[0] = c[0]; str[1] = c[1];
  } 
  else {
    str[0] = c[1]; str[1] = c[0]   ;
  }
}

/******************************************************************************
*
*  copy 4 bytes to big endian
*
******************************************************************************/

void copy_4_bytes(char *str, void *source)
{
  char *c;
  c = source;
  if (my_endian==1) {
    str[0] = c[0]; str[1] = c[1]; str[2] = c[2]; str[3] = c[3];
  } 
  else {
    str[0] = c[3]; str[1] = c[2]; str[2] = c[1]; str[3] = c[0];
  }
}

/******************************************************************************
*
*  read diffraction pattern
*
******************************************************************************/

int read_diffpat(char *fname, float fmax, float dyn)
{
  FILE  *infile;
  float *dist, x, y, f000;
  char  line[255];
  char  str[255];
  int   cont, input_endian, i, j, k, sz, nnx, nny, nnz, ind1, ind2, n;

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
        n = sscanf(line+2,"%s %d",str,&dim);
        if      (str[0]=='B') input_endian=1;
        else if (str[0]=='L') input_endian=0;
        else error("file header corrupt!"); 
      } else  if (line[1]=='D') {
        if (dim!=sscanf(line+2,"%d %d %d",&nnx,&nny,&nnz))
          error("file header corrupt (dimension)!");
      } else  if (line[1]=='S') {
        if (dim!=sscanf(line+2,"%f %f %f",&ddx,&ddy,&ddz))
          error("file header corrupt (dimension)!");
      } else  if (line[1]=='E') {
	cont=0;
      }
    }
  } while (cont);
  
  /* read diffraction pattern */
  sz   = nnx * nny * nnz;
  dist = (float *) malloc( sz * sizeof(float) );
  if (NULL==dist) error("Cannot allocate data for diffraction pattern");
  if (sz!=fread(dist,sizeof(float),sz,infile))
    error("Cannot read diffraction pattern data");
  fclose(infile);

  /* do endian swap if necessary */
  if (input_endian!=endian()) {
    for (i=0; i<sz; i++) swap_4_bytes(dist+i);
  }

  /* covert to more convenient data layout */
  nx = nnx / 2; dimx = 2 * nx + 1;
  ny = nny / 2; dimy = 2 * ny + 1;
  nz = nnz - 1; dimz = 2 * nz + 1;
  size = dimx * dimy * dimz;

  /* allocate new diffraction pattern */
  diffpat = (float *) calloc( size, sizeof(float) );
  if (NULL==diffpat) error("Cannot allocate data for diffraction pattern");

  /* if pseudo-3d (single slice), do just a transposition */
  if (nz==0) {
    for (i=0; i<nnx; i++)
      for (j=0; j<nny; j++)
        diffpat[i * nny + j] = dist[j * nnx + i];
  } 
  else {
    /* +++ quadrant */
    for (i=0; i<=nx; i++)
      for (j=0; j<=ny; j++)
        for (k=0; k<=nz; k++) {
          ind1 = ((i + nx) * dimy + (j + ny)) * dimz + (k + nz);
          ind2 = (i * nny + j) * nnz + k;
          diffpat[ind1] = dist[ind2];
        }

    /* -++ quadrant */
    for (i=1; i<=nx; i++)
      for (j=0; j<=ny; j++)
        for (k=0; k<=nz; k++) {
          ind1 = ((-i + nx) * dimy + (j + ny)) * dimz + (k + nz);
          ind2 = ((nnx - i) * nny + j) * nnz + k;
          diffpat[ind1] = dist[ind2];
        }

    /* +-+ quadrant */
    for (i=0; i<=nx; i++)
      for (j=1; j<=ny; j++)
        for (k=0; k<=nz; k++) {
          ind1 = ((i + nx) * dimy + (-j + ny)) * dimz + (k + nz);
          ind2 = (i * nny + (nny - j)) * nnz + k;
          diffpat[ind1] = dist[ind2];
        }

    /* --+ quadrant */
    for (i=1; i<=nx; i++)
      for (j=1; j<=ny; j++)
        for (k=0; k<=nz; k++) {
          ind1 = ((-i + nx) * dimy + (-j + ny)) * dimz + (k + nz);
          ind2 = ((nnx - i) * nny + (nny - j)) * nnz + k;
          diffpat[ind1] = dist[ind2];
        }

    /* lower half */
    for (i=-nx; i<=nx; i++)
      for (j=-ny; j<=ny; j++)
        for (k=1; k<=nz; k++) {
          ind1 = ((-i + nx) * dimy + (-j + ny)) * dimz + (-k + nz);
          ind2 = (( i + nx) * dimy + ( j + ny)) * dimz + ( k + nz);
          diffpat[ind1] = diffpat[ind2];
        }

    /* free temporary memory */
    free(dist);
  }
  printf("f000: %e\n", diffpat[(nx * dimy + ny) * dimz + nz] );

  /* cut at fmax */
  if (fmax>0) {
    for (i=0; i<size; i++) {
      diffpat[i] = MIN(diffpat[i],fmax);
    }
  }

  /* logaritmic scale */
  if (dyn>0) {
    float x, y;
    /*
    y = (float) log( (double) (diffpat[(nx * dimy + ny) * dimz + nz] / dyn) );
    */
    y = (float) log( (double) dyn );
    for (i=0; i<size; i++) {
      x = (float) (log( (double) diffpat[i] ) - y); 
      diffpat[i] = MAX(x,0.0);
    }
  }

  return dim;

}

/******************************************************************************
*
*  project diffraction pattern on the three main axis
*
******************************************************************************/

void axis_projections_3d()
{
  float *histx, *histy, *histz;
  int   t, i, j, k, l;

  histx = (float *) calloc( dimx, sizeof(float) );
  histy = (float *) calloc( dimy, sizeof(float) );
  histz = (float *) calloc( dimz, sizeof(float) );
  if ((histx==NULL) || (histy==NULL) || (histx==NULL))
    error("out of memory");

  for (i=0; i<dimx; i++)
    for (j=0; j<dimy; j++)
      for (k=0; k<dimz; k++) {
        l = (i * dimy + j) * dimz + k;
        histx[i] += diffpat[l];
        histy[j] += diffpat[l];
        histz[k] += diffpat[l];
      }

  printf("# Projection on x-Axis:\n");
  for (i=0; i<dimx; i++) printf("%d %e\n", i-nx, histx[i]);

  printf("\n\n# Projection on y-Axis:\n");
  for (i=0; i<dimy; i++) printf("%d %e\n", i-ny, histy[i]);

  printf("\n\n# Projection on z-Axis:\n");
  for (i=0; i<dimz; i++) printf("%d %e\n", i-nz, histz[i]);

}

/******************************************************************************
*
*  write pgm file
*
******************************************************************************/

void write_pgm(char *outfile, char *pix, int dx, int dy)
{
  int  len = dx * dy;
  char fname[255];
  FILE *out;

  sprintf(fname,"%s.pgm",outfile);
  if (NULL==(out=fopen(fname,"w"))) 
    error("Cannot open pgm-file.");
  fprintf(out,"P5\n%d %d\n255\n", dx, dy);
  if (len!=fwrite(pix,sizeof(char),len,out))
    error("Cannot write pgm-file.");
  fclose(out);

}

/******************************************************************************
*
*  write 2D virvo file
*
******************************************************************************/

void write_virvo_2d(char *outfile, char *pix, int dx, int dy, 
                    float sx, float sy, float sz)
{
  unsigned short us;
  unsigned int   ui;
  unsigned char  uc;
  float fl;
  int   len = dx * dy;
  char  header[48], *str, fname[255]; 
  FILE  *out;

  /* make volume file header */
  my_endian=endian();
  str = header;
  sprintf(str,"%s","VIRVO-XVF");     str +=9;   /* file type */
  us =  48; copy_2_bytes(str,&us);   str +=2;   /* header size */
  ui =  dx; copy_4_bytes(str,&ui);   str +=4;   /* dim_x */
  ui =  dy; copy_4_bytes(str,&ui);   str +=4;   /* dim_y */
  ui =   1; copy_4_bytes(str,&ui);   str +=4;   /* dim_z */
  ui =   1; copy_4_bytes(str,&ui);   str +=4;   /* number of frames */
  uc =   8; *str = uc;               str++;     /* bits per voxel */
  fl =  sx; copy_4_bytes(str,&fl);   str +=4;   /* x-length of voxel */
  fl =  sy; copy_4_bytes(str,&fl);   str +=4;   /* y-length of voxel */
  fl =  sz; copy_4_bytes(str,&fl);   str +=4;   /* z-length of voxel */
  fl = 1.0; copy_4_bytes(str,&fl);   str +=4;   /* secs per frame */
  us =   0; copy_2_bytes(str,&us);   str +=2;   /* number of transf. func. */
  us =   0; copy_2_bytes(str,&us);   str +=2;   /* type of transf. func. */

  /* write virvo file */
  sprintf(fname,"%s.xvf",outfile);
  if (NULL==(out=fopen(fname,"w"))) 
    error("Cannot open volume file.");
  if (48!=fwrite(header,sizeof(char),48,out))
    error("Cannot write volume file header.");
  if (len!=fwrite(pix,sizeof(char),len,out))
    error("Cannot write volume file data.");
  fclose(out);

}

void write_virvo_2d_16bit(char *outfile, unsigned short *pix, int dx, int dy, 
                    float sx, float sy, float sz)
{
  unsigned short us;
  unsigned int   ui;
  unsigned char  uc;
  float fl;
  int   len = dx * dy * 2;
  char  header[48], *str, fname[255]; 
  FILE  *out;

  /* make volume file header */
  my_endian=endian();
  str = header;
  sprintf(str,"%s","VIRVO-XVF");     str +=9;   /* file type */
  us =  48; copy_2_bytes(str,&us);   str +=2;   /* header size */
  ui =  dx; copy_4_bytes(str,&ui);   str +=4;   /* dim_x */
  ui =  dy; copy_4_bytes(str,&ui);   str +=4;   /* dim_y */
  ui =   1; copy_4_bytes(str,&ui);   str +=4;   /* dim_z */
  ui =   1; copy_4_bytes(str,&ui);   str +=4;   /* number of frames */
  uc =  16; *str = uc;               str++;     /* bits per voxel */
  fl =  sx; copy_4_bytes(str,&fl);   str +=4;   /* x-length of voxel */
  fl =  sy; copy_4_bytes(str,&fl);   str +=4;   /* y-length of voxel */
  fl =  sz; copy_4_bytes(str,&fl);   str +=4;   /* z-length of voxel */
  fl = 1.0; copy_4_bytes(str,&fl);   str +=4;   /* secs per frame */
  us =   0; copy_2_bytes(str,&us);   str +=2;   /* number of transf. func. */
  us =   0; copy_2_bytes(str,&us);   str +=2;   /* type of transf. func. */

  /* write virvo file */
  sprintf(fname,"%s.xvf",outfile);
  if (NULL==(out=fopen(fname,"w"))) 
    error("Cannot open volume file.");
  if (48!=fwrite(header,sizeof(char),48,out))
    error("Cannot write volume file header.");
  if (len!=fwrite(pix,sizeof(char),len,out))
    error("Cannot write volume file data.");
  fclose(out);

}


/******************************************************************************
*
*  write pseudo-2D diffraction pattern
*
******************************************************************************/

void write_2d(char *fname, float *dist, int dx, int dy, 
              float sx, float sy, float sz)
{
  char c;
  int  n = dx * dy;
  FILE *out;

  if (NULL==(out=fopen(fname,"w"))) 
    error("Cannot open output file.");
  if (endian()) c='B'; 
  else          c='L';
  fprintf(out,"#F %c 3 0 1\n",c);
  fprintf(out,"#D %d %d 1\n", dx, dy);
  fprintf(out,"#S %f %f %f\n", sx, sy, sz);
  fprintf(out,"#E\n");
  if (n!=fwrite(dist,sizeof(float),n,out))
    error("Cannot write output file");
  fclose(out);
}


/******************************************************************************
*
*  write pictures of xy-slices
*
******************************************************************************/

void dist_xy(char *infile, int min, int max, int mode)
{
  float *hist, inv, fmax;
  char  *pix;
  int   dimxy, i, j, k, l;
  FILE  *out;
  char  fname[255]; 

  /* allocate data */
  dimxy  = dimx  * dimy;
  hist = (float *) calloc( dimxy, sizeof(float) );
  pix  = (char  *) calloc( dimxy, sizeof(char ) );
  if ((hist==NULL) || (pix==NULL)) error("out of memory");

  /* make distribution */
  for (i=0; i<dimx; i++)
    for (j=0; j<dimy; j++)
      for (k=min; k<max; k++) {
        l = (i * dimy + j) * dimz + k + nz;
        hist[j * dimx + i] += diffpat[l];
      }

  /* normalize distribution */
  if (mode<2) {
    fmax = 0.0;
    for (i=0; i<dimxy; i++) fmax = MAX(fmax,hist[i]);
    inv = 255/fmax;
    if (mode==1)
      for (i=0; i<dimxy; i++) pix[i] = (char) (      MIN(hist[i],fmax) * inv);
    else
      for (i=0; i<dimxy; i++) pix[i] = (char) (255 - MIN(hist[i],fmax) * inv);
  }

  /* make file name */
  if (nz==0) {
    sprintf(fname,"%s",infile);
  }
  else if (max-min==1) {
    sprintf(fname,"%s.xy.%d",infile,min);
  }
  else {
    sprintf(fname,"%s.xy.%d-%d",infile,min,max-1);
  }

  /* write distribution */
  if (mode==1) {
    write_virvo_2d(fname, pix, dimx, dimy, ddx, ddy, ddz);
  }
  else if (mode==2) {
    write_2d(fname, hist, dimx, dimy, ddx, ddy, ddz);
  } 
  else {
    write_pgm(fname, pix, dimx, dimy);
  }
}

/******************************************************************************
*
*  write pictures of xz-slices
*
******************************************************************************/

void dist_xz(char *infile, int min, int max, int mode)
{
  float *hist, inv, fmax;
  char  *pix;
  int   dimxz, i, j, k, l;
  FILE  *out;
  char  fname[255]; 

  /* allocate data */
  dimxz  = dimx  * dimz;
  hist = (float *) calloc( dimxz, sizeof(float) );
  pix  = (char  *) calloc( dimxz, sizeof(char ) );
  if ((hist==NULL) || (pix==NULL)) error("out of memory");

  /* make distribution */
  for (i=0; i<dimx; i++)
    for (j=min; j<max; j++)
      for (k=0; k<dimz; k++) {
        l = (i * dimy + j + ny) * dimz + k;
        hist[k * dimx + i] += diffpat[l];
      }

  /* normalize distribution */
  if (mode<2) {
    fmax = 0.0;
    for (i=0; i<dimxz; i++) fmax = MAX(fmax,hist[i]);
    inv = 255/fmax;
    if (mode==1)
      for (i=0; i<dimxz; i++) pix[i] = (char) (      MIN(hist[i],fmax) * inv);
    else
      for (i=0; i<dimxz; i++) pix[i] = (char) (255 - MIN(hist[i],fmax) * inv);
  }

  /* make file name */
  if (nz==0) {
    sprintf(fname,"%s",infile);
  }
  else if (max-min==1) {
    sprintf(fname,"%s.xz.%d",infile,min);
  }
  else {
    sprintf(fname,"%s.xz.%d-%d",infile,min,max-1);
  }

  /* write distribution */
  if (mode==1) {
    write_virvo_2d(fname, pix, dimx, dimz, ddx, ddz, ddy);
  }
  else if (mode==2) {
    write_2d(fname, hist, dimx, dimz, ddx, ddz, ddy);
  }
  else {
    write_pgm(fname, pix, dimx, dimz);
  }
}

/******************************************************************************
*
*  write pictures of yz-slices
*
******************************************************************************/

void dist_yz(char *infile, int min, int max, int mode)
{
  float *hist, inv, fmax;
  char  *pix;
  int   dimyz, i, j, k, l;
  FILE  *out;
  char  fname[255]; 

  /* allocate data */
  dimyz = dimy  * dimz;
  hist  = (float *) calloc( dimyz, sizeof(float) );
  pix   = (char  *) calloc( dimyz, sizeof(char ) );
  if ((hist==NULL) || (pix==NULL)) error("out of memory");

  /* make distribution */
  for (i=min; i<max; i++)
    for (j=0; j<dimy; j++)
      for (k=0; k<dimz; k++) {
        l = ((i + nx) * dimy + j) * dimz + k;
        hist[k * dimy + j] += diffpat[l];
      }

  /* normalize distribution */
  if (mode<2) {
    fmax = 0.0;
    for (i=0; i<dimyz; i++) fmax = MAX(fmax,hist[i]);
    inv = 255/fmax;
    if (mode==1)
      for (i=0; i<dimyz; i++) pix[i] = (char) (      MIN(hist[i],fmax) * inv);
    else
      for (i=0; i<dimyz; i++) pix[i] = (char) (255 - MIN(hist[i],fmax) * inv);
  }

  /* make file name */
  if (nz==0) {
    sprintf(fname,"%s",infile);
  }
  else if (max-min==1) {
    sprintf(fname,"%s.yz.%d",infile,min);
  }
  else {
    sprintf(fname,"%s.yz.%d-%d",infile,min,max-1);
  }

  /* write distribution */
  if (mode==1) {
    write_virvo_2d(fname, pix, dimy, dimz, ddy, ddz, ddx);
  }
  else if (mode==2) {
    write_2d(fname, hist, dimy, dimz, ddy, ddz, ddx);
  }
  else {
    write_pgm(fname, pix, dimy, dimz);
  }
}

/******************************************************************************
*
*  3D virvo volume data
*
******************************************************************************/

void virvo_picture_3d(char *infile, int min_x, int min_y, int min_z,
                                    int max_x, int max_y, int max_z)
{
  float max=0.0, fl;
  unsigned char *vol;
  int   len, dx, dy, dz, t, i, j, k, l, m;
  unsigned short us;
  unsigned int   ui;
  unsigned char  uc;
  char  fname[255], header[48], *str; 
  FILE  *out;

  min_x = MAX(min_x,0); max_x = MIN(max_x,dimx); dx = max_x - min_x;
  min_y = MAX(min_y,0); max_y = MIN(max_y,dimy); dy = max_y - min_y;
  min_z = MAX(min_z,0); max_z = MIN(max_z,dimz); dz = max_z - min_z;
  len   = dx * dy * dz;

  /* make volume file header */
  my_endian=endian();
  str = header;
  sprintf(str,"%s","VIRVO-XVF");    str +=9;   /* file type */
  us = 48; copy_2_bytes(str,&us);   str +=2;   /* header size */
  ui = dx; copy_4_bytes(str,&ui);   str +=4;   /* dim_x */
  ui = dy; copy_4_bytes(str,&ui);   str +=4;   /* dim_y */
  ui = dz; copy_4_bytes(str,&ui);   str +=4;   /* dim_z */
  ui =  1; copy_4_bytes(str,&ui);   str +=4;   /* number of frames */
  uc =  8; *str = uc;               str++;     /* bits per voxel */
  fl = ddx; copy_4_bytes(str,&fl);  str +=4;   /* x-length of voxel */
  fl = ddy; copy_4_bytes(str,&fl);  str +=4;   /* y-length of voxel */
  fl = ddz; copy_4_bytes(str,&fl);  str +=4;   /* z-length of voxel */
  fl = 1.0; copy_4_bytes(str,&fl);  str +=4;   /* secs per frame */
  us =   0; copy_2_bytes(str,&us);  str +=2;   /* number of transf. func. */
  us =   0; copy_2_bytes(str,&us);  str +=2;   /* type of transf. func. */

  vol  = (unsigned char  *) calloc( len, sizeof(char ) );
  if (vol==NULL) error("out of memory");

  /* compute renormalization factor */
  for (i=min_x; i<max_x; i++)
    for (j=min_y; j<max_y; j++)
      for (k=min_z; k<max_z; k++) {
        l = (i * dimy + j) * dimz + k;
        max = MAX(max,diffpat[l]);
      }
  max = 255/max;

  /* compute volume data */
  for (i=min_x; i<max_x; i++)
    for (j=min_y; j<max_y; j++)
      for (k=min_z; k<max_z; k++) {
        l = (i * dimy + j) * dimz + k;
        m = ((k-min_z) * dy + (j-min_y)) * dx + (i-min_x);
        vol[m] = (unsigned char) (diffpat[l]*max);
      }

  /* write volume data */
  sprintf(fname,"%s.xvf",infile);
  if (NULL==(out=fopen(fname,"w"))) 
    error("Cannot open volume file.");
  if (48!=fwrite(header,sizeof(char),48,out))
    error("Cannot write volume file header.");
  if (len!=fwrite(vol,sizeof(char),len,out))
    error("Cannot write volume file data.");
  fclose(out);

}

/******************************************************************************
*
*  usage
*
******************************************************************************/

void usage(char *progname)
{
  printf("Usage:  %s -i infile\n",progname);
  printf("        %s [-v] [-c fmax] [-l dyn] infile dir min max\n",progname);
  printf("        %s [-c fmax] [-l dyn] infile", progname);
  printf(" xmin ymin zmin xmax ymax zmax\n");
  exit(1);
}

/******************************************************************************
*
*  main
*
******************************************************************************/

int main(int argc, char **argv) 
{
  int dim, dir=0, info=0, mode=0, min, max;
  int min_x, min_y, min_z, max_x, max_y, max_z;
  float fmax=0.0, dyn=0;
  char *progname, *infile;

  /* parse command line options */
  progname = strdup(argv[0]);
  while ((argc > 1) && (argv[1][0] =='-')) {
    if (argv[1][1]=='l') {
      dyn   = atof(argv[2]);
      argc -= 2;
      argv += 2;
    }
    else if (argv[1][1]=='c') {
      fmax  = atof(argv[2]);
      argc -= 2;
      argv += 2;
    }
    else if (argv[1][1]=='i') {
      info  = 1;
      argc -= 1;
      argv += 1;
    }
    else if (argv[1][1]=='v') {
      mode  = 1;
      argc -= 1;
      argv += 1;
    }
    else if (argv[1][1]=='r') {
      mode  = 2;
      argc -= 1;
      argv += 1;
    }
    else usage(progname);
  }

  /* scan further arguments */
  if (argc>1) infile = strdup(argv[1]);
  else usage(progname);
 
  if (argc==5) {
    sscanf(argv[2], "%d", &dir);
    sscanf(argv[3], "%d", &min); 
    sscanf(argv[4], "%d", &max);
  }
  else if (argc==8) {
    sscanf(argv[2], "%d", &min_x);
    sscanf(argv[3], "%d", &min_y); 
    sscanf(argv[4], "%d", &min_z);
    sscanf(argv[5], "%d", &max_x);
    sscanf(argv[6], "%d", &max_y); 
    sscanf(argv[7], "%d", &max_z);
  }
  else usage(progname);

  /* read atom distribution */
  if (3 != read_diffpat(infile,fmax,dyn))
    error("Implemented for 3D data only!");

  /* write summary or pictures */
  if (info) {
    axis_projections_3d();
  } else {
    if (dir==3) {
      dist_xy(infile,min,max,mode);
    } else if (dir==2) { 
      dist_xz(infile,min,max,mode);
    } else if (dir==1) {
      dist_yz(infile,min,max,mode);
    } else if (dir==0) {
      virvo_picture_3d(infile,min_x,min_y,min_z,max_x,max_y,max_z);
    }
  }
  return 0;
}
