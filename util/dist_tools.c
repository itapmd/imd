
#include "dist_tools.h"

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
*  do endian swap for a two byte quantity
*
******************************************************************************/

void swap_2_bytes(void *data)
{
  char c, *str;
  str = (char *) data;
  c = str[0]; str[0] = str[1]; str[1] = c;
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
  if (endian()==1) {
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
  if (endian()==1) {
    str[0] = c[0]; str[1] = c[1]; str[2] = c[2]; str[3] = c[3];
  } 
  else {
    str[0] = c[3]; str[1] = c[2]; str[2] = c[1]; str[3] = c[0];
  }
}

/******************************************************************************
*
*  rgb color of float value
*
******************************************************************************/

void rgbcolor( unsigned char *dat, float val)
{
  int i;
  float v;
  float rtab[5] = { 0.02, 0.03, 0.02, 0.23, 0.45 };
  float gtab[5] = { 0.02, 0.23, 0.45, 0.23, 0.02 };
  float btab[5] = { 0.45, 0.23, 0.02, 0.03, 0.02 };

  /*  Original table from RUS - we use more saturation
  float rtab[5] = { 0.10, 0.05, 0.10, 0.75, 0.75 };
  float gtab[5] = { 0.20, 0.75, 0.50, 0.75, 0.05 };
  float btab[5] = { 0.50, 0.75, 0.25, 0.05, 0.05 }; */
  
  /* Get RBG values from linear interpolation */
  i = (int) (val * 3.9999);
  v = 4.0 * val - i;
  dat[0] = (unsigned char) (255 * (rtab[i] + v * (rtab[i+1] - rtab[i])));
  dat[1] = (unsigned char) (255 * (gtab[i] + v * (gtab[i+1] - gtab[i])));
  dat[2] = (unsigned char) (255 * (btab[i] + v * (btab[i+1] - btab[i])));
}

/******************************************************************************
*
*  initialize float distribution
*
******************************************************************************/

void init_float_dist( float_dist_t *fl )
{
  fl->max_n=0;
  fl->max_len=0;
  fl->min=NULL;
  fl->max=NULL;
  fl->dat=NULL;
  fl->cont=NULL;
}

/******************************************************************************
*
*  allocate float distribution
*
******************************************************************************/

void alloc_float_dist( float_dist_t *fl )
{
  /* allocate or enlarge arrays */
  if (fl->n * fl->nbin > fl->max_len) {
    fl->max_len = fl->n * fl->nbin;
    fl->dat = (float *) realloc(fl->dat, fl->max_len * sizeof(float) );
  }
  if (fl->n > fl->max_n) {
    fl->max_n = fl->n;
    fl->min  = (float *) realloc( fl->min,  fl->max_n  * sizeof(float) );
    fl->max  = (float *) realloc( fl->max,  fl->max_n  * sizeof(float) );
    fl->cont = (char **) realloc( fl->cont, fl->max_n  * sizeof(char*) );
  }
  if ((NULL==fl->dat) || (NULL==fl->min) || 
      (NULL==fl->max) || (NULL==fl->cont))
    error("Cannot allocate arrays");
}

/******************************************************************************
*
*  initialize byte distribution
*
******************************************************************************/

void init_byte_dist( byte_dist_t *bt )
{
  bt->max_len=0;
  bt->dat=NULL;
}

/******************************************************************************
*
*  allocate byte distribution
*
******************************************************************************/

void alloc_byte_dist( byte_dist_t *bt )
{
  /* allocate or enlarge array */
  if (bt->n * bt->nbin > bt->max_len) {
    bt->max_len = bt->n * bt->nbin;
    bt->dat = (unsigned char *) 
      realloc(bt->dat, bt->max_len * sizeof(unsigned char) );
    if (NULL==bt->dat) error("Cannot allocate byte array");
  }
}

/******************************************************************************
*
*  read float distribution
*
******************************************************************************/

int read_float_dist( float_dist_t *dist, char *fname )
{
  FILE *infile;
  char line[255], str[255], *contents, *token;
  int  i, j, k, n, ix, iy, iz;
  int  cont=1, ascii=0, input_endian, n_coord;

  /* open file */
  if (NULL==(infile=fopen(fname,"r"))) return -1;

  /* read file header */
  do {
    fgets(line,255,infile);
    if (line[0]!='#') {
      error("file header corrupt!");
    } else {
      if (line[1]=='F') {
        n = sscanf(line+2, "%s %d %d %d", 
                   str, &(dist->dim), &n_coord, &(dist->n) );
        if ((n_coord!=0) && (n_coord!=dist->dim))
          error("strange number of bin coordinates!");
        if (n<4) error("file header corrupt (format)!");
        if      (str[0]=='B') input_endian = 1;
        else if (str[0]=='L') input_endian = 0;
        else if (str[0]=='A')        ascii = 1;
        else error("file header corrupt (format)!"); 
      } else  if (line[1]=='C') {
        contents = strdup( line );
      } else  if (line[1]=='D') {
        if (dist->dim != sscanf(line+2, "%d %d %d",
                                &(dist->dimx), &(dist->dimy), &(dist->dimz)))
          error("file header corrupt (dimension)!");
        if (dist->dim<3) dist->dimz = 1;
        if (dist->dim<2) dist->dimy = 1;
      } else  if (line[1]=='S') {
        if (dist->dim != sscanf(line+2, "%f %f %f",
                                &(dist->sx), &(dist->sy), &(dist->sz)))
          error("file header corrupt (pixel/voxel size)!");
        if (dist->dim<3) dist->dimz = 1.0;
        if (dist->dim<2) dist->dimy = 1.0;
      } else  if (line[1]=='E') {
	cont=0;
      }
    }
  } while (cont);
  
  /* compute array size */
  dist->nbin = dist->dimx * dist->dimy * dist->dimz;
  dist->len =  dist->nbin * dist->n;

  /* (re)allocate distribution */
  alloc_float_dist( dist );

  /* set contents entries */
  token = strtok(contents," \t\r\n");
  for (i=0; i<n_coord; i++) token = strtok(NULL," \t\r\n");
  for (i=0; i<dist->n; i++) dist->cont[i] = strdup( strtok(NULL," \t\r\n") );

  /* read distribution */
  if (ascii) {
    if (n_coord==2) { 
      for (i=0; i<dist->nbin; i++) {
        fscanf(infile, "%d %d", &ix, &iy);
        j = IDX(dist,ix,iy,0,0);
        for (k=0; k<dist->n; k++) fscanf(infile, "%f", dist->dat + j++);
      }
    } else if (n_coord==3) {
      for (i=0; i<dist->nbin; i++) {
        fscanf(infile, "%d %d %d", &ix, &iy, &iz);
        j = IDX(dist,ix,iy,iz,0);
        for (k=0; k<dist->n; k++) fscanf(infile, "%f", dist->dat + j++);
      }
    } else {
      for (i=0; i<dist->len; i++) {
        fscanf(infile, "%f", dist->dat + i);
      }
    }
  } else {
    if (dist->len != fread(dist->dat, sizeof(float), dist->len, infile))
      error("Cannot read distribution");
    if (input_endian != endian()) {
      for (i=0; i<dist->len; i++) swap_4_bytes(dist->dat + i);
    }
  }

  /* close input file */
  fclose(infile);

  /* compute minima and maxima */
  for (i=0; i<dist->n; i++) {
    dist->min[i] = dist->dat[i];
    dist->max[i] = dist->dat[i];
  }
  for (i=0; i<dist->nbin; i++) {
    for (j=0; j<dist->n; j++) {
      k = i * dist->n + j;
      dist->min[j] = MIN( dist->min[j], dist->dat[k] );
      dist->max[j] = MAX( dist->max[j], dist->dat[k] );
    }
  }
  return 0;
}

/******************************************************************************
*
*  write float distribution (binary only)
*
******************************************************************************/

void write_float_dist( float_dist_t *dist, char *fname )
{
  FILE *outfile;
  char c;
  int  i;

  /* open file */
  if (NULL==(outfile=fopen(fname,"w"))) error("Cannot write distribution");

  /* write file header */
  if (endian()) c='B'; else c='L';
  fprintf(outfile, "#F %c %d 0 %d\n", c, dist->dim, dist->n);
  fprintf(outfile, "#C");
  for (i=0; i<dist->n; i++) fprintf(outfile, " %s", dist->cont[i]);
  fprintf(outfile, "\n");
  if (2==dist->dim) {
    fprintf(outfile, "#D %d %d\n",    dist->dimx, dist->dimy);
    fprintf(outfile, "#S %e %e\n",    dist->sx,   dist->sy  );
  }
  else if (3==dist->dim) {
    fprintf(outfile, "#D %d %d %d\n", dist->dimx, dist->dimy, dist->dimz);
    fprintf(outfile, "#S %e %e %e\n", dist->sx,   dist->sy,   dist->sz  );
  }
  else error("Strange dimension of distribution!");
  fprintf(outfile, "#E\n");

  /* write data */
  if (dist->len!=fwrite(dist->dat, sizeof(float), dist->len, outfile))
    error("Cannot write distribution");
  fclose(outfile);
}

/******************************************************************************
*
*  cut and convert float to scalar8 distribution
*
******************************************************************************/

void float2scalar8_dist( float_dist_t *fl, byte_dist_t *bt, float min, 
  float max, int n, int llx, int lly, int llz, int urx, int ury, int urz)
{
  int ix, iy, iz, j, k;
  float tmp, sc;

  if ((llx<0) || (urx>fl->dimx) || 
      (lly<0) || (ury>fl->dimy) || 
      (llz<0) || (urz>fl->dimz)) error("cut region not contained in dist");
  if (n >= fl->n) error("too few components");

  bt->dim  = fl->dim;
  bt->dimx = urx - llx;
  bt->dimy = ury - lly;
  bt->dimz = urz - llz;
  bt->len  = bt->dimx * bt->dimy * bt->dimz;
  bt->nbin = bt->len;
  bt->sx   = fl->sx;
  bt->sy   = fl->sy;
  bt->sz   = fl->sz;
  bt->n    = 1;

  /* (re)allocate distribution */
  alloc_byte_dist( bt );

  /* most volume renderers need the x-direction running the fastest */
  sc = 255.99 / (max - min);
  k  = 0;
  for (iz=llz; iz<urz; iz++)
    for (iy=lly; iy<ury; iy++)
      for (ix=llx; ix<urx; ix++) {
        tmp = MAX( min, ELM(fl,ix,iy,iz,n) );
        tmp = MIN( tmp, max ) * sc;
        bt->dat[k++] = (unsigned char) tmp;
      }
}

/******************************************************************************
*
*  cut and convert float to scalar12 distribution
*
******************************************************************************/

void float2scalar12_dist( float_dist_t *fl, byte_dist_t *bt, float min, 
  float max, int n, int llx, int lly, int llz, int urx, int ury, int urz)
{
  int ix, iy, iz, i, j, k;
  unsigned short *us;
  float tmp, sc;

  if ((llx<0) || (urx>fl->dimx) || 
      (lly<0) || (ury>fl->dimy) || 
      (llz<0) || (urz>fl->dimz)) error("cut region not contained in dist");
  if (n >= fl->n) error("too few components");

  bt->dim  = fl->dim;
  bt->dimx = urx - llx;
  bt->dimy = ury - lly;
  bt->dimz = urz - llz;
  bt->nbin = bt->dimx * bt->dimy * bt->dimz;
  bt->n    = 2; /* we have two bytes per voxel */
  bt->len  = bt->nbin * bt->n;
  bt->sx   = fl->sx;
  bt->sy   = fl->sy;
  bt->sz   = fl->sz;

  /* (re)allocate distribution */
  alloc_byte_dist( bt );

  /* most volume renderers need the x-direction running the fastest */
  sc = 4095.99 / (max - min);
  us = (unsigned short *) bt->dat;
  k  = 0;
  for (iz=llz; iz<urz; iz++)
    for (iy=lly; iy<ury; iy++)
      for (ix=llx; ix<urx; ix++) {
        tmp = MAX( min, ELM(fl,ix,iy,iz,n) );
        tmp = MIN( tmp, max ) * sc;
        us[k++] = (unsigned short) tmp;
      }

  /* endian swap if necessary - we need little endian */
  if (endian()) 
    for (i=0; i<k; i++) swap_2_bytes( (void *) (us+i) );

}

/******************************************************************************
*
*  cut and convert 2D float to ppm distribution
*
******************************************************************************/

void float2ppm_dist( float_dist_t *fl, byte_dist_t *bt, float min, float max,
  int n, int llx, int lly, int urx, int ury, 
  void (*color)(unsigned char*, float))
{
  int ix, iy, iz, j, k;
  float tmp, sc;

  if (fl->dim!=2) error("Need 2D distribution!");
  if ((llx<0) || (urx>fl->dimx) || 
      (lly<0) || (ury>fl->dimy)) error("cut region not contained in dist");
  if (n >= fl->n) error("too few components");

  bt->dim  = fl->dim;
  bt->dimx = urx - llx;
  bt->dimy = ury - lly;
  bt->nbin = bt->dimx * bt->dimy;
  bt->len  = bt->nbin * 3;
  bt->sx   = fl->sx;
  bt->sy   = fl->sy;
  bt->n    = 3;

  /* (re)allocate distribution */
  alloc_byte_dist( bt );

  sc = 1.0 / (max - min);
  k  = 0;
  for (iy=ury-1; iy>=lly; iy--)
    for (ix=llx; ix<urx; ix++) {
      tmp = MAX( min, ELM(fl,ix,iy,0,n) );
      tmp = sc * (MIN( tmp, max ) - min);
      color(bt->dat + k, tmp); 
      k+=3;
    }
}

/******************************************************************************
*
*  write 2D ppm distribution to ppm file
*
******************************************************************************/

void write_ppm( byte_dist_t *ppm, char *fname )
{
  FILE *outfile;
  if (ppm->dim!=2) error("2D distribution required");
  if (NULL==(outfile=fopen(fname,"w"))) error("Cannot open ppm file");
  fprintf( outfile, "P6\n%d %d\n255\n", ppm->dimx, ppm->dimy );
  fwrite( ppm->dat, sizeof(unsigned char), ppm->len, outfile ); 
  fclose( outfile);
}
