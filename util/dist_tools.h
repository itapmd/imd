
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX(a,b)        ((a) > (b) ? (a) : (b))
#define MIN(a,b)        ((a) < (b) ? (a) : (b))
#define IDX(d,x,y,z,j)  ((((x)*(d)->dimy+(y))*(d)->dimz+(z))*(d)->n+(j))
#define ELM(d,x,y,z,j)  ((d)->dat[ IDX(d,x,y,z,j) ])


typedef struct {
  int   dim, dimx, dimy, dimz, n, max_n, nbin, len, max_len;
  float sx, sy, sz;
  float *min, *max, *dat;
  char  **cont;
} float_dist_t;

typedef struct {
  int   dim, dimx, dimy, dimz, n, nbin, len, max_len;
  float sx, sy, sz;
  unsigned char *dat;
} byte_dist_t;

/* utility functions */
void error(char *msg);
int  endian(void);
void swap_4_bytes(void*);
void swap_2_bytes(void*);
void copy_2_bytes(char*,void*);
void copy_4_bytes(char*,void*);
void rgbcolor( unsigned char *dat, float val);

/* initialization and allocation */
void  init_float_dist(float_dist_t *fl);
void alloc_float_dist(float_dist_t *bt);
void   init_byte_dist( byte_dist_t *fl);
void  alloc_byte_dist( byte_dist_t *bt);

/* conversion of distributions */
void float2scalar8_dist(float_dist_t *fl, byte_dist_t *bt, float min, 
  float max, int n, int llx, int lly, int llz, int urx, int ury, int urz);
void float2scalar12_dist(float_dist_t *fl, byte_dist_t *bt, float min, 
  float max, int n, int llx, int lly, int llz, int urx, int ury, int urz);
void float2ppm_dist(float_dist_t *fl, byte_dist_t *bt, float min, float max,
  int n, int llx, int lly, int urx, int ury,
  void (*color)(unsigned char*, float));

/* reading and writing distributions */
int  read_float_dist( float_dist_t *dist, char *fname);
void write_ppm( byte_dist_t *rgb, char *fname );
void write_float_dist( float_dist_t *dist, char *fname );
