
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
*  common routines and global variables, prototyps, makros, typedefs etc.
*  for various IMD utility programs
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

/******************************************************************************
*
*  some usful macros 
*
******************************************************************************/

/* number of slots in the histograms */
#ifndef SLOTS
#define SLOTS 1000 
#endif

/* maximal number of neighbors for CONN */
#ifndef MAXNEIGH
#define MAXNEIGH 10
#endif

/* number of slots in the histograms */
#ifdef TWOD
#define DIM 2
#else
#define DIM 3
#endif

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SQR(a) (a)*(a)
#define TRUNC (int)

#if defined(__GNUC__)
#define INLINE inline
#else
#define INLINE
#endif

/* avoid p % q, which is terribly slow */
/* on SGI, inline doesn't really work :-( */
#if defined(alpha)
#pragma inline(MOD)
#elif defined(sgi)
#pragma inline global (MOD)
#endif
INLINE static int MOD(int p, int q)
{
  int stmp=p;
  while (stmp>=q) stmp-=q;
  return stmp;
}

/* scalar product for vectors */
#ifdef TWOD
#define SPROD(a,b) (((a).x * (b).x) + ((a).y * (b).y))
#else
#define SPROD(a,b) (((a).x * (b).x) + ((a).y * (b).y) + ((a).z * (b).z))
#endif

/* Dynamically allocated 5D arrray -- half vector version */
#define PTR_5D_V(var,n,a,b,c,d,dim) (((var) + \
				      ((n)*(dim.i)*(dim.j)*(dim.k)*(dim.l)) + \
				      ((a)*(dim.j)*(dim.k)*(dim.l)) + \
				      ((b)*(dim.k)*(dim.l)) + \
				      ((c)*(dim.l)) + \
				      (d)))

/* Dynamically allocated 4D arrray -- half vector version */
#define PTR_4D_V(var,i,j,k,l,dim) (((var) + \
                                 ((i)*(dim.x)*(dim.y)*(dim.z)) + \
                                 ((j)*(dim.y)*(dim.z)) + \
                                 ((k)*(dim.z)) + \
                                 (l)))

/* Dynamically allocated 3D arrray -- half vector version */
#define PTR_3D_V(var,i,j,k,dim) (((var) + \
                                 ((i)*(dim.y)*(dim.z)) + \
                                 ((j)*(dim.z)) + \
                                 (k)))

/* Dynamically allocated 3D arrray -- full vector version */
#define PTR_3D_VV(var,coord,dim) (((var) + \
                                 ((coord.x)*(dim.y)*(dim.z)) + \
                                 ((coord.y)*(dim.z)) + \
                                 (coord.z)))

/* Dynamically allocated 2D arrray -- half vector version */
#define PTR_2D_V(var,i,j,dim) (((var) + \
                               ((i)*(dim.y)) + \
                                (j)))
                                

/* Dynamically allocated 2D arrray -- full vector version */
#define PTR_2D_VV(var,coord,dim) (((var) + \
                                 ((coord.x)*(dim.y)) + \
                                  (coord.y)))

#ifdef TWOD
#define PTR_VV   PTR_2D_VV
#else
#define PTR_VV   PTR_3D_VV
#endif

/* allocation increment */
#ifdef STRAIN
#define CSTEP   2
#else
#define CSTEP   10
#endif

/* Tolerance values for imd_stress */
#if defined(STRESS) || defined(ELCO)
#define NUM        200
#define TOL        1.0e-6
#define TOL2       1.0e-10
#define TOL_VERT2  1.0e-6
#define TOL_DIST2  2.0e-11
#endif

#ifdef ELCO
#define I(a,b)   [(((a)*neigh_len) + (b))]
#define J(a,b,c) [(((((a)*neigh_len) + (b))*neigh_len) + (c))]
#endif

#ifdef PS
#define ABS(a) ((a)>=0) ? (a) : (-(a))
#define NSTEP 1
#define PI 3.1415927
#define PIN 0.017453293
#define Max -1.0e10
#define Min 1.0e10
#define PT 0.03527777778
#define TOL 0.00001
#endif

/*****************************************************************************
*
*  include files
*
*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


/*****************************************************************************
*
*  typedefs
*
*****************************************************************************/

#ifdef SINGLE
typedef float real;
#else
typedef double real;
#endif

#ifdef TWOD
typedef struct {real x; real y; } vektor;
typedef struct {int  x; int  y; } ivektor;
#else
typedef struct {real x; real y; real z; } vektor;
typedef struct {int  x; int  y;  int z; } ivektor;
#endif

typedef struct {int  x; int  y; } ivektor2d;
typedef struct {int  x; int  y; int z; } ivektor3d;
typedef struct {int  i; int  x; int  y; int z; } ivektor4d;
typedef struct {int n; int  i; int  j; int  k; int l; } ivektor5d;

#ifdef PS
typedef struct {real x; real y; } vektor2d;
typedef struct {real x; real y; real z;} vektor3d;
#endif

#ifdef ELCO
typedef struct { real xx, xy, xz, yx, yy, yz, zx, zy, zz; } tensor;

typedef struct { real c11, c12, c13, c14, c15, c16, c22, c23, c24, 
  c25, c26, c33, c34, c35, c36, c44, c45, c46, c55, c56, c66;
} tensor6d; 
#endif

#ifdef COVALENT
/* Neighbor table for Tersoff potential */
typedef struct {
#ifndef RING
  real        *dist;
#endif
  short       *typ;
  void        **cl;
  int         *num;
  int         n;
  int         n_max;
} neightab;

typedef neightab* neighptr;
#endif

/* Basic Data Type - The Cell */
typedef struct { 
  vektor *ort;
#ifdef STRAIN
  vektor *dsp;
  vektor *strain;
  vektor *strain_offdia;
  short  *empty;
#elif defined(STRESS)
  vektor *stress;
  vektor *stress_offdia;
  real   *vol;
#else
  int    *sorte;
#endif
#if defined(CONN) || defined(ELCO) || defined(COORD)
  int    *nummer;
#endif
#ifdef COORD
  real    *coord;
#endif
#ifdef COVALENT
  neightab *neightab_array;
#endif
#ifdef ELCO
  tensor   *stress;
  tensor6d *elco;
  real     *vol;
#endif
#ifdef RING
  int      *del;
  int      *hops;
  int      *sp_hops;
  int      *color;
  int      *sp_color;
  neightab * perm_neightab_array;
#endif
#ifdef PS
  real *enc;
#endif
  int     n;
  int     n_max;
} cell;

#ifdef COVALENT
typedef cell* cellptr;
#endif

/* String used for Filenames etc. */
typedef char str255[255];

#if defined(STRESS) || defined(ELCO)
#ifndef TWOD
typedef vektor vektorstr[NUM]; 
#endif
#endif

#ifdef RING
typedef struct atom_ {
  cell *cl;
  int  num;
  int  status;
} atom;

typedef struct queue_elmt {
  cell              *cl;
  int               num;
  struct queue_elmt *next;
} Queue_elmt;
#endif

#ifdef PS
typedef struct list_elmt {
  cell *cl;
  int  num;
  real  d;
  struct list_elmt *next;
} List_elmt;
#endif

/* for parameter reading */
typedef enum ParamType {
  PARAM_STR, PARAM_STRPTR,
  PARAM_INT, PARAM_INT_COPY,
  PARAM_INTEGER, PARAM_INTEGER_COPY,
  PARAM_REAL, PARAM_REAL_COPY
} PARAMTYPE;


/*****************************************************************************
*
*  prototypes
*
*****************************************************************************/

int  main (int argc, char **argv);
void read_parameters(int argc, char **argv);
void read_atoms(str255 infilename);
void usage(void);
void getparamfile(char *infile);
int  getparam(char *param_name, void *param, PARAMTYPE ptype, 
              int pnum_min, int pnum_max);
ivektor cell_coord(vektor pos);
void init_cells(void);
void error(char *mesg);
void alloc_cell(cell *cl, int count);
void do_work(void (*do_cell_pair)(cell *p, cell *q, vektor pbc));
void do_cell_pair(cell *p, cell *q, vektor pbc);
void write_data(void);

#ifdef ANGLE
void do_angle(cell *p, cell *q, cell *qq, vektor pbc, vektor ppbc);
void calc_angles(void);
#endif

#ifdef STRAIN
void read_displacement(str255 infilename);
void calc_strain(void);
#endif

#if defined(STRESS) || defined(ELCO)
void voronoi(void);
#ifdef TWOD
void do_voronoi_2d(void);
#else
void do_voronoi_3d(void);
#endif
void sort(void);
#endif

#ifdef STRESS
void read_stress(str255 infilename);
void write_stress(int restart);
#endif

#ifdef TORSION
void do_angle(cell *p, cell *q, cell *r, cell *s, 
              vektor pbc_q, vektor pbc_r, vektor pbc_s);
void calc_angles(void);
#endif

#ifdef CONN
void make_numbers(void);
#endif

#ifdef COVALENT
void do_neighbour_tables(cell *p, cell *q, vektor pbc);
#endif

#ifdef TERSOFF
void init_tersoff(void);
#elif STIWEB
void init_stiweb(void);
#elif KEATING
void init_keating(void);
#endif

#ifdef ELCO
void init_elco(void);

#ifdef TERSOFF
void do_elco_tersoff(void);
#elif STIWEB
void do_elco_stiweb(void);
#elif KEATING
void do_elco_keating(void);
#endif
void write_stress(void);
void write_elco(void);
#endif

#ifdef RING
void search_rings(void);
void go_forward(void);
void compute_hops(void);
int sp_ring(void);
void update_neighbour_tables(cell *p, int i);
void write_data(void);
Queue_elmt *queue_create(cell *cl, int num);
Queue_elmt *queue_enqueue(Queue_elmt *queue, cell *cl, int num);
void queue_dequeue(Queue_elmt **qptr, cell **clptr, int *numptr); 
#endif

#ifdef PS
void init_ps(void);
void sort_atoms(void);
void printcolor(FILE *out, real type, real zval, int bond);
void setlinew(FILE *out, real linew);
real setradius(int type);
vektor2d transform(vektor2d vekt, real scl, vektor2d trans);
void write_profile(void);
void write_settings(void);
vektor maxvektor(vektor v1, vektor v2);
vektor minvektor(vektor v1, vektor v2);
vektor zrotate(vektor pos, real angle);
void draw_colorencoding(FILE *out, vektor2d position, vektor2d size); 
#ifndef TWOD
void draw_picture3d(void);
real fdistance(vektor pos);
real pdistance(vektor pos);
real scale(vektor pos, real sclfct);
vektor2d proj(vektor pos, real sclfct, vektor2d trans);
vektor xrotate(vektor pos, real angle);
vektor yrotate(vektor pos, real angle);
#else
void draw_picture2d(void);
vektor2d scl(vektor pos, real sclfct, vektor2d trans);
#endif
#endif

/*****************************************************************************
*
*  global variables
*
*****************************************************************************/

cell *cell_array;    /* Array of Cells */
ivektor cell_dim;    /* Dimension of above */
int natoms;          /* Total number of atoms */
int ntypes=0;        /* Total number of different atom types */
int vtypes=0;        /* Number of virtual types */
int use_vtypes=0;    /* flag for using virtual types */
str255 progname;     /* Name of current executable argv[0] */
int curline;         /* Number of current line for parameter reading */
str255 error_msg;    /* string for error message */

/* The simulation box and its inverse */
vektor box_x, tbox_x;
vektor box_y, tbox_y;
#ifndef TWOD
vektor box_z, tbox_z;
#endif
vektor height;
#ifdef TWOD
#ifdef PS
ivektor pbc_dirs = {0,0};
#else 
ivektor pbc_dirs = {1,1};    /* directions with pbc - default is PBC */
#endif
#else
#ifdef PS
ivektor pbc_dirs = {0,0,0};
#else
ivektor pbc_dirs = {1,1,1};  /* directions with pbc - default is PBC */
#endif
#endif

/* Filenames */
str255 infilename;    /* Input File */
str255 outfilename;   /* Output File */
char *paramfilename;  /* parameter file */

/* Bookkeeping */
int  restart = -1;
int  avpos   = -1;
real r2_cut;
real r_min = 1.0;     /* default value */
#ifdef PAIR
real r_max = 5.0;     /* default value */
#else
real r_max = 1.0;     /* default value */
#endif
int  nangles = 0;
#ifdef CONN
real *r_cut;
int  n_min=0, n_max=0;
#endif
#ifdef STRAIN
real r_cell = 1.0;    /* default value */
#endif
#if defined(STRESS) || defined(ELCO)
int  atomcount = 0;   /* statistics */ 
int  maxneigh  = 0, sumneigh  = 0;
int  maxvert   = 0, sumvert   = 0;
int  maxedges  = 0, sumedges  = 0;
#ifndef TWOD
int  maxfaces  = 0, sumfaces  = 0;
#endif
#endif

/* The histograms */
#ifndef RING
real      *histogram;
#endif
#ifdef PAIR
ivektor3d hist_dim;
#endif
#ifdef TORSION
ivektor5d hist_dim;
#endif
#ifdef ANGLE
ivektor4d hist_dim;
#endif

#if defined(STRESS) || defined(ELCO)
int     neighnum;        /* Number of neighbour atoms */
vektor  *candcoord;      /* Coordinates of neighbour atoms */
real    *canddist2;      /* Distance squared of neighbour atoms */ 
#ifndef TWOD
real    volume;          /* Volume of Voronoi cell, 3d */
#else
real    area;            /* Area of Voronoi cell, 2d */
#endif
#endif

#ifdef COORD
real      *numbers;
int       use_unity = -1;
int       local     = -1;
int       global    = -1;
int       c_max     = 10;
ivektor2d num_dim;
int       bonds     = 0; /* Total number of bonds */
#endif

#ifdef CONN
/* Connection matrix */
ivektor2d cm_dim = {0,MAXNEIGH};
int *cm, *nn, *tp, *num, *ind;
vektor *pos;
#endif

#ifdef ELCO
tensor6d c = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
tensor sigma = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
real   r_cell = -1.0;
real   vol, epot = 0.0, bulkm = 0.0, dbulkm_dp = 0.0;
int    stresstens = 0, moduli = 0, all_moduli = 0;
#endif

#ifdef COVALENT
int  neigh_len = 6;
int  nbonds = 0;
#endif
#ifdef TERSOFF
real ter_r_cut[10][10], ter_r2_cut[10][10], ters_r_cut[55] ;
real ter_r0[10][10], ters_r0[55];
real ter_a[10][10], ters_a[55];
real ter_b[10][10], ters_b[55];
real ter_la[10][10], ters_la[55];
real ter_mu[10][10], ters_mu[55];
real ters_ga[10];
real ters_n[10];
real ter_c2[10], ters_c[10];
real ter_d2[10], ters_d[10];
real ters_h[10];
real ter_chi[10][10], ters_chi[45];
real ter_om[10][10], ters_om[45];
#elif STIWEB
real stiweb_a[55];
real sw_a[10][10];
real stiweb_b[55];
real sw_b[10][10];
real stiweb_p[55];
real sw_p[10][10];
real stiweb_q[55];
real sw_q[10][10];
real stiweb_a1[55];
real sw_a1[10][10];
real stiweb_de[55];
real sw_de[10][10];
real sw_2_a1[10][10];
real stiweb_a2[55];
real sw_a2[10][10];
real stiweb_la[550];
real sw_la[10][10][10];
real stiweb_ga[55];
real sw_ga[10][10];
real sw_r_cut[10][10];
#elif KEATING
real keating_alpha[55];
real keat_alpha[10][10];
real keating_d[55];
real keat_d[10][10];
real keating_r_cut[55];
real keat_r_cut[10][10];
real keat_r2_cut[10][10];
real keating_beta[550];
real keat_beta[10][10][10];
#endif

#ifdef RING
int *histogram;
int total_rings = 0; /* Number of sp rings found */
int first = 1;       /* =1: First generation of neighbour tables */
int max_length = 10; /* Maximum length of paths, default 10 */
Queue_elmt *queue;
int queue_length;
Queue_elmt *sp_queue;
int sp_queue_length;
atom *stack;
int stack_end;
#endif

#ifdef PS
List_elmt *atomlist;
int    natoms_drawn = 0;
int    nbonds_drawn = 0;
real   dist_min = Min, dist_max = Max;
real   r_atom = -1.0; 
real   r_bond = -1.0;
long   atomcolor = 1;
real   atombrightness = 0.0;
real   bondbrightness = -0.3;
long    bondcolor = 0;
real   ambient = 1.0;
long   radii = 0;
long   eng   = 0;
long   kineng = 0;
long   color_enc = 0;
int    colorencoding = 0;
real   linew = -1.0;
real   alinew = 1.0;
real   blinew = 1.0;
real   slinew = 1.0;
real   llinew = -1.0;
real   shading = 45.0;
real   scaling = 1.0;
real   frame = 0.0;
real   frame_linew = 2.0;
real   rframe = -1.0;
real   rframe_linew = 1.0;
real   sframe = -1.0;
real   sframe_linew = 1.0;
real   frame_border = 9.0;
real   frame_border_linew = 0.0;
vektor2d translation = {0.0, 0.0};
real   zoom = 1.0;
real   ps_width = -1.0;
real   ps_height = -1.0;
int    axes = -1;
real   spect = -1.0;
real   backgrd = 0.0;
real   needle = -1.0;
real   depth  = 0.0;
real   user_scale = 1.0;
real   wireframe = -1.0;
real   crop = -1.0;
int    settings = -1;
vektor ursprung;
vektor cunit;
#ifndef TWOD
vektor maxl={Max,Max,Max}, minl={Min,Min,Min}; 
vektor real_maxl={Max,Max,Max}, real_minl={Min,Min,Min};
#else
vektor maxl={Max,Max}, minl={Min,Min}; 
vektor real_maxl={Max,Max}, real_minl={Min,Min};
#endif
real   maxeng=Max, mineng=Min;
real   angx = 0.0, angy = 0.0, angz = 0.0;
#ifndef TWOD
real   zeta = Max;
vektor foc = { 0.0, 0.0, Max};
#endif
real   r_cut_vec[55] = { -1.0 }; 
real   ters_r_cut[55] = { -1.0 };
real   r_cut[10][10];
real   r2cut[10][10];
vektor3d color_array[9];
vektor3d enc_array[8];
real   shade_array[9];
real   radii_array[9];
int    atomcol[9] = {0,1,2,3,4,5,6,7,8};
int    bondcol[9] = {0,1,2,3,4,5,6,7,8};
int    rad[9];
int    enc[8] = {1,2,3,4,5,6,7,8};
int    enc_len;
#ifndef TWOD
vektor unitv[3] = {{1,0,0},{0,1,0},{0,0,1}};
#else
vektor unitv[2] = {{1,0},{0,1}};
#endif
#endif

/******************************************************************************
*
*  Common functions
*
******************************************************************************/

/******************************************************************************
*
*  To determine the cell into which a given particle belongs, we
*  have to transform the cartesian coordinates of the particle into
*  the coordinate system spanned by the vectors of the box edges. 
*  This yields coordinates in the interval [0..1] that are
*  multiplied by cell_dim to get the cell's index.
*
******************************************************************************/

#ifdef TWOD

/* compute transformation matrix */
void make_box(void)
{
  real volume;

  /* volume */
  volume = box_x.x * box_y.y - box_x.y * box_y.x;
  if (0==volume) error("Box Edges are parallel.");

  /* compute tbox_j such that SPROD(box_i,tbox_j) == delta_ij */
  tbox_x.x =   box_y.y / volume;
  tbox_x.y = - box_y.x / volume;
  tbox_y.x = - box_x.y / volume;
  tbox_y.y =   box_x.x / volume;

  /* squares of the box heights perpendicular to the faces */
  height.x = 1.0 / SPROD(tbox_x,tbox_x);
  height.y = 1.0 / SPROD(tbox_y,tbox_y);
}

#else

/* vector product */ 
vektor vec_prod(vektor u, vektor v)
{
  vektor w;
  w.x = u.y * v.z - u.z * v.y;
  w.y = u.z * v.x - u.x * v.z;
  w.z = u.x * v.y - u.y * v.x;
  return w;
}

/* compute transformation matrix */
void make_box(void)
{
  real volume;

  /* compute tbox_j such that SPROD(box_i,tbox_j) == delta_ij */
  /* first unnormalized */
  tbox_x = vec_prod( box_y, box_z );
  tbox_y = vec_prod( box_z, box_x );
  tbox_z = vec_prod( box_x, box_y );

  /* volume */
  volume = SPROD( box_x, tbox_x );
  if (0==volume) error("Box Edges are parallel.");

#ifdef ELCO
  vol = volume;
#endif

  /* normalization */
  tbox_x.x /= volume;  tbox_x.y /= volume;  tbox_x.z /= volume;
  tbox_y.x /= volume;  tbox_y.y /= volume;  tbox_y.z /= volume;
  tbox_z.x /= volume;  tbox_z.y /= volume;  tbox_z.z /= volume;

  /* squares of the box heights perpendicular to the faces */
  height.x = 1.0 / SPROD(tbox_x,tbox_x);
  height.y = 1.0 / SPROD(tbox_y,tbox_y);
  height.z = 1.0 / SPROD(tbox_z,tbox_z);
}

#endif


/******************************************************************************
*
*  Compute the size of the cells, and initialize the cell array.
*  There must be at least two cells in each direction.
*
******************************************************************************/

void init_cells(void)
{
  vektor cell_scale;
  cell   *p;
  int    nmax, i;

  make_box();

  /* scaling factors box/cell */
  cell_scale.x = sqrt( (double)(r2_cut / height.x) );
  cell_scale.y = sqrt( (double)(r2_cut / height.y) );
#ifndef TWOD
  cell_scale.z = sqrt( (double)(r2_cut / height.z) );
#endif

#ifdef TWOD
  printf("Minimal cell size: \n\t ( %f %f ) \n\t ( %f %f ) \n",
	 box_x.x * cell_scale.x, box_x.y * cell_scale.x, 
	 box_y.x * cell_scale.y, box_y.y * cell_scale.y); 
#else
  printf("Minimal cell size: \n");
  printf("\t ( %f %f %f ) \n\t ( %f %f %f ) \n\t ( %f %f %f )\n",
     box_x.x * cell_scale.x, box_x.y * cell_scale.x, box_x.z * cell_scale.x,
     box_y.x * cell_scale.y, box_y.y * cell_scale.y, box_y.z * cell_scale.y,
     box_z.x * cell_scale.z, box_z.y * cell_scale.z, box_z.z * cell_scale.z);
#endif

  /* set up cell array */
  cell_dim.x = (int) ( 1.0 / cell_scale.x );
  cell_dim.y = (int) ( 1.0 / cell_scale.y );
#ifndef TWOD
  cell_dim.z = (int) ( 1.0 / cell_scale.z );
#endif
  
  /* If an integer number of cells does not fit exactly into the box, the
     cells are enlarged accordingly */
  cell_scale.x = 1.0 / cell_dim.x;
  cell_scale.y = 1.0 / cell_dim.y;
#ifndef TWOD
  cell_scale.z = 1.0 / cell_dim.z;
#endif

#ifdef TWOD
  printf("Actual cell size: \n\t ( %f %f ) \n\t ( %f %f ) \n",
	 box_x.x * cell_scale.x, box_x.y * cell_scale.x, 
	 box_y.x * cell_scale.y, box_y.y * cell_scale.y); 
  printf("Cell array dimensions: %d %d\n", cell_dim.x,cell_dim.y);
  cell_array = (cell *) malloc(cell_dim.x * cell_dim.y * sizeof(cell));
#else
  printf("Actual cell size: \n");
  printf("\t ( %f %f %f ) \n\t ( %f %f %f ) \n\t ( %f %f %f )\n",
     box_x.x * cell_scale.x, box_x.y * cell_scale.x, box_x.z * cell_scale.x,
     box_y.x * cell_scale.y, box_y.y * cell_scale.y, box_y.z * cell_scale.y,
     box_z.x * cell_scale.z, box_z.y * cell_scale.z, box_z.z * cell_scale.z);
  printf("Cell array dimensions: %d %d %d\n",cell_dim.x,cell_dim.y,cell_dim.z);
  cell_array = (cell *) malloc(cell_dim.x*cell_dim.y*cell_dim.z*sizeof(cell));
#endif
  if (NULL == cell_array) error("Cannot allocate memory for cells");

  nmax = cell_dim.x * cell_dim.y;
#ifndef TWOD
  nmax = nmax * cell_dim.z;
#endif

  /* initialize cells */
  for (i=0; i<nmax; ++i) {
    p = cell_array + i;
    p->n_max = 0;
    alloc_cell(p, CSTEP);
  }

}


/*****************************************************************************
*
*  cell_coord gives the (integral) cell_coorinates of a position
*
*****************************************************************************/

ivektor cell_coord(vektor ort)
{
  ivektor coord;

  /* Map positions to boxes */
  coord.x = (int) TRUNC( cell_dim.x * SPROD(ort,tbox_x) );
  coord.y = (int) TRUNC( cell_dim.y * SPROD(ort,tbox_y) );
#ifndef TWOD
  coord.z = (int) TRUNC( cell_dim.z * SPROD(ort,tbox_z) );
#endif

  /* Roundoff errors put atoms slightly out of the simulation cell */
  coord.x = coord.x <   0         ?             0 : coord.x;
  coord.x = coord.x >= cell_dim.x ? cell_dim.x -1 : coord.x;
  coord.y = coord.y <   0         ?             0 : coord.y;
  coord.y = coord.y >= cell_dim.y ? cell_dim.y -1 : coord.y;
#ifndef TWOD
  coord.z = coord.z <   0         ?             0 : coord.z;
  coord.z = coord.z >= cell_dim.z ? cell_dim.z -1 : coord.z;
#endif
  return coord;
}


/******************************************************************************
*
*  allocate memory for a cell
*
******************************************************************************/

void alloc_cell(cell *cl, int count)
{
  int i;
#ifdef COVALENT
  neightab *neigh;
#endif

  /* initialize if it is the first call */
  if (cl->n_max == 0) {
    cl->ort           = NULL;
#ifdef STRAIN
    cl->dsp           = NULL;
    cl->strain        = NULL;
    cl->strain_offdia = NULL;
    cl->empty         = NULL;
#elif defined(STRESS)
    cl->stress        = NULL;
    cl->stress_offdia = NULL;
    cl->vol           = NULL;
#else
    cl->sorte         = NULL;
#endif
#if defined(CONN) || defined(ELCO) || defined(COORD)
    cl->nummer        = NULL;
#endif
#ifdef COORD
    cl->coord         = NULL;
#endif
#ifdef COVALENT
    cl->neightab_array= NULL;
#endif
#ifdef ELCO
    cl->stress        = NULL;
    cl->elco          = NULL;
    cl->vol           = NULL;
#endif
#ifdef RING
    cl->del            = NULL;
    cl->hops           = NULL;
    cl->sp_hops        = NULL;
    cl->color          = NULL;
    cl->sp_color       = NULL;
    cl->perm_neightab_array = NULL;
#endif
#ifdef PS
    cl->enc            = NULL;
#endif
    cl->n             = 0;
  }

  /* (re)allocate */
  cl->ort    = (vektor *) realloc(cl->ort,    count * sizeof(vektor));
#ifdef STRAIN
  cl->dsp    = (vektor *) realloc(cl->dsp,    count * sizeof(vektor));
  cl->strain = (vektor *) realloc(cl->strain, count * sizeof(vektor));
  cl->strain_offdia = (vektor *) realloc(cl->strain_offdia,
                                              count * sizeof(vektor));
  cl->empty  = (short  *) realloc(cl->empty,  count * sizeof(short));
#elif defined(STRESS)
  cl->stress = (vektor *) realloc(cl->stress, count * sizeof(vektor));
  cl->stress_offdia = (vektor *) realloc(cl->stress_offdia,
                                              count * sizeof(vektor));
  cl->vol    = (real   *) realloc(cl->vol,    count * sizeof(real));
#else
  cl->sorte  = (int    *) realloc(cl->sorte,  count * sizeof(int));
#endif
#if defined(CONN) || defined(ELCO) || defined(COORD)
  cl->nummer = (int    *) realloc(cl->nummer, count * sizeof(int));
#endif
#ifdef COORD
  cl->coord  = (real   *) realloc(cl->coord,  ntypes * count * sizeof(real));
#endif
#ifdef ELCO
  cl->stress = (tensor *) realloc(cl->stress, count * sizeof(tensor));
  cl->elco   = (tensor6d *) realloc(cl->elco, count * sizeof(tensor6d));
  cl->vol    = (real   *) realloc(cl->vol,    count * sizeof(real));
#endif
#ifdef RING
  cl->del      = (int    *) realloc(cl->del,      count * sizeof(int));
  cl->hops     = (int    *) realloc(cl->hops,     count * sizeof(int));
  cl->sp_hops  = (int    *) realloc(cl->sp_hops,  count * sizeof(int));
  cl->color    = (int    *) realloc(cl->color,    count * sizeof(int));
  cl->sp_color = (int    *) realloc(cl->sp_color, count * sizeof(int));
#endif
#ifdef PS
  cl->enc     = (real   *) realloc(cl->enc,     count * sizeof(real) ); 
#endif
#ifdef COVALENT
  cl->neightab_array = (neightab *) realloc( cl->neightab_array, 
					      count * sizeof(neightab));
   if (NULL == cl->neightab_array) 
      error("Cannot allocate neighbor tables");

   /* Allocate memory for neighbour tables */
   for (i = cl->n_max; i < count; ++i) {
     neigh = cl->neightab_array + i;

     neigh->n     = 0;
     neigh->n_max = neigh_len;
#ifndef RING
     neigh->dist  = (real *)  malloc( neigh_len * 3 * sizeof(real) );
#endif
     neigh->typ   = (short *) malloc( neigh_len * sizeof(short) );
     neigh->cl    = (void **) malloc( neigh_len * sizeof(cellptr) );
     neigh->num   = (int *)   malloc( neigh_len * sizeof(int) );

     if (
#ifndef RING
       (neigh->dist==NULL) || 
#endif
       (neigh->typ==NULL) || 
       (neigh->cl  ==NULL) || (neigh->num==NULL) )
       error("Cannot allocate memory for neighbor table");
   }
#endif
#ifdef RING
  cl->perm_neightab_array = (neightab *) realloc( cl->perm_neightab_array, 
					      count * sizeof(neightab));
   if (NULL == cl->perm_neightab_array) 
      error("Cannot allocate permanent neighbor tables");

   /* Allocate memory for permanent neighbour tables */
   for (i = cl->n_max; i < count; ++i) {
     neigh = cl->perm_neightab_array + i;

     neigh->n     = 0;
     neigh->n_max = neigh_len;
     neigh->typ   = (short *) malloc( neigh_len * sizeof(short) );
     neigh->cl    = (void **) malloc( neigh_len * sizeof(cellptr) );
     neigh->num   = (int *)   malloc( neigh_len * sizeof(int) );

     if ((neigh->typ==NULL) || (neigh->cl  ==NULL) || 
	 (neigh->num==NULL) )
       error("Cannot allocate memory for permanent neighbor table");
   }
#endif

  /* check if it worked */
  if ( (NULL==cl->ort)
#ifdef STRAIN
    || (NULL==cl->dsp)
    || (NULL==cl->strain)
    || (NULL==cl->strain_offdia)
    || (NULL==cl->empty)
#elif defined(STRESS)
    || (NULL==cl->stress)
    || (NULL==cl->stress_offdia)
    || (NULL==cl->vol)
#else   
    || (NULL==cl->sorte)
#endif
#if defined(CONN) || defined(ELCO) || defined(COORD)
    || (NULL==cl->nummer)
#endif
#ifdef COORD
    || (NULL==cl->coord)
#endif
#ifdef ELCO
    || (NULL==cl->stress)
    || (NULL==cl->elco)
    || (NULL==cl->vol)
#endif
#ifdef RING
    || (NULL==cl->del)
    || (NULL==cl->hops)
    || (NULL==cl->sp_hops)
    || (NULL==cl->color)
    || (NULL==cl->sp_color)
#endif
#ifdef PS
   || (NULL==cl->enc)
#endif
  ) error("Cannot allocate memory for cell.");
  cl->n_max  = count;
}


/******************************************************************************
*
*   read_parameters reads the command line parameters, 
*   and the parameter file given on the command line.
*
******************************************************************************/

void read_parameters(int argc,char **argv)
{
  str255 fname;
  FILE *testfile;

/* Check for Restart, process options */

  strcpy(progname,argv[0]);
  while ((argc > 1) && (argv[1][0] =='-')) {
    switch (argv[1][1]) {
      /* r - restart */
    case 'r':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          restart = atoi(argv[2]);
          --argc;
          ++argv;
        }
      }
      else restart = atoi(&argv[1][2]);
      break;
     /* A - use avpos outfiles */
    case 'A':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          avpos = atoi(argv[2]);
          --argc;
          ++argv;
        }
      }
      else avpos = atoi(&argv[1][2]);
      break;  
#ifdef PAIR
      /* a - minimum radius */
    case 'a':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          r_min = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else r_min = atof(&argv[1][2]);
      break;
#endif
#ifndef CONN
      /* e - maximum radius */
    case 'e':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          r_max = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else r_max = atof(&argv[1][2]);
#ifdef STRAIN
      r_cell = r_max;
#endif
      break;
#endif
#ifdef COORD
      /* u - use cutoff-function equal to 1 */
    case 'u':
      use_unity = 1;
      break;
     /* l - write local coordination numbers */
    case 'l':
      local  = 1;
      break;
     /* g - write global coordination numbers */
    case 'g':
      global = 1;
      break; 
      /* C - largest coordination number */
    case 'C':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          c_max = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else c_max = atof(&argv[1][2]);
      break;
#endif 
#if defined(STRAIN) || defined(ELCO)
      /* c - cell width */
    case 'c':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          r_cell = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else r_cell = atof(&argv[1][2]);
      break;
#endif
#ifdef ELCO
      /* s - Ouput local stress tensor */
    case 's':
      stresstens = 1;
      break;
      /* m - Output local elastic moduli c11, c12, c44 */
    case 'm':
      moduli = 1;
      break;
      /* M - Output all local elastic moduli */
    case 'M':
      all_moduli = 1;
      break;
#endif
#ifdef RING
      /* l - Maximal length of rings */
    case 'l':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          max_length = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else max_length = atof(&argv[1][2]);
      break;
      /* n - Size of neighbour tables */
    case 'n':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          neigh_len = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else neigh_len = atof(&argv[1][2]);
      break;
#endif
#ifdef PS
      /* a - draw coordinate axes */
    case 'a':
     if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          axes = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else axes = atof(&argv[1][2]);
      break;
      /* b - brightness of bonds */
    case 'b':
     if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          bondbrightness = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else bondbrightness = atof(&argv[1][2]);
      break;
      /* B - bond radius */
    case 'B':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          r_bond = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else r_bond = atof(&argv[1][2]);
      break;
      /* c - crop boundary strip */
    case 'c':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          crop = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else crop = atof(&argv[1][2]);
      break;
      /* C - atom colors */
    case 'C':
     if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          atomcolor = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else atomcolor = atof(&argv[1][2]);
      break;
#ifndef TWOD
      /* d - shading with depth */
    case 'd':
     if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          depth = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else depth = atof(&argv[1][2]);
      break;
#endif
      /* D - atom radii */
    case 'D':
     if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          radii = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else radii = atof(&argv[1][2]);
      break; 
      /* E - Use color encoding for energy */
    case 'E':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          eng = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else eng = atof(&argv[1][2]);
      break;
      /* f - draw boxes */
    case 'f':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          frame = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else frame = atof(&argv[1][2]);
      break;
      /* F - linewidth of bounding box */
    case 'F':
     if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          frame_linew = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else frame_linew = atof(&argv[1][2]);
      break; 
      /* g - background color */
    case 'g':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          backgrd = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else backgrd = atof(&argv[1][2]);
      break;
       /* G - colorencoding of column G */ 
    case 'G':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          color_enc = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else color_enc = atof(&argv[1][2]);
      break;
      /* h - height of picture */
    case 'h':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          ps_height = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else ps_height = atof(&argv[1][2]);
      break;
      /* H - linewidth of atoms, bonds, and 
       bond separators */
    case 'H':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          linew = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else linew = atof(&argv[1][2]);
      break;
      /* i - x translation */
    case 'i':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          translation.x = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else translation.x = atof(&argv[1][2]);
      break;
      /* I - Brigthness of atom shading */
    case 'I':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          ambient = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else ambient = atof(&argv[1][2]);
      break;
      /* j - y translation */
    case 'j':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          translation.y = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else translation.y = atof(&argv[1][2]);
      break;
      /* J - brightness of atoms */
    case 'J':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          atombrightness = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else atombrightness = atof(&argv[1][2]);
      break;
     /* k - show color encoding of energy */
    case 'k':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          spect = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else spect = atof(&argv[1][2]);
      break;
     /* K write file with parameter settings */
    case 'K':
      settings = 1;
      break;
     /* l linewidth of bond borders */
    case 'l':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          blinew = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else blinew = atof(&argv[1][2]);
      break;
     /* L - line width of atom borders */
    case 'L':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          alinew = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else alinew = atof(&argv[1][2]);
      break;
      /* m - zoom */
    case 'm':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          zoom = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else zoom = atof(&argv[1][2]);
      break;
      /* M - Scale factor with depth */
    case 'M':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          user_scale = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else user_scale = atof(&argv[1][2]);
      break;
#ifndef TWOD
      /* n - draw lightshadow */
    case 'n':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          llinew = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else llinew = atof(&argv[1][2]);
      break;
     /* N - bonds like needles */
    case 'N':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          needle = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else needle = atof(&argv[1][2]);
      break;
#endif
     /* o - draw simulation box */
    case 'o':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          sframe = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else sframe = atof(&argv[1][2]);
      break;
     /* O - linewidth of simulation box */
    case 'O':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          sframe_linew = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else sframe_linew = atof(&argv[1][2]);
      break;
     /* P - linewidth of bond separators */
    case 'P':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          slinew = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else slinew = atof(&argv[1][2]);
      break;
     /* q - color of border of frames */
    case 'q':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          frame_border = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else frame_border = atof(&argv[1][2]);
      break;
     /* Q - linewidth of border of frames */
    case 'Q':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          frame_border_linew = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else frame_border_linew = atof(&argv[1][2]);
      break;
      /* R - atom radius */
    case 'R':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          r_atom = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else r_atom = atof(&argv[1][2]);
      break;
      /* s - scale factor */
    case 's':
     if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          scaling = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else scaling = atof(&argv[1][2]);
      break; 
      /* S - shading of atoms */
    case 'S':
     if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          shading = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else shading = atof(&argv[1][2]);
      break; 
      /* t - kinetic energy */
    case 't':
     if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          kineng = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else kineng = atof(&argv[1][2]);
      break; 
#ifndef TWOD
      /* T - Locatio of projection plane */
    case 'T':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          zeta = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else zeta = atof(&argv[1][2]);
      break;
#endif
      /* u - draw real bounding box */
    case 'u':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          rframe = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else rframe = atof(&argv[1][2]);
      break;
      /* U - draw real bounding box */
    case 'U':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          rframe_linew = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else rframe_linew = atof(&argv[1][2]);
      break;
      /* V - color of bonds */
    case 'V':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          bondcolor = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else bondcolor = atof(&argv[1][2]);
      break;
      /* w - width of picture */
    case 'w':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          ps_width = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else ps_width = atof(&argv[1][2]);
      break;
     /* W - draw wireframe */
    case 'W':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          wireframe = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else wireframe = atof(&argv[1][2]);
      break;
#ifndef TWOD 
      /* x - Rotation */
    case 'x':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          angx = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else angx = atof(&argv[1][2]);
      break;
      /* y - Rotation */
    case 'y':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          angy = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else angy = atof(&argv[1][2]);
      break;
#endif
      /* z - Rotation */
    case 'z':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          angz = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else angz = atof(&argv[1][2]);
      break;
#ifndef TWOD
      /* X - Center of projection */
    case 'X':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          foc.x = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else foc.x = atof(&argv[1][2]);
      break;
     /* Y - Center of projection */
    case 'Y':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          foc.y = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else foc.y = atof(&argv[1][2]);
      break;
      /* Z - Center of projection */
    case 'Z':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          foc.z = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else foc.z = atof(&argv[1][2]);
      break;
#endif
#endif
    case 'p':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          paramfilename = strdup(argv[2]);
          --argc;
          ++argv;
        }
      }
      else paramfilename = strdup(&argv[1][2]);
      break;
    case 'v':
      use_vtypes=1;
      break;
    default:
      printf("Illegal option %s \n",argv[1]);
      usage();
      exit(-1);
    }
    ++argv;
    --argc;
  }

  getparamfile(paramfilename);
  if (use_vtypes==1) ntypes = vtypes;

  /* Get restart parameters if restart */
  if (-1 != restart) {
    sprintf(fname,"%s.%d.itr",outfilename,restart);
    testfile = fopen(fname,"r");
    if (NULL==testfile) { 
      sprintf(fname,"%s.%05d.itr",outfilename,restart);
    } else {
      fclose(testfile);
    }
    getparamfile(fname);
  }
  else if  (-1 != avpos) {
    sprintf(fname,"%s.%d.avp.itr",outfilename,avpos);
    testfile = fopen(fname,"r");
    if (NULL==testfile) { 
      sprintf(fname,"%s.%05d.avp.itr",outfilename,avpos);
    } else {
      fclose(testfile);
    }
    getparamfile(fname);
  }

#ifdef STRAIN
  if (-1 == restart) {
    sprintf(infilename, "%s.dsp", outfilename);
  }
  else {
    sprintf(infilename,"%s.%u.dsp",outfilename,restart);
    testfile = fopen(infilename,"r");
    if (NULL==testfile) { 
      sprintf(infilename,"%s.%05d.dsp",outfilename,restart);
    } else {
      fclose(testfile);
    }
  }
#elif defined(STRESS)
  if (-1 == restart) {
    sprintf(infilename, "%s.press",outfilename);
  }
  else {
    sprintf(infilename,"%s.%u.press",outfilename,restart);
    testfile = fopen(infilename,"r");
    if (NULL==testfile) { 
      sprintf(infilename,"%s.%05d.press",outfilename,restart);
    } else {
      fclose(testfile);
    }
  }
#else
  if (-1 != restart) {
    sprintf(infilename,"%s.%u.chkpt",outfilename,restart);
    testfile = fopen(infilename,"r");
    if (NULL==testfile) { 
      sprintf(infilename,"%s.%05d.chkpt",outfilename,restart);
    } else {
      fclose(testfile);
    }
  }
  else if (-1 != avpos) {
    sprintf(infilename,"%s.%u.avp",outfilename,avpos);
    testfile = fopen(infilename,"r");
    if (NULL==testfile) { 
      sprintf(infilename,"%s.%05d.avp",outfilename,avpos);
    } else {
      fclose(testfile);
    }
  }
#endif

#ifdef STRAIN
  /* r_cell >= r_max is required */
  if (r_cell < r_max) error("Cell smaller than cutoff radius!");
  r2_cut = SQR(r_cell);
#endif
#if defined(STRESS) || defined (RING)
  r2_cut = SQR(r_max);  
#endif

} 


/*****************************************************************************
*
*  read one parameter (routine from imd_param.c)
*
*****************************************************************************/

int getparam(char *param_name, void *param, PARAMTYPE ptype, 
             int pnum_min, int pnum_max)
{
  static char errmsg[256];
  char *str;
  int i;
  int numread;

  numread = 0;
  if (ptype == PARAM_STR) {
    str = strtok(NULL," \t\n");
    if (str == NULL) {
      fprintf(stderr,"parameter for %s missing in line %u\n",
              param_name,curline);
      error("string expected");
    }
    else strncpy((char *)param,str,pnum_max);
    numread++;
  }
  else if (ptype == PARAM_STRPTR) {
    str = strtok(NULL," \t\n");
    if (str == NULL) {
      fprintf(stderr,"parameter for %s missing in line %u\n",
              param_name,curline);
      error("string expected");
    }
    else *((char**)param) = strdup(str);
    numread++;
  }
  else if (ptype == PARAM_INT) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," \t\n");
      if (str == NULL) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"integer vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      else ((int*)param)[i] = atoi(str);
      numread++;
    }
    for (i=pnum_min; i<pnum_max; i++) {
      if ((str = strtok(NULL," \t\n")) != NULL) {
        ((int*)param)[i] = atoi(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_INT_COPY) {
    int ival = 0;
    for (i=0; i<pnum_max; i++) {
      str = strtok(NULL," \t\n");
      if (str != NULL) {
        ival = atoi(str);
        numread++; /* return number of parameters actually read */
      }
      else if (i<pnum_min) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"integer vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      ((int*)param)[i] = ival;
    }
  }
  else if (ptype == PARAM_INTEGER) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," \t\n");
      if (str == NULL) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"integer vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      else ((int *)param)[i] = atoi(str);
      numread++;
    }
    for (i=pnum_min; i<pnum_max; i++) {
      if ((str = strtok(NULL," \t\n")) != NULL) {
        ((int *)param)[i] = atoi(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_INTEGER_COPY) {
    int ival = 0;
    for (i=0; i<pnum_max; i++) {
      str = strtok(NULL," \t\n");
      if (str != NULL) {
        ival = atoi(str);
        numread++; /* return number of parameters actually read */
      }
      else if (i<pnum_min) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"integer vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      ((int *)param)[i] = (int)ival;
    }
  }
  else if (ptype == PARAM_REAL) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," \t\n");
      if (str == NULL) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"real vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      else ((real*)param)[i] = atof(str);
      numread++;
    }
    for (i=pnum_min; i<pnum_max; i++) {
      if ((str = strtok(NULL," \t\n")) != NULL) {
        ((real*)param)[i] = atof(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_REAL_COPY) {
    real rval = 0;
    for (i=0; i<pnum_max; i++) {
      str = strtok(NULL," \t\n");
      if (str != NULL) {
        rval = atof(str);
        numread++; /* return number of parameters actually read */
      }
      else if (i<pnum_min) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"real vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      ((real*)param)[i] = rval;
    }
  }
  return numread;
} /* getparam */


/*****************************************************************************
*
*  read tag-based parameter file
*  lines beginning with comment characters '#' or blank lines are skipped
*
*****************************************************************************/

void getparamfile(char *paramfname)
{
  FILE *pf;
  char buffer[1024];
  char *token;
  char *res;

  curline = 0;
  pf = fopen(paramfname,"r");
  if (NULL == pf) {
    sprintf(error_msg,"Cannot open parameter file %s",paramfname);
    error(error_msg);
  }

  do {
    res=fgets(buffer,1024,pf);
    if (NULL == res) break; /* probably EOF reached */
    curline++;
    token = strtok(buffer," \t\n");
    if (NULL == token) continue; /* skip blank lines */
    if (token[0]=='#') continue; /* skip comments */

    if (strcasecmp(token,"coordname")==0) {
      /* file name for atom coordinate input data */
      getparam("coordname",infilename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"outfiles")==0) {
      /* output file basename */
      getparam("outfiles",outfilename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"box_x")==0) {
      /* 'x' or first vector for box */
      getparam("box_x",&box_x,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"box_y")==0) {
      /* 'y' or second vector for box */
      getparam("box_y",&box_y,PARAM_REAL,DIM,DIM);
    }
#ifndef TWOD
    else if (strcasecmp(token,"box_z")==0) {
      /* 'z' or third vector for box */
      getparam("box_z",&box_z,PARAM_REAL,DIM,DIM);
    }
#endif
#ifndef PS
    else if (strcasecmp(token,"pbc_dirs")==0) {
      /* directions with periodic boundary conditions */
      getparam("pbc_dirs",&pbc_dirs,PARAM_INT,DIM,DIM);
    }
#endif
    else if (strcasecmp(token,"ntypes")==0) {
      /* number of atom types */
      getparam("ntypes",&ntypes,PARAM_INT,1,1);
      vtypes = MAX(vtypes,ntypes);
    }
    else if (strcasecmp(token,"total_types")==0) {
      /* number of virtual atom types */
      getparam("total_types",&vtypes,PARAM_INT,1,1);
    }
#ifdef CONN
    else if (strcasecmp(token,"r_cut")==0) {
      /* cutoff radii */
      int nn;
      if (use_vtypes) nn = SQR(vtypes);
      else            nn = SQR(ntypes);
      r_cut = (real *) calloc(nn,sizeof(real));
      if (NULL == r_cut) error("cannot allocate r_cut");
      getparam("r_cut",r_cut,PARAM_REAL,nn,nn);
    }
#endif
#ifdef COVALENT
    else if (strcasecmp(token,"neigh_len")==0) {
      /* number of neighbors */
      getparam("neigh_len",&neigh_len,PARAM_INT,1,1);
    }
#endif
#if defined(TERSOFF) || defined(PS)
    else if (strcasecmp(token,"ters_r_cut")==0) {     
      getparam("ters_r_cut",ters_r_cut,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
#endif
#ifdef TERSOFF
    else if (strcasecmp(token,"ters_r0")==0) {     
      getparam("ters_r0",ters_r0,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"ters_a")==0) {     
      getparam("ters_a",ters_a,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"ters_b")==0) {     
      getparam("ters_b",ters_b,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"ters_la")==0) {     
      getparam("ters_la",ters_la,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"ters_mu")==0) {     
      getparam("ters_mu",ters_mu,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"ters_chi")==0) {     
      getparam("ters_chi",ters_chi,PARAM_REAL,ntypes*(ntypes-1)/2,ntypes*(ntypes-1)/2);
    }
    else if (strcasecmp(token,"ters_om")==0) {     
      getparam("ters_om",ters_om,PARAM_REAL,ntypes*(ntypes-1)/2,ntypes*(ntypes-1)/2);
    }
    else if (strcasecmp(token,"ters_ga")==0) {     
      getparam("ters_ga",ters_ga,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"ters_n")==0) {     
      getparam("ters_n",ters_n,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"ters_c")==0) {     
      getparam("ters_c",ters_c,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"ters_d")==0) {     
      getparam("ters_d",ters_d,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"ters_h")==0) {     
      getparam("ters_h",ters_h,PARAM_REAL,ntypes,ntypes);
    }
#elif STIWEB
    else if (strcasecmp(token,"stiweb_a")==0) {     
      getparam("stiweb_a",stiweb_a,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"stiweb_b")==0) {     
      getparam("stiweb_b",stiweb_b,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"stiweb_p")==0) {     
      getparam("stiweb_p",stiweb_p,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"stiweb_q")==0) {     
      getparam("stiweb_q",stiweb_q,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"stiweb_a1")==0) {     
      getparam("stiweb_a1",stiweb_a1,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"stiweb_de")==0) {     
      getparam("stiweb_de",stiweb_de,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"stiweb_a2")==0) {     
      getparam("stiweb_a2",stiweb_a2,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"stiweb_ga")==0) {     
      getparam("stiweb_ga",stiweb_ga,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"stiweb_la")==0) {     
      getparam("stiweb_la",stiweb_la,PARAM_REAL,ntypes*ntypes*(ntypes+1)/2,ntypes*ntypes*(ntypes+1)/2);
    }
#elif KEATING
    else if (strcasecmp(token,"keating_alpha")==0) {     
      getparam("keating_alpha",keating_alpha,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"keating_d")==0) {     
      getparam("keating_d",keating_d,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"keating_r_cut")==0) {     
      getparam("keating_r_cut",keating_r_cut,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"keating_beta")==0) {     
      getparam("keating_beta",keating_beta,PARAM_REAL,ntypes*ntypes*(ntypes+1)/2,ntypes*ntypes*(ntypes+1)/2);
    }
#endif
#ifdef PS
    else if (strcasecmp(token,"r_cut")==0) {
      /* Cutoff parameter */
      getparam("r_cut",&r_cut_vec,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
#endif
  } while (!feof(pf));
  fclose(pf);

#ifdef STRAIN
  if ( box_x.x < r_cell || box_y.y < r_cell 
#ifndef TWOD
       || box_z.z < r_cell 
#endif
       )
    error("Cell size greater than box!");
#endif

} /* getparamfile */


/******************************************************************************
*
*  error -- Complain and abort
*
******************************************************************************/

void error(char *mesg)
{
  printf("Error: %s\n",mesg);
  exit(2);
}


/******************************************************************************
*
*  read_atoms - reads atoms and velocities into the cell-array
*
*  The file format is flat ascii, one atom per line, lines beginning
*  with '#' denote comments. Each line consists of
*
*  number type mass x y [z] vx vy [vz] e [rest]
*
*  where
*
*  number   is an arbitrary integer number assigned to each atom
*  type     is the atom's index to the potenital table
*  mass     is the mass of the atom
*  x,y,z    are the atom's coordinates
*  vx,vy,vz are the atoms's velocities
*  e        is the atom's potential energy
*  rest     is ignored until end of line
*
******************************************************************************/

void read_atoms(str255 infilename)
{
  FILE *infile;
  char buf[512];
  int i, p, s, n;
  vektor pos, v;
  real m, e, c;
#ifdef PS
  real tmp = 0.0;
#endif
  cell *to;
  ivektor cellc;

  /* we first try the old checkpoint name, then the new */
  infile = fopen(infilename,"r");
  if ((NULL == infile) && (restart != -1)) {
    infilename = strcat(infilename,".chkpt");
    infile = fopen(infilename,"r");
  }
  if (NULL==infile) {
    sprintf(error_msg,"Cannot open atoms file %s",infilename);
    error(error_msg);
  }

  natoms=0;
  
  /* Read the input file line by line */
  while (!feof(infile)) {

    buf[0] = (char) NULL;
    fgets(buf,sizeof(buf),infile);
    while ('#'==buf[1]) fgets(buf,sizeof(buf),infile); /* eat comments */

#ifdef SINGLE
#ifdef TWOD
    p = sscanf(buf,"%d %d %f %f %f %f %f %f %f",
	       &n,&s,&m,&pos.x,&pos.y,&v.x,&v.y,&e,&c);
#else
    p = sscanf(buf,"%d %d %f %f %f %f %f %f %f %f %f",
	       &n,&s,&m,&pos.x,&pos.y,&pos.z,&v.x, &v.y, &v.z,&e,&c);
#endif
#else
#ifdef TWOD
    p = sscanf(buf,"%d %d %lf %lf %lf %lf %lf %lf %lf",
	       &n,&s,&m,&pos.x,&pos.y,&v.x,&v.y,&e,&c);
#else
    p = sscanf(buf,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	       &n,&s,&m,&pos.x,&pos.y,&pos.z,&v.x, &v.y, &v.z,&e,&c);
#endif
#endif

    if (p>0) {

#ifdef PS
      /* Determine bounding box before rotation */
      real_minl = minvektor(real_minl, pos);
      real_maxl = maxvektor(real_maxl, pos);

      /* Rotations of atom positions */
#ifndef TWOD
      if ( angx != 0.0 ) 
	pos = xrotate(pos, angx);
      if ( angy != 0.0 ) 
	pos = yrotate(pos, angy);
#endif
      if ( angz != 0.0 ) 
	pos = zrotate(pos, angz);
#endif

      /* compute target cell */
      cellc = cell_coord(pos);
      to = PTR_VV(cell_array,cellc,cell_dim);
      /* enlarge it if necessary */
      if (to->n >= to->n_max) alloc_cell(to,to->n_max+CSTEP);
      /* put the data */
      to->ort   [to->n] = pos;
#if (!defined(STRAIN) && !defined(STRESS))
      to->sorte [to->n] = MOD(s,ntypes);
#endif
#if defined(CONN) || defined(ELCO) || defined(COORD)
      to->nummer[to->n] = n;
#endif
#ifdef CONN
      n_min = MIN(n_min,n);
      n_max = MAX(n_max,n);
#endif
#ifdef COORD
      /* Initialization */
      for ( i=0; i<ntypes; i++) 
	to->coord[to->n * ntypes + i] = 0.0;
#endif
#ifdef RING
      to->del[to->n] = 0;
#endif

#ifdef PS
      /* Determine bounding box after rotation */
      minl = minvektor( minl, pos);
      maxl = maxvektor( maxl, pos);

#ifndef TWOD
      if ( p > 6 && ( eng != 0 || kineng != 0 || color_enc != 0 ) ) {
#else
      if ( p > 5 && ( eng != 0 || kineng != 0 || color_enc != 0 ) ) {
#endif
	if ( eng != 0 )
	  tmp = e;
	else if ( kineng != 0 ) 
	  tmp = 0.5 * m * SPROD(v,v);
	else if ( color_enc != 0 )
	  tmp = c;

	maxeng = MAX(maxeng, tmp);
	mineng = MIN(mineng, tmp);

	to->enc[to->n] = tmp;
      }
#endif

      to->n++;
      natoms++;
    }
  }
  fclose(infile);  
}


/******************************************************************************
*
*  loop over all cell pairs
*
******************************************************************************/

void do_work(void (*do_cell_pair)(cell *p, cell *q, vektor pbc))
{
  cell *p,*q;
  int i,j,k;
  int l,m,n;
  int r,s,t;
  vektor pbc;

  /* for each cell */
  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
#ifndef TWOD
      for (k=0; k < cell_dim.z; ++k)
#endif
      {
#ifdef TWOD
        p = PTR_2D_V(cell_array,i,j  ,cell_dim);
#else
        p = PTR_3D_V(cell_array,i,j,k,cell_dim);
#endif
        /* For half of the neighbours of this cell */
        for (l=0; l <= 1; ++l)
          for (m=-l; m <= 1; ++m)
#ifndef TWOD
            for (n=(l==0 ? -m  : -l ); n <= 1; ++n)
#endif
            {
              /* Calculate Indicies of Neighbour */
              r = i+l;  pbc.x = 0;
              s = j+m;  pbc.y = 0;
#ifndef TWOD
              t = k+n;  pbc.z = 0;
#endif

              /* deal with periodic boundary conditions if necessary */
              if (r<0) {
                if (pbc_dirs.x==1) {
                  r = cell_dim.x-1; 
                  pbc.x -= box_x.x;      
                  pbc.y -= box_x.y;
#ifndef TWOD
                  pbc.z -= box_x.z;
#endif
                } else continue;
              }
              if (s<0) {
                if (pbc_dirs.y==1) {
                  s = cell_dim.y-1;
                  pbc.x -= box_y.x;      
                  pbc.y -= box_y.y;
#ifndef TWOD
                  pbc.z -= box_y.z;
#endif
                } else continue;
              }
#ifndef TWOD
              if (t<0) {
                if (pbc_dirs.z==1) {
                  t = cell_dim.z-1;
                  pbc.x -= box_z.x;      
                  pbc.y -= box_z.y;
                  pbc.z -= box_z.z;
                } else continue;
              }
#endif
              if (r>cell_dim.x-1) {
                if (pbc_dirs.x==1) {
                  r = 0; 
                  pbc.x += box_x.x;      
                  pbc.y += box_x.y;
#ifndef TWOD
                  pbc.z += box_x.z;
#endif
                } else continue;
              }
              if (s>cell_dim.y-1) {
                if (pbc_dirs.y==1) {
                  s = 0; 
                  pbc.x += box_y.x;      
                  pbc.y += box_y.y;
#ifndef TWOD
                  pbc.z += box_y.z;
#endif
                } else continue;
              }
#ifndef TWOD
              if (t>cell_dim.z-1) {
                if (pbc_dirs.z==1) {
                  t = 0; 
                  pbc.x += box_z.x;      
                  pbc.y += box_z.y;
                  pbc.z += box_z.z;
                } else continue;
              }
#endif

	      /* Neighbour cell (note that p==q ist possible) */
#ifdef TWOD
              q = PTR_2D_V(cell_array,r,s,cell_dim);
#else
              q = PTR_3D_V(cell_array,r,s,t,cell_dim);
#endif
              /* Do the work */
              do_cell_pair(p,q,pbc);
            }
      }
}




