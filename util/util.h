
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

/* Dynamically allocated 2D arrray -- sort of */
#define PTR_2D(var,i,j,dim_j) \
  (((var) + ((i)*(dim_j)) + (j)))

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
/* Conversion from eV/A^3 to GPa */
#define CONV 160.21767
#endif

#ifdef PAIR_POT
#define PSTEP 50
#endif

#ifdef EAM
#define EAM_RHO(cell,i) ((cell)->eam_rho_h[i])
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
#include "param.h"

/*****************************************************************************
*
*  typedefs
*
*****************************************************************************/

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
  real    *poteng;
#endif
#ifdef COVALENT
  neightab *neightab_array;
#endif
#ifdef ELCO
  tensor   *stress;
  tensor6d *elco;
  real     *vol;
#ifdef EAM
  tensor   *eam_stress;
  real     *eam_press;
  real     *eam_bulkm;
  real     *eam_dbulkm;
#endif
#endif
#ifdef EAM
  real     *eam_rho_h;
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

#ifdef PAIR_POT
typedef struct {
  real *begin;      /* first value in the table */
  real *end;        /* last value in the table (followed by extra zeros) */
  real *step;       /* table increment */
  real *invstep;    /* inverse of increment */
  int  maxsteps;    /* physical length of the table */
  real *table;      /* the actual data */
} pot_table_t;
#endif

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

#ifdef CONN
void make_numbers(void);
#endif

#ifdef COORD
void init_coord(void);
#endif

#ifdef COVALENT
void do_neighbour_tables(cell *p, cell *q, vektor pbc);
#endif

#ifdef EAM
void init_eam(void);
#endif

#ifdef ELCO
void init_elco(void);
#ifdef PAIR_POT
void do_elco_pair(cell *p, cell *q, vektor pbc);
#endif
#ifdef TERSOFF
void do_elco_tersoff(void);
#elif STIWEB
void do_elco_stiweb(void);
#elif KEATING
void do_elco_keating(void);
#endif
#ifdef EAM
void do_elco_eam(cell *p, cell *q, vektor pbc);
#endif
void write_stress(void);
void write_elco(void);
#endif

#ifdef PAIR_POT
void init_pair(void);
void read_pot_table( pot_table_t *pt, char *filename, int ncols );
void read_pot_table1(pot_table_t *pt, int ncols, char *filename, FILE *infile);
void read_pot_table2(pot_table_t *pt, int ncols, char *filename, FILE *infile);
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

#ifdef TERSOFF
void init_tersoff(void);
#elif STIWEB
void init_stiweb(void);
#elif KEATING
void init_keating(void);
#endif

#ifdef TORSION
void do_angle(cell *p, cell *q, cell *r, cell *s, 
              vektor pbc_q, vektor pbc_r, vektor pbc_s);
void calc_angles(void);
#endif

/*****************************************************************************
*
*  global variables
*
*****************************************************************************/

#ifdef MAIN 
#define EXTERN 
#define INIT(data) = data

#define nulltensor6d  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
#define nulltensor    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
#define nullivektor2d {0,0}
#define  nullvektor2d {0.0,0.0}
#define nullivektor   {0,0,0}
#define einsivektor2d {1,1}
#define einsivektor   {1,1,1}
#define minvektor2d   {Min,Min}
#define minvektor3d   {Min,Min,Min}
#define maxvektor2d   {Max,Max}
#define maxvektor3d   {Max,Max,Max}
#define nullmaxvektor {0.0,0.0,Max}
#define figures0to8   {0,1,2,3,4,5,6,7,8}
#define figures1to8   {1,2,3,4,5,6,7,8}
#define unitvektor2d  {{1,0},{0.1}}
#define unitvektor3d  {{1,0,0},{0,1,0},{0,0,1}}

#else
#define EXTERN extern
#define INIT(data)
#endif

EXTERN cell *cell_array;       /* Array of Cells */
EXTERN ivektor cell_dim;       /* Dimension of above */
EXTERN int natoms;             /* Total number of atoms */
EXTERN int ntypes INIT(0);     /* Total number of different atom types */
EXTERN int vtypes INIT(0);     /* Number of virtual types */
EXTERN int use_vtypes INIT(0); /* flag for using virtual types */
EXTERN str255 progname;        /* Name of current executable argv[0] */
EXTERN int curline;          /* Number of current line for parameter reading */
EXTERN str255 error_msg;       /* string for error message */
/* The simulation box and its inverse */
EXTERN vektor box_x, tbox_x;
EXTERN vektor box_y, tbox_y;
#ifndef TWOD
EXTERN vektor box_z, tbox_z;
#endif
EXTERN vektor height;

#ifdef TWOD
#ifdef PS
EXTERN ivektor pbc_dirs INIT(nullivektor2d);
#else 
EXTERN ivektor pbc_dirs INIT(einsivektor2d);  /* directions with pbc - default is PBC */
#endif
#else
#ifdef PS
EXTERN ivektor pbc_dirs INIT(nullivektor);
#else
EXTERN ivektor pbc_dirs INIT(einsivektor);  /* directions with pbc - default is PBC */
#endif
#endif

/* Filenames */
EXTERN str255 infilename;    /* Input File */
EXTERN str255 outfilename;   /* Output File */
EXTERN char *paramfilename;  /* parameter file */
/* Bookkeeping */
EXTERN int  restart INIT(-1);
EXTERN int  avpos   INIT(-1);
EXTERN real r2_cut;
#if defined(ANGLE) || defined(COORD)
EXTERN real r_min INIT(0.0);     /* default value */
#else
EXTERN real r_min INIT(1.0);     /* default value */
#endif
#ifdef PAIR
EXTERN real r_max INIT(5.0);     /* default value */
#else
EXTERN real r_max INIT(1.0);     /* default value */
#endif
#if defined(ANGLE) || defined(TORSION)
EXTERN int  nangles INIT(0);
#endif

#ifdef CONN
EXTERN real *r_cut;
EXTERN int  n_min INIT(0);
EXTERN int  n_max INIT(0);
#endif

#ifdef STRAIN
real r_cell INIT(1.0);    /* default value */
#endif

#if defined(STRESS) || defined(ELCO)
EXTERN int  atomcount INIT(0);   /* statistics */ 
EXTERN int  maxneigh  INIT(0);
EXTERN int  sumneigh  INIT(0);
EXTERN int  maxvert   INIT(0);
EXTERN int  sumvert   INIT(0);
EXTERN int  maxedges  INIT(0);
EXTERN int  sumedges  INIT(0);
#ifndef TWOD
EXTERN int maxfaces  INIT(0);
EXTERN int sumfaces  INIT(0);
#endif
#endif

/* The histograms */
#ifndef RING
EXTERN real      *histogram;
#endif
#if defined(ANGLE) || defined(PAIR) || defined(TORSION)
EXTERN int slots INIT(1000);
#endif
#ifdef ANGLE
EXTERN ivektor4d hist_dim;
#endif
#ifdef PAIR
EXTERN ivektor3d hist_dim;
#endif
#ifdef TORSION
EXTERN ivektor5d hist_dim;
#endif

#if defined(STRESS) || defined(ELCO)
EXTERN int     neighnum;        /* Number of neighbour atoms */
EXTERN vektor  *candcoord;      /* Coordinates of neighbour atoms */
EXTERN real    *canddist2;      /* Distance squared of neighbour atoms */ 
#ifndef TWOD
EXTERN real    volume;          /* Volume of Voronoi cell, 3d */
#else
EXTERN real    area;            /* Area of Voronoi cell, 2d */
#endif
#endif

#ifdef COORD
EXTERN real      *numbers;
EXTERN int       use_unity INIT(-1);
EXTERN int       local     INIT(-1);
EXTERN int       global    INIT(-1);
EXTERN int       c_max     INIT(10);
EXTERN ivektor2d num_dim;
EXTERN int       bonds     INIT(0); /* Total number of bonds */
EXTERN real      r_cut_vec[55] INIT({-1.0});
EXTERN real      r_cut[10][10];
EXTERN int       write_poteng INIT(0);
#endif

#ifdef CONN
/* Connection matrix */
EXTERN ivektor2d cm_dim ;
EXTERN int maxneigh INIT(10);
EXTERN int *cm, *nn, *tp, *num, *ind;
EXTERN vektor *pos;

#endif

#ifdef ELCO
EXTERN tensor6d c INIT(nulltensor6d);
EXTERN tensor sigma INIT(nulltensor);
EXTERN real   r_cell INIT(-1.0);
EXTERN real   vol;
EXTERN real   epot INIT(0.0); 
EXTERN real   press INIT(0.0); 
EXTERN real   bulkm INIT(0.0);
EXTERN real   dbulkm_dp INIT(0.0);
EXTERN int    stresstens INIT(0); 
EXTERN int    moduli INIT(0); 
EXTERN int    all_moduli INIT(0);
#endif

#ifdef COVALENT
EXTERN int  neigh_len INIT(6);
EXTERN int  nbonds    INIT(0);
#endif
#ifdef TERSOFF
EXTERN real ter_r_cut[10][10], ter_r2_cut[10][10], ters_r_cut[55] ;
EXTERN real ter_r0[10][10], ters_r0[55];
EXTERN real ter_a[10][10], ters_a[55];
EXTERN real ter_b[10][10], ters_b[55];
EXTERN real ter_la[10][10], ters_la[55];
EXTERN real ter_mu[10][10], ters_mu[55];
EXTERN real ters_ga[10];
EXTERN real ters_n[10];
EXTERN real ter_c2[10], ters_c[10];
EXTERN real ter_d2[10], ters_d[10];
EXTERN real ters_h[10];
EXTERN real ter_chi[10][10], ters_chi[45];
EXTERN real ter_om[10][10], ters_om[45];
#elif STIWEB
EXTERN real stiweb_a[55];
EXTERN real sw_a[10][10];
EXTERN real stiweb_b[55];
EXTERN real sw_b[10][10];
EXTERN real stiweb_p[55];
EXTERN real sw_p[10][10];
EXTERN real stiweb_q[55];
EXTERN real sw_q[10][10];
EXTERN real stiweb_a1[55];
EXTERN real sw_a1[10][10];
EXTERN real stiweb_de[55];
EXTERN real sw_de[10][10];
EXTERN real sw_2_a1[10][10];
EXTERN real stiweb_a2[55];
EXTERN real sw_a2[10][10];
EXTERN real stiweb_la[550];
EXTERN real sw_la[10][10][10];
EXTERN real stiweb_ga[55];
EXTERN real sw_ga[10][10];
EXTERN real sw_r_cut[10][10];
#elif KEATING
EXTERN real keating_alpha[55];
EXTERN real keat_alpha[10][10];
EXTERN real keating_d[55];
EXTERN real keat_d[10][10];
EXTERN real keating_r_cut[55];
EXTERN real keat_r_cut[10][10];
EXTERN real keat_r2_cut[10][10];
EXTERN real keating_beta[550];
EXTERN real keat_beta[10][10][10];
#endif

#ifdef PAIR_POT
EXTERN pot_table_t pair_pot;         /* potential data structure */
EXTERN str255 potfilename;           /* Potential */
#endif

#ifdef EAM
EXTERN pot_table_t embed_pot;         /* embedding energy table  */
EXTERN pot_table_t rho_h_tab;         /* electron transfer table */
EXTERN str255 eam_emb_E_filename;     /* embedding energy file   */
EXTERN str255 eam_at_rho_filename;    /* electron transfer file  */
#endif

#ifdef RING
EXTERN int *histogram;
EXTERN int total_rings INIT(0); /* Number of sp rings found */
EXTERN int first       INIT(1); /* =1: First generation of neighbour tables */
EXTERN int max_length  INIT(10); /* Maximum length of paths, default 10 */
EXTERN Queue_elmt *queue;
EXTERN int queue_length;
EXTERN Queue_elmt *sp_queue;
EXTERN int sp_queue_length;
EXTERN atom *stack;
EXTERN int stack_end;
#endif
#ifdef PS
EXTERN List_elmt *atomlist;
EXTERN int    natoms_drawn INIT(0);
EXTERN int    nbonds_drawn INIT(0);
EXTERN real   dist_min INIT(Min);
EXTERN real   dist_max INIT(Max);
EXTERN real   r_atom   INIT(-1.0); 
EXTERN real   r_bond   INIT(-1.0);
EXTERN long   atomcolor INIT(1);
EXTERN real   atombrightness INIT(0.0);
EXTERN real   bondbrightness INIT(-0.3);
EXTERN long   bondcolor INIT(0);
EXTERN real   ambient INIT(1.0);
EXTERN long   radii   INIT(0);
EXTERN long   eng     INIT(0);
EXTERN long   kineng  INIT(0);
EXTERN long   color_enc INIT(0);
EXTERN int    colorencoding INIT(0);
EXTERN real   linew  INIT(-1.0);
EXTERN real   alinew INIT(1.0);
EXTERN real   blinew INIT(1.0);
EXTERN real   slinew INIT(1.0);
EXTERN real   llinew INIT(-1.0);
EXTERN real   shading INIT(45.0);
EXTERN real   scaling INIT(1.0);
EXTERN real   frame INIT(0.0);
EXTERN real   frame_linew INIT(2.0);
EXTERN real   rframe INIT(-1.0);
EXTERN real   rframe_linew INIT(1.0);
EXTERN real   sframe INIT(-1.0);
EXTERN real   sframe_linew INIT(1.0);
EXTERN real   frame_border INIT(9.0);
EXTERN real   frame_border_linew INIT(0.0);
EXTERN vektor2d translation INIT(nullvektor2d);
EXTERN real   zoom INIT(1.0);
EXTERN real   ps_width INIT(-1.0);
EXTERN real   ps_height INIT(-1.0);
EXTERN int    axes INIT(-1);
EXTERN real   spect INIT(-1.0);
EXTERN real   backgrd INIT(0.0);
EXTERN real   needle INIT(-1.0);
EXTERN real   depth  INIT(0.0);
EXTERN real   user_scale INIT(1.0);
EXTERN real   wireframe INIT(-1.0);
EXTERN real   crop INIT(-1.0);
EXTERN int    settings INIT(-1);
EXTERN vektor ursprung;
EXTERN vektor cunit;
#ifndef TWOD
EXTERN vektor maxl INIT(maxvektor3d); 
EXTERN vektor minl INIT(minvektor3d); 
EXTERN vektor real_maxl INIT(maxvektor3d); 
EXTERN vektor real_minl INIT(minvektor3d);
#else
EXTERN vektor maxl INIT(maxvektor2d); 
EXTERN vektor minl INIT(minvektor2d); 
EXTERN vektor real_maxl INIT(maxvektor2d); 
EXTERN vektor real_minl INIT(minvektor2d);
#endif
EXTERN real   maxeng INIT(Max); 
EXTERN real   mineng INIT(Min);
EXTERN real   angx INIT(0.0); 
EXTERN real   angy INIT(0.0); 
EXTERN real   angz INIT(0.0);
#ifndef TWOD
EXTERN real   zeta INIT(Max);
EXTERN vektor foc  INIT(nullmaxvektor);
#endif
EXTERN real   r_cut_vec[55] INIT({-1.0}); 
EXTERN real   ters_r_cut[55] INIT({ -1.0});
EXTERN real   r_cut[10][10];
EXTERN real   r2cut[10][10];
EXTERN vektor3d color_array[9];
EXTERN vektor3d enc_array[8];
EXTERN real   shade_array[9];
EXTERN real   radii_array[9];
EXTERN int    atomcol[9] INIT(figures0to8);
EXTERN int    bondcol[9] INIT(figures0to8);
EXTERN int    rad[9];
EXTERN int    enc[8] INIT(figures1to8);
EXTERN int    enc_len;
#ifndef TWOD
EXTERN vektor unitv[3] INIT(unitvektor3d);
#else
EXTERN vektor unitv[2] INIT(unitvektor2d);
#endif
#endif

#ifdef MAIN

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

#ifdef COVALENT
  neightab *neigh;
  int i;
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
    cl->poteng        = NULL;
#endif
#ifdef COVALENT
    cl->neightab_array= NULL;
#endif
#ifdef ELCO
    cl->stress        = NULL;
    cl->elco          = NULL;
    cl->vol           = NULL;
#ifdef EAM
    cl->eam_stress    = NULL;
    cl->eam_press     = NULL;
    cl->eam_bulkm     = NULL;
    cl->eam_dbulkm    = NULL;
#endif
#endif
#ifdef EAM
    cl->eam_rho_h     = NULL;
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
  cl->poteng = (real   *) realloc(cl->poteng, ntypes * count * sizeof(real));
#endif
#ifdef ELCO
  cl->stress = (tensor *) realloc(cl->stress, count * sizeof(tensor));
  cl->elco   = (tensor6d *) realloc(cl->elco, count * sizeof(tensor6d));
  cl->vol    = (real   *) realloc(cl->vol,    count * sizeof(real));
#ifdef EAM
  cl->eam_stress = (tensor *) realloc(cl->eam_stress, count * sizeof(tensor));
  cl->eam_press  = (real *)   realloc(cl->eam_press, count * sizeof(real));
  cl->eam_bulkm  = (real *)   realloc(cl->eam_bulkm, count * sizeof(real));
  cl->eam_dbulkm = (real *)   realloc(cl->eam_dbulkm, count * sizeof(real));
#endif
#endif
#ifdef EAM
  cl->eam_rho_h = (real *) realloc(cl->eam_rho_h, count * sizeof(real));
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
    || (NULL==cl->poteng)
#endif
#ifdef ELCO
    || (NULL==cl->stress)
    || (NULL==cl->elco)
    || (NULL==cl->vol)
#ifdef EAM
    || (NULL==cl->eam_stress)
    || (NULL==cl->eam_press)
    || (NULL==cl->eam_bulkm)
    || (NULL==cl->eam_dbulkm)
#endif
#endif
#ifdef EAM
    || (NULL==cl->eam_rho_h)
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
#if defined(ANGLE) || defined(PAIR) || defined(COORD)
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
#if defined(ANGLE) || defined(PAIR) || defined(TORSION)
      /* n - number of slots in histogram */
    case 'n':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          slots = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else slots = atof(&argv[1][2]);
      break;
#endif
#ifdef CONN
      /* n - number of maximal neighbours */
    case 'n':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          maxneigh = atof(argv[2]);
          --argc;
          ++argv;
        }
      }
      else maxneigh = atof(&argv[1][2]);
      cm_dim.x = 0;
      cm_dim.y = maxneigh;
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
     /* E - write potential energy */
    case 'E':
      write_poteng = 1;
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

  pf = fopen(paramfname,"r");
  if (NULL == pf) {
    sprintf(error_msg,"Cannot open parameter file %s",paramfname);
    error(error_msg);
  }

  do {
    res=fgets(buffer,1024,pf);
    if (NULL == res) break; /* probably EOF reached */
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
#ifdef PAIR_POT
    else if (strcasecmp(token,"potfile")==0) {
      getparam("potfile",potfilename,PARAM_STR,1,255);
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
#ifdef EAM
    else if (strcasecmp(token,"core_potential_file")==0) {
      /* EAM:Filename for the tabulated Core-Core Potential (r^2) */
      getparam("core_potential_file",potfilename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"embedding_energy_file")==0) {
      /* EAM:Filename for the tabulated Embedding Enery(rho_h) */
      getparam("embedding_energy_file",eam_emb_E_filename,PARAM_STR,1,255);
    }
   else if (strcasecmp(token,"atomic_e-density_file")==0) {
      /* EAM:Filename for the tabulated atomic electron density(r_ij^2) */
      getparam("atomic_e-density_file",eam_at_rho_filename,PARAM_STR,1,255);
    }
#endif
#if defined(COORD) || defined(PS)
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
  int p, s, n;
  vektor pos, v;
  real m, e, c;
#ifdef COORD
  int i;
#endif
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

#ifdef TWOD
    p = sscanf(buf,"%d %d %lf %lf %lf %lf %lf %lf %lf",
	       &n,&s,&m,&pos.x,&pos.y,&v.x,&v.y,&e,&c);
#else
    p = sscanf(buf,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	       &n,&s,&m,&pos.x,&pos.y,&pos.z,&v.x, &v.y, &v.z,&e,&c);
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
      /* Read potential energy */
      to->poteng[to->n] = e;
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

#endif /* MAIN */
