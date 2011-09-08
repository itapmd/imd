
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2004 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
*  Global variables, prototyps, makros, typedefs etc.
*  for IMD configuration utility programs
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

/******************************************************************************
*
*  some usful makros 
*
******************************************************************************/

#define  UTIL_H 1

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

#ifdef CNA
#define MAX_NEIGH 12
#define MAX_BONDS 24
#define MAX_TYPES 25
#define PI 3.1415927
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

#if defined(PS) || defined(CNA)
#define NSTEP 1
#endif

#ifdef PS
#define ABS(a) ((a)>=0) ? (a) : (-(a))
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
#include <float.h>
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
#ifdef SLIP
    vektor *dsp;
    vektor *slip;
#endif
#ifdef STRAIN
  vektor *dsp;
  vektor *strain;
  vektor *strain_offdia;
  short  *empty;
#elif defined(STRESS)
  vektor *stress;
  vektor *stress_offdia;
  real   *vol;
  int    *sorte;
  real   *masse;
  int    *nummer;
#else
  int    *sorte;
#endif
#if defined(CONN) || defined(ELCO) || defined(COORD) || defined(CNA) ||defined(REMAT)
  int    *nummer;
#endif
#ifdef REMAT
  real   *masse;
  int *nnn;
  int *flag;
#endif
#ifdef CNA
 short  *cna;
  short  *mark;
  real   *masse;
  real   *ppc;
#endif
#ifdef OVEC
  vektor *ovec;
#endif
#if defined(COORD) || defined(CNA)
     real    *coord;
    //short    *coord;
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
void read_command_line(int argc, char **argv);
void read_arg_bool(int   *argcptr, char ***argvptr, int    *parptr);
void read_arg_int(int    *argcptr, char ***argvptr, int    *parptr);
void read_arg_long(int   *argcptr, char ***argvptr, long   *parptr);
void read_arg_real(int   *argcptr, char ***argvptr, real   *parptr);
void read_arg_string(int *argcptr, char ***argvptr, char  **parptr);
void read_arg_vektor(int *argcptr, char ***argvptr, vektor *parptr);
void read_arg_ivektor4d(int *argcptr, char ***argvptr, ivektor4d *parptr);
void read_parameters(void);
void getparamfile(char *infile);
void read_atoms(str255 infilename);
void usage(void);
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

#ifdef CNA
void init_cna(void);
void do_cna(void);
//void do_cna(cell *p, cell *q, vektor pbc);
void domino(int start, int end, int listlength, int *max_chain, int *chain);
void sort_pair_types(void);
void write_atoms(void);
void write_statistics(void);
void write_rdf(void);
int  atom_in_pbox(vektor pos);
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
void write_elco_select(void);
#endif

#ifdef REMAT
void calc_n_toclose_nn(int n_nnn_crit);
void write_atoms(void);

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
#define nullvektor2d  {0.0,0.0}
#define nullivektor   {0,0,0}
#define nullivektor4d {0,0,0,0}
#define einsivektor2d {1,1}
#define einsivektor   {1,1,1}
#define minvektor2d   {DBL_MAX,DBL_MAX}
#define minvektor3d   {DBL_MAX,DBL_MAX,DBL_MAX}
#define maxvektor2d   {-DBL_MAX,-DBL_MAX}
#define maxvektor3d   {-DBL_MAX,-DBL_MAX,-DBL_MAX}
#ifdef PS
#define nullmaxvektor {0.0,0.0,Max}
#endif
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
EXTERN int curline;            /* Number of current line for parameter reading */
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

EXTERN int box_from_header INIT(1);

/* Filenames */
EXTERN str255 infilename;    /* Input File */
EXTERN str255 outfilename;   /* Output File */
EXTERN char *paramfilename;  /* parameter file */
/* Bookkeeping */
EXTERN int  restart INIT(-1);
EXTERN int  avpos   INIT(-1);
EXTERN real r2_cut;
#if defined(ANGLE) || defined(COORD) || defined(CNA) 
EXTERN real r_min INIT(0.0);     /* default value */
#else
EXTERN real r_min INIT(1.0);     /* default value */
#endif
#ifdef PAIR
EXTERN real r_max INIT(5.0);     /* default value */
#elif CNA
EXTERN real r_max INIT(-1.0);
#else
EXTERN real r_max INIT(1.0);     /* default value */
#endif
#if defined(ANGLE) || defined(TORSION)
EXTERN int  nangles INIT(0);
#endif

#ifdef CNA
EXTERN ivektor4d pair_type INIT(nullivektor4d);
EXTERN int pair_type_short INIT(0);
EXTERN int writeatoms INIT(0);
EXTERN int rdf INIT(0);
EXTERN int nearestneighbour INIT(0);
EXTERN int cna_pairs INIT(0);
EXTERN int bondlist[MAX_BONDS][3];
EXTERN ivektor4d type_list[MAX_TYPES];
EXTERN int type_count[MAX_TYPES];
EXTERN int type_sort[MAX_TYPES];
EXTERN int type_list_length INIT(0);
EXTERN int *rdf_tab[MAX_TYPES];
EXTERN real rcut;
//void do_cna(cell *p, cell *q, vektor pbc);
EXTERN int  cna_stat[5] INIT({0});
#endif

#ifdef CNA
/* use of partial box */
#ifndef TWOD
EXTERN vektor ll INIT(maxvektor3d);
EXTERN vektor ur INIT(minvektor3d);
#else
EXTERN vektor ll INIT(maxvektor2d);
EXTERN vektor ur INIT(minvektor2d);
#endif
EXTERN int    use_ll INIT(0);
EXTERN int    use_ur INIT(0);
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
#if defined(ANGLE) || defined(PAIR) || defined(TORSION) || defined(CNA)
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

#if defined(CNA) || defined(COORD)
EXTERN real      r_cut_vec[55] INIT({-1.0});
EXTERN real      r_cut[10][10];
#endif

#ifdef COORD
EXTERN real      *numbers;
EXTERN int       use_unity INIT(-1);
EXTERN int       local     INIT(-1);
EXTERN int       global    INIT(-1);
EXTERN int       c_max     INIT(10);
EXTERN ivektor2d num_dim;
EXTERN int       bonds     INIT(0); /* Total number of bonds */
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
EXTERN int    cindex INIT(-1);
EXTERN int    select_moduli INIT(0);
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

#ifdef REMAT
EXTERN int   n_nnn[12];
EXTERN int   ndeletedatoms;
EXTERN real  r_crit;

#endif
