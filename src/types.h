
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
* types.h -- Data types for IMD
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/


/* double precision is the default */
#ifdef DOUBLE
typedef double real;
#define REAL MPI_DOUBLE
#else 
typedef float real;
#define REAL MPI_FLOAT
#endif

/* Crays use 64bit ints. Thats too much just to enumerate the atoms */
#if defined(CRAY) || defined(t3e)
typedef short int shortint;
typedef short int integer;
#define SHORT   MPI_SHORT
#define INTEGER MPI_SHORT
/* on alphas, all data should be 32 bit aligned */
#elif defined(alpha)
typedef       int shortint;
typedef       int integer;
#define SHORT   MPI_INT
#define INTEGER MPI_INT
#else
typedef short int shortint;
typedef       int integer;
#define SHORT   MPI_SHORT
#define INTEGER MPI_INT
#endif

/* 2d Vector real */
typedef struct {real x; real y; } vektor2d;
/* 2d Vector integer */
typedef struct {int x; int y;   } ivektor2d;

/* 3d Vector real */
typedef struct {real x; real y; real z; } vektor3d;
/* 3d Vector integer */
typedef struct {int x; int y; int z; } ivektor3d;

/* 4d Vector real for temp. input */
typedef struct {real x; real y; real z; real z2; } vektor4d;
typedef struct {int x; int y; int z; int z2; } ivektor4d; 

#ifdef TWOD
typedef struct {real xx, yy, xy; } sym_tensor;
#else
typedef struct {real xx, yy, zz, yz, zx, xy; } sym_tensor;
#endif

/*

Dirty trick to share as much code as possible between 2D and 3D
incarnations of IMD.

The generic 'vektor' type is either 2D or 3D. The vektorNd types are
specific about the dimension. I use the vektorNd version whenever
possible to avoid confusion.

*/

#ifdef TWOD
typedef vektor2d    vektor;
typedef ivektor2d  ivektor;
#else
typedef vektor3d    vektor;
typedef ivektor3d  ivektor;
#endif

#ifdef COVALENT
/* per particle neighbor table for COVALENT */
typedef struct {
    real        *dist;
    shortint    *typ;
    void        **cl;
    integer     *num;
    int         n;
    int         n_max;
} neightab;

typedef neightab* neighptr;
#endif

/* Basic Data Type - The Cell */

typedef struct {
  real        *ort;
#ifndef MONOLJ
  integer     *nummer;   
  shortint    *sorte;
  real        *masse;
  real        *pot_eng;
#endif
#ifdef EAM2                 /* EAM2: variable for the host electron density */
  real        *eam2_rho_h;
#endif
#ifdef DISLOC
  real        *Epot_ref;
  real        *ort_ref;
#endif
#ifdef AVPOS
  real        *sheet;
  real        *avpos;
  real        *av_epot;
#endif
#ifdef ORDPAR
  shortint    *nbanz;
#endif
#ifdef REFPOS
  real        *refpos;
#endif
#ifdef NVX
  real        *heatcond;
#endif
#ifdef STRESS_TENS
  sym_tensor  *presstens;
#endif
  real        *impuls;
  real        *kraft;
#ifdef COVALENT
  neightab    **neigh;
#endif
#ifdef UNIAX
  real        *traeg_moment;
  real        *achse;
  real        *shape;
  real        *pot_well;
  real        *dreh_impuls;
  real        *dreh_moment;
#endif
  int         n;
  int         n_max;
} cell;

typedef struct {
  integer np, nq;
  signed char ipbc[4];
} pair;

typedef cell* cellptr;

/* Buffer for messages */

typedef struct { real    *data;
		 int         n;
		 int     n_max;
	       } msgbuf;


/* String used for Filenames etc. */
typedef char str255[255];

/* data structure to store a potential table or a function table */ 
typedef struct {
  real *begin;      /* first value in the table */
  real *end;        /* last value in the table (followed by extra zeros) */
  real *step;       /* table increment */
  real *invstep;    /* inverse of increment */
  int  maxsteps;    /* physical length of the table */
  real *table;      /* the actual data */
} pot_table_t;

/* data structure for histograms */
typedef struct {
#ifdef STRESS_TENS
#ifdef SHOCK
  float *press_histxx, minpxx, maxpxx;
  float *press_histyy, minpyy, maxpyy;
  float *kin_histxx,   minkxx, maxkxx;
  float *kin_histxxu,  minkxxu,maxkxxu;
  float *kin_histyy,   minkyy, maxkyy;
#ifndef TWOD
  float *press_histzz, minpzz, maxpzz;
  float *kin_histzz,   minkzz, maxkzz;
#endif
#else /* not SHOCK */
  float *press_histxx, minpxx, maxpxx;
  float *press_histyy, minpyy, maxpyy;
  float *press_histxy, minpxy, maxpxy;
#ifndef TWOD
  float *press_histzz, minpzz, maxpzz;
  float *press_histzx, minpzx, maxpzx;
  float *press_histyz, minpyz, maxpyz;
#endif
#endif /* SHOCK */
#endif /* STRESS_TENS */
  float *kin_hist,     minkin, maxkin;
  float *pot_hist,     minpot, maxpot;
  integer *num_hist;
  int     size;
  ivektor dim;
  vektor  ll, ur;
  real    binvol;
} hist_t;

/* data structure for timers */
typedef struct {
#ifdef MPI                /* with MPI_Wtime */
  double start;           /* time when timer was started */
#elif defined(USE_RUSAGE) /* with getrusage */ 
  struct rusage start;    /* time when timer was started */
#else                     /* with times */
  struct tms start;       /* time when timer was started */
#endif
  double total;           /* accumulation of (stop_time - start_time) */
} imd_timer;
