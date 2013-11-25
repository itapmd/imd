
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2011 Institute for Theoretical and Applied Physics,
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
/* on NEC SX, shorts are not vectorizable */
#elif defined(ALPHA) || defined(SX)
typedef int shortint;
typedef int integer;
#define SHORT   MPI_INT
#define INTEGER MPI_INT
#else
typedef int shortint;
typedef int integer;
#define SHORT   MPI_INT
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

typedef struct {real yz, zx, xy; } off_tensor;

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

#if defined(COVALENT) || defined(NNBR_TABLE)
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

#ifdef NYETENSOR
typedef struct {
	real lcm[3][3];		  /* Latticae corresponding matrix */
	real nyeTensor[3][3]; /* Nye tensor */
	real bv[3];			  /* Burgers vector */
	real ls[3];           /* Line sense / Line direction of dislocation */
} nyeTensorInfo;
typedef struct {
	real n[14][3];		  /* Reference nearest neighbors in perfect crystal */
} neighPerfType;
#endif

/* data structure for KIM */
#ifdef KIM
typedef struct {
  /* pointer to the KIM_API object */
  void *pkim;

  /* values set in init_kim_info() */
  int model_using_half;
  int model_using_cluster;
  int model_using_Rij;

  /* values set in model_has_flags(), called by kim_init() */
  int model_has_forces;
  int model_has_energy;
  int model_has_particleEnergy;
  int model_has_process_dEdr;

  /* values set in kim_init(), after call to string_init(_) */
  int ind_coordinates;
  int ind_numberOfParticles;
  int ind_numberContributingParticles;
  int ind_numberParticleTypes;
  int ind_particleTypes;
  int ind_get_neigh;
  int ind_neighObject;
  int ind_cutoff;
  int ind_energy;
  int ind_particleEnergy;
  int ind_forces;
  int ind_process_dEdr;

  /* index of the current cell */
  int cell_ind;
  int *cell_list;
  int *cell_offset;
  int cell_atom_ind;
  int *cell_index_atom;

  /* pointer for the particle types mapping */
  int *kim_particle_codes;
  int iterator_position;
} imd_kim_t;
#endif /* KIM */

#ifdef BBOOST
/* per particle neighbor table for BBOOST */
typedef struct {
    real        *distref1;   /* array containg info on bondlength */
    real        *distref2;
    integer     *numref1;    /* array containg the neighbor info */
    integer     *numref2;
    int         nbondsref1;  /* number of bonds this atoms had in reference pos */
    int         nbondsref2;
} bb_neightab;

typedef bb_neightab* bb_neighptr;
#endif



/* Basic Data Type - The Cell */

typedef struct {
    real        *ort;
#ifdef BBOOST
    real        *bb_refposone; /* reference positions for bondboost */
    real        *bb_refpostwo;
    real        *bb_oldpos; /* position before boosting */
    bb_neightab **bb_neigh;
#endif
#ifndef MONOLJ
  integer     *nummer;
  shortint    *sorte;
  shortint    *vsorte;
  real        *masse;
  real        *pot_eng;
#endif
#ifdef EAM2
  real        *eam_rho;     /* host electron density */
  real        *eam_dF;      /* derivative of embedding energy */
#ifdef EEAM
  real        *eeam_p_h;    /* energy modification term */
  real        *eeam_dM;     /* derivative of energy modification */
#endif
#endif
#ifdef ADP
  real        *adp_mu;
  sym_tensor  *adp_lambda;
#endif
#ifdef VARCHG
  real        *charge;      /* individual charge for each particle */
#endif
#ifdef SM
  real *chi_sm;              /* electronegativity */
  real *z_sm;                /* effective core charge */
  real *j_sm;               /* repulsiveness parameter */
  real *v_sm;                /* integral */

  real *b_sm;                /* Ax=b */
  real *x_sm;                /* Ax=b */
  real *r_sm;                /* residuum Ax=b */
  real *d_sm;                /* conjugate directions Ax=b */
  real *s_sm;                /* auxiliary variable Ax=b */
  real *q_sm;                /* initial value */

#endif
#ifdef DIPOLE
  real        *dp_E_stat;    /* electric field at atom location */
  real        *dp_E_ind;     /* induced field at atom location */
  real        *dp_E_old_1;   /* old field from previous steps */
  real        *dp_E_old_2;   /* old field from previous steps */
  real        *dp_E_old_3;   /* old field from previous steps */
  real        *dp_p_stat;    /* static dipoles from Short-Range interaction */
  real        *dp_p_ind;     /* induced dipoles */
#endif /* DIPOLE */
#ifdef CG
  real        *h;           /* Conjugated Gradient: search vektor */
  real        *g;           /* Conjugated Gradient: old forces */
  real        *old_ort;     /* CG: old locations, needed for linmin */
#endif
#ifdef DAMP
  real        *damp_f; /* damping function for that atom, position dependent */
#endif
#ifdef DISLOC
  real        *Epot_ref;
  real        *ort_ref;
#endif
#ifdef CNA
  long        *mark;
#endif
#ifdef AVPOS
  real        *sheet;
  real        *avpos;
  real        *av_epot;
#endif
#ifdef NNBR
  shortint    *nbanz;
#endif
#ifdef REFPOS
  real        *refpos;
#endif
#ifdef HC
  real        *hcaveng;
#endif
#ifdef STRESS_TENS
  sym_tensor  *presstens;
#ifdef AVPOS
  sym_tensor  *avpresstens;
#endif
#endif
#ifdef SHOCK
  real        *pxavg;
#endif
  real        *impuls;
  real        *kraft;
#if defined(COVALENT) || defined(NNBR_TABLE)
  neightab    **neigh;
#endif
#ifdef NBLIST
  real        *nbl_pos;
#endif
#ifdef UNIAX
  real        *achse;
  real        *dreh_impuls;
  real        *dreh_moment;
#endif
#if defined(VEC) && defined(MPI)
  integer     *ind;
#endif
  int         n;
  int         n_max;
#ifdef VEC
  int         n_buf;
#endif
#ifdef TTM
  ivektor fd_cell_idx;
#endif
#ifdef ADA
  char *adaType;
  char *hopsToDefect;
#endif
#ifdef NYETENSOR
  nyeTensorInfo **nyeTens;
#endif

} cell;

typedef cell* cellptr;

#ifdef VEC
typedef struct {
  int *ind;
  int n;
  int n_max;
} minicell;
#else
typedef cell minicell;
#endif

typedef struct {
  char format;
  int  endian;
  int  n_number;
  int  n_type;
  int  n_mass;
  int  n_pos;
  int  n_vel;
  int  n_data;
  int  n_items;
#ifdef REFPOS
  int  n_refpos_x;
#endif
#ifdef DISLOC
  int  n_x_ref;
  int  n_Epot_ref;
#endif
#ifdef VARCHG
  int  n_charge;
#endif
} header_info_t;

typedef union {
  integer i;
  float   f;
} i_or_f;

typedef union {
#ifdef DOUBLE
  integer i[2];
#else
  integer i;
#endif
  real    r;
} i_or_r;

typedef union {
  integer i[2];
  double  d;
} i_or_d;

typedef struct {
  integer np, nq;
#if ! defined NBLIST || defined BBOOST
  signed char ipbc[4];
#endif
} pair;

#ifdef NBLIST
typedef struct {
  integer np, nq[NNBCELL];
} cell_nbrs_t;
#endif

/* Buffer for messages */
typedef struct {
  real *data;
  int  n;
  int  n_max;
} msgbuf;

/* String used for Filenames etc. */
typedef char str255[255];

#ifdef TTM

#define FE_VACUUM 0
#define FE_MATTER 1

/* structure for FD lattice elements */
typedef struct
{
  int natoms; /* number of atoms in encompassed MD cells */
  cellptr * md_cellptrs; /* array of pointers to MD cells */
  real temp; /* electron temperature */
  real xi; /* damping parameter for coupling to MD system */
  real md_temp; /* avg. temperature of the MD cells */
  real source; /* thermal power to be coupled into electronic system;
  source term for pdeq, needs to be updated every timestep if time-dependent
  source is desired (for example, look at laser_rescale_ttm() in imd_laser.c) */
  vektor3d v_com; /* velocity of the center of mass of MD cells */
} ttm_Element;

#endif /*TTM*/

/* data structure to store a potential table or a function table */
typedef struct {
  real *begin;      /* first value in the table */
  real *end;        /* last value in the table (followed by extra zeros) */
  real *step;       /* table increment */
  real *invstep;    /* inverse of increment */
  int  *len;        /* length of the individual columns */
  int  ncols;       /* number of columns in the table */
  int  maxsteps;    /* physical length of the table */
  real *table;      /* the actual data */
#ifdef SPLINE
  real *table2;     /* second derivatives for spine interpolation */
#endif
} pot_table_t;

#ifdef LINPOT
/* data structure to store a potential table or a function table */
typedef struct {
  real *begin;      /* first value in the table */
  real *end;        /* last value in the table */
  real *step;       /* table increment */
  real *invstep;    /* inverse of increment */
  int  *len;        /* length of the individual columns */
  int  ncols;       /* number of columns in the table */
  real **table;     /* the actual data */
} lin_pot_table_t;
#endif

/* data structure for timers */
typedef struct {
#ifdef MPI                  /* with MPI_Wtime */
#ifdef MPE
  int mpe_id_begin;         /* MPE event id (begin) */
  int mpe_id_end;           /* MPE event id (end) */
  int mpe_flag;             /* whether to do MPE logging */
#endif
  double start;             /* time when timer was started */
#elif defined(USE_WALLTIME) /* with gettimeofday */
  struct timeval start;     /* time when timer was started */
#elif defined(OMP)          /* with omp_get_wtime */
  double start;             /* time when timer was started */
#elif defined(USE_RUSAGE)   /* with getrusage */
  struct rusage start;      /* time when timer was started */
#else                       /* with times */
  struct tms start;         /* time when timer was started */
#endif
  double total;             /* accumulation of (stop_time - start_time) */
} imd_timer;

