/******************************************************************************
*
* types.h -- Data types for the imd Package
*
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/


/* Single precision is the default */
#ifdef DOUBLE
typedef double real;
#define MPI_REAL MPI_DOUBLE
#else 
typedef float real;
#define MPI_REAL MPI_FLOAT
#endif

/* Crays use 64bit ints. Thats too much just to enumerate the atoms */
#if defined(CRAY) || defined(t3e)
typedef short int shortint;
typedef short int integer;
#define SHORT   MPI_SHORT
#define INTEGER MPI_SHORT
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

#ifdef TTBP
/* per particle neighbor table for TTBP */
typedef struct {
    real        *dist;
    shortint    *typ;
    int         n;
    int         n_max;
} neightab;

typedef neightab* neighptr;
#endif

/* Basic Data Type - The Cell */
typedef struct {
    real        *ort;
#ifndef MONOLJ
    shortint    *sorte;
    real        *masse;
    real        *pot_eng;
#ifdef DISLOC
    real        *Epot_ref;
    real        *ort_ref;
#endif
#ifdef REFPOS
    real        *refpos;
#endif
#ifdef TRANSPORT
    real        *heatcond;
#endif
#ifdef STRESS_TENS
    real        *presstens;
    real        *presstens_offdia;
#endif
    integer     *nummer;
#endif
    real        *impuls;
    real        *kraft;
#ifdef TTBP
    neightab    **neigh;
#endif
    int         n;
    int         n_max;
} cell;

typedef struct {
  integer np, nq;
  signed char ipbc[4];
} pair;

/* Buffer for messages */

typedef struct { real    *data;
		 int         n;
		 int     n_max;
	       } msgbuf;

/* String used for Filenames etc. */
typedef char str255[255];






