
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
* makros.h -- Some useful makros for IMD
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#if defined(__GNUC__)
#define INLINE inline
#else
#define INLINE
#endif

/* avoid p % q, which is terribly slow */
/* on SGI, inline doesn't really work :-( */
#if !defined(MONO) && !defined(BINARY)
#if defined(t3e)
#pragma _CRI inline(MOD)
#elif defined(ALPHA)
#pragma inline(MOD)
#elif defined(sgi)
#pragma inline global (MOD)
#endif
INLINE static int MOD(shortint p, int q)
{
  int stmp=p;
  while (stmp>=q) stmp-=q;
  return stmp;
}
#endif

/* Sometimes we use array where we should use vectors but... */
#define X(i) [DIM*(i)  ]
#define Y(i) [DIM*(i)+1]
#define Z(i) [DIM*(i)+2]

#ifdef MONOLJ
#define SORTE(cell,i)           0
#define VSORTE(cell,i)          0
#define NUMMER(cell,i)          0
#define MASSE(cell,i)           1.0
#define POTENG(cell,i)          0.0
#else
#if defined(MONO)
#define SORTE(cell,i)           0
#elif defined(BINARY)
#define SORTE(cell,i)           (((cell)->sorte[i]) & 0x1)
#else
#define SORTE(cell,i)           MOD((cell)->sorte[i],ntypes)
#endif
#define VSORTE(cell,i)          ((cell)->sorte[i])
#define NUMMER(cell,i)          ((cell)->nummer[i])
#define MASSE(cell,i)           ((cell)->masse[i])
#define POTENG(cell,i)          ((cell)->pot_eng[i])
#endif

#define ORT(cell,i,sub)         ((cell)->ort sub(i))
#define KRAFT(cell,i,sub)       ((cell)->kraft sub(i))
#define IMPULS(cell,i,sub)      ((cell)->impuls sub(i))

#ifdef EAM2
#define EAM_RHO(cell,i)         ((cell)->eam2_rho_h[i])
#endif
#ifdef CG
#define CG_G(cell,i,sub)        ((cell)->g sub(i))
#define CG_H(cell,i,sub)        ((cell)->h sub(i))
#define OLD_ORT(cell,i,sub)     ((cell)->old_ort sub(i))
#endif
#ifdef DISLOC
#define EPOT_REF(cell,i)        ((cell)->Epot_ref[i])
#define ORT_REF(cell,i,sub)     ((cell)->ort_ref sub(i))
#endif
#ifdef AVPOS
#define AV_POS(cell,i,sub)      ((cell)->avpos sub(i))
#define SHEET(cell,i,sub)       ((cell)->sheet sub(i))
#define AV_EPOT(cell,i)         ((cell)->av_epot[i])
#endif
#ifdef ORDPAR
#define NBANZ(cell,i)           (cell)->nbanz[(i)]
#endif
#ifdef REFPOS
#define REF_POS(cell,i,sub)     ((cell)->refpos sub(i))
#endif
#ifdef NVX
#define HEATCOND(cell,i)        ((cell)->heatcond[i])
#endif
#ifdef STRESS_TENS
#define PRESSTENS(cell,i,sub)   ((cell)->presstens[i].sub)
#endif
#ifdef COVALENT
#define NEIGH(cell,i)           ((cell)->neigh[i])
#endif
#ifdef UNIAX
#define TRAEG_MOMENT(cell,i)    ((cell)->traeg_moment[i])
#define ACHSE(cell,i,sub)       ((cell)->achse sub(i))
#define SHAPE(cell,i,sub)       ((cell)->shape sub(i))
#define POT_WELL(cell,i,sub)    ((cell)->pot_well sub(i))
#define DREH_IMPULS(cell,i,sub) ((cell)->dreh_impuls sub(i))
#define DREH_MOMENT(cell,i,sub) ((cell)->dreh_moment sub(i))
#endif

#ifdef MPI
#define CELLS(k) cells[k]
#else
#define CELLS(k) k
#endif

/* Max gibt den groesseren von zwei Werten */
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* Min gibt den kleineren  von zwei Werten */
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/* Sqr quadriert sein Argument */
#if defined(__GNUC__) || defined(__SASC)
inline static real SQR(real x)
{
  return x*x;
}
#else
#define SQR(a) ((a)*(a))
#endif

/* Abs berechnet den Betrag einer Zahl */
#define ABS(a) ((a) >0 ? (a) : -(a))

/* How many dimension are there? */
#ifdef TWOD
#define DIM 2
#else
#define DIM 3
#endif

#ifdef EWALD
#define I(a,b,c) [(((2*(a)+1)*(b))+(c))]
#define J(a,b,c) [(((a)*(b))+(c))]   
#endif

/* Skalarprodukt */
/* Vectors */
#define SPROD3D(a,b) (((a).x * (b).x) + ((a).y * (b).y) + ((a).z * (b).z))
#define SPROD2D(a,b) (((a).x * (b).x) + ((a).y * (b).y))
/* Arrays */
#define SPRODN3D(a,b) \
  (((a)[0] * (b)[0]) + ((a)[1] * (b)[1]) + ((a)[2] * (b)[2]))
#define SPRODN2D(a,b) (((a)[0] * (b)[0]) + ((a)[1] * (b)[1]))
/* Mixed Arrray, Vector */
#define SPRODX3D(a,b) \
   (((a)[0] * (b).x) + ((a)[1] * (b).y) + ((a)[2] * (b).z))
#define SPRODX2D(a,b) (((a)[0] * (b).x) + ((a)[1] * (b).y))
                           

#ifdef TWOD
#define SPROD(a,b)         SPROD2D(a,b)
#define SPRODN(a,b)        SPRODN2D(a,b)
#define SPRODX(a,b)        SPRODX2D(a,b)
#else
#define SPROD(a,b)         SPROD3D(a,b)
#define SPRODN(a,b)        SPRODN3D(a,b)
#define SPRODX(a,b)        SPRODX3D(a,b)
#endif

/* Dynamically allocated 3D arrray -- sort of */
#define PTR_3D(var,i,j,k,dim_i,dim_j,dim_k) \
  (((var) + ((i)*(dim_j)*(dim_k)) + ((j)*(dim_k)) + (k)))

/* Dynamically allocated 3D arrray -- half vector version */
#define PTR_3D_V(var,i,j,k,dim) \
  (((var) + ((i)*(dim.y)*(dim.z)) + ((j)*(dim.z)) + (k)))

/* Dynamically allocated 3D arrray -- full vector version */
#define PTR_3D_VV(var,coord,dim) \
  (((var) + ((coord.x)*(dim.y)*(dim.z)) + ((coord.y)*(dim.z)) + (coord.z)))

/* Dynamically allocated 2D arrray -- sort of */
#define PTR_2D(var,i,j,dim_i,dim_j) \
  (((var) + ((i)*(dim_j)) + (j)))

/* Dynamically allocated 2D arrray -- half vector version */
#define PTR_2D_V(var,i,j,dim) \
  (((var) + ((i)*(dim.y)) + (j)))

/* Dynamically allocated 2D arrray -- full vector version */
#define PTR_2D_VV(var,coord,dim) \
  (((var) + ((coord.x)*(dim.y)) + (coord.y)))

#ifdef TWOD
#define PTR     PTR_2D
#define PTR_V   PTR_2D_V
#define PTR_VV  PTR_2D_VV
#else
#define PTR     PTR_3D
#define PTR_V   PTR_3D_V
#define PTR_VV  PTR_3D_VV
#endif

/* simulation ensembles */
#define ENS_EMPTY     0
#define ENS_NVE       1
#define ENS_MIK       2
#define ENS_NVT       3
#define ENS_NPT_ISO   4
#define ENS_NPT_AXIAL 5
#define ENS_FRAC      8
#define ENS_SLLOD     9
#define ENS_NVX      11
#define ENS_STM      13
#define ENS_FTG      14
#define ENS_CG       15

/* output formats for distributions */
#define DIST_FORMAT_BINARY       1
#define DIST_FORMAT_ASCII_COORD  2
#define DIST_FORMAT_ASCII        3

/* All the logic in this program */
#define TRUE         1
#define FALSE        0

/* Some constants for Message passing, should all have unique values */
#define CELL_TAG   100
#define BUFFER_TAG 200
#define OUTBUF_TAG 300
#define INBUF_TAG  400
#define AT_BUF_TAG 500

/* some systems have different versions of trunc and floor float and double */
#if sgi || t3e
#ifdef DOUBLE
#define FLOOR floor
#else
#define FLOOR floorf
#endif
#else
#define FLOOR floor
#endif

#define TRUNC (int)

#ifndef M_PI
#define M_PI 4.0*atan(1.0)
#endif

/* Definition of the value that should be minimized */
#ifdef CGE   /* completely based on energy, no use of gradient information */
#define CGVAL tot_pot_energy
#endif
#ifdef CGEF   /* minimization of epot, but uses gradient information */
#define CGVAL tot_pot_energy
#endif
#ifdef CGF   /* minimization of total forces */
#define CGVAL fnorm
#endif
