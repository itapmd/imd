/******************************************************************************
*
* config.h -- Configuration rules and constants for the imd Package
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

/******************************************************************************
*
* Misc. Settings
*
******************************************************************************/

/* Double Precision is default */
#ifndef SINGLE
#define DOUBLE
#endif

/* define SAVEMEM for small cell structure            */
/* define AR for using actio-reactio principle        */
/* define USE_RUSAGE for using getrusage()-routines   */
/* define USE_CLOCK to use clock() instead of times() */

#ifdef MONOLJ
/* #define SAVEMEM
   #undef AR */
#endif

/******************************************************************************
*
* Statistical Ensembles
*
******************************************************************************/

/* Microcanonical */

#ifdef NVE

#endif

/* Constant Temperature & Pressure */

#ifdef NPT
#define DYNAMIC_CELLS
#endif

#ifdef NPT_axial
#define P_AXIAL
#endif

/* Canonical */

#ifdef NVT

#endif

/* Stadium damping for fracture */

#ifdef FRAC
#define NOPBC
#endif

/* Pulling on borders */

#ifdef PULL
#define NOPBC
#endif

/* Shear */

#ifdef MIKSHEAR
#define NOPBC
#endif

#ifdef HOMSHEAR
#define NOPBC
#endif

#ifdef MSQD
#define REFPOS
#endif

#ifdef CORRELATE
#define REFPOS
#endif

/******************************************************************************
*
* Architectures
*
******************************************************************************/

/* Single CPU, Workstation */

#ifdef RISC

#endif

/* Parallel Machine using Message Passing Interface */

#ifdef MPI

#endif

#ifdef CRAY 

#endif

#ifdef SX4
#undef USE_SOCKETS
#endif

/* Metacomputing */

#ifdef PACX
/* #define MPI_Ssend MPI_Send */
#endif

/******************************************************************************
*
* Constants
*
******************************************************************************/

/* Memory allocation increments for potential and cells */

#define PSTEP 50
#define CSTEP 10












