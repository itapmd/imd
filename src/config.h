
/******************************************************************************
*
* config.h -- Configuration rules and constants for the imd Package
*
******************************************************************************/

/******************************************************************************
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
#define SAVEMEM
#undef AR
#define NODBG_DIST
#endif

/******************************************************************************
*
* Statistical Ensembles
*
******************************************************************************/

#ifdef NPT_axial
#define P_AXIAL
#endif

/* heat transport */
#ifdef NVX
#define RNEMD
#endif
#if defined(NVX) || defined(RNEMD)
#define TRANSPORT
#endif

#ifdef MSQD
#define REFPOS
#endif

#ifdef CORRELATE
#define REFPOS
#endif

#ifdef SHOCK
#define STRESS_TENS
#endif

#ifdef EAM
#undef AR
#endif

#ifdef EAM2
#undef AR
#endif

#if defined(TTBP) || defined(TERSOFF)
#define COVALENT
#endif

#ifdef COVALENT
#ifndef AR
#define AR
#endif
#endif

#ifdef UNIAX
#undef AR
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


/******************************************************************************
*
* Constants
*
******************************************************************************/

/* memory allocation increment for potential */
#define PSTEP 50

/* security margin for buffer sizes */
#define CSTEP 10

/* interval for checking MPI buffers */
#define BUFSTEP 100

/* maximum number of data items an atom can have */
#ifdef TWOD
#define MAX_ATOM_SIZE 20
#else
#ifdef UNIAX
#define MAX_ATOM_SIZE 40
#else
#define MAX_ATOM_SIZE 25
#endif
#endif

/* number of data items by which an MPI buffer is increased */
#define BUFFER_SIZE_INC 1024 

/* buffer size for serial read_atoms */
#define INPUT_BUF_SIZE 16384
