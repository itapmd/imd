
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

/* double precision is the default */
#ifndef SINGLE
#define DOUBLE
#endif
/* the world record switch */
#ifdef MONOLJ
#define NODBG_DIST
#endif

/* define USE_RUSAGE for using getrusage()-routines   */
/* define USE_CLOCK to use clock() instead of times() */

/******************************************************************************
*
* Interactions
*
******************************************************************************/

/* backwards compatibility for potential files without headers */
#ifdef EAM2
#define DEFAULT_POTFILE_TYPE 2
#else
#define DEFAULT_POTFILE_TYPE 1
#endif

/* we always need PAIR_POT, unless TERSOFF, UNIAX, or MONOLJ are defined */
/* note that  PAIR_POT is the default if no interaction is specified     */
#if !(defined(TERSOFF) || defined(UNIAX) || defined(MONOLJ))
#ifndef PAIR_POT
#define PAIR_POT
#endif
#endif

/* for EAM2, TTBP, and EWALD, we also need PAIR_POT */
#if (defined(EAM2) || defined(TTBP) || defined(EWALD))
#ifndef PAIR_POT
#define PAIR_POT 
#endif
#endif

/* shortcut for covalent interactions */
#if defined(TTBP) || defined(TERSOFF)
#define COVALENT
#endif

/* check whether we use AR across CPU borders */
#ifdef MPI

/* AR is the default. We could make the default machine dependent */
#define AR  

/* UNIAX must not use AR - not implemented */
#ifdef UNIAX
#undef AR
#endif

/* for COVALENT, AR *must* be set */
#ifdef COVALENT
#ifndef AR
#define AR
#endif
#endif

#endif /* MPI */

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

/******************************************************************************
*
* Architectures
*
******************************************************************************/

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

/* output buffer size (in bytes) */
#define OUTPUT_BUF_SIZE 131072
