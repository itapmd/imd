
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

/* we always need PAIR_POT, unless TERSOFF, KEATING, STIWEB, UNIAX, MONOLJ, */
/* SPRING, or PAIR_PRE are defined */
/* note that  PAIR_POT is the default if no interaction is specified     */
#if !(defined(TERSOFF) || defined(KEATING) || defined(STIWEB) || defined(UNIAX) || defined(MONOLJ) || defined(SPRING) || defined(PAIR_PRE))
#ifndef PAIR_POT
#define PAIR_POT
#endif
#endif

/* for EAM2, TTBP, and PAIR, we also need PAIR_POT */
#if (defined(EAM2) || defined(TTBP) || defined(PAIR))
#ifndef PAIR_POT
#define PAIR_POT 
#endif
#endif

/* PAIR is default */
#if !(defined(TERSOFF) || defined(KEATING) || defined(STIWEB) || defined(EWALD) || defined(SPRING) || defined(PAIR_PRE))
#ifndef PAIR
#define PAIR
#endif
#endif

/* If predefined potentials are used, no potential table is read in unless */
/* the option pair is used explicitly */
#if defined(PAIR_PRE) && !defined(PAIR)
#ifdef PAIR_POT
#undef PAIR_POT
#endif
#endif

/* shortcut for covalent interactions */
#if defined(KEATING)|| defined(TTBP) || defined(TERSOFF) || defined(STIWEB)
#define COVALENT
#endif

/* check whether we use AR across CPU borders */
#ifdef MPI

/* AR is the default. We could make the default machine dependent */
#define AR  

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
#ifdef CG
#define MAX_ATOM_SIZE 35
#else
#define MAX_ATOM_SIZE 25
#endif
#endif
#endif

/* number of data items by which an MPI buffer is increased */
#define BUFFER_SIZE_INC 1024 

/* buffer size for serial read_atoms */
#define INPUT_BUF_SIZE 16384

/* output buffer size (in bytes) */
#define OUTPUT_BUF_SIZE 131072

