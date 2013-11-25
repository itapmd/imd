
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

#ifdef CBE
#define NBL
#define SINGLE
#define MEMALIGN
#endif

/* double precision is the default */
#ifndef SINGLE
#define DOUBLE
#endif

/* the world record switch */
#ifdef MONOLJ
#define NODBG_DIST
#endif

/* KIM only works with neighbor lists enabled */
#ifdef KIM
#define NBL
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

/* we always need PAIR, unless MEAM, KEATING, UNIAX, FCS, or EWALD */
/* Note that PAIR is the default, if no interaction is specified */
#if !(defined(MEAM) || defined(KEATING) || defined(UNIAX) || defined(USEFCS) || defined(EWALD) || defined(TERSOFFMOD))
#ifndef PAIR
#define PAIR
#endif
#endif

/* ADP also implies EAM2 */
#if defined(ADP) && !defined(EAM2)
#define EAM2
#endif

/* for EAM2 and TTBP, we also need PAIR */
#if (defined(EAM2) || defined(TTBP))
#ifndef PAIR
#define PAIR
#endif
#endif

/* VARCHG and DIPOLE require COULOMB */
#if (defined(VARCHG) || defined(DIPOLE))
#ifndef COULOMB
#define COULOMB
#endif
#endif

/* for SM we also need EWALD (unless we have NBL) and VARCHG */
#if defined SM
#if !defined(EWALD) && !defined(NBL)
#define EWALD
#endif
#ifndef VARCHG
#define VARCHG
#endif
#endif

/* default short-range potential for DIPOLE is MORSE */
#if (defined(DIPOLE) && !defined(BUCK))
#define MORSE
#endif

/* shortcut for covalent interactions */
#if defined(MEAM) || defined(KEATING) || defined(TTBP) || defined(TERSOFF) || defined(TERSOFFMOD) || defined(STIWEB) || defined(TERNBCC) || defined(XT) || defined(CNA)
#define COVALENT
#endif

/* Enabling creation of nearest neighbor tables, using the COVALENT tables */
/* Combining covalent interaction with ADA or NYETENSOR may by risky*/
#ifdef ADA
#define NNBR_TABLE
#define ADA_ACKLAND_CONFIG   0
#define ADA_FCC_CONFIG 1
#define ADA_BCC_CONFIG 2
#endif


/* percentage of r2_cut, in which generated potentials are corrected
   so that the forces go continuously to zero */
#define POT_TAIL 0.05

/* with VEC or MPI or NBL we always use buffer cells */
#if defined(MPI) || defined(VEC) || defined(NBL)
#define BUFCELLS
#endif

#ifdef VEC
/* #define MULTIPOT */
#define N_POT_TAB 128
#define NBLIST
#endif

#ifdef NBL
#define NBLIST
#endif

#ifdef BUFCELLS

/* AR is the default. We could make the default machine dependent */
#if !defined(CBE) || defined(CBE2)
#define AR
#endif

#ifdef AR
#  define NNBCELL 14
#else
#  define NNBCELL 27
#endif

/* for COVALENT, AR *must* be set */
#ifdef COVALENT
#ifndef AR
#define AR
#endif
#endif

#endif /* BUFCELLS */

#ifdef KIM
#undef PAIR
#undef AR
#undef NNBCELL
#define NNBCELL 27
#endif

/******************************************************************************
*
* Statistical Ensembles and Integrators
*
******************************************************************************/

#if defined(CG) || defined(MIK) || defined(GLOK) || defined(DEFORM)
#ifndef FNORM
#define FNORM
#endif
#endif

/* relaxation integrators */
#if defined(MIK) || defined(GLOK) || defined(CG)
#define RELAX
#endif

#if defined(AND) || defined(NVT) || defined(NPT) || defined(FRAC) || defined(FINNIS) || defined(STM) ||defined(BER)

#define TEMPCONTROL
#endif

/* GLOK is NVE with additional features */
#ifdef GLOK
#ifndef NVE
#define NVE
#endif
#ifndef MIX
#define MIX
#endif
#endif

#ifdef NPT_axial
#define P_AXIAL
#endif

/* heat transport */
#ifdef HC
#ifndef STRESS_TENS
#define STRESS_TENS
#endif
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

#ifdef NMOLDYN
#define REFPOS
#endif

#ifdef FEFL
#define REFPOS
#endif

/******************************************************************************
*
* Architectures
*
******************************************************************************/

#if defined(MPI) && (MPI_VERSION==2)
#define MPI2
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

#if defined(COVALENT) || defined(NNBR_TABLE)
#ifdef MEAM
#define NEIGH_LEN_INIT 12
#else
#define NEIGH_LEN_INIT 4
#endif
#define NEIGH_LEN_INC  2
#endif

#ifdef CNA
#define MAX_NEIGH 12
#define MAX_BONDS 24
#define MAX_TYPES 25
#endif

/* pressure relaxation */
#define RELAX_FULL     0
#define RELAX_AXIAL    1
#define RELAX_ISO      2

/* total buffer size in MB for serial read_atoms */
#define INPUT_BUF_SIZE 64

/* output buffer size in MB */
#define OUTPUT_BUF_SIZE 32

/* maximum number of items per atom in config file (plus number and type) */
#define MAX_ITEMS_CONFIG 16
