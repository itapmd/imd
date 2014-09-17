
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2012 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* constants.h 
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

/* How many dimension are there? */
#ifdef TWOD
#define DIM 2
#else
#define DIM 3
#endif

/* simulation ensembles */
#define ENS_EMPTY     0
#define ENS_NVE       1
#define ENS_MIK       2
#define ENS_NVT       3
#define ENS_NPT_ISO   4
#define ENS_NPT_AXIAL 5
#define ENS_GLOK      6
#define ENS_FRAC      8
#define ENS_SLLOD     9
#define ENS_NVX      11
#define ENS_STM      13
#define ENS_FTG      14
#define ENS_CG       15
#define ENS_FINNIS   16
#define ENS_TTM      17

/* FCS methods */
#define FCS_METH_EMPTY  0
#define FCS_METH_DIRECT 1
#define FCS_METH_PEPC   2
#define FCS_METH_FMM    3
#define FCS_METH_P3M    4
#define FCS_METH_P2NFFT 5
#define FCS_METH_VMG    6
#define FCS_METH_PP3MG  7

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
#define ANNOUNCE_TAG 600

/* Definition of the value that should be minimized */
#define CGE  0 /* completely based on energy, no use of gradient information */
#define CGEF 1 /* minimization of epot, but uses gradient information */

#ifdef LOADBALANCE
#define LB_REAL_CELL 1
#define LB_BUFFER_CELL 2
#define LB_EMPTY_CELL 3
#define LB_NON_PBC_BUFFER_CELL 4
#define LB_CELL_SUBLEVELS 11 /*each cell is divided into sublevels for finer approximations of cpu domains*/

#define LB_SEND_FORCE   1
#define LB_SEND_CELL    2
#define LB_DONT_SEND    0

#endif
