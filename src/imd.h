
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
* imd.h -- Header file for all modules of IMD
*
******************************************************************************/

/* C stuff */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* support for timers */
#ifndef MPI
#include <unistd.h>
#ifdef USE_RUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#else
#include <sys/times.h>
#include <sys/types.h>
#endif
#endif

/* Machine specific headers */
#ifdef MPI
#include <mpi.h>
#endif
#ifdef OMP
#include <omp.h>
#endif

/* Konfiguration */
#include "config.h"

/* Data types */
#include "types.h"

/* Some makros */
#include "makros.h"

/* Socket headers */
#ifdef USE_SOCKETS
#include "sockets.h"
#endif

/* Function Prototypes */
#include "prototypes.h"

/* Global Variables */
#include "globals.h"
