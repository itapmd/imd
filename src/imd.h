/******************************************************************************
*
* imd.h -- Header file for all Modules of the imd Package
*
*
******************************************************************************/


/* C stuff */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#ifndef USE_CLOCK
#include <sys/times.h>
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
