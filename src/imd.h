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
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#ifndef USE_CLOCK
#include <sys/times.h>
#endif
#ifndef HAVE_NO_STDLIB
#include <stdlib.h>
#else
#include <malloc.h>
#endif
#ifdef sgi
#include <bstring.h>
#endif
#ifdef OMP
#include <omp.h>
#endif
/* #include <sys/stat.h> */
/* #include <sys/types.h> */


/* Machine specific headers */
#ifdef MPI
#include <mpi.h>
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
