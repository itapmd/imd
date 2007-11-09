/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2007 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

/******************************************************************************
*
*  Main file for PyIMD, the pythonized version of IMD
*
*  As PyIMD requires no main(), this file contains only some 
*  interfacing routines, and provides the global variables
* 
******************************************************************************/

#define MAIN
#include "imd.h"

/* get an element of a real array */
real RealGetElm( real *p, int i) {
  return p[i];
}

/* set an element of a real array */
void RealSetElm( real *p, int i, real val) {
  p[i] = val;
}

/* get an element of an int array */
int IntGetElm( int *p, int i) {
  return p[i];
}

/* set an element of an int array */
void IntSetElm( int *p, int i, int val) {
  p[i] = val;
}
