
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
* imd_deform.c -- deform sample
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

#ifdef HOMDEF   /* homogeneous deformation with pbc */

/*****************************************************************************
*
* expand_sample()
*
*****************************************************************************/

void expand_sample(void)
{
  int k;
  
  /* Apply expansion */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (k=0; k<ncells; ++k) {
    int i;
    cell *p;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) {
      p->ort X(i) *= expansion.x;
      p->ort Y(i) *= expansion.y;
#ifndef TWOD
      p->ort Z(i) *= expansion.z;
#endif
    }
  }

  /* new box size */
#ifdef TWOD
  box_x.x *= expansion.x;  box_y.x *= expansion.x;
  box_x.y *= expansion.y;  box_y.y *= expansion.y;
#else
  box_x.x *= expansion.x;  box_x.y *= expansion.y;  box_x.z *= expansion.z;
  box_y.x *= expansion.x;  box_y.y *= expansion.y;  box_y.z *= expansion.z;
  box_z.x *= expansion.x;  box_z.y *= expansion.y;  box_z.z *= expansion.z;
#endif
  make_box();

} /* expand sample */


/*****************************************************************************
*
* shear_sample()
*
*****************************************************************************/

void shear_sample(void)
{
  int k;

  /* Apply shear */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (k=0; k<ncells; ++k) {
    int i;
    cell *p;
    real tmport[2];
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) {
      tmport[0]  = shear_factor.x * p->ort Y(i);
      tmport[1]  = shear_factor.y * p->ort X(i);
      p->ort X(i) += tmport[0];
      p->ort Y(i) += tmport[1];
    }
  }

  /* new box size */
  box_y.x += shear_factor.x * box_y.y;
  box_x.y += shear_factor.y * box_x.x;
  make_box();

} /* shear sample */


/*****************************************************************************
* 
* lin_deform()
*
*****************************************************************************/

void lin_deform(void)
{
   int k;
   real tmpbox[3];
   
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (k=0; k<ncells; ++k) {
    int i;
    cell *p;
    real tmport[3];
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) {
       /* transform atom positions */
      tmport[0] =   lindef_x.x * p->ort X(i) + lindef_x.y * p->ort Y(i)
#ifndef TWOD
                  + lindef_x.z * p->ort Z(i)
#endif
                                                 ;
      tmport[1] =   lindef_y.x * p->ort X(i) + lindef_y.y * p->ort Y(i)
#ifndef TWOD
                  + lindef_y.z * p->ort Z(i)
#endif
                                                 ;
#ifndef TWOD
      tmport[2] =   lindef_z.x * p->ort X(i) + lindef_z.y * p->ort Y(i)
                  + lindef_z.z * p->ort Z(i) ;
#endif

      p->ort X(i) += tmport[0];
      p->ort Y(i) += tmport[1];
#ifndef TWOD
      p->ort Z(i) += tmport[2];
#endif
    }
  }

  /* transform first box vector */
  tmpbox[0] =   lindef_x.x * box_x.x + lindef_x.y * box_x.y
#ifndef TWOD
              + lindef_x.z * box_x.z
#endif
                                             ;
  tmpbox[1] =   lindef_y.x * box_x.x + lindef_y.y * box_x.y
#ifndef TWOD
              + lindef_y.z * box_x.z
#endif
                                             ;
#ifndef TWOD
  tmpbox[2] =   lindef_z.x * box_x.x + lindef_z.y * box_x.y
              + lindef_z.z * box_x.z ;
#endif

  box_x.x += tmpbox[0];
  box_x.y += tmpbox[1];
#ifndef TWOD
  box_x.z += tmpbox[2];
#endif
  
  /* transform second box vector */
  tmpbox[0] =   lindef_x.x * box_y.x + lindef_x.y * box_y.y
#ifndef TWOD
              + lindef_x.z * box_y.z
#endif
                                             ;
  tmpbox[1] =   lindef_y.x * box_y.x + lindef_y.y * box_y.y
#ifndef TWOD
              + lindef_y.z * box_y.z
#endif
                                             ;
#ifndef TWOD
  tmpbox[2] =   lindef_z.x * box_y.x + lindef_z.y * box_y.y
              + lindef_z.z * box_y.z ;
#endif

  box_y.x += tmpbox[0];
  box_y.y += tmpbox[1];
#ifndef TWOD
  box_y.z += tmpbox[2];
#endif
  
  /* transform third box vector */
#ifndef TWOD
  tmpbox[0] =   lindef_x.x * box_z.x + lindef_x.y * box_z.y
              + lindef_x.z * box_z.z ;

  tmpbox[1] =   lindef_y.x * box_z.x + lindef_y.y * box_z.y
              + lindef_y.z * box_z.z ;

  tmpbox[2] =   lindef_z.x * box_z.x + lindef_z.y * box_z.y
              + lindef_z.z * box_z.z ;

  box_z.x += tmpbox[0];
  box_z.y += tmpbox[1];
  box_z.z += tmpbox[2];
#endif

  /* apply box changes */
  make_box();
 
} /* lin_deform */


#endif /* HOMDEF */


#ifdef DEFORM

/*****************************************************************************
*
* deform_sample()
*
*****************************************************************************/

void deform_sample(void) 
{
  int k;
  /* loop over all atoms */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (k=0; k<ncells; ++k) {
    int i;
    cell *p;
    int sort;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) {
      sort = p->sorte[i];
      /* move particles with virtual types  */
      p->ort X(i) += (deform_shift + sort)->x;
      p->ort Y(i) += (deform_shift + sort)->y;
#ifndef TWOD
      p->ort Z(i) += (deform_shift + sort)->z;
#endif
    }
  }
}

#endif /* DEFORM */

