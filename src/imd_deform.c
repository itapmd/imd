
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2004 Institute for Theoretical and Applied Physics,
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
  for (k=0; k<NCELLS; ++k) {
    int i;
    cell *p;
    p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      ORT(p,i,X) *= expansion.x;
      ORT(p,i,Y) *= expansion.y;
#ifndef TWOD
      ORT(p,i,Z) *= expansion.z;
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
  real tmpbox;

  /* Apply shear */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (k=0; k<NCELLS; ++k) {
    int i;
    cell *p;
    real tmport[2];
    p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      tmport[0]   = shear_factor.x * ORT(p,i,Y);
      tmport[1]   = shear_factor.y * ORT(p,i,X);
      ORT(p,i,X) += tmport[0];
      ORT(p,i,Y) += tmport[1];
    }
  }

  /* new box size */
  tmpbox = box_x.x;
  box_x.x += shear_factor.x * box_x.y;
  box_x.y += shear_factor.y * tmpbox;

  tmpbox = box_y.x;
  box_y.x += shear_factor.x * box_y.y;
  box_y.y += shear_factor.y * tmpbox;

#ifndef TWOD
  tmpbox = box_z.x;
  box_z.x += shear_factor.x * box_z.y;
  box_z.y += shear_factor.y * tmpbox;
#endif

  make_box();

} /* shear sample */


/*****************************************************************************
* 
* lin_deform()
*
*****************************************************************************/

#ifdef TWOD
void lin_deform(vektor dx, vektor dy,            real scale)
#else
void lin_deform(vektor dx, vektor dy, vektor dz, real scale)
#endif
{
   int k;
   real tmpbox[3];
   
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (k=0; k<NCELLS; ++k) {
    int i;
    cell *p;
    real tmport[3];
    p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      /* transform atom positions */
#ifdef TWOD
      tmport[0] = dx.x * ORT(p,i,X) + dx.y * ORT(p,i,Y);
      tmport[1] = dy.x * ORT(p,i,X) + dy.y * ORT(p,i,Y);
#else
      tmport[0] = dx.x * ORT(p,i,X) + dx.y * ORT(p,i,Y) + dx.z * ORT(p,i,Z);
      tmport[1] = dy.x * ORT(p,i,X) + dy.y * ORT(p,i,Y) + dy.z * ORT(p,i,Z);
      tmport[2] = dz.x * ORT(p,i,X) + dz.y * ORT(p,i,Y) + dz.z * ORT(p,i,Z);
#endif
      ORT(p,i,X) += scale * tmport[0];
      ORT(p,i,Y) += scale * tmport[1];
#ifndef TWOD
      ORT(p,i,Z) += scale * tmport[2];
#endif
    }
  }

  /* transform first box vector */
  tmpbox[0] = scale * SPROD(dx,box_x);
  tmpbox[1] = scale * SPROD(dy,box_x);
#ifndef TWOD
  tmpbox[2] = scale * SPROD(dz,box_x);
#endif

  box_x.x += tmpbox[0];
  box_x.y += tmpbox[1];
#ifndef TWOD
  box_x.z += tmpbox[2];
#endif
  
  /* transform second box vector */
  tmpbox[0] = scale * SPROD(dx,box_y);
  tmpbox[1] = scale * SPROD(dy,box_y);
#ifndef TWOD
  tmpbox[2] = scale * SPROD(dz,box_y);
#endif

  box_y.x += tmpbox[0];
  box_y.y += tmpbox[1];
#ifndef TWOD
  box_y.z += tmpbox[2];
#endif

  /* transform third box vector */
#ifndef TWOD
  tmpbox[0] = scale * SPROD(dx,box_z);
  tmpbox[1] = scale * SPROD(dy,box_z);
  tmpbox[2] = scale * SPROD(dz,box_z);

  box_z.x += tmpbox[0];
  box_z.y += tmpbox[1];
  box_z.z += tmpbox[2];
#endif
 
  /* apply box changes */
  make_box();
 
} /* lin_deform */

/*****************************************************************************
*
* relax_pressure()
*
*****************************************************************************/

void relax_pressure()
{
#ifdef TWOD
  vektor dx = {0.0, 0.0}, dy = {0.0, 0.0};
#else
  vektor dx = {0.0, 0.0, 0.0}, dy = {0.0, 0.0, 0.0}, dz = {0.0, 0.0, 0.0};
#endif

  dx.x = pressure / bulk_module;
  dy.y = pressure / bulk_module;
#ifndef TWOD
  dz.z = pressure / bulk_module;
#endif

#ifdef STRESS_TENS
  if ((relax_mode == RELAX_FULL) || (relax_mode == RELAX_AXIAL)) {
    calc_tot_presstens();
    dx.x += (tot_presstens.xx - pressure) / shear_module;
    dy.y += (tot_presstens.yy - pressure) / shear_module;
#ifndef TWOD
    dz.z += (tot_presstens.zz - pressure) / shear_module;
#endif
  }
  if (relax_mode == RELAX_FULL) {
    dx.y  = dy.x = tot_presstens.xy / shear_module;
#ifndef TWOD
    dy.z  = dz.y = tot_presstens.yz / shear_module;
    dz.x  = dx.z = tot_presstens.zx / shear_module;
#endif
  }
#endif /* STRESS_TENS */

#ifdef TWOD
  lin_deform(dx, dy,     relax_rate);
#else
  lin_deform(dx, dy, dz, relax_rate);
#endif
}

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
  for (k=0; k<NCELLS; ++k) {
    int i;
    cell *p;
    int sort;
    vektor ort;
    real shear;

    p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      sort = VSORTE(p,i);

      if ( *(shear_def + sort) == 1 ) {
	  ort.x = ORT(p,i,X) - (deform_base + sort)->x;
	  ort.y = ORT(p,i,Y) - (deform_base + sort)->y;
#ifndef TWOD
	  ort.z = ORT(p,i,Z) - (deform_base + sort)->z;
#endif
	  shear = SPROD( *(deform_shear + sort), ort);
	}
      else
	shear = 1.0;

      /* move particles with virtual types  */
      ORT(p,i,X) += shear * deform_size * (deform_shift + sort)->x;
      ORT(p,i,Y) += shear * deform_size * (deform_shift + sort)->y;
#ifndef TWOD
      ORT(p,i,Z) += shear * deform_size * (deform_shift + sort)->z;
#endif
    }
  }
}

#endif /* DEFORM */

