
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

#ifdef DAMP /* deform the stadium correspondingly */
#ifdef TWOD
  center.x   *= expansion.x;  center.y   *= expansion.y; 
  stadium.x  *= expansion.x;  stadium.y  *= expansion.y;
  stadium2.x *= expansion.x;  stadium2.y *= expansion.y;
#else
  center.x   *= expansion.x;  center.y   *= expansion.y;  center.z   *= expansion.z;
  stadium.x  *= expansion.x;  stadium.y  *= expansion.y;  stadium.z  *= expansion.z;
  stadium2.x *= expansion.x;  stadium2.y *= expansion.y;  stadium2.z *= expansion.z;
#endif
#endif
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

#ifdef DAMP /* deform the stadium correspondingly */
  tmpbox = center.x;
  center.x += shear_factor.x * center.y;
  center.y += shear_factor.y * tmpbox;

  tmpbox = stadium.x;
  stadium.x += shear_factor.x * stadium.y;
  stadium.y += shear_factor.y * tmpbox;

  tmpbox = stadium2.x;
  stadium2.x += shear_factor.x * stadium2.y;
  stadium2.y += shear_factor.y * tmpbox;

#endif
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
 
#ifdef DAMP /* deform the stadium correspondingly */

  tmpbox[0] = scale * SPROD(dx,center);
  tmpbox[1] = scale * SPROD(dy,center);
#ifndef TWOD
  tmpbox[2] = scale * SPROD(dz,center);
#endif

  center.x += tmpbox[0];
  center.y += tmpbox[1];
#ifndef TWOD
  center.z += tmpbox[2];
#endif

  tmpbox[0] = scale * SPROD(dx,stadium);
  tmpbox[1] = scale * SPROD(dy,stadium);
#ifndef TWOD
  tmpbox[2] = scale * SPROD(dz,stadium);
#endif

  stadium.x += tmpbox[0];
  stadium.y += tmpbox[1];
#ifndef TWOD
  stadium.z += tmpbox[2];
#endif


  tmpbox[0] = scale * SPROD(dx,stadium2);
  tmpbox[1] = scale * SPROD(dy,stadium2);
#ifndef TWOD
  tmpbox[2] = scale * SPROD(dz,stadium2);
#endif

  stadium2.x += tmpbox[0];
  stadium2.y += tmpbox[1];
#ifndef TWOD
  stadium2.z += tmpbox[2];
#endif

#endif


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
  sym_tensor pt;
  real pp;

#ifdef STRESS_TENS

  /* here, we support relaxation to arbitrary external pressure tensor */
  calc_tot_presstens();
  pt.xx = tot_presstens.xx / volume - presstens_ext.xx;
  pt.yy = tot_presstens.yy / volume - presstens_ext.yy;
  pt.xy = tot_presstens.xy / volume - presstens_ext.xy;
#ifndef TWOD
  pt.zz = tot_presstens.zz / volume - presstens_ext.zz;
  pt.yz = tot_presstens.yz / volume - presstens_ext.yz;
  pt.zx = tot_presstens.zx / volume - presstens_ext.zx;
  pp = (pt.xx + pt.yy + pt.zz) / 3.0;
#else
  pp = (pt.xx + pt.yy) / 2.0;
#endif
  if ((relax_mode == RELAX_FULL) || (relax_mode == RELAX_AXIAL)) {
    dx.x = pp / bulk_module + (pt.xx - pp) / shear_module;
    dy.y = pp / bulk_module + (pt.yy - pp) / shear_module;
#ifndef TWOD
    dz.z = pp / bulk_module + (pt.zz - pp) / shear_module;
#endif
  } else {
    dx.x = pp / bulk_module;
    dy.y = pp / bulk_module;
#ifndef TWOD
    dz.z = pp / bulk_module;
#endif
  }
  if (relax_mode == RELAX_FULL) {
    dx.y  = dy.x = pt.xy / shear_module;
#ifndef TWOD
    dy.z  = dz.y = pt.yz / shear_module;
    dz.x  = dx.z = pt.zx / shear_module;
#endif
  }

#else  /* not STRESS_TENS */

  /* here, we support only relaxation to scalar pressure zero */
  dx.x = pressure / bulk_module;
  dy.y = pressure / bulk_module;
#ifndef TWOD
  dz.z = pressure / bulk_module;
#endif

#endif /* not STRESS_TENS */

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

