
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2010 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_fefl.c -- compute free energy using the Frenkel-Ladd method
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/******************************************************************************
*
* compute forces due to harmonic force wrt to reference knofiguration
*
******************************************************************************/

/* tot_pot_energy, KRAFT: global */

void calc_fefl(void)
{
  int k, i, n;
  real tmpvec1[8], tmpvec2[8];
  real tmp_virial, tmp_tot_pot_ene;

#ifdef P_AXIAL
  vektor tmp_vir_vect;
#endif

  vektor d;
  real   dd,ee,ff;
  int    p_typ; 

  tmp_tot_pot_ene = 0.0;
  tmp_virial     = 0.0;
#ifdef P_AXIAL
  tmp_vir_vect.x = 0.0;
  tmp_vir_vect.y = 0.0;
#ifndef TWOD
  tmp_vir_vect.z = 0.0;
#endif
#endif

  for (k=0; k<NCELLS; ++k) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      
      /* Elongation from reference position */
      d.x =   ORT(p,i,X)-REF_POS(p,i,X); 
      d.y =   ORT(p,i,Y)-REF_POS(p,i,Y);
#ifndef TWOD
      d.z =   ORT(p,i,Z)-REF_POS(p,i,Z);
      dd  =   SPROD(d,d);
#else
      dd  =   SPROD2D(d,d);
#endif
      
      p_typ = SORTE(p,i);

      /* Addition pair energies and Einstein energy */
      ee  =   0.5 * spring_rate[p_typ] * dd;
      POTENG(p,i) = (1-lambda) * POTENG(p,i) + lambda * ee;
      tmp_tot_pot_ene += ee;

      /* Addition pair forces and Einstein forces */      
      ff  =   spring_rate[p_typ];
      KRAFT(p,i,X) = (1-lambda) * KRAFT(p,i,X) + lambda * ff * d.x; 
      KRAFT(p,i,Y) = (1-lambda) * KRAFT(p,i,Y) + lambda * ff * d.y; 
#ifndef TWOD
      KRAFT(p,i,Z) = (1-lambda) * KRAFT(p,i,Z) + lambda * ff * d.z; 
#endif      
      /* Sign correct ? */
#ifdef P_AXIAL
      tmp_vir_vect.x += d.x * d.x * ff;
      tmp_vir_vect.y += d.y * d.y * ff;
#ifndef TWOD
      tmp_vir_vect.z += d.z * d.z * ff;
#endif
#else
      tmp_virial     += dd * ff;
#endif

    }
  }

  /* Add total pair energy and total pair virial to Einstein energy and virial */ 
  tot_pot_energy = (1-lambda) * tot_pot_energy + lambda * tmp_tot_pot_ene;
#ifdef P_AXIAL
  vir_xx = (1-lambda) * vir_xx + lambda * tmp_vir_vect.x;
  vir_yy = (1-lambda) * vir_yy + lambda * tmp_vir_vect.y;
#ifndef TWOD
  vir_zz = (1-lambda) * vir_zz + lambda * tmp_vir_vect.z;
#endif
#else
  virial = (1-lambda) * virial + lambda * tmp_virial;
#endif

  /* tot_pot_energy muss noch korrigiert werden */

#ifdef MPI
#ifdef TWOD
  /* sum up results of different CPUs */
  tmpvec1[0] = tot_pot_energy;
  tmpvec1[1] = virial;
  tmpvec1[2] = vir_xx;
  tmpvec1[3] = vir_yy;
  /* Die Nichtdiagonalkomponenten werden tatsächlich nirgendwo in IMD berechnet! 
     Sie sollten sich auch bei allen Zentralkräften effektiv wegheben */
  /* tmpvec1[4] = vir_xy; */

  MPI_Allreduce( tmpvec1, tmpvec2, 5, REAL, MPI_SUM, cpugrid); 
#else
  /* tmpvec1[5] = vir_yz;
     tmpvec1[6] = vir_zx; */
  tmpvec1[7] = vir_zz;

  MPI_Allreduce( tmpvec1, tmpvec2, 8, REAL, MPI_SUM, cpugrid); 
#endif
#ifdef TWOD
  tot_pot_energy = tmpvec2[0];
  virial         = tmpvec2[1];
  vir_xx         = tmpvec2[2];
  vir_yy         = tmpvec2[3];
  /*  vir_xy         = tmpvec2[4]; */
#else
  /*  vir_yz         = tmpvec2[5];
      vir_zx         = tmpvec2[6]; */
  vir_zz         = tmpvec2[7];
#endif
#endif
}
