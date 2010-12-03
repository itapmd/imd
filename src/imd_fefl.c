
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
  real tmp_tot_ein_ene;
  real tmp;

  vektor d;
  real   dd,ee,ff;
  int    p_typ; 

  tot_harm_energy = 0.0;
  tmp_tot_ein_ene = 0.0;

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

      /* Addition of Einstein energy */
      ee  =   spring_rate[p_typ] * dd;
      tmp_tot_ein_ene += ee;

      /* Addition of pair forces and Einstein forces */      
      ff  =  -spring_rate[p_typ];
      KRAFT(p,i,X) = (1-lambda) * KRAFT(p,i,X) + lambda * ff * d.x; 
      KRAFT(p,i,Y) = (1-lambda) * KRAFT(p,i,Y) + lambda * ff * d.y; 
#ifndef TWOD
      KRAFT(p,i,Z) = (1-lambda) * KRAFT(p,i,Z) + lambda * ff * d.z; 
#endif      
    }
  }

  /* Add total Einstein energy */ 
  tot_harm_energy += tmp_tot_ein_ene;

#ifdef MPI
  MPI_Allreduce( &tot_harm_energy, &tmp, 1, REAL, MPI_SUM, cpugrid);
  tot_harm_energy = tmp;
#endif
}
