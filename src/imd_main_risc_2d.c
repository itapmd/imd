
/******************************************************************************
*
* imd_main_risc_2d.c -- main loop, risc specific part, two dimensions
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"


/*****************************************************************************
*
*  calc_forces()
*
*****************************************************************************/

void calc_forces(void)
{
  int n, k;

  /* clear global accumulation variables */
  tot_pot_energy = 0.0;
  virial = 0.0;
  vir_x  = 0.0;
  vir_y  = 0.0;

  /* clear per atom accumulation variables */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (k=0; k<ncells; ++k) {
    int  i;
    cell *p;
    p = cell_array + k;
    for (i=0; i<p->n; ++i) {
      p->kraft X(i) = 0.0;
      p->kraft Y(i) = 0.0;
      p->pot_eng[i] = 0.0;
#ifdef ORDPAR
      p->nbanz[i] = 0;
#endif
#ifdef NVX
      p->heatcond[i] = 0.0;
#endif     
#ifdef STRESS_TENS
      p->presstens X(i) = 0.0;
      p->presstens Y(i) = 0.0;
      p->presstens_offdia[i] = 0.0;
#endif      
    }
  }

  /* compute forces for all pairs of cells */
  for (n=0; n<4; ++n ) {
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tot_pot_energy,virial,vir_x,vir_y)
#endif
    for (k=0; k<npairs[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n]+k;
      pbc.x = P->ipbc[0] * box_x.x + P->ipbc[1] * box_y.x;
      pbc.y = P->ipbc[0] * box_x.y + P->ipbc[1] * box_y.y;
      do_forces(cell_array + P->np, cell_array + P->nq, pbc);
    }
  }
}



