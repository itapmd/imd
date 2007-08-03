#ifndef CALC_WP_C__
#define CALC_WP_C__

#include "imd_cbe.h"

/* This is basically the calc_wp rotuine cut & pasted from imd_forces_cbe */

/* Macro for the calculation of the Lennard-Jones potential */
#define LJ(pot,grad,r2)  {                                               \
   int const col4 = col*4;                                               \
   flt const tmp2 = pt.lj_sig[col4] / r2;                                \
   flt const tmp6 = tmp2 * tmp2 * tmp2;                                  \
   pot  = pt.lj_eps[col4] * tmp6 * (tmp6 - 2.0) - pt.lj_shift[col4];     \
   grad = - 12.0 * pt.lj_eps[col4] * tmp6 * (tmp6 - 1.0) / r2;           \
}


/******************************************************************************
*
*  calc_wp
*
******************************************************************************/

void calc_wp(wp_t *wp)
{
  typedef flt fltvec[4];
  fltvec d, f; 
  flt r2, *pi, *pj, *fi, *fj, pot, grad;
  int i, c1, col, inc=pt.ntypes * pt.ntypes, n, j;

  wp->totpot = wp->virial = 0.0;
  for (i=0; i<4*wp->n2; i++) { 
      wp->force[i] = 0.0;
  }

  /* fetch sizes and pointers */
  /* fetch data at wp->pos, wp->typ, wp->ti */

  for (i=0; i<wp->n1; i++) {

    short *ttb;

    /* fetch data at wp->tb + wp->ti[2*i] */
    ttb = wp->tb + wp->ti[2*i];

    c1 = pt.ntypes * wp->typ[i];
    fi = wp->force + 4*i;
    pi = wp->pos   + 4*i;

    for (n=0; n<wp->ti[2*i+1]; n++) {

      j    = ttb[n];
      pj   = wp->pos + 4*j;
      d[0] = pj[0] - pi[0];
      d[1] = pj[1] - pi[1];
      d[2] = pj[2] - pi[2];
      r2   = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
      col  = c1 + wp->typ[j];

      if (r2 <= pt.r2cut[4*col]) {

  	  /* Calculation of potential has been moved to a macro */
   	  LJ(pot,grad,r2)

          wp->totpot += pot;
          pot        *= 0.5;   /* avoid double counting */
          wp->virial -= r2  * grad;

          f[0] = d[0] * grad;
          f[1] = d[1] * grad;
          f[2] = d[2] * grad;

          fi[0] += f[0];
          fi[1] += f[1];
          fi[2] += f[2];
          fi[3] += pot;

          fj = wp->force + 4*j;
          fj[0] -= f[0];
          fj[1] -= f[1];
          fj[2] -= f[2];
          fj[3] += pot;

      }
    }
  }
  /* return wp->force, wp->totpot, wp->virial */
}



#endif
