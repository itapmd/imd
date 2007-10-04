
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2007 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_cbe_calc_ppu.c -- calc_wp for the PPU
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"
#include "imd_cbe.h"

/* Macro for the calculation of the Lennard-Jones potential */
#define LJ(pot,grad,r2)  {                                               \
   int   const col4 = col*4;                                             \
   float const tmp2 = pt.lj_sig[col4] / r2;                              \
   float const tmp6 = tmp2 * tmp2 * tmp2;                                \
   pot  = pt.lj_eps[col4] * tmp6 * (tmp6 - 2.0) - pt.lj_shift[col4];     \
   grad = - 12.0 * pt.lj_eps[col4] * tmp6 * (tmp6 - 1.0) / r2;           \
}

#ifndef ON_PPU

/******************************************************************************
*
*  Lennard-Jones interactions on PPU - dummy routine, real one runs on SPU
*
******************************************************************************/

void calc_wp(wp_t * const wp) {}

/******************************************************************************
*
*  Neighbor tables on PPU - dummy routine, real one runs on SPU
*
******************************************************************************/

void calc_tb(wp_t * const wp) {}


#else  /* ON_PPU */

#ifdef CBE_DIRECT

/******************************************************************************
*
*  Lennard-Jones interactions on PPU - CBE_DIRECT version
*
******************************************************************************/

void calc_wp(wp_t *wp)
{
  cell_dta_t *p, *q;
  float d[4], f[4]; 
  float r2, *pi, *pj, *fi, *fj, pot, grad;
  int   i, c1, col, inc=pt.ntypes * pt.ntypes, n, j, m;

  wp->totpot = wp->virial = 0.0;
  p = wp->cell_dta;

  for (m=0; m<NNBCELL; m++) {

    q = wp->cell_dta + m;
    if (0==q->n) continue;

    for (i=0; i<p->n; i++) {

      short *ttb = q->tb + q->ti[2*i];

      c1 = pt.ntypes * p->typ[i];
      fi = p->force + 4*i;
      pi = p->pos   + 4*i;

      for (n=0; n<q->ti[2*i+1]; n++) {

        j    = ttb[n];
        pj   = q->pos + 4*j;
        d[0] = pj[0] - pi[0];
        d[1] = pj[1] - pi[1];
        d[2] = pj[2] - pi[2];
        r2   = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
        col  = c1 + q->typ[j];

        if (r2 <= pt.r2cut[4*col]) {

    	  /* Calculation of potential has been moved to a macro */
   	  LJ(pot,grad,r2)

#ifdef AR
          wp->totpot += pot;
          pot        *= 0.5;   /* avoid double counting */
          wp->virial -= r2  * grad;
#else
          pot        *= 0.5;   /* avoid double counting */
          wp->totpot += pot;
          wp->virial -= r2  * grad * 0.5;
#endif

          f[0] = d[0] * grad;
          f[1] = d[1] * grad;
          f[2] = d[2] * grad;

          fi[0] += f[0];
          fi[1] += f[1];
          fi[2] += f[2];
          fi[3] += pot;
#ifdef AR
          fj = q->force + 4*j;
          fj[0] -= f[0];
          fj[1] -= f[1];
          fj[2] -= f[2];
          fj[3] += pot;
#endif
	}
      }
    }
  }
}


/******************************************************************************
*
*  calc_tb: Neighbor tables on PPU -- CBE_DIRECT version
*
******************************************************************************/

static void calc_tb(wp_t * const wp) {

  cell_dta_t *p = wp->cell_dta;
  float cellsz = wp->totpot;
  int inc_int   = 128 / sizeof(int  ) - 1;
  int inc_short = 128 / sizeof(short) - 1;
  int at_inc = (2*p->n + inc_int) & (~inc_int);
  int m, nn=0;

  /* for each neighbor cell */
  for (m=0; m<NNBCELL; m++) {

    cell_dta_t *q = wp->cell_dta + m;

    int l=0, n=0, i;

    q->tb = p->tb + nn;

    if (0==q->n) continue;

    /* for each atom in cell */
    for (i=0; i<p->n; i++) {

      int    jstart, j, rr;
      vektor d1;

      d1.x = p->pos[4*i  ];
      d1.y = p->pos[4*i+1];
      d1.z = p->pos[4*i+2];

      /* for each neighboring atom */
      q->ti[l] = n;
#ifdef AR
      jstart = (m==0) ? i+1 : 0;
#else
      jstart = 0;
#endif
      for (j=jstart; j<q->n; j++) {
        vektor d;
        real   r2;
#ifndef AR
        if ((m==0) && (i==j)) continue;
#endif
        d.x = q->pos[4*j  ] - d1.x;
        d.y = q->pos[4*j+1] - d1.y;
        d.z = q->pos[4*j+2] - d1.z;
        r2  = SPROD(d,d);
        if (r2 < cellsz) q->tb[n++] = j;
      }
      q->ti[l+1] = n - q->ti[l];
      l += 2;

      /* if n is not divisible by 4, pad with copies of q->n */
      rr = n % 4;
      if (rr>0) for (j=rr; j<4; j++) q->tb[n++] = q->n;
    }

    /* enlarge n to next 128 byte boundary */
    n = (n + inc_short) & (~inc_short); 

    nn += n;
    if (nn > nb_max) {
      wp->flag = -1;
      return;
    }
  }
  wp->flag = nn;
}


#else  /* not CBE_DIRECT */

/******************************************************************************
*
*  Lennard-Jones interactions on PPU - indirect version
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

#endif  /* CBE_DIRECT */

#endif  /* ON_PPU */
