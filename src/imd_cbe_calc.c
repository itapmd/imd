
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
* imd_cbe_calc.c
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include <stdio.h>
#include "imd_cbe.h"

#ifdef __SPU__

#include <spu_intrinsics.h>

/******************************************************************************
*
*  Lennard-Jones interaction work package for the SPU
*
******************************************************************************/

void calc_wp(wp_t *wp)
{
  vector float f00  = spu_splats( (float)   0.0  );
  vector float f001 = spu_splats( (float)   0.001);
  vector float f05n = spu_splats( (float)  -0.5  );
  vector float f05l = {0.0, 0.0, 0.0, 0.5};
  vector float f10  = spu_splats( (float)   1.0  );
  vector float f20  = spu_splats( (float)   2.0  );
  vector float f12  = spu_splats( (float) -12.0  );
  vector float vir  = spu_splats( (float)   0.0  );
  vector float r2cut   = pt.r2cut   [0];
  vector float ljsig   = pt.lj_sig  [0];
  vector float ljeps   = pt.lj_eps  [0];
  vector float ljshift = pt.lj_shift[0];
  vector signed   int  i00 = spu_splats( (int) 0 );
  vector unsigned char s0  = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19};
  vector unsigned char s1  = {0,1,2,3,4,5,6,7,8,9,10,11,20,21,22,23};
  vector unsigned char s2  = {0,1,2,3,4,5,6,7,8,9,10,11,24,25,26,27};
  vector unsigned char s3  = {0,1,2,3,4,5,6,7,8,9,10,11,28,29,30,31};
  vector unsigned char ss0 = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
  vector unsigned char ss1 = {4,5,6,7,4,5,6,7,4,5,6,7,4,5,6,7};
  vector unsigned char ss2 = {8,9,10,11,8,9,10,11,8,9,10,11,8,9,10,11};
  vector unsigned char ss3 = {12,13,14,15,12,13,14,15,12,13,14,15,12,13,14,15};
  vector float d0, d1, d2, d3, d20, d21, d22, d23, r2, r2i;
  vector float tmp2, tmp3, tmp4, tmp6, tmp7, pot;
  vector float grad, pots, ff0, ff1, ff2, ff3;
  vector float d2a, d2b, ffa, ffb, dummy=f00, *pi, *fi, ff0s, ff1s, ff2s, ff3s;
  vector unsigned int ms, ms1, ms2, ms3;
  vector signed int tj;
  int    i, c1, n, j0, j1, j2, j3;

  /* clear accumulation variables */
  wp->totpot = 0.0;
  wp->virial = 0.0;
  for (i=0; i<wp->n2; i++) wp->force[i] = f00;

  /* fetch sizes and pointers
     fetch data at wp->pos, wp->typ, wp->ti */
  for (i=0; i<wp->n1; i++) {

    short *ttb;

    /* fetch data at wp->tb + wp->ti[2*i] */
    ttb = wp->tb + wp->ti[2*i];

#ifndef MONO
    c1 = 2 * wp->typ[i];
#endif

    fi = wp->force + i;
    pi = wp->pos   + i;

    /* clear accumulation variables */
    pots = f00;  ff0s = f00;  ff1s = f00;  ff2s = f00;  ff3s = f00;

    /* we treat four neighbors at a time, so that we can vectorize */
    /* ttb is padded with copies of i, which have to be masked     */
    for (n = 0; n < wp->ti[2*i+1]; n += 4) {

      /* indices of neighbors */
      j0  = ttb[n  ];
      j1  = ttb[n+1];
      j2  = ttb[n+2];
      j3  = ttb[n+3];

#ifndef MONO
      /* if not MONO, we assume up to two atom types */
      /* mask for type dependent selections */
      tj = spu_promote( wp->typ[j0],     0 );
      tj = spu_insert ( wp->typ[j1], tj, 1 );
      tj = spu_insert ( wp->typ[j2], tj, 2 );
      tj = spu_insert ( wp->typ[j3], tj, 3 );
      ms = spu_cmpeq( tj, i00 );
      r2cut   = spu_sel( pt.r2cut   [c1+1], pt.r2cut   [c1], ms );
      ljsig   = spu_sel( pt.lj_sig  [c1+1], pt.lj_sig  [c1], ms );
      ljeps   = spu_sel( pt.lj_eps  [c1+1], pt.lj_eps  [c1], ms );
      ljshift = spu_sel( pt.lj_shift[c1+1], pt.lj_shift[c1], ms );
#endif

      /* distance vectors */
      d0  = spu_sub( wp->pos[j0], *pi );
      d1  = spu_sub( wp->pos[j1], *pi );
      d2  = spu_sub( wp->pos[j2], *pi );
      d3  = spu_sub( wp->pos[j3], *pi );

      d20 = spu_mul( d0, d0 );
      d21 = spu_mul( d1, d1 );
      d22 = spu_mul( d2, d2 );
      d23 = spu_mul( d3, d3 );

      d20 = spu_add( d20, spu_rlqwbyte( d20, 8 ) );
      d21 = spu_add( d21, spu_rlqwbyte( d21, 8 ) );
      d22 = spu_add( d22, spu_rlqwbyte( d22, 8 ) );
      d23 = spu_add( d23, spu_rlqwbyte( d23, 8 ) );

      d20 = spu_add( d20, spu_rlqwbyte( d20, 4 ) );
      d21 = spu_add( d21, spu_rlqwbyte( d21, 4 ) );
      d22 = spu_add( d22, spu_rlqwbyte( d22, 4 ) );
      d23 = spu_add( d23, spu_rlqwbyte( d23, 4 ) );

      d2a = spu_sel( d20, d21, spu_maskw(7) );
      d2b = spu_sel( d22, d23, spu_maskw(1) );
      r2  = spu_sel( d2a, d2b, spu_maskw(3) );

      /* compute inverses of r2 */
      ms1 = spu_cmpgt( r2cut, r2 );  /* cutoff mask */
      ms2 = spu_cmpgt( r2, f001 );   /* mask zeros */
      r2  = spu_sel( f10, r2, ms2  );
      /* first estimate inverse, then sharpen it */
      r2i = spu_re( r2 ); 
      r2i = spu_mul( r2i, spu_sub( f20, spu_mul( r2i, r2 ) ) ); 

      /* mask unwanted values */
      ms3 = spu_and( ms1, ms2 );
      r2  = spu_sel( f00, r2,  ms3 );
      r2i = spu_sel( f00, r2i, ms3 );

      /* compute LJ interaction */
      tmp2 = spu_mul( r2i,  ljsig );
      tmp3 = spu_mul( tmp2, ljeps );
      tmp4 = spu_mul( tmp2, tmp2 );
      tmp6 = spu_mul( tmp4, tmp2 );
      tmp7 = spu_mul( tmp4, tmp3 );
      pot  = spu_msub( tmp7, spu_sub(tmp6, f20), spu_sel(f00, ljshift, ms3) );
      grad = spu_mul( spu_msub(tmp7, tmp6, tmp7), spu_mul(f12, r2i) );

      /* add up potential energy */
      pots = spu_add( pots, pot );
      pot  = spu_mul( pot,  f05n );  /* avoid double counting */

      /* add to total virial */
      vir  = spu_madd( r2, grad, vir );

      /* the forces */
      ff0  = spu_mul( d0, spu_shuffle( grad, dummy, ss0 ) );
      ff1  = spu_mul( d1, spu_shuffle( grad, dummy, ss1 ) );
      ff2  = spu_mul( d2, spu_shuffle( grad, dummy, ss2 ) );
      ff3  = spu_mul( d3, spu_shuffle( grad, dummy, ss3 ) );

      /* add forces and potential on first particle */
      ff0s = spu_add( ff0s, ff0 );
      ff1s = spu_add( ff1s, ff1 );
      ff2s = spu_add( ff2s, ff2 );
      ff3s = spu_add( ff3s, ff3 );

      /* add forces and potential on second particle */
      wp->force[j0] = spu_sub( wp->force[j0], spu_shuffle( ff0, pot, s0 ) ); 
      wp->force[j1] = spu_sub( wp->force[j1], spu_shuffle( ff1, pot, s1 ) ); 
      wp->force[j2] = spu_sub( wp->force[j2], spu_shuffle( ff2, pot, s2 ) ); 
      wp->force[j3] = spu_sub( wp->force[j3], spu_shuffle( ff3, pot, s3 ) );

    }

    /* add contribution to total poteng */
    pots = spu_add( pots, spu_rlqwbyte( pots, 8 ) );
    pots = spu_add( pots, spu_rlqwbyte( pots, 4 ) );
    wp->totpot += spu_extract( pots, 0 );

    /* add force of first particle */
    ffa = spu_add( ff0s, ff1s );
    ffb = spu_add( ff2s, ff3s );
    *fi = spu_add( *fi,  ffa  );
    *fi = spu_add( *fi,  ffb  );

    /* add potential of first particle */
    *fi = spu_madd( pots, f05l, *fi );

  } 

  /* set contribution to total virial */
  vir = spu_add( vir, spu_rlqwbyte( vir, 8 ) );
  vir = spu_add( vir, spu_rlqwbyte( vir, 4 ) );
  wp->virial = - spu_extract( vir, 0 );

  /* return forces, poteng, totpot, virial */

}

#else

/******************************************************************************
*
*  Lennard-Jones interaction work package for the PPU
*
******************************************************************************/

/* Macro for the calculation of the Lennard-Jones potential */
#define LJ(pot,grad,r2)  {                                               \
   int const col4 = col*4;                                               \
   flt const tmp2 = pt.lj_sig[col4] / r2;                                \
   flt const tmp6 = tmp2 * tmp2 * tmp2;                                  \
   pot  = pt.lj_eps[col4] * tmp6 * (tmp6 - 2.0) - pt.lj_shift[col4];     \
   grad = - 12.0 * pt.lj_eps[col4] * tmp6 * (tmp6 - 1.0) / r2;           \
}

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
