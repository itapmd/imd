/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/*****************************************************************************
*
*  Macros for potential table access. We use macros in order to enforce
*  inlining, which is hard to achieve with functions. We need inlining
*  because these access functions are called very frequently.
*
*  $Revision$
*  $Date$
*
******************************************************************************/

#define POS_TRUNC(x) ((int)(x))

/*****************************************************************************
*
*  Evaluate potential table with quartic interpolation. 
*  Note: we need (1/r)(dV/dr) = 2 * dV/dr^2 --> use with equidistant r^2 
*  col is p_typ * ntypes + q_typ
*
******************************************************************************/

#define PAIR_INT4(phi, dphi, ddphi, dddphi, pt, col, inc, r2, is_short)      \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, p3, p4, *ptr;                            \
  real fac0, fac1, fac2, fac3, fac4, dfac0, dfac1, dfac2, dfac3, dfac4;      \
  real ddfac0, ddfac1, ddfac2, ddfac3, ddfac4;				     \
  real dddfac0, dddfac1, dddfac2, dddfac3, dddfac4;                          \
  int  k;                                                                    \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  /* we need one extra value at the lower end for interpolation */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col] - (pt).step[col];                              \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  k     = POS_TRUNC(r2a * istep);                                            \
  chi   = (r2a - k * (pt).step[col]) * istep;                                \
                                                                             \
  /* factors for the interpolation */                                        \
  fac0 =    1.0/24.0 * chi * (chi-1.0) * (chi-2.0) * (chi-3.0);              \
  fac1 =  -(1.0/6.0) * (chi*chi-1.0) * (chi-2.0) * (chi-3.0);                \
  fac2 =        0.25 * chi * (chi+1.0) * (chi-2.0) * (chi-3.0);              \
  fac3 =  -(1.0/6.0) * chi * (chi*chi-1.0) * (chi-3.0);                      \
  fac4 =  (1.0/24.0) * chi * (chi*chi-1.0) * (chi-2.0);                      \
                                                                             \
  /* factors for the interpolation of the first derivative */                \
  dfac0 =   1.0/24.0 * ( (2.0*chi-1.0)*(chi-2.0)*(chi-3.0)                   \
      		       + chi*(chi-1.0)*(2.0*chi-5.0) );                      \
  dfac1 = -(1.0/6.0) * ( 2.0*chi*(chi-2.0)*(chi-3.0)                         \
			 + (chi*chi-1.0)*(2.0*chi-5.0) );                    \
  dfac2 =        0.25 * ( (2.0*chi+1.0)*(chi-2.0)*(chi-3.0)                  \
  		       + chi*(chi+1.0)*(2.0*chi-5.0) );                      \
  dfac3 =  -(1.0/6.0) * ( (3.0*chi*chi-1.0)*(chi-3.0)                        \
       		       + chi*(chi*chi-1.0) );                                \
  dfac4 =    1.0/24.0 * ( (3.0*chi*chi-1.0)*(chi-2.0)                        \
       		       + chi*(chi*chi-1.0) );                                \
                                                                             \
  /* factors for the interpolation of the second derivative */               \
  ddfac0 = (1.0/12.0) * ( 6.0*chi*chi - 18.0*chi + 11.0 );                   \
  ddfac1 = -(1.0/3.0) * ( 6.0*chi*chi - 15.0*chi + 5.0 );                    \
  ddfac2 = (1.0/2.0)  * ( 6.0*chi*chi - 12.0*chi + 1.0 );                    \
  ddfac3 = -(1.0/3.0) * ( 6.0*chi*chi - 9.0*chi - 1.0 );                     \
  ddfac4 = (1.0/12.0) * ( 6.0*chi*chi - 6.0*chi - 1.0 );                     \
                                                                             \
  /* factors for the interpolation of the third derivative */                \
  dddfac0 = (1.0/2.0) * ( 2.0 * chi - 3.0 );                                 \
  dddfac1 =             - 4.0 * chi + 5.0;                                   \
  dddfac2 =       6.0  * ( chi - 1.0 );                                      \
  dddfac3 =             - 4.0 * chi + 3.0;                                   \
  dddfac4 = (1.0/2.0)  * ( 2.0 * chi - 1.0 );                                \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k, (col), (inc));                                 \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr; ptr += (inc);                                                  \
  p3  = *ptr; ptr += (inc);                                                  \
  p4  = *ptr;                                                                \
                                                                             \
  /* potential energy */                                                     \
  phi = fac0 * p0 + fac1 * p1 + fac2 * p2 + fac3 * p3 + fac4 * p4;           \
	                                                                     \
  /* first derivative, dV/d(r^2) */                                          \
  dphi = istep *                                                             \
    (dfac0 * p0 + dfac1 * p1 + dfac2 * p2 + dfac3 * p3 + dfac4 * p4);        \
                                                                             \
  /* second derivative, d^2V/d(r^2)^2 */                                     \
  ddphi = istep * istep *                                                    \
    (ddfac0 * p0 + ddfac1 * p1 + ddfac2 * p2 + ddfac3 * p3 + ddfac4 * p4);   \
                                                                             \
  /* third derivative d^3V/d(r^2)^3 */                                       \
  dddphi = istep * istep * istep *                                           \
    (dddfac0 * p0 + dddfac1 * p1 + dddfac2 * p2                              \
    + dddfac3 * p3 + dddfac4 * p4);                                          \
									     \
}

/*****************************************************************************
*
*  Evaluate the derivatives of a function with quartic interpolation. 
*
******************************************************************************/

#define DERIV_FUNC(dphi, ddphi, dddphi, pt, col, inc, r2, is_short)          \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, p3, p4, *ptr;                            \
  real dfac0, dfac1, dfac2, dfac3, dfac4;                                    \
  real ddfac0, ddfac1, ddfac2, ddfac3, ddfac4;				     \
  real dddfac0, dddfac1, dddfac2, dddfac3, dddfac4;                          \
  int  k;                                                                    \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  /* we need one extra value at the lower end for interpolation */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col] - (pt).step[col];                              \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  k     = POS_TRUNC(r2a * istep);                                            \
  chi   = (r2a - k * (pt).step[col]) * istep;                                \
                                                                             \
  /* factors for the interpolation of the first derivative */                \
  dfac0 =  (1.0/24.0) * ( (2.0*chi-1.0)*(chi-2.0)*(chi-3.0)                  \
      		       + chi*(chi-1.0)*(2.0*chi-5.0) );                      \
  dfac1 =  -(1.0/6.0) * ( 2.0*chi*(chi-2.0)*(chi-3.0)                        \
			 + (chi*chi-1.0)*(2.0*chi-5.0) );                    \
  dfac2 =        0.25 * ( (2.0*chi+1.0)*(chi-2.0)*(chi-3.0)                  \
  		       + chi*(chi+1.0)*(2.0*chi-5.0) );                      \
  dfac3 =  -(1.0/6.0) * ( (3.0*chi*chi-1.0)*(chi-3.0)                        \
       		       + chi*(chi*chi-1.0) );                                \
  dfac4 =    1.0/24.0 * ( (3.0*chi*chi-1.0)*(chi-2.0)                        \
       		       + chi*(chi*chi-1.0) );                                \
                                                                             \
  /* factors for the interpolation of the second derivative */               \
  ddfac0 = (1.0/12.0) * ( 6.0*chi*chi - 18.0*chi + 11.0 );                   \
  ddfac1 = -(1.0/3.0) * ( 6.0*chi*chi - 15.0*chi + 5.0 );                    \
  ddfac2 = (1.0/2.0)  * ( 6.0*chi*chi - 12.0*chi + 1.0 );                    \
  ddfac3 = -(1.0/3.0) * ( 6.0*chi*chi - 9.0*chi - 1.0 );                     \
  ddfac4 = (1.0/12.0) * ( 6.0*chi*chi - 6.0*chi - 1.0 );                     \
  /* factors for the interpolation of the third derivative */                \
  dddfac0 = (1.0/2.0) * ( 2.0*chi - 3.0 );                                   \
  dddfac1 =             - 4.0*chi + 5.0;                                     \
  dddfac2 =       6.0  * ( chi - 1.0 );                                      \
  dddfac3 =             - 4.0*chi + 3.0;                                     \
  dddfac4 = (1.0/2.0)  * ( 2.0*chi - 1.0 );                                  \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k, (col), (inc));                                 \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr; ptr += (inc);                                                  \
  p3  = *ptr; ptr += (inc);                                                  \
  p4  = *ptr;                                                                \
	                                                                     \
  /* first derivative, dV/d(r^2) */                                          \
  dphi = istep *                                                             \
    (dfac0 * p0 + dfac1 * p1 + dfac2 * p2 + dfac3 * p3 + dfac4 * p4);        \
                                                                             \
  /* second derivative, d^2V/d(r^2)^2 */                                     \
  ddphi = istep * istep *                                                    \
    (ddfac0 * p0 + ddfac1 * p1 + ddfac2 * p2 + ddfac3 * p3 + ddfac4 * p4);   \
                                                                             \
  /* third derivative d^3V/d(r^2)^3 */                                       \
  dddphi = istep * istep * istep *                                           \
    (dddfac0 * p0 + dddfac1 * p1 + dddfac2 * p2                              \
    + dddfac3 * p3 + dddfac4 * p4);                                          \
}

#ifdef EAM

/*****************************************************************************
*
*  Evaluate tabulated function with cubic interpolation. 
*
******************************************************************************/

#define VAL_FUNC(val, pt, col, inc, r2, is_short)                            \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, p3;                                      \
  real fac0, fac1, fac2, fac3, *ptr;                                         \
  int  k;                                                                    \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  /* we need one extra value at the lower end for interpolation */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col] - (pt).step[col];                              \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  k     = POS_TRUNC(r2a * istep);                                            \
  chi   = (r2a - k * (pt).step[col]) * istep;                                \
                                                                             \
  /* factors for the interpolation */                                        \
  fac0 = -(1.0/6.0) * chi * (chi-1.0) * (chi-2.0);                           \
  fac1 =        0.5 * (chi*chi-1.0) * (chi-2.0);                             \
  fac2 =       -0.5 * chi * (chi+1.0) * (chi-2.0);                           \
  fac3 =  (1.0/6.0) * chi * (chi*chi-1.0);                                   \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k, (col), (inc));                                 \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr; ptr += (inc);                                                  \
  p3  = *ptr;                                                                \
                                                                             \
  /* the function value */                                                   \
  val = fac0 * p0 + fac1 * p1 + fac2 * p2 + fac3 * p3;                       \
}

#endif
