
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

/* 4-point (cubic) or 3-point (quadratic) interpolation of function tables */
#ifdef FOURPOINT
#define   PAIR_INT   PAIR_INT3
#define   VAL_FUNC   VAL_FUNC3
#define DERIV_FUNC DERIV_FUNC3
#else
#define   PAIR_INT   PAIR_INT2
#define   VAL_FUNC   VAL_FUNC2
#define DERIV_FUNC DERIV_FUNC2
#endif

/*****************************************************************************
*
*  Evaluate monolj potential 
*
******************************************************************************/

#define PAIR_INT_MONOLJ(pot, grad, r2)                                        \
{                                                                             \
  real sig_d_rad2, sig_d_rad6, sig_d_rad12;                                   \
                                                                              \
  sig_d_rad2  = 2.0 / (r2);                                                   \
  sig_d_rad6  = sig_d_rad2 * sig_d_rad2 * sig_d_rad2;                         \
  sig_d_rad12 = sig_d_rad6 * sig_d_rad6;                                      \
                                                                              \
  grad = - 6 * sig_d_rad2 * ( sig_d_rad12 - sig_d_rad6 );                     \
  pot  = sig_d_rad12 - 2.0 * sig_d_rad6 - monolj_shift;                       \
}

/*****************************************************************************
*
*  Evaluate potential table with quadratic interpolation. 
*  Returns the potential value and twice the derivative.
*  Note: we need (1/r)(dV/dr) = 2 * dV/dr^2 --> use with equidistant r^2 
*  col is p_typ * ntypes + q_typ
*
******************************************************************************/

#define PAIR_INT2(pot, grad, pt, col, inc, r2, is_short)                     \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;                           \
  int  k;                                                                    \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  k     = (int) (r2a * istep);                                               \
  chi   = (r2a - k * (pt).step[col]) * istep;                                \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k, (col), (pt).maxsteps, (inc));                  \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr;                                                                \
  dv  = p1 - p0;                                                             \
  d2v = p2 - 2 * p1 + p0;                                                    \
                                                                             \
  /* potential and twice the derivative */                                   \
  pot  = p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;                        \
  grad = 2 * istep * (dv + (chi - 0.5) * d2v);                               \
}

/*****************************************************************************
*
*  Evaluate potential table with cubic interpolation. 
*  Returns the potential value and twice the derivative.
*  Note: we need (1/r)(dV/dr) = 2 * dV/dr^2 --> use with equidistant r^2 
*  col is p_typ * ntypes + q_typ
*
******************************************************************************/

#define PAIR_INT3(pot, grad, pt, col, inc, r2, is_short)                     \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, p3, *ptr;                                \
  real fac0, fac1, fac2, fac3, dfac0, dfac1, dfac2, dfac3;                   \
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
  k     = (int) (r2a * istep);                                               \
  chi   = (r2a - k * (pt).step[col]) * istep;                                \
                                                                             \
  /* factors for the interpolation */                                        \
  fac0 = -(1.0/6.0) * chi * (chi-1.0) * (chi-2.0);                           \
  fac1 =        0.5 * (chi*chi-1.0) * (chi-2.0);                             \
  fac2 =       -0.5 * chi * (chi+1.0) * (chi-2.0);                           \
  fac3 =  (1.0/6.0) * chi * (chi*chi-1.0);                                   \
                                                                             \
  /* factors for the interpolation of the derivative */                      \
  dfac0 = -(1.0/6.0) * ((3.0*chi-6.0)*chi+2.0);                              \
  dfac1 =        0.5 * ((3.0*chi-4.0)*chi-1.0);                              \
  dfac2 =       -0.5 * ((3.0*chi-2.0)*chi-2.0);                              \
  dfac3 =    1.0/6.0 * (3.0*chi*chi-1.0);                                    \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k, (col), (pt).maxsteps, (inc));                  \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr; ptr += (inc);                                                  \
  p3  = *ptr;                                                                \
                                                                             \
  /* potential energy */                                                     \
  pot = fac0 * p0 + fac1 * p1 + fac2 * p2 + fac3 * p3;                       \
                                                                             \
  /* twice the derivative */                                                 \
  grad = 2 * istep * (dfac0 * p0 + dfac1 * p1 + dfac2 * p2 + dfac3 * p3);    \
}

/*****************************************************************************
*
*  Evaluate tabulated function with quadratic interpolation. 
*
******************************************************************************/

#define VAL_FUNC2(val, pt, col, inc, r2, is_short)                           \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;                           \
  int  k;                                                                    \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  k     = (int) (r2a * istep);                                               \
  chi   = (r2a - k * (pt).step[col]) * istep;                                \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k, (col), (pt).maxsteps, (inc));                  \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr;                                                                \
  dv  = p1 - p0;                                                             \
  d2v = p2 - 2 * p1 + p0;                                                    \
                                                                             \
  /* the function value */                                                   \
  val = p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;                         \
}

/*****************************************************************************
*
*  Evaluate tabulated function with cubic interpolation. 
*
******************************************************************************/

#define VAL_FUNC3(val, pt, col, inc, r2, is_short)                           \
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
  k     = (int) (r2a * istep);                                               \
  chi   = (r2a - k * (pt).step[col]) * istep;                                \
                                                                             \
  /* factors for the interpolation */                                        \
  fac0 = -(1.0/6.0) * chi * (chi-1.0) * (chi-2.0);                           \
  fac1 =        0.5 * (chi*chi-1.0) * (chi-2.0);                             \
  fac2 =       -0.5 * chi * (chi+1.0) * (chi-2.0);                           \
  fac3 =  (1.0/6.0) * chi * (chi*chi-1.0);                                   \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k, (col), (pt).maxsteps, (inc));                  \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr; ptr += (inc);                                                  \
  p3  = *ptr;                                                                \
                                                                             \
  /* the function value */                                                   \
  val = fac0 * p0 + fac1 * p1 + fac2 * p2 + fac3 * p3;                       \
}

/*****************************************************************************
*
*  Evaluate the derivative of a function with quadratic interpolation. 
*  Returns *twice* the derivative.
*
******************************************************************************/

#define DERIV_FUNC2(grad, pt, col, inc, r2, is_short)                        \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;                           \
  int  k;                                                                    \
                                                                             \
  /* check for distances shorter than minimal distance in table */           \
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a   = 0;                                                               \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  k     = (int) (r2a * istep);                                               \
  chi   = (r2a - k * (pt).step[col]) * istep;                                \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k, (col), (pt).maxsteps, (inc));                  \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr;                                                                \
  dv  = p1 - p0;                                                             \
  d2v = p2 - 2 * p1 + p0;                                                    \
                                                                             \
  /* twice the derivative */                                                 \
  grad = 2 * istep * (dv + (chi - 0.5) * d2v);                               \
}

/*****************************************************************************
*
*  Evaluate the derivative of a function with cubic interpolation. 
*  Returns *twice* the derivative.
*
******************************************************************************/

#define DERIV_FUNC3(grad, pt, col, inc, r2, is_short)                        \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, p3, *ptr;                                \
  real dfac0, dfac1, dfac2, dfac3;                                           \
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
  k     = (int) (r2a * istep);                                               \
  chi   = (r2a - k * (pt).step[col]) * istep;                                \
                                                                             \
  /* factors for the interpolation of the 1. derivative */                   \
  dfac0 = -(1.0/6.0) * ((3.0*chi-6.0)*chi+2.0);                              \
  dfac1 =        0.5 * ((3.0*chi-4.0)*chi-1.0);                              \
  dfac2 =       -0.5 * ((3.0*chi-2.0)*chi-2.0);                              \
  dfac3 =    1.0/6.0 * (3.0*chi*chi-1.0);                                    \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k, (col), (pt).maxsteps, (inc));                  \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr; ptr += (inc);                                                  \
  p3  = *ptr;                                                                \
                                                                             \
  /* twice the derivative */                                                 \
  grad = 2 * istep * (dfac0 * p0 + dfac1 * p1 + dfac2 * p2 + dfac3 * p3);    \
}
