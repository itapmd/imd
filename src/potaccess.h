
/*****************************************************************************
*
*  Macros for potential table access. We use macros in order to enforce
*  inlining, which hard to achieve with functions. We need inlining
*  because these access functions are called very frequently.
*
*  $Revision$
*  $Date$
*
******************************************************************************/


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
*  Returned are potential value V and gradient (dV/dr)/r.
*  col is p_typ * ntypes + q_typ
*  pot, grad, r2_short and is_short must be variables, not expressions.
*
******************************************************************************/

#define PAIR_INT2(pot, grad, pt, col, inc, r2, r2_short, is_short)          \
{                                                                           \
  real r2a, istep, chi, p0, p1, p2, dv, d2v;                                \
  real *ptr;                                                                \
  int  k;                                                                   \
                                                                            \
  /* check for distances shorter than minimal distance in table */          \
  r2a = (r2) - (pt).begin[col];                                             \
  if (r2a < 0) {                                                            \
    r2_short = MIN(r2_short,(r2));                                          \
    r2a = 0;                                                                \
    is_short = 1;                                                           \
  }                                                                         \
                                                                            \
  /* Indices into potential table */                                        \
  istep = (pt).invstep[col];                                                \
  k     = (int) (r2a * istep);                                              \
  chi   = (r2a - k * (pt).step[col]) * istep;                               \
                                                                            \
  /* make potential table access as efficient as possible - use pointers */ \
  ptr = PTR_2D((pt).table, k, (col), (pt).maxsteps, (inc));                 \
  p0  = *ptr; ptr += (inc);                                                 \
  p1  = *ptr; ptr += (inc);                                                 \
  p2  = *ptr;                                                               \
  dv  = p1 - p0;                                                            \
  d2v = p2 - 2 * p1 + p0;                                                   \
                                                                            \
  /* norm of gradient and potential energy */                               \
  grad = 2 * istep * (dv + (chi - 0.5) * d2v);                              \
  pot  = p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;                       \
}
