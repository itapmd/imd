
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2006 Institute for Theoretical and Applied Physics,
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

/* 4-point (cubic), spline, or 3-point (quadratic) interpolation of 
   function tables */
#if   defined(FOURPOINT)
#define   PAIR_INT   PAIR_INT3
#define   VAL_FUNC   VAL_FUNC3
#define DERIV_FUNC DERIV_FUNC3
#elif defined(SPLINE)
#define   PAIR_INT   PAIR_INT_SP
#define   VAL_FUNC   VAL_FUNC_SP
#define DERIV_FUNC DERIV_FUNC_SP
#else
#define   PAIR_INT   PAIR_INT2
#define   VAL_FUNC   VAL_FUNC2
#define DERIV_FUNC DERIV_FUNC2
#endif

/* compensate for non-standard cast in icc (which is faster) */
/* this works only for a positive argument */
#ifdef RCD
#define POS_TRUNC(x) ((int)(-0.5+(x)))
#else
#define POS_TRUNC(x) ((int)(x))
#endif

/*****************************************************************************
*
*  Evaluate harmonic crystal
*
******************************************************************************/

#define ATOM_INT_EC(pot, grad, p_typ, r2)		        	\
{									\
  pot += 0.5 * spring_rate[p_typ] * r2;			        	\
  grad += spring_rate[p_typ];				        	\
}

/*****************************************************************************
*
*  Evaluate Lennard-Jones potential
*
******************************************************************************/

#define PAIR_INT_LJ(pot, grad, p_typ, q_typ, r2)			\
{									\
  real sig_d_rad2, sig_d_rad6, sig_d_rad12;				\
									\
  sig_d_rad2  = lj_sigma[p_typ][q_typ] * lj_sigma[p_typ][q_typ] / (r2);	\
  sig_d_rad6  = sig_d_rad2 * sig_d_rad2 * sig_d_rad2;			\
  sig_d_rad12 = sig_d_rad6 * sig_d_rad6;				\
									\
  pot += lj_epsilon[p_typ][q_typ] * ( sig_d_rad12 - 2.0 * sig_d_rad6 )	\
    - lj_shift[p_typ][q_typ];						\
  grad += - 12.0 * lj_epsilon[p_typ][q_typ] / (r2)			\
    * ( sig_d_rad12 - sig_d_rad6 );					\
}

#ifdef VEC

#define PAIR_INT_LJ_VEC(pot, grad, col, r2)		                   \
{									   \
  real sig_d_rad2, sig_d_rad6, sig_d_rad12;				   \
									   \
  sig_d_rad2  = lj_sigma2_vec[col] / (r2);                                 \
  sig_d_rad6  = sig_d_rad2 * sig_d_rad2 * sig_d_rad2;			   \
  sig_d_rad12 = sig_d_rad6 * sig_d_rad6;				   \
									   \
  pot  = lj_epsilon_vec[col] * (sig_d_rad12 - 2.0 * sig_d_rad6);           \
  grad = - 12.0 * lj_epsilon_vec[col] * (sig_d_rad12 - sig_d_rad6) / (r2); \
}

/*****************************************************************************
*
*  Evaluate Lennard-Jones potential
*
******************************************************************************/

#define PAIR_INT_LJG(pot, grad, p_typ, q_typ, r2)			\
{									\
  real sig_d_rad2, sig_d_rad6, sig_d_rad12, expo, dr;			\
									\
  sig_d_rad2  = lj_sigma[p_typ][q_typ] * lj_sigma[p_typ][q_typ] / (r2);	\
  sig_d_rad6  = sig_d_rad2 * sig_d_rad2 * sig_d_rad2;			\
  sig_d_rad12 = sig_d_rad6 * sig_d_rad6;				\
                                                                        \
  dr = ( sqrt(r2)-ljg_r0[p_typ][q_typ] ) / ljg_sig[p_typ][q_typ];       \
  expo = exp ( - 0.5 * dr * dr );					\
									\
  pot += lj_epsilon[p_typ][q_typ] * ( sig_d_rad12 - 2.0 * sig_d_rad6 )	\
    - lj_shift[p_typ][q_typ]                                            \
    - ljg_eps[p_typ][q_typ] * expo;					\
  grad += - 12.0 * lj_epsilon[p_typ][q_typ] / (r2)			\
    * ( sig_d_rad12 - sig_d_rad6 )                                      \
    - ljg_eps[p_typ][q_typ] * dr * expo / ljg_sig[p_typ][q_typ];	\

/*****************************************************************************
*
*  Evaluate potential table with quadratic interpolation. 
*  Returns the potential value and twice the derivative.
*  Note: we need (1/r)(dV/dr) = 2 * dV/dr^2 --> use with equidistant r^2 
*  col is p_typ * ntypes + q_typ
*
******************************************************************************/

#define PAIR_INT_VEC(pot, grd, pt, col, inc, r2)                             \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;                           \
  long kk;                                                                   \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2 * istep;                                                        \
  kk    = (long) (r2a);                                                      \
  chi   = r2a - kk;                                                          \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, kk, (col), (pt).maxsteps, (inc));                 \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr;                                                                \
  dv  = p1 - p0;                                                             \
  d2v = p2 - 2 * p1 + p0;                                                    \
                                                                             \
  /* potential and twice the derivative */                                   \
  pot = p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;                         \
  grd = 2 * istep * (dv + (chi - 0.5) * d2v);                                \
}

#define VAL_FUNC_VEC(pot, pt, col, inc, r2)                                  \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;                           \
  long kk;                                                                   \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2 * istep;                                                        \
  kk    = (long) (r2a);                                                      \
  chi   = r2a - kk;                                                          \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, kk, (col), (pt).maxsteps, (inc));                 \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr;                                                                \
  dv  = p1 - p0;                                                             \
  d2v = p2 - 2 * p1 + p0;                                                    \
                                                                             \
  /* potential value */                                                      \
  pot = p0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;                         \
}

#define DERIV_FUNC_VEC(grd, pt, col, inc, r2)                                \
{                                                                            \
  real r2a, istep, chi, p0, p1, p2, dv, d2v, *ptr;                           \
  long kk;                                                                   \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2 * istep;                                                        \
  kk    = (long) (r2a);                                                      \
  chi   = r2a - kk;                                                          \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, kk, (col), (pt).maxsteps, (inc));                 \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr;                                                                \
  dv  = p1 - p0;                                                             \
  d2v = p2 - 2 * p1 + p0;                                                    \
                                                                             \
  /* twice the derivative */                                                 \
  grd = 2 * istep * (dv + (chi - 0.5) * d2v);                                \
}

#endif

/*****************************************************************************
*
*  Evaluate Morse potential
*
******************************************************************************/

#define PAIR_INT_MORSE(pot, grad, p_typ, q_typ, r2)			    \
{									    \
  real r, exppot, cexppot;						    \
									    \
  r       = sqrt((r2));							    \
  exppot  = exp( - morse_alpha[p_typ][q_typ]				    \
		 * ( r - morse_sigma[p_typ][q_typ] ) );			    \
  cexppot = 1.0 - exppot;						    \
									    \
  pot  += morse_epsilon[p_typ][q_typ] * ( cexppot * cexppot - 1.0 )	    \
    - morse_shift[p_typ][q_typ];					    \
  grad += 2.0 * morse_alpha[p_typ][q_typ] * morse_epsilon[p_typ][q_typ] / r \
    * exppot * cexppot;							    \
}

/*****************************************************************************
*
*  Evaluate Buckingham potential 
*
******************************************************************************/

#define PAIR_INT_BUCK(pot, grad, p_typ, q_typ, r2)			   \
{									   \
  real rinv, rinv2, powpot, exppot, invs2;				   \
									   \
  rinv   = buck_sigma[p_typ][q_typ] / sqrt((r2));			   \
  rinv2  = rinv * rinv;							   \
  powpot  = buck_c[p_typ][q_typ] * rinv2 * rinv2 * rinv2;		   \
  exppot = buck_a[p_typ][q_typ] * exp ( - 1.0 / rinv );			   \
  invs2   = 1.0 / ( buck_sigma[p_typ][q_typ] * buck_sigma[p_typ][q_typ] ); \
									   \
  pot += exppot - powpot - buck_shift[p_typ][q_typ];			   \
  grad += ( - exppot * rinv + 6 * powpot * rinv2 ) * invs2;		   \
}

/*****************************************************************************
*
*  Evaluate pair potential for Keating potential 
*
******************************************************************************/

#define PAIR_INT_KEATING(pot, grad, p_typ, q_typ, r2)                         \
{                                                                             \
  real tmp_pot, tmp;                                                          \
                                                                              \
  tmp_pot = 3.0 * keat_alpha[p_typ][q_typ]                                    \
          / ( 8.0 * keat_d[p_typ][q_typ] * keat_d[p_typ][q_typ] );            \
  tmp = r2 - keat_d[p_typ][q_typ] * keat_d[p_typ][q_typ];                     \
                                                                              \
  pot = tmp_pot * tmp * tmp;                                                  \
  grad = 4.0 * tmp_pot * tmp;                                                 \
}

/*****************************************************************************
*
*  Evaluate pair potential for Stillinger-Weber potential 
*
******************************************************************************/

#define PAIR_INT_STIWEB(pot, grad, p_typ, q_typ, r2)                          \
{                                                                             \
  real radius, phi_r, phi_a, inv_c, inv_r, f_cut;                             \
                                                                              \
  radius = sqrt((r2));                                                        \
  phi_r  = sw_a[p_typ][q_typ] * pow( radius, - sw_p[p_typ][q_typ] );          \
  phi_a  = - sw_b[p_typ][q_typ] * pow( radius, - sw_q[p_typ][q_typ] );        \
  inv_c  = 1.0 / ( (radius) - sw_a1[p_typ][q_typ] );                          \
  inv_r  = 1.0 / (radius);                                                    \
  f_cut  = exp( sw_de[p_typ][q_typ] * inv_c );                                \
                                                                              \
  pot  = ( phi_r + phi_a ) * f_cut;                                           \
  grad = ( - pot * sw_de[p_typ][q_typ] * inv_c * inv_c                        \
         - f_cut * inv_r * ( sw_p[p_typ][q_typ] * phi_r                       \
                           + sw_q[p_typ][q_typ] * phi_a ) ) * inv_r;          \
                                                                              \
}


/*****************************************************************************
*
*  Evaluate potential table with linear interpolation. 
*  Returns the potential value and twice the derivative.
*  Note: we need (1/r)(dV/dr) = 2 * dV/dr^2 --> use with equidistant r^2 
*  col is p_typ * ntypes + q_typ
*
******************************************************************************/

#define PAIR_INT_LIN(pot, grad, pt, col, inc, r2, is_short)                  \
{                                                                            \
  real r2a, chi, *t;                                                         \
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
  r2a   = r2a *  (pt).invstep[col];                                          \
  k     = POS_TRUNC(r2a);                                                    \
  chi   = r2a - k;                                                           \
                                                                             \
  /* potential and twice the derivative */                                   \
  t    = (pt).table[col];                                                    \
  pot  = t[2*k  ]*(1.0-chi) + t[2*k+2]*chi;                                  \
  grad = t[2*k+1]*(1.0-chi) + t[2*k+3]*chi;                                  \
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
  r2a   = r2a * istep;                                                       \
  k     = POS_TRUNC(r2a);                                                    \
  /* k     = MIN( POS_TRUNC(r2a), (pt).len[col]-3 ); */                      \
  chi   = r2a - k;                                                           \
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
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2a * istep;                                                       \
  k     = MAX( POS_TRUNC(r2a), 1 );                                          \
  /* k     = MIN( MAX( POS_TRUNC(r2a), 1 ), (pt).len[col]-3 ); */            \
  chi   = r2a - k;                                                           \
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
  ptr = PTR_2D((pt).table, k-1, (col), (pt).maxsteps, (inc));                \
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
*  Evaluate potential table with spline interpolation. 
*  Returns the potential value and twice the derivative.
*  Note: we need (1/r)(dV/dr) = 2 * dV/dr^2 --> use with equidistant r^2 
*  col is p_typ * ntypes + q_typ
*
******************************************************************************/

#define PAIR_INT_SP(pot, grad, pt, col, inc, r2, is_short)                   \
{                                                                            \
  real r2a, a, b, a2, b2, istep, step, st6, p1, p2, d21, d22;                \
  int k;                                                                     \
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
  step  = (pt).step[col];                                                    \
  r2a   = r2a * istep;                                                       \
  k     = POS_TRUNC(r2a);                                                    \
  /* k     = MIN( POS_TRUNC(r2a), (pt).len[col]-2 ); */                      \
  b     = r2a - k;                                                           \
  a     = 1.0 - b;                                                           \
                                                                             \
  /* intermediate values */                                                  \
  k     = k * (inc) + (col);                                                 \
  p1    = (pt).table[k];                                                     \
  d21   = (pt).table2[k];                                                    \
  k    += (inc);                                                             \
  p2    = (pt).table[k];                                                     \
  d22   = (pt).table2[k];                                                    \
  a2    = a * a - 1;                                                         \
  b2    = b * b - 1;                                                         \
  st6   = step / 6;                                                          \
                                                                             \
  /* potential and twice the derivative */                                   \
  pot  = a * p1 + b * p2 + (a * a2 * d21 + b * b2 * d22) * st6 * step;       \
  grad = 2*((p2 - p1) * istep + ((3*b2 + 2) * d22 - (3*a2 + 2) * d21) * st6);\
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
  r2a   = r2a * istep;                                                       \
  k     = POS_TRUNC(r2a);                                                    \
  /* k     = MIN( POS_TRUNC(r2a), (pt).len[col]-3 ); */                      \
  chi   = r2a - k;                                                           \
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
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2a * istep;                                                       \
  k     = MAX( POS_TRUNC(r2a), 1 );                                          \
  /* k     = MIN( MAX( POS_TRUNC(r2a), 1 ), (pt).len[col]-3 ); */            \
  chi   = r2a - k;                                                           \
                                                                             \
  /* factors for the interpolation */                                        \
  fac0 = -(1.0/6.0) * chi * (chi-1.0) * (chi-2.0);                           \
  fac1 =        0.5 * (chi*chi-1.0) * (chi-2.0);                             \
  fac2 =       -0.5 * chi * (chi+1.0) * (chi-2.0);                           \
  fac3 =  (1.0/6.0) * chi * (chi*chi-1.0);                                   \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k-1, (col), (pt).maxsteps, (inc));                \
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
*  Evaluate tabulated function with spline interpolation. 
*
******************************************************************************/

#define VAL_FUNC_SP(val, pt, col, inc, r2, is_short)                         \
{                                                                            \
  real r2a, a, b, a2, b2, istep, step, st6, p1, p2, d21, d22;                \
  int k;                                                                     \
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
  step  = (pt).step[col];                                                    \
  r2a   = r2a * istep;                                                       \
  k     = POS_TRUNC(r2a);                                                    \
  /* k     = MIN( POS_TRUNC(r2a), (pt).len[col]-2 ); */                      \
  b     = r2a - k;                                                           \
  a     = 1.0 - b;                                                           \
                                                                             \
  /* intermediate values */                                                  \
  k     = k * (inc) + (col);                                                 \
  p1    = (pt).table[k];                                                     \
  d21   = (pt).table2[k];                                                    \
  k    += (inc);                                                             \
  p2    = (pt).table[k];                                                     \
  d22   = (pt).table2[k];                                                    \
  a2    = a * a - 1;                                                         \
  b2    = b * b - 1;                                                         \
  st6   = step / 6;                                                          \
                                                                             \
  /* the function value */                                                   \
  val  = a * p1 + b * p2 + (a * a2 * d21 + b * b2 * d22) * st6 * step;       \
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
  r2a   = r2a * istep;                                                       \
  k     = POS_TRUNC(r2a);                                                    \
  /* k     = MIN( POS_TRUNC(r2a), (pt).len[col]-3 ); */                      \
  chi   = r2a - k;                                                           \
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
  r2a = MIN((r2),(pt).end[col]);                                             \
  r2a = r2a - (pt).begin[col];                                               \
  if (r2a < 0) {                                                             \
    r2a = 0;                                                                 \
    is_short = 1;                                                            \
  }                                                                          \
                                                                             \
  /* indices into potential table */                                         \
  istep = (pt).invstep[col];                                                 \
  r2a   = r2a * istep;                                                       \
  k     = MAX( POS_TRUNC(r2a), 1 );                                          \
  /* k     = MIN( MAX( POS_TRUNC(r2a), 1 ), (pt).len[col]-3 ); */            \
  chi   = r2a - k;                                                           \
                                                                             \
  /* factors for the interpolation of the 1. derivative */                   \
  dfac0 = -(1.0/6.0) * ((3.0*chi-6.0)*chi+2.0);                              \
  dfac1 =        0.5 * ((3.0*chi-4.0)*chi-1.0);                              \
  dfac2 =       -0.5 * ((3.0*chi-2.0)*chi-2.0);                              \
  dfac3 =    1.0/6.0 * (3.0*chi*chi-1.0);                                    \
                                                                             \
  /* intermediate values */                                                  \
  ptr = PTR_2D((pt).table, k-1, (col), (pt).maxsteps, (inc));                \
  p0  = *ptr; ptr += (inc);                                                  \
  p1  = *ptr; ptr += (inc);                                                  \
  p2  = *ptr; ptr += (inc);                                                  \
  p3  = *ptr;                                                                \
                                                                             \
  /* twice the derivative */                                                 \
  grad = 2 * istep * (dfac0 * p0 + dfac1 * p1 + dfac2 * p2 + dfac3 * p3);    \
}

/*****************************************************************************
*
*  Evaluate the derivative of a function with spline interpolation. 
*  Returns *twice* the derivative.
*
******************************************************************************/

#define DERIV_FUNC_SP(grad, pt, col, inc, r2, is_short)                      \
{                                                                            \
  real r2a, a, b, a2, b2, istep, step, st6, p1, p2, d21, d22;                \
  int k;                                                                     \
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
  step  = (pt).step[col];                                                    \
  r2a   = r2a * istep;                                                       \
  k     = POS_TRUNC(r2a);                                                    \
  /* k     = MIN( POS_TRUNC(r2a), (pt).len[col]-2 ); */                      \
  b     = r2a - k;                                                           \
  a     = 1.0 - b;                                                           \
                                                                             \
  /* intermediate values */                                                  \
  k     = k * (inc) + (col);                                                 \
  p1    = (pt).table[k];                                                     \
  d21   = (pt).table2[k];                                                    \
  k    += (inc);                                                             \
  p2    = (pt).table[k];                                                     \
  d22   = (pt).table2[k];                                                    \
  a2    = a * a - 1;                                                         \
  b2    = b * b - 1;                                                         \
  st6   = step / 6;                                                          \
                                                                             \
  /* twice the derivative */                                                 \
  grad = 2*((p2 - p1) * istep + ((3*b2 + 2) * d22 - (3*a2 + 2) * d21) * st6);\
}

