
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2006 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* Routines for CG integrator - adapted from Numerical Recipies
*
* $Revision$
* $Date$
*
******************************************************************************/

#include <math.h>
// these are now parameters:
/* #define ITMAX 100 */
//#define ZEPS 1.0e-10
//#define GLIMIT 100.0  /* was  100.0 */
#define CGOLD 0.3819660
#define GOLD 1.618034
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? FABS(a) : -FABS(a))

/* brent from Num. Rec. modified for our cg */
/* fb should be the minimal value */

int brent(real ax, real bx, real cx, real fa, real fb, real fc, real *alphamin)
{
  int iter,ITMAX;
  real a, b, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
  real e = 0.0, d = 0.0, tol = linmin_tol;
  ITMAX=linmin_maxsteps;

 
  /* sorting is already done in mnbrak */
  /*   a = (ax < cx) ? ax : cx; */
  /*   b = (ax > cx) ? ax : cx; */
  a=ax; b=bx;
  x = w = v = bx;
  
  // fb = fonedim(b); 
 
  fw = fv = fx = fb;
  if ((cg_infolevel>0) && (0==myid)) {
      printf("in brent starting with  a %.12e b %.12e fb %.12e\n",a,b,fb);
      fflush(stdout);
  }

  for (iter = 1; iter <= ITMAX; iter++) {
    xm = 0.5 * (a+b);
    tol1 = tol * FABS(x) + cg_zeps;
    tol2 = 2.0 * tol1;

    if ((cg_infolevel>0) && (0==myid)) {
      printf("in brent, iter %d x %e v %e w %e fx %e fv %e fw %e\n",
             iter, x,v,w,fx,fv,fw);
      fflush(stdout);
    }

    if (FABS(x-xm) <= (tol2-0.5*(b-a))) {
      *alphamin = x;
      return iter;
    }
    if (FABS(e) > tol1) {
      r = (x-w) * (fx-fv);
      q = (x-v) * (fx-fw);
      p = (x-v)*q - (x-w)*r;
      q = 2.0 * (q-r);
      if (q > 0.0) p = -p;
      q = FABS(q);
      etemp = e;
      e = d;
      if (FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) {
        e = (x >= xm) ? (a-x) : (b-x);
        d = CGOLD * e;
      }
      else {
        d = p/q;
        u = x+d;
        if ((u-a < tol2) || (b-u < tol2)) d = SIGN(tol1,xm-x);
      }
    } 
    else {
      e = (x >= xm) ? (a-x) : (b-x);
      d = CGOLD * e;
    }
    u  = ((FABS(d) >= tol1) ? (x+d) : (x+SIGN(tol1,d)));
    fu = fonedim(u);
    if (fu <= fx) {
      if (u >= x) a = x; else b = x;
      SHFT( v, w, x, u)
      SHFT(fv,fw,fx,fu)
    } 
    else {
      if (u < x) a = u; else b = u;
      if ((fu <= fw) || (w == x)) {
        v = w; fv = fw;
        w = u; fw = fu;
      } 
      else if ((fu <= fv) || (v == x) || (v == w)) {
        v = u; fv = fu;
      }
    }
  }
  error("Too many iterations in brent");
  return 0;
}


/* bracketing of minimum from Num Rec. */
 
int mnbrak(real *ax, real *bx, real *cx, real *fa, real *fb, real *fc)
{
  real ulim, u, r, q, fu, dum;
  int  ctf;

  /* fb is calculated in linmin, fa is the old value at alpha_a=0 */
  if (*fb > *fa) {
    SHFT(dum, *ax, *bx, dum)
    SHFT(dum, *fb, *fa, dum)
  }

  *cx = (*bx) + GOLD * (*bx-*ax);
  if( (*bx-*ax) <= 1e-4) {
    *cx = (*bx) + 0.05;
  }

  /* changed to num rec., make sure i'm not going in the wrong direction */
  if (*cx <0.0) {
    *cx = (*bx) - 1.0/GOLD * (*bx-*ax);
  }
  if (*cx >=1.0) {
    *cx = 1.0;
  }
  
  *fc = fonedim(*cx);
  ctf = 1;

  while (*fb >= *fc) {  // Peters code: >=
    r = (*bx-*ax) * (*fb-*fc);
    q = (*bx-*cx) * (*fb-*fa);
    u = (*bx) - ((*bx-*cx)*q - (*bx-*ax)*r) / 
          (2.0 * SIGN( MAX(FABS(q-r), TINY), q-r));
    ulim = (*bx) + cg_glimit * (*cx-*bx);
    if ((*bx-u) * (u-*cx) > 0.0) {
      fu = fonedim(u);
      ctf++;
      if (fu < *fc) {
        *ax = (*bx);
        *bx = u;
        *fa = (*fb);
        *fb = fu;
	return ctf;
      } 
      else if (fu > *fb) {
        *cx = u;
        *fc = fu;
        return ctf;
      }
      u  = (*cx) + GOLD * (*cx-*bx);
      fu = fonedim(u);
      ctf++;
    } 
    else if ((*cx-u) * (u-ulim) > 0.0) {
      fu = fonedim(u);
      ctf++;
      if (fu < *fc) {
        SHFT(*bx, *cx, u, *cx + GOLD * (*cx - *bx))
        SHFT(*fb, *fc, fu, fonedim(u))
        ctf++;
      }
    } 
    else if ((u-ulim) * (ulim-*cx) >= 0.0) {
      u  = ulim;
      fu = fonedim(u);
      ctf++;
    } 
    else {
      u  = (*cx) + GOLD * (*cx-*bx);
      fu = fonedim(u);
      ctf++;
    }
    SHFT(*ax, *bx, *cx,  u)
    SHFT(*fa, *fb, *fc, fu)
  }
  return ctf;
}

#undef GOLD
//#undef GLIMIT
#undef TINY
//#undef SHFT

#undef CGOLD
//#undef ZEPS

/* #undef ITMAX */
