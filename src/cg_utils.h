#include <math.h>
#define NRANSI
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


#define GOLD 1.618034
#define GLIMIT 100.0/* was  100.0 */
#define TINY 1.0e-20

/* brent from Num. Rec. modified for our cg */
/* fb should be the minimal value */
int brent(real ax, real bx, real cx, real fb,real *alphamin)
{
	int iter;
	real a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	real e=0.0;

	real tol = linmin_tol;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;

/*	fw=fv=fx=(*f)(x); */
	fw=fv=fx=fb;

	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*alphamin=x;
			/* return fx; */
			/* the ort should be in the minimum position, with the correct Epot & fnorm */
			return (iter); 
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));

		fu=fonedim(u);

		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	error("Too many iterations in brent");

}



/* braketing of minimum from Num Rec. */
 
int mnbrak(real *ax, real *bx, real *cx, real *fa, real *fb, real *fc)
{
	real ulim,u,r,q,fu,dum;
	int ctf;
/*	*fa=(*func)(*ax);   = old_cgval */
/*	*fb=(*func)(*bx);   is already initialised */ 
	
	/*  printf("start mnbrak\n");fflush(stdout); */

	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);

	/* *fc=(*func)(*cx);*/
	*fc = fonedim(*cx);
	ctf=1;
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
		    /* fu=(*func)(u); */
		    fu=fonedim(u);
		    ctf++;
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return(ctf);
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return(ctf);
			}
			u=(*cx)+GOLD*(*cx-*bx);
			/*fu=(*func)(u);*/
			fu=fonedim(u);
			ctf++;
		} else if ((*cx-u)*(u-ulim) > 0.0) {
		    /* fu=(*func)(u); */
		    fu=fonedim(u);
		    ctf++;
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				  /* SHFT(*fb,*fc,fu,(*func)(u))*/
				SHFT(*fb,*fc,fu,fonedim(u))
				  ctf++;
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			/*	fu=(*func)(u);*/
			fu=fonedim(u);
			ctf++;
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			/* fu=(*func)(u); */
			fu=fonedim(u);
			ctf++;
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	     }
	/*   printf("in mnbrack: ctf = %d  \n",ctf);fflush(stdout); */
	return(ctf);
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software X!05.W4z4'>4. */


#undef ITMAX
#undef CGOLD
#undef ZEPS
/* (C) Copr. 1986-92 Numerical Recipes Software X!05.W4z4'>4. */







