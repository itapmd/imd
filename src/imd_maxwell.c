/******************************************************************************
*
* maxwell.c -- initialize velocity with a maxwell distribution
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/*
*
* Converted directely from an example in the book of Allen and Tildesley
*
*/

void maxwell(real temp)
               

/* *******************************************************************
   ** TRANSLATIONAL VELOCITIES FROM MAXWELL-BOLTZMANN DISTRIBUTION  **
   **                                                               **
   ** THE DISTRIBUTION IS DETERMINED BY TEMPERATURE AND (UNIT) MASS.**
   ** THIS ROUTINE IS GENERAL, AND CAN BE USED FOR ATOMS, LINEAR    **
   ** MOLECULES, AND NON-LINEAR MOLECULES.                          **
   **                                                               **
   ** ROUTINE REFERENCED:                                           **
   **                                                               **
   ** REAL FUNCTION GAUSS ( DUMMY )                                 **
   **    RETURNS A UNIFORM RANDOM NORMAL VARIATE FROM A             **
   **    DISTRIBUTION WITH ZERO MEAN AND UNIT VARIANCE.             **
   ******************************************************************* */
 
 
{ 
   int         k, natom;
   vektor      tot_impuls;
   static long dummy = 0;

#ifdef UNIAX
   real xisq;
   real xi0;
   real xi1;
   real xi2;
   real dot;
   real norm;
   real osq;
#endif

   real        TEMP;
   real scale, xx;
   int num, nhalf;

   TEMP = temp;
   tot_impuls.x = 0.0;
   tot_impuls.y = 0.0;
#ifndef TWOD
   tot_impuls.z = 0.0;
#endif

   natom = 0;

#ifdef TRANSPORT
   nhalf = tran_nlayers / 2;
   scale = tran_nlayers / box_x.x;
#endif

   /* Temperatur setzen */
   for (k=0; k<ncells; ++k) {

      int i;
      cell *p;
      p = cell_array + CELLS(k);

      for (i=0; i<p->n; ++i) {

#ifdef TRANSPORT
         /* which layer? */
         num = scale * p->ort X(i);
         if (num < 0)             num = 0;
         if (num >= tran_nlayers) num = tran_nlayers-1;

         if (num == 0) {
            TEMP = tran_Tleft;
         } else if (num == nhalf) {
            TEMP = tran_Tright;
         } else if (num < nhalf) {
            xx = p->ort X(i) - box_x.x / tran_nlayers;
            TEMP = tran_Tleft + (tran_Tright - tran_Tleft) * 
                               xx * tran_nlayers / (box_x.x * (nhalf-1));
         } else {
            xx = box_x.x - p->ort X(i) - box_x.x / tran_nlayers;
            TEMP = tran_Tleft + (tran_Tright - tran_Tleft) * 
                               xx * tran_nlayers / (box_x.x * (nhalf-1));
	 }
#endif

         if (NUMMER(p,i) >= 0) {
            p->impuls X(i) = gasdev( &seed ) * sqrt(TEMP * MASSE(p,i));
            p->impuls Y(i) = gasdev( &seed ) * sqrt(TEMP * MASSE(p,i));
#ifndef TWOD
            p->impuls Z(i) = gasdev( &seed ) * sqrt(TEMP * MASSE(p,i));
#endif
            ++natom;
            tot_impuls.x += p->impuls X(i);
            tot_impuls.y += p->impuls Y(i);
#ifndef TWOD
            tot_impuls.z += p->impuls Z(i);
#endif
         } else {
            p->impuls X(i) = 0.0;
            p->impuls Y(i) = 0.0;
#ifndef TWOD
            p->impuls Z(i) = 0.0;
#endif
         }

#ifdef UNIAX

	   /* angular velocities for uniaxial molecules */

	   /* choose a random vector in space */

	   xisq = 1.0;

	   while ( xisq >= 1.0 ) 
	     {
	       xi1 = 2.0 * ran1( &dummy ) - 1.0 ;
	       xi2 = 2.0 * ran1( &dummy ) - 1.0 ;
	       xisq = xi1 * xi1 + xi2 * xi2 ;
	     }

	   xi0 = sqrt( 1.0 - xisq ) ;

	   p->dreh_impuls X(i) = 2.0 * xi1 * xi0 ;
	   p->dreh_impuls Y(i) = 2.0 * xi2 * xi0 ;
	   p->dreh_impuls Z(i) = 1.0 - 2.0 * xisq ;

	   /* constrain the vector to be perpendicular 
	      to the molecule */

	   dot = SPRODN(p->dreh_impuls,i,p->achse,i);

	   p->dreh_impuls X(i) -= dot * p->achse X(i) ; 
	   p->dreh_impuls Y(i) -= dot * p->achse Y(i) ; 
	   p->dreh_impuls Z(i) -= dot * p->achse Z(i) ; 

	   /* renormalize vector */	   

	   osq = SPRODN(p->dreh_impuls,i,p->dreh_impuls,i);
	   norm = sqrt( osq );

	   p->dreh_impuls X(i) /= norm ;
	   p->dreh_impuls Y(i) /= norm ;
	   p->dreh_impuls Z(i) /= norm ;

	   /* choose the magnitude of the angular momentum */

	   osq = - 2.0 * p->traeg_moment[i] * TEMP
	     * log( ran1( &dummy ) ) ;
	   norm = sqrt( osq );

	   p->dreh_impuls X(i) *= norm ;
	   p->dreh_impuls Y(i) *= norm ;
	   p->dreh_impuls Z(i) *= norm ;

           if (NUMMER(p,i) >= 0) {
              p->dreh_impuls X(i) = 0.0;
              p->dreh_impuls Y(i) = 0.0;
              p->dreh_impuls Z(i) = 0.0;
           }

#endif /* UNIAX */

#ifdef SHOCK
	 if (shock_mode == 1) {
	   if ( p->ort X(i) < shock_strip ) 
	       p->impuls X(i) += shock_speed * MASSE(p,i);
	 }
	 if (shock_mode == 2) {
	   if ( p->ort X(i) < box_x.x*0.5 ) 
	     p->impuls X(i) += shock_speed * MASSE(p,i);
	   else
	     p->impuls X(i) -= shock_speed * MASSE(p,i);
	 }
#endif
      }
   }

   /* CPU could be empty */
   if (0==natom) return;

   tot_impuls.x /= natom;
   tot_impuls.y /= natom;
#ifndef TWOD
   tot_impuls.z /= natom;
#endif

   /* Temperatur setzen */
   for (k=0; k<ncells; ++k) {
      int i;
      cell *p;
      p = cell_array + CELLS(k);
      for (i=0; i<p->n; ++i) { 
         if (NUMMER(p,i) >= 0) {
            p->impuls X(i) -= tot_impuls.x;
            p->impuls Y(i) -= tot_impuls.y;
#ifndef TWOD
            p->impuls Z(i) -= tot_impuls.z;
#endif
         }
      }
   }
} 
 
 
/* gasdev and ran1 are routines from numerical recipes to generate a
   gaussian random distribution

*/
 

float gasdev(long *idum)
{
	float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software X!05.W4z4'>4. */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software X!05.W4z4'>4. */

