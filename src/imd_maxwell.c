/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2004 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* maxwell.c -- initialize velocity with a maxwell distribution
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

real gaussian(const real);

/*
*
* Converted directely from an example in the book of Allen and Tildesley
*
*/

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

void maxwell(real temp)
{ 
   int         k;
   vektor      tot_impuls;
   int         nactive_x, nactive_y, nactive_z;
   static long dummy = 0;
   int slice;

#ifdef UNIAX
   real xisq;
   real xi0;
   real xi1;
   real xi2;
   real dot;
   real norm;
   real osq;
#endif
   real   TEMP;
   real   scale, xx, tmp;
   int    num, nhalf, typ;

   TEMP = temp;
   tot_impuls.x = 0.0;   nactive_x = 0;
   tot_impuls.y = 0.0;   nactive_y = 0;
#ifndef TWOD
   tot_impuls.z = 0.0;   nactive_z = 0;
#endif

#ifdef NVX
   nhalf = tran_nlayers / 2;
   scale = tran_nlayers / box_x.x;
#endif


   /* set temperature */
   for (k=0; k<NCELLS; ++k) {

      int i;
      cell *p;
      vektor *rest;

      p = CELLPTR(k);

      for (i=0; i<p->n; ++i) {

#ifdef NVX
         /* which layer? */
         num = scale * ORT(p,i,X);
         if (num < 0)             num = 0;
         if (num >= tran_nlayers) num = tran_nlayers-1;

         if (num == 0) {
            TEMP = tran_Tleft;
         } else if (num == nhalf) {
            TEMP = tran_Tright;
         } else if (num < nhalf) {
            xx = ORT(p,i,X) - box_x.x / tran_nlayers;
            TEMP = tran_Tleft + (tran_Tright - tran_Tleft) * 
                               xx * tran_nlayers / (box_x.x * (nhalf-1));
         } else {
            xx = box_x.x - ORT(p,i,X) - box_x.x / tran_nlayers;
            TEMP = tran_Tleft + (tran_Tright - tran_Tleft) * 
                               xx * tran_nlayers / (box_x.x * (nhalf-1));
	 }
#endif

#ifdef FTG
	  /* calc slice and set TEMP  */
	 tmp = ORT(p,i,X) / box_x.x;
	 slice = (int) nslices * tmp;
	 if (slice<0)        slice = 0;
	 if (slice>=nslices) slice = nslices -1;;
	 
	 TEMP=  Tleft + (Tright-Tleft)*(slice-nslices_Left+1) /
	   (real) (nslices-nslices_Left-nslices_Right+1);
    
	 if (slice>=nslices-nslices_Right) TEMP = Tright;
	 if (slice<nslices_Left)           TEMP=  Tleft;
#endif
         
	 tmp  = sqrt(TEMP * MASSE(p,i));
         rest = restrictions + VSORTE(p,i);
#ifndef RIGID
         IMPULS(p,i,X) = gaussian(tmp) * rest->x;
         IMPULS(p,i,Y) = gaussian(tmp) * rest->y;
#ifndef TWOD
         IMPULS(p,i,Z) = gaussian(tmp) * rest->z;
#endif
#else
	 /* superatoms get velocity zero */
	 if (superatom[VSORTE(p,i)]>-1) {
	   IMPULS(p,i,X) = 0.0;
	   IMPULS(p,i,Y) = 0.0;
#ifndef TWOD
	   IMPULS(p,i,Z) = 0.0;
#endif
	 }
#endif
         nactive_x += (int) rest->x;
         nactive_y += (int) rest->y;
#ifndef TWOD
         nactive_z += (int) rest->z;
#endif
         tot_impuls.x += IMPULS(p,i,X);
         tot_impuls.y += IMPULS(p,i,Y);
#ifndef TWOD
         tot_impuls.z += IMPULS(p,i,Z);
#endif

#ifdef UNIAX

         /* angular velocities for uniaxial molecules */

         /* choose a random vector in space */

         do {
           xi1 = 2.0 * drand48() - 1.0 ;
           xi2 = 2.0 * drand48() - 1.0 ;
           xisq = xi1 * xi1 + xi2 * xi2 ;
         } while ( xisq >= 1.0 );

         xi0 = sqrt( 1.0 - xisq ) ;

         DREH_IMPULS(p,i,X) = 2.0 * xi1 * xi0 ;
         DREH_IMPULS(p,i,Y) = 2.0 * xi2 * xi0 ;
         DREH_IMPULS(p,i,Z) = 1.0 - 2.0 * xisq ;

        /* constrain the vector to be perpendicular to the molecule */

        dot = SPRODN( &DREH_IMPULS(p,i,X), &ACHSE(p,i,X) );

        DREH_IMPULS(p,i,X) -= dot * ACHSE(p,i,X) ; 
        DREH_IMPULS(p,i,Y) -= dot * ACHSE(p,i,Y) ; 
        DREH_IMPULS(p,i,Z) -= dot * ACHSE(p,i,Z) ; 

        /* renormalize vector */	   

        osq = SPRODN( &DREH_IMPULS(p,i,X), &DREH_IMPULS(p,i,X) );
        norm = sqrt( osq );

        DREH_IMPULS(p,i,X) /= norm ;
        DREH_IMPULS(p,i,Y) /= norm ;
        DREH_IMPULS(p,i,Z) /= norm ;

        /* choose the magnitude of the angular momentum */

        osq = - 2.0 * uniax_inert * TEMP * log( drand48() ) ;
        norm = sqrt( osq );

        DREH_IMPULS(p,i,X) *= norm ;
        DREH_IMPULS(p,i,Y) *= norm ;
        DREH_IMPULS(p,i,Z) *= norm ;

#endif /* UNIAX */

#ifdef SHOCK
	/* plate against bulk */
	 if (shock_mode == 1) {
	   if ( ORT(p,i,X) < shock_strip ) 
	       IMPULS(p,i,X) += shock_speed * MASSE(p,i);
	 }
	 /* two halves against one another */
	 if (shock_mode == 2) {
	   if ( ORT(p,i,X) < box_x.x*0.5 ) 
	     IMPULS(p,i,X) += shock_speed * MASSE(p,i);
	   else
	     IMPULS(p,i,X) -= shock_speed * MASSE(p,i);
	 }
	 /* bulk against wall */
	 if (shock_mode == 3) IMPULS(p,i,X) += shock_speed * MASSE(p,i);
#endif
      }
   }

#ifdef CLONE

   /* compute the total momentum afresh */
   tot_impuls.x = 0.0;
   tot_impuls.y = 0.0;
#ifndef TWOD
   tot_impuls.z = 0.0;
#endif

   /* set velocities of clones equal */
   for (k=0; k<NCELLS; k++) {

      int i, j;
      cell *p;
      p = CELLPTR(k);

      for (i=0; i<p->n; i+=nclones) {

        tot_impuls.x += nclones * IMPULS(p,i,X);
        tot_impuls.y += nclones * IMPULS(p,i,Y);
#ifndef TWOD
        tot_impuls.z += nclones * IMPULS(p,i,Z);
#endif
	for (j=1; j<nclones; j++) {
          IMPULS(p,i+j,X) = IMPULS(p,i,X);
          IMPULS(p,i+j,Y) = IMPULS(p,i,Y);
#ifndef TWOD
          IMPULS(p,i+j,Z) = IMPULS(p,i,Z);
#endif
        }
      }
   }

#endif /* CLONE */

   tot_impuls.x = nactive_x == 0 ? 0.0 : tot_impuls.x / nactive_x;
   tot_impuls.y = nactive_y == 0 ? 0.0 : tot_impuls.y / nactive_y;
#ifndef TWOD
   tot_impuls.z = nactive_z == 0 ? 0.0 : tot_impuls.z / nactive_z;
#endif

   /* correct center of mass momentum */
   for (k=0; k<NCELLS; ++k) {
      int i;
      cell *p;
      vektor *rest;
      p = CELLPTR(k);
      for (i=0; i<p->n; ++i) {
         rest = restrictions + VSORTE(p,i);
         IMPULS(p,i,X) -= tot_impuls.x * rest->x;
         IMPULS(p,i,Y) -= tot_impuls.y * rest->y;
#ifndef TWOD
         IMPULS(p,i,Z) -= tot_impuls.z * rest->z;
#endif
      }
   }

} 

/* The following is adapted from the GNU Scientific Library, Version 1.1.
   See http://www.gnu.org/software/gsl/ 
*/

/* Polar (Box-Mueller) method; See Knuth v2, 3rd ed, p122 */

real gaussian(const real sigma)
{
  double x, y, r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      x = -1 + 2 * drand48();
      y = -1 + 2 * drand48();

      /* see if it is in the unit circle */
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  /* Box-Muller transform */
  return (real) (sigma * y * sqrt (-2.0 * log (r2) / r2));
}
