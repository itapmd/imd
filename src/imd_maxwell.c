/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
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
   ivektor     nactive_vec;
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
   tot_impuls.x = 0.0;   nactive_vec.x = 0;
   tot_impuls.y = 0.0;   nactive_vec.y = 0;
#ifndef TWOD
   tot_impuls.z = 0.0;   nactive_vec.z = 0;
#endif

#ifdef NVX
   nhalf = tran_nlayers / 2;
   scale = tran_nlayers / box_x.x;
#endif

   /* set temperature */
   for (k=0; k<ncells; ++k) {

      int i;
      cell *p;
      p = cell_array + CELLS(k);

      for (i=0; i<p->n; ++i) {

#ifdef NVX
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

#ifdef FTG
	  /* calc slice and set TEMP  */
	 tmp = p->ort X(i)/box_x.x;
	 slice = (int) nslices * tmp;
	 if (slice<0)        slice = 0;
	 if (slice>=nslices) slice = nslices -1;;
	 
	 TEMP=  Tleft + (Tright-Tleft)*(slice-nslices_Left+1) /
	   (real) (nslices-nslices_Left-nslices_Right+1);
    
	 if(slice>=nslices-nslices_Right)  TEMP = Tright;
	 if(slice<nslices_Left)            TEMP=  Tleft;
#endif
         
	 tmp = sqrt(TEMP * MASSE(p,i));
         typ = VSORTE(p,i);
         p->impuls X(i) = gaussian(tmp) * (restrictions + typ)->x;
         p->impuls Y(i) = gaussian(tmp) * (restrictions + typ)->y;
#ifndef TWOD
         p->impuls Z(i) = gaussian(tmp) * (restrictions + typ)->z;
#endif
         nactive_vec.x += (int) (restrictions + typ)->x;
         nactive_vec.y += (int) (restrictions + typ)->y;
#ifndef TWOD
         nactive_vec.z += (int) (restrictions + typ)->z;
#endif
         tot_impuls.x += p->impuls X(i);
         tot_impuls.y += p->impuls Y(i);
#ifndef TWOD
         tot_impuls.z += p->impuls Z(i);
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

         p->dreh_impuls X(i) = 2.0 * xi1 * xi0 ;
         p->dreh_impuls Y(i) = 2.0 * xi2 * xi0 ;
         p->dreh_impuls Z(i) = 1.0 - 2.0 * xisq ;

        /* constrain the vector to be perpendicular to the molecule */

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

        osq = - 2.0 * p->traeg_moment[i] * TEMP * log( drand48() ) ;
        norm = sqrt( osq );

        p->dreh_impuls X(i) *= norm ;
        p->dreh_impuls Y(i) *= norm ;
        p->dreh_impuls Z(i) *= norm ;

#endif /* UNIAX */

#ifdef SHOCK
	/* plate against bulk */
	 if (shock_mode == 1) {
	   if ( p->ort X(i) < shock_strip ) 
	       p->impuls X(i) += shock_speed * MASSE(p,i);
	 }
	 /* two halves against one another */
	 if (shock_mode == 2) {
	   if ( p->ort X(i) < box_x.x*0.5 ) 
	     p->impuls X(i) += shock_speed * MASSE(p,i);
	   else
	     p->impuls X(i) -= shock_speed * MASSE(p,i);
	 }
	 /* bulk against wall */
	 if (shock_mode == 3) p->impuls X(i) += shock_speed * MASSE(p,i);
#endif
      }
   }

   tot_impuls.x = nactive_vec.x == 0 ? 0.0 : tot_impuls.x / nactive_vec.x;
   tot_impuls.y = nactive_vec.y == 0 ? 0.0 : tot_impuls.y / nactive_vec.y;
#ifndef TWOD
   tot_impuls.z = nactive_vec.z == 0 ? 0.0 : tot_impuls.z / nactive_vec.z;
#endif

   /* correct center of mass momentum */
   for (k=0; k<ncells; ++k) {
      int i;
      cell *p;
      p = cell_array + CELLS(k);
      for (i=0; i<p->n; ++i) {
         typ = VSORTE(p,i);
         p->impuls X(i) -= tot_impuls.x * (restrictions + typ)->x;
         p->impuls Y(i) -= tot_impuls.y * (restrictions + typ)->y;
#ifndef TWOD
         p->impuls Z(i) -= tot_impuls.z * (restrictions + typ)->z;
#endif
      }
   }
} 

/* The following is adapted from the GNU Scientific Library at
   http://www.gnu.org/software/gsl/ 
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
