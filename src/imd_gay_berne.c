#ifdef UNIAX

/********************************************************/
/*   interaction between two identical                  */
/*   uniaxial Gay-Berne particles.                      */
/********************************************************/

#include "imd.h"

void gay_berne ( vektor r12, vektor e1, vektor e2, real rsqr, 
		 vektor s1, vektor w1, 
		 real *pot12, vektor *force12, 
		 vektor *torque12, vektor *torque21 )
{
  const real huge = 1.0e+30 ;
  const real tiny = 1.0e-10 ;
  const real crit = 1.0e-03 ;

  /* Gay-Berne exponents */

  const real mu = 2.0 ;
  const real nu = 1.0 ;
  const real omu = 1.0 / mu ;
  
  /* uniaxial Gay-Berne parameters */

  /* anisotropy parameters */

  real chi  ;  
  real chip ;

  /* minimum contact distance */

  real sig0 ;

  /* energy unit */

  real eps0 ;

  /* local variables */

  real rr, orr, re1, re2, e12 ;
  real eps, epsp, well, aid, scalar, sig ;
  real aux, aux2, aux6, aux12 ;
  real rep, rem, repp, remp ;
  real pot_rr, pot_re1, pot_re2, pot_e12 ;

  /******************************************************/

  /* check for uniaxiality */

  if ( s1.x != s1.y ) {
    error("not a uniaxial molecule!\n"); }
  
  if ( w1.x != w1.y ) {
    error("not a uniaxial molecule!\n"); }

  /* Gay-Berne anisotropy parameters */

  chi = ( pow(s1.z,2.0)  - pow(s1.x,2.0) )
    / ( pow(s1.z,2.0) + pow(s1.x,2.0) ) ;

  chip = ( pow(w1.x,omu) - pow(w1.z,omu) )
    / ( pow(w1.x,omu) + pow(w1.z,omu) ) ;

  /* minimum contact distance */

  if ( s1.x <= s1.z ) 
    sig0 = s1.x ;
  else
    sig0 = s1.z ;

  /* energy unit */

  eps0 = 1.0 ;

  /******************************************************/

  /* setup initial values: a steep repulsive barrier */

  *pot12 = huge ;

  force12->x = - huge ;
  force12->y = - huge ;
  force12->z = - huge ;

  torque12->x = - huge ;
  torque12->y = - huge ;
  torque12->z = - huge ;

  torque21->x = - huge ;
  torque21->y = - huge ;
  torque21->z = - huge ;

  /******************************************************/

  /* distance between Gay-Berne centers */

  rr = sqrt( rsqr ) ;
  orr = 1.0 / ( rr + tiny ) ;

  /* arguments of Gay-Berne potential */

  re1 = orr * SPROD(r12,e1) ;
  re2 = orr * SPROD(r12,e2) ;
  e12 = SPROD(e1,e2) ;

  /******************************************************/

  /* Gay-Berne epsilon function */

  eps = eps0 / sqrt( 1.0 - pow( chi * e12, 2.0 ) + tiny ) ;

  /* Gay-Berne epsilon prime function */

  epsp = 1.0 - 0.5 * chip * 
    ( pow( re1 + re2, 2.0 ) / ( 1.0 + chip * e12 + tiny ) +
      pow( re1 - re2, 2.0 ) / ( 1.0 - chip * e12 + tiny ) ) ;
      
  /* Gay-Berne well */

  well = pow( eps, nu ) * pow( epsp, mu ) ;

  /* Gay-Berne sigma function */

  aid = 1.0 - 0.5 * chi * 
    ( pow( re1 + re2, 2.0 ) / ( 1.0 + chi * e12 + tiny ) +
      pow( re1 - re2, 2.0 ) / ( 1.0 - chi * e12 + tiny ) ) ;
  
  sig = sig0 / sqrt( aid + tiny ) ;

  /* Gay-Berne 12-6 expressions */

  scalar = rr - sig + sig0 ;

  if ( scalar > crit ) 
    {
      aux = sig0 / scalar ;
      aux2 = aux * aux ;
      aux6 = aux2 * aux2 * aux2 ;
      aux12 = aux6 * aux6 ;

      /* pair potential energy */

      *pot12 = 4.0 * well * ( aux12 - aux6 ) ;

      /***************************************************/
      /* forces and torques between pair of molecules */

      rep = ( re1 + re2 ) / ( 1.0 + chi * e12 + tiny ) ;
      rem = ( re1 - re2 ) / ( 1.0 - chi * e12 + tiny ) ;

      repp = ( re1 + re2 ) / ( 1.0 + chip * e12 + tiny ) ;
      remp = ( re1 - re2 ) / ( 1.0 - chip * e12 + tiny ) ;

      /* derivative of potential energy 
	 with respect to molecular distance r(12) */
      
      pot_rr = well * ( - 12.0 * aux12 + 6.0 * aux6 )
	* aux / sig0 ;

      /* derivative of potential energy 
	 with respect to cos(theta1) */

      pot_re1 = - well * mu * chip * ( repp + remp )
	* ( aux12 - aux6 ) / ( epsp + tiny )
	+ 3.0 * well * chi * ( rep + rem ) * aux
	* ( 2.0 * aux12 - aux6 ) 
	/ sqrt( pow( aid + tiny, 3.0 ) ) ;

      /* derivative of potential energy 
	 with respect to cos(theta2) */
      
      pot_re2 = - well * mu * chip * ( repp - remp )
	* ( aux12 - aux6 ) / ( epsp + tiny )
	+ 3.0 * well * chi * ( rep - rem ) * aux
	* ( 2.0 * aux12 - aux6 ) 
	/ sqrt( pow( aid + tiny, 3.0 ) ) ;

      /* derivative of potential energy 
	 with respect to cos(phi12) */

      pot_e12 = well * nu * pow( chi, 2.0 ) * e12
	* ( aux12 - aux6 ) 
	/ ( 1.0 - pow( chi * e12, 2.0 ) + tiny )
	- 0.5 * well * mu * pow( chip, 2.0 )
	* ( aux12 - aux6 ) 
	* ( - pow( re1 + re2, 2.0 ) /
	    pow( 1.0 + chip * e12 + tiny, 2.0 )
	    + pow( re1 - re2, 2.0 ) /
	    pow( 1.0 - chip * e12 + tiny, 2.0 ) )
	/ ( epsp + tiny )
	+ 0.5 * well * pow( chi, 2.0 ) * aux
	* ( 6.0 * aux12 - 3.0 * aux6 ) 
	* ( - pow( re1 + re2, 2.0 ) /
	    pow( 1.0 + chi * e12 + tiny, 2.0 )
	    + pow( re1 - re2, 2.0 ) /
	    pow( 1.0 - chi * e12 + tiny, 2.0 ) )
	/ sqrt( pow( aid + tiny , 3.0 ) ) ;

      /* force on molecule 1 due to molecule 2 */

      force12->x = - 4.0 * 
	( orr * pot_rr * r12.x + pot_re1 * orr *
	  ( e1.x - orr * re1 * r12.x ) + pot_re2 * orr *
	  ( e2.x - orr * re2 * r12.x ) ) ;

      force12->y = - 4.0 * 
	( orr * pot_rr * r12.y + pot_re1 * orr *
	  ( e1.y - orr * re1 * r12.y ) + pot_re2 * orr *
	  ( e2.y - orr * re2 * r12.y ) ) ;

      force12->z = - 4.0 * 
	( orr * pot_rr * r12.z + pot_re1 * orr *
	  ( e1.z - orr * re1 * r12.z ) + pot_re2 * orr *
	  ( e2.z - orr * re2 * r12.z ) ) ;
      
      /* torque on molecule 1 due to molecule 2 */
      
      torque12->x = - 4.0 *
	( orr * pot_re1 * ( e1.y * r12.z - e1.z * r12.y )
	  + pot_e12 * ( e1.y * e2.z - e1.z * e2.y ) ) ;

      torque12->y = - 4.0 *
	( orr * pot_re1 * ( e1.z * r12.x - e1.x * r12.z )
	  + pot_e12 * ( e1.z * e2.x - e1.x * e2.z ) ) ;

      torque12->z = - 4.0 *
	( orr * pot_re1 * ( e1.x * r12.y - e1.y * r12.x )
	  + pot_e12 * ( e1.x * e2.y - e1.y * e2.x ) ) ;

      /* torque on molecule 2 due to molecule 1 */
      
      torque21->x = - 4.0 * 
	( orr * pot_re2 * ( e2.y * r12.z - e2.z * r12.y )
	  + pot_e12 * ( e2.y * e1.z - e2.z * e1.y ) ) ;

      torque21->y = - 4.0 * 
	( orr * pot_re2 * ( e2.z * r12.x - e2.x * r12.z )
	  + pot_e12 * ( e2.z * e1.x - e2.x * e1.z ) ) ;

      torque21->z = - 4.0 * 
	( orr * pot_re2 * ( e2.x * r12.y - e2.y * r12.x )
	  + pot_e12 * ( e2.x * e1.y - e2.y * e1.x ) ) ;
    }
}

#endif
