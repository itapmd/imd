/*****************************************************************************
*
*  imd_forces_ttbp.c -- force loop for three body forces
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/******************************************************************************
*
*  three body forces for TTBP, using neighbor tables computed in do_forces
*
******************************************************************************/

void do_forces_ttbp(cell *p)
{
  vektor force,d,dd,tmp_d,dtheta_dr;
  int    i,j,k,s_k;					       /* counter */
  int    p_typ,j_typ,k_typ;				         /* sorte */
  neightab  *neigh;                     /* neighbor table of one particle */
  real   *tmp_ptr;                                               /* dummy */
  real   dE_dr;
  real   pot_zwi,tmp_pot;				           /* pot */
  real   radius,radius2,rradius,rradius2;		      /* distance */
  real   *s_potptr;				                   /* pot */
  real   s_chi,s_pot_k0,s_pot_k1,s_pot_k2;                         /* pot */
  real   s_dv,s_d2v;                                               /* pot */
  real   s_pot_grad_ik,s_pot_zwi_ik;		                   /* pot */
  real   s_pot_grad_ij,s_pot_zwi_ij,tmp_f2;	                   /* pot */
  real   theta_grad;                /* actual value of theta(jik) in grad */
  real   theta_rad;               /* actual value of theta(jik) in radian */
  real   cos_theta;                      	    /* cosine(theta(jik)) */
  real   dE_dtheta;	                                   /* dE / dtheta */
  real   d_acos;	                        /* d cos(theta) / d theta */
  real   tmp_factor;	                             
  real   tmp_denom;	         /* inverse of denominator radius*rradius */
  real	 tmp_sp;	                    /* SPROD of vectors ij and ik */
  real   tmp_pi  = 3.141592654/180.0;                              /* ... */
  real   tmp_sin;                                                /* dummy */
  real   tmp_virial = 0.0;
#ifdef P_AXIAL
  vektor tmp_vir_vect;
  tmp_vir_vect.x = 0.0;
  tmp_vir_vect.y = 0.0;
  tmp_vir_vect.z = 0.0;
#endif

  /* The angle theta is the angle at atom i formed by the atoms j-i-k */

  /* For each atom in cell */
  for (i=0; i<p->n; ++i) {

    p_typ    = p->sorte[i];
    neigh    = p->neigh[i];

    /* first neighbor in pair */
    for (j=0; j<neigh->n-1; ++j) {

      j_typ    = neigh->typ[j];
      tmp_ptr  = &neigh->dist[j*3];
      d.x      = *tmp_ptr; ++tmp_ptr;
      d.y      = *tmp_ptr; ++tmp_ptr;
      d.z      = *tmp_ptr;

      /* Calculate distance ij */
      radius2  = SPROD(d,d);
      /* Check for distances, shorter than minimal distance */
      if (radius2 <= ttbp_r2_0) radius2 = ttbp_r2_0;
      radius = sqrt(radius2);

      /* ttbp smooth function ij: potential table */
      s_k      = (int) ((radius2 - ttbp_r2_0) * ttbp_inv_r2_step);
      s_chi    = (radius2 - ttbp_r2_0 - s_k * ttbp_r2_step) * ttbp_inv_r2_step;
      s_potptr = PTR_3D_V(ttbp_potential, s_k, p_typ, j_typ, ttbp_pot_dim);
      s_pot_k0 = *s_potptr; s_potptr += ttbp_pot_dim.y * ttbp_pot_dim.z;
      s_pot_k1 = *s_potptr; s_potptr += ttbp_pot_dim.y * ttbp_pot_dim.z;
      s_pot_k2 = *s_potptr;
      s_dv     = s_pot_k1 - s_pot_k0;
      s_d2v    = s_pot_k2 - 2 * s_pot_k1 + s_pot_k0;
      s_pot_grad_ij = 2 * ttbp_inv_r2_step * ( s_dv + (s_chi - 0.5) * s_d2v );
      s_pot_zwi_ij = s_pot_k0 + s_chi * s_dv + 0.5 * s_chi * (s_chi-1) * s_d2v;

      /* second neighbor in pair */
      for (k=j+1; k<neigh->n; ++k) {

        k_typ    = neigh->typ[k];
        tmp_ptr  = &neigh->dist[k*3];
        dd.x     = *tmp_ptr; ++tmp_ptr;
        dd.y     = *tmp_ptr; ++tmp_ptr;
        dd.z     = *tmp_ptr;

        /* Calculate distance ik */
        rradius2  = SPROD(dd,dd);
        /* Check for distances, shorter than minimal distance */
        if (rradius2 <= ttbp_r2_0) rradius2 = ttbp_r2_0;
	rradius  = sqrt(rradius2);

        /* ttbp smooth function ik: potential table */
        s_k   = (int) ((rradius2 - ttbp_r2_0) * ttbp_inv_r2_step);
        s_chi = (rradius2 - ttbp_r2_0 - s_k * ttbp_r2_step) * ttbp_inv_r2_step;
        s_potptr = PTR_3D_V(ttbp_potential, s_k, p_typ, k_typ , ttbp_pot_dim);
        s_pot_k0 = *s_potptr; s_potptr += ttbp_pot_dim.y * ttbp_pot_dim.z;
        s_pot_k1 = *s_potptr; s_potptr += ttbp_pot_dim.y * ttbp_pot_dim.z;
        s_pot_k2 = *s_potptr;
        s_dv     = s_pot_k1 - s_pot_k0;
        s_d2v    = s_pot_k2 - 2 * s_pot_k1 + s_pot_k0;
        s_pot_grad_ik = 2 * ttbp_inv_r2_step * ( s_dv + (s_chi-0.5) * s_d2v );
        s_pot_zwi_ik = s_pot_k0 + s_chi * s_dv + 0.5*s_chi * (s_chi-1) * s_d2v;

	/* Calculate SP of vectors ij and ik */
	tmp_sp    = SPROD(d,dd);
	/* Calculate cosine(theta(jik)) */
	cos_theta = tmp_sp/radius/rradius;

        /*  incase some security is needed:
        if (cos_theta >  1.0) { cos_theta =  1.0; };
        if (cos_theta < -1.0) { cos_theta = -1.0; };
        */

	/* Calculate theta(jik) */
	theta_rad = acos(cos_theta);
	/* calculate f2 = f(ij)*f(ik) */
 	tmp_f2    = s_pot_zwi_ij * s_pot_zwi_ik;       /* smooth */

        /* FOURIER potential term */
        pot_zwi = ttbp_constant[p_typ]*(ttbp_c0[p_typ] +
                                        ttbp_c1[p_typ] * cos_theta +
                                        ttbp_c2[p_typ] * cos(2*theta_rad) +
                                        ttbp_c3[p_typ] * cos(3*theta_rad));
	tmp_pot		= pot_zwi ;
 	pot_zwi         = pot_zwi * tmp_f2 ;	       /* smooth */
	p->pot_eng[i]  += pot_zwi;
	tot_pot_energy += pot_zwi;

	/* Forces ... the horror starts 
	 * F = - d E / d r  
	 *   = - d E / d Theta  *  d Theta / d r
	 *   = - d E / d Theta  *  d_acos *  dtheta_dr 
	 * smoothing 
	 * F = - d (E*Fij*Fik) / d r
	 *   = - d E / d r * Fij * Fik   - E * d Fij 
         *                    / d r * Fik  - E * Fij * d Fik / d r
	 */

        dE_dtheta = - ttbp_constant[p_typ] *
                      (  ttbp_c1[p_typ] * sin(theta_rad) + 
                       2*ttbp_c2[p_typ] * sin(2*theta_rad) + 
                       3*ttbp_c3[p_typ] * sin(3*theta_rad));

        tmp_sin = sin(theta_rad);
        /*  in case some security is needed
        if ( (tmp_sin <  0.0001) && (tmp_sin >= 0.0) ) tmp_sin =  0.0001;
        if ( (tmp_sin > -0.0001) && (tmp_sin <= 0.0) ) tmp_sin = -0.0001;
        */
        d_acos          = - 1.0 / tmp_sin;
	tmp_denom	=   1.0 / (radius*rradius);
	tmp_factor	=   tmp_denom * d_acos;

	dtheta_dr.x	= tmp_factor * ( d.x  * ( tmp_sp / radius2  - 1) +
			                 dd.x * ( tmp_sp / rradius2 - 1) ) ;
	dtheta_dr.y	= tmp_factor * ( d.y  * ( tmp_sp / radius2  - 1) +
			                 dd.y * ( tmp_sp / rradius2 - 1) ) ;
	dtheta_dr.z	= tmp_factor * ( d.z  * ( tmp_sp / radius2  - 1) +
			                 dd.z * ( tmp_sp / rradius2 - 1) ) ;

	force.x  = - dE_dtheta * dtheta_dr.x * tmp_f2
	           + tmp_pot   * (  s_pot_grad_ij * d.x   * s_pot_zwi_ik
	                          + s_pot_grad_ik * dd.x  * s_pot_zwi_ij);
	force.y  = - dE_dtheta * dtheta_dr.y * tmp_f2
	           + tmp_pot   * (  s_pot_grad_ij * d.y   * s_pot_zwi_ik
	                          + s_pot_grad_ik * dd.y  * s_pot_zwi_ij);
	force.z  = - dE_dtheta * dtheta_dr.z * tmp_f2
	           + tmp_pot   * (  s_pot_grad_ij * d.z   * s_pot_zwi_ik
	                          + s_pot_grad_ik * dd.z  * s_pot_zwi_ij);
	p->kraft X(i)  += force.x;
	p->kraft Y(i)  += force.y;
	p->kraft Z(i)  += force.z;

        tmp_d.x         = p->ort X(i);
    	tmp_d.y         = p->ort Y(i);
    	tmp_d.z         = p->ort Z(i);
#ifdef P_AXIAL
        tmp_vir_vect.x -= tmp_d.x * force.x;
        tmp_vir_vect.y -= tmp_d.y * force.y;
        tmp_vir_vect.z -= tmp_d.z * force.z;
#else
        tmp_virial     -= SPROD(tmp_d,force);
#endif
      }
    }
  }

#ifdef P_AXIAL
  vir_vect.x += tmp_vir_vect.x;
  vir_vect.y += tmp_vir_vect.y;
  virial     += tmp_vir_vect.x;
  virial     += tmp_vir_vect.y;
  vir_vect.z += tmp_vir_vect.z;
  virial     += tmp_vir_vect.z;
#else
  virial     += tmp_virial;
#endif

}


