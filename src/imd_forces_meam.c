
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2009 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/*****************************************************************************
*
*  imd_forces_meam.c -- force loops for MEAM potential
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"
#include "potaccess.h"

void do_forces2(cell *p, real *Epot, real *Virial, 
                real *Vir_xx, real *Vir_yy, real *Vir_zz,
                real *Vir_yz, real *Vir_zx, real *Vir_xy)
{
  static vektor *d    = NULL, *dfc  = NULL, *ds   = NULL;
  static vektor *dfl1 = NULL, *dfl2 = NULL, *dfl3 = NULL; 
  static real *r = NULL, *invr = NULL, *r2 = NULL, *invr2 = NULL;
  static real *cos = NULL, *fc = NULL, *s = NULL;
  static real *rho_a0 = NULL, *rho_a1 = NULL, *rho_a2 = NULL, *rho_a3 = NULL;
  static real *fl1  = NULL, *fl2  = NULL, *fl3  = NULL;
  static int  curr_len = 0;
  cell     *jcell, *kcell;
  int      i, j, k, l, m, jnum, knum, p_typ, j_typ, k_typ;
  neightab *neigh;
  vektor   d_jk;
  real     r2_jk;
  real     f, df, phi, dphi, pot_zwi;
  int      col1, col2, is_short=0, idummy=0, inc = ntypes * ntypes;
  real     x_ik, x_jk;
  vektor   dx_ik_k, dx_jk_k, dx_ik_j, dx_jk_j;
  real     tmp, tmp1, tmp2, tmp3, tmp4;
  vektor   tmpv;
  real     *tmpptr;
  real     c, c_red, invdelta_c, s_kij, dgc;
  real     l1, l2, l3, dl1, dl2, dl3;
  vektor   force_j, force_k;
  real     g, gamma, egamma;
  real     pref1, pref2;
  real     rho, drho, rho_0, rho2_1, rho2_2, rho2_3;
  vektor   drho_0, drho2_1, drho2_2, drho2_3;
  real     epsilon = 1e-10;
  real     meam_t1_av, meam_t2_av, meam_t3_av;
  real     t1, t2, t3;
  real     tmp_virial = 0.0;
#ifdef P_AXIAL
  vektor   tmp_vir_vect = {0.0, 0.0, 0.0};
#endif

  if (curr_len < neigh_len) {
    d       = (vektor *) realloc( d,      neigh_len * sizeof(vektor) );
    r       = (real   *) realloc( r,      neigh_len * sizeof(real) );
    invr    = (real   *) realloc( invr,   neigh_len * sizeof(real) );
    r2      = (real   *) realloc( r2,     neigh_len * sizeof(real) );
    invr2   = (real   *) realloc( invr2,  neigh_len * sizeof(real) );
    fc      = (real   *) realloc( fc,     neigh_len * sizeof(real) );
    dfc     = (vektor *) realloc( dfc,    neigh_len * sizeof(vektor) );
    rho_a0  = (real   *) realloc( rho_a0, neigh_len * sizeof(real) );
    rho_a1  = (real   *) realloc( rho_a1, neigh_len * sizeof(real) );
    rho_a2  = (real   *) realloc( rho_a2, neigh_len * sizeof(real) );
    rho_a3  = (real   *) realloc( rho_a3, neigh_len * sizeof(real) );
    fl1     = (real   *) realloc( fl1,    neigh_len * sizeof(real) );
    fl2     = (real   *) realloc( fl2,    neigh_len * sizeof(real) );
    fl3     = (real   *) realloc( fl3,    neigh_len * sizeof(real) );
    dfl1    = (vektor *) realloc( dfl1,   neigh_len * sizeof(vektor) );
    dfl2    = (vektor *) realloc( dfl2,   neigh_len * sizeof(vektor) );
    dfl3    = (vektor *) realloc( dfl3,   neigh_len * sizeof(vektor) );
    s       = (real   *) realloc( s,      neigh_len * sizeof(real) );
    ds      = (vektor *) realloc( ds,  neigh_len * neigh_len * sizeof(vektor));
    cos     = (real   *) realloc( cos, neigh_len * neigh_len * sizeof(real) );
    if ( d==NULL || r==NULL || fc==NULL || dfc==NULL || rho_a0==NULL 
        || rho_a1==NULL || rho_a2==NULL || rho_a3==NULL || fl1==NULL 
        || fl2==NULL || fl3==NULL || dfl1==NULL || dfl2==NULL || dfl3==NULL 
        || s==NULL || ds==NULL || invr==NULL || r2==NULL || invr2==NULL 
        || cos==NULL )
      error("Cannot allocate memory for neighbour data!");
    curr_len = neigh_len;
  }

  /* for each atom in cell */
  for (i=0; i<p->n; ++i) {

    p_typ = SORTE(p,i);
    neigh = NEIGH(p,i);

    /* first loop: compute distances */

    tmpptr = neigh->dist;
    /* for each neighbor of i */
    for (j=0; j<neigh->n; ++j) {

      d[j].x   = *tmpptr++;
      d[j].y   = *tmpptr++;
      d[j].z   = *tmpptr++;
      r2[j]    = SPROD(d[j],d[j]);
      invr2[j] = 1.0 / r2[j];
      r[j]     = sqrt(r2[j]);
      invr[j]  = 1.0 / r[j];
    }

    /* second loop: compute screening function, cutoff function,
       and atomic electron density functions */

    /* for each neighbor of i */
    for (j=0; j<neigh->n; ++j) {

      j_typ  = NSORTE(neigh,j);

      s[j] = 1.0;

      if ( meam_cmax[0][0][0] == meam_cmin[0][0][0] ) {
	/* switch off screening function */
	for (m=0; m<neigh->n; ++m) {
	  ds I(j,m) = nullvek;
	}
      }
      else {

	/* initializations */
	ds I(j,j) = nullvek;

	/* for each neighbor of i other than j */
	for (k=0; k<neigh->n; ++k) {

	  cos I(j,k) = SPROD(d[j],d[k]) * invr[j] * invr[k];

	  if ( k!=j ) {

	    ds I(j,k) = nullvek;

	    /* location of atom k */	    
	    d_jk.x = d[k].x - d[j].x;
	    d_jk.y = d[k].y - d[j].y;
	    d_jk.z = d[k].z - d[j].z;

	    if ( cos I(j,k) > epsilon && SPROD(d_jk,d[j]) < -epsilon ) {

	      r2_jk  = SPROD(d_jk,d_jk);
		  
	      /* compute x_ik and x_jk */
	      x_ik = r2[k] * invr2[j];
	      x_jk = r2_jk * invr2[j];
	      
	      /* compute C_kij */
	      tmp1 = x_ik - x_jk;
	      tmp2 = 1.0 - SQR(tmp1);
	      k_typ = NSORTE(neigh,k);

	      invdelta_c = 1.0 / ( meam_cmax[k_typ][p_typ][j_typ] 
				  - meam_cmin[k_typ][p_typ][j_typ] );
	      c = ( 2.0 * ( x_ik + x_jk ) - SQR(tmp1) - 1.0 ) / tmp2;
	      c_red = ( c - meam_cmin[k_typ][p_typ][j_typ] ) * invdelta_c;
	    }
	    else
	      continue;

	    if ( c_red >= 1.0 )
	      /* S_kij=1 */
	      continue;
	    else if ( c_red <= epsilon ) {
	      /* S_kij=0 => S_ij=0 */
	      s[j] = 0.0;
	      break;
	    }
	    else {
	      s_kij = SQR( 1.0 - SQR(SQR( 1.0 - c_red )));
	      s[j] *= s_kij;
	      
	      /* compute derivatives of x_ik and x_jk */
	      tmp = 2.0 * invr2[j];
	      dx_ik_k.x = tmp * d[k].x; 
	      dx_ik_k.y = tmp * d[k].y; 
	      dx_ik_k.z = tmp * d[k].z; 
		  
	      dx_jk_k.x = tmp * d_jk.x; 
	      dx_jk_k.y = tmp * d_jk.y; 
	      dx_jk_k.z = tmp * d_jk.z; 
	      
	      dx_ik_j.x = -tmp * x_ik * d[j].x;
	      dx_ik_j.y = -tmp * x_ik * d[j].y;
	      dx_ik_j.z = -tmp * x_ik * d[j].z;
	      
	      dx_jk_j.x = -tmp * ( x_jk * d[j].x + d_jk.x );
	      dx_jk_j.y = -tmp * ( x_jk * d[j].y + d_jk.y );
	      dx_jk_j.z = -tmp * ( x_jk * d[j].z + d_jk.z );
	      
	      dgc = 8.0 * SQR( 1.0 - c_red ) * ( 1.0 - c_red ) 
		* ( 1.0 - SQR(SQR( 1.0 - c_red )) ) * invdelta_c;
	      
	      /* Compute derivatives of ln(S_ij) */
	      tmp = 2.0 * dgc / ( s_kij * tmp2 );
	      tmp1 *= c - 1.0;
	      ds IX(j,k) = tmp * ( dx_ik_k.x + dx_jk_k.x 
				  + tmp1 * ( dx_ik_k.x - dx_jk_k.x ) );
	      ds IY(j,k) = tmp * ( dx_ik_k.y + dx_jk_k.y 
				  + tmp1 * ( dx_ik_k.y - dx_jk_k.y ) );
	      ds IZ(j,k) = tmp * ( dx_ik_k.z + dx_jk_k.z 
				  + tmp1 * ( dx_ik_k.z - dx_jk_k.z ) );
	      
	      ds IX(j,j) += tmp * ( dx_ik_j.x + dx_jk_j.x 
				   + tmp1 * ( dx_ik_j.x - dx_jk_j.x ) );
	      ds IY(j,j) += tmp * ( dx_ik_j.y + dx_jk_j.y 
				   + tmp1 * ( dx_ik_j.y - dx_jk_j.y ) );
	      ds IZ(j,j) += tmp * ( dx_ik_j.z + dx_jk_j.z 
				   + tmp1 * ( dx_ik_j.z - dx_jk_j.z ) );
	    } /* c_red */
	  } /* k!=j */
	} /* k*/

      } /* end: screening function */

      if ( s[j] > 0.0 ) {

	/* compute f_cut and the derivative of ln(f_cut) */
	if ( r[j] <= meam_rcut[p_typ][j_typ]-meam_deltar[p_typ][j_typ] ) {
	  fc[j]  = 1.0;
	  dfc[j] = nullvek;
	}
	else {
	  tmp     = 1.0 
	    - ( meam_rcut[p_typ][j_typ] - r[j] ) / meam_deltar[p_typ][j_typ];
	  tmp1    = 1.0 - SQR(SQR(tmp));
	  fc[j]    = SQR(tmp1);
	  
	  tmp2 = - 8.0 / ( tmp1 * meam_deltar[p_typ][j_typ] * r[j] ) 
	    * tmp * SQR(tmp);
	  dfc[j].x = tmp2 * d[j].x;
	  dfc[j].y = tmp2 * d[j].y;
	  dfc[j].z = tmp2 * d[j].z;
	}
	
	/* compute atomic electron density functions */
	tmp1 = s[j] * fc[j] * meam_f0[j_typ];
	tmp2 = r[j] * invmeam_r0[j_typ] - 1.0;

	if ( have_eldensity_file ) {
	  col2 = p_typ * ntypes + j_typ;
	  VAL_FUNC(rho_a0[j], el_density, col2, inc, r2[j], is_short);
	  rho_a0[j] *= s[j];
	}
	else       
	  rho_a0[j] = tmp1 * exp( - meam_beta0[j_typ] * tmp2 );
	
	rho_a1[j] = tmp1 * exp( - meam_beta1[j_typ] * tmp2 ); 
	rho_a2[j] = tmp1 * exp( - meam_beta2[j_typ] * tmp2 ); 
	rho_a3[j] = tmp1 * exp( - meam_beta3[j_typ] * tmp2 ); 

      } /* s[j] > 0 */
    } /* j */

    /* third loop: compute local electron densities and some temporary sums */

    rho_0  = 0.0;
    rho2_1 = 0.0;
    rho2_2 = 0.0;
    rho2_3 = 0.0;

    if ( meam_t_average ) {
      meam_t1_av = 0.0;
      meam_t2_av = 0.0;
      meam_t3_av = 0.0;
    }

    /* for each neighbor of i */
    for (j=0; j<neigh->n; ++j) {

      if ( s[j] > 0.0 ) {

	fl1[j]    = 0.0;
	fl2[j]    = 0.0;
	fl3[j]    = 0.0;
	dfl1[j] = nullvek;
	dfl2[j] = nullvek;
	dfl3[j] = nullvek;
	
	rho_0  += rho_a0[j];

	if ( meam_t_average ) {
	  j_typ = NSORTE(neigh,j);
	  meam_t1_av += meam_t1[j_typ] * rho_a0[j];
	  meam_t2_av += meam_t2[j_typ] * rho_a0[j];
	  meam_t3_av += meam_t3[j_typ] * rho_a0[j];
	}	

	/* for each neighbor of i */
	for (k=0; k<neigh->n; ++k) {

	  if ( s[k] > 0.0 ) { 

	    /* Legendre polynomials */
	    l1 = cos I(j,k); 
	    l2 = SQR(l1) - 1.0 / 3.0;
	    l3 = l1 * ( SQR(l1) - 0.6 );

	    fl1[j] += rho_a1[k] * l1;
	    fl2[j] += rho_a2[k] * l2;
	    fl3[j] += rho_a3[k] * l3;

	    /* derivatives of Legendre polynomials */
	    dl1 = 1.0;
	    dl2 = 2.0 * l1;
	    dl3 = 3.0 * SQR(l1) - 0.6;

	    tmp = l1 * invr[j];
	    
	    tmpv.x = d[k].x * invr[k] - tmp * d[j].x;
	    tmpv.y = d[k].y * invr[k] - tmp * d[j].y;
	    tmpv.z = d[k].z * invr[k] - tmp * d[j].z;

	    tmp1 = rho_a1[k] * dl1;
	    tmp2 = rho_a2[k] * dl2;
	    tmp3 = rho_a3[k] * dl3;

	    dfl1[j].x += tmp1 * tmpv.x;
	    dfl1[j].y += tmp1 * tmpv.y;
	    dfl1[j].z += tmp1 * tmpv.z;
	    dfl2[j].x += tmp2 * tmpv.x;
	    dfl2[j].y += tmp2 * tmpv.y;
	    dfl2[j].z += tmp2 * tmpv.z;
	    dfl3[j].x += tmp3 * tmpv.x;
	    dfl3[j].y += tmp3 * tmpv.y;
	    dfl3[j].z += tmp3 * tmpv.z;
	  }
	}

	rho2_1 += fl1[j] * rho_a1[j];
	rho2_2 += fl2[j] * rho_a2[j];
	rho2_3 += fl3[j] * rho_a3[j];
      }
    } /* j */

    if ( meam_t_average ) {
      t1 = meam_t1_av / rho_0;
      t2 = meam_t2_av / rho_0;
      t3 = meam_t3_av / rho_0;
    }
    else {
      t1 = meam_t1[p_typ];
      t2 = meam_t2[p_typ];
      t3 = meam_t3[p_typ];
    }

    if ( rho_0 == 0.0 ) {
      gamma  = 0.0;
      egamma = 1.0;
      g      = 1.0;
      f      = 0.0;
      df     = 0.0;
    }
    else {
      gamma = ( t1 * rho2_1 + t2 * rho2_2 + t3 * rho2_3 ) / SQR(rho_0);
      egamma = exp( - gamma );
      g = 2.0 / ( 1.0 + egamma );
      
      /* the local electron density */
      rho = rho_0 * g; 

      if ( have_embed_potfile ) {
	PAIR_INT(f, df, embed_pot, p_typ, ntypes, rho, idummy);
	df *= 0.5;
      }
      if ( have_pre_embed_pot ) {
	/* compute embedding function analytically */
	tmp1 = rho * invmeam_rho0[p_typ];
	tmp2 = log(tmp1);
	tmp3 =  meam_e[p_typ] * meam_a[p_typ];
	f  = tmp3 * tmp1 * tmp2;
	df = tmp3 * invmeam_rho0[p_typ] * ( tmp2 + 1.0 );
      }

      /* update potential energy */
      *Epot       += f; 
      POTENG(p,i) += f;
    }


    /* fourth loop: compute forces */
    
    if ( rho_0 != 0.0 ) {
      tmp = egamma / ( 1.0 + egamma );
      pref1 = - df * g * ( 1.0 - 2.0 * gamma * tmp );
      pref2 = - df * g * tmp / rho_0;
    }
    else {
      pref1 = 0.0;
      pref2 = 0.0;
    }

    /* for each neighbor of i */
    for (j=0; j<neigh->n; ++j) {

      if ( s[j] > 0.0 ) {

	j_typ = NSORTE(neigh,j);
	jcell = NZELLE(neigh,j);
	jnum  = NNUMMER(neigh,j);      

	force_j = nullvek;

	/* compute pair interaction, contribution from atom j */
	
	col2   = p_typ * ntypes + j_typ;
	PAIR_INT(phi, dphi, pair_pot, col2, inc, r2[j], is_short);
	pot_zwi = 0.5 * s[j] * phi;

	/* update potential energy */
	*Epot       += pot_zwi;
	POTENG(p,i) += pot_zwi;

	/* force on particle j */
	force_j.x += - pot_zwi * ds IX(j,j) - 0.5 * s[j] * dphi * d[j].x;
	force_j.y += - pot_zwi * ds IY(j,j) - 0.5 * s[j] * dphi * d[j].y;
	force_j.z += - pot_zwi * ds IZ(j,j) - 0.5 * s[j] * dphi * d[j].z;

	/* compute interaction from embedding part,
	   contribution from atom j */

	tmp3 = invmeam_r0[j_typ] * invr[j];
	tmp4 = 2.0 * invr[j];
	tmpv.x = ds IX(j,j) + dfc[j].x;
	tmpv.y = ds IY(j,j) + dfc[j].y;
	tmpv.z = ds IZ(j,j) + dfc[j].z;
	
	if ( have_eldensity_file ) {
	  col1 = j_typ * ntypes + p_typ;
	  DERIV_FUNC(drho, el_density, col1, inc, r2[j], is_short);	
	  drho_0.x = rho_a0[j] * ds IX(j,j) + s[j] * drho * d[j].x; 
	  drho_0.y = rho_a0[j] * ds IY(j,j) + s[j] * drho * d[j].y; 
	  drho_0.z = rho_a0[j] * ds IZ(j,j) + s[j] * drho * d[j].z; 
	}
	else {
	  tmp = meam_beta0[j_typ] * tmp3;
	  drho_0.x = rho_a0[j] * ( tmpv.x - tmp * d[j].x );
	  drho_0.y = rho_a0[j] * ( tmpv.y - tmp * d[j].y );
	  drho_0.z = rho_a0[j] * ( tmpv.z - tmp * d[j].z );
	}

	tmp = meam_beta1[j_typ] * tmp3;
	tmp1 = rho_a1[j] * fl1[j];
	tmp2 = rho_a1[j] * tmp4;
	drho2_1.x = tmp1 * ( tmpv.x - tmp * d[j].x ) + tmp2 * dfl1[j].x;
	drho2_1.y = tmp1 * ( tmpv.y - tmp * d[j].y ) + tmp2 * dfl1[j].y;
	drho2_1.z = tmp1 * ( tmpv.z - tmp * d[j].z ) + tmp2 * dfl1[j].z;

	tmp = meam_beta2[j_typ] * tmp3;
	tmp1 = rho_a2[j] * fl2[j];
	tmp2 = rho_a2[j] * tmp4;
	drho2_2.x = tmp1 * ( tmpv.x - tmp * d[j].x ) + tmp2 * dfl2[j].x;
	drho2_2.y = tmp1 * ( tmpv.y - tmp * d[j].y ) + tmp2 * dfl2[j].y;
	drho2_2.z = tmp1 * ( tmpv.z - tmp * d[j].z ) + tmp2 * dfl2[j].z;
	
	tmp = meam_beta3[j_typ] * tmp3;
	tmp1 = rho_a3[j] * fl3[j];
	tmp2 = rho_a3[j] * tmp4;
	drho2_3.x = tmp1 * ( tmpv.x - tmp * d[j].x ) + tmp2 * dfl3[j].x;
	drho2_3.y = tmp1 * ( tmpv.y - tmp * d[j].y ) + tmp2 * dfl3[j].y;
	drho2_3.z = tmp1 * ( tmpv.z - tmp * d[j].z ) + tmp2 * dfl3[j].z;

	force_j.x += pref1 * drho_0.x 
	  + pref2 * ( t1 * drho2_1.x + t2 * drho2_2.x + t3 * drho2_3.x );
	force_j.y += pref1 * drho_0.y 
	  + pref2 * ( t1 * drho2_1.y + t2 * drho2_2.y + t3 * drho2_3.y );
	force_j.z += pref1 * drho_0.z 
	  + pref2 * ( t1 * drho2_1.z + t2 * drho2_2.z + t3 * drho2_3.z );

	KRAFT(jcell,jnum,X) += force_j.x;
	KRAFT(jcell,jnum,Y) += force_j.y;
	KRAFT(jcell,jnum,Z) += force_j.z;

	KRAFT(p,i,X) -= force_j.x;      
	KRAFT(p,i,Y) -= force_j.y;      
	KRAFT(p,i,Z) -= force_j.z;      

#ifdef P_AXIAL
        tmp_vir_vect.x += d[j].x * force_j.x;
        tmp_vir_vect.y += d[j].y * force_j.y;
        tmp_vir_vect.z += d[j].z * force_j.z;
#else
        tmp_virial     += SPROD(d[j],force_j);
#endif

#ifdef STRESS_TENS
        if (do_press_calc) {
          PRESSTENS(p,i,xx) += d[j].x * force_j.x;
          PRESSTENS(p,i,yy) += d[j].y * force_j.y;
          PRESSTENS(p,i,zz) += d[j].z * force_j.z;
          PRESSTENS(p,i,yz) += d[j].y * force_j.z;
          PRESSTENS(p,i,zx) += d[j].z * force_j.x;
          PRESSTENS(p,i,xy) += d[j].x * force_j.y;
	}
#endif

	if ( s[j] < 1.0 ) { 

	  /* for each neighbor of i other than j */
	  for (k=0; k<neigh->n; ++k) 
	    if ( k!=j ) {
	      
	      kcell = NZELLE(neigh,k);
	      knum  = NNUMMER(neigh,k);	

	      force_k = nullvek;
		
	      /* compute pair interaction, contribution from atom k */
	      force_k.x -= pot_zwi * ds IX(j,k);
	      force_k.y -= pot_zwi * ds IY(j,k);
	      force_k.z -= pot_zwi * ds IZ(j,k);

	      /* compute interaction from embedding part,
		 contribution from atom k */
	      drho_0.x = rho_a0[j] * ds IX(j,k);
	      drho_0.y = rho_a0[j] * ds IY(j,k);  
	      drho_0.z = rho_a0[j] * ds IZ(j,k);  
	      
	      if ( s[k] > 0.0 ) {
		tmp = rho_a1[j] * fl1[k];
		drho2_1.x = tmp * ds IX(j,k);
		drho2_1.y = tmp * ds IY(j,k);
		drho2_1.z = tmp * ds IZ(j,k);
		
		tmp = rho_a2[j] * fl2[k];
		drho2_2.x = tmp * ds IX(j,k); 
		drho2_2.y = tmp * ds IY(j,k); 
		drho2_2.z = tmp * ds IZ(j,k); 
		
		tmp = rho_a3[j] * fl3[k];
		drho2_3.x = tmp * ds IX(j,k); 
		drho2_3.y = tmp * ds IY(j,k); 
		drho2_3.z = tmp * ds IZ(j,k); 
	      }
	      else {
		drho2_1 = nullvek;
		drho2_2 = nullvek;
		drho2_3 = nullvek;
	      }

	      force_k.x += pref1 * drho_0.x 
		+ pref2 * ( t1 * drho2_1.x + t2 * drho2_2.x + t3 * drho2_3.x );
	      force_k.y += pref1 * drho_0.y 
		+ pref2 * ( t1 * drho2_1.y + t2 * drho2_2.y + t3 * drho2_3.y );
	      force_k.z += pref1 * drho_0.z 
		+ pref2 * ( t1 * drho2_1.z + t2 * drho2_2.z + t3 * drho2_3.z );

	      KRAFT(kcell,knum,X) += force_k.x;
	      KRAFT(kcell,knum,Y) += force_k.y;
	      KRAFT(kcell,knum,Z) += force_k.z;
	      
	      KRAFT(p,i,X) -= force_k.x;
	      KRAFT(p,i,Y) -= force_k.y;
	      KRAFT(p,i,Z) -= force_k.z;

#ifdef P_AXIAL
        tmp_vir_vect.x += d[k].x * force_k.x;
        tmp_vir_vect.y += d[k].y * force_k.y;
        tmp_vir_vect.z += d[k].z * force_k.z;

#else
        tmp_virial     += SPROD(d[k],force_k);
#endif

#ifdef STRESS_TENS
        if (do_press_calc) {
          PRESSTENS(p,i,xx) += d[k].x * force_k.x;
          PRESSTENS(p,i,yy) += d[k].y * force_k.y;
          PRESSTENS(p,i,zz) += d[k].z * force_k.z;
          PRESSTENS(p,i,yz) += d[k].y * force_k.z;
          PRESSTENS(p,i,zx) += d[k].z * force_k.x;
          PRESSTENS(p,i,xy) += d[k].x * force_k.y;
	}
#endif
	    } /* k */
	} /* s[j]<1 */
      }
    } /* j */
  } /* i */  

  /* print warning if short distance occurred */
  if (is_short==1) fprintf(stderr, "\n Short distance!\n");

#ifdef P_AXIAL
  *Vir_xx += tmp_vir_vect.x;
  *Vir_yy += tmp_vir_vect.y;
  *Virial += tmp_vir_vect.x;
  *Virial += tmp_vir_vect.y;
#ifndef TWOD
  *Vir_zz += tmp_vir_vect.z;
  *Virial += tmp_vir_vect.z;
#endif
#else
  *Virial += tmp_virial;
#endif 

}


void init_meam(void) {

  int  i, j, k, n;
  int have_common_cmax = 0;
  int have_common_cmin = 0;

  if ( ntypes == 1 || meam_cmax_lin[1] == 0.0 )
      have_common_cmax = 1;
  if ( ntypes == 1 || meam_cmin_lin[1] == 0.0 )
      have_common_cmin = 1;

  n = 0;
  for (i=0; i<ntypes; i++) {
      
      if (meam_r0[i]==0.0)
	  error("meam_r0 is zero!");
      else
	  invmeam_r0[i] = 1.0 / meam_r0[i];
      
      if (meam_rho0[i]==0.0)
	  error("meam_rho0 is zero!");
      else
	  invmeam_rho0[i] = 1.0 / meam_rho0[i];
      
      for (j=i; j<ntypes; j++) {
	  meam_deltar[i][j]  = meam_deltar[j][i]  = meam_deltar_lin[n];
	  meam_rcut[i][j]    = meam_rcut[j][i]    = meam_rcut_lin[n];
	  meam_r2_cut[i][j]  = meam_r2_cut[j][i]  = SQR(meam_rcut[i][j]); 
	  n++;
      }
  }
  
  n = 0;
  for (k=0; k<ntypes; k++)
      for (i=0; i<ntypes; i++)
	  for (j=i; j<ntypes; j++) {
	      if ( have_common_cmax ) 
		  meam_cmax[k][i][j] = meam_cmax[k][j][i] = meam_cmax_lin[0];
	      else
		  meam_cmax[k][i][j] = meam_cmax[k][j][i] = meam_cmax_lin[n];
	      
	      if ( have_common_cmin )
		  meam_cmin[k][i][j] = meam_cmin[k][j][i] = meam_cmin_lin[0];
	      else
		  meam_cmin[k][i][j] = meam_cmin[k][j][i] = meam_cmin_lin[n];
	      
	      n++;
	  }
  
  /* update neighbor table cutoff */
  if (NULL==neightab_r2cut) {
      neightab_r2cut = (real *) calloc( ntypes * ntypes, sizeof(real) );
      if (NULL==neightab_r2cut) 
	  error("cannot allocate memory for neightab_r2cut");
  }
  for (i=0; i<ntypes; i++)
      for (j=0; j<ntypes; j++)
	  neightab_r2cut[i*ntypes+j] = 
	MAX( neightab_r2cut[i*ntypes+j], meam_r2_cut[i][j] );
}


