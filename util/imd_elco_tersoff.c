
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
* imd_elco_tersoff 
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "util.h"

/******************************************************************************
*
*  init_tersoff
*
******************************************************************************/

void init_tersoff(void) {

  int i, j, n = 0;
  real tmp;

  /* Parameters for more than one atom type */
  for (i=0; i<ntypes; i++) {
    ter_c2[i] = ters_c[i] * ters_c[i];
    ter_d2[i] = ters_d[i] * ters_d[i];
    for (j=i; j<ntypes; j++) {
      ter_r_cut[i][j]  = ter_r_cut[j][i]  = ters_r_cut[n];
      ter_r2_cut[i][j] = ter_r2_cut[j][i] = ter_r_cut[i][j] * ter_r_cut[i][j];
      ter_r0[i][j]     = ter_r0[j][i]     = ters_r0[n];
      ter_a[i][j]      = ter_a[j][i]      = ters_a[n];
      ter_b[i][j]      = ter_b[j][i]      = ters_b[n];
      ter_la[i][j]     = ter_la[j][i]     = ters_la[n];
      ter_mu[i][j]     = ter_mu[j][i]     = ters_mu[n];
      ++n;      
    }
  }

  for (i=0; i<ntypes; i++) ter_chi[i][i] = 1.0;
  if ( ntypes>1 ) {
    for (i=0; i<(ntypes-1); i++)
      for (j=(i+1); j<ntypes; j++) {
        ter_chi[i][j] = ter_chi[j][i] 
                      = ters_chi[i * ( 2 * ntypes - i - 3 ) / 2 + j - 1]; 
      }
  }

  for (i=0; i<ntypes; i++) ter_om[i][i] = 1.0;
  if ( ntypes>1 ) {
    for (i=0; i<(ntypes-1); i++)
      for (j=(i+1); j<ntypes; j++) {
        ter_om[i][j] = ter_om[j][i] 
                     = ters_om[i * ( 2 * ntypes - i - 3 ) / 2 + j - 1]; 
      }
  }
 
  if( r_cell != -1.0)
    r2_cut = SQR(r_cell);
  else { 
    tmp = 0.0;
    for (i=0; i<ntypes; ++i)
      for (j=0; j<ntypes; ++j)
	tmp = MAX( tmp, ter_r2_cut[i][j] );
    r2_cut = MAX(r2_cut,tmp);
  }

}

/******************************************************************************
*
*  do_elco_tersoff -- computes stress tensor and elastic constants using the 
*             Tersoff potential
*
******************************************************************************/

void do_elco_tersoff(void)
{
  cell   *p;
  int    ic, jc, kc;
  real   *r, *fc, *dfc, *ddfc, *dddfc;
  vektor *d;
  neightab *neigh;
  vektor dcos_j, dcos_k, dzeta_j, dphi_j;
  tensor ddcos_jj, ddcos_jk, ddcos_kk, ddzeta_jj, ddphi_jj, *ddzeta_kk;
  real   *dzeta_bk, *ddzeta_bkl, *dddzeta_bklm;
  real   *db_k, *ddb_lk, *dddb_lkm;
  vektor *dzeta_k;
  tensor *ddphi_jk, *ddphi_lk, *ddzeta_jk;
  cell   *jcell;
  int    i, j, k, l, m, p_typ, j_typ, k_typ, jnum;
  real   *tmpptr;
  real   pot_zwi, tmp_grad;
  real   cos_theta, cut_tmp, cut_tmp_j;
  real   zeta, g_theta, b_ij;
  real   tmp_jj, tmp_jk, tmp_kk, tmp_j2, tmp_k2;
  real   phi_r, phi_a;
  real   dphi_r, dphi_a, ddphi_r, ddphi_a, dddphi_r, dddphi_a;
  real   tmp_1, tmp_2, tmp_3, tmp_4, tmp_5, tmp_6;
  real   tmp1_zeta, tmp2_zeta, tmp3_zeta, tmp4_zeta, tmp5_zeta;
  real   tmp1_phi, tmp2_phi, tmp3_phi, tmp4_phi, tmp5_phi;
  real   tmp6_phi, tmp7_phi, tmp8_phi, tmp9_phi;
  real   tmp_b1, tmp_b2, tmp_b3;
  real   tmp;

  d            = (vektor *) malloc( neigh_len                         * sizeof(vektor) );
  r            = (real *)   malloc( neigh_len                         * sizeof(real)   );
  fc           = (real *)   malloc( neigh_len                         * sizeof(real)   );
  dfc          = (real *)   malloc( neigh_len                         * sizeof(real)   );
  ddfc         = (real *)   malloc( neigh_len                         * sizeof(real)   );
  dddfc        = (real *)   malloc( neigh_len                         * sizeof(real)   );
  dzeta_k      = (vektor *) malloc( neigh_len                         * sizeof(vektor) );
  dzeta_bk     = (real   *) malloc( neigh_len                         * sizeof(real) );
  ddzeta_bkl   = (real   *) malloc( neigh_len                         * sizeof(real) );
  dddzeta_bklm = (real   *) malloc( neigh_len                         * sizeof(real) );
  db_k         = (real   *) malloc( neigh_len                         * sizeof(real) );
  ddb_lk       = (real   *) malloc( neigh_len * neigh_len             * sizeof(real) );
  dddb_lkm     = (real   *) malloc( neigh_len * neigh_len * neigh_len * sizeof(real) );
  ddphi_jk     = (tensor *) malloc( neigh_len                         * sizeof(tensor) );
  ddphi_lk     = (tensor *) malloc( neigh_len * neigh_len             * sizeof(tensor) );
  ddzeta_jk    = (tensor *) malloc( neigh_len                         * sizeof(tensor) );
  ddzeta_kk    = (tensor *) malloc( neigh_len                         * sizeof(tensor) );

  if ((d==NULL) || (r==NULL) || (fc==NULL) || (dfc==NULL) 
      || (ddfc==NULL) || (dddfc==NULL)
      || (dzeta_k==NULL) || (dzeta_bk==NULL) || (ddzeta_bkl==NULL) 
      || (ddzeta_jk==NULL) || (ddzeta_kk==NULL) 
      || (db_k==NULL) || (ddb_lk==NULL) || (dddb_lkm==NULL) 
      || (ddphi_jk==NULL) || (ddphi_lk==NULL) )
    error("cannot allocate memory for temporary neighbor data");

  /*     k
          \
           \
	    i----j  */

  /* Initializations */
  for (i=0; i<neigh_len; i++)
    for (j=0; j<neigh_len; ++j) 
    {
      ddphi_lk I(i,j) .xx = 0.0;
      ddphi_lk I(i,j) .xy = 0.0;
      ddphi_lk I(i,j) .yy = 0.0;
      ddphi_lk I(i,j) .yz = 0.0;
      ddphi_lk I(i,j) .zx = 0.0;
      ddphi_lk I(i,j) .zz = 0.0;
    }

  /* For each cell */
  for (ic=0; ic < cell_dim.x; ++ic)
    for (jc=0; jc < cell_dim.y; ++jc)
      for (kc=0; kc < cell_dim.z; ++kc)
      {
        p = PTR_3D_V(cell_array,ic,jc,kc,cell_dim);

	/* For each atom in cell */
	for (i=0; i<p->n; ++i) 
	{

	  p_typ   = p->sorte[i];
	  neigh   = &p->neightab_array[i];

	  /* Construct some data for all neighbors */
	  tmpptr = neigh->dist;
	  for (j=0; j<neigh->n; ++j) 
	  {

	    /* Type, distance vector, radii */
	    j_typ   = neigh->typ[j];
	    d[j].x  = *tmpptr++;
	    d[j].y  = *tmpptr++;
	    d[j].z  = *tmpptr++;
	    r[j]    = sqrt(SPROD(d[j],d[j]));

	    /* Cutoff function and its derivatives */
	    cut_tmp   = M_PI / ( ter_r_cut[p_typ][j_typ] - ter_r0[p_typ][j_typ] );
	    cut_tmp_j = cut_tmp * ( r[j] - ter_r0[p_typ][j_typ] );
	    if ( r[j] < ter_r0[p_typ][j_typ] ) {
	      fc[j]    = 1.0; 
	      dfc[j]   = 0.0;
	      ddfc[j]  = 0.0;
	      dddfc[j] = 0.0;
	    }
	    else if ( r[j] > ter_r_cut[p_typ][j_typ] ) {
	      fc[j]    = 0.0;
	      dfc[j]   = 0.0;
	      ddfc[j]  = 0.0;
	      dddfc[j] = 0.0;
	    }
	    else {
	      fc[j]    =   0.5 * ( 1.0 + cos( cut_tmp_j ) );
	      dfc[j]   = - 0.5 * cut_tmp * sin( cut_tmp_j );
	      ddfc[j]  = - 0.5 * cut_tmp * cut_tmp * cos( cut_tmp_j );
	      dddfc[j] = - cut_tmp * cut_tmp * dfc[j];
	    }      
	  } /* j */

          /*********************************************************************/

	  /* For each neighbor of i */
	  for (j=0; j<neigh->n; ++j)
	  {

	    j_typ  = neigh->typ[j];
	    jcell  = (cell *) neigh->cl [j];
	    jnum   = neigh->num[j];
	    tmp_jj = 1 / ( r[j] * r[j] );

	    /* Initializations */
	    zeta         = 0.0;     
	    dzeta_j.x    = 0.0; dzeta_j.y    = 0.0; dzeta_j.z    = 0.0;
	    ddzeta_jj.xx = 0.0; ddzeta_jj.xy = 0.0; ddzeta_jj.yy = 0.0;
	    ddzeta_jj.yz = 0.0; ddzeta_jj.zx = 0.0; ddzeta_jj.zz = 0.0;

	    /* For each neighbor of i other than j */
	    for (k=0; k<neigh->n; ++k) if (k!=j) 
	    {

	      k_typ = neigh->typ[k];

	      /* Angular term */
	      tmp_jk    = 1 / ( r[j] * r[k] );  
	      cos_theta = SPROD(d[j],d[k]) * tmp_jk;
	      tmp_1     = ters_h[p_typ] - cos_theta;
	      tmp_2     = 1 / ( ter_d2[p_typ] + tmp_1 * tmp_1 );
	      g_theta   = 1 + ter_c2[p_typ]/ter_d2[p_typ] - ter_c2[p_typ] * tmp_2;
	      
	      /* zeta */
	      zeta  += fc[k] * ter_om[p_typ][k_typ] * g_theta; 

	      /* tmp variables */
	      tmp_j2 = cos_theta / ( r[j] * r[j] );
	      tmp_k2 = cos_theta / ( r[k] * r[k] );
	      tmp_kk = 1 /  ( r[k] * r[k] );
	      
	      /* Derivatives of cos(theta) */
	      dcos_j.x = tmp_jk * d[k].x - tmp_j2 * d[j].x;
	      dcos_j.y = tmp_jk * d[k].y - tmp_j2 * d[j].y;
	      dcos_j.z = tmp_jk * d[k].z - tmp_j2 * d[j].z;
	      
	      dcos_k.x = tmp_jk * d[j].x - tmp_k2 * d[k].x;
	      dcos_k.y = tmp_jk * d[j].y - tmp_k2 * d[k].y;
	      dcos_k.z = tmp_jk * d[j].z - tmp_k2 * d[k].z;
	      

	      ddcos_jj.xx = - 2 * tmp_jk * tmp_jj * d[j].x * d[k].x 
		            + 3 * tmp_j2 * tmp_jj * d[j].x * d[j].x - tmp_j2;

	      ddcos_jj.xy = - tmp_jk * tmp_jj * ( d[j].x * d[k].y + d[j].y * d[k].x )
		            + 3 * tmp_j2 * tmp_jj * d[j].x * d[j].y;

	      ddcos_jj.yy = - 2 * tmp_jk * tmp_jj * d[j].y * d[k].y 
		            + 3 * tmp_j2 * tmp_jj * d[j].y * d[j].y - tmp_j2;

	      ddcos_jj.yz = - tmp_jk * tmp_jj * ( d[j].y * d[k].z + d[j].z * d[k].y )
		            + 3 * tmp_j2 * tmp_jj * d[j].y * d[j].z;

	      ddcos_jj.zx = - tmp_jk * tmp_jj * ( d[j].z * d[k].x + d[j].x * d[k].z )
		            + 3 * tmp_j2 * tmp_jj * d[j].z * d[j].x;

	      ddcos_jj.zz = - 2 * tmp_jk * tmp_jj * d[j].z * d[k].z 
		            + 3 * tmp_j2 * tmp_jj * d[j].z * d[j].z - tmp_j2;


	      ddcos_jk.xx = - tmp_jj * tmp_jk * d[j].x * d[j].x 
		            + tmp_k2 * tmp_jj * d[j].x * d[k].x
		            - tmp_kk * tmp_jk * d[k].x * d[k].x + tmp_jk;

	      ddcos_jk.xy = - tmp_jj * tmp_jk * d[j].x * d[j].y 
		            + tmp_k2 * tmp_jj * d[j].x * d[k].y
		            - tmp_kk * tmp_jk * d[k].x * d[k].y;

	      ddcos_jk.xz = - tmp_jj * tmp_jk * d[j].x * d[j].z 
		            + tmp_k2 * tmp_jj * d[j].x * d[k].z
		            - tmp_kk * tmp_jk * d[k].x * d[k].z;

	      ddcos_jk.yx = - tmp_jj * tmp_jk * d[j].y * d[j].x 
		            + tmp_k2 * tmp_jj * d[j].y * d[k].x
		            - tmp_kk * tmp_jk * d[k].y * d[k].x;

	      ddcos_jk.yy = - tmp_jj * tmp_jk * d[j].y * d[j].y 
		            + tmp_k2 * tmp_jj * d[j].y * d[k].y
		            - tmp_kk * tmp_jk * d[k].y * d[k].y + tmp_jk;

	      ddcos_jk.yz = - tmp_jj * tmp_jk * d[j].y * d[j].z 
		              + tmp_k2 * tmp_jj * d[j].y * d[k].z
		              - tmp_kk * tmp_jk * d[k].y * d[k].z;

	      ddcos_jk.zx = - tmp_jj * tmp_jk * d[j].z * d[j].x 
		            + tmp_k2 * tmp_jj * d[j].z * d[k].x
		            - tmp_kk * tmp_jk * d[k].z * d[k].x;

	      ddcos_jk.zy = - tmp_jj * tmp_jk * d[j].z * d[j].y 
		            + tmp_k2 * tmp_jj * d[j].z * d[k].y
		            - tmp_kk * tmp_jk * d[k].z * d[k].y;

	      ddcos_jk.zz = - tmp_jj * tmp_jk * d[j].z * d[j].z 
		            + tmp_k2 * tmp_jj * d[j].z * d[k].z
		            - tmp_kk * tmp_jk * d[k].z * d[k].z + tmp_jk;


	      ddcos_kk.xx = - 2 * tmp_jk * tmp_kk * d[j].x * d[k].x 
		            + 3 * tmp_k2 * tmp_kk * d[k].x * d[k].x - tmp_k2;

	      ddcos_kk.xy = - tmp_jk * tmp_kk * ( d[j].x * d[k].y + d[k].x * d[j].y ) 
		            + 3 * tmp_k2 * tmp_kk * d[k].x * d[k].y;

	      ddcos_kk.yy = - 2 * tmp_jk * tmp_kk * d[j].y * d[k].y 
		            + 3 * tmp_k2 * tmp_kk * d[k].y * d[k].y - tmp_k2;

	      ddcos_kk.yz = - tmp_jk * tmp_kk * ( d[j].y * d[k].z + d[k].y * d[j].z ) 
		            + 3 * tmp_k2 * tmp_kk * d[k].y * d[k].z;

	      ddcos_kk.zx = - tmp_jk * tmp_kk * ( d[j].z * d[k].x + d[k].z * d[j].x ) 
		            + 3 * tmp_k2 * tmp_kk * d[k].z * d[k].x;

	      ddcos_kk.zz = - 2 * tmp_jk * tmp_kk * d[j].z * d[k].z 
		            + 3 * tmp_k2 * tmp_kk * d[k].z * d[k].z - tmp_k2;

	      /* tmp variables for derivatives of zeta */
	      tmp_3     = 2 * ter_c2[p_typ] * tmp_1 * tmp_2 * tmp_2 
		            * fc[k] * ter_om[p_typ][k_typ];
	      tmp_grad  = dfc[k] / r[k] * g_theta * ter_om[p_typ][k_typ];
	      tmp1_zeta = ter_om[p_typ][k_typ] * fc[k] * 2 * ter_c2[p_typ] * tmp_2 * tmp_2;
	      tmp2_zeta = tmp_2 * ( ter_d2[p_typ] - 3 * tmp_1 * tmp_1 );
	      tmp3_zeta = - ter_om[p_typ][k_typ] * dfc[k] / r[k] 
		          * 2 * ter_c2[p_typ] * tmp_1 * tmp_2 * tmp_2;
	      tmp4_zeta = ter_om[p_typ][k_typ] * tmp_kk * ( ddfc[k] - dfc[k] / r[k] ) * g_theta;
	      tmp5_zeta = ter_om[p_typ][k_typ] * dfc[k] / r[k] * g_theta;

	      /* First and second derivatives of zeta */
	      dzeta_bk[k]  = dfc[k] * g_theta;          /* For B' */

	      dzeta_j.x   -= tmp_3 * dcos_j.x;
	      dzeta_j.y   -= tmp_3 * dcos_j.y;
	      dzeta_j.z   -= tmp_3 * dcos_j.z;

	      dzeta_k[k].x = tmp_grad * d[k].x - tmp_3 * dcos_k.x;
	      dzeta_k[k].y = tmp_grad * d[k].y - tmp_3 * dcos_k.y;
	      dzeta_k[k].z = tmp_grad * d[k].z - tmp_3 * dcos_k.z;

	      ddzeta_bkl[k] = ddfc[k] * g_theta; /* For B' */

	      ddzeta_jj.xx += tmp1_zeta * ( tmp2_zeta * dcos_j.x * dcos_j.x
		                                - tmp_1 * ddcos_jj.xx );

	      ddzeta_jj.xy += tmp1_zeta * ( tmp2_zeta * dcos_j.x * dcos_j.y
		                                - tmp_1 * ddcos_jj.xy );

	      ddzeta_jj.yy += tmp1_zeta * ( tmp2_zeta * dcos_j.y * dcos_j.y
		                                - tmp_1 * ddcos_jj.yy );

	      ddzeta_jj.yz += tmp1_zeta * ( tmp2_zeta * dcos_j.y * dcos_j.z
		                                - tmp_1 * ddcos_jj.yz );

	      ddzeta_jj.zx += tmp1_zeta * ( tmp2_zeta * dcos_j.z * dcos_j.x
		                                - tmp_1 * ddcos_jj.zx );

	      ddzeta_jj.zz += tmp1_zeta * ( tmp2_zeta * dcos_j.z * dcos_j.z
		                                - tmp_1 * ddcos_jj.zz );

		
	      ddzeta_jk[k].xx = tmp1_zeta * ( tmp2_zeta * dcos_j.x * dcos_k.x
		                                  - tmp_1 * ddcos_jk.xx )
		                              + tmp3_zeta * d[k].x * dcos_j.x; 

	      ddzeta_jk[k].xy = tmp1_zeta * ( tmp2_zeta * dcos_j.x * dcos_k.y
		                                  - tmp_1 * ddcos_jk.xy )
		                              + tmp3_zeta * d[k].y * dcos_j.x; 

	      ddzeta_jk[k].xz = tmp1_zeta * ( tmp2_zeta * dcos_j.x * dcos_k.z
		                                  - tmp_1 * ddcos_jk.xz )
		                              + tmp3_zeta * d[k].z * dcos_j.x; 

	      ddzeta_jk[k].yx = tmp1_zeta * ( tmp2_zeta * dcos_j.y * dcos_k.x
		                                  - tmp_1 * ddcos_jk.yx )
		                              + tmp3_zeta * d[k].x * dcos_j.y; 

	      ddzeta_jk[k].yy = tmp1_zeta * ( tmp2_zeta * dcos_j.y * dcos_k.y
						  - tmp_1 * ddcos_jk.yy )
		                              + tmp3_zeta * d[k].y * dcos_j.y; 

	      ddzeta_jk[k].yz = tmp1_zeta * ( tmp2_zeta * dcos_j.y * dcos_k.z
		                                  - tmp_1 * ddcos_jk.yz )
		                              + tmp3_zeta * d[k].z * dcos_j.y; 

	      ddzeta_jk[k].zx = tmp1_zeta * ( tmp2_zeta * dcos_j.z * dcos_k.x
		                                  - tmp_1 * ddcos_jk.zx )
		                              + tmp3_zeta * d[k].x * dcos_j.z; 

	      ddzeta_jk[k].zy = tmp1_zeta * ( tmp2_zeta * dcos_j.z * dcos_k.y
		                                  - tmp_1 * ddcos_jk.zy )
		                              + tmp3_zeta * d[k].y * dcos_j.z; 

	      ddzeta_jk[k].zz = tmp1_zeta * ( tmp2_zeta * dcos_j.z * dcos_k.z
		                                  - tmp_1 * ddcos_jk.zz )
		                              + tmp3_zeta * d[k].z * dcos_j.z; 


	      ddzeta_kk[k].xx = + tmp4_zeta * d[k].x * d[k].x
		            + tmp3_zeta * 2 * d[k].x * dcos_k.x
		  + tmp1_zeta * ( tmp2_zeta * dcos_k.x * dcos_k.x - tmp_1 * ddcos_kk.xx )
		  + tmp5_zeta;

	      ddzeta_kk[k].xy = tmp4_zeta * d[k].x * d[k].y 
		            + tmp3_zeta * ( d[k].x * dcos_k.y + d[k].y * dcos_k.x )
		+ tmp1_zeta * ( tmp2_zeta * dcos_k.x * dcos_k.y - tmp_1 * ddcos_kk.xy ) ;

	      ddzeta_kk[k].yy = + tmp4_zeta * d[k].y * d[k].y 
		            + tmp3_zeta * 2 * d[k].y * dcos_k.y
		  + tmp1_zeta * ( tmp2_zeta * dcos_k.y * dcos_k.y - tmp_1 * ddcos_kk.yy )
		  + tmp5_zeta;

	      ddzeta_kk[k].yz = tmp4_zeta * d[k].y * d[k].z 
		            + tmp3_zeta * ( d[k].y * dcos_k.z + d[k].z * dcos_k.y )
		+ tmp1_zeta * ( tmp2_zeta * dcos_k.y * dcos_k.z - tmp_1 * ddcos_kk.yz ) ;

	      ddzeta_kk[k].zx = tmp4_zeta * d[k].z * d[k].x 
		            + tmp3_zeta * ( d[k].z * dcos_k.x + d[k].x * dcos_k.z )
		+ tmp1_zeta * ( tmp2_zeta * dcos_k.z * dcos_k.x - tmp_1 * ddcos_kk.zx ) ;

	      ddzeta_kk[k].zz = + tmp4_zeta * d[k].z * d[k].z 
		            + tmp3_zeta * 2 * d[k].z * dcos_k.z
		  + tmp1_zeta * ( tmp2_zeta * dcos_k.z * dcos_k.z - tmp_1 * ddcos_kk.zz )
		  + tmp5_zeta;

	      /* Third derivative of zeta, for B' */
	      dddzeta_bklm[k] = dddfc[k] * g_theta;

	    } /* k */

	    phi_r  = 0.5 * ter_a[p_typ][j_typ] * exp( - ter_la[p_typ][j_typ] * r[j] );
	    phi_a  = 0.5 * ter_b[p_typ][j_typ] * exp( - ter_mu[p_typ][j_typ] * r[j] );
	    tmp_4  = pow( ters_ga[p_typ] * zeta, ters_n[p_typ] );

	    b_ij  = ter_chi[p_typ][j_typ] * pow( 1 + tmp_4, - 1 / ( 2 * ters_n[p_typ] ) );

	    pot_zwi  = phi_r - b_ij * phi_a;

	    epot += fc[j] * pot_zwi;

	    if ( zeta == 0.0 )   /* only one neighbor of i */
	      tmp_5 = 0.0;
	    else
	      tmp_5 = - b_ij * fc[j] * phi_a * tmp_4 / ( 2 * zeta * ( 1 + tmp_4 ) );

	    tmp_6   = - ( fc[j] * ( - phi_r * ter_la[p_typ][j_typ] 
			+ phi_a * ter_mu[p_typ][j_typ] * b_ij )  + dfc[j] * pot_zwi ) / r[j];

	    /* First derivatives of phi */
	    dphi_r   = - ter_la[p_typ][j_typ] * phi_r;
	    dphi_a   = - ter_mu[p_typ][j_typ] * phi_a;

	    dphi_j.x = tmp_6 * d[j].x + tmp_5 * dzeta_j.x;
	    dphi_j.y = tmp_6 * d[j].y + tmp_5 * dzeta_j.y;
	    dphi_j.z = tmp_6 * d[j].z + tmp_5 * dzeta_j.z;

	    /* tmp variables for derivatives of phi */
	    tmp1_phi = tmp_jj * ( ddfc[j] - dfc[j] / r[j] ) * pot_zwi;
	    tmp2_phi = dfc[j] / r[j] * pot_zwi;
	    tmp3_phi = 2 * dfc[j] * tmp_jj * (- phi_r * ter_la[p_typ][j_typ] 
					      + phi_a * ter_mu[p_typ][j_typ] * b_ij);
	    tmp4_phi = fc[j] * ter_la[p_typ][j_typ] * phi_r / r[j];
	    tmp5_phi = - fc[j] * b_ij * ter_mu[p_typ][j_typ] * phi_a / r[j];
	    tmp6_phi = ter_la[p_typ][j_typ] / r[j] + tmp_jj;
	    tmp7_phi = ter_mu[p_typ][j_typ] / r[j] + tmp_jj;
	    if ( zeta == 0.0 ) {
	      tmp8_phi = 0.0;
	      tmp9_phi = 0.0;
	    }
	    else {
	      tmp8_phi = phi_a / r[j] * ( dfc[j] - fc[j] * ter_mu[p_typ][j_typ] ) 
	               * b_ij * tmp_4 / ( 2.0 * zeta * ( 1.0 + tmp_4 ) );
	      tmp9_phi = - tmp_5 * ( ters_n[p_typ] - 1.0 - 1.5 * tmp_4 ) / ( zeta * ( 1.0 + tmp_4 ) );
	    }

	    /* tmp variables for computation of B' */
	    if ( zeta == 0.0 ) {   /* only one neighbor of i */
	      tmp_b1 = 0.0;
	      tmp_b2 = 0.0;
	      tmp_b3 = 0.0;
	    }
	    else {
	      tmp_b1 = - b_ij * tmp_4 / ( 2.0 * zeta * ( 1.0 + tmp_4 ) );

	      tmp_b2 = tmp_b1 * ( ters_n[p_typ] - 1.0 - 1.5 * tmp_4 ) / ( zeta * ( 1.0 + tmp_4 ) );

	      tmp_b3 = tmp_b1 * ( - 2.0 - ters_n[p_typ] * ters_n[p_typ] + 3.0 * ters_n[p_typ]  
		 + tmp_4 * ( 9.0 * ters_n[p_typ] - 11 ) / 2.0 
		 + tmp_4 * tmp_4 * ( 2.0 * ters_n[p_typ] - 15.0 ) / 4.0 ) 
	      / (zeta * zeta * ( 1.0 + tmp_4 ) * ( 1.0 + tmp_4 ) );
	    }

	    /* Second derivatives of phi */
	    ddphi_r = - ter_la[p_typ][j_typ] * dphi_r;
	    ddphi_a = - ter_mu[p_typ][j_typ] * dphi_a;

	    ddphi_jj.xx = tmp1_phi * d[j].x * d[j].x + tmp2_phi
	                + tmp3_phi * d[j].x * d[j].x
	      + tmp4_phi * ( d[j].x * d[j].x * tmp6_phi - 1 )
	      + tmp5_phi * ( d[j].x * d[j].x * tmp7_phi - 1 ) 
	      + 2 * tmp8_phi * dzeta_j.x * d[j].x 
	      + tmp9_phi * dzeta_j.x * dzeta_j.x - tmp_5 * ddzeta_jj.xx;

	    ddphi_jj.xy = tmp1_phi* d[j].x * d[j].y 
	               + tmp3_phi * d[j].x * d[j].y 
	      + ( tmp4_phi * tmp6_phi + tmp5_phi * tmp7_phi ) * d[j].x * d[j].y 
	      + tmp8_phi * ( dzeta_j.x * d[j].y + dzeta_j.y * d[j].x ) 
	      + tmp9_phi * dzeta_j.x * dzeta_j.y - tmp_5 * ddzeta_jj.xy;

	    ddphi_jj.yy = tmp1_phi * d[j].y * d[j].y + tmp2_phi
	                + tmp3_phi * d[j].y * d[j].y
	      + tmp4_phi * ( d[j].y * d[j].y * tmp6_phi - 1 )
	      + tmp5_phi * ( d[j].y * d[j].y * tmp7_phi - 1 ) 
	      + 2 * tmp8_phi * dzeta_j.y * d[j].y 
	      + tmp9_phi * dzeta_j.y * dzeta_j.y - tmp_5 * ddzeta_jj.yy;

	    ddphi_jj.yz = tmp1_phi* d[j].y * d[j].z 
	               + tmp3_phi * d[j].y * d[j].z 
	      + ( tmp4_phi * tmp6_phi + tmp5_phi * tmp7_phi ) * d[j].y * d[j].z 
	      + tmp8_phi * ( dzeta_j.y * d[j].z + dzeta_j.z * d[j].y ) 
	      + tmp9_phi * dzeta_j.y * dzeta_j.z - tmp_5 * ddzeta_jj.yz;

	    ddphi_jj.zx = tmp1_phi* d[j].z * d[j].x 
	               + tmp3_phi * d[j].z * d[j].x 
	      + ( tmp4_phi * tmp6_phi + tmp5_phi * tmp7_phi ) * d[j].z * d[j].x 
	      + tmp8_phi * ( dzeta_j.z * d[j].x + dzeta_j.x * d[j].z ) 
	      + tmp9_phi * dzeta_j.z * dzeta_j.x - tmp_5 * ddzeta_jj.zx;

	    ddphi_jj.zz = tmp1_phi * d[j].z * d[j].z + tmp2_phi
	                + tmp3_phi * d[j].z * d[j].z
	      + tmp4_phi * ( d[j].z * d[j].z * tmp6_phi - 1 )
	      + tmp5_phi * ( d[j].z * d[j].z * tmp7_phi - 1 ) 
	      + 2 * tmp8_phi * dzeta_j.z * d[j].z 
	      + tmp9_phi * dzeta_j.z * dzeta_j.z - tmp_5 * ddzeta_jj.zz;

	    for (k=0; k<neigh->n; ++k) if (k!=j) 
	    {

	      ddphi_jk[k].xx = tmp8_phi * dzeta_k[k].x * d[j].x
		+ tmp9_phi * dzeta_j.x * dzeta_k[k].x - tmp_5 * ddzeta_jk[k].xx;

	      ddphi_jk[k].xy = tmp8_phi * dzeta_k[k].y * d[j].x
		+ tmp9_phi * dzeta_j.x * dzeta_k[k].y - tmp_5 * ddzeta_jk[k].xy;

	      ddphi_jk[k].xz = tmp8_phi * dzeta_k[k].z * d[j].x
		+ tmp9_phi * dzeta_j.x * dzeta_k[k].z - tmp_5 * ddzeta_jk[k].xz;

	      ddphi_jk[k].yx = tmp8_phi * dzeta_k[k].x * d[j].y
		+ tmp9_phi * dzeta_j.y * dzeta_k[k].x - tmp_5 * ddzeta_jk[k].yx;

	      ddphi_jk[k].yy = tmp8_phi * dzeta_k[k].y * d[j].y
		+ tmp9_phi * dzeta_j.y * dzeta_k[k].y - tmp_5 * ddzeta_jk[k].yy;

	      ddphi_jk[k].yz = tmp8_phi * dzeta_k[k].z * d[j].y
		+ tmp9_phi * dzeta_j.y * dzeta_k[k].z - tmp_5 * ddzeta_jk[k].yz;

	      ddphi_jk[k].zx = tmp8_phi * dzeta_k[k].x * d[j].z
		+ tmp9_phi * dzeta_j.z * dzeta_k[k].x - tmp_5 * ddzeta_jk[k].zx;

	      ddphi_jk[k].zy = tmp8_phi * dzeta_k[k].y * d[j].z
		+ tmp9_phi * dzeta_j.z * dzeta_k[k].y - tmp_5 * ddzeta_jk[k].zy;

	      ddphi_jk[k].zz = tmp8_phi * dzeta_k[k].z * d[j].z
		+ tmp9_phi * dzeta_j.z * dzeta_k[k].z - tmp_5 * ddzeta_jk[k].zz;

	      /* First derivative of b_ij , for B' */	      
	      db_k[k] = tmp_b1 * dzeta_bk[k];      

	      for(l=0; l<neigh->n; ++l) if (l!=j) 
	      {

		ddphi_lk I(l,k) .xx = tmp9_phi * dzeta_k[l].x * dzeta_k[k].x;
		ddphi_lk I(l,k) .xy = tmp9_phi * dzeta_k[l].x * dzeta_k[k].y;
		ddphi_lk I(l,k) .yy = tmp9_phi * dzeta_k[l].y * dzeta_k[k].y;
		ddphi_lk I(l,k) .yz = tmp9_phi * dzeta_k[l].y * dzeta_k[k].z;
		ddphi_lk I(l,k) .zx = tmp9_phi * dzeta_k[l].z * dzeta_k[k].x;
		ddphi_lk I(l,k) .zz = tmp9_phi * dzeta_k[l].z * dzeta_k[k].z;

		/* Second derivative of b_ij, for B' */
		ddb_lk I(l,k) = tmp_b2 * dzeta_bk[k] * dzeta_bk[l];
					       
		if ( l==k ) {
 
		  ddphi_lk I(l,k) .xx += - tmp_5 * ddzeta_kk[k].xx; 
		  ddphi_lk I(l,k) .xy += - tmp_5 * ddzeta_kk[k].xy; 
		  ddphi_lk I(l,k) .yy += - tmp_5 * ddzeta_kk[k].yy; 
		  ddphi_lk I(l,k) .yz += - tmp_5 * ddzeta_kk[k].yz; 
		  ddphi_lk I(l,k) .zx += - tmp_5 * ddzeta_kk[k].zx; 
		  ddphi_lk I(l,k) .zz += - tmp_5 * ddzeta_kk[k].zz;

		  /* Second derivative of b_ij, for B' */
		  ddb_lk I(l,k) += - tmp_b1 * ddzeta_bkl[k];
		  
		}

		for(m=0; m<neigh->n; ++m) if (m!=j) 
		{
		  dddb_lkm J(l,k,m) = tmp_b3 * dzeta_bk[l] * dzeta_bk[k] * dzeta_bk[m];

		  if ( k==l )
		    dddb_lkm J(l,k,m) += tmp_b2 * dzeta_bk[m] * ddzeta_bkl[k];
		  if ( k==m )
		    dddb_lkm J(l,k,m) += tmp_b2 * dzeta_bk[l] * ddzeta_bkl[k];
		  if ( l==m ) {
		    dddb_lkm J(l,k,m) += tmp_b2 * dzeta_bk[k] * ddzeta_bkl[l];

		    if ( k==l )
		      dddb_lkm J(l,k,m) += - tmp_b1 * dddzeta_bklm[l];
		  }
		}

	      } /* l */
	    } /* k */

	    /* Third derivatives of phi */
	    dddphi_r = - ter_la[p_typ][j_typ] * ddphi_r;
	    dddphi_a = - ter_mu[p_typ][j_typ] * ddphi_a;

	    /* Compute stress and elastic constants, contribution from j */
	    tmp = 0.5 * d[j].x * dphi_j.x;
	    p->stress[i].xx        -= tmp;
	    jcell->stress[jnum].xx -= tmp;
	    sigma.xx               -= 2 * tmp;
	    tmp = 0.5 * d[j].y * dphi_j.y;
	    p->stress[i].yy        -= tmp;
	    jcell->stress[jnum].yy -= tmp;
	    sigma.yy               -= 2 * tmp;
	    tmp = 0.5 * d[j].z * dphi_j.z;
	    p->stress[i].zz        -= tmp;
	    jcell->stress[jnum].zz -= tmp;
	    sigma.zz               -= 2 * tmp;
	    tmp = 0.25 * ( d[j].y * dphi_j.z + d[j].z * dphi_j.y );
	    p->stress[i].yz        -= tmp;
	    jcell->stress[jnum].yz -= tmp;
	    sigma.yz               -= 2 * tmp;
	    tmp = 0.25 * ( d[j].z * dphi_j.x + d[j].x * dphi_j.z );
	    p->stress[i].zx        -= tmp;
	    jcell->stress[jnum].zx -= tmp;
	    sigma.zx               -= 2 * tmp;
	    tmp = 0.25 * ( d[j].x * dphi_j.y + d[j].y * dphi_j.x );
	    p->stress[i].xy        -= tmp;
	    jcell->stress[jnum].xy -= tmp;
	    sigma.xy               -= 2 * tmp;

	    tmp = 0.5   * ddphi_jj.xx * d[j].x * d[j].x;
	    p->elco[i].c11        += tmp;
	    jcell->elco[jnum].c11 += tmp;
	    c.c11                 += 2 * tmp;
	    tmp = 0.5   * ddphi_jj.xy * d[j].x * d[j].y;
	    p->elco[i].c12        += tmp;
	    jcell->elco[jnum].c12 += tmp;
	    c.c12                 += 2 * tmp;
	    tmp = 0.5   * ddphi_jj.zx * d[j].x * d[j].z;
	    p->elco[i].c13        += tmp;
	    jcell->elco[jnum].c13 += tmp;
	    c.c13                 += 2 * tmp;
	    tmp = 0.5   * ddphi_jj.yy * d[j].y * d[j].y;
	    p->elco[i].c22        += tmp;
	    jcell->elco[jnum].c22 += tmp;
	    c.c22                 += 2 * tmp;
	    tmp = 0.5   * ddphi_jj.yz * d[j].y * d[j].z;
	    p->elco[i].c23        += tmp;
	    jcell->elco[jnum].c23 += tmp;
	    c.c23                 += 2 * tmp;
	    tmp = 0.5   * ddphi_jj.zz * d[j].z * d[j].z;
	    p->elco[i].c33        += tmp;
	    jcell->elco[jnum].c33 += tmp;
	    c.c33                 += 2 * tmp;
	    tmp = 0.125 * ( ddphi_jj.xx * d[j].y * d[j].y 
	                  + ddphi_jj.yy * d[j].x * d[j].x
	              + 2 * ddphi_jj.xy * d[j].x * d[j].y );
	    p->elco[i].c44        += tmp;
	    jcell->elco[jnum].c44 += tmp;
	    c.c44                 += 2 * tmp;
	    tmp = 0.125 * ( ddphi_jj.yz * d[j].z * d[j].x
	                  + ddphi_jj.zz * d[j].y * d[j].x
	                  + ddphi_jj.xy * d[j].z * d[j].z
	                  + ddphi_jj.zx * d[j].y * d[j].z );
	    p->elco[i].c45        += tmp;
	    jcell->elco[jnum].c45 += tmp;
	    c.c45                 += 2 * tmp;
	    tmp = 0.125 * ( ddphi_jj.xy * d[j].y * d[j].z
	                  + ddphi_jj.yy * d[j].x * d[j].z
	                  + ddphi_jj.zx * d[j].y * d[j].y
	                  + ddphi_jj.yz * d[j].x * d[j].y );
	    p->elco[i].c46        += tmp;
	    jcell->elco[jnum].c46 += tmp;
	    c.c46                 += 2 * tmp;
	    tmp = 0.125 * ( ddphi_jj.zz * d[j].x * d[j].x 
	                  + ddphi_jj.xx * d[j].z * d[j].z
	              + 2 * ddphi_jj.zx * d[j].z * d[j].x );
	    p->elco[i].c55        += tmp;
	    jcell->elco[jnum].c55 += tmp;
	    c.c55                 += 2 * tmp;
	    tmp = 0.125 * ( ddphi_jj.zx * d[j].x * d[j].y
	                    + ddphi_jj.xx * d[j].z * d[j].y
	                    + ddphi_jj.yz * d[j].x * d[j].x
	                    + ddphi_jj.xy * d[j].z * d[j].x );
	    p->elco[i].c56        += tmp;
	    jcell->elco[jnum].c56 += tmp;
	    c.c56                 += 2 * tmp;
	    tmp = 0.125 * ( ddphi_jj.xx * d[j].y * d[j].y 
	                  + ddphi_jj.yy * d[j].x * d[j].x
	              + 2 * ddphi_jj.xy * d[j].x * d[j].y );
	    p->elco[i].c66        += tmp;
	    jcell->elco[jnum].c66 += tmp;
	    c.c66                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.xy * d[j].x * d[j].z
	                 + ddphi_jj.zx * d[j].x * d[j].y );
	    p->elco[i].c14        += tmp;
	    jcell->elco[jnum].c14 += tmp;
	    c.c14                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.zx * d[j].x * d[j].x
	                 + ddphi_jj.xx * d[j].z * d[j].x );
	    p->elco[i].c15        += tmp;
	    jcell->elco[jnum].c15 += tmp;
	    c.c15                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.xx * d[j].x * d[j].y
	                 + ddphi_jj.xy * d[j].x * d[j].x );
	    p->elco[i].c16        += tmp;
	    jcell->elco[jnum].c16 += tmp;
	    c.c16                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.yy * d[j].y * d[j].z
	                 + ddphi_jj.yz * d[j].y * d[j].y );
	    p->elco[i].c24        += tmp;
	    jcell->elco[jnum].c24 += tmp;
	    c.c24                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.yz * d[j].y * d[j].x
	                 + ddphi_jj.xy * d[j].y * d[j].z );
	    p->elco[i].c25        += tmp;
	    jcell->elco[jnum].c25 += tmp;
	    c.c25                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.xy * d[j].y * d[j].y
	                 + ddphi_jj.yy * d[j].x * d[j].y );
	    p->elco[i].c26        += tmp;
	    jcell->elco[jnum].c26 += tmp;
	    c.c26                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.yz * d[j].z * d[j].z
	                 + ddphi_jj.zz * d[j].y * d[j].z );
	    p->elco[i].c34        += tmp;
	    jcell->elco[jnum].c34 += tmp;
	    c.c34                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.zz * d[j].z * d[j].x
	                 + ddphi_jj.zx * d[j].z * d[j].z );
	    p->elco[i].c35        += tmp;
	    jcell->elco[jnum].c35 += tmp;
	    c.c35                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.zx * d[j].z * d[j].y
	                    + ddphi_jj.yz * d[j].z * d[j].x );
	    p->elco[i].c36        += tmp;
	    jcell->elco[jnum].c36 += tmp;
	    c.c36                 += 2 * tmp;	    

	    /* B', contribution from j */
	    dbulkm_dp += ( dddfc[j] * ( phi_r - b_ij * phi_a ) 
	      + 3.0 * ddfc[j] * ( dphi_r - b_ij * dphi_a )
	      + 3.0 * dfc[j] * ( ddphi_r - b_ij * ddphi_a )
	      + fc[j] * ( dddphi_r - b_ij * dddphi_a ) ) 
	      * r[j] * r[j] * r[j];

	    /* Compute stress and elastic constants, contributions from k */
	    for (k=0; k<neigh->n; ++k) if (k!=j) 
	    {
	     tmp = 0.5 * d[k].x * tmp_5 * dzeta_k[k].x; 
	     p->stress[i].xx        -= tmp;
	     jcell->stress[jnum].xx -= tmp;
	     sigma.xx               -= 2 * tmp;
	     tmp = 0.5 * d[k].y * tmp_5 * dzeta_k[k].y;
	     p->stress[i].yy        -= tmp;
	     jcell->stress[jnum].yy -= tmp;
	     sigma.yy               -= 2 * tmp;
	     tmp = 0.5 * d[k].z * tmp_5 * dzeta_k[k].z;
	     p->stress[i].zz        -= tmp;
	     jcell->stress[jnum].zz -= tmp;
	     sigma.zz               -= 2 * tmp;
	     tmp = 0.25 * ( d[k].y * tmp_5 * dzeta_k[k].z  
			  + d[k].z * tmp_5 * dzeta_k[k].y );
	     p->stress[i].yz        -= tmp;
	     jcell->stress[jnum].yz -= tmp;
	     sigma.yz               -= 2 * tmp;
	     tmp = 0.25 * ( d[k].z * tmp_5 * dzeta_k[k].x  
			  + d[k].x * tmp_5 * dzeta_k[k].z );
	     p->stress[i].zx        -= tmp;
	     jcell->stress[jnum].zx -= tmp;
	     sigma.zx               -= 2 * tmp;
	     tmp = 0.25 * ( d[k].x * tmp_5 * dzeta_k[k].y  
		          + d[k].y * tmp_5 * dzeta_k[k].x );
	     p->stress[i].xy        -= tmp;
	     jcell->stress[jnum].xy -= tmp;
	     sigma.xy               -= 2 * tmp;

	     tmp = ddphi_jk[k].xx * d[j].x * d[k].x;
	     p->elco[i].c11        += tmp;
	     jcell->elco[jnum].c11 += tmp;
	     c.c11                 += 2 * tmp;
	     tmp = 0.5 * ( ddphi_jk[k].xy * d[j].x * d[k].y 
			   + ddphi_jk[k].yx * d[j].y * d[k].x );
	     p->elco[i].c12        += tmp;
	     jcell->elco[jnum].c12 += tmp;
	     c.c12                 += 2 * tmp;
	     tmp = 0.5 * ( ddphi_jk[k].xz * d[j].x * d[k].z 
			   + ddphi_jk[k].zx * d[j].z * d[k].x );
	     p->elco[i].c13        += tmp;
	     jcell->elco[jnum].c13 += tmp;
	     c.c13                 += 2 * tmp;
	     tmp = ddphi_jk[k].yy * d[j].y * d[k].y;
	     p->elco[i].c22        += tmp;
	     jcell->elco[jnum].c22 += tmp;
	     c.c22                 += 2 * tmp;
	     tmp = 0.5 * ( ddphi_jk[k].yz * d[j].y * d[k].z 
			   + ddphi_jk[k].zy * d[j].z * d[k].y );
	     p->elco[i].c23        += tmp;
	     jcell->elco[jnum].c23 += tmp;
	     c.c23                 += 2 * tmp;
	     tmp = ddphi_jk[k].zz * d[j].z * d[k].z;
	     p->elco[i].c33        += tmp;
	     jcell->elco[jnum].c33 += tmp;
	     c.c33                 += 2 * tmp;
	     tmp = 0.25 * ( ddphi_jk[k].yy * d[j].z * d[k].z
			    + ddphi_jk[k].zz * d[j].y * d[k].y
			    + ddphi_jk[k].yz * d[j].z * d[k].y
			    + ddphi_jk[k].zy * d[j].y * d[k].z );
	     p->elco[i].c44        += tmp;
	     jcell->elco[jnum].c44 += tmp;
	     c.c44                 += 2 * tmp;
	     tmp = 0.125 * ( ddphi_jk[k].yz * d[j].z * d[k].x
			     + ddphi_jk[k].zy * d[j].x * d[k].z
			     + ddphi_jk[k].zz * ( d[j].y * d[k].x + d[j].x * d[j].y )
			     + 2 * ddphi_jk[k].yx * d[j].z * d[k].z
			     + ddphi_jk[k].zx * d[j].y * d[k].z
			     + ddphi_jk[k].xz * d[j].z * d[k].y );
	     p->elco[i].c45        += tmp;
	     jcell->elco[jnum].c45 += tmp;
	     c.c45                 += 2 * tmp;	     
	     tmp = 0.125 * ( ddphi_jk[k].xy * d[j].y * d[k].z
			     + ddphi_jk[k].yx * d[j].z * d[k].y
			     + ddphi_jk[k].yy * ( d[j].x * d[k].z + d[j].z * d[j].x )
			     + 2 * ddphi_jk[k].xz * d[j].y * d[k].y
			     + ddphi_jk[k].yz * d[j].x * d[k].y
			     + ddphi_jk[k].zy * d[j].y * d[k].x );
	     p->elco[i].c46        += tmp;
	     jcell->elco[jnum].c46 += tmp;
	     c.c46                 += 2 * tmp;
	     tmp = 0.25 * ( ddphi_jk[k].zz * d[j].x * d[k].x 
			    + ddphi_jk[k].xx * d[j].z * d[k].z
			    + ddphi_jk[k].zx * d[j].x * d[k].z
			    + ddphi_jk[k].xz * d[j].z * d[k].x );
	     p->elco[i].c55        += tmp;
	     jcell->elco[jnum].c55 += tmp;
	     c.c55                 += 2 * tmp;
	     tmp = 0.125 * ( ddphi_jk[k].zx * d[j].x * d[k].y
			     + ddphi_jk[k].xz * d[j].y * d[k].x
			     + ddphi_jk[k].xx * ( d[j].z * d[k].y + d[j].y * d[j].z )
			     + 2 * ddphi_jk[k].zy * d[j].x * d[k].x
			     + ddphi_jk[k].xy * d[j].z * d[k].x
			     + ddphi_jk[k].yx * d[j].x * d[k].z );
	     p->elco[i].c56        += tmp;
	     jcell->elco[jnum].c56 += tmp;
	     c.c56                 += 2 * tmp;
	     tmp = 0.25 * ( ddphi_jk[k].xx * d[j].y * d[k].y 
			    + ddphi_jk[k].yy * d[j].x * d[k].x
			    + ddphi_jk[k].xy * d[j].y * d[k].x
			    + ddphi_jk[k].yx * d[j].x * d[k].y );
	     p->elco[i].c66        += tmp;
	     jcell->elco[jnum].c66 += tmp;
	     c.c66                 += 2 * tmp;
	     tmp = 0.25 * ( ddphi_jk[k].xy * d[j].x * d[k].z
			    + ddphi_jk[k].yx * d[j].z * d[k].x
			    + ddphi_jk[k].xz * d[j].x * d[k].y
			    + ddphi_jk[k].zx * d[j].y * d[k].x );
	     p->elco[i].c14        += tmp;
	     jcell->elco[jnum].c14 += tmp;
	     c.c14                 += 2 * tmp;
	     tmp = 0.25 * ( ( ddphi_jk[k].xz + ddphi_jk[k].zx ) * d[j].x * d[k].x
			    + ddphi_jk[k].xx * ( d[j].x * d[k].z + d[j].z * d[k].x ) );
	     p->elco[i].c15        += tmp;
	     jcell->elco[jnum].c15 += tmp;
	     c.c15                 += 2 * tmp;
	     tmp = 0.25 * ( ddphi_jk[k].xx * ( d[j].x * d[k].y + d[j].y * d[k].x )
			    + ( ddphi_jk[k].xy + ddphi_jk[k].yx ) * d[j].x * d[k].x );
	     p->elco[i].c16        += tmp;
	     jcell->elco[jnum].c16 += tmp;
	     c.c16                 += 2 * tmp;
	     tmp = 0.25 * ( ddphi_jk[k].yy * ( d[j].y * d[k].z + d[j].z * d[k].y )
			    + ( ddphi_jk[k].yz + ddphi_jk[k].zy ) * d[j].y * d[k].y );
	     p->elco[i].c24        += tmp;
	     jcell->elco[jnum].c24 += tmp;
	     c.c24                 += 2 * tmp;
	     tmp = 0.25 * ( ddphi_jk[k].yz * d[j].y * d[k].x
			    + ddphi_jk[k].zy * d[j].x * d[k].y
			    + ddphi_jk[k].yx * d[j].y * d[k].z
			    + ddphi_jk[k].xy * d[j].z * d[k].y );
	     p->elco[i].c25        += tmp;
	     jcell->elco[jnum].c25 += tmp;
	     c.c25                 += 2 * tmp;
	     tmp = 0.25 * ( ( ddphi_jk[k].yx + ddphi_jk[k].xy ) * d[j].y * d[k].y
			    + ddphi_jk[k].yy * ( d[j].y * d[k].x + d[j].x * d[k].y ) );
	     p->elco[i].c26        += tmp;
	     jcell->elco[jnum].c26 += tmp;
	     c.c26                 += 2 * tmp;
	     tmp = 0.25 * ( ( ddphi_jk[k].zy + ddphi_jk[k].yz ) * d[j].z * d[k].z
			    + ddphi_jk[k].zz * ( d[j].z * d[k].y + d[j].y * d[k].z ) );
	     p->elco[i].c34        += tmp;
	     jcell->elco[jnum].c34 += tmp;
	     c.c34                 += 2 * tmp;
	     tmp = 0.25 * ( ddphi_jk[k].zz * ( d[j].z * d[k].x + d[j].x * d[k].z )
			    + ( ddphi_jk[k].zx + ddphi_jk[k].xz ) * d[j].z * d[k].z );
	     p->elco[i].c35        += tmp;
	     jcell->elco[jnum].c35 += tmp;
	     c.c35                 += 2 * tmp;
	     tmp = 0.25 * ( ddphi_jk[k].zx * d[j].z * d[k].y
			    + ddphi_jk[k].xz * d[j].y * d[k].z
			    + ddphi_jk[k].zy * d[j].z * d[k].x
			    + ddphi_jk[k].yz * d[j].x * d[k].z );
	     p->elco[i].c36        += tmp;
	     jcell->elco[jnum].c36 += tmp;
	     c.c36                 += 2 * tmp;

	     /* B', contribution from k */
	     dbulkm_dp -= 3.0 * db_k[k] 
	       * ( ddfc[j] * phi_a  + 2.0 * dfc[j] * dphi_a + fc[j] * ddphi_a )
	       * r[j] * r[j] * r[k];

		for (l=0; l<neigh->n; ++l) if (l!=j) 
		{
		  tmp = 0.5 * ddphi_lk I(l,k) .xx * d[l].x * d[k].x;
		  p->elco[i].c11        += tmp;
		  jcell->elco[jnum].c11 += tmp;
		  c.c11                 += 2 * tmp;
		  tmp = 0.5 * ddphi_lk I(l,k) .xy * d[l].x * d[k].y;
		  p->elco[i].c12        += tmp;
		  jcell->elco[jnum].c12 += tmp;
		  c.c12                 += 2 * tmp;
		  tmp = 0.5 * ddphi_lk I(l,k) .zx * d[l].z * d[k].x;
		  p->elco[i].c13        += tmp;
		  jcell->elco[jnum].c13 += tmp;
		  c.c13                 += 2 * tmp;
		  tmp = 0.5 * ddphi_lk I(l,k) .yy * d[l].y * d[k].y;
		  p->elco[i].c22        += tmp;
		  jcell->elco[jnum].c22 += tmp;
		  c.c22                 += 2 * tmp;
		  tmp = 0.5 * ddphi_lk I(l,k) .yz * d[l].y * d[k].z;
 		  p->elco[i].c23        += tmp;
		  jcell->elco[jnum].c23 += tmp;
		  c.c23                 += 2 * tmp;
		  tmp = 0.5 * ddphi_lk I(l,k) .zz * d[l].z * d[k].z;
 		  p->elco[i].c33        += tmp;
		  jcell->elco[jnum].c33 += tmp;
		  c.c33                 += 2 * tmp;
		  tmp = 0.125 * ( ddphi_lk I(l,k) .xx * d[l].y * d[k].y
		                + ddphi_lk I(l,k) .yy * d[l].x * d[k].x
		            + 2 * ddphi_lk I(l,k) .xy * d[l].y * d[k].x );
 		  p->elco[i].c44        += tmp;
		  jcell->elco[jnum].c44 += tmp;
		  c.c44                 += 2 * tmp;
		  tmp = 0.125 * ( ddphi_lk I(l,k) .yz * d[l].z * d[k].x
		                + ddphi_lk I(l,k) .zz * d[l].y * d[k].x 
		                + ddphi_lk I(l,k) .yx * d[l].z * d[k].z
		                + ddphi_lk I(l,k) .zx * d[l].y * d[k].z );
 		  p->elco[i].c45        += tmp;
		  jcell->elco[jnum].c45 += tmp;
		  c.c45                 += 2 * tmp;		  
		  tmp = 0.125 * ( ddphi_lk I(l,k) .xy * d[l].y * d[k].z
		                + ddphi_lk I(l,k) .yy * d[l].x * d[k].z 
		                + ddphi_lk I(l,k) .xz * d[l].y * d[k].y
		                + ddphi_lk I(l,k) .yz * d[l].x * d[k].y );
 		  p->elco[i].c46        += tmp;
		  jcell->elco[jnum].c46 += tmp;
		  c.c46                 += 2 * tmp;	
		  tmp = 0.125 * ( ddphi_lk I(l,k) .zz * d[l].x * d[k].x
		                + ddphi_lk I(l,k) .xx * d[l].z * d[k].z
		            + 2 * ddphi_lk I(l,k) .zx * d[l].x * d[k].z );
 		  p->elco[i].c55        += tmp;
		  jcell->elco[jnum].c55 += tmp;
		  c.c55                 += 2 * tmp;	
		  tmp = 0.125 * ( ddphi_lk I(l,k) .zx * d[l].x * d[k].y
		                + ddphi_lk I(l,k) .xx * d[l].z * d[k].y 
		                + ddphi_lk I(l,k) .zy * d[l].x * d[k].x
		                + ddphi_lk I(l,k) .xy * d[l].z * d[k].x );
 		  p->elco[i].c56        += tmp;
		  jcell->elco[jnum].c56 += tmp;
		  c.c56                 += 2 * tmp;	
		  tmp = 0.125 * ( ddphi_lk I(l,k) .xx * d[l].y * d[k].y
		                + ddphi_lk I(l,k) .yy * d[l].x * d[k].x
		            + 2 * ddphi_lk I(l,k) .xy * d[l].y * d[k].x );
 		  p->elco[i].c66        += tmp;
		  jcell->elco[jnum].c66 += tmp;
		  c.c66                 += 2 * tmp;			  
		  tmp = 0.25  * ( ddphi_lk I(l,k) .xy * d[l].x * d[k].z
		                 + ddphi_lk I(l,k) .zx * d[l].y * d[k].x );
 		  p->elco[i].c14        += tmp;
		  jcell->elco[jnum].c14 += tmp;
		  c.c14                 += 2 * tmp;	
		  tmp = 0.25 * ( ddphi_lk I(l,k) .zx * d[l].x * d[k].x
		               + ddphi_lk I(l,k) .xx * d[l].x * d[k].z );
 		  p->elco[i].c15        += tmp;
		  jcell->elco[jnum].c15 += tmp;
		  c.c15                 += 2 * tmp;	
		  tmp = 0.25  * ( ddphi_lk I(l,k) .xx * d[l].x * d[k].y
		                + ddphi_lk I(l,k) .xy * d[l].x * d[k].x );
 		  p->elco[i].c16        += tmp;
		  jcell->elco[jnum].c16 += tmp;
		  c.c16                 += 2 * tmp;	
		  tmp = 0.25  * ( ddphi_lk I(l,k) .yy * d[l].y * d[k].z 
		                + ddphi_lk I(l,k) .yz * d[l].y * d[k].y );
		  p->elco[i].c24        += tmp;
		  jcell->elco[jnum].c24 += tmp;
		  c.c24                 += 2 * tmp;	
		  tmp = 0.25  * ( ddphi_lk I(l,k) .yz * d[l].y * d[k].x
		                + ddphi_lk I(l,k) .xy * d[l].z * d[k].y );
		  p->elco[i].c25        += tmp;
		  jcell->elco[jnum].c25 += tmp;
		  c.c25                 += 2 * tmp;	
		  tmp = 0.25  * ( ddphi_lk I(l,k) .xy * d[l].y * d[k].y
		                + ddphi_lk I(l,k) .yy * d[l].y * d[k].x );
		  p->elco[i].c26        += tmp;
		  jcell->elco[jnum].c26 += tmp;
		  c.c26                 += 2 * tmp;	
		  tmp = 0.25  * ( ddphi_lk I(l,k) .zx * d[l].z * d[k].z
		                + ddphi_lk I(l,k) .zz * d[l].z * d[k].y );
		  p->elco[i].c34        += tmp;
		  jcell->elco[jnum].c34 += tmp;
		  c.c34                 += 2 * tmp;	
		  tmp = 0.25  * ( ddphi_lk I(l,k) .zz * d[l].z * d[k].x
		                + ddphi_lk I(l,k) .zx * d[l].z * d[k].z );
		  p->elco[i].c35        += tmp;
		  jcell->elco[jnum].c35 += tmp;
		  c.c35                 += 2 * tmp;	
		  tmp = 0.25  * ( ddphi_lk I(l,k) .zx * d[l].z * d[k].y
		                + ddphi_lk I(l,k) .yz * d[l].x * d[k].z );
		  p->elco[i].c36        += tmp;
		  jcell->elco[jnum].c36 += tmp;
		  c.c36                 += 2 * tmp;	

		  /* B', contribution from l */
		  dbulkm_dp -= 3.0 * ddb_lk I(l,k) 
		    * ( dfc[j] * phi_a + fc[j] * dphi_a  )
		    * r[j] * r[k] * r[l];
		  
		  /* B' contribution from m */
		  for (m=0; m<neigh->n; ++m) if (m!=j) 
		    dbulkm_dp -= fc[j] * dddb_lkm J(l,k,m) * phi_a 
		      * r[k] * r[l] * r[m];

		} /* l */

	    } /* k */
            
	  } /* neighbor j */

	} /* i */

      } /* loop over cells */ 

}
