
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
* imd_elco_stiweb
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "util.h"

/******************************************************************************
*
*  init_stiweb
*
******************************************************************************/

void init_stiweb(void) {

  int  i, j, k, n, m;
  real tmp;

  /* Parameters for more than one atom type */
  n = 0; m = 0;
  for (i=0; i<ntypes; i++) 
    for (j=i; j<ntypes; j++) {
      sw_a[i][j]  = sw_a[j][i]  = stiweb_a[n];
      sw_b[i][j]  = sw_b[j][i]  = stiweb_b[n];
      sw_p[i][j]  = sw_p[j][i]  = stiweb_p[n];
      sw_q[i][j]  = sw_q[j][i]  = stiweb_q[n];
      sw_a1[i][j] = sw_a1[j][i] = stiweb_a1[n];
      sw_de[i][j] = sw_de[j][i] = stiweb_de[n];
      sw_ga[i][j] = sw_ga[j][i] = stiweb_ga[n];
      sw_a2[i][j] = sw_a2[j][i] = stiweb_a2[n];
      n++;
      for (k=0; k<ntypes; k++) {
	sw_la[k][i][j] = sw_la[k][j][i] = stiweb_la[m];
	m++;
      }
    }

  for (i=0; i<ntypes; ++i)
    for (j=0; j<ntypes; ++j) 
      sw_2_a1[i][j] = sw_a1[i][j] * sw_a1[i][j];

  for (i=0; i<ntypes; ++i)
    for (j=0; j<ntypes; ++j) 
      sw_r_cut[i][j] = MAX(sw_a1[i][j], sw_a2[i][j] );
 
  if( r_cell != -1.0)
    r2_cut = SQR(r_cell);
  else { 
    tmp = 0.0;
    for (i=0; i<ntypes; ++i)
      for (j=0; j<ntypes; ++j) {
	tmp = MAX( tmp, sw_a1[i][j] );
	tmp = MAX( tmp, sw_a2[i][j] );
      }
    r2_cut = MAX(r2_cut,tmp*tmp);
  }
}

/******************************************************************************
*
*  do_elco_stiweb -- computes stress tensor and elastic constants using the 
*             Stillinger-Weber potential
*
******************************************************************************/

void do_elco_stiweb(void)
{
  cell   *p;
  int    ic, jc, kc;
  real   *r, *fc, *dfc;
  vektor *d;
  neightab *neigh;
  vektor dcos_j, dcos_k;
  tensor ddcos_jj, ddcos_jk, ddcos_kk;
  vektor dphi_j, dphi_k;
  tensor ddphi_jj, ddphi_jk, ddphi_kk;
  cell   *jcell, *kcell;
  int    i, j, k, p_typ, j_typ, k_typ, jnum, knum;
  real   *tmpptr;
  real   tmp_r, tmp_jj;
  real   inv_c, inv_r, tmp_frc1, tmp_frc2, tmp_frc;
  real   pot, f_cut;
  real   tmp_phi1, tmp_jk, tmp_kk, tmp_j2, tmp_k2;
  real   tmp_cos, tmp_l, tmp_j, tmp_k, cos_theta;
  real   phi_r, phi_a;
  real   tmp, tmp2_j, tmp2_k, tmp2_l;
  real   tmp1_b, tmp2_b, ddphi_jj_rr;
  real   tmp3_b, tmp_b_j, tmp_b_k, ddphi_jj_b, ddphi_kk_b, ddphi_jk_b;
  real   tmp_p, tmp_q, dddphi_jjj_rrr;
  real   tmp_aj, tmp_ak, dddphi_jjj, dddphi_jjk, dddphi_kkj, dddphi_kkk;

  d    = (vektor *) malloc( neigh_len * sizeof(vektor) );
  r    = (real *)   malloc( neigh_len * sizeof(real)   );
  fc   = (real *)   malloc( neigh_len * sizeof(real)   );
  dfc  = (real *)   malloc( neigh_len * sizeof(real)   );

  if ((d==NULL) || (r==NULL) || (fc==NULL) || (dfc==NULL))
    error("cannot allocate memory for temporary neighbor data");

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

	    /* cutoff function for three-body term */
	    tmp_r  = 1.0 / ( r[j] - sw_a2[p_typ][j_typ] );
	    fc[j]  = exp( sw_ga[p_typ][j_typ] * tmp_r );
	    dfc[j] = - sw_ga[p_typ][j_typ] * tmp_r * tmp_r;
	     
	  } /* j */

	  /*************************************************************/

	  /* For each neighbor of i */
	  for (j=0; j<neigh->n; ++j)
	  {

	    j_typ = neigh->typ[j];
	    jcell = (cell *) neigh->cl [j];
	    jnum  = neigh->num[j];
	    tmp_jj = 1 / ( r[j] * r[j] );

	    /* Potential energy of pair potential part */
	    phi_r  =   0.5 * sw_a[p_typ][j_typ] * pow( r[j], - sw_p[p_typ][j_typ] );
	    phi_a  = - 0.5 * sw_b[p_typ][j_typ] * pow( r[j], - sw_q[p_typ][j_typ] );
	    inv_c  = 1.0 / ( r[j] - sw_a1[p_typ][j_typ] );
	    inv_r  = 1.0 / r[j];
	    f_cut  = exp( sw_de[p_typ][j_typ] * inv_c );

	    pot  = ( phi_r + phi_a ) * f_cut;

	    epot += pot;

	    /* First derivative of pair potential part */
	    tmp_frc1 = pot * sw_de[p_typ][j_typ] * inv_c * inv_c;
	    tmp_frc2 = f_cut * inv_r * ( sw_p[p_typ][j_typ] * phi_r 
					 + sw_q[p_typ][j_typ] * phi_a );
	    tmp_frc = - inv_r * ( tmp_frc1 + tmp_frc2 );

	    dphi_j.x = tmp_frc * d[j].x;
	    dphi_j.y = tmp_frc * d[j].y;
	    dphi_j.z = tmp_frc * d[j].z;

	    /* Compute stress for pair potential part */
	    tmp = 0.5 * dphi_j.x * d[j].x;
	    p->stress[i].xx        += tmp;
	    jcell->stress[jnum].xx += tmp;
	    sigma.xx               += 2 * tmp;
	    tmp = 0.5 * dphi_j.y * d[j].y;
	    p->stress[i].yy        += tmp;
	    jcell->stress[jnum].yy += tmp;
	    sigma.yy               += 2 * tmp;
	    tmp = 0.5 * dphi_j.z * d[j].z;
	    p->stress[i].zz        += tmp;
	    jcell->stress[jnum].zz += tmp;
	    sigma.zz               += 2 * tmp;
	    tmp = 0.5 * dphi_j.y * d[j].z;
	    p->stress[i].yz        += tmp;
	    jcell->stress[jnum].yz += tmp;
	    sigma.yz               += 2 * tmp;
	    tmp = 0.5 * dphi_j.z * d[j].x;
	    p->stress[i].zx        += tmp;
	    jcell->stress[jnum].zx += tmp;
	    sigma.zx               += 2 * tmp;
	    tmp = 0.5 * dphi_j.x * d[j].y;
	    p->stress[i].xy        += tmp;
	    jcell->stress[jnum].xy += tmp;
	    sigma.xy               += 2 * tmp;

	    /* Second derivatives of pair potential part */
	    tmp_phi1 = ( ( tmp_frc1 + 2 * tmp_frc2 ) 
	      * ( sw_de[p_typ][j_typ] * inv_c * inv_c + inv_r ) 
	        + 2 * tmp_frc1 * inv_c + f_cut * inv_r * inv_r 
	      * ( sw_p[p_typ][j_typ] * sw_p[p_typ][j_typ] * phi_r 
		+ sw_q[p_typ][j_typ] * sw_q[p_typ][j_typ] * phi_a ) ) * inv_r * inv_r;

	    ddphi_jj.xx = tmp_phi1 * d[j].x * d[j].x + tmp_frc;
	    ddphi_jj.xy = tmp_phi1 * d[j].x * d[j].y;
	    ddphi_jj.yy = tmp_phi1 * d[j].y * d[j].y + tmp_frc;
	    ddphi_jj.yz = tmp_phi1 * d[j].y * d[j].z;
	    ddphi_jj.zx = tmp_phi1 * d[j].z * d[j].x;
	    ddphi_jj.zz = tmp_phi1 * d[j].z * d[j].z + tmp_frc;

	    /* For computation of Bulk modulus */
	    tmp1_b = sw_de[p_typ][j_typ] * r[j] * inv_c * inv_c;
	    tmp2_b = 2.0 * tmp1_b * r[j] * inv_c;
	    tmp_p = sw_p[p_typ][j_typ];
	    tmp_q = sw_q[p_typ][j_typ];

	    ddphi_jj_rr = ( tmp2_b + tmp_p 
	      + ( tmp1_b + tmp_p ) * ( tmp1_b + tmp_p ) ) * phi_r * f_cut
	     + ( tmp2_b + tmp_q 
	      + ( tmp1_b + tmp_q ) * ( tmp1_b + tmp_q ) ) * phi_a * f_cut;

	    bulkm += ddphi_jj_rr;

	    /* For computation of B' */
	    tmp3_b = 3.0 * tmp2_b * r[j] * inv_c;
	    
	    dddphi_jjj_rrr = - ( tmp3_b + 2.0 * tmp_p 
	      + 3.0 * ( tmp1_b + tmp_p ) * ( tmp2_b + tmp_p )
	      + ( tmp1_b + tmp_p ) * ( tmp1_b + tmp_p ) * ( tmp1_b + tmp_p ) ) * phi_r * f_cut
	     - ( tmp3_b + 2.0 * tmp_q 
	      + 3.0 * ( tmp1_b + tmp_q ) * ( tmp2_b + tmp_q )
	      + ( tmp1_b + tmp_q ) * ( tmp1_b + tmp_q ) * ( tmp1_b + tmp_q ) ) * phi_a * f_cut;

	    dbulkm_dp += dddphi_jjj_rrr;

	    /* Compute elastic constants for pair potential part */
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
	    tmp = 0.125 * ( ddphi_jj.yy * d[j].z * d[j].z 
	                  + ddphi_jj.zz * d[j].y * d[j].y
	              + 2 * ddphi_jj.yz * d[j].y * d[j].z );
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

	    for (k=j+1; k<neigh->n; ++k) {

	      k_typ = neigh->typ[k];
	      kcell = (cell *) neigh->cl [k];
	      knum  = neigh->num[k];

	      tmp_jk    = 1 / ( r[j] * r[k] );  
	      cos_theta = SPROD(d[j],d[k]) * tmp_jk;
	      tmp_cos = cos_theta + 1.0 / 3.0;
	      tmp_j2 = cos_theta * tmp_jj;
	      tmp_kk = 1 /  ( r[k] * r[k] );
	      tmp_k2 = cos_theta * tmp_kk;

	      /* Potential energy of three-body part */
	      epot += sw_la[p_typ][j_typ][k_typ] * fc[j] * fc[k] * tmp_cos * tmp_cos;

	      /* First derivatives of cos(theta) */
	      dcos_j.x = tmp_jk * d[k].x - tmp_j2 * d[j].x;
	      dcos_j.y = tmp_jk * d[k].y - tmp_j2 * d[j].y;
	      dcos_j.z = tmp_jk * d[k].z - tmp_j2 * d[j].z;
	      
	      dcos_k.x = tmp_jk * d[j].x - tmp_k2 * d[k].x;
	      dcos_k.y = tmp_jk * d[j].y - tmp_k2 * d[k].y;
	      dcos_k.z = tmp_jk * d[j].z - tmp_k2 * d[k].z;
	
	      /* First derivatives of three-body potential part */	      
	      tmp_l = sw_la[p_typ][j_typ][k_typ] * fc[j] * fc[k] * tmp_cos;
	      tmp_j = tmp_cos * dfc[j] / r[j];
	      tmp_k = tmp_cos * dfc[k] / r[k];

	      dphi_j.x = tmp_l * ( tmp_j * d[j].x + 2.0 * dcos_j.x );
	      dphi_j.y = tmp_l * ( tmp_j * d[j].y + 2.0 * dcos_j.y );
	      dphi_j.z = tmp_l * ( tmp_j * d[j].z + 2.0 * dcos_j.z );

	      dphi_k.x = tmp_l * ( tmp_k * d[k].x + 2.0 * dcos_k.x );
	      dphi_k.y = tmp_l * ( tmp_k * d[k].y + 2.0 * dcos_k.y );
	      dphi_k.z = tmp_l * ( tmp_k * d[k].z + 2.0 * dcos_k.z );     

	      /* Compute stress for three-body potential part */	      
	      tmp = 0.5 * ( dphi_j.x * d[j].x + dphi_k.x * d[k].x );
	      p->stress[i].xx        += tmp;
	      jcell->stress[jnum].xx += 0.5 * tmp;
	      kcell->stress[knum].xx += 0.5 * tmp;
	      sigma.xx               += 2 * tmp;
	      tmp = 0.5 * ( dphi_j.y * d[j].y + dphi_k.y * d[k].y );
	      p->stress[i].yy        += tmp;
	      jcell->stress[jnum].yy += 0.5 * tmp;
	      kcell->stress[knum].yy += 0.5 * tmp;
	      sigma.yy               += 2 * tmp;
	      tmp = 0.5 * ( dphi_j.z * d[j].z + dphi_k.z * d[k].z );
	      p->stress[i].zz        += tmp;
	      jcell->stress[jnum].zz += 0.5 * tmp;
	      kcell->stress[knum].zz += 0.5 * tmp;
	      sigma.zz               += 2 * tmp;
	      tmp = 0.25 * ( dphi_j.y * d[j].z + dphi_k.y * d[k].z
			   + dphi_j.z * d[j].y + dphi_k.z * d[k].y );
	      p->stress[i].yz        += tmp;
	      jcell->stress[jnum].yz += 0.5 * tmp;
	      kcell->stress[knum].yz += 0.5 * tmp;
	      sigma.yz               += 2 * tmp;
	      tmp = 0.25 * ( dphi_j.z * d[j].x + dphi_k.z * d[k].x
			   + dphi_j.x * d[j].z + dphi_k.x * d[k].z );
	      p->stress[i].zx        += tmp;
	      jcell->stress[jnum].zx += 0.5 * tmp;
	      kcell->stress[knum].zx += 0.5 * tmp;
	      sigma.zx               += 2 * tmp;
	      tmp = 0.25 * ( dphi_j.x * d[j].y + dphi_k.x * d[k].y
			   + dphi_j.y * d[j].x + dphi_k.y * d[k].x );
	      p->stress[i].xy        += tmp;
	      jcell->stress[jnum].xy += 0.5 * tmp;
	      kcell->stress[knum].xy += 0.5 * tmp;
	      sigma.xy               += 2 * tmp;
	      
	      /* Second derivatives of cos(theta) */ 
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

	      /* Second derivatives of three-body potential part */
	      tmp_j = dfc[j] / r[j];
	      tmp2_j = - ( 3 * r[j] - sw_a2[p_typ][j_typ] ) * tmp_jj / ( r[j] - sw_a2[p_typ][j_typ] );
	      tmp_k = dfc[k] / r[k];
	      tmp2_k = - ( 3 * r[k] - sw_a2[p_typ][k_typ] ) * tmp_kk / ( r[k] - sw_a2[p_typ][k_typ] );
	      tmp2_l = 2 * sw_la[p_typ][j_typ][k_typ] * fc[j] * fc[k];

	      ddphi_jj.xx = tmp_j * ( dphi_j.x * d[j].x   
	        + tmp_l * ( tmp2_j * d[j].x * d[j].x + 4 * d[j].x * dcos_j.x + 1.0 ) )
	        + 2 * tmp_l * ddcos_jj.xx
	       	+ tmp2_l * dcos_j.x * dcos_j.x;
 
	      ddphi_jj.xy = tmp_j * ( dphi_j.y * d[j].x  
	        + tmp_l * ( tmp2_j * d[j].x * d[j].y + 2 * ( d[j].x * dcos_j.y + d[j].y * dcos_j.x ) ) )
	        + 2 * tmp_l * ddcos_jj.xy
	       	+ tmp2_l * dcos_j.x * dcos_j.y;

	      ddphi_jj.yy = tmp_j * ( dphi_j.y * d[j].y   
	        + tmp_l * ( tmp2_j * d[j].y * d[j].y + 4 * d[j].y * dcos_j.y + 1.0 ) )
	        + 2 * tmp_l * ddcos_jj.yy
	       	+ tmp2_l * dcos_j.y * dcos_j.y;

	      ddphi_jj.yz = tmp_j * ( dphi_j.z * d[j].y  
	        + tmp_l * ( tmp2_j * d[j].y * d[j].z + 2 * ( d[j].y * dcos_j.z + d[j].z * dcos_j.y ) ) )
	        + 2 * tmp_l * ddcos_jj.yz
	       	+ tmp2_l * dcos_j.y * dcos_j.z; 

	      ddphi_jj.zx = tmp_j * ( dphi_j.x * d[j].z  
	        + tmp_l * ( tmp2_j * d[j].z * d[j].x + 2 * ( d[j].z * dcos_j.x + d[j].x * dcos_j.z ) ) )
	        + 2 * tmp_l * ddcos_jj.zx
	       	+ tmp2_l * dcos_j.z * dcos_j.x; 

	      ddphi_jj.zz = tmp_j * ( dphi_j.z * d[j].z   
	        + tmp_l * ( tmp2_j * d[j].z * d[j].z + 4 * d[j].z * dcos_j.z + 1.0 ) )
	        + 2 * tmp_l * ddcos_jj.zz
	       	+ tmp2_l * dcos_j.z * dcos_j.z; 

	      ddphi_jk.xx = tmp_j * d[j].x * dphi_k.x + 2 * tmp_l * ( tmp_k * dcos_j.x * d[k].x + ddcos_jk.xx )
		+ tmp2_l * dcos_j.x * dcos_k.x;

	      ddphi_jk.xy = tmp_j * d[j].x * dphi_k.y + 2 * tmp_l * ( tmp_k * dcos_j.x * d[k].y + ddcos_jk.xy )
		+ tmp2_l * dcos_j.x * dcos_k.y;

	      ddphi_jk.xz = tmp_j * d[j].x * dphi_k.z + 2 * tmp_l * ( tmp_k * dcos_j.x * d[k].z + ddcos_jk.xz )
		+ tmp2_l * dcos_j.x * dcos_k.z;

	      ddphi_jk.yx = tmp_j * d[j].y * dphi_k.x + 2 * tmp_l * ( tmp_k * dcos_j.y * d[k].x + ddcos_jk.yx )
		+ tmp2_l * dcos_j.y * dcos_k.x;

	      ddphi_jk.yy = tmp_j * d[j].y * dphi_k.y + 2 * tmp_l * ( tmp_k * dcos_j.y * d[k].y + ddcos_jk.yy )
		+ tmp2_l * dcos_j.y * dcos_k.y;

	      ddphi_jk.yz = tmp_j * d[j].y * dphi_k.z + 2 * tmp_l * ( tmp_k * dcos_j.y * d[k].z + ddcos_jk.yz )
		+ tmp2_l * dcos_j.y * dcos_k.z;

	      ddphi_jk.zx = tmp_j * d[j].z * dphi_k.x + 2 * tmp_l * ( tmp_k * dcos_j.z * d[k].x + ddcos_jk.zx )
		+ tmp2_l * dcos_j.z * dcos_k.x;

	      ddphi_jk.zy = tmp_j * d[j].z * dphi_k.y + 2 * tmp_l * ( tmp_k * dcos_j.z * d[k].y + ddcos_jk.zy )
		+ tmp2_l * dcos_j.z * dcos_k.y;

	      ddphi_jk.zz = tmp_j * d[j].z * dphi_k.z + 2 * tmp_l * ( tmp_k * dcos_j.z * d[k].z + ddcos_jk.zz )
		+ tmp2_l * dcos_j.z * dcos_k.z;

	      ddphi_kk.xx = tmp_k * ( d[k].x * dphi_k.x   
	        + tmp_l * ( tmp2_k * d[k].x * d[k].x + 4 * d[k].x * dcos_k.x + 1.0 ) )
	        + 2 * tmp_l * ddcos_kk.xx
	       	+ tmp2_l * dcos_k.x * dcos_k.x;

	      ddphi_kk.xy = tmp_k * ( d[k].x * dphi_k.y   
	        + tmp_l * ( tmp2_k * d[k].x * d[k].y + 2 * ( d[k].x * dcos_k.y + d[k].y * dcos_k.x ) ) )
	        + 2 * tmp_l * ddcos_kk.xy
	       	+ tmp2_l * dcos_k.x * dcos_k.y;

	      ddphi_kk.yy = tmp_k * ( d[k].y * dphi_k.y   
	        + tmp_l * ( tmp2_k * d[k].y * d[k].y + 4 * d[k].y * dcos_k.y + 1.0 ) )
	        + 2 * tmp_l * ddcos_kk.yy
	       	+ tmp2_l * dcos_k.y * dcos_k.y; 

	      ddphi_kk.yz = tmp_k * ( d[k].y * dphi_k.z   
	        + tmp_l * ( tmp2_k * d[k].y * d[k].z + 2 * ( d[k].y * dcos_k.z + d[k].z * dcos_k.y ) ) )
	        + 2 * tmp_l * ddcos_kk.yz
	       	+ tmp2_l * dcos_k.y * dcos_k.z;

	      ddphi_kk.zx = tmp_k * ( d[k].z * dphi_k.x   
	        + tmp_l * ( tmp2_k * d[k].z * d[k].x + 2 * ( d[k].z * dcos_k.x + d[k].x * dcos_k.z ) ) )
	        + 2 * tmp_l * ddcos_kk.zx
	       	+ tmp2_l * dcos_k.z * dcos_k.x; 

	      ddphi_kk.zz = tmp_k * ( d[k].z * dphi_k.z   
	        + tmp_l * ( tmp2_k * d[k].z * d[k].z + 4 * d[k].z * dcos_k.z + 1.0 ) )
	        + 2 * tmp_l * ddcos_kk.zz
	       	+ tmp2_l * dcos_k.z * dcos_k.z; 

	      /* For computation of bulk modulus */
	      tmp3_b = sw_la[p_typ][j_typ][k_typ] * fc[j] * fc[k] * tmp_cos * tmp_cos;
	      tmp_b_j = sw_ga[p_typ][j_typ] / ( r[j] - sw_a2[p_typ][j_typ] );
	      tmp_b_k = sw_ga[p_typ][k_typ] / ( r[k] - sw_a2[p_typ][k_typ] );
	      tmp_aj = sw_a2[p_typ][j_typ];
	      tmp_ak = sw_a2[p_typ][k_typ];

	      ddphi_jj_b = tmp_b_j / ( ( r[j] - tmp_aj ) * ( r[j] - tmp_aj ) )
		* ( 2.0 + tmp_b_j ) * tmp3_b;
	      ddphi_kk_b = tmp_b_k / ( r[k] - tmp_ak ) / ( r[k] - tmp_ak )
		* ( 2.0 + tmp_b_k ) * tmp3_b;
	      ddphi_jk_b = tmp_b_j / ( r[j] - tmp_aj )
		* tmp_b_k / ( r[k] - tmp_ak ) * tmp3_b;

	      bulkm += ddphi_jj_b * r[j] * r[j] 
		+ 2 * ddphi_jk_b * r[j] * r[k] 
		+ ddphi_kk_b * r[k] * r[k];

	      /* For computation of B' */	      
	      dddphi_jjj = - tmp_b_j / ( ( r[j] - tmp_aj ) * ( r[j] - tmp_aj ) * ( r[j] - tmp_aj ) )
		* ( 6.0 + tmp_b_j * ( 6.0 + tmp_b_j ) ) * tmp3_b;

	      dddphi_kkk = - tmp_b_k / ( ( r[k] - tmp_ak ) * ( r[k] - tmp_ak ) * ( r[k] - tmp_ak ) )
		* ( 6.0 + tmp_b_k * ( 6.0 + tmp_b_k ) ) * tmp3_b;

	      dddphi_jjk = - tmp_b_j * tmp_b_k / ( ( r[j] - tmp_aj ) * ( r[j] - tmp_aj ) * ( r[k] - tmp_ak ) )
		* ( tmp_b_j + 2.0 ) * tmp3_b;

	      dddphi_kkj = - tmp_b_j * tmp_b_k / ( ( r[j] - tmp_aj ) * ( r[k] - tmp_ak ) * ( r[k] - tmp_ak ) )
		* ( tmp_b_k + 2.0 ) * tmp3_b;

	      dbulkm_dp += dddphi_jjj * r[j] * r[j] * r[j]
                + 3.0 * dddphi_jjk * r[j] * r[j] * r[k]
                + 3.0 * dddphi_kkj * r[k] * r[k] * r[j]
		+  dddphi_kkk * r[k] * r[k] * r[k];

	      /* Compute elastic constants of the three-body potential part */
	      tmp = 0.5 * 
		      ( ddphi_jj.xx * d[j].x * d[j].x 
		  + 2 * ddphi_jk.xx * d[j].x * d[k].x 
		      + ddphi_kk.xx * d[k].x * d[k].x );
	      p->elco[i].c11        += tmp;
	      jcell->elco[jnum].c11 += 0.5 * tmp;
	      kcell->elco[knum].c11 += 0.5 * tmp;	      
	      c.c11                 += 2 * tmp;
	      tmp = 0.5 * 
		      ( ddphi_jj.xy * d[j].x * d[j].y 
		      + ddphi_jk.xy * d[j].x * d[k].y
		      + ddphi_jk.yx * d[k].x * d[j].y
		      + ddphi_kk.xy * d[k].x * d[k].y ); 
	      p->elco[i].c12        += tmp;
	      jcell->elco[jnum].c12 += 0.5 * tmp;
	      kcell->elco[knum].c12 += 0.5 * tmp;	      
	      c.c12                 += 2 * tmp;
	      tmp = 0.5 * 
		      ( ddphi_jj.zx * d[j].x * d[j].z 
		      + ddphi_jk.xz * d[j].x * d[k].z
		      + ddphi_jk.zx * d[k].x * d[j].z
		      + ddphi_kk.zx * d[k].x * d[k].z );   
	      p->elco[i].c13        += tmp;
	      jcell->elco[jnum].c13 += 0.5 * tmp;
	      kcell->elco[knum].c13 += 0.5 * tmp;	      
	      c.c13                 += 2 * tmp;
	      tmp = 0.5 * 
		      ( ddphi_jj.yy * d[j].y * d[j].y 
		  + 2 * ddphi_jk.yy * d[j].y * d[k].y 
		      + ddphi_kk.yy * d[k].y * d[k].y );  
	      p->elco[i].c22        += tmp;
	      jcell->elco[jnum].c22 += 0.5 * tmp;
	      kcell->elco[knum].c22 += 0.5 * tmp;	      
	      c.c22                 += 2 * tmp;
	      tmp = 0.5 * 
		      ( ddphi_jj.yz * d[j].y * d[j].z 
		      + ddphi_jk.yz * d[j].y * d[k].z
		      + ddphi_jk.zy * d[k].y * d[j].z
		      + ddphi_kk.yz * d[k].y * d[k].z );
	      p->elco[i].c23        += tmp;
	      jcell->elco[jnum].c23 += 0.5 * tmp;
	      kcell->elco[knum].c23 += 0.5 * tmp;	      
	      c.c23                 += 2 * tmp;
	      tmp = 0.5 * 
		      ( ddphi_jj.zz * d[j].z * d[j].z 
		  + 2 * ddphi_jk.zz * d[j].z * d[k].z 
		      + ddphi_kk.zz * d[k].z * d[k].z ); 
	      p->elco[i].c33        += tmp;
	      jcell->elco[jnum].c33 += 0.5 * tmp;
	      kcell->elco[knum].c33 += 0.5 * tmp;	      
	      c.c33                 += 2 * tmp;
	      tmp = 0.125 * 
		      ( ddphi_jj.yy * d[j].z * d[j].z 
		  + 2 * ddphi_jk.yy * d[j].z * d[k].z 
		      + ddphi_kk.yy * d[k].z * d[k].z 
		  + 2 * ddphi_jj.yz * d[j].y * d[j].z 
		   + 2 * ddphi_jk.yz * d[j].z * d[k].y + 2 * ddphi_jk.zy * d[j].y * d[k].z 
		   + 2 * ddphi_kk.yz * d[k].y * d[k].z 
		  + ddphi_jj.zz * d[j].y * d[j].y 
		   + 2 * ddphi_jk.zz * d[j].y * d[k].y 
		   + ddphi_kk.zz * d[k].y * d[k].y ); 
	      p->elco[i].c44        += tmp;
	      jcell->elco[jnum].c44 += 0.5 * tmp;
	      kcell->elco[knum].c44 += 0.5 * tmp;	      
	      c.c44                 += 2 * tmp;
	      tmp = 0.125 *
		( ddphi_jj.yz * d[j].z * d[j].x 
		  + ddphi_jk.yz * d[j].z * d[k].x
		   + ddphi_jk.zy * d[j].x * d[k].z
		   + ddphi_kk.yz * d[k].z * d[k].x 
		  + ddphi_jj.zz * d[j].y * d[j].x 
		  + ddphi_jk.zz * ( d[j].y * d[k].x + d[j].x * d[k].y ) 
		   + ddphi_kk.zz * d[k].y * d[k].x
		  + ddphi_jj.xy * d[j].z * d[j].z 
		  + ( ddphi_jk.yx + ddphi_jk.xy ) * d[j].z * d[k].z 
		   + ddphi_kk.xy * d[k].z * d[k].z
		  + ddphi_jj.zx * d[j].y * d[j].z 
		   + ddphi_jk.zx * d[j].y * d[k].z
		   + ddphi_jk.xz * d[j].z * d[k].y 
		   + ddphi_kk.zx * d[k].y * d[k].z ); 
	      p->elco[i].c45        += tmp;
	      jcell->elco[jnum].c45 += 0.5 * tmp;
	      kcell->elco[knum].c45 += 0.5 * tmp;	      
	      c.c45                 += 2 * tmp;              
	      tmp = 0.125 *
		( ddphi_jj.xy * d[j].y * d[j].z 
		  + ddphi_jk.xy * d[j].y * d[k].z
		   + ddphi_jk.yx * d[j].z * d[k].y
		   + ddphi_kk.xy * d[k].y * d[k].z 
		  + ddphi_jj.yy * d[j].x * d[j].z 
		  + ddphi_jk.yy * ( d[j].x * d[k].z + d[j].z * d[k].x ) 
		   + ddphi_kk.yy * d[k].x * d[k].z
		  + ddphi_jj.zx * d[j].y * d[j].y 
		  + ( ddphi_jk.xz + ddphi_jk.zx ) * d[j].y * d[k].y 
		   + ddphi_kk.zx * d[k].y * d[k].y
		  + ddphi_jj.yz * d[j].x * d[j].y 
		   + ddphi_jk.yz * d[j].x * d[k].y
		   + ddphi_jk.zy * d[j].y * d[k].x 
		   + ddphi_kk.yz * d[k].x * d[k].y ); 
	      p->elco[i].c46        += tmp;
	      jcell->elco[jnum].c46 += 0.5 * tmp;
	      kcell->elco[knum].c46 += 0.5 * tmp;	      
	      c.c46                 += 2 * tmp;              
	      tmp = 0.125 * 
		      ( ddphi_jj.zz * d[j].x * d[j].x 
		  + 2 * ddphi_jk.zz * d[j].x * d[k].x 
		      + ddphi_kk.zz * d[k].x * d[k].x 
		  + 2 * ddphi_jj.zx * d[j].z * d[j].x 
		   + 2 * ddphi_jk.zx * d[j].x * d[k].z + 2 * ddphi_jk.xz * d[j].z * d[k].x 
		   + 2 * ddphi_kk.zx * d[k].z * d[k].x 
		  + ddphi_jj.xx * d[j].z * d[j].z 
		   + 2 * ddphi_jk.xx * d[j].z * d[k].z 
		   + ddphi_kk.xx * d[k].z * d[k].z ); 
	      p->elco[i].c55        += tmp;
	      jcell->elco[jnum].c55 += 0.5 * tmp;
	      kcell->elco[knum].c55 += 0.5 * tmp;	      
	      c.c55                 += 2 * tmp;
	      tmp = 0.125 *
		( ddphi_jj.zx * d[j].x * d[j].y 
		  + ddphi_jk.zx * d[j].x * d[k].y
		   + ddphi_jk.xz * d[j].y * d[k].x
		   + ddphi_kk.zx * d[k].x * d[k].y 
		  + ddphi_jj.xx * d[j].z * d[j].y 
		  + ddphi_jk.xx * ( d[j].z * d[k].y + d[j].y * d[k].z ) 
		   + ddphi_kk.xx * d[k].z * d[k].y
		  + ddphi_jj.yz * d[j].x * d[j].x 
		  + ( ddphi_jk.zy + ddphi_jk.yz ) * d[j].x * d[k].x 
		   + ddphi_kk.yz * d[k].x * d[k].x
		  + ddphi_jj.xy * d[j].z * d[j].x 
		   + ddphi_jk.xy * d[j].z * d[k].x
		   + ddphi_jk.yx * d[j].x * d[k].z 
		   + ddphi_kk.xy * d[k].z * d[k].x ); 
	      p->elco[i].c56        += tmp;
	      jcell->elco[jnum].c56 += 0.5 * tmp;
	      kcell->elco[knum].c56 += 0.5 * tmp;	      
	      c.c56                 += 2 * tmp;              
	      tmp = 0.125 * 
		      ( ddphi_jj.xx * d[j].y * d[j].y 
		  + 2 * ddphi_jk.xx * d[j].y * d[k].y 
		      + ddphi_kk.xx * d[k].y * d[k].y 
		  + 2 * ddphi_jj.xy * d[j].x * d[j].y 
		   + 2 * ddphi_jk.xy * d[j].y * d[k].x + 2 * ddphi_jk.yx * d[j].x * d[k].y 
		   + 2 * ddphi_kk.xy * d[k].x * d[k].y 
		  + ddphi_jj.yy * d[j].x * d[j].x 
		   + 2 * ddphi_jk.yy * d[j].x * d[k].x 
		   + ddphi_kk.yy * d[k].x * d[k].x ); 
	      p->elco[i].c66        += tmp;
	      jcell->elco[jnum].c66 += 0.5 * tmp;
	      kcell->elco[knum].c66 += 0.5 * tmp;	      
	      c.c66                 += 2 * tmp;
	      tmp = 0.25 * 
		 ( ddphi_jj.xy * d[j].x * d[j].z 
		  + ddphi_jk.xy * d[j].x * d[k].z
 		  + ddphi_jk.yx * d[j].z * d[k].x
		  + ddphi_kk.xy * d[k].x * d[k].z 
		 + ddphi_jj.zx * d[j].x * d[j].y 
		  + ddphi_jk.xz * d[j].x * d[k].y
		  + ddphi_jk.zx * d[j].y * d[k].x 
		  + ddphi_kk.zx * d[k].x * d[k].y ); 
	      p->elco[i].c14        += tmp;
	      jcell->elco[jnum].c14 += 0.5 * tmp;
	      kcell->elco[knum].c14 += 0.5 * tmp;	      
	      c.c14                 += 2 * tmp;
	      tmp = 0.25 * 
		 ( ddphi_jj.zx * d[j].x * d[j].x 
		  + ( ddphi_jk.xz + ddphi_jk.zx ) * d[j].x * d[k].x
		  + ddphi_kk.zx * d[k].x * d[k].x 
		 + ddphi_jj.xx * d[j].x * d[j].z 
	          + ddphi_jk.xx * ( d[j].x * d[k].z + d[j].z * d[k].x )
		  + ddphi_kk.xx * d[k].x * d[k].z ); 
	      p->elco[i].c15        += tmp;
	      jcell->elco[jnum].c15 += 0.5 * tmp;
	      kcell->elco[knum].c15 += 0.5 * tmp;	      
	      c.c15                 += 2 * tmp;
	      tmp = 0.25 * 
		 ( ddphi_jj.xx * d[j].x * d[j].y 
	          + ddphi_jk.xx * ( d[j].x * d[k].y + d[j].y * d[k].x )
		  + ddphi_kk.xx * d[k].x * d[k].y 
		 + ddphi_jj.xy * d[j].x * d[j].x 
	          + ( ddphi_jk.xy + ddphi_jk.yx ) * d[j].x * d[k].x 
		  + ddphi_kk.xy * d[k].x * d[k].x ); 
	      p->elco[i].c16        += tmp;
	      jcell->elco[jnum].c16 += 0.5 * tmp;
	      kcell->elco[knum].c16 += 0.5 * tmp;	      
	      c.c16                 += 2 * tmp;
	      tmp = 0.25 * 
		 ( ddphi_jj.yy * d[j].y * d[j].z 
	          + ddphi_jk.yy * ( d[j].y * d[k].z + d[j].z * d[k].y )
		  + ddphi_kk.yy * d[k].y * d[k].z 
		 + ddphi_jj.yz * d[j].y * d[j].y 
	          + ( ddphi_jk.yz + ddphi_jk.zy ) * d[j].y * d[k].y 
		  + ddphi_kk.yz * d[k].y * d[k].y ); 
	      p->elco[i].c24        += tmp;
	      jcell->elco[jnum].c24 += 0.5 * tmp;
	      kcell->elco[knum].c24 += 0.5 * tmp;	      
	      c.c24                 += 2 * tmp;
	      tmp = 0.25 * 
		 ( ddphi_jj.yz * d[j].y * d[j].x 
		  + ddphi_jk.yz * d[j].y * d[k].x
 		  + ddphi_jk.zy * d[j].x * d[k].y
		  + ddphi_kk.yz * d[k].y * d[k].x 
		 + ddphi_jj.xy * d[j].y * d[j].z 
		  + ddphi_jk.yx * d[j].y * d[k].z
		  + ddphi_jk.xy * d[j].z * d[k].y 
		  + ddphi_kk.xy * d[k].y * d[k].z ); 
	      p->elco[i].c25        += tmp;
	      jcell->elco[jnum].c25 += 0.5 * tmp;
	      kcell->elco[knum].c25 += 0.5 * tmp;	      
	      c.c25                 += 2 * tmp;
	      tmp = 0.25 * 
		 ( ddphi_jj.xy * d[j].y * d[j].y 
		  + ( ddphi_jk.yx + ddphi_jk.xy ) * d[j].y * d[k].y
		  + ddphi_kk.xy * d[k].y * d[k].y
		 + ddphi_jj.yy * d[j].y * d[j].x 
	          + ddphi_jk.yy * ( d[j].y * d[k].x + d[j].x * d[k].y )
		  + ddphi_kk.yy * d[k].y * d[k].x ); 
	      p->elco[i].c26        += tmp;
	      jcell->elco[jnum].c26 += 0.5 * tmp;
	      kcell->elco[knum].c26 += 0.5 * tmp;	      
	      c.c26                 += 2 * tmp;
	      tmp = 0.25 * 
		 ( ddphi_jj.yz * d[j].z * d[j].z 
		  + ( ddphi_jk.zy + ddphi_jk.yz ) * d[j].z * d[k].z
		  + ddphi_kk.yz * d[k].z * d[k].z
		 + ddphi_jj.zz * d[j].z * d[j].y 
	          + ddphi_jk.zz * ( d[j].z * d[k].y + d[j].y * d[k].z )
		  + ddphi_kk.zz * d[k].z * d[k].y ); 
	      p->elco[i].c34        += tmp;
	      jcell->elco[jnum].c34 += 0.5 * tmp;
	      kcell->elco[knum].c34 += 0.5 * tmp;
	      c.c34                 += 2 * tmp;
	      tmp = 0.25 * 
		 ( ddphi_jj.zz * d[j].z * d[j].x 
	          + ddphi_jk.zz * ( d[j].z * d[k].x + d[j].x * d[k].z )
		  + ddphi_kk.zz * d[k].z * d[k].x 
		 + ddphi_jj.zx * d[j].z * d[j].z 
	          + ( ddphi_jk.zx + ddphi_jk.xz ) * d[j].z * d[k].z 
		  + ddphi_kk.zx * d[k].z * d[k].z ); 
	      p->elco[i].c35        += tmp;
	      jcell->elco[jnum].c35 += 0.5 * tmp;
	      kcell->elco[knum].c35 += 0.5 * tmp;	      
	      c.c35                 += 2 * tmp;
	      tmp = 0.25 * 
		 ( ddphi_jj.zx * d[j].z * d[j].y 
		  + ddphi_jk.zx * d[j].z * d[k].y
 		  + ddphi_jk.xz * d[j].y * d[k].z
		  + ddphi_kk.zx * d[k].z * d[k].y 
		 + ddphi_jj.yz * d[j].z * d[j].x 
		  + ddphi_jk.zy * d[j].z * d[k].x
		  + ddphi_jk.yz * d[j].x * d[k].z 
		  + ddphi_kk.yz * d[k].z * d[k].x ); 
	      p->elco[i].c36        += tmp;
	      jcell->elco[jnum].c36 += 0.5 * tmp;
	      kcell->elco[knum].c36 += 0.5 * tmp;	      
	      c.c36                 += 2 * tmp;

	    } /* k */

	  } /* j */

	}

      }

}
