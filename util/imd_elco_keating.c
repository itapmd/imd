
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
* imd_elco_keating
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "util.h"

/******************************************************************************
*
*  init_keating
*
******************************************************************************/

void init_keating(void) {

  int  i, j, k, n, m;
  real tmp;

  /* Parameters for more than one atom type */
  n = 0; m = 0;
  for (i=0; i<ntypes; i++) 
    for (j=i; j<ntypes; j++) {
      keat_alpha[i][j]  = keat_alpha[j][i]  = keating_alpha[n];
      keat_d[i][j]      = keat_d[j][i]      = keating_d[n];
      keat_r_cut[i][j]  = keat_r_cut[j][i]  = keating_r_cut[n];
      keat_r2_cut[i][j] = keat_r2_cut[j][i] = SQR( keat_r_cut[i][j] );
      n++; 
      for (k=0; k<ntypes; k++) {
	keat_beta[k][i][j] = keat_beta[k][j][i] = keating_beta[m];
	m++;
      }
    }

  tmp = 0.0;
  for (i=0; i<ntypes; ++i)
    for (j=0; j<ntypes; ++j)
      tmp = MAX( tmp, keat_r2_cut[i][j] );
  r2_cut = MAX(r2_cut, tmp*tmp);

}	  

/******************************************************************************
*
*  do_elco_keating -- computes stress tensor and elastic constants using the 
*             Keating potential
*
******************************************************************************/

void do_elco_keating(void)
{

  int    ic, jc, kc, i, j, k;
  cell   *p, *jcell, *kcell;
  int    jnum, knum;
  int    p_typ, j_typ, k_typ;
  neightab *neigh;
  real   *tmpptr, *r2;
  vektor *d;
  real   tmp_d2, tmp_sp, tmp_pre, tmp_pot, pot, tmp_frc, tmp, tmp_phi1, tmp_phi;
  vektor dphi_j, dphi_k;
  tensor ddphi_jj, ddphi_jk, ddphi_kk;
  real   ddphi_jj_rr, dddphi_jjj_rrr;

  d    = (vektor *) malloc( neigh_len * sizeof(vektor) );
  r2   = (real *)   malloc( neigh_len * sizeof(real)   );

  if ((d==NULL) || (r2==NULL))
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
	  neigh   = p->neightab_array + i;

	  /* Construct some data for all neighbors */
	  tmpptr = neigh->dist;
	  for (j=0; j<neigh->n; ++j) 
	  {
	    /* Type, distance vector, radii */
	    j_typ   = neigh->typ[j];
	    d[j].x  = *tmpptr++;
	    d[j].y  = *tmpptr++;
	    d[j].z  = *tmpptr++;
	    r2[j]    = SPROD(d[j],d[j]);
     
	  } /* j */

	  /*************************************************************/

	  /* For each neighbor of i */
	  for (j=0; j<neigh->n; ++j)
	  {

	    j_typ = neigh->typ[j];
	    jcell = (cell *) neigh->cl[j];
	    jnum  = neigh->num[j];

	    /* Potential energy of pair potential part */
	    tmp_d2  = keat_d[p_typ][j_typ] * keat_d[p_typ][j_typ];
	    tmp_pre = 3.0 * keat_alpha[p_typ][j_typ] / ( 16.0 * tmp_d2 );
	    tmp_pot = r2[j] - tmp_d2;

	    pot  = tmp_pre * tmp_pot * tmp_pot;

	    epot += pot;

	    /* First derivative of pair potential part */
	    tmp_frc = 4.0 * tmp_pre * tmp_pot;

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
	    tmp_phi1 = 8.0 * tmp_pre;

	    ddphi_jj.xx = tmp_phi1 * d[j].x * d[j].x + tmp_frc;
	    ddphi_jj.xy = tmp_phi1 * d[j].x * d[j].y;
	    ddphi_jj.yy = tmp_phi1 * d[j].y * d[j].y + tmp_frc;
	    ddphi_jj.yz = tmp_phi1 * d[j].y * d[j].z;
	    ddphi_jj.zx = tmp_phi1 * d[j].z * d[j].x;
	    ddphi_jj.zz = tmp_phi1 * d[j].z * d[j].z + tmp_frc;

	    /* For computation of Bulk modulus */
	    ddphi_jj_rr = 8.0 * tmp_pre * ( 2.0 * r2[j] + tmp_pot );

	    bulkm += 0.5 * ddphi_jj_rr * r2[j];

	    /* For computation of B' */    
	    dddphi_jjj_rrr = 6.0 * tmp_phi1 * r2[j] * r2[j];

	    dbulkm_dp += 0.5 * dddphi_jjj_rrr;

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
	      kcell = (cell *) neigh->cl[k];
	      knum  = neigh->num[k];

	      /* Potential energy of three-body part */
	      tmp_sp  = SPROD(d[j],d[k]);
	      tmp_d2  = keat_d[p_typ][j_typ] * keat_d[p_typ][k_typ];
	      tmp_pre = 3.0 * keat_beta[p_typ][j_typ][k_typ] / ( 8.0 * tmp_d2 );
	      tmp_pot = tmp_sp + tmp_d2 / 3.0;

	      epot   += tmp_pre * tmp_pot * tmp_pot; 

	      /* First derivatives of three-body potential part */
	      tmp_frc  = 2.0 * tmp_pre * tmp_pot;

	      dphi_j.x = tmp_frc * d[k].x;
	      dphi_j.y = tmp_frc * d[k].y;
	      dphi_j.z = tmp_frc * d[k].z;

	      dphi_k.x = tmp_frc * d[j].x;
	      dphi_k.y = tmp_frc * d[j].y;
	      dphi_k.z = tmp_frc * d[j].z;     

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

	      /* Second derivatives of three-body potential part */
	      tmp_phi = 2.0 * tmp_pre;

	      ddphi_jj.xx = tmp_phi * d[k].x * d[k].x; 
	      ddphi_jj.xy = tmp_phi * d[k].x * d[k].y;
	      ddphi_jj.yy = tmp_phi * d[k].y * d[k].y;
	      ddphi_jj.yz = tmp_phi * d[k].y * d[k].z;
	      ddphi_jj.zx = tmp_phi * d[k].z * d[k].x;
	      ddphi_jj.zz = tmp_phi * d[k].z * d[k].z;

	      ddphi_jk.xx = tmp_phi * ( d[j].x * d[k].x + tmp_pot );
	      ddphi_jk.xy = tmp_phi * d[j].y * d[k].x;
	      ddphi_jk.xz = tmp_phi * d[j].z * d[k].x;
	      ddphi_jk.yx = tmp_phi * d[j].x * d[k].y;
	      ddphi_jk.yy = tmp_phi * ( d[j].y * d[k].y + tmp_pot );
	      ddphi_jk.yz = tmp_phi * d[j].z * d[k].y;
	      ddphi_jk.zx = tmp_phi * d[j].x * d[k].z;
	      ddphi_jk.zy = tmp_phi * d[j].y * d[k].z;
	      ddphi_jk.zz = tmp_phi * ( d[j].z * d[k].z + tmp_pot );

	      ddphi_kk.xx = tmp_phi * d[j].x * d[j].x; 
	      ddphi_kk.xy = tmp_phi * d[j].x * d[j].y;
	      ddphi_kk.yy = tmp_phi * d[j].y * d[j].y;
	      ddphi_kk.yz = tmp_phi * d[j].y * d[j].z;
	      ddphi_kk.zx = tmp_phi * d[j].z * d[j].x;
	      ddphi_kk.zz = tmp_phi * d[j].z * d[j].z;

	      /* For computation of bulk modulus */
	      bulkm += 4.0 * tmp_pre * tmp_sp * ( 2.0 * tmp_sp + tmp_pot );

	      /* For computation of B' */	      
	      dbulkm_dp += 24.0 * tmp_pre * tmp_sp * tmp_sp;

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
