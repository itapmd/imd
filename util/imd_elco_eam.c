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
* imd_elco_eam
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "util.h"
#include "potaccess.h"

/******************************************************************************
*
*  init_eam, Initializations
*
******************************************************************************/

void init_eam(void) 
{
  /* read the tabulated embedding energy function */
  read_pot_table(&embed_pot,eam_emb_E_filename,ntypes);

  /* read the tabulated electron density function */
  read_pot_table(&rho_h_tab,eam_at_rho_filename,ntypes*ntypes);
}

/******************************************************************************
*
*  do_elco_eam -- computes stress tensor and elastic constants using the 
*             EAM potential
*
******************************************************************************/

void do_elco_eam(cell *p, cell *q, vektor pbc)
{
  int i, j, jstart;
  vektor tmp_d, d;
  int same_cell;
  int p_typ, q_typ, col1, col2;
  real r2, eam_energy, pref1, pref2, tmp;
  real f_i_strich, f_i_zweistrich, f_i_dreistrich;
  real f_j_strich, f_j_zweistrich, f_j_dreistrich;
  real rho_j_strich, rho_j_zweistrich, rho_j_dreistrich;
  real rho_i_strich, rho_i_zweistrich, rho_i_dreistrich;
  int is_short=0, inc=ntypes*ntypes;
  real xx, xy, yy, yz, zz, zx;

  /* for each atom in first cell */
  for (i=0; i<p->n; ++i) {

    tmp_d.x = p->ort[i].x - pbc.x;
    tmp_d.y = p->ort[i].y - pbc.y;
#ifndef TWOD
    tmp_d.z = p->ort[i].z - pbc.z;
#endif
    p_typ   = p->sorte[i];

#ifdef TWOD
    same_cell = ((p==q) && (pbc.x==0) && (pbc.y==0));
#else
    same_cell = ((p==q) && (pbc.x==0) && (pbc.y==0) && (pbc.z==0));
#endif

    if (same_cell) {
      PAIR_INT4(eam_energy, f_i_strich, f_i_zweistrich, f_i_dreistrich, 
	       embed_pot, p_typ, ntypes, EAM_RHO(p,i), is_short)

      epot += eam_energy;

      xx = p->eam_stress[i].xx;
      xy = p->eam_stress[i].xy;
      yy = p->eam_stress[i].yy;
#ifndef TWOD
      yz = p->eam_stress[i].yz;
      zz = p->eam_stress[i].zz;
      zx = p->eam_stress[i].zx;
#endif

      p->stress[i].xx += f_i_strich * xx;
      sigma.xx        += f_i_strich * xx;
      p->stress[i].xy += f_i_strich * xy;
      sigma.xy        += f_i_strich * xy;
      p->stress[i].yy += f_i_strich * yy;
      sigma.yy        += f_i_strich * yy;
#ifndef TWOD
      p->stress[i].yz += f_i_strich * yz;
      sigma.yz        += f_i_strich * yz;
      p->stress[i].zz += f_i_strich * zz;
      sigma.zz        += f_i_strich * zz;
      p->stress[i].zx += f_i_strich * zx;
      sigma.zx        += f_i_strich * zx;
#endif

      p->elco[i].c11  += f_i_zweistrich * xx * xx;
      c.c11           += f_i_zweistrich * xx * xx;
      p->elco[i].c12  += f_i_zweistrich * xx * yy;
      c.c12           += f_i_zweistrich * xx * yy;
      p->elco[i].c66  += f_i_zweistrich * xy * xy;
      c.c66           += f_i_zweistrich * xy * xy;
      p->elco[i].c22  += f_i_zweistrich * yy * yy;
      c.c22           += f_i_zweistrich * yy * yy;
#ifndef TWOD
      p->elco[i].c13  += f_i_zweistrich * xx * zz;
      c.c13           += f_i_zweistrich * xx * zz;
      p->elco[i].c55  += f_i_zweistrich * zx * zx;
      c.c55           += f_i_zweistrich * zx * zx;
      p->elco[i].c23  += f_i_zweistrich * yy * zz;
      c.c23           += f_i_zweistrich * yy * zz;
      p->elco[i].c44  += f_i_zweistrich * yz * yz;
      c.c44           += f_i_zweistrich * yz * yz;
      p->elco[i].c33  += f_i_zweistrich * zz * zz;
      c.c33           += f_i_zweistrich * zz * zz;
      p->elco[i].c45  += f_i_zweistrich * yz * zx;
      c.c45           += f_i_zweistrich * yz * zx;
      p->elco[i].c46  += f_i_zweistrich * yz * xy;
      c.c46           += f_i_zweistrich * yz * xy;
      p->elco[i].c25  += f_i_zweistrich * yy * zx;
      c.c25           += f_i_zweistrich * yy * zx;
      p->elco[i].c56  += f_i_zweistrich * zx * xy;
      c.c56           += f_i_zweistrich * zx * xy;
      p->elco[i].c14  += f_i_zweistrich * xx * yz;
      c.c14           += f_i_zweistrich * xx * yz;
      p->elco[i].c15  += f_i_zweistrich * xx * zx;
      c.c15           += f_i_zweistrich * xx * zx;
#endif
      p->elco[i].c16  += f_i_zweistrich * xx * xy;
      c.c16           += f_i_zweistrich * xx * xy;
      p->elco[i].c26  += f_i_zweistrich * yy * xy;
      c.c26           += f_i_zweistrich * yy * xy;
#ifndef TWOD
      p->elco[i].c24  += f_i_zweistrich * yy * yz;
      c.c24           += f_i_zweistrich * yy * yz;
      p->elco[i].c34  += f_i_zweistrich * zz * yz;
      c.c34           += f_i_zweistrich * zz * yz;
      p->elco[i].c35  += f_i_zweistrich * zz * zx;
      c.c35           += f_i_zweistrich * zz * zx;
      p->elco[i].c36  += f_i_zweistrich * zz * xy;
      c.c36           += f_i_zweistrich * zz * xy;
#endif
      press           += f_i_strich * p->eam_press[i];

      bulkm           += f_i_zweistrich * p->eam_press[i] * p->eam_press[i]
	                   + f_i_strich * p->eam_bulkm[i];

      dbulkm_dp       += f_i_dreistrich 
	                 * p->eam_press[i] * p->eam_press[i] * p->eam_press[i]
			   + 3.0 * f_i_zweistrich 
			     * p->eam_press[i] * p->eam_bulkm[i]
			       + f_i_strich * p->eam_dbulkm[i];
    } 
    else {      
      DERIV_FUNC(f_i_strich, f_i_zweistrich, f_i_dreistrich, 
		 embed_pot, p_typ, ntypes, EAM_RHO(p,i), is_short)
    }

    jstart = (same_cell ? i+1 : 0);

    /* for each atom in neighbouring cell */
    for (j=jstart; j<q->n; ++j) {

      /* calculate distance */ 
      d.x = q->ort[j].x - tmp_d.x;
      d.y = q->ort[j].y - tmp_d.y;
#ifndef TWOD
      d.z = q->ort[j].z - tmp_d.z;
#endif
      q_typ = q->sorte[j];
      r2    = SPROD(d,d);
      col1  = q_typ * ntypes + p_typ;
      col2  = p_typ * ntypes + q_typ;

      if ((r2 < rho_h_tab.end[col1]) || (r2 < rho_h_tab.end[col2])) {

        DERIV_FUNC(f_j_strich, f_j_zweistrich, f_j_dreistrich, 
		   embed_pot, q_typ, ntypes, EAM_RHO(q,j), is_short);

        DERIV_FUNC(rho_i_strich, rho_i_zweistrich, rho_i_dreistrich, 
		   rho_h_tab, col1, inc, r2, is_short);

        if ( col1 == col2 ) {
          rho_j_strich     = rho_i_strich;
          rho_j_zweistrich = rho_i_zweistrich;
	} else {

          DERIV_FUNC(rho_j_strich, rho_j_zweistrich, rho_j_dreistrich, 
		     rho_h_tab, col2, inc, r2, is_short);
	}

	pref1 = 2.0 * 
	  ( f_i_strich * rho_j_zweistrich + f_j_strich * rho_i_zweistrich );

	pref2 = f_i_strich * rho_j_strich + f_j_strich * rho_i_strich;

	tmp = d.x * d.x * d.x * d.x * pref1 + d.x * d.x * pref2;
	p->elco[i].c11 += tmp;
	q->elco[j].c11 += tmp;
	c.c11          += 2.0 * tmp;
	tmp = d.x * d.x * d.y * d.y * pref1;
	p->elco[i].c12 += tmp;
	q->elco[j].c12 += tmp;
	c.c12          += 2.0 * tmp;
	tmp += 0.5 * ( d.x * d.x + d.y * d.y ) * pref2;
	p->elco[i].c66 += tmp;
	q->elco[j].c66 += tmp;
	c.c66          += 2.0 * tmp;
	tmp = d.y * d.y * d.y * d.y * pref1 + d.y * d.y * pref2;
	p->elco[i].c22 += tmp;
	q->elco[j].c22 += tmp;
	c.c22          += 2.0 * tmp;
#ifndef TWOD
	tmp = d.x * d.x * d.z * d.z * pref1;
	p->elco[i].c13 += tmp;
	q->elco[j].c13 += tmp;
	c.c13          += 2.0 * tmp;
	tmp += 0.5 * ( d.x * d.x + d.z * d.z ) * pref2; 
	p->elco[i].c55 += tmp;
	q->elco[j].c55 += tmp;
	c.c55          += 2.0 * tmp;
	tmp = d.y * d.y * d.z * d.z * pref1;
	p->elco[i].c23 += tmp;
	q->elco[j].c23 += tmp;
	c.c23          += 2.0 * tmp;
	tmp += 0.5 * ( d.y * d.y + d.z * d.z ) * pref2; 
	p->elco[i].c44 += tmp;
	q->elco[j].c44 += tmp;
	c.c44          += 2.0 * tmp;
	tmp = d.z * d.z * d.z * d.z * pref1 + d.z * d.z * pref2;
	p->elco[i].c33 += tmp;
	q->elco[j].c33 += tmp;
	c.c33          += 2.0 * tmp;
	tmp = d.x * d.y * d.z * d.z * pref1 + 0.25 * d.x * d.y * pref2;
	p->elco[i].c45 += tmp;
	q->elco[j].c45 += tmp;
	c.c45          += 2.0 * tmp;
	tmp = d.x * d.y * d.y * d.z * pref1 + 0.25 * d.x * d.z * pref2;
	p->elco[i].c46 += tmp;
	q->elco[j].c46 += tmp;
	c.c46          += 2.0 * tmp;
	tmp -= 0.25 * d.x * d.z * pref2;
	p->elco[i].c25 += tmp;
	q->elco[j].c25 += tmp;
	c.c25          += 2.0 * tmp;
	tmp = d.x * d.x * d.y * d.z * pref1 + 0.25 * d.y * d.z * pref2;
	p->elco[i].c56 += tmp;
	q->elco[j].c56 += tmp;
	c.c56          += 2.0 * tmp;
	tmp -= 0.25 * d.y * d.z * pref2;
	p->elco[i].c14 += tmp;
	q->elco[j].c14 += tmp;
	c.c14          += 2.0 * tmp;
	tmp = d.x * d.x * d.x * d.z * pref1 + 0.5 * d.x * d.z * pref2;
	p->elco[i].c15 += tmp;
	q->elco[j].c15 += tmp;
	c.c15          += 2.0 * tmp;
#endif
	tmp = d.x * d.x * d.x * d.y * pref1 + 0.5 * d.x * d.y * pref2;
	p->elco[i].c16 += tmp;
	q->elco[j].c16 += tmp;
	c.c16          += 2.0 * tmp;
	tmp = d.x * d.y * d.y * d.y * pref1 + 0.5 * d.x * d.y * pref2;
	p->elco[i].c26 += tmp;
	q->elco[j].c26 += tmp;
	c.c26          += 2.0 * tmp;
#ifndef TWOD
	tmp = d.y * d.y * d.y * d.z * pref1 + 0.5 * d.y * d.z * pref2;
	p->elco[i].c24 += tmp;
	q->elco[j].c24 += tmp;
	c.c24          += 2.0 * tmp;
	tmp = d.y * d.z * d.z * d.z * pref1 + 0.5 * d.y * d.z * pref2;
	p->elco[i].c34 += tmp;
	q->elco[j].c34 += tmp;
	c.c34          += 2.0 * tmp;
	tmp = d.x * d.z * d.z * d.z * pref1 + 0.5 * d.x * d.z * pref2;
	p->elco[i].c35 += tmp;
	q->elco[j].c35 += tmp;
	c.c35          += 2.0 * tmp;
	tmp = d.x * d.y * d.z * d.z * pref1;
	p->elco[i].c36 += tmp;
	q->elco[j].c36 += tmp;
	c.c36          += 2.0 * tmp;
#endif
      }
    }
  }
}
