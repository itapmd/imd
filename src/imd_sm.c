
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2011 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/*****************************************************************************
*
* imd_sm.c -- Minimization calculations of the Streitz and Mintmire model
*
******************************************************************************/

/*****************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"
#include "potaccess.h"

/*****************************************************************************
*
* Read in nuclear attraction potential, coulomb repulsive potential, and 
* tabulated error function 
*
******************************************************************************/

void init_sm(void)
{
  read_pot_table(&na_pot_tab,na_pot_filename,ntypes*ntypes,1);
  read_pot_table(&cr_pot_tab,cr_pot_filename,ntypes*ntypes,1);
  read_pot_table(&erfc_r_tab,erfc_filename,ntypes*ntypes,1);
}


/*****************************************************************************
*
* Compute the electronegativity
*
******************************************************************************/

void do_electronegativity(void)
{
  int i, j, k, r, s, jstart;
  vektor d, tmp_d;
  real r2;
  int is_short=0;
  int col1, col2, inc=ntypes*ntypes;
  int q_typ, p_typ;
  cell *p, *q;
  real na_pot_p, na_pot_q, cr_pot;
  real chi_sm_p, chi_sm_q, z_sm_p, z_sm_q;

  /* Inititalisierung aller Zellen */
  for(r=0; r<ncells; ++r)
    {
      p = CELLPTR(r);
      for (i=0; i<p->n; ++i) CHI_SM(p,i) = 0;
    }
	  
  /* Loop over all pairs of cells */
  for(r=0; r<ncells; ++r)
    for(s=r; s<ncells; ++s)
      {
	p = CELLPTR(r);
	q = CELLPTR(s);
	
	/* for each atom in first cell */
	for (i=0; i<p->n; ++i) 
	  {
	    p_typ   = SORTE(p,i);
	    chi_sm_p = chi_0[p_typ];
	    z_sm_p = z_es[p_typ];
	    
	    /* For each atom in second cell */
	    jstart = (p==q ? i+1 : 0);
	    
	    for (j=jstart; j<q->n; ++j) {
	      
	      /* calculate distance */ 
	      d.x = ORT(p,i,X) - ORT(q,j,X);
	      d.y = ORT(p,i,Y) - ORT(q,j,Y);
	      d.z = ORT(p,i,Z) - ORT(q,j,Z);
	      
	      q_typ = SORTE(q,j);
	      chi_sm_q = chi_0[q_typ];
	      z_sm_q = z_es[q_typ];

	      r2    = SPROD(d,d);
	      col1  = q_typ * ntypes + p_typ;
	      col2  = p_typ * ntypes + q_typ;
	      
#ifdef DEBUG
	      if (0==r2) { char msgbuf[256];
		sprintf(msgbuf, "Distance 1 is zero in module sm between particles %d and %d!\n",
			NUMMER(p,i), NUMMER(q,j));
		error(msgbuf);
	      }
#endif
	      
	      /* compute electronegativity */
	      if (r2 < na_pot_tab.end[col2]){
		VAL_FUNC(na_pot_p, na_pot_tab, col2, inc, r2, is_short);
	      }
	      if (r2 < na_pot_tab.end[col1]){
		VAL_FUNC(na_pot_q, na_pot_tab, col1, inc, r2, is_short);
	      }
	      if (r2 < cr_pot_tab.end[col2]){
		VAL_FUNC(cr_pot, cr_pot_tab, col2, inc, r2, is_short);
	      }
	      chi_sm_p += z_sm_q*(na_pot_p-cr_pot);
	      CHI_SM(p,i) += chi_sm_p;

	      chi_sm_q += z_sm_p*(na_pot_q-cr_pot);
	      CHI_SM(q,j) += chi_sm_q;
	    }
	  }
      }
}

/*****************************************************************************
*
* Computes the real space part of v_i = V_ijq_j
*
******************************************************************************/

void do_v_real(void)
{

  int i, j, k, r, s, jstart;
  vektor d, tmp_d;
  real r2;
  int is_short=0;
  int col1, col2, inc=ntypes*ntypes;
  int q_typ, p_typ;
  cell *p, *q;
  real erfc_r, cr_pot;
  real j_sm_p, j_sm_q, v_sm_p, v_sm_q;

  /* Inititalisierung aller Zellen */
  for(r=0; r<ncells; ++r)
    {
      p = CELLPTR(r);
      for (i=0; i<p->n; ++i) V_SM(p,i) = 0;
    }
	  
  /* Loop over all pairs of cells */
  for(r=0; r<ncells; ++r)
    for(s=r; s<ncells; ++s)
      {
	p = CELLPTR(r);
	q = CELLPTR(s);
	
	/* for each atom in first cell */
	for (i=0; i<p->n; ++i) 
	  {
	    p_typ   = SORTE(p,i);
	    j_sm_p  = j_0[p_typ];
	    v_sm_p = Q_SM(p,i)*(j_sm_p-ew_vorf);

	    /* For each atom in second cell */
	    jstart = (p==q ? i+1 : 0);
	    
	    /* for each atom in neighbouring cell */
	    for (j=jstart; j<q->n; ++j) {
	      
	      /* calculate distance */ 
	      d.x = ORT(p,i,X) - ORT(q,j,X);
	      d.y = ORT(p,i,Y) - ORT(q,j,Y);
	      d.z = ORT(p,i,Z) - ORT(q,j,Z);

	      /* minimum image convention? */
	      
	      q_typ = SORTE(q,j);
	      j_sm_q  = j_0[q_typ];
	      v_sm_q = Q_SM(q,i)*(j_sm_q-ew_vorf);

	      r2    = SPROD(d,d);
	      col1  = q_typ * ntypes + p_typ;
	      col2  = p_typ * ntypes + q_typ;
	      
#ifdef DEBUG
	      if (0==r2) { char msgbuf[256];
		sprintf(msgbuf, "Distance 2 is zero in module between particles sm %d and %d!\n",
			NUMMER(p,i), NUMMER(q,j));
		error(msgbuf);
	      }
#endif
	      
	      /* compute real space term of v_i */
	      if (r2 < erfc_r_tab.end[col2]){
		VAL_FUNC(erfc_r, erfc_r_tab, col2, inc, r2, is_short);
	      }
	      if (r2 < cr_pot_tab.end[col2]){
		VAL_FUNC(cr_pot, cr_pot_tab, col2, inc, r2, is_short);
	      }
	      v_sm_p += Q_SM(q,j)*(erfc_r+cr_pot);
	      V_SM(p,i) += v_sm_p;

	      v_sm_q += Q_SM(p,j)*(erfc_r+cr_pot);
	      V_SM(q,j) += v_sm_q;
	    }
	  }
      }
}

/******************************************************************************
*
*  do_v_kspace
*
*  computes the fourier part of the Ewald sum
*
******************************************************************************/

void do_v_kspace(void)
{
  int    i, j, k, l, m, n, c;
  int    cnt, typ, offx, offy, offz;
  int    px, py, pz, mx, my, mz;
  real   tmp, tmp_virial=0.0, sum_cos, sum_sin;
  real   kforce, kpot;
  real v_k;

  /* Compute exp(ikr) recursively */
  px = (ew_nx+1) * natoms;  mx = (ew_nx-1) * natoms;
  py = (ew_ny+1) * natoms;  my = (ew_ny-1) * natoms;
  pz = (ew_nz+1) * natoms;  mz = (ew_nz-1) * natoms;
  cnt = 0;
  for (c=0; c<ncells; c++) {
    cell *p = CELLPTR(c);
    for (i=0; i<p->n; i++) {
      tmp = twopi * SPRODX(ORT,p,i,tbox_x);
      coskx[px+cnt] =  cos(tmp);
      coskx[mx+cnt] =  coskx[px+cnt];
      sinkx[px+cnt] =  sin(tmp);
      sinkx[mx+cnt] = -sinkx[px+cnt]; 
      tmp = twopi * SPRODX(ORT,p,i,tbox_y);
      cosky[py+cnt] =  cos(tmp);
      cosky[my+cnt] =  cosky[py+cnt];
      sinky[py+cnt] =  sin(tmp);
      sinky[my+cnt] = -sinky[py+cnt];
      tmp = twopi * SPRODX(ORT,p,i,tbox_z);
      coskz[pz+cnt] = cos(tmp);
      sinkz[pz+cnt] = sin(tmp);
      cnt++;
    }
  }

  for (j=2; j<=ew_nx; j++) {
    int pp, qq, mm, ee, i;
    pp  = (ew_nx+j  ) * natoms;
    qq  = (ew_nx+j-1) * natoms;
    mm  = (ew_nx-j  ) * natoms;
    ee  = (ew_nx  +1) * natoms;
    for (i=0; i<natoms; i++) {
      coskx[pp+i] =   coskx[qq+i] * coskx[ee+i] - sinkx[qq+i] * sinkx[ee+i];
      coskx[mm+i] =   coskx[pp+i];
      sinkx[pp+i] =   coskx[qq+i] * sinkx[ee+i] - sinkx[qq+i] * coskx[ee+i];
      sinkx[mm+i] = - sinkx[pp+i];
    }
  }

  for (j=2; j<=ew_ny; j++) {
    int pp, qq, mm, ee, i;
    pp  = (ew_ny+j  ) * natoms;
    qq  = (ew_ny+j-1) * natoms;
    mm  = (ew_ny-j  ) * natoms;
    ee  = (ew_ny  +1) * natoms;
    for (i=0; i<natoms; i++) {
      cosky[pp+i] =   cosky[qq+i] * cosky[ee+i] - sinky[qq+i] * sinky[ee+i];
      cosky[mm+i] =   cosky[pp+i];
      sinky[pp+i] =   cosky[qq+i] * sinky[ee+i] - sinky[qq+i] * cosky[ee+i];
      sinky[mm+i] = - sinky[pp+i];
    }
  }

  for (j=2; j<=ew_nz; j++) {
    int pp, qq, ee, i;
    pp  = (ew_nz+j  ) * natoms;
    qq  = (ew_nz+j-1) * natoms;
    ee  = (ew_nz  +1) * natoms;
    for (i=0; i<natoms; i++) {
      coskz[pp+i] =   coskz[qq+i] * coskz[ee+i] - sinkz[qq+i] * sinkz[ee+i];
      sinkz[pp+i] =   coskz[qq+i] * sinkz[ee+i] - sinkz[qq+i] * coskz[ee+i];
    }
  }

  /* Loop over all reciprocal vectors */
  for (k=0; k<ew_totk; k++) {

    cnt     = 0; 
    offx    = ew_ivek[k].x * natoms;
    offy    = ew_ivek[k].y * natoms;
    offz    = ew_ivek[k].z * natoms;
    sum_cos = 0.0;
    sum_sin = 0.0;
    v_k     = 0.0;

    /* Compute exp(ikr) and sums thereof */
    for (c=0; c<ncells; c++) {

      cell *p = CELLPTR(c);

      for (i=0; i<p->n; i++) {

        coskr[cnt] =   coskx[offx+cnt] * cosky[offy+cnt] * coskz[offz+cnt]
                     - sinkx[offx+cnt] * sinky[offy+cnt] * coskz[offz+cnt]
                     - sinkx[offx+cnt] * cosky[offy+cnt] * sinkz[offz+cnt] 
                     - coskx[offx+cnt] * sinky[offy+cnt] * sinkz[offz+cnt];

        sinkr[cnt] = - sinkx[offx+cnt] * sinky[offy+cnt] * sinkz[offz+cnt]
                     + sinkx[offx+cnt] * cosky[offy+cnt] * coskz[offz+cnt]
                     + coskx[offx+cnt] * sinky[offy+cnt] * coskz[offz+cnt]
                     + coskx[offx+cnt] * cosky[offy+cnt] * sinkz[offz+cnt];

        sum_cos += Q_SM(p,i) * coskr[cnt];
        sum_sin += Q_SM(p,i) * sinkr[cnt];
        cnt++;
      }
    }

    /* updates */
    cnt = 0;
    for (c=0; c<ncells; c++) {

      cell *p = CELLPTR(c);

      for (i=0; i<p->n; i++) {

	/* update fourier part of v_i */
        v_k   = ew_expk[k] * (sinkr[cnt] * sum_sin + coskr[cnt] * sum_cos);
        V_k_SM(p,i)  += v_k;   
	/* the total vector v_i */
	V_SM(p,i) += v_k;
        cnt++;
      }
    }
  }  /* k */
}
/*****************************************************************************
*
* Conjugate gradient algorithm for solving the system Ax=b
*
******************************************************************************/

void do_cg(void)
{
  int k, kstep=0, kstepmax=1000;
  real beta, alpha, rho[kstepmax];
  real dad, tolerance, tolerance2, norm_b, epsilon_cg=0.001;
  
  /* compute V_SM and V_K_SM */
  do_v_real();
  do_v_kspace();
      
  /*  norm_b +=SPRODN(B,p,i,B,p,i); */
  norm_b=0.0;
  for (k=0; k<NCELLS; ++k) {
    int  i;
    cell *p;
    p = CELLPTR(k);
    /* loop over all particles */
    for (i=0; i<p->n; ++i) {      
      norm_b += B_SM(p,i)*B_SM(p,i);

      /* initial values */
      X_SM(p,i) = Q_SM(p,i);
      R_SM(p,i) = B_SM(p,i)-V_SM(p,i);
    }
  }
  
  tolerance = epsilon_cg*SQRT(norm_b);
  tolerance2 = SQR(tolerance);
  
  rho[kstep] = 0.0;
  /*  rho[kstep] += SPRODN(R,p,i,R,p,i); */
  for (k=0; k<NCELLS; ++k) {
    int  i;
    cell *p;
    p = CELLPTR(k);
    /* loop over all particles */
    for (i=0; i<p->n; ++i) {      
      rho[kstep] += R_SM(p,i)*R_SM(p,i);
    }
  }
  
  while ((rho[kstep] > tolerance2) && (kstep < kstepmax))
    {
      kstep++;

      if (kstep == 1)
	for (k=0; k<NCELLS; ++k) {
	  int  i;
	  cell *p;
	  p = CELLPTR(k);
	  /* loop over all particles */
	  for (i=0; i<p->n; ++i) {      
	    Q_SM(p,i) = D_SM(p,i) = R_SM(p,i);
	  }
	}
      else
	{
	  beta = rho[kstep-1]/rho[kstep-2];	  
	  for (k=0; k<NCELLS; ++k) {
	    int  i;
	    cell *p;
	    p = CELLPTR(k);
	    /* loop over all particles */
	    for (i=0; i<p->n; ++i) {      
	      Q_SM(p,i) = D_SM(p,i) = R_SM(p,i) + beta* D_SM(p,i);

	    }
	  }
	}
      
      /* update V_SM and V_K_SM */
      do_v_real();
      do_v_kspace();
      
      dad = 0.0;
      /* dad += SPRODN(D,p,i,V,p,i); */
      for (k=0; k<NCELLS; ++k) {
	int  i;
	cell *p;
	p = CELLPTR(k);
	/* loop over all particles */
	for (i=0; i<p->n; ++i) {      
	  dad += D_SM(p,i)*V_SM(p,i);
	}
      }
      
      alpha = rho[kstep-1]/dad;

      for (k=0; k<NCELLS; ++k) {
	int  i;
	cell *p;
	p = CELLPTR(k);
	/* loop over all particles */
	for (i=0; i<p->n; ++i) {      
	  X_SM(p,i) += alpha*D_SM(p,i);
	  R_SM(p,i) -= alpha*V_SM(p,i);
	}
      }	  

      rho[kstep] = 0.0;
      /* rho[kstep]+= SPRODN(R,p,i,R,p,i); */
      for (k=0; k<NCELLS; ++k) {
	int  i;
	cell *p;
	p = CELLPTR(k);
	/* loop over all particles */
	for (i=0; i<p->n; ++i) {      
	  rho[kstep] += R_SM(p,i)*R_SM(p,i);
	}
      }
    }

}

/*****************************************************************************
*
* Charge update
*
******************************************************************************/

void do_charge_update(void)
{
  
  int k, typ;
  real sum1, sum2, potchem;
  real q_Al, q_O, q_tot;
  int n_O, n_Al;
  
  /* Update electronegativity since coordinates have changed */

  do_electronegativity();

  /* Solving the first linear system V_ij s_j = -chi_i */
  
  /* loop over all cells */
  for (k=0; k<NCELLS; ++k) {
    int  i;
    cell *p;
    p = CELLPTR(k);
    
    /* loop over all particles */
    for (i=0; i<p->n; ++i) {
      
      B_SM(p,i) = -CHI_SM(p,i);
      /* Initial value of the charges */
      Q_SM(p,i) = CHARGE(p,i);
    }
  }
  
  do_cg();
  
  /* Sum up for getting charges */
  sum1=0.0;

  /* loop over all cells */
  for (k=0; k<NCELLS; ++k) {
    int  i;
    cell *p;
    p = CELLPTR(k);
    
    /* loop over all particles */
    for (i=0; i<p->n; ++i) {
      sum1 += X_SM(p,i);
      S_SM(p,i) = X_SM(p,i);
    }
  }
  
  /* Solving the second linear system V_ij t_j = -1 */

  /* loop over all cells */
  for (k=0; k<NCELLS; ++k) {
    int  i;
    cell *p;
    p = CELLPTR(k);
    
    /* loop over all particles */
    for (i=0; i<p->n; ++i) {
      
      B_SM(p,i) = -1;
      /* Initial value of the charges */
      Q_SM(p,i) = 0.0;
    }
  }

  do_cg();
  
  /* Sum up for getting charges */
  sum2=0.0;
  
  /* loop over all cells */
  for (k=0; k<NCELLS; ++k) {
    int  i;
    cell *p;
    p = CELLPTR(k);
    
    /* loop over all particles */
    for (i=0; i<p->n; ++i) {
      sum2 += X_SM(p,i);
    }
  }
  
  /* Calculate chemical potential */

  potchem = sum1/sum2;
      
  /* Calculate new charges of the atoms */

  q_tot = 0.0;
  q_Al = 0.0;
  q_O = 0.0;
  n_Al = 0;
  n_O = 0;

  /* loop over all cells */
  for (k=0; k<NCELLS; ++k) {
    int  i;
    cell *p;
    p = CELLPTR(k);
    
    /* loop over all particles */
    for (i=0; i<p->n; ++i) {
      
      CHARGE(p,i) = S_SM(p,i)-potchem*X_SM(p,i);
      q_tot += CHARGE(p,i);
      typ = SORTE(p,i);
      if (typ == 0) {
	q_Al += CHARGE(p,i);
	n_Al++;
      }
      else {
	q_O += CHARGE(p,i);
	n_O++;
      }
    }
  }

  if ((n_Al != 0) || (n_O != 0))
    {
      q_Al/n_Al;
      q_O/n_O;
    }
  
  printf("Sums: sum1 = %e sum2 = %e\n", sum1, sum2);
  printf("Total charge: qtot = %e\n", q_tot);
#ifndef DEBUG
  printf("Total charge and number of Al: qAl = %e\n", q_Al/n_Al);
  printf("Total charge and number of O: qO = %e\n", q_O/n_O);
#else
  printf("Average charge of Al: qAl = %e %d\n", q_Al,n_Al);
  printf("Average charge of O: qO = %e %d\n", q_O,n_O);
#endif
}

