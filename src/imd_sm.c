
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
*
*
******************************************************************************/

#include "imd.h"
#include "potaccess.h"

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
  real na_pot, cr_pot;
  real chi_sm, z_sm;

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
	    chi_sm = chi_0[p_typ];
	    
	    /* For each atom in second cell */
	    jstart = (p==q ? i : 0);
	    
	    for (j=jstart; j<q->n; ++j) {
	      
	      /* calculate distance */ 
	      d.x = ORT(p,i,X) - ORT(q,j,X);
	      d.y = ORT(p,i,Y) - ORT(q,j,Y);
	      d.z = ORT(p,i,Z) - ORT(q,j,Z);
	      
	      q_typ = SORTE(q,j);
	      z_sm = z_es[q_typ];
	      r2    = SPROD(d,d);
	      col1  = q_typ * ntypes + p_typ;
	      col2  = p_typ * ntypes + q_typ;
	      
#ifdef DEBUG
	      if (0==r2) { char msgbuf[256];
		sprintf(msgbuf, "Distance is zero between particles %d and %d!\n",
			NUMMER(p,i), NUMMER(q,j));
		error(msgbuf);
	      }
#endif
	      
	      /* compute electronegativity */
	      if (r2 < na_pot_tab.end[col2]){
		VAL_FUNC(na_pot, na_pot_tab, col2, inc, r2, is_short);
	      }
	      if (r2 < cr_pot_tab.end[col2]){
		VAL_FUNC(cr_pot, cr_pot_tab, col2, inc, r2, is_short);
	      }
	      chi_sm += z_sm*(na_pot-cr_pot);
	      CHI_SM(p,i) = chi_sm;
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
  real q_0, j_sm, v_sm;

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
	    j_sm  = j_0[p_typ];
	    /* Initial value of the charges */
	    CHARGE(p,i)= q_0;
	    v_sm = CHARGE(p,i)*(j_sm-ew_vorf);

	    /* For each atom in second cell */
	    jstart = (p==q ? i : 0);
	    
	    /* for each atom in neighbouring cell */
	    for (j=jstart; j<q->n; ++j) {
	      
	      /* calculate distance */ 
	      d.x = ORT(p,i,X) - ORT(q,j,X);
	      d.y = ORT(p,i,Y) - ORT(q,j,Y);
	      d.z = ORT(p,i,Z) - ORT(q,j,Z);
	      
	      q_typ = SORTE(q,j);
	      CHARGE(q,j) = q_0;
	      r2    = SPROD(d,d);
	      col1  = q_typ * ntypes + p_typ;
	      col2  = p_typ * ntypes + q_typ;
	      
#ifdef DEBUG
	      if (0==r2) { char msgbuf[256];
		sprintf(msgbuf, "Distance is zero between particles %d and %d!\n",
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
	      v_sm += CHARGE(q,j)*(erfc_r+cr_pot);
	      V_SM(p,i) = v_sm;
	    }
	  }
      }
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
  
  /*  norm_b +=SPRODN(B,p,i,B,p,i); */
  norm_b=0.0;
  for (k=0; k<NCELLS; ++k) {
    int  i, sort;
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
    int  i, sort;
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
	  int  i, sort;
	  cell *p;
	  p = CELLPTR(k);
	  /* loop over all particles */
	  for (i=0; i<p->n; ++i) {      
	    D_SM(p,i) = R_SM(p,i);
	  }
	}
      else
	{
	  beta = rho[kstep-1]/rho[kstep-2];	  
	  for (k=0; k<NCELLS; ++k) {
	    int  i, sort;
	    cell *p;
	    p = CELLPTR(k);
	    /* loop over all particles */
	    for (i=0; i<p->n; ++i) {      
	      D_SM(p,i) = R_SM(p,i) + beta* D_SM(p,i);
	    }
	  }
	}
      
      dad = 0.0;
      /* dad += SPRODN(D,p,i,V,p,i); */
      for (k=0; k<NCELLS; ++k) {
	int  i, sort;
	cell *p;
	p = CELLPTR(k);
	/* loop over all particles */
	for (i=0; i<p->n; ++i) {      
	  dad += D_SM(p,i)*V_SM(p,i);
	}
      }
      
      alpha = rho[kstep-1]/dad;

      for (k=0; k<NCELLS; ++k) {
	int  i, sort;
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
	int  i, sort;
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
  
  /* Compute V_SM, V_k_SM */
  do_v_real();
  do_forces_ewald_fourier();

  /* Solving the first linear system V_ij s_j = -chi_i */
  
  /* loop over all cells */
  for (k=0; k<NCELLS; ++k) {
    int  i, sort;
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
    int  i, sort;
    cell *p;
    p = CELLPTR(k);
    
    /* loop over all particles */
    for (i=0; i<p->n; ++i) {
      S_SM(p,i) = X_SM(p,i);
      sum1 += S_SM(p,i);
    }
  }
  
  /* Solving the second linear system V_ij t_j = -1 */

  /* loop over all cells */
  for (k=0; k<NCELLS; ++k) {
    int  i, sort;
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
    int  i, sort;
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
    int  i, sort;
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
  
  printf("Total charge: qtot = %e", q_tot);
  printf("Average charge of Al: qAl = %e", q_Al/n_Al);
  printf("Average charge of O: qO = %e", q_O/n_O);
}

