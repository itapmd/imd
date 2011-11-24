
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
#ifdef DEBUG
  printf("init_sm\n");
#ifndef COULOMB
  printf("COULOMB not defined\n");
#else
  printf("COULOMB defined!!\n");
#endif
  printf("computing initial charges ew_nmax=%d\n",ew_nmax);
#endif

#ifdef NBLIST
#ifdef MPI
  setup_buffers();  /* setup MPI buffers */
#endif
  /* fill the buffer cells for first charge update */
  if ((!sm_fixed_charges) && ((charge_update_steps > 0) && steps % charge_update_steps == 0)) {
    send_cells(copy_cell,pack_cell,unpack_cell);
    charge_update_sm();
  }
#else
  if ((!sm_fixed_charges) && ((charge_update_steps > 0) && steps % charge_update_steps == 0))
    do_charge_update();
#endif
}


/*****************************************************************************
*
* Compute the electronegativity
*
******************************************************************************/

void do_electronegativity(void)
{
#ifdef DEBUG
  printf("do_electronegativity\n");
#endif
  int i, j, k, r, s, jstart;
  vektor d, tmp_d;
  real r2;
  int is_short=0;
  int col1, col2, inc=ntypes*ntypes;
  int q_typ, p_typ;
  cell *p, *q;
  real tmp_sprod, tmp_boxl;
  real na_pot_p, na_pot_q, cr_pot;
  real z_sm_p, z_sm_q;

  /* Initialization of all cells and compuation of diagonal part */
  for(r=0; r<ncells; ++r)
    {
      p = CELLPTR(r);
      for (i=0; i<p->n; ++i) 
	{ 
	  p_typ   = SORTE(p,i);
	  CHI_SM(p,i) = sm_chi_0[p_typ];
	}
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
	    z_sm_p = sm_Z[p_typ];
	    
	    /* For each atom in second cell */
	    jstart = (p==q ? i+1 : 0);
	    
	    for (j=jstart; j<q->n; ++j) {
	      
	      /* calculate distance */ 
	      d.x = ORT(p,i,X) - ORT(q,j,X);
	      d.y = ORT(p,i,Y) - ORT(q,j,Y);
	      d.z = ORT(p,i,Z) - ORT(q,j,Z);

	      /* Apply the minimum image convention ? */
	      tmp_sprod = SPROD(d,box_x);
	      tmp_boxl  = 0.5 * SPROD(box_x,box_x);
	      if ( tmp_sprod > tmp_boxl )
		{
		  d.x -= box_x.x;
		  d.y -= box_x.y;
		  d.z -= box_x.z;
		}
	      if ( tmp_sprod < -tmp_boxl ) 
		{
		  d.x += box_x.x;
		  d.y += box_x.y;
		  d.z += box_x.z;
		}
	      
	      tmp_sprod = SPROD(d,box_y);
	      tmp_boxl  = 0.5 * SPROD(box_y,box_y);
	      if ( tmp_sprod > tmp_boxl ) 
		{
		  d.x -= box_y.x;
		  d.y -= box_y.y;
		  d.z -= box_y.z;
		}
	      if ( tmp_sprod < -tmp_boxl ) 
		{
		  d.x += box_y.x;
		  d.y += box_y.y;
		  d.z += box_y.z;
		}
	      
	      tmp_sprod = SPROD(d,box_z);
	      tmp_boxl  = 0.5 * SPROD(box_z,box_z);
	      if ( tmp_sprod > tmp_boxl ) 
		{
		  d.x -= box_z.x;
		  d.y -= box_z.y;
		  d.z -= box_z.z;
		}
	      if ( tmp_sprod < -tmp_boxl ) 
		{
		  d.x += box_z.x;
		  d.y += box_z.y;
		  d.z += box_z.z;
		}
	      
	      q_typ = SORTE(q,j);
	      z_sm_q = sm_Z[q_typ];

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
	      if (r2 < na_pot_tab.end[col2]){  /* AABB */
		VAL_FUNC(na_pot_p, na_pot_tab, col2, inc, r2, is_short);
	      }
	      if (r2 < na_pot_tab.end[col1]){
		VAL_FUNC(na_pot_q, na_pot_tab, col1, inc, r2, is_short);
	      }
	      if (r2 < cr_pot_tab.end[col2]){  /* ABBC */
		VAL_FUNC(cr_pot, cr_pot_tab, col2, inc, r2, is_short);
	      }
	      CHI_SM(p,i) += z_sm_q*(na_pot_p-cr_pot)*coul_eng;
	      CHI_SM(q,j) += z_sm_p*(na_pot_q-cr_pot)*coul_eng;
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

#ifdef DEBUG
  printf("do_v_real\n");
#endif

  int i, j, k, r, s, jstart;
  vektor d, tmp_d;
  real r2;
  int is_short=0;
  int col1, col2, inc=ntypes*ntypes;
  int q_typ, p_typ;
  cell *p, *q;
  real tmp_sprod, tmp_boxl;
  real erfc_r=0, cr_pot=0;
  real j_sm_p, j_sm_q, v_sm_p, v_sm_q;

  /* Inititalisation of all cells and computation of diagonal part */
  for(r=0; r<ncells; ++r)
    {
      p = CELLPTR(r);
      for (i=0; i<p->n; ++i) 
	{
	  p_typ   = SORTE(p,i);
	  j_sm_p  = sm_J_0[p_typ];
	  /* definition */
	  V_SM(p,i) = Q_SM(p,i)*(j_sm_p-2*ew_vorf*coul_eng);
	}
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

	    /* For each atom in second cell */
	    jstart = (p==q ? i+1 : 0);
	    
	    /* for each atom in neighbouring cell */
	    for (j=jstart; j<q->n; ++j) {
	      
	      /* calculate distance */ 
	      d.x = ORT(p,i,X) - ORT(q,j,X);
	      d.y = ORT(p,i,Y) - ORT(q,j,Y);
	      d.z = ORT(p,i,Z) - ORT(q,j,Z);

	      /* Apply the minimum image convention ? */
	      tmp_sprod = SPROD(d,box_x);
	      tmp_boxl  = 0.5 * SPROD(box_x,box_x);
	      if ( tmp_sprod > tmp_boxl )
		{
		  d.x -= box_x.x;
		  d.y -= box_x.y;
		  d.z -= box_x.z;
		}
	      if ( tmp_sprod < -tmp_boxl ) 
		{
		  d.x += box_x.x;
		  d.y += box_x.y;
		  d.z += box_x.z;
		}
	      
	      tmp_sprod = SPROD(d,box_y);
	      tmp_boxl  = 0.5 * SPROD(box_y,box_y);
	      if ( tmp_sprod > tmp_boxl ) 
		{
		  d.x -= box_y.x;
		  d.y -= box_y.y;
		  d.z -= box_y.z;
		}
	      if ( tmp_sprod < -tmp_boxl ) 
		{
		  d.x += box_y.x;
		  d.y += box_y.y;
		  d.z += box_y.z;
		}
	      
	      tmp_sprod = SPROD(d,box_z);
	      tmp_boxl  = 0.5 * SPROD(box_z,box_z);
	      if ( tmp_sprod > tmp_boxl ) 
		{
		  d.x -= box_z.x;
		  d.y -= box_z.y;
		  d.z -= box_z.z;
		}
	      if ( tmp_sprod < -tmp_boxl ) 
		{
		  d.x += box_z.x;
		  d.y += box_z.y;
		  d.z += box_z.z;
		}

	      r2    = SPROD(d,d);
	      /* redundant due to symmetry p<->q
		 col1  = q_typ * ntypes + p_typ; */
	      col2  = p_typ * ntypes + q_typ;
	      
#ifdef DEBUG
	      if (0==r2) { char msgbuf[256];
		sprintf(msgbuf, "Distance 2 is zero in module between particles sm %d and %d!\n",
			NUMMER(p,i), NUMMER(q,j));
		error(msgbuf);
	      }
#endif
	      
	      /* compute real space term of v_i */
	      if (r2 < erfc_r_tab.end[col2]){ /* AAAA */
		VAL_FUNC(erfc_r, erfc_r_tab, col2, inc, r2, is_short);
	      }
	      if (r2 < cr_pot_tab.end[col2]){ /* ABBC */
		VAL_FUNC(cr_pot, cr_pot_tab, col2, inc, r2, is_short);
	      }
	      /* definition */
	      V_SM(p,i) += Q_SM(q,j)*(erfc_r+cr_pot)*coul_eng;
	      V_SM(q,j) += Q_SM(p,j)*(erfc_r+cr_pot)*coul_eng;
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
#ifdef DEBUG
  printf("do_v_kspace\n");
#endif
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
  real totpot, tmp;

  /* if (myid==0) printf("start do_cg\n"); */

  /* update V_SM */
#ifdef NBLIST
  calc_sm_pot();
#else
  do_v_real();
  do_v_kspace();
#endif

  norm_b     = 0.0;
  totpot     = 0.0;
  rho[kstep] = 0.0;
  for (k=0; k<ncells; ++k) {
    int  i;
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      /* initial values */
      X_SM(p,i)   = Q_SM(p,i);
      R_SM(p,i)   = B_SM(p,i)-V_SM(p,i);
      totpot     += V_SM(p,i);
      norm_b     += B_SM(p,i)*B_SM(p,i);
      rho[kstep] += R_SM(p,i)*R_SM(p,i);
    }
  }
#ifdef MPI
  MPI_Allreduce(&norm_b,   &tmp, 1, REAL, MPI_SUM, cpugrid); norm_b=tmp;
  MPI_Allreduce(&totpot,   &tmp, 1, REAL, MPI_SUM, cpugrid); totpot=tmp;
  MPI_Allreduce(rho+kstep, &tmp, 1, REAL, MPI_SUM, cpugrid); rho[kstep]=tmp;
#endif
  
  tolerance = epsilon_cg*SQRT(norm_b);
  tolerance2 = SQR(tolerance);

#ifdef DEBUG
  printf("tolerance after kstep %d: %e, %e, totpot: %e\n", 
         kstep,tolerance,SQRT(rho[kstep])/SQRT(norm_b),totpot);
#endif
  /*
  printf("a rank %d, kstep %d, tol %f\n", myid,kstep,rho[kstep]);
  */
  while ((rho[kstep] > tolerance2) && (kstep < kstepmax)) {

    kstep++;

    beta = (kstep == 1) ? 0.0 : rho[kstep-1]/rho[kstep-2];
    for (k=0; k<ncells; ++k) {
      int  i;
      cell *p = CELLPTR(k);
      for (i=0; i<p->n; ++i) {
        Q_SM(p,i) = D_SM(p,i) = R_SM(p,i) + beta * D_SM(p,i);
      }
    }

    /* update V_SM */
#ifdef NBLIST
    calc_sm_pot();
#else
    do_v_real();
    do_v_kspace();
#endif

    dad = 0.0;
    totpot = 0.0;
    for (k=0; k<ncells; ++k) {
      int  i;
      cell *p = CELLPTR(k);
      for (i=0; i<p->n; ++i) {
        dad    += D_SM(p,i)*V_SM(p,i);
        totpot += V_SM(p,i);
      }
    }
#ifdef MPI
    MPI_Allreduce(&dad,    &tmp, 1, REAL, MPI_SUM, cpugrid); dad   =tmp;
    MPI_Allreduce(&totpot, &tmp, 1, REAL, MPI_SUM, cpugrid); totpot=tmp;
#endif

    alpha = rho[kstep-1]/dad;
    rho[kstep] = 0.0;
    for (k=0; k<ncells; ++k) {
      int  i;
      cell *p = CELLPTR(k);
      for (i=0; i<p->n; ++i) {      
        X_SM(p,i)  += alpha*D_SM(p,i);
        R_SM(p,i)  -= alpha*V_SM(p,i);
        rho[kstep] += R_SM(p,i)*R_SM(p,i);
      }
    }
#ifdef MPI
    MPI_Allreduce(rho+kstep, &tmp, 1, REAL, MPI_SUM, cpugrid); rho[kstep]=tmp;
#endif
      
#ifdef DEBUG
    printf("tolerance after kstep %d: %e, %e, totpot: %e\n", 
           kstep,tolerance,SQRT(rho[kstep])/SQRT(norm_b),totpot);
#endif
    /*
    printf("rank %d, kstep %d, tol %f\n", myid,kstep,rho[kstep]);
    */
  }
}

/*****************************************************************************
*
* Charge update
*
******************************************************************************/

void do_charge_update(void)
{
#ifdef DEBUG
  printf("do_charge_update\n");
#endif
  
  int k, typ;
  real sum1, sum2, potchem;
  real q_Al, q_O, q_tot, tmp;
  
  /* Update electronegativity since coordinates have changed */
#ifdef NBLIST
  calc_sm_chi();
#else
  do_electronegativity();
#endif

  /* Solving the first linear system V_ij s_j = -chi_i */
  
  for (k=0; k<ncells; ++k) {
    int  i;
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      B_SM(p,i) = -CHI_SM(p,i);
      /* Initial value of the charges */
      Q_SM(p,i) = CHARGE(p,i);
    }
  }
  
#ifdef DEBUG
  printf("do_cg %d\n",1);
#endif
  do_cg();
  
  /* Sum up for getting charges */
  sum1=0.0;
  for (k=0; k<ncells; ++k) {
    int  i;
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      sum1     += X_SM(p,i);
      S_SM(p,i) = X_SM(p,i);
    }
  }
  
  /* Solving the second linear system V_ij t_j = -1 */

  for (k=0; k<ncells; ++k) {
    int  i;
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      B_SM(p,i) = -1;
      /* Initial value of the charges */
      Q_SM(p,i) = 0.0;
    }
  }

#ifdef DEBUG
  printf("do_cg %d\n",2);
#endif
  do_cg();
  
  /* Sum up for getting charges */
  sum2=0.0;
    for (k=0; k<ncells; ++k) {
    int  i;
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      sum2 += X_SM(p,i);
    }
  }
  
  /* Calculate chemical potential */
  potchem = sum1/sum2;
      
  /* Calculate new charges of the atoms */
  q_tot = 0.0; q_Al = 0.0; q_O = 0.0;
  for (k=0; k<ncells; ++k) {
    int  i;
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      CHARGE(p,i) = S_SM(p,i)-potchem*X_SM(p,i);
      q_tot += CHARGE(p,i); 
      typ = SORTE(p,i);
      if (typ == 0) {
	q_Al += CHARGE(p,i);
      }
      else {
	q_O += CHARGE(p,i);
      }
    }
  }
#ifdef MPI
  MPI_Allreduce(&sum1,  &tmp,  1, REAL,    MPI_SUM, cpugrid); sum1 =tmp;
  MPI_Allreduce(&sum2,  &tmp,  1, REAL,    MPI_SUM, cpugrid); sum2 =tmp;
  MPI_Allreduce(&q_tot, &tmp,  1, REAL,    MPI_SUM, cpugrid); q_tot=tmp;
  MPI_Allreduce(&q_Al,  &tmp,  1, REAL,    MPI_SUM, cpugrid); q_Al =tmp;
  MPI_Allreduce(&q_O,   &tmp,  1, REAL,    MPI_SUM, cpugrid); q_O  =tmp;
#endif
  
  if (0==myid) {
    printf("Sums: sum1 = %e sum2 = %e\n", sum1, sum2);
    printf("Total charge: qtot = %e\n", q_tot);
#ifndef DEBUG
    printf("Average charge of Al: qAl = %e\n", q_Al/num_sort[0]);
    printf("Average charge of O: qO = %e\n",   q_O /num_sort[1]);
#else
    printf("Total charge and number of Al: qAl = %e %d %e\n", 
           q_Al, num_sort[0], q_Al/num_sort[0]);
    printf("Total charge and number of O: qO = %e %d %e\n", 
           q_O, num_sort[1], q_O/num_sort[1]);
#endif
  }
}

/*****************************************************************************
*
* Charge update a la Elsener et al. (Mod. Sim. Mat. Sci. 16, 0250006 (2008))
*
*   We solve the linear system
*
*      [ V 1 ] [ Q_SM ] = [ -CHI_SM ]
*      [ 1 0 ] [  Q_0 ]   [    0    ]
*
*   The matrix-vector product V*Q_SM ist stored in V_SM (by calc_sm_pot). 
*   Q_SM stores the subsequent charge corrections (with Q_0 the negative
*   of the chemical potential), R_SM the residuals of the system above.
*
******************************************************************************/

void charge_update_sm(void) {

  real tmpvec1[3], tmpvec2[3], *tmpvec;
  real r_old, r_new, alpha, beta;
  real Q_0, V_0, R_0, sm_tol = 1e-5;
  int  i, k, itr=0, sm_max_itr = 10;

#ifdef MPI
  tmpvec = tmpvec2;
#else
  tmpvec = tmpvec1;
#endif

  /* assign initial charges */
  for (k=0; k<ncells; ++k) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      Q_SM(p,i) = CHARGE(p,i);
    }
  }
#ifdef NBLIST
  calc_sm_chi();
  calc_sm_pot();
#else
  do_electronegativity();
  do_v_real();
  do_v_kspace();
#endif

  /* first residuals and first charge correction; we choose the 
     initial chemical potential such that the residuals are minimal */
  tmpvec1[0] = tmpvec1[1] = tmpvec1[2] = 0.0; 
  for (k=0; k<ncells; ++k) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      Q_SM(p,i) = R_SM(p,i) = -CHI_SM(p,i) - V_SM(p,i);
      tmpvec1[0] += CHARGE(p,i);
      tmpvec1[1] += R_SM(p,i);
      tmpvec1[2] += SQR(R_SM(p,i));
    }
  }
#ifdef MPI
  MPI_Allreduce(tmpvec1, tmpvec2, 3, REAL, MPI_SUM, cpugrid);
#endif
  R_0 = Q_0 = -tmpvec[0];
  tmpvec[1] /= natoms;
  r_new = r_old = (tmpvec[2] + SQR(R_0)) / natoms - SQR(tmpvec[1]);
  for (k=0; k<ncells; ++k) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      Q_SM(p,i) -= tmpvec[1];
      R_SM(p,i) -= tmpvec[1];
    }
  }

  if (myid==0) printf("itr: %d, r_new: %f\n", itr, r_new);

  /* now the iteration starts ... */
  while ((r_new > sm_tol) && (itr++ < sm_max_itr)) {

#ifdef NBLIST
    calc_sm_pot();
#else
    do_v_real();
    do_v_kspace();
#endif

    /* reduction loop for size of correction*/
    tmpvec1[0] = tmpvec1[1] = 0.0; 
    for (k=0; k<ncells; ++k) {
      cell *p = CELLPTR(k);
      for (i=0; i<p->n; ++i) {
        V_SM(p,i)  += Q_0;
        tmpvec1[0] += Q_SM(p,i);
        tmpvec1[1] += Q_SM(p,i) * V_SM(p,i);
      }
    }
#ifdef MPI
    MPI_Allreduce(tmpvec1, tmpvec2, 2, REAL, MPI_SUM, cpugrid);
#endif
    V_0 = tmpvec[0];
    alpha = natoms * r_old / (tmpvec[1] + V_0 * Q_0);

    /* new residuals, corrected charge */
    tmpvec1[0] = 0.0;
    for (k=0; k<ncells; ++k) {
      cell *p = CELLPTR(k);
      for (i=0; i<p->n; ++i) {
        CHARGE(p,i) += alpha * Q_SM(p,i);
        R_SM(p,i)   -= alpha * V_SM(p,i);
        tmpvec1[0]  += SQR(R_SM(p,i));
      }
    }
#ifdef MPI
    MPI_Allreduce(tmpvec1, tmpvec2, 1, REAL, MPI_SUM, cpugrid);
#endif
    R_0  -= alpha * V_0;
    r_new = (tmpvec[0] + SQR(R_0)) / natoms;

    if (myid==0) printf("itr: %d, r_new: %f\n", itr, r_new);

    /* stop if already close enough */
    if ((r_new < sm_tol) || (itr >= sm_max_itr)) break;

    beta  = r_new / r_old; 
    r_old = r_new;

    /* preparation vor the next charge correction */
    for (k=0; k<ncells; ++k) {
      cell *p = CELLPTR(k);
      for (i=0; i<p->n; ++i) {
        Q_SM(p,i)  = R_SM(p,i) + beta * Q_SM(p,i);
      }
    }
    Q_0 = R_0 + beta * Q_0;
  }

  /* print average charges for each atom type */
  tmpvec1[0] = tmpvec1[1] = 0.0;
  for (k=0; k<ncells; ++k) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      tmpvec1[SORTE(p,i)] += CHARGE(p,i);
    }
  }
#ifdef MPI
  MPI_Allreduce( tmpvec1, tmpvec2, 2, REAL, MPI_SUM, cpugrid);
#endif
  if (0==myid) {
    printf("Total charge:        qtot = %e\n", (tmpvec[0]+tmpvec[1])/natoms);
    printf("Average charge of Al: qAl = %e\n", tmpvec[0] / num_sort[0]);
    printf("Average charge of O:   qO = %e\n", tmpvec[1] / num_sort[1]);
  }

}

