
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2011 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_forces_nbl.c -- force loop with neighbor lists
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"
#include "potaccess.h"

/*****************************************************************************

  Neighbor list format

  Normally, each atom is identified by its cell number k, and the number
  i within that cell. Here, we want to identify it with a single number n. 
  To translate between the two labellings, two auxiliary arrays cl_off
  and cl_num are used, such that the following relations hold:
  n = cl_off[k] + i, k = cl_num[n]. Note that this numbering also
  includes buffer cell atoms. It runs consecutively over the atoms
  in cell 0, cell 1, etc.
  
  The neighbors of atom i are stored in the array tb, in the index range
  tl[i] .. tl[i+1]-1. Unlike the numbers of neighbor atoms, i enumerates 
  only the real atoms, not buffer atoms, and refers to a different numbering 
  scheme. It runs first over the atoms in the first inner cell, then those 
  in the second inner cell, etc.

******************************************************************************/

#define NBLMINLEN 100000

int  *tl=NULL, *tb=NULL, *cl_off=NULL, *cl_num=NULL, nb_max=0;


/******************************************************************************
*
*  deallocate (largest part of) neighbor list
*
******************************************************************************/

void deallocate_nblist(void)
{
#if defined(DEBUG) || defined(TIMING)
  if (myid==0)
    printf("Size of neighbor table: %d MB\n", 
           nb_max * sizeof(int) / SQR(1024) );
#endif
  if (tb) free(tb);
  tb = NULL;
  have_valid_nbl = 0;
}

/******************************************************************************
*
*  estimate_nblist_size
*
******************************************************************************/

int estimate_nblist_size(void)
{
  int  c, tn=1;

  /* for all cells */
  for (c=0; c<ncells2; c++) {

    int i, c1 = cnbrs[c].np;
    cell *p   = cell_array + c1;

    /* for each atom in cell */
    for (i=0; i<p->n; i++) {

      int    m;
      vektor d1;

      d1.x = ORT(p,i,X);
      d1.y = ORT(p,i,Y);
#ifndef TWOD
      d1.z = ORT(p,i,Z);
#endif

      /* for each neighboring atom */
      for (m=0; m<14; m++) {   /* this is not TWOD ready! */

        int    c2, jstart, j;
        real   r2;
        cell   *q;

        c2 = cnbrs[c].nq[m];
        if (c2<0) continue;
        if (c2==c1) jstart = i+1;
        else        jstart = 0;
        q = cell_array + c2;
#ifdef ia64
#pragma ivdep
#endif
        for (j=jstart; j<q->n; j++) {
          vektor d;
          d.x = ORT(q,j,X) - d1.x;
          d.y = ORT(q,j,Y) - d1.y;
#ifndef TWOD
          d.z = ORT(q,j,Z) - d1.z;
#endif
          r2  = SPROD(d,d);
          if (r2 < cellsz) tn++;
        }
      }
    }
  }
#ifdef MPI
  /* printf ("myid: %d nb-list size: %d\n",myid,tn);fflush(stdout); */
#endif
  return tn;
}

/******************************************************************************
*
*  make_nblist
*
******************************************************************************/

void make_nblist(void)
{
  static int at_max=0, pa_max=0, ncell_max=0;
  int  c, i, k, n, tn, at, cc;

  /* update reference positions */
  for (k=0; k<ncells; k++) {
    cell *p = cell_array + cnbrs[k].np;
#ifdef ia64
#pragma ivdep,swp
#endif
    for (i=0; i<p->n; i++) {
      NBL_POS(p,i,X) = ORT(p,i,X);
      NBL_POS(p,i,Y) = ORT(p,i,Y);
#ifndef TWOD
      NBL_POS(p,i,Z) = ORT(p,i,Z);
#endif
    }
  }

  /* (re)allocate cl_off */
  if (nallcells > ncell_max) {
    cl_off = (int *) realloc( cl_off, nallcells * sizeof(int) );
    if (cl_off==NULL) error("cannot allocate neighbor table");
    ncell_max = nallcells;
  }

  /* count atom numbers (including buffer atoms) */
  at=0;
  for (k=0; k<nallcells; k++) {
    cell *p = cell_array + k;
    cl_off[k] = at;
    at += p->n;
  }

  /* (re-)allocate neighbor table */
  if (at >= at_max) {
    free(tl);
    free(cl_num);
    at_max = (int) (1.1 * at);
    tl     = (int *) malloc(at_max * sizeof(int));
    cl_num = (int *) malloc(at_max * sizeof(int));
  }
  if (NULL==tb) {
    if (0==last_nbl_len) 
      nb_max = MAX(((int) (nbl_size * estimate_nblist_size())), (NBLMINLEN));
    else
      nb_max = (int) (nbl_size * last_nbl_len);
    tb = (int *) malloc(nb_max * sizeof(int));
  }
  else if (last_nbl_len * sqrt(nbl_size) > nb_max) {
    free(tb);
    nb_max = (int  ) (nbl_size * last_nbl_len);
    tb     = (int *) malloc(nb_max * sizeof(int));
  }
  if ((tl==NULL) || (tb==NULL) || (cl_num==NULL)) 
    error("cannot allocate neighbor table");

  /* set cl_num */
  n=0;
  for (k=0; k<nallcells; k++) {
    cell *p = cell_array + k;
    for (i=0; i<p->n; i++) cl_num[n++] = k;
  }

  /* for all cells */
  n=0; tn=0; tl[0]=0;
  for (c=0; c<ncells2; c++) {

    int  c1 = cnbrs[c].np;
    cell *p = cell_array + c1;

    /* for each atom in cell */
    for (i=0; i<p->n; i++) {

      int    m;
      vektor d1;

      d1.x = ORT(p,i,X);
      d1.y = ORT(p,i,Y);
#ifndef TWOD
      d1.z = ORT(p,i,Z);
#endif

      /* for each neighboring atom */
      for (m=0; m<14; m++) {   /* this is not TWOD ready! */
        int  c2, jstart, j;
        cell *q;
        c2 = cnbrs[c].nq[m];
        if (c2<0) continue;
        if (c2==c1) jstart = i+1;
        else        jstart = 0;
        q = cell_array + c2;
#ifdef ia64
#pragma ivdep
#endif
        for (j=jstart; j<q->n; j++) {
          vektor d;
          real   r2;
          d.x = ORT(q,j,X) - d1.x;
          d.y = ORT(q,j,Y) - d1.y;
#ifndef TWOD
          d.z = ORT(q,j,Z) - d1.z;
#endif
          r2  = SPROD(d,d);
          if (r2 < cellsz) {
            tb[tn++] = cl_off[c2] + j;
          }
        }
      }
      pa_max = MAX(pa_max,tn-tl[n]);
      tl[++n] = tn;
      if (tn > nb_max-2*pa_max) {
        error("neighbor table full - increase nbl_size");
      }
    }
  }
  last_nbl_len   = tn;
  have_valid_nbl = 1;
  nbl_count++;
}

/******************************************************************************
*
*  calc_forces
*
******************************************************************************/

void calc_forces(int steps)
{
  int  i, b, k, n=0, is_short=0, idummy=0;
  real tmpvec1[8], tmpvec2[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

#ifdef DIPOLE
  static int dp_E_calc=0; 	/* Number of field iterations */
  int dp_it=0;			/* Number of dipole iterations */
  int dp_p_calc=0;		/* Calculate dipoles or keep them */
  /* TODO: Communicate! */
  int dp_converged=0;
  real dp_sum_old, dp_sum=1.;
  real max_diff=10.;
  real *dp_E_shift;
  dp_p_calc = ((dp_fix-1 + dp_fix*dp_E_calc)>0 ) ? 0 : 1;
#endif

  if (0==have_valid_nbl) {
#ifdef MPI
    /* check message buffer size */
    if (0 == nbl_count % BUFSTEP) setup_buffers();
#endif
    /* update cell decomposition */
    fix_cells();
  }

  /* fill the buffer cells */
  send_cells(copy_cell,pack_cell,unpack_cell);

  /* make new neighbor lists */
  if (0==have_valid_nbl) make_nblist();

  /* clear global accumulation variables */
  tot_pot_energy = 0.0;
#ifdef SM
  tot_sm_es_energy = 0.0;
#endif
  virial = 0.0;
  vir_xx = 0.0;
  vir_yy = 0.0;
  vir_xy = 0.0;
  vir_zz = 0.0;
  vir_yz = 0.0;
  vir_zx = 0.0;
  nfc++;

  /* clear per atom accumulation variables, also in buffer cells */
  for (k=0; k<nallcells; k++) {
    cell *p = cell_array + k;
#ifdef ia64
#pragma ivdep,swp
#endif
    for (i=0; i<p->n; i++) {
      KRAFT(p,i,X) = 0.0;
      KRAFT(p,i,Y) = 0.0;
#ifndef TWOD
      KRAFT(p,i,Z) = 0.0;
#endif
#if defined(STRESS_TENS)
      PRESSTENS(p,i,xx) = 0.0;
      PRESSTENS(p,i,yy) = 0.0;
      PRESSTENS(p,i,xy) = 0.0;
#ifndef TWOD
      PRESSTENS(p,i,zz) = 0.0;
      PRESSTENS(p,i,yz) = 0.0;
      PRESSTENS(p,i,zx) = 0.0;
#endif
#endif     
#ifndef MONOLJ
      POTENG(p,i) = 0.0;
#endif
#ifdef NNBR
      NBANZ(p,i) = 0;
#endif
#ifdef CNA
      if (cna) MARK(p,i) = 0;
#endif
#ifdef COVALENT
      NEIGH(p,i)->n = 0;
#endif
#ifdef EAM2
      EAM_RHO(p,i) = 0.0;
#ifdef EEAM
      EAM_P(p,i) = 0.0;
#endif
#endif
#ifdef ADP
      ADP_MU    (p,i,X)  = 0.0;
      ADP_MU    (p,i,Y)  = 0.0;
      ADP_MU    (p,i,Z)  = 0.0;
      ADP_LAMBDA(p,i,xx) = 0.0;
      ADP_LAMBDA(p,i,yy) = 0.0;
      ADP_LAMBDA(p,i,xy) = 0.0;
      ADP_LAMBDA(p,i,zz) = 0.0;
      ADP_LAMBDA(p,i,yz) = 0.0;
      ADP_LAMBDA(p,i,zx) = 0.0;
#endif
#ifdef DIPOLE
      DP_E_STAT(p,i,X) = 0.0;
      DP_E_STAT(p,i,Y) = 0.0;
      DP_E_STAT(p,i,Z) = 0.0;
      DP_P_STAT(p,i,X) = 0.0;
      DP_P_STAT(p,i,Y) = 0.0;
      DP_P_STAT(p,i,Z) = 0.0;
      DP_E_IND(p,i,X)  = 0.0;
      DP_E_IND(p,i,Y)  = 0.0;
      DP_E_IND(p,i,Z)  = 0.0;
      if ( dp_p_calc ) {
	DP_P_IND(p,i,X)  = 0.0;
	DP_P_IND(p,i,Y)  = 0.0;
	DP_P_IND(p,i,Z)  = 0.0;
      }
#endif /* dipole */
    }
  }

  /* clear total forces */
#ifdef RIGID
  if ( nsuperatoms>0 ) 
    for(i=0; i<nsuperatoms; i++) {
      superforce[i].x = 0.0;
      superforce[i].y = 0.0;
#ifndef TWOD
      superforce[i].z = 0.0;
#endif
    }
#endif

#ifdef EWALD
  if (steps==0) {
    ewald_time.total = 0.0;
    imd_start_timer( &ewald_time );
  }
#endif

  /* pair interactions - for all atoms */
  n=0;
  for (k=0; k<ncells; k++) {
    cell *p = cell_array + cnbrs[k].np;
    for (i=0; i<p->n; i++) {

#ifdef STRESS_TENS
#ifdef TWOD
      sym_tensor pp = {0.0,0.0,0.0};
#else
      sym_tensor pp = {0.0,0.0,0.0,0.0,0.0,0.0};
#endif
#endif
#ifdef ADP
      real       tmp;
      vektor     mu = {0.0,0.0,0.0};
      sym_tensor la = {0.0,0.0,0.0,0.0,0.0,0.0};
#endif
#ifdef COULOMB
      real       phi, grphi, chg;
#endif
#ifdef DIPOLE
      real       tmp;
      vektor     Estat = {0.0,0.0,0.0};
      vektor     pstat = {0.0,0.0,0.0};
#endif
#ifdef TWOD
      vektor d1, ff = {0.0,0.0};
#else
      vektor d1, ff = {0.0,0.0,0.0};
#endif
      real   ee = 0.0;
      real   eam_r = 0.0, eam_p = 0.0;
      int    m, it, nb = 0;

      d1.x = ORT(p,i,X);
      d1.y = ORT(p,i,Y);
#ifndef TWOD
      d1.z = ORT(p,i,Z);
#endif
      it   = SORTE(p,i);

      /* loop over neighbors */
#ifdef ia64
#pragma ivdep
#endif
      for (m=tl[n]; m<tl[n+1]; m++) {

        vektor d, force;
        cell   *q;
        real   pot, grad, r2, rho_h;
        int    c, j, jt, col, col2, inc = ntypes * ntypes;

        c = cl_num[ tb[m] ];
        j = tb[m] - cl_off[c];
        q = cell_array + c;

        d.x = ORT(q,j,X) - d1.x;
        d.y = ORT(q,j,Y) - d1.y;
#ifndef TWOD
        d.z = ORT(q,j,Z) - d1.z;
#endif
        r2  = SPROD(d,d);
        jt  = SORTE(q,j);
        col = it * ntypes + jt;
        col2= jt * ntypes + it;

        /* compute pair interactions */
#if defined(PAIR) || defined(KEATING)
        /* PAIR and KEATING are mutually exclusive */
#if defined(PAIR)
        if (r2 <= pair_pot.end[col])
#elif defined(KEATING)
        if (r2 < keat_r2_cut[it][jt]) 
#endif
	{
#if defined(PAIR)
#ifdef LINPOT
          PAIR_INT_LIN(pot, grad, pair_pot_lin, col, inc, r2, is_short);
#else
	  PAIR_INT(pot, grad, pair_pot, col, inc, r2, is_short);
#endif
#elif defined(KEATING)
          PAIR_INT_KEATING(pot, grad, it, jt, r2);
#endif

          tot_pot_energy += pot;
          force.x = d.x * grad;
          force.y = d.y * grad;
#ifndef TWOD
          force.z = d.z * grad;
#endif
          KRAFT(q,j,X) -= force.x;
          KRAFT(q,j,Y) -= force.y;
#ifndef TWOD
          KRAFT(q,j,Z) -= force.z;
#endif
          ff.x         += force.x;
          ff.y         += force.y;
#ifndef TWOD
          ff.z         += force.z;
#endif

#ifdef FLAGEDATOMS
	  if(VSORTE(q,j) == flagedatomstype && VSORTE(p,i) == flagedatomstype)
	    {
	      //	      printf("Atom nr %d of type %d interacting with %d: Pair forces : %e %e %e\n",
	      printf("%d %d %d %e %e %e %e %e %e\n",
		     NUMMER(p,i),VSORTE(p,i),NUMMER(q,j),d.x, d.y, d.z, force.x,force.y,force.z);
	      fflush(stdout);
	    }
#endif

#ifndef MONOLJ
          pot *= 0.5;   /* avoid double counting */
#ifdef NNBR
          if (r2 < nb_r2_cut[col ]) nb++;
          if (r2 < nb_r2_cut[col2]) NBANZ(q,j)++;
#endif
#ifdef ORDPAR
          if (r2 < op_r2_cut[col ]) ee          += op_weight[col ] * pot;
          if (r2 < op_r2_cut[col2]) POTENG(q,j) += op_weight[col2] * pot;
#else
          ee          += pot;
          POTENG(q,j) += pot;
#endif
#endif
#ifdef P_AXIAL
          vir_xx -= d.x * force.x;
          vir_yy -= d.y * force.y;
#ifndef TWOD
          vir_zz -= d.z * force.z;
#endif
#else
          virial -= r2  * grad;
#endif

#ifdef STRESS_TENS
          if (do_press_calc) {
            /* avoid double counting of the virial */
            force.x *= 0.5;
            force.y *= 0.5;
#ifndef TWOD
            force.z *= 0.5;
#endif
            pp.xx             -= d.x * force.x;
            PRESSTENS(q,j,xx) -= d.x * force.x;
            pp.yy             -= d.y * force.y;
            PRESSTENS(q,j,yy) -= d.y * force.y;
            pp.xy             -= d.x * force.y;
            PRESSTENS(q,j,xy) -= d.x * force.y;
#ifndef TWOD
            pp.zz             -= d.z * force.z;
            PRESSTENS(q,j,zz) -= d.z * force.z;
            pp.yz             -= d.y * force.z;
            PRESSTENS(q,j,yz) -= d.y * force.z;
            pp.zx             -= d.z * force.x;
            PRESSTENS(q,j,zx) -= d.z * force.x;
#endif
	  }
#endif
        }

#endif /* PAIR || KEATING */

#ifdef EAM2
        /* compute host electron density */
        if (r2 < rho_h_tab.end[col])  {
          VAL_FUNC(rho_h, rho_h_tab, col, inc, r2, is_short);
          eam_r += rho_h;
#ifdef EEAM
          eam_p += rho_h*rho_h; 
#endif
        }
        if (it==jt) {
          if (r2 < rho_h_tab.end[col]) {
            EAM_RHO(q,j) += rho_h;
#ifdef EEAM
            EAM_P(q,j) += rho_h*rho_h;
#endif
          } 
        } else {
          if (r2 < rho_h_tab.end[col2]) {
            VAL_FUNC(rho_h, rho_h_tab, col2, inc, r2, is_short);
            EAM_RHO(q,j) += rho_h; 
#ifdef EEAM
            EAM_P(q,j) += rho_h*rho_h; 
#endif
          }
        }
#endif

#ifdef ADP
        /* compute adp_mu */
        if (r2 < adp_upot.end[col])  {
          VAL_FUNC(pot, adp_upot, col, inc, r2, is_short);
          tmp = pot * d.x;  mu.x += tmp;  ADP_MU(q,j,X) -= tmp;
          tmp = pot * d.y;  mu.y += tmp;  ADP_MU(q,j,Y) -= tmp;
          tmp = pot * d.z;  mu.z += tmp;  ADP_MU(q,j,Z) -= tmp;
        }
        /* compute adp_lambda */
        if (r2 < adp_wpot.end[col])  {
          VAL_FUNC(pot, adp_wpot, col, inc, r2, is_short);
          tmp = pot * d.x * d.x;  la.xx += tmp;  ADP_LAMBDA(q,j,xx) += tmp;
          tmp = pot * d.y * d.y;  la.yy += tmp;  ADP_LAMBDA(q,j,yy) += tmp;
          tmp = pot * d.z * d.z;  la.zz += tmp;  ADP_LAMBDA(q,j,zz) += tmp;
          tmp = pot * d.y * d.z;  la.yz += tmp;  ADP_LAMBDA(q,j,yz) += tmp;
          tmp = pot * d.z * d.x;  la.zx += tmp;  ADP_LAMBDA(q,j,zx) += tmp;
          tmp = pot * d.x * d.y;  la.xy += tmp;  ADP_LAMBDA(q,j,xy) += tmp;
        }
#endif /* ADP */

#ifdef COULOMB
#ifdef VARCHG
        chg = CHARGE(p,i) * CHARGE(q,j);
#else
	chg = charge[it]  * charge[jt];
#endif
	if (r2 < ew_r2_cut) {	
	  if (SQR(chg)>0.) {
#ifdef SM
            real cr_pot=0.0, cr_gr=0.0, na_pot_p=0.0, na_pot_q=0.0, na_gr_p=0.0, na_gr_q=0.0, sm_es_energy=0.0;
            real z_sm_p = sm_Z[it] * CHARGE(q,j) * coul_eng;
            real z_sm_q = sm_Z[jt] * CHARGE(p,i) * coul_eng;
#endif
	    /* Constant electric field from charges */
	    /* Coulomb potential is in column 0 */
            int incr = coul_table.ncols;
	    PAIR_INT(phi, grphi, coul_table, 0, incr, r2, is_short);

	    /* Coulomb Energy */
	    pot     = chg * phi;
#ifdef SM
	    sm_es_energy  = chg * phi;
#endif
	    grad    = chg * grphi;
#ifdef SM
            /* Coulomb repulsion potential */
            if (r2 < cr_pot_tab.end[col]) {
              PAIR_INT(cr_pot, cr_gr, cr_pot_tab, col, inc, r2, is_short);
              pot  += cr_pot * (chg * coul_eng - z_sm_p - z_sm_q);
	      sm_es_energy  += cr_pot * (chg * coul_eng - z_sm_p - z_sm_q);
              grad += cr_gr * (chg * coul_eng - z_sm_p - z_sm_q);
            }
            /* nuclear attraction potential */
            if (r2 < na_pot_tab.end[col]) {
              PAIR_INT(na_pot_p, na_gr_p, na_pot_tab, col, inc, r2, is_short);
            }
            if (r2 < na_pot_tab.end[col2]) {
              PAIR_INT(na_pot_q, na_gr_q, na_pot_tab, col2, inc, r2,is_short);
            }
              pot  += z_sm_q * na_pot_p + z_sm_p * na_pot_q;
              sm_es_energy += z_sm_q * na_pot_p + z_sm_p * na_pot_q;
              grad += z_sm_q * na_gr_p + z_sm_p * na_gr_q;

#endif

#ifdef SM
	    tot_sm_es_energy += sm_es_energy;	
#endif

	    tot_pot_energy += pot;
	    force.x = d.x * grad;
	    force.y = d.y * grad;
	    force.z = d.z * grad;

#ifdef EXTF
	    real chg_single;
#ifdef VARCHG
	    chg_single = CHARGE(p,i);
#else
	    chg_single = charge[it];
#endif
	    force.x += chg_single * extf.x; 
	    force.y += chg_single * extf.y; 
	    force.z += chg_single * extf.z; 
#endif /* EXTF */

	    KRAFT(q,j,X) -= force.x;
	    KRAFT(q,j,Y) -= force.y;
	    KRAFT(q,j,Z) -= force.z;
	    ff.x         += force.x;
	    ff.y         += force.y;
	    ff.z         += force.z;
            pot          *= 0.5;   /* avoid double counting */
	    ee           += pot;
	    POTENG(q,j)  += pot;
#ifdef P_AXIAL
	    vir_xx -= d.x * force.x;
	    vir_yy -= d.y * force.y;
	    vir_zz -= d.z * force.z;
#else
	    virial -= r2  * grad;
#endif

#ifdef STRESS_TENS
	    if (do_press_calc) {
	      /* avoid double counting of the virial */
	      force.x *= 0.5;
	      force.y *= 0.5;
	      force.z *= 0.5;
	      pp.xx             -= d.x * force.x;
	      PRESSTENS(q,j,xx) -= d.x * force.x;
	      pp.yy             -= d.y * force.y;
	      PRESSTENS(q,j,yy) -= d.y * force.y;
	      pp.xy             -= d.x * force.y;
	      PRESSTENS(q,j,xy) -= d.x * force.y;
	      pp.zz             -= d.z * force.z;
	      PRESSTENS(q,j,zz) -= d.z * force.z;
	      pp.yz             -= d.y * force.z;
	      PRESSTENS(q,j,yz) -= d.y * force.z;
	      pp.zx             -= d.z * force.x;
	      PRESSTENS(q,j,zx) -= d.z * force.x;
	    }
#endif

#ifdef DIPOLE
#ifdef VARCHG
	    /* Field for Dipole calculation */
	    if (dp_p_calc) {
#ifdef SM
	      Estat.x += d.x * grphi * CHARGE(q,j)
		+ d.x * coul_eng * (CHARGE(q,j)-sm_Z[jt])*na_gr_q;
	      Estat.y += d.y * grphi * CHARGE(q,j)
		+ d.y * coul_eng * (CHARGE(q,j)-sm_Z[jt])*na_gr_q;
	      Estat.z += d.z * grphi * CHARGE(q,j)
		+ d.z * coul_eng * (CHARGE(q,j)-sm_Z[jt])*na_gr_q;
	      DP_E_STAT(q,j,X) -= d.x * grphi * CHARGE(p,i)
		+ d.x * coul_eng * (CHARGE(p,i)-sm_Z[it])*na_gr_p;
	      DP_E_STAT(q,j,Y) -= d.y * grphi * CHARGE(p,i)
		+ d.y * coul_eng * (CHARGE(p,i)-sm_Z[it])*na_gr_p;
	      DP_E_STAT(q,j,Z) -= d.z * grphi * CHARGE(p,i)
		+ d.z * coul_eng * (CHARGE(p,i)-sm_Z[it])*na_gr_p;
#else
	      Estat.x += d.x * grphi * CHARGE(q,j);
	      Estat.y += d.y * grphi * CHARGE(q,j);
	      Estat.z += d.z * grphi * CHARGE(q,j);
	      DP_E_STAT(q,j,X) -= d.x * grphi * CHARGE(p,i);
	      DP_E_STAT(q,j,Y) -= d.y * grphi * CHARGE(p,i);
	      DP_E_STAT(q,j,Z) -= d.z * grphi * CHARGE(p,i);
#endif
            }
#else
	    /* Field for Dipole calculation */
	    if (dp_p_calc) {
	      Estat.x += d.x * grphi * charge[jt];
	      Estat.y += d.y * grphi * charge[jt];
	      Estat.z += d.z * grphi * charge[jt];
	      DP_E_STAT(q,j,X) -= d.x * grphi * charge[it];
	      DP_E_STAT(q,j,Y) -= d.y * grphi * charge[it];
	      DP_E_STAT(q,j,Z) -= d.z * grphi * charge[it];
	    
#ifdef EXTF
	      Estat.x += extf.x;
	      Estat.y += extf.y;
	      Estat.z += extf.z;
	      DP_E_STAT(q,j,X) += extf.x;
	      DP_E_STAT(q,j,Y) += extf.y;
	      DP_E_STAT(q,j,Z) += extf.z;
#endif
	    }
#endif
#endif
	  }
#ifdef DIPOLE
#ifdef VARCHG
	  /* calculate short-range dipoles field */
	  /* short-range fn.: 3rd column ff. */
	  if (dp_p_calc) {
	    col=(it <= jt) ?
	      it * ntypes + jt - ((it * (it + 1))/2)
	      : jt * ntypes + it - ((jt * (jt + 1))/2);
	    VAL_FUNC(pot,coul_table,2+col, 2+ntypepairs, r2, is_short);
	    tmp = pot*CHARGE(q,j)*dp_alpha[it];
	    if (SQR(tmp)>0) {
	      pstat.x -= tmp * d.x;
	      pstat.y -= tmp * d.y;
	      pstat.z -= tmp * d.z;
	    }
	    tmp = pot*CHARGE(p,i)*dp_alpha[jt];
	    if (SQR(tmp)>0){
	      DP_P_STAT(q,j,X) += tmp * d.x;
	      DP_P_STAT(q,j,Y) += tmp * d.y;
	      DP_P_STAT(q,j,Z) += tmp * d.z;
	    }
	  }
#else
	  /* calculate short-range dipoles field */
	  /* short-range fn.: 3rd column ff. */
	  if (dp_p_calc) {
	    col=(it <= jt) ?
	      it * ntypes + jt - ((it * (it + 1))/2)
	      : jt * ntypes + it - ((jt * (jt + 1))/2);
	    VAL_FUNC(pot,coul_table,2+col, 2+ntypepairs, r2, is_short);
	    tmp = pot*charge[jt]*dp_alpha[it];
	    if (SQR(tmp)>0) {
	      pstat.x -= tmp * d.x;
	      pstat.y -= tmp * d.y;
	      pstat.z -= tmp * d.z;
#ifdef EXTF
	      pstat.x += dp_alpha[it] * extf.x;
	      pstat.y += dp_alpha[it] * extf.y;
	      pstat.z += dp_alpha[it] * extf.z;
#endif
	    }
	    tmp = pot*charge[it]*dp_alpha[jt];
	    if (SQR(tmp)>0){
	      DP_P_STAT(q,j,X) += tmp * d.x;
	      DP_P_STAT(q,j,Y) += tmp * d.y;
	      DP_P_STAT(q,j,Z) += tmp * d.z;
#ifdef EXTF
	      DP_P_STAT(q,j,X) += dp_alpha[jt] * extf.x;
	      DP_P_STAT(q,j,Y) += dp_alpha[jt] * extf.y;
	      DP_P_STAT(q,j,Z) += dp_alpha[jt] * extf.z;
#endif
	    }
	  }
#endif
#endif /* DIPOLE */
	}
#endif /* COULOMB */

#ifdef COVALENT
        /* make neighbor tables for covalent systems */
        if (r2 < neightab_r2cut[col]) {

          neightab *neigh;

          /* update neighbor table of particle i */
          neigh = NEIGH(p,i);
          if (neigh->n_max <= neigh->n) {
            increase_neightab( neigh, neigh->n_max + NEIGH_LEN_INC );
          }
          neigh->typ[neigh->n] = jt;
          neigh->cl [neigh->n] = q;
          neigh->num[neigh->n] = j;
          neigh->dist[3*neigh->n  ] = d.x;
          neigh->dist[3*neigh->n+1] = d.y;
          neigh->dist[3*neigh->n+2] = d.z;
          neigh->n++;

          /* update neighbor table of particle j */
          neigh = NEIGH(q,j);
          if (neigh->n_max <= neigh->n) {
            increase_neightab( neigh, neigh->n_max + NEIGH_LEN_INC );
          }
          neigh->typ[neigh->n] = it;
          neigh->cl [neigh->n] = p;
          neigh->num[neigh->n] = i;
          neigh->dist[3*neigh->n  ] = -d.x;
          neigh->dist[3*neigh->n+1] = -d.y;
          neigh->dist[3*neigh->n+2] = -d.z;
          neigh->n++;
        }
#endif  /* COVALENT */

      }
      KRAFT(p,i,X) += ff.x;
      KRAFT(p,i,Y) += ff.y;
      KRAFT(p,i,Z) += ff.z;
#ifndef MONOLJ
      POTENG(p,i)  += ee;
#endif
#ifdef EAM2
      EAM_RHO(p,i) += eam_r;
#ifdef EEAM
      EAM_P(p,i)   += eam_p;
#endif
#endif
#ifdef ADP
      ADP_MU    (p,i,X)  += mu.x;
      ADP_MU    (p,i,Y)  += mu.y;
      ADP_MU    (p,i,Z)  += mu.z;
      ADP_LAMBDA(p,i,xx) += la.xx;
      ADP_LAMBDA(p,i,yy) += la.yy;
      ADP_LAMBDA(p,i,zz) += la.zz;
      ADP_LAMBDA(p,i,yz) += la.yz;
      ADP_LAMBDA(p,i,zx) += la.zx;
      ADP_LAMBDA(p,i,xy) += la.xy;
#endif
#ifdef DIPOLE
      if (dp_p_calc) {
	DP_E_STAT(p,i,X)   += Estat.x;
	DP_E_STAT(p,i,Y)   += Estat.y;
	DP_E_STAT(p,i,Z)   += Estat.z;
	DP_P_STAT(p,i,X)   += pstat.x;
	DP_P_STAT(p,i,Y)   += pstat.y;
	DP_P_STAT(p,i,Z)   += pstat.z;
	/* Field Extrapolation */
	if (dp_E_calc>2) {
	  DP_E_IND(p,i,X) = 3.*DP_E_OLD_1(p,i,X) - 3.*DP_E_OLD_2(p,i,X) +
	    DP_E_OLD_3(p,i,X);
	  DP_E_IND(p,i,Y) = 3.*DP_E_OLD_1(p,i,Y) - 3.*DP_E_OLD_2(p,i,Y) +
	    DP_E_OLD_3(p,i,Y);
	  DP_E_IND(p,i,Z) = 3.*DP_E_OLD_1(p,i,Z) - 3.*DP_E_OLD_2(p,i,Z) +
	    DP_E_OLD_3(p,i,Z);
	  DP_E_OLD_3(p,i,X) = 0.;
	  DP_E_OLD_3(p,i,Y) = 0.;
	  DP_E_OLD_3(p,i,Z) = 0.;
	} else {
	  DP_E_IND(p,i,X) = DP_E_OLD_1(p,i,X);
	  DP_E_IND(p,i,Y) = DP_E_OLD_1(p,i,Y);
	  DP_E_IND(p,i,Z) = DP_E_OLD_1(p,i,Z);
	}
      }
#endif
#ifdef STRESS_TENS
      if (do_press_calc) {
        PRESSTENS(p,i,xx) += pp.xx;
        PRESSTENS(p,i,yy) += pp.yy;
        PRESSTENS(p,i,xy) += pp.xy;
#ifndef TWOD
        PRESSTENS(p,i,zz) += pp.zz;
        PRESSTENS(p,i,yz) += pp.yz;
        PRESSTENS(p,i,zx) += pp.zx;
#endif
      }
#endif
#ifdef NNBR
      NBANZ(p,i)    += nb;
#endif
      n++;
    }
#ifdef DIPOLE
    if (dp_p_calc) {
      dp_E_shift = p->dp_E_old_3;
      p->dp_E_old_3 = p->dp_E_old_2;
      p->dp_E_old_2 = p->dp_E_old_1;
      p->dp_E_old_1 = dp_E_shift;
    }
#endif /* DIPOLE */
  }
  if (is_short) fprintf(stderr,"Short distance, pair, step %d!\n",steps);

#ifdef EWALD
  if (steps==0) {
    imd_stop_timer( &ewald_time );
  }
#endif

#ifdef COVALENT

  /* complete neighbor tables for covalent systems */
  for (k=ncells; k<ncells2; k++) {
    cell *p = cell_array +cnbrs[k].np;
    for (i=0; i<p->n; i++) {

      vektor d1;
      int    m, it;

      d1.x = ORT(p,i,X);
      d1.y = ORT(p,i,Y);
      d1.z = ORT(p,i,Z);
      it   = SORTE(p,i);

      /* loop over neighbors */
      for (m=tl[n]; m<tl[n+1]; m++) {

        int    c, j, jt, col, inc = ntypes * ntypes; 
        vektor d;
        real   r2;
        cell   *q;

        c = cl_num[ tb[m] ];
        j = tb[m] - cl_off[c];
        q = cell_array + c;

        d.x = ORT(q,j,X) - d1.x;
        d.y = ORT(q,j,Y) - d1.y;
        d.z = ORT(q,j,Z) - d1.z;
        r2  = SPROD(d,d);
        jt  = SORTE(q,j);
        col = it * ntypes + jt;

        /* make neighbor tables */
        if (r2 <= neightab_r2cut[col]) {

          neightab *neigh;

          /* update neighbor table of particle i */
          neigh = NEIGH(p,i);
          if (neigh->n_max <= neigh->n) {
            increase_neightab( neigh, neigh->n_max + NEIGH_LEN_INC );
          }
          neigh->typ[neigh->n] = jt;
          neigh->cl [neigh->n] = q;
          neigh->num[neigh->n] = j;
          neigh->dist[3*neigh->n  ] = d.x;
          neigh->dist[3*neigh->n+1] = d.y;
          neigh->dist[3*neigh->n+2] = d.z;
          neigh->n++;

          /* update neighbor table of particle j */
          /* we do not need a neighbor table in buffer cells
          neigh = NEIGH(q,j);
          if (neigh->n_max <= neigh->n) {
            increase_neightab( neigh, neigh->n_max + NEIGH_LEN_INC );
          }
          neigh->typ[neigh->n] = it;
          neigh->cl [neigh->n] = p;
          neigh->num[neigh->n] = i;
          neigh->dist[3*neigh->n  ] = -d.x;
          neigh->dist[3*neigh->n+1] = -d.y;
          neigh->dist[3*neigh->n+2] = -d.z;
          neigh->n++;
          */
        }
      }
      n++;
    }
  }

#ifndef CNA
  /* second force loop for covalent systems */
  for (k=0; k<ncells; ++k) {
    do_forces2(cell_array + cnbrs[k].np,
               &tot_pot_energy, &virial, &vir_xx, &vir_yy, &vir_zz,
                                         &vir_yz, &vir_zx, &vir_xy);
  }
#endif

#endif /* COVALENT */

#ifdef EAM2

  /* collect host electron density */
  send_forces(add_rho,pack_rho,unpack_add_rho);

  /* compute embedding energy and its derivative */
  for (k=0; k<ncells; k++) {
    cell *p = CELLPTR(k);
    real pot, tmp, tr;
#ifdef ia64
#pragma ivdep,swp
#endif
    for (i=0; i<p->n; i++) {
      PAIR_INT( pot, EAM_DF(p,i), embed_pot, SORTE(p,i), 
                ntypes, EAM_RHO(p,i), idummy);
      POTENG(p,i)    += pot;
      tot_pot_energy += pot;
#ifdef EEAM
      PAIR_INT( pot, EAM_DM(p,i), emod_pot, SORTE(p,i), 
                ntypes, EAM_P(p,i), idummy);
      POTENG(p,i)    += pot;
      tot_pot_energy += pot;
#endif
#ifdef ADP
      tr  = (ADP_LAMBDA(p,i,xx) + ADP_LAMBDA(p,i,yy) + ADP_LAMBDA(p,i,zz))/3.0;
      tmp = ADP_LAMBDA(p,i,xx) - tr; pot  = SQR(tmp);
      tmp = ADP_LAMBDA(p,i,yy) - tr; pot += SQR(tmp);
      tmp = ADP_LAMBDA(p,i,zz) - tr; pot += SQR(tmp);
      tmp = ADP_LAMBDA(p,i,yz);      pot += SQR(tmp) * 2.0;
      tmp = ADP_LAMBDA(p,i,zx);      pot += SQR(tmp) * 2.0;
      tmp = ADP_LAMBDA(p,i,xy);      pot += SQR(tmp) * 2.0;
      tmp = ADP_MU    (p,i,X);       pot += SQR(tmp);
      tmp = ADP_MU    (p,i,Y);       pot += SQR(tmp);
      tmp = ADP_MU    (p,i,Z);       pot += SQR(tmp);
      pot *= 0.5;
      POTENG(p,i)    += pot;
      tot_pot_energy += pot;
#endif
    }
  }

  /* distribute derivative of embedding energy */
  send_cells(copy_dF,pack_dF,unpack_dF);

  /* EAM interactions - for all atoms */
  n=0;
  for (k=0; k<ncells; k++) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; i++) {

#ifdef STRESS_TENS
      sym_tensor pp = {0.0,0.0,0.0,0.0,0.0,0.0};
#endif
#ifdef ADP
      sym_tensor la1;
      vektor mu1;
#endif
      vektor d1, ff = {0.0,0.0,0.0};
      int m, it;

      d1.x = ORT(p,i,X);
      d1.y = ORT(p,i,Y);
      d1.z = ORT(p,i,Z);
#ifdef ADP
      mu1.x  = ADP_MU    (p,i,X);
      mu1.y  = ADP_MU    (p,i,Y);
      mu1.z  = ADP_MU    (p,i,Z);
      la1.xx = ADP_LAMBDA(p,i,xx);
      la1.yy = ADP_LAMBDA(p,i,yy);
      la1.zz = ADP_LAMBDA(p,i,zz);
      la1.yz = ADP_LAMBDA(p,i,yz);
      la1.zx = ADP_LAMBDA(p,i,zx);
      la1.xy = ADP_LAMBDA(p,i,xy);
#endif
      it   = SORTE(p,i);

      /* loop over neighbors */
#ifdef ia64
#pragma ivdep,swp
#endif
      for (m=tl[n]; m<tl[n+1]; m++) {

        vektor d, force = {0.0,0.0,0.0};
        real   r2;
        int    c, j, jt, col1, col2, inc = ntypes * ntypes, have_force=0;
        cell   *q;

        c = cl_num[ tb[m] ];
        j = tb[m] - cl_off[c];
        q = cell_array + c;

        d.x  = ORT(q,j,X) - d1.x;
        d.y  = ORT(q,j,Y) - d1.y;
        d.z  = ORT(q,j,Z) - d1.z;
        r2   = SPROD(d,d);
        jt   = SORTE(q,j);
        col1 = jt * ntypes + it;
        col2 = it * ntypes + jt;

        if ((r2 < rho_h_tab.end[col1]) || (r2 < rho_h_tab.end[col2])) {

          real pot, grad, rho_i_strich, rho_j_strich, rho_i, rho_j;

          /* take care: particle i gets its rho from particle j.    */
          /* This is tabulated in column it*ntypes+jt.              */
          /* Here we need the giving part from column jt*ntypes+it. */

          /* rho_strich_i(r_ij) */
#ifndef EEAM
          DERIV_FUNC(rho_i_strich, rho_h_tab, col1, inc, r2, is_short);
#else
          /* rho_strich_i(r_ij) and rho_i(r_ij) */
          PAIR_INT(rho_i, rho_i_strich, rho_h_tab, col1, inc, r2, is_short);
#endif

          /* rho_strich_j(r_ij) */
          if (col1==col2) {
            rho_j_strich = rho_i_strich;
#ifdef EEAM
            rho_j = rho_i;
#endif
          } else {
#ifndef EEAM
            DERIV_FUNC(rho_j_strich, rho_h_tab, col2, inc, r2, is_short);
#else
            PAIR_INT(rho_j, rho_j_strich, rho_h_tab, col2, inc, r2, is_short);
#endif
	  }

          /* put together (dF_i and dF_j are by 0.5 too big) */
          grad = 0.5 * (EAM_DF(p,i)*rho_j_strich + EAM_DF(q,j)*rho_i_strich);
#ifdef EEAM
          /* 0.5 times 2 from derivative simplified to 1 */
          grad += (EAM_DM(p,i) * rho_j * rho_j_strich +
                   EAM_DM(q,j) * rho_i * rho_i_strich);
#endif

          /* store force in temporary variable */
          force.x = d.x * grad;
          force.y = d.y * grad;
          force.z = d.z * grad;
          have_force=1;
        }

#ifdef ADP
        /* forces due to dipole distortion */
        if (r2 < adp_upot.end[col1]) {
          vektor mu;
          real pot, grad, tmp;
          PAIR_INT(pot, grad, adp_upot, col1, inc, r2, is_short);
          mu.x = mu1.x - ADP_MU(q,j,X);
          mu.y = mu1.y - ADP_MU(q,j,Y);
          mu.z = mu1.z - ADP_MU(q,j,Z);
          tmp  = SPROD(mu,d) * grad;
          force.x += mu.x * pot + tmp * d.x;
          force.y += mu.y * pot + tmp * d.y;
          force.z += mu.z * pot + tmp * d.z;
          have_force=1;
        }
        /* forces due to quadrupole distortion */
        if (r2 < adp_wpot.end[col1]) {
          sym_tensor la;
          vektor v;
          real pot, grad, nu, f1, f2;
          PAIR_INT(pot, grad, adp_wpot, col1, inc, r2, is_short);
          la.xx = la1.xx + ADP_LAMBDA(q,j,xx);
          la.yy = la1.yy + ADP_LAMBDA(q,j,yy);
          la.zz = la1.zz + ADP_LAMBDA(q,j,zz);
          la.yz = la1.yz + ADP_LAMBDA(q,j,yz);
          la.zx = la1.zx + ADP_LAMBDA(q,j,zx);
          la.xy = la1.xy + ADP_LAMBDA(q,j,xy);
          v.x = la.xx * d.x + la.xy * d.y + la.zx * d.z;
          v.y = la.xy * d.x + la.yy * d.y + la.yz * d.z;
          v.z = la.zx * d.x + la.yz * d.y + la.zz * d.z;
          nu  = (la.xx + la.yy + la.zz) / 3.0;
          f1  = 2.0 * pot;
          f2  = (SPROD(v,d) - nu * r2) * grad - nu * f1; 
          force.x += f1 * v.x + f2 * d.x;
          force.y += f1 * v.y + f2 * d.y;
          force.z += f1 * v.z + f2 * d.z;
          have_force=1;
        }
#endif

#ifdef FLAGEDATOMS
	  if(VSORTE(q,j) == flagedatomstype && VSORTE(p,i) == flagedatomstype)
	    {
	      //	      printf("Atom nr %d of type %d interacting with %d: Embed forces : %e %e %e\n",
	      //     NUMMER(p,i),VSORTE(p,i),NUMMER(q,j),force.x,force.y,force.z);
	      printf("%d %d %d %e %e %e %e %e %e\n",
		     NUMMER(p,i),VSORTE(p,i),NUMMER(q,j),d.x, d.y, d.z, force.x,force.y,force.z);
	      fflush(stdout);
	    }
#endif
        /* accumulate forces */
        if (have_force) {
          KRAFT(q,j,X) -= force.x;
          KRAFT(q,j,Y) -= force.y;
          KRAFT(q,j,Z) -= force.z;
          ff.x         += force.x;
          ff.y         += force.y;
          ff.z         += force.z;
#ifdef P_AXIAL
          vir_xx       -= d.x * force.x;
          vir_yy       -= d.y * force.y;
          vir_zz       -= d.z * force.z;
#else
          virial       -= SPROD(d,force);
#endif

#ifdef STRESS_TENS
          if (do_press_calc) {
            /* avoid double counting of the virial */
            force.x *= 0.5;
            force.y *= 0.5;
            force.z *= 0.5;
 
            pp.xx -= d.x * force.x;
            pp.yy -= d.y * force.y;
            pp.zz -= d.z * force.z;
            pp.yz -= d.y * force.z;
            pp.zx -= d.z * force.x;
            pp.xy -= d.x * force.y;

            PRESSTENS(q,j,xx) -= d.x * force.x;
            PRESSTENS(q,j,yy) -= d.y * force.y;
            PRESSTENS(q,j,zz) -= d.z * force.z;
            PRESSTENS(q,j,yz) -= d.y * force.z;
            PRESSTENS(q,j,zx) -= d.z * force.x;
            PRESSTENS(q,j,xy) -= d.x * force.y;
          }
#endif
        }
      }
      KRAFT(p,i,X) += ff.x;
      KRAFT(p,i,Y) += ff.y;
      KRAFT(p,i,Z) += ff.z;
#ifdef STRESS_TENS
      if (do_press_calc) {
        PRESSTENS(p,i,xx) += pp.xx;
        PRESSTENS(p,i,yy) += pp.yy;
        PRESSTENS(p,i,zz) += pp.zz;
        PRESSTENS(p,i,yz) += pp.yz;
        PRESSTENS(p,i,zx) += pp.zx;
        PRESSTENS(p,i,xy) += pp.xy;
      }
#endif
      n++;
    }
  }
  if (is_short) fprintf(stderr, "\n Short distance, EAM, step %d!\n",steps);

#endif /* EAM2 */

#ifdef COULOMB
  /* coulomb self energy part */
  for (k=0; k<ncells; k++) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; i++) {
      real chg = CHARGE(p,i);
#ifdef SM
      int  t   = SORTE(p,i);
      real pot = (0.5*(2*ew_vorf * coul_eng - sm_J_0[t]) * chg - sm_chi_0[t]) * chg;
#else
      real pot = ew_vorf * SQR(chg) * coul_eng;
#endif
#ifdef SM
      tot_sm_es_energy -= pot;
#endif
      tot_pot_energy -= pot;
      POTENG(p,i)    -= pot;
    }
  }
#endif

#ifdef DIPOLE

  /* collect constant electric fields and dipole moments */
  send_forces(add_dipole,pack_dipole,unpack_add_dipole);
  
  if (dp_p_calc) {
    while (dp_converged==0) {
      dp_sum_old=dp_sum;
      dp_sum=0.0;
      n=0;
      /* Set field, dipoles */
      for (k=0; k<ncells; k++) { 
	cell *p = CELLPTR(k);
#ifdef ia64 
#pragma ivdep,swp
#endif
	for (i=0; i<p->n; i++) {
	  int it;
	  vektor Etot;
	  it   = SORTE(p,i);
	  if (SQR(dp_alpha[it])>0) {
	    if (dp_it) {
	      Etot.x = dp_mix * DP_E_IND(p,i,X)
		+ (1.-dp_mix) * DP_E_OLD_1(p,i,X)
		+ DP_E_STAT(p,i,X);
	      Etot.y = dp_mix * DP_E_IND(p,i,Y)
		+ (1.-dp_mix) * DP_E_OLD_1(p,i,Y)
		+ DP_E_STAT(p,i,Y);
	      Etot.z = dp_mix * DP_E_IND(p,i,Z)
		+ (1.-dp_mix) * DP_E_OLD_1(p,i,Z)
		+ DP_E_STAT(p,i,Z);
	    } else {
	      Etot.x = DP_E_IND(p,i,X) + DP_E_STAT(p,i,X);
	      Etot.y = DP_E_IND(p,i,Y) + DP_E_STAT(p,i,Y);
	      Etot.z = DP_E_IND(p,i,Z) + DP_E_STAT(p,i,Z);
	    }
	    DP_P_IND(p,i,X) = Etot.x * dp_alpha[it] + DP_P_STAT(p,i,X);
	    DP_P_IND(p,i,Y) = Etot.y * dp_alpha[it] + DP_P_STAT(p,i,Y);
	    DP_P_IND(p,i,Z) = Etot.z * dp_alpha[it] + DP_P_STAT(p,i,Z);
	    
/* #ifdef DEBUG */
/* 	    printf("%d\t%d\t%g\t%g\t%g\t", NUMMER(p,i),it,ORT(p,i,X), */
/* 		   ORT(p,i,Y),ORT(p,i,Z)); */
/* 	    printf("%10.6f\t%10.6f\t%10.6f\t",Etot.x,Etot.y,Etot.z); */
/* 	    printf("%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t", */
/* 		   DP_E_IND(p,i,X), DP_E_IND(p,i,Y),DP_E_IND(p,i,Z), */
/* 		   DP_P_IND(p,i,X), DP_P_IND(p,i,Y),DP_P_IND(p,i,Z)); */
/* 	    printf("%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\n", */
/* 		   DP_P_STAT(p,i,X), DP_P_STAT(p,i,Y),DP_P_STAT(p,i,Z), */
/* 		   DP_E_STAT(p,i,X), DP_E_STAT(p,i,Y),DP_E_STAT(p,i,Z)); */
/* 	  printf("%d\t1\t%10.6f\t%10.6f\t%10.6f\t%10.6f\n",NUMMER(p,i),Etot.x, */
/* 		 DP_E_IND(p,i,X), DP_P_IND(p,i,X), DP_P_STAT(p,i,X)); */
/* 	  printf("%d\t2\t%10.6f\t%10.6f\t%10.6f\t%10.6f\n",NUMMER(p,i),Etot.y, */
/* 		 DP_E_IND(p,i,Y), DP_P_IND(p,i,Y), DP_P_STAT(p,i,Y)); */
/* 	  printf("%d\t3\t%10.6f\t%10.6f\t%10.6f\t%10.6f\n",NUMMER(p,i),Etot.z, */
/* 		 DP_E_IND(p,i,Z), DP_P_IND(p,i,Z), DP_P_STAT(p,i,Z)); */
/* #endif */
	    
	    DP_E_OLD_1(p,i,X) = DP_E_IND(p,i,X);
	    DP_E_OLD_1(p,i,Y) = DP_E_IND(p,i,Y);
	    DP_E_OLD_1(p,i,Z) = DP_E_IND(p,i,Z);
	    /* Self Term */
/* 	    DP_E_IND(p,i,X) = DP_P_IND(p,i,X) * dp_self; */
/* 	    DP_E_IND(p,i,Y) = DP_P_IND(p,i,Y) * dp_self; */
/* 	    DP_E_IND(p,i,Z) = DP_P_IND(p,i,Z) * dp_self; */
	    DP_E_IND(p,i,X) = 0.;
	    DP_E_IND(p,i,Y) = 0.;
	    DP_E_IND(p,i,Z) = 0.;
	  }
	}
      }
      /* Distribute dipole moments */
      send_cells(copy_pind,pack_pind,unpack_pind);

      /* compute DIPOLE strength  */
      for (k=0; k<ncells; k++) { 
	cell *p = CELLPTR(k); 
#ifdef ia64 
#pragma ivdep,swp
#endif
	for (i=0; i<p->n; i++) {
	  int    m, it, nb = 0;
	  vektor d1, pi, Eind={0.0,0.0,0.0};
	  it   = SORTE(p,i);

	  if ( SQR(dp_alpha[it])>0) { 

	    d1.x = ORT(p,i,X);
	    d1.y = ORT(p,i,Y);
	    d1.z = ORT(p,i,Z);

      
	    pi.x = DP_P_IND(p,i,X);
	    pi.y = DP_P_IND(p,i,Y);
	    pi.z = DP_P_IND(p,i,Z);
      

	    /* loop over neighbors */
#ifdef ia64
#pragma ivdep,swp
#endif
	    for (m=tl[n]; m<tl[n+1]; m++) {

	      vektor d;
	      real   r2;
	      int    c, j, jt, col1, col2, inc = ntypes * ntypes, have_force=0;
	      cell   *q;
	  
	      c = cl_num[ tb[m] ];
	      j = tb[m] - cl_off[c];
	      q = cell_array + c;

	      d.x  = ORT(q,j,X) - d1.x;
	      d.y  = ORT(q,j,Y) - d1.y;
	      d.z  = ORT(q,j,Z) - d1.z;
	      r2   = SPROD(d,d);
	      jt   = SORTE(q,j);
	      col1 = jt * ntypes + it;
	      col2 = it * ntypes + jt;
	  
	      /* Dipole-Dipole */
	      if ( (SQR(dp_alpha[jt])>0) && (r2 < ew_r2_cut )) {
		real pot,tmp;
		vektor pj;
		pj.x=DP_P_IND(q,j,X);
		pj.y=DP_P_IND(q,j,Y);
		pj.z=DP_P_IND(q,j,Z);
		/* smooth r^3 cutoff */
		VAL_FUNC(pot,coul_table,1,2+ntypepairs,r2,is_short);
		pot *= coul_eng;
		tmp=SPROD(pj,d);
		Eind.x += pot* ((3.0/r2)*tmp*d.x - pj.x);
		Eind.y += pot* ((3.0/r2)*tmp*d.y - pj.y);
		Eind.z += pot* ((3.0/r2)*tmp*d.z - pj.z);
		tmp=SPROD(pi,d);
		DP_E_IND(q,j,X) += pot * ((3.0/r2)*tmp*d.x -pi.x);
		DP_E_IND(q,j,Y) += pot * ((3.0/r2)*tmp*d.y -pi.y);
		DP_E_IND(q,j,Z) += pot * ((3.0/r2)*tmp*d.z -pi.z);
	      }	  
	    }
	    DP_E_IND(p,i,X) += Eind.x;
	    DP_E_IND(p,i,Y) += Eind.y;
	    DP_E_IND(p,i,Z) += Eind.z;
	  }
	  n++;
	}
      }
      /* Collect Electric fields */
      send_forces(add_field,pack_field,unpack_add_field);
      n=0;
      for (k=0; k<ncells; k++) {
	cell *p = CELLPTR(k);
	for (i=0; i<p->n; i++) {
	  int it;
	  real dp_tmp;
	  it   = SORTE(p,i);
	  if(SQR(dp_alpha[it])>0){
	    dp_sum += SQR(dp_alpha[it]*(DP_E_OLD_1(p,i,X)-DP_E_IND(p,i,X)));
	    dp_sum += SQR(dp_alpha[it]*(DP_E_OLD_1(p,i,Y)-DP_E_IND(p,i,Y)));
	    dp_sum += SQR(dp_alpha[it]*(DP_E_OLD_1(p,i,Y)-DP_E_IND(p,i,Y)));
	  }
	}
      }
      dp_sum /= 3.0*natoms;
      dp_sum=sqrt(dp_sum);
#ifdef DEBUG
      printf("#dipole deviation at step %d: %g\n",steps,dp_sum);
#endif /*DEBUG */
      if ((dp_sum > max_diff) || ( dp_it>50)) { 
	fprintf(stderr, "\n Convergence Error, dipole, step %d: ", \
			steps);
	fprintf(stderr,"dp_sum = %g, dp_it=%d \n",dp_sum,dp_it);
	n=0;
	for (k=0; k<ncells; k++) {
	  cell *p = CELLPTR(k);
	  for (i=0; i<p->n; i++) {
	    int it;
	    real dp_tmp;
	    it   = SORTE(p,i);
	    if(SQR(dp_alpha[it])>0){
	      DP_P_IND(p,i,X) = dp_alpha[it]*DP_E_STAT(p,i,X)+DP_P_STAT(p,i,X);
	      DP_P_IND(p,i,Y) = dp_alpha[it]*DP_E_STAT(p,i,Y)+DP_P_STAT(p,i,Y);
	      DP_P_IND(p,i,Z) = dp_alpha[it]*DP_E_STAT(p,i,Z)+DP_P_STAT(p,i,Z);
	      DP_E_IND(p,i,X) = DP_E_STAT(p,i,X);
	      DP_E_IND(p,i,Y) = DP_E_STAT(p,i,Y);
	      DP_E_IND(p,i,Z) = DP_E_STAT(p,i,Z);
	    }
	    n++;
	  }
	}
	send_cells(copy_pind,pack_pind,unpack_pind);
	dp_converged=1;
	dp_sum=1.;
      }
/*       if  (fabs(dp_sum)-dp_sum_old) < dp_tol) /\* reasonable? *\/ */
      if (dp_sum < dp_tol)
	dp_converged=1;
      dp_it++;

    } /* Dipole iteration */
  }
  /* DIPOLE interactions - for all atoms */

#ifdef DEBUG
/* Don't use this unless you are prepared for massive output */
//  printf("No,type,pos,force,p_ind,E_ind,p_stat,E_stat\n");
#endif

  n=0;
  for (k=0; k<ncells; k++) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; i++) {

#ifdef STRESS_TENS
      sym_tensor pp = {0.0,0.0,0.0,0.0,0.0,0.0};
#endif
      vektor d1, pi, ff = {0.0,0.0,0.0};
      real   ee = 0.0;
      real   dp_energy = 0.0;
      int    m, it, nb = 0;
       
      d1.x = ORT(p,i,X);
      d1.y = ORT(p,i,Y);
      d1.z = ORT(p,i,Z);
      it   = SORTE(p,i);

      
      pi.x = DP_P_IND(p,i,X);
      pi.y = DP_P_IND(p,i,Y);
      pi.z = DP_P_IND(p,i,Z);

#ifdef VARCHG
real pconst_i=SQR(CHARGE(p,i))+SQR(dp_alpha[it]);
#else
real pconst_i=SQR(charge[it])+SQR(dp_alpha[it]);
#endif

      if (pconst_i>0) { 
	/* loop over neighbors */
#ifdef ia64
#pragma ivdep,swp
#endif
	for (m=tl[n]; m<tl[n+1]; m++) {

	  vektor d, force = {0.0,0.0,0.0};
	  real   r2;
	  int    c, j, jt, col1, col2, inc = ntypes * ntypes, have_force=0;
	  cell   *q;
	  
	  c = cl_num[ tb[m] ];
	  j = tb[m] - cl_off[c];
	  q = cell_array + c;

	  d.x  = ORT(q,j,X) - d1.x;
	  d.y  = ORT(q,j,Y) - d1.y;
	  d.z  = ORT(q,j,Z) - d1.z;
	  r2   = SPROD(d,d);
	  jt   = SORTE(q,j);
	  col1 = jt * ntypes + it;
	  col2 = it * ntypes + jt;
	  
	  if (r2 < ew_r2_cut ) {	
	    /* XXX */
	    real pot, grad, val, dval, pdotd, tmp, pdotp, proj_i, proj_j;
	    real valsr,dvalsr;
	    vektor pj;

	    pot=0.0; grad=0.0;

	    pj.x=DP_P_IND(q,j,X);
	    pj.y=DP_P_IND(q,j,Y);
	    pj.z=DP_P_IND(q,j,Z);

	    /* smooth cutoff function is 2nd in coul_table */
	    PAIR_INT(val, dval, coul_table, 1, 2+ntypepairs, r2, is_short);
	    
	    val  *= coul_eng;
	    dval *= coul_eng;
	    /* short-range function is 3rd in coul_table */
	    col1=(it <= jt) ?
	      it * ntypes + jt - ((it * (it + 1))/2)
	      :jt * ntypes + it - ((jt * (jt + 1))/2);
	    PAIR_INT(valsr,dvalsr,coul_table,2+col1,2+ntypepairs,
		     r2,is_short);
	    /* Dipole at i -Charge at j interaction */
#ifdef VARCHG
	      if (SQR(dp_alpha[it])*SQR(CHARGE(q,j))>0) {

	      pdotd=SPROD(pi,d);
	      pot += val*CHARGE(q,j)*pdotd;
	      force.x += CHARGE(q,j) * (dval * pdotd * d.x + val * pi.x);
	      force.y += CHARGE(q,j) * (dval * pdotd * d.y + val * pi.y);
	      force.z += CHARGE(q,j) * (dval * pdotd * d.z + val * pi.z);
	      
	      /* short range interaction */
	      pot += CHARGE(q,j)*pdotd*valsr;
	      force.x += CHARGE(q,j) * (dvalsr*pdotd *d.x +valsr*pi.x);
	      force.y += CHARGE(q,j) * (dvalsr*pdotd *d.y +valsr*pi.y);
	      force.z += CHARGE(q,j) * (dvalsr*pdotd *d.z +valsr*pi.z);

	      have_force=1;
	    }
#else
	    if (SQR(dp_alpha[it])*SQR(charge[jt])>0) {

	      pdotd=SPROD(pi,d);
	      pot += val*charge[jt]*pdotd;
	      force.x += charge[jt] * (dval * pdotd * d.x + val * pi.x);
	      force.y += charge[jt] * (dval * pdotd * d.y + val * pi.y);
	      force.z += charge[jt] * (dval * pdotd * d.z + val * pi.z);
	      
	      /* short range interaction */
	      pot += charge[jt]*pdotd*valsr;
	      force.x += charge[jt] * (dvalsr*pdotd *d.x +valsr*pi.x);
	      force.y += charge[jt] * (dvalsr*pdotd *d.y +valsr*pi.y);
	      force.z += charge[jt] * (dvalsr*pdotd *d.z +valsr*pi.z);

	      have_force=1;
	    }
#endif
	    /* Dipole at j -Charge at i interaction */
#ifdef VARCHG
	    if (SQR(dp_alpha[jt])*SQR(CHARGE(p,i))>0) {
	      /* smooth cutoff function is 2nd in coul_table */

	      pdotd=SPROD(pj,d);
	      pot -= val*CHARGE(p,i)*pdotd;
	      force.x -= CHARGE(p,i) * (dval * pdotd * d.x + val * pj.x);
	      force.y -= CHARGE(p,i) * (dval * pdotd * d.y + val * pj.y);
	      force.z -= CHARGE(p,i) * (dval * pdotd * d.z + val * pj.z);

	      /* short range interaction */
	      pot -= CHARGE(p,i)*pdotd*valsr;
	      force.x -= CHARGE(p,i) * (dvalsr*pdotd *d.x +valsr*pj.x);
	      force.y -= CHARGE(p,i) * (dvalsr*pdotd *d.y +valsr*pj.y);
	      force.z -= CHARGE(p,i) * (dvalsr*pdotd *d.z +valsr*pj.z);

	      have_force=1;
	    }
#else
	    if (SQR(dp_alpha[jt])*SQR(charge[it])>0) {
	      /* smooth cutoff function is 2nd in coul_table */

	      pdotd=SPROD(pj,d);
	      pot -= val*charge[it]*pdotd;
	      force.x -= charge[it] * (dval * pdotd * d.x + val * pj.x);
	      force.y -= charge[it] * (dval * pdotd * d.y + val * pj.y);
	      force.z -= charge[it] * (dval * pdotd * d.z + val * pj.z);

	      /* short range interaction */
	      pot -= charge[it]*pdotd*valsr;
	      force.x -= charge[it] * (dvalsr*pdotd *d.x +valsr*pj.x);
	      force.y -= charge[it] * (dvalsr*pdotd *d.y +valsr*pj.y);
	      force.z -= charge[it] * (dvalsr*pdotd *d.z +valsr*pj.z);

	      have_force=1;
	    }
#endif
	    /* Dipole-Dipole interaction */
	    if (SQR(dp_alpha[it])*SQR(dp_alpha[jt])>0) {
	      pdotp=SPROD(pi,pj);
	      proj_i=SPROD(pi,d);
	      proj_j=SPROD(pj,d);
	      tmp = pdotp - 3.0 * proj_i * proj_j /r2; 
	      pot += val * tmp;
	      force.x -= -dval * tmp * d.x 
		- val * ( 6.0 / SQR(r2) *proj_i *proj_j * d.x 
			  - 3.0 / r2 * (proj_i * pj.x + proj_j * pi.x));
	      force.y -= -dval * tmp * d.y 
		- val * ( 6.0 / SQR(r2) *proj_i *proj_j * d.y 
			  - 3.0 / r2 * (proj_i * pj.y + proj_j * pi.y));
	      force.z -= -dval * tmp * d.z 
		- val * ( 6.0 / SQR(r2) *proj_i *proj_j * d.z 
			  - 3.0 / r2 * (proj_i * pj.z + proj_j * pi.z));
	      have_force = 1;
	    }

	    if (have_force) {
	      tot_pot_energy += pot;

	      KRAFT(q,j,X) -= force.x;
	      KRAFT(q,j,Y) -= force.y;
	      KRAFT(q,j,Z) -= force.z;
	      ff.x         += force.x;
	      ff.y         += force.y;
	      ff.z         += force.z;
              pot          *= 0.5;   /* avoid double counting */
	      ee           += pot;
	      POTENG(q,j)  += pot;
	      
#ifdef P_AXIAL
	      vir_xx -= d.x * force.x;
	      vir_yy -= d.y * force.y;
	      vir_zz -= d.z * force.z;
#else
	      virial       -= SPROD(d,force);
#endif
	    
#ifdef STRESS_TENS
	      if (do_press_calc) {
		/* avoid double counting of the virial */
		force.x *= 0.5;
		force.y *= 0.5;
		force.z *= 0.5;
		pp.xx             -= d.x * force.x;
		PRESSTENS(q,j,xx) -= d.x * force.x;
		pp.yy             -= d.y * force.y;
		PRESSTENS(q,j,yy) -= d.y * force.y;
		pp.xy             -= d.x * force.y;
		PRESSTENS(q,j,xy) -= d.x * force.y;
		pp.zz             -= d.z * force.z;
		PRESSTENS(q,j,zz) -= d.z * force.z;
		pp.yz             -= d.y * force.z;
		PRESSTENS(q,j,yz) -= d.y * force.z;
		pp.zx             -= d.z * force.x;
		PRESSTENS(q,j,zx) -= d.z * force.x;
	      }
#endif
	    }
	  }	  
        }	
        /* Dipole self energy */
	dp_energy = (SQR(dp_alpha[it])>0)?0.5*SPROD(pi,pi)/dp_alpha[it]:0.;
	tot_pot_energy += dp_energy; /* avoid double counting */
	ee += dp_energy;

	POTENG(p,i)  += ee;
	KRAFT(p,i,X) += ff.x;
	KRAFT(p,i,Y) += ff.y;
	KRAFT(p,i,Z) += ff.z;
#ifdef STRESS_TENS
	if (do_press_calc) {
	  PRESSTENS(p,i,xx) += pp.xx;
	  PRESSTENS(p,i,yy) += pp.yy;
	  PRESSTENS(p,i,zz) += pp.zz;
	  PRESSTENS(p,i,yz) += pp.yz;
	  PRESSTENS(p,i,zx) += pp.zx;
	  PRESSTENS(p,i,xy) += pp.xy;
	}
#endif
      }
      n++;
#ifdef DEBUG
/* Don't use this unless you are prepared for massive output*/

/*       printf("%d\t%d\t%10.6f\t%10.6f\t%10.6f\t",  */
/* 	     NUMMER(p,i),it,d1.x,d1.y,d1.z); */
/*       printf("%10.6f\t%10.6f\t%10.6f\t",  */
/* 	     KRAFT(p,i,X),KRAFT(p,i,Y),KRAFT(p,i,Z)); */
/*       printf("%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t", */
/* 	     pi.x,pi.y,pi.z,DP_E_IND(p,i,X),DP_E_IND(p,i,Y),DP_E_IND(p,i,Z)); */
/*       printf("%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\n", */
/* 		   DP_P_STAT(p,i,X), DP_P_STAT(p,i,Y),DP_P_STAT(p,i,Z), */
/* 		   DP_E_STAT(p,i,X), DP_E_STAT(p,i,Y),DP_E_STAT(p,i,Z)); */

#endif
    }
    /* Shift Field history */
    dp_E_shift = p->dp_E_old_1;
    p->dp_E_old_1 = p->dp_E_ind;
    p->dp_E_ind = dp_E_shift;
    
  }
  
  if (is_short) fprintf(stderr, "\n Short distance, dipole, step %d!\n",\
			steps);
  dp_E_calc++; 			/* increase field calc counter */
#endif /* DIPOLE */

  /* EWALD is only partially parallelized */
#ifdef EWALD 
#ifdef MPI
  if ((ew_nmax >= 0) || (ew_kcut > 0)) 
    error("option EWALD is only partially parallelized");
#endif
  do_forces_ewald(steps);
#endif 

#ifdef MPI
  /* sum up results of different CPUs */
  tmpvec1[0]     = tot_pot_energy;
  tmpvec1[1]     = virial;
  tmpvec1[2]     = vir_xx;
  tmpvec1[3]     = vir_yy;
  tmpvec1[4]     = vir_zz;
  tmpvec1[5]     = vir_xy;
  tmpvec1[6]     = vir_yz;
  tmpvec1[7]     = vir_zx;
  MPI_Allreduce( tmpvec1, tmpvec2, 8, REAL, MPI_SUM, cpugrid); 
  tot_pot_energy = tmpvec2[0];
  virial         = tmpvec2[1];
  vir_xx         = tmpvec2[2];
  vir_yy         = tmpvec2[3];
  vir_zz         = tmpvec2[4];
  vir_xy         = tmpvec2[5];
  vir_yz         = tmpvec2[6];
  vir_zx         = tmpvec2[7];
#endif

  /* add forces back to original cells/cpus */
  send_forces(add_forces,pack_forces,unpack_forces);

}

/******************************************************************************
*
*  check_nblist
*
******************************************************************************/

void check_nblist()
{
  real   r2, max1=0.0, max2;
  vektor d;
  int    k;

  /* compare with reference positions */
  for (k=0; k<NCELLS; k++) {
    int  i;
    cell *p = CELLPTR(k);
#ifdef ia64
#pragma ivdep,swp
#endif
    for (i=0; i<p->n; i++) {
      d.x = ORT(p,i,X) - NBL_POS(p,i,X);
      d.y = ORT(p,i,Y) - NBL_POS(p,i,Y);
#ifndef TWOD
      d.z = ORT(p,i,Z) - NBL_POS(p,i,Z);
#endif
      r2 = SPROD(d,d);
      if (r2 > max1) max1 = r2;
    }
  }

#ifdef MPI
  MPI_Allreduce( &max1, &max2, 1, REAL, MPI_MAX, cpugrid); 
#else
  max2 = max1;
#endif
  if (max2 > SQR(0.5*nbl_margin)) have_valid_nbl = 0;
}


#ifdef SM

/******************************************************************************
*
*  calc_sm_pot
*
******************************************************************************/

void calc_sm_pot()
{
  int i, k, n=0, m, is_short=0, inc=ntypes*ntypes;

  /* fill the buffer cells */
  send_cells(copy_sm_charge,pack_sm_charge,unpack_sm_charge);

  /* clear per atom accumulation variables, also in buffer cells */
  for (k=0; k<nallcells; k++) {
    cell *p = cell_array + k;
    for (i=0; i<p->n; i++) {
      V_SM(p,i) = 0.0;
    }
  }

  /* pair interactions - for all atoms */
  n=0;
  for (k=0; k<ncells; k++) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; i++) {

      vektor d1;
      real   phi, ch_i, pot = 0.0;
      int    p_typ;

      d1.x  = ORT(p,i,X);
      d1.y  = ORT(p,i,Y);
      d1.z  = ORT(p,i,Z);
      ch_i  = Q_SM(p,i);
      p_typ = SORTE(p,i);

      /* loop over neighbors */
      for (m=tl[n]; m<tl[n+1]; m++) {

        vektor d;
        real   r2, ch_j;
        int    c  = cl_num[ tb[m] ];
        int    j  = tb[m] - cl_off[c];
        cell   *q = cell_array + c;
        int    col2;

        d.x  = ORT(q,j,X) - d1.x;
        d.y  = ORT(q,j,Y) - d1.y;
        d.z  = ORT(q,j,Z) - d1.z;
        r2   = SPROD(d,d);
        ch_j = Q_SM(q,j);

        if (SQR(ch_i * ch_j) > 0.0) {
          if (r2 < ew_r2_cut) {	
            /* Coulomb potential is in column 0, contains already coul_eng */
            int incr = coul_table.ncols;
            VAL_FUNC(phi, coul_table, 0, incr, r2, is_short);
            pot        += phi * ch_j;
            V_SM(q,j)  += phi * ch_i;
          }
          col2 = p_typ * ntypes + SORTE(q,j);
          if (r2 < cr_pot_tab.end[col2]) {
            VAL_FUNC(phi, cr_pot_tab, col2, inc, r2, is_short);
            pot        += phi * ch_j * coul_eng;
            V_SM(q,j)  += phi * ch_i * coul_eng;
          }
        }
      }
      V_SM(p,i) += pot;
      n++;
    }
  }
  if (is_short) fprintf(stderr,"Short distance in calc_sm_pot!\n");

  /* contribution of coulomb self energy */
  for (k=0; k<ncells; k++) {
    real tmp = ew_vorf * 2;
    cell *p  = CELLPTR(k);
    for (i=0; i<p->n; i++) {
      int typ = SORTE(p,i);
      V_SM(p,i) -= (tmp*coul_eng - sm_J_0[typ]) * Q_SM(p,i);
    }
  }

  /* add SM potentials back to original cells/cpus */
  send_forces(add_sm_pot,pack_sm_pot,unpack_add_sm_pot);
}

/******************************************************************************
*
*  calc_sm_chi
*
******************************************************************************/

void calc_sm_chi()
{
  int i, k, n=0, m, is_short=0;

  if (0==have_valid_nbl) {
#ifdef MPI
    /* check message buffer size */
    if (0 == nbl_count % BUFSTEP) setup_buffers();
#endif
    /* update cell decomposition */
    fix_cells();
  }

  /* make new neighbor lists */
  if (0==have_valid_nbl) make_nblist();

  /* clear per atom accumulation variables, also in buffer cells */
  for (k=0; k<nallcells; k++) {
    cell *p = cell_array + k;
    for (i=0; i<p->n; i++) {
      CHI_SM(p,i) = 0.0;
    }
  }

  /* pair interactions - for all atoms */
  n=0;
  for (k=0; k<ncells; k++) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; i++) {

      vektor d1;
      real   z_sm_p, ch_i, chi_tmp = 0.0;
      int    p_typ;

      p_typ  = SORTE(p,i);
      z_sm_p = sm_Z[p_typ]*coul_eng;

      d1.x = ORT(p,i,X);
      d1.y = ORT(p,i,Y);
      d1.z = ORT(p,i,Z);
      ch_i = CHARGE(p,i);

      /* loop over neighbors */
      for (m=tl[n]; m<tl[n+1]; m++) {

        vektor d;
        real   r2, z_sm_q, ch_j;
        int    c  = cl_num[ tb[m] ];
        int    j  = tb[m] - cl_off[c];
        cell   *q = cell_array + c;
        int    q_typ, col1, col2, inc=ntypes*ntypes;;

        q_typ = SORTE(q,j);
        z_sm_q = sm_Z[q_typ]*coul_eng;
        col1  = q_typ * ntypes + p_typ;
        col2  = p_typ * ntypes + q_typ;

        d.x  = ORT(q,j,X) - d1.x;
        d.y  = ORT(q,j,Y) - d1.y;
        d.z  = ORT(q,j,Z) - d1.z;
        r2   = SPROD(d,d);
        ch_j = CHARGE(q,j);

        /* compute electronegativity */
        if (SQR(ch_i * ch_j) > 0.0) {
          real na_pot_p, na_pot_q, cr_pot; 
          na_pot_p = na_pot_q = cr_pot = 0.0; 
          if (r2 < na_pot_tab.end[col2]) {
            VAL_FUNC(na_pot_p, na_pot_tab, col2, inc, r2, is_short);
          }
          if (r2 < na_pot_tab.end[col1]) {
            VAL_FUNC(na_pot_q, na_pot_tab, col1, inc, r2, is_short);
          }
          if (r2 < cr_pot_tab.end[col2]) {
            VAL_FUNC(cr_pot, cr_pot_tab, col2, inc, r2, is_short);
          }
          chi_tmp     += z_sm_q * (na_pot_p - cr_pot);
          CHI_SM(q,j) += z_sm_p * (na_pot_q - cr_pot);
        }
      }
      CHI_SM(p,i) += chi_tmp;
      n++;
    }
  }
  if (is_short) fprintf(stderr,"Short distance in calc_sm_chi!\n");

  /* add chi back to original cells/cpus */
  send_forces(add_sm_chi,pack_sm_chi,unpack_add_sm_chi);

  /* add sm_chi_0 */
  for (k=0; k<ncells; k++) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; i++) {
      int t = SORTE(p,i);
      CHI_SM(p,i) += sm_chi_0[t];
    }
  }

}

#endif  /* SM */
