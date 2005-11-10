/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2005 Institute for Theoretical and Applied Physics,
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

int  *tl=NULL, *tb=NULL, *cl_off=NULL, *cl_num=NULL;

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
  return tn;
}

/******************************************************************************
*
*  make_nblist
*
******************************************************************************/

void make_nblist(void)
{
  static int at_max=0, nb_max=0, pa_max=0;
  int  c, i, k, n, tn, at, cc;

#ifdef MPI
  if (0 == nbl_count % BUFSTEP) setup_buffers();
#endif

  /* update cell decomposition */
  do_boundaries();
  fix_cells();

  /* fill the buffer cells */
  send_cells(copy_cell,pack_cell,unpack_cell);

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

  /* allocate cl_off */
  if (nbl_count==0) {
    cl_off = (int *) malloc( nallcells * sizeof(int) );
    if (cl_off==NULL) error("cannot allocate neighbor table");
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
  if (nbl_count==0) {
    nb_max = (int  ) (nbl_size * estimate_nblist_size());
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
      if (tn > nb_max-2*pa_max) 
        error("neighbor table full - increase nbl_size");
    }
  }
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

  /* fill the buffer cells */
  send_cells(copy_cell,pack_cell,unpack_cell);

  /* clear global accumulation variables */
  tot_pot_energy = 0.0;
  virial = 0.0;
  vir_xx = 0.0;
  vir_yy = 0.0;
  vir_xy = 0.0;
  vir_zz = 0.0;
  vir_yz = 0.0;
  vir_zx = 0.0;
  nfc++;

  /* clear per atom accumulation variables */
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
#ifdef STRESS_TENS
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
#ifdef ORDPAR
      NBANZ(p,i) = 0;
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
#ifdef NVX
      HEATCOND(p,i) = 0.0;
#endif     
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
#ifdef TWOD
      vektor d1, ff = {0.0,0.0};
#else
      vektor d1, ff = {0.0,0.0,0.0};
#endif
      real   ee = 0.0, hc = 0.0;
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

        /* compute pair interactions */
#if defined(PAIR) || defined(KEATING) || defined(STIWEB)
        /* PAIR, KEATING, and STIWEB are mutually exclusive */
#if defined(PAIR)
        if (r2 <= pair_pot.end[col]) {
#ifdef LINPOT
          PAIR_INT_LIN(pot, grad, pair_pot_lin, col, inc, r2, is_short)
#else
          PAIR_INT(pot, grad, pair_pot, col, inc, r2, is_short)
#endif
#elif defined(KEATING)
        if (r2 < keat_r2_cut[it][jt]) {
          PAIR_INT_KEATING(pot, grad, it, jt, r2)
#elif defined(STIWEB)  
        if (r2 < sw_2_a1[it][jt]) {
	  PAIR_INT_STIWEB(pot, grad, it, jt, r2)
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

#ifndef MONOLJ
#ifdef EAM2
          pot *= 0.5;   /* avoid double counting */
#endif
#ifdef ORDPAR
          if (r2 < op_r2_cut[it][jt]) {
            ee          += op_weight[it][jt] * pot;
            POTENG(q,j) += op_weight[jt][it] * pot;
            nb++;
            NBANZ(q,j)++;
          }
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
#ifdef NVX
          hc            += pot - r2 * grad;
          HEATCOND(q,j) += pot - r2 * grad;
#endif
        }

#endif /* PAIR || KEATING || STIWEB */

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
          col2 = jt * ntypes + it;
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
#ifdef ORDPAR
      NBANZ(p,i)    += nb;
#endif
#ifdef NVX
      HEATCOND(p,i) += hc;
#endif
      n++;
    }
  }
  if (is_short) printf("short distance!\n");

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

  /* second force loop for covalent systems */
  for (k=0; k<ncells; ++k) {
    do_forces2(cell_array + cnbrs[k].np,
               &tot_pot_energy, &virial, &vir_xx, &vir_yy, &vir_zz,
                                         &vir_yz, &vir_zx, &vir_xy);
  }

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
  if (is_short) fprintf(stderr, "\n Short distance!\n");

#endif /* EAM2 */

  /* EWALD is not parallelized */
#if defined(EWALD) && !defined(MPI) 
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
  if (max2 > SQR(0.5*nbl_margin)) make_nblist();
}
