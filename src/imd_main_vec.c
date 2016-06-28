
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2006 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_main_vec_3d.c -- force loop, vector specific part, three dimensions
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"
#include "potaccess.h"

#ifndef SX
#define restrict
#endif

/* for VEC2, the following are no longer used */
#define N0 140000     /* at least atoms.n */
#define MC 150        /* maximum number of particles in cell */
#define BS 100
#define N1 MC*14
#define N2 MC*4

/******************************************************************************
*
*  calc_forces, using Verlet neighbor lists
*
*  We use buffer cells on the surface of the cell array, even for 
*  the serial version. Atoms in the buffer cells have positions with
*  periodic boundary conditions applied, where applicable.
*
*  We can use two versions of calc_forces:
*
*     VEC2:  vectorizing LLC algorithm as in Grest, Dünwald, Kremer
*     VEC3:  keep one atom fixed, run over neighbors
*
******************************************************************************/

#if !defined(VEC2) && !defined(VEC3)
#define VEC2
#endif

#ifdef VEC2

#ifdef MONO
#define COL(i)  0
#define COL2(i) 0
#else
#define COL(i)  col[i]
#define COL2(i) col2[i]
#endif

/* int li[N0*N2], lj[N0*N2], col[N0*N2], lbeg[N1], lend[N1], ltot; */ 
long *li=NULL, *lj=NULL, *col=NULL, *col2=NULL, *lbeg=NULL, *lend=NULL, ltot; 
long *pair_filt=NULL;
long *eam_filt0=NULL, *eam_filt1=NULL, *eam_filt2=NULL, *eam_filt3=NULL;
real *rr=NULL, *grad=NULL;
vektor *dd=NULL, *ff=NULL;

/******************************************************************************
*
*  make_nblist
*
*  This version implements the vectorizing LLC algorithm from
*  Grest, Dünweg & Kremer, Comp. Phys. Comm. 55, 269-285 (1989).
*
******************************************************************************/

void make_nblist()
{
  static long max_atoms=0, max_pairs=0, max_cell_max=0, max_sub_len=0;
  static long *li1 = NULL, *lj1 = NULL;
  long   i, j, k, l, m, n, len1, len2, max_cell, max_len=0;
  /* int    li1[N0], lj1[N0]; */

  /* update reference positions */
  for (i=0; i<DIM*atoms.n; i++) atoms.nbl_pos[i] = atoms.ort[i];

  /* get maximum cell occupancy */
  max_cell = 0;
  for (i=0; i<nallcells; i++) {
    n = (cell_array+i)->n;
    if (max_cell < n) max_cell = n;
  }

  /* allocate neighbor table arrays */
  if (max_atoms < atoms.n) {
    free(li1);
    free(lj1);
    max_atoms = (long) (atoms.n * 1.2);
    li1 = (long *) malloc( max_atoms * sizeof(long) );
    lj1 = (long *) malloc( max_atoms * sizeof(long) );
    if ((li1==NULL) || (lj1==NULL))
      error("cannot allocate temporary pair arrays");
  }
  if (max_pairs < 4 * max_cell * atoms.n) {
    free(li);
    free(lj);
    free(col);
    max_pairs = (long) (4 * max_cell * atoms.n * 1.2);
    li   = (long *) malloc( max_pairs * sizeof(long) );
    lj   = (long *) malloc( max_pairs * sizeof(long) );
    col  = (long *) malloc( max_pairs * sizeof(long) );
    if ((li==NULL) || (lj==NULL) || (col==NULL))
      error("cannot allocate temporary pair arrays");
#if defined(EAM2) || defined(NNBR)
    free(col2);
    col2 = (long *) malloc( max_pairs * sizeof(long) );
    if (col2==NULL)
      error("cannot allocate temporary pair arrays");
#endif
  }
  if (max_cell_max < max_cell) {
    free(lbeg);
    free(lend);
    max_cell_max = (long) (max_cell * 1.2);
    lbeg = (long *) malloc( 14 * max_cell_max * sizeof(long) );
    lend = (long *) malloc( 14 * max_cell_max * sizeof(long) );
    if ((lbeg==NULL) || (lend==NULL))
      error("cannot allocate temporary arrays");
  }

  len2=0;
  ltot=0;
  for (n=0; n<14; n++) {
    for (k=0; k<max_cell; k++) {

      lbeg[ltot] = len2;

#ifdef FTRACE
      ftrace_region_begin("nbl_get_indices");
#endif

      /* make list of index pairs */
      len1 = 0;
      for (i=0; i<max_cell; i++) {
        j = i + k;
        if (j>=max_cell) j -= max_cell;
        if ((n==0) && (j<=i)) continue;  /* twice the same cell */
        for (m=0; m<ncells; m++) {
          minicell *ci, *cj;
          long cjn = cnbrs[m].nq[n];
          if (cjn==-1) continue;
          ci = cell_array + cnbrs[m].np;
          cj = cell_array + cjn;
          if ( (i >= ci->n) || (j >= cj->n) ) continue;
          li1[len1] = ci->ind[i];
          lj1[len1] = cj->ind[j];
          len1++;
	}
      }
      if (len1 > max_atoms) printf("max_atoms %d %d\n", max_atoms, len1 );

#ifdef FTRACE
      ftrace_region_end  ("nbl_get_indices");
      ftrace_region_begin("nbl_calc_distances");
#endif

      /* reduce pair list */
      for (l=0; l<len1; l++) {
        real * restrict d1, * restrict d2, dx, dy, dz, r2;
        d1 = atoms.ort + DIM * li1[l];
        d2 = atoms.ort + DIM * lj1[l];
        dx = d2[0]-d1[0];
        dy = d2[1]-d1[1];
        dz = d2[2]-d1[2];
        r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < cellsz) {
          li[len2] = li1[l];
          lj[len2] = lj1[l];
          len2++;
        }
      }

#ifdef FTRACE
      ftrace_region_end("nbl_calc_distances");
#endif

      lend[ltot]=len2;
      max_len = MAX( max_len, lend[ltot] - lbeg[ltot] );
      if (lend[ltot] > lbeg[ltot]) ltot++;
    }
  }
  if (max_pairs < len2) error("neighbor table overflow"); 

  /* allocate arrays for filter flags */
  if (max_len > max_sub_len) {
    max_sub_len = (long) (1.2 * max_len);
    free(dd);
    free(ff);
    free(rr);
    free(grad);
    free(pair_filt);
#ifdef EAM2
    free(eam_filt0);
    free(eam_filt1);
    free(eam_filt2);
#ifndef MONO
    free(eam_filt3);
#endif
#endif
    rr   = (real   *) malloc( max_sub_len * sizeof(real  ) );
    grad = (real   *) malloc( max_sub_len * sizeof(real  ) );
    dd   = (vektor *) malloc( max_sub_len * sizeof(vektor) );
    ff   = (vektor *) malloc( max_sub_len * sizeof(vektor) );
    pair_filt = (long *) malloc( max_sub_len * sizeof(long) );
    if ((dd==NULL) || (rr==NULL) || (grad==NULL) || 
        (pair_filt==NULL) || (ff==NULL))
      error("cannot allocate temporary pair arrays");
#ifdef EAM2
    eam_filt0 = (long *) malloc( max_sub_len * sizeof(long) );
    eam_filt1 = (long *) malloc( max_sub_len * sizeof(long) );
    eam_filt2 = (long *) malloc( max_sub_len * sizeof(long) );
#ifdef MONO
    eam_filt3 = eam_filt0;
#else
    eam_filt3 = (long *) malloc( max_sub_len * sizeof(long) );
#endif
    if ((eam_filt0==NULL) || (eam_filt1==NULL) || (eam_filt2==NULL) ||
        (eam_filt3==NULL))
      error("cannot allocate temporary pair arrays");
#endif
  }

#ifndef MONO

#ifdef FTRACE
  ftrace_region_begin("nbl_make_col");
#endif

  /* precompute potential column for each atom pair */
  for (l=0; l<lend[ltot-1]; l++) {
    col [l] = SORTE(&atoms,li[l]) * ntypes + SORTE(&atoms,lj[l]);
#if defined(EAM2) || defined(NNBR)
    col2[l] = SORTE(&atoms,li[l]) + ntypes * SORTE(&atoms,lj[l]);
#endif
  }

#ifdef FTRACE
  ftrace_region_end("nbl_make_col");
#endif

#endif /* not MONO */

  have_valid_nbl = 1;
  nbl_count++;
}

/******************************************************************************
*
*  calc_forces
*
*  This version implements the vectorizing LLC algorithm from
*  Grest, Dünweg & Kremer, Comp. Phys. Comm. 55, 269-285 (1989).
*
******************************************************************************/

void calc_forces(int steps)
{
  static int len=0;
  long   k, m, is_short=0, idummy=0;
  real   tmpvec1[5], tmpvec2[5] = {0.0,0.0,0.0,0.0,0.0};
  static real *kraft1 = NULL, *pot_eng1 = NULL, *eam_rho1 = NULL;
  static shortint *nbanz1 = NULL;
  static sym_tensor *presstens1 = NULL;

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
  virial = 0.0;
  vir_xx = 0.0;
  vir_yy = 0.0;
  vir_zz = 0.0;
  nfc++;

  /* (re)allocate temporary force and energy arrays */
  if (len < atoms.n_buf) {
    free(pot_eng1);
    free(eam_rho1);
    free(kraft1);
    len = (int) (atoms.n_buf * 1.2);
    pot_eng1 = (real *) malloc(       len * sizeof(real) );
    kraft1   = (real *) malloc( DIM * len * sizeof(real) );
    if ((pot_eng1 == NULL) || (kraft1 == NULL))
      error("cannot allocate temporary force and energy arrays");
#ifdef EAM2
    eam_rho1 = (real *) malloc( len * sizeof(real) );
    if (eam_rho1 == NULL)
      error("cannot allocate temporary host electron density array");
#endif
#ifdef NNBR
    nbanz1 = (shortint *) malloc( len * sizeof(shortint) );
    if (nbanz1 == NULL)
      error("cannot allocate temporary coordination number array");
#endif
#ifdef STRESS_TENS
    presstens1 = (sym_tensor *) malloc( len * sizeof(sym_tensor) );
    if (presstens1 == NULL)
      error("cannot allocate temporary pressure tensor array");
#endif
  }

  /* clear per atom accumulation variables */
  for (k=0; k<DIM*atoms.n_buf; k++) atoms.kraft  [k] = 0.0;
  for (k=0; k<    atoms.n_buf; k++) atoms.pot_eng[k] = 0.0;
  for (k=0; k<DIM*atoms.n_buf; k++)      kraft1  [k] = 0.0;
  for (k=0; k<    atoms.n_buf; k++)      pot_eng1[k] = 0.0;
#ifdef EAM2
  for (k=0; k<    atoms.n_buf; k++) atoms.eam_rho[k] = 0.0;
  for (k=0; k<    atoms.n_buf; k++)      eam_rho1[k] = 0.0;
#endif
#ifdef NNBR
  for (k=0; k<    atoms.n_buf; k++) atoms.nbanz[k]   = 0;
  for (k=0; k<    atoms.n_buf; k++)      nbanz1[k]   = 0;
#endif
#ifdef STRESS_TENS
  for (k=0; k<atoms.n_buf; k++) {
    atoms.presstens[k].xx = 0.0;
    atoms.presstens[k].yy = 0.0;
    atoms.presstens[k].zz = 0.0;
    atoms.presstens[k].yz = 0.0;
    atoms.presstens[k].zx = 0.0;
    atoms.presstens[k].xy = 0.0;
  }
  for (k=0; k<atoms.n_buf; k++) {
    presstens1[k].xx = 0.0;
    presstens1[k].yy = 0.0;
    presstens1[k].zz = 0.0;
    presstens1[k].yz = 0.0;
    presstens1[k].zx = 0.0;
    presstens1[k].xy = 0.0;
  }
#endif

  /* loop over independent pair sublists */
  for (k=0; k<ltot; k++) {

    long s, pfmax=0, ef0max=0, ef1max=0, ef2max=0;

#ifdef FTRACE
    ftrace_region_begin("calc_forces1");
#endif

#ifdef SX
#pragma cdir vector,nodep
#endif
    /* compute distances */
    for (m=0; m<lend[k]-lbeg[k]; m++) {

      const long l = lbeg[k] + m;
      real *d1, *d2;
      d1 = atoms.ort + DIM * li[l];
      d2 = atoms.ort + DIM * lj[l];
      dd[m].x = d2[0]-d1[0];
      dd[m].y = d2[1]-d1[1];
      dd[m].z = d2[2]-d1[2];
      rr[m] = SQR(dd[m].x) + SQR(dd[m].y) + SQR(dd[m].z);

      /* filter interacting pairs */
      if (rr[m] <= pair_pot.end[COL(l)]) {
        pair_filt[pfmax++] = m;
      }

#ifdef EAM2
      /* filter interacting pairs for EAM */
      if (COL(l)==COL2(l)) {
        if (rr[m] < rho_h_tab.end[COL(l) ]) {
          eam_filt0[ef0max++] = m;
        }
      }
#ifndef MONO
      else {
        if (rr[m] < rho_h_tab.end[COL(l) ]) {
          eam_filt1[ef1max++] = m;
        }
        if (rr[m] < rho_h_tab.end[COL2(l)]) {
          eam_filt2[ef2max++] = m;
        }
      }
#endif /* not MONO */
#endif /* EAM2 */

    }

#ifdef FTRACE
    ftrace_region_end("calc_forces1");
    ftrace_region_begin("calc_forces2");
#endif

    /* compute pair forces */
#ifdef SX
#pragma cdir vector,nodep
#endif
    for (s=0; s<pfmax; s++) {

      const long m = pair_filt[s], l = lbeg[k] + m;
      const long i = li[l], j = lj[l];
      const long inc = ntypes * ntypes;
      real  ee, grad, r2c;
      vektor f;

#ifdef LJ
      PAIR_INT_LJ_VEC(ee, grad, COL(l), rr[m])
#else
      r2c = MAX(0.0, rr[m] - pair_pot.begin[COL(l)]);
      PAIR_INT_VEC(ee, grad, pair_pot, COL(l), inc, r2c)
#endif

      /* store force in temporary variables */
      f.x             = dd[m].x * grad;
      f.y             = dd[m].y * grad;
      f.z             = dd[m].z * grad;
#ifdef P_AXIAL
      vir_xx         -= dd[m].x * f.x;
      vir_yy         -= dd[m].y * f.y;
      vir_zz         -= dd[m].z * f.z;
#else
      virial         -= rr[m]   * grad;
#endif
      tot_pot_energy += ee;
#ifdef EAM2
      ee *= 0.5;   /* avoid double counting */
#endif

      /* store force to first particle */
      kraft1 X(i) += f.x;
      kraft1 Y(i) += f.y;
      kraft1 Z(i) += f.z;
      pot_eng1[i] += ee;

      /* store force to second particle */
      atoms.kraft X(j) -= f.x;
      atoms.kraft Y(j) -= f.y;
      atoms.kraft Z(j) -= f.z;
      atoms.pot_eng[j] += ee;

#ifdef STRESS_TENS
      /* we need to keep a copy; avoid double counting of the virial */
      /* if (do_press_calc) {...} would prevent vectorisation */
      ff[m].x = 0.5 * f.x;
      ff[m].y = 0.5 * f.y;
      ff[m].z = 0.5 * f.z;
#endif
    }

#ifdef FTRACE
    ftrace_region_end("calc_forces2");
#endif

#ifdef NNBR
    for (s=0; s<pfmax; s++) {
      const long m = pair_filt[s], l = lbeg[k] + m;
      const long i = li[l], j = lj[l];
      if (rr[m] < nb_r2_cut[COL(l) ]) nbanz1[i]++;
      if (rr[m] < nb_r2_cut[COL2(l)]) atoms.nbanz[j]++;
    }
#endif

#ifdef STRESS_TENS

#ifdef FTRACE
    ftrace_region_begin("calc_forces3");
#endif

    if (do_press_calc) {
#ifdef SX
#pragma cdir vector,nodep
#endif
      for (s=0; s<pfmax; s++) {

        const long m = pair_filt[s], l = lbeg[k] + m;
        const long i = li[l], j = lj[l];

        presstens1     [i].xx -= dd[m].x * ff[m].x;
        atoms.presstens[j].xx -= dd[m].x * ff[m].x;
        presstens1     [i].yy -= dd[m].y * ff[m].y;
        atoms.presstens[j].yy -= dd[m].y * ff[m].y;
        presstens1     [i].zz -= dd[m].z * ff[m].z;
        atoms.presstens[j].zz -= dd[m].z * ff[m].z;
        presstens1     [i].yz -= dd[m].y * ff[m].z;
        atoms.presstens[j].yz -= dd[m].y * ff[m].z;
        presstens1     [i].zx -= dd[m].z * ff[m].x;
        atoms.presstens[j].zx -= dd[m].z * ff[m].x;
        presstens1     [i].xy -= dd[m].x * ff[m].y;
        atoms.presstens[j].xy -= dd[m].x * ff[m].y;
      }
    }

#ifdef FTRACE
    ftrace_region_end("calc_forces3");
#endif

#endif

#ifdef EAM2

    /* compute host electron density */

#ifdef FTRACE
    ftrace_region_begin("eam_rho");
#endif

#ifdef SX
#pragma cdir vector,nodep
#endif
    for (s=0; s<ef0max; s++) {
      const long m   = eam_filt0[s], l = lbeg[k] + m;
      const long inc = ntypes * ntypes;
      real  r2c, rho_h;
      r2c = MAX(0.0, rr[m] - rho_h_tab.begin[COL(l)]);
      VAL_FUNC_VEC(rho_h, rho_h_tab, COL(l), inc, r2c);
      atoms.eam_rho[ li[l] ] += rho_h; 
      eam_rho1     [ lj[l] ] += rho_h;
    }

#ifndef MONO

#ifdef SX
#pragma cdir vector,nodep
#endif
    for (s=0; s<ef1max; s++) {
      const long m   = eam_filt1[s], l = lbeg[k] + m;
      const long inc = ntypes * ntypes;
      real  r2c, rho_h;
      r2c = MAX(0.0, rr[m] - rho_h_tab.begin[COL(l)]);
      VAL_FUNC_VEC(rho_h, rho_h_tab, COL(l), inc, r2c);
      atoms.eam_rho[ li[l] ] += rho_h; 
    }

#ifdef SX
#pragma cdir vector,nodep
#endif
    for (s=0; s<ef2max; s++) {
      const long m   = eam_filt2[s], l = lbeg[k] + m;
      const long inc = ntypes * ntypes;
      real  r2c, rho_h;
      r2c = MAX(0.0, rr[m] - rho_h_tab.begin[COL2(l)]);
      VAL_FUNC_VEC(rho_h, rho_h_tab, COL2(l), inc, r2c);
      eam_rho1[ lj[l] ] += rho_h; 
    }

#endif /* not MONO */

#ifdef FTRACE
    ftrace_region_end("eam_rho");
#endif

#endif /* EAM2 */

  }
  if (is_short) printf("short distance!\n");

#ifdef EAM2

  /* add up eam_rho on first and second particles */
  for (k=0; k<atoms.n_buf; k++) atoms.eam_rho[k] += eam_rho1[k];

  /* collect host electron density */
  send_forces(add_rho,pack_rho,unpack_add_rho);

#ifdef FTRACE
  ftrace_region_begin("eam_energy");
#endif

  /* compute embedding energy and its derivative */
#ifdef SX
#pragma cdir vector,nodep
#endif
  for (m=0; m<atoms.n; m++) {
    long st = SORTE(&atoms,m);
    real pot, rho = MAX( 0.0, atoms.eam_rho[m] - embed_pot.begin[st] );
    PAIR_INT_VEC( pot, atoms.eam_dF[m], embed_pot, st, ntypes, rho);
    atoms.pot_eng[m] += pot;
    tot_pot_energy   += pot;
  }

#ifdef FTRACE
  ftrace_region_end("eam_energy");
#endif

  /* distribute derivative of embedding energy */
  send_cells(copy_dF,pack_dF,unpack_dF);

  /* EAM forces - loop over independent pair sublists */
  for (k=0; k<ltot; k++) {

    long s, ef0max=0, ef1max=0, ef2max=0, ef3max=0;

#ifdef FTRACE
    ftrace_region_begin("eam_forces1");
#endif

#ifdef SX
#pragma cdir vector,nodep
#endif
    for (m=0; m<lend[k]-lbeg[k]; m++) {

      const long l = lbeg[k] + m;
      real *d1, *d2;
      d1 = atoms.ort + DIM * li[l];
      d2 = atoms.ort + DIM * lj[l];
      dd[m].x = d2[0]-d1[0];
      dd[m].y = d2[1]-d1[1];
      dd[m].z = d2[2]-d1[2];
      rr[m]   = SQR(dd[m].x) + SQR(dd[m].y) + SQR(dd[m].z);
      grad[m] = 0.0;
    }

#ifdef FTRACE
    ftrace_region_end("eam_forces1");
    ftrace_region_begin("eam_forces2");
#endif

    /* filter interacting pairs */
#ifdef MONO

#ifdef SX
#pragma cdir vector,nodep
#endif
    for (m=0; m<lend[k]-lbeg[k]; m++) {
      if (rr[m] < rho_h_tab.end[0]) {
        eam_filt0[ef0max++] = m;
      }
    }
    ef3max = ef0max;

#else  /* not MONO */

#ifdef SX
#pragma cdir vector,nodep
#endif
    for (m=0; m<lend[k]-lbeg[k]; m++) {
      const long l = lbeg[k] + m;
      if (COL(l)==COL2(l)) {
        if (rr[m] < rho_h_tab.end[COL2(l)]) {
          eam_filt0[ef0max++] = m;
        }
      }
      else {
        if (rr[m] < rho_h_tab.end[COL2(l)]) {
          eam_filt1[ef1max++] = m;
        }
        if (rr[m] < rho_h_tab.end[COL(l) ]) {
          eam_filt2[ef2max++] = m;
        }
      }
      if ((rr[m]<rho_h_tab.end[COL(l)]) || (rr[m]<rho_h_tab.end[COL2(l)])) {
        eam_filt3[ef3max++] = m;
      }
    }
#endif /* not MONO */

#ifdef FTRACE
    ftrace_region_end("eam_forces2");
    ftrace_region_begin("eam_forces3");
#endif

    /* compute potential derivatives */
#ifdef SX
#pragma cdir vector,nodep
#endif
    for (s=0; s<ef0max; s++) {
      const long inc = ntypes * ntypes;
      const long m = eam_filt0[s], l = lbeg[k] + m;
      real  r2c, rho_i_str;
      r2c = MAX( 0.0, rr[m] - rho_h_tab.begin[COL2(l)] );
      DERIV_FUNC_VEC(rho_i_str, rho_h_tab, COL2(l), inc, r2c);
      grad[m] = 0.5 * (atoms.eam_dF[li[l]] + atoms.eam_dF[lj[l]]) * rho_i_str;
    }

#ifndef MONO

#ifdef SX
#pragma cdir vector,nodep
#endif
    for (s=0; s<ef1max; s++) {
      const long inc = ntypes * ntypes;
      const long m = eam_filt1[s], l = lbeg[k] + m;
      real  r2c, rho_i_str;
      r2c = MAX( 0.0, rr[m] - rho_h_tab.begin[COL2(l)] );
      DERIV_FUNC_VEC(rho_i_str, rho_h_tab, COL2(l), inc, r2c);
      grad[m] += 0.5 * atoms.eam_dF[lj[l]] * rho_i_str;
    }

#ifdef SX
#pragma cdir vector,nodep
#endif
    for (s=0; s<ef2max; s++) {
      const long inc = ntypes * ntypes;
      const long m = eam_filt2[s], l = lbeg[k] + m;
      real  r2c, rho_j_str;
      r2c = MAX( 0.0, rr[m] - rho_h_tab.begin[COL(l)] );
      DERIV_FUNC_VEC(rho_j_str, rho_h_tab, COL(l), inc, r2c);
      grad[m] += 0.5 * atoms.eam_dF[li[l]] * rho_j_str;
    }

#endif /* not MONO */

#ifdef FTRACE
    ftrace_region_end("eam_forces3");
    ftrace_region_begin("eam_forces4");
#endif

    /* compute EAM forces */
#ifdef SX
#pragma cdir vector,nodep
#endif
    for (s=0; s<ef3max; s++) {

      const long m = eam_filt3[s], l = lbeg[k] + m;
      const long i = li[l], j = lj[l]; 
      vektor f;

      /* store force in temporary variable */
      f.x     = dd[m].x * grad[m];
      f.y     = dd[m].y * grad[m];
      f.z     = dd[m].z * grad[m];
#ifdef P_AXIAL
      vir_xx -= dd[m].x * f.x;
      vir_yy -= dd[m].y * f.y;
      vir_zz -= dd[m].z * f.z;
#else
      virial -= rr[m]   * grad[m];
#endif

      /* store force to first particle */
      kraft1 X(i) += f.x;
      kraft1 Y(i) += f.y;
      kraft1 Z(i) += f.z;

      /* store force to second particle */
      atoms.kraft X(j) -= f.x;
      atoms.kraft Y(j) -= f.y;
      atoms.kraft Z(j) -= f.z;

#ifdef STRESS_TENS
      /* we need to keep a copy; avoid double counting of the virial */
      /* if (do_press_calc) {...} would prevent vectorisation */
      ff[m].x = 0.5 * f.x;
      ff[m].y = 0.5 * f.y;
      ff[m].z = 0.5 * f.z;
#endif

    }

#ifdef FTRACE
    ftrace_region_end("eam_forces4");
#endif

#ifdef STRESS_TENS

#ifdef FTRACE
    ftrace_region_begin("eam_forces5");
#endif

    if (do_press_calc) {
#ifdef SX
#pragma cdir vector,nodep
#endif
      for (s=0; s<ef3max; s++) {

        const long m = eam_filt3[s], l = lbeg[k] + m;
        const long i = li[l], j = lj[l]; 

        presstens1     [i].xx -= dd[m].x * ff[m].x;
        atoms.presstens[j].xx -= dd[m].x * ff[m].x;
        presstens1     [i].yy -= dd[m].y * ff[m].y;
        atoms.presstens[j].yy -= dd[m].y * ff[m].y;
        presstens1     [i].zz -= dd[m].z * ff[m].z;
        atoms.presstens[j].zz -= dd[m].z * ff[m].z;
        presstens1     [i].yz -= dd[m].y * ff[m].z;
        atoms.presstens[j].yz -= dd[m].y * ff[m].z;
        presstens1     [i].zx -= dd[m].z * ff[m].x;
        atoms.presstens[j].zx -= dd[m].z * ff[m].x;
        presstens1     [i].xy -= dd[m].x * ff[m].y;
        atoms.presstens[j].xy -= dd[m].x * ff[m].y;
      }
    }

#ifdef FTRACE
    ftrace_region_end("eam_forces5");
#endif

#endif  /* STRESS_TENS */

  }
  if (is_short) printf("short distance!\n");

#endif  /* EAM2 */

#ifdef FTRACE
  ftrace_region_begin("store_forces");
#endif

  /* add up forces on first and second particles */
  for (k=0; k<DIM*atoms.n_buf; k++) atoms.kraft  [k] += kraft1  [k];
  for (k=0; k<    atoms.n_buf; k++) atoms.pot_eng[k] += pot_eng1[k];
#ifdef NNBR
  for (k=0; k<    atoms.n_buf; k++) atoms.nbanz[k]   += nbanz1[k];
#endif
#ifdef STRESS_TENS
  for (k=0; k<atoms.n_buf; k++) {
    atoms.presstens[k].xx += presstens1[k].xx;
    atoms.presstens[k].yy += presstens1[k].yy;
    atoms.presstens[k].zz += presstens1[k].zz;
    atoms.presstens[k].yz += presstens1[k].yz;
    atoms.presstens[k].zx += presstens1[k].zx;
    atoms.presstens[k].xy += presstens1[k].xy;
  }
#endif

#ifdef FTRACE
    ftrace_region_end("store_forces");
#endif

#ifdef MPI
  /* sum up results of different CPUs */
  tmpvec1[0] = tot_pot_energy;
  tmpvec1[1] = virial;
  tmpvec1[2] = vir_xx;
  tmpvec1[3] = vir_yy;
  tmpvec1[4] = vir_zz;
  MPI_Allreduce( tmpvec1, tmpvec2, 5, REAL, MPI_SUM, cpugrid); 
  tot_pot_energy = tmpvec2[0];
  virial         = tmpvec2[1];
  vir_xx         = tmpvec2[2];
  vir_yy         = tmpvec2[3];
  vir_zz         = tmpvec2[4];
#endif

  /* add forces back to original cells/cpus */
  send_forces(add_forces,pack_forces,unpack_forces);

}

#endif /* VEC2 */

#ifdef VEC3

#ifdef MONO
#define COL 0
#else
#define COL col
#endif

int tab[N2*N0], tlen[N0];

/******************************************************************************
*
*  make_nblist
*
*  This version keeps one particle fixed, and treats all its neighbors.
*
******************************************************************************/

void make_nblist()
{
  int  c, i, tb[N1];

  /* update reference positions */
  for (i=0; i<DIM*atoms.n; i++) atoms.nbl_pos[i] = atoms.ort[i];

  /* for all cells */
  for (c=0; c<ncells; c++) {

    minicell *p;
    p = cell_array + cnbrs[c].np;

    /* for each atom in cell */
    for (i=0; i<p->n; i++) {

      real *d1;
      int  *tt, tl, tn, ii, jj, j, k, m;
      minicell *q;

      /* indices of particles from the same cell */
      tt = tb - (i+1);
      tl =    - (i+1);
      for (j = i+1; j < p->n; j++) tt[j] = p->ind[j];
      tt += p->n;
      tl += p->n;

      /* indices of particles from neighbor cells */
      for (m=1; m<14; m++) {
        if (cnbrs[c].nq[m]<0) continue;
        q = cell_array + cnbrs[c].nq[m];
        for (j = 0; j < q->n; j++) tt[j] = q->ind[j];
        tt += q->n;
        tl += q->n;
      }

      /* narrow neighbor table */
      tn = 0;
      ii = p->ind[i];
      tt = tab + N2*ii;
      d1 = atoms.ort + DIM * ii;
#ifdef SX
#pragma cdir vector,nodep
#endif
      for (k=0; k<tl; k++) {
        real *d2, dx, dy, dz, r2;
        d2 = atoms.ort + DIM * tb[k];
        dx = d2[0]-d1[0];
        dy = d2[1]-d1[1];
        dz = d2[2]-d1[2];
        r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < cellsz) {
          tt[tn] = tb[k];          
          tn++;
        }
      }
      tlen[ii] = tn;
    }
  }
  have_valid_nbl = 1;
  nbl_count++;
}

/******************************************************************************
*
*  calc_forces
*
*  This version keeps one particle fixed, and treats all its neighbors.
*
******************************************************************************/

void calc_forces(int steps)
{
  int    i, b, k, is_short=0, nblocks, blocksize, flag[BS*N2];
  int    it[BS*N2], lj[BS*N2], lbeg[BS], lend[BS];
  vektor dd[BS*N2], ff[BS*N2], fi[N0];
  real   rr[BS*N2], ee[BS*N2], ei[N0];
  real   tmpvec1[2], tmpvec2[2] = {0.0, 0.0};

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
  virial = 0.0;

  /* clear per atom accumulation variables */
  for (k=0; k<DIM*atoms.n_buf; k++) atoms.kraft  [k] = 0.0;
  for (k=0; k<    atoms.n_buf; k++) atoms.pot_eng[k] = 0.0;
  /* clear temporary storage */
  for (k=0; k<atoms.n; k++) {
    fi[k].x          = 0.0;
    fi[k].y          = 0.0;
    fi[k].z          = 0.0;
    ei[k]            = 0.0;
  }

  /* to save memory, we treat atoms in blocks */
  nblocks   = atoms.n / BS + 1;
  blocksize = atoms.n / nblocks + 1;

  for (b=0; b<nblocks; b++) {

    int i, l, k, imin, imax, len;

    imin = b*blocksize;
    imax = MIN(atoms.n, imin+blocksize);

#ifdef FTRACE
    ftrace_region_begin("calc_distances");
#endif

    /* compute distances */
    len=0;
    for (i=imin; i<imax; i++) {

      int  ityp, j;
      real d1[3], *d2, dx, dy, dz, r2;

      lbeg[i-imin] = len;
      d1[0] = atoms.ort X(i);
      d1[1] = atoms.ort Y(i);
      d1[2] = atoms.ort Z(i);
      ityp  = SORTE(&atoms,i);
#ifdef SX
#pragma cdir vector,loopcnt=256
#endif
      for (k=0; k<tlen[i]; k++) {
        j  = tab[N2*i+k];
        d2 = atoms.ort + DIM * j;
        dd[len].x = d2[0]-d1[0];
        dd[len].y = d2[1]-d1[1];
        dd[len].z = d2[2]-d1[2];
        rr[len] = dd[len].x*dd[len].x+dd[len].y*dd[len].y+dd[len].z*dd[len].z;
        lj[len] = j;
#ifndef MONO
        it[len] = ityp;
#endif
        len++;
      }
      lend[i-imin] = len;
    }

#ifdef FTRACE
    ftrace_region_end  ("calc_distances");
    ftrace_region_begin("calc_forces");
#endif

    /* compute forces */
    for (l=0; l<len; l++) {

      real grad;
      int  col, inc = ntypes * ntypes; 

      flag[l] = 0;
#ifndef MONO
      col = it[l] * ntypes + SORTE(&atoms,lj[l]);
#endif

      /* compute pair interactions */
      if (rr[l] <= pair_pot.end[COL]) {

        /* beware: we must not use k as index in the macro's arguments! */
#ifdef MULTIPOT
        int pp = l % N_POT_TAB;
        PAIR_INT(ee[l], grad, pair_pot_ar[pp], COL, inc, rr[l], is_short)
#else
#ifdef LINPOT
        PAIR_INT_LIN(ee[l], grad, pair_pot_lin, COL, inc, rr[l], is_short)
#else
#ifdef LJ
        PAIR_INT_LJ_VEC(ee[l], grad, COL, rr[l])
#else
        PAIR_INT(ee[l], grad, pair_pot, COL, inc, rr[l], is_short)
#endif
#endif
#endif
        /* store force in temporary variables */
        ff[l].x = dd[l].x * grad;
        ff[l].y = dd[l].y * grad;
        ff[l].z = dd[l].z * grad;
        virial -= rr[l]   * grad;
        flag[l] = 1;
      }
    }

#ifdef FTRACE
    ftrace_region_end  ("calc_forces");
    ftrace_region_begin("store_forces");
#endif

    for (i=imin; i<imax; i++) {
#ifdef SX
#pragma cdir vector,nodep,loopcnt=256
#endif
      for (k=lbeg[i-imin]; k<lend[i-imin]; k++) {
        if (flag[k]) {  /* is there a better way to deal with this? */
          int j = lj[k];
          atoms.kraft X(j) -= ff[k].x;
          atoms.kraft Y(j) -= ff[k].y;
          atoms.kraft Z(j) -= ff[k].z;
          atoms.pot_eng[j] += ee[k];
          fi[i].x          += ff[k].x;
          fi[i].y          += ff[k].y;
          fi[i].z          += ff[k].z;
          ei[i]            += ee[k];
        }
      }
    }

#ifdef FTRACE
    ftrace_region_end("store_forces");
#endif

  }  /* loop over blocks */
  if (is_short) printf("short distance!\n");

  /* store the remaining forces */
  for (i=0; i<atoms.n; i++) {
    atoms.kraft X(i) += fi[i].x;
    atoms.kraft Y(i) += fi[i].y;
    atoms.kraft Z(i) += fi[i].z;
    atoms.pot_eng[i] += ei[i];
    tot_pot_energy   += ei[i];
  }

#ifdef MPI
  /* sum up results of different CPUs */
  tmpvec1[0] = tot_pot_energy;
  tmpvec1[1] = virial;
  MPI_Allreduce( tmpvec1, tmpvec2, 2, REAL, MPI_SUM, cpugrid); 
  tot_pot_energy = tmpvec2[0];
  virial         = tmpvec2[1];
#endif

  /* add forces back to original cells/cpus */
  send_forces(add_forces,pack_forces,unpack_forces);

}

#endif /* VEC3 */

/******************************************************************************
*
*  check_nblist
*
******************************************************************************/

void check_nblist()
{
  int  i;
  real dx, dy, dz, r2, max1=0.0, max2;

  for (i=0; i<atoms.n; i++) {
    dx = atoms.ort X(i) - atoms.nbl_pos X(i);
    dy = atoms.ort Y(i) - atoms.nbl_pos Y(i);
    dz = atoms.ort Z(i) - atoms.nbl_pos Z(i);
    r2 = dx*dx + dy*dy +dz*dz;
    if (r2 > max1) max1 = r2;
  }

#ifdef MPI
  MPI_Allreduce( &max1, &max2, 1, REAL, MPI_MAX, cpugrid); 
#else
  max2 = max1;
#endif
  if (max2 > SQR(0.5*nbl_margin)) have_valid_nbl = 0;
}
