
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2004 Institute for Theoretical and Applied Physics,
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
#define COL(i) 0
#else
#define COL(i) col[i]
#endif

/* int li[N0*N2], lj[N0*N2], col[N0*N2], lbeg[N1], lend[N1], ltot; */ 
int *li=NULL, *lj=NULL, *col=NULL, *lbeg=NULL, *lend=NULL, ltot; 

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
  static int max_atoms=0, max_pairs=0, max_cell_max=0;
  long   i, j, k, l, m, n, len1, len2, max_cell;
  int    *li1 = NULL, *lj1 = NULL;
  /* int    li1[N0], lj1[N0]; */

#ifdef MPI
  if ((steps == steps_min) || (0 == steps % BUFSTEP)) setup_buffers();
#endif

  /* update cell decomposition */
  do_boundaries();
  fix_cells();

  /* update reference positions */
  for (i=0; i<DIM*atoms.n; i++) atoms.nbl_pos[i] = atoms.ort[i];

  /* fill the buffer cells */
  send_cells(copy_cell,pack_cell,unpack_cell);

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
    max_atoms = (int) (atoms.n * 1.2);
    li1 = (int *) malloc( max_atoms * sizeof(int) );
    lj1 = (int *) malloc( max_atoms * sizeof(int) );
    if ((li1==NULL) || (lj1==NULL))
      error("cannot allocate temporary pair arrays");
  }
  if (max_pairs < 4 * max_cell * atoms.n) {
    free(li);
    free(lj);
    free(col);
    max_pairs = (int) (4 * max_cell * atoms.n * 1.2);
    li  = (int *) malloc( max_pairs * sizeof(int) );
    lj  = (int *) malloc( max_pairs * sizeof(int) );
    col = (int *) malloc( max_pairs * sizeof(int) );
    if ((li==NULL) || (lj==NULL) || (col==NULL))
      error("cannot allocate temporary pair arrays");
  }
  if (max_cell_max < max_cell) {
    free(lbeg);
    free(lend);
    max_cell_max = (int) (max_cell * 1.2);
    lbeg = (int *) malloc( 14 * max_cell_max * sizeof(int) );
    lend = (int *) malloc( 14 * max_cell_max * sizeof(int) );
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
          int cjn = cnbrs[m].nq[n];
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
      if (lend[ltot] > lbeg[ltot]) ltot++;
    }
  }
  if (max_pairs < len2) error("neighbor table overflow"); 

#ifndef MONO

#ifdef FTRACE
  ftrace_region_begin("nbl_make_col");
#endif

  /* precompute potential column for each atom pair */
  for (l=0; l<lend[ltot-1]; l++) {
    col[l] = SORTE(&atoms,li[l]) * ntypes + SORTE(&atoms,lj[l]);
  }

#ifdef FTRACE
  ftrace_region_end("nbl_make_col");
#endif

#endif /* not MONO */

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
  int    k, m, is_short=0;
  real   rr, ee, tmpvec1[2], tmpvec2[2] = {0.0, 0.0};
  vektor dd, ff;
  static real *kraft1 = NULL, *pot_eng1 = NULL;

#ifdef MPI
  if ((steps == steps_min) || (0 == steps % BUFSTEP)) setup_buffers();
#endif

  /* fill the buffer cells */
  send_cells(copy_cell,pack_cell,unpack_cell);

  /* clear global accumulation variables */
  tot_pot_energy = 0.0;
  virial = 0.0;
  nfc++;

  /* (re)allocate temporary force and energy arrays */
  if (len < atoms.n_buf) {
    free(pot_eng1);
    free(kraft1);
    len = (int) (atoms.n_buf * 1.2);
    pot_eng1 = (real *) malloc(       len * sizeof(real) );
    kraft1   = (real *) malloc( DIM * len * sizeof(real) );
    if ((pot_eng1 == NULL) || (kraft1 == NULL))
      error("cannot allocate temporary force and energy arrays");
  }

  /* clear per atom accumulation variables */
  for (k=0; k<DIM*atoms.n_buf; k++) atoms.kraft  [k] = 0.0;
  for (k=0; k<    atoms.n_buf; k++) atoms.pot_eng[k] = 0.0;
  for (k=0; k<DIM*atoms.n_buf; k++)      kraft1  [k] = 0.0;
  for (k=0; k<    atoms.n_buf; k++)      pot_eng1[k] = 0.0;

  /* loop over independent pair sublists */
  for (k=0; k<ltot; k++) {

#ifdef FTRACE
    ftrace_region_begin("calc_forces");
#endif

#ifdef SX
#pragma vdir vector,nodep
#endif
    /* compute forces */
    for (m=0; m<lend[k]-lbeg[k]; m++) {

      real grad, *d1, *d2;
      int  i, j, l = lbeg[k] + m; 
      const int inc = ntypes * ntypes;
      d1   = atoms.ort + DIM * li[l];
      d2   = atoms.ort + DIM * lj[l];
      dd.x = d2[0]-d1[0];
      dd.y = d2[1]-d1[1];
      dd.z = d2[2]-d1[2];
      rr   = SPROD(dd,dd);

      /* compute pair interactions */
      if (rr <= pair_pot.end[COL(l)]) {

        /* beware: we must not use k as index in the macro's arguments! */
#ifdef MULTIPOT
        int pp = m % N_POT_TAB;
        PAIR_INT(ee, grad, pair_pot_ar[pp], COL(l), inc, rr, is_short)
#else
#ifdef LINPOT
        PAIR_INT_LIN(ee, grad, pair_pot_lin, COL(l), inc, rr, is_short)
#else
#ifdef LJ
        PAIR_INT_LJ_VEC(ee, grad, COL(l), rr)
#else
        PAIR_INT(ee, grad, pair_pot, COL(l), inc, rr, is_short)
#endif
#endif
#endif
        /* store force in temporary variables */
        ff.x            = dd.x * grad;
        ff.y            = dd.y * grad;
        ff.z            = dd.z * grad;
        virial         -= rr   * grad;
        tot_pot_energy += ee;

        /* store force to first particle */
        i = li[l];
        kraft1 X(i) += ff.x;
        kraft1 Y(i) += ff.y;
        kraft1 Z(i) += ff.z;
        pot_eng1[i] += ee;

        /* store force to second particle */
        j = lj[l];
        atoms.kraft X(j) -= ff.x;
        atoms.kraft Y(j) -= ff.y;
        atoms.kraft Z(j) -= ff.z;
        atoms.pot_eng[j] += ee;
      }
    }

#ifdef FTRACE
    ftrace_region_end("calc_forces");
#endif

  }
  if (is_short) printf("short distance!\n");

#ifdef FTRACE
  ftrace_region_begin("store_forces");
#endif

  /* add up forces on first and second particles */
  for (k=0; k<DIM*atoms.n_buf; k++) atoms.kraft  [k] += kraft1  [k];
  for (k=0; k<    atoms.n_buf; k++) atoms.pot_eng[k] += pot_eng1[k];

#ifdef FTRACE
    ftrace_region_end("store_forces");
#endif

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

#ifdef MPI
  if ((steps == steps_min) || (0 == steps % BUFSTEP)) setup_buffers();
#endif

  /* update cell decomposition */
  do_boundaries();
  fix_cells();

  /* update reference positions */
  for (i=0; i<DIM*atoms.n; i++) atoms.nbl_pos[i] = atoms.ort[i];

  /* fill the buffer cells */
  send_cells(copy_cell,pack_cell,unpack_cell);

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
#pragma vdir vector,nodep
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

#ifdef MPI
  if ((steps == steps_min) || (0 == steps % BUFSTEP)) setup_buffers();
#endif

  /* fill the buffer cells */
  send_cells(copy_cell,pack_cell,unpack_cell);

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
#pragma vdir vector,loopcnt=256
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
#pragma vdir vector,nodep,loopcnt=256
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
  if (max2 > SQR(0.5*nbl_margin)) make_nblist();
}
