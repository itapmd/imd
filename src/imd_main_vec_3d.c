
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

/* we can use several versions of calc_forces:

   VERSION1:  straight forward first attempt
   VERSION2:  vectorizing LLC algorithm as in Grest, Dünwald, Kremer
   VERSION3:  refinement of VERSION1
   VERSION4:  intermediate between VERSION1 and VERSION3

   on SX6, only VERSION2 currently works; strange: VERSION1 used to work!
   on IA64, VERSION3 is surprisingly fast!
*/

#define VERSION3

/******************************************************************************
*
*  calc_forces 
*
*  We use buffer cells on the surface of the cell array, even for 
*  the serial version. Atoms in the buffer cells have positions with
*  periodic boundary conditions applied, where applicable.
*
******************************************************************************/

#ifdef VERSION1

#define MC 40      /* maximum atoms per minicell */
#define N1 MC*14
#define N2 MC*4
#define N0 17000   /* at least atoms.n */

void calc_forces(int steps)
{
  long   i, j, k, l, is_short=0, max, len;
  int    tab[N0*N1], tlen[N0], itypes[N0*N2];
  int    li[N0*N2], lj[N0*N2], lbeg[N0], lend[N0];
  vektor dd[N0*N2], ff[N0*N2], fi[N0];
  real   rr[N0*N2], ee[N0*N2], ei[N0];
  real   tmpvec1[2], tmpvec2[2] = {0.0, 0.0}, tmp;

#ifdef MPI
  if ((steps == steps_min) || (0 == steps % BUFSTEP)) setup_buffers();
#endif

  /* fill the buffer cells */
  send_cells(copy_cell,pack_cell,unpack_cell);

  /* clear global accumulation variables */
  tot_pot_energy = 0.0;
  virial = 0.0;

  /* clear per atom accumulation variables */
  for (k=0; k<atoms.n_buf; k++) {
    atoms.kraft X(k) = 0.0;
    atoms.kraft Y(k) = 0.0;
    atoms.kraft Z(k) = 0.0;
    atoms.pot_eng[k] = 0.0;
  }

#ifdef FTRACE
  ftrace_region_begin("get_indices");
#endif

  /* course neighbor table */
  for (k=0; k<atoms.n; k++) tlen[k] = 0;
  for (k=0; k<npairs[0]; k++) {

    pair *P;
    minicell *p, *q;
    int i, j, ii, jj, jstart;

    P = pairs[0] + k;
    p = cell_array + P->np;
    q = cell_array + P->nq;

    /* for each atom in first cell */
    for (i=0; i<p->n; i++) {
      int *tt;
      ii = p->ind[i];
      jstart = ((p==q) ? i+1 : 0);
      tt = tab + N1*ii + tlen[ii] - jstart;

      /* for each atom in neighbouring cell */
      for (j = jstart; j < q->n; j++) {
        tt[j] = q->ind[j];
      }
      tlen[ii] += q->n - jstart;
    }
  }

#ifdef FTRACE
  ftrace_region_end  ("get_indices");
  ftrace_region_begin("calc_distances");
#endif

  /* narrow neighbor table */
  len=0;
  lbeg[0]=0;
  for (i=0; i<atoms.n; i++) {
    int j, ityp;
    real *d1, *d2, dx, dy, dz, r2;
    d1 = atoms.ort + DIM * i;
    ityp = VSORTE(&atoms,i);
    for (k=0; k<tlen[i]; k++) {
      j  = tab[N1*i+k];
      d2 = atoms.ort + DIM * j;
      dx = d2[0]-d1[0];
      dy = d2[1]-d1[1];
      dz = d2[2]-d1[2];
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < cellsz) {
        li[len]   = i;
        lj[len]   = j;
        rr[len]   = r2;
        dd[len].x = dx;
        dd[len].y = dy;
        dd[len].z = dz;
        itypes[len] = ityp;
        len++;
      }
    }
    lend[i  ] = len;
    lbeg[i+1] = len;
  }

#ifdef FTRACE
  ftrace_region_end  ("calc_distances");
  ftrace_region_begin("calc_forces");
#endif

  /* compute forces */
  for (l=0; l<len; l++) {

    real grad;
    int  col, inc = ntypes * ntypes; 

    col = itypes[l] * ntypes + VSORTE(&atoms,lj[l]);
    /* col = VSORTE(&atoms,li[l]) * ntypes + VSORTE(&atoms,lj[l]); */

    /* compute pair interactions */
    if (rr[l] <= pair_pot.end[col]) {

      /* beware: we must not use k as index in the macro's arguments! */
#ifdef MULTIPOT
      int pp = l % N_POT_TAB;
      PAIR_INT(ee[l], grad, pair_pot_ar[pp], col, inc, rr[l], is_short)
#else
#ifdef LINPOT
      PAIR_INT_LIN(ee[l], grad, pair_pot_lin, col, inc, rr[l], is_short)
#else
      PAIR_INT(ee[l], grad, pair_pot, col, inc, rr[l], is_short)
#endif
#endif

      /* store force in temporary variables */
      ff[l].x = dd[l].x * grad;
      ff[l].y = dd[l].y * grad;
      ff[l].z = dd[l].z * grad;
      virial -= rr[l]   * grad;
    }
  }
  if (is_short) printf("short distance!\n");

#ifdef FTRACE
  ftrace_region_end  ("calc_forces");
  ftrace_region_begin("store_forces");
#endif

  /* accumulate forces */
  for (i=0; i<atoms.n; i++) {
    fi[i].x = 0.0;
    fi[i].y = 0.0;
    fi[i].z = 0.0;
    ei[i]   = 0.0;
  }
  for (i=0; i<atoms.n; i++) {
#ifdef SX
#pragma vdir vector,nodep
#endif
    for (k=lbeg[i]; k<lend[i]; k++) {
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
  for (i=0; i<atoms.n; i++) {
    atoms.kraft X(i) += fi[i].x;
    atoms.kraft Y(i) += fi[i].y;
    atoms.kraft Z(i) += fi[i].z;
    atoms.pot_eng[i] += ei[i];
    tot_pot_energy   += ei[i];
  }

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

#endif /* VERSION1 */

/******************************************************************************
*
*  calc_forces 
*
*  We use buffer cells on the surface of the cell array, even for 
*  the serial version. Atoms in the buffer cells have positions with
*  periodic boundary conditions applied, where applicable.
*
*  This version implements the vectorizing LLC algorithm from
*  Grest, Dünweg & Kremer, Comp. Phys. Comm. 55, 269-285 (1989).
*
******************************************************************************/

#ifdef VERSION2

#define N0 17000   /* at least atoms.n */

void calc_forces(int steps)
{
  long   i, j, k, l, m, n, is_short=0, len1, len2, max_cell;

  int    li1[N0], lj1[N0], li2[N0], lj2[N0], ll[N0];
  vektor dd[N0], ff[N0];
  real   rr[N0], ee[N0];

  real   tmpvec1[2], tmpvec2[2] = {0.0, 0.0}, tmp;

#ifdef MPI
  if ((steps == steps_min) || (0 == steps % BUFSTEP)) setup_buffers();
#endif

  /* fill the buffer cells */
  send_cells(copy_cell,pack_cell,unpack_cell);

  /* get maximum cell occupancy */
  max_cell = 0;
  for (i=0; i<nallcells; i++) {
    n = (cell_array+i)->n;
    if (max_cell < n) max_cell = n;
  }

  /* clear global accumulation variables */
  tot_pot_energy = 0.0;
  virial = 0.0;

  /* clear per atom accumulation variables */
  for (k=0; k<atoms.n_buf; k++) {
    atoms.kraft X(k) = 0.0;
    atoms.kraft Y(k) = 0.0;
    atoms.kraft Z(k) = 0.0;
    atoms.pot_eng[k] = 0.0;
  }

  for (n=0; n<14; n++) {

    for (k=0; k<max_cell; k++) {

#ifdef FTRACE
      ftrace_region_begin("get_indices");
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

#ifdef FTRACE
      ftrace_region_end  ("get_indices");
      ftrace_region_begin("calc_distances");
#endif

      /* reduce pair list */
      len2=0;
      for (l=0; l<len1; l++) {
        real *d1, *d2, dx, dy, dz, r2;
        d1 = atoms.ort + DIM * li1[l];
        d2 = atoms.ort + DIM * lj1[l];
        dx = d2[0]-d1[0];
        dy = d2[1]-d1[1];
        dz = d2[2]-d1[2];
        r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < cellsz) {
          li2[len2]   = li1[l];
          lj2[len2]   = lj1[l];
          rr [len2]   = r2;
          dd [len2].x = dx;
          dd [len2].y = dy;
          dd [len2].z = dz;
          len2++;
        }
      }

#ifdef FTRACE
      ftrace_region_end  ("calc_distances");
      ftrace_region_begin("calc_forces");
#endif

#ifdef SX
#pragma vdir vector,nodep
#endif
      /* compute forces */
      for (l=0; l<len2; l++) {

        real grad;
        int  j, col, inc = ntypes * ntypes; 

        col = VSORTE(&atoms,li2[l]) * ntypes + VSORTE(&atoms,lj2[l]);

        /* compute pair interactions */
        if (rr[l] <= pair_pot.end[col]) {

          /* beware: we must not use k as index in the macro's arguments! */
          PAIR_INT(ee[l], grad, pair_pot, col, inc, rr[l], is_short)

          /* store force in temporary variables */
          ff[l].x         = dd[l].x * grad;
          ff[l].y         = dd[l].y * grad;
          ff[l].z         = dd[l].z * grad;
          virial         -= rr[l]   * grad;
          tot_pot_energy += ee[l];

          /* store force to second particle */
          j = lj2[l];
          atoms.kraft X(j) -= ff[l].x;
          atoms.kraft Y(j) -= ff[l].y;
          atoms.kraft Z(j) -= ff[l].z;
          atoms.pot_eng[j] += ee[l];
        }
      }
#ifdef FTRACE
      ftrace_region_end  ("calc_forces");
      ftrace_region_begin("store_forces");
#endif

#ifdef SX
#pragma vdir vector,nodep
#endif
      /* accumulate remaining forces */
      for (l=0; l<len2; l++) {
        int i = li2[l];
        atoms.kraft X(i) += ff[l].x;
        atoms.kraft Y(i) += ff[l].y;
        atoms.kraft Z(i) += ff[l].z;
        atoms.pot_eng[i] += ee[l];
      }

#ifdef FTRACE
      ftrace_region_end("store_forces");
#endif

    }
  }
  if (is_short) printf("short distance!\n");

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

#endif /* VERSION2 */


/******************************************************************************
*
*  calc_forces 
*
*  We use buffer cells on the surface of the cell array, even for 
*  the serial version. Atoms in the buffer cells have positions with
*  periodic boundary conditions applied, where applicable.
*
******************************************************************************/

#ifdef VERSION3

/*
#define DEBUG
*/

#define MC 40      /* maximum atoms per minicell */
#define N1 MC*14
#define N2 MC*4
#define N0 17000   /* at least atoms.n */

void calc_forces(int steps)
{
  long   c, i, j, k, l, m, is_short=0, max, len;
  int    tab[MC*N1], tlen[N0], itypes[MC*N2];
  int    li[MC*N2], lj[MC*N2], lbeg[MC], lend[MC];
  vektor dd[MC*N2], ff[MC*N2], fi[N0];
  real   rr[MC*N2], ee[MC*N2], ei[N0];
  real   tmpvec1[2], tmpvec2[2] = {0.0, 0.0}, tmp;

#ifdef MPI
  if ((steps == steps_min) || (0 == steps % BUFSTEP)) setup_buffers();
#endif

  /* fill the buffer cells */
  send_cells(copy_cell,pack_cell,unpack_cell);

  /* clear global accumulation variables */
  tot_pot_energy = 0.0;
  virial = 0.0;

  /* clear per atom accumulation variables */
  for (k=0; k<atoms.n_buf; k++) {
    atoms.kraft X(k) = 0.0;
    atoms.kraft Y(k) = 0.0;
    atoms.kraft Z(k) = 0.0;
    atoms.pot_eng[k] = 0.0;
    /* temporary storage */
    fi[k].x          = 0.0;
    fi[k].y          = 0.0;
    fi[k].z          = 0.0;
    ei[k]            = 0.0;
  }

  /* for all cells */
  for (c=0; c<ncells; c++) {

    minicell *p;

    p = cell_array + cnbrs[c].np;

#ifdef FTRACE
    ftrace_region_begin("get_indices");
#endif

    /* for each atom in cell */
    for (i=0; i<p->n; i++) {

      int *tt, ii;
      minicell *q;

      ii = p->ind[i];

      /* indices of particles from the same cell */
      tt       = tab + N1*i - (i+1);
      tlen[i]  =            - (i+1);
      for (j = i+1; j < p->n; j++) tt[j] = p->ind[j];
      tt      += p->n;
      tlen[i] += p->n;

      /* indices of particles from neighbor cells */
      for (m=1; m<14; m++) {
        q = cell_array + cnbrs[c].nq[m];
#ifdef SX
#pragma vdir vector,nodep
#endif
        for (j = 0; j < q->n; j++) tt[j] = q->ind[j];
        tt      += q->n;
        tlen[i] += q->n;
      }
    }
#ifdef DEBUG
    for (i=0; i<tlen[0]; i++) printf("%d ",tab[i]);
    printf("\n%d\n",tlen[0]);
#endif

#ifdef FTRACE
    ftrace_region_end  ("get_indices");
    ftrace_region_begin("calc_distances");
#endif

    /* narrow neighbor table */
    len=0;
    lbeg[0]=0;
    for (i=0; i<p->n; i++) {
      int jj, ii, ityp;
      real *d1, *d2, dx, dy, dz, r2;
      ii = p->ind[i];
      d1 = atoms.ort + DIM * ii;
      ityp = VSORTE(&atoms,ii);
      for (k=0; k<tlen[i]; k++) {
        jj = tab[N1*i+k];
        d2 = atoms.ort + DIM * jj;
        dx = d2[0]-d1[0];
        dy = d2[1]-d1[1];
        dz = d2[2]-d1[2];
        r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < cellsz) {
          li[len]   = ii;
          lj[len]   = jj;
          rr[len]   = r2;
          dd[len].x = dx;
          dd[len].y = dy;
          dd[len].z = dz;
          itypes[len] = ityp;
          len++;
        }
      }
      lend[i  ] = len;
      lbeg[i+1] = len;
    }
#ifdef DEBUG
    for (i=0; i<lend[0]; i++) printf("%d %d %f\n",li[i],lj[i],rr[i]);
    printf("%d\n",lend[0]);
#endif

#ifdef FTRACE
    ftrace_region_end  ("calc_distances");
    ftrace_region_begin("calc_forces");
#endif

    /* compute forces */
    for (l=0; l<len; l++) {

      real grad;
      int  col, inc = ntypes * ntypes; 

      col = itypes[l] * ntypes + VSORTE(&atoms,lj[l]);
      /* col = VSORTE(&atoms,li[l]) * ntypes + VSORTE(&atoms,lj[l]); */

      /* compute pair interactions */
      if (rr[l] <= pair_pot.end[col]) {

        /* beware: we must not use k as index in the macro's arguments! */
#ifdef MULTIPOT
        int pp = l % N_POT_TAB;
        PAIR_INT(ee[l], grad, pair_pot_ar[pp], col, inc, rr[l], is_short)
#else
#ifdef LINPOT
        PAIR_INT_LIN(ee[l], grad, pair_pot_lin, col, inc, rr[l], is_short)
#else
        PAIR_INT(ee[l], grad, pair_pot, col, inc, rr[l], is_short)
#endif
#endif

        /* store force in temporary variables */
        ff[l].x = dd[l].x * grad;
        ff[l].y = dd[l].y * grad;
        ff[l].z = dd[l].z * grad;
        virial -= rr[l]   * grad;
      }
    }
    if (is_short) printf("short distance!\n");

#ifdef DEBUG
    for (i=0; i<lend[0]; i++) 
      printf("%f %f %f %f\n",ff[i].x,ff[i].y,ff[i].z,ee[i]);
    printf("%d\n",lend[0]);
#endif


#ifdef FTRACE
    ftrace_region_end  ("calc_forces");
    ftrace_region_begin("store_forces");
#endif

    for (l=0; l<p->n; l++) {
#ifdef SX
#pragma vdir vector,nodep
#endif
      for (k=lbeg[l]; k<lend[l]; k++) {
        int i=li[k], j = lj[k];
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

#ifdef FTRACE
    ftrace_region_end("store_forces");
#endif

  }  /* loop over all cells */

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

#endif /* VERSION3 */


#ifdef VERSION4

#define MC 40      /* maximum atoms per minicell */
#define N1 MC*14
#define N2 MC*4
#define N0 17000   /* at least atoms.n */

void calc_forces(int steps)
{
  long   i, j, k, l, m, is_short=0, max, len;
  int    tab[N0*N1], tlen[N0], itypes[N0*N2];
  int    li[N0*N2], lj[N0*N2], lbeg[N0], lend[N0];
  vektor dd[N0*N2], ff[N0*N2], fi[N0];
  real   rr[N0*N2], ee[N0*N2], ei[N0];
  real   tmpvec1[2], tmpvec2[2] = {0.0, 0.0}, tmp;

#ifdef MPI
  if ((steps == steps_min) || (0 == steps % BUFSTEP)) setup_buffers();
#endif

  /* fill the buffer cells */
  send_cells(copy_cell,pack_cell,unpack_cell);

  /* clear global accumulation variables */
  tot_pot_energy = 0.0;
  virial = 0.0;

  /* clear per atom accumulation variables */
  for (k=0; k<atoms.n_buf; k++) {
    atoms.kraft X(k) = 0.0;
    atoms.kraft Y(k) = 0.0;
    atoms.kraft Z(k) = 0.0;
    atoms.pot_eng[k] = 0.0;
  }

#ifdef FTRACE
  ftrace_region_begin("get_indices");
#endif

  /* course neighbor table */
  for (k=0; k<ncells; k++) {

    minicell *p, *q;
    int *tt, ii;

    p = cell_array + cnbrs[k].np;

    /* for each atom in first cell */
    for (i=0; i<p->n; i++) {

      ii = p->ind[i];

      /* indices of particles from the same cell */
      tt       = tab + N1*ii - (i+1);
      tlen[ii] =             - (i+1);
      for (j = i+1; j < p->n; j++) tt[j] = p->ind[j];
      tt       += p->n;
      tlen[ii] += p->n;

      /* indices of particles from neighbor cells */
      for (m=1; m<14; m++) {
        q = cell_array + cnbrs[k].nq[m];
        for (j = 0; j < q->n; j++) tt[j] = q->ind[j];
        tt       += q->n;
        tlen[ii] += q->n;
      }
    }
  }

#ifdef FTRACE
  ftrace_region_end  ("get_indices");
  ftrace_region_begin("calc_distances");
#endif

  /* narrow neighbor table */
  len=0;
  lbeg[0]=0;
  for (i=0; i<atoms.n; i++) {
    int j, ityp;
    real *d1, *d2, dx, dy, dz, r2;
    d1 = atoms.ort + DIM * i;
    ityp = VSORTE(&atoms,i);
    for (k=0; k<tlen[i]; k++) {
      j  = tab[N1*i+k];
      d2 = atoms.ort + DIM * j;
      dx = d2[0]-d1[0];
      dy = d2[1]-d1[1];
      dz = d2[2]-d1[2];
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < cellsz) {
        li[len]   = i;
        lj[len]   = j;
        rr[len]   = r2;
        dd[len].x = dx;
        dd[len].y = dy;
        dd[len].z = dz;
        itypes[len] = ityp;
        len++;
      }
    }
    lend[i  ] = len;
    lbeg[i+1] = len;
  }

#ifdef FTRACE
  ftrace_region_end  ("calc_distances");
  ftrace_region_begin("calc_forces");
#endif

  /* compute forces */
  for (l=0; l<len; l++) {

    real grad;
    int  col, inc = ntypes * ntypes; 

    col = itypes[l] * ntypes + VSORTE(&atoms,lj[l]);
    /* col = VSORTE(&atoms,li[l]) * ntypes + VSORTE(&atoms,lj[l]); */

    /* compute pair interactions */
    if (rr[l] <= pair_pot.end[col]) {

      /* beware: we must not use k as index in the macro's arguments! */
#ifdef MULTIPOT
      int pp = l % N_POT_TAB;
      PAIR_INT(ee[l], grad, pair_pot_ar[pp], col, inc, rr[l], is_short)
#else
#ifdef LINPOT
      PAIR_INT_LIN(ee[l], grad, pair_pot_lin, col, inc, rr[l], is_short)
#else
      PAIR_INT(ee[l], grad, pair_pot, col, inc, rr[l], is_short)
#endif
#endif

      /* store force in temporary variables */
      ff[l].x = dd[l].x * grad;
      ff[l].y = dd[l].y * grad;
      ff[l].z = dd[l].z * grad;
      virial -= rr[l]   * grad;
    }
  }
  if (is_short) printf("short distance!\n");

#ifdef FTRACE
  ftrace_region_end  ("calc_forces");
  ftrace_region_begin("store_forces");
#endif

  /* accumulate forces */
  for (i=0; i<atoms.n; i++) {
    fi[i].x = 0.0;
    fi[i].y = 0.0;
    fi[i].z = 0.0;
    ei[i]   = 0.0;
  }
  for (i=0; i<atoms.n; i++) {
#ifdef SX
#pragma vdir vector,nodep
#endif
    for (k=lbeg[i]; k<lend[i]; k++) {
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
  for (i=0; i<atoms.n; i++) {
    atoms.kraft X(i) += fi[i].x;
    atoms.kraft Y(i) += fi[i].y;
    atoms.kraft Z(i) += fi[i].z;
    atoms.pot_eng[i] += ei[i];
    tot_pot_energy   += ei[i];
  }

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

#endif /* VERSION4 */
