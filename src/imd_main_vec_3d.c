
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

/******************************************************************************
*
*  calc_forces 
*
*  We use buffer cells on the surface of the cell array, even for 
*  the serial version. Atoms in the buffer cells have positions with
*  periodic boundary conditions applied, where applicable.
*
******************************************************************************/

#define MC 40      /* maximum atoms per minicell */
#define N1 MC*14
#define N2 MC*4
#define N0 17000   /* at least atoms.n */

void calc_forces(int steps)
{
  int    i, j, k, l, is_short=0, max, len;
  int    tab[N0*N1], tlen[N0];
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

  /* narrow neighbor table */
  len=0;
  lbeg[0]=0;
  for (i=0; i<atoms.n; i++) {
    int j, col, inc=ntypes*ntypes;
    real *d1, *d2, dx, dy, dz, r2;
    d1 = atoms.ort + DIM * i;
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
        len++;
      }
    }
    lend[i  ] = len;
    lbeg[i+1] = len;
  }

  /* compute forces */
  for (l=0; l<len; l++) {

    real grad;
    int  col, inc = ntypes * ntypes; 

    col = SORTE(&atoms,li[l]) * ntypes + SORTE(&atoms,lj[l]);

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
