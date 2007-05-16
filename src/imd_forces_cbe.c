
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2007 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_forces_cbe.c -- force loop for CBE
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

typedef float flt;

typedef struct {
  int   k, n1, n1_max, n2, n2_max, len, len_max;
  flt   totpot, virial;
  flt   *pos, *force;
  int   *typ, *tb, *ta, *te;
} wp_t;

typedef struct {
  int ntypes;
  flt *r2cut, *lj_sig, *lj_eps, *lj_shift;
} pt_t;

int  *ta=NULL, *te=NULL, *tb=NULL, *tb_off=NULL, *at_off=NULL, nb_max=0;
wp_t *wp=NULL;
pt_t pt;

/******************************************************************************
*
*  collect potential data for SPE
*
******************************************************************************/

void mk_pt(void)
{
  int i, j, k, col;
  flt r2, r6;

  pt.ntypes   = ntypes;
  pt.r2cut    = (flt *) malloc( 4 * ntypes * ntypes * sizeof(flt) );
  pt.lj_sig   = (flt *) malloc( 4 * ntypes * ntypes * sizeof(flt) );
  pt.lj_eps   = (flt *) malloc( 4 * ntypes * ntypes * sizeof(flt) );
  pt.lj_shift = (flt *) malloc( 4 * ntypes * ntypes * sizeof(flt) );
  if ((NULL==pt.r2cut) || (NULL==pt.lj_sig) || (NULL==pt.lj_eps) || 
      (NULL==pt.lj_shift)) error("cannot allocate potential package");
  for (i=0; i<ntypes; i++)
    for (j=0; j<ntypes; j++) {
      col = i*ntypes + j;
      for (k=0; k<4; k++) { /* make vectors: 4 copies of each value */
        pt.r2cut  [4*col+k] = (flt) SQR(r_cut_lin   [col]);
        pt.lj_sig [4*col+k] = (flt) SQR(lj_sigma_lin[col]);
        pt.lj_eps [4*col+k] = (flt) lj_epsilon_lin  [col] ;
        r2 = pt.lj_sig[4*col+k] / pt.r2cut[4*col+k];
        r6 = r2 * r2 * r2;
        pt.lj_shift[4*col+k] = -pt.lj_eps[4*col+k] * r6 * (r6 - 2.0);
      }
    }
}

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
      d1.z = ORT(p,i,Z);

      /* for each neighboring atom */
      for (m=0; m<14; m++) {   /* this is not TWOD ready! */

        int    c2, jstart, j;
        real   r2;
        cell   *q;

        c2 = cnbrs[c].nq[m];
        if (c2<0) continue;
        jstart = (c2==c1) ? i+1: 0;
        q = cell_array + c2;
        for (j=jstart; j<q->n; j++) {
          vektor d;
          d.x = ORT(q,j,X) - d1.x;
          d.y = ORT(q,j,Y) - d1.y;
          d.z = ORT(q,j,Z) - d1.z;
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
  static int at_max=0, pa_max=0, ncell_max=0, tn_max=0, cl_max=0;
  static int *nb=NULL, *ty=NULL;
  int  c, i, k, n, tn, at, cc, atm, l, t, j, nn;

  /* update reference positions */
  for (k=0; k<ncells; k++) {
    cell *p = cell_array + cnbrs[k].np;
    for (i=0; i<p->n; i++) {
      NBL_POS(p,i,X) = ORT(p,i,X);
      NBL_POS(p,i,Y) = ORT(p,i,Y);
      NBL_POS(p,i,Z) = ORT(p,i,Z);
    }
  }

  /* (re-)allocate offsets */
  if (ncells > cl_max) {
    free(at_off);
    free(tb_off);
    cl_max = ncells;
    at_off = (int *) malloc( (cl_max+1) * sizeof(int) );
    tb_off = (int *) malloc( (cl_max+1) * sizeof(int) );
    if ((NULL==at_off) || (NULL==tb_off))
      error("cannot allocate neighbor table");
  }

  /* count atom numbers */
  at=0; atm=0;
  for (k=0; k<ncells; k++) {
    cell *p = CELLPTR(k);
    at_off[k] = at;
    at += p->n * ntypes;
    atm = MAX( atm, p->n );
  }

  /* allocate intermediate neighbor array */
  if (15*atm > tn_max) {
    free(nb);
    free(ty);
    tn_max = 16*atm;
    nb = (int *) malloc( tn_max * sizeof(int) );
    ty = (int *) malloc( tn_max * sizeof(int) );
    if ((NULL==nb) || (NULL==ty)) error("cannot allocate neighbor table");
  }

  /* (re-)allocate neighbor table */
  if (at >= at_max) {
    free(ta);
    free(te);
    at_max = (int) (1.1 * at);
    ta = (int *) malloc(at_max * sizeof(int));
    te = (int *) malloc(at_max * sizeof(int));
  }
  if (NULL==tb) {
    if (0==last_nbl_len) nb_max = (int) (nbl_size * estimate_nblist_size());
    else                 nb_max = (int) (nbl_size * last_nbl_len);
    tb = (int *) malloc(nb_max * sizeof(int));
  }
  else if (last_nbl_len * sqrt(nbl_size) > nb_max) {
    free(tb);
    nb_max = (int  ) (nbl_size * last_nbl_len);
    tb     = (int *) malloc(nb_max * sizeof(int));
  }
  if ((ta==NULL) || (te==NULL) || (tb==NULL)) 
    error("cannot allocate neighbor table");

  /* for all cells */
  n=0; l=0;
  for (c=0; c<ncells; c++) {

    int  c1 = cnbrs[c].np;
    cell *p = cell_array + c1;

    at_off[c] = l;
    tb_off[c] = n;

    /* for each atom in cell */
    for (i=0; i<p->n; i++) {

      int    m, off;
      vektor d1;

      d1.x = ORT(p,i,X);
      d1.y = ORT(p,i,Y);
      d1.z = ORT(p,i,Z);

      /* for each neighboring atom */
      off = 0; tn = 0;
      for (m=0; m<14; m++) {   /* this is not TWOD ready! */
        int  c2, jstart, j;
        cell *q;
        c2 = cnbrs[c].nq[m];
        if (c2<0) continue;
        jstart = (c2==c1) ? i+1 : 0;
        q = cell_array + c2;
        for (j=jstart; j<q->n; j++) {
          vektor d;
          real   r2;
          d.x = ORT(q,j,X) - d1.x;
          d.y = ORT(q,j,Y) - d1.y;
          d.z = ORT(q,j,Z) - d1.z;
          r2  = SPROD(d,d);
          if (r2 < cellsz) {
            nb[tn  ] = off + j;
            ty[tn++] = SORTE(q,j);
            if (tn>tn_max) error("tn overflow");
          }
        }
        off += q->n;
      }
      for (t=0; t<ntypes; t++) {
        ta[l] = n - tb_off[c];
        for (j=0; j<tn; j++) 
          if (t==ty[j]) tb[n++] = nb[j];
        te[l] = n - tb_off[c];
        /* enlarge n to next multiple of 4 */
        nn = n % 4;
        n = nn ? n + 4 - nn : n;
        pa_max = MAX(pa_max,te[l]-ta[l]);
        if (n > nb_max-2*pa_max) {
          error("neighbor table full - increase nbl_size");
        }
        l++;
      }
    }
    
  }

  tb_off[ncells] = n;
  last_nbl_len   = n;
  have_valid_nbl = 1;
  nbl_count++;
}

/******************************************************************************
*
*  make_wp
*
******************************************************************************/

void make_wp(int k, wp_t *wp)
{
  cell_nbrs_t *c = cnbrs + k;
  int   m, j, i, n, l, t, min;
  flt   *pp;
  int   *tt;

  /* store dimensions */
  wp->k  = k;
  wp->n1 = cell_array[c->np].n;
  wp->n2 = 0;
  for (m=0; m<14; m++) {
    j = c->nq[m];
    if (j > -1) wp->n2 += cell_array[j].n;
  }

  /* allocate or enlarge wp */
  if (wp->n2 > wp->n2_max) {
    if (wp->n2_max > 0) {
      free(wp->pos);
      free(wp->force);
      free(wp->typ);
    }
    wp->n2_max = wp->n2 + 50;
    wp->pos    = (flt *) malloc( 4 * wp->n2_max * sizeof(flt) );
    wp->force  = (flt *) malloc( 4 * wp->n2_max * sizeof(flt) );
    wp->typ    = (int *) malloc(     wp->n2_max * sizeof(int) );
    if ((NULL==wp->pos) || (NULL==wp->force) || (NULL==wp->typ))
      error("cannot allocate workpackage");
  }

  /* copy position and type */
  pp = wp->pos;
  tt = wp->typ;
  n = 0; l = 0;
  for (m=0; m<14; m++) {
    j = c->nq[m];
    if (j == -1) continue;
    cell  *q = cell_array + j;
    for (i=0; i<q->n; i++) {
      pp[l++] = (flt) ORT(q,i,X);
      pp[l++] = (flt) ORT(q,i,Y);
      pp[l++] = (flt) ORT(q,i,Z);
      pp[l++] = (flt) 0.0;
      tt[n++] = SORTE(q,i);
    }
  }
  
  /* set pointers to neighbor tables */
  wp->ta  = ta + at_off[k];
  wp->te  = te + at_off[k];
  wp->tb  = tb + tb_off[k];
  wp->len = tb_off[k+1] - tb_off[k];

}

/******************************************************************************
*
*  calc_wp
*
******************************************************************************/

void calc_wp(wp_t *wp, int *is_short)
{
  flt d[4], f[4], r2, *pi, *pj, *fi, *fj, pot, grad;
  int i, l, c1, t, col, inc=pt.ntypes * pt.ntypes, n, j;

  wp->totpot = 0.0;
  wp->virial = 0.0;
  for (i=0; i<4*wp->n2; i++) wp->force[i] = 0.0;

  l = 0;
  for (i=0; i<wp->n1; i++) {

    c1 = pt.ntypes * wp->typ[i];
    fi = wp->force + 4*i;
    pi = wp->pos   + 4*i;

    for (t=0; t<pt.ntypes; t++) {

      col = c1 + t;

      for (n=wp->ta[l]; n<wp->te[l]; n++) {

        j    = wp->tb[n];
        pj   = wp->pos + 4*j;
        d[0] = pj[0] - pi[0];
        d[1] = pj[1] - pi[1];
        d[2] = pj[2] - pi[2];
        r2   = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];

        if (r2 <= pt.r2cut[4*col]) {

          flt tmp2, tmp6;
          tmp2 = pt.lj_sig[4*col] / r2;
          tmp6 = tmp2 * tmp2 * tmp2;
          pot  = pt.lj_eps[4*col] * tmp6 * (tmp6 - 2.0) - pt.lj_shift[4*col];
          grad = - 12.0 * pt.lj_eps[col] * tmp6 * (tmp6 - 1.0) / r2;

          //PAIR_INT(pot, grad, pair_pot, col, inc, r2, *is_short)

          wp->totpot += pot;
          pot        *= 0.5;   /* avoid double counting */
          wp->virial -= r2  * grad;

          f[0] = d[0] * grad;
          f[1] = d[1] * grad;
          f[2] = d[2] * grad;

          fi[0] += f[0];
          fi[1] += f[1];
          fi[2] += f[2];
          fi[3] += pot;

          fj = wp->force + 4*j;
          fj[0] -= f[0];
          fj[1] -= f[1];
          fj[2] -= f[2];
          fj[3] += pot;

	}
      }
      l++;
    }
  }
}

/******************************************************************************
*
*  store_wp
*
******************************************************************************/

void store_wp(wp_t *wp)
{
  int m, j, i, n=0;

  /* copy force and potential energy to cell array */
  for (m=0; m<14; m++) {
    j = cnbrs[wp->k].nq[m];
    if (j == -1) continue;
    cell *q = cell_array + j;
    for (i=0; i<q->n; i++) {
      KRAFT (q,i,X) += wp->force[n++];
      KRAFT (q,i,Y) += wp->force[n++];
      KRAFT (q,i,Z) += wp->force[n++];
      POTENG(q,i  ) += wp->force[n++];
    }
  }
  tot_pot_energy += wp->totpot;
  virial         += wp->virial;
}

/******************************************************************************
*
*  calc_forces
*
******************************************************************************/

void calc_forces(int steps)
{
  int  k, i, off, is_short=0;
  real tmpvec1[2], tmpvec2[2] = {0.0, 0.0};

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
  nfc++;

  /* clear per atom accumulation variables */
  for (k=0; k<nallcells; k++) {
    cell *p = cell_array + k;
    for (i=0; i<p->n; i++) {
      KRAFT(p,i,X) = 0.0;
      KRAFT(p,i,Y) = 0.0;
      KRAFT(p,i,Z) = 0.0;
      POTENG(p,i)  = 0.0;
    }
  }

  /* allocate wp if necessary */
  if (NULL==wp) {
    wp = (wp_t *) malloc( sizeof(wp_t) );
    if (NULL==wp) error("cannot allocate workpackage");
    wp->n2_max = 0;
  }

  /* compute interactions */
  for (k=0; k<ncells; k++) {
    make_wp(k, wp);
    calc_wp(wp, &is_short);
    store_wp(wp);
  }
  if (is_short) printf("short distance!\n");

#ifdef MPI
  /* sum up results of different CPUs */
  tmpvec1[0]     = tot_pot_energy;
  tmpvec1[1]     = virial;
  MPI_Allreduce( tmpvec1, tmpvec2, 2, REAL, MPI_SUM, cpugrid); 
  tot_pot_energy = tmpvec2[0];
  virial         = tmpvec2[1];
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
    for (i=0; i<p->n; i++) {
      d.x = ORT(p,i,X) - NBL_POS(p,i,X);
      d.y = ORT(p,i,Y) - NBL_POS(p,i,Y);
      d.z = ORT(p,i,Z) - NBL_POS(p,i,Z);
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
