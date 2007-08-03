
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

/* CBE declarations */
#include "imd_cbe.h"

int   *ti=NULL, *tb_off=NULL, *at_off=NULL, nb_max=0;
short *tb=NULL;
wp_t  *wp=NULL;
pt_t  pt;

/******************************************************************************
*
*  collect potential data for SPE
*
******************************************************************************/

void mk_pt(void)
{
  int i, j, k, col;
  flt r2, r6;

  /* Allocation size common to all mallocs */
  size_t const nelem = 4*ntypes*ntypes;
  size_t const asze  = nelem*(sizeof (flt));

  /* Set ntypes */
  pt.ntypes   = ntypes;

  /* Allocation modified by Frank Pister */
  /*
  pt.r2cut    = (flt *) malloc_aligned(asze, 16);
  pt.lj_sig   = (flt *) malloc_aligned(asze, 16);
  pt.lj_eps   = (flt *) malloc_aligned(asze, 16);
  pt.lj_shift = (flt *) malloc_aligned(asze, 16);
  */

  /* Allocate one array which is 4 times the size of every single array.
     As every individual array contains a number of items (nelem, see above)
     which is dividable by 4, the large array's number of elements is 
     dividable by 16.
     So its size is a multiple of 16bytes which is required for DMA.
     Besides that, we only need one call the memory allocation routine.

     In the DMAs it is assumed that pt.r2cut point to a block of memory
     aligned appropriatly and the the remaining pointers in pt_t point to
     some location inside that block.

     Note that, using this allocation scheme, lj_xxx must not be free'd!
   */
  pt.r2cut    = (flt*)(malloc_aligned(4*asze, 16,16));
  pt.lj_sig   = pt.r2cut  + nelem;
  pt.lj_eps   = pt.lj_sig + nelem;
  pt.lj_shift = pt.lj_eps + nelem;
  


  if ((NULL==pt.r2cut) || (NULL==pt.lj_sig) || (NULL==pt.lj_eps) || 
      (NULL==pt.lj_shift)) error("cannot allocate potential package");
  for (i=0; i<ntypes; i++)
    for (j=0; j<ntypes; j++) {
      col = i*ntypes + j;
      for (k=0; k<4; k++) { /* make vectors: 4 copies of each value */
        pt.r2cut  [4*col+k] = (flt) r2_cut      [i][j];
        pt.lj_sig [4*col+k] = (flt) SQR(lj_sigma[i][j]);
        pt.lj_eps [4*col+k] = (flt) lj_epsilon  [i][j] ;
        r2 = pt.lj_sig[4*col+k] / pt.r2cut[4*col+k];
        r6 = r2 * r2 * r2;
        pt.lj_shift[4*col+k] = pt.lj_eps[4*col+k] * r6 * (r6 - 2.0);
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
           (int)(nb_max * sizeof(int) / SQR(1024)) );
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
  int c, tn=1;
  int inc_short = 128 / sizeof(short);

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
      /* enlarge tn to next 128 byte boundary */
      tn = ((tn + inc_short - 1) / inc_short) * inc_short; 
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
  static int at_max=0, pa_max=0, ncell_max=0, cl_max=0;
  int  c, i, k, n, nn, at, cc, atm, l, t, j;
  int inc_int   = 128 / sizeof(int);
  int inc_short = 128 / sizeof(short);

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
  at=0;
  for (k=0; k<ncells; k++) {
    cell *p = CELLPTR(k);
    at_off[k] = at;
    /* next block on 128 byte boundary */
    at += ((2*p->n + inc_int - 1) / inc_int) * inc_int;
  }
  at_off[ncells] = at;

  /* (re-)allocate neighbor table */
  if (at >= at_max) {
    free(ti);
    at_max = (int) (1.1 * at);
    ti = (int *) malloc_aligned(at_max * sizeof(int), 128,16);
  }
  if (NULL==tb) {
    if (0==last_nbl_len) nb_max = (int) (nbl_size * estimate_nblist_size());
    else                 nb_max = (int) (nbl_size * last_nbl_len);
    tb = (short *) malloc_aligned(nb_max* sizeof(short), 0,0);
  }
  else if (last_nbl_len * sqrt(nbl_size) > nb_max) {
    free(tb);
    nb_max = (int    ) (nbl_size * last_nbl_len);
    tb     = (short *) malloc_aligned(nb_max* sizeof(short), 0,0);
  }
  if ((ti==NULL) || (tb==NULL)) 
    error("cannot allocate neighbor table");

  /* for all cells */
  nn=0;
  for (c=0; c<ncells; c++) {

    int   c1 = cnbrs[c].np;
    short *ttb;
    cell  *p = cell_array + c1;

    l = at_off[c];
    tb_off[c] = nn;
    ttb = tb + nn;
    n = 0;

    /* for each atom in cell */
    for (i=0; i<p->n; i++) {

      int    m, off;
      vektor d1;

      d1.x = ORT(p,i,X);
      d1.y = ORT(p,i,Y);
      d1.z = ORT(p,i,Z);

      /* for each neighboring atom */
      off = 0;
      ti[l] = n;
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
          if (r2 < cellsz) ttb[n++] = off + j;
        }
        off += q->n;
      }
      ti[l+1] = n - ti[l];
      l += 2;
      /* enlarge n to next 128 byte boundary */
      n = ((n + inc_short - 1) / inc_short) * inc_short; 
    }
    pa_max = MAX(pa_max,n);
    if (n + nn + 2*pa_max > nb_max) {
      error("neighbor table full - increase nbl_size");
    }
    nn += n;
  }

  tb_off[ncells] = nn;
  last_nbl_len   = nn;
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
    wp->pos    = (flt *) malloc_aligned(4 * wp->n2_max * sizeof(flt), 0,0);
    wp->force  = (flt *) malloc_aligned(4 * wp->n2_max * sizeof(flt), 0,0);
    wp->typ    = (int *) malloc_aligned(    wp->n2_max * sizeof(int), 0,0);
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
  wp->ti  = ti + at_off[k];
  wp->tb  = tb + tb_off[k];
  wp->len = tb_off[k+1] - tb_off[k];
  wp->n1_max  = 0;  /* max size for ti - nothing allocated here */
  wp->len_max = 0;  /* max size for tb - nothing allocated here */
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
  int  k, i, off;
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

    /* Debugging output */
    /* printf("wp on PPU after make_wp:\n");  wp_out(printf, wp); */

    /* calc_wp(wp); */
    schedtospu(wp); 

    /* Debugging output */
    /* printf("wp on PPU before store_wp:\n");  wp_out(printf, wp); */

    store_wp(wp);
  }

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
