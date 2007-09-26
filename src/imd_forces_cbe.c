
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


/* PPU intrinsics */
#include <ppu_intrinsics.h>

/* IMD declarations */
#include "imd.h"
/* CBE declarations */
#include "imd_cbe.h"

/* Global variables only needed in this file */
static int *ti=NULL;
static int *tb_off=NULL;
static int *at_off=NULL;
static int nb_max=0;
static short *tb=NULL;
#ifdef CBE_DIRECT
static int n_max=0, len_max=0;
#ifdef ON_PPU
static cell_dta_t *cell_dta=NULL;
#else
static cell_ea_t  *cell_dta=NULL;
#endif
#endif

/* Potential which is need by calc_wp */
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

  /* Allocate one array which is 4 times the size of a single array */
  pt.r2cut    = (flt*)(malloc_aligned(4*asze, 128, 16));
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

#ifdef CBE_DIRECT

/******************************************************************************
*
*  estimate_nblist_size
*
******************************************************************************/

int estimate_nblist_size(void)
{
  int c, tn=1;
  int const inc_short = 128 / sizeof(short) - 1;

  /* for all cells */
  for (c=0; c<ncells; c++) {

    int m, c1 = cnbrs[c].np, tn2 = 1, i;
    cell *p   = cell_array + c1;

    /* for each neighboring cell */
    for (m=0; m<NNBCELL; m++) {

      /* for each atom in first cell */
      for (i=0; i<p->n; i++) {

        int    c2, jstart, j;
        real   r2;
        cell   *q;
        vektor d1;

        d1.x = ORT(p,i,X);
        d1.y = ORT(p,i,Y);
        d1.z = ORT(p,i,Z);

        c2 = cnbrs[c].nq[m];
        if (c2<0) continue;
#ifdef AR
        jstart = (c2==c1) ? i+1: 0;
#else
        jstart = 0;
#endif
        q = cell_array + c2;

        /* for each atom in neighbor cell */
        for (j=jstart; j<q->n; j++) {
          vektor d;
          d.x = ORT(q,j,X) - d1.x;
          d.y = ORT(q,j,Y) - d1.y;
          d.z = ORT(q,j,Z) - d1.z;
          r2  = SPROD(d,d);
          if (r2 < cellsz) tn2++;
        }
        /* enlarge tn2 to next value divisible by 4 */
        tn2 = (tn2 + 3) & (~3); 
      }
      /* enlarge tn to next 128 byte boundary */
      tn2 = (tn2 + inc_short) & (~inc_short); 
    }
    tn  = MAX(tn2,tn);
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
  int c, k, nn, at, i;
  int inc_int   = 128 / sizeof(int)   - 1;
  int inc_short = 128 / sizeof(short) - 1;

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
    free(cell_dta);
    cl_max   = ncells;
#ifdef ON_PPU
    cell_dta = (cell_dta_t *) malloc( (nallcells+1) * sizeof(cell_dta_t) );
#else
    cell_dta = (cell_ea_t  *) malloc( (nallcells+1) * sizeof(cell_ea_t)  );
#endif
    at_off = (int *) malloc( (cl_max+1)           * sizeof(int) );
    tb_off = (int *) malloc( (cl_max+1) * NNBCELL * sizeof(int) );
    if ((NULL==at_off) || (NULL==tb_off) || (NULL==cell_dta))
      error("cannot allocate neighbor table");
  }

  /* initialize cell_dta */
  for (k=0; k<nallcells; k++) {
    cell *p = cell_array + k;
#ifdef ON_PPU
    cell_dta[k].n     = p->n;
    cell_dta[k].pos   = p->ort;
    cell_dta[k].force = p->kraft;
    cell_dta[k].typ   = p->sorte;
    cell_dta[k].ti    = NULL;
    cell_dta[k].tb    = NULL;
#else
    cell_dta[k].n = p->n;
    PTR2EA( p->ort,   cell_dta[k].pos_ea   );
    PTR2EA( p->kraft, cell_dta[k].force_ea );
    PTR2EA( p->sorte, cell_dta[k].typ_ea   );
    PTR2EA( NULL,     cell_dta[k].ti_ea    );
    PTR2EA( NULL,     cell_dta[k].tb_ea    );
#endif
  }

  /* count atom numbers */
  at=0; n_max=0;
  for (k=0; k<ncells; k++) {
    cell *p = CELLPTR(k);
    at_off[k] = at;
    /* next block on 128 byte boundary */
    at += ((2*p->n + inc_int) & (~inc_int)) * NNBCELL;
    n_max = MAX(n_max,p->n);
  }
  at_off[ncells] = at;
  n_max++;

  /* (re-)allocate neighbor table */
  if (at >= at_max) {
    free(ti);
    at_max = (int) (1.1 * at * NNBCELL);
    ti = (int *) malloc_aligned(at_max * sizeof(int), 128, 16);
  }
  if (NULL==tb) {
    if (0==last_nbl_len) 
      nb_max = (int) (nbl_size * ncells * estimate_nblist_size());
    else                 
      nb_max = (int) (nbl_size * last_nbl_len);
    tb = (short *) malloc_aligned( nb_max * sizeof(short), 0, 0);
  }
  else if (last_nbl_len * sqrt(nbl_size) > nb_max) {
    free(tb);
    nb_max = (int    ) (nbl_size * last_nbl_len);
    tb     = (short *) malloc_aligned( nb_max * sizeof(short), 0, 0);
  }
  if ((ti==NULL) || (tb==NULL)) 
    error("cannot allocate neighbor table");

  /* for all cells */
  nn=0; len_max=0;
  for (c=0; c<ncells; c++) {

    int   c1 = cnbrs[c].np, m;
    cell  *p = cell_array + c1;
    int   at_inc = (2*p->n + inc_int) & (~inc_int);

    /* for each neighbor cell */
    for (m=0; m<NNBCELL; m++) {

      int   c2 = cnbrs[c].nq[m];
      int   l  = at_off[c] + m * at_inc, n = 0, i;
      cell  *q;
      short *ttb = tb + nn;

      tb_off[c*NNBCELL+m] = nn;

      c2 = cnbrs[c].nq[m];
      if (c2<0) continue;
      q = cell_array + c2;

      /* for each atom in cell */
      for (i=0; i<p->n; i++) {

        int    jstart, j, rr;
        vektor d1;

        d1.x = ORT(p,i,X);
        d1.y = ORT(p,i,Y);
        d1.z = ORT(p,i,Z);

        /* for each neighboring atom */
        ti[l] = n;
#ifdef AR
        jstart = (m==0) ? i+1 : 0;
#else
        jstart = 0;
#endif
        for (j=jstart; j<q->n; j++) {
          vektor d;
          real   r2;
#ifndef AR
          if ((m==0) && (i==j)) continue;
#endif
          d.x = ORT(q,j,X) - d1.x;
          d.y = ORT(q,j,Y) - d1.y;
          d.z = ORT(q,j,Z) - d1.z;
          r2  = SPROD(d,d);
          if (r2 < cellsz) ttb[n++] = j;
        }
        ti[l+1] = n - ti[l];
        l += 2;

        /* if n is not divisible by 4, pad with copies of q->n */
        rr = n % 4;
        if (rr>0) for (j=rr; j<4; j++) ttb[n++] = q->n;
      }

      /* enlarge n to next 128 byte boundary */
      n = (n + inc_short) & (~inc_short); 

      nn += n;
      pa_max = MAX(pa_max,n);
      if (nn + 2*pa_max > nb_max) {
        error("neighbor table full - increase nbl_size");
      }
      len_max = MAX( len_max, nn - tb_off[c*NNBCELL+m] );
    }
  }

  tb_off[ncells*NNBCELL] = nn;
  last_nbl_len   = nn;
  have_valid_nbl = 1;
  nbl_count++;
}


/******************************************************************************
*
*  make_wp
*
******************************************************************************/

static wp_t* make_wp(int k, wp_t *wp)
{
  cell_nbrs_t *c = cnbrs + k;
  cell *p = CELLPTR(k);
  int  inc_int = 128 / sizeof(int) - 1, m, n;
  int  at_inc  = (2*p->n + inc_int) & (~inc_int);

  wp->k = k;
  wp->n_max = n_max;      /* allocate for this many atoms     */
  wp->len_max = len_max;  /* allocate for this many neighbors */
  wp->ti_len = at_inc;
  for (m=0; m<NNBCELL; m++) {
    n = c->nq[m];
    if (n<0) {
      wp->cell_dta[m].n = 0;
      continue;
    }
    wp->cell_dta[m] = cell_dta[n];
#ifdef ON_PPU
    wp->cell_dta[m].ti = ti + at_off[k] + m * at_inc;
    wp->cell_dta[m].tb = tb + tb_off[k * NNBCELL + m];
#else
    wp->cell_dta[m].len = tb_off[k*NNBCELL+m+1] - tb_off[k*NNBCELL+m];
    PTR2EA( ti + at_off[k] + m * at_inc,  wp->cell_dta[m].ti_ea );
    PTR2EA( tb + tb_off[k * NNBCELL + m], wp->cell_dta[m].tb_ea );
#endif
  }


  /* Debugging output */
  /*
  fprintf(stdout, "wp k=%d created for CBE_DIRECT with n_max=%d, len_max=%d\n",
          k, wp->n_max, wp->len_max);
  fflush(stdout);
  */

  return wp;
}


/******************************************************************************
*
*  store_wp
*
******************************************************************************/

static INLINE_ void store_wp(wp_t const* const wp)
{
  tot_pot_energy += wp->totpot;
  virial         += wp->virial;
}


#else  /* not CBE_DIRECT */

/******************************************************************************
*
*  estimate_nblist_size
*
******************************************************************************/

int estimate_nblist_size(void)
{
  int c, tn=1;
  int inc_short = 128 / sizeof(short) - 1;

  /* for all cells */
  for (c=0; c<ncells; c++) {

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
      tn = (tn + inc_short) & (~inc_short); 
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
  int  c, i, k, n, nn, at, l, t, j;
  int inc_int   = 128 / sizeof(int)   - 1;
  int inc_short = 128 / sizeof(short) - 1;

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
    at += (2*p->n + inc_int) & (~inc_int);
  }
  at_off[ncells] = at;

  /* (re-)allocate neighbor table */
  if (at >= at_max) {
    free(ti);
    at_max = (int) (1.1 * at);
    ti = (int *) malloc_aligned(at_max * sizeof(int), 128,16);
  }
  if (NULL==tb) {
    if (0==last_nbl_len) 
      nb_max = (int) (nbl_size * ncells * estimate_nblist_size());
    else                 
      nb_max = (int) (nbl_size * last_nbl_len);
    tb = (short *) malloc_aligned(nb_max* sizeof(short), 0,0);
  }
  else if (last_nbl_len * sqrt(nbl_size) > nb_max) {
    free(tb);
    nb_max = (int    ) (nbl_size * last_nbl_len);
    tb     = (short *) malloc_aligned( nb_max * sizeof(short), 0, 0);
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

      int    m, off, rr;
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

      /* if n is not divisible by 4, pad with copies of i */
      rr = n % 4;
      if (rr>0) for (j=rr; j<4; j++) ttb[n++] = i;

      /* enlarge n to next 128 byte boundary */
      n = (n + inc_short) & (~inc_short); 

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

static wp_t* make_wp(int k, wp_t *wp)
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

  /* Return ptr to the wp just created */
  return wp;
}


/******************************************************************************
*
*  store_wp
*
******************************************************************************/

static void store_wp(wp_t const* const wp)
{
  /* Some indices */
  int m, j, i, n=0;

  /* Deref. pointer */
  flt const* const force = wp->force;

  /* copy force and potential energy to cell array */
  for (m=0; m<14; m++) {
    j = cnbrs[wp->k].nq[m];
    if (j == -1) continue;
    cell *q = cell_array + j;
    for (i=0; i<q->n; i++) {
      KRAFT (q,i,X) += force[n++];
      KRAFT (q,i,Y) += force[n++];
      KRAFT (q,i,Z) += force[n++];
      POTENG(q,i  ) += force[n++];
    }
  }
  tot_pot_energy += wp->totpot;
  virial         += wp->virial;
}

#endif  /* CBE_DIRECT */


/******************************************************************************
*
*  calc_forces_spu
*
******************************************************************************/

static void calc_forces_spu(int const steps)
{
  /* The following code uses mailboxes to synchronize between PPU & SPUs.
     To do so, it uses direct (memory mapped) access to the control area
     of the SPUs.
     Writing/reading to the members of the control block struct assume
     that writes and reads to unsigneds are atomic.
  */
  /* Furthermore, we need the following bit masks: */
  enum { OUTMBX_CNT_MASK  = 0x000000ffu,
         INMBX_CNT_MASK   = 0x0000ff00u,
         OUTMBX_CNT_SHIFT = 0u,
         INMBX_CNT_SHIFT  = 8u
       };

  /* Some indices */
  int k, kdone, i, ispu;
  int ispumax = num_spus;


  /* The following array keeps track of the states of the SPUs */
  typedef enum { IDLE=0u, WORKING=1u } Tspustate;
  Tspustate spustate[N_SPU_THREADS_MAX];

  /* Some statistics: spuload[k] contains the number of wps which have
     been sheduled to SPU k  */
  unsigned spuload[N_SPU_THREADS_MAX];


  /* Assume that all SPUs are idle and are awaiting a mailbox message 
     when entering this routine */
  for ( ispu=0;    (ispu<ispumax);   ++ispu ) {
      spustate[ispu]=IDLE;
      spuload[ispu]=0;
  }

  /* fprintf(stdout, "nallcells=%i  ncells=%i\n", nallcells, ncells); fflush(stdout); */

  /* While there is still work to be scheduled or there are still 
     results to be picked up. */
  for (  k=kdone=0;       (kdone<ncells) || (k<ncells);      )  {

      /* fprintf(stdout, "k=%u  kdone=%u\n", k,kdone); */

      /* Index for wp/exch */
      unsigned iwp1;
      for (  ispu=iwp1=0;     (ispu<ispumax);     ++ispu, iwp1+=N_BUFLEV  ) {

           /* Control area used to access the mailbox and ptr. to state */
  	   spe_spu_control_area_p const pctl = cbe_spucontrolarea[ispu];
           Tspustate* const pstat = spustate+ispu;

           /* Work package */
           wp_t* const pwp = cbe_wp+iwp1;
#if ! defined(CBE_DIRECT)
           /* Exchange buffer needed if no direct access enabled */
           exch_t* const pexch = cbe_exch+iwp1;
#endif


           /* fprintf(stdout, "SPU %u is %s\n", ispu, ((WORKING==*pstat) ? "working" : "notworking"));  */
 

           /* Idle SPU and still work to schedule */
           if (  (IDLE==*pstat) && (k<ncells)  ) {
   	       /* Create a work package... */
   	       make_wp(k, pwp);
#if ! defined (CBE_DIRECT)
               /* ...and also create an exch buffer */
               create_exch(pwp, pexch);
#endif

 
               /* Signal SPU by writing to its mailbox when space is available */
               __lwsync();  /* synchronize memory before releasing it */

               while ( 0u == ((pctl->SPU_Mbox_Stat) & INMBX_CNT_MASK) ) {}
               pctl->SPU_In_Mbox = WPSTRT1;

	       /* Mark SPU as working */
	       ++k;
  	       *pstat=WORKING;

               /* Update statistics */
               ++(spuload[ispu]);

               /* This SPU is working, so we can move on to the next one */
               continue;
           }


           /* Has some wp been scheduled to current SPU and are there
              still results to be picked up? */
           if (  (WORKING==*pstat) && (kdone<ncells) ) {

	        /* Check wether SPU has finished by reading its mailbox */
	        if ( 0u != ((pctl->SPU_Mbox_Stat) & OUTMBX_CNT_MASK) ) {

		    /* Data available, so read it */
		    unsigned const spumsg = pctl->SPU_Out_Mbox;

                    /* fprintf(stdout, "%u from spu %u\n", spumsg, ispu); */

                    __lwsync();  /* synchronize memory before accessing it */

                    if ( WPDONE1 == spumsg ) {
		        /* Copy results from exch buffer */
#if ! defined(CBE_DIRECT)
		        pwp->totpot = pexch->totpot;
                        pwp->virial = pexch->virial;
#endif

                        /* Store the wp */
                        store_wp(pwp);

  		        /* SPU's done with its work, so it's idle again 
                           and one more wp is finished */
                        ++kdone;
		        *pstat=IDLE;                      
                    }
                }
                else {
		   /* This SPU is working and has not yet finished, so move
                      on to the next one */
 		   continue;
                }
           }
      }
  }

  /* Print SPU summary */
  /*
  for (  ispu=0;   (ispu<ispumax);    ++ispu  ) {
      fprintf(stdout, "%u  ", (unsigned)(spustate[ispu]));
  }
  fputc('\n',stdout);
  for (  ispu=0;   (ispu<ispumax);    ++ispu ) {
      fprintf(stdout, "%u  ", spuload[ispu]);
  }
  fputc('\n',stdout);
  */
}


/******************************************************************************
*
*  calc_forces_ppu
*
******************************************************************************/

static void calc_forces_ppu(int const steps)
{
  /* (Global) work package */
  wp_t* const wp = cbe_wp;
  int k;

  /* compute interactions */
  for (k=0; k<ncells; ++k ) {
    make_wp(k, wp);
    calc_wp(wp);
    store_wp(wp);
  }
}


/******************************************************************************
*
*  calc_forces
*
******************************************************************************/

void calc_forces(int const steps)
{
  int  k, i;
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

#if defined(ON_PPU) || defined(AR)
  /* clear per atom accumulation variables */
  for (  k=0;    k<nallcells;  ++k ) {
    cell *p = cell_array + k;
    for (i=0; i<p->n; i++) {
      KRAFT(p,i,X) = 0.0;
      KRAFT(p,i,Y) = 0.0;
      KRAFT(p,i,Z) = 0.0;
      POTENG(p,i)  = 0.0;
    }
  }
#endif

  /* compute the forces - either on PPU or on SPUs*/
#ifdef ON_PPU
  calc_forces_ppu(steps); 
#else
  calc_forces_spu(steps);
#endif

#ifdef MPI
  /* sum up results of different CPUs */
  tmpvec1[0]     = tot_pot_energy;
  tmpvec1[1]     = virial;
  MPI_Allreduce( tmpvec1, tmpvec2, 2, REAL, MPI_SUM, cpugrid); 
  tot_pot_energy = tmpvec2[0];
  virial         = tmpvec2[1];
#endif

#ifdef AR
  /* add forces back to original cells/cpus */
  send_forces(add_forces,pack_forces,unpack_forces);
#endif

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
