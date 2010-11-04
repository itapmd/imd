
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2008 Institute for Theoretical and Applied Physics,
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

/* Needed for multithreaded "SPU driver" */
#include <pthread.h>

/* PPU intrinsics */
#include <ppu_intrinsics.h>

/* IMD declarations */
#include "imd.h"
/* CBE declarations */
#include "imd_cbe.h"

/* Potential data structure */
pt_t ALIGNED_(16,pt);

/* extra data for integrator on SPU */
#ifdef SPU_INT
mv_t ALIGNED_(16,mv);
#endif

/* Global variables only needed in this file */
static int *ti= NULL;
static int *tb_off=NULL;
static int *at_off=NULL;
static short *tb=NULL;
static int nb_max=0;
static int n_max=0, len_max=0, len_max2=0;
#ifdef ON_PPU
static cell_dta_t *cell_dta=NULL;
#else
static cell_ea_t  *cell_dta=NULL;
#endif


/* Constants to access the SPU mailboxes */
/* Shifts */
enum MBXSHIFT {
    OUTMBX_CNT_SHIFT = (unsigned)0,  /* outbox */
    INMBX_CNT_SHIFT  = (unsigned)8   /* inbox  */
};
/* Masks */
enum MBXMASK { 
    OUTMBX_CNT_MASK  = (unsigned)0x000000ff,  /* outbox */
    INMBX_CNT_MASK   = (unsigned)0x0000ff00   /* inbox */
};


/* Read Nvpmbx values per mbox */
static void mboxes2tuples(spe_spu_control_area_p const* const cfrst, 
                          unsigned const Nmbx, unsigned const Nvpmbx,
                          unsigned* const res,
                          unsigned* const cnt)
{
    /* Iterators */
    register unsigned* icnt;
    register spe_spu_control_area_p const* ic;
    register unsigned* ires;
    /* Mailbox index, #of reads remaining */
    register unsigned imbx;
    register unsigned nrem;

    for ( icnt=cnt, imbx=Nmbx;   (imbx>0u);    ++icnt, --imbx ) {
        (*icnt)=0u;
    }

    for ( nrem=Nmbx*Nvpmbx;    (nrem>0u);    ) {
        for ( ic=cfrst, icnt=cnt, ires=res, imbx=Nmbx; (imbx>0u); 
              ++ic, ires+=Nvpmbx,  ++icnt, --imbx ) {
	   /* Have to read? */
	   if ( (*icnt)<Nvpmbx ) {
	       /* Mask */
	       enum { OMSK = (unsigned)0x000000ff };

	       /* Can read from mbox ? */
	       if ( 0u != (((*ic)->SPU_Mbox_Stat) & OMSK)  ) {
 		   *(ires+(*icnt)) = (*ic)->SPU_Out_Mbox;
                   ++(*icnt);
                   --nrem;
               }
           }
        }
    }
}


static void value2mboxes(spe_spu_control_area_p const* const cfrst, 
                         unsigned const Nmbx,
                         unsigned const v,
                         unsigned* const cnt)
{
    /* Iterators */
    register unsigned* icnt;
    register spe_spu_control_area_p const* ic;
    /* Mailbox index, #of writes remaining */
    register unsigned imbx;
    register unsigned nrem;

    for ( imbx=Nmbx, icnt=cnt;     (imbx>0u);      ++icnt, --imbx  ) {
        (*icnt) = 1u;
    }

    for ( nrem=Nmbx;    (nrem>0u);    ) {
        for (ic=cfrst, icnt=cnt, imbx=Nmbx; (imbx>0u); ++ic, ++icnt, --imbx) {
	   /* Have to write? */
	   if ( 1u == (*icnt) ) {
	       /* Mask */
	       enum { IMSK = (unsigned)0x0000ff00 };

	       /* Can write to mbox ? */
	       if ( 0u != (((*ic)->SPU_Mbox_Stat) & IMSK)  ) {
		   (*ic)->SPU_In_Mbox =  v;
                   *icnt = 0u;
                   --nrem;
               }
           }
        }
    }
}


static void tuples2mboxes(spe_spu_control_area_p const* const cfrst, 
                          unsigned const Nmbx, unsigned const Nvpmbx,
                          unsigned const* const val,
                          unsigned      * const cnt)
{
    /* Iterators */
    register unsigned* icnt;
    register spe_spu_control_area_p const* ic;
    register unsigned const* ival = val;
    /* Mailbox index, #of writes remaining */
    register unsigned imbx;
    register unsigned nrem;

    for ( icnt=cnt, imbx=Nmbx;   (imbx>0u);    ++icnt, --imbx ) {
        (*icnt)=0u;
    }

    for ( nrem=Nmbx*Nvpmbx;    (nrem>0u);    ) {
        for ( ic=cfrst, ival=val, icnt=cnt, imbx=Nmbx;     (imbx>0u);    ++ic, ival+=Nvpmbx, ++icnt, --imbx ) {
	   /* Have to write? */
	   if ( (*icnt)<Nvpmbx ) {
	       /* Mask */
	       enum { IMSK = (unsigned)0x0000ff00 };

	       /* Can write to mbox ? */
	       if ( 0u != (((*ic)->SPU_Mbox_Stat) & IMSK)  ) {
		   (*ic)->SPU_In_Mbox =  *(ival+(*icnt));
                   ++(*icnt);
                   --nrem;
               }
           }
        }
    }
}


/******************************************************************************
*
*  collect potential data for SPE
*
******************************************************************************/

#ifdef LJ

void mk_pt(void)
{
  int   i, j, k, col;
  float r2, r6;

  /* Allocation size common to all mallocs */
  size_t const nelem = 4*ntypes*ntypes;
  size_t const asze  = nelem*(sizeof(float));

  /* Set ntypes */
  pt.ntypes   = ntypes;

  /* Allocate one array which is 4 times the size of a single array */
  pt.r2cut    = (float *) malloc_aligned(4*asze, 128, 16);
  pt.lj_sig   = pt.r2cut  + nelem;
  pt.lj_eps   = pt.lj_sig + nelem;
  pt.lj_shift = pt.lj_eps + nelem;
  PTR2EA( pt.r2cut, pt.ea );
#ifdef SPU_INT
  PTR2EA( &mv, pt.mv_ea );
#endif

  if (NULL==pt.r2cut) error("cannot allocate potential package");
  for (i=0; i<ntypes; i++)
    for (j=0; j<ntypes; j++) {
      col = i*ntypes + j;
      for (k=0; k<4; k++) { /* make vectors: 4 copies of each value */
        pt.r2cut  [4*col+k] = (float) r2_cut      [i][j];
        pt.lj_sig [4*col+k] = (float) SQR(lj_sigma[i][j]);
        pt.lj_eps [4*col+k] = (float) lj_epsilon  [i][j] ;
        r2 = pt.lj_sig[4*col+k] / pt.r2cut[4*col+k];
        r6 = r2 * r2 * r2;
        pt.lj_shift[4*col+k] = pt.lj_eps[4*col+k] * r6 * (r6 - 2.0);
      }
    }
}

#else

/* cubic spline interpolation */
#define PAIR_INT(pot,grad,r2,col)  {                                         \
                                                                             \
  float r2a, a, b, a2, b2, istep, step, st6, *y;                             \
  int   k, col4=col*4;                                                       \
                                                                             \
  /* indices into potential table */                                         \
  istep = pt.invstep[col4];                                                  \
  step  = pt.step[col4];                                                     \
  r2a   = MIN(r2,pt.r2cut[col4]);                                            \
  r2a   = r2a * istep;                                                       \
  r2a   = r2a - pt.begin[col4];                                              \
  k     = MAX(0,(int) r2a);                                                  \
  b     = r2a - k;                                                           \
  a     = 1.0 - b;                                                           \
  a2    = a * a - 1;                                                         \
  b2    = b * b - 1;                                                         \
  st6   = step / 6;                                                          \
  y     = pt.tab[col] + 4 * k;                                               \
                                                                             \
  /* potential and twice the derivative */                                   \
  pot  = a * y[0] + b * y[2] + (a * a2 * y[1] + b * b2 * y[3]) * st6 * step; \
  grad = 2*((y[2] - y[0])*istep + ((3*b2+2) * y[3] - (3*a2+2) * y[1])* st6); \
}

void mk_pt(void)
{
  float *u    = NULL, r2;
  int   n     = cbe_pot_steps, off = 0, i, j, k, m;
  int   inc   = ntypes * ntypes;
  int   hsize = 4 * ntypes * ntypes;
  int   psize = 2 * (n+1) * ntypes * (ntypes+1);
  int   tsize = 4 * hsize + psize;

  /* allocate data */
  pt.data = (float *) malloc_aligned( tsize * sizeof(float), 128, 16 );
  u       = (float *) malloc( (n+1) * sizeof(float) );
  if ((NULL==pt.data) || (NULL==u)) error("cannot allocate potential package");

  /* fill in header of potential table */
  PTR2EA( pt.data, pt.ea );
#ifdef SPU_INT
  PTR2EA( &mv, pt.mv_ea );
#endif
  pt.ntypes  = ntypes;
  pt.nsteps  = n;
  pt.begin   = pt.data + off;  off += hsize;
  pt.r2cut   = pt.data + off;  off += hsize;
  pt.step    = pt.data + off;  off += hsize;
  pt.invstep = pt.data + off;  off += hsize;
  for (i=0; i<ntypes; i++)
    for (j=i; j<ntypes; j++) {
      pt.tab[i*ntypes+j] = pt.tab[j*ntypes+i] = pt.data + off;  
      off += 4*(n+1);
    }
  for (k=0; k<ntypes*ntypes; k++) {

    /* determine where to start tabulation */
    m = 0;
    while (pair_pot.table[m*inc+k] > cbe_pot_max) m++;
    r2 = pair_pot.begin[k] + m * pair_pot.step[k];

    /* fill in potential table header */
    for (i=0; i<4; i++) {
      int col = 4*k+i;
      pt.r2cut  [col] = pair_pot.end[k];
      pt.step   [col] = (pt.r2cut[col] - r2) / n;
      pt.invstep[col] = 1.0 / pt.step[col];
      pt.begin  [col] = r2 * pt.invstep[col];
    }
  }

  /* loop over columns */
  for (i=0; i<ntypes; i++)
    for (j=i; j<ntypes; j++) {

      int   col  = i*ntypes + j, idummy;
      float *y   = pt.tab[col], *y2 = y + 1;
      float step = pt.step[col];
      float r2, p, qn, un, q0, fdummy;

      /* fill in the potential values */
      for (k=0; k<n; k++) {
        r2 = (pt.begin[col] + k) * step;
	pair_int2(y+4*k, &fdummy, &idummy, &pair_pot, col, inc, r2);
      }
      y[4*n] = 0.0;  /* last potential value is zero */

      /* set first derivative at the left end correctly */
      y2[0] = -0.5;
      r2    = pt.begin[col] * step; 
      pair_int2( &fdummy, &q0, &idummy, &pair_pot, col, inc, r2);
      u[0]  = (3.0 / step) * ( (y[4]-y[0]) / step - q0 * 0.5 );

      for (k=1; k<n; k++) {
        p = 0.5 * y2[4*(k-1)] + 2.0;
        y2[4*k] = -0.5 / p;
        u[k] = (y[4*(k+1)] - 2*y[4*k] + y[4*(k-1)]) / step; 
        u[k] = (6.0 * u[k] / (2*step) - 0.5 * u[k-1]) / p;
      }

      /* set first derivative zero at right end */
      qn = 0.5;
      un = (3.0 / step) * (y[4*(n-1)] - y[4*n]) / step;

      y2[4*n] = (un - qn * u[n-1]) / (qn * y2[4*(n-1)] + 1.0);
      for (k=n-1; k>=0; k--) 
        y2[4*k] = y2[4*k] * y2[4*(k+1)] + u[k];

      /* fill in second copies of y and y2 */
      for (k=0; k<=4*(n-1); k+=4) {
        y [k+2] = y [k+4];
        y2[k+2] = y2[k+4];
      }

      /* for security, we continue the last interpolation polynomial */
      y [4*n+2] = 2*y [4*n]-y [4*(n-1)]+SQR(step)*y2[4*n];
      y2[4*n+2] = 2*y2[4*n]-y2[4*(n-1)];

    }

  free(u);

  if (1==debug_potential) {

    /* potential as interpolated from SPU potential table */
    for (i=0; i<ntypes*ntypes; i++) {
      FILE *out;
      char fname[128];
      float r2, pot, grad;
      float start = pt.begin[4*i] * pt.step[4*i];
      float dr2 = (pt.r2cut[4*i] - start) / 10000;
      sprintf(fname, "cbe_pot.%d", i);
      out = fopen(fname,"w");
      for (k=0; k<=10000; k++) {
        r2 = start + k * dr2;
        PAIR_INT(pot,grad,r2,i);
        fprintf(out,"%f %f %f\n", r2, pot, grad);
      }
      fclose(out);
    }

    /* SPU potential table */
    for (i=0; i<ntypes*ntypes; i++) {
      FILE *out;
      char fname[128];
      float r2, pot, grad, *f = pt.tab[i];
      sprintf(fname, "pt.%d", i);
      out = fopen(fname,"w");
      for (k=0; k<=pt.nsteps; k++) {
        r2 = (pt.begin[4*i] + k) * pt.step[4*i];
        fprintf(out,"%f %f %f %f %f\n",r2,f[4*k],f[4*k+1],f[4*k+2],f[4*k+3]);
      }
      fclose(out);
    }

  }
}

#endif


/******************************************************************************
*
*  prepare extra data for integrator on SPU
*
******************************************************************************/

#ifdef SPU_INT

void mk_mv(void)
{
  int i, t;
#ifndef NVT
  float eta=0.0;
#endif

  float fric  =        1.0 - eta * timestep / 2.0;
  float ifric = 1.0 / (1.0 + eta * timestep / 2.0);
  float fac1  = fric * ifric;
  float fac2  = timestep * ifric;

  mv.ts  [0] = mv.ts  [1] = mv.ts  [2] = timestep; mv.ts  [3] = 0.0;
  mv.fac1[0] = mv.fac1[1] = mv.fac1[2] = fac1;     mv.fac1[3] = 0.0;
  mv.fac2[0] = mv.fac2[1] = mv.fac2[2] = fac2;     mv.fac2[3] = 0.0;
#ifdef MIK
  if (ensemble == ENS_MIK) 
    mv.mikmask[0] = mv.mikmask[1] = mv.mikmask[2] = 0.0;  mv.mikmask[3] = 0.0;
  else
    mv.mikmask[0] = mv.mikmask[1] = mv.mikmask[2] = 1e10; mv.mikmask[3] = 0.0;
#endif
  for (i=0; i<4; i++) {
    for (t=0; t<ntypes; t++) mv.imass[4*t+i] = 1.0 / masses[t];
  }
  for (t=0; t<vtypes; t++) {
    mv.restr[4*t  ] = restrictions[t].x;
    mv.restr[4*t+1] = restrictions[t].y;
    mv.restr[4*t+2] = restrictions[t].z;
    mv.restr[4*t+3] = 1.0;  /* here, the energy is stored! */
#ifdef FBC
    mv.fbc_f[4*t  ] = fbc_forces[t].x;
    mv.fbc_f[4*t+1] = fbc_forces[t].y;
    mv.fbc_f[4*t+2] = fbc_forces[t].z;
    mv.fbc_f[4*t+3] = 0.0;  /* here, the energy is stored! */
#endif
  }
}


#endif

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
*  make_tb
*
******************************************************************************/

static void make_tb(int k, wp_t *wp)
{
  int m;
  /* cell *p = CELLPTR(k); */
  cell_nbrs_t const * const c = cnbrs + k;
  int const inc_int = 128 / sizeof(int) - 1;
  int const at_inc  = (2*(CELLPTR(k))->n + inc_int) & (~inc_int);

  wp->totpot = cellsz;          /* we abuse totpot for cellsz here */
  wp->nb_max = nb_max;
  wp->k = k;
  wp->n_max = n_max;            /* allocate for this many atoms     */
  wp->len_max = n_max * n_max;  /* allocate for this many neighbors */
  wp->ti_len = at_inc;
  wp->flag = 1;


  for (  m=0;    EXPECT_FALSE_(m<NNBCELL);  ++m ) {
    int const n = c->nq[m];
    if (n<0) {
      wp->cell_dta[m].n = 0;
      continue;
    }
    wp->cell_dta[m] = cell_dta[n];
#ifdef ON_PPU
    wp->cell_dta[m].ti = ti + at_off[k] + m * at_inc;
#else
    PTR2EA( ti + at_off[k] + m * at_inc,  wp->cell_dta[m].ti_ea );
#endif
  }
#ifdef ON_PPU
  wp->cell_dta[0].tb = tb + tb_off[k * NNBCELL];
#else
  PTR2EA( tb + tb_off[k * NNBCELL], wp->cell_dta[0].tb_ea );
#endif
}


static void make_tb_args(argbuf_t* parg, unsigned n, int k)
{
     for ( ;     n>0u;    ++parg, ++k, --n ) {
         make_tb(k, ((wp_t*)parg));
     }

     /* Sync memory */
     __lwsync();
}


/******************************************************************************
*
*  store_tb
*
******************************************************************************/

static void store_tb(wp_t const* const wp)
{
  register int off, off2=0, m;
  int const k = wp->k;

  /* Some "base" pointer */
  register int* const tb_off_Nk = tb_off + NNBCELL*k;

  if (wp->flag<0) {
       error("error flag in store_tb");
  }

  for ( m=1;  EXPECT_TRUE_(m<NNBCELL);   ++m ) {
#ifdef ON_PPU
    off = (int) (wp->cell_dta[m].tb       - wp->cell_dta[m-1].tb      );
#else
    /* off = (int) (wp->cell_dta[m].tb_ea[1] - wp->cell_dta[m-1].tb_ea[1]);
       off = off / sizeof(short);
     */
    off = (int)((wp->cell_dta[m].tb_ea - wp->cell_dta[m-1].tb_ea) / (sizeof (short)));
    
#endif
    /* tb_off[k*NNBCELL+m] = tb_off[k*NNBCELL+m-1] + off;  */
    tb_off_Nk[m] = tb_off_Nk[m-1] + off;
    len_max = MAX( len_max, off ); 
    if ( EXPECT_TRUE_(m>1) ) { 
       len_max2 = MAX( len_max2, off );
    }
    off2 += off;
  }
  len_max = MAX( len_max, wp->flag - off2 ); 
  last_nbl_len = MAX( last_nbl_len, wp->flag );
}


static void store_tb_args(argbuf_t* parg, unsigned n)
{
     __lwsync();
     for (  ;      n>0u;      ++parg, --n ) {
         store_tb((wp_t const*)parg);
     }
}


/******************************************************************************
*
*  calc_tb_ppu
*
******************************************************************************/

static void calc_tb_ppu(void)
{
  /* (Global) work package */
  wp_t* const wp = cbe_wp_begin;
  int k;

  len_max = 0;
  last_nbl_len = 0;

  /* compute neighbor tables */
  for (k=0; k<ncells; ++k ) {
    make_tb(k, wp);
    calc_tb(wp);
    store_tb(wp);
  }
}


/******************************************************************************
*
*  calc_tb_spu
*
******************************************************************************/

static void calc_tb_spu(void)
{
  /* Needed for SPU mailboxes reading/writing */
  unsigned tmpcnt[N_SPU_THREADS_MAX], mbxval[N_SPU_THREADS_MAX*2];
  unsigned ival;
  unsigned const num_spus2 = num_spus << 1u;


  len_max = 0;
  len_max2 = 0;
  last_nbl_len = 0;



  /* flag 1: calculate neighbor tables */
  /* do_work_spu(DOTB);  */


  value2mboxes(cbe_spucontrolarea, num_spus, DOTB,   tmpcnt);

  do_work_spu_mbuf(make_tb_args, store_tb_args);

  value2mboxes(cbe_spucontrolarea, num_spus, 0u,  tmpcnt);


  /* Get respones from SPUs */
  /* fprintf(stdout, "Getting result locations from SPUs\n"); fflush(stdout); */
  mboxes2tuples(cbe_spucontrolarea, num_spus, 1,  mbxval,  tmpcnt);


  /*
  fprintf(stdout,
         "calc_tb_spu: n_max = %d, len_max = %d, len_max2 = %d\n", 
         n_max, len_max, len_max2);
  fflush(stdout);
  */
}


/******************************************************************************
*
*  make_nblist
*
******************************************************************************/

void make_nblist(void)
{
  static int at_max=0, cl_max=0;
  int c, k, at, i;
  int inc_int   = 128 / sizeof(int)   - 1;
  int inc_short = 128 / sizeof(short) - 1;

  /* update reference positions */
  for (k=0; k<ncells; k++) {
    cell *p = cell_array + cnbrs[k].np;
    for (i=0; i<p->n; i++) {
      NBL_POS(p,i,X) = ORT(p,i,X);
      NBL_POS(p,i,Y) = ORT(p,i,Y);
      NBL_POS(p,i,Z) = ORT(p,i,Z);
      ORT    (p,i,W) = 0.0;
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
#ifdef SPU_INT
    PTR2EA( p->impuls,cell_dta[k].imp_ea   );
#endif
    PTR2EA( NULL,     cell_dta[k].ti_ea    );
    PTR2EA( NULL,     cell_dta[k].tb_ea    );
#endif
  }

  /* count atom numbers */
  at=0; n_max=0;
  for (  k=0; EXPECT_TRUE_(k<ncells);   ++k  ) {
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
  if ( EXPECT_FALSE_(NULL==tb) ) {
    nb_max = 
        (0==last_nbl_len)  ?
            ((int)(nbl_size * estimate_nblist_size())) :
	    ((int)(nbl_size * last_nbl_len));

    nb_max = (nb_max + inc_short) & (~inc_short);
    tb = (short *) malloc_aligned( nb_max * (cl_max+1) * sizeof(short), 0, 0);
  }
  else if (last_nbl_len * sqrt(nbl_size) > nb_max) {
    free(tb);
    nb_max = (int) (nbl_size * last_nbl_len);
    nb_max = (nb_max + inc_short) & (~inc_short);
    tb = (short *) malloc_aligned( nb_max * (cl_max+1) * sizeof(short), 0, 0);
  }
  if (  EXPECT_FALSE_((ti==NULL) || (tb==NULL)) ) 
    error("cannot allocate neighbor table");
  for (  c=0; EXPECT_TRUE_(c<=ncells);   c++ /* What else? :-) */  ) {
      tb_off[c*NNBCELL] = c * nb_max; 
  }

#ifndef NBL_ON_PPU

#ifdef ON_PPU
  calc_tb_ppu();
#else
  calc_tb_spu();
#endif

#else

  /* for all cells */
  last_nbl_len = 0; len_max=0;
  for (  c=0;  EXPECT_TRUE_(c<ncells);   ++c ) {

    int   c1 = cnbrs[c].np, m;
    cell  *p = cell_array + c1;
    int   at_inc = (2*p->n + inc_int) & (~inc_int);
    int   nn = tb_off[c*NNBCELL];

    /* for each neighbor cell */
    for (  m=0; EXPECT_TRUE_(m<NNBCELL);   ++m  ) {

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
      if (  EXPECT_FALSE_(nn - tb_off[c*NNBCELL] > nb_max) ) {
        error("neighbor table full - increase nbl_size");
      }
      len_max = MAX( len_max, nn - tb_off[c*NNBCELL+m] );
    }
    last_nbl_len = MAX( last_nbl_len, nn - tb_off[c*NNBCELL] );
  }

#endif  /* NBL_ON_PPU */

  have_valid_nbl = 1;
  nbl_count++;
}


/******************************************************************************
*
*  make_wp
*
******************************************************************************/

static void make_wp(int k, wp_t *wp)
{
  int m;
  cell_nbrs_t const * const c = cnbrs + k;
  int const inc_int = 128 / sizeof(int) - 1;
  int const at_inc  = (2*(CELLPTR(k))->n + inc_int) & (~inc_int);

  wp->k = k;
  wp->n_max = n_max;      /* allocate for this many atoms     */
  wp->len_max = len_max;  /* allocate for this many neighbors */
  wp->ti_len = at_inc;
  wp->flag = 2;

  for ( m=0;   EXPECT_FALSE_(m<NNBCELL);    ++m ) {
    int const n = c->nq[m];
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
#ifndef ON_PPU
  wp->cell_dta[NNBCELL-1].len = len_max2;  /* this is not optimal... */
#endif
}


/* Create narg (or argbufsze if it is less) wp in buffer starting
   at parg, with serial number starting k
 */
static void make_wp_args(argbuf_t* parg, unsigned n, int k)
{
     for ( ;    EXPECT_TRUE_(n>0u);      ++k, --n, ++parg ) {
         make_wp(k, ((wp_t*)parg));
     }

     /* Sync memory */
     __lwsync();
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

static void store_wp_args(argbuf_t* parg, unsigned narg)
{
    /* Sync memory */
     __lwsync();

     for (  ;     EXPECT_TRUE_(narg>0u);     ++parg, --narg ) {
         tot_pot_energy += ((wp_t const*)parg)->totpot;
         virial         += ((wp_t const*)parg)->virial;
     }
}

/* Dummy function: Do not store anyting at all */
static void store_wp_none(argbuf_t* unused1, unsigned unused2)
{}


/******************************************************************************
*
*  do_forces_spu
*
******************************************************************************/

/* (Slightly generalized) multibuffered version */
void do_work_spu_mbuf(void (* const mkargs)(argbuf_t*, unsigned n, int k0),
                      void (* const stargs)(argbuf_t*, unsigned n)
                     )
 {

    /* Number of elements in "low buffer" (starting at 0) and
       number of elements in "high buffer" (starting at NHALF) */
    unsigned  int nlo[N_SPU_THREADS_MAX], nhi[N_SPU_THREADS_MAX];

    /* Number of SPUs, indices */
    register unsigned int const ispumax=num_spus;
    register unsigned int ispu, iarg;

    /* Half the number of buffers (rounded down) */
    unsigned int const nbuf_half = num_bufs/2;
    /* Buffer sizes/offsets */
    unsigned int const narg_lo = nbuf_half;
    unsigned int const narg_hi = num_bufs-nbuf_half;
    unsigned int const offs_lo = 0u;
    unsigned int const offs_hi = nbuf_half;

    /* Workpackage index (== number of workpackages scheduled) */
    register int k;

    /* Number of WPs done & number of WPs remaining to be scheduled
       nwp is the initial number of WPs
     */
    register unsigned nfrom,nto;

    /* Initialization */
    for ( ispu=0;   EXPECT_TRUE_(ispu<ispumax);   ++ispu ) {
        /* Scheduling: No work has been scheduled yet */
        nlo[ispu]=nhi[ispu]=0u;
    }

    /* While not all WPs have been completed yet */
    for ( k=0, nto=nfrom=ncells;; ) {

        /* Iterate over all SPUs */
        for (ispu=iarg=0; EXPECT_TRUE_(ispu<ispumax); ++ispu, iarg+=N_ARGBUF) {
	    /* Number of slots which must be available in inbox... */
	    enum { IBXCNT=2u };
            /* ...shifted to inbox count position */
  	    enum { IBX2=((unsigned)(IBXCNT<<INMBX_CNT_SHIFT)) };

            /* Pointers to SPU argument buffers / WP buffers */
            register argbuf_t* const parg = cbe_arg_begin + iarg;
            register unsigned int narg;
          
  	    /* Control area for current SPU, used for mailboxing */
   	    spe_spu_control_area_p const pctl  = cbe_spucontrolarea[ispu];


            /* If there is still work to be scheduled, use current SPU
               in case it has a free buffer.
             */
           
  	    /* Use lower buffer? */
   	    if ( EXPECT_TRUE_(nto>0u) && (0u==nlo[ispu]) &&
                 (((pctl->SPU_Mbox_Stat) & INMBX_CNT_MASK) >= IBX2)  ) {
                /* 1st msg. to SPU: number of arguments */
                pctl->SPU_In_Mbox = (narg = (narg_lo<nto ? narg_lo : nto));

                /* Create arguments */
                mkargs(parg+offs_lo, narg, k);

                /* 2nd msg. to SPU: offset */
                pctl->SPU_In_Mbox = offs_lo;


                /* Mark low buffer as being filled with ncrea elements */
                nto -= (nlo[ispu]=narg);
                k   += narg;
            }

            /* Use higher buffer? */
            if ( EXPECT_TRUE_(nto>0u) && (0u==nhi[ispu]) &&
                 (((pctl->SPU_Mbox_Stat) & INMBX_CNT_MASK) >= IBX2)   ) {
                /* 1st msg. to SPU: number of arguments */
                pctl->SPU_In_Mbox  =  (narg = (narg_hi<nto ? narg_hi : nto));

                /* Create arguments */
                mkargs(parg+offs_hi, narg, k);

                /* 2nd msg. to SPU: offset */
                pctl->SPU_In_Mbox = offs_hi;

		/* Mark high buffer as being filled with ncrea elements */
		nto -= (nhi[ispu]=narg);
                k   += narg;
            }

	    /* Can we read from Mbox?
               That is: has SPU finished working on some WP? */
            if ( 0u != ((pctl->SPU_Mbox_Stat) & OUTMBX_CNT_MASK) ) {
		/* Read Mbox message which is just the offset */
                unsigned const offs = pctl->SPU_Out_Mbox;

                /* Number of WPs / offset */
                unsigned int* const pn = ((offs_lo==offs) ? nlo : nhi);

                /* Get number of arguments & store them back */
                stargs(parg+offs,  (narg=pn[ispu]));


                /* Update number of WPs to be fetched, leave loop
                   if we are finished */
                if ( EXPECT_FALSE_(0u == (nfrom-=narg)) ) {
	  	    return;
                }

                /* Buffer is free again */
                pn[ispu]=0u;
            }
           
        }
    }

}


void do_work_spu(int const flag)
{
  /* The following code uses mailboxes to synchronize between PPU & SPUs.
     To do so, it uses direct (memory mapped) access to the control area
     of the SPUs.
     Writing/reading to the members of the control block struct assume
     that writes and reads to unsigneds are atomic.
  */

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

  /* fprintf(stdout, "nallcells=%i  ncells=%i\n", nallcells, ncells); 
     fflush(stdout); */

  /* While there is still work to be scheduled or there are still 
     results to be picked up. */
  for (  k=kdone=0;       (kdone<ncells) || (k<ncells);      )  {

      /* fprintf(stdout, "k=%u  kdone=%u\n", k,kdone); */

      /* Index for wp */
      unsigned iwp1;
      for (  ispu=iwp1=0;     (ispu<ispumax);     ++ispu, iwp1+=N_ARGBUF  ) {

           /* Control area used to access the mailbox and ptr. to state */
  	   spe_spu_control_area_p const pctl  = cbe_spucontrolarea[ispu];
           Tspustate*             const pstat = spustate+ispu;

           wp_t*   const pwp   = ((wp_t*)(cbe_arg_begin+iwp1));

           /* fprintf(stdout, "SPU %u is %s\n", ispu, 
                      ((WORKING==*pstat) ? "working" : "notworking"));  */
 
           /* Idle SPU and still work to schedule */
           if (  (IDLE==*pstat) && (k<ncells)  ) {
   	       /* Create a work package... */
	       switch (flag) {
	         case DOTB: make_tb(k, pwp); break;  /* neighbor tables */
	         case DOWP: make_wp(k, pwp); break;  /* force calc */
	         default: error("unknown flag in do_work_spu");
               }
 
               /* Signal SPU by writing to its mailbox when space 
                  is available */
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

                        /* Store the wp */
   	                switch (flag) {
	                  case DOTB: store_tb(pwp); break; /* neighbor tab */
	                  case DOWP: store_wp(pwp); break; /* force calc */
	                  default: error("unknown flag in do_work_spu");
                        }

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
  wp_t* const wp = cbe_wp_begin;
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
  unsigned int tmpcnt[N_SPU_THREADS_MAX], mbxval[N_SPU_THREADS_MAX*2];
  unsigned int ival, iarg;
  unsigned const num_spus2 = num_spus << 1u;

  int  k, i;
  float tmpvec1[7], tmpvec2[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

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
  if ( 0==have_valid_nbl ) { 
      make_nblist();
  }

  /* clear global accumulation variables */
  tot_pot_energy = virial = tot_kin_energy = 0.0;
#ifdef NVT
  E_kin_2 = 0.0;
#endif
#ifdef FNORM
  fnorm = 0.0;
#endif
#ifdef GLOK
  PxF = pnorm = 0.0;
#endif
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

#ifdef SPU_INT
  /* prepare extra data for integrator on SPU */
  mk_mv();
#endif

  /* compute the forces - either on PPU or on SPUs*/
#ifdef ON_PPU
  calc_forces_ppu(steps); 
#else
  /* flag 2: calculate pair forces */
  /* do_work_spu(DOWP);   */

  /* Tell all SPUs to go to "mode WP" */
  value2mboxes(cbe_spucontrolarea, num_spus,  DOWP,  tmpcnt);

  /* Distribute work */
  do_work_spu_mbuf(make_wp_args, store_wp_none);

  /* Tell SPUs to stop processing WPs */
  value2mboxes(cbe_spucontrolarea, num_spus, 0u,    tmpcnt);

  /* Get responses from SPU mailboxes */
  /* fprintf(stdout, "Getting result locations from SPUs\n"); fflush(stdout);*/
  mboxes2tuples(cbe_spucontrolarea, num_spus, 1,  mbxval,  tmpcnt);


#ifndef PPUSUM

  /* Sum up energies */
  __lwsync();
  for ( ival=iarg=0u;  EXPECT_TRUE_(ival<num_spus);  ++ival, iarg+=N_ARGBUF ) {
      /* Energy/virial WP is located at offset mbxval[ival] 
         in argument buffer of SPU # ival. */
      register wp_t const* const pwp = 
         (wp_t const*)(cbe_arg_begin+iarg+mbxval[ival]);
      tot_pot_energy += pwp->totpot;
      virial         += pwp->virial;
#ifdef SPU_INT
      tot_kin_energy += pwp->totkin;
#ifdef NVT
      E_kin_2        += pwp->E_kin_2;
#endif
#ifdef FNORM
      fnorm          += pwp->fnorm;
#endif
#ifdef GLOK
      PxF            += pwp->PxF;
      pnorm          += pwp->pnorm;
#endif
#endif  /* SPU_INT */
  }

#endif

#endif


#ifdef MPI
  /* sum up results of different CPUs */
  tmpvec1[0]     = tot_pot_energy;
  tmpvec1[1]     = virial;
  tmpvec1[2]     = tot_kin_energy;
  tmpvec1[3]     = E_kin_2;
  tmpvec1[4]     = fnorm;
  tmpvec1[5]     = PxF;
  tmpvec1[6]     = pnorm;

  MPI_Allreduce( tmpvec1, tmpvec2, 7, MPI_FLOAT, MPI_SUM, cpugrid);

  tot_pot_energy = tmpvec2[0];
  virial         = tmpvec2[1];
  tot_kin_energy = tmpvec2[2];
  E_kin_2        = tmpvec2[3];
  fnorm          = tmpvec2[4];
  PxF            = tmpvec2[5];
  pnorm          = tmpvec2[6];
#endif

#ifdef AR
  /* add forces back to original cells/cpus */
  send_forces(add_forces,pack_forces,unpack_forces);
#endif

#ifdef NVT
  /* time evolution of constraints */
  if (ensemble == ENS_NVT) {
    float ttt = nactive * temperature;
    eta += timestep * (E_kin_2 / ttt - 1.0) * isq_tau_eta;
  }
#endif

#ifdef GLOK
  PxF /= (sqrtf(fnorm) * sqrtf(pnorm));
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
