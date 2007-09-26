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
* imd_cbe_calc_spu.c  -- calc_wp on SPU
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

/* ISO C std. headers */
#include <stdio.h>
/* SPU functions/macrso, DMA */
#include <spu_intrinsics.h>
#include <spu_mfcio.h>

/* The IMD configuration (macro CBE_DIRECT) */
#include "config.h"
/* DMA typedefs...*/
#include "imd_cbe.h"




static INLINE_ int iceil128(int const x) {
    return (x+127) & (~127);
}

static INLINE_ int iceil16(int const x) {
    return (x+15) & (~15);
}




#ifdef ON_PPU

/******************************************************************************
*
*  Lennard-Jones interactions on SPU - dummy routine, real one runs on PPU
*
******************************************************************************/

void calc_wp(wp_t *wp) {}

#else  /* not ON_PPU */

#ifdef CBE_DIRECT

/******************************************************************************
*
*  Lennard-Jones interactions on SPU - CBE_DIRECT version
*
******************************************************************************/





/* Helper routine to set all elements in [p,p+k) to x */
static INLINE_ void fvecset(register vector float* p, register int k,
                            register vector float const x)
{
     for (   ; (k>0);    --k, ++p ) {
         (*p) = x;
     }
}


/*  Allocates Pointers inside *buf with at least the following sizes:
      *pos, *force:  wp->n_max * 4 * sizeof(float)
      *typ, *ti:     wp->n_max * sizeof(int), wp->n_max * 2 * sizeof(int)
      *tb:           wp->len_max * sizeof(short)

      If pf0 is non-NULL the force array is filled with *pf0
*/
static void allocate_buf(cell_dta_t* const buf,
                         int const n_max, int const len_max,
                         vector float const* const pf0
                        )
{
    /* Current length of workspace */
    extern unsigned wrkspc_len;

    /* pos & force have the same size */
    int const spos = iceil128(  n_max *     sizeof(flt[4]));
    /* Sizes of the othe members */
    int const styp = iceil128(  n_max *     sizeof(int));
    int const sti  = iceil128(  n_max * 2 * sizeof(int));
    int const stb  = iceil128(len_max *     sizeof(short));

    /* Total number of bytes needed */
    int const stot = spos + spos + styp + sti + stb;


    /* Enough (work) space? */
    if ( EXPECT_TRUE(wrkspc_len>=stot) ) {
        /* Type "pointer to byte" & "raw memory pointer"
           (used for buffer mgmt) */
        typedef unsigned char* Pbyt;
        typedef void*          Praw;

        /* Current address of workspace */
        extern void*  wrkspc_adr;

        /* Set pointers */
        buf->pos   = (Praw)(wrkspc_adr);
        buf->force = (Praw)(((Pbyt)(buf->pos))   + spos);
        buf->typ   = (Praw)(((Pbyt)(buf->force)) + spos);
        buf->ti    = (Praw)(((Pbyt)(buf->typ))   + styp);
        buf->tb    = (Praw)(((Pbyt)(buf->ti))    + sti);

        /* Update */
        wrkspc_adr  = ((Pbyt)(wrkspc_adr)) + stot;
        wrkspc_len -= stot;

        /* Set forces to zero? */
        if ( pf0 ) {
	    fvecset((vector float*)(buf->force), n_max,  *pf0);
	}
    }
    else {
       fprintf(stderr, "Could not allocate %d bytes in allocate_buf (only %u bytes available)\n", stot, wrkspc_len);

       buf->pos   = 0;
       buf->force = 0;
       buf->typ   = 0;
       buf->ti    = 0;
       buf->tb    = 0;
    }
}






static void mdma64_rec(register unsigned char* const p,
                       register unsigned long long const lea, register unsigned const sze,
                       register unsigned const tag, register unsigned const cmd)
{
    /* Maxium transfer size possible */
    enum { maxsze  = (unsigned)(16*1024)          };
    enum { maxstep = (unsigned long long)(maxsze) };

    /* High/Low part of EA */
    register unsigned const hi = (unsigned)(lea>>32u);
    register unsigned const lo = (unsigned)lea;

    /* Only one DMA needed? End recusrion */
    if ( sze<=maxsze ) {
        spu_mfcdma64(p, hi,lo, sze,  tag,cmd);
        return;
    }

    /* Start a DMA */
    spu_mfcdma64(p, hi,lo, maxsze,  tag,cmd);
    /* Tail recursive call... */
    mdma64_rec(p+maxsze,  lea+maxstep,  sze-maxsze,  tag,cmd);
}









/* DMA more than 16k using multiple DMAs with the same tag */
void mdma64_iter(register unsigned char* p,
                 register unsigned long long lea, register unsigned remsze,
                 register unsigned const tag, register unsigned const cmd)
{
    /* Maximum number of bytes per DMA (which is a multiple of 16) */
    enum { maxsze  = (unsigned)(16*1024)        };
    enum { maxstep = (unsigned long long)maxsze };

    for(;;) {
        /* Split EA into low & hi part */
        register unsigned const hi = (unsigned)(lea>>32u);
        register unsigned const lo = (unsigned)lea;

        /* We're done if not more than maxsze bytes are to be xfered */
        if ( remsze<=maxsze ) {
  	   spu_mfcdma64(p, hi,lo,  remsze, tag, cmd);
 	   return;
        }

        /* Xfer one chunk of maximum size  */
        spu_mfcdma64(p,  hi,lo,  maxsze, tag, cmd);

        /* Update addresses and number of bytes */
        remsze -= maxsze;
        p      += maxsze;
        lea    += maxstep;
    }
}



/* DMA more than 16K using multiple DMAs */
INLINE_ void mdma64(void* const p,
                    unsigned const* const ea, unsigned const size,
                    unsigned const tag, unsigned const cmd
                   )
{
    typedef unsigned long long  T64;

    /* Debugging output */
    /*
    fprintf(stdout, "DMAing %u bytes between %p and (0x%x,0x%x) with tag %u\n", 
            size, p,  ea[0], ea[1],  tag);
    fflush(stdout);
    */

    mdma64_iter(p, (((T64)(ea[0]))<<32u)+((T64)(ea[1])), size,  tag, cmd);
}


/* One DMA only (added for syntactical reasons only) */
INLINE_ void dma64(void* const p,
                   unsigned const* const ea, unsigned const size,
                   unsigned const tag, unsigned const cmd
                  )
{
    spu_mfcdma64(p,  ea[0],ea[1],  size, tag,cmd);
}


  


/* Init a DMA get from the EAs specified in *addr to the LS addresses
   specified ton buf.
 */
static INLINE_ 
  void init_fetch(cell_dta_t* const buf, const int ti_len, 
                  cell_ea_t const* const addr, unsigned const tag)
{
    /* Start 4 seperate DMAs (each of which may be larger than 16K as 
       mdma64 is used).
       Note that list DMA is not possible here, as the ptr. members of *buf
       are aligned to 128-byte boundary, but in list DMA, LS addresses are
       automatically rounded up to the next 16-byte boundary.
     */
    mdma64(buf->pos, addr->pos_ea, iceil16(addr->n*sizeof(flt[4])),
           tag, MFC_GET_CMD);

    mdma64(buf->typ, addr->typ_ea, iceil16(addr->n*sizeof(int)),
           tag, MFC_GET_CMD);
  
    mdma64(buf->ti,  addr->ti_ea,  iceil16(ti_len*sizeof(int)),
           tag, MFC_GET_CMD);

    mdma64(buf->tb,  addr->tb_ea,  iceil16(addr->len*sizeof(short)),
           tag, MFC_GET_CMD);


    /* Do we really have to DMA the force?
    mdma64(buf->force, addr->force_ea, iceil16(addr->n*sizeof(flt[4])),
           tag, MFC_GET_CMD);
    */


    /* Also copy n */
    buf->n = addr->n;
}




/* void wait_fetch( cell_dta_t *buf) {} */
static INLINE_ void wait_fetch(unsigned const tag) {
    spu_writech(MFC_WrTagMask, (1u<<tag));
    spu_mfcstat(MFC_TAG_UPDATE_ALL);
}





/* Start a DMA back to main memory without waiting for it */
static INLINE_
  void return_forces(cell_dta_t const* const buf, cell_ea_t const* const addr)
{
    /* Contains some tag which may be used in return_forces */
    extern unsigned const forces_tag;

    /* Debugging output */
    /* fprintf(stdout, "DMAing back forces with tag %u\n", forces_tag); fflush(stdout); */

    mdma64(buf->force, addr->force_ea, iceil16(addr->n*sizeof(flt[4])),
           forces_tag, MFC_PUT_CMD);
}





/* Macro for the calculation of the Lennard-Jones potential */
#define LJ(pot,grad,r2)  {                                               \
   int   const col4 = col*4;                                             \
   float const tmp2 = lj_sig[col4] / r2;                                 \
   float const tmp6 = tmp2 * tmp2 * tmp2;                                \
   pot  = lj_eps[col4] * tmp6 * (tmp6 - 2.0) - lj_shift[col4];           \
   grad = - 12.0 * lj_eps[col4] * tmp6 * (tmp6 - 1.0) / r2;              \
}






void calc_wp(wp_t *wp)
{
  /* Some DMA tag constants for the 3 buf's 
     btagx is the tag for DMAs to buf[x]
   */
  enum { btag0=0u, btag1=1u, btag2=2u };

  /* Tags for the buffer ptrs.  */
  unsigned qtag, next_qtag, old_qtag;
  cell_dta_t buf[3], *p, *q, *next_q, *old_q;

  vector float f00  = spu_splats( (float)   0.0  );
  vector float f001 = spu_splats( (float)   0.001);
#ifdef AR
  vector float f05n = spu_splats( (float)  -0.5  );
  vector float f05l = {0.0, 0.0, 0.0, 0.5};
#else
  vector float f05  = spu_splats( (float)   0.5  );
  vector float f1l  = {0.0, 0.0, 0.0, 1.0};
#endif
  vector float f10  = spu_splats( (float)   1.0  );
  vector float f20  = spu_splats( (float)   2.0  );
  vector float f12  = spu_splats( (float) -12.0  );
  vector float vir  = spu_splats( (float)   0.0  );
  vector float r2cut   = pt.r2cut   [0];
  vector float ljsig   = pt.lj_sig  [0];
  vector float ljeps   = pt.lj_eps  [0];
  vector float ljshift = pt.lj_shift[0];
  vector signed   int  i00 = spu_splats( (int) 0 );
#ifdef AR
  vector unsigned char s0  = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19};
  vector unsigned char s1  = {0,1,2,3,4,5,6,7,8,9,10,11,20,21,22,23};
  vector unsigned char s2  = {0,1,2,3,4,5,6,7,8,9,10,11,24,25,26,27};
  vector unsigned char s3  = {0,1,2,3,4,5,6,7,8,9,10,11,28,29,30,31};
#endif
  vector unsigned char ss0 = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
  vector unsigned char ss1 = {4,5,6,7,4,5,6,7,4,5,6,7,4,5,6,7};
  vector unsigned char ss2 = {8,9,10,11,8,9,10,11,8,9,10,11,8,9,10,11};
  vector unsigned char ss3 = {12,13,14,15,12,13,14,15,12,13,14,15,12,13,14,15};
  vector float d0, d1, d2, d3, d20, d21, d22, d23, r2, r2i;
  vector float tmp2, tmp3, tmp4, tmp6, tmp7, pot;
  vector float grad, pots, ff0, ff1, ff2, ff3;
  vector float d2a, d2b, ffa, ffb, dummy=f00, ff0s, ff1s, ff2s, ff3s;
  vector float *pi, *fi, *qpos;
#ifdef AR
  vector float *qforce;
#endif
  vector unsigned int ms, ms1, ms2, ms3;
  vector signed int tj;
  int    i, c1, n, m, j0, j1, j2, j3;

  /* fprintf(stdout, "calc_wp (DIRECT) starts on SPU\n");  fflush(stdout); */

  /* The initial fetch to buffer q */
  allocate_buf(&buf[0], wp->n_max, wp->len_max, &f00);
  q = p = &(buf[0]);
  qtag  = btag0;
  init_fetch(q, wp->ti_len, wp->cell_dta, qtag);

  /* Allocate the remaining buffers... */
  allocate_buf(&buf[1], wp->n_max, wp->len_max, NULL);
  allocate_buf(&buf[2], wp->n_max, wp->len_max, NULL);

  /* ...and set pointers to them */
  next_q      = &(buf[1]);
  next_qtag   = btag1;

  /* Set scalars to zero */
  wp->totpot = 0.0;
  wp->virial = 0.0;

  m = 0;
  do {

    wait_fetch(qtag);

    do {
      m++;
    } while ((wp->cell_dta[m].n==0) && (m<NNBCELL));

    if (m<NNBCELL) {
      init_fetch( next_q, wp->ti_len, wp->cell_dta + m, next_qtag );
    }

    /* positions in q as vector float */
    qpos = (vector float *) q->pos;

    for (i=0; i<p->n; i++) {

      short const * const ttb = q->tb + q->ti[2*i];

#ifndef MONO
      c1 = 2 * p->typ[i];
#endif

      fi = ((vector float *) p->force) + i;
      pi = ((vector float *) p->pos)   + i;

      /* clear accumulation variables */
      pots = f00;  ff0s = f00;  ff1s = f00;  ff2s = f00;  ff3s = f00;

      /* we treat four neighbors at a time, so that we can vectorize */
      /* ttb is padded with copies of i, which have to be masked     */
      qpos[q->n] = *pi;      /* needed for the padding values in ttb */
      for (n=0; n<q->ti[2*i+1]; n+=4) {

        /* indices of neighbors */
        j0  = ttb[n  ];
        j1  = ttb[n+1];
        j2  = ttb[n+2];
        j3  = ttb[n+3];

#ifndef MONO
        /* if not MONO, we assume up to two atom types */
        /* mask for type dependent selections */
        tj = spu_promote( q->typ[j0],     0 );
        tj = spu_insert ( q->typ[j1], tj, 1 );
        tj = spu_insert ( q->typ[j2], tj, 2 );
        tj = spu_insert ( q->typ[j3], tj, 3 );
        ms = spu_cmpeq( tj, i00 );
        r2cut   = spu_sel( pt.r2cut   [c1+1], pt.r2cut   [c1], ms );
        ljsig   = spu_sel( pt.lj_sig  [c1+1], pt.lj_sig  [c1], ms );
        ljeps   = spu_sel( pt.lj_eps  [c1+1], pt.lj_eps  [c1], ms );
        ljshift = spu_sel( pt.lj_shift[c1+1], pt.lj_shift[c1], ms );
#endif

        /* distance vectors */
        d0  = spu_sub( qpos[j0], *pi );
        d1  = spu_sub( qpos[j1], *pi );
        d2  = spu_sub( qpos[j2], *pi );
        d3  = spu_sub( qpos[j3], *pi );

        /* clear the 4th component */
        d0 = spu_sel( d0, f00, spu_maskw(1) );
        d1 = spu_sel( d1, f00, spu_maskw(1) );
        d2 = spu_sel( d2, f00, spu_maskw(1) );
        d3 = spu_sel( d3, f00, spu_maskw(1) );

        d20 = spu_mul( d0, d0 );
        d21 = spu_mul( d1, d1 );
        d22 = spu_mul( d2, d2 );
        d23 = spu_mul( d3, d3 );

        d20 = spu_add( d20, spu_rlqwbyte( d20, 8 ) );
        d21 = spu_add( d21, spu_rlqwbyte( d21, 8 ) );
        d22 = spu_add( d22, spu_rlqwbyte( d22, 8 ) );
        d23 = spu_add( d23, spu_rlqwbyte( d23, 8 ) );

        d20 = spu_add( d20, spu_rlqwbyte( d20, 4 ) );
        d21 = spu_add( d21, spu_rlqwbyte( d21, 4 ) );
        d22 = spu_add( d22, spu_rlqwbyte( d22, 4 ) );
        d23 = spu_add( d23, spu_rlqwbyte( d23, 4 ) );

        d2a = spu_sel( d20, d21, spu_maskw(7) );
        d2b = spu_sel( d22, d23, spu_maskw(1) );
        r2  = spu_sel( d2a, d2b, spu_maskw(3) );

        /* compute inverses of r2 */
        ms1 = spu_cmpgt( r2cut, r2 );  /* cutoff mask */
        ms2 = spu_cmpgt( r2, f001 );   /* mask zeros */
        r2  = spu_sel( f10, r2, ms2  );
        /* first estimate inverse, then sharpen it */
        r2i = spu_re( r2 ); 
        r2i = spu_mul( r2i, spu_sub( f20, spu_mul( r2i, r2 ) ) ); 

        /* mask unwanted values */
        ms3 = spu_and( ms1, ms2 );
        r2  = spu_sel( f00, r2,  ms3 );
        r2i = spu_sel( f00, r2i, ms3 );

        /* compute LJ interaction */
        tmp2 = spu_mul( r2i,  ljsig );
        tmp3 = spu_mul( tmp2, ljeps );
        tmp4 = spu_mul( tmp2, tmp2 );
        tmp6 = spu_mul( tmp4, tmp2 );
        tmp7 = spu_mul( tmp4, tmp3 );
        pot  = spu_msub( tmp7, spu_sub(tmp6, f20), spu_sel(f00, ljshift, ms3));
        grad = spu_mul( spu_msub(tmp7, tmp6, tmp7), spu_mul(f12, r2i) );

#ifdef AR
        /* add up potential energy */
        pots = spu_add( pots, pot );
        pot  = spu_mul( pot,  f05n );  /* avoid double counting */
#else
        /* add up potential energy */
        pots = spu_madd( pot, f05, pots );  /* avoid double counting */
#endif
        /* add to total virial */
        vir  = spu_madd( r2, grad, vir );

        /* the forces */
        ff0  = spu_mul( d0, spu_shuffle( grad, dummy, ss0 ) );
        ff1  = spu_mul( d1, spu_shuffle( grad, dummy, ss1 ) );
        ff2  = spu_mul( d2, spu_shuffle( grad, dummy, ss2 ) );
        ff3  = spu_mul( d3, spu_shuffle( grad, dummy, ss3 ) );

        /* add forces and potential on first particle */
        ff0s = spu_add( ff0s, ff0 );
        ff1s = spu_add( ff1s, ff1 );
        ff2s = spu_add( ff2s, ff2 );
        ff3s = spu_add( ff3s, ff3 );

#ifdef AR
        /* add forces and potential on second particle */
        qforce = (vector float *) q->force;
        qforce[j0] = spu_sub( qforce[j0], spu_shuffle( ff0, pot, s0 ) ); 
        qforce[j1] = spu_sub( qforce[j1], spu_shuffle( ff1, pot, s1 ) ); 
        qforce[j2] = spu_sub( qforce[j2], spu_shuffle( ff2, pot, s2 ) ); 
        qforce[j3] = spu_sub( qforce[j3], spu_shuffle( ff3, pot, s3 ) );
#endif
      }

      /* add contribution to total poteng */
      pots = spu_add( pots, spu_rlqwbyte( pots, 8 ) );
      pots = spu_add( pots, spu_rlqwbyte( pots, 4 ) );
      wp->totpot += spu_extract( pots, 0 );

      /* add force of first particle */
      ffa = spu_add( ff0s, ff1s );
      ffb = spu_add( ff2s, ff3s );
      *fi = spu_add( *fi,  ffa  );
      *fi = spu_add( *fi,  ffb  );

      /* add potential of first particle */
#ifdef AR
      *fi = spu_madd( pots, f05l, *fi );
#else
      *fi = spu_madd( pots, f1l, *fi );
#endif
    }

#ifdef AR
    /* write back forces in *q - locking required! */
#endif

    /*
    if ((0==wp->k) && (NNBCELL==m)) {
      for (i=0; i<p->n; i++)
        printf("%d %d %e %e %e %e %e %e %e\n", 
          i, p->typ[i], p->pos[4*i], p->pos[4*i+1], p->pos[4*i+2],
          p->force[4*i], p->force[4*i+1], p->force[4*i+2], p->force[4*i+3]);
    }
    */

    if (q == p) {  
       old_q    = &(buf[2]);
       old_qtag = btag2;
    }
    else {
       old_q    = q;
       old_qtag = qtag;
    }
 
    q    = next_q;
    qtag = next_qtag;

    next_q    = old_q;
    next_qtag = old_qtag;

  } while (m<NNBCELL);

  /* set contribution to total virial */
#ifndef AR
  vir = spu_mul( vir, f05 );  /* avoid double counting */
#endif
  vir = spu_add( vir, spu_rlqwbyte( vir, 8 ) );
  vir = spu_add( vir, spu_rlqwbyte( vir, 4 ) );
  wp->virial = - spu_extract( vir, 0 );

  /* Init DMA back to main memory. No need to wait for it to complete here,
     as we will do that in the main (work) loop. */
  return_forces(p, wp->cell_dta);
}

#else  /* not CBE_DIRECT */

/******************************************************************************
*
*  Lennard-Jones interactions on SPU - indirect version
*
******************************************************************************/

void calc_wp(wp_t *wp)
{
  vector float f00  = spu_splats( (float)   0.0  );
  vector float f001 = spu_splats( (float)   0.001);
  vector float f05n = spu_splats( (float)  -0.5  );
  vector float f05l = {0.0, 0.0, 0.0, 0.5};
  vector float f10  = spu_splats( (float)   1.0  );
  vector float f20  = spu_splats( (float)   2.0  );
  vector float f12  = spu_splats( (float) -12.0  );
  vector float vir  = spu_splats( (float)   0.0  );
  vector float r2cut   = pt.r2cut   [0];
  vector float ljsig   = pt.lj_sig  [0];
  vector float ljeps   = pt.lj_eps  [0];
  vector float ljshift = pt.lj_shift[0];
  vector signed   int  i00 = spu_splats( (int) 0 );
  vector unsigned char s0  = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19};
  vector unsigned char s1  = {0,1,2,3,4,5,6,7,8,9,10,11,20,21,22,23};
  vector unsigned char s2  = {0,1,2,3,4,5,6,7,8,9,10,11,24,25,26,27};
  vector unsigned char s3  = {0,1,2,3,4,5,6,7,8,9,10,11,28,29,30,31};
  vector unsigned char ss0 = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
  vector unsigned char ss1 = {4,5,6,7,4,5,6,7,4,5,6,7,4,5,6,7};
  vector unsigned char ss2 = {8,9,10,11,8,9,10,11,8,9,10,11,8,9,10,11};
  vector unsigned char ss3 = {12,13,14,15,12,13,14,15,12,13,14,15,12,13,14,15};
  vector float d0, d1, d2, d3, d20, d21, d22, d23, r2, r2i;
  vector float tmp2, tmp3, tmp4, tmp6, tmp7, pot;
  vector float grad, pots, ff0, ff1, ff2, ff3;
  vector float d2a, d2b, ffa, ffb, dummy=f00, *fi, ff0s, ff1s, ff2s, ff3s;
  vector float const* pi;
  vector unsigned int ms, ms1, ms2, ms3;
  vector signed int tj;
  int    i, c1, n, j0, j1, j2, j3;

  /* clear accumulation variables */
  wp->totpot = 0.0;
  wp->virial = 0.0;
  for (i=0; i<wp->n2; i++) wp->force[i] = f00;

  /* fetch sizes and pointers
     fetch data at wp->pos, wp->typ, wp->ti */
  for (i=0; i<wp->n1; i++) {

    short const *ttb;

    /* fetch data at wp->tb + wp->ti[2*i] */
    ttb = wp->tb + wp->ti[2*i];

#ifndef MONO
    c1 = 2 * wp->typ[i];
#endif

    fi = wp->force + i;
    pi = wp->pos   + i;

    /* clear accumulation variables */
    pots = f00;  ff0s = f00;  ff1s = f00;  ff2s = f00;  ff3s = f00;

    /* we treat four neighbors at a time, so that we can vectorize */
    /* ttb is padded with copies of i, which have to be masked     */
    for (n = 0; n < wp->ti[2*i+1]; n += 4) {

      /* indices of neighbors */
      j0  = ttb[n  ];
      j1  = ttb[n+1];
      j2  = ttb[n+2];
      j3  = ttb[n+3];

#ifndef MONO
      /* if not MONO, we assume up to two atom types */
      /* mask for type dependent selections */
      tj = spu_promote( wp->typ[j0],     0 );
      tj = spu_insert ( wp->typ[j1], tj, 1 );
      tj = spu_insert ( wp->typ[j2], tj, 2 );
      tj = spu_insert ( wp->typ[j3], tj, 3 );
      ms = spu_cmpeq( tj, i00 );
      r2cut   = spu_sel( pt.r2cut   [c1+1], pt.r2cut   [c1], ms );
      ljsig   = spu_sel( pt.lj_sig  [c1+1], pt.lj_sig  [c1], ms );
      ljeps   = spu_sel( pt.lj_eps  [c1+1], pt.lj_eps  [c1], ms );
      ljshift = spu_sel( pt.lj_shift[c1+1], pt.lj_shift[c1], ms );
#endif

      /* distance vectors */
      d0  = spu_sub( wp->pos[j0], *pi );
      d1  = spu_sub( wp->pos[j1], *pi );
      d2  = spu_sub( wp->pos[j2], *pi );
      d3  = spu_sub( wp->pos[j3], *pi );

      d20 = spu_mul( d0, d0 );
      d21 = spu_mul( d1, d1 );
      d22 = spu_mul( d2, d2 );
      d23 = spu_mul( d3, d3 );

      d20 = spu_add( d20, spu_rlqwbyte( d20, 8 ) );
      d21 = spu_add( d21, spu_rlqwbyte( d21, 8 ) );
      d22 = spu_add( d22, spu_rlqwbyte( d22, 8 ) );
      d23 = spu_add( d23, spu_rlqwbyte( d23, 8 ) );

      d20 = spu_add( d20, spu_rlqwbyte( d20, 4 ) );
      d21 = spu_add( d21, spu_rlqwbyte( d21, 4 ) );
      d22 = spu_add( d22, spu_rlqwbyte( d22, 4 ) );
      d23 = spu_add( d23, spu_rlqwbyte( d23, 4 ) );

      d2a = spu_sel( d20, d21, spu_maskw(7) );
      d2b = spu_sel( d22, d23, spu_maskw(1) );
      r2  = spu_sel( d2a, d2b, spu_maskw(3) );

      /* compute inverses of r2 */
      ms1 = spu_cmpgt( r2cut, r2 );  /* cutoff mask */
      ms2 = spu_cmpgt( r2, f001 );   /* mask zeros */
      r2  = spu_sel( f10, r2, ms2  );
      /* first estimate inverse, then sharpen it */
      r2i = spu_re( r2 ); 
      r2i = spu_mul( r2i, spu_sub( f20, spu_mul( r2i, r2 ) ) ); 

      /* mask unwanted values */
      ms3 = spu_and( ms1, ms2 );
      r2  = spu_sel( f00, r2,  ms3 );
      r2i = spu_sel( f00, r2i, ms3 );

      /* compute LJ interaction */
      tmp2 = spu_mul( r2i,  ljsig );
      tmp3 = spu_mul( tmp2, ljeps );
      tmp4 = spu_mul( tmp2, tmp2 );
      tmp6 = spu_mul( tmp4, tmp2 );
      tmp7 = spu_mul( tmp4, tmp3 );
      pot  = spu_msub( tmp7, spu_sub(tmp6, f20), spu_sel(f00, ljshift, ms3) );
      grad = spu_mul( spu_msub(tmp7, tmp6, tmp7), spu_mul(f12, r2i) );

      /* add up potential energy */
      pots = spu_add( pots, pot );
      pot  = spu_mul( pot,  f05n );  /* avoid double counting */

      /* add to total virial */
      vir  = spu_madd( r2, grad, vir );

      /* the forces */
      ff0  = spu_mul( d0, spu_shuffle( grad, dummy, ss0 ) );
      ff1  = spu_mul( d1, spu_shuffle( grad, dummy, ss1 ) );
      ff2  = spu_mul( d2, spu_shuffle( grad, dummy, ss2 ) );
      ff3  = spu_mul( d3, spu_shuffle( grad, dummy, ss3 ) );

      /* add forces and potential on first particle */
      ff0s = spu_add( ff0s, ff0 );
      ff1s = spu_add( ff1s, ff1 );
      ff2s = spu_add( ff2s, ff2 );
      ff3s = spu_add( ff3s, ff3 );

      /* add forces and potential on second particle */
      wp->force[j0] = spu_sub( wp->force[j0], spu_shuffle( ff0, pot, s0 ) ); 
      wp->force[j1] = spu_sub( wp->force[j1], spu_shuffle( ff1, pot, s1 ) ); 
      wp->force[j2] = spu_sub( wp->force[j2], spu_shuffle( ff2, pot, s2 ) ); 
      wp->force[j3] = spu_sub( wp->force[j3], spu_shuffle( ff3, pot, s3 ) );

    }

    /* add contribution to total poteng */
    pots = spu_add( pots, spu_rlqwbyte( pots, 8 ) );
    pots = spu_add( pots, spu_rlqwbyte( pots, 4 ) );
    wp->totpot += spu_extract( pots, 0 );

    /* add force of first particle */
    ffa = spu_add( ff0s, ff1s );
    ffb = spu_add( ff2s, ff3s );
    *fi = spu_add( *fi,  ffa  );
    *fi = spu_add( *fi,  ffb  );

    /* add potential of first particle */
    *fi = spu_madd( pots, f05l, *fi );

  } 

  /* set contribution to total virial */
  vir = spu_add( vir, spu_rlqwbyte( vir, 8 ) );
  vir = spu_add( vir, spu_rlqwbyte( vir, 4 ) );
  wp->virial = - spu_extract( vir, 0 );

  /* return forces, poteng, totpot, virial */

}

#endif  /* not CBE_DIRECT */

#endif  /* not ON_PPU */
