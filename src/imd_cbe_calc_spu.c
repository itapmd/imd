
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
#ifdef SPU_INT
#include <vector/sum_across_float3.h>
#endif
#include "config.h"
#include "imd_cbe.h"


#ifdef ON_PPU

/******************************************************************************
*
*  Pair interactions on SPU - dummy routine, real one runs on PPU
*
******************************************************************************/

void calc_wp_direct(wp_t *wp, cell_dta_t* const buf, unsigned const otag) {}

/******************************************************************************
*
*  Neighbor tables on SPU - dummy routine, real one runs on PPU
*
******************************************************************************/

void calc_tb_direct(wp_t *wp, cell_dta_t* const buf, unsigned const otag) {}

#else  /* not ON_PPU */


/* Allocate memory for 3 buffers pointed to by pbuf:
   Allocates pointers inside *pbuf with at least the following sizes:
      *pos, *force:  wp->n_max * 4 * sizeof(float)
      *typ, *ti:     wp->n_max * sizeof(int), wp->n_max * 2 * sizeof(int)
      *tb:           wp->len_max * sizeof(short)

   Only pbuf[0]->force is allocated (and its components set to zero)
   pbuf[1]->force & pbuf[2]->force are set to NULL as they won't be needed

 */


/* Init a DMA get from the EAs specified in *addr to the LS addresses
   specified ton buf.
 */
static INLINE_ 
  void init_fetch_wp(cell_dta_t* const buf, const int ti_len, 
                  cell_ea_t const* const addr, unsigned const tag, 
                  int const flag)
{
    /* Start 4 seperate DMAs.
       Note that list DMA is not possible here, as the ptr. members of *buf
       are aligned to 128-byte boundary, but in list DMA, LS addresses are
       automatically rounded up to the next 16-byte boundary.
     */
    dma64(buf->pos, addr->pos_ea, iceil16(addr->n*sizeof(vector float)),
           tag, MFC_GET_CMD);

#ifndef MONO
    dma64(buf->typ, addr->typ_ea, iceil16(addr->n*sizeof(int)),
           tag, MFC_GET_CMD);
#endif

    dma64(buf->ti,  addr->ti_ea,  iceil16(ti_len*sizeof(int)),
           tag, MFC_GET_CMD);

    /* here we possibly need a multi-DMA */
    mdma64(buf->tb,  addr->tb_ea,  iceil16(addr->len*sizeof(short)),
           tag, MFC_GET_CMD);

#ifdef SPU_INT
    if ( EXPECT_FALSE_(flag) )
      dma64(buf->imp, addr->imp_ea, iceil16(addr->n*sizeof(vector float)),
            tag, MFC_GET_CMD);
#endif

    /* Also copy n */
    buf->n = addr->n;
}



static INLINE_ 
  void init_fetch_tb(cell_dta_t* const buf, 
                     cell_ea_t const* const addr, unsigned const tag)
{
    /* we need only the positions here */
    dma64(buf->pos, addr->pos_ea, iceil16(addr->n*sizeof(vector float)),
           tag, MFC_GET_CMD);

    /* Also copy n */
    buf->n = addr->n;
}



static void wait_tag(unsigned const tag) {
    wait_dma((1u<<tag),  MFC_TAG_UPDATE_ALL);
}



/* Start a DMA back to main memory without waiting for it */
static INLINE_
  void init_return(cell_dta_t const* const buf, cell_ea_t const* const addr,
                   unsigned const tag)
{
    dma64(buf->force, addr->force_ea, iceil16(addr->n*sizeof(vector float)),
           tag, MFC_PUT_CMD);
#ifdef SPU_INT
    dma64(buf->pos,   addr->pos_ea,   iceil16(addr->n*sizeof(vector float)),
           tag, MFC_PUT_CMD);
    dma64(buf->imp,   addr->imp_ea,   iceil16(addr->n*sizeof(vector float)),
           tag, MFC_PUT_CMD);
#endif
}


/* Start a DMA back to main memory without waiting for it */
static INLINE_
  void init_return_tb(cell_dta_t const* const buf, cell_ea_t const* const addr,
                      int const ti_len, int const tb_len, unsigned const tag)
{
  dma64 (buf->ti, addr->ti_ea, iceil16(ti_len*sizeof(int)  ), tag,MFC_PUT_CMD);
  /* here we possibly need a multi-DMA */
  mdma64(buf->tb, addr->tb_ea, iceil16(tb_len*sizeof(short)), tag,MFC_PUT_CMD);
}


/******************************************************************************
*
*  Integrator for SPU -- can do NVE, NVT, MIK, GLOK
*
******************************************************************************/

#ifdef SPU_INT

void vecprint(char *str, vector float v)
{
  printf("%s %f %f %f %f\n", str, spu_extract(v,0), spu_extract(v,1),
         spu_extract(v,2), spu_extract(v,3) );
  fflush(stdout);
}

void move_atoms_spu(wp_t *wp, cell_dta_t *p)
{
  int i;
  vector float tmp;
  vector float Ekin1 = spu_splats( (float)0.0 );
  vector float Ekin2 = spu_splats( (float)0.0 );
#ifdef FNORM
  vector float fnorm = spu_splats( (float)0.0 );
#endif
#ifdef GLOK
  vector float PxF   = spu_splats( (float)0.0 );
  vector float pnorm = spu_splats( (float)0.0 );
#endif

  for (i=0; EXPECT_TRUE_(i<p->n); ++i) {

    vector float * const fce = ((vector float *) p->force) + i;
    vector float * const pos = ((vector float *) p->pos)   + i;
    vector float * const imp = ((vector float *) p->imp)   + i;
    int sort = p->typ[i];

#ifdef FBC
    *fce = spu_add( *fce, mv.fbc_f[sort] ); /* external fbc force */
#endif
    *fce = spu_mul( *fce, mv.restr[sort] ); /* restrictions */

#ifdef FNORM
    fnorm = spu_madd( *fce, *fce, fnorm);   /* add to force norm */
#endif

    /* 2 x old kinetic energy */
    Ekin1 = spu_madd( spu_mul( *imp, *imp ), mv.imass[sort], Ekin1 );
#ifdef NVT
    tmp   = spu_mul( *imp, mv.fac1 );
    *imp  = spu_madd( *fce, mv.fac2, tmp );  /* updated momentum */
#else
    *imp  = spu_madd( *fce, mv.ts, *imp );   /* updated momentum */
#endif

#ifdef MIK
    tmp = spu_splats( _dot_product3( *imp, *fce ) );
    /* for MIK, mikmask is 0.0, else +infinity */
    msk  = spu_cmpgt( mv.mikmask, tmp );
    *imp = spu_sel( *imp, spu_splats( (float)0.0 ), msk );
#endif

#ifdef GLOK
    tmp   = spu_mul( *imp, mv.imass[sort] );
    PxF   = spu_madd( tmp, *fce, PxF );
    pnorm = spu_madd( tmp, tmp, pnorm );
#endif

    /* 2 x new kinetic energy */
    Ekin2 = spu_madd( spu_mul( *imp, *imp ), mv.imass[sort], Ekin2 );

    tmp  = spu_mul( mv.ts, mv.imass[sort] );
    *pos = spu_madd( tmp, *imp, *pos );          /* new position */

  } 

  /* set accumulation variables */
  Ekin2 = spu_mul( Ekin2, spu_splats((float)0.5) );
  Ekin1 = spu_madd( Ekin1, spu_splats((float)0.5), Ekin2 );
  Ekin1 = spu_mul( Ekin1, spu_splats((float)0.5) );

#ifdef PPUSUM

  wp->totkin    = _sum_across_float3( Ekin1 );
#ifdef NVT
  wp->E_kin_2   = _sum_across_float3( Ekin2 );  /* this is needed for NVT */
#endif
#ifdef FNORM
  wp->fnorm     = _sum_across_float3( fnorm );
#endif
#ifdef GLOK
  wp->PxF       = _sum_across_float3( PxF   );
  wp->pnorm     = _sum_across_float3( pnorm );
#endif

#else  /* not PPUSUM */

  psum.totkin  += _sum_across_float3( Ekin1 );
#ifdef NVT
  psum.E_kin_2 += _sum_across_float3( Ekin2 );  /* this is needed for NVT */
#endif
#ifdef FNORM
  psum.fnorm   += _sum_across_float3( fnorm );
#endif
#ifdef GLOK
  psum.PxF     += _sum_across_float3( PxF   );
  psum.pnorm   += _sum_across_float3( pnorm );
#endif

#endif  /* PPUSUM */

}

#endif


/******************************************************************************
*
*  Pair interactions -- Lennard-Jones or tabulated
*
******************************************************************************/

/* This routine will use tags 0,1,2 for DMA as well as otag */
void calc_wp_direct(wp_t *wp,
                    /*
                    void* const tempbuf, unsigned const tempbuf_len,
                    void* const resbuf,  unsigned const resbuf_len,
                    */
                    cell_dta_t* const buf,
                    unsigned const otag)
{
  unsigned tag=0u;  /* tag for incoming data  */

  cell_dta_t  *p, *q, *next_q, *old_q;


  vector float const f00  = spu_splats( (float)   0.0  );
  vector float const f001 = spu_splats( (float)   0.001);
#ifdef AR
  vector float const f05n = spu_splats( (float)  -0.5  );
  vector float const f05l = {0.0, 0.0, 0.0, 0.5};
#else
  vector float const f05  = spu_splats( (float)   0.5  );
  vector float const f1l  = {0.0, 0.0, 0.0, 1.0};
#endif
  vector float const f10  = spu_splats( (float)   1.0  );
  vector float const f20  = spu_splats( (float)   2.0  );
#ifdef LJ
  vector float const f12  = spu_splats( (float) -12.0  );
#else
  vector float const f30  = spu_splats( (float)   3.0  );
  vector float const f6i  = spu_splats( (float) (1.0/6.0) );
#endif
  vector float vir  = f00;
  vector float pott = f00;
#ifdef LJ
  vector float r2cut   = pt.r2cut   [0];
  vector float ljsig   = pt.lj_sig  [0];
  vector float ljeps   = pt.lj_eps  [0];
  vector float ljshift = pt.lj_shift[0];
#else
  vector float r2cut   = pt.r2cut   [0];
  vector float ptbeg   = pt.begin   [0];
  vector float ptstep  = pt.step    [0];
  vector float ptistep = pt.invstep [0];
  vector unsigned int ptcol = pt.tab[0];
#endif
  vector signed int i00   = spu_splats( (int) 0 );
#ifdef AR
  vector unsigned char const s0  = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19};
  vector unsigned char const s1  = {0,1,2,3,4,5,6,7,8,9,10,11,20,21,22,23};
  vector unsigned char const s2  = {0,1,2,3,4,5,6,7,8,9,10,11,24,25,26,27};
  vector unsigned char const s3  = {0,1,2,3,4,5,6,7,8,9,10,11,28,29,30,31};
#endif
  vector unsigned char const ss0 = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
  vector unsigned char const ss1 = {4,5,6,7,4,5,6,7,4,5,6,7,4,5,6,7};
  vector unsigned char const ss2 = {8,9,10,11,8,9,10,11,8,9,10,11,8,9,10,11};
  vector unsigned char const ss3 = {12,13,14,15,12,13,14,15,12,13,14,15,12,13,14,15};
#ifndef LJ
  vector unsigned char const vpat1 = { 0x00, 0x01, 0x02, 0x03, 
                                       0x04, 0x05, 0x06, 0x07,
                                       0x10, 0x11, 0x12, 0x13,  
                                       0x14, 0x15, 0x16, 0x17 };
  vector unsigned char const vpat2 = { 0x08, 0x09, 0x0a, 0x0b,
                                       0x0c, 0x0d, 0x0e, 0x0f,
                                       0x18, 0x19, 0x1a, 0x1b,
                                       0x1c, 0x1d, 0x1e, 0x1f };
  vector unsigned char const vpat3 = { 0x00, 0x01, 0x02, 0x03,
                                       0x10, 0x11, 0x12, 0x13,
                                       0x08, 0x09, 0x0a, 0x0b,
                                       0x18, 0x19, 0x1a, 0x1b };
  vector unsigned char const vpat4 = { 0x04, 0x05, 0x06, 0x07,
                                       0x14, 0x15, 0x16, 0x17,
                                       0x0c, 0x0d, 0x0e, 0x0f,
                                       0x1c, 0x1d, 0x1e, 0x1f };
#endif
  vector float d0, d1, d2, d3, d20, d21, d22, d23, r2, tmp2, tmp3, tmp4;
  vector float pot, grad, pots, ff0, ff1, ff2, ff3;
  vector float d2a, d2b, ffa, ffb, dummy=f00, ff0s, ff1s, ff2s, ff3s;
#ifdef LJ
  vector float r2i, tmp6, tmp7;
#else
  vector float tmp0, tmp1, tmp5, a, a2, a3, b, b2, b3;
  vector float pt0, pt1, pt2, pt3, p1, p2, y1, y2, st6;
#endif
  vector float *qpos;
#ifdef AR
  vector float *qforce;
#endif

  vector unsigned int ms1, ms2;
#ifndef MONO
  vector unsigned int ms;
  vector signed   int tj;
#endif
#ifdef LJ
  vector unsigned int ms3;
#else
  vector unsigned int kk;
#endif
  int    i, c1, n, m;




  /* The initial fetch to buffer q */
  p = q = buf;
  init_fetch_wp(q, wp->ti_len, wp->cell_dta, tag, 1);

  /* Set pointers to the remaining buffers */
  next_q = buf+1;

  /* Set forces to zero using qpos to iterate over elements */
  for ( qpos=(vector float*)(buf[0].force), i=wp->n_max;    EXPECT_TRUE_(i>0);   --i, ++qpos ) {
    (*qpos) = f00;
  }
#ifdef OLD
  /* Set scalars to zero */
  wp->totpot = wp->virial = 0;
#endif
  m = 0;
  do {

    wait_tag(tag);
    do {
        ++m;
    } while ((wp->cell_dta[m].n==0) && (m<NNBCELL));

    if (m<NNBCELL) {
        init_fetch_wp(next_q, wp->ti_len, wp->cell_dta + m, tag, 0);
    }

    /* positions in q as vector float */
    qpos = (vector float *) q->pos;

    for ( i=0;     EXPECT_TRUE_(i<p->n);    ++i ) {

      vector float*       const fi = ((vector float *) p->force) + i;
      vector float const* const pi = ((vector float *) p->pos)   + i;

      short const * const ttb = q->tb + q->ti[2*i];

#ifndef MONO
      c1 = 2 * p->typ[i];
#endif

      /* clear accumulation variables */
      pots = f00;  ff0s = f00;  ff1s = f00;  ff2s = f00;  ff3s = f00;

      /* we treat four neighbors at a time, so that we can vectorize */
      /* ttb is padded with copies of i, which have to be masked     */
      qpos[q->n] = *pi;      /* needed for the padding values in ttb */
      for ( n=0;   EXPECT_TRUE_(n<q->ti[2*i+1]);   n+=4 ) {

        /* indices of neighbors */
        register short int const* const ttbn = ttb+n;
        register int const j0 = ttbn[0];
        register int const j1 = ttbn[1];
        register int const j2 = ttbn[2];
        register int const j3 = ttbn[3];

#ifndef MONO
        /* if not MONO, we assume up to two atom types */
        /* mask for type dependent selections */
        tj = spu_promote( q->typ[j0],     0 );
        tj = spu_insert ( q->typ[j1], tj, 1 );
        tj = spu_insert ( q->typ[j2], tj, 2 );
        tj = spu_insert ( q->typ[j3], tj, 3 );
        ms = spu_cmpeq( tj, i00 );
#ifdef LJ
        r2cut   = spu_sel( pt.r2cut   [c1+1], pt.r2cut   [c1], ms );
        ljsig   = spu_sel( pt.lj_sig  [c1+1], pt.lj_sig  [c1], ms );
        ljeps   = spu_sel( pt.lj_eps  [c1+1], pt.lj_eps  [c1], ms );
        ljshift = spu_sel( pt.lj_shift[c1+1], pt.lj_shift[c1], ms );
#else
        r2cut   = spu_sel( pt.r2cut  [c1+1], pt.r2cut  [c1], ms );
        ptbeg   = spu_sel( pt.begin  [c1+1], pt.begin  [c1], ms );
        ptstep  = spu_sel( pt.step   [c1+1], pt.step   [c1], ms );
        ptistep = spu_sel( pt.invstep[c1+1], pt.invstep[c1], ms );
        ptcol   = spu_sel( pt.tab    [c1+1], pt.tab    [c1], ms );
#endif
#endif

        /* distance vectors */
        d0  = spu_sub( qpos[j0], *pi );
        d1  = spu_sub( qpos[j1], *pi );
        d2  = spu_sub( qpos[j2], *pi );
        d3  = spu_sub( qpos[j3], *pi );

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

#ifdef LJ
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
#else
        ms1 = spu_cmpgt( r2cut, r2 );  /* cutoff mask */
        ms2 = spu_cmpgt( r2, f001 );   /* mask zeros */
        r2  = spu_sel( r2cut, r2, spu_and( ms1, ms2 ) );
        r2  = spu_msub(r2, ptistep, ptbeg );
        kk  = spu_convtu( r2, 0 );
        b   = spu_sub( r2, spu_convtf(kk,0) );
        a   = spu_sub( f10, b );
        a2  = spu_msub( a, a, f10 );
        b2  = spu_msub( b, b, f10 );
        a3  = spu_mul( a, a2 );
        b3  = spu_mul( b, b2 );
        st6 = spu_mul( ptstep, f6i ) ;

        kk  = spu_add( kk, ptcol );
        pt0 = pt.data[spu_extract(kk,0)];
        pt1 = pt.data[spu_extract(kk,1)];
        pt2 = pt.data[spu_extract(kk,2)];
        pt3 = pt.data[spu_extract(kk,3)];

        tmp0 = spu_shuffle(pt0, pt2, vpat1);
        tmp1 = spu_shuffle(pt1, pt3, vpat1);
        tmp2 = spu_shuffle(pt0, pt2, vpat2);
        tmp3 = spu_shuffle(pt1, pt3, vpat2);

        p1 = spu_shuffle(tmp0, tmp1, vpat3);
        y1 = spu_shuffle(tmp0, tmp1, vpat4);
        p2 = spu_shuffle(tmp2, tmp3, vpat3);
        y2 = spu_shuffle(tmp2, tmp3, vpat4);

	tmp1 = spu_madd( a, p1, spu_mul(b,p2) );
	tmp2 = spu_madd( a3, y1, spu_mul(b3,y2) );
        tmp3 = spu_mul( st6, ptstep );
        pot  = spu_madd( tmp2, tmp3, tmp1 );

        tmp1 = spu_mul( spu_sub(p2, p1), ptistep );
        tmp2 = spu_madd( f30, b2, f20 );
        tmp3 = spu_madd( f30, a2, f20 );
        tmp4 = spu_msub( tmp2, y2, spu_mul(tmp3, y1) );
        tmp5 = spu_madd( tmp4, st6, tmp1 );
        grad = spu_mul( f20, tmp5 );
#endif

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
      pott = spu_add( pots, pott );
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

    old_q  = (q == p) ? &(buf[2]) : q;
    q      = next_q;
    next_q = old_q;

  } while (EXPECT_TRUE_(m<NNBCELL));

  /* set contribution to total virial */
#ifndef AR
  vir = spu_mul( vir, f05 );  /* avoid double counting */
#endif
  vir = spu_add( vir, spu_rlqwbyte( vir, 8 ) );
  vir = spu_add( vir, spu_rlqwbyte( vir, 4 ) );

#ifdef PPUSUM
  wp->virial   = -spu_extract( vir,  0 );
  wp->totpot   =  spu_extract( pott, 0 );
#else
  psum.virial -=  spu_extract( vir,  0 );
  psum.totpot +=  spu_extract( pott, 0 );
#endif

#ifdef SPU_INT
  move_atoms_spu(wp, p);
#endif

  /* Init DMA back to main memory. No need to wait for it to complete here,
     as we will do that in the main (work) loop. */
  init_return(p, wp->cell_dta, otag);
}


/******************************************************************************
*
*  calc_tb: Neighbor tables on SPU
*
******************************************************************************/

void calc_tb_direct(wp_t *wp,
                    /*
                    void* const tempbuf, unsigned const tempbuf_len,
                    void* const resbuf,  unsigned const resbuf_len,
                    */
                    cell_dta_t* const buf,
                    unsigned const otag)
{
  vector float d, *qpos;
  float const cellsz = wp->totpot;
  int const inc_short = 128 / sizeof(short) - 1;
  int m, next_m, nn=0;

  cell_dta_t *p, *q, *next_q, *old_q;


  /* The initial fetch to buffer q */
  p = q = buf;
  init_fetch_tb(q, wp->cell_dta, otag);

  /* Set pointers to the remaining buffers */
  next_q = buf+1;

  /* for each neighbor cell */
  m = 0;
  do {

    int l=0, n=0, i;

    wait_tag(otag);

    wp->cell_dta[m].tb_ea = wp->cell_dta[0].tb_ea + nn * sizeof(short); 

    next_m = m;
    do {
      next_m++;
      if (m<NNBCELL) {
        wp->cell_dta[next_m].tb_ea = wp->cell_dta[next_m-1].tb_ea; 
      }
    } while ((wp->cell_dta[next_m].n==0) && (next_m<NNBCELL));

    if (next_m<NNBCELL) {
      init_fetch_tb(next_q, wp->cell_dta + next_m, otag);
    }

    /* positions in q as vector float */
    qpos = (vector float *) q->pos;

    /* for each atom in cell */
    for ( i=0; EXPECT_TRUE_(i<p->n);  ++i ) {

      vector float* const pi = ((vector float *) p->pos) + i;
      int    jstart, j, rr;

      /* for each neighboring atom */
      q->ti[l] = n;
#ifdef AR
      jstart = (m==0) ? i+1 : 0;
#else
      jstart = 0;
#endif
      for ( j=jstart;   EXPECT_TRUE_(j<q->n);   ++j ) {
        float r2;
#ifndef AR
        if ((m==0) && (i==j)) continue;
#endif
        d  = spu_sub( qpos[j], *pi );
        d  = spu_mul( d, d );
        r2 = spu_extract(d,0) + spu_extract(d,1) + spu_extract(d,2); 
        if (r2 < cellsz) q->tb[n++] = j;
      }
      q->ti[l+1] = n - q->ti[l];
      l += 2;

      /* if n is not divisible by 4, pad with copies of q->n */
      rr = n % 4;
      if (rr>0) for (j=rr; j<4; j++) q->tb[n++] = q->n;
    }

    /* enlarge n to next 128 byte boundary */
    n = (n + inc_short) & (~inc_short); 

    nn += n;
    if ( EXPECT_FALSE_(nn > wp->nb_max) ) {
      wp->flag = -1;
      printf("nb_max = %d, nn = %d\n", wp->nb_max, nn);
      return;
    }

    init_return_tb(q, wp->cell_dta + m, wp->ti_len, n, otag);

    m      = next_m;
    old_q  = (q == p) ? &(buf[2]) : q;
    q      = next_q;
    next_q = old_q;

  } while (EXPECT_TRUE_(m<NNBCELL));

  wp->flag = nn;

}

#endif  /* not ON_PPU */
