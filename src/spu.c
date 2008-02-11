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
* spu.c -- main program on SPU
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

/* ISO C std. headers */
/* #include <stddef.h> */
#include <stdio.h>

/* SPU headers */
#include <spu_mfcio.h>
#include <spu_internals.h>

/* SDK Headers containing inline functions/macros */
/*
#include <cond_init.h>
#include <cond_wait.h>
#include <cond_signal.h>
*/

/* Local/program specific headers */
#include "imd_cbe.h"

/* potential data structure */
pt_t ALIGNED_(16,pt);


/* Suffix (kilobyte) */
enum { Ki=1024u };

/* Maximum decrementer value */
/* static unsigned int const DecMax = -1; */
#define DecMax ((unsigned int)(-1))


/* Argument types for main function */
typedef unsigned long long  ui64_t;


/* Maximum */
#define MAX(a,b) ((a)>(b) ? (a) : (b))

static INLINE_ unsigned uimax(unsigned const a, unsigned const b) {
    return MAX(a,b);
}


#if defined(CBE_DIRECT)  /* "Direct" case */


/* Set result buffer (force) */
static INLINE_ void set_res_buffers(cell_dta_t* const pbuf, int const n_max,
                                    void* const resbuf, unsigned const resbuf_len) {
    /* Size of vector arrays (e.g. force) */
    int const svec = iceil128(n_max * sizeof(vector float));

    /* forces arrays will only be used for pbuf[0] */
    /* pbuf[1].force=pbuf[2].force=0; */

    if ( EXPECT_TRUE_(resbuf_len >= svec) ) {
        /* At the moment, only pbuf[0].forces will store results */
        pbuf[0].force = (float*)resbuf;
    }
    else {
        /* fprintf(stderr, "set_buffers: Can not allocate %d bytes for result buffers\n", svec); */
       pbuf[0].force = 0;
    }

}

/* Set temp. buffer memory (which will be needed during calculation only) */
static void set_tmp_buffers(cell_dta_t* const pbuf,
                            int const n_max, int const len_max,
                            void* const tmpbuf, unsigned const tmpbuf_len)
{
    /* Size of vector arrays (e.g. position) */
    int const svec = iceil128(  n_max *     sizeof(vector float));
    /* Sizes of the scalar members */
    int const styp = iceil128(  n_max *     sizeof(int));
    int const sti  = iceil128(  n_max * 2 * sizeof(int));
    int const stb  = iceil128(len_max *     sizeof(short));
    int const stb2 = iceil128(len_max/2 *   sizeof(short));

    /* Total number of bytes needed for temp. storage */
    int const stottmp = 3*(svec+styp)+2*(sti+stb);


    /* fprintf(stdout, "set_tmp_buffers: n_max=%d, len_max=%d\n", n_max, len_max); fflush(stdout); */

    /* Enoug temp space available? */
    if ( EXPECT_TRUE_(tmpbuf_len >= stottmp) ) {
        /* Allocation starts here */
        register unsigned char* loc = (unsigned char*)tmpbuf;

        /* Pos arrays */
        pbuf[0].pos = (float*)loc;
        loc+=svec;
        pbuf[1].pos = (float*)loc;
        loc+=svec;
        pbuf[2].pos = (float*)loc;
        loc+=svec;

        /* typ arrays */
        pbuf[0].typ = (int*)loc;
        loc+=styp;
        pbuf[1].typ = (int*)loc;
        loc+=styp;
        pbuf[2].typ = (int*)loc;
        loc+=styp;

        /* ti arrays */
        pbuf[0].ti  = (int*)loc;
        loc+=sti;
        pbuf[1].ti  = (int*)loc;
        loc+=sti;
        pbuf[2].ti  = pbuf[0].ti;

        /* tb arrays */
        pbuf[0].tb  = (short*)loc;
        loc+=stb;
        pbuf[1].tb  = (short*)loc;
        loc+=stb2;  /* this one can be smaller */
        pbuf[2].tb  = pbuf[0].tb;
    }
    else {
        /*  fprintf(stderr, "set_buffers: Can not allocate %d bytes for temp buffers\n", 
	    stottmp); */
    }

}





/* Init get (always the same for wp & tb) */
static void init_get(void* ibuf, void* obuf, ea_t iea, unsigned const itag)
{
    /* Start DMA to input buffer usinga fence command, as the input buffer
       may have been used for a previous put command to system-memory */
    dma64(ibuf, iea, (sizeof (argbuf_t)), itag, MFC_GETF_CMD);
}



/* Temp. DMA buffer */
static unsigned char ALIGNED_(128, direct_tmp[80*Ki]);

/* Buffer "control block" for calc_direct routines */
static cell_dta_t direct_buf[3];

/* Last parameters used to setup tmp. buffers in bufcb */
static int last_n_max=0u;
static int last_len_max=0u;


/* Size of ouput buffers (used in main) */
enum { DIRECT_OBUFSZE=((unsigned)(30*Ki)) };



static void wp_direct_spusum(void* i, void* o,  ea_t unused, unsigned otag)
{
    /* Size of buffer for the results */
    enum { RESSZE=DIRECT_OBUFSZE-(sizeof (argbuf_t)) };

    /* Store ptr. to buffers in register */
    register cell_dta_t* const B = direct_buf;

    /* Ptrs. to input/summatiuon workpackages */
    register wp_t* const pwp  = (wp_t*)i;
    register wp_t* const psum = (wp_t*)o;

    /* Some length parameters from WP */
    register int const n_max   = pwp->n_max;
    register int const len_max = pwp->len_max;


    /* Set buffers */
    set_tmp_buffers(B, n_max, len_max,  &direct_tmp,  (sizeof direct_tmp));
    set_res_buffers(B, n_max,           ((argbuf_t*)o)+1, RESSZE);

    /* Calc. forces */
    calc_wp_direct(pwp, B, otag);


    /* Add up forces/virial in result buffer */
    psum->virial += pwp->virial;
    psum->totpot += pwp->totpot;
}


static void wp_direct(void* i, void* o, ea_t oea, unsigned otag)
{
    /* Store ptr. to buffers in register */
    register cell_dta_t* const B = direct_buf;

    /* Current work package */
    register wp_t* const pwp = ((wp_t*)i);

    /* Some length parameters from WP */
    register int const n_max   = pwp->n_max;
    register int const len_max = pwp->len_max;




    /* Set buffers  */
    set_tmp_buffers(B, n_max, len_max,  &direct_tmp, (sizeof direct_tmp));
    set_res_buffers(B, n_max,           o,           DIRECT_OBUFSZE);

    /* Calc forces */
    calc_wp_direct(pwp, B,   otag);


    /* Init outbound DMA: */
    dma64(i,  oea,  (sizeof (argbuf_t)),  otag, MFC_PUT_CMD);
}








static void tb_direct(void* i, void* o, ea_t oea, unsigned otag)
{
    /* Current WP */
    register wp_t* const pwp = ((wp_t*)i);

    /* Store ptr. to buffers in register */
    register cell_dta_t* const B = direct_buf;

    /* Some length parametrs form WP */
    register int const n_max   = pwp->n_max;
    register int const len_max = pwp->len_max;



    /* Set buffers  */
    set_tmp_buffers(B,  n_max, len_max,  &direct_tmp, (sizeof direct_tmp));
    set_res_buffers(B,  n_max,           o,           DIRECT_OBUFSZE);

    /* Calc. NBL */
    calc_tb_direct(pwp, B,  otag);

    /* DMA WP back */
    dma64(i, oea,  (sizeof(argbuf_t)),  otag, MFC_PUT_CMD);
}














#else   /* "Indirect" case */



/* Given an exch_t, build the corresponding wp_t, that is:
   Copy the scalar members
   Allocate memory for all the array members (using memory pointed to by a).
   Fetch all the input array data to LS using DMA (possibly a list DMA) with
   the specified tag
*/
static void start_create_wp(exch_t const* const exch,
                            void* const ibuf, unsigned const ibufsze,
                            void* const obuf, unsigned const obufsze,
                            unsigned const itag,
                            wp_t* const wp)
{
     /* Number of bytes needed for the pos/force vectors */
     unsigned const nvecbytes = (exch->n2) * sizeof(float[4]);

     /* Number of bytes remaining */
     unsigned       nrem;
     unsigned char* apos;

     /* Boundaries used for rounding up */
     register vector unsigned int const bndv = { 15,127,  0,0 };
     
     /* Vector of sizes:
        szev[0]=to be transferd, szev[1]=to be allocated (remaing components
        are unused). Access the scalars using the following macros
     */
     vector unsigned int szev;

#define ADV(v)  spu_extract((v),1)
#define XFR(v)  spu_extract((v),0)
 
     /* Number of bytes to advance/allocate in LS */
     register unsigned nadv;



     /* First allocate input buffers */
     nrem=ibufsze;
     apos=ibuf;

     /* pos */
     nadv=ADV(szev=CEILV(nvecbytes, bndv));
     if ( EXPECT_TRUE_(nrem>=nadv)  ) {
         mdma64(apos,  exch->pos, XFR(szev),  itag, MFC_GET_CMD);
         wp->pos = (void*)apos;
         nrem-=nadv;
         apos+=nadv;
     }
     else {
         fprintf(stderr, "Could not allocate %u bytes for pos\n", nadv);
         wp->pos = 0;
     }

     /* typ */
     nadv=ADV(szev=CEILV(exch->n2*sizeof(int),  bndv));
     if ( EXPECT_TRUE_(nrem>=nadv) ) {
         mdma64(apos,  exch->typ, XFR(szev),  itag, MFC_GET_CMD);
         wp->typ=(void*)apos;
         nrem-=nadv;
         apos+=nadv;
     }
     else {
         fprintf(stderr, "Could not allocate %u bytes for typ\n", nadv);
         wp->typ = 0;
     }

     /* ti */
     nadv = ADV(szev=CEILV(2*exch->n2*sizeof(int), bndv));
     if ( EXPECT_TRUE_(nrem>=nadv) ) {
         mdma64(apos,  exch->ti, XFR(szev),  itag, MFC_GET_CMD);
         wp->ti=(void*)apos;
         nrem-=nadv;
         apos+=nadv;
     }
     else {
         fprintf(stderr, "Could not allocate %u bytes for ti\n", nadv);
         wp->ti = 0;
      }

     /* tb */
     nadv=ADV(szev=CEILV(exch->len*sizeof(short), bndv));
     if ( EXPECT_TRUE_(nrem>=nadv) ) {
         mdma64(apos, exch->tb,  XFR(szev),   itag, MFC_GET_CMD);
         wp->tb=(void*)apos;
         nrem-=nadv;
         apos+=nadv;
     }
     else {
         fprintf(stderr, "Could not allocate %u bytes for tb\n", nadv);
         wp->tb = 0;
     }


     /* Now allocate output buffers */
     nrem = obufsze;
     apos = obuf;

     /* force (no DMA needed) */
     nadv=ADV(szev=CEILV(nvecbytes, bndv));
     if ( EXPECT_TRUE_(nrem>=nadv) ) {
          wp->force = (void*)apos;
          nrem -= nadv;
          apos += nadv;

     }
     else {
         fprintf(stderr, "Could not allocate %u bytes for forces\n", nadv);
         wp->force=0;
     }


#undef XFR
#undef ADV



     /* Now copy all scalar members (but totpot & virial) from exch to wp */
     wp->k       = exch->k;

     wp->n1      = exch->n1;
     wp->n1_max  = exch->n1_max;

     wp->n2      = exch->n2;
     wp->n2_max  = exch->n2_max;

     wp->len     = exch->len;
     wp->len_max = exch->len_max;
}




/* Some buffer size for "Direct case" */
enum { IBUFSZE=((unsigned)(70*Ki)), OBUFSZE=((unsigned)(16*Ki)) };

static void init_get(void* ibuf, void* obuf, 
                     ea_t iea,  unsigned const itag)
{
    /* Input exch_t & wp_t are DMAed in and are placed in the input buffer */
    exch_t* const pexchin = (exch_t*)ibuf;
    wp_t*   const pwp     = (wp_t*)(((argbuf_t*)ibuf)+1);
    

    /* Number of bytes remaining in the input/output buffers */
    enum { IREM=IBUFSZE-(sizeof (argbuf_t[2])),
           OREM=OBUFSZE-(sizeof (argbuf_t[1]))   };

    /* Tag/mask for short DMA of controlblock (exch_t) */
    enum { cbtag=(unsigned)5          };
    enum { cbmsk=(unsigned)(1<<cbtag) };

    /* First get the controll block (exch_t) from ea */
    dma64(pexchin, iea,  (sizeof (argbuf_t)),  cbtag, MFC_GET_CMD);
    wait_dma(cbmsk,  MFC_TAG_UPDATE_ALL);

    /* Now start creating the work package */
    start_create_wp(pexchin,
                    ((argbuf_t*)ibuf)+2, IREM,  ((argbuf_t*)obuf)+1, OREM,
                    itag, pwp);
}




static void init_put(void* obuf,  ea_t oea, unsigned otag)
{
    /* The exch to be DMAed out is the first element of output buffer */
    exch_t* pexchout = (exch_t*)obuf;

    /* DMA back the force (LS address is stored in exch.pos_ea[0]) */
    mdma64((void*)(pexchout->pos[0]), pexchout->force,
           uiceil16((pexchout->n2) * sizeof(float[4])),
           otag, MFC_PUT_CMD
          );

    /* DMA back the exch_t containing virial & energy */
    dma64(pexchout, oea,  (sizeof (argbuf_t)),  otag, MFC_PUT_CMD);
}



static void wrk(void* ibuf, void* obuf,  unsigned otag)
{
   /* The workpackage is located in input buffer followed by the input exch */
   wp_t*         const pwp = (wp_t*)ibuf;
   exch_t const* const pin = (exch_t*)(((argbuf_t*)ibuf)+1);

   /* Output exch_t is the first element in output buffer 
      (followed by the forces)                      */
   exch_t*       const pout = (exch_t*)obuf;



   /* The main work to be done */
   calc_wp(pwp);


   /* Copy some members from input buffers (exch & wp) */
   /* The scalars calculated by calc_wp */
   pout->totpot = pwp->totpot;
   pout->virial = pwp->virial;

   /* Serial number */
   pout->k      = pwp->k;

   /* n2 is also needed as it determines the length of the force array 
      to be DMAed back */
   pout->n2     = pwp->n2;

   /* Also copy EA of force */
   pout->force[0] = pin->force[0];
   pout->force[1] = pin->force[1];

   /* Store LS pointer to force as unsigned in exch.pos_ea[0] */
   pout->pos[0] = (unsigned)(pwp->force);
}



#endif /* CBE_DIRECT */













/* C.f. Prog. Handbook p. 687  */

/* Double bufferd workloop for a list of given eas */
/* After this function returns there might still be some
   outbound DMAs pending. The unsigned returned by this function is
   a bit mask of the tags of pending DMAs.
 */
static unsigned workloop_ea_stat(
   /* List of n input & output EAs */
   register ea_t const* const ieas, register ea_t const* const oeas, register unsigned int n,
   /* Input/output buffers */
   register void* const* const ib,  register void* const* const ob,
   /* Tags & masks */
   register unsigned int const* const t, register unsigned int const* const m,
   /* Functions to get, to process and to write back data */
   register void (* const getf)(void* ibuf, void* obuf, ea_t const iea, unsigned it),
   register void (* const wrkf)(void* ibuf, void* obuf, ea_t const oea, unsigned ot, tick_t wtime)
)
{
    /* Next buffer & EA indices */
    register unsigned int nxtbf, nxtea;
    /* Current buffer & EA indices */
    register unsigned int curbf=0u, curea=0u;


    /* Timing */
    register tick_t dt=0, t0=0;
    tick_t tstart[2], tend[2];

    getf(ib[curbf],ob[curbf], ieas[curea], t[curbf]);
    for (;;) {
        /* Finsish work loop? */
        if ( EXPECT_FALSE_(0u == (--n)) ) {
  	    break;
        }

        /* Next buffer/EA indices: 0,1 -> 1,0 ("swap" 0 and 1) */
        nxtbf=curbf^1u;   
        nxtea=curea+1u;

        /* Get next buffer */
        getf(ib[nxtbf],ob[nxtbf],   ieas[nxtea], t[nxtbf]);

         
        /* Wait for current input/output buffers to become available */
        TDIFF( {
        spu_writech(MFC_WrTagMask,   m[curbf]);
        spu_writech(MFC_WrTagUpdate, MFC_TAG_UPDATE_ALL);
        (void)(spu_readch(MFC_RdTagStat));
        }, t0,dt,  ticks, tick_diff
        )

        /* Work on current buffer */
        wrkf(ib[curbf],ob[curbf],     oeas[curea], t[curbf], dt);

        /* Move on to next buffer/EA */
        curbf=nxtbf;
        curea=nxtea;
    }


    /* Wait for last buffer to become available */
    TDIFF( {
    spu_writech(MFC_WrTagMask,   m[curbf]);
    spu_writech(MFC_WrTagUpdate, MFC_TAG_UPDATE_ALL);
    (void)(spu_readch(MFC_RdTagStat));
    },  t0,dt,  ticks,tick_diff
    )

    /* Work on buffer */
    wrkf(ib[curbf],ob[curbf],  oeas[curea], t[curbf], dt);

    /* Return the mask client code should wait for */
    return m[curbf];
}


static unsigned workloop_ea(
   /* List of n input & output EAs */
   register ea_t const* const ieas, register ea_t const* const oeas, register unsigned int n,
   /* Input/output buffers */
   register void* const* const ib,  register void* const* const ob,
   /* Tags & masks */
   register unsigned int const* const t, register unsigned int const* const m,
   /* Functions to get, to process and to write back data */
   register void (* const getf)(void* ibuf, void* obuf, ea_t const iea, unsigned it),
   register void (* const wrkf)(void* ibuf, void* obuf, ea_t const oea, unsigned ot)
)
{
    /* Next buffer & EA indices */
    register unsigned int nxtbf, nxtea;
    /* Current buffer & EA indices */
    register unsigned int curbf=0u, curea=0u;


    /* Timing */
    register tick_t dt=0, t0=0;

    getf(ib[curbf],ob[curbf], ieas[curea], t[curbf]);
    for (;;) {
        /* Finsish work loop? */
        if ( EXPECT_FALSE_(0u == (--n)) ) {
  	    break;
        }

        /* Next buffer/EA indices: 0,1 -> 1,0 ("swap" 0 and 1) */
        nxtbf=curbf^1u;   
        nxtea=curea+1u;

        /* Get next buffer */
        getf(ib[nxtbf],ob[nxtbf],   ieas[nxtea], t[nxtbf]);

         
        /* Wait for current input/output buffers to become available */
        spu_writech(MFC_WrTagMask,   m[curbf]);
        spu_writech(MFC_WrTagUpdate, MFC_TAG_UPDATE_ALL);
        (void)(spu_readch(MFC_RdTagStat));

        /* Work on current buffer */
        wrkf(ib[curbf],ob[curbf],     oeas[curea], t[curbf]);

        /* Move on to next buffer/EA */
        curbf=nxtbf;
        curea=nxtea;
    }


    /* Wait for last buffer to become available */
    spu_writech(MFC_WrTagMask,   m[curbf]);
    spu_writech(MFC_WrTagUpdate, MFC_TAG_UPDATE_ALL);
    (void)(spu_readch(MFC_RdTagStat));


    /* Work on buffer */
    wrkf(ib[curbf],ob[curbf],  oeas[curea], t[curbf]);

    /* Return the mask client code should wait for */
    return m[curbf];
}


/* Start initialization (pt only, at the moment) */
static void start_init(ea_t const envea, unsigned int const itag)
{
  /* Buffer for potential data */
#ifdef LJ
  enum { memsize=16u };   /* for ntypes <= 2 */
#else
  enum { memsize=384u };  /* for ntypes <= 2 and <= 120 potential lines */
  int i, j;
#endif
  static vector float ALIGNED_(128, pbuf[memsize]);
  unsigned int off=0, inc;

  /* Get the potential control block - 16 bytes is enough */
  dma64( &pt, envea, 16, itag, MFC_GET_CMD);
  wait_dma((1u<<itag), MFC_TAG_UPDATE_ALL);

  /* Get the potential data */
  dma64(pbuf, pt.ea, sizeof(pbuf), itag, MFC_GET_CMD);

  /* Set some pointers to the data */
#ifdef LJ
  inc = pt.ntypes * pt.ntypes;
  pt.r2cut    = pbuf;            off += inc;
  pt.lj_sig   = pt.r2cut + off;  off += inc;
  pt.lj_eps   = pt.r2cut + off;  off += inc;
  pt.lj_shift = pt.r2cut + off;  off += inc;
#else
  inc = pt.ntypes * pt.ntypes;
  pt.data    = pbuf;
  pt.begin   = pt.data + off;  off += inc;
  pt.r2cut   = pt.data + off;  off += inc;
  pt.step    = pt.data + off;  off += inc;
  pt.invstep = pt.data + off;  off += inc;
  for (i=0; i<pt.ntypes; i++)
    for (j=i; j<pt.ntypes; j++) {
      pt.tab[i*pt.ntypes+j] = spu_splats( off );
      pt.tab[j*pt.ntypes+i] = spu_splats( off );
      off += pt.nsteps + 1;
    }
#endif

  if ( EXPECT_FALSE_(off > memsize) ) {
    printf("memory overflow in potential data\n");
  }

}



/* Make EA list in eas (npairs elements),
   with initial value eastart
*/
static INLINE_ void arg_eas(register ea_t cur64, register ea_t const stride64,
                            register ea_t* p64, register unsigned int k)
{
    for (  ;       EXPECT_TRUE_(k>0u);      --k,  ++p64, cur64+=stride64  ) {
        (*p64)=cur64;
    }
}



/* argp "points" to an array of lengt N_ARGBUF of argbuf_t */
int main(ui64_t const id, ui64_t const argp, ui64_t const envp)
{
    /* Some tags/mask used for. misc purposes */
    enum { misctag = ((unsigned)1) };
    enum { miscmsk = ((unsigned)(1u<<misctag)) };

    /* Stride in argbuf EA space */
    enum { EAstride=((unsigned)(sizeof (argbuf_t))) };

    /* Effective addresses (pairs of unsigneds) of the arguments
       in main memory */
    ea_t ALIGNED_(8, eabuf[N_ARGBUF]), ALIGNED_(16, env[N_ENVEA]);

    /* Some timing stuff */
    register tick_t t0, dt;

    /* Get the environment list */
    mdma64(env,  envp,  (sizeof (ea_t[N_ENVEA])),   misctag, MFC_GET_CMD);
    (void)(wait_dma(miscmsk, MFC_TAG_UPDATE_ALL));

    /* Fetch env and start creating the corresponding pt_t (Using DMA) */
    start_init(env[0],  misctag);

    /* Set EAs of the buffers in main memory */
    arg_eas(argp,  EAstride, eabuf,N_ARGBUF);


    /* Wait for init. DMA to finish */
    (void)(wait_dma(miscmsk, MFC_TAG_UPDATE_ALL));


    /* _cond_wait(env[1], env[2]); */


    /* Set decrementer to max. value for further timing */
    spu_writech(SPU_WrDec, DecMax);

    for (;;) {
        /* Read work mode / flag */
        unsigned int const mode = spu_read_in_mbox();

	/* Output buffers */
        static unsigned char ALIGNED_(128, obuf0[DIRECT_OBUFSZE]);
        static unsigned char ALIGNED_(128, obuf1[DIRECT_OBUFSZE]);
        static void* const ob[] = { &obuf0, &obuf1 };


        /* mode==0: No more work   */
        if ( EXPECT_FALSE_(0u==mode) ) {
 	    break;
        }
        else {
 	   /* "Functors passed to workloop" */
 	   register void (* getfunc)(void* i, void* o, ea_t iea, unsigned it);
           register void (* wrkfunc)(void* i, void* o, ea_t oea, unsigned ot);

           /* Input values from mailbox */
           unsigned offs,narg;


           /* Mode dependent setup */
           switch (mode) {
	       case DOWP : {
		  getfunc=init_get;
                  wrkfunc=wp_direct_spusum;
                  /* Set accumulation variables */
                  ((wp_t*)&obuf0)->virial = ((wp_t*)&obuf1)->virial =
		  ((wp_t*)&obuf0)->totpot = ((wp_t*)&obuf1)->totpot = 0.0f;
                  break;
               }
	       case DOTB : {
		  getfunc=init_get;
                  wrkfunc=tb_direct;
                  break;
               }
	       default: {
		  break;
               }
           };

  	   for(;;) {
               /* Input buffers */
               static argbuf_t ALIGNED_(128, ibuf0);
               static argbuf_t ALIGNED_(128, ibuf1);
               /* Pointers to the beginning of the buffers */
               static void* const ib[] = { &ibuf0, &ibuf1 };
               
               /* Tags & the corresponding masks to be used */
               static unsigned int const tags[] = { 0,1 };
               static unsigned int const masks[] = { 1<<0, 1<<1 };

               /* Mask of pending DMAs */
               unsigned omsk;

               /* Input/output EAs */
               ea_t const* eas;

               /* Read #of WPs, if number is >0 also read offset */
               narg=spu_read_in_mbox();
               if ( EXPECT_FALSE_(0u==narg) ) {
		  break;
               }
               offs=spu_read_in_mbox();


               /* Debugging output */
               /* fprintf(stdout, "%u WPs at offset %u in mode %u\n", 
                  narg,offs,mode); fflush(stdout); */

               /* Process the arguments */
               eas = eabuf+offs;
               omsk = workloop_ea(eas,eas,narg, ib,ob, tags,masks, 
                                  getfunc,wrkfunc);

               /* Wait for all outbound DMAs */
               wait_dma(omsk, MFC_TAG_UPDATE_ALL);

               /* Write back offset */
               spu_write_out_mbox(offs);
           }

           /* Mode dependent Postprocessing */
           if ( DOWP==mode ) {
	      /* The resulting wp is DMAed back with the following tag/mask */
	      enum { t=(unsigned)0       };
              enum { m=(unsigned)(1u<<t) };

              /* Offset of result buffer */
              enum { RESOFFS  = (unsigned)0 };

	      /* Total energy/virial in obuf0 */
	      ((wp_t*)&obuf0)->totpot += ((wp_t const*)&obuf1)->totpot;
              ((wp_t*)&obuf0)->virial += ((wp_t const*)&obuf1)->virial;

              /* DMA results back */
              dma64(&obuf0, eabuf[RESOFFS], sizeof(argbuf_t), t, MFC_PUT_CMD);
              wait_dma(m, MFC_TAG_UPDATE_ALL);


              /* Tell PPU where to pick up accu */
              spu_write_out_mbox(RESOFFS);
           }
           if ( DOTB==mode ) {
	      enum { RESOFFS=((unsigned)1) };
	      /* Just acknoledge */
	      spu_write_out_mbox(RESOFFS);
           }
        }
    }

    return 0;
}
