
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
#include <malloc_align.h>
#include <free_align.h>

/*
#include <cond_init.h>
#include <cond_wait.h>
#include <cond_signal.h>
*/

/* Local/program specific headers */
#include "config.h"
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
                            mem_buf_t *tmpbuf)
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

    /* (re-)allocate tmpbuf */
    if ( EXPECT_FALSE_(tmpbuf->len < stottmp) ) {
      _free_align(tmpbuf->data);
      tmpbuf->data = (unsigned char *) _malloc_align( stottmp, 7 );
      if ( EXPECT_TRUE_((int)tmpbuf->data) ) {
        tmpbuf->len = stottmp;
      }
      else {
        printf("tmpbuf allocation failed, %d bytes requested\n", stottmp);
        tmpbuf->len = 0;
      }
    }

    /* Enough temp space available? */
    if ( EXPECT_TRUE_(tmpbuf->len >= stottmp) ) {
        /* Allocation starts here */
        register unsigned char* loc = (unsigned char *) tmpbuf->data;

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
}





/* Init get (always the same for wp & tb) */
static void init_get(void* ibuf, void* obuf, ea_t iea, unsigned const itag)
{
    /* Start DMA to input buffer usinga fence command, as the input buffer
       may have been used for a previous put command to system-memory */
    dma64(ibuf, iea, (sizeof (argbuf_t)), itag, MFC_GETF_CMD);
}



/* Temp. DMA buffer */
static mem_buf_t tmpbuf = {NULL, 0};

/* Buffer "control block" for calc_direct routines */
static cell_dta_t direct_buf[3];

/* Last parameters used to setup tmp. buffers in bufcb */
static int last_n_max=0u;
static int last_len_max=0u;


/* Size of ouput buffers (used in main) */
enum { DIRECT_OBUFSZE=((unsigned)(5*Ki)) };



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
    set_tmp_buffers(B, n_max, len_max, &tmpbuf );
    set_res_buffers(B, n_max, ((argbuf_t*)o)+1, RESSZE);

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
    set_tmp_buffers(B, n_max, len_max, &tmpbuf );
    set_res_buffers(B, n_max, o, DIRECT_OBUFSZE);

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
    set_tmp_buffers(B, n_max, len_max, &tmpbuf );
    set_res_buffers(B, n_max, o, DIRECT_OBUFSZE);

    /* Calc. NBL */
    calc_tb_direct(pwp, B,  otag);

    /* DMA WP back */
    dma64(i, oea,  (sizeof(argbuf_t)),  otag, MFC_PUT_CMD);
}


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
