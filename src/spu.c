/* ISO C std. headers */
/* #include <stddef.h> */
#include <stdio.h>



/* SPU headers */
#include <spu_mfcio.h>
#include <spu_internals.h>

/* Local/program specific headers */
#include "imd_cbe.h"






/* Suffix (kilobyte) */
enum { Ki=1024u };



/* Argument types for main function */
typedef unsigned long long ui64_t;


/* Take *p as a 64bit (in general ty) quantity & add n to it */
#define INC(ty,p,n)  ((*((ty*)(p))) += (n))
#define INC64(p,n)   INC(ui64_t, (p),(n))





/* First/past the end byte position of member m in structure like object s */
#define MEMB_FRST_POS(s,m) (((unsigned char const*)&((s).m)) - ((unsigned char const*)&(s)))
#define MEMB_LAST_POS(s,m) (MEMB_FRST_POS((s),m) + (sizeof ((s).m)))

/* Maximum */
#define MAX(a,b) ((a)>(b) ? (a) : (b))

static INLINE_ unsigned uimax(unsigned const a, unsigned const b) {
    return MAX(a,b);
}




#if defined(CBE_DIRECT)  /* "Direct" case */

/* Size of ouputbuffer */
enum { DIRECT_OBUFSZE=((unsigned)(30*Ki)) };



/**/
static void init_get_direct(void* ibuf, void* obuf,
                            unsigned const* iea, unsigned const itag)
{
    /* Start DMA to input buffer */
    dma64(ibuf, iea,  (sizeof (argbuf_t)),  itag,  MFC_GET_CMD);
}

static void init_put_direct(void* obuf,
                            unsigned const* oea,  unsigned otag)
{
    /* Results are now in xfer buffer, so transfer it (at least partially) */
    dma64(obuf, oea,  (sizeof (argbuf_t)),  otag, MFC_PUT_CMD);
}


static INLINE_ void vmemcpy(void const* const psrc0, void* const pdst0, unsigned const n0)
{
    typedef vector unsigned char Tv;
    register Tv const *psrc = (Tv const*)psrc0;
    register Tv       *pdst = (Tv*)pdst;
    register unsigned  nv   = n0/(sizeof (Tv));

    for (  ;   (nv>0u);    --nv, ++pdst, ++psrc ) {
        *pdst = *psrc;
    }
}

static void wrk_direct(void* i, void* o,  unsigned otag)
{
    /* Temp buffer */
    static unsigned char ALIGNED_(128, tmp[90*Ki]);

    /* Input work package */
    wp_t* const pwp  = (wp_t*)i;

    /* Output work package... */
    wp_t* const pout = (wp_t*)o;
    /* followed by some memory fore forces which are also stored in 
       output buffer */
    void*       pres = ((argbuf_t*)o)+1;

    /* Number of bytes remaining in result buffer */
    enum { RESSZE = (unsigned)(DIRECT_OBUFSZE-(sizeof (argbuf_t))) };


    /* fprintf(stdout, "wrk_direct: %u bytes remaining for results\n",  ressze); fflush(stdout);  */


    /* The main work to be done: Calc. the input WP */
    switch ( pwp->flag ) {
   	   case 1: {
 	       /* NBL */
   	       calc_tb_direct(pwp,  &tmp, (sizeof tmp),  pres, RESSZE, otag);

               /* Copy the entire WP */
               (*pout) = (*pwp);

               break;
	   }
   	   case 2: {
               /* Forces, virial & energy */
  	       calc_wp_direct(pwp,  &tmp, (sizeof tmp),  pres, RESSZE, otag);

               /* Copy some scalar to output WP */
               pout->k    = pwp->k;
               pout->flag = pwp->flag;

               pout->virial = pwp->virial;
               pout->totpot = pwp->totpot;

               break;
           }
 	   default: {
	       /* Error */
               fprintf(stderr, "unknown SPU task: \n");
               break;
          }
    }
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
     vector unsigned szev;

#define ADV(v)  spu_extract((v),1)
#define XFR(v)  spu_extract((v),0)
 
     /* Number of bytes to advance/allocate in LS */
     register unsigned nadv;



     /* First allocate input buffers */
     nrem=ibuflen;
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
     nrem = obuflen;
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
         pres_ea[0] = pres_ea[1]  =  pres_sze[0]=0;
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




/**/
enum { IBUFSZE=((unsigned)(70*Ki)), OBUFSZE=((unsigned)(16*Ki)) };

static void init_get(void* ibuf, void* obuf, 
                     unsigned const* iea,  unsigned const itag)
{
    /* Input exch_t & wp_t are DMAed in and are placed in the input buffer */
    exch_t* const pexchin = (exch_t*)ibuf;
    wp_t*   const pwp     = (wp_t)(((argbuf_t*)ibuf)+1);
    

    /* Number of bytes remaining in the input/output buffers */
    enum { IREM=IBUFSZE-(sizeof (argbuf_t[2])),
           OREM=OBUFSZE-(sizeof (argbuf_t[1]))   };

    /* Tag/mask for short DMA of controlblock (exch_t) */
    enum { cbtag=(unsigned)5          };
    enum { cbmsk=(unsigned)(1<<cbtag) };

    /* First get the controll block (exch_t) from ea */
    dma64(exchin, ea,  (sizeof (argbuf_t)),  cbtag, MFC_GET_CMD)
    wait_dma(cbmsk,  MFC_STATUS_UPDATE_ALL);

    /* Now start creating the work package */
    start_create_wp(pexchin,
                    ((argbuf_t*)ibuf)+2, IREM,  ((argbuf_t*)obuf)+1, OREM,
                    itag, pwp);
}




static void init_put(void* obuf,  unsigned const* oea, unsigned otag)
{
    /* The exch to be DMAed out is the first element of output buffer */
    exch_t* pexchout = (exch_t*)obuf;

    /* DMA back the force (LS address is stored in exch.pos_ea[0]) */
    mdma64((void*)(pexchout->pos_ea[0]), pexchout->force_ea,
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
   exch_t const* const pin = (exch_t*)(((argbuf_t)ibuf)+1);

   /* Output exch_t is the first element in output buffer 
      (followed by the forces)                      */
   exch_t*       const pout = (exch_t*)obuf;



   /* The main work to be done */
   calc_wp(pwp);


   /* Copy some members from input buffers (exch & wp) */
   /* The scalars calculated by calc_wp */
   pout->totpot = wp->totpot;
   pout->virial = wp->virial;

   /* Serial number */
   pout->k      = wp->k;

   /* n2 is also needed as it determines the length of the force array 
      to be DMAed back */
   pout->n2     = wp->n2;

   /* Also copy EA of force */
   pout->force_ea[0] = pin->force_ea[0];
   pout->force_ea[1] = pin->force_ea[1];

   /* Store LS pointer to force as unsigned in exch.pos_ea[0] */
   pout->pos_ea[0] = (unsigned)(wp->force);
}



#endif /* CBE_DIRECT */





/* Double bufferd workloop for list of given eas */
static int workloop_ea(/* List of n input & output EAs */
                       unsigned const* const ineas,
                       unsigned const* const outeas,
                       unsigned const n,
                       /* Input/output buffers */
                       void* const ib0, void* const ib1,
                       void* const ob0, void* const ob1,
                       /* Functions to get, to process and to write back data */
                       void (* const getf)(void* ibuf, void* obuf,
                                           unsigned const* iea, unsigned it),
                       void (* const wrkf)(void* ibuf, void* obuf, unsigned ot),
                       void (* const putf)(void* obuf, unsigned const* oea, unsigned ot)
		      )
{
    /* DMA Tags/masks */
    enum { t0=(unsigned)0,       t1=(unsigned)1       };
    enum { m0=(unsigned)(1<<t0), m1=(unsigned)(1<<t1) };
    enum { m01=(unsigned)(m0|m1) };



    /* Number of WPs remaining to be DMAed in, number of WPs done */
    register unsigned nrem  = n;
    register unsigned ndone = 0;

    /* Pointer to eff. addr. pairs */
    register unsigned const* iea0 = ineas;
    register unsigned const* iea1 = iea0+2;
    register unsigned const* oea0 = outeas;
    register unsigned const* oea1 = oea0+2;



    /* As we use b0 & b1 for DMA, we may assume that they are (at least)
       16-byte aligned. */


    /* Initial DMA to buffer 0 (if there are WPs remaining) */
    if ( EXPECT_TRUE_(nrem>0) ) {
        getf(ib0,ob0,  iea0,t0);
        --nrem;
    }
    else {
        /* Return right now, if there are no WPs at all */
        return 0;
    }

    for(;;) {
        /* Stride for walking thru EA list  */
        enum { eastride=((unsigned)(2*2)) };


        /* Init DMA to buffer 1, if there are WPs remaining */
        if ( EXPECT_TRUE_(nrem>0) ) {
           getf(ib1,ob1, iea1,t1);
           --nrem;
        }

        /* Wait for input and ouput buffers 0 to be available */
        wait_dma(m0, MFC_TAG_UPDATE_ALL);
        wrkf(ib0,ob0,        t0);
        putf(ob0,      oea0, t0);

        /* All WPs done? */
        if ( EXPECT_FALSE_(++ndone == n) ) {
            break;
        }

        /* Next EAs for buffer 0 */
  	iea0+=eastride;
        oea0+=eastride;

        /* Again, init DMA to buffer 0, if there are WPs remaining
         */
        if ( EXPECT_TRUE_(nrem>0) ) {
  	    getf(ib0,ob0, iea0,t0);
            --nrem;
        }

        /* Wait for, process & start DMAing back buffer 1*/
        wait_dma(m1, MFC_TAG_UPDATE_ALL);
        wrkf(ib1,ob1,       t1);
        putf(ob1,     oea1, t1);

        /* All WPs done? */
        if ( EXPECT_FALSE_(++ndone == n) ) {
  	   break;
        }

        /* Next EA for buffer 1 */
        iea1+=eastride;
        oea1+=eastride;
    }

    /* There may be still some (outbound) DMA which is not yet completed */
    wait_dma(m01, MFC_TAG_UPDATE_ALL);

    return 0;
}




static int workloop(unsigned const* const ea0,
                    void* const b, unsigned const bsze,
                    void (* const getf)(void* buf, unsigned sze,
                                        unsigned const* ea, unsigned it),
                    void (* const wrkf)(void* buf, unsigned sze,
                                        unsigned it, unsigned ot),
                    void (* const putf)(void* buf, unsigned sze,
                                        unsigned const* ea, unsigned ot)
     )
{
    /**/
    for (;;) {
        /* Tag & mask */
        enum { itag=((unsigned)0),          otag=((unsigned)1) };
        enum { imsk=((unsigned)(1u<<itag)), omsk=((unsigned)(1u<<otag)) };

        /* result from inbox (command or number of wps to be processed) */
        unsigned ibx;

        /* Wait for message from PPU */
        if ( EXPECT_FALSE_(SPUEXIT == (ibx = spu_read_in_mbox())) ) {
  	    break;
        }

        getf(b,bsze,  ea0, itag);
        wait_dma(imsk, MFC_TAG_UPDATE_ALL);

        wrkf(b,bsze, itag,otag);

        putf(b,bsze, ea0, otag);
	wait_dma(omsk, MFC_TAG_UPDATE_ALL);

        /* Done calculating & transferring wp */
        spu_write_out_mbox(WPDONE1);
    }

    return 0;
}









/* Create potential data by DMA */
static void start_create_pt(void* const pbuf, unsigned const szebuf,
                            env_t const* const env, pt_t* const p,
                            unsigned const tag)
{
     /* There are 4 pointers/EAs in pt_t/env_t
        Every pointer in env points to nelem items of type flt, so
        there are 4*nelem itmes to be DMAed in total
     */
     unsigned const nelem   = (env->ntypes)*(env->ntypes);
     /* Total number of bytes/allocation units needed */
     unsigned const nszetot = 4 * nelem * (sizeof (flt[4]));

     /* Enough memory available? */
     if ( EXPECT_TRUE_(szebuf>=nszetot) ) {
         /* Start DMAing the memory pointed to by pt.r2cut to LS, also getting data
            for lj_sig,... which follows after that in main mem. */

        mdma64(pbuf, env->r2cut, nszetot,  tag, MFC_GET_CMD);

        /* Set the pointers in the resulting pt_t (see modified allocation
           scheme in mk_pt)  & copy scalar member while DMA is running
	*/
        p->r2cut    = ((vector float *)pbuf);
        p->lj_sig   = p->r2cut  + nelem;
        p->lj_eps   = p->lj_sig + nelem;
        p->lj_shift = p->lj_eps + nelem;
        p->ntypes   = env->ntypes;
     }
     else {
        fprintf(stderr, "Could not allocate %u bytes for env.\n", nszetot);
        p->r2cut    = 0;
        p->lj_sig   = 0;
        p->lj_eps   = 0;
        p->lj_shift = 0;
     }
}




/* Start initialization (pt only, at the moment) */
static void start_init(unsigned const* const envea, unsigned const itag)
{
    /* Buffer for env. data */
    static unsigned char ALIGNED_(16, envdata[300]);

    /* env control block  */
    envbuf_t ALIGNED_(16, envbuf);

    /* 1st, DMA the control block containing the pointers */
    mdma64(&envbuf,  envea, (sizeof envbuf),  itag, MFC_GET_CMD);
    (void)wait_dma((1u<<itag), MFC_TAG_UPDATE_ALL);

    /* Start DMA of pt */
    start_create_pt(envdata, (sizeof envdata), ((env_t*)&envbuf), &pt,  itag);
}











/* Make EA list in eas (npairs elements),
   with initial value eastart
*/
static void arg_eas(unsigned const* const eastart,
                    unsigned long long const stride64,
                    unsigned* eas, unsigned const npairs)
{
    /* Number of pairs remaining */
    register unsigned k=npairs;
  
    /* Treat 32-bit pairs as 64bit words */
    typedef unsigned long long T64;
    register T64 cur64 = *((T64 const*)eastart);
    register T64* p64 = (T64*)eas;

    /**/
    for (  ;      (k>0);     --k, ++p64, cur64+=stride64  ) {
        (*p64)=cur64;
    }
}






/* Type of effective address of argument */
typedef union argea {
    unsigned ea32[2];
    ui64_t   ea64;
} argea_t;


/* argp "points" to an array of lengt N_ARGBUF of argbuf_t */
int main(ui64_t const id, argea_t const argp, argea_t const envp)
{
    /* (Working) buffer memory */
    static argbuf_t ALIGNED_(128, ibuf0);
    static argbuf_t ALIGNED_(128, ibuf1);
    static unsigned char ALIGNED_(128, obuf0[DIRECT_OBUFSZE]);
    static unsigned char ALIGNED_(128, obuf1[DIRECT_OBUFSZE]);

    /* Buffer effective addresses: An even number (at least N_ARGBUF) of pairs
       of unsigneds */
    unsigned int ALIGNED_(8, eabuf[N_ARGBUF<<1u]);



    /* Set decrementer to max. value for timing */
    spu_writech(SPU_WrDec, (~(0u)));




    /* Fetch env and start creating the corresponding pt_t (Using DMA) */
    start_init(envp.ea32, 0);

    /* Set EAs of the buffers in main memory */
    arg_eas(argp.ea32,  (sizeof(argbuf_t)), eabuf,N_ARGBUF);

    /* Wait for init. DMA */
    (void)wait_dma((1u<<0u), MFC_TAG_UPDATE_ALL);

    /* fprintf(stdout, "Starting work loop on SPU\n"); fflush(stdout); */

    /*
    workloop(argp.ea32,  &buf0, (sizeof buf0),
#if defined(CBE_DIRECT)
             init_get_direct, wrk_direct, init_put_direct
#else
             init_get,        wrk,        init_put
#endif
            );
    */


    for (;;) {
        /* Read mailbox */
        unsigned const offs = spu_read_in_mbox();
        unsigned const narg = spu_read_in_mbox();

        /* Use the EAs staring at offs*2 for input/output */
        unsigned const* const eas = eabuf + (offs<<1u);

        /* fprintf(stdout, "SPU %u: %u arguments at offset %u\n", ((unsigned)id), narg, offs); fflush(stdout); */

        workloop_ea(eas,eas, narg,
                    &ibuf0, &ibuf1,    &obuf0, &obuf1,
#if defined(CBE_DIRECT)
             init_get_direct, wrk_direct, init_put_direct
#else
             init_get,        wrk,        init_put
#endif
                   );

	spu_write_out_mbox(offs);
    }



    /* We should not arrive here as we only leave the main program if
       SPUEXIT is passed via mbox */
    return 1;
}








/* SPU (local) copy of pt which is initialized (once) via DMA */
pt_t pt;
