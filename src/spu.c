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







/* First/past the end byte position of member m in structure like object s */
#define MEMB_FRST_POS(s,m) (((unsigned char const*)&((s).m)) - ((unsigned char const*)&(s)))
#define MEMB_LAST_POS(s,m) (MEMB_FRST_POS((s),m) + (sizeof ((s).m)))

/* Maximum */
#define MAX(a,b) ((a)>(b) ? (a) : (b))

static INLINE_ unsigned uimax(unsigned const a, unsigned const b) {
    return MAX(a,b);
}




#if defined(CBE_DIRECT)  /* "Direct" case */





/* Init get (always the same for wp & tb) */
static void init_get(void* ibuf, void* obuf, unsigned const* iea, unsigned const itag)
{
    /* Start DMA to input buffer */
    /* dma64(ibuf, iea,  (sizeof (argbuf_t)),  itag,  MFC_GET_CMD); */
    spu_mfcdma64(ibuf, *iea,*(iea+1), (sizeof (argbuf_t)), itag, MFC_GET_CMD);
}


/* Init put (also the same for wp & tb) */
static void init_put(void* obuf, unsigned const* oea,  unsigned otag)
{
    /* Results are now in xfer buffer, so transfer it (at least partially) */
    /* dma64(obuf, oea,  (sizeof (argbuf_t)),  otag, MFC_PUT_CMD); */
    spu_mfcdma64(obuf, *oea,*(oea+1), (sizeof (argbuf_t)), otag, MFC_PUT_CMD);
}



/* Init put of exactly noting */
static void init_put_none(void* obuf, unsigned const* ea, unsigned otag)
{}




/* Temp. DMA buffer */
static unsigned char ALIGNED_(128, direct_tmp[80*Ki]);


/* Size of ouput buffers */
enum { DIRECT_OBUFSZE=((unsigned)(30*Ki)) };

/* Number of bytes remaining in result buffer */
enum { RESSZE = (unsigned)(DIRECT_OBUFSZE-(sizeof (argbuf_t))) };



static void wp_direct_spusum(void* i, void* o, unsigned otag)
{
    /* Ptrs. to input/summatiuon workpackages */
    wp_t* const pwp  = (wp_t*)i;
    wp_t* const psum = (wp_t*)o;


    /* Calculate to forces */
    calc_wp_direct(pwp,  &direct_tmp,(sizeof direct_tmp), 
                         ((argbuf_t*)o)+1,RESSZE,
                   otag);


    /* Add up forces/virial in result buffer */
    psum->virial += pwp->virial;
    psum->totpot += pwp->totpot;
}

static void wp_direct(void* i, void* o, unsigned otag)
{
    /* Ptrs. to input/output workpackages */
    wp_t* const pwp  = (wp_t*)i;
    wp_t* const pout = (wp_t*)o;


    /* Calculate to forces */
    calc_wp_direct(pwp,  direct_tmp,(sizeof direct_tmp),
                         ((argbuf_t*)o)+1,RESSZE,
                   otag);

    /* Copy some variables (scalars) to output buffer */
    pout->k    = pwp->k;
    pout->flag = pwp->flag;

    pout->virial = pwp->virial;
    pout->totpot = pwp->totpot;
}


/* Vectorized version of memcpy:
   copy nvec vectors from dst to src. Callers must make sure that the
   objects pointed to be dst and src are aligned to a 16-byte boundary.
 */
static void vmemcpy(void* const dst, void const* const src, unsigned const nvec)
{
    /* Vector type */
    typedef __vector unsigned char  Tvec;
    /* Iterators/number of vectors remaining */
    register Tvec const* isrc = ((Tvec const*)src);
    register Tvec      * idst = ((Tvec      *)dst);
    register unsigned    nrem = nvec;

    /* Until all vectors have been copied */
    for (;;) {
        /* Leave loop if there are no more vectors remaining
           (which is quite unlikely) */
        if ( EXPECT_FALSE_(0u==nrem) ) {
  	    return;
        }
        /* Copy, update pointers and number of vectors */
        (*(idst++)) = (*(isrc++));
        --nrem;
    }
}




static void tb_direct(void* i, void* o, unsigned otag)
{
    /* Some sizes */
    enum { SVEC=(sizeof(vector unsigned char)), SWP=(sizeof (wp_t)) };
    /* Number of vectors needed for one wp_t */
    enum { NVEC=(unsigned)(((SWP+SVEC-1) & (~(SVEC-1))) / SVEC) };

    /* Ptrs. to input/output workpackages */
    wp_t* const pwpin  = (wp_t*)i;
    wp_t* const pwpout = (wp_t*)o;

    /* Calculate NBL */
    calc_tb_direct(pwpin,  &direct_tmp,(sizeof direct_tmp),
                           ((argbuf_t*)o)+1, RESSZE,
                   otag);

    /* Here, we have to copy the entire WP to the output buffer */
    /* (*pwpout) = (*pwpin); */
    vmemcpy(pwpout,pwpin, NVEC);    
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




/**/
enum { IBUFSZE=((unsigned)(70*Ki)), OBUFSZE=((unsigned)(16*Ki)) };

static void init_get(void* ibuf, void* obuf, 
                     unsigned const* iea,  unsigned const itag)
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




static void init_put(void* obuf,  unsigned const* oea, unsigned otag)
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








/* Double bufferd workloop for list of given eas */
/* After this function returns there might still be some
   outbound DMAs pending. The unsigned returned by this function is
   a bit mask of the tags of pending DMAs.
 */
static unsigned workloop_ea(
   /* List of n input & output EAs */
   unsigned const* const ineas, unsigned const* const outeas, unsigned const n,
   /* Input/output buffers */
   void* const ib0, void* const ib1, void* const ob0, void* const ob1,
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


    /* Initial DMA to buffer 0 (if there are WPs remaining) */
    if ( EXPECT_TRUE_(nrem>0) ) {
        getf(ib0,ob0,  iea0,  t0);
        --nrem;
    }
    else {
        /* Return right now, if there are no WPs at all.
           Also, there's not outbound DMA to wait for */
        return 0;
    }

    for(;;) {
        /* Init DMA to buffer 1, if there are WPs remaining */
        if ( EXPECT_TRUE_(nrem>0) ) {
           getf(ib1,ob1,  iea1, t1);
           --nrem;
        }

        /* Wait for input and ouput buffers 0 to be available */
        /* wait_dma(m0, MFC_TAG_UPDATE_ALL); */
        spu_writech(MFC_WrTagMask,   m0);
        spu_writech(MFC_WrTagUpdate, MFC_TAG_UPDATE_ALL);
        spu_readch(MFC_RdTagStat);

        wrkf(ib0,ob0,        t0);
        putf(ob0,      oea0, t0);

        /* All WPs done? */
        if ( EXPECT_FALSE_(++ndone == n) ) {
            break;
        }

        /* Next EAs for buffer 0 */
  	iea0+=4;
        oea0+=4;


        /* Again, init DMA to buffer 0, if there are WPs remaining
         */
        if ( EXPECT_TRUE_(nrem>0) ) {
  	    getf(ib0,ob0, iea0, t0);
            --nrem;
        }

        /* Wait for, process & start DMAing back buffer 1*/
        /* wait_dma(m1, MFC_TAG_UPDATE_ALL); */
        spu_writech(MFC_WrTagMask,   m1);
        spu_writech(MFC_WrTagUpdate, MFC_TAG_UPDATE_ALL);
        spu_readch(MFC_RdTagStat);

        wrkf(ib1,ob1,        t1);
        putf(ob1,     oea1,  t1);

        /* All WPs done? */
        if ( EXPECT_FALSE_(++ndone == n) ) {
  	   break;
        }

        /* Next EA for buffer 1 */
        iea1+=4;
        oea1+=4;

    }

    /* Return the mask client code should wait for */
    return m01;
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
    /* Buffer for env. data, must be static potential data will be
       accessed after pt has been created. */
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

    for (  ;       (k>0);      --k,  ++p64,  cur64+=stride64  ) {
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
    /* Effective addresses (pairs of unsigneds) of the arguments
       in main memory */
    unsigned int ALIGNED_(8, eabuf[N_ARGBUF*2]);



    /* Set decrementer to max. value for timing */
    spu_writech(SPU_WrDec, (~(0u)));


    /* Fetch env and start creating the corresponding pt_t (Using DMA) */
    start_init(envp.ea32, 0);

    /* Set EAs of the buffers in main memory */
    arg_eas(argp.ea32,  (sizeof(argbuf_t)), eabuf,N_ARGBUF);


    /* Wait for init. DMA */
    (void)wait_dma((1u<<0u), MFC_TAG_UPDATE_ALL);





    for (;;) {
        /* Read work mode / flag */
        unsigned const mode = spu_read_in_mbox();

	/* Output buffers */
        static unsigned char ALIGNED_(128, obuf0[DIRECT_OBUFSZE]);
        static unsigned char ALIGNED_(128, obuf1[DIRECT_OBUFSZE]);

        /* mode==0: No more work   */
        if ( EXPECT_FALSE_(0==mode) ) {
 	    return 0;
        }
        else {
 	   /* "Functors passed to workloop" */
 	   void (* getfunc)(void* i, void* o, unsigned const* iea, unsigned it);
           void (* putfunc)(void* o, unsigned const* oea, unsigned ot);
           void (* wrkfunc)(void* i, void* o, unsigned ot);

           /* Input values from mailbox */
           unsigned offs,narg;


           /* Mode dependent setup */
           switch (mode) {
	       case DOWP : {
		  getfunc=init_get;
                  putfunc=init_put_none;
                  wrkfunc=wp_direct_spusum;
                  /* Set accumulation variables */
                  ((wp_t*)&obuf0)->virial = ((wp_t*)&obuf1)->virial =
		  ((wp_t*)&obuf0)->totpot = ((wp_t*)&obuf1)->totpot = 0;
                  break;
               }
	       case DOTB : {
		  getfunc=init_get;
                  putfunc=init_put;
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

               /* Mask of pending DMAs */
               unsigned omsk;

               /* Input/output EAs */
               unsigned const* eas;

               /* Read #of WPs, if number is >0 also read offset */
               if ( EXPECT_FALSE_(0 == (narg=spu_read_in_mbox())) ) {
		  break;
               }
               offs=spu_read_in_mbox();


               /* Debugging output */
               /* fprintf(stdout, "%u WPs at offset %u in mode %u\n", narg,offs,mode); fflush(stdout); */

               /* Process the arguments */
               eas = eabuf + (offs<<1u);
               omsk = workloop_ea(eas,eas,narg,
                                  &ibuf0,&ibuf1, &obuf0,&obuf1,
                                  getfunc,wrkfunc,putfunc);


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
              enum { RESOFFS2 = (unsigned)(RESOFFS*2u) };

	      /* Total energy/virial in obuf0 */
	      ((wp_t*)&obuf0)->totpot += ((wp_t const*)&obuf1)->totpot;
              ((wp_t*)&obuf0)->virial += ((wp_t const*)&obuf1)->virial;

              /* DMA results back */
              dma64(&obuf0, eabuf+RESOFFS2,  sizeof(argbuf_t),
                    t, MFC_PUT_CMD);
              wait_dma(m, MFC_TAG_UPDATE_ALL);


              /* Tell PPU where to pick up accu */
              spu_write_out_mbox(RESOFFS);
           }
           if ( DOTB==mode ) {
	      enum { RESOFFS=(unsigned)1 };
	      /* Just acknoledge */
	      spu_write_out_mbox(RESOFFS);
           }
        }
    }





    /* We should not arrive here as we only leave the main program if
       SPUEXIT is passed via mbox */
    return 1;
}








/* SPU (local) copy of pt which is initialized (once) via DMA */
pt_t pt;
