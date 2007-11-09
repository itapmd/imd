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
static void init_get(void* ibuf, void* obuf, ea_t iea, unsigned const itag)
{
    /* Start DMA to input buffer usinga fence command, as the input buffer
       may have been used for a previous put command to system-memory */
    /* dma64(ibuf, iea,  (sizeof (argbuf_t)),  itag,  MFC_GETF_CMD); */
    dma64(ibuf, iea, (sizeof (argbuf_t)), itag, MFC_GETF_CMD);
}



/* Temp. DMA buffer */
static unsigned char ALIGNED_(128, direct_tmp[80*Ki]);


/* Size of ouput buffers */
enum { DIRECT_OBUFSZE=((unsigned)(30*Ki)) };



static void wp_direct_spusum(void* i, void* o,  ea_t unused, unsigned otag)
{
    /* Ptrs. to input/summatiuon workpackages */
    wp_t* const pwp  = (wp_t*)i;
    wp_t* const psum = (wp_t*)o;

    /* Size of buffer for the results */
    enum { RESSZE=DIRECT_OBUFSZE-(sizeof (argbuf_t)) };

    /* Calculate to forces */
    calc_wp_direct(pwp,  &direct_tmp,(sizeof direct_tmp), 
                   ((argbuf_t*)o)+1, RESSZE,
                   otag);


    /* Add up forces/virial in result buffer */
    psum->virial += pwp->virial;
    psum->totpot += pwp->totpot;
}

static void wp_direct(void* i, void* o, ea_t oea, unsigned otag)
{
    /* Calculate to forces */
    calc_wp_direct(((wp_t*)i),  direct_tmp,(sizeof direct_tmp),
                   o, DIRECT_OBUFSZE,
                   otag);


    /* Init outbound DMA: */
    dma64(i,  oea,  (sizeof (argbuf_t)),  otag, MFC_PUT_CMD);
}



/* Vectorized version of memcpy:
   copy nvec vectors from dst to src. Callers must make sure that the
   objects pointed to be dst and src are aligned to a 16-byte boundary.
 */
static void* vmemcpy(void* const dst, void const* const src, unsigned const nvec)
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
  	    break;
        }
        /* Copy, update pointers and number of vectors */
        (*(idst++)) = (*(isrc++));
        --nrem;
    }

    /* Retrun address of destination buffer */
    return dst;
}




static void tb_direct(void* i, void* o, ea_t oea, unsigned otag)
{
    /* Calculate NBL */
    calc_tb_direct(((wp_t*)i),  &direct_tmp,(sizeof direct_tmp),
                   o, DIRECT_OBUFSZE,
                   otag);

    /* We need a fence here */
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




/**/
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








/* Double bufferd workloop for list of given eas */
/* After this function returns there might still be some
   outbound DMAs pending. The unsigned returned by this function is
   a bit mask of the tags of pending DMAs.
 */
static unsigned workloop_ea_(
   /* List of n input & output EAs */
   ea_t const* const ineas, ea_t const* const outeas, unsigned const n,
   /* Input/output buffers */
   void* const* const ib,  void* const* const ob,
   /* Tags */
   unsigned const* const t,
   /* Functions to get, to process and to write back data */
   void (* const getf)(void* ibuf, void* obuf, ea_t const iea, unsigned it),
   void (* const wrkf)(void* ibuf, void* obuf, ea_t const oea, unsigned ot)
)
{
    /* Number of WPs remaining to be DMAed in, number of WPs done */
    register unsigned nget = n;
    register unsigned nwrk = n;

    /* Pointer to eff. addr. pairs */
    typedef ea_t const*  Pea;
    Pea iea[2], oea[2];

    /* Masks (qwords) */
    unsigned m[2];

    /* Set masks */
    m[0]=1u<<t[0];
    m[1]=1u<<t[1];

    /* Set input/output EAs */
    iea[1]=(iea[0]=ineas )+1;
    oea[1]=(oea[0]=outeas)+1;

    /* Initial DMA to buffer 0 (if there are WPs remaining) */
    if ( EXPECT_TRUE_(nget>0u) ) {
        getf(ib[0],ob[0],  *(iea[0]),t[0]);
        --nget;
    }
    else {
        /* Return right now, if there are no WPs at all.
           Also, there's no outbound DMA to wait for */
        return 0;
    }




    for(;;) {
        /* Init DMA to buffer 1, if there are WPs remaining */
        if ( EXPECT_TRUE_(nget>0u) ) {
           getf(ib[1],ob[1],  *(iea[1]),t[1]);
           --nget;
        }

        /* Wait for input and ouput buffers 0 to be available */
        /* wait_dma(m0, MFC_TAG_UPDATE_ALL); */
        spu_writech(MFC_WrTagMask,   m[0]);
        spu_writech(MFC_WrTagUpdate, MFC_TAG_UPDATE_ALL);
        spu_readch(MFC_RdTagStat);

        wrkf(ib[0],ob[0],    *(oea[0]), t[0]);

        /* All WPs done? */
        if ( EXPECT_FALSE_(0u == --nwrk) ) {
            break;
        }

        /* Next EAs for buffer 0 */
  	iea[0]+=2;
        oea[0]+=2;


        /* Again, init DMA to buffer 0, if there are WPs remaining
         */
        if ( EXPECT_TRUE_(nget>0u) ) {
  	    getf(ib[0],ob[0],  *(iea[0]),t[0]);
            --nget;
        }

        /* Wait for, process & start DMAing back buffer 1*/
        /* wait_dma(m1, MFC_TAG_UPDATE_ALL); */
        spu_writech(MFC_WrTagMask,   m[1]);
        spu_writech(MFC_WrTagUpdate, MFC_TAG_UPDATE_ALL);
        spu_readch(MFC_RdTagStat);

        wrkf(ib[1],ob[1],  *(oea[1]),t[1]);

        /* All WPs done? */
        if ( EXPECT_FALSE_(0u == --nwrk) ) {
  	   break;
        }

        /* Next EA for buffer 1 */
        iea[1]+=2;
        oea[1]+=2;

    }

    /* Return the mask client code should wait for */
    return (m[0]|m[1]);
}




/* C.f. Prog. Handbook p. 687  */

/* Double bufferd workloop for list of given eas */
/* After this function returns there might still be some
   outbound DMAs pending. The unsigned returned by this function is
   a bit mask of the tags of pending DMAs.
 */
static unsigned workloop_ea(
   /* List of n input & output EAs */
   ea_t const* const ieas, ea_t const* const oeas, unsigned const nea,
   /* Input/output buffers */
   void* const* const ib,  void* const* const ob,
   /* Tags */
   unsigned const* const t,
   /* Functions to get, to process and to write back data */
   void (* const getf)(void* ibuf, void* obuf, ea_t const iea, unsigned it),
   void (* const wrkf)(void* ibuf, void* obuf, ea_t const oea, unsigned ot)
)
{
    /* Masks */
    unsigned m[2];

    /* Ptrs. to EA pairs / buffers */
    typedef ea_t const* Pea;
    Pea iea[2], oea[2];


    /* Indices */
    register unsigned cur = 0u;
    register unsigned nxt = cur^1u;
    register unsigned n   = nea;


    /* Set masks */
    m[cur]=1u<<t[cur];
    m[nxt]=1u<<t[nxt];


    /* Set EAs / buffers */
    iea[nxt]=iea[cur]=ieas;
    oea[nxt]=oea[cur]=oeas;

    /* Get buffer 0 */
    getf(ib[cur],ob[cur],  *(iea[cur]),t[cur]);


    for (;;) {
        /* Are we done? */
        if ( EXPECT_FALSE_(0u == (--n)) ) {
  	    break;
        }

        /* Next index: cur=0,1 -> nxt=1,0 */
        nxt=cur^1u;   
        /* Next EA pairs */    
        iea[nxt] = iea[cur]+1;
        oea[nxt] = oea[cur]+1;
        /* Get next buffer */
        getf(ib[nxt],ob[nxt],   *(iea[nxt]), t[nxt]);

         
        /* Wait for current input/output buffers to become available */
        spu_writech(MFC_WrTagMask,   m[cur]);
        spu_writech(MFC_WrTagUpdate, MFC_TAG_UPDATE_ALL);
        (void)(spu_readch(MFC_RdTagStat));

        /* Work on current buffer */
        wrkf(ib[cur],ob[cur],     *(oea[cur]), t[cur]);

        /* Move on to next buffer */
        cur=nxt;
    }


    /* Wait for last buffer to become available */
    spu_writech(MFC_WrTagMask,   m[cur]);
    spu_writech(MFC_WrTagUpdate, MFC_TAG_UPDATE_ALL);
    (void)(spu_readch(MFC_RdTagStat));

    /* Work on buffer */
    wrkf(ib[cur],ob[cur],  *(oea[cur]), t[cur]);

    /* Return the mask client code should wait for */
    return m[cur];
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
static void start_init(ea_t const envea, unsigned const itag)
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
static void arg_eas(ea_t const eastart,
                    ea_t const stride64,
                    ea_t* eas, unsigned const neas)
{
    /* Number of pairs remaining */
    register unsigned k=neas;
  
    register ea_t  cur64 = eastart;
    register ea_t* p64 = eas;

    for (  ;       (k>0);      --k,  ++p64, cur64+=stride64  ) {
        (*p64)=cur64;
    }
}






/* Type of effective address of argument */
typedef union argea {
    unsigned ea32[2];
    ea_t     ea64;
} argea_t;





/* argp "points" to an array of lengt N_ARGBUF of argbuf_t */
int main(ui64_t const id, argea_t const argp, argea_t const envp)
{
    /* Effective addresses (pairs of unsigneds) of the arguments
       in main memory */
    ea_t ALIGNED_(8, eabuf[N_ARGBUF]);

    /* Tags (as vectors) */
    qword qwtag[2];

    /* Set decrementer to max. value for timing */
    spu_writech(SPU_WrDec, (~(0u)));


    /* Fetch env and start creating the corresponding pt_t (Using DMA) */
    start_init(envp.ea64, 0);

    /* Set EAs of the buffers in main memory */
    arg_eas(argp.ea64,  (sizeof(argbuf_t)), eabuf,N_ARGBUF);

    /* Create tag qwords */
    qwtag[0] = si_from_uint(0u);
    qwtag[1] = si_from_uint(1u);


    /* Wait for init. DMA */
    (void)wait_dma((1u<<0u), MFC_TAG_UPDATE_ALL);





    for (;;) {
        /* Read work mode / flag */
        unsigned const mode = spu_read_in_mbox();

	/* Output buffers */
        static unsigned char ALIGNED_(128, obuf0[DIRECT_OBUFSZE]);
        static unsigned char ALIGNED_(128, obuf1[DIRECT_OBUFSZE]);
        static void* const ob[] = { &obuf0, &obuf1 };


        /* mode==0: No more work   */
        if ( EXPECT_FALSE_(0==mode) ) {
 	    return 0;
        }
        else {
 	   /* "Functors passed to workloop" */
 	   void (* getfunc)(void* i, void* o, ea_t iea, unsigned it);
           void (* wrkfunc)(void* i, void* o, ea_t oea, unsigned ot);

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
               
               /* Tags to be used */
               static unsigned int const tags[] = { 0,1 };

               /* Mask of pending DMAs */
               unsigned omsk;

               /* Input/output EAs */
               ea_t const* eas;

               /* Read #of WPs, if number is >0 also read offset */
               if ( EXPECT_FALSE_(0 == (narg=spu_read_in_mbox())) ) {
		  break;
               }
               offs=spu_read_in_mbox();


               /* Debugging output */
               /* fprintf(stdout, "%u WPs at offset %u in mode %u\n", narg,offs,mode); fflush(stdout); */

               /* Process the arguments */
               eas = eabuf+offs;
               omsk = workloop_ea(eas,eas,narg, ib,ob, tags, getfunc,wrkfunc);


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
              dma64(&obuf0, eabuf[RESOFFS],  sizeof(argbuf_t),
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
