/* ISO C std. headers */
/* #include <stddef.h> */
#include <stdio.h>



/* SPU headers */
#include <spu_mfcio.h>
#include <spu_internals.h>

/* Local/program specific headers */
#include "imd_cbe.h"













#if ! defined(CBE_DIRECT)  /* Only needed in "indirect case" */
/* Given an exch_t, build the corresponding wp_t, that is:
   Copy the scalar members
   Allocate memory for all the array members (using memory pointed to by a).
   Fetch all the input array data to LS using DMA (possibly a list DMA) with
   the specified tag
*/
static void start_create_wp(exch_t const* const exch,
                            void* const tmpbuf, unsigned const tmplen,
                            void* const resbuf, unsigned const reslen,
                            unsigned const tag,
                            wp_t* const wp,
                            unsigned* pres_ea, unsigned* pres_sze)
{
     /* Number of bytes needed for the pos/force vectors */
     unsigned const nvecbytes = (exch->n2) * sizeof(flt[4]);

     /* Number of bytes remaining */
     unsigned       nrem;
     unsigned char* apos;

     /* Boundaries used for rounding up */
     register vector unsigned const bndv = { 15,127,  0,0 };
     
     /* Vector of sizes:
        szev[0]=to be transferd, szev[1]=to be allocated (remaing components
        are unused). Access the scalars using the following macros
     */
     vector unsigned szev;

#define ADV(v)  spu_extract((v),1)
#define XFR(v)  spu_extract((v),0)
 
     /* Number of bytes to advance/allocate in LS */
     register unsigned nadv;



     /* First allocate temp. mem. */
     nrem=tmplen;
     apos=tmpbuf;

     /* pos */
     nadv=ADV(szev=CEILV(nvecbytes, bndv));
     if ( EXPECT_TRUE(nrem>=nadv)  ) {
         mdma64(apos,  exch->pos, XFR(szev),  tag, MFC_GET_CMD);
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
     if ( EXPECT_TRUE(nrem>=nadv) ) {
         mdma64(apos,  exch->typ, XFR(szev),  tag, MFC_GET_CMD);
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
     if ( EXPECT_TRUE(nrem>=nadv) ) {
         mdma64(apos,  exch->ti, XFR(szev),  tag, MFC_GET_CMD);
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
     if ( EXPECT_TRUE(nrem>=nadv) ) {
         mdma64(apos, exch->tb,  XFR(szev),   tag, MFC_GET_CMD);
         wp->tb=(void*)apos;
         nrem-=nadv;
         apos+=nadv;
     }
     else {
         fprintf(stderr, "Could not allocate %u bytes for tb\n", nadv);
         wp->tb = 0;
     }



     /* Now, allocate result buffer */
     nrem = reslen;
     apos = resbuf;

     /* force (no DMA needed) */
     nadv=ADV(szev=CEILV(nvecbytes, bndv));
     if ( EXPECT_TRUE(nrem>=nadv) ) {
          wp->force = (void*)apos;
          nrem -= nadv;
          apos += nadv;

          /* Output EA and size for force */
          pres_ea[0] = exch->force[0];
          pres_ea[1] = exch->force[1];
          pres_sze[0] = XFR(szev);
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
     wp->n2_max  = exch->n2;
     wp->len     = exch->len;
     wp->len_max = exch->len_max;
}




/* DMA forces back to main memory */
static void start_DMA_results(unsigned const* const exch_ea,
                              unsigned const exch_out_sze,
                              wp_t const* const wp, exch_t* const exch,
                              unsigned const force_sze,
                              unsigned const tag
                             )
{
     /* First DMA the force array directly back to main memory... */
     /* DMA64(wp->force, res_sze[0],  exch->force,  tag, MFC_PUT_CMD); */
     mdma64(wp->force, exch->force,  force_sze, tag, MFC_PUT_CMD);

     /* In the meantime:
        Copy the scalars which have been updated in wp back to exch
        Leaving all the other scalar members unchanged
     */
     exch->totpot = wp->totpot;
     exch->virial = wp->virial;

     /* ...then DMA the updated exch_t */
     /* DMA64(exch, (sizeof (*exch)), exch_ea, tag, MFC_PUT_CMD); */
     mdma64(exch, exch_ea,  exch_out_sze,  tag, MFC_PUT_CMD);
}


#endif /* CBE_DIRECT */






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
     if ( EXPECT_TRUE(szebuf>=nszetot) ) {
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
static void start_init(unsigned const* const envea,
                       unsigned const tag, unsigned const cbtag)
{
    /* Buffer for env. data */
    static unsigned char ALIGNED_(16, envdata[300]);

    /* env control block  */
    envbuf_t ALIGNED_(16, envbuf);

    /* 1st, DMA the control block containing the pointers */
    mdma64(&envbuf,  envea, (sizeof envbuf),  cbtag, MFC_GET_CMD);
    spu_writech(MFC_WrTagMask, (1u<<cbtag));
    spu_mfcstat(MFC_TAG_UPDATE_ALL);

    /* Start DMA of pt */
    start_create_pt(envdata, (sizeof envdata), ((env_t*)envbuf), &pt,  tag);
}









/* First/past the end byte position of member m in structur object s */
#define MEMB_FRST_POS(s,m) (((unsigned char const*)&((s).m)) - ((unsigned char const*)&(s)))
#define MEMB_LAST_POS(s,m) (MEMB_FRST_POS((s),m) + (sizeof ((s).m)))

/* Maximum */
#define MAX(a,b) ((a)>(b) ? (a) : (b))




/* Argument types for main function */
typedef unsigned long long ui64_t;

/* 64bit EA split into 2 32bit parts (hi: ea32[0],  lo: ea32[1])
   used for argument & environment "pointer"
 */
typedef union {
    ui64_t   ea64;
    unsigned ea32[2];
} arg_ea_t, env_ea_t;


/* Advance argument pointer a by  n bytes */
static INLINE_ arg_ea_t advance_arg_ea(arg_ea_t a, ui64_t const n) {
    a.ea64+=n;
    return a;
}



/* argp0 "points" to an array of exch_t */
int main(ui64_t const id, arg_ea_t const argp0, env_ea_t const envp)
{
    /* EA of 2nd buffer  */
    arg_ea_t const argp1 = advance_arg_ea(argp0, sizeof (argbuf_t));


    /* Some tags used in the main routine to DMA workpackages/control blocks
       Tags > 2 will be used in order not to interfer with DMAs initiated
       in calc_wp_direct
     */

    /* Input tags */
    enum { itag0=3u,    itag1=4u    };
    /* Output tags are the same as input tags */
    enum { otag0=itag0, otag1=itag1 };
    /* Control block tag */
    enum { cbtag=5u };

    /* The corresponding masks */
    enum { imsk0=1u<<itag0, imsk1=1u<<itag1,
           omsk0=1u<<otag0, omsk1=1u<<otag1,
           cbmsk=1u<<cbtag };

    /* Suffix (kilobyte) */
    enum { Ki=1024u };

    /* wp_t or exch_t are buffered here */
    argbuf_t ALIGNED_(128, argbuf0), ALIGNED_(128, argbuf1);

    /* Buffers for the results */
    unsigned char ALIGNED_(128,resbuf0[17*Ki]), ALIGNED_(128, resbuf1[1*Ki]);



    /* Number of argument bytes which must be DMAed in.. */
    unsigned const arg_in_sze = (sizeof argbuf0);
    /* ..and out */
    unsigned const arg_out_sze = uiceil16(
                                 MAX(MEMB_LAST_POS(*((exch_t*)argbuf0), totpot),
                                     MEMB_LAST_POS(*((exch_t*)argbuf0), virial)
				    )
                                 );



    /* Some type info */
    /* printf("Sizes on SPU:\n");  sizeinfo(printf);   fflush(stdout); */


    /* Set decrementer to max. value */
    spu_writech(SPU_WrDec, (~(0u)));


    /* Fetch env and start creating the corresponding pt_t (Using DMA) */
    start_init(envp.ea32, itag0, cbtag);


    /* fprintf(stdout, "Starting work loop on SPU\n"); fflush(stdout); */




    /* Work loop */
    for(;;) {
       /* Buffer memory for DMAs (allocated "dynamically" from the following
          array) initiated by calc_wp */
       unsigned char ALIGNED_(128, calc_temp[90*Ki]); 


#if ! defined(CBE_DIRECT)
        /* Extra wp's are needed as only exch_t are passed which are
           not usable directly */
        wp_t wp0, wp1;

        /* EAs & sizes of the results, which at the moment is 
          the force vector only */
        enum { Nres=1u };
        unsigned res_ea[Nres<<1u], res_sze[Nres];
#endif




        /* Wait for signal from PPU */
        if ( EXPECT_FALSE(SPUEXIT==spu_read_in_mbox()) ) {
  	    return 0;
        }



#if defined(CBE_DIRECT)
        /* DMA the wp_t directly */
        mdma64(argbuf0, argp0.ea32, arg_in_sze,  itag0, MFC_GET_CMD);
#else
        /* Fetch exch controll block via DMA & wait for it to complete */
        mdma64(argbuf0, argp0.ea32, arg_in_sze,  cbtag,  MFC_GET_CMD);
        spu_writech(MFC_WrTagMask,  cbmsk);
        spu_mfcstat(MFC_TAG_UPDATE_ALL);

        /* Now, given the exch, start DMAing the array memebers in wp */
        start_create_wp(((exch_t const*)argbuf0),
                        calc_temp, (sizeof calc_temp),
                        resbuf0,   (sizeof resbuf0),
                        itag0,
                        &wp0,  res_ea, res_sze);
#endif




        /* Wait for for inbound DMA */
        spu_writech(MFC_WrTagMask,  imsk0);
        spu_mfcstat(MFC_TAG_UPDATE_ALL);


        /* Debugging output */
        /* fprintf(stdout, "Got wp k=%d\n", wp0.k);  fflush(stdout); */

        /* The main work to be done. */
#if defined(CBE_DIRECT)
        calc_wp_direct(((wp_t*)argbuf0),
                       calc_temp,(sizeof calc_temp), resbuf0,(sizeof resbuf0),
                       otag0);
#else
        calc_wp(&wp0);
#endif


#if defined(CBE_DIRECT)
        /* Forecs are DMAed back in calc_wp in "direct case", so here, we just
           DMA back the controll block containing virial & energy */
        mdma64(argbuf0, argp0.ea32, arg_out_sze,  otag0, MFC_PUT_CMD);
#else
        /* DMA results back to main memory  */
        start_DMA_results(argp0.ea32, arg_out_sze,
                          &wp0, ((exch_t*)argbuf0),
                          res_sze[0],
                          otag0
                         );
#endif  /* CBE_DIRECT */


        /* Wait for outbound DMA (including the DMA of the forces initiated
           in calc_wp) */
        spu_writech(MFC_WrTagMask, omsk0);
        spu_mfcstat(MFC_TAG_UPDATE_ALL);


        /* Done calculating & transferring wp */
        spu_write_out_mbox(WPDONE1);
    }


    /* We should not arrive here as we only leave the main program if
       SPUEXIT is passed via mbox */
    return 1;
}




/* SPU (local) copy of pt which is initialized (once) via DMA */
pt_t pt;
