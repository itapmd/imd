/* ISO C std. headers */
/* #include <stddef.h> */
/* #include <stdio.h> */



/* SPU headers */
#include <spu_mfcio.h>

/* Local/program specific headers */
#include "imd_cbe.h"



/* Is often needed */
typedef vector float  vecflt;


/* SPU (local) copy of pt which is initialized (once) via DMA and
   is needed by calc_wp
*/
pt_t pt;


/* Argument types for main function */
typedef unsigned long long ui64_t;

typedef union {
    ui64_t   ea64;
    unsigned ea32[(sizeof (ui64_t)) / (sizeof (unsigned))];
} Targ, Tenv;













/* Componentwise round up x to next multiples of r
   Components of boundary vector r must be a power of 2 minus 1
 */

/* Macro (on the SPU, x may be a scalar) 
   (Replace spu_ -> vec_ to get AltiVec)
 */
#define UICEILV(x,r)  spu_andc(spu_add((r),(x)), (r))

static INLINE_ vector unsigned uiceilv(register vector unsigned const x,
                                       register vector unsigned const r)
{
    /* Just call the macro */
    return UICEILV(x,r);
}











/* Given an exch_t, build the corresponding wp_t, that is:
   Copy the scalar members
   Allocate memory for all the array members (using memory pointed to by a).
   Fetch all the input array data to LS using DMA (possibly a list DMA) with
   the specified tag
*/
static void start_create_wp(void* const pbuf, unsigned const sizebuf,
                            exch_t const* const exch, wp_t* const wp,  unsigned const tag,
                            unsigned* pres_ea, unsigned* pres_sze)
{
     /* Number of bytes needed for the vectors */
     unsigned const nvecbytes = exch->n2*sizeof(vecflt);

     /* Number of bytes remaining */
     unsigned nrem = sizebuf;
     char*    apos = pbuf;

     /* Boundaries used for rounding up */
     register vector unsigned const bndv = { 15,127, 1023,4095 };
     
     /* Vector of sizes:
        szev[0]=to be transferd, szev[1]=to be allocated (remaing components
        are unused). Access the scalars using the following macros
     */
     vector unsigned szev;

#define ADV(v)  spu_extract((v),1)
#define XFR(v)  spu_extract((v),0)
 
     /* Number of bytes to advance in LS */
     register unsigned nadv;




     /* pos */
     nadv=ADV(szev=UICEILV(nvecbytes, bndv));
     if ( EXPECT_TRUE(nrem>=nadv)  ) {
         spu_mfcdma64((wp->pos=(void*)apos),  (exch->pos)[0], (exch->pos)[1],
                      XFR(szev),  tag, MFC_GET_CMD
		     );
         nrem-=nadv;
         apos+=nadv;
     }
     else {
         wp->pos = 0;
     }

     /* typ */
     nadv=ADV(szev=UICEILV(exch->n2*sizeof(int),  bndv));
     if ( EXPECT_TRUE(nrem>=nadv) ) {
         spu_mfcdma64((wp->typ=(void*)apos),  (exch->typ)[0], (exch->typ)[1],
                      XFR(szev),  tag, MFC_GET_CMD
                     );
         nrem-=nadv;
         apos+=nadv;
     }
     else {
         wp->typ = 0;
     }

     /* ti */
     nadv = ADV(szev=UICEILV(2*exch->n2*sizeof(int), bndv));
     if ( EXPECT_TRUE(nrem>=nadv) ) {
         spu_mfcdma64((wp->ti=(void*)apos),  (exch->ti)[0], (exch->ti)[1],
                      XFR(szev),   tag, MFC_GET_CMD
                     );
         nrem-=nadv;
         apos+=nadv;
     }
     else {
         wp->ti = 0;
     }

     /* tb */
     nadv=ADV(szev=UICEILV(exch->len*sizeof(short), bndv));
     if ( EXPECT_TRUE(nrem>=nadv) ) {
         spu_mfcdma64((wp->tb=(void*)apos), (exch->tb)[0], (exch->tb)[1],
                       XFR(szev),   tag, MFC_GET_CMD
                      );
         nrem-=nadv;
         apos+=nadv;
     }
     else {
         wp->tb = 0;
     }




     /* force (no DMA needed) */
     nadv=ADV(szev=UICEILV(nvecbytes, bndv));
     if ( EXPECT_TRUE(nrem>=nadv) ) {
          wp->force = ((void*)apos);
          nrem -= nadv;
          apos += nadv;

          /* Output EA and size for force */
          pres_ea[0]  = exch->force[0];
          pres_ea[1]  = exch->force[1];
          pres_sze[0] = XFR(szev);
     }
     else {
         pres_ea[0] = pres_ea[1]  =  pres_sze[0]=0;
         wp->force=0;
     }


#undef XFR
#undef ADV



     /* Now copy the scalar members from exch to wp */
     wp->k       = exch->k;
     wp->n1      = exch->n1;
     wp->n1_max  = exch->n1_max;
     wp->n2      = exch->n2;
     wp->n2_max  = exch->n2;
     wp->len     = exch->len;
     wp->len_max = exch->len_max;
}





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
     unsigned const nszetot = 4 * nelem * (sizeof (vecflt));

     /* Enough memory available? */
     if ( EXPECT_TRUE(szebuf>=nszetot) ) {
         /* Start DMAing the memory pointed to by pt.r2cut to LS, also getting data
            for lj_sig,... which follows after that in main mem. */
 
        /* DMA64((*a), nszetot, env->r2cut,  tag, MFC_GET_CMD); */
        spu_mfcdma64(pbuf, (env->r2cut)[0], (env->r2cut)[1], nszetot,
                     tag, MFC_GET_CMD);

        /* Set the pointers in the resulting pt_t (see modified allocation scheme
           in mk_pt)
	*/
        p->r2cut    = ((vecflt *)pbuf);
        p->lj_sig   = p->r2cut  + nelem;
        p->lj_eps   = p->lj_sig + nelem;
        p->lj_shift = p->lj_eps + nelem;

	/* Copy the scalar member */
        p->ntypes   = env->ntypes;
     }
}


static void start_init(void* const pbuf, unsigned const szebuf,
                       unsigned const* const envea,
                       pt_t* const ppt, unsigned const tag)
{
    /* Tag used for control block */
    enum { cbtag=5u };
    /* First DMA env which is only neede here */
    env_t ALIGNED_(16, env);


    spu_mfcdma64(&env,  envea[0],envea[1], (sizeof env),  cbtag, MFC_GET_CMD);
    spu_writech(MFC_WrTagMask, (1u<<cbtag));
    spu_mfcstat(MFC_TAG_UPDATE_ALL);

    /* Start DMA of pt */
    start_create_pt(pbuf, szebuf,  &env, ppt,  tag);
}



/* DMA forces back to main memory */
static void start_DMA_results(unsigned const* const exch_ea, unsigned const tag,
                              wp_t const* const wp, exch_t* const exch,
                              unsigned const* const res_sze
                             )
{
     /* First DMA the force array directly back to main memory... */
     /* DMA64(wp->force, res_sze[0],  exch->force,  tag, MFC_PUT_CMD); */
     spu_mfcdma64(wp->force, (exch->force)[0], (exch->force)[1],  res_sze[0],
                  tag, MFC_PUT_CMD
                 );

     /* In the meantime:
        Copy the scalars which have been updated in wp back to exch
        Leaving all the other scalar members unchanged
     */
     exch->totpot = wp->totpot;
     exch->virial = wp->virial;

     /* ...then DMA the updated exch_t */
     /* DMA64(exch, (sizeof (*exch)), exch_ea, tag, MFC_PUT_CMD); */
     spu_mfcdma64(exch, exch_ea[0], exch_ea[1],  (sizeof (*exch)),
                  tag, MFC_PUT_CMD
                 );
}








int main(ui64_t const id, Targ const argp, Tenv const envp)
{
    /* Unsigned integer constant, e.g. needed for tags/masks */
    typedef unsigned const Tuc;

    /* Some tags and the correspondning masks */
    enum { itag0=1u, itag1=2u,
           otag0=3u, otag1=4u,
           cbtag=5u };
    enum { imsk0=1u<<itag0, imsk1=1u<<itag1,
           omsk0=1u<<otag0, omsk1=1u<<otag1,
           cbmsk=1u<<cbtag };

 
    /* Local memory for DMAs (allocated "dynamically" on the stack) */
    unsigned char ALIGNED_(128, dmabuf0[50u*1024u]), 
               /* ALIGNED_(128, dmabuf1[40u*1024u]), */
                  ALIGNED_( 16, envbuf[512]);
 

    /* Set decrementer to max. value */
    spu_writech(SPU_WrDec, (~(0u)));


    /* Fetch env and start creating the corresponding pt_t (Using DMA) */
    start_init(&envbuf, (sizeof envbuf), envp.ea32,  &pt, itag0);


    /* fprintf(stdout, "Starting work loop on SPU\n"); fflush(stdout); */

    /* Work loop */
    for(;;) {
        /* Local copies of work packages on the SPU & the exch control block */
        wp_t wp0, wp1;
        exch_t ALIGNED_(16, exch0), ALIGNED_(16, exch1);

        /* EAs & sizes of the results, which at the moment is 
          the force vector only */
        enum { Nres=1u };
        unsigned res_ea[Nres<<1u], res_sze[Nres];



        /* Wait for signal from PPU */
        if ( EXPECT_FALSE(SPUEXIT==spu_read_in_mbox()) ) {
  	    return 0;
        }

        /* Fetch exch controll block: Start DMA & wait for it to complete */
        /* DMA64(&(exch[0]), (sizeof exch[0]), argp.ea32, 5u, MFC_GET_CMD); */
        spu_mfcdma64(&exch0,  (argp.ea32)[0],(argp.ea32)[1], (sizeof exch0),
                     cbtag,  MFC_GET_CMD);
        spu_writech(MFC_WrTagMask,  cbmsk);
        spu_mfcstat(MFC_TAG_UPDATE_ALL);


        /* fprintf(stdout, "Got control block\n"); fflush(stdout); */



        /* Start the creation (via DMA) of a wp_t using the exch_t  & 
           wait for it to finish 
	*/
        start_create_wp(dmabuf0, (sizeof dmabuf0), &exch0, &wp0,
                        itag0,  res_ea,res_sze);
        spu_writech(MFC_WrTagMask,  imsk0);
        spu_mfcstat(MFC_TAG_UPDATE_ALL);


        /* The main work to be done */
        calc_wp(&wp0);


        /* DMA results back to main memory & wait for it to complete */
        start_DMA_results(argp.ea32, otag0, &wp0, &exch0, res_sze);
        spu_writech(MFC_WrTagMask, omsk0);
        spu_mfcstat(MFC_TAG_UPDATE_ALL);


        /* Done calculating & transferring wp */
        spu_write_out_mbox(WPDONE1);
    }

    /* We should not arrive here as we only leave the main program if
       SPUEXIT is passed via mbox */
    return 1;
}
