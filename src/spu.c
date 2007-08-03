/* ISO C std. headers */
#include <stdio.h>

#if ! defined __SPU__
#include <stddef.h>
#endif

/* SPU headers */
#include <spu_mfcio.h>

/* Local/program specific headers */
#include "imd_cbe.h"

/* SPU (local) copy of pt which is initialized via DMA and is 
   needed by calc_wp
*/
pt_t pt;

/* Type of effective address used */
typedef unsigned long long  Tui64;
typedef ea_t                Tea;

/* Number of allocation units of size ausze needed for a block of size sze */
unsigned nalloc_block(unsigned const sze, unsigned const ausze)
{   
    /* Number of units */
    unsigned const nau = sze/ausze;

    /* Have to add one more alloction unit? */
    return (((nau*ausze) == sze) ?  nau : (1u+nau));
}

/* Number of allocation units of size "ausze" needed for 
   nelem items of size itmsze
*/
unsigned nalloc_items(unsigned const nelem, unsigned const itmsze,
                      unsigned const ausze)
{
     return nalloc_block(itmsze*nelem, ausze);
}

/* DMA sze bytes between LS addr. p and effective address ea in main mem */
void (DMA64)(void* const p, unsigned const sze,
             unsigned const* const ea,
             unsigned const tag, unsigned const cmd)
{
     spu_mfcdma64(p, ea[0], ea[1],  sze, tag, cmd);
}

/* List DMA with a list as returned by eas2les */
void (LDMA64)(void* const p, unsigned const* const list, unsigned const Nelem,
              unsigned const tag, unsigned const cmd)
{
    /* Size (in bytes) of one list element */
    enum { elemsze = 2u*(sizeof (unsigned)) };

    /*Call LDMA64 defined above with common high part of EAs taken
      from the usual position in the array
    */
    spu_mfcdma64(p, list[Nelem<<1u], ((unsigned)(list)), (Nelem*elemsze),
                 tag,cmd);
}


/* Allocate n items of size itmsze copying EAs
   The number of bytes actually allocated (which is a multiple of sizeof(Tau))
   is written to *isze.
*/
void* (alloc)(unsigned const n,  unsigned const itmsze,
              void** a, unsigned* const nrem,  unsigned const ausze,
              unsigned const* ieasrc,  unsigned** ieadst,  unsigned** isze)
{
     /* Number of allocation units/bytes needed */
     unsigned const nalloc = nalloc_items(n,itmsze, ausze);
     unsigned const nbytes = nalloc*ausze;

     /* Enough memory? */
     if ( EXPECT_TRUE((*nrem)>=nalloc) ) {
         /* Allocation starts here */
         void* const allocstart = *a;

         /* Copy EA, updating the iterators */
          *(*ieadst) = *( ieasrc);
         ++(*ieadst);
         ++( ieasrc);
          *(*ieadst) = *( ieasrc);
         ++(*ieadst);

         /* Copy size (in bytes), updating iterator */
         *(*isze) = nbytes;
         ++(*isze);

         /* Update allocation pointer/number of units remaining */
         *a     = ((unsigned char*)(*a)) + nbytes;
         *nrem -= nalloc;

         /* Return ptr. to beginning of the block just allocated */
         return allocstart;
     }

     /* Error */
     return 0;
}

/* Given an exch_t build the corresponding wp_t, that is:
   Copy the scalar members
   Allocate memory for all the array members (using memory pointed to by a).
   Fetch all the input array data to LS using DMA (possibly a list DMA) with
   the specified tag
*/
void start_create_wp(void** a, unsigned* nrem, unsigned const ausze,
                     exch_t const* const exch, wp_t* const wp,  unsigned const tag,
                     unsigned* pres_ea, unsigned* pres_sze)
{
     /* We need 4 EAs/sizes for the DMA + one high part of the EAs
        as well as some iterators over the EA/size arrays
     */
     enum {  Nea=4u };
     unsigned ALIGNED_(8, ea[(Nea<<1u)+1u]), sze[Nea];
     unsigned *pea=ea, *psze=sze;


     /* Get pointers to LS memory allocated from a, writing the size (in
        bytes to psze), copying the corresponding effective addresses to
        ea buffer
     */
     void* const astart=(*a);
#define ALLOCMEMB(mb,ty,nele)                                               \
   if ( EXPECT_FALSE(0 == ((wp->mb)=(ty*)alloc((nele), (sizeof (ty)), a,nrem, ausze, (exch->mb), &pea, &psze))) ) {   \
        return;                                                             \
   }
     ALLOCMEMB(pos, flt,   4*exch->n2)
     ALLOCMEMB(typ, int,     exch->n2)
     ALLOCMEMB(ti,  int,   2*exch->n2)
     ALLOCMEMB(tb,  short,   exch->len)
#undef ALLOCMEMB


     /* Start DMA to fetch the arrays */
     if ( EXPECT_TRUE(0 != eas2les(ea,sze,Nea)) ) {
         /* List DMA, using the list just built */
         LDMA64(astart, ea,Nea, tag, MFC_GETL_CMD);
     }
     else {
         /* 4 single DMAs for the 4 members */
#define DMAMEMB(member,nn) DMA64(wp->member,sze[nn], exch->member, tag, MFC_GET_CMD)
         DMAMEMB(pos, 0);
         DMAMEMB(typ, 1);
         DMAMEMB(ti,  2);
         DMAMEMB(tb,  3);
#undef DMAMEMB
     }
     


     /* Also allocate LS memory for force, copying its ea and alloction size */
     wp->force = (flt*)(alloc(4*exch->n2, (sizeof (flt)),
                              a,nrem,ausze, exch->force, &pres_ea, &pres_sze
                             )
                       );


     /* Now copy the scalar members from exch to wp */
#define CPYMEMB(member) ((wp->member)=(exch->member))
     CPYMEMB(k);
     CPYMEMB(n1);
     CPYMEMB(n1_max);
     CPYMEMB(n2);
     CPYMEMB(n2_max);
     CPYMEMB(len);
     CPYMEMB(len_max);

     /* No need to copy totpot & virial as it is initialized in calc_wp
        CPYMEMB(totpot);
        CPYMEMB(virial);
     */
#undef CPYMEMB
}



void start_create_pt(void** a, unsigned* nrem, unsigned const ausze,
                     env_t const* const env, pt_t* const p,  unsigned const tag)
{
     /* There are 4 pointers/EAs in pt_t/env_t
        Every pointer in env points to nelem items of type flt, so
        there are 4*nelem itmes to be DMAed in total
     */
     unsigned const nelem   = 4*(env->ntypes)*(env->ntypes);
     /* Total number of bytes/allocation units needed */
     unsigned const nszetot = 4*nelem*(sizeof (flt));
     unsigned const nau = nalloc_block(nszetot, ausze);

     /* Enough memory available? */
     if ( EXPECT_TRUE((*nrem)>=nau) ) {
         /* Start DMAing the memory pointed to by pt.r2cut to LS, also getting data
            for lj_sig,... which follows after that in main mem. */
 
        DMA64((*a), nszetot, env->r2cut,  tag, MFC_GET_CMD);

        /* Set the pointers in the resulting pt_t (see modified allocation scheme
           in mk_pt)
	*/
        p->r2cut    = ((flt*)(*a));
        p->lj_sig   = p->r2cut  + nelem;
        p->lj_eps   = p->lj_sig + nelem;
        p->lj_shift = p->lj_eps + nelem;

	/* Copy the scalar member */
        p->ntypes   = env->ntypes;

        /* Update alloction pointer/remaing size */
        (*a)   = ((unsigned char*)(*a)) + (nau*ausze);
        *nrem -= nau;
     }
}


/* DMA forces back to main memory */
void start_DMA_results(unsigned const* const exch_ea, unsigned const tag,
                       wp_t const* const wp, exch_t* const exch,
                       unsigned const* const res_sze
                      )
{
     /* First DMA the force array directly back to main memory... */
     DMA64(wp->force, res_sze[0],  exch->force,  tag, MFC_PUT_CMD);

     /* In the meantime:
        Copy the scalars which have been updated in wp back to exch
        Leaving all the other scalar members unchanged
     */
     exch->totpot = wp->totpot;
     exch->virial = wp->virial;

     /* ...then DMA the updated exch_t */
     DMA64(exch, (sizeof (*exch)), exch_ea, tag, MFC_PUT_CMD);
}

/* Argument type for main function */
typedef union {
    unsigned long long ea64;
    unsigned           ea32[2];
} Targ, Tenv;


int main(Tui64 const id, Targ const argp, Tenv const envp)
{
    /* Unsigned integer constant, e.g. needed for tags/masks */
    typedef unsigned const Tuc;

    /* Some tags and masks */
    static Tuc intag[]  = {      1u,       2u  };
    static Tuc inmsk[]  = { (1u<<1u), (1u<<2u) };
    static Tuc outtag[] = {      3u,       4u  };
    static Tuc outmsk[] = { (1u<<3u), (1u<<4u) };


    /* Local copies of work packages on the SPU */
    wp_t wp[1];

    /* Controllblocks used for DMA */
    env_t  ALIGNED_(16, env);
    exch_t ALIGNED_(16, exch[1]);



    /* Local memory for DMAs */
    enum { allocsize=16u };
    unsigned char ALIGNED_(16, DMAmem[50u*1024u]);
    /* Number of alloction units remaining, current point of allocation */
    unsigned DMArem = ASIZE(DMAmem);
    void* DMAloc = AFRST(DMAmem);

    /* EAs & sizes of the results, which at the moment is 
       the force vector only */
    enum { Nres=1u };
    unsigned res_ea[Nres<<1u], res_sze[Nres];

    /* printf("SPU main started\n"); */

    /* Fetch env and start creating the corresponding pt_t */
    DMA64(&env, (sizeof env), envp.ea32,  5u, MFC_GET_CMD);
    spu_writech(MFC_WrTagMask, (1u<<5u));
    spu_mfcstat(MFC_TAG_UPDATE_ALL);
    start_create_pt(&DMAloc,&DMArem, allocsize,  &env,&pt,  intag[0]);

    /* printf("Starting work loop on SPU\n"); */

    /* Work loop */
    for(;;) {
        /* Fetch exch controll block: Start DMA & wait for it to complete */
        DMA64(&(exch[0]), (sizeof exch[0]), argp.ea32, 5u, MFC_GET_CMD);
        spu_writech(MFC_WrTagMask,  (1u<<5u));
        spu_mfcstat(MFC_TAG_UPDATE_ALL);

        /* Start the creation (via DMA) of a wp_t using the exch_t  & 
           wait for it to finish 
	*/
        start_create_wp(&DMAloc,&DMArem, allocsize, &(exch[0]),&(wp[0]),  intag[0],
                        res_ea,res_sze);
        spu_writech(MFC_WrTagMask,  inmsk[0]);
        spu_mfcstat(MFC_TAG_UPDATE_ALL);

        /* The main work to be done */
        calc_wp(&(wp[0]));

        /* DMA results back to main memory & wait for it to complete */
        start_DMA_results(argp.ea32, outtag[0], &(wp[0]), &(exch[0]), res_sze);
        spu_writech(MFC_WrTagMask, outmsk[0]);
        spu_mfcstat(MFC_TAG_UPDATE_ALL);


        /* At the moment, just stop working here */
        break;
    }

    return 0;
}
