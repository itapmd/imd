
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
* imd_cbe_util.c -- utility functions for CBE
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

/* IMD config */
#include "config.h"

/* CBE specific stuff */
#include "imd_cbe.h"

/* Headers for SPU */
#if defined(__SPU__)

#include <spu_mfcio.h>

#endif

/* Headers for PPU */

#if defined(__PPU__) 

/* ISO C std. headers */
#include <limits.h>
#include <stddef.h>
/* Hosted */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* POSIX std. headers */
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <pthread.h>

/* CBE Headers */
#include <libspe2.h>

/* SDK Headers containing inline functions/macros */
/* libsync */
/* conditional variables */
#include <cond_init.h>
#include <cond_wait.h>
#include <cond_signal.h>
#include <cond_broadcast.h>
/* Mutexes */
#include <mutex_init.h>
#include <mutex_lock.h>
/* Atomic operations */
#include <atomic_inc_return.h>
#include <atomic_set.h>

#endif


/******************************************************************************
*
*  SPU specific part
*
******************************************************************************/

#if defined(__SPU__)

/* 16k are maximum which may be DMAed in one call to spu_mfc... */
enum { DMAMAX   = ((unsigned)(16*1024)) };
enum { DMAMAXea = ((ea_t)DMAMAX)        };


void mdma64_rec(register unsigned char* const p,
                register ea_t const lea, register unsigned const sze,
                register unsigned const tag, register unsigned const cmd)
{
    /* Only one DMA needed? End recusrion */
    if ( EXPECT_TRUE_(sze<=DMAMAX) ) {
        dma64(p, lea, sze,  tag,cmd);
        return;
    }

    /* Start a DMA */
    dma64(p, lea,  DMAMAX,  tag,cmd);
    /* Tail recursive call... */
    mdma64_rec(p+DMAMAX,  lea+DMAMAX,  sze-DMAMAX,  tag,cmd);
}


/* DMA more than 16k using multiple DMAs with the same tag */
void mdma64_iter(register unsigned char* p,
                 register ea_t lea, register unsigned remsze,
                 register unsigned const tag, register unsigned const cmd)
{
    for(;;) {
        /* We're done if not more than maxsze bytes are to be xfered */
        if ( EXPECT_TRUE_(remsze<=DMAMAX) ) {
  	   dma64(p, lea,  remsze, tag, cmd);
 	   return;
        }

        /* Xfer one chunk of maximum size  */
        dma64(p,  lea,  DMAMAX, tag, cmd);

        /* fputs(stdout,"More than one DMA needed in mdma64_iter\n"); */

        /* Update addresses and number of bytes */
        remsze -= DMAMAX;
        p      += DMAMAX;
        lea    += DMAMAX;
    }
}


/* Wait for DMA specified by mask m using type for update */
INLINE_ unsigned wait_dma(unsigned const m, unsigned const type) {
    spu_writech(MFC_WrTagMask,   m);
    spu_writech(MFC_WrTagUpdate, type);
    return spu_readch(MFC_RdTagStat);
}


/* (Just a) wrapper for DMA with less than 16K */
INLINE_ void (dma64)(register void* const p,
                            register ea_t const ea, register unsigned int const size,
                            register unsigned int const tag, register unsigned int const cmd
                  )
{
    si_wrch(MFC_EAH,   si_from_uint((unsigned int)(ea>>32u)));
    si_wrch(MFC_EAL,   si_from_uint((unsigned int)ea));
    si_wrch(MFC_LSA,   si_from_ptr(p));
    si_wrch(MFC_Size,  si_from_uint(size));
    si_wrch(MFC_TagID, si_from_uint(tag));
    si_wrch(MFC_Cmd,   si_from_uint(cmd));
}


/* DMA more than 16K using multiple DMAs */
INLINE_ void mdma64(void* const p, ea_t const ea, unsigned const size,
                           unsigned const tag, unsigned const cmd
                          )
{
    /* Byte type */
    typedef unsigned char    Tbyte;

    /* Defined somwhere in imd_cbe_util.c */
    extern void mdma64_rec(register Tbyte* const, register ea_t const, register unsigned const,
                           register unsigned const, register unsigned const);
    extern void mdma64_iter(register Tbyte*,  register ea_t, register unsigned,
                            register unsigned const, register unsigned const);



    mdma64_iter(((Tbyte*)p),  ea, size,  tag,cmd);
}

#endif /* SPU specific part */


/******************************************************************************
*
*  PPU specific part
*
******************************************************************************/

#if defined(__PPU__)

/******************************************************************************
*
*  Get time base from /proc/cpuinfo
*
******************************************************************************/

unsigned long tbfreq(void) {

  char *str, line[128];
  FILE *inp;
  unsigned long timebase = 0u;

  inp = fopen("/proc/cpuinfo", "r");
  if (NULL==inp) return 0u;
  while (!feof(inp)) {
    str=fgets(line,128,inp);
    if (NULL==str) break;
    if (strncmp(line, "timebase", 8) == 0) {
      str = strstr(line,":") + 1;
      sscanf(str, "%lu", &timebase);
    }
  }
  fclose(inp);
  return timebase;
}


/******************************************************************************
*
*  Memory management
*
******************************************************************************/

/* Return the page size */
size_t pagesize(void)
{
    return sysconf(_SC_PAGESIZE);
}


/* Helper function: Replace zero argument by page size */
static size_t repl0pgsze(size_t const n)
{
    return ((0==n) ? (pagesize()) : n);
}

/* Allocate at least s bytes of memory which are aligned ton an alig0
   boundary. The requested size is rounded up towards the next higher 
   multiple of mult0 */
void* (malloc_aligned)(size_t const s0, size_t const alig0, size_t const mult0)
{
    /* Use pagesize for mult if it is zero */
    size_t const mult = repl0pgsze(mult0);

    /* s0=qu*mult + some remainder which we do not need
       k*mult <= s0, but never k*mult>s0
     */
    size_t const qu = s0/mult;


    /* Try to allocate aligend memory returning the staring address.
       pagesize is used for the alignment parameter if alig0 is zero.
       Also, the amount of memory requested is rounded up to the next
       multiple of mult which is larger than s0.
     */
    void* res;
    return ((0 == posix_memalign(&res,
                                 (repl0pgsze(alig0)),
                                 ((qu*mult ==s0) ?  s0 : ((qu+1)*mult))
                                )
            )
               ? res : NULL
           );
}


/* The same as above allocation nelem elements of size elsze each */
void* (calloc_aligned)(size_t const nelem, size_t const elsze, 
                       size_t const alig0, size_t const mult0)
{
     return malloc_aligned((nelem*elsze), alig0, mult0);
}


/******************************************************************************
*
*  Various declarations
*
******************************************************************************/

/* Swaps intger variables a and b (of same type)*/
#define SWAP(a,b) { (a)^=(b); (b)^=(a); }
/* Sets x to zero without even mentioning zero :-) */
#define CLEAR(x)  ((x)^=(x))


/* Some local abreviations for two libspe2 data types
   which are often used here */
typedef spe_program_handle_t Thndl;
typedef spe_context_ptr_t    Tctxt;
/* The same for the pthread type */
typedef pthread_t            Tthr;


/* Arrays of buffers:  */

/* Environment */
/* cbe_env_begin is passed as environment address to every SPU thread */
static ea_t ALIGNED_(128, enveabuf[N_ENVEA]);
ea_t* const cbe_envea_begin = enveabuf;


/* Arguments */
static argbuf_t ALIGNED_(128, abuf[N_SPU_THREADS_MAX * N_ARGBUF]);
/* cbe_arg_begin + N_ARGBUF*k is passed as argument address to SPU thread k */
argbuf_t* const cbe_arg_begin = abuf;


/* Additional workpackage buffer */
static wp_t wbuf[N_SPU_THREADS_MAX * N_ARGBUF]; 
wp_t* const cbe_wp_begin = wbuf;


/* Stop info for the SPU threads */
static spe_stop_info_t sinfo[N_SPU_THREADS_MAX];
spe_stop_info_t const* const cbe_stopinfo = sinfo; 


/* SPU contexts.
   Those must be available to the client code as they will be needed for
   mailbox communication */
static spe_context_ptr_t c[N_SPU_THREADS_MAX];
spe_context_ptr_t const* const cbe_spucontext = c;



/* Control areas */
static spe_spu_control_area_p contrarea[N_SPU_THREADS_MAX];
spe_spu_control_area_p* const cbe_spucontrolarea = contrarea;


/* Reservation granule type large enough to hold locks/mutexes, etc...
   unsigned is used as base type as some libsync routines require the 
   adresses to be 32bit aligned.
 */
typedef unsigned int uirgran_t[CBE_GRANULE_NINT];

/* Mutexes/Condition variables used by libsync functions */
static uirgran_t ALIGNED_(16, lsmutx[128/(sizeof (uirgran_t))]);
static uirgran_t ALIGNED_(16, lscvar[128/(sizeof (uirgran_t))]);

/* Number of mutexes/cond. vars. */
enum {
   NMUTX = (unsigned)((sizeof lsmutx)/(sizeof (lsmutx[0]))),
   NCVAR = (unsigned)((sizeof lscvar)/(sizeof (lscvar[0])))
};


/* Minimum of m or the number of usable SPEs */
int min_usable_spes(int const m, int const cpu_node) {
    /* Get number of usable SPEs from OS */
    int const u = spe_cpu_info_get(SPE_COUNT_USABLE_SPES, cpu_node);

    if ( -1 == u ) {
        /* Error getting number */
        return -1;
    }

    /* Now return the minum */
    return ((m<u) ? m : u);
}


/* Argument struct for the pthreads which pass control to the SPU */
typedef struct spu_pthr_arg {
     /* Argument & environment passed to the SPU thread*/
     void *spu_arg;
     void *spu_env;

     /* Context to switch to */
     spe_context_ptr_t spu_ctxt;

     /* Entry point for SPU code */
     /* unsigned spu_entry; */

     /* Location of stop info */
     spe_stop_info_t* spu_stopinfo_loc;
} spu_pthr_arg_t;



/* Trivial pthread function which just switches to SPU context and
   returns a pointer to the SPU program's stop info. */
static void* spu_pthr(void* p0)
{
    /* Check the argument ptr. and use it in case it's valid */
    if ( p0 ) {
        /* SPU instruction counter, initialized to start of SPU program */
        unsigned spu_entry = SPE_DEFAULT_ENTRY;

        /* Get the (constant) argument pointed to by p0 */
        spu_pthr_arg_t const* const p = p0;


        /* The last message from PPU */
        /* fprintf(stdout, "About to switch to SPU context\n");  
           fflush(stdout); */

        if ( spe_context_run(p->spu_ctxt,  &spu_entry, 0u,
                             p->spu_arg, p->spu_env,
                             p->spu_stopinfo_loc
                            ) 
           ) {
    	      return p->spu_stopinfo_loc;
        }
    }

    /* Error */
    return 0;
}



/* Number of threads initialized */
static unsigned int nspus=0;

/* Get number of threads initialized */
unsigned cbe_get_nspus() {
    return nspus;
}


/* SPU multithreading:  IDs of the POSIX threads */
static Tthr  t[N_SPU_THREADS_MAX];


/* The "real" initialization routine performing the actual work */
int cbe_init0(int const nspu_req, int const cpu_node)
{
    /* The env - obsolete */
    /* static env_t ALIGNED_(16,env); */

    /* Arguments for the SPU threads */
    static spu_pthr_arg_t pthrargs[N_SPU_THREADS_MAX];

    /* Number of threads started */
    unsigned int nstrt=0;

    /* Number of SPUs which may be initialized */
    int const imax = min_usable_spes(((nspu_req<N_SPU_THREADS_MAX) ?
                                     nspu_req : N_SPU_THREADS_MAX), cpu_node
                                    );

    /* Some indices/iterators */
    int i;
    Tctxt* pctxt;
    spu_pthr_arg_t* parg;


    /* No SPUs available */
    if ( imax < 1 ) {
        nspus=0;
        return -1;
    }


    /* Set the environement EAs */
    cbe_envea_begin[0] = EA_CAST(&pt);
    cbe_envea_begin[1] = EA_CAST(lscvar);
    cbe_envea_begin[2] = EA_CAST(lsmutx);


    /* init condition variable(s) / mutexes */
    for ( i=0;  i<NCVAR;  ++i ) {
        _cond_init(EA_CAST(lscvar+i));
    }
    for ( i=0;  i<NMUTX;  ++i ) {
        _mutex_init(EA_CAST(lsmutx+i));
    }
    _mutex_lock(EA_CAST(lsmutx+0));

    /* Set "ground state" for some publically accessible values */
    for ( i=0;   i<N_SPU_THREADS_MAX ;  ++i  ) {
        /* context & control area */
        c[i]=0;
        contrarea[i]=0;
    }

    /* Sync memory before creating threads */
    /* __lwsync(); */

    for ( i=0, pctxt=c, parg=pthrargs;    (i<imax);   )  {
        /* Set the arguments for the pthread */ 
        /* SPU thread arguments, also set wp fields to initial value
           such that the wp gets allocated upon first call to make_wp */
        parg->spu_arg = cbe_arg_begin + i*N_ARGBUF;

        /* Environment */
        parg->spu_env = cbe_envea_begin;

        /* Stop info of SPU program is written here */
        parg->spu_stopinfo_loc=sinfo+i;

        /* Create context & load program into it */
        if ( 0 != ((*pctxt) = (parg->spu_ctxt) = 
               spe_context_create(SPE_MAP_PS,0)) ) {
   	    /* This handle must be defined somewhere else and the program
               must contain the calc_wp code. */
            extern spe_program_handle_t hndle_cbe_calc;

            /* Get control area for the context just created */
            contrarea[i] = spe_ps_area_get((*pctxt), SPE_CONTROL_AREA);
            /*  fprintf(stdout, "Control area %i at %p\n", i, contrarea[i]);
                fflush(stdout);  */

            if ( 0 == spe_program_load((*pctxt), &hndle_cbe_calc)  ) {
                /* Now start the a pthread */
                if ( 0 == pthread_create((t+i), NULL, &spu_pthr, parg)  ) {
		    /* (Another) thread started successfully */
  		    ++nstrt;

                    /* Move on to next thread */
                    ++i;
                    ++parg;
                    ++pctxt;

                    /* fprintf(stdout, "pthread %u started!\n", nstrt); */
                }
                else {
		    /* Could not start thread, so destroy the (unused!)
                       SPU context */
		    spe_context_destroy(*pctxt);
                    (*pctxt)=0;
                    contrarea[i]=0;
                    fprintf(stderr, "Could not create pthread\n");
                }
   	    } else {
   	         fprintf(stderr, "Could not load handle\n");
                 return(-1);
            }
        }
        else {
      	   fprintf(stderr, "Could not create context\n");
        }
    }


    /* Store the number of threads actually started */
    nspus=nstrt;

    /* OK, CBE initialized */
    return 0;
}


/* Init stuff for calc_threads. Try to start SPU threads
   This function may only be called once. It's not thread safe
   The minimum of nspu_req and the number of SPUs available (according
   cpu_node) specifys the number of threads actually started.
   The number of threas started is returned.
*/
int cbe_init(int const nspu_req, int const cpu_node) {
    /* Prevent init routine from being called more than once
       using the following flag.
     */
    static int flag[CBE_GRANULE_NINT] = { 0 };
    if ( 0 == _atomic_inc_return(EA_CAST(flag)) ) {
        /* 1st. Init */
        return cbe_init0(nspu_req,cpu_node);
    }

    /* Alread initialized */
    _atomic_set(EA_CAST(flag), 1);
    return 1;
}


void cbe_shutdown(void) {
    unsigned i;
    spe_context_ptr_t* pctxt;

    for ( i=0, pctxt=c;    i<N_SPU_THREADS_MAX;    ++i, ++pctxt ) {
        spe_context_destroy(*pctxt);
        (*pctxt)=0;
        contrarea[i]=0;
    }
}


#endif /* PPU specific part */
