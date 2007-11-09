/* IMD config (e.g. CBE_DIRECT macro) */
#include "config.h"
/* CBE specific stuff */
#include "imd_cbe.h"


/********* Common (SPU && PPU) part +++++++++*/


unsigned* eas2les(unsigned ea_pairs[], unsigned const sizes[], unsigned const N){
    /* Pointer/index */
    register unsigned const* psrc;
    register unsigned* pdst;
    register unsigned k;

    /* Pointer to hi part: Store the common high part of all EAs
       right after the list elements */
    unsigned* const phi = ea_pairs+(N<<1u);

    /* Check wether all hi parts are the same and return -1 if
       they are not
    */
    for ( *phi=*(psrc=ea_pairs), psrc+=2, k=N;    (k>1);     psrc+=2, --k ) {
        if ( (*phi) != (*psrc) ) {
	  /* There are EAs in the list which are not equal */
  	   /* DBGPRINTF(stdout, "%u %u\n", (*phi), (*psrc)); */
  	   return 0;
        }
    }


    /* DBGPRINTF(stdout, "Common high part of all EAs: %u\n", (*phi)); */

    /* All hi parts are the same (their value has been save in hi),
       such that they can now be overwritten by the size arguments */
    for ( psrc=sizes, pdst=ea_pairs,  k=N;         k>0;    pdst+=2, ++psrc, --k ) {
        /* Copy 16 bits of the size argument, leaving the "Reserved" part
           untouched. Also, the stall-notify bit will not be set  */
        enum { msk=(1u<<17u)-1u };
        (*pdst) = (*psrc);
        (*pdst) &= msk;

        /* DBGPRINTF(stdout, "DMA-List element: %u bytes @ (0x%x,0x%x)\n",  (*pdst), (*phi), *(pdst+1)); */
    }

    /* Everything went OK */
    return ea_pairs;
}



/* Print sizes of some types */
int sizeinfo(int (* const of)(char const[],...)) {
    return
      of("sizeof(wp_t)       = %u\n"
         "sizeof(exch_t)     = %u\n"
#if defined(CBE_DIRECT)
         "sizeof(cell_dta_t) = %u\n"
         "sizeof(cell_ea_t)  = %u\n"
#endif
         ,

         (unsigned)(sizeof(wp_t)),
         (unsigned)(sizeof(exch_t))
#if defined(CBE_DIRECT)
         ,
         (unsigned)(sizeof(cell_dta_t)),
         (unsigned)(sizeof(cell_ea_t))
#endif
	);
}




/* Print vector */
int vecout(int (* const of)(char const[],...),   /* General output function*/
           int (* const elemout)(int (* const of2)(char const[],...), void const*), /* Perelement output function using general function */
           void const* (* const nxt)(void const*),  /* Get next element */
           void const* s, unsigned n, /* n elements starting at s */
           char const sep[]
           )
{
    /* The result: number of characters written */
    int res=0;

    if ( n ) {
       /* 0=!n that is: There are some elements to be printed */
       /* Iterate over elements */
       for( ;  s;  s=nxt(s) ) {
           if ( --n ) {
   	       /* Print element followed be a separator */
 	       int rc;
   	       if ( (rc=elemout(of,s)) < 0 ) {
		   return rc;
               }
               if ( (rc=of("%s", sep)) < 0 ) {
		   return rc;
               }
           }
           else {
	      /* Last element */
   	      int rc;
  	      if ( (rc=elemout(of,s)) < 0 ) {
	  	  return rc;
              }
              break;
           }
       }
    }

    /* Return total number of characters written */
    return res;
}



/* Helper output functors for various types */
static void const* fltnxt(void const* p) {
     return ((flt const*)p)+1;
}

static int fltout(int (*of)(char const[],...), void const* pflt) {
     return of("%f", *((flt const*)pflt));
}


static void const* intnxt(void const* p) {
     return ((int const*)p)+1;
}

static int intout(int (*of)(char const[],...), void const* pint) {
     return of("%d", *((int const*)pint));
}


/* Advance iterators */
static void const* shtnxt(void const* p) {
     return ((short const*)p)+1;
}

static int shtout(int (*of)(char const[],...), void const* pshrt) {
     return of("%hd", *((short const*)pshrt));
}





/* Check wether effective address a is aligned to boundary specified by b.
   This macro assumes that an effective address is an array of unsigneds
   with the lowest part at the end of the array. e must be a array and must
   not be a pointer.

   If needed, this macro might be moved to the .h header in order to make
   it availably to client code.
 */
/* #define EA_ALIGNED(e,b) (0 == e[((sizeof (e))/(sizeof (e)[0]))-1] % (b)) */
#define EA_ALIGNED(e,b) (0 == ((e)%(b)) )


/* Return if p is not aligned to 16 byte boundary which is the minimum
   requirement for DMAs  */
#define RET(ea,alig) { if ( ! EA_ALIGNED((ea),(alig)) ) { return 0; } }


/* Check wether all EAs in *p are aligned to boundary specified by a*/
int (exch_aligned)(exch_t const* const p, unsigned const a) {
    /* Return "error" if parameters are not OK */
    if ( 0==a ) {
        return 0;
    }

    /* Check all EA members */
#if defined(CBE_DIRECT)
    {
    int i=0;
    cell_ea_t const* pdta = p->cell_dta;
    for (  ;   (i<NNBCELL);    ++i, ++pdta ) {
        RET(pdta->pos_ea,   a);
        RET(pdta->force_ea, a);
        RET(pdta->typ_ea,   a);
        RET(pdta->ti_ea,    a);
        RET(pdta->tb_ea,    a);
    }
    }
#else
    RET(p->pos,   a);
    RET(p->force, a);
    RET(p->typ,   a);
    RET(p->ti,    a);
    RET(p->tb,    a);
#endif

    /* All EAs we checked are valid if we arrive here. */
    return 1;
}


/* Check wether all EAs in *p are aligned to boundary specified by a*/
int (env_aligned)(env_t const* const p, unsigned const a) {
    /* Return "error" if parameters are not OK */
    if ( 0==a ) {
        return 0;
    }


    /* Check all EA members */
    RET(p->r2cut,a);
    RET(p->lj_sig,a);
    RET(p->lj_eps,a);
    RET(p->lj_shift,a);

    /* All EAs we checked are valid if we arrive here. */
    return 1;
}

#undef RET
#undef EA_ALIGNED














/* Needed in the output of effective addresses  */
/* #define EA(ptr) (*(ptr)), (*((ptr)+1)) */
#define EA(ea)   ((unsigned)(((ea_t)(ea))>>32u)),  ((unsigned)((ea_t)(ea)))
#define EAFMT "(0x%x,0x%x)"

/* Output of env_t *e using printf-like function of */
int (env_out)(int (*of)(char const[],...), env_t const* const e)
{
    return
    of("ntypes   = %d\n"
       "r2cut    = " EAFMT "\n"
       "lj_sig   = " EAFMT "\n"
       "lj_eps   = " EAFMT "\n"
       "lj_shift = " EAFMT "\n",


        (e->ntypes),
        EA(e->r2cut), EA(e->lj_sig), EA(e->lj_eps), EA(e->lj_shift)
      );
}




#if defined(CBE_DIRECT)   /* "Direct " */

/* exch_t & wp_t are the same */
int (wp_out)(int (*of)(char const[],...), exch_t const* const w) {
    /* To be done */
    return -1;
}

int (exch_out)(int (*of)(char const[],...), exch_t const* const e) {
    return wp_out(of,e);
}

#else  /* "Indirect" */

/* Output of work package wp using printf-like functions of */
int (wp_out)(int (*of)(char const[],...), wp_t const* const e)
{
   /* Seperator */
   static char const sep[] = " ";

   /* The result */
   int res=0;

   /* Scalars */
   res += of("k       = %d\n", e->k);
   res += of("n1      = %d\n", e->n1);
   res += of("n1_max  = %d\n", e->n1_max);
   res += of("n2      = %d\n", e->n2);
   res += of("n2_max  = %d\n", e->n2_max);
   res += of("len     = %d\n", e->len);
   res += of("len_max = %d\n", e->len_max);
   res += of("totpot  = %f\n", e->totpot);
   res += of("virial  = %f\n", e->virial);

   res += of("pos     = ");
   res += vecout(of, fltout, fltnxt,  e->pos,    4*(e->n2), sep);
   res += of("\n");

   
   res += of("force   = ");
   res += vecout(of, fltout, fltnxt,  e->force,  4*(e->n2), sep);
   res += of("\n");

   res += of("typ     = ");
   res += vecout(of, intout, intnxt, e->typ,    (e->n2),   sep);
   res += of("\n");

   res += of("tb      = ");
   res += vecout(of, shtout, shtnxt, e->tb,      (e->len),  sep);
   res += of("\n");

   res += of("ti      = ");
   res += vecout(of, intout, intnxt, e->ti,      2*(e->n2), sep);
   res += of("\n");

   return res;
}

/* Output of *e unsing printf-like function of */
int (exch_out)(int (*of)(char const[],...), exch_t const* const e)
{
   return
   of("k       = %d\n"
      "n1      = %d\n"
      "n1_max  = %d\n"
      "n2      = %d\n"
      "n2_max  = %d\n"
      "len     = %d\n"
      "len_max = %d\n"
      "totpot  = %f\n"
      "virial  = %f\n"
      "pos     = " EAFMT "\n"
      "force   = " EAFMT "\n"
      "typ     = " EAFMT "\n"
      "tb      = " EAFMT "\n"
      "ti      = " EAFMT "\n",

      e->k, e->n1, e->n1_max, e->n2, e->n2_max, e->len, e->len_max,
      e->totpot, e->virial,
      EA(e->pos), EA(e->force), EA(e->typ),
      EA(e->tb), EA(e->ti)
     );
}

#endif  /* CBE_DIRECT */



#undef EA
#undef EAFMT












/* Output of pt_t *e using  printf-like function of */
int (pt_out)(int (*of)(char const[],...), pt_t const* const e)
{
     /* Sepeartor */
     static char const sep[]=" ";

     /* Number of elements in each array */
     unsigned const nout = 4*(e->ntypes)*(e->ntypes);

     /* The result */
     int res=0;

     /* Scalar */
     res += of("ntypes   = %d\n", (e->ntypes));

     /* Arrays */
     res += of("r2cut   = ");
     res += vecout(of, fltout,fltnxt, e->r2cut,    nout, sep);
     res += of("\n");

     res += of("lj_sig   = ");
     res += vecout(of, fltout,fltnxt, e->lj_sig,   nout, sep);
     res += of("\n");

     res += of("lj_eps   = ");
     res += vecout(of, fltout,fltnxt, e->lj_eps,   nout, sep);
     res += of("\n");

     res += of("lj_shift = ");
     res += vecout(of, fltout,fltnxt, e->lj_shift, nout, sep);
     res += of("\n");

     return res;
}





























#if defined(__SPU__)  /************** SPU part ************/

/* SPU headers */
#include <spu_mfcio.h>




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

        /* Update addresses and number of bytes */
        remsze -= DMAMAX;
        p      += DMAMAX;
        lea    += DMAMAX;
    }
}


#endif /* SPU part */
























#if defined(__PPU__)   /************ PPU part **************/

/* ISO C std. headers */
#include <limits.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/* POSIX std. headers */
#include <unistd.h>
#include <pthread.h>
/* CBE Headers */
#include <libspe2.h>




/* Secure version of strcpy/strcat for array A as dst. argument
   Those can't be implemented as functions. Even thoug A appears twice
   in the macro, the second occurrance is as an argument to the sizeof operator
 */
#define STRACPY(A,ct) strncpy((A), ct, sizeof(A))
#define STRACAT(A,ct) strncat((A), ct, sizeof(A))
/* Macro version of fget which deduces the size of the buffer from the argument*/
#define FGETA(A,f) fgets(A, (sizeof (A)), f)



/* Function to check a string:
   return s if it points to an non-empty string, dflt is returned otherwise
 */

/* Macro */
#define strchk(s,dflt) (s ? (('\0' != (*s)) ? s : dflt) :  dflt)

/* Real function using that macro */
static char const* (strchk)(char const s[], char const dflt[])
{
    return strchk(s,dflt);
}



/* The default timebase line in cpuinfo */
static char const tbase[] = "timebase";

/* Read timebase from open file */
timebase_t timebase_file(FILE* const f, char const grepval[])
{
    /* Only read from valid file */
    if ( f ) {
        char const* const tbkey = strchk(grepval, tbase);
        /* Input line by line */
        while ( ! feof(f) ) {
   	    /* Key an separator are read first with the following formats */
  	    char key[200];
            static char const fmt1[] = "%200s%*s";

    	    /* Read 1st part of line, which is key followed by a separator */
            if (  EOF !=  fscanf(f, fmt1,  key) ) {
  	          /* Get the rest of the line which is the value */
	          char val[0x100];
                  fgets(val, (sizeof val),  f);

                  /* Found requested key? */
                  if ( 0 == strcmp(key,tbkey) ) {
		       /* The result */
                       unsigned long const res=strtoul(val, NULL, 10);
                       return ((ULONG_MAX!=res) ? res : 0);
                  }
	    }
        }
    }

    /* Error */
    return 0;
}


#undef STRACAT
#undef STRACPY
#undef FGETA





/* CBE Time base functions  */

/* Helper function which also closes the file f  */
static timebase_t timebase_file_close(FILE* const f, char const line[])
{
  unsigned const res = timebase_file(f,line);
  if ( f ) {
      (void)(fclose(f));
  }
  return res;
}




/* Path to cpuinfo */
static char const cipath[] = "/proc/cpuinfo";


/* Read time base from path */
timebase_t timebase_path(char const path[], char const line[])
{
     /* File to be opened for reading */
     static char const mode[] = "r";
     return timebase_file_close(fopen(strchk(path,cipath), mode), line);
}



#undef strchk


/* Just read the usual timebase entry from the usual file */
unsigned long timebase_proc(void)
{
    return timebase_path(cipath,tbase);
}










/* Allocation of aligned memory */

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
static envbuf_t ALIGNED_(128, ebuf[N_ENVBUF]);
/* cbe_env_begin is passed as environment address to every SPU thread */
envbuf_t* const cbe_env_begin = ebuf;


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

























/* Create env_t from pt_t */
env_t* (create_env)(pt_t const* const p, env_t* const e)
{
     /* The int member may just be copied */
     e->ntypes = p->ntypes;


     /* The ptrs have to be cast to effective addresses */
     PTR2EA(p->r2cut,    e->r2cut);
     PTR2EA(p->lj_sig,   e->lj_sig);
     PTR2EA(p->lj_eps,   e->lj_eps);
     PTR2EA(p->lj_shift, e->lj_shift);


     /* Return ptr. to the env_t */
     return e;
}



/* Create exch_t *e from wp_t *wp  */
exch_t* (create_exch)(wp_t const* const wp, exch_t* const e)
{
#if defined(CBE_DIRECT)  /* "Direct" */
    /* "Direct case", that is wp_t==exch_t, such that we can just
        copy/assign the members */
   (*e) = (*wp);
#else /* "Indirect" */
    /* Some types may be copied 1:1 */
    e->k       = wp->k;

    e->n1      = wp->n1;
    e->n1_max  = wp->n1_max;

    e->n2      = wp->n2;
    e->n2_max  = wp->n2_max;

    e->len     = wp->len;
    e->len_max = wp->len_max;

    /* Do not copy totpot & virial */

    /* Pointers have to be cast to integers (ea_t), asserting that they
       are aligned (at least) to a byte boundary  */
    PTR2EA(wp->pos,    e->pos);
    PTR2EA(wp->force,  e->force);
    PTR2EA(wp->typ,    e->typ);
    PTR2EA(wp->ti,     e->ti);
    PTR2EA(wp->tb,     e->tb);

#endif
    /* Return ptr. to the exch_t */
    return e;
}










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
   return a pointer to the SPU programs stop info. */
static void* spu_pthr(void* p0)
{
    /* Check the argument ptr. and use it in case it's valid */
    if ( p0 ) {
        /* SPU instruction counter, initialized to start of SPU program */
        unsigned spu_entry = SPE_DEFAULT_ENTRY;

        /* Get the (constant) argument pointed to be ptr */
        spu_pthr_arg_t const* const p = p0;


        /* The last message from PPU */
        /* fprintf(stdout, "About to switch to SPU context\n");  fflush(stdout); */

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
static unsigned nspus=0;

/* Get number of threads initialized */
unsigned cbe_get_nspus() {
    return nspus;
}


/* SPU multithreading:  IDs of the POSIX threads */
static Tthr  t[N_SPU_THREADS_MAX];


/* Init stuff for calc_threads. Try to start SPU threads
   This function may only be called once. It's not thread safe
   The minimum of nspu_req and the number of SPUs available (according
   cpu_node) specifys the number of threads actually started.
   The number of threas started is returned.
*/
int cbe_init(int const nspu_req, int const cpu_node)
{
    /* Arguments for the SPU threads */
    static spu_pthr_arg_t pthrargs[N_SPU_THREADS_MAX];

    /* Number of threads started */
    unsigned nstrt=0;

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
        return nspus=0;
    }


    /* Copy pt to env_t
       It is assumend, that pt has been setup before we arrve here */
    create_env(&pt, ((env_t*)cbe_env_begin));


#if ! defined(CBE_DIRECT)
    for ( i=0;  i<(N_SPU_THREADS_MAX * N_ARGBUF);  ++i ) {
       wbuf[i].n2_max = 0;
    }
#endif

    /* Set "ground state" for some publically accessible values */
    for ( i=0;   i<N_SPU_THREADS_MAX ;  ++i  ) {
        /* context & control area */
        c[i]=0;
        contrarea[i]=0;
    }

    for ( i=0, pctxt=c, parg=pthrargs;    (i<imax);   )  {
        /* Set the arguments for the pthread */ 
        /* SPU thread arguments, also set wp fields to initial value
           such that the wp gets allocated upon first call to make_wp */
        parg->spu_arg = cbe_arg_begin + i*N_ARGBUF;

        /* Environment */
        parg->spu_env = cbe_env_begin;

        /* Stop info of SPU program is written here */
        parg->spu_stopinfo_loc=sinfo+i;

        /* Create context & load program into it */
        if ( 0 != ((*pctxt) = (parg->spu_ctxt) = spe_context_create(SPE_MAP_PS,0)) ) {
   	    /* This handle must be defined somewhere else and the program
               must contain the calc_wp code. */
            extern spe_program_handle_t hndle_cbe_calc;

            /* Get control area for the context just created */
            contrarea[i] = spe_ps_area_get((*pctxt), SPE_CONTROL_AREA);
            /*  fprintf(stdout, "Control area %i at %p\n", i, contrarea[i]); fflush(stdout);  */

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
                 exit(1);
            }
        }
        else {
      	   fprintf(stderr, "Could not create context\n");
        }
    }


    /* Return & store the number of threads actually started */
    return nspus=nstrt;
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





#endif /* PPU part */
