/* The header file containing the declarations */
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





/* Check wether effective address a is aligned to boundary specified by b.
   This macro assumes that an effective address is an array of unsigneds
   with the lowest part at the end of the array. e must be a array and must
   not be a pointer.

   If needed, this macro might be moved to the .h header in order to make
   it availably to client code.
 */
#define EA_ALIGNED(e,b) (0 == e[((sizeof e)/(sizeof e[0]))-1] % (b))

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
    RET(p->pos,a);
    RET(p->force,a);
    RET(p->typ,a);
    RET(p->ti,a);
    RET(p->tb,a);

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
#define EA(ptr) (*(ptr)), (*((ptr)+1))
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

#undef EA
#undef EAFMT





/* Print vector */
int vecout(int (*of)(char const[],...),   /* General output function*/
           int (*elemout)(int (*of2)(char const[],...), void const*), /* Perelement output function using general function */
           void const* (*nxt)(void const*),  /* Get next element */
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



static void const* shtnxt(void const* p) {
     return ((short const*)p)+1;
}

static int shtout(int (*of)(char const[],...), void const* pshrt) {
     return of("%hd", *((short const*)pshrt));
}




/* Output of pt_t *e using  printf-like function of */
void (pt_out)(int (*of)(char const[],...), pt_t const* const e)
{
     /* Sepeartor */
     static char const sep[]=" ";

     /* Number of elements in each array */
     unsigned const nout = 4*(e->ntypes)*(e->ntypes);

     /* Scalar */
     of("ntypes   = %d\n", (e->ntypes));

     /* Arrays */
     of("r2cut   = ");
     vecout(of, fltout,fltnxt, e->r2cut,    nout, sep);
     of("\n");

     of("lj_sig   = ");
     vecout(of, fltout,fltnxt, e->lj_sig,   nout, sep);
     of("\n");

     of("lj_eps   = ");
     vecout(of, fltout,fltnxt, e->lj_eps,   nout, sep);
     of("\n");

     of("lj_shift = ");
     vecout(of, fltout,fltnxt, e->lj_shift, nout, sep);
     of("\n");
}

/* Output of work package wp using printf-like functions of */
void (wp_out)(int (*of)(char const[],...), wp_t const* const e)
{
   /* Seperator */
   static char const sep[] = " ";

   /* Scalars */
   of("k       = %d\n", e->k);
   of("n1      = %d\n", e->n1);
   of("n1_max  = %d\n", e->n1_max);
   of("n2      = %d\n", e->n2);
   of("n2_max  = %d\n", e->n2_max);
   of("len     = %d\n", e->len);
   of("len_max = %d\n", e->len_max);
   of("totpot  = %f\n", e->totpot);
   of("virial  = %f\n", e->virial);

   of("pos     = ");
   vecout(of, fltout, fltnxt,  e->pos,    4*(e->n2), sep);
   of("\n");

   
   of("force   = ");
   vecout(of, fltout, fltnxt,  e->force,  4*(e->n2), sep);
   of("\n");

   of("typ     = ");
   vecout(of, intout, intnxt, e->typ,    (e->n2),   sep);
   of("\n");

   of("tb      = ");
   vecout(of, shtout, shtnxt, e->tb,      (e->len),  sep);
   of("\n");

   of("ti      = ");
   vecout(of, intout, intnxt, e->ti,      2*(e->n2), sep);
   of("\n");
}



























#if defined(__SPU__)  /************** SPU part ************/
/* 
  Empty
*/
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







/* Some local abreviations for two libspe2 data types
   which are often used here */
typedef spe_program_handle_t Thndl;
typedef spe_context_ptr_t    Tctxt;
/* The same for the pthread type */
typedef pthread_t            Tthr;
















/* Work packages */
static wp_t w[N_SPU_THREADS_MAX * N_BUFLEV];
wp_t* const cbe_wp = w;


/* DMA exchange buffers, aligned 16byte boundary */
static exch_t ALIGNED_(16, exch[N_SPU_THREADS_MAX * N_BUFLEV]);
exch_t* const cbe_exch = exch;


/* DMA env. buffer, aligned to 16byte boundary */
static env_t ALIGNED_(16, env[1]);
env_t* const cbe_env = env;



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







/* Schedule work (package) to SPU  */
void schedtospu(wp_t* const pwp)
{
     calc_wp(pwp);
}











/* Helper function which converts pointers to EAs */
static unsigned* cpyea(void const* const ptr, unsigned* dstfrst, unsigned* dstlast)
{
    /* Last unsigned in ptr. representation */
    unsigned const* const srcfrst = (unsigned const*)(&ptr);
    unsigned const*       srclast = srcfrst + ((sizeof ptr)/(sizeof *dstfrst));


    /* Copy from source */
    while (  (dstfrst!=dstlast)  &&   (srcfrst!=srclast) ) {
        *(--dstlast)=*(--srclast);
    }
    /* Zero fill */
    while ( dstfrst != dstlast ) {
         *(--dstlast)=0u;
    }


    /* Return ptr. to beginning of buffer */
    return dstlast;
}




/* Copy pointer and non-pointer members as well as global variables */
/* In the following dst must be something link structure. or pointer-> 
   member must be the name of a member of the corresponding struct. */
#define CPY(src,dst,member)         ((dst).member)=((src).member)
#define CPYMEMBPTR(src,dst,member)  (cpyea(((src).member), AFRST((dst).member), ALAST((dst).member)))
#define CPYGLB(dst,global)          (((dst).global)=global)
#define CPYGLOBPTR(dst, global)     (cpyea(global, AFRST((dst).global),ALAST((dst).global)))


/* Create env_t from pt_t */
env_t* (create_env)(pt_t const* const p, env_t* const e)
{
     /* The int member may just be copied */
     CPY(*p,*e, ntypes);

     /* The ptrs have to be cast to effective addresses */
     CPYMEMBPTR(*p,*e, r2cut); 
     CPYMEMBPTR(*p,*e, lj_sig);
     CPYMEMBPTR(*p,*e, lj_eps);
     CPYMEMBPTR(*p,*e, lj_shift);

     /* Return ptr. to the env_t */
     return e;
}



/* Create exch_t *e from wp_t *wp  */
exch_t* (create_exch)(wp_t const* const wp, exch_t* const e)
{
    /* Some types may be copied 1:1 */
    CPY(*wp,*e, k);
    CPY(*wp,*e, n1);
    CPY(*wp,*e, n1_max);
    CPY(*wp,*e, n2);
    CPY(*wp,*e, n2_max);
    CPY(*wp,*e, len);
    CPY(*wp,*e, len_max);
    CPY(*wp,*e, totpot);
    CPY(*wp,*e, virial);


    /* Pointers have to be cast to integers (ea_t), asserting that they
       are aligned (at least) to a byte boundary  */
    CPYMEMBPTR(*wp,*e, pos);
    CPYMEMBPTR(*wp,*e, force);
    CPYMEMBPTR(*wp,*e, typ);
    CPYMEMBPTR(*wp,*e, ti);
    CPYMEMBPTR(*wp,*e, tb);

    /* Return ptr. to the exch_t */
    return e;
}


/* Copy pointer and non-pointer members as well as global variables */
/* In the following dst must be something link structure. or pointer-> 
   member must be the name of a member of the corresponding struct. */
#undef CPY
#undef CPYMEMBPTR
#undef CPYGLB
#undef CPYGLOBPTR







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


/* SPU multithreading: */
/* IDs of the POSIX threads */
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
    int const imax = min_usable_spes(((nspu_req < N_SPU_THREADS_MAX) ? nspu_req : N_SPU_THREADS_MAX), cpu_node);

    /* Some indices/iterators */
    int i;
    Tctxt* pctxt;
    spu_pthr_arg_t* parg;


    /* No SPUs available */
    if ( imax <1 ) {
        return nspus=0;
    }


    /* Copy pt to env_t
       It is assumend, that pt has been setup before we arrve here */
    create_env(&pt, cbe_env+0);


    /* Set "ground state" for some publically accessible values */
    for ( i=0;   i<(N_SPU_THREADS_MAX * N_BUFLEV);  ++i  ) {
        w[i].n2_max = 0;
        c[i]=0;
        contrarea[i]=0;
    }

    for ( i=0, pctxt=c, parg=pthrargs;    (i<imax);   )  {
        /* Set the arguments for the pthread */ 
        /* SPU thread arguments, also set wp fields to initial value
           such that the wp gets allocated upon first call to make_wp */
        parg->spu_arg=cbe_exch + (i*N_BUFLEV);

        /* Environment */
        parg->spu_env=cbe_env+0;

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
