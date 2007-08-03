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


























/************** SPU part ************/
#if defined(__SPU__)
  /* Empty */
#endif /* SPU part */























/************ PPU part **************/
#if ! defined(__SPU__)

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




/* Some local abreviations for two libspe2 data types
   which are often used here */
typedef spe_program_handle_t Thndl;
typedef spe_context_ptr_t    Tctxt;
/* The same for the pthread type */
typedef pthread_t            Tthr;





/* Number of SPU threads to be used  */
enum {
   N_SPU_THREADS_MAX = 
      (unsigned)
#if defined (SPU_THREADS_MAX)
      (SPU_THREADS_MAX) /* Use the #defined value */
#else
      (6) /* Use default which is also OK for the PS3 */
#endif
};


/* Check k */
static int (chkk)(int const k) {
    return ((k>=0) && (k<N_SPU_THREADS_MAX));
}

/* Implemented as a macro, in which k appears twice */
/* #define chkk(k) ((k>=0) && (k<N_SPU_THREADS_MAX)) */




/* Get wp for SPU k */
static wp_t wp[N_SPU_THREADS_MAX];

/* Get a buffer for the work package */
wp_t* cbe_get_wp(int const k)
{
    return ((chkk(k))  ?  (wp+k)  :  0);
}



/* Get exch_t for SPU k */
exch_t* cbe_get_argp(int const k, wp_t const* const pwp)
{
    /* Ptr to buffer for SPU k */
    exch_t* res;

    if ( chkk(k) ) {
        /* This must be aligned to 16byte boundary */
        static exch_t ALIGNED_(16, exch[N_SPU_THREADS_MAX]);
        res=exch+k;
    }
    else {
        res=0;
    }

    /* Init result with data from pwp (if supplied) */
    if ( pwp ) {
        if ( res ) {
            (void)(create_exch(pwp,res));
            /* fprintf(stdout, "exch_t created.\n"); */
        }
    }
    return res;
}


/* Get env_t for SPU k */
env_t* cbe_get_envp(int const k, pt_t const* const ppt)
{
    /* Ptr to buffer for SPU k */
    env_t* res;

    if ( chkk(k) ) {
         /* This must be aligned to 16byte boundary */
         static env_t ALIGNED_(16, env[N_SPU_THREADS_MAX]);
         res=env+k;
    }
    else {
        res=0;
    }

    /* Init result with data from ppt (if supplied) */
    if ( ppt ) {
        if ( res ) {
            (void)(create_env(ppt,res));
            /* fprintf(stdout, "env_t created.\n"); */
        }
    }
    return res;
}


/* Get instruction pointer for SPU k */

static unsigned spuentry[N_SPU_THREADS_MAX] = { SPE_DEFAULT_ENTRY };

unsigned* cbe_get_spuentry(int const k, unsigned const entry0)
{
    /* Instr. pointer for SPU k */
    unsigned* const res = ((chkk(k)) ? (spuentry+k) : 0);
    
    /* Init result with entry0 */
    if ( res ) {
       (*res) = entry0;
    }
    return res;
}

/* Same as spuentry but with instruction counter set to default start value */
unsigned* cbe_get_spustart(int const k)
{
     return cbe_get_spuentry(k, SPE_DEFAULT_ENTRY);
}



/* Get stop info buffer for SPU k */
spe_stop_info_t* cbe_get_stopinfop(int const k)
{
    /* Global stopinfo buffers */
    if ( chkk(k) ) {
        static spe_stop_info_t sinfo[N_SPU_THREADS_MAX];
        return sinfo+k;
    }

    return 0;
}




/* Program handle to be used. Returns NULL if valid handle can't be found */
static Thndl* cbe_get_handle(int const k)
{
    /* The same handle for all SPU as they are all supposed to do the same. */
    if ( chkk(k) ) {
        extern spe_program_handle_t h;
        return &h;
    }


    /* Error: Index out of range */
    return 0;
}





/* ...for the SPU pthreads */
static Tthr  t[N_SPU_THREADS_MAX];
static Tctxt c[N_SPU_THREADS_MAX];


/* My NULL-ptrs */
#define CTXT0 ((Tctxt)(0))
#define HNDL0 ((Thndl*)(0))

/* Get context k */
spe_context_ptr_t cbe_get_context(int const k)
{
    return (chkk(k)) ? (c[k]) : (CTXT0);
}




/* CBE specfic initialization.
   This function must be called prior to using any other function or
   object from imd_cbe_util.h if not stated other wise
  */
void cbe_init(void)
{
     /* Indices, iterators */
     int i;
     Tctxt* pctxt;



     /* Debugging message */
     /* fprintf(stdout, "cbe_init()\n"); */

     for ( i=0, pctxt=c;      i<N_SPU_THREADS_MAX;      ++i, ++pctxt ) {
          /* Set wp fields to zero, set instruction counters */
          wp[i].n2_max=0;
          spuentry[i]=SPE_DEFAULT_ENTRY;

          /* Create a new context */
          if (  CTXT0 != (*pctxt = spe_context_create(0,0))  ) {
  	       /* Context created! Now, get the corresponding handle */
  	       spe_program_handle_t* const hndl = cbe_get_handle(i);
               if ( HNDL0 != hndl ) {
		   /* Got handle. Now load the program into context,
                      then move on to
                      next context */
		   if ( 0 == spe_program_load(*pctxt,hndl) ) {
		       continue;
                   }
               }
               /* Could not get handle or could not load program into ctxt */
  	       spe_context_destroy(*pctxt);
               *pctxt=CTXT0;
          }
     }
}

/* CBE specific cleanup */
void cbe_shutdown(void)
{
     /* Indices, iterators */
     int i;
     spe_context_ptr_t* pctxt;

#if DBGOUT
     /* Debugging message */
     fprintf(stdout, "cbe_shutdown()\n");
#endif

     /* Destroy contexts */
     for ( i=N_SPU_THREADS_MAX, pctxt=c;     i>0;     --i, ++pctxt ) {
          spe_context_destroy(*pctxt);
          *pctxt=CTXT0;
     }
}

#undef CTXT0
#undef HNDL0









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



/* Secure version of strcpy/strcat for array A as dst. argument
   Those can't be implemented as functions. Even thoug A appears twice
   in the macro, the second occurrance is as an argument to the sizeof operator
 */
#define STRACPY(A,ct) strncpy((A), ct, sizeof(A))
#define STRACAT(A,ct) strncat((A), ct, sizeof(A))
/* Macro version of fget which deduces the size of the buffer from the argument*/
#define FGETA(A,f) fgets(A, (sizeof (A)), f)


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











/* Schedule workpackage to SPU */


/* Schedule work (package) to SPU  */
void schedtospu(wp_t* const pwp)
{
     /* The SPU to be used */
     int const spunr = 0;

     /* Location of "controlblocks" for DMA */
     exch_t* const parg = cbe_get_argp(spunr, pwp);
     env_t* const  penv = cbe_get_envp(spunr, &pt);

     if ( ! (parg && penv) ) {
         fprintf(stderr, "Could not get buffers for argument/environment!\n");
         exit(2);
     }

     /* Create the control blocks & make sure the EAs contained therein
        are aligned */
     if ( !  exch_aligned(parg, 16) ) {
        fprintf(stderr, "Effective address in wp is not aligned appropriatly.\n");
        exit(2);
     }
     /* Debugging output */
     /* exch_out(printf, parg); */

     if ( ! env_aligned(penv, 16) ) {
        fprintf(stderr, "Effective address in pt is not aligned appropriatly.\n");
        exit(2);
     }
     /* Debugging output */
     /* env_out(printf, penv); */


     /* Debugging output */
     /* fprintf(stdout, "Switching to SPU context in calc_wp_on_SPU\n"); fflush(stdout); */

     /* Change to SPU context passing the addresses of the control blocks
       as argp & envp. */
     spe_context_run(cbe_get_context(spunr),  cbe_get_spustart(spunr),
                     0, parg, penv, cbe_get_stopinfop(spunr)
                    );

     /* Debugging output */
     /* fprintf(stdout, "Back from SPU context\n");  fflush(stdout); */



     /* To update wp in main mem, only the componenet totpot & virial 
        have to be copied from the exch_t, as the force has been DMAed
        directly */
     pwp->totpot = parg->totpot;
     pwp->virial = parg->virial;
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
         *(--dstlast)=0;
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






#endif /* PPU part */
