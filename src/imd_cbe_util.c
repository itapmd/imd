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




















































/* Helper function for conversion functors */
static unsigned int add_digit(unsigned int const n,
                              unsigned int B, char const d)
{
    /* n shifted to the left by one position (in its base B representation) */
    unsigned int const nB = B*n;

    /* Append digit (if it is a valid one) */
    if ( '0' == d ) { return nB;    }
    if ( '1' == d ) { return nB+ 1; }
    if ( '2' == d ) { return nB+ 2; }
    if ( '3' == d ) { return nB+ 3; }
    if ( '4' == d ) { return nB+ 4; }
    if ( '5' == d ) { return nB+ 5; }
    if ( '6' == d ) { return nB+ 6; }
    if ( '7' == d ) { return nB+ 7; }
    if ( '8' == d ) { return nB+ 8; }
    if ( '9' == d ) { return nB+ 9; }
    if ( 'a' == d ) { return nB+10; }
    if ( 'b' == d ) { return nB+11; }
    if ( 'c' == d ) { return nB+12; }
    if ( 'd' == d ) { return nB+13; }
    if ( 'e' == d ) { return nB+14; }
    if ( 'f' == d ) { return nB+15; }

    /* Invalid digit, so just return the original number */
    return n;
}

/* dec,hex,bin,oct may be used as conversion functor U in str2ui */
unsigned int dec(unsigned int n, char d) {
    return add_digit(n, 10u, d);
}

unsigned int hex(unsigned int n, char d) {
    return add_digit(n, 0x10u, d);
}

unsigned int bin(unsigned int n, char d) {
    return add_digit(n,  2u, d);
}

unsigned int oct(unsigned int n, char d) {
    return add_digit(n, 010u, d);
}


/* Convert digit string str to unsigned using conversion functor U */
unsigned int addstr2ui(unsigned int n,
                       char const* const sbeg, char const* const send,
                       unsigned (* const U)(unsigned cur, char digit)
                      )
{
    if ( (0!=sbeg) && (0!=send) ) {
        register char const* p=sbeg;
   
        for(   ;   (send!=p);    ++p ) {
            n = U(n,*p);
        }
    }

    /* Return updated (increased) n  */
    return n;
}


/* Find first subsequence [subfrst,sublast) in [tfrst,tlast)
   returning a pointer to the start of that subsequence
 */
char const* strmatch(register char const* const tfrst, register char const* const tlast, 
                     register char const* const subfrst, register char const* const sublast)
{
    /* Iterator over t range */
    register char const* it=tfrst;
    for (  ;  (tlast!=it);  ) {
        if ( (*it) == (*subfrst)  ) {
   	    /* Found 1st matching char. */
  	    register char const* it2 = it;
  	    register char const* isub2 = subfrst;
            /* 1st character matched, so move on to next one */
            ++it2;
            ++isub2;

            for ( ;   ((sublast!=isub2) && (tlast!=it2));   ++it2, ++isub2 ) {
	        if ( ! ((*it2) == (*isub2)) ) {
		   /* At least one character didn't match */
		   goto next; /* We can't use break here. */
                }
            }
            /* Found a matching sequence */
            return it;
        }
        /* Just move on to next item in sequence t. */
   	next:
        ++it;
    }

    /* No match found */
    return tlast;
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



#endif /* SPU part */
























#if defined(__PPU__)   /************ PPU part **************/

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












/* Get timebase frequency from contents of cpuinfo file */
unsigned long tbfreq_string(char const* const sbeg, char const* const send) {
    if ( (0!=sbeg) && (0!=send) ) {
        /* (Pointers to) constant characters are often needed in this block,
           so we introduce the following typedefs */
        typedef char const Tuc, *Puc;


        /* Patterns/strings to be searched for */
        /*                  01234567  */
        static Tuc key[] = "timebase";
        Puc const keypos = strmatch(sbeg, send,  key, key+(sizeof(key)-1));

        if ( keypos != send ) {
            static Tuc col[] = ":";
            Puc const colpos = strmatch(keypos+1, send,  col, col+(sizeof(col)-1));

            if ( colpos != send ) {
	        static Tuc eol[] = "\n";
  	        Puc const eolpos = strmatch(colpos+1, send,  eol, eol+((sizeof(eol)-1)));
                if ( send != eolpos ) {
		    return addstr2ui(0u,  colpos, eolpos,  dec);
	        }
            }
        }   
    }

    /* Could not convert, return (non sense) frequency 0 */
    return 0u;
}


/* Read cpuinfo from fd into buffer buf of length bufsze parse it and return
   the timebase freqeuncy. */
unsigned long tbfreq_fd(int const fd,  char* const buf, unsigned const bufsze) {    /* Got valid handle? */
    if ( -1 != fd ) {
       /* Read data from fd into the following buffer, storing the
          number of characters read */
       int const nread = read(fd, buf,bufsze);

       /* Successfully read the entire file? */
       if ( (-1!=nread) && (nread<=bufsze-1) ) {
  	    /* Terminate string just read and pass it to the parsing function
               to get the final result */
     	    return tbfreq_string(buf, buf+nread);
       }
    }

    /* Return 0 to indicate an error */
    return 0u;
}


/* Get time base frequency from OS (as suggested by the CBE programming handbook) */
unsigned long tbfreq(void) {
    /* Open cpuinfo file */
    static char const cipath[] = "/proc/cpuinfo";
    int const fd = open(cipath, O_RDONLY);

    /* Input buffer (hopefully large enough ) */
    char buf[4*1024];

    /* Get timebase frequency reading fd. */
    unsigned int const res = tbfreq_fd(fd, buf, (sizeof buf));

    /* Close file before returning */
    (void)(close(fd));

    /* Return the result */
    return res;
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
/* static envbuf_t ALIGNED_(128, ebuf[N_ENVBUF]); */
/* cbe_env_begin is passed as environment address to every SPU thread */
/* envbuf_t* const cbe_env_begin = ebuf; */

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
static unsigned int nspus=0;

/* Get number of threads initialized */
unsigned cbe_get_nspus() {
    return nspus;
}


/* SPU multithreading:  IDs of the POSIX threads */
static Tthr  t[N_SPU_THREADS_MAX];




/* The "real" initialization routine perfomring the actual work */
int cbe_init0(int const nspu_req, int const cpu_node)
{
    /* The env */
    static env_t ALIGNED_(16,env);

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
    cbe_envea_begin[0] = EA_CAST(&env);
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



    /* Copy pt to env_t
       It is assumend that pt has been setup before we arrive here */
    create_env(&pt, &env);




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





#endif /* PPU part */
