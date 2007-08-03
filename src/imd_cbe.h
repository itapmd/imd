#ifndef IMD_CBE_UTIL_H__
#define IMD_CBE_UTIL_H__


/********* Common (SPU & PPU) part +++++++++*/

#if defined (__cplusplus)
extern "C" {
#endif


/* The ALIGNED_(b,v) macro is used to align
   varibale v to a boundary specified by b.
   E.g. Use it to specify 16-byte alignment (which is needed for DMA)
   as follows:   buffer_type ALIGNED_(16, buf);
 */

/* __attribute__ is define on both PPU & SPU if gcc is used */
#if defined(__GNUC__)
#    define HAVE_ATTRIBUTE_ALIGNED
#endif

/* Some attributes used to specify alignment */
#if defined(HAVE_ATTRIBUTE_ALIGNED)
#   /* Define ALIGNED to used the gcc attribute if it is not already defined */
#   if ! defined(ALIGNED_)
#       define ALIGNED_(alg, var) var __attribute__((aligned(alg)))
#   endif
#else
#  /* Just define ALIGNED to be the variable */
#  /* a NOOP */
#  define ALIGNED_(alg, var) var
#endif




/* Some branch prediction macros */

/* spu-gcc has __builtin_expect */
#if defined(__SPU__)
#     define HAVE_BUILTIN_EXPECT
#endif

/* Define EXPECT_TRUE & EXPECT_FALSE if they have not yet been */
/* defined somewhere else */

/* Use the following macros in if-statements */
#if defined(HAVE_BUILTIN_EXPECT)
#   if ! defined(EXPECT)
#       define EXPECT(x,v)   __builtin_expect((x),(v))
#   endif
#else
#   /* The macros are just replicate the argument if __builtin_expect */
#   /*  is not available */
#   if ! defined(EXPECT)
#        define EXPECT(x,v)  (x)
#   endif
#endif

/* Define EXPECT_TRUE/FALSE in terms of EXPECT defined above */
# if ! defined(EXPECT_TRUE)
#    define EXPECT_TRUE(x)   EXPECT((x), 1)
#endif
#if ! defined(EXPECT_FALSE)
#    define EXPECT_FALSE(x)  EXPECT((x), 0)
#endif





/* Some casts (as macros).  Use with care and use sparingly!
   CAST(t,x) casts expression x to type t.
   ACAST(t,x) casts the address of x to a pointer to t (that is t*)
 */
#if defined (__cplusplus)
#    define  CAST(typ, expr)  (reinterpret_cast<typ> (  expr))
#    define ACAST(typ, expr)  (reinterpret_cast<typ*>(&(expr)))
#else
#    define  CAST(typ, expr)  ((typ )  (expr))
#    define ACAST(typ, expr)  ((typ*)(&(expr)))
#endif





/* Some utility "functions" for arrays defined as macros  */
/* Number of elements in array a (May be evaluated at compile time) */
#define ASIZE(a) ((sizeof (a))/(sizeof (a[0])))

/* Get pointer (itartor) to first element and the 
   past-the-end-iterator of the array a
 */
#define AFRST(a) (&(a[0]))
#define ALAST(a) ((a)+ASIZE(a))

/* Warning: Do not use pointers as arguments */


/* An effective address is just an array of unsigneds 
   ea_t[0] (the first component) is the high part
   ea_t[1] (the second component) is the low part
 */


/* The floating point type */
typedef float flt;

typedef unsigned ea_t[2];
typedef ea_t ea32array_t;


/* The work package type */
typedef struct {
  int   k, n1, n2, len;           /* package number and actual sizes  */
  int   n1_max, n2_max, len_max;  /* allocated size (not transferred) */
  flt   totpot, virial;
  flt   *pos, *force;             /* length: 4*n2, 4*n2 */
  int   *typ, *ti;                /* length: n2, 2*n2   */
  short *tb;                      /* length: len        */
} wp_t;

/* Basically the same structure as the wp_t but with ptrs. replaced
   by effective addresses and aligned to 16 or even 128-bytes boundary */
typedef struct exch {
    /* Non-ptr. variables which may be passed "as is" */
    int   k, n1, n2, len;           /* package number and actual sizes  */
    int   n1_max, n2_max, len_max;  /* allocated size (not transferred) */
    flt   totpot, virial;

    ea32array_t   pos, force;             /* length: 4*n2, 4*n2 */

    /* "Pointers" to integer arrays
        typ "points" to an array of n2 items of type int
        tb "points" to len items of type short
     */
    ea32array_t  typ, ti;                /* length: n2, 2*n2   */
    ea32array_t  tb;                     /* length: len        */


    /* Padding */
    unsigned char pad[4];
} exch_t;



/* The potential type */
typedef struct {
  int ntypes;
  flt *r2cut, *lj_sig, *lj_eps, *lj_shift;
} pt_t;

/* This should be defined in spu.c or imd_force_cbe */
extern pt_t pt;

/* ptr<->ea_t mapping is similar to wp_t */
/* This type is supposed to be passed as envirnoment pointer */
typedef struct env {
    /* Copy of ntypes */
    int ntypes;

    /* Each of the following "points" to arrays with
       4*ntypes*ntypes items of type flt */
    ea32array_t r2cut, lj_sig, lj_eps, lj_shift; 



    /* Padding */
    unsigned char pad[12];
} env_t;


/* The main calculation routine */
void calc_wp(wp_t *wp);

/* Transform effective address pairs to DMA list elements:

   The hi part of effective addresses is written to ea_pairs[2*N] if
   all hi parts in ea_pairs are the same. In this case the hi parts will be
   overwritten with the corresponding sizes. This way, eas[0]...eas[2*N-1] may
   be used as list in list DMAs.

   NULL is returned to indicate that there are some hi parts which differ and 
   thus list DMA is not possible.

   Otherwise, a ptr. to the  beginning of the list (that is ea_pairs)
   is returned.

   eas_pairs[2*N] is written to in any case.
  

   This function may cause undefined behaviour if N<1
 */
unsigned* eas2les(unsigned ea_pairs[], unsigned const sizes[], unsigned const N);





/* Check wether all EAs in *p are aligned to boundary specified by a*/
int (exch_aligned)(exch_t const* const p, unsigned const a);

/* Check wether all EAs in *p are aligned to boundary specified by a*/
int (env_aligned)(env_t const* const p, unsigned const a);




/* Output of env_t *e using printf-like function of */
int (env_out)(int (*of)(char const[],...), env_t const* const e);

/* Output of *e unsing printf-like function of */
int (exch_out)(int (*of)(char const[],...), exch_t const* const e);


/* Print vector */
int vecout(int (*of)(char const[],...),   /* General output function*/
           int (*elemout)(int (*of2)(char const[],...), void const*), /* Perelement output function using general function */
           void const* (*nxt)(void const*),  /* Get next element */
           void const* s, unsigned n, /* n elements starting at s */
           char const sep[]
           );


/* Output of pt_t *e using  printf-like function of */
void (pt_out)(int (*of)(char const[],...), pt_t const* const e);

/* Output of work package wp using printf-like functions of */
void (wp_out)(int (*of)(char const[],...), wp_t const* const e);


#if defined (__cplusplus)
}
#endif
















/************** SPU part ************/
#if defined(__SPU__)
#endif /* SPU part */

#if defined (__cplusplus)
extern "C" {
#endif

#if defined (__cplusplus)
}
#endif

#endif /* SPU */



















/************ PPU part **************/
#if ! defined(__SPU__)

/* SPU runtime management */
#include <libspe2.h>

#if defined (__cplusplus)
extern "C" {
#endif

/* Get a buffer for the work package */
wp_t* cbe_get_wp(int const k);

/* Get exch_t for SPU k */
exch_t* cbe_get_argp(int const k, wp_t const* const pwp);

/* Get env_t for SPU k */
env_t* cbe_get_envp(int const k, pt_t const* const ppt);

/* Get instruction pointer for SPU k */
unsigned* cbe_get_spuentry(int const k, unsigned const entry0);

/* Same as spuentry but with instruction counter set to default start value */
unsigned* cbe_get_spustart(int const k);

/* Get stop info buffer for SPU k */
spe_stop_info_t* cbe_get_stopinfop(int const k);

/* Get context k */
spe_context_ptr_t cbe_get_context(int const k);


/* CBE specfic initialization.
   This function must be called prior to using any other function or
   object from imd_cbe_util.h if not stated other wise
  */
void cbe_init(void);

/* CBE specific cleanup */
void cbe_shutdown(void);






/* CBE time base utilities */

/* The used to represent the time base */
typedef unsigned long timebase_t;

/* Read timebase from open file */
timebase_t timebase_file(FILE* const, char const line[]);

/* Read timebase from open file */
timebase_t timebase_path(char const path[], char const line[]);

/* Just read the usual timebase entry from the usual file */
timebase_t timebase_proc(void);




/* Memory allocation/alignment */

/* Return the page size */
size_t pagesize(void);

/* Allocate at least s bytes of memory which are aligned ton an alig0
   boundary. The requested size is rounded up towards the next higher 
   multiple of mult0 */
void* (malloc_aligned)(size_t const s0, size_t const alig0, size_t const mult0);

/* The same as above allocation nelem elements of size elsze each */
void* (calloc_aligned)(size_t const nelem, size_t const elsze, 
                       size_t const alig0, size_t const mult0);




/* Schedule work (package) to SPU  */
void schedtospu(wp_t* const pwp);




/* Create env_t from pt_t */
env_t* (create_env)(pt_t const* const p, env_t* const e);

/* Create exch_t *e from wp_t *wp  */
exch_t* (create_exch)(wp_t const* const wp, exch_t* const e);



#if defined (__cplusplus)
}
#endif


#endif /* PPU part */
