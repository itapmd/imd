#ifndef IMD_CBE_UTIL_H__
#define IMD_CBE_UTIL_H__


/********* Common (SPU & PPU) part +++++++++*/

#if defined (__cplusplus)
extern "C" {
#endif







/* Vector type in general */
#define VEC(ty) __vector ty

/* The 4 component floating point vector
   This vector type must match the flt type declared in imd_cbe.h
   The FLTVEC macro is needed for syntactical reasons when
   vector literals are used.  Whenever possible, the real typedef fltvec
   should be preferrd over the macro FLTVEC
   */
#define VECFLT VEC(float)
typedef VECFLT vecflt;


/* Mask vector type */
#define VECMSK  VEC(unsigned)
typedef VECMSK  vecmsk;



/* Initializer list for some vector literals */
#if defined(__GNUC__)
#   /* Unfortantly, gcc does not support the AltiVec list format */
#   /* but used curly braces instead */
#   define VECINIT2(a,b)     {a,b}
#   define VECINIT4(a,b,c,d) {a,b,c,d}
#else
#   define VECINIT2(a,b)     (a,b)
#   define VECINIT4(a,b,c,d) (a,b,c,d)
#endif










/* Apparently, vector types can't be initialized using the AltiVec syntax
   as the initializer expression in fact do generate instructions which
   just set the vector components to some value. So we initialize vector
   using the following union which contains on array member of the same size
   as the vector (in terms of elements). Note that the array member is
   written to at initialization time, but subsequently the vector member
   is used. Strictly speaking, this is undefined according to the ISO C
   standard, but is one way to initialize a vector according to the
   AltiVec API.
   Please also note that we can use the usual array syntax for initialization.
*/
#define VECUNION(ty) union {                         \
    ty      arr[(sizeof(VEC(ty))) / (sizeof(ty))];   \
    VEC(ty) vec;                                     \
}





#if ! defined(INLINE_)
#    if defined(__cplusplus)
#         define INLINE_ inline
#    else
#        if defined (__GNUC__)
#             define INLINE_ __inline
#        else
#             define INLINE_
#        endif /* GNUC */
#    endif  /* C++  */
#endif /* INLINE_ */

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
#ifdef __SPU__
  vector float *pos, *force;      /* length: n2, n2 */
#else
  flt   *pos, *force;             /* length: 4*n2, 4*n2 */
#endif
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
#ifdef __SPU__
  vector float *r2cut, *lj_sig, *lj_eps, *lj_shift;
#else
  flt *r2cut, *lj_sig, *lj_eps, *lj_shift;
#endif
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


/* Tokens used to synchronize between PPU / SPU */
typedef enum sync_token { MBXNONE=0u,
                          WPCREA1=1u, WPDONE1=2u, WPSTRT1=4u,
                          SPUEXIT } sync_token_t;




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
#if defined(__PPU__)

/* SPU runtime management */
#include <libspe2.h>

#if defined (__cplusplus)
extern "C" {
#endif


/* Max. number of SPU threads which may be managed by the
   following utility functions   */
enum {
   N_SPU_THREADS_MAX = 
      (unsigned)
#if defined (SPU_THREADS_MAX)
      (SPU_THREADS_MAX) /* Use the #defined value */
#else
      (32) /* Use a default which should be large enough */
#endif
};


/* Number of buffer levels, that is number of exch_t/wp_t buffers
   per SPU thread
 */
enum { N_BUFLEV=2u };


/* Will be passed as arguments to SPU programs there are at least 
   N_BUFLEV * N_SPU_THREADS_MAX elements in the array.
   cbe_exch+0*N_BUFLEV is passed to SPU 0, cbe_exch+1*N_BUFLEV is passed to
   SPU 1,   cbe_exch+k*N_BUFLEV is passed to SPU k...
*/
extern exch_t* const cbe_exch;

/* Work packages */
extern wp_t* const cbe_wp;

/* To be passed as environment to SPU programs there is at least
   one element in the array */
extern env_t* const cbe_env;


/* Stop info buffer for SPU k (readonly) */
extern spe_stop_info_t const* const cbe_stopinfo;

/* Context for SPU k (readonly) */
extern spe_context_ptr_t const* const cbe_spucontext;


/* Pointers to control areas the SPUs */
typedef spe_spu_control_area_t volatile*  spe_spu_control_area_p;
extern spe_spu_control_area_p* const cbe_spucontrolarea;



/* Get number of threads initialized by cbe_init */
unsigned cbe_get_nspus();


/* Init stuff for calc_threads. Try to start SPU threads
   This function may only be called once. It's not thread safe
   The minimum of nspu_req and the number of SPUs available (according
   cpu_node) specifys the number of threads actually started.
   The number of threas started is returned.
*/
unsigned cbe_init(unsigned const nspu_req, int const cpu_node);



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









/* Create env_t from pt_t */
env_t* (create_env)(pt_t const* const p, env_t* const e);

/* Create exch_t *e from wp_t *wp  */
exch_t* (create_exch)(wp_t const* const wp, exch_t* const e);



#if defined (__cplusplus)
}
#endif


#endif /* PPU part */
