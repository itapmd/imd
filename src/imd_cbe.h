#ifndef IMD_CBE_H__
#define IMD_CBE_H__

/* We need config here, as we have to test the CBE_DIRECT in the 
   DMA typedefs */
#include "config.h"


/********* Common (SPU & PPU) part +++++++++*/

#if defined (__cplusplus)
extern "C" {
#endif







/* Specify wether a (member) variable is to be const on PPU or SPU side
   (but not neccessarily on the other side)  */
#if defined(__SPU__)
#define ppu_const
#define spu_const const
#endif

#if defined(__PPU__)
#define ppu_const const
#define spu_const
#endif

#if ! defined(INLINE_)
#    if defined(__cplusplus)
#         define INLINE_ inline
#    else
#        if defined (__GNUC__)
#             define INLINE_ __inline
#        else
#             define INLINE_
#        endif /* __GNUC__ */
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
#define ASIZE(a) ((sizeof (a))/(sizeof ((a)[0])))

/* Get pointer (itartor) to first element and the 
   past-the-end-iterator of the array a
 */
#define AFRST(a) (&((a)[0]))
#define ALAST(a) ((a)+ASIZE(a))






/* An effective address is just an array of unsigneds 
   ea_t[0] (the first component) is the high part
   ea_t[1] (the second component) is the low part
 */
typedef unsigned  ea_t[2];



/* The floating point type */
typedef float flt;








#ifdef CBE_DIRECT  /* "Direct " */

/* Cell data within a work package;
   On SPU, allocate for n_max atoms and len_max neighbor indices,
   namely 4*n_max floats (pos, force), n_max ints (typ), 2*n_max ints (ti), 
   and len_max shorts (tb), all rounded up to 128 bytes.
   n_max and len_max are taken from wp_t; 
   n is set to the actual n from cell_ea_t. */
typedef struct {
  int   n;  /* actual number of atoms */

  /* "Output data" on SPU */
  float  *force;

  /* "Input data" on SPU */
  float  *pos;
  int    *typ;
  int    *ti;
  short  *tb;
} cell_dta_t;


/* Effective addresses of cell data within a work package;
   The cell contains n atoms and len neighbor indices;
   Transfer 4*n floats (pos, force), n ints (typ), 2*n ints (ti),
   and len shorts (tb), all rounded up to 128 bytes */
typedef struct {
  int  n, len;
  ea_t pos_ea, force_ea, typ_ea, ti_ea, tb_ea;
} cell_ea_t;


/* The work package type; essentially consists of a list of cells. 
   k       = package number, 
   n_max   = upper bound on number of atoms in a cell,
   len_max = upper bound on number of neighbor indices in a cell */
typedef struct wp {
  float totpot;
  float virial;

  int nb_max;
  int flag;

  int k;
  int n_max;
  int len_max;
  int ti_len; 

#ifdef ON_PPU
  cell_dta_t cell_dta[NNBCELL];
#else
  /* Control blocks (not to be changed on SPU after DMA) */
  cell_ea_t  cell_dta[NNBCELL];
#endif

} wp_t;


/* The wp_t may be used in DMAs directly */
typedef wp_t exch_t;


#else  /* "Indirect" */


/* The work package type */
typedef struct wp {
  int   k, n1, n2, len;           /* package number and actual sizes  */
  int   n1_max, n2_max, len_max;  /* allocated size (not transferred) */

  flt   totpot, virial;

#if defined(__SPU__)
  vector float const *pos;      /* length: n2 */
  vector float       *force;    /*         n2 */
#else
  flt       *pos;             /* length: 4*n2 */ 
  flt const *force;	      /*         4*n2 */
#endif

  int    *typ;  /* length: n2, 2*n2   */
  int    *ti;                
  short  *tb;   /* length: len        */
} wp_t;




/* Basically the same structure as the wp_t but with ptrs. replaced
   by effective addresses and aligned to 16 or even 128-bytes boundary
*/
typedef struct exch {
    /* Non-ptr. variables which may be passed "as is" */
    int k;
    int n1;
    int n2;
    int len;           /* package number and actual sizes  */
    int n1_max;
    int n2_max;
    int len_max;  /* allocated size (not transferred) */

    flt totpot;
    flt virial;

    ea_t pos;
    ea_t force;       /* length: 4*n2, 4*n2 */

    /* "Pointers" to integer arrays
        typ "points" to an array of n2 items of type int
        tb "points" to len items of type short
     */
    ea_t typ;
    ea_t ti;            /* length: n2, 2*n2   */
    ea_t tb;            /* length: len        */


    /* Padding */
    unsigned char pad[4];
} exch_t;



#endif   /* CBE_DIRECT */







/* The potential type */
typedef struct {
  int ntypes;
#if defined(__SPU__)
  vector float const* r2cut;
  vector float const* lj_sig;
  vector float const* lj_eps;
  vector float const* lj_shift;
#else
  flt* r2cut;
  flt* lj_sig;
  flt* lj_eps;
  flt* lj_shift;
#endif
} pt_t;

/* This should be defined in spu.c or imd_force_cbe */
extern  pt_t pt;






/* ptr<->ea_t mapping is similar to wp_t */
/* This type is supposed to be passed as envirnoment pointer */
typedef struct env {
    /* Copy of ntypes */
    int ntypes;

    /* Each of the following "points" to arrays with
       4*ntypes*ntypes items of type flt */
    ea_t r2cut, lj_sig, lj_eps, lj_shift; 



    /* Padding */
    unsigned char pad[12];
} env_t;


/* Tokens used to synchronize between PPU / SPU */
enum {
  MBXNONE=0u,
  WPCREA1=1u, WPDONE1=2u, WPSTRT1=4u,
  SPUEXIT
};




/* The main calculation routine(s) */
void do_work_spu(unsigned const flag);
void calc_wp(wp_t *wp);
void calc_tb(wp_t *wp);


#if defined(CBE_DIRECT)
void calc_wp_direct(wp_t*,
                    void* const, unsigned const, void* const, unsigned const,
                    unsigned const otag);

void calc_tb_direct(wp_t*,
                    void* const, unsigned const, void* const, unsigned const,
                    unsigned const otag);

#endif




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

/* Output of *e unsing printf-like function of */
int (exch_out)(int (*of)(char const[],...), exch_t const* const e);



/* Check wether all EAs in *p are aligned to boundary specified by a*/
int (env_aligned)(env_t const* const p, unsigned const a);

/* Output of env_t *e using printf-like function of */
int (env_out)(int (*of)(char const[],...), env_t const* const e);




/* Print vector */
int vecout(int (* const of)(char const[],...),   /* General output function*/
           int (* const elemout)(int (* const of2)(char const[],...), void const*), /* Per element output function using general function */
           void const* (* const nxt)(void const*),  /* Get next element */
           void const* s, unsigned n, /* n elements starting at s */
           char const sep[]
           );


/* Output of pt_t *e using  printf-like function of */
int (pt_out)(int (*of)(char const[],...), pt_t const* const e);

/* Output of work package wp using printf-like functions of */
int (wp_out)(int (*of)(char const[],...), wp_t const* const e);




/* There shall be N_ARGBUF arguement buffers for each SPU thread.
   All SPU threads shall receive a private argument address which does
   not equal an argument address passed to an other SPU.

   There shall be N_ENVBUF environment buffers which are shared by all
   SPU threads. All SPU shall receive the same environment address
 */
enum { N_ARGBUF=2, N_ENVBUF=2 };

/* Buffer padding (must be a multiple of 16, should be 128) */
enum { BUFPAD  = 128u };
enum { BUFPAD1 = BUFPAD-1u };

/* The types used as buffers for arguments and environment */
typedef unsigned char argbuf_t[(sizeof(exch_t)+ BUFPAD1) & (~BUFPAD1)];
typedef unsigned char envbuf_t[(sizeof(env_t) + BUFPAD1) & (~BUFPAD1)];


#if defined (__cplusplus)
}
#endif

















#if defined(__SPU__)          /************** SPU part ************/


/* SPU specific headers */
#include <spu_mfcio.h>
#include <spu_intrinsics.h>

#if defined (__cplusplus)
extern "C" {
#endif


/* DMA more than 16K using multiple DMAs */
static INLINE_ void mdma64(void* const p,
                           unsigned const* const ea, unsigned const size,
                           unsigned const tag, unsigned const cmd
                          )
{
    /* 64 bit effective address */
    typedef unsigned long long  T64;

    /* Defined somwhere in imd_cbe_util.c */
    extern void mdma64_rec(register unsigned char* const p, register T64 const,
                           register unsigned const,
                           register unsigned const, register unsigned const);
    extern void mdma64_iter(register unsigned char*,  register T64,
                            register unsigned,
                            register unsigned const, register unsigned const);



    mdma64_iter((unsigned char*)p, (((T64)(ea[0]))<<32u)+((T64)(ea[1])),
                size,  tag, cmd);
}



/* (Just a) wrapper for DMA with less than 16K
   This has the same signature as the function above
 */
static INLINE_ void dma64(void* const p,
                          unsigned const* const ea, unsigned const size,
                          unsigned const tag, unsigned const cmd
                  )
{
    spu_mfcdma64(p,  ea[0],ea[1],  size, tag,cmd);
}







/* Generic macro for rounding up vector components in terms of
   SPU vector commands.
   On the SPU, x may be a scalar which is converted to a vector first
 */
#define CEILV(x,r) spu_andc(spu_add((r),(x)), (r))


static INLINE_
vector unsigned uiceilv(vector unsigned const v, vector unsigned const m)
{
     return CEILV(v,m);
}





#if defined (__cplusplus)
}
#endif

#endif /* SPU */




















#if defined(__PPU__)           /************ PPU part **************/

/* WORDSIZE & size_t */
#include <limits.h>
#include <stddef.h>


/* SPU runtime management */
#include <libspe2.h>



#if defined (__cplusplus)
extern "C" {
#endif



/* Map some vector commands */
#define spu_andc  vec_andc
#define spu_add   vec_add



/* PTR2EA converts ptr to an effective address, and places it at *ea */

/* Assuming a pointer is 32 bits wide only:
   Copy ptr to lower part, setting higher part to zero
   ea must point to an array of (at least) 2 objects which are (at least)
   32 bits wide each.
   The 1st component is set to zero, the pointer is copied to the second one
*/
#if 32 == __WORDSIZE
#  /* Warning: ea is used twice inside the macro! */
#  define PTR2EA(ptr,ea) { *(ea)=0;  *((void const**)((ea)+1)) = (void const*)(ptr); }
#  define EA2PTR(ea)     (*((void**)((ea)+1)))
#endif

/* Assuming a pointer is 64 bits wide:
   ea must point to an object which is (at least) 64 bits wide
   The pointer is just copied there.
*/
#if 64 == __WORDSIZE
#  define PTR2EA(ptr,ea) {  *((void const**)(ea)) = (void const*)(ptr); }
#  define EA2PTR(ea)     (*((void**)(ea)))
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






/* Pointers to the begining of the arrays of buffer
   (to make the static arrays defined above publically available) */

/* cbe_env_begin is passed as environment address to every SPU thread */
extern envbuf_t* const cbe_env_begin;
/* cbe_arg_begin + N_ARGBUF*k is passed as argument address to SPU thread k */
extern argbuf_t* const cbe_arg_begin;


/* Additional work package buffer */
extern wp_t* const cbe_wp_begin;





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
int cbe_init(int const nspu_req, int const cpu_node);



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










/* Common part again. 
   "Generic" declarations follow which need declarations from PPU/SPU parts
 */







/* Some rounding functions */
static INLINE_ int iceil128(int const x) {
    return (x+127) & (~127);
}

static INLINE_ int iceil16(int const x) {
    return (x+15) & (~15);
}


static INLINE_ unsigned uiceil128(unsigned const x) {
    return (x+127u) & (~127u);
}

static INLINE_ unsigned uiceil16(unsigned const x) {
    return (x+15u) & (~15u);
}




#endif /* IMD_CBE_H */
