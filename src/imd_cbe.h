#ifndef IMD_CBE_H__
#define IMD_CBE_H__

/* We need config here, as we have to test the CBE_DIRECT in the 
   DMA typedefs */
#include "config.h"


/********* Common (SPU & PPU) part +++++++++*/

#if defined (__cplusplus)
extern "C" {
#endif



/* ticks() is a qunatity of tick_t, a type which may different on PPU and SPU.
   ticks() changes at the same frequency on both the PPU and the SPU,
   however it is unspecified wether they count upwards or downwords.
   tick64_t is 64 bits wide on both the PPU and the SPU. A ticks_t may be
   converted to a tick64_t, but not vice versa. Both types are
   unsigned integer types.
   tick_diff(x,y) always returns the time (in ticks) elapsed between
   the statements x=ticks() and a following y=ticks().
 */
#if defined(__SPU__)
#define set_ticks(m) (spu_writech(SPU_WrDec, (m)))
#define ticks() (spu_readch(SPU_RdDec))
#define tick_diff(tfst,tsnd) ((tfst)-(tsnd))  /* tsnd<=tfst on SPU (where time ticks backwards) */
typedef unsigned int       tick_t;
typedef unsigned long long tick64_t;
#endif
#if defined(__PPU__)
#define set_ticks(dummy) /* This is ignored on the SPU, where the tb req must not be written by user software */
#define ticks __mftb
#define tick_diff(tfst,tsnd) ((tsnd)-(tfst))  /* tsnd>=tfist on PPU  */
typedef unsigned long long  tick_t, tick64_t;
#endif


/* Get start timestamp (t0), end timestamp (t1) of block cmd using
   function T */
#define TSTMPS(cmd, t0,t1, T) { (t0)=(T()); { cmd } (t1)=(T()); }

/* Run command cmd as a block,
   getting initial timestamp t0, and the time dt it took to execute.
   Functor T returns current time, Functor D calculates time difference
 */
#define TDIFF(cmd, t0,dt,  T,D) { TSTMPS(cmd, t0,dt, T)  (dt)=(D(t0,dt)); }






/* Inlining */
#if defined(__cplusplus)  /* C++ */
#    define INLINE_ inline
#else                     /* C */
#    if defined (__GNUC__) /* gcc */
#        define INLINE_ __inline
#    elif (defined(__IBMC__) || defined(__IBMCPP__)) /* xlc */
#        if defined(__IBM_GCC_INLINE)
#            define INLINE_ __inline__
#        else
#            define INLINE_ __inline
#        endif
#    else
#        error "INLINE_ macro empty"
#        define INLINE_
#    endif
#endif


/* The ALIGNED_(b,v) macro is used to align
   varibale v to a boundary specified by b.
   E.g. Use it to specify 16-byte alignment (which is needed for DMA)
   as follows:   buffer_type ALIGNED_(16, buf);
 */
#if defined(__GNUC__)  /* gcc */
#   define ALIGNED_(alg,var) var __attribute__((aligned(alg)))  
#elif (defined(__IBMC__) || defined(__IBMCPP__))  /* xlc */
#   define ALIGNED_(alg,var) var __attribute__((aligned(alg)))
#else
#   /* Just define ALIGNED to be the variable */
#   /* a NOOP */
#   error "ALIGNED_ macro was defined, but may not align objects."
#   define ALIGNED_(alg,var) var
#endif




/* Some branch prediction macros */
#if defined (__GNUC__) /* gcc */
#    define EXPECT_(x,v)  __builtin_expect((x), (v))
#elif (defined(__IBMC__) || defined(__IBMCPP__) )  /* xlc */
#    define EXPECT_(x,v)  __builtin_expect((x), (v))
#else
#    error "EXPECT macro was defined but will not have branch predicting effect"
#    define EXPECT_(x,v)  (x)
#endif


#define EXPECT_TRUE_(x)   EXPECT_((x), 1)
#define EXPECT_FALSE_(x)  EXPECT_((x), 0)






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















/* Effective address type: unsigned long long is 64 bits wide on
   the SPU and on the PPU (both in 32bit & 64bit mode). */
typedef unsigned long long int  ea_t;



/* A reservation granule is 128 bytes large on the CBE (which is also the PPE
   cache line size).
   According to the PowerPC architecture, this size is implementation dependent
   but must be pow(2,n) with n>=4.
   An array unsigned/signed int[CBE_GRANULE_NINT] is of that size.
 */
enum { CBE_GRANULE_NBYTE = ((unsigned)128) };
enum { CBE_GRANULE_NINT  = (unsigned)(CBE_GRANULE_NBYTE/(sizeof (int))) };


/* The floating point type */
typedef float flt;








#ifdef CBE_DIRECT  /* "Direct " */

/* Symbolic constants for the flag field */
enum { DOTB=((int)1), DOWP=((int)2) };

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

  /* Some timing values */
  tick_t timings[1];
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


    /* Padding, no longer needed */
    /* unsigned char pad[4]; */
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
} env_t;





/* Constants passed via mailbox */
enum { WPDONE1=(unsigned)0, WPSTRT1=(unsigned)1, SPUEXIT=(unsigned)2 };


/* There shall be N_ARGBUF arguement buffers for each SPU thread.
   All SPU threads shall receive a private argument address which does
   not equal an argument address passed to an other SPU.

   There shall be N_ENVBUF environment buffers which are shared by all
   SPU threads. All SPU shall receive the same environment address
 */
enum { /* Max. number of argument buffers per SPU */
       N_ARGBUF=((unsigned)32),
       /* Max. number of env. buffers */
       N_ENVBUF=((unsigned)2),
       /* Max number of EAs in env. list */
       N_ENVEA=(unsigned)(128u/(sizeof(ea_t)))
};

/* Buffer padding (must be a multiple of 16, should be a multiple 128) */
enum { BUFPAD = ((unsigned)128) };



/* x padded to the next boundary... */
/* ...given a mask of some lower 1-bits */
#define CEILONES(x,m) (((x)+(m)) & (~(m)))
/* ...given a power of two */
#define CEILPOW2(x,p) CEILONES((x), ((p)-1))
/* ...given a bit postion e */
#define CEILBPOS(x,e) CEILPOW2((x), (1<<(e)))


/* Maximum size of the two types a and b */
#define MAXSIZE2(a,b) (((sizeof(a))>(sizeof(b))) ? (sizeof(a)) : (sizeof(b)))


typedef union {
    /* Dummy data */

    /* Arguments/work packages */
    wp_t   w;
    exch_t e;

    /* Or simply an EA list */
    ea_t ea[(MAXSIZE2(wp_t,exch_t))/(sizeof (ea_t))];

    /* Padding */
    unsigned char pad[CEILPOW2(MAXSIZE2(wp_t,exch_t),BUFPAD)];
} argbuf_t;


typedef union {
    /* Dummy data: */

    /* Potential type */
    pt_t  p;
    env_t e;

    /* Or simply an EA list */
    ea_t ea[(MAXSIZE2(pt_t,env_t))/(sizeof (ea_t))];

    /* Padding */
    unsigned char pad[CEILPOW2(MAXSIZE2(pt_t,env_t),BUFPAD)];
} envbuf_t;














/* The main calculation routine(s) */
void calc_wp(wp_t *wp);
void calc_tb(wp_t *wp);

#if defined(CBE_DIRECT)
void calc_wp_direct(wp_t*,
                    /* void* const, unsigned const, void* const, unsigned const, */
                    cell_dta_t* const pbuffers,
                    unsigned const otag);

void calc_tb_direct(wp_t*,
                    /* void* const, unsigned const, void* const, unsigned const, */
                    cell_dta_t* const pbuffers,
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



/* Unsigned integer type which may hold a pointer 
  (and vice versa)*/
typedef unsigned int  uintaddr_t;



/* Wait for DMA specified by mask m using type for update */
static INLINE_ unsigned wait_dma(unsigned const m, unsigned const type) {
    spu_writech(MFC_WrTagMask,   m);
    spu_writech(MFC_WrTagUpdate, type);
    return spu_readch(MFC_RdTagStat);
}





/* (Just a) wrapper for DMA with less than 16K */
static INLINE_ void (dma64)(register void* const p,
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
static INLINE_ void mdma64(void* const p, ea_t const ea, unsigned const size,
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











/* Generic macro for rounding up vector components in terms of
   SPU vector commands.
   On the SPU, x may be a scalar which is converted to a vector first
 */
#define CEILV(x,r) spu_andc(spu_add((r),(x)), (r))


static INLINE_
vector unsigned int uiceilv(vector unsigned int const v, vector unsigned int const m)
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

/* Some PPU macros are needed to implement other macros */
#include <ppu_intrinsics.h>

/* SPU runtime management */
#include <libspe2.h>




#if defined (__cplusplus)
extern "C" {
#endif



/* Map some vector commands */
#define spu_andc  vec_andc
#define spu_add   vec_add



/* Pointer size in bits (defined as a macro such that it may be
   used further down in an #if directive) */
#if defined(__GNUC__)  /* gcc */
#   if (64 == __WORDSIZE)
#       define PPU_PTRBITS_ 64
#   elif (32 == __WORDSIZE)
#      define PPU_PTRBITS_ 32
#   else
#      error "Can't determine pointer size"
#   endif
#elif (defined(__IBMC__) || defined(__IBMCPP__))  /* xlc */
#   if defined(__64bit__)
#        define PPU_PTRBITS_ (8*8) 
#   else
#        define PPU_PTRBITS_ (4*8)
#   endif
#else
#    error "Can't determine pointer size"
#endif

/* The same value as above as a compile time constant */
enum { PPU_PTRBITS = (unsigned)(PPU_PTRBITS_) };






/* Unsigned integer type which may used to represent a pointer
   (and vice versa). This type is not larger than the ea_t.
*/
typedef
#if (32 == PPU_PTRBITS_)
    unsigned int
#elif (64 == PPU_PTRBITS_)
    ea_t
#endif
uintaddr_t;

/* Conversion from pointer to effective address is just a cast */
#define EA_CAST(ptr)   ((uintaddr_t)(ptr))
#define PTR2EA(ptr,ea) { (void)((ea)=EA_CAST(ptr)); }




/* Max. number of SPU threads which may be managed by the
   following utility functions   */
enum {
   N_SPU_THREADS_MAX =  (unsigned)
#if defined (SPU_THREADS_MAX)
      (SPU_THREADS_MAX) /* Use the #defined value */
#else
      (32) /* Use a default which should be large enough */
#endif
};






/* Pointers to the begining of the arrays of buffer
   (to make the static arrays defined above publically available) */

/* cbe_env_begin is passed as environment address to every SPU thread */
/* extern envbuf_t* const cbe_env_begin; */

/* cbe_envea_begin is passed as environment address to every SPU thread 
   and points to an array of EAs ("pointers") */
extern ea_t* const cbe_envea_begin;

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
*/
int cbe_init(int const nspu_req, int const cpu_node);



/* CBE specific cleanup */
void cbe_shutdown(void);







/* CBE time base utilities */


/* Get timebase frequency from contents of cpuinfo file */
unsigned long tbfreq_string(char const* const sfrst, char const* const slast);

/* Read cpuinfo from fd into buffer buf of length bufsze parse it and return
   the timebase freqeuncy. */
unsigned long tbfreq_fd(int const fd,  char* const buf, unsigned const bufsze);

/* Get time base frequency from OS (as suggested by the CBE programming handbook) */
unsigned long tbfreq(void);






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




/* Output of pt_t *e using  printf-like function of */
int (pt_out)(int (*of)(char const[],...), pt_t const* const e);

/* Output of work package wp using printf-like functions of */
int (wp_out)(int (*of)(char const[],...), wp_t const* const e);





/* The SPU work scheduling routines */
void do_work_spu(int const flag);


#if defined(CBE_DIRECT)
/* Multibuffered version of do_work_spu */
void do_work_spu_mbuf(void (* const mkf)(argbuf_t*, unsigned, int),
                      void (* const stf)(argbuf_t*, unsigned)
                     );
#endif



#if defined (__cplusplus)
}
#endif


#endif /* PPU part */





/* Common part again. 
   "Generic" declarations follow which need declarations from PPU/SPU parts
 */







#endif /* IMD_CBE_H */
