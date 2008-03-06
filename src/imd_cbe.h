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
* imd_cbe.h -- header file for CBE
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#ifndef IMD_CBE_H__
#define IMD_CBE_H__

#include "config.h"

#ifdef __SPU__

/* SPU specific headers */
#include <spu_mfcio.h>
#include <spu_intrinsics.h>

#endif

#ifdef __PPU__

/* WORDSIZE & size_t */
#include <limits.h>
#include <stddef.h>

/* Some PPU macros are needed to implement other macros */
#include <ppu_intrinsics.h>

/* SPU runtime management */
#include <libspe2.h>

#endif


#if defined (__cplusplus)
extern "C" {
#endif


/********* Common (SPU & PPU) part +++++++++*/


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



/* Symbolic constants for the flag field */
enum { DOTB=((int)1), DOWP=((int)2) };

/* Cell data within a work package;
   On SPU, allocate for n_max atoms and len_max neighbor indices,
   namely 4*n_max floats (pos, force, imp), n_max ints (typ), 
   2*n_max ints (ti), and len_max shorts (tb), all rounded up to 128 bytes.
   n_max and len_max are taken from wp_t; 
   n is set to the actual n from cell_ea_t. */
typedef struct {
  int   n;  /* actual number of atoms */
  float  *force;
  float  *pos;
#ifdef SPU_INT
  float  *imp;
#endif
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
#ifdef SPU_INT
  ea_t imp_ea;
#endif
} cell_ea_t;


/* The work package type; essentially consists of a list of cells. 
   k       = package number, 
   n_max   = upper bound on number of atoms in a cell,
   len_max = upper bound on number of neighbor indices in a cell */
typedef struct wp {
  float totpot;
  float virial;
#ifdef SPU_INT
  float totkin;
#ifdef NVT
  float E_kin_2;
#endif
#ifdef FNORM
  float fnorm;
#endif
#ifdef GLOK
  float PxF, pnorm;
#endif
#endif
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



#ifdef LJ

/* potential type for Lennard-Jones */
typedef struct {
  int ntypes;
  ea_t ea;
#ifdef SPU_INT
  ea_t mv_ea;
#endif
#if defined(__SPU__)
  vector float const* r2cut;
  vector float const* lj_sig;
  vector float const* lj_eps;
  vector float const* lj_shift;
#else
  float *r2cut;
  float *lj_sig;
  float *lj_eps;
  float *lj_shift;
#endif
} pt_t;

#else

/* potential type for tabulated pair potentials */
typedef struct {
  int   ntypes;      /* number of atom types */
  int   nsteps;      /* number of tabulations steps */
  ea_t  ea;          /* effective address of data section */
#ifdef SPU_INT
  ea_t  mv_ea;       /* address of extra data for integrator */
#endif
#ifdef __SPU__
  vector float const *begin;    /* first value in the table */
  vector float const *r2cut;    /* last value in the table */
  vector float const *step;     /* table increment */
  vector float const *invstep;  /* inverse of increment */
  vector float const *data;     /* the actual data */
  vector unsigned int tab[4];   /* where columns start (max. 2x2 columns) */
#else
  float *begin;      /* first value in the table */
  float *r2cut;      /* last value in the table */
  float *step;       /* table increment */
  float *invstep;    /* inverse of increment */
  float *data;       /* the actual data */
  float *tab[4];     /* pointers to tabulation data (max. 2x2 columns) */
#endif
} pt_t;

#endif

/* This is defined in spu.c and imd_forces_cbe.c */
extern  pt_t pt;

#ifdef SPU_INT

/* extra data for integrator on SPU */
/* we assume here ntypes <= 2 and total_types <= 6 */
typedef struct {
#ifdef __SPU__
  vector float ts, fac1, fac2, mikmask, imass[2], restr[6], fbc_f[6];
#else
  float ts[4], fac1[4], fac2[4], mikmask[4], imass[8], restr[24], fbc_f[24];
#endif
} mv_t;

/* This is defined in spu.c and imd_forces_cbe.c */
extern  mv_t mv;

#endif

/* summation WP */
extern wp_t psum;

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


/* data structure for memory buffer */
typedef struct {
  unsigned int  len;      /* size of buffer in bytes */
  unsigned char *data;    /* data section of buffer  */
} mem_buf_t;


typedef union {
    /* Dummy data */

    /* Arguments/work packages */
    wp_t   w;

    /* Or simply an EA list */
    ea_t ea[sizeof(wp_t)/sizeof(ea_t)];

    /* Padding */
    unsigned char pad[CEILPOW2(sizeof(wp_t),BUFPAD)];

} argbuf_t;


/* The main calculation routine(s) */
void calc_wp(wp_t *wp);
void calc_tb(wp_t *wp);

void calc_wp_direct(wp_t*,
                    /* void* const, unsigned const, void* const, unsigned const, */
                    cell_dta_t* const pbuffers,
                    unsigned const otag);

void calc_tb_direct(wp_t*,
                    /* void* const, unsigned const, void* const, unsigned const, */
                    cell_dta_t* const pbuffers,
                    unsigned const otag);



#if defined(__SPU__)          /************** SPU part ************/

/* Unsigned integer type which may hold a pointer 
  (and vice versa)*/
typedef unsigned int  uintaddr_t;

INLINE_ unsigned wait_dma(unsigned const m, unsigned const type);
INLINE_ void dma64(register void* const p,
    register ea_t const ea, register unsigned int const size,
    register unsigned int const tag, register unsigned int const cmd);
INLINE_ void mdma64(void* const p, ea_t const ea, unsigned const size,
                           unsigned const tag, unsigned const cmd);

/* Generic macro for rounding up vector components in terms of
   SPU vector commands.
   On the SPU, x may be a scalar which is converted to a vector first
 */
#define CEILV(x,r) spu_andc(spu_add((r),(x)), (r))

#endif /* SPU */





#if defined(__PPU__)           /************ PPU part **************/


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



/* Unsigned integer type which may be used to represent a pointer
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
/* Get time base frequency from OS (see CBE programming handbook) */
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


/* The SPU work scheduling routines */
void do_work_spu(int const flag);


/* Multibuffered version of do_work_spu */
void do_work_spu_mbuf(void (* const mkf)(argbuf_t*, unsigned, int),
                      void (* const stf)(argbuf_t*, unsigned)
                     );


#endif /* PPU part */


#if defined (__cplusplus)
}
#endif


#endif /* IMD_CBE_H */
