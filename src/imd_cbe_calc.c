#if defined (__PPU__)
/* Vector intrinsics on the PPU (AltiVec) */
#include <altivec.h>
#include <stdio.h>
#endif

#if defined(__SPU__)
/* Vector intrinsics on the SPU */
#include <spu_intrinsics.h>
/* #include <stdio.h> */
#endif


#include "imd_cbe.h"



/* --- Common part --- */

/* Macro for the calculation of the Lennard-Jones potential */
#define LJ(pot,grad,r2)  {                                                 \
   int const col4     = col<<2u;                                           \
   flt const tmp2     = pt.lj_sig[col4] / r2;                              \
   flt const tmp6     = tmp2 * tmp2 * tmp2;                                \
   flt const ptljtmp6 = pt.lj_eps[col4] * tmp6;                            \
   grad = ptljtmp6 * (tmp6-1.0) * (-12.0/r2);                          \
   pot  = ptljtmp6 * (tmp6-2.0)                 - pt.lj_shift[col4];   \
}



/* --- (Plain old) Scalar version --- */
static void calc_wp_scalar(wp_t* const wp)
{
  flt d[4], f[4], r2, *pi, *pj, *fi, *fj, pot, grad;
  int i, c1, col, inc=pt.ntypes * pt.ntypes, n, j;

  wp->totpot = wp->virial = 0.0;
  for (i=0; i<4*wp->n2; i++) { 
      wp->force[i] = 0.0;
  }

  /* fetch sizes and pointers */
  /* fetch data at wp->pos, wp->typ, wp->ti */

  for (i=0; i<wp->n1; i++) {

    short *ttb;

    /* fetch data at wp->tb + wp->ti[2*i] */
    ttb = wp->tb + wp->ti[2*i];

    c1 = pt.ntypes * wp->typ[i];
    fi = wp->force + 4*i;
    pi = wp->pos   + 4*i;

    for (n=0; n<wp->ti[2*i+1]; n++) {

      j    = ttb[n];
      pj   = wp->pos + 4*j;
      d[0] = pj[0] - pi[0];
      d[1] = pj[1] - pi[1];
      d[2] = pj[2] - pi[2];
      r2   = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
      col  = c1 + wp->typ[j];

      if (r2 <= pt.r2cut[4*col]) {

  	  /* Calculation of potential has been moved to a macro */
   	  LJ(pot,grad,r2)

          wp->totpot += pot;
          pot        *= 0.5;   /* avoid double counting */
          wp->virial -= r2  * grad;

          f[0] = d[0] * grad;
          f[1] = d[1] * grad;
          f[2] = d[2] * grad;

          fi[0] += f[0];
          fi[1] += f[1];
          fi[2] += f[2];
          fi[3] += pot;

          fj = wp->force + 4*j;
          fj[0] -= f[0];
          fj[1] -= f[1];
          fj[2] -= f[2];
          fj[3] += pot;

      }
    }
  }
  /* return wp->force, wp->totpot, wp->virial */
}











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

























/******************************************************************************
*
*  calc_wp
*
******************************************************************************/



#if defined(__PPU__)   /* --- PPU part --- */

/* (Trivial) mapping which makes some SPU intrinsics available as
   AltiVec intrinsics on the PPU */
#define spu_add   vec_add
#define spu_sub   vec_sub
#define spu_madd  vec_madd
#define spu_sel   vec_sel



/* Choose scalar or vector version */
#define CALC_DISPATCH_PPU   calc_wp_scalar

void calc_wp(wp_t* wp) {
     CALC_DISPATCH_PPU(wp);
}

#endif /* PPU */









#if defined(__SPU__)     /* -- SPU Part --- */

/* (Trivial) mappings which make some AltiVec intrinsics
   available in SPU code */
#define vec_add   spu_add
#define vec_sub   spu_sub
#define vec_madd  spu_madd
#define vec_sel   spu_sel




/* Sum over 1st three components of vector float v */
static INLINE_ float vf3sum(register VEC(float) const v) {
     return spu_extract(spu_add(spu_add(spu_rlqwbyte(v,4), v), 
                                spu_rlqwbyte(v,8)
                               ),  
                        0
                       );
}

static INLINE_ float vf3abs2(register VEC(float) const v) {
     return vf3sum(spu_mul(v,v));
}







static void calc_wp_vec1(wp_t* const wp)
{
   /* User short (16bit) as index */
   typedef unsigned short  index_t;

   /* The null vector initialized at compile (load?) time using a
      union */
   static VECUNION(float) const unull = { { 0, 0, 0, 0 } };


   /* Load some (vector) registers with some special vectors */
   register vecflt const vnull = unull.vec;


  /* Ptrs to positon & force vectors:
     Thses casts are well defined as both pointers in wp are correctly aligned
     (16-byte boundary as required by DMA).
     Posistions are constant (readonly), forces will be updated
   */
  register vecflt const* const pos   = (vecflt const*)(wp->pos);
  register vecflt *      const force = (vecflt*)(wp->force);
  register int const* const    typ   = wp->typ;

  /* Iterator for 1st particle */
  register index_t i;
  index_t const imax = wp->n1;


  /* Total negative virial, potential */
  register flt virial_, totpot;

  /* Set output variables (potential, virial & all forces) to zero  */
  virial_ = totpot = 0;
  for ( i=0;      i<(wp->n2);         ++i ) {
      force[i] = vnull;
  }


  /* For all particles in wp (1st particle) */
  for (  i=0;        i<imax;           ++i ) {

    /* Some pointers into tables/lists */
    int   const* const pti2i = (wp->ti)+(i<<1u);
    register short const* const ttb   = (wp->tb)+(*pti2i);
    index_t const c1 = pt.ntypes * typ[i];

    /* Index in NBL, maximum value of that index */
    register index_t n;
    index_t const nmax = *(pti2i+1);


    for ( n=0;        n<nmax;      ++n  ) {

      /* Index of 2nd particle is looked up in the NBL */
      register index_t const j    = ttb[n];
      register index_t const col  = c1+typ[j];

      /* Distance vector between particle i and particle j  (and vice versa) */
      register vecflt const di = vec_sub(pos[j],  pos[i]);
      register vecflt const dj = vec_sub(vnull, di);


      /* Distance squared */
      flt const r2 = vf3abs2(di);


      /* No interaction if distance is to large */
      if (  r2 <= pt.r2cut[col<<2u] ) {

          /* Gradient & potential */
          /* scalar values */
          flt grad, pot;
          /* vector containing these scalar values at some position */
          register vecflt vgrad;
          register vecflt vpot;
          /* New values for fi & fj */
          register vecflt newfi;
          register vecflt newfj;

  	  /*  Given r2, calculate Potential & Gradient */
   	  LJ(pot, grad, r2)
 
	  /* Update virial & potential */
          totpot   += pot;
          virial_  += r2  * grad;
          /* pot           *= 0.5f; */  /* avoid double counting */

          /* Generate vectors containing gradient & potential
             vgrad = (grad, grad, grad,       0)
             vpot  = (   0,    0,    0, 0.5*pot)
           */
          vpot  = spu_insert(0.5*pot,  vnull,            3);
          vgrad = spu_insert(0.0,      spu_splats(grad), 3);

          /* Update forces */
          newfi = vec_add(vec_madd(di, vgrad, vpot), force[i]);
	  newfj = vec_add(vec_madd(dj, vgrad, vpot), force[j]);
          /* Store updated forces */
          force[i]=newfi;
          force[j]=newfj;
      }
    }
  }
  /* return wp->force, wp->totpot, wp->virial */
  wp->virial = -virial_;
  wp->totpot = totpot;
}





#define CALC_DISPATCH_SPU   calc_wp_vec1

void calc_wp(wp_t* wp) {
     unsigned dt;
     dt = spu_readch(SPU_RdDec);
     CALC_DISPATCH_SPU(wp);
     dt -= spu_readch(SPU_RdDec);

     /* fprintf(stdout, "calc_wp on SPU took %u decrementer ticks.\n", dt); */
}

#endif  /* SPU */
