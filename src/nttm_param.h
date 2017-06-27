/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2007 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* nttm_parameter.h -- Material parameters for nttm simulations
*
******************************************************************************/

/******************************************************************************
* $Revision: 1.0 $
* $Date: 2013/09/25 17:43:09 $
******************************************************************************/

#define SILICON

#if defined(SILICON)

/* inverse photon energy [1/eV] cor. to 775 nm */
#define NTTM_HNU_INV 0.625081

/* two-photon  absorption coefficient [A^-1] (DOI:10.1007/s00339-011-6573-z) */
#define NTTM_BETA 0.0

/* room temperature 298 K in eV */
#define NTTM_T_ROOM 0.0256796

/* Auger recombination coefficient [A^6/time unit] (DOI:10.1007/s00339-011-6573-z) */
#define NTTM_GAMMA 3668.59

/* cutoff temperature for band gab = melting temperature */
#define NTTM_T_CUT 0.145374
/* cutoff start for band gap = T_cut - 500 K */
#define NTTM_T_0 0.1022875161

/* equillibrium carrier density [A^-3] */
#define NTTM_N(T) 1e-12//0.2070627924*pow((T)/CTTM_T_ROOM,1.5)*2.5e-5*exp(-0.5*1.16/(T))

/* carrier-phonon relaxation coefficient [1/time unit] (DOI:10.1007/s00339-011-6573-z) */
#define NTTM_TAU_INVERSE(N) 1.0/(23.5744855386*(1.0+(N)*(N)*2777777.777778))
#define NTTM_TAU(N) (23.5744855386*(1.0+(N)*(N)*2777777.777778))

/* cutoff function for energy gap */

#define NTTM_FC(T) ( (T < NTTM_T_0)?1.0:((T > NTTM_T_CUT)?0.0 :\
    0.5 * (1.0 + 1.125*cos(M_PI*(T-NTTM_T_0)/(NTTM_T_CUT-NTTM_T_0))-\
    0.125*cos(3.0*M_PI*(T-NTTM_T_0)/(NTTM_T_CUT-NTTM_T_0)))) )

#define NTTM_DFC(T) ( (T<NTTM_T_0 || T>NTTM_T_CUT)?1.0 :\
    -0.5*M_PI/(NTTM_T_CUT-NTTM_T_0)*(1.125*sin(M_PI*(T-NTTM_T_0)/(NTTM_T_CUT-NTTM_T_0))-\
    0.375*sin(3.0*M_PI*(T-NTTM_T_0)/(NTTM_T_CUT-NTTM_T_0))) )

//#define NTTM_FC(R,T)                                                      \
//{                                                                         \
//    if (T < NTTM_T_0) R = 1.0;                                            \
//    else if (T > NTTM_T_CUT) R = 0.0;                                     \
//    else{                                                                 \
//      double tmp = M_PI * (T - NTTM_T_0) / (NTTM_T_CUT - NTTM_T_0);       \
//      R = 0.5 * (1.0 + 1.125*cos(tmp) - 0.125*cos(3.0*tmp));              \
//    }                                                                     \
//}

//#define NTTM_DFC(R,T)                                                          \
//{                                                                              \
//    if (T < NTTM_T_0) R = 0.0;                                                 \
//    else if (T > NTTM_T_CUT) R = 0.0;                                          \
//    else{                                                                      \
//      double tmp = M_PI * (T - NTTM_T_0) / (NTTM_T_CUT - NTTM_T_0);            \
//      R = - 0.5 * M_PI / (NTTM_T_CUT - NTTM_T_0) * (1.125*sin(tmp) - 0.375*sin(3.0*tmp));\
//    }                                                                          \
//}

/* band-gap energy [eV] (DOI:10.1007/s00339-011-6573-z) */
#define NTTM_EGAP(T,N) NTTM_FC(T)*(1.16 - 8.14639*T*T/(T+0.093) - 1.5*pow(N,1.0/3.0))
#define NTTM_EGAP_DN(N) -0.5*NTTM_FC(T)*pow(N,-2.0/3.0)
#define NTTM_EGAP_DT(T,N) NTTM_FC(T)*8.14639*T/(T+0.093)*(T/(T+0.093)-2.)+\
    NTTM_DFC(T)*(1.16-8.14639*T*T/(T+0.093)-1.5*pow(N,1.0/3.0))


/* impact ionization coefficient [1/time unit] (DOI:10.1007/s00339-011-6573-z) */
#define NTTM_THETA(Eg,Tc) 0.0003665*exp(-1.5*(Eg)/(Tc))

/* carrier thermal conductivity coefficients: K(T) = A*T + B (DOI:10.1007/s00339-011-6573-z) */
//#define NTTM_K_A 6.1007777
//#define NTTM_K_B -0.04099459
#define NTTM_K_A 9.8462375711
#define NTTM_K_B -0.823026763812

#define NTTM_K_TCUT 0.09287549436115936
#define NTTM_K_BIGA NTTM_K_A/(10.0*pow(NTTM_K_TCUT,9.0))

#define NTTM_K(T) (T < NTTM_K_TCUT) ? NTTM_K_BIGA*pow(T,10.0) : NTTM_K_A*T+NTTM_K_B
#define NTTM_DK(T) (T < NTTM_K_TCUT) ? 10*NTTM_K_BIGA*pow(T,9.0) : NTTM_K_A

#define NTTM_EPS_RE(N) (13.794 - 889.347*N)
#define NTTM_EPS_IM(N) (0.0594 + 712.931*N)

#define NTTM_F_SQ(N) 0.5*(sqrt(NTTM_EPS_RE(N)*NTTM_EPS_RE(N) + NTTM_EPS_IM(N)*NTTM_EPS_IM(N)) + NTTM_EPS_RE(N))
#define NTTM_G_SQ(N) 0.5*(sqrt(NTTM_EPS_RE(N)*NTTM_EPS_RE(N) + NTTM_EPS_IM(N)*NTTM_EPS_IM(N)) - NTTM_EPS_RE(N))

#define NTTM_ALPHA_PRE 0.00166442
#define NTTM_ALPHA(N) NTTM_ALPHA_PRE*sqrt(NTTM_G_SQ(N))

#define NTTM_REFLECTIVITY(N) ((sqrt(NTTM_F_SQ(N))-1.)*(sqrt(NTTM_F_SQ(N))-1.)+NTTM_G_SQ(N))\
                             /((sqrt(NTTM_F_SQ(N))+1.)*(sqrt(NTTM_F_SQ(N))+1.)+NTTM_G_SQ(N))

/* free-carrier absorption cross-section [A^2] (DOI:10.1007/s00339-011-6573-z) */
#define NTTM_BIGTHETA(Tl) 0.051*(Tl)/NTTM_T_ROOM


#elif defined(GERMANIUM)

#endif


