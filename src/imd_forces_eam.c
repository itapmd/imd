/*****************************************************************************
*
* imd_forces -- force loop
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/* do forces, risc version; EAM - jh 5/98 */

/* ----------------------------------------------------------------------
*  EAM: Finnis/Sinclair Phil. Mag. A (1984) Vol 50 no 1 p 45 
*		  
*   tot_energy 	= 1/2 * sum_i^N      sum_{j.ne.i}^N   V_{i,j}
*		  - A * sum_i^N sqrt(sum_{j.ne.i}^N phi_{i,j})
*		= pair potential
*		  - density correction ("cohesive function") 
*		= sum_i^N [ 1/2         sum_{j.ne.i}^N   V_{i,j}
*		             - A * sqrt(sum_{j.ne.i}^N phi_{i,j}) ]
*   The cutoffs for the pair potential and for the density correction
*   are neither identical nor do they exhibit a constant relationship
*   between (r2_cut is not always smaller than eam_r2_cut).
*   Therefore, two 'if' statements are needed; one for each cutoff        
*
*   do_forces_eam_1: calculation of the pair potential and its derivative
*   		     and save the eam neighbors
*   do_forces_eam_2: calculation of the cohesive potential function
*   		     and calculation of the derivative of the cohesive 
*   	             potential function
* -------------------------------------------------------------------- */

/* -------------------------------------------- */
void do_forces_eam_1(cell *p, cell *q, vektor pbc)

/* Part 1: calc of pair potential and its derivative 
*  and save the eam neighbors */

{
  int i,j,k;
  vektor d;
  vektor tmp_d;
  vektor force;
  real radius2;
  vektor tmp_forces;
  real tmp_energy;
  real tmp_virial;
#ifdef P_AXIAL
  vektor tmp_vir_vect;
#endif
  real pot_zwi;
  real pot_grad;
  int jstart,jend;
  int c;
#ifdef MONOLJ
  real sig_d_rad2,sig_d_rad6,sig_d_rad12;
#else
  real r2_short = r2_end;
  real pot_k0,pot_k1,pot_k2;
  real dv, d2v;
  int q_typ,p_typ;
  real *potptr;
  real chi;
#endif
  real *qptr;

  int  eam_k,eam_pni,eam_pnj;     /* dummy, atom number */
  real eam_phi_ij;        /* interaction between i and j */
  real eam_r_ij;          /* distance, needed for eam forces */

  tmp_virial     = 0.0;
#ifdef P_AXIAL
  tmp_vir_vect.x = 0.0;
  tmp_vir_vect.y = 0.0;
#ifndef TWOD
  tmp_vir_vect.z = 0.0;
#endif
#endif

  /* For each atom in first cell */
  for (i = 0;i < p->n; ++i) {

    eam_pni = p->nummer[i];

    /* Some compilers don't find the expressions that are invariant 
       to the inner loop. I'll have to define my own temp variables. */
    tmp_d.x = p->ort X(i) - pbc.x;
    tmp_d.y = p->ort Y(i) - pbc.y;
#ifndef TWOD
    tmp_d.z = p->ort Z(i) - pbc.z;
#endif
#ifndef MONOLJ
    p_typ   = p->sorte[i];
#endif
    
    /* atom i interacts with j */
    jstart = (p==q ? i+1 : 0);

    for (j = jstart; j < q->n; ++j) {

      /* Calculate distance  */
      d.x      = q->ort X(j) - tmp_d.x;
      d.y      = q->ort Y(j) - tmp_d.y;
#ifndef TWOD
      d.z      = q->ort Z(j) - tmp_d.z;
#endif
      radius2  = SPROD(d,d);
      eam_r_ij = sqrt(radius2);

 #ifndef NODBG_DIST
      if (0==radius2) { char msgbuf[256];
        sprintf(msgbuf,"Pair distance is zero: i=%d, j=%d\n",i,j);
        error(msgbuf);
      }
#else
      if (0==radius2) error("Pair distance is zero.");
#endif

      /* 1. Cutoff: Calc of pair interaction and forces */
      if (radius2 <= r2_cut) {

      	/* Check for distances, shorter than minimal distance in pot. table */
	if (radius2 <= r2_0) {
	  radius2 = r2_0; 
	}; 

	/* Indices into potential table */
	k      = (int) ((radius2 - r2_0) * inv_r2_step);
	q_typ  = q->sorte[j];
	/* Abweichung von k (k ist int, chi ist real);
	   entspricht dem Ausdruck (x-x0) in der Taylorreihe */
	chi    = (radius2 - r2_0 - k * r2_step) * inv_r2_step;
	
	/* A single access to the potential table involves two multiplications 
	   We use an intermediate pointer to avoid this as much as possible.
	   Note: This relies on layout of the pot-table in memory!!! */
	/* k0, k1, k2: Potentialwerte an drei aufeinanderfolgenden k */
 	potptr = PTR_3D_V(potential, k, p_typ, q_typ , pot_dim);
	pot_k0 = *potptr; potptr += pot_dim.y + pot_dim.z;
	pot_k1 = *potptr; potptr += pot_dim.y + pot_dim.z;
	pot_k2 = *potptr;

	/* dv: 	1. Abl. des Pots in Einheiten der Einheitsschrittweite k
	   d2v:	2. Abl. des Pots in Einheiten der Einheitsschrittweite k */
	dv     = pot_k1 - pot_k0;
	d2v    = pot_k2 - 2 * pot_k1 + pot_k0;

  	/* - Calculation of the total energy 'tot_pot_energy'
	   - Norm of Gradient: in x direction (y,z analog)
	 	d
	   	-- pot_zwi  = x * pot_grad
	       	dx           			            
	     (pot_grad includes negative sign!) 				*/ 
	pot_grad = 2 * inv_r2_step * ( dv + (chi - 0.5) * d2v );
        /* Potential energy of atom: 
	   Taylorreihe um k0 mit Restglied -1/2 * chi * d2v */
	pot_zwi  =  pot_k0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;
        
	/* Store forces in temp */
	force.x  = d.x * pot_grad;
	force.y  = d.y * pot_grad;
#ifndef TWOD
	force.z  = d.z * pot_grad;
#endif

        /* Accumulate forces  */
	p->kraft X(i)  += force.x;
	p->kraft Y(i)  += force.y;
	q->kraft X(j)  -= force.x;
	q->kraft Y(j)  -= force.y;
#ifndef TWOD
	p->kraft Z(i)  += force.z;
	q->kraft Z(j)  -= force.z;
#endif
#ifndef MONOLJ
	q->pot_eng[j]  += pot_zwi;
	p->pot_eng[i]  += pot_zwi;
#endif
	tot_pot_energy += 2*pot_zwi;

#ifdef P_AXIAL
        tmp_vir_vect.x -= d.x * d.x * pot_grad;
        tmp_vir_vect.y -= d.y * d.y * pot_grad;
#ifndef TWOD
        tmp_vir_vect.z -= d.z * d.z * pot_grad;
#endif
#else
        tmp_virial     -= radius2 * pot_grad;  
#endif
	/* negativ, da pot_grad gleich force !! */
#ifdef STRESS_TENS
	p->presstens X(i) -= d.x * d.x * pot_grad;
	p->presstens Y(i) -= d.y * d.y * pot_grad;
	q->presstens X(j) -= d.x * d.x * pot_grad;
	q->presstens Y(j) -= d.y * d.y * pot_grad;
#ifndef TWOD
	p->presstens Z(i) -= d.z * d.z * pot_grad;
	q->presstens Z(j) -= d.z * d.z * pot_grad;
#endif
#endif


      }; /* if 1. Cutoff */

      /* 2. Cutoff: Calc of Finnis/Sinclair EAM cohesive function */

      if (radius2 <= eam_r2_cut) {

      	/* Check for distances, shorter than minimal distance */
	if (eam_r_ij <= eam_r_0) {
	  eam_r_ij = eam_r_0; 
	}; 

	eam_pnj    = q->nummer[j];	
	/* EAM-CF:  fixed i and j */
	eam_phi_ij = (eam_r_ij-eam_r_cut)*(eam_r_ij-eam_r_cut);

        /* save #j in the #i array field k > 0  (-> calc of forces) */
        eam_ij[eam_pni*eam_len]         += 1;
        eam_k                            = eam_ij[eam_pni*eam_len];
        eam_ij[eam_pni*eam_len+eam_k]    = q->nummer[j]; 
        eam_dij_x[eam_pni*eam_len+eam_k] = d.x; 
        eam_dij_y[eam_pni*eam_len+eam_k] = d.y; 
        eam_dij_z[eam_pni*eam_len+eam_k] = d.z; 
        /* EAM-CF: fixed i: sum of j's */
        eam_rho[eam_pni] += eam_phi_ij;

        /* save #i in the #j array field k > 0  (-> calc of forces) */
        eam_ij[eam_pnj*eam_len]         += 1;
        eam_k                            = eam_ij[eam_pnj*eam_len];
        eam_ij[eam_pnj*eam_len+eam_k]    = p->nummer[i]; 
        eam_dij_x[eam_pnj*eam_len+eam_k] = -d.x; 
        eam_dij_y[eam_pnj*eam_len+eam_k] = -d.y; 
        eam_dij_z[eam_pnj*eam_len+eam_k] = -d.z; 
        /* EAM-CF: fixed j: sum of i's */
        eam_rho[eam_pnj] += eam_phi_ij;

      }; /* if 2. Cutoff */

    }; /* for j */

  }; /* for i */

#ifndef MONOLJ
  /* A little security */
  if (r2_short < r2_0) 
    printf("\n Short distance! r2: %f\n",r2_short);
#endif

#ifdef P_AXIAL
  vir_vect.x += tmp_vir_vect.x;
  vir_vect.y += tmp_vir_vect.y;
  virial     += tmp_vir_vect.x;
  virial     += tmp_vir_vect.y;
#ifndef TWOD
  vir_vect.z += tmp_vir_vect.z;
  virial     += tmp_vir_vect.z;
#endif
#else
  virial     += tmp_virial;
#endif

  return;
} /* do_forces_eam_1 */

/* -------------------------------------------- */
void do_forces_eam_2(cell *p, cell *q, vektor pbc)

/* Part 2: calc of cohesive potential function 
*
*  Part 3: The cohesive function force on atom i is the sum of:
*   - the change in energy due to change in density at the site i itself 
*     as atom i moves (force of eam_rho[i])
*   and
*   -  change in energy due to change in density at all the other sites  
*      as atom i moves (force of \sum_j eam_rho[j]; rij < eam_r_cut)	 */

{
  int i,j,k;
  vektor tmp_d;
  vektor d;
  vektor force;
#ifdef P_AXIAL
  vektor tmp_vir_vect;
#endif
  real tmp_virial;
  real radius2;
  real pot_zwi;
  real pot_grad;
  int q_typ, p_typ;
  int jstart,jend;

  int  eam_k;             /* dummy */
  int  eam_pni;           /* dummy */
  real eam_r_ik;          /* distance, needed for eam forces */
  real eam_tmp_i;         /* dummy */
  real eam_tmp_k;         /* dummy */
  real eam_pot_grad;      /* gradient des potentials */
  real eam_cf;            /* cohesive function */
  real eam_cf_i;          /* dummy */

  eam_cf = 0.;

  tmp_virial     = 0.0;
#ifdef P_AXIAL
  tmp_vir_vect.x = 0.0;
  tmp_vir_vect.y = 0.0;
#ifndef TWOD
  tmp_vir_vect.z = 0.0;
#endif
#endif

  /* For each atom in first cell: calculate force */
  for (i = 0;i < p->n; ++i) {

    eam_pni    = p->nummer[i];

    /* EAM-CF: sum of all i energies */
    eam_cf_i   = sqrt(eam_rho[eam_pni]);
    eam_cf    += eam_cf_i;

    eam_tmp_i  = 1.0/eam_cf_i;
    jstart     = 1;
    jend       = eam_ij[eam_pni*eam_len];

    /* pair potential and cohesive function for atom i */
    p->pot_eng[i] -= eam_A * eam_cf_i; 

    /* interaction of i with selected j's (called: eam_k)
       (selection in first for statement of i)            */

    for (j = jstart; j <= jend; ++j) {

      /* distance */
      d.x      = eam_dij_x[eam_pni*eam_len+j]; 
      d.y      = eam_dij_y[eam_pni*eam_len+j]; 
      d.z      = eam_dij_z[eam_pni*eam_len+j]; 
      radius2  = SPROD(d,d);
      if (0==radius2) error("EAM-Force Distance is zero.");
      eam_r_ik = sqrt(radius2);

      /* number of atom k */
      eam_k    	    = eam_ij[eam_pni*eam_len+j];  
      eam_tmp_k     = 1.0/sqrt(eam_rho[eam_k]);
      eam_pot_grad  = -eam_A*(eam_tmp_i+eam_tmp_k)*(eam_r_ik-eam_r_cut)/eam_r_ik;
      
      /* store forces in temp */
      force.x       = d.x * eam_pot_grad;
      force.y       = d.y * eam_pot_grad;
      force.z       = d.z * eam_pot_grad;

      /* Accumulate forces */
      p->kraft X(i) += force.x;
      p->kraft Y(i) += force.y;
      p->kraft Z(i) += force.z;

#ifdef P_AXIAL
      tmp_vir_vect.x -= d.x * d.x * eam_pot_grad;
      tmp_vir_vect.y -= d.y * d.y * eam_pot_grad;
#ifndef TWOD
      tmp_vir_vect.z -= d.z * d.z * eam_pot_grad;
#endif
#else
      tmp_virial     -= radius2 * eam_pot_grad;  
#endif
	/* negativ, da pot_grad gleich force !! */
#ifdef STRESS_TENS
      p->presstens X(i) -= d.x * d.x * eam_pot_grad;
      p->presstens Y(i) -= d.y * d.y * eam_pot_grad;
      q->presstens X(j) -= d.x * d.x * eam_pot_grad;
      q->presstens Y(j) -= d.y * d.y * eam_pot_grad;
#ifndef TWOD
      p->presstens Z(i) -= d.z * d.z * eam_pot_grad;
      q->presstens Z(j) -= d.z * d.z * eam_pot_grad;
#endif
#endif


    }; /* for j */

  }; /* for i */

  /* total energy: pair potential + cohesive function */
  tot_pot_energy  -= eam_A * eam_cf;

#ifdef P_AXIAL
  vir_vect.x += tmp_vir_vect.x;
  vir_vect.y += tmp_vir_vect.y;
  virial     += tmp_vir_vect.x;
  virial     += tmp_vir_vect.y;
#ifndef TWOD
  vir_vect.z += tmp_vir_vect.z;
  virial     += tmp_vir_vect.z;
#endif
#else
  virial     += tmp_virial;
#endif

  return;
} /* do_forces_eam_2 */

