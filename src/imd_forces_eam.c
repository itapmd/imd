
/*****************************************************************************
*
* imd_forces_eam.c -- force loop for Finnis-Sinclair EAM
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"


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
*
*  Modified (7/99):  calculation pot_k0 and pot_k1 (+ -> *)
*                    sort  = 0: only for EAM element type
*		     sort >= 1: only for element type described by pair 
*                               interaction
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
  real r2_short=cellsz;
  int jstart,jend;
  int column, is_short=0;
  int q_typ,p_typ;
  real *qptr;

  integer  eam_k,eam_pni,eam_pnj;     /* dummy, atom number */
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
  for (i=0; i < p->n; ++i) {

    eam_pni = p->nummer[i];
    if (eam_pni < 0 ) {
      eam_pni = -eam_pni;
    }
    /* Some compilers don't find the expressions that are invariant 
       to the inner loop. I'll have to define my own temp variables. */
    tmp_d.x = p->ort X(i) - pbc.x;
    tmp_d.y = p->ort Y(i) - pbc.y;
#ifndef TWOD
    tmp_d.z = p->ort Z(i) - pbc.z;
#endif
    p_typ   = SORTE(p,i);

    /* atom i interacts with j */
    jstart = (p==q ? i+1 : 0);

    for (j = jstart; j < q->n; ++j) {

      /* Calculate distance  */
      d.x      = q->ort X(j) - tmp_d.x;
      d.y      = q->ort Y(j) - tmp_d.y;
#ifndef TWOD
      d.z      = q->ort Z(j) - tmp_d.z;
#endif
      q_typ    = SORTE(q,j);
      column   = p_typ * ntypes + q_typ;
      radius2  = SPROD(d,d);
      eam_r_ij = sqrt(radius2);

#ifndef NODBG_DIST
      if (0==radius2) { char msgbuf[256];
        sprintf(msgbuf,"Pair distance is zero: i=%d (#%d), j=%d (#%d)\n",
                i,eam_pni,j,eam_pnj);
        error(msgbuf);
      }
#else
      if (0==radius2) error("Pair distance is zero.");
#endif

      /* 1. Cutoff: pair potential */
      if (radius2 <= pair_pot.end[column]) {

      	/* Check for distances, shorter than minimal distance in pot. table */
        if (1==pair_int2(&pot_zwi,&pot_grad,&pair_pot,column,radius2)) {
          r2_short = MIN(r2_short,radius2);
	  is_short=1; 
	}

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
        if ((NUMMER(p,i)>=0) || (NUMMER(q,j)>=0))
          tot_pot_energy += pot_zwi;

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
	p->presstens_offdia[i] -= d.x * d.y * pot_grad;
	q->presstens_offdia[j] -= d.x * d.y * pot_grad;
#else
	p->presstens Z(i) -= d.z * d.z * pot_grad;
	q->presstens Z(j) -= d.z * d.z * pot_grad;
	p->presstens_offdia X(i) -= d.y * d.z * pot_grad;
	p->presstens_offdia Y(i) -= d.z * d.x * pot_grad;
	q->presstens_offdia X(j) -= d.y * d.z * pot_grad;
	q->presstens_offdia Y(j) -= d.z * d.x * pot_grad;
	p->presstens_offdia Z(i) -= d.x * d.y * pot_grad;
	q->presstens_offdia Z(j) -= d.x * d.y * pot_grad;
#endif
#endif


      }; /* if 1. Cutoff */

      /* Calc of EAM cohesive function ONLY for atom sort "0" */
      /*                           ... not nice, but it works */

      if (p_typ == 0 && q_typ == 0) {
      
        /* 2. Cutoff: Calc of Finnis/Sinclair EAM cohesive function */

        if (radius2 <= eam_r2_cut) {

      	  /* Check for distances, shorter than minimal distance */
	  if (eam_r_ij <= eam_r_0) {
	    eam_r_ij = eam_r_0; 
  	  }; 

	  eam_pnj    = q->nummer[j];	
	  if (eam_pnj < 0 ) {
	    eam_pnj = -eam_pnj;
	  }
	  /* EAM-CF:  fixed i and j */
	  eam_phi_ij = (eam_r_ij-eam_r_cut)*(eam_r_ij-eam_r_cut);

          /* save #j in the #i array field k > 0  (-> calc of forces) */
          eam_ij[eam_pni*eam_len]         += (real) 1;
          eam_k                            = (integer) eam_ij[eam_pni*eam_len];
          if (eam_k >= eam_len) {
            error("neighbor table too small, increase eam_len");
	  }
          eam_ij[eam_pni*eam_len+eam_k]    = (real) eam_pnj; 
          eam_dij_x[eam_pni*eam_len+eam_k] = d.x; 
          eam_dij_y[eam_pni*eam_len+eam_k] = d.y; 
          eam_dij_z[eam_pni*eam_len+eam_k] = d.z; 
          /* EAM-CF: fixed i: sum of j's */
          eam_rho[eam_pni] += eam_phi_ij;

          /* save #i in the #j array field k > 0  (-> calc of forces) */
          eam_ij[eam_pnj*eam_len]         += (real) 1;
          eam_k                            = (integer) eam_ij[eam_pnj*eam_len];
          if (eam_k >= eam_len) {
            error("neighbor table too small, increase eam_len");
	  }
          eam_ij[eam_pnj*eam_len+eam_k]    = (real) eam_pni; 
          eam_dij_x[eam_pnj*eam_len+eam_k] = -d.x; 
          eam_dij_y[eam_pnj*eam_len+eam_k] = -d.y; 
          eam_dij_z[eam_pnj*eam_len+eam_k] = -d.z; 
          /* EAM-CF: fixed j: sum of i's */
          eam_rho[eam_pnj] += eam_phi_ij;

        } /* if 2. Cutoff */

      } /* p_typ, q_typ ==0 */

    } /* for j */

  } /* for i */

  /* A little security */
  if (is_short==1) printf("\n Short distance! r2: %f\n",r2_short);

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

} /* do_forces_eam_1 */

/* -------------------------------------------- */
void do_forces_eam_2(cell *p)

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
  int q_typ, p_typ;
  int jstart,jend;

  integer  eam_k;         /* dummy */
  integer  eam_pni;       /* dummy */
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

    p_typ      = SORTE(p,i);

    /* force only for atom sort "0" */
    /* ... not nice, but it works */
    if (p_typ == 0) {

      eam_pni    = p->nummer[i];
      if ( eam_pni < 0 ) {
        eam_pni = -eam_pni;
      }

      /* EAM-CF: sum of all i energies */
      eam_cf_i   = sqrt(eam_rho[eam_pni]);
      if (NUMMER(p,i)>0) eam_cf += eam_cf_i;

      eam_tmp_i  = 1.0/eam_cf_i;
      jstart     = 1;
      jend       = (integer) eam_ij[eam_pni*eam_len];

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
        eam_k         = (integer) eam_ij[eam_pni*eam_len+j];  
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
#ifdef TWOD
	p->presstens_offdia[i] -= d.x * d.y * eam_pot_grad;
	q->presstens_offdia[j] -= d.x * d.y * eam_pot_grad;
#else
        p->presstens Z(i) -= d.z * d.z * eam_pot_grad;
        q->presstens Z(j) -= d.z * d.z * eam_pot_grad;
	p->presstens_offdia X(i) -= d.y * d.z * eam_pot_grad;
	p->presstens_offdia Y(i) -= d.z * d.x * eam_pot_grad;
	q->presstens_offdia X(j) -= d.y * d.z * eam_pot_grad;
	q->presstens_offdia Y(j) -= d.z * d.x * eam_pot_grad;
	p->presstens_offdia Z(i) -= d.x * d.y * eam_pot_grad;
	q->presstens_offdia Z(j) -= d.x * d.y * eam_pot_grad;
#endif
#endif

      } /* for j */

    } /* p_typ == 0 */

  } /* for i */

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

} /* do_forces_eam_2 */

