
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2004 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
*  do_forces for ASYMPOT, and second force loop for EAM2
*
******************************************************************************/

/******************************************************************************
*  $Revision$
*  $Date$
******************************************************************************/

#include "imd.h"
#include "potaccess.h"

/******************************************************************************
*
*  special version of do_forces for asymmetric core potentials
*
******************************************************************************/

#ifdef ASYMPOT

void do_forces(cell *p, cell *q, vektor pbc, real *Epot, real *Virial, 
               real *Vir_xx, real *Vir_yy, real *Vir_zz,
               real *Vir_yz, real *Vir_zx, real *Vir_xy)
{
  int i,j,k;
  vektor d;
  vektor tmp_d;
  vektor force;
  real r2, rho_h;
  real tmp_virial;
#ifdef P_AXIAL
  vektor tmp_vir_vect;
#endif
  real pot_zwi, pot_grad;
  int col1, col2, is_short=0, inc = ntypes * ntypes;
  int jstart, q_typ, p_typ;
  real *qptr, *pfptr, *qfptr, *qpdptr, *ppdptr, *qpoptr, *ppoptr;
  
  tmp_virial     = 0.0;
#ifdef P_AXIAL
  tmp_vir_vect.x = 0.0;
  tmp_vir_vect.y = 0.0;
  tmp_vir_vect.z = 0.0;
#endif
    
  /* for each atom in first cell */
  for (i=0; i<p->n; ++i) {

    tmp_d.x = ORT(p,i,X) - pbc.x;
    tmp_d.y = ORT(p,i,Y) - pbc.y;
    tmp_d.z = ORT(p,i,Z) - pbc.z;

    p_typ  = SORTE(p,i);
#ifdef TWOD
    jstart = (((p==q) && (pbc.x==0) && (pbc.y==0))               ? i+1 : 0);
#else
    jstart = (((p==q) && (pbc.x==0) && (pbc.y==0) && (pbc.z==0)) ? i+1 : 0);
#endif
    qptr   = &ORT(q,jstart,X);

    /* for each atom in neighbouring cell */
    for (j = jstart; j < q->n; ++j) {

      /* calculate distance */
      d.x = *qptr - tmp_d.x; ++qptr;
      d.y = *qptr - tmp_d.y; ++qptr;
      d.z = *qptr - tmp_d.z; ++qptr;

      q_typ = SORTE(q,j);
      col1  = p_typ * ntypes + q_typ;
      col2  = q_typ * ntypes + p_typ;
      r2    = SPROD(d,d);

#ifdef DEBUG
      if (0==r2) { char msgbuf[256];
        sprintf(msgbuf, "Distance is zero between particles %d and %d!\n",
                NUMMER(p,i), NUMMER(q,j));
        error(msgbuf);
      }
#endif

      /* compute pair interactions, first on particle i */
      if (r2 <= pair_pot.end[col1]) {
        PAIR_INT(pot_zwi, pot_grad, pair_pot, col1, inc, r2, is_short)

        /* store force in temporary variable */
        force.x = d.x * pot_grad;
        force.y = d.y * pot_grad;
        force.z = d.z * pot_grad;

        /* accumulate forces */
        pfptr = &KRAFT(p,i,X);
        *pfptr     += force.x; 
        *(++pfptr) += force.y; 
        *(++pfptr) += force.z; 

        /* the first half of the pot. energy of this bond */
        pot_zwi     *= 0.5;
        *Epot       += pot_zwi;
        POTENG(p,i) += pot_zwi;

        /* for the virial, we take the mean forces on the two particles */
        force.x       *= 0.5;
        force.y       *= 0.5;
        force.z       *= 0.5;
        pot_grad      *= 0.5;

#ifdef P_AXIAL
        tmp_vir_vect.x -= d.x * force.x;
        tmp_vir_vect.y -= d.y * force.y;
        tmp_vir_vect.z -= d.z * force.z;
#else
        tmp_virial     -= r2 * pot_grad;
#endif

#ifdef STRESS_TENS
        if (do_press_calc) {
          PRESSTENS(p,i,xx) -= d.x * force.x;
          PRESSTENS(p,i,yy) -= d.y * force.y;
          PRESSTENS(p,i,zz) -= d.z * force.z;
          PRESSTENS(p,i,yz) -= d.y * force.z;
          PRESSTENS(p,i,zx) -= d.z * force.x;
          PRESSTENS(p,i,xy) -= d.x * force.y;
	}
#endif
      }

      /* compute pair interactions, now on particle j */
      if (r2 <= pair_pot.end[col2]) {
        if (col1!=col2) {
          PAIR_INT(pot_zwi, pot_grad, pair_pot, col2, inc, r2, is_short);
	}

        /* store force in temporary variable */
        force.x = d.x * pot_grad;
        force.y = d.y * pot_grad;
        force.z = d.z * pot_grad;

        /* accumulate forces */
        qfptr = &KRAFT(q,j,X);
        *qfptr     -= force.x; 
        *(++qfptr) -= force.y; 
        *(++qfptr) -= force.z; 

        /* the second half of the pot. energy of this bond */
        pot_zwi     *= 0.5;
        *Epot       += pot_zwi;
        POTENG(q,j) += pot_zwi;

        /* for the virial, we take the mean forces on the two particles */
        force.x       *= 0.5;
        force.y       *= 0.5;
        force.z       *= 0.5;
        pot_grad      *= 0.5;

#ifdef P_AXIAL
        tmp_vir_vect.x -= d.x * force.x;
        tmp_vir_vect.y -= d.y * force.y;
        tmp_vir_vect.z -= d.z * force.z;
#else
        tmp_virial     -= r2 * pot_grad;
#endif

#ifdef STRESS_TENS
        if (do_press_calc) {
          PRESSTENS(q,j,xx) -= d.x * force.x;
          PRESSTENS(q,j,yy) -= d.y * force.y;
          PRESSTENS(q,j,zz) -= d.z * force.z;
          PRESSTENS(q,j,yz) -= d.y * force.z;
          PRESSTENS(q,j,zx) -= d.z * force.x;
          PRESSTENS(q,j,xy) -= d.x * force.y;
	}
#endif
      }

      /* compute host electron density */
      if (r2 < rho_h_tab.end[col1])  {
        VAL_FUNC(rho_h, rho_h_tab, col1, inc, r2, is_short);
        EAM_RHO(p,i) += rho_h; 
#ifdef EEAM
        EAM_P(p,i) += rho_h*rho_h;
#endif
      }
      if (r2 < rho_h_tab.end[col2]) {
        if (col1!=col2) {
          VAL_FUNC(rho_h, rho_h_tab, col2, inc, r2, is_short);
        }
        EAM_RHO(q,j) += rho_h; 
#ifdef EEAM
        EAM_P(q,j) += rho_h*rho_h; 
#endif
      }

    } /* for j */
  } /* for i */

  if (is_short==1) printf("\n Short distance!\n");

#ifdef P_AXIAL
  *Vir_xx += tmp_vir_vect.x;
  *Vir_yy += tmp_vir_vect.y;
  *Virial += tmp_vir_vect.x;
  *Virial += tmp_vir_vect.y;
#ifndef TWOD
  *Vir_zz += tmp_vir_vect.z;
  *Virial += tmp_vir_vect.z;
#endif
#else
  *Virial += tmp_virial;
#endif 

}

#endif /* ASYMPOT */

/******************************************************************************
*
*  compute embedding energy and its derivative for all atoms
*
******************************************************************************/

void do_embedding_energy(void)
{
  int k;

#ifdef _OPENMP
#pragma omp parallel for schedule(runtime) reduction(+:tot_pot_energy)
#endif
  for (k=0; k<NCELLS; k++) {
    int  i, idummy=0;
    real pot;
    cell *p;
    p = CELLPTR(k);
    for (i=0; i<p->n; i++) {
      PAIR_INT( pot, EAM_DF(p,i), embed_pot, SORTE(p,i), 
                ntypes, EAM_RHO(p,i), idummy);
      POTENG(p,i)    += pot;
      tot_pot_energy += pot;
#ifdef EEAM
      PAIR_INT( pot, EAM_DM(p,i), emod_pot, SORTE(p,i), 
                ntypes, EAM_P(p,i), idummy);
      POTENG(p,i)    += pot;
      tot_pot_energy += pot;
#endif
    }
  }
}


/******************************************************************************
*
*  second force loop, calculates the force and the energy 
*  caused by the embedding electron density
*  uses Phi(r2), Rho(r2), F(rho) and its derivatives
*  also used for EEAM
*
******************************************************************************/

void do_forces_eam2(cell *p, cell *q, vektor pbc, real *Virial, 
                    real *Vir_xx, real *Vir_yy, real *Vir_zz,
                    real *Vir_yz, real *Vir_zx, real *Vir_xy)
{
  int i,j,k,same_cell;
  vektor d, tmp_d, force;
  real r2;
  int  is_short=0, idummy=0;
  int  jstart, q_typ, p_typ;
  int  col1, col2, inc=ntypes*ntypes;
  real *qptr, *pfptr, *qfptr, *qpdptr, *ppdptr, *qpoptr, *ppoptr;
  real tmp_virial=0.0;
#ifdef P_AXIAL
  vektor tmp_vir_vect = {0.0, 0.0, 0.0};
#endif
  real eam2_force, rho_i_strich, rho_j_strich;
#ifdef EEAM
  real rho_i, rho_j;
#endif

  /* for each atom in first cell */
  for (i=0; i<p->n; ++i) {

    tmp_d.x = ORT(p,i,X) - pbc.x;
    tmp_d.y = ORT(p,i,Y) - pbc.y;
    tmp_d.z = ORT(p,i,Z) - pbc.z;

    p_typ   = SORTE(p,i);

#ifdef TWOD
    same_cell = ((p==q) && (pbc.x==0) && (pbc.y==0));
#else
    same_cell = ((p==q) && (pbc.x==0) && (pbc.y==0) && (pbc.z==0));
#endif

    jstart = (same_cell ? i+1 : 0);
    qptr   = &ORT(q,jstart,X);

    /* for each atom in neighbouring cell */
    for (j=jstart; j<q->n; ++j) {

      /* calculate distance */ 
      d.x = *qptr - tmp_d.x; ++qptr;
      d.y = *qptr - tmp_d.y; ++qptr;
      d.z = *qptr - tmp_d.z; ++qptr;

      q_typ = SORTE(q,j);
      r2    = SPROD(d,d);
      col1  = q_typ * ntypes + p_typ;
      col2  = p_typ * ntypes + q_typ;

      if ((r2 < rho_h_tab.end[col1]) || (r2 < rho_h_tab.end[col2])) {

        /* take care: particle i gets its rho from particle j.
           This is tabulated in column p_typ*ntypes+q_typ.
           Here we need the giving part from column q_typ*ntypes+p_typ.
         */

        /* rho_strich_i(r_ij) */
#ifndef EEAM
        DERIV_FUNC(rho_i_strich, rho_h_tab, col1, inc, r2, is_short);
#else
        /* rho_strich_i(r_ij) and rho_i(r_ij) */
        PAIR_INT(rho_i, rho_i_strich, rho_h_tab, col1, inc, r2, is_short);
#endif

        /* rho_strich_j(r_ij) */
        if (col1==col2) {
          rho_j_strich = rho_i_strich;
#ifdef EEAM
          rho_j = rho_i;
#endif
	} else {
#ifndef EEAM
          DERIV_FUNC(rho_j_strich, rho_h_tab, col2, inc, r2, is_short);
#else
          PAIR_INT(rho_j, rho_j_strich, rho_h_tab, col2, inc, r2, is_short);
#endif
	}

        /* put together (dF_i and dF_j are by 0.5 too big) */
        eam2_force = 0.5 * (EAM_DF(p,i)*rho_j_strich+EAM_DF(q,j)*rho_i_strich);
#ifdef EEAM
        /* 0.5 times 2 from derivative simplified to 1 */
        eam2_force += (EAM_DM(p,i) * rho_j * rho_j_strich +
                     + EAM_DM(q,j) * rho_i * rho_i_strich);
#endif

        /* store force in temporary variable */
        force.x = d.x * eam2_force;
        force.y = d.y * eam2_force;
        force.z = d.z * eam2_force;

        /* accumulate forces */
        pfptr = &KRAFT(p,i,X);
        qfptr = &KRAFT(q,j,X);
        *pfptr     += force.x; 
        *qfptr     -= force.x; 
        *(++pfptr) += force.y; 
        *(++qfptr) -= force.y; 
        *(++pfptr) += force.z; 
        *(++qfptr) -= force.z; 

#ifdef P_AXIAL
        tmp_vir_vect.x -= d.x * force.x;
        tmp_vir_vect.y -= d.y * force.y;
        tmp_vir_vect.z -= d.z * force.z;
#else
        tmp_virial     -= r2  * eam2_force;
#endif

#ifdef STRESS_TENS
        if (do_press_calc) {
          /* avoid double counting of the virial */
          force.x *= 0.5;
          force.y *= 0.5;
          force.z *= 0.5;
 
          PRESSTENS(p,i,xx) -= d.x * force.x;
          PRESSTENS(p,i,yy) -= d.y * force.y;
          PRESSTENS(p,i,zz) -= d.z * force.z;
          PRESSTENS(p,i,yz) -= d.y * force.z;
          PRESSTENS(p,i,zx) -= d.z * force.x;
          PRESSTENS(p,i,xy) -= d.x * force.y;

          PRESSTENS(q,j,xx) -= d.x * force.x;
          PRESSTENS(q,j,yy) -= d.y * force.y;
          PRESSTENS(q,j,zz) -= d.z * force.z;
          PRESSTENS(q,j,yz) -= d.y * force.z;
          PRESSTENS(q,j,zx) -= d.z * force.x;
          PRESSTENS(q,j,xy) -= d.x * force.y;
	}
#endif

      } /* if in the cutoff range */
    } /* for j */
  } /* for i */

  /* print warning if short distance occurred */
  if (is_short==1) fprintf(stderr, "\n Short distance!\n");

#ifdef P_AXIAL
  *Vir_xx += tmp_vir_vect.x;
  *Vir_yy += tmp_vir_vect.y;
  *Virial += tmp_vir_vect.x;
  *Virial += tmp_vir_vect.y;
#ifndef TWOD
  *Vir_zz += tmp_vir_vect.z;
  *Virial += tmp_vir_vect.z;
#endif
#else
  *Virial += tmp_virial;
#endif 

} /* do_forces_eam2 */
