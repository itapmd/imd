
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
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

    tmp_d.x = p->ort X(i) - pbc.x;
    tmp_d.y = p->ort Y(i) - pbc.y;
    tmp_d.z = p->ort Z(i) - pbc.z;

    p_typ  = SORTE(p,i);
    jstart = (p==q ? i+1 : 0);
    qptr   = q->ort + DIM * jstart;

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
        pfptr = p->kraft + DIM * i;
        *pfptr     += force.x; 
        *(++pfptr) += force.y; 
        *(++pfptr) += force.z; 

        /* the first half of the pot. energy of this bond */
        pot_zwi       *= 0.5;
        *Epot         += pot_zwi;
        p->pot_eng[i] += pot_zwi;

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
          p->presstens[i].xx -= d.x * force.x;
          p->presstens[i].yy -= d.y * force.y;
          p->presstens[i].zz -= d.z * force.z;
          p->presstens[i].yz -= d.y * force.z;
          p->presstens[i].zx -= d.z * force.x;
          p->presstens[i].xy -= d.x * force.y;
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
        qfptr = q->kraft + DIM * j;
        *qfptr     -= force.x; 
        *(++qfptr) -= force.y; 
        *(++qfptr) -= force.z; 

        /* the second half of the pot. energy of this bond */
        pot_zwi       *= 0.5;
        *Epot         += pot_zwi;
        q->pot_eng[j] += pot_zwi;

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
          q->presstens[j].xx -= d.x * force.x;
          q->presstens[j].yy -= d.y * force.y;
          q->presstens[j].zz -= d.z * force.z;
          q->presstens[j].yz -= d.y * force.z;
          q->presstens[j].zx -= d.z * force.x;
          q->presstens[j].xy -= d.x * force.y;
	}
#endif
      }

      /* compute host electron density */
      if (r2 < rho_h_tab.end[col1])  {
        VAL_FUNC(rho_h, rho_h_tab, col1, inc, r2, is_short);
        p->eam2_rho_h[i] += rho_h; 
      }
      if (r2 < rho_h_tab.end[col2]) {
        if (col1!=col2) {
          VAL_FUNC(rho_h, rho_h_tab, col2, inc, r2, is_short);
        }
        q->eam2_rho_h[j] += rho_h; 
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
*  second force loop, calculates the force and the energy 
*  caused by the embedding electron density
*  uses Phi(r2), Rho(r2), F(rho) and its derivatives
*
******************************************************************************/

void do_forces_eam2(cell *p, cell *q, vektor pbc, real *Epot, real *Virial, 
                    real *Vir_xx, real *Vir_yy, real *Vir_zz,
                    real *Vir_yz, real *Vir_zx, real *Vir_xy)
{
  int i,j,k;
  vektor d;
  vektor tmp_d;
  vektor force;
  real r2;
  int  is_short=0, idummy=0;
  int  jstart, q_typ, p_typ;
  int  col1, col2, inc=ntypes*ntypes;
  real *qptr, *pfptr, *qfptr, *qpdptr, *ppdptr, *qpoptr, *ppoptr;
  real tmp_virial=0.0;
#ifdef P_AXIAL
  vektor tmp_vir_vect = {0.0, 0.0, 0.0};
#endif
  real eam2_energy, eam2_force;
  real f_i_strich, f_j_strich;
  real rho_i_strich, rho_j_strich;

  /* for each atom in first cell */
  for (i=0; i<p->n; ++i) {

    tmp_d.x = p->ort X(i) - pbc.x;
    tmp_d.y = p->ort Y(i) - pbc.y;
    tmp_d.z = p->ort Z(i) - pbc.z;

    p_typ   = SORTE(p,i);

    if (p==q) {
      /* f_i and f_i_strich */
      PAIR_INT(eam2_energy, f_i_strich, embed_pot, p_typ, 
               ntypes, p->eam2_rho_h[i], idummy);
      /* add energy only once per particle (p==q) */
      p->pot_eng[i]  += eam2_energy;
      *Epot          += eam2_energy;
    } else {      
      /* only f_i_strich */
      DERIV_FUNC(f_i_strich, embed_pot, p_typ, 
                 ntypes, p->eam2_rho_h[i], idummy);
    }

    jstart = (p==q ? i+1 : 0);
    qptr   = q->ort + DIM * jstart;

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

        /* f_j_strich(rho_h_j) */
        DERIV_FUNC(f_j_strich, embed_pot, q_typ, 
                   ntypes, q->eam2_rho_h[j], idummy);

        /* take care: particle i gets its rho from particle j.
           This is tabulated in column p_typ*ntypes+q_typ.
           Here we need the giving part from column q_typ*ntypes+p_typ.
         */

        /* rho_strich_i(r_ij) */
        DERIV_FUNC(rho_i_strich, rho_h_tab, col1, inc, r2, is_short);

        /* rho_strich_j(r_ij) */
        if (col1==col2) {
          rho_j_strich = rho_i_strich;
	} else {
          DERIV_FUNC(rho_j_strich, rho_h_tab, col2, inc, r2, is_short);
	}

        /* put together (f_i_strich and f_j_strich are by 0.5 too big) */
        eam2_force = 0.5 * (f_i_strich*rho_j_strich+f_j_strich*rho_i_strich);

        /* store force in temporary variable */
        force.x = d.x * eam2_force;
        force.y = d.y * eam2_force;
        force.z = d.z * eam2_force;

        /* accumulate forces */
        pfptr = p->kraft + DIM * i;
        qfptr = q->kraft + DIM * j;
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
 
          p->presstens[i].xx -= d.x * force.x;
          p->presstens[i].yy -= d.y * force.y;
          p->presstens[i].zz -= d.z * force.z;
          p->presstens[i].yz -= d.y * force.z;
          p->presstens[i].zx -= d.z * force.x;
          p->presstens[i].xy -= d.x * force.y;

          q->presstens[j].xx -= d.x * force.x;
          q->presstens[j].yy -= d.y * force.y;
          q->presstens[j].zz -= d.z * force.z;
          q->presstens[j].yz -= d.y * force.z;
          q->presstens[j].zx -= d.z * force.x;
          q->presstens[j].xy -= d.x * force.y;
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
