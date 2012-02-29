
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2011 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/*****************************************************************************
*
* imd_forces -- force loop
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"
#include "potaccess.h"

/******************************************************************************
*
*  do_forces, version for scalar processors
*
*  computes the forces between atoms in two given cells
*
******************************************************************************/

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
  int col, col2, is_short=0, inc = ntypes * ntypes;
  int jstart, q_typ, p_typ;
  real *qptr, *pfptr, *qfptr, *qpdptr, *ppdptr, *qpoptr, *ppoptr;
  
  tmp_virial     = 0.0;
#ifdef P_AXIAL
  tmp_vir_vect.x = 0.0;
  tmp_vir_vect.y = 0.0;
#ifndef TWOD
  tmp_vir_vect.z = 0.0;
#endif
#endif
    
  /* for each atom in first cell */
  for (i=0; i<p->n; ++i) {

    tmp_d.x = ORT(p,i,X) - pbc.x;
    tmp_d.y = ORT(p,i,Y) - pbc.y;
#ifndef TWOD
    tmp_d.z = ORT(p,i,Z) - pbc.z;
#endif

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
#ifndef TWOD
      d.z = *qptr - tmp_d.z; ++qptr;
#endif

      q_typ = SORTE(q,j);
      col   = p_typ * ntypes + q_typ;
      col2  = q_typ * ntypes + p_typ;
      r2    = SPROD(d,d);

#ifdef DEBUG
      if (0==r2) { char msgbuf[256];
        sprintf(msgbuf, "Distance is zero between particles %d and %d!\n",
                NUMMER(p,i), NUMMER(q,j));
        error(msgbuf);
      }
#endif

      /* compute pair interactions */
#if defined(PAIR) || defined(KEATING)
      /* PAIR and KEATING are mutually exclusive */
#if defined(PAIR)
      if (r2 <= pair_pot.end[col]) {
#ifdef LINPOT
        PAIR_INT_LIN(pot_zwi, pot_grad, pair_pot_lin, col, inc, r2, is_short)
#else
        PAIR_INT(pot_zwi, pot_grad, pair_pot, col, inc, r2, is_short)
#endif
#elif defined(KEATING)
      if (r2 < keat_r2_cut[p_typ][q_typ]) {
	PAIR_INT_KEATING(pot_zwi, pot_grad, p_typ, q_typ, r2)
#endif

        /* store force in temporary variable */
        force.x = d.x * pot_grad;
        force.y = d.y * pot_grad;
#ifndef TWOD
        force.z = d.z * pot_grad;
#endif

        /* accumulate forces */
        pfptr = &KRAFT(p,i,X);
        qfptr = &KRAFT(q,j,X);
        *pfptr     += force.x; 
        *qfptr     -= force.x; 
        *(++pfptr) += force.y; 
        *(++qfptr) -= force.y; 
#ifndef TWOD
        *(++pfptr) += force.z; 
        *(++qfptr) -= force.z; 
#endif
        *Epot      += pot_zwi;

#ifndef MONOLJ
        pot_zwi *= 0.5;   /* avoid double counting */
#ifdef NNBR
        if (r2 < nb_r2_cut[col ]) NBANZ(p,i)++;
        if (r2 < nb_r2_cut[col2]) NBANZ(q,j)++;
#endif
#ifdef ORDPAR
        if (r2 < op_r2_cut[col ]) POTENG(p,i) += op_weight[col ] * pot_zwi;
        if (r2 < op_r2_cut[col2]) POTENG(q,j) += op_weight[col2] * pot_zwi;
#else
        POTENG(p,i) += pot_zwi;
        POTENG(q,j) += pot_zwi;
#endif
#endif

#ifdef P_AXIAL
        tmp_vir_vect.x -= d.x * force.x;
        tmp_vir_vect.y -= d.y * force.y;
#ifndef TWOD
        tmp_vir_vect.z -= d.z * force.z;
#endif
#else
        tmp_virial     -= r2 * pot_grad;
#endif

#ifdef STRESS_TENS
        if (do_press_calc) {
          /* avoid double counting of the virial */
          force.x *= 0.5;
          force.y *= 0.5;
#ifndef TWOD
          force.z *= 0.5;
#endif
          PRESSTENS(p,i,xx) -= d.x * force.x;
          PRESSTENS(q,j,xx) -= d.x * force.x;
          PRESSTENS(p,i,yy) -= d.y * force.y;
          PRESSTENS(q,j,yy) -= d.y * force.y;
          PRESSTENS(p,i,xy) -= d.x * force.y;
          PRESSTENS(q,j,xy) -= d.x * force.y;
#ifndef TWOD
          PRESSTENS(p,i,zz) -= d.z * force.z;
          PRESSTENS(q,j,zz) -= d.z * force.z;
          PRESSTENS(p,i,yz) -= d.y * force.z;
          PRESSTENS(q,j,yz) -= d.y * force.z;
          PRESSTENS(p,i,zx) -= d.z * force.x;
          PRESSTENS(q,j,zx) -= d.z * force.x;
#endif
	}
#endif
      }
#endif /* PAIR || KEATING */

#ifdef EAM2
      /* compute host electron density */
      if (r2 < rho_h_tab.end[col])  {
        VAL_FUNC(rho_h, rho_h_tab, col,  inc, r2, is_short);
        EAM_RHO(p,i) += rho_h; 
#ifdef EEAM
        EAM_P(p,i) += rho_h*rho_h; 
#endif
      }
      if (p_typ==q_typ) {
        if (r2 < rho_h_tab.end[col]) 
          {EAM_RHO(q,j) += rho_h;
#ifdef EEAM
           EAM_P(q,j) += rho_h*rho_h;
#endif
          }
      } else {
        col2 = q_typ * ntypes + p_typ;
        if (r2 < rho_h_tab.end[col2]) {
          VAL_FUNC(rho_h, rho_h_tab, col2, inc, r2, is_short);
          EAM_RHO(q,j) += rho_h; 
#ifdef EEAM
          EAM_P(q,j) += rho_h*rho_h; 
#endif
        }
      }
#endif

#ifdef COVALENT
      /* make neighbor tables for covalent systems */
      if (r2 <= neightab_r2cut[col]) {

        neightab *neigh;
        real  *tmp_ptr;

        /* update neighbor table of particle i */
        neigh = NEIGH(p,i);
        if (neigh->n_max <= neigh->n) {
          increase_neightab( neigh, neigh->n_max + NEIGH_LEN_INC );
        }
        neigh->typ[neigh->n] = q_typ;
        neigh->cl [neigh->n] = q;
        neigh->num[neigh->n] = j;
        tmp_ptr  = &neigh->dist[3*neigh->n];
        *tmp_ptr = d.x; ++tmp_ptr; 
        *tmp_ptr = d.y; ++tmp_ptr; 
        *tmp_ptr = d.z;
        neigh->n++;

        /* update neighbor table of particle j */
        neigh = NEIGH(q,j);
        if (neigh->n_max <= neigh->n) {
          increase_neightab( neigh, neigh->n_max + NEIGH_LEN_INC );
        }
        neigh->typ[neigh->n] = p_typ;
        neigh->cl [neigh->n] = p;
        neigh->num[neigh->n] = i;
        tmp_ptr  = &neigh->dist[3*neigh->n];
        *tmp_ptr = -d.x; ++tmp_ptr; 
        *tmp_ptr = -d.y; ++tmp_ptr; 
        *tmp_ptr = -d.z;
        neigh->n++;
      }
#endif  /* COVALENT */

    } /* for j */
  } /* for i */

#ifdef DEBUG
  if (is_short==1) printf("\n Short distance!\n");
#endif
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
