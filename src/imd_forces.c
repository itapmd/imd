
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

/******************************************************************************
*
*  do_forces, version for scalar processors
*
*  computes the forces between atoms in two given cells
*
******************************************************************************/

void do_forces(cell *p, cell *q, vektor pbc)

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
  
  tmp_virial     = 0.0;
#ifdef P_AXIAL
  tmp_vir_vect.x = 0.0;
  tmp_vir_vect.y = 0.0;
#ifndef TWOD
  tmp_vir_vect.z = 0.0;
#endif
#endif
    
  /* For each atom in first cell */
  for (i=0; i<p->n; ++i) {

    /* For each atom in neighbouring cell */
    /* Some compilers don't find the expressions that are invariant 
       to the inner loop. I'll have to define my own temp variables. */

    tmp_d.x = p->ort X(i) - pbc.x;
    tmp_d.y = p->ort Y(i) - pbc.y;
#ifndef TWOD
    tmp_d.z = p->ort Z(i) - pbc.z;
#endif
    p_typ   = SORTE(p,i);

    jstart = (p==q ? i+1 : 0);
    qptr   = q->ort + DIM * jstart;
    
    for (j = jstart; j < q->n; ++j) {

      q_typ = SORTE(q,j);
      
      /* Calculate distance  */
      d.x = *qptr - tmp_d.x;
      ++qptr;
      d.y = *qptr - tmp_d.y;
      ++qptr;
#ifndef TWOD
      d.z = *qptr - tmp_d.z;
      ++qptr;
#endif

      column  = p_typ * ntypes + q_typ;
      radius2 = SPROD(d,d);

      if (0==radius2) { char msgbuf[256];
#ifdef TWOD
        sprintf(msgbuf,
                "Distance is zero: nrs=%d %d\norte: %f %f, %f %f\n",
                NUMMER(p,i),NUMMER(q,i),
                p->ort X(i),p->ort Y(i),
                q->ort X(j),q->ort Y(j));
        error(msgbuf);
#else
        sprintf(msgbuf,
                "Distance is zero: nrs=%d %d\norte: %f %f %f, %f %f %f\n",
                NUMMER(p,i),NUMMER(q,i),
                p->ort X(i),p->ort Y(i),p->ort Z(i),
                q->ort X(j),q->ort Y(j),q->ort Z(j));
        error(msgbuf);
#endif
      }

#ifndef TERSOFF
#ifdef MONOLJ
      if (radius2 <= monolj_r2_cut) {
        pair_int_monolj(&pot_zwi,&pot_grad,radius2);
#else
      /* 1. Cutoff: pair potential */
      if (radius2 <= pair_pot.end[column]) {

      	/* Check for distances, shorter than minimal distance in pot. table */
        if (1==pair_int2(&pot_zwi,&pot_grad,&pair_pot,column,radius2)) {
          r2_short = MIN(r2_short,radius2);
	  is_short=1; 
	}
#endif  /* MONOLJ */
        
        /* Store forces in temp */
        force.x = d.x * pot_grad;
        force.y = d.y * pot_grad;
#ifndef TWOD
        force.z = d.z * pot_grad;
#endif

        /* Accumulate forces */
        p->kraft X(i) += force.x;
        p->kraft Y(i) += force.y;
        q->kraft X(j) -= force.x;
        q->kraft Y(j) -= force.y;
#ifndef TWOD
        p->kraft Z(i) += force.z;
        q->kraft Z(j) -= force.z;
#endif

#ifndef MONOLJ
#ifdef ORDPAR
        if (radius2 < op_r2_cut[p_typ][q_typ]) {
           p->pot_eng[i] += pot_zwi * op_weight[p_typ][q_typ];         
           q->pot_eng[j] += pot_zwi * op_weight[q_typ][p_typ];         
#ifndef TWOD
	   p->nbanz[i]++;
	   q->nbanz[j]++;
#endif /* TWOD or NOT TWOD (Hamlet, Act 3, Scene 1) */
        }
#else
        q->pot_eng[j] += pot_zwi;
        p->pot_eng[i] += pot_zwi;
#endif
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
#ifdef TWOD
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
#ifdef NVX
        p->heatcond[i] += pot_zwi - radius2 * pot_grad;
        q->heatcond[j] += pot_zwi - radius2 * pot_grad;
#endif
      }  /* if */
#endif /* TERSOFF */

#ifdef TTBP
      /* 2. Cutoff: make neighbor tables for TTBP */
      if (radius2 <= smooth_pot.end[column]) {
#endif
#ifdef TERSOFF
      /* 2. Cutoff: make neighbor tables for TERSOFF */ 
      if (radius2 <= ter_r2_cut[p_typ][q_typ]) {
#endif 
#ifdef COVALENT
        neightab *neigh;
        real  *tmp_ptr;

        /* update neighbor table of particle i */
        neigh = p->neigh[i];
        if (neigh->n_max <= neigh->n) {
          error("neighbor table too small, increase neigh_len");
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
        neigh = q->neigh[j];
        if (neigh->n_max <= neigh->n) {
          error("neighbor table too small, increase neigh_len");
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

#ifndef TERSOFF
#ifndef MONOLJ
  if (is_short==1) printf("\n Short distance! r2: %f\n",r2_short);
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
#endif /* TERSOFF */ 

}
