
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

void do_forces(cell *p, cell *q, vektor pbc, real *Epot, 
               real *Virial, real *Vir_x, real *Vir_y, real *Vir_z)
{
  int i,j,k;
  vektor d;
  vektor tmp_d;
  vektor force;
  real radius2;
  real tmp_virial;
#ifdef P_AXIAL
  vektor tmp_vir_vect;
#endif
  real pot_zwi;
  real pot_grad;
  int jstart,jend;
  int col, is_short=0, inc = ntypes * ntypes;
  int q_typ,p_typ;
  real *qptr, *pfptr, *qfptr, *qpdptr, *ppdptr, *qpoptr, *ppoptr;
  
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
      d.x = *qptr - tmp_d.x; ++qptr;
      d.y = *qptr - tmp_d.y; ++qptr;
#ifndef TWOD
      d.z = *qptr - tmp_d.z; ++qptr;
#endif

      col  = p_typ * ntypes + q_typ;
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
        PAIR_INT_MONOLJ(pot_zwi,pot_grad,radius2)
#else
      /* 1. Cutoff: pair potential */
      if (radius2 <= pair_pot.end[col]) {
        PAIR_INT2(pot_zwi,pot_grad,pair_pot,col,inc,radius2,is_short)
#endif  /* MONOLJ */
        
        /* Store forces in temp */
        force.x = d.x * pot_grad;
        force.y = d.y * pot_grad;
#ifndef TWOD
        force.z = d.z * pot_grad;
#endif

        /* Accumulate forces */
        pfptr = p->kraft + DIM * i;
        qfptr = q->kraft + DIM * j;
        *pfptr     += force.x; 
        *qfptr     -= force.x; 
        *(++pfptr) += force.y; 
        *(++qfptr) -= force.y; 
#ifndef TWOD
        *(++pfptr) += force.z; 
        *(++qfptr) -= force.z; 
#endif

#ifndef MONOLJ
#ifdef ORDPAR
        if (radius2 < op_r2_cut[p_typ][q_typ]) {
	  p->pot_eng[i] += op_weight[p_typ][q_typ]*pot_zwi;
	  q->pot_eng[j] += op_weight[q_typ][p_typ]*pot_zwi;
	  p->nbanz[i]++;
	  q->nbanz[j]++;
        }
#else
        q->pot_eng[j] += pot_zwi;
        p->pot_eng[i] += pot_zwi;
#endif
#endif
        *Epot   += pot_zwi;

#ifdef P_AXIAL
        tmp_vir_vect.x -= d.x * force.x;
        tmp_vir_vect.y -= d.y * force.y;
#ifndef TWOD
        tmp_vir_vect.z -= d.z * force.z;
#endif
#else
        tmp_virial     -= radius2 * pot_grad;  
#endif

#ifdef STRESS_TENS
        ppdptr = p->presstens + DIM * i;
        qpdptr = q->presstens + DIM * j;
        *ppdptr     -= d.x * force.x;
        *qpdptr     -= d.x * force.x;
        *(++ppdptr) -= d.y * force.y;
        *(++qpdptr) -= d.y * force.y;
#ifdef TWOD
#ifndef SHOCK
        p->presstens_offdia[i] -= d.x * force.y;
        q->presstens_offdia[j] -= d.x * force.y;
#endif
#else
        *(++ppdptr) -= d.z * force.z;
        *(++qpdptr) -= d.z * force.z;
#ifndef SHOCK
        ppoptr = p->presstens_offdia + DIM * i;
        qpoptr = q->presstens_offdia + DIM * j;
        *ppoptr     -= d.y * force.z;
        *qpoptr     -= d.y * force.z;
        *(++ppoptr) -= d.z * force.x;
        *(++qpoptr) -= d.z * force.x;
        *(++ppoptr) -= d.x * force.y;
        *(++qpoptr) -= d.x * force.y;
#endif
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
      if (radius2 <= smooth_pot.end[col]) {
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
  if (is_short==1) printf("\n Short distance!\n");
#endif

#ifdef P_AXIAL
  *Vir_x  += tmp_vir_vect.x;
  *Vir_y  += tmp_vir_vect.y;
  *Virial += tmp_vir_vect.x;
  *Virial += tmp_vir_vect.y;
#ifndef TWOD
  *Vir_z  += tmp_vir_vect.z;
  *Virial += tmp_vir_vect.z;
#endif
#else
  *Virial += tmp_virial;
#endif 
#endif /* TERSOFF */ 

}
