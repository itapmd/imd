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

#if defined(PVPCRAY) || defined(SX4)
#include "imd_forces_vec.c"
#else

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
#ifndef MONOLJ
    p_typ   = p->sorte[i];
#endif

    jstart = (p==q ? i+1 : 0);
    qptr   = q->ort + DIM * jstart;
    
    for (j = jstart; j < q->n; ++j) {
      
      /* Calculate distance  */
      d.x = *qptr - tmp_d.x;
      ++qptr;
      d.y = *qptr - tmp_d.y;
      ++qptr;
#ifndef TWOD
      d.z = *qptr - tmp_d.z;
      ++qptr;
#endif

      radius2 = SPROD(d,d);

#ifndef NODBG_DIST
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
#else
      if (0==radius2) error("Distance is zero.");
#endif

      if (radius2 <= r2_cut) {

        /* Calculate force, if distance smaller than cutoff */ 
#ifdef MONOLJ
        sig_d_rad2  = 2.0 / radius2;
        sig_d_rad6  = sig_d_rad2 * sig_d_rad2 * sig_d_rad2;
        sig_d_rad12 = sig_d_rad6 * sig_d_rad6;

        pot_grad = - 6 * sig_d_rad2 * ( sig_d_rad12 - sig_d_rad6 );
        pot_zwi  = sig_d_rad12 - 2.0 * sig_d_rad6;
#else
      	/* Check for distances, shorter than minimal distance in pot. table */
	if (radius2 <= r2_0) {
	  radius2 = r2_0; 
	}

	/* Indices into potential table */
	k     = (int) ((radius2 - r2_0) * inv_r2_step);
	q_typ = q->sorte[j];
	chi = (radius2 - r2_0 - k * r2_step) * inv_r2_step;
	
	/* A single access to the potential table involves two multiplications 
	   We use a intermediate pointer to aviod this as much as possible.
	   Note: This relies on layout of the pot-table in memory!!! */

	potptr = PTR_3D_V(potential, k, p_typ, q_typ , pot_dim);
	pot_k0 = *potptr; potptr += pot_dim.y * pot_dim.z;
	pot_k1 = *potptr; potptr += pot_dim.y * pot_dim.z;
	pot_k2 = *potptr;

	dv  = pot_k1 - pot_k0;
	d2v = pot_k2 - 2 * pot_k1 + pot_k0;

	/* Norm of Gradient */
	pot_grad = 2 * inv_r2_step * ( dv + (chi - 0.5) * d2v );
        /* Potential energy of atom */
	pot_zwi =  pot_k0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;
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

#ifdef ORDPAR
        if (radius2 < op_r2_cut[p_typ][q_typ]) {
           p->pot_eng[i] += pot_zwi * op_weight[p_typ][q_typ];         
           q->pot_eng[j] += pot_zwi * op_weight[q_typ][p_typ];         
        }
#else
#ifndef MONOLJ
	q->pot_eng[j] += pot_zwi;
	p->pot_eng[i] += pot_zwi;
#endif
#endif
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
#ifdef TRANSPORT
        p->heatcond[i] += pot_zwi - radius2 * pot_grad;
        q->heatcond[j] += pot_zwi - radius2 * pot_grad;
#endif

#if (defined(TTBP) && !defined(TWOD))
	/* 2. Cutoff: make neighbor tables for TTBP */
        if (radius2 <= ttbp_r2_cut[p_typ][q_typ]) {

          neightab *neigh;
          real  *tmp_ptr;

          /* update neighbor table of particle i */
          neigh = p->neigh[i];
          if (neigh->n_max == neigh->n) {
            error("neighbor table too small, increase ttbp_len");
          }
          neigh->typ[neigh->n] = q_typ;
          tmp_ptr  = &neigh->dist[3*neigh->n];
	  *tmp_ptr = d.x; ++tmp_ptr; 
	  *tmp_ptr = d.y; ++tmp_ptr; 
	  *tmp_ptr = d.z;
          neigh->n++;

          /* update neighbor table of particle j */
          neigh = q->neigh[j];
          if (neigh->n_max == neigh->n) {
            error("neighbor table too small, increase ttbp_len");
          }
          neigh->typ[neigh->n] = p_typ;
          tmp_ptr  = &neigh->dist[3*neigh->n];
	  *tmp_ptr = -d.x; ++tmp_ptr; 
	  *tmp_ptr = -d.y; ++tmp_ptr; 
	  *tmp_ptr = -d.z;
          neigh->n++;

        }
#endif  /* TTPB */

      } /* if */

    } /* for j */

  } /* for i */

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
} /* do_forces */

#endif /* do_forces for scalar processors */
