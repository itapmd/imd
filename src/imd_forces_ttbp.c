/*****************************************************************************
*
*  imd_forces_ttbp.c -- force loop for three body forces
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/******************************************************************************
*
*  three body forces for TTBP, using neighbor tables computed in do_forces
*
******************************************************************************/

void do_forces_ttbp(cell *p)
{
  static real   *r2 = NULL, *r = NULL, *pot = NULL, *grad = NULL;
  static vektor *d  = NULL;
  neightab *neigh;
  vektor force_j, force_k;
  cell   *jcell, *kcell;
  int    i, j, k, p_typ, j_typ, knum, jnum;
  real   *tmpptr, *potptr;
  real   pot_zwi, tmp_pot, tmp_grad, tmp, tmp_j, tmp_k;
  real   chi, dv, d2v, pot_k0, pot_k1, pot_k2, r2j;
  real   tmp_f2, cos_theta, tmp_sp;
  real   tmp_virial = 0.0;
#ifdef P_AXIAL
  vektor tmp_vir_vect = {0.0, 0.0, 0.0};
#endif

  d    = (vektor *) realloc( d,    ttbp_len * sizeof(vektor) );
  r2   = (real *)   realloc( r2,   ttbp_len * sizeof(real)   );
  r    = (real *)   realloc( r,    ttbp_len * sizeof(real)   );
  pot  = (real *)   realloc( pot,  ttbp_len * sizeof(real)   );
  grad = (real *)   realloc( grad, ttbp_len * sizeof(real)   );
  if ((d==NULL) || (r2==NULL) || (r==NULL) || (pot==NULL) || (grad==NULL))
    error("cannot allocate memory for temporary neighbor data");

  /*           j
              /|
             / |
            i  |
             \ |
              \|
               k       */

  /* For each atom in cell */
  for (i=0; i<p->n; ++i) {

    p_typ   = p->sorte[i];
    neigh   = p->neigh[i];

    /* construct some data for all neighbors */
    tmpptr = neigh->dist;
    for (j=0; j<neigh->n; ++j) {

      /* type, distance vector, radii */
      j_typ   = neigh->typ[j];
      d[j].x  = *tmpptr++;
      d[j].y  = *tmpptr++;
      d[j].z  = *tmpptr++;
      r2[j]   = SPROD(d[j],d[j]);
      r[j]    = sqrt(r2[j]);

      /* smoothing potential */
      r2j     = MAX( r2[j], ttbp_r2_0);
      k       = (int) ((r2j - ttbp_r2_0) * ttbp_inv_r2_step);
      chi     = (r2j - ttbp_r2_0 - k * ttbp_r2_step) * ttbp_inv_r2_step;
      potptr  = PTR_3D_V(ttbp_potential, k, p_typ, j_typ, ttbp_pot_dim);
      pot_k0  = *potptr; potptr += ttbp_pot_dim.y * ttbp_pot_dim.z;
      pot_k1  = *potptr; potptr += ttbp_pot_dim.y * ttbp_pot_dim.z;
      pot_k2  = *potptr;
      dv      = pot_k1 - pot_k0;
      d2v     = pot_k2 - 2 * pot_k1 + pot_k0;
      grad[j] = 2 * ttbp_inv_r2_step * ( dv + (chi - 0.5) * d2v );
      pot[j]  = pot_k0 + chi * dv + 0.5 * chi * (chi-1) * d2v;

    }

    /* for each pair of neighbors */
    for (j=0; j<neigh->n-1; ++j)
      for (k=j+1; k<neigh->n; ++k) {

        /* FOURIER potential term */
        tmp_sp    = SPROD(d[j],d[k]);
        cos_theta = tmp_sp / (r[j] * r[k]);
        tmp       = cos_theta + 1.0 / ttbp_sp[p_typ];
        tmp_pot   = ttbp_constant[p_typ] * tmp * tmp;
        tmp_grad  = ttbp_constant[p_typ] * 2 * tmp;

        /* smoothing potential, total potential */
        tmp_f2          = pot[j] * pot[k];
        pot_zwi         = tmp_pot * tmp_f2;
        if (NUMMER(p,i)>=0) tot_pot_energy += pot_zwi;

        /* forces */
        tmp   = -tmp_f2 * tmp_grad / (r[j] * r[k]);
        tmp_j = tmp * tmp_sp / r2[j] + tmp_pot * grad[j] * pot[k];
        tmp_k = tmp * tmp_sp / r2[k] + tmp_pot * grad[k] * pot[j];

        force_j.x = tmp_j * d[j].x - tmp * d[k].x;
        force_j.y = tmp_j * d[j].y - tmp * d[k].y;
        force_j.z = tmp_j * d[j].z - tmp * d[k].z;

        force_k.x = tmp_k * d[k].x - tmp * d[j].x;
        force_k.y = tmp_k * d[k].y - tmp * d[j].y;
        force_k.z = tmp_k * d[k].z - tmp * d[j].z;

        /* update force on particle i */
        p->kraft X(i) += force_j.x + force_k.x;
        p->kraft Y(i) += force_j.y + force_k.y;
        p->kraft Z(i) += force_j.z + force_k.z;
        p->pot_eng[i] += pot_zwi;

        /* update force on particle j */
        jcell = neigh->cl [j];
        jnum  = neigh->num[j];
        jcell->kraft X(jnum) -= force_j.x;
        jcell->kraft Y(jnum) -= force_j.y;
        jcell->kraft Z(jnum) -= force_j.z;
        jcell->pot_eng[jnum] += pot_zwi;

        /* update force on particle k */
        kcell = neigh->cl [k];
        knum  = neigh->num[k];
        kcell->kraft X(knum) -= force_k.x;
        kcell->kraft Y(knum) -= force_k.y;
        kcell->kraft Z(knum) -= force_k.z;
        kcell->pot_eng[knum] += pot_zwi;

#ifdef P_AXIAL
        tmp_vir_vect.x += d[j].x * force_j.x + d[k] * force_k.x;
        tmp_vir_vect.y += d[j].y * force_j.y + d[k] * force_k.y;
        tmp_vir_vect.z += d[j].z * force_j.z + d[k] * force_k.z;
#else
        tmp_virial     += SPROD(d[j],force_j) + SPROD(d[k],force_k);
#endif

    } /* neighbor pairs */

  } /* i */

#ifdef P_AXIAL
  vir_vect.x += tmp_vir_vect.x;
  virial     += tmp_vir_vect.x;
  vir_vect.y += tmp_vir_vect.y;
  virial     += tmp_vir_vect.y;
  vir_vect.z += tmp_vir_vect.z;
  virial     += tmp_vir_vect.z;
#else
  virial     += tmp_virial;
#endif

}


/******************************************************************************
*
*  do_neightab - compute neighbor table
*
******************************************************************************/

void do_neightab(cell *p, cell *q, vektor pbc)
{
  int i, j, k;
  int jstart, jend;
  int q_typ, p_typ;
  vektor d, tmp_d;
  real *qptr, radius2, r2_short = r2_end;

  /* For each atom in first cell */
  for (i=0; i<p->n; ++i) {

    tmp_d.x = p->ort X(i) - pbc.x;
    tmp_d.y = p->ort Y(i) - pbc.y;
    tmp_d.z = p->ort Z(i) - pbc.z;
    p_typ   = p->sorte[i];

    jstart = (p==q ? i+1 : 0);
    qptr   = q->ort + DIM * jstart;
    
    /* For each atom in neighbouring cell */
    for (j = jstart; j < q->n; ++j) {

      q_typ = q->sorte[j];
      
      /* Calculate distance  */
      d.x = *qptr++ - tmp_d.x;
      d.y = *qptr++ - tmp_d.y;
      d.z = *qptr++ - tmp_d.z;

      radius2 = SPROD(d,d);

      if (0==radius2) { char msgbuf[256];
        sprintf(msgbuf,
                "Distance is zero: nrs=%d %d\norte: %f %f %f, %f %f %f\n",
                NUMMER(p,i),NUMMER(q,i),
                p->ort X(i),p->ort Y(i),p->ort Z(i),
                q->ort X(j),q->ort Y(j),q->ort Z(j));
        error(msgbuf);
      }

      /* make neighbor tables for TTBP */
      if (radius2 <= ttbp_r2_cut[p_typ][q_typ]) {

        neightab *neigh;
        real  *tmp_ptr;

        /* update neighbor table of particle i */
        neigh = p->neigh[i];
        if (neigh->n_max <= neigh->n) {
          error("neighbor table too small, increase ttbp_len");
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
        /* we do not need a neighbor table in buffer cells
        neigh = q->neigh[j];
        if (neigh->n_max <= neigh->n) {
          error("neighbor table too small, increase ttbp_len");
        }
        neigh->typ[neigh->n] = p_typ;
        neigh->cl [neigh->n] = p;
        neigh->num[neigh->n] = i;
        tmp_ptr  = &neigh->dist[3*neigh->n];
        *tmp_ptr = -d.x; ++tmp_ptr; 
        *tmp_ptr = -d.y; ++tmp_ptr; 
        *tmp_ptr = -d.z;
        neigh->n++;
        */
      }
    } /* for j */
  } /* for i */

#ifndef MONOLJ
  if (r2_short < r2_0) printf("\n Short distance! r2: %f\n",r2_short);
#endif

}

