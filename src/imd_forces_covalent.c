
/*****************************************************************************
*
*  imd_forces_covalent.c -- force loops for many-body forces
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"
#include "potaccess.h"

#ifdef TTBP

/******************************************************************************
*
*  three body forces for TTBP, using neighbor tables computed in do_forces
*
******************************************************************************/

void do_forces2(cell *p)
{
  static real   *r2 = NULL, *r = NULL, *pot = NULL, *grad = NULL;
  static vektor *d  = NULL;
  neightab *neigh;
  vektor force_j, force_k;
  cell   *jcell, *kcell;
  int    i, j, k, p_typ, j_typ, knum, jnum, col, inc = ntypes * ntypes;
  int    is_short=0;
  real   *tmpptr;
  real   pot_zwi, tmp_pot, tmp_grad, tmp, tmp_j, tmp_k, r2_short=0.0;
  real   tmp_f2, cos_theta, tmp_sp;
  real   tmp_virial = 0.0;
#ifdef P_AXIAL
  vektor tmp_vir_vect = {0.0, 0.0, 0.0};
#endif

  d    = (vektor *) realloc( d,    neigh_len * sizeof(vektor) );
  r2   = (real *)   realloc( r2,   neigh_len * sizeof(real)   );
  r    = (real *)   realloc( r,    neigh_len * sizeof(real)   );
  pot  = (real *)   realloc( pot,  neigh_len * sizeof(real)   );
  grad = (real *)   realloc( grad, neigh_len * sizeof(real)   );
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

    p_typ   = SORTE(p,i);
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
      col = p_typ * ntypes + j_typ;
      PAIR_INT2(pot[j],grad[j],smooth_pot,col,inc,r2[j],r2_short,is_short)
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
        tot_pot_energy += pot_zwi;

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
        tmp_vir_vect.x += d[j].x * force_j.x + d[k].x * force_k.x;
        tmp_vir_vect.y += d[j].y * force_j.y + d[k].y * force_k.y;
        tmp_vir_vect.z += d[j].z * force_j.z + d[k].z * force_k.z;
#else
        tmp_virial     += SPROD(d[j],force_j) + SPROD(d[k],force_k);
#endif

    } /* neighbor pairs */

  } /* i */

#ifdef P_AXIAL
  vir_x  += tmp_vir_vect.x;
  virial += tmp_vir_vect.x;
  vir_y  += tmp_vir_vect.y;
  virial += tmp_vir_vect.y;
  vir_z  += tmp_vir_vect.z;
  virial += tmp_vir_vect.z;
#else
  virial += tmp_virial;
#endif

}

#endif /* TTBP */

#ifdef TERSOFF

/******************************************************************************
*
*  forces for Tersoff potential, using neighbor tables computed in do_forces
*
******************************************************************************/

void do_forces2(cell *p)
{
  static real   *r = NULL, *fc = NULL, *dfc = NULL;
  static vektor *d  = NULL;
  neightab *neigh;
  vektor dcos_j, dcos_k, dzeta_i, dzeta_j, force_j;
  static vektor *dzeta_k = NULL; 
  cell   *jcell, *kcell;
  int    i, j, k, p_typ, j_typ, k_typ, knum, jnum;
  real   *tmpptr;
  real   pot_zwi, tmp_grad;
  real   cos_theta, cut_tmp, cut_tmp_j;
  real   zeta, g_theta, b_ij;
  real   tmp_jk, tmp_j2, tmp_k2;
  real   phi_r, phi_a;
  real   tmp_1, tmp_2, tmp_3, tmp_4, tmp_5, tmp_6;
  real   tmp_virial = 0.0;
#ifdef P_AXIAL
  vektor tmp_vir_vect = {0.0, 0.0, 0.0};
#endif

  d    = (vektor *) realloc( d,    neigh_len * sizeof(vektor) );
  r    = (real *)   realloc( r,    neigh_len * sizeof(real)   );
  fc   = (real *)   realloc( fc,   neigh_len * sizeof(real)   );
  dfc  = (real *)   realloc( dfc,  neigh_len * sizeof(real)   );
  dzeta_k = (vektor *) realloc( dzeta_k,  neigh_len * sizeof(vektor) );
  if ((d==NULL) || (r==NULL) || (fc==NULL) || (dfc==NULL) || (dzeta_k==NULL))
    error("cannot allocate memory for temporary neighbor data");

  /*     k
          \
           \
	    i----j  */


  /* For each atom in cell */
  for (i=0; i<p->n; ++i) {

    p_typ   = SORTE(p,i);
    neigh   = p->neigh[i];

    /* construct some data for all neighbors */
    tmpptr = neigh->dist;
    for (j=0; j<neigh->n; ++j) {

      /* type, distance vector, radii */
      j_typ   = neigh->typ[j];
      d[j].x  = *tmpptr++;
      d[j].y  = *tmpptr++;
      d[j].z  = *tmpptr++;
      r[j]    = sqrt(SPROD(d[j],d[j]));

      /* cutoff function */
      cut_tmp   = M_PI / ( ter_r_cut[p_typ][j_typ] - ter_r0[p_typ][j_typ] );
      cut_tmp_j = cut_tmp * ( r[j] - ter_r0[p_typ][j_typ] );
      if ( r[j] < ter_r0[p_typ][j_typ] ) {
	fc[j]   = 1.0; 
	dfc[j]  = 0.0;
      }
      else if ( r[j] > ter_r_cut[p_typ][j_typ] ) {
	fc[j]   = 0.0;
	dfc[j]  = 0.0;
      }
      else {
	fc[j]   = 0.5 * ( 1.0 + cos( cut_tmp_j ) );
	dfc[j]  = - 0.5 * cut_tmp * sin( cut_tmp_j );
      }      
    } /* j */

    /* for each neighbor of i */
    for (j=0; j<neigh->n; ++j){

      j_typ = neigh->typ[j];

      zeta = 0.0;     
      dzeta_i.x = 0.0; dzeta_i.y = 0.0; dzeta_i.z = 0.0;
      dzeta_j.x = 0.0; dzeta_j.y = 0.0; dzeta_j.z = 0.0;
      
      /* for each neighbor of i other than j */
      for (k=0; k<neigh->n; ++k) 
	if (k!=j) {

	k_typ = neigh->typ[k];
  
	/* angular term */
	tmp_jk    = 1 / ( r[j] * r[k] );  
        cos_theta = SPROD(d[j],d[k]) * tmp_jk;
	tmp_1     = ters_h[p_typ] - cos_theta;
	tmp_2     = 1 / ( ter_d2[p_typ] + tmp_1 * tmp_1 );
	g_theta   = 1 + ter_c2[p_typ]/ter_d2[p_typ] - ter_c2[p_typ] * tmp_2;

	/* zeta */
	zeta  += fc[k] * ter_om[p_typ][k_typ] * g_theta; 

        tmp_j2 = cos_theta / ( r[j] * r[j] );
	tmp_k2 = cos_theta / ( r[k] * r[k] );

	/* derivatives of cos(theta), 
	   note that dcos_i + dcos_j + dcos_k = 0 */
        dcos_j.x = tmp_jk * d[k].x - tmp_j2 * d[j].x;
	dcos_j.y = tmp_jk * d[k].y - tmp_j2 * d[j].y;
	dcos_j.z = tmp_jk * d[k].z - tmp_j2 * d[j].z;

	dcos_k.x = tmp_jk * d[j].x - tmp_k2 * d[k].x;
	dcos_k.y = tmp_jk * d[j].y - tmp_k2 * d[k].y;
	dcos_k.z = tmp_jk * d[j].z - tmp_k2 * d[k].z;

	tmp_3    = 2 * ter_c2[p_typ] * tmp_1 * tmp_2 * tmp_2 * fc[k] * ter_om[p_typ][k_typ];
	tmp_grad = dfc[k] / r[k] * g_theta * ter_om[p_typ][k_typ];

	/* derivatives of zeta; dzeta_i is not the full derivative */
	dzeta_k[k].x = tmp_grad * d[k].x - tmp_3 * dcos_k.x;
	dzeta_k[k].y = tmp_grad * d[k].y - tmp_3 * dcos_k.y;
	dzeta_k[k].z = tmp_grad * d[k].z - tmp_3 * dcos_k.z;

	dzeta_i.x   -= dzeta_k[k].x;
	dzeta_i.y   -= dzeta_k[k].y;
	dzeta_i.z   -= dzeta_k[k].z;

	dzeta_j.x   -= tmp_3 * dcos_j.x;
	dzeta_j.y   -= tmp_3 * dcos_j.y;
	dzeta_j.z   -= tmp_3 * dcos_j.z;
 
      } /* k */

      phi_r  = 0.5 * ter_a[p_typ][j_typ] * exp(- ter_la[p_typ][j_typ] * r[j]);
      phi_a  = 0.5 * ter_b[p_typ][j_typ] * exp(- ter_mu[p_typ][j_typ] * r[j]);
      tmp_4  = pow( ters_ga[p_typ] * zeta, ters_n[p_typ] );

      b_ij  = ter_chi[p_typ][j_typ] * pow( 1 + tmp_4, -1 / ( 2 * ters_n[p_typ] ));

      pot_zwi         = phi_r - b_ij * phi_a;
      tot_pot_energy += fc[j] * pot_zwi;

      if ( zeta == 0.0 )   /* only one neighbor of i */
	tmp_5 = 0.0;
      else
        tmp_5 = - b_ij * fc[j] * phi_a * tmp_4 / ( 2 * zeta * ( 1 + tmp_4 ) );
      tmp_6   = ( fc[j] * (- phi_r * ter_la[p_typ][j_typ] + phi_a * ter_mu[p_typ][j_typ] * b_ij)  + dfc[j] * pot_zwi ) / r[j];

      /* tmp force on particle j */
      force_j.x = - tmp_6 * d[j].x + tmp_5 * dzeta_j.x;
      force_j.y = - tmp_6 * d[j].y + tmp_5 * dzeta_j.y;
      force_j.z = - tmp_6 * d[j].z + tmp_5 * dzeta_j.z;

      /* update force on particle k */
      for (k=0; k<neigh->n; ++k) 
	if (k!=j) {
        kcell = neigh->cl [k];
        knum  = neigh->num[k];
        kcell->kraft X(knum) += tmp_5 * dzeta_k[k].x;
        kcell->kraft Y(knum) += tmp_5 * dzeta_k[k].y;
        kcell->kraft Z(knum) += tmp_5 * dzeta_k[k].z;
#ifdef P_AXIAL
	tmp_vir_vect.x += tmp_5 * d[k].x * dzeta_k[k].x;
	tmp_vir_vect.y += tmp_5 * d[k].y * dzeta_k[k].y;
	tmp_vir_vect.z += tmp_5 * d[k].z * dzeta_k[k].z;
#else
	tmp_virial     += tmp_5 * SPROD(d[k],dzeta_k[k]);
#endif
#ifdef STRESS_TENS
	kcell->presstens X(knum)        -= d[k].x * tmp_5 * dzeta_k[k].x;
	kcell->presstens Y(knum)        -= d[k].y * tmp_5 * dzeta_k[k].y;
	kcell->presstens Z(knum)        -= d[k].z * tmp_5 * dzeta_k[k].z;
	kcell->presstens_offdia X(knum) -= ( d[k].y * tmp_5 * dzeta_k[k].z + d[k].z * tmp_5 * dzeta_k[k].y ) / 2;
	kcell->presstens_offdia Y(knum) -= ( d[k].z * tmp_5 * dzeta_k[k].x + d[k].x * tmp_5 * dzeta_k[k].z ) / 2;
	kcell->presstens_offdia Z(knum) -= ( d[k].x * tmp_5 * dzeta_k[k].y + d[k].y * tmp_5 * dzeta_k[k].x ) / 2;
#endif 
	}
      
      /* update force on particle j */
      jcell = neigh->cl [j];
      jnum  = neigh->num[j];
      jcell->kraft X(jnum) += force_j.x;
      jcell->kraft Y(jnum) += force_j.y;
      jcell->kraft Z(jnum) += force_j.z;
      jcell->pot_eng[jnum] += fc[j] * pot_zwi;

      /* update force on particle i */
      p->kraft X(i) += tmp_5 * dzeta_i.x - force_j.x;
      p->kraft Y(i) += tmp_5 * dzeta_i.y - force_j.y;
      p->kraft Z(i) += tmp_5 * dzeta_i.z - force_j.z;
      p->pot_eng[i] += fc[j] * pot_zwi;

#ifdef P_AXIAL
      tmp_vir_vect.x += d[j].x * force_j.x;
      tmp_vir_vect.y += d[j].y * force_j.y;
      tmp_vir_vect.z += d[j].z * force_j.z;
#else
      tmp_virial     += SPROD(d[j],force_j);
#endif
#ifdef STRESS_TENS
      jcell->presstens X(jnum)        -= d[j].x * force_j.x;
      jcell->presstens Y(jnum)        -= d[j].y * force_j.y;
      jcell->presstens Z(jnum)        -= d[j].z * force_j.z;
      jcell->presstens_offdia X(jnum) -= ( d[j].y * force_j.z + d[j].z * force_j.y ) / 2;
      jcell->presstens_offdia Y(jnum) -= ( d[j].z * force_j.x + d[j].x * force_j.z ) / 2;
      jcell->presstens_offdia Z(jnum) -= ( d[j].x * force_j.y + d[j].y * force_j.x ) / 2;
#endif

    } /* neighbor j */

  } /* i */

#ifdef P_AXIAL
  vir_x  += tmp_vir_vect.x;
  virial += tmp_vir_vect.x;
  vir_y  += tmp_vir_vect.y;
  virial += tmp_vir_vect.y;
  vir_z  += tmp_vir_vect.z;
  virial += tmp_vir_vect.z;
#else
  virial += tmp_virial;
#endif

}

#endif /* TERSOFF */

/******************************************************************************
*
*  do_neightab - compute neighbor table
*
******************************************************************************/

void do_neightab(cell *p, cell *q, vektor pbc)
{
  int i, j, k;
  int jstart, jend;
  int q_typ, p_typ, column;
  vektor d, tmp_d;
  real *qptr, radius2;

  /* For each atom in first cell */
  for (i=0; i<p->n; ++i) {

    tmp_d.x = p->ort X(i) - pbc.x;
    tmp_d.y = p->ort Y(i) - pbc.y;
    tmp_d.z = p->ort Z(i) - pbc.z;
    p_typ   = SORTE(p,i);

    jstart = (p==q ? i+1 : 0);
    qptr   = q->ort + DIM * jstart;
    
    /* For each atom in neighbouring cell */
    for (j = jstart; j < q->n; ++j) {

      q_typ = SORTE(q,j);
      
      /* Calculate distance  */
      d.x = *qptr++ - tmp_d.x;
      d.y = *qptr++ - tmp_d.y;
      d.z = *qptr++ - tmp_d.z;

      column  = p_typ * ntypes + q_typ;
      radius2 = SPROD(d,d);

      if (0==radius2) { char msgbuf[256];
        sprintf(msgbuf,
                "Distance is zero: nrs=%d %d\norte: %f %f %f, %f %f %f\n",
                NUMMER(p,i),NUMMER(q,i),
                p->ort X(i),p->ort Y(i),p->ort Z(i),
                q->ort X(j),q->ort Y(j),q->ort Z(j));
        error(msgbuf);
      }
#ifdef TTBP
      /* make neighbor tables for TTBP */
      if (radius2 <= smooth_pot.end[column])
#endif
#ifdef TERSOFF
      /* make neighbor tables for TERSOFF */
      if (radius2 <= ter_r2_cut[p_typ][q_typ])
#endif
      {        
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
        /* we do not need a neighbor table in buffer cells
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
        */
      }
    } /* for j */
  } /* for i */

}


#ifdef TERSOFF
void init_tersoff(void) {

  int i, j, n = 0;
  real tmp;

  /* parameters for more than one atom type */
  for (i=0; i<ntypes; i++) {
    ter_c2[i] = ters_c[i] * ters_c[i];
    ter_d2[i] = ters_d[i] * ters_d[i];
    for (j=i; j<ntypes; j++) {
      ter_r_cut[i][j]  = ter_r_cut[j][i]  = ters_r_cut[n];
      ter_r2_cut[i][j] = ter_r2_cut[j][i] = ter_r_cut[i][j] * ter_r_cut[i][j];
      ter_r0[i][j]     = ter_r0[j][i]     = ters_r0[n];
      ter_a[i][j]      = ter_a[j][i]      = ters_a[n];
      ter_b[i][j]      = ter_b[j][i]      = ters_b[n];
      ter_la[i][j]     = ter_la[j][i]     = ters_la[n];
      ter_mu[i][j]     = ter_mu[j][i]     = ters_mu[n];
      ++n;      
    }
  }
  for (i=0; i<ntypes; i++) ter_chi[i][i] = 1.0;
  if ( ntypes>1 ) {
    for (i=0; i<(ntypes-1); i++)
      for (j=(i+1); j<ntypes; j++) {
        ter_chi[i][j] = ter_chi[j][i] 
                      = ters_chi[i * ( 2 * ntypes - i - 3 ) / 2 + j - 1]; 
      }
  }
  for (i=0; i<ntypes; i++) ter_om[i][i] = 1.0;
  if ( ntypes>1 ) {
    for (i=0; i<(ntypes-1); i++)
      for (j=(i+1); j<ntypes; j++) {
        ter_om[i][j] = ter_om[j][i] 
                     = ters_om[i * ( 2 * ntypes - i - 3 ) / 2 + j - 1]; 
      }
  }  
  tmp = 0.0;
  for (i=0; i<ntypes; ++i)
    for (j=0; j<ntypes; ++j)
      tmp = MAX( tmp, ter_r2_cut[i][j] );
  cellsz = MAX(cellsz,tmp);
}

#endif 
