
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2010 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

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

#ifdef KEATING

/******************************************************************************
*
*  three body forces for Keating potential, 
*        using neighbor tables computed in do_forces
*
******************************************************************************/

void do_forces2(cell *p, real *Epot, real *Virial, 
                real *Vir_xx, real *Vir_yy, real *Vir_zz,
                real *Vir_yz, real *Vir_zx, real *Vir_xy)
{
  static vektor *d  = NULL;
  static int    curr_len = 0;
  neightab *neigh;
  vektor force_j, force_k;
  cell   *jcell, *kcell;
  int    i, j, k, p_typ, j_typ, k_typ, knum, jnum;
  real   *tmpptr;
  real   pot_zwi, tmp_sp, tmp_d2, tmp;
  real   tmp_virial = 0.0;
#ifdef P_AXIAL
  vektor tmp_vir_vect = {0.0, 0.0, 0.0};
#endif

  if (curr_len < neigh_len) {
    d = (vektor *) realloc( d, neigh_len * sizeof(vektor) );
    if ( d==NULL )
      error("cannot allocate memory for temporary neighbor data");
    curr_len = neigh_len;
  }

  /* For each atom in cell */
  for (i=0; i<p->n; ++i) {

    p_typ   = SORTE(p,i);
    neigh   = NEIGH(p,i);

    /* construct some data for all neighbors */
    tmpptr = neigh->dist;
    for (j=0; j<neigh->n; ++j) {

      /* type, distance vector */
      j_typ   = neigh->typ[j];
      d[j].x  = *tmpptr++;
      d[j].y  = *tmpptr++;
      d[j].z  = *tmpptr++;
    }

    /* for each pair of neighbors */
    for (j=0; j<neigh->n-1; ++j)
      for (k=j+1; k<neigh->n; ++k) {

	j_typ = neigh->typ[j];
	k_typ = neigh->typ[k];

        /* potential term */
        tmp_sp  = SPROD(d[j],d[k]);
	tmp_d2  = keat_d[p_typ][j_typ] * keat_d[p_typ][k_typ];
	pot_zwi = 3.0 * keat_beta[p_typ][j_typ][k_typ] / ( 8.0 * tmp_d2 )
	        * SQR( tmp_sp + tmp_d2 / 3.0 ); 

        *Epot   += pot_zwi;
        pot_zwi /= 3.0;  /* avoid triple counting */

        /* forces */
        tmp   = 3.0 * keat_beta[p_typ][j_typ][k_typ] / ( 4.0 * tmp_d2)
	      * ( SPROD(d[j],d[k]) + tmp_d2 / 3.0 );

        force_j.x = tmp * d[k].x;
        force_j.y = tmp * d[k].y;
        force_j.z = tmp * d[k].z;

        force_k.x = tmp * d[j].x;
        force_k.y = tmp * d[j].y;
        force_k.z = tmp * d[j].z;
        /* update force on particle i */
        KRAFT(p,i,X) += force_j.x + force_k.x;
        KRAFT(p,i,Y) += force_j.y + force_k.y;
        KRAFT(p,i,Z) += force_j.z + force_k.z;
        POTENG(p,i)  += pot_zwi;

        /* update force on particle j */
        jcell = (cell *) neigh->cl [j];
        jnum  = neigh->num[j];
        KRAFT(jcell,jnum,X) -= force_j.x;
        KRAFT(jcell,jnum,Y) -= force_j.y;
        KRAFT(jcell,jnum,Z) -= force_j.z;
        POTENG(jcell,jnum)  += pot_zwi;

        /* update force on particle k */
        kcell = (cell *) neigh->cl [k];
        knum  = neigh->num[k];
        KRAFT(kcell,knum,X) -= force_k.x;
        KRAFT(kcell,knum,Y) -= force_k.y;
        KRAFT(kcell,knum,Z) -= force_k.z;
        POTENG(kcell,knum)  += pot_zwi;

#ifdef P_AXIAL
        tmp_vir_vect.x += d[j].x * force_j.x + d[k].x * force_k.x;
        tmp_vir_vect.y += d[j].y * force_j.y + d[k].y * force_k.y;
        tmp_vir_vect.z += d[j].z * force_j.z + d[k].z * force_k.z;
#else
        tmp_virial     += SPROD(d[j],force_j) + SPROD(d[k],force_k);
#endif

#ifdef STRESS_TENS /* Distribute stress among atoms */
        if (do_press_calc) {
          tmp = 0.25 * ( force_j.x * d[j].x + force_k.x * d[k].x );
          PRESSTENS(p,i,xx)        += 2.0 * tmp;
          PRESSTENS(jcell,jnum,xx) += tmp;
          PRESSTENS(kcell,knum,xx) += tmp;
          tmp = 0.25 * ( force_j.y * d[j].y + force_k.y * d[k].y );
          PRESSTENS(p,i,yy)        += 2.0 * tmp;
          PRESSTENS(jcell,jnum,yy) += tmp;
          PRESSTENS(kcell,knum,yy) += tmp;
          tmp = 0.25 * ( force_j.z * d[j].z + force_k.z * d[k].z );
          PRESSTENS(p,i,zz)        += 2.0 * tmp;
          PRESSTENS(jcell,jnum,zz) += tmp;
          PRESSTENS(kcell,knum,zz) += tmp;
          tmp = 0.125 * ( force_j.y * d[j].z + force_k.y * d[k].z 
                        + force_j.z * d[j].y + force_k.z * d[k].y );
          PRESSTENS(p,i,yz)        += 2.0 * tmp;
          PRESSTENS(jcell,jnum,yz) += tmp;
          PRESSTENS(kcell,knum,yz) += tmp;	
          tmp =  0.125 * ( force_j.z * d[j].x + force_k.z * d[k].x 
                         + force_j.x * d[j].z + force_k.x * d[k].z );
          PRESSTENS(p,i,zx)        += 2.0 * tmp;
          PRESSTENS(jcell,jnum,zx) += tmp;
          PRESSTENS(kcell,knum,zx) += tmp;	
          tmp =  0.125 * ( force_j.x * d[j].y + force_k.x * d[k].y 
                         + force_j.y * d[j].x + force_k.y * d[k].x );
          PRESSTENS(p,i,xy)        += 2.0 * tmp;
          PRESSTENS(jcell,jnum,xy) += tmp;
          PRESSTENS(kcell,knum,xy) += tmp;
        }
#endif

    } /* neighbor pairs */

  } /* i */

#ifdef P_AXIAL
  *Vir_xx += tmp_vir_vect.x;
  *Virial += tmp_vir_vect.x;
  *Vir_yy += tmp_vir_vect.y;
  *Virial += tmp_vir_vect.y;
  *Vir_zz += tmp_vir_vect.z;
  *Virial += tmp_vir_vect.z;
#else
  *Virial += tmp_virial;
#endif

}

#endif

#ifdef TTBP

/******************************************************************************
*
*  three body forces for TTBP, using neighbor tables computed in do_forces
*
******************************************************************************/

void do_forces2(cell *p, real *Epot, real *Virial, 
                real *Vir_xx, real *Vir_yy, real *Vir_zz,
                real *Vir_yz, real *Vir_zx, real *Vir_xy)
{
  static real   *r2 = NULL, *r = NULL, *pot = NULL, *grad = NULL;
  static vektor *d  = NULL;
  static int    curr_len = 0;
  neightab *neigh;
  vektor force_j, force_k;
  cell   *jcell, *kcell;
  int    i, j, k, p_typ, j_typ, k_typ, knum, jnum, col, inc = ntypes * ntypes;
  int    is_short=0;
  real   *tmpptr;
  real   pot_zwi, tmp_pot, tmp_grad, tmp, tmp_j, tmp_k;
  real   tmp_f2, cos_theta, tmp_sp;
  real   tmp_virial = 0.0;
#ifdef P_AXIAL
  vektor tmp_vir_vect = {0.0, 0.0, 0.0};
#endif

  if (curr_len < neigh_len) {
    d    = (vektor *) realloc( d,    neigh_len * sizeof(vektor) );
    r2   = (real *)   realloc( r2,   neigh_len * sizeof(real)   );
    r    = (real *)   realloc( r,    neigh_len * sizeof(real)   );
    pot  = (real *)   realloc( pot,  neigh_len * sizeof(real)   );
    grad = (real *)   realloc( grad, neigh_len * sizeof(real)   );
    if ((d==NULL) || (r2==NULL) || (r==NULL) || (pot==NULL) || (grad==NULL))
      error("cannot allocate memory for temporary neighbor data");
    curr_len = neigh_len;
  }

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
    neigh   = NEIGH(p,i);

    /* shortcut for types without 3-body interactions */
    if (ttbp_constant[p_typ] == 0.0) continue;

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
      PAIR_INT(pot[j],grad[j],smooth_pot,col,inc,r2[j],is_short)
    }

    /* for each pair of neighbors */
    for (j=0; j<neigh->n-1; ++j)
      for (k=j+1; k<neigh->n; ++k) {
	if(ttbp_vas == 0 || (r[j] < ttbp_cut[0] && r[k] < ttbp_cut[0])) {
	k_typ   = neigh->typ[k];
        /* FOURIER potential term */
        tmp_sp    = SPROD(d[j],d[k]);
        cos_theta = tmp_sp / (r[j] * r[k]);
        tmp       = cos_theta + 1.0 / ttbp_sp[p_typ];
#ifdef XT
        tmp_pot   = ttbp_constant[p_typ] * g(cos_theta);
        tmp_grad  = ttbp_constant[p_typ] * dg(cos_theta);
#else
        tmp_pot   = B[j_typ][p_typ][k_typ] * ttbp_constant[p_typ] * tmp * tmp;
        tmp_grad  = B[j_typ][p_typ][k_typ] * ttbp_constant[p_typ] * 2 * tmp;
#endif
        /* smoothing potential, total potential */
        tmp_f2   = pot[j]  * pot[k];
        pot_zwi  = tmp_pot * tmp_f2;
        *Epot   += pot_zwi;
        pot_zwi /= 3.0;  /* avoid triple counting */

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
        KRAFT(p,i,X) += force_j.x + force_k.x;
        KRAFT(p,i,Y) += force_j.y + force_k.y;
        KRAFT(p,i,Z) += force_j.z + force_k.z;
        POTENG(p,i)  += pot_zwi;

        /* update force on particle j */
        jcell = (cell *) neigh->cl [j];
        jnum  = neigh->num[j];
        KRAFT(jcell,jnum,X) -= force_j.x;
        KRAFT(jcell,jnum,Y) -= force_j.y;
        KRAFT(jcell,jnum,Z) -= force_j.z;
        POTENG(jcell,jnum)   += pot_zwi;

        /* update force on particle k */
        kcell = (cell *) neigh->cl [k];
        knum  = neigh->num[k];
        KRAFT(kcell,knum,X) -= force_k.x;
        KRAFT(kcell,knum,Y) -= force_k.y;
        KRAFT(kcell,knum,Z) -= force_k.z;
        POTENG(kcell,knum)  += pot_zwi;

#ifdef P_AXIAL
        tmp_vir_vect.x += d[j].x * force_j.x + d[k].x * force_k.x;
        tmp_vir_vect.y += d[j].y * force_j.y + d[k].y * force_k.y;
        tmp_vir_vect.z += d[j].z * force_j.z + d[k].z * force_k.z;
#else
        tmp_virial     += SPROD(d[j],force_j) + SPROD(d[k],force_k);
#endif

#ifdef STRESS_TENS /* Distribute stress among atoms */
        if (do_press_calc) {
          tmp = 0.25 * ( force_j.x * d[j].x + force_k.x * d[k].x );
          PRESSTENS(p,i,xx)        += 2.0 * tmp;
          PRESSTENS(jcell,jnum,xx) += tmp;
          PRESSTENS(kcell,knum,xx) += tmp;
          tmp = 0.25 * ( force_j.y * d[j].y + force_k.y * d[k].y );
          PRESSTENS(p,i,yy)        += 2.0 * tmp;
          PRESSTENS(jcell,jnum,yy) += tmp;
          PRESSTENS(kcell,knum,yy) += tmp;
          tmp = 0.25 * ( force_j.z * d[j].z + force_k.z * d[k].z );
          PRESSTENS(p,i,zz)        += 2.0 * tmp;
          PRESSTENS(jcell,jnum,zz) += tmp;
          PRESSTENS(kcell,knum,zz) += tmp;
          tmp = 0.125 * ( force_j.y * d[j].z + force_k.y * d[k].z 
                        + force_j.z * d[j].y + force_k.z * d[k].y );
          PRESSTENS(p,i,yz)        += 2.0 * tmp;
          PRESSTENS(jcell,jnum,yz) += tmp;
          PRESSTENS(kcell,knum,yz) += tmp;	
          tmp =  0.125 * ( force_j.z * d[j].x + force_k.z * d[k].x 
                         + force_j.x * d[j].z + force_k.x * d[k].z );
          PRESSTENS(p,i,zx)        += 2.0 * tmp;
          PRESSTENS(jcell,jnum,zx) += tmp;
          PRESSTENS(kcell,knum,zx) += tmp;	
          tmp =  0.125 * ( force_j.x * d[j].y + force_k.x * d[k].y 
                         + force_j.y * d[j].x + force_k.y * d[k].x );
          PRESSTENS(p,i,xy)        += 2.0 * tmp;
          PRESSTENS(jcell,jnum,xy) += tmp;
          PRESSTENS(kcell,knum,xy) += tmp;
        }
#endif

	} 
    } /* neighbor pairs */

  } /* i */

  if (is_short==1) printf("\n Short distance!\n");

#ifdef P_AXIAL
  *Vir_xx += tmp_vir_vect.x;
  *Virial += tmp_vir_vect.x;
  *Vir_yy += tmp_vir_vect.y;
  *Virial += tmp_vir_vect.y;
  *Vir_zz += tmp_vir_vect.z;
  *Virial += tmp_vir_vect.z;
#else
  *Virial += tmp_virial;
#endif

}

#endif /* TTBP */

#ifdef STIWEB

/******************************************************************************
*
*  forces for Stillinger-Weber potential, 
*             using neighbor tables computed in do_forces
*
******************************************************************************/

void do_forces2(cell *p, real *Epot, real *Virial, 
                real *Vir_xx, real *Vir_yy, real *Vir_zz,
                real *Vir_yz, real *Vir_zx, real *Vir_xy)
{
  static vektor  *d = NULL;
  static real    *r = NULL, *fc = NULL, *dfc = NULL;
  static int     curr_len = 0;
  neightab *neigh;
  int      i, j, k, p_typ, k_typ, j_typ, jnum, knum, col;
  vektor   force_j, force_k;
  cell     *jcell, *kcell;
  real     *tmpptr;
  real     tmp_r, tmp_sp, cos_theta, tmp, pot_zwi;
  real     tmp_grad1, tmp_grad2, tmp_jj, tmp_kk, tmp_jk;
  real     tmp_1, tmp_2;
  real     tmp_virial = 0.0;
#ifdef P_AXIAL
  vektor   tmp_vir_vect = {0.0, 0.0, 0.0};
#endif

  if (curr_len < neigh_len) {
    d    = (vektor *) realloc( d,    neigh_len * sizeof(vektor) );
    r    = (real *)   realloc( r,    neigh_len * sizeof(real)   );
    fc   = (real *)   realloc( fc,   neigh_len * sizeof(real)   );
    dfc  = (real *)   realloc( dfc,  neigh_len * sizeof(real)   );
    if ( (d==NULL) || (r==NULL) || (fc==NULL) || (dfc==NULL) )
      error("cannot allocate memory for temporary neighbor data");
    curr_len = neigh_len;
  }

  /* For each atom in cell */
  for (i=0; i<p->n; ++i) {

    p_typ   = SORTE(p,i);
    neigh   = NEIGH(p,i);
    
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
      tmp_r = r[j] - sw_a2[p_typ][j_typ];
      if (tmp_r < -0.01 * sw_ga[p_typ][j_typ]) {
        tmp_r  = 1.0 / tmp_r;
        fc [j] = exp( sw_ga[p_typ][j_typ] * tmp_r );
        dfc[j] = - fc[j] * sw_ga[p_typ][j_typ] * tmp_r * tmp_r / r[j];
      } else {
        fc [j] = 0.0;
        dfc[j] = 0.0;
      }      
    } /* j */

    /* for each pair of neighbors */
    for (j=0; j<neigh->n-1; ++j)
      for (k=j+1; k<neigh->n; ++k) {

	j_typ = neigh->typ[j];
	jcell = (cell *) neigh->cl [j];
	jnum  = neigh->num[j];
	k_typ = neigh->typ[k];
	kcell = (cell *) neigh->cl [k];
	knum  = neigh->num[k];

        /* shortcut for types without 3-body interactions */
        if (sw_la[p_typ][j_typ][k_typ] == 0.0) continue;

	/* potential term */
	tmp_sp    = SPROD(d[j],d[k]);
	cos_theta = tmp_sp / (r[j] * r[k]);
	tmp       = cos_theta + 1.0 / 3.0;
#ifdef TERNBCC
	pot_zwi   = sw_la[p_typ][j_typ][k_typ] * fc[j] * fc[k] * g(cos_theta);
#else
	pot_zwi   = sw_la[p_typ][j_typ][k_typ] * fc[j] * fc[k] * tmp * tmp;
#endif

	/* total potential */
	*Epot   += pot_zwi;
        pot_zwi /= 3.0;  /* avoid triple counting */

	/* forces */
#ifdef TERNBCC
	tmp_grad1  = sw_la[p_typ][j_typ][k_typ] * fc[j] * fc[k] * dg(cos_theta);
	tmp_grad2  = sw_la[p_typ][j_typ][k_typ] * g(cos_theta);
#else
	tmp_grad1  = sw_la[p_typ][j_typ][k_typ] * fc[j] * fc[k] * 2 * tmp;
	tmp_grad2  = sw_la[p_typ][j_typ][k_typ] * tmp * tmp;
#endif
	tmp_jj = 1.0 / ( r[j] * r[j] );
	tmp_kk = 1.0 / ( r[k] * r[k] );
	tmp_jk = 1.0 / ( r[j] * r[k] );
	tmp_1 = tmp_grad2 * dfc[j] * fc[k] - tmp_grad1 * cos_theta * tmp_jj;
	tmp_2 = tmp_grad1 * tmp_jk;

	force_j.x = tmp_1 * d[j].x + tmp_2 * d[k].x;
       	force_j.y = tmp_1 * d[j].y + tmp_2 * d[k].y;
	force_j.z = tmp_1 * d[j].z + tmp_2 * d[k].z;

	tmp_1 = tmp_grad2 * dfc[k] * fc[j] - tmp_grad1 * cos_theta * tmp_kk;
	force_k.x = tmp_1 * d[k].x + tmp_2 * d[j].x;
	force_k.y = tmp_1 * d[k].y + tmp_2 * d[j].y;
	force_k.z = tmp_1 * d[k].z + tmp_2 * d[j].z;

	/* update force on particle i */
	KRAFT(p,i,X) += force_j.x + force_k.x;
	KRAFT(p,i,Y) += force_j.y + force_k.y;
	KRAFT(p,i,Z) += force_j.z + force_k.z;
	POTENG(p,i)  += pot_zwi;
	
	/* update force on particle j */
	KRAFT(jcell,jnum,X) -= force_j.x;
	KRAFT(jcell,jnum,Y) -= force_j.y;
	KRAFT(jcell,jnum,Z) -= force_j.z;
	POTENG(jcell,jnum)  += pot_zwi;

	/* update force on particle k */
	KRAFT(kcell,knum,X) -= force_k.x;
	KRAFT(kcell,knum,Y) -= force_k.y;
	KRAFT(kcell,knum,Z) -= force_k.z;
	POTENG(kcell,knum)  += pot_zwi;
#ifdef P_AXIAL
	tmp_vir_vect.x += d[j].x * force_j.x + d[k].x * force_k.x;
	tmp_vir_vect.y += d[j].y * force_j.y + d[k].y * force_k.y;
	tmp_vir_vect.z += d[j].z * force_j.z + d[k].z * force_k.z;
#else
	tmp_virial     += SPROD(d[j],force_j) + SPROD(d[k],force_k);
#endif

#ifdef STRESS_TENS /* Distribute stress among atoms */
        if (do_press_calc) {
          tmp = 0.25 * ( force_j.x * d[j].x + force_k.x * d[k].x );
          PRESSTENS(p,i,xx)        += 2.0 * tmp;
          PRESSTENS(jcell,jnum,xx) += tmp;
          PRESSTENS(kcell,knum,xx) += tmp;
          tmp = 0.25 * ( force_j.y * d[j].y + force_k.y * d[k].y );
          PRESSTENS(p,i,yy)        += 2.0 * tmp;
          PRESSTENS(jcell,jnum,yy) += tmp;
          PRESSTENS(kcell,knum,yy) += tmp;
          tmp = 0.25 * ( force_j.z * d[j].z + force_k.z * d[k].z );
          PRESSTENS(p,i,zz)        += 2.0 * tmp;
          PRESSTENS(jcell,jnum,zz) += tmp;
          PRESSTENS(kcell,knum,zz) += tmp;
          tmp = 0.125 * ( force_j.y * d[j].z + force_k.y * d[k].z 
                        + force_j.z * d[j].y + force_k.z * d[k].y );
          PRESSTENS(p,i,yz)        += 2.0 * tmp;
          PRESSTENS(jcell,jnum,yz) += tmp;
          PRESSTENS(kcell,knum,yz) += tmp;	
          tmp =  0.125 * ( force_j.z * d[j].x + force_k.z * d[k].x 
                         + force_j.x * d[j].z + force_k.x * d[k].z );
          PRESSTENS(p,i,zx)        += 2.0 * tmp;
          PRESSTENS(jcell,jnum,zx) += tmp;
          PRESSTENS(kcell,knum,zx) += tmp;	
          tmp =  0.125 * ( force_j.x * d[j].y + force_k.x * d[k].y 
                         + force_j.y * d[j].x + force_k.y * d[k].x );
          PRESSTENS(p,i,xy)        += 2.0 * tmp;
          PRESSTENS(jcell,jnum,xy) += tmp;
          PRESSTENS(kcell,knum,xy) += tmp;
        }
#endif

      } /* neighbor pairs */
  
  } /* i */

#ifdef P_AXIAL
  *Vir_xx += tmp_vir_vect.x;
  *Virial += tmp_vir_vect.x;
  *Vir_yy += tmp_vir_vect.y;
  *Virial += tmp_vir_vect.y;
  *Vir_zz += tmp_vir_vect.z;
  *Virial += tmp_vir_vect.z;
#else
  *Virial += tmp_virial;
#endif

}

#endif /* STIWEB */

#if defined (TERNBCC) || defined(XT)

/******************************************************************************
*
*  g(cos_theta) for ternary BCC potential
*
******************************************************************************/

real g(real cos_theta){
  real gc;
  if (                         cos_theta < -5. / 6. ) 
    gc=             (cos_theta + 1. ) * (cos_theta + 1. );
  if (cos_theta >= -5. / 6. && cos_theta < -3. / 6. ) 
    gc= 1. / 18. - (cos_theta + 2. / 3. ) * (cos_theta + 2. / 3. );
  if (cos_theta >= -3. / 6. && cos_theta < -1. /6. )
    gc=            (cos_theta + 1. / 3. ) * (cos_theta + 1. / 3. );
  if (cos_theta >= -1. / 6. && cos_theta <  1. /6. )
    gc= 1. / 18. -  cos_theta * cos_theta ;
  if (cos_theta >=  1. / 6. )
    gc=            (cos_theta - 1. / 3. ) * (cos_theta - 1. / 3. );
  return gc;
/*  return (cos(3./2.*3.14159265*cos_theta)); */
/* return cos_theta + 1. /3. ; */
}

/******************************************************************************
*
*  dg(cos_theta) for ternary BCC potential
*
******************************************************************************/

real dg(real cos_theta){
  real dgc;
  if (                         cos_theta < -5. / 6. ) 
    dgc=             cos_theta + 1. ;
  if (cos_theta >= -5. / 6. && cos_theta < -3. / 6. ) 
    dgc=            -(cos_theta + 2. /3. );
  if (cos_theta >= -3. / 6. && cos_theta < -1. / 6. )
    dgc=             cos_theta + 1. /3. ;
  if (cos_theta >= -1. / 6. && cos_theta <  1. / 6. )
    dgc=            -cos_theta;
  if (cos_theta >=  1. / 6. )
    dgc=             cos_theta - 1. / 3. ;
  return dgc;
/*  return (-3./4.*3.14159265*sin(3.*3.14159625*cos_theta)); */
/* return cos_theta + 1. / 3.; */
}
#endif

#ifdef TERSOFF

/******************************************************************************
*
*  forces for Tersoff potential, using neighbor tables computed in do_forces
*
******************************************************************************/

void do_forces2(cell *p, real *Epot, real *Virial, 
                real *Vir_xx, real *Vir_yy, real *Vir_zz,
                real *Vir_yz, real *Vir_zx, real *Vir_xy)
{
  static real   *r = NULL, *fc = NULL, *dfc = NULL;
  static vektor *d = NULL;
  static int    curr_len = 0;
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
  real   tmp, tmp_1, tmp_2, tmp_3, tmp_4, tmp_5, tmp_6;
  real   tmp_virial = 0.0;
#ifdef P_AXIAL
  vektor tmp_vir_vect = {0.0, 0.0, 0.0};
#endif

  if (curr_len < neigh_len) {
    d    = (vektor *) realloc( d,    neigh_len * sizeof(vektor) );
    r    = (real *)   realloc( r,    neigh_len * sizeof(real)   );
    fc   = (real *)   realloc( fc,   neigh_len * sizeof(real)   );
    dfc  = (real *)   realloc( dfc,  neigh_len * sizeof(real)   );
    dzeta_k = (vektor *) realloc( dzeta_k,  neigh_len * sizeof(vektor) );
    if ((d==NULL) || (r==NULL) || (fc==NULL) || (dfc==NULL) || (dzeta_k==NULL))
      error("cannot allocate memory for temporary neighbor data");
    curr_len = neigh_len;
  }

  /*     k
          \
           \
	    i----j  */

  /* For each atom in cell */
  for (i=0; i<p->n; ++i) {

    p_typ   = SORTE(p,i);
    neigh   = NEIGH(p,i);

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
    for (j=0; j<neigh->n; ++j) {

      j_typ = neigh->typ[j];

      /* shortcut for types without 3-body interactions */
      if (ter_b[p_typ][j_typ] == 0.0) continue;

      jcell = (cell *) neigh->cl [j];
      jnum  = neigh->num[j];

      zeta = 0.0;     
      dzeta_i.x = 0.0; dzeta_i.y = 0.0; dzeta_i.z = 0.0;
      dzeta_j.x = 0.0; dzeta_j.y = 0.0; dzeta_j.z = 0.0;
      
      /* for each neighbor of i other than j */
      for (k=0; k<neigh->n; ++k) {
 
	if (k==j) continue;

	k_typ = neigh->typ[k];
  
	/* angular term */
	tmp_jk    = 1 / ( r[j] * r[k] );  
        cos_theta = SPROD(d[j],d[k]) * tmp_jk;
#ifdef TERSOFF2
	tmp_1     = ter_h[p_typ][j_typ] - cos_theta;
	tmp_2     = 1 / ( ter_d2[p_typ][j_typ] + tmp_1 * tmp_1 );
	g_theta   = 1 + ter_c2[p_typ][j_typ]/ter_d2[p_typ][j_typ] 
	  - ter_c2[p_typ][j_typ] * tmp_2;
#else
	tmp_1     = ters_h[p_typ] - cos_theta;
	tmp_2     = 1 / ( ter_d2[p_typ] + tmp_1 * tmp_1 );
	g_theta   = 1 + ter_c2[p_typ]/ter_d2[p_typ] - ter_c2[p_typ] * tmp_2;
#endif
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

#ifdef TERSOFF2
	tmp_3    = 2 * ter_c2[p_typ][j_typ] * tmp_1 * tmp_2 * tmp_2 
	  * fc[k] * ter_om[p_typ][k_typ];
#else
	tmp_3    = 2 * ter_c2[p_typ] * tmp_1 * tmp_2 * tmp_2 
	  * fc[k] * ter_om[p_typ][k_typ];
#endif
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

      phi_a = 0.5 * ter_b[p_typ][j_typ] * exp(- ter_mu[p_typ][j_typ] * r[j]);
#ifdef TERSOFF2
      tmp_4 = pow( ter_ga[p_typ][j_typ] * zeta, ter_n[p_typ][j_typ] );

      b_ij  = pow( 1 + tmp_4, -1 / ( 2 * ter_n[p_typ][j_typ] ));
#else
      tmp_4 = pow( ters_ga[p_typ] * zeta, ters_n[p_typ] );

      b_ij  = pow( 1 + tmp_4, -1 / ( 2 * ters_n[p_typ] ));
#endif
      pot_zwi  = - b_ij * phi_a;
      *Epot   += fc[j] * pot_zwi;

      if ( zeta == 0.0 )   /* only one neighbor of i */
	tmp_5 = 0.0;
      else
        tmp_5 = - b_ij * fc[j] * phi_a * tmp_4 / ( 2 * zeta * ( 1 + tmp_4 ) );
      tmp_6   = ( fc[j] * phi_a * ter_mu[p_typ][j_typ] * b_ij  
                  + dfc[j] * pot_zwi ) / r[j];

      /* tmp force on particle j */
      force_j.x = - tmp_6 * d[j].x + tmp_5 * dzeta_j.x;
      force_j.y = - tmp_6 * d[j].y + tmp_5 * dzeta_j.y;
      force_j.z = - tmp_6 * d[j].z + tmp_5 * dzeta_j.z;

      /* update force on particle k */
      for (k=0; k<neigh->n; ++k) 
	if (k!=j) {
        kcell = (cell *) neigh->cl [k];
        knum  = neigh->num[k];
        KRAFT(kcell,knum,X) += tmp_5 * dzeta_k[k].x;
        KRAFT(kcell,knum,Y) += tmp_5 * dzeta_k[k].y;
        KRAFT(kcell,knum,Z) += tmp_5 * dzeta_k[k].z;

#ifdef P_AXIAL
	tmp_vir_vect.x += tmp_5 * d[k].x * dzeta_k[k].x;
	tmp_vir_vect.y += tmp_5 * d[k].y * dzeta_k[k].y;
	tmp_vir_vect.z += tmp_5 * d[k].z * dzeta_k[k].z;
#else
	tmp_virial     += tmp_5 * SPROD(d[k],dzeta_k[k]);
#endif
#ifdef STRESS_TENS
        if (do_press_calc) {
	  tmp = 0.5 * d[k].x * tmp_5 * dzeta_k[k].x;
	  PRESSTENS(p,i,xx)        -= tmp;
  	  PRESSTENS(jcell,jnum,xx) -= tmp;
	  tmp = 0.5 * d[k].y * tmp_5 * dzeta_k[k].y;
	  PRESSTENS(p,i,yy)        -= tmp;
  	  PRESSTENS(jcell,jnum,yy) -= tmp;
	  tmp = 0.5 * d[k].z * tmp_5 * dzeta_k[k].z;
	  PRESSTENS(p,i,zz)        -= tmp;
  	  PRESSTENS(jcell,jnum,zz) -= tmp;
	  tmp = 0.25 * ( d[k].y * tmp_5 * dzeta_k[k].z + 
                         d[k].z * tmp_5 * dzeta_k[k].y );
	  PRESSTENS(p,i,yz)        -= tmp;
  	  PRESSTENS(jcell,jnum,yz) -= tmp;
	  tmp = 0.25 * ( d[k].z * tmp_5 * dzeta_k[k].x + 
                         d[k].x * tmp_5 * dzeta_k[k].z );
	  PRESSTENS(p,i,zx)        -= tmp;
  	  PRESSTENS(jcell,jnum,zx) -= tmp;	  
	  tmp = 0.25 * ( d[k].x * tmp_5 * dzeta_k[k].y + 
			 d[k].y * tmp_5 * dzeta_k[k].x );
	  PRESSTENS(p,i,xy)        -= tmp;
  	  PRESSTENS(jcell,jnum,xy) -= tmp;
	}
#endif 
      }
      
      pot_zwi *= 0.5;  /* avoid double counting */

      /* update force on particle j */
      KRAFT(jcell,jnum,X) += force_j.x;
      KRAFT(jcell,jnum,Y) += force_j.y;
      KRAFT(jcell,jnum,Z) += force_j.z;
      POTENG(jcell,jnum)  += fc[j] * pot_zwi;

      /* update force on particle i */
      KRAFT(p,i,X) += tmp_5 * dzeta_i.x - force_j.x;
      KRAFT(p,i,Y) += tmp_5 * dzeta_i.y - force_j.y;
      KRAFT(p,i,Z) += tmp_5 * dzeta_i.z - force_j.z;
      POTENG(p,i)  += fc[j] * pot_zwi;

#ifdef P_AXIAL
      tmp_vir_vect.x += d[j].x * force_j.x;
      tmp_vir_vect.y += d[j].y * force_j.y;
      tmp_vir_vect.z += d[j].z * force_j.z;
#else
      tmp_virial     += SPROD(d[j],force_j);
#endif
#ifdef STRESS_TENS
      if (do_press_calc) {
	tmp = 0.5 * d[j].x * force_j.x; 
	PRESSTENS(p,i,xx)        -= tmp;
	PRESSTENS(jcell,jnum,xx) -= tmp;
	tmp = 0.5 * d[j].y * force_j.y;
	PRESSTENS(p,i,yy)        -= tmp;
	PRESSTENS(jcell,jnum,yy) -= tmp;
	tmp = 0.5 * d[j].z * force_j.z;
	PRESSTENS(p,i,zz)        -= tmp;
	PRESSTENS(jcell,jnum,zz) -= tmp;
	tmp = 0.25 * ( d[j].y * force_j.z + d[j].z*force_j.y );
	PRESSTENS(p,i,yz)        -= tmp;
	PRESSTENS(jcell,jnum,yz) -= tmp;
	tmp = 0.25 * ( d[j].z * force_j.x + d[j].x*force_j.z );
	PRESSTENS(p,i,zx)        -= tmp;
	PRESSTENS(jcell,jnum,zx) -= tmp;
	tmp = 0.25 * ( d[j].x * force_j.y + d[j].y*force_j.x );
	PRESSTENS(p,i,xy)        -= tmp;
	PRESSTENS(jcell,jnum,xy) -= tmp;
      }
#endif

    } /* neighbor j */

  } /* i */

#ifdef P_AXIAL
  *Vir_xx += tmp_vir_vect.x;
  *Virial += tmp_vir_vect.x;
  *Vir_yy += tmp_vir_vect.y;
  *Virial += tmp_vir_vect.y;
  *Vir_zz += tmp_vir_vect.z;
  *Virial += tmp_vir_vect.z;
#else
  *Virial += tmp_virial;
#endif

}

#endif /* TERSOFF */



#ifdef TERSOFFMOD

/******************************************************************************
*
*  forces for modified Tersoff potential, 
*  using neighbor tables computed in do_forces
*
******************************************************************************/

void do_forces2(cell *p, real *Epot, real *Virial, 
                real *Vir_xx, real *Vir_yy, real *Vir_zz,
                real *Vir_yz, real *Vir_zx, real *Vir_xy)
{

  static real *r = NULL, *fc = NULL, *dfc = NULL;
  static vektor *er = NULL;
  static int curr_len = 0;

  neightab *neigh;
  cell *jcell, *kcell;

  real *tmpptr;
  real tmp, tmp1, tmp2, tmp3, tmp4;
  real cos_theta, g, dg_cos, zeta, dzeta_cos, dzeta_ij, dzeta_ik, b, phi;

  vektor dcos_j, dcos_k, gradi_zeta, gradj_zeta, force_j;
  static vektor *gradk_zeta = NULL;
  
  real tmp_virial = 0.0;
#ifdef P_AXIAL
  vektor tmp_vir_vect = {0.0, 0.0, 0.0};
#endif


  // memory allocation
  if (curr_len < neigh_len) {
    er   = (vektor *) realloc( er,   neigh_len * sizeof(vektor) );
    r    = (real *)   realloc( r,    neigh_len * sizeof(real)   );
    fc   = (real *)   realloc( fc,   neigh_len * sizeof(real)   );
    dfc  = (real *)   realloc( dfc,  neigh_len * sizeof(real)   );
    gradk_zeta = (vektor *) realloc( gradk_zeta,  neigh_len * sizeof(vektor) );
    if ((er==NULL) || (r==NULL) || (fc==NULL) || (dfc==NULL) || (gradk_zeta==NULL))
      error("cannot allocate memory for temporary neighbor data");
    curr_len = neigh_len;
  }

  int i, j, k, p_typ, j_typ, k_typ, knum, jnum, i_id, j_id;
  
  /* For each atom in cell */
  for (i=0; i<p->n; ++i) {

    p_typ   = SORTE(p,i);
    neigh   = NEIGH(p,i);
    i_id = NUMMER(p,i);

    /* construct some data for all neighbors */
    tmpptr = neigh->dist;

    for (j=0; j<neigh->n; ++j) {

      /* type, distance vector, radii */
      j_typ   = neigh->typ[j];
      er[j].x  = *tmpptr++;
      er[j].y  = *tmpptr++;
      er[j].z  = *tmpptr++;
      r[j]    = sqrt(SPROD(er[j],er[j]));
      er[j].x /= r[j]; er[j].y /= r[j]; er[j].z /= r[j];

      /* cutoff function */
      tmp1   = M_PI / ( ter_r_cut[p_typ][j_typ] - ter_r0[p_typ][j_typ] );
      tmp2   = tmp1 * ( r[j] - ter_r0[p_typ][j_typ] );
      if ( r[j] < ter_r0[p_typ][j_typ] ) {
        fc[j]   = 1.0; 
        dfc[j]  = 0.0;
      }else if ( r[j] > ter_r_cut[p_typ][j_typ] ) {
        fc[j]   = 0.0;
        dfc[j]  = 0.0;
      }else {
        //fc[j]   = 0.5 * ( 1.0 + cos( tmp2 ) );
        //dfc[j]  = - 0.5 * tmp1 * sin( tmp2 );
        fc[j]   = 0.5 * ( 1.0 + 1.125*cos(tmp2) - 0.125*cos(3.0*tmp2) );
        dfc[j]  = - 0.5 * tmp1 * ( 1.125*sin(tmp2) - 0.375*sin(3.0*tmp2) );
      }
    } /* j */

    /* for each neighbor of i */
    for (j=0; j<neigh->n; ++j) {

      j_typ = neigh->typ[j];

      /* shortcut for types without 3-body interactions */
      //if (ter_b[p_typ][j_typ] == 0.0) continue;

      jcell = (cell *) neigh->cl [j];
      jnum  = neigh->num[j];
      j_id = NUMMER(jcell,jnum);

      zeta  = 0.0; dzeta_ij = 0.0;
      gradi_zeta.x = 0.0; gradi_zeta.y = 0.0; gradi_zeta.z = 0.0;
      gradj_zeta.x = 0.0; gradj_zeta.y = 0.0; gradj_zeta.z = 0.0;

      /* for each neighbor of i other than j */
      for (k=0; k<neigh->n; ++k) {
        if (k==j) continue;

        //k_typ = neigh->typ[k];

        cos_theta = SPROD(er[j],er[k]);

        /* derivatives of cos(theta), 
          note that dcos_i + dcos_j + dcos_k = 0 */
        dcos_j.x = (er[k].x - er[j].x*cos_theta)/r[j];
        dcos_j.y = (er[k].y - er[j].y*cos_theta)/r[j];
        dcos_j.z = (er[k].z - er[j].z*cos_theta)/r[j];

        dcos_k.x = (er[j].x - er[k].x*cos_theta)/r[k];
        dcos_k.y = (er[j].y - er[k].y*cos_theta)/r[k];
        dcos_k.z = (er[j].z - er[k].z*cos_theta)/r[k];
        /* -------------------------------------*/
#ifdef TERSOFFMOD2

        /* angular-dependent term g_jik and its derivative */

        /* argument */
        tmp1 = ter_h[p_typ][j_typ] - cos_theta;
        /* denominator */
        tmp2 = 1.0/(ter_c3[p_typ][j_typ]+tmp1*tmp1);
        /* exp. term */
        tmp3 = ter_c4[p_typ][j_typ] * exp(-ter_c5[p_typ][j_typ]*tmp1*tmp1);
        /* g_ijk */
        g = ter_c1[p_typ][j_typ] + ter_c2[p_typ][j_typ]*tmp1*tmp1*tmp2*(1.0+tmp3);
        /* dg/dcos */
        dg_cos = 2.0*ter_c2[p_typ][j_typ]*tmp1*tmp2*
              (ter_c5[p_typ][j_typ]*tmp1*tmp1*tmp3 - ter_c3[p_typ][j_typ]*tmp2*(1.0+tmp3));


        /* zeta term and its derivatives */

        tmp1 = r[j] - r[k];
        tmp2 = ter_alpha[p_typ][j_typ]*ter_beta[p_typ][j_typ]*pow(tmp1, ter_beta[p_typ][j_typ] - 1);
        tmp3 = exp(ter_alpha[p_typ][j_typ] * pow(tmp1, ter_beta[p_typ][j_typ]));

        dzeta_ik  = (dfc[k]-fc[k]*tmp2)*g*tmp3;
        
#else

        /* angular-dependent term g_jik and its derivative */

        /* argument */
        tmp1 = ters_h[p_typ] - cos_theta;
        /* denominator */
        tmp2 = 1.0/(ters_c3[p_typ]+tmp1*tmp1);
        /* exp. term */
        tmp3 = ters_c4[p_typ] * exp(-ters_c5[p_typ]*tmp1*tmp1);
        /* g_ijk */
        g = ters_c1[p_typ] + ters_c2[p_typ]*tmp1*tmp1*tmp2*(1.0+tmp3);
        /* dg/dcos */
        dg_cos = 2.0*ters_c2[p_typ]*tmp1*tmp2*
              (ters_c5[p_typ]*tmp1*tmp1*tmp3 - ters_c3[p_typ]*tmp2*(1.0+tmp3));


        /* zeta term and its derivatives */

        tmp1 = r[j] - r[k];
        tmp2 = ters_alpha[p_typ]*ters_beta[p_typ]*pow(tmp1, ters_beta[p_typ] - 1);
        tmp3 = exp(ters_alpha[p_typ] * pow(tmp1, ters_beta[p_typ]));

        dzeta_ik  = (dfc[k]-fc[k]*tmp2)*g*tmp3;
#endif
        tmp4 = fc[k]*g*tmp3;

        zeta  += tmp4;
        dzeta_ij += tmp4*tmp2;
        dzeta_cos = fc[k]*dg_cos*tmp3;

        gradk_zeta[k].x = dzeta_cos*dcos_k.x + dzeta_ik*er[k].x;
        gradk_zeta[k].y = dzeta_cos*dcos_k.y + dzeta_ik*er[k].y;
        gradk_zeta[k].z = dzeta_cos*dcos_k.z + dzeta_ik*er[k].z;

        gradi_zeta.x   -= gradk_zeta[k].x;
        gradi_zeta.y   -= gradk_zeta[k].y;
        gradi_zeta.z   -= gradk_zeta[k].z;

        gradj_zeta.x   += dzeta_cos*dcos_j.x;
        gradj_zeta.y   += dzeta_cos*dcos_j.y;
        gradj_zeta.z   += dzeta_cos*dcos_j.z;
 
      }

      gradj_zeta.x   += dzeta_ij*er[j].x;
      gradj_zeta.y   += dzeta_ij*er[j].y;
      gradj_zeta.z   += dzeta_ij*er[j].z;

#ifdef TERSOFFMOD2
      tmp1 = pow(zeta, ter_eta[p_typ][j_typ]);
      b = pow(1. + tmp1, -ter_delta[p_typ][j_typ]);

      tmp2 = 0.5*b*ter_b[p_typ][j_typ]*exp(-ter_mu[p_typ][j_typ]*r[j]);

      if ( zeta == 0.0 )   /* only one neighbor of i */
        tmp3 = 0.0;
      else
        tmp3 = tmp2*fc[j]*ter_eta[p_typ][j_typ]*ter_delta[p_typ][j_typ]*tmp1/((1.+tmp1)*zeta);
#else
      tmp1 = pow(zeta, ters_eta[p_typ]);
      b = pow(1. + tmp1, -ters_delta[p_typ]);

      tmp2 = 0.5*b*ter_b[p_typ][j_typ]*exp(-ter_mu[p_typ][j_typ]*r[j]);

      if ( zeta == 0.0 )   /* only one neighbor of i */
        tmp3 = 0.0;
      else
        tmp3 = tmp2*fc[j]*ters_eta[p_typ]*ters_delta[p_typ]*tmp1/((1.+tmp1)*zeta);
#endif

      /* attractive part */
      phi = -tmp2;
      tmp4 = -tmp2*(dfc[j]-ter_mu[p_typ][j_typ]*fc[j]);

      /* repulsive part */
      tmp1 = 0.5*ter_a[p_typ][j_typ]*exp(-ter_la[p_typ][j_typ]*r[j]);
      phi += tmp1;
      tmp4 += tmp1*(dfc[j]-ter_la[p_typ][j_typ]*fc[j]);

      *Epot += fc[j] * phi;
      phi *= 0.5;  /* avoid double counting */

      force_j.x = -tmp4*er[j].x - tmp3*gradj_zeta.x;
      force_j.y = -tmp4*er[j].y - tmp3*gradj_zeta.y;
      force_j.z = -tmp4*er[j].z - tmp3*gradj_zeta.z;

      /* update force on particle j */
      KRAFT(jcell,jnum,X) += force_j.x;
      KRAFT(jcell,jnum,Y) += force_j.y;
      KRAFT(jcell,jnum,Z) += force_j.z;
      POTENG(jcell,jnum)  += fc[j]*phi; 
      /* update force on particle i */
      KRAFT(p,i,X) -= tmp3*gradi_zeta.x + force_j.x;
      KRAFT(p,i,Y) -= tmp3*gradi_zeta.y + force_j.y;
      KRAFT(p,i,Z) -= tmp3*gradi_zeta.z + force_j.z;
      POTENG(p,i)  += fc[j]*phi;
      
#ifdef P_AXIAL
      tmp_vir_vect.x += r[j]*er[j].x * force_j.x;
      tmp_vir_vect.y += r[j]*er[j].y * force_j.y;
      tmp_vir_vect.z += r[j]*er[j].z * force_j.z;
#else
      tmp_virial     += r[j]*SPROD(er[j],force_j);
#endif
#ifdef STRESS_TENS
      if (do_press_calc) {
        tmp = 0.5 * r[j]*er[j].x * force_j.x; 
        PRESSTENS(p,i,xx)        += tmp;
        PRESSTENS(jcell,jnum,xx) += tmp;
        tmp = 0.5 * r[j]*er[j].y * force_j.y;
        PRESSTENS(p,i,yy)        += tmp;
        PRESSTENS(jcell,jnum,yy) += tmp;
        tmp = 0.5 * r[j]*er[j].z * force_j.z;
        PRESSTENS(p,i,zz)        += tmp;
        PRESSTENS(jcell,jnum,zz) += tmp;
        tmp = 0.25 * r[j]*( er[j].y * force_j.z + er[j].z*force_j.y );
        PRESSTENS(p,i,yz)        += tmp;
        PRESSTENS(jcell,jnum,yz) += tmp;
        tmp = 0.25 * r[j]*( er[j].z * force_j.x + er[j].x*force_j.z );
        PRESSTENS(p,i,zx)        += tmp;
        PRESSTENS(jcell,jnum,zx) += tmp;
        tmp = 0.25 * r[j]*( er[j].x * force_j.y + er[j].y*force_j.x );
        PRESSTENS(p,i,xy)        += tmp;
        PRESSTENS(jcell,jnum,xy) += tmp;
      }
#endif


      /* update force on particle k */
      for (k=0; k<neigh->n; ++k) 
        if (k!=j) {
          kcell = (cell *) neigh->cl [k];
          knum  = neigh->num[k];
          KRAFT(kcell,knum,X) -= tmp3 * gradk_zeta[k].x;
          KRAFT(kcell,knum,Y) -= tmp3 * gradk_zeta[k].y;
          KRAFT(kcell,knum,Z) -= tmp3 * gradk_zeta[k].z;

#ifdef P_AXIAL
          tmp_vir_vect.x -= tmp3 * r[k] * er[k].x * gradk_zeta[k].x;
          tmp_vir_vect.y -= tmp3 * r[k] * er[k].y * gradk_zeta[k].y;
          tmp_vir_vect.z -= tmp3 * r[k] * er[k].z * gradk_zeta[k].z;
#else
          tmp_virial     -= tmp3 * r[k] * SPROD(er[k],gradk_zeta[k]);
#endif
#ifdef STRESS_TENS
          if (do_press_calc) {
            tmp = 0.5 * r[k] * er[k].x * tmp3 * gradk_zeta[k].x;
            PRESSTENS(p,i,xx)        -= tmp;
            PRESSTENS(jcell,jnum,xx) -= tmp;
            tmp = 0.5 * r[k] * er[k].y * tmp3 * gradk_zeta[k].y;
            PRESSTENS(p,i,yy)        -= tmp;
            PRESSTENS(jcell,jnum,yy) -= tmp;
            tmp = 0.5 * r[k] * er[k].z * tmp3 * gradk_zeta[k].z;
            PRESSTENS(p,i,zz)        -= tmp;
            PRESSTENS(jcell,jnum,zz) -= tmp;
            tmp = 0.25 *r[k]*tmp3*( er[k].y * gradk_zeta[k].z + er[k].z * gradk_zeta[k].y );
            PRESSTENS(p,i,yz)        -= tmp;
            PRESSTENS(jcell,jnum,yz) -= tmp;
            tmp = 0.25 *r[k]*tmp3*( er[k].z * gradk_zeta[k].x + er[k].x * gradk_zeta[k].z );
            PRESSTENS(p,i,zx)        -= tmp;
            PRESSTENS(jcell,jnum,zx) -= tmp;    
            tmp = 0.25 *r[k]*tmp3*( er[k].x * gradk_zeta[k].y + er[k].y * gradk_zeta[k].x );
            PRESSTENS(p,i,xy)        -= tmp;
            PRESSTENS(jcell,jnum,xy) -= tmp;
          }
#endif 
        }
    }
  }

#ifdef P_AXIAL
  *Vir_xx += tmp_vir_vect.x;
  *Virial += tmp_vir_vect.x;
  *Vir_yy += tmp_vir_vect.y;
  *Virial += tmp_vir_vect.y;
  *Vir_zz += tmp_vir_vect.z;
  *Virial += tmp_vir_vect.z;
#else
  *Virial += tmp_virial;
#endif

}

#endif /* TERSOFFMOD */

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

    tmp_d.x = ORT(p,i,X) - pbc.x;
    tmp_d.y = ORT(p,i,Y) - pbc.y;
    tmp_d.z = ORT(p,i,Z) - pbc.z;
    p_typ   = SORTE(p,i);

#ifdef TWOD
    jstart = (((p==q) && (pbc.x==0) && (pbc.y==0))               ? i+1 : 0);
#else
    jstart = (((p==q) && (pbc.x==0) && (pbc.y==0) && (pbc.z==0)) ? i+1 : 0);
#endif
    qptr   = &ORT(q,jstart,X);
    
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
                ORT(p,i,X), ORT(p,i,Y), ORT(p,i,Z),
                ORT(q,j,X), ORT(q,j,Y), ORT(q,j,Z) );
        error(msgbuf);
      }

      /* make neighbor tables for covalent systems */
      if (radius2 <= neightab_r2cut[column]) {        
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
        /* we do not need a neighbor table in buffer cells
        neigh = q->neigh[j];
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
        */
      }
    } /* for j */
  } /* for i */

}

#ifdef KEATING

/******************************************************************************
*
*  init_keating
*
******************************************************************************/

void init_keating(void) {

  int  i, j, k, n, m;
  real tmp;

  /* parameters for more than one atom type */
  n = 0; m = 0;
  for (i=0; i<ntypes; i++) 
    for (j=i; j<ntypes; j++) {
      keat_alpha[i][j]  = keat_alpha[j][i]  = keating_alpha[n];
      keat_d[i][j]      = keat_d[j][i]      = keating_d[n];
      keat_r_cut[i][j]  = keat_r_cut[j][i]  = keating_r_cut[n];
      keat_r2_cut[i][j] = keat_r2_cut[j][i] = SQR( keat_r_cut[i][j] );
      n++; 
      for (k=0; k<ntypes; k++) {
	keat_beta[k][i][j] = keat_beta[k][j][i] = keating_beta[m];
	m++;
      }
    }
  tmp = 0.0;
  for (i=0; i<ntypes; ++i)
    for (j=0; j<ntypes; ++j)
      tmp = MAX( tmp, keat_r2_cut[i][j] );
  cellsz = MAX(cellsz,tmp);

  /* update neighbor table cutoff */
  if (NULL==neightab_r2cut) {
    neightab_r2cut = (real *) calloc( ntypes * ntypes, sizeof(real) );
    if (NULL==neightab_r2cut) 
      error("cannot allocate memory for neightab_r2cut");
  }
  for (i=0; i<ntypes; i++)
    for (j=0; j<ntypes; j++)
      neightab_r2cut[i*ntypes+j] = 
	MAX( neightab_r2cut[i*ntypes+j], keat_r2_cut[i][j] );
}

#endif /* KEATING */

#ifdef TTBP

/******************************************************************************
*
*  init_ttbp
*
******************************************************************************/

void init_ttbp(void) {
  int  i;
  /* update neighbor table cutoff */
  if (NULL==neightab_r2cut) {
    neightab_r2cut = (real *) calloc( ntypes * ntypes, sizeof(real) );
    if (NULL==neightab_r2cut) 
      error("cannot allocate memory for neightab_r2cut");
  }
  for (i=0; i<ntypes*ntypes; i++)
    neightab_r2cut[i] = MAX( neightab_r2cut[i], smooth_pot.end[i] );
  if (ttbp_vas == 0){
    for (i=0; i<8; ++i)
      ttbp_constant2[i] = 1.;
  }
  B[0][0][0] = ttbp_constant2[0];
  B[0][0][1] = ttbp_constant2[1];
  B[0][1][0] = ttbp_constant2[2];
  B[0][1][1] = ttbp_constant2[3];
  B[1][0][0] = ttbp_constant2[4];
  B[1][0][1] = ttbp_constant2[5];
  B[1][1][0] = ttbp_constant2[6];
  B[1][1][1] = ttbp_constant2[7];
}

#endif /* TTBP */

#ifdef STIWEB

/******************************************************************************
*
*  init_stiweb
*
******************************************************************************/

void init_stiweb(void) {

  int  i, j;

  /* update neighbor table cutoff */
  if (NULL==neightab_r2cut) {
    neightab_r2cut = (real *) calloc( ntypes * ntypes, sizeof(real) );
    if (NULL==neightab_r2cut) 
      error("cannot allocate memory for neightab_r2cut");
  }
  for (i=0; i<ntypes; i++)
    for (j=0; j<ntypes; j++) {
      real tmp = SQR(sw_a2[i][j]);
      neightab_r2cut[i*ntypes+j] = MAX( neightab_r2cut[i*ntypes+j], tmp );
      cellsz = MAX( cellsz, tmp );
    }
}

#endif /* STIWEB */

#ifdef TERSOFF

/******************************************************************************
*
*  init_tersoff
*
******************************************************************************/

void init_tersoff(void) {

  int i, j, n = 0;
  real tmp;

  /* parameters for more than one atom type */
  for (i=0; i<ntypes; i++) {
#ifndef TERSOFF2
    ter_c2[i] = ters_c[i] * ters_c[i];
    ter_d2[i] = ters_d[i] * ters_d[i];
#endif
    for (j=i; j<ntypes; j++) {
      ter_r_cut[i][j]  = ter_r_cut[j][i]  = ters_r_cut[n];
      ter_r2_cut[i][j] = ter_r2_cut[j][i] = ter_r_cut[i][j] * ter_r_cut[i][j];
      ter_r0[i][j]     = ter_r0[j][i]     = ters_r0[n];
      ter_a[i][j]      = ter_a[j][i]      = ters_a[n];
      ter_b[i][j]      = ter_b[j][i]      = ters_b[n];
      ter_la[i][j]     = ter_la[j][i]     = ters_la[n];
      ter_mu[i][j]     = ter_mu[j][i]     = ters_mu[n];
#ifdef TERSOFF2
      ter_ga[i][j]     = ter_ga[j][i]     = ters_ga[n];
      ter_n[i][j]      = ter_n[j][i]      = ters_n[n];
      ter_c[i][j]      = ter_c[j][i]      = ters_c[n];
      ter_c2[i][j]     = ter_c2[j][i]     = ter_c[i][j] * ter_c[i][j];
      ter_d[i][j]      = ter_d[j][i]      = ters_d[n];
      ter_d2[i][j]     = ter_d2[j][i]     = ter_d[i][j] * ter_d[i][j];
      ter_h[i][j]      = ter_h[j][i]      = ters_h[n];
#endif
      ++n;      
    }
  }

  /* absorb ters_chi into ter_b */
  if ( ntypes>1 ) {
    for (i=0; i<(ntypes-1); i++)
      for (j=(i+1); j<ntypes; j++) {
        ter_b[i][j] *= ters_chi[i * ( 2 * ntypes - i - 3 ) / 2 + j - 1]; 
        ter_b[j][i] *= ters_chi[i * ( 2 * ntypes - i - 3 ) / 2 + j - 1]; 
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

  /* update neighbor table cutoff */
  if (NULL==neightab_r2cut) {
    neightab_r2cut = (real *) calloc( ntypes * ntypes, sizeof(real) );
    if (NULL==neightab_r2cut) 
      error("cannot allocate memory for neightab_r2cut");
  }
  for (i=0; i<ntypes; i++)
    for (j=0; j<ntypes; j++)
      neightab_r2cut[i*ntypes+j] = 
	MAX( neightab_r2cut[i*ntypes+j], ter_r2_cut[i][j] );
}

#endif /* TERSOFF */

#ifdef TERSOFFMOD

/******************************************************************************
*
*  init_tersoffmod
*
******************************************************************************/

void init_tersoffmod(void) {

	int i, j, n = 0;

	for (i=0; i<ntypes; i++) {
		for (j=i; j<ntypes; j++) {
			ter_r_cut[i][j]  = ter_r_cut[j][i]  = ters_r_cut[n];
			ter_r2_cut[i][j] = ter_r2_cut[j][i] = ter_r_cut[i][j] * ter_r_cut[i][j];
			ter_r0[i][j]     = ter_r0[j][i]     = ters_r0[n];
			ter_a[i][j]      = ter_a[j][i]      = ters_a[n];
			ter_b[i][j]      = ter_b[j][i]      = ters_b[n];
			ter_la[i][j]     = ter_la[j][i]     = ters_la[n];
			ter_mu[i][j]     = ter_mu[j][i]     = ters_mu[n];
#ifdef TERSOFFMOD2
			ter_eta[i][j]    = ter_eta[j][i]    = ters_eta[n];
			ter_delta[i][j]  = ter_delta[j][i]  = ters_delta[n];
			ter_alpha[i][j]  = ter_alpha[j][i]  = ters_alpha[n];
			ter_beta[i][j]   = ter_beta[j][i]   = ters_beta[n];
			ter_c1[i][j]     = ter_c1[j][i]     = ters_c1[n];
			ter_c2[i][j]     = ter_c2[j][i]     = ters_c2[n];
			ter_c3[i][j]     = ter_c3[j][i]     = ters_c3[n];
			ter_c4[i][j]     = ter_c4[j][i]     = ters_c4[n];
			ter_c5[i][j]     = ter_c5[j][i]     = ters_c5[n];
			ter_h[i][j]      = ter_h[j][i]      = ters_h[n];
#endif
			++n;      
		}
	}


	real tmp = 0.0;
	for (i=0; i<ntypes; ++i)
		for (j=0; j<ntypes; ++j)
			tmp = MAX( tmp, ter_r2_cut[i][j] );
	cellsz = MAX(cellsz,tmp);

	/* update neighbor table cutoff */
	if (NULL==neightab_r2cut) {
		neightab_r2cut = (real *) calloc( ntypes * ntypes, sizeof(real) );
		if (NULL==neightab_r2cut) 
			error("cannot allocate memory for neightab_r2cut");
	}
	for (i=0; i<ntypes; i++)
		for (j=0; j<ntypes; j++)
			neightab_r2cut[i*ntypes+j] = 
				MAX( neightab_r2cut[i*ntypes+j], ter_r2_cut[i][j] );
}

#endif /* TERSOFFMOD */

#ifdef NNBR_TABLE
/**************************************************************/
/* computing a complete nearest neighbor table                */
/* which is used in processing steps (ADA, NYETENSOR, ..)     */
/* if already computed by COVALENT, does nothing              */
/* uses 'nbr_r2cut' as the threshold (constant for all atoms) */
/**************************************************************/
void do_neightab_complete() {
#ifndef COVALENT
	vektor pbc = {0.0, 0.0, 0.0}; /* atoms in buffer cells have pbc applied */
	int i, k, x,y,z;
	minicell *cell;
	if (nnbr_done == 0) { /*Check if already computed in this step*/
#ifdef MPI
	if (1==parallel_output) fix_cells();
#endif
		/*null the old list*/
		for (k = 0; k < nallcells; k++) {
			minicell *p = cell_array + k;
			for (i = 0; i < p->n; i++) {
				NEIGH(p, i)->n = 0;
			}
		}

		sync_cells(copy_cell,pack_cell,unpack_cell);

		/*
		 * constant cell offset to 14 neighbors
		 * the complete neighbor builder creates even neighbor-lists
		 * between two buffers cells to ensure correctness for
		 * neighbors of neighbors (provided the is cutoff chosen large enough)
		 */
		int cOff[14][3];
		cOff[0][0]  = 0;	cOff[0][1]  =  0;	cOff[0][2]  =  0;
		cOff[1][0]  = 0;	cOff[1][1]  =  0;	cOff[1][2]  =  1;
		cOff[2][0]  = 0;	cOff[2][1]  =  1;	cOff[2][2]  = -1;
		cOff[3][0]  = 0;	cOff[3][1]  =  1;	cOff[3][2]  =  0;
		cOff[4][0]  = 0;	cOff[4][1]  =  1;	cOff[4][2]  =  1;
		cOff[5][0]  = 1;	cOff[5][1]  = -1;	cOff[5][2]  = -1;
		cOff[6][0]  = 1;	cOff[6][1]  = -1;	cOff[6][2]  =  0;
		cOff[7][0]  = 1;	cOff[7][1]  = -1;	cOff[7][2]  =  1;
		cOff[8][0]  = 1;	cOff[8][1]  =  0;	cOff[8][2]  = -1;
		cOff[9][0]  = 1;	cOff[9][1]  =  0;	cOff[9][2]  =  0;
		cOff[10][0] = 1;	cOff[10][1] =  0;	cOff[10][2] =  1;
		cOff[11][0] = 1;	cOff[11][1] =  1;	cOff[11][2] = -1;
		cOff[12][0] = 1;	cOff[12][1] =  1;	cOff[12][2] =  0;
		cOff[13][0] = 1;	cOff[13][1] =  1;	cOff[13][2] =  1;

		for (x=0;x<cell_dim.x; x++){
			for (y=0;y<cell_dim.y; y++){
				for (z=0;z<cell_dim.z; z++){
					cell = PTR_3D_V(cell_array, x, y, z, cell_dim);
					for (i=0; i<14; i++){
						int x2, y2, z2;
						x2 = x + cOff[i][0];
						y2 = y + cOff[i][1];
						z2 = z + cOff[i][2];

						if (x2 >= 0 && x2 < cell_dim.x && y2 >= 0 && y2 < cell_dim.y && z2 >= 0 && z2 < cell_dim.z){
#ifndef NBLIST
							ivektor offset;
#ifdef LOADBALANCE
							offset = lb_cell_offset;
							offset.x+=1;
							offset.y+=1;
							offset.z+=1;
#else
							offset.x = (cell_dim.x-2) * my_coord.x;
							offset.y = (cell_dim.y-2) * my_coord.y;
							offset.z = (cell_dim.z-2) * my_coord.z;
#endif
							/*Handling of the cells along periodic boundaries*/
							/*Test if the absolute cell index is along the boundaries*/
							pbc.x = 0.0; pbc.y = 0.0; pbc.z = 0.0;

							if (pbc_dirs.x) {
								if (x+offset.x==1 && x2+offset.x==0)
									pbc.x = -(box_x.x + box_y.x + box_z.x);
								else if (x+offset.x==0 && x2+offset.x==1)
									pbc.x = +(box_x.x + box_y.x + box_z.x);
								else if (x+offset.x==global_cell_dim.x+1 && x2+offset.x==global_cell_dim.x)
									pbc.x = -(box_x.x + box_y.x + box_z.x);
								else if (x+offset.x==global_cell_dim.x && x2+offset.x==global_cell_dim.x+1)
									pbc.x = +(box_x.x + box_y.x + box_z.x);
							}

							if (pbc_dirs.y) {
								if (y+offset.y==1 && y2+offset.y==0)
									pbc.y = -(box_x.y + box_y.y + box_z.y);
								else if (y+offset.y==0 && y2+offset.y==1)
									pbc.y = +(box_x.y + box_y.y + box_z.y);
								else if (y+offset.y==global_cell_dim.y+1 && y2+offset.y==global_cell_dim.y)
									pbc.y = -(box_x.y + box_y.y + box_z.y);
								else if (y+offset.y==global_cell_dim.y && y2+offset.y==global_cell_dim.y+1)
									pbc.y = +(box_x.y + box_y.y + box_z.y);
							}

							if (pbc_dirs.z) {
								if (z+offset.z==1 && z2+offset.z==0)
									pbc.z = -(box_x.z + box_y.z + box_z.z);
								else if (z+offset.z==0 && z2+offset.z==1)
									pbc.z = +(box_x.z + box_y.z + box_z.z);
								else if (z+offset.z==global_cell_dim.z+1 && z2+offset.z==global_cell_dim.z)
									pbc.z = -(box_x.z + box_y.z + box_z.z);
								else if (z+offset.z==global_cell_dim.z && z2+offset.z==global_cell_dim.z+1)
									pbc.z = +(box_x.z + box_y.z + box_z.z);
							}
#endif

							do_neightab2(cell, PTR_3D_V(cell_array, x2, y2, z2, cell_dim), pbc);
						}
					}
				}
			}
		}
		nnbr_done = 1;
	}
#endif /*COVALENT*/
}



/********************************************************************/
/* inserts neighboring atoms in the neighbor table of the two cells */
/* called within do_neightab_complete()						    	*/
/********************************************************************/
void do_neightab2(cell *p, cell *q, vektor pbc)
{
  int i, j;
  int jstart;
  int q_typ, p_typ;
  vektor d, tmp_d;
  real *qptr, radius2;

  /* For each atom in first cell */
  for (i=0; i<p->n; ++i) {

    tmp_d.x = ORT(p,i,X) - pbc.x;
    tmp_d.y = ORT(p,i,Y) - pbc.y;
    tmp_d.z = ORT(p,i,Z) - pbc.z;
	p_typ   = SORTE(p,i);

    jstart = (((p==q) && (pbc.x==0) && (pbc.y==0) && (pbc.z==0)) ? i+1 : 0);
    qptr   = &ORT(q,jstart,X);

    /* For each atom in neighboring cell */
    for (j = jstart; j < q->n; ++j) {
      /* Calculate distance  */
      d.x = *qptr++ - tmp_d.x;
      d.y = *qptr++ - tmp_d.y;
      d.z = *qptr++ - tmp_d.z;

      radius2 = SPROD(d,d);

      /* make neighbor tables*/
      if (radius2 <= ada_nbr_r2cut) {
        neightab *neigh;
        real  *tmp_ptr;

        /* update neighbor table of particle i */
        neigh = NEIGH(p,i);
        if (neigh->n_max <= neigh->n) {
          increase_neightab( neigh, neigh->n_max + NEIGH_LEN_INC );
        }
        neigh->typ[neigh->n] = SORTE(q,j);
        neigh->cl [neigh->n] = q;
        neigh->num[neigh->n] = j;
        tmp_ptr  = &neigh->dist[DIM*neigh->n];
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
        tmp_ptr  = &neigh->dist[DIM*neigh->n];
        *tmp_ptr = -d.x; ++tmp_ptr;
        *tmp_ptr = -d.y; ++tmp_ptr;
        *tmp_ptr = -d.z;
        neigh->n++;
      }
    } /* for j */
  } /* for i */

}

#endif
