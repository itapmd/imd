
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2006 Institute for Theoretical and Applied Physics,
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

        *Epot  += pot_zwi;

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
  int    i, j, k, p_typ, j_typ, knum, jnum, col, inc = ntypes * ntypes;
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

        /* FOURIER potential term */
        tmp_sp    = SPROD(d[j],d[k]);
        cos_theta = tmp_sp / (r[j] * r[k]);
        tmp       = cos_theta + 1.0 / ttbp_sp[p_typ];
        tmp_pot   = ttbp_constant[p_typ] * tmp * tmp;
        tmp_grad  = ttbp_constant[p_typ] * 2 * tmp;

        /* smoothing potential, total potential */
        tmp_f2   = pot[j]  * pot[k];
        pot_zwi  = tmp_pot * tmp_f2;
        *Epot   += pot_zwi;

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
      tmp_r  = 1.0 / ( r[j] - sw_a2[p_typ][j_typ] );
      fc[j]  = exp( sw_ga[p_typ][j_typ] * tmp_r );
      dfc[j] = - fc[j] * sw_ga[p_typ][j_typ] * tmp_r * tmp_r / r[j];
      
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

	pot_zwi   = sw_la[p_typ][j_typ][k_typ] * fc[j] * fc[k] * tmp * tmp;

	/* total potential */
	*Epot   += pot_zwi;

	/* forces */
	tmp_grad1  = sw_la[p_typ][j_typ][k_typ] * fc[j] * fc[k] * 2 * tmp;
	tmp_grad2  = sw_la[p_typ][j_typ][k_typ] * tmp * tmp;
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
