
/******************************************************************************
*
*  force loops for EAM2
*
******************************************************************************/

/******************************************************************************
*  $Revision$
*  $Date$
******************************************************************************/

#include "imd.h"
#include "potaccess.h"

/* Personal Debug Switch 
#define DEBUG_INFO 1000 
#define NR 94
*/

/*#define ATOMNR 100  quick hack to print the position of a special atom
		      better (later) Nr. in parameterfile... */ 

/******************************************************************************
*
*  do_forces for EAM2
*
*  first loop, calculates forces and energy for the core-core potential, 
*  and the embedding electron density at the atoms' sites
*
*  doesn't use 'actio = reactio' !!!!
*
******************************************************************************/

void do_forces(cell *p, cell *q, vektor pbc, real *Epot, 
               real *Virial, real *Vir_x, real *Vir_y, real *Vir_z)     
{
  int i,j,k;
  vektor d;
  vektor tmp_d;
  vektor force;
  real r2, rho_h;
  real tmp_virial = 0.0;
#ifdef P_AXIAL
  vektor tmp_vir_vect = {0.0, 0.0, 0.0};
#endif
  real pot_zwi, pot_grad;
  int  col, inc = ntypes * ntypes, is_short=0;
  int  q_typ,p_typ;
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
    qptr    = q->ort;
    
    for (j=0; j<q->n; ++j) {

      /* Calculate distance  */
      d.x = tmp_d.x - *qptr;  ++qptr;
      d.y = tmp_d.y - *qptr;  ++qptr;
#ifndef TWOD
      d.z = tmp_d.z - *qptr;  ++qptr;
#endif

      /* don't compute self interaction */
      if (!(p==q && i==j)) {

        q_typ = SORTE(q,j);
        col   = p_typ * ntypes + q_typ;
        r2    = SPROD(d,d);

#ifdef DEBUG_INFO
	/* check if distance == 0 */
	if (0==r2) { char msgbuf[256];
          sprintf(msgbuf, "Distance is zero between particles %d and %d!\n",
                  NUMMER(p,i), NUMMER(q,j));
          error(msgbuf);
	}
#endif

        if (r2 < core_pot.end[col]) {

          PAIR_INT(pot_zwi, pot_grad, core_pot, col, inc, r2, is_short);

	  /* Store forces in temp */
	  force.x = d.x * pot_grad;
	  force.y = d.y * pot_grad;
#ifndef TWOD
	  force.z = d.z * pot_grad;
#endif
	  
          /* Accumulate forces due to core-core potential */
	  p->kraft X(i) -= force.x;
	  p->kraft Y(i) -= force.y;
#ifndef TWOD
	  p->kraft Z(i) -= force.z;
#endif

	  p->pot_eng[i] += pot_zwi * 0.5;  /* avoid double counting */
	  *Epot         += pot_zwi * 0.5;  /* each pair occurs twice */
          /* virial and pressure tensor are updated twice for each pair */
          pot_grad      *= 0.5;

#ifdef DEBUG_INFO
	  if (p->nummer[i]==NR){
	    printf("d.x: %lf d.y: %lf d.z: %lf -> r: %lf\n",
                   d.x,d.y,d.z,eam2_r);
	    printf("-Phi'(r): %lf -> force: %.16lf %.16lf %.16lf\n",
		   eam2_r* pot_grad,force.x, force.y,force.z);
	    printf("Ges. Kraft: f.x: %.16lf f.y:%.16lf f.z:%.16lf\n",
		   p->kraft X(i),p->kraft Y(i),p->kraft Z(i));
	    printf("p_typ: %d, q_typ: %d  p->pot_eng[i]: %.12lf\n",
                   p_typ,q_typ,p->pot_eng[i]);
	  }
#endif
	  
#ifdef P_AXIAL
	  tmp_vir_vect.x += d.x * d.x * pot_grad;
	  tmp_vir_vect.y += d.y * d.y * pot_grad;
#ifndef TWOD
	  tmp_vir_vect.z += d.z * d.z * pot_grad;
#endif
#else
	  tmp_virial     += r2  * pot_grad;  
#endif

#ifdef STRESS_TENS
	  p->presstens X(i) += d.x * d.x * pot_grad;
	  p->presstens Y(i) += d.y * d.y * pot_grad;
#ifdef TWOD
	  p->presstens_offdia[i] += d.x * d.y * pot_grad;
#else
	  p->presstens Z(i) += d.z * d.z * pot_grad;
	  p->presstens_offdia X(i) += d.y * d.z * pot_grad;
	  p->presstens_offdia Y(i) += d.z * d.x * pot_grad;
	  p->presstens_offdia Z(i) += d.x * d.y * pot_grad;
#endif
#endif
	} /* core potential */

        if (r2 < rho_h_tab.end[col]) {
          VAL_FUNC(rho_h, rho_h_tab, col, inc, r2, is_short);
          p->eam2_rho_h[i] += rho_h; 
        }

      } /* if ! q==p.. */

#ifdef DEBUG_INFO
      else{
	if (p->nummer[i]==NR)
	  printf("++++++ Ort von Teilchen %d : x: %lf y: %lf z:%lf \n",
		 p->nummer[i],tmp_d.x,tmp_d.y,tmp_d.z);
      }
#endif

    } /* for j */
  } /* for i */

  /* print warning if short distance occurred */
  if (is_short==1) fprintf(stderr, "\n Short distance!\n");
  
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

} /* do_forces */ 


/******************************************************************************
*
*  second force loop, calculates the force and the energy 
*  caused by the embedding electron density
*  uses Phi(r2), Rho(r2), F(rho) and it's derivatives
*
******************************************************************************/

void do_forces_eam2(cell *p, cell *q, vektor pbc, real *Epot, 
                    real *Virial, real *Vir_x, real *Vir_y, real *Vir_z)
{
  int i,j,k;
  vektor d;
  vektor tmp_d;
  vektor force;
  real r2;
  int  is_short=0, idummy=0;
  real tmp_virial=0.0;
  int q_typ,p_typ;
  real *qptr; 
#ifdef P_AXIAL
  vektor tmp_vir_vect = {0.0, 0.0, 0.0};
#endif
  real dummy, eam2_energy, eam2_force;
  real f_i_strich, f_j_strich;
  real rho_i_strich, rho_j_strich;

  /* For each atom in first cell */
  for (i=0; i<p->n; ++i) {
    tmp_d.x = p->ort X(i) - pbc.x;
    tmp_d.y = p->ort Y(i) - pbc.y;
#ifndef TWOD
    tmp_d.z = p->ort Z(i) - pbc.z;
#endif
    p_typ   = SORTE(p,i);

    if (p==q) {
      /* f_i and f_i_strich -- TEST CUTOFF?? */
      PAIR_INT(eam2_energy, f_i_strich, embed_pot, p_typ, 
               ntypes, p->eam2_rho_h[i], idummy);
      /* add energy only once per particle (p==q) */
      p->pot_eng[i]  += eam2_energy;
      *Epot          += eam2_energy;
    } else {      
      /* only f_i_strich -- TEST CUTOFF?? */
      DERIV_FUNC(f_i_strich, embed_pot, p_typ, 
                 ntypes, p->eam2_rho_h[i], idummy);
    }

    /* go over ALL (except q==p && i==j) atoms, no "Actio = Reactio" */
    qptr = q->ort;
    for (j=0; j<q->n; ++j) {

      /* calculate distance -- before the next if! */ 
      d.x = tmp_d.x - *qptr;	++qptr;
      d.y = tmp_d.y - *qptr;	++qptr;
#ifndef TWOD
      d.z = tmp_d.z - *qptr;	++qptr;
#endif

      /* don't compute 'selfinteraction' */
      if (!((q==p) && (i==j))) {

        q_typ = SORTE(q,j);
        r2    = SPROD(d,d);

        if ((r2<rho_h_tab.end[p_typ*ntypes+q_typ]) || 
            (r2<rho_h_tab.end[q_typ*ntypes+p_typ])) {

          /* f_j_strich(rho_h_j) -- TEST CUTOFF??? */
          DERIV_FUNC(f_j_strich, embed_pot, q_typ, 
                     ntypes, q->eam2_rho_h[j], idummy);

          /* Take care: particle i gets its rho from particle j.
             This is tabulated in column p_typ*ntypes+q_typ.
             Here we need the giving part from column q_typ*ntypes+p_typ.
           */

          /* rho_strich_i(r_ij) */
          DERIV_FUNC(rho_i_strich, rho_h_tab, q_typ*ntypes+p_typ, 
                     ntypes*ntypes, r2, is_short);

          /* rho_strich_j(r_ij) */
          DERIV_FUNC(rho_j_strich, rho_h_tab, p_typ*ntypes+q_typ, 
                     ntypes*ntypes, r2, is_short);

          /* put together (f_i_strich and f_j_strich are by 0.5 too big) */
	  eam2_force = -0.5*(f_i_strich*rho_j_strich+f_j_strich*rho_i_strich);

          p->kraft X(i) += d.x  * eam2_force;
          p->kraft Y(i) += d.y  * eam2_force;
#ifndef TWOD
          p->kraft Z(i) += d.z  * eam2_force;
#endif

#ifdef DEBUG_INFO
          if (p->nummer[i]==NR){
            printf("****TTTT p_typ: %d, q_typ: %d\n",p_typ,q_typ);
            printf("**** r=%.12lf thisrho:%.14lf other_rho:%.14lf\n",
                   r2, p->eam2_rho_h[i], q->eam2_rho_h[j]);
            printf("**** f_i`: %.14lf rho_j`: %.14lf  f_j`: %.12lf rho_i`: %.14lf\n",
	           f_i_strich, rho_j_strich, f_j_strich, rho_i_strich);
            printf("**** eam2_force: %.14lf f.x:%.16lf f.y:%.16lf f.z:%.16lf \n",
		   eam2_force, d.x*eam2_force, d.y*eam2_force, d.z*eam2_force);
            printf("**** d.x:%lf d.y:%lf d.z:%lf\n",d.x,d.y,d.z);
            printf("**** Total Force :f.x: %.16lf f.y:%.16lf f.z:%.16lf\n",
		   p->kraft X(i),p->kraft Y(i),p->kraft Z(i));
	  }
#endif

          /* virial and pressure tensor are updated twice for each pair */
          eam2_force *= 0.5;
#ifdef P_AXIAL
          tmp_vir_vect.x -= d.x * d.x * eam2_force;
          tmp_vir_vect.y -= d.y * d.y * eam2_force;
#ifndef TWOD
          tmp_vir_vect.z -= d.z * d.z * eam2_force;
#endif
#else
          tmp_virial     -= r2 * eam2_force;  
#endif

#ifdef STRESS_TENS
          p->presstens X(i) -= d.x * d.x * eam2_force;
          p->presstens Y(i) -= d.y * d.y * eam2_force;
#ifdef TWOD
          p->presstens_offdia[i] -= d.x * d.y * eam2_force;
#else
          p->presstens Z(i) -= d.z * d.z * eam2_force;
          p->presstens_offdia X(i) -= d.y * d.z * eam2_force;
          p->presstens_offdia  Y(i) -= d.z * d.x * eam2_force;
          p->presstens_offdia Z(i) -= d.x * d.y * eam2_force;
#endif
#endif
	} /* if in the cutoff range */
      } /* if not the same particle */    
    } /* for j */
  } /* for i */

  /* print warning if short distance occurred */
  if (is_short==1) fprintf(stderr, "\n Short distance!\n");

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

} /* do_forces_eam2 */
