
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

/*****************************************************************************
*
* imd_forces -- force loop
* this works for uniaxial molecules!
*
******************************************************************************/

#include "imd.h"

/* do_forces_uniax */

void do_forces(cell *p, cell *q, vektor pbc, real *Epot, real *Virial, 
               real *Vir_xx, real *Vir_yy, real *Vir_zz,
               real *Vir_yz, real *Vir_zx, real *Vir_xy)
{
  int i, j ;
  int jstart ;

  real tmp_virial ;
  vektor tmp_vir_vect ; 
  vektor tmp_r12 ;

  vektor r12 ;
  vektor e1 ;
  vektor e2 ;
  real rsqr ;
  vektor s1 ;
  vektor w1 ;

  real pot12 ;
  vektor force12 ;
  vektor torque12 ;
  vektor torque21 ;

  /* actual pair virial and virial components */
  
  tmp_virial     = 0.0;
  tmp_vir_vect.x = 0.0;
  tmp_vir_vect.y = 0.0;
  tmp_vir_vect.z = 0.0;
    
  /* For each atom in first cell */
  for (i = 0;i < p->n; ++i) {
    /* For each atom in neighbouring cell */
    /* Some compilers don't find the expressions that are invariant 
       to the inner loop. I'll have to define my own temp variables. */

    tmp_r12.x = p->ort X(i) - pbc.x;
    tmp_r12.y = p->ort Y(i) - pbc.y;
    tmp_r12.z = p->ort Z(i) - pbc.z;

    e1.x = p->achse X(i) ;
    e1.y = p->achse Y(i) ;
    e1.z = p->achse Z(i) ;

    s1.x = p->shape X(i) ;
    s1.y = p->shape Y(i) ;
    s1.z = p->shape Z(i) ;

    w1.x = p->pot_well X(i) ;
    w1.y = p->pot_well Y(i) ;
    w1.z = p->pot_well Z(i) ;

#ifdef TWOD
    jstart = (((p==q) && (pbc.x==0) && (pbc.y==0))               ? i+1 : 0);
#else
    jstart = (((p==q) && (pbc.x==0) && (pbc.y==0) && (pbc.z==0)) ? i+1 : 0);
#endif
    
    for (j = jstart; j < q->n; ++j) {
      
      /* Calculate distance */

      r12.x = tmp_r12.x - q->ort X(j) ;
      r12.y = tmp_r12.y - q->ort Y(j) ;
      r12.z = tmp_r12.z - q->ort Z(j) ;

      rsqr = SPROD(r12,r12);

      e2.x = q->achse X(j) ;
      e2.y = q->achse Y(j) ;
      e2.z = q->achse Z(j) ;

#ifndef NODBG_DIST
      if (0==rsqr) 
	{ 
	  char msgbuf[256];
	  sprintf(msgbuf,"Distance is zero: i=%d, j=%d\n",i,j);
	  error(msgbuf);
	}
#else
      if (0==rsqr) error("Distance is zero.");
#endif

      if (rsqr <= uniax_r2_cut) {

	/* calculate interactions, if distance smaller than cutoff radius */ 

	gay_berne( r12, e1, e2, rsqr, s1, w1, &pot12,
		   &force12, &torque12, &torque21) ;
	
        /* accumulate forces */

	p->kraft X(i) += force12.x;
	p->kraft Y(i) += force12.y;
	p->kraft Z(i) += force12.z;

	q->kraft X(j) -= force12.x;
	q->kraft Y(j) -= force12.y;
	q->kraft Z(j) -= force12.z;

        /* accumulate torques */

	p->dreh_moment X(i) += torque12.x;
	p->dreh_moment Y(i) += torque12.y;
	p->dreh_moment Z(i) += torque12.z;

	q->dreh_moment X(j) += torque21.x;
	q->dreh_moment Y(j) += torque21.y;
	q->dreh_moment Z(j) += torque21.z;

	p->pot_eng[i] += pot12;
	q->pot_eng[j] += pot12;
        *Epot         += pot12;

        tmp_vir_vect.x += r12.x * force12.x ;
        tmp_vir_vect.y += r12.y * force12.y ;
        tmp_vir_vect.z += r12.z * force12.z ;
	tmp_virial += r12.x * force12.x
	  + r12.y * force12.y + r12.z * force12.z ;

      }; /* if */

    }; /* for j */

  }; /* for i */

  *Vir_xx += tmp_vir_vect.x;
  *Vir_yy += tmp_vir_vect.y;
  *Vir_zz += tmp_vir_vect.z;
  *Virial += tmp_virial ;

} /* do_forces_uniax */

