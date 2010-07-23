
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2010 Institute for Theoretical and Applied Physics,
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

    tmp_r12.x = ORT(p,i,X) - pbc.x;
    tmp_r12.y = ORT(p,i,Y) - pbc.y;
    tmp_r12.z = ORT(p,i,Z) - pbc.z;

    e1.x = ACHSE(p,i,X);
    e1.y = ACHSE(p,i,Y);
    e1.z = ACHSE(p,i,Z);

#ifdef TWOD
    jstart = (((p==q) && (pbc.x==0) && (pbc.y==0))               ? i+1 : 0);
#else
    jstart = (((p==q) && (pbc.x==0) && (pbc.y==0) && (pbc.z==0)) ? i+1 : 0);
#endif
    
    for (j = jstart; j < q->n; ++j) {
      
      /* Calculate distance */

      r12.x = tmp_r12.x - ORT(q,j,X);
      r12.y = tmp_r12.y - ORT(q,j,Y);
      r12.z = tmp_r12.z - ORT(q,j,Z);

      rsqr = SPROD(r12,r12);

      e2.x = ACHSE(q,j,X);
      e2.y = ACHSE(q,j,Y);
      e2.z = ACHSE(q,j,Z);

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

	gay_berne( r12, e1, e2, rsqr, uniax_sig, uniax_eps, &pot12,
		   &force12, &torque12, &torque21) ;
	
        /* accumulate forces */

	KRAFT(p,i,X) += force12.x;
	KRAFT(p,i,Y) += force12.y;
	KRAFT(p,i,Z) += force12.z;

	KRAFT(q,j,X) -= force12.x;
	KRAFT(q,j,Y) -= force12.y;
	KRAFT(q,j,Z) -= force12.z;

        /* accumulate torques */

	DREH_MOMENT(p,i,X) += torque12.x;
	DREH_MOMENT(p,i,Y) += torque12.y;
	DREH_MOMENT(p,i,Z) += torque12.z;

	DREH_MOMENT(q,j,X) += torque21.x;
	DREH_MOMENT(q,j,Y) += torque21.y;
	DREH_MOMENT(q,j,Z) += torque21.z;

        *Epot       += pot12;
        pot12       *= 0.5;   /* avoid double counting */
	POTENG(p,i) += pot12;
	POTENG(q,j) += pot12;

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

