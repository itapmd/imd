/*****************************************************************************
*
* imd_forces_eam2_risc -- force loops
* uses some parts of do_forces, so please keep up to date
* version that uses tables in r2
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

/* Personal Debug Switch 
#define DEBUG_INFO 1000 
#define NR 94
*/

/*#define ATOMNR 100  quick hack to print the position of a special atom
		      better (later) Nr. in parameterfile... */ 

#include "imd.h"
/* #include "makros.h" */
#include <math.h>

/* not supported yet by eam2*********
   #if defined(PVPCRAY) || defined(SX4) 
   #include "imd_forces_vec.c" 
   #else 
*************************************/

/******************************************************************************
*
*  eam2_do_forces1, version for scalar processors
*
*  computes the forces between atoms in two given cells
*
*  first loop, calculates forces and energy for the 
*  core-core potential, and the embedding electron density at the atoms' sites
*
*    doesn't use 'actio = reactio' !!!!
******************************************************************************/



void eam2_do_forces1(cell *p, cell *q, vektor pbc)     
{
  int i,j,k;
  vektor d;
  vektor tmp_d;
  vektor force;
  real radius2;
  vektor tmp_forces;
  real tmp_energy;

#ifdef P_AXIAL
  vektor tmp_vir_vect;
#endif
  real pot_zwi;
  real pot_grad;
  int jstart,jend;
  int c;
  real pot_k0,pot_k1,pot_k2;
  int q_typ,p_typ;

  real *qptr;  

  real tmp_virial     = 0.0;
#ifdef P_AXIAL
  tmp_vir_vect.x = 0.0;
  tmp_vir_vect.y = 0.0;
#ifndef TWOD
  tmp_vir_vect.z = 0.0;
#endif
#endif
  
  /* EAM2 variables */
  real eam2_r;
  real eam2_x;
  real eam2_dr,eam2_dr_inv, eam2_r_cut,eam2_r0;
  int  eam2_k, eam2_nsteps;

  /* makes the interpolation readable */
  real eam2_interpol_fac1=0.0,eam2_interpol_fac2=0.0;
  real eam2_interpol_fac3=0.0,eam2_interpol_fac4=0.0;
  real *eam2_rho_at_potptr;
  real *eam2_phi_potptr;
  real eam2_rho=0.0;
  real eam2_rho_at_k1=0.0,eam2_rho_at_k2=0.0;
  real eam2_rho_at_k3=0.0,eam2_rho_at_k4=0.0;
  real eam2_phi_k1=0.0,eam2_phi_k2=0.0;
  real eam2_phi_k3=0.0,eam2_phi_k4=0.0;
  real eam2_dphi_interpol_fac1=0.0,eam2_dphi_interpol_fac2=0.0;
  real eam2_dphi_interpol_fac3=0.0,eam2_dphi_interpol_fac4=0.0; 


  /* For each atom in first cell */
  for (i = 0;i < p->n; ++i) {

    /* For each atom in neighbouring cell */
    /* Some compilers don't find the expressions that are invariant 
       to the inner loop. I'll have to define my own temp variables. */

    tmp_d.x = p->ort X(i) - pbc.x;
    tmp_d.y = p->ort Y(i) - pbc.y;
#ifndef TWOD
    tmp_d.z = p->ort Z(i) - pbc.z;
#endif
    p_typ   = SORTE(p,i);
   
    jstart =0;                             /* we go over all particles */
    qptr   = q->ort + DIM * jstart;
    
    for (j = jstart; j < q->n; ++j) {
       q_typ = SORTE(q,j);
      /* Calculate distance  */
      d.x = tmp_d.x - *qptr;
      ++qptr;
      d.y = tmp_d.y - *qptr;
      ++qptr;
#ifndef TWOD
      d.z = tmp_d.z - *qptr;
      ++qptr;
#endif

      if(!(p==q && i==j)){         /* don't compute self interaction !*/
	
	radius2 = SPROD(d,d);

	/* check if Distance == 0 */
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


	/* eam2_r  = sqrt(radius2);                  
	   this is history, but i didn't wanted to change all 
	   variable names to ...r2, sorry, 
	   maybe i'll find the time to clean up...
	 */  
	eam2_r = radius2;

	eam2_r_cut=(*PTR_2D(eam2_phi_r_end,p_typ,q_typ,ntypes,ntypes));
	eam2_dr = *PTR_2D(eam2_phi_r_step,p_typ,q_typ,ntypes,ntypes); 
	if (eam2_r < eam2_r_cut-(2.0*eam2_dr)) {                              

	  /* interpolation: 4 Point Lagrange */
	  
	  /* get information about the table, depending on the atomtype */
	  eam2_r0 = *PTR_2D(eam2_phi_r_begin,p_typ,q_typ,ntypes,ntypes);
	  eam2_dr_inv = 1.0/eam2_dr;
	  eam2_nsteps = (int)( (eam2_r_cut - eam2_r0)*eam2_dr_inv +0.5 ); /*rounding?*/
      
	  /* one quick hack to treat the borders of table */
	  eam2_x= (eam2_r - eam2_r0)*eam2_dr_inv;
	  eam2_k= (int) (eam2_x);                           
	  /* eam2_k= MIN(eam2_k, eam2_nsteps-2); eliminated because if (eam2_r < .. */ 
	  eam2_k= MAX(eam2_k,1);
	  eam2_x= eam2_x-eam2_k;
	  eam2_x= MIN(eam2_x,2.0);

	  /* factors for the interpolation */
	  eam2_interpol_fac1 = (-1.0/6.0)*eam2_x*(eam2_x-1.0)*(eam2_x-2.0);
	  eam2_interpol_fac2 = 0.5*(eam2_x*eam2_x-1.0)*(eam2_x-2.0);
	  eam2_interpol_fac3 = (-0.5)*eam2_x*(eam2_x+1.0)*(eam2_x-2.0);
	  eam2_interpol_fac4 = (1.0/6.0)*eam2_x*(eam2_x*eam2_x-1.0);

	  /* factors for the interpolation of the 1. derivative */
	  eam2_dphi_interpol_fac1 = (-1.0/6.0) *((3*eam2_x-6.0)*eam2_x+2.0);           
	  eam2_dphi_interpol_fac2 = 0.5*((3.0*eam2_x-4.0)*eam2_x-1.0);                  
	  eam2_dphi_interpol_fac3 = (-0.5)*((3.0*eam2_x-2.0)*eam2_x-2.0);
	  eam2_dphi_interpol_fac4 = (1.0/6.0)*(3.0*eam2_x*eam2_x-1.0);
	  
	  /* intermediate pointers */
	  eam2_phi_potptr=PTR_3D(eam2_phi,eam2_k-1,p_typ,q_typ,eam2_max_phi_r_steps,ntypes,ntypes);
	  eam2_phi_k1 = *eam2_phi_potptr; eam2_phi_potptr += ntypes * ntypes;
	  eam2_phi_k2 = *eam2_phi_potptr; eam2_phi_potptr += ntypes * ntypes;
	  eam2_phi_k3 = *eam2_phi_potptr; eam2_phi_potptr += ntypes * ntypes;
	  eam2_phi_k4 = *eam2_phi_potptr;

	  /* interpolation of potential energy, Factor 0.5 depending on the table */ 
	  pot_zwi = 0.5* ( eam2_interpol_fac1 * eam2_phi_k1 +
	                   eam2_interpol_fac2 * eam2_phi_k2 +
			   eam2_interpol_fac3 * eam2_phi_k3 +
			   eam2_interpol_fac4 * eam2_phi_k4  );
	  
	  /* 1. derivative of interpolationfunction gives us the derivative of the
	     potential energy. the table is in r2, so (see Allen Tildesley p.145)
	     with   w = r * d/dr V(r)
	            F = - \vec r /r * d/dr V(r) = - w/r**2 * \vec r
	     we want to calculate with tables in s=r**2:
	     w(s)/s =2* d/ds V(s) */

	  pot_grad = -  eam2_dr_inv*2.0*( eam2_dphi_interpol_fac1 *eam2_phi_k1 +          
			                  eam2_dphi_interpol_fac2 *eam2_phi_k2 +
				          eam2_dphi_interpol_fac3 *eam2_phi_k3 +
				          eam2_dphi_interpol_fac4 *eam2_phi_k4  );
	  

	  /* Store forces in temp */
	  force.x = d.x * pot_grad;
	  force.y = d.y * pot_grad;
#ifndef TWOD
	  force.z = d.z * pot_grad;
#endif
	  
        /* Accumulate forces due to core-core potential */
	  p->kraft X(i) += force.x;
	  p->kraft Y(i) += force.y;
#ifndef TWOD
	  p->kraft Z(i) += force.z;
#endif

	  p->pot_eng[i]  += pot_zwi; 
	  tot_pot_energy += pot_zwi;

#ifdef DEBUG_INFO
	  if (p->nummer[i]==NR){
	    printf("d.x: %lf d.y: %lf d.z: %lf -> r: %lf\n",d.x,d.y,d.z,eam2_r);
	    printf("-Phi'(r): %lf -> force: %.16lf %.16lf %.16lf\n",
		   eam2_r* pot_grad,force.x, force.y,force.z);
	    printf("Ges. Kraft: f.x: %.16lf f.y:%.16lf f.z:%.16lf\n",
		   p->kraft X(i),p->kraft Y(i),p->kraft Z(i));
	    printf("p_typ: %d, q_typ: %d  p->pot_eng[i]: %.12lf\n",p_typ,q_typ,p->pot_eng[i]);
	  }
#endif
	  /* pure analogy programming */
	  
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
#ifdef TWOD
	  p->presstens_offdia[i] -= d.x * d.y * pot_grad;
#else
	  p->presstens Z(i) -= d.z * d.z * pot_grad;
	  p->presstens_offdia X(i) -= d.y * d.z * pot_grad;
	  p->presstens_offdia Y(i) -= d.z * d.x * pot_grad;
	  p->presstens_offdia Z(i) -= d.x * d.y * pot_grad;
#endif
#endif
	}; /* if radius2 <= r2_cut */

	eam2_r_cut = *PTR_2D(eam2_r_end,p_typ,q_typ,ntypes,ntypes);
	eam2_dr = *PTR_2D(eam2_r_step,p_typ,q_typ,ntypes,ntypes);
	if(eam2_r < eam2_r_cut-(2.0*eam2_dr))
	  {
	  /* here comes the EAM2  related stuff:
	   * give each atom the electron density at the atoms position 
	   * interpolation: 4 Point Lagrange */

	    /* get informatio nabaout the table, depending on type */
	    eam2_r0    = *PTR_2D(eam2_r_begin,p_typ,q_typ,ntypes,ntypes);
	    eam2_dr_inv =1.0/eam2_dr;
	    eam2_nsteps = (int)((eam2_r_cut - eam2_r0 )*eam2_dr_inv +0.5);/*rounding?*/

	    /* one quick hack to treat the boarders of table */
	    eam2_x= (eam2_r - eam2_r0)*eam2_dr_inv;
	    eam2_k= (int) (eam2_x);                           
	    /* eam2_k= MIN(eam2_k, eam2_nsteps-2);*/                           
	    eam2_k= MAX(eam2_k,1);
	    eam2_x= eam2_x-eam2_k;
	    eam2_x= MIN(eam2_x,2.0);
	
	    /* here comes the interpolation...*/
	    eam2_interpol_fac1 = (-1./6.0) *eam2_x*(eam2_x-1.0)*(eam2_x-2.0);
	    eam2_interpol_fac2 = 0.5*(eam2_x*eam2_x-1.0)*(eam2_x-2.0);
	    eam2_interpol_fac3 = (-0.5)*eam2_x*(eam2_x+1.0)*(eam2_x-2.0);
	    eam2_interpol_fac4 = (1.0/6.0)*eam2_x*(eam2_x*eam2_x-1.0);

	    /* version with intermediate pointers */
	    eam2_rho_at_potptr = PTR_3D(eam2_rho_at, eam2_k-1, p_typ, q_typ, 
					eam2_max_r_steps, ntypes, ntypes);
	    eam2_rho_at_k1 = *eam2_rho_at_potptr; eam2_rho_at_potptr += ntypes * ntypes;
	    eam2_rho_at_k2 = *eam2_rho_at_potptr; eam2_rho_at_potptr += ntypes * ntypes;
	    eam2_rho_at_k3 = *eam2_rho_at_potptr; eam2_rho_at_potptr += ntypes * ntypes;
	    eam2_rho_at_k4 = *eam2_rho_at_potptr;

	    eam2_rho = eam2_interpol_fac1 * eam2_rho_at_k1 +
	               eam2_interpol_fac2 * eam2_rho_at_k2 +
                       eam2_interpol_fac3 * eam2_rho_at_k3 +
	               eam2_interpol_fac4 * eam2_rho_at_k4 ;

	    p->eam2_rho_h[i] +=eam2_rho; 

	  } /* if(eam2_r<eam2_r_cut-2*eam_dr) */
      } /* if ! q==p.. */

#ifdef DEBUG_INFO
      else{
	if (p->nummer[i]==NR)
	  printf("++++++ Ort von Teilchen %d : x: %lf y: %lf z:%lf \n",
		 p->nummer[i],tmp_d.x,tmp_d.y,tmp_d.z);
      }
#endif

    }; /* for j */
  }; /* for i */
  
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


} /* eam2_do_forces1 */



/* second force loop, calculates the force and the energy 
   caused by the embedding electron density
   uses Phi(r2), Rho(r2), F(rho) and it's derivatives
*/

void eam2_do_forces2(cell *p, cell *q, vektor pbc)
{
  int i,j,k;
  vektor d;
  vektor tmp_d;
  vektor force;
  real radius2;
  vektor tmp_forces;
  real tmp_energy;
  real tmp_virial=0.0;
  int jstart,jend;
  real r2_short = cellsz;
  int q_typ,p_typ;
  real *qptr; 
#ifdef P_AXIAL
  vektor tmp_vir_vect;
#endif
  /* EAM2 variables */
  real eam2_r;
  real eam2_energy;
  real eam2_x,eam2_rho2_x;
  int  eam2_k,eam2_rho2_k;
  real eam2_this_rho, eam2_other_rho;
  real eam2_r0=0.0, eam2_rho0=0.0;
  real eam2_rho_x, eam2_drho,eam2_drho_inv;
  real eam2_rho20=0.0;
  real eam2_drho2,eam2_drho2_inv;
  int  eam2_rho_k, eam2_rho_nsteps,eam2_rho2_nsteps;
  real eam2_dr_pq,eam2_dr_qp,eam2_dr_inv;
  int  eam2_r_nsteps;
  real rho_r_cut_pq,rho_r_cut_qp;
  /* makes the interpolation readable, don't really need different variables! */
  real eam2_f_interpol_fac1=0.0,eam2_f_interpol_fac2=0.0;
  real eam2_f_interpol_fac3=0.0,eam2_f_interpol_fac4=0.0;
  real eam2_df_interpol_fac1=0.0,eam2_df_interpol_fac2=0.0;
  real eam2_df_interpol_fac3=0.0,eam2_df_interpol_fac4=0.0;
  real eam2_df2_interpol_fac1=0.0,eam2_df2_interpol_fac2=0.0;
  real eam2_df2_interpol_fac3=0.0,eam2_df2_interpol_fac4=0.0;
  real eam2_drho_at_interpol_fac1=0.0,eam2_drho_at_interpol_fac2=0.0;
  real eam2_drho_at_interpol_fac3=0.0,eam2_drho_at_interpol_fac4=0.0;
  real *eam2_f_i_potptr, *eam2_rho_at_potptr;
  real eam2_f_i_k1=0.0,eam2_f_i_k2=0.0,eam2_f_i_k3=0.0,eam2_f_i_k4=0.0;
  real eam2_rho_at_k1=0.,eam2_rho_at_k2=0.,eam2_rho_at_k3=0.,eam2_rho_at_k4=0.;

  real eam2_f_i_strich,eam2_f_j_strich,eam2_rho_i_strich,eam2_rho_j_strich;
  real eam2_force;

#ifdef P_AXIAL
  tmp_vir_vect.x = 0.0;
  tmp_vir_vect.y = 0.0;
#ifndef TWOD
  tmp_vir_vect.z = 0.0;
#endif
#endif

  /* For each atom in first cell */
  for (i = 0;i < p->n; ++i) {
    tmp_d.x = p->ort X(i) - pbc.x;
    tmp_d.y = p->ort Y(i) - pbc.y;
#ifndef TWOD
    tmp_d.z = p->ort Z(i) - pbc.z;
#endif
    p_typ   = SORTE(p,i);

    /*go over ALL (except q==p && i==j) atoms, no "Actio = Reactio" */
    jstart = 0;                     
    qptr   = q->ort + DIM * jstart;
    eam2_this_rho = p->eam2_rho_h[i];

   /* get information about the Embedding energy function table (eam2_f_i)*/ 
    eam2_drho       = *(eam2_rho_step+p_typ); 
    eam2_drho_inv   = 1.0/eam2_drho;
    eam2_rho0       = *(eam2_rho_begin+p_typ);
    eam2_rho_nsteps = (int)((*(eam2_rho_end+p_typ)-eam2_rho0)*eam2_drho_inv +.5);


    /* f_i_strich(rho_h_i)  *******************************************************/

    eam2_rho_x= (eam2_this_rho - eam2_rho0)*eam2_drho_inv;                      
    eam2_rho_k= (int) (eam2_rho_x);                           
    eam2_rho_k= MIN(eam2_rho_k, eam2_rho_nsteps - 2);                           
    eam2_rho_k= MAX(eam2_rho_k,1);
    eam2_rho_x=eam2_rho_x - eam2_rho_k;
    eam2_rho_x= MIN(eam2_rho_x,2.0);
    
    eam2_f_i_potptr = PTR_2D(eam2_f_i, eam2_rho_k-1, p_typ, eam2_max_rho_steps, ntypes);
    eam2_f_i_k1 = *eam2_f_i_potptr; eam2_f_i_potptr += ntypes;
    eam2_f_i_k2 = *eam2_f_i_potptr; eam2_f_i_potptr += ntypes;
    eam2_f_i_k3 = *eam2_f_i_potptr; eam2_f_i_potptr += ntypes;
    eam2_f_i_k4 = *eam2_f_i_potptr;

    /* factors for the interpolation of the 1. derivative of f_i(this_rho) */
    eam2_df_interpol_fac1 = (-1.0/6.0) *((3.0*eam2_rho_x-6.0)*eam2_rho_x+2.0);           
    eam2_df_interpol_fac2 = 0.5*((3.0*eam2_rho_x-4.0)*eam2_rho_x-1.0);                
    eam2_df_interpol_fac3 = (-0.5)*((3.0*eam2_rho_x-2.0)*eam2_rho_x-2.0);
    eam2_df_interpol_fac4 = (1.0/6.0)*(3.0*eam2_rho_x*eam2_rho_x-1.0);
    
    
    eam2_f_i_strich = eam2_drho_inv*( eam2_df_interpol_fac1 *eam2_f_i_k1 +     
				      eam2_df_interpol_fac2 *eam2_f_i_k2 +
				      eam2_df_interpol_fac3 *eam2_f_i_k3 +
				      eam2_df_interpol_fac4 *eam2_f_i_k4 );

	   
    for (j = jstart; j < q->n; ++j) 
      {
	/* Calculate distance, this has to happen before the if, else *qptr isn't actual!*/
	d.x = tmp_d.x - *qptr;
	++qptr;
	d.y = tmp_d.y - *qptr;
	++qptr;
#ifndef TWOD
	d.z = tmp_d.z - *qptr;
	++qptr;
#endif
	if((q==p && i==j)) 
	  {
	    /* here comes the energy, energy is calculated only once per particle! */ 
	    
	    /* get table info */
	    eam2_drho       = *(eam2_rho_step+p_typ); 
	    eam2_drho_inv   = 1.0/eam2_drho;
	    eam2_rho0       = *(eam2_rho_begin+p_typ);
	    eam2_rho_nsteps = (int)((*(eam2_rho_end+p_typ)-eam2_rho0)*eam2_drho_inv + .5);

	    
	    /* handle boarders */
	    eam2_rho_x= (eam2_this_rho - eam2_rho0)*eam2_drho_inv;
	    eam2_rho_k= (int) (eam2_rho_x);                           
	    eam2_rho_k= MIN(eam2_rho_k, eam2_rho_nsteps - 2);          
	    eam2_rho_k= MAX(eam2_rho_k,1);
	    eam2_rho_x= eam2_rho_x - eam2_rho_k;
	    eam2_rho_x= MIN(eam2_rho_x,2.0);
	    

	    eam2_f_i_potptr = PTR_2D(eam2_f_i, eam2_rho_k-1, p_typ, eam2_max_rho_steps, ntypes); 
	    eam2_f_i_k1 = *eam2_f_i_potptr; eam2_f_i_potptr += ntypes;   
	    eam2_f_i_k2 = *eam2_f_i_potptr; eam2_f_i_potptr += ntypes;  
	    eam2_f_i_k3 = *eam2_f_i_potptr; eam2_f_i_potptr += ntypes;  
	    eam2_f_i_k4 = *eam2_f_i_potptr; 
		
	    eam2_f_interpol_fac1 = (-1.0/6.0) *eam2_rho_x*(eam2_rho_x-1.0)*(eam2_rho_x-2.0);
	    eam2_f_interpol_fac2 = 0.5*(eam2_rho_x*eam2_rho_x-1.0)*(eam2_rho_x-2.0);
	    eam2_f_interpol_fac3 = (-0.5)*eam2_rho_x*(eam2_rho_x+1.0)*(eam2_rho_x-2.0);
	    eam2_f_interpol_fac4 = (1.0/6.0)*eam2_rho_x*(eam2_rho_x*eam2_rho_x-1.0);
	    
	    eam2_energy = eam2_f_interpol_fac1 * eam2_f_i_k1 +
	                  eam2_f_interpol_fac2 * eam2_f_i_k2 +
	                  eam2_f_interpol_fac3 * eam2_f_i_k3 +
	                  eam2_f_interpol_fac4 * eam2_f_i_k4 ;
	    
	    p->pot_eng[i]  += eam2_energy;
	    tot_pot_energy += eam2_energy;

#ifdef ATOMNR
	   if (p->nummer[i]==ATOMNR){
	     printf("%.8lf %.8lf %.8lf %lf %lf %lf\n", tmp_d.x, tmp_d.y, tmp_d.z,(p->impuls X(i))/p->masse[i],(p->impuls Y(i))/p->masse[i],(p->impuls Z(i))/p->masse[i]);
	    	  }
#endif 
	    
#ifdef DEBUG_INFO
	  if (p->nummer[i]==NR){
	    printf("!!!!!!!! rho_h: %.12lf eam2_energy: %.12lf\t",eam2_this_rho,eam2_energy);
	    printf("p_typ: %d, p->pot_eng[i]: %.12lf\n",p_typ,p->pot_eng[i]);
	  }
#endif
	  }
	else{                       /* don't compute 'selfinteraction' */
	  q_typ = SORTE(q,j);
	  radius2 = SPROD(d,d);
	
	  /* eam2_r  = sqrt(radius2); sorry, i didn't changed the variable names to ...r2 */ 
	  eam2_r = radius2;

	  eam2_other_rho = q->eam2_rho_h[j];

	  rho_r_cut_pq = *PTR_2D(eam2_r_end,p_typ,q_typ,ntypes,ntypes);
	  rho_r_cut_qp = *PTR_2D(eam2_r_end,q_typ,p_typ,ntypes,ntypes); 
	  eam2_dr_pq   = *PTR_2D(eam2_r_step,p_typ,q_typ,ntypes,ntypes);
	  eam2_dr_qp   = *PTR_2D(eam2_r_step,q_typ,p_typ,ntypes,ntypes);
	  
	  if(eam2_r<rho_r_cut_pq-(2.0*eam2_dr_pq) || eam2_r<rho_r_cut_qp-(2.0*eam2_dr_qp))
	    {
	      
	      /* now we go for the forces... */
	      /* use the usual 4 Point Lagrange interpolation */ 
	      
	      /* f_j_strich(rho_h_j)  *******************************************************/
	      
	      /* get information about the Embedding energy function table (eam2_f_i)*/ 
	      eam2_drho2       = *(eam2_rho_step+q_typ); 
	      eam2_drho2_inv   = 1.0/eam2_drho2;
	      eam2_rho20       = *(eam2_rho_begin+q_typ);
	      eam2_rho2_nsteps = (int)((*(eam2_rho_end+q_typ)-eam2_rho20)*eam2_drho2_inv + .5);

	      
	      /* treat the boarders of table */
	      eam2_rho2_x= (eam2_other_rho - eam2_rho20)*eam2_drho2_inv;
	      eam2_rho2_k= (int) (eam2_rho2_x);                           
	      eam2_rho2_k= MIN(eam2_rho2_k, eam2_rho2_nsteps - 2);                           
	      eam2_rho2_k= MAX(eam2_rho2_k,1);
	      eam2_rho2_x= eam2_rho2_x - eam2_rho2_k;
	      eam2_rho2_x= MIN(eam2_rho2_x,2.0);
	  
	      eam2_f_i_potptr = PTR_2D(eam2_f_i, eam2_rho2_k-1, q_typ, eam2_max_rho_steps, ntypes);
	      eam2_f_i_k1 = *eam2_f_i_potptr; eam2_f_i_potptr += ntypes;
	      eam2_f_i_k2 = *eam2_f_i_potptr; eam2_f_i_potptr += ntypes;
	      eam2_f_i_k3 = *eam2_f_i_potptr; eam2_f_i_potptr += ntypes;
	      eam2_f_i_k4 = *eam2_f_i_potptr;

	      /* factors for the interpolation of the 1. derivative of f_j(other_rho) */
	      eam2_df2_interpol_fac1 = (-1.0/6.0) *((3.0*eam2_rho2_x-6.0)*eam2_rho2_x+2.0);           
	      eam2_df2_interpol_fac2 = 0.5*((3.0*eam2_rho2_x-4.0)*eam2_rho2_x-1.0);                   
	      eam2_df2_interpol_fac3 = (-0.5)*((3.0*eam2_rho2_x-2.0)*eam2_rho2_x-2.0);
	      eam2_df2_interpol_fac4 = (1.0/6.0)*(3.0*eam2_rho2_x*eam2_rho2_x-1.0);
	      
	      eam2_f_j_strich=eam2_drho_inv*( eam2_df2_interpol_fac1 *eam2_f_i_k1 +          
					      eam2_df2_interpol_fac2 *eam2_f_i_k2 +
					      eam2_df2_interpol_fac3 *eam2_f_i_k3 +
					      eam2_df2_interpol_fac4 *eam2_f_i_k4 );

		    
	      /* rho_at_strich_i(r_ij) ****************************************************/
	      /* is of course tabulated in r**2, this gives a factor 2/r, see first force loop
		 for details
		 take care: particle i gets its rho from particle j:
	                      PTR_3D(eam2_rho_at,eam2_k-1,p_typ,q_typ,eam2_max_r_steps,ntypes,ntypes);
			      means in the table: 0 0   NI Rho file
			                          0 1 -> 1 gives his charge to 0 so: AL Rho file
                                                  1 0 -> 0 gives to 1 so :           NI Rho file
						  1 1   AL Rho file
			    BUT here i needs to have the his proper Rho-function evaluated! 
			    that means the giving part of p_type:
			      PTR_3D(eam2_rho_at,eam2_k-1,p_typ,q_typ,eam2_max_r_steps,ntypes,ntypes);
	      */

	      /* get information about the atomic electron density function table (eam2_rho_at)*/ 
	     
	      eam2_dr_inv   = 1.0/eam2_dr_qp;
	      eam2_r0       = *PTR_2D(eam2_r_begin,q_typ,p_typ,ntypes,ntypes);
	      eam2_r_nsteps = (int)((rho_r_cut_qp-eam2_r0)*eam2_dr_inv +.5);
	      
	      eam2_x= (eam2_r - eam2_r0)*eam2_dr_inv;
	      eam2_k= (int) (eam2_x);                           
	      /* eam2_k= MIN(eam2_k, eam2_r_nsteps-2);  */     
	      eam2_k= MAX(eam2_k,1);
	      eam2_x= eam2_x-eam2_k;
	      eam2_x= MIN(eam2_x,2.0);
	      
	      /* factors for the interpolation of the 1. derivative of rho_at(r)  */
	      eam2_drho_at_interpol_fac1 = (-1.0/6.0) *((3.0*eam2_x-6.0)*eam2_x+2.0);           
	      eam2_drho_at_interpol_fac2 = 0.5*((3.0*eam2_x-4.0)*eam2_x-1.0);         
	      eam2_drho_at_interpol_fac3 = (-0.5)*((3.0*eam2_x-2.0)*eam2_x-2.0);
	      eam2_drho_at_interpol_fac4 = (1.0/6.0)*(3.0*eam2_x*eam2_x-1.0);

	      eam2_rho_at_potptr = PTR_3D(eam2_rho_at,eam2_k-1,q_typ,p_typ,eam2_max_r_steps,ntypes,ntypes);
	      eam2_rho_at_k1 = *eam2_rho_at_potptr; eam2_rho_at_potptr += ntypes * ntypes;
	      eam2_rho_at_k2 = *eam2_rho_at_potptr; eam2_rho_at_potptr += ntypes * ntypes;
	      eam2_rho_at_k3 = *eam2_rho_at_potptr; eam2_rho_at_potptr += ntypes * ntypes;
	      eam2_rho_at_k4 = *eam2_rho_at_potptr;
	      
	      eam2_rho_i_strich=2.0*eam2_dr_inv*( eam2_drho_at_interpol_fac1 *eam2_rho_at_k1 +     
					          eam2_drho_at_interpol_fac2 *eam2_rho_at_k2 +
					          eam2_drho_at_interpol_fac3 *eam2_rho_at_k3 +
					          eam2_drho_at_interpol_fac4 *eam2_rho_at_k4 );
		
	      
	      /* rho_at_strich_j(r_ij) the same (same r) just with inverted p_type, q_type */ 
	      /* get information about the atomic electron density function table (eam2_rho_at)*/ 
	      eam2_dr_inv   = 1.0/eam2_dr_pq;
	      eam2_r0       = *PTR_2D(eam2_r_begin,p_typ,q_typ,ntypes,ntypes);
	      eam2_r_nsteps = (int)( (rho_r_cut_pq - eam2_r0 )*eam2_dr_inv +.5);

	      
	      eam2_x= (eam2_r - eam2_r0)*eam2_dr_inv;
	      eam2_k= (int) (eam2_x);                           
	      /* eam2_k= MIN(eam2_k, eam2_r_nsteps-2); */       
	      eam2_k= MAX(eam2_k,1);
	      eam2_x= eam2_x-eam2_k;
	      eam2_x= MIN(eam2_x,2.0);
	      
	      /* factors for the interpolation of the 1. derivative of rho_at(r)  */
	      eam2_drho_at_interpol_fac1 = (-1.0/6.0) *((3.0*eam2_x-6.0)*eam2_x+2.0);           
	      eam2_drho_at_interpol_fac2 = 0.5*((3.0*eam2_x-4.0)*eam2_x-1.0);         
	      eam2_drho_at_interpol_fac3 = (-0.5)*((3.0*eam2_x-2.0)*eam2_x-2.0);
	      eam2_drho_at_interpol_fac4 = (1.0/6.0)*(3.0*eam2_x*eam2_x-1.0);

	      eam2_rho_at_potptr = PTR_3D(eam2_rho_at,eam2_k-1,p_typ,q_typ,eam2_max_r_steps,ntypes,ntypes);
	      eam2_rho_at_k1 = *eam2_rho_at_potptr; eam2_rho_at_potptr += ntypes * ntypes;
	      eam2_rho_at_k2 = *eam2_rho_at_potptr; eam2_rho_at_potptr += ntypes * ntypes;
	      eam2_rho_at_k3 = *eam2_rho_at_potptr; eam2_rho_at_potptr += ntypes * ntypes;
	      eam2_rho_at_k4 = *eam2_rho_at_potptr;
		    
	      eam2_rho_j_strich=2.0*eam2_dr_inv*( eam2_drho_at_interpol_fac1 *eam2_rho_at_k1 +     
					      eam2_drho_at_interpol_fac2 *eam2_rho_at_k2 +
					      eam2_drho_at_interpol_fac3 *eam2_rho_at_k3 +
					      eam2_drho_at_interpol_fac4 *eam2_rho_at_k4 );

		    
	      /* put together ************************************************************/ 
	      eam2_force = -(eam2_f_i_strich*eam2_rho_j_strich + eam2_f_j_strich*eam2_rho_i_strich);
	    
	      p->kraft X(i) += d.x  * eam2_force;
	      p->kraft Y(i) += d.y  * eam2_force;
#ifndef TWOD
	      p->kraft Z(i) += d.z  * eam2_force;
#endif

#ifdef DEBUG_INFO
	      if (p->nummer[i]==NR){
		printf("****TTTT p_typ: %d, q_typ: %d\n",p_typ,q_typ);
		printf("**** r=%.12lf thisrho:%.14lf other_rho:%.14lf\n",
		       eam2_r,eam2_this_rho,eam2_other_rho);
		printf("**** f_i`: %.14lf rho_j`: %.14lf  f_j`: %.12lf rho_i`: %.14lf\n",
		       eam2_f_i_strich,eam2_rho_j_strich,eam2_f_j_strich,eam2_rho_i_strich);
		printf("**** eam2_force: %.14lf f.x:%.16lf f.y:%.16lf f.z:%.16lf \n"
		       ,eam2_force, d.x*eam2_force, d.y*eam2_force, d.z*eam2_force);
		printf("**** d.x:%lf d.y:%lf d.z:%lf\n",d.x,d.y,d.z);
		printf("**** Total Force :f.x: %.16lf f.y:%.16lf f.z:%.16lf\n",
		       p->kraft X(i),p->kraft Y(i),p->kraft Z(i));
				    }
#endif

		    /* Anlogy-programming of the following ... */
		   		    
#ifdef P_AXIAL
		    tmp_vir_vect.x -= d.x * d.x * eam2_force;
		    tmp_vir_vect.y -= d.y * d.y * eam2_force;
#ifndef TWOD
		    tmp_vir_vect.z -= d.z * d.z * eam2_force;
#endif
#else
		    tmp_virial     -= radius2 * eam2_force;  
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

		    
/* TODO: 
 * treatment  of ORDPAR, NVX 
   refer to do_forces, imd_forces_eam,
 */
	    } /* if in the cutoff range */
	} /* if not the same particle */
	    
      }; /* for j */
  }; /* for i */

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

} /* eam2_do_forces2 */









