
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_integrate -- various md integrators
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/*#define ATNR 2826  easier debuging... */

/*****************************************************************************
*
* Basic NVE Integrator
*
*****************************************************************************/

#if defined(NVE) || defined(EPITAX)

void move_atoms_nve(void)
{
  int k;
  real tmpvec1[3], tmpvec2[3];
  static int count = 0;
#if defined(NVE) || !defined(EPITAX)
  tot_kin_energy = 0.0;
#endif
  fnorm = 0.0;
  PxF   = 0.0;

  /* loop over all cells */
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tot_kin_energy,fnorm,PxF)
#endif
  for (k=0; k<ncells; ++k) {

    int  i, sort;
    cell *p;
    real kin_energie_1, kin_energie_2, tmp;
#ifdef UNIAX    
    real rot_energie_1, rot_energie_2;
    real dot, norm;
    vektor cross;
#endif
    p = cell_array + CELLS(k);

#ifdef PVPCRAY
#pragma ivdep
#endif
#ifdef SX4
#pragma vdir vector,nodep
#endif
    for (i=0; i<p->n; ++i) {

#if defined(EPITAX) && !defined(NVE) 
      /* beam atoms are always integrated by NVE */
      if ( (NUMMER(p,i)>epitax_sub_n) && (p->pot_eng[i]>(epitax_ctrl * epitax_poteng_min)) ) {
#endif

        kin_energie_1 = SPRODN(p->impuls,i,p->impuls,i);
#ifdef UNIAX
        rot_energie_1 = SPRODN(p->dreh_impuls,i,p->dreh_impuls,i);
#endif

        sort = VSORTE(p,i);
#ifdef FBC
        /* give virtual particles their extra force */
	p->kraft X(i) += (fbc_forces + sort)->x;
	p->kraft Y(i) += (fbc_forces + sort)->y;
#ifndef TWOD
	p->kraft Z(i) += (fbc_forces + sort)->z;
#endif

#endif
	/* and set their force (->momentum) in restricted directions to 0 */
	p->kraft X(i) *= (restrictions + sort)->x;
	p->kraft Y(i) *= (restrictions + sort)->y;
#ifndef TWOD
	p->kraft Z(i) *= (restrictions + sort)->z;
#endif

#ifdef FNORM
	fnorm +=  SPRODN(p->kraft,i,p->kraft,i) / MASSE(p,i);
#endif

	p->impuls X(i) += timestep * p->kraft X(i);
        p->impuls Y(i) += timestep * p->kraft Y(i);
#ifndef TWOD
        p->impuls Z(i) += timestep * p->kraft Z(i);
#endif

	/* "Globale Konvergenz": like mik, just with the global 
           force and momentum vectors */
#ifdef GLOK              
	PxF   +=  SPRODN(p->impuls,i,p->kraft,i);
#endif

#ifdef UNIAX
        dot = 2.0 * SPRODN(p->dreh_impuls,i,p->achse,i);
        p->dreh_impuls X(i) += timestep * p->dreh_moment X(i)
                               - dot * p->achse X(i);
        p->dreh_impuls Y(i) += timestep * p->dreh_moment Y(i)
                               - dot * p->achse Y(i);
        p->dreh_impuls Z(i) += timestep * p->dreh_moment Z(i)
                               - dot * p->achse Z(i);
#endif

        kin_energie_2 = SPRODN(p->impuls,i,p->impuls,i);
#ifdef UNIAX
        rot_energie_2 = SPRODN(p->dreh_impuls,i,p->dreh_impuls,i);
#endif
        /* sum up kinetic energy on this CPU */
        tot_kin_energy += (kin_energie_1 + kin_energie_2) / (4 * MASSE(p,i));
#ifdef UNIAX
        tot_kin_energy += (rot_energie_1 + rot_energie_2) 
                                    / (4 * p->traeg_moment[i]);	  
#endif	  

        /* new positions */
        tmp = timestep / MASSE(p,i);
        p->ort X(i) += tmp * p->impuls X(i);
        p->ort Y(i) += tmp * p->impuls Y(i);
#ifndef TWOD
        p->ort Z(i) += tmp * p->impuls Z(i);
#endif

#ifdef SHOCK
	if (shock_mode == 3) {
	  if (p->ort X(i) > box_x.x) p->impuls X(i) = -p->impuls X(i);
	}
#endif

        /* new molecular axes */
#ifdef UNIAX
        cross.x = p->dreh_impuls Y(i) * p->achse Z(i)
                - p->dreh_impuls Z(i) * p->achse Y(i);
        cross.y = p->dreh_impuls Z(i) * p->achse X(i)
                - p->dreh_impuls X(i) * p->achse Z(i);
        cross.z = p->dreh_impuls X(i) * p->achse Y(i)
                - p->dreh_impuls Y(i) * p->achse X(i);

        p->achse X(i) += timestep * cross.x / p->traeg_moment[i];
        p->achse Y(i) += timestep * cross.y / p->traeg_moment[i];
        p->achse Z(i) += timestep * cross.z / p->traeg_moment[i];

        norm = sqrt( SPRODN(p->achse,i,p->achse,i) );
	    
        p->achse X(i) /= norm;
        p->achse Y(i) /= norm;
        p->achse Z(i) /= norm;
#endif    

#ifdef STRESS_TENS
#ifdef SHOCK
	/* plate against bulk */
        if (shock_mode == 1) {
          if ( p->ort X(i) < shock_strip )
            p->presstens[i].xx += (p->impuls X(i) - shock_speed * MASSE(p,i)) 
                    * (p->impuls X(i) - shock_speed * MASSE(p,i)) / MASSE(p,i);
          else
	    p->presstens[i].xx += p->impuls X(i) * p->impuls X(i)/MASSE(p,i);
        }
	/* two halves against one another */
        if (shock_mode == 2) {
          if ( p->ort X(i) < box_x.x*0.5 )
            p->presstens[i].xx += (p->impuls X(i) - shock_speed * MASSE(p,i)) 
                    * (p->impuls X(i) - shock_speed * MASSE(p,i)) / MASSE(p,i);
          else
            p->presstens[i].xx += (p->impuls X(i) + shock_speed * MASSE(p,i)) 
                    * (p->impuls X(i) + shock_speed * MASSE(p,i)) / MASSE(p,i);
        }
	/* bulk against wall */
        if (shock_mode == 3) p->presstens[i].xx += (p->impuls X(i) - 
	  shock_speed * MASSE(p,i)) * (p->impuls X(i) - shock_speed * 
           MASSE(p,i)) / MASSE(p,i);
#else
        p->presstens[i].xx += p->impuls X(i) * p->impuls X(i)/MASSE(p,i);
#endif
        p->presstens[i].yy += p->impuls Y(i) * p->impuls Y(i)/MASSE(p,i);
#ifndef TWOD
        p->presstens[i].zz += p->impuls Z(i) * p->impuls Z(i)/MASSE(p,i);
        p->presstens[i].yz += p->impuls Y(i) * p->impuls Z(i)/MASSE(p,i);
        p->presstens[i].zx += p->impuls Z(i) * p->impuls X(i)/MASSE(p,i);
#endif
        p->presstens[i].xy += p->impuls X(i) * p->impuls Y(i)/MASSE(p,i);
#endif /* STRESS_TENS */
    }
#if defined(EPITAX) && !defined(NVE)
    }
#endif
  }

#ifdef MPI
  /* add up results from different CPUs */
  tmpvec1[0] = tot_kin_energy;
  tmpvec1[1] = fnorm;
  tmpvec1[2] = PxF;

  MPI_Allreduce( tmpvec1, tmpvec2, 3, REAL, MPI_SUM, cpugrid);

  tot_kin_energy = tmpvec2[0];
  fnorm          = tmpvec2[1];
  PxF            = tmpvec2[2];
#endif
  
#ifdef AND
  /* Andersen Thermostat -- Initialize the velocities now and then */
  ++count;
  if ((tmp_interval!=0) && (0==count%tmp_interval)) maxwell(temperature);
#endif


}

#else

void move_atoms_nve(void) 
{
  if (myid==0)
  error("the chosen ensemble NVE is not supported by this binary");
}

#endif


/*****************************************************************************
*
*  NVE Integrator with microconvergence relaxation
*
*****************************************************************************/

#ifdef MIK

void move_atoms_mik(void)
{
  int k;
  real tmpvec1[2], tmpvec2[2];
  static int count = 0;
  tot_kin_energy = 0.0;
  fnorm = 0.0;

#ifdef AND
  /* Andersen Thermostat -- Initialize the velocities now and then */
  ++count;
  if ((tmp_interval!=0) && (0==count%tmp_interval)) maxwell(temperature);
#endif

  /* loop over all cells */
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tot_kin_energy,fnorm)
#endif
  for (k=0; k<ncells; ++k) {

    int  i, sort;
    cell *p;
    real kin_energie_1, kin_energie_2, tmp;
    p = cell_array + CELLS(k);

#ifdef PVPCRAY
#pragma ivdep
#endif
#ifdef SX4
#pragma vdir vector,nodep
#endif
    for (i=0; i<p->n; ++i) {

#ifdef EPITAX
      /* only substrate atoms are integrated by MIK */
      if ( (NUMMER(p,i)<=epitax_sub_n) || (p->pot_eng[i]<=(epitax_ctrl * epitax_poteng_min)) ) {
#endif

        kin_energie_1 = SPRODN(p->impuls,i,p->impuls,i);

	sort = VSORTE(p,i);
#ifdef FBC
        /* give virtual particles their extra force */
	p->kraft X(i) += (fbc_forces + sort)->x;
	p->kraft Y(i) += (fbc_forces + sort)->y;
#ifndef TWOD
	p->kraft Z(i) += (fbc_forces + sort)->z;
#endif
#endif /* FBC */

	/* and set their force (->momentum) in restricted directions to 0 */
	p->kraft X(i) *= (restrictions + sort)->x;
	p->kraft Y(i) *= (restrictions + sort)->y;
#ifndef TWOD
	p->kraft Z(i) *= (restrictions + sort)->z;
#endif
	

#ifdef FNORM
	fnorm +=  SPRODN(p->kraft,i,p->kraft,i) / MASSE(p,i);
#endif

        p->impuls X(i) += timestep * p->kraft X(i);
        p->impuls Y(i) += timestep * p->kraft Y(i);
#ifndef TWOD
        p->impuls Z(i) += timestep * p->kraft Z(i);
#endif

	/* Mikroconvergence Algorithm - set velocity zero if a*v < 0 */
	if (0.0 > SPRODN(p->impuls,i,p->kraft,i)) {
          p->impuls X(i) = 0.0;
          p->impuls Y(i) = 0.0;
#ifndef TWOD
          p->impuls Z(i) = 0.0;
#endif
        } else { /* new positions */
          tmp = timestep / MASSE(p,i);
          p->ort X(i) += tmp * p->impuls X(i);
          p->ort Y(i) += tmp * p->impuls Y(i);
#ifndef TWOD
          p->ort Z(i) += tmp * p->impuls Z(i);
#endif
        }

#ifdef ATNR
        if(p->nummer[i]==ATNR) {
          printf("impuls: %.16f  %.16f  %.16f \n",
                 p->impuls X(i),p->impuls Y(i),p->impuls Z(i));
          printf("ort:  %.16f  %.16f  %.16f \n",
                 p->ort X(i),p->ort Y(i),p->ort Z(i));
          fflush(stdout);
        }
#endif
        kin_energie_2 =  SPRODN(p->impuls,i,p->impuls,i);

        /* sum up kinetic energy on this CPU */ 
        tot_kin_energy += (kin_energie_1 + kin_energie_2) / (4.0 * MASSE(p,i));

#ifdef STRESS_TENS
        p->presstens[i].xx += p->impuls X(i) * p->impuls X(i)/MASSE(p,i);
        p->presstens[i].yy += p->impuls Y(i) * p->impuls Y(i)/MASSE(p,i);
#ifndef TWOD
        p->presstens[i].zz += p->impuls Z(i) * p->impuls Z(i)/MASSE(p,i);
        p->presstens[i].yz += p->impuls Y(i) * p->impuls Z(i)/MASSE(p,i);
        p->presstens[i].zx += p->impuls Z(i) * p->impuls X(i)/MASSE(p,i);
#endif
        p->presstens[i].xy += p->impuls X(i) * p->impuls Y(i)/MASSE(p,i);
#endif
    }
#ifdef EPITAX
    }
#endif
  }

#ifdef MPI
  /* add up results from different CPUs */
  tmpvec1[0] = tot_kin_energy;
  tmpvec1[1] = fnorm;

  MPI_Allreduce( tmpvec1, tmpvec2, 2, REAL, MPI_SUM, cpugrid);

  tot_kin_energy = tmpvec2[0];
  fnorm          = tmpvec2[1];
#endif

}

#else

void move_atoms_mik(void) 
{
  if (myid==0)
  error("the chosen ensemble MIK is not supported by this binary");
}

#endif


/*****************************************************************************
*
* NVT Integrator with Nose Hoover Thermostat 
*
*****************************************************************************/

#ifdef NVT

void move_atoms_nvt(void)

{
  int k;
  real tmpvec1[3], tmpvec2[3], ttt;
  real E_kin_1 = 0.0, E_kin_2 = 0.0;
  real reibung, eins_d_reib;
  real E_rot_1 = 0.0, E_rot_2 = 0.0;
#ifdef UNIAX
  real reibung_rot,  eins_d_reib_rot;
#endif
  fnorm = 0.0;

  reibung         =        1.0 - eta * timestep / 2.0;
  eins_d_reib     = 1.0 / (1.0 + eta * timestep / 2.0);
#ifdef UNIAX
  reibung_rot     =        1.0 - eta_rot * timestep / 2.0;
  eins_d_reib_rot = 1.0 / (1.0 + eta_rot * timestep / 2.0);
#endif
   
#ifdef _OPENMP
#pragma omp parallel for reduction(+:E_kin_1,E_kin_2,E_rot_1,E_rot_2,fnorm)
#endif
  for (k=0; k<ncells; ++k) {

    int i;
    int sort;
    cell *p;
    real tmp;
#ifdef UNIAX
    real dot, norm ;
    vektor cross ;
#endif
    p = cell_array + CELLS(k);

    for (i=0; i<p->n; ++i) {

#ifdef EPITAX
      /* only substrate atoms are integrated by NVT */
      if ( (NUMMER(p,i)<=epitax_sub_n) || (p->pot_eng[i]<=(epitax_ctrl * epitax_poteng_min)) ) {
#endif

        /* twice the old kinetic energy */
        E_kin_1 +=  SPRODN(p->impuls,i,p->impuls,i) / MASSE(p,i);
#ifdef UNIAX
        E_rot_1 +=  SPRODN(p->dreh_impuls,i,p->dreh_impuls,i) 
                                                  / p->traeg_moment[i];
#endif

	sort = VSORTE(p,i);
#ifdef FBC
        /* give virtual particles their extra force */
	p->kraft X(i) += (fbc_forces + sort)->x;
	p->kraft Y(i) += (fbc_forces + sort)->y;
#ifndef TWOD
	p->kraft Z(i) += (fbc_forces + sort)->z;
#endif

#endif
	p->kraft X(i) *= (restrictions + sort)->x;
	p->kraft Y(i) *= (restrictions + sort)->y;
#ifndef TWOD
	p->kraft Z(i) *= (restrictions + sort)->z;
#endif
#ifdef FNORM
	fnorm +=  SPRODN(p->kraft,i,p->kraft,i) / MASSE(p,i);
#endif

	p->impuls X(i) = (p->impuls X(i) * reibung + timestep * p->kraft X(i)) 
                           * eins_d_reib * (restrictions + sort)->x;
        p->impuls Y(i) = (p->impuls Y(i) * reibung + timestep * p->kraft Y(i)) 
                           * eins_d_reib * (restrictions + sort)->y;
#ifndef TWOD
        p->impuls Z(i) = (p->impuls Z(i) * reibung + timestep * p->kraft Z(i)) 
                           * eins_d_reib * (restrictions + sort)->z;
#endif

#ifdef SLLOD
	p->impuls X(i) += epsilon * p->impuls Y(i) / MASSE(p,i);
#endif

#ifdef UNIAX
        /* new angular momenta */
        dot = 2.0 * SPRODN(p->dreh_impuls,i,p->achse,i);

        p->dreh_impuls X(i) = eins_d_reib_rot
            * ( p->dreh_impuls X(i) * reibung_rot
                + timestep * p->dreh_moment X(i) - dot * p->achse X(i) );
        p->dreh_impuls Y(i) = eins_d_reib_rot
            * ( p->dreh_impuls Y(i) * reibung_rot
                + timestep * p->dreh_moment Y(i) - dot * p->achse Y(i) );
        p->dreh_impuls Z(i) = eins_d_reib_rot
            * ( p->dreh_impuls Z(i) * reibung_rot
                + timestep * p->dreh_moment Z(i) - dot * p->achse Z(i) );
#endif

        /* twice the new kinetic energy */ 
        E_kin_2 += SPRODN(p->impuls,i,p->impuls,i) / MASSE(p,i);
#ifdef UNIAX
        E_rot_2 += SPRODN(p->dreh_impuls,i,p->dreh_impuls,i) 
                                                     / p->traeg_moment[i];
#endif

        /* new positions */
        tmp = timestep / MASSE(p,i);
        p->ort X(i) += tmp * p->impuls X(i);
        p->ort Y(i) += tmp * p->impuls Y(i);
#ifndef TWOD
        p->ort Z(i) += tmp * p->impuls Z(i);
#endif

#ifdef SLLOD
	p->ort X(i) += epsilon * p->ort Y(i);
#endif /* SLLOD */

#ifdef UNIAX
        cross.x = p->dreh_impuls Y(i) * p->achse Z(i)
                - p->dreh_impuls Z(i) * p->achse Y(i);
        cross.y = p->dreh_impuls Z(i) * p->achse X(i)
                - p->dreh_impuls X(i) * p->achse Z(i);
        cross.z = p->dreh_impuls X(i) * p->achse Y(i)
                - p->dreh_impuls Y(i) * p->achse X(i);

        p->achse X(i) += timestep * cross.x / p->traeg_moment[i];
        p->achse Y(i) += timestep * cross.y / p->traeg_moment[i];
        p->achse Z(i) += timestep * cross.z / p->traeg_moment[i];

        norm = sqrt( SPRODN(p->achse,i,p->achse,i) );

        p->achse X(i) /= norm;
        p->achse Y(i) /= norm;
        p->achse Z(i) /= norm;
#endif

#ifdef STRESS_TENS
        p->presstens[i].xx += p->impuls X(i) * p->impuls X(i)/MASSE(p,i);
        p->presstens[i].yy += p->impuls Y(i) * p->impuls Y(i)/MASSE(p,i);
#ifndef TWOD
        p->presstens[i].zz += p->impuls Z(i) * p->impuls Z(i)/MASSE(p,i);
        p->presstens[i].yz += p->impuls Y(i) * p->impuls Z(i)/MASSE(p,i);
        p->presstens[i].zx += p->impuls Z(i) * p->impuls X(i)/MASSE(p,i);
#endif
        p->presstens[i].xy += p->impuls X(i) * p->impuls Y(i)/MASSE(p,i);
#endif
    }
#ifdef EPITAX
    }
#endif
  }
  
#ifdef SLLOD
  /* new box size */
  box_y.x += epsilon * box_y.y;
  make_box();
#endif

#ifdef UNIAX
  tot_kin_energy = ( E_kin_1 + E_kin_2 + E_rot_1 + E_rot_2 ) / 4.0;
#else
  tot_kin_energy = ( E_kin_1 + E_kin_2 ) / 4.0;
#endif

#ifdef MPI
  /* add up results from different CPUs */
  tmpvec1[0] = tot_kin_energy;
  tmpvec1[1] = E_kin_2;
  tmpvec1[2] = E_rot_2;

  MPI_Allreduce( tmpvec1, tmpvec2, 3, REAL, MPI_SUM, cpugrid);

  tot_kin_energy = tmpvec2[0];
  E_kin_2        = tmpvec2[1];
  E_rot_2        = tmpvec2[2];
#endif

  /* time evolution of constraints */
  ttt  = nactive * temperature;
  eta += timestep * (E_kin_2 / ttt - 1.0) * isq_tau_eta;
#ifdef UNIAX
  ttt  = nactive_rot * temperature;
  eta_rot += timestep * (E_rot_2 / ttt - 1.0) * isq_tau_eta_rot;
#endif
  
}

#else

void move_atoms_nvt(void) 
{
  if (myid==0)
  error("the chosen ensemble NVT is not supported by this binary");
}

#endif


#ifdef NPT

/******************************************************************************
*
*  compute initial dynamical pressure
*
******************************************************************************/

void calc_dyn_pressure(void)
{
  int  k;
  real tmpvec1[5], tmpvec2[5];

  /* initialize data */
  dyn_stress_x = 0.0;
  dyn_stress_y = 0.0;
  dyn_stress_z = 0.0;
  Ekin_old     = 0.0;
  Erot_old     = 0.0;

  /* loop over all cells */
#ifdef _OPENMP
#pragma omp parallel for reduction(+:dyn_stress_x,dyn_stress_y,dyn_stress_z,Ekin_old,Erot_old)
#endif
  for (k=0; k<ncells; ++k) {

    int i;
    cell *p;
    real tmp;
    p = cell_array + CELLS(k);

    /* loop over atoms in cell */
    for (i=0; i<p->n; ++i) {
      tmp = 1.0 / MASSE(p,i);
      dyn_stress_x += p->impuls X(i) * p->impuls X(i) * tmp;
      dyn_stress_y += p->impuls Y(i) * p->impuls Y(i) * tmp;
#ifndef TWOD
      dyn_stress_z += p->impuls Z(i) * p->impuls Z(i) * tmp;
#endif
#ifdef UNIAX
      Erot_old += SPRODN(p->dreh_impuls,i,p->dreh_impuls,i)/p->traeg_moment[i];
#endif
    }
  }

  Ekin_old  = dyn_stress_x + dyn_stress_y;
#ifndef TWOD
  Ekin_old += dyn_stress_z;
#endif

#ifdef MPI
  /* add up results from different CPUs */
  tmpvec1[0]   = dyn_stress_x;
  tmpvec1[1]   = dyn_stress_y;
  tmpvec1[2]   = dyn_stress_z;
  tmpvec1[3]   = Ekin_old;
  tmpvec1[4]   = Erot_old;

  MPI_Allreduce( tmpvec1, tmpvec2, 5, REAL, MPI_SUM, cpugrid);

  dyn_stress_x = tmpvec2[0];
  dyn_stress_y = tmpvec2[1];
  dyn_stress_z = tmpvec2[2];
  Ekin_old     = tmpvec2[3];
  Erot_old     = tmpvec2[4];
#endif

}

#endif /* NPT */

/******************************************************************************
*
* NPT Integrator with Nose Hoover Thermostat
*
******************************************************************************/

#ifdef NPT_iso

void move_atoms_npt_iso(void)
{
  int  k;
  real Ekin_new = 0.0, Erot_new = 0.0;
  real pfric, pifric, rfric, rifric;
  real tmpvec1[5], tmpvec2[5], ttt;
  real reib, ireib;

  fnorm    = 0.0;
#ifdef UNIAX
  pressure = (0.6 * (Ekin_old + Erot_old) + virial) / (DIM * volume);
#else
  pressure = (Ekin_old + virial) / (DIM * volume) ;
#endif

  /* time evolution of xi */
  xi_old.x = xi.x;
  xi.x += timestep * (pressure-pressure_ext.x) * volume * isq_tau_xi / nactive;

  /* some constants used later on */
  pfric  =        1.0 - (xi_old.x + eta) * timestep / 2.0;
  pifric = 1.0 / (1.0 + (xi.x     + eta) * timestep / 2.0);
  rfric  =        1.0 + (xi.x          ) * timestep / 2.0;
  rifric = 1.0 / (1.0 - (xi.x          ) * timestep / 2.0);
#ifdef UNIAX
  reib  =        1.0 - eta_rot * timestep / 2.0;
  ireib = 1.0 / (1.0 + eta_rot * timestep / 2.0);
#endif

  /* loop over all cells */
#ifdef _OPENMP
#pragma omp parallel for reduction(+:Ekin_new,Erot_new,fnorm)
#endif
  for (k=0; k<ncells; ++k) {

    int i;
    cell *p;
    real tmp;
#ifdef UNIAX
    real dot, norm ;
    vektor cross ;
#endif
    p = cell_array + CELLS(k);

    for (i=0; i<p->n; ++i) {

#ifdef FNORM
      fnorm +=  SPRODN(p->kraft,i,p->kraft,i) / MASSE(p,i);
#endif

      /* new momenta */
      p->impuls X(i) = (pfric*p->impuls X(i)+timestep*p->kraft X(i))*pifric;
      p->impuls Y(i) = (pfric*p->impuls Y(i)+timestep*p->kraft Y(i))*pifric;
#ifndef TWOD
      p->impuls Z(i) = (pfric*p->impuls Z(i)+timestep*p->kraft Z(i))*pifric;
#endif

#ifdef UNIAX
      /* new angular momenta */
      dot = 2.0 * SPRODN(p->dreh_impuls,i,p->achse,i);
      p->dreh_impuls X(i) = ireib * ( p->dreh_impuls X(i) * reib 
              + timestep * p->dreh_moment X(i) - dot * p->achse X(i) );
      p->dreh_impuls Y(i) = ireib * ( p->dreh_impuls Y(i) * reib 
              + timestep * p->dreh_moment Y(i) - dot * p->achse Y(i) );
      p->dreh_impuls Z(i) = ireib * ( p->dreh_impuls Z(i) * reib 
              + timestep * p->dreh_moment Z(i) - dot * p->achse Z(i) );
#endif

      /* twice the new kinetic energy */ 
      Ekin_new += SPRODN(p->impuls,i,p->impuls,i) / MASSE(p,i);
#ifdef UNIAX
      Erot_new += SPRODN(p->dreh_impuls,i,p->dreh_impuls,i)/p->traeg_moment[i];
#endif

      /* new positions */
      tmp = timestep / MASSE(p,i);
      p->ort X(i) = (rfric * p->ort X(i) + p->impuls X(i) * tmp) * rifric;
      p->ort Y(i) = (rfric * p->ort Y(i) + p->impuls Y(i) * tmp) * rifric;
#ifndef TWOD
      p->ort Z(i) = (rfric * p->ort Z(i) + p->impuls Z(i) * tmp) * rifric;
#endif

#ifdef UNIAX
      /* new molecular axes */
      cross.x = p->dreh_impuls Y(i) * p->achse Z(i)
              - p->dreh_impuls Z(i) * p->achse Y(i);
      cross.y = p->dreh_impuls Z(i) * p->achse X(i)
              - p->dreh_impuls X(i) * p->achse Z(i);
      cross.z = p->dreh_impuls X(i) * p->achse Y(i)
              - p->dreh_impuls Y(i) * p->achse X(i);

      p->achse X(i) += timestep * cross.x / p->traeg_moment[i];
      p->achse Y(i) += timestep * cross.y / p->traeg_moment[i];
      p->achse Z(i) += timestep * cross.z / p->traeg_moment[i];

      norm = sqrt( SPRODN(p->achse,i,p->achse,i) );

      p->achse X(i) /= norm;
      p->achse Y(i) /= norm;
      p->achse Z(i) /= norm;
#endif

#ifdef STRESS_TENS
        p->presstens[i].xx += p->impuls X(i) * p->impuls X(i)/MASSE(p,i);
        p->presstens[i].yy += p->impuls Y(i) * p->impuls Y(i)/MASSE(p,i);
#ifndef TWOD
        p->presstens[i].zz += p->impuls Z(i) * p->impuls Z(i)/MASSE(p,i);
        p->presstens[i].yz += p->impuls Y(i) * p->impuls Z(i)/MASSE(p,i);
        p->presstens[i].zx += p->impuls Z(i) * p->impuls X(i)/MASSE(p,i);
#endif
        p->presstens[i].xy += p->impuls X(i) * p->impuls Y(i)/MASSE(p,i);
#endif
    }
  }

#ifdef MPI
  /* add up results from all CPUs */
  tmpvec1[0] = Ekin_new;
  tmpvec1[1] = Erot_new;
  tmpvec1[2] = fnorm;

  MPI_Allreduce( tmpvec1, tmpvec2, 3, REAL, MPI_SUM, cpugrid);

  Ekin_new = tmpvec2[0];
  Erot_new = tmpvec2[1];
  fnorm    = tmpvec2[2];
#endif

#ifdef UNIAX
  tot_kin_energy = ( Ekin_old + Ekin_new + Erot_old + Erot_new) / 4.0;
#else
  tot_kin_energy = ( Ekin_old + Ekin_new ) / 4.0;
#endif

  /* time evolution of eta */
  ttt      = nactive * temperature;
  eta     += timestep * (Ekin_new / ttt - 1.0) * isq_tau_eta;
  Ekin_old = Ekin_new;
#ifdef UNIAX
  ttt      = nactive_rot * temperature;
  eta_rot += timestep * (Erot_new / ttt - 1.0) * isq_tau_eta_rot;
  Erot_old = Erot_new;
#endif

  /* time evolution of box size */
  ttt = (1.0 + xi.x * timestep / 2.0) / (1.0 - xi.x * timestep / 2.0);
  if (ttt<0) error("box size has become negative!");
  box_x.x *= ttt;
  box_x.y *= ttt;
  box_y.x *= ttt;
  box_y.y *= ttt;
#ifndef TWOD
  box_x.z *= ttt;
  box_y.z *= ttt;
  box_z.x *= ttt;
  box_z.y *= ttt;
  box_z.z *= ttt;
#endif  
  make_box();
}

#else

void move_atoms_npt_iso(void) 
{
  if (myid==0)
  error("the chosen ensemble NPT_ISO is not supported by this binary");
}

#endif


/******************************************************************************
*
*  NPT Integrator with Nose Hoover Thermostat
*
******************************************************************************/

#ifdef NPT_axial

void move_atoms_npt_axial(void)
{
  int k;
  real Ekin_new = 0.0, ttt, tmpvec1[5], tmpvec2[5];
  vektor pfric, pifric, rfric, rifric, tvec;

  fnorm    = 0.0;
  stress_x = (dyn_stress_x + vir_xx) / volume;  dyn_stress_x = 0.0;
  stress_y = (dyn_stress_y + vir_yy) / volume;  dyn_stress_y = 0.0;
#ifndef TWOD
  stress_z = (dyn_stress_z + vir_zz) / volume;  dyn_stress_z = 0.0;
#endif

  /* time evolution of xi */
  ttt  = timestep * volume * isq_tau_xi / nactive;
  xi_old.x = xi.x;  xi.x += ttt * (stress_x - pressure_ext.x);
  xi_old.y = xi.y;  xi.y += ttt * (stress_y - pressure_ext.y);
#ifndef TWOD
  xi_old.z = xi.z;  xi.z += ttt * (stress_z - pressure_ext.z);
#endif

  /* some constants used later on */
  pfric.x  =        1.0 - (xi_old.x + eta) * timestep / 2.0;
  pifric.x = 1.0 / (1.0 + (xi.x     + eta) * timestep / 2.0);
  rfric.x  =        1.0 + (xi.x          ) * timestep / 2.0;
  rifric.x = 1.0 / (1.0 - (xi.x          ) * timestep / 2.0);
  pfric.y  =        1.0 - (xi_old.y + eta) * timestep / 2.0;
  pifric.y = 1.0 / (1.0 + (xi.y     + eta) * timestep / 2.0);
  rfric.y  =        1.0 + (xi.y          ) * timestep / 2.0;
  rifric.y = 1.0 / (1.0 - (xi.y          ) * timestep / 2.0);
#ifndef TWOD
  pfric.z  =        1.0 - (xi_old.z + eta) * timestep / 2.0;
  pifric.z = 1.0 / (1.0 + (xi.z     + eta) * timestep / 2.0);
  rfric.z  =        1.0 + (xi.z          ) * timestep / 2.0;
  rifric.z = 1.0 / (1.0 - (xi.z          ) * timestep / 2.0);
#endif

  /* loop over all cells */
#ifdef _OPENMP
#pragma omp parallel for reduction(+:Ekin,dyn_stress_x,dyn_stress_y,dyn_stress_z,fnorm)
#endif
  for (k=0; k<ncells; ++k) {

    int i;
    cell *p;
    real tmp;
    p = cell_array + CELLS(k);

    /* loop over atoms in cell */
    for (i=0; i<p->n; ++i) {

#ifdef FNORM
      fnorm +=  SPRODN(p->kraft,i,p->kraft,i) / MASSE(p,i);
#endif

      /* new momenta */
      p->impuls X(i) = (pfric.x * p->impuls X(i)
                                   + timestep * p->kraft X(i)) * pifric.x;
      p->impuls Y(i) = (pfric.y * p->impuls Y(i)
                                   + timestep * p->kraft Y(i)) * pifric.y;
#ifndef TWOD
      p->impuls Z(i) = (pfric.z * p->impuls Z(i)
                                   + timestep * p->kraft Z(i)) * pifric.z;
#endif

      /* new stress tensor (dynamic part only) */
      tmp = 1.0 / MASSE(p,i);
      dyn_stress_x += p->impuls X(i) * p->impuls X(i) * tmp;
      dyn_stress_y += p->impuls Y(i) * p->impuls Y(i) * tmp;
#ifndef TWOD
      dyn_stress_z += p->impuls Z(i) * p->impuls Z(i) * tmp;
#endif

      /* twice the new kinetic energy */ 
      Ekin_new += SPRODN(p->impuls,i,p->impuls,i) / MASSE(p,i);
	  
      /* new positions */
      tmp = timestep / MASSE(p,i);
      p->ort X(i) = (rfric.x * p->ort X(i) + p->impuls X(i) * tmp) * rifric.x;
      p->ort Y(i) = (rfric.y * p->ort Y(i) + p->impuls Y(i) * tmp) * rifric.y;
#ifndef TWOD
      p->ort Z(i) = (rfric.z * p->ort Z(i) + p->impuls Z(i) * tmp) * rifric.z;
#endif

#ifdef STRESS_TENS
        p->presstens[i].xx += p->impuls X(i) * p->impuls X(i)/MASSE(p,i);
        p->presstens[i].yy += p->impuls Y(i) * p->impuls Y(i)/MASSE(p,i);
#ifndef TWOD
        p->presstens[i].zz += p->impuls Z(i) * p->impuls Z(i)/MASSE(p,i);
        p->presstens[i].yz += p->impuls Y(i) * p->impuls Z(i)/MASSE(p,i);
        p->presstens[i].zx += p->impuls Z(i) * p->impuls X(i)/MASSE(p,i);
#endif
        p->presstens[i].xy += p->impuls X(i) * p->impuls Y(i)/MASSE(p,i);
#endif
    }
  }

#ifdef MPI
  /* add up results from different CPUs */
  tmpvec1[0]   = Ekin_new;
  tmpvec1[1]   = fnorm;
  tmpvec1[2]   = dyn_stress_x;
  tmpvec1[3]   = dyn_stress_y;
  tmpvec1[4]   = dyn_stress_z;
  MPI_Allreduce( tmpvec1, tmpvec2, 5, REAL, MPI_SUM, cpugrid);
  Ekin_new     = tmpvec2[0];
  fnorm        = tmpvec2[1];
  dyn_stress_x = tmpvec2[2];
  dyn_stress_y = tmpvec2[3];
  dyn_stress_z = tmpvec2[4];
#endif

  /* time evolution of eta */
  tot_kin_energy = ( Ekin_old + Ekin_new ) / 4.0;
  ttt      = nactive * temperature;
  eta     += timestep * (Ekin_new / ttt - 1.0) * isq_tau_eta;
  Ekin_old = Ekin_new;

  /* time evolution of box size */
  tvec.x   = (1.0 + xi.x * timestep / 2.0) / (1.0 - xi.x * timestep / 2.0);
  tvec.y   = (1.0 + xi.y * timestep / 2.0) / (1.0 - xi.y * timestep / 2.0);
  if ((tvec.x<0) || (tvec.y<0)) error("box size has become negative!");
  box_x.x *= tvec.x;
  box_x.y *= tvec.x;
  box_y.x *= tvec.y;
  box_y.y *= tvec.y;
#ifndef TWOD
  tvec.z = (1.0 + xi.z * timestep / 2.0) / (1.0 - xi.z * timestep / 2.0);
  if (tvec.z<0) error("box size has become negative!");
  box_x.z *= tvec.x;
  box_y.z *= tvec.y;
  box_z.x *= tvec.z;
  box_z.y *= tvec.z;
  box_z.z *= tvec.z;
#endif
  make_box();
}

#else

void move_atoms_npt_axial(void) 
{
  if (myid==0)
  error("the chosen ensemble NPT_AXIAL is not supported by this binary");
}

#endif


/*****************************************************************************
*
*  NVE Integrator with stadium damping and fixed borders 
*  for fracture studies
*
*****************************************************************************/

#ifdef FRAC

void move_atoms_frac(void)

{
  int k;
  real tmpvec1[6], tmpvec2[6], ttt;
  real E_kin_1        = 0.0, E_kin_2        = 0.0; 
  real E_kin_damp1    = 0.0, E_kin_damp2    = 0.0;
  real E_kin_stadium1 = 0.0, E_kin_stadium2 = 0.0;
  real reibung, reibung_y, eins_d_reib, eins_d_reib_y;
  real epsilontmp, eins_d_epsilontmp;

  real f; /* stadium function: the bath tub !!!!*/

  fnorm     = 0.0;
  sum_f     = 0.0;
  n_stadium = 0;

  if(expansionmode==1)
      dotepsilon = dotepsilon0 / (1.0 + dotepsilon0 * steps * timestep);
      

  /* loop over all atoms */
#ifdef _OPENMP
#pragma omp parallel for reduction(+:E_kin_1,E_kin_2,E_kin_damp1,E_kin_damp2,E_kin_stadium1,E_kin_stadium2,sum_f,n_stadium,fnorm)
#endif
  for (k=0; k<ncells; ++k) {

    int i;
    int sort;
    cell *p;
    real tmp,tmp1,tmp2;

    p = cell_array + CELLS(k);

    for (i=0; i<p->n; ++i) {
	
	/* dampingmode == 2 -> global viscous damping !!!! */
	if(dampingmode == 2){ 
	    f = 1.0; 
	} else {
	    /* Calculate stadium function f */
	    tmp1 = SQR((p->ort X(i)-center.x)/box_x.x);
	    tmp2 = SQR((p->ort Y(i)-center.y)/box_y.y);
	    f = (tmp1+tmp2-SQR(stadium.x/box_x.x))/\
		(.25- SQR(stadium.x/box_x.x));
	}
	

	if (f<= 0.0) {
	    f = 0.0;
	    n_stadium += DIM;
	}
	if (f>1.0) f = 1.0;

	sort = VSORTE(p,i);

        /* add up f considering the restriction vector  */
#ifdef TWOD
	sum_f+= f * ( (restrictions + sort)->x + 
		      (restrictions + sort)->y   )/2.0;
#else
	sum_f+= f * ( (restrictions + sort)->x + 
		      (restrictions + sort)->y +  
		      (restrictions + sort)->z  )/3.0;
#endif
	
       	/* twice the old kinetic energy */
	E_kin_1 +=  SPRODN(p->impuls,i,p->impuls,i) / MASSE(p,i);	
	if (f == 0.0)
	    E_kin_stadium1 +=      SPRODN(p->impuls,i,p->impuls,i) / MASSE(p,i);
	if (f >  0.0)
	    E_kin_damp1    +=  f * SPRODN(p->impuls,i,p->impuls,i) / MASSE(p,i);
	
#ifdef FBC
        /* give virtual particles their extra force */
	p->kraft X(i) += (fbc_forces + sort)->x;
	p->kraft Y(i) += (fbc_forces + sort)->y;
#ifndef TWOD
	p->kraft Z(i) += (fbc_forces + sort)->z;
#endif
#endif

	p->kraft X(i) *= (restrictions + sort)->x;
	p->kraft Y(i) *= (restrictions + sort)->y;
#ifndef TWOD
	p->kraft Z(i) *= (restrictions + sort)->z;
#endif

#ifdef FNORM
	fnorm +=  SPRODN(p->kraft,i,p->kraft,i) / MASSE(p,i);
#endif

	reibung       =        1.0 -  gamma_damp * f * timestep / 2.0;
	eins_d_reib   = 1.0 / (1.0 +  gamma_damp * f * timestep / 2.0);
	reibung_y     =        1.0 - (gamma_damp * f + dotepsilon) * 
	                        timestep / 2.0;
	eins_d_reib_y = 1.0 / (1.0 + (gamma_damp * f + dotepsilon) * 
			        timestep / 2.0);
	
        /* new momenta */
	p->impuls X(i) = (p->impuls X(i)   * reibung   + timestep * p->kraft X(i)) 
                           * eins_d_reib   * (restrictions + sort)->x;
        p->impuls Y(i) = (p->impuls Y(i)   * reibung_y + timestep * p->kraft Y(i)) 
                           * eins_d_reib_y * (restrictions + sort)->y;
#ifndef TWOD
        p->impuls Z(i) = (p->impuls Z(i)   * reibung   + timestep * p->kraft Z(i)) 
                           * eins_d_reib   * (restrictions + sort)->z;
#endif                  

	/* twice the new kinetic energy */ 
	E_kin_2 +=  SPRODN(p->impuls,i,p->impuls,i) / MASSE(p,i);	
	
	if (f == 0.0)
	    E_kin_stadium2 +=      SPRODN(p->impuls,i,p->impuls,i) / MASSE(p,i);
	if (f > 0.0)
	    E_kin_damp2    +=  f * SPRODN(p->impuls,i,p->impuls,i) / MASSE(p,i);


	/* new positions */
        tmp = timestep / MASSE(p,i);
	epsilontmp =               1.0 + dotepsilon * timestep / 2.0;
	eins_d_epsilontmp = 1.0 / (1.0 - dotepsilon * timestep / 2.0);

        p->ort X(i) +=  tmp * p->impuls X(i);
        p->ort Y(i)  = (tmp * p->impuls Y(i) + epsilontmp * p->ort Y(i))
	                * eins_d_epsilontmp;

#ifndef TWOD
        p->ort Z(i) +=  tmp * p->impuls Z(i);
#endif

#ifdef STRESS_TENS
        p->presstens[i].xx += p->impuls X(i) * p->impuls X(i)/MASSE(p,i);
        p->presstens[i].yy += p->impuls Y(i) * p->impuls Y(i)/MASSE(p,i);
#ifndef TWOD
        p->presstens[i].zz += p->impuls Z(i) * p->impuls Z(i)/MASSE(p,i);
        p->presstens[i].yz += p->impuls Y(i) * p->impuls Z(i)/MASSE(p,i);
        p->presstens[i].zx += p->impuls Z(i) * p->impuls X(i)/MASSE(p,i);
#endif
        p->presstens[i].xy += p->impuls X(i) * p->impuls Y(i)/MASSE(p,i);
#endif
    }
  }

  tot_kin_energy = ( E_kin_1        + E_kin_2         ) / 4.0;
  E_kin_stadium  = ( E_kin_stadium1 + E_kin_stadium2  ) / 4.0;
  E_kin_damp     = ( E_kin_damp1    + E_kin_damp2     ) / 4.0;


#ifdef MPI
  /* add up results from different CPUs */
  tmpvec1[0] = tot_kin_energy;
  tmpvec1[1] = E_kin_stadium;
  tmpvec1[2] = E_kin_damp;
  tmpvec1[3] = E_kin_damp2;
  tmpvec1[4] = n_stadium;
  tmpvec1[5] = sum_f;

  MPI_Allreduce( tmpvec1, tmpvec2, 6, REAL, MPI_SUM, cpugrid);

  tot_kin_energy = tmpvec2[0];
  E_kin_stadium  = tmpvec2[1];
  E_kin_damp     = tmpvec2[2];
  E_kin_damp2    = tmpvec2[3];
  n_stadium      = tmpvec2[4];
  sum_f          = tmpvec2[5];
#endif

  ttt   = 2.0 * temperature * sum_f;

  /* time evolution of constraints */
  /* dampingmode: 0 -> viscous damping (default); 
                  1 -> Nose-Hoover; 
		  2 -> global viscous damping */

  if(dampingmode == 1){
      gamma_damp += timestep * (E_kin_damp2 / ttt - 1.0) * gamma_bar;
  } else {
      gamma_damp  =            (1.0 - ttt / E_kin_damp2) * gamma_bar;
  }
      
}

#else

void move_atoms_frac(void) 
{
  if (myid==0)
  error("the chosen ensemble FRAC is not supported by this binary");
}

#endif

#ifdef STM

/*****************************************************************************
*
*  NVT Integrator with Stadium 
*
*****************************************************************************/

void move_atoms_stm(void)

{
  int k;
  /* we handle 2 ensembles ensindex = 0 -> NVT ;ensindex = 1 -> NVE */
  int ensindex = 0;
  real kin_energie_1[2] = {0.0,0.0}, kin_energie_2[2] = {0.0,0.0};
  real tmpvec1[5], tmpvec2[5], ttt;
  n_stadium = 0;

  /* loop over all atoms */
#ifdef _OPENMP
#pragma omp parallel for reduction(+:kin_energie_1[0],kin_energie_1[1],kin_energie_2[0],kin_energie_2[2],n_stadium)
#endif
  for (k=0; k<ncells; ++k) {

    int i;
    cell *p;
    real reibung, eins_d_reib;
    real tmp;
    vektor d;
    int sort=0;

    p = cell_array + CELLS(k);

    for (i=0; i<p->n; ++i) {

        /* Check if outside or inside the ellipse: */	
        tmp = SQR((p->ort X(i)-center.x)/stadium.x) +
              SQR((p->ort Y(i)-center.y)/stadium.y) - 1;
        if (tmp <= 0) {
          /* We are inside the ellipse: */
          reibung = 1.0;
          eins_d_reib = 1.0;
	  n_stadium += DIM;
	  ensindex = 1;
        } else {
          reibung     =      1 - eta * timestep / 2.0;
          eins_d_reib = 1 / (1 + eta * timestep / 2.0);
	  ensindex = 0;
        }

        /* twice the old kinetic energy */
        kin_energie_1[ensindex] +=  SPRODN(p->impuls,i,p->impuls,i) / MASSE(p,i);

        /* new momenta */
	sort = VSORTE(p,i);
	p->impuls X(i) = (p->impuls X(i)*reibung + timestep * p->kraft X(i))
                          * eins_d_reib * (restrictions + sort)->x;
        p->impuls Y(i) = (p->impuls Y(i)*reibung + timestep * p->kraft Y(i))
                          * eins_d_reib * (restrictions + sort)->y;

        /* twice the new kinetic energy */ 
        kin_energie_2[ensindex] +=  SPRODN(p->impuls,i,p->impuls,i) / MASSE(p,i);

        /* new positions */
        tmp = timestep * MASSE(p,i);
        p->ort X(i) += tmp * p->impuls X(i);
        p->ort Y(i) += tmp * p->impuls Y(i);
    }
  }
  
  tot_kin_energy  = (kin_energie_1[0] + kin_energie_2[0]) / 4.0;
  E_kin_stadium   = (kin_energie_1[1] + kin_energie_2[1]) / 4.0;
#ifdef MPI
  /* add up results from all CPUs */
  tmpvec1[0] = tot_kin_energy;
  tmpvec1[1] = kin_energie_2[0];
  tmpvec1[2] = E_kin_stadium;
  tmpvec1[3] = kin_energie_2[1];
  tmpvec1[4] = (real)n_stadium;
  MPI_Allreduce( tmpvec1, tmpvec2, 5, REAL, MPI_SUM, cpugrid);

  tot_kin_energy     = tmpvec2[0];
  kin_energie_2[0]   = tmpvec2[1];
  E_kin_stadium = tmpvec2[2];
  kin_energie_2[1]   = tmpvec2[3];
  n_stadium              = (int)tmpvec2[4];
#endif

  /* Zeitentwicklung der Parameter */
  ttt  = (nactive - n_stadium) * temperature;
  eta += timestep * (kin_energie_2[0] / ttt - 1.0) * isq_tau_eta;
}

#else

void move_atoms_stm(void) 
{
  if (myid==0)
  error("the chosen ensemble STM is not supported by this binary");
}

#endif

/******************************************************************************
*
*  NVX Integrator for heat conductivity
* 
******************************************************************************/

#ifdef NVX

void move_atoms_nvx(void)

{
  int  k;
  real Ekin_1, Ekin_2;
  real Ekin_left = 0.0, Ekin_right = 0.0;
  int  natoms_left = 0, natoms_right = 0;
  real px, vol, real_tmp;
  int  num, nhalf, int_tmp;  
  real scale, rescale, Rescale;
  vektor tot_impuls_left, tot_impuls_right, vectmp;
  real inv_mass_left=0.0, inv_mass_right=0.0;
 
  tot_kin_energy = 0.0;
  tot_impuls_left.x  = 0.0;
  tot_impuls_right.x = 0.0;
  tot_impuls_left.y  = 0.0;
  tot_impuls_right.y = 0.0;
#ifndef TWOD
  tot_impuls_left.z  = 0.0;
  tot_impuls_right.z = 0.0;
#endif

  nhalf = tran_nlayers / 2;
  scale = tran_nlayers / box_x.x;

  /* loop over all atoms */
  for (k=0; k<ncells; ++k) {

    int i;
    cell *p;
    p = cell_array + CELLS(k);

    for (i=0; i<p->n; ++i) {

      Ekin_1 = SPRODN(p->impuls,i,p->impuls,i) / MASSE(p,i); 
      px = p->impuls X(i);

      /* new momenta */
      p->impuls X(i) += timestep * p->kraft X(i); 
      p->impuls Y(i) += timestep * p->kraft Y(i); 
#ifndef TWOD
      p->impuls Z(i) += timestep * p->kraft Z(i); 
#endif

      /* twice the new kinetic energy */ 
      Ekin_2 = SPRODN(p->impuls,i,p->impuls,i) / MASSE(p,i); 
      px = (px + p->impuls X(i)) / 2.0;

      /* new positions */
      p->ort X(i) += timestep * p->impuls X(i) / MASSE(p,i);
      p->ort Y(i) += timestep * p->impuls Y(i) / MASSE(p,i);
#ifndef TWOD
      p->ort Z(i) += timestep * p->impuls Z(i) / MASSE(p,i);
#endif

      tot_kin_energy += (Ekin_1 + Ekin_2) / 4.0;

      /* which layer */
      num = scale * p->ort X(i);
      if (num < 0)             num = 0;
      if (num >= tran_nlayers) num = tran_nlayers-1;

      /* temperature control and heat conductivity */
      if (num == 0) {
        Ekin_left += Ekin_2;
        inv_mass_left     += 1/(MASSE(p,i));
        tot_impuls_left.x += p->impuls X(i);
        tot_impuls_left.y += p->impuls Y(i);
#ifndef TWOD
        tot_impuls_left.z += p->impuls Z(i);
#endif
        natoms_left++;
      } else if  (num == nhalf) {
        Ekin_right += Ekin_2;
        inv_mass_right     += 1/(MASSE(p,i));
        tot_impuls_right.x += p->impuls X(i);
        tot_impuls_right.y += p->impuls Y(i);
#ifndef TWOD
        tot_impuls_right.z += p->impuls Z(i);
#endif
        natoms_right++;
      } else if (num < nhalf) {
        heat_cond += (p->heatcond[i] + (Ekin_1 + Ekin_2)/2.0 )
                                  * px / MASSE(p,i);
      } else {
        heat_cond -= (p->heatcond[i] + (Ekin_1 + Ekin_2)/2.0 )
                                  * px / MASSE(p,i);
      }
    }
  }

#ifdef MPI
  /* Add up results from all cpus */
  MPI_Allreduce( &tot_kin_energy, &real_tmp, 1, REAL, MPI_SUM, cpugrid);
  tot_kin_energy                 = real_tmp;
  MPI_Allreduce( &Ekin_left,      &real_tmp, 1, REAL, MPI_SUM, cpugrid);
  Ekin_left                      = real_tmp;
  MPI_Allreduce( &Ekin_right,     &real_tmp, 1, REAL, MPI_SUM, cpugrid);
  Ekin_right                     = real_tmp;
  MPI_Allreduce( &inv_mass_left,  &real_tmp, 1, REAL, MPI_SUM, cpugrid);
  inv_mass_left                  = real_tmp;
  MPI_Allreduce( &inv_mass_right, &real_tmp, 1, REAL, MPI_SUM, cpugrid);
  inv_mass_right                 = real_tmp;
  MPI_Allreduce( &tot_impuls_left,&vectmp, DIM, REAL, MPI_SUM, cpugrid);
  tot_impuls_left                = vectmp;
  MPI_Allreduce(&tot_impuls_right,&vectmp, DIM, REAL, MPI_SUM, cpugrid);
  tot_impuls_right               = vectmp;
  MPI_Allreduce( &natoms_left,    &int_tmp,  1, MPI_INT,  MPI_SUM, cpugrid);
  natoms_left                    = int_tmp;
  MPI_Allreduce( &natoms_right,   &int_tmp,  1, MPI_INT,  MPI_SUM, cpugrid);
  natoms_right                   = int_tmp;
#endif

  inv_mass_left      /= 2.0;
  inv_mass_right     /= 2.0;
  tot_impuls_left.x  /= natoms_left;
  tot_impuls_right.x /= natoms_right;
  tot_impuls_left.y  /= natoms_left;
  tot_impuls_right.y /= natoms_right;
#ifndef TWOD
  tot_impuls_left.z  /= natoms_left;
  tot_impuls_right.z /= natoms_right;
#endif

  /* rescale factors for momenta */
  real_tmp = Ekin_left
             - inv_mass_left * SPROD(tot_impuls_left,tot_impuls_left);
  rescale = sqrt( DIM * tran_Tleft * natoms_left / real_tmp  );
  real_tmp = Ekin_right
             - inv_mass_right * SPROD(tot_impuls_right,tot_impuls_right);
  Rescale = sqrt( DIM * tran_Tright * natoms_right / real_tmp  );

  for (k=0; k<ncells; ++k) {

    int i;
    cell *p;
    p = cell_array + CELLS(k);

    for (i=0; i<p->n; ++i) {
	    
      /* which layer? */
      num = scale * p->ort X(i);
      if (num < 0)             num = 0;
      if (num >= tran_nlayers) num = tran_nlayers-1;

      /* rescale momenta */
      if (num == 0) {
        p->impuls X(i) = (p->impuls X(i)-tot_impuls_left.x)*rescale;
        p->impuls Y(i) = (p->impuls Y(i)-tot_impuls_left.y)*rescale;
#ifndef TWOD
        p->impuls Z(i) = (p->impuls Z(i)-tot_impuls_left.z)*rescale;
#endif
      } else if (num == nhalf) {
        p->impuls X(i) = (p->impuls X(i)-tot_impuls_right.x)*Rescale;
        p->impuls Y(i) = (p->impuls Y(i)-tot_impuls_right.y)*Rescale;
#ifndef TWOD
        p->impuls Z(i) = (p->impuls Z(i)-tot_impuls_right.z)*Rescale;
#endif
      }

#ifdef STRESS_TENS
        p->presstens[i].xx += p->impuls X(i) * p->impuls X(i)/MASSE(p,i);
        p->presstens[i].yy += p->impuls Y(i) * p->impuls Y(i)/MASSE(p,i);
#ifndef TWOD
        p->presstens[i].zz += p->impuls Z(i) * p->impuls Z(i)/MASSE(p,i);
        p->presstens[i].yz += p->impuls Y(i) * p->impuls Z(i)/MASSE(p,i);
        p->presstens[i].zx += p->impuls Z(i) * p->impuls X(i)/MASSE(p,i);
#endif
        p->presstens[i].xy += p->impuls X(i) * p->impuls Y(i)/MASSE(p,i);
#endif
    }
  }
#ifdef RNEMD
  heat_transfer += tran_Tleft *natoms_left  * DIM/2 - Ekin_left/2;  /* hot  */ 
  heat_transfer -= tran_Tright*natoms_right * DIM/2 - Ekin_right/2; /* cold */
#endif
}

#else

void move_atoms_nvx(void) 
{
  if (myid==0) 
  error("the chosen ensemble NVX is not supported by this binary");
}

#endif

