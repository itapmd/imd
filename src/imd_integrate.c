
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

/*****************************************************************************
*
* Basic NVE Integrator
*
*****************************************************************************/

#if defined(NVE) || defined(EPITAX) || defined(CG)

void move_atoms_nve(void){
  int k;
  real tmpvec1[4], tmpvec2[4];
  static int count = 0;
#ifdef CLONE
  int clones;
#endif

#if defined(NVE) || !defined(EPITAX)
  tot_kin_energy = 0.0;
#endif
  fnorm   = 0.0;
  PxF     = 0.0;
  omega_E = 0.0;

  /* loop over all cells */
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tot_kin_energy,fnorm,omega_E,PxF)
#endif

  for (k=0; k<ncells; ++k) { /* loop over all cells */

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

    
    for (i=0; i<p->n; ++i) { /* loop over all atoms in the cell */


#if defined(EPITAX) && !defined(NVE) 
        /* beam atoms are always integrated by NVE */
        if ( (NUMMER(p,i) <= epitax_sub_n) && 
             (POTENG(p,i) <= epitax_ctrl * epitax_poteng_min) ) continue;
#endif

        kin_energie_1 = SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) );
#ifdef UNIAX
        rot_energie_1 = SPRODN( &DREH_IMPULS(p,i,X), &DREH_IMPULS(p,i,X) );
#endif

        sort = VSORTE(p,i);
#ifdef FBC
        /* give virtual particles their extra force */
	KRAFT(p,i,X) += (fbc_forces + sort)->x;
	KRAFT(p,i,Y) += (fbc_forces + sort)->y;
#ifndef TWOD
	KRAFT(p,i,Z) += (fbc_forces + sort)->z;
#endif
#endif

	/* and set their force (->momentum) in restricted directions to 0 */
	KRAFT(p,i,X) *= (restrictions + sort)->x;
	KRAFT(p,i,Y) *= (restrictions + sort)->y;
#ifndef TWOD
	KRAFT(p,i,Z) *= (restrictions + sort)->z;
#endif

#ifdef FNORM
	fnorm   += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) );
#endif
#ifdef EINSTEIN
	omega_E += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) ) / MASSE(p,i);
#endif

	IMPULS(p,i,X) += timestep * KRAFT(p,i,X);
        IMPULS(p,i,Y) += timestep * KRAFT(p,i,Y);
#ifndef TWOD
        IMPULS(p,i,Z) += timestep * KRAFT(p,i,Z);
#endif

	/* "Globale Konvergenz": like mik, just with the global 
           force and momentum vectors */
#ifdef GLOK              
	PxF += SPRODN( &IMPULS(p,i,X), &KRAFT(p,i,X) );
#endif

#ifdef UNIAX
        dot = 2.0 * SPRODN( &DREH_IMPULS(p,i,X), &ACHSE(p,i,X) );
        DREH_IMPULS(p,i,X) += timestep * DREH_MOMENT(p,i,X)
                               - dot * ACHSE(p,i,X);
        DREH_IMPULS(p,i,Y) += timestep * DREH_MOMENT(p,i,Y)
                               - dot * ACHSE(p,i,Y);
        DREH_IMPULS(p,i,Z) += timestep * DREH_MOMENT(p,i,Z)
                               - dot * ACHSE(p,i,Z);
#endif

        kin_energie_2 = SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) );
#ifdef UNIAX
        rot_energie_2 = SPRODN( &DREH_IMPULS(p,i,X), &DREH_IMPULS(p,i,X) ); 
#endif
        /* sum up kinetic energy on this CPU */
        tot_kin_energy += (kin_energie_1 + kin_energie_2) / (4 * MASSE(p,i));
#ifdef UNIAX
        tot_kin_energy += (rot_energie_1 + rot_energie_2) / 
                                                    (4 * TRAEG_MOMENT(p,i) );
#endif	  

        /* new positions */
        tmp = timestep / MASSE(p,i);
        ORT(p,i,X) += tmp * IMPULS(p,i,X);
        ORT(p,i,Y) += tmp * IMPULS(p,i,Y);
#ifndef TWOD
        ORT(p,i,Z) += tmp * IMPULS(p,i,Z);
#endif

#ifdef SHOCK
	if (shock_mode == 3) {
	  if (ORT(p,i,X) > box_x.x) IMPULS(p,i,X) = -IMPULS(p,i,X);
	}
#endif

        /* new molecular axes */
#ifdef UNIAX
        cross.x = DREH_IMPULS(p,i,Y) * ACHSE(p,i,Z)
                - DREH_IMPULS(p,i,Z) * ACHSE(p,i,Y);
        cross.y = DREH_IMPULS(p,i,Z) * ACHSE(p,i,X)
                - DREH_IMPULS(p,i,X) * ACHSE(p,i,Z);
        cross.z = DREH_IMPULS(p,i,X) * ACHSE(p,i,Y)
                - DREH_IMPULS(p,i,Y) * ACHSE(p,i,X);

        ACHSE(p,i,X) += timestep * cross.x / TRAEG_MOMENT(p,i);
        ACHSE(p,i,Y) += timestep * cross.y / TRAEG_MOMENT(p,i);
        ACHSE(p,i,Z) += timestep * cross.z / TRAEG_MOMENT(p,i);

        norm = sqrt( SPRODN( &ACHSE(p,i,X), &ACHSE(p,i,X) ) );
	    
        ACHSE(p,i,X) /= norm;
        ACHSE(p,i,Y) /= norm;
        ACHSE(p,i,Z) /= norm;
#endif    

#ifdef STRESS_TENS
#ifdef SHOCK
	/* plate against bulk */
        if (shock_mode == 1) {
          if ( ORT(p,i,X) < shock_strip )
            PRESSTENS(p,i,xx) += (IMPULS(p,i,X) - shock_speed * MASSE(p,i)) 
                    * (IMPULS(p,i,X) - shock_speed * MASSE(p,i)) / MASSE(p,i);
          else
	    PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
        }
	/* two halves against one another */
        if (shock_mode == 2) {
          if ( ORT(p,i,X) < box_x.x*0.5 )
            PRESSTENS(p,i,xx) += (IMPULS(p,i,X) - shock_speed * MASSE(p,i)) 
                    * (IMPULS(p,i,X) - shock_speed * MASSE(p,i)) / MASSE(p,i);
          else
            PRESSTENS(p,i,xx) += (IMPULS(p,i,X) + shock_speed * MASSE(p,i)) 
                    * (IMPULS(p,i,X) + shock_speed * MASSE(p,i)) / MASSE(p,i);
        }
	/* bulk against wall */
        if (shock_mode == 3) PRESSTENS(p,i,xx) += (IMPULS(p,i,X) - 
	  shock_speed * MASSE(p,i)) * (IMPULS(p,i,X) - shock_speed * 
           MASSE(p,i)) / MASSE(p,i);
#else
        PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
#endif
        PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
        PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
        PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
#endif /* STRESS_TENS */
      	
#if defined(EPITAX) && !defined(NVE)
    }   
#endif
    
#ifdef CLONE
	for(clones=1;clones<=nclones;clones++)
	    { 	/* and now do the same to all the clones */

		kin_energie_1 = SPRODN( &IMPULS(p,i+clones,X), &IMPULS(p,i+clones,X) );

		/* for everything which has to do with summing up something 
		   we use the original atom */
#ifdef FNORM
		fnorm   += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) );
#endif
#ifdef EINSTEIN
		omega_E += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) ) / MASSE(p,i);
#endif

		/* needed for kin_energy_1 */
		IMPULS(p,i+clones,X) = IMPULS(p,i,X);
		IMPULS(p,i+clones,Y) = IMPULS(p,i,Y);
		IMPULS(p,i+clones,Z) = IMPULS(p,i,Z);

#ifdef GLOK              
		PxF += SPRODN( &IMPULS(p,i,X), &KRAFT(p,i,X) );
#endif

		kin_energie_2 = SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) );
		tot_kin_energy += (kin_energie_1 + kin_energie_2) / (4 * MASSE(p,i));

		/* for new positions and other properties 
		   of the clones use the impuls of original atom */
		tmp = timestep / MASSE(p,i);
		ORT(p,i+clones,X) += tmp * IMPULS(p,i,X);
		ORT(p,i+clones,Y) += tmp * IMPULS(p,i,Y);
#ifndef TWOD
		ORT(p,i+clones,Z) += tmp * IMPULS(p,i,Z);
#endif

	       
#ifdef STRESS_TENS
		PRESSTENS(p,i+clones,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
		PRESSTENS(p,i+clones,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
		PRESSTENS(p,i+clones,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
		PRESSTENS(p,i+clones,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
		PRESSTENS(p,i+clones,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
		PRESSTENS(p,i+clones,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
#endif 
	    }
    i+=nclones;
  

#endif

    }
  }

#ifdef MPI
  /* add up results from different CPUs */
  tmpvec1[0] = tot_kin_energy;
  tmpvec1[1] = fnorm;
  tmpvec1[2] = PxF;
  tmpvec1[3] = omega_E;

  MPI_Allreduce( tmpvec1, tmpvec2, 4, REAL, MPI_SUM, cpugrid);

  tot_kin_energy = tmpvec2[0];
  fnorm          = tmpvec2[1];
  PxF            = tmpvec2[2];
  omega_E        = tmpvec2[3];
#endif
  
#ifdef AND
  /* Andersen Thermostat -- Initialize the velocities now and then */
  ++count;
  if ((tempintv!=0) && (0==count%tempintv)) maxwell(temperature);
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
  fnorm   = 0.0;

#ifdef CLONE
  int clones;
#endif

#ifdef AND
  /* Andersen Thermostat -- Initialize the velocities now and then */
  ++count;
  if ((tempintv!=0) && (0==count%tempintv)) maxwell(temperature);
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
        if ( (NUMMER(p,i) > epitax_sub_n) && 
             (POTENG(p,i) > epitax_ctrl * epitax_poteng_min) ) continue;
#endif

        kin_energie_1 = SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) );

	sort = VSORTE(p,i);
#ifdef FBC
        /* give virtual particles their extra force */
	KRAFT(p,i,X) += (fbc_forces + sort)->x;
	KRAFT(p,i,Y) += (fbc_forces + sort)->y;
#ifndef TWOD
	KRAFT(p,i,Z) += (fbc_forces + sort)->z;
#endif
#endif /* FBC */

	/* and set their force (->momentum) in restricted directions to 0 */
	KRAFT(p,i,X) *= (restrictions + sort)->x;
	KRAFT(p,i,Y) *= (restrictions + sort)->y;
#ifndef TWOD
	KRAFT(p,i,Z) *= (restrictions + sort)->z;
#endif
	
#ifdef FNORM
	fnorm += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) );
#endif

        IMPULS(p,i,X) += timestep * KRAFT(p,i,X);
        IMPULS(p,i,Y) += timestep * KRAFT(p,i,Y);
#ifndef TWOD
        IMPULS(p,i,Z) += timestep * KRAFT(p,i,Z);
#endif

	/* Mikroconvergence Algorithm - set velocity zero if a*v < 0 */
	if (0.0 > SPRODN( &IMPULS(p,i,X), &KRAFT(p,i,X) ) ) {
          IMPULS(p,i,X) = 0.0;
          IMPULS(p,i,Y) = 0.0;
#ifndef TWOD
          IMPULS(p,i,Z) = 0.0;
#endif
        } else { /* new positions */
          tmp = timestep / MASSE(p,i);
          ORT(p,i,X) += tmp * IMPULS(p,i,X);
          ORT(p,i,Y) += tmp * IMPULS(p,i,Y);
#ifndef TWOD
          ORT(p,i,Z) += tmp * IMPULS(p,i,Z);
#endif
        }
        kin_energie_2 = SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) );

        /* sum up kinetic energy on this CPU */ 
        tot_kin_energy += (kin_energie_1 + kin_energie_2) / (4.0 * MASSE(p,i));

#ifdef STRESS_TENS
        PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
        PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
        PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
        PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
#endif

#ifdef CLONE
	for(clones=1;clones<=nclones;clones++)
	    { 	/* and now do the same to all the clones */

        kin_energie_1 = SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) );


#ifdef FNORM
	fnorm += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) );
#endif

        IMPULS(p,i+clones,X)  =       IMPULS(p,i,X) ;
        IMPULS(p,i+clones,Y)  =       IMPULS(p,i,Y) ;
        IMPULS(p,i+clones,Z)  =       IMPULS(p,i,Z) ;

	tmp = timestep / MASSE(p,i);
	ORT(p,i+clones,X) += tmp * IMPULS(p,i,X);
	ORT(p,i+clones,Y) += tmp * IMPULS(p,i,Y);
	ORT(p,i+clones,Z) += tmp * IMPULS(p,i,Z);
        
        kin_energie_2 = SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) );

        /* sum up kinetic energy on this CPU */ 
        tot_kin_energy += (kin_energie_1 + kin_energie_2) / (4.0 * MASSE(p,i));

#ifdef STRESS_TENS
        PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
        PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
        PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
        PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
#endif

	    }
    i+=nclones;
#endif
    }
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
  real tmpvec1[5], tmpvec2[5], ttt;
  real E_kin_1 = 0.0, E_kin_2 = 0.0;
  real reibung, eins_d_reib;
  real E_rot_1 = 0.0, E_rot_2 = 0.0;
#ifdef UNIAX
  real reibung_rot,  eins_d_reib_rot;
#endif
#ifdef CLONE
  int clones;
#endif
  fnorm   = 0.0;
  omega_E = 0.0;

  reibung         =        1.0 - eta * timestep / 2.0;
  eins_d_reib     = 1.0 / (1.0 + eta * timestep / 2.0);
#ifdef UNIAX
  reibung_rot     =        1.0 - eta_rot * timestep / 2.0;
  eins_d_reib_rot = 1.0 / (1.0 + eta_rot * timestep / 2.0);
#endif
   
#ifdef _OPENMP
#pragma omp parallel for reduction(+:E_kin_1,E_kin_2,E_rot_1,E_rot_2,fnorm,omega_E)
#endif
  for (k=0; k<ncells; ++k) {  /* loop over cells */

    int i;
    int sort;
    cell *p;
    real tmp;
#ifdef UNIAX
    real dot, norm ;
    vektor cross ;
#endif
    p = cell_array + CELLS(k);

    for (i=0; i<p->n; ++i) {  /* loop over atoms */

#ifdef EPITAX
        /* only substrate atoms are integrated by NVT */
        if ( (NUMMER(p,i) > epitax_sub_n) && 
             (POTENG(p,i) > epitax_ctrl * epitax_poteng_min) ) continue;
#endif

        /* twice the old kinetic energy */
        E_kin_1 += SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / MASSE(p,i);
#ifdef UNIAX
        E_rot_1 += SPRODN( &DREH_IMPULS(p,i,X), &DREH_IMPULS(p,i,X) ) / 
                                                       TRAEG_MOMENT(p,i);
#endif

	sort = VSORTE(p,i);
#ifdef FBC
        /* give virtual particles their extra force */
	KRAFT(p,i,X) += (fbc_forces + sort)->x;
	KRAFT(p,i,Y) += (fbc_forces + sort)->y;
#ifndef TWOD
	KRAFT(p,i,Z) += (fbc_forces + sort)->z;
#endif

#endif
	KRAFT(p,i,X) *= (restrictions + sort)->x;
	KRAFT(p,i,Y) *= (restrictions + sort)->y;
#ifndef TWOD
	KRAFT(p,i,Z) *= (restrictions + sort)->z;
#endif
#ifdef FNORM
	fnorm   += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) );
#endif
#ifdef EINSTEIN
	omega_E += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) ) / MASSE(p,i);
#endif

	IMPULS(p,i,X) = (IMPULS(p,i,X) * reibung + timestep * KRAFT(p,i,X)) 
                           * eins_d_reib * (restrictions + sort)->x;
        IMPULS(p,i,Y) = (IMPULS(p,i,Y) * reibung + timestep * KRAFT(p,i,Y)) 
                           * eins_d_reib * (restrictions + sort)->y;
#ifndef TWOD
        IMPULS(p,i,Z) = (IMPULS(p,i,Z) * reibung + timestep * KRAFT(p,i,Z)) 
                           * eins_d_reib * (restrictions + sort)->z;
#endif

#ifdef UNIAX
        /* new angular momenta */
        dot = 2.0 * SPRODN( &DREH_IMPULS(p,i,X), &ACHSE(p,i,X) );

        DREH_IMPULS(p,i,X) = eins_d_reib_rot
            * ( DREH_IMPULS(p,i,X) * reibung_rot
                + timestep * DREH_MOMENT(p,i,X) - dot * ACHSE(p,i,X) );
        DREH_IMPULS(p,i,Y) = eins_d_reib_rot
            * ( DREH_IMPULS(p,i,Y) * reibung_rot
                + timestep * DREH_MOMENT(p,i,Y) - dot * ACHSE(p,i,Y) );
        DREH_IMPULS(p,i,Z) = eins_d_reib_rot
            * ( DREH_IMPULS(p,i,Z) * reibung_rot
                + timestep * DREH_MOMENT(p,i,Z) - dot * ACHSE(p,i,Z) );
#endif

        /* twice the new kinetic energy */ 
        E_kin_2 += SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / MASSE(p,i);
#ifdef UNIAX
        E_rot_2 += SPRODN( &DREH_IMPULS(p,i,X), &DREH_IMPULS(p,i,X) ) / 
                                                       TRAEG_MOMENT(p,i);
#endif

        /* new positions */
        tmp = timestep / MASSE(p,i);
        ORT(p,i,X) += tmp * IMPULS(p,i,X);
        ORT(p,i,Y) += tmp * IMPULS(p,i,Y);
#ifndef TWOD
        ORT(p,i,Z) += tmp * IMPULS(p,i,Z);
#endif

#ifdef UNIAX
        cross.x = DREH_IMPULS(p,i,Y) * ACHSE(p,i,Z)
                - DREH_IMPULS(p,i,Z) * ACHSE(p,i,Y);
        cross.y = DREH_IMPULS(p,i,Z) * ACHSE(p,i,X)
                - DREH_IMPULS(p,i,X) * ACHSE(p,i,Z);
        cross.z = DREH_IMPULS(p,i,X) * ACHSE(p,i,Y)
                - DREH_IMPULS(p,i,Y) * ACHSE(p,i,X);

        ACHSE(p,i,X) += timestep * cross.x / TRAEG_MOMENT(p,i);
        ACHSE(p,i,Y) += timestep * cross.y / TRAEG_MOMENT(p,i);
        ACHSE(p,i,Z) += timestep * cross.z / TRAEG_MOMENT(p,i);

        norm = sqrt( SPRODN( &ACHSE(p,i,X), &ACHSE(p,i,X) ));

        ACHSE(p,i,X) /= norm;
        ACHSE(p,i,Y) /= norm;
        ACHSE(p,i,Z) /= norm;
#endif

#ifdef STRESS_TENS
        PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
        PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
        PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
        PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
#endif
    

#ifdef CLONE
	for(clones=1;clones<=nclones;clones++)
	    { 	/* and now do the same to all the clones */

        /* twice the old kinetic energy */
        E_kin_1 += SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / MASSE(p,i);

#ifdef FNORM
	fnorm   += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) );
#endif
#ifdef EINSTEIN
	omega_E += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) ) / MASSE(p,i);
#endif

	IMPULS(p,i+clones,X) = IMPULS(p,i,X);
	IMPULS(p,i+clones,Y) = IMPULS(p,i,Y);
	IMPULS(p,i+clones,Z) = IMPULS(p,i,Z);

        /* twice the new kinetic energy */ 
        E_kin_2 += SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / MASSE(p,i);

        /* new positions */
        tmp = timestep / MASSE(p,i);
        ORT(p,i+clones,X) += tmp * IMPULS(p,i,X);
        ORT(p,i+clones,Y) += tmp * IMPULS(p,i,Y);
        ORT(p,i+clones,Z) += tmp * IMPULS(p,i,Z);

#ifdef STRESS_TENS
        PRESSTENS(p,i+clones,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
        PRESSTENS(p,i+clones,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
        PRESSTENS(p,i+clones,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i+clones,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i+clones,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
        PRESSTENS(p,i+clones,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
#endif
    }    
	i+=nclones;
    
#endif
	   
    } 
}
  
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
  tmpvec1[3] = fnorm;
  tmpvec1[4] = omega_E;

  MPI_Allreduce( tmpvec1, tmpvec2, 5, REAL, MPI_SUM, cpugrid);

  tot_kin_energy = tmpvec2[0];
  E_kin_2        = tmpvec2[1];
  E_rot_2        = tmpvec2[2];
  fnorm          = tmpvec2[3];
  omega_E        = tmpvec2[4];
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


/*****************************************************************************
*
* NVT Integrator with Nose Hoover Thermostat and some shearing (?)
*
*****************************************************************************/

#ifdef SLLOD

void move_atoms_sllod(void)

{
  int k;
  real tmpvec1[5], tmpvec2[5], ttt;
  real E_kin_1 = 0.0, E_kin_2 = 0.0;
  vektor reibung, eins_d_reib;
  real E_rot_1 = 0.0, E_rot_2 = 0.0;
#ifdef UNIAX
  real reibung_rot,  eins_d_reib_rot;
#endif
  fnorm   = 0.0;

#ifdef TWOD
  reibung.x         =        1.0 - (eta+shear_rate.x) * timestep / 2.0;
  eins_d_reib.x     = 1.0 / (1.0 + (eta+shear_rate.x) * timestep / 2.0);
  reibung.y         =        1.0 - (eta+shear_rate.y) * timestep / 2.0;
  eins_d_reib.y     = 1.0 / (1.0 + (eta+shear_rate.y) * timestep / 2.0);
#else
  reibung.x         =        1.0 - (eta+shear_rate.z+shear_rate2.y) * timestep / 2.0;
  eins_d_reib.x     = 1.0 / (1.0 + (eta+shear_rate.z+shear_rate2.y) * timestep / 2.0);
  reibung.y         =        1.0 - (eta+shear_rate.x+shear_rate2.z) * timestep / 2.0;
  eins_d_reib.y     = 1.0 / (1.0 + (eta+shear_rate.x+shear_rate2.z) * timestep / 2.0);
  reibung.z         =        1.0 - (eta+shear_rate.y+shear_rate2.x) * timestep / 2.0;
  eins_d_reib.z     = 1.0 / (1.0 + (eta+shear_rate.y+shear_rate2.x) * timestep / 2.0);
#endif
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

        /* twice the old kinetic energy */
        E_kin_1 += SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / MASSE(p,i);
#ifdef UNIAX
        E_rot_1 += SPRODN( &DREH_IMPULS(p,i,X), &DREH_IMPULS(p,i,X) ) / 
                                                       TRAEG_MOMENT(p,i);
#endif

	sort = VSORTE(p,i);
#ifdef FBC
        /* give virtual particles their extra force */
	KRAFT(p,i,X) += (fbc_forces + sort)->x;
	KRAFT(p,i,Y) += (fbc_forces + sort)->y;
#ifndef TWOD
	KRAFT(p,i,Z) += (fbc_forces + sort)->z;
#endif

#endif
	KRAFT(p,i,X) *= (restrictions + sort)->x;
	KRAFT(p,i,Y) *= (restrictions + sort)->y;
#ifndef TWOD
	KRAFT(p,i,Z) *= (restrictions + sort)->z;
#endif
#ifdef FNORM
	fnorm += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) );
#endif

	IMPULS(p,i,X) = (IMPULS(p,i,X) * reibung.x + timestep * KRAFT(p,i,X)) 
                           * eins_d_reib.x * (restrictions + sort)->x;
        IMPULS(p,i,Y) = (IMPULS(p,i,Y) * reibung.y + timestep * KRAFT(p,i,Y)) 
                           * eins_d_reib.y * (restrictions + sort)->y;
#ifndef TWOD
        IMPULS(p,i,Z) = (IMPULS(p,i,Z) * reibung.z + timestep * KRAFT(p,i,Z)) 
                           * eins_d_reib.z * (restrictions + sort)->z;
#endif

#ifdef UNIAX
        /* new angular momenta */
        dot = 2.0 * SPRODN( &DREH_IMPULS(p,i,X), &ACHSE(p,i,X) );

        DREH_IMPULS(p,i,X) = eins_d_reib_rot
            * ( DREH_IMPULS(p,i,X) * reibung_rot
                + timestep * DREH_MOMENT(p,i,X) - dot * ACHSE(p,i,X) );
        DREH_IMPULS(p,i,Y) = eins_d_reib_rot
            * ( DREH_IMPULS(p,i,Y) * reibung_rot
                + timestep * DREH_MOMENT(p,i,Y) - dot * ACHSE(p,i,Y) );
        DREH_IMPULS(p,i,Z) = eins_d_reib_rot
            * ( DREH_IMPULS(p,i,Z) * reibung_rot
                + timestep * DREH_MOMENT(p,i,Z) - dot * ACHSE(p,i,Z) );
#endif

        /* twice the new kinetic energy */ 
        E_kin_2 += SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / MASSE(p,i);
#ifdef UNIAX
        E_rot_2 += SPRODN( &DREH_IMPULS(p,i,X), &DREH_IMPULS(p,i,X) ) / 
                                                       TRAEG_MOMENT(p,i);
#endif

        /* new positions */
        tmp = timestep / MASSE(p,i);
        ORT(p,i,X) += tmp * IMPULS(p,i,X);
        ORT(p,i,Y) += tmp * IMPULS(p,i,Y);
#ifndef TWOD
        ORT(p,i,Z) += tmp * IMPULS(p,i,Z);
	/* sllod specific */
        ORT(p,i,X) += shear_rate.z  * ORT(p,i,Y);
        ORT(p,i,X) += shear_rate2.y * ORT(p,i,Z);
        ORT(p,i,Y) += shear_rate.x  * ORT(p,i,Z);
        ORT(p,i,Y) += shear_rate2.z * ORT(p,i,X);
        ORT(p,i,Z) += shear_rate.y  * ORT(p,i,X);
        ORT(p,i,Z) += shear_rate2.x * ORT(p,i,Y);
#else
        ORT(p,i,X) += shear_rate.x * ORT(p,i,Y);
        ORT(p,i,Y) += shear_rate.y * ORT(p,i,X);
#endif
#ifdef UNIAX
        cross.x = DREH_IMPULS(p,i,Y) * ACHSE(p,i,Z)
                - DREH_IMPULS(p,i,Z) * ACHSE(p,i,Y);
        cross.y = DREH_IMPULS(p,i,Z) * ACHSE(p,i,X)
                - DREH_IMPULS(p,i,X) * ACHSE(p,i,Z);
        cross.z = DREH_IMPULS(p,i,X) * ACHSE(p,i,Y)
                - DREH_IMPULS(p,i,Y) * ACHSE(p,i,X);

        ACHSE(p,i,X) += timestep * cross.x / TRAEG_MOMENT(p,i);
        ACHSE(p,i,Y) += timestep * cross.y / TRAEG_MOMENT(p,i);
        ACHSE(p,i,Z) += timestep * cross.z / TRAEG_MOMENT(p,i);

        norm = sqrt( SPRODN( &ACHSE(p,i,X), &ACHSE(p,i,X) );

        ACHSE(p,i,X) /= norm;
        ACHSE(p,i,Y) /= norm;
        ACHSE(p,i,Z) /= norm;
#endif

#ifdef STRESS_TENS
        PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
        PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
        PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
        PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
#endif
    }
  }
  
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
  tmpvec1[3] = fnorm;

  MPI_Allreduce( tmpvec1, tmpvec2, 4, REAL, MPI_SUM, cpugrid);

  tot_kin_energy = tmpvec2[0];
  E_kin_2        = tmpvec2[1];
  E_rot_2        = tmpvec2[2];
  fnorm          = tmpvec2[3];
#endif

  /* adjusting the box */
#ifdef TWOD
  box_x.y += shear_rate.y*box_y.y;
  box_y.x += shear_rate.x*box_x.x;
#else
  box_y.x += shear_rate.z  * box_y.y;
  box_z.x += shear_rate2.y * box_z.z;
  box_z.y += shear_rate.x  * box_z.z;
  box_x.y += shear_rate2.z * box_x.x;
  box_x.z += shear_rate.y  * box_x.x;
  box_y.z += shear_rate2.x * box_y.y;
#endif
  make_box();

  /* time evolution of constraints */
  ttt  = nactive * temperature;
  eta += timestep * (E_kin_2 / ttt - 1.0) * isq_tau_eta;
#ifdef UNIAX
  ttt  = nactive_rot * temperature;
  eta_rot += timestep * (E_rot_2 / ttt - 1.0) * isq_tau_eta_rot;
#endif
  
}

#else

void move_atoms_sllod(void) 
{
  if (myid==0)
  error("the chosen ensemble SLLOD is not supported by this binary");
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
      dyn_stress_x += IMPULS(p,i,X) * IMPULS(p,i,X) * tmp;
      dyn_stress_y += IMPULS(p,i,Y) * IMPULS(p,i,Y) * tmp;
#ifndef TWOD
      dyn_stress_z += IMPULS(p,i,Z) * IMPULS(p,i,Z) * tmp;
#endif
#ifdef UNIAX
      Erot_old += SPRODN( &DREH_IMPULS(p,i,X), &DREH_IMPULS(p,i,X) ) /
                                                        TRAEG_MOMENT(p,i);
#endif
    }
  }

  /* twice the kinetic energy */
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
  real tmpvec1[4], tmpvec2[4], ttt;
  real reib, ireib;

  fnorm    = 0.0;
  omega_E  = 0.0;
#ifdef UNIAX
  pressure = (0.6 * (Ekin_old + Erot_old) + virial) / (DIM * volume);
#else
  pressure = (Ekin_old + virial) / (DIM * volume) ;
#endif
#ifdef CLONE
  int clones;
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
#pragma omp parallel for reduction(+:Ekin_new,Erot_new,fnorm,omega_E)
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
      fnorm   += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) );
#endif
#ifdef EINSTEIN
      omega_E += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) ) / MASSE(p,i);
#endif

      /* new momenta */
      IMPULS(p,i,X) = (pfric*IMPULS(p,i,X)+timestep*KRAFT(p,i,X))*pifric;
      IMPULS(p,i,Y) = (pfric*IMPULS(p,i,Y)+timestep*KRAFT(p,i,Y))*pifric;
#ifndef TWOD
      IMPULS(p,i,Z) = (pfric*IMPULS(p,i,Z)+timestep*KRAFT(p,i,Z))*pifric;
#endif

#ifdef UNIAX
      /* new angular momenta */
      dot = 2.0 * SPRODN( &DREH_IMPULS(p,i,X), &ACHSE(p,i,X) );
      DREH_IMPULS(p,i,X) = ireib * ( DREH_IMPULS(p,i,X) * reib 
              + timestep * DREH_MOMENT(p,i,X) - dot * ACHSE(p,i,X) );
      DREH_IMPULS(p,i,Y) = ireib * ( DREH_IMPULS(p,i,Y) * reib 
              + timestep * DREH_MOMENT(p,i,Y) - dot * ACHSE(p,i,Y) );
      DREH_IMPULS(p,i,Z) = ireib * ( DREH_IMPULS(p,i,Z) * reib 
              + timestep * DREH_MOMENT(p,i,Z) - dot * ACHSE(p,i,Z) );
#endif

      /* twice the new kinetic energy */ 
      Ekin_new += SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / MASSE(p,i);
#ifdef UNIAX
      Erot_new += SPRODN( &DREH_IMPULS(p,i,X), &DREH_IMPULS(p,i,X) ) /
                                                      TRAEG_MOMENT(p,i);
#endif

      /* new positions */
      tmp = timestep / MASSE(p,i);
      ORT(p,i,X) = (rfric * ORT(p,i,X) + IMPULS(p,i,X) * tmp) * rifric;
      ORT(p,i,Y) = (rfric * ORT(p,i,Y) + IMPULS(p,i,Y) * tmp) * rifric;
#ifndef TWOD
      ORT(p,i,Z) = (rfric * ORT(p,i,Z) + IMPULS(p,i,Z) * tmp) * rifric;
#endif

#ifdef UNIAX
      /* new molecular axes */
      cross.x = DREH_IMPULS(p,i,Y) * ACHSE(p,i,Z)
              - DREH_IMPULS(p,i,Z) * ACHSE(p,i,Y);
      cross.y = DREH_IMPULS(p,i,Z) * ACHSE(p,i,X)
              - DREH_IMPULS(p,i,X) * ACHSE(p,i,Z);
      cross.z = DREH_IMPULS(p,i,X) * ACHSE(p,i,Y)
              - DREH_IMPULS(p,i,Y) * ACHSE(p,i,X);

      ACHSE(p,i,X) += timestep * cross.x / TRAEG_MOMENT(p,i);
      ACHSE(p,i,Y) += timestep * cross.y / TRAEG_MOMENT(p,i);
      ACHSE(p,i,Z) += timestep * cross.z / TRAEG_MOMENT(p,i);

      norm = sqrt( SPRODN( &ACHSE(p,i,X), &ACHSE(p,i,X) ) );

      ACHSE(p,i,X) /= norm;
      ACHSE(p,i,Y) /= norm;
      ACHSE(p,i,Z) /= norm;
#endif

#ifdef STRESS_TENS
      PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
      PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
        PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
        PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
#endif

#ifdef CLONE
	for(clones=1;clones<=nclones;clones++)
	    { 	/* and now do the same to all the clones */
		
#ifdef FNORM
      fnorm   += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) );
#endif
#ifdef EINSTEIN
      omega_E += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) ) / MASSE(p,i);
#endif

      /* new momenta */
      IMPULS(p,i+clones,X) =       IMPULS(p,i,X);
      IMPULS(p,i+clones,Y) =       IMPULS(p,i,Y);
      IMPULS(p,i+clones,Z) =       IMPULS(p,i,Z);

      /* twice the new kinetic energy */ 
      Ekin_new += SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / MASSE(p,i);

      /* new positions */
      tmp = timestep / MASSE(p,i);
      ORT(p,i+clones,X) = (rfric * ORT(p,i+clones,X) + IMPULS(p,i,X) * tmp) * rifric;
      ORT(p,i+clones,Y) = (rfric * ORT(p,i+clones,Y) + IMPULS(p,i,Y) * tmp) * rifric;
      ORT(p,i+clones,Z) = (rfric * ORT(p,i+clones,Z) + IMPULS(p,i,Z) * tmp) * rifric;

#ifdef STRESS_TENS
        PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
        PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
        PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
        PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
#endif

	    }
	i+=nclones;

#endif
    }
  }

#ifdef MPI
  /* add up results from all CPUs */
  tmpvec1[0] = Ekin_new;
  tmpvec1[1] = Erot_new;
  tmpvec1[2] = fnorm;
  tmpvec1[3] = omega_E;

  MPI_Allreduce( tmpvec1, tmpvec2, 4, REAL, MPI_SUM, cpugrid);

  Ekin_new = tmpvec2[0];
  Erot_new = tmpvec2[1];
  fnorm    = tmpvec2[2];
  omega_E  = tmpvec2[3];
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
  real Ekin_new = 0.0, ttt, tmpvec1[6], tmpvec2[6];
  vektor pfric, pifric, rfric, rifric, tvec;

  fnorm    = 0.0;
  omega_E  = 0.0;
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
#pragma omp parallel for reduction(+:Ekin,dyn_stress_x,dyn_stress_y,dyn_stress_z,fnorm,omega_E)
#endif
  for (k=0; k<ncells; ++k) {

    int i;
    cell *p;
    real tmp;
    p = cell_array + CELLS(k);

    /* loop over atoms in cell */
    for (i=0; i<p->n; ++i) {

      tmp = 1.0 / MASSE(p,i);
#ifdef FNORM
      fnorm   += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) );
#endif
#ifdef EINSTEIN
      omega_E += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) ) * tmp;
#endif
#ifdef STRESS_TENS
      PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) * tmp;
      PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) * tmp;
#ifndef TWOD
      PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) * tmp;
      PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) * tmp;
      PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) * tmp;
#endif
      PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) * tmp;
#endif

      /* new momenta */
      IMPULS(p,i,X) = (pfric.x * IMPULS(p,i,X)
                                   + timestep * KRAFT(p,i,X)) * pifric.x;
      IMPULS(p,i,Y) = (pfric.y * IMPULS(p,i,Y)
                                   + timestep * KRAFT(p,i,Y)) * pifric.y;
#ifndef TWOD
      IMPULS(p,i,Z) = (pfric.z * IMPULS(p,i,Z)
                                   + timestep * KRAFT(p,i,Z)) * pifric.z;
#endif

      /* new stress tensor (dynamic part only) */
      dyn_stress_x += IMPULS(p,i,X) * IMPULS(p,i,X) * tmp;
      dyn_stress_y += IMPULS(p,i,Y) * IMPULS(p,i,Y) * tmp;
#ifndef TWOD
      dyn_stress_z += IMPULS(p,i,Z) * IMPULS(p,i,Z) * tmp;
#endif

      /* twice the new kinetic energy */ 
      Ekin_new += SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) * tmp;
	  
      /* new positions */
      tmp *= timestep;
      ORT(p,i,X) = (rfric.x * ORT(p,i,X) + IMPULS(p,i,X) * tmp) * rifric.x;
      ORT(p,i,Y) = (rfric.y * ORT(p,i,Y) + IMPULS(p,i,Y) * tmp) * rifric.y;
#ifndef TWOD
      ORT(p,i,Z) = (rfric.z * ORT(p,i,Z) + IMPULS(p,i,Z) * tmp) * rifric.z;
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
  tmpvec1[5]   = omega_E;
  MPI_Allreduce( tmpvec1, tmpvec2, 6, REAL, MPI_SUM, cpugrid);
  Ekin_new     = tmpvec2[0];
  fnorm        = tmpvec2[1];
  dyn_stress_x = tmpvec2[2];
  dyn_stress_y = tmpvec2[3];
  dyn_stress_z = tmpvec2[4];
  omega_E      = tmpvec2[5];
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
  real tmpvec1[7], tmpvec2[7], ttt;
  real E_kin_1        = 0.0, E_kin_2        = 0.0; 
  real E_kin_damp1    = 0.0, E_kin_damp2    = 0.0;
  real E_kin_stadium1 = 0.0, E_kin_stadium2 = 0.0;
  real reibung, reibung_y, eins_d_reib, eins_d_reib_y;
  real epsilontmp, eins_d_epsilontmp;

  real tmp, f; /* stadium function: the bath tub !!!!*/

  fnorm     = 0.0;
  sum_f     = 0.0;
  n_stadium = 0;

#ifdef CLONE
  int clones;
#endif

  if(expansionmode==1)
      dotepsilon = dotepsilon0 / (1.0 + dotepsilon0 * steps * timestep);
      

  /* loop over all atoms */
#ifdef _OPENMP
#pragma omp parallel for reduction(+:E_kin_1,E_kin_2,E_kin_damp1,E_kin_damp2,E_kin_stadium1,E_kin_stadium2,sum_f,n_stadium,fnorm)
#endif
  for (k=0; k<ncells; ++k){ 

    int i;
    int sort;
    cell *p;
    real tmp,tmp1,tmp2;

    p = cell_array + CELLS(k);

    for (i=0; i<p->n; ++i){
	
	/* if half axis in x-direction is zero: global viscous damping ! */
	if(stadium.x <= 0.0){ 
	    f = 1.0; 
	} else {
	    /* Calculate stadium function f */
	    tmp1 = SQR((ORT(p,i,X)-center.x)/(2.0*stadium2.x));
	    tmp2 = SQR((ORT(p,i,Y)-center.y)/(2.0*stadium2.y));
	    f    = (tmp1+tmp2-SQR(stadium.x/(2.0*stadium2.x)))/\
		(.25- SQR(stadium.x/(2.0*stadium2.x)));
	
	}
	
	if (f<= 0.0) {
	    f = 0.0;
	    n_stadium += DIM;
	    /* what about the restrictions?? */
#ifdef CLONE
	    n_stadium += nclones*DIM; /* also the clones aren't in f */
#endif
	}
	if (f>1.0) f = 1.0;

       /* we smooth the stadium function: to get a real bath tub !*/  
	f    = .5 * (1 + sin(-M_PI/2.0 + M_PI*f));

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
        tmp = SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / MASSE(p,i);
	E_kin_1 += tmp;
	if (f == 0.0) E_kin_stadium1 +=      tmp;
	if (f >  0.0) E_kin_damp1    +=  f * tmp;
	
#ifdef FBC
        /* give virtual particles their extra force */
	KRAFT(p,i,X) += (fbc_forces + sort)->x;
	KRAFT(p,i,Y) += (fbc_forces + sort)->y;
#ifndef TWOD
	KRAFT(p,i,Z) += (fbc_forces + sort)->z;
#endif
#endif

	KRAFT(p,i,X) *= (restrictions + sort)->x;
	KRAFT(p,i,Y) *= (restrictions + sort)->y;
#ifndef TWOD
	KRAFT(p,i,Z) *= (restrictions + sort)->z;
#endif

#ifdef FNORM
	fnorm   += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) );
#endif

	reibung       =        1.0 -  gamma_damp * f * timestep / 2.0;
	eins_d_reib   = 1.0 / (1.0 +  gamma_damp * f * timestep / 2.0);
	reibung_y     =        1.0 - (gamma_damp * f + dotepsilon) * 
	                        timestep / 2.0;
	eins_d_reib_y = 1.0 / (1.0 + (gamma_damp * f + dotepsilon) * 
			        timestep / 2.0);
	
        /* new momenta */
	IMPULS(p,i,X) = (IMPULS(p,i,X)  * reibung   + timestep * KRAFT(p,i,X))
                           * eins_d_reib   * (restrictions + sort)->x;
        IMPULS(p,i,Y) = (IMPULS(p,i,Y)  * reibung_y + timestep * KRAFT(p,i,Y))
                           * eins_d_reib_y * (restrictions + sort)->y;
#ifndef TWOD
        IMPULS(p,i,Z) = (IMPULS(p,i,Z)  * reibung   + timestep * KRAFT(p,i,Z))
                           * eins_d_reib   * (restrictions + sort)->z;
#endif                  

	/* twice the new kinetic energy */
        tmp = SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / MASSE(p,i);
	E_kin_2 += tmp;
	if (f == 0.0) E_kin_stadium2 +=      tmp;
	if (f >  0.0) E_kin_damp2    +=  f * tmp;


	/* new positions */
        tmp = timestep / MASSE(p,i);
	epsilontmp =               1.0 + dotepsilon * timestep / 2.0;
	eins_d_epsilontmp = 1.0 / (1.0 - dotepsilon * timestep / 2.0);

        ORT(p,i,X) +=  tmp * IMPULS(p,i,X);
        ORT(p,i,Y)  = (tmp * IMPULS(p,i,Y) + epsilontmp * ORT(p,i,Y))
	                * eins_d_epsilontmp;

#ifndef TWOD
        ORT(p,i,Z) +=  tmp * IMPULS(p,i,Z);
#endif

#ifdef STRESS_TENS
        PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
        PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
        PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
        PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
#endif

#ifdef CLONE
	for(clones=1;clones<=nclones;clones++)
	    { 	/* and now do the same to all the clones */
	
		/* f and sort stays the same */

        /* add up f considering the restriction vector  */

	sum_f+= f * ( (restrictions + sort)->x + 
		      (restrictions + sort)->y +  
		      (restrictions + sort)->z  )/3.0;
	
       	/* twice the old kinetic energy */
	/* this has to be redone here, because we need 
		   to work also with the old impuls */
        tmp = SPRODN( &IMPULS(p,i+clones,X), &IMPULS(p,i+clones,X) ) / MASSE(p,i);
	E_kin_1 += tmp;
	if (f == 0.0) E_kin_stadium1 +=      tmp;
	if (f >  0.0) E_kin_damp1    +=  f * tmp;


#ifdef FNORM
	fnorm   += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) );
#endif

        /* new momenta */
	IMPULS(p,i+clones,X) = IMPULS(p,i,X);
	IMPULS(p,i+clones,Y) = IMPULS(p,i,Y);
	IMPULS(p,i+clones,Z) = IMPULS(p,i,Z);

	/* twice the new kinetic energy */
        tmp = SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / MASSE(p,i);
	E_kin_2 += tmp;
	if (f == 0.0) E_kin_stadium2 +=      tmp;
	if (f >  0.0) E_kin_damp2    +=  f * tmp;


	/* new positions */
        tmp = timestep / MASSE(p,i);
	/* allready known */
/* 	epsilontmp =               1.0 + dotepsilon * timestep / 2.0; */
/* 	eins_d_epsilontmp = 1.0 / (1.0 - dotepsilon * timestep / 2.0); */

        ORT(p,i+clones,X) +=  tmp * IMPULS(p,i,X);
        ORT(p,i+clones,Y)  = (tmp * IMPULS(p,i,Y) + epsilontmp * ORT(p,i,Y))
	                * eins_d_epsilontmp;
	ORT(p,i+clones,Z) +=  tmp * IMPULS(p,i,Z);

#ifdef STRESS_TENS
        PRESSTENS(p,i+clones,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
        PRESSTENS(p,i+clones,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
        PRESSTENS(p,i+clones,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i+clones,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i+clones,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
        PRESSTENS(p,i+clones,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
#endif

	

	    }
    i+=nclones;
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

  ttt   = DIM * temperature * sum_f;

  /* time evolution of constraints */
  /* dampingmode: 0 -> viscous damping (default); 
                  1 -> Nose-Hoover; */

  if(dampingmode == 1){
      gamma_damp += timestep * (E_kin_damp2 / ttt - 1.0) * gamma_bar;
  } else {
      
      if( E_kin_damp2 != 0.0){
	  gamma_damp  =        (1.0 - ttt / E_kin_damp2) * gamma_bar;
      } else {
	  gamma_damp  =  0.0;
      }

  }
      
}

#else

void move_atoms_frac(void) 
{
  if (myid==0)
  error("the chosen ensemble FRAC is not supported by this binary");
}

#endif

/*****************************************************************************
*
*  Integrator for fracture studies with temperature gradient 
*  
*
*****************************************************************************/

#ifdef FTG

void move_atoms_ftg(void)

{
  int j, k;
  static real *E_kin_1     = NULL; static real *E_kin_2     = NULL;
  static real *ftgtmpvec1  = NULL; static real *ftgtmpvec2  = NULL;
  static int  *iftgtmpvec1 = NULL; static int  *iftgtmpvec2 = NULL;

  real tmp,tmp1,tmp2;
  real ttt;
  real reibung, reibung_y, eins_d_reib, eins_d_reib_y;
  real epsilontmp, eins_d_epsilontmp;
  int slice;
  real gamma_tmp;

#ifdef CLONE
  int clones;
#endif

  /* alloc vector versions of E_kin and  ftgtmpvect*/
  if (NULL==E_kin_1) {
    E_kin_1=(real*) malloc(nslices*sizeof(real));
    if (NULL==E_kin_1) 
      error("Cannot allocate memory for E_kin_1 vector\n");
  }
  if (NULL==E_kin_2) {
    E_kin_2=(real*) malloc(nslices*sizeof(real));
    if (NULL==E_kin_2) 
      error("Cannot allocate memory for E_kin_2 vector\n");
  }
  if (NULL==ftgtmpvec1) {
    ftgtmpvec1=(real*) malloc(nslices*sizeof(real));
    if (NULL==ftgtmpvec1) 
      error("Cannot allocate memory for ftgtmpvec1 vector\n");
  }
  if (NULL==ftgtmpvec2) {
    ftgtmpvec2=(real*) malloc(nslices*sizeof(real));
    if (NULL==ftgtmpvec2) 
      error("Cannot allocate memory for ftgtmpvec2 vector\n");
  }
  if (NULL==iftgtmpvec1) {
    iftgtmpvec1=(int*) malloc(nslices*sizeof(int));
    if (NULL==iftgtmpvec1) 
      error("Cannot allocate memory for iftgtmpvec1 vector\n");
  }
  if (NULL==iftgtmpvec2) {
    iftgtmpvec2=(int*) malloc(nslices*sizeof(int));
    if (NULL==iftgtmpvec2) 
      error("Cannot allocate memory for iftgtmpvec2 vector\n");
  }

  for (j=0; j<nslices; j++) {
    *(E_kin_1   +j) = 0.0;
    *(E_kin_2   +j) = 0.0;
    *(ninslice  +j) = 0;
  }

  fnorm     = 0.0;

  if(expansionmode==1)
      dotepsilon = dotepsilon0 / (1.0 + dotepsilon0 * steps * timestep);
      

  /* loop over all atoms */

  for (k=0; k<ncells; ++k) {

    int i;
    int sort;
    cell *p;

    p = cell_array + CELLS(k);

    for (i=0; i<p->n; ++i) {
	
      sort = VSORTE(p,i);

      /* calc slice */
      tmp = ORT(p,i,X)/box_x.x; 
      slice = (int) (nslices *tmp);
      if (slice<0)        slice = 0;
      if (slice>=nslices) slice = nslices-1;;      
     
      /* if half axis in y-direction is given: local viscous damping !!! */
      if(stadium.y != 0.0){ 

	/* calc desired temperature */
	temperature = Tleft + (Tright-Tleft) * (nslices*tmp - nslices_Left)
	  /(nslices - nslices_Left - nslices_Right);
	if (temperature < Tleft ) temperature = Tleft;
	if (temperature > Tright) temperature = Tright;
	
	/* calc kinetic "temperature" for actual atom */
	tmp  = SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / MASSE(p,i);
#ifdef TWOD
	tmp2 = ( (restrictions + sort)->x + 
		 (restrictions + sort)->y   );
#else
	tmp2 = ( (restrictions + sort)->x + 
		 (restrictions + sort)->y +  
		 (restrictions + sort)->z  );
#endif
	if(tmp2!=0) tmp /= (real)tmp2;
	
	/* calc damping factor form position */
	gamma_tmp = (fabs(ORT(p,i,Y)-center.y) - stadium.y)/
	  (stadium2.y-stadium.y);
	if ( gamma_tmp < 0.0)  gamma_tmp = 0.0;
	if ( gamma_tmp > 1.0)  gamma_tmp = 1.0;
	
	/* smooth the gamma_tmp funktion*/
	gamma_tmp = .5 * (1 + sin(-M_PI/2.0 + M_PI*gamma_tmp));
	
	/* to share the code with the non local version we overwrite 
	 the gamma values every timestep */
	*(gamma_ftg+slice)  = (gamma_min + gamma_bar * gamma_tmp) 
	  * (tmp-temperature) 
	  / sqrt(SQR(tmp) + SQR(temperature/delta_ftg));     
      } 

      
      /* add up degrees of freedom  considering restriction vector  */
#ifdef TWOD
      *(ninslice + slice) += ( (restrictions + sort)->x + 
			       (restrictions + sort)->y   );
#else
      *(ninslice + slice) += ( (restrictions + sort)->x + 
			       (restrictions + sort)->y +  
			       (restrictions + sort)->z  );
#endif
      
      /* twice the old kinetic energy */
      *(E_kin_1 + slice) += 
        SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / MASSE(p,i);

#ifdef FBC
        /* give virtual particles their extra force */
	KRAFT(p,i,X) += (fbc_forces + sort)->x;
	KRAFT(p,i,Y) += (fbc_forces + sort)->y;
#ifndef TWOD
	KRAFT(p,i,Z) += (fbc_forces + sort)->z;
#endif
#endif

	KRAFT(p,i,X) *= (restrictions + sort)->x;
	KRAFT(p,i,Y) *= (restrictions + sort)->y;
#ifndef TWOD
	KRAFT(p,i,Z) *= (restrictions + sort)->z;
#endif

#ifdef FNORM
	fnorm   += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) );
#endif

	reibung       =        1.0 -  *(gamma_ftg + slice) * timestep / 2.0;
	eins_d_reib   = 1.0 / (1.0 +  *(gamma_ftg + slice) * timestep / 2.0);
	reibung_y     =        1.0 - (*(gamma_ftg + slice) + dotepsilon) * 
	                        timestep / 2.0;
	eins_d_reib_y = 1.0 / (1.0 + (*(gamma_ftg + slice) + dotepsilon) * 
			        timestep / 2.0);
	
        /* new momenta */
	IMPULS(p,i,X) = (IMPULS(p,i,X)  * reibung   + timestep * KRAFT(p,i,X))
                           * eins_d_reib   * (restrictions + sort)->x;
        IMPULS(p,i,Y) = (IMPULS(p,i,Y)  * reibung_y + timestep * KRAFT(p,i,Y))
                           * eins_d_reib_y * (restrictions + sort)->y;
#ifndef TWOD
        IMPULS(p,i,Z) = (IMPULS(p,i,Z)  * reibung   + timestep * KRAFT(p,i,Z))
                           * eins_d_reib   * (restrictions + sort)->z;
#endif                  

	/* twice the new kinetic energy */ 
	*(E_kin_2 + slice) += 
          SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / MASSE(p,i);
	
	/* new positions */
        tmp = timestep / MASSE(p,i);
	epsilontmp =               1.0 + dotepsilon * timestep / 2.0;
	eins_d_epsilontmp = 1.0 / (1.0 - dotepsilon * timestep / 2.0);

        ORT(p,i,X) +=  tmp * IMPULS(p,i,X);
        ORT(p,i,Y)  = (tmp * IMPULS(p,i,Y) + epsilontmp * ORT(p,i,Y))
	                * eins_d_epsilontmp;

#ifndef TWOD
        ORT(p,i,Z) +=  tmp * IMPULS(p,i,Z);
#endif

#ifdef STRESS_TENS
        PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
        PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
        PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
        PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
#endif

#ifdef CLONE
	for(clones=1;clones<=nclones;clones++)
	    { 	/* and now do the same to all the clones */
	    }
    i+=nclones;
#endif  

    }
  }

  tot_kin_energy = 0.0; 
  for (j=0; j<nslices; j++){
    tot_kin_energy += ( *(E_kin_1 + j) + *(E_kin_2 + j)) / 4.0;
    *(E_kin_ftg+j)  = ( *(E_kin_1 + j) + *(E_kin_2 + j)) / 4.0;
  }

#ifdef DEBUG
  printf("%d: ", myid);
  for (j=0;j<nslices;j++)
    printf("%3.6f ", *(E_kin_2 +j));
  printf("\n");
#endif

#ifdef MPI
  /* add up results from different CPUs */
  
  for (j=0; j<nslices; j++) 
    *(ftgtmpvec1 + j) = *(E_kin_ftg + j);
  MPI_Allreduce( ftgtmpvec1, ftgtmpvec2, nslices, REAL, MPI_SUM, cpugrid);
  for (j=0; j<nslices; j++) 
    *(E_kin_ftg + j) = *(ftgtmpvec2 + j);
  
  for (j=0; j<nslices; j++) 
    *(ftgtmpvec1 + j) = *(E_kin_2 + j);
  MPI_Allreduce( ftgtmpvec1, ftgtmpvec2, nslices, REAL, MPI_SUM, cpugrid);
  for (j=0; j<nslices; j++) 
    *(E_kin_2 + j) = *(ftgtmpvec2 + j);
  
  tmp1 = tot_kin_energy;
  MPI_Allreduce( &tmp1, &tmp2, 1, REAL, MPI_SUM, cpugrid);
  tot_kin_energy = tmp2;
  
  for (j=0; j<nslices; j++) 
    *(iftgtmpvec1 +j) = *(ninslice   + j);
  MPI_Allreduce( iftgtmpvec1, iftgtmpvec2, nslices, MPI_INT, MPI_SUM, cpugrid);
  for (j=0; j<nslices; j++) 
    *(ninslice + j) = *(iftgtmpvec2 +j); 
#endif


  for (j=0; j<nslices; j++) {
    temperature =  Tleft + (Tright-Tleft)*(j-nslices_Left+1) /
      (real) (nslices-nslices_Left-nslices_Right+1);
    
    if(j>=nslices-nslices_Right)  temperature = Tright;
    if(j<nslices_Left)            temperature = Tleft;
    
    ttt   = temperature * *(ninslice+j);
    
    /* time evolution of constraints */
    /* dampingmode: 0 -> viscous damping (default); 
       1 -> Nose-Hoover; */
    if(0.0 == ttt){
      *(gamma_ftg+j)  = 0.0;
    } 
    else if (dampingmode == 1) {
      *(gamma_ftg+j) += timestep * ( *(E_kin_2+j) / ttt - 1.0) * gamma_bar;
    } 
    else if (dampingmode == 0) {
      *(gamma_ftg+j)  =            (1.0 - ttt /  *(E_kin_2+j)) * gamma_bar;    
    }
  }


}

#else

void move_atoms_ftg(void) 
{
  if (myid==0)
  error("the chosen ensemble FTG is not supported by this binary");
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
        tmp = SQR((ORT(p,i,X)-center.x)/stadium.x) +
              SQR((ORT(p,i,Y)-center.y)/stadium.y) - 1;
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
        kin_energie_1[ensindex] += 
          SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / MASSE(p,i);

        /* new momenta */
	sort = VSORTE(p,i);
	IMPULS(p,i,X) = (IMPULS(p,i,X)*reibung + timestep * KRAFT(p,i,X))
                          * eins_d_reib * (restrictions + sort)->x;
        IMPULS(p,i,Y) = (IMPULS(p,i,Y)*reibung + timestep * KRAFT(p,i,Y))
                          * eins_d_reib * (restrictions + sort)->y;

        /* twice the new kinetic energy */ 
        kin_energie_2[ensindex] +=
          SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / MASSE(p,i);

        /* new positions */
        tmp = timestep * MASSE(p,i);
        ORT(p,i,X) += tmp * IMPULS(p,i,X);
        ORT(p,i,Y) += tmp * IMPULS(p,i,Y);
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

      Ekin_1 = SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / MASSE(p,i);
      px = IMPULS(p,i,X);

      /* new momenta */
      IMPULS(p,i,X) += timestep * KRAFT(p,i,X); 
      IMPULS(p,i,Y) += timestep * KRAFT(p,i,Y); 
#ifndef TWOD
      IMPULS(p,i,Z) += timestep * KRAFT(p,i,Z); 
#endif

      /* twice the new kinetic energy */ 
      Ekin_2 = SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / MASSE(p,i);
      px = (px + IMPULS(p,i,X)) / 2.0;

      /* new positions */
      ORT(p,i,X) += timestep * IMPULS(p,i,X) / MASSE(p,i);
      ORT(p,i,Y) += timestep * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
      ORT(p,i,Z) += timestep * IMPULS(p,i,Z) / MASSE(p,i);
#endif

      tot_kin_energy += (Ekin_1 + Ekin_2) / 4.0;

      /* which layer */
      num = scale * ORT(p,i,X);
      if (num < 0)             num = 0;
      if (num >= tran_nlayers) num = tran_nlayers-1;

      /* temperature control and heat conductivity */
      if (num == 0) {
        Ekin_left += Ekin_2;
        inv_mass_left     += 1/(MASSE(p,i));
        tot_impuls_left.x += IMPULS(p,i,X);
        tot_impuls_left.y += IMPULS(p,i,Y);
#ifndef TWOD
        tot_impuls_left.z += IMPULS(p,i,Z);
#endif
        natoms_left++;
      } else if  (num == nhalf) {
        Ekin_right += Ekin_2;
        inv_mass_right     += 1/(MASSE(p,i));
        tot_impuls_right.x += IMPULS(p,i,X);
        tot_impuls_right.y += IMPULS(p,i,Y);
#ifndef TWOD
        tot_impuls_right.z += IMPULS(p,i,Z);
#endif
        natoms_right++;
      } else if (num < nhalf) {
        heat_cond += (HEATCOND(p,i) + (Ekin_1 + Ekin_2) / 2) * px / MASSE(p,i);
      } else {
        heat_cond -= (HEATCOND(p,i) + (Ekin_1 + Ekin_2) / 2) * px / MASSE(p,i);
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
      num = scale * ORT(p,i,X);
      if (num < 0)             num = 0;
      if (num >= tran_nlayers) num = tran_nlayers-1;

      /* rescale momenta */
      if (num == 0) {
        IMPULS(p,i,X) = (IMPULS(p,i,X) - tot_impuls_left.x) * rescale;
        IMPULS(p,i,Y) = (IMPULS(p,i,Y) - tot_impuls_left.y) * rescale;
#ifndef TWOD
        IMPULS(p,i,Z) = (IMPULS(p,i,Z) - tot_impuls_left.z) * rescale;
#endif
      } else if (num == nhalf) {
        IMPULS(p,i,X) = (IMPULS(p,i,X) - tot_impuls_right.x) * Rescale;
        IMPULS(p,i,Y) = (IMPULS(p,i,Y) - tot_impuls_right.y) * Rescale;
#ifndef TWOD
        IMPULS(p,i,Z) = (IMPULS(p,i,Z) - tot_impuls_right.z) * Rescale;
#endif
      }

#ifdef STRESS_TENS
        PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
        PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
        PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
        PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
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


/*****************************************************************************
*
* Move the atoms for the  Conjugated Gradient
*
*****************************************************************************/

#if defined(CG) 

void move_atoms_cg(real alpha)
{
  int k;
#ifdef CLONE
  int clones;
#endif
  /* loop over all cells */
  for (k=0; k<ncells; ++k) {

    int  i, sort;
    cell *p;

    p = cell_array + CELLS(k);

    for (i=0; i<p->n; ++i) {
        /* CG:  move atoms in search direction for linmin */
        ORT(p,i,X) = OLD_ORT(p,i,X) + alpha * CG_H(p,i,X);
        ORT(p,i,Y) = OLD_ORT(p,i,Y) + alpha * CG_H(p,i,Y);
#ifndef TWOD
        ORT(p,i,Z) = OLD_ORT(p,i,Z) + alpha * CG_H(p,i,Z);
#endif
#ifdef CLONE /* i hope that the clone atoms will experience the same force 
		as the org. atoms, so i'm not changing anything else... */
	for(clones=1;clones<=nclones;clones++)
	    {
		ORT(p,i,X) = OLD_ORT(p,i,X) + alpha * CG_H(p,i+clones,X);
		ORT(p,i,Y) = OLD_ORT(p,i,Y) + alpha * CG_H(p,i+clones,Y);
#ifndef TWOD
		ORT(p,i,Z) = OLD_ORT(p,i,Z) + alpha * CG_H(p,i+clones,Z);
#endif
	    }
    i+=nclones;
#endif
    }
  }
}
#else

void move_atoms_cg(real alpha) 
{
  if (myid==0)
  error("the chosen ensemble CG is not supported by this binary");
}

#endif

