
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2011 Institute for Theoretical and Applied Physics,
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

#if defined(NVE) || defined(EPITAX)

void move_atoms_nve(void)
{
  int k;
  static int count = 0;
  real tmp_f_max2=0.0;
  real tmp_x_max2=0.0;
#ifdef DAMP
  real tmp1, tmp2, tmp3, f, maxax, maxax2;
#endif
#ifdef BER
  real cc;
#endif
#ifdef BER /* reproduce Ju Li's implementation of the Berendsen Thermostat */
     cc = 1. - timestep/ tauber * ( ( 2.0 * tot_kin_energy / nactive+8.6174101569719990e-06) / (temperature+8.6174101569719990e-06) - 1. );
     /* this would be the standard formula */ 
   //cc = 1.+timestep/tauber*((temperature+8.6174101569719990e-06)/(2.0*tot_kin_energy/nactive+8.6174101569719990e-06)- 1. );

   
      if (cc < 0.5) cc = 0.5;
      else if (cc > 2.0) cc = 2.0;
      cc=sqrt(cc);
#endif
      
  /* epitax may call this routine for other ensembles,
     in which case we do not reset tot_kin_energy */
  if ((ensemble==ENS_NVE) || (ensemble==ENS_GLOK)) tot_kin_energy = 0.0;

  fnorm   = 0.0;
  xnorm   = 0.0;
  pnorm   = 0.0;
  PxF     = 0.0;
  omega_E = 0.0;
#ifdef DAMP
  n_damp  = 0;
  tot_kin_energy_damp = 0.0;
#endif

#ifdef SHOCK
  if (do_press_calc) calc_pxavg();
#endif




  
  /* loop over all cells */
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tot_kin_energy,fnorm,omega_E,PxF,pnorm)
#endif
  for (k=0; k<NCELLS; ++k) { /* loop over all cells */

    int  i,j, sort;
    cell *p;
    real kin_energie_1, kin_energie_2, tmp;
#ifdef UNIAX    
    real rot_energie_1, rot_energie_2;
    real dot, norm;
    vektor cross;
#endif
#ifdef RIGID
    int satom;
    real relmass;
#endif

#ifdef DAMP
    real kin_energie_damp_1,kin_energie_damp_2,tmp2,rampedtemp,zeta_finnis;
#endif

    p = CELLPTR(k);

#ifdef CLONE
    for (i=0; i<p->n; i+=nclones)
      for (j=1; j<nclones; j++) {
        KRAFT(p,i+j,X)  = KRAFT(p,i,X);
        KRAFT(p,i+j,Y)  = KRAFT(p,i,Y);
#ifndef TWOD
        KRAFT(p,i+j,Z)  = KRAFT(p,i,Z);
#endif
        IMPULS(p,i+j,X) = IMPULS(p,i,X);
        IMPULS(p,i+j,Y) = IMPULS(p,i,Y);
#ifndef TWOD
        IMPULS(p,i+j,Z) = IMPULS(p,i,Z);
#endif
      }
#endif /* CLONE */

#ifdef SX
#pragma vdir vector,nodep
#endif

    for (i=0; i<p->n; ++i) { /* loop over all atoms in the cell */

#ifdef EPITAX 
      /* beam atoms are always integrated by NVE */
      if ( (ensemble != ENS_NVE) &&
           (NUMMER(p,i) <= epitax_sub_n) && 
           (POTENG(p,i) <= epitax_ctrl * epitax_poteng_min) ) continue;
#endif
#ifndef DAMP
      kin_energie_1 = SPRODN(IMPULS,p,i,IMPULS,p,i);
#endif

#ifdef UNIAX
      rot_energie_1 = SPRODN(DREH_IMPULS,p,i,DREH_IMPULS,p,i);
#endif

      sort = VSORTE(p,i);

#ifdef RIGID
      if ( superatom[sort] > -1 ) {

        satom   = superatom[sort];
        relmass = MASSE(p,i) / supermass[satom];

        if ( (superrestrictions + satom)->x )
          KRAFT(p,i,X) = superforce[satom].x * relmass; 
        if ( (superrestrictions + satom)->y )
          KRAFT(p,i,Y) = superforce[satom].y * relmass;
#ifndef TWOD
        if ( (superrestrictions + satom)->z )
          KRAFT(p,i,Z) = superforce[satom].z * relmass;
#endif
      }
#endif

#if defined(FBC) && !defined(RIGID) 
      /* give virtual particles their extra force */
      KRAFT(p,i,X) += (fbc_forces + sort)->x;
      KRAFT(p,i,Y) += (fbc_forces + sort)->y;
#ifndef TWOD
      KRAFT(p,i,Z) += (fbc_forces + sort)->z;
#endif
#endif

#if defined(FBC) && defined(BEND)
      /* give virtual particles their extra force */
      KRAFT(p,i,X) += (bend_forces + sort)->x;
      KRAFT(p,i,Y) += (bend_forces + sort)->y;
#ifndef TWOD
      KRAFT(p,i,Z) += (bend_forces + sort)->z;
#endif
#endif


      

      /* and set their force (->momentum) in restricted directions to 0 */
      KRAFT(p,i,X) *= (restrictions + sort)->x;
      KRAFT(p,i,Y) *= (restrictions + sort)->y;
#ifndef TWOD
      KRAFT(p,i,Z) *= (restrictions + sort)->z;
#endif

#ifdef FNORM
      fnorm += SPRODN(KRAFT,p,i,KRAFT,p,i);
      /* determine the biggest force component */
      tmp_f_max2 = MAX(SQR(KRAFT(p,i,X)),tmp_f_max2);
      tmp_f_max2 = MAX(SQR(KRAFT(p,i,Y)),tmp_f_max2);
#ifndef TWOD
      tmp_f_max2 = MAX(SQR(KRAFT(p,i,Z)),tmp_f_max2);
#endif
#endif
#ifdef EINSTEIN
      omega_E += SPRODN(KRAFT,p,i,KRAFT,p,i) / MASSE(p,i);
#endif

#ifndef DAMP /*  Normal NVE */
      IMPULS(p,i,X) += timestep * KRAFT(p,i,X);
      IMPULS(p,i,Y) += timestep * KRAFT(p,i,Y);
#ifndef TWOD
      IMPULS(p,i,Z) += timestep * KRAFT(p,i,Z);
#endif
#else  /* Damping layers */

      /* use a local thermostat: Finnis
         We fix a temperature gradient from temp to zero
         in the same way the damping constant is ramped up.
         The mean temperature in the damping layers
         has no meaning, only temperature of the inner part is output */
      /* the stadium function for each atom could also be calculated 
         in forces_nbl to save time */
      /* it is the users responsability that stadium.i / stadium2.i
         is equal for all i */

      maxax = MAX(MAX(stadium.x,stadium.y),stadium.z);
      maxax2 = MAX(MAX(stadium2.x,stadium2.y),stadium2.z);

      /* Calculate stadium function f */
      tmp1 = (stadium2.x==0) ? 0 : SQR((ORT(p,i,X)-center.x)/(2.0*stadium2.x));
      tmp2 = (stadium2.y==0) ? 0 : SQR((ORT(p,i,Y)-center.y)/(2.0*stadium2.y));
      tmp3 = (stadium2.z==0) ? 0 : SQR((ORT(p,i,Z)-center.z)/(2.0*stadium2.z));

      f    = (tmp1+tmp2+tmp3-SQR(maxax / (2.0*maxax2)) ) /
                (0.25 - SQR(maxax / (2.0*maxax2)) );
      if      (f<= 0.0) f = 0.0;
      else if (f>1.0)   f = 1.0;

      /* we smooth the stadium function: to get a real bath tub !*/
      DAMPF(p,i) = .5 * (1 + sin(-M_PI/2.0 + M_PI*f));


      if (DAMPF(p,i) == 0.0) { /* take care of possible rounding errors ? */

        kin_energie_1 = SPRODN(IMPULS,p,i,IMPULS,p,i);

        IMPULS(p,i,X) += timestep * KRAFT(p,i,X);
        IMPULS(p,i,Y) += timestep * KRAFT(p,i,Y);
#ifndef TWOD
        IMPULS(p,i,Z) += timestep * KRAFT(p,i,Z);
#endif
        kin_energie_2 = SPRODN(IMPULS,p,i,IMPULS,p,i);
        tot_kin_energy += (kin_energie_1 + kin_energie_2) / (4 * MASSE(p,i));
      }
      else {
        kin_energie_damp_1 = SPRODN(IMPULS,p,i,IMPULS,p,i);

        tmp   = kin_energie_damp_1 / MASSE(p,i); /* local temp */
        tmp2  = (restrictions + sort)->x + (restrictions + sort)->y;
#ifndef TWOD
        tmp2 += (restrictions + sort)->z;
#endif
        n_damp += tmp2;

        if (tmp2 != 0) tmp /= tmp2;

        /* to account for restricted mobilities */
	rampedtemp  = (tmp2 !=0) ? 

          (tmp2/3.0 * damptemp * (1.0 - DAMPF(p,i))) : 0.0;

        if(rampedtemp !=0.0) {
          zeta_finnis = zeta_0 * (tmp-rampedtemp)
             / sqrt(SQR(tmp) + SQR(rampedtemp*delta_finnis)) * DAMPF(p,i);
        }
	else zeta_finnis = zeta_0;

        /* new momenta */
        IMPULS(p,i,X) += (-1.0*IMPULS(p,i,X) * zeta_finnis + KRAFT(p,i,X)) 
                         * timestep * (restrictions + sort)->x ;
        IMPULS(p,i,Y) += (-1.0*IMPULS(p,i,Y) * zeta_finnis + KRAFT(p,i,Y)) 
                         * timestep * (restrictions + sort)->y;
#ifndef TWOD
        IMPULS(p,i,Z) += (-1.0*IMPULS(p,i,Z) * zeta_finnis + KRAFT(p,i,Z)) 
                         * timestep * (restrictions + sort)->z;
#endif
        kin_energie_damp_2 =  SPRODN(IMPULS,p,i,IMPULS,p,i);
        tot_kin_energy_damp += (kin_energie_damp_1 + kin_energie_damp_2) 
                               / (4.0 * MASSE(p,i)) ;
      }
#endif /* DAMP */

#if defined (GLOK)|| defined(MIX)

      /* "Global Convergence": */ 
      /* like mik, just with the global force and momentum vectors */
      /* change to velocity norm, change names later... */
      PxF   += SPRODN(IMPULS,p,i,KRAFT,p,i)/MASSE(p,i);
      pnorm += SPRODN(IMPULS,p,i,IMPULS,p,i)/MASSE(p,i)/MASSE(p,i);

#ifdef MIX
      /* global version of MIX  with adaptive mixing */  
      /* 'turn' the velocities a little bit more along the forces... */
      IMPULS(p,i,X) = (1.0-mix)*IMPULS(p,i,X) + mix * KRAFT(p,i,X)
                      * mixforcescalefac * MASSE(p,i);
      IMPULS(p,i,Y) = (1.0-mix)*IMPULS(p,i,Y) + mix * KRAFT(p,i,Y)
                      * mixforcescalefac * MASSE(p,i);
      IMPULS(p,i,Z) = (1.0-mix)*IMPULS(p,i,Z) + mix * KRAFT(p,i,Z)
                      * mixforcescalefac * MASSE(p,i);
#endif   

#endif /* GLOK || MIX */

#ifdef UNIAX
      dot = 2.0 * SPRODN(DREH_IMPULS,p,i,ACHSE,p,i);
      DREH_IMPULS(p,i,X) += timestep * DREH_MOMENT(p,i,X)
                             - dot * ACHSE(p,i,X);
      DREH_IMPULS(p,i,Y) += timestep * DREH_MOMENT(p,i,Y)
                             - dot * ACHSE(p,i,Y);
      DREH_IMPULS(p,i,Z) += timestep * DREH_MOMENT(p,i,Z)
                             - dot * ACHSE(p,i,Z);
#endif

#ifndef DAMP
      kin_energie_2 = SPRODN(IMPULS,p,i,IMPULS,p,i);
#endif
#ifdef UNIAX
      rot_energie_2 = SPRODN(DREH_IMPULS,p,i,DREH_IMPULS,p,i); 
#endif
#ifndef DAMP
      tot_kin_energy += (kin_energie_1 + kin_energie_2) / (4 * MASSE(p,i));
#endif
#ifdef UNIAX
      tot_kin_energy += (rot_energie_1 + rot_energie_2) / (4 * uniax_inert);
#endif	  

#ifdef BER
      //  if(steps%16==0)
      {
          IMPULS(p,i,X) *= cc;
          IMPULS(p,i,Y) *= cc;
#ifndef TWOD
          IMPULS(p,i,Z) *= cc;
#endif
      }
#endif
      
      /* new positions */
      tmp = timestep / MASSE(p,i);
      ORT(p,i,X) += tmp * IMPULS(p,i,X);
      ORT(p,i,Y) += tmp * IMPULS(p,i,Y);
#ifndef TWOD
      ORT(p,i,Z) += tmp * IMPULS(p,i,Z);
#endif

#ifdef RELAXINFO
      xnorm += tmp * tmp * SPRODN(IMPULS,p,i,IMPULS,p,i);
      /* determine the biggest force component */
      tmp_x_max2 =  MAX(SQR(tmp*IMPULS(p,i,X)),tmp_x_max2);
      tmp_x_max2 =  MAX(SQR(tmp*IMPULS(p,i,Y)),tmp_x_max2);
#ifndef TWOD
      tmp_x_max2 =  MAX(SQR(tmp*IMPULS(p,i,Z)),tmp_x_max2);
#endif
#endif

#ifdef SHOCK
      if (shock_mode == 3 && steps > 0 ) {
        if (ORT(p,i,X) > box_x.x) {
          IMPULS(p,i,X) = -IMPULS(p,i,X);
          ORT(p,i,X) = 2 * box_x.x - ORT(p,i,X);
        }
      }
      if (shock_mode == 4) {
        real rand = shock_speed_l * timestep * steps ;
        if (ORT(p,i,X) < rand ) {
          IMPULS(p,i,X) = -IMPULS(p,i,X) + 2 * shock_speed_l * MASSE(p,i);
          ORT(p,i,X) = 2 * rand - ORT(p,i,X);
        }
        if (ORT(p,i,X) > box_x.x - rand ) {
          IMPULS(p,i,X) = -IMPULS(p,i,X) - 2 * shock_speed_r * MASSE(p,i);
          ORT(p,i,X) = 2 * ( box_x.x - rand ) - ORT(p,i,X);
        }
      }
#endif

#ifdef UNIAX
      /* new molecular axes */
      cross.x = DREH_IMPULS(p,i,Y) * ACHSE(p,i,Z)
              - DREH_IMPULS(p,i,Z) * ACHSE(p,i,Y);
      cross.y = DREH_IMPULS(p,i,Z) * ACHSE(p,i,X)
              - DREH_IMPULS(p,i,X) * ACHSE(p,i,Z);
      cross.z = DREH_IMPULS(p,i,X) * ACHSE(p,i,Y)
              - DREH_IMPULS(p,i,Y) * ACHSE(p,i,X);

      ACHSE(p,i,X) += timestep * cross.x / uniax_inert;
      ACHSE(p,i,Y) += timestep * cross.y / uniax_inert;
      ACHSE(p,i,Z) += timestep * cross.z / uniax_inert;

      norm = SQRT( SPRODN(ACHSE,p,i,ACHSE,p,i) );
	    
      ACHSE(p,i,X) /= norm;
      ACHSE(p,i,Y) /= norm;
      ACHSE(p,i,Z) /= norm;
#endif    

#ifdef STRESS_TENS
      if (do_press_calc) {
#ifdef SHOCK
        PRESSTENS(p,i,xx) += (IMPULS(p,i,X) - PXAVG(p,i)) * 
                             (IMPULS(p,i,X) - PXAVG(p,i))  / MASSE(p,i); 
        PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
        PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,zx) += (IMPULS(p,i,X) - PXAVG(p,i)) * 
                              IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,xy) += (IMPULS(p,i,X) - PXAVG(p,i)) *
                              IMPULS(p,i,Y) / MASSE(p,i);
#else /* not SHOCK */
        PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
        PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
        PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
        PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
#endif /* SHOCK */
      }
#endif /* STRESS_TENS */
    }
  }

#ifdef MPI
  { /* add up results from different CPUs */
    int nc = 0;
    double tmpvec1[8], tmpvec2[8];

    tmpvec1[nc++] = tot_kin_energy;
    tmpvec1[nc++] = fnorm;
    tmpvec1[nc++] = PxF;
    tmpvec1[nc++] = omega_E;
    tmpvec1[nc++] = pnorm;
    tmpvec1[nc++] = xnorm;
#ifdef DAMP
    tmpvec1[nc++] = tot_kin_energy_damp;
    tmpvec1[nc++] = n_damp;
#endif

    MPI_Allreduce( tmpvec1, tmpvec2, nc, REAL, MPI_SUM, cpugrid);

    nc = 0;
    tot_kin_energy      = tmpvec2[nc++];
    fnorm               = tmpvec2[nc++];
    PxF                 = tmpvec2[nc++];
    omega_E             = tmpvec2[nc++];
    pnorm               = tmpvec2[nc++];
    xnorm               = tmpvec2[nc++];
#ifdef DAMP
    tot_kin_energy_damp = tmpvec2[nc++];
    n_damp              = tmpvec2[nc++];
#endif
  }
#ifdef FNORM
  MPI_Allreduce( &tmp_f_max2, &f_max2, 1, REAL, MPI_MAX, cpugrid);
#endif
#ifdef RELAXINFO
  MPI_Allreduce( &tmp_x_max2, &x_max2, 1, REAL, MPI_MAX, cpugrid);
#endif

#else /* not MPI */

#ifdef FNORM
  f_max2 = tmp_f_max2;
#endif
#ifdef RELAXINFO
  x_max2 = tmp_x_max2;
#endif

#endif /* MPI */

#if defined (GLOK) || defined (MIX)
  PxF /= (SQRT(fnorm) * SQRT(pnorm));
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
* Basic NVE Integrator with TTM (two temperature model), stripped of most other options
*
*****************************************************************************/

#if defined(TTM)

void move_atoms_ttm(void)
{
  int k;
/*  static int count = 0;*/
  real tmpvec1[8], tmpvec2[8];

  tot_kin_energy = 0.0;
#ifdef DEBUG
  E_ph_auf_local = 0.0;
#endif /*DEBUG*/
  omega_E = 0.0;

  /* loop over all cells */
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tot_kin_energy,omega_E)
#endif
  for (k=0; k<NCELLS; ++k) { /* loop over all cells */

    int  i,j, sort;
    int fd_i, fd_j, fd_k;
    double fd_xi;
    cell *p;
    real kin_energie_1, kin_energie_2, tmp;

    p = CELLPTR(k);

    
    fd_i=p->fd_cell_idx.x;
    fd_j=p->fd_cell_idx.y;
    fd_k=p->fd_cell_idx.z;

    /* get coupling constant, if fd cell is inactive, no coupling */
    fd_xi=(l1[fd_i][fd_j][fd_k].natoms>=fd_min_atoms)?(l1[fd_i][fd_j][fd_k].xi):(0.0);
    
    for (i=0; i<p->n; ++i) { /* loop over all atoms in the cell */
#ifdef DEBUG
      double delta_E_atom;
#endif

      kin_energie_1 = SPRODN(IMPULS,p,i,IMPULS,p,i);
      sort = VSORTE(p,i);


#if defined(FBC)    /* give virtual particles their extra force */
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

#ifdef EINSTEIN
      omega_E += SPRODN(KRAFT,p,i,KRAFT,p,i) / MASSE(p,i);
#endif

#ifdef DEBUG
      /* dE = dt * v * F(el->ph)
       * mit F(el->ph) = xi*m*v_therm
       * ******************************* */
      delta_E_atom = 1/MASSE(p,i) * timestep
	             * (  IMPULS(p,i,X) * fd_xi * MASSE(p,i) * ( IMPULS(p,i,X)/MASSE(p,i) - l1[fd_i][fd_j][fd_k].v_com.x)
			+ IMPULS(p,i,Y) * fd_xi * MASSE(p,i) * ( IMPULS(p,i,Y)/MASSE(p,i) - l1[fd_i][fd_j][fd_k].v_com.y)
#ifndef TWOD
			+ IMPULS(p,i,Z) * fd_xi * MASSE(p,i) * ( IMPULS(p,i,Z)/MASSE(p,i) - l1[fd_i][fd_j][fd_k].v_com.z)
#endif /*TWOD*/
		       );
      E_ph_auf_local += delta_E_atom;
#endif /*DEBUG*/


      /* TTM: p += (F + xi*m*v_therm) * dt  
       * ********************************** */
      
      IMPULS(p,i,X) += timestep * ( KRAFT(p,i,X) + fd_xi * MASSE(p,i) * ( IMPULS(p,i,X)/MASSE(p,i) - l1[fd_i][fd_j][fd_k].v_com.x) );
      IMPULS(p,i,Y) += timestep * ( KRAFT(p,i,Y) + fd_xi * MASSE(p,i) * ( IMPULS(p,i,Y)/MASSE(p,i) - l1[fd_i][fd_j][fd_k].v_com.y) );
#ifndef TWOD
      IMPULS(p,i,Z) += timestep * ( KRAFT(p,i,Z) + fd_xi * MASSE(p,i) * ( IMPULS(p,i,Z)/MASSE(p,i) - l1[fd_i][fd_j][fd_k].v_com.z) );
#endif


      kin_energie_2 = SPRODN(IMPULS,p,i,IMPULS,p,i);

      tot_kin_energy += (kin_energie_1 + kin_energie_2) / (4 * MASSE(p,i));


      /* new positions */
      tmp = timestep / MASSE(p,i);
      ORT(p,i,X) += tmp * IMPULS(p,i,X);
      ORT(p,i,Y) += tmp * IMPULS(p,i,Y);
#ifndef TWOD
      ORT(p,i,Z) += tmp * IMPULS(p,i,Z);
#endif

#ifdef STRESS_TENS
      if (do_press_calc) {
        PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
        PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
        PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
        PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
      }
#endif /* STRESS_TENS */
    }
  }

#ifdef MPI
  /* add up results from different CPUs */
  tmpvec1[0] = tot_kin_energy;
  tmpvec1[3] = omega_E;
#ifdef DEBUG
  tmpvec1[7] = E_ph_auf_local;
#endif /*DEBUG*/

  /*  MPI_Allreduce( tmpvec1, tmpvec2, 5, REAL, MPI_SUM, cpugrid); */
  MPI_Allreduce( tmpvec1, tmpvec2, 8, REAL, MPI_SUM, cpugrid); 


  tot_kin_energy = tmpvec2[0];
  omega_E        = tmpvec2[3];
#ifdef DEBUG
  E_ph_auf      += tmpvec2[7];
#endif /*DEBUG*/

#endif /* MPI */


}
#else

void move_atoms_ttm(void) 
{
  if (myid==0)
  error("the chosen ensemble TTM is not supported by this binary");
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
  real tmpvec1[3], tmpvec2[3];
  real tmp_f_max2=0.0;
  real tmp_x_max2=0.0;
  real mass = 0.006084;


  static int count = 0;
 

  /* implementation of adaptive mik time step */


 tot_kin_energy = 0.0;
  fnorm   = 0.0;
  xnorm   = 0.0;
  //  pnorm   = 0.0;

#ifdef AND
  /* Andersen Thermostat -- Initialize the velocities now and then */
  ++count;
  if ((tempintv!=0) && (0==count%tempintv)) maxwell(temperature);
#endif

  /* loop over all cells */
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tot_kin_energy,fnorm)
#endif
  for (k=0; k<NCELLS; ++k) {

    int  i, j, sort;
    cell *p;
    real kin_energie_1, kin_energie_2, tmp;
#ifdef RIGID
    int satom;
    real relmass;
#endif

    p = CELLPTR(k);

#ifdef CLONE
    for (i=0; i<p->n; i+=nclones)
      for (j=1; j<nclones; j++) {
        KRAFT(p,i+j,X)  = KRAFT(p,i,X);
        KRAFT(p,i+j,Y)  = KRAFT(p,i,Y);
#ifndef TWOD
        KRAFT(p,i+j,Z)  = KRAFT(p,i,Z);
#endif
        IMPULS(p,i+j,X) = IMPULS(p,i,X);
        IMPULS(p,i+j,Y) = IMPULS(p,i,Y);
#ifndef TWOD
        IMPULS(p,i+j,Z) = IMPULS(p,i,Z);
#endif
      }
#endif /* CLONE */

#ifdef SX
#pragma vdir vector,nodep
#endif
    for (i=0; i<p->n; ++i) {

#ifdef EPITAX
        /* only substrate atoms are integrated by MIK */
        if ( (NUMMER(p,i) > epitax_sub_n) && 
             (POTENG(p,i) > epitax_ctrl * epitax_poteng_min) ) continue;
#endif

        kin_energie_1 = SPRODN(IMPULS,p,i,IMPULS,p,i);

	sort = VSORTE(p,i);

#ifdef RIGID
	if ( superatom[sort] > -1 ) {

	  satom   = superatom[sort];
	  relmass = MASSE(p,i) / supermass[satom];

	  if ( (superrestrictions + satom)->x )
	    KRAFT(p,i,X) = superforce[satom].x * relmass; 
	  if ( (superrestrictions + satom)->y )
	    KRAFT(p,i,Y) = superforce[satom].y * relmass;
#ifndef TWOD
	  if ( (superrestrictions + satom)->z )
	    KRAFT(p,i,Z) = superforce[satom].z * relmass;
#endif
	}
#endif

#if defined(FBC) && !defined(RIGID)
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
	fnorm += SPRODN(KRAFT,p,i,KRAFT,p,i);
        /* determine the biggest force component */
        tmp_f_max2 = MAX(SQR(KRAFT(p,i,X)),tmp_f_max2);
        tmp_f_max2 = MAX(SQR(KRAFT(p,i,Y)),tmp_f_max2);
#ifndef TWOD
        tmp_f_max2 = MAX(SQR(KRAFT(p,i,Z)),tmp_f_max2);
#endif
#endif
        IMPULS(p,i,X) += timestep * KRAFT(p,i,X);
        IMPULS(p,i,Y) += timestep * KRAFT(p,i,Y);
#ifndef TWOD
        IMPULS(p,i,Z) += timestep * KRAFT(p,i,Z);
#endif

	/* Mikroconvergence Algorithm - set velocity zero if a*v < 0 */
	if (0.0 > SPRODN(IMPULS,p,i,KRAFT,p,i) ) {
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
#ifdef RELAXINFO
	//	pnorm += SPRODN(IMPULS,p,i,IMPULS,p,i)/MASSE(p,i)/MASSE(p,i);
	xnorm   += tmp * tmp* SPRODN(IMPULS,p,i,IMPULS,p,i);
        /* determine the biggest force component */
        tmp_x_max2 =  MAX(SQR(tmp*IMPULS(p,i,X)),tmp_x_max2);
        tmp_x_max2 =  MAX(SQR(tmp*IMPULS(p,i,Y)),tmp_x_max2);
#ifndef TWOD
        tmp_x_max2 =  MAX(SQR(tmp*IMPULS(p,i,Z)),tmp_x_max2);
#endif
#endif
        kin_energie_2 = SPRODN(IMPULS,p,i,IMPULS,p,i);

        /* sum up kinetic energy on this CPU */ 
        tot_kin_energy += (kin_energie_1 + kin_energie_2) / (4.0 * MASSE(p,i));

#ifdef STRESS_TENS
        if (do_press_calc) {
          PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
          PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
          PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
          PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
          PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
          PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
        }
#endif
    }
  }

#ifdef MPI
  /* add up results from different CPUs */
  tmpvec1[0] = tot_kin_energy;
  tmpvec1[1] = fnorm;
  tmpvec1[2] = xnorm;

  MPI_Allreduce( tmpvec1, tmpvec2, 3, REAL, MPI_SUM, cpugrid);

  tot_kin_energy = tmpvec2[0];
  fnorm          = tmpvec2[1];
  xnorm          = tmpvec2[2];

#ifdef FNORM
  MPI_Allreduce( &tmp_f_max2, &f_max2, 1, REAL, MPI_MAX, cpugrid);
#endif
#ifdef RELAXINFO
  MPI_Allreduce( &tmp_x_max2, &x_max2, 1, REAL, MPI_MAX, cpugrid);
#endif
#else
#ifdef FNORM
  f_max2 = tmp_f_max2;
#endif
#ifdef RELAXINFO
  x_max2 = tmp_x_max2;
#endif

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
  real tmpvec1[6], tmpvec2[6], ttt;

  real E_kin_1 = 0.0, E_kin_2 = 0.0;
  real reibung, eins_d_reib;
  real E_rot_1 = 0.0, E_rot_2 = 0.0;
  real tmp_f_max2=0.0;
  real tmp_x_max2=0.0;
#ifdef UNIAX
  real reibung_rot,  eins_d_reib_rot;
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
  for (k=0; k<NCELLS; ++k) {  /* loop over cells */

    int i, j, sort;
    cell *p;
    real tmp;
#ifdef UNIAX
    real dot, norm ;
    vektor cross ;
#endif
#ifdef RIGID
    int satom;
    real relmass;
#endif

    p = CELLPTR(k);

#ifdef CLONE
    for (i=0; i<p->n; i+=nclones)
      for (j=1; j<nclones; j++) {
        KRAFT(p,i+j,X)  = KRAFT(p,i,X);
        KRAFT(p,i+j,Y)  = KRAFT(p,i,Y);
#ifndef TWOD
        KRAFT(p,i+j,Z)  = KRAFT(p,i,Z);
#endif
        IMPULS(p,i+j,X) = IMPULS(p,i,X);
        IMPULS(p,i+j,Y) = IMPULS(p,i,Y);
#ifndef TWOD
        IMPULS(p,i+j,Z) = IMPULS(p,i,Z);
#endif
      }
#endif /* CLONE */

    for (i=0; i<p->n; ++i) {  /* loop over atoms */

#ifdef EPITAX
        /* only substrate atoms are integrated by NVT */
        if ( (NUMMER(p,i) > epitax_sub_n) && 
             (POTENG(p,i) > epitax_ctrl * epitax_poteng_min) ) continue;
#endif

        /* twice the old kinetic energy */
        E_kin_1 += SPRODN(IMPULS,p,i,IMPULS,p,i) / MASSE(p,i);
#ifdef UNIAX
        E_rot_1 += SPRODN(DREH_IMPULS,p,i,DREH_IMPULS,p,i) / uniax_inert;
#endif

	sort = VSORTE(p,i);

#ifdef RIGID
	if ( superatom[sort] > -1 ) {

	  satom   = superatom[sort];
	  relmass = MASSE(p,i) / supermass[satom];

	  if ( (superrestrictions + satom)->x )
	    KRAFT(p,i,X) = superforce[satom].x * relmass; 
	  if ( (superrestrictions + satom)->y )
	    KRAFT(p,i,Y) = superforce[satom].y * relmass;
#ifndef TWOD
	  if ( (superrestrictions + satom)->z )
	    KRAFT(p,i,Z) = superforce[satom].z * relmass;
#endif
	}
#endif
#if defined(FBC) && !defined(RIGID)
        /* give virtual particles their extra force */
	KRAFT(p,i,X) += (fbc_forces + sort)->x;
	KRAFT(p,i,Y) += (fbc_forces + sort)->y;
#ifndef TWOD
	KRAFT(p,i,Z) += (fbc_forces + sort)->z;
#endif
#endif
        
#if defined(FBC) && defined(BEND)
      /* give virtual particles their extra force */
      KRAFT(p,i,X) += (bend_forces + sort)->x;
      KRAFT(p,i,Y) += (bend_forces + sort)->y;
#ifndef TWOD
      KRAFT(p,i,Z) += (bend_forces + sort)->z;
#endif
#endif


      KRAFT(p,i,X) *= (restrictions + sort)->x;
      KRAFT(p,i,Y) *= (restrictions + sort)->y;
#ifndef TWOD
      KRAFT(p,i,Z) *= (restrictions + sort)->z;
#endif
#ifdef FNORM

	fnorm   += SPRODN(KRAFT,p,i,KRAFT,p,i);
        /* determine the biggest force component */
        tmp_f_max2 = MAX(SQR(KRAFT(p,i,X)),tmp_f_max2);
        tmp_f_max2 = MAX(SQR(KRAFT(p,i,Y)),tmp_f_max2);

#ifndef TWOD
        tmp_f_max2 = MAX(SQR(KRAFT(p,i,Z)),tmp_f_max2);
#endif
#endif
#ifdef EINSTEIN
	omega_E += SPRODN(KRAFT,p,i,KRAFT,p,i) / MASSE(p,i);
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
        dot = 2.0 * SPRODN(DREH_IMPULS,p,i,ACHSE,p,i);

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
        E_kin_2 += SPRODN(IMPULS,p,i,IMPULS,p,i) / MASSE(p,i);
#ifdef UNIAX
        E_rot_2 += SPRODN(DREH_IMPULS,p,i,DREH_IMPULS,p,i) / uniax_inert;
#endif

        /* new positions */
        tmp = timestep / MASSE(p,i);
        ORT(p,i,X) += tmp * IMPULS(p,i,X);
        ORT(p,i,Y) += tmp * IMPULS(p,i,Y);
#ifndef TWOD
        ORT(p,i,Z) += tmp * IMPULS(p,i,Z);
#endif
#ifdef RELAXINFO
	xnorm   += tmp * tmp* SPRODN(IMPULS,p,i,IMPULS,p,i);
        /* determine the biggest force component */
        tmp_x_max2 =  MAX(SQR(tmp*IMPULS(p,i,X)),tmp_x_max2);
        tmp_x_max2 =  MAX(SQR(tmp*IMPULS(p,i,Y)),tmp_x_max2);
#ifndef TWOD
        tmp_x_max2 =  MAX(SQR(tmp*IMPULS(p,i,Z)),tmp_x_max2);
#endif
#endif
#ifdef UNIAX
        cross.x = DREH_IMPULS(p,i,Y) * ACHSE(p,i,Z)
                - DREH_IMPULS(p,i,Z) * ACHSE(p,i,Y);
        cross.y = DREH_IMPULS(p,i,Z) * ACHSE(p,i,X)
                - DREH_IMPULS(p,i,X) * ACHSE(p,i,Z);
        cross.z = DREH_IMPULS(p,i,X) * ACHSE(p,i,Y)
                - DREH_IMPULS(p,i,Y) * ACHSE(p,i,X);

        ACHSE(p,i,X) += timestep * cross.x / uniax_inert;
        ACHSE(p,i,Y) += timestep * cross.y / uniax_inert;
        ACHSE(p,i,Z) += timestep * cross.z / uniax_inert;

        norm = SQRT( SPRODN(ACHSE,p,i,ACHSE,p,i) );

        ACHSE(p,i,X) /= norm;
        ACHSE(p,i,Y) /= norm;
        ACHSE(p,i,Z) /= norm;
#endif

#ifdef STRESS_TENS
        if (do_press_calc) {
          PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
          PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
          PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
          PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
          PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
          PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
	}
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
  tmpvec1[5] = xnorm;

  MPI_Allreduce( tmpvec1, tmpvec2, 6, REAL, MPI_SUM, cpugrid);

  tot_kin_energy = tmpvec2[0];
  E_kin_2        = tmpvec2[1];
  E_rot_2        = tmpvec2[2];
  fnorm          = tmpvec2[3];
  omega_E        = tmpvec2[4];
  xnorm          = tmpvec2[5];

#ifdef FNORM
  MPI_Allreduce( &tmp_f_max2, &f_max2, 1, REAL, MPI_MAX, cpugrid);
#endif
#ifdef RELAXINFO
  MPI_Allreduce( &tmp_x_max2, &x_max2, 1, REAL, MPI_MAX, cpugrid);
#endif
#else
#ifdef FNORM
  f_max2 = tmp_f_max2;
#endif
#ifdef RELAXINFO
  x_max2 = tmp_x_max2;
#endif
#endif /* MPI */

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
  real tmpvec1[4], tmpvec2[4], ttt;
  real E_kin_1 = 0.0, E_kin_2 = 0.0;
  vektor reibung, eins_d_reib;
  real E_rot_1 = 0.0, E_rot_2 = 0.0;
 real tmp_f_max2=0.0;
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
  for (k=0; k<NCELLS; ++k) {

    int i;
    int sort;
    cell *p;
    real tmp;
#ifdef UNIAX
    real dot, norm ;
    vektor cross ;
#endif
    p = CELLPTR(k);

    for (i=0; i<p->n; ++i) {

        /* twice the old kinetic energy */
        E_kin_1 += SPRODN(IMPULS,p,i,IMPULS,p,i) / MASSE(p,i);
#ifdef UNIAX
        E_rot_1 += SPRODN(DREH_IMPULS,p,i,DREH_IMPULS,p,i) / uniax_inert;
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
	fnorm += SPRODN(KRAFT,p,i,KRAFT,p,i);
        /* determine the biggest force component */
        tmp_f_max2 = MAX(SQR(KRAFT(p,i,X)),tmp_f_max2);
        tmp_f_max2 = MAX(SQR(KRAFT(p,i,Y)),tmp_f_max2);
#ifndef TWOD
        tmp_f_max2 = MAX(SQR(KRAFT(p,i,Z)),tmp_f_max2);
#endif
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
        dot = 2.0 * SPRODN(DREH_IMPULS,p,i,ACHSE,p,i);

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
        E_kin_2 += SPRODN(IMPULS,p,i,IMPULS,p,i) / MASSE(p,i);
#ifdef UNIAX
        E_rot_2 += SPRODN(DREH_IMPULS,p,i,DREH_IMPULS,p,i) / uniax_inert;
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

        ACHSE(p,i,X) += timestep * cross.x / uniax_inert;
        ACHSE(p,i,Y) += timestep * cross.y / uniax_inert;
        ACHSE(p,i,Z) += timestep * cross.z / uniax_inert;

        norm = SQRT( SPRODN(ACHSE,p,i,ACHSE,p,i) );

        ACHSE(p,i,X) /= norm;
        ACHSE(p,i,Y) /= norm;
        ACHSE(p,i,Z) /= norm;
#endif

#ifdef STRESS_TENS
        if (do_press_calc) {
          PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
          PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
          PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
          PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
          PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
          PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
	}
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

#ifdef FNORM
  MPI_Allreduce( &tmp_f_max2, &f_max2, 1, REAL, MPI_MAX, cpugrid);
#endif

#else
#ifdef FNORM
  f_max2 = tmp_f_max2;
#endif

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
  for (k=0; k<NCELLS; ++k) {

    int i;
    cell *p;
    real tmp;
    p = CELLPTR(k);

    /* loop over atoms in cell */
    for (i=0; i<p->n; ++i) {
      tmp = 1.0 / MASSE(p,i);
      dyn_stress_x += IMPULS(p,i,X) * IMPULS(p,i,X) * tmp;
      dyn_stress_y += IMPULS(p,i,Y) * IMPULS(p,i,Y) * tmp;
#ifndef TWOD
      dyn_stress_z += IMPULS(p,i,Z) * IMPULS(p,i,Z) * tmp;
#endif
#ifdef UNIAX
      Erot_old += SPRODN(DREH_IMPULS,p,i,DREH_IMPULS,p,i) / uniax_inert;
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
  real tmpvec1[5], tmpvec2[5], ttt;
  real reib, ireib;
  real tmp_f_max2=0.0;
  static real d_pressure;

  PxF     = 0.0;  
#ifdef RELAXINFO
  real tmp_x_max2=0.0;
  real tmppnorm, tmpfnorm,tmportx,tmporty,tmportz;
  xnorm   = 0.0;
  pnorm   = 0.0;
#endif
  
  if (steps == steps_min) {
    calc_dyn_pressure();
    if (isq_tau_xi==0.0) xi.x = 0.0;
  }

  fnorm    = 0.0;
  omega_E  = 0.0;
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
#pragma omp parallel for reduction(+:Ekin_new,Erot_new,fnorm,omega_E)
#endif
  for (k=0; k<NCELLS; ++k) {

    int i, j;
    cell *p;
    real tmp;
#ifdef UNIAX
    real dot, norm ;
    vektor cross ;
#endif
    p = CELLPTR(k);

#ifdef CLONE
    for (i=0; i<p->n; i+=nclones)
      for (j=1; j<nclones; j++) {
        KRAFT(p,i+j,X)  = KRAFT(p,i,X);
        KRAFT(p,i+j,Y)  = KRAFT(p,i,Y);
#ifndef TWOD
        KRAFT(p,i+j,Z)  = KRAFT(p,i,Z);
#endif
        IMPULS(p,i+j,X) = IMPULS(p,i,X);
        IMPULS(p,i+j,Y) = IMPULS(p,i,Y);
#ifndef TWOD
        IMPULS(p,i+j,Z) = IMPULS(p,i,Z);
#endif
      }
#endif /* CLONE */

    for (i=0; i<p->n; ++i) {

#ifdef FNORM
      fnorm   += SPRODN(KRAFT,p,i,KRAFT,p,i);
     /* determine the biggest force component */
      tmp_f_max2 = MAX(SQR(KRAFT(p,i,X)),tmp_f_max2);
      tmp_f_max2 = MAX(SQR(KRAFT(p,i,Y)),tmp_f_max2);
#ifndef TWOD
      tmp_f_max2 = MAX(SQR(KRAFT(p,i,Z)),tmp_f_max2);
#endif
#endif
#ifdef EINSTEIN
      omega_E += SPRODN(KRAFT,p,i,KRAFT,p,i) / MASSE(p,i);
#endif
#ifdef RELAXINFO
      PxF   += SPRODN(IMPULS,p,i,KRAFT,p,i)/MASSE(p,i);
      pnorm += SPRODN(IMPULS,p,i,IMPULS,p,i)/MASSE(p,i)/MASSE(p,i);
#endif

      /* new momenta */
      IMPULS(p,i,X) = (pfric*IMPULS(p,i,X)+timestep*KRAFT(p,i,X))*pifric;
      IMPULS(p,i,Y) = (pfric*IMPULS(p,i,Y)+timestep*KRAFT(p,i,Y))*pifric;
#ifndef TWOD
      IMPULS(p,i,Z) = (pfric*IMPULS(p,i,Z)+timestep*KRAFT(p,i,Z))*pifric;
#endif

#ifdef UNIAX
      /* new angular momenta */
      dot = 2.0 * SPRODN(DREH_IMPULS,p,i,ACHSE,p,i);
      DREH_IMPULS(p,i,X) = ireib * ( DREH_IMPULS(p,i,X) * reib 
              + timestep * DREH_MOMENT(p,i,X) - dot * ACHSE(p,i,X) );
      DREH_IMPULS(p,i,Y) = ireib * ( DREH_IMPULS(p,i,Y) * reib 
              + timestep * DREH_MOMENT(p,i,Y) - dot * ACHSE(p,i,Y) );
      DREH_IMPULS(p,i,Z) = ireib * ( DREH_IMPULS(p,i,Z) * reib 
              + timestep * DREH_MOMENT(p,i,Z) - dot * ACHSE(p,i,Z) );
#endif

      /* twice the new kinetic energy */ 
      Ekin_new += SPRODN(IMPULS,p,i,IMPULS,p,i) / MASSE(p,i);
#ifdef UNIAX
      Erot_new += SPRODN(DREH_IMPULS,p,i,DREH_IMPULS,p,i) / uniax_inert;
#endif

      /* new positions */
#ifdef RELAXINFO
      tmportx=ORT(p,i,X);
      tmporty=ORT(p,i,Y);
      tmportz=ORT(p,i,Z);
#endif
      tmp = timestep / MASSE(p,i);
      ORT(p,i,X) = (rfric * ORT(p,i,X) + IMPULS(p,i,X) * tmp) * rifric;
      ORT(p,i,Y) = (rfric * ORT(p,i,Y) + IMPULS(p,i,Y) * tmp) * rifric;
#ifndef TWOD
      ORT(p,i,Z) = (rfric * ORT(p,i,Z) + IMPULS(p,i,Z) * tmp) * rifric;
#endif

#ifdef RELAXINFO
      xnorm +=  SQR(ORT(p,i,X)-tmportx) + SQR(ORT(p,i,Y)-tmporty)+SQR(ORT(p,i,Z)-tmportz);
      /* determine the biggest displacement component */
      tmp_x_max2 =  MAX( SQR(ORT(p,i,X)-tmportx),tmp_x_max2);
      tmp_x_max2 =  MAX( SQR(ORT(p,i,Y)-tmporty),tmp_x_max2);
#ifndef TWOD
      tmp_x_max2 =  MAX( SQR(ORT(p,i,Z)-tmportz),tmp_x_max2);
#endif
#endif
      
#ifdef UNIAX
      /* new molecular axes */
      cross.x = DREH_IMPULS(p,i,Y) * ACHSE(p,i,Z)
              - DREH_IMPULS(p,i,Z) * ACHSE(p,i,Y);
      cross.y = DREH_IMPULS(p,i,Z) * ACHSE(p,i,X)
              - DREH_IMPULS(p,i,X) * ACHSE(p,i,Z);
      cross.z = DREH_IMPULS(p,i,X) * ACHSE(p,i,Y)
              - DREH_IMPULS(p,i,Y) * ACHSE(p,i,X);

      ACHSE(p,i,X) += timestep * cross.x / uniax_inert;
      ACHSE(p,i,Y) += timestep * cross.y / uniax_inert;
      ACHSE(p,i,Z) += timestep * cross.z / uniax_inert;

      norm = SQRT( SPRODN(ACHSE,p,i,ACHSE,p,i) );

      ACHSE(p,i,X) /= norm;
      ACHSE(p,i,Y) /= norm;
      ACHSE(p,i,Z) /= norm;
#endif

#ifdef STRESS_TENS
      if (do_press_calc) {
        PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
        PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
        PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
        PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
      }
#endif
    }
  }

#ifdef MPI
  /* add up results from all CPUs */
  tmpvec1[0] = Ekin_new;
  tmpvec1[1] = Erot_new;
  tmpvec1[2] = fnorm;
  tmpvec1[3] = omega_E;
  tmpvec1[4] = PxF;

  MPI_Allreduce( tmpvec1, tmpvec2, 5, REAL, MPI_SUM, cpugrid);

  Ekin_new = tmpvec2[0];
  Erot_new = tmpvec2[1];
  fnorm    = tmpvec2[2];
  omega_E  = tmpvec2[3];
  PxF      = tmpvec2[4];
#ifdef FNORM
  MPI_Allreduce( &tmp_f_max2, &f_max2, 1, REAL, MPI_MAX, cpugrid);
#endif
#ifdef RELAXINFO  
    MPI_Allreduce( &tmp_x_max2, &x_max2, 1, REAL, MPI_MAX, cpugrid);
#endif
#else
#ifdef FNORM
  f_max2 = tmp_f_max2;
#endif
#ifdef RELAXINFO  
  x_max2=tmp_x_max2;
#endif

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

  /* increment external pressure */
  if (steps == steps_min) {
    if (use_curr_pressure==1) {
      pressure_ext.x = pressure;
      use_curr_pressure = 0;
    }
    d_pressure = (pressure_end.x - pressure_ext.x) / (steps_max - steps_min);
  }
  pressure_ext.x += d_pressure;

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
  int k, sort;
  real tmp_f_max2=0.0;
  real Ekin_new = 0.0, ttt, tmpvec1[6], tmpvec2[6];
  vektor pfric, pifric, rfric, rifric, tvec;
  static vektor d_pressure;

  if (steps == steps_min) {
    calc_dyn_pressure();
    if (isq_tau_xi==0.0) {
      xi.x = 0.0;
      xi.y = 0.0;
#ifndef TWOD
      xi.z = 0.0;
#endif
    }
    xi.x *= relax_dirs.x;
    xi.y *= relax_dirs.y;
#ifndef TWOD
    xi.z *= relax_dirs.z;
#endif
  }

  fnorm    = 0.0;
  omega_E  = 0.0;
  stress_x = (dyn_stress_x + vir_xx) / volume;  dyn_stress_x = 0.0;
  stress_y = (dyn_stress_y + vir_yy) / volume;  dyn_stress_y = 0.0;
#ifndef TWOD
  stress_z = (dyn_stress_z + vir_zz) / volume;  dyn_stress_z = 0.0;
#endif

  /* time evolution of xi */
  ttt  = timestep * volume * isq_tau_xi / nactive;
  xi_old.x = xi.x;  xi.x += ttt * (stress_x - pressure_ext.x) * relax_dirs.x;
  xi_old.y = xi.y;  xi.y += ttt * (stress_y - pressure_ext.y) * relax_dirs.y;
#ifndef TWOD
  xi_old.z = xi.z;  xi.z += ttt * (stress_z - pressure_ext.z) * relax_dirs.z;
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
  for (k=0; k<NCELLS; ++k) {

    int i, sort;
    cell *p;
    real tmp;
    p = CELLPTR(k);

    /* loop over atoms in cell */
    for (i=0; i<p->n; ++i) {

      tmp = 1.0 / MASSE(p,i);
#ifdef FNORM
      fnorm   += SPRODN(KRAFT,p,i,KRAFT,p,i);
     /* determine the biggest force component */
      tmp_f_max2 = MAX(SQR(KRAFT(p,i,X)),tmp_f_max2);
      tmp_f_max2 = MAX(SQR(KRAFT(p,i,Y)),tmp_f_max2);
#ifndef TWOD
      tmp_f_max2 = MAX(SQR(KRAFT(p,i,Z)),tmp_f_max2);
#endif
#endif
#ifdef EINSTEIN
      omega_E += SPRODN(KRAFT,p,i,KRAFT,p,i) * tmp;
#endif
#ifdef STRESS_TENS
      if (do_press_calc) {
        PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) * tmp;
        PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) * tmp;
#ifndef TWOD
        PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) * tmp;
        PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) * tmp;
        PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) * tmp;
#endif
        PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) * tmp;
      }
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

    sort = VSORTE(p,i);

  	/* and set their force (->momentum) in restricted directions to 0 */
    IMPULS(p,i,X) *= (restrictions + sort)->x;
   	IMPULS(p,i,Y) *= (restrictions + sort)->y;
#ifndef TWOD
   	IMPULS(p,i,Z) *= (restrictions + sort)->z;
#endif


      /* new stress tensor (dynamic part only) */
      dyn_stress_x += IMPULS(p,i,X) * IMPULS(p,i,X) * tmp;
      dyn_stress_y += IMPULS(p,i,Y) * IMPULS(p,i,Y) * tmp;
#ifndef TWOD
      dyn_stress_z += IMPULS(p,i,Z) * IMPULS(p,i,Z) * tmp;
#endif

      /* twice the new kinetic energy */ 
      Ekin_new += SPRODN(IMPULS,p,i,IMPULS,p,i) * tmp;
	  
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
#ifdef FNORM
  MPI_Allreduce( &tmp_f_max2, &f_max2, 1, REAL, MPI_MAX, cpugrid);
#endif

#else
#ifdef FNORM
  f_max2 = tmp_f_max2;
#endif

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

  /* increment external pressure */
  if (steps == steps_min) {
    if (use_curr_pressure==1) {
      pressure_ext.x = stress_x;
      pressure_ext.y = stress_y;
#ifndef TWOD
      pressure_ext.z = stress_z;
#endif
      use_curr_pressure = 0;
    }
    d_pressure.x = (pressure_end.x-pressure_ext.x) / (steps_max-steps_min);
    d_pressure.y = (pressure_end.y-pressure_ext.y) / (steps_max-steps_min);
#ifndef TWOD
    d_pressure.z = (pressure_end.z-pressure_ext.z) / (steps_max-steps_min);
#endif
  }
  pressure_ext.x += d_pressure.x;
  pressure_ext.y += d_pressure.y;
#ifndef TWOD
  pressure_ext.z += d_pressure.z;
#endif

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
  int  k;
  real tmpvec1[6], tmpvec2[6], ttt;
  real tmp_f_max2=0.0;
  real E_kin_1        = 0.0, E_kin_2        = 0.0; 
  real E_kin_damp1    = 0.0, E_kin_damp2    = 0.0;
  real E_kin_stadium1 = 0.0, E_kin_stadium2 = 0.0;
  real reibung, reibung_y, eins_d_reib, eins_d_reib_y;
  real epsilontmp, eins_d_epsilontmp;

  real tmp, f; /* stadium function: the bath tub !!!!*/

  fnorm     = 0.0;
  sum_f     = 0.0;
  n_stadium = 0;

  if(expansionmode==1)
      dotepsilon = dotepsilon0 / (1.0 + dotepsilon0 * steps * timestep);

  /* loop over all cells */
#ifdef _OPENMP
#pragma omp parallel for reduction(+:E_kin_1,E_kin_2,E_kin_damp1,E_kin_damp2,E_kin_stadium1,E_kin_stadium2,sum_f,n_stadium,fnorm)
#endif
  for (k=0; k<NCELLS; ++k){ 

    int i, j, sort;
    cell *p;
    real tmp,tmp1,tmp2;

    p = CELLPTR(k);

#ifdef CLONE
    for (i=0; i<p->n; i+=nclones)
      for (j=1; j<nclones; j++) {
        KRAFT(p,i+j,X)  = KRAFT(p,i,X);
        KRAFT(p,i+j,Y)  = KRAFT(p,i,Y);
#ifndef TWOD
        KRAFT(p,i+j,Z)  = KRAFT(p,i,Z);
#endif
        IMPULS(p,i+j,X) = IMPULS(p,i,X);
        IMPULS(p,i+j,Y) = IMPULS(p,i,Y);
#ifndef TWOD
        IMPULS(p,i+j,Z) = IMPULS(p,i,Z);
#endif
      }
#endif /* CLONE */

    /* loop over all atoms in the cell */
    for (i=0; i<p->n; ++i) {
	
	/* if half axis in x-direction is zero: global viscous damping ! */
	if(stadium.x <= 0.0) { 
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
        tmp = SPRODN(IMPULS,p,i,IMPULS,p,i) / MASSE(p,i);
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
        fnorm   += SPRODN(KRAFT,p,i,KRAFT,p,i);
        /* determine the biggest force component */
        tmp_f_max2 = MAX(SQR(KRAFT(p,i,X)),tmp_f_max2);
        tmp_f_max2 = MAX(SQR(KRAFT(p,i,Y)),tmp_f_max2);
#ifndef TWOD
        tmp_f_max2 = MAX(SQR(KRAFT(p,i,Z)),tmp_f_max2);
#endif
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
        tmp = SPRODN(IMPULS,p,i,IMPULS,p,i) / MASSE(p,i);
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
        if (do_press_calc) {
          PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
          PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
          PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
          PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
          PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
          PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
	}
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
   real tmp_f_max2=0.0;
  static real *E_kin_1     = NULL; static real *E_kin_2     = NULL;
  static real *ftgtmpvec1  = NULL; static real *ftgtmpvec2  = NULL;
  static int  *iftgtmpvec1 = NULL; static int  *iftgtmpvec2 = NULL;

  real tmp,tmp1,tmp2;
  real ttt;
  real reibung, reibung_y, eins_d_reib, eins_d_reib_y;
  real epsilontmp, eins_d_epsilontmp;
  int slice;
  real gamma_tmp;

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

  fnorm = 0.0;

  if(expansionmode==1)
      dotepsilon = dotepsilon0 / (1.0 + dotepsilon0 * steps * timestep);
      
  /* loop over all cells */
  for (k=0; k<NCELLS; ++k) {

    int i, j, sort;
    cell *p;

    p = CELLPTR(k);

#ifdef CLONE
    for (i=0; i<p->n; i+=nclones)
      for (j=1; j<nclones; j++) {
        KRAFT(p,i+j,X)  = KRAFT(p,i,X);
        KRAFT(p,i+j,Y)  = KRAFT(p,i,Y);
#ifndef TWOD
        KRAFT(p,i+j,Z)  = KRAFT(p,i,Z);
#endif
        IMPULS(p,i+j,X) = IMPULS(p,i,X);
        IMPULS(p,i+j,Y) = IMPULS(p,i,Y);
#ifndef TWOD
        IMPULS(p,i+j,Z) = IMPULS(p,i,Z);
#endif
      }
#endif /* CLONE */

    /* loop over all atoms in cell */
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
	tmp  = SPRODN(IMPULS,p,i,IMPULS,p,i) / MASSE(p,i);
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
	gamma_tmp = (FABS(ORT(p,i,Y)-center.y) - stadium.y)/
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
      *(E_kin_1 + slice) += SPRODN(IMPULS,p,i,IMPULS,p,i) / MASSE(p,i);

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
        fnorm   += SPRODN(KRAFT,p,i,KRAFT,p,i);
        /* determine the biggest force component */
        tmp_f_max2 = MAX(SQR(KRAFT(p,i,X)),tmp_f_max2);
        tmp_f_max2 = MAX(SQR(KRAFT(p,i,Y)),tmp_f_max2);
#ifndef TWOD
        tmp_f_max2 = MAX(SQR(KRAFT(p,i,Z)),tmp_f_max2);
#endif
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
	*(E_kin_2 + slice) +=  SPRODN(IMPULS,p,i,IMPULS,p,i) / MASSE(p,i);
	
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
        if (do_press_calc) {
          PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
          PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
          PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
          PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
          PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
          PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
	}
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

/*****************************************************************************
*
*  Integrator with local temperature (Finnis)
*  
*
*****************************************************************************/

#ifdef FINNIS

void move_atoms_finnis(void)
{
  int j, k;

  real E_kin_1, E_kin_2;
  real tmp, tmp1, tmp2;
  real ttt;
  real temperature_at;
  real zeta_finnis;
  real tmp_f_max2=0.0;
  fnorm = 0.0;

  /* loop over all cells */
  for (k=0; k<NCELLS; ++k) {

    int i, j, sort;
    vektor *rest;
    cell *p;

    p = CELLPTR(k);

#ifdef CLONE
    for (i=0; i<p->n; i+=nclones)
      for (j=1; j<nclones; j++) {
        KRAFT(p,i+j,X)  = KRAFT(p,i,X);
        KRAFT(p,i+j,Y)  = KRAFT(p,i,Y);
#ifndef TWOD
        KRAFT(p,i+j,Z)  = KRAFT(p,i,Z);
#endif
        IMPULS(p,i+j,X) = IMPULS(p,i,X);
        IMPULS(p,i+j,Y) = IMPULS(p,i,Y);
#ifndef TWOD
        IMPULS(p,i+j,Z) = IMPULS(p,i,Z);
#endif
      }
#endif /* CLONE */

    /* loop over all atoms in the cell */
    for (i=0; i<p->n; ++i) {
	
      sort = VSORTE(p,i);
      rest = restrictions + sort;

      /* calc kinetic "temperature" for actual atom */
      tmp  = SPRODN(IMPULS,p,i,IMPULS,p,i) / MASSE(p,i);
#ifdef TWOD
      tmp2 = rest->x + rest->y;
#else
      tmp2 = rest->x + rest->y + rest->z;
#endif
      if (tmp2 != 0) tmp /= tmp2;
      /* to account for restricted mobilities and to avoid singularities */
      temperature_at = (tmp2 !=0) ? (tmp2/3.0 * temperature) : (1e-10); 

      /* to share the code with the non local version we overwrite 
	 the zeta values every timestep */
      zeta_finnis = zeta_0 * (tmp-temperature_at) 
	/ sqrt(SQR(tmp) + SQR(temperature_at*delta_finnis));     

      /* twice the old kinetic energy */
      E_kin_1 += SPRODN(IMPULS,p,i,IMPULS,p,i) / MASSE(p,i);

#ifdef FBC
      /* give virtual particles their extra force */
      KRAFT(p,i,X) += (fbc_forces + sort)->x;
      KRAFT(p,i,Y) += (fbc_forces + sort)->y;
#ifndef TWOD
      KRAFT(p,i,Z) += (fbc_forces + sort)->z;
#endif
#endif

      KRAFT(p,i,X) *= rest->x;
      KRAFT(p,i,Y) *= rest->y;
#ifndef TWOD
      KRAFT(p,i,Z) *= rest->z;
#endif

#ifdef FNORM
      fnorm   += SPRODN(KRAFT,p,i,KRAFT,p,i);
      /* determine the biggest force component */
      tmp_f_max2 = MAX(SQR(KRAFT(p,i,X)),tmp_f_max2);
      tmp_f_max2 = MAX(SQR(KRAFT(p,i,Y)),tmp_f_max2);
#ifndef TWOD
      tmp_f_max2 = MAX(SQR(KRAFT(p,i,Z)),tmp_f_max2);
#endif
#endif

      /* new momenta */
      IMPULS(p,i,X) += (-1.0*IMPULS(p,i,X) * zeta_finnis + KRAFT(p,i,X)) * timestep 
                       * rest->x;
      IMPULS(p,i,Y) += (-1.0*IMPULS(p,i,Y) * zeta_finnis + KRAFT(p,i,Y)) * timestep 
                       * rest->y;
#ifndef TWOD
      IMPULS(p,i,Z) += (-1.0*IMPULS(p,i,Z) * zeta_finnis + KRAFT(p,i,Z)) * timestep 
                       * rest->z;
#endif                  

      /* twice the new kinetic energy */ 
      E_kin_2 += SPRODN(IMPULS,p,i,IMPULS,p,i) / MASSE(p,i);

      /* new positions */
      tmp = timestep / MASSE(p,i);
      ORT(p,i,X) +=  tmp * IMPULS(p,i,X);
      ORT(p,i,Y) +=  tmp * IMPULS(p,i,Y);
#ifndef TWOD
      ORT(p,i,Z) +=  tmp * IMPULS(p,i,Z);
#endif

#ifdef STRESS_TENS
      if (do_press_calc) {
        PRESSTENS(p,i,xx) += IMPULS(p,i,X) * IMPULS(p,i,X) / MASSE(p,i);
        PRESSTENS(p,i,yy) += IMPULS(p,i,Y) * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
        PRESSTENS(p,i,zz) += IMPULS(p,i,Z) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,yz) += IMPULS(p,i,Y) * IMPULS(p,i,Z) / MASSE(p,i);
        PRESSTENS(p,i,zx) += IMPULS(p,i,Z) * IMPULS(p,i,X) / MASSE(p,i);
#endif
        PRESSTENS(p,i,xy) += IMPULS(p,i,X) * IMPULS(p,i,Y) / MASSE(p,i);
      }
#endif
    }
  }

  tot_kin_energy = ( E_kin_1 + E_kin_2 ) / 4.0;

#ifdef MPI
  /* add up results from different CPUs */
  tmp1 = tot_kin_energy;
  MPI_Allreduce( &tmp1, &tmp2, 1, REAL, MPI_SUM, cpugrid);
  tot_kin_energy = tmp2;
#endif

}

#else

void move_atoms_finnis(void) 
{
  if (myid==0)
  error("the chosen ensemble FINNIS is not supported by this binary");
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
  real tmp_f_max2=0.0;
  /* we handle 2 ensembles ensindex = 0 -> NVT ;ensindex = 1 -> NVE */
  int ensindex = 0;
  real kin_energie_1[2] = {0.0,0.0}, kin_energie_2[2] = {0.0,0.0};
  real tmpvec1[5], tmpvec2[5], ttt;
  n_stadium = 0;

  /* loop over all atoms */
#ifdef _OPENMP
#pragma omp parallel for reduction(+:kin_energie_1[0],kin_energie_1[1],kin_energie_2[0],kin_energie_2[2],n_stadium)
#endif
  for (k=0; k<NCELLS; ++k) {

    int i;
    cell *p;
    real reibung, eins_d_reib;
    real tmp;
    vektor d;
    int sort=0;

    p = CELLPTR(k);

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
        kin_energie_1[ensindex] += SPRODN(IMPULS,p,i,IMPULS,p,i) / MASSE(p,i);

        /* new momenta */
	sort = VSORTE(p,i);
	IMPULS(p,i,X) = (IMPULS(p,i,X)*reibung + timestep * KRAFT(p,i,X))
                          * eins_d_reib * (restrictions + sort)->x;
        IMPULS(p,i,Y) = (IMPULS(p,i,Y)*reibung + timestep * KRAFT(p,i,Y))
                          * eins_d_reib * (restrictions + sort)->y;

        /* twice the new kinetic energy */ 
        kin_energie_2[ensindex] += SPRODN(IMPULS,p,i,IMPULS,p,i) / MASSE(p,i);

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
  E_kin_stadium      = tmpvec2[2];
  kin_energie_2[1]   = tmpvec2[3];
  n_stadium          = (int)tmpvec2[4];
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
*  NVX Integrator for heat conductivity (direct method)
* 
******************************************************************************/

#ifdef NVX

void move_atoms_nvx(void)

{
  int  k, num, nhalf;  
  real Ekin_1, Ekin_2, Ekin_right, Ekin_left, delta_E, Evec1[3], Evec2[3];
  real scale, xx, rescale_left, rescale_right;

  Evec1[0] = 0.0;
  Evec1[1] = 0.0;
  Evec1[2] = 0.0;

  nhalf   = hc_nlayers / 2;
  scale   = hc_nlayers / box_x.x;
  delta_E = hc_heatcurr * 2 * box_y.y * box_z.z * timestep;

  /* loop over all atoms */
  for (k=0; k<NCELLS; ++k) {

    int i;
    cell *p = CELLPTR(k);

    for (i=0; i<p->n; ++i) {

      /* twice the old kinetic energy */ 
      Ekin_1 = SPRODN(IMPULS,p,i,IMPULS,p,i) / MASSE(p,i);

      /* new momenta */
      IMPULS(p,i,X) += timestep * KRAFT(p,i,X); 
      IMPULS(p,i,Y) += timestep * KRAFT(p,i,Y); 
#ifndef TWOD
      IMPULS(p,i,Z) += timestep * KRAFT(p,i,Z); 
#endif

      /* twice the new kinetic energy */ 
      Ekin_2 = SPRODN(IMPULS,p,i,IMPULS,p,i) / MASSE(p,i);

      /* new positions */
      ORT(p,i,X) += timestep * IMPULS(p,i,X) / MASSE(p,i);
      ORT(p,i,Y) += timestep * IMPULS(p,i,Y) / MASSE(p,i);
#ifndef TWOD
      ORT(p,i,Z) += timestep * IMPULS(p,i,Z) / MASSE(p,i);
#endif
      Evec1[0] += (Ekin_1 + Ekin_2) / 4.0;

      /* kinetic energy of layers 0 and nhalf */
      xx = ORT(p,i,X);
      if (xx<0.0) xx += box_x.x;
      num = (int) (scale * xx);
      if      (num >= hc_nlayers) num -= hc_nlayers;
      if      (num == 0    ) Evec1[1] += Ekin_2;
      else if (num == nhalf) Evec1[2] += Ekin_2;
    }
  }

#ifdef MPI
  /* Add up results from all cpus */
  MPI_Allreduce( Evec1, Evec2, 3, REAL, MPI_SUM, cpugrid);
  tot_kin_energy = Evec2[0];
  Ekin_left      = Evec2[1];
  Ekin_right     = Evec2[2];
#else
  tot_kin_energy = Evec1[0];
  Ekin_left      = Evec1[1];
  Ekin_right     = Evec1[2];
#endif

  /* rescale factors for momenta */
  rescale_left  = SQRT( 1.0 - delta_E / Ekin_left  );
  rescale_right = SQRT( 1.0 + delta_E / Ekin_right );

  /* rescale the momenta */
  for (k=0; k<NCELLS; ++k) {

    int i;
    cell *p = CELLPTR(k);

    for (i=0; i<p->n; ++i) {
	    
      /* which layer? */
      xx = ORT(p,i,X);
      if (xx<0.0) xx += box_x.x;
      num = (int) (scale * xx);
      if (num >= hc_nlayers) num -= hc_nlayers;

      /* rescale momenta */
      if (num == 0) {
        IMPULS(p,i,X) *= rescale_left;
        IMPULS(p,i,Y) *= rescale_left;
#ifndef TWOD
        IMPULS(p,i,Z) *= rescale_left;
#endif
      } else if (num == nhalf) {
        IMPULS(p,i,X) *= rescale_right;
        IMPULS(p,i,Y) *= rescale_right;
#ifndef TWOD
        IMPULS(p,i,Z) *= rescale_right;
#endif
      }
    }
  }
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
* Move the atoms for the Conjugated Gradient relaxator
*
*****************************************************************************/

#ifdef CG 

void move_atoms_cg(real alpha)
{
  int k;
  real tmp_x_max2 = 0.0;
  real tmpvec1[1], tmpvec2[1];
  /* loop over all cells */
  xnorm=0;
  for (k=0; k<NCELLS; ++k) {

    int  i, j, sort;
    cell *p;

    p = CELLPTR(k);

#ifdef CLONE
    for (i=0; i<p->n; i+=nclones)
      for (j=1; j<nclones; j++) {
        KRAFT(p,i+j,X)  = KRAFT(p,i,X);
        KRAFT(p,i+j,Y)  = KRAFT(p,i,Y);
#ifndef TWOD
        KRAFT(p,i+j,Z)  = KRAFT(p,i,Z);
#endif
      }
#endif /* CLONE */

    for (i=0; i<p->n; ++i) {
      /* CG:  move atoms in search direction for linmin */
      ORT(p,i,X) = OLD_ORT(p,i,X) + alpha * CG_H(p,i,X);
      ORT(p,i,Y) = OLD_ORT(p,i,Y) + alpha * CG_H(p,i,Y);
#ifndef TWOD
      ORT(p,i,Z) = OLD_ORT(p,i,Z) + alpha * CG_H(p,i,Z);
#endif

#ifdef RELAXINFO
      xnorm   += alpha * alpha *  SPRODN(CG_H,p,i,CG_H,p,i);
      /* determine the biggest force component */
      tmp_x_max2 = MAX( alpha * alpha *SQR(CG_H(p,i,X)),tmp_x_max2);
      tmp_x_max2 = MAX( alpha * alpha *SQR(CG_H(p,i,Y)),tmp_x_max2);
#ifndef TWOD
      tmp_x_max2 = MAX( alpha * alpha *SQR(CG_H(p,i,Z)),tmp_x_max2);
#endif
#endif
    }
  }

#ifdef RELAXINFO
#ifdef MPI
  tmpvec1[0] = xnorm;
  MPI_Allreduce( tmpvec1, tmpvec2, 1, REAL, MPI_SUM, cpugrid); 
  xnorm  = tmpvec2[0];
  MPI_Allreduce( &tmp_x_max2, &x_max2, 1, REAL, MPI_MAX, cpugrid);
#else
  x_max2 = tmp_x_max2;
#endif
#endif

  if ((cg_infolevel>0) && (0==myid)) {
    printf("moveatoms, alpha %.12e , xmax %.12e\n", alpha, SQRT(x_max2));
    fflush(stdout);
  }

}

#else

void move_atoms_cg(real alpha) 
{
  if (myid==0)
    error("the chosen ensemble CG is not supported by this binary");
}

#endif


/* steepest descent step, needed for ACG , could also be used just 
   to do a steepest descent */

#if defined(CG) || defined(SD) 

void move_atoms_sd(real alpha)
{
  int k;
  real tmp_x_max2 = 0.0;
  real tmpvec1[1], tmpvec2[1];
  /* loop over all cells */
  xnorm=0;
  for (k=0; k<NCELLS; ++k) {

    int  i, j, sort;
    cell *p;

    p = CELLPTR(k);

#ifdef CLONE
    for (i=0; i<p->n; i+=nclones)
      for (j=1; j<nclones; j++) {
        KRAFT(p,i+j,X)  = KRAFT(p,i,X);
        KRAFT(p,i+j,Y)  = KRAFT(p,i,Y);
#ifndef TWOD
        KRAFT(p,i+j,Z)  = KRAFT(p,i,Z);
#endif
      }
#endif /* CLONE */

    for (i=0; i<p->n; ++i) {
      /* CG:  move atoms in force direction for linmin */
      ORT(p,i,X) = OLD_ORT(p,i,X) + alpha * KRAFT(p,i,X);
      ORT(p,i,Y) = OLD_ORT(p,i,Y) + alpha * KRAFT(p,i,Y);
#ifndef TWOD
      ORT(p,i,Z) = OLD_ORT(p,i,Z) + alpha * KRAFT(p,i,Z);
#endif

#ifdef RELAXINFO
      xnorm   += alpha * alpha *  SPRODN(KRAFT,p,i,KRAFT,p,i);
     /* determine the biggest force component */
      tmp_x_max2 = MAX( alpha * alpha *SQR(KRAFT(p,i,X)),tmp_x_max2);
      tmp_x_max2 = MAX( alpha * alpha *SQR(KRAFT(p,i,Y)),tmp_x_max2);
#ifndef TWOD
      tmp_x_max2 = MAX( alpha * alpha *SQR(KRAFT(p,i,Z)),tmp_x_max2);
#endif

#endif

    }
  }

#ifdef RELAXINFO
#ifdef MPI
  tmpvec1[0] = xnorm;
  MPI_Allreduce( tmpvec1, tmpvec2, 1, REAL, MPI_SUM, cpugrid); 
  xnorm          = tmpvec2[0];
  MPI_Allreduce( &tmp_x_max2, &x_max2, 1, REAL, MPI_MAX, cpugrid);
#else
  x_max2 = tmp_x_max2;
#endif
#endif

  if ((cg_infolevel>0) && (0==myid)) {
    printf("moveatoms, alpha %.12e , xmax %.12e\n", alpha, SQRT(x_max2) );
    fflush(stdout);
  }

}

#else

void move_atoms_sd(real alpha) 
{
  if (myid==0)
    error("the chosen ensemble CG or SD is not supported by this binary");
}

#endif

#ifdef SHOCK

/*****************************************************************************
*
* Calculate average momentum in x direction
*
*****************************************************************************/

void calc_pxavg(void) 
{
  integer *num_1, *num_2;
  real *dat_1, *dat_2, scale;
  int i, k;

  /* backup if dist_ur is not set */
  if (0.0==dist_ur.x) {
    dist_ur.x = box_x.x;
    dist_ur.y = box_y.y;
#ifndef TWOD
    dist_ur.z = box_z.z;
#endif
  }

#ifdef MPI
  dat_1 = (real    *) malloc( dist_dim.x * sizeof(real  ) );
  num_1 = (integer *) malloc( dist_dim.x * sizeof(integer) );
  dat_2 = (real    *) malloc( dist_dim.x * sizeof(real  ) );
  num_2 = (integer *) malloc( dist_dim.x * sizeof(integer) );
  if ((NULL==dat_1) || (NULL==num_1) || (NULL==dat_2) || (NULL==num_2))
    error("Cannot allocate distribution data.");
#else
  dat_1 = (real    *) malloc( dist_dim.x * sizeof(real  ) );
  num_1 = (integer *) malloc( dist_dim.x * sizeof(integer) );
  dat_2 = dat_1;
  num_2 = num_1;
  if ((NULL==dat_1) || (NULL==num_1))
    error("Cannot allocate distribution data.");
#endif

  /* the bins are orthogonal slices in space */
  scale = dist_dim.x / (dist_ur.x - dist_ll.x);

  /* clear distributions */
  for (i=0; i<dist_dim.x; i++) {
    dat_1[i] = 0.0;
    num_1[i] = 0;
  }

  /* loop over all atoms */
  for (k=0; k<NCELLS; ++k) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      int n = (int) (scale * (ORT(p,i,X) - dist_ll.x));
      if ((n < 0) || (n >= dist_dim.x)) continue;
      num_1[n]++;
      dat_1[n] += IMPULS(p,i,X);
    }
  }

#ifdef MPI  /* add up results form different CPUs */
  MPI_Allreduce( dat_1, dat_2, dist_dim.x, REAL,    MPI_SUM, cpugrid);
  MPI_Allreduce( num_1, num_2, dist_dim.x, INTEGER, MPI_SUM, cpugrid);
#endif

  /* normalize distribution */
  for (i=0; i<dist_dim.x; i++) {
    if (num_2[i] > 0) dat_2[i] /= num_2[i];
  }

  /* loop over all atoms */
  for (k=0; k<NCELLS; ++k) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      int n = (int) (scale * (ORT(p,i,X) - dist_ll.x) - 0.5);
      if (n < 0) {
        PXAVG(p,i) = dat_2[0];
      } 
      else if (n >= dist_dim.x-1) {
        PXAVG(p,i) = dat_2[dist_dim.x-1];
      }
      else if (num_2[n]>0){
        real chi = (ORT(p,i,X) - n / scale) * scale; 
        PXAVG(p,i) = dat_2[n] * (1-chi) + chi * dat_2[n+1];
      }
      else { 
	PXAVG(p,i) = dat_2[n];
      }
    }
  }

#ifdef MPI
  free(dat_1);
  free(num_1);
  free(dat_2);
  free(num_2);
#else
  free(dat_1);
  free(num_1);
#endif

}

#endif
