/******************************************************************************
*
* imd_integrate -- various md integrators
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"


/*****************************************************************************
*
* Basic NVE Verlet Integrator
*
*****************************************************************************/

#ifdef NVE

void move_atoms_nve(void)

{
  cell *p;
  int i;
  real kin_energie_1,kin_energie_2;
  int r,s,t;  
  real tmp = 0;
  int iend;

  tot_kin_energy = 0;

  /* loop over all atoms */
  for ( r = cellmin.x; r < cellmax.x; ++r )
    for ( s = cellmin.y; s < cellmax.y; ++s )
#ifndef TWOD
      for ( t = cellmin.z; t < cellmax.z; ++t )
#endif
	{

#ifndef TWOD
	p = PTR_3D_V(cell_array, r, s, t, cell_dim);
#else
	p = PTR_2D_V(cell_array, r, s,    cell_dim);
#endif
	iend = p->n;

#ifdef PVPCRAY
#pragma ivdep
#endif
#ifdef SX4
#pragma vdir vector,nodep
#endif
	for (i = 0;i < iend; ++i) {

	  kin_energie_1 =  SPRODN(p->impuls,i,p->impuls,i);

	  /* new momenta */
	  p->impuls X(i) += timestep * p->kraft X(i);
	  p->impuls Y(i) += timestep * p->kraft Y(i);
#ifndef TWOD
	  p->impuls Z(i) += timestep * p->kraft Z(i);
#endif

	  kin_energie_2 =  SPRODN(p->impuls,i,p->impuls,i);

          /* sum up kinetic energy on this CPU */
          tmp += (kin_energie_1 + kin_energie_2) / ( 4 * MASSE(p,i) );

          /* new positions - do not move atoms with negative numbers */
	  if (NUMMER(p,i) >= 0) {
	    p->ort X(i) += timestep * p->impuls X(i) / MASSE(p,i);
	    p->ort Y(i) += timestep * p->impuls Y(i) / MASSE(p,i);
#ifndef TWOD
            p->ort Z(i) += timestep * p->impuls Z(i) / MASSE(p,i);
#endif
	  } else {
	    p->impuls X(i) = 0.0;
	    p->impuls Y(i) = 0.0;
#ifndef TWOD
	    p->impuls Z(i) = 0.0;
#endif
          }

#ifdef STRESS_TENS
#ifdef SHOCK
	  if (shock_mode == 2) {
                if ( p->ort X(i) < box_x.x*0.5 )
		  p->presstens X(i) += ( p->impuls X(i) - shock_speed * MASSE(p,i)) * (p->impuls X(i) - shock_speed *MASSE(p,i))/ MASSE(p,i);
		else
		  p->presstens X(i) += ( p->impuls X(i) + shock_speed * MASSE(p,i)) * (p->impuls X(i) + shock_speed *MASSE(p,i))/ MASSE(p,i);
	      }
#else
          p->presstens X(i) += p->impuls X(i) * p->impuls X(i) / MASSE(p,i);
#endif
          p->presstens Y(i) += p->impuls Y(i) * p->impuls Y(i) / MASSE(p,i);
          p-> presstens_offdia Z(i) += p->impuls X(i) * p->impuls Y(i) /  MASSE(p,i);
#ifndef TWOD
          p->presstens Z(i) += p->impuls Z(i) * p->impuls Z(i) / MASSE(p,i);
          p-> presstens_offdia X(i) += p->impuls Y(i) * p->impuls Z(i) /  MASSE(p,i);
          p-> presstens_offdia Y(i) += p->impuls Z(i) * p->impuls X(i) /  MASSE(p,i);
#endif
#endif
        }
      }

#ifdef MPI
  /* Add kinetic energy form all cpus */
  MPI_Allreduce( &tmp, &tot_kin_energy, 1, MPI_REAL, MPI_SUM, cpugrid);
#else
  tot_kin_energy = tmp;
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
* NVE Verlet Integrator with microconvergence relaxation 
*
*****************************************************************************/

#ifdef MIK

void move_atoms_mik(void)

{
  cell *p;
  int i,iend;
  real kin_energie_1,kin_energie_2;
  int r,s,t;  
  real tmp = 0.0; 
  real tmp2;

  tot_kin_energy = 0.0;

  /* loop over all atoms */
  for ( r = cellmin.x; r < cellmax.x; ++r )
    for ( s = cellmin.y; s < cellmax.y; ++s )
#ifndef TWOD
      for ( t = cellmin.z; t < cellmax.z; ++t )
#endif
	{

#ifndef TWOD
	p = PTR_3D_V(cell_array, r, s, t, cell_dim);
#else
	p = PTR_2D_V(cell_array, r, s,    cell_dim);
#endif
	iend = p->n;

#ifdef PVPCRAY
#pragma ivdep
#endif
#ifdef SX4
#pragma vdir vector,nodep
#endif
	for (i = 0;i < iend; ++i) {

	  kin_energie_1 =  SPRODN(p->impuls,i,p->impuls,i);

	  /* Neue Impulse */
	  p->impuls X(i) += timestep * p->kraft X(i);
	  p->impuls Y(i) += timestep * p->kraft Y(i);
#ifndef TWOD
	  p->impuls Z(i) += timestep * p->kraft Z(i);
#endif

	  /* Mikroconvergence Algorithm - set velocity zero if a*v < 0 */
	  /* Do not move atoms in strip or with negative numbers */
	  if ((0 > SPRODN(p->impuls,i,p->kraft,i)) || (NUMMER(p,i)<0)) {
	    p->impuls X(i) = 0.0;
	    p->impuls Y(i) = 0.0;
#ifndef TWOD
	    p->impuls Z(i) = 0.0;
#endif
          } else { /* new positions */
            tmp2 = timestep / MASSE(p,i);
            p->ort X(i) += tmp2 * p->impuls X(i);
            p->ort Y(i) += tmp2 * p->impuls Y(i);
#ifndef TWOD
            p->ort Z(i) += tmp2 * p->impuls Z(i);
#endif
 	  }

	  kin_energie_2 =  SPRODN(p->impuls,i,p->impuls,i);

       	  /* sum up kinetic energy on this CPU */ 
          tmp += (kin_energie_1 + kin_energie_2) / ( 4.0 * MASSE(p,i) );

#ifdef STRESS_TENS
          p-> presstens X(i) += p->impuls X(i) * p->impuls X(i) / MASSE(p,i);
          p-> presstens Y(i) += p->impuls Y(i) * p->impuls Y(i) / MASSE(p,i);
          p-> presstens_offdia Z(i) += p->impuls X(i) * p->impuls Y(i) /  MASSE(p,i);
#ifndef TWOD
          p-> presstens Z(i) += p->impuls Z(i) * p->impuls Z(i) / MASSE(p,i);
          p-> presstens_offdia X(i) += p->impuls Y(i) * p->impuls Z(i) /  MASSE(p,i);
          p-> presstens_offdia Y(i) += p->impuls Z(i) * p->impuls X(i) /  MASSE(p,i);
#endif
#endif
	}
      }

#ifdef MPI
  /* Add kinetic energy form all cpus */
  MPI_Allreduce( &tmp, &tot_kin_energy, 1, MPI_REAL, MPI_SUM, cpugrid);
#else  
  tot_kin_energy = tmp;
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
* Minimizer:
*
* A sort of ... Steepest Descent - JH 1999: 02/03/...
*
*****************************************************************************/

#ifdef MSD

void move_atoms_msd(void)

{
  cell *p;
  int i,iend;
  real tmp_msd  = 0.0;
  real tmp_step = 0.0;
  int r,s,t;  
  real tmp = 0.0; 				/* tot_kin_energy = 0 ! */

  /* loop over all atoms */
  for ( r = cellmin.x; r < cellmax.x; ++r )
    for ( s = cellmin.y; s < cellmax.y; ++s )
#ifndef TWOD
      for ( t = cellmin.z; t < cellmax.z; ++t )
#endif
        {

#ifndef TWOD
        p = PTR_3D_V(cell_array, r, s, t, cell_dim);
#else
        p = PTR_2D_V(cell_array, r, s,    cell_dim);
#endif

        tmp_step = timestep; 
        iend = p->n;
        for (i = 0;i < iend; ++i) {
          
          /* new positions - do not move atoms with negative numbers */
          if (NUMMER(p,i) >= 0) {
            tmp_msd = tmp_step / sqrt(SPRODN(p->kraft,i,p->kraft,i)) ;
            p->ort X(i) += tmp_msd * p->kraft X(i);
            p->ort Y(i) += tmp_msd * p->kraft Y(i);
#ifndef TWOD
            p->ort Z(i) += tmp_msd * p->kraft Z(i);
#endif 
          }

          /* new momenta */
          p->impuls X(i) = 0.0;
          p->impuls Y(i) = 0.0;
#ifndef TWOD
          p->impuls Z(i) = 0.0;
#endif
        }
      }

#ifdef MPI
  /* Add kinetic energy form all cpus */
  MPI_Allreduce( &tmp, &tot_kin_energy, 1, MPI_REAL, MPI_SUM, cpugrid);
#else  
  tot_kin_energy = tmp;
#endif

}

#else  

void move_atoms_msd(void) 
{
  if (myid==0)
  error("the chosen ensemble MSD is not supported by this binary");
}

#endif 

/*****************************************************************************
*
* NVT Verlet Integrator with Nose Hoover Thermostat
*
*****************************************************************************/

#ifdef NVT

void move_atoms_nvt(void)

{
  cell *p;
  int i;
  real kin_energie_1,kin_energie_2;
  int r,s,t;  
  real reibung;
  real eins_d_reibung;
  real tmp,tmp2;

  tot_kin_energy = 0;
  kin_energie_1  = 0;
  kin_energie_2  = 0;
  reibung        =      1 - eta * timestep / 2.0;
  eins_d_reibung = 1 / (1 + eta * timestep / 2.0);
   
  /* loop over all atoms */
  for ( r = cellmin.x; r < cellmax.x; ++r )
    for ( s = cellmin.y; s < cellmax.y; ++s )
#ifndef TWOD
      for ( t = cellmin.z; t < cellmax.z; ++t )
#endif
      {
        
#ifndef TWOD
	p = PTR_3D_V(cell_array, r, s, t, cell_dim);
#else
	p = PTR_2D_V(cell_array, r, s,    cell_dim);
#endif

	for (i = 0;i < p->n; ++i) {

          /* twice the old kinetic energy */
	  kin_energie_1 +=  SPRODN(p->impuls,i,p->impuls,i) / MASSE(p,i);

	  /* new momenta */
	  p->impuls X(i) = (p->impuls X(i)*reibung + timestep * p->kraft X(i)) 
	                   * eins_d_reibung;
	  p->impuls Y(i) = (p->impuls Y(i)*reibung + timestep * p->kraft Y(i)) 
	                   * eins_d_reibung;
#ifndef TWOD
	  p->impuls Z(i) = (p->impuls Z(i)*reibung + timestep * p->kraft Z(i)) 
	                   * eins_d_reibung;
#endif

	  /* twice the new kinetic energy */ 
	  kin_energie_2  +=  SPRODN(p->impuls,i,p->impuls,i) / MASSE(p,i);
	  
          /* new positions - do not move atoms with negative numbers */
	  if (NUMMER(p,i) >= 0) {
	    p->ort X(i) += timestep * p->impuls X(i) / MASSE(p,i);
	    p->ort Y(i) += timestep * p->impuls Y(i) / MASSE(p,i);
#ifndef TWOD
	    p->ort Z(i) += timestep * p->impuls Z(i) / MASSE(p,i);
#endif
	  } else {
	    p->impuls X(i) = 0.0;
	    p->impuls Y(i) = 0.0;
#ifndef TWOD
	    p->impuls Z(i) = 0.0;
#endif
	  }

#ifdef STRESS_TENS
          p->presstens X(i) += p->impuls X(i) * p->impuls X(i) / MASSE(p,i);
          p->presstens Y(i) += p->impuls Y(i) * p->impuls Y(i) / MASSE(p,i);
          p-> presstens_offdia Z(i) += p->impuls X(i) * p->impuls Y(i) /  MASSE(p,i);
#ifndef TWOD
          p->presstens Z(i) += p->impuls Z(i) * p->impuls Z(i) / MASSE(p,i);
          p-> presstens_offdia X(i) += p->impuls Y(i) * p->impuls Z(i) /  MASSE(p,i);
          p-> presstens_offdia Y(i) += p->impuls Z(i) * p->impuls X(i) /  MASSE(p,i);
#endif
#endif
        }
      }
  
  tot_kin_energy = (kin_energie_1 + kin_energie_2) / 4.0;

#ifdef MPI
  /* add kinetic energy from all cpus */
  MPI_Allreduce( &tot_kin_energy, &tmp,  1, MPI_REAL, MPI_SUM, cpugrid);
  MPI_Allreduce( &kin_energie_2,  &tmp2, 1, MPI_REAL, MPI_SUM, cpugrid);

  tot_kin_energy = tmp;
  kin_energie_2  = tmp2;
#endif

  /* update parameters */
  tmp  = DIM * natoms * temperature;
  eta += timestep * (kin_energie_2 / tmp - 1.0) * isq_tau_eta;
  
}

#else

void move_atoms_nvt(void) 
{
  if (myid==0)
  error("the chosen ensemble NVT is not supported by this binary");
}

#endif


/******************************************************************************
*
* NPT Verlet Integrator with Nose Hoover Thermostat
*
******************************************************************************/

#ifdef NPT_iso

void move_atoms_npt_iso(void)

{
  cell *p;
  int  i, r, s, t;
  real Ekin_old, Ekin_new;
  real tmp, tmp1, tmp2;
  real fric, ifric;

  Ekin_old         = 0.0;
  Ekin_new         = 0.0;
  box_size.x      += 2.0 * timestep * xi.x;  /* relative box size change */  
  actual_shrink.x *= box_size.x;
  fric             = 1.0 - (xi.x + eta) * timestep / 2.0;
  ifric            = 1.0 / ( 1.0 + (xi.x + eta) * timestep / 2.0 );

  /* loop over all cells */
  for ( r = cellmin.x; r < cellmax.x; ++r )
    for ( s = cellmin.y; s < cellmax.y; ++s )
#ifndef TWOD
      for ( t = cellmin.z; t < cellmax.z; ++t )
#endif
      {

#ifndef TWOD
	p = PTR_3D_V(cell_array, r, s, t, cell_dim);
#else
	p = PTR_2D_V(cell_array, r, s,    cell_dim);
#endif

        /* loop over atoms in cell */
	for (i = 0;i < p->n; ++i) {
      
	  /* twice the old kinetic energy */ 
	  Ekin_old += SPRODN(p->impuls,i,p->impuls,i) / MASSE(p,i);

	  /* new momenta */
	  p->impuls X(i) = (fric*p->impuls X(i)+timestep*p->kraft X(i))*ifric;
	  p->impuls Y(i) = (fric*p->impuls Y(i)+timestep*p->kraft Y(i))*ifric;
#ifndef TWOD
	  p->impuls Z(i) = (fric*p->impuls Z(i)+timestep*p->kraft Z(i))*ifric;
#endif

	  /* twice the new kinetic energy */ 
	  Ekin_new += SPRODN(p->impuls,i,p->impuls,i) / MASSE(p,i);
	  
	  /* new positions */
          tmp = p->impuls X(i) * (1.0 + box_size.x) / (2.0 * MASSE(p,i));
	  p->ort X(i) = box_size.x * ( p->ort X(i) + timestep * tmp );
          tmp = p->impuls Y(i) * (1.0 + box_size.x) / (2.0 * MASSE(p,i));
	  p->ort Y(i) = box_size.x * ( p->ort Y(i) + timestep * tmp );
#ifndef TWOD
          tmp = p->impuls Z(i) * (1.0 + box_size.x) / (2.0 * MASSE(p,i));
	  p->ort Z(i) = box_size.x * ( p->ort Z(i) + timestep * tmp );
#endif

#ifdef STRESS_TENS
          p-> presstens X(i) += p->impuls X(i) * p->impuls X(i) / MASSE(p,i);
          p-> presstens Y(i) += p->impuls Y(i) * p->impuls Y(i) / MASSE(p,i);
          p-> presstens_offdia Z(i) += p->impuls X(i) * p->impuls Y(i) /  MASSE(p,i);
#ifndef TWOD
          p-> presstens Z(i) += p->impuls Z(i) * p->impuls Z(i) / MASSE(p,i);
          p-> presstens_offdia X(i) += p->impuls Y(i) * p->impuls Z(i) /  MASSE(p,i);
          p-> presstens_offdia Y(i) += p->impuls Z(i) * p->impuls X(i) /  MASSE(p,i);
#endif
#endif
	}
      }

#ifdef MPI
  /* add data from all cpus */
  MPI_Allreduce( &Ekin_old, &tmp,  1, MPI_REAL, MPI_SUM, cpugrid);
  MPI_Allreduce( &Ekin_new, &tmp1, 1, MPI_REAL, MPI_SUM, cpugrid);

  Ekin_old = tmp;
  Ekin_new = tmp1;
#endif

  tot_kin_energy  = (Ekin_old + Ekin_new) / 4.0;
  pressure        = (2.0 * tot_kin_energy + virial) / (DIM * volume);

  /* udate parameters */
  tmp  = DIM * natoms * temperature;
  eta += timestep * (Ekin_new / tmp - 1.0) * isq_tau_eta;

  tmp = xi_old.x + timestep * 2.0 * (pressure - pressure_ext.x) * volume
                          * isq_tau_xi / (natoms * temperature);
  xi_old.x = xi.x;
  xi.x = tmp;

  /* new box size */
  box_x.x *= box_size.x;  ibox_x.x /= box_size.x;  tbox_x.x /= box_size.x;
  box_x.y *= box_size.x;  ibox_x.x /= box_size.x;  tbox_x.x /= box_size.x;
  box_y.x *= box_size.x;  ibox_y.x /= box_size.x;  tbox_y.x /= box_size.x;
  box_y.y *= box_size.x;  ibox_y.x /= box_size.x;  tbox_y.x /= box_size.x;
#ifndef TWOD
  box_x.z *= box_size.x;  ibox_x.z /= box_size.x;  tbox_x.z /= box_size.x;
  box_y.z *= box_size.x;  ibox_y.z /= box_size.x;  tbox_y.z /= box_size.x;
  box_z.x *= box_size.x;  ibox_z.x /= box_size.x;  tbox_z.x /= box_size.x;
  box_z.y *= box_size.x;  ibox_z.y /= box_size.x;  tbox_z.y /= box_size.x;
  box_z.z *= box_size.x;  ibox_z.z /= box_size.x;  tbox_z.z /= box_size.x;
#endif  

  /* old box_size relative to the current one, which is set to 1.0 */
  box_size.x = 1.0 / box_size.x;

  /* check whether the box has not changed too much */
  if (actual_shrink.x < limit_shrink.x) {
     cells_too_small = 1;
  }
  if ((actual_shrink.x < limit_shrink.x) 
       || (actual_shrink.x > limit_growth.x)) {
     revise_cell_division = 1;
  }

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
* NPT Verlet Integrator with Nose Hoover Thermostat
*
******************************************************************************/

#ifdef NPT_axial

void move_atoms_npt_axial(void)

{
  cell *p;
  int i, r, s, t;
  real Ekin, tmp, xi_tmp;
  vektor fric, ifric, tmpvec;

  /* initialize data, and compute new box size */
  Ekin             = 0.0;
  stress.x         = 0.0;
  stress.y         = 0.0;
  box_size.x      += 2.0 * timestep * xi.x;  /* relative box size change */  
  box_size.y      += 2.0 * timestep * xi.y;  /* relative box size change */  
  actual_shrink.x *= box_size.x;
  actual_shrink.y *= box_size.y;
  fric.x           = 1.0 - (xi.x + eta) * timestep / 2.0;
  fric.y           = 1.0 - (xi.y + eta) * timestep / 2.0;
  ifric.x          = 1.0 / ( 1.0 + (xi.x + eta) * timestep / 2.0 );
  ifric.y          = 1.0 / ( 1.0 + (xi.y + eta) * timestep / 2.0 );
#ifndef TWOD
  stress.z         = 0.0;
  box_size.z      += 2.0 * timestep * xi.z;  /* relative box size change */  
  actual_shrink.z *= box_size.z;
  fric.z           = 1.0 - (xi.z + eta) * timestep / 2.0;
  ifric.z          = 1.0 / ( 1.0 + (xi.z + eta) * timestep / 2.0 );
#endif

  /* loop over all cells */
  for ( r = cellmin.x; r < cellmax.x; ++r )
    for ( s = cellmin.y; s < cellmax.y; ++s )
#ifndef TWOD
      for ( t = cellmin.z; t < cellmax.z; ++t )
#endif
	{
#ifndef TWOD
	p = PTR_3D_V(cell_array, r, s, t, cell_dim);
#else
	p = PTR_2D_V(cell_array, r, s,    cell_dim);
#endif

        /* loop over atoms in cell */
	for (i = 0;i < p->n; ++i) {

          /* contribution of old p's to stress */
          stress.x += p->impuls X(i) * p->impuls X(i) / MASSE(p,i);
          stress.y += p->impuls Y(i) * p->impuls Y(i) / MASSE(p,i);
#ifndef TWOD
          stress.z += p->impuls Z(i) * p->impuls Z(i) / MASSE(p,i);
#endif

	  /* new momenta */
	  p->impuls X(i) = (fric.x * p->impuls X(i)
                                   + timestep * p->kraft X(i)) * ifric.x;
	  p->impuls Y(i) = (fric.y * p->impuls Y(i)
                                   + timestep * p->kraft Y(i)) * ifric.y;
#ifndef TWOD
	  p->impuls Z(i) = (fric.z * p->impuls Z(i)
                                   + timestep * p->kraft Z(i)) * ifric.z;
#endif

	  /* new kinetic energy, and contribution of new p's to stress */
          tmp       = p->impuls X(i) * p->impuls X(i) / MASSE(p,i);
          stress.x += tmp;
          Ekin     += tmp;
          tmp       = p->impuls Y(i) * p->impuls Y(i) / MASSE(p,i);
          stress.y += tmp;
          Ekin     += tmp;
#ifndef TWOD
          tmp       = p->impuls Z(i) * p->impuls Z(i) / MASSE(p,i);
          stress.z += tmp;
          Ekin     += tmp;
#endif
	  
	  /* new positions */
          tmp = p->impuls X(i) * (1.0 + box_size.x) / (2.0 * MASSE(p,i));
	  p->ort X(i) = box_size.x * (p->ort X(i) + timestep * tmp);
          tmp = p->impuls Y(i) * (1.0 + box_size.y) / (2.0 * MASSE(p,i));
	  p->ort Y(i) = box_size.y * (p->ort Y(i) + timestep * tmp);
#ifndef TWOD
          tmp = p->impuls Z(i) * (1.0 + box_size.z) / (2.0 * MASSE(p,i));
	  p->ort Z(i) = box_size.z * (p->ort Z(i) + timestep * tmp);
#endif

#ifdef STRESS_TENS
          p-> presstens X(i) += p->impuls X(i) * p->impuls X(i) / MASSE(p,i);
          p-> presstens Y(i) += p->impuls Y(i) * p->impuls Y(i) / MASSE(p,i);
          p-> presstens_offdia Z(i) += p->impuls X(i) * p->impuls Y(i) /  MASSE(p,i);
#ifndef TWOD
          p-> presstens Z(i) += p->impuls Z(i) * p->impuls Z(i) / MASSE(p,i);
          p-> presstens_offdia X(i) += p->impuls Y(i) * p->impuls Z(i) /  MASSE(p,i);
          p-> presstens_offdia Y(i) += p->impuls Z(i) * p->impuls X(i) /  MASSE(p,i);
#endif
#endif

	}
      }

#ifdef MPI
  /* Add kinetic energy form all cpus */
  MPI_Allreduce( &Ekin,   &tmp,      1, MPI_REAL, MPI_SUM, cpugrid);
  MPI_Allreduce( &stress, &tmpvec, DIM, MPI_REAL, MPI_SUM, cpugrid);

  Ekin     = tmp;
  stress.x = tmpvec.x;
  stress.y = tmpvec.y;
#ifndef TWOD
  stress.z = tmpvec.z;
#endif
#endif

  tot_kin_energy  = stress.x + stress.y;
  stress.x        = (stress.x / 2.0 + vir_vect.x) / volume;
  stress.y        = (stress.y / 2.0 + vir_vect.y) / volume;
#ifndef TWOD
  tot_kin_energy += stress.z;
  stress.z        = (stress.z / 2.0 + vir_vect.z) / volume;
#endif
  tot_kin_energy /= 4.0;

  /* update parameters */
  tmp  = DIM * natoms * temperature;
  eta += timestep * (Ekin / tmp - 1.0) * isq_tau_eta;

  tmp  = timestep * 2.0 * volume * isq_tau_xi / (natoms * temperature);
  
  xi_tmp   = xi_old.x + tmp * (stress.x - pressure_ext.x);
  xi_old.x = xi.x;
  xi.x     = xi_tmp;
  xi_tmp   = xi_old.y + tmp * (stress.y - pressure_ext.y);
  xi_old.y = xi.y;
  xi.y     = xi_tmp;
#ifndef TWOD
  xi_tmp   = xi_old.z + tmp * (stress.z - pressure_ext.z);
  xi_old.z = xi.z;
  xi.z     = xi_tmp;
#endif

  /* new box size (box is rectangular) */
  box_x.x *= box_size.x;  ibox_x.x /= box_size.x;  tbox_x.x /= box_size.x;
  box_y.y *= box_size.y;  ibox_y.y /= box_size.y;  tbox_y.y /= box_size.y;
#ifndef TWOD
  box_z.z *= box_size.z;  ibox_z.z /= box_size.z;  tbox_z.z /= box_size.z;
#endif  

  /* old box size relative to new one */
  box_size.x  = 1.0 / box_size.x;
  box_size.y  = 1.0 / box_size.y;
#ifndef TWOD
  box_size.z  = 1.0 / box_size.z;
#endif  

  /* check whether box has not changed too much */
  if  ((actual_shrink.x < limit_shrink.x) 
    || (actual_shrink.y < limit_shrink.y) 
#ifndef TWOD
    || (actual_shrink.z < limit_shrink.z) 
#endif
     ) {
     cells_too_small = 1;
  };
  if  ((actual_shrink.x < limit_shrink.x) || (actual_shrink.x > limit_growth.x)
    || (actual_shrink.y < limit_shrink.y) || (actual_shrink.y > limit_growth.y)
#ifndef TWOD
    || (actual_shrink.z < limit_shrink.z) || (actual_shrink.z > limit_growth.z)
#endif
     ) {
     revise_cell_division = 1;
  };

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
* NVT Verlet Integrator with Andersen Thermostat for NVT 
*
*****************************************************************************/

#ifdef AND

void move_atoms_and(void)

{
  cell *p;
  int i;
  real kin_energie_1,kin_energie_2;
  int r,s,t;  
  real tmp = 0;
  static int count = 0;

  tot_kin_energy = 0;

  /* loop over all atoms */
  for ( r = cellmin.x; r < cellmax.x; ++r )
    for ( s = cellmin.y; s < cellmax.y; ++s )
#ifndef TWOD
      for ( t = cellmin.z; t < cellmax.z; ++t )
#endif
	{

#ifndef TWOD
	p = PTR_3D_V(cell_array, r, s, t, cell_dim);
#else
	p = PTR_2D_V(cell_array, r, s,    cell_dim);
#endif

	for (i = 0;i < p->n; ++i) {

	  kin_energie_1 =  SPRODN(p->impuls,i,p->impuls,i);

	  /* new momenta */
	  p->impuls X(i) += timestep * p->kraft X(i);
	  p->impuls Y(i) += timestep * p->kraft Y(i);
#ifndef TWOD
	  p->impuls Z(i) += timestep * p->kraft Z(i);
#endif

	  kin_energie_2 =  SPRODN(p->impuls,i,p->impuls,i);

	  /* sum up kinetic energy */ 
	  tot_kin_energy += (kin_energie_1 + kin_energie_2) / (4 * MASSE(p,i));
	  
          /* new positions - do not move atoms  with negative numbers */
	  if (NUMMER(p,i) >= 0) {
	    p->ort X(i) += timestep * p->impuls X(i) / MASSE(p,i);
	    p->ort Y(i) += timestep * p->impuls Y(i) / MASSE(p,i);
#ifndef TWOD
	    p->ort Z(i) += timestep * p->impuls Z(i) / MASSE(p,i);
#endif
	  } else {
	    p->impuls X(i) = 0.0;
	    p->impuls Y(i) = 0.0;
#ifndef TWOD
	    p->impuls Z(i) = 0.0;
#endif
	  }

#ifdef STRESS_TENS
          p-> presstens X(i) += p->impuls X(i) * p->impuls X(i) / MASSE(p,i);
          p-> presstens Y(i) += p->impuls Y(i) * p->impuls Y(i) / MASSE(p,i);
          p-> presstens_offdia Z(i) += p->impuls X(i) * p->impuls Y(i) /  MASSE(p,i);
#ifndef TWOD
          p-> presstens Z(i) += p->impuls Z(i) * p->impuls Z(i) / MASSE(p,i);
          p-> presstens_offdia X(i) += p->impuls Y(i) * p->impuls Z(i) /  MASSE(p,i);
          p-> presstens_offdia Y(i) += p->impuls Z(i) * p->impuls X(i) /  MASSE(p,i);
#endif
#endif
	}
      }

#ifdef MPI
  /* Add kinetic energy for all CPUs */
  MPI_Allreduce( &tot_kin_energy, &tmp, 1, MPI_REAL, MPI_SUM, cpugrid);
  tot_kin_energy = tmp;
#endif

  /* Andersen Thermostat -- Initialize the velocities now and then */
  ++count;
  if ((tmp_interval!=0) && (0==count%tmp_interval)) maxwell(temperature);

}

#else

void move_atoms_and(void) 
{
  if (myid==0)
  error("the chosen ensemble AND is not supported by this binary");
}

#endif


/*****************************************************************************
*
*  Monte Carlo "integrator" 
*
*****************************************************************************/

#ifdef MC

void move_atoms_mc(void)
{
   one_mc_step();
}

#else

void move_atoms_mc(void)
{
  if (myid==0)
  error("the chosen ensemble MC is not supported by this binary");
}

#endif


/*****************************************************************************
*
* Verlet Integrator with stadion damping and fixed borders 
* for fracture studies 
*
*****************************************************************************/

#ifdef FRAC

void move_atoms_frac(void)

{
  cell *p;
  int i;
  real kin_energie_1,kin_energie_2;
  int r,s,t;  
  vektor stad, c_halbe;
  real gamma;
  real tmp = 0;
  static int count = 0;
  tot_kin_energy = 0;

  stad.x = 2 / stadion.x;
  stad.y = 2 / stadion.y;

  c_halbe.x = box_x.x / 2.0;
  c_halbe.y = box_y.y / 2.0;

  /* loop over all atoms */
  for ( r = cellmin.x; r < cellmax.x; ++r )
    for ( s = cellmin.y; s < cellmax.y; ++s )
#ifndef TWOD
      for ( t = cellmin.z; t < cellmax.z; ++t )
#endif
	{

#ifndef TWOD
	p = PTR_3D_V(cell_array, r, s, t, cell_dim);
#else
	p = PTR_2D_V(cell_array, r, s,    cell_dim);
#endif

	for (i = 0;i < p->n; ++i) {

          kin_energie_1 =  SPRODN(p->impuls,i,p->impuls,i);

          /* Calculate stadium function */
          gamma = gamma_bar * 
          ( SQR( SQR( (p->ort X(i) - c_halbe.x) / 2*stadion.x)
                 + SQR( (p->ort Y(i) - c_halbe.y) / 2*stadion.y)) -
            2* ( SQR( (p->ort X(i) - c_halbe.x) / 2*stadion.x)
                 - SQR( (p->ort Y(i) - c_halbe.y) / 2*stadion.y))
                    + 1);

              if (gamma_cut <= gamma) gamma = gamma_cut;
              if ((SQR((p->ort X(i) - c_halbe.x)  * stad.x)
                  + SQR((p->ort Y(i) - c_halbe.y) * stad.y)) < 1.0) gamma=0;
                  
	  /* new momenta */
	  p->impuls X(i) = 
	    (p->kraft X(i) * timestep 
	      + p->impuls X(i) * (1-(1/(2 * MASSE(p,i))) * gamma * timestep))
	      / (1 + (1/(2 * MASSE(p,i))) * gamma * timestep);
	  p->impuls Y(i) = 
	    (p->kraft Y(i) * timestep 
	      + p->impuls Y(i) * (1-(1/(2 * MASSE(p,i))) * gamma * timestep))
	      / (1 + (1/(2 * MASSE(p,i))) * gamma * timestep);
#ifndef TWOD
	  p->impuls Z(i) = 
	    (p->kraft Z(i) * timestep 
	      + p->impuls Z(i) * (1-(1/(2 * MASSE(p,i))) * gamma * timestep))
	      / (1 + (1/(2 * MASSE(p,i))) * gamma * timestep);
#endif

	  kin_energie_2 =  SPRODN(p->impuls,i,p->impuls,i);

	  /* sum up kinetic energy */ 
	  tot_kin_energy += (kin_energie_1 + kin_energie_2) / (4*MASSE(p,i));
	  
	  /* new positions - do not move atoms with negative numbers */
	  if (NUMMER(p,i) >= 0) {
	    p->ort X(i) += timestep * p->impuls X(i) / MASSE(p,i);
	    p->ort Y(i) += timestep * p->impuls Y(i) / MASSE(p,i);
#ifndef TWOD
	    p->ort Z(i) += timestep * p->impuls Z(i) / MASSE(p,i);
#endif
	  } else {
	    p->impuls X(i) = 0.0;
	    p->impuls Y(i) = 0.0;
#ifndef TWOD
	    p->impuls Z(i) = 0.0;
#endif
	  }

#ifdef STRESS_TENS
          p-> presstens X(i) += p->impuls X(i) * p->impuls X(i) / MASSE(p,i);
          p-> presstens Y(i) += p->impuls Y(i) * p->impuls Y(i) / MASSE(p,i);
          p-> presstens_offdia Z(i) += p->impuls X(i) * p->impuls Y(i) /  MASSE(p,i);
#ifndef TWOD
          p-> presstens Z(i) += p->impuls Z(i) * p->impuls Z(i) / MASSE(p,i);
          p-> presstens_offdia X(i) += p->impuls Y(i) * p->impuls Z(i) /  MASSE(p,i);
          p-> presstens_offdia Y(i) += p->impuls Z(i) * p->impuls X(i) /  MASSE(p,i);
#endif
#endif
	}
      }

#ifdef MPI
  /* Add kinetic energy for all CPUs */
  MPI_Allreduce( &tot_kin_energy, &tmp, 1, MPI_REAL, MPI_SUM, cpugrid);
  tot_kin_energy = tmp;
#endif

}

#else

void move_atoms_frac(void) 
{
  if (myid==0)
  error("the chosen ensemble FRAC is not supported by this binary");
}

#endif


/******************************************************************************
*
*  NVX  Integrator for heat conductivity
* 
******************************************************************************/

#ifdef NVX

void move_atoms_nvx(void)

{
  cell *p;
  int	i;
  real  Ekin_1, Ekin_2;
  real  Ekin_left = 0.0, Ekin_right = 0.0;
  integer natoms_left = 0, natoms_right = 0;
  real  px, vol, real_tmp;
  integer num, nhalf, r, s, t, int_tmp;  
  real	scale, rescale, Rescale;
  vektor tot_impuls_left, tot_impuls_right, vectmp;
  real inv_mass_left=0.0, inv_mass_right=0.0;
 
  heat_cond = 0.0;
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
  for ( r = cellmin.x; r < cellmax.x; ++r )
    for ( s = cellmin.y; s < cellmax.y; ++s )
#ifndef TWOD
      for ( t = cellmin.z; t < cellmax.z; ++t )
#endif
      {
#ifdef TWOD
        p = PTR_2D_V(cell_array, r, s,    cell_dim);
#else
        p = PTR_3D_V(cell_array, r, s, t, cell_dim);
#endif
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
  MPI_Allreduce( &tot_kin_energy, &real_tmp, 1, MPI_REAL, MPI_SUM, cpugrid);
  tot_kin_energy                 = real_tmp;
  MPI_Allreduce( &heat_cond,      &real_tmp, 1, MPI_REAL, MPI_SUM, cpugrid);
  heat_cond                      = real_tmp;
  MPI_Allreduce( &Ekin_left,      &real_tmp, 1, MPI_REAL, MPI_SUM, cpugrid);
  Ekin_left                      = real_tmp;
  MPI_Allreduce( &Ekin_right,     &real_tmp, 1, MPI_REAL, MPI_SUM, cpugrid);
  Ekin_right                     = real_tmp;
  MPI_Allreduce( &inv_mass_left,  &real_tmp, 1, MPI_REAL, MPI_SUM, cpugrid);
  inv_mass_left                  = real_tmp;
  MPI_Allreduce( &inv_mass_right, &real_tmp, 1, MPI_REAL, MPI_SUM, cpugrid);
  inv_mass_right                 = real_tmp;
  MPI_Allreduce( tot_impuls_left,  vectmp, DIM, MPI_REAL, MPI_SUM, cpugrid);
  tot_impuls_left                = vectmp;
  MPI_Allreduce( tot_impuls_right, vectmp, DIM, MPI_REAL, MPI_SUM, cpugrid);
  tot_impuls_right               = vectmp;
  MPI_Allreduce( &natoms_left,    &int_tmp,  1, INTEGER,  MPI_SUM, cpugrid);
  natoms_left                    = int_tmp;
  MPI_Allreduce( &natoms_right,   &int_tmp,  1, INTEGER,  MPI_SUM, cpugrid);
  natoms_right                   = int_tmp;
#endif

#ifdef TWOD
  vol = box_x.x * box_y.y           * (tran_nlayers-2) / tran_nlayers;
#else
  vol = box_x.x * box_y.y * box_z.z * (tran_nlayers-2) / tran_nlayers;
#endif
  heat_cond /= (2 * vol * (tran_Tleft -  tran_Tright));

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

  for ( r = cellmin.x; r < cellmax.x; ++r )
    for ( s = cellmin.y; s < cellmax.y; ++s )
#ifndef TWOD
      for ( t = cellmin.z; t < cellmax.z; ++t )
#endif
        {
#ifdef TWOD
	  p = PTR_2D_V(cell_array, r, s,    cell_dim);
#else
	  p = PTR_3D_V(cell_array, r, s, t, cell_dim);
#endif
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
	    p-> presstens X(i) += p->impuls X(i) * p->impuls X(i) / MASSE(p,i);
	    p-> presstens Y(i) += p->impuls Y(i) * p->impuls Y(i) / MASSE(p,i);
            p-> presstens_offdia Z(i) += p->impuls X(i) * p->impuls Y(i) /  MASSE(p,i);
#ifndef TWOD
	    p-> presstens Z(i) += p->impuls Z(i) * p->impuls Z(i) / MASSE(p,i);
            p-> presstens_offdia X(i) += p->impuls Y(i) * p->impuls Z(i) /  MASSE(p,i);
            p-> presstens_offdia Y(i) += p->impuls Z(i) * p->impuls X(i) /  MASSE(p,i);
#endif
#endif
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
