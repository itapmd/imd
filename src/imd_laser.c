/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2008 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_laser.c -- Laser heating of sample by different rescaling algorithms
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/
#include <assert.h>
#include "imd.h"


void init_laser()
{
  laser_p_peak=laser_mu*laser_sigma_e/laser_sigma_t/sqrt(2*M_PI);
  laser_sigma_t_squared=laser_sigma_t*laser_sigma_t;

  if (0==myid) {
      printf( "Parameter laser_rescale_mode is %d\n", laser_rescale_mode );
      printf( "Parameter laser_delta_temp is  %1.10f\n", laser_delta_temp );
#ifndef TWOD
      printf( "Laser irradiates from direction (%d, %d, %d)\n", laser_dir.x,
	      laser_dir.y, laser_dir.z);
#else
      printf( "Laser irradiates from direction (%d, %d)\n", laser_dir.x,
	      laser_dir.y);
#endif /*TWOD*/

      if (laser_mu==0.0) {
	printf( "Absorption length is infinite.\n" );
      } else {
        printf( "Absorption length is %1.10f\n", 1.0/laser_mu );
      }
      printf( "Laser energy density is %1.10f\n", laser_sigma_e);
      printf( "Laser pulse duration (sigma) is %1.10f\n", laser_sigma_t);
      printf( "Time t_0 of laser pulse is %1.10f\n",laser_t_0);
      printf( "(%1.10f time steps after start of simulation)\n", 
             laser_t_0/timestep);
  }
}

void laser_rescale_dummy()
{
  /* do nothing */
}

/* Routine to generate random 3D unit vector */
void rand_uvec_3d(real* x, real* y, real *z)
{

  real cosph, sinph, theta, costh, sinth;

  cosph = 1.0 - 2.0 * drand48();
  sinph = sqrt( 1.0 - cosph * cosph );
  theta = 2.0 * M_PI * drand48();
  costh = cos(theta);
  sinth = sin(theta);

  *x = cosph;
  *y = costh * sinph;
  *z = sinth * sinph;
}


/* Routine to generate random 2D unit vector */
void rand_uvec_2d(real* x, real* y)
{
  real phi;

  phi   = drand48() * 2 * M_PI;
  *x = cos(phi);
  *y = sin(phi);
}

real laser_calc_depth( cell *p, int i )
{
  real depth;
  depth = laser_dir.x * ORT(p,i,X)
        + laser_dir.y * ORT(p,i,Y)
#ifndef TWOD
        + laser_dir.z * ORT(p,i,Z)
#endif
        ;
  depth -= laser_offset;

  if (depth < 0.0) depth=0.0; /* we don't want negative depths */

  return depth;
}


void laser_rescale_1()
{ /* Rescale all atom velocities so that exactly the same amount of
     kinetic energy is added to the atoms at a given depth.
     The direction of the momentum vectors remains unchanged, except for
     zero-velocity atoms, which will be given a random direction.*/

  real exp_gauss_time_etc, gauss_time_squared;
  int k;
  real cosph, sinph, theta, costh, sinth;

  gauss_time_squared = timestep * steps - laser_t_0;
  gauss_time_squared *= gauss_time_squared;
  exp_gauss_time_etc = exp(-gauss_time_squared/laser_sigma_t_squared/2.0)
                       * laser_p_peak * timestep * laser_atom_vol;

  for (k=0; k<NCELLS; k++) {
    cell *p;
    p = CELLPTR(k);
    int i;
    for (i=0; i < p->n; i++) {
      real p_0_square; /* square of initial atom momentum */
      real de; /* Kinetic energy to be added to the atom */
      real depth; /* Distance from origin along laser_dir */
      real tmpx, tmpy, tmpsqr;
#ifndef TWOD
      real tmpz;
#endif
      real scale_p; /* Scale factor for atom momentum */

      p_0_square = SPRODN(IMPULS,p,i,IMPULS,p,i);
      depth = laser_calc_depth(p,i);

      de = exp(-laser_mu*depth) * exp_gauss_time_etc;

      if ( p_0_square == 0.0 ) { /* we need a direction for the momentum. */
#ifndef TWOD
        /* find random 3d unit vector */
        rand_uvec_3d(&tmpx, &tmpy, &tmpz);
#else
        /* find random 2d unit vector */
        rand_uvec_2d(&tmpx, &tmpy);      
#endif
        scale_p = sqrt( de * 2.0 * MASSE(p,i) );
        IMPULS(p,i,X) = tmpx * scale_p;
        IMPULS(p,i,Y) = tmpy * scale_p;
#ifndef TWOD
        IMPULS(p,i,Z) = tmpz * scale_p;
#endif
      } else { /* rescale present momentum */

        scale_p = sqrt( de * 2.0 * MASSE(p,i) / p_0_square + 1.0 );
        IMPULS(p,i,X) *= scale_p;
        IMPULS(p,i,Y) *= scale_p;
#ifndef TWOD
        IMPULS(p,i,Z) *= scale_p;
#endif

      }
    }
  }
} /* void laser_rescale_1()*/


void laser_rescale_2()
/* Instead of just rescaling, add the momentum increment in random direction.
   Only thereafter rescale so the right amount of energy gets absorbed.
   Involves several square roots, probably slower
*/ 
{

  real exp_gauss_time_etc, gauss_time_squared;
  int k;

  gauss_time_squared = timestep * steps - laser_t_0;
  gauss_time_squared *= gauss_time_squared;
  exp_gauss_time_etc = exp(-gauss_time_squared/laser_sigma_t_squared)
                       * laser_p_peak * timestep * laser_atom_vol;

  for (k=0; k<NCELLS; k++) {
    cell *p;
    p = CELLPTR(k);
    int i;
    for (i=0; i < p->n; i++) {
      real p_0_square; /* square of initial atom momentum */
      real p_0;
      real de; /* Kinetic energy to be added to the atom */
      real dp;
      real depth; /* Distance from origin along laser_dir */
      real tmpx, tmpy, tmpsqr;
#ifndef TWOD
      real tmpz;
#endif
      real scale_p; /* Scale factor for atom momentum */

      p_0_square = SPRODN(IMPULS,p,i,IMPULS,p,i);
      p_0=sqrt(p_0_square);

      depth = laser_calc_depth(p,i);

      de = exp(-laser_mu*depth) * exp_gauss_time_etc;
    
      dp = sqrt( p_0_square + 2*de*MASSE(p,i) ) - p_0;
    

      /* Momentum increment is to point in a random direction... */
#ifndef TWOD
      rand_uvec_3d(&tmpx, &tmpy, &tmpz);
#else
      rand_uvec_2d(&tmpx, &tmpy);
#endif
      IMPULS(p,i,X) += tmpx * dp;
      IMPULS(p,i,Y) += tmpy * dp;
#ifndef TWOD
      IMPULS(p,i,Z) += tmpz * dp;
#endif

      /* rescale present momentum to an absolute value of dp + p_0 */
      scale_p = (p_0 + dp) / SQRT( SPRODN(IMPULS,p,i,IMPULS,p,i) );
      IMPULS(p,i,X) *= scale_p;
      IMPULS(p,i,Y) *= scale_p;
#ifndef TWOD
      IMPULS(p,i,Z) *= scale_p;
#endif
    
    }
  }
} /* void laser_rescale_2()*/


void laser_rescale_3()
/* 3rd Method: 1. Rescale
               2. Add a small random component orthogonal to resulting vector
*/
{

  /* doesn't do anything yet, will probably be removed soon */

} /* void laser_rescale_3()*/

#ifdef TTM
real ttm_calc_depth( int i, int j, int k )
{
  real depth;
  depth = ( fd_h.x * ((i-1) + my_coord.x*(local_fd_dim.x-2)) ) * laser_dir.x
        + ( fd_h.y * ((j-1) + my_coord.y*(local_fd_dim.y-2)) ) * laser_dir.y
        + ( fd_h.z * ((k-1) + my_coord.z*(local_fd_dim.z-2)) ) * laser_dir.z;
#ifdef DEBUG
  assert(depth >= 0.0); /* we don't want negative depths here */
#endif
  depth -= laser_offset;
  return (depth<0)?0:depth;
}

/* This is laser_rescale_mode == 4 */
void laser_rescale_ttm()
{
  /* This function just writes the exponential source term into the FD net.
   * Heating of the electrons occurs later,
   * in the numerical solution of the pdeq.                              */

  int i,j,k;
  real exp_gauss_time_etc, gauss_time_squared, depth;
  
  gauss_time_squared = timestep * steps - laser_t_0;
  gauss_time_squared *= gauss_time_squared;
  exp_gauss_time_etc = exp(-gauss_time_squared/laser_sigma_t_squared/2.0)
    * laser_p_peak ;/*  * fd_h.x*fd_h.y*fd_h.z not needed */

  /* loop over all FD cells, excluding ghost layers */
  
  for(i=1; i<local_fd_dim.x-1; ++i)
  {
    for(j=1; j<local_fd_dim.y-1; ++j)
    {
      for(k=1; k<local_fd_dim.z-1; ++k)
      {
	depth = ttm_calc_depth(i,j,k); 
        l2[i][j][k].source = l1[i][j][k].source
	                   = exp(-laser_mu*depth) * exp_gauss_time_etc; 
      }
    }
  }
}
#endif


