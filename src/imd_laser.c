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

/* Various laser profiles */
#ifdef LASERYZ
#include "imd_laser_profiles.c"
#endif

double calc_laser_atom_vol(double deltax, int leftside, int rightside, int* xdens_1)
{
  /* The volume per atom (laser_atom_vol) is calculated. 
     
     We calculate this, because in huge samples, the volume per atom inside the
     irradiated volume can differ from the volume per atom of the whole sample, i.e.
     for the second laser pulse, the surface already starts melting, and a different
     local volume per atom value is obvious */
  
  int l=0; 
  int totald = 0;              /* total number of atoms inside the irradiated volume */

  double intensityfrac = 0.01; /* to which fraction the intensity should decay */
  double xpenetrate = 0.0;     /* penetration depth of laser at which intensity is at obove percentage */

  

#ifdef LASERYZ
#ifdef DEBUG
  for (l=leftside; l<= rightside; l++)
  {
    totald += xdens_1[l];
  }
  
  double voltemp = ((double)rightside*deltax-((double)leftside+0.5)*deltax)*box_y.y*box_z.z/(double)totald;
  printf("laser_atom_vol: %f , leftside: %i rightside: %i leftside_x: %f, totaldens: %i \n", voltemp,  leftside, rightside, ((double)(leftside)+0.5)*deltax, totald);
  
  totald = 0;
#endif
#endif
 xpenetrate = -log(intensityfrac)/laser_mu;

  if ( ( xpenetrate + ((double)(leftside) + 0.5)*deltax ) < (double)rightside*deltax)
  {
    
     /* In this case, the intesity decays to less than 1% INSIDE the irradiated volume, 
	we have to find  the 'new' rightside dcell number to calculate the proper number 
	density */
	
    /* calculates the dcell-number of the most inner irradiated volume */
    
    rightside = (int) (( xpenetrate + ((double)(leftside) + 0.5)*deltax )/deltax);
    
    /* add up the atoms inside the irradiated region */
    
    for (l=leftside; l<= rightside; l++)
    {
      totald += xdens_1[l];
    }

#ifdef LASERYZ
#ifdef DEBUG    
    printf("laser_atom_vol: %f , leftside: %i rightside: %i leftside_x: %f, totaldens: %i \n", xpenetrate*box_y.y*box_z.z/(double)totald,  leftside, rightside, ((double)(leftside)+0.5)*deltax, totald);
#endif
#endif
    
    /* return the new value for laser_atom_vol */
    return xpenetrate*box_y.y*box_z.z/(double)totald;
  
  }
  
  else
    
  {
    /* In this case, the intensity decays to less than 1% OUTSIDE the actual volume, so 
       we can computer the number density over the total volume, which is from leftside 
       to the rightside */
    
 
    /* add up the atoms inside the irradiated region */
    for (l=leftside; l<= rightside; l++)
    {
      totald += xdens_1[l];
    }
    
#ifdef LASERYZ
#ifdef DEBUG    
    printf("laser_atom_vol: %f , leftside: %i rightside: %i leftside_x: %f, totaldens: %i \n",((double)rightside-((double)leftside+0.5))*deltax*box_y.y*box_z.z/(double)totald ,  leftside, rightside, ((double)(leftside)+0.5)*deltax, totald);
#endif
#endif
    
    /*return the new value for laser_atom_vol */
    return (double)(rightside-leftside-0.5)*deltax*box_y.y*box_z.z/(double)totald;
    
  }

}


double get_surface()
{
  /* Calculates a 1D (anlong the x-axes) density distribution, to 'find'
     the surface. 1D is sufficient here, because laser radiation is only
     availible in (1,0,0) direction (yet) */
  
  /* One density-cell is chosen to be ~ 2.5 A, but I depends on the system*/
  
  double deltax = 2.5;
  int ndcells = (int)(box_x.x/deltax);
  int *xdens_1, *xdens_2 = NULL;
  
  int l,k;
  int rightside = ndcells;
  int leftside = 0;

#ifdef MPI2
  MPI_Alloc_mem( ndcells * sizeof(int), MPI_INFO_NULL, &xdens_1 );
  MPI_Alloc_mem( ndcells * sizeof(int), MPI_INFO_NULL, &xdens_2 );
#elif defined(MPI)
  xdens_1 = (int *) malloc( ndcells * sizeof(int) );
  xdens_2 = (int *) malloc( ndcells * sizeof(int) );
#else
  xdens_1 = (int *) malloc( ndcells * sizeof(int) );
  xdens_2 = (int *) malloc( ndcells * sizeof(int) );
#endif

  for (l=0; l<(ndcells); l++) {
    xdens_1[l] = 0;
    xdens_2[l] = 0;
   }

  /* sum over the density-cells */
 
  for (l=0; l<(ndcells); l++) {
    
    /* check if an atoms x-coordiante is in the l-th dcell */
    
    for (k=0; k<NCELLS; k++) {
      cell *p;
      p = CELLPTR(k);
      int i;
      
      /* if so, add 1 atom to the xdens of the l-th dcell */
      
      for (i=0; i < p->n; i++) {
	if( (ORT(p,i,X) > (double)l * deltax ) && (ORT(p,i,X) < (double)(l+1) * deltax ) )
	  xdens_2[l] += 1;
	
      }
    }
  } 

  /* add the resultes for the different CPUs */
  
#ifdef MPI
  MPI_Reduce( xdens_2, xdens_1, ndcells, MPI_INT, MPI_SUM, 0, cpugrid);
#else
  for(l=0; l<ndcells;l++)
  {
    xdens_1[l] = xdens_2[l];
  }
#endif
 
  if(myid==0){
  
    /* find the right side of the sample, which is the most outer dcell with density != 0 */
    
    for (l=ndcells-1;l>0;l--)
    {
      if( (xdens_1[l] == 0) && (xdens_1[l-1]  != 0) )
      {
	rightside = l-1;
	break;
      }
      
    }


#ifdef DEBUG
    for (l=0; l<ndcells; l++)
      printf("num(%i): %i, intervall: %f - %f \n", l, xdens_1[l], l* deltax, (l+1)* deltax); 
#endif

    /* find the actual surface (not clusters or vapour)
       by starting from the very right side and look for two adjacent dcells with density 0 */
    
    for (l=rightside; l>0; l--){
      if((xdens_1[l]==0)&&(xdens_1[l-1]==0))
	break;
    }
    
    leftside = l+1;
    
    /* check if there are only a few atoms inside the actual top 2 guessed 
       'surface dcells'; if so adjusts the surface slightly; do also for
       the rightside */
    
    if ((xdens_1[leftside]<500))
    {
      if ((xdens_1[leftside+1]<500))
	leftside = l+3;
      else
	leftside = l+2;
    }
    
    
    if ((xdens_1[rightside]<500))
    {
      if((xdens_1[rightside-1]<500))
	rightside-=2;
      else
	rightside-=1;
    }
       
  }

  /* free arrays*/
#ifdef MPI2
  MPI_Free_mem(xdens_1);
  MPI_Free_mem(xdens_2);
#elif defined(MPI)
  free(xdens_1);
  free(xdens_2);
#else
  free(xdens_1);
  free(xdens_2);
#endif
  
  laser_atom_vol=calc_laser_atom_vol(deltax, leftside, rightside, xdens_1);
   
  
  /* send the new value for laser_atom_vol to the other CPUs */
#ifdef MPI
  MPI_Bcast( &laser_atom_vol,     1, MPI_DOUBLE,  0, MPI_COMM_WORLD);
#endif

#ifdef PDECAY
  double samplesize = (rightside - leftside)*deltax;
    
  ramp_start = (1.0 - ramp_fraction) * samplesize + (double)(leftside+0.5)*deltax;
  ramp_end = (double)(rightside) * deltax;

#ifdef MPI
  MPI_Bcast( &ramp_start,   1, REAL,  0, MPI_COMM_WORLD);
  MPI_Bcast( &ramp_end,     1, REAL,  0, MPI_COMM_WORLD);
#endif

#endif


  /* return the new actual surface x-coordiante (parameter laser_offset) */
  return (double)(leftside+0.5)*deltax;
  
}



void init_laser()
{

 int k;

#ifdef LASERYZ
  laser_p_peak0=laser_mu*laser_sigma_e/laser_sigma_t/sqrt(2*M_PI);
  //printf("ppeak: %f mu: %f se: %f", laser_p_peak, laser_mu, laser_sigma_e);
#else  
   laser_p_peak=laser_mu*laser_sigma_e/laser_sigma_t/sqrt(2*M_PI); 
#endif
  laser_sigma_t_squared=laser_sigma_t*laser_sigma_t;

  if (laser_t_1>0){
#ifdef LASERYZ
    laser_p_peak1=laser_mu*laser_sigma_e1/laser_sigma_t1/sqrt(2*M_PI);
    //printf("ppeak: %f mu: %f se: %f", laser_p_peak, laser_mu, laser_sigma_e);
#else  
    laser_p_peak1=laser_mu*laser_sigma_e1/laser_sigma_t1/sqrt(2*M_PI); 
#endif
    laser_sigma_t1_squared=laser_sigma_t1*laser_sigma_t1;
  }

#ifdef LASERYZ
  laser_sigma_w0 = 1.0 / ( laser_sigma_w0 * laser_sigma_w0 );
#endif

  /* get the surface coordinate and write it into laser_offset */
  
  laser_offset = get_surface();
 
#ifdef MPI
  MPI_Bcast( &laser_offset,     1, REAL,  0, MPI_COMM_WORLD);
#endif

#ifdef LASERYZ
    /* if beam-coordinates are given in realtive values*/
    /* rescale them to absolute values */
    if( (laser_sigma_w_y > 0.0) && (laser_sigma_w_y < 1.0) ){
      laser_sigma_w_y *= box_y.y;
    }
    if( (laser_sigma_w_z > 0.0) && (laser_sigma_w_z < 1.0) ){
      laser_sigma_w_z *= box_z.z;
    }   
#endif
  
  if(myid==0){  
    printf("laser offset is set to: %f \n", laser_offset);
    printf("laser_atom_vol is set to: %f \n", laser_atom_vol);
    
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

#ifdef LASERYZ
      printf( "Laser beam diameter is %3.2f A. \n", 2.0 / ( sqrt(laser_sigma_w0 ) ) );
      printf( "Laser beam hits the target at y = %3.2f A z = %3.2f A.\n", laser_sigma_w_y, laser_sigma_w_z );
      printf( "Laser offset is %f.\n", laser_offset);
      printf( "Laser atom_vol is %f. \n", laser_atom_vol);
#endif
      printf( "Laser energy density is %1.10f\n", laser_sigma_e);
      printf( "Laser pulse duration (sigma) is %1.10f\n", laser_sigma_t);
      printf( "Time t_0 of laser pulse is %1.10f\n",laser_t_0);
      printf( "(%1.10f time steps after start of simulation)\n", 
	      laser_t_0/timestep);

      if (laser_t_1>0) {
      printf( "Laser energy density of the second pulse is %1.10f\n", laser_sigma_e1);
      printf( "Laser pulse duration (sigma) of the second pulse is %1.10f\n", laser_sigma_t1);
      printf( "Time t_1 of the second laser pulse is %1.10f\n",laser_t_1);
      printf( "(%1.10f time steps after start of simulation)\n", 
	      laser_t_1/timestep);
      }

#ifdef PDECAY
  if (0==myid) {
    printf( "Parameter pdecay_mode is %i \n", pdecay_mode );
    printf( "Parameter xipdecay is %f \n", xipdecay );
    printf( "Parameter ramp_start is %f A \n", ramp_start );
    printf( "Parameter ramp_end is %f A \n", ramp_end );
  }
#endif
      
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

  real exp_gauss_time_etc, gauss_time_squared, gauss_time_squared1;
  int k;
  real cosph, sinph, theta, costh, sinth;

  gauss_time_squared = timestep * steps - laser_t_0;
  gauss_time_squared *= gauss_time_squared;
  exp_gauss_time_etc = exp(-gauss_time_squared/laser_sigma_t_squared/2.0)
                       * laser_p_peak * timestep * laser_atom_vol;
  if (laser_t_1>0) {
    gauss_time_squared1 = timestep * steps - laser_t_1;
    gauss_time_squared1 *= gauss_time_squared1;
    exp_gauss_time_etc += exp(-gauss_time_squared1/laser_sigma_t1_squared/2.0)
                       * laser_p_peak1 * timestep * laser_atom_vol;
  }


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

#ifdef PDECAY
	if( ORT(p,i,X) > ramp_start )
	  {
	    switch ( pdecay_mode){
	    case 0:
	      {
		/* y= m * x + b */
		double m = 1 / ( (ramp_end - ramp_start ) ) ;
		double b = - ramp_start / ( ramp_end - ramp_start );
		
		double xia = ORT(p,i,X) *m + b ;
	
		IMPULS(p,i,X) *= 1.0 - ( ORT(p,i,X) *m + b)  ;
		
		if(steps==-1)
		  printf(" %f %f \n",  ORT(p,i,X),1.0 - ( ORT(p,i,X) * m + b) );
		break;
	      }
	    case 1:	      
	      {
		double a = 1.0 / ( ramp_start - ramp_end);
		a *= a;
	
		IMPULS(p,i,X) *= a * ( ORT(p,i,X) - ramp_end ) * ( ORT(p,i,X) - ramp_end );

		if(steps==-1)
		  printf("%f %f \n",  ORT(p,i,X), a * ( ORT(p,i,X) - ramp_end ) * ( ORT(p,i,X) - ramp_end ));
		break;
	      }
	      break;
	    case 2:
	      {
		double m = 1 / ( (ramp_end - ramp_start ) ) ;
		double b = - ramp_start / ( ramp_end - ramp_start );
		
		KRAFT(p,i,X) -=  ( IMPULS(p,i,X)/MASSE(p,i)) * ( ORT(p,i,X) *m + b  ) * xipdecay;
		
		
		if(steps==-1)
		  printf(" %f %f \n",  ORT(p,i,X),ORT(p,i,X) *m + b);
		break;
	      }
	      break;
	    case 3:
	      {
		double a = 1.0 / ( ramp_end - ramp_start);
		a *= a;
			
		 KRAFT(p,i,X) -=  ( IMPULS(p,i,X)/MASSE(p,i)) * xipdecay * a * ( ORT(p,i,X) - ramp_start ) * ( ORT(p,i,X) - ramp_start );

		if(steps==-1)
		  printf("%f %f \n",  ORT(p,i,X), a * ( ORT(p,i,X) - ramp_start ) * ( ORT(p,i,X) - ramp_start ) );
		break;
	      }
	      break;
	    default:
	      error("Illegal value for parameter pdecay_mode.\n");
	      break;
	    }
	  }
#endif 


#ifdef LASERYZ  /* spacial dependence of laser beam */    
    
     
      double x = ORT(p,i,X);
      double y = ORT(p,i,Y);
      double z = ORT(p,i,Z);
   
      
     
      de = exp( -laser_mu*depth ) * exp_gauss_time_etc * laser_intensity_profile(x,y,z);     
  
    

#else
      de = exp(-laser_mu*depth) * exp_gauss_time_etc;
       
#endif         

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

  real exp_gauss_time_etc, gauss_time_squared, gauss_time_squared1;
  int k;

  gauss_time_squared = timestep * steps - laser_t_0;
  gauss_time_squared *= gauss_time_squared;
  exp_gauss_time_etc = exp(-gauss_time_squared/laser_sigma_t_squared/2.0)
                       * laser_p_peak * timestep * laser_atom_vol;
  if (laser_t_1>0) {
    gauss_time_squared1 = timestep * steps - laser_t_1;
    gauss_time_squared1 *= gauss_time_squared1;
    exp_gauss_time_etc += exp(-gauss_time_squared1/laser_sigma_t1_squared/2.0)
                       * laser_p_peak1 * timestep * laser_atom_vol;
  }

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

#ifdef LASERYZ  /* spacial dependence of laser beam */    
    
    
      double x = ORT(p,i,X);
      double y = ORT(p,i,Y);
      double z = ORT(p,i,Z);
      de = exp( -laser_mu*depth ) * exp_gauss_time_etc * laser_intensity_profile(x,y,z);   

#else
      de = exp(-laser_mu*depth) * exp_gauss_time_etc;
    
#endif
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
  real exp_gauss_time_etc, gauss_time_squared, gauss_time_squared1, depth;
  
  gauss_time_squared = timestep * steps - laser_t_0;
  gauss_time_squared *= gauss_time_squared;
  exp_gauss_time_etc = exp(-gauss_time_squared/laser_sigma_t_squared/2.0)
    * laser_p_peak ;/*  * fd_h.x*fd_h.y*fd_h.z not needed */
  if (laser_t_1>0) {
    gauss_time_squared1 = timestep * steps - laser_t_1;
    gauss_time_squared1 *= gauss_time_squared1;
    exp_gauss_time_etc += exp(-gauss_time_squared1/laser_sigma_t1_squared/2.0)
      * laser_p_peak1;/*  * fd_h.x*fd_h.y*fd_h.z not needed */
  }

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

#ifdef PDECAY

  for (k=0; k<NCELLS; k++) {
    cell *p;
    p = CELLPTR(k);
    int i;
    for (i=0; i < p->n; i++) {
	if( ORT(p,i,X) > ramp_start )
	  {
	    switch ( pdecay_mode){
	    case 0:
	      {
		double m = 1 / ( (ramp_end - ramp_start ) ) ;
		double b = - ramp_start / ( ramp_end - ramp_start );
		
		double xia = ORT(p,i,X) *m + b ;
	
		IMPULS(p,i,X) *= 1.0- ( ORT(p,i,X) *m + b)  ;
		
		if(steps==-1)
		  printf(" %f %f \n",  ORT(p,i,X),1.0 - ( ORT(p,i,X) * m + b) );
		break;
	      }
	    case 1:	      
	      {
		double a = 1.0 / ( ramp_start - ramp_end);
		a *= a;
	
		IMPULS(p,i,X) *= a * ( ORT(p,i,X) - ramp_end ) * ( ORT(p,i,X) - ramp_end );

		if(steps==-1)
		  printf("%f %f \n",  ORT(p,i,X), a * ( ORT(p,i,X) - ramp_end ) * ( ORT(p,i,X) - ramp_end ));
		break;
	      }
	      break;
	    case 2:
	      {
		double m = 1 / ( (ramp_end - ramp_start ) ) ;
		double b = - ramp_start / ( ramp_end - ramp_start );
		
		KRAFT(p,i,X) -=  ( IMPULS(p,i,X)/MASSE(p,i)) * ( ORT(p,i,X) *m + b  ) * xipdecay;
		
		
		if(steps==-1)
		  printf(" %f %f \n",  ORT(p,i,X),ORT(p,i,X) *m + b);
		break;
	      }
	      break;
	    case 3:
	      {
		double a = 1.0 / ( ramp_end - ramp_start);
		a *= a;
			
		 KRAFT(p,i,X) -=  ( IMPULS(p,i,X)/MASSE(p,i)) * xipdecay * a * ( ORT(p,i,X) - ramp_start ) * ( ORT(p,i,X) - ramp_start );

		if(steps==-1)
		  printf("%f %f \n",  ORT(p,i,X), a * ( ORT(p,i,X) - ramp_start ) * ( ORT(p,i,X) - ramp_start ) );
		break;
	      }
	      break;
	    default:
	      error("Illegal value for parameter pdecay_mode.\n");
	      break;
	    }
	  }

    }
  }
#endif 

}
#endif


