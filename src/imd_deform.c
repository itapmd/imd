/******************************************************************************
*
* imd_shear_hom.c -- shear sample homogeneously
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

#ifdef HOM

/*****************************************************************************
*
* shear_sample()
*
*****************************************************************************/



void shear_sample(void)
{

  int i,r,s,t;
  cell *p;
  vektor2d d,u;
  real umax,umin;
  real tmp_umax,tmp_umin;
  real xmax,xmin;
  real tmp_xmax,tmp_xmin;
  real sclx;
  real theta;
  real radius;
  real amue;
  real kappa;
  int flag=0;
  
  
#ifdef MPI
  if (0==myid)
#endif
      printf("Shearing sample.\n");
  
      /* Seek for Min/Max u and x */
  for ( r = cellmin.x; r < cellmax.x; ++r )
      for ( s = cellmin.y; s < cellmax.y; ++s )
#ifndef TWOD
          for ( t = cellmin.z; t < cellmax.z; ++t )
#endif
          {
        
#ifdef TWOD
            p = PTR_2D_V(cell_array, r, s, cell_dim);
#else
            p = PTR_3D_V(cell_array, r, s, t, cell_dim);
#endif
            for (i = 0;i < p->n; ++i) {

/* Catch first value */
              if (0==flag) {
                flag = -1;
                xmax = p->ort X(i);
                xmin = p->ort X(i);
              } else {
                xmax = MAX(xmax,p->ort X(i));
                xmin = MIN(xmin,p->ort X(i));
              };
            }
          }

#ifdef MPI
  MPI_Allreduce(&xmax, &tmp_xmax, 1, MPI_REAL, MPI_MAX, cpugrid);
  MPI_Allreduce(&xmin, &tmp_xmin, 1, MPI_REAL, MPI_MIN, cpugrid);

  xmax = tmp_xmax;
  xmin = tmp_xmin;
#endif

#ifdef MPI
  if (0==myid)
#endif
      printf("Min X: %f\nMax X: %f\n",
             xmin,xmax);

  
          /* Apply shear */
  for ( r = cellmin.x; r < cellmax.x; ++r )
      for ( s = cellmin.y; s < cellmax.y; ++s )
#ifndef TWOD
          for ( t = cellmin.z; t < cellmax.z; ++t )
#endif
          {
        
#ifdef TWOD
            p = PTR_2D_V(cell_array, r, s, cell_dim);
#else
            p = PTR_3D_V(cell_array, r, s, t, cell_dim);
#endif
            for (i = 0;i < p->n; ++i)
              p->ort Y(i) += shear_max * p->ort X(i) / xmax;
                            
          }
#ifdef MPI
  if (0==myid)
#endif
  return;
} /*shear_sample*/

/*****************************************************************************
*
* expand_sample()
*
*****************************************************************************/

void expand_sample(void)

{

  int i,r,s,t;
  cell *p;
  vektor2d d,u;
  real umax,umin;
  real tmp_umax,tmp_umin;
  real xmax,xmin;
  real tmp_xmax,tmp_xmin;
  real sclx;
  real theta;
  real radius;
  real amue;
  real kappa;
  int flag=0;
  
  
#ifdef MPI
  if (0==myid)
#endif
      printf("Expanding sample.\n");
  
          /* Apply field */
  for ( r = cellmin.x; r < cellmax.x; ++r )
      for ( s = cellmin.y; s < cellmax.y; ++s )
#ifndef TWOD
          for ( t = cellmin.z; t < cellmax.z; ++t )
#endif
          {
        
#ifdef TWOD
            p = PTR_2D_V(cell_array, r, s, cell_dim);
#else
            p = PTR_3D_V(cell_array, r, s, t, cell_dim);
#endif
            for (i = 0;i < p->n; ++i) {

              p->ort X(i) *= expansion;
              p->ort Y(i) *= expansion;
	      box_x.x *= expansion;
	      box_y.y *= expansion;
#ifndef TWOD
              p->ort Z(i) *= expansion;
	      box_z.z *= expansion;
#endif
                            
            }
          }
#ifdef MPI
  if (0==myid)
#endif
      printf("Expansion done.\n");
  return;
} /* expand sample */

#endif /* HOM */


#ifdef PULL

/*****************************************************************************
*
* deform_atoms()
*
*****************************************************************************/

void deform_atoms(void) {

  cell *p;
  int i;
  int r,s,t;
  real tot_kin_energy_prev;

  tot_kin_energy_prev = tot_kin_energy;

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
	  if ((dnoshsteps > annealsteps) && 
	      ((tot_kin_energy_prev < ekin_threshold) || 
	       (dnoshsteps > maxdnoshsteps))) {
	    /* which direction of pulling? */
	    if (p->ort X(i) <= strip_width) {
	      p->ort X(i) -= strip_shift.x;
	      p->ort Y(i) -= strip_shift.y;
#ifndef TWOD
	      p->ort Z(i) -= strip_shift.z;
#endif
	    }
	    if (p->ort X(i) >= box_x.x - strip_width) {
	      p->ort X(i) += strip_shift.x;
	      p->ort Y(i) += strip_shift.y;
#ifndef TWOD
	      p->ort Z(i) += strip_shift.z;
#endif        
	    }
	    dnoshsteps = 0;
	  } else {
	    dnoshsteps++;
	  }
	  p->impuls X(i) = 0.0;
	  p->impuls Y(i) = 0.0;
#ifndef TWOD
	  p->impuls Z(i) = 0.0;
#endif
	} /* i - loop */
      }

} /* deform_atoms */

#endif /* PULL */


