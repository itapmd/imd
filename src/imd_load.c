/******************************************************************************
*
* imd_load.c -- load sample for crack
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

/*****************************************************************************
*
* load_sample()
*
*****************************************************************************/

#include "imd.h"

void load_sample(void)

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
  
  
  kappa = ( kel + 2 * mue ) / kel;
  amue  = 8 * mue;
  
#ifdef MPI
  if (0==myid)
#endif
      printf("Loading sample.\n");
  
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

                  /* Lock atoms in strip */
              
              if ((p->ort X(i) < strip_width) || (p->ort X(i) > (box_x.x - strip_width))
                  && p->nummer[i] > 0)
                  p->nummer[i] = - p->nummer[i];

              
                    
              d.x = p->ort X(i) - tip.x;
              d.y = p->ort Y(i) - tip.y;
              
              radius = sqrt(SPROD2D(d,d));
              theta = atan2( d.x, d.y );
              
              u.x = kcrit * sqrt(2*radius) / amue *
                  (((2*kappa - 1) * sin(0.5*theta)) - sin(1.5*theta));
              u.y = kcrit * sqrt(2*radius) / amue *
                  (((2*kappa + 1) * cos(0.5*theta)) - cos(1.5*theta));

/* Catch first value */
              if (0==flag) {
                flag = -1;
                umax = u.x;
                umin = u.x;
                xmax = p->ort X(i);
                xmin = p->ort X(i);
              } else {
                umax = MAX(umax,u.x);
                umin = MIN(umin,u.x);
                xmax = MAX(xmax,p->ort X(i));
                xmin = MIN(xmin,p->ort X(i));
              };
            }
          }

#ifdef MPI
  MPI_Allreduce(&umax, &tmp_umax, 1, MPI_REAL, MPI_MAX, cpugrid);
  MPI_Allreduce(&umin, &tmp_umin, 1, MPI_REAL, MPI_MIN, cpugrid);
  MPI_Allreduce(&xmax, &tmp_xmax, 1, MPI_REAL, MPI_MAX, cpugrid);
  MPI_Allreduce(&xmin, &tmp_xmin, 1, MPI_REAL, MPI_MIN, cpugrid);

  umax = tmp_umax;
  umin = tmp_umin;
  xmax = tmp_xmax;
  xmin = tmp_xmin;
#endif

#ifdef MPI
  if (0==myid)
#endif
      printf("Min X: %f\nMax X: %f\nMin u: %f\nMax u: %f\n",
             xmin,xmax,umin,umax);

  
          /* Apply load */
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

              d.x = p->ort X(i) - tip.x;
              d.y = p->ort Y(i) - tip.y;

              radius = sqrt(SPROD2D(d,d));
              theta = atan2( d.x, d.y );
              
              u.x = kcrit * sqrt(2*radius) / amue *
                  (((2*kappa - 1) * sin(0.5*theta)) - sin(1.5*theta));
              u.y = kcrit * sqrt(2*radius) / amue *
                  (((2*kappa + 1) * cos(0.5*theta)) - cos(1.5*theta));

              p->ort X(i) += u.x;
              p->ort Y(i) += u.y;
              
              if (p->ort X(i) > tip.x )
                  u.x = (p->ort X(i) - tip.x) / (xmax - tip.x) * (umax - u.x);
              else
                  u.x = (p->ort X(i) - tip.x) / (xmin - tip.x) * (umin - u.x);
              
              p->ort X(i) += u.x; 
              
            }
          }
#ifdef MPI
  if (0==myid)
#endif
      printf("Loading done.\n");
  return;
}






