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

/*****************************************************************************
*
* shear_sample()
*
*****************************************************************************/

#include "imd.h"

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
            for (i = 0;i < p->n; ++i) {

              p->ort Y(i) += shear_max * p->ort X(i) / xmax;
                            
            }
          }
  #ifdef MPI
  if (0==myid)
#endif
      printf("Shear done.\n");
  return;
}






