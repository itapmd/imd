/******************************************************************************
*
* imd_pacx.c -- Replacement of MPI routines for pacx
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/******************************************************************************
*
* compute cpu coordinates from rank
*
******************************************************************************/

ivektor my_cart_coords(int myid)

#ifndef TWOD
{
  int xfactor;
  int yfactor;

  xfactor=cpu_dim.y*cpu_dim.z;
  yfactor=cpu_dim.z;
  my_coord.x=myid/xfactor;
  my_coord.y=(myid-my_coord.x*xfactor)/yfactor;
  my_coord.z=myid-my_coord.x*xfactor-my_coord.y*yfactor;
  return(my_coord);
}
#else
{
  int xfactor;
  
  xfactor=cpu_dim.y;
  
  my_coord.x=myid/xfactor;
  my_coord.y=myid-my_coord.x*xfactor;
  return(my_coord);
}
#endif

/******************************************************************************
*
* compute rank from cpu coordinates
*
******************************************************************************/

void my_cart_rank(ivektor my_coord)

#ifndef TWOD
{
  int xfactor;
  int yfactor;
  int rank;
  int *ort;

  xfactor=cpu_dim.y*cpu_dim.z;
  yfactor=cpu_dim.z;

  /* Berechne den Rang int rank */

  rank=my_coord.z+yfactor*my_coord.y+xfactor*my_coord.x;

  /* Berechne den Pointer int *ort */

  ort=PTR_3D_VV(cpu_ranks, my_coord, cpu_dim);

  /* Schreibe an die Adresse int *ort den Wert int rank hin */

  *ort=rank;
}
#else
{
  int xfactor;
  int rank;
  int *ort;

  xfactor=cpu_dim.y;

  rank=my_coord.y+xfactor*my_coord.x;
  ort=PTR_2D_VV(cpu_ranks, my_coord, cpu_dim);
  *ort=rank;
}
#endif

