/***************************************************************************
 *  heat conductivity 
 *
 * $RCSfile$
 * $Revision$
 * $Date$
 ***************************************************************************/

#include "imd.h"

/******************************************************************************
*
* write_temp 
*
******************************************************************************/

void write_temp(int steps)

{
  FILE  *outtemp;
  str255  fnametemp;
  cell *p;
  real scale;
  integer num, nhalf;
  int i, r, s, t;
  static real    *temp_hist, *temp_hist_1 = NULL, *temp_hist_2 = NULL;
  static integer *num_hist,   *num_hist_1 = NULL,  *num_hist_2 = NULL;
  
  /* the temp bins are orthogonal boxes in space */
  nhalf = tran_nlayers / 2;
  scale = tran_nlayers / box_x.x;

  /* allocate histogram arrays */
  if (NULL==temp_hist_1) {
    temp_hist_1 = (real *) malloc( (nhalf+1) * sizeof(real) );
    if (NULL==temp_hist_1) 
      error("Cannot allocate temperature  array.");
  }  
  if (NULL==num_hist_1) {
    num_hist_1 = (integer *) malloc( (nhalf+1) * sizeof(integer) );
    if (NULL==num_hist_1) 
      error("Cannot allocate temperature  array.");
  }
#ifdef MPI
  if (NULL==temp_hist_2) {
    temp_hist_2 = (real *) malloc( (nhalf+1) * sizeof(real) );
    if ((NULL==temp_hist_2) && (myid==0)) 
      error("Cannot allocate temperature  array.");
  }  
  if (NULL==num_hist_2) {
    num_hist_2 = (integer *) malloc( (nhalf+1) * sizeof(integer) );
    if ((NULL==num_hist_2) && (myid==0)) 
      error("Cannot allocate temperature  array.");
  }
#endif

  for (i = 0; i <= nhalf; i++) {
    temp_hist_1[i] = 0.0;
     num_hist_1[i] = 0;
#ifdef MPI
    temp_hist_2[i] = 0.0;
     num_hist_2[i] = 0;
#endif
  }

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

              /* which layer? */
              num = scale * p->ort X(i);
              if (num > nhalf) num = tran_nlayers - num;
              if (num < 0)     num = 0;

              temp_hist_1[num] += SPRODN(p->impuls,i,p->impuls,i)/p->masse[i];
              num_hist_1[num]++;
	   }
        }

#ifdef MPI
  /* add up results form different CPUs */
  MPI_Allreduce( temp_hist_1, temp_hist_2, 
                 nhalf+1, MPI_REAL, MPI_SUM, cpugrid);
  temp_hist = temp_hist_2;
  MPI_Allreduce( num_hist_1, num_hist_2, 
                 nhalf+1, INTEGER, MPI_SUM, cpugrid);
  num_hist  = num_hist_2;
#else
  temp_hist = temp_hist_1;
  num_hist  = num_hist_1;
#endif

  /* write temperature distribution */
  if (myid==0) {
    sprintf(fnametemp,"%s.tempdist",outfilename);
    outtemp = fopen(fnametemp,"a");
    if (NULL == outtemp) error("Cannot open temperatur file.");

    fprintf(outtemp,"%10.4e", steps * timestep);
    for (i = 0; i <= nhalf; i++) {
      if (num_hist[i] > 0) temp_hist[i] /= (2*num_hist[i]);
#ifndef TWOD   
      temp_hist[i] *= (2.0/DIM);
#endif
      fprintf(outtemp," %10.4e", temp_hist[i] );
    }
    fprintf(outtemp,"\n");
    fclose(outtemp);
  }

}







































































