/***************************************************************************
 * uniaxial histograms
 *
 * $RCSfile$
 * $Revision$
 * $Date$
 ***************************************************************************/

#include "imd.h"

#ifdef TRANSPORT
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
              if (num < 0)             num = 0;
              if (num >= tran_nlayers) num = tran_nlayers-1;
              if (num > nhalf) num = tran_nlayers - num;

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
    /* the variable heat_cond was formerly written to the .eng file */
    fprintf(outtemp," %10.4e", heat_cond);
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

#endif

#ifdef STRESS_TENS
/******************************************************************************
*
* write_press 
*
******************************************************************************/

void write_press(int steps)

{
  FILE  *outpress;
  str255  fnamepress;
  cell *p;
  real scalex, scaley, layvol, laydens;
  integer numx, numy, press_nlayers_tot;
  int fzhlr, i, r, s, t;
#ifndef TWOD
  int scalez, numz;
#endif
  static real  *press_histxx, *press_histxx_1 = NULL, *press_histxx_2 = NULL;
  static real  *press_histyy, *press_histyy_1 = NULL, *press_histyy_2 = NULL;
  static real  *press_histzz, *press_histzz_1 = NULL, *press_histzz_2 = NULL;
  static real  *kin_hist, *kin_hist_1 = NULL, *kin_hist_2 = NULL;
  static real  *pot_hist, *pot_hist_1 = NULL, *pot_hist_2 = NULL;
  static integer *num_hist,   *num_hist_1 = NULL,  *num_hist_2 = NULL;
  
  /* the press bins are orthogonal boxes in space */
  scalex = press_nlayers.x / box_x.x;
  scaley = press_nlayers.y / box_y.y;
  press_nlayers_tot = press_nlayers.x*press_nlayers.y;
#ifndef TWOD
  scalez = press_nlayers.z / box_z.z;
  press_nlayers_tot = press_nlayers.x*press_nlayers.y*press_nlayers.z;
#endif
  layvol = volume / press_nlayers_tot;

  /* allocate histogram arrays */
  if (NULL==press_histxx_1) {
    press_histxx_1 = (real *) malloc( press_nlayers_tot * sizeof(real) );
    if (NULL==press_histxx_1) 
      error("Cannot allocate xx pressure tensor array.");
  }  
  if (NULL==press_histyy_1) {
    press_histyy_1 = (real *) malloc( press_nlayers_tot * sizeof(real) );
    if (NULL==press_histyy_1) 
      error("Cannot allocate yy pressure tensor array.");
  }
#ifndef TWOD  
  if (NULL==press_histzz_1) {
    press_histzz_1 = (real *) malloc( press_nlayers_tot * sizeof(real) );
    if (NULL==press_histzz_1) 
      error("Cannot allocate zz pressure tensor array.");
  }
#endif  
  if (NULL==kin_hist_1) {
    kin_hist_1 = (real *) malloc( press_nlayers_tot * sizeof(real) );
    if (NULL==kin_hist_1) 
      error("Cannot allocate kinetic energy array.");
  }  
  if (NULL==pot_hist_1) {
    pot_hist_1 = (real *) malloc( press_nlayers_tot * sizeof(real) );
    if (NULL==pot_hist_1) 
      error("Cannot allocate potential energy array.");
  }  
  if (NULL==num_hist_1) {
    num_hist_1 = (integer *) malloc( press_nlayers_tot * sizeof(integer) );
    if (NULL==num_hist_1) 
      error("Cannot allocate number count array.");
  }
#ifdef MPI
  if (NULL==press_histxx_2) {
    press_histxx_2 = (real *) malloc( press_nlayers_tot * sizeof(real) );
    if ((NULL==press_histxx_2) && (myid==0)) 
      error("Cannot allocate xx pressure tensor  array.");
  }  
  if (NULL==press_histyy_2) {
    press_histyy_2 = (real *) malloc( press_nlayers_tot * sizeof(real) );
    if ((NULL==press_histyy_2) && (myid==0)) 
      error("Cannot allocate yy pressure tensor  array.");
  }
#ifndef TWOD
  if (NULL==press_histzz_2) {
    press_histzz_2 = (real *) malloc( press_nlayers_tot * sizeof(real) );
    if ((NULL==press_histzz_2) && (myid==0)) 
      error("Cannot allocate zz pressure tensor  array.");
  }
#endif
  if (NULL==kin_hist_2) {
    kin_hist_2 = (real *) malloc( press_nlayers_tot * sizeof(real) );
    if (NULL==kin_hist_2) 
      error("Cannot allocate kinetic energy array.");
  }  
  if (NULL==pot_hist_2) {
    pot_hist_2 = (real *) malloc( press_nlayers_tot * sizeof(real) );
    if (NULL==pot_hist_2) 
      error("Cannot allocate potential energy array.");
  }  
  if (NULL==num_hist_2) {
    num_hist_2 = (integer *) malloc( press_nlayers_tot * sizeof(integer) );
    if ((NULL==num_hist_2) && (myid==0)) 
      error("Cannot allocate number count array.");
  }
#endif

  for (i = 0; i < press_nlayers_tot; i++) {
    press_histxx_1[i] = 0.0;
    press_histyy_1[i] = 0.0;
#ifndef TWOD
    press_histzz_1[i] = 0.0;
#endif
    kin_hist_1[i] = 0;
    pot_hist_1[i] = 0;
    num_hist_1[i] = 0;
#ifdef MPI
    press_histxx_2[i] = 0.0;
    press_histyy_2[i] = 0.0;
#ifndef TWOD
    press_histzz_2[i] = 0.0;
#endif
    kin_hist_2[i] = 0;
    pot_hist_2[i] = 0;
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
	      numx = scalex * p->ort X(i);
              if (numx < 0)             numx = 0;
              if (numx >= press_nlayers.x)      numx = press_nlayers.x-1;
	      numy = scaley * p->ort Y(i);
              if (numy < 0)             numy = 0;
              if (numy >= press_nlayers.y)      numy = press_nlayers.y-1;
#ifndef TWOD
	      numz = scalez * p->ort Z(i);
              if (numz < 0)             numz = 0;
              if (numz >= press_nlayers.z)      numz = press_nlayers.z-1;
#endif

#ifdef TWOD
              press_histxx_1[numx*press_nlayers.y+numy] += p->presstens X(i);
              press_histyy_1[numx*press_nlayers.y+numy] += p->presstens Y(i);
	      kin_hist_1[numx*press_nlayers.y+numy] += SPRODN(p->impuls,i,p->impuls,i) / (2*p->masse[i]);
	      pot_hist_1[numx*press_nlayers.y+numy] += p->pot_eng[i];
              num_hist_1[numx*press_nlayers.y+numy]++;
#else
              press_histxx_1[numx*press_nlayers.y*press_nlayers.z+numy*press_nlayers.z+numz] += p->presstens X(i);
              press_histyy_1[numx*press_nlayers.y*press_nlayers.z+numy*press_nlayers.z+numz] += p->presstens Y(i);
              press_histzz_1[numx*press_nlayers.y*press_nlayers.z+numy*press_nlayers.z+numz] += p->presstens Z(i);
	      kin_hist_1[numx*press_nlayers.y*press_nlayers.z+numy*press_nlayers.z+numz] += SPRODN(p->impuls,i,p->impuls,i) / (2*p->masse[i]);
	      pot_hist_1[numx*press_nlayers.y*press_nlayers.z+numy*press_nlayers.z+numz] += p->pot_eng[i];
              num_hist_1[numx*press_nlayers.y*press_nlayers.z+numy*press_nlayers.z+numz]++;
#endif
	   }
        }

#ifdef MPI
  /* add up results form different CPUs */
  MPI_Allreduce( press_histxx_1, press_histxx_2, 
                 DIM * press_nlayers., MPI_REAL, MPI_SUM, cpugrid);
  press_histxx = press_histxx_2;
  MPI_Allreduce( press_histyy_1, press_histyy_2, 
                 DIM * press_nlayers., MPI_REAL, MPI_SUM, cpugrid);
  press_histxx = press_histyy_2;
  MPI_Allreduce( press_histzz_1, press_histzz_2, 
                 DIM * press_nlayers.+, MPI_REAL, MPI_SUM, cpugrid);
  press_histzz = press_hist_2;
  MPI_Allreduce( kin_hist_1, kin_hist_2, 
                 DIM * press_nlayers.+, MPI_REAL, MPI_SUM, cpugrid);
  press_histzz = kin_hist_2;
  MPI_Allreduce( pot_hist_1, pot_hist_2, 
                 DIM * press_nlayers.+, MPI_REAL, MPI_SUM, cpugrid);
  press_histzz = pot_hist_2;
  MPI_Allreduce( num_hist_1, num_hist_2, 
                 press_nlayers., INTEGER, MPI_SUM, cpugrid);
  num_hist  = num_hist_2;
#else
  press_histxx = press_histxx_1;
  press_histyy = press_histyy_1;
  press_histzz = press_histzz_1;
  kin_hist = kin_hist_1;
  pot_hist = pot_hist_1;
  num_hist  = num_hist_1;
#endif

  /* write pressure tensor distribution */
  if (myid==0) {

    fzhlr = steps / press_interval;

    sprintf(fnamepress,"%s.%u.pressdist",outfilename,fzhlr);
    outpress = fopen(fnamepress,"w");
    if (NULL == outpress) error("Cannot open pressure tensor file.");

    /*    fprintf(outpress,"%10.4e\n", steps * timestep); */
    for (i = 0; i < press_nlayers_tot; i++) {
      if (num_hist[i] > 0) {
	laydens = num_hist[i] / layvol;
	kin_hist[i] /= num_hist[i];
	pot_hist[i] /= num_hist[i];
	press_histxx[i] /= num_hist[i];
	press_histyy[i] /= num_hist[i];
#ifndef TWOD
	press_histzz[i] /= num_hist[i];
#endif
      }
#ifdef DEBUG /* print histbox coordinate */
#ifndef TWOD
      fprintf(outpress,"%d %d %d %10.4e %10.4e %10.4e %10e %10.4e %10.4e\n", 
	      (i-i%(press_nlayers.y*press_nlayers.z))/(press_nlayers.y*press_nlayers.z)%press_nlayers.x,(i-i%press_nlayers.z)/press_nlayers.z%press_nlayers.y, i%press_nlayers.z,press_histxx[i], press_histyy[i], press_histzz[i], laydens, 
	      kin_hist[i], pot_hist[i] );
#else
      fprintf(outpress,"%d %d %10.4e %10.4e %10.4e %10.4e %10.4e\n", 
	      (i-i%press_nlayers.y)/5, i%press_nlayers.y, press_histxx[i], press_histyy[i], laydens,
	      kin_hist[i], pot_hist[i] );
#endif
#endif /* DEBUG */

#ifndef TWOD
      fprintf(outpress,"%10.4e %10.4e %10.4e %10e %10.4e %10.4e\n", 
	      press_histxx[i], press_histyy[i], press_histzz[i], laydens, 
	      kin_hist[i], pot_hist[i] );
#else
      fprintf(outpress,"%10.4e %10.4e %10.4e %10.4e %10.4e\n", 
	      press_histxx[i], press_histyy[i], laydens,
	      kin_hist[i], pot_hist[i] );
#endif    
    }
    fprintf(outpress,"\n");
    fclose(outpress);
  }

}
#endif






