
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_transport.c -- heat transport (and pressure histograms)
*
******************************************************************************/

/***************************************************************************
 * $Revision$
 * $Date$
 ***************************************************************************/

#include "imd.h"


/******************************************************************************
*
*  rnmend_heat_exchange -- exchange velocities of two particles
*
*  current limitation: one particle type, serial only
*
******************************************************************************/

void rnemd_heat_exchange()
{
  int  k;
  real swap;

  cell *mincell, *maxcell;
  int  minatom, maxatom;
  real minEkin=tot_kin_energy, maxEkin=0.0;

  int  nhalf = tran_nlayers / 2;
  real scale = tran_nlayers / box_x.x;

  /* find hot and cold particles to be swapped */
  for (k=0; k<ncells; ++k) {

    int  i,num;
    cell *p;
    real tmp;
    p = cell_array + CELLS(k);

    for (i=0; i<p->n; ++i) {

      tmp = SPRODN(p->impuls,i,p->impuls,i) / (2 * MASSE(p,i));

      /* which layer? */
      num = scale * p->ort X(i);
      if (num < 0)             num = 0;
      if (num >= tran_nlayers) num = tran_nlayers-1;

      /* minimum in hot layer */
      if ((minEkin > tmp) && (num==0)) {
        minEkin = tmp;
        mincell = p;
        minatom = i;
      } 
      /* maximum in cold layer */
      if ((maxEkin < tmp) && (num==nhalf)) {
        maxEkin = tmp;
        maxcell = p;
        maxatom = i;
      }
    }
  }

  /* swap the velocities */
  swap = maxcell->impuls X(maxatom);
  maxcell->impuls X(maxatom) = mincell->impuls X(minatom);
  mincell->impuls X(minatom) = swap;
  swap = maxcell->impuls Y(maxatom);
  maxcell->impuls Y(maxatom) = mincell->impuls Y(minatom);
  mincell->impuls Y(minatom) = swap;
#ifndef TWOD
  swap = maxcell->impuls Z(maxatom);
  maxcell->impuls Z(maxatom) = mincell->impuls Z(minatom);
  mincell->impuls Z(minatom) = swap;
#endif

  /* accumulate heat transfer */
  heat_transfer += maxEkin - minEkin;

}


/******************************************************************************
*
* write_temp_dist
*
******************************************************************************/

void write_temp_dist(int steps)
{
  FILE  *outtemp;
  str255  fnametemp;
  real scale, vol;
  int  num, nlayer, k, i, nhalf;
  static real    *temp_hist, *temp_hist_1 = NULL, *temp_hist_2 = NULL;
  static integer *num_hist,   *num_hist_1 = NULL,  *num_hist_2 = NULL;
  int numb = 0, numb_tmp;
  real a, SxiTi=0.0, Sxi=0.0, STi=0.0, Sxisq=0.0, temp, xx, tmp;
  
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
  for (k=0; k<ncells; ++k) {

    int  i;
    cell *p;
    p = cell_array + CELLS(k);

    for (i=0; i<p->n; ++i) {

      /* which layer? */
      xx = p->ort X(i);
      num = scale * xx;
      if (num < 0)             num = 0;
      if (num >= tran_nlayers) num = tran_nlayers-1;
      if (num > nhalf) {
        num = tran_nlayers - num;
        xx  = box_x.x - xx + box_x.x / tran_nlayers; 
      }
      temp = SPRODN(p->impuls,i,p->impuls,i) / (2*MASSE(p,i));
      temp_hist_1[num] += temp;
      num_hist_1[num]++;

      /* fit temperature gradient */
      if ((num!=0) && (num!=nhalf)) {
        Sxi += xx;
        STi += temp;
        SxiTi += temp * xx;
        Sxisq += xx * xx;
        numb++;
      }
    }
  }

#ifdef MPI
  /* add up results form different CPUs */
#ifdef NVX
  MPI_Allreduce( &heat_cond, &tmp, 1, REAL, MPI_SUM, cpugrid);
  heat_cond = tmp;
#endif
  MPI_Reduce( temp_hist_1, temp_hist_2, 
              nhalf+1, REAL, MPI_SUM, 0, cpugrid);
  temp_hist = temp_hist_2;
  MPI_Reduce( num_hist_1, num_hist_2, 
              nhalf+1, INTEGER,  MPI_SUM, 0, cpugrid);
  num_hist  = num_hist_2;
  MPI_Reduce( &Sxi, &tmp, 1, REAL, MPI_SUM, 0, cpugrid );
  Sxi = tmp;
  MPI_Reduce( &STi, &tmp, 1, REAL, MPI_SUM, 0, cpugrid );
  STi = tmp;
  MPI_Reduce( &SxiTi, &tmp, 1, REAL, MPI_SUM, 0, cpugrid );
  SxiTi = tmp;
  MPI_Reduce( &Sxisq, &tmp, 1, REAL, MPI_SUM, 0, cpugrid );
  Sxisq = tmp;
  MPI_Reduce( &numb, &numb_tmp, 1, MPI_INT, MPI_SUM, 0, cpugrid );
  numb = numb_tmp; 
#else
  temp_hist = temp_hist_1;
  num_hist  = num_hist_1;
#endif


  /* write temperature distribution */
  if (myid==0) {

    Sxi /= numb;
    STi /= numb;
    SxiTi /= numb;
    Sxisq /= numb;
    a = - (SxiTi - Sxi * STi) / (Sxisq - Sxi * Sxi);

    sprintf(fnametemp,"%s.tempdist",outfilename);
    outtemp = fopen(fnametemp,"a");
    if (NULL == outtemp) error("Cannot open temperatur file.");

    /* write current time */
    fprintf(outtemp,"%10.4e", steps * timestep);

#ifdef RNEMD
    /* write heat current density determined from heat transfer */
    /* heat flows away in two directions -> twice the cross section */
    heat_transfer /= 2 * box_y.y * tran_interval * timestep;
#ifndef TWOD
    heat_transfer /= box_z.z;
#endif
    fprintf(outtemp," %10.4e", heat_transfer);
    heat_transfer = 0.0;
#endif

#ifdef NVX
#ifdef TWOD
    vol = box_x.x * box_y.y           * (tran_nlayers-2) / tran_nlayers;
#else
    vol = box_x.x * box_y.y * box_z.z * (tran_nlayers-2) / tran_nlayers;
#endif
    heat_cond /= (vol * tran_interval);
    /* write heat current density determined by Gillan-Evans-Algorithm */
    fprintf(outtemp," %10.4e", heat_cond);
    heat_cond = 0.0;
#endif

    /* write measured temperature gradient */
    fprintf(outtemp," %10.4e", a);

#ifdef NVX
    /* write nominal temperature gradient */
    fprintf(outtemp," %10.4e", 2 * (tran_Tleft - tran_Tright) / box_x.x );
#endif

    /* write temperature histogram */
    for (i = 0; i <= nhalf; i++) {
      if (num_hist[i] > 0) temp_hist[i] /= num_hist[i];
#ifndef TWOD
      temp_hist[i] *= (2.0/DIM);
#endif
      fprintf(outtemp," %10.4e", temp_hist[i] );
    }

    fprintf(outtemp,"\n");
    fclose(outtemp);
  }
}
