
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
*  imd_histogram.c -- energy and pressure histograms
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/******************************************************************************
*
*  make_histograms
*
******************************************************************************/

void make_histograms(hist_t *hist)
{
  cell *p;
  real scalex, scaley, scalez, Ekin, tmp;
  int  num, numx, numy, numz, size;
  int  i, j, k;
  static float *press_histxx_1 = NULL, *press_histxx_2 = NULL;
  static float *press_histyy_1 = NULL, *press_histyy_2 = NULL;
  static float *press_histxy_1 = NULL, *press_histxy_2 = NULL;
#ifndef TWOD
  static float *press_histzz_1 = NULL, *press_histzz_2 = NULL;
  static float *press_histzx_1 = NULL, *press_histzx_2 = NULL;
  static float *press_histyz_1 = NULL, *press_histyz_2 = NULL;
#endif
#ifdef SHOCK
  static float *kin_histxx_1   = NULL, *kin_histxx_2   = NULL;
  static float *kin_histxxu_1  = NULL, *kin_histxxu_2  = NULL;
  static float *kin_histyy_1   = NULL, *kin_histyy_2   = NULL;
#ifndef TWOD
  static float *kin_histzz_1   = NULL, *kin_histzz_2   = NULL;
#endif
#endif
  static float *kin_hist_1     = NULL, *kin_hist_2     = NULL;
  static float *pot_hist_1     = NULL, *pot_hist_2     = NULL;
  static integer *num_hist_1   = NULL, *num_hist_2     = NULL;

  /* the bins are orthogonal boxes in space */
  scalex = hist->dim.x / box_x.x;
  scaley = hist->dim.y / box_y.y;
  size   = hist->dim.x * hist->dim.y;
#ifndef TWOD
  scalez = hist->dim.z / box_z.z;
  size  *= hist->dim.z;
#endif
  hist->binvol = volume / size;
  hist->size   = size;

  /* allocate histogram arrays */

#ifdef STRESS_TENS
  if (NULL==press_histxx_1) {
    press_histxx_1 = (float *) malloc( size * sizeof(float) );
    if (NULL==press_histxx_1) 
      error("Cannot allocate xx pressure tensor array.");
  }  
  if (NULL==press_histyy_1) {
    press_histyy_1 = (float *) malloc( size * sizeof(float) );
    if (NULL==press_histyy_1) 
      error("Cannot allocate yy pressure tensor array.");
  }
#ifndef SHOCK
  if (NULL==press_histxy_1) {
    press_histxy_1 = (float *) malloc( size * sizeof(float) );
    if (NULL==press_histxy_1) 
      error("Cannot allocate xy pressure tensor array.");
  }
#endif
#ifndef TWOD  
  if (NULL==press_histzz_1) {
    press_histzz_1 = (float *) malloc( size * sizeof(float) );
    if (NULL==press_histzz_1) 
      error("Cannot allocate zz pressure tensor array.");
  }
#ifndef SHOCK
  if (NULL==press_histzx_1) {
    press_histzx_1 = (float *) malloc( size * sizeof(float) );
    if (NULL==press_histzx_1) 
      error("Cannot allocate zx pressure tensor array.");
  }
  if (NULL==press_histyz_1) {
    press_histyz_1 = (float *) malloc( size * sizeof(float) );
    if (NULL==press_histyz_1) 
      error("Cannot allocate yz pressure tensor array.");
  }
#endif /* SHOCK */
#endif /* not TWOD */
#endif /* STRESS_TENS */

#ifdef SHOCK
  if (NULL==kin_histxx_1) {
    kin_histxx_1 = (float *) malloc( size * sizeof(float) );
    if (NULL==kin_histxx_1) 
      error("Cannot allocate xx kinetic energy tensor array.");
  }  
  if (NULL==kin_histyy_1) {
    kin_histyy_1 = (float *) malloc( size * sizeof(float) );
    if (NULL==kin_histyy_1) 
      error("Cannot allocate yy kinetic energy tensor array.");
  }
  if (NULL==kin_histxxu_1) {
    kin_histxxu_1 = (float *) malloc( size * sizeof(float) );
    if (NULL==kin_histxxu_1) 
      error("Cannot allocate xxu kinetic energy tensor array.");
  }
#ifndef TWOD  
  if (NULL==kin_histzz_1) {
    kin_histzz_1 = (float *) malloc( size * sizeof(float) );
    if (NULL==kin_histzz_1) 
      error("Cannot allocate zz pressure tensor array.");
  }
#endif
#endif /* SHOCK */

  if (NULL==kin_hist_1) {
    kin_hist_1 = (float *) malloc( size * sizeof(float) );
    if (NULL==kin_hist_1) 
      error("Cannot allocate kinetic energy array.");
  }
  if (NULL==pot_hist_1) {
    pot_hist_1 = (float *) malloc( size * sizeof(float) );
    if (NULL==pot_hist_1) 
      error("Cannot allocate potential energy array.");
  }  
  if (NULL==num_hist_1) {
    num_hist_1 = (integer *) malloc( size * sizeof(integer) );
    if (NULL==num_hist_1) 
      error("Cannot allocate number count array.");
  }

#ifdef MPI

#ifdef STRESS_TENS
  if (NULL==press_histxx_2) {
    press_histxx_2 = (float *) malloc( size * sizeof(float) );
    if ((NULL==press_histxx_2) && (myid==0)) 
      error("Cannot allocate xx pressure tensor  array.");
  }  
  if (NULL==press_histyy_2) {
    press_histyy_2 = (float *) malloc( size * sizeof(float) );
    if ((NULL==press_histyy_2) && (myid==0)) 
      error("Cannot allocate yy pressure tensor  array.");
  }
#ifndef SHOCK
  if (NULL==press_histxy_2) {
    press_histxy_2 = (float *) malloc( size * sizeof(float) );
    if ((NULL==press_histxy_2) && (myid==0)) 
      error("Cannot allocate xy pressure tensor  array.");
  }
#endif
#ifndef TWOD
  if (NULL==press_histzz_2) {
    press_histzz_2 = (float *) malloc( size * sizeof(float) );
    if ((NULL==press_histzz_2) && (myid==0)) 
      error("Cannot allocate zz pressure tensor  array.");
  }
#ifndef SHOCK
  if (NULL==press_histzx_2) {
    press_histzx_2 = (float *) malloc( size * sizeof(float) );
    if ((NULL==press_histzx_2) && (myid==0)) 
      error("Cannot allocate zx pressure tensor  array.");
  }
  if (NULL==press_histyz_2) {
    press_histyz_2 = (float *) malloc( size * sizeof(float) );
    if ((NULL==press_histyz_2) && (myid==0)) 
      error("Cannot allocate yz pressure tensor  array.");
  }
#endif /* SHOCK */
#endif /* not TWOD */
#endif /* STRESS_TENS */

#ifdef SHOCK
  if (NULL==kin_histxx_2) {
    kin_histxx_2 = (float *) malloc( size * sizeof(float) );
    if (NULL==kin_histxx_2) 
      error("Cannot allocate xx kinetic energy tensor array.");
  }  
  if (NULL==kin_histyy_2) {
    kin_histyy_2 = (float *) malloc( size * sizeof(float) );
    if (NULL==kin_histyy_2) 
      error("Cannot allocate yy kinetic energy tensor array.");
  }
  if (NULL==kin_histxxu_2) {
    kin_histxxu_2 = (float *) malloc( size * sizeof(float) );
    if (NULL==kin_histxxu_2) 
      error("Cannot allocate xxu kinetic energy tensor array.");
  }
#ifndef TWOD  
  if (NULL==kin_histzz_2) {
    kin_histzz_2 = (float *) malloc( size * sizeof(float) );
    if (NULL==kin_histzz_2) 
      error("Cannot allocate zz pressure tensor array.");
  }
#endif
#endif /* SHOCK */

  if (NULL==kin_hist_2) {
    kin_hist_2 = (float *) malloc( size * sizeof(float) );
    if (NULL==kin_hist_2) 
      error("Cannot allocate kinetic energy array.");
  }  
  if (NULL==pot_hist_2) {
    pot_hist_2 = (float *) malloc( size * sizeof(float) );
    if (NULL==pot_hist_2) 
      error("Cannot allocate potential energy array.");
  }  
  if (NULL==num_hist_2) {
    num_hist_2 = (integer *) malloc( size * sizeof(integer) );
    if ((NULL==num_hist_2) && (myid==0)) 
      error("Cannot allocate number count array.");
  }

#endif /* MPI */

  /* clear histograms */
  for (i = 0; i < size; i++) {

#ifdef STRESS_TENS
#ifdef SHOCK
    press_histxx_1[i] = 0.0;
    press_histyy_1[i] = 0.0;
    kin_histxx_1  [i] = 0.0;
    kin_histyy_1  [i] = 0.0;
    kin_histxxu_1 [i] = 0.0;
#ifndef TWOD
    press_histzz_1[i] = 0.0;
    kin_histzz_1  [i] = 0.0;
#endif
#else /* not SHOCK */
    press_histxx_1[i] = 0.0;
    press_histyy_1[i] = 0.0;
    press_histxy_1[i] = 0.0;
#ifndef TWOD
    press_histzz_1[i] = 0.0;
    press_histzx_1[i] = 0.0;
    press_histyz_1[i] = 0.0;
#endif
#endif /* SHOCK */
#endif /* STRESS_TENS */
    kin_hist_1[i] = 0; 
    pot_hist_1[i] = 0;
    num_hist_1[i] = 0;

#ifdef MPI
#ifdef STRESSTENS
#ifdef SHOCK
    press_histxx_2[i] = 0.0;
    press_histyy_2[i] = 0.0;
    kin_histxx_2  [i] = 0.0;
    kin_histyy_2  [i] = 0.0;
    kin_histxxu_2 [i] = 0.0;
#ifndef TWOD
    press_histzz_2[i] = 0.0;
    kin_histzz_2  [i] = 0.0;
#endif
#else /* not SHOCK */
    press_histxx_2[i] = 0.0;
    press_histyy_2[i] = 0.0;
    press_histxy_2[i] = 0.0;
#ifndef TWOD
    press_histzz_2[i] = 0.0;
    press_histzx_2[i] = 0.0;
    press_histyz_2[i] = 0.0;
#endif
#endif /* SHOCK */
#endif /* STRESS_TENS */
    kin_hist_2[i] = 0;
    pot_hist_2[i] = 0;
    num_hist_2[i] = 0;
#endif /* MPI */

  }

  /* loop over all atoms */
  for (k=0; k<ncells; ++k) {

    p = cell_array + CELLS(k);

    for (i = 0; i < p->n; ++i) {
      /* which bin? */
      numx = scalex * p->ort X(i);
      if (numx < 0)             numx = 0;
      if (numx >= hist->dim.x)  numx = hist->dim.x-1;
      numy = scaley * p->ort Y(i);
      if (numy < 0)             numy = 0;
      if (numy >= hist->dim.y)  numy = hist->dim.y-1;
      num = numx * hist->dim.y + numy;
#ifndef TWOD
      numz = scalez * p->ort Z(i);
      if (numz < 0)             numz = 0;
      if (numz >= hist->dim.z)  numz = hist->dim.z-1;
      num = num * hist->dim.z + numz;
#endif

#ifdef STRESS_TENS

#ifdef SHOCK
      press_histxx_1[num] += p->presstens X(i);
      press_histyy_1[num] += p->presstens Y(i);
      kin_histxx_1  [num] += SQR(p->impuls X(i)) / (2*MASSE(p,i));
      kin_histyy_1  [num] += SQR(p->impuls Y(i)) / (2*MASSE(p,i));
      /* average v_xx - u_p  relative to moving pistons */
      tmp = shock_speed * MASSE(p,i);
      /* plate against bulk */
      if (shock_mode == 1) {
        if ( p->ort X(i) < shock_strip ) 
          kin_histxxu_1[num] += SQR(p->impuls X(i) - tmp) / (2*MASSE(p,i));
        else
          kin_histxxu_1[num] += SQR(p->impuls X(i)) / (2*MASSE(p,i));
      }
      /* two halves against one another */
      if (shock_mode == 2) {
        if ( p->ort X(i) < box_x.x*0.5 )
          kin_histxxu_1[num] += SQR(p->impuls X(i) - tmp) / (2*MASSE(p,i));
        else
          kin_histxxu_1[num] += SQR(p->impuls X(i) + tmp) / (2*MASSE(p,i));
      }
      /* bulk against wall */
      if (shock_mode == 3) 
          kin_histxxu_1[num] += SQR(p->impuls X(i) - tmp) / (2*MASSE(p,i));

#ifndef TWOD
      press_histzz_1[num] += p->presstens Z(i);
      kin_histzz_1  [num] += SQR(p->impuls Z(i)) / (2*MASSE(p,i));
#endif

#else /* not SHOCK */

      press_histxx_1[num] += p->presstens X(i);
      press_histyy_1[num] += p->presstens Y(i);
#ifdef TWOD
      press_histxy_1[num] += p->presstens_offdia[i];
#else
      press_histzz_1[num] += p->presstens Z(i);
      press_histxy_1[num] += p->presstens_offdia Z(i);
      press_histzx_1[num] += p->presstens_offdia Y(i);
      press_histyz_1[num] += p->presstens_offdia X(i);
#endif
#endif /* SHOCK */

#endif /* STRESS_TENS */

      Ekin = SPRODN(p->impuls,i,p->impuls,i) / (2* MASSE(p,i));
      kin_hist_1[num] += Ekin; 
#ifdef DISLOC
      if (Epot_diff==1)
        pot_hist_1[num] += p->pot_eng[i] - p->Epot_ref[i];
      else
#endif
#if defined(ORDPAR) && !defined(TWOD)
      pot_hist_1[num] += (p->nbanz[i]==0)?0:p->pot_eng[i]/p->nbanz[i];
#else
      pot_hist_1[num] += POTENG(p,i);
#endif
      num_hist_1[num]++;

    }
  }

  
#ifdef MPI  /* add up results form different CPUs */

#ifdef STRESS_TENS
  MPI_Reduce(press_histxx_1,press_histxx_2,size,MPI_FLOAT,MPI_SUM,0,cpugrid);
  hist->press_histxx = press_histxx_2;
  MPI_Reduce(press_histyy_1,press_histyy_2,size,MPI_FLOAT,MPI_SUM,0,cpugrid);
  hist->press_histyy = press_histyy_2;
#ifndef SHOCK
  MPI_Reduce(press_histxy_1,press_histxy_2,size,MPI_FLOAT,MPI_SUM,0,cpugrid);
  hist->press_histxy = press_histxy_2;
#endif
#ifndef TWOD
  MPI_Reduce(press_histzz_1,press_histzz_2,size,MPI_FLOAT,MPI_SUM,0,cpugrid);
  hist->press_histzz = press_histzz_2;
#ifndef SHOCK
  MPI_Reduce(press_histzx_1,press_histzx_2,size,MPI_FLOAT,MPI_SUM,0,cpugrid);
  hist->press_histzx = press_histzx_2;
  MPI_Reduce(press_histyz_1,press_histyz_2,size,MPI_FLOAT,MPI_SUM,0,cpugrid);
  hist->press_histyz = press_histyz_2;
#endif
#endif
#endif /* STRESS_TENS */

#ifdef SHOCK
  MPI_Reduce(kin_histxx_1, kin_histxx_2, size,MPI_FLOAT,MPI_SUM,0,cpugrid);
  hist->kin_histxx  = kin_histxx_2;
  MPI_Reduce(kin_histyy_1, kin_histyy_2, size,MPI_FLOAT,MPI_SUM,0,cpugrid);
  hist->kin_histyy  = kin_histyy_2;
  MPI_Reduce(kin_histxxu_1,kin_histxxu_2,size,MPI_FLOAT,MPI_SUM,0,cpugrid);
  hist->kin_histxxu = kin_histxxu_2;
#ifndef TWOD
  MPI_Reduce(kin_histzz_1, kin_histzz_2, size,MPI_FLOAT,MPI_SUM,0,cpugrid);
  hist->kin_histzz  = kin_histzz_2;
#endif
#endif /* SHOCK */

  MPI_Reduce( kin_hist_1, kin_hist_2, size, MPI_FLOAT, MPI_SUM, 0, cpugrid);
  hist->kin_hist     = kin_hist_2;
  MPI_Reduce( pot_hist_1, pot_hist_2, size, MPI_FLOAT, MPI_SUM, 0, cpugrid);
  hist->pot_hist     = pot_hist_2;
  MPI_Reduce( num_hist_1, num_hist_2, size, INTEGER,   MPI_SUM, 0, cpugrid);
  hist->num_hist     = num_hist_2;

#else /* not MPI */

#ifdef STRESS_TENS
#ifdef SHOCK
  hist->press_histxx = press_histxx_1;
  hist->press_histyy = press_histyy_1;
  hist->kin_histxx   = kin_histxx_1;
  hist->kin_histyy   = kin_histyy_1;
  hist->kin_histxxu  = kin_histxxu_1;
#ifndef TWOD
  hist->press_histzz = press_histzz_1;
  hist->kin_histzz   = kin_histzz_1;
#endif
#else /* not SHOCK */
  hist->press_histxx = press_histxx_1;
  hist->press_histyy = press_histyy_1;
  hist->press_histxy = press_histxy_1;
#ifndef TWOD
  hist->press_histzz = press_histzz_1;
  hist->press_histzx = press_histzx_1;
  hist->press_histyz = press_histyz_1;
#endif
#endif /* SHOCK */
#endif /* STRESS_TENS */
  hist->kin_hist     = kin_hist_1; 
  hist->pot_hist     = pot_hist_1;
  hist->num_hist     = num_hist_1;

#endif /* not MPI */

  /* normalize histograms */
  if (myid==0) {
    for (i=0; i < size; i++) {
      if (hist->num_hist[i] > 0) {
#ifdef STRESS_TENS
#ifdef SHOCK
	hist->press_histxx[i] /= hist->num_hist[i];
	hist->press_histyy[i] /= hist->num_hist[i];
	hist->kin_histxx  [i] /= hist->num_hist[i];
	hist->kin_histyy  [i] /= hist->num_hist[i];
	hist->kin_histxxu [i] /= hist->num_hist[i];
#ifndef TWOD
	hist->press_histzz[i] /= hist->num_hist[i];
	hist->kin_histzz  [i] /= hist->num_hist[i];
#endif
#else /* not SHOCK */
	hist->press_histxx[i] /= hist->num_hist[i];
	hist->press_histyy[i] /= hist->num_hist[i];
	hist->press_histxy[i] /= hist->num_hist[i];
#ifndef TWOD
	hist->press_histzz[i] /= hist->num_hist[i];
	hist->press_histzx[i] /= hist->num_hist[i];
	hist->press_histyz[i] /= hist->num_hist[i];
#endif
#endif /* SHOCK */
#endif /* STRESS_TENS */
	hist->kin_hist[i]     /= hist->num_hist[i];
	hist->pot_hist[i]     /= hist->num_hist[i];
      }
    }
  } 

  /* compute minima and maxima of histograms */
  if (myid==0) {

    j=0;
    while (hist->num_hist[j]==0) j++;

#ifdef STRESS_TENS

#ifdef SHOCK
    hist->maxpxx  = hist->press_histxx[j];
    hist->minpxx  = hist->press_histxx[j];
    hist->maxpyy  = hist->press_histyy[j];
    hist->minpyy  = hist->press_histyy[j];
    hist->maxkxx  = hist->kin_histxx  [j];
    hist->minkxx  = hist->kin_histxx  [j];
    hist->maxkyy  = hist->kin_histyy  [j];
    hist->minkyy  = hist->kin_histyy  [j];
    hist->maxkxxu = hist->kin_histxxu [j];
    hist->minkxxu = hist->kin_histxxu [j];
#ifndef TWOD
    hist->maxpzz  = hist->press_histzz[j];
    hist->minpzz  = hist->press_histzz[j];
    hist->maxkzz  = hist->kin_histzz  [j];
    hist->minkzz  = hist->kin_histzz  [j];
#endif
#else /* not SHOCK */
    hist->maxpxx  = hist->press_histxx[j];
    hist->minpxx  = hist->press_histxx[j];
    hist->maxpyy  = hist->press_histyy[j];
    hist->minpyy  = hist->press_histyy[j];
    hist->maxpxy  = hist->press_histxy[j];
    hist->minpxy  = hist->press_histxy[j];
#ifndef TWOD
    hist->maxpzz  = hist->press_histzz[j];
    hist->minpzz  = hist->press_histzz[j];
    hist->maxpyz  = hist->press_histyz[j];
    hist->minpyz  = hist->press_histyz[j];
    hist->maxpzx  = hist->press_histzx[j];
    hist->minpzx  = hist->press_histzx[j];
#endif
#endif /* SHOCK */

#endif /* STRESS_TENS */

    hist->minpot  = hist->pot_hist[j];
    hist->maxpot  = hist->pot_hist[j];
    hist->minkin  = hist->kin_hist[j];
    hist->maxkin  = hist->kin_hist[j];

    for (i=j+1; i<size; i++) {
      if (hist->num_hist[i]>0) {

#ifdef STRESS_TENS
#ifdef SHOCK
        hist->maxpxx  = MAX( hist->maxpxx,  hist->press_histxx[i] );
        hist->minpxx  = MIN( hist->minpxx,  hist->press_histxx[i] );
        hist->maxpyy  = MAX( hist->maxpyy,  hist->press_histyy[i] );
        hist->minpyy  = MIN( hist->minpyy,  hist->press_histyy[i] );
        hist->maxkxx  = MAX( hist->maxkxx,  hist->kin_histxx  [i] );
        hist->minkxx  = MIN( hist->minkxx,  hist->kin_histxx  [i] );
        hist->maxkyy  = MAX( hist->maxkyy,  hist->kin_histyy  [i] );
        hist->minkyy  = MIN( hist->minkyy,  hist->kin_histyy  [i] );
        hist->maxkxxu = MAX( hist->maxkxxu, hist->kin_histxxu [i] );
        hist->minkxxu = MIN( hist->minkxxu, hist->kin_histxxu [i] );
#ifndef TWOD
        hist->maxpzz  = MAX( hist->maxpzz,  hist->press_histzz[i] );
        hist->minpzz  = MIN( hist->minpzz,  hist->press_histzz[i] );
#endif
#else /* not SHOCK */
        hist->maxpxx  = MAX( hist->maxpxx,  hist->press_histxx[i] );
        hist->minpxx  = MIN( hist->minpxx,  hist->press_histxx[i] );
        hist->maxpyy  = MAX( hist->maxpyy,  hist->press_histyy[i] );
        hist->minpyy  = MIN( hist->minpyy,  hist->press_histyy[i] );
        hist->maxpxy  = MAX( hist->maxpxy,  hist->press_histxy[i] );
        hist->minpxy  = MIN( hist->minpxy,  hist->press_histxy[i] );
#ifndef TWOD
        hist->maxpzz  = MAX( hist->maxpzz,  hist->press_histzz[i] );
        hist->minpzz  = MIN( hist->minpzz,  hist->press_histzz[i] );
        hist->maxpyz  = MAX( hist->maxpyz,  hist->press_histyz[i] );
        hist->minpyz  = MIN( hist->minpyz,  hist->press_histyz[i] );
        hist->maxpzx  = MAX( hist->maxpzx,  hist->press_histzx[i] );
        hist->minpzx  = MIN( hist->minpzx,  hist->press_histzx[i] );
#endif
#endif /* SHOCK */
#endif /* STRESS_TENS */

        hist->maxpot  = MAX( hist->maxpot,  hist->pot_hist[i] );
        hist->minpot  = MIN( hist->minpot,  hist->pot_hist[i] );
        hist->maxkin  = MAX( hist->maxkin,  hist->kin_hist[i] );
        hist->minkin  = MIN( hist->minkin,  hist->kin_hist[i] );
      }
    }
  }

}


/******************************************************************************
*
* write_distrib write spatial distribution of potential and kinetice energy
*
******************************************************************************/

void write_distrib(int steps)
{
  FILE   *outpot, *outkin, *outminmax;
  str255 fnamepot, fnamekin, fnameminmax;
  hist_t hist;
  int    fzhlr, i, r, s, t, count_pot, count_kin;

  hist.dim = dist_dim;
  make_histograms(&hist);

  if (0==myid) {

    fzhlr = steps / dis_interval;
    sprintf(fnamepot,"%s.%u.pot.dist",outfilename,fzhlr);
    sprintf(fnamekin,"%s.%u.kin.dist",outfilename,fzhlr);
    sprintf(fnameminmax,"%s.minmax.dist",outfilename);

    outpot = fopen(fnamepot,"w");
    if (NULL == outpot) error("Cannot open pot distrib file.");
    outkin = fopen(fnamekin,"w");
    if (NULL == outkin) error("Cannot open kin distrib file.");

    outminmax = fopen(fnameminmax, "a");
    fprintf(outminmax, "%d %f %f %f %f\n", 
            fzhlr, hist.minpot, hist.maxpot, hist.minkin, hist.maxkin);
    fclose(outminmax);

    if (dist_binary_io) {
      count_pot=fwrite(hist.pot_hist, sizeof(float), hist.size, outpot);
      count_kin=fwrite(hist.kin_hist, sizeof(float), hist.size, outkin);
      if ((count_pot!=hist.size) || (count_kin!=hist.size)) {
        char msg[255];
        sprintf(msg,"dist write incomplete - cnt_pot = %d, cnt_kin = %d", 
                count_pot, count_kin );
        warning(msg);
      }
    } else {
      i=0;
      for (r=0; r<hist.dim.x; r++) {
        for (s=0; s<hist.dim.y; s++) {
#ifndef TWOD
          for (t=0; t<hist.dim.z; t++) {
#endif
            fprintf(outpot,"%f\n", hist.pot_hist[i]);
            fprintf(outkin,"%f\n", hist.kin_hist[i]);
            i++;
#ifndef TWOD
          }
          fprintf(outpot,"\n");
          fprintf(outkin,"\n");
#endif
	}
        fprintf(outpot,"\n");
        fprintf(outkin,"\n");
      }
    }
    fclose(outpot);
    fclose(outkin);
  }
}


#ifdef STRESS_TENS

#ifndef SHOCK

/******************************************************************************
*
* write_press_dist (distributions of energy and stress tensor)
*
******************************************************************************/

void write_press_dist(int steps)
{
  FILE   *outfile;
  str255 fname;
  hist_t hist;
  int    fzhlr, i, r, s, t;

  hist.dim = press_dim;
  make_histograms(&hist);

  if (0==myid) {

    fzhlr = steps / press_interval;
    sprintf(fname,"%s.%u.pressdist",outfilename,fzhlr);
    outfile = fopen(fname,"w");
    if (NULL == outfile) error("Cannot open pressure tensor file.");

    i=0;
    for (r=0; r<hist.dim.x; r++) {
      for (s=0; s<hist.dim.y; s++) {
#ifdef TWOD
        fprintf(outfile, "%e %e %e %e %e %e\n", 
	        hist.press_histxx[i], hist.press_histyy[i],
                hist.press_histxy[i], hist.num_hist[i] / hist.binvol,
                hist.kin_hist[i], hist.pot_hist[i] );
        i++;
#else
        for (t=0; t<hist.dim.z; t++) {
          fprintf(outfile, "%e %e %e %e %e %e %e %e %e\n", 
	          hist.press_histxx[i], hist.press_histyy[i], 
                  hist.press_histzz[i], hist.press_histyz[i], 
                  hist.press_histzx[i], hist.press_histxy[i],
                  hist.num_hist[i] / hist.binvol, 
                  hist.kin_hist[i], hist.pot_hist[i] );
          i++;
	}
        fprintf(outfile,"\n");
#endif
      }
      fprintf(outfile,"\n");
    }
    fclose(outfile);
  }
}


#else /* SHOCK */


/******************************************************************************
*
* write_press_dist_shock (distribution of stress tensor for shock waves)
*
******************************************************************************/

void write_press_dist_shock(int steps)
{
  FILE   *outfile;
  str255 fname;
  hist_t hist;
  int    fzhlr, i;

  /* we make only a 1D histogram */
  hist.dim.x = press_dim.x;
  hist.dim.y = 1;
#ifndef TWOD
  hist.dim.z = 1;
#endif

  make_histograms(&hist);

  /* write pressure tensor distribution */

  if (myid==0) {

    fzhlr = steps / press_interval;

    sprintf(fname,"%s.%u.pressdist",outfilename,fzhlr);
    outfile = fopen(fname,"w");
    if (NULL == outfile) error("Cannot open pressure tensor file.");

    for (i = 0; i < hist.size; i++) {
#ifndef TWOD
      fprintf(outfile,
        "%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10e\n", 
        hist.press_histxx[i], hist.press_histyy[i], hist.press_histzz[i], 
        hist.num_hist[i] / hist.binvol, 
        hist.kin_histxx[i], hist.kin_histxxu[i], hist.kin_histyy[i],
        hist.kin_histzz[i], hist.pot_hist[i] );
#else
      fprintf(outfile,
        "%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n", 
        hist.press_histxx[i], hist.press_histyy[i], 
        hist.num_hist[i] / hist.binvol,
        hist.kin_histxx[i], hist.kin_histxxu[i], hist.kin_histyy[i], 
        hist.pot_hist[i] ); 
#endif
    }
    fprintf(outfile,"\n");
    fclose(outfile);
  }
}

#endif /* SHOCK */

#endif /* STRESS_TENS */
