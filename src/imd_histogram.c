
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
  real scalex, scaley, scalez, Ekin, tmp, binvol;
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
  scalex = hist->dim.x / (hist->ur.x - hist->ll.x);
  scaley = hist->dim.y / (hist->ur.y - hist->ll.y);
  binvol = 1.0 / (scalex * scaley);
  size   = hist->dim.x * hist->dim.y;
#ifndef TWOD
  scalez = hist->dim.z / (hist->ur.z - hist->ll.z);
  binvol = binvol / scalez;
  size  *= hist->dim.z;
#endif
  hist->binvol = binvol;
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
  for (i=0; i<size; i++) {

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

    for (i=0; i<p->n; ++i) {

      /* which bin? */
      numx = scalex * (p->ort X(i) - hist->ll.x);
      if ((numx < 0) || (numx >= hist->dim.x)) continue;
      numy = scaley * (p->ort Y(i) - hist->ll.y);
      if ((numy < 0) || (numy >= hist->dim.y)) continue;
      num = numx * hist->dim.y + numy;
#ifndef TWOD
      numz = scalez * (p->ort Z(i) - hist->ll.z);
      if ((numz < 0) || (numz >= hist->dim.z)) continue;
      num = num * hist->dim.z + numz;
#endif

#ifdef STRESS_TENS

#ifdef SHOCK
      press_histxx_1[num] += p->presstens[i].xx;
      press_histyy_1[num] += p->presstens[i].yy;
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
      press_histzz_1[num] += p->presstens[i].zz;
      kin_histzz_1  [num] += SQR(p->impuls Z(i)) / (2*MASSE(p,i));
#endif

#else /* not SHOCK */

      press_histxx_1[num] += p->presstens[i].xx;
      press_histyy_1[num] += p->presstens[i].yy;
#ifndef TWOD
      press_histzz_1[num] += p->presstens[i].zz;
      press_histzx_1[num] += p->presstens[i].zx;
      press_histyz_1[num] += p->presstens[i].yz;
#endif
      press_histxy_1[num] += p->presstens[i].xy;
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
*  write header for histograms of potential and kinetic energy
*
******************************************************************************/

void write_distrib_header(FILE *out, hist_t *hist, char *type)
{
  char c;
  int  n_coord;
  time_t now;
  vektor s;

  /* format line -- format dim n_coord n_data */
  if (dist_binary_io)
    c = is_big_endian ? 'B' : 'L';
  else
    c = 'A';
  n_coord = dist_has_coords ? DIM : 0;
  fprintf(out, "#F %c %d %d 1\n", c, DIM, n_coord);

  /* contents line */
  if (dist_has_coords)
#ifdef TWO
    fprintf(out, "#C x y %s\n",   type);
#else
    fprintf(out, "#C x y z %s\n", type);
#endif
  else
    fprintf(out, "#C %s\n",       type);

  /* dimension line */
#ifdef TWOD
  fprintf(out, "#D %d %d\n",    hist->dim.x, hist->dim.y);
#else
  fprintf(out, "#D %d %d %d\n", hist->dim.x, hist->dim.y, hist->dim.z);
#endif

  /* bin size line */
  s.x = (hist->ur.x - hist->ll.x) / hist->dim.x;
  s.y = (hist->ur.y - hist->ll.y) / hist->dim.y;
#ifdef TWOD
  fprintf(out, "#S %e %e\n",    s.x, s.y);
#else
  s.z = (hist->ur.z - hist->ll.z) / hist->dim.z;
  fprintf(out, "#S %e %e %e\n", s.x, s.y, s.z);
#endif

  /* endheader line */
  time(&now);
  fprintf(out, "## Generated by %s on %s", progname, ctime(&now));
  fprintf(out, "#E\n");

}

/******************************************************************************
*
*  write histograms of potential and kinetic energy
*
******************************************************************************/

void write_distrib(int steps)
{
  FILE   *outpot, *outkin, *outminmax;
  str255 fnamepot, fnamekin, fnameminmax;
  hist_t hist;
  int    fzhlr, i, r, s, t, count_pot, count_kin;
  unsigned char c;
  char d;

  hist.dim = dist_dim;
  hist.ll  = hist_ll;
  hist.ur  = hist_ur;
  make_histograms(&hist);

  if (0==myid) {

    fzhlr = steps / dis_interval;
    if (virvo_io==1) {
      sprintf(fnamepot,"%s.%u.pot.rvf",outfilename,fzhlr);
      sprintf(fnamekin,"%s.%u.kin.rvf",outfilename,fzhlr);
    } else {
      sprintf(fnamepot,"%s.%u.pot.dist",outfilename,fzhlr);
      sprintf(fnamekin,"%s.%u.kin.dist",outfilename,fzhlr);
      sprintf(fnameminmax,"%s.minmax.dist",outfilename);
    }

    outpot = fopen(fnamepot,"w");
    if (NULL == outpot) error("Cannot open pot distrib file.");
    outkin = fopen(fnamekin,"w");
    if (NULL == outkin) error("Cannot open kin distrib file.");

    outminmax = fopen(fnameminmax, "a");
    if ((virvo_io==0)&&(use_header)) {
      write_distrib_header(outpot, &hist, "Epot");
      write_distrib_header(outkin, &hist, "Ekin");
      fprintf(outminmax, 
              "# contents count min_Epot max_Epot min_Ekin max_Ekin\n");
    }
    if (virvo_io==1) {
      c = (unsigned char)((hist.dim.x & 65280)>>8);
      fprintf(outpot, "%c", c);
      fputc(c, outkin);
      c = (unsigned char)(hist.dim.x & 255); 
      fprintf(outpot, "%c", c);
      fputc(c, outkin);
      c = (unsigned char)((hist.dim.y & 65280)>>8);
      fprintf(outpot, "%c", c);
      fputc(c, outkin);
      c = (unsigned char)(hist.dim.y & 255); 
      fprintf(outpot, "%c", c);
      fputc(c, outkin);
#ifndef TWOD
      c = (unsigned char)((hist.dim.z & 65280)>>8);
      fprintf(outpot, "%c", c);
      fputc(c, outkin);
      c = (unsigned char)(hist.dim.z & 255); 
      fprintf(outpot, "%c", c);
      fputc(c, outkin);
#endif
    }

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
	    if (dist_has_coords) {
	      fprintf(outpot, "%d %d ", r, s);
	      fprintf(outkin, "%d %d ", r, s);
#ifndef TWOD
	      fprintf(outpot, "%d ", t);
	      fprintf(outkin, "%d ", t);
#endif
	    }
	    if (virvo_io==1) {
	      c = (unsigned char)(256.0*(hist.pot_hist[i]-hist.minpot)/(hist.maxpot-hist.minpot));
	      fputc(c, outpot);
	      c = (unsigned char)(256.0*(hist.kin_hist[i]-hist.minkin)/(hist.maxkin-hist.minkin));
	      fputc(c, outkin);
	      i++;
	    } else {
	      if (norm_hist) {
		fprintf(outpot,"%f\n", 
			(hist.num_hist[i])?hist.pot_hist[i]/hist.num_hist[i]:-10000.0);
		fprintf(outkin,"%f\n", 
			(hist.num_hist[i])?hist.kin_hist[i]/hist.num_hist[i]:0.0);
	      } else {
		fprintf(outpot,"%f\n", hist.pot_hist[i]);
		fprintf(outkin,"%f\n", hist.kin_hist[i]);
	      }
	      i++;
	    }
#ifndef TWOD
	  }
	  if (virvo_io==0) {
	    fprintf(outpot,"\n");
	    fprintf(outkin,"\n");
	  }
#endif
	}
	if (virvo_io==0) {
	  fprintf(outpot,"\n");
	  fprintf(outkin,"\n");
	}
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
*  write header for histograms of pressure tensor
*
******************************************************************************/

void write_press_dist_header(FILE *out, hist_t *hist)
{
  char c;
  int n_coord;
  time_t now;
  vektor s;

  /* format line -- format dim n_coord n_data */
  /*
  if (dist_binary_io)
    c = is_big_endian ? 'B' : 'L';
  else
  */
    c = 'A';
  n_coord = dist_has_coords ? DIM : 0;
  fprintf(out, "#F %c %d %d %d\n", c, n_coord, (DIM*(DIM+1))/2);

  /* contents line */
#ifdef TWOD
  if (dist_has_coords)
    fprintf(out, "#C x y P_xx P_yy P_xy density Ekin Epot\n");
  else
    fprintf(out, "#C P_xx P_yy P_xy density Ekin Epot\n");
#else
  if (dist_has_coords)
    fprintf(out, "#C x y z P_xx P_yy P_zz P_yz P_zy P_xy density Ekin Epot\n");
  else
    fprintf(out, "#C P_xx P_yy P_zz P_yz P_zx P_xy density Ekin Epot\n");
#endif

  /* dimension line */
#ifdef TWOD
  fprintf(out, "#D %d %d\n",    hist->dim.x, hist->dim.y);
#else
  fprintf(out, "#D %d %d %d\n", hist->dim.x, hist->dim.y, hist->dim.z);
#endif

  /* bin size line */
  s.x = (hist->ur.x - hist->ll.x) / hist->dim.x;
  s.y = (hist->ur.y - hist->ll.y) / hist->dim.y;
#ifdef TWOD
  fprintf(out, "#S %e %e\n",    s.x, s.y);
#else
  s.z = (hist->ur.z - hist->ll.z) / hist->dim.z;
  fprintf(out, "#S %e %e %e\n", s.x, s.y, s.z);
#endif

  /* endheader line */
  time(&now);
  fprintf(out, "## Generated by %s on %s", progname, ctime(&now));
  fprintf(out, "#E\n");

}

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
  hist.ur  = hist_ur;
  hist.ll  = hist_ll;
  make_histograms(&hist);

  if (0==myid) {

    fzhlr = steps / press_interval;
    sprintf(fname,"%s.%u.pressdist",outfilename,fzhlr);
    outfile = fopen(fname,"w");
    if (NULL == outfile) error("Cannot open pressure tensor file.");
    write_press_dist_header(outfile, &hist);

    i=0;
    for (r=0; r<hist.dim.x; r++) {
      for (s=0; s<hist.dim.y; s++) {
#ifdef TWOD
        if (dist_has_coords) fprintf(outfile, "%d %d ", r, s);
        fprintf(outfile, "%e %e %e %e %e %e\n", 
	        hist.press_histxx[i], hist.press_histyy[i],
q                hist.press_histxy[i], hist.num_hist[i] / hist.binvol,
                hist.kin_hist[i], hist.pot_hist[i] );
        i++;
#else
        for (t=0; t<hist.dim.z; t++) {
          if (dist_has_coords) fprintf(outfile, "%d %d %d ", r, s, t);
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
*  write header for histograms of pressure tensor (shock mode)
*
******************************************************************************/

void write_press_dist_shock_header(FILE *out, hist_t *hist)
{
  char c;
  int n_coord;
  time_t now;

  /* format line -- format dim n_coord n_data */
  /*
  if (dist_binary_io)
    c = is_big_endian ? 'B' : 'L';
  else
  */
    c = 'A';
  n_coord = dist_has_coords ? 1 : 0;
  fprintf(out, "#F %c 1 %d %d\n", c, 1, 2*DIM+3);

  /* contents line */
#ifdef TWOD
  if (dist_has_coords) fprintf(out, "#C x P_xx P_yy density");
  else                 fprintf(out, "#C P_xx P_yy density");
  fprintf(out, " Ekin_xx Ekin_xxu Ekin_yy Epot\n");
#else
  if (dist_has_coords) fprintf(out, "#C x P_xx P_yy P_zz density");
  else                 fprintf(out, "#C P_xx P_yy P_zz density");
  fprintf(out, " Ekin_xx Ekin_xxu Ekin_yy Ekin_zz Epot\n");
#endif

  /* dimension line - we make only a 1-dim histogram */
  fprintf(out, "#D %d\n", hist->dim.x);

  /* bin size line */
  fprintf(out, "#S %e\n", (hist->ur.x - hist->ll.x) / hist->dim.x);

  /* endheader line */
  time(&now);
  fprintf(out, "## Generated by %s on %s", progname, ctime(&now));
  fprintf(out, "#E\n");

}

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
  hist.ur    = hist_ur;
  hist.ll    = hist_ll;
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
    write_press_dist_shock_header(outfile, &hist);

    for (i = 0; i < hist.size; i++) {
      if (dist_has_coords) fprintf(outfile, "%d ", i);      
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

#ifdef ATDIST

/******************************************************************************
*
*  initialize atoms distribution array
*
******************************************************************************/
  
void init_atoms_dist()
{
  int size, i;

  /* compute array size */
  atoms_dist_size  = atoms_dist_dim.x * atoms_dist_dim.y;
#ifndef TWOD
  atoms_dist_size *= atoms_dist_dim.z;
#endif
  size = atoms_dist_size * ntypes;

  /* backup if pic_ur is not set */
  if (0.0==pic_ur.x) {
    pic_ur.x = box_x.x;
    pic_ur.y = box_y.y;
#ifndef TWOD
    pic_ur.z = box_z.z;
#endif
  }

  /* the bins are orthogonal boxes in space */
  atoms_dist_scale.x = atoms_dist_dim.x / (pic_ur.x - pic_ll.x);
  atoms_dist_scale.y = atoms_dist_dim.y / (pic_ur.y - pic_ll.y);
#ifndef TWOD
  atoms_dist_scale.z = atoms_dist_dim.z / (pic_ur.z - pic_ll.z);
#endif

  /* allocate histogram array */
  if (NULL==atoms_dist) {
    atoms_dist = (float *) malloc( size * sizeof(float) );
    if (NULL==atoms_dist) error("Cannot allocate atoms distribution array.");
  }

  /* initialize histogram array */
  for (i=0; i<size; i++) atoms_dist[i]=0.0;
}

/******************************************************************************
*
*  update atoms distribution array
*
******************************************************************************/
  
void update_atoms_dist()
{
  int  num, numx, numy, numz;
  cell *p;
  int  i, k;

  /* loop over all atoms */
  for (k=0; k<ncells; ++k) {
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) {
      /* continue if atom is not inside selected box */
      if ((p->ort X(i) < pic_ll.x) || (p->ort X(i) > pic_ur.x) ||
#ifndef TWOD
          (p->ort Z(i) < pic_ll.z) || (p->ort Z(i) > pic_ur.z) || 
#endif
          (p->ort Y(i) < pic_ll.y) || (p->ort Y(i) > pic_ur.y)) continue;
      /* which bin? */
      numx = atoms_dist_scale.x * (p->ort X(i) - pic_ll.x);
      if (numx < 0)                   numx = 0;
      if (numx >= atoms_dist_dim.x)   numx = atoms_dist_dim.x-1;
      numy = atoms_dist_scale.y * (p->ort Y(i) - pic_ll.y);
      if (numy < 0)                   numy = 0;
      if (numy >= atoms_dist_dim.y)   numy = atoms_dist_dim.y-1;
      num = numx * atoms_dist_dim.y + numy;
#ifndef TWOD
      numz = atoms_dist_scale.z * (p->ort Z(i) - pic_ll.z);
      if (numz < 0)                   numz = 0;
      if (numz >= atoms_dist_dim.z)   numz = atoms_dist_dim.z-1;
      num = num  * atoms_dist_dim.z + numz;
#endif
      num = SORTE(p,i) * atoms_dist_size + num;
      atoms_dist[num] += 1.0;
    }
  }
}

/******************************************************************************
*
*  write atoms distribution array
*
******************************************************************************/
  
void write_atoms_dist()
{
  int  num, numx, numy, numz, size;
  int  i, k;
  char c;
  str255 fname;
  FILE *out;

  if (myid == 0) {
    sprintf(fname,"%s.atoms_dist",outfilename);
    out = fopen(fname,"w");
    if (NULL == out) error("Cannot open atoms distribution file.");

    c = is_big_endian ? 'B' : 'L';
    fprintf(out,"#F %s %d 0 %d\n", c, DIM, ntypes);
    fprintf(out,"#C");
    for (i=0; i<ntypes; i++) fprintf(out," density_%d",i);
    fprintf(out,"\n");
#ifdef TWOD
    fprintf(out,"#D %d %d\n", atoms_dist_dim.x, atoms_dist_dim.y);
    fprintf(out,"#S %f %f\n", 1.0/atoms_dist_scale.x, 1-0/atoms_dist_scale.y);
#else
    fprintf(out,"#D %d %d %d\n",
      atoms_dist_dim.x, atoms_dist_dim.y, atoms_dist_dim.z);
    fprintf(out,"#S %f %f %f\n", 
      1.0/atoms_dist_scale.x, 1.0/atoms_dist_scale.y, 1.0/atoms_dist_scale.z);
#endif
    fprintf(out,"#E\n");

    size = atoms_dist_size * ntypes;
    if (size!=fwrite(atoms_dist,sizeof(float),size,out))
      error("Cannot write atoms_dist");

    fclose(out);
  }
}

#endif /* ATDIST */
