
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2004 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
*  imd_distrib.c -- distributions of various quantities
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/******************************************************************************
*
*  write distributions
*
******************************************************************************/

void write_distrib(int steps)
{
  char contents[255];
  static dist_t dist;
  int flag, fzhlr;

  /* backup if dist_ur is not set */
  if (0.0==dist_ur.x) {
    dist_ur.x = box_x.x;
    dist_ur.y = box_y.y;
#ifndef TWOD
    dist_ur.z = box_z.z;
#endif
  }

  dist.dim = dist_dim;
  dist.ll  = dist_ll;
  dist.ur  = dist_ur;

  flag = 1;
  fzhlr = steps / dist_int;

  if (dist_Ekin_flag) {
    make_distrib_select( &dist, 1, &flag, dist_Ekin_fun);
    if (myid==0) write_distrib_select( &dist, dist_Ekin_flag, 1, fzhlr, 
                                       "Ekin", "Ekin");
  }
  if (dist_Epot_flag) {
    make_distrib_select( &dist, 1, &flag, dist_Epot_fun);
    if (myid==0) write_distrib_select( &dist, dist_Epot_flag, 1, fzhlr, 
                                       "Epot", "Epot");
  }

#ifdef STRESS_TENS
  if (dist_press_flag) {
    make_distrib_select( &dist, 1, &flag, dist_press_fun);
    if (myid==0) write_distrib_select( &dist, dist_press_flag, 1, fzhlr, 
                                       "press", "press");
  }
  if (dist_presstens_flag) {
#ifdef TWOD
    sprintf(contents, "P_xx P_yy P_xy");
#else
    sprintf(contents, "P_xx P_yy P_zz P_yz P_zx P_xy");
#endif
    make_distrib_select( &dist, DIM*(DIM+1)/2, &flag, dist_presstens_fun);
    if (myid==0) write_distrib_select( &dist, dist_presstens_flag, 
                              DIM*(DIM+1)/2, fzhlr, "presstens", contents);
  }
#endif /* STRESS_TENS */

#ifdef SHOCK
  if (dist_Ekin_long_flag) {
    make_distrib_select( &dist, 1, &flag, dist_Ekin_long_fun);
    if (myid==0) write_distrib_select( &dist, dist_Ekin_long_flag, 1, fzhlr, 
                                       "Ekin_long", "Ekin_long");
  }
  if (dist_Ekin_trans_flag) {
    make_distrib_select( &dist, 1, &flag, dist_Ekin_trans_fun);
    if (myid==0) write_distrib_select( &dist, dist_Ekin_trans_flag, 1, fzhlr, 
                                       "Ekin_trans", "Ekin_trans");
  }
  if (dist_shock_shear_flag) {
    make_distrib_select( &dist, 1, &flag, dist_shock_shear_fun);
    if (myid==0) write_distrib_select( &dist, dist_shock_shear_flag, 1, fzhlr, 
                                       "shock_shear", "shock_shear");
  }
#endif /* SHOCK */

}

/******************************************************************************
*
*  selection function for kinetic energy
*
******************************************************************************/

void dist_Ekin_fun(float *dat, cell *p, int i)
{
  *dat += SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / (2 * MASSE(p,i));
}

/******************************************************************************
*
*  selection function for potential energy
*
******************************************************************************/

void dist_Epot_fun(float *dat, cell *p, int i)
{
#ifdef DISLOC
  if (Epot_diff==1)
    *dat += POTENG(p,i) - EPOT_REF(p,i);
  else
#endif
    *dat += POTENG(p,i);
}

#ifdef STRESS_TENS

/******************************************************************************
*
*  selection function for scalar pressure
*
******************************************************************************/

void dist_press_fun(float *dat, cell *p, int i)
{
#ifdef TWOD
  *dat += (PRESSTENS(p,i,xx) + PRESSTENS(p,i,yy)) / 2.0;
#else
  *dat += (PRESSTENS(p,i,xx) + PRESSTENS(p,i,yy) + PRESSTENS(p,i,zz)) / 3.0;
#endif
}

/******************************************************************************
*
*  selection function for pressure tensor
*
******************************************************************************/

void dist_presstens_fun(float *dat, cell *p, int i)
{
  int k=0;
  dat[k++] += PRESSTENS(p,i,xx);
  dat[k++] += PRESSTENS(p,i,yy);
#ifndef TWOD
  dat[k++] += PRESSTENS(p,i,zz);
  dat[k++] += PRESSTENS(p,i,yz);
  dat[k++] += PRESSTENS(p,i,zx);
#endif
  dat[k++] += PRESSTENS(p,i,xy);
}

#endif /* STRESS_TENS */

#ifdef SHOCK

/******************************************************************************
*
*  selection function for various shock quantities
*
******************************************************************************/

#ifdef STRESS_TENS
/* shear stress */
void dist_shock_shear_fun(float *dat, cell *p, int i)
{
  *dat += (PRESSTENS(p,i,xx)-(PRESSTENS(p,i,yy)+PRESSTENS(p,i,zz))/2.0)/2.0;
}
#endif

/* longitudinal kinetic energy */
void dist_Ekin_trans_fun(float *dat, cell *p, int i)
{
  *dat += (SQR(IMPULS(p,i,Y)) + SQR(IMPULS(p,i,Z))) / (4 * MASSE(p,i));
}

/* transversal kinetic energy */
void dist_Ekin_long_fun(float *dat, cell *p, int i)
{
  float tmp;

  /* average v_xx - u_p  relative to moving pistons */
  tmp = shock_speed * MASSE(p,i);
  /* plate against bulk */
  if (shock_mode == 1) {
    if ( ORT(p,i,X) < shock_strip ) 
      *dat += SQR(IMPULS(p,i,X) - tmp) / (2*MASSE(p,i));
    else
      *dat += SQR(IMPULS(p,i,X)) / (2*MASSE(p,i));
  }
  /* two halves against one another */
  if (shock_mode == 2) {
    if ( ORT(p,i,X) < box_x.x*0.5 )
      *dat += SQR(IMPULS(p,i,X) - tmp) / (2*MASSE(p,i));
    else
      *dat += SQR(IMPULS(p,i,X) + tmp) / (2*MASSE(p,i));
  }
  /* bulk against wall */
  if (shock_mode == 3) 
    *dat += SQR(IMPULS(p,i,X) - tmp) / (2*MASSE(p,i));
}

#endif /* SHOCK */


/******************************************************************************
*
*  make distribution of selected variables
*
******************************************************************************/

void make_distrib_select(dist_t *dist, int n, int *flag, 
                         void (*fun)(float*, cell*, int))
{
  cell *p;
  real scalex, scaley, scalez, Ekin, tmp;
  int  num, numx, numy, numz, size;
  int  i, j, k;
  static float   *dat_1 = NULL, *dat_2 = NULL;
  static float   *min   = NULL, *max   = NULL;
  static integer *num_1 = NULL, *num_2 = NULL;

  /* the bins are orthogonal boxes in space */
  scalex = dist->dim.x / (dist->ur.x - dist->ll.x);
  scaley = dist->dim.y / (dist->ur.y - dist->ll.y);
  size   = dist->dim.x * dist->dim.y;
#ifndef TWOD
  scalez = dist->dim.z / (dist->ur.z - dist->ll.z);
  size  *= dist->dim.z;
#endif
  dist->size   = size;

  /* allocate distribution arrays */
  if (NULL==dat_1) {
    dat_1 = (float   *) realloc( dat_1, n * size * sizeof(float  ) );
    num_1 = (integer *) realloc( num_1,     size * sizeof(integer) );
    min   = (float   *) realloc( min,   n *        sizeof(float  ) );
    max   = (float   *) realloc( max,   n *        sizeof(float  ) );
    if ((NULL==dat_1) || (NULL==num_1) || (NULL==min) || (NULL==max))
      error("Cannot allocate distribution data.");
  }
#ifdef MPI
  if (NULL==dat_2) {
    dat_2 = (float   *) realloc( dat_2, n * size * sizeof(float  ) );
    num_2 = (integer *) realloc( num_2,     size * sizeof(integer) );
    if ((NULL==dat_2) || (NULL==num_2))
      error("Cannot allocate distribution data.");
  }
#endif

  /* clear distributions */
  for (i=0; i<n*size; i++) {
    dat_1[i] = 0.0;
  }
  for (i=0; i<size; i++) {
    num_1[i] = 0;
  }

  /* loop over all atoms */
  for (k=0; k<NCELLS; ++k) {
    p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      /* which bin? */
      numx = scalex * (ORT(p,i,X) - dist->ll.x);
      if ((numx < 0) || (numx >= dist->dim.x)) continue;
      numy = scaley * (ORT(p,i,Y) - dist->ll.y);
      if ((numy < 0) || (numy >= dist->dim.y)) continue;
      num = numx * dist->dim.y + numy;
#ifndef TWOD
      numz = scalez * (ORT(p,i,Z) - dist->ll.z);
      if ((numz < 0) || (numz >= dist->dim.z)) continue;
      num = num * dist->dim.z + numz;
#endif
      num_1[num]++;
      num *= n;
      (*fun)(dat_1 + num, p, i);
    }
  }

#ifdef MPI  /* add up results form different CPUs */
  MPI_Reduce( dat_1, dat_2, n * size, MPI_FLOAT, MPI_SUM, 0, cpugrid);
  dist->dat = dat_2;
  if (*flag) { /* communicate the same num only once */
    MPI_Reduce( num_1, num_2, size, INTEGER,   MPI_SUM, 0, cpugrid);
    dist->num = num_2;
    *flag = 0;
  }
#else
  dist->dat = dat_1; 
  dist->num = num_1;
#endif
  dist->min = min;
  dist->max = max;

  /* normalize distributions */
  if (myid==0) {
    for (i=0; i < size; i++) {
      if (dist->num[i] > 0) {
        for (j=0; j<n; j++) dist->dat[n*i+j] /= dist->num[i];
      }
    }
  } 

  /* compute minima and maxima of distributions */
  if (myid==0) {
    j=0;
    while (dist->num[j]==0) j++;
    for (k=0; k<n; k++) {
      min[k] = dist->dat[j];
      max[k] = dist->dat[j];
      for (i=j+1; i<size; i++) {
        if (dist->num[i]>0) {
          dist->min[k]  = MIN( min[k],  dist->dat[n*i+k] );
          dist->max[k]  = MAX( max[k],  dist->dat[n*i+k] );
	}
      }
    }
  }

}

/******************************************************************************
*
*  write distribution of selected variables
*
******************************************************************************/

void write_distrib_select( dist_t *dist, int mode, int n, int fzhlr,
                           char *suffix, char *cont)
{
  FILE *outfile;
  char fname[255];
  int i, j, count, r, s, t;

  /* write minmax */
  sprintf(fname, "%s.minmax.%s", outfilename, suffix);
  outfile = fopen(fname, "a");
  if (NULL == outfile) error("Cannot open minmax file.");
  fprintf( outfile, "%d ", fzhlr );
  for (i=0; i<n; i++) fprintf(outfile, " %e %e", dist->min[i], dist->max[i]);
  fprintf(outfile, "\n");
  fclose(outfile);

  /* open distribution file, write header */
  sprintf(fname, "%s.%u.%s", outfilename, fzhlr, suffix);
  outfile = fopen(fname, "w");
  if (NULL == outfile) error("Cannot open distribution file.");
  write_distrib_header(outfile, dist, mode, n, cont);

  /* write distribution */
  if (mode==DIST_FORMAT_BINARY) {
    count = fwrite(dist->dat, sizeof(float), n*dist->size, outfile); 
    if (count != n*dist->size) warning("distribution write incomplete!");
  } 
  else if ((mode==DIST_FORMAT_ASCII) || (mode==DIST_FORMAT_ASCII_COORD)) {
    i=0;
    for (r=0; r<dist->dim.x; r++)
      for (s=0; s<dist->dim.y; s++)
#ifndef TWOD
        for (t=0; t<dist->dim.z; t++)
#endif
	{
	  if (mode==DIST_FORMAT_ASCII_COORD) {
#ifdef TWOD
            fprintf(outfile, "%d %d", r, s);
#else
            fprintf(outfile, "%d %d %d", r, s, t);
#endif
          }
          for (j=0; j<n; j++) fprintf(outfile," %e", dist->dat[i++]);
          fprintf(outfile, "\n");
        }
  }
  else error("unknown distribution output format");
  fclose(outfile);
}

/******************************************************************************
*
*  write header of distribution files 
*
******************************************************************************/

void write_distrib_header(FILE *out, dist_t *dist, int mode, int n, char *cont)
{
  char c;
  int  n_coord;
  time_t now;
  vektor s;

  /* format line -- format dim n_coord n_data */
  if (mode==DIST_FORMAT_BINARY)
    c = is_big_endian ? 'B' : 'L';
  else
    c = 'A';
  n_coord = (mode==DIST_FORMAT_ASCII_COORD) ? DIM : 0;
  fprintf(out, "#F %c %d %d %d\n", c, DIM, n_coord, n);

  /* contents line */
  if (mode==DIST_FORMAT_ASCII_COORD)
#ifdef TWO
    fprintf(out, "#C x y %s\n",   cont);
#else
    fprintf(out, "#C x y z %s\n", cont);
#endif
  else
    fprintf(out, "#C %s\n",       cont);

  /* dimension line */
#ifdef TWOD
  fprintf(out, "#D %d %d\n",    dist->dim.x, dist->dim.y);
#else
  fprintf(out, "#D %d %d %d\n", dist->dim.x, dist->dim.y, dist->dim.z);
#endif

  /* bin size line */
  s.x = (dist->ur.x - dist->ll.x) / dist->dim.x;
  s.y = (dist->ur.y - dist->ll.y) / dist->dim.y;
#ifdef TWOD
  fprintf(out, "#S %e %e\n",    s.x, s.y);
#else
  s.z = (dist->ur.z - dist->ll.z) / dist->dim.z;
  fprintf(out, "#S %e %e %e\n", s.x, s.y, s.z);
#endif

  /* endheader line */
  time(&now);
  fprintf(out, "## Generated by %s on %s", progname, ctime(&now));
  fprintf(out, "#E\n");

}

#ifdef ATDIST

/******************************************************************************
*
*  initialize atoms distribution array
*
******************************************************************************/
  
void init_atdist()
{
  int size, i;

  /* compute array size */
  atdist_size  = atdist_dim.x * atdist_dim.y;
#ifndef TWOD
  atdist_size *= atdist_dim.z;
#endif
  size = atdist_size * ntypes;

  /* atdist_ll and atdist_ur must be set */
  if (0.0==atdist_ur.x) {
    error("atdist_ll and atdist_ur must be set");
  }

  /* the bins are orthogonal boxes in space */
  atdist_scale.x = atdist_dim.x / (atdist_ur.x - atdist_ll.x);
  atdist_scale.y = atdist_dim.y / (atdist_ur.y - atdist_ll.y);
#ifndef TWOD
  atdist_scale.z = atdist_dim.z / (atdist_ur.z - atdist_ll.z);
#endif

  /* allocate distribution array */
  if (NULL==atdist) {
    atdist = (float *) malloc( size * sizeof(float) );
    if (NULL==atdist) error("Cannot allocate atoms distribution array.");
  }

  /* initialize distribution array */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i=0; i<size; i++) atdist[i]=0.0;
}

/******************************************************************************
*
*  update atoms distribution array
*
******************************************************************************/
  
void update_atdist()
{
  int  num, numx, numy, numz;
  cell *p;
  int  i, k, ix, iy, iz;
  real x, y, z, t, co, si;
  static vektor tot_velocity = {0.0,0.0,0.0};
  static int step_count=0, max_count=100;
  int count=0;

  co = cos(atdist_phi);
  si = sin(atdist_phi);

  /* correct center of mass velocity */
  for (k=0; k<NCELLS; ++k) {
    int i;
    cell *p;
    p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      if (SORTE(p,i)>0) {
        tot_velocity.x += IMPULS(p,i,X) / MASSE(p,i);
        tot_velocity.y += IMPULS(p,i,Y) / MASSE(p,i);
        tot_velocity.z += IMPULS(p,i,Z) / MASSE(p,i);
        count++;
      }
    }
  }
  step_count++;
  if (step_count % max_count) {
    for (k=0; k<NCELLS; ++k) {
      int i;
      cell *p;
      p = CELLPTR(k);
      for (i=0; i<p->n; ++i) {
        IMPULS(p,i,X) -= tot_velocity.x * MASSE(p,i) / count;
        IMPULS(p,i,Y) -= tot_velocity.y * MASSE(p,i) / count;
        IMPULS(p,i,Z) -= tot_velocity.z * MASSE(p,i) / count;
      }
    }
    tot_velocity.x = 0.0;
    tot_velocity.y = 0.0;
    tot_velocity.z = 0.0;
    count=0;
  } else {
    step_count++;
  }

  /* loop over all atoms */
  for (k=0; k<NCELLS; ++k) {
    p = CELLPTR(k);
    for (i=0; i<p->n; ++i) 

      /* periodic continuation */
      for (ix=atdist_per_ll.x; ix<=atdist_per_ur.x; ix++)
#ifndef TWOD
        for (iz=atdist_per_ll.z; iz<=atdist_per_ur.z; iz++)
#endif
          for (iy=atdist_per_ll.y; iy<=atdist_per_ur.y; iy++) {
#ifdef TWOD
            x = ORT(p,i,X) + ix * box_x.x + iy * box_y.x;
            y = ORT(p,i,Y) + ix * box_x.y + iy * box_y.y;
#else
            x = ORT(p,i,X) + ix * box_x.x + iy * box_y.x + iz * box_z.x;
            y = ORT(p,i,Y) + ix * box_x.y + iy * box_y.y + iz * box_z.y;
            z = ORT(p,i,Z) + ix * box_x.z + iy * box_y.z + iz * box_z.z;
#endif
            t =  co * x + si * y;
            y = -si * x + co * y;
            x = t;

            /* continue if atom is not inside selected box */
            if ((x < atdist_ll.x) || (x > atdist_ur.x) ||
#ifndef TWOD
                (z < atdist_ll.z) || (z > atdist_ur.z) || 
#endif
                (y < atdist_ll.y) || (y > atdist_ur.y)) continue;

            /* which bin? */
            numx = atdist_scale.x * (x - atdist_ll.x);
            if (numx < 0)                   numx = 0;
            if (numx >= atdist_dim.x)   numx = atdist_dim.x-1;
            numy = atdist_scale.y * (y - atdist_ll.y);
            if (numy < 0)                   numy = 0;
            if (numy >= atdist_dim.y)   numy = atdist_dim.y-1;
            num = numx * atdist_dim.y + numy;
#ifndef TWOD
            numz = atdist_scale.z * (z - atdist_ll.z);
            if (numz < 0)                   numz = 0;
            if (numz >= atdist_dim.z)   numz = atdist_dim.z-1;
            num = num  * atdist_dim.z + numz;
#endif
            num = SORTE(p,i) * atdist_size + num;
            atdist[num] += 1.0;
	  }
  }
}

/******************************************************************************
*
*  write atoms distribution array
*
******************************************************************************/
  
void write_atdist()
{
  int  num, numx, numy, numz, size;
  int  i, k;
  char c;
  str255 fname;
  FILE *out;

  if (myid == 0) {
    sprintf(fname,"%s.atdist",outfilename);
    out = fopen(fname,"w");
    if (NULL == out) error("Cannot open atoms distribution file.");

    c = is_big_endian ? 'B' : 'L';
    fprintf(out,"#F %c %d 0 %d\n", c, DIM, ntypes);
    fprintf(out,"#C");
    for (i=0; i<ntypes; i++) fprintf(out," density_%d",i);
    fprintf(out,"\n");
#ifdef TWOD
    fprintf(out,"#D %d %d\n", atdist_dim.x, atdist_dim.y);
    fprintf(out,"#S %f %f\n", 1.0 / atdist_scale.x, 1.0 / atdist_scale.y);
#else
    fprintf(out,"#D %d %d %d\n",
      atdist_dim.x, atdist_dim.y, atdist_dim.z);
    fprintf(out,"#S %f %f %f\n", 
      1.0 / atdist_scale.x, 1.0 / atdist_scale.y, 1.0 / atdist_scale.z);
#endif
    fprintf(out,"#E\n");

    size = atdist_size * ntypes;
    if (size!=fwrite(atdist, sizeof(float), size, out))
      error("Cannot write atoms distribution");

    fclose(out);
  }
}

#endif /* ATDIST */

#ifdef DIFFPAT

/******************************************************************************
*
*  initialize atoms distribution array
*
******************************************************************************/
  
void init_diffpat()
{
  int size, i;

#ifdef OMP
  fftwf_init_threads();
  fftwf_plan_with_nthreads(omp_get_max_threads());
#endif

#ifdef TIMING
  time_fft.total      = 0.0;
  time_fft_plan.total = 0.0;
#endif

  /* compute array size */
  diffpat_size  = diffpat_dim.x * diffpat_dim.y * 2 * (diffpat_dim.z / 2 + 1);

  /* diffpat_ll and diffpat_ur must be set */
  if (0.0==diffpat_ur.x) {
    error("diffpat_ll and diffpat_ur must be set");
  }

  /* diffpat_weight must be set */
  if (0.0==diffpat_weight[0]) {
    error("diffpat_weight must be set");
  }

  /* the bins are orthogonal boxes in space */
  diffpat_scale.x = diffpat_dim.x / (diffpat_ur.x - diffpat_ll.x);
  diffpat_scale.y = diffpat_dim.y / (diffpat_ur.y - diffpat_ll.y);
  diffpat_scale.z = diffpat_dim.z / (diffpat_ur.z - diffpat_ll.z);

  /* allocate arrays */
  if (NULL==diffdist) {
    diffdist = (float *) fftwf_malloc( diffpat_size * sizeof(float) );
    diffpat  = (float *) malloc( (diffpat_size / 2) * sizeof(float) );
    if ((NULL==diffdist) || (NULL==diffpat))
      error("Cannot allocate diffraction pattern array.");
  }

  /* make fftw plan */
#ifdef TIMING
  imd_start_timer(&time_fft_plan);
#endif
  if ((diffpat_end - diffpat_start) % diffpat_int > 50)
    diffpat_plan = fftwf_plan_dft_r2c_3d(
      diffpat_dim.x, diffpat_dim.y, diffpat_dim.z,
      diffdist, (fftwf_complex *) diffdist, FFTW_MEASURE);
  else
    diffpat_plan = fftwf_plan_dft_r2c_3d(
      diffpat_dim.x, diffpat_dim.y, diffpat_dim.z,
      diffdist, (fftwf_complex *) diffdist, FFTW_ESTIMATE);
#ifdef TIMING
  imd_stop_timer(&time_fft_plan);
  printf("Time for FFT plan: %f\n", time_fft_plan.total);
#endif

  /* initialize arrays */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i=0; i<diffpat_size;   i++) diffdist[i]=0.0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i=0; i<diffpat_size/2; i++) diffpat [i]=0.0;

}

/******************************************************************************
*
*  update atoms distribution array
*
******************************************************************************/
  
void update_diffpat(int steps)
{
  int   num, numx, numy, numz, k, i;
  real  x, y, z;
  int   dimz  = 2 * (diffpat_dim.z / 2 + 1);
  int   dimz2 =     (diffpat_dim.z / 2 + 1);
  fftwf_complex *dist_out = (fftwf_complex *) diffdist;

  /* loop over all atoms */
  for (k=0; k<NCELLS; ++k) {
    cell  *p;
    p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      x = ORT(p,i,X);
      y = ORT(p,i,Y);
      z = ORT(p,i,Z);
      /* continue if atom is not inside selected box */
      if ((x < diffpat_ll.x) || (x > diffpat_ur.x) ||
          (z < diffpat_ll.z) || (z > diffpat_ur.z) || 
          (y < diffpat_ll.y) || (y > diffpat_ur.y)) continue;
      /* which bin? */
      numx = diffpat_scale.x * (x - diffpat_ll.x);
      if (numx < 0)              numx = 0;
      if (numx >= diffpat_dim.x) numx = diffpat_dim.x-1;
      numy = diffpat_scale.y * (y - diffpat_ll.y);
      if (numy < 0)              numy = 0;
      if (numy >= diffpat_dim.y) numy = diffpat_dim.y-1;
      numz = diffpat_scale.z * (z - diffpat_ll.z);
      if (numz < 0)              numz = 0;
      if (numz >= diffpat_dim.z) numz = diffpat_dim.z-1;
      num = (numx * diffpat_dim.y + numy) * dimz + numz;
      diffdist[num] += diffpat_weight[ SORTE(p,i) ];
    }
  }

  /* increment diffraction pattern */
  if (0==steps%diffpat_int) {
#ifdef TIMING
    imd_start_timer(&time_fft);
#endif
    fftwf_execute(diffpat_plan);
#ifdef TIMING
    imd_stop_timer(&time_fft);
#endif
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i=0; i<diffpat_size/2; i++)
      diffpat[i] += (float)(SQR(dist_out[i][0])+SQR(dist_out[i][1]));
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i=0; i<diffpat_size; i++) diffdist[i]=0.0;
  }

}

/******************************************************************************
*
*  write diffraction pattern
*
******************************************************************************/

void write_diffpat()
{
  int  num, numx, numy, numz, len;
  int  dimz2 = diffpat_dim.z / 2 + 1;
  real pi, ddx, ddy, ddz;
  char c;
  str255 fname;
  FILE *out;

  if (myid == 0) {

    /* open file */
    sprintf(fname,"%s.diffpat",outfilename);
    if (NULL==(out=fopen(fname,"w"))) error("Cannot open output file");

    pi  = 4 * atan( (double) 1.0 );
    ddx = 2 * pi * diffpat_scale.x / diffpat_dim.x;
    ddy = 2 * pi * diffpat_scale.y / diffpat_dim.y;
    ddz = 2 * pi * diffpat_scale.z / diffpat_dim.z;

    /* write file header */
    if (endian()) c='B'; else c='L';
    fprintf(out, "#F %c 3 0 1\n", c );
    fprintf(out, "#C Fourier\n");
    fprintf(out, "#D %d %d %d\n", diffpat_dim.x, diffpat_dim.y, dimz2);
    fprintf(out, "#S %e %e %e\n", ddx, ddy, ddz);
    fprintf(out, "#E\n");

    /* write data */
    len = diffpat_size / 2;
    if (len!=fwrite(diffpat, sizeof(float), len, out))
      error("Cannot write distribution");
    fclose(out);

#ifdef TIMING
    printf("Time for FFT: %f\n", time_fft.total);
#endif

  }
}

#endif /* DIFFPAT */
