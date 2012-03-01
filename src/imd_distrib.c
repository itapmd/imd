
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2012 Institute for Theoretical and Applied Physics,
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

float *dat_1 = NULL, *dat_2 = NULL;
float *num_1 = NULL, *num_2 = NULL;
int   dist_size;

/******************************************************************************
*
*  write distributions
*
******************************************************************************/

void write_distrib(int steps)
{
  char contents[255];
  int fzhlr, n, i, j, k;

  is_big_endian = endian();

  /* backup if dist_ur is not set */
  if (0.0==dist_ur.x) {
    dist_ur.x = box_x.x;
    dist_ur.y = box_y.y;
#ifndef TWOD
    dist_ur.z = box_z.z;
#endif
  }

  dist_size  = dist_dim.x * dist_dim.y;
#ifndef TWOD
  dist_size *= dist_dim.z;
#endif
#ifdef BG
  n = 1; /* here we write presstens components in separate files */
#else
  n = dist_presstens_flag ? DIM*(DIM+1)/2 : 1; 
#endif
#ifndef MPI
  dist_chunk_size = dist_size;
#endif
#if defined(BGL) && (defined(TIMING) || defined(DEBUG))
  if (myid==0) 
    printf("%d MB free before distribution allocation\n", get_free_mem());
#endif
#if defined(BG) && defined(NBLIST)
  deallocate_nblist();
#endif
#if defined(BGL) && defined(NBLIST) && (defined(TIMING) || defined(DEBUG))
  if (myid==0) 
    printf("%d MB free after nblist deallocation\n", get_free_mem());
#endif

  /* allocate distribution arrays */
#ifdef MPI2
  MPI_Alloc_mem( n *       dist_size * sizeof(float), MPI_INFO_NULL, &dat_1 );
  MPI_Alloc_mem(           dist_size * sizeof(float), MPI_INFO_NULL, &num_1 );
  MPI_Alloc_mem( n * dist_chunk_size * sizeof(float), MPI_INFO_NULL, &dat_2 );
  num_2 = dat_1;
#elif defined(MPI)
  dat_1 = (float *) malloc( n *       dist_size * sizeof(float) );
  num_1 = (float *) malloc(           dist_size * sizeof(float) );
  dat_2 = (float *) malloc( n * dist_chunk_size * sizeof(float) );
  num_2 = dat_1;
#else
  dat_1 = (float *) malloc( n * dist_size * sizeof(float) );
  num_1 = (float *) malloc(     dist_size * sizeof(float) );
  dat_2 = dat_1;
  num_2 = num_1;
#endif
  if ((NULL==dat_1) || (NULL==num_1) || (NULL==dat_2) || (NULL==num_2))
    error("Cannot allocate distribution data.");

#if defined(BGL) && (defined(TIMING) || defined(DEBUG))
  if (myid==0) 
    printf("%d MB free after distribution allocation\n", get_free_mem());
#endif

  fzhlr = steps / dist_int;

  /* make number density distribution */
  make_distrib_density();

  if (dist_Ekin_flag) {
    make_write_distrib_select( 1, dist_Ekin_fun, 
      dist_Ekin_flag, fzhlr, "Ekin", "Ekin");
  }
  if (dist_Epot_flag) {
    make_write_distrib_select( 1, dist_Epot_fun, 
      dist_Epot_flag, fzhlr, "Epot", "Epot");
  }

#ifdef STRESS_TENS
  if (dist_press_flag) {
    make_write_distrib_select( 1, dist_press_fun, 
      dist_press_flag, fzhlr, "press", "press");
  }
  if (dist_presstens_flag) {
#ifdef TWOD
    sprintf(contents, "P_xx P_yy P_xy");
#else
    sprintf(contents, "P_xx P_yy P_zz P_yz P_zx P_xy");
#endif
#ifdef BG
    /* to save memory, we write each componend in separate file */
    make_write_distrib_select( 1, dist_presstens_xx_fun, 
      dist_presstens_flag, fzhlr, "presstens_xx", "presstens_xx" );
    make_write_distrib_select( 1, dist_presstens_yy_fun,
      dist_presstens_flag, fzhlr, "presstens_yy", "presstens_yy" );
#ifndef TWOD
    make_write_distrib_select( 1, dist_presstens_zz_fun,
      dist_presstens_flag, fzhlr, "presstens_zz", "presstens_zz" );
    make_write_distrib_select( 1, dist_presstens_yz_fun,
      dist_presstens_flag, fzhlr, "presstens_yz", "presstens_yz" );
    make_write_distrib_select( 1, dist_presstens_zx_fun,
      dist_presstens_flag, fzhlr, "presstens_zx", "presstens_zx" );
#endif
    make_write_distrib_select( 1, dist_presstens_xy_fun,
      dist_presstens_flag, fzhlr, "presstens_xy", "presstens_xy" );
#else
    make_write_distrib_select( DIM*(DIM+1)/2, dist_presstens_fun,
      dist_presstens_flag, fzhlr, "presstens", contents);
#endif /* BG */
  }
#endif /* STRESS_TENS */

#ifdef SHOCK
  if (dist_vxavg_flag) {
    make_write_distrib_select( 1, dist_vxavg_fun,
      dist_vxavg_flag, fzhlr, "vxavg", "vxavg");
  }
  if (dist_Ekin_long_flag) {
    make_write_distrib_select( 1, dist_Ekin_long_fun,
      dist_Ekin_long_flag, fzhlr, "Ekin_long", "Ekin_long");
  }
  if (dist_Ekin_trans_flag) {
    make_write_distrib_select( 1, dist_Ekin_trans_fun,
      dist_Ekin_trans_flag, fzhlr, "Ekin_trans", "Ekin_trans");
  }
  if (dist_Ekin_comp_flag) {
    make_write_distrib_select( 1, dist_Ekin_comp_fun,
      dist_Ekin_comp_flag, fzhlr, "Ekin_comp", "Ekin_comp");
  }
  if (dist_shock_shear_flag) {
    make_write_distrib_select( 1, dist_shock_shear_fun,
      dist_shock_shear_flag, fzhlr, "shock_shear", "shock_shear");
  }
  if (dist_shear_aniso_flag) {
    make_write_distrib_select( 1, dist_shear_aniso_fun,
      dist_shear_aniso_flag, fzhlr, "shear_aniso", "shear_aniso");
  }
  if (dist_pressoff_flag) {
    make_write_distrib_select( 1, dist_pressxy_fun,
      dist_pressoff_flag, fzhlr, "pressxy", "pressxy");
    make_write_distrib_select( 1, dist_pressyz_fun,
      dist_pressoff_flag, fzhlr, "pressyz", "pressyz");
    make_write_distrib_select( 1, dist_presszx_fun,
      dist_pressoff_flag, fzhlr, "presszx", "presszx");
  }
#endif /* SHOCK */

  /* write density distribution */
  if ((myid==0) && (dist_dens_flag)) {
    write_distrib_density(dist_dens_flag, fzhlr);
  }

  /* free distribution arrays */
#ifdef MPI2
  MPI_Free_mem(dat_1);
  MPI_Free_mem(num_1);
  MPI_Free_mem(dat_2);
#elif defined(MPI)
  free(dat_1);
  free(num_1);
  free(dat_2);
#else
  free(dat_1);
  free(num_1);
#endif

#if defined(BGL) && (defined(TIMING) || defined(DEBUG))
  if (myid==0) 
    printf("%d MB free after distribution deallocation\n", get_free_mem());
#endif

}

/******************************************************************************
*
*  selection function for kinetic energy
*
******************************************************************************/

void dist_Ekin_fun(float *dat, cell *p, int i)
{
  *dat += SPRODN(IMPULS,p,i,IMPULS,p,i) / (2 * MASSE(p,i));
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
*  selection functions for pressure tensor
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

void dist_presstens_xx_fun(float *dat, cell *p, int i)
{
  *dat += PRESSTENS(p,i,xx);
}

void dist_presstens_yy_fun(float *dat, cell *p, int i)
{
  *dat += PRESSTENS(p,i,yy);
}

#ifndef TWOD

void dist_presstens_zz_fun(float *dat, cell *p, int i)
{
  *dat += PRESSTENS(p,i,zz);
}

void dist_presstens_yz_fun(float *dat, cell *p, int i)
{
  *dat += PRESSTENS(p,i,yz);
}

void dist_presstens_zx_fun(float *dat, cell *p, int i)
{
  *dat += PRESSTENS(p,i,zx);
}

#endif

void dist_presstens_xy_fun(float *dat, cell *p, int i)
{
  *dat += PRESSTENS(p,i,xy);
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
void dist_shear_aniso_fun(float *dat, cell *p, int i)
{
  *dat += PRESSTENS(p,i,yy)-PRESSTENS(p,i,zz);
}
void dist_pressxy_fun(float *dat, cell *p, int i)
{
  *dat += PRESSTENS(p,i,xy);
}
void dist_pressyz_fun(float *dat, cell *p, int i)
{
  *dat += PRESSTENS(p,i,yz);
}
void dist_presszx_fun(float *dat, cell *p, int i)
{
  *dat += PRESSTENS(p,i,zx);
}
#endif

/* transversal kinetic energy */
void dist_Ekin_trans_fun(float *dat, cell *p, int i)
{
  *dat += (SQR(IMPULS(p,i,Y)) + SQR(IMPULS(p,i,Z))) / (4 * MASSE(p,i));
}

/* difference kinetic energy */
void dist_Ekin_comp_fun(float *dat, cell *p, int i)
{
  *dat += (SQR(IMPULS(p,i,Y)) - SQR(IMPULS(p,i,Z))) / (2 * MASSE(p,i));
}

/* longitudinal kinetic energy */
void dist_Ekin_long_fun(float *dat, cell *p, int i)
{
  *dat += SQR(IMPULS(p,i,X) - PXAVG(p,i)) / (2*MASSE(p,i));
}

/* average sample velocity */
void dist_vxavg_fun(float *dat, cell *p, int i)
{
  *dat += IMPULS(p,i,X) / MASSE(p,i);
}

#endif /* SHOCK */

/******************************************************************************
*
*  make density distribution
*
******************************************************************************/

void make_distrib_density(void)
{
  cell *p;
  real scalex, scaley, scalez;
  int  num, numx, numy, numz;
  int  i, j, k, m, chunk_size;

  /* the bins are orthogonal boxes in space */
  scalex = dist_dim.x / (dist_ur.x - dist_ll.x);
  scaley = dist_dim.y / (dist_ur.y - dist_ll.y);
#ifndef TWOD
  scalez = dist_dim.z / (dist_ur.z - dist_ll.z);
#endif

  /* clear distribution */
  for (i=0; i<dist_size; i++) num_2[i] = 0.0;

  /* loop over all atoms */
  for (k=0; k<NCELLS; ++k) {
    p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      /* which bin? */
      numx = scalex * (ORT(p,i,X) - dist_ll.x);
      if ((numx < 0) || (numx >= dist_dim.x)) continue;
      numy = scaley * (ORT(p,i,Y) - dist_ll.y);
      if ((numy < 0) || (numy >= dist_dim.y)) continue;
      num = numx * dist_dim.y + numy;
#ifndef TWOD
      numz = scalez * (ORT(p,i,Z) - dist_ll.z);
      if ((numz < 0) || (numz >= dist_dim.z)) continue;
      num = num * dist_dim.z + numz;
#endif
      num_2[num] += 1.0;
    }
  }

  /* add up results form different CPUs */
#ifdef MPI
#ifdef BG
  /* doing it in several chunks saves memory */
  for (m = 0; m < dist_size; m += dist_chunk_size) { 
    chunk_size = MIN( dist_size - m, dist_chunk_size );
    MPI_Reduce(num_2+m, num_1+m, chunk_size, MPI_FLOAT, MPI_SUM, 0, cpugrid);
  }
#else
  MPI_Reduce( num_2, num_1, dist_size, MPI_FLOAT, MPI_SUM, 0, cpugrid);
#endif
#endif

}

/******************************************************************************
*
*  write density distribution
*
******************************************************************************/

void write_distrib_density(int mode, int fzhlr)
{
  FILE  *outfile;
  char  fname[255];
  float fac, max=0.0, min=1e10;
  int   i, j, count, r, s, t;
  real  vol;

  /* open distribution file, write header */
  sprintf(fname, "%s.%u.%s", outfilename, fzhlr, "dens");
  outfile = fopen(fname, "w");
  if (NULL == outfile) error("Cannot open distribution file.");
  write_distrib_header(outfile, mode, 1, "dens");

  /* compute density, minima and maxima */
  vol = (dist_ur.x - dist_ll.x) * (dist_ur.y - dist_ll.y);
#ifndef TWOD
  vol *= (dist_ur.z - dist_ll.z);
#endif
  fac = dist_size / vol;
  for (i=0; i<dist_size; i++) { 
    dat_1[i] = num_1[i] * fac;
    max = MAX( max, dat_1[i] );
    min = MIN( min, dat_1[i] );
  }

  /* write distribution */
  if (mode==DIST_FORMAT_BINARY) {
    count = fwrite(dat_1, sizeof(float), dist_size, outfile); 
    if (count != dist_size) warning("distribution write incomplete!");
  } 
  else if ((mode==DIST_FORMAT_ASCII) || (mode==DIST_FORMAT_ASCII_COORD)) {
    i=0;
    for (r=0; r<dist_dim.x; r++)
      for (s=0; s<dist_dim.y; s++)
#ifndef TWOD
        for (t=0; t<dist_dim.z; t++)
#endif
	{
	  if (mode==DIST_FORMAT_ASCII_COORD) {
#ifdef TWOD
            fprintf(outfile, "%d %d ", r, s);
#else
            fprintf(outfile, "%d %d %d ", r, s, t);
#endif
          }
          fprintf(outfile,"%e\n", dat_1[i++]);
        }
  }
  else error("unknown distribution output format");
  fclose(outfile);

  /* write minmax */
  sprintf(fname, "%s.minmax.%s", outfilename, "dens");
  outfile = fopen(fname, "a");
  if (NULL == outfile) error("Cannot open minmax file.");
  fprintf(outfile, "%d %e %e\n", fzhlr, min, max);
  fclose(outfile);

}

/******************************************************************************
*
*  make and write distribution of selected variables
*
******************************************************************************/

void make_write_distrib_select(int n, void (*fun)(float*, cell*, int),
                               int mode, int fzhlr, char *suffix, char *cont)
{
  cell  *p;
  real  scalex, scaley, scalez;
  int   num, numx, numy, numz;
  int   i, j, k, count, r, s, t, m, chunk_size;
  float max[6], min[6];
  FILE  *outfile=NULL;
  char  fname[255];

  /* the bins are orthogonal boxes in space */
  scalex = dist_dim.x / (dist_ur.x - dist_ll.x);
  scaley = dist_dim.y / (dist_ur.y - dist_ll.y);
#ifndef TWOD
  scalez = dist_dim.z / (dist_ur.z - dist_ll.z);
#endif

  /* clear distribution */
  for (i=0; i<n*dist_size; i++) dat_1[i] = 0.0;

  /* loop over all atoms */
  for (k=0; k<NCELLS; ++k) {
    p = CELLPTR(k);
    for (i=0; i<p->n; ++i) {
      /* which bin? */
      numx = scalex * (ORT(p,i,X) - dist_ll.x);
      if ((numx < 0) || (numx >= dist_dim.x)) continue;
      numy = scaley * (ORT(p,i,Y) - dist_ll.y);
      if ((numy < 0) || (numy >= dist_dim.y)) continue;
      num = numx * dist_dim.y + numy;
#ifndef TWOD
      numz = scalez * (ORT(p,i,Z) - dist_ll.z);
      if ((numz < 0) || (numz >= dist_dim.z)) continue;
      num = num * dist_dim.z + numz;
#endif
      (*fun)(dat_1 + n * num, p, i);
    }
  }

  /* open distribution file, write header */
  if (myid==0) {
    sprintf(fname, "%s.%u.%s", outfilename, fzhlr, suffix);
    outfile = fopen(fname, "w");
    if (NULL == outfile) error("Cannot open distribution file.");
    write_distrib_header(outfile, mode, n, cont);
    for (k=0; k<n; k++) {
      min[k] =  1e+10;
      max[k] = -1e+10;
    }
  }

  /* collect and write data in handy chunks */
  r = s = t = 0;
  for (m = 0; m < dist_size; m += dist_chunk_size) { 

    chunk_size = MIN( dist_size - m, dist_chunk_size );

    /* add up results form different CPUs */
#ifdef MPI
    MPI_Reduce(dat_1+n*m, dat_2, n*chunk_size, MPI_FLOAT, MPI_SUM, 0, cpugrid);
#endif

    if (myid==0) {

      /* normalize distribution, compute minima and maxima */
      for (i=0; i<chunk_size; i++) {
        if (num_1[m+i] > 0.0) {
          for (k=0; k<n; k++) {
            dat_2[n*i+k] /= num_1[m+i];
            min[k] = MIN( min[k], dat_2[n*i+k] );
            max[k] = MAX( max[k], dat_2[n*i+k] );
          }
        }
      }

      /* write distribution */
      if (mode==DIST_FORMAT_BINARY) {
        count = fwrite(dat_2, sizeof(float), n*chunk_size, outfile); 
        if (count != n*chunk_size) warning("distribution write incomplete!");
      } 
      else if ((mode==DIST_FORMAT_ASCII) || (mode==DIST_FORMAT_ASCII_COORD)) {
	for (i=0; i<chunk_size; i++) {
          if (mode==DIST_FORMAT_ASCII_COORD) {
#ifdef TWOD
            fprintf(outfile, "%d %d", r, s++);
            if (s==dist_dim.y) { s=0; r++; }
#else
            fprintf(outfile, "%d %d %d", r, s, t++);
            if (t==dist_dim.z) { t=0; s++; }
            if (s==dist_dim.y) { s=0; r++; }
#endif
          }
#ifdef BG
          fprintf(outfile," %e\n", dat_2[i]);
#else
          for (k=0; k<n; k++) fprintf(outfile," %e", dat_2[n*i+k]);
          fprintf(outfile, "\n");
#endif
        }
      }
      else error("unknown distribution output format");
    }
  }

  /* close distribution file */
  if (myid==0) fclose(outfile);

  /* write minmax */
  if (myid==0) {
    sprintf(fname, "%s.minmax.%s", outfilename, suffix);
    outfile = fopen(fname, "a");
    if (NULL == outfile) error("Cannot open minmax file.");
    fprintf( outfile, "%d ", fzhlr );
    for (i=0; i<n; i++) fprintf(outfile, " %e %e", min[i], max[i]);
    fprintf(outfile, "\n");
    fclose(outfile);
  }

}

/******************************************************************************
*
*  write header of distribution files 
*
******************************************************************************/

void write_distrib_header(FILE *out, int mode, int n, char *cont)
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
  fprintf(out, "#D %d %d\n",    dist_dim.x, dist_dim.y);
#else
  fprintf(out, "#D %d %d %d\n", dist_dim.x, dist_dim.y, dist_dim.z);
#endif

  /* bin size line */
  s.x = (dist_ur.x - dist_ll.x) / dist_dim.x;
  s.y = (dist_ur.y - dist_ll.y) / dist_dim.y;
#ifdef TWOD
  fprintf(out, "#S %e %e\n",    s.x, s.y);
#else
  s.z = (dist_ur.z - dist_ll.z) / dist_dim.z;
  fprintf(out, "#S %e %e %e\n", s.x, s.y, s.z);
#endif

  /* endheader line */
  time(&now);
  fprintf(out, "## Generated on %s", ctime(&now) ); 
  fprintf(out, "## by %s (version of %s)\n", progname, DATE);
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
#ifdef CM_HACK
  static vektor tot_velocity = {0.0,0.0,0.0};
  static int step_count=0, max_count=100;
  int count=0;
#endif

  co = cos(atdist_phi);
  si = sin(atdist_phi);

#ifdef CM_HACK
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
#endif /* CM_HACK */

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

  is_big_endian = endian();

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
  imd_init_timer( &time_fft,      0, NULL, NULL );
  imd_init_timer( &time_fft_plan, 0, NULL, NULL );
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
