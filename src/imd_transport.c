
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2011 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_transport.c -- heat transport
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

#ifdef HC

/******************************************************************************
*
*  add up microscopic heat current
*
******************************************************************************/

void do_heat_cond(int steps)
{
  int    i, k;
  real   e, tmp;
  vektor pp, vv, tmpv;
  static real fac = 0.0;

#ifdef TWOD
  hc.x = 0.0; hc.y = 0.0;
#else
  hc.x = 0.0; hc.y = 0.0; hc.z = 0.0;
#endif

  for (k=0; k<NCELLS; ++k) { /* loop over all cells */

    cell *p = CELLPTR(k);

    for (i=0; i<p->n; ++i) { /* loop over all atoms in the cell */

      /* momenta at the time of the force computation */
      pp.x = IMPULS(p,i,X) + 0.5 * timestep * KRAFT(p,i,X); 
      pp.y = IMPULS(p,i,Y) + 0.5 * timestep * KRAFT(p,i,Y); 
#ifndef TWOD
      pp.z = IMPULS(p,i,Z) + 0.5 * timestep * KRAFT(p,i,Z); 
#endif
      /* current total energy */ 
      e = SPROD(pp,pp) / (2*MASSE(p,i));
      if (steps < hc_start) fac += e; /* add up kin. Energy for av. temp. */
      e += POTENG(p,i);

      /* average total energy of each atom */
      if   (steps == hc_av_start) HCAVENG(p,i)  = e;
      else if (steps <  hc_start) HCAVENG(p,i) += e;
      else if (steps == hc_start) HCAVENG(p,i) /= (hc_start - hc_av_start);

      /* add up microscopic heat conductivity */ 
      /* at this point, PRESSTENS consists only of the virial part */
      if (steps >= hc_start) {
#ifdef TWOD
        vv.x = PRESSTENS(p,i,xx) * pp.x + PRESSTENS(p,i,xy) * pp.y;
        vv.y = PRESSTENS(p,i,xy) * pp.x + PRESSTENS(p,i,yy) * pp.y;
#else
        vv.x = PRESSTENS(p,i,xx) * pp.x + PRESSTENS(p,i,xy) * pp.y
                                        + PRESSTENS(p,i,zx) * pp.z;
        vv.y = PRESSTENS(p,i,xy) * pp.x + PRESSTENS(p,i,yy) * pp.y
                                        + PRESSTENS(p,i,yz) * pp.z;
        vv.z = PRESSTENS(p,i,zx) * pp.x + PRESSTENS(p,i,yz) * pp.y
                                        + PRESSTENS(p,i,zz) * pp.z;
#endif
        e -= HCAVENG(p,i);  /* subtract average energy */
        hc.x += (pp.x * e + 0.5 * vv.x) / MASSE(p,i);
        hc.y += (pp.y * e + 0.5 * vv.y) / MASSE(p,i);
#ifndef TWOD
        hc.z += (pp.z * e + 0.5 * vv.z) / MASSE(p,i);
#endif
      }
    }
  }
#ifdef MPI
  if (steps == hc_start) {
    MPI_Allreduce( &fac, &tmp, 1, REAL, MPI_SUM, cpugrid);
    fac = tmp;
  }
  MPI_Allreduce( &hc, &tmpv, DIM, REAL, MPI_SUM, cpugrid);
  hc = tmpv;
#endif

  /* scaling factor: 1 / (sqrt(volume) * temperature) */
  if (steps == hc_start) {
    real temp = 2 * fac / (DIM * natoms * (hc_start - hc_av_start));
    fac = 1.0 / (SQRT(volume) * temp);
  }

  if ((myid==0) && (hc_int > 0) && (steps >= hc_start) && 
      ((steps - hc_start) % hc_int == 0)) {
#ifdef TWOD
    hc.x *= fac;  hc.y *= fac;
#else
    hc.x *= fac;  hc.y *= fac;  hc.z *= fac;
#endif
    write_heat_current(steps);
  }

}

#endif /* HC */

#ifdef NVX

/******************************************************************************
*
*  update and write temperature profile
*
******************************************************************************/

void write_temp_dist(int steps)
{
  FILE   *outtemp;
  str255 fnametemp;
  real   scale;
  int    i, k, nhalf;
  static real    *temp_hist_1 = NULL, *temp_hist_2 = NULL;
  static integer  *num_hist_1 = NULL,  *num_hist_2 = NULL;
#ifdef MPI
  static real grad_fit_1[5], grad_fit_2[5];
#else
  static real grad_fit_1[5], *grad_fit_2;
#endif
  
  /* the temp bins are orthogonal boxes in space */
  nhalf = hc_nlayers / 2;
  scale = hc_nlayers / box_x.x;

  /* allocate and clear histogram arrays */
  if (NULL==temp_hist_1) {
    temp_hist_1 = (real    *) malloc( (nhalf+1) * sizeof(real   ) );
     num_hist_1 = (integer *) malloc( (nhalf+1) * sizeof(integer) );
#ifdef MPI
    temp_hist_2 = (real    *) malloc( (nhalf+1) * sizeof(real   ) );
     num_hist_2 = (integer *) malloc( (nhalf+1) * sizeof(integer) );
#else
    temp_hist_2 = temp_hist_1;
     num_hist_2 =  num_hist_1;
     grad_fit_2 =  grad_fit_1;
#endif
    if ( (NULL==temp_hist_1) || (NULL==num_hist_1) 
#ifdef MPI
      || (NULL==temp_hist_2) || (NULL==num_hist_2)
#endif
      ) error("Cannot allocate temperature array.");
    for (i=0; i<=nhalf; i++) {
      temp_hist_1[i] = 0.0;
       num_hist_1[i] = 0;
    }
    for (i=0; i<5; i++) grad_fit_1[i] = 0.0;
    /* write header of temperature gradient file */
    if ((0==myid) && (0==imdrestart)) {
      sprintf(fnametemp, "%s.hcgrad", outfilename);
      outtemp = fopen(fnametemp, "w");
      if (NULL == outtemp) error("Cannot open temperature gradient file.");
      fprintf(outtemp, "# count gradT deltaT kappa kappa[W/mK]\n");
      fclose(outtemp);
    }
    /* write header of temperature profile file */
    if ((0==myid) && (0==imdrestart)) {
      sprintf(fnametemp, "%s.hcprof", outfilename);
      outtemp = fopen(fnametemp, "w");
      if (NULL == outtemp) error("Cannot open temperature profile file.");
      fprintf(outtemp, "# %d %14.4e\n", nhalf+1, hc_heatcurr);
      fclose(outtemp);
    }
  }

  /* construct temperature profile */
  for (k=0; k<NCELLS; ++k) {

    cell *p = CELLPTR(k);

    for (i=0; i<p->n; ++i) {

      real xx, temp;
      int  num;

      /* which layer? */
      xx = ORT(p,i,X);
      if (xx<0.0) xx += box_x.x;
      num = (int) (scale * xx);
      if (num >= hc_nlayers) num -= hc_nlayers;
      if (num > nhalf) {
        num = hc_nlayers - num;
        xx  = box_x.x - xx + box_x.x / hc_nlayers; 
      }
      temp = SPRODN(IMPULS,p,i,IMPULS,p,i) / (2*MASSE(p,i));
      temp_hist_1[num] += temp;
       num_hist_1[num]++;

      /* data for temperature gradient fitting */
      if ((num>2) && (num<nhalf-2)) {
        grad_fit_1[0] += xx;
        grad_fit_1[1] += temp;
        grad_fit_1[2] += temp * xx;
        grad_fit_1[3] += xx * xx;
        grad_fit_1[4] += 1.0;
      }
    }
  }

  /* write temperature profile, and clear arrays afterwards */
  if (0 == steps % hc_int) {

#ifdef MPI
    /* add up results form different CPUs */
    MPI_Reduce(temp_hist_1, temp_hist_2, nhalf+1, REAL,    MPI_SUM,0,cpugrid);
    MPI_Reduce(num_hist_1,  num_hist_2,  nhalf+1, INTEGER, MPI_SUM,0,cpugrid);
    MPI_Reduce(grad_fit_1,  grad_fit_2,        5, REAL,    MPI_SUM,0,cpugrid);
#endif

    /* write temperature profile */
    if (myid==0) {

      real Sxi, STi, SxiTi, Sxi2, a, fact, kappa;

      /* get estimate of temperature gradient */
      Sxi   = grad_fit_2[0] / grad_fit_2[4];
      STi   = grad_fit_2[1] / grad_fit_2[4];
      SxiTi = grad_fit_2[2] / grad_fit_2[4];
      Sxi2  = grad_fit_2[3] / grad_fit_2[4];
      a     = (SxiTi - Sxi * STi) / (Sxi2 - Sxi * Sxi);

      /* write temperature gradient file */
      sprintf(fnametemp, "%s.hcgrad", outfilename);
      outtemp = fopen(fnametemp, "a");
      if (NULL == outtemp) error("Cannot open temperature gradient file.");
      kappa = hc_heatcurr / a;
      /* conversion factor of kappa to SI units */
      fact  = 1.6022e-19 / (1.0179e-14 * 1e-10 * 11605);
      fprintf(outtemp, "%d %10.4e %10.4e %10.4e %10.4e\n", 
              hc_count++, a, 0.5*a*box_x.x, kappa, fact*kappa);
      fclose(outtemp);

      /* open temperature profile file */
      sprintf(fnametemp, "%s.hcprof", outfilename);
      outtemp = fopen(fnametemp, "a");
      if (NULL == outtemp) error("Cannot open temperature profile file.");

      /* write temperature profile */
      fprintf(outtemp, "\n");
      for (i=0; i<=nhalf; i++) {
        if (num_hist_2[i] > 0) temp_hist_2[i] /= num_hist_2[i];
#ifndef TWOD
        temp_hist_2[i] *= (2.0/DIM);
#endif
        fprintf(outtemp, "%10.4e %10.4e\n", (i+0.5) / scale, temp_hist_2[i] );
      }

      /* close file */
      fprintf(outtemp,"\n");
      fclose(outtemp);
    }

    /* clear arrays */
    for (i=0; i<=nhalf; i++) {
      temp_hist_1[i] = 0.0;  
       num_hist_1[i] = 0;
    }
    for (i=0; i<5; i++) grad_fit_1[i] = 0.0;

  }
}

#endif /* NVX */
