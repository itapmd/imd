
/***************************************************************************
 * $RCSfile$
 * $Revision$
 * $Date$
 ***************************************************************************/

#include "imd.h"

#ifdef RNEMD

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

#endif


#ifdef TRANSPORT

/******************************************************************************
*
* write_temp_dist
*
******************************************************************************/

void write_temp_dist(int steps)

{
  FILE  *outtemp;
  str255  fnametemp;
  cell *p;
  real scale;
  int  num, nlayer;
  int i, r, s, t, nhalf;
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

              temp_hist_1[num] += SPRODN(p->impuls,i,p->impuls,i)/ MASSE(p,i);
              num_hist_1[num]++;
	   }
        }

#ifdef MPI
  /* add up results form different CPUs */
  MPI_Reduce( temp_hist_1, temp_hist_2, 
              nhalf+1, MPI_REAL, MPI_SUM, 0, cpugrid);
  temp_hist = temp_hist_2;
  MPI_Reduce( num_hist_1, num_hist_2, 
              nhalf+1, INTEGER,  MPI_SUM, 0, cpugrid);
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
#ifdef NVX
    /* the variable heat_cond was formerly written to the .eng file */
    fprintf(outtemp," %10.4e", heat_cond);
#elif RNEMD
    fprintf(outtemp," %10.4e", heat_transfer);
#endif
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

#ifdef SHOCK

/******************************************************************************
*
* write_press_dist_shock (distribution of stress tensor for shock waves)
*
******************************************************************************/

void write_press_dist_shock(int steps)
{
  FILE  *outpress;
  str255  fnamepress;
  cell *p;
  real scalex, scaley, layvol, laydens;
  int  num, numx, numy, press_dim_tot;
  int fzhlr, i, r, s, t;
#ifndef TWOD
  int scalez, numz;
#endif
  static real  *press_histxx, *press_histxx_1 = NULL, *press_histxx_2 = NULL;
  static real  *press_histyy, *press_histyy_1 = NULL, *press_histyy_2 = NULL;
#ifndef TWOD
  static real  *press_histzz, *press_histzz_1 = NULL, *press_histzz_2 = NULL;
#endif
  static real  *kin_histxx, *kin_histxx_1 = NULL, *kin_histxx_2 = NULL;
  static real  *kin_histxxu, *kin_histxxu_1 = NULL, *kin_histxxu_2 = NULL;
  static real  *kin_histyy, *kin_histyy_1 = NULL, *kin_histyy_2 = NULL;
#ifndef TWOD
  static real  *kin_histzz, *kin_histzz_1 = NULL, *kin_histzz_2 = NULL;
#endif
  static real  *pot_hist, *pot_hist_1 = NULL, *pot_hist_2 = NULL;
  static integer *num_hist,   *num_hist_1 = NULL,  *num_hist_2 = NULL;
  
  /* the press bins are orthogonal boxes in space */

  press_dim.y=1;
  scalex = press_dim.x / box_x.x;
  scaley = press_dim.y / box_y.y;
  press_dim_tot = press_dim.x*press_dim.y;
#ifndef TWOD
  press_dim.z=1;
  scalez = press_dim.z / box_z.z;
  press_dim_tot = press_dim.x*press_dim.y*press_dim.z;
#endif
  layvol = volume / press_dim_tot;

  /* allocate histogram arrays */
  if (NULL==press_histxx_1) {
    press_histxx_1 = (real *) malloc( press_dim_tot * sizeof(real) );
    if (NULL==press_histxx_1) 
      error("Cannot allocate xx pressure tensor array.");
  }  
  if (NULL==press_histyy_1) {
    press_histyy_1 = (real *) malloc( press_dim_tot * sizeof(real) );
    if (NULL==press_histyy_1) 
      error("Cannot allocate yy pressure tensor array.");
  }
#ifndef TWOD  
  if (NULL==press_histzz_1) {
    press_histzz_1 = (real *) malloc( press_dim_tot * sizeof(real) );
    if (NULL==press_histzz_1) 
      error("Cannot allocate zz pressure tensor array.");
  }
#endif  
  if (NULL==kin_histxx_1) {
    kin_histxx_1 = (real *) malloc( press_dim.x * sizeof(real) );
    if (NULL==kin_histxx_1) 
      error("Cannot allocate pressure tensor array.");
  }  
  if (NULL==kin_histxxu_1) {
    kin_histxxu_1 = (real *) malloc( press_dim.x * sizeof(real) );
    if (NULL==kin_histxxu_1) 
      error("Cannot allocate pressure tensor array.");
  }  
  if (NULL==kin_histyy_1) {
    kin_histyy_1 = (real *) malloc( press_dim.x * sizeof(real) );
    if (NULL==kin_histyy_1) 
      error("Cannot allocate pressure tensor array.");
  }  
#ifndef TWOD
  if (NULL==kin_histzz_1) {
    kin_histzz_1 = (real *) malloc( press_dim.x * sizeof(real) );
    if (NULL==kin_histzz_1) 
      error("Cannot allocate pressure tensor array.");
  }  
#endif
  if (NULL==pot_hist_1) {
    pot_hist_1 = (real *) malloc( press_dim_tot * sizeof(real) );
    if (NULL==pot_hist_1) 
      error("Cannot allocate potential energy array.");
  }  
  if (NULL==num_hist_1) {
    num_hist_1 = (integer *) malloc( press_dim_tot * sizeof(integer) );
    if (NULL==num_hist_1) 
      error("Cannot allocate number count array.");
  }
#ifdef MPI
  if (NULL==press_histxx_2) {
    press_histxx_2 = (real *) malloc( press_dim_tot * sizeof(real) );
    if ((NULL==press_histxx_2) && (myid==0)) 
      error("Cannot allocate xx pressure tensor  array.");
  }  
  if (NULL==press_histyy_2) {
    press_histyy_2 = (real *) malloc( press_dim_tot * sizeof(real) );
    if ((NULL==press_histyy_2) && (myid==0)) 
      error("Cannot allocate yy pressure tensor  array.");
  }
#ifndef TWOD
  if (NULL==press_histzz_2) {
    press_histzz_2 = (real *) malloc( press_dim_tot * sizeof(real) );
    if ((NULL==press_histzz_2) && (myid==0)) 
      error("Cannot allocate zz pressure tensor  array.");
  }
#endif
  if (NULL==kin_histxx_2) {
    kin_histxx_2 = (real *) malloc( press_dim.x * sizeof(real) );
    if (NULL==kin_histxx_2) 
      error("Cannot allocate pressure tensor array.");
  }  
  if (NULL==kin_histxxu_2) {
    kin_histxxu_2 = (real *) malloc( press_dim.x * sizeof(real) );
    if (NULL==kin_histxxu_2) 
      error("Cannot allocate pressure tensor array.");
  }  
  if (NULL==kin_histyy_2) {
    kin_histyy_2 = (real *) malloc( press_dim.x * sizeof(real) );
    if (NULL==kin_histyy_2) 
      error("Cannot allocate pressure tensor array.");
  }  
#ifndef TWOD
  if (NULL==kin_histzz_2) {
    kin_histzz_2 = (real *) malloc( press_dim.x * sizeof(real) );
    if (NULL==kin_histzz_2) 
      error("Cannot allocate pressure tensor array.");
  }  
#endif
  if (NULL==pot_hist_2) {
    pot_hist_2 = (real *) malloc( press_dim_tot * sizeof(real) );
    if (NULL==pot_hist_2) 
      error("Cannot allocate potential energy array.");
  }  
  if (NULL==num_hist_2) {
    num_hist_2 = (integer *) malloc( press_dim_tot * sizeof(integer) );
    if ((NULL==num_hist_2) && (myid==0)) 
      error("Cannot allocate number count array.");
  }
#endif

  for (i = 0; i < press_dim.x; i++) {
    kin_histxx_1[i] = 0;
    kin_histxxu_1[i] = 0;
    kin_histyy_1[i] = 0;
#ifndef TWOD
    kin_histzz_1[i] = 0;
#endif
#ifdef MPI
    kin_histxx_2[i] = 0;
    kin_histxxu_2[i] = 0;
    kin_histyy_2[i] = 0;
#ifndef TWOD
    kin_histzz_2[i] = 0;
#endif
#endif
  }
  for (i = 0; i < press_dim_tot; i++) {
    press_histxx_1[i] = 0.0;
    press_histyy_1[i] = 0.0;
#ifndef TWOD
    press_histzz_1[i] = 0.0;
#endif
    pot_hist_1[i] = 0;
    num_hist_1[i] = 0;
#ifdef MPI
    press_histxx_2[i] = 0.0;
    press_histyy_2[i] = 0.0;
#ifndef TWOD
    press_histzz_2[i] = 0.0;
#endif
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
              if (numx >= press_dim.x)      numx = press_dim.x-1;
	      numy = scaley * p->ort Y(i);
              if (numy < 0)             numy = 0;
              if (numy >= press_dim.y)      numy = press_dim.y-1;
              num = numx * press_dim.y + numy;
#ifndef TWOD
	      numz = scalez * p->ort Z(i);
              if (numz < 0)             numz = 0;
              if (numz >= press_dim.z)      numz = press_dim.z-1;
              num = num  * press_dim.z + numz;
#endif
#ifdef TWOD
              press_histxx_1[num] += p->presstens X(i);
              press_histyy_1[num] += p->presstens Y(i);

              /* average v_xx relative to center of mass*/

              kin_histxx_1[numx] += p->impuls X(i) * p->impuls X(i) 
                / (2* MASSE(p,i));

              /* average v_xx - u_p  relative to moving pistons */
              if (shock_mode == 2) {
                if ( p->ort X(i) < box_x.x*0.5 )
                  kin_histxxu_1[numx]+=(p->impuls X(i)-shock_speed* MASSE(p,i)) 
                    *(p->impuls X(i)-shock_speed* MASSE(p,i))/(2* MASSE(p,i));
                else
                  kin_histxxu_1[numx]+=(p->impuls X(i)+shock_speed* MASSE(p,i)) 
                    *(p->impuls X(i)+shock_speed* MASSE(p,i))/(2* MASSE(p,i));
              }
              if (shock_mode == 1) {
                if ( p->ort X(i) < shock_strip ) 
		  kin_histxxu_1[numx]+=
		    (p->impuls X(i)-shock_speed* MASSE(p,i)) 
		    *(p->impuls X(i)-shock_speed* MASSE(p,i))
		    /(2* MASSE(p,i));
		else
		  kin_histxxu_1[numx]+=p->impuls X(i)*p->impuls X(i)
		    /(2* MASSE(p,i));
              }

              kin_histyy_1[numx] += p->impuls Y(i) * p->impuls Y(i) 
                / (2* MASSE(p,i));

	      /* ENDE Johannes Stuff */

	      pot_hist_1[num] += p->pot_eng[i];
              num_hist_1[num]++;
#else
              press_histxx_1[num] += p->presstens X(i);
              press_histyy_1[num] += p->presstens Y(i);
              press_histzz_1[num] += p->presstens Z(i);

              /* average v_xx relative to center of mass*/

              kin_histxx_1[numx] += p->impuls X(i) * p->impuls X(i) 
                / (2* MASSE(p,i));

              /* average v_xx - u_p  relative to moving pistons */
              if (shock_mode == 2) {
                if ( p->ort X(i) < box_x.x*0.5 )
                  kin_histxxu_1[numx]+=(p->impuls X(i)-shock_speed* MASSE(p,i)) 
                    *(p->impuls X(i)-shock_speed* MASSE(p,i))/(2* MASSE(p,i));
                else
                  kin_histxxu_1[numx]+=(p->impuls X(i)+shock_speed* MASSE(p,i)) 
                    *(p->impuls X(i)+shock_speed* MASSE(p,i))/(2* MASSE(p,i));
              }
              if (shock_mode == 1) {
                if ( p->ort X(i) < shock_strip ) 
		  kin_histxxu_1[numx]+=
		    (p->impuls X(i)-shock_speed* MASSE(p,i)) 
		    *(p->impuls X(i)-shock_speed* MASSE(p,i))
		    /(2* MASSE(p,i));
		else
		  kin_histxxu_1[numx]+=p->impuls X(i)*p->impuls X(i)
		    /(2* MASSE(p,i));
              }

              kin_histyy_1[numx] += p->impuls Y(i) * p->impuls Y(i) 
                / (2* MASSE(p,i));
#ifndef TWOD
              kin_histzz_1[numx] += p->impuls Z(i) * p->impuls Z(i) 
                / (2* MASSE(p,i));
#endif

	      pot_hist_1[num] += p->pot_eng[i];
              num_hist_1[num]++;
#endif
	   }
        }

#ifdef MPI
  /* add up results form different CPUs */
  MPI_Reduce( press_histxx_1, press_histxx_2, 
              press_dim_tot, MPI_REAL, MPI_SUM, 0, cpugrid);
  press_histxx = press_histxx_2;
  MPI_Reduce( press_histyy_1, press_histyy_2, 
              press_dim_tot, MPI_REAL, MPI_SUM, 0, cpugrid);
  press_histyy = press_histyy_2;
#ifndef TWOD
  MPI_Reduce( press_histzz_1, press_histzz_2, 
              press_dim_tot, MPI_REAL, MPI_SUM, 0, cpugrid);
  press_histzz = press_histzz_2;
#endif
  MPI_Reduce( pot_hist_1, pot_hist_2, 
              press_dim_tot, MPI_REAL, MPI_SUM, 0, cpugrid);
  pot_hist = pot_hist_2;
  MPI_Reduce( num_hist_1, num_hist_2, 
              press_dim_tot, INTEGER, MPI_SUM, 0, cpugrid);
  num_hist = num_hist_2;
  MPI_Reduce( kin_histxx_1, kin_histxx_2, 
              press_dim.x, MPI_REAL, MPI_SUM, 0, cpugrid);
  kin_histxx = kin_histxx_2;

  MPI_Reduce( kin_histxxu_1, kin_histxxu_2, 
              press_dim.x, MPI_REAL, MPI_SUM, 0, cpugrid);
  kin_histxxu = kin_histxxu_2;

  MPI_Reduce( kin_histyy_1, kin_histyy_2, 
              press_dim.x, MPI_REAL, MPI_SUM, 0, cpugrid);
  kin_histyy = kin_histyy_2;
#ifndef TWOD
  MPI_Reduce( kin_histzz_1, kin_histzz_2, 
              press_dim.x, MPI_REAL, MPI_SUM, 0, cpugrid);
  kin_histzz = kin_histzz_2;
#endif
  MPI_Reduce( pot_hist_1, pot_hist_2, 
              press_dim.x, MPI_REAL, MPI_SUM, 0, cpugrid);
  pot_hist = pot_hist_2;
  MPI_Reduce( num_hist_1, num_hist_2, 
              press_dim.x, INTEGER, MPI_SUM, 0, cpugrid);
  num_hist  = num_hist_2;
#else
  press_histxx = press_histxx_1;
  press_histyy = press_histyy_1;
#ifndef TWOD
  press_histzz = press_histzz_1;
#endif
  kin_histxx = kin_histxx_1;
  kin_histxxu = kin_histxxu_1;
  kin_histyy = kin_histyy_1;
#ifndef TWOD
  kin_histzz = kin_histzz_1;
#endif
  pot_hist = pot_hist_1;
  num_hist = num_hist_1;
#endif

  /* write pressure tensor distribution */

  if (myid==0) {

    fzhlr = steps / press_interval;

    sprintf(fnamepress,"%s.%u.pressdist",outfilename,fzhlr);
    outpress = fopen(fnamepress,"w");
    if (NULL == outpress) error("Cannot open pressure tensor file.");
#ifdef SHOCK
    for (i = 0; i < press_dim.x; i++) {
      if (num_hist[i] > 0) {
	laydens = num_hist[i] / layvol;
	pot_hist[i] /= num_hist[i];
	press_histxx[i] /= num_hist[i];
	press_histyy[i] /= num_hist[i];
#ifndef TWOD
	press_histzz[i] /= num_hist[i];
#endif
        kin_histxx[i] /= num_hist[i];
        kin_histxxu[i] /= num_hist[i];
        kin_histyy[i] /= num_hist[i];
#ifndef TWOD
        kin_histzz[i] /= num_hist[i];
#endif
      }

#ifndef TWOD
      fprintf(outpress,"%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10e\n", 
              press_histxx[i], press_histyy[i], press_histzz[i], laydens, 
              kin_histxx[i], kin_histxxu[i], kin_histyy[i],
              kin_histzz[i], pot_hist[i] );
#else
      fprintf(outpress,"%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n", 
              press_histxx[i], press_histyy[i], laydens,
              kin_histxx[i], kin_histxxu[i], kin_histyy[i], 
              pot_hist[i] ); 
#endif
      }
#endif

    fprintf(outpress,"\n");
    fclose(outpress);

  }
}

#else /* not SHOCK */

/******************************************************************************
*
* write_press_dist (distribution of stress tensor)
*
******************************************************************************/

void write_press_dist(int steps)
{
  FILE  *outpress;
  str255  fnamepress;
  cell *p;
  real scalex, scaley, scalez, layvol, laydens, Ekin;
  int  num, numx, numy, numz, press_dim_tot;
  int fzhlr, i, r, s, t;
  static real  *press_histxx, *press_histxx_1 = NULL, *press_histxx_2 = NULL;
  static real  *press_histyy, *press_histyy_1 = NULL, *press_histyy_2 = NULL;
  static real  *press_histxy, *press_histxy_1 = NULL, *press_histxy_2 = NULL;
#ifndef TWOD
  static real  *press_histzz, *press_histzz_1 = NULL, *press_histzz_2 = NULL;
  static real  *press_histzx, *press_histzx_1 = NULL, *press_histzx_2 = NULL;
  static real  *press_histyz, *press_histyz_1 = NULL, *press_histyz_2 = NULL;
#endif
  static real  *kin_hist, *kin_hist_1 = NULL, *kin_hist_2 = NULL;
  static real  *pot_hist, *pot_hist_1 = NULL, *pot_hist_2 = NULL;
  static integer *num_hist,   *num_hist_1 = NULL,  *num_hist_2 = NULL;
  
  /* the press bins are orthogonal boxes in space */
  scalex = press_dim.x / box_x.x;
  scaley = press_dim.y / box_y.y;
  press_dim_tot = press_dim.x*press_dim.y;
#ifndef TWOD
  scalez = press_dim.z / box_z.z;
  press_dim_tot *= press_dim.z;
#endif
  layvol = volume / press_dim_tot;

  /* allocate histogram arrays */
  if (NULL==press_histxx_1) {
    press_histxx_1 = (real *) malloc( press_dim_tot * sizeof(real) );
    if (NULL==press_histxx_1) 
      error("Cannot allocate xx pressure tensor array.");
  }  
  if (NULL==press_histyy_1) {
    press_histyy_1 = (real *) malloc( press_dim_tot * sizeof(real) );
    if (NULL==press_histyy_1) 
      error("Cannot allocate yy pressure tensor array.");
  }
  if (NULL==press_histxy_1) {
    press_histxy_1 = (real *) malloc( press_dim_tot * sizeof(real) );
    if (NULL==press_histxy_1) 
      error("Cannot allocate xy pressure tensor array.");
  }
#ifndef TWOD  
  if (NULL==press_histzz_1) {
    press_histzz_1 = (real *) malloc( press_dim_tot * sizeof(real) );
    if (NULL==press_histzz_1) 
      error("Cannot allocate zz pressure tensor array.");
  }
  if (NULL==press_histzx_1) {
    press_histzx_1 = (real *) malloc( press_dim_tot * sizeof(real) );
    if (NULL==press_histzx_1) 
      error("Cannot allocate zx pressure tensor array.");
  }
  if (NULL==press_histyz_1) {
    press_histyz_1 = (real *) malloc( press_dim_tot * sizeof(real) );
    if (NULL==press_histyz_1) 
      error("Cannot allocate yz pressure tensor array.");
  }
#endif  
  if (NULL==kin_hist_1) {
    kin_hist_1 = (real *) malloc( press_dim_tot * sizeof(real) );
    if (NULL==kin_hist_1) 
      error("Cannot allocate kinetic energy array.");
  }
  if (NULL==pot_hist_1) {
    pot_hist_1 = (real *) malloc( press_dim_tot * sizeof(real) );
    if (NULL==pot_hist_1) 
      error("Cannot allocate potential energy array.");
  }  
  if (NULL==num_hist_1) {
    num_hist_1 = (integer *) malloc( press_dim_tot * sizeof(integer) );
    if (NULL==num_hist_1) 
      error("Cannot allocate number count array.");
  }
#ifdef MPI
  if (NULL==press_histxx_2) {
    press_histxx_2 = (real *) malloc( press_dim_tot * sizeof(real) );
    if ((NULL==press_histxx_2) && (myid==0)) 
      error("Cannot allocate xx pressure tensor  array.");
  }  
  if (NULL==press_histyy_2) {
    press_histyy_2 = (real *) malloc( press_dim_tot * sizeof(real) );
    if ((NULL==press_histyy_2) && (myid==0)) 
      error("Cannot allocate yy pressure tensor  array.");
  }
  if (NULL==press_histxy_2) {
    press_histxy_2 = (real *) malloc( press_dim_tot * sizeof(real) );
    if ((NULL==press_histxy_2) && (myid==0)) 
      error("Cannot allocate xy pressure tensor  array.");
  }
#ifndef TWOD
  if (NULL==press_histzz_2) {
    press_histzz_2 = (real *) malloc( press_dim_tot * sizeof(real) );
    if ((NULL==press_histzz_2) && (myid==0)) 
      error("Cannot allocate zz pressure tensor  array.");
  }
  if (NULL==press_histzx_2) {
    press_histzx_2 = (real *) malloc( press_dim_tot * sizeof(real) );
    if ((NULL==press_histzx_2) && (myid==0)) 
      error("Cannot allocate zx pressure tensor  array.");
  }
  if (NULL==press_histyz_2) {
    press_histyz_2 = (real *) malloc( press_dim_tot * sizeof(real) );
    if ((NULL==press_histyz_2) && (myid==0)) 
      error("Cannot allocate yz pressure tensor  array.");
  }
#endif
  if (NULL==kin_hist_2) {
    kin_hist_2 = (real *) malloc( press_dim_tot * sizeof(real) );
    if (NULL==kin_hist_2) 
      error("Cannot allocate kinetic energy array.");
  }  
  if (NULL==pot_hist_2) {
    pot_hist_2 = (real *) malloc( press_dim_tot * sizeof(real) );
    if (NULL==pot_hist_2) 
      error("Cannot allocate potential energy array.");
  }  
  if (NULL==num_hist_2) {
    num_hist_2 = (integer *) malloc( press_dim_tot * sizeof(integer) );
    if ((NULL==num_hist_2) && (myid==0)) 
      error("Cannot allocate number count array.");
  }
#endif /* MPI */

  /* clear histograms */
  for (i = 0; i < press_dim_tot; i++) {
    press_histxx_1[i] = 0.0;
    press_histyy_1[i] = 0.0;
    press_histxy_1[i] = 0.0;
#ifndef TWOD
    press_histzz_1[i] = 0.0;
    press_histzx_1[i] = 0.0;
    press_histyz_1[i] = 0.0;
#endif
    kin_hist_1[i] = 0; 
    pot_hist_1[i] = 0;
    num_hist_1[i] = 0;
#ifdef MPI
    press_histxx_2[i] = 0.0;
    press_histyy_2[i] = 0.0;
    press_histxy_2[i] = 0.0;
#ifndef TWOD
    press_histzz_2[i] = 0.0;
    press_histzx_2[i] = 0.0;
    press_histyz_2[i] = 0.0;
#endif
    kin_hist_2[i] = 0;
    pot_hist_2[i] = 0;
    num_hist_2[i] = 0;
#endif /* MPI */
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
           for (i = 0; i < p->n; ++i) {
              /* which layer? */
	      numx = scalex * p->ort X(i);
              if (numx < 0)             numx = 0;
              if (numx >= press_dim.x)      numx = press_dim.x-1;
	      numy = scaley * p->ort Y(i);
              if (numy < 0)             numy = 0;
              if (numy >= press_dim.y)      numy = press_dim.y-1;
              num = numx * press_dim.y + numy;
#ifndef TWOD
	      numz = scalez * p->ort Z(i);
              if (numz < 0)             numz = 0;
              if (numz >= press_dim.z)      numz = press_dim.z-1;
              num = num  * press_dim.z + numz;
#endif

              Ekin = SPRODN(p->impuls,i,p->impuls,i) / (2* MASSE(p,i));
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
	      kin_hist_1[num] += Ekin; 
	      pot_hist_1[num] += p->pot_eng[i];
              num_hist_1[num]++;
	   }
        }

#ifdef MPI
  /* add up results form different CPUs */
  MPI_Reduce( press_histxx_1, press_histxx_2, 
              press_dim_tot, MPI_REAL, MPI_SUM, 0, cpugrid);
  press_histxx = press_histxx_2;
  MPI_Reduce( press_histyy_1, press_histyy_2, 
              press_dim_tot, MPI_REAL, MPI_SUM, 0, cpugrid);
  press_histyy = press_histyy_2;
  MPI_Reduce( press_histxy_1, press_histxy_2, 
              press_dim_tot, MPI_REAL, MPI_SUM, 0, cpugrid);
  press_histxy = press_histxy_2;
#ifndef TWOD
  MPI_Reduce( press_histzz_1, press_histzz_2, 
              press_dim_tot, MPI_REAL, MPI_SUM, 0, cpugrid);
  press_histzz = press_histzz_2;
  MPI_Reduce( press_histzx_1, press_histzx_2, 
              press_dim_tot, MPI_REAL, MPI_SUM, 0, cpugrid);
  press_histzx = press_histzx_2;
  MPI_Reduce( press_histyz_1, press_histyz_2, 
              press_dim_tot, MPI_REAL, MPI_SUM, 0, cpugrid);
  press_histyz = press_histyz_2;
#endif
  MPI_Reduce( kin_hist_1, kin_hist_2, 
              press_dim_tot, MPI_REAL, MPI_SUM, 0, cpugrid);
  kin_hist     = kin_hist_2;
  MPI_Reduce( pot_hist_1, pot_hist_2, 
              press_dim_tot, MPI_REAL, MPI_SUM, 0, cpugrid);
  pot_hist     = pot_hist_2;
  MPI_Reduce( num_hist_1, num_hist_2, 
              press_dim_tot, INTEGER, MPI_SUM, 0, cpugrid);
  num_hist     = num_hist_2;
#else /* not MPI */
  press_histxx = press_histxx_1;
  press_histyy = press_histyy_1;
  press_histxy = press_histxy_1;
#ifndef TWOD
  press_histzz = press_histzz_1;
  press_histzx = press_histzx_1;
  press_histyz = press_histyz_1;
#endif
  kin_hist     = kin_hist_1; 
  pot_hist     = pot_hist_1;
  num_hist     = num_hist_1;
#endif /* not MPI */

  /* write pressure tensor distribution */
  if (myid==0) {

    fzhlr = steps / press_interval;
    sprintf(fnamepress,"%s.%u.pressdist",outfilename,fzhlr);
    outpress = fopen(fnamepress,"w");
    if (NULL == outpress) error("Cannot open pressure tensor file.");

    for (i = 0; i < press_dim_tot; i++) {
      if (num_hist[i] > 0) {
	laydens = num_hist[i] / layvol;
	kin_hist[i]     /= num_hist[i];
	pot_hist[i]     /= num_hist[i];
	press_histxx[i] /= num_hist[i];
	press_histyy[i] /= num_hist[i];
	press_histxy[i] /= num_hist[i];
#ifndef TWOD
	press_histzz[i] /= num_hist[i];
	press_histzx[i] /= num_hist[i];
	press_histyz[i] /= num_hist[i];
#endif
      }
    }
      
#ifdef TWOD
    i=0;
    for (r = 0; r < press_dim.x; r++) {
      for (s = 0; s < press_dim.y; s++) {
        fprintf(outpress,
                "%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n", 
	        press_histxx[i], press_histyy[i], press_histxy[i],
                laydens, kin_hist[i], pot_hist[i] );
        i++;
      }
      fprintf(outpress,"\n");
    }
#else
    i=0;
    for (r = 0; r < press_dim.x; r++) {
      for (s = 0; s < press_dim.y; s++) {
        for (t = 0; t < press_dim.z; t++) {
          fprintf(outpress,
              "%10.4e %10.4e %10.4e %10e %10.4e %10.4e %10e %10.4e %10.4e\n", 
	      press_histxx[i], press_histyy[i], press_histzz[i], 
	      press_histyz[i], press_histzx[i], press_histxy[i],
              laydens, kin_hist[i], pot_hist[i] );
          i++;
	}
        fprintf(outpress,"\n");
      }
      fprintf(outpress,"\n");
    }
#endif
    fprintf(outpress,"\n");
    fclose(outpress);
  }
}

#endif /* not SHOCK */


/******************************************************************************
*
* write_press_atoms (pressure tensor for each atom)
*
******************************************************************************/

void write_press_atoms(int steps)
{
  FILE  *outpress;
  str255  fnamepress;
  cell *p;
  int fzhlr, i, r, s, t;

  fzhlr = steps / press_interval;

#ifdef MPI
  sprintf(fnamepress,"%s.%u.pressdist.%u",outfilename,fzhlr,myid);
#else
  sprintf(fnamepress,"%s.%u.pressdist",outfilename,fzhlr);
#endif
  outpress = fopen(fnamepress,"w");
  if (NULL == outpress) error("Cannot open pressure tensor file.");

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
        for (i = 0; i < p->n; ++i) {
#ifdef TWOD
          fprintf(outpress,"%10.4e %10.4e %10.4e %10.4e %10.4e\n", 
                  p->ort X(i),p->ort Y(i),
                  p->presstens X(i),p->presstens Y(i),
                  p->presstens_offdia[i]);
#else
          fprintf(outpress,"%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n", 
                  p->ort X(i),p->ort Y(i),p->ort Z(i),
                  p->presstens X(i),p->presstens Y(i),p->presstens Z(i),
                  p->presstens_offdia X(i),p->presstens_offdia Y(i),
                  p->presstens_offdia Z(i));
#endif
        }
      }
  fclose(outpress);
}

#endif /* STRESS_TENS */









