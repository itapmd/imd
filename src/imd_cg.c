/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2008 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_cg.c -- conjugate gradient stuff
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"
#include "cg_utils.h"

/*****************************************************************************
*
*  reset cg integrator
*
*****************************************************************************/

void reset_cg(void)
{
  /* initialisations: h, g, f_max */
  calc_forces(0);
  calc_fnorm_g_h();
  //cg_poteng = tot_pot_energy / natoms; calculate with larger numbers
  cg_poteng = tot_pot_energy  ;
  old_cg_poteng = cg_poteng;
}

/*****************************************************************************
*
*  one cg line minimization
*
*****************************************************************************/
#if defined (CG)&& !defined(ACG)
void cg_step(int steps)
{
  /* minimization in one 'direction' */
  if (cg_reset_int>0) {
    if (0==steps%cg_reset_int) reset_cg();
  }
  cg_poteng = tot_pot_energy ;
  old_cg_poteng = cg_poteng;

  //  if ((cg_infolevel>0) && (0==myid)) { printf("fnorm old = %e\n", SQRT( fnorm / nactive ) );}
  
  linmin();
  
  // if ((cg_infolevel>0) && (0==myid)) { printf("fnorm new = %e\n", SQRT( fnorm / nactive ) );}
  
  cg_calcgamma(); /* calc gg, dgg */

  /* sets old_ort = ort, h, g, needs gamma and gets f_max2 */ 
  set_hg();
  
}
#endif

#ifdef ACG
void acg_step(int steps)
{
  /* minimization in one 'direction' */
  if (cg_reset_int>0) {
    if (0==steps%cg_reset_int) reset_cg();
  }
  cg_poteng = tot_pot_energy ;
  old_cg_poteng = cg_poteng;

  findalpha();
  
  cg_calcgamma(); /* calc gg, dgg */

  /* sets old_ort = ort, h, g, needs gamma and gets f_max2 */ 
  set_hg();

  acg_alpha *= acg_incfac;
}
#endif

/******************************************************************************
 *
 * linmin : line minimization for cg
 *
 *****************************************************************************/

int linmin()
{
  real alpha_a, alpha_b, alpha_c, fa, fb, fc,dum;
  real alphamin;
  static real old_alphamin = 0.1;
  real f_max;
  int  iter1=0, iter2=0;
  
  /* can we do the scaling in an other way, additional lowest scale to move? */
  f_max = sqrt(f_max2);
  /* alpha_b = linmin_dmax/sqrt(f_max2);  got short distances: new scaling*/  

  /* scaling: alpha_b should be 0.01 and alpha_c 0.02 */
  /* this can be achived by setting linmin_dmax to 0.001 */

  if(f_max > 100*linmin_dmax) {
    alpha_b = linmin_dmax/f_max;
  }
  // else if (f_max<linmin_dmin) {
    //alpha_b = linmin_dmax *linmin_dmin/f_max;
  else if (old_alphamin<linmin_dmin) {

    //     alpha_b=2.0 * linmin_dmin/old_alphamin;
    //alpha_b=2.0 * 0.1;
    // alpha_b=2.0 * linmin_dmin/old_alphamin;
    alpha_b=0.05;

  }
  else {
    alpha_b = 2.0* old_alphamin;  
    //    alpha_b = linmin_dmax;
  }
  
  
  // this is anyhow set in mnbrak
  alpha_c = 2.0 * alpha_b;
  alpha_a = 0.0;

  fa = old_cg_poteng;
  
  fb = fonedim(alpha_b);

  /* decide which method to take to braket a mimimum, */
  /* at the moment only mbrak, later zbrak? */
  /* call by reference Num Rec. p297 */
  
  if ((cg_infolevel>0) && (0==myid)) {
    printf("befor mnbrak alpha_a %le alpha_b %le alpha_c %le %.12lf %.12lf\n",
           alpha_a,alpha_b,alpha_c,fa,fb);
    fflush(stdout);
  }

  iter1 = mnbrak(&alpha_a,&alpha_b,&alpha_c,&fa,&fb,&fc); 


  if (alpha_a > alpha_b) {
     SHFT(dum,alpha_a,alpha_b,dum)
     SHFT(dum,fa,fb,dum)
  }
  if (alpha_b > alpha_c) {
    SHFT(dum,alpha_b,alpha_c,dum)
      SHFT(dum,fb,fc,dum)
  }
  if (alpha_a > alpha_b) {
    SHFT(dum,alpha_a,alpha_b,dum)
      SHFT(dum,fa,fb,dum)
  }
  
  if ((cg_infolevel>0) && (0==myid)) {
    printf("after mnbrak alpha_a %le alpha_b %le alpha_c %le %.12lf %.12lf %.12lf\n",
           alpha_a,alpha_b,alpha_c,fa,fb,fc);
    fflush(stdout);
  }

  if(iter1 <1) {
    printf("error in mnbrak %lf %lf %lf %.12lf %.12lf %.12lf\n",
           alpha_a,alpha_b,alpha_c,fa,fb,fc);
    fflush(stdout);
  }

/*   if ((cg_infolevel>0) && (0==myid)) { */
/*   printf("after mnbrak alpha_a %le alpha_b %le alpha_c %le %.12lf %.12lf %.12lf\n", */
/*            alpha_a,alpha_b,alpha_c,fa,fb,fc); */
/*   fflush(stdout); */
/*   } */


  /* decide which method to take to search the mimimum */
  if (cg_mode == CGEF) { /* not implemented yet */
    /* in ort should be the coord. of min pos. */
    /* iter2 = dbrent (alpha_a,alpha_b,alpha_c,fa,fb,fc); */ 
  }
  else {
    iter2 = brent(alpha_a,alpha_b,alpha_c,fa,fb,fc,&alphamin);
  }
  old_alphamin = alphamin;

  /* info message */
  if ((cg_infolevel>0) && (0==myid)) {
    printf("iter1= %d iter2 = %d alphamin = %e f_max = %e fnorm %e  epot %.12e \n",
           iter1,iter2,alphamin,f_max,sqrt(fnorm/nactive),tot_pot_energy );
    fflush(stdout); 
  }
  return (iter1 + iter2);
}


#ifdef ACG

/******************************************************************************
 *
 * findalpha: simplified line minimization for acg from Kai Nordlund
 *
 *****************************************************************************/

/* int findalpha() */
/* { */

/*   int  iter=0; */
/*   real fb; */

/*   fb =old_cg_poteng; */
/*   while (fb >= old_cg_poteng) */
/*     { */
/*       fb = fonedim(acg_alpha); */
/*       iter++; */
/*       if (fb>=old_cg_poteng) */
/* 	{ */
/* 	  acg_alpha = (acg_alpha<=1e-10 || iter>50) ? acg_alpha *acg_incfac : acg_alpha * acg_decfac; */
/* 	} */
/*       if (iter > 100) */
/* 	{ */
/* 	  printf("iter > 100"); */
/* 	  reset_cg(); */
/* 	  break; */
/* 	} */
/*     } */

/*   return (iter); */
/* } */

int findalpha()
{

  int  iter=0;
  real fb;

  fb =old_cg_poteng;
  while (fb >= old_cg_poteng)
    {
	if (iter <4)     // this is ACG as in Kain Nordlunds script
	    {
		fb = fonedim(acg_alpha);
		iter++;
		if (fb>=old_cg_poteng)
		    {
			acg_alpha =  acg_alpha * acg_decfac;
		    }
	    }
	else if (iter>50)  // bail-out by Kai Nordlund
	    {
		write_eng_file(nfc+iter);
		write_cgconfig(nfc+iter);

		if (myid==0) {
		  /*printf("ACG not able to find step minimizing fonedim\n");*/
		  printf("nfc = %d epot = %22.16f\n", nfc, old_cg_poteng );
		  error("ACG not able to find step minimizing fonedim");
		}

	    }

	else          
	    {
	      if (iter==4)
			acg_alpha=acg_init_alpha;
	      fb = fonedim_sd(acg_alpha);
	      iter++;
	      if (fb>=old_cg_poteng)
		{
		  acg_alpha =  acg_alpha * acg_decfac;
		}
	    }
    
	
    }
  return (iter);
}

#endif







/* ondimensional function, used by brent */
/* sets the global variables epot,fnorm corresponding to alpha */
real fonedim(real alpha)  
{
  move_atoms_cg(alpha);
#ifdef NBLIST
  check_nblist();
#else
  fix_cells();
#endif
  calc_forces(1);
  calc_fnorm();
  // return tot_pot_energy / natoms;
  return tot_pot_energy ;
}


/* ondimensional function, used by brent */
/* sets the global variables epot,fnorm corresponding to alpha */
real fonedim_sd(real alpha)  
{
  move_atoms_sd(alpha);
#ifdef NBLIST
  check_nblist();
#else
  fix_cells();
#endif
  calc_forces(1);
  calc_fnorm();
  // return tot_pot_energy / natoms;
  return tot_pot_energy ;
}


/* calculates the forcenorm - the forces have to be calculated before! */
void calc_fnorm(void)
{
  int k;
  real tmp_fnorm=0.0;
  fnorm = 0.0;

  /* loop over all cells */
  for (k=0; k<NCELLS; ++k) {

    int  i, sort;
    cell *p;
    p = CELLPTR(k);

    /* loop over all particles */
    for (i=0; i<p->n; ++i) {

      sort = VSORTE(p,i);

#ifdef FBC
      /* give virtual particles their extra force */
      KRAFT(p,i,X) += (fbc_forces + sort)->x;
      KRAFT(p,i,Y) += (fbc_forces + sort)->y;
#ifndef TWOD
      KRAFT(p,i,Z) += (fbc_forces + sort)->z;
#endif
#endif

      /* and set their force in restricted directions to 0 */ 
      KRAFT(p,i,X) *= (restrictions + sort)->x;
      KRAFT(p,i,Y) *= (restrictions + sort)->y;
#ifndef TWOD
      KRAFT(p,i,Z) *= (restrictions + sort)->z;
#endif

      tmp_fnorm +=  SPRODN(KRAFT,p,i,KRAFT,p,i);
    }
  }

#ifdef MPI
  /* add up results from different CPUs */
  MPI_Allreduce( &tmp_fnorm, &fnorm, 1, REAL, MPI_SUM, cpugrid);
#else
  fnorm = tmp_fnorm;
#endif

}

void calc_fnorm_g_h(void)
{
  int k;
  real tmp_fnorm=0.0, tmp_f_max2=0.0;
  fnorm = 0.0;
  
  /* loop over all cells */
  for (k=0; k<NCELLS; ++k) {

    int  i, sort;
    cell *p;
    p = CELLPTR(k);

    /* loop over all particles */
    for (i=0; i<p->n; ++i) {

      sort = VSORTE(p,i);
#ifdef FBC
      /* give virtual particles their extra force */
      KRAFT(p,i,X) += (fbc_forces + sort)->x;
      KRAFT(p,i,Y) += (fbc_forces + sort)->y;
#ifndef TWOD
      KRAFT(p,i,Z) += (fbc_forces + sort)->z;
#endif
#endif

      /* and set their force (->momentum) in restricted directions to 0 */ 
      KRAFT(p,i,X) *= (restrictions + sort)->x;
      KRAFT(p,i,Y) *= (restrictions + sort)->y;
#ifndef TWOD
      KRAFT(p,i,Z) *= (restrictions + sort)->z;
#endif

      tmp_fnorm +=  SPRODN(KRAFT,p,i,KRAFT,p,i);
      /* initialise old_ort */
      OLD_ORT(p,i,X) = ORT(p,i,X);
      OLD_ORT(p,i,Y) = ORT(p,i,Y);
#ifndef TWOD
      OLD_ORT(p,i,Z) = ORT(p,i,Z);
#endif

      /* initialise search vectors */
      CG_H(p,i,X) = KRAFT(p,i,X);
      CG_H(p,i,Y) = KRAFT(p,i,Y);
#ifndef TWOD
      CG_H(p,i,Z) = KRAFT(p,i,Z);
#endif

      CG_G(p,i,X) = KRAFT(p,i,X);
      CG_G(p,i,Y) = KRAFT(p,i,Y);
#ifndef TWOD
      CG_G(p,i,Z) = KRAFT(p,i,Z);
#endif

      /* determine the biggest force component */
      tmp_f_max2 = MAX(SQR(KRAFT(p,i,X)),tmp_f_max2);
      tmp_f_max2 = MAX(SQR(KRAFT(p,i,Y)),tmp_f_max2);
#ifndef TWOD
      tmp_f_max2 = MAX(SQR(KRAFT(p,i,Z)),tmp_f_max2);
#endif
    }
  }

#ifdef MPI
  /* find the maximum from all CPUs */
  MPI_Allreduce( &tmp_f_max2, &f_max2, 1, REAL, MPI_MAX, cpugrid);
  /* add up results from different CPUs */
  MPI_Allreduce( &tmp_fnorm, &fnorm, 1, REAL, MPI_SUM, cpugrid);
#else
  fnorm = tmp_fnorm;
  f_max2 = tmp_f_max2;
#endif

}

/* calculation of gamma, dgg,gg */
void cg_calcgamma(void)
{
  int k;
  vektor tmpvec;
  real tmp_f_max2 = 0.0, tmp_gg = 0.0, tmp_dgg = 0.0;

  /* loop over all cells */
  for (k=0; k<NCELLS; ++k) {

    int  i, sort;
    cell *p;
    real  tmp;

    p = CELLPTR(k);

    for (i=0; i<p->n; ++i) {

      sort = VSORTE(p,i);
      tmp_gg +=  SPRODN(CG_G,p,i,CG_G,p,i);

      /* which method to construct the conjugated direction */
      if (cg_fr == 1) { /* Fletcher-Reeves */
        tmp_dgg +=  SPRODN(KRAFT,p,i,KRAFT,p,i);
      }
      else { /* Polak-Ribiere */
        tmpvec.x = KRAFT(p,i,X) - CG_G(p,i,X);
        tmpvec.y = KRAFT(p,i,Y) - CG_G(p,i,Y);
#ifndef TWOD
        tmpvec.z = KRAFT(p,i,Z) - CG_G(p,i,Z);
#endif
        tmp_dgg += SPRODX(KRAFT,p,i,tmpvec);
      }
    }
  }

#ifdef MPI
  MPI_Allreduce( &tmp_gg,  &gg,  1, REAL, MPI_SUM, cpugrid);
  MPI_Allreduce( &tmp_dgg, &dgg, 1, REAL, MPI_SUM, cpugrid);
#else
  gg  = tmp_gg;
  dgg = tmp_dgg;
#endif
  cg_gamma = dgg/gg;
  // if ((cg_infolevel>0) && (0==myid)) 
  //{
  //printf("cg_gamma %.12e   gg %.12e\n",cg_gamma,gg);fflush(stdout);
  //}

}

void set_hg(void)
{
  int k;
  real tmp_f_max2 = 0.0;

  /* loop over all cells */
  for (k=0; k<NCELLS; ++k) {

    int  i;
    cell *p;
    p = CELLPTR(k);

    for (i=0; i<p->n; ++i) {

      OLD_ORT(p,i,X) = ORT(p,i,X);
      OLD_ORT(p,i,Y) = ORT(p,i,Y);
#ifndef TWOD
      OLD_ORT(p,i,Z) = ORT(p,i,Z);
#endif
      CG_G(p,i,X) = KRAFT(p,i,X);
      CG_G(p,i,Y) = KRAFT(p,i,Y);
#ifndef TWOD
      CG_G(p,i,Z) = KRAFT(p,i,Z);
#endif
      CG_H(p,i,X) = KRAFT(p,i,X) + cg_gamma * CG_H(p,i,X);
      CG_H(p,i,Y) = KRAFT(p,i,Y) + cg_gamma * CG_H(p,i,Y);
#ifndef TWOD
      CG_H(p,i,Z) = KRAFT(p,i,Z) + cg_gamma * CG_H(p,i,Z);
#endif
      /* determine the biggest force component*/
      tmp_f_max2 = MAX(SQR(KRAFT(p,i,X)),tmp_f_max2);
      tmp_f_max2 = MAX(SQR(KRAFT(p,i,Y)),tmp_f_max2);
#ifndef TWOD
      tmp_f_max2 = MAX(SQR(KRAFT(p,i,Z)),tmp_f_max2);
#endif
    }
  }

#ifdef MPI
  MPI_Allreduce( &tmp_f_max2, &f_max2, 1, REAL, MPI_MAX, cpugrid);
#else
  f_max2 = tmp_f_max2;
#endif

}
