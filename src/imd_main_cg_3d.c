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
* imd_main__cg_3d.c -- main loop for Conjugated gradient, three dimensions
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
*  main_loop
*
*****************************************************************************/

void main_loop(void)
{

  real tmp_pot_energy;

  int ctf=0;
  int astep=0;
  int linminsteps;

#ifdef FBC
  int l;
  vektor nullv={0.0,0.0,0.0};
  vektor temp_df;
  vektor *fbc_df;
  fbc_df = (vektor *) malloc(vtypes*DIM*sizeof(real));
  if (NULL==fbc_df)
    error("Can't allocate memory for fbc_df\n");
#endif


  if (0==myid) printf( "Starting simulation %d\n", simulation );


#ifdef DEFORM
  deform_int = 0;
#endif

/********************** outer simulation loop ****************************/
/* steps now no 'time' steps, one step equals a total relaxation         */

  for (steps=steps_min; steps <= steps_max; ++steps) {


  /******* changing the boundary / initial  conditions after relaxation ***/
 
#ifdef FBC
    temp_df.x = 0.0;
    temp_df.y = 0.0;
    temp_df.z = 0.0;  
    for (l=0; l<vtypes; l++) *(fbc_df+l) = temp_df;

    for (l=0;l<vtypes;l++) {
	(fbc_df+l)->x = (fbc_dforces+l)->x;
	(fbc_df+l)->y = (fbc_dforces+l)->y;
	(fbc_df+l)->z = (fbc_dforces+l)->z;
    } 
#endif /* FBC */

#ifdef HOMDEF
    if ((exp_interval > 0)  expand_sample();
    if ((hom_interval > 0)   shear_sample();
    if ((lindef_interval >0)   lin_deform();
#endif

#ifdef DEFORM
	if(steps>0)     /* first relax undeformed sample */
	deform_sample();
#endif
  /* initialize with random noise ("temperature"), if necessary *
   * -> let the system evolve anealsteps with nve                        */
     if (do_maxwell) maxwell(temperature);
     for(astep=0;astep<annealsteps;astep++)
    {
	calc_forces(astep);
	ctf++;
	move_atoms_nve(); 
	do_boundaries();    
	fix_cells();  
	
    }

    /***** Now do the Conjugate gradient ************/

    /* initialisations: h,g, Fmax */
    calc_forces(steps);
    calc_fnorm_g_h();

    printf(" fmax= %lf \n",sqrt(fmax2));fflush(stdout);
    old_cgval = CGVAL;

    for (cgsteps = 0 ; cgsteps < cg_maxsteps; cgsteps++)
    {
	/*  printf (" cgstep %d, CGVAL = %lf \n",cgsteps,CGVAL);fflush(stdout); */
	/* minimization in one 'direction' */
	
	linminsteps=linmin();
	ctf += linminsteps;
	/*  printf (" cgstep %d, ctf %d linminsteps %d, CGVAL = %lf \n",cgsteps,ctf,linminsteps,CGVAL);fflush(stdout); */
        /* Convergence test: change of Epot or smaller than fnorm*/
#if defined(CGE) || defined(CGEF)
	if (SQR(CGVAL - old_cgval) <= SQR(cg_threshold))  
	    break;
#else
	if (CGVAL <= cg_threshold)
	    break;
#endif
	old_cgval = CGVAL;

	cg_calcgamma();         /* calc gg, dgg */
	
        /* sets ort = old_ort, h, g, needs gamma and gets fmax2*/ 
	set_hg();
	
	do_boundaries();    
	fix_cells();
    
        /* overwrites the checkpoint file after each cgstep */  
	write_config_select(0,"cgchkpt",write_atoms_config,write_header_config);  
	write_eng_file(ctf);
    }
    
    /* write 'relaxed' config */

   
    write_config(steps);


#ifdef USE_SOCKETS
    if ((socket_int>0) && (0==steps%socket_int)) check_socket(steps);
#endif

#ifdef TIMING
    imd_stop_timer(&time_io);
#endif

#ifdef FBC
    /* fbc_forces is already initialised with beginforces */
    for (l=0; l<vtypes; l++){ 
      /* fbc_df=0 if MIK && ekin> ekin_threshold */
      (fbc_forces+l)->x += (fbc_df+l)->x;  
      (fbc_forces+l)->y += (fbc_df+l)->y;
      (fbc_forces+l)->z += (fbc_df+l)->z;
    } 
#endif


  }

  /* clean up the current phase, and clear restart flag */
  restart=0;
  if (0==myid) {
    write_itr_file(-1, steps_max,"");
    printf( "End of simulation %d\n", simulation );
  } 
  
  steps = ctf; /* for easier comparison of the time spent */
}


/******************************************************************************
*
* do_boundaries
*
* Apply periodic boundaries to all atoms
* Could change so that only cells on surface to some work
*
******************************************************************************/

void do_boundaries(void)
{
  int k;

  /* for each cell in bulk */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (k=0; k<ncells; ++k) {

    int l,i;
    cell *p;

    p = cell_array + CELLS(k);

    /* PBC in x direction */
    if (pbc_dirs.x==1)
    for (l=0; l<p->n; ++l) {
      i = -FLOOR(SPRODX(p->ort,l,tbox_x));
      p->ort X(l)    += i * box_x.x;
      p->ort Y(l)    += i * box_x.y;
      p->ort Z(l)    += i * box_x.z;
#ifdef MSQD
      p->refpos X(l) += i * box_x.x;
      p->refpos Y(l) += i * box_x.y;
      p->refpos Z(l) += i * box_x.z;
#endif
#ifdef AVPOS
      p->sheet X(l)  -= i * box_x.x;
      p->sheet Y(l)  -= i * box_x.y;
      p->sheet Z(l)  -= i * box_x.z;
#endif
    }

    /* PBC in y direction */
    if (pbc_dirs.y==1)
    for (l=0; l<p->n; ++l) {
      i = -FLOOR(SPRODX(p->ort,l,tbox_y));
      p->ort X(l)    += i * box_y.x;
      p->ort Y(l)    += i * box_y.y;
      p->ort Z(l)    += i * box_y.z;
#ifdef MSQD
      p->refpos X(l) += i * box_y.x;
      p->refpos Y(l) += i * box_y.y;
      p->refpos Z(l) += i * box_y.z;
#endif
#ifdef AVPOS
      p->sheet X(l)  -= i * box_y.x;
      p->sheet Y(l)  -= i * box_y.y;
      p->sheet Z(l)  -= i * box_y.z;
#endif
    }

    /* PBC in z direction */
    if (pbc_dirs.z==1)
    for (l=0; l<p->n; ++l) {
      i = -FLOOR(SPRODX(p->ort,l,tbox_z));
      p->ort X(l)    += i * box_z.x;
      p->ort Y(l)    += i * box_z.y;
      p->ort Z(l)    += i * box_z.z;
#ifdef MSQD
      p->refpos X(l) += i * box_z.x;
      p->refpos Y(l) += i * box_z.y;
      p->refpos Z(l) += i * box_z.z;
#endif
#ifdef AVPOS
      p->sheet X(l)  -= i * box_z.x;
      p->sheet Y(l)  -= i * box_z.y;
      p->sheet Z(l)  -= i * box_z.z;
#endif
    }

  }
}




/*******************************************************************************
 *
 *                  linmin : linear minimization for cg                        *
 *
 ******************************************************************************/

int linmin()
{
    real alpha_a, alpha_b, alpha_c, fa,fb,fc;
    real alphamin;
    real fmax;
    int iter1, iter2;
    real linmin_dmin = 0.00001;          /* will be later parameter */
/* initialisation */
    alpha_a = 0.0;                     
    iter1=0; 
    iter2=0;

/* can we do the scaling in an other way, additional lowest scale to move ? */
    fmax = sqrt (fmax2);
    /*    alpha_b = linmin_dmax/sqrt(fmax2);  got short distances: new scaling*/  
    
    /* scaling: alpha_b should be 0.01 and alpha_c 0.02 */
    if(fmax > 1.0) 
      {
	alpha_b = 0.01/fmax;
	alpha_c = 0.02/fmax;
      }
    else if(fmax<0.001)
      {
	alpha_b = 0.01 *0.001/fmax;
	alpha_c = 0.02 *0.001/fmax;
      }
    else
      {
	alpha_b = 0.01;
	alpha_c = 0.02;
      }
    fa = old_cgval;
    fb = fonedim(alpha_b);

/*      printf("ID: %d before mnbrak: fmax= %lf alpha_a= %lf alpha_b=%lf fa= %lf fb=%lf \n",myid,fmax,alpha_a,alpha_b,fa,fb); */
/*      fflush(stdout);      */
 
/* decide which method to take to braket a mimimum, at the moment only mbrak, later zbrak? */
    iter1 = mnbrak (&alpha_a,&alpha_b,&alpha_c,&fa,&fb,&fc); /* call by reference Num Rec. p297 */
   /*   if(iter1 <0) */
/*      { */
/*  	printf("error in mnbrak %lf %lf %lf %.12lf %.12lf %.12lf\n",&alpha_a,&alpha_b,&alpha_c,&fa,&fb,&fc);fflush(stdout); */
/*      } */
    /*  printf("after mnbrak \n");fflush(stdout); */
/* decide which method to take to search the mimimum */
#ifdef CGEF /* not implemented yet */
    iter2 = dbrent (alpha_a,alpha_b,alpha_c,fa,fb,fc); /* in ort should be the coord. of min pos. */
#else
    iter2 =  brent (alpha_a,alpha_b,alpha_c,fb,&alphamin);
#endif
  /*    printf("ID: %d  in linmin: iter1= %d iter2 =%d \n",myid,iter1,iter2);fflush(stdout);  */
    return (iter1 + iter2);
}


/* ondimensional function, used by brent */
real fonedim ( real alpha)  /* sets the global variables epot,fnorm corresponding to alpha */
{
    move_atoms_cg(alpha);
    calc_forces(0);       /* why does calc_forces needs steps ? */
    calc_fnorm();
    return (CGVAL);
}



/* calculates the forcenorm
   the forces have to be calculated before!   */
void  calc_fnorm(void)
{
  int k;
  real tmp_fnorm;

  fnorm = 0.0;

  /* loop over all cells */
#ifdef _OPENMP
#pragma omp parallel for reduction(+:fnorm)
#endif
  for (k=0; k<ncells; ++k) {

    int  i, sort;
    cell *p;
    real  tmp;

    p = cell_array + CELLS(k);

/* loop over all particles */
    for (i=0; i<p->n; ++i) {

	sort = VSORTE(p,i);

#ifdef FBC
        /* give virtual particles their extra force */
	p->kraft X(i) += (fbc_forces + sort)->x;
	p->kraft Y(i) += (fbc_forces + sort)->y;
#ifndef TWOD
	p->kraft Z(i) += (fbc_forces + sort)->z;
#endif
#endif

  	/* and set their force (->momentum) in restricted directions to 0 */ 
	p->kraft X(i) *= (restrictions + sort)->x;
	p->kraft Y(i) *= (restrictions + sort)->y;
#ifndef TWOD
	p->kraft Z(i) *= (restrictions + sort)->z;
#endif

	tmp_fnorm +=  SPRODN(p->kraft,i,p->kraft,i);
    }
  }
#ifdef MPI
/* add up results from different CPUs */
  MPI_Allreduce( &tmp_fnorm, &fnorm, 1, MPI_REAL, MPI_SUM, cpugrid);
#else
  fnorm = tmp_fnorm;
#endif

}

void calc_fnorm_g_h(void)
{
  int k;
  real tmp_fnorm;
  real tmp_fmax2;
  
  static int count = 0;
  fnorm = 0.0;
  tmp_fmax2 = 0.0;
  
  /* loop over all cells */
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tot_kin_energy,fnorm)
#endif
  for (k=0; k<ncells; ++k) {

    int  i, sort;
    cell *p;
    real  tmp;

    p = cell_array + CELLS(k);

/* loop over all particles */
    for (i=0; i<p->n; ++i) {

	sort = VSORTE(p,i);
#ifdef FBC
        /* give virtual particles their extra force */
	p->kraft X(i) += (fbc_forces + sort)->x;
	p->kraft Y(i) += (fbc_forces + sort)->y;
#ifndef TWOD
	p->kraft Z(i) += (fbc_forces + sort)->z;
#endif
#endif

  	/* and set their force (->momentum) in restricted directions to 0 */ 
	p->kraft X(i) *= (restrictions + sort)->x;
	p->kraft Y(i) *= (restrictions + sort)->y;
#ifndef TWOD
	p->kraft Z(i) *= (restrictions + sort)->z;
#endif

	tmp_fnorm +=  SPRODN(p->kraft,i,p->kraft,i);
	/* initialise old_ort */
	p->old_ort X(i) = p->ort X(i);
	p->old_ort Y(i) = p->ort Y(i);
	p->old_ort Z(i) = p->ort Z(i);

	/* initialise search vectors */
	p->h X(i) = p->kraft X(i);
	p->h Y(i) = p->kraft Y(i);
#ifndef TWOD
	p->h Z(i) = p->kraft Z(i);
#endif

	p->g X(i) = p->kraft X(i);
	p->g Y(i) = p->kraft Y(i);
#ifndef TWOD
	p->g Z(i) = p->kraft Z(i);
#endif

/* determine the biggest force component */
	tmp_fmax2 = MAX(SQR(p->kraft X(i)),tmp_fmax2);
	tmp_fmax2 = MAX(SQR(p->kraft Y(i)),tmp_fmax2);
	tmp_fmax2 = MAX(SQR(p->kraft Z(i)),tmp_fmax2);
    }
  }
#ifdef MPI
/* find the maximum from all CPUs */
  MPI_Allreduce( &tmp_fmax2, &fmax2, 1, MPI_REAL, MPI_MAX, cpugrid);
/* add up results from different CPUs */
  MPI_Allreduce( &tmp_fnorm, &fnorm, 1, MPI_REAL, MPI_SUM, cpugrid);
#else
  fnorm = tmp_fnorm;
  fmax2 = tmp_fmax2;
#endif

  return;  
}

/* calculation of gamma, dgg,gg */
void cg_calcgamma(void)
{
  int k;
  real tmpvec1[3], tmpvec2[3];
  vektor tmpvec;
  real tmp_fmax2,tmp_gg,tmp_dgg;

  static int count = 0;
  tmp_fmax2 = 0.0;
  tmp_dgg = 0.0;
  tmp_gg = 0.0;

  /* loop over all cells */
  for (k=0; k<ncells; ++k) 
    {

    int  i, sort;
    cell *p;
    real  tmp;

    p = cell_array + CELLS(k);

    for (i=0; i<p->n; ++i) 
      {

	sort = VSORTE(p,i);
	tmp_gg +=  SPRODN(p->g,i,p->g,i);

/* which method to construct the conjugated direction */
#ifdef FR /* Fletcher-Reeves */
	tmp_dgg +=  SPRODN(p->kraft,i,p->kraft,i);
#else  /* Polak-Ribiere */
	tmpvec.x = p->kraft X(i) - p->g X(i);
	tmpvec.y = p->kraft Y(i) - p->g Y(i);
	tmpvec.z = p->kraft Z(i) - p->g Z(i);
	tmp_dgg  += SPRODX(p->kraft,i,tmpvec); 
#endif
/* do we need this here again? is already in set_hg ! */
/* determine the biggest force component */
	/*  tmp_fmax2 = MAX(SQR(p->kraft X(i)),tmp_fmax2); */
/*  	tmp_fmax2 = MAX(SQR(p->kraft Y(i)),tmp_fmax2); */
/*  	tmp_fmax2 = MAX(SQR(p->kraft Z(i)),tmp_fmax2); */
      }
    }

#ifdef MPI
/* do we need this here again? is already in set_hg ! */
/* find the maximum from all CPUs */
  /*  MPI_Allreduce( &tmp_fmax2, &fmax2, 1, MPI_REAL, MPI_MAX, cpugrid); */
/* add up results from different CPUs */
  MPI_Allreduce( &tmp_gg,&gg , 1, REAL, MPI_SUM, cpugrid);
  MPI_Allreduce( &tmp_dgg,&dgg , 1, REAL, MPI_SUM, cpugrid);
#else
  fmax2 = tmp_fmax2;
  gg = tmp_gg;
  dgg = tmp_dgg;
#endif
  cg_gamma = dgg/gg;

}

void set_hg(void)
{
  int k;
  real tmpvec1[3], tmpvec2[3];
  vektor tmpvec;
  real tmp_fmax2,tmp_gg,tmp_dgg;

  static int count = 0;
  tmp_fmax2 = 0.0;
  tmp_dgg = 0.0;
  tmp_gg = 0.0;

  /* loop over all cells */
  for (k=0; k<ncells; ++k) {

    int  i, sort;
    cell *p;
    real  tmp;

    p = cell_array + CELLS(k);

    for (i=0; i<p->n; ++i) {

	p->old_ort X(i) = p->ort X(i);
	p->old_ort Y(i) = p->ort Y(i);
	p->old_ort Z(i) = p->ort Z(i);

	p->g X(i) = p->kraft X(i);
	p->g Y(i) = p->kraft Y(i);
	p->g Z(i) = p->kraft Z(i);

	p->h X(i) = p->kraft X(i) + cg_gamma * p->h X(i);
	p->h Y(i) = p->kraft Y(i) + cg_gamma * p->h Y(i);
	p->h Z(i) = p->kraft Z(i) + cg_gamma * p->h Z(i);

/* determine the biggest force component*/
	tmp_fmax2 = MAX(SQR(p->kraft X(i)),tmp_fmax2);
	tmp_fmax2 = MAX(SQR(p->kraft Y(i)),tmp_fmax2);
	tmp_fmax2 = MAX(SQR(p->kraft Z(i)),tmp_fmax2);
    }
  }

#ifdef MPI
/*    printf("ID: %d  tmp_fmax2: %lf\n",myid,tmp_fmax2);fflush(stdout); */
   MPI_Allreduce( &tmp_fmax2, &fmax2, 1, MPI_REAL, MPI_MAX, cpugrid);
/*    printf("ID: %d  fmax2: %lf\n",myid,fmax2);fflush(stdout);  */
/* add up results from different CPUs */
  MPI_Allreduce( &tmp_gg,&gg , 1, REAL, MPI_SUM, cpugrid);
  MPI_Allreduce( &tmp_dgg,&dgg , 1, REAL, MPI_SUM, cpugrid);
#else
  fmax2 = tmp_fmax2;
  gg = tmp_gg;
  dgg = tmp_dgg;
#endif
  cg_gamma = dgg/gg;

  return;
}




















































