/******************************************************************************
*
* imd_misc.c -- Some Misc. Routines for the imd package
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/* personal debug switch
#define DEBUGINFO 1000 
*/
/******************************************************************************
*
*  Usage -- educate users
*
******************************************************************************/

void usage(void)

{ 

#ifdef MPI
  fprintf(stderr,"%s [-r<nnn>] [-p paramter-file]\n",progname); 
  MPI_Abort(MPI_COMM_WORLD, 1);
#else
  fprintf(stderr,"%s [-r<nnn>] [-p paramter-file]\n",progname); 
#endif
  exit(1); 

}

/******************************************************************************
*
* error -- Complain and abort
*
******************************************************************************/

void error(char *mesg)

{
  
#ifdef MPI
  fprintf(stderr,"CPU %d: %s\n",myid,mesg);
  MPI_Abort(MPI_COMM_WORLD, 1);
#else
  fprintf(stderr,"Error: %s\n",mesg);
#endif
  exit(2);
}


#ifndef MONOLJ

/******************************************************************************
*
* read_potential - reads potential into the potential-array
*
* The format of the potential file is:
* 
* r**2 V00 V01 V02 ... V10 V11 V12 ... VNN
*
* N is the number of different atom types
*
* Note that it is assumed that Vij == Vji and that the r**2 are 
* aequidistant.
*
* Adds: 30/6/98 J. Hahn
*	add function 'read_ttbp_potential '
*	(distance smoothing of three body potential)
*
******************************************************************************/

/* Should somehow handle comments */

/* ------------------------ */
/* classical pair potential */
/* ------------------------ */

void read_potential(str255 potfilename)

{
  real last_r =  0;
  double r;
  FILE *infile;
  int  npot = 0;
  int  nmax;
  int  i,j,k;
  double val;
  real delta;
  int tablesize;
  int c;

#ifdef MPI
  if (0 == myid) { /* Read Potential only on master processor */
#endif

  infile = fopen(potfilename,"r");
  if (NULL==infile) error("Can't open potential file.");

#ifndef STATIC_POT
  tablesize = PSTEP*ntypes*ntypes*sizeof(real);
  potential = (real *) malloc(tablesize);
  if (NULL==potential) error("Can't allocate memory for potential.");
  nmax = PSTEP;
#endif
 
  r2_0   = -1;
  r2_end = -1;

  /* input loop */
  while (!feof(infile)) {

#ifndef STATIC_POT
    /* still some space left? */ 
    if (((npot%PSTEP) == 0) && (npot>0)) {
      tablesize = (npot+PSTEP)*ntypes*ntypes*sizeof(real);
      potential = (real *)
	realloc(potential,tablesize);
      if (NULL==potential) error("Can't extend memory for potential.");
      nmax = npot+PSTEP;
    };
#endif
    /*  read in potential */
    if ( 1 != fscanf(infile,"%lf",&r) ) break;
    
    for (i=0;i<ntypes;++i)
      for (j=0;j<ntypes;++j) {
	if ( 1 != fscanf(infile,"%lf", &val)) 
           error("Line incomplete in potential file.");
#ifdef STATIC_POT
	potential[i][j][npot] = val;
#else
	*PTR_3D(potential,npot,i,j,nmax,ntypes,ntypes) = val;
#endif
      };
    /* Catch first value */
    if (r2_0  == -1) r2_0   = r;
    if (r2_end== -1) r2_end = r;
    
    /* Get minimum and maximum  */
    r2_0   = MIN(r2_0  ,r);
    r2_end = MAX(r2_end,r);
    
    
    /* Check if the r really go up */
    if (r<last_r) error("Value of r is not incremental.");
    
    last_r = r;
    ++npot;
    
  };

  fclose(infile);

  r2_step = (r2_end - r2_0) / npot;

  printf("Read potential %s with %d lines.\n",potfilename,npot);
  printf("Starts at r2_0 : %f, r_0  : %f\n",r2_0,sqrt(r2_0));
  printf("Ends at r2_end : %f, r_end: %f\n",r2_end,sqrt(r2_end));
  printf("Step is r2_step: %f\n",r2_step);

  /* Shift potential to zero */

  for (i=0; i<ntypes; ++i)
    for (j=0; j<ntypes; ++j) {
#ifdef STATIC_POT
      delta = potential[i][j][(npot-1)];
#else
      delta = *PTR_3D(potential,(npot-1),i,j,nmax,ntypes,ntypes);
#endif
      printf("Potential %1d%1d shifted by %f\n",i,j,delta);
      for (k=0; k<npot; ++k) 
#ifdef STATIC_POT
        potential[i][j][k] -=delta;
#else
	*PTR_3D(potential,k,i,j,nmax,ntypes,ntypes) -= delta;
#endif
    };

  printf("\n");

  /* The Interpolation uses k+1 and k+2, so we add zeros at end of table */

  for(k=1; k<=5; ++k){
#ifndef STATIC_POT
    /* still some space left? */ 
    if (((npot%PSTEP) == 0) && (npot>0)) {
      tablesize = (npot+PSTEP)*ntypes*ntypes*sizeof(real);
      potential = (real *)
	realloc(potential,tablesize);
      if (NULL==potential) error("Can't extend memory for potential.");
      nmax = npot+PSTEP;
    };
#endif
    for (i=0;i<ntypes;++i)
      for (j=0;j<ntypes;++j)	
#ifdef STATIC_POT
	potential[i][j][npot] = 0.0;
#else
	*PTR_3D(potential,npot,i,j,nmax,ntypes,ntypes) = 0.0;
#endif
    ++npot;
  };

#ifdef STATIC_POT
  if ((MAXPOTLEN<npot) || (MAXATOMTYPES<ntypes)) 
    error("Messed up potenital table.");
#endif


#ifdef MPI
  };

  /* Broadcast potential table to clients */

  MPI_Bcast( &nmax,      1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &tablesize, 1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &r2_step  , 1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &r2_end   , 1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &r2_0     , 1, MPI_REAL, 0, MPI_COMM_WORLD);

#ifdef STATIC_POT

  MPI_Bcast( potential, MAXPOTLEN*MAXATOMTYPES*MAXATOMTYPES
	    , MPI_REAL, 0, MPI_COMM_WORLD);

#else
  if (0!=myid) potential = (real *) malloc(tablesize);
  if (NULL==potential) error("Can't allocate memory for potential on client.");

  MPI_Bcast( potential, tablesize / sizeof(real), MPI_REAL, 0, MPI_COMM_WORLD);
#endif
#endif

  /* Set cutoff */
  r2_cut = r2_end + r2_step;

  /* Save dimensions of potential array */
#ifdef STATIC_POT
  pot_dim.x = MAXPOTLEN;
  pot_dim.y = MAXATOMTYPES;
  pot_dim.z = MAXATOMTYPES;
#else
  pot_dim.x = nmax;
  pot_dim.y = ntypes;
  pot_dim.z = ntypes;
#endif
  
  inv_r2_step = 1 / r2_step;
  
} /* end of pair potential */



#ifdef EAM2
/* input of the tabulated function values for... */

void eam2_read_core_pot(str255 core_pot_filename)

     /* numstep is the number of datapoints!!!
	required layout of the Core-Core-Potential file:
	really ugly, but easy to code and allows different
	stepsizes etc. can be improved latter.
	requires the values to belong to aequidistant increasing r^2.

	at the beginning ntypes*ntypes times a line of
	        r_begin r_end r_step
	then comes the values of the Core Potential:
	in form of one(!) column, first for atom pair  00, then 01 and so on
	after that a linefeed would be nice (allows easier splitting of files
	for gnuplot), and then the next column: 1 and so on (ntypes times) 
     */

{
  FILE *infile;
  int eam2_phi_tablesize,eam2_phi_r_infosize;
  int i,j,k;
  real r_begin,r_end,r_step,val=0;
  real numstep;
  int number_of_steps;

#ifdef MPI
  if (0 == myid) { /* Read Potential only on master processor */
#endif

  infile = fopen(core_pot_filename,"r");
  if (NULL==infile) error("Can't open core_potential_file.");

  /* reading the information about the function table */
  eam2_phi_r_infosize=ntypes*ntypes*sizeof(real);
  eam2_phi_r_begin=(real *) malloc(eam2_phi_r_infosize);
  if (NULL==eam2_phi_r_begin) 
    error("Can't allocate memory for eam2_phi_r_begin.");
  eam2_phi_r_end=(real *) malloc(eam2_phi_r_infosize);
  if (NULL==eam2_phi_r_end) 
    error("Can't allocate memory for eam2_phi_r_end.");
  eam2_phi_r_step=(real *) malloc(eam2_phi_r_infosize);
  if (NULL==eam2_phi_r_step) 
    error("Can't allocate memory for eam2_phi_r_step.");

  for(i=0;i<ntypes;i++){
    for(j=0;j<ntypes;j++){
      if ( 1 != fscanf(infile,"%lf", &r_begin)) 
        error("Info line in core_potential_file corrupt.");
      if ( 1 != fscanf(infile,"%lf", &r_end)) 
        error("Info line in core_potential_file corrupt.");
      if ( 1 != fscanf(infile,"%lf", &r_step)) 
        error("Info line in core_potential_file corrupt.");
      numstep=1+(r_end - r_begin)/r_step;
      number_of_steps = (int) (numstep+0.5);
      *PTR_2D(eam2_phi_r_begin,i,j,ntypes,ntypes) =r_begin ;
      *PTR_2D(eam2_phi_r_end,i,j,ntypes,ntypes) =r_end ;
      *PTR_2D(eam2_phi_r_step,i,j,ntypes,ntypes) =r_step ;
      if (number_of_steps>eam2_max_phi_r_steps) {
        eam2_max_phi_r_steps=number_of_steps; 
      }
      /* beware of stuff like numstep = 499.999999 */
      if (fabs(number_of_steps - numstep) >= 0.1) {
        fprintf(stderr,"WARNING: core_potential_file:\
                       numstep: %lf rounded to:%lf \n",
                       numstep, number_of_steps);
      }
    }
  }
  
  /* reading the functiontables */
  eam2_phi_tablesize = eam2_max_phi_r_steps*ntypes*ntypes*sizeof(real);
  eam2_phi = (real *) malloc(eam2_phi_tablesize);
  if (NULL==eam2_phi) error("Can't allocate memory for eam2_phi.");
  
  /* input loop */
  for (i=0;i<ntypes;++i)
    for (j=0;j<ntypes;++j) {
      r_end  =*PTR_2D(eam2_phi_r_end,i,j,ntypes,ntypes);
      r_begin=*PTR_2D(eam2_phi_r_begin,i,j,ntypes,ntypes);
      r_step =*PTR_2D(eam2_phi_r_step,i,j,ntypes,ntypes);
      numstep=1+(r_end - r_begin)/r_step;
      number_of_steps = (int)(numstep+0.5);
      for(k=0;k<number_of_steps;k++){
	if ( 1 != fscanf(infile,"%lf", &val)) 
	    error("wrong format of core_potential_file.");
	*PTR_3D(eam2_phi,k,i,j,eam2_max_phi_r_steps,ntypes,ntypes) = val;
#ifdef DEBUGINFO
	printf("\n Phi_r: i=%d j=%d k=%d val=%lf",i,j,k,val);
#endif
	}
    };

  fclose(infile);
  printf("Read tabulated Core Potential Function %s for %d atoms.\n",
         core_pot_filename ,ntypes);
  printf("eam2_max_phi_r_steps: %d\n",eam2_max_phi_r_steps);
#ifdef MPI
  };
  /* Broadcast table and table information to clients */
  MPI_Bcast( &eam2_max_phi_r_steps,      1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &eam2_phi_tablesize, 1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &eam2_phi_r_infosize, 1, MPI_INT,  0, MPI_COMM_WORLD);
  if (0!=myid) eam2_phi  = (real *) malloc(eam2_phi_tablesize);
  if (NULL==eam2_phi) 
    error("Can't allocate memory for eam2_phi table on client.");
  if (0!=myid) eam2_phi_r_begin= (real *) malloc(eam2_phi_r_infosize);
  if (NULL==eam2_phi_r_begin) 
    error("Can't allocate memory for eam2_phi_r_begin on client.");
  if (0!=myid) eam2_phi_r_end = (real *) malloc(eam2_phi_r_infosize);
  if (NULL==eam2_phi_r_end) 
    error("Can't allocate memory for eam2_phi_r_end on client.");
  if (0!=myid) eam2_phi_r_step = (real *) malloc(eam2_phi_r_infosize);
  if (NULL==eam2_phi_r_step) 
    error("Can't allocate memory for eam2_phi_r_step on client.");
  MPI_Bcast( eam2_phi, eam2_phi_tablesize / sizeof(real), 
             MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_phi_r_begin, eam2_phi_r_infosize / sizeof(real), 
             MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_phi_r_end, eam2_phi_r_infosize / sizeof(real), 
             MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_phi_r_step, eam2_phi_r_infosize / sizeof(real), 
             MPI_REAL, 0, MPI_COMM_WORLD);
#endif
} /* end of eam2_read_core_pot */


void eam2_read_embedding_energy(str255 emb_E_filename)

     /* numstep = number of datapoints!!!
	required layout of the embedding energy file:
	really ugly, but easy to code and allows different
	stepsizes etc. can be improved latter.
	requires the values to belong to aequidistant increasing rho_h.

	at the beginning ntypes times a line of
	        rho_begin rho_end rho_step
	then comes the values of the embedding energy:
	in form of one(!) column, first for taom type 0
	after that a linefeed would be nice (allows easier splitting of files
	for gnuplot), and then the next column: 1 and so on (ntypes times) 
     */
{
  FILE *infile;
  int eam2_f_i_tablesize,eam2_rho_infosize;
  int i,k;
  real rho_begin,rho_end,rho_step,val=0;
  real numstep;
  int number_of_steps;

#ifdef MPI
  if (0 == myid) { /* Read Potential only on master processor */
#endif

  infile = fopen(emb_E_filename,"r");
  if (NULL==infile) error("Can't open embedding_energy_file.");

  /* reading the information about the function table */
  eam2_rho_infosize=ntypes*sizeof(real);
  eam2_rho_begin=(real *) malloc(eam2_rho_infosize);
  if (NULL==eam2_rho_begin) error("Can't allocate memory for eam2_rho_begin.");
  eam2_rho_end=(real *) malloc(eam2_rho_infosize);
  if (NULL==eam2_rho_end) error("Can't allocate memory for eam2_rho_end.");
  eam2_rho_step=(real *) malloc(eam2_rho_infosize);
  if (NULL==eam2_rho_step) error("Can't allocate memory for eam2_rho_step.");

  for(i=0;i<ntypes;i++){
    if ( 1 != fscanf(infile,"%lf", &rho_begin))  
      error("Info line in embedding_energy_file corrupt.");
    if ( 1 != fscanf(infile,"%lf",&rho_end)) 
      error("Info line in embedding_energy_file corrupt.");
    if ( 1 != fscanf(infile,"%lf",&rho_step)) 
      error("Info line in embedding_energy_file corrupt.");

    numstep=1+(rho_end - rho_begin)/rho_step;
    number_of_steps = (int)(numstep+0.5);  

    *(eam2_rho_begin+i) = rho_begin;
    *(eam2_rho_end+i)   = rho_end; 
    *(eam2_rho_step+i)  = rho_step; 
    if (number_of_steps>eam2_max_rho_steps) { 
      eam2_max_rho_steps=number_of_steps;
    }

    /* beware of stuff like numstep = 499.999999 */
    if (fabs(number_of_steps - numstep) >= 0.1) {
        fprintf(stderr,"WARNING: embedding_energy_file:\
                       numstep: %lf rounded to:%lf \n",
                       numstep, number_of_steps);
      }
  }
  
  /* reading the functiontables */
  eam2_f_i_tablesize = eam2_max_rho_steps*ntypes*sizeof(real);
  eam2_f_i = (real *) malloc(eam2_f_i_tablesize);
  if (NULL==eam2_f_i) error("Can't allocate memory for eam2_f_i .");
  
  /* input loop */
  for (i=0;i<ntypes;i++){
      rho_end  =*(eam2_rho_end+i);
      rho_begin=*(eam2_rho_begin+i);
      rho_step =*(eam2_rho_step+i);
      numstep=1+(rho_end - rho_begin)/rho_step;
      number_of_steps = (int)(numstep+0.5);
      for(k=0;k<number_of_steps;k++){
	if ( 1 != fscanf(infile,"%lf", &val)) 
	  error("wrong format of embedding_energy_file.");
	*PTR_2D(eam2_f_i,k,i,eam2_max_rho_steps,ntypes) = val;
#ifdef DEBUGINFO
	printf("F_rho: i=%d  k=%d val=%lf\n",i,k,val);
#endif
      }
    };

  fclose(infile);
  printf("Read tabulated Embedding Energy Function %s for %d atoms.\n",
         emb_E_filename,ntypes);
  printf("eam2_max_rho_steps = %d\n",eam2_max_rho_steps);

#ifdef MPI
  };
  /* Broadcast table and table information to clients */
  MPI_Bcast( &eam2_max_rho_steps,      1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &eam2_f_i_tablesize, 1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &eam2_rho_infosize, 1, MPI_INT,  0, MPI_COMM_WORLD);
  if (0!=myid) eam2_f_i  = (real *) malloc(eam2_f_i_tablesize);
  if (NULL==eam2_f_i) 
    error("Can't allocate memory for embedding energy table on client.");
  if (0!=myid) eam2_rho_begin= (real *) malloc(eam2_rho_infosize);
  if (NULL==eam2_rho_begin) 
    error("Can't allocate memory for eam2_rho_begin on client.");
  if (0!=myid) eam2_rho_end = (real *) malloc(eam2_rho_infosize);
  if (NULL==eam2_rho_end) 
    error("Can't allocate memory for eam2_rho_end on client.");
  if (0!=myid) eam2_rho_step = (real *) malloc(eam2_rho_infosize);
  if (NULL==eam2_rho_step) 
    error("Can't allocate memory for eam2_rho_step on client.");
  MPI_Bcast( eam2_f_i, eam2_f_i_tablesize / sizeof(real), 
             MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_rho_begin, eam2_rho_infosize / sizeof(real), 
             MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_rho_end, eam2_rho_infosize / sizeof(real), 
             MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_rho_step, eam2_rho_infosize / sizeof(real), 
             MPI_REAL, 0, MPI_COMM_WORLD);
#endif
} /* end of eam2_readenergy */


void eam2_read_atomic_rho(str255 at_rho_filename)

     /* numstep =number of datapoints!!!
	required layout of the atomic e-density file:
	really ugly, but easy to code and allows different
	stepsizes etc. can be improved latter.
	requires the values to belong to aequidistant increasing r^2.

	at the beginning ntypes*ntypes times a line of
	        r_begin r_end r_step
	then comes the values of the electron desity at atomtype i from type j:
	in form of one(!) column, first for 00 after that a linefeed 
        would be nice (allows easier splitting of files for gnuplot), 
        and then the next column: 01 and so on (ntypes*ntypes times) 
     */
{
  FILE *infile;
  int eam2_rho_at_tablesize,eam2_r_infosize;
  int i,j,k;
  real r_begin,r_end,r_step,val=0;
  real numstep;
  int number_of_steps;

#ifdef MPI
  if (0 == myid) { /* Read Potential only on master processor */
#endif

  infile = fopen(at_rho_filename,"r");
  if (NULL==infile) error("Can't open atomic_e-density_file.");

  /* reading the information about the function table */
  eam2_r_infosize=ntypes*ntypes*sizeof(real);
  eam2_r_begin=(real *) malloc(eam2_r_infosize);
  if (NULL==eam2_r_begin) error("Can't allocate memory for eam2_r_begin.");
  eam2_r_end=(real *) malloc(eam2_r_infosize);
  if (NULL==eam2_r_end) error("Can't allocate memory for eam2_r_end.");
  eam2_r_step=(real *) malloc(eam2_r_infosize);
  if (NULL==eam2_r_step) error("Can't allocate memory for eam2_r_step.");

  for (i=0;i<ntypes;i++) {
    for (j=0;j<ntypes;j++) {
      if ( 1 != fscanf(infile,"%lf", &r_begin)) 
        error("Info line in atomic_e-density_file corrupt.");
      if ( 1 != fscanf(infile,"%lf", &r_end)) 
        error("Info line in atomic_e-density_file corrupt.");
      if ( 1 != fscanf(infile,"%lf", &r_step)) 
        error("Info line in atomic_e-density_file corrupt.");
      numstep=1+(r_end - r_begin)/r_step;
      number_of_steps = (int) (numstep+0.5);
      *PTR_2D(eam2_r_begin,i,j,ntypes,ntypes) =r_begin ;
      *PTR_2D(eam2_r_end,i,j,ntypes,ntypes) =r_end ;
      *PTR_2D(eam2_r_step,i,j,ntypes,ntypes) =r_step ;
      if(number_of_steps>eam2_max_r_steps) {
        eam2_max_r_steps=number_of_steps; 
      }

      /* what should i do: numstep=499.99999932*/
      if (fabs(number_of_steps - numstep) >= 0.1) {
        fprintf(stderr,"WARNING: atomic_e-density_file:\
                       numstep: %lf rounded to:%lf \n",
		numstep, number_of_steps);
      }
    }
  }

  /* reading the functiontables */
  eam2_rho_at_tablesize = eam2_max_r_steps*ntypes*ntypes*sizeof(real);
  eam2_rho_at = (real *) malloc(eam2_rho_at_tablesize);
  if (NULL==eam2_rho_at) error("Can't allocate memory for  eam2_rho_at.");
  
  /* input loop */
  for (i=0;i<ntypes;++i)
    for (j=0;j<ntypes;++j) {
      r_end  =*PTR_2D(eam2_r_end,i,j,ntypes,ntypes);
      r_begin=*PTR_2D(eam2_r_begin,i,j,ntypes,ntypes);
      r_step =*PTR_2D(eam2_r_step,i,j,ntypes,ntypes);
      numstep=1+(r_end - r_begin)/r_step;
      number_of_steps = (int) (numstep+0.5);
      for(k=0;k<number_of_steps;k++){
	if ( 1 != fscanf(infile,"%lf", &val)) 
	    error("wrong format of atomic_e-density_file.");
	*PTR_3D(eam2_rho_at,k,i,j,eam2_max_r_steps,ntypes,ntypes) = val;
#ifdef DEBUGINFO
	printf("Rho_r: i=%d j=%d k=%d val=%lf\n",i,j,k,val);
#endif
	}
    };

  fclose(infile);
  printf("Read tabulated atomic e-density Function %s for %d atoms.\n",at_rho_filename ,ntypes);
  printf("eam2_max_r_steps: %d\n",eam2_max_r_steps);
#ifdef MPI
  };
  /* Broadcast table and table information to clients */
  MPI_Bcast( &eam2_max_r_steps,      1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &eam2_rho_at_tablesize, 1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &eam2_r_infosize, 1, MPI_INT,  0, MPI_COMM_WORLD);
  if (0!=myid) eam2_rho_at  = (real *) malloc(eam2_rho_at_tablesize);
  if (NULL==eam2_rho_at) 
    error("Can't allocate memory for eam2_rho_at table on client.");
  if (0!=myid) eam2_r_begin= (real *) malloc(eam2_r_infosize);
  if (NULL==eam2_r_begin) 
    error("Can't allocate memory for eam2_r_begin on client.");
  if (0!=myid) eam2_r_end = (real *) malloc(eam2_r_infosize);
  if (NULL==eam2_r_end) 
    error("Can't allocate memory for eam2_r_end on client.");
  if (0!=myid) eam2_r_step = (real *) malloc(eam2_r_infosize);
  if (NULL==eam2_r_step) 
    error("Can't allocate memory for eam2_r_step on client.");
  MPI_Bcast( eam2_rho_at, eam2_rho_at_tablesize / sizeof(real), 
             MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_r_begin, eam2_r_infosize / sizeof(real), 
             MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_r_end, eam2_r_infosize / sizeof(real), 
             MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_r_step, eam2_r_infosize / sizeof(real), 
             MPI_REAL, 0, MPI_COMM_WORLD);
#endif
} /* end of eam2_read_atomic_rho */

#endif /* EAM2 function tables */


#ifdef TTBP

/* ----------------------------------------- */
/* smoothing of three body potential in TTBP */
/* ----------------------------------------- */

void read_ttbp_potential(str255 ttbp_potfilename)

{

  real last_r = 0.0;
  double r;
  FILE *infile;
  int  npot = 0;
  int  ttbp_nmax;
  int  i,j,k;
  double val;
  real delta;
  int ttbp_tablesize;
  int c;

#ifdef MPI
  if (0 == myid) { /* Read Potential only on master processor */
#endif

  infile = fopen(ttbp_potfilename,"r");
  if (NULL==infile) error("Can't open potential file.");

#ifndef STATIC_POT
  ttbp_tablesize = PSTEP*ntypes*ntypes*sizeof(real);
  ttbp_potential = (real *) malloc(ttbp_tablesize);
  if (NULL==ttbp_potential) error("Can't allocate memory for TTBP potential.");
  ttbp_nmax = PSTEP;
#endif
 
  ttbp_r2_0   = -1;
  ttbp_r2_end = -1;

  /* zero ttbp cut */
  for (i=0;i<ntypes;++i) 
    for (j=0;j<ntypes;++j) {
      ttbp_r2_cut[i][j] = -1;
  }

  /* input loop */
  while (!feof(infile)) {

#ifndef STATIC_POT
    /* still some space left? */ 
    if (((npot%PSTEP) == 0) && (npot>0)) {
      ttbp_tablesize = (npot+PSTEP)*ntypes*ntypes*sizeof(real);
      ttbp_potential = (real *)	realloc(ttbp_potential,ttbp_tablesize);
      if (NULL==ttbp_potential) 
        error("Cannot extend memory for TTBP potential.");
      ttbp_nmax = npot+PSTEP;
    }
#endif
    /*  read in potential */
    if ( 1 != fscanf(infile,"%lf",&r) ) break;
    
    for (i=0;i<ntypes;++i) 
      for (j=0;j<ntypes;++j) {
	if ( 1 != fscanf(infile,"%lf", &val)) 
          error("Line incomplete in potential file.");
#ifdef STATIC_POT
	ttbp_potential[i][j][npot] = val;
#else
	*PTR_3D(ttbp_potential,npot,i,j,ttbp_nmax,ntypes,ntypes) = val;
#endif
        if ((ttbp_r2_cut[i][j]==-1) && (val==0)) ttbp_r2_cut[i][j] = r;
      };
    /* Catch first value */
    if (ttbp_r2_0  == -1) ttbp_r2_0   = r;
    if (ttbp_r2_end== -1) ttbp_r2_end = r;
    
    /* Get minimum and maximum  */
    ttbp_r2_0   = MIN(ttbp_r2_0  ,r);
    ttbp_r2_end = MAX(ttbp_r2_end,r);

    /* Check if the r really go up */
    if (r<last_r) error("TTBP: Value of r is not incremental.");
    
    last_r = r;
    ++npot;
    
  };

  fclose(infile);

  ttbp_r2_step = (ttbp_r2_end - ttbp_r2_0) / npot;

  printf("TTBP Read potential %s with %d lines.\n",ttbp_potfilename,npot);
  printf("TTBP Starts at r2_0 : %f, r_0  : %f\n",ttbp_r2_0,sqrt(ttbp_r2_0));
  printf("TTBP Ends at r2_end : %f, r_end: %f\n",ttbp_r2_end,sqrt(ttbp_r2_end));
  printf("TTBP Step is r2_step: %f\n",ttbp_r2_step);

    for (i=0;i<ntypes;++i) 
      for (j=0;j<ntypes;++j) {
        printf("TTBP Smooth potential cutoff for %d and %d is: %f\n", 
               i,j,sqrt(ttbp_r2_cut[i][j]));
    }

  /* Shift potential to zero */

  for (i=0; i<ntypes; ++i)
    for (j=0; j<ntypes; ++j) {
#ifdef STATIC_POT
      delta = ttbp_potential[i][j][(npot-1)];
#else
      delta = *PTR_3D(ttbp_potential,(npot-1),i,j,ttbp_nmax,ntypes,ntypes);
#endif
      printf("TTBP Potential %1d%1d shifted by %f\n",i,j,delta);
      for (k=0; k<npot; ++k) 
#ifdef STATIC_POT
        ttbp_potential[i][j][k] -=delta;
#else
	*PTR_3D(ttbp_potential,k,i,j,ttbp_nmax,ntypes,ntypes) -= delta;
#endif
    };

  printf("\n");

  /* The Interpolation uses k+1 and k+2, so we add zeros at end of table */

  for (k=1; k<=5; ++k){
#ifndef STATIC_POT
    /* still some space left? */ 
    if (((npot%PSTEP) == 0) && (npot>0)) {
      ttbp_tablesize = (npot+PSTEP)*ntypes*ntypes*sizeof(real);
      ttbp_potential = (real *)
	realloc(ttbp_potential,ttbp_tablesize);
      if (NULL==ttbp_potential) 
        error("Can't extend memory for ttbp potential.");
      ttbp_nmax = npot+PSTEP;
    }
#endif
    for (i=0;i<ntypes;++i)
      for (j=0;j<ntypes;++j)	
#ifdef STATIC_POT
	ttbp_potential[i][j][npot] = 0.0;
#else
	*PTR_3D(ttbp_potential,npot,i,j,ttbp_nmax,ntypes,ntypes) = 0.0;
#endif
    ++npot;
  };

#ifdef STATIC_POT
  if ((MAXPOTLEN<npot) || (MAXATOMTYPES<ntypes)) 
    error("Messed up potenital table.");
#endif


#ifdef MPI
  }

  /* Broadcast potential table to clients */

  MPI_Bcast( &ttbp_nmax     , 1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &ttbp_tablesize, 1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &ttbp_r2_step  , 1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ttbp_r2_end   , 1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ttbp_r2_0     , 1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ttbp_r2_cut,   100, MPI_REAL, 0, MPI_COMM_WORLD);

#ifdef STATIC_POT

  MPI_Bcast( ttbp_potential, MAXPOTLEN*MAXATOMTYPES*MAXATOMTYPES,
	     MPI_REAL, 0, MPI_COMM_WORLD);

#else
  if (0!=myid) ttbp_potential = (real *) malloc(ttbp_tablesize);
  if (NULL==ttbp_potential) 
    error("Can't allocate memory for ttbp_potential on client.");

  MPI_Bcast( ttbp_potential, ttbp_tablesize / sizeof(real), 
             MPI_REAL, 0, MPI_COMM_WORLD);
#endif
#endif /* MPI */

  /* Set cutoff */
  for (i=0;i<ntypes;++i) 
    for (j=0;j<ntypes;++j) {
      ttbp_r2_cut[i][j] = ttbp_r2_cut[i][j] + ttbp_r2_step;
    };

  /* Save dimensions of potential array */
#ifdef STATIC_POT
  ttbp_pot_dim.x = MAXPOTLEN;
  ttbp_pot_dim.y = MAXATOMTYPES;
  ttbp_pot_dim.z = MAXATOMTYPES;
#else
  ttbp_pot_dim.x = ttbp_nmax;
  ttbp_pot_dim.y = ntypes;
  ttbp_pot_dim.z = ntypes;
#endif
  
  ttbp_inv_r2_step = 1 / ttbp_r2_step;
  
} /* end of smooth TTBP potential */

#endif /* TTBP */

#endif /* not MONOLJ */


/******************************************************************************
*
* generate_atoms - generate atoms for initial configuration
*
* filenames starting with a dot don't specify a file to read from
* but a crystal structure to generate as an initial configuration:
*
* .fcc      -- generates fcc structure
* .nacl     -- generates binary nacl structure (atom type 0 and 1)
* .tiqc     -- generates a truncated icosahedra quasicrystal
* .hex      -- generates 2D hexagonal crystal
* .lav      -- generates a cubic Laves structure A15 (MgCu2)
*
* The lattice constant of the crystal structures (fcc and nacl) is 2.0.
*
******************************************************************************/

void generate_atoms(str255 mode)
{
  cell *to;
  int addnumber, tmp, ninc;
  int r, s, t, i, k;
#ifdef MPI
  MPI_Status status;
#endif

  if ((num_sort=calloc(ntypes,sizeof(int)))==NULL) {
      error("cannot allocate memory for num_sort\n");
  };

#ifdef DISLOC
  calc_Epot_ref=1;
#endif

  /* We always have to initialize the velocities on generated structures */
  do_maxwell=1;

#ifdef TWOD
  if (0 == strcmp(mode,".hex")) {          /* hexagonal crystal */
    init_hex();
    init_cells();
    generate_hex();
#else /* 3D */
  if (0 == strcmp(mode,".fcc")) {          /* fcc */
    init_cubic();
    init_cells();
    generate_fcc(0);
  } else if (0 == strcmp(mode,".nacl")) {  /* NaCl */
    init_cubic();
    init_cells();
    generate_fcc(1);
  } else if (0 == strcmp(mode,".lav")) {   /* Laves */
    init_cubic();
    init_cells();
    generate_lav();
#ifdef QUASI
  } else if (0 == strcmp(mode,".tiqc")) {  /* quasicrystal */
    init_qc();
    init_cells();
    generate_qc();
#endif
#endif /* 3D */
  } else error("Filename with dot specifies unknown structure.");

#ifdef MPI
#ifndef MONOLJ
  /* Get numbering of atoms right even across CPUs */
  ninc = 0;
  if (0==myid)
     MPI_Send( &natoms, 1, MPI_INT, myid + 1, 0, cpugrid );
  else if ((0<myid) && (myid<(num_cpus-1))) {
     MPI_Recv( &tmp,    1, MPI_INT, myid - 1, 0, cpugrid, &status );
     ninc = tmp;
     tmp += natoms;
     MPI_Send( &tmp,    1, MPI_INT, myid + 1, 0, cpugrid );

  }
  else if (myid==(num_cpus-1)) {
     MPI_Recv( &tmp,    1, MPI_INT, myid - 1, 0, cpugrid, &status );
     ninc = tmp;
  }
  /* loop over all atoms, fix numbers */
  for (k=0; k<ncells; ++k) {
    int i;
    cell *p;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) p->nummer[i] += ninc;
  }
#endif /* MONOLJ */

  MPI_Allreduce( &natoms,  &addnumber, 1, MPI_INT, MPI_SUM, cpugrid);
  natoms = addnumber;
  MPI_Allreduce( &nactive, &addnumber, 1, MPI_INT, MPI_SUM, cpugrid);
  nactive = addnumber;
  for (i=0; i<ntypes; i++) {
    MPI_Allreduce(&num_sort[i], &addnumber, 1, MPI_INT, MPI_SUM, cpugrid);
    num_sort[i] = addnumber;
  }
#endif /* MPI */

  if (0== myid) {
    printf("Generated %s structure with %d atoms.\n",mode,natoms);
    addnumber=num_sort[0];
    printf("num_sort = [ %u",num_sort[0]);
    for (i=1; i<ntypes; i++) {
      printf(", %u",num_sort[i]);
      addnumber+=num_sort[i];
    };
    printf(" ],  total = %u\n",addnumber);
  }
  return;
}


#ifdef TWOD

/* initialize for hexagonal crystal */
void init_hex(void)
{
  if (box_param.x != 0) {
    box_x.x = box_param.x * sqrt(3.0);
    box_x.y = 0.0;
    box_y.x = 0.0;
    box_y.y = box_param.y;
  } else { /* backward compatibility */
    box_param.x = (int) box_x.x;
    box_param.y = (int) box_y.y;
    box_x.x *= sqrt(3.0);
  }
}

/* generate hexagonal crystal */
void generate_hex()
{
  cell    *input;
  vektor  min, max;
  ivektor cellc;
  int     to_cpu;
  int     i, j, typ, sign;
  real    x, y;

#ifdef MPI
  min.x =  my_coord.x      * 2 * box_param.x / cpu_dim.x;
  max.x = (my_coord.x + 1) * 2 * box_param.x / cpu_dim.x;
  min.y =  my_coord.y      * 2 * box_param.y / cpu_dim.y;
  max.y = (my_coord.y + 1) * 2 * box_param.y / cpu_dim.y;
#else
  min.x = 0; max.x = 2 * box_param.x;
  min.y = 0; max.y = 2 * box_param.y;
#endif

  /* Set up 1 atom input cell */
  input = (cell *) malloc(sizeof(cell));
  if (0==input) error("Can't allocate input cell.");
  input->n_max=0;
  alloc_cell(input, 1);
  input->masse[0] = 1.0;

  natoms = 0;
  nactive = 0;

  for (i=min.x ; i<max.x; i++)
    for (j=min.y; j<max.y; j++) {

      sign = 1;
      typ  = (i+j) % 2;
      if (typ > 0) continue;

      x = (i+0.5) * sqrt(3.0) * 0.5;
      y = (j+0.5) * 0.5;

#if defined(FRAC) || defined(DEFORM)
      /* leave boundary open if necessary */
      if ((x < strip_width/2) || (x > box_x.x - strip_width/2) ||
          (y < strip_width/2) || (y > box_y.y - strip_width/2)) continue;
      /* boundary atoms in x-direction get negative numbers */
      if ((x < strip_width)   || (x > box_x.x - strip_width)) sign = -1;
#endif

      natoms++;
      if (sign>0) nactive++;

      input->n = 1;
      input->ort X(0) = x;
      input->ort Y(0) = y;
      input->nummer[0] = sign * natoms;
      input->sorte[0] = typ;
      num_sort[typ]++;
      cellc = cell_coord(x,y);
#ifdef MPI
      to_cpu = cpu_coord(cellc);
      if (to_cpu==myid) {
        cellc = local_cell_coord(input->ort X(0), input->ort Y(0));
        move_atom(cellc, input, 0);
      };
#else
      move_atom(cellc, input, 0);
#endif
  }
} 

#endif /* TWOD */


#ifndef TWOD

/* initialize for cubic crystals */
void init_cubic(void)
{
  if (box_param.x != 0) {
    box_x.x = box_param.x; box_x.y = 0.0;         box_x.z = 0.0;
    box_y.x = 0.0;         box_y.y = box_param.y; box_y.z = 0.0;
    box_z.x = 0.0;         box_z.y = 0.0;         box_z.z = box_param.z;
  }
}

/* generate fcc or NaCl crystal */
void generate_fcc(int maxtyp)
{
  cell    *input;
  vektor  min, max;
  ivektor cellc;
  int     to_cpu;
  int     x, y, z, typ, sign;

#ifdef MPI
  min.x =  my_coord.x      * box_x.x / cpu_dim.x;
  max.x = (my_coord.x + 1) * box_x.x / cpu_dim.x;
  min.y =  my_coord.y      * box_y.y / cpu_dim.y;
  max.y = (my_coord.y + 1) * box_y.y / cpu_dim.y;
  min.z =  my_coord.z      * box_z.z / cpu_dim.z;
  max.z = (my_coord.z + 1) * box_z.z / cpu_dim.z;

#else
  min.x = 0; max.x = box_x.x;
  min.y = 0; max.y = box_y.y;
  min.z = 0; max.z = box_z.z;
#endif

  /* Set up 1 atom input cell */
  input = (cell *) malloc(sizeof(cell));
  if (0==input) error("Can't allocate input cell.");
  input->n_max=0;
  alloc_cell(input, 1);
#ifndef MONOLJ
  input->masse[0] = 1.0;
#endif

  natoms=0;
  nactive=0;

  for (x=min.x ; x<max.x; x++)
    for (y=min.y; y<max.y; y++)
      for (z=min.z; z<max.z; z++) {
 
        sign = 1;
        typ  = (x+y+z) % 2;

#if defined(FRAC) || defined(DEFORM)
        /* leave boundary open if necessary */
        if ((x+0.5 < strip_width/2) || (x+0.5 > box_x.x-strip_width/2) ||
            (y+0.5 < strip_width/2) || (y+0.5 > box_y.y-strip_width/2) ||
            (z+0.5 < strip_width/2) || (z+0.5 > box_z.z-strip_width/2)) continue;
        /* boundary atoms in x-direction get negatiave numbers */
        if ((x+0.5 < strip_width)   || (x+0.5 > box_x.x-strip_width)) sign=-1;
#endif

        /* if fcc, only atoms of type 0 */
        if (typ > maxtyp) continue;

	++natoms;
        if (sign>0) nactive++;

        input->n = 1;
        input->ort X(0) = x + 0.5;
        input->ort Y(0) = y + 0.5;
        input->ort Z(0) = z + 0.5;
        cellc = cell_coord(input->ort X(0), input->ort Y(0), input->ort Z(0));
#ifndef MONOLJ
        input->nummer[0] = sign * natoms;
        input->sorte[0] = typ;
        num_sort[typ]++;
#else
        num_sort[0]++;
#endif

#ifdef MPI
	to_cpu = cpu_coord(cellc);
        if (to_cpu==myid) {
	  cellc = local_cell_coord(input->ort X(0), input->ort Y(0), 
                                   input->ort Z(0));
	  move_atom(cellc, input, 0);
	};
#else
	move_atom(cellc, input, 0);
#endif

  };
} 

/* generate a cubic Laves structure crystal */

void generate_lav()
{
  cell    *input;
  vektor  min, max, lmin, lmax, rmax, rmin;
  ivektor cellc;
  int     to_cpu;
  real     x, y, z, co;
  real     px[24],py[24],pz[24];
  int     i,j,k,l,typ,sign,pa[24];

  co = sqrt(2.)/8.;

  px[0] = 0; py[0] = 0; pz[0] = 0; pa[0] = 0;
  px[1] = 2; py[1] = 2; pz[1] = 0; pa[1] = 0;
  px[2] = 2; py[2] = 0; pz[2] = 2; pa[2] = 0;
  px[3] = 0; py[3] = 2; pz[3] = 2; pa[3] = 0;
  px[4] = 3; py[4] = 3; pz[4] = 3; pa[4] = 1;
  px[5] = 5; py[5] = 5; pz[5] = 5; pa[5] = 1;

  px[6] = 4; py[6] = 4; pz[6] = 0; pa[6] = 0;
  px[7] = 6; py[7] = 6; pz[7] = 0; pa[7] = 0;
  px[8] = 6; py[8] = 4; pz[8] = 2; pa[8] = 0;
  px[9] = 4; py[9] = 6; pz[9] = 2; pa[9] = 0;
  px[10] = 7; py[10] = 7; pz[10] = 3; pa[10] = 1;
  px[11] = 1; py[11] = 1; pz[11] = 5; pa[11] = 1;

  px[12] = 4; py[12] = 0; pz[12] = 4; pa[12] = 0;
  px[13] = 6; py[13] = 2; pz[13] = 4; pa[13] = 0;
  px[14] = 6; py[14] = 0; pz[14] = 6; pa[14] = 0;
  px[15] = 4; py[15] = 2; pz[15] = 6; pa[15] = 0;
  px[16] = 7; py[16] = 3; pz[16] = 7; pa[16] = 1;
  px[17] = 1; py[17] = 5; pz[17] = 1; pa[17] = 1;

  px[18] = 0; py[18] = 4; pz[18] = 4; pa[18] = 0;
  px[19] = 2; py[19] = 6; pz[19] = 4; pa[19] = 0;
  px[20] = 2; py[20] = 4; pz[20] = 6; pa[20] = 0;
  px[21] = 0; py[21] = 6; pz[21] = 6; pa[21] = 0;
  px[22] = 3; py[22] = 7; pz[22] = 7; pa[22] = 1;
  px[23] = 5; py[23] = 1; pz[23] = 1; pa[23] = 1;

#ifdef MPI

  max.x = (int)(box_x.x/sqrt(8.)+0.5);
  max.y = (int)(box_y.y/sqrt(8.)+0.5);
  max.z = (int)(box_z.z/sqrt(8.)+0.5);

  if (myid==0) printf("Zahl der Atome: %d\n",(int)(max.x*max.y*max.z)*24);

  rmin.x = my_coord.x      * box_x.x / cpu_dim.x;
  rmax.x = (my_coord.x + 1) * box_x.x / cpu_dim.x;
  rmin.y = my_coord.y      * box_y.y / cpu_dim.y;
  rmax.y = (my_coord.y + 1) * box_y.y / cpu_dim.y;
  rmin.z = my_coord.z      * box_z.z / cpu_dim.z;
  rmax.z = (my_coord.z + 1) * box_z.z / cpu_dim.z;

  lmin.x = (my_coord.x -1)     * box_x.x / cpu_dim.x;
  lmax.x = (my_coord.x + 2) * box_x.x / cpu_dim.x;
  lmin.y = (my_coord.y -1)     * box_y.y / cpu_dim.y;
  lmax.y = (my_coord.y + 2) * box_y.y / cpu_dim.y;
  lmin.z = (my_coord.z -1)     * box_z.z / cpu_dim.z;
  lmax.z = (my_coord.z + 2) * box_z.z / cpu_dim.z;

  min.x = (int)(lmin.x/sqrt(8.)+0.5);
  max.x = (int)(lmax.x/sqrt(8.)+0.5);
  min.y = (int)(lmin.y/sqrt(8.)+0.5);
  max.y = (int)(lmax.y/sqrt(8.)+0.5);
  min.z = (int)(lmin.z/sqrt(8.)+0.5);
  max.z = (int)(lmax.z/sqrt(8.)+0.5);

  /*printf("%d %f %f %f %f %f %f\n",myid,min.x,max.x,min.y,max.y,min.z,max.z);*/

#else

  min.x=0;
  min.y=0;
  min.z=0;
  max.x=(int)(box_x.x/sqrt(8.)+0.5);
  max.y=(int)(box_y.y/sqrt(8.)+0.5);
  max.z=(int)(box_z.z/sqrt(8.)+0.5);

  printf("Zahl der Atome: %d\n",(int)(max.x*max.y*max.z)*24);

#endif

  /* Set up 1 atom input cell */
  input = (cell *) malloc(sizeof(cell));
  if (0==input) error("Can't allocate input cell.");
  input->n_max=0;
  alloc_cell(input, 1);
#ifndef MONOLJ
  input->masse[0] = 1.0;
#endif

  natoms=0;
  nactive=0;

  for (i=min.x ; i<max.x; i++)
    for (j=min.y; j<max.y; j++)
      for (k=min.z; k<max.z; k++) 

	for (l=0;l<24;l++) {
	  x=(px[l]+8.*i)/sqrt(8.);
	  y=(py[l]+8.*j)/sqrt(8.);
	  z=(pz[l]+8.*k)/sqrt(8.);
	  typ=pa[l];
          sign=1;

/*	  if ((x+co < rmin.x) || (x+co > rmax.x) ||
	      (y+co < rmin.y) || (y+co > rmax.y) ||
	      (z+co < rmin.z) || (z+co > rmax.z)) continue; */

#if defined(FRAC) || defined(DEFORM)
	  /* leave boundary open if necessary */
	  if ((x+co < strip_width/2) || (x+co > box_x.x-strip_width/2) ||
	      (y+co < strip_width/2) || (y+co > box_y.y-strip_width/2) ||
	      (z+co < strip_width/2) || (z+co > box_z.z-strip_width/2)) continue;
	  /* boundary atoms in x-direction get negative numbers */
	  if ((x+co < strip_width) || (x+co > box_x.x-strip_width)) sign=-1;
#endif
	  
	  ++natoms;
          if (sign>0) nactive++;

	  input->n = 1;
	  input->ort X(0) = x + co;
	  input->ort Y(0) = y + co;
	  input->ort Z(0) = z + co;
	  cellc = cell_coord(input->ort X(0), input->ort Y(0), 
                                              input->ort Z(0));
#ifndef MONOLJ
	  input->nummer[0] = sign * natoms;
	  input->sorte[0] = typ;
	  num_sort[typ]++;
#endif
	  
#ifdef MPI
	  to_cpu = cpu_coord(cellc);
	  if (to_cpu==myid) {
	    cellc = local_cell_coord(input->ort X(0), input->ort Y(0), 
				     input->ort Z(0));
	    move_atom(cellc, input, 0);
	  }
#else
	  move_atom(cellc, input, 0);
#endif
	  
	}
} 


#endif /* not TWOD */


