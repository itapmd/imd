
/******************************************************************************
*
* imd_misc.c -- Some Misc. Routines for the imd package
*
******************************************************************************/

/******************************************************************************
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
  if (myid==0) {
    fprintf(stderr,"%s [-r<nnn>] [-p paramter-file]\n",progname); 
    fflush(stderr);
  }
#ifdef MPI
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  exit(1); 
}

/******************************************************************************
*
*  print a warning
*
******************************************************************************/

void warning(char *msg)
{ 
  if (myid==0) {
    fprintf(stderr,"WARNING: %s\n",msg);
    fflush(stderr);
  }
}

/******************************************************************************
*
* error -- Complain and abort
*
******************************************************************************/

void error(char *msg)
{
#ifdef MPI
  fprintf(stderr,"Error on CPU %d: %s\n",myid,msg);
  fflush(stderr);
  MPI_Abort(MPI_COMM_WORLD, 1);
#else
  fprintf(stderr,"Error: %s\n",msg);
  fflush(stderr);
#endif
  exit(2);
}

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
        char msg[255];
        sprintf(msg,"core_potential_file: numstep: %lf rounded to: %lf",
                numstep, number_of_steps);
        warning(msg);
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
             REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_phi_r_begin, eam2_phi_r_infosize / sizeof(real), 
             REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_phi_r_end, eam2_phi_r_infosize / sizeof(real), 
             REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_phi_r_step, eam2_phi_r_infosize / sizeof(real), 
             REAL, 0, MPI_COMM_WORLD);
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
      char msg[255];
      sprintf(msg,"embedding_energy_file: numstep: %lf rounded to: %lf",
              numstep, number_of_steps);
      warning(msg);
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
             REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_rho_begin, eam2_rho_infosize / sizeof(real), 
             REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_rho_end, eam2_rho_infosize / sizeof(real), 
             REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_rho_step, eam2_rho_infosize / sizeof(real), 
             REAL, 0, MPI_COMM_WORLD);
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
        char msg[255];
        sprintf(msg,"atomic_e-density_file: numstep: %lf rounded to: %lf",
		numstep, number_of_steps);
        warning(msg);
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
             REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_r_begin, eam2_r_infosize / sizeof(real), 
             REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_r_end, eam2_r_infosize / sizeof(real), 
             REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_r_step, eam2_r_infosize / sizeof(real), 
             REAL, 0, MPI_COMM_WORLD);
#endif
} /* end of eam2_read_atomic_rho */

#endif /* EAM2 function tables */

