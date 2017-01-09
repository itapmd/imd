
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2011 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#define MAIN

/* include file also declares global variables and has function prototypes */
#include "imd.h"

/* Performance Application Programming Interface */
#ifdef PAPI
#include <papi.h>
#endif

/* Main module of the IMD Package;
   Various versions are built by conditional compilation */

int main(int argc, char **argv)
{
  int start;
  int simulation = 1, finished = 0;

#ifdef PAPI
  float rtime, ptime, mflops;
  long_long flpins;
#endif

#if defined(MPI) || defined(NEB)
  MPI_Init(&argc, &argv);
  init_mpi();
#endif

  /* initialize timers */
  imd_init_timer( &time_total,      0, NULL,        NULL    );
  imd_init_timer( &time_setup,      1, "setup",     "white" );
  imd_init_timer( &time_main,       0, NULL,        NULL    );
  imd_init_timer( &time_output,     1, "output",    "cyan"  );
  imd_init_timer( &time_input,      1, "input",     "orange");
  imd_init_timer( &time_integrate,  1, "integrate", "green" );
  imd_init_timer( &time_forces,     1, "forces",    "yellow");

  read_command_line(argc,argv);

  time(&tstart);

  /* start some timers (after starting MPI!) */
  imd_start_timer(&time_total);
  imd_start_timer(&time_setup);

  /* read parameters for first simulation phase */
  finished = read_parameters(paramfilename, simulation);

  /* initialize all potentials */
  setup_potentials();

#ifdef TIMING
  imd_start_timer(&time_input);
#endif
  
#ifdef NEB
  read_atoms_neb(infilename);
#else
  /* filenames starting with _ denote internal generation of configuration */
  if ('_' == infilename[0]) {
    generate_atoms(infilename);
  }
  else {
    read_atoms(infilename);
  }
#endif

#ifdef TIMING
  imd_stop_timer(&time_input);
#endif

  start = steps_min;  /* keep starting step number */

  /* write .eng file header */
  imd_write_eng_headers();
  
  /* initialize modules */
  imd_init_modules();

#ifdef PAPI
  PAPI_flops(&rtime,&ptime,&flpins,&mflops);
#endif

  imd_stop_timer(&time_setup);

  /* first phase of the simulation */
  if (steps_min <= steps_max) {
    if ((0==myid) && (0==myrank))
      printf("Starting simulation %d\n", simulation);
    if (main_loop(simulation)) finished = 1;
  }

  /* execute further phases of the simulation */
  while (finished==0) {
    simulation++;
    finished = read_parameters(paramfilename, simulation);
    if (steps_min <= steps_max) {
      make_box();  /* make sure the box size is still ok */
      if ((0==myid) && (0==myrank))
        printf("Starting simulation %d\n", simulation);
      if (main_loop(simulation)) finished = 1;
    }
  }

  imd_stop_timer(&time_total);

#ifdef PAPI
  PAPI_flops(&rtime,&ptime,&flpins,&mflops);
#endif

  /* write execution time summary */
  if ((0 == myid) && (0 == myrank)){

    time(&tend);
    steps_max -= start;
    if (steps_max < 0)
      error("Start step greater than maximum number of steps");

    printf("Done simulation.\n\n");
    printf("%s\n\n",progname);
    printf("started at %s", ctime(&tstart));
    printf("finished at %s\n", ctime(&tend));

#ifdef PAPI
    printf("Achieved %f megaflops per CPU\n\n", mflops);
#endif
    
    imd_print_summary();
  }

  imd_cleanup();

  /* kill MPI */
#if defined(MPI) || defined(NEB)
  shutdown_mpi();
#endif

  return 0;
}

/*****************************************************************************/

void imd_write_eng_headers( void )
{
  if ((imdrestart==0) && (eng_int>0)) write_eng_file_header();
#ifdef EXTPOT
  if ((imdrestart==0) && (eng_int>0)) write_fext_header();
#endif
#ifdef RELAX
  if ((imdrestart==0) && (eng_int>0))
  {
      if ( (ensemble==ENS_MIK) || (ensemble==ENS_GLOK) || (ensemble==ENS_CG) )
          write_ssdef_header();
  }
#endif
}

void imd_init_modules( void )
{

#ifdef EPITAX
  if (0 == myid)
    printf("EPITAX: Largest substrate atom number: %d\n", epitax_sub_n);
  epitax_number = epitax_sub_n;
  epitax_level  = substrate_level();
  check_boxheight();
#endif

#ifdef REFPOS
  init_refpos();
#endif

#ifdef BBOOST
  init_bboost();
#endif
  
#ifdef BEND
  init_bend();
#endif

#ifdef ZAPP
  init_zapp();
#endif

#ifdef KIM
  init_kim();
#endif
  
#ifdef CNA
  init_cna();
#endif

#ifdef ADA
  init_ada();
#endif

#ifdef NYETENSOR
  init_NyeTensor();
#endif

#ifdef LASER
  init_laser();
#endif

#ifdef CYCLE
  init_cycle();
#endif

#ifdef TTM
  init_ttm();
#endif

#ifdef SM
  init_sm();
#endif

#ifdef LOADBALANCE
  lb_computeVariance();
  if (myid == 0 && (0==myrank)){
    printf("LOAD BALANCING: Initial max_load %f\n", lb_maxLoad);
    printf("LOAD BALANCING: Initial min_load %f\n", lb_minLoad);
    printf("LOAD BALANCING: Initial variance %f\n",  lb_loadVariance);
  }

  if (lb_preRuns > 0){
      fix_cells();
      setup_buffers();
      for (i=0; i<lb_preRuns;i++){
          int success = balanceLoad(lb_getLoad, lb_getCenterOfGravity, 0, 0);
          lb_computeVariance();
          write_lb_file(-lb_preRuns+i, 1);
          if (!success) break;
      }
      if (lb_writeStatus)
            write_lb_status(0);
      if (myid == 0 && (0==myrank)){
        printf("LOAD BALANCING: After initial runs max_load %f\n", lb_maxLoad);
        printf("LOAD BALANCING: After initial runs min_load %f\n", lb_minLoad);
        printf("LOAD BALANCING: After initial runs variance %f\n", lb_loadVariance);
      }
  }
#endif
  
}

void imd_cleanup( void )
{

  /* close open files */
  if (NULL!= eng_file) fclose(eng_file);

#ifdef EXTPOT
  if (NULL!= ind_file) fclose(ind_file);
#endif

#ifdef LOADBALANCE
  if (NULL!= lblog_file) fclose(lblog_file);
#endif
  
#ifdef RELAX
  if (NULL!= ssdef_file) fclose(ssdef_file);
#endif

  if (NULL!=msqd_file) fclose(msqd_file);
  
  /* empty all cells */
  int k;
  for (k=0; k<nallcells; k++) cell_array[k].n = 0;

  /* free potential tables */
#if defined(PAIR)
  free_pot_table(&pair_pot);
#endif
#ifdef TTBP
  free_pot_table(&smooth_pot);
#endif
#ifdef EAM2
  free_pot_table(&embed_pot);
  free_pot_table(&rho_h_tab);
#ifdef EEAM
  free_pot_table(&emod_pot);
#endif
#endif
#ifdef SM
  free_pot_table(&na_pot_tab);
  free_pot_table(&cr_pot_tab);
  free_pot_table(&erfc_r_tab);
#endif
#ifdef KIM
  destroy_kim();
#endif

}

void imd_print_summary( void )
{

#ifdef NBLIST
  printf("Neighbor list update every %d steps on average\n\n",
         steps_max / MAX(nbl_count,1));
#endif

#ifdef EPITAX
  printf("EPITAX: %d atoms created.\n", nepitax);
#endif
    
#ifdef NEB
  real Emax=-999999;
  real Emin=999999;
  int maxi=0;

  printf ("NEB:\n # Image Epot\n");
  for(i=0;i<neb_nrep;i++){
      if ( neb_epot_im[i] < Emin)
          Emin=neb_epot_im[i];
      if ( neb_epot_im[i] > Emax){
          maxi=i;
          Emax=neb_epot_im[i];
      }
      printf(" %d %lf\n",i, neb_epot_im[i]);
  }
  printf ("Saddlepoint: %d Activation Energy: %lf \n",maxi,Emax-Emin);
#endif
  
  int num_threads = 1;
#ifdef OMP
  int num_threads = omp_get_max_threads();
#endif

#ifdef MPI
  printf("Did %d steps with %ld atoms on %d physical CPUs.\n",
         steps_max+1, natoms, num_cpus * num_threads / hyper_threads);
#ifdef OMP
  printf("Used %d processes with %d threads each.\n", num_cpus, num_threads);
#endif
#else
#ifdef OMP
  printf("Did %d steps with %ld atoms on %d physical CPUs.\n",
         steps_max+1, natoms, num_threads / hyper_threads);
#else
  printf("Did %d steps with %ld atoms.\n", steps_max+1, natoms);
#endif
#endif

  printf("Used %f seconds cputime,\n", time_total.total);
  printf("%f seconds excluding setup time,\n", time_main.total);
  real tmp =  ((num_cpus * num_threads / hyper_threads) * time_main.total /
               (steps_max+1)) / natoms;
  printf("%e cpuseconds per step and atom\n", tmp);
  printf("(inverse is %e).\n\n", 1.0/tmp);

#ifdef TIMING
  printf("Output time:   %e seconds or %.1f %% of main loop\n",
         time_output.total, 100*time_output.total /time_main.total);
  printf("Input  time:   %e seconds or %.1f %% of main loop\n",
         time_input.total,100*time_input.total/time_main.total);
  printf("Force  time:   %e seconds or %.1f %% of main loop\n",
         time_forces.total,100*time_forces.total/time_main.total);
#endif
  
  fflush(stdout);
  fflush(stderr);

}
