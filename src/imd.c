
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
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
#ifdef PAPI
#include <papi.h>
#endif

/* Main module of the IMD Package;
   Various versions are built by conditional compilation */

int main(int argc, char **argv)
{
  int start, num_threads;
  time_t tstart, tend;
  real tmp;
#ifdef PAPI
  float rtime, ptime, mflops;
  long_long flpins;
#endif
  
#ifdef MPI
  init_mpi(&argc,argv);
#endif

  is_big_endian = endian();
  read_command_line(argc,argv);

  /* loop for online visualization */
  do {

  time(&tstart);

  /* reset timers */
  time_total.total = 0.0;
  time_setup.total = 0.0;
  time_main.total  = 0.0;
  time_io.total    = 0.0;

  /* start some timers (after starting MPI!) */
  imd_start_timer(&time_total);
  imd_start_timer(&time_setup);

  /* read parameters for first simulation phase */
  read_parameters();

  /* initialize random number generator */
  srand48(seed);

#ifdef PAIR
  /* read pair potential file - also used for TTBP, EAM2, TERSOFF, EWALD, .. */
  if (have_potfile)
    read_pot_table(&pair_pot,potfilename,ntypes*ntypes,1);
  /* initialize analytically defined potentials */
  if (have_pre_pot) init_pre_pot();
#ifdef LINPOT
  make_lin_pot_table(pair_pot, &pair_pot_lin);
#endif
#endif

#ifdef TTBP
  /* read TTBP smoothing potential file */
  read_pot_table(&smooth_pot,ttbp_potfilename,ntypes*ntypes,1);
#endif

#ifdef EAM2
  /* read the tabulated embedding energy function */
  read_pot_table(&embed_pot,eam2_emb_E_filename,ntypes,0);
  /* read the tabulated electron density function */
  read_pot_table(&rho_h_tab,eam2_at_rho_filename,ntypes*ntypes,1);
#endif

#ifdef MEAM
  if (have_potfile)
    read_pot_table(&pair_pot,potfilename,ntypes*ntypes,1);
  /* read the tabulated embedding energy function */
  if (have_embed_potfile)
    read_pot_table(&embed_pot,meam_emb_E_filename,ntypes,0);
  /* read the tabulated electron density function */
  if (have_eldensity_file)
    read_pot_table(&el_density,meam_eldensity_filename,ntypes,1);    
  init_meam();
#endif

#ifdef KEATING
  init_keating();
#endif
#ifdef STIWEB
  init_stiweb();
#endif
#ifdef TERSOFF
  init_tersoff();
#endif

  /* filenames starting with _ denote internal 
     generation of the intitial configuration */
  if ('_' == infilename[0]) {
    if (0 == myid) { 
      printf("Generating atoms: %s.\n", infilename);fflush(stdout);
    }
    generate_atoms(infilename);
  }
  else {
    if (0 == myid) {
      printf("Reading atoms.\n");fflush(stdout);
    }
    make_box();
    read_atoms(infilename);
  }
  if (0 == myid) printf("Done reading atoms.\n");

#ifdef EPITAX
  if (0 == myid) 
    printf("EPITAX: Largest substrate atom number: %d\n", epitax_sub_n);
  epitax_number = epitax_sub_n;
  epitax_level  = substrate_level();
  check_boxheight();
#endif

#ifdef EWALD
  init_ewald();
#endif

  /* initialize socket I/O */
#ifdef USE_SOCKETS
  if ((myid == 0) && (socket_int > 0)) init_socket();
#endif

  start = steps_min;  /* keep starting step number */

  /* write .eng file header */
  if ((imdrestart==0) && (eng_int>0)) write_eng_file_header();

  imd_stop_timer(&time_setup);
  imd_start_timer(&time_main);
#ifdef PAPI
  PAPI_flops(&rtime,&ptime,&flpins,&mflops);
#endif

  /* first phase of the simulation */
  if (steps_min <= steps_max) main_loop();

  /* execute further phases of the simulation */
  while (finished==0) {
    int tmp;
    simulation++;
    if (0==myid) {
      getparamfile( paramfilename, simulation );
      sprintf( itrfilename,"%s-final.itr", outfilename );
      tmp = finished;
      getparamfile( itrfilename, 1 );
      finished = tmp;
    }
#ifdef MPI
    broadcast_params();
#endif

#ifdef KEATING
    init_keating();
#endif
#ifdef STIWEB
    init_stiweb();
#endif
#ifdef TERSOFF
    init_tersoff();
#endif

    if (steps_min <= steps_max) {
      make_box();  /* make sure the box size is still ok */
      main_loop();
    }
  }

  imd_stop_timer(&time_main);
  imd_stop_timer(&time_total);
#ifdef PAPI
  PAPI_flops(&rtime,&ptime,&flpins,&mflops);
#endif

  /* write execution time summary */
  if (0 == myid) {
    if (NULL!= eng_file) fclose( eng_file);
    if (NULL!=msqd_file) fclose(msqd_file);
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

#ifdef EPITAX
    if (0 == myid) printf("EPITAX: %d atoms created.\n", nepitax);
#endif

#ifdef OMP
    num_threads = omp_get_max_threads();
#else
    num_threads = 1;
#endif

#ifdef MPI
    printf("Did %d steps with %d atoms on %d CPUs.\n", 
           steps_max+1, natoms, num_cpus*num_threads);
#ifdef OMP
    printf("Used %d processes with %d threads each.\n", num_cpus, num_threads);
#endif
#else
#ifdef OMP
    printf("Did %d steps with %d atoms on %d CPUs.\n", 
           steps_max+1, natoms, num_threads);
#else
    printf("Did %d steps with %d atoms.\n", steps_max+1, natoms);
#endif
#endif

    printf("Used %f seconds cputime,\n", time_total.total);
    printf("%f seconds excluding setup time,\n", time_main.total);
    tmp =  (num_cpus * num_threads * time_main.total / (steps_max+1)) / natoms;
    printf("%e cpuseconds per step and atom\n", tmp);
    printf("(inverse is %e).\n\n", 1.0/tmp);

#ifdef TIMING
    printf("Input/Output time:   %e seconds or %.1f %% of main loop\n",
           time_io.total,100*time_io.total/time_main.total);
#endif
  
     fflush(stdout);
     fflush(stderr);
  }

  /* if looping, clean up old simulation */
  if (loop) {

    int k;

    /* empty all cells */
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
#endif

    volume_init = 0.0;

    /* force init_cells */
    max_height.x = 0.0;
    max_height.y = 0.0;
#ifndef TWOD
    max_height.z = 0.0;
#endif

    /* close open files */
    if (myid==0) {
      if ( eng_file!=NULL) fclose( eng_file);  eng_file=NULL;
      if (msqd_file!=NULL) fclose(msqd_file); msqd_file=NULL;
    }
  }

  } while (loop);

  /* kill MPI */
#ifdef MPI
  shutdown_mpi();
#endif

  exit(0);

}
