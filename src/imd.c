/******************************************************************************
*
* imd -- The ITAP Molecular Dynamics Package
*
* Author: J. Stadler, ITAP, University of Stuttgart
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/


#define MAIN

/* Include file also declares global Variables and has Function Prototypes */
#include "imd.h"
#ifdef USE_RUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#endif

/* Main module of the IMD Package
   Various versions are built by conditional compilation */

void main(int argc, char **argv)

{
  real tmp;
    
  int start;
  time_t tstart, tend;
  real total_rtime = (real)0;
#ifndef USE_CLOCK
  struct tms utstart,utend;
#else
  clock_t utend;
#endif
#ifdef USE_RUSAGE
  /* getrusage() liefert auf manchen Systemen genauere Zeiten als times() */
  struct rusage start_ru, end_ru;
  struct timeval start_tval,end_tval;
#endif
  real total_time;
    
  time(&tstart);
#ifndef USE_CLOCK
  times(&utstart);
#endif
#ifdef USE_RUSAGE
  getrusage(RUSAGE_SELF,&start_ru);
#endif

  /* Fire up MPI */
#ifdef MPI
  init_mpi(argc,argv);
  time_start = MPI_Wtime();
  time_last = time_start;
  time_comm = 0.0;
  time_io   = 0.0;
  time_setup= 0.0;
  time_comm_force=0.0;
  time_comm_ar=0.0;
  time_calc_local=0.0;
  time_calc_nonlocal=0.0;
#endif

  /* read command line and parameters for first simulation phase */
  read_parameters(argc,argv);

#ifdef MPI
  broadcast_params();
#endif

  /* Initialize socket i/o */
#ifdef USE_SOCKETS
#ifdef MPI
  if (myid == 0)
#endif /* MPI */
    init_client();
#endif /* USE_SOCKETS */

#ifndef MONOLJ
  /* Read Potential from file */
  read_potential(potfilename);
#ifdef TTBP
  /* Read TTBP Potential from file */
  read_ttbp_potential(ttbp_potfilename);
#endif
#endif

  /* Filenames starting with denote internal 
  generation of the intitial configuration */
  if ('.' == infilename[0]) {
    if (0 == myid) { 
      printf("Generating atoms: %s.\n", infilename);fflush(stdout);
    }
    generate_atoms(infilename);
  }
  else {
    if (0 == myid) {
      printf("Reading atoms.\n");fflush(stdout);
    }
    init_cells();
    read_atoms(infilename);
  }

#ifdef MPI
  if (0 == myid) 
#endif
  printf("Done reading atoms.\n");

  start = steps_min;  /* keep starting step number */

  /* first phase of the simulation */
  if (steps_min <= steps_max) main_loop();

  /* execute further phases of the simulation */
  while (finished==0) {
    simulation++;
    steps_min = steps_max;
#ifdef MPI
    if (0==myid)
#endif
    getparamfile( paramfilename, simulation );
#ifdef MPI
    broadcast_params();
#endif
    if (steps_min <= steps_max) {
      init_cells();  /* a new cell division might be necessary or useful */
      main_loop();
    };
  };
  
/* Write execution time summary */

#ifdef MPI
  if (0 == myid) {
    time_stop = MPI_Wtime();
#endif
#ifndef USE_CLOCK
    times(&utend);
#endif
    time(&tend);
#ifdef USE_RUSAGE
    memset(&end_ru,0,sizeof(struct rusage));
    getrusage(RUSAGE_SELF,&end_ru);
#endif

    steps_max -= start;

#ifdef USE_RUSAGE
    total_rtime =  (real)(utend.tms_utime - utstart.tms_utime)/(real)CLK_TCK;
    total_time = (double)end_ru.ru_utime.tv_sec+(double)end_ru.ru_utime.tv_usec/1E6-
      ((double)start_ru.ru_utime.tv_sec + (double)start_ru.ru_utime.tv_usec/1E6);
#else
#ifndef USE_CLOCK
    total_time =  (real)(utend.tms_utime - utstart.tms_utime)/(real)CLK_TCK;
#else
    total_time = (real)clock()/(real)CLOCKS_PER_SEC;
#endif
#endif
  
    printf("Done simulation.\n\n");
    printf("%s\n\n",progname);
    printf("started at  %s", ctime(&tstart));
    printf("finished at %s\n", ctime(&tend));
#ifdef USE_RUSAGE
/* Zum Testen eingefuegt -- MH */
#ifndef __linux__ /* not yet implemented in the Linux 2.0.x kernel */
    printf("Maximum resident size RSS %lu Bytes\n",
           (unsigned long)end_ru.ru_maxrss);
    printf("Code  IXRSS=%lu Bytes\n",(unsigned long)end_ru.ru_ixrss);
    printf("Data  IDRSS=%lu Bytes\n",(unsigned long)end_ru.ru_idrss);
    printf("Stack ISRSS=%lu Bytes\n",(unsigned long)end_ru.ru_isrss);
#endif
#endif
    printf("Used %f seconds cputime,\ndid %d steps with %d atoms\nfor %.3e cpuseconds per step and atom.\n",
	   total_time, 
           steps_max, 
           natoms, 
           total_time/ (steps_max * natoms));
    printf("Inverse is %.3e steps * atom per cpusecond.\n",
	   ((float)steps_max*(float)natoms)/total_time);
    printf("Time for pair interaction t_pair is: %.3e seconds\n\n",
           total_time / ((float)steps_max*(float)natoms*
                         (4.0/3.0)*M_PI*pow(r2_cut,3.0/2.0)*(natoms/volume)));
#ifdef MPI

    time_comm += time_comm_force + time_comm_ar;
    time_calc += time_calc_local + time_calc_nonlocal;

    total_time = time_stop - time_start;
    printf("MPI timings:\n");
    printf("Total time: %e seconds.\n",total_time);
    printf("Comp. time: %e seconds or %.1f %% of total.\n",
           time_calc,100*time_calc/total_time);
    printf("I/O   time: %e seconds or %.1f %% of total.\n",
           time_io,100*time_io/total_time);
    printf("Comm. time: %e seconds or %.1f %% of total.\n",
           time_comm,100*time_comm/total_time);
    printf("Setup time: %e seconds or %.1f %% of total.\n\n",
           time_setup,100*time_setup/total_time);
    printf("Comm. time force   : %e seconds.\n",time_comm_force);
    printf("Comm. time act=rea : %e seconds.\n",time_comm_ar   );
    printf("Comp. time local   : %e seconds.\n",time_calc_local);
    printf("Comp. time nonlocal: %e seconds.\n",time_calc_nonlocal);

    printf("\nExcluding Setup time, using MPI Timings:\n\n");
    total_time -= time_setup;
    printf("Used %f seconds cputime,\ndid %d steps with %d atoms\nfor %.3e cpuseconds per step and atom.\n",
	   total_time, 
           steps_max, 
           natoms, 
           total_time/ (steps_max*natoms));
    printf("Inverse is %.3e steps * atom per cpusecond.\n",
	   ((float)steps_max*(float)natoms)/total_time);
    printf("Time for pair interaction t_pair is: %.3e seconds\n",
           total_time / (steps_max*natoms*
                         (4.0/3.0)*M_PI*pow(r2_cut,3.0/2.0)*(natoms/volume)));
    printf("Time for pair interaction per cpu t_pair(P=1) is: %.3e seconds\n",
    total_time*num_cpus/ ((float)steps_max*(float)natoms*
                          (4.0/3.0)*M_PI*pow(r2_cut,3.0/2.0)*(natoms/volume)));
    printf("Comp. time: %e seconds or %.1f %% of total.\n",
           time_calc,100*time_calc/total_time);
    printf("I/O   time: %e seconds or %.1f %% of total.\n",
           time_io,100*time_io/total_time);
    printf("Comm. time: %e seconds or %.1f %% of total.\n",
           time_comm,100*time_comm/total_time);
  };

  /* Kill MPI */
  shutdown_mpi();
#endif

  exit(0);

}
