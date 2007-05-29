
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2007 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_time.c -- IMD timer routines
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"
#ifdef BGL
#include <rts.h>
double bgl_clockspeed=1.0e-6/700.0;
#endif

/******************************************************************************
*
*  initialize IMD timer
*
******************************************************************************/

void imd_init_timer(imd_timer *timer, int flag, char *description, char *color)
{
#if defined(MPI) && defined(MPE)
  timer->mpe_flag = flag;
  if (flag) {
    MPE_Log_get_state_eventIDs( &(timer->mpe_id_begin), &(timer->mpe_id_end));
    MPE_Describe_state( timer->mpe_id_begin, timer->mpe_id_end, 
                        description, color );
  }
#endif
  timer->total = 0.0;
}

/******************************************************************************
*
*  set timer.start to the current time
*
******************************************************************************/

void imd_start_timer(imd_timer *timer)
{
#ifdef MPI
#ifdef MPE
  if (timer->mpe_flag) MPE_Log_event( timer->mpe_id_begin, 0, NULL );
#endif
#ifdef BGL
  timer->start = rts_get_timebase() * bgl_clockspeed;
#else
  timer->start = MPI_Wtime();
#endif
#elif defined(USE_WALLTIME)
  gettimeofday(&(timer->start),NULL);
#elif defined(OMP)
  timer->start = omp_get_wtime();
#elif defined(USE_RUSAGE)
  getrusage(RUSAGE_SELF,&(timer->start));
#else
  times(&(timer->start)); 
#endif
}

/******************************************************************************
*
*  add (current time - timer.start) to timer.total
*
******************************************************************************/

void imd_stop_timer(imd_timer *timer)
{
#ifdef MPI
#ifdef MPE
  if (timer->mpe_flag) MPE_Log_event( timer->mpe_id_end, 0, NULL );
#endif
#ifdef BGL
  timer->total += rts_get_timebase() * bgl_clockspeed - timer->start;
#else
  timer->total += MPI_Wtime() - timer->start;
#endif
#elif defined(USE_WALLTIME)
  struct timeval now;
  gettimeofday(&now,NULL);
  timer->total += 
    (double)(now.tv_sec  - timer->start.tv_sec) +
    (double)(now.tv_usec - timer->start.tv_usec)/1E6;
#elif defined(OMP)
  timer->total += omp_get_wtime() - timer->start;
#elif defined(USE_RUSAGE)
  struct rusage now;
  getrusage(RUSAGE_SELF,&now);
  timer->total += 
    (double)(now.ru_utime.tv_sec  - timer->start.ru_utime.tv_sec) +
    (double)(now.ru_utime.tv_usec - timer->start.ru_utime.tv_usec)/1E6;
#else
  struct tms now;
  double inv_tick = 0.0;
  if (inv_tick==0.0) inv_tick = 1.0 / (double) sysconf(_SC_CLK_TCK);
  times(&now);
  timer->total += (double)(now.tms_utime - timer->start.tms_utime)*inv_tick;
#endif
}
