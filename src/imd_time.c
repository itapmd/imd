
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
* imd_time.c -- IMD timer routines
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/******************************************************************************
*
*  set timer.start to the current time
*
******************************************************************************/

void imd_start_timer(imd_timer *timer)
{
#ifdef MPI
  timer->start = MPI_Wtime();
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
  timer->total += MPI_Wtime() - timer->start;
#elif defined(USE_RUSAGE)
  timeval now;
  getrusage(RUSAGE_SELF,&now);
  timer->total += 
    (double)(now.ru_utime.tv_sec  - timer->start.ru_utime.tv_sec) +
    (double)(now.ru_utime.tv_usec - timer->start.ru_utime.tv_usec)/1E6;
#else
  struct tms now;
  double inv_tick = 0.0;
  if (inv_tick==0.0) inv_tick = 1.0 / (double) sysconf(_SC_CLK_TCK);
  times(&now);
  timer->total += (double)(now.tms_utime - timer->start.tms_utime) * inv_tick;
#endif
}




