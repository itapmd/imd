//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifndef _VVSTOPWATCH_H_
#define _VVSTOPWATCH_H_

#ifdef WIN32
  #include <time.h>
#else
  #include <sys/time.h>
#endif

/** System independent implementation of a stop watch.
  The stop watch can be started, stopped, resetted, and its time
  can be read. Stopping does not reset the watch. Initially, the
  stop watch is resetted. Resetting stops the watch. Reading the time
  does not stop the watch. <P>

  Example usage:<PRE>

  vvStopwatch sw = new vvStopwatch(); // create new stop watch instance
  sw->start();                        // reset counter
  // *** do something ***
  float time1 = sw->getTime();        // get an intermediate time but don't stop counting
  // *** do something ***
  sw->stop();                         // stop counting
  float time2 = sw->getTime();        // get another intermediate time
  sw->start();                        // continue counting
  // *** do something ***
  float time3 = sw->getTime();        // get the final time
  sw->reset();                        // stop and reset the watch for next usage
  cout << "The total time measured is: " << time3 << " seconds." << endl;
  delete sw;                          // remove stop watch from memory
  </PRE>

  @author Juergen Schulze (schulze@hlrs.de)
*/
class vvStopwatch
{
  private:
#ifdef WIN32
    clock_t baseTime;         ///< system time when stop watch was triggered last
    LARGE_INTEGER baseTimeQP; ///< base time when using QueryPerformance API
    LARGE_INTEGER freq;       ///< frequency if QueryPerformance API is used
    bool useQueryPerformance; ///< true=use QueryPerformance API
#else
    timeval baseTime;         ///< system time when stop watch was triggered last
#endif
    float accTime;            ///< additional passed time to current measurement (accumulated time)
    bool  running;            ///< true = stop watch is running

    float getCurrentMeasurement();

  public:
    vvStopwatch();
    void  start();
    void  stop();
    void  reset();
    float getTime();
};

#endif

//============================================================================
// End of File
//============================================================================
