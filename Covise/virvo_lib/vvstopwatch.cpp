//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifdef WIN32
  #include <windows.h>    // required for QueryPerformance API
#endif
#include "vvstopwatch.h"

//#define VV_STANDALONE      // uncomment for demonstration

//----------------------------------------------------------------------------
/// Constructor. Initializes time variable with zero.
vvStopwatch::vvStopwatch()
{
#ifdef WIN32
  baseTime = 0;
  baseTimeQP.QuadPart = 0;
  useQueryPerformance = (QueryPerformanceFrequency(&freq)) ? true : false;
#else
  baseTime.tv_sec  = 0;
  baseTime.tv_usec = 0;
#endif
  accTime = 0.0f;
  running = false;
}

//----------------------------------------------------------------------------
/// Start or restart measurement but don't reset counter.
void vvStopwatch::start()
{
#ifdef WIN32
  if (useQueryPerformance) QueryPerformanceCounter(&baseTimeQP);
  else baseTime = clock(); 
#elif __linux
  struct timezone tz;
  gettimeofday(&baseTime, &tz);
#else
  void* v = NULL;
  gettimeofday(&baseTime, v);
#endif
  running = true;
}

//----------------------------------------------------------------------------
/** Get the current stop watch time. This does not stop or reset the counter.
  @return time difference in seconds
*/
float vvStopwatch::getTime()
{
  if (running)
    return (accTime + getCurrentMeasurement());
  else
    return accTime;
}

//----------------------------------------------------------------------------
/// Stop measurement but don't reset counter.
void vvStopwatch::stop()
{ 
  accTime += getCurrentMeasurement();
  running = false;
}

//----------------------------------------------------------------------------
/// Reset counter and stop measurement.
void vvStopwatch::reset()
{
  accTime = 0.0f;
  running = false;
}

//----------------------------------------------------------------------------
/// Return the time passed since the last start command [seconds].
float vvStopwatch::getCurrentMeasurement()
{
  float dt;           // measured time difference [seconds]

  if (!running) return 0.0f;

#ifdef WIN32
  if (useQueryPerformance)
  {
    LARGE_INTEGER now;
    QueryPerformanceCounter(&now);
    dt = float((float(now.QuadPart) - float(baseTimeQP.QuadPart)) / float(freq.QuadPart));
  }
  else 
  {
    clock_t now = clock();
    dt = (float)(now - baseTime) / (float)CLOCKS_PER_SEC;
  }
#else

  #ifdef __linux
    struct timezone dummy;
  #else
    void* dummy = NULL;
  #endif
  timeval now;          // current system time
  
  gettimeofday(&now, &dummy);
  time_t sec  = now.tv_sec  - baseTime.tv_sec;
  long   usec = now.tv_usec - baseTime.tv_usec;
  dt   = (float)sec + (float)usec / 1000000.0f;

#endif

  return dt;
}

//============================================================================
// Functions for STANDALONE mode
//============================================================================

#ifdef VV_STANDALONE

#ifdef _STANDARD_C_PLUS_PLUS
  #include <iostream>
  using std::cerr;
  using std::endl;
#else
  #include <iostream.h>
#endif

int main(int argc, char** argv)
{
  vvStopwatch* watch;
  char input[128];

  watch = new vvStopwatch();

  cerr << "Input something to start stopwatch: " << endl;
  cin >> input;
  watch->start();

  cerr << "Input something to stop watch: " << endl;
  cin >> input;
  watch->stop();
  cerr << "Current time: " << watch->getTime() << endl;

  cerr << "Input something to continue taking time: " << endl;
  cin >> input;
  watch->start();

  cerr << "Input something to stop watch: " << endl;
  cin >> input;
  watch->stop();
  cerr << "Total time: " << watch->getTime() << endl;

  delete watch;

  return 0;
}

#endif

//============================================================================
// End of File
//============================================================================
