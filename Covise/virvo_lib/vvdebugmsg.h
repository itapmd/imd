//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifndef _VVDEBUGMSG_H_
#define _VVDEBUGMSG_H_

/** Debug message manager. 
    Allows the programmer to output debug messages only when required.
    @author Juergen Schulze-Doebold
*/
class vvDebugMsg
{
  public:
    /// Valid level types
    enum LevelType    
    {
      NO_MESSAGES   =  0,    ///< no messages are printed     
      FEW_MESSAGES  =  1,    ///< only the most important messages are printed
      MOST_MESSAGES =  2,    ///< all other messages which don't appear frequently are also printed
      ALL_MESSAGES  =  3     ///< also messages which appear frequently are printed    
    };

  private:
    static const char* DEBUG_TEXT;      ///< string to be printed at debug message
    static LevelType debugLevel;        ///< current debug level

  public:
    static void setDebugLevel(LevelType);
    static void setDebugLevel(int);
    static LevelType getDebugLevel();
    static void msg(int, const char*);
    static void msg(int, const char*, int);
    static void msg(int, const char*, int, int);
    static void msg(int, const char*, int, int, int);
    static void msg(int, const char*, int, int, int, int);
    static void msg(int, const char*, float);
    static void msg(int, const char*, float, float);
    static void msg(int, const char*, float, float, float);
    static void msg(int, const char*, float, float, float, float);
    static void msg(int, const char*, const char*);
    static bool isActive(int);
};

#endif

//****************************************************************************
// End of File
//****************************************************************************
