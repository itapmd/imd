//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifdef _STANDARD_C_PLUS_PLUS
  #include <iostream>
  using std::cerr;
  using std::endl;
#else
  #include <iostream.h>
#endif

#include "vvdebugmsg.h"

const char* vvDebugMsg::DEBUG_TEXT = "###Debug message: ";    ///< String printed before each debug message
vvDebugMsg::LevelType vvDebugMsg::debugLevel = NO_MESSAGES;     ///< Default debugging mode

/// Setter method for debugLevel
void vvDebugMsg::setDebugLevel(LevelType newLevel)
{
  debugLevel = newLevel;
  if (debugLevel!=NO_MESSAGES)
    cerr << "vvDebugMsg level set to: " << debugLevel << endl;
}

/// Setter method for debugLevel
void vvDebugMsg::setDebugLevel(int newLevel)

{
  debugLevel = (LevelType)newLevel;
  if (debugLevel!=NO_MESSAGES)
    cerr << "vvDebugMsg level set to: " << debugLevel << endl;
}

/// Getter method for debugLevel
vvDebugMsg::LevelType vvDebugMsg::getDebugLevel()
{
  return debugLevel;
}

/// Print an information string for debugging.
/// @param level level number of this message
/// @param text  information string being printed if current debug level is lower or
///              equal to given one
void vvDebugMsg::msg(int level, const char* text)
{
  if (level <= debugLevel)
    cerr << DEBUG_TEXT << text << endl; 
}

/// Print debug information and an integer.
/// @param level  level number of this message
/// @param text   information string being printed if current debug level is lower or
///               equal to given one
/// @param number integer string to be displayed after the message string
void vvDebugMsg::msg(int level, const char* text, int number)
{
  if (level <= debugLevel)
    cerr << DEBUG_TEXT << text << number << endl;
}

/// Print debug string and two integers.
/// @param level  level number of this message
/// @param text   information string being printed if current debug level is lower or
///               equal to given one
/// @param n1     first integer string to be displayed
/// @param n2     second integer string to be displayed
void vvDebugMsg::msg(int level, const char* text, int n1, int n2)
{
  if (level <= debugLevel)
    cerr << DEBUG_TEXT << text << n1 << ", " << n2 << endl;
}

/// Print debug string and three integers.
/// @param level  level number of this message
/// @param text   information string being printed if current debug level is lower or
///               equal to given one
/// @param n1     first integer string to be displayed
/// @param n2     second integer string to be displayed
/// @param n3     third integer string to be displayed
void vvDebugMsg::msg(int level, const char* text, int n1, int n2, int n3)
{
  if (level <= debugLevel)
    cerr << DEBUG_TEXT << text << n1 << ", " << n2 << ", " << n3 << endl;
}

/// Print debug string and four integers.
/// @param level  level number of this message
/// @param text   information string being printed if current debug level is lower or
///               equal to given one
/// @param n1     first integer string to be displayed
/// @param n2     second integer string to be displayed
/// @param n3     third integer string to be displayed
/// @param n4     fourth integer string to be displayed
void vvDebugMsg::msg(int level, const char* text, int n1, int n2, int n3, int n4)
{
  if (level <= debugLevel)
    cerr << DEBUG_TEXT << text << n1 << ", " << n2 << ", " << n3 << ", " << n4 << endl;
}

/// Print debug string and a float number.
/// @param level  level number of this message
/// @param text   information string being printed if current debug level is lower or
///               equal to given one
/// @param number float number to be displayed
void vvDebugMsg::msg(int level, const char* text, float number)
{
  if (level <= debugLevel)
    cerr << DEBUG_TEXT << text << number << endl;
}

/// Print debug string and two float numbers.
/// @param level  level number of this message
/// @param text   information string being printed if current debug level is lower or
///               equal to given one
/// @param n1     first float number to be displayed
/// @param n2     second float number to be displayed
void vvDebugMsg::msg(int level, const char* text, float n1, float n2)
{
  if (level <= debugLevel)
    cerr << DEBUG_TEXT << text << n1 << ", " << n2 << endl;
}

/// Print debug string and three float numbers.
/// @param level  level number of this message
/// @param text   information string being printed if current debug level is lower or
///               equal to given one
/// @param n1     first float number to be displayed
/// @param n2     second float number to be displayed
/// @param n3     third float number to be displayed
void vvDebugMsg::msg(int level, const char* text, float n1, float n2, float n3)
{
  if (level <= debugLevel)
    cerr << DEBUG_TEXT << text << n1 << ", " << n2 << ", " << n3 << endl;
}

/// Print debug string and four float numbers.
/// @param level  level number of this message
/// @param text   information string being printed if current debug level is lower or
///               equal to given one
/// @param n1     first float number to be displayed
/// @param n2     second float number to be displayed
/// @param n3     third float number to be displayed
/// @param n4     fourth float number to be displayed
void vvDebugMsg::msg(int level, const char* text, float n1, float n2, float n3, float n4)
{
  if (level <= debugLevel)
    cerr << DEBUG_TEXT << text << n1 << ", " << n2 << ", " << n3 << ", " << n4 << endl;
}

/// Print debug string and an additional string.
/// @param level  level number of this message
/// @param text   information string being printed if current debug level is lower or
///               equal to given one
/// @param str    additional string to be displayed
void vvDebugMsg::msg(int level, const char* text, const char* str)
{
  if (level <= debugLevel)
    cerr << DEBUG_TEXT << text << str << endl;
}

/// Check for active debug messages on the given debug level.
/// @param level debug level of this request
//  @return Return true if debug info shall be printed
bool vvDebugMsg::isActive(int level)
{
  if (level <= debugLevel) return true;
  else return false;
}

//****************************************************************************
// End of File
//****************************************************************************
