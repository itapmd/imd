//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifndef _VVGLTOOLS_H_
#define _VVGLTOOLS_H_

//============================================================================
// Class Definitions
//============================================================================

/** Collection of OpenGL raleted tools. 
    Consists of static helper functions which are project independent.
    @author Juergen Schulze-Doebold
*/
class vvGLTools
{
  public:
    enum DisplayStyle           /// string display style for extensions display
    {
      CONSECUTIVE = 0,          ///< using entire line length
      ONE_BY_ONE  = 1           ///< one extension per line
    };
    static const char* checkGLstate();
    static bool isGLextensionSupported(const char*);
    static void displayOpenGLextensions(DisplayStyle);
    static void checkOpenGLextensions();
};

#endif

//============================================================================
// End of File
//============================================================================
