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

#ifdef WIN32
  #include <windows.h>
#endif

#include <string.h>
#include <GL/gl.h>
#include "vvgltools.h"

//============================================================================
// Method Definitions
//============================================================================

//----------------------------------------------------------------------------
/** Check OpenGL for errors.
    @return error string if there was an error, otherwise return NULL
*/
const char* vvGLTools::checkGLstate()
{
  GLenum e;           // error value
  const char* msg;    // error message
  
  e = glGetError();
  switch (e)
  {
    case GL_INVALID_ENUM:
      msg =  "GL_INVALID_ENUM";
      break;
    case GL_INVALID_VALUE:
      msg =  "GL_INVALID_VALUE";
      break;
    case GL_INVALID_OPERATION:
      msg =  "GL_INVALID_OPERATION";
      break;
    case GL_STACK_OVERFLOW:
      msg =  "GL_STACK_OVERFLOW";
      break;
    case GL_STACK_UNDERFLOW:
      msg =  "GL_STACK_UNDERFLOW";
      break;
    case GL_OUT_OF_MEMORY:
      msg =  "GL_OUT_OF_MEMORY";
      break;
#ifdef GL_SGI_texture_color_table
    case GL_TABLE_TOO_LARGE_EXT:
      msg =  "GL_TABLE_TOO_LARGE_EXT";
      break;
#elif defined(GL_TABLE_TOO_LARGE)
    case GL_TABLE_TOO_LARGE:
      msg =  "GL_TABLE_TOO_LARGE";
      break;
#endif
    case GL_NO_ERROR: 
      msg = NULL;
      break;
    default: 
      msg = "unknown"; 
      break;
  }
  return msg;
}

//----------------------------------------------------------------------------
/** Checks OpenGL for a specific extension.
    @param extension OpenGL extension to check for (e.g. "GL_EXT_bgra")
    @return true if extension is supported
*/
bool vvGLTools::isGLextensionSupported(const char* extension)
{
  const GLubyte *extensions = NULL;
  const GLubyte *start;
  GLubyte *where, *terminator;

  // Check requested extension name for existence and for spaces:
  where = (GLubyte*)strchr(extension, ' ');
  if (where || *extension=='\0') return false;

  // Get extensions string from OpenGL:
  extensions = glGetString(GL_EXTENSIONS);
  if (extensions=='\0') return false;
  
  // Parse OpenGL extensions string:
  start = extensions;
  for (;;) 
  {
    where = (GLubyte*)strstr((const char*)start, extension);
    if (!where) return false;
    terminator = where + strlen(extension);
    if (where==start || *(where - 1)==' ')
      if (*terminator==' ' || *terminator=='\0')
        return true;
    start = terminator;
  }
}

//----------------------------------------------------------------------------
/** Display the OpenGL extensions which are supported by the system at 
  run time.
  @param style display style
*/
void vvGLTools::displayOpenGLextensions(DisplayStyle style)
{
  char* extensions = NULL;    // OpenGL extensions string
  char* extCopy;              // local copy of extensions string for modifications
  int i;

  extensions = (char*)glGetString(GL_EXTENSIONS);

  switch (style)
  {
    default:
    case CONSECUTIVE: 
      cerr << extensions << endl;
      break;
    case ONE_BY_ONE: 
      extCopy = new char[strlen(extensions) + 1];
      strcpy(extCopy, extensions);
      for (i=0; i<(int)strlen(extCopy); ++i)
        if (extCopy[i] == ' ') extCopy[i] = '\n';
      cerr << extCopy << endl;
      delete[] extCopy;
      break;
  }
}

//----------------------------------------------------------------------------
/** Check for some specific OpenGL extensions.
  Displays the status of volume rendering related extensions, each on a separate line.
*/
void vvGLTools::checkOpenGLextensions()
{
  char* status[3] = {"supported", "not found"};

  cerr << "GL_EXT_texture3D.............";
  cerr << ((vvGLTools::isGLextensionSupported("GL_EXT_texture3D")) ? status[0] : status[1]) << endl;

  cerr << "GL_SGI_texture_color_table...";
  cerr << ((vvGLTools::isGLextensionSupported("GL_SGI_texture_color_table")) ? status[0] : status[1]) << endl;

  cerr << "GL_EXT_paletted_texture......";
  cerr << ((vvGLTools::isGLextensionSupported("GL_EXT_paletted_texture")) ? status[0] : status[1]) << endl;

  cerr << "GL_EXT_blend_equation........";
  cerr << ((vvGLTools::isGLextensionSupported("GL_EXT_blend_equation")) ? status[0] : status[1]) << endl;
}


//============================================================================
// End of File
//============================================================================
