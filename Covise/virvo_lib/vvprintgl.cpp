//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifdef WIN32
  #include <windows.h>
#else 
  #include <string.h>
  #include <stdarg.h>
#endif
#include <stdio.h>	
#include <GL/gl.h>
#ifdef VV_GLUT
  #include <glut.h>
#endif
#include "vvtoolshed.h"
#include "vvprintgl.h"
#include "vvdebugmsg.h"

/** Constructor.
  Builds the bitmap font.
  @param hDC Windows device context
*/
vvPrintGL::vvPrintGL()               
{
  vvDebugMsg::msg(3, "vvPrintGL::vvPrintGL()");

#ifdef VV_GLUT   

// Glut provides its own simple font handling which does not 
// need to be initialized explicitly.

#elif WIN32
  HFONT font;     // Windows Font ID
  HDC   hDC;      // Windows device context

  base = glGenLists(96);                // Storage For 96 Characters

  font = CreateFont(-24,                // Height Of Font
          0,                            // Width Of Font
          0,                            // Angle Of Escapement
          0,                            // Orientation Angle
          FW_BOLD,                      // Font Weight
          FALSE,                        // Italic
          FALSE,                        // Underline
          FALSE,                        // Strikeout
          ANSI_CHARSET,                 // Character Set Identifier
          OUT_TT_PRECIS,                // Output Precision
          CLIP_DEFAULT_PRECIS,          // Clipping Precision
          ANTIALIASED_QUALITY,          // Output Quality
          FF_DONTCARE|DEFAULT_PITCH,    // Family And Pitch
          "Courier New");               // Font Name

  hDC = wglGetCurrentDC();                                       
  SelectObject(hDC, font);              // Selects The Font We Want
                                        
  wglUseFontBitmaps(hDC, 32, 96, base); // Builds 96 Characters Starting At Character 32

#else
/*  TODO: make sure it works!
  Display     *dpy = __glutDisplay;
  XFontStruct *fontInfo;      // storage for our font. 

  base = glGenLists(96);      // storage for 96 characters. 

  // Load the font. What fonts you have is    
  // system dependent, but on my system they are     
  // in /usr/X11R6/lib/X11/fonts/, with fonts.alias and    
  // fonts.dir explaining what fonts the .pcf.gz files     
  // are. In any case, one of these 2 fonts should be     
  // on your system - or you won't see any text:

  fontInfo = XLoadQueryFont(dpy, "-adobe-helvetica-medium-r-normal--18-*-*-*-p-*-iso8859-1");
  if (fontInfo == NULL) 
  {
    fontInfo = XLoadQueryFont(dpy, "fixed");
    if (fontInfo == NULL) 
      cout << "no X font found" << endl;
  }

  // Start at character 32 (space), get 96 characters (a few characters past z), and 
  // store them starting at base:
  glXUseXFont(fontInfo->fid, 32, 96, base);

  // Free that font's info now since we've got the display lists:
  XFreeFont(dpy, fontInfo);
*/
#endif
}

/** Destructor.
  Deletes the bitmap font.
*/
vvPrintGL::~vvPrintGL()   
{
  vvDebugMsg::msg(3, "vvPrintGL::~vvPrintGL()");

#ifndef VV_GLUT
  glDeleteLists(base, 96);   // Delete All 96 Characters
#endif
}

//----------------------------------------------------------------------------
/** Print a text string to the current screen position.
  The current OpenGL drawing color is used.
  @param x,y text position [0..1]
  @param fmt printf compatible argument format
*/
void vvPrintGL::print(float x, float y, const char *fmt, ...)        
{
  vvDebugMsg::msg(3, "vvPrintGL::print()");

  va_list ap;                       // Pointer To List Of Arguments
  char    text[1024];               // Holds Our String
                                    
  if (fmt == NULL)                  // If There's No Text
    return;                         // Do Nothing
                                    
  va_start(ap, fmt);                // Parses The String For Variables
    vsprintf(text, fmt, ap);        // And Converts Symbols To Actual Numbers
  va_end(ap);                       // Results Are Stored In Text
  
  saveGLState();

#ifdef VV_GLUT

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
  glRasterPos2f(x, y);                  // set text position
  for (uint i=0; i<strlen(text); ++i)   // write characters one-by-one
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, text[i]);    

#else

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
  glRasterPos2f(x, y);                  // set text position
  glListBase(base - 32);            // Sets The Base Character to 32
  glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);  // Draws The Display List Text

#endif

  restoreGLState();
}

//----------------------------------------------------------------------------
/// Save GL state for text display.
void vvPrintGL::saveGLState()
{
  vvDebugMsg::msg(3, "vvPrintGL::saveGLState()");

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();

  glPushAttrib(GL_CURRENT_BIT | GL_TRANSFORM_BIT | GL_LIST_BIT);
}

//----------------------------------------------------------------------------
/// Restore GL state to previous state.
void vvPrintGL::restoreGLState()
{
  vvDebugMsg::msg(3, "vvPrintGL::restoreGLState()");

  glPopAttrib();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}



//============================================================================
// End of File
//============================================================================
