//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifndef _VVPRINTGL_H_
#define _VVPRINTGL_H_

#include <float.h>

//============================================================================
// Class Definition
//============================================================================

/** This class allows 2D billboard type text to be printed on an OpenGL
    canvas. There are separate implementations for use with and without
    GLUT. Define the variable VV_GLUT for the compiler to enable GLUT text.

    @author Juergen Schulze-Doebold
*/
class vvPrintGL
{
  private:
#ifndef VV_REMOTE_RENDERING
#ifndef VV_GLUT
    GLuint	base;				// Base Display List For The Font Set
#endif
    // GL state variables:
    GLint glsRasterPos[4];            ///< stores GL_CURRENT_RASTER_POSITION
    GLfloat glsColor[4];              ///< stores GL_CURRENT_COLOR
#endif    
    void saveGLState();
    void restoreGLState();

  public:
    vvPrintGL();
    virtual ~vvPrintGL();
    void print(float, float, const char *, ...);
};

#endif

//============================================================================
// End of File
//============================================================================
