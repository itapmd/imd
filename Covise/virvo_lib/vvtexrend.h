//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifndef _VVTEXREND_H_
#define _VVTEXREND_H_

#ifdef WIN32
  #include <windows.h>
  #include <GL/gl.h>
  #include <vvglext.h>
#elif __linux
  #include <GL/gl.h>
#else
  #include <GL/gl.h>
#endif
#include "vvvoldesc.h"
#include "vvrenderer.h"

//============================================================================
// Class Definitions
//============================================================================

/** Volume rendering engine using a texture-based algorithm.
  Textures can be drawn as planes or spheres. In planes mode a rendering
  quality can be given (determining the number of texture slices used), and
  the texture normal can be set according to the application's needs.<P>
  The data points are located at the grid as follows:<BR>
  The outermost data points reside at the very edge of the drawn object,
  the other values are evenly distributed inbetween.
  @author Juergen Schulze-Doebold (schulze@hlrs.de)
  @see vvRenderer
*/
class vvTexRend : public vvRenderer
{
  public:
    enum ErrorType      /// Error Codes
    {
      OK,               ///< no error
      TRAM_ERROR,       ///< not enough texture memory
      NO3DTEX           ///< 3D textures not supported on this hardware
    };
    enum RenderType     /// Rendering method
    {
      BEST,             ///< choose best according to display and data types
      SLICES2D,         ///< render slices parallel to xy axis plane using 2D textures
      CUBIC2D,          ///< render slices parallel to all axis planes using 2D textures
      CUBIC3D,          ///< render slices parallel to all axis planes using a 3D texture
      PLANAR3D,         ///< render planar slices using a 3D texture
      SPHERIC3D,        ///< render spheres originating at viewer using a 3D texture
      INDEX3D           ///< similar to PLANAR3D, but using indexed colors
    };

  private:
    enum AxisType               /// names for coordinate axes
    { X_AXIS, Y_AXIS, Z_AXIS }; 
    uchar rgbaConv[4096][4];    ///< density to RGBA conversion table (max. 12 bit density supported)
    int   texels[3];            ///< width, height and depth of volume, including empty space [texels]
    float texMin[3];            ///< minimum texture value of object [0..1] (to prevent border interpolation)
    float texMax[3];            ///< maximum texture value of object [0..1] (to prevent border interpolation)
    int   textures;             ///< number of textures stored in TRAM
    GLuint* texNames;           ///< names of texture slices stored in TRAM
    RenderType methodReq;       ///< rendering method requested by application
    RenderType methodUsed;      ///< rendering method currently used
    bool  extTex3d;             ///< true = 3D texturing supported
    bool  extColLUT;            ///< true = SGI texture color lookup table supported
    bool  extPalTex;            ///< true = OpenGL 1.2 paletted textures supported
    bool  extBlendEq;           ///< true = blending equation supported
		vvVector3 viewDir; 					///< user's current viewing direction [object coordinates]
		vvVector3 objDir; 	   			///< direction from viewer to object [object coordinates]
    bool interpolation;         ///< interpolation mode: true=linear interpolation (default), false=nearest neighbor

#ifdef WIN32
    PFNGLTEXIMAGE3DEXTPROC glTexImage3DEXT;
    PFNGLCOLORTABLESGIPROC glColorTableSGI;
    PFNGLCOLORTABLEEXTPROC glColorTableEXT;
    PFNGLBLENDEQUATIONEXTPROC glBlendEquationEXT;
#else
    typedef void (glTexImage3DEXT_type)(GLenum, GLint, GLenum, GLsizei, GLsizei, GLsizei, GLint, GLenum, GLenum, const GLvoid*);
    typedef void (glColorTableSGI_type)(GLenum, GLenum, GLsizei, GLenum, GLenum, const GLvoid*);
    typedef void (glColorTableEXT_type)(GLenum, GLenum, GLsizei, GLenum, GLenum, const GLvoid*);
    glTexImage3DEXT_type* glTexImage3DEXT;
    glColorTableSGI_type* glColorTableSGI;
    glColorTableEXT_type* glColorTableEXT;
#endif

    // GL state variables:
    GLboolean glsCulling;             ///< stores GL_CULL_FACE
    GLboolean glsBlend;               ///< stores GL_BLEND
    GLboolean glsColorMaterial;       ///< stores GL_COLOR_MATERIAL
    GLint glsBlendSrc;                ///< stores glBlendFunc(source,...)
    GLint glsBlendDst;                ///< stores glBlendFunc(...,destination)
    GLboolean glsTexColTable;         ///< stores GL_TEXTURE_COLOR_TABLE_SGI
    GLboolean glsLighting;            ///< stores GL_LIGHTING
    GLint glsMatrixMode;              ///< stores GL_MATRIX_MODE
    
    void removeTextures();
    void makeTextures();
    void makeTexturesSlices2D();
    ErrorType makeTexturesCubic2D();
    ErrorType makeTextures3D();
    void setGLenvironment();
    void unsetGLenvironment();
    void renderTex3DSpherical(vvMatrix4*);
    void renderTex3DPlanar(vvMatrix4*);   
    void renderTex3DFaces(AxisType, float, float, float);
    void renderTex2DSlices(float);
    void renderTex2DCubic(AxisType, float, float, float);
    RenderType findBestRenderingMethod();
    void makeColorLUT(float);
    int  getLUTSize();
    
  public:
    vvTexRend(vvVolDesc*, RenderType=BEST);
    virtual ~vvTexRend();
    void  renderVolumeGL();
    void  updateTransferFunction();
    void  updateVolumeData();
    void  activateClippingPlane();
    void  deactivateClippingPlane();
    void  setNumLights(int);
    void  setVoxelSize(vvVector3*);
    bool  instantClassification();
    void  setViewingDirection(vvVector3*);
    void  setObjectDirection(vvVector3*);
    void  setParameter(float, ParameterType, char* = NULL);
    float getParameter(ParameterType, char* = NULL);
    static bool isSupported(RenderType);
    RenderType getMethodUsed();
};

#endif

//============================================================================
// End of File
//============================================================================
