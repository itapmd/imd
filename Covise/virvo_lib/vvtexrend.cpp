//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifdef _STANDARD_C_PLUS_PLUS
  #include <iostream>
  #include <iomanip>
  using std::cerr;
  using std::endl;
  using std::setprecision;
#else
  #include <iostream.h>
  #include <iomanip.h>
#endif

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#ifdef WIN32
  #include <windows.h>
  #include <GL/gl.h>
  #include <vvglext.h>
#elif __linux
  #define GL_GLEXT_PROTOTYPES 1
  #include <GL/gl.h>
  #include "vvglext.h"
  #include <string.h>
  #include <values.h>
#else
  #include <GL/gl.h>
  #include "vvglext.h"
#endif
#ifndef WIN32
  #include <dlfcn.h>
  #include "vvdynlib.h"
#endif
#include "vvvecmath.h"
#include "vvdebugmsg.h"
#include "vvtoolshed.h"
#include "vvgltools.h"
#include "vvsphere.h"
#include "vvtexrend.h"
#include "vvstopwatch.h"

//----------------------------------------------------------------------------
/** Constructor.
  @param vd  volume description
  @param m   rendering method (default: choose best)
*/
vvTexRend::vvTexRend(vvVolDesc* vd, RenderType m) : vvRenderer(vd)
{
  vvDebugMsg::msg(1, "vvTexRend::vvTexRend()");

  rendererType = TEXREND;
  texNames = NULL;
  quality = 1.0f;
  viewDir.zero();
	objDir.zero();
  interpolation = true;

  // Find out which OpenGL extensions are supported:
  extTex3d    = vvGLTools::isGLextensionSupported("GL_EXT_texture3D");
  extColLUT   = vvGLTools::isGLextensionSupported("GL_SGI_texture_color_table");
  extPalTex   = vvGLTools::isGLextensionSupported("GL_EXT_paletted_texture");
  extBlendEq  = vvGLTools::isGLextensionSupported("GL_EXT_blend_equation");

#ifdef WIN32
  glTexImage3DEXT = (PFNGLTEXIMAGE3DEXTPROC)wglGetProcAddress("glTexImage3DEXT");
  if (glTexImage3DEXT==NULL) extTex3d = false;
  glColorTableSGI = (PFNGLCOLORTABLESGIPROC)wglGetProcAddress("glColorTableSGI");
  if (glColorTableSGI==NULL) extColLUT = false;
  glColorTableEXT = (PFNGLCOLORTABLEEXTPROC)wglGetProcAddress("glColorTableEXT");
  if (glColorTableEXT==NULL) extPalTex = false;
  glBlendEquationEXT = (PFNGLBLENDEQUATIONEXTPROC)wglGetProcAddress("glBlendEquationEXT");
  if (glBlendEquationEXT==NULL) extBlendEq = false;
#else
  VV_SHLIB_HANDLE dynLib = vvDynLib::open("libGL.so", RTLD_NOW);
  glTexImage3DEXT = (glTexImage3DEXT_type*)vvDynLib::sym(dynLib, "glTexImage3DEXT");
  if (glTexImage3DEXT==NULL) extTex3d = false;
  glColorTableSGI = (glColorTableSGI_type*)vvDynLib::sym(dynLib, "glColorTableSGI");
  if (glColorTableSGI==NULL) extColLUT = false;
  glColorTableEXT = (glColorTableEXT_type*)vvDynLib::sym(dynLib, "glColorTableEXT");
  if (glColorTableEXT==NULL) extPalTex = false;
  vvDynLib::close(dynLib);
#endif

  methodReq = m;   // save rendering method requested by application
  methodUsed = SLICES2D;  // default to simplest method
  vvDebugMsg::msg(1, "Requested rendering method: ", m);

  textures = 0;

  // Reset rendering method because a new volume is loaded:
  methodUsed = (methodReq==BEST) ? (findBestRenderingMethod()) : methodReq;

  // Determine volume size in world space:
  size.e[0] = vd->dist[0] * (float)vd->vox[0];
  size.e[1] = vd->dist[1] * (float)vd->vox[1];
  size.e[2] = vd->dist[2] * (float)vd->vox[2];

  updateTransferFunction();
  if (methodUsed==INDEX3D) makeTextures();	// wasn't done by updateTransferFunction()
}

//----------------------------------------------------------------------------
/// Destructor
vvTexRend::~vvTexRend()
{
  vvDebugMsg::msg(1, "vvTexRend::~vvTexRend()");
  removeTextures();
}

//----------------------------------------------------------------------------
/// Chooses the best rendering method depending on the graphics hardware's 
/// capabilities. 
vvTexRend::RenderType vvTexRend::findBestRenderingMethod()
{
  vvDebugMsg::msg(1, "vvTexRend::findBestRenderingMethod()");
  
  if (extTex3d)
  {
    if (vd->bpv<3 && (extColLUT || extPalTex)) return INDEX3D;
    else return PLANAR3D;
  }
  else return SLICES2D;
}

//----------------------------------------------------------------------------
/// Remove all textures from texture memory.
void vvTexRend::removeTextures()
{
  vvDebugMsg::msg(1, "vvTexRend::removeTextures()");

  if (textures>0)
  {
    glDeleteTextures(textures, texNames);
    delete[] texNames;
    texNames = NULL;
    textures = 0;
  }
}

//----------------------------------------------------------------------------
/// Generate textures.
void vvTexRend::makeTextures()
{
  bool done = false;

  vvDebugMsg::msg(2, "vvTexRend::makeTextures()");

  // Compute texture dimensions (must be power of 2):
  texels[0]  = vvToolshed::getTextureSize(vd->vox[0]);
  texels[1] = vvToolshed::getTextureSize(vd->vox[1]);
  texels[2] = vvToolshed::getTextureSize(vd->vox[2]);

  while (!done)
  {
    switch (methodUsed)
    {
      default:
      case SLICES2D:  
        makeTexturesSlices2D();
        done = true;
        break;

      case CUBIC2D:   
        if (makeTexturesCubic2D() != OK && methodReq == BEST)
          methodUsed = SLICES2D;
        else
          done = true;
        break;

      case CUBIC3D:   
      case PLANAR3D:  
      case SPHERIC3D: 
        if (makeTextures3D() != OK && methodReq == BEST)
          methodUsed = CUBIC2D;
        else
          done = true;
        break;

      case INDEX3D:
        if (makeTextures3D() != OK && methodReq == BEST)
          methodUsed = CUBIC2D;
        else
          done = true;
        break;
    }
  }
}

//----------------------------------------------------------------------------
/// Generate 2D textures for slices 2D mode.
void vvTexRend::makeTexturesSlices2D()
{
  GLint glWidth;      // return value from OpenGL call
  uchar* rgbaSlice;   // RGBA slice for texture memory
  int rawVal[4];      // raw values for R,G,B,A
  int sliceSize;      // number of voxels in a raw slice
  int texSize;        // texture size in bytes
  int sliceOffset, heightOffset;  // offsets in raw volume data
  int frames;         // sequence timesteps
  int index, srcIndex;
  int s, x, y, f, c;
  uchar* raw;         // raw volume data
  bool accommodated = true;  // false if a texture cannot be accommodated in TRAM

  vvDebugMsg::msg(1, "vvTexRend::makeTexturesSlices2D()");

  removeTextures();   // first remove previously generated textures from TRAM
  
  frames = vd->frames;

  textures  = vd->vox[2] * frames;
  texSize   = texels[0] * texels[1] * 4;   // 4 for RGBA texels
  
  vvDebugMsg::msg(1, "Texture width  = ", texels[0]);
  vvDebugMsg::msg(1, "Texture height = ", texels[1]);
  vvDebugMsg::msg(1, "Texture slices = ", vd->vox[2]);
  vvDebugMsg::msg(1, "Timesteps      = ", frames);
  
  // Generate texture data:
  rgbaSlice = new uchar[texSize];   
  memset(rgbaSlice, 0, texSize);   // initialize with 0's for invisible borders

  // Generate texture names:
  assert(texNames==NULL);
  texNames = new GLuint[textures];
  glGenTextures(textures, texNames);  

  // Generate texture contents:
  sliceSize = vd->getSliceSize();   // buffer for speed
  vvDebugMsg::msg(2, "Transferring 2D textures to TRAM. Total size [KB]: ", 
    frames * vd->vox[2] * texSize / 1024);
  for (f=0; f<frames; ++f)
  {
    raw = vd->getRaw(f);
    for (s=0; s<vd->vox[2]; ++s)
    {
      sliceOffset = s * sliceSize;
      for (y=0; y<vd->vox[1]; ++y)
      {
        heightOffset = y * vd->vox[0] * vd->bpv;
        switch (vd->bpv)
        {
          case 1: // 1 bpv, interpreted as 8 bit density value
          case 2: // 2 bpv, interpreted as 12 bit density value
            for (x=0; x<vd->vox[0]; ++x)
            {
              if (vd->bpv == 1)
                rawVal[0] = raw[x + sliceOffset + heightOffset];
              else
              {
                srcIndex = (x+x) + sliceOffset + heightOffset;
                rawVal[0] = ((int)raw[srcIndex] << 8) | (int)raw[srcIndex + 1];
                rawVal[0] >>= 4;    // make 12 bit index into LUT
              }
              index = x + y * texels[0];
              rgbaSlice[4 * index    ] = rgbaConv[rawVal[0]][0];
              rgbaSlice[4 * index + 1] = rgbaConv[rawVal[0]][1];
              rgbaSlice[4 * index + 2] = rgbaConv[rawVal[0]][2];
              rgbaSlice[4 * index + 3] = rgbaConv[rawVal[0]][3];
            }
            break;

          case 3: // 3 bpv, interpreted as RGB color
					case 4: // 4 bpv, interpreted as RGB and density
            for (x=0; x<vd->vox[0]; ++x)
            {
              index = x + y * texels[0];
              for (c=0; c<vd->bpv; ++c)
                rawVal[c] = raw[sliceOffset + vd->bpv * x + heightOffset + c];

              for (c=0; c<3; ++c)
                rgbaSlice[4 * index + c] = (uchar)rawVal[c];

							if (vd->bpv==3)
  	            rgbaSlice[4 * index + 3] = (uchar)(((int)rgbaConv[rawVal[0]][0] + // Alpha is mean of sum of RGB conversion table results
    	                                      (int)rgbaConv[rawVal[1]][1] + 
      	                                    (int)rgbaConv[rawVal[2]][2]) / 3);
							else		// if bpv==4
								rgbaSlice[4 * index + 3] = rgbaConv[rawVal[3]][3];
            }
            break;
						
          default: break;
        }
      }
      glBindTexture(GL_TEXTURE_2D, texNames[vd->vox[2] * f + s]);    
      glPixelStorei(GL_UNPACK_ALIGNMENT,1);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, (interpolation) ? GL_LINEAR : GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, (interpolation) ? GL_LINEAR : GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
      
      // Only load texture if it can be accommodated:
      glTexImage2D(GL_PROXY_TEXTURE_2D, 0, GL_RGBA,
        texels[0], texels[1], 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
      glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &glWidth); 
      if (glWidth!=0) 
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texels[0], texels[1], 
          0, GL_RGBA, GL_UNSIGNED_BYTE, rgbaSlice);
      else
        accommodated = false;
    }
  }

  if (accommodated==false)
    cerr << "Insufficient texture memory for 2D textures." << endl;

  delete[] rgbaSlice;
}

//----------------------------------------------------------------------------
/// Generate 2D textures for cubic 2D mode.
vvTexRend::ErrorType vvTexRend::makeTexturesCubic2D()
{
  const int TEXELSIZE = 4;  // bytes needed for each texel
  GLint glWidth;            // return value from OpenGL call
  uchar* rgbaSlice[3];      // RGBA slice data for texture memory for each principal axis
  int rawVal[4];            // raw values for R,G,B,A
  int rawSliceSize;         // number of bytes in a slice of the raw data array
  int rawLineSize;          // number of bytes in a row of the raw data array
  int texSize[3];           // size of a 2D texture in bytes for each principal axis
  int frames;               // sequence timesteps
  uchar* raw;               // raw volume data
  int texIndex=0;           // index of current texture
  int texSliceIndex;        // index of current voxel in texture
  int tw[3], th[3];         // current texture width and height for each principal axis
  int rw[3], rh[3], rs[3];  // raw data width, height, slices for each principal axis
  int rawStart[3];          // starting offset into raw data, for each principal axis
  int rawStepW[3];          // raw data step size for texture row, for each principal axis
  int rawStepH[3];          // raw data step size for texture column, for each principal axis
  int rawStepS[3];          // raw data step size for texture slices, for each principal axis
  uchar* rawVoxel;          // current raw data voxel
  int i, s, w, h, f, c;      
  bool accommodated = true;  // false if a texture cannot be accommodated in TRAM

  vvDebugMsg::msg(1, "vvTexRend::makeTexturesCubic2D()");

  removeTextures();   // first remove previously generated textures from TRAM
  
  frames = vd->frames;
  textures = (vd->vox[0] + vd->vox[1] + vd->vox[2]) * frames;   // total number of textures
  
  vvDebugMsg::msg(1, "Total number of 2D textures:    ", textures);
  vvDebugMsg::msg(1, "Total size of 2D textures [KB]: ", 
    frames * 3 * texels[0] * texels[1] * texels[2] / 1024);
  
  // Generate texture names:
  assert(texNames==NULL);
  texNames = new GLuint[textures];
  glGenTextures(textures, texNames);  

  // Initialize texture sizes:
  th[1] = tw[2] = texels[0];
  tw[0] = th[2] = texels[1];
  tw[1] = th[0] = texels[2];
  for (i=0; i<3; ++i)
    texSize[i] = tw[i] * th[i] * TEXELSIZE; 

  // Initialize raw data sizes:
  rs[0] = rh[1] = rw[2] = vd->vox[0];
  rw[0] = rs[1] = rh[2] = vd->vox[1];
  rh[0] = rw[1] = rs[2] = vd->vox[2];
  rawSliceSize = vd->getSliceSize();
  rawLineSize  = vd->vox[0] * vd->bpv;

  // Initialize raw data access info:
  rawStart[0] = (vd->vox[2]) * rawSliceSize - vd->bpv; 
  rawStepW[0] = -rawLineSize;
  rawStepH[0] = -rawSliceSize;
  rawStepS[0] = -vd->bpv;
  rawStart[1] = (vd->vox[2] - 1) * rawSliceSize; 
  rawStepW[1] = -rawSliceSize;
  rawStepH[1] = vd->bpv;
  rawStepS[1] = rawLineSize;
  rawStart[2] = (vd->vox[1] - 1) * rawLineSize;; 
  rawStepW[2] = vd->bpv;
  rawStepH[2] = -rawLineSize;
  rawStepS[2] = rawSliceSize;

  // Generate texture data arrays:
  for (i=0; i<3; ++i)
  {
    rgbaSlice[i] = new uchar[texSize[i]];   
  }

  // Generate texture data:
  for (f=0; f<frames; ++f)
  {
    raw = vd->getRaw(f);  // points to beginning of frame in raw data
    for (i=0; i<3; ++i)   // generate textures for each principal axis
    {
      memset(rgbaSlice[i], 0, texSize[i]);   // initialize with 0's for invisible empty regions
  
      // Generate texture contents:
      for (s=0; s<rs[i]; ++s)   // loop thru texture and raw data slices
      {
        for (h=0; h<rh[i]; ++h)   // loop thru raw data rows
        {
          // Set voxel to starting position in raw data array:
          rawVoxel = raw + rawStart[i] + s * rawStepS[i] + h * rawStepH[i]; 

          for (w=0; w<rw[i]; ++w)   // loop thru raw data columns
          {
            texSliceIndex = TEXELSIZE * (w + h * tw[i]);
            switch (vd->bpv)
            {
              case 1: // 1 bpv, interpreted as 8 bit density value
              case 2: // 2 bpv, interpreted as 12 bit density value
                if (vd->bpv == 1)
                  rawVal[0] = int(*rawVoxel);
                else
                {
                  rawVal[0] = (int(*rawVoxel) << 8) | int(*(rawVoxel+1));
                  rawVal[0] >>= 4;    // make 12 bit LUT index
                }
                for (c=0; c<TEXELSIZE; ++c)
                  rgbaSlice[i][texSliceIndex + c] = rgbaConv[rawVal[0]][c];
                break;
						    
              case 3: // 3 bpv, interpreted as RGB color
					    case 4: // 4 bpv, interpreted as RGB and density
                for (c=0; c<vd->bpv; ++c)   // fetch raw value from memory
                  rawVal[c] = *(rawVoxel + c);

                for (c=0; c<3; ++c) // duplicate RGB to texture memory
                  rgbaSlice[i][texSliceIndex + c] = (uchar)rawVal[c];

							  if (vd->bpv==3) // if volume is RGB: compute alpha from RGB
  	              rgbaSlice[i][texSliceIndex + 3] = 
                    (uchar)(((int)rgbaConv[rawVal[0]][0] + // Alpha is mean of sum of RGB conversion table results
            	               (int)rgbaConv[rawVal[1]][1] + 
      	                     (int)rgbaConv[rawVal[2]][2]) / 3);
							  else		// if volume is RGBA: lookup alpha in LUT
								  rgbaSlice[i][texSliceIndex + 3] = rgbaConv[rawVal[3]][3];
                break;
						    
              default: break;
            }
            rawVoxel += rawStepW[i];
          }
        }
        glBindTexture(GL_TEXTURE_2D, texNames[texIndex]);    
        glPixelStorei(GL_UNPACK_ALIGNMENT,1);
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, (interpolation) ? GL_LINEAR : GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, (interpolation) ? GL_LINEAR : GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

        // Only load texture if it can be accommodated:
        glTexImage2D(GL_PROXY_TEXTURE_2D, 0, GL_RGBA,
          tw[i], th[i], 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
        glGetTexLevelParameteriv(GL_PROXY_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &glWidth); 
        if (glWidth!=0) 
          glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, tw[i], th[i], 
            0, GL_RGBA, GL_UNSIGNED_BYTE, rgbaSlice[i]);
        else
          accommodated = false;

        ++texIndex;
      }
    }
  }

  if (accommodated==false)
    cerr << "Insufficient texture memory for cubic 2D textures." << endl;

  assert(texIndex == textures);  
  for (i=0; i<3; ++i)
    delete[] rgbaSlice[i];
  return OK;
}

//----------------------------------------------------------------------------
/// Generate 3D textures.
vvTexRend::ErrorType vvTexRend::makeTextures3D()
{
  ErrorType err = OK;
  GLint glWidth;          // return value from OpenGL call
  uchar* texData;         // data for texture memory
  uchar* raw;             // raw volume data
  int rawVal[4];          // raw values for R,G,B,A
  int sliceSize;          // number of voxels in a raw slice
  int texLineOffset;      // index into currently processed line of texture data array [texels]
  int texOffset;          // currently processed texel [texels]
  int texSize;            // texture size in bytes
  int frames;             // number of time steps in sequence
  int rawSliceOffset;     // slice offset in raw data
  int heightOffset, srcIndex;
  int s, x, y, f, c;      
  bool accommodated = true;  // false if a texture cannot be accommodated in TRAM

  vvDebugMsg::msg(1, "vvTexRend::makeTextures3D()");

  if (!extTex3d) return NO3DTEX;

  removeTextures();   // first remove previously generated textures from TRAM
  frames = vd->frames;

  // Compute texture dimensions (must be power of 2):
  texels[0] = vvToolshed::getTextureSize(vd->vox[0]);
  texels[1] = vvToolshed::getTextureSize(vd->vox[1]);
  texels[2] = vvToolshed::getTextureSize(vd->vox[2]);
  textures  = frames;
  if (methodUsed==INDEX3D)
  {
    texSize = texels[0] * texels[1] * texels[2] * 2;
    vvDebugMsg::msg(2, "Using paletted textures.");
  }
  else
  {
    texSize = texels[0] * texels[1] * texels[2] * 4;
    vvDebugMsg::msg(2, "Using RGBA textures.");
  }
  
  vvDebugMsg::msg(1, "3D Texture width     = ", texels[0]);
  vvDebugMsg::msg(1, "3D Texture height    = ", texels[1]);
  vvDebugMsg::msg(1, "3D Texture depth     = ", texels[2]);
  vvDebugMsg::msg(1, "3D Texture size (KB) = ", texSize / 1024);

  texData = new uchar[texSize];   // R+G+B+A = 4
  memset(texData, 0, texSize);    // initialize with 0's

  sliceSize = vd->getSliceSize();   // buffer for speed

  vvDebugMsg::msg(2, "Creating texture names. # of names: ", frames);

  assert(texNames==NULL);
  texNames = new GLuint[frames];
  glGenTextures(frames, texNames);  // generate texture name

  // Generate texture contents:
  vvDebugMsg::msg(2, "Transferring textures to TRAM. Total size [KB]: ", 
    frames * texSize / 1024);
  for (f=0; f<frames; ++f)
  {
    raw = vd->getRaw(f);
    for (s=0; s<vd->vox[2]; ++s)
    {
      rawSliceOffset = s * sliceSize;
      for (y=0; y<vd->vox[1]; ++y)
      {
        heightOffset = y * vd->vox[0] * vd->bpv;
        texLineOffset = (vd->vox[1]-1-y) * texels[0] + (vd->vox[2]-1-s) * texels[0] * texels[1];
        switch (vd->bpv)
        {
          case 1: // 1 bpv, interpreted as 8 bit density value
          case 2: // 2 bpv, interpreted as 12 bit density value
            for (x=0; x<vd->vox[0]; ++x)
            {
              srcIndex = vd->bpv * x + rawSliceOffset + heightOffset;
              if (vd->bpv == 1)
                rawVal[0] = (int)raw[srcIndex];
              else
              {
                rawVal[0] = ((int)raw[srcIndex] << 8) | (int)raw[srcIndex + 1];
                rawVal[0] >>= 4;    // make 12 bit index into LUT
              }
              texOffset = x + texLineOffset;
              if (methodUsed==INDEX3D)
              {
                if (vd->bpv==2) rawVal[0] >>= 4;  // 12 bit mode?
                if (extColLUT)
                  texData[2*texOffset] = texData[2*texOffset+1] = (uchar)rawVal[0];
                else  // extPalTex
                  texData[texOffset] = (uchar)rawVal[0];
              }
              else
                for (c=0; c<4; ++c)
                  texData[4 * texOffset + c] = rgbaConv[rawVal[0]][c];
            }
            break;
          case 3: // 3 bpv, interpreted as RGB color
					case 4: // 4 bpv, interpreted as RGB and density
            if (methodUsed!=INDEX3D)    // RGB mode -> no palette
            {
              for (x=0; x<vd->vox[0]; ++x)
              {
                texOffset = x + texLineOffset;
                for (c=0; c<vd->bpv; ++c)
                  rawVal[c] = (int)raw[rawSliceOffset + heightOffset + x * vd->bpv + c];

	              for (c=0; c<3; ++c)
                  texData[4 * texOffset + c] = (uchar)rawVal[c];
									
								if (vd->bpv==3)
      	          texData[4 * texOffset + 3] = (uchar)(((int)rgbaConv[rawVal[0]][0] + 
        	                                  (int)rgbaConv[rawVal[1]][1] + 
          	                                (int)rgbaConv[rawVal[2]][2]) / 3);
								else		// if bpv==4
									texData[4 * texOffset + 3] = rgbaConv[rawVal[3]][3];
              }
            }
            break;
          default: break;
        }
      }
    }

    if (methodUsed == SPHERIC3D)
    {
      // Set edge values to 0 for spheric textures, because textures
      // may exceed texel volume:
      for (s=0; s<vd->vox[2]; ++s) 
      {
        for (y=0; y<vd->vox[1]; ++y) 
        {
          for (x=0; x<vd->vox[0]; ++x) 
          {
	          if ((s == 0) || (s==vd->vox[2]-1) ||
	              (y == 0) || (y==vd->vox[1]-1) ||
	              (x == 0) || (x==vd->vox[0]-1)) 
            {
	            texOffset = x + y * texels[0] + s * texels[0] * texels[1];
              if (methodUsed==INDEX3D)
              {
                if (extColLUT)
                  texData[2 * texOffset] = texData[2 * texOffset + 1] = 0;
                else
                  texData[texOffset] = 0;
              }
              else
              {
	              texData[4 * texOffset    ] = texData[4 * texOffset + 1] =
	              texData[4 * texOffset + 2] = texData[4 * texOffset + 3] = 0;
              }
	          }
          }
        }
      }
    }
    
    glBindTexture(GL_TEXTURE_3D_EXT, texNames[f]);    
    glPixelStorei(GL_UNPACK_ALIGNMENT,1);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameteri(GL_TEXTURE_3D_EXT, GL_TEXTURE_MAG_FILTER, (interpolation) ? GL_LINEAR : GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D_EXT, GL_TEXTURE_MIN_FILTER, (interpolation) ? GL_LINEAR : GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D_EXT, GL_TEXTURE_WRAP_R_EXT, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_3D_EXT, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_3D_EXT, GL_TEXTURE_WRAP_T, GL_CLAMP);

    if (methodUsed==INDEX3D)
    {
      // Only load texture if it can be accommodated:
      if (extColLUT)
        glTexImage3DEXT(GL_PROXY_TEXTURE_3D_EXT, 0, GL_LUMINANCE_ALPHA,
          texels[0], texels[1], texels[2], 0, GL_LUMINANCE_ALPHA, GL_UNSIGNED_BYTE, NULL);
      else  // for extPalTex
        glTexImage3DEXT(GL_PROXY_TEXTURE_3D_EXT, 0, GL_COLOR_INDEX8_EXT,
          texels[0], texels[1], texels[2], 0, GL_COLOR_INDEX , GL_UNSIGNED_BYTE, NULL);        

      glGetTexLevelParameteriv(GL_PROXY_TEXTURE_3D_EXT, 0, GL_TEXTURE_WIDTH, &glWidth); 
      if (glWidth!=0) 
      {
        if (extColLUT)
          glTexImage3DEXT(GL_TEXTURE_3D_EXT, 0, GL_LUMINANCE_ALPHA, texels[0], texels[1], texels[2], 0,
                          GL_LUMINANCE_ALPHA, GL_UNSIGNED_BYTE, texData);
        else  // for extPalTex
          glTexImage3DEXT(GL_TEXTURE_3D_EXT, 0, GL_COLOR_INDEX8_EXT, texels[0], texels[1], texels[2], 0,
                          GL_COLOR_INDEX, GL_UNSIGNED_BYTE, texData);
      }
      else
        accommodated = false;
    }
    else
    {
      // Only load texture if it can be accommodated:
      glTexImage3DEXT(GL_PROXY_TEXTURE_3D_EXT, 0, GL_RGBA,
        texels[0], texels[1], texels[2], 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
      glGetTexLevelParameteriv(GL_PROXY_TEXTURE_3D_EXT, 0, GL_TEXTURE_WIDTH, &glWidth); 
      if (glWidth!=0) 
        glTexImage3DEXT(GL_TEXTURE_3D_EXT, 0, GL_RGBA, texels[0], texels[1], texels[2], 0,
                        GL_RGBA, GL_UNSIGNED_BYTE, texData);
      else
        accommodated = false;
    }
  }

  if (accommodated==false)
  {
    cerr << "Insufficient texture memory for 3D texture(s)." << endl;
    err = TRAM_ERROR;
  }

  delete[] texData;
  return err;
}

//----------------------------------------------------------------------------
/// Update transfer function from volume description.
void vvTexRend::updateTransferFunction()
{
  int i, c;
  float* rgba;
  int lutEntries;

  vvDebugMsg::msg(1, "vvTexRend::updateTransferFunction()");

  lutEntries = getLUTSize();
  rgba = new float[4 * lutEntries];

  // Generate arrays from pins:
  vd->tf.pins.setColorModel(vvPinList::RGB_MODEL);
  vd->tf.pins.computeFunction(lutEntries, vvPinList::RGBA, rgba);

  // Copy RGBA values to internal array:
  for (i=0; i<lutEntries; ++i)
		for (c=0; c<4; ++c)
    {
      if (methodUsed==INDEX3D && extPalTex)
	      rgbaConv[i][c] = (uchar)(rgba[c * lutEntries + i] * 254.0f);    // TODO: really 254???
      else
        rgbaConv[i][c] = (uchar)(rgba[c * lutEntries + i] * 255.0f);
    }
  delete[] rgba;

  if (methodUsed==INDEX3D)  
    makeColorLUT(1.0f); // generate color lookup table
  else 
    makeTextures();  // if no color table is used, textures must be reloaded
}
  
//----------------------------------------------------------------------------
// see parent in vvRenderer
void vvTexRend::updateVolumeData()
{
  makeTextures();
}

//----------------------------------------------------------------------------
/// Set GL environment for texture rendering.
void vvTexRend::setGLenvironment()
{
  vvDebugMsg::msg(3, "vvTexRend::setGLenvironment()");

  // Save current GL state:
  glGetBooleanv(GL_CULL_FACE, &glsCulling);
  glGetBooleanv(GL_BLEND, &glsBlend);
  glGetBooleanv(GL_COLOR_MATERIAL, &glsColorMaterial);
  glGetIntegerv(GL_BLEND_SRC, &glsBlendSrc);
  glGetIntegerv(GL_BLEND_DST, &glsBlendDst);
  glGetBooleanv(GL_LIGHTING, &glsLighting);
  glGetIntegerv(GL_MATRIX_MODE, &glsMatrixMode);
  if (extColLUT)
    glGetBooleanv(GL_TEXTURE_COLOR_TABLE_SGI, &glsTexColTable);

  // Set new GL state:
  glDisable(GL_CULL_FACE);       
	glDisable(GL_LIGHTING); 
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_BLEND);
  if (extColLUT)
    glDisable(GL_TEXTURE_COLOR_TABLE_SGI);  
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glMatrixMode(GL_TEXTURE);
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);

  bool mip = false;
  if (extBlendEq && mip)
    glBlendEquationEXT(GL_MAX_EXT);   // set maximum intensity projection
   
  vvDebugMsg::msg(3, "vvTexRend::setGLenvironment() done");
}

//----------------------------------------------------------------------------
/// Unset GL environment for texture rendering.
void vvTexRend::unsetGLenvironment()
{
  vvDebugMsg::msg(3, "vvTexRend::unsetGLenvironment()");

  if (glsCulling==(uchar)true) glEnable(GL_CULL_FACE);
  else glDisable(GL_CULL_FACE);

  if (glsBlend==(uchar)true) glEnable(GL_BLEND);
  else glDisable(GL_BLEND);

  if (glsColorMaterial==(uchar)true) glEnable(GL_COLOR_MATERIAL);
  else glDisable(GL_COLOR_MATERIAL);

  if (extColLUT)
  {
    if (glsTexColTable==(uchar)true) glEnable(GL_TEXTURE_COLOR_TABLE_SGI);
    else glDisable(GL_TEXTURE_COLOR_TABLE_SGI);
  }
  
  if (glsLighting==(uchar)true) glEnable(GL_LIGHTING);
  else glDisable(GL_LIGHTING);
  
  glBlendFunc(glsBlendSrc, glsBlendDst);
  glMatrixMode(glsMatrixMode);
  vvDebugMsg::msg(3, "vvTexRend::unsetGLenvironment() done");
}

//----------------------------------------------------------------------------
/** Render a volume entirely if probeSize=0 or a cubic sub-volume of size probeSize.
  @param mv        model-view matrix
*/
void vvTexRend::renderTex3DPlanar(vvMatrix4* mv)
{
  vvMatrix4 invMV;        	// inverse of model-view matrix
  vvMatrix4 pm;             // OpenGL projection matrix
  vvVector3 size2;        	// half object size
  vvVector3 objPos;       	// object midpoint position in world space
  vvVector3 p1, p2;       	// points on volume edges in textures space
  vvVector3 isect[6];     	// intersection points, maximum of 6 allowed when intersecting a plane and a volume [object space]
  vvVector3 texcoord[6];  	// intersection points in texture coordinate space [0..1]
  vvVector3 farthest;     	// volume vertex farthest from the viewer
  vvVector3 delta;        	// distance vector between textures [object space]
  vvVector3 normal;       	// normal vector of textures
  vvVector3 temp;         	// temporary vector
  vvVector3 origin;       	// origin (0|0|0) transformed to object space
  vvVector3 eye;          	// user's eye position [object space]
  vvVector3 normClipPoint;	// normalized point on clipping plane
  vvVector3 clipPosObj; 		// clipping plane position in object space w/o position
  vvVector3 probePosObj;    // probe midpoint [object space]
  vvVector3 probeSizeObj;   // probe size [object space]
  vvVector3 probeTexels;    // number of texels in each probe dimension
  vvVector3 probeMin, probeMax;  // probe min and max coordinates [object space]
  vvVector3 texSize;        // size of 3D texture [object space]
  vvVector3 pos;            // volume location
  float     diagonal;  	  	// probe diagonal [object space]
	float     maxSize;  			// maximum volume edge length [object space]
  float     maxDist;      	// maximum length of texture drawing path
  int       optNumSlices;   // optimum number of texture slices for 1:1 slices to sample rate
  int       numSlices;     	// number of texture slices drawn
  int       i;              // general counter

  vvDebugMsg::msg(3, "vvTexRend::renderTex3DPlanar()");
	
  if (!extTex3d) return;	// needs 3D texturing extension

  // Determine texture object dimensions and half object size as a shortcut:
  for (i=0; i<3; ++i)
  {
    texSize.e[i] = size.e[i] * (float)texels[i] / (float)vd->vox[i];
    size2.e[i]   = 0.5f * size.e[i];
  }
  pos.copy(&vd->position);

  // Calculate inverted modelview matrix:
  invMV.copy(mv);
  invMV.invert();

  // Find eye position:
  getEyePosition(&eye);
  eye.multiply(&invMV);

	if (probeSize > 0.0f)
	{
	  // Convert probe midpoint coordinates to object space w/o position:
	  probePosObj.copy(&probePos);
	  probePosObj.sub(&pos); 		// eliminate object position from probe position

  	// Compute probe min/max coordinates in object space:
		maxSize = ts_max(size[0], size[1], size[2]);
  	for (i=0; i<3; ++i)
  	{
    	probeMin[i] = probePosObj[i] - (probeSize * maxSize) / 2.0f;
    	probeMax[i] = probePosObj[i] + (probeSize * maxSize) / 2.0f;
  	}

  	// Constrain probe boundaries to volume data area:
  	for (i=0; i<3; ++i)
  	{
    	if (probeMin[i] > size2[i] || probeMax[i] < -size2[i]) 
    	{
				vvDebugMsg::msg(3, "probe outside of volume");
      	return;
    	}
    	if (probeMin[i] < -size2[i]) probeMin[i] = -size2[i];
    	if (probeMax[i] >  size2[i]) probeMax[i] =  size2[i];
			probePosObj[i] = (probeMax[i] + probeMin[i]) / 2.0f;
  	}

		// Compute probe edge lengths:
  	for (i=0; i<3; ++i)   
    	probeSizeObj[i] = probeMax[i] - probeMin[i];
	}
	else		// probe mode off
	{
		probeSizeObj.copy(&size);		
		probeMin.set(-size2[0], -size2[1], -size2[2]);
		probeMax.copy(&size2);
		probePosObj.zero();
	}
  
  // Compute length of probe diagonal [object space]:
  diagonal = (float)sqrt(
    probeSizeObj[0] * probeSizeObj[0] +
    probeSizeObj[1] * probeSizeObj[1] +
    probeSizeObj[2] * probeSizeObj[2]);

  // Initialize texture counters: draw twice as many slices as there are samples
  // in the diagonal.
	if (probeSize>0.0f)
	{
  	probeTexels.zero();
    for (i=0; i<3; ++i)
	    probeTexels[i] = texels[i] * probeSizeObj[i] / texSize.e[i];
	}
	else		// probe mode off
	{
		probeTexels.set((float)vd->vox[0], (float)vd->vox[1], (float)vd->vox[2]);
	}

	optNumSlices = (int)(sqrt((double)(
    probeTexels[0] * probeTexels[0] + 
    probeTexels[1] * probeTexels[1] + 
		probeTexels[2] * probeTexels[2])));  // probe diagonal
	numSlices = (int)(quality * optNumSlices);
  if (numSlices < 1)  // make sure that at least one slice is drawn
    numSlices = 1;

  vvDebugMsg::msg(3, "Number of textures rendered: ", numSlices);

  // Use alpha correction in indexed mode: adapt alpha values to number of textures:
  if (instantClassification())
  {
    float diagonalVoxels;
    diagonalVoxels = sqrtf(float(vd->vox[0] * vd->vox[0] + 
                                 vd->vox[1] * vd->vox[1] + 
                                 vd->vox[2] * vd->vox[2]));
    makeColorLUT(diagonalVoxels / float(numSlices));
  }

  // Get projection matrix:
  getProjectionMatrix(&pm);

  // Compute normal vector of textures using the follwing strategy:
  // For orthographic projections or if viewDir is (0|0|0) use 
	// (0|0|1) as the normal vector.
  // Otherwise use objDir as the normal.
	// Exception: if user's eye is inside object and probe mode is off, 
	// then use viewDir as the normal.
  if (clipMode)
  {
    normal.copy(&clipNormal);
  }
  else if (pm.isProjOrtho() || (viewDir.e[0]==0.0f && viewDir.e[1]==0.0f && viewDir.e[2]==0.0f))
  {
    // Draw slices parallel to projection plane:
    normal.set(0.0f, 0.0f, 1.0f);   // (0|0|1) is normal on projection plane
    normal.multiply(&invMV);
    origin.zero();
    origin.multiply(&invMV);
    normal.sub(&origin);
  }
  else if (probeSize<=0.0f && isInVolume(&eye))
  {
    // Draw slices perpendicular to viewing direction:
    normal.copy(&viewDir);   
    normal.negate();    // viewDir points away from user, the normal should point towards them
  }
  else
  {
    // Draw slices perpendicular to line eye->object:
		normal.copy(&objDir);
		normal.negate();
  }

  // Compute distance vector between textures:
  normal.normalize();
  delta.copy(&normal);
  delta.scale(diagonal / ((float)numSlices));

  // Compute farthest point to draw texture at:
	farthest.copy(&delta);
	farthest.scale((float)(numSlices - 1) / -2.0f);
  farthest.add(&probePosObj);
  
  if (clipMode)     // clipping plane present?
  {
    // Adjust numSlices and set farthest point so that textures are only
    // drawn up to the clipPoint. (Delta may not be changed
    // due to the automatic opacity correction.)
    // First find point on clipping plane which is on a line perpendicular
    // to clipping plane and which traverses the origin:
    temp.copy(&delta);
    temp.scale(-0.5f);
    farthest.add(&temp);    		// add a half delta to farthest
		clipPosObj.copy(&clipPoint);
		clipPosObj.sub(&pos);
		temp.copy(&probePosObj);
		temp.add(&normal);
    normClipPoint.isectPlaneLine(&normal, &clipPosObj, &probePosObj, &temp);
    maxDist = farthest.distance(&normClipPoint);
    numSlices = (int)(maxDist / delta.length()) + 1;
    temp.copy(&delta);
    temp.scale((float)(1 - numSlices));
    farthest.copy(&normClipPoint);
    farthest.add(&temp);
    if (singleSlice)
    {
      // Compute slice position:
      delta.scale((float)(numSlices-1));
      farthest.add(&delta);
      numSlices = 1;

      // Make slice opaque if possible:
      if (instantClassification())
        makeColorLUT(0.0f);
    }
  }
    
  // Translate object by its position:
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glTranslatef(pos.e[0], pos.e[1], pos.e[2]);

  vvVector3 texPoint;     // arbitrary point on current texture
  int isectCnt;           // intersection counter
  int j,k;              	// counters
  int drawn = 0;      		// counter for drawn textures

  // Volume render a 3D texture:
  glEnable(GL_TEXTURE_3D_EXT);
  glBindTexture(GL_TEXTURE_3D_EXT, texNames[vd->getCurrentFrame()]);

  texPoint.copy(&farthest);
  for (i=0; i<numSlices; ++i)    // loop thru all drawn textures
  {
    // Search for intersections between texture plane (defined by texPoint and 
    // normal) and texture object (0..1):
    isectCnt = isect->isectPlaneCuboid(&normal, &texPoint, &probeMin, &probeMax);

    texPoint.add(&delta);
    if (isectCnt<3) continue;   // at least 3 intersections needed for drawing
    
    // Put the intersecting 3 to 6 vertices in cyclic order to draw adjacent
    // and non-overlapping triangles:
    isect->cyclicSort(isectCnt, &normal); 

    // Generate vertices in texture coordinates:
    for (j=0; j<isectCnt; ++j)
    {
      for (k=0; k<3; ++k)
      {
        texcoord[j][k] = (isect[j][k] + size2.e[k]) / size.e[k]; 
        texcoord[j][k] = texcoord[j][k] * (texMax[k] - texMin[k]) + texMin[k];
      }
    }

    glBegin(GL_TRIANGLE_FAN);
    glColor4f(1.0, 1.0, 1.0, 1.0); 
    glNormal3f(normal[0], normal[1], normal[2]);
    ++drawn;    
    for (j=0; j<isectCnt; ++j)
    {
      glTexCoord3f(texcoord[j][0], texcoord[j][1], texcoord[j][2]);   
      glVertex3f(isect[j][0], isect[j][1], isect[j][2]); 
    }
    glEnd();
  }
  vvDebugMsg::msg(3, "Number of textures drawn: ", drawn);
  glDisable(GL_TEXTURE_3D_EXT);

  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

//----------------------------------------------------------------------------
/** Render the volume using a 3D texture (needs 3D texturing extension).
  Spherical slices are surrounding the observer.
  @param view       model-view matrix
*/
void vvTexRend::renderTex3DSpherical(vvMatrix4* view)
{
  vvVector3 temp;        // temporary vector
  int i, k;
  int ix, iy, iz;
  float  spacing;      // texture spacing
  float  maxDist = 0.0;
  float  minDist = 0.0;
  vvVector3 texSize;    // size of 3D texture [object space]
  vvVector3 texSize2;   // half size of 3D texture [object space]
  vvVector3 volumeVertices[8];
  vvMatrix4 invView;
  int optNumSlices;
  int numShells;
  
  if (!extTex3d) return;

  vvDebugMsg::msg(3, "vvTexRend::renderTex3DSpherical()");

	optNumSlices = (int)(sqrt((double)(
    vd->vox[0] * vd->vox[0] + 
    vd->vox[1] * vd->vox[1] + 
		vd->vox[2] * vd->vox[2])));  // diagonal
	numShells = (int)(quality * optNumSlices);
  if (numShells < 1)  // make sure that at least one shell is drawn
    numShells = 1;

  // Determine texture object dimensions:
  for (i=0; i<3; ++i)
  {
    texSize.e[i]  = size.e[i] * (float)texels[i] / (float)vd->vox[i];
    texSize2.e[i] = 0.5f * texSize.e[i];
  }

  invView.copy(view);
  invView.invert();

  // generates the vertices of the cube (volume) in world coordinates
  int vertexIdx = 0;
  for (ix=0; ix<2; ++ix) 
    for (iy=0; iy<2; ++iy) 
      for (iz=0; iz<2; ++iz) 
      {
	      volumeVertices[vertexIdx].e[0] = (float)ix;
	      volumeVertices[vertexIdx].e[1] = (float)iy;
	      volumeVertices[vertexIdx].e[2] = (float)iz;
	      // transfers vertices to world coordinates:
	      for (k=0; k<3; ++k)
	        volumeVertices[vertexIdx].e[k] = 
	          (volumeVertices[vertexIdx].e[k] * 2.0f - 1.0f) * texSize2.e[k];
	      volumeVertices[vertexIdx].multiply(view);
	      vertexIdx++;
      }

  // Determine maximal and minimal distance of the volume from the eyepoint:
  maxDist = minDist = volumeVertices[0].length();
  for (i = 1; i<7; i++) 
  {
    float dist = volumeVertices[i].length();
    if (dist > maxDist)  maxDist = dist;
    if (dist < minDist)  minDist = dist;
  }

  maxDist *= 1.4f;
  minDist *= 0.5f;

  // transfer the eyepoint to the object coordinates of the volume
  // to check whether the camera is inside the volume:
  vvVector3 eye(0.0,0.0,0.0);
  eye.multiply(&invView);
  int inside = 1;
  for (k=0; k<3; ++k)
  {
    if (eye.e[k] < -texSize2.e[k] || eye.e[k] > texSize2.e[k])
      inside = 0;
  }
  if (inside != 0) 
    minDist = 0.0f;

  // Determine texture spacing:
  spacing = (maxDist-minDist) / (float)(numShells-1); 

  vvSphere shell;
  shell.subdivide();
  shell.subdivide();
  shell.setVolumeDim(&texSize);
  shell.setViewMatrix(view);
  float offset[3];
  for (k=0; k<3; ++k)
    offset[k] = -(0.5f - (texMin[k] + texMax[k]) / 2.0f);
  shell.setTextureOffset(offset);

  float  radius;

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  // Volume render a 3D texture:
  glEnable(GL_TEXTURE_3D_EXT);
  glBindTexture(GL_TEXTURE_3D_EXT, texNames[0]);

  // Enable clipping plane if appropriate:
  if (clipMode) activateClippingPlane();

  radius = maxDist;
  for (i=0; i<numShells; ++i)    // loop thru all drawn textures
  {    
    shell.setRadius(radius);
    shell.calculateTexCoords();
    shell.performCulling();
    glColor4f(1.0, 1.0, 1.0, 1.0);
    shell.render();
    radius -= spacing;
  }
  deactivateClippingPlane();
  glDisable(GL_TEXTURE_3D_EXT);
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

//----------------------------------------------------------------------------
/** Render the volume using a 3D texture (needs 3D texturing extension).
  The slices are drawn parallel to the volume faces.
  @param principal  principal viewing axis
  @param zx,zy,zz   z coordinates of transformed base vectors
*/
void vvTexRend::renderTex3DFaces(AxisType principal, float zx, float zy, float zz)
{
  vvVector3 normal;                            // normal vector for slices
  vvVector3 texTL, texTR, texBL, texBR;        // texture coordinates (T=top etc.)
  vvVector3 texSpacing;                        // spacing for texture coordinates
  vvVector3 objTL, objTR, objBL, objBR;        // object coordinates in world space
  vvVector3 objSpacing;                        // texture spacing in world space
  vvVector3 pos;                               // volume location
  float     size2[3];                          // half object size
  int       numTextures;                       // number of texture slices drawn
	int       optNumTextures;    						  	 // optimum number of texture slices drawn
  int i;

  vvDebugMsg::msg(3, "vvTexRend::renderTex3DFaces()");

  if (!extTex3d) return;

  pos.copy(&vd->position);

  // Initialize texture counter:
  optNumTextures = ts_max(texels[0], texels[1], texels[2]);
  numTextures = (int)(quality * optNumTextures);  
  if (numTextures < 1)  // make sure that at least one slice is drawn
    numTextures = 1;

  vvDebugMsg::msg(3, "Number of textures needed:   ", optNumTextures);
  vvDebugMsg::msg(3, "Number of textures drawn:    ", numTextures);

  // Generate half object size as a shortcut:
  for (i=0; i<3; ++i)
    size2[i] = 0.5f * size.e[i];

  // Initialize parameters upon principal viewing direction:
  switch (principal)
  {
    case X_AXIS: // zx>0 -> draw left to right
      // Coordinate system:
      //     z
      //     |__y
      //   x/
      objTL.set(-size2[0],-size2[1], size2[2]);
      objTR.set(-size2[0], size2[1], size2[2]);
      objBL.set(-size2[0],-size2[1],-size2[2]);
      objBR.set(-size2[0], size2[1],-size2[2]);

      texTL.set(texMin[0], texMin[1], texMax[2]);
      texTR.set(texMin[0], texMax[1], texMax[2]);
      texBL.set(texMin[0], texMin[1], texMin[2]);
      texBR.set(texMin[0], texMax[1], texMin[2]);

      objSpacing.set(2.0f * size2[0] / (float)(numTextures - 1), 0.0f, 0.0f);
      texSpacing.set((texMax[0] - texMin[0]) / (float)(numTextures - 1), 0.0f, 0.0f);
      normal.set(1.0f, 0.0f, 0.0f);
      if (zx<0)     // reverse order?
      {
        normal.e[0]     = -normal.e[0];
        objTL.e[0]      = objTR.e[0] = objBL.e[0] = objBR.e[0] = size2[0];
        texTL.e[0]      = texTR.e[0] = texBL.e[0] = texBR.e[0] = texMax[0];
        objSpacing.e[0] = -objSpacing.e[0];
        texSpacing.e[0] = -texSpacing.e[0];
      }
      break; 
      
    case Y_AXIS:  // zy>0 -> draw bottom to top
      // Coordinate system:
      //     x
      //     |__z
      //   y/
      objTL.set( size2[0],-size2[1],-size2[2]);
      objTR.set( size2[0],-size2[1], size2[2]);
      objBL.set(-size2[0],-size2[1],-size2[2]);
      objBR.set(-size2[0],-size2[1], size2[2]);

      texTL.set(texMax[0], texMin[1], texMin[2]);
      texTR.set(texMax[0], texMin[1], texMax[2]);
      texBL.set(texMin[0], texMin[1], texMin[2]);
      texBR.set(texMin[0], texMin[1], texMax[2]);

      objSpacing.set(0.0f, 2.0f * size2[1] / (float)(numTextures - 1), 0.0f);
      texSpacing.set(0.0f, (texMax[1] - texMin[1]) / (float)(numTextures - 1), 0.0f);
      normal.set(0.0f, 1.0f, 0.0f);
      if (zy<0)     // reverse order?
      {
        normal.e[1]     = -normal.e[1];
        objTL.e[1]      = objTR.e[1] = objBL.e[1] = objBR.e[1] = size2[1];
        texTL.e[1]      = texTR.e[1] = texBL.e[1] = texBR.e[1] = texMax[1];
        objSpacing.e[1] = -objSpacing.e[1];
        texSpacing.e[1] = -texSpacing.e[1];
      }
      break; 
      
    case Z_AXIS: // zz>0 -> draw back to front
      // Coordinate system:
      //     y
      //     |__x
      //   z/
      objTL.set(-size2[0], size2[1],-size2[2]);
      objTR.set( size2[0], size2[1],-size2[2]);
      objBL.set(-size2[0],-size2[1],-size2[2]);
      objBR.set( size2[0],-size2[1],-size2[2]);

      texTL.set(texMin[0], texMax[1], texMin[2]);
      texTR.set(texMax[0], texMax[1], texMin[2]);
      texBL.set(texMin[0], texMin[1], texMin[2]);
      texBR.set(texMax[0], texMin[1], texMin[2]);

      objSpacing.set(0.0f, 0.0f, size.e[2] / (float)(numTextures - 1));
      texSpacing.set(0.0f, 0.0f, (texMax[2]-texMin[2]) / (float)(numTextures - 1));
      normal.set(0.0f, 0.0f, 1.0f);
      if (zz<0)     // reverse order?
      {
        normal.e[2]     = -normal.e[2];
        objTL.e[2]      = objTR.e[2] = objBL.e[2] = objBR.e[2] = size2[2];
        texTL.e[2]      = texTR.e[2] = texBL.e[2] = texBR.e[2] = texMax[2];
        objSpacing.e[2] = -objSpacing.e[2];
        texSpacing.e[2] = -texSpacing.e[2];
      }
      break; 
      
    default: break;
  }

  // Translate object by its position:
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glTranslatef(pos.e[0], pos.e[1], pos.e[2]);

  // Enable clipping plane if appropriate:
  if (clipMode) activateClippingPlane();

  // Volume render a 3D texture:
  glEnable(GL_TEXTURE_3D_EXT);
  glBindTexture(GL_TEXTURE_3D_EXT, texNames[vd->getCurrentFrame()]);
  for (i=0; i<numTextures; ++i)
  {
    glBegin(GL_QUADS);
      glColor4f(1.0, 1.0, 1.0, 1.0); 
      glNormal3f(normal.e[0], normal.e[1], normal.e[2]);
      glTexCoord3f(texTL.e[0], texTL.e[1], texTL.e[2]); 
        glVertex3f(objTL.e[0], objTL.e[1], objTL.e[2]);   
      glTexCoord3f(texBL.e[0], texBL.e[1], texBL.e[2]); 
        glVertex3f(objBL.e[0], objBL.e[1], objBL.e[2]);   
      glTexCoord3f(texBR.e[0], texBR.e[1], texBR.e[2]); 
        glVertex3f(objBR.e[0], objBR.e[1], objBR.e[2]);   
      glTexCoord3f(texTR.e[0], texTR.e[1], texTR.e[2]); 
        glVertex3f(objTR.e[0], objTR.e[1], objTR.e[2]);   
    glEnd();
    texTL.add(&texSpacing);
    texBL.add(&texSpacing);
    texBR.add(&texSpacing);
    texTR.add(&texSpacing);
    objTL.add(&objSpacing);
    objBL.add(&objSpacing);
    objBR.add(&objSpacing);
    objTR.add(&objSpacing);
  }
  deactivateClippingPlane();
  glDisable(GL_TEXTURE_3D_EXT);
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

//----------------------------------------------------------------------------
/** Render the volume using 2D textures (OpenGL 1.1 compatible).
  @param zz  z coordinate of transformed z base vector
*/
void vvTexRend::renderTex2DSlices(float zz)
{
  vvVector3 normal;          // normal vector for slices
  vvVector3 size2;           // half texture size
  vvVector3 pos;             // object location
  float     texSpacing;      // spacing for texture coordinates
  float     zPos;            // texture z position
  float     texStep;         // step between texture indices
  float     texIndex;        // current texture index
  int       numTextures;     // number of textures drawn
  int       i;

  vvDebugMsg::msg(3, "vvTexRend::renderTex2DSlices()");

  // Translate object by its position:
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  pos.copy(&vd->position);
  glTranslatef(pos.e[0], pos.e[1], pos.e[2]);

  // Enable clipping plane if appropriate:
  if (clipMode) activateClippingPlane();

  // Generate half object size as shortcut:
  size2.e[0] = 0.5f * size.e[0];
  size2.e[1] = 0.5f * size.e[1];
  size2.e[2] = 0.5f * size.e[2];

  numTextures = (int)(quality * vd->vox[2]);
  if (numTextures < 1) numTextures = 1;

  normal.set(0.0f, 0.0f, 1.0f);
  zPos = -size2.e[2];
  if (numTextures>1)  // prevent division by zero
  {
    texSpacing = size.e[2] / (float)(numTextures - 1);
    texStep = (float)(vd->vox[2] - 1) / (float)(numTextures - 1);
  }
  else
  {
    texSpacing = 0.0f;
    texStep = 0.0f;
    zPos = 0.0f;
  }
  
  texIndex = float(vd->getCurrentFrame()) * float(vd->vox[2]);    // offset for current time step
  if (zz>0.0f) // draw textures back to front?
  {
    texIndex += float(vd->vox[2] - 1);
    texStep  = -texStep;
  }
  else    // draw textures front to back
  {
    zPos        = -zPos;
    texSpacing  = -texSpacing;
    normal.e[2] = -normal.e[2];
  }
  
  // Volume rendering with multiple 2D textures:
  glEnable(GL_TEXTURE_2D);

  for (i=0; i<numTextures; ++i)
  { 
    glBindTexture(GL_TEXTURE_2D, texNames[vvToolshed::round(texIndex)]); 

    glBegin(GL_QUADS);
      glColor4f(1.0, 1.0, 1.0, 1.0); 
      glNormal3f(normal.e[0], normal.e[1], normal.e[2]);
      glTexCoord2f(texMin[0], texMin[1]); glVertex3f(-size2.e[0],  size2.e[1], zPos);
      glTexCoord2f(texMin[0], texMax[1]); glVertex3f(-size2.e[0], -size2.e[1], zPos);
      glTexCoord2f(texMax[0], texMax[1]); glVertex3f( size2.e[0], -size2.e[1], zPos);
      glTexCoord2f(texMax[0], texMin[1]); glVertex3f( size2.e[0],  size2.e[1], zPos); 
    glEnd();

    zPos += texSpacing;
    texIndex += texStep;
  }

  glDisable(GL_TEXTURE_2D);
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  deactivateClippingPlane();
  vvDebugMsg::msg(3, "Number of textures stored: ", vd->vox[2]);
  vvDebugMsg::msg(3, "Number of textures drawn:  ", numTextures);
}

//----------------------------------------------------------------------------
/** Render the volume using 2D textures, switching to the optimum
    texture set to prevent holes. 
  @param principal  principal viewing axis
  @param zx,zy,zz   z coordinates of transformed base vectors
*/
void vvTexRend::renderTex2DCubic(AxisType principal, float zx, float zy, float zz)
{
  vvVector3 normal;                         // normal vector for slices
  vvVector3 texTL, texTR, texBL, texBR;     // texture coordinates (T=top etc.)
  vvVector3 objTL, objTR, objBL, objBR;     // object coordinates in world space
  vvVector3 texSpacing;                     // distance between textures
  vvVector3 pos;                            // object location
  float  size2[3];                          // half object sizes
  float  texStep;                           // step size for texture names
  float  texIndex;                          // textures index
  int    numTextures;                       // number of textures drawn
  int    frameTextures;                     // number of textures per frame
  int    i;

  vvDebugMsg::msg(3, "vvTexRend::renderTex2DCubic()");

  // Translate object by its position:
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  pos.copy(&vd->position);
  glTranslatef(pos.e[0], pos.e[1], pos.e[2]);

  // Enable clipping plane if appropriate:
  if (clipMode) activateClippingPlane();

  // Initialize texture parameters:
  frameTextures = vd->vox[0] + vd->vox[1] + vd->vox[2];
  numTextures = (int)(quality * (float)ts_max(vd->vox[0], vd->vox[1], vd->vox[2]));  
  if (numTextures < 2)  numTextures = 2;   // make sure that at least one slice is drawn to prevent division by zero

  // Generate half object size as a shortcut:
  for (i=0; i<3; ++i)
    size2[i] = 0.5f * size.e[i];

  // Initialize parameters upon principal viewing direction:
  switch (principal)
  {
    case X_AXIS: // zx>0 -> draw left to right
      // Coordinate system:
      //     z
      //     |__y
      //   x/
      objTL.set(-size2[0],-size2[1], size2[2]);
      objTR.set(-size2[0], size2[1], size2[2]);
      objBL.set(-size2[0],-size2[1],-size2[2]);
      objBR.set(-size2[0], size2[1],-size2[2]);

      texTL.set(texMin[1], texMax[2], 0.0f);
      texTR.set(texMax[1], texMax[2], 0.0f);
      texBL.set(texMin[1], texMin[2], 0.0f);
      texBR.set(texMax[1], texMin[2], 0.0f);

      texSpacing.set(size.e[0] / float(numTextures - 1), 0.0f, 0.0f);
      texStep = -1.0f * float(vd->vox[0] - 1) / float(numTextures - 1);
      normal.set(1.0f, 0.0f, 0.0f);
      texIndex = float(vd->getCurrentFrame() * frameTextures);
      if (zx<0)     // reverse order? draw right to left
      {
        normal.e[0]     = -normal.e[0];
        objTL.e[0]      = objTR.e[0] = objBL.e[0] = objBR.e[0] = size2[0];
        texSpacing.e[0] = -texSpacing.e[0];
        texStep         = -texStep;
      }
      else
      {
        texIndex += float(vd->vox[0] - 1);
      }
      break; 
      
    case Y_AXIS:  // zy>0 -> draw bottom to top
      // Coordinate system:
      //     x
      //     |__z
      //   y/
      objTL.set( size2[0],-size2[1],-size2[2]);
      objTR.set( size2[0],-size2[1], size2[2]);
      objBL.set(-size2[0],-size2[1],-size2[2]);
      objBR.set(-size2[0],-size2[1], size2[2]);

      texTL.set(texMin[2], texMax[0], 0.0f);
      texTR.set(texMax[2], texMax[0], 0.0f);
      texBL.set(texMin[2], texMin[0], 0.0f);
      texBR.set(texMax[2], texMin[0], 0.0f);
      
      texSpacing.set(0.0f, size.e[1] / float(numTextures - 1), 0.0f);
      texStep = -1.0f * float(vd->vox[1] - 1) / float(numTextures - 1);
      normal.set(0.0f, 1.0f, 0.0f);
      texIndex = float(vd->getCurrentFrame() * frameTextures + vd->vox[0]);
      if (zy<0)     // reverse order? draw top to bottom
      {
        normal.e[1]     = -normal.e[1];
        objTL.e[1]      = objTR.e[1] = objBL.e[1] = objBR.e[1] = size2[1];
        texSpacing.e[1] = -texSpacing.e[1];
        texStep         = -texStep;
      }
      else
      {
        texIndex += float(vd->vox[1] - 1);
      }
      break; 

    case Z_AXIS: // zz>0 -> draw back to front
    default:
      // Coordinate system:
      //     y
      //     |__x
      //   z/
      objTL.set(-size2[0], size2[1],-size2[2]);
      objTR.set( size2[0], size2[1],-size2[2]);
      objBL.set(-size2[0],-size2[1],-size2[2]);
      objBR.set( size2[0],-size2[1],-size2[2]);

      texTL.set(texMin[0], texMax[1], 0.0f);
      texTR.set(texMax[0], texMax[1], 0.0f);
      texBL.set(texMin[0], texMin[1], 0.0f);
      texBR.set(texMax[0], texMin[1], 0.0f);
      
      texSpacing.set(0.0f, 0.0f, size.e[2] / float(numTextures - 1));
      normal.set(0.0f, 0.0f, 1.0f);
      texStep = -1.0f * float(vd->vox[2] - 1) / float(numTextures - 1);
      texIndex = float(vd->getCurrentFrame() * frameTextures + vd->vox[0] + vd->vox[1]);
      if (zz<0)     // reverse order? draw front to back
      {
        normal.e[2]     = -normal.e[2];
        objTL.e[2]      = objTR.e[2] = objBL.e[2] = objBR.e[2] = size2[2];
        texSpacing.e[2] = -texSpacing.e[2];
        texStep         = -texStep;
      }
      else   // draw back to front
      {
        texIndex += float(vd->vox[2] - 1);
      }
      break; 
  }

  // Volume render a 2D texture:
  glEnable(GL_TEXTURE_2D);
  for (i=0; i<numTextures; ++i)
  {
    glBindTexture(GL_TEXTURE_2D, texNames[vvToolshed::round(texIndex)]);

    glBegin(GL_QUADS);
      glColor4f(1.0, 1.0, 1.0, 1.0); 
      glNormal3f(normal.e[0], normal.e[1], normal.e[2]);
      glTexCoord2f(texTL.e[0], texTL.e[1]); glVertex3f(objTL.e[0], objTL.e[1], objTL.e[2]);   
      glTexCoord2f(texBL.e[0], texBL.e[1]); glVertex3f(objBL.e[0], objBL.e[1], objBL.e[2]);   
      glTexCoord2f(texBR.e[0], texBR.e[1]); glVertex3f(objBR.e[0], objBR.e[1], objBR.e[2]);   
      glTexCoord2f(texTR.e[0], texTR.e[1]); glVertex3f(objTR.e[0], objTR.e[1], objTR.e[2]);   
    glEnd();
    objTL.add(&texSpacing);
    objBL.add(&texSpacing);
    objBR.add(&texSpacing);
    objTR.add(&texSpacing);
    
    texIndex += texStep;
  }
  glDisable(GL_TEXTURE_2D);
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  deactivateClippingPlane();
}

//----------------------------------------------------------------------------
/** Render the volume onto currently selected drawBuffer.
 Viewport size in world coordinates is -1.0 .. +1.0 in both x and y direction
*/
void vvTexRend::renderVolumeGL()
{
  AxisType principal;                       // principal viewing direction
  vvMatrix4 mv;                             // current modelview matrix
  vvVector3 origin(0.0f, 0.0f, 0.0f);       // zero vector
  vvVector3 xAxis(1.0f, 0.0f, 0.0f);        // vector in x axis direction
  vvVector3 yAxis(0.0f, 1.0f, 0.0f);        // vector in y axis direction
  vvVector3 zAxis(0.0f, 0.0f, 1.0f);        // vector in z axis direction
	vvVector3 probeSizeObj; 									// probe size [object space]
  vvStopwatch* sw = NULL;                   // stop watch for performance measurements
  float zx, zy, zz;                         // base vector z coordinates
	float maxEdge;													  // maximum edge length
  int i;

  vvDebugMsg::msg(3, "vvTexRend::renderVolumeGL()");

  if (timing)
  {
    sw = new vvStopwatch();
    sw->start();
  }

  setGLenvironment();

  // Determine texture object extensions:
  for (i=0; i<3; ++i)
  {
    texMin[i] = 0.5f / (float)texels[i];
    texMax[i] = (float)vd->vox[i] / (float)texels[i] - texMin[i];
  }

  // Find principle viewing direction:

  // Get OpenGL modelview matrix:
  getModelviewMatrix(&mv);

  // Transform 4 point vectors with the modelview matrix:
  origin.multiply(&mv);
  xAxis.multiply(&mv);
  yAxis.multiply(&mv);
  zAxis.multiply(&mv);

  // Generate coordinate system base vectors from those vectors:
  xAxis.sub(&origin);
  yAxis.sub(&origin);
  zAxis.sub(&origin);

  xAxis.normalize();
  yAxis.normalize();
  zAxis.normalize();

  // Only z component of base vectors is needed:
  zx = xAxis.e[2];
  zy = yAxis.e[2];
  zz = zAxis.e[2];

  if (fabs(zx) > fabs(zy))
  {
    if (fabs(zx) > fabs(zz)) principal = X_AXIS;
    else principal = Z_AXIS;
  }
  else 
  {
    if (fabs(zy) > fabs(zz)) principal = Y_AXIS;
    else principal = Z_AXIS;
  }

  // Draw all boundaries:
  if (boundaries || probeSize>0.0f)
  {
    drawBoundingBox(&size, &vd->position, boundColor, false);
		if (probeSize>0.0f)
		{
			maxEdge = ts_max(size[0], size[1], size[2]);
			probeSizeObj.set(maxEdge * probeSize, maxEdge * probeSize, maxEdge * probeSize); 
	    drawBoundingBox(&probeSizeObj, &probePos, probeColor, false);
		}
  }

  switch (methodUsed)
  {
    default:
    case SLICES2D:  renderTex2DSlices(zz); break;
    case CUBIC2D:   renderTex2DCubic(principal, zx, zy, zz); break;
    case CUBIC3D:   renderTex3DFaces(principal, zx, zy, zz); break;
    case SPHERIC3D: renderTex3DSpherical(&mv); break;
    case INDEX3D:   
      if (extColLUT)
        glEnable(GL_TEXTURE_COLOR_TABLE_SGI);
    // fall through
    case PLANAR3D:  
			renderTex3DPlanar(&mv);
      break;
  }

  // Draw visible boundaries:
  if (boundaries || probeSize>0.0f)
  {
    drawBoundingBox(&size, &vd->position, boundColor, true);
		if (probeSize>0.0f)
		{
			maxEdge = ts_max(size[0], size[1], size[2]);
			probeSizeObj.set(maxEdge * probeSize, maxEdge * probeSize, maxEdge * probeSize); 
	    drawBoundingBox(&probeSizeObj, &probePos, probeColor, true);
		}
  }

  vvRenderer::renderVolumeGL();

  if (timing)
  {
    cerr << "Rendering time [ms]: " << setprecision(7) << ((float)sw->getTime() * 1000.0f) << endl;
    delete sw;
  }

  unsetGLenvironment();
}

//----------------------------------------------------------------------------
/** Activate the previously set clipping plane.
    Clipping plane parameters have to be set with setClippingPlane().
*/
void vvTexRend::activateClippingPlane()
{
  GLdouble planeEq[4];    // plane equation
  vvVector3 clipNormal2;  // clipping normal pointing to opposite direction
  float thickness;        // thickness of single slice clipping plane

  vvDebugMsg::msg(3, "vvTexRend::activateClippingPlane()");

  // Generate OpenGL compatible clipping plane parameters:
  // normal points into oppisite direction
  planeEq[0] = -clipNormal[0];
  planeEq[1] = -clipNormal[1];
  planeEq[2] = -clipNormal[2];
  planeEq[3] = clipNormal.dot(&clipPoint);
  glClipPlane(GL_CLIP_PLANE0, planeEq);
  glEnable(GL_CLIP_PLANE0);

  // Generate second clipping plane in single slice mode:
  if (singleSlice)
  {
    thickness = size[0] / 100.0f;
    clipNormal2.copy(&clipNormal);
    clipNormal2.negate();
    planeEq[0] = -clipNormal2[0];
    planeEq[1] = -clipNormal2[1];
    planeEq[2] = -clipNormal2[2];
    planeEq[3] = clipNormal2.dot(&clipPoint) + thickness;
    glClipPlane(GL_CLIP_PLANE1, planeEq);
    glEnable(GL_CLIP_PLANE1);
  }
}

//----------------------------------------------------------------------------
/** Deactivate the clipping plane.
*/
void vvTexRend::deactivateClippingPlane()
{
  vvDebugMsg::msg(3, "vvTexRend::deactivateClippingPlane()");
  glDisable(GL_CLIP_PLANE0);
  if (singleSlice) glDisable(GL_CLIP_PLANE1);
}

//----------------------------------------------------------------------------
/** Set number of lights in the scene. 
  Fixed material characteristics are used with each setting.
  @param numLights  number of lights in scene (0=ambient light only)
*/
void vvTexRend::setNumLights(int numLights)
{
  float ambient[]  = {0.5f, 0.5f, 0.5f, 1.0f};
  float pos0[] = {0.0f, 10.0f, 10.0f, 0.0f};
  float pos1[] = {0.0f, -10.0f, -10.0f, 0.0f};

  vvDebugMsg::msg(1, "vvTexRend::setNumLights()");

  // Generate light source 1:
  glEnable(GL_LIGHT0);
  glLightfv(GL_LIGHT0, GL_POSITION, pos0);
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);

  // Generate light source 2:
  glEnable(GL_LIGHT1);
  glLightfv(GL_LIGHT1, GL_POSITION, pos1);
  glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);

  // At least 2 lights:
  if (numLights >= 2)
    glEnable(GL_LIGHT1);
  else 
    glDisable(GL_LIGHT1);

  // At least one light:
  if (numLights >= 1) 
    glEnable(GL_LIGHT0);
  else // no lights selected
    glDisable(GL_LIGHT0);
}

//----------------------------------------------------------------------------
/// Set voxel size and resize object accordingly.
void vvTexRend::setVoxelSize(vvVector3* s)
{
  vvDebugMsg::msg(3, "vvTexRend::setVoxelSize()");
  for (int i=0; i<3; ++i)
    size.e[i] *= s->e[i] / vd->dist[i];
  vvRenderer::setVoxelSize(s);
}

//----------------------------------------------------------------------------
/// @return true if classification is done in no time
bool vvTexRend::instantClassification()
{
  vvDebugMsg::msg(3, "vvTexRend::instantClassification()");
  if (methodUsed==INDEX3D) return true;
  else return false;
}

//----------------------------------------------------------------------------
/// Returns the number of entries in the RGBA lookup table.
int vvTexRend::getLUTSize()
{
  vvDebugMsg::msg(1, "vvTexRend::getLUTSize()");
  return (vd->bpv==2 && methodUsed!=INDEX3D) ? 4096 : 256;
}

//----------------------------------------------------------------------------
/** Make a color lookup table.
 Important observation: glColorTableSGI can have a maximum width of
 1024 RGBA entries on IR2 graphics!
 @param dist  slice distance relative to 3D texture sample point distance 
              (1.0 for 1:1 mapping, 0.0 for opaque representation).
*/
void vvTexRend::makeColorLUT(float dist)
{
  vvDebugMsg::msg(1, "Generating texture color lookup table. Slice distance = ", dist);

  uchar tmpConv[4096][4];     // temporary RGBA conversion table
  float aOrig;                // original alpha channel value [0..1]
  float aCorr;                // corrected alpha channel value [0..1]
  int lutEntries;             // number of entries in the RGBA lookup table
  int i;

  lutEntries = getLUTSize();
  if (extColLUT || extPalTex)
  {
    if (dist==1)    // is slice distance equal to sampling distance?
    {
      if (extColLUT)
      {
        glColorTableSGI(GL_TEXTURE_COLOR_TABLE_SGI, GL_RGBA, 
          lutEntries, GL_RGBA, GL_UNSIGNED_BYTE, &rgbaConv[0][0]);
      }
      else
      {
        glColorTableEXT(GL_TEXTURE_3D_EXT, GL_RGBA, 
          lutEntries, GL_RGBA, GL_UNSIGNED_BYTE, &rgbaConv[0][0]);
      }
    }
    else
    {
      // Copy LUT entries and modify alpha channel:
      memcpy(&tmpConv[0][0], &rgbaConv[0][0], lutEntries * 4 * sizeof(uchar));
      for (i=0; i<lutEntries; ++i)
      {
        if (dist<=0.0)    // for 0 distance draw opaque slices
          aCorr = 1.0f;
        else
        {
          aOrig = (float)rgbaConv[i][3] / 255.0f;
          aCorr = 1.0f - powf(1.0f - aOrig, dist);
        }
        tmpConv[i][3] = (uchar)(aCorr * 255.0f);
      }
      
      if (extColLUT)
      {
        glColorTableSGI(GL_TEXTURE_COLOR_TABLE_SGI, GL_RGBA, 
          lutEntries, GL_RGBA, GL_UNSIGNED_BYTE, &tmpConv[0][0]);
      }
      else
      {
        glColorTableEXT(GL_TEXTURE_3D_EXT, GL_RGBA, 
          lutEntries, GL_RGBA, GL_UNSIGNED_BYTE, &tmpConv[0][0]);
      }
    }
  }
}

//----------------------------------------------------------------------------
/** Set user's viewing direction.
  This information is needed to correctly orientate the texture slices
  in 3D texturing mode if the user is inside the volume.
  @param vd  viewing direction in object coordinates
*/
void vvTexRend::setViewingDirection(vvVector3* vd)
{
  vvDebugMsg::msg(3, "vvTexRend::setViewingDirection()");
  viewDir.copy(vd);
}

//----------------------------------------------------------------------------
/** Set the direction from the viewer to the object.
  This information is needed to correctly orientate the texture slices
  in 3D texturing mode if the viewer is outside of the volume.
  @param vd  object direction in object coordinates
*/
void vvTexRend::setObjectDirection(vvVector3* vd)
{
  vvDebugMsg::msg(3, "vvTexRend::setObjectDirection()");
  objDir.copy(vd);
}

//----------------------------------------------------------------------------
// see parent
void vvTexRend::setParameter(float newValue, ParameterType param, char*)
{
  vvDebugMsg::msg(3, "vvTexRend::setParameter()");
  switch (param)
  {
    case vvRenderer::VV_SLICEINT:
      interpolation = (newValue == 0.0f) ? false : true;
      makeTextures();
      break;
    default: break;
  }
}

//----------------------------------------------------------------------------
// see parent
float vvTexRend::getParameter(ParameterType param, char*)
{
  vvDebugMsg::msg(3, "vvTexRend::getParameter()");

  switch (param)
  {
    case vvRenderer::VV_SLICEINT:
      return (interpolation) ? 1.0f : 0.0f;
      break;
    default: return vvRenderer::getParameter(param);
  }
}

//----------------------------------------------------------------------------
/** Get information on hardware support for rendering modes.
  This routine cannot guarantee that a requested method works on any
  dataset, but it will at least work for datasets which fit into
  texture memory.
  @param rt rendering type to get information about
  @return true if the requested rendering type is supported by 
    the system's OpenGL hardware.
*/
bool vvTexRend::isSupported(RenderType rt)
{
  vvDebugMsg::msg(3, "vvTexRend::isSupported()");
  switch (rt)
  {
    case BEST:
    // fall thru
    case SLICES2D:
    // fall thru
    case CUBIC2D: return true; 
      break; 

    case CUBIC3D:   
    // fall thru
    case PLANAR3D:  
    // fall thru
    case SPHERIC3D: 
        return vvGLTools::isGLextensionSupported("GL_EXT_texture3D");

    case INDEX3D:
        return (vvGLTools::isGLextensionSupported("GL_EXT_texture3D") && 
               (vvGLTools::isGLextensionSupported("GL_SGI_texture_color_table") ||
                vvGLTools::isGLextensionSupported("GL_EXT_paletted_texture")));
    
    default: return false;
  }
}

//----------------------------------------------------------------------------
/** Return the currently used rendering method.
  The returned method does not need to be the same as the method requested
  in the constructor.
*/
vvTexRend::RenderType vvTexRend::getMethodUsed()
{
  vvDebugMsg::msg(3, "vvTexRend::getMethodUsed()");
  return methodUsed;
}

//============================================================================
// End of File
//============================================================================
