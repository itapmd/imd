//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifndef _VVRENDERER_H_
#define _VVRENDERER_H_

#include "vvvecmath.h"
#include "vvvoldesc.h"
#include "vvrenderer.h"

//============================================================================
// Class Definition
//============================================================================

/** Abstract volume rendering class.
  @author Juergen Schulze-Doebold (schulze@hlrs.de)
  @see vvSoftVR
  @see vvTexRend
  @see vvVolTex
  @see vvRenderVP
*/
class vvRenderer
{
  public:
    enum RendererType   /// Current renderer
    {
      TEXREND = 0,      ///< texture based renderer
      SOFTPAR,          ///< software based renderer for parallel projection
      SOFTPER,          ///< software based renderer for perspective projection
      VOLPACK,          ///< Phil Lacroute's VolPack renderer
      SIMIAN,           ///< Joe Kniss's Simian renderer
      IMGREND,          ///< 2D image renderer
      BLOWTEX,          ///< blown up textures renderer
      UNKNOWN           ///< unknown renderer
    };
    enum ParameterType  ///  Names for rendering parameters
    {   // important notice: if number assignments are changed, make sure to change them in Java code, too!
      VV_ALPHA     = 0, ///< alpha value for Simian widgets
      VV_EYE       = 1, ///< left or right eye for stereo
      VV_SLICEINT  = 2, ///< inter-slice interpolation on/off
      VV_WARPINT   = 3, ///< warp interpolation on/off
      VV_COMPRESS  = 4, ///< RLE data compression on/off
      VV_MULTIPROC = 5, ///< multi-processiong on/off
      VV_PREINT    = 6, ///< pre-integration on/off
      VV_SLICEBUF  = 7, ///< slice buffer mode
      VV_LOOKUP    = 8, ///< bilinear interpolation at lookup in pre-integration table on/off
      VV_OPCORR    = 9  ///< opacity correction on/off
    };

  private:
    vvMatrix4   modelview;      ///< copy of current modelview matrix
    vvMatrix4   projection;     ///< copy of current projection matrix

  protected:
    RendererType rendererType;  ///< currently used renderer type
    vvVolDesc*  vd;             ///< volume description
    vvVector3   size;           ///< volume size (width, height, depth) in world coordinates
    bool        orientation;    ///< true = display object orientation
    bool        palette;        ///< true = display transfer function palette
    bool        boundaries;     ///< true = display volume boundaries
		float       boundColor[3];  ///< boundary color (R,G,B in [0..1])
		float       probeColor[3];  ///< probe boundary color (R,G,B in [0..1])
    vvVector3   clipPoint;      ///< point on clipping plane
    vvVector3   clipNormal;     ///< clipping plane normal
    bool        clipMode;       ///< true = clipping plane enabled, false=disabled
    bool        singleSlice;    ///< true = use single slice in clipping mode
    float       probeSize;      ///< edge length of probe (fraction of entire object) [0..1], 0.0=probe mode off
    vvVector3   probePos;       ///< world space coordinates of probe midpoint [mm]
    float       quality;        ///< rendering image quality (0=minimum, 1=sampling rate, >1=oversampling)
    bool        timing;         ///< true = measure and display rendering times

    void        init();         ///< initialization routine
    
  // Class Methods:
  public:                 // public methods will be inherited as public
    vvRenderer(vvVolDesc*);
    virtual ~vvRenderer();

    // Static methods:
    static float adaptQuality(float, float, float, float);

    // Public methods that should be redefined by subclasses:
    virtual RendererType getRendererType();
    virtual void  makeHistogram(int, float*);
    virtual void  renderVolumeGL();
    virtual void  renderVolumeRGB(int, int, uchar*);
    virtual void  updateTransferFunction();
    virtual void  updateVolumeData();
    virtual void  setClippingPlane(vvVector3*, vvVector3*);
    virtual void  setClippingMode(bool);
    virtual bool  getClippingMode();
    virtual void  setSingleSliceMode(bool);
    virtual bool  getSingleSliceMode();
    virtual int   getNumFrames();
    virtual int   getCurrentFrame();
    virtual void  setCurrentFrame(int);
    virtual void  setBoundariesMode(bool);
    virtual void  setTimingMode(bool);
    virtual bool  getTimingMode();
    virtual void  setOrientationMode(bool);
		virtual void  setBoundariesColor(float, float, float);
    virtual void  getBoundariesColor(float*, float*, float*);
		virtual void  setProbeColor(float, float, float);
    virtual void  getProbeColor(float*, float*, float*);
    virtual void  setPaletteMode(bool);
    virtual void  setLights(int, bool);
		virtual void  setQuality(float);
    virtual float getQuality();
    virtual void  setPosition(vvVector3*);
    virtual void  getPosition(vvVector3*);
    virtual void  setObjectSize(vvVector3*);
    virtual void  getObjectSize(vvVector3*);
    virtual void  setVoxelSize(vvVector3*);
    virtual void  getVoxelSize(vvVector3*);
    virtual void  resizeEdgeMax(float);
    virtual void  renderCoordinates();
    virtual void  renderPalette();
    virtual void  drawBoundingBox(vvVector3*, vvVector3*, float*, bool);
    virtual bool  instantClassification();
    virtual void  setViewingDirection(vvVector3*);
    virtual void  setObjectDirection(vvVector3*);
		virtual void  setProbePosition(vvVector3*);
		virtual void  getProbePosition(vvVector3*);
		virtual void  setProbeSize(float);
		virtual float getProbeSize();
    virtual void  getModelviewMatrix(vvMatrix4*);
    virtual void  getProjectionMatrix(vvMatrix4*);
    virtual void  setModelviewMatrix(vvMatrix4*);
    virtual void  setProjectionMatrix(vvMatrix4*);
    virtual void  getEyePosition(vvVector3*);
    virtual bool  isInVolume(vvVector3*);
	  virtual float getAlphaValue(float, float, float);
    virtual void  setParameter(float, ParameterType, char* = NULL);
    virtual float getParameter(ParameterType, char* = NULL);
    virtual void  profileStart();
    virtual void  profileStop();
};

#endif

//============================================================================
// End of File
//============================================================================
