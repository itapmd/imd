//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifndef _VVIMGREND_H_
#define _VVIMGREND_H_

#ifdef WIN32
  #include <winsock2.h>
  #include <windows.h>
  #include <GL/gl.h>
  #ifndef VV_REMOTE_RENDERING
    #include <vvglext.h>
  #endif
#elif __linux
  #ifndef VV_REMOTE_RENDERING
    #include <GL/gl.h>
    #include "vvglext.h"
  #endif
#else
  #ifndef VV_REMOTE_RENDERING
    #include <GL/gl.h>
  #endif
#endif
#include "vvrenderer.h"
#include "vvimage.h"
#include "vvsocketio.h"
#include "vvvecmath.h"
#include "vvsoftimg.h"
#include "vvvoldesc.h"
#include "vvstopwatch.h"

//============================================================================
// Class Definitions
//============================================================================

/** Rendering engine for rendering 2D images which are projected on a 
  2D texture which is modified by a warp matrix. The 2D images are generated
  by a remote rendering engine: VRemote. This mechanism allows the rendering
  to take place on a remote parallel computer, while the resulting images
  can be displayed locally.
  @author Juergen Schulze (schulze@hlrs.de)
  @see vvRenderer
*/
class vvImgRend : public vvRenderer
{
  public:
    enum ActionType     /// Action types for rendering, depend on user input
    {
      RENDER=0,         ///< request for a new rendered image
      SET_TF,           ///< the transfer function has changed
      SET_FRAME,        ///< use another animation frame
      SET_INTERPOL,     ///< set interpolation mode
      SET_QUALITY,      ///< set image quality
      SET_SIZE,         ///< set object size
      SET_POSITION,     ///< set object position
      SET_EYE,          ///< set current eye (left or right)
      REQUEST_IMAGE,    ///< request most recently rendered intermediate image from remote renderer
      CLOSE_CONNECTION, ///< close socket connection
      PROFILE_START,    ///< start profiling
      PROFILE_STOP      ///< stop profiling
    };

    vvImgRend(vvVolDesc*, int, const char* = NULL);
    virtual ~vvImgRend();
    virtual void  renderVolumeGL();
    virtual void  updateTransferFunction();
    virtual void  setCurrentFrame(int);
    virtual void  setQuality(float);
    virtual void  setObjectSize(vvVector3*);
    virtual void  setPosition(vvVector3*);
    virtual void  resizeEdgeMax(float);
    virtual bool  instantClassification();
    virtual void  profileStart();
    virtual void  profileStop();
    virtual void  setParameter(float, ParameterType, char* = NULL);
    virtual float getParameter(ParameterType, char* = NULL);

  private:
    static const float TIMEOUT; ///< timeout for socket connection
    vvMatrix4 iwWarp;         ///< matrix to warp image with before rendering
    vvImage* image;           ///< intermediate image in vvImage format
    vvSoftImg* intImg;        ///< intermediate image to display
    vvSocketIO* socket;       ///< socket for connection to remote renderer
    vvStopwatch* decodeTime;  ///< time for decoding the intermediate image
    vvStopwatch* drawTime;    ///< time for drawing the intermediate image
    vvStopwatch* idleTime;    ///< idle time
    vvStopwatch* tcpTime;     ///< TCP transfer time
    bool firstImage;          ///< true before the first rendering command
    int  prevFrame;           ///< previous frame index
    float prevQuality;        ///< previous quality
    bool interpolation;       ///< interpolation mode: true=linear interpolation (default), false=nearest neighbor

    int createSocketConnection(vvSocketIO**, int, const char*);
};

#endif

//============================================================================
// End of File
//============================================================================
