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
  using std::ios;
#else
  #include <iostream.h>
#endif

#include <assert.h>
#include "vvdebugmsg.h"
#include "vvimgrend.h"

const float vvImgRend::TIMEOUT = 600.0f;

//----------------------------------------------------------------------------
/** Constructor.
  @param vd           volume description
  @param port         port number for socket connection
*/
vvImgRend::vvImgRend(vvVolDesc* vd, int port, const char* host) : vvRenderer(vd)
{
  vvDebugMsg::msg(1, "vvImgRend::vvImgRend()");

  rendererType = IMGREND;
  intImg  = new vvSoftImg();
  image   = new vvImage();
  socket  = NULL;
  decodeTime = NULL;
  drawTime   = NULL;
  idleTime   = NULL;
  tcpTime    = NULL;
  firstImage = true;
  prevFrame = -1;
  prevQuality = -1.0f;
  interpolation = true;

  assert(createSocketConnection(&socket, port, host)==0);
 
  assert(socket->putVolume(vd) == vvSocket::VV_OK);

  updateTransferFunction();
}

//----------------------------------------------------------------------------
/// Destructor
vvImgRend::~vvImgRend()
{
  vvDebugMsg::msg(1, "vvImgRend::~vvImgRend()");
  
  // Send ID for close connection:
  int id = CLOSE_CONNECTION;
  socket->putData(&id, 1, vvSocketIO::VV_INT);
  
  // Clear memory:
  delete socket;
  delete image;
  delete intImg;
  delete decodeTime;
  delete drawTime;
  delete idleTime;
  delete tcpTime;
}

//----------------------------------------------------------------------------
/** Render 2D image to current OpenGL viewport, but don't use the OpenGL
  modelview matrix but the class's warp matrix.
*/
void vvImgRend::renderVolumeGL()
{
  vvMatrix4 mv;
  vvMatrix4 pm;
  vvVector3 pos;          // object location
  int id;

  vvDebugMsg::msg(1, "vvImgRend::renderVolumeGL()");

  // Don't draw the first time of rendering, so that afterwards 
  // always the image from the last call is rendered (pipelining).
  if (!firstImage)
  {
    id = REQUEST_IMAGE;
    if (tcpTime) tcpTime->start();
    if (socket->putData(&id, 1, vvSocketIO::VV_INT) != vvSocket::VV_OK)
    {
      cerr << "Cannot send image request command to socket." << endl;
      return;
    }
    if (tcpTime) tcpTime->stop();

    // Receive warp matrix:
    if (idleTime) idleTime->start();
    if (socket->getMatrix(&iwWarp) != vvSocket::VV_OK)
    {
      cerr << "Cannot get warp matrix from socket." << endl;
      return;
    }
    if (idleTime) idleTime->stop();

    // Receive intermediate image:
    if (tcpTime) tcpTime->start();
    if (socket->getImage(image) != vvSocket::VV_OK)
    {
      cerr << "Cannot get intermediate image from socket." << endl;
      return;
    }
    if (tcpTime) tcpTime->stop();

    // Decode intermediate image:
    if (decodeTime) decodeTime->start();
    image->decode();
    if (decodeTime) decodeTime->stop();
    intImg->setImageData(image->getWidth(), image->getHeight(), image->getImagePtr());

    if (drawTime) drawTime->start();

    // Draw back part of bounding box:
    pos.copy(&vd->position);
    if (boundaries) 
      drawBoundingBox(&size, &pos, boundColor, false);   // draw back boundaries

    // Draw warped intermediate image:
    intImg->warpTex(&iwWarp);

    // Draw front part of bounding box:
    if (boundaries) 
      drawBoundingBox(&size, &pos, boundColor, true);   // draw front boundaries

    vvRenderer::renderVolumeGL();         // draw coordinate axes

    if (drawTime) drawTime->stop();
  }
  else firstImage = false;

  // Prepare next image:

  // Send ID for rendering:
  if (tcpTime) tcpTime->start();
  id = RENDER;
  if (socket->putData(&id, 1, vvSocketIO::VV_INT) != vvSocket::VV_OK)
  {
    cerr << "Cannot send rendering command to socket." << endl;
    return;
  }
  if (tcpTime) tcpTime->stop();

  // Send modelview matrix:
  getModelviewMatrix(&mv);
  if (tcpTime) tcpTime->start();
  if (socket->putMatrix(&mv) != vvSocket::VV_OK)
  {
    cerr << "Cannot send modelview matrix to socket." << endl;
    return;
  }
  if (tcpTime) tcpTime->stop();

  // Send projection matrix:
  getProjectionMatrix(&pm);
  if (tcpTime) tcpTime->start();
  if (socket->putMatrix(&pm) != vvSocket::VV_OK)
  {
    cerr << "Cannot send projection matrix to socket." << endl;
    return;
  }
  if (tcpTime) tcpTime->stop();
  
  if (vvDebugMsg::isActive(3))
  {
    mv.print("mv sent:");
    pm.print("pm sent:");
  }
}

//----------------------------------------------------------------------------
// see parent for comments
void vvImgRend::updateTransferFunction()
{
  float* pinArray;    // serialized pins (float)
  int pinArraySize;   // number of floats in pin array
  int id = SET_TF;

  vvDebugMsg::msg(1, "vvImgRend::updateTransferFunction()");

  // Send ID for new transfer function:
  if (socket->putData(&id, 1, vvSocketIO::VV_INT) != vvSocket::VV_OK)
  {
    cerr << "Cannot send TF ID to socket" << endl;
    return;
  }

  // Serialize transfer function:
  pinArraySize = vvPin::ELEMENTS_PER_PIN * vd->tf.pins.count();
  pinArray = new float[pinArraySize];
  vd->tf.pins.setColorModel(vvPinList::RGB_MODEL);
  vd->tf.pins.makeArray(pinArray);

  // Send transfer function:
  if (socket->putData(&pinArraySize, 1, vvSocketIO::VV_INT) != vvSocket::VV_OK)
  {
    cerr << "Cannot send pin array size to socket" << endl;
    return;
  }

  if (socket->putData((void*)pinArray, pinArraySize, vvSocketIO::VV_FLOAT) != vvSocket::VV_OK)
  {
    cerr << "Cannot send TF to socket" << endl;
    return;
  }
  delete[] pinArray;
}

//----------------------------------------------------------------------------
// see parent for comments
void vvImgRend::setCurrentFrame(int newFrame)
{
  vvDebugMsg::msg(1, "vvImgRend::setCurrentFrame()");

  int id = SET_FRAME;

  if (newFrame==prevFrame) return;
  prevFrame = newFrame;

  if (socket->putData(&id, 1, vvSocketIO::VV_INT) != vvSocket::VV_OK)
  {
    cerr << "Cannot send frame ID to socket" << endl;
    return;
  }

  if (socket->putData(&newFrame, 1, vvSocketIO::VV_INT) != vvSocket::VV_OK)
  {
    cerr << "Cannot send frame to socket" << endl;
    return;
  }
}

//----------------------------------------------------------------------------
/** Establish a socket connection as a server.
  @param sock socket structure to use
  @param port port number to use for the connection
  @return 0 if ok
*/
int vvImgRend::createSocketConnection(vvSocketIO** sock, int port, const char* host)
{
  vvDebugMsg::msg(1, "vvImgRend::createSocketConnection()");

  if (host)
  {
    *sock = new vvSocketIO(port, (char*)host, vvSocketIO::VV_TCP);
    (*sock)->set_sock_param(TIMEOUT, TIMEOUT);
    cerr << "Opening client socket." << endl;
  }
  else
  {
    *sock = new vvSocketIO(port, vvSocketIO::VV_TCP);
    (*sock)->set_sock_param(TIMEOUT, TIMEOUT);
    cerr << "Opening server socket." << endl;
  }
 
  switch ((*sock)->init())
  {
    case vvSocket::VV_OK:            
      cerr<<"Socket was opened successfully"<<endl;  
      return 0;
    case vvSocket::VV_CREATE_ERROR:  
      cerr<<"Socket cannot be created"<<endl; 
      return -1;
    case vvSocket::VV_TIMEOUT_ERROR: 
      cerr<<"Timeout"<<endl; 
      return -2; 
    default:                        
      cerr<<"Cannot open socket"<<endl; 
      return -3; 
  }        
}

//----------------------------------------------------------------------------
// see parent for comment
void vvImgRend::setParameter(float newValue, ParameterType param, char*)
{
  int iMode;
  int id, eye;

  vvDebugMsg::msg(3, "vvImgRend::setParameter()");
  switch (param)
  {
    case vvRenderer::VV_SLICEINT:
      id = SET_INTERPOL;
      interpolation = (newValue == 0.0f) ? false : true;
      if (socket->putData(&id, 1, vvSocketIO::VV_INT) != vvSocket::VV_OK)
      {
        cerr << "Cannot send data to socket" << endl;
        return;
      }

      iMode = (interpolation) ? 1 : 0;
      if (socket->putData(&iMode, 1, vvSocketIO::VV_INT) != vvSocket::VV_OK)
      {
        cerr << "Cannot send data to socket" << endl;
        return;
      }
      break;
    case vvRenderer::VV_EYE:
      eye = int(newValue);
      id = SET_EYE;
      if (tcpTime) tcpTime->start();
      if (socket->putData(&id, 1, vvSocketIO::VV_INT) != vvSocket::VV_OK)
      {
        cerr << "Cannot send eye ID to socket." << endl;
        return;
      }
      if (socket->putData(&eye, 1, vvSocketIO::VV_INT) != vvSocket::VV_OK)
      {
        cerr << "Cannot send eye info to socket." << endl;
        return;
      }
      if (tcpTime) tcpTime->stop();
      break;
    default: break;
  }
}

//----------------------------------------------------------------------------
// see parent for comment
float vvImgRend::getParameter(ParameterType param, char*)
{
  vvDebugMsg::msg(3, "vvImgRend::getParameter()");

  switch (param)
  {
    case vvRenderer::VV_SLICEINT:
      return (interpolation) ? 1.0f : 0.0f;
      break;
    default: return vvRenderer::getParameter(param);
  }
}

//----------------------------------------------------------------------------
/** Additionally to the parent function, the quality is sent to the
  remote renderer. If the quality hasn't changed since the last call,
  nothing is done.
  @param newQuality new quality value (>0)
*/
void vvImgRend::setQuality(float newQuality)
{
  vvDebugMsg::msg(3, "vvImgRend::setQuality()", newQuality);

  int id = SET_QUALITY;

  if (prevQuality == newQuality) return;  // omit unnecessary network traffic
  
  if (socket->putData(&id, 1, vvSocketIO::VV_INT) != vvSocket::VV_OK)
  {
    cerr << "Cannot send quality ID to socket" << endl;
    return;
  }

  if (socket->putData(&newQuality, 1, vvSocketIO::VV_FLOAT) != vvSocket::VV_OK)
  {
    cerr << "Cannot send quality to socket" << endl;
    return;
  }

  prevQuality = newQuality;
  vvRenderer::setQuality(newQuality);
}

//----------------------------------------------------------------------------
/** Set volume size (width, height, depth) in world coordinates and
  pass the new size to the remote renderer.
*/
void vvImgRend::setObjectSize(vvVector3* s)
{
  vvDebugMsg::msg(1, "vvImgRend::setObjectSize()");
  
  int id = SET_SIZE;

  vvRenderer::setObjectSize(s);

  if (socket->putData(&id, 1, vvSocketIO::VV_INT) != vvSocket::VV_OK)
  {
    cerr << "Cannot send size ID to socket" << endl;
    return;
  }

  if (socket->putData(s->e, 3, vvSocketIO::VV_FLOAT) != vvSocket::VV_OK)
  {
    cerr << "Cannot send size to socket" << endl;
    return;
  }
}

//----------------------------------------------------------------------------
/** Set volume position in world coordinates and
  pass the new position to the remote renderer.
*/
void vvImgRend::setPosition(vvVector3* p)
{
  vvDebugMsg::msg(1, "vvImgRend::setPosition()");
  
  int id = SET_POSITION;

  vvRenderer::setPosition(p);

  if (socket->putData(&id, 1, vvSocketIO::VV_INT) != vvSocket::VV_OK)
  {
    cerr << "Cannot send position ID to socket" << endl;
    return;
  }

  if (socket->putData(p->e, 3, vvSocketIO::VV_FLOAT) != vvSocket::VV_OK)
  {
    cerr << "Cannot send position to socket" << endl;
    return;
  }
}

//----------------------------------------------------------------------------
/// Resize volume so that the longest edge becomes the length of len
void vvImgRend::resizeEdgeMax(float len)
{
  vvDebugMsg::msg(1, "vvImgRend::resizeEdgeMax()");
  vvRenderer::resizeEdgeMax(len);
  setObjectSize(&size);
}

//----------------------------------------------------------------------------
/// This renderer allows instant classification.
bool vvImgRend::instantClassification()
{
  vvDebugMsg::msg(1, "vvImgRend::instantClassification()");
  return true;
}

//----------------------------------------------------------------------------
/// Start benchmarking.
void vvImgRend::profileStart()
{
  vvDebugMsg::msg(1, "vvImgRend::profileStart()");
  
  int id = PROFILE_START;

  if (socket->putData(&id, 1, vvSocketIO::VV_INT) != vvSocket::VV_OK)
  {
    cerr << "Cannot send profile start ID to socket" << endl;
    return;
  }
  decodeTime = new vvStopwatch();
  drawTime   = new vvStopwatch();
  idleTime   = new vvStopwatch();
  tcpTime    = new vvStopwatch();
}

//----------------------------------------------------------------------------
/// Stop benchmarking.
void vvImgRend::profileStop()
{
  vvDebugMsg::msg(1, "vvImgRend::profileStop()");
  
  int id = PROFILE_STOP;
  
  if (socket->putData(&id, 1, vvSocketIO::VV_INT) != vvSocket::VV_OK)
  {
    cerr << "Cannot send profile stop ID to socket" << endl;
    return;
  }
  
  // Print profiling results:
  cerr.setf(ios::fixed, ios::floatfield);
  cerr.precision(3);
  cerr << "Intermediate image size [pixels].................." << intImg->width << " x " << intImg->height << endl;
  cerr << "Total image decoding time [sec]..................." << decodeTime->getTime() << endl;
  cerr << "Total drawing time [sec].........................." << drawTime->getTime() << endl;
  cerr << "Total idle time [sec]............................." << idleTime->getTime() << endl;
  cerr << "Total TCP transfer time [sec]....................." << tcpTime->getTime() << endl;

  // Free stopwatch memory:
  delete decodeTime;
  delete drawTime;
  delete idleTime;
  delete tcpTime;
  decodeTime = NULL;
  drawTime   = NULL;
  idleTime   = NULL;
  tcpTime    = NULL;
}

//============================================================================
// End of File
//============================================================================
