//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#include <stdlib.h>
#ifdef WIN32
  #include <windows.h>
#endif
#ifndef VV_REMOTE_RENDERING
  #include <GL/gl.h>
#endif
#include <string.h>
#include <assert.h>
#include "vvrenderer.h"
#include "vvdebugmsg.h"
#include "vvprintgl.h"

/*ARGSUSED*/

//----------------------------------------------------------------------------
/** Constructor, called when program starts and when the rendering method changes
  @param voldesc volume descriptor, should be deleted when not needed anymore
*/
vvRenderer::vvRenderer(vvVolDesc* voldesc)
{
  vvDebugMsg::msg(1, "vvRenderer::vvRenderer()");
  assert(voldesc!=NULL);
  vd = voldesc;
  init();
}

//----------------------------------------------------------------------------
/// Initialization routine for class variables.
void vvRenderer::init()
{
  vvDebugMsg::msg(1, "vvRenderer::init()");

  rendererType = UNKNOWN;
  orientation = false;
  palette = false;
  boundaries = false;
  quality = 1.0f;
  boundColor[0] = boundColor[1] = boundColor[2] = 1.0f;
  probeColor[0] = 1.0f; probeColor[1] = probeColor[2] = 0.0f;
  probeSize = 0.0f;
  probePos.zero();
  clipMode = false;
  clipPoint.set(0.0f, 0.0f, 0.0f);
  clipNormal.set(0.0f, 0.0f, 1.0f);
  timing = false;

  // Determine volume size in world space:
  size.e[0] = vd->dist[0] * (float)vd->vox[0];
  size.e[1] = vd->dist[1] * (float)vd->vox[1];
  size.e[2] = vd->dist[2] * (float)vd->vox[2];
}

//----------------------------------------------------------------------------
/** Destructor, called when program ends or when rendering method changes.
   Clear up all dynamically allocated memory space here.
*/
vvRenderer::~vvRenderer()
{
  vvDebugMsg::msg(1, "vvRenderer::~vvRenderer()");
}

//----------------------------------------------------------------------------
/** Adapt image quality to frame rate.
  The method assumes a reciprocal dependence of quality and frame rate.
  The threshold region is always entered from below.
  @param quality    current image quality (>0, the bigger the better)
  @param curFPS     current frame rate [fps]
  @param desFPS     desired frame rate [fps]
  @param threshold  threshold to prevent flickering [fraction of frame rate difference]
  @return the adapted quality
*/
float vvRenderer::adaptQuality(float quality, float curFPS, float desFPS, float threshold)
{
  const float MIN_QUALITY = 0.001f;
  const float MAX_QUALITY = 10.0f;

  vvDebugMsg::msg(1, "vvRenderer::adaptQuality()");
  
  if (curFPS>0.0f && curFPS<1000.0f && desFPS>0.0f && threshold>=0.0f)  // plausibility check
  {
    if (curFPS>desFPS || ((desFPS - curFPS) / desFPS) > threshold)
    {
      quality *= curFPS / desFPS;
      if (quality < MIN_QUALITY) quality = MIN_QUALITY;
      else if (quality > MAX_QUALITY) quality = MAX_QUALITY;
    }
  }
  return quality;
}

//----------------------------------------------------------------------------
/** Returns the type of the renderer used.
*/
vvRenderer::RendererType vvRenderer::getRendererType()
{
  return rendererType;
}

//----------------------------------------------------------------------------
/** Returns a voxel value histogram array. The histogram is computed for
  all values in volume animation.
  @param buckets histogram resolution
  @param hist _allocated_ memory for 'buckets' floating point values
  @return histogram values in 'hist'
*/
void vvRenderer::makeHistogram(int buckets, float* hist)
{
  int* count;
  
  count = new int[buckets];
  vd->makeHistogram(-1, buckets, count);
  vd->normalizeHistogram(buckets, count, hist, vvVolDesc::VV_LOGARITHMIC);
  delete[] count;
}

//----------------------------------------------------------------------------
/** Core display rendering routine. 
  Should render the volume to the currently selected draw buffer. This parent
  method renders the coordinate axes and the palette, if the respective
  modes are set.
*/
void vvRenderer::renderVolumeGL()
{
  vvDebugMsg::msg(3, "vvRenderer::renderVolumeGL()");

  // Draw legend if requested:
  if (orientation==true) renderCoordinates();
  if (palette==true) renderPalette();
}

//----------------------------------------------------------------------------
/** Copy the currently displayed renderer image to a memory buffer and
  resize the image if necessary.
  Important: OpenGL canvas size should be larger than rendererd image size
  for optimum quality!
  @param w,h image size in pixels
  @param data _allocated_ memory space providing w*h*3 bytes of memory space
  @return memory space to which volume was rendered. This need not be the same
          as data, if internal space is used.
*/
void vvRenderer::renderVolumeRGB(int w, int h, uchar* data)
{
#ifndef VV_REMOTE_RENDERING
  GLint glsDrawBuffer;  // state of draw buffer
  uchar* screenshot;
  int viewPort[4];    // x, y, width, height of viewport
  int x, y;
  int srcIndex, dstIndex, srcX, srcY;
  int offX, offY;   // offsets in source image to maintain aspect ratio
  int srcWidth, srcHeight;  // actually used area of source image

  vvDebugMsg::msg(3, "vvRenderer::renderVolumeRGB(), size: ", w, h);
  glGetIntegerv(GL_VIEWPORT, viewPort);

  screenshot = new uchar[viewPort[2] * viewPort[3] * 3];
  glGetIntegerv(GL_DRAW_BUFFER, &glsDrawBuffer);  // save draw buffer
  glDrawBuffer(GL_FRONT);                         // set draw buffer to front in order to read image data   
  glReadPixels(0, 0, viewPort[2], viewPort[3], GL_RGB, GL_UNSIGNED_BYTE, screenshot); // read image data
  glDrawBuffer(glsDrawBuffer);    // restore draw buffer

  // Maintain aspect ratio:
  if ((float)viewPort[2] / (float)viewPort[3] > (float)w / (float)h)  // movie image more narrow than OpenGL window?
  {
    srcHeight = viewPort[3];
    srcWidth = srcHeight * w / h;
    offX = (viewPort[2] - srcWidth) / 2;
    offY = 0;
  }
  else    // movie image wider than OpenGL window
  {
    srcWidth = viewPort[2];
    srcHeight = h * srcWidth / w;
    offX = 0;
    offY = (viewPort[3] - srcHeight) / 2;
  }

  // Now resample image data:
  for (y=0; y<h; ++y)
    for (x=0; x<w; ++x)
    {
      dstIndex = 3 * (x + (h - y - 1) * w);
      srcX = offX + srcWidth * x / w;
      srcY = offY + srcHeight * y / h;
      srcIndex = 3 * (srcX + srcY * viewPort[2]);
      memcpy(&data[dstIndex], &screenshot[srcIndex], 3);
    }
  delete[] screenshot;
#else
  // prevent unreferenced parameter warnings:
  w = w;
  h = h;
  data = data;  
#endif
}

//----------------------------------------------------------------------------
/// Update transfer function in renderer
void vvRenderer::updateTransferFunction()
{
  vvDebugMsg::msg(1, "vvRenderer::updateTransferFunction()");
}

//----------------------------------------------------------------------------
/** Update volume data in renderer.
  This function is called when the volume data in vvVolDesc were modified.
*/
void vvRenderer::updateVolumeData()
{
  vvDebugMsg::msg(1, "vvRenderer::updateVolumeData()");
}

//----------------------------------------------------------------------------
/** Set clipping plane parameters.
  Parameters must be given in object coordinates. 
  The clipping plane will be fixed in object space.
  The clipping plane will only be used for clipping the volume. I will
  be switched off before the volume render command returns.
  @param point   arbitrary point on clipping plane
  @param normal  normal vector of clipping plane, direction points 
                 to clipped area. If length of normal is 0, the clipping plane
                 parameters are not changed.
*/
void vvRenderer::setClippingPlane(vvVector3* point, vvVector3* normal)
{
  vvDebugMsg::msg(3, "vvRenderer::setClippingPlane()");
  if (normal->length()==0.0f) return; // ignore call if normal is zero
  clipPoint.copy(point);
  clipNormal.copy(normal);
  clipNormal.normalize();
}

//----------------------------------------------------------------------------
/** Set clipping plane mode.
  @param newMode true for clipping plane on
*/
void vvRenderer::setClippingMode(bool newMode)
{
  vvDebugMsg::msg(3, "vvRenderer::setClippingMode()");
  clipMode = newMode;
}

//----------------------------------------------------------------------------
/** Get clipping plane mode.
  @return true for clipping plane on
*/
bool vvRenderer::getClippingMode()
{
  vvDebugMsg::msg(3, "vvRenderer::getClippingMode()");
  return clipMode;
}

//----------------------------------------------------------------------------
/** Set single slice clipping mode.
  @param newMode true for single slice clipping on. Clip mode must be on to
                 take this parameter into effect.
*/
void vvRenderer::setSingleSliceMode(bool newMode)
{
  vvDebugMsg::msg(3, "vvRenderer::setSingleSliceMode()");
  singleSlice = newMode;
}

//----------------------------------------------------------------------------
/** Get single slice clipping mode.
  @return true for single slice clipping on 
*/
bool vvRenderer::getSingleSliceMode()
{
  vvDebugMsg::msg(3, "vvRenderer::getSingleSliceMode()");
  return singleSlice;
}

//----------------------------------------------------------------------------
/** Returns the number of animation frames. 
*/
int vvRenderer::getNumFrames()
{
  vvDebugMsg::msg(3, "vvRenderer::getNumFrames()");
  return vd->frames;
}

//----------------------------------------------------------------------------
/** Returns index of current animation frame. 
  (first frame = 0, <0 if undefined)
*/
int vvRenderer::getCurrentFrame()
{
  vvDebugMsg::msg(3, "vvRenderer::getCurrentFrame()");
  return vd->getCurrentFrame();
}

//----------------------------------------------------------------------------
/** Set new frame index.
  @param index  new frame index (0 for first frame)
*/
void vvRenderer::setCurrentFrame(int index)
{
  vvDebugMsg::msg(3, "vvRenderer::setCurrentFrame()");
  if (index == vd->getCurrentFrame()) return;
  if (index < 0) index = 0;
  if (index >= vd->frames) index = vd->frames - 1;
  vvDebugMsg::msg(3, "New frame index: ", index);
  vd->setCurrentFrame(index);
}

//----------------------------------------------------------------------------
/** Set boundaries mode.
  @param bMode  true if boundaries are to be shown
*/
void vvRenderer::setBoundariesMode(bool bMode)
{
  vvDebugMsg::msg(3, "vvRenderer::setBoundariesMode()");
  boundaries = bMode;
}

//----------------------------------------------------------------------------
/** Set time measurement mode.
  @param newMode  true if computation times are to be displayed in text window
*/
void vvRenderer::setTimingMode(bool newMode)
{
  vvDebugMsg::msg(1, "vvRenderer::setTimingMode()");
  timing = newMode;
}

//----------------------------------------------------------------------------
/** Get time measurement mode.
  @return true if computation times are to be displayed in text window
*/
bool vvRenderer::getTimingMode()
{
  vvDebugMsg::msg(1, "vvRenderer::getTimingMode()");
  return timing;
}

//----------------------------------------------------------------------------
/** Set light parameters.
  @param numLights     number of lights in scene (0=ambient light only)
  @param stickyLights  true if lights stick to object
*/
void vvRenderer::setLights(int numLights, bool stickyLights)
{
  vvDebugMsg::msg(1, "vvRenderer::setLights()");

  numLights = numLights;
  stickyLights = stickyLights;
}

//----------------------------------------------------------------------------
/** Set rendering quality. 
  Lower quality means higher frame rates.
  @param q  image quality: range from 0.0 (highest frame rate) to infinity.
                           1.0 is 1:1 voxels:slices ratio
*/
void vvRenderer::setQuality(float q)
{
  vvDebugMsg::msg(3, "vvRenderer::setQuality()");

  if (q < 0.0f) quality = 0.0f;
  else quality = q;
}

//----------------------------------------------------------------------------
/// Returns current rendering quality.
float vvRenderer::getQuality()
{
  vvDebugMsg::msg(3, "vvRenderer::getQuality()");
  return quality;
}

//----------------------------------------------------------------------------
/// Set volume size (width, height, depth) in world coordinates
void vvRenderer::setObjectSize(vvVector3* s)
{
  vvDebugMsg::msg(1, "vvRenderer::setObjectSize()");
  size.copy(s);
}

//----------------------------------------------------------------------------
/// Get volume size (width, height, depth) in world coordinates
void vvRenderer::getObjectSize(vvVector3* s)
{
  vvDebugMsg::msg(1, "vvRenderer::getObjectSize()");
  s->copy(&size);
}

//----------------------------------------------------------------------------
/// Resize volume so that the longest edge becomes the length of len
void vvRenderer::resizeEdgeMax(float len)
{
  float maxLen;       // maximum edge length

  vvDebugMsg::msg(1, "vvRenderer::resizeEdgeMax()");

  // Determine volume dimensions in world space:
  maxLen = (float)ts_max((float)vd->vox[0] * vd->dist[0], 
                         (float)vd->vox[1] * vd->dist[1], 
                         (float)vd->vox[2] * vd->dist[2]);
  size.e[0] = vd->dist[0] * (float)vd->vox[0] / maxLen;
  size.e[1] = vd->dist[1] * (float)vd->vox[1] / maxLen;
  size.e[2] = vd->dist[2] * (float)vd->vox[2] / maxLen;
  size.scale(len);
}

//----------------------------------------------------------------------------
/// Set voxel dimensions (width, height, depth)
void vvRenderer::setVoxelSize(vvVector3* voxSize)
{
  vvDebugMsg::msg(1, "vvRenderer::setVoxelSize()");
  vd->dist[0] = voxSize->e[0];
  vd->dist[1] = voxSize->e[1];
  vd->dist[2] = voxSize->e[2];
}

//----------------------------------------------------------------------------
/// Get voxel dimensions (width, height, depth)
void vvRenderer::getVoxelSize(vvVector3* voxSize)
{
  vvDebugMsg::msg(1, "vvRenderer::getVoxelSize()");
  voxSize->e[0] = vd->dist[0];
  voxSize->e[1] = vd->dist[1];
  voxSize->e[2] = vd->dist[2];
}

//----------------------------------------------------------------------------
/** Set orientation display mode.
  @param oMode  true if axis orientation is to be shown
*/
void vvRenderer::setOrientationMode(bool oMode)
{
  vvDebugMsg::msg(1, "vvRenderer::setOrientationMode()");
  orientation = oMode;
}

//----------------------------------------------------------------------------
/** Set boundaries color.
	@param red,green,blue   RGB values [0..1] of boundaries color values
*/
void vvRenderer::setBoundariesColor(float red, float green, float blue)
{
  vvDebugMsg::msg(1, "vvRenderer::setBoundariesColor()");
	red   = ts_clamp(red,   0.0f, 1.0f);
	green = ts_clamp(green, 0.0f, 1.0f);
	blue  = ts_clamp(blue,  0.0f, 1.0f);
  boundColor[0] = red;
  boundColor[1] = green;
  boundColor[2] = blue;
}

//----------------------------------------------------------------------------
/** Get boundaries color.
	@param red,green,blue   RGB values [0..1] of boundaries color values
*/
void vvRenderer::getBoundariesColor(float* red, float* green, float* blue)
{
  vvDebugMsg::msg(1, "vvRenderer::getBoundariesColor()");
  *red   = boundColor[0];
  *green = boundColor[1];
  *blue  = boundColor[2];
}

//----------------------------------------------------------------------------
/** Set probe boundaries color.
	@param red,green,blue   RGB values [0..1] of probe boundaries color values
*/
void vvRenderer::setProbeColor(float red, float green, float blue)
{
  vvDebugMsg::msg(3, "vvRenderer::setProbeColor()");
	red   = ts_clamp(red,   0.0f, 1.0f);
	green = ts_clamp(green, 0.0f, 1.0f);
	blue  = ts_clamp(blue,  0.0f, 1.0f);
  probeColor[0] = red;
  probeColor[1] = green;
  probeColor[2] = blue;
}

//----------------------------------------------------------------------------
/** Get probe boundaries color.
	@param red,green,blue   RGB values [0..1] of probe boundaries color values
*/
void vvRenderer::getProbeColor(float* red, float* green, float* blue)
{
  vvDebugMsg::msg(1, "vvRenderer::getProbeColor()");
  *red   = probeColor[0];
  *green = probeColor[1];
  *blue  = probeColor[2];
}

//----------------------------------------------------------------------------
/** Set palette display mode.
  @param pMode  true if transfer colors palette is to be shown
*/
void vvRenderer::setPaletteMode(bool pMode)
{
  vvDebugMsg::msg(1, "vvRenderer::setPaletteMode()");
  palette = pMode;
}

//----------------------------------------------------------------------------
/** Render axis coordinates in bottom right corner.
  Arrows are of length 1.0<BR>
  Colors: x-axis=red, y-axis=green, z-axis=blue
*/
void vvRenderer::renderCoordinates()
{
#ifndef VV_REMOTE_RENDERING
  vvMatrix4 mv;                 // current modelview matrix 
  vvVector3 column;             // column vector
  GLboolean glsLighting;        // stores GL_LIGHTING
  int viewPort[4];              // x, y, width, height of viewport
  float aspect;                 // viewport aspect ratio
  float half[2];                // half viewport dimensions (x,y)
  int i;

  // Save lighting mode:
  glGetBooleanv(GL_LIGHTING, &glsLighting);
	glDisable(GL_LIGHTING); 

  // Get viewport parameters: 
  glGetIntegerv(GL_VIEWPORT, viewPort);
  aspect = (float)viewPort[2] / (float)viewPort[3];

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();

  half[0] = (aspect < 1.0f) ? 1.0f : aspect;
  half[1] = (aspect > 1.0f) ? 1.0f : (1.0f / aspect);
  glOrtho(-half[0], half[0], -half[1], half[1], 10.0f, -10.0f);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();

  // Compute modelview matrix:
  getModelviewMatrix(&mv);
  mv.killTrans();
  for (i=0; i<3; ++i)   // normalize base vectors to remove scaling
  {
    mv.getColumn(i, &column);
    column.normalize();
    mv.setColumn(i, &column);
  }    
  mv.translate(0.8f * half[0], -0.8f * half[1], 0.0f);
  mv.scale(0.2f, 0.2f, 0.2f);
  setModelviewMatrix(&mv);

  // Draw axis cross:
  glBegin(GL_LINES);
    glColor3f(1.0f, 0.0f, 0.0f);   // red
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(1.0f, 0.0f, 0.0f);

    glColor3f(0.0f, 1.0f, 0.0f);   // green
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 1.0f, 0.0f);
    
    glColor3f(0.0f, 0.0f, 1.0f);   // blue
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
  glEnd();

  // Draw arrows:
  glBegin(GL_TRIANGLES);
    glColor3f(1.0f, 0.0f, 0.0f);   // red
    glVertex3f(1.0f, 0.0f, 0.0f);
    glVertex3f(0.8f, 0.0f,-0.2f);
    glVertex3f(0.8f, 0.0f, 0.2f);

    glColor3f(0.0f, 1.0f, 0.0f);   // green
    glVertex3f(0.0f, 1.0f, 0.0f);
    glVertex3f(-0.2f, 0.8f, 0.0f);
    glVertex3f(0.2f, 0.8f, 0.0f);

    glColor3f(0.0f, 0.0f, 1.0f);   // blue
    glVertex3f(0.0f, 0.0f, 1.0f);
    glVertex3f(-0.2f,-0.0f, 0.8f);
    glVertex3f(0.2f, 0.0f, 0.8f);
  glEnd();
  
  // Restore matrix states:
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);

  // Restore lighting mode:
  if (glsLighting==(uchar)true) glEnable(GL_LIGHTING);
  else glDisable(GL_LIGHTING);
#endif
}

//----------------------------------------------------------------------------
/// Render transfer function palette at left border.
void vvRenderer::renderPalette()
{
#ifndef VV_REMOTE_RENDERING
  const int WIDTH = 10;         // palette width [pixels]
  GLfloat viewport[4];          // OpenGL viewport information (position and size)
	GLfloat glsRasterPos[4];	    // current raster position (glRasterPos)
  float* colors;                // palette colors
  uchar* image;                 // palette image
  int w, h;                     // width and height of palette
  int x, y, c;

  vvDebugMsg::msg(1, "vvRenderer::renderPalette()");

  if (vd->bpv > 2) return;    // palette only makes sense with scalar data

  // Get viewport size:
  glGetFloatv(GL_VIEWPORT, viewport);
  if (viewport[2]<=0 || viewport[3]<=0) return;   // safety first

  // Save matrix states:
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(-1.0f, 1.0f, -1.0f, 1.0f, 10.0f, -10.0f);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  // Store raster position:  
	glGetFloatv(GL_CURRENT_RASTER_POSITION, glsRasterPos);

  // Compute palette image:
  w = WIDTH;
  h = (int)viewport[3];
  colors = new float[h * 3];
  vd->tf.pins.setColorModel(vvPinList::RGB_MODEL);
  vd->tf.pins.computeFunction(h, vvPinList::RGB, colors); 
  image = new uchar[w * h * 3];
  for (x=0; x<w; ++x)
    for (y=0; y<h; ++y)
      for (c=0; c<3; ++c)
        image[c + 3 * (x + w * y)] = (uchar)(colors[c * h + y] * 255.0f);
    
  // Draw palette:
  glRasterPos2f(-1.0f,-1.0f);   // pixmap origin is bottom left corner of output window
  glDrawPixels(w, h, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)image);
  delete[] image;
  delete[] colors;

  // Display min and max values:
  vvPrintGL* printGL;
  printGL = new vvPrintGL();
  printGL->print(-0.90f,  0.9f,  "%-9.2f", vd->realMax);
  printGL->print(-0.90f, -0.95f, "%-9.2f", vd->realMin);
  delete printGL;
 
  // Restore raster position:  
	glRasterPos4fv(glsRasterPos);

  // Restore matrix states:
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
#endif
}

//----------------------------------------------------------------------------
/** Render a bounding box.
  Volume vertex names:<PRE>
      4____ 7        y
     /___ /|         |
   0|   3| |         |___x
    | 5  | /6       /
    |/___|/        z
    1    2
  </PRE>
	@param oSize  boundary box size [object coordinates]
	@param oPos   position of boundary box center [object coordinates]
	@param color  bounding box color (R,G,B) [0..1], array of 3 floats expected
  @param front  true if front boundaries are to be drawn, otherwise back boundaries are drawn
*/
void vvRenderer::drawBoundingBox(vvVector3* oSize, vvVector3* oPos, float* color, bool front)
{
#ifndef VV_REMOTE_RENDERING
  vvVector3 vertvec[8]; // vertex vectors in object space
  vvVector3 projvert[8];// projected vertex vectors
  vvVector3 edge[2];    // edge vectors
  vvVector3 normal;     // volume face normal
  vvMatrix4 mv;
  vvMatrix4 pm;
  GLboolean glsLighting;    // stores GL_LIGHTING
  float vertices[8][3] =    // volume vertices
    {{-0.5, 0.5, 0.5},
     {-0.5,-0.5, 0.5},
     { 0.5,-0.5, 0.5},
     { 0.5, 0.5, 0.5},
     {-0.5, 0.5,-0.5},
     {-0.5,-0.5,-0.5},
     { 0.5,-0.5,-0.5},
     { 0.5, 0.5,-0.5}};
  int faces[6][4] =         // volume faces. first 3 values are used for normal compuation
    {{7, 3, 2, 6},
     {0, 3, 7, 4},
     {2, 3, 0, 1},
     {4, 5, 1, 0},
     {1, 5, 6, 2},
     {6, 5, 4, 7}};
	int i,j;

  vvDebugMsg::msg(3, "vvRenderer::drawBoundingBox()");

  // Save lighting state:
  glGetBooleanv(GL_LIGHTING, &glsLighting);

  // Translate boundaries by volume position:
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();   // save modelview matrix
  glTranslatef(oPos->e[0], oPos->e[1], oPos->e[2]);

  // Get OpenGL matrices:
  getModelviewMatrix(&mv);
  getProjectionMatrix(&pm);
  
  // Create vertex vectors:
  mv.multiplyPost(&pm);
  for (i=0; i<8; ++i)
  {
    vertvec[i].set(vertices[i][0], vertices[i][1], vertices[i][2]);
    vertvec[i].scale(oSize);
    projvert[i].copy(&vertvec[i]);
    projvert[i].multiply(&mv);
  }

	// Set box color:
  glColor4f(color[0], color[1], color[2], 1.0); 

  // Disable lighting:
  glDisable(GL_LIGHTING);

  // Compute face normals and draw faces if normal points to user::
  for (i=0; i<6; ++i)
  {
    // Compute visibility of each face by computing the z coordinate
    // of the projected normals:
    edge[0].copy(&projvert[faces[i][1]]);
    edge[0].sub( &projvert[faces[i][0]]);
    edge[1].copy(&projvert[faces[i][2]]);
    edge[1].sub( &projvert[faces[i][1]]);
    normal.copy(&edge[0]);
    normal.cross(&edge[1]);

    // Draw borders of desired faces:
    if ((normal[2]>0.0f && front) || (normal[2]<0.0f && (!front)))
    {
      glBegin(GL_LINE_STRIP);  
      for (j=0; j<4; ++j)    
        glVertex3f(vertvec[faces[i][j]][0], vertvec[faces[i][j]][1], vertvec[faces[i][j]][2]);
      glVertex3f(vertvec[faces[i][0]][0], vertvec[faces[i][0]][1], vertvec[faces[i][0]][2]);
      glEnd();
    }
  }

  glPopMatrix();    // restore modelview matrix

  // Restore lighting state:  
  if (glsLighting==(uchar)true) glEnable(GL_LIGHTING);
  else glDisable(GL_LIGHTING);
#else
  front = front;
  color = color;
  oPos = oPos;
  oSize = oSize;
#endif
}

//----------------------------------------------------------------------------
/** Find out if classification can be done in real time.
  @return true if updateTransferFunction() can be processed immediately
          (eg. for color indexed textures), otherwise false is returned.
*/
bool vvRenderer::instantClassification()
{
  vvDebugMsg::msg(1, "vvRenderer::instantClassification()");
  return false;
}

//----------------------------------------------------------------------------
/// Set volume position.
void vvRenderer::setPosition(vvVector3* p)
{
  vvDebugMsg::msg(3, "vvRenderer::setPosition()");
  vd->position.copy(p);
}

//----------------------------------------------------------------------------
/// Get volume position.
void vvRenderer::getPosition(vvVector3* p)
{
  vvDebugMsg::msg(3, "vvRenderer::getPosition()");
  p->copy(&vd->position);
}

//----------------------------------------------------------------------------
/** Set the direction in which the user is currently viewing. 
  The vector originates in the user's eye and points along the
  viewing direction.
*/
void vvRenderer::setViewingDirection(vvVector3*)
{
  vvDebugMsg::msg(3, "vvRenderer::setViewingDirection()");
}

//----------------------------------------------------------------------------
/// Set the direction from the user to the object. 
void vvRenderer::setObjectDirection(vvVector3*)
{
  vvDebugMsg::msg(3, "vvRenderer::setObjectDirection()");
}

//----------------------------------------------------------------------------
/** Set the probe position.
  @param pos  position [object space]
*/
void vvRenderer::setProbePosition(vvVector3* pos)
{
  vvDebugMsg::msg(3, "vvRenderer::setProbePosition()");
  probePos.copy(pos);
}

//----------------------------------------------------------------------------
/** Get the probe position.
  @param pos  returned position [object space]
*/
void vvRenderer::getProbePosition(vvVector3* pos)
{
  vvDebugMsg::msg(3, "vvRenderer::getProbePosition()");
  pos->copy(&probePos);
}

//----------------------------------------------------------------------------
/** Set the probe size.
  @param newSize  probe size. 0.0 turns off probe draw mode
*/
void vvRenderer::setProbeSize(float newSize)
{
  vvDebugMsg::msg(3, "vvRenderer::setProbeSize()");
  probeSize = newSize;
}

//----------------------------------------------------------------------------
/** Get the probe size.
  @return probe size (0.0 = probe mode off)
*/
float vvRenderer::getProbeSize()
{
  vvDebugMsg::msg(3, "vvRenderer::getProbeSize()");
  return probeSize;
}

//----------------------------------------------------------------------------
/** Get the current modelview matrix.
  @param a matrix which will be set to the current modelview matrix
*/
void vvRenderer::getModelviewMatrix(vvMatrix4* mv)
{
#ifndef VV_REMOTE_RENDERING
  GLfloat glmatrix[16];   // OpenGL compatible matrix

  vvDebugMsg::msg(3, "vvRenderer::getModelviewMatrix()");
  glGetFloatv(GL_MODELVIEW_MATRIX, glmatrix);
  mv->getGL((float*)glmatrix);
#else
  mv->copy(&modelview);  
#endif
}

//----------------------------------------------------------------------------
/** Get the current projection matrix.
  @param a matrix which will be set to the current projection matrix
*/
void vvRenderer::getProjectionMatrix(vvMatrix4* pm)
{
  vvDebugMsg::msg(3, "vvRenderer::getProjectionMatrix()");

#ifndef VV_REMOTE_RENDERING
  GLfloat glmatrix[16];   // OpenGL compatible matrix
  glGetFloatv(GL_PROJECTION_MATRIX, glmatrix);
  pm->getGL((float*)glmatrix);
#else
  pm->copy(&projection);
#endif
}

//----------------------------------------------------------------------------
/** Set the OpenGL modelview matrix.
  @param new OpenGL modelview matrix
*/
void vvRenderer::setModelviewMatrix(vvMatrix4* mv)
{
  vvDebugMsg::msg(3, "vvRenderer::setModelviewMatrix()");

#ifndef VV_REMOTE_RENDERING
  GLfloat glmatrix[16];   // OpenGL compatible matrix
  mv->makeGL((float*)glmatrix);
  glMatrixMode(GL_MODELVIEW);
  glLoadMatrixf(glmatrix);
#endif
  modelview.copy(mv);
}

//----------------------------------------------------------------------------
/** Set the OpenGL projection matrix.
  @param new OpenGL projection matrix
*/
void vvRenderer::setProjectionMatrix(vvMatrix4* pm)
{
  vvDebugMsg::msg(3, "vvRenderer::setProjectionMatrix()");

#ifndef VV_REMOTE_RENDERING
  GLfloat glmatrix[16];   // OpenGL compatible matrix
  pm->makeGL((float*)glmatrix);
  glMatrixMode(GL_PROJECTION);
  glLoadMatrixf(glmatrix);
  glMatrixMode(GL_MODELVIEW);
#endif
  projection.copy(pm);
}

//----------------------------------------------------------------------------
/** Compute user's eye position.
  @param eye  vector to receive eye position [world space]
*/
void vvRenderer::getEyePosition(vvVector3* eye)
{
  vvMatrix4 invPM;            // inverted projection matrix
  vvVector4 projEye;          // eye x PM

  vvDebugMsg::msg(3, "vvRenderer::getEyePosition()");

  getProjectionMatrix(&invPM);
  invPM.invert();
  projEye.set(0.0f, 0.0f, -1.0f, 0.0f);
  projEye.multiply(&invPM);
  eye->copy(&projEye);
}

//----------------------------------------------------------------------------
/** Find out if user is inside of volume.
  @param point  point to test [object coordinates]
  @return true if the given point is inside or on the volume boundaries.
*/
bool vvRenderer::isInVolume(vvVector3* point)
{
  vvVector3 size2;    // half object size
  vvVector3 pos;      // object location
  int i;

  vvDebugMsg::msg(3, "vvRenderer::isInVolume()");

  pos.copy(&vd->position);
  size2.copy(&size);
  size2.scale(0.5f);
  for (i=0; i<3; ++i)
  {
    if (point->e[i] < (pos.e[i] - size2.e[i])) return false;
    if (point->e[i] > (pos.e[i] + size2.e[i])) return false;
  }
  return true;
}

//----------------------------------------------------------------------------
/** Gets the alpha value nearest to the point specified in x, y and z 
  coordinates with consideration of the alpha transfer function.
  @param x,y,z  point to test [object coordinates]
  @return normalized alpha value or -1.0 if point is outside of volume
*/
float vvRenderer::getAlphaValue(float x, float y, float z)
{
  vvVector3 size2;        // half object size
  vvVector3 point(x,y,z); // point to get alpha value at
  vvVector3 pos;
  float index;            // floating point index value into alpha TF [0..1]
  int vp[3];              // position of nearest voxel to x/y/z [voxel space]
  int i;
  uchar* ptr;

  vvDebugMsg::msg(3, "vvRenderer::getAlphaValue()");

  size2.copy(&size);
  size2.scale(0.5f);
  pos.copy(&vd->position);

  for (i=0; i<3; ++i)
  {
    if (point.e[i] < (pos.e[i] - size2.e[i])) return -1.0f;
    if (point.e[i] > (pos.e[i] + size2.e[i])) return -1.0f;

    vp[i] = int(float(vd->vox[i]) * (point.e[i] - pos.e[i] + size2.e[i]) / size.e[i]);
	  vp[i] = ts_clamp(vp[i], 0, vd->vox[i]-1);
  }

  vp[1] = vd->vox[1] - vp[1] - 1;
  vp[2] = vd->vox[2] - vp[2] - 1;
  ptr = vd->getRaw(getCurrentFrame()) + vd->bpv * (vp[0] + vp[1] * vd->vox[0] + vp[2] * vd->vox[0] * vd->vox[1]);

  // Determine index into alpha LUT:
  switch (vd->bpv)
  {
    default:
    case 1: index = float(*ptr) / 255.0f; break;
	  case 2: index = (float(*(ptr+1)) + float(int(*ptr) << 8)) / 65535.0f; break;
	  case 3: index = (float(*ptr) + float(*(ptr+1)) + float(*(ptr+2))) / (3.0f * 255.0f); break;
	  case 4: index = float(*(ptr+3)) / 255.0f; break;
  }

  // Determine alpha value:  
  return vd->tf.pins.getAlpha(index);
}

//----------------------------------------------------------------------------
/** Set a new value for a parameter. This function can be interpreted
    differently by every actual renderer implementation.
  @param newValue   new value
  @param param      parameter to change
  @param objName    name of the object to change (default: NULL)
*/
void vvRenderer::setParameter(float, ParameterType, char*)
{
  vvDebugMsg::msg(3, "vvRenderer::setParameter()");
}

//----------------------------------------------------------------------------
/** Get a parameter value.
  @param param    parameter to get value of
  @param objName  name of the object to get value of (default: NULL)
*/
float vvRenderer::getParameter(ParameterType, char*)
{
  vvDebugMsg::msg(3, "vvRenderer::getParameter()");
  return -VV_FLT_MAX;
}

//----------------------------------------------------------------------------
/// Start benchmarking.
void vvRenderer::profileStart()
{
  vvDebugMsg::msg(1, "vvRenderer::profileStart()");
}

//----------------------------------------------------------------------------
/// Stop benchmarking.
void vvRenderer::profileStop()
{
  vvDebugMsg::msg(1, "vvRenderer::profileStop()");
}

//============================================================================
// End of File
//============================================================================
