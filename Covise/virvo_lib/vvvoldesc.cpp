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

#ifdef WIN32
  #include <float.h>
#endif
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "vvdebugmsg.h"
#include "vvtoolshed.h"
#include "vvvecmath.h"
#include "vvvoldesc.h"
#include "vvsllist.h"

#ifdef __sun
  #define logf log
#endif

//============================================================================
// Class vvVolDesc
//============================================================================

//----------------------------------------------------------------------------
/// Constructor.
vvVolDesc::vvVolDesc()
{
  vvDebugMsg::msg(1, "vvVolDesc::vvVolDesc(1)");
  initialize();
}

//----------------------------------------------------------------------------
/** Constructor with filename initialization.
 @param fn  volume file name
*/
vvVolDesc::vvVolDesc(const char* fn)
{
  vvDebugMsg::msg(1, "vvVolDesc::vvVolDesc(2)");
  initialize();
  setFilename(fn);
}

//----------------------------------------------------------------------------
/** Constructor with volume data initialization.
  The volume data will not be replicated and not be deleted upon deletion of
  vvVolDesc, so it _must_ be deleted by the caller!
 @param fn  volume file name (use "COVISE" if source is COVISE)
 @param w   width in pixels
 @param h   height in pixels
 @param s   number of volume slices
 @param f   number of animation frames
 @param b   number of bytes per voxel
 @param d   pointer to pointer array of raw voxel data
*/
vvVolDesc::vvVolDesc(const char* fn, int w, int h, int s, int f, int b, 
  uchar** d)
{
  uchar* data;
  int i;
  
  vvDebugMsg::msg(1, "vvVolDesc::vvVolDesc(3)");
  initialize();
  setFilename(fn);
  vox[0] = w;
  vox[1] = h;
  vox[2] = s;
  bpv    = b;
  frames = f;
  
  for (i=0; i<f; ++i)
  {
    // Replicate data because later on Covise conversion happens:
    data = new uchar[getFrameSize()];
    memcpy(data, d[i], getFrameSize());
    addFrame(data, DELETE_DATA);
  }
  if (strcmp("COVISE", fn)==0)				// convert data if source is COVISE
	  convertCoviseToVirvo();
}

//----------------------------------------------------------------------------
/** Constructor with volume data initialization 
    in floating point format for density only data.
 @param fn     file name (use "COVISE" if source is COVISE)
 @param w,h,s  width, height, slices
 @param f      number of animation frames
 @param d      pointer array to raw voxel data (must be deleted by caller!)
*/
vvVolDesc::vvVolDesc(const char* fn, int w, int h, int s, int f, float** d)
{
  uchar* data;
  int i;
  float frameMin, frameMax;     // minimum and maximum value in one frame

  vvDebugMsg::msg(1, "vvVolDesc::vvVolDesc(4)");
  initialize();
  vox[0] = w;
  vox[1] = h;
  vox[2] = s;
  bpv    = 1;
  frames = f;
  realMin = VV_FLT_MAX;
  realMax = -realMin;
  for (i=0; i<f; ++i)
  {
    vvToolshed::getMinMaxIgnore(d[i], getFrameVoxels(), VV_FLT_MAX, &frameMin, &frameMax);
    if (frameMin < realMin) realMin = frameMin;
    if (frameMax > realMax) realMax = frameMax;
  }
  for (i=0; i<f; ++i)
  {
    data = new uchar[getFrameSize()];
    vvToolshed::convertFloat2UCharClampZero(d[i], data, getFrameVoxels(), 
      realMin, realMax, VV_FLT_MAX);
    addFrame(data, DELETE_DATA);
  }
  if (strcmp("COVISE", fn)==0)				// convert data if source is COVISE
  	convertCoviseToVirvo();
}

//----------------------------------------------------------------------------
/** Constructor with volume data initialization for one timestep
 in floating point format for RGB data
 @param fn     file name (use "COVISE" if source is COVISE)
 @param w,h,s  width, height, slices
 @param f      number of animation frames
 @param r,g,b  pointer arrays to raw voxel data for red, green, blue [0.0f..1.0f]
               (must be deleted by caller!)
*/
vvVolDesc::vvVolDesc(const char* fn, int w, int h, int s, int f, float** r, 
  float** g, float** b)
{
  uchar* data;
  float* c[3];
  int i,j,k;
  int frameSize;    // number of bytes per volume animation frame
  int compSize;     // size of one color component

  vvDebugMsg::msg(1, "vvVolDesc::vvVolDesc(5)");
  initialize();
  vox[0] = w;
  vox[1] = h;
  vox[2] = s;
  bpv    = 3;
  frames = f;
  frameSize = getFrameSize();
  compSize  = getFrameVoxels();

  // Convert float to uchar:
  for (k=0; k<f; ++k)
  {
    c[0] = r[k];
    c[1] = g[k];
    c[2] = b[k];
    data = new uchar[frameSize];
    for (j=0; j<3; ++j)
      for (i=0; i<compSize; ++i)
        data[i * 3 + j] = (uchar)(255.0f * c[j][i]);
    addFrame(data, DELETE_DATA);
  }
  if (strcmp("COVISE", fn)==0)				// convert data if source is COVISE
	  convertCoviseToVirvo();
}

//----------------------------------------------------------------------------
/** Constructor with volume data initialization for a 2D image in RGB format.
 @param fn   image file name
 @param w,h  width and height of image in pixels
 @param d    pointer to raw voxel data (must be deleted by caller!)
*/
vvVolDesc::vvVolDesc(const char* fn, int w, int h, uchar* d)
{
  uchar* data;

  vvDebugMsg::msg(1, "vvVolDesc::vvVolDesc(6)");
  initialize();
  setFilename(fn);
  vox[0] = w;
  vox[1] = h;
  vox[2] = 1;
  bpv    = 3;
  frames = 1;
  data = new uchar[getFrameSize()];
  memcpy(data, d, getFrameSize());
  addFrame(data, DELETE_DATA);
}

//----------------------------------------------------------------------------
/** Copy constructor.
 Copies all vvVolDesc data including transfer functions and raw voxel data.
 @param v  source volume description
 @param f  frame index to copy (-1 for all frames, -2 for no copying of raw data)
*/
vvVolDesc::vvVolDesc(vvVolDesc* v, int f)
{
  int i;
  
  vvDebugMsg::msg(3, "vvVolDesc::vvVolDesc(7)");
  initialize();
  setFilename(v->filename);
  vox[0] = v->vox[0];
  vox[1] = v->vox[1];
  vox[2] = v->vox[2];
  bpv    = v->bpv;
  frames = 0;
  currentFrame = (f==-1) ? v->currentFrame : 0;
  for (i=0; i<3; ++i)
    dist[i] = v->dist[i];
  dt = v->dt;

  // Copy icon:
  iconSize = v->iconSize;
  if (iconSize > 0)
  {
    iconData = new uchar[iconSize * 3];
    memcpy(iconData, v->iconData, iconSize * 3);
  }

  // Copy transfer functions:
  tf.copy(&v->tf);

  if (f==-1) 
  {
    for (i=0; i<v->frames; ++i)
    {
      copyFrame(v->getRaw(i));      
      ++frames;
    }
  }
  else if (f!=-2)
  {
    copyFrame(v->getRaw(f));
    ++frames;
  }
}

//----------------------------------------------------------------------------
/// Destructor
vvVolDesc::~vvVolDesc()
{
  vvDebugMsg::msg(1, "vvVolDesc::~vvVolDesc()");
  removeSequence();
  delete[] filename;
  delete[] iconData;
}

//----------------------------------------------------------------------------
/** Initialization of class attributes. 
  May only be called by constructors!
*/
void vvVolDesc::initialize()
{
  vvDebugMsg::msg(2, "vvVolDesc::initialize()");
  filename = NULL;
  iconData = NULL;
  simulate = false;
  setDefaults();
  removeSequence();
  currentFrame = 0;
  delete[] filename;
  filename = NULL;
  tf.clear();
  iconSize = 0;
  if (iconData!=NULL)
  {
    delete[] iconData;
    iconData = NULL;
  }
}

//----------------------------------------------------------------------------
/// Set default values for all serializable attributes.
void vvVolDesc::setDefaults()
{
  vvDebugMsg::msg(2, "vvVolDesc::setDefaults()");

  vox[0] = vox[1] = vox[2] = frames = 0;
  frames = 0;
  bpv = 1;
  stype = VV_INTEGER;
  dist[0] = dist[1] = dist[2] = 1.0f;
  dt = 1.0f;
  realMin = 0.0f;
  realMax = 1.0f;
  position.zero();
}

//----------------------------------------------------------------------------
/// Remove all frames of animation sequence.
void vvVolDesc::removeSequence()
{
  vvDebugMsg::msg(2, "vvVolDesc::removeSequence()");
  if (raw.isEmpty()) return;
  raw.removeAll();
}

//----------------------------------------------------------------------------
/** Get slice size.
 @return slice size in number of bytes
*/
int vvVolDesc::getSliceSize()
{
  return vox[0] * vox[1] * bpv;
}

//----------------------------------------------------------------------------
/** Get frame size.
 @return frame size in number of bytes
*/
int vvVolDesc::getFrameSize()
{
  return vox[0] * vox[1] * vox[2] * bpv;
}

//----------------------------------------------------------------------------
/** Get movie size.
 @return movie size in bytes (movie = timely sequence of all volume frames)
*/
int vvVolDesc::getMovieSize()
{
  return vox[0] * vox[1] * vox[2] * frames * bpv;
}

//----------------------------------------------------------------------------
/// Get number of voxels in a slice.
int vvVolDesc::getSliceVoxels()
{
  return vox[0] * vox[1];
}

//----------------------------------------------------------------------------
/// Get number of voxels in a frame.
int vvVolDesc::getFrameVoxels()
{
  return vox[0] * vox[1] * vox[2];
}

//----------------------------------------------------------------------------
/// Get number of voxels in the volume movie.
int vvVolDesc::getMovieVoxels()
{
  return vox[0] * vox[1] * vox[2] * frames;
}

//----------------------------------------------------------------------------
/** Merge two volume datasets. 
 The data will be moved to the target volume.<BR>
 If the volumes have the same data type (bpv), two cases are considered:
 <OL>
  <LI>The volumes have equal frame sizes
      => the new volume sequence is concatenated to the current volume sequence.</LI>
  <LI>Both volumes consist of one frame, and they share the slice sizes
      => the new volume frame is added at the back of the current volume</LI>
 </OL>
 The source volume sequence is moved to the end of the current sequence.<BR>
 No problems occur if either one of the sequences is empty.<BR>
 @param src new volume sequence. Will end up with no volume data.
 @return OK if successful
*/
vvVolDesc::ErrorType vvVolDesc::merge(vvVolDesc* src)
{
  uchar* newRaw;
  uchar* rd;
  int f, i;

  vvDebugMsg::msg(2, "vvVolDesc::merge()");
  if (src->frames==0)   // is source src empty?
    return OK;

  if (bpv != src->bpv && frames != 0)   // are data types the same?
    return TYPE_ERROR;
  
  if (frames==0)    // is the target VD empty?
  {
    // Copy all volume data to this volume:
    vox[0] = src->vox[0];
    vox[1] = src->vox[1];
    vox[2] = src->vox[2];
    frames = src->frames;
    bpv    = src->bpv;
    stype  = src->stype;
    dt     = src->dt;
    raw.merge(&src->raw);
    for (i=0; i<3; ++i)
      dist[i] = src->dist[i];
    realMin = src->realMin;
    realMax = src->realMax;
    position.copy(&src->position);
    currentFrame = src->currentFrame;
    tf.copy(&src->tf);

    // Delete sequence information from source src:
    src->bpv = src->vox[0] = src->vox[1] = src->vox[2] = src->frames = src->currentFrame = 0;
    return OK;
  }
  
  if (vox[0] == src->vox[0] && vox[1] == src->vox[1] && 
    vox[2] == src->vox[2] && vox[2] > 1)    // do volume dimensions match?
  { 
    // Append all volume time steps to target volume:
    raw.merge(&src->raw);
    frames = raw.count();

    // Delete sequence information from source src:
    src->bpv = src->vox[0] = src->vox[1] = src->vox[2] = src->frames = src->currentFrame = 0;
    return OK;
  }
  
  if (vox[0] == src->vox[0] && vox[1] == src->vox[1] && 
    frames == src->frames)    // do slice dimensions and animation length match?
  {
    // Append all slices of each animation step to target volume:
    for (f=0; f<frames; ++f)
    {
      raw.makeCurrent(f);
      rd = raw.getData();
      newRaw = new uchar[getFrameSize() + src->getFrameSize()];
      memcpy(newRaw, rd, getFrameSize()); // copy current frame to new raw data array      
      memcpy(newRaw + getFrameSize(), src->getRaw(f), src->getFrameSize());  // copy source frame to new raw data array
      raw.remove();
      if (f==0) raw.insertBefore(newRaw, true);
      else raw.insertAfter(newRaw, true);
    }
    vox[2] += src->vox[2];    // update target slice number
    src->removeSequence();    // delete copied frames from source
    return OK;
  }
  
  return TYPE_ERROR;   // more sophisticated merging methods can be implemented later on
}


//----------------------------------------------------------------------------
/** Returns a pointer to the raw data of a specific frame.
  @param frame  index of desired frame (0 for first frame) if frame does not
                exist, NULL will be returned        
*/
uchar* vvVolDesc::getRaw(int frame)
{
  if (frame<0 || frame>=frames) return NULL;    // frame does not exist
  raw.makeCurrent(frame);
  return raw.getData();
}

//----------------------------------------------------------------------------
/// Return the pointer to the raw data of the current frame.
uchar* vvVolDesc::getRaw()
{
  return getRaw(currentFrame);
}

//----------------------------------------------------------------------------
/** Adds a new frame to the animation sequence. The data has to be in
    the appropriate format, according to the vvVolDesc::bpv setting.
    The data itself is not copied to the sequence. 
    It has to be deleted by the caller if deleteData is false,
    or is deleted automatically if deleteData is true.
    <BR>
    The frames variable is not adjusted, this must be done separately.
    <BR>
    The data bytes are stored in the following order: 
    Starting point is the top left voxel on the 
    front slice. For each voxel, bpv bytes are stored with
    the most significant byte first. The voxels are stored
    in English writing order as in a book (first left to right, 
    then top to bottom, then front to back).
    After the front slice is stored, the second one is stored in the
    same order of voxels. Last stored is the right bottom voxel of
    the back slice.
  @param ptr          pointer to raw data
  @param deleteData   data deletion type: delete or don't delete when not used anymore
*/
void vvVolDesc::addFrame(uchar* ptr, DeleteType deleteData)
{
  raw.append(ptr, (deleteData==DELETE_DATA) ? true : false);
}

/// Return the number of frames actually stored.
int vvVolDesc::getStoredFrames()
{
  return raw.count();
}

//----------------------------------------------------------------------------
/** Copies frame data from memory and adds them as a new frame 
    to the animation sequence. 
  @param ptr  pointer to raw source data, _must_ be deleted by the caller!
*/
void vvVolDesc::copyFrame(uchar* ptr)
{
  uchar* newData;

  vvDebugMsg::msg(3, "vvVolDesc::copyFrame()");
  newData = new uchar[getFrameSize()];
  memcpy(newData, ptr, getFrameSize());
  raw.append(newData, true);
}

//----------------------------------------------------------------------------
/** Normalize the values of a histogram to the range [0..1].
  @param buckets number of entries in arrays
  @param count   absolute histogram values
  @param normalized pointer to _allocated_ array with buckets entries
  @param type    type of normalization
  @return normalized values in 'normalized'
*/
void vvVolDesc::normalizeHistogram(int buckets, int* count, float* normalized, NormalizationType type)
{
  int max = 0;          // maximum counter value
  int i;

  vvDebugMsg::msg(2, "vvVolDesc::normalizeHistogram()");
  assert(count);

  // Find maximum counter value:
  for (i=0; i<buckets; ++i)   
    max = ts_max(max, count[i]);
  
  // Normalize counter values:
  for (i=0; i<buckets; ++i)   
  {
    if (count[i] < 1) normalized[i] = 0.0f;
    else if (type==VV_LOGARITHMIC)
    {
      normalized[i] = logf(float(count[i])) / logf(float(max));
    }
    else normalized[i] = float(count[i]) / float(max);
  }
}

//----------------------------------------------------------------------------
/** Generate voxel value histogram array.
  Counts:
  <ul>
    <li>number of data value occurrences for 8 or 16 bit per voxel data</li>
    <li>number of occurrences of each color component for RGB data</li>
    <li>number of density value occurrences for 32 bit per voxel data</li>
  </ul>
  @param frame   frame to generate histogram for (0 is first frame, -1 for all frames)
  @param buckets number of counters to use for histogram computation
  @param count   _allocated_ array of 'buckets' entries
  @return histogram values in 'count'
*/
void vvVolDesc::makeHistogram(int frame, int buckets, int* count)
{
  uchar* raw;           // raw voxel data
  int c, f, i;          // counters
  int numVoxels;        // number of voxels per frame
  int voxVal;           // voxel value
  int bucket;           // bucket
  int index;            // index into raw voxel data
  float valPerBucket;   // scalar data values per bucket

  vvDebugMsg::msg(2, "vvVolDesc::makeHistogram()");

  memset(count, 0, buckets * sizeof(int)); // initialize counter array

  numVoxels = getFrameVoxels();
  valPerBucket = float(getAlphaRange()) / float(buckets);
  for (f=0; f<frames; ++f)
  {
    if (frame>-1 && frame!=f) continue;     // only compute histogram for a specific frame
    raw = getRaw(f);
    for (i=index=0; i<numVoxels; ++i)  // count each voxel value
    {
      switch (bpv)
      {
        case 1: 
          voxVal = raw[index]; 
          bucket = int(float(voxVal) / valPerBucket);
          ++count[bucket];
          break;
        case 2: 
          voxVal = int(raw[index] << 8) | int(raw[index + 1]);
          bucket = int(float(voxVal) / valPerBucket);
          ++count[bucket];
          break;
        case 3:
          for (c=0; c<3; ++c)
          {
            voxVal = (int)raw[index + c];
            bucket = int(float(voxVal) / valPerBucket);
            ++count[bucket];
          }
          break;
        case 4:
          voxVal = raw[index + 3];
          bucket = int(float(voxVal) / valPerBucket);
          ++count[bucket];
          break;
        default: break;
      }
      index += bpv;
    }
  }
}

//----------------------------------------------------------------------------
/** Generates RGB texture values for the histogram.
 Texture values are returned as 3 bytes per texel, bottom to top, 
 ordered RGBRGBRGB...
 @param twidth,theight  width and height of texture [texels]
 @param alpha           true=generate alpha channel
 @param data            pointer to _pre-allocated_ memory space providing 
                        width * height * [3|4] bytes
*/
void vvVolDesc::makeHistogramTexture(int twidth, int theight, bool alpha, uchar* data)
{
  int    buckets;                   // resolution of histogram
  uchar  bg[4] = {  0,  0,  0,  0}; // background color (RGB)
  uchar  fg[4] = {255,255,255,255}; // foreground color (RGB)
  uchar* tex;                       // pointer to current texel
  float* hist;                      // histogram values (float)
  int*   count;                     // histogram values (integer)
  int    bpt;                       // bytes per texel: 3 for RGB, 4 for RGBA
  int    x, y;                      // current texel position
  int    barHeight;                 // histogram bar height [texels]
  int    histIndex;                 // index in histogram array
  int    c;                         // color component counter

  vvDebugMsg::msg(2, "vvVolDesc::makeHistogramTexture()");
  assert(data!=NULL);

  buckets = twidth;
  count = new int[buckets];
  hist = new float[buckets];
  makeHistogram(-1, buckets, count);
  normalizeHistogram(buckets, count, hist, VV_LOGARITHMIC);

  // Fill histogram texture with background values:
  if (alpha) bpt = 4;
  else       bpt = 3;
  tex = data;
  for (y=0; y<theight; ++y)
    for (x=0; x<twidth; ++x)
      for (c=0; c<bpt; ++c)
      {
        *tex = bg[c]; 
        ++tex;
      }

  // Draw histogram bars:
  for (x=0; x<twidth; ++x)
  {
    // Find histogram index for current x position in texture:
    histIndex = (int)(((float)x / (float)(twidth-1)) * (float)(buckets-1));
    ts_clamp(histIndex, 0, buckets-1);

    // Find height of histogram bar:
    barHeight = (int)(hist[histIndex] * (float)(theight-1));
    for (y=0; y<barHeight; ++y)
    {
      for (c=0; c<bpt; ++c)
        data[bpt * (x + y * twidth)] = fg[c];
    }
  }

  delete[] hist;
  delete[] count;
}

//----------------------------------------------------------------------------
/** Set file name.
  @param fn  file name including path (e.g. "c:\data\volumes\testvol.xvf")
*/
void vvVolDesc::setFilename(const char* fn)
{
  delete[] filename;
	if (fn==NULL) 
		filename = NULL;
  else
	{
		filename = new char[strlen(fn) + 1];
	  strcpy(filename, fn);
	}
  return;
}

//----------------------------------------------------------------------------
/** Get file name.
  @return pointer to file name. Don't delete this pointer.
*/
const char* vvVolDesc::getFilename()
{
  return filename;
}

//----------------------------------------------------------------------------
/** Set current frame.
  @param f  current frame
*/
void vvVolDesc::setCurrentFrame(int f)
{
  if (f>=0 && f<frames) currentFrame = f;
}

//----------------------------------------------------------------------------
/// Get current frame.
int vvVolDesc::getCurrentFrame()
{
  return currentFrame;
}

//----------------------------------------------------------------------------
/** Get range of alpha values.
 @return the number of distinct alpha values
*/
int vvVolDesc::getAlphaRange()
{
  return (bpv==2) ? 65536 : 256;
}

//----------------------------------------------------------------------------
/** Set volume size.
  @param w,h,s  number of voxels for width, height, slices
  @param b      byte per voxel
*/
void vvVolDesc::setVoxelSize(int w, int h, int s, int b)
{
  vox[0] = w;
  vox[1] = h;
  vox[2] = s;
  bpv = b;
}

//----------------------------------------------------------------------------
/** Set real world volume size.
  @param size world space volume size [mm] -> adjusts voxel distance values (dist[]).
              Requires correctly set voxel size (vox[])
*/
void vvVolDesc::setRealSize(vvVector3* size)
{
  for (int i=0; i<3; ++i)
    if (vox[i] > 1)   // prevent division by zero. don't change dist if vox=1
      dist[i] = size->e[i] / (vox[i] - 1);
}

//----------------------------------------------------------------------------
/** Get real world volume size.
  @param size return value: world space volume size [mm]
*/
void vvVolDesc::getRealSize(vvVector3* size)
{
  for (int i=0; i<3; ++i)
    size->e[i] = dist[i] * (vox[i] - 1);
}

//----------------------------------------------------------------------------
/** Converts the voxel format.
  The strategy depends on the number of bytes per voxel both in the source
  and in the destination volume:<PRE>
  Source  Destination  Strategy
  -----------------------------
  8 bit   16 bit       Shift source value left by 8 bit.
  8 bit   RGB          Make gray RGB value from 8 bit value.
  8 bit   RGBA         Make gray RGB value from 8 bit value and use 8 bit value for alpha.
  16 bit  8 bit        Shift source value right by 8 bit.
  16 bit  RGB          Shift source value right by 8 bit and continue as with 8 bit source.
  16 bit  RGBA         Shift source value right by 8 bit and continue as with 8 bit source.
  RGB     8 bit        Compute averate value of R,G,B and store as 8 bit value.
  RGB     16 bit       Compute averate value of R,G,B and store as 16 bit value.
  RGB     RGBA         Compute averate value of R,G,B and store as alpha.
  RGBA    8 bit        Use alpha as 8 bit value.
  RGBA    16 bit       Use alpha as 16 bit value.
  RGBA    RGB          Drop alpha.
  </PRE>
  bpv=1,2, or 3: first voxel byte is used<BR>
  bpv=4: last voxel byte (alpha) is used

  @param newBPV  new number of bytes per voxel
  @param verbose true = print progress info
*/
void vvVolDesc::convertVoxelFormat(int newBPV, bool verbose)
{
  uchar* newRaw;
  uchar* rd;
  uchar* src;
  uchar* dst;
  uchar srcVox[4];   // source voxel value
  uchar dstVox[4];   // destination voxel value
  float scalar, r, g, b;
  int newSliceSize;
  int x, y, z, f;
  uchar average;

  vvDebugMsg::msg(2, "vvVolDesc::convertVoxelFormat()");

  // Verify input parameters:
  if (bpv==newBPV) return;            // already done?
  if (newBPV<1 || newBPV>4) return;   // ignore invalid values

  // Check for simulation:
  if (simulate) 
  {
    bpv = newBPV;
    return;
  }

  newSliceSize = vox[0] * vox[1] * newBPV;
  if (verbose) vvToolshed::initProgress(vox[2] * frames);
  raw.first();
  for (f=0; f<frames; ++f)
  {
    rd = raw.getData();
    newRaw = new uchar[newSliceSize * vox[2]];
    src = rd;
    dst = newRaw;
    for (z=0; z<vox[2]; ++z)
    {
      for (y=0; y<vox[1]; ++y)
        for (x=0; x<vox[0]; ++x)
        {
          memcpy(srcVox, src, bpv);

          // Perform actual conversion:
          switch (bpv)  // switch by source voxel type
          {
            case 1: // 8 bit source
              // fall through
            case 2: // 16 bit source
              switch (newBPV)   // switch by destination voxel type
              {
                case 2: dstVox[1] = 0;
                  // fall through
                case 1: dstVox[0] = srcVox[0]; 
                  break;
                case 4: 
                  // fall through
                case 3:
                  if (bpv==1)
                    scalar = (float)srcVox[0] / 255.0f;
                  else 
                    scalar = ((float)srcVox[0] * 256.0f + (float)srcVox[1]) / 65535.0f;
                  tf.pins.getRGB(scalar, &r, &g, &b);
                  dstVox[0] = (uchar)(r * 255.0f);
                  dstVox[1] = (uchar)(g * 255.0f);
                  dstVox[2] = (uchar)(b * 255.0f);
                  if (newBPV==4)
                    dstVox[3] = srcVox[0];
                  break;
                default: break;
              }
              break;
            case 3: // RGB source
              average = (uchar)(((int)srcVox[0] + (int)srcVox[1] + (int)srcVox[2]) / 3);
              switch (newBPV)   // switch by destination voxel type
              {
                case 2: dstVox[1] = 0;
                  // fall through
                case 1: dstVox[0] = average; break;
                case 4: dstVox[0] = srcVox[0]; 
                        dstVox[1] = srcVox[1]; 
                        dstVox[2] = srcVox[2];
                        dstVox[3] = average;
                        break;
                default: break;
              }
              break;
            case 4: // RGBA source
              switch (newBPV)   // switch by destination voxel type
              {
                case 2: dstVox[1] = 0;
                  // fall through
                case 1: dstVox[0] = srcVox[3]; break;
                case 3: dstVox[0] = srcVox[0]; 
                        dstVox[1] = srcVox[1]; 
                        dstVox[2] = srcVox[2];
                        break;
                default: break;
              }
              break;
            default: break;
          }
          memcpy(dst, dstVox, newBPV);
          src += bpv;
          dst += newBPV;
        }
      if (verbose) vvToolshed::printProgress(z + vox[2] * f);
    }
    raw.remove();
    if (f==0) raw.insertBefore(newRaw, true);
    else raw.insertAfter(newRaw, true);
    raw.next();
  }      
  bpv = newBPV;
}

//----------------------------------------------------------------------------
/** Convert the number of modalities.
  @param newBPV  new number of modalities
  @param verbose true = print progress info
*/
void vvVolDesc::convertModalities(int newBPV, bool verbose)
{
  uchar* newRaw;
  uchar* rd;
  uchar* src;
  uchar* dst;
  uchar srcVox[4];   // source voxel value
  uchar dstVox[4];   // destination voxel value
  int newSliceSize;
  int x, y, z, f, i;

  vvDebugMsg::msg(2, "vvVolDesc::convertModalities()");

  // Verify input parameters:
  if (stype!=VV_MULTIMODAL)
  {
    cerr << "Storage type must be multi-modal to change modalities." << endl;
    return;
  }
  if (bpv==newBPV) return;            // already done?
  if (newBPV<1 || newBPV>4) return;   // ignore invalid values

  // Check for simulation:
  if (simulate) 
  {
    bpv = newBPV;
    return;
  }

  newSliceSize = vox[0] * vox[1] * newBPV;
  if (verbose) vvToolshed::initProgress(vox[2] * frames);
  raw.first();
  for (f=0; f<frames; ++f)
  {
    rd = raw.getData();
    newRaw = new uchar[newSliceSize * vox[2]];
    src = rd;
    dst = newRaw;
    for (z=0; z<vox[2]; ++z)
    {
      for (y=0; y<vox[1]; ++y)
      {
        for (x=0; x<vox[0]; ++x)
        {
          memcpy(srcVox, src, bpv);

          // Perform actual conversion:
          for (i=0; i<newBPV; ++i)
          {
            if (i < bpv) dstVox[i] = srcVox[i];
            else dstVox[i] = 0;
          }
          memcpy(dst, dstVox, newBPV);
          src += bpv;
          dst += newBPV;
        }
      }
      if (verbose) vvToolshed::printProgress(z + vox[2] * f);
    }
    raw.remove();
    if (f==0) raw.insertBefore(newRaw, true);
    else raw.insertAfter(newRaw, true);
    raw.next();
  }      
  bpv = newBPV;
}

//----------------------------------------------------------------------------
/** Bit-shifts all bytes of each voxel.
  @param bits number of bits to shift (>0 = right, <0 = left)
  @param verbose true = print progress info
*/
void vvVolDesc::bitShiftData(int bits, bool verbose)
{
  int  x, y, z, b, f;
  long pixel;
  long byte;
  int  shift;
  uchar* rd;
  int sliceSize;
  int offset;
  
  vvDebugMsg::msg(2, "vvVolDesc::bitShiftData()");
  if (bpv>sizeof(long)) return;   // shift only works up to sizeof(long) byte per pixel
  if (bits==0) return;
  if (simulate) return;   // this function does not affect the volume header values  

  sliceSize = getSliceSize();
  shift = ts_max(bits, -bits);    // find absolute value
  if (verbose) vvToolshed::initProgress(vox[2] * frames);
  raw.first();
  for (f=0; f<frames; ++f)
  {
    rd = raw.getData();
    raw.next();
    for (z=0; z<vox[2]; ++z)
    {
      for (y=0; y<vox[1]; ++y)
        for (x=0; x<vox[0]; ++x)
        {
          offset = x * bpv + y * vox[0] * bpv + z * sliceSize;
          pixel = 0;
          for (b=0; b<bpv; ++b)
          {
            byte = (int)rd[offset + b];
            byte = byte << ((bpv-b-1) * 8);
            pixel += byte;
          }
          if (bits>0)
            pixel = pixel >> shift;
          else
            pixel = pixel << shift;
          for (b=0; b<bpv; ++b)
          {
            byte = pixel >> ((bpv-b-1) * 8);
            rd[offset + b] = (uchar)(byte & 0xFF); 
          }
        }
      if (verbose) vvToolshed::printProgress(z + vox[2] * f);
    }
  }
}

//----------------------------------------------------------------------------
/** Invert all bits of each scalar voxel value.
*/
void vvVolDesc::invert()
{
  int  x, y, z, b, f;
  uchar* ptr;    // pointer to currently worked on byte
  uchar* rd;
  
  vvDebugMsg::msg(2, "vvVolDesc::invert()");
  
  if (simulate) return;  // this function does not affect the volume header values  

  raw.first();
  for (f=0; f<frames; ++f)
  {
    rd = raw.getData();
    raw.next();
    ptr = &rd[0];
    for (z=0; z<vox[2]; ++z)
      for (y=0; y<vox[1]; ++y)
        for (x=0; x<vox[0]; ++x)
          for (b=0; b<bpv; ++b)
          {
            *ptr = (uchar)(~(*ptr));
            ++ptr;
          }
  }
}

//----------------------------------------------------------------------------
/// Convert 24 bit RGB to 8 bit RGB332.
void vvVolDesc::convertRGB24toRGB8()
{
  int x, y, z, f, i;
  uchar* newRaw;
  int pixel;
  int color;
  int shift[3] = {3, 3, 2};  // 3 bit red, 3 bit green, 2 bit blue
  uchar* rd;
  int newSliceSize;
  int oldSliceSize;

  vvDebugMsg::msg(2, "vvVolDesc::convertRGB24toRGB8()");
  if (bpv!=3) return;   // cannot work on non-24bit-modes
  if (simulate)
  {
    bpv = 1;
    return;
  }

  oldSliceSize = getSliceSize();
  newSliceSize = vox[0] * vox[1];
  raw.first();
  for (f=0; f<frames; ++f)
  {
    rd = raw.getData();
    newRaw = new uchar[vox[0] * vox[1] * vox[2]];
    for (z=0; z<vox[2]; ++z)
      for (y=0; y<vox[1]; ++y)
        for (x=0; x<vox[0]; ++x)
        {
          pixel = 0;
          for (i=0; i<3; ++i)
          {
            color = rd[x * bpv + y * vox[0] * bpv + z * oldSliceSize + i];
            pixel &= 0xFF00;
            pixel += color;
            pixel <<= shift[i];
          }
          newRaw[x + y * vox[0] + z * newSliceSize] = (uchar)(pixel >> 8);
        }
    raw.remove();
    if (f==0) raw.insertBefore(newRaw, true);
    else raw.insertAfter(newRaw, true);
    raw.next();
  }
  bpv = 1;
}

//----------------------------------------------------------------------------
/** Flip voxel data along a coordinate axis.<PRE>
  Coordinate system:
      y    
      |_x
    z/     
  </PRE>
  Strategy:
  <UL>
  <LI>X Axis: Copy each line of voxels to a buffer and reversely copy voxels 
      back to the data set one-by-one.</LI>
  <LI>Y Axis: Copy each line of voxels to a buffer and exchange the
      line with the destination line.</LI>
  <LI>Z Axis: Copy each slice of voxels to a buffer and exchange it
      with the destination slice.</LI>
  </UL>
  @param axis axis along which to flip the volume data
*/
void vvVolDesc::flip(AxisType axis)
{
  uchar* rd;
  uchar* voxelData;   // temporary buffer for voxel data
  uchar* dst;         // destination pointer
  uchar* src;         // source pointer
  int lineSize;
  int sliceSize;
  int x, y, z, f;
  
  vvDebugMsg::msg(2, "vvVolDesc::flip()");
  if (simulate) return;   // this function does not affect the volume header values  

  lineSize = vox[0] * bpv;
  sliceSize = getSliceSize();
  if (axis==Z_AXIS)
    voxelData = new uchar[sliceSize];
  else
    voxelData = new uchar[lineSize];
  raw.first();
  for (f=0; f<frames; ++f)
  {
    rd = raw.getData();
    switch (axis)
    {
      case X_AXIS:
        dst = rd;
        for (z=0; z<vox[2]; ++z)
          for (y=0; y<vox[1]; ++y)
          {
            memcpy((void*)voxelData, (void*)dst, lineSize);
            src = voxelData + (vox[0]-1) * bpv;
            for (x=0; x<vox[0]; ++x)
            {
              memcpy(dst, src, bpv);
              dst += bpv;
              src -= bpv;
            }
          }
        break;
      case Y_AXIS:
        for (z=0; z<vox[2]; ++z)
          for (y=0; y<vox[1]/2; ++y)
          {
            src = rd + y * lineSize + z * sliceSize;
            dst = rd + (vox[1] - y - 1) * lineSize + z * sliceSize;
            memcpy((void*)voxelData, (void*)dst, lineSize);
            memcpy((void*)dst, (void*)src, lineSize);
            memcpy((void*)src, (void*)voxelData, lineSize);
          }
        break;
      case Z_AXIS:
        for (z=0; z<vox[2]/2; ++z)
        {
          dst = rd + z * sliceSize;
          src = rd + (vox[2]-z-1) * sliceSize;
          memcpy((void*)voxelData, (void*)dst, sliceSize);
          memcpy((void*)dst, (void*)src, sliceSize);
          memcpy((void*)src, (void*)voxelData, sliceSize);
        }
        break;
      default: break;
    }
    raw.next();
  }
  delete[] voxelData;
}

//----------------------------------------------------------------------------
/** Rotate voxel data about a coordinate axis.<PRE>
  Coordinate system:
      y    
      |_x
    z/     
  </PRE>
  @param axis axis about which to rotate the volume data
  @param dir  direction into which to rotate when looking at the origin 
              from the positive half of the chosen coordinate axis (-1=left, 1=right)
*/
void vvVolDesc::rotate(AxisType axis, int dir)
{
  uchar* rd;
  uchar* dst;         // destination pointer
  uchar* src;         // source pointer
  uchar* newRaw;      // new volume data
  int frameSize;
  int newWidth, newHeight, newSlices;   // dimensions of rotated volume
  int x, y, z, f;
  int xpos, ypos, zpos;
  
  vvDebugMsg::msg(2, "vvVolDesc::rotate()");
  if (dir!=-1 && dir!=1) return;        // validate direction

  // Compute the new volume size:
  switch (axis)
  {
    case X_AXIS:
      newWidth  = vox[0];
      newHeight = vox[2];
      newSlices = vox[1];
      break;
    case Y_AXIS:
      newWidth  = vox[2];
      newHeight = vox[1];
      newSlices = vox[0];
      break;
    case Z_AXIS:
      newWidth  = vox[1];
      newHeight = vox[0];
      newSlices = vox[2];
      break;
    default:  // no change
      newWidth  = vox[0];
      newHeight = vox[1];
      newSlices = vox[2];
      break;
  }

  if (simulate)
  {
    vox[0] = newWidth;
    vox[1] = newHeight;
    vox[2] = newSlices;
    return;
  }

  frameSize = getFrameSize();
  raw.first();
  for (f=0; f<frames; ++f)
  {
    rd = raw.getData();
    newRaw = new uchar[frameSize];
    src = rd;
    switch (axis)
    {
      case X_AXIS:
        for (y=0; y<newHeight; ++y)
        {
          if (dir>0) ypos = y;
          else       ypos = newHeight - 1 - y;
          for (z=0; z<newSlices; ++z)
          {
            if (dir>0) zpos = newSlices - 1 - z;
            else       zpos = z;
            for (x=0; x<newWidth; ++x)
            { 
              dst = newRaw + bpv * (x + ypos * newWidth + zpos * newWidth * newHeight);
              memcpy((void*)dst, (void*)src, bpv);
              src += bpv;
            }
          }
        }
        break;
      case Y_AXIS:
        for (x=0; x<newWidth; ++x)
        {
          if (dir>0) xpos = x;
          else       xpos = newWidth - 1 - x;
          for (y=0; y<newHeight; ++y)
            for (z=0; z<newSlices; ++z)
            { 
              if (dir>0) zpos = newSlices - 1 - z;
              else       zpos = z;
              dst = newRaw + bpv * (xpos + y * newWidth + zpos * newWidth * newHeight);
              memcpy((void*)dst, (void*)src, bpv);
              src += bpv;
            }
        }
        break;
      case Z_AXIS:
        for (z=0; z<newSlices; ++z)
          for (x=0; x<newWidth; ++x)
          {
            if (dir>0) xpos = newWidth - 1 - x;
            else       xpos = x;
            for (y=0; y<newHeight; ++y)
            { 
              if (dir>0) ypos = y;
              else       ypos = newHeight - 1 - y;
              dst = newRaw + bpv * (xpos + ypos * newWidth + z * newWidth * newHeight);
              memcpy((void*)dst, (void*)src, bpv);
              src += bpv;
            }
          }
        break;
      default: break;
    }
    raw.remove();
    if (f==0) raw.insertBefore(newRaw, true);
    else raw.insertAfter(newRaw, true);
    raw.next();
  }
  vox[0] = newWidth;
  vox[1] = newHeight;
  vox[2] = newSlices;
}

//----------------------------------------------------------------------------
/** Convert 24 bit RGB planar (RRRR..., GGGG..., BBBB...)
   to 24 bit RGB interleaved (RGB, RGB, RGB, ...).
   Planar format must be repeated in each volume frame.
*/
void vvVolDesc::convertRGBPlanarToRGBInterleaved()
{
  uchar* raw;
  uchar* tmpData;
  int i, c, f;
  int voxels;
  int frameSize;
  
  vvDebugMsg::msg(2, "vvVolDesc::convertRGBPlanarToRGBInterleaved()");
  if (bpv != 3) return;   // this routine only works on RGB volumes
  if (simulate) return; 
  
  frameSize = getFrameSize();
  voxels = getFrameVoxels();
  tmpData = new uchar[frameSize];
  for (f=0; f<frames; ++f)
  {
    raw = getRaw(f);
    for (i=0; i<voxels; ++i)
      for (c=0; c<3; ++c)
        tmpData[i * 3 + c] = raw[c * voxels + i];  
    memcpy(raw, tmpData, frameSize);
  }
  delete[] tmpData;
}

//----------------------------------------------------------------------------
/** Toggle endianness of bytes in voxels.
  <UL>
    <LI>bpv=1: nothing to be done</LI>
    <LI>bpv=2: high and low byte are swapped</LI>
    <LI>bpv=3: original byte order: RGB,  order after toggle: BGR</LI>
    <LI>bpv=4: original byte order: RGBA, order after toggle: ABGR</LI>
  </UL>
*/
void vvVolDesc::toggleEndianness()
{
  uchar* rd;
  int    x, y, z, f;
  int    sliceSize;
  int    rowOffset, sliceOffset, voxelOffset;
  uchar  buffer;
  
  vvDebugMsg::msg(2, "vvVolDesc::toggleEndianness()");
  if (bpv<2) return;   // for 1 bpv nothing is to be done
  if (simulate) return; 
  
  sliceSize = getSliceSize();
  raw.first();
  for (f=0; f<frames; ++f)
  {
    rd = raw.getData();
    raw.next();
    for (z=0; z<vox[2]; ++z)
    {
      sliceOffset = z * sliceSize;
      for (y=0; y<vox[1]; ++y)
      {
        rowOffset = sliceOffset + y * vox[0] * bpv;
        for (x=0; x<vox[0]; ++x)
        {
          voxelOffset = x * bpv + rowOffset;

          // Swap first and last byte of voxel:
          buffer = rd[voxelOffset];
          rd[voxelOffset] = rd[voxelOffset + bpv - 1];
          rd[voxelOffset + bpv - 1] = buffer;

          // For 32 bit voxels also swap middle bytes:
          if (bpv==4)
          {
            buffer = rd[voxelOffset + 1];
            rd[voxelOffset + 1] = rd[voxelOffset + 2];
            rd[voxelOffset + 2] = buffer;
          }
        }
      }
    } 
  }
}

//----------------------------------------------------------------------------
/** Toggle sign of data values.
  <UL>
    <LI>bpv=1: invert most significant bit</LI>
    <LI>bpv=2: invert most significant bit</LI>
    <LI>bpv=3: invert most significant bits of each byte</LI>
    <LI>bpv=4: invert most significant bits of each byte</LI>
  </UL>
*/
void vvVolDesc::toggleSign()
{
  uchar* rd;
  int    sliceSize;
  int    f, i;
  
  vvDebugMsg::msg(2, "vvVolDesc::toggleSign()");
  if (simulate) return; 
  
  sliceSize = getSliceSize();
  raw.first();
  for (f=0; f<frames; ++f)
  {
    rd = raw.getData();
    raw.next();
    for (i=0; i<sliceSize; i+=2)
    {
      rd[i] ^= 0x80;      // bitwise exclusive or with 10000000b to negate MSB
      if (bpv!=2)   // don't negate low byte of 16 bit data but every other type
      {
        rd[i+1] ^= 0x80;      // bitwise exclusive or with 10000000b to negate MSB
      }
    } 
  }
}

//----------------------------------------------------------------------------
/** Returns true if a specific voxel byte is nonzero in any voxel of the volume.
  @param b  byte index (0=first byte) [0..3]
*/
bool vvVolDesc::isByteUsed(int b)
{
  uchar* rd;
  int    x, y, z, f;
  int    sliceSize;
  int    rowOffset, sliceOffset, voxelOffset;
  
  vvDebugMsg::msg(2, "vvVolDesc::isByteUsed()");
  if (b<0 || b>=bpv) return false;    // check for valid byte index
  
  sliceSize = getSliceSize();
  raw.first();
  for (f=0; f<frames; ++f)
  {
    rd = raw.getData();
    raw.next();
    for (z=0; z<vox[2]; ++z)
    {
      sliceOffset = z * sliceSize;
      for (y=0; y<vox[1]; ++y)
      {
        rowOffset = sliceOffset + y * vox[0] * bpv;
        for (x=0; x<vox[0]; ++x)
        {
          voxelOffset = x * bpv + rowOffset;
          if (rd[voxelOffset + b] != 0) return true;
        }
      }
    }
  }
  return false;
}

//----------------------------------------------------------------------------
/** Crop each volume of the animation to a sub-volume.
  @param x,y,z  coordinates of top-left-front corner of sub-volume
  @param w,h,s  width, height, and number of slices of sub-volume
*/
void vvVolDesc::crop(int x, int y, int z, int w, int h, int s)
{
  int j, i, f;
  uchar* newRaw;
  uchar* rd;
  int xmin, xmax, ymin, ymax, zmin, zmax;
  int newWidth, newHeight, newSlices;
  int newSliceSize;
  int oldSliceSize;
  uchar *src, *dst;

  vvDebugMsg::msg(2, "vvVolDesc::crop()");

  // Find minimum and maximum values for crop:
  xmin = ts_max(0, ts_min(x, x + w - 1));  
  ymin = ts_max(0, ts_min(y, y + h - 1));  
  zmin = ts_max(0, ts_min(z, z + s - 1));  
  xmax = ts_min(vox[0] -1, ts_max(x, x + w - 1));  
  ymax = ts_min(vox[1]-1, ts_max(y, y + h - 1));  
  zmax = ts_min(vox[2]-1, ts_max(z, z + s - 1));  

  // Set new volume dimensions:
  newWidth  = xmax - xmin + 1;
  newHeight = ymax - ymin + 1;
  newSlices = zmax - zmin + 1;

  // Check for simulation:
  if (simulate)
  {
    vox[0] = newWidth;
    vox[1] = newHeight;
    vox[2] = newSlices;
    return;
  }

  // Now cropping can be done:
  oldSliceSize = getSliceSize();
  newSliceSize = newWidth * newHeight * bpv;
  raw.first();
  for (f=0; f<frames; ++f)
  {
    rd = raw.getData();
    newRaw = new uchar[newSliceSize * newSlices];
    for (j=0; j<newSlices; ++j)
      for (i=0; i<newHeight; ++i)
      {
        src = rd + (j + zmin) * oldSliceSize + 
              (i + ymin) * vox[0] * bpv + xmin * bpv;
        dst = newRaw + j * newSliceSize + i * newWidth * bpv;
        memcpy(dst, src, newWidth * bpv);
      }
    raw.remove();
    if (f==0) raw.insertBefore(newRaw, true);
    else raw.insertAfter(newRaw, true);
    raw.next();
  }
  vox[0] = newWidth;
  vox[1] = newHeight;
  vox[2] = newSlices;
}

//----------------------------------------------------------------------------
/** Remove time steps at start and end of volume animation.
  @param start first time step to keep in the sequence [0..timesteps-1]
  @param steps number of steps to keep
*/
void vvVolDesc::cropTimesteps(int start, int steps)
{
  int i;

  if (simulate)
  {
    frames = steps;
    return;
  }

  raw.first();
  
  // Remove steps before the desired range:
  for (i=0; i<start; ++i)
    raw.remove();
  
  // Remove steps after the desired range:
  raw.last();
  for (i=0; i<frames-start-steps; ++i)
    raw.remove();

  frames = raw.count();
}

//----------------------------------------------------------------------------
/** Resize each volume of the animation. The real voxel size parameters
  are adjusted accordingly.
  @param w,h,s   new width, height, and number of slices
  @param ipt     interpolation type to use for resampling
  @param verbose true = verbose mode
*/
void vvVolDesc::resize(int w, int h, int s, InterpolationType ipt, bool verbose)
{
  uchar* newRaw;            // pointer to new volume data
  uchar* rd;
  int newSliceSize, newFrameSize;
  int oldSliceVoxels;
  uchar *src, *dst;
  int ix, iy, iz;           // integer source voxel coordinates
  float fx, fy, fz;         // floating point source voxel coordinates
  uchar interpolated[4];    // interpolated voxel values
  int f, x, y, z;

  vvDebugMsg::msg(2, "vvVolDesc::resize()");

  // Validate resize parameters:
  if (w<=0 || h<=0 || s<=0) return;
  if (w==vox[0] && h==vox[1] && s==vox[2]) return; // already done

  // Check for simulation:
  if (simulate)
  {
    vox[0] = w;
    vox[1] = h;
    vox[2] = s;
    return;
  }

  // Now resizing can be done:
  oldSliceVoxels = getSliceVoxels();
  newSliceSize = w * h * bpv;
  newFrameSize = newSliceSize * s;
  if (verbose) vvToolshed::initProgress(s * frames);
  raw.first();
  for (f=0; f<frames; ++f)
  {
    rd = raw.getData();
    newRaw = new uchar[newFrameSize];
    dst = newRaw;

    // Traverse destination data:
    for (z=0; z<s; ++z)
    {
      for (y=0; y<h; ++y)
        for (x=0; x<w; ++x)
        {
          // Compute source coordinates of current destination voxel:
          if (ipt==TRILINEAR)   // trilinear interpolation
          {
            // Compute source coordinates of current destination voxel:
            fx = (float)x / (float)(w-1) * (float)(vox[0]-1);
            fy = (float)y / (float)(h-1) * (float)(vox[1]-1);
            fz = (float)z / (float)(s-1) * (float)(vox[2]-1);
            trilinearInterpolation(f, fx, fy, fz, interpolated);

            // Copy interpolated voxel data to destination voxel:
            memcpy(dst, interpolated, bpv);
          }
          else    // nearest neighbor interpolation
          {
            // Compute source coordinates of current destination voxel:
            if (w>1) ix = x * (vox[0]-1)  / (w-1);
            else     ix = 0;
            if (h>1) iy = y * (vox[1]-1) / (h-1);
            else     iy = 0;
            if (s>1) iz = z * (vox[2]-1) / (s-1);
            else     iz = 0;
            ix = ts_clamp(ix, 0, vox[0]-1);
            iy = ts_clamp(iy, 0, vox[1]-1);
            iz = ts_clamp(iz, 0, vox[2]-1);

            // Copy source voxel data to destination voxel:
            src = rd + bpv * (ix + iy * vox[0] + iz * oldSliceVoxels);
            memcpy(dst, src, bpv);
          }
          dst += bpv;
        }
      if (verbose) vvToolshed::printProgress(z + s * f);
    }
    raw.remove();
    if (f==0) raw.insertBefore(newRaw, true);
    else raw.insertAfter(newRaw, true);
    raw.next();
  }
  // Adjust voxel size:
  dist[0] *= float(vox[0]) / float(w);
  dist[1] *= float(vox[1]) / float(h);
  dist[2] *= float(vox[2]) / float(s);

  // Use new size:
  vox[0] = w;
  vox[1] = h;
  vox[2] = s;
}

//----------------------------------------------------------------------------
/** Shift each volume of the animation by a number of voxels.
  Rotary boundary conditions are applied.
  Strategy: The volume is shifted slicewise so that only one
  volume slice has to be duplicated in memory.
  @param sx,sy,sz  shift amounts. Positive values shift into the
                   positive axis direction.
*/
void vvVolDesc::shift(int sx, int sy, int sz)
{
  uchar* rd;
  uchar* newRaw;    // shifted volume data
  uchar* src;
  uchar* dst;
  int lineSize, sliceSize, frameSize; 
  int sval[3];      // shift amount
  int f, x, y, z, i;

  vvDebugMsg::msg(2, "vvVolDesc::shift()");

  if (simulate) return; // this function does not affect the volume header

  // Consider rotary boundary conditions and make shift values positive:
  if (sx==0 && sy==0 && sz==0) return;

  sval[0] =  sx % vox[0];
  sval[1] = -sy % vox[1];
  sval[2] = -sz % vox[2];
  if (sval[0]<0) sval[0] += vox[0];
  if (sval[1]<0) sval[1] += vox[1];
  if (sval[2]<0) sval[2] += vox[2];

  // Now shifting starts:
  lineSize  = vox[0] * bpv;
  sliceSize = getSliceSize();
  frameSize = getFrameSize();
  raw.first();
  for (f=0; f<frames; ++f)
  {
    rd = raw.getData();

    for (i=0; i<3; ++i)
    {
      if (sval[i] > 0)
      {
        newRaw = new uchar[frameSize];
        switch (i)
        {
          case 0:   // x shift
            src = rd;
            for (z=0; z<vox[2]; ++z)
              for (y=0; y<vox[1]; ++y)
                for (x=0; x<vox[0]; ++x)
                {
                  dst = newRaw + z * sliceSize + y * lineSize + ((x + sval[0]) % vox[0]) * bpv;
                  memcpy(dst, src, bpv);
                  src += bpv;
                }
            break;
          case 1:   // y shift
            src = rd;
            for (z=0; z<vox[2]; ++z)
              for (y=0; y<vox[1]; ++y)
              {
                dst = newRaw + z * sliceSize + ((y + sval[1]) % vox[1]) * lineSize;
                memcpy(dst, src, lineSize);
                src += lineSize;
              }
            break;
          default:
          case 2:   // z shift
            src = rd;
            for (z=0; z<vox[2]; ++z)
            {
              dst = newRaw + ((z + sval[2]) % vox[2]) * sliceSize;
              memcpy(dst, src, sliceSize);
              src += sliceSize;
            }
            break;
        }
        raw.remove();
        if (f==0) raw.insertBefore(newRaw, true);
        else raw.insertAfter(newRaw, true);
      }
    }   
    raw.next();
  }
}

//----------------------------------------------------------------------------
/** Convert COVISE volume data to Virvo format.
 The difference is that in COVISE the innermost and the outermost
 voxel loops are inverted.
*/
void vvVolDesc::convertCoviseToVirvo()
{
  uchar* raw;
  uchar* tmpData;
  int z, y, x, f;
  int frameSize;
	int srcIndex; 	// index into COVISE volume array
	uchar* ptr;
  
  vvDebugMsg::msg(2, "vvVolDesc::convertCoviseToVirvo()");
  if (simulate) return;

  frameSize = getFrameSize();
  tmpData = new uchar[frameSize];
  for (f=0; f<frames; ++f)
  {
    raw = getRaw(f);
		ptr = tmpData;
		for (z=0; z<vox[2]; ++z)
			for (y=0; y<vox[1]; ++y)
				for (x=0; x<vox[0]; ++x)
				{
					srcIndex = bpv * ((vox[2]-z-1)+ (vox[1]-y-1) * vox[2] + x * vox[1] * vox[2]);
					memcpy(ptr,	raw + srcIndex, bpv);
					ptr += bpv; 	// skip to next voxel
				}
    memcpy(raw, tmpData, frameSize);
	}
  delete[] tmpData;
  
#ifdef __linux
  if (bpv==4) toggleEndianness();   // RGBA data are transferred as packed colors, which are integers
#endif
}

//----------------------------------------------------------------------------
/** Convert Virvo volume data to COVISE format.
 The difference is that the innermost and the outermost
 voxel loops are inverted.
*/
void vvVolDesc::convertVirvoToCovise() 
{
  uchar* raw;
  uchar* tmpData;
  uchar* ptr;
  int    z, y, x, f;
  int    frameSize;
  int    dstIndex; 	// index into COVISE volume array
  
  vvDebugMsg::msg(2, "vvVolDesc::convertVirvoToCovise()");
  if (simulate) return;

  frameSize = getFrameSize();
  tmpData = new uchar[frameSize];
  for (f=0; f<frames; ++f)
  {
    raw = getRaw(f);
  	ptr = raw;
  	for (z=0; z<vox[2]; ++z)
	  	for (y=0; y<vox[1]; ++y)
		  	for (x=0; x<vox[0]; ++x)
			  {
  				dstIndex = bpv * (x * vox[1] * vox[2] + (vox[1]-y-1) * vox[2] + (vox[2]-z-1));
	  			memcpy(tmpData + dstIndex, ptr, bpv);
					ptr += bpv; 	// skip to next voxel
			  }
    memcpy(raw, tmpData, frameSize);
  }
  delete[] tmpData;
}

//----------------------------------------------------------------------------
/** Convert Virvo volume data to OpenGL format.
 The OpenGL format is as follows: counting starts at the backmost slice bottom left,
 it continues to the right, then up and then to the front.
*/
void vvVolDesc::convertVirvoToOpenGL()
{
  uchar* raw;
  uchar* tmpData;
	uchar* ptr;
  int    z, y, x, f;
  int    frameSize;
	int    dstIndex; 	// index into OpenGL volume array
  
  vvDebugMsg::msg(2, "vvVolDesc::convertVirvoToOpenGL()");
  if (simulate) return;

  frameSize = getFrameSize();
  tmpData = new uchar[frameSize];
  for (f=0; f<frames; ++f)
  {
    raw = getRaw(f);
  	ptr = raw;
		for (z=0; z<vox[2]; ++z)
			for (y=0; y<vox[1]; ++y)
				for (x=0; x<vox[0]; ++x)
				{
					dstIndex = bpv * (x + (vox[1] - y - 1) * vox[0] + (vox[2] - z - 1) * vox[0] * vox[1]);
					memcpy(tmpData + dstIndex, ptr, bpv);
					ptr += bpv;
				}
    memcpy(raw, tmpData, frameSize);
	}
  delete[] tmpData;
}

//----------------------------------------------------------------------------
/** Convert OpenGL volume data to Virvo format.
 The OpenGL format is as follows: counting starts at the backmost slice bottom left,
 it continues to the right, then up and then to the front.
*/
void vvVolDesc::convertOpenGLToVirvo()
{
  uchar* raw;
  uchar* tmpData;
	uchar* ptr;
  int    z, y, x, f;
  int    frameSize;
	int    dstIndex; 	// index into Virvo volume array
  
  vvDebugMsg::msg(2, "vvVolDesc::convertOpenGLToVirvo()");
  if (simulate) return;

  frameSize = getFrameSize();
  tmpData = new uchar[frameSize];
  for (f=0; f<frames; ++f)
  {
    raw = getRaw(f);
  	ptr = raw;
		for (z=0; z<vox[2]; ++z)
			for (y=0; y<vox[1]; ++y)
				for (x=0; x<vox[0]; ++x)
				{
					dstIndex = bpv * (x + (vox[1] - y - 1) * vox[0] + (vox[2] - z - 1)  * vox[0] * vox[1]);
					memcpy(tmpData + dstIndex, ptr, bpv);
					ptr += bpv;
				}
    memcpy(raw, tmpData, frameSize);
	}
  delete[] tmpData;
}

//----------------------------------------------------------------------------
/// Return the size of the icon (size = width = height).
int vvVolDesc::getIconSize()
{
  return iconSize;
}

//----------------------------------------------------------------------------
/** Return the icon's RGB data (uchar array of iconSize*iconSize*3 bytes).
 The caller must not delete this array!
*/
uchar* vvVolDesc::getIconData()
{
  return iconData;
}

//----------------------------------------------------------------------------
/** Creates an icon from passed RGB data.
 The caller must delete the passed rgb array.
 @param s    icon size in pixels (size = width = height), if 0 icon will be deleted
 @param rgb  icon image data array (size * size * 3 bytes expected)
*/
void vvVolDesc::makeIcon(int s, const uchar* rgb)
{
  vvDebugMsg::msg(2, "vvVolDesc::makeIcon()");

  if (s<=0) 
  {
    iconSize = 0;
    if (iconData!=NULL) delete[] iconData;
    iconData = NULL;
  }
  else
  {
    iconSize = s;
    iconData = new uchar[iconSize * iconSize * 3];
    memcpy(iconData, rgb, iconSize * iconSize * 3);
  }  
}


//----------------------------------------------------------------------------
/** Convert flat data to a sphere.
  The sphere coordinate equations are taken from Bronstein page 154f.
  The routine traverses the destination voxels, reversely computes
  their locations on the planar data set, and copies these voxel
  values. The source volume is mapped on a sphere which can be 
  imagined to be located 'behind' the source volume: Up remains
  up, the 'seam' is located at the front of the sphere. Lower
  z values become located a more remote regions from the sphere 
  center.<P>
  The used coordinate system is:<PRE>
        z
       /
      /___x
     | 
     |
     y
  </PRE>
  R is the distance from the current voxel to the volume center.<BR>
  Phi is the angle in the x/z axis, starting at (0|0|-1), going
  past (1|0|0) and (0|0|1) and back to (0|0|-1).<BR>
  Theta is the angle from the position vector of the current
  voxel to the vector (0|-1|0).
  @param outer   outer sphere diameter [voxels]
  @param inner   inner sphere diameter [voxels]
  @param ipt     interpolation type
  @param verbose true = verbose mode
*/
void vvVolDesc::makeSphere(int outer, int inner, InterpolationType ipt, bool verbose)
{
  uchar* rd;              // raw data of current source frame
  vvVector3 center;       // sphere center position [voxel space]
  vvVector3 v;            // currently processed voxel coordinates [voxel space]
  uchar* newRaw;          // raw data of current destination frame
  float dist;             // distance from voxel to sphere center
  float phi;              // angle in x/z plane
  float theta;            // angle to y axis
  uchar *src, *dst;       // source and destination volume data
  float radius;           // outer radius of sphere [voxels]
  float core;             // core radius [voxels]
  int newFrameSize;       // sphere volume frame size [voxels]
  int sliceVoxels;        // number of voxels per slice in source volume
  int f, x, y, z;         // loop counters
  float sx, sy, sz;       // coordinates in source volume
  float ringSize;         // precomputed ring size
  uchar interpolated[4];  // interpolated voxel values

  vvDebugMsg::msg(2, "vvVolDesc::makeSphere()");
  if (outer<1 || inner<0) return;
  if (simulate)
  {
    vox[0] = vox[1] = vox[2] = outer;
    return;
  }

  newFrameSize = outer * outer * outer * bpv;
  if (outer>1)
    radius = (float)(outer-1) / 2.0f;
  else
    radius = 0.5f;
  if (inner>1)
    core = (float)(inner-1) / 2.0f;
  else
    core = 0.0f;
  ringSize = radius - core;
  if (ringSize<1.0f) ringSize = 1.0f; // prevent division by zero later on
  center.set(radius, radius, radius);
  sliceVoxels = vox[0] * vox[1];
  if (verbose) vvToolshed::initProgress(outer * frames);
  raw.first();
  for (f=0; f<frames; ++f)
  {
    rd = raw.getData();
    newRaw = new uchar[newFrameSize];
    dst = newRaw;

    // Traverse destination data:
    for (z=0; z<outer; ++z)
    {
      for (y=0; y<outer; ++y)
        for (x=0; x<outer; ++x)
        {
          // Compute sphere coordinates of current destination voxel:
          v.set((float)x, (float)y, (float)z);
          v.sub(&center);
          v[1] = -v[1];   // adapt to vvVecmath coordinate system
          v[2] = -v[2];   // adapt to vvVecmath coordinate system
          v.getSpherical(&dist, &phi, &theta);

          // Map sphere coordinates to planar coordinates in source volume:
          if (dist<=radius && dist>=core)   // is voxel within source volume?
          {
            // Compute source coordinates of current destination voxel:
            sx = phi / (2.0f * VV_PI) * (float)vox[0];
            sy = theta / VV_PI * (float)vox[1];
            sz = (float)vox[2] - 1.0f - ((dist-core) / ringSize * (float)vox[2]);
            if (ipt==TRILINEAR)  // trilinear interpolation
            {
              trilinearInterpolation(f, sx, sy, sz, interpolated);
              memcpy(dst, interpolated, bpv);
            }
            else    // nearest neighbor
            {
              sx = ts_clamp(sx, 0.0f, (float)(vox[0]-1));
              sy = ts_clamp(sy, 0.0f, (float)(vox[1]-1));
              sz = ts_clamp(sz, 0.0f, (float)(vox[2]-1));
              src = rd + bpv * ((int)sx + (int)sy * vox[0] + (int)sz * sliceVoxels);
              memcpy(dst, src, bpv);
            }
          }
          else
            memset(dst, 0, bpv);    // outside of sphere

          dst += bpv;
        }
      if (verbose) vvToolshed::printProgress(z + outer * f);
    }
    raw.remove();
    if (f==0) raw.insertBefore(newRaw, true);
    else raw.insertAfter(newRaw, true);
    raw.next();
  }
  vox[0] = vox[1] = vox[2] = outer;
}

//----------------------------------------------------------------------------
/** Display one line of volume information data.
  @param desc   description text, NULL for none
*/
void vvVolDesc::printInfoLine(const char* desc)
{
  if (desc!=NULL)
    cerr << desc << " ";

  cerr << frames;
  if (frames!=1) cerr << " frames, ";
  else cerr << " frame, ";
  cerr << vox[0] << " x " << vox[1] << " x " << vox[2];
  if (frames > 1) cerr << " voxels per frame, ";
  else cerr << " voxels, ";
  cerr << bpv;
  if (bpv > 1) cerr << " bytes";
  else cerr << " byte";
  cerr << " per voxel" << endl;
}

//----------------------------------------------------------------------------
/// Display verbose volume information from file header.
void vvVolDesc::printVolumeInfo()
{
  cerr << "File name:                         " << filename << endl;
  cerr << "Volume size [voxels]:              " << vox[0] << " x " << vox[1] << " x " << vox[2] << endl;
  cerr << "Number of frames:                  " << frames << endl;
  cerr << "Bytes per voxel:                   " << bpv << endl;
  cerr << "Storage type:                      " << int(stype) << endl;
  cerr << "Voxels per frame:                  " << getFrameVoxels() << endl;
  cerr << "Bytes per frame:                   " << getFrameSize() << endl;
  cerr << "Sample distances:                  " << setprecision(3) << dist[0] << " x " << dist[1] << " x " << dist[2] << endl;
  cerr << "Time step duration [s]:            " << setprecision(3) << dt << endl;
  cerr << "Physical data range:               " << realMin << " to " << realMax << endl;
  cerr << "Object location [mm]:              " << position[0] << ", " << position[1] << ", " << position[2] << endl;
  cerr << "Icon stored:                       " << ((iconSize > 0) ? "yes" : "no") << endl;
}

//----------------------------------------------------------------------------
/// Display statistics about the volume data.
void vvVolDesc::printStatistics()
{
  int scalarMin, scalarMax;

  findMinMax(&scalarMin, &scalarMax);
  cerr << "Scalar value range:                " << scalarMin << " to " << scalarMax << endl;
  cerr << "Number of different data values:   " << findNumUsed() << endl;
  cerr << "Zero voxels in first frame:        " << setprecision(4) << 
    (100.0f * findNumValue(0,0) / getFrameVoxels()) << " %" << endl;
  cerr << "Transparent voxels in first frame: " << setprecision(4) << 
    (100.0f * findNumTransparent(0) / getFrameVoxels()) << " %" << endl;
}

//----------------------------------------------------------------------------
/** Display histogram as text. Omit non-occuring values.
  @param frame   frame to compute histogram for (-1 for all frames)
  @param buckets number of buckets to use for values
*/
void vvVolDesc::printHistogram(int frame)
{
  int* hist;
  int i;
  int buckets;
  
  buckets = getAlphaRange();
  hist = new int[buckets];
  makeHistogram(frame, buckets, hist);
  for (i=0; i<buckets; ++i)
  {
    if (hist[i] > 0) cerr << i << ": " << hist[i] << endl;
  }
  delete[] hist;
}

//----------------------------------------------------------------------------
/** Display the voxel data as hex numbers.
  @param frame number of frame in dataset
  @param slice number of slice to display
  @param width number of voxels to display per row (starting left, 0 for all voxels)
  @param height number of voxels to display per column (starting at top, 0 for all voxels))
*/
void vvVolDesc::printVoxelData(int frame, int slice, int width, int height)
{
  uchar* raw;   // pointer to raw volume data
  int x,y;      // voxel counters
  int nx,ny;    // number of voxels to display per row/column
  int val;
  
  raw = getRaw(frame);
  raw += slice * getSliceSize();  // move pointer to beginning of slice to display
  if (width<=0) nx = vox[0];
  else nx = ts_min(width, vox[0]);
  if (height<=0) ny = vox[1];
  else ny = ts_min(height, vox[1]);
  cerr << "Voxel data of frame " << frame << ":" << endl;
  for (y=0; y<ny; ++y)
  {
    cerr << "Row " << y << ": ";
    for (x=0; x<nx; ++x)
    {
      switch (bpv)
      {
        default:
        case 1:
          val = int(raw[bpv * (x + y * vox[0])]);
          cerr << val << " ";
          break;
        case 2:
          val = 256 * int(raw[bpv * (x + y * vox[0])]) + 
                int(raw[1 + bpv * (x + y * vox[0])]);
          cerr << val << " ";
          break;
        case 3:
          cerr << "{" << int(raw[bpv * (x + y * vox[0])]) << "," <<
                         int(raw[1 + bpv * (x + y * vox[0])]) << "," <<
                         int(raw[2 + bpv * (x + y * vox[0])]) << "} ";
          break;
        case 4:
          cerr << "{" << int(raw[bpv * (x + y * vox[0])]) << "," <<
                         int(raw[1 + bpv * (x + y * vox[0])]) << "," <<
                         int(raw[2 + bpv * (x + y * vox[0])]) << "," <<
                         int(raw[3 + bpv * (x + y * vox[0])]) << "} ";
          break;
      }
    }
    cerr << endl;
  }
}

//----------------------------------------------------------------------------
/** Perform a trilinear interpolation.
  Neighboring vertex indices:<PRE>
      4____ 7
     /___ /| 
   0|   3| | 
    | 5  | /6
    |/___|/  
    1    2
  </PRE>
  @param f      volume animation frame
  @param x,y,z  coordinates for which the interpolation must be performed
  @param result interpolated voxel data. This must be a pointer to 
                an _allocated_ memory space of bpv bytes.
*/
void vvVolDesc::trilinearInterpolation(int f, float x, float y, float z, uchar* result)
{
  uchar* rd;
  uchar* neighbor[8];   // pointers to neighboring voxels
  float val[8];         // neighboring voxel values
  int tfl[3];           // coordinates of neighbor 0 (top-front-left)
  float dist[3];        // distance to neighbor 0
  int lineSize, sliceSize;
  int i, j, passes;
  int interpolated;   // interpolated value

  vvDebugMsg::msg(3, "vvVolDesc::trilinearInterpolation()");

  if (vox[0]<2 || vox[1]<2 || vox[2]<2) return;

  // Check for valid frame index:
  if (f<0 || f>=frames) return; 
  
  // Constrain voxel position to the valid region:
  x = ts_clamp(x, 0.0f, (float)(vox[0]-1));
  y = ts_clamp(y, 0.0f, (float)(vox[1]-1));
  z = ts_clamp(z, 0.0f, (float)(vox[2]-1));

  // Compute coordinates of neighbor 0:
  tfl[0] = (int)x;
  tfl[1] = (int)y;
  tfl[2] = (int)z;

  // Compute distance to neighbor 0:
  if (tfl[0] < (vox[0]-1))
    dist[0] = x - (float)tfl[0];
  else   // border values need special treatment
  {
    --tfl[0];
    dist[0] = 1.0f;
  }
  if (tfl[1] < (vox[1]-1))
    dist[1] = y - (float)tfl[1];
  else   // border values need special treatment
  {
    --tfl[1];
    dist[1] = 1.0f;
  }
  if (tfl[2] < (vox[2]-1))
    dist[2] = z - (float)tfl[2];
  else   // border values need special treatment
  {
    --tfl[2];
    dist[2] = 1.0f;
  }

  raw.makeCurrent(f);
  rd = raw.getData();   // get pointer to voxel data

  // Compute pointers to neighboring voxels:
  sliceSize = vox[0] * vox[1] * bpv;
  lineSize  = vox[0] * bpv;
  neighbor[0] = rd + bpv * tfl[0] + tfl[1] * lineSize + tfl[2] * sliceSize;
  neighbor[1] = neighbor[0] + lineSize;
  neighbor[2] = neighbor[1] + bpv;
  neighbor[3] = neighbor[0] + bpv;
  neighbor[4] = neighbor[0] + sliceSize;
  neighbor[5] = neighbor[1] + sliceSize;
  neighbor[6] = neighbor[2] + sliceSize;
  neighbor[7] = neighbor[3] + sliceSize;

  // Prepare passes, depending on voxel data type:
  switch (bpv)
  {
    case 1: // fall through
    case 2: passes = 1; break;
    case 3: passes = 3; break;
    case 4: passes = 4; break;
    default: return;
  }

  for (j=0; j<passes; ++j)
  {
    // Get neighboring voxel values:
    for (i=0; i<8; ++i)
    {
      switch (bpv)
      {
        default:
        case 1: val[i] = (float)(*neighbor[i]); break;
        case 2: val[i] = ((float)(*neighbor[i])) * 256.0f + (float)(*(neighbor[i]+1)); break;
        case 3: // fall through
        case 4: val[i] = (float)(*(neighbor[i] + j)); break;
      }
    }

    // Trilinearly interpolate values:
    interpolated = (int)(
      val[0] * (1.0f - dist[0]) * (1.0f - dist[1]) * (1.0f - dist[2]) +
      val[1] * (1.0f - dist[0]) * dist[1]          * (1.0f - dist[2]) +
      val[2] * dist[0]          * dist[1]          * (1.0f - dist[2]) +
      val[3] * dist[0]          * (1.0f - dist[1]) * (1.0f - dist[2]) +
      val[4] * (1.0f - dist[0]) * (1.0f - dist[1]) * dist[2] +
      val[5] * (1.0f - dist[0]) * dist[1]          * dist[2] +
      val[6] * dist[0]          * dist[1]          * dist[2] +
      val[7] * dist[0]          * (1.0f - dist[1]) * dist[2] );

    switch (bpv)
    {
      default:
      case 1: result[0] = (uchar)interpolated; break;
      case 2: result[0] = (uchar)(interpolated / 256);
              result[1] = (uchar)(interpolated % 256); break;
      case 3: // fall through 
      case 4: result[j] = (uchar)interpolated; break;
    }
  }
}

//----------------------------------------------------------------------------
/** Draws a 3D line into each animation frame of the dataset.
  Only works when bpv == 1!
  @param p1x,p1y,p1z   line start point
  @param p2x,p2y,p2z   line end point
  @param val           value of line voxels [0..255]
*/
void vvVolDesc::drawLine(int p1x, int p1y, int p1z, int p2x, int p2y, int p2z, 
  int val)
{
  uchar* raw;
  int f;        // frame counter

  vvDebugMsg::msg(3, "vvVolDesc::drawLine()");

  if (bpv!=1) return;
  for (f=0; f<frames; ++f)
  {
    raw = getRaw(f);
    vvToolshed::draw3DLine(p1x, p1y, p1z, p2x, p2y, p2z, (uchar)val, 
      raw, vox[0], vox[1], vox[2]);
  }
}

//----------------------------------------------------------------------------
/** Serializes all volume attributes to a memory buffer.
  Serialization is performed for volume size, number of voxels in each axis,
  number of time steps, real voxel sizes, etc.<BR>
  This serialization does NOT include the voxel data and the transfer functions!<BR>
  Since the serialization buffer must be allocated before calling the function,
  its size can be found by calling the function with the NULL parameter or no
  parameter at all. Here's an example:<PRE>
  int num_bytes = serializeAttributes();
  uchar* buffer = new uchar[num_bytes];
  serializeAttributes(buffer);
  </PRE><BR>
  The @see deserializeAttributes command can be used to deserialize the
  buffer values.<P>
  Here is the exact description of the buffer values in the serialized
  buffer:
  <PRE>
   Length          Data Type        VolDesc Attribute
 ---------------------------------------------------------------
   3 x 4 bytes     unsigned int     vox[0..2]
   4 bytes         unsigned int     frames
   1 byte          unsigned char    bpv (bytes per voxel)
   3 x 4 bytes     float            dist[0..2]
   4 bytes         float            dt
   2 x 4 bytes     float            realMin, realMax
   3 x 4 bytes     float            pos
   1 byte          unsigned char    storage type
  </PRE>
  @param buffer pointer to _allocated_ memory for serialized attributes
  @return number of bytes required for serialization buffer
*/
int vvVolDesc::serializeAttributes(uchar* buffer)
{
  uchar* ptr;     // pointer to current serialization buffer element

  vvDebugMsg::msg(3, "vvVolDesc::serializeAttributes()");

  if (buffer != NULL)
  {
    ptr = buffer;
    ptr += vvToolshed::write32 (ptr, vox[0]);    
    ptr += vvToolshed::write32 (ptr, vox[1]);    
    ptr += vvToolshed::write32 (ptr, vox[2]);    
    ptr += vvToolshed::write32 (ptr, frames);    
    ptr += vvToolshed::write8    (ptr, uchar(bpv));
    ptr += vvToolshed::writeFloat(ptr, dist[0]);   
    ptr += vvToolshed::writeFloat(ptr, dist[1]);   
    ptr += vvToolshed::writeFloat(ptr, dist[2]);   
    ptr += vvToolshed::writeFloat(ptr, dt);        
    ptr += vvToolshed::writeFloat(ptr, realMin);   
    ptr += vvToolshed::writeFloat(ptr, realMax);   
    ptr += vvToolshed::writeFloat(ptr, position.e[0]);  
    ptr += vvToolshed::writeFloat(ptr, position.e[1]);  
    ptr += vvToolshed::writeFloat(ptr, position.e[2]);  
    ptr += vvToolshed::write8    (ptr, uchar(stype));
    assert(ptr - buffer == SERIAL_ATTRIB_SIZE);
  }
  return SERIAL_ATTRIB_SIZE;
}

//----------------------------------------------------------------------------
/** Deserializes all volume attributes from a memory buffer.
  The @see serializeAttributes command can be used to serialize the
  attributes.
  @param buffer   pointer to _allocated_ memory for serialized attributes
  @param bufSize  size of buffer [bytes]. Values smaller than the default
                  size are allowed and only fill the values up to the 
                  passed value. The remaining values will be set to default values.
*/
void vvVolDesc::deserializeAttributes(uchar* buffer, int bufSize)
{
  uchar* ptr;       // pointer to current serialization buffer element

  vvDebugMsg::msg(3, "vvVolDesc::deserializeAttributes()");
  assert(buffer!=NULL);

  // Set default values for all serializable attributes:
  setDefaults();

  ptr = buffer;
  if (ptr+4 - buffer <= bufSize)
    vox[0] = vvToolshed::read32(ptr);  
  else return;
  ptr += 4;
  if (ptr+4 - buffer <= bufSize)
    vox[1] = vvToolshed::read32(ptr);  
  else return;
  ptr += 4;
  if (ptr+4 - buffer <= bufSize)
    vox[2] = vvToolshed::read32(ptr);  
  else return;
  ptr += 4;
  if (ptr+4 - buffer <= bufSize)
    frames = vvToolshed::read32(ptr);  
  else return;
  ptr += 4;
  if (ptr+1 - buffer <= bufSize)
    bpv = vvToolshed::read8(ptr);
  else return;
  ptr += 1;
  if (ptr+4 - buffer <= bufSize)
    dist[0] = vvToolshed::readFloat(ptr);  
  else return;
  ptr += 4;
  if (ptr+4 - buffer <= bufSize)
    dist[1]  = vvToolshed::readFloat(ptr);  
  else return;
  ptr += 4;
  if (ptr+4 - buffer <= bufSize)
    dist[2] = vvToolshed::readFloat(ptr);  
  else return;
  ptr += 4;
  if (ptr+4 - buffer <= bufSize)
    dt = vvToolshed::readFloat(ptr);  
  else return;
  ptr += 4;
  if (ptr+4 - buffer <= bufSize)
    realMin = vvToolshed::readFloat(ptr);  
  else return;
  ptr += 4;
  if (ptr+4 - buffer <= bufSize)
    realMax = vvToolshed::readFloat(ptr);  
  else return;
  ptr += 4;
  if (ptr+4 - buffer <= bufSize)
    position.e[0] = vvToolshed::readFloat(ptr);  
  else return;
  ptr += 4;
  if (ptr+4 - buffer <= bufSize)
    position.e[1] = vvToolshed::readFloat(ptr);  
  else return;
  ptr += 4;
  if (ptr+4 - buffer <= bufSize)
    position.e[2] = vvToolshed::readFloat(ptr);  
  else return;
  ptr += 4;
  if (ptr+1 - buffer <= bufSize)
    stype = StorageType(vvToolshed::read8(ptr));
  else return;
  ptr += 1;
}

//----------------------------------------------------------------------------
/** Set the data of one slice to new values.
  @param newData new voxel data for one specific slice 
                 (bpv*width*height voxels expected). The data must be arranged
                 in the same way as the volume data, i.e. all bytes of the voxel
                 must be in a row. The data will be copied to the volume dataset,
                 so it must be freed by the _caller_!
  @param slice   slice index (first slice = 0 = default)
  @param frame   frame index (first frame = 0 = default)
*/
void vvVolDesc::setSliceData(uchar* newData, int slice, int frame)
{
  uchar* dst;       // pointer to beginning of slice
  int sliceSize;    // shortcut for speed
  
  if (frames>0 && vox[2]>0)   // make sure at least one slice is stored
  {
    sliceSize = getSliceSize();
    dst = getRaw(frame) + slice * sliceSize;
    memcpy(dst, newData, sliceSize);
  }
}

//----------------------------------------------------------------------------
/** Extract a slice from a volume. Z direction slices are drawn just as they
  are in the volume. For X direction slices, the height is the same as for 
  Z slices, and the width is the depth of the volume. For Y direction slices,
  the width is the same as for the volume, and the height is the depth of
  the volume. The slices are drawn as they are seen from the positive end of
  the respective axis.
  @param frame  animation frame to use [0..frames-1]
  @param axis   axis for which to generate slice data (0=x, 1=y, 2=z)
  @param slice  slice index to create, relative to slicing axis (>=0) 
  @param dst    _allocated_ space for sliceWidth * sliceHeight * bpv bytes
*/
void vvVolDesc::extractSlice(int frame, int axis, int slice, uchar* dst)
{
  uchar* raw;       // raw volume data of current frame
  int sliceSize;    // bytes per volume slice (z-axis view)
  int lineSize;     // bytes per voxel line
  int i, j;

  sliceSize = vox[0] * vox[1] * bpv;
  lineSize  = vox[0] * bpv;
  raw = getRaw(frame);
  switch (axis)
  {
    case 0: // x axis
      for (j=0; j<vox[1]; ++j)
        for (i=0; i<vox[2]; ++i)
        {
          memcpy(dst, raw + i * sliceSize + j * lineSize + (vox[0] - slice - 1) * bpv, bpv);
          dst += bpv;
        }
      break;
    case 1: // y axis
      for (i=0; i<vox[2]; ++i)
        memcpy(dst + i * lineSize, raw + (vox[2] - i - 1) * sliceSize + slice * lineSize, lineSize);
      break;
    case 2: // z axis
    default:
      memcpy(dst, raw + slice * sliceSize, sliceSize);
      break;
  }
}

//----------------------------------------------------------------------------
/** Corrects an image with interlaced slices. The first slice will remain
  the first slice, the second slice will be taken from halfway into the
  dataset.
*/
void vvVolDesc::deinterlace()
{
  uchar* rd;
  uchar* volBuf;   // temporary buffer for voxel data of one time step
  uchar* dst;      // destination pointer
  uchar* src;      // source pointer
  int sliceSize;
  int frameSize;
  int z, f;
  
  vvDebugMsg::msg(2, "vvVolDesc::deinterlace()");
  if (simulate) return;   // this function does not affect the volume header values  

  sliceSize = getSliceSize();
  frameSize = getFrameSize();
  volBuf = new uchar[frameSize];
  assert(volBuf);
  raw.first();
  for (f=0; f<frames; ++f)
  {
    rd = raw.getData();
    memcpy(volBuf, rd, frameSize);  // make backup copy of volume
    for (z=0; z<vox[2]; ++z)
    {
      dst = rd + z * sliceSize;
      if ((z % 2) == 0)  // even slice number?
        src = volBuf + (z/2) * sliceSize;
      else
        src = volBuf + (z/2 + (vox[2]+1)/2) * sliceSize;
      
      // Swap source and destination slice:
      memcpy((void*)dst, (void*)src, sliceSize);
    }
    raw.next();
  }
  delete[] volBuf;
}
   
//----------------------------------------------------------------------------
/** Find the minimum and maximum scalar value. 
  For RGB data the min/max values of any component is returned, for RGBA data
  the min/max values of the alpha component is returned.
  @param scalarMin,scalarMax  minimum and maximum scalar values in volume animation
*/
void vvVolDesc::findMinMax(int* scalarMin, int* scalarMax)
{
  int f;
  int mi, ma;

  vvDebugMsg::msg(2, "vvVolDesc::findMinMax()");

  *scalarMin = (bpv==2) ? 65535 : 255;
  *scalarMax = 0;
  for (f=0; f<frames; ++f)
  {
    if (bpv==2) vvToolshed::getMinMax16bitBE(getRaw(f), getFrameVoxels(), &mi, &ma);
    else if (bpv==4) vvToolshed::getMinMaxAlpha(getRaw(f), getFrameVoxels(), &mi, &ma);
    else vvToolshed::getMinMax(getRaw(f), getFrameSize(), &mi, &ma);
    if (mi < *scalarMin) *scalarMin = mi;
    if (ma > *scalarMax) *scalarMax = ma;
  }
}

//----------------------------------------------------------------------------
/** Find the number of voxels with a specific value.<br>
  For RGB volumes, all three elements must be zero to be considered.
  For RGBA volumes, the alpha component must be zero to be considered.
  @param frame frame index to look at (first frame = 0)
  @param val   value to look for
*/
int vvVolDesc::findNumValue(int frame, int val)
{
  uchar* raw;
  int i;
  int num = 0;
  int frameSize;
  int sval;

  vvDebugMsg::msg(2, "vvVolDesc::findNumValue()");

  frameSize = getFrameSize();
  
  // Search volume:
  raw = getRaw(frame);
  for (i=0; i<frameSize; i+=bpv)
  {
    switch (bpv)
    {
      case 1:
        if (raw[i] == val) ++num;
        break;
      case 2:
        sval = int(raw[i]) * 256 + raw[i+1];
        if (sval == val) ++num;
        break;
      case 3:
        if (raw[i]==0 && raw[i+1]==0 && raw[i+2]==0) ++num;
        break;
      case 4:
        if (raw[3]==0) ++num;
        break;
      default: break;
    }
  }
  return num;
}

//----------------------------------------------------------------------------
/** Find the number of different data values used in a dataset.
*/
int vvVolDesc::findNumUsed()
{
  int f, i;
  bool* used;   // true = scalar value occurs in array
  uchar* raw;     
  int numValues;
  int numUsed = 0;
  int value;
  int frameSize;

  vvDebugMsg::msg(2, "vvVolDesc::findNumUsed()");

  frameSize = getFrameSize();
  numValues = (bpv==2) ? 65536 : 256;
  used = new bool[numValues];

  // Initialize occurrance array:
  for (i=0; i<numValues; ++i)
    used[i] = false;

  // Fill occurrance array:
  for (f=0; f<frames; ++f)
  {
    raw = getRaw(f);
    switch (bpv)
    {
      case 2:
        for (i=0; i<frameSize; i+=2)
        {
          value = (int(raw[i]) << 8) | int(raw[i+1]);
          used[value] = true;
        }
        break;
      default:
        for (i=0; i<frameSize; ++i)
          used[raw[i]] = true;
        break;
    }
  }

  // Count number of 'true' entries in occurrance array:
  for (i=0; i<numValues; ++i)
    if (used[i]==true) ++numUsed;

  delete[] used;
  return numUsed;
}

//----------------------------------------------------------------------------
/** Find the number of transparent data values in a dataset.
  If there is no transfer function defined in the VD,
  a linear ramp is assumed.
  @param frame frame index to look at (first frame = 0)
*/
int vvVolDesc::findNumTransparent(int frame)
{
  float* rgba = NULL;
  uchar* raw;
  int i;
  int numTransparent = 0;
  int frameSize;
  int lutEntries = 0;
  int scalar;
  bool noTF;      // true = no TF present in file

  vvDebugMsg::msg(2, "vvVolDesc::findNumTransparent()");

  frameSize = getFrameSize();
  noTF = tf.pins.isEmpty();
  
  if (!noTF)
  {
    switch(bpv)
    {
      case 1:  lutEntries = 256;   break;
      case 2:  lutEntries = 65536; break;
      case 3:  lutEntries = 0;     break;
      case 4:  lutEntries = 256;   break;
      default: lutEntries = 0;     break;
    }

    rgba = new float[4 * lutEntries];

    // Generate arrays from pins:
    tf.pins.setColorModel(vvPinList::RGB_MODEL);
    tf.pins.computeFunction(lutEntries, vvPinList::RGBA, rgba);
  }

  // Search volume:
  raw = getRaw(frame);
  switch (bpv)
  {
    case 1:
    case 2:
    case 4:
      for (i=bpv-1; i<frameSize; i += bpv)
      {
        if (bpv==2) scalar = (int(raw[i-1]) << 8) | int(raw[i]);
        else scalar = raw[i];
        if (noTF)
        {
          if (scalar==0) 
            ++numTransparent;
        }
        else
        {
          if (rgba[3 * lutEntries + scalar]==0.0f) 
            ++numTransparent;
        }
      }
      break;
    default:
      numTransparent = -1;
      break;
  }

  if (!noTF) delete[] rgba;

  return numTransparent;
}

//----------------------------------------------------------------------------
/** 
  Expand the data range to occupy the entire value range. This is useful to take
  advantage of the entire possible value range. For example, if the scanned
  16 bit data occupy only values between 50 and 1000, they will be mapped to
  a range of 0 to 65535.
  @param verbose true = display progress
*/
void vvVolDesc::expandDataRange(bool verbose)
{
  int smin, smax;   // scalar minimum and maximum

  vvDebugMsg::msg(2, "vvVolDesc::expandDataRange()");

  findMinMax(&smin, &smax);
  zoomDataRange(smin, smax, verbose);
}

//----------------------------------------------------------------------------
/** Zoom the data range to zoom in on a subset of the value range.
  @param low,high bottom and top limit of data range to zoom to [scalar value]
  @param verbose true = display progress
*/
void vvVolDesc::zoomDataRange(int low, int high, bool verbose)
{
  uchar* raw;       // raw volume data
  int frameSize;    // bytes per frame
  int f,i;
  int ival, irange;
  float fmin, fmax, fval, frange;
  float oldRealMin; // old real minimum
  
  vvDebugMsg::msg(2, "vvVolDesc::zoomDataRange()");

  if (simulate) return;

  frameSize = getFrameSize();
  fmin = float(low);
  fmax = float(high);
  frange = fmax - fmin;

  // Compute new real world range:
  irange = (bpv==2) ? 65535 : 255;
  oldRealMin = realMin;
  realMin = fmin / float(irange) * (realMax - oldRealMin) + oldRealMin;
  realMax = fmax / float(irange) * (realMax - oldRealMin) + oldRealMin;
  
  // Perform the actual expansion:
  if (verbose) vvToolshed::initProgress(frames);
  for (f=0; f<frames; ++f)
  {
    raw = getRaw(f);
    for (i=0; i<frameSize; ++i)
    {
      switch (bpv)
      {
        case 1:
        case 3:
        case 4:
          fval = float(*raw);
          ival = int(255.0f * (fval-fmin) / frange);
          ival = ts_clamp(ival, 0, 255);
          *raw = uchar(ival);
          ++raw;
          break;
        case 2:
          ival = (int(*raw) << 8) | int(*(raw+1));
          fval = float(ival);
          ival = int(65535.0f * (fval - fmin) / frange);
          ival = ts_clamp(ival, 0, 65535);
          *raw = uchar(ival >> 8);
          ++raw;
          *raw = uchar(ival & 0xFF);
          ++raw;
          ++i;  // additionally increase i
          break;
        default: break;
      }    
    }
    if (verbose) vvToolshed::printProgress(f);
  }
  if (verbose) cerr << endl;
}

//----------------------------------------------------------------------------
/** Blend two volumes together. Both volumes must have the same size and
  the same data type. If the number of time steps is different,
  the sequence with less steps will be repeated. The result will have
  the same number of frames as the source volume.
  @param blendVD volume to blend
  @param alpha   alpha value applied to voxels of blend volume [0..1],
                 source volume voxels are multiplied by (1-alpha)
  @param verbose true = display progress
*/
void vvVolDesc::blend(vvVolDesc* blendVD, float alpha, bool verbose)
{
  uchar* raw;       // raw volume data
  uchar* rawBlend;  // raw volume data of file to blend
  int frameSize;    // bytes per frame
  int f, i, fBlend;
  float val1, val2;
  int blended;      // result from blending operation

  vvDebugMsg::msg(2, "vvVolDesc::blend()");
  
  if (bpv != blendVD->bpv || vox[0] != blendVD->vox[0] ||
      vox[1] != blendVD->vox[1] || vox[2] != blendVD->vox[2] ||
      dist[0] != blendVD->dist[0] || dist[1] != blendVD->dist[1] ||
      dist[2] != blendVD->dist[2])
  {
    cerr << "Cannot blend: volumes to blend must have same size, bytes per voxel, and voxel distances" << endl;
    return;
  }

  frameSize = getFrameSize();
  fBlend = 0;

  if (verbose) vvToolshed::initProgress(frames);
  for (f=0; f<frames; ++f)
  {
    raw = getRaw(f);
    rawBlend = blendVD->getRaw(fBlend);
    for (i=0; i<frameSize; ++i)
    {
      switch(bpv)
      {
        case 1:
        case 3:
        case 4:
          val1 = (1.0f - alpha) * float(raw[i]);
          val2 = alpha * float(rawBlend[i]);
          blended = ts_clamp(int(val1 + val2), 0, 255);
          raw[i] = uchar(blended);
          break;
        case 2:
          if (i%2 == 0)   // only act on every other byte because of 16 bit values
          {
            val1 = (1.0f - alpha) * float(int(raw[i]) * 256 + int(raw[i+1]));
            val2 = alpha * float(int(rawBlend[i]) * 256 + int(rawBlend[i+1]));
            blended = ts_clamp(int(val1 + val2), 0, 65535);
            raw[i] = uchar(blended >> 8);
            raw[i+1] = uchar(blended & 0xff);
          }
          break;
      }
    }
    if (verbose) vvToolshed::printProgress(f);
    ++fBlend;
    if (fBlend>=blendVD->frames) fBlend = 0;
  }
  if (verbose) cerr << endl;
}

///// EOF /////
