//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifdef WIN32
  #include <windows.h>
  #include <float.h>
#else
  #include <unistd.h>
#endif
#ifdef __linux
  #include <values.h>
#endif
#include <stdio.h>
#include <string.h>
#ifdef _STANDARD_C_PLUS_PLUS
#include <iostream>
#include <iomanip>
using std::cerr;
using std::setw;
#else
#include <iostream.h>
#include <iomanip.h>
#endif
#include <float.h>
#include <math.h>
#include <assert.h>
#include "vvtoolshed.h"

#ifdef __sun
  #define ceilf ceil
  #define powf pow 
  #define atanf atan
  #define sqrtf sqrt
  #define sinf sin
  #define cosf cos
#endif

//#define VV_STANDALONE      // define to perform self test

int vvToolshed::progressSteps = 0;

//============================================================================
// Method Definitions
//============================================================================

//----------------------------------------------------------------------------
/** Case insensitive string comparison
    @param str1,str2 pointers to strings that are being compared
    @return
      <UL>
        <LI> 0 if equal
        <LI>-1 if str1<str2
        <LI> 1 if str1>str2
      </UL>
*/
int vvToolshed::strCompare(const char* str1, const char* str2)
{
  #ifdef WIN32
    return stricmp(str1, str2);
  #else
    return strcasecmp(str1, str2);
  #endif
}

//----------------------------------------------------------------------------
/** Case insensitive string comparison with a number of characters.
    @param str2,str2 pointers to strings that are being compared
    @param n = number of characters to compare
    @return the same values as in #strCompare(const char*, const char*)
*/  
int vvToolshed::strCompare(const char* str1, const char* str2, int n)
{
  #ifdef WIN32
    return strnicmp(str1, str2, n);
  #else
    return strncasecmp(str1, str2, n);
  #endif
}

//----------------------------------------------------------------------------
/** Case insensitive string suffix comparison.
    @param str    pointer to string
    @param suffix pointer to suffix
    @return true if suffix is the suffix of str
*/
bool vvToolshed::isSuffix(const char* str, const char* suffix)
{
  if (vvToolshed::strCompare(str + strlen(str) - strlen(suffix), suffix) == 0)
    return true;
  else return false;
}

//----------------------------------------------------------------------------
/** Convert HSB color model to RGB.
    @param a hue (0..360) (becomes red)
    @param b saturation (0..1) (becomes green)
    @param c value = brightness (0..1) (becomes blue)
    @return RGB values in a,b,c
*/
void vvToolshed::HSBtoRGB(float* a, float* b, float* c)
{
  float red, green, blue;

  HSBtoRGB(*a, *b, *c, &red, &green, &blue);
  *a = red;
  *b = green;
  *c = blue;
}

//----------------------------------------------------------------------------
/** Convert HSB color model to RGB.
    @param h hue (0..360)
    @param s saturation (0..1)
    @param v value = brightness (0..1)
    @return RGB values (0..1) in r,g,b
*/
void vvToolshed::HSBtoRGB(float h, float s, float v, float* r, float* g, float* b)
{
  float f, p, q, t;
  int i;

  // Clamp values to their valid ranges:
  h = ts_clamp(h, 0.0f, 1.0f);
  s = ts_clamp(s, 0.0f, 1.0f);
  v = ts_clamp(v, 0.0f, 1.0f);

  // Convert hue:
  if (h == 1.0f) h = 0.0f;
  h *= 360.0f;

  if (s==0.0f)    // grayscale value?
  {
    *r = v;
    *g = v;
    *b = v;
  }
  else
  {
    h /= 60.0;
    i = (int)floor(h);
    f = h - i;
    p = v * (1.0f - s);
    q = v * (1.0f - (s * f));
    t = v * (1.0f - (s * (1.0f - f)));
    switch (i)
    {
      case 0: *r = v; *g = t; *b = p; break;
      case 1: *r = q; *g = v; *b = p; break;
      case 2: *r = p; *g = v; *b = t; break;
      case 3: *r = p; *g = q; *b = v; break;
      case 4: *r = t; *g = p; *b = v; break;
      case 5: *r = v; *g = p; *b = q; break;
    }
  }
}

//----------------------------------------------------------------------------
/** Convert RGB colors to HSB model
    @param a red [0..360] (becomes hue)
    @param b green [0..1] (becomes saturation)
    @param c blue [0..1]  (becomes brightness)
    @return HSB in a,b,c
*/
void vvToolshed::RGBtoHSB(float* a, float* b, float* c)
{
  float h,s,v;
  RGBtoHSB(*a, *b, *c, &h, &s, &v);
  *a = h;
  *b = s;
  *c = v;
}

//----------------------------------------------------------------------------
/** Convert RGB colors to HSB model.
    @param r,g,b RGB values [0..1]
    @return h = hue [0..1], s = saturation [0..1], v = value = brightness [0..1]
*/
void vvToolshed::RGBtoHSB(float r, float g, float b, float* h, float* s, float* v)
{
  float max, min, delta;
  
  // Clamp input values to valid range:
  r = ts_clamp(r, 0.0f, 1.0f);
  g = ts_clamp(g, 0.0f, 1.0f);
  b = ts_clamp(b, 0.0f, 1.0f);

  max = ts_max(r, ts_max(g, b));
  min = ts_min(r, ts_min(g, b));
  *v = max;
  *s = (max != 0.0f) ? ((max - min) / max) :0.0f;
  if (*s == 0.0f) *h = 0.0f;
  else
  {
    delta = max - min;
    if (r==max)
      *h = (g - b) / delta;
    else if (g==max)
      *h = 2.0f + (b - r) / delta;
    else if (b==max)
      *h = 4.0f + (r - g) / delta;
    *h *= 60.0f;
    if (*h < 0.0f)
      *h += 360.0f;
  }
  *h /= 360.0f;
}

//----------------------------------------------------------------------------
/** Copies the tail string after the last occurrence of a given character.
    Example: str="c:\ local\ testfile.dat", c='\' => suffix="testfile.dat"
    @param suffix <I>allocated</I> space for the found string
    @param str    source string
    @param c      character after which to copy characters
    @return result in suffix, empty string if c was not found in str
*/
void vvToolshed::strcpyTail(char* suffix, const char* str, char c)
{
  int i, j;

  // Search for c in pathname:
  i = strlen(str) - 1;
  while (i>=0 && str[i]!=c)
    --i;

  // Extract tail string:
  if (i<0)    // c not found?
    strcpy(suffix, str);
  else
  {
    for (j=i+1; j<(int)strlen(str); ++j)
      suffix[j-i-1] = str[j];
    suffix[j-i-1] = '\0';
  }
}

//----------------------------------------------------------------------------
/** Copies the head string before the first occurrence of a given character.
    Example: str="c:\ local\ testfile.dat", c='.' => head="c:\ local\ testfile"
    @param head  <I>allocated</I> space for the found string
    @param str    source string
    @param c      character before which to copy characters
    @return result in head, empty string if c was not found in str
*/
void vvToolshed::strcpyHead(char* head, const char* str, char c)
{
  int i = 0;

  if (strchr(str, c) == NULL) 
  {
    head[0] = '\0';
    return;
  }
  while (str[i] != c)
  {
    head[i] = str[i];
    ++i;
  }
  head[i] = '\0';
}

//----------------------------------------------------------------------------
/** Removes leading and trailing spaces from a string.
    Example: str="  hello " => str="hello"
    @param str    string to trim
    @return result in str
*/
void vvToolshed::strTrim(char* str)
{
  int i;

  // Trim trailing spaces:
  for (i=strlen(str)-1; i>0; --i)
  {
    if (str[i]==' ') str[i] = '\0';
    else break;
  }
  if (str[0]=='\0') return;   // done

  // Trim leading spaces:
  i=0;
  while (str[i]==' ')
  {
    ++i;
  }
  if (i==0) return;   // done
  strcpy(str, str+i);
}

//----------------------------------------------------------------------------
/** Extracts a filename from a given path.
    Directory elements have to be separated by '/' or '\' depending on OS.
    @param filename <I>allocated</I> space for filename (e.g. "testfile.dat")
    @param pathname file including entire path (e.g. "/usr/local/testfile.dat")
    @return result in filename
*/
void vvToolshed::extractFilename(char* filename, const char* pathname)
{
#ifdef WIN32
  strcpyTail(filename, pathname, '\\');
#else
  strcpyTail(filename, pathname, '/');
#endif
}

//----------------------------------------------------------------------------
/** Extracts a directory name from a given path.
    Directory elements have to be separated by '/' or '\' depending on OS.
    @param dirname  <I>allocated</I> space for directory name (e.g. "/usr/local/" or "c:\user\")
    @param pathname file including entire path (e.g. "/usr/local/testfile.dat" or "c:\user\testfile.dat")
    @return result in dirname
*/
void vvToolshed::extractDirname(char* dirname, const char* pathname)
{
  int i, j;
  char delim;   // path delimiter

#ifdef WIN32
  delim = '\\';
#else
  delim = '/';
#endif

  // Search for '\' or '/' in pathname:
  i = strlen(pathname) - 1;
  while (i>=0 && pathname[i]!=delim)
    --i;

  // Extract preceding string:
  if (i<0)    // delimiter not found?
    strcpy(dirname, "");
  else
  {
    for (j=0; j<=i; ++j)
      dirname[j] = pathname[j];
    dirname[j] = '\0';
  }
}

//----------------------------------------------------------------------------
/** Extracts an extension from a given path or filename.
    @param extension <I>allocated</I> space for extension (e.g. "dat")
    @param pathname  file including entire path (e.g. "/usr/local/testfile.dat")
    @return result in extension
*/
void vvToolshed::extractExtension(char* extension, const char* pathname)
{
  strcpyTail(extension, pathname, '.');
}

//----------------------------------------------------------------------------
/** Extracts the base file name from a given path or filename, excluding 
    the '.' delimiter.
    @param basename  <I>allocated</I> memory space for basename (e.g. "testfile").
                     Memory must be allocated for at least strlen(pathname)+1 chars!
    @param pathname  file including entire path (e.g. "/usr/local/testfile.dat")
    @return result in basename
*/
void vvToolshed::extractBasename(char* basename, const char* pathname)
{
  int i;

  extractFilename(basename, pathname);

  // Search for '.' in pathname:
  i = strlen(basename) - 1;
  while (i>=0 && basename[i]!='.')
    --i;

  if (i>0) basename[i] = '\0';    // convert '.' to '\0' to terminate string
}

//----------------------------------------------------------------------------
/** Remove the extension from a path string. If no '.' is present in path
    string, the path is removed without changes.
    @param basepath  <I>allocated</I> space for path without extension 
                     (e.g., "/usr/local/testfile")
    @param pathname  file including entire path (e.g. "/usr/local/testfile.dat")
    @return result in basepath
*/
void vvToolshed::extractBasePath(char* basepath, const char* pathname)
{
  int i, j;

  // Search for '.' in pathname:
  i = strlen(pathname) - 1;
  while (i>=0 && pathname[i]!='.')
    --i;

  // Extract tail string:
  if (i<0)    // '.' not found?
  {
    strcpy(basepath, pathname);
  }
  else
  {
    for (j=0; j<i; ++j)
      basepath[j] = pathname[j];
    basepath[j] = '\0';
  }
}

//----------------------------------------------------------------------------
/** Replaces a file extension with a new one, overwriting the old one.
    If the pathname does not have an extension yet, the new extension will be
    added.
    @param newPath      _allocated space_ for resulting path name with new 
                        extension (e.g. "/usr/local/testfile.txt")
    @param newExtension new extension without '.' (e.g. "txt")
    @param pathname     file including entire path (e.g. "/usr/local/testfile.dat")
    @return result in newPath
*/
void vvToolshed::replaceExtension(char* newPath, const char* newExtension, const char* pathname)
{
  char* pointPos;
  int baseNameLen;    // length of base file name, including point

  pointPos = strrchr(pathname, '.');
  if (pointPos==NULL) // is there a point in pathname?
  {
    // No point, so just add new extension:
    strcpy(newPath, pathname);
    strcat(newPath, ".");
    strcat(newPath, newExtension);
  }
  else
  {
    baseNameLen = pointPos-pathname+1;
    memcpy(newPath, pathname, baseNameLen); // copy everything before the point, including the point
    newPath[baseNameLen] = '\0';
    strcat(newPath, newExtension);
  }
}

//----------------------------------------------------------------------------
/** Increases the filename (filename must include an extension!)
 @return true if successful, false if filename couldn't be 
         increased
*/
bool vvToolshed::increaseFilename(char* filename)
{
  bool done = false;
  int i;
  char ext[256];
  
  extractExtension(ext, filename);
  if (strlen(ext)==0) i=strlen(filename) - 1;
  else i = strlen(filename) - strlen(ext) - 2;
  while (!done)
  {
    if (i<0 || filename[i]<'0' || filename[i]>'9') 
      return false;

    if (filename[i] == '9')   // overflow?
    {
      filename[i] = '0';
      --i;
    }
    else 
    {
      ++filename[i];
      done = 1;
    } 
  }
  return true;
}

//----------------------------------------------------------------------------
/** Draws a line in a 3D volume dataset using Bresenham's algorithm.
    The volume must consist of one byte per voxel!
    Both line end points must lie within the volume. The Coordinate system is:
    <PRE>
           y  
           |__ x
          / 
         z
    </PRE>
    The volume data is arranged like this:
    <UL>
      <LI>origin is top left front
      <LI>width in positive x direction
      <LI>height in negative y direction
      <LI>slices in negative z direction
    </UL>
    @param x0,y0,z0  line starting point in voxels
    @param x1,y1,z1  line end point in voxels
    @param scalar    line color
    @param data      pointer to raw volume data
    @param w,h,s     width/height/slices of volume data array in voxels
*/
void vvToolshed::draw3DLine(int x0, int y0, int z0, int x1, int y1, int z1, 
                          uchar scalar, uchar* data, int w, int h, int s)
{
  int xd, yd, zd;
  int x, y, z;
  int ax, ay, az;
  int sx, sy, sz;
  int dx, dy, dz;

  x0 = ts_clamp(x0, 0, w-1);
  x1 = ts_clamp(x1, 0, w-1);
  y0 = ts_clamp(y0, 0, h-1);
  y1 = ts_clamp(y1, 0, h-1);
  z0 = ts_clamp(z0, 0, s-1);
  z1 = ts_clamp(z1, 0, s-1);

  dx = x1 - x0;
  dy = y1 - y0;
  dz = z1 - z0;

  ax = ts_abs(dx) << 1;
  ay = ts_abs(dy) << 1;
  az = ts_abs(dz) << 1;

  sx = ts_zsgn(dx);
  sy = ts_zsgn(dy);
  sz = ts_zsgn(dz);

  x = x0;
  y = y0;
  z = z0;

  if (ax >= ts_max(ay, az))            // x is dominant
  {
    yd = ay - (ax >> 1);
    zd = az - (ax >> 1);
    for (;;)
    {
      data[z * w * h + y * w + x] = scalar;
      if (x == x1) return;
      if (yd >= 0)
      {
        y += sy;
        yd -= ax;
      }
      if (zd >= 0)
      {
        z += sz;
        zd -= ax;
      }
      x += sx;
      yd += ay;
      zd += az;
    }
  }
  else if (ay >= ts_max(ax, az))            // y is dominant 
  {
    xd = ax - (ay >> 1);
    zd = az - (ay >> 1);
    for (;;)
    {
      data[z * w * h + y * w + x] = scalar;
      if (y == y1) return;
      if (xd >= 0)
      {
        x += sx;
        xd -= ay;
      }
      if (zd >= 0)
      {
        z += sz;
        zd -= ay;
      }
      y += sy;
      xd += ax;
      zd += az;
    }
  }
  else if (az >= ts_max(ax, ay))            // z is dominant 
  {
    xd = ax - (az >> 1);
    yd = ay - (az >> 1);
    for (;;)
    {
      data[z * w * h + y * w + x] = scalar;
      if (z == z1) return;
      if (xd >= 0)
      {
        x += sx;
        xd -= az;
      }
      if (yd >= 0)
      {
        y += sy;
        yd -= az;
      }
      z += sz;
      xd += ax;
      yd += ay;
    }
  }
}

//----------------------------------------------------------------------------
/** Draws a line in a 2D image dataset using Bresenham's algorithm.
    Both line end points must lie within the image. The coordinate system is:
    <PRE>
           y  
           |__ x
    </PRE>
    The image data is arranged like this:
    <UL>
      <LI>origin is top left
      <LI>width is in positive x direction
      <LI>height is in negative y direction
    </UL>
    @param x0/y0  line starting point in pixels
    @param x1/y1  line end point in pixels
    @param color  line color, 32 bit value: bits 0..7=first color, 
                  8..15=second color etc.
    @param data   pointer to raw image data
    @param bpp    byte per pixel (e.g. 3 for 24 bit RGB), range: [1..4]
    @param w/h    width/height of image data array in pixels
*/
void vvToolshed::draw2DLine(int x0, int y0, int x1, int y1, 
                          uint color, uchar* data, int bpp, int w, int h)
{
  int xd, yd;
  int x, y;
  int ax, ay;
  int sx, sy;
  int dx, dy;
  int i;
  uchar col[4];     // color components; 0=most significant byte

  assert(bpp <= 4);

  col[0] = (uchar)((color >> 24) & 0xff);
  col[1] = (uchar)((color >> 16) & 0xff);
  col[2] = (uchar)((color >> 8)  & 0xff);
  col[3] = (uchar)(color & 0xff);

  x0 = ts_clamp(x0, 0, w-1);
  x1 = ts_clamp(x1, 0, w-1);
  y0 = ts_clamp(y0, 0, h-1);
  y1 = ts_clamp(y1, 0, h-1);

  dx = x1 - x0;
  dy = y1 - y0;

  ax = ts_abs(dx) << 1;
  ay = ts_abs(dy) << 1;

  sx = ts_zsgn(dx);
  sy = ts_zsgn(dy);

  x = x0;
  y = y0;

  if (ax >= ay)            // x is dominant
  {
    yd = ay - (ax >> 1);
    for (;;)
    {
      for (i=0; i<bpp; ++i)
        data[bpp * (y * w + x) + i] = col[i];
      if (x == x1) return;
      if (yd >= 0)
      {
        y += sy;
        yd -= ax;
      }
      x += sx;
      yd += ay;
    }
  }
  else if (ay >= ax)            // y is dominant 
  {
    xd = ax - (ay >> 1);
    for (;;)
    {
      for (i=0; i<bpp; ++i)
        data[bpp * (y * w + x) + i] = col[i];
      if (y == y1) return;
      if (xd >= 0)
      {
        x += sx;
        xd -= ay;
      }
      y += sy;
      xd += ax;
    }
  }
}

//----------------------------------------------------------------------------
/** Compute texture hardware compatible numbers.
    @param imgSize  the image size [pixels]
    @return the closest power-of-2 value that is greater than or equal to imgSize.
*/
int vvToolshed::getTextureSize(int imgSize)
{
  return (int)powf(2.0f, (float)ceil(log((float)imgSize) / log(2.0f)));
}

//----------------------------------------------------------------------------
/** Checks if a file exists.
    @param filename file name to check for
    @return true if file exists
*/
bool vvToolshed::isFile(const char* filename)
{
  FILE* fp;

  fp = fopen(filename, "rb");
  if (fp==NULL) return false;
  fclose(fp);
  return true;
}

//----------------------------------------------------------------------------
/** Figures out the size of a file in bytes.
    @param  filename file name including path
    @return file size in bytes or -1 on error
*/
long vvToolshed::getFileSize(const char* filename)
{
  FILE* fp;
  long size;

  fp = fopen(filename, "rb");
  if (fp==NULL) return -1;
  if (fseek(fp, 0L, SEEK_END) != 0)
  {
    fclose(fp);
    return -1;
  }
  size = ftell(fp);
  fclose(fp);
  return size;
}

//----------------------------------------------------------------------------
/** Find the minimum and the maximum values in an uchar data array.
    @param data        source data array
    @param elements    number of bytes in source array
    @return minimum and maximum in min and max
*/
void vvToolshed::getMinMax(const uchar* data, int elements, int* min, int* max)
{
  int i;

  *min = 255;
  *max = 0;

  for (i=0; i<elements; ++i)
  {
    if (data[i] > *max) *max = data[i];
    if (data[i] < *min) *min = data[i];
  }
}

//----------------------------------------------------------------------------
/** Find the minimum and the maximum values in a 16 bit big endian data array.
    @param data        source data array
    @param elements    number of 16 bit elements in source array
    @return minimum and maximum in min and max
*/
void vvToolshed::getMinMax16bitBE(const uchar* data, int elements, int* min, int* max)
{
  int i;
  int value;
  int bytes;

  *min = 65535;
  *max = 0;
  bytes = 2 * elements;

  for (i=0; i<bytes; i+=2)
  {
    value = (int(data[i]) << 8) | int(data[i+1]);
    if (value > *max) *max = value;
    if (value < *min) *min = value;
  }
}

//----------------------------------------------------------------------------
/** Find the minimum and the maximum values in an RGBA dataset, only 
    considering the alpha component.
    @param data        source data array
    @param elements    number of RGBA elements in source array
    @return minimum and maximum of the alpha component in min and max
*/
void vvToolshed::getMinMaxAlpha(const uchar* data, int elements, int* min, int* max)
{
  int i;
  int bytes;

  *min = 255;
  *max = 0;
  bytes = 4 * elements;

  for (i=3; i<bytes; i+=4)
  {
    if (data[i] > *max) *max = data[i];
    if (data[i] < *min) *min = data[i];
  }
}

//----------------------------------------------------------------------------
/** Find the minimum and the maximum values in a float data array.
    @param data        source data array
    @param elements    number of elements in source array
    @return minimum and maximum in min and max
*/
void vvToolshed::getMinMax(const float* data, int elements, float* min, float* max)
{
  int i;

  *min = FLT_MAX;
  *max = -(*min);

  for (i=0; i<elements; ++i)
  {
    if (data[i] > *max) *max = data[i];
    if (data[i] < *min) *min = data[i];
  }
}

//----------------------------------------------------------------------------
/** Find the minimum and the maximum values in a data array and specify
    a value which is to be ignored, i.e., it does not change the determined
    minimum or maximum values.
    @param data        source data array
    @param elements    number of elements in source array
    @param ignore      value which is to be ignored (e.g, FLT_MAX)
    @return minimum and maximum in min and max
*/
void vvToolshed::getMinMaxIgnore(const float* data, int elements, float ignore, 
  float* min, float* max)
{
  int i;

  *min = FLT_MAX;
  *max = -(*min);

  for (i=0; i<elements; ++i)
  {
    if (data[i] != ignore)
    {
      if (data[i] > *max) *max = data[i];
      if (data[i] < *min) *min = data[i];
    }
  }
}

//----------------------------------------------------------------------------
/** Convert a sequence of uchar values to float.
    Be sure to have the floatArray allocated before this call!
    @param ucharArray  source array
    @param floatArray  destination array
    @param elements    number of uchar array elements to convert
    @return result in floatArray
*/
void vvToolshed::convertUChar2Float(const uchar* ucharArray, float* floatArray, int elements)
{
  int i;
  for (i=0; i<elements; ++i)
    floatArray[i] = (float)((double)ucharArray[i] / 255.0);
}

//----------------------------------------------------------------------------
/** Convert a sequence of float values to uchar.
    The uchar values cover a range of 0.0 to 1.0 of the float values.
    Be sure to have the ucharArray allocated before this call!
    @param floatArray  source array
    @param ucharArray  destination array
    @param elements    number of float array elements to convert
    @return result in ucharArray
*/
void vvToolshed::convertFloat2UChar(const float* floatArray, 
    uchar* ucharArray, int elements)
{
  for (int i=0; i<elements; ++i)
    ucharArray[i] = (uchar)(255.0 * floatArray[i]);
}

//----------------------------------------------------------------------------
/** Convert a sequence of float values to uchar.
    The uchar values will cover the range defined by the minimum and 
    maximum float values.
    Be sure to have the ucharArray allocated before this call!
    @param floatArray  source array
    @param ucharArray  destination array
    @param elements    number of float array elements to convert
    @param min,max     minimum and maximum float values which will be assigned
                       to 0 and 255 respectively.
    @return result in ucharArray
*/
void vvToolshed::convertFloat2UCharClamp(const float* floatArray, 
    uchar* ucharArray, int elements, float min, float max)
{
  int i;

  if (min>=max) 
    memset(ucharArray, 0, elements);
  else
  {
    for (i=0; i<elements; ++i)
      ucharArray[i] = (uchar)(255.0f * (floatArray[i] - min) / (max - min));
  }
}

//----------------------------------------------------------------------------
/** Convert a sequence of float values to 16 bit values.
    The 16 bit values will cover the range defined by the minimum and 
    maximum float values.
    Be sure to have the ucharArray allocated before this call (requires
    elements * 2 bytes)!
    @param floatArray  source array
    @param ucharArray  destination array (16 bit values)
    @param elements    number of float array elements to convert
    @param min,max     minimum and maximum float values which will be assigned
                       to 0 and 65535 respectively.
    @return result in ucharArray
*/
void vvToolshed::convertFloat2ShortClamp(const float* floatArray, 
    uchar* ucharArray, int elements, float min, float max)
{
  int i;
  int shortValue;

  if (min>=max) 
  {
    memset(ucharArray, 0, elements);
  }
  else
  {
    for (i=0; i<elements; ++i)
    {
      shortValue = int(65535.0f * (floatArray[i] - min) / (max - min));
      ucharArray[2*i]   = (uchar)((shortValue >> 8) & 255);
      ucharArray[2*i+1] = (uchar)(shortValue & 255);
    }
  }
}

//----------------------------------------------------------------------------
/** Convert a sequence of float values to uchar.
    The uchar values will cover the range defined by the maximum 
    and minimum float values.
    The uchar value of 0 is set for all float values of the 'zero' value,
    and only for them.
    The minimum float value which is not 'zero' becomes an uchar value of 1.
    Be sure to have the ucharArray allocated before this call!
    @param floatArray  source array
    @param ucharArray  destination array
    @param elements    number of float array elements to convert
    @param min,max     minimum and maximum float values which will be assigned
                       to 1 and 255 respectively.
    @param zero        float value which will become 0 in uchar array
    @return result in ucharArray
*/
void vvToolshed::convertFloat2UCharClampZero(const float* floatArray, 
    uchar* ucharArray, int elements, float min, float max, float zero)
{
  int i;
  
  if (min>=max)       
    memset(ucharArray, 0, elements);
  else
    for (i=0; i<elements; ++i)
    {
      if (floatArray[i] == zero) 
        ucharArray[i] = (uchar)0;
      else
      {
        if (min==max)
          ucharArray[i] = (uchar)1;
        else
          ucharArray[i] = (uchar)((uchar)(254.0f * (floatArray[i] - min) / (max - min)) + (uchar)1);
      }
    }
}

//----------------------------------------------------------------------------
/** Compute the largest prime factor which is not the number itself.
    @param number number to examine (>1)
    @return largest prime factor, -1 on error
*/
int vvToolshed::getLargestPrimeFactor(int number)
{
  int remainder;
  int factor = 2, largest = 1;

  if (number < 2) return -1;
  remainder = number;
  while (factor < remainder/2)
  {
    if ((remainder % factor) == 0)
    {
      remainder /= factor;
      largest = factor;
    }
    else 
      ++factor;
  }
  if (largest==1) return 1;
  else return ts_max(remainder, largest);
}

//----------------------------------------------------------------------------
/** Round the float value to the nearest integer value.
    @param fval value to round
    @return rounded value
*/
int vvToolshed::round(float fval)
{
  return (int)floor(fval + 0.5f);
}

//----------------------------------------------------------------------------
/** Initialize progress display.
    @param total total number of progress steps
*/
void vvToolshed::initProgress(int total)
{
  progressSteps = total;
  cerr << "     ";
}

//----------------------------------------------------------------------------
/** Print progress.
  Format: 'xxx %'.
  @param current current progress step
*/
void vvToolshed::printProgress(int current)
{
  int percent, i;

  if (progressSteps<2) percent = 100;
  else
    percent = 100 * current / (progressSteps - 1);
  for (i=0; i<5; ++i)
    cerr << (char)8;      // ASCII 8 = backspace (BS)
  cerr << setw(3) << percent << " %";
}

//----------------------------------------------------------------------------
/** Run length encode (RLE) a sequence in memory.
  Encoding scheme: X is first data chunk (unsigned char).<UL>
  <LI>if X<128:  copy next X+1 chunks (literal run)</LI>
  <LI>if X>=128: repeat next chunk X-127 times (replicate run)</LI></UL>
  @param out  destination position in memory (must be _allocated_!)
  @param in   source location in memory
  @param size number of bytes to encode
  @param symbol_size  bytes per chunk
  @param space  number of bytes allocated for destination array.
                Encoding process is stopped when this number is reached.
  @return number of bytes written to destination memory or -1 if there is not
          enough destination memory, -2 if an invalid data size was passed       
  @see decodeRLE
  @author Michael Poehnl
*/
int vvToolshed::encodeRLE(uchar* out, uchar* in, int size, int symbol_size, int space) 
{
  int same_symbol=1;
  int diff_symbol=0;
  int src=0;
  int dest=0;
  bool same;
  int i;

  if ((size % symbol_size) != 0)
  {
    return -2;
  }

  while (src < (size - symbol_size))
  {
    same = true;
    for (i=0; i<symbol_size; i++)
    {
      if (in[src+i] != in[src+symbol_size+i])
      {
        same = false;
        break;
      }
    }
    if (same)
    {
      if (same_symbol == 129)
      {
        // TODO: hier testen ob dest noch innerhalb von space ist
        out[dest] = (uchar)(126+same_symbol);      
        dest += symbol_size+1;
        same_symbol = 1;  
      }
      else
      {
        same_symbol++;
        if (diff_symbol > 0)
        {
          // TODO: hier testen ob dest noch innerhalb von space ist
          out[dest] = (uchar)(diff_symbol-1);
          dest += 1+symbol_size*diff_symbol;
          diff_symbol=0;  
        }
        if (same_symbol == 2)
        {
          if ((dest+1+symbol_size) > space)
          {
            return -1;
          }                  
          memcpy(&out[dest+1], &in[src], symbol_size);
        }
      }
    }
    else
    {
      if (same_symbol > 1)
      {
        // TODO: hier testen ob dest noch innerhalb von space ist
        out[dest] = (uchar)(126+same_symbol);      
        dest += symbol_size+1;
        same_symbol = 1;           
      }
      else
      {
        if ((dest+1+diff_symbol*symbol_size+symbol_size) > space)
        {
          return -1;
        }                              
        memcpy(&out[dest+1+diff_symbol*symbol_size], &in[src], symbol_size);           
        diff_symbol++;
        if (diff_symbol == 128)
        {
          // TODO: hier testen ob dest noch innerhalb von space ist
          out[dest] = (uchar)(diff_symbol-1);
          dest += 1+symbol_size*diff_symbol;
          diff_symbol=0;             
        }
      }
    }
    src += symbol_size;
  }
  if (same_symbol > 1)
  {
    // TODO: hier testen ob dest noch innerhalb von space ist
    out[dest] = (uchar)(126+same_symbol);      
    dest += symbol_size+1;
  }
  else
  {
    if ((dest+1+diff_symbol*symbol_size+symbol_size) > space)
    {
      return -1;
    }                   
    memcpy(&out[dest+1+diff_symbol*symbol_size], &in[src], symbol_size);           
    diff_symbol++;            
    out[dest] = (uchar)(diff_symbol-1);
    dest += 1+symbol_size*diff_symbol;
  }
  return dest;  
}   

//----------------------------------------------------------------------------
/** Decode a run length encoded (RLE) sequence. 
  Data chunks of any byte aligned size can be processed.
  @param out  destination position in memory (_allocated_ space for max bytes)
  @param in   source location in memory
  @param size number of bytes in source array to decode
  @param symbol_size  bytes per chunk (e.g., to encode 24 bit RGB data, use bpc=3)
  @param space  number of allocated bytes in destination memory (for range checking)
  @return number of bytes written to destination memory. If max would
          have been exceeded, -1 is returned
  @see encodeRLE
  @author Michael Poehnl
*/
int vvToolshed::decodeRLE(uchar* out, uchar* in, int size, int symbol_size, int space)
{   
  int src=0;
  int dest=0;
  int i, length;

  while (src < size)
  {   
    length = (int)in[src];
    if (length > 127)
    {
      for(i=0; i<(length - 126); i++)
      {
        if ((dest + symbol_size) > space)
        {
          return -1;
        }                      
        memcpy(&out[dest], &in[src+1], symbol_size);
        dest += symbol_size;
      }
      src += 1+symbol_size;
    }
    else
    {
      length++;
      if ((dest + length*symbol_size) > space)
      {
        return -1;
      }                                     
      memcpy(&out[dest], &in[src+1], symbol_size*length);
      dest += length*symbol_size; 	         
      src += 1+symbol_size*length;
    }
  }
  return dest;
}
   
//----------------------------------------------------------------------------
/** Run length encode (RLE) a sequence of 8 bit values in memory.
  Encoding scheme: X is first data byte (unsigned char).<UL>
  <LI>if X<128:  copy next X+1 bytes (literal run)</LI>
  <LI>if X>=128: repeat next byte X-127 times (replicate run)</LI></UL>
  @param dst  destination position in memory (must be _allocated_!)
  @param src  source location in memory
  @param len  number of bytes to encode
  @param max  number of bytes allocated in destination memory.
              Encoding process is stopped when this number is reached
  @return number of bytes written to destination memory
  @see decodeRLEFast
*/
int vvToolshed::encodeRLEFast(uchar* dst, uchar* src, int len, int max)
{
  int offset;     // start position of currently processed run in source array
  int index;      // index in source array
  int out;        // index in destination array
  int i;          // counter
  uchar cur;      // currently processed data byte

  offset = out = 0;
  while (offset < len) 
  {
    index = offset;
    cur = src[index++];   // load first data byte from source array
    while (index<len && index-offset<128 && src[index]==cur)
      index++;    // search for replicate run
    if (index-offset==1)  // generate literal run
    {
      // Failed to "replicate" the current byte. See how many to copy.
      // Avoid a replicate run of only 2-pixels after a literal run. There
      // is no gain in this, and there is a risk of loss if the run after
      // the two identical pixels is another literal run. So search for
      // 3 identical pixels.
      while (index<len && index-offset<128 && (src[index]!=src[index-1] || 
             index>1 && src[index]!=src[index-2]))
        index++;
       // Check why this run stopped. If it found two identical pixels, reset
       // the index so we can add a run. Do this twice: the previous run
       // tried to detect a replicate run of at least 3 pixels. So we may be
       // able to back up two pixels if such a replicate run was found.
       while (index<len && src[index]==src[index-1])
         index--;
       if (out < max)
         dst[out++] = (uchar)(index - offset - 1);
       for (i=offset; i<index; i++)
         if (out < max) dst[out++] = src[i];
    } 
    else    // generate replicate run
    {
      if (out < max)
        dst[out++] = (uchar)(index - offset + 127);
      if (out < max)
        dst[out++] = cur;
    }
    offset = index;
  } 
  return out;
}

//----------------------------------------------------------------------------
/** Decode a run length encoded (RLE) sequence of 8 bit values.
  @param dst  destination position in memory (must be _allocated_!)
  @param src  source location in memory
  @param len  number of bytes to decode
  @param max  number of allocated bytes in destination memory (for range checking)
  @return number of bytes written to destination memory. If max would
          have been exceeded, max+1 is returned
  @see encodeRLEFast
*/
int vvToolshed::decodeRLEFast(uchar* dst, uchar* src, int len, int max)
{
  int count;    // RLE counter
  int out=0;      // counter for written output bytes

  while (len > 0) 
  {
    count = (int)*src++;
    if (count > 127)      // replicate run?
    {
      count -= 127;           // remove bias
      if (out+count <= max)    // don't exceed allocated memory array
        memset(dst, *src++, count);
      else
      {
        if (out < max)    // write as much as possible
          memset(dst, *src++, max-out);
        return max+1;
      }
      len -= 2;
    }
    else                  // literal run
    { 
      ++count;    // remove bias
      if (out+count <= max)   // don't exceed allocated memory array
        memcpy(dst, src, count);
      else
      {
        if (out < max)    // write as much as possible
          memcpy(dst, src, max-out);
        return max+1;
      }
      src += count;
      len -= count + 1;
    }
    dst += count;
    out += count;
  }
  return out;
}

//----------------------------------------------------------------------------
/** Get the number of system processors.
  @return number of processors
*/
int vvToolshed::getNumProcessors()
{
#ifdef WIN32
  SYSTEM_INFO sysinfo;
  GetSystemInfo(&sysinfo);
  return sysinfo.dwNumberOfProcessors;
#else
  char hostName[128];
  gethostname(hostName, 128);
  // TODO: find out how to do that in Unix
  if      (strCompare(hostName, "vision")==0) return 16;
  else if (strCompare(hostName, "visky") ==0) return 4;
  else if (strCompare(hostName, "visit") ==0) return 2;
  else return 1;
#endif
}

//----------------------------------------------------------------------------
/** Returns RGBA texture values for a hue/saturation color chooser.
 Texture values are returned in data as 4 bytes per texel, bottom to top, 
 ordered RGBARGBARGBA...
 @param width,height   width and height of texture [pixels]
 @param brightness     brightness of color values [0..1]
 @param data           pointer to allocated memory space providing width * height * 4 bytes
*/
void vvToolshed::makeColorBoardTexture(int width, int height, float brightness, 
    uchar* data)
{
  float h, s, v;   // hue, saturation, value
  float r, g, b;   // RGB
  float nx, ny;    // x and y normalized to range [0..1]
	float dx, dy; 	 // distance from center
  int   i = 0;     // index of current texel element
  int   x, y;      // current texel position

  for (y=0; y<height; ++y)
    for (x=0; x<width; ++x)
    {
      nx = (float)x / (float)(width-1);
      ny = (float)y / (float)(height-1);
			dx = 2.0f * nx - 1.0f;
			dy = 2.0f * ny - 1.0f;
			if ( (dx * dx + dy * dy) > 1.0f) 
			{
				// Outer area is black:
        data[i++] = (uchar)0;    // red
        data[i++] = (uchar)0;    // green 
        data[i++] = (uchar)0;    // blue
        data[i++] = (uchar)255;  // alpha
			}
			else 
			{
				v = brightness; 	// circle area has requested brightness
     	  convertXY2HS(nx, ny, &h, &s);
        vvToolshed::HSBtoRGB(h, s, v, &r, &g, &b);
        data[i++] = (uchar)(r * 255.0f);  // red
        data[i++] = (uchar)(g * 255.0f);  // green 
        data[i++] = (uchar)(b * 255.0f);  // blue
        data[i++] = (uchar)255;           // alpha
      }
		}
}

//----------------------------------------------------------------------------
/** The given x|y coordinates of the mouse are translated to hue and 
    saturation values. 
 Mouse coordinate 0|0 is bottom left, 1|1 is top right. 
 Hue and saturation values are in range [0..1].
*/
void vvToolshed::convertXY2HS(float x, float y, float* hue, float* saturation)
{ 
  float dx, dy; 	// distance from center of x/y area

  // Determine hue:
  dx = x - 0.5f;
  dy = y - 0.5f;

  if (dx==0.0f)
  {
    if (dy>=0.0f) *hue = 0.0f;
    else *hue = 180.0f;
  }
  else if (dy==0.0f)
  {
    if (dx>0.0f) *hue = 90.0f;
    else *hue = 270.0f;
  }
  else
  {  
    if      (dx>0.0f && dy>0.0f) *hue = atanf(dx / dy);
    else if (dx>0.0f && dy<0.0f) *hue = TS_PI * 0.5f + atanf(-dy / dx);
    else if (dx<0.0f && dy<0.0f) *hue = TS_PI + atanf(-dx / -dy);
    else                         *hue = TS_PI * 1.5f + atanf(dy / -dx);
  }
  *hue /= (2.0f * TS_PI);
  *hue = ts_clamp(*hue, 0.0f, 1.0f);

  // Determine saturation:
	dx *= 2.0f;
	dy *= 2.0f;
  *saturation = sqrtf(dx * dx + dy * dy);
  *saturation = ts_clamp(*saturation, 0.0f, 1.0f);
}

//----------------------------------------------------------------------------
/** The given hue and saturation values are converted to mouse x|y coordinates. 
 Hue and saturation values are in range [0..1].
 Mouse coordinate 0|0 is bottom left, 1|1 is top right. 
*/
void vvToolshed::convertHS2XY(float hue, float saturation, float* x, float* y)
{
  float angle;    // angle of point xy position within color circle
  float dx, dy;   // point position relative to circle midpoint

  angle = hue * 2.0f * TS_PI;
  dx = 0.5f * saturation * sinf(angle);
  dy = 0.5f * saturation * cosf(angle);
  *x = dx + 0.5f;
  *y = dy + 0.5f;
  *x = ts_clamp(*x, 0.0f, 1.0f);
  *y = ts_clamp(*y, 0.0f, 1.0f);
}

//----------------------------------------------------------------------------
/** Read an unsigned char value from a file.
*/
uchar vvToolshed::read8(FILE* src)
{
  uchar val;
  
  fread(&val, 1, 1, src);
  return val;
}

//----------------------------------------------------------------------------
/** Write an unsigned char value to a file.
  @return number of bytes written
*/
int vvToolshed::write8(FILE* dst, uchar val)
{
  fwrite(&val, 1, 1, dst);
  return 1;
}

//----------------------------------------------------------------------------
/** Read a big endian unsigned short value system independently from a file
  (most significant byte first).
*/
ushort vvToolshed::read16BE(FILE* src)
{
  uchar buf[2];
  int val;
  
  fread(buf, 2, 1, src);
  val = (int)buf[0] * (int)256 + (int)buf[1];
  return (ushort)val;
}

//----------------------------------------------------------------------------
/** Read a little endian unsigned short value system independently from a file
  (least significant byte first).
*/
ushort vvToolshed::read16LE(FILE* src)
{
  uchar buf[2];
  int val;
  
  fread(buf, 2, 1, src);
  val = (int)buf[0] + (int)buf[1] * (int)256;
  return (ushort)val;
}

//----------------------------------------------------------------------------
/** Write a big endian unsigned short value system independently to a file
  (most significant byte first).
  @return number of bytes written
*/
int vvToolshed::write16BE(FILE* fp, ushort val)
{
  uchar buf[2];

  buf[0] = (uchar)(val >> 8); 
  buf[1] = (uchar)(val & 0xFF);
  fwrite(buf, 2, 1, fp);
  return 2;
}

//----------------------------------------------------------------------------
/** Write a little endian unsigned short value system independently to a file
  (least significant byte first).
  @return number of bytes written
*/
int vvToolshed::write16LE(FILE* fp, ushort val)
{
  uchar buf[2];

  buf[0] = (uchar)(val & 0xFF);
  buf[1] = (uchar)(val >> 8); 
  fwrite(buf, 2, 1, fp);
  return 2;
}

//----------------------------------------------------------------------------
/** Read a big endian unsigned long value system independently from a file.
 Read four bytes in a row in unix-style (most significant byte first).
*/
ulong vvToolshed::read32BE(FILE* src)
{
  uchar buf[4];
  ulong val;
  
  fread(buf, 4, 1, src);
  val = (ulong)buf[0] * (ulong)16777216 + (ulong)buf[1] * (ulong)65536 + 
        (ulong)buf[2] * (ulong)256 + (ulong)buf[3];
  return (ulong)val;
}

//----------------------------------------------------------------------------
/** Read a little endian unsigned long value system independently from a file.
 Read four bytes in a row in unix-style (least significant byte first).
*/
ulong vvToolshed::read32LE(FILE* src)
{
  uchar buf[4];
  ulong val;
  
  fread(buf, 4, 1, src);
  val = (ulong)buf[3] * (ulong)16777216 + (ulong)buf[2] * (ulong)65536 + 
        (ulong)buf[1] * (ulong)256 + (ulong)buf[0];
  return (ulong)val;
}

//----------------------------------------------------------------------------
/** Write a big endian unsigned long value system independently to a file. 
  Write four bytes in a row in unix-style (most significant byte first).
  @return number of bytes written
*/
int vvToolshed::write32BE(FILE* fp, ulong val)
{
  uchar buf[4];

  buf[0] = (uchar)(val  >> 24); 
  buf[1] = (uchar)((val >> 16) & 0xFF);
  buf[2] = (uchar)((val >> 8)  & 0xFF);
  buf[3] = (uchar)(val & 0xFF);
  fwrite(buf, 4, 1, fp);
  return 4;
}

//----------------------------------------------------------------------------
/** Write a little endian unsigned long value system independently to a file. 
  Write four bytes in a row in unix-style (least significant byte first).
  @return number of bytes written
*/
int vvToolshed::write32LE(FILE* fp, ulong val)
{
  uchar buf[4];

  buf[0] = (uchar)(val & 0xFF);
  buf[1] = (uchar)((val >> 8)  & 0xFF);
  buf[2] = (uchar)((val >> 16) & 0xFF);
  buf[3] = (uchar)(val  >> 24); 
  fwrite(buf, 4, 1, fp);
  return 4;
}

//----------------------------------------------------------------------------
/** Read a 32 bit float value system independently from a file.
 Read four bytes in a row in unix-style (most significant byte first).
*/
float vvToolshed::readFloat(FILE* src)
{
  uchar buf[4];
  uchar tmp;
  
  fread(buf, 4, 1, src);
  if (getEndianness()==LITTLE_END)  
  { 
    // Reverse byte order:
    tmp = buf[0]; buf[0] = buf[3]; buf[3] = tmp;
    tmp = buf[1]; buf[1] = buf[2]; buf[2] = tmp;
  }
  return *((float*)buf);
}

//----------------------------------------------------------------------------
/** Write a 32 bit float value system independently to a file.
  Write four bytes in a row in unix-style (most significant byte first).
  @return number of bytes written
*/
int vvToolshed::writeFloat(FILE* fp, float val)
{
  uchar* buf;
  uchar tmp;
  
  if (getEndianness()==LITTLE_END)  
  { 
    // Reverse byte order:
    buf = (uchar*)&val;
    tmp = buf[0]; buf[0] = buf[3]; buf[3] = tmp;
    tmp = buf[1]; buf[1] = buf[2]; buf[2] = tmp;
  }

  fwrite(&val, 4, 1, fp);
  return 4;
}

//----------------------------------------------------------------------------
/** Read an unsigned char value from a buffer.
*/
uchar vvToolshed::read8(uchar* src)
{
  return *src;
}

//----------------------------------------------------------------------------
/** Write an unsigned char value to a buffer.
  @return number of bytes written
*/
int vvToolshed::write8(uchar* src, uchar val)
{
  *src = val;
  return sizeof(uchar);
}

//----------------------------------------------------------------------------
/** Read a big endian unsigned short value system independently from a buffer
  (most significant byte first).
*/
ushort vvToolshed::read16BE(uchar* src)
{
  int val;
  
  val = (int)src[0] * (int)256 + (int)src[1];
  return (ushort)val;
}

//----------------------------------------------------------------------------
/** Read a little endian unsigned short value system independently from a buffer
  (least significant byte first).
*/
ushort vvToolshed::read16LE(uchar* src)
{
  int val;
  
  val = (int)src[0] + (int)src[1] * (int)256;
  return (ushort)val;
}

//----------------------------------------------------------------------------
/** Write a big endian unsigned short value system independently to a buffer
  (most significant byte first).
  @param buf pointer to 2 bytes of _allocated_ memory
  @return number of bytes written
*/
int vvToolshed::write16BE(uchar* buf, ushort val)
{
  buf[0] = (uchar)(val >> 8); 
  buf[1] = (uchar)(val & 0xFF);
  return sizeof(ushort);
}

//----------------------------------------------------------------------------
/** Write a little endian unsigned short value system independently to a buffer
  (least significant byte first). 
  @param buf pointer to 2 bytes of _allocated_ memory
  @return number of bytes written
*/
int vvToolshed::write16LE(uchar* buf, ushort val)
{
  buf[0] = (uchar)(val & 0xFF);
  buf[1] = (uchar)(val >> 8); 
  return sizeof(ushort);
}

//----------------------------------------------------------------------------
/** Read a big endian unsigned long value system independently from a buffer.
 Read four bytes in a row in unix-style (most significant byte first).
*/
ulong vvToolshed::read32BE(uchar* buf)
{
  ulong val;
  
  val = (ulong)buf[0] * (ulong)16777216 + (ulong)buf[1] * (ulong)65536 + 
        (ulong)buf[2] * (ulong)256 + (ulong)buf[3];
  return (ulong)val;
}

//----------------------------------------------------------------------------
/** Read a little endian unsigned long value system independently from a buffer.
 Read four bytes in a row in unix-style (least significant byte first).
*/
ulong vvToolshed::read32LE(uchar* buf)
{
  ulong val;
  
  val = (ulong)buf[3] * (ulong)16777216 + (ulong)buf[2] * (ulong)65536 + 
        (ulong)buf[1] * (ulong)256 + (ulong)buf[0];
  return (ulong)val;
}

//----------------------------------------------------------------------------
/** Write a big endian unsigned long value system independently to a buffer. 
  Write four bytes in a row in unix-style (most significant byte first).
  @return number of bytes written
*/
int vvToolshed::write32BE(uchar* buf, ulong val)
{
  buf[0] = (uchar)(val  >> 24); 
  buf[1] = (uchar)((val >> 16) & 0xFF);
  buf[2] = (uchar)((val >> 8)  & 0xFF);
  buf[3] = (uchar)(val & 0xFF);
  return sizeof(ulong);
}

//----------------------------------------------------------------------------
/** Write a little endian unsigned long value system independently to a buffer. 
  Write four bytes in a row in unix-style (least significant byte first).
  @return number of bytes written
*/
int vvToolshed::write32LE(uchar* buf, ulong val)
{
  buf[0] = (uchar)(val & 0xFF);
  buf[1] = (uchar)((val >> 8)  & 0xFF);
  buf[2] = (uchar)((val >> 16) & 0xFF);
  buf[3] = (uchar)(val  >> 24); 
  return sizeof(ulong);
}

//----------------------------------------------------------------------------
/** Read a 32 bit float value system independently from a buffer.
 Read four bytes in a row in unix-style (most significant byte first).
*/
float vvToolshed::readFloat(uchar* buf)
{
  float  fval;
  uchar* ptr;
  uchar  tmp;
  
  assert(sizeof(float)==4);
  memcpy(&fval, buf, 4);
  if (getEndianness()==LITTLE_END)  
  { 
    // Reverse byte order:
    ptr = (uchar*)&fval;
    tmp = ptr[0]; ptr[0] = ptr[3]; ptr[3] = tmp;
    tmp = ptr[1]; ptr[1] = ptr[2]; ptr[2] = tmp;
  }
  return fval;
}

//----------------------------------------------------------------------------
/** Write a 32 bit float value system independently to a buffer.
  Write four bytes in a row in unix-style (most significant byte first).
  @return number of bytes written
*/
int vvToolshed::writeFloat(uchar* buf, float val)
{
  uchar tmp;
  
  assert(sizeof(float)==4);
  memcpy(buf, &val, 4);
  if (getEndianness()==LITTLE_END)  
  { 
    // Reverse byte order:
    tmp = buf[0]; buf[0] = buf[3]; buf[3] = tmp;
    tmp = buf[1]; buf[1] = buf[2]; buf[2] = tmp;
  }
  return sizeof(float);
}

//----------------------------------------------------------------------------
/** Make a float array system independent:
  convert each four byte float to unix-style (most significant byte first).
  @param numValues number of float values in array
  @param array     pointer to float array (size of array must be 4*numValues!)
*/
void vvToolshed::makeArraySystemIndependent(int numValues, float* array)
{
  uchar* buf;   // array pointer in uchar format
  int i;
  uchar tmp;    // temporary byte value from float array, needed for swapping

  assert(sizeof(float) == 4);
  if (getEndianness()==BIG_END)  return;    // nothing needs to be done

  buf = (uchar*)array;
  for (i=0; i<numValues; ++i)
  {
    // Reverse byte order:
    tmp = buf[0]; buf[0] = buf[3]; buf[3] = tmp;
    tmp = buf[1]; buf[1] = buf[2]; buf[2] = tmp;
    buf += 4;
  }
}

//----------------------------------------------------------------------------
/** Make a system independent float array system dependent:
  convert each four byte float value back to system style.
  @param numValues number of float values in array
  @param array     pointer to float array (size of array must be 4*numValues!)
*/
void vvToolshed::makeArraySystemDependent(int numValues, float* array)
{
  // Swapping bytes is the same as above, therefore use the same code:
  makeArraySystemIndependent(numValues, array);
}

//----------------------------------------------------------------------------
/** Returns the current system's endianness.
*/
vvToolshed::EndianType vvToolshed::getEndianness()
{
  float one = 1.0f;   // memory representation of 1.0 on big endian machines: 3F 80 00 00
  uchar* ptr;

  ptr = (uchar*)&one;
  if (*ptr == 0x3f) return BIG_END;
  else              return LITTLE_END;
}

//----------------------------------------------------------------------------
/** Suspend process for a specific time. If milliseconds are not available
  on a specific system type, seconds are used (e.g., on Cray systems).
  @param msec suspension time [milliseconds]
*/
void vvToolshed::sleep(int msec)
{
#ifdef WIN32  
  Sleep(msec);
#elif CRAY
  sleep(msec / 1000);
#else
  usleep(msec * 1000);
#endif
}

//----------------------------------------------------------------------------
/// Main function for standalone test mode.
#ifdef VV_STANDALONE
int main(int, char**)
{
#ifdef WIN32
  char* pathname={"c:\\user\\testfile.dat"};
#else
  char* pathname={"/usr/local/testfile.dat"};
#endif
  char  teststring[256];

  cout << "ts_max(2,9)  = " << ts_max(2,9)  << endl;
  cout << "ts_min(2,9)  = " << ts_min(2,9)  << endl;
  cout << "ts_abs(-7)   = " << ts_abs(-7)   << endl;
  cout << "ts_sgn(-9.1) = " << ts_sgn(-9.1) << endl;
  cout << "ts_zsgn(0.0) = " << ts_zsgn(0.0) << endl;
  cout << "ts_zsgn(-2)  = " << ts_zsgn(-2)  << endl;
  cout << "ts_clamp(1.2, 1.0, 2.0)  = " << ts_clamp(1.2f, 1.0f, 2.0f)  << endl;
  cout << "ts_clamp(-0.5, 1.0, 2.0)  = " << ts_clamp(0.5f, 1.0f, 2.0f)  << endl;
  cout << "ts_clamp(2.1, 1.0, 2.0)  = " << ts_clamp(2.1f, 1.0f, 2.0f)  << endl;

  cout << "isSuffix(" << pathname << "), 'Dat' = ";
  if (vvToolshed::isSuffix(pathname, "Dat") == true)
    cout << "true" << endl;
  else
    cout << "false" << endl;

  cout << "isSuffix(" << pathname << "), 'data' = ";
  if (vvToolshed::isSuffix(pathname, "data") == true)
    cout << "true" << endl;
  else
    cout << "false" << endl;

  vvToolshed::extractFilename(teststring, pathname);
  cout << "extractFilename(" << pathname << ") = " << teststring << endl;

  vvToolshed::extractDirname(teststring, pathname);
  cout << "extractDirname(" << pathname << ") = " << teststring << endl;

  vvToolshed::extractExtension(teststring, pathname);
  cout << "extractExtension(" << pathname << ") = " << teststring << endl;

  vvToolshed::extractBasename(teststring, pathname);
  cout << "extractBasename(" << pathname << ") = " << teststring << endl;

  cout << "getTextureSize(84) = " << vvToolshed::getTextureSize(84) << endl;

  char* testData = {"ABABACACACABABABACABCD"};
  char encoded[100];
  char decoded[100];
  int len;
  int bpc = 2;
  cout << "Unencoded: " << testData << endl;
  len = vvToolshed::encodeRLE((uchar*)encoded, (uchar*)testData, strlen(testData), bpc, 100);
  len = vvToolshed::decodeRLE((uchar*)decoded, (uchar*)encoded, len, bpc, 100);
  decoded[len] = '\0';
  cout << "Decoded:   " << decoded << endl;

  return 1;
}

#endif

//============================================================================
// End of File
//============================================================================
