//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifdef WIN32
  #include <float.h>
#endif
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "vvdebugmsg.h"
#include "vvvecmath.h"
#include "vvtransfunc.h"

//============================================================================
// Class vvTransFunc
//============================================================================

//----------------------------------------------------------------------------
/// Constructor
vvTransFunc::vvTransFunc()
{
  nextBufferEntry = 0;
  bufferUsed      = 0;
  setDiscreteColors(0);
}

//----------------------------------------------------------------------------
/// Destructor
vvTransFunc::~vvTransFunc()
{
}

//----------------------------------------------------------------------------
/** Set default color values in the transfer function.
 The previous color pins are deleted, alpha pins are not affected.
*/
void vvTransFunc::setDefaultColors(int index)
{
  vvDebugMsg::msg(2, "vvTransFunc::setDefaultColors()");
  
  pins.clearColor();
  pins.setColorModel(vvPinList::RGB_MODEL);
  switch (index)
  {
    case 0:       // bright colors
    default:
      // Set RGBA table to bright colors (range: blue->green->red):
      pins.add(new vvPin(vvPin::COLOR, 0.0f,  0.0f, 0.0f, 1.0f));
      pins.add(new vvPin(vvPin::COLOR, 0.33f, 0.0f, 1.0f, 1.0f));
      pins.add(new vvPin(vvPin::COLOR, 0.67f, 1.0f, 1.0f, 0.0f));
      pins.add(new vvPin(vvPin::COLOR, 1.0f,  1.0f, 0.0f, 0.0f));
      break;

    case 1:       // hue gradient
      // Set RGBA table to maximum intensity and value HSB colors:
      pins.add(new vvPin(vvPin::COLOR, 0.0f, 1.0f, 0.0f, 0.0f));
      pins.add(new vvPin(vvPin::COLOR, 0.2f, 1.0f, 1.0f, 0.0f));
      pins.add(new vvPin(vvPin::COLOR, 0.4f, 0.0f, 1.0f, 0.0f));
      pins.add(new vvPin(vvPin::COLOR, 0.6f, 0.0f, 1.0f, 1.0f));
      pins.add(new vvPin(vvPin::COLOR, 0.8f, 0.0f, 0.0f, 1.0f));
      pins.add(new vvPin(vvPin::COLOR, 1.0f, 1.0f, 0.0f, 1.0f));
      break;

    case 2:       // grayscale ramp
      // Set RGBA table to grayscale ramp (range: black->white). 
      // Set only a black pin to 0.0 because 1.0 is automatically
      // white:
      pins.add(new vvPin(vvPin::COLOR, 0.0f, 0.0f, 0.0f, 0.0f));
      break;

    case 3:       // white
      // Set RGBA table to all white values by not adding
      // any pin to the default.
      break;
  }
}

//----------------------------------------------------------------------------
/// Returns the number of default color schemes.
int vvTransFunc::getNumDefaultColors()
{
	return 4;
}

//----------------------------------------------------------------------------
/** Set default alpha values in the transfer function.
 The previous alpha pins are deleted, color pins are not affected.
*/
void vvTransFunc::setDefaultAlpha(int index)
{
  vvDebugMsg::msg(2, "vvTransFunc::setDefaultAlpha()");

  pins.clearAlpha();
  switch (index)
  {
    case 0:     // ascending (0->1)
    default:
      pins.add(new vvPin(vvPin::ALPHA_RAMP, 0.5f,  1.0f, 1.0f));
      break;
    case 1:     // descending (1->0)
      pins.add(new vvPin(vvPin::ALPHA_RAMP, 0.5f, -1.0f, 1.0f));
      break;
    case 2:     // opaque (all 1)
      pins.add(new vvPin(vvPin::ALPHA_HAT, 0.5f, 1.0f, 1.0f, 1.0f));
      break;
  }
}

//----------------------------------------------------------------------------
/// Returns the number of default alpha schemes.
int vvTransFunc::getNumDefaultAlpha()
{
	return 3;
}

//----------------------------------------------------------------------------
/// Returns the number of non-empty preset pin lists.
int vvTransFunc::getNumUsedPresets()
{
	int i, used = 0;
  for (i=0; i<getNumPresets(); ++i)
    if (!presets[i].isEmpty()) ++used;
  return used;
}

//----------------------------------------------------------------------------
/** Returns true if a default color scheme is set.
  @param index number of default color scheme
*/
bool vvTransFunc::isDefaultColors(int index)
{
  vvTransFunc test1, test2;

  test1.setDefaultColors(index);
  test2.copy(this);
  test2.pins.clearAlpha();
	if (test1.pins.isEqual(&test2.pins))
    return true;
  else return false;
}

//----------------------------------------------------------------------------
/** Returns true if a default alpha scheme is set.
  @param index number of default alpha scheme
*/
bool vvTransFunc::isDefaultAlpha(int index)
{
  vvTransFunc test1, test2;

  test1.setDefaultAlpha(index);
  test2.copy(this);
  test2.pins.clearColor();
	if (test1.pins.isEqual(&test2.pins))
    return true;
  else return false;
}

//----------------------------------------------------------------------------
/** Returns true if no transfer function is set (no pins in list).
*/
bool vvTransFunc::isEmpty()
{
  bool empty = true;

  if (!pins.isEmpty()) empty = false;
  for (int i=0; i<getNumPresets(); ++i)
    if (!presets[i].isEmpty()) empty = false;
  return empty;
}

//----------------------------------------------------------------------------
/** Returns RGBA texture values for a color preview bar.
@param num    number of pixels to generate (>0)
@param colors pointer to _allocated_ memory providing space for num x 2 x 4 bytes.
              Byte quadruple 0 will be considered to correspond with scalar 
              value 0.0, quadruple num-1 will be considered to correspond with 
              scalar value 1.0. The resulting RGBA values are stored in the
              following order: RGBARGBARGBA...
*/
void vvTransFunc::makeColorBar(int num, uchar* colors)
{
  float* components;  // component values
  int c, x;

  vvDebugMsg::msg(3, "vvTransFunc::makeColorBar()");

  assert(colors != NULL);
  if (num<1) return;    // already done

  // Compute color components:
  components = new float[num * 4 * 2];    // four bytes per pixel, 2 lines
  assert(components != NULL);
  pins.setColorModel(vvPinList::RGB_MODEL);
  pins.computeFunction(num, vvPinList::RGBA, &components[0]);

  // Rearrange values and convert to uchar:
  for (x=0; x<num; ++x)
	{
    for (c=0; c<3; ++c)
    {
      colors[x * 4 + c] = colors[x * 4 + c + num * 4] = 
        uchar(components[c * num + x] * 255.0f);
    }
		colors[x * 4 + 3] = (uchar)255;
		colors[x * 4 + 3 + num * 4] = uchar(components[3 * num + x] * 255.0f);
	}
  delete[] components;
}

//----------------------------------------------------------------------------
/** Create a look-up table of 8-bit integer values from current transfer
  function.
  @param entries number of LUT entries (typically 256 or 4096, depending on bpv)
  @param lut     _allocated_ space with space for entries*4 bytes
*/
void vvTransFunc::make8bitLUT(int entries, uchar* lut)
{
  float* rgba;    // temporary LUT in floating point format
  int i, c;       

  vvDebugMsg::msg(1, "vvTransFunc::make8bitLUT()");

  rgba = new float[4 * entries];

  // Generate arrays from pins:
  pins.setColorModel(vvPinList::RGB_MODEL);
  pins.computeFunction(entries, vvPinList::RGBA, rgba);

  // Copy RGBA values to internal array:
  for (i=0; i<entries; ++i)
		for (c=0; c<4; ++c)
    {
      *lut = uchar(rgba[c * entries + i] * 255.0f);
      ++lut;
    }

  delete[] rgba;
}

//----------------------------------------------------------------------------
/** Create RGBA texture values for a bitmap display of the transfer function.
@param width,height size of texture [pixels]
@param texels       pointer to _allocated_ memory providing space for width x height x 4 bytes.
                    The resulting RGBA values are stored in the following order: RGBARGBARGBA...
                    The first byte contains the bottom left pixel of the texture.
@param colors       true = include color bar in texture, false = only display opacity information
                    If colors is true, the top of the resulting texture displays an opaque color bar,
                    below is a translucent color bar. Both bars account for 1/10th of the vertical
                    image size.
*/
void vvTransFunc::makeFunctionTexture(int width, int height, uchar* texels, bool colors)
{
  const uchar ALPHA_COLOR[4] = {170,170,170,255};   // alpha function color
  const uchar BG_COLOR[4] = {0,0,0,255};            // background color
  const int BPT = 4;                  // bytes per texel
  const float COLOR_FRACTION = 0.1f;  // fraction of each color bar from total vertical image size
  uchar* colorBar;                    // pointer to color bar values which are returned from makeColorBar()
  uchar* colorBarOffset;              // pointer to first color bar texel (bottom left corner of color bar area)
  float xval;                         // floating point x value [0..1]
  int x,y;                            // image coordinates [pixels]
  int colorBarHeight;                 // height of each color bar [texels]
  int barHeight;                      // height of an alpha bar [pixels]
  int i;

  vvDebugMsg::msg(3, "vvTransFunc::makeFunctionTexture()");
  assert(texels != NULL && width>0 && height>0);

  // Clear destination texture:
  for (i=0; i<width*height; ++i)
    memcpy(texels + i * BPT, BG_COLOR, BPT);

  // Generate color bar:
  if (colors)
  {
    colorBar = new uchar[width * 2 * BPT];
    makeColorBar(width, colorBar);
    colorBarHeight = int(float(height) * COLOR_FRACTION);
    colorBarOffset = texels + (height - 2 * colorBarHeight) * width * BPT;
    for (i=0; i<colorBarHeight; ++i)    // replicate bottom line of color bar
      memcpy(colorBarOffset + i * width * BPT, colorBar, width * BPT);
    colorBarOffset += colorBarHeight * width * BPT;
    for (i=0; i<colorBarHeight; ++i)    // replicate bottom line of color bar
      memcpy(colorBarOffset + i * width * BPT, colorBar + width * BPT, width * BPT);
    delete colorBar;
  }
  else colorBarHeight = 0;

  // Generate opacity function:
  for (x=0; x<width; ++x)
  {
    xval = float(x) / float(width-1);
    barHeight = int(float(height - 2 * colorBarHeight) * pins.computeY(xval, vvPinList::ALPHA));
    for (y=0; y<barHeight; ++y)
      memcpy(texels + BPT * (x + y * width), ALPHA_COLOR, BPT);
  }
}

//----------------------------------------------------------------------------
/** Copy active transfer function.
 The source transfer function instance may be deleted after copying.
 @param src pointer to source transfer function instance
*/
void vvTransFunc::copy(vvTransFunc* src)
{
  vvDebugMsg::msg(2, "vvTransFunc::copy()");
  pins.copy(&src->pins);
}

//----------------------------------------------------------------------------
/** Copy preset transfer functions.
 The source transfer function instance may be deleted after copying.
 @param src pointer to source transfer function instance
*/
void vvTransFunc::copyPresets(vvTransFunc* src)
{
  int i;

  vvDebugMsg::msg(2, "vvTransFunc::copyPresets()");
  for (i=0; i<NUM_PRESETS; ++i)
    presets[i].copy(&src->presets[i]);
}

//----------------------------------------------------------------------------
/** Stores the current transfer function to a preset slot.
 @param index slot to save function to [0..getNumPresets()]
*/
void vvTransFunc::storePreset(int index)
{
  vvDebugMsg::msg(2, "vvTransFunc::storePreset()");
  if (index<0 || index>=getNumPresets()) return;
  presets[index].copy(&pins);
}

//----------------------------------------------------------------------------
/** Restores the transfer function from a preset slot.
  Range of 'index': [0..getNumPresets()]
*/
void vvTransFunc::restorePreset(int index)
{
  vvDebugMsg::msg(2, "vvTransFunc::restorePreset()");
	pins.copy(&presets[index]);
}

//----------------------------------------------------------------------------
/// Returns the number of available preset slots.
int vvTransFunc::getNumPresets()
{
  vvDebugMsg::msg(2, "vvTransFunc::getNumPresets()");
  return NUM_PRESETS;
}

//----------------------------------------------------------------------------
/// Delete all pins from all transfer functions.
void vvTransFunc::clear()
{
  int i;

  vvDebugMsg::msg(2, "vvTransFunc::clear()");
  pins.clear();
  for (i=0; i<NUM_PRESETS; ++i)
    presets[i].clear();
}

//----------------------------------------------------------------------------
/// Store the current pin list in the undo ring buffer.
void vvTransFunc::putUndoBuffer()
{
  int prevBufferEntry;  	  // previous buffer entry
  
  vvDebugMsg::msg(1, "vvTransFunc::putUndoBuffer()");

  // Check if this TF is already in ring buffer:
  if (nextBufferEntry > 0) prevBufferEntry = nextBufferEntry - 1;
  else prevBufferEntry = BUFFER_SIZE - 1;
  if (bufferUsed>0 && pins.isEqual(&buffer[prevBufferEntry])) return;
  
  // Add TF to ring buffer:
  buffer[nextBufferEntry].copy(&pins);
  if (bufferUsed < BUFFER_SIZE) ++bufferUsed;
  if (nextBufferEntry < BUFFER_SIZE-1) ++nextBufferEntry;
  else nextBufferEntry = 0;
}

//----------------------------------------------------------------------------
/** Restore the latest element from the undo ring buffer to the current pin list.
  If the ring buffer is empty, nothing happens.
*/
void vvTransFunc::getUndoBuffer()
{
  int bufferEntry;

  vvDebugMsg::msg(1, "vvTransFunc::getUndoBuffer()");

  if (bufferUsed==0) return;    // ring buffer is empty
  if (nextBufferEntry > 0) bufferEntry = nextBufferEntry - 1;
  else bufferEntry = BUFFER_SIZE - 1;
  pins.copy(&buffer[bufferEntry]);
  nextBufferEntry = bufferEntry;
  --bufferUsed;
}

//----------------------------------------------------------------------------
/// Clear the undo ring buffer.
void vvTransFunc::clearUndoBuffer()
{
  vvDebugMsg::msg(1, "vvTransFunc::clearUndoBuffer()");

  bufferUsed      = 0;
  nextBufferEntry = 0;
}

//----------------------------------------------------------------------------
/** Set the number of discrete colors to use for color interpolation.
  @param numColors number of discrete colors (use 0 for smooth colors)
*/
void vvTransFunc::setDiscreteColors(int numColors)
{
  vvDebugMsg::msg(1, "vvTransFunc::setDiscreteColors()");

  if (numColors<0) numColors = 0;
  discreteColors = numColors;
  pins.setDiscreteColors(discreteColors);
  for (int i=0; i<NUM_PRESETS; ++i)
    presets[i].setDiscreteColors(discreteColors);  
}

//----------------------------------------------------------------------------
/** @return the number of discrete colors used for color interpolation.
            0 means smooth colors.
*/
int vvTransFunc::getDiscreteColors()
{
  vvDebugMsg::msg(1, "vvTransFunc::getDiscreteColors()");
  return discreteColors;
}

//============================================================================
// Class vvPin
//============================================================================

int vvPin::idCounter = 0;   // ID counter starts at 0 and counts up

//----------------------------------------------------------------------------
/** Constructor
 @param t     pin type (see enum PinType)
 @param xpos  x position of pin [0..1]
 @param a,b,c value components<BR>
*/
vvPin::vvPin(PinType t, float xpos, float a, float b, float c)
{
  id = idCounter;
  ++idCounter;
  type = t;
  x    = ts_clamp(xpos, 0.0f, 1.0f);
  if (t==vvPin::COLOR)
    v[0] = ts_clamp(a, 0.0f, 1.0f);
  else if (t==vvPin::ALPHA_BLANK)
    v[0] = ts_clamp(a, 0.0f, VV_FLT_MAX);
  else if (t==vvPin::ALPHA_HAT)
    v[0] = (float)fabs(a);
  else
    v[0] = a;
  v[1] = ts_clamp(b, 0.0f, 1.0f);
  v[2] = ts_clamp(c, 0.0f, 1.0f);
  next = NULL;
}

//----------------------------------------------------------------------------
/// Copy constructor
vvPin::vvPin(vvPin* pin)
{
  id   = pin->id;
  type = pin->type;
  x    = pin->x;
  for (int c=0; c<3; ++c)
    v[c] = pin->v[c];
  next = NULL;
}

//----------------------------------------------------------------------------
/// @return pin ID
int vvPin::getID()
{
  return id;
}

//============================================================================
// Class vvPinList and descendents
//============================================================================

//----------------------------------------------------------------------------
/// Constructor
vvPinList::vvPinList()
{
  root = NULL;
  colorModel = RGB_MODEL;
  discreteColors = 0;
}

//----------------------------------------------------------------------------
/** Destructor.
 Deletes all pins from the list.
*/
vvPinList::~vvPinList()
{
  clear();
}

//----------------------------------------------------------------------------
/** Add a vvPinList to the list.
  Insertion is done at the correct position sorted by x values.
  @param newPin  pointer to new pin (caller must _not_ delete this object)
*/
void vvPinList::add(vvPin* newPin)
{
  vvPin* pin;
  vvPin* nextPin;

  vvDebugMsg::msg(3, "vvPinList::add()");

  if (root==NULL)     // list is empty
  {
    root = newPin;
  }
  else if (newPin->x < root->x) // insert pin as new root
  {
    pin = root;
    root = newPin;
    root->next = pin;
  }
  else      // regular insertion
  {
    pin = root;
    for (;;)
    {
      if (pin->next==NULL) // append as last pin?
      {
        nextPin = newPin;
        pin->next = nextPin;
        newPin->next = NULL;
        return;
      }
      else if (newPin->x <= pin->next->x)  // insert here
      {
        nextPin = pin->next;
        pin->next = newPin;
        pin->next->next = nextPin;
        return;
      }
      else 
        pin = pin->next;
    }
  }
}

//----------------------------------------------------------------------------
/** Remove a pin from the list.
 @param p  pin index (0 for root pin)
*/
void vvPinList::remove(int p)
{
  vvPin* pin;
  vvPin* prev;    // previous pin
  int i;

  vvDebugMsg::msg(3, "vvPinList::remove()");

  if (root==NULL || p<0) return;     // no pins in list or invalid index
  pin = prev = root;
  if (p==0)   // remove root pin?
  { 
    root = root->next;
    delete pin;
  }
  else
  {
    for (i=0; i<p; ++i)
    {
      prev = pin;
      pin = pin->next;
      if (pin == NULL) return;    // return if pin number doesn't exist
    }    
    prev->next = pin->next;
    delete pin;
  }
}
//----------------------------------------------------------------------------
/** Get vvPin pointer with index p.
  @return NULL for invalid indices
*/
vvPin* vvPinList::getPin(int p)
{
  vvPin* pin;
  int i;

  vvDebugMsg::msg(2, "vvPinList::remove()");

  if (root==NULL || p<0) return NULL;     // no pins in list or invalid index
  pin = root;
  for (i=0; i<p; ++i)
  {
    pin = pin->next;
    if (pin == NULL) return NULL;    // return if pin number doesn't exist
  }    
  return pin;
}

//----------------------------------------------------------------------------
/** Remove a pin from the list.
  @param p  pointer to pin
*/
void vvPinList::remove(vvPin* p)
{
  vvPin* pin;
  vvPin* nextPin;
  vvPin* prevPin;

  vvDebugMsg::msg(2, "vvPinList::remove()");
  pin = prevPin = root;
  while (pin != NULL)
  {
    nextPin = pin->next;    
    if (pin==p)
    {
      if (pin==root) root = nextPin;
      else prevPin->next = nextPin;
      delete pin;
      return;
    }
    prevPin = pin;
    pin = nextPin;
  }
}

//----------------------------------------------------------------------------
/** Move a pin to another x position.
 @param p       pointer to pin
 @param newPos  new x position
*/
void vvPinList::move(vvPin* p, float newPos)
{
  vvPin* pin;
  vvPin* nextPin;
  vvPin* prevPin;

  vvDebugMsg::msg(3, "vvPinList::move()");
  pin = prevPin = root;
  while (pin != NULL)
  {
    nextPin = pin->next;    
    if (pin==p)
    {
      // Remove pin from current list position:
      if (pin==root) root = nextPin;
      else prevPin->next = nextPin;
      
      // Change the x-coordinate:
      p->x = newPos;

      // Re-insert pin into list:
      add(p);
      return;
    }
    prevPin = pin;
    pin = nextPin;
  }
}

//----------------------------------------------------------------------------
/** Move all alpha pins relative to their current position.
  @param dist value to offset all alpha pins
*/
void vvPinList::moveAlphaRel(float dist)
{
  vvDebugMsg::msg(2, "vvPinList::moveAlphaRel()");

  vvPin* pin;

  if (root==NULL) return;     // no pins in list
  pin = root;
  while (pin != NULL)
  {
    if (pin->type==vvPin::ALPHA_HAT || 
				pin->type==vvPin::ALPHA_RAMP ||
        pin->type==vvPin::ALPHA_BLANK) 
    {
      pin->x += dist;
    }
    pin = pin->next;
  }
}

//----------------------------------------------------------------------------
/** Delete all pins of given pin type from the list.
  @param pt  pin type to delete
*/
void vvPinList::clearPins(vvPin::PinType pt)
{
  vvPin* pin;
  vvPin* nextPin;
  vvPin* prevPin;

  vvDebugMsg::msg(3, "vvPinList::clearPins()");
  pin = prevPin = root;
  while (pin != NULL)
  {
    nextPin = pin->next;    
    if (pin->type == pt)
    {
      if (pin==root) root = nextPin;
      else prevPin->next = nextPin;
      delete pin;
    }
    else prevPin = pin;
    pin = nextPin;
  }
}

//----------------------------------------------------------------------------
/// Delete all pins of all pin types from the list.
void vvPinList::clear()
{
  vvPin* pin;
  vvPin* nextPin;

  vvDebugMsg::msg(3, "vvPinList::clear()");
  pin = root;
  while (pin != NULL)
  {
    nextPin = pin->next;    
    delete pin;
    pin = nextPin;
  }
  root = NULL;
}

//----------------------------------------------------------------------------
/// Delete all color pins from the list
void vvPinList::clearColor()
{
  clearPins(vvPin::COLOR);
}

//----------------------------------------------------------------------------
/// Delete all alpha pins from the list
void vvPinList::clearAlpha()
{
  clearPins(vvPin::ALPHA_HAT);
  clearPins(vvPin::ALPHA_RAMP);
  clearPins(vvPin::ALPHA_BLANK);
}

//----------------------------------------------------------------------------
/// Returns the total number of pins in the list
int vvPinList::count()
{
  vvPin* pin;
  int i=0;

  vvDebugMsg::msg(3, "vvPinList::count()");
  if (root==NULL) return 0;     // no pins in list
  pin = root;
  while (pin != NULL)
  {
    ++i;
    pin = pin->next;
  }
  return i;
}

//----------------------------------------------------------------------------
/// Returns true if there are no pins of any type in the pin list.
bool vvPinList::isEmpty()
{
  if (root==NULL) return true;
  else return false;
}

//----------------------------------------------------------------------------
/// Returns the total number of color pins in the list
int vvPinList::countColor()
{
  vvPin* pin;
  int i=0;

  vvDebugMsg::msg(2, "vvPinList::countColor()");
  if (root==NULL) return 0;     // no pins in list
  pin = root;
  while (pin != NULL)
  {
    if (pin->type==vvPin::COLOR) ++i;
    pin = pin->next;
  }
  return i;
}

//----------------------------------------------------------------------------
/// Returns the total number of alpha pins in the list
int vvPinList::countAlpha()
{
  vvPin* pin;
  int i=0;

  vvDebugMsg::msg(2, "vvPinList::countAlpha()");
  if (root==NULL) return 0;     // no pins in list
  pin = root;
  while (pin != NULL)
  {
    if (pin->type==vvPin::ALPHA_HAT || 
				pin->type==vvPin::ALPHA_RAMP ||
        pin->type==vvPin::ALPHA_BLANK) 
      ++i;
    pin = pin->next;
  }
  return i;
}

//----------------------------------------------------------------------------
/** Compares two pin lists.
  @return true if lists are equal, otherwise false
*/
bool vvPinList::isEqual(vvPinList* otherList)
{
	int i;
	vvPin *ownPin, *otherPin;
	
  vvDebugMsg::msg(2, "vvPinList::isEqual()");
	if (root==NULL && otherList->root==NULL) return true; // both lists are empty
	ownPin = root;
	otherPin = otherList->root;
	while (ownPin != NULL && otherPin != NULL)
	{
		// Check all pin attributes for equality and return
		// as soon as a value differs:
		if (ownPin->type != otherPin->type) return false;
		if (ownPin->x != otherPin->x) return false;
	  for (i=0; i<3; ++i)
		  if (ownPin->v[i] != otherPin->v[i]) return false;
		ownPin = ownPin->next;
		otherPin = otherPin->next;
	}
	if (ownPin==NULL && otherPin==NULL) return true;
	else return false;
}

//----------------------------------------------------------------------------
/** Computes and returns the y value for a given x value,
    interpolating linearly between pins.
 @param x   x value for which to compute y value
 @param ct  component type for which to compute y value (RED, GREEN, BLUE, ALPHA)
*/
float vvPinList::computeY(float x, CompType ct)
{
  vvPin*  pin;              // current pin
  float   alpha;            // current alpha value
  float   maxAlpha;         // accumulated alpha value
  float   centerX;          // x position of range center for discrete colors
  float   rangeWidth;       // width of each discrete colors range
  int     currentRange;     // range index in which the current x value is located

  assert(ct==RED || ct==GREEN || ct==BLUE || ct==ALPHA);
  if (root==NULL || x<0.0f || x>1.0f)       // no pins in list
    return 0.0f;

  if (ct==ALPHA)
  {
    pin = root;
    maxAlpha = 0.0f;
    
    // Find maximum alpha value among alpha pins:
    while (pin != NULL)
    {
      switch (pin->type)
      {
        case vvPin::ALPHA_HAT:
          alpha = computeAlphaHatY(pin, x);
          break;
        case vvPin::ALPHA_RAMP:
          alpha = computeAlphaRampY(pin, x);
          break;
        case vvPin::ALPHA_BLANK:      // blank pins are dominant
          if (x >= (pin->x - pin->v[0] / 2.0f) && x <= (pin->x + pin->v[0] / 2.0f))
            return 0.0f;
          else alpha = 0.0f;
          break;
        default: alpha = 0.0f; break;
      }
      maxAlpha = ts_max(alpha, maxAlpha);
      pin = pin->next;
    }
    return maxAlpha;    
  }
  else    // red, green, or blue
  {
    if (discreteColors<=0)
      return computeLinearColorY(x, ct);     // using smooth colors
    else
    {
      // Using discrete colors: determine color range which the
      // value belongs to and look up its center value for the 
      // color:
      rangeWidth = 1.0f / discreteColors;     
      currentRange = int(x * discreteColors);
      if (currentRange>=discreteColors)     // constrain range to valid ranges
        currentRange = discreteColors - 1;
      centerX = currentRange * rangeWidth + (rangeWidth / 2.0f);
      return computeLinearColorY(centerX, ct);
    }
  }
}

//----------------------------------------------------------------------------
/** Computes and returns the y color value for a given x value,
    interpolating linearly between pins.
 @param x   x value for which to compute y value
 @param ct  component type for which to compute y value (RED, GREEN, BLUE)
*/
float vvPinList::computeLinearColorY(float x, CompType ct)
{
  vvPin*  pin;              // current pin
  int     index;            // component index
  float   prevX, prevY;     // previous x and y values

  assert(ct==RED || ct==GREEN || ct==BLUE);
  if (root==NULL || x<0.0f || x>1.0f)       // no pins in list
    return 0.0f;
  // Find component index:
  if      (ct==RED)   index=0;
  else if (ct==GREEN) index=1;
  else                index=2;
  prevX = 0.0f;
  prevY = 1.0f;
  pin = root;
  for (;;)
  {
    if (pin==NULL)    // done?
      return interpolateLinear(prevX, prevY, 1.0f, 1.0f, x);

    else if (pin->type==vvPin::COLOR)
    {
      if (x < pin->x)    // is x in current interval?
        return interpolateLinear(prevX, prevY, pin->x, pin->v[index], x);
      else if (x == pin->x)
        return pin->v[index];
      else
      {
        prevX = pin->x;
        prevY = pin->v[index];
      }
    }

    pin = pin->next;
  }
}

//----------------------------------------------------------------------------
/** Returns H, S, and B values for an x value.
 @param x      x value to get HSB values for [0..1]
 @param h,s,v  pointers to float variables to receive HSB values [0..1] 
*/
void vvPinList::getHSB(float x, float* h, float* s, float* v)
{
  float r,g,b;

  vvDebugMsg::msg(3, "vvPinList::getHSB()");
  r = computeY(x, RED);
  g = computeY(x, GREEN);
  b = computeY(x, BLUE);
  vvToolshed::RGBtoHSB(r, g, b, h, s, v);
}

//----------------------------------------------------------------------------
/** Returns red, green, and blue values for an x value.
 @param x      x value to get HSB values for [0..1]
 @param r,g,b  pointers to float variables to receive RGB values [0..1] 
*/
void vvPinList::getRGB(float x, float* r, float* g, float* b)
{
  vvDebugMsg::msg(3, "vvPinList::getRGB()");
  *r = computeY(x, RED);
  *g = computeY(x, GREEN);
  *b = computeY(x, BLUE);
}

//----------------------------------------------------------------------------
/** Returns the alpha value for an x value.
 @param x      x value to get the alpha value for [0..1]
 @return alpha value at x position [0..1]
*/
float vvPinList::getAlpha(float x)
{
  vvDebugMsg::msg(3, "vvPinList::getAlpha()");
  return computeY(x, ALPHA);
}

//----------------------------------------------------------------------------
/** Compute y value for an x value of an alpha hat.
  The hat consists of two flanks and a plateau.
*/
float vvPinList::computeAlphaHatY(vvPin* pin, float x)
{
  float y;                    // wanted y value
  float platBegin, platEnd;   // start and end x values of plateau
  
  assert(pin->type==vvPin::ALPHA_HAT);

  // Initialization:
  pin->v[0] = (float)fabs(pin->v[0]);   // slope must be positive
  platBegin = pin->x - pin->v[2] / 2.0f;
  platEnd   = pin->x + pin->v[2] / 2.0f;

  // Find point on hat edges:
  if (x <= pin->x)    // x is left of hat center
    y  = interpolateLinear(platBegin, pin->v[1], pin->v[0], x);
  else                // x is right of hat center
    y  = interpolateLinear(platEnd, pin->v[1], -pin->v[0], x);

  return ts_clamp(y, 0.0f, pin->v[1]);  // constrain to min/max values
}

//----------------------------------------------------------------------------
/// Compute y value for an x value of an alpha ramp.
float vvPinList::computeAlphaRampY(vvPin* pin, float x)
{
  float y;      // wanted y value

  assert(pin->type==vvPin::ALPHA_RAMP);
  y  = interpolateLinear(pin->x, pin->v[1] / 2.0f, pin->v[0], x);
  return ts_clamp(y, 0.0f, pin->v[1]);
}

//----------------------------------------------------------------------------
/** Interpolate linearly between two x/y value pairs.
 @param x1, y1  one x/y value pair
 @param x2, y2  another x/y value pair
 @param x       x value for which to find interpolated y value
 @return interpolated y value
*/
float vvPinList::interpolateLinear(float x1, float y1, float x2, float y2, float x)
{
  if (x1==x2) return ts_max(y1, y2);    // on equal x values: return maximum value
  if (x1 > x2)    // make x1 less than x2
  {
    ts_swap(x1, x2);
    ts_swap(y1, y2);
  }
  return (y2 - y1) * (x - x1) / (x2 - x1) + y1;  
}

//----------------------------------------------------------------------------
/** Interpolate linearly given a x/y value pair and a slope.
 @param x1, y1  point on straight line (x/y value pair)
 @param slope   slope of line
 @param x       x value for which to find interpolated y value
 @return interpolated y value
*/
float vvPinList::interpolateLinear(float x1, float y1, float slope, float x)
{
  return (y1 + slope * (x - x1));
}

//----------------------------------------------------------------------------
/** Copy a pin list from another pin list.
  @param pl  source pin list
*/
void vvPinList::copy(vvPinList* pl)
{
  vvPin* pin;

  vvDebugMsg::msg(3, "vvPinList::copy()");  
  clear();
  colorModel = pl->colorModel;
  pin = pl->root;
  while (pin != NULL)
  {
    add(new vvPin(pin));    
    pin = pin->next;
  }
}

//----------------------------------------------------------------------------
/** Make a float array from a pin list.
 @param pins _allocated_ array of pins (ELEMENTS_PER_PIN float values per pin)
*/
void vvPinList::makeArray(float* pins)
{
  vvPin* pin = root;
  int i, c;

  vvDebugMsg::msg(2, "vvPinList::makeArray()");
  i = 0;
  while (pin != NULL)
  {
    pins[i++] = (float)pin->type;
    pins[i++] = pin->x;
    for (c=0; c<3; ++c)
      pins[i++] = pin->v[c];
    pin = pin->next;
  }
}

//----------------------------------------------------------------------------
/** Make a pin list of pins from a float array of pin attributes.
 @param num   number of pins in list
 @param pins  array of pins: ELEMENTS_PER_PIN float values in range [0..1] per pin
*/
void vvPinList::makePins(int num, const float* pins)
{
  vvPin::PinType type;
  float x,v[3];
  int c;
  
  vvDebugMsg::msg(3, "vvPinList::makePins()");
  clear();
  for (int i=0; i<num; ++i)
  { 
    type = (vvPin::PinType)((int)pins[vvPin::ELEMENTS_PER_PIN * i]);
    x = pins[vvPin::ELEMENTS_PER_PIN * i + 1];
    for (c=0; c<3; ++c)
      v[c] = pins[vvPin::ELEMENTS_PER_PIN * i + c + 2];
    add(new vvPin(type, x, v[0], v[1], v[2]));
  }
}

//----------------------------------------------------------------------------
/** Compute function values to a float array. 
 Order of components: RRR..GGG..BBB.. (for RGB type array)
 @param entries  number of array entries to compute
 @param ct       component type: determines which components are to be written
                 to the array (e.g. "RED | GREEN | BLUE" for RGB)
 @param array    _allocated_ float array in which to store computed values.
                 values range is from 0.0 to 1.0. 
                 Space for 'entries' float values must be provided.
*/
void vvPinList::computeFunction(int entries, int ct, float* array)
{
  int c, i, index;
  CompType comp[4] = {RED, GREEN, BLUE, ALPHA};

  vvDebugMsg::msg(3, "vvPinList::computeFunction()");

  for (i=0; i<entries; ++i)
  {
    index = 0;
    for (c=0; c<4; ++c)
      if ((ct & comp[c]) != 0)
      {
        array[index * entries + i] = computeY((float)i / (float)(entries-1), comp[c]);
        ++index;
      }
  }
}

//----------------------------------------------------------------------------
/// Returns the currently selected color model
vvPinList::ColorModelType vvPinList::getColorModel()
{
  vvDebugMsg::msg(3, "vvPinList::getColorModel()");
  return colorModel;
}

//----------------------------------------------------------------------------
/** Set a new color model for color pins.
  @param cm  new color model: RGB_MODEL or HSB_MODEL
*/
void vvPinList::setColorModel(ColorModelType cm)
{
  vvPin* pin;

  vvDebugMsg::msg(3, "vvPinList::setColorModel()");

  if (colorModel==cm || root==NULL) return;    // nothing needs to be done

  pin = root;
  while (pin != NULL)
  {
    if (pin->type==vvPin::COLOR)
    {
      if (cm==HSB_MODEL)
        vvToolshed::RGBtoHSB(&pin->v[0], &pin->v[1], &pin->v[2]);
      else
        vvToolshed::HSBtoRGB(&pin->v[0], &pin->v[1], &pin->v[2]);
    }
    pin = pin->next;
  }
  colorModel = cm;
}

//----------------------------------------------------------------------------
/** Set the number of discrete colors to use for color interpolation.
  @param numColors number of discrete colors (use 0 for smooth colors)
*/
void vvPinList::setDiscreteColors(int numColors)
{
  vvDebugMsg::msg(3, "vvPinList::setDiscreteColors()");
  if (numColors<0) numColors = 0;
  discreteColors = numColors;
}

//============================================================================
// End of File
//============================================================================
