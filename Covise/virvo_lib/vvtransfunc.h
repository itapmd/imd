//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifndef _VVTRANSFUNC_H_
#define _VVTRANSFUNC_H_

#include "vvtoolshed.h"

/** A single transfer function pin.
  @author Juergen Schulze-Doebold (schulze@hlrs.de)
  @see vvTransFunc
*/
class vvPin                       
{
  private:
    static int idCounter;       ///< pin ID counter
    int id;                     ///< unique pin identifier

  public:
    enum ElementsType           /// elements per pin: type, x, v[0], v[1], v[2]
    { ELEMENTS_PER_PIN = 5 };  
    enum PinType                /// type of pin
    {
      COLOR       = 0,          ///< color pin (RGB)
      ALPHA_HAT   = 1,          ///< alpha hat
      ALPHA_RAMP  = 2,          ///< alpha ramp
      ALPHA_BLANK = 3           ///< alpha blank (this pin type dominates all other alpha pins)
    };
    PinType type;               ///< pin type
    float x;                    ///< scalar value for which pin is valid [0..1]<BR>
                                ///< Important: this value may not be changed while 
                                ///< a pin is in the pin list. Remove pin from the list
                                ///< and add as a new pin instead, because the pin list
                                ///< must always remain sorted!<BR>
    float v[3];                 ///< storage for up to 3 pin property values:
                                ///< <UL>
                                ///<   <LI>color pin:<UL> 
                                ///<       <LI>0 = red/hue [0..1]</LI>
                                ///<       <LI>1 = green/saturation [0..1]</LI>
                                ///<       <LI>2 = blue/brightness [0..1]</LI></UL></LI>
                                ///<   <LI>alpha hat:<UL> 
                                ///<       <LI>0 = flank slope [0..oo]</LI>
                                ///<       <LI>1 = maximum value [0..1]</LI>
                                ///<       <LI>2 = width at top [0..1] (0=peak)</LI></UL></LI>
                                ///<   <LI>alpha ramp:<UL> 
                                ///<       <LI>0 = flank slope [-oo..oo]</LI>
                                ///<       <LI>1 = maximum value [0..1]</LI></UL></LI>
                                ///<   <LI>alpha blank: 0 = width of transparent area</LI>
                                ///< </UL>
    vvPin* next;                ///< pointer to next pin (NULL if no more pins in list)

    vvPin(PinType, float, float, float=0.0f, float=0.0f);
    vvPin(vvPin*);
    int getID();
};

/** A list of pins for RGBA transfer functions.
  @author Juergen Schulze-Doebold (schulze@hlrs.de)
  @see vvPin
  @see vvTransFunc
*/
class vvPinList                   
{         
  public:
    enum ColorModelType         /// color models
    {
      RGB_MODEL = 0,
      HSB_MODEL = 1
    };
  protected:  
    vvPin* root;                  ///< pointer to root pin (NULL if no pins stored)
    ColorModelType colorModel;    ///< currently selected color model        
    int discreteColors;           ///< number of discrete colors, copied from vvTransfunc

    float interpolateLinear(float, float, float, float, float);
    float interpolateLinear(float, float, float, float);
    void  clearPins(vvPin::PinType);

  public:
    enum CompType     /// component type: to be or'ed for composite forms: RED | GREEN | BLUE = RGB
    {
      RED   = 1,
      GREEN = 2,
      BLUE  = 4,
      ALPHA = 8,
      RGB   = RED | GREEN | BLUE,
      RGBA  = RGB | ALPHA
    };

    vvPinList();
    virtual ~vvPinList();
    void  add(vvPin*);
    vvPin*  getPin(int);
    void  remove(int);
    void  remove(vvPin*);
    void  move(vvPin*, float);
    void  moveAlphaRel(float);
    void  clear();
    void  clearColor();
    void  clearAlpha();
    int   count();
    int   countColor();
    int   countAlpha();
		bool  isEqual(vvPinList*);
    bool  isEmpty();
    float computeY(float, CompType);
    void  getHSB(float, float*, float*, float*);
    void  getRGB(float, float*, float*, float*);
    float getAlpha(float);
    float computeAlphaHatY(vvPin*, float);
    float computeAlphaRampY(vvPin*, float);
    float computeLinearColorY(float, CompType);
    void  copy(vvPinList*);
    void  makeArray(float*);
    void  makePins(int, const float*);
    void  computeFunction(int, int, float*);
    void  setColorModel(ColorModelType);
    void  setDiscreteColors(int);
    ColorModelType getColorModel();
};

/** Description of a transfer function.
  @author Juergen Schulze-Doebold (schulze@hlrs.de)
  @see vvPin
  @see vvPinList
*/
class vvTransFunc               
{
  private:
    enum PresetsType                /// number of preset transfer functions
    { NUM_PRESETS = 8 };  
    enum BufferLengthType           /// number of elements in ring buffer
    { BUFFER_SIZE = 20 };             
    vvPinList buffer[BUFFER_SIZE];  ///< ring buffer which can be used to implement Undo functionality
    int nextBufferEntry;            ///< index of next ring buffer entry to use for storage
    int bufferUsed;                 ///< number of ring buffer entries used
    int discreteColors;             ///< number of discrete colors to use for color interpolation (0 for smooth colors)

  public: 
    vvPinList pins;                 ///< pin list for current transfer function
    vvPinList presets[NUM_PRESETS]; ///< pin lists for preset transfer functions

    vvTransFunc();
    virtual ~vvTransFunc();
    void setDefaultColors(int=0);
		int  getNumDefaultColors();
    void setDefaultAlpha(int=0);
		int  getNumDefaultAlpha();
    int  getNumUsedPresets();
    bool isDefaultColors(int=0);
    bool isDefaultAlpha(int=0);
    bool isEmpty();
    void makeColorBar(int, uchar*);
    void make8bitLUT(int, uchar*);
    void makeFunctionTexture(int, int, uchar*, bool);
    void copy(vvTransFunc*);
    void copyPresets(vvTransFunc*);
    void storePreset(int);
    void restorePreset(int);
    int  getNumPresets();
    void clear();
    void putUndoBuffer();
    void getUndoBuffer();
    void clearUndoBuffer();
    void setDiscreteColors(int);
    int  getDiscreteColors();
};

#endif

//============================================================================
// End of File
//============================================================================
