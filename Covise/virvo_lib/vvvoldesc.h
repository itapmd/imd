//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifndef _VVVOLDESC_H_
#define _VVVOLDESC_H_

#include <stdlib.h>
#include "vvtoolshed.h"
#include "vvvecmath.h"
#include "vvtransfunc.h"
#include "vvsllist.h"

//============================================================================
// Class Definition
//============================================================================

/** Volume description.
  The volume description contains basically all the elements which describe
  the volume data. Most of this information can be saved to a file using
  the appropriate file format.
  @author Juergen Schulze-Doebold (schulze@hlrs.de)
*/
class vvVolDesc               
{
  public: 
    enum ErrorType                ///  error Codes
    {
      OK,                         ///< no error
      TYPE_ERROR                  ///< data type error
    };
    enum AxisType                 ///  names for coordinate axes
    { X_AXIS, Y_AXIS, Z_AXIS }; 
    enum InterpolationType        ///  interpolation types to use for resampling
    {                           
      NEAREST,                    ///< nearest neighbor
      TRILINEAR                   ///< trilinear
    };
    enum DeleteType               /// types for data deletion
    {
      DELETE_DATA,                ///< delete data when not used anymore
      NO_DELETE                   ///< don't delete data because it will be deleted by the caller
    };
    enum NormalizationType        /// type of normalization
    { 
      VV_LINEAR,                  ///< linear
      VV_LOGARITHMIC              ///< logarithmic
    };
    enum StorageType              /// voxel data storage type
    {
      VV_INTEGER = 0,             ///< storage of integer values, depends on bpv: 1=8bit, 2=16bit, 3=RGB, 4=RGBA
      VV_FLOAT,                   ///< storage of floating point values, depends on bpv: 4=single precision, 8=double precision
      VV_MULTIMODAL               ///< storage of multi-modal data with one byte per modality, value of bpv gives number of modalities
    };
    enum                          ///< size of serialization buffer for attributes [bytes]
    { SERIAL_ATTRIB_SIZE = 3*4 + 4 + 1 + 3*4 + 4 + 2*4 + 3*4 + 1}; // 54
    
    int   vox[3];                 ///< width, height and number of slices of volume [voxels] (0 if no volume loaded)
    int   frames;                 ///< number of animation frames in movie (0 if no volume data stored)
    int   bpv;                    ///< byte per voxel (default = 1):<UL>
                                  ///< <LI>1 = 1 uchar density value</LI>
                                  ///< <LI>2 = 16 bit unsigned density value (12 bit data must be located in 12 most significant bits, the rest set to zero)</LI>
                                  ///< <LI>3 = 3 uchar RGB values</LI>
                                  ///< <LI>4 = 3 uchar RGB values + 1 uchar density value</LI></UL>
    StorageType stype;            ///< voxel data storage type
    float dist[3];                ///< Distance between sampling points in x/y/z direction [mm]
    float dt;                     ///< Length of an animation time step [seconds]
    float realMin;                ///< physical equivalent of minimum scalar value (only applies to 8 or 16 bit data)
    float realMax;                ///< physical equivalent of maximum scalar value (only applies to 8 or 16 bit data)
    vvVector3 position;           ///< location of volume center [mm]
    vvTransFunc tf;               ///< transfer functions
    bool  simulate;               ///< true = only simulate conversion routines: don't actually modify volume data, modify volume header only

    // Constructors and destructors:
    vvVolDesc();                    
    vvVolDesc(const char*);
    vvVolDesc(const char*, int, int, int, int, int, uchar**);
    vvVolDesc(const char*, int, int, int, int, float**);
    vvVolDesc(const char*, int, int, int, int, float**, float**, float**);
    vvVolDesc(const char*, int, int, uchar*);
    vvVolDesc(vvVolDesc*, int);
    virtual ~vvVolDesc();    

    // Getters and setters:
    int    getSliceSize();
    int    getFrameSize();
    int    getMovieSize();
    int    getSliceVoxels();
    int    getFrameVoxels();
    int    getMovieVoxels();
    uchar* getRaw(int);
    uchar* getRaw();
    void   setVoxelSize(int, int, int, int);
    void   setRealSize(vvVector3*);
    void   getRealSize(vvVector3*);
    const char* getFilename();
    void   setFilename(const char*);
    int    getIconSize();
    uchar* getIconData();
    void   setCurrentFrame(int);
    int    getCurrentFrame();

    // Conversion routines:
    void   convertVoxelFormat(int, bool=false);
    void   convertModalities(int, bool=false);
    void   bitShiftData(int, bool=false);
    void   invert();
    void   convertRGB24toRGB8();
    void   flip(AxisType);
    void   rotate(AxisType, int);
    void   convertRGBPlanarToRGBInterleaved();
    void   toggleEndianness();
    void   crop(int, int, int, int, int, int);
    void   cropTimesteps(int, int);
    void   resize(int, int, int, InterpolationType, bool=false);
    void   shift(int, int, int);
    void   convertCoviseToVirvo();
    void   convertVirvoToCovise();
    void   convertVirvoToOpenGL();
    void   convertOpenGLToVirvo();
    void   makeSphere(int, int, InterpolationType, bool=false);
    void   expandDataRange(bool = false);
    void   zoomDataRange(int, int, bool = false);
    void   toggleSign();
    void   blend(vvVolDesc*, float, bool=false);

    // Other routines:
    ErrorType merge(vvVolDesc*);
    void   addFrame(uchar*, DeleteType);
    int    getStoredFrames();
    void   copyFrame(uchar*);
    void   removeSequence();
    void   makeHistogram(int, int, int*);
    void   normalizeHistogram(int, int*, float*, NormalizationType);
    void   makeHistogramTexture(int, int, bool, uchar*);
    bool   isByteUsed(int);
    int    getAlphaRange();
    void   makeIcon(int, const uchar*);
    void   printInfoLine(const char* = NULL);
    void   printVolumeInfo();
    void   printStatistics();
    void   printVoxelData(int, int, int=0, int=0);
    void   printHistogram(int);
    void   trilinearInterpolation(int, float, float, float, uchar*);
    void   drawLine(int, int, int, int, int, int, int);
    int    serializeAttributes(uchar* = NULL);
    void   deserializeAttributes(uchar*, int=SERIAL_ATTRIB_SIZE);
    void   setSliceData(uchar*, int=0, int=0);
    void   extractSlice(int, int, int, uchar*);
    void   deinterlace();
    void   findMinMax(int*, int*);
    int    findNumValue(int, int);
    int    findNumUsed();
    int    findNumTransparent(int);
    
  private:
    char*  filename;          ///< name of volume data file, including extension, excluding path ("" if undefined)
    int    currentFrame;      ///< current animation frame
    int    iconSize;          ///< size = width and height of icon [pixels] (0 if no icon stored)
    uchar* iconData;          ///< icon image data as RGB triples (RGB, RGB, ...), starting top left, 
                              ///< then going right, then down
    vvSLList<uchar*> raw;     ///< pointer list to raw volume data
    
    void initialize(); 
    void setDefaults();
};

#endif

//============================================================================
// End of File
//============================================================================
