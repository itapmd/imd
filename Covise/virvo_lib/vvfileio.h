//****************************************************************************
// Project Affiliation: Virvo (Virtual Reality Volume Renderer)
// Copyright:           (c) 2002 Juergen Schulze-Doebold. All rights reserved.
// Author's E-Mail:     schulze@hlrs.de
// Institution:         University of Stuttgart, Supercomputing Center (HLRS)
//****************************************************************************

#ifndef _VVFILEIO_H_
#define _VVFILEIO_H_

#include <assert.h>
#include "vvvoldesc.h"

/** File load and save routines for volume data.
  The following file formats are supported:<BR>
  RVF, XVF, VF, AVF, 3D TIFF, Visible Human, raw, RGB, TGA, PGM, PPM<BR>
  When a 2D image file is loaded, the loader looks for a numbered sequence
  automatically, for instance: file001.rvf, file002.rvf, ...
  @author Juergen Schulze (schulze@hlrs.de)
*/
class vvFileIO
{
  public:
    class ParticleTimestep
    {
      public:
        int numParticles;   // number of particles
        float* pos[3];      // x/y/z positions
        float* val;         // values
        float min,max;      // minimum and maximum scalar values        
        ParticleTimestep(int np)
        {
          numParticles = np;
          min = max = 0;
          for(int i=0; i<3; ++i) 
          {
            pos[i] = new float[numParticles];
            assert(pos[i]);
          }
          val = new float[numParticles];
          assert(val);
        }
        ~ParticleTimestep()
        {
          for(int i=0; i<3; ++i) delete[] pos[i];
          delete[] val;
        }
    };
    enum ErrorType      /// Error Codes
    {
      OK,               ///< no error
      PARAM_ERROR,      ///< parameter error
      FILE_ERROR,       ///< file IO error
      FILE_EXISTS,      ///< file exists error
      FILE_NOT_FOUND,   ///< file not found error
      DATA_ERROR,       ///< data format error
      FORMAT_ERROR,     ///< file format error (e.g. no valid TIF file)
      VD_ERROR          ///< volume descriptor (vvVolDesc) error
    };
    enum LoadType        /// Load options
    {
      ALL_DATA = 0xFFFF, ///< load all data
      HEADER   = 0x0001, ///< load header
      ICON     = 0x0002, ///< load icon
      RAW_DATA = 0x0004, ///< load volume raw data
      TRANSFER = 0x0008  ///< load transfer functions
    };

    vvFileIO();  
    ErrorType saveVolumeData(vvVolDesc*, bool, int = ALL_DATA);
    ErrorType loadVolumeData(vvVolDesc*, int = ALL_DATA);
    ErrorType loadDicomFile(vvVolDesc*, int* = NULL, int* = NULL, float* = NULL);
    ErrorType loadRawFile(vvVolDesc*, int, int, int, int, int);
    ErrorType loadXB7File(vvVolDesc*,int=128,int=8,bool=true);
    ErrorType loadCPTFile(vvVolDesc*,int=128,int=8,bool=true);
    void      setCompression(bool);

  private:
    char xvfID[10];     ///< XVF file ID
    char nrrdID[9];     ///< nrrd file ID
    int  sections;      ///< bit coded list of file sections to load
    bool compression;   ///< true = compression on (default)

    int       readASCIIint(FILE*);
    void      writePinList(FILE*, vvPinList*);
    ErrorType loadWLFile(vvVolDesc*);
    ErrorType loadASCFile(vvVolDesc*);
    ErrorType saveRVFFile(vvVolDesc*);
    ErrorType loadRVFFile(vvVolDesc*);
    ErrorType saveXVFFile(vvVolDesc*);
    ErrorType loadXVFFile(vvVolDesc*);
    ErrorType loadAVFFile(vvVolDesc*);
    ErrorType saveVFFile(vvVolDesc*);
    ErrorType loadVFFile(vvVolDesc*);
    ErrorType loadTIFFile(vvVolDesc*);
    ErrorType loadRawFile(vvVolDesc*);
    ErrorType saveRawFile(vvVolDesc*);
    ErrorType loadRGBFile(vvVolDesc*);
    ErrorType loadTGAFile(vvVolDesc*);
    ErrorType loadPXMRawImage(vvVolDesc*);
    ErrorType loadVHDAnatomicFile(vvVolDesc*);
    ErrorType loadVHDMRIFile(vvVolDesc*);
    ErrorType loadVHDCTFile(vvVolDesc*);
    ErrorType loadVMRFile(vvVolDesc*);
    ErrorType loadVTCFile(vvVolDesc*);
    ErrorType loadNrrdFile(vvVolDesc*);
    ErrorType saveNrrdFile(vvVolDesc*);
    ErrorType loadXIMGFile(vvVolDesc*);
    ErrorType computeDefaultVolume(vvVolDesc*, int);
    ErrorType savePXMSlices(vvVolDesc*, bool);
};

#endif

//============================================================================
// End of File
//============================================================================
