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
  using std::setw;
  using std::hex;
  using std::dec;
#else
  #include <iostream.h>
  #include <iomanip.h>
#endif

#include <math.h>
#include <limits.h>
#include <ctype.h>
#include "vvfileio.h"
#include "vvvffile.h"
#include "vvtoolshed.h"
#include "vvdebugmsg.h"
#include "vvtokenizer.h"
#include "vvdicom.h"
#include "vvarray.h"

#ifdef __sun
  #define powf pow
#endif

//----------------------------------------------------------------------------
/// Constructor
vvFileIO::vvFileIO()
{
  vvDebugMsg::msg(1, "vvFileIO::vvFileIO()");
  assert(sizeof(float) == 4);
  strcpy(xvfID, "VIRVO-XVF");
  strcpy(nrrdID, "NRRD0001");
  sections = ALL_DATA;
  compression = true;
}

//----------------------------------------------------------------------------
/** Read the next ASCII integer character from a file. 
  Ignores all other characters (as well as CR, LF etc.).
  @return 0 if no integer string could be found
*/
int vvFileIO::readASCIIint(FILE* src)
{
  char intString[64];     // string to hold integer value
  int  i;                 // index to current letter in integer string
  char c;                 // most recently read character
  bool done = false;      // true if integer reading done
  bool found = false;     // true if ASCII integer reading has begun

  i = 0;
  while (!done && !feof(src))
  {
    c = (char)fgetc(src);
    if (c>='0' && c<='9')
    {
      found = true;
      intString[i] = c;
      ++i;
    }
    else if (found==true) done = true;
  }
  intString[i] = '\0';
  return atoi(intString);
}

//----------------------------------------------------------------------------
/** Writes a pin list to a file.
*/
void vvFileIO::writePinList(FILE* fp, vvPinList* pins)
{
  float* pinArray;    // serialized pins
  int numPins;        // number of pins in the list
  int i;
  
  numPins = pins->count();
  pinArray = new float[vvPin::ELEMENTS_PER_PIN * numPins];
  pins->setColorModel(vvPinList::RGB_MODEL);
  pins->makeArray(pinArray);
  for (i=0; i<vvPin::ELEMENTS_PER_PIN * numPins; ++i)
    vvToolshed::writeFloat(fp, pinArray[i]);
  delete[] pinArray;
}

//----------------------------------------------------------------------------
/// Loader for voxel files in worldlines format.
vvFileIO::ErrorType vvFileIO::loadWLFile(vvVolDesc* vd)
{
  FILE* fp;
  int val[3];                     // scalar voxel value ([0] = x, [1] = y, [2] = z)
  int i;                          // index counter
  int col = 0;                    // voxel scalar value
  int max[3]={1, 1, 1};           // maximum values in each dimension
  bool done;                      // true if done
  uchar* raw;                     // raw volume data

  vvDebugMsg::msg(1, "vvFileIO::loadWLFile()");

  if (vd->getFilename()==NULL) return FILE_ERROR;
  fp = fopen(vd->getFilename(), "r");
  if (fp==NULL) return FILE_ERROR;

  vd->removeSequence();

  // Find maximum dimension values:  
  done = false;
  while (!done)
  {
    for (i=0; i<3; ++i)
    {
      if (feof(fp)) done = true;
      else 
      {
        val[i] = readASCIIint(fp);    
        if ((val[i]+1) > max[i]) max[i] = val[i] + 1;   // adjust maximum value. +1 needs to be there to include index 0 
      }
    }
  }
  
  vd->vox[0] = max[0];
  vd->vox[1] = max[1];
  vd->vox[2] = max[2];
  vd->frames = 1;
  
  // Allocate memory for volume:    
  raw = new uchar[vd->getFrameSize()];

  // Initialize volume data:  
  for (i=0; i<vd->getFrameSize(); ++i)
    raw[i] = (uchar)0;

  // Read volume data:
  fseek(fp, 0, SEEK_SET);
  done = false;
  while (!done)
  {
    for (i=0; i<3; ++i)
    {
      if (feof(fp)) done = true;
      else 
      {
        val[i] = readASCIIint(fp);    
        if (val[i] > max[i])  // safety check
        {
          vvDebugMsg::msg(1, "Error: Voxel coordinate exceeds limit."); 
          fclose(fp);
          delete[] raw;
          return FILE_ERROR;
        }
      }
    }
    raw[val[0] + val[1] * max[0] + val[2] * max[0] * max[1]] = uchar(col);
    ++col;
    if (col > 255) col = 0;
  }
  
  vd->addFrame(raw, vvVolDesc::DELETE_DATA);
  fclose(fp);  
  return OK;
}

//----------------------------------------------------------------------------
/** Loader for voxel file in ASCII format.
  3 integers in first line indicating size of 3 dimensions (XSIZE, YSIZE, ZSIZE),
  then quadruples of integers indicating visible voxel position (0..XSIZE-1 etc.) 
  and opacity (0..255).<BR>
  1 byte opacity value per voxel, 0=transparent, 255=opaque
*/
vvFileIO::ErrorType vvFileIO::loadASCFile(vvVolDesc* vd)
{
  FILE* fp;
  int x, y, z, op, i;
  uchar* raw;   // raw volume data
  
  vvDebugMsg::msg(1, "vvFileIO::loadASCFile()");

  if (vd->getFilename()==NULL) return FILE_ERROR;
  fp = fopen(vd->getFilename(), "rt");
  if (fp==NULL)
  {
    vvDebugMsg::msg(1, "Error: Cannot open file.");
    return FILE_ERROR;
  }

  vd->removeSequence();  
  
  fscanf(fp, "%d %d %d", &x, &y, &z);
  if (x>0 && y>0 && z>0) 
  {
    vd->vox[0] = x;
    vd->vox[1] = y;
    vd->vox[2] = z;
    vd->frames = 1;
  }
  else    // invalid volume dimensions
  {
    fclose(fp);
    vvDebugMsg::msg(1, "Error: Invalid header in ASC file.");
    return FILE_ERROR;
  }
  
  raw = new uchar[vd->getFrameSize()];
  
  for (i=0; i<vd->getFrameSize(); ++i)
    raw[i] = (uchar)0;    // initialize with opacity 0
  
  while (!feof(fp))
  {
    fscanf(fp, "%d %d %d %d", &x, &y, &z, &op);
    if (x>vd->vox[0]-1 || y>vd->vox[1]-1 || z>vd->vox[2]-1 || x<0 || y<0 || z<0 || op<0 || op>255)
    {
      vvDebugMsg::msg(1, "Error: Invalid value in ASC file.");
      fclose(fp);
      delete[] raw;
      return FILE_ERROR;
    }
    raw[x + y * vd->vox[0] + z * vd->vox[0] * vd->vox[1]] = (uchar)op;
  }  
  fclose(fp);

  vd->addFrame(raw, vvVolDesc::DELETE_DATA);
  return OK;
}

//----------------------------------------------------------------------------
/** Save current frame to a .RVF (raw volume data) file.
 Only volume dimensions and 8 bit raw data is saved in this format.
*/
vvFileIO::ErrorType vvFileIO::saveRVFFile(vvVolDesc* vd)
{
  FILE* fp;             // volume file pointer
  int frameSize;        // size of a frame in bytes
  uchar* raw;           // raw volume data
  vvVolDesc* v;           // temporary volume description
  
  vvDebugMsg::msg(1, "vvFileIO::saveRVFFile()");

  // Save volume data:
  v = new vvVolDesc(vd, vd->getCurrentFrame());   // copy current frame to a new VD
  if (vd->bpv!=1) 
  {
    cerr << "Converting data to 1 bpv" << endl;
    v->convertVoxelFormat(1);
  }
  frameSize = v->getFrameSize();
  raw = v->getRaw();  // save only current frame of loaded sequence
  if (frameSize==0 || raw==NULL) 
  {
    delete v;
    return VD_ERROR;
  }

  if ( (fp = fopen(v->getFilename(), "wb")) == NULL) // now open file to write
  {
    vvDebugMsg::msg(1, "Error: Cannot open file to write.");
    delete v;
    return FILE_ERROR;
  }

  vvToolshed::write16(fp, (ushort)v->vox[0]);
  vvToolshed::write16(fp, (ushort)v->vox[1]);
  vvToolshed::write16(fp, (ushort)v->vox[2]);

  // Write volume data:
  if ((int)fwrite(raw, 1, frameSize, fp) != frameSize)
  {
    vvDebugMsg::msg(1, "Error: Cannot write voxel data to file.");
    fclose(fp);
    delete v;
    return FILE_ERROR;
  }
  
  fclose(fp);
  delete v;
  return OK;
}

//----------------------------------------------------------------------------
/** Loader for voxel file in rvf (raw volume file) format.
  File specification (byte order: most significant first = big endian):<P>
  Header: 3 x 16 bit for width, height, slices<BR>
  Data: width x height x slices bytes 8 bit voxel data;
        order: top left front first, then to right,
        then to bottom, then to back
*/
vvFileIO::ErrorType vvFileIO::loadRVFFile(vvVolDesc* vd)
{
  FILE* fp;             // volume file pointer
  int frameSize;        // size of a frame in bytes
  uchar* raw;           // raw volume data
   
  vvDebugMsg::msg(1, "vvFileIO::loadRVFFile()");

  if (vd->getFilename()==NULL) return FILE_ERROR;
  if ( (fp = fopen(vd->getFilename(), "rb")) == NULL) 
  {
    vvDebugMsg::msg(1, "Error: Cannot open file.");
    return FILE_ERROR;
  }

  vd->removeSequence();

  // Read header:
  vd->vox[0] = vvToolshed::read16(fp);
  vd->vox[1] = vvToolshed::read16(fp);
  vd->vox[2] = vvToolshed::read16(fp);
  vd->frames = 1;
  vd->bpv    = 1;

  // Create new data space for volume data:
  if ((sections & RAW_DATA) != 0)
  {
    frameSize = vd->getFrameSize();
    raw = new uchar[frameSize];
  
    // Load volume data:
    if ((int)fread(raw, 1, frameSize, fp) != frameSize)
    {
      vvDebugMsg::msg(1, "Error: Insufficient voxel data in RVF file.");
      fclose(fp);
      delete[] raw;
      return FILE_ERROR;
    }
  
    vd->addFrame(raw, vvVolDesc::DELETE_DATA);
  }
  fclose(fp);
  return OK;
}

//----------------------------------------------------------------------------
/** Save volume data to a .XVF (extended volume data) file.
 File specification (byte order: most significant first = big endian):
 <PRE>
 Header: 
   Offset Bytes Data Type       Description
 ---------------------------------------------------------------
      0   9     char            file ID string: "VIRVO-XVF"
      9   2     unsigned short  offset to beginning of data area, from top of file [bytes]
     11   2 x 4 unsigned int    width and height of volume [voxels]
     19   4     unsigned int    number of slices per time step
     23   4     unsigned int    number of frames in volume animation (time steps)
     27   1     unsigned char   bytes per voxel (for details see vvvoldesc.h)
     28   3 x 4 float           real world voxel size (width, height, depth) [mm]
     40   4     float           length of a time step in the volume animation [seconds]
     44   2 x 4 float           physical data range covered by voxel data (minimum, maximum)
     52   3 x 4 float           real world location of volume center (x,y,z) [mm]
     64   1     unsigned char   compression type (0=uncompressed, 1=RLE)
     65   2     unsigned short  number of transfer functions
     67   2     unsigned short  type of transfer function: 0 = 4 x 256 Byte, 
                                1 = list of control pins
     69   1     unsigned char   storage type (for details see vvvoldesc.h)                             
     
 Data area:           
   Data starts at "offset to data area".
   Voxel order: voxel at top left front first, then to right, 
   then to bottom, then to back, then frames. All bytes of each voxel
   are stored successively.
   In RLE encoding mode, a 4 byte big endian value precedes each frame, 
   telling the number of RLE encoded bytes that will follow. If this 
   value is zero, the frame is stored without encoding.

 Now follow the transfer functions. 
 Each transfer function of type 0 consists of:
 - Zero terminated description string
 - Transfer function data in RGBA format: 
   First all R's, then all G's, etc.
   Each R/G/B/A entry is coded as one unsigned byte.
   The table length depends on the number of bits per voxel: 
    8 bits per voxel => 4 * 256 bytes (RGBA)
   16 bits per voxel => 4 * 4096 bytes (RGBA)
   24 bits per voxel => 3 * 256 bytes (alpha conversion for RGB)
   32 bits per voxel => 1 * 256 bytes (alpha conversion for density)
 - Transfer function in pin format:
   Each pin consists of 5 float values. The list is terminated by
   5 -1.0 values.

 Hints: 
   The big endian hexadecimal representations of some important floating point values are:
    1.0 = 3F 80 00 00
   -1.0 = BF 80 00 00
 </PRE><P>

  Up to now, the following header sizes (=offset to data area) were used in files:<BR>
  28: ID string + offset + width + height + slices + frames + bpv<BR>
  48: + dist + dt + num Transf + type Transf<BR>
  69: + realMin + realMax + position + compression
  70: + storage type
*/
vvFileIO::ErrorType vvFileIO::saveXVFFile(vvVolDesc* vd)
{
  const int HEADER_SIZE = 70;             // header size [bytes]
  uchar serialized[vvVolDesc::SERIAL_ATTRIB_SIZE];  // space for serialized volume data
  vvPinList* pins;                        // pin list
  FILE* fp;                               // volume file pointer
  uchar* raw;                             // raw volume data
  uchar* encoded = NULL;                  // encoded volume data
  char tfName[16];                        // transfer function name
  int f;                                  // counter for frames
  int i, j;                               // counters
  int index;                              // preset index
  int frames;                             // volume animation frames
  int frameSize;                          // frame size
  int numTF;                              // number of transfer functions to save
  int encodedSize;                        // number of bytes in encoded array
    
  vvDebugMsg::msg(1, "vvFileIO::saveXVFFile()");

  // Prepare variables:
  frames = vd->frames;
  if (frames==0) return VD_ERROR;
  frameSize = vd->getFrameSize();

  // Open file:
  if ( (fp = fopen(vd->getFilename(), "wb")) == NULL) // now open file to write
  {
    vvDebugMsg::msg(1, "Error: Cannot open file to write.");
    return FILE_ERROR;
  }

  // Write header:
  fputs(xvfID, fp);
  vvToolshed::write16(fp, HEADER_SIZE);   // total header size in bytes
  vd->serializeAttributes(serialized);
  fwrite(serialized, vvVolDesc::SERIAL_ATTRIB_SIZE, 1, fp);

  vvToolshed::write8(fp, uchar((compression) ? 1 : 0));     // compression type
  numTF = vd->tf.getNumUsedPresets() + 1;
  if (numTF==1 && vd->tf.pins.isEmpty()) numTF = 0;   // don't save any TF if none is defined
  vvToolshed::write16(fp, (ushort)numTF);   // number of transfer functions
  vvToolshed::write16(fp, 1);       // type of transfer functions (1 = pin list)

  // Write volume data frame by frame:
  if (compression==1) encoded = new uchar[frameSize];
  for (f=0; f<frames; ++f)
  {
    raw = vd->getRaw(f);
    if (raw==NULL) 
    { 
      vvDebugMsg::msg(1, "Error: no data available for frame", f);
      fclose(fp);
      delete[] encoded;
      return VD_ERROR;
    }
    if (compression)
    {
      encodedSize = vvToolshed::encodeRLE(encoded, raw, frameSize, vd->bpv, frameSize);
      if (encodedSize>=0)  // compression possible?
      {
        vvToolshed::write32(fp, encodedSize);   // write length of encoded frame
        if ((int)fwrite(encoded, 1, encodedSize, fp) != encodedSize)
        {
          vvDebugMsg::msg(1, "Error: Cannot write voxel data to file.");
          fclose(fp);
          delete[] encoded;
          return FILE_ERROR;
        }
      }
      else  // no compression possible -> store unencoded
      {
        vvToolshed::write32(fp, 0);   // write zero to mark as unencoded
        if ((int)fwrite(raw, 1, frameSize, fp) != frameSize)
        {
          vvDebugMsg::msg(1, "Error: Cannot write voxel data to file.");
          fclose(fp);
          delete[] encoded;
          return FILE_ERROR;
        }
      }
    }
    else    // no compression
    {
      if ((int)fwrite(raw, 1, frameSize, fp) != frameSize)
      {
        vvDebugMsg::msg(1, "Error: Cannot write voxel data to file.");
        fclose(fp);
        delete[] encoded;
        return FILE_ERROR;
      }
    }
  }
  delete[] encoded;

  // Write transfer function(s):
  for (i=index=0; i<numTF; ++i)
  {
    if (i==0) pins = &vd->tf.pins;
    else 
    {
      pins = &vd->tf.presets[index];
      ++index;
    }

    if (!pins->isEmpty())
    {
      // Write TF name:
      sprintf(tfName, "TF %d", i+1);
      fputs(tfName, fp);    // default description string
      fputc(0, fp);         // string must be zero-terminated

      // Write pins:
      writePinList(fp, pins);
      for (j=0; j<vvPin::ELEMENTS_PER_PIN; ++j)
        vvToolshed::writeFloat(fp, -1.0f);
    }
  }

  // Clean up:
  fclose(fp);
  return OK;
}

//----------------------------------------------------------------------------
/** Save volume data to a Nrrd file (Gordon Kindlmann's proprietary format).
  See http://www.cs.utah.edu/~gk/teem/nrrd/ for more information.
  This file format cannot save transfer functions in the Virvo format.
*/
vvFileIO::ErrorType vvFileIO::saveNrrdFile(vvVolDesc* vd)
{
  char buf[256];                          // text buffer
  FILE* fp;                               // volume file pointer
  uchar* raw;                             // raw volume data
  int frameSize;                          // frame size
    
  vvDebugMsg::msg(1, "vvFileIO::saveNrrdFile()");
  
  if (vd->bpv>2) 
  {
    cerr << "Can only save 8 and 16 bit per voxel data in nrrd format." << endl;
    return FORMAT_ERROR;
  }

  // Prepare variables:
  if (vd->frames==0) return VD_ERROR;
  if (vd->frames>1) 
    cerr << "The nrrd writer will only write the first animation frame." << endl;
  frameSize = vd->getFrameSize();

  // Open file:
  if ( (fp = fopen(vd->getFilename(), "wb")) == NULL) // now open file to write
  {
    vvDebugMsg::msg(1, "Error: Cannot open file to write.");
    return FILE_ERROR;
  }

  // Write magic string:
  fputs(nrrdID, fp);
  fputc('\n', fp);

  // Write content:
  vvToolshed::extractBasename(buf, vd->getFilename());  
  fprintf(fp, "content: %s\n", buf);
  
  // Write data type:
  fputs("type: ", fp);
  switch (vd->bpv)
  {
    case 1: fputs("unsigned char\n", fp); break;
    case 2: fputs("unsigned short\n", fp); break;
    default: assert(0); break;
  }
  
  // Write dimension:
  fputs("dimension: 3\n", fp);

  // Write sizes:
  fprintf(fp, "sizes: %d %d %d\n", vd->vox[0], vd->vox[1], vd->vox[2]);
  
  // Write encoding:
  fputs("encoding: raw\n", fp);
 
  // Write spacings:
  fprintf(fp, "spacings: %f %f %f\n", vd->dist[0], vd->dist[1], vd->dist[2]);
  
  // Write endianness:
  fputs("endian: big\n", fp);

  // Write scalar data:
  fputc('\n', fp);
  raw = vd->getRaw(0);
  if (raw==NULL) 
  { 
    cerr << "Error: no data available for frame 0" << endl;
    fclose(fp);
    return VD_ERROR;
  }
  if ((int)fwrite(raw, 1, frameSize, fp) != frameSize)
  {
    vvDebugMsg::msg(1, "Error: Cannot write voxel data to file.");
    fclose(fp);
    return FILE_ERROR;
  }

  // Clean up:
  fclose(fp);
  return OK;
}

//----------------------------------------------------------------------------
/** Loader for voxel file in xvf (extended volume file) format.
 File format: see saveXVFFile()
*/
vvFileIO::ErrorType vvFileIO::loadXVFFile(vvVolDesc* vd)
{
  FILE* fp;             // volume file pointer
  vvPinList* pins;      // pin list
  float v[vvPin::ELEMENTS_PER_PIN];           // pin components
  uchar serialized[vvVolDesc::SERIAL_ATTRIB_SIZE];  // space for serialized volume data
  char tfName[257];     // transfer function name
  int f;                // counter for frames
  int c, i;             // counters
  int frameSize;        // size of a frame in bytes
  uchar* raw;           // raw volume data
  uchar* encoded = NULL;// encoded volume data
  int headerSize;       // total header size in bytes, including ID string
  int ctype;            // compression type
  int tnum;             // number of transfer functions
  int ttype;            // type of transfer function
  vvPin::PinType type;  // pin type
  bool done;
  int serializedSize;   // size of serialized part of header
  int encodedSize;      // size of encoded data array
  int offset;           // byte offset into file
   
  vvDebugMsg::msg(1, "vvFileIO::loadXVFFile()");

  if (vd->getFilename()==NULL) return FILE_ERROR;
	
  if ( (fp = fopen(vd->getFilename(), "rb")) == NULL) 
  {
    vvDebugMsg::msg(1, "Error: Cannot open file.");
    return FILE_ERROR;
  }

  vd->removeSequence();   // delete previous volume sequence

  // Read header:
  for (i=0; i<(int)strlen(xvfID); ++i)
  {
    if (fgetc(fp) != xvfID[i])
    {
      vvDebugMsg::msg(1, "Error: Invalid file ID string. Expected: ", xvfID);
      fclose(fp);
      return DATA_ERROR;
    }
  } 
  headerSize = vvToolshed::read16(fp);   // total header size in bytes
  if (headerSize<=28)   // early files didn't have transfer function info
    serializedSize = headerSize - strlen(xvfID);
  else
    serializedSize = headerSize - strlen(xvfID) - 7;

  fread(serialized, serializedSize, 1, fp);
  vd->deserializeAttributes(serialized, serializedSize);

  // Allow either bit or byte per voxel in vd->bpv:
  if (vd->bpv==8 || vd->bpv==16 || vd->bpv==24 || vd->bpv==32)
    vd->bpv /= 8;
	assert(vd->bpv>=1 && vd->bpv<=4);

  if (headerSize <= 28)   // early file format?
  {
    ctype = 0;
    tnum  = 0;
    ttype = 0;
  }
  else
  {
    ctype  = vvToolshed::read8(fp);
    tnum   = vvToolshed::read16(fp);
    ttype  = vvToolshed::read16(fp);
  }
  frameSize = vd->getFrameSize();

  // Print file information:
  if (vvDebugMsg::isActive(1))
  {  
    cerr << "XVF header size:                  " << headerSize << endl;
    cerr << "XVF compression type:             " << ctype << endl;
    cerr << "XVF number of transfer functions: " << tnum << endl;
    cerr << "XVF type of transfer function(s): " << ttype << endl;
  }

  // Load volume data:
  if ((sections & RAW_DATA) != 0)
  {
    fseek(fp, headerSize, SEEK_SET);
    if (ctype==1) encoded = new uchar[frameSize];
    for (f=0; f<vd->frames; ++f)
    {
      raw = new uchar[frameSize];   // create new data space for volume data
      switch (ctype)
      {
        default:
        case 0: // no compression
          if ((int)fread(raw, 1, frameSize, fp) != frameSize)
          {
            vvDebugMsg::msg(1, "Error: Insuffient voxel data in file.");
            fclose(fp);
            delete[] raw;
            delete[] encoded;
            return DATA_ERROR;
          }
          break;
        case 1: // RLE encoding
          encodedSize = vvToolshed::read32(fp);
          if (encodedSize>0)
          {
            if ((int)fread(encoded, 1, encodedSize, fp) != encodedSize)
            {
              vvDebugMsg::msg(1, "Error: Insuffient voxel data in file.");
              fclose(fp);
              delete[] raw;
              delete[] encoded;
              return DATA_ERROR;
            }
            if (vvToolshed::decodeRLE(raw, encoded, encodedSize, vd->bpv, frameSize) < 0)
            {
              vvDebugMsg::msg(1, "Error: Decoding exceeds frame size.");
              fclose(fp);
              delete[] raw;
              delete[] encoded;
              return DATA_ERROR;
            }
          } 
          else  // no encoding
          {
            if ((int)fread(raw, 1, frameSize, fp) != frameSize)
            {
              vvDebugMsg::msg(1, "Error: Insuffient voxel data in file.");
              fclose(fp);
              delete[] raw;
              delete[] encoded;
              return DATA_ERROR;
            }
          }
          break;          
      }
      vd->addFrame(raw, vvVolDesc::DELETE_DATA);
    }
    delete[] encoded;
  }

  // Read transfer function(s):
  if ((sections & TRANSFER) != 0)
  {
    offset = headerSize;
    fseek(fp, offset, SEEK_SET);
    for (f=0; f<vd->frames; ++f)
    {
      switch (ctype)
      {
        default:
        case 0: offset += frameSize;
                break;
        case 1: encodedSize = vvToolshed::read32(fp);
                offset += 4;    // add 4 bytes for encoded size value
                if (encodedSize==0) offset += frameSize;
                else offset += encodedSize;
                break;
      }
      fseek(fp, offset, SEEK_SET);
    }
    for (i=0; i<tnum; ++i)
    { 
      // Read zero terminated TF name:
      c = 0;
      do 
      {
        tfName[c] = char(fgetc(fp));
        if (feof(fp)) break;
        ++c;
      } while (tfName[c-1] != 0);

      if (!feof(fp))
      {
        if (ttype != 1) break;   // only accept pin lists

        if (i==0) pins = &vd->tf.pins;
        else      pins = &vd->tf.presets[i-1];
        pins->clear();

        done = false;
        while (done==false)
        {
          for (c=0; c<vvPin::ELEMENTS_PER_PIN; ++c)
            v[c] = vvToolshed::readFloat(fp);
          if (v[0]==-1.0f || feof(fp)) done = true;
          else
          {
            type = (vvPin::PinType)((int)v[0]);
            pins->add(new vvPin(type, v[1], v[2], v[3], v[4]));
          }
        } 
      }
    }  
  }

  // Clean up:
  fclose(fp);
  return OK;
}

//----------------------------------------------------------------------------
/** Loader for voxel file in avf (ASCII volume file) format.
 File specification:
 <PRE>
 Header: 

   In the header, several lines give information about the data format.
   Each consists of an identifier and a value, separated by whitespace.
   Each line can contain one identifier and one value. 
   This file format cannot store transfer functions.
   At any point in the file comments starting with '#' are allowed, they
   cover the rest of the current line.

   The following abbreviations are used:  
   <int>            for integer values
   <float>          for floating point values
   <OPT1|OPT2|OPT3> for an options list

   The following lines are required:
   WIDTH <int>      the width of the volume [voxels]
   HEIGHT <int>     the height of the volume [voxels]
   SLICES <int>     the number of slices in the volume [voxels]
   
   The following lines are optional. 
   If they are missing, default values are used:
   FRAMES <int>     the number of data sets contained in the file
                    (default: 1)
   MIN <float>      the minimum data value, smaller values will be constrained
                    to this value (default: 0.0)
   MAX <float>      the maximum data value, larger values will be constrained
                    to this value (default: 1.0)
   FORMAT <SCALAR8|SCALAR16|RGB|RGBA>  voxel data format (default: SCALAR8):                  
                    SCALAR8  = scalar values, passed as 1 number, stored as 8 bit integers
                    SCALAR16 = scalar values, passed as 1 number, stored as 16 bit integers
                    RGB      = color values, consisting of a red, a green, and 
                               a blue color component, passed as 3 numbers in a row,
                               stored as 3x8 bit
                    RGBA     = color values, consisting of a red, a green, 
                               a blue, and an opacity (alpha) value,
                               passed as 4 numbers in a row, stored as 4x8 bit
   XDIST <float>    the sample distance in x direction (-> width) [mm] 
                    (default: 1.0)
   YDIST <float>    the sample distance in y direction (-> height) [mm]
                    (default: 1.0)
   ZDIST <float>    the sample distance in z direction (-> slices) [mm]
                    (default: 1.0)
   TIME <float>     the length of each time step for transient data [s]
                    (default: 1.0)

 Data area:         
   
   The data area starts right after the header. The voxel data values
   are listed, separated by whitespace and/or end-of-line markers.
   Both float and integer values are accepted.
   Voxel order: voxel at top left front first, then to right, 
   then to bottom, then to back, then frames. 
   All elements of each voxel (RGB, RGBA) are stored consecutively.

 Sample file:

 WIDTH   4
 HEIGHT  3
 SLICES  2
 FRAMES  1
 MIN     0.0
 MAX     1.0
 FORMAT  SCALAR8    # 8 bit data
 XDIST   1.0
 YDIST   1.0
 ZDIST   1.0
 TIME    1.0
 0.9 0.9 0.9 0.9
 0.9 0.2 0.3 0.9
 0.9 0.2 0.4 0.9
 0.8 0.8 0.8 0.8
 0.8 0.1 0.1 0.8
 0.8 0.0 0.0 0.8
</PRE>
*/
vvFileIO::ErrorType vvFileIO::loadAVFFile(vvVolDesc* vd)
{
  vvTokenizer* tokenizer; // ASCII file tokenizer
  vvTokenizer::TokenType ttype;  // currently processed token type
  FILE* fp;               // volume file pointer
  uchar* raw;             // raw volume data
  float cval;             // data value constrained to min/max
  int   vval;             // data value in raw data range
  int f;                  // counter for frames
  int i, x, y, z, v;      // counters
  int frameSize;          // size of a frame in bytes
  int identifier;         // ID of string identifier in file header
  int numPerVoxel;        // numbers per voxel in ASCII file(8/16bit: 1, 24bit: 3, 32bit: 4)
  bool done;          
  bool error;
   
  vvDebugMsg::msg(1, "vvFileIO::loadAVFFile()");

  if (vd->getFilename()==NULL) return FILE_ERROR;
	
  if ( (fp = fopen(vd->getFilename(), "rb")) == NULL) 
  {
    vvDebugMsg::msg(1, "Error: Cannot open file.");
    return FILE_ERROR;
  }

  // Delete previous volume sequence
  vd->removeSequence();   

  // Set default values:
  vd->vox[0]  = 0;
  vd->vox[1]  = 0;
  vd->vox[2]  = 0;
  vd->frames  = 1;
  vd->bpv     = 1;
  vd->dist[0] = 1.0f;
  vd->dist[1] = 1.0f;
  vd->dist[2] = 1.0f;
  vd->dt      = 1.0f;
  vd->realMin = 0.0f;
  vd->realMax = 1.0f;

  // Read header data:
  tokenizer = new vvTokenizer(fp);
  tokenizer->setCommentCharacter('#');
  tokenizer->setEOLisSignificant(false);
  tokenizer->setCaseConversion(vvTokenizer::VV_UPPER);
  tokenizer->setParseNumbers(true);
  tokenizer->setWhitespaceCharacter('=');
  done = error = false;
  while (!done)
  {
    // Read identifier:
    ttype = tokenizer->nextToken();
    if (ttype != vvTokenizer::VV_WORD)
    {
      done = true;
      tokenizer->pushBack();
      continue;
    }
    else if (strcmp(tokenizer->sval, "WIDTH")==0)
      identifier = 0;
    else if (strcmp(tokenizer->sval, "HEIGHT")==0)
      identifier = 1;
    else if (strcmp(tokenizer->sval, "SLICES")==0)
      identifier = 2;
    else if (strcmp(tokenizer->sval, "FRAMES")==0)
      identifier = 3;
    else if (strcmp(tokenizer->sval, "MIN")==0)
      identifier = 4;
    else if (strcmp(tokenizer->sval, "MAX")==0)
      identifier = 5;
    else if (strcmp(tokenizer->sval, "XDIST")==0)
      identifier = 6;
    else if (strcmp(tokenizer->sval, "YDIST")==0)
      identifier = 7;
    else if (strcmp(tokenizer->sval, "ZDIST")==0)
      identifier = 8;
    else if (strcmp(tokenizer->sval, "TIME")==0)
      identifier = 9;
    else if (strcmp(tokenizer->sval, "FORMAT")==0)
      identifier = 10;
    else 
    {
      done = error = true;
      continue;
    }
    
    // Read assigned value:
    ttype = tokenizer->nextToken();
    if (ttype == vvTokenizer::VV_NUMBER)
    {
      switch (identifier)
      {
        case 0: vd->vox[0] = (int)tokenizer->nval; break;
        case 1: vd->vox[1] = (int)tokenizer->nval;break;
        case 2: vd->vox[2] = (int)tokenizer->nval;break;
        case 3: vd->frames = (int)tokenizer->nval;break;
        case 4: vd->realMin = tokenizer->nval; break;
        case 5: vd->realMax = tokenizer->nval; break;
        case 6: vd->dist[0] = tokenizer->nval;break;
        case 7: vd->dist[1] = tokenizer->nval;break;
        case 8: vd->dist[2] = tokenizer->nval;break;
        case 9: vd->dt = tokenizer->nval;break;
        default: break;
      }        
    }
    else if (ttype == vvTokenizer::VV_WORD)
    {
      if (identifier == 10)
      {
        if (strcmp(tokenizer->sval, "SCALAR8")==0)
          vd->bpv = 1;
        else if (strcmp(tokenizer->sval, "SCALAR16")==0)
          vd->bpv = 2;
        else if (strcmp(tokenizer->sval, "RGB")==0)
          vd->bpv = 3;
        else if (strcmp(tokenizer->sval, "RGBA")==0)
          vd->bpv = 4;
        else 
          error = done = true;
      }
    }
    else
      error = done = true;
  }  
  if (error)
  {
    cerr << "Read error in line " << tokenizer->getLineNumber() << " of file " <<
      vd->getFilename() << endl;
    delete tokenizer;
    fclose(fp);
    return DATA_ERROR;
  }
  
  // Check for consistence:
  if (vd->vox[0]<=0 || vd->vox[1]<=0 || vd->vox[2]<=0 || 
    vd->frames<=0 || vd->realMin>=vd->realMax)
  {
    vvDebugMsg::msg(1, "Error: Invalid file information in header");
    delete tokenizer;
    fclose(fp);
    return DATA_ERROR;
  }

  // Load voxel data:
  numPerVoxel = (vd->bpv == 2) ? 1 : vd->bpv;
  frameSize = vd->getFrameSize();
  if ((sections & RAW_DATA) != 0)
  {
    for (f=0; f<vd->frames; ++f)
    {
      raw = new uchar[frameSize];   // create new data space for volume data
      i = 0;
      for (z=0; z<vd->vox[2]; ++z)
        for (y=0; y<vd->vox[1]; ++y)
          for (x=0; x<vd->vox[0]; ++x)
            for (v=0; v<numPerVoxel; ++v)
            {
              ttype = tokenizer->nextToken();
              if (ttype != vvTokenizer::VV_NUMBER)
              {
                cerr << "Syntax error in line " << tokenizer->getLineNumber() << endl;
                delete tokenizer;
                fclose(fp);
                return DATA_ERROR;
              }
              cval = (tokenizer->nval - vd->realMin) / (vd->realMax - vd->realMin);
              cval = ts_clamp(cval, 0.0f, 1.0f);   // constrain value to 0..1
              if (vd->bpv==2)     // 16 bit value?
              {
                vval = int(65535.0f * cval);
                raw[i++] = (uchar)(vval >> 8); 
                raw[i++] = (uchar)(vval & 0xFF);
              }
              else raw[i++] = (uchar)(255.0f * cval);
            }
      vd->addFrame(raw, vvVolDesc::DELETE_DATA);
    }
  }

  // Clean up:
  delete tokenizer;
  fclose(fp);
  return OK;
}

//----------------------------------------------------------------------------
/** Saver for vf file format
  @return error code, OK if successful
*/
vvFileIO::ErrorType vvFileIO::saveVFFile(vvVolDesc* vd)
{
  vvvffile* vff;
  int i,j;
  int code, para, lengt;
  unsigned char* coded_data = NULL;
  vvPin* curPin;      // current pin
  
  vvDebugMsg::msg(1, "vvFileIO::saveVFFile()");

  vff = new vvvffile((char*)vd->getFilename());
  vff->setDataFile();
  vff->setArrayDimensions(vd->vox[0], vd->vox[1], vd->vox[2], vd->bpv);
  vff->setVoxelsize(vd->dist[0], vd->dist[1], vd->dist[2]);
  for (i=0; i<vd->frames; ++i)
  {
    code = vff->findBestTypeOfCoding(vd->getRaw(i), &para, &lengt);
    coded_data = new uchar[lengt];
    vff->encodeData(vd->getRaw(i), code, para, coded_data);
    vff->addToDataFile(lengt, (uchar)code, coded_data);
    delete[] coded_data; 
  }

  vff->initPins();
  
  // Add current pin list:
  for (i=0; i<vd->tf.pins.count(); ++i)
  {
    curPin = vd->tf.pins.getPin(i);
    vff->addPin(1, curPin->type + 1, curPin->v[0], curPin->v[1], 
      curPin->v[2], curPin->x, NULL);
  }

  // Add preset pin lists:
  for (j=0; j<vd->tf.getNumPresets(); ++j)
    for (i=0; i<vd->tf.presets[j].count(); ++i)
    {
      curPin = vd->tf.presets[j].getPin(i);
      vff->addPin(j+2, curPin->type + 1, curPin->v[0], curPin->v[1], 
        curPin->v[2], curPin->x, NULL);
    }

  vff->setValues(220, 30, vd->frames, (char*)vd->getFilename()); 
  vff->setArrayDimensions(vd->vox[0], vd->vox[1], vd->vox[2], vd->bpv);
  vff->setHotspot(vd->vox[0] / 2, vd->vox[1] / 2, vd->vox[2] / 2);
  vff->setFileDescription("n/a");
  if (vd->getIconSize() > 0)
    vff->setIcon(2, vd->getIconData());
  else
    vff->setIcon(0, NULL);
  vff->setVoxelsize(vd->dist[0], vd->dist[1], vd->dist[2]);
  vff->writeFile();    
  delete vff;
  return OK;
}

//----------------------------------------------------------------------------
/** Loader for xb7 files. These files are used in SFB 382, project C15 (Prof.
  Herrmann, Stefan Luding, Stefan Miller, Sean McNamara).<br> 
  xb7-Format:<br>
  Series of snapshots at times t in blocks of N+1 lines each (header + N particles) 
  with 8 numbers, separated by whitespace.<p>
  
  line 1:  N  t  x1 y1 z1  x2 y2 z2 <br>
  line 2 - line N+1:  x y z  vx vy vz  r  i <p>

  N: number of particles<br>
  t: time<br>
  x1/y1/z1-x2/y2/z2  size of simulation box (hier x1=y1=z1=0, x2=y2=z2)<br>
  x/y/z: coordinates (particle center is within box, periodic boundaries)<br>
  vx/vy/vz: speed<br>
  r: diameter<br>  
  i: collision frequency<p> 

  An example file with 3 particles is:
  <pre>
  3 40.9594 0 0 0 0.102355 0.102355 0.102355	
  0.0914023 0.0886842 0.0880599	-0.000187777 -4.58716e-05 -0.000202447	0.0005	219 
  0.0183272 0.0727637 0.0348822	4.57354e-05 -0.000339601 0.000512404	0.0005	259 
  0.0955405 0.0885498 0.00429593	-0.000176341 -0.000405909 -0.000278665	0.0005	1487 
  </pre>
  For time dependent data, the voxel grid is sized to the simulation box of the first
  time step. In subsequent time steps particles which are not within the first box will
  be discarded. Also, the value range fitting will be done only with the particles of the
  first time step.
  @param vd volume description
  @param maxEdgeLength volume will be shaped like the simulation box with this maximum edge length [voxels]
  @param densityParam parameter to use for voxel density: 0=x, 1=y, 2=z, 3=vx, 4=vy, 5=vz, 6=r, 7=i, 8=sqrt(vx*vx+vy*vy+vz*vz)
  @param useGlobalMinMax true if min and max scalar values should be determined globally, false for within time steps only
*/
vvFileIO::ErrorType vvFileIO::loadXB7File(vvVolDesc* vd, int maxEdgeLength, int densityParam, bool useGlobalMinMax)
{
  vvTokenizer* tokenizer;         // ASCII file tokenizer
  vvTokenizer::TokenType ttype;   // currently processed token type
  vvSLList<ParticleTimestep*> timesteps;  // particle storage for all time steps
  FILE* fp;                       // volume file pointer
  uchar* raw;                     // raw volume data
  float boxMin[3];                // simulation box min values
  float boxMax[3];                // simulation box max values
  float boxSize[3];               // simulation box size
  float maxBoxSize;               // maximum box size
  float globalMin,globalMax;      // real min and max values over all time steps
  float minVal,maxVal;            // real min and max values
  float param[9];                 // atom parameters of one line
  float val;                      // particle value
  int numParticles=0;             // number of particles in current time step
  int numTimesteps;               // number of time steps in file
  int frameSize;                  // number of bytes per frame
  int iVal;                       // integer density
  int iPos[3];                    // position of particle in volume
  int index;
  int i,j,t;
  bool error = false;
  
  vvDebugMsg::msg(1, "vvFileIO::loadXB7File()");

  if (vvToolshed::isFile(vd->getFilename())==false)
    return FILE_NOT_FOUND;
	
  if ((fp = fopen(vd->getFilename(), "rb")) == NULL) 
  {
    vvDebugMsg::msg(1, "Error: Cannot open file.");
    return FILE_ERROR;
  }

  cerr << "Reading XB7 file: edge length=" << maxEdgeLength << ", density parameter=" << densityParam <<
    ", global min/max=" << useGlobalMinMax << endl;

  // Delete previous volume sequence
  vd->removeSequence();   

  // Set default values:
  vd->bpv = 2;
  vd->dt = 0.1f;

  // Initialize tokenizer:
  tokenizer = new vvTokenizer(fp);
  tokenizer->setEOLisSignificant(true);
  tokenizer->setCaseConversion(vvTokenizer::VV_LOWER);
  tokenizer->setParseNumbers(true);

  // Read all time step data but don't create volume yet:
  for (i=0; i<3; ++i)
  {
    boxMin[i] =  VV_FLT_MAX;
    boxMax[i] = -VV_FLT_MAX;
  }
  globalMin =  VV_FLT_MAX; 
  globalMax = -VV_FLT_MAX;
  for(;;)
  {
    // Parse header:  
    for(i=0; i<8; ++i)
    {
      // Read identifier:
      ttype = tokenizer->nextToken();
      if (ttype != vvTokenizer::VV_NUMBER)
      {
        error = true;
        break;
      }
      switch (i)
      {
        case 0: numParticles = int(tokenizer->nval); break;
        case 1: break;    // ignore time value
        case 2: if (tokenizer->nval < boxMin[0]) boxMin[0] = tokenizer->nval; break;
        case 3: if (tokenizer->nval < boxMin[1]) boxMin[1] = tokenizer->nval; break;
        case 4: if (tokenizer->nval < boxMin[2]) boxMin[2] = tokenizer->nval; break;
        case 5: if (tokenizer->nval > boxMax[0]) boxMax[0] = tokenizer->nval; break;
        case 6: if (tokenizer->nval > boxMax[1]) boxMax[1] = tokenizer->nval; break;
        case 7: if (tokenizer->nval > boxMax[2]) boxMax[2] = tokenizer->nval; break;
        default: break;
      }        
    }
    ttype = tokenizer->nextToken();
    if (ttype!=vvTokenizer::VV_EOL) error = true;
    if (error)
    {
      cerr << "Parse error in line " << tokenizer->getLineNumber() << endl;
      delete tokenizer;
      timesteps.removeAll();
      fclose(fp);
      return DATA_ERROR;
    }
    
    // Create new time step:
    timesteps.append(new ParticleTimestep(numParticles), true);  
    cerr << "Reading " << numParticles << " particles" << endl;
    
    timesteps.getData()->min = VV_FLT_MAX;
    timesteps.getData()->max = -VV_FLT_MAX;

    // Load particles, but don't create the volume yet:
    for (i=0; i<numParticles; ++i)
    {
      // Read an entire line including new line character:
      for (j=0; j<9; ++j)
      {
        if (i==numParticles-1 && j==8) continue;  // ignore last EOL
        ttype = tokenizer->nextToken();
        if ((j<8 && ttype != vvTokenizer::VV_NUMBER) || (j==8 && ttype!=vvTokenizer::VV_EOL))
        {
          cerr << "Parse error in line " << tokenizer->getLineNumber() << endl;
          delete tokenizer;
          timesteps.removeAll();
          fclose(fp);
          return DATA_ERROR;
        }
        param[j] = tokenizer->nval;
      }
      for (j=0; j<3; ++j)   // memorize particle position
      {
        timesteps.getData()->pos[j][i] = param[j];
      }
      if (densityParam<8) // memorize particle value
      {
        val = param[densityParam];
      }
      else  // speed
      {
        val = sqrtf(param[3] * param[3] + param[4] * param[4] + param[5] * param[5]);
      }
      timesteps.getData()->val[i] = val;
      if (val < timesteps.getData()->min) timesteps.getData()->min = val;
      if (val > timesteps.getData()->max) timesteps.getData()->max = val;
    }
    cerr << "Timestep: scalar min,max: " << timesteps.getData()->min << "," << timesteps.getData()->max << endl;
    if (timesteps.getData()->min < globalMin) globalMin = timesteps.getData()->min;
    if (timesteps.getData()->max > globalMax) globalMax = timesteps.getData()->max;
    
    // Look for another time step:
    do
    {
      ttype = tokenizer->nextToken();
    } while (ttype != vvTokenizer::VV_EOF && ttype != vvTokenizer::VV_NUMBER);
    if (ttype==vvTokenizer::VV_EOF) break;
    else tokenizer->pushBack();
  }
  delete tokenizer;
  fclose(fp);
  numTimesteps = timesteps.count();
  cerr << numTimesteps << " time steps read" << endl;
  cerr << "Global min,max: " << globalMin << "," << globalMax << endl;

  // Check for consistency:
  if (boxMin[0] >= boxMax[0] || boxMin[1] >= boxMax[1] || boxMin[2] >= boxMax[2])
  {
    cerr << "Error: invalid box size in header: " << boxMin[0] << " " << boxMin[1] << " " << boxMin[2] << 
      " " << boxMax[0] << " " << boxMax[1] << " " << boxMax[2] << endl;
    timesteps.removeAll();
    return DATA_ERROR;
  }
  vd->realMin = globalMin;
  vd->realMax = globalMax;

  // Now that all particles are read from all time steps, the volumes can be generated.
  
  // Compute header values:
  for(i=0; i<3; ++i)
  {
    boxSize[i] = boxMax[i] - boxMin[i];
  }
  maxBoxSize = ts_max(boxSize[0], boxSize[1], boxSize[2]);
  for(i=0; i<3; ++i)
  {
    vd->vox[i] = int(float(maxEdgeLength) * boxSize[i] / maxBoxSize);
    vd->vox[i] = ts_clamp(vd->vox[i], 1, maxEdgeLength);
  }
  for (i=0; i<3; ++i)
  {
    vd->dist[i] = boxSize[i] / float(vd->vox[i]);
  }

  frameSize = vd->getFrameSize();
  timesteps.first();
  for(t=0; t<numTimesteps; ++t)
  {
    raw = new uchar[frameSize];
    assert(raw);
    memset(raw, 0, frameSize);
    numParticles = timesteps.getData()->numParticles;

    for (i=0; i<numParticles; ++i)
    {
      for (j=0; j<3; ++j)
      {
        iPos[j] = int(float(vd->vox[j] - 1) * (timesteps.getData()->pos[j][i] - boxMin[j]) / (boxMax[j] - boxMin[j]));
        iPos[j] = ts_clamp(iPos[j], 0, vd->vox[j] - 1);
      }
      // Allow values from 1 to MAX_16BIT:
      if (useGlobalMinMax)
      {
        minVal = globalMin;
        maxVal = globalMax;
      }
      else
      {
        minVal = timesteps.getData()->min;
        maxVal = timesteps.getData()->max;
      }
      if (maxVal > minVal) iVal = int(65534.0f * (timesteps.getData()->val[i] - minVal) / (maxVal - minVal)) + 1;
      else iVal = 65535;
      iVal = ts_clamp(iVal, 1, 65535);
      index = 2 * (iPos[0] + iPos[1] * vd->vox[0] + iPos[2] * vd->vox[0] * vd->vox[1]);
      raw[index] = uchar(iVal >> 8);
      raw[index + 1] = uchar(iVal & 0xff);
    }
    vd->addFrame(raw, vvVolDesc::DELETE_DATA);
    timesteps.next();
  }

  vd->frames = vd->getStoredFrames();
  assert(vd->frames == timesteps.count());
  timesteps.removeAll();
  return OK;
}

//----------------------------------------------------------------------------
/// Loader for vf and xb7 file formats.
vvFileIO::ErrorType vvFileIO::loadVFFile(vvVolDesc* vd)
{
  vvvffile* vff;
  uchar* raw;   // raw volume data
  int i,j;
  uchar* iconData;
  int iconSize; // icon size [bytes]
  int iconEdge; // icon edge length [pixels]
  vvPin* newPin;
  int result;
  int type, description_length;
  float float1, float2, float3, xposition;
  bool allDone;
  vvPin::PinType pinType;    // pin type

  vvDebugMsg::msg(1, "vvFileIO::loadVFFile()");
  if (vd->getFilename()==NULL) return FILE_ERROR;

  vd->removeSequence();

  vff = new vvvffile((char*)vd->getFilename());
  vd->frames = vff->getNumberOfDataArrays();
  vff->getArrayDimensions(&vd->vox[0], &vd->vox[1], &vd->vox[2], &vd->bpv);
  vff->getVoxelsize(&vd->dist[0], &vd->dist[1], &vd->dist[2]);

  // Read icon from file:
  if ((sections & ICON) != 0)
  {
    iconSize = vff->getIconSize();
    iconData = new uchar[iconSize];
    vff->getIcon(iconData);
    if (iconData!=NULL)
    {
      iconEdge = (int)sqrt((double)(iconSize/3));
      vd->makeIcon(iconEdge, iconData);
    }
  }

  // Read raw data arrays from file:
  if ((sections & RAW_DATA) != 0)
  {
    // Safety check:
    if (vd->getFrameSize()==0 || vd->frames==0)
    {
      delete vff;
      return vvFileIO::DATA_ERROR;
    }

    for (i=0; i<vd->frames; ++i)
    {
      // Create space for one frame:
      raw = new uchar[vd->getFrameSize()];

      // Read data array:
      vff->readDataArray(&raw, i+1);
      vd->addFrame(raw, vvVolDesc::DELETE_DATA);
    }
  }

  // Read transfer functions:
  if ((sections & TRANSFER) != 0)
  {
    // Add pins to current list:
    vd->tf.pins.clear();
    for (i=0;; ++i)
    {
      result = vff->getPin(1, i+1, &type, &float1, &float2, &float3, &xposition, &description_length);
      if (result==1)
        break;
      pinType = (vvPin::PinType)(type-1);
      newPin = new vvPin(pinType, xposition, float1, float2, float3);
      vd->tf.pins.add(newPin);
    }

    // Add to preset pin lists:
    allDone = false;
    for (j=0; allDone==false; ++j)
    {
      vd->tf.presets[j].clear();
      for (i=0;; ++i)
      {
        result = vff->getPin(j+1, i+1, &type, &float1, &float2, &float3, &xposition, &description_length);
        if (result==1)
        {
          if (i==0) allDone = true;
          break;
        }
        pinType = (vvPin::PinType)(type-1);
        newPin = new vvPin(pinType, xposition, float1, float2, float3);
        vd->tf.presets[j].add(newPin);
      } 
    } 
  }

  // Clean up:
  delete vff;
  return OK;
}

//----------------------------------------------------------------------------
// TODO: parse header!
/** Loader for IMD checkpoint files with particle data. 
  These files are used at ITAP in SFB 382, project C14 
  (Prof. Trebin, Gunther Schaaf, Franz Gaehler, Silvia Hocker).<br> 
  cpt-Format:<br>
  This is an example file with three atoms:
  <PRE>
  #F A 1 1 1 3 3 1
  #C number type mass x y z vx vy vz Epot
  #X      3.0766609395000000e+02 0.0000000000000000e+00 0.0000000000000000e+00
  #Y      0.0000000000000000e+00 1.0442535916900000e+02 0.0000000000000000e+00
  #Z      0.0000000000000000e+00 0.0000000000000000e+00 1.4357751050999990e+01
  ## Generated by /hwwt3e/rus/ita/pof30/bin/cray-t3e/imd_mpi_nve_stress_ordpar_efilter on Thu Jun 27 23:56:46 2002
  #E
  18368 0     1.000000    12.311201    48.337746     1.031030    -0.032552     0.047432    -0.014428   -17.623187
  18800 0     1.000000    12.310159    48.341527     3.080848    -0.024594     0.040695    -0.009033   -17.630691
  15772 1     1.000000    10.766565    47.946747     2.054420    -0.009312    -0.063240    -0.027128   -21.059210
  </PRE>
  For time dependent data, multiple files are stored on disk with increasing filename numbers, 
  e.g., timestep001.cpt, timestep002.cpt, etc.
  @param vd volume description
  @param maxEdgeLength volume will be shaped like the simulation box with this maximum edge length [voxels]
  @param densityParam value index to use for voxel density, use -1 for speed (vx*vx+vy*vy+vz*vz)
  @param useGlobalMinMax true if min and max scalar values should be determined globally, false for within time steps only
*/
vvFileIO::ErrorType vvFileIO::loadCPTFile(vvVolDesc* vd, int maxEdgeLength, int densityParam, bool useGlobalMinMax)
{
  vvTokenizer* tokenizer;           // ASCII file tokenizer
  vvTokenizer::TokenType ttype;     // currently processed token type
  vvSLList<ParticleTimestep*> timesteps; // particle storage for all time steps
  vvArray<float> particles;         // densities of current time step
  vvArray<float> xpos;              // x positions of current time step
  vvArray<float> ypos;              // y positions of current time step
  vvArray<float> zpos;              // z positions of current time step
  FILE* fp;                         // volume file pointer
  uchar* raw;                       // raw volume data
  char* filename;                   // current particles file name
  float boxMin[3];                  // simulation box min values
  float boxMax[3];                  // simulation box max values
  float boxSize[3];                 // simulation box size
  float maxBoxSize;                 // maximum box size
  float globalMin,globalMax;        // density min and max values over all time steps
  float minVal,maxVal;              // density min and max values of current time step
  float val;                        // particle value
  int numParticles=0;               // number of particles in current time step
  int numTimesteps;                 // number of time steps in file
  int frameSize;                    // number of bytes per frame
  int iVal;                         // integer density
  int iPos[3];                      // position of particle in volume
  int index;
  int i,j,t;
  float speed[3];
  
  vvDebugMsg::msg(1, "vvFileIO::loadCPTFile()");

  filename = new char[strlen(vd->getFilename()) + 1];
  strcpy(filename, vd->getFilename());

  if (vvToolshed::isFile(filename)==false)
    return FILE_NOT_FOUND;
  	
  cerr << "Checkpoint reader parameters: edge length=" << maxEdgeLength << ", density parameter=" << densityParam <<
    ", global min/max=" << useGlobalMinMax << endl;

  // Delete previous volume sequence
  vd->removeSequence();   

  // Set default values:
  vd->bpv = 2;
  vd->dt = 0.1f;

  // Initialize variables:
  for (i=0; i<3; ++i)
  {
    boxMin[i] =  VV_FLT_MAX;
    boxMax[i] = -VV_FLT_MAX;
  }
  globalMin =  VV_FLT_MAX; 
  globalMax = -VV_FLT_MAX;
  
  // Loop thru time steps:
  for(;;)   
  {
    if ((fp = fopen(filename, "rb")) == NULL) 
    {
      vvDebugMsg::msg(1, "Error: Cannot open file: ", filename);
      return FILE_ERROR;
    }

    // Initialize tokenizer:
    tokenizer = new vvTokenizer(fp);
    tokenizer->setEOLisSignificant(true);
    tokenizer->setCaseConversion(vvTokenizer::VV_LOWER);
    tokenizer->setParseNumbers(true);
    tokenizer->setCommentCharacter('#');

    particles.clear();
    xpos.clear();
    ypos.clear();
    zpos.clear();
    minVal =  VV_FLT_MAX;
    maxVal = -VV_FLT_MAX;
    
    // Load particles, but don't create the volume yet:
    for (;;)    // loop thru particles in one time step
    {
      // Read an entire line of numbers:
      for(i=0; tokenizer->nextToken() == vvTokenizer::VV_NUMBER; ++i)
      {
        // Memorize position and adjust simulation box:
        if (i>=0 && i<=2)
        {
          if (i==0) xpos.append(tokenizer->nval);
          else if (i==1) ypos.append(tokenizer->nval);
          else if (i==2) zpos.append(tokenizer->nval);
          if (tokenizer->nval < boxMin[i]) boxMin[i] = tokenizer->nval;       
          if (tokenizer->nval > boxMax[i]) boxMax[i] = tokenizer->nval;       
        }
        
        // Memorize density value:
        if (densityParam==-1)
        {
          if (i>=3 && i<=5) speed[i-3] = tokenizer->nval;
          if (i==5) val = sqrtf(speed[0] * speed[0] + speed[1] * speed[1] + speed[2] * speed[2]);
        }
        else if (i==densityParam) val = tokenizer->nval;
      }
      particles.append(val);
      if (val < minVal) minVal = val;
      if (val > maxVal) maxVal = val;

      // Look for another particle:
      do
      {
        ttype = tokenizer->nextToken();
      } while (ttype != vvTokenizer::VV_EOF && ttype != vvTokenizer::VV_NUMBER);
      if (ttype==vvTokenizer::VV_EOF) break;
      else tokenizer->pushBack();
    }
    delete tokenizer;
    fclose(fp);
    cerr << "Timestep: scalar min,max: " << minVal << "," << maxVal << endl;
    if (minVal < globalMin) globalMin = minVal;
    if (maxVal > globalMax) globalMax = maxVal;

    // Create new time step and copy data to it:    
    timesteps.append(new ParticleTimestep(particles.count()));
    timesteps.getData()->max = maxVal;
    timesteps.getData()->min = minVal;
    memcpy(timesteps.getData()->val, particles.getArrayPtr(), particles.count());
    memcpy(timesteps.getData()->pos[0], xpos.getArrayPtr(), xpos.count());
    memcpy(timesteps.getData()->pos[1], ypos.getArrayPtr(), ypos.count());
    memcpy(timesteps.getData()->pos[2], zpos.getArrayPtr(), zpos.count());
    
    // Look for another time step:
    if (!vvToolshed::increaseFilename(filename)) break;
    if (vvToolshed::isFile(filename)==false) break;
  }
  numTimesteps = timesteps.count();
  cerr << numTimesteps << " time steps read" << endl;
  cerr << "Global: scalar min,max: " << globalMin << "," << globalMax << endl;

  vd->realMin = globalMin;
  vd->realMax = globalMax;

  // Now that all particles are read from all time steps, the volumes can be generated.
  
  // Compute header values:
  for(i=0; i<3; ++i)
  {
    boxSize[i] = boxMax[i] - boxMin[i];
  }
  maxBoxSize = ts_max(boxSize[0], boxSize[1], boxSize[2]);
  for(i=0; i<3; ++i)
  {
    vd->vox[i] = int(float(maxEdgeLength) * boxSize[i] / maxBoxSize);
    vd->vox[i] = ts_clamp(vd->vox[i], 1, maxEdgeLength);
  }
  for (i=0; i<3; ++i)
  {
    vd->dist[i] = boxSize[i] / float(vd->vox[i]);
  }
  
  frameSize = vd->getFrameSize();
  timesteps.first();
  for(t=0; t<numTimesteps; ++t)
  {
    raw = new uchar[frameSize];
    assert(raw);
    memset(raw, 0, frameSize);
    numParticles = timesteps.getData()->numParticles;

    for (i=0; i<numParticles; ++i)
    {
      for (j=0; j<3; ++j)
      {
        iPos[j] = int(float(vd->vox[j] - 1) * (timesteps.getData()->pos[j][i] - boxMin[j]) / (boxMax[j] - boxMin[j]));
        iPos[j] = ts_clamp(iPos[j], 0, vd->vox[j] - 1);
      }
      // Allow values from 1 to MAX_16BIT:
      if (useGlobalMinMax)
      {
        minVal = globalMin;
        maxVal = globalMax;
      }
      else
      {
        minVal = timesteps.getData()->min;
        maxVal = timesteps.getData()->max;
      }
      if (maxVal > minVal) iVal = int(65534.0f * (timesteps.getData()->val[i] - minVal) / (maxVal - minVal)) + 1;
      else iVal = 65535;
      iVal = ts_clamp(iVal, 1, 65535);
      index = 2 * (iPos[0] + iPos[1] * vd->vox[0] + iPos[2] * vd->vox[0] * vd->vox[1]);
      raw[index] = uchar(iVal >> 8);
      raw[index + 1] = uchar(iVal & 0xff);
    }
    vd->addFrame(raw, vvVolDesc::DELETE_DATA);
    timesteps.next();
  }

  vd->frames = vd->getStoredFrames();
  assert(vd->frames == timesteps.count());
  timesteps.removeAll();
  delete[] filename;
  return OK;
}

//----------------------------------------------------------------------------
/// Loader for voxel file in tif (Tagged Image File) format.
vvFileIO::ErrorType vvFileIO::loadTIFFile(vvVolDesc* vd)
{
  const ushort ENDIANNESS = 19789;    // TIF endianness for big-endian (Unix) style
  const ushort MAGICNUMBER = 42;      // TIF magic number
  FILE* fp;             // volume file pointer
  ushort endian, magic; // file format test values
  ulong ifdpos;         // position of first IFD
  int numEntries;       // number of entries in IFD
  int i;                // counter
  ushort tag;           // IFD-tag
  ushort dataType;      // IFD data type
  ulong  numData;       // IFD: number of data values
  ulong  value;         // IFD data value
  ulong  nextIFD;       // pointer to next IFD
  ushort tileWidth=0;   // tile width in voxels
  ushort tileHeight=0;  // tile height in voxels
  ulong  tileOffset=0;  // tile offset in file
  ulong* tilePos;       // array of tile positions
  ulong  numTiles=0;    // total number of tiles in file
  int    numTilesX;     // number of tiles horizontally
  int    numTilesY;     // number of tiles vertically
  int    tpx, tpy, tpz; // tile starting position in volume data space
  int    y;             // counter for voxel lines
  int    offset;        // volume data offset to first byte of tile
  ErrorType err = OK;   // error
  uchar* raw;           // raw volume data
   
  vvDebugMsg::msg(1, "vvFileIO::loadTIFFile()");

  if (vd->getFilename()==NULL) return FILE_ERROR;
  if ( (fp = fopen(vd->getFilename(), "rb")) == NULL) 
  {
    vvDebugMsg::msg(1, "Error: Cannot open file.");
    return FILE_ERROR;
  }

  // Check file format:
  endian = vvToolshed::read16(fp);
  magic = vvToolshed::read16(fp);
  if (endian != ENDIANNESS || magic != MAGICNUMBER)
  {
    fclose(fp);
    return FORMAT_ERROR;
  }

  vd->removeSequence();
  
  // Find and process first IFD:
  ifdpos = vvToolshed::read32(fp);
  fseek(fp, ifdpos, SEEK_SET);
  numEntries = vvToolshed::read16(fp);

  vvDebugMsg::msg(2, "TIFF IFD Tags: ", numEntries);
  for (i=0; i<numEntries; ++i)    // process all IFD entries
  {
    tag      = vvToolshed::read16(fp);
    dataType = vvToolshed::read16(fp);
    numData  = vvToolshed::read32(fp);
    value    = vvToolshed::read32(fp);

    if (dataType==3) value = value >> 16;   // 16 bit values are left aligned

    if (vvDebugMsg::isActive(3))
    {
      cerr << "Tag: " << hex << setw(4) << tag << ", Data Type: " << dataType << 
        ", Data Entries: " << setw(3) << numData << ", Value: " << setw(8) 
        << value << dec << endl;
    }

    switch (tag)
    {
      case 0x100: vd->vox[0] = value; break;
      case 0x101: vd->vox[1] = value; break;
      case 0x102: if (value != 8) err = DATA_ERROR; break; // must be 8 bit per voxel
      case 0x103: if (value != 1) err = DATA_ERROR; break; // must be uncompressed
      case 0x115: if (value != 1) err = DATA_ERROR; break; // must have one component only (density)
      case 0x142: tileWidth  = (ushort)value; break;
      case 0x143: tileHeight = (ushort)value; break;
      case 0x144: numTiles = numData; tileOffset = value; break;
      case 0x80e5: vd->vox[2] = value; break;
      default: break;
    }
  }

  nextIFD = vvToolshed::read32(fp);    // check for further IFDs
  if (nextIFD==0) vvDebugMsg::msg(3, "No more IFDs in file.");
  else vvDebugMsg::msg(1, "There are more IFDs in the file.");

  if (vd->getFrameSize()==0 || err!=OK)  // check for plausibility
  {
    vvDebugMsg::msg(1, "Error: Invalid volume dimensions or file error.");
    fclose(fp);
    return DATA_ERROR;
  }  

  // Allocate memory for volume data:
  raw = new uchar[vd->getFrameSize()];

  // Load tile offsets:
  tilePos = new ulong[numTiles];
  fseek(fp, tileOffset, SEEK_SET);
  for (i=0; i<(int)numTiles; ++i)
    tilePos[i] = vvToolshed::read32(fp);    

  // Compute tiles distribution (in z direction there are as many tiles as slices):
  numTilesX = int((double)vd->vox[0] / (double)tileWidth) + 1;
  numTilesY = int((double)vd->vox[1] / (double)tileHeight) + 1;
  
  // Load volume data:
  for (i=0; i<(int)numTiles; ++i)
  {
    fseek(fp, tilePos[i], SEEK_SET);
    tpx = i % numTilesX;
    tpy = (i / numTilesX) % numTilesY;
    tpz = i / (numTilesX * numTilesY);
    offset = tpx * tileWidth + tpy * tileHeight * vd->vox[0] + tpz * vd->getSliceSize();
    for (y=0; y<tileHeight; ++y)
    {
      if (fread(raw + offset + y * vd->vox[0], 1, tileWidth, fp) != tileWidth)
      {
        vvDebugMsg::msg(1, "Error: Too short file for volume data in TIF file.");
        fclose(fp);
        delete[] raw;
        return DATA_ERROR;
      }
    }
  }

  delete[] tilePos;
  fclose(fp);
  vd->addFrame(raw, vvVolDesc::DELETE_DATA);
  vd->frames = 1;
  return OK;
}

//----------------------------------------------------------------------------
/// Loads an rgb image file as a one-sliced volume.
vvFileIO::ErrorType vvFileIO::loadRGBFile(vvVolDesc* vd)
{
  const int DIMENSIONS_OFFSET = 6;
  const int DATA_OFFSET = 512;
  FILE* fp;
  uint read;
  uchar* rawData;

  vvDebugMsg::msg(1, "vvFileIO::loadRGBFile()");
  if ( (fp=fopen(vd->getFilename(), "rb")) == NULL)
  {
    vvDebugMsg::msg(1, "Error: Cannot open RGB file.");
    return FILE_ERROR;
  }

  // Check magic number:
  fseek(fp, 0, SEEK_SET);
  if (vvToolshed::read16(fp) != 474)
  {
    vvDebugMsg::msg(1, "Error: Invalid magic number in RGB file.");
    fclose(fp);
    return FORMAT_ERROR;
  }

  // Read dimensions:
  fseek(fp, DIMENSIONS_OFFSET, SEEK_SET);
  vd->vox[0] = vvToolshed::read16(fp);
  vd->vox[1] = vvToolshed::read16(fp);
  vd->bpv    = 1;
  vd->vox[2] = 1;
  
  // Read data:
  fseek(fp, DATA_OFFSET, SEEK_SET);
  rawData = new uchar[vd->getSliceSize()]; 
  read = fread(rawData, vd->getSliceSize(), 1, fp);
  if (read != 1)
  {
    vvDebugMsg::msg(1, "Error: RGB file corrupt.");
    fclose(fp);
    delete[] rawData;
    return FILE_ERROR;
  }
  
  vd->addFrame(rawData, vvVolDesc::DELETE_DATA);  
  ++vd->frames;
  fclose(fp);
  return OK;
}

//----------------------------------------------------------------------------
/** Loads a TGA image file as a one-sliced volume.
  Only uncompressed, non-indexed image formats are supported.
*/
vvFileIO::ErrorType vvFileIO::loadTGAFile(vvVolDesc* vd)
{
  int offset_imageSpec = 8;
  int offset_idBlock = 18;
  int offset_data;

  FILE* fp;
  uint read;
  uchar* rawData;
  
  // TGA header
  uchar idLength;
  uchar colorMapType;
  uchar imageType;
  ushort colorMapOrigin;
  ushort colorMapLength;
  uchar colorMapEntrySize;


  ushort imageOriginX;
  ushort imageOriginY;
  ushort imageWidth;
  ushort imageHeigth;
  uchar imagePixelSize;
  uchar imageDescriptorByte;

  char * idBlock;

  uchar aux;
  int i;
  int total;

  vvDebugMsg::msg(1, "FileIO::loadTGAFile()");
  if ( (fp=fopen(vd->getFilename(), "rb")) == NULL)
  {
    vvDebugMsg::msg(1, "Error: Cannot open TGA file.");
    return FILE_ERROR;
  }
  
  vvDebugMsg::msg(1, "TGA file header:");
  // read ID block length
  fread(&idLength, 1, 1, fp);
  vvDebugMsg::msg(1, " ID block length: ", idLength);

  // read color map type
  fread(&colorMapType, 1, 1, fp);
  vvDebugMsg::msg(1, " Color map type: ", colorMapType);

  // read image type
  fread(&imageType, 1, 1, fp);
  vvDebugMsg::msg(1, " Image type: ", imageType);

  // read color map infos
  if ( 0 != colorMapType) 
  {
	  colorMapOrigin = vvToolshed::read16(fp, vvToolshed::VV_LITTLE_END);
	  colorMapLength = vvToolshed::read16(fp, vvToolshed::VV_LITTLE_END);
	  fread(&colorMapEntrySize, 1, 1, fp);
	  vvDebugMsg::msg(1, " Color map origin: ", colorMapOrigin);
	  vvDebugMsg::msg(1, " Color map length: ", colorMapLength);
	  vvDebugMsg::msg(1, " Color map entry size: ", colorMapEntrySize);
  }

  // read image specification block
  fseek(fp, offset_imageSpec, SEEK_SET);
  imageOriginX = vvToolshed::read16(fp, vvToolshed::VV_LITTLE_END);
  imageOriginY = vvToolshed::read16(fp, vvToolshed::VV_LITTLE_END);
  imageWidth   = vvToolshed::read16(fp, vvToolshed::VV_LITTLE_END);
  imageHeigth  = vvToolshed::read16(fp, vvToolshed::VV_LITTLE_END);
  fread(&imagePixelSize, 1, 1, fp);
  fread(&imageDescriptorByte, 1, 1, fp);

  vvDebugMsg::msg(1, " Origin X: ", imageOriginX);
  vvDebugMsg::msg(1, " Origin Y: ", imageOriginY);
  vvDebugMsg::msg(1, " Width: ", imageWidth);
  vvDebugMsg::msg(1, " Height: ", imageHeigth);
  vvDebugMsg::msg(1, " Pixel size: ", imagePixelSize);
  vvDebugMsg::msg(1, " Descriptor byte: ", imageDescriptorByte);
  vvDebugMsg::msg(1, "  Number of attribute bits per pixel: ", (imageDescriptorByte & 7));
  vvDebugMsg::msg(1, "  Reserved: ", (imageDescriptorByte & 8)>>3);
  vvDebugMsg::msg(1, "  Screen origin bit: ", (imageDescriptorByte & 16)>>4);
  vvDebugMsg::msg(1, "  Data storage interleaving flag: ", (imageDescriptorByte & 96)>>5);
  
  // read ID block
  fseek(fp, offset_idBlock, SEEK_SET);
  if (0 < idLength) 
  {
  	idBlock = new char[idLength+1];
	  fread(idBlock, 1, idLength, fp);
	  if (NULL != idBlock) 
    {
  		idBlock[idLength]='\0';
		  vvDebugMsg::msg(1, " Image ID block: ", idBlock);
		  delete(idBlock);
	  }
  }
  
  if (2!=imageType) 
  {
    vvDebugMsg::msg(1, "Error: Image type not supported,");
    vvDebugMsg::msg(1, "please use uncompressed RGB(A) only!");
    fclose(fp);
  	return FILE_ERROR;
  }

  // assign image params to volume
  vd->vox[0] = imageWidth;
  vd->vox[1] = imageHeigth;
  vd->vox[2] = 1;
  vd->bpv    = imagePixelSize/8;
  vd->frames = 1;

  total = vd->getSliceSize();

  // compute data block offset
  offset_data = offset_idBlock + idLength;

	// ### TODO: include color map offset (colormaps not supported anyways)

  // Read data:
  fseek(fp, offset_data, SEEK_SET);
  rawData = new uchar[total]; 

  read = fread(rawData, total, 1, fp);
  if (read != 1)
	{
    vvDebugMsg::msg(1, "Error: TGA file corrupt.");
	  fclose(fp);
	  delete[] rawData;
	  return FILE_ERROR;
	}

	// Byte per pixel value of 3 or 4 implies that the image is RGB(A). 
  // However TGA stores it as BGR(A) so we'll have to swap R and B.
  if (vd->bpv >= 3)
	  for (i=0; i < total; i += vd->bpv) 
    {
		  aux = rawData[i];
		  rawData[i] = rawData[i+2];
		  rawData[i+2] = aux;
	  }
	
  vd->addFrame(rawData, vvVolDesc::DELETE_DATA);  
  vd->flip(vvVolDesc::Y_AXIS);

  fclose(fp);
  return OK;
}

//----------------------------------------------------------------------------
/** Loads a raw volume file w/o knowing its structure
 Several automatic detection algorithms are tried.
 If filename is of format: 
 <PRE>
    <filename|width|x|height|x|slices|.dat> 
    (example: "cthead256x256x128.dat")
 </PRE>
 a volume size of |width| x |height| x |slices| voxels is tried
 in addition to the automatic detection algorithms.
*/
vvFileIO::ErrorType vvFileIO::loadRawFile(vvVolDesc* vd)
{
  const int NUM_ALGORITHMS = 4;     // number of different size detection algorithms
  char* ptr;                        // pointer to current character
  long lSize;
  char filename[1024];              // buffer for filename
  int size, voxels;
  int width, height, slices, bpv;   // volume characteristics
  int remainder;
  int cubRoot;                      // cubic root
  int sqrRoot;                      // square root
  int attempt;
  int factor;
  int i;

  vvDebugMsg::msg(1, "vvFileIO::loadRawFile(0)");

  lSize = vvToolshed::getFileSize(vd->getFilename());
  if (lSize <= 0) return FILE_ERROR;
  if (lSize > INT_MAX) return FORMAT_ERROR;
  size = (int)lSize;

  for (attempt=0; attempt<NUM_ALGORITHMS; ++attempt) // try different ways to find the volume dimensions
  {
    for (bpv=1; bpv<=4; ++bpv)
    {
      if ((size % bpv) != 0) continue;
      else voxels = size / bpv;
      cubRoot = (int)powf((float)voxels, 1.0f/3.0f);
      width = height = slices = 0; 
      switch (attempt)
      {
        case 0:         // check for dimensions given in filename
          vvToolshed::extractFilename(filename, vd->getFilename());

          // Search for beginning of size information string:
          ptr = strchr(filename, '.');
          *ptr = '\0';
          i = 0;
          while (ptr > filename && i < 3) 
          {
            --ptr;
            if (!isdigit(*ptr)) 
            {
              switch (i)
              {
                case 0: slices = atoi(ptr+1); break;
                case 1: height = atoi(ptr+1); break;
                case 2: width  = atoi(ptr+1); break;
                default: break;              
              }            
              ++i;
              *ptr = '\0';    // convert delimiters to string terminators
            }
            else if (ptr==filename && i==2) width = atoi(ptr);
          }
          break;

        case 1: 
        default:        // Check for cubic volume:
          width = height = slices = cubRoot;
          break;

        case 2:         // Check for slices being powers of 2:
          width = vvToolshed::getTextureSize(cubRoot);
          sqrRoot = (int)sqrt((double)voxels / (double)width);
          height = vvToolshed::getTextureSize(sqrRoot);
          slices = voxels / width / height;
          break;

        case 3:         // Check for square slices and slice edge length greater than volume depth:
          width = slices = 1;
          remainder = size;
          while ((factor = vvToolshed::getLargestPrimeFactor(remainder)) > 1 && width < cubRoot)
          {
            if ((remainder % (factor*factor)) == 0)   // is factor contained twice?
            {
              width *= factor;
              remainder /= (factor * factor);
            }
            else 
            {
              slices *= factor;
              remainder /= factor;
            }
          }
          slices *= remainder;
          height = width;
          break;
      }
      if (bpv * width * height * slices == size)
        return loadRawFile(vd, width, height, slices, bpv, 0);
    }    
  }
  return FORMAT_ERROR;
}

//----------------------------------------------------------------------------
/** Saves a raw data file: no header information is written.
 The number of bytes per voxel is similar to memory format.
 Only one frame of a time dependent dataset is written.
*/
vvFileIO::ErrorType vvFileIO::saveRawFile(vvVolDesc* vd)
{
  FILE* fp;             // volume file pointer
  int frameSize;        // size of a frame in bytes
  uchar* raw;           // raw volume data
  
  vvDebugMsg::msg(1, "vvFileIO::saveRawFile()");

  // Save volume data:
  frameSize = vd->getFrameSize();
  raw = vd->getRaw();  // save only current frame of loaded sequence
  if (frameSize==0 || raw==NULL) 
  {
    return VD_ERROR;
  }

  if ( (fp = fopen(vd->getFilename(), "wb")) == NULL) // now open file to write
  {
    vvDebugMsg::msg(1, "Error: Cannot open file to write.");
    return FILE_ERROR;
  }

  // Write volume data:
  if ((int)fwrite(raw, 1, frameSize, fp) != frameSize)
  {
    vvDebugMsg::msg(1, "Error: Cannot write voxel data to file.");
    fclose(fp);
    return FILE_ERROR;
  }
  
  fclose(fp);
  return OK;
}

//----------------------------------------------------------------------------
/** Loads a raw volume file of which the structure is known
  @param w      width
  @param h      height
  @param s      slices (use 1 for 2D image files)
  @param b      byte per voxel
  @param header header size in bytes (= number of bytes to skip at beginning of file)
*/
vvFileIO::ErrorType vvFileIO::loadRawFile(vvVolDesc* vd, int w, int h, int s, int b, int header)
{
  FILE* fp;
  uint read;
  uchar* rawData;

  if (b<1 || b>4) return FORMAT_ERROR;
  
  vvDebugMsg::msg(1, "vvFileIO::loadRawFile(1)");
  if ( (fp=fopen(vd->getFilename(), "rb")) == NULL)
  {
    vvDebugMsg::msg(1, "Error: Cannot open raw file.");
    return FILE_ERROR;
  }

  vd->vox[0] = w;
  vd->vox[1] = h;
  vd->vox[2] = s;
  vd->bpv    = b;

  fseek(fp, header, SEEK_SET);    // skip header
  rawData = new uchar[vd->getFrameSize()]; 
  read = fread(rawData, vd->getFrameSize(), 1, fp);
  if (read != 1)
  {
    vvDebugMsg::msg(1, "Error: Raw file corrupt.");
    fclose(fp);
    delete[] rawData;
    return FILE_ERROR;
  }

  fclose(fp);
  vd->addFrame(rawData, vvVolDesc::DELETE_DATA);
  ++vd->frames;
  return OK;
}

//----------------------------------------------------------------------------
/// Loads a PGM or PPM binary image file.
vvFileIO::ErrorType vvFileIO::loadPXMRawImage(vvVolDesc* vd)
{
  const int BUFSIZE = 128;
  FILE* fp;
  uint read;
  uchar* rawData;
  char buf[3][BUFSIZE];
  bool isPGM;         // true=PGM, false=PPM

  vvDebugMsg::msg(1, "vvFileIO::loadPXMRawImage()");
  if ( (fp=fopen(vd->getFilename(), "rb")) == NULL)
  {
    vvDebugMsg::msg(1, "Error: Cannot open PGM/PPM file.");
    return FILE_ERROR;
  }

  // Read magic number:
  fgets(buf[0], BUFSIZE, fp);
  if (vvToolshed::strCompare("P5", buf[0], 2) == 0)
    isPGM = true;
  else if (vvToolshed::strCompare("P6", buf[0], 2) == 0)
    isPGM = false;
  else
  {
    fclose(fp);
    vvDebugMsg::msg(1, "Error: Wrong magic number in PGM/PPM file. Use binary format.");
    return DATA_ERROR;
  }

  // Read width and height:
  do
  {
    fgets(buf[0], BUFSIZE, fp);
  } while (buf[0][0] == '#');
  sscanf(buf[0], "%s %s", buf[1], buf[2]);
  vd->vox[0] = atoi(buf[1]);
  vd->vox[1] = atoi(buf[2]);

  // Read maxval:
  do
  {
    fgets(buf[0], BUFSIZE, fp);
  } while (buf[0][0] == '#');

  // Read image data:
  vd->vox[2] = 1;
  if (isPGM) vd->bpv = 1;
  else       vd->bpv = 3;
  rawData = new uchar[vd->getFrameSize()]; 
  read = fread(rawData, vd->getFrameSize(), 1, fp);
  if (read != 1)
  {
    vvDebugMsg::msg(1, "Error: PGM/PPM file corrupt.");
    fclose(fp);
    delete[] rawData;
    return DATA_ERROR;
  }

  fclose(fp);
  vd->addFrame(rawData, vvVolDesc::DELETE_DATA);
  ++vd->frames;
  return OK;
}

//----------------------------------------------------------------------------
/** Loads a DICOM 3.0 image file
 (DICOM = Digital Imaging COmmunications in Medicine)
 @param vd        volume description
 @param dcmSeq    DICOM sequence ID (NULL if not required)
 @param dcmSlice  DICOM slice ID (NULL if not required)
 @param dcmSPos   DICOM slice location (NULL if not required)
*/
vvFileIO::ErrorType vvFileIO::loadDicomFile(vvVolDesc* vd, int* dcmSeq, int* dcmSlice, float* dcmSPos)
{
  vvDicom* dicomReader;
  vvDicomProperties prop;
  int i;

  dicomReader = new vvDicom(&prop);
  if (!dicomReader->readDicomFile((char*)vd->getFilename()))
  {
    delete dicomReader;
    vvDebugMsg::msg(1, "Error: Cannot open Dicom file.");
    return FILE_ERROR;
  }
  
  if (vvDebugMsg::isActive(1)) prop.print();

  // Make sure variables are tested for NULL because they might be default:
  if (dcmSeq   != NULL) *dcmSeq   = prop.sequence;
  if (dcmSlice != NULL) *dcmSlice = prop.image;
  if (dcmSPos  != NULL) *dcmSPos  = prop.slicePos;

  vd->vox[0] = prop.width;
  vd->vox[1] = prop.height;
  vd->vox[2] = 1;
  for (i=0; i<3; ++i)
  {
    vd->dist[i] = prop.dist[i];
  }
  vd->bpv = prop.bpp;
  vd->addFrame(prop.raw, vvVolDesc::DELETE_DATA);
  ++vd->frames;
  
  // Make big endian data:
  if (prop.littleEndian) vd->toggleEndianness();
  
  // Shift bits so that most significant used bit is leftmost:
  vd->bitShiftData(prop.highBit - (prop.bpp * 8 - 1));
  
  // Make unsigned data:
  if (prop.isSigned) vd->toggleSign();
  
  delete dicomReader;
  return OK;
}

//----------------------------------------------------------------------------
/** Loads a Visible Human anatomic (photo) slice file.
*/
vvFileIO::ErrorType vvFileIO::loadVHDAnatomicFile(vvVolDesc* vd)
{
  ErrorType err;
  
  err = loadRawFile(vd, 2048, 1216, 1, 3, 0);
  if (err != OK) return err;
  vd->convertRGBPlanarToRGBInterleaved();
  vd->crop(300, 100, 0, 1400, 950, 1);    // images are 2048 x 1216 but contain unnecessary information in the border region
  return OK;
}

//----------------------------------------------------------------------------
/** Loads a Visible Human MRI slice file.
*/
vvFileIO::ErrorType vvFileIO::loadVHDMRIFile(vvVolDesc* vd)
{
  ErrorType err;

  err = loadRawFile(vd, 256, 256, 1, 2, 7900);
  if (err != OK) return err;
  vd->bitShiftData(-4);        // image is (about) 12 bit, shift it to be in correct 16 bit representation
  return OK;
}

//----------------------------------------------------------------------------
/** Loads a Visible Human CT slice file.
*/
vvFileIO::ErrorType vvFileIO::loadVHDCTFile(vvVolDesc* vd)
{
  ErrorType err;

  err = loadRawFile(vd, 512, 512, 1, 2, 3416);
  if (err != OK) return err;
  vd->bitShiftData(-4);        // image is 12 bit, shift to appear as 16 bit
  return OK;
}

//----------------------------------------------------------------------------
/** Loader for BrainVoyager VMR files.
  VMR files contain anatomical 3D data stored as bytes. 
  The format is very simple consisting of a small header prior to the 
  actual data. The only important point is that you understand how the 
  three axes are ordered (see below).<P>
  VMR header:<BR>
  <PRE>
  BYTES	 DATA TYPE	            DESCRIPTION
  2      16 bit little endian  DimX, dimension of X axis
  2	     16 bit little endian  DimY, dimension of Y axis
  2	     16 bit little endian  DimZ, dimension of Z axis
  </PRE>
  Each data element (intensity value) is represented in 1 byte. 
  The data is organized in three loops: DimZ, DimY, DimX
*/
vvFileIO::ErrorType vvFileIO::loadVMRFile(vvVolDesc* vd)
{
  FILE* fp;             // volume file pointer
  int frameSize;        // size of a frame in bytes
  uchar* raw;           // raw volume data
   
  vvDebugMsg::msg(1, "vvFileIO::loadVMRFile()");

  if (vd->getFilename()==NULL) return FILE_ERROR;
  if ( (fp = fopen(vd->getFilename(), "rb")) == NULL) 
  {
    vvDebugMsg::msg(1, "Error: Cannot open file.");
    return FILE_ERROR;
  }

  vd->removeSequence();

  // Read header:
  vd->vox[0] = vvToolshed::read16(fp, vvToolshed::VV_LITTLE_END);
  vd->vox[1] = vvToolshed::read16(fp, vvToolshed::VV_LITTLE_END);
  vd->vox[2] = vvToolshed::read16(fp, vvToolshed::VV_LITTLE_END);
  vd->frames = 1;
  vd->bpv    = 1;

  // Create new data space for volume data:
  if ((sections & RAW_DATA) != 0)
  {
    frameSize = vd->getFrameSize();
    raw = new uchar[frameSize];
  
    // Load volume data:
    if ((int)fread(raw, 1, frameSize, fp) != frameSize)
    {
      vvDebugMsg::msg(1, "Error: Insufficient voxel data in VMR file.");
      fclose(fp);
      delete[] raw;
      return FILE_ERROR;
    }
  
    vd->addFrame(raw, vvVolDesc::DELETE_DATA);
  }
  fclose(fp);
  return OK;
}

//----------------------------------------------------------------------------
/** Loader for BrainVoyager VTC files.
  A VTC file contains the functional data (time series) of one experimental 
  run in a 3D format, i.e. in Talairach space. The binary file contains a 
  variable-length header followed by the actual 4D data.<P>
  Header:
  <PRE>
  BYTES DATA TYPE   DESCRIPTION
  2     short int   version number
  N     byte	 	    name of FMR file whose STC data has been transformed
  M     byte        name of the linked protocol (PRT) file
  2     short int   NrOfVolumes (number of volumes, measurements, time points)
  2     short int   VTC-resolution, i.e. 3 -> one voxel = 3 x 3 x 3 mm
  2     short int   XStart
  2     short int   XEnd
  2     short int   YStart
  2     short int   YEnd
  2     short int   ZStart
  2     short int   ZEnd
  2     short int	 	Hemodynamic delay, simple shift value
  4     float	 	    TR [ms]
  4     float       Hemodynamic function, delta parameter
  4     float       Hemodynamic function, tau parameter
  2     short int   Segment size, used for time course separation
  2     short int   Segment offset, used for time course separation
  </PRE>
  The order of voxels in the file is:<BR>
  timesteps, voxels/line, lines, slices
*/
vvFileIO::ErrorType vvFileIO::loadVTCFile(vvVolDesc* vd)
{
  FILE* fp;             // volume file pointer
  uchar* raw;           // raw volume data
  uchar* buf;           // buffer for data from file
  uchar* bufPtr;        // pointer to data in file buffer
  uchar bak;            // backup value
  uchar** frameRaw;     // pointer to beginning of raw data in each animation frame
  int frameSize;        // size of a frame in bytes
  int vtcSliceSize;     // slice size of VTC file (contains all time step data)
  int start, end;
  int version;          // file version
  int f,i,x,y;
   
  vvDebugMsg::msg(1, "vvFileIO::loadVTCFile()");

  if (vd->getFilename()==NULL) return FILE_ERROR;
  if ( (fp = fopen(vd->getFilename(), "rb")) == NULL) 
  {
    vvDebugMsg::msg(1, "Error: Cannot open file.");
    return FILE_ERROR;
  }

  vd->removeSequence();

  // Read header:
  vd->bpv = 2;
  version = vvToolshed::read16(fp, vvToolshed::VV_LITTLE_END);   // read version number
  while (fgetc(fp)!=0);       // ignore FMR file name string
  while (fgetc(fp)!=0);       // ignore PRT file name string
  vd->frames = vvToolshed::read16(fp, vvToolshed::VV_LITTLE_END);
  vd->dist[0] = vd->dist[1] = vd->dist[2] = vvToolshed::read16(fp, vvToolshed::VV_LITTLE_END);
  for (i=0; i<3; ++i)
  {
    start      = vvToolshed::read16(fp, vvToolshed::VV_LITTLE_END);
    end        = vvToolshed::read16(fp, vvToolshed::VV_LITTLE_END);
    vd->vox[i] = int((end - start) / vd->dist[0]);
  }
  if (version==2)   // the following parameters are only in the header if version equals 2
  {
    // Ignore the extra header information:
    vvToolshed::read16(fp, vvToolshed::VV_LITTLE_END);
    vvToolshed::readFloat(fp, vvToolshed::VV_LITTLE_END);
    vvToolshed::readFloat(fp, vvToolshed::VV_LITTLE_END);
    vvToolshed::readFloat(fp, vvToolshed::VV_LITTLE_END);
    vvToolshed::read16(fp, vvToolshed::VV_LITTLE_END);
  }

  if ((sections & RAW_DATA) != 0)
  {
    frameSize = vd->getFrameSize();

    // First allocate space for entire volume animation:
    for (f=0; f<vd->frames; ++f)
    {
      raw = new uchar[frameSize];   // create new data space for volume data
      assert(raw);
      vd->addFrame(raw, vvVolDesc::DELETE_DATA);  // add uninitialized data to volume
    }
  
    // Now we can fill the frames with the data from disk:
    vtcSliceSize = vd->getMovieSize() / vd->vox[2];
    buf = new uchar[vtcSliceSize];
    frameRaw = new uchar*[vd->frames];
  
    for (f=0; f<vd->frames; ++f)    // store pointers to frame data in array for speed
      frameRaw[f] = vd->getRaw(f);

    for (i=0; i<vd->vox[2]; ++i)
    {
      if ((int)fread(buf, 1, vtcSliceSize, fp) != vtcSliceSize)
      {
        vvDebugMsg::msg(1, "Error: Insufficient voxel data in file.");
        fclose(fp);
        vd->removeSequence();
        return FILE_ERROR;
      }
      bufPtr = buf;

      // Copy data from buffer to actual volume storage:
      for (y=0; y<vd->vox[1]; ++y)
        for (x=0; x<vd->vox[0]; ++x)
          for (f=0; f<vd->frames; ++f)
          {
            // Swap the bytes because they are stored as little endian:            
            bak = *bufPtr;
            *bufPtr = *(bufPtr+1);
            *(bufPtr+1) = bak;

            memcpy(frameRaw[f], bufPtr, 2);
            frameRaw[f] += 2;
            bufPtr += 2;
          }
    }
    delete[] frameRaw;
    delete[] buf;
  }
  fclose(fp);
  return OK;
}

//----------------------------------------------------------------------------
/** Loader for voxel file in nrrd (teem volume file) format.
*/
vvFileIO::ErrorType vvFileIO::loadNrrdFile(vvVolDesc* vd)
{
  FILE* fp;               // volume file pointer
  vvTokenizer* tokenizer; // stream tokenizer
  vvTokenizer::TokenType ttype = vvTokenizer::VV_NOTHING;  // token type
  vvTokenizer::TokenType prevTT = vvTokenizer::VV_NOTHING; // previous token type
  uchar* raw;               // raw volume data
  int f, i;                 // counters
  int frameSize;            // size of a frame in bytes
  int dimension = 0;        // dimension of the volume dataset
  bool bigEnd = true;       // true = big endian
   
  vvDebugMsg::msg(1, "vvFileIO::loadNrrdFile()");

  if (vd->getFilename()==NULL) return FILE_ERROR;
	
  if ( (fp = fopen(vd->getFilename(), "rb")) == NULL) 
  {
    vvDebugMsg::msg(1, "Error: Cannot open file.");
    return FILE_ERROR;
  }
  vd->removeSequence();   // delete previous volume sequence
  
  for (i=0; i<(int)strlen(nrrdID); ++i)
  {
    if (fgetc(fp) != nrrdID[i])
    {
      cerr << "Error: Invalid file ID string." << endl;
      fclose(fp);
      return DATA_ERROR;
    }
  } 

  // Create tokenizer:
  tokenizer = new vvTokenizer(fp);
  tokenizer->setCommentCharacter('#');
  tokenizer->setEOLisSignificant(true);
  tokenizer->setCaseConversion(vvTokenizer::VV_LOWER);
  tokenizer->setParseNumbers(true);
  tokenizer->setWhitespaceCharacter(':');
  
  // Parse header:
  vd->vox[2] = 1;
  vd->frames = 1;   // default values
  do
  {
    prevTT = ttype;
    ttype = tokenizer->nextToken();
    if (ttype == vvTokenizer::VV_EOF || ttype == vvTokenizer::VV_NUMBER) 
    {
      cerr << "Invalid nrrd file format." << endl;
      delete tokenizer;
      fclose(fp);
      return FORMAT_ERROR;
    }
    else if (ttype == vvTokenizer::VV_EOL)
    { // do nothing
    }
    else if (strcmp(tokenizer->sval, "content")==0)
    {
      // ignore content information
      tokenizer->nextLine();
    }
    else if (strcmp(tokenizer->sval, "type")==0)
    {
      ttype = tokenizer->nextToken();
      if (ttype != vvTokenizer::VV_WORD || strcmp(tokenizer->sval, "unsigned")!=0)
        cerr << "unknown type" << endl;
      else
      {
        ttype = tokenizer->nextToken();
        if (strcmp(tokenizer->sval, "char")==0) vd->bpv = 1;
        else if (strcmp(tokenizer->sval, "short")==0) vd->bpv = 2;
        else cerr << "unknown type" << endl;
      }
    }
    else if (strcmp(tokenizer->sval, "dimension")==0)
    {
      ttype = tokenizer->nextToken();
      if (ttype == vvTokenizer::VV_NUMBER)
      {
        dimension = int(tokenizer->nval);
        if (dimension < 1 || dimension > 4) cerr << "dimension must be 1 to 4" << endl;
      }
      else cerr << "invalid dimension" << endl;      
    }
    else if (strcmp(tokenizer->sval, "sizes")==0)
    {
      for (i=0; i<dimension; ++i)
      {
        ttype = tokenizer->nextToken();
        if (i==0)
        {
          // Guess if first entry is number of modalities or width.
          // Assume number of modalities if first size is 2, 3, or 4.
          switch (int(tokenizer->nval))
          {
            case 2:  
            case 3:  
            case 4: vd->bpv = int(tokenizer->nval); 
                    vd->stype = vvVolDesc::VV_MULTIMODAL;
                    break;
            default: vd->vox[0] = int(tokenizer->nval); break;
          }
        }
        else
        {
          if (vd->bpv>=2 && vd->bpv<=4) vd->vox[i-1] = int(tokenizer->nval);
          else if (i==3) vd->frames = int(tokenizer->nval);
          else vd->vox[i] = int(tokenizer->nval);
        }
      }
    }
    else if (strcmp(tokenizer->sval, "spacings")==0)
    {
      bool multiModal = false;
      for (i=0; i<dimension; ++i)
      {
        ttype = tokenizer->nextToken();
        if (i==0 && ttype==vvTokenizer::VV_WORD)  // if first value is NaN, expect multi-modal data
        {
          vd->dt = 0.0f;
          multiModal = true;
        }
        else if (i>0 && multiModal)  // still multi-modal data
        {
          vd->dist[i-1] = tokenizer->nval;
        }
        else    // only one modality
        {
          if (i==3) vd->dt = tokenizer->nval;
          else vd->dist[i] = tokenizer->nval;
        }
      }
    }
    else if (strcmp(tokenizer->sval, "endian")==0)
    {
      ttype = tokenizer->nextToken();
      if (strcmp(tokenizer->sval, "little") == 0) bigEnd = false;
      else bigEnd = true;
    }
    else if (strcmp(tokenizer->sval, "encoding")==0)
    {
      ttype = tokenizer->nextToken();
      if (strcmp(tokenizer->sval, "raw") != 0) 
      {
        cerr << "Can only process raw data." << endl;
        delete tokenizer;
        fclose(fp);
        return FORMAT_ERROR;
      }
    }
    else
    {
      tokenizer->nextLine();
    }
  } while (ttype != vvTokenizer::VV_EOL || prevTT != vvTokenizer::VV_EOL);  // stop when two EOL in a row
  delete tokenizer;

  frameSize = vd->getFrameSize();

  // Load volume data:
  if ((sections & RAW_DATA) != 0)
  {
    for (f=0; f<vd->frames; ++f)
    {
      raw = new uchar[frameSize];   // create new data space for volume data
      if ((int)fread(raw, 1, frameSize, fp) != frameSize)
      {
        vvDebugMsg::msg(1, "Error: Insuffient voxel data in file.");
        fclose(fp);
        delete[] raw;
        return DATA_ERROR;
      }
      vd->addFrame(raw, vvVolDesc::DELETE_DATA);
    }
  }
  
  if (!bigEnd) vd->toggleEndianness();

  // Clean up:
  fclose(fp);
  return OK;
}

//----------------------------------------------------------------------------
/** Load XIMG image file. This file format was created by General Electric.
*/
vvFileIO::ErrorType vvFileIO::loadXIMGFile(vvVolDesc* vd)
{
  FILE* fp;
  uchar* rawData;
  uint read;
  char magic[4];    // XIMG magic number
  int offset;       // offset to data area
  int compression;  // compression format
  int i;

  vvDebugMsg::msg(1, "vvFileIO::loadXIMGFile()");
  if ( (fp=fopen(vd->getFilename(), "rb")) == NULL)
  {
    vvDebugMsg::msg(1, "Error: Cannot open XIMG file.");
    return FILE_ERROR;
  }

  // Read magic number:
  for (i=0; i<4; ++i)
  {
    magic[i] = vvToolshed::read8(fp);
  }
  if (magic[0]!='I' || magic[1]!='M' || magic[2]!='G' || magic[3]!='F')
  {
    fclose(fp);
    cerr << "Wrong magic number in XIMG file." << endl;
    return DATA_ERROR;
  }
  
  // Read offset to data area:
  offset = vvToolshed::read32(fp);
  
  // Read image size:
  vd->vox[0] = vvToolshed::read32(fp);
  vd->vox[1] = vvToolshed::read32(fp);
  
  // Read bpv:
  vd->bpv = vvToolshed::read32(fp) / 8;
  
  // Read compression:
  compression = vvToolshed::read32(fp);
  if (compression!=1) 
  {
    fclose(fp);
    cerr << "Compression type must be 'rectangular'." << endl;
    return DATA_ERROR;
  }

  // Read image data:
  fseek(fp, offset, SEEK_SET);
  vd->vox[2] = 1;
  rawData = new uchar[vd->getFrameSize()]; 
  read = fread(rawData, vd->getFrameSize(), 1, fp);
  if (read != 1)
  {
    vvDebugMsg::msg(1, "Error: XIMG file corrupt.");
    fclose(fp);
    delete[] rawData;
    return DATA_ERROR;
  }

  fclose(fp);
  vd->addFrame(rawData, vvVolDesc::DELETE_DATA);
  ++vd->frames;
  return OK;
}

//----------------------------------------------------------------------------
/** Saves all slices of the current volume as PPM or PGM images.
 A numbered suffix of four digits will be added to the file names.
 @param overwrite   true to overwrite existing files
*/
vvFileIO::ErrorType vvFileIO::savePXMSlices(vvVolDesc* vd, bool overwrite)
{
  FILE* fp;
  ErrorType err = OK;
  int digits;             // number of digits used for file numbering
  char** filenames;       // list of filenames
  int len;                // filename length
  int i, j, k;
  int sliceSize;
  int tmpSliceSize = 0;
  char buffer[1024];
  uchar* slice;           // original slice data
  uchar* tmpSlice = NULL; // temporary slice data

  vvDebugMsg::msg(1, "vvFileIO::savePXMSlices()");

  if (vd->frames<1 || vd->vox[2]<1) return DATA_ERROR;
  
  // Generate file names:
  digits = 1 + int(log((double)vd->vox[2]) / log(10.0));
  filenames = new char*[vd->vox[2]];
  len = strlen(vd->getFilename());
  for (i=0; i<vd->vox[2]; ++i)
  {
    filenames[i] = new char[len + digits + 2];   // add 2 for '-' and '\0'
    vvToolshed::extractDirname(buffer, vd->getFilename());
    strcpy(filenames[i], buffer);
    vvToolshed::extractBasename(buffer, vd->getFilename());
    strcat(filenames[i], buffer);
    if (vd->vox[2] > 1)
    {
      sprintf(buffer, "-%0*d.", digits, i);
      strcat(filenames[i], buffer);
    }
    else
      strcat(filenames[i], ".");
    if (vd->bpv<=2) strcat(filenames[i], "pgm");
    else strcat(filenames[i], "ppm");
  }

  // Check files for existence:
  if (!overwrite)
    for (i=0; i<vd->vox[2]; ++i)
      if (vvToolshed::isFile(filenames[i]))    // check if file exists
      {
        vvDebugMsg::msg(1, "Error - file exists: ", filenames[i]);
        err = FILE_EXISTS;
      }

  // Write files:
  sliceSize = vd->getSliceSize();
  if (vd->bpv==2 || vd->bpv==4)
  {
    tmpSliceSize = vd->vox[0] * vd->vox[1] * (vd->bpv-1);
    tmpSlice = new uchar[tmpSliceSize];
  }
  for (i=0; i<vd->vox[2] && err==OK; ++i)
  {
    // Open file to write:
    if ( (fp = fopen(filenames[i], "wb")) == NULL)
    {
      err = FILE_ERROR;
      continue;
    }

    // Write header:
    if (vd->bpv <= 2) fprintf(fp, "P5\n");   // grayscale
    else fprintf(fp, "P6\n");                // RGB

    fprintf(fp, "%d %d\n", vd->vox[0], vd->vox[1]);  // write dimensions
    fprintf(fp, "%d\n", 255);   // write maximum value

    // Write data:
    slice = vd->getRaw() + i * sliceSize;
    switch (vd->bpv)
    {
      case 1:
      case 3:
        if ((int)fwrite(slice, sliceSize, 1, fp) != 1)
          err = FILE_ERROR;
        break;
      case 2:
      case 4:
        for (j=0; j<vd->getSliceVoxels(); ++j)
          for (k=0; k<vd->bpv-1; ++k)
            tmpSlice[j * (vd->bpv-1) + k] = slice[j * vd->bpv + k];
        if ((int)fwrite(tmpSlice, tmpSliceSize, 1, fp) != 1)
          err = FILE_ERROR;
        break;
      default: break;
    }

    fclose(fp);
  }

  // Free memory:    
  if (vd->bpv==2 || vd->bpv==4)
    delete[] tmpSlice;
  for (i=0; i<vd->vox[2]; ++i)
    delete[] filenames[i];
  delete[] filenames;

  return err;
}

//----------------------------------------------------------------------------
/** Save volume data to a volume file. 
  The file format is determined from the filename extension.
  @param vd        volume description 
  @param overwrite true to overwrite existing file
  @param sec       bit encoded list of file sections to be saved (if present in file).
                   This value defaults to saving all data to a file.
  @return NO_ERROR if successful
*/
vvFileIO::ErrorType vvFileIO::saveVolumeData(vvVolDesc* vd, bool overwrite, int sec)
{
  vvDebugMsg::msg(1, "vvFileIO::saveVolumeData(), file name: ", vd->getFilename());

  if (vd==NULL) return PARAM_ERROR;                    // volume description missing

  if (vd->getFilename()==NULL) return PARAM_ERROR;
  if (strlen(vd->getFilename()) < 3) return PARAM_ERROR;    // filename too short

  if (!overwrite && vvToolshed::isFile(vd->getFilename()))    // check if file exists
  {
    vvDebugMsg::msg(1, "Error: File exists:", vd->getFilename());
    return FILE_EXISTS;
  }

  sections = sec;

  if (vvToolshed::isSuffix(vd->getFilename(), ".rvf"))
    return saveRVFFile(vd);

  if (vvToolshed::isSuffix(vd->getFilename(), ".xvf"))
    return saveXVFFile(vd);

  if (vvToolshed::isSuffix(vd->getFilename(), ".dat"))
    return saveRawFile(vd);

  if (vvToolshed::isSuffix(vd->getFilename(), ".nrd"))
    return saveNrrdFile(vd);

  if (vvToolshed::isSuffix(vd->getFilename(), ".ppm") ||
      vvToolshed::isSuffix(vd->getFilename(), ".pgm"))
    return savePXMSlices(vd, overwrite);

  vvDebugMsg::msg(1, "Error in saveVolumeData: unknown extension");
  return PARAM_ERROR;
}

//----------------------------------------------------------------------------
/** Load volume data from a volume file.
  If filename is undefined, compute default volume.
  @param vd   volume description 
  @param sec  bit encoded list of file sections to be loaded (if present in file).
              This value defaults to loading all data in file.
*/
vvFileIO::ErrorType vvFileIO::loadVolumeData(vvVolDesc* vd, int sec)
{
  vvDebugMsg::msg(1, "vvFileIO::loadVolumeData()");

  ErrorType err = OK;
  char* suffix;

  if (vd==NULL) return PARAM_ERROR;                    // volume description missing

  if (vd->getFilename()==NULL || strlen(vd->getFilename()) == 0)
  {
    vd->setFilename("default.xvf");
    return computeDefaultVolume(vd, 1); 
  }

  if (vvToolshed::isFile(vd->getFilename())==false)
    return FILE_NOT_FOUND;

  sections = sec;

  suffix = new char[strlen(vd->getFilename())+1];
  vvToolshed::extractExtension(suffix, vd->getFilename());

  // Load files according to extension:
  if (vvToolshed::strCompare(suffix, "wl") == 0)
    err = loadWLFile(vd);

  else if (vvToolshed::strCompare(suffix, "rvf") == 0)
    err = loadRVFFile(vd);

  else if (vvToolshed::strCompare(suffix, "xvf") == 0)
    err = loadXVFFile(vd);

  else if (vvToolshed::strCompare(suffix, "avf") == 0)
    err = loadAVFFile(vd);

  else if (vvToolshed::strCompare(suffix, "xb7") == 0)
    err = loadXB7File(vd);

  else if (vvToolshed::strCompare(suffix, "vf")  == 0)
    err = loadVFFile(vd);

  else if (vvToolshed::strCompare(suffix, "asc") == 0)
    err = loadASCFile(vd);

  else if (vvToolshed::strCompare(suffix, "tga") == 0) 
    err = loadTGAFile(vd);

  else if (vvToolshed::strCompare(suffix, "tif") == 0 || 
           vvToolshed::strCompare(suffix, "tiff") == 0)
    err = loadTIFFile(vd);

  else if (vvToolshed::strCompare(suffix, "fro") == 0 ||    // VHD CT data
           vvToolshed::strCompare(suffix, "fre") == 0)
    err = loadVHDCTFile(vd);

  else if (vvToolshed::strCompare(suffix, "pd") == 0 ||     // VHD MRI data
           vvToolshed::strCompare(suffix, "t1") == 0 ||
           vvToolshed::strCompare(suffix, "t2") == 0 ||
           vvToolshed::strCompare(suffix, "loc") == 0)
    err = loadVHDMRIFile(vd);

  else if (vvToolshed::strCompare(suffix, "rgb") == 0)      // SGI RGB file
    err = loadRGBFile(vd);

  else if (vvToolshed::strCompare(suffix, "pgm") == 0 ||    // PGM file
           vvToolshed::strCompare(suffix, "ppm") == 0)      // PPM file
    err = loadPXMRawImage(vd);

  else if (vvToolshed::strCompare(suffix, "raw") == 0)      // VHD anatomic
    err = loadVHDAnatomicFile(vd);

  else if (vvToolshed::strCompare(suffix, "dat") == 0)      // DAT file = raw volume data w/o header information
    err = loadRawFile(vd);

  else if (vvToolshed::strCompare(suffix, "dcm") == 0 ||    // DICOM file
           vvToolshed::strCompare(suffix, "dcom") == 0)
    err = loadDicomFile(vd);

  else if (vvToolshed::strCompare(suffix, "vmr") == 0)      // VMR file = BrainVoyager anatomical 3D data
    err = loadVMRFile(vd);

  else if (vvToolshed::strCompare(suffix, "vtc") == 0)      // VTC file = BrainVoyager functional data (time series)
    err = loadVTCFile(vd);

  else if (vvToolshed::strCompare(suffix, "nrd") == 0 ||    // NRRD file = Teem nrrd volume file
           vvToolshed::strCompare(suffix, "nrrd") == 0)    
    err = loadNrrdFile(vd);

  else if (vvToolshed::strCompare(suffix, "ximg") == 0)     // XIMG = General Electric MRI file
    err = loadXIMGFile(vd);

  // Unknown extension error:
  else
  {
    vvDebugMsg::msg(1, "Cannot load volume: unknown extension. File name:", vd->getFilename());
    err = PARAM_ERROR;
  }
  
  return err;
}

//----------------------------------------------------------------------------
/** Computes a volume data set with given sizes
  @param b  byte per voxel
*/
vvFileIO::ErrorType vvFileIO::computeDefaultVolume(vvVolDesc* vd, int b)
{
  const int DEFAULT_WIDTH  = 32;
  const int DEFAULT_HEIGHT = 32;
  const int DEFAULT_SLICES = 32;
  const int DEFAULT_FRAMES = 8;
  const int EQUATION = 0;             // 0 is default
  int x, y, z, f;     // counters
  uchar* rd;          // raw volume data
	int index;
  
  vvDebugMsg::msg(1, "vvFileIO::computeDefaultVolume()");

  vd->vox[0] = DEFAULT_WIDTH;
  vd->vox[1] = DEFAULT_HEIGHT;
  vd->vox[2] = DEFAULT_SLICES;
  vd->setFilename("default.xvf");
  vd->bpv = b;
  vd->dt = 0.1f;

  for (f=0; f<DEFAULT_FRAMES; ++f)
  {
    rd = new uchar[vd->getFrameSize()];
    for (z=0; z<vd->vox[2]; ++z)
      for (y=0; y<vd->vox[1]; ++y)
        for (x=0; x<vd->vox[0]; ++x)
				{
					switch (vd->bpv)
					{
						case 1:
              switch (EQUATION)
              {
                case 0:
                default:
		              rd[x + y * vd->vox[0] + z * vd->getSliceSize()] = 
  		              (uchar)((x*(y+1)*(z+1)*(f+1)) % 256);
                  break;
                case 1:   
                  rd[x + y * vd->vox[0] + z * vd->getSliceSize()] = (uchar)((z * 8) & 0xFF);
                  break;
                case 2:
                  rd[x + y * vd->vox[0] + z * vd->getSliceSize()] = (uchar)0;
                  break;
              }
							break;
						case 4:
						default:
							index = x * vd->bpv + y * vd->vox[0] * vd->bpv + z * vd->getSliceSize();
		          rd[index]   = (uchar)((x*(y+1)*(z+1)*(f+1)) % 256);
		          rd[index+1] = (uchar)(x & 0xFF);
		          rd[index+2] = (uchar)(y & 0xFF);
		          rd[index+3] = (uchar)(x & 0xFF);
						  break;
					}
				}
    vd->addFrame(rd, vvVolDesc::DELETE_DATA);
    ++vd->frames;
  }
  return OK;
}

//----------------------------------------------------------------------------
/** Set compression mode for data compression in files.
  This parameter is only used if the file type supports it.
  @param newCompression true = compression on
*/
void vvFileIO::setCompression(bool newCompression)
{
  compression = newCompression;
}

//============================================================================
// End of File
//============================================================================
