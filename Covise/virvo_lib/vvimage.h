//****************************************************************************
// Filename:          vvimage.h
// Author:            Michael Poehnl
// Institution:       University of Stuttgart, Supercomputing Center
// History:           08-01-2002  Creation date
//****************************************************************************


#ifndef _VVIMAGE_H_
#define _VVIMAGE_H_

#include <string.h>
#include <assert.h>
#include "vvtoolshed.h"
#ifdef VV_XVID
#include "vvxvid.h"
#endif
//----------------------------------------------------------------------------  
/**This class provides different encoding and decoding types for RGB images. <BR>

Supported code types:
- no encoding (code type 0)
- Run Lenght Encoding over the whole image (code type 1)
- Run Lenght Encoding over a quadratic part of the image (code type 2).
  Therefore start and end pixels for width an height must be specified.
  The rest of the image is interpreted as background and the pixels get the
  value 0,0,0,0. (picture width from 0 - width-1, picture height from 0 - height-1)
- Xvid Encoding (code type 3). For this type the VV_XVID Flag must be set.<BR>

Here is an example code fragment for encoding and decoding an image with
800 x 600 pixels :<BR>
<PRE>

//Create a new image class instance 
vvImage* im = new vvImage(600, 800, (char *)imagepointer);

//Encode with normal RLE
if(im->encode(1) < 0)
{   
    delete im;
    return -1;
}

//Or encode with RLE but only the lower half of the image
if(im->encode(2, 0, 799, 300, 599 ) < 0)
{   
    delete im;
    return -1;
}
//Decode the image 
if(im->decode())
    return -1;

delete im;
</PRE> 
*/
class vvImage
{
public:
        
        vvImage(short, short, uchar*);
        vvImage();
        virtual ~vvImage();      
        int encode(short, short sh=-1, short eh=-1, short sw=-1, short ew=-1);
        int decode();
        void setNewImage(short, short, uchar*); 
        void setHeight(short);
        void setWidth(short);
        void setCodeType(short);
        void setSize(int);
        void setXviDSize(int);
        void setImagePtr(uchar*);
        void setKeyframe(int);
        void setNewImagePtr(uchar*);
        void setXvidStyle(int);
        void setXvidQuant(int);
        short getHeight();
        short getWidth();
        short getCodeType();
        int getSize();
        int getXviDSize();
        int getKeyframe();
        uchar* getImagePtr();    
        uchar* getCodedImage();
        uchar* getXviDCodedImage();
        int alloc_mem();
        
private:
        
        enum Type
        {
                VV_SERVER,
                VV_CLIENT    
        };
        
        Type t;  
        short height;
        short width;
        short codetype;
        int size;
        int xvidsize;
        int keyframe;
        uchar* imageptr;
        uchar* codedimage;
        uchar* xvidimageptr;
        uchar* xvidcodedimage;  
        uchar* tmpimage;  
        int xvidstyle;
        int xvidquant;
#ifdef VV_XVID
        vvXviD* xvidEncoder;
        vvXviD* xvidDecoder;
#endif
        
        int spec_RLC_encode(int, short, short, int dest=0);
        int spec_RLC_decode(int, short, int src=0);
        int gen_RLC_encode(uchar*, uchar*, int, int, int);
        int gen_RLC_decode(uchar*, uchar*, int, int, int);
        void put_diff(short&, int&);
        void put_same(short&, int&);
        int xvidEncode();
        int xvidDecode();
}; 



#endif
