//****************************************************************************
// Filename:          vvimage.cpp
// Author:            Michael Poehnl
// Institution:       University of Stuttgart, Supercomputing Center
// History:           08-01-2002  Creation date
//****************************************************************************
#ifdef _STANDARD_C_PLUS_PLUS
#include <iostream>
using std::cerr;
using std::endl;
#else
#include <iostream.h>
#endif

#include "vvimage.h"
#include "vvdebugmsg.h"
#ifdef VV_XVID
#include "vvxvid.h"
#endif

//----------------------------------------------------------------------------
/** Constructor for initialization with an image
    @param h   picture height
    @param w   picture width 
    @param image   pointer to the image
*/
vvImage::vvImage(short h, short w, uchar* image)
        : height(h), width(w), imageptr(image)
{   
        vvDebugMsg::msg(3, "vvImage::vvImage(): ", w, h);
        xvidsize = 0;
        size = height*width*4;   
        codetype = 0;
        codedimage = new uchar[size];
        xvidimageptr = new uchar[width*height*6];
        xvidcodedimage = new uchar[width*height*6]; 
        tmpimage = new uchar[width*height];
        t = VV_SERVER;
        xvidstyle = 0;
        xvidquant = 1; 
#ifdef VV_XVID
        xvidEncoder = new vvXviD();
        xvidEncoder->create_enc(w, h);
        xvidDecoder = new vvXviD();
        xvidDecoder->create_dec(w, h);
#endif
}


//----------------------------------------------------------------------------
/** Constructor for an empty image
*/
vvImage::vvImage()
{
        height = 0;
        width = 0;
        size = 0;
        xvidsize = 0;
        imageptr = 0;
        codedimage = 0;
        xvidimageptr = 0;
        xvidcodedimage = 0;
        tmpimage = 0;
        codetype = 0;
        t = VV_CLIENT;
        xvidstyle = 0;
        xvidquant = 1;  
#ifdef VV_XVID
        xvidEncoder = new vvXviD();
        xvidDecoder = new vvXviD();  
#endif
}



//----------------------------------------------------------------------------
/** Destructor
*/ 
vvImage::~vvImage()
{
        if (t == VV_CLIENT)
        {
                if (imageptr == codedimage)         
                {
                        if (imageptr != 0)
                                delete[] imageptr;
                }
                else 
                {
                        if (imageptr != 0)
                                delete[] imageptr;
                        if (codedimage != 0)     
                                delete[] codedimage;
                }
        }
        else     
                delete[] codedimage;
        if (xvidimageptr !=0)
                delete[] xvidimageptr;  
        if (xvidcodedimage != 0)
                delete[] xvidcodedimage;
        if (tmpimage != 0)
                delete[] tmpimage;   
#ifdef VV_XVID
        if (xvidEncoder != 0)
                delete xvidEncoder;
        if (xvidDecoder != 0)
                delete xvidDecoder;
#endif
}

//----------------------------------------------------------------------------
/**Reinitializes an image object with new width and height
@param h   picture height
@param w   picture width 
@param image   pointer to the image
*/
void vvImage::setNewImage(short h, short w, uchar* image)
{
        vvDebugMsg::msg(3, "vvImage::setNewImage(): ", w, h);
#ifdef VV_XVID
        if (width != w || height != h)
        {
                xvidEncoder->create_enc(w, h);
                xvidDecoder->create_dec(w, h);
        }
#endif
        height =h;
        width = w;
        imageptr = image;
        size = height*width*4;
        codetype = 0;
        if (codedimage != 0)
                delete[] codedimage;
        codedimage = new uchar[size];  
        if (xvidimageptr !=0)
                delete [] xvidimageptr;
        if (xvidcodedimage != 0)
                delete [] xvidcodedimage; 
        if (tmpimage != 0)
                delete [] tmpimage;     
        xvidimageptr = new uchar[width*height*6];
        xvidcodedimage = new uchar[width*height*6];
        tmpimage = new uchar[width*height];  
}

//----------------------------------------------------------------------------
/**Sets the image pointer to a new image which has the same height and
width as the old one.
@param image   pointer to the image
*/
void vvImage::setNewImagePtr(uchar* image)
{
        imageptr = image;
        size = height*width*4;
        codetype = 0;  
}


//----------------------------------------------------------------------------
/**Encodes an image
@param ct   codetype to use (see detailed description of the class)
@param sw   start pixel relating to width
@param ew   end pixel relating to width
@param sh   start pixel relating to height
@param eh   end pixel relating to height
@return size of encoded image in bytes, or -1 on error
*/
int vvImage::encode(short ct, short sw, short ew, short sh, short eh)
{
        short realheight, realwidth;
        int start;
        float cr;
        
        if (size <= 0)
        {
                vvDebugMsg::msg(1, "Illegal image parameters ");
                return -1;
        }   
        switch(ct)
        {
        case 0:cr=1;break;
        case 1:
        {
                if (spec_RLC_encode(0, height, width))
                {        
                        vvDebugMsg::msg(1, "No compression possible");
                        codetype = 0;
                }
                else
                        codetype = 1;
                cr = (float)size / (height*width*4);
        }break;
        case 2:
        {
                codetype = 2;  
                if(sh<0 || eh<0  || sw<0 || ew<0 ||
                   (realheight=short(eh-sh+1))<=0 || (realwidth=short(ew-sw+1))<=0 ||
                   eh > height-1 || ew > width-1)
                {
                        vvDebugMsg::msg(1,"Wrong usage vvImage::encode()");
                        return -1;
                }
                start = (sh)*width*4 + (sw)*4;
                vvToolshed::write32(&codedimage[0],(ulong)start);
                vvToolshed::write16(&codedimage[4],realwidth);
                if (spec_RLC_encode(start, realheight, realwidth, 6))
                {        
                        vvDebugMsg::msg(1,"No compression possible");
                        codetype = 0;
                } 
                cr = (float)size / (height*width*4); 
        }break;    
#ifdef VV_XVID
        case 3:
        {
                int i;
                codetype = 3;
                for (i=0; i<width*height; ++i)
                        memcpy(&xvidimageptr[i * 3], &imageptr[i * 4], 3); 
                for (i=0; i<width*height; ++i)
                        memcpy(&tmpimage[i], &imageptr[i * 4 +3], 1);                 
                imageptr = tmpimage; 
                if (xvidEncode())
                {
                        vvDebugMsg::msg(1,"Error: xvidEncode()");
                        return -1;
                }     
                if ( (size = gen_RLC_encode(imageptr, codedimage, width*height, 1, width*height*4)) < 0)
                {
                        vvDebugMsg::msg(1,"Error: gen_RLC_encode()");
                        return -1;
                }
                imageptr = codedimage;       
                cr = (float)(size+xvidsize) / (height*width*4);            
        }break;
#endif
        default:
                vvDebugMsg::msg(1,"Unknown encoding type ", ct );
                return -1;
        }   
        vvDebugMsg::msg(2, "compression rate: ", cr);
        vvDebugMsg::msg(3, "image encoding succeeded");
        return size;
}

//----------------------------------------------------------------------------
/** Decodes an image
*/
int vvImage::decode()
{
        short  realwidth;
        int start;
        
        switch(codetype)
        {
        case 0: imageptr = codedimage;break;
        case 1:
        {
                spec_RLC_decode(0, width);
        }break;
        case 2: 
        {
                memset(imageptr, 0, height*width*4);
                start = (int)vvToolshed::read32(&codedimage[0]);
                realwidth = vvToolshed::read16(&codedimage[4]);
                spec_RLC_decode(start, realwidth, 6);        
        }break;
#ifdef VV_XVID
        case 3:
        {       
                int i;
                if (xvidDecode())
                {
                        vvDebugMsg::msg(1,"Error: xvidDecode()");
                        return -1;
                }
                for (i=0; i<width*height; ++i)
                        memcpy(&imageptr[i * 4], &xvidimageptr[i * 3], 3);              
                if (gen_RLC_decode(codedimage, tmpimage, size, 1, width*height))
                {
                        vvDebugMsg::msg(1,"Error: gen_RLC_decode()");
                        return -1;
                }         
                for (i=0; i<width*height; ++i)
                        memcpy(&imageptr[i * 4 + 3], &tmpimage[i], 1); 
                size = width*height*4;
        }break;
#endif
        default:
                vvDebugMsg::msg(1,"No encoding type with that identifier");
                return -1;
        }         
        codetype = 0;
        vvDebugMsg::msg(3, "image decoding succeeded");
        return 0;
}
//----------------------------------------------------------------------------
/** Sets the image height.
*/
void vvImage::setHeight(short h)
{
        height = h;
}

//----------------------------------------------------------------------------
/** Sets the image width.
*/
void vvImage::setWidth(short w)
{
        width = w;
}

//----------------------------------------------------------------------------
/** Sets the code type.
*/
void vvImage::setCodeType(short ct)
{
        codetype = ct;
}

//----------------------------------------------------------------------------
/** Sets the image size.
*/
void vvImage::setSize(int s)
{
        size = s;
}

//----------------------------------------------------------------------------
/** Sets the XviD image size.
*/
void vvImage::setXviDSize(int s)
{
        xvidsize = s;
}

//----------------------------------------------------------------------------
/** Sets the image pointer.
*/
void vvImage::setImagePtr(uchar* image)
{
        imageptr = image;    
}

//----------------------------------------------------------------------------
/** Sets the style of Xvid encoding
*/
void vvImage::setXvidStyle(int s)
{
        if ( (s<0) || (s>6) )
        {
                vvDebugMsg::msg(1, "XvidStyle hast to be between 0 and 6, using 0 now");
                xvidstyle = 0;
        }
        else
                xvidstyle = s;
}

//----------------------------------------------------------------------------
/** Sets the value for the Xvid quantizer
*/
void vvImage::setXvidQuant(int q)
{
        if ( (q<1) || (q>31) )
        {
                vvDebugMsg::msg(1,"XvidQuant has to be between 1 and 31, using 1 now");
                xvidquant = 1;
        }
        else
                xvidquant = q;  
}
//----------------------------------------------------------------------------
/** Sets an key frame
@param k 
*/
void vvImage::setKeyframe(int k)
{
        keyframe = k;
}

//----------------------------------------------------------------------------
/**Returns the image height
*/
short vvImage::getHeight()
{
        return height;
}

//----------------------------------------------------------------------------
/** Returns the image width
 */
short vvImage::getWidth()
{
        return width;
}

//----------------------------------------------------------------------------
/** Returns the code type
 */
short vvImage::getCodeType()
{
        return codetype;
}

//----------------------------------------------------------------------------
/** Returns the image size in bytes
*/
int vvImage::getSize()
{
        return size;
}

//----------------------------------------------------------------------------
/** Returns the Xvid image size in bytes
*/
int vvImage::getXviDSize()
{
        return xvidsize;
}
//----------------------------------------------------------------------------
/** Returns the key frame
*/
int vvImage::getKeyframe()
{
        return keyframe;
}

//----------------------------------------------------------------------------
/** Returns the pointer to the image
*/
uchar* vvImage::getImagePtr()
{
        return imageptr;
}

//----------------------------------------------------------------------------
/** Returns the pointer to the encoded image
*/
uchar* vvImage::getCodedImage()
{
        return codedimage;
}
//----------------------------------------------------------------------------
/** Returns the pointer to the encoded XviD image
*/
uchar* vvImage::getXviDCodedImage()
{
        return xvidcodedimage;
}

//----------------------------------------------------------------------------
/**Does the Run Length Encoding for a defined cutout of an image. 
@param start   start pixel for RLE encoding
@param h   height of pixel square to encode
@param w   width of pixel square to encode
@param dest   start writing in coded image at position dest
*/
int vvImage::spec_RLC_encode(int start, short h, short w, int dest)
{
        short samePixel=1;
        short diffPixel=0;
        int src;
        int l,m;
        
        for ( int i=0; i < h; i++)
        {
                src = start + i*width*4;   
                for ( int j=0; j < w; j++)
                {    
                        if (j == (w-1))
                        {   
                                l=1;
                                if (i == (h-1))
                                        m=0;
                                else
                                        m =1;
                        }  
                        else
                        {
                                m=1;
                                l=0;
                        }     
                        if (imageptr[src] == imageptr[m*(src+4+l*(width-w)*4)] &&   
                            imageptr[src+1] == imageptr[m*(src+5+l*(width-w)*4)] && 
                            imageptr[src+2] == imageptr[m*(src+6+l*(width-w)*4)] &&
                            imageptr[src+3] == imageptr[m*(src+7+l*(width-w)*4)] )  
                        {             
                                if(samePixel == 129)
                                        put_same(samePixel, dest);
                                else
                                {
                                        samePixel++;
                                        if(diffPixel > 0 )
                                                put_diff(diffPixel, dest);             
                                        if(samePixel == 2)
                                        {
                                                if ((dest+5) > size)
                                                        return -1;
                                                memcpy(&codedimage[dest+1], &imageptr[src], 4);       
                                        }               
                                }
                        }
                        else
                        {
                                if (samePixel > 1)
                                        put_same(samePixel, dest);                  
                                else 
                                {
                                        if ((dest+5+4*diffPixel) > size)
                                                return -1;           
                                        memcpy(&codedimage[dest+1+diffPixel*4], &imageptr[src], 4);           
                                        diffPixel++;
                                        if(diffPixel == 128)      
                                                put_diff(diffPixel, dest);                   
                                }           
                        }
                        src += 4;
                }
        }
        if (samePixel > 1)
        {
                samePixel--;
                put_same(samePixel, dest);
        }
        else if (diffPixel > 0)
                put_diff(diffPixel, dest);
        imageptr = codedimage;
        size = dest;
        return 0;
}
//----------------------------------------------------------------------------
/** Does the Run Length Decoding for a cutout of an image
@param start   start pixel where the decoded pixel square is
written
@param w   width of pixel square to decode
@param src   start position of encoded pixels in coded image
*/
int vvImage::spec_RLC_decode(int start, short w, int src)
{
        int dest;
        int length;
        
        dest = start;
        while (src < size)
        {
                length = (int)codedimage[src];
                if (length > 127)
                {
                        for(int i=0; i<(length - 126); i++)
                        {
                                if (((dest-start-4*w)% (4*width)) == 0 && dest != start)
                                        dest += (width-w)*4;
                                memcpy(&imageptr[dest], &codedimage[src+1], 4);
                                dest += 4;    
                        }
                        src += 5;
                }
                else
                {
                        length++;
                        for(int i=0; i<(length); i++)
                        {
                                if (((dest-start-4*w)% (4*width)) == 0 && dest != start)
                                        dest += (width-w)*4;	     
                                memcpy(&imageptr[dest], &codedimage[src+1+i*4], 4);      	 
                                dest +=4;
                        }
                        src += 1+4*length;
                }
        }
        size = height*width*4;
        return 0;
}


//----------------------------------------------------------------------------
/** Writes a RLE encoded set of same pixels.
@param sP   number of same pixels
@param d   destination in coded image where to write
*/
void vvImage::put_same(short& sP, int& d)
{
        codedimage[d] = (uchar)(126+sP);      
        d += 5;
        sP = 1;
}   

//----------------------------------------------------------------------------
/** Writes a RLE encoded set of different pixels.
@param dP   number of different pixels
@param d   destination in coded image where to write
*/
void vvImage::put_diff(short& dP, int& d)
{
        codedimage[d] = (uchar)(dP-1);
        d += 1+4*dP;
        dP=0;
}   

//----------------------------------------------------------------------------
/** Allocates momory for a new image
*/
int vvImage::alloc_mem()
{
        vvDebugMsg::msg(3, "vvImage::alloc_mem(): ", width, height);
        
#ifdef VV_XVID
        xvidEncoder->create_enc(width, height);
        xvidDecoder->create_dec(width, height);
#endif

        if (imageptr == codedimage)         
        { 
                if (imageptr != 0)
                        delete[] imageptr;
        }
        else 
        {
                if (imageptr != 0)
                        delete[] imageptr;
                if (codedimage != 0)
                        delete[] codedimage;   
        }
        if (xvidimageptr !=0)
                delete [] xvidimageptr;
        if (xvidcodedimage != 0)
                delete [] xvidcodedimage;     
        if (tmpimage != 0)
                delete [] tmpimage;       
        if (codetype != 0)
        {
                imageptr = new uchar[height*width*4];
                if (!imageptr)
                        return -1;
        }
        if (codetype == 3)
        {
                xvidimageptr = new uchar[height*width*6];
                if (!xvidimageptr)
                        return -1;        
                xvidcodedimage = new uchar[height*width*6];
                if (!xvidcodedimage)
                        return -1;      
                tmpimage = new uchar[height*width];
                if (!tmpimage)
                        return -1;                       
        }          
        codedimage = new uchar[height*width*4];
        if (!codedimage)
                return -1;   
        
        return 0;
}

//----------------------------------------------------------------------------
/** Does the Xvid encoding
 */
int vvImage::xvidEncode()
{
#ifdef VV_XVID
        assert(xvidEncoder); 
        if (xvidEncoder->enc_frame(xvidimageptr, xvidcodedimage, &xvidsize, &keyframe, xvidstyle, xvidquant))
        {
                vvDebugMsg::msg(1,"Error xvidEncode()");
                return -1;
        }
        vvDebugMsg::msg(3, "encoded Xvid image size: ", xvidsize);
#endif
        return 0;
}

//----------------------------------------------------------------------------
/** Does the Xvid decoding
 */
int vvImage::xvidDecode()
{
#ifdef VV_XVID
        int newsize;
        
        assert(xvidDecoder);
        if (xvidDecoder->dec_frame(xvidcodedimage, xvidimageptr,  xvidsize, &newsize))
        {
                vvDebugMsg::msg(1,"Error xvidDecode()");
                return -1;
        }
#endif
        return 0;
}

//----------------------------------------------------------------------------
/** general function for the RLC encoding
*/
int vvImage::gen_RLC_encode(uchar* in, uchar* out, int size, int symbol_size, int space) // size=total size in byte
{
        int same_symbol=1;
        int diff_symbol=0;
        int src=0;
        int dest=0;
        bool same;
        int i;
        
        if ((size % symbol_size) != 0)
        {
                vvDebugMsg::msg(1,"No RLC encoding possible with this parameters");
                return -1;
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
                                out[dest] = (uchar)(126+same_symbol);      
                                dest += symbol_size+1;
                                same_symbol = 1;  
                        }
                        else
                        {
                                same_symbol++;
                                if (diff_symbol > 0)
                                {
                                        out[dest] = (uchar)(diff_symbol-1);
                                        dest += 1+symbol_size*diff_symbol;
                                        diff_symbol=0;  
                                }
                                if (same_symbol == 2)
                                {
                                        if ((dest+1+symbol_size) > space)
                                        {
                                                vvDebugMsg::msg(1,"Not enough memory to encode");
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
                                out[dest] = (uchar)(126+same_symbol);      
                                dest += symbol_size+1;
                                same_symbol = 1;           
                        }
                        else
                        {
                                if ((dest+1+diff_symbol*symbol_size+symbol_size) > space)
                                {
                                        vvDebugMsg::msg(1,"Not enough memory to encode");
                                        return -1;
                                }                              
                                memcpy(&out[dest+1+diff_symbol*symbol_size], &in[src], symbol_size);           
                                diff_symbol++;
                                if (diff_symbol == 128)
                                {
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
                out[dest] = (uchar)(126+same_symbol);      
                dest += symbol_size+1;
        }
        else
        {
                if ((dest+1+diff_symbol*symbol_size+symbol_size) > space)
                {
                        vvDebugMsg::msg(1,"Not enough memory to encode");
                        return -1;
                }                   
                memcpy(&out[dest+1+diff_symbol*symbol_size], &in[src], symbol_size);           
                diff_symbol++;            
                out[dest] = (uchar)(diff_symbol-1);
                dest += 1+symbol_size*diff_symbol;
        }
        if (dest > size)
        {
                vvDebugMsg::msg(1,"No compression possible with RLC !!!");
        }
        return dest;  
}   

//----------------------------------------------------------------------------
/** general function for the RLC decoding
*/   
int vvImage::gen_RLC_decode(uchar* in, uchar* out, int size, int symbol_size, int space)
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
                                        vvDebugMsg::msg(1,"Not enough memory to decode");
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
                                vvDebugMsg::msg(1,"Not enough memory to decode");
                                return -1;
                        }                                     
                        memcpy(&out[dest], &in[src+1], symbol_size*length);
                        dest += length*symbol_size; 	         
                        src += 1+symbol_size*length;
                }
        }
        return 0;
}
