//****************************************************************************
// Filename:          vvxvid.h.cpp
// Author:            Michael Poehnl
// Institution:       University of Stuttgart, Supercomputing Center
// History:           08-01-2002  Creation date
//****************************************************************************

#ifndef _VVXVID_H_
#define _VVXVID_H_

#ifdef VV_XVID
#include "xvid.h"

/**This class is the interface to the XviD library (xvidcore-0.9.0)
   It is used by the vvImage class for the XviD encoding of RGB frames. <BR>

   @author Michael Poehnl
*/
class vvXviD
{
public:      
        vvXviD( float fr=25.0f, int min_q=1, int max_q=31, int br=900, int max_k=250);
        ~vvXviD();
        int create_enc(int w, int h);
        int create_dec(int w, int h);
        int enc_frame(unsigned char* src, unsigned char* dst, int* enc_size, int* key, int style, int quant);
        int dec_frame(unsigned char* src, unsigned char* dst, int src_size, int* dst_size);
        void set_framerate(float fr);
        void set_quantizer(int min_q, int max_q);
        void set_bitrate(int br);
        void set_max_key_interval(int max_k);
        
private:
        float framerate;
        int min_quantizer;
        int max_quantizer;
        int bitrate;
        int max_key_interval; 
        int stride;
        void *enc_handle;
        void *dec_handle;
        int del_enc();
        int del_dec();   
        
        
};

#endif
#endif   

