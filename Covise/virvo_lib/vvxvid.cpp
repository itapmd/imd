//****************************************************************************
// Filename:          vvxvid.cpp
// Author:            Michael Poehnl
// Institution:       University of Stuttgart, Supercomputing Center
// History:           19-12-2002  Creation date
//****************************************************************************
#ifdef VV_XVID

#include "vvxvid.h"
#include "vvdebugmsg.h"
#ifndef NULL
#ifdef __GNUG__
#define NULL (__null)
#else
#define NULL (0)
#endif
#endif

//----------------------------------------------------------------------------
/** Destructor
*/
static int const motion_presets[7] = {
	0,                                                        // Q 0
	PMV_EARLYSTOP16,                                          // Q 1
	PMV_EARLYSTOP16,                                          // Q 2
	PMV_EARLYSTOP16 | PMV_HALFPELREFINE16,                    // Q 3
	PMV_EARLYSTOP16 | PMV_HALFPELREFINE16,                    // Q 4
	PMV_EARLYSTOP16 | PMV_HALFPELREFINE16 | PMV_EARLYSTOP8 |  // Q 5
	PMV_HALFPELREFINE8,
	PMV_EARLYSTOP16 | PMV_HALFPELREFINE16 | PMV_EXTSEARCH16 | // Q 6
	PMV_USESQUARES16 | PMV_EARLYSTOP8 | PMV_HALFPELREFINE8
};

//----------------------------------------------------------------------------
/** Destructor
*/
static int const general_presets[7] = {
	XVID_H263QUANT,	                              // Q 0
	XVID_MPEGQUANT,                               // Q 1
	XVID_H263QUANT,                               // Q 2
	XVID_H263QUANT | XVID_HALFPEL,                // Q 3
	XVID_H263QUANT | XVID_HALFPEL | XVID_INTER4V, // Q 4
        XVID_H263QUANT | XVID_HALFPEL | XVID_INTER4V, // Q 5
	XVID_H263QUANT | XVID_HALFPEL | XVID_INTER4V  // Q 6
};

//----------------------------------------------------------------------------
/** Constructor
@param fr  framerate
@param min_q  lower bound for quantizer
@param max_q  upper bound for quantizer
@param br  target bitrate
@param max_k  maximum key interval
*/
vvXviD::vvXviD(float fr, int min_q, int max_q, int br, int max_k)
        :framerate(fr), min_quantizer(min_q), max_quantizer(max_q), bitrate(br), max_key_interval(max_k)
{
        enc_handle = 0;
        dec_handle = 0; 
        stride = 0; 
}

//----------------------------------------------------------------------------
/** Destructor
*/
vvXviD::~vvXviD()
{
        if ( enc_handle != 0)
        {
                if (del_enc())
                        vvDebugMsg::msg(1,"error: del_enc()");
        }
        if (dec_handle != 0)
        {
                if (del_dec())
                        vvDebugMsg::msg(1,"error: del_dec()");
        }   
}

//----------------------------------------------------------------------------
/** Creates an XviD encoder
@param w  width of frames 
@param h  height of frames 
@return   0 for success, != 0 for error
*/
int vvXviD::create_enc(int w, int h)
{   
	int xerr;
	XVID_INIT_PARAM xinit;
	XVID_ENC_PARAM xparam;
	
        if (enc_handle != 0)
                if (del_enc())
                        vvDebugMsg::msg(1,"error: del_enc()");
        xinit.cpu_flags = XVID_CPU_FORCE;
        xvid_init(0, 0, &xinit, 0);
        if (xinit.api_version != API_VERSION)
        {
                vvDebugMsg::msg(1,"Wrong Xvid library version");
                return -1;
        }      
	xparam.width = w;
	xparam.height = h;
	if ((framerate - (int)framerate) < 0.001)
	{
		xparam.fincr = 1;
		xparam.fbase = (int)framerate;
	}
	else
	{
		xparam.fincr = 1000;
		xparam.fbase = (int)(1000 * framerate);
	}
        xparam.rc_reaction_delay_factor = -1;  //default values
        xparam.rc_averaging_period = -1;     //default values
        xparam.rc_buffer = -1;               //default values
        xparam.rc_bitrate = bitrate*1000;
        xparam.min_quantizer = min_quantizer;
        xparam.max_quantizer = max_quantizer;
        xparam.max_key_interval = max_key_interval;
        xerr = xvid_encore(0, XVID_ENC_CREATE, &xparam, 0);
        enc_handle=xparam.handle;
        if (xerr == 0)
                vvDebugMsg::msg(3, "XviD Encoder created");
        return xerr;   
}

//----------------------------------------------------------------------------
/** Creates an XviD decoder
@param w  width of frames
@param h  height of frames 
@return   0 for success, != 0 for error
*/
int vvXviD::create_dec(int w, int h)
{   
	int xerr;
	XVID_INIT_PARAM xinit;
	XVID_DEC_PARAM xparam;  
        
        if (dec_handle != 0)
                if (del_dec())
                        vvDebugMsg::msg(1,"error: del_dec()");
        xinit.cpu_flags = XVID_CPU_FORCE;
        xinit.cpu_flags = 0;
        xvid_init(NULL, 0, &xinit, NULL); 
        if (xinit.api_version != API_VERSION)
        {
                vvDebugMsg::msg(1,"Wrong Xvid library version");
                return -1;
        }      
	xparam.width = w;
        stride = w;
	xparam.height = h; 
        xerr = xvid_decore(NULL, XVID_DEC_CREATE, &xparam, NULL);
        dec_handle = xparam.handle;
        if (xerr == 0)
                vvDebugMsg::msg(3, "XviD Decoder created");   
        return xerr;      
}

//----------------------------------------------------------------------------
/** Encodes a frame
@param src  pointer to frame to encode 
@param dst  pointer to destination
@param enc_size  OUT, size of encoded frame
@param key  OUT, encoded frame is a keyframe if != 0
@param style  style of encoding [0 .. 6] 
@param quant  quantizer value to use for that frame [1 .. 31]
@return   0 for success, != 0 for error
*/
int vvXviD::enc_frame(unsigned char* src, unsigned char* dst, int* enc_size, int* key, int style, int quant)
{
	int xerr;
	XVID_ENC_FRAME xframe;
        
        if (style < 0 || style > 6)
                style = 0; 
        if (quant < 1 || quant > 31)
                quant = 1;    
        xframe.bitstream = dst;
	xframe.length = -1; 	// this is written by the routine
        xframe.image = src;
        xframe.colorspace = XVID_CSP_RGB24;
        xframe.intra = -1;
        xframe.quant = quant;
        xframe.motion = motion_presets[style];
        xframe.general = general_presets[style];
        xframe.quant_intra_matrix = xframe.quant_inter_matrix = 0;
        xerr = xvid_encore(enc_handle, XVID_ENC_ENCODE, &xframe, 0);
        *key = xframe.intra;
        *enc_size = xframe.length;
        if (xerr == 0)
                vvDebugMsg::msg(3, "frame encoded in XviD style");      
        return xerr;
}

//----------------------------------------------------------------------------
/** Decodes a frame
@param src  pointer to encoded frame  
@param dst  pointer to destination
@param src_size  size of encoded frame
@param dst_size  OUT, size of decoded frame
@return   0 for success, != 0 for error
*/
int vvXviD::dec_frame(unsigned char* src, unsigned char* dst, int src_size, int* dst_size)
{
        int xerr;
        XVID_DEC_FRAME xframe;
        
        xframe.bitstream = src;
        xframe.length = src_size;
        xframe.stride = stride;
        xframe.image = dst;
        xframe.colorspace = XVID_CSP_RGB24 | XVID_CSP_VFLIP;
        xerr = xvid_decore(dec_handle, XVID_DEC_DECODE, &xframe, 0);
        *dst_size =  xframe.length;
        if (xerr == 0)
                vvDebugMsg::msg(3, "XviD frame decoded");         
        return xerr;  
}

//----------------------------------------------------------------------------
/** Sets the target framerate 
@param fr  framerate
*/
void vvXviD::set_framerate(float fr)
{
        framerate = fr;
}

//----------------------------------------------------------------------------
/** Sets the quantizer bounds
@param min_q  lower quantizer bound 
@param max_q  upper quantizer bound
*/
void vvXviD::set_quantizer(int min_q, int max_q)
{
        min_quantizer = min_q;
        max_quantizer = max_q;
}

//----------------------------------------------------------------------------
/** Sets the target bitrate
@param br  bitrate
*/
void vvXviD::set_bitrate(int br)
{
        bitrate = br;
}

//----------------------------------------------------------------------------
/** Sets the maximum interval for key frames
@param max_k maximum interval in frames
*/
void vvXviD::set_max_key_interval(int max_k)
{
        max_key_interval = max_k;
}

//----------------------------------------------------------------------------
/** Deletes the encoder
@return   0 for success, != 0 for error
*/
int vvXviD::del_enc()
{
	int xerr;
        
        xerr = xvid_encore(enc_handle, XVID_ENC_DESTROY, 0, 0);
        enc_handle = 0;
        return xerr;
}   

//----------------------------------------------------------------------------
/** Deletes the decoder
@return   0 for success, != 0 for error
*/
int vvXviD::del_dec()
{
        int xerr;
        
        xerr = xvid_decore(dec_handle, XVID_DEC_DESTROY, 0, 0);
        dec_handle = 0;
        return xerr;
}


#endif
