#ifndef _IMG_RENDER_ACTION_
#define _IMG_RENDER_ACTION_

// own headers

#include "imginfo.h"
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/fields/SoSFBool.h>
#include <stdio.h>
#include <GL/gl.h>

typedef unsigned char *ucptr;
typedef unsigned char *ImgRenderActionCB(const ImageInfo &info, int UserData);

/* declaration of extern variables */ 
extern unsigned char TwoDSimulation; /* Boolean for 2d-Simulation*/

// ---------------------------------------------------------------

class ImgViewer;

class ImgRenderAction : public SoGLRenderAction
{ private:
     const SoSFBool    *Stereo;
     SoSFBool    *leftImg;
     const SoSFBool    *interactive;
     ImgViewer         *viewer;
     ImgRenderActionCB *CallBack;
     unsigned char     *Image[2];
     int               CallBackArg;

  public:
    ImgRenderAction(ImgViewer              *view,
		    SoSFBool               *Ste,
		    SoSFBool               *left,
		    SoSFBool               *interac,
		    const SbViewportRegion &viewportRegion)/*, 
                    SbBool                 useCurrentGLValues=FALSE)*/
     : SoGLRenderAction(viewportRegion/*,useCurrentGLValues*/) 
       { Stereo=Ste;
         leftImg=left;
	 interactive=interac;
	 viewer=view;
	 Image[0]=NULL;
	 Image[1]=NULL;
	 CallBack=NULL;}
     void DrawImage(void);
     void getImgSize(int *image_width, int *image_height);
     void DrawImage2D(void);
     void RotateImage(int axis, float angle, unsigned char *image);
     void setImgRenderCB(ImgRenderActionCB *cb, int cba=0)
      { CallBack=cb; CallBackArg=cba; }
     virtual ~ImgRenderAction() {;}
    
    virtual void	apply(SoNode *node);
    //  {
	//if (interactive->getValue()){  
          //if ((global.viewer->DrawStyle() == VIEW_SAME_AS_STILL) && (!TwoDSimulation)){
            //DrawImage();
          //}
          //else{
            //SoGLRenderAction::apply(node);
          //}
        //}
	//else{
          //if(TwoDSimulation){
            //DrawImage2D(); // SoGLRenderAction::apply(node); 
          //}
          //else{
            //DrawImage(); // SoGLRenderAction::apply(node); 
          //}
	//}
      //}

    virtual void	apply(SoPath *path)
      { SoGLRenderAction::apply(path);
	if (interactive->getValue()){  
	  SoGLRenderAction::apply(path);
        }
	else{
	  DrawImage();
        }
      }

    virtual void	apply(const SoPathList &pathList,
			      SbBool obeysRules = FALSE)
       {
	 if (interactive->getValue()){
	   SoGLRenderAction::apply(pathList,obeysRules);
         }
	 else{
	   DrawImage();
         }
       }
};


#endif
