/*                         
                           (C) 1997
              Computer Centre University of Stuttgart
                         Allmandring 30
                       D-70550 Stuttgart
                            Germany

THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,
EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY
WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

Permission to use, copy, modify and distribute this software is hereby 
granted without fee, provided that the above and this permission notice 
appear in all copies of this software and that the authorship is 
acknowledged.

*/
/* program name ..............  VolRend  
   file name .................  ImgRenderAction.C
   author ....................  Roland Niemeier
   version ...................  0.9beta0
   date of last change .......  07/27/97
   description : interface for direct volume visualization of scalar fields. 
                 The application development is intended for 
                 the interactive visualization of volume data sets. 
                 It is based on Open Inventor, volpack, xforms.
   purpose of ImgRenderAction.C: 
                 gets volume and camera values, makes a callback
                 to the renderer and displays the image(s)
*/

// Open Inventor, GL and math includes
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/SbLinear.h>
#include <Inventor/nodes/SoTransform.h>
#include <GL/gl.h>
#include <math.h>

// specific includes 
#include "ImgRenderAction.h"	// defines class ImgRenderAction
#include "ImgViewer.h" 		// defines class ImgViewer
#include "ImgExaminerViewer.h" 	// defines class ImgExaminerViewer
#include "global.h" 		// GlobalInfoTyp structure
unsigned char *get_image_pointer(void);

// declaration of extern variables
extern Boolean DrawFlag; // Boolean for Rendering/Drawing if a volume is loaded
/* two color modes: grayscale or rgb */
enum ColorModes {
LUMINANCE_MODE = 0,
RGB_MODE = 1
};
extern ColorModes CurrentColorMode; 
extern GlobalInfoTyp global; // pointer to viewer, etc. 



// definition of file local and global variables
const int IMAGE_HEIGHT_OFFSET=32; // for down decoration bar
const int IMAGE_WIDTH_OFFSET=60;  // for left and right decoration bars
extern const int COLOR_CHANNELS = 3;
extern unsigned char FORCEDRAW = 0;
static SbVec3f translation; // translation vector 

// definition of functions (methods) 

//virtual void apply(SoNode *node)
void ImgRenderAction::apply(SoNode *node)
{
  if (interactive->getValue()){  
    if ((global.viewer->getDrawStyle(ImgViewer::INTERACTIVE) == ImgViewer::VIEW_SAME_AS_STILL) && (!TwoDSimulation)){
      DrawImage();
    }
    else{
      SoGLRenderAction::apply(node);
    }
  }
  else{
    if(TwoDSimulation){
      DrawImage2D(); // SoGLRenderAction::apply(node); 
    }
    else{
      DrawImage(); // SoGLRenderAction::apply(node); 
    }
  }
}

void ImgRenderAction::getImgSize(int *SizeX, int *SizeY)
/************************************************************************/
/* Returns the current size of the image		                */
/************************************************************************/
{
    SbVec2s Size=viewer->getSize();
    *SizeX=(int)Size[0];
    *SizeY=(int)Size[1];
    if (global.viewer->isDecoration()){
      *SizeX=*SizeX-IMAGE_WIDTH_OFFSET;
      *SizeY=*SizeY-IMAGE_HEIGHT_OFFSET;
    }
}
void ImgRenderAction::DrawImage2D(void)
/************************************************************************/
/* Draws the image of a 2D Simulation					*/
/************************************************************************/
{
    SbVec2s Size=viewer->getSize();
    int SizeX=Size[0];
    int SizeY=Size[1];
    glDisable(GL_DEPTH_TEST);
    glViewport(0, 0, SizeX, SizeY);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, SizeX, 0, SizeY, -1e30, 1e30);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glRasterPos2i(0, 0);
    if (global.viewer->isDecoration()){
      SizeX=SizeX-IMAGE_WIDTH_OFFSET;
      SizeY=SizeY-IMAGE_HEIGHT_OFFSET;
    }
    printf("SizeX and SizeY: %d %d\n",SizeX,SizeY);
    Image[0] = get_image_pointer();
    if (DrawFlag)
    {
        switch (CurrentColorMode){
          case LUMINANCE_MODE: 
            printf("Drawing greyscale Image\n"); 
            glDrawPixels(SizeX, SizeY,
                GL_LUMINANCE,GL_UNSIGNED_BYTE, (const GLvoid *)Image[0]);
            break;
          case RGB_MODE:
            printf("Drawing rgb Image\n"); 
            glDrawPixels(SizeX, SizeY,
                GL_RGB,GL_UNSIGNED_BYTE, (const GLvoid *)Image[0]);
            break;
          default: fprintf(stderr,"Error: Unknown CurrentColorMode");
            exit(1);
        }
    }
}

void ImgRenderAction::DrawImage(void)
/************************************************************************/
/* Makes a CallBack to the volume renderer if				*/
/* camera position, object orientation or image size is changed 	*/
/* otherwise redraws the image						*/
/************************************************************************/
{
    /* viewing matrices (stereo) from previous call */
    static SbMatrix oldmat[2];  
    /* eye position from previous call */
    static SbVec3f  oldeye[2];  
    /* sizes of images from previous call */
    static float oldSizeX, oldSizeY; 
    /* structure to store characteristics for rendering */
    static ImageInfo ImgInfo;
    /* current camera position and orientation */
    SoCamera *camera=global.viewer->getCamera();
    SbRotation rotation, scale_rotation;
    SbVec3f scale_translation;
    /* viewing matrices */
    SbMatrix mat[2], projection_mat;
    SbViewVolume view_volume = camera->getViewVolume(0.0);
    /* camera eye position */
    SbVec3f eye[2]; 
    /* image width and height */
    SbVec2s Size=viewer->getSize();
    int SizeX=Size[0];
    int SizeY=Size[1];

    /* if stereo is used, which image should be drawn: left or right?  */
    int PicNo = (Stereo->getValue()) ? leftImg->getValue() : 0 ;

    /* store the current camera position into eye[PicNo] */
    eye[PicNo] =camera->position.getValue();
    //mat[PicNo] =camera->orientation.getValue();
    /* store the viewing matrices */
    view_volume.getMatrices(mat[PicNo],projection_mat);

    /* select translation and rotation parts of viewing matrix */
    mat[PicNo].getTransform(translation,rotation,
                     scale_translation,scale_rotation);

    /* variables for starting a new rendering if decoration is turned on/off */
    static SbBool old_decoration_flag = global.viewer->isDecoration();
    SbBool redraw_if_decoration_changes;
    if (old_decoration_flag != global.viewer->isDecoration()){
      old_decoration_flag = global.viewer->isDecoration();
      redraw_if_decoration_changes = TRUE;
    }
    else {
      redraw_if_decoration_changes = FALSE;
    }
    if (global.viewer->isDecoration()){
      ImgInfo.SizeX=SizeX-IMAGE_WIDTH_OFFSET;
      ImgInfo.SizeY=SizeY-IMAGE_HEIGHT_OFFSET;
    }
    else{
      ImgInfo.SizeX=SizeX;
      ImgInfo.SizeY=SizeY;
    }


    /* perspective or orthographic camera */
    if (camera->isOfType(SoPerspectiveCamera::getClassTypeId())){
      global.Camera_Type = 1;
    }
    else global.Camera_Type = 0;


    if (    (eye[PicNo]!=oldeye[PicNo])|| (mat[PicNo]!=oldmat[PicNo])
	      || (SizeX!=oldSizeX) || (SizeY!=oldSizeY) 
              || redraw_if_decoration_changes || FORCEDRAW) 
      {
        //printf("needs a new rendering\n");
        FORCEDRAW = 0;
        /* prepare ImgInfo for Callback */
	ImgInfo.perspect[0] = mat[PicNo][0][0];
	ImgInfo.perspect[1] = mat[PicNo][0][1];
	ImgInfo.perspect[2] = mat[PicNo][0][2];
	ImgInfo.perspect[3] = 0.0;
	ImgInfo.perspect[4] = mat[PicNo][1][0];
	ImgInfo.perspect[5] = mat[PicNo][1][1];
	ImgInfo.perspect[6] = mat[PicNo][1][2];
	ImgInfo.perspect[7] = 0.0;
	ImgInfo.perspect[8] = mat[PicNo][2][0];
	ImgInfo.perspect[9] = mat[PicNo][2][1];
	ImgInfo.perspect[10] = mat[PicNo][2][2];
	ImgInfo.perspect[11] = 0.0;
	ImgInfo.perspect[12] = translation[0];
	ImgInfo.perspect[13] = translation[1];
	ImgInfo.perspect[14] = translation[2];
	ImgInfo.perspect[15] = 1.0;
	ImgInfo.eye[0] = eye[PicNo][0];
	ImgInfo.eye[1] = eye[PicNo][1];
	ImgInfo.eye[2] = eye[PicNo][2];
	


        /* calling the Renderer */
	if (CallBack  && DrawFlag)
	  Image[PicNo]=CallBack(ImgInfo,PicNo);
//	  Image[PicNo]=CallBack(ImgInfo,CallBackArg);

        /* copy eye positions and viewing matrices for next call */
	oldeye[PicNo]=eye[PicNo];
	oldmat[PicNo]=mat[PicNo];
	oldSizeX=SizeX;
	oldSizeY=SizeY;

        glDisable(GL_DEPTH_TEST);
        glViewport(0, 0, SizeX, SizeY);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0, SizeX, 0, SizeY, -1e30, 1e30);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        glRasterPos2i(0, 0);
      }

    if (Image[PicNo] && DrawFlag)
      /* draw the image */
      {
        //printf("Drawing picture number %d\n",PicNo);
        switch (CurrentColorMode){
          case LUMINANCE_MODE: 
            glDrawPixels(ImgInfo.SizeX, ImgInfo.SizeY,
                GL_LUMINANCE,GL_UNSIGNED_BYTE, (const GLvoid *)Image[PicNo]);
            break;
          case RGB_MODE: 
            glDrawPixels(ImgInfo.SizeX, ImgInfo.SizeY,
                GL_RGB,GL_UNSIGNED_BYTE, (const GLvoid *)Image[PicNo]);
            break;
          default: fprintf(stderr,"Error: Unknown CurrentColorMode");
            exit(1);
        }
       // some of the following statements may be useful to switch buffers for
       // a new rendering or to force a rendering
       //global.viewer->setBufferingType(DOUBLE);
       //glDrawBuffer(GL_FRONT);
       //glDrawBuffer(global.viewer->isDoubleBuffer() ? GL_BACK : GL_FRONT);
       //global.viewer->getSceneManager()->render(FALSE, FALSE);
       //if (global.viewer->isDoubleBuffer()){
       //  Display *display = global.viewer->getDisplay();
       //  Drawable d = DefaultRootWindow(display);
       //  glXSwapBuffers(display,d);}
      }
}

void ImgRenderAction::RotateImage(int axis, float angle, unsigned char *)
//void ImgRenderAction::RotateImage(int axis, float angle, unsigned char *image)
/************************************************************************/
/* rotates camera and calls renderer					*/
/* if stereo for both images 						*/
/* before calling ImgRenderAction::DrawImage for drawing		*/
/************************************************************************/
{

  /* for the use of axis vectors */
  SbVec3f x_axis, y_axis, z_axis;
  x_axis[0] = 1; x_axis[1] = 0; x_axis[2] = 0;
  y_axis[0] = 0; y_axis[1] = 1; y_axis[2] = 0;
  z_axis[0] = 0; z_axis[1] = 0; z_axis[2] = 1;

  /* rotation increment */
  angle=(angle/180.0)*acos(-1);

  /* rotation matrices around x, y, or z-axis */
  SbRotation rot_mat[2];
  switch (axis){
    case 1: 
      rot_mat[0].setValue(x_axis,angle);
      break;
    case 2:
      rot_mat[0].setValue(y_axis,angle);
      break;
    case 3:
      // different conventions of z-axis in volpack and inventor 
      rot_mat[0].setValue(z_axis,-angle);
      break;
    default:
      fprintf(stderr,"Error: switch default in RotateImage !?\n");
    exit(1);
  }

  /* get the old camera orientation and position */
  SoCamera *camera=global.viewer->getCamera();
  SbRotation rot = camera->orientation.getValue();
  SbVec3f eye = camera->position.getValue(), rot_eye;

  /* set the new rotated camera orientation and position */
  rot_mat[0] *= rot;
  rot_mat[0].multVec(translation,rot_eye);
  rot_eye = -rot_eye;
  camera->orientation.setValue(rot_mat[0]);
  camera->position.setValue(rot_eye);
  global.viewer->setCamera(camera);

  /* calculate the left and right eye positions if stereo is active */
  if (Stereo->getValue()){
    rot_mat[1] = rot_mat[0];
    // rotate the right eye camera by + stereoOffset/2 degrees
    rot_mat[0] = SbRotation(SbVec3f(0, 1, 0), + global.viewer->getStereoOffset() * M_PI / 360.0) * rot_mat[0];
    // rotate the left eye camera by - stereoOffset/2 degrees
    rot_mat[1] *= rot;
    rot_mat[1] = SbRotation(SbVec3f(0, 1, 0), - global.viewer->getStereoOffset() * M_PI / 360.0) * rot_mat[1];
  }

  /* prepare for new renderings */
  ImageInfo ImgInfo;
  SbVec2s Size=viewer->getSize();
  SbMatrix Matrot;
  int SizeX=Size[0], SizeY=Size[1];
  if (global.viewer->isDecoration()){
    ImgInfo.SizeX=SizeX-IMAGE_WIDTH_OFFSET;
    ImgInfo.SizeY=SizeY-IMAGE_HEIGHT_OFFSET;
  }
  else{
    ImgInfo.SizeX=SizeX;
    ImgInfo.SizeY=SizeY;
  }
    rot_mat[0].getValue(Matrot);
    ImgInfo.perspect[0] = Matrot[0][0];
    ImgInfo.perspect[1] = Matrot[1][0];
    ImgInfo.perspect[2] = Matrot[2][0];
    ImgInfo.perspect[3] = 0.0;
    ImgInfo.perspect[4] = Matrot[0][1];
    ImgInfo.perspect[5] = Matrot[1][1];
    ImgInfo.perspect[6] = Matrot[2][1];
    ImgInfo.perspect[7] = 0.0;
    ImgInfo.perspect[8] = Matrot[0][2];
    ImgInfo.perspect[9] = Matrot[1][2];
    ImgInfo.perspect[10] = Matrot[2][2];
    ImgInfo.perspect[11] = 0.0;
    ImgInfo.perspect[12] = -translation[0];
    ImgInfo.perspect[13] = -translation[1];
    ImgInfo.perspect[14] = -translation[2];
    ImgInfo.perspect[15] = 1.0;
    ImgInfo.eye[0] = rot_eye[0];
    ImgInfo.eye[1] = rot_eye[1];
    ImgInfo.eye[2] = rot_eye[2];

  if (CallBack){
    Image[0]=CallBack(ImgInfo,0);
      if (Stereo->getValue()){
	rot_mat[1].multVec(translation,rot_eye);
	rot_eye = -rot_eye;
	rot_mat[1].getValue(Matrot);
	ImgInfo.perspect[0] = Matrot[0][0];
	ImgInfo.perspect[1] = Matrot[1][0];
	ImgInfo.perspect[2] = Matrot[2][0];
	ImgInfo.perspect[3] = 0.0;
	ImgInfo.perspect[4] = Matrot[0][1];
	ImgInfo.perspect[5] = Matrot[1][1];
	ImgInfo.perspect[6] = Matrot[2][1];
	ImgInfo.perspect[7] = 0.0;
	ImgInfo.perspect[8] = Matrot[0][2];
	ImgInfo.perspect[9] = Matrot[1][2];
	ImgInfo.perspect[10] = Matrot[2][2];
	ImgInfo.perspect[11] = 0.0;
	ImgInfo.perspect[12] = -translation[0];
	ImgInfo.perspect[13] = -translation[1];
	ImgInfo.perspect[14] = -translation[2];
	ImgInfo.perspect[15] = 1.0;
	ImgInfo.eye[0] = rot_eye[0];
	ImgInfo.eye[1] = rot_eye[1];
	ImgInfo.eye[2] = rot_eye[2];
        Image[1]=CallBack(ImgInfo,1);
      }
 }

 /* finally draw the already rendered images */
 if (DrawFlag){
   global.viewer->actualRedraw();
 }
}
