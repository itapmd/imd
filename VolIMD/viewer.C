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
   file name .................  viewer.C
   author ....................  Roland Niemeier
   version ...................  0.9beta0
   date of last change .......  07/25/97
   description : interface for direct volume visualization of scalar fields. 
                 The application development is intended for 
                 the interactive visualization of volume data sets. 
                 It is based on Open Inventor, VolPack, Xforms.
   purpose of viewer.C: Collection of routines for the Open Inventor viewer 
*/

// Nodes
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoRotation.h>
#include <Inventor/nodes/SoScale.h>
#include <Inventor/nodes/SoTransformSeparator.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/nodes/SoOrthographicCamera.h>

// Basic types
#include <Inventor/SbLinear.h>

// Fields
#include <Inventor/fields/SoSFRotation.h>
#include <Inventor/fields/SoSFVec3f.h>
#include <Inventor/fields/SoSFImage.h>

// Actions and other Open Inventor tools

#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/manips/SoTrackballManip.h>
#include <Inventor/Xt/SoXtComponent.h>
#include <Inventor/Xt/SoXt.h>
#include <Inventor/SoDB.h>
#include <Inventor/SoInput.h>

#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/Xt/SoXtDirectionalLightEditor.h>
#include <Inventor/SoPath.h>
#include <Inventor/draggers/SoTranslate1Dragger.h>
#include <Inventor/draggers/SoTranslate2Dragger.h>


// Basic C includes

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <iostream.h>
#include <new.h>
#include <bstring.h>
#include <fstream.h>
#include <stdlib.h>
#include <string.h>

// specific includes 

#include "ImgExaminerViewer.h"
#include "ImgRenderAction.h"
#include "global.h"
#include "render.h"

SoSeparator *root;
SoDirectionalLight *light_source = new SoDirectionalLight;
extern unsigned char FORCEDRAW;
// ---------------------------------------------------------------

void WindowCloseCB(void *, SoXtComponent *)
/************************************************************************/
/* Exit all					       			*/
/* if stereo is used set monitor back  		                        */
/************************************************************************/
{
   //system("/usr/gfx/setmon 72HZ");
   exit(1);
}
// ---------------------------------------------------------------

void SetLightCB(void *,const SoDirectionalLight *directional_light)
/************************************************************************/
/* Creates an Open Inventor light editor and attaches it to 		*/
/* the current_light							*/
/************************************************************************/
{
SbVec3f direction = directional_light->direction.getValue();
float intensity = directional_light->intensity.getValue();
SbVec3f color = directional_light->color.getValue();
RendererSetLightCB((double)color[0]*intensity,
                   (double)color[1]*intensity,
                   (double)color[2]*intensity, 
                   (double)direction[0],
                   (double)direction[1],
                  -(double)direction[2]);
}

// ---------------------------------------------------------------

void ModifyLightCB(double r_light, double g_light, double b_light,
                   double view_x, double view_y, double view_z)
/************************************************************************/
/* Creates an Open Inventor light editor and attaches it to 		*/
/* the current_light							*/
/************************************************************************/
{
   /* initialize the light direction and color */
   root->ref();
   
//   SoDirectionalLight *light_source = new SoDirectionalLight;
   root->addChild(light_source);
   view_z = -view_z;
   light_source->direction.setValue(view_x, view_y, view_z);
   printf("Setting light direction to %f %f %f\n",view_x, view_y, view_z);
   light_source->color.setValue(r_light, g_light, b_light);

   SoPath *light_path = new SoPath;
   light_path->append(light_source);

   SoXtDirectionalLightEditor *ltEditor = new SoXtDirectionalLightEditor();
   ltEditor->attach(light_path);
   ltEditor->setTitle("Light Editor");
   ltEditor->show();
   ltEditor->addLightChangedCallback(SetLightCB,0);
}
// ---------------------------------------------------------------

ImgExaminerViewer *CreateViewer(Widget myWindow)
/************************************************************************/
/* Create the viewer							*/
/************************************************************************/
{
   ImgExaminerViewer *myViewer;
   myViewer = new ImgExaminerViewer(myWindow,"Viewer");
   myViewer->setDoubleBuffer(FALSE);
   myViewer->setBufferingType(ImgExaminerViewer::BUFFER_INTERACTIVE);
   myViewer->setFeedbackVisibility(FALSE);   
   myViewer->setStereoViewing(FALSE);
//   system("/usr/gfx/setmon STR_TOP");
   
   myViewer->setSize(SbVec2s(400,400)); 
   //myViewer->setDecoration(FALSE);
   //myViewer->setPopupMenuEnabled(FALSE);
   //myViewer->setTitle("VolVR");
   myViewer->setWindowCloseCallback(WindowCloseCB,NULL);
   //myViewer->setBackgroundColor(SbColor(1.0,1.0,1.0));   
   return myViewer;
}

// ---------------------------------------------------------------

//SoSeparator *buildRoot(SoNode *scene, SoPerspectiveCamera *&camera)
SoSeparator *buildRoot(SoNode *scene, SoOrthographicCamera *&camera)
/************************************************************************/
/* create some open inventor nodes   					*/
/************************************************************************/
{
  SoSeparator *root = new SoSeparator;
  root->ref();
  
  SoTransformSeparator *transSep = new SoTransformSeparator;
  root->addChild(transSep);
  
  SoScale *scale = new SoScale;
  scale->scaleFactor.setValue(1.0,1.0,1.0);
  transSep->addChild(scale);
  //root->addChild(new SoTrackballManip);

//  camera = new SoPerspectiveCamera;
  camera = new SoOrthographicCamera;
  transSep->addChild(camera);
  
  root->addChild(scene);
  
  root->unrefNoDelete();
  return root;
}


// ---------------------------------------------------------------

SoSeparator *buildObject()
/************************************************************************/
/* add a trackball for interactive mode					*/
/************************************************************************/
{
   SoSeparator *root = new SoSeparator;
   root->ref();
   root->addChild(new SoTrackballManip);
   root->addChild(new SoDirectionalLight);
   /*
   // for drawing some arrows along the coordinate axis
   SoSeparator *xaxisSeparator = new SoSeparator;
   SoSeparator *yaxisSeparator = new SoSeparator;
   SoSeparator *zaxisSeparator = new SoSeparator;
   SoTransform *xaxisTransform = new SoTransform;
   SoTransform *yaxisTransform = new SoTransform;
   SoTransform *zaxisTransform = new SoTransform;
   SoMaterial  *xaxisMaterial  = new SoMaterial;
   SoMaterial  *yaxisMaterial  = new SoMaterial;
   SoMaterial  *zaxisMaterial  = new SoMaterial;
   yaxisTransform->rotation.setValue(0,0,1,1.);
   zaxisTransform->rotation.setValue(0,1,0,1.);
   xaxisSeparator->addChild(xaxisTransform);
   yaxisSeparator->addChild(yaxisTransform);
   zaxisSeparator->addChild(zaxisTransform);
   xaxisMaterial->diffuseColor.setValue(1.0, 0.0, 0.0);   // Red
   xaxisMaterial->specularColor.setValue(1.0, 0.0, 0.0);   // Red
   xaxisMaterial->ambientColor.setValue(1.0, 0.0, 0.0);   // Red
   yaxisMaterial->diffuseColor.setValue(0.0, 1.0, 0.0);   // Green
   yaxisMaterial->specularColor.setValue(0.0, 1.0, 0.0);   // Green
   yaxisMaterial->ambientColor.setValue(0.0, 1.0, 0.0);   // Green
   zaxisMaterial->diffuseColor.setValue(0.0, 0.0, 1.0);   // Blue
   zaxisMaterial->specularColor.setValue(0.0, 0.0, 1.0);   // Blue
   zaxisMaterial->ambientColor.setValue(0.0, 0.0, 1.0);   // Blue
   xaxisSeparator->addChild(xaxisMaterial);
   yaxisSeparator->addChild(yaxisMaterial);
   zaxisSeparator->addChild(zaxisMaterial);
   SoTranslate1Dragger *xaxis = new SoTranslate1Dragger;
   SoTranslate1Dragger *yaxis = new SoTranslate1Dragger;
   SoTranslate1Dragger *zaxis = new SoTranslate1Dragger;
   xaxisSeparator->addChild(xaxis);
   yaxisSeparator->addChild(yaxis);
   zaxisSeparator->addChild(zaxis);

   //yaxis->translation.setValue(0.0,0.0,1.0);
   root->addChild(xaxisSeparator);
   root->addChild(yaxisSeparator);
   root->addChild(zaxisSeparator);
   */
   root->unrefNoDelete();
   return root;
}

// ---------------------------------------------------------------

void initViewer(ImgExaminerViewer* &myViewer, Widget form, Widget menu_bar)
/************************************************************************/
/* is called by main, initializes viewers characteristics		*/
/************************************************************************/
{
   // Initialize Inventor. This returns a main window to use.
   // If unsuccessful, exit.


   //SoTransform *cameraTransform = new SoTransform;
   //SoPerspectiveCamera *camera;
   SoOrthographicCamera *camera;
   
   myViewer=CreateViewer(form);
   
   //SoSeparator 
   root = buildRoot(buildObject(),camera);
   root->ref();
   //root->addChild(cameraTransform);
   myViewer->setSceneGraph(root);
   myViewer->setCamera(camera);
   myViewer->viewAll();
   myViewer->setDrawStyle(ImgExaminerViewer::INTERACTIVE,
                          ImgExaminerViewer::VIEW_LINE);

   myViewer->getImgRenderAction()->setImgRenderCB(ImgRenderCB);

   // place the adjustable top menu bar
   Arg args[6];
   int n = 0;
   XtSetArg(args[n], XmNtopAttachment, XmATTACH_WIDGET); n++;
   XtSetArg(args[n], XmNtopWidget, menu_bar); n++;
   XtSetArg(args[n], XmNleftAttachment, XmATTACH_FORM); n++;
   XtSetArg(args[n], XmNrightAttachment, XmATTACH_FORM); n++;
   XtSetArg(args[n], XmNbottomAttachment, XmATTACH_FORM); n++;
   XtSetValues(myViewer->getWidget(), args, n);

   // show the viewer
   myViewer->show();
   // initialize HomePosition
   myViewer->saveHomePosition();
}

// ---------------------------------------------------------------
