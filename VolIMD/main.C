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
   file name .................  main.C
   author ....................  Roland Niemeier
   version ...................  0.9beta0
   date of last change .......  07/25/96
   description : interface for direct volume visualization of scalar fields. 
                 The application development is intended for 
                 the interactive visualization of volume data sets. 
                 It is based on Open Inventor, VolPack, Xforms.
   purpose of main.C: initialization and call of main loop
*/

// std, X11, Xmotif, Inventor includes

#include <stdlib.h>
#include <X11/Intrinsic.h>
#include <Xm/MainW.h>
#include <Xm/Form.h>
#include <Inventor/Xt/SoXt.h>

// specific includes

#include "ImgExaminerViewer.h"
#include "render.h"
#include "viewer.h"
#include "menu.h"
#include "global.h"


static String fallback[] = 
  {
     "*background: grey70",
     "*fontList:-*-helvetica-bold-r-*-*-12-*=charset", 
     "*geometry: 460x466",
     "*Quit*foreground: red", 
     NULL
  };
  
Widget toplevel;
GlobalInfoTyp global;

// ---------------------------------------------------------------

void forms_initialization(int *argc, char **argv);
void main(int argc, char **argv)

/************************************************************************/
/* Initializes widgets, graphics, connections, etc.			*/
/* enters main loop							*/
/************************************************************************/
{
   XtAppContext app_context;
   Arg args[1];
   
   // ----- Initialize X Application and create a form
   toplevel=XtVaAppInitialize(&app_context,"MyClass",NULL,0,
                                   &argc,argv,fallback,NULL,0);
   SoXt::init(toplevel);

   Widget form = XtCreateManagedWidget("form",xmFormWidgetClass,toplevel,args,0);
   Widget MenuBar = createMenuBar(form, &global);

   forms_initialization(&argc, argv);
   initialize_volpack();
   initViewer(global.viewer,form,MenuBar);
   SoXt::show(form);
   SoXt::show(toplevel);
   SoXt::mainLoop();      // Main Inventor event loop
}

