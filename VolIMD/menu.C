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
   file name .................  menu.C
   author ....................  Roland Niemeier
   version ...................  0.9beta0
   date of last change .......  07/25/97
   description : interface for direct volume visualization of scalar fields. 
                 The application development is intended for 
                 the interactive visualization of volume data sets. 
                 It is based on Open Inventor, volpack, xforms.
   purpose of menu.C: Collection of menu bar callback functions and popup menus
                      that do not use xforms. 
*/
// Xstuff
#include <X11/Intrinsic.h>
#include <Xm/DialogS.h>
#include <Xm/Xm.h>
#include <X11/StringDefs.h>
#include <Xm/Label.h>
#include <Xm/Form.h>
#include <Xm/PushB.h>
#include <Xm/PushBG.h>
#include <Xm/RowColumn.h>
#include <Xm/CascadeB.h>
#include <Xm/CascadeBG.h>
#include <Xm/MessageB.h>
#include <Xm/ToggleB.h>
#include <Xm/ToggleBG.h>
#include <Xm/SeparatoG.h>
#include <Xm/Text.h>
#include <Xm/FileSB.h>
// Inventor stuff
//#include <Inventor/Xt/SoXt.h>
//#include <Inventor/Xt/viewers/SoXtExaminerViewer.h>
//#include <Inventor/nodes/SoCone.h>
// standard unix
#include <stdio.h>
#include <stdlib.h>
// for starting rendering functions
#include "ImgRenderAction.h"
#include "global.h"
#include "render.h"
#include "volpack.h"
#include "tokens.h"
#include "guiforms.h"
#include "client.h";

// declaration of extern variables

extern vpContext  *vpc;	/* rendering context */
extern Widget toplevel; /* toplevel widget */
extern GlobalInfoTyp global; /* viewer pointer, etc. */
extern int endianByteSwap; /* FG - this is an int */
enum ColorModes { /* either grayscale or rgb mode */
  LUMINANCE_MODE = 0,
  RGB_MODE = 1
};
extern ColorModes CurrentColorMode;
/* number of rotations if menu rotate is selected */
extern int number_of_rotations; 
/* angle increment between two pictures if menu rotate is selected */
extern double rotation_increment;
extern unsigned char KineticEnergyFile;
// definition of variables for local (routines in this file) and extern usage
extern unsigned char TwoDSimulation=0;
static unsigned char HistogramDisplayed = 0;
static unsigned char portOnThisMachine = 1;
extern unsigned char usePreviousLimits = 0;
// local structures

enum MenuButtonTypes {
  M_SEPARATOR,
  M_PUSH_BUTTON,
  M_TOGGLE_BUTTON,
  M_RADIO_BUTTON
};

typedef struct
{
  int number;
  WidgetList widget;
}
RadioButtonGrTyp;

typedef struct
{
   Widget widget;
   int id;
   float userData;
   GlobalInfoTyp *globalInfo;
} ButtonDataTyp;

typedef struct
{
  char	         *name;
  char           *resname;
  int		  id;
  int		  buttonType;	// PUSH, TOGGLE, RADIO
  XtCallbackProc  callback;     // who to call ?
  float 	  userData;
  char		 *accelerator;	// e.g. "Alt <Key> p" or "Ctrl <Key> u"
  char		 *accelText;	// text that appears in the menu item
} MenuButtonTyp;

typedef struct
{
    char		*name;
    int			 id;
    MenuButtonTyp	*subMenu;
    int			 subItemCount;
    XtCallbackProc	 callback;
} MenuTyp;

// functions defined in this file

void FileCB (Widget,XtPointer,XtPointer);
void ModifyCB(Widget,XtPointer,XtPointer);
void RotateCB(Widget,XtPointer,XtPointer);
void SwitchCB  (Widget,XtPointer,XtPointer);
void SimulationCB(Widget , XtPointer client_data, XtPointer);
void Create_FileBox(int identifier);
/* motif file selection boxes */
void open_fileCB  (Widget, int, XmSelectionBoxCallbackStruct *);
void openClassified_fileCB(Widget, int type, XmSelectionBoxCallbackStruct *cbs);
void save_imageCB  (Widget, int, XmSelectionBoxCallbackStruct *);

// local variables

unsigned char file_selected = 0; /* boolean for file selection */
char *filename; /* name of selected file */
Widget FselBox; /* motif file selection box widget*/
enum FILEBOX {OKAY, CANCEL}; /* only two buttons for File Selection Box */

// =================================================================

#include "menu.inc"
/************************************************************************/
/* Defines layout and CallBack Functions of menu bar			*/
/************************************************************************/

// =================================================================

void FileCB(Widget,
            XtPointer client_data, /*normally an integer used to decide
                                 which menu was selected. */
            XtPointer)
/************************************************************************/
/* Callback for PullDownMenu File					*/
/************************************************************************/
{
  ButtonDataTyp *button = (ButtonDataTyp *) client_data;
// Case Open and classify ...
  if (button->id == M_FILE_OPEN)
    {
      Create_FileBox(M_FILE_OPEN);
    }
// Case Open classified ...
  if (button->id == M_FILE_OPEN_CLASSIFIED)
    {
      Create_FileBox(M_FILE_OPEN_CLASSIFIED);
    }
// Case Build Octree ...
  if (button->id == M_FILE_BUILD_OCTREE)
    {
      Create_FileBox(M_FILE_OPEN);

      if (file_selected){
        change_filenames(filename);
        build_raw_volume();      
        build_octree();
        // and render new volume (rotation by degree 0 is used)
        multi_rotate_image(1,0.0,1);
        file_selected = 0;
      }      
    }
// Case Build Raw Volume  ...
  if (button->id == M_FILE_BUILD_RAW_VOLUME)
    {
      build_raw_volume();      
    }
// Case Save image ...
  if (button->id == M_FILE_SAVE_IMAGE)
    {
      Create_FileBox(M_FILE_SAVE_IMAGE);
    }
// Case Select image sequence name ... for storing a sequence of images
  if (button->id == M_FILE_SELECT_SEQUENCE)
    {
      Create_FileBox(M_FILE_SELECT_SEQUENCE);
    }
// Case Open classified volume sequence ... for loading a sequence of classified volumes
  if (button->id == M_FILE_LOAD_SEQUENCE)
    {
      Create_FileBox(M_FILE_LOAD_SEQUENCE);
    }
// Case Open unclassified volume sequence ... for loading a sequence of unclassified volumes
  if (button->id == M_FILE_UNCLASSIFIED_SEQUENCE)
    {
      Create_FileBox(M_FILE_UNCLASSIFIED_SEQUENCE);
    }
// Case Quit
  if (button->id == M_FILE_QUIT) 
    {
      //sendQuit(button->globalInfo->soc);
      exit(0);
    }
}

// ==================================================================

void ModifyCB(Widget,XtPointer client_data,XtPointer)
/************************************************************************/
/* Callback for PullDownMenu Modify					*/
/* Modifies render characteristics with sliders, xyplots, etc.		*/
/* Cases: Material, Lighting, Classification, Rotation Increment	*/
/************************************************************************/
{
  ButtonDataTyp *button = (ButtonDataTyp *) client_data;
  int option = (int) button->id;
  //int option = (int) button->userData;
  //sendModify(button->globalInfo->soc,option);
  switch (option){
    case M_MODIFY_MATERIAL_PROPERTIES:
      modify_material_properties();
      break;
    case M_MODIFY_MATERIAL_NUMBERS:
      modify_material_numbers();
      break;
    case M_MODIFY_LIGHT_PROPERTIES:
      modify_light_properties();
      break;
    case M_MODIFY_LIGHT_NUMBERS:
      modify_light_numbers();
      break;
    case M_MODIFY_CLASSIFICATION:
      // classification case popups an active xyplot and needs xforms
      activate_xyplot();
      break;
    case M_MODIFY_HISTOGRAM_BOUNDARIES:
      // for modifying the rotation parameters (angle increment & rotation number)
      activate_histogram_slider();
      break;
    case M_MODIFY_ROTATION_INCREMENT:
      // for modifying the rotation parameters (angle increment & rotation number)
      activate_rotation_slider();
      break;
    case M_MODIFY_SCALING:
      // for modifying the scaling values (relative distances between voxels) 
      activate_scaling_sliders();
      break;
    case M_MODIFY_DEPTH_CUEING:
      // for modifying the depth cueing parameters front factor and fog density
      activate_depth_cueing_slider();
      break;
    case M_MODIFY_MINIMUM_OPACITY:
      // for modifying the minimum opacity
      activate_minimum_opacity_slider();
      break;
    case M_MODIFY_MAXIMUM_OPACITY:
      // for modifying the maximum opacity
      activate_maximum_opacity_slider();
      break;
    default:
      fprintf(stderr,"Bug in ModifyCB\n");
      exit(1);
  }
}

// ==================================================================

void SwitchCB(Widget , XtPointer client_data, XtPointer)
/************************************************************************/
/* Callback for PullDownMenu Switch					*/
/* Toggles some volpack state variables				 	*/
/* further action until now not implemented				*/
/************************************************************************/
	
{
  int vp_option, done = 0;
  ButtonDataTyp *button = (ButtonDataTyp *) client_data;
  switch ((int)button->id) {
    case M_SWITCH_SHADOW:
      vp_option = VP_SHADOW;
      vpToggle(vp_option);
      fprintf(stderr,"Sorry, not implemented!\n");
      break;
    case M_SWITCH_DEPTH_CUE:
      vp_option = VP_DEPTH_CUE;
      vpToggle(vp_option);
      break;
    case M_SWITCH_LIGHT_BOTH:
      vp_option = VP_LIGHT_BOTH_SIDES;
      vpToggle(vp_option);
      multi_rotate_image(1,0.0,1);
      break;
    case M_SWITCH_REVERSE_SURFACE:
      vp_option = VP_REVERSE_SURFACE_SIDES;
      vpToggle(vp_option);
      fprintf(stderr,"Sorry, not implemented!\n");
      break;
    case M_SWITCH_COLOR_MODE:
      switch_color_mode();
      done = 1;
      break;
    case M_SWITCH_USE_OCTREE:
      switch_octree_mode();
      done = 1;
      break;
    case M_SWITCH_SAVE_SEQUENCE:
      switch_image_store_mode();
      done = 1;
      break;
    case M_SWITCH_EDITOR:
      switch_light_editor();
      done = 1;
      break;
    case M_SWITCH_FILENAME:
      KineticEnergyFile = (KineticEnergyFile) ? 0 : 1;
      done = 1;
      break;
    case M_SWITCH_TWOD:
      TwoDSimulation = (TwoDSimulation) ? 0 : 1;
      done = 1;
      break;
    case M_SWITCH_HISTOGRAM:
      if (HistogramDisplayed){
        remove_histogram();
      }
      else{
        calculate_histogram();
      }
      HistogramDisplayed = HistogramDisplayed ? 0 : 1;
      break;
    default:
      printf("Case default in SwitchCB\n");
      break;
  }
}
// ==================================================================

void SimulationCB(Widget , XtPointer client_data, XtPointer)
/************************************************************************/
/* Callback for PullDownMenu Simulation					*/
/* contact the simulation and get some data from simulation	 	*/
/************************************************************************/
{
  ButtonDataTyp *button = (ButtonDataTyp *) client_data;
  char token;
  switch ((int)button->id) {
    case M_SIM_GET_DISTRIBUTION:
      token = T_DISTRIBUTION;
      if (portOnThisMachine){
        connect_client(token);
      }
      else{
        connect_server(token);
      }
      break;
    case M_SIM_GET_PICTURE:
      token = T_PICTURE;
      if (portOnThisMachine){
        connect_client(token);
      }
      else{
        connect_server(token);
      }
      break;
    case M_SIM_TOGGLE_STORAGE:
      token = T_TOGGLE_STORAGE;
      printf("deferred implementation ...\n");
      break;
    case M_SIM_CHANGE_PORT:
      token = T_PORT_ADDRESS;
      ask_port_address();
      //printf("deferred implementation ...\n");
      break;
    case M_SIM_CHANGE_NAME:
      ask_server_name();
      break;
    case M_SIM_PORT_MACHINE:
      portOnThisMachine = portOnThisMachine ? 0 : 1;
      break;
    case M_SIM_TOGGLE_LIMITS:
      usePreviousLimits = usePreviousLimits ? 0 : 1; 
      break;
    case M_SIM_BYTE_SWAP:
      endianByteSwap = endianByteSwap ? 0 : 1; 
      break;
    case M_SIM_TERMINATE:
      token = T_QUIT;
      connect_client(token);
      break;
    default:
      printf("Case default in SimulationCB !?\n");
      break;
  }
}
	
// ==================================================================

void RotateCB(Widget, XtPointer client_data, XtPointer)
/************************************************************************/
/* Callback for PullDownMenu Rotate					*/
/* Rotate the viewing coordinates around the axis and render volume 	*/
/* Cases: Around x-axis, y-axis, z-axis					*/
/************************************************************************/
{
  ButtonDataTyp *button = (ButtonDataTyp *) client_data;
  int option = (int) button->userData;
  // option: 1 = x-axis, 2 = y-axis, 3 = z-axis;
  multi_rotate_image(option, (float)rotation_increment, number_of_rotations);
}

// ==================================================================

void openClassified_fileCB(Widget, int type, XmSelectionBoxCallbackStruct *cbs)
/************************************************************************/
/* CB initiated by open in file pulldown menu				*/
/* file selection box ( OK, FILTER, CANCEL)	                        */
/* until now only implementation for loading classified volumes		*/
/************************************************************************/
{
//Widget dialog;
//unsigned char *image;
char *filename;

    switch(type) {
	case OKAY:
	    XmStringGetLtoR(cbs->value, XmSTRING_DEFAULT_CHARSET, 
			 &filename);
            printf("file %s selected\n",filename);
            //dialog = create_dialog(FselBox);
            //XtManageChild(dialog);
            load_classified_volume(filename);
            // and render new volume (rotation by degree 0.000001 is used
            // because display of volume data only if camera position 
            // is changed)
	    multi_rotate_image(1,0.000001,1);
            //image = rotate_image(1, 0.0);
            //global.viewer->getImgRenderAction()->RotateImage(1,0,image);
	    break;
	case CANCEL:
	    break;
    }
    XtUnmanageChild(FselBox);
}


void open_fileCB(Widget, int type, XmSelectionBoxCallbackStruct *cbs)
/************************************************************************/
/* CB initiated by open in file pulldown menu				*/
/* file selection box ( OK, FILTER, CANCEL)	                        */
/* until now only implementation for loading classified volumes		*/
/************************************************************************/
{
//Widget dialog;
//unsigned char *image;

    switch(type) {
	case OKAY:
	    XmStringGetLtoR(cbs->value, XmSTRING_DEFAULT_CHARSET, 
			 &filename);
            printf("file %s selected\n",filename);
            //dialog = create_dialog(FselBox);
            //XtManageChild(dialog);
            XtUnmanageChild(FselBox);
            XtDestroyWidget(FselBox);
            file_selected = 1;
            if (file_selected){
            load_and_classify_volume(filename);
            //and render new volume (rotation by degree 0 is used)
            multi_rotate_image(1,0.0,1);
            file_selected = 0;
            }
	    break;
	case CANCEL:
            XtUnmanageChild(FselBox);
	    break;
    }
 //   XtUnmanageChild(FselBox);
}

// ==================================================================

void save_imageCB(Widget, int type, XmSelectionBoxCallbackStruct *cbs)
/************************************************************************/
/* CB initiated by save in file pulldown menu				*/
/* file selection box ( OK, FILTER, CANCEL)	                        */
/* for saving images							*/
/************************************************************************/
{
char *filename;

    switch(type) {
	case OKAY:
	    XmStringGetLtoR(cbs->value, XmSTRING_DEFAULT_CHARSET, 
			 &filename);
            printf("file %s selected\n",filename);
            StorePGM(filename);
	    break;
	case CANCEL:
	    break;
    }
    XtUnmanageChild(FselBox);
}

// ==================================================================

void save_sequenceCB(Widget, int type, XmSelectionBoxCallbackStruct *cbs)
/************************************************************************/
/* CB initiated by save in file pulldown menu				*/
/* file selection box ( OK, FILTER, CANCEL)	                        */
/* for saving image sequences						*/
/************************************************************************/
{
char *filename;

    switch(type) {
	case OKAY:
	    XmStringGetLtoR(cbs->value, XmSTRING_DEFAULT_CHARSET, 
			 &filename);
            printf("file %s selected\n",filename);
            set_image_sequence_name(filename);
	    break;
	case CANCEL:
	    break;
    }
    XtUnmanageChild(FselBox);
}

// ==================================================================

void load_sequenceCB(Widget, int type, XmSelectionBoxCallbackStruct *cbs)
/************************************************************************/
/* CB initiated by save in file pulldown menu				*/
/* file selection box ( OK, FILTER, CANCEL)	                        */
/* for loading classified volume files					*/
/************************************************************************/
{
char *filename;

    switch(type) {
	case OKAY:
	    XmStringGetLtoR(cbs->value, XmSTRING_DEFAULT_CHARSET, 
			 &filename);
            printf("file %s selected\n",filename);
            XtUnmanageChild(FselBox);
            load_and_render_sequence(filename);
	    break;
	case CANCEL:
            XtUnmanageChild(FselBox);
	    break;
    }
    //XtUnmanageChild(FselBox);
}

// ==================================================================

void load_unclassified_sequenceCB(Widget, int type, XmSelectionBoxCallbackStruct *cbs)
/************************************************************************/
/* CB initiated by save in file pulldown menu				*/
/* file selection box ( OK, FILTER, CANCEL)	                        */
/* for loading unclassified volume files				*/
/************************************************************************/
{
char *filename;

    switch(type) {
	case OKAY:
	    XmStringGetLtoR(cbs->value, XmSTRING_DEFAULT_CHARSET, 
			 &filename);
            printf("file %s selected\n",filename);
            XtUnmanageChild(FselBox);
            load_unclassified_sequence(filename);
	    break;
	case CANCEL:
            XtUnmanageChild(FselBox);
	    break;
    }
    //XtUnmanageChild(FselBox);
}

// ==================================================================

void Create_FileBox(int identifier)
/************************************************************************/
/* Create a file selection box                                          */
/************************************************************************/
{
Widget mshell;
/****** create a dialog shell widget ******/
    mshell = XmCreateDialogShell(toplevel, "FileSelection", NULL, 0);


    //XtManageChild(mshell);
/****** create a file selection box ******/ 
    switch(identifier){
      case M_FILE_OPEN:  
        FselBox = XtVaCreateWidget("FileSel", 
		    xmFileSelectionBoxWidgetClass, mshell, 
                    XmNpattern,	XmStringCreateSimple("*.dat"),
                    NULL,0);
        XtManageChild(FselBox);
        XtUnmanageChild(XmFileSelectionBoxGetChild(FselBox,
		 XmDIALOG_HELP_BUTTON));
        XtAddCallback(FselBox,  XmNokCallback,	
		    (XtCallbackProc)open_fileCB, (XtPointer)OKAY);
        XtAddCallback(FselBox,  XmNcancelCallback,  
		    (XtCallbackProc)open_fileCB, (XtPointer)CANCEL);
    
        break;
      case M_FILE_OPEN_CLASSIFIED:  
        FselBox = XtVaCreateWidget("FileSel", 
		    xmFileSelectionBoxWidgetClass, mshell, 
                    XmNpattern,	XmStringCreateSimple("*.cv"),
                    NULL,0);
        XtManageChild(FselBox);
        XtUnmanageChild(XmFileSelectionBoxGetChild(FselBox,
		 XmDIALOG_HELP_BUTTON));
        XtAddCallback(FselBox,  XmNokCallback,	
		    (XtCallbackProc)openClassified_fileCB, (XtPointer)OKAY);
        XtAddCallback(FselBox,  XmNcancelCallback,  
		    (XtCallbackProc)openClassified_fileCB, (XtPointer)CANCEL);
    
        break;
      case M_FILE_SAVE_IMAGE:
        switch (CurrentColorMode){
          case LUMINANCE_MODE:
            FselBox = XtVaCreateWidget("FileSel", 
		    xmFileSelectionBoxWidgetClass, mshell, 
                    XmNpattern,	XmStringCreateSimple("*.pgm"),
                    NULL,0);
            break;
          case RGB_MODE:
            FselBox = XtVaCreateWidget("FileSel", 
		    xmFileSelectionBoxWidgetClass, mshell, 
                    XmNpattern,	XmStringCreateSimple("*.ppm"),
                    NULL,0);
            break;
          default: fprintf(stderr,"Error: Unknown CurrentColorMode");
                      exit(1);
       }
        XtManageChild(FselBox);
        XtUnmanageChild(XmFileSelectionBoxGetChild(FselBox,
		 XmDIALOG_HELP_BUTTON));
        XtAddCallback(FselBox,  XmNokCallback,	
		    (XtCallbackProc)save_imageCB, (XtPointer)OKAY);
        XtAddCallback(FselBox,  XmNcancelCallback,  
		    (XtCallbackProc)save_imageCB, (XtPointer)CANCEL);
        break;
      case M_FILE_SELECT_SEQUENCE:
        FselBox = XtVaCreateWidget("FileSel", 
		    xmFileSelectionBoxWidgetClass, mshell, 
                    XmNpattern,	XmStringCreateSimple("*.0"),
                    NULL,0);
        XtManageChild(FselBox);
        XtUnmanageChild(XmFileSelectionBoxGetChild(FselBox,
		 XmDIALOG_HELP_BUTTON));
        XtAddCallback(FselBox,  XmNokCallback,	
		    (XtCallbackProc)save_sequenceCB, (XtPointer)OKAY);
        XtAddCallback(FselBox,  XmNcancelCallback,  
		    (XtCallbackProc)save_sequenceCB, (XtPointer)CANCEL);
        
        break;
      case M_FILE_LOAD_SEQUENCE:
        FselBox = XtVaCreateWidget("FileSel", 
		    xmFileSelectionBoxWidgetClass, mshell, 
                    XmNpattern,	XmStringCreateSimple("*.cv.0"),
                    NULL,0);
        XtManageChild(FselBox);
        XtUnmanageChild(XmFileSelectionBoxGetChild(FselBox,
		 XmDIALOG_HELP_BUTTON));
        XtAddCallback(FselBox,  XmNokCallback,	
		    (XtCallbackProc)load_sequenceCB, (XtPointer)OKAY);
        XtAddCallback(FselBox,  XmNcancelCallback,  
		    (XtCallbackProc)load_sequenceCB, (XtPointer)CANCEL);
        
        break;
      case M_FILE_UNCLASSIFIED_SEQUENCE:
        FselBox = XtVaCreateWidget("FileSel", 
		    xmFileSelectionBoxWidgetClass, mshell, 
                    XmNpattern,	XmStringCreateSimple("*.dat.0"),
                    NULL,0);
        XtManageChild(FselBox);
        XtUnmanageChild(XmFileSelectionBoxGetChild(FselBox,
		 XmDIALOG_HELP_BUTTON));
        XtAddCallback(FselBox,  XmNokCallback,	
		    (XtCallbackProc)load_unclassified_sequenceCB, (XtPointer)OKAY);
        XtAddCallback(FselBox,  XmNcancelCallback,  
		    (XtCallbackProc)load_unclassified_sequenceCB, (XtPointer)CANCEL);
        
        break;
      default:
        fprintf(stderr,"Unmanaged case in Create_FileBox");
        break;
    }
}

// ==================================================================

void RadioCB(Widget myWidget, XtPointer client_data, XtPointer)
/************************************************************************/
/* toggles the selection of the group                                   */
/************************************************************************/
	
{
  RadioButtonGrTyp *radio = (RadioButtonGrTyp *) client_data;
  int i;
  for (i=0;i<radio->number;i++)
    XmToggleButtonSetState(radio->widget[i],(radio->widget[i]==myWidget),FALSE);
}

// ==================================================================

Widget createMenuBar(Widget parent, GlobalInfoTyp *global)
/************************************************************************/
/* creates the menu bar according to menu.inc       		        */
/************************************************************************/
{
    Arg			args[12];
    int			i, j, n, id;
    WidgetList		buttons, subButtons;
    int			itemCount, subItemCount;
    WidgetClass	    	widgetClass;
    String  	    	callbackReason;
    Widget 		topbarMenuWidget;
    Widget              radio[256];
    int                 numRadio;
    RadioButtonGrTyp   *RadioGroup=NULL;
    ButtonDataTyp      *ButtonData;

    // create topbar menu
    topbarMenuWidget = XmCreateMenuBar(parent, "menuBar", NULL, 0);

    itemCount = XtNumber(Menu);
    buttons   = (WidgetList) XtMalloc(itemCount * sizeof(Widget));
    ButtonData = (ButtonDataTyp*) XtMalloc(M_ITEMS * sizeof(ButtonDataTyp));

    for (i = 0; i < itemCount; i++) {
	// Make Topbar menu button
	Widget subMenu = XmCreatePulldownMenu(topbarMenuWidget, NULL, NULL, 0);
	
	// set button data for main menu buttons
	id = Menu[i].id;
	ButtonData[id].widget = subMenu;
	ButtonData[id].id = id;
	ButtonData[id].userData = 0.0;
	ButtonData[id].globalInfo=global;

        // add menu mapping callback
        if (Menu[i].callback)    
    	  XtAddCallback(subMenu, XmNmapCallback,Menu[i].callback,(XtPointer) &ButtonData[id]);

    	XtSetArg(args[0], XmNsubMenuId, subMenu);
	buttons[i] = XtCreateWidget(Menu[i].name,
	    xmCascadeButtonGadgetClass, topbarMenuWidget, args, 1);

	// make submenu buttons
	subItemCount = Menu[i].subItemCount;
	subButtons = (WidgetList) XtMalloc(subItemCount * sizeof(Widget));
	
	numRadio=0;
	
	for (j = 0; j < subItemCount; j++) {
	
	    if (Menu[i].subMenu[j].id == M_SEPARATOR)
	      {
		 subButtons[j] = XtCreateWidget(NULL, xmSeparatorGadgetClass, 
		      subMenu, NULL, 0);
		 if (numRadio)                  // Radio Groups end at separators
		   {
		     RadioGroup->number=numRadio;
		     RadioGroup->widget = (Widget *) XtMalloc(numRadio*sizeof(Widget)) ;
		     int k;
		     for (k=0;k<numRadio;k++)
		       RadioGroup->widget[k]=radio[k];
		     numRadio=0;
		   }
	      }
	    else {

		id = Menu[i].subMenu[j].id;
	        ButtonData[id].widget = subMenu;
		ButtonData[id].id = id;
		ButtonData[id].userData = Menu[i].subMenu[j].userData;
		ButtonData[id].globalInfo=global;

    	    	switch (Menu[i].subMenu[j].buttonType) {
		    case M_PUSH_BUTTON:
			widgetClass = xmPushButtonWidgetClass;
			callbackReason = XmNactivateCallback;
			n = 0;
		    	break;
		    case M_TOGGLE_BUTTON:
			widgetClass = xmToggleButtonWidgetClass;
			callbackReason = XmNvalueChangedCallback;
			n = 0;
		    	break;
		    case M_RADIO_BUTTON:
		        if (numRadio==0)        // If no group exists
			    RadioGroup = (RadioButtonGrTyp*)
			                         XtMalloc(sizeof(RadioButtonGrTyp));
			widgetClass = xmToggleButtonWidgetClass;
			callbackReason = XmNvalueChangedCallback;
			XtSetArg(args[0], XmNindicatorType, XmONE_OF_MANY);
			n = 1;
		    	break;
		    default:
			fprintf(stderr,"%d,%d,%s\n",__LINE__,__FILE__,"Menu INTERNAL ERROR: bad buttonType");
		    	break;
		}
		
		// check for keyboard accelerator
		char *accel = Menu[i].subMenu[j].accelerator;
		char *accelText = Menu[i].subMenu[j].accelText;
		if (accel != NULL) {
		    XtSetArg(args[n], XmNaccelerator, accel); n++;
		    
		    if (accelText != NULL) {
			XmString xmstr = XmStringCreate(accelText,
					 XmSTRING_DEFAULT_CHARSET);
			XtSetArg(args[n], XmNacceleratorText, xmstr); n++;
			//??? can we ever free the xmstr?
		    }
		}

	        XtSetArg(args[n], XmNlabelString, XmStringCreate(Menu[i].subMenu[j].name,
		                                  XmSTRING_DEFAULT_CHARSET));
		n++;
		
		subButtons[j] = XtCreateWidget(
		    Menu[i].subMenu[j].resname,
		    widgetClass,
		    subMenu, args, n);
		    
		// --- RADIOBUTTON special
		if (Menu[i].subMenu[j].buttonType == M_RADIO_BUTTON)
		  {
		    XtAddCallback(subButtons[j],callbackReason,RadioCB,(XtPointer) RadioGroup);
		    radio[numRadio] = subButtons[j];
		    numRadio++;
		  }
		  
		id = Menu[i].subMenu[j].id;
		ButtonData[id].widget = subButtons[j];
		ButtonData[id].id = id; 
		XtAddCallback(subButtons[j],callbackReason,Menu[i].subMenu[j].callback,
		              (XtPointer) &ButtonData[id]);
	    }

	}
        if (numRadio)                  // still an open Radio group ?
          {
            RadioGroup->number=numRadio;
            RadioGroup->widget = (Widget *) XtMalloc (numRadio*sizeof(Widget));
            int k;
            for (k=0;k<numRadio;k++)
               RadioGroup->widget[k]=radio[k];
            numRadio=0;
          }
	XtManageChildren(subButtons, subItemCount);
	XtFree((char *)subButtons);
    }
    XtManageChildren(buttons, itemCount);
    XtFree((char *)buttons);
    
    //
    // layout the menu bar
    //
    n = 0;
    XtSetArg(args[n], XmNtopAttachment, XmATTACH_FORM); n++;
    XtSetArg(args[n], XmNleftAttachment, XmATTACH_FORM); n++;
    XtSetArg(args[n], XmNrightAttachment, XmATTACH_FORM); n++;
    XtSetValues(topbarMenuWidget, args, n);
    
    XtManageChild(topbarMenuWidget);
    return topbarMenuWidget;
}



