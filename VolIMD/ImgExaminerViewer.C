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

//**************************************************************************
//
// * Description    : Inventor examiner viewer 
//                   
// * Class(es)      : ImgExaminerViewer
//
// * inherited from : ImgFullViewer ImgViewer
//
// * Author  : Dirk Rantzau / Andreas Werner
//
// * History : 29.03.94 V 1.0
//
//**************************************************************************

#include <X11/Intrinsic.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>

#include <Xm/LabelG.h>
#include <Xm/Text.h>
#include <Xm/ToggleB.h>
#include <Xm/ToggleBG.h>
#include <Xm/Form.h>
#include <Sgm/ThumbWheel.h>

#include <Inventor/SoDB.h>
#include <Inventor/SoInput.h>
#include <Inventor/sensors/SoTimerSensor.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/nodes/SoScale.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSwitch.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/projectors/SbSphereSheetProjector.h>
#include <Inventor/actions/SoSearchAction.h>
#include <Inventor/Xt/SoXt.h>
#include <Inventor/Xt/SoXtCursors.h>
#include <Inventor/Xt/SoXtResource.h>
#include <Inventor/Xt/SoXtIcons.h>
#include "ImgExaminerViewer.h"
#include "InvPixmapButton.h"
#include <GL/gl.h>


/*
 *  Defines
 */

enum ViewerModes {
    IDLE_MODE, 
    SPIN_MODE, 
    TRANS_MODE, 
    DOLLY_MODE, 
    ROLL_MODE, 
    ROLL_ACTIVE_MODE, 
    SEEK_MODE, 
};

// list of custom push buttons
enum {
    CAM_PUSH = 0,
    PUSH_NUM,
};


// size of the rotation buffer, which is used to animate the spinning ball.
#define ROT_BUFF_SIZE 3

#define ANIM_RATE (1/60.0)  // animation frame rate

/*
 * Macros
 */

#define KEY_PRESS_CHECK_VALID_MACRO \
    if (mode != IDLE_MODE && mode != SPIN_MODE && \
	mode != TRANS_MODE && mode != DOLLY_MODE) \
	break;

#define KEY_RELEASE_SWITCH_MODE_MACRO(ke) \
{ \
    int newMode = ((ke)->state & Button1Mask) ? \
	( ((ke)->state & Button2Mask) ? DOLLY_MODE : SPIN_MODE ) : \
	( ((ke)->state & Button2Mask) ? TRANS_MODE : IDLE_MODE ); \
    switchMode(newMode); \
}

//
// The point of interest geometry description
//
char *ImgExaminerViewer::geometryBuffer = "\
#Inventor V2.0 ascii\n\
\
Separator { \
    PickStyle { style UNPICKABLE } \
    LightModel { model BASE_COLOR } \
    MaterialBinding { value PER_PART } \
    DrawStyle { lineWidth 2 } \
    Coordinate3 { point [0 0 0, 1 0 0, 0 1 0, 0 0 1] } \
    BaseColor { rgb [1 0 0, 0 1 0, 0 0 1] } \
    IndexedLineSet { coordIndex [1, 0, 2, -1, 0, 3] } \
     \
    LightModel { model PHONG } \
    MaterialBinding { value OVERALL } \
    Complexity { value .1 } \
    Separator { \
    	Material { \
	    diffuseColor    [ 0.5 0 0 ] \
	    emissiveColor   [ 0.5 0 0 ] \
	} \
	Translation { translation 1 0 0 } \
    	RotationXYZ { axis Z angle -1.570796327 } \
    	Cone { bottomRadius .2 height .3 } \
    } \
    Separator { \
    	Material { \
	    diffuseColor    [ 0 0.5 0 ] \
	    emissiveColor   [ 0 0.5 0 ] \
	} \
	Translation { translation 0 1 0 } \
    	Cone { bottomRadius .2 height .3 } \
    } \
    Material { \
	diffuseColor    [ 0 0 0.5 ] \
	emissiveColor   [ 0 0 0.5 ] \
    } \
    Translation { translation 0 0 1 } \
    RotationXYZ { axis X angle 1.570796327 } \
    Cone { bottomRadius .2 height .3 } \
} ";


static char *thisClassName = "ImgExaminerViewer";

////////////////////////////////////////////////////////////////////////
//
// Public constructor - build the widget right now
//
ImgExaminerViewer::ImgExaminerViewer(
    Widget parent,
    const char *name, 
    SbBool buildInsideParent, 
    ImgFullViewer::BuildFlag b, 
    ImgViewer::Type t)
	: ImgFullViewer(
	    parent,
	    name, 
	    buildInsideParent, 
	    b, 
	    t, 
	    FALSE) // tell GLWidget not to build just yet  
//
////////////////////////////////////////////////////////////////////////
{
    // In this case, render area is what the app wants, so buildNow = TRUE
    constructorCommon(TRUE);
}

////////////////////////////////////////////////////////////////////////
//
// SoEXTENDER constructor - the subclass tells us whether to build or not
//
ImgExaminerViewer::ImgExaminerViewer(
    Widget parent,
    const char *name, 
    SbBool buildInsideParent, 
    ImgFullViewer::BuildFlag b, 
    ImgViewer::Type t, 
    SbBool buildNow)
	: ImgFullViewer(
	    parent,
	    name, 
	    buildInsideParent,
	    b,  
	    t, 
	    FALSE) // tell GLWidget not to build just yet  
//
////////////////////////////////////////////////////////////////////////
{
    // In this case, render area may be what the app wants, 
    // or it may want a subclass of render area. Pass along buildNow
    // as it was passed to us.
    constructorCommon(buildNow);
}

////////////////////////////////////////////////////////////////////////
//
// Called by the constructors
//
// private
//
void
ImgExaminerViewer::constructorCommon(SbBool buildNow)
//
////////////////////////////////////////////////////////////////////////
{
    // init local vars
    addVisibilityChangeCallback(visibilityChangeCB, this);
    mode = IDLE_MODE;
    createdCursors = FALSE;
    spinCursor = panCursor = dollyCursor = rollCursor = seekCursor = 0;
    firstBuild = TRUE;
    setSize( SbVec2s(500, 390) );  // default size
    setClassName(thisClassName);
    
    // feedback vars
    feedbackFlag = FALSE;
    feedbackRoot = NULL;
    feedbackSwitch = NULL;
    feedbackSize = 20.0;
    feedbackSizeWheel = NULL;
    
    // init animation variables
    animationEnabled = TRUE;
    animatingFlag = FALSE;
    rotBuffer = new SbRotation[ROT_BUFF_SIZE];
    animationSensor = new
	SoTimerSensor(ImgExaminerViewer::animationSensorCB, this);
    animationSensor->setInterval(ANIM_RATE);

    // init the projector class
    SbViewVolume vv;
    vv.ortho(-1, 1, -1, 1, -10, 10);
    sphereSheet = new SbSphereSheetProjector;
    sphereSheet->setViewVolume( vv );
    sphereSheet->setSphere( SbSphere( SbVec3f(0, 0, 0), .7) );
    
    // assign decoration names
    setPopupMenuString("Examiner Viewer");
    setBottomWheelString("Roty");
    setLeftWheelString("Rotx");
    setPrefSheetString("Examiner Viewer Preference Sheet");
    
    for (int i=0; i<PUSH_NUM; i++)
	buttonList[i] = NULL;
    

    
    // Build the widget tree, and let SoXtComponent know about our base widget.
    if (buildNow) {
	Widget w = buildWidget(getParentWidget());
	setBaseWidget(w);
    }
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Destructor.
//
// Use: public

ImgExaminerViewer::~ImgExaminerViewer()
//
////////////////////////////////////////////////////////////////////////
{
    delete animationSensor;
    
    for (int i=0; i<PUSH_NUM; i++)
	delete buttonList[i];
    
    delete sphereSheet;
    if (feedbackRoot)
      feedbackRoot->unref();
 
    // free the viewer cursors
    if (getDisplay()) {
	Display *display = getDisplay();
	if (spinCursor) XFreeCursor(display, spinCursor);
	if (panCursor) XFreeCursor(display, panCursor);
	if (dollyCursor) XFreeCursor(display, dollyCursor);
	if (rollCursor) XFreeCursor(display, rollCursor);
	if (seekCursor) XFreeCursor(display, seekCursor);
    }
 
    
    delete [] rotBuffer;
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Shows/hide the point of rotation feedback geometry
//
// Use: public
void
ImgExaminerViewer::setFeedbackVisibility(SbBool flag)
//
////////////////////////////////////////////////////////////////////////
{
    if (flag == feedbackFlag || camera == NULL) {
	feedbackFlag = flag;
	return;
    }
    
    //
    // find the camera parent to insert/remove the feedback root
    //
    SoSearchAction sa;
    sa.setNode(flag ? (SoNode *)camera : (SoNode *)feedbackRoot);
    sa.apply(sceneRoot);
    SoFullPath *path = (SoFullPath *) sa.getPath();
    if (!path) {

	return;
    }
    SoGroup *parent = (SoGroup *)path->getNode(path->getLength() - 2);
    
    feedbackFlag = flag;
    
    // make sure the feedback has been built
    if (!feedbackRoot)
	createFeedbackNodes();
    
    //
    // inserts/remove the feedback axis group
    //
    if (feedbackFlag) {
	
	// return if geometry is already there (this should be an error !)
	if (parent->findChild(feedbackRoot) >= 0)
	    return;
	
	// place the feedback right after the headlight if the headlight
	// is turned on. Otherwise place it right after the camera.
	if ( isHeadlight() )
	    parent->insertChild(feedbackRoot, parent->findChild(camera) + 2);
	else
	    parent->insertChild(feedbackRoot, parent->findChild(camera) + 1);
	
	// make sure the feedback switch is turned to the correct state now
	// that the feedback root has been inserted in the scene
	feedbackSwitch->whichChild.setValue(viewingFlag ? SO_SWITCH_ALL : SO_SWITCH_NONE);
    }
    else
	parent->removeChild(feedbackRoot);
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Sets the feedback size.
//
// Use: public
void
ImgExaminerViewer::setFeedbackSize(int newSize)
//
////////////////////////////////////////////////////////////////////////
{
    if (feedbackSize == newSize)
	return;
    
    // assign new value and redraw (since it is not a field in the scene)
    feedbackSize = newSize;
    if (isFeedbackVisible() && isViewing())
	redraw();
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Remove the extra geometry when doing a viewAll.
//
// Use: virtual public
void
ImgExaminerViewer::viewAll()
//
////////////////////////////////////////////////////////////////////////
{
    // stop spinning
    if ( isAnimating() )
    	stopAnimating();
    
    // temporarily remove the feedback geometry
    if (feedbackFlag && isViewing() && feedbackSwitch)
	feedbackSwitch->whichChild.setValue( SO_SWITCH_NONE );
    
    // call the base class
    ImgFullViewer::viewAll();
    
    // now add the geometry back in
    if (feedbackFlag && isViewing() && feedbackSwitch)
	feedbackSwitch->whichChild.setValue( SO_SWITCH_ALL );
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Call the parent class and stop animation if any
//
// Use: virtual public
void
ImgExaminerViewer::resetToHomePosition()
//
////////////////////////////////////////////////////////////////////////
{
    // stop spinning
    if ( isAnimating() )
    	stopAnimating();
    
    // call the base class
    ImgFullViewer::resetToHomePosition();
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Call the parent class and insert/remove the feedback root
//
// Use: virtual public
void
ImgExaminerViewer::setCamera(SoCamera *newCamera)
//
////////////////////////////////////////////////////////////////////////
{
    if (camera == newCamera)
	return;
    
    // set the right thumbwheel label and toggle button image based on 
    // the camera type
    if (newCamera != NULL && (camera == NULL || 
	newCamera->getTypeId() != camera->getTypeId())) {
	if (newCamera->isOfType(SoOrthographicCamera::getClassTypeId())) {
	    if (buttonList[CAM_PUSH])
		buttonList[CAM_PUSH]->setIcon(so_xt_ortho_bits, 
		    so_xt_icon_width, so_xt_icon_height);
	    setRightWheelString("Zoom");
	}
	else {
	    if (buttonList[CAM_PUSH])
		buttonList[CAM_PUSH]->setIcon(so_xt_persp_bits, 
		    so_xt_icon_width, so_xt_icon_height);
	    setRightWheelString("Dolly");
	}
    }
    
    // detach feedback which depends on camera
    if ( feedbackFlag ) {
	setFeedbackVisibility(FALSE);
	feedbackFlag = TRUE;  // can later be turned on
    }
    
    // call parent class
    ImgFullViewer::setCamera(newCamera);
    
    // attach feedback back on
    if ( feedbackFlag ) {
	feedbackFlag = FALSE; // enables routine to be called
	setFeedbackVisibility(TRUE);
    }
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Call the base class and sets the correct cursors on the window,
//  plus insert/remove the feedback geometry.
//
// Use: virtual public
void
ImgExaminerViewer::setViewing(SbBool flag)
//
////////////////////////////////////////////////////////////////////////
{
    if (flag == viewingFlag)
	return;
    
    // call the parent class
    ImgFullViewer::setViewing(flag);
    
    // set the right cursor on the window
    Widget w = getRenderAreaWidget();
    Window window = (w != NULL) ? XtWindow(w) : NULL;
    if (window != NULL) {
	if (!createdCursors)
	    defineCursors();
	if (isViewing())
	    XDefineCursor(XtDisplay(w), window, spinCursor);
	else
	    XUndefineCursor(XtDisplay(w), window);
    }
    
    // stops animation if any
    if ( isAnimating() )
	stopAnimating();
    
    // show/hide the feedback geometry based on the viewing state
    if (feedbackFlag && feedbackSwitch)
	feedbackSwitch->whichChild.setValue(viewingFlag ? SO_SWITCH_ALL : SO_SWITCH_NONE);

}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Process the given event to do viewing stuff
//
// Use: virtual protected
void
ImgExaminerViewer::processEvent(XAnyEvent *xe)
//
////////////////////////////////////////////////////////////////////////
{
    if ( processCommonEvents(xe) )
	return;
    
    if (!createdCursors) {
	defineCursors();
	Widget w = getRenderAreaWidget();
	XDefineCursor(XtDisplay(w), XtWindow(w), spinCursor);
    }
    
    XButtonEvent    *be;
    XMotionEvent    *me;
    XKeyEvent	    *ke;
    
    SbVec2s raSize = getGlxSize();
    
    switch(xe->type) {
    case ButtonPress:
	be = (XButtonEvent *)xe;
	locator[0] = be->x;
	locator[1] = raSize[1] - be->y;
	if (be->button == Button1) {
	    
	    if (mode != SEEK_MODE)
		interactiveCountInc();
	    stopAnimating();
	    
	    switch (mode) {
		case IDLE_MODE:	switchMode(SPIN_MODE); break;
		case TRANS_MODE: switchMode(DOLLY_MODE); break;
		case ROLL_MODE:	switchMode(ROLL_ACTIVE_MODE); break;
		case SEEK_MODE: seekToPoint(locator); break;
	    }
	}
	else if (be->button == Button2) {
	    
	    if (mode != SEEK_MODE)
		interactiveCountInc();
	    stopAnimating();
	    
	    switch (mode) {
		case IDLE_MODE: switchMode(TRANS_MODE); break;
		case SPIN_MODE: switchMode(DOLLY_MODE); break;
	    }
	}
	break;
	
    case ButtonRelease:
	be = (XButtonEvent *)xe;
	locator[0] = be->x;
	locator[1] = raSize[1] - be->y;
	if (be->button == Button1) {
	    switch (mode) {
		case SPIN_MODE:
		    // check if we need to start spinning
		    if (animationEnabled && lastMotionTime == be->time) {
			animatingFlag = TRUE;
			computeAverage = TRUE;
			animationSensor->schedule();
			interactiveCountInc();
		    }
		    switchMode(IDLE_MODE);
		    break;
		case DOLLY_MODE: switchMode(TRANS_MODE); break;
		case ROLL_ACTIVE_MODE: switchMode(ROLL_MODE); break;
	    }
	    
	    if (mode != SEEK_MODE)
		interactiveCountDec();
	}
	else if (be->button == Button2) {
	    switch (mode) {
		case TRANS_MODE: switchMode(IDLE_MODE); break;
		case DOLLY_MODE: switchMode(SPIN_MODE); break;
	    }
	    
	    if (mode != SEEK_MODE)
		interactiveCountDec();
	}
	break;
	
    case MotionNotify:
	me = (XMotionEvent *)xe;
	switch (mode) {
	    case SPIN_MODE:
		lastMotionTime = me->time;
		spinCamera(SbVec2f(me->x/float(raSize[0]), (raSize[1] - me->y)/float(raSize[1])));
		break;
	    case TRANS_MODE:
		translateCamera(SbVec2f(me->x/float(raSize[0]), (raSize[1] - me->y)/float(raSize[1])));
		break;
	    case DOLLY_MODE:
		dollyCamera( SbVec2s(me->x, raSize[1] - me->y) );
		break;
	    case ROLL_ACTIVE_MODE:
		rollCamera( SbVec2s(me->x, raSize[1] - me->y) );
		break;
	}
	break;
	
    case KeyPress:
	ke = (XKeyEvent *)xe;
	locator[0] = ke->x;
	locator[1] = raSize[1] - ke->y;
	switch ( XLookupKeysym(ke, 0) ) {
	    case XK_Control_L:
	    case XK_Control_R:
		KEY_PRESS_CHECK_VALID_MACRO
		if (ke->state & Button1Mask)
		    switchMode(ROLL_ACTIVE_MODE);
		else
		    switchMode(ROLL_MODE);
		break;
	}
	break;
	
    case KeyRelease:
	ke = (XKeyEvent *)xe;
	locator[0] = ke->x;
	locator[1] = raSize[1] - ke->y;
	switch ( XLookupKeysym(ke, 0) ) {
	    case XK_Control_L:
	    case XK_Control_R:
		if (mode != ROLL_MODE && mode != ROLL_ACTIVE_MODE)
		    break;
		KEY_RELEASE_SWITCH_MODE_MACRO(ke)
		break;
	}
	break;
    }
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    switches to the specified viewer mode
//
// Use: private
void
ImgExaminerViewer::switchMode(int newMode)
//
////////////////////////////////////////////////////////////////////////
{
    // assing new mode
    SbBool redrawNeeded = FALSE;
    int prevMode = mode;
    mode = newMode;
    
    // needed to define new cursors
    Widget w = getRenderAreaWidget();
    Display *display = XtDisplay(w);
    Window window = XtWindow(w);
    SbVec2s raSize = getGlxSize();
    
    // check the old viewer mode for redraw need
    switch (prevMode) {
	case ROLL_MODE:
	case ROLL_ACTIVE_MODE:
	    redrawNeeded = TRUE;
	    break;
    }
    
    // switch to new viewer mode
    switch (newMode) {
	case IDLE_MODE:
	    XDefineCursor(display, window, spinCursor);
	    break;
	    
	case SPIN_MODE:
	    XDefineCursor(display, window, spinCursor);
	    
	    // set the sphere sheet starting point
	    sphereSheet->project(
		SbVec2f(locator[0]/float(raSize[0]), locator[1]/float(raSize[1])) );
	    
	    // reset the animation queue
	    firstIndex = 0;
	    lastIndex = -1;
	    break;
	    
	case TRANS_MODE:
	    XDefineCursor(display, window, panCursor);
	    
	    {
	    // Figure out the focal plane
	    SbMatrix mx;
	    mx = camera->orientation.getValue();
	    SbVec3f forward(-mx[2][0], -mx[2][1], -mx[2][2]);
	    SbVec3f fp = camera->position.getValue() + 
		forward * camera->focalDistance.getValue();
	    focalplane = SbPlane(forward, fp);
	    
	    // map mouse starting position onto the panning plane
	    SbViewVolume    cameraVolume;
	    SbLine	    line;
	    cameraVolume = camera->getViewVolume();
	    cameraVolume.projectPointToLine(
		SbVec2f(locator[0]/float(raSize[0]), locator[1]/float(raSize[1])), line);
	    focalplane.intersect(line, locator3D);
	    }
	    break;
	    
	case DOLLY_MODE:
	    XDefineCursor(display, window, dollyCursor);
	    break;
	    
	case ROLL_MODE:
	    XDefineCursor(display, window, rollCursor);
	    redrawNeeded = TRUE;
	    break;
	    
	case ROLL_ACTIVE_MODE:
	    XUndefineCursor(display, window);
	    redrawNeeded = TRUE;
	    break;
    }
    
    if (redrawNeeded)
	redraw();
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    draws viewer feedback during a render area redraw of the scene.
//
// Use: virtual public
void
ImgExaminerViewer::actualRedraw()
//
////////////////////////////////////////////////////////////////////////
{
    // place the feedback at the focal point
    // ??? we really only need to do this when the camera changes
    if (isViewing() && feedbackFlag && camera != NULL && feedbackRoot) {
	
	// adjust the position to be at the focal point
	SbMatrix mx;
	mx = camera->orientation.getValue();
	SbVec3f forward(-mx[2][0], -mx[2][1], -mx[2][2]);
	feedbackTransNode->translation = camera->position.getValue() + 
	    camera->focalDistance.getValue() * forward;
	
	// adjust the size to be constant on the screen
	float height;
	if (camera->isOfType(SoPerspectiveCamera::getClassTypeId())) {
	    float angle = ((SoPerspectiveCamera *)camera)->heightAngle.getValue();
	    height = camera->focalDistance.getValue() * ftan(angle/2);
	}
	else if (camera->isOfType(SoOrthographicCamera::getClassTypeId()))
	    height = ((SoOrthographicCamera *)camera)->height.getValue() / 2;
	
	// ??? getGlxSize[1] == 0 the very first time, so return in that case
	// ??? else the redraws are 3 times slower from now on !! (alain)
	if (getGlxSize()[1] != 0) {
	    float size = 2.0 * height * feedbackSize / float (getGlxSize()[1]);
	    feedbackScaleNode->scaleFactor.setValue(size, size, size);
	}
    }
    
    // have the base class draw the scene
    ImgFullViewer::actualRedraw();
    
    // now draw the viewer extra feedback
    if (isViewing() && (mode == ROLL_MODE || mode == ROLL_ACTIVE_MODE)) {
	
	setFeedbackOrthoProjection(getGlxSize());
	
	if (mode == ROLL_ACTIVE_MODE)
	    drawViewerRollFeedback(getGlxSize()/2, locator);
	else
	    drawViewerCrossFeedback(getGlxSize()/2);
	
	// now restore state
	restoreGLStateAfterFeedback();
    }
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Enable/disable the animation feature of the Examiner
//
// Use: public
void
ImgExaminerViewer::setAnimationEnabled(SbBool flag)
//
////////////////////////////////////////////////////////////////////////
{
    if (animationEnabled == flag)
	return;
    
    animationEnabled = flag;
    if ( !animationEnabled && isAnimating())
        stopAnimating();
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Stops ongoing animation (if any)
//
// Use: public
void
ImgExaminerViewer::stopAnimating()
//
////////////////////////////////////////////////////////////////////////
{
    if (animatingFlag) {
	animatingFlag = FALSE;
	animationSensor->unschedule();
	interactiveCountDec();
    }
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Call the base class and sets the correct cursors on the window
//
// Use: virtual protected
void
ImgExaminerViewer::setSeekMode(SbBool flag)
//
////////////////////////////////////////////////////////////////////////
{
    if ( !isViewing() )
	return;
    
    // stop spinning
    if (isAnimating())
    	stopAnimating();
    
    // call the base class
    ImgFullViewer::setSeekMode(flag);
    
    mode = isSeekMode() ? SEEK_MODE : IDLE_MODE;
    
    // set the right cursor now
    Widget w = getRenderAreaWidget();
    if (w != NULL && XtWindow(w) != NULL) {
	if (!createdCursors)
	    defineCursors();
	if (isSeekMode())
	    XDefineCursor(XtDisplay(w), XtWindow(w), seekCursor);
	else
	    XDefineCursor(XtDisplay(w), XtWindow(w), spinCursor);
    }
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Redefine this routine to add some viewer specific stuff.
//
// Use: virtual protected
void
ImgExaminerViewer::createPrefSheet()
//
////////////////////////////////////////////////////////////////////////
{
    // create the preference sheet shell and form widget
    Widget shell, form;
    createPrefSheetShellAndForm(shell, form);
    
    // create all of the parts
    Widget widgetList[20];
    int num = 0;
    //RN widgetList[num++] = createSeekPrefSheetGuts(form);
    widgetList[num++] = createZoomPrefSheetGuts(form);
    widgetList[num++] = ImgFullViewer::createClippingPrefSheetGuts(form);
    widgetList[num++] = createStereoPrefSheetGuts(form);
    widgetList[num++] = createExamPrefSheetGuts(form);
    
    layoutPartsAndMapPrefSheet(widgetList, num, form, shell);
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Creates the viewer extra pref sheet stuff
//
// Use: private
Widget
ImgExaminerViewer::createExamPrefSheetGuts(Widget parent)
//
////////////////////////////////////////////////////////////////////////
{
    // 
    // changed 11.04.94 due to bus error in the COVISE renderer
    // when adding the point of rotation axis
    // D. Rantzau
    // ( commented out by /// )

    // this change was removed because customers
    // wanted the possibility of a display of the rotation axis 
    // R. Niemeier 03.04.97
    
    Widget toggles[2], labels[2];
///    Widget toggles[1], labels[1];
    Arg args[12];
    int n;
    
    // create a form to hold everything together
    Widget form = XtCreateWidget("", xmFormWidgetClass, parent, NULL, 0);
    
    // create all the parts
    n = 0;
    XtSetArg(args[n], XmNset, animationEnabled); n++;
    XtSetArg(args[n], XmNspacing, 0); n++;
    XtSetArg(args[n], XmNhighlightThickness, 0); n++;
    toggles[0] = XtCreateWidget("", xmToggleButtonGadgetClass, form, args, n);
    labels[0] = XtCreateWidget("Enable spin animation", 
	xmLabelGadgetClass, form, NULL, 0);
    XtAddCallback(toggles[0], XmNvalueChangedCallback, 
	(XtCallbackProc) ImgExaminerViewer::animPrefSheetToggleCB, 
	(XtPointer) this);
    
    n = 0;
    XtSetArg(args[n], XmNset, feedbackFlag); n++;
    XtSetArg(args[n], XmNspacing, 0); n++;
    XtSetArg(args[n], XmNhighlightThickness, 0); n++;
    toggles[1] = XtCreateWidget("", xmToggleButtonGadgetClass, form, args, n);
    labels[1] = XtCreateWidget("Show point of rotation axes", 
	xmLabelGadgetClass, form, NULL, 0);
    XtAddCallback(toggles[1], XmNvalueChangedCallback, 
	(XtCallbackProc) ImgExaminerViewer::feedbackPrefSheetToggleCB, 
	(XtPointer) this);
    
    // layout
    n = 0;
    XtSetArg(args[n], XmNleftAttachment,    XmATTACH_WIDGET); n++;
    XtSetArg(args[n], XmNleftWidget,	    toggles[0]); n++;
    XtSetArg(args[n], XmNtopAttachment,	    XmATTACH_OPPOSITE_WIDGET); n++;
    XtSetArg(args[n], XmNtopWidget,	    toggles[0]); n++;
    XtSetArg(args[n], XmNbottomAttachment,  XmATTACH_OPPOSITE_WIDGET); n++;
    XtSetArg(args[n], XmNbottomWidget,	    toggles[0]); n++;
    XtSetValues(labels[0], args, n);
    
    n = 0;
    XtSetArg(args[n], XmNtopAttachment,	    XmATTACH_WIDGET); n++;
    XtSetArg(args[n], XmNtopWidget,	    toggles[0]); n++;
    XtSetArg(args[n], XmNtopOffset,	    10); n++;
    XtSetValues(toggles[1], args, n);
    
    n = 0;
    XtSetArg(args[n], XmNleftAttachment,    XmATTACH_WIDGET); n++;
    XtSetArg(args[n], XmNleftWidget,	    toggles[1]); n++;
    XtSetArg(args[n], XmNtopAttachment,	    XmATTACH_OPPOSITE_WIDGET); n++;
    XtSetArg(args[n], XmNtopWidget,	    toggles[1]); n++;
    XtSetArg(args[n], XmNbottomAttachment,  XmATTACH_OPPOSITE_WIDGET); n++;
    XtSetArg(args[n], XmNbottomWidget,	    toggles[1]); n++;
    XtSetValues(labels[1], args, n);
    
    // manage children
    XtManageChildren(toggles, 2);
    XtManageChildren(labels, 2);
///     XtManageChildren(toggles, 1);
///     XtManageChildren(labels, 1);
    
    if (feedbackFlag && camera)
	toggleFeedbackWheelSize(toggles[1]);
    
    return form;
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Brings the viewer help card (called by "?" push button)
//
// Use: virtual protected
void
ImgExaminerViewer::openViewerHelpCard()
//
////////////////////////////////////////////////////////////////////////
{
    // tell the base class to open the file for us
    openHelpCard("MyExaminerViewer.help");
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Call the base class and stop animating
//
// Use: virtual protected

void
ImgExaminerViewer::bottomWheelStart()
//
////////////////////////////////////////////////////////////////////////
{
    ImgFullViewer::bottomWheelStart();
    stopAnimating();
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Call the base class and stop animating
//
// Use: virtual protected

void
ImgExaminerViewer::leftWheelStart()
//
////////////////////////////////////////////////////////////////////////
{
    ImgFullViewer::bottomWheelStart();
    stopAnimating();
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Rotates the object around the screen x axis (called by thumb wheel).
//
// Use: virtual protected

void
ImgExaminerViewer::bottomWheelMotion(float newVal)
//
////////////////////////////////////////////////////////////////////////
{
    // get rotation and apply to camera
    SbVec3f axis(0, 1, 0);
    SbRotation rot(axis, bottomWheelVal - newVal);
    rotateCamera(rot);
    
    bottomWheelVal = newVal;
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Rotates the object around the screen y axis (called by thumb wheel).
//
// Use: virtual protected

void
ImgExaminerViewer::leftWheelMotion(float newVal)
//
////////////////////////////////////////////////////////////////////////
{
    // get rotation and apply to camera
    SbVec3f axis(1, 0, 0);
    SbRotation rot(axis, newVal - leftWheelVal);
    rotateCamera(rot);
    
    leftWheelVal = newVal;
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//   	Moves the camera closer/further away from the plane of interest
//  (perspective camera case), else change the camera height (orthographic
//  camera case).
//
// Use: virtual protected

void
ImgExaminerViewer::rightWheelMotion(float newVal)
//
////////////////////////////////////////////////////////////////////////
{
    if (camera == NULL)
	return;
    
    if (camera->isOfType(SoOrthographicCamera::getClassTypeId())) {
	// change the ortho camera height
	SoOrthographicCamera *cam = (SoOrthographicCamera *) camera;
	cam->height = cam->height.getValue() * powf(2.0, newVal - rightWheelVal);
    }
    else {
	// shorter/grow the focal distance given the wheel rotation
	float focalDistance = camera->focalDistance.getValue();;
	float newFocalDist = focalDistance;
	newFocalDist *= powf(2.0, newVal - rightWheelVal);
	
	// finally reposition the camera
	SbMatrix mx;
	mx = camera->orientation.getValue();
	SbVec3f forward(-mx[2][0], -mx[2][1], -mx[2][2]);
	camera->position = camera->position.getValue() + 
			   (focalDistance - newFocalDist) * forward;
	camera->focalDistance = newFocalDist;
    }
    
    rightWheelVal = newVal;
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    This routine is used to define cursors (can only be called after
//  window has been realized).
//
// Use: private

void
ImgExaminerViewer::defineCursors()
//
////////////////////////////////////////////////////////////////////////
{
    XColor foreground;
    Pixmap source;
    Display *display = getDisplay();
    Drawable d = DefaultRootWindow(display);
    
    // set color
    foreground.red = 65535;
    foreground.green = foreground.blue = 0;
    
    // spin cursor
    source = XCreateBitmapFromData(display, d, 
	so_xt_curved_hand_bits, so_xt_curved_hand_width, so_xt_curved_hand_height);
    spinCursor = XCreatePixmapCursor(display, source, source, 
	&foreground, &foreground, so_xt_curved_hand_x_hot, so_xt_curved_hand_y_hot);
    XFreePixmap(display, source);
    
    // panning cursor
    source = XCreateBitmapFromData(display, d, 
	so_xt_flat_hand_bits, so_xt_flat_hand_width, so_xt_flat_hand_height);
    panCursor = XCreatePixmapCursor(display, source, source, 
	&foreground, &foreground, so_xt_flat_hand_x_hot, so_xt_flat_hand_y_hot);
    XFreePixmap(display, source);
    
    // dolly cursor
    source = XCreateBitmapFromData(display, d, 
	so_xt_pointing_hand_bits, so_xt_pointing_hand_width, so_xt_pointing_hand_height);
    dollyCursor = XCreatePixmapCursor(display, source, source, 
	&foreground, &foreground, so_xt_pointing_hand_x_hot, so_xt_pointing_hand_y_hot);
    XFreePixmap(display, source);
    
    // rolling cursor
    source = XCreateBitmapFromData(display, d, 
	so_xt_roll_bits, so_xt_roll_width, so_xt_roll_height);
    rollCursor = XCreatePixmapCursor(display, source, source, 
	&foreground, &foreground, so_xt_roll_x_hot, so_xt_roll_y_hot);
    XFreePixmap(display, source);
    
    // seek cursor
    source = XCreateBitmapFromData(display, d, 
	so_xt_target_bits, so_xt_target_width, so_xt_target_height);
    seekCursor = XCreatePixmapCursor(display, source, source, 
	&foreground, &foreground, so_xt_target_x_hot, so_xt_target_y_hot);
    XFreePixmap(display, source);
    
    createdCursors = TRUE;
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Rotates the camera around pt of interest by given rotation
//
// Use: private

void
ImgExaminerViewer::rotateCamera(const SbRotation &rot)
//
////////////////////////////////////////////////////////////////////////
{
    if (camera == NULL)
	return;
    
    // get center of rotation
    SbRotation camRot = camera->orientation.getValue();
    float radius = camera->focalDistance.getValue();
    SbMatrix mx;
    mx = camRot;
    SbVec3f forward( -mx[2][0], -mx[2][1], -mx[2][2]);
    SbVec3f center = camera->position.getValue()
	+ radius * forward;
    
    // apply new rotation to the camera
    camRot = rot * camRot;
    camera->orientation = camRot;
    
    // reposition camera to look at pt of interest
    mx = camRot;
    forward.setValue( -mx[2][0], -mx[2][1], -mx[2][2]);
    camera->position = center - radius * forward;
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Rolls the camera around it's forward direction given the new mouse
//  location.
//
// Use: private

void
ImgExaminerViewer::rollCamera(const SbVec2s &newLocator)
//
////////////////////////////////////////////////////////////////////////
{
    if (camera == NULL)
	return;
    
    SbVec2s center = getGlxSize()/2;
    SbVec2s p1, p2;
    float angle;
    
    // get angle of rotation
    p1 = locator - center;
    p2 = newLocator - center;
    // checking needed so that NaN won't occur
    angle = (p2[0]==0 && p2[1]==0) ? 0 : atan2(p2[1], p2[0]);
    angle -= (p1[0]==0 && p1[1]==0) ? 0 : atan2(p1[1], p1[0]);
    
    // now find the rotation and rotate camera
    SbVec3f axis(0, 0, -1);
    SbRotation rot;
    rot.setValue(axis, angle);
    camera->orientation = rot * camera->orientation.getValue();
    
    locator = newLocator;
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Moves the camera into the plane defined by the camera forward vector
//  and the focal point to follow the new mouse location.
//
// Use: private

void
ImgExaminerViewer::translateCamera(const SbVec2f &newLocator)
//
////////////////////////////////////////////////////////////////////////
{
    if (camera == NULL)
	return;
    
    // map new mouse location into the camera focal plane
    SbViewVolume    cameraVolume;
    SbLine	    line;
    SbVec3f	    newLocator3D;
    cameraVolume = camera->getViewVolume();
    cameraVolume.projectPointToLine(newLocator, line);
    focalplane.intersect(line, newLocator3D);
    
    // move the camera by the delta 3D position amount
    camera->position = camera->position.getValue() + 
	(locator3D - newLocator3D);
    
    // You would think we would have to set locator3D to
    // newLocator3D here.  But we don't, because moving the camera
    // essentially makes locator3D equal to newLocator3D in the
    // transformed space, and we will project the next newLocator3D in
    // this transformed space.
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    rotates the camera using the sheet sphere projector, given the new
//  mouse location.
//
// Use: private

void
ImgExaminerViewer::spinCamera(const SbVec2f &newLocator)
//
////////////////////////////////////////////////////////////////////////
{
    // find rotation and rotate camera
    SbRotation rot;
    sphereSheet->projectAndGetRotation(newLocator, rot);
    rot.invert();

    rotateCamera(rot);
    
    // save rotation for animation
    lastIndex = ((lastIndex+1) % ROT_BUFF_SIZE);
    rotBuffer[lastIndex] = rot;
    
    // check if queue is full
    if (((lastIndex+1) % ROT_BUFF_SIZE) == firstIndex)
	firstIndex = ((firstIndex+1) % ROT_BUFF_SIZE);
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Moves the camera forward/backward based on the new mouse potion.
//  (perspective camera), else change the camera height (orthographic
//  camera case).
//
// Use: private

void
ImgExaminerViewer::dollyCamera(const SbVec2s &newLocator)
//
////////////////////////////////////////////////////////////////////////
{
    if (camera == NULL)
	return;
    
    // moving the mouse up/down will move the camera futher/closer.
    // moving the camera sideway will not move the camera at all
    float d = (newLocator[1] - locator[1]) / 40.0;
    
    if (camera->isOfType(SoOrthographicCamera::getClassTypeId())) {
	// change the ortho camera height
	SoOrthographicCamera *cam = (SoOrthographicCamera *) camera;
	cam->height = cam->height.getValue() * powf(2.0, d);
    }
    else {
	// shorter/grow the focal distance given the mouse move
	float focalDistance = camera->focalDistance.getValue();;
	float newFocalDist = focalDistance * powf(2.0, d);
	
	// finally reposition the camera
	SbMatrix mx;
	mx = camera->orientation.getValue();
	SbVec3f forward(-mx[2][0], -mx[2][1], -mx[2][2]);
	camera->position = camera->position.getValue() + 
			   (focalDistance - newFocalDist) * forward;
	camera->focalDistance = newFocalDist;
    }
    
    locator = newLocator;
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Routine which animates the ball spinning (called by sensor).
//
// Use: private

void
ImgExaminerViewer::doSpinAnimation()
//
////////////////////////////////////////////////////////////////////////
{
    //
    // check if average rotation needs to be computed
    //
    
    if (computeAverage) {
	float averageAngle, angle;
	SbVec3f averageAxis, axis;
	
	// get number of samples
	int num = (((lastIndex - firstIndex) + 1 + 
	    ROT_BUFF_SIZE) % ROT_BUFF_SIZE);
	
	// check for not enough samples
	if (num < 2) {
	    stopAnimating();
	    return;
	}
	
	// get average axis of rotation
	// ??? right now only take one sample
	rotBuffer[firstIndex].getValue(averageAxis, angle);
	
	// get average angle of rotation
	averageAngle = 0;
	for (int i=0; i<num; i++) {
	    int n = (firstIndex + i) % ROT_BUFF_SIZE;
	    rotBuffer[n].getValue(axis, angle);
	    averageAngle += angle;
	}
	averageAngle /= float(num);
	
	averageRotation.setValue(averageAxis, averageAngle);
	computeAverage = FALSE;
    }
    
    //
    // rotate camera by average rotation
    //
    rotateCamera(averageRotation);
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Show/hide the pref sheet feedback size wheel and label.
//
// Use: private

void
ImgExaminerViewer::toggleFeedbackWheelSize(Widget toggle)
//
////////////////////////////////////////////////////////////////////////
{
    if ( feedbackFlag ) {
	Widget parent = XtParent(toggle);
	Arg args[12];
	int n;
	
	// create the label/thumb/text/label in the toggle parent
	feedbackLabel[0] = XtCreateWidget("axes size:", 
	    xmLabelGadgetClass, parent, NULL, 0);
	
	n = 0;
	XtSetArg(args[n], XmNvalue, 0); n++;
	XtSetArg(args[n], SgNangleRange, 0); n++;
	XtSetArg(args[n], SgNunitsPerRotation, 360); n++;
	XtSetArg(args[n], SgNshowHomeButton, FALSE); n++;
	XtSetArg(args[n], XmNhighlightThickness, 0); n++;
	XtSetArg(args[n], XmNorientation, XmHORIZONTAL); n++;
	feedbackSizeWheel = SgCreateThumbWheel(parent, NULL, args, n);
	
	XtAddCallback(feedbackSizeWheel, XmNvalueChangedCallback, 
	    (XtCallbackProc) ImgExaminerViewer::feedbackSizeWheelCB, (XtPointer) this);
	XtAddCallback(feedbackSizeWheel, XmNdragCallback, 
	    (XtCallbackProc) ImgExaminerViewer::feedbackSizeWheelCB, (XtPointer) this);
	feedbackSizeWheelVal = 0;
	
	n = 0;
	char str[15];
	sprintf(str, "%d", int(feedbackSize));
	XtSetArg(args[n], XmNvalue, str); n++;
	XtSetArg(args[n], XmNhighlightThickness, 1); n++;
	XtSetArg(args[n], XmNcolumns, 3); n++;
	feedbackField = XtCreateWidget("", xmTextWidgetClass, 
	    parent, args, n);
	
	XtAddCallback(feedbackField, XmNactivateCallback, 
	    (XtCallbackProc) ImgExaminerViewer::feedbackSizeFieldCB,
	    (XtPointer) this);
	
	feedbackLabel[1] = XtCreateWidget("pixels", 
	    xmLabelGadgetClass, parent, NULL, 0);
	
	// layout
	n = 0;
	XtSetArg(args[n], XmNleftAttachment,	XmATTACH_FORM); n++;
	XtSetArg(args[n], XmNleftOffset,	20); n++;
	XtSetArg(args[n], XmNtopAttachment,	XmATTACH_WIDGET); n++;
	XtSetArg(args[n], XmNtopWidget,		toggle); n++;
	XtSetArg(args[n], XmNtopOffset,		5); n++;
	XtSetValues(feedbackLabel[0], args, n);
	
	n = 0;
	XtSetArg(args[n], XmNleftAttachment,	XmATTACH_WIDGET); n++;
	XtSetArg(args[n], XmNleftWidget,	feedbackLabel[0]); n++;
	XtSetArg(args[n], XmNleftOffset,	5); n++;
	XtSetArg(args[n], XmNtopAttachment,	XmATTACH_OPPOSITE_WIDGET); n++;
	XtSetArg(args[n], XmNtopWidget,		feedbackLabel[0]); n++;
	XtSetValues(feedbackSizeWheel, args, n);
	
	n = 0;
	XtSetArg(args[n], XmNleftAttachment,	XmATTACH_WIDGET); n++;
	XtSetArg(args[n], XmNleftWidget,	feedbackSizeWheel); n++;
	XtSetArg(args[n], XmNleftOffset,	3); n++;
	XtSetArg(args[n], XmNtopAttachment,	XmATTACH_OPPOSITE_WIDGET); n++;
	XtSetArg(args[n], XmNtopWidget,		feedbackSizeWheel); n++;
	XtSetArg(args[n], XmNtopOffset,		-5); n++;
	XtSetValues(feedbackField, args, n);
	
	n = 0;
	XtSetArg(args[n], XmNleftAttachment,	XmATTACH_WIDGET); n++;
	XtSetArg(args[n], XmNleftWidget,	feedbackField); n++;
	XtSetArg(args[n], XmNleftOffset,	5); n++;
	XtSetArg(args[n], XmNbottomAttachment,	XmATTACH_OPPOSITE_WIDGET); n++;
	XtSetArg(args[n], XmNbottomWidget,	feedbackLabel[0]); n++;
	XtSetValues(feedbackLabel[1], args, n);
	
	// manage children
	XtManageChild(feedbackLabel[0]);
	XtManageChild(feedbackSizeWheel);
	XtManageChild(feedbackField);
	XtManageChild(feedbackLabel[1]);
    }
    else {
	// destroys the widgets
	XtDestroyWidget(feedbackLabel[1]);
	XtDestroyWidget(feedbackField);
	XtDestroyWidget(feedbackSizeWheel);
	XtDestroyWidget(feedbackLabel[0]);
    }
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Build the parent class widget, then register this widget.
//
// Use: protected
Widget
ImgExaminerViewer::buildWidget(Widget parent)
//
////////////////////////////////////////////////////////////////////////
{
    // Create the root widget and register it with a class name
    Widget w = ImgFullViewer::buildWidget(parent);
    
    // If first build, get resource values
    if (firstBuild) {
	// Full viewer registered the widget for us
	SoXtResource xr(w);
	SbBool boolvar;
	short val;
	
	if (xr.getResource("spinAnimation", "SpinAnimation", boolvar))
	    setAnimationEnabled(boolvar);
	if (xr.getResource("pointOfRotationAxes",
                           "PointOfRotationAxes", boolvar))
	    setFeedbackVisibility(boolvar);
	if (xr.getResource("axesSize", "AxesSize", val))
	    feedbackSize = val;
	
	firstBuild = FALSE;
    }
    
    return w;
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	add our own button to the existing list
//
// Use: virtual protected
void
ImgExaminerViewer::createViewerButtons(Widget parent)
//
////////////////////////////////////////////////////////////////////////
{
     // get the default buttons
    ImgFullViewer::createViewerButtons(parent);
    
    // allocate our buttons
    buttonList[CAM_PUSH] = new InvPixmapButton(parent, FALSE);
    buttonList[CAM_PUSH]->setIcon(so_xt_persp_bits, so_xt_icon_width, so_xt_icon_height);
    Widget w = buttonList[CAM_PUSH]->getWidget();
    XtAddCallback(w, XmNactivateCallback,
	(XtCallbackProc) ImgExaminerViewer::camPushCB, (XtPointer) this);
    
    // add this button to the list...
    viewerButtonWidgets->append(w);
}
////////////////////////////////////////////////////////////////////////
//
// Description:
//	read the point of interest geometry, which will be placed
//  right after the camera node (or right after the headlight
//  if the headlight is turned on).
//
// Use: private
void
ImgExaminerViewer::createFeedbackNodes()
//
////////////////////////////////////////////////////////////////////////
{
    // make sure we havn't built this yet...
    if (feedbackRoot)
	return;
    
    feedbackRoot	= new SoSeparator(1);
    feedbackSwitch	= new SoSwitch(3);
    feedbackTransNode	= new SoTranslation;
    feedbackScaleNode	= new SoScale;
    feedbackRoot->ref();
    feedbackRoot->addChild( feedbackSwitch );
    feedbackSwitch->addChild( feedbackTransNode );
    feedbackSwitch->addChild( feedbackScaleNode );
    SoInput in;
    in.setBuffer((void *)geometryBuffer, (size_t) strlen(geometryBuffer));
    SoNode *node;
    SbBool ok = SoDB::read(&in, node);
    if (ok && node != NULL)
	feedbackSwitch->addChild(node);

}
//
// redefine those generic virtual functions
//
const char *
ImgExaminerViewer::getDefaultWidgetName() const
{ return thisClassName; }

const char *
ImgExaminerViewer::getDefaultTitle() const
{ return "Examiner Viewer"; }

const char *
ImgExaminerViewer::getDefaultIconTitle() const
{ return "Examiner Viewer"; }



//
////////////////////////////////////////////////////////////////////////
// static callbacks stubs
////////////////////////////////////////////////////////////////////////
//

void
ImgExaminerViewer::camPushCB(Widget, ImgExaminerViewer *v, void *)
{ v->toggleCameraType(); }


void
ImgExaminerViewer::animationSensorCB(void *v, SoSensor *)
{ ((ImgExaminerViewer *) v)->doSpinAnimation(); }

void
ImgExaminerViewer::animPrefSheetToggleCB(Widget toggle, 
    ImgExaminerViewer *v, void *)
{
    v->setAnimationEnabled( XmToggleButtonGetState(toggle) );
}

void
ImgExaminerViewer::feedbackPrefSheetToggleCB(Widget toggle, 
    ImgExaminerViewer *v, void *)
{
    // show/hide the feedback
    v->setFeedbackVisibility( XmToggleButtonGetState(toggle) );
    
    // show/hide the extra size setting
    v->toggleFeedbackWheelSize(toggle);
}

void
ImgExaminerViewer::feedbackSizeWheelCB(Widget, ImgExaminerViewer *v, XtPointer *d)
{
    static SbBool firstDrag = TRUE;
    SgThumbWheelCallbackStruct *data = (SgThumbWheelCallbackStruct *) d;
    
    if (data->reason == XmCR_DRAG) {
	// for the first move, invoke the start callbacks
	if (firstDrag) {
	    v->interactiveCountInc();
	    firstDrag = FALSE;
	}
	
	// grow/shrink the feedback based on the wheel rotation
	v->feedbackSize *= powf(12, (data->value - v->feedbackSizeWheelVal) / 360.0);
	v->feedbackSizeWheelVal = data->value;
	
	// update the text field
	char str[15];
	sprintf(str, "%d", int(v->feedbackSize));
	XmTextSetString(v->feedbackField, str);
	
	// redraw since the wheel size isn't a field in the scene
	if (v->isViewing())
	    v->redraw();
    }
    else {
	// reason = XmCR_VALUE_CHANGED, invoke the finish callbacks
	v->interactiveCountDec();
	firstDrag = TRUE;
    }
}

void
ImgExaminerViewer::feedbackSizeFieldCB(Widget field, ImgExaminerViewer *v, void *)
{
    // get text value from the field
    char *str = XmTextGetString(field);
    int val;
    if ( sscanf(str, "%d", &val) && val > 0)
	v->setFeedbackSize(val);
    else
	val = int(v->feedbackSize);
    free(str);
    
    // reformat text field
    char valStr[15];
    sprintf(valStr, "%d", val);
    XmTextSetString(field, valStr);
    
    // make the text field loose the focus
    XmProcessTraversal(SoXt::getShellWidget(field), XmTRAVERSE_CURRENT);
}

// called when the viewer becomes visible/hidden - when hidden, make
// sure to temporary stop any ongoing animation (and restart it as soon
// as we become visible).
//
void
ImgExaminerViewer::visibilityChangeCB(void *pt, SbBool visible)
{
    ImgExaminerViewer *p = (ImgExaminerViewer *)pt;
    
    // only do this if we are/were spinning....
    if (! p->animatingFlag)
	return;
    
    if (visible) {
	// we now are visible again so reschedule the timer sensor
	p->animationSensor->schedule();
    }
    else
	// if hidden, unschedule the timer sensor, but don't change the
	// animatingFlag var to let us know we need to turn it back on
	// when we become visible....
	p->animationSensor->unschedule();
}
