/* $Log$
 * Revision 1.1  1994/04/12  13:39:31  zrfu0125
 * Initial revision
 * */


//**************************************************************************
//
// * Description    : Inventor viewer base class
//                   
// * Class(es)      : ImgViewer
//
// * inherited from : none
//
// * Author  : Dirk Rantzau / Andreas Werner
//
// * Last Modification : 25/04/96 by Roland Niemeier 
//                       (adaption of ImgViewer::setCurrentDrawStyle
//                        for volume rendering)
//
//**************************************************************************


#include <stdio.h>
#include <math.h>

#include <X11/Intrinsic.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>

#ifdef __sgi
#include <X11/extensions/SGIStereo.h>

// XSGISetStereoBuffer is not being defined, see bug 196992
extern "C" {
extern Status XSGISetStereoBuffer(
#if NeedFunctionPrototypes
	Display *,
	Window,
	int
#endif
);
}
#endif

#include <GL/gl.h>
#include <GL/glu.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/SbBox.h>
#include <Inventor/SbLinear.h>
#include <Inventor/SbPList.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/SoDB.h>
#include <Inventor/SoPath.h>
#include <Inventor/SoPath.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/Xt/SoXtClipboard.h>
#include "ImgViewer.h"
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/actions/SoSearchAction.h>
#include <Inventor/errors/SoDebugError.h>
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoComplexity.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoLightModel.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/nodes/SoResetTransform.h>
#include <Inventor/nodes/SoRotation.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSwitch.h>
#include <Inventor/nodes/SoScale.h>
#include <Inventor/sensors/SoFieldSensor.h>
#include <Inventor/sensors/SoTimerSensor.h>

#include "ImgRenderAction.h"

#define ANIM_RATE (1/60.0)  // animation frame rate
#define MINIMUM_NEAR_PLANE 0.001 // min near plane value relative to far plane


////////////////////////////////////////////////////////////////////////
//
//  Constructor.
//
// Use: protected

ImgViewer::ImgViewer(
    Widget parent,
    const char *name, 
    SbBool buildInsideParent,
    ImgViewer::Type t, 
    SbBool buildNow) 
	: SoXtRenderArea(
	    parent,
	    name, 
	    buildInsideParent, 
	    TRUE,   // getMouseInput
	    TRUE,   // getKeyboardInput
	    FALSE)  // buildNow
//
////////////////////////////////////////////////////////////////////////
{
    // init local vars
    addVisibilityChangeCallback(visibilityChangeCB, this);
    type = t;
    camera = NULL;
    cameraType = SoPerspectiveCamera::getClassTypeId();
    createdCamera = FALSE;
    viewingFlag = TRUE;
    interactiveFlag = FALSE;
    startCBList = new SoCallbackList;
    finishCBList = new SoCallbackList;
    interactiveCount = 0;
    bufferType = isDoubleBuffer() ? BUFFER_DOUBLE : BUFFER_SINGLE;
    stereoOffset = 3.0;
#ifdef __sgi
    useSGIStereoExt = FALSE;
#endif
    
    // init auto clipping stuff
    autoClipFlag = TRUE;
    // ??? we don't have a valid size yet....(alain)
    bboxAction = new SoGetBoundingBoxAction(SbVec2s(1,1));
    
    // copy/paste support
    clipboard = NULL;
    
    // init seek animation variables
    seekDistance = 50.0;
    seekDistAsPercentage = TRUE;
    seekModeFlag = FALSE;
    detailSeekFlag = TRUE;
    seekAnimTime = 2.0;
    seekAnimationSensor = new SoTimerSensor(ImgViewer::seekAnimationSensorCB, this);
    seekAnimationSensor->setInterval(ANIM_RATE);

    //
    // build the small internal graph (nodes used for draw style stuff)
    //
    sceneRoot		= new SoSeparator(4);
    drawStyleSwitch	= new SoSwitch(5);
    drawStyleNode   	= new SoDrawStyle;
    lightModelNode  	= new SoLightModel;
    colorNode		= new SoBaseColor;
    complexityNode	= new SoComplexity;
    sceneGraph 	    	= NULL;
    
    // note: we cannot setSceneGraph on the renderArea in the constructor
    // since it calls virtual functions, and all of our members aren't
    // initialized yet. We'll call it the first time our setSceneGraph
    // is called.
    sceneRoot->ref();
    sceneRoot->renderCaching.setValue(SoSeparator::OFF); // no caching there
    sceneRoot->renderCulling.setValue(SoSeparator::OFF); // no culling there
    sceneRoot->addChild(drawStyleSwitch);
    drawStyleSwitch->addChild(drawStyleNode);
    drawStyleSwitch->addChild(lightModelNode);
    drawStyleSwitch->addChild(colorNode);
    drawStyleSwitch->addChild(complexityNode);
    
    // set the draw style vars and fields that don't change
    stillDrawStyle = interactiveDrawStyle = VIEW_VOLPACK_SHADER;
    interactiveDrawStyle = VIEW_SAME_AS_STILL;
    drawStyleSwitch->whichChild = SO_SWITCH_NONE;
    
    drawStyleNode->setOverride(TRUE); // only use style field
    drawStyleNode->pointSize = 3.0;
    drawStyleNode->lineWidth.setIgnored(TRUE);
    drawStyleNode->linePattern.setIgnored(TRUE);
    
    lightModelNode->setOverride(TRUE);
    colorNode->setOverride(TRUE);
    
    complexityNode->setOverride(TRUE);
    complexityNode->textureQuality = 0; // always turn texture off under switch
    complexityNode->value = 0.15;
    
    addStartCallback(ImgViewer::drawStyleStartCallback);
    addFinishCallback(ImgViewer::drawStyleFinishCallback);
    
    //
    // headlightGroup - we have a rotation which keeps the headlight
    // moving whenever the camera moves,  and a reset xform so
    // that the rest of the scene is not affected by the first rot.
    // these leaves the direction field in the headlight open for the
    // user to edit, allowing for the direction to change w.r.t. the camera.
    //
    headlightGroup  = new SoGroup;
    headlightRot    = new SoRotation;
    headlightNode   = new SoDirectionalLight;
    headlightGroup->ref();
    headlightGroup->addChild(headlightRot);
    headlightGroup->addChild(headlightNode);
    headlightGroup->addChild(new SoResetTransform);
    headlightNode->direction.setValue(SbVec3f(.2, -.2, -.9797958971));
    headlightFlag = TRUE;
    
    // setup sensor on camera for headlight
    headlightSensor = new SoFieldSensor(ImgViewer::headlightSensorCB, this);
    headlightSensor->setPriority(1);  // make it happen before redraws
    
    // Build the widget tree, and let SoXtComponent know about our base widget.
    if (buildNow) {
	Widget w = buildWidget(getParentWidget());
	setBaseWidget(w);
    }
    
    
    // ##### Add-on: use own Render class ######## AW ###############

    StereoOn  = FALSE;
    LeftImage = FALSE;
    
    render = new ImgRenderAction (this,&StereoOn,&LeftImage,
                                  &Interactive,
                                  getViewportRegion()/*,
				  getGLRenderAction()->isInheriting()*/);
				  
    SoXtRenderArea::setGLRenderAction(render);
    
    // ###############################################################
}

////////////////////////////////////////////////////////////////////////
//
//    Destructor.
//
// Use: protected

ImgViewer::~ImgViewer()
//
////////////////////////////////////////////////////////////////////////
{
    // detach everything
    if ( sceneGraph != NULL )
	setSceneGraph(NULL);
    sceneRoot->unref();
    
    // delete everything
    delete headlightSensor;
    delete seekAnimationSensor;
    delete clipboard;
    delete bboxAction;
    delete startCBList;
    delete finishCBList;
    headlightGroup->unref();
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//  	Set a new user supplied scene graph.
//
// use: virtual public
//
void
ImgViewer::setSceneGraph(SoNode *newScene)
//
////////////////////////////////////////////////////////////////////////
{
    // if we haven't already given the render area a scene graph sceneRoot,
    // give it the scene graph now. This is a one shot deal, which
    // cannot be done in the constructor.
    if (SoXtRenderArea::getSceneGraph() == NULL)
	SoXtRenderArea::setSceneGraph(sceneRoot);
    
    //
    // detach everything that depends on the old sceneGraph
    //
    if ( sceneGraph != NULL ) {
	setCamera(NULL);
	sceneRoot->removeChild(sceneGraph);
    }
    
    sceneGraph = newScene;
    
    //
    // now assign the new sceneGraph, find or create the new camera 
    // and attach things back.
    //
    if ( sceneGraph != NULL ) {
	sceneRoot->addChild(sceneGraph);
	
	// search for first camera in the scene
	SoSearchAction sa;
	sa.setType(SoCamera::getClassTypeId());
	sa.setSearchingAll(FALSE); // don't look under off switches
	sa.apply(sceneGraph);
	
	SoCamera *newCamera = NULL;
	if (sa.getPath())
	    newCamera = (SoCamera *)((SoFullPath *)sa.getPath())->getTail();
	
	// if no camera found create one of the right kind...
	if ( newCamera == NULL ) {
	    
	    if (cameraType == SoPerspectiveCamera::getClassTypeId())
		newCamera = new SoPerspectiveCamera;
	    else if (cameraType == SoOrthographicCamera::getClassTypeId())
		newCamera = new SoOrthographicCamera;
	    else {
#ifdef DEBUG
		SoDebugError::post("ImgViewer::setSceneGraph",
		    "unknown camera type!");
#endif
		// ??? what should we do here ?
		newCamera = new SoPerspectiveCamera;
	    }
	    createdCamera = TRUE;
	    
	    if (type == ImgViewer::BROWSER)
		// add camera after drawstyle stuff
		sceneRoot->insertChild(newCamera, 1);
	    else {
		// check to make sure scene starts with at least a group node
		if ( sceneGraph->isOfType(SoGroup::getClassTypeId()) )
		    ((SoGroup *)sceneGraph)->insertChild(newCamera, 0);
		else {
		    // make scene start with a group node
		    SoGroup *group = new SoGroup;
		    group->addChild(newCamera);
		    group->addChild(sceneGraph);
		    sceneRoot->addChild(group);
		    sceneRoot->removeChild(sceneGraph);
		    sceneGraph = group;
		}
	    }
	    
	    newCamera->viewAll(sceneGraph, SbViewportRegion(getGlxSize()),4);
	}
	
	setCamera(newCamera);
    }
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//  	Return the user supplied scene graph.
//
// use: public, virtual

SoNode *
ImgViewer::getSceneGraph()
//
////////////////////////////////////////////////////////////////////////
{
    return sceneGraph;
}


////////////////////////////////////////////////////////////////////////
//
// Description:
//    Sets the camera to use.
//
// Use: virtual public
void
ImgViewer::setCamera(SoCamera *newCamera)
//
////////////////////////////////////////////////////////////////////////
{
    // check for trivual return
    if (camera == newCamera)
	return;
    
    //
    // detach everything that depended on the old camera
    //
    if ( camera != NULL ) {
	
        if (headlightFlag) {
	    setHeadlight(FALSE);
	    headlightFlag = TRUE;  // can later be turned on
        }
    	
	if (viewingFlag) {
	    setViewing(FALSE);
	    viewingFlag = TRUE;  // can later be turned on
	}
	
	// remove the camera if we created one outside of the
	// scene graph.
    	if (createdCamera && type == ImgViewer::BROWSER) {
	    if (sceneRoot->findChild(camera) >= 0)
    		sceneRoot->removeChild(camera);
	    createdCamera = FALSE;
	}
	
    	camera->unref();
    }
    
    camera = newCamera;
    
    //
    // attach everything that depends on the new camera
    //
    if ( camera != NULL) {
	camera->ref();
	
	if (headlightFlag) {
	    headlightFlag = FALSE;  // enables the routine to be called
	    setHeadlight(TRUE);
	}
	
	if (viewingFlag) {
	    viewingFlag = FALSE;  // enables the routine to be called
	    setViewing(TRUE);
	}
	
	saveHomePosition();
    }
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Set the camera type to create.
//
// Use: virtual public
void
ImgViewer::setCameraType(SoType type)
//
////////////////////////////////////////////////////////////////////////
{
    if (type.isDerivedFrom(SoPerspectiveCamera::getClassTypeId()) ||
	type.isDerivedFrom(SoOrthographicCamera::getClassTypeId()))
	cameraType = type;
#ifdef DEBUG
    else
	SoDebugError::post("ImgViewer::setCameraType",
			"unknown camera type!");
#endif
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    See the whole scene from the camera
//
// Use: public
void
ImgViewer::viewAll()
//
////////////////////////////////////////////////////////////////////////
{
    if ( camera != NULL )
	camera->viewAll(sceneGraph,SbViewportRegion(getGlxSize()));
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Sets the viewing mode.
//
// Use: virtual public
void
ImgViewer::setViewing(SbBool flag)
//
////////////////////////////////////////////////////////////////////////
{
    if (flag == viewingFlag)
	return;
    
    viewingFlag = flag;
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Sets the auto clipping mode
//
// Use: public
void
ImgViewer::setAutoClipping(SbBool flag)
//
////////////////////////////////////////////////////////////////////////
{
    if (autoClipFlag == flag)
	return;
    
    autoClipFlag = flag;
    
    // cause a redraw to correctly place the near and far plane now that
    // auto clipping is on.
    if (autoClipFlag)
	redraw();
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	sets stereo mode
//
// Use: virtual public
void
ImgViewer::setStereoViewing(SbBool flag)
//
////////////////////////////////////////////////////////////////////////
{
    if (flag == isStereoViewing())
	return;
    
    // First, check to see if the OpenGL stereo visual can be created
    setStereoBuffer(flag);
    
    // #####  New field: SbBool StereoOn = 1/0  #####  AW
    StereoOn = flag;
    
#ifdef __sgi
    // since OpenGL stereo failed, see if the SGI extension will work
    // by checking whether the X server supports it....
    int first_event, first_error;
    if (flag != isStereoViewing() && 
	XSGIStereoQueryExtension(getDisplay(), &first_event, &first_error)) {
	
	if (flag) {
	    // make sure the current window will support stereo
	    // ??? if we havn't been managed yet, just assume this visual
	    // ??? will support stereo viewing (see bug 
	    if (! getNormalWindow())
		useSGIStereoExt = TRUE;
	    else if (XSGIQueryStereoMode(getDisplay(), getNormalWindow()) !=
		X_STEREO_UNSUPPORTED)
		// stereo will be turned on in the rendering....
		useSGIStereoExt = TRUE;
	}
	else {
	    // turn stereo off on the window
	    useSGIStereoExt = FALSE;
	    
	    // clear the left/right buffers to prevent gost images from
	    // the other view...(until the user resets the monitor with setmon)
	    if (isRGBMode()) {
		SbColor color = getBackgroundColor();
		glClearColor(color[0], color[1], color[2], 0);
	    }
	    else
		glClearIndex(getBackgroundIndex());
	    
	    glDrawBuffer(GL_FRONT_AND_BACK);
	    
	    XSGISetStereoBuffer(getDisplay(), getNormalWindow(), STEREO_BUFFER_LEFT);
	    XSync(getDisplay(), False);
	    glClear(GL_COLOR_BUFFER_BIT);
	    
	    XSGISetStereoBuffer(getDisplay(), getNormalWindow(), STEREO_BUFFER_RIGHT);
	    XSync(getDisplay(), False);
	    glClear(GL_COLOR_BUFFER_BIT);
	    
	    glDrawBuffer( isDoubleBuffer() ? GL_BACK : GL_FRONT);
	}
	
	// now cause a redraw to see the affect since we havn't changed
	// the actual visual (unlike OpenGL)
	if (flag == isStereoViewing())
	    redraw();
    }
#endif
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	gets stereo mode
//
// Use: virtual public
SbBool
ImgViewer::isStereoViewing()
//
////////////////////////////////////////////////////////////////////////
{
#ifdef __sgi
    return (isStereoBuffer() || useSGIStereoExt);
#else
    // done in SoXtGLWidget
    return isStereoBuffer();
#endif
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    adds a directional light to the scene graph that is a headlight
//    positioned over the left shoulder of the camera.  It has a sensor
//    on the camera, so that it always is pointing in the same direction
//    as the camera.  The sensor is not delayed, so that it is "always"
//    accurate.
//
// Use: public
//
void
ImgViewer::setHeadlight(SbBool flag)
//
////////////////////////////////////////////////////////////////////////
{
    // check for trivual return
    if (camera == NULL || flag == headlightFlag) {
	headlightFlag = flag;
	return;
    }
    
    //
    // find the camera parent to insert/remove the headlight
    //
    SoSearchAction sa;
    sa.setNode(flag ? (SoNode *)camera: (SoNode *)headlightGroup);
    sa.apply(sceneRoot);
    SoFullPath *path = (SoFullPath *) sa.getPath();
    if (!path) {
#if DEBUG
	SoDebugError::post("ImgViewer::setHeadlight",
			    "ERROR: cannot find camera in graph");
#endif
	return;
    }
    SoGroup *parent = 
	    (SoGroup *)path->getNode(path->getLength() - 2);
    
    headlightFlag = flag;
    
    //
    // inserts/remove the headlight group node
    //
    if (headlightFlag) {


        //return if headlight is already there (this should be an error!)
	if (parent->findChild(headlightGroup) >= 0)
	  return;

	parent->insertChild(headlightGroup, parent->findChild(camera) + 1);

	if (isVisible())
	    activate();
    }
    else {
	parent->removeChild(headlightGroup);
	deactivate();
    }
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//   	Sets the drawing style.
//
// Use: public
void
ImgViewer::setDrawStyle(
	    ImgViewer::DrawType type, ImgViewer::DrawStyle style)
//
////////////////////////////////////////////////////////////////////////
{
    if (type == STILL) {
	if (stillDrawStyle == style)
	    return;
	if (style == VIEW_SAME_AS_STILL) {
#ifdef DEBUG
	    SoDebugError::post("ImgViewer::setDrawStyle", 
		"illegal VIEW_SAME_AS_STILL draw style passed for STILL !");
#endif
	    return; 
	}
	stillDrawStyle = style;
	
	if (! interactiveFlag || interactiveDrawStyle == VIEW_SAME_AS_STILL)
	    setCurrentDrawStyle(style);
    }
    else {
	// else it type == INTERACTIVE
	
	if (interactiveDrawStyle == style)
	    return;
	interactiveDrawStyle = style;
	
	if (interactiveFlag) {
	    if (interactiveDrawStyle == VIEW_SAME_AS_STILL)
		setCurrentDrawStyle(stillDrawStyle);
	    else
		setCurrentDrawStyle(style);
	}
    }
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//   	Sets the current drawing style. This only changes the nodes to
//  match what is passed (called from multiple places) and doesn't
//  affect the current state.
//
// Use: private
void
ImgViewer::setCurrentDrawStyle(ImgViewer::DrawStyle style)
//
////////////////////////////////////////////////////////////////////////
{
    if (style != VIEW_VOLPACK_SHADER)
	drawStyleSwitch->whichChild = SO_SWITCH_ALL;
    
    switch(style) {
	case VIEW_VOLPACK_SHADER:
	    drawStyleSwitch->whichChild = SO_SWITCH_NONE;
	    break;
	case VIEW_CALLBACK_SHADER:
	    drawStyleSwitch->whichChild = SO_SWITCH_NONE;
	    break;
	    
       //only VIEW_VOLPACK_SHADER and VIEW_CALLBACK_SHADER are currently used
	case VIEW_HIDDEN_LINE:
	    // texture is always off under the switch node.
	    // List only stuff common to both rendering passes
	    // (the rest is done when rendering)
	    drawStyleNode->style.setIgnored(FALSE);
	    drawStyleNode->pointSize.setIgnored(TRUE);
	    lightModelNode->model = SoLightModel::BASE_COLOR;
	    lightModelNode->model.setIgnored(FALSE);
	    complexityNode->type.setIgnored(TRUE);
	    complexityNode->value.setIgnored(TRUE);
	    break;
	    
        //only VIEW_VOLPACK_SHADER and VIEW_CALLBACK_SHADER are currently used
	case VIEW_NO_TEXTURE:
	//case VIEW_LOW_COMPLEXITY:
	    // texture is always off under the switch node
	    //drawStyleNode->style.setIgnored(TRUE);
	    //drawStyleNode->pointSize.setIgnored(TRUE);
	    //lightModelNode->model.setIgnored(TRUE);
	    //colorNode->rgb.setIgnored(TRUE);
	    //complexityNode->type.setIgnored(TRUE);
	    //complexityNode->value.setIgnored(style != VIEW_LOW_COMPLEXITY);
	    //break;
	    
	case VIEW_LOW_RES_LINE:
	    // texture is always off under the switch node
	    drawStyleNode->style.setIgnored(TRUE);
	    drawStyleNode->pointSize.setIgnored(TRUE);
	    lightModelNode->model.setIgnored(TRUE);
	    colorNode->rgb.setIgnored(TRUE);
	    complexityNode->type.setIgnored(TRUE);
	    complexityNode->value.setIgnored(style != VIEW_LOW_RES_LINE);
	    break;
	case VIEW_LINE:
	case VIEW_POINT:
	case VIEW_LOW_RES_POINT:
	    // texture is always off under the switch node
	    drawStyleNode->style = (style == VIEW_LINE || style == VIEW_LOW_RES_LINE) ? 
		SoDrawStyle::LINES : SoDrawStyle::POINTS;
	    drawStyleNode->style.setIgnored(FALSE);
	    drawStyleNode->pointSize.setIgnored(style != VIEW_POINT && style != VIEW_LOW_RES_POINT);
	    lightModelNode->model = SoLightModel::BASE_COLOR;
	    lightModelNode->model.setIgnored(FALSE);
	    colorNode->rgb.setIgnored(TRUE);
	    
	    // Force a lower complexity for the low res draw styles
	    // ??? this only works if the object didn't have
	    // ??? something lower in the first place...
	    if (style == VIEW_LOW_RES_LINE || style == VIEW_LOW_RES_POINT) {
		complexityNode->type = SoComplexity::OBJECT_SPACE;
		complexityNode->type.setIgnored(FALSE);
		complexityNode->value.setIgnored(FALSE);
	    }
	    else {
		complexityNode->type.setIgnored(TRUE);
		complexityNode->value.setIgnored(TRUE);
	    }
	    break;
	    
	case VIEW_BBOX:
	    // texture is always off under the switch node
	    drawStyleNode->style = SoDrawStyle::LINES;
	    drawStyleNode->style.setIgnored(FALSE);
	    drawStyleNode->pointSize.setIgnored(TRUE);
	    lightModelNode->model = SoLightModel::BASE_COLOR;
	    lightModelNode->model.setIgnored(FALSE);
	    colorNode->rgb.setIgnored(TRUE);
	    complexityNode->type = SoComplexity::BOUNDING_BOX;
	    complexityNode->type.setIgnored(FALSE);
	    complexityNode->value.setIgnored(TRUE);
	    break;
	    
	case VIEW_SAME_AS_STILL:
#ifdef DEBUG
	    SoDebugError::post("ImgViewer::setCurrentDrawStyle", "VIEW_SAME_AS_STILL was passed !");
#endif
	    break;
    }
    
    setZbufferState();
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//   	Gets the drawing style.
//
// Use: public
ImgViewer::DrawStyle
ImgViewer::getDrawStyle(ImgViewer::DrawType type)
//
////////////////////////////////////////////////////////////////////////
{
    return (type == STILL ? stillDrawStyle : interactiveDrawStyle);
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//   	redefine this routine to also correctly set the buffering type
//
// Use: virtual public
void
ImgViewer::setNormalVisual(XVisualInfo *vis)
//
////////////////////////////////////////////////////////////////////////
{
    // call parent class
    SoXtRenderArea::setNormalVisual(vis);
    
    // now update the buffering type
    if (isDoubleBuffer())
	setBufferingType(BUFFER_DOUBLE);
    else
	setBufferingType(BUFFER_SINGLE);
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//   	Sets the buffering style.
//
// Use: public
void
ImgViewer::setBufferingType(ImgViewer::BufferType type)
//
////////////////////////////////////////////////////////////////////////
{
    if (bufferType == type)
    	return;
    
    // remove interactive callback
    if (bufferType == BUFFER_INTERACTIVE) {
	removeStartCallback(ImgViewer::bufferStartCallback);
	removeFinishCallback(ImgViewer::bufferFinishCallback);
    }
    
    bufferType = type;
    
    switch(bufferType) {
	case BUFFER_SINGLE:
	    setDoubleBuffer(FALSE);
	    break;
	case BUFFER_DOUBLE:
	    setDoubleBuffer(TRUE);
	    break;
	case BUFFER_INTERACTIVE:
	    setDoubleBuffer(FALSE);
	    addStartCallback(ImgViewer::bufferStartCallback);
	    addFinishCallback(ImgViewer::bufferFinishCallback);
	    break;
    }
    
  
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Externally set the viewer into/out off seek mode (default OFF). Actual
//  seeking will not happen until the viewer decides to (ex: mouse click).
//
//  Note: setting the viewer out of seek mode while the camera is being
//  animated will stop the animation to the current location.
//
// use: virtual protected

void
ImgViewer::setSeekMode(SbBool flag)
//
////////////////////////////////////////////////////////////////////////
{
    if (!isViewing())
	return;
    
    seekModeFlag = flag;
    
    // check if seek is being turned off while seek animation is happening
    if ( !seekModeFlag && seekAnimationSensor->isScheduled() ) {
	seekAnimationSensor->unschedule();
	interactiveCountDec();
    }
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Adjust the camera clipping planes before a redraw.
//
// use: virtual public.

void
ImgViewer::actualRedraw()
//
////////////////////////////////////////////////////////////////////////
{
    if ( isAutoClipping() && ! isStereoViewing())
	adjustCameraClippingPlanes();
    
    //
    // Check to see if we are in stereo mode, if so draw the scene
    // twice with the camera offseted between the two views, else
    // do a simple redraw.
    //
    
    if ( isStereoViewing() && camera != NULL) {
	
	// Check the camera type, since stereo is different:
	//
	// Ortho camera: setereo is accomplished by simply rorating
	// the camera (around the point of interest) by 6 degree. 
	//
	// Perspective camera: we translate the camera and rotate
	// them to look at the same point of interest (idealy we also would
	// make sure the plane of convergence is exactly the same for
	// both perspective views, unfortunatly we cannot do this with
	// the current symetric view volumes).
	//
	
	// save the camera original values to restore the camera after
	// both views are rendered.
	SbVec3f	    camOrigPos = camera->position.getValue();
	SbRotation  camOrigRot = camera->orientation.getValue();
	float	    camOrigAspect = camera->aspectRatio.getValue();
	int	    camOrigViewportMapping = camera->viewportMapping.getValue();
	
	// get the camera focal point
	SbMatrix mx;
	mx = camOrigRot;
	SbVec3f forward( -mx[2][0], -mx[2][1], -mx[2][2]);
	float radius = camera->focalDistance.getValue();
	SbVec3f center = camOrigPos + radius * forward;
	
#ifdef __sgi
	//
	// if we are splitting the screen in half (loose vertical resolution)
	// then change the aspect ratio to squish the objects to make them 
	// look square again through the stereo glasses
	//
	if (useSGIStereoExt) {
	    SbVec2s windowSize = getGlxSize();
	    camera->aspectRatio = 0.5 * windowSize[0] / (float) windowSize[1];
	    camera->viewportMapping = SoCamera::LEAVE_ALONE;
	}
#endif
	
	//
	// change the camera for the LEFT eye view, and render
	//
	
	// ##### LeftImage-Field setzen ##### AW
	LeftImage = TRUE;
		
#ifdef __sgi
	if (useSGIStereoExt) {
	    XSGISetStereoBuffer(getDisplay(), getNormalWindow(), STEREO_BUFFER_LEFT);
	    XSync(getDisplay(), False);
	}
	else
#endif
	    glDrawBuffer( isDoubleBuffer() ? GL_BACK_LEFT : GL_FRONT_LEFT);
	// rotate the camera by - stereoOffset/2 degrees
	camera->orientation = 
	    SbRotation(SbVec3f(0, 1, 0), - stereoOffset * M_PI / 360.0) * camOrigRot;

        fprintf(stderr,"Changed camera orientation\n");
	
	// reposition camera to look at pt of interest
	mx = camera->orientation.getValue();
	forward.setValue( -mx[2][0], -mx[2][1], -mx[2][2]);
	camera->position = center - radius * forward;

        fprintf(stderr,"Changed camera position\n");
        
	
	adjustCameraClippingPlanes();
	doRendering();
	
	//
	// change the camera for the RIGHT eye view, and render
	//
	
	// ##### LeftImage-Field setzen ##### AW
	LeftImage = FALSE;
		
#ifdef __sgi
	if (useSGIStereoExt) {
	    XSGISetStereoBuffer(getDisplay(), getNormalWindow(), STEREO_BUFFER_RIGHT);
	    XSync(getDisplay(), False);
	}
	else
#endif
	    glDrawBuffer( isDoubleBuffer() ? GL_BACK_RIGHT : GL_FRONT_RIGHT);
	// rotate the camera by + stereoOffset/2 degrees
	camera->orientation = 
	    SbRotation(SbVec3f(0, 1, 0), stereoOffset * M_PI / 360.0) * camOrigRot;
	
	// reposition camera to look at pt of interest
	mx = camera->orientation.getValue();
	forward.setValue( -mx[2][0], -mx[2][1], -mx[2][2]);
	camera->position = center - radius * forward;
	
	adjustCameraClippingPlanes();
	doRendering();
	
	
	//
	// reset the camera original values now that we are done rendering
	// the stereo views.
	camera->enableNotify(FALSE); // don't cause a redraw
	camera->position = camOrigPos;
	camera->orientation = camOrigRot;
#ifdef __sgi
	if (useSGIStereoExt) {
	    camera->aspectRatio = camOrigAspect;
	    camera->viewportMapping = camOrigViewportMapping;
	}
#endif
	camera->enableNotify(TRUE);
	
#ifdef __sgi
	if (! useSGIStereoExt)
#endif
	// restore to draw to both buffer (viewer feedback)
	    glDrawBuffer( isDoubleBuffer() ? GL_BACK : GL_FRONT);
    }
    //
    // else not stereo viewing, so do the regular rendering....
    //
    else
	doRendering();
}


////////////////////////////////////////////////////////////////////////
//
// Description:
//	Do a multiple pass rendering if necessary, else simply call
//  SoXtRenderAre::actualRedraw() method.
//
// use: private
void
ImgViewer::doRendering()
//
////////////////////////////////////////////////////////////////////////
{
    // check if we need two pass rendering
    if ( (stillDrawStyle == VIEW_HIDDEN_LINE && 
	 (! interactiveFlag || interactiveDrawStyle == VIEW_SAME_AS_STILL)) ||
	(interactiveFlag && interactiveDrawStyle == VIEW_HIDDEN_LINE) ) {
	
	// ??? what do we do about highlights ??
	
	// the smaller the near clipping plane is relative to the far
	// plane, the smaller the zbuffer offset needs to be (because
	// the granularity will be pretty big). The closer the clipping
	// planes are relative to each other, the bigger the zbuffer offset
	// needs to be (because the zbuffer granularity will be small).
	// The scale factor was found empirically to work best with the
	// current settings of near/far.
	float zOffset = camera->nearDistance.getValue() / 
	    (40 * camera->farDistance.getValue());
	
	//
	// render the first pass as solid, using the background color
	// for the object base color.
	//
	
	drawStyleNode->style = SoDrawStyle::FILLED;
	colorNode->rgb = getBackgroundColor();
	colorNode->rgb.setIgnored(FALSE);
	
	// ??? this should match the SoXtRenderArea::actualRedraw()
	// ??? method exactly (apart for not clearing the z-buffer)
	glDepthRange(zOffset, 1); // enable wireframe to be draw on top
	//getSceneManager()->render(FALSE, FALSE);
	getSceneManager()->render(isClearBeforeRender(), TRUE);
	
	//
	// render the second pass as wireframe
	// (the first pass rendered the objects solid with base color
	// set to the background color to set the zbuffer values)
	//
	
	drawStyleNode->style = SoDrawStyle::LINES;
	colorNode->rgb.setIgnored(TRUE);
	
	// ??? this should match the SoXtRenderArea::actualRedraw()
	// ??? method exactly (apart for not clearing the color and z-buffer)
	glDepthRange(0,1-zOffset); // enable wireframe to be draw on top
	//getSceneManager()->render(FALSE, FALSE);
	
	//glDepthRange(0, 1); // restore the range
    }
    else
	// ??? this should match the SoXtRenderArea::actualRedraw()
	// ??? method exactly (apart for not clearing the z-buffer)
	getSceneManager()->render(isClearBeforeRender(), ! isZbufferOff());

}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Attach the sensor to the camera (if necessary).
//
// use: private

void
ImgViewer::activate()
//
////////////////////////////////////////////////////////////////////////
{
    // attach sensor to camera for headlight
    if (camera && headlightFlag && headlightSensor->getAttachedField() == NULL){
	headlightSensor->attach(&camera->orientation);
	
	// forces things to be in sync whith camera
	headlightSensor->schedule();
    }
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Detach the sensor.
//
// use: private

void
ImgViewer::deactivate()
//
////////////////////////////////////////////////////////////////////////
{
    headlightSensor->detach();
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Returns TRUE if the zbuffer should be off (based on the viewer
//  draw styles).
//
// Use: private

SbBool
ImgViewer::isZbufferOff()
//
////////////////////////////////////////////////////////////////////////
{
    DrawStyle style = (interactiveFlag ? interactiveDrawStyle : stillDrawStyle);
    if (interactiveFlag && interactiveDrawStyle == VIEW_SAME_AS_STILL)
	style = stillDrawStyle;
    
    // for these draw styles, turn the zbuffer off
    return (style == VIEW_LOW_RES_LINE || style == VIEW_LOW_RES_POINT 
	|| style == VIEW_BBOX);
}


////////////////////////////////////////////////////////////////////////
//
// Description:
//    Sets the zbuffer state on the current window. This is called whenever
//  the windows changes (called by SoXtGLWidget::widgetChanged()) or when
//  the viewer draw style changes.
//
// Use: private

void
ImgViewer::setZbufferState()
//
////////////////////////////////////////////////////////////////////////
{
    if (getNormalWindow() == NULL)
	return;
    
    glXMakeCurrent(getDisplay(), getNormalWindow(), getNormalContext());
    
    if (isZbufferOff())
	glDisable(GL_DEPTH_TEST);
    else
	glEnable(GL_DEPTH_TEST);
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    saves the camera values for later restore.
//
// Use: virtual public

void
ImgViewer::saveHomePosition()
//
////////////////////////////////////////////////////////////////////////
{
    if (camera == NULL)
	return;
    
    origPosition = camera->position.getValue();
    origOrientation = camera->orientation.getValue();
    origNearDistance = camera->nearDistance.getValue();
    origFarDistance = camera->farDistance.getValue();
    origFocalDistance = camera->focalDistance.getValue();
    
    // save camera height (changed by zooming)
    if (camera->isOfType(SoPerspectiveCamera::getClassTypeId()))
	origHeight = ((SoPerspectiveCamera *)camera)->heightAngle.getValue();
    else if (camera->isOfType(SoOrthographicCamera::getClassTypeId()))
	origHeight = ((SoOrthographicCamera *)camera)->height.getValue();
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    reset the camera to it's saved values.
//
// Use: virtual public

void
ImgViewer::resetToHomePosition()
//
////////////////////////////////////////////////////////////////////////
{
    if (camera == NULL)
	return;
    
    camera->position = origPosition;
    camera->orientation = origOrientation;
    camera->nearDistance = origNearDistance;
    camera->farDistance = origFarDistance;
    camera->focalDistance = origFocalDistance;
    
    // restore camera height (changed by zooming)
    if (camera->isOfType(SoPerspectiveCamera::getClassTypeId()))
	((SoPerspectiveCamera *)camera)->heightAngle.setValue(origHeight);
    else if (camera->isOfType(SoOrthographicCamera::getClassTypeId()))
	((SoOrthographicCamera *)camera)->height.setValue(origHeight);
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Process a common set of events which are shared accross all 
//  viewers. Returning TRUE if the event was processed.
//
// Use: protected
SbBool
ImgViewer::processCommonEvents(XAnyEvent *xe)
//
////////////////////////////////////////////////////////////////////////
{
    // check if the application wants to handle the event itself
    // instead of giving it to the viewer. This can be used to disable
    // some simple viewer functionality (like the Arrow keys or Esc key).
    // ??? this is a simple work around for bug #113991 - Xt translation
    // ??? tables would be better than dealing with events directly.
    if (isViewing() && invokeAppCB(xe))
	return TRUE;
    
    // check for special key which turns viewing on/off
    // regardless of whether viewing is off
    KeySym key = ! XK_Escape;
    if (xe->type == KeyPress) {
	key = XLookupKeysym((XKeyEvent *)xe, 0);
	if (key == XK_Escape) {
	    // check to see if the app should have a first short
	    if (! isViewing() && invokeAppCB(xe))
		return TRUE;
	    
	    // else toggle the viewing mode...
	    setViewing( !isViewing() );
	    return TRUE;
	}
    }
    
    // send the event to the scene graph if viewing is off
    if ( !isViewing() ) {
	// don't send that special Esc event twice...
	if (key != XK_Escape)
	    SoXtRenderArea::processEvent(xe);
	return TRUE;
    }
    
    // if no camera discard events
    if (camera == NULL)
	return TRUE;
    
    SbBool handled = TRUE;
    
    switch(xe->type) {
	case KeyPress:
	    switch ( key ) {
		case XK_Home:
		    resetToHomePosition();
		    break;
		case XK_s:
		    setSeekMode( !isSeekMode() );
		    // ??? this is kind of a hack, but it is needed
		    // ??? until a better solution is found
		    if ( isSeekMode() && interactiveCount != 0 ) {
			interactiveCount = 0;
			finishCBList->invokeCallbacks(this);
		    }
		    break;
		case XK_Left:
		case XK_Up:
		case XK_Right:
		case XK_Down:
		    arrowKeyPressed(key);
		    break;
		default:
		    handled = FALSE;
		    break;
	    }
	    break;
	    
	default:
	    handled = FALSE;
	    break;
    }
    
    return handled;
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Increment the intercative viewing counter.
//
// Use: protected.

void
ImgViewer::interactiveCountInc()
//
////////////////////////////////////////////////////////////////////////
{
    interactiveCount++;
    
    if (interactiveCount == 1)
      {
        Interactive=1;
	startCBList->invokeCallbacks(this);
      }
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Decrement the intercative viewing counter.
//
// Use: protected.

void
ImgViewer::interactiveCountDec()
//
////////////////////////////////////////////////////////////////////////
{
    if (interactiveCount > 0) {
	interactiveCount--;
	if (interactiveCount == 0)
          {
	    finishCBList->invokeCallbacks(this);
            Interactive=0;
          }
    }
}

////////////////////////////////////////////////////////////////////////
//
// Description:
// This routine is used by subclasses to initiate the seek animation. Given a
// screen mouse location, this routine will return the picked point
// and the normal at that point. It will also schedule the sensor to animate
// if necessary. The routine retuns TRUE if something got picked...
//
// Note: if detailSeek is on, the point and normal correspond to the exact 
//	 3D location under the cursor.
//	 if detailSeek if off, the object bbox center and the camera 
//	 orientation are instead returned.
//
// Use: protected.

SbBool
ImgViewer::seekToPoint(const SbVec2s &mouseLocation)
//
////////////////////////////////////////////////////////////////////////
{
    if (camera == NULL) {
	setSeekMode(FALSE);
	return FALSE;
    }
    
    // do the picking
    SoRayPickAction pick =  SoRayPickAction(SbViewportRegion(getGlxSize()));
    pick.setPoint(mouseLocation);
    pick.setRadius(1.0);
    pick.setPickAll(FALSE); // pick only the closest object
    pick.apply(sceneRoot);
    
    // makes sure something got picked
    SoPickedPoint *pp = pick.getPickedPoint();
    if ( pp == NULL ) {
	setSeekMode(FALSE);
	return FALSE;
    }
    
    //
    // Get picked point and normal if detailtSeek
    //
    if (detailSeekFlag) {
	
	seekPoint = pp->getPoint();
	seekNormal = pp->getNormal();
	
	// check to make sure normal points torward the camera, else
	// flip the normal around
	if ( seekNormal.dot(camera->position.getValue() - seekPoint) < 0 )
	    seekNormal.negate();
    }
    //
    // else get object bounding box as the seek point and the camera
    // orientation as the normal.
    //
    else {
	// get center of object's bounding box
	SoGetBoundingBoxAction bba = SoGetBoundingBoxAction(SbViewportRegion(getGlxSize()));
	bba.apply(pp->getPath());
	SbBox3f bbox = bba.getBoundingBox();
	seekPoint = bbox.getCenter();
	
	// keep the camera oriented the same way
	SbMatrix mx;
	mx = camera->orientation.getValue();
	seekNormal.setValue(mx[2][0], mx[2][1], mx[2][2]);
    }
    
    
    //
    // now check if animation sensor needs to be scheduled
    //
    
    computeSeekVariables = TRUE;
    if (seekAnimTime <= ANIM_RATE) {
	
	// jump to new location, no animation needed
	interpolateSeekAnimation(1.0);
    }
    else {
	// get animation starting time
	seekStartTime.setToTimeOfDay();
	
	// schedule sensor and call viewer start callbacks
	if ( ! seekAnimationSensor->isScheduled() ) {
	    seekAnimationSensor->schedule();
	    interactiveCountInc();
	}
    }
    
    return TRUE;    // successfull
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//  
// Use: virtual protected

void
ImgViewer::interpolateSeekAnimation(float t)
//
////////////////////////////////////////////////////////////////////////
{
    //
    // check if camera new and old position/orientation have already
    // been computed.
    //
    if (computeSeekVariables) {
	
	// save camera starting point
	oldCamPosition = camera->position.getValue();
	oldCamOrientation = camera->orientation.getValue();
	
	// now compute new camera position and rotation
	newCamPosition = seekPoint + 
	    seekNormal * camera->focalDistance.getValue();
	if ( isDetailSeek() ) {
	    
	    // get the camera forward view vector
	    SbMatrix mx;
	    mx = oldCamOrientation;
	    SbVec3f viewVector(-mx[2][0], -mx[2][1], -mx[2][2]);
	    
	    // calculate camera's orientation change
	    SbRotation changeOrient(viewVector, -seekNormal);
	    
	    newCamOrientation = oldCamOrientation * changeOrient;
	}
	else
	    newCamOrientation = oldCamOrientation;
	
	computeSeekVariables = FALSE;
    }
    
    
    //
    // Now position the camera according to the animation time
    //
    
    // use and ease-in ease-out approach
    float cos_t = 0.5 - 0.5 * fcos(t * M_PI);
    
    // get camera new rotation
    camera->orientation = SbRotation::slerp(oldCamOrientation, newCamOrientation, cos_t);
    
    // get camera new position
    camera->position = oldCamPosition + (newCamPosition - oldCamPosition) * cos_t;
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	Adjust the camera clipping planes based on the scene bounding 
//  box. (called before every redraws)
//
// use: private

void
ImgViewer::adjustCameraClippingPlanes()
//
////////////////////////////////////////////////////////////////////////
{
    if (camera == NULL)
	return;
    
    // get the scene bounding box
    bboxAction->setViewportRegion(SbViewportRegion(getGlxSize()));
    bboxAction->apply(sceneRoot);

    SbXfBox3f xfbbox = bboxAction->getXfBoundingBox();
    
    // get camera transformation and apply to xfbbox
    // to align the bounding box to camera space.
    // This will enable us to simply use the z values of the
    // transformed bbox for near and far plane values.
    SbMatrix mx;
    mx.setTranslate(- camera->position.getValue());
    xfbbox.transform(mx);
    mx = camera->orientation.getValue().inverse();
    xfbbox.transform(mx);
    
    // get screen align bbox and figure the near and far plane values
    SbBox3f bbox = xfbbox.project();
    // take negative value and opposite to what one might think 
    // because the camera points down the -Z axis
    float far = - bbox.getMin()[2];
    float near = - bbox.getMax()[2];
    
    // scene is behind the camera so don't change the planes
    if (far < 0)
	return;
    
    // check for minimum near plane value (Value will be negative 
    // when the camera is inside the bounding box).
    // Note: there needs to be a minimum near value for perspective
    // camera because of zbuffer resolution problem (plus the values
    // has to be positive). There is no such restriction for
    // an Orthographic camera (see behind you).
    if (!camera->isOfType(SoOrthographicCamera::getClassTypeId())) {
	if (near < (MINIMUM_NEAR_PLANE * far))
	    near = MINIMUM_NEAR_PLANE * far;
    }
    
    // give the near and far distances a little bit of slack in case
    // the object lies against the bounding box, otherwise the object
    // will be poping in and out of view (ex cube is the same as it's bbox)
    near *= 0.999;
    far *= 1.001;
    
    // finally assign camera plane values
    if (camera->nearDistance.getValue() != near)
	camera->nearDistance = near;
    if (camera->farDistance.getValue() != far)
	camera->farDistance = far;
}

////////////////////////////////////////////////////////////////////////
//
//  Copy the camera onto the clipboard.
//
//  Use: private
//
void
ImgViewer::copyView(Time eventTime)
//
////////////////////////////////////////////////////////////////////////
{
    if (camera == NULL)
	return;
    
    if (clipboard == NULL)
    	clipboard = new SoXtClipboard(getWidget());
    
    clipboard->copy(camera, eventTime);
}

////////////////////////////////////////////////////////////////////////
//
//  Retrieve the selection from the X server and paste it when it
//  arrives (in our pasteDone callback).
//
//  Use: private
//
void
ImgViewer::pasteView(Time eventTime)
//
////////////////////////////////////////////////////////////////////////
{
    if (clipboard == NULL)
    	clipboard = new SoXtClipboard(getWidget());
    
    clipboard->paste(eventTime, ImgViewer::pasteDoneCB, this);
}


////////////////////////////////////////////////////////////////////////
//
//  This is called by Xt when the data is ready to be pasted.
//
//  Use: static, private
//
void 
ImgViewer::pasteDoneCB(void *userData, SoPathList *pathList)
//
////////////////////////////////////////////////////////////////////////
{
    SoCamera *newCamera = NULL;
    
    // search for a camera in the paste data
    for (int i = 0; i < pathList->getLength(); i++) {
	SoFullPath *fullP = (SoFullPath *) (*pathList)[i];
	if (fullP->getTail()->isOfType(SoCamera::getClassTypeId())) {
	    newCamera = (SoCamera *) fullP->getTail();
	    break;
	}
    }

    if (newCamera != NULL)
	((ImgViewer *) userData)->changeCameraValues(newCamera);
    
    // We delete the callback data when done with it.
    delete pathList;
}

////////////////////////////////////////////////////////////////////////
//
//  Change the values of our camera to newCamera.
//??? animate from old values to new?
//
//  Use: virtual, protected
//
void 
ImgViewer::changeCameraValues(SoCamera *newCamera)
//
////////////////////////////////////////////////////////////////////////
{
    if (camera == NULL)
	return;
    
    // only paste cameras of the same type
    if (camera->getTypeId() != newCamera->getTypeId())
	return;

    // give our camera the values of the new camera
    camera->position	    = newCamera->position;
    camera->orientation	    = newCamera->orientation;
    camera->nearDistance    = newCamera->nearDistance;
    camera->farDistance	    = newCamera->farDistance;
    camera->focalDistance   = newCamera->focalDistance;

    // get the height or heightAngle
    if (camera->isOfType(SoPerspectiveCamera::getClassTypeId()))
	((SoPerspectiveCamera *)camera)->heightAngle = 
		((SoPerspectiveCamera *)newCamera)->heightAngle;
    else
	((SoOrthographicCamera *)camera)->height = 
		((SoOrthographicCamera *)newCamera)->height;
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Toggles the current camera type (perspective <--> orthographic)
//
//  Use: virtual protected
//
void
ImgViewer::toggleCameraType()
//
////////////////////////////////////////////////////////////////////////
{
    if (camera == NULL)
	return;
    
    // create the camera of the opposite kind and compute the wanted height
    // or heightAngle of the new camera.
    SoCamera *newCam;
    if (camera->isOfType(SoPerspectiveCamera::getClassTypeId())) {
	float angle = ((SoPerspectiveCamera *)camera)->heightAngle.getValue();
	float height = camera->focalDistance.getValue() * ftan(angle/2);
	newCam = new SoOrthographicCamera;
	((SoOrthographicCamera *)newCam)->height = 2 * height;
    }
    else if (camera->isOfType(SoOrthographicCamera::getClassTypeId())) {
	float height = ((SoOrthographicCamera *)camera)->height.getValue() / 2;
	float angle = fatan(height / camera->focalDistance.getValue());
	newCam = new SoPerspectiveCamera;
	((SoPerspectiveCamera *)newCam)->heightAngle = 2 * angle;
    }
    else {
#ifdef DEBUG
	SoDebugError::post("ImgViewer::toggleCameraType", "unknown camera type!");
#endif
	return;
    }
    
    newCam->ref();
    
    // copy common stuff from the old to the new camera
    newCam->viewportMapping = camera->viewportMapping.getValue();
    newCam->position = camera->position.getValue();
    newCam->orientation = camera->orientation.getValue();
    newCam->aspectRatio = camera->aspectRatio.getValue();
    newCam->focalDistance = camera->focalDistance.getValue();
    
    // search for the old camera and replace it by the new camera
    SoSearchAction sa;
    sa.setNode(camera);
    sa.apply(sceneRoot);
    SoFullPath *fullCamPath = (SoFullPath *) sa.getPath();
    if (fullCamPath) {
	SoGroup *parent = (SoGroup *)fullCamPath->getNode(fullCamPath->getLength() - 2);
	parent->insertChild(newCam, parent->findChild(camera));
	SoCamera *oldCam = camera;
	setCamera(newCam);
	
	// remove the old camera if it is still there (setCamera() might
	// have removed it) and set the created flag to true (for next time)
	if (parent->findChild(oldCam) >= 0)
	    parent->removeChild(oldCam);
	createdCamera = TRUE;
    }
#ifdef DEBUG
    else
	SoDebugError::post("ImgViewer::toggleCameraType", "camera not found!");
#endif
    
    newCam->unref();
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Called by the processCommonEvent routine whenever the arrow keys
//  are pressed. Translate the camera in the viewing plane in the arrow
//  direction half a screen at a time.
//
//  Use: private
//
void
ImgViewer::arrowKeyPressed(KeySym key)
//
////////////////////////////////////////////////////////////////////////
{
    if (camera == NULL)
	return;
    
    // get the camera near plane height value
    float dist;
    if (camera->isOfType(SoPerspectiveCamera::getClassTypeId())) {
	float angle = ((SoPerspectiveCamera *)camera)->heightAngle.getValue();
	float length = camera->nearDistance.getValue();
	dist = length * ftan(angle);
    }
    else if (camera->isOfType(SoOrthographicCamera::getClassTypeId()))
	dist = ((SoOrthographicCamera *)camera)->height.getValue();
    dist /= 2.0;
    
    // get camera right/left/up/down direction
    SbMatrix mx;
    mx = camera->orientation.getValue();
    SbVec3f dir;
    switch(key) {
	case XK_Up:
	    dir.setValue(mx[1][0], mx[1][1], mx[1][2]);
	    break;
	case XK_Down:
	    dir.setValue(-mx[1][0], -mx[1][1], -mx[1][2]); 
	    break;
	case XK_Right:
	    dir.setValue(mx[0][0], mx[0][1], mx[0][2]);
	    dist *= camera->aspectRatio.getValue();
	    break;
	case XK_Left:
	    dir.setValue(-mx[0][0], -mx[0][1], -mx[0][2]);
	    dist *= camera->aspectRatio.getValue();
	    break;
    }
    
    // finally reposition the camera
    camera->position = camera->position.getValue() + dist * dir;
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//	This routine is called when the widget/window is changed (single
//  to double, or for a new visual) so set the right zbuffer state (based
//  on the current draw style) for that new window.
//
//	HACK: redefines this routine from SoXtRenderArea to also check to make
//  sure that MyViewer::setBufferingType() is updated if users call 
//  SoXtGLWidget::setDoubleBuffer() instead.
//	Idealy setDoubleBuffer() would be virtual and redefined by the
//  viewers to call setBufferingType() instead - but that would brake
//  2.0.1 binary compatibility from 2.0 - see bug 198545. Instead we
//  know that widgetChanged() will be called at least when setDoubleBuffer()
//  is called.
//
//  Use: virtual protected
//
void
ImgViewer::widgetChanged(Widget w)
//
////////////////////////////////////////////////////////////////////////
{
    // set the correct zbuffer state on the new window
    if (w != NULL)
	setZbufferState();
    
    // call the parent class
    SoXtRenderArea::widgetChanged(w);
    
    //
    // HACK - now check to see if we need to update our setBufferingType() method
    // based on what the SoXtGLWidget::getDoubleBuffer() mode is.
    // (this won't busy loop because setDoubleBuffer() and setBufferingType()
    // will return if the state doesn't change).
    //
    
    switch(bufferType) {
	case BUFFER_SINGLE:
	    if (isDoubleBuffer()) // GLX is double, we think single
		setBufferingType(BUFFER_DOUBLE);
	    break;
	case BUFFER_DOUBLE:
	    if (! isDoubleBuffer()) // GLX is single, we think double
		setBufferingType(BUFFER_SINGLE);
	    break;
	case BUFFER_INTERACTIVE:
	    // don't do anything - this will toggle single/double
	    break;
    }
}



//
////////////////////////////////////////////////////////////////////////
// static callbacks stubs
////////////////////////////////////////////////////////////////////////
//



void
ImgViewer::bufferStartCallback(void *, ImgViewer *v)
{
    v->setDoubleBuffer(TRUE);
}

void
ImgViewer::bufferFinishCallback(void *, ImgViewer *v)
{
    v->setDoubleBuffer(FALSE);
}

void
ImgViewer::drawStyleStartCallback(void *, ImgViewer *v)
{
    v->interactiveFlag = TRUE;  // must happen first
    
    if (v->interactiveDrawStyle == v->stillDrawStyle ||
	v->interactiveDrawStyle == VIEW_SAME_AS_STILL)
	return;
    
    v->setCurrentDrawStyle(v->interactiveDrawStyle);
}

void
ImgViewer::drawStyleFinishCallback(void *, ImgViewer *v)
{
    v->interactiveFlag = FALSE;  // must happen first
    
    if (v->interactiveDrawStyle == v->stillDrawStyle ||
	v->interactiveDrawStyle == VIEW_SAME_AS_STILL)
	return;
    
    v->setCurrentDrawStyle(v->stillDrawStyle);
}

////////////////////////////////////////////////////////////////////////
//
//  Called whenever the camera changes values. The headlight is updated
//  accordingly to follow the camera new position.
//
//  Use: static private
//
void
ImgViewer::headlightSensorCB(void *p, SoSensor *)
//
////////////////////////////////////////////////////////////////////////
{
    ImgViewer *v = (ImgViewer *)p;
    if (!v->isVisible())
	return;
    
    // now update the headlight
    v->headlightRot->rotation.setValue(
	v->camera->orientation.getValue());
}

////////////////////////////////////////////////////////////////////////
//
//	Called whenever the seek animation sensor fires. Finds the amount 
//  of time since we started the seek and call the subclasses routine
//  to do the correct interpolation.
//
//  Use: static private
//
void
ImgViewer::seekAnimationSensorCB(void *p, SoSensor *)
//
////////////////////////////////////////////////////////////////////////
{
    ImgViewer *v = (ImgViewer *)p;
    
    // get the time difference
    SbTime time = SbTime::getTimeOfDay();
    SbTime timeDiff = time - v->seekStartTime;
    float t = float(timeDiff.getValue() / v->seekAnimTime);
    if (t > 1.0)
    	t = 1.0;
    
    // check if this gonna be the last one
    if ((1.0 - t) < 0.0001)
	t = 1.0;
    
    // call subclasses to interpolate the animation
    v->interpolateSeekAnimation(t);
    
    // stops seek if this was the last interval
    if (t == 1.0)
	v->setSeekMode(FALSE);
}

//
// called whenever the component becomes visibble or not
//
void
ImgViewer::visibilityChangeCB(void *pt, SbBool visible)
{
    ImgViewer *p = (ImgViewer *)pt;
    
    if (visible)
	p->activate();
    else
	p->deactivate();
}




//
////////////////////////////////////////////////////////////////////////
// viewer feedback convenience routines
////////////////////////////////////////////////////////////////////////
//



/*
 * Defines
 */

// color used in feedback
#define DARK_COLOR	glColor3ub(90, 90, 90)
#define LIGHT_COLOR	glColor3ub(230, 230, 230)

#define LINE_THIN   3	// line thickness used in feedback
#define	LINE_THICK  (LINE_THIN + 2)
#define CROSS 	    8	// size of cross hair at screen center for Roll
#define RADIUS	    15	// radius of center circle (in pix) for Roll
#define ANGLE_LEN   14   // angular size in degrees of Roll anchor


/*
 * Globals
 */

#define ARROW_SIZE  6.0	// size in pix of arrow head

// anchor arrow head description
static float arrow_data[3][3] = {
    -ARROW_SIZE, 0, 0,
    0, 2*ARROW_SIZE, 0,
    ARROW_SIZE, 0, 0
};


/*
 * Macros
 */

#define	DRAW_ARROW_MACRO    \
    DARK_COLOR;	\
    glBegin(GL_POLYGON);    \
    glVertex3fv(arrow_data[0]);	\
    glVertex3fv(arrow_data[1]);	\
    glVertex3fv(arrow_data[2]);	\
    glEnd();	\
    LIGHT_COLOR;	\
    glLineWidth(1); \
    glBegin(GL_LINE_LOOP);	\
    glVertex3fv(arrow_data[0]);	\
    glVertex3fv(arrow_data[1]);	\
    glVertex3fv(arrow_data[2]);	\
    glEnd();


////////////////////////////////////////////////////////////////////////
//
// Description:
//    Sets the default ortho projection when doing viewer feedback. The 
//  zbuffer/lighting is automatically turned off.
//
//  Use: static protected
//
void
ImgViewer::setFeedbackOrthoProjection(const SbVec2s &size)
//
////////////////////////////////////////////////////////////////////////
{
    // push the gl state to revert it back later....
    glPushAttrib(GL_DEPTH_BUFFER_BIT | GL_LIGHTING_BIT | GL_LINE_BIT);
    
    // ??? should we worry about restoring this matrix later ?
    glViewport(0, 0, size[0], size[1]);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, size[0], 0, size[1], -1, 1);
    
    // disable zbuffer and lighting....
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);

}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    restores the state that was changed when setFeedbackOrthoProjection()
//  is called.
//
//  Use: static protected
//
void
ImgViewer::restoreGLStateAfterFeedback()
//
////////////////////////////////////////////////////////////////////////
{
     // restore the gl state that were saved in setFeedbackOrthoProjection()
    glPopAttrib();
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    Draws a simple 2 colored cross at the given location.
//
//  Use: static protected
//
void
ImgViewer::drawViewerCrossFeedback(SbVec2s loc)
//
////////////////////////////////////////////////////////////////////////
{
    LIGHT_COLOR;
    glLineWidth(4);
    glBegin(GL_LINES);
    glVertex2s(loc[0]-CROSS, loc[1]);
    glVertex2s(loc[0]+CROSS, loc[1]);
    glVertex2s(loc[0], loc[1]-CROSS);
    glVertex2s(loc[0], loc[1]+CROSS);
    glEnd();
    
    DARK_COLOR;
    glLineWidth(2);
    glBegin(GL_LINES);
    glVertex2s(loc[0]-CROSS+1, loc[1]);
    glVertex2s(loc[0]+CROSS-1, loc[1]);
    glVertex2s(loc[0], loc[1]-CROSS+1);
    glVertex2s(loc[0], loc[1]+CROSS-1);
    glEnd();
}

////////////////////////////////////////////////////////////////////////
//
// Description:
//    draws the anchor roll feedback given the point of rotation and the
//  current mouse location.
//
//  Use: static protected
//
void
ImgViewer::drawViewerRollFeedback(SbVec2s center, SbVec2s loc)
//
////////////////////////////////////////////////////////////////////////
{
    // get angle and distance of mouse from center of rotation
    float ang, dist;
    float vx = loc[0] - center[0];
    float vy = loc[1] - center[1];
    if (vx==0 && vy==0) {
	ang = 0;
	dist = 0;
    }
    else {
	ang = atan2(vy, vx) * 180.0 / M_PI;
	dist = fsqrt(vx*vx + vy*vy);
    }
    float cirAng = -ang + 90; // gluPartialDisk() angle is REALLY backward !!
    
    static GLUquadricObj *quad = NULL;
    if (! quad)	quad = gluNewQuadric();
    
    // draw all of the circles (first inner, then outer)
    glTranslatef(center[0], center[1], 0);
    LIGHT_COLOR;
    gluDisk(quad, RADIUS, RADIUS+LINE_THICK, 20, 2);
    gluPartialDisk(quad, dist-2, dist+LINE_THICK-2, 20, 2, cirAng - ANGLE_LEN, 2 * ANGLE_LEN);
    DARK_COLOR;
    gluDisk(quad, RADIUS+1, RADIUS+LINE_THICK-1, 20, 2);
    gluPartialDisk(quad, dist-1, dist+LINE_THICK-3, 20, 2, cirAng - ANGLE_LEN, 2 * ANGLE_LEN);
    glTranslatef(-center[0], -center[1], 0); // undo the translation
    
    // draw connecting line from center to outer circle
    glLineWidth(LINE_THICK);
    LIGHT_COLOR;
    glBegin(GL_LINES);
    glVertex2s(center[0], center[1]);
    glVertex2s(loc[0], loc[1]);
    glEnd();
    glLineWidth(LINE_THIN);
    DARK_COLOR;
    glBegin(GL_LINES);
    glVertex2s(center[0], center[1]);
    glVertex2s(loc[0], loc[1]);
    glEnd();
    
    // draw the CCW arrow
    glPushMatrix();
    glTranslatef(center[0], center[1], 0);
    glRotatef(ang+ANGLE_LEN, 0, 0, 1);
    glTranslatef(dist, 0, 0);
    DRAW_ARROW_MACRO
    glPopMatrix();
    
    // draw the CW arrow
    glPushMatrix();
    glTranslatef(center[0], center[1], 0);
    glRotatef(ang-ANGLE_LEN, 0, 0, 1);
    glTranslatef(dist, 0, 0);
    glScalef(1, -1, 1);
    DRAW_ARROW_MACRO
    glPopMatrix();
}

