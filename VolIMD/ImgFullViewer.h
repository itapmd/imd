#ifndef  _INV_FULL_VIEWER_
#define  _INV_FULL_VIEWER_

/* $Id$ */

/* $Log$
 * Revision 1.1  1994/04/12  13:39:31  zrfu0125
 * Initial revision
 * */



#include <X11/Xm/Xm.h>
#include "ImgViewer.h"
#include <Inventor/SbPList.h>

// classes used
class	SoFieldSensor;
class	SoXtResource;
class	InvPixmapButton;

//////////////////////////////////////////////////////////////////////////////
//
//  Class: ImgFullViewer
//
//	The ImgFullViewer component class is the abstract base class for all
//  viewers which include a decoration around the rendering area. The 
//  decoration is made of thumbwheels, sliders and push/toggle buttons. There
//  is also a popup menu which includes all of the viewing options and methods.
//
//////////////////////////////////////////////////////////////////////////////

// C-api: abstract
// C-api: prefix=SoXtFullVwr
class ImgFullViewer : public ImgViewer {
  public:
    // This specifies what should be build by default in the constructor
    enum BuildFlag {
	BUILD_NONE	    = 0x00, 
	BUILD_DECORATION    = 0x01, 
	BUILD_POPUP	    = 0x02, 
	BUILD_ALL	    = 0xff, 
    };
    
    //
    // Show/hide the viewer component trims (default ON)
    //
    // C-api: name=setDecor
    void    	    setDecoration(SbBool onOrOff);
    // C-api: name=isDecor
    SbBool  	    isDecoration()	    { return decorationFlag; }
    
    //
    // Enable/disable the popup menu (default enabled)
    //
    // C-api: name=SetPopupEnabled
    void    	    setPopupMenuEnabled(SbBool trueOrFalse);
    // C-api: name=IsPopupEnabled
    SbBool  	    isPopupMenuEnabled()    { return popupEnabled; }
    
    //
    // Add/remove push buttons for the application, which will be placed 
    // in the left hand side decoration trim.
    // The add() method appends the button to the end of the list, while 
    // insert() places the button at the specified index (starting at 0).
    //
    // C-api: name=getAppPushBtnParent
    // returns the parent widget, which is needed when creating new buttons
    // NOTE that this will be NULL if the decoration is NOT created in the
    // constructor (see the BuildFlag) until it is shown.
    Widget	    getAppPushButtonParent() const { return appButtonForm; }
    // C-api: name=addAppPushBtn
    void	    addAppPushButton(Widget newButton);
    // C-api: name=insertAppPushBtn
    void	    insertAppPushButton(Widget newButton, int index);
    // C-api: name=removeAppPushBtn
    void	    removeAppPushButton(Widget oldButton);
    // C-api: name=findAppPushBtn
    int		    findAppPushButton(Widget oldButton)
				{ return appButtonList->find(oldButton); }
    // C-api: name=lengthAppPushBtn
    int		    lengthAppPushButton()
				{ return appButtonList->getLength(); }
    
    // C-api: name=getRAWidget
    Widget	    getRenderAreaWidget()   { return raWidget; }
    
    // redefine these from the base class
    virtual void    setViewing(SbBool onOrOff);
    virtual void    setCamera(SoCamera *cam);
    float	    getCameraZoom();
    virtual void    hide();
    
  protected:
    // Constructor/Destructor
    ImgFullViewer(
	Widget parent,
	const char *name, 
	SbBool buildInsideParent, 
	ImgFullViewer::BuildFlag flag, 
	ImgViewer::Type type, 
	SbBool buildNow);
    ~ImgFullViewer();
    
    // general decoration vars
    SbBool	    decorationFlag;
    Widget	    mgrWidget; // form which manages all other widgets
    Widget	    raWidget;  // render area widget
    Widget	    leftTrimForm, bottomTrimForm, rightTrimForm;
    
    // zoom slider vars
    Widget  	    zoomForm, zoomField, zoomSlider;
    SoFieldSensor   *zoomSensor; // attached to camera zoom field
    
    // thumb wheel vars
    Widget	    rightWheel, bottomWheel, leftWheel;
    char    	    *rightWheelStr, *bottomWheelStr, *leftWheelStr;
    float   	    rightWheelVal, bottomWheelVal, leftWheelVal;
    Widget  	    rightWheelLabel, bottomWheelLabel, leftWheelLabel;
    
    // widget list of viewer buttons
    SbPList	    *viewerButtonWidgets;
    
    // The button widget should be used as the parent widget
    // when creating new buttons
    Widget	    getButtonWidget() const { return appButtonForm; }
    
    // popup menu vars
    SbBool  	    popupEnabled;
    Widget  	    popupWidget, *popupToggleWidgets;
    Widget  	    *drawStyleWidgets, bufferStyleWidgets[3];
    char    	    *popupTitle;
    
    //
    // Build/destroy decoration routines
    //
    Widget 		buildWidget(Widget parent);
    void		buildLeftWheel(Widget parent);
    
    virtual void    	buildDecoration(Widget parent);
    virtual Widget    	buildLeftTrim(Widget parent);
    virtual Widget    	buildBottomTrim(Widget parent);
    virtual Widget    	buildRightTrim(Widget parent);
    virtual Widget    	buildZoomSlider(Widget parent);
    Widget		buildAppButtons(Widget parent);
    Widget		buildViewerButtons(Widget parent);
    virtual void	createViewerButtons(Widget parent);
    
    //
    // popup menu build routines
    //
    virtual void  	buildPopupMenu();
    virtual void    	destroyPopupMenu();
    void		setPopupMenuString(const char *name);
    Widget		buildFunctionsSubmenu(Widget popup);
    Widget		buildDrawStyleSubmenu(Widget popup);
    
    //
    // Preference sheet build routines
    //
    void	    setPrefSheetString(const char *name);
    virtual void    createPrefSheet();
    void	    createPrefSheetShellAndForm(Widget &shell, Widget &form);
    void	    createDefaultPrefSheetParts(Widget widgetList[], 
				int &num, Widget form);
    void	    layoutPartsAndMapPrefSheet(Widget widgetList[], 
				int num, Widget form, Widget shell);
    Widget	    createSeekPrefSheetGuts(Widget parent);
    Widget	    createSeekDistPrefSheetGuts(Widget parent);
    Widget	    createZoomPrefSheetGuts(Widget parent);
    Widget	    createClippingPrefSheetGuts(Widget parent);
    Widget	    createStereoPrefSheetGuts(Widget parent);
    
    // Subclasses SHOULD redefine these to do viewer specific tasks
    // since these do nothing in the base class. those routines are 
    // called by the thumb wheels whenever they rotate.
    virtual void    	rightWheelMotion(float);
    virtual void    	bottomWheelMotion(float);
    virtual void    	leftWheelMotion(float);
    
    // Subclasses can redefine these to add viewer functionality. They
    // by default call start and finish interactive viewing methods.
    // (defined in the base class).
    virtual void    	rightWheelStart();
    virtual void    	bottomWheelStart();
    virtual void    	leftWheelStart();
    virtual void    	rightWheelFinish();
    virtual void    	bottomWheelFinish();
    virtual void    	leftWheelFinish();
    
    // Subclasses SHOULD set those wheel string to whatever functionality
    // name they are redefining the thumb wheels to do. Default names are
    // "Motion X, "Motion Y" and "Motion Z" for bottom, left and right wheels.
    void		setBottomWheelString(const char *name);
    void		setLeftWheelString(const char *name);
    void		setRightWheelString(const char *name);
    
    // Subclasses SHOULD redefine the openViewerHelpCard() routine to bring 
    // their help card (called by push button and popup menu entry).
    // They can simply call SoXtComponent::openHelpCard() to open their
    // specific help card.
    virtual void	openViewerHelpCard();
    
    // redefine this to load the popup menu colormap
    virtual void    afterRealizeHook();

 private:
    SbBool	    firstBuild; // set FALSE after buildWidget called once
    
    // app button vars
    Widget	    appButtonForm;
    SbPList	    *appButtonList;
    void	    doAppButtonLayout(int start);

    // this is called the first time the widget is built
    void	    getResources(SoXtResource *xr);
    
    // popup menu callbacks
    static void     popMenuCallback(Widget, ImgFullViewer *, XEvent *, Boolean *);
    static void     drawStyleMenuPick(Widget, int id, void *);
    static void     bufferStyleMenuPick(Widget, int id, void *);
    static void     menuPick(Widget, int id, XmAnyCallbackStruct *);
    static void	    menuDisplay(Widget, ImgFullViewer *, XtPointer *);
    
    // zoom slider vars and routines
    SbVec2f	    zoomSldRange;
    void    	    setCameraZoom(float zoom);
    void	    setZoomSliderPosition(float zoom);
    void	    setZoomFieldString(float zoom);
    static void	    zoomSliderCB(Widget, ImgFullViewer *, XtPointer *);
    static void	    zoomFieldCB(Widget, ImgFullViewer *, XtPointer *);
    static void	    zoomSensorCB(void *, SoSensor *);
    
    // push buttons vars and callbacks
    InvPixmapButton  *buttonList[10];
    static void	    pushButtonCB(Widget, int id, void *);
    
    // pref sheet variables
    Widget	    prefSheetShellWidget;
    char	    *prefSheetStr;
    static void	    prefSheetDestroyCB(Widget, ImgFullViewer *, void *);
    
    // seek pref sheet callbacks
    static void	    seekPrefSheetFieldCB(Widget, ImgFullViewer *, void *);
    static void	    seekPrefSheetToggle1CB(Widget, Widget, void *);
    static void	    seekPrefSheetToggle2CB(Widget, Widget, void *);
    
    // zoom pref sheet callbacks
    static void	    zoomPrefSheetMinFieldCB(Widget, ImgFullViewer *, void *);
    static void	    zoomPrefSheetMaxFieldCB(Widget, ImgFullViewer *, void *);
    
    // seek dist pref sheet vars and callbacks
    int		    seekDistWheelVal;
    Widget	    seekDistField;
    static void	    seekDistWheelCB(Widget, ImgFullViewer *v, XtPointer *d);
    static void	    seekDistFieldCB(Widget, ImgFullViewer *, void *);
    static void	    seekDistPercPrefSheetToggleCB(Widget, Widget, void *);
    static void	    seekDistAbsPrefSheetToggleCB(Widget, Widget, void *);
    
    // clipping plane pref sheet callbacks and variables
    Widget	    clipWheelForm;
    int		    clipNearWheelVal, clipFarWheelVal;
    Widget	    clipNearField, clipFarField;
    static void	    clipPrefSheetToggleCB(Widget, Widget, void *);
    static void	    clipNearWheelCB(Widget, ImgFullViewer *v, XtPointer *d);
    static void	    clipFarWheelCB(Widget, ImgFullViewer *v, XtPointer *d);
    static void	    clipFieldCB(Widget, ImgFullViewer *, void *);
    
    // stereo viewing pref sheet callbacks and variables
    Widget	    stereoWheelForm, stereoField, stereoLabel;
    int		    stereoWheelVal;
    static void	    stereoPrefSheetToggleCB(Widget, Widget, void *);
    static void	    stereoWheelCB(Widget, ImgFullViewer *, XtPointer *);
    static void	    stereoFieldCB(Widget, ImgFullViewer *, void *);
    
    // error dialog destroy callback
    static void	    destroyDialogCB(Widget, void *, void *);
    
    // zoom slider sensor activation 
    static void visibilityChangeCB(void *pt, SbBool visible);
    void        activate();		// connects the sensor
    void        deactivate();		// disconnects the sensor
    
    // callback functions for thumbwheels
    static  void    rightWheelCB(Widget, ImgFullViewer *v, XtPointer *d);
    static  void    bottomWheelCB(Widget, ImgFullViewer *v, XtPointer *d);
    static  void    leftWheelCB(Widget, ImgFullViewer *v, XtPointer *d);
};


#endif  /* _INV_FULL_VIEWER_ */
