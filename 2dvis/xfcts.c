#include <stdio.h>
#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>
#include <X11/Xaw/Dialog.h>
#include <X11/Shell.h>
#include <X11/Xaw/Label.h>
#include <X11/Xaw/Paned.h>
#include <X11/Xaw/Box.h>
#include <X11/Xaw/Command.h>
#include <X11/Xaw/Toggle.h> 
#include "globals.h"
#include "xfcts.h"

extern Widget toplevel;
extern Display *dpy;
int radiovalue;

struct drag_struct {
  GC   gc;
  GC   xor_gc;
  char *object;
} t_drag_client;

XtEventHandler t_press_event(Widget w,struct drag_struct *drag,XButtonEvent *ev) {
  int mousex,button;
  button=ev->button;
  if (button!=1) return;
  mousex=ev->x;
  temperature=mousex*.005;
}

/* some local identifiers */
static int ReturnCode;         /* used by standard dialogs */

static XtCallbackProc rgback(Widget w,struct radio_struct *rs,XtPointer call) {
       *(rs->variable) = rs->value;
       radiovalue = rs->value;
}

/*
* callback for standard dialog
* just copies the button-code which is passed as ClientData to a
* public identifier
*/
static void CB_OK(Widget W, XtPointer ClientData, XtPointer CallData) {
  ReturnCode = (int) ClientData;
}

/* creates a dialog box as a child of ... (see arguments) */
Widget CreateDialogBox(char *LabelText,
		       DialogButtonTypePtr Buttons, 
		       size_t Count, char *title) {
  Widget  dialog, popup;
  
  /* create the popup shell */
  popup = XtVaCreatePopupShell(title, transientShellWidgetClass,
			       toplevel,
			       XtNinput, True,
			       XtNallowShellResize, True,
			       NULL);
  dialog = XtVaCreateManagedWidget("dialog", dialogWidgetClass,
				   popup,
				   XtNlabel, LabelText,
				   NULL);

  /* now we add the buttons and install the accelerators */
  AddButtons(dialog, NULL, Buttons, Count);
  while (Count--)
    XtInstallAccelerators(dialog, Buttons[Count].W);
  return(popup);
}

/* creates a dialog box as a child of ... (see arguments) */
Widget CreateRadioBox(char *LabelText,char *title,DialogButtonTypePtr Buttons, 
		       size_t Count,String *names, 
		       int *values,int *variable,int initial) {
  Widget  dialog, popup;
  struct    radio_struct *rs;
  Widget    togs[100];
  int  n,index,i,ntogs,len;
  char *init_name;
  Widget    box,group;
  Arg  wargs[10];

  /* create the popup shell */
  popup = XtVaCreatePopupShell(title, transientShellWidgetClass,
			       toplevel,
			       XtNinput, True,
			       XtNallowShellResize, True,
			       NULL);
  dialog = XtVaCreateManagedWidget("dialog", panedWidgetClass,
				   popup,
				   XtNlabel, LabelText,
				   NULL);
  box = XtCreateManagedWidget("box",boxWidgetClass,dialog,
			      NULL, 0);


  ntogs = 0;

  while(names[ntogs] != 0) {
    togs[ntogs] = XtCreateWidget(names[ntogs],
				 toggleWidgetClass,box,NULL,0);
    n=0;
    XtSetArg(wargs[n], XtNmin, 25); n++;
    XtSetArg(wargs[n], XtNmax, 75); n++;
    XtSetArg(wargs[n], XtNpreferredPaneSize, 50); n++;
    XtSetValues(togs[ntogs],wargs,n);

    
    if(ntogs == 0)
      group = togs[0];

    XawToggleChangeRadioGroup(togs[ntogs],group);
    
    rs = (struct radio_struct *) malloc(sizeof *rs);
    rs->variable = variable;
    rs->value = values[ntogs];

    XtAddCallback(togs[ntogs],XtNcallback,rgback,rs);

    ntogs++;
  }

  XtManageChildren(togs,ntogs);

  index = -1;
  for(i=0; i < ntogs; i++)
    if(values[i] == initial)
      index = i;
  
  if(index  <  0)
    index = 0;
  
  *variable = values[index];
  
  len = strlen(names[index]) + 1;
  init_name = (char *)malloc(len);
  strcpy(init_name,names[index]);

  XtSetArg(wargs[0],XtNradioData,init_name);
  XtSetValues(togs[index],wargs,1);
  
  XawToggleSetCurrent(group,init_name);

  /* now we add the buttons and install the accelerators */
  AddButtons(dialog, NULL, Buttons, Count);
  while (Count--)
    XtInstallAccelerators(dialog, Buttons[Count].W);
  return(popup);
}

/* creates a dialog box as a child of ... (see arguments) */
Widget CreateSlideBox(char *LabelText,
		       DialogButtonTypePtr Buttons, 
		       size_t Count, char *title) {
  Widget  dialog, popup, box,rect;
  int n;
  Arg wargs[3];

  /* create the popup shell */
  popup = XtVaCreatePopupShell(title, transientShellWidgetClass,
			       toplevel,
			       XtNinput, True,
			       XtNallowShellResize, True,
			       NULL);
  dialog = XtVaCreateManagedWidget("dialog", panedWidgetClass,
				   popup,
				   XtNlabel, LabelText,
				   NULL);
  box = XtCreateManagedWidget("box",boxWidgetClass,dialog,
			      NULL, 0);
  rect = XtCreateManagedWidget("rect",boxWidgetClass,box,
			      NULL, 0);

  n = 0;
  XtSetArg(wargs[n], XtNheight, 20); n++;
  XtSetArg(wargs[n], XtNwidth, 200); n++;
  XtSetValues(rect, wargs, n);

  XtAddEventHandler(popup, ButtonPressMask, FALSE,
		    t_press_event, &t_drag_client);

  AddButtons(dialog, NULL, Buttons, Count);
  while (Count--)
    XtInstallAccelerators(dialog, Buttons[Count].W);
  return(popup);
}

/* adds the buttons to a dialog */
Widget AddButtons(Widget Parent, Widget Top,
		  DialogButtonTypePtr Buttons, size_t Count) {
  int             i;
  
  for (i = 0; i < Count; i++) {
    /* skip button if there's no label,
     * used for dialogs without default button
     */
    if (!Buttons[i].Label)
      continue;
    Buttons[i].W=XtVaCreateManagedWidget(Buttons[i].Name,
					 commandWidgetClass,
					 Parent,
					 XtNlabel, Buttons[i].Label,
					 XtNfromHoriz, i ? Buttons[i-1].W : NULL,
					 XtNfromVert, Top,
					 XtNresizable, True,
					 NULL);
    
    XtAddCallback(Buttons[i].W,
		  XtNcallback, Buttons[i].Callback, Buttons[i].ClientData);
  }
  return(Buttons[Count-1].W);
}

/* prints out a question returns the answer */
int QuestionDialog(char *MessageText) {
  Widget popup;
  static DialogButtonType buttons[] =
    { { "defaultButton", "  OK  ", CB_OK, (XtPointer) OK_BUTTON, NULL },{ "", "  CANCEL  ", CB_OK, (XtPointer) CANCEL_BUTTON, NULL }};
  
  /* create dialog box */
  popup = CreateDialogBox(MessageText, buttons, ENTRIES(buttons),MessageText);
  StartDialog(popup);
  /* wait for dialog to complete */
  DialogEventLoop(&ReturnCode);
  EndDialog(popup);
  return ReturnCode;
}

/* color encoding dialog */
int SaveImageDialog(char *MessageText) {
  Widget popup;
  int var;
  String names[] = {
       " gif  ",
       " vrml ",
       " pdb",
       0
  };
  int values[] = {
       1,
       2,
       3,
  };
  static DialogButtonType buttons[] =
  { { "defaultButton", "     OK     ", CB_OK, (XtPointer) OK_BUTTON, NULL }};

  /* create radio box */
  popup = CreateRadioBox(MessageText,MessageText, buttons, ENTRIES(buttons),names,values,&var,savimg_mode+1);
  StartDialog(popup);
  /* wait for radio box to complete */
  DialogEventLoop(&ReturnCode);
  EndDialog(popup);

  return radiovalue;
}

/* color encoding dialog */
int ColorEncodingDialog(char *MessageText) {
  Widget popup;
  int var;
  String names[] = {
       "  nothing  ",
       "    type   ",
       "pot. energy",
       "kin. energy",
       "#neighbours",
       "   number  ",
       0
  };
  int values[] = {
       1,
       2,
       3,
       4,
       5,
       6
  };
  static DialogButtonType buttons[] =
  { { "defaultButton", "     OK     ", CB_OK, (XtPointer) OK_BUTTON, NULL }};

  /* create radio box */
  popup = CreateRadioBox(MessageText,MessageText, buttons, ENTRIES(buttons),names,values,&var,col_mode+1);
  StartDialog(popup);
  /* wait for radio box to complete */
  DialogEventLoop(&ReturnCode);
  EndDialog(popup);

  return radiovalue;
}

/* size encoding dialog */
int SizeEncodingDialog(char *MessageText) {
  Widget popup;
  int var;
  String names[] = {
       "nothing",
       "  type ",
       0
  };
  int values[] = {
       1,
       2
  };
  static DialogButtonType buttons[] =
  { { "defaultButton", "     OK     ", CB_OK, (XtPointer) OK_BUTTON, NULL }};

  /* create radio box */
  popup = CreateRadioBox(MessageText,MessageText, buttons, ENTRIES(buttons),names,values,&var,size_mode+1);
  StartDialog(popup);
  /* wait for radio box to complete */
  DialogEventLoop(&ReturnCode);
  EndDialog(popup);

  return radiovalue;
}

/* Bond mode dialog */
int BondTypeDialog(char *MessageText) {
  Widget popup;
  int var;
  String names[] = {
       " static",
       "dynamic",
       0
  };
  int values[] = {
       1,
       2,
  };
  static DialogButtonType buttons[] =
  { { "defaultButton", "     OK     ", CB_OK, (XtPointer) OK_BUTTON, NULL }};

  /* create radio box */
  popup = CreateRadioBox(MessageText,MessageText, buttons, ENTRIES(buttons),names,values,&var,bond_type+1);
  StartDialog(popup);
  /* wait for radio box to complete */
  DialogEventLoop(&ReturnCode);
  EndDialog(popup);

  return radiovalue;
}

/* Atom mode dialog */
int AtomModeDialog(char *MessageText) {
  Widget popup;
  int var;
  String names[] = {
       "   off    ",
       "draw atoms",
       0
  };
  int values[] = {
       1,
       2,
  };
  static DialogButtonType buttons[] =
  { { "defaultButton", "     OK     ", CB_OK, (XtPointer) OK_BUTTON, NULL }};

  /* create radio box */
  popup = CreateRadioBox(MessageText,MessageText, buttons, ENTRIES(buttons),names,values,&var,atom_mode+1);
  StartDialog(popup);
  /* wait for radio box to complete */
  DialogEventLoop(&ReturnCode);
  EndDialog(popup);

  return radiovalue;
}

/* Bond mode dialog */
int BondModeDialog(char *MessageText) {
  Widget popup;
  int var;
  String names[] = {
       "off",
       "just count",
       "draw bonds",
       0
  };
  int values[] = {
       1,
       2,
       3,
  };
  static DialogButtonType buttons[] =
  { { "defaultButton", "     OK     ", CB_OK, (XtPointer) OK_BUTTON, NULL }};

  /* create radio box */
  popup = CreateRadioBox(MessageText,MessageText, buttons, ENTRIES(buttons),names,values,&var,bond_mode+1);
  StartDialog(popup);
  /* wait for radio box to complete */
  DialogEventLoop(&ReturnCode);
  EndDialog(popup);

  return radiovalue;
}

/* color encoding dialog */
int QuasiSwitchDialog(char *MessageText) {
  Widget popup;
  int var;
  String names[] = {
       "    periodic   ",
       "quasip. 2 types",
       "quasip. 3 types",
       "genburgers     ",
       0
  };
  int values[] = {
       1,
       2,
       3,
       4,
  };
  static DialogButtonType buttons[] =
  { { "defaultButton", "     OK     ", CB_OK, (XtPointer) OK_BUTTON, NULL }};

  /* create radio box */
  popup = CreateRadioBox(MessageText,MessageText, buttons, ENTRIES(buttons),names,values,&var,qp+1);
  StartDialog(popup);
  /* wait for radio box to complete */
  DialogEventLoop(&ReturnCode);
  EndDialog(popup);

  return radiovalue;
}

/* prints out a message */
void TemperatureDialog(char *MessageText) {
  Widget popup;
  static DialogButtonType buttons[] =
  { { "defaultButton", "     OK     ", CB_OK, (XtPointer) OK_BUTTON, NULL }};
  
  /* create dialog box */
  popup = CreateSlideBox(MessageText, buttons, ENTRIES(buttons),MessageText);
  StartDialog(popup);
  /* wait for dialog to complete */
  DialogEventLoop(&ReturnCode);
  EndDialog(popup);
}

/* prints out a message */
void MessageDialog(char *MessageText) {
  Widget popup;
  static DialogButtonType buttons[] =
  { { "defaultButton", "     OK     ", CB_OK, (XtPointer) OK_BUTTON, NULL }};
  
  /* create dialog box */
  popup = CreateDialogBox(MessageText, buttons, ENTRIES(buttons),MessageText);
  StartDialog(popup);
  /* wait for dialog to complete */
  DialogEventLoop(&ReturnCode);
  EndDialog(popup);
}

/*
* pops up a dialog
* saves output window in pixmap
* first realize the dialog and
* focus keyboard events from the application to it
*/
void StartDialog(Widget Popup) {
  /*  SaveOutputWindow();*/
  XtRealizeWidget(Popup);
  CenterDialog(Popup, True, True);
  XtPopup(Popup, XtGrabExclusive);
}

/* removes a dialog from screen and from memory
 * pixmap is released by next dialog */
void EndDialog(Widget Popup) {
  XtPopdown(Popup);
  XtDestroyWidget(Popup);
  XFlush(dpy);
  XSync(dpy, False);
}

/* waits until returncode is differnt from -1
 * see also CB_OK() */
int DialogEventLoop(int *Code) {
  XEvent  event;
  
  for(*Code = -1; *Code == -1;) {
    XtNextEvent(&event);
    XtDispatchEvent(&event);
  }
  return(*Code);
}


/* centers the widget around the current cursor position */
static void CenterDialog(Widget Popup, Boolean CenterX, Boolean CenterY) {
  Window root, child;
  int    root_x, root_y,child_x, child_y;
  unsigned int    mask;
  Dimension               width, height;

  /* get current pointer position relativ to it's parent */
  XQueryPointer(dpy, XtWindow(Popup),
		&root, &child,
		&root_x, &root_y,
		&child_x, &child_y,
		&mask);
  
  /* get the dialogs size */
  XtVaGetValues(Popup,
		XtNheight, &height,
		XtNwidth, &width,
		NULL);  
  
  /* make sure position is inside our screen */
  if (CenterX)
    if ((root_x -= (width/2)) < 0)
      root_x = 0;
    else
      if (root_x > WidthOfScreen(XtScreen(Popup)) -width)
	root_x =  WidthOfScreen(XtScreen(Popup)) -width;
  
  if (CenterY)
    if ((root_y -= (height/2)) < 0)
      root_y = 0;
    else
      if (root_y > HeightOfScreen(XtScreen(Popup)) -height)
	root_y =  HeightOfScreen(XtScreen(Popup)) -height;
  
  XtVaSetValues(Popup,
		XtNx, root_x,
		XtNy, root_y,
		NULL);
}

