#ifndef _INV_PIXMAP_BUTTON_
#define _INV_PIXMAP_BUTTON_


/* $Id$ */

/* $Log$
 * Revision 1.1  1994/04/12  13:39:31  zrfu0125
 * Initial revision
 * */




//**************************************************************************
//
// * Description    : Pixmap button class 
//                   
// * Class(es)      : InvPixmapButton
//
// * inherited from : none
//
// * Author  : Dirk Rantzau
//
// * History : 29.03.94 V 1.0
//
//**************************************************************************


#include <X11/Intrinsic.h>
#include <Inventor/SbBasic.h>

class InvPixmapButton {
  public:
    InvPixmapButton(Widget parent, SbBool selectable);
    ~InvPixmapButton();
    
    // return the motif push button
    Widget	getWidget()	    { return widget; }
    
    // set the icon to use for the pixmap
    void	setIcon(char *icon, int width, int height);
    
    // Highlight the pixmap to show it it selected
    void	select(SbBool onOrOff);
    SbBool	isSelected()	    { return selectFlag; }
    
  private:
    Widget	parent, widget;
    SbBool	selectFlag, selectable;
    Pixmap	normalPixmap, armPixmap, selectPixmap;
};

#endif // _INV_PIXMAP_BUTTON_
