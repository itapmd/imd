#ifndef LIMITS_DIALOG_H
#define LIMITS_DIALOG_H

class coPanel;
class coPushButton;
class coPopupHandle;
class coToggleButton;
class coValuePoti;
class IMDPlugin;

/** Dialog window to set filter values.
*/
class LimitsDialog: public coButtonActor, public coValuePotiActor
{
  public:
    coPanel* panel;
    coPopupHandle* handle;
    coToggleButton** checkbox;
    coValuePoti** potiMin;
    coValuePoti** potiMax;

    LimitsDialog(const char*, IMDPlugin*);
    ~LimitsDialog();
    void  setMin(float m);
    void  setMax(float m);
    float getMin();
    float getMax();
    void  setVisible(bool);
    void  update();

  protected:
    IMDPlugin* myPlugin;
    coLabel** label;
    coPushButton* applyDisp;
    coPushButton* applyTrans;
    
    void buttonEvent(coButton*);
    void potiValueChanged(float, float, coValuePoti*, int =-1);
};

#endif
