#include <coVRModuleSupport.h>
#include <coPanel.h>
#include <coPopupHandle.h>
#include <coFlatPanelGeometry.h>
#include <coFlatButtonGeometry.h>
#include "coRectButtonGeometry.h"
#include <coValuePoti.h>
#include <coButton.h>
#include <coLabel.h>
#include <coFrame.h>
#include "LimitsDialog.h"
#include "IMDPlugin.h"

/// Constructor
LimitsDialog::LimitsDialog(const char* title, IMDPlugin* plugin)
{
  char labelString[ATOMS_FILT_SIZE][16] = {"Sorte","X pos","Y pos","Z pos","Ekin","Epot","Nbanz"};
  float checkboxOffset[ATOMS_FILT_SIZE] = {3,6,6,6,3,6,3};
  float x[4] = {0, 20, 80, 130};
  float dy = 35;
  float y;
  float z = 1;
  float checkboxSize = 12;
  float labelSize = 17;
  float potiSize = 0.5;
  float buttonsY = 0;
  float buttonSize[2] = {70, 25};
  int i;

  myPlugin = plugin;
  panel = new coPanel(new coFlatPanelGeometry(coUIElement::BLACK));
  handle = new coPopupHandle(title);

  // Create parameter arrays:
  checkbox = new coToggleButton*[ATOMS_FILT_SIZE];
  label = new coLabel*[ATOMS_FILT_SIZE];
  potiMin = new coValuePoti*[ATOMS_FILT_SIZE];
  potiMax = new coValuePoti*[ATOMS_FILT_SIZE];
  
  // Loop thru all parameters:
  y = ATOMS_FILT_SIZE * dy;
  for (i=0; i<ATOMS_FILT_SIZE; ++i)   
  {
    // Create widgets:
    checkbox[i] = new coToggleButton(new coFlatButtonGeometry("UI/haken"),this);
    checkbox[i]->setPos(x[0], y+checkboxOffset[i], z);
    checkbox[i]->setSize(checkboxSize);
    checkbox[i]->setState(false, false);
    label[i] = new coLabel();
    label[i]->setString(labelString[i]);
    label[i]->setPos(x[1], y+labelSize/5, z);
    label[i]->setFontSize(labelSize);
    potiMin[i] = new coValuePoti("min", this, "Volume/valuepoti-bg");
    potiMin[i]->setPos(x[2], y-potiSize, z);
    potiMin[i]->setSize(potiSize);
    potiMin[i]->setMin(0);
    potiMin[i]->setMax(1);
    potiMin[i]->setValue(0);
    if (i==0 || i==6) potiMin[i]->setInteger(true);
    potiMax[i] = new coValuePoti("max", this, "Volume/valuepoti-bg");
    potiMax[i]->setPos(x[3], y-potiSize, z);
    potiMax[i]->setSize(potiSize);
    potiMax[i]->setMin(0);
    potiMax[i]->setMax(1);
    potiMax[i]->setValue(0);
    if (i==0 || i==6) potiMax[i]->setInteger(true);
    y -= dy;

    // Add widgets to panel:
    panel->addElement(checkbox[i]);
    panel->addElement(label[i]);
    panel->addElement(potiMin[i]);
    panel->addElement(potiMax[i]);
  }  

  applyDisp = new coPushButton(new coRectButtonGeometry(buttonSize[0],buttonSize[1],"IMD/display"),this);
  applyDisp->setPos(5, buttonsY);
  applyTrans = new coPushButton(new coRectButtonGeometry(buttonSize[0],buttonSize[1],"IMD/transfer"),this);
  applyTrans->setPos(85, buttonsY);
  panel->addElement(applyDisp);
  panel->addElement(applyTrans);

  panel->setScale(3); // set size of panel and contents
  panel->resize();

  coFrame *frame = new coFrame("UI/Frame");
  frame->addElement(panel);
  handle->addElement(frame);
}

/// Destructor
LimitsDialog::~LimitsDialog()
{
  int i;
  
  for (i=0; i<ATOMS_FILT_SIZE; ++i)
  {
    delete checkbox[i];
    delete label[i];
    delete potiMin[i];
    delete potiMax[i];
  }
  delete[] checkbox;
  delete[] label;
  delete[] potiMin;
  delete[] potiMax;
}

void LimitsDialog::setMin(float m)
{
}

void LimitsDialog::setMax(float m)
{
}

float LimitsDialog::getMin()
{
  return 0;
}

float LimitsDialog::getMax()
{
  return 1;
}

void LimitsDialog::setVisible(bool visible)
{
  handle->setVisible(visible);
}

void LimitsDialog::update()
{
  handle->update();
}

void LimitsDialog::potiValueChanged(float oldValue, float newValue, 
  coValuePoti* poti, int context)
{
}

void LimitsDialog::buttonEvent(coButton* button)
{
  if (button==applyDisp)
  {
    myPlugin->updateDisplay();
  }
  else if (button==applyTrans)
  {
    myPlugin->updateTransfer();
  }
}
