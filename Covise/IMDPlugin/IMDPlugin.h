// **************************************************************************
//
// Description:   IMD plugin for simulation steering with IMD (ITAP)
//
// Author:        Juergen Schulze
//
// Creation Date: 2003-01-09
//
// **************************************************************************

#ifndef IMD_PLUGIN_H
#define IMD_PLUGIN_H

#include "coMenu.h"
#include <vvsocketio.h>
#include <vvtoolshed.h>
#include <vvarray.h>
#include "AtomScenegraph.h"
#include "LimitsDialog.h"

#define PROTOCOL_VERSION_MAJOR  0
#define PROTOCOL_VERSION_MINOR  1

#define VIS_INIT               10
#define VIS_INIT_ATOMS         15
#define VIS_WRITE_ATOMS        20
#define VIS_WRITE_DISTRIB      30
#define VIS_CHANGE_PARAMS      40
#define VIS_RESTART            50
#define VIS_QUIT               99 
#define VIS_WRITE_QUIT        100
#define VIS_PARAM_DEFORM        1

#define ATOMS_FLAG_SIZE 6
#define ATOMS_FILT_SIZE 7

class coButtonMenuItem;
class coSubMenuItem;
class coPotiMenuItem;
class coCheckboxMenuItem;
class coCheckboxGroup;

/** Data for one atom.
*/
class IMDAtom
{
  public:
    float sorte;
    float pos[3];     ///< x/y/z position
    float impuls[3];  ///< x/y/z impulse
    float ekin;       ///< speed
    float epot;       ///< potential energy
    float nbanz;      ///< number of neighbors

    IMDAtom();
    void print();
};

/** Plugin to visualize data from the IMD molecule simulation software which
  was created at ITAP (University of Stuttgart).
  To use you will have to activate the plugin and start the IMD program
  in interactive mode on some remote machine. Then you can establish the
  connection from the IMD menu in VR.
  The following covise.config entries are supported:
  <pre>
  # TCP connection parameters:
  IMDPluginConfig
  {
    Host vision.rus.uni-stuttgart.de
    Port 31050
    AtomTriangles 20
    MaxAtomsRange 10000
  }
  
  # Default IMDPlugin filter menu position:
  Filters
  {
    MenuPosition 500.0 500.0 -400.0
    MenuOrientation -45.0 0.0 0.0
    MenuSize 1.0
  }
  </pre>
*/
class IMDPlugin : public coMenuListener
{
  private:
    typedef struct 
    {
      int sorte;
      int ort;
      int impuls;
      int Ekin;
      int Epot;
      int nbanz;
    } atoms_flag_t;

    typedef struct 
    {
      float sorte;
      float x;
      float y;
      float z;
      float Ekin;
      float Epot;
      float nbanz;
    } atoms_filt_t;

    coMenu* imdMenu;
    coMenu* limitsMenu;
    coMenu* simMenu;
    coMenu* colorMenu;
    coMenu* prefMenu;
    coMenu* flagsMenu;
    
    // Main menu items:
    coButtonMenuItem* restartItem;
    coButtonMenuItem* initItem;
    coButtonMenuItem* quitItem;
    coButtonMenuItem* limitsItem;
    coSubMenuItem* flagsItem;
    coSubMenuItem* simItem;
    coSubMenuItem* colorItem;
    coSubMenuItem* prefItem;
    
    // Simulation menu items:
    coPotiMenuItem* dehnItem;
    coButtonMenuItem* updateSimItem;
    
    // Color menu items:
    coCheckboxGroup* colorGroup;
    coCheckboxMenuItem* sorteItem;
    coCheckboxMenuItem* ekinItem;
    coCheckboxMenuItem* epotItem;
    coCheckboxMenuItem* nbanzItem;
    
    // Preferences menu items:
    coPotiMenuItem* sizeItem;
    coCheckboxMenuItem* animItem;
    coPotiMenuItem* animTimeItem;
    coPotiMenuItem* numAtomsItem;
    coCheckboxMenuItem* verboseItem;
    
    // Send flags items:
    coCheckboxMenuItem** sendFlagsItem;
        
    vvSocketIO* socket;
    AtomScenegraph* atomScenegraph;
    float minValue[ATOMS_FILT_SIZE];
    float maxValue[ATOMS_FILT_SIZE];

    vvSocket::EndianType endian;        ///< endianness for data transferred over network
    vvArray<IMDAtom*> atoms;            ///< atom data
    atoms_flag_t stockFlags;            ///< flags determining which atom components can be transferred (are in stock)
    const char* remoteHost;             ///< remote host used for TCP connection
    int   remotePort;                   ///< port used for TCP connection
    int   coloring;                     ///< 0=sorte, 1=ekin, 2=epot, 3=nbanz
    LimitsDialog* limits;               ///< dialog window for parameter limits
    int maxNumAtoms;                    ///< maximum number of atoms to constrain display to
  
    void menuEvent(coMenuItem*);
    void menuReleaseEvent(coMenuItem*);
    vvSocketIO* createSocket();
    void initConnection();
    void restartConnection();
    void quitConnection();
    void listAtoms();
    void changeParams();
    void createRGBColor(float, float*, float*, float*);
    float spaceBottom(float);
    float spaceTop(float);
    bool isFiltered(IMDAtom*);
    void getSimulationParams();

  public:
    IMDPlugin(coVRModule *m);
    ~IMDPlugin();
    void preFrame();
    void updateDisplay();
    void updateTransfer();
};

#endif

// EOF
