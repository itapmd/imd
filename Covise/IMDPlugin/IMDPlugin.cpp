// **************************************************************************
//
// Description:   IMD plugin for simulation steering with IMD (ITAP)
//
// Author:        Juergen Schulze
//
// Creation Date: 2003-01-09
//
// **************************************************************************

#include <coVRModuleSupport.h>
#include <coRowMenu.h>
#include <coButtonMenuItem.h>
#include <coCheckboxMenuItem.h>
#include <coSubMenuItem.h>
#include <coPotiMenuItem.h>
#include <Performer/pf/pfGeode.h>
#include <Performer/pf/pfDCS.h>
#include <Performer/pr/pfMaterial.h>
#include "IMDPlugin.h"
#include "AtomScenegraph.h"
#include "LimitsDialog.h"

IMDPlugin* plugin = NULL;

int coVRInit(coVRModule *m)
{
  if (plugin==NULL) plugin = new IMDPlugin(m);
  return 0;
}

// REQUIRED!
void coVRDelete(coVRModule *m)
{
  delete plugin;
}


// coVRPreFrame()
// called before each frame
void coVRPreFrame()
{
  plugin->preFrame();
}

/// Constructor
IMDPlugin::IMDPlugin(coVRModule *m)
{
  coSubMenuItem* pinboardEntry;
  const char* remotePortString;
  const char* maxAtomsString;
  int i;

  cerr << "IMDPlugin started" << endl;

  // Initialize variables:
  endian = vvSocket::VV_LITTLE_END;
  socket = NULL;
  atomScenegraph = NULL;
  coloring = 0;
  remotePort = 31050;
  remoteHost = NULL;
  sendFlagsItem = new coCheckboxMenuItem*[ATOMS_FLAG_SIZE];
  
  // Read atom constraint:
  maxAtomsString = CoviseConfig::getEntry("IMDPluginConfig.MaxAtomsRange");
  if (maxAtomsString) maxNumAtoms = atoi(maxAtomsString);
  else maxNumAtoms = 10000;
  cerr << "Atoms constraint maximum value: " << maxNumAtoms << endl;
  
  // Initialize flags:
  stockFlags.sorte  = 0;
  stockFlags.ort    = 0;
  stockFlags.impuls = 0;
  stockFlags.Ekin   = 0;
  stockFlags.Epot   = 0;
  stockFlags.nbanz  = 0;
  
  // Initialize min/max values and filters:
  for (i=0; i<ATOMS_FILT_SIZE; ++i)
  {
    minValue[i]  = 0.0f;
    maxValue[i]  = 1.0f;
  }

  // Create limits dialog:
  limits = new LimitsDialog("Filters", this); 
  // The name passed to the LimitsDialog constructor determines the 
  // section in covise.config to use for menu appearance settings.
  limits->setVisible(false);

  // Create simulation parameters menu:
  simMenu  = new coRowMenu("Simulation Parameters");
  dehnItem = new coPotiMenuItem("Dehnungsparameter", -5.0, 5.0, 0.0);
  updateSimItem = new coButtonMenuItem("Get value");
  dehnItem->setMenuListener(this);
  dehnItem->setInteger(false);
  updateSimItem->setMenuListener(this);
  simMenu->add(dehnItem);
  simMenu->add(updateSimItem);
  
  // Create color assignment menu:
  colorMenu = new coRowMenu("Color Assignment");
  colorGroup = new coCheckboxGroup();
  ekinItem  = new coCheckboxMenuItem("Ekin", false, colorGroup);
  epotItem  = new coCheckboxMenuItem("Epot", false, colorGroup);
  nbanzItem = new coCheckboxMenuItem("nbanz", false, colorGroup);
  sorteItem = new coCheckboxMenuItem("Sorte", true, colorGroup);
  sorteItem->setMenuListener(this);
  ekinItem->setMenuListener(this);
  epotItem->setMenuListener(this);
  nbanzItem->setMenuListener(this);
  colorMenu->add(sorteItem);
  colorMenu->add(ekinItem);
  colorMenu->add(epotItem);
  colorMenu->add(nbanzItem);

  // Create preferences menu:
  prefMenu = new coRowMenu("Preferences");
  animItem      = new coCheckboxMenuItem("Auto animate", false);
  animTimeItem  = new coPotiMenuItem("Update period (sec)", 0.05, 20.0, 1.0);
  sizeItem      = new coPotiMenuItem("Atom diameter", 0.1, 1.0, 0.3);
  numAtomsItem  = new coPotiMenuItem("Max num atoms", 0.0, float(maxNumAtoms), 0.0);
  verboseItem   = new coCheckboxMenuItem("Verbose mode", false);
  numAtomsItem->setInteger(true);
  animItem->setMenuListener(this);
  animTimeItem->setMenuListener(this);
  sizeItem->setMenuListener(this);
  prefMenu->add(animItem);
  prefMenu->add(animTimeItem);
  prefMenu->add(sizeItem);
  prefMenu->add(numAtomsItem);
  prefMenu->add(verboseItem);
  
  // Create send flags menu:
  flagsMenu = new coRowMenu("Send Flags");
  sendFlagsItem[0] = new coCheckboxMenuItem("Sorte", true);
  sendFlagsItem[1] = new coCheckboxMenuItem("Ort", true);
  sendFlagsItem[2] = new coCheckboxMenuItem("Impuls", true);
  sendFlagsItem[3] = new coCheckboxMenuItem("Ekin", true);
  sendFlagsItem[4] = new coCheckboxMenuItem("Epot", true);
  sendFlagsItem[5] = new coCheckboxMenuItem("nbanz", true); 
  for (i=0; i<ATOMS_FLAG_SIZE; ++i)
  {
    flagsMenu->add(sendFlagsItem[i]);
  }
  
  // Create main menu:
  imdMenu     = new coRowMenu("VisIMD");
  initItem    = new coButtonMenuItem("Get parameters");
  limitsItem  = new coButtonMenuItem("Parameter limits...");
  flagsItem   = new coSubMenuItem("Send flags");
  simItem     = new coSubMenuItem("Simulation settings");
  colorItem   = new coSubMenuItem("Color assignment");
  prefItem    = new coSubMenuItem("Preferences");
  restartItem = new coButtonMenuItem("Restart simulation");
  quitItem    = new coButtonMenuItem("Quit simulation");
  initItem->setMenuListener(this);
  limitsItem->setMenuListener(this);
  flagsItem->setMenu(flagsMenu);
  simItem->setMenu(simMenu);
  colorItem->setMenu(colorMenu);
  prefItem->setMenu(prefMenu);  
  restartItem->setMenuListener(this);
  quitItem->setMenuListener(this);
  imdMenu->add(initItem);
  imdMenu->add(limitsItem);
  imdMenu->add(flagsItem);
  imdMenu->add(simItem);
  imdMenu->add(colorItem);
  imdMenu->add(prefItem);
  imdMenu->add(restartItem);
  imdMenu->add(quitItem);

  pinboardEntry = new coSubMenuItem("IMD");
  cover->getMenu()->add(pinboardEntry);
  pinboardEntry->setMenu(imdMenu);

  // Process covise.config entries:
  remoteHost = CoviseConfig::getEntry("IMDPluginConfig.Host");
  remotePortString = CoviseConfig::getEntry("IMDPluginConfig.Port");
  if (remotePortString) remotePort = atoi(remotePortString);
  if (remoteHost)
  {
    cerr << "Connection mode: client" << endl;
    cerr << "Remote host:     " << remoteHost << endl;
  }
  else
  {
    cerr << "Connection mode: server" << endl;
  }
  cerr << "Remote port:     " << remotePort << endl;
   
  cerr << "IMDPlugin initialized" << endl;
}


/// Destructor
IMDPlugin::~IMDPlugin()
{
  delete atomScenegraph;
  delete socket;

  // Delete menu components:
  delete updateSimItem;
  delete verboseItem;
  delete animItem;
  delete dehnItem;
  delete quitItem;
  delete simItem;
  delete initItem;
  delete simMenu;
  delete limitsMenu;
  delete prefMenu;
  delete imdMenu;
  
  delete limits;
  delete[] sendFlagsItem;
}  

/// Process release events
void IMDPlugin::menuReleaseEvent(coMenuItem* item)
{
  if (item==sizeItem)
  {
    updateDisplay();
  }
  else if (item==dehnItem)
  {
    changeParams();
  }
}

/// Process VR menu events
void IMDPlugin::menuEvent(coMenuItem* item)
{
  if (item==initItem)
  {
    initConnection();
  }
  else if (item==quitItem)
  {
    quitConnection();
  }
  else if (item==updateSimItem)
  {
    getSimulationParams();
  }
  else if (item==restartItem)
  {
    restartConnection();
  }
  else if (item==limitsItem)
  {
    limits->setVisible(true);
  }
  else if (item==sorteItem)
  {
    coloring = 0;
    updateDisplay();
  }
  else if (item==ekinItem)
  {
    coloring = 1;
    updateDisplay();
  }
  else if (item==epotItem)
  {
    coloring = 2;
    updateDisplay();
  }
  else if (item==nbanzItem)
  {
    coloring = 3;
    updateDisplay();
  }
}

vvSocketIO* IMDPlugin::createSocket()
{
  vvSocketIO* sock;

  if (remoteHost==NULL) // server connection
  {
    cerr << "Opening server connection" << endl;
    sock = new vvSocketIO(remotePort, vvSocketIO::VV_TCP);
  }
  else      // client connection
  {
    cerr << "Opening client connection" << endl;
    sock = new vvSocketIO(remotePort, (char*)remoteHost, vvSocketIO::VV_TCP);
  }
  sock->set_sock_param(300.0f, 300.0f, 65535, 0);
  if (sock->init() != vvSocketIO::VV_OK)
  {
    delete sock;
    return NULL;
  }
  else return sock;
}

/// Send initialization command to simulation
void IMDPlugin::initConnection()
{
  uchar buffer[4];
  char endianString[2][7] = {"little", "big"};
  int i, buf;
  static bool firstInit = true;

  // Create socket:
  if (socket==NULL) socket = createSocket();
  if (socket==NULL)
  {
    cerr << "Cannot establish connection" << endl;
    return;
  }  
  else cerr << "Connection established" << endl;


  // Initialize connection:
  if (socket->write8(VIS_INIT) != vvSocketIO::VV_OK)
  {
    cerr << "Cannot write VIS_INIT to socket" << endl;
    return;
  }
  if (socket->read_data(buffer, 4) == vvSocketIO::VV_OK)
  {
    cerr << "Initialization data received" << endl;
    if (buffer[0]!=PROTOCOL_VERSION_MAJOR) cerr << "unexpected PROTOCOL_VERSION_MAJOR: " << int(buffer[0]) << endl;
    if (buffer[1]!=PROTOCOL_VERSION_MINOR) cerr << "unexpected PROTOCOL_VERSION_MINOR: " << int(buffer[1]) << endl;
    endian = (buffer[2]==1) ? vvSocket::VV_BIG_END : vvSocket::VV_LITTLE_END;
    cerr << "Endianness: " << endianString[buffer[2]] << endl;
    if (buffer[3]!=3) cerr << "unexpected DIM (should be 3): " << int(buffer[3]) << endl;
  }
  else
  {
    cerr << "Cannot read data from socket" << endl;
    return;
  }
  
  // Initialize atom parameters:
  if (socket->write8(VIS_INIT_ATOMS) != vvSocketIO::VV_OK) assert(0);
  for (i=0; i<ATOMS_FLAG_SIZE; ++i)
  {
    buf = socket->read32(endian);
    switch(i)
    {
      case 0: stockFlags.sorte  = buf; break;
      case 1: stockFlags.ort    = buf; break;
      case 2: stockFlags.impuls = buf; break;
      case 3: stockFlags.Ekin   = buf; break;
      case 4: stockFlags.Epot   = buf; break;
      case 5: stockFlags.nbanz  = buf; break;
      default: break;
    }
  }
   
  if (verboseItem->getState())
  {
    cerr << "Atom flags received:" << endl;
    cerr << "sorte  = " << stockFlags.sorte << endl;
    cerr << "ort    = " << stockFlags.ort << endl;
    cerr << "impuls = " << stockFlags.impuls << endl;
    cerr << "Ekin   = " << stockFlags.Ekin << endl;
    cerr << "Epot   = " << stockFlags.Epot << endl;
    cerr << "nbanz  = " << stockFlags.nbanz << endl;
    cerr << endl;
  }

  // Update send and filter flag checkboxes:
  if (stockFlags.sorte==0)  
  {
    sendFlagsItem[0]->setState(false);
    limits->checkbox[0]->setState(false);
  }
  if (stockFlags.ort==0)
  {
    sendFlagsItem[1]->setState(false);
    limits->checkbox[1]->setState(false);
    limits->checkbox[2]->setState(false);
    limits->checkbox[3]->setState(false);
  }
  if (stockFlags.impuls==0) 
  {
    sendFlagsItem[2]->setState(false);
  }
  if (stockFlags.Ekin==0)
  {
    sendFlagsItem[3]->setState(false);
    limits->checkbox[4]->setState(false);
  }
  if (stockFlags.Epot==0)
  {
    sendFlagsItem[4]->setState(false);
    limits->checkbox[5]->setState(false);
  }
  if (stockFlags.nbanz==0)
  {
    sendFlagsItem[5]->setState(false);
    limits->checkbox[6]->setState(false);
  }
  
  for (i=0; i<ATOMS_FILT_SIZE; ++i)
  {
    minValue[i] = socket->readFloat(endian);
  }
  for (i=0; i<ATOMS_FILT_SIZE; ++i)
  {
    maxValue[i] = socket->readFloat(endian);
  }
  if (verboseItem->getState()) 
  {
    cerr << "Value ranges received:" << endl;
    cerr << "sorte: " << minValue[0] << " to " << maxValue[0] << endl;
    cerr << "x:     " << minValue[1] << " to " << maxValue[1] << endl;
    cerr << "y:     " << minValue[2] << " to " << maxValue[2] << endl;
    cerr << "z:     " << minValue[3] << " to " << maxValue[3] << endl;
    cerr << "Ekin:  " << minValue[4] << " to " << maxValue[4] << endl;
    cerr << "Epot:  " << minValue[5] << " to " << maxValue[5] << endl;
    cerr << "nbanz: " << minValue[6] << " to " << maxValue[6] << endl;
    cerr << endl;
  }
  
  // Update min and max filter values:
  for (i=0; i<ATOMS_FILT_SIZE; ++i)
  {
    if (i==0 || i==6)  // treat sorte and nbanz differently: prevent min=max
    {
      limits->potiMin[i]->setMin(minValue[i]);
      limits->potiMax[i]->setMin(minValue[i]);
      // Integer poti: prevent min=max
      if (minValue[i] == maxValue[i]) 
      {
        limits->potiMin[i]->setMax(minValue[i] + 1);
        limits->potiMax[i]->setMax(minValue[i] + 1);
      }
      else
      {
        limits->potiMin[i]->setMax(maxValue[i]);
        limits->potiMax[i]->setMax(maxValue[i]);
      }
      if (firstInit)
      {
        limits->potiMin[i]->setValue(minValue[i]);
        limits->potiMax[i]->setValue(maxValue[i]);
      }
    }
    else
    {
      float min = spaceBottom(minValue[i]);
      float max = spaceTop(maxValue[i]);
      limits->potiMin[i]->setMin(min);
      limits->potiMin[i]->setMax(max);
      limits->potiMax[i]->setMin(min);
      limits->potiMax[i]->setMax(max);
      if (firstInit)
      {
        limits->potiMin[i]->setValue(minValue[i]);
        limits->potiMax[i]->setValue(maxValue[i]);
      }
    }
  }

  firstInit = false;
}

/// Send new parameter set to simulation
void IMDPlugin::changeParams()
{
  int changeFlag = 1;
  int step;
  float newValue;

  if (!socket) 
  {
    cerr << "no connection" << endl;
    return;
  }

  // Send request:
  cerr << "Sending new deform paramter: " << dehnItem->getValue() << endl;
  socket->write8(VIS_CHANGE_PARAMS);
  socket->write32(VIS_PARAM_DEFORM, endian);
  socket->write32(changeFlag, endian);
  socket->writeFloat(dehnItem->getValue(), endian);

  // Receive reply:
  step = socket->read32(endian);
  cerr << "Time step: " << step << endl;
  
  newValue = socket->readFloat(endian);
  cerr << "New value: " << newValue << endl;
}

/// Update atoms by requesting them from IMD.
void IMDPlugin::updateTransfer()
{
  atoms_flag_t sendFlags;   // flags determining which atom components should be sent
  atoms_flag_t filterFlags; // flags determining which atom components should be filtered
  IMDAtom* newAtom;         // atom currently being read 
  int numAtoms;             // number of atoms that are going to be sent in a row, 0 terminates time step
  int floatsPerAtom;        // number of floating point values transferred per atom
  int timestep;             // current time step
  int i,j;

  if (verboseItem->getState()) cerr << "IMDPlugin::updateTransfer()" << endl;
  if (!socket) 
  { 
    cerr << "no connection" << endl;
    return;
  }

  // Read send flags from menu:
  sendFlags.sorte  = sendFlagsItem[0]->getState();
  sendFlags.ort    = sendFlagsItem[1]->getState();
  sendFlags.impuls = sendFlagsItem[2]->getState();
  sendFlags.Ekin   = sendFlagsItem[3]->getState();
  sendFlags.Epot   = sendFlagsItem[4]->getState();
  sendFlags.nbanz  = sendFlagsItem[5]->getState();
  
  // Read filter flags from menu:
  filterFlags.sorte  = limits->checkbox[0]->getState();
  if (limits->checkbox[1]->getState() || 
      limits->checkbox[2]->getState() || 
      limits->checkbox[3]->getState())
  {
    filterFlags.ort = true;
  }
  else filterFlags.ort = false;
  filterFlags.impuls = 0;
  filterFlags.Ekin   = limits->checkbox[4]->getState();
  filterFlags.Epot   = limits->checkbox[5]->getState();
  filterFlags.nbanz  = limits->checkbox[6]->getState();

  // Verify flags:
  if (stockFlags.sorte==0)  sendFlags.sorte  = filterFlags.sorte  = 0;
  if (stockFlags.ort==0)    sendFlags.ort    = filterFlags.ort    = 0;
  if (stockFlags.impuls==0) sendFlags.impuls = filterFlags.impuls = 0;
  if (stockFlags.Ekin==0)   sendFlags.Ekin   = filterFlags.Ekin   = 0;
  if (stockFlags.Epot==0)   sendFlags.Epot   = filterFlags.Epot   = 0;
  if (stockFlags.nbanz==0)  sendFlags.nbanz  = filterFlags.nbanz  = 0;

  // Send command:
  cerr << "Sending VIS_WRITE_ATOMS" << endl;
  socket->write8(VIS_WRITE_ATOMS);
  cerr << "VIS_WRITE_ATOMS sent" << endl;

  // Send send flags:
  socket->write32(sendFlags.sorte,  endian);
  socket->write32(sendFlags.ort,    endian);
  socket->write32(sendFlags.impuls, endian);
  socket->write32(sendFlags.Ekin,   endian);
  socket->write32(sendFlags.Epot,   endian);
  socket->write32(sendFlags.nbanz,  endian);

  if (verboseItem->getState()) 
  {
    cerr << "Send flags sent:" << endl;
    cerr << "sorte:  " << sendFlags.sorte << endl;
    cerr << "ort:    " << sendFlags.ort << endl;
    cerr << "impuls: " << sendFlags.impuls << endl;
    cerr << "Ekin:   " << sendFlags.Ekin << endl;
    cerr << "Epot:   " << sendFlags.Epot << endl;
    cerr << "nbanz:  " << sendFlags.nbanz << endl;
    cerr << endl;
  }

  // Send filter flags:
  socket->write32(filterFlags.sorte,  endian);
  socket->write32(filterFlags.ort,    endian);
  socket->write32(filterFlags.impuls, endian);
  socket->write32(filterFlags.Ekin,   endian);
  socket->write32(filterFlags.Epot,   endian);
  socket->write32(filterFlags.nbanz,  endian);

  if (verboseItem->getState()) 
  {
    cerr << "Filter flags sent:" << endl;
    cerr << "sorte:  " << filterFlags.sorte << endl;
    cerr << "ort:    " << filterFlags.ort << endl;
    cerr << "impuls: " << filterFlags.impuls << endl;
    cerr << "Ekin:   " << filterFlags.Ekin << endl;
    cerr << "Epot:   " << filterFlags.Epot << endl;
    cerr << "nbanz:  " << filterFlags.nbanz << endl;
    cerr << endl;
  }

  // Send min filter values:
  for (i=0; i<ATOMS_FILT_SIZE; ++i)
  {
    socket->writeFloat(limits->potiMin[i]->getValue(), endian);
  }

  // Send max filter values:
  for (i=0; i<ATOMS_FILT_SIZE; ++i)
  {
    socket->writeFloat(limits->potiMax[i]->getValue(), endian);
  }

  if (verboseItem->getState()) 
  {
    cerr << "Filter values sent:" << endl;
    cerr << "sorte: " << limits->potiMin[0]->getValue() << " to " << limits->potiMax[0]->getValue() << endl;
    cerr << "x:     " << limits->potiMin[1]->getValue() << " to " << limits->potiMax[1]->getValue() << endl;
    cerr << "y:     " << limits->potiMin[2]->getValue() << " to " << limits->potiMax[2]->getValue() << endl;
    cerr << "z:     " << limits->potiMin[3]->getValue() << " to " << limits->potiMax[3]->getValue() << endl;
    cerr << "Ekin:  " << limits->potiMin[4]->getValue() << " to " << limits->potiMax[4]->getValue() << endl;
    cerr << "Epot:  " << limits->potiMin[5]->getValue() << " to " << limits->potiMax[5]->getValue() << endl;
    cerr << "nbanz: " << limits->potiMin[6]->getValue() << " to " << limits->potiMax[6]->getValue() << endl;
    cerr << endl;
  }

  // Receive time step and number of atoms:
  timestep = socket->read32(endian);
  floatsPerAtom = socket->read32(endian);
  cerr << "Current time step: " << timestep << endl;
  if (verboseItem->getState()) cerr << "Floating point values per atom: " << floatsPerAtom << endl;

  // Receive atom configuration:
  atoms.clear();
  for (;;)    // stop only on break
  {
    numAtoms = socket->read32(endian);
    if (verboseItem->getState()) cerr << "Number of atoms to receive: " << numAtoms << endl;
    if (numAtoms<=0) break;   // done

    // Read atom parameters:
    for (i=0; i<numAtoms; ++i)
    { 
      newAtom = new IMDAtom();

      if (sendFlags.sorte) 
      {
        newAtom->sorte = socket->readFloat(endian);
      }
      if (sendFlags.ort)
      {
        for (j=0; j<3; ++j)
        {
          newAtom->pos[j] = socket->readFloat(endian);
        }
      }
      if (sendFlags.impuls) 
      {
        for (j=0; j<3; ++j)
        {
          newAtom->impuls[j] = socket->readFloat(endian);
        }
      }
      if (sendFlags.Ekin) 
      {
        newAtom->ekin = socket->readFloat(endian);
      }
      if (sendFlags.Epot) 
      {
        newAtom->epot = socket->readFloat(endian);
      }
      if (sendFlags.nbanz)
      {
        newAtom->nbanz = socket->readFloat(endian);
      }
      atoms.append(newAtom);
    }
  }
  cerr << "Number of atoms received: " << atoms.count() << endl;
//  if (verboseItem->getState()) listAtoms();
  updateDisplay();
}

/// Send quit command to simulation
void IMDPlugin::quitConnection()
{
  uchar buffer = VIS_WRITE_QUIT;
 
  if (!socket)
  {
    cerr << "no connection" << endl;
    return;
  }
  socket->write8(buffer);
  delete socket;
  socket = NULL;
}

/// Send restart command to simulation
void IMDPlugin::restartConnection()
{
  uchar buffer = VIS_RESTART;
 
  if (!socket)
  {
    cerr << "no connection" << endl;
    return;
  }
  socket->write8(buffer);
}

/// Display all stored atoms on cerr.
void IMDPlugin::listAtoms()
{
  int numAtoms;
  int i;

  numAtoms = atoms.count();
  cerr << "Number of atoms in list: " << numAtoms << endl;
  for (i=0; i<numAtoms; ++i)
  {
    cerr << "Atom " << i << ": ";
    atoms[i]->print();
    cerr << endl;
  }
}

/// Update graphical atom representation.
void IMDPlugin::updateDisplay()
{
  int numAtoms;
  float color, red, green, blue;
  int i;

  if (verboseItem->getState()) cerr << "IMDPlugin::updateDisplay()" << endl;
  numAtoms = atoms.count();
  if (atomScenegraph) delete atomScenegraph;
  atomScenegraph = new AtomScenegraph();
  
  // Constrain number of atoms to selected value:
  if (int(numAtomsItem->getValue()) > 0)
  {
    if (numAtoms > int(numAtomsItem->getValue())) 
    {
      numAtoms = int(numAtomsItem->getValue());
      cerr << "Number of atoms constrained to: " << numAtoms << endl;
    }
  }
  
  // Actually attach atoms to scenegraph:
  for(i=0; i<numAtoms; ++i)
  {
    if (!isFiltered(atoms[i]))
    {
      switch(coloring)
      {
        default:
        case 0:  
          if (minValue[0]==maxValue[0]) color = 0.5f;
          else color = (atoms[i]->sorte - minValue[0]) / 
                       (maxValue[0] - minValue[0]);
          break;
        case 1:
          if (minValue[4]==maxValue[4]) color = 0.5f;
          else color = (atoms[i]->ekin - minValue[4]) / 
                       (maxValue[4] - minValue[4]);
          break;
        case 2:
          if (minValue[5]==maxValue[5]) color = 0.5f;
          else color = (atoms[i]->epot - minValue[5]) / 
                       (maxValue[5] - minValue[5]);
          break;
        case 3:
          if (minValue[6]==maxValue[6]) color = 0.5f;
          else color = (atoms[i]->nbanz - minValue[6]) / 
                       (maxValue[6] - minValue[6]);
          break;
      }
      color = coClamp(color, 0.0f, 1.0f);
      createRGBColor(color, &red, &green, &blue);
      atomScenegraph->addAtom(atoms[i]->pos[0], atoms[i]->pos[1], 
        atoms[i]->pos[2], sizeItem->getValue(), red, green, blue);
    }
  }
}

/// Called before each frame
void IMDPlugin::preFrame()
{
  static double prevTime = pfGetTime();
  
  if (animItem->getState())
  {
    if ((pfGetTime() - prevTime) >= animTimeItem->getValue())
    {
      updateTransfer();
      prevTime = pfGetTime();
    }
  }

  limits->update();
}

/** Create an RGB color from a scalar color.
  Use color gradient: blue, green, yellow, red
  @param color scalar color [0..1]
  @param r,g,b RGB color [0..1]
*/
void IMDPlugin::createRGBColor(float color, float* r, float* g, float* b)
{
  // blue   = 0,0,1
  // green  = 0,1,0
  // yellow = 1,1,0
  // red    = 1,0,0
  
  // red:
  if (color < 0.33) *r = 0.0;
  else if (color > 0.67) *r = 1.0;
  else *r = (color - 0.33) * 3.0;

  // green:  
  if (color < 0.33) *g = color * 3.0;
  else if (color > 0.67) *g = 1.0 - ((color - 0.67) * 3.0);
  else *g = 1.0;

  // blue:
  *b = (color > 0.33) ? 0.0 : (1.0 - color * 3.0);
}

/// Returns a value that is 10% below the passed value.
float IMDPlugin::spaceBottom(float val)
{
  return(val - 0.1 * fabs(val));
}

/// Returns a value that is 10% above the passed value.
float IMDPlugin::spaceTop(float val)
{
  return(val + 0.1 * fabs(val));
}

/** Apply the current filter settings to the atom.
  @return true if atom is to be filtered, false if it is to be displayed
*/
bool IMDPlugin::isFiltered(IMDAtom* atom)
{
  bool filter = false;
  
  if (limits->checkbox[0]->getState() && 
      (atom->sorte < limits->potiMin[0]->getValue() || 
       atom->sorte > limits->potiMax[0]->getValue())) 
    filter = true;
  if (limits->checkbox[1]->getState() && 
      (atom->pos[0] < limits->potiMin[1]->getValue() || 
       atom->pos[0] > limits->potiMax[1]->getValue())) 
    filter = true;
  if (limits->checkbox[2]->getState() && 
      (atom->pos[1] < limits->potiMin[2]->getValue() || 
       atom->pos[1] > limits->potiMax[2]->getValue())) 
    filter = true;
  if (limits->checkbox[3]->getState() && 
      (atom->pos[2] < limits->potiMin[3]->getValue() || 
       atom->pos[2] > limits->potiMax[3]->getValue())) 
    filter = true;
  if (limits->checkbox[4]->getState() && 
      (atom->ekin < limits->potiMin[4]->getValue() || 
       atom->ekin > limits->potiMax[4]->getValue())) 
    filter = true;
  if (limits->checkbox[5]->getState() && 
      (atom->epot < limits->potiMin[5]->getValue() || 
       atom->epot > limits->potiMax[5]->getValue())) 
    filter = true;
  if (limits->checkbox[6]->getState() && 
      (atom->nbanz < limits->potiMin[6]->getValue() || 
       atom->nbanz > limits->potiMax[6]->getValue())) 
    filter = true;
    
  return filter;
}

/// Receive simulation parameters.
void IMDPlugin::getSimulationParams()
{
  int changeFlag = 0;
  int step;
  float dehnValue;

  if (!socket) 
  {
    cerr << "no connection" << endl;
    return;
  }

  // Send request:
  cerr << "Requesting deform paramter" << endl;
  socket->write8(VIS_CHANGE_PARAMS);
  socket->write32(VIS_PARAM_DEFORM, endian);
  socket->write32(changeFlag, endian);

  // Receive reply:
  step = socket->read32(endian);
  cerr << "Time step: " << step << endl;
  
  dehnValue = socket->readFloat(endian);
  cerr << "Dehnungsparameter: " << dehnValue << endl;
  if (dehnValue < dehnItem->getMin()) dehnItem->setMin(dehnValue);
  if (dehnValue > dehnItem->getMax()) dehnItem->setMax(dehnValue);
  dehnItem->setValue(dehnValue);
}

/// Constructor
IMDAtom::IMDAtom()
{
  sorte = 0.0f;
  pos[0] = pos[1] = pos[2] = 0.0f;
  impuls[0] = impuls[1] = impuls[2] = 0.0f;
  ekin = 0.0f;
  epot = 0.0f;
  nbanz = 0;
}

/// Print atom parameters on cerr.
void IMDAtom::print()
{
  cerr.setf(ios::fixed, ios::floatfield);
  cerr.precision(5);
  cerr << sorte << " " << pos[0] << " " << pos[1] << " " << pos[2] << " " << 
    impuls[0] << " " << impuls[1] << " " << impuls[2] << " " << 
    ekin << " " << epot << " " << nbanz;
}

// EOF
