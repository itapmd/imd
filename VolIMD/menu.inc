/* program name ..............  VolRend  
   file name .................  menu.inc
   author ....................  Roland Niemeier
   version ...................  1.0beta
   date of last change .......  05/22/97
   description : interface for direct volume visualization of scalar fields. 
                 The application development is intended for 
                 the interactive visualization of volume data sets. 
                 It is based on Open Inventor, volpack, xforms.
   purpose of menu.inc: defines layout and callback functions of menu bar
*/

enum MenuEntries {
  M_FILE = 0,
  M_FILE_OPEN_CLASSIFIED,
  M_FILE_OPEN,
  M_FILE_BUILD_RAW_VOLUME,
  M_FILE_BUILD_OCTREE,
  M_FILE_SAVE_IMAGE,
  M_FILE_SELECT_SEQUENCE,
  M_FILE_LOAD_SEQUENCE,
  M_FILE_UNCLASSIFIED_SEQUENCE,
  M_FILE_QUIT,
  
  M_MODIFY,
  M_MODIFY_MATERIAL_PROPERTIES,
  M_MODIFY_MATERIAL_NUMBERS,
  M_MODIFY_LIGHT_PROPERTIES,
  M_MODIFY_LIGHT_NUMBERS,
  M_MODIFY_CLASSIFICATION,
  M_MODIFY_HISTOGRAM_BOUNDARIES,
  M_MODIFY_ROTATION_INCREMENT,
  M_MODIFY_SCALING,
  M_MODIFY_DEPTH_CUEING,
  M_MODIFY_MINIMUM_OPACITY,
  M_MODIFY_MAXIMUM_OPACITY,

  M_SWITCH,
  M_SWITCH_SHADOW,
  M_SWITCH_DEPTH_CUE,
  M_SWITCH_LIGHT_BOTH,
  M_SWITCH_REVERSE_SURFACE,
  M_SWITCH_SAVE_SEQUENCE,
  M_SWITCH_COLOR_MODE,
  M_SWITCH_USE_OCTREE,
  M_SWITCH_EDITOR,
  M_SWITCH_FILENAME,
  M_SWITCH_TWOD,
  M_SWITCH_HISTOGRAM,
  
  M_ROTATE,
  M_ROTATE_AROUND_X,
  M_ROTATE_AROUND_Y,
  M_ROTATE_AROUND_Z,

  M_SIMULATION,
  M_SIM_GET_DISTRIBUTION,
  M_SIM_GET_PICTURE,
  M_SIM_TOGGLE_LIMITS,
  M_SIM_TOGGLE_STORAGE,
  M_SIM_DRAW_HISTO,
  M_SIM_CHANGE_PORT,
  M_SIM_PORT_MACHINE,
  M_SIM_BYTE_SWAP,
  M_SIM_CHANGE_NAME,
  M_SIM_TERMINATE,
 
  M_ITEMS
};

// ===========================================================================================
//
//    Pulldown menus:
//
//     Label, ResourceName, Id, Typ, Callback,       Accelerator, AcceleratorLabel, Userdata
//
//     char*, char*       , int,int, XtCallbackProc, char*,       char*,            int
//
// ===========================================================================================

static MenuButtonTyp FileMenu[] =
{
  {"Open classified file ...",	"Open classified",	M_FILE_OPEN_CLASSIFIED,	M_PUSH_BUTTON,	FileCB,	0 },
  {"Open unclassified file ...",	"Open",	M_FILE_OPEN,	M_PUSH_BUTTON,	FileCB,	1 },
//  {"Build unclassified volume ...",	"Build raw volume",	M_FILE_BUILD_RAW_VOLUME,	M_PUSH_BUTTON,	FileCB,	2 },
  {"Build octree for fast classification ...",	"Build octree",	M_FILE_BUILD_OCTREE,	M_PUSH_BUTTON,	FileCB,	3 },
  {"Save image ...",	"Save",	M_FILE_SAVE_IMAGE,	M_PUSH_BUTTON,	FileCB,	4 },
  {"Select image sequence name ...",	"Save",	M_FILE_SELECT_SEQUENCE,	M_PUSH_BUTTON,	FileCB,	5 },
  {"Open classified volume sequence ...",	"Load",	M_FILE_LOAD_SEQUENCE,	M_PUSH_BUTTON,	FileCB,	6 },
  {"Open unclassified volume sequence ...",	"Load",	M_FILE_UNCLASSIFIED_SEQUENCE,	M_PUSH_BUTTON,	FileCB,	7 },
  {"Quit",	"Quit",		M_FILE_QUIT,	M_PUSH_BUTTON,	FileCB,	7,"Alt <Key> q","Alt-q" },
};

static MenuButtonTyp ModifyMenu[] =
{
  {"material properties ...",	"M1",	 M_MODIFY_MATERIAL_PROPERTIES,	 M_PUSH_BUTTON, ModifyCB,  1 },	
  //{"material numbers...",	"M2",	 M_MODIFY_MATERIAL_NUMBERS,	 M_PUSH_BUTTON, ModifyCB,  2 },	
  {"light properties ...",	"M3",	 M_MODIFY_LIGHT_PROPERTIES,	 M_PUSH_BUTTON, ModifyCB, 3 },	
  {"light numbers ...",	"M4",	 M_MODIFY_LIGHT_NUMBERS,	 M_PUSH_BUTTON, ModifyCB, 4 },	
  {"classification ...","M5",	 M_MODIFY_CLASSIFICATION,M_PUSH_BUTTON, ModifyCB, 5 },
  {"histogram values ...",	"M6",	 M_MODIFY_HISTOGRAM_BOUNDARIES,	 M_PUSH_BUTTON, ModifyCB,  6 },	
  {"rotation values ...",	"M6",	 M_MODIFY_ROTATION_INCREMENT,	 M_PUSH_BUTTON, ModifyCB,  7 },	
  {"scaling ...",	"M7",	 M_MODIFY_SCALING,	 M_PUSH_BUTTON, ModifyCB,  8 },	
  {"depth cueing values ...",	"M8",	 M_MODIFY_DEPTH_CUEING,	 M_PUSH_BUTTON, ModifyCB,  9 },	
  {"minimum opacity ...",	"M9",	 M_MODIFY_MINIMUM_OPACITY,	 M_PUSH_BUTTON, ModifyCB,  10 },	
  {"maximum opacity ...",	"M10",	 M_MODIFY_MAXIMUM_OPACITY,	 M_PUSH_BUTTON, ModifyCB,  11 },	
};

static MenuButtonTyp SwitchMenu[] =
{
//  {"shadow",		"SHADOW",	M_SWITCH_SHADOW, M_TOGGLE_BUTTON,	 SwitchCB, 0},
  {"depth cueing",	"DEPTH_CUE",	M_SWITCH_DEPTH_CUE, M_TOGGLE_BUTTON,	 SwitchCB, 0},
  {"light both sides",	"LIGHT_BOTH",	M_SWITCH_LIGHT_BOTH, M_TOGGLE_BUTTON,	 SwitchCB, 0},
  {"save sequence",   "SAVE_SEQUENCE",   M_SWITCH_SAVE_SEQUENCE, M_TOGGLE_BUTTON,SwitchCB, 0},
  {"RGB Colors", "COLOR_MODE", M_SWITCH_COLOR_MODE, M_TOGGLE_BUTTON,SwitchCB, 0},
  {"use octree", "Use the octree", M_SWITCH_USE_OCTREE, M_TOGGLE_BUTTON,SwitchCB, 0},
  {"Inventor light editor", "Select Editor", M_SWITCH_EDITOR, M_TOGGLE_BUTTON,SwitchCB, 0},
  {"Kinetic Energy", "Select Filename", M_SWITCH_FILENAME, M_TOGGLE_BUTTON,SwitchCB, 0},
  {"2D Simulation", "Select 2d", M_SWITCH_TWOD, M_TOGGLE_BUTTON,SwitchCB, 0},
  {"Histogram", "histo", M_SWITCH_HISTOGRAM, M_TOGGLE_BUTTON,SwitchCB, 0},
};

static MenuButtonTyp RotateMenu[] =
{
  {"Around x-axis", "RotX", M_ROTATE_AROUND_X, M_PUSH_BUTTON, RotateCB, 1 },	
  {"Around y-axis", "RotY", M_ROTATE_AROUND_Y, M_PUSH_BUTTON, RotateCB, 2 },	
  {"Around z-axis", "RotZ", M_ROTATE_AROUND_Z, M_PUSH_BUTTON, RotateCB, 3 },	
};

static MenuButtonTyp SimulationMenu[] =
{
  {"get distribution", "getDist", M_SIM_GET_DISTRIBUTION, M_PUSH_BUTTON, SimulationCB, 1 },	
  {"get picture", "getPic", M_SIM_GET_PICTURE, M_PUSH_BUTTON, SimulationCB, 2 },	
  {"change port address", "port address", M_SIM_CHANGE_PORT, M_PUSH_BUTTON, SimulationCB, 3 },	
// FG - the following three are not needed for IMD
//  {"change server name", "server name", M_SIM_CHANGE_NAME, M_PUSH_BUTTON, SimulationCB, 4 },	
//  {"port on this machine", "toggleServer", M_SIM_PORT_MACHINE, M_TOGGLE_BUTTON, SimulationCB, 5 },	
//  {"endian byte swap", "toggleServer", M_SIM_BYTE_SWAP, M_TOGGLE_BUTTON, SimulationCB, 6 },	
  {"use previous limits", "toggleStore", M_SIM_TOGGLE_LIMITS, M_TOGGLE_BUTTON, SimulationCB, 7 },	
  {"store results", "toggleStore", M_SIM_TOGGLE_STORAGE, M_TOGGLE_BUTTON, SimulationCB, 8 },	
  {"terminate", "terminate sim", M_SIM_TERMINATE, M_PUSH_BUTTON, SimulationCB, 9 },	
};

// ===========================================================================================
//
//    Main Menu:
//
//     Label, Id,  PulldownMenu,   NumberOfEntries, Callback
//
//     char*, int, MenuButtonTyp*, int,             XtCallbackProc
//
// ===========================================================================================

static MenuTyp Menu[] =
{
  {"File",        M_FILE,   FileMenu,    XtNumber(FileMenu)   },
  {"Modify",  M_MODIFY,  ModifyMenu,   XtNumber(ModifyMenu)  },
  {"Switch",  M_SWITCH,   SwitchMenu,    XtNumber(SwitchMenu)   },
  {"Rotate",  M_ROTATE,   RotateMenu,    XtNumber(RotateMenu)   },
  {"Simulation",  M_SIMULATION,   SimulationMenu,    XtNumber(SimulationMenu)   },
};

