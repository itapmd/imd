/*****************************************************************************
*
* read in parameter files (tag based)                               MH 260298
* 
* $Revision$
* $Date$
*
******************************************************************************/

#include "imd.h"

/* the following are needed for gettimeofday */
#include <sys/time.h>
#include <unistd.h>

#if defined(__GNUC__) && defined(__STRICT_ANSI__)
extern char *strdup(char *);
#endif

#ifdef __WATCOMC__
#define strcasecmp strcmpi
#endif

/* To do yet: improve checking of bad input files ! */
/* e.g. for forgotten starttemp/endtemp              */
/* clean up prototypes.h 
                           */
typedef enum ParamType {
  PARAM_STR, PARAM_STRPTR,
  PARAM_INT, PARAM_INT_COPY,
  PARAM_INTEGER, PARAM_INTEGER_COPY,
  PARAM_REAL, PARAM_REAL_COPY
} PARAMTYPE;

int curline; /* number of current line */

/*****************************************************************************
*
* Parameter aus Zeile auslesen / get parameter from line
*
*****************************************************************************/

/* Parameter:
   param_name ... Parametername (fuer Fehlermeldungen)
   param ........ Adresse der Variable fuer den Parameter
   ptype ........ Parametertyp
                  folgende Werte sind zulaessig:
                  PARAM_STR : String, deklariert als char[]
                  PARAM_STRPTR : String, deklariert als Zeiger auf char*
                  PARAM_INT : Integer-Wert(e)
                  PARAM_INT_COPY : Integer-Wert(e), kopierend
                  PARAM_REAL : Real-Wert(e)
                  PARAM_REAL_COPY : Real-Wert(e), kopierend
                  
   pnum_min ..... Minimale Anzahl der einzulesenden Werte
                  (Falls weniger Werte gelesen werden koennen als verlangt,
                  wird ein Fehler gemeldet).
   pnum_max ..... Maximale Anzahl der einzulesenden Werte
                  (Die nicht kopierenden Routinen lesen hoechstens
                  pnum_max Werte aus der uebergebenen Zeile ein,
                  weitere Werte werden ignoriert. Falls weniger als
                  pnum_max Werte vorhanden sind, wird das Lesen
                  abgebrochen, es wird kein Fehler gemeldet,
                  wenn mindestens pnum_min Werte abgesaettigt wurden.
                  Die kopierenden Routinen melden ebenfalls keinen
                  Fehler, wenn mindestens pnum_min Werte abgesaettigt
                  wurden. Falls weniger als pnum_max Werte vorhanden sind,
                  werden die restlichen Werte mit Kopien des zuletzt gelesenen
                  Werts aufgefuellt.

  Resultat:
  nichtkopierende Routinen: Die Anzahl der gelesenen Werte wird zurueckgegeben.
  kopierende Routinen: Die Anzahl der tatsaechlich gelesenen Werte wird
                       zurueckgegeben. Resultat = pnum_max - Anzahl der Kopien
*/

int getparam(char *param_name, void *param, PARAMTYPE ptype, 
             int pnum_min, int pnum_max)
{
  static char errmsg[256];
  char *str;
  int i;
  int numread;

  numread = 0;
  if (ptype == PARAM_STR) {
    str = strtok(NULL," \t\r\n");
    if (str == NULL) {
      sprintf(errmsg,"parameter for %s missing in line %u\nstring expected",
              param_name,curline);
      error(errmsg);
    }
    else strncpy((char *)param,str,pnum_max);
    numread++;
  }
  else if (ptype == PARAM_STRPTR) {
    str = strtok(NULL," \t\r\n");
    if (str == NULL) {
      sprintf(errmsg,"parameter for %s missing in line %u\nstring expected",
              param_name,curline);
      error(errmsg);
    }
    else *((char**)param) = strdup(str);
    numread++;
  }
  else if (ptype == PARAM_INT) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," \t\r\n");
      if (str == NULL) {
        sprintf(errmsg,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg+strlen(errmsg),"integer vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      else ((int*)param)[i] = atoi(str);
      numread++;
    }
    for (i=pnum_min; i<pnum_max; i++) {
      if ((str = strtok(NULL," \t\r\n")) != NULL) {
        ((int*)param)[i] = atoi(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_INT_COPY) {
    int ival = 0;
    for (i=0; i<pnum_max; i++) {
      str = strtok(NULL," \t\r\n");
      if (str != NULL) {
        ival = atoi(str);
        numread++; /* return number of parameters actually read */
      }
      else if (i<pnum_min) {
        sprintf(errmsg,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg+strlen(errmsg),"integer vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      ((int*)param)[i] = ival;
    }
  }
  else if (ptype == PARAM_INTEGER) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," \t\r\n");
      if (str == NULL) {
        sprintf(errmsg,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg+strlen(errmsg),"integer vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      else ((integer*)param)[i] = atoi(str);
      numread++;
    }
    for (i=pnum_min; i<pnum_max; i++) {
      if ((str = strtok(NULL," \t\r\n")) != NULL) {
        ((integer*)param)[i] = atoi(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_INTEGER_COPY) {
    int ival = 0;
    for (i=0; i<pnum_max; i++) {
      str = strtok(NULL," \t\r\n");
      if (str != NULL) {
        ival = atoi(str);
        numread++; /* return number of parameters actually read */
      }
      else if (i<pnum_min) {
        sprintf(errmsg,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg+strlen(errmsg),"integer vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      ((integer*)param)[i] = (integer)ival;
    }
  }
  else if (ptype == PARAM_REAL) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," \t\r\n");
      if (str == NULL) {
        sprintf(errmsg,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg+strlen(errmsg),"real vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      else ((real*)param)[i] = atof(str);
      numread++;
    }
    for (i=pnum_min; i<pnum_max; i++) {
      if ((str = strtok(NULL," \t\r\n")) != NULL) {
        ((real*)param)[i] = atof(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_REAL_COPY) {
    real rval = 0;
    for (i=0; i<pnum_max; i++) {
      str = strtok(NULL," \t\r\n");
      if (str != NULL) {
        rval = atof(str);
        numread++; /* return number of parameters actually read */
      }
      else if (i<pnum_min) {
        sprintf(errmsg,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg+strlen(errmsg),"real vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      ((real*)param)[i] = rval;
    }
  }
  return numread;
} /* getparam */


/*****************************************************************************
*
* read in parameter file in new format (tag based) with name <paramfname>
*
* lines beginning with comment characters '#' or blank lines are skipped
*
*****************************************************************************/

void getparamfile(char *paramfname, int sim)
{
  FILE *pf;
  char buffer[1024];
  char *token;
  char *res;
  str255 tmpstr;
  int tmp;
  real rtmp;
#ifdef MC
  int mc_seed_int;
#endif


#ifdef TWOD
  vektor3d tempforce;
  vektor nullv={0.0,0.0};
  vektor3d tempvek;
  vektor einsv={1.0,1.0};
  vektor3d tempshift;
#else 
  vektor4d tempforce;
  vektor nullv={0.0,0.0,0.0};
  vektor4d tempvek;
  vektor einsv={1.0,1.0,1.0};
  vektor4d tempshift;
#endif
  vektor force;
  vektor vek;
  vektor shift;
  int k;


  int i;

  curline = 0;
  pf = fopen(paramfname,"r");
  if (NULL == pf) {
    sprintf(tmpstr,"Could not open parameter file with name %s\n",paramfname);
    error(tmpstr);
  }

  /* set the random number generator seed to the */
  /* negative of the current time in seconds */
  /* this will be superseeded by a fixed value from the parameter file */
  { 
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv,&tz);
    seed = (long) -tv.tv_sec;
  }

  do {
    res=fgets(buffer,1024,pf);
    if (NULL == res) { finished=1; break; }; /* probably EOF reached */
    curline++;
    token = strtok(buffer," \t\r\n");
    if (NULL == token) continue; /* skip blank lines */
    if (token[0]=='#') continue; /* skip comments */

    if (strcasecmp(token,"simulation")==0) {
      /* file name for atom coordinate input data */
      getparam("simulation",&tmp,PARAM_INT,1,1);
      if (sim < tmp) break;
    }
    else if (strcasecmp(token,"seed")==0) {
      /* seed for random number generator in maxwell */
      int tmp;
      getparam("seed",&tmp,PARAM_INT,1,1);
      seed = (long) tmp;
      if (seed > 0) seed = -seed;
    }
    else if (strcasecmp(token,"coordname")==0) {
      /* file name for atom coordinate input data */
      getparam("coordname",infilename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"outfiles")==0) {
      /* output file basename */
      getparam("outfiles",outfilename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"potfile")==0) {
      /* filename for potential data */
      getparam("potfile",potfilename,PARAM_STR,1,255);
    }
#ifdef DISLOC
    else if (strcasecmp(token,"reffile")==0) {
      /* filename for reference configuration */
      getparam("reffile",reffilename,PARAM_STR,1,255);
    }
#endif
    else if (strcasecmp(token,"ensemble")==0) {
      /* ensemble */
      getparam("ensemble",tmpstr,PARAM_STR,1,255);
      if (strcasecmp(tmpstr,"nve")==0) {
        ensemble = ENS_NVE;
        move_atoms = move_atoms_nve;
      }
      else if (strcasecmp(tmpstr,"mik")==0) {
        ensemble = ENS_MIK;
        move_atoms = move_atoms_mik;
      }
       else if (strcasecmp(tmpstr,"nvt")==0) {
        ensemble = ENS_NVT;
        move_atoms = move_atoms_nvt;
      }
      else if (strcasecmp(tmpstr,"nvx")==0) {
        ensemble = ENS_NVX;
        move_atoms = move_atoms_nvx;
      }
      else if (strcasecmp(tmpstr,"npt_iso")==0) {
        ensemble = ENS_NPT_ISO;
        move_atoms = move_atoms_npt_iso;
      }
      else if (strcasecmp(tmpstr,"npt_axial")==0) {
        ensemble = ENS_NPT_AXIAL;
        move_atoms = move_atoms_npt_axial;
      }
      else if (strcasecmp(tmpstr,"and")==0) {
        error("please use nve ensemble with option and");
      }
      else if (strcasecmp(tmpstr,"mc")==0) {
        ensemble = ENS_MC;
        move_atoms = move_atoms_mc;
      }
      else if (strcasecmp(tmpstr,"frac")==0) {
        ensemble = ENS_FRAC;
        move_atoms = move_atoms_frac;
      }
      else if (strcasecmp(tmpstr,"stm")==0) {
        ensemble = ENS_STM;
        move_atoms = move_atoms_stm;
      } else {
        error("unknown ensemble");
      }
    }
    else if (strcasecmp(token,"maxsteps")==0) {
      /* number of steps for total simulation */
      getparam("maxsteps",&steps_max,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"startstep")==0) {
      /* (re)starting step for the simulation */
      getparam("startstep",&steps_min,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"checkpt_int")==0) {
      /* number of steps between checkpoints / period for checkpoints */
      getparam("checkpt_int",&rep_interval,PARAM_INT,1,1);
    }

    else if (strcasecmp(token,"eng_int")==0) {
      /* energy data output interval */
      getparam("eng_int",&eng_interval,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dist_int")==0) {
      /* number of steps between energy dist. writes */
      getparam("dist_int",&dis_interval,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dist_binary_io")==0) {
      /* binary io flag for energy distributions */
      getparam("dist_binary_io",&dist_binary_io,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"pic_int")==0) {
      /* number of steps between picture writes */
      getparam("pic_int",&pic_interval,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"onl_int")==0) {
      /* number of steps between online visualization */
      getparam("onl_int",&onl_interval,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dist_dim")==0) {
      getparam("dist_dim",&dist_dim,PARAM_INT,DIM,DIM);
    }
    else if (strcasecmp(token,"pbc_dirs")==0) {
      /* directions with periodic boundary conditions */
      getparam("pbc_dirs",&pbc_dirs,PARAM_INT,DIM,DIM);
    }
#ifdef EFILTER
    else if (strcasecmp(token,"ef_checkpt_int")==0) {
      /* number of steps between energy filtered checkpoints */
      getparam("ef_checkpt_int",&efrep_interval,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"e_pot_lower")==0) {
      /* lower end of energy window */
      getparam("e_pot_lower",&lower_e_pot,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"e_pot_upper")==0) {
      /* upper end of energy window */
      getparam("e_pot_upper",&upper_e_pot,PARAM_REAL,1,1);
    }
#endif
    else if (strcasecmp(token,"total_types")==0) {
      /* TOTAL nuber of atoms: ntypes + virtualtypes */
      getparam("total_types",&vtypes,PARAM_INT,1,1);
      restrictions=(vektor*)realloc(restrictions,vtypes*DIM*sizeof(real));
      if (NULL==restrictions)
	error("Cannot allocate memory for restriction vectors\n");
      for(k=0; k<vtypes; k++)
       *(restrictions+k) = einsv;
#ifdef FBC
      /* Allocation & Initialisation of fbc_forces */
      fbc_forces = (vektor *) malloc(vtypes*DIM*sizeof(real));
      if (NULL==fbc_forces)
	error("Cannot allocate memory for fbc_forces\n");
      for(k=0; k<vtypes; k++)
       *(fbc_forces+k) = nullv;

      fbc_beginforces = (vektor *) malloc(vtypes*DIM*sizeof(real));
      if (NULL==fbc_beginforces)
	error("Cannot allocate memory for fbc_beginforces\n");
      for(k=0; k<vtypes; k++)
       *(fbc_beginforces+k) = nullv;
#ifdef MIK
      fbc_dforces = (vektor *) malloc(vtypes*DIM*sizeof(real));
      if (NULL==fbc_dforces)
	error("Cannot allocate memory for fbc_dforces\n");
      for(k=0; k<vtypes; k++)
       *(fbc_dforces+k) = nullv; 
#else
      fbc_endforces = (vektor *) malloc(vtypes*DIM*sizeof(real));
      if (NULL==fbc_endforces)
	error("Cannot allocate memory for fbc_endforces\n");
      for(k=0; k<vtypes; k++)
       *(fbc_endforces+k) = nullv;
#endif
#endif /*FBC*/ 
#ifdef DEFORM
      /* Allocation & Initialisation of deform_shift */
      deform_shift = (vektor *) malloc(vtypes*DIM*sizeof(real));
      if (NULL==deform_shift)
	error("Cannot allocate memory for deform_shift\n");
      for(k=0; k<vtypes; k++)
       *(deform_shift+k) = nullv;
#endif
    }

#ifdef FBC
    else if (strcasecmp(token,"extra_startforce")==0) {
      /* extra force for virtual types */
      /* format: type force.x force.y (force.z) read in a temp. vektor */
      getparam("extra_startforce",&tempforce,PARAM_REAL,DIM+1,DIM+1);
      if (tempforce.x>vtypes-1)
       error("Force defined for non existing virtual atom type\n");
      force.x = tempforce.y;
      force.y = tempforce.z;
#ifndef TWOD
      force.z = tempforce.z2;
#endif
      *(fbc_beginforces+(int)(tempforce.x)) = force;
      *(fbc_forces+(int)(tempforce.x)) = force; 
    }

#ifdef SMIK
    else if (strcasecmp(token,"relaxsteps")==0) {
      /* steps nve integration befor mik  */
      getparam("relaxsteps",&relaxsteps,PARAM_INT,1,1);
    }
#endif

#ifdef MIK
    else if (strcasecmp(token,"fbc_ekin_threshold")==0) {
      /* epsilon criterium to increment extra force*/
      getparam("fbc_ekin_threshold",&fbc_ekin_threshold,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"fbc_annealsteps")==0) {
      /* max nr of steps before shears */
      getparam("fbc_annealsteps",&fbc_annealsteps,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"fbc_waitsteps")==0) {
      /* max nr of steps between shears */
      getparam("fbc_waitsteps",&fbc_waitsteps,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"extra_dforce")==0) {
      /* extra force increment for virtual types */
      /* format: type force.x force.y (force.z) read in a temp. vektor */
      getparam("extra_dforce",&tempforce,PARAM_REAL,DIM+1,DIM+1);
      if (tempforce.x>vtypes-1)
       error("Force increment defined for non existing virtual atom type\n");
      force.x = tempforce.y;
      force.y = tempforce.z;
#ifndef TWOD
      force.z = tempforce.z2;
#endif
      *(fbc_dforces+(int)(tempforce.x)) = force;
    }
#else
    else if (strcasecmp(token,"extra_endforce")==0) {
      /* extra force for virtual types */
      /* format: type force.x force.y (force.z) read in a temp. vektor */
      getparam("extra_endforce",&tempforce,PARAM_REAL,DIM+1,DIM+1);
      if (tempforce.x>vtypes-1)
       error("Force defined for non existing virtual atom type\n");
      force.x = tempforce.y;
      force.y = tempforce.z;
#ifndef TWOD
      force.z = tempforce.z2;
#endif
      *(fbc_endforces+(int)(tempforce.x)) = force;
    }
#endif
#endif /* FBC */

    else if (strcasecmp(token,"restrictionvector")==0) {
      /* restrictions for virtual types */
      /* format: type  1 1 (1) (=all directions ok) read in a temp. vektor */
      getparam("restrictionvector",&tempvek,PARAM_REAL,DIM+1,DIM+1);
      if (tempvek.x>vtypes-1)
       error("Restriction defined for non existing virtual atom type\n");
      vek.x = tempvek.y;
      vek.y = tempvek.z;
#ifndef TWOD
      vek.z = tempvek.z2;
#endif
      *(restrictions+(int)(tempvek.x)) = vek;
    }
    else if (strcasecmp(token,"box_x")==0) {
      /* 'x' or first vector for box */
      getparam("box_x",&box_x,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"box_y")==0) {
      /* 'y' or second vector for box */
      getparam("box_y",&box_y,PARAM_REAL,DIM,DIM);
    }
#ifndef TWOD
    else if (strcasecmp(token,"box_z")==0) {
      /* 'z' or third vector for box */
      getparam("box_z",&box_z,PARAM_REAL,DIM,DIM);
    }
#endif
    else if (strcasecmp(token,"box_param")==0) {
      /*  box parameters for generated structures */
      getparam("box_param",&box_param,PARAM_INT,DIM,DIM);
    }
    else if (strcasecmp(token,"timestep")==0) {
      /* size of timestep (in MD units) */
      getparam("timestep",&timestep,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"ntypes")==0) {
      /* number of atom types */
      getparam("ntypes",&ntypes,PARAM_INT,1,1);
      /*if there are no virtual atoms*/
      if (vtypes==0) vtypes=ntypes;
      restrictions=(vektor*)realloc(restrictions,vtypes*DIM*sizeof(real));
      if (NULL==restrictions)
	error("Cannot allocate memory for restriction vectors\n");
      for(k=0; k<ntypes; k++)
       *(restrictions+k) = einsv;
    }
    else if (strcasecmp(token,"starttemp")==0) {
      /* temperature at start of sim. */
      getparam("starttemp",&temperature,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"use_current_temp")==0) {
      /* set imposed temperature to current system temperature */
      use_curr_temp = 1;
    }
#if defined(AND) || defined(NVT) || defined(NPT) || defined(STM)
    else if (strcasecmp(token,"endtemp")==0) {
      /* temperature at end of sim. */
      getparam("endtemp",&end_temp,PARAM_REAL,1,1);
    }
#endif
#if defined(STM)
    else if (strcasecmp(token,"stadium")==0) {
      getparam("stadium",&stadium,PARAM_REAL,2,2);
    }
    else if (strcasecmp(token,"center")==0) {
      getparam("center",&center,PARAM_REAL,2,2);
    }
#endif
#ifdef MONOLJ
    else if (strcasecmp(token,"r_cut")==0) {
      /* cutoff radius */
      getparam("r_cut",&monolj_r2_cut,PARAM_REAL,1,1);
      monolj_r2_cut = SQR(monolj_r2_cut);
      cellsz = MAX(cellsz,monolj_r2_cut);
    }
    else if (strcasecmp(token,"cellsize")==0) {
      /* minimal cell diameter */
      getparam("cellsize",&rtmp,PARAM_REAL,1,1);
      cellsz = MAX(cellsz,SQR(rtmp));
    }
#endif
    else if (strcasecmp(token,"initsize")==0) {
      /* initial cell size */
      getparam("initsize",&initsz,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"incrsize")==0) {
      /* initial cell size */
      getparam("incrsize",&incrsz,PARAM_INT,1,1);
    }
#ifdef AND
    else if (strcasecmp(token,"tempintv")==0) {
      /* temperature interval */
      getparam("tempintv",&tmp_interval,PARAM_INT,1,1);
    }
#endif
#if defined(NVT) || defined(NPT) || defined(STM)
    else if (strcasecmp(token,"eta")==0) {
      /* eta variable for NVT or NPT thermostat */
      getparam("eta",&eta,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"tau_eta")==0) {
      /* time constant tau_eta for thermostat */
      getparam("tau_eta",&inv_tau_eta,PARAM_REAL,1,1);
      if (inv_tau_eta == (real)0) {
        error("tau_eta is zero.\n");
      }
      inv_tau_eta = 1.0/inv_tau_eta;
    }
    else if (strcasecmp(token,"isq_tau_eta")==0) {
      /* inverse of square of time constant tau_eta for thermostat */
      getparam("isq_tau_eta",&inv_tau_eta,PARAM_REAL,1,1);
      inv_tau_eta = sqrt(inv_tau_eta);
    }
    else if (strcasecmp(token,"inv_tau_eta")==0) {
      /* inverse of time constant tau_eta for thermostat */
      getparam("inv_tau_eta",&inv_tau_eta,PARAM_REAL,1,1);
    }
#ifdef UNIAX
    else if (strcasecmp(token,"eta_rot")==0) {
      /* eta variable of rotational motion for NVT or NPT thermostat */
      getparam("eta_rot",&eta_rot,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"tau_eta_rot")==0) {
      /* time constant tau_eta for thermostat of rotational motion */
      getparam("tau_eta_rot",&inv_tau_eta_rot,PARAM_REAL,1,1);
      if (inv_tau_eta_rot == (real)0) {
        error("tau_eta_rot is zero.\n");
      }
      inv_tau_eta_rot = 1.0/inv_tau_eta_rot;
    }
    else if (strcasecmp(token,"isq_tau_eta_rot")==0) {
      /* squared inverse of time constant for thermostat of rot. motion */
      getparam("isq_tau_eta_rot",&inv_tau_eta_rot,PARAM_REAL,1,1);
      inv_tau_eta_rot = sqrt(inv_tau_eta_rot);
    }
    else if (strcasecmp(token,"inv_tau_eta_rot")==0) {
      /* inverse of time constant for thermostat of rotational motion */
      getparam("inv_tau_eta_rot",&inv_tau_eta_rot,PARAM_REAL,1,1);
    }
#endif
#endif
#ifdef MC
    else if (strcasecmp(token,"mc_beta")==0) {
      /* Monte Carlo inverse temperature */
      getparam("mc_beta",&mc_beta,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"mc_len")==0) {
      /* Monte Carlo mean jump length */
      getparam("mc_len",&mc_len,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"mc_seed")==0) {
      /* seed of Monte Carlo random number generator */
      getparam("mc_seed",&mc_seed_int,PARAM_INT,1,1);
      mc_seed = (long)mc_seed_int;
    }
#endif
#ifdef FRAC
    else if (strcasecmp(token,"dampstadium")==0) { /* Damping stadium */
      getparam("dampstadium",&stadium,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"dampgamma")==0) { /* Damping factor gamma */
      getparam("dampgamma",&gamma_bar,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"dampcutoff")==0) { /* Damping cutoff */
      getparam("dampcutoff",&gamma_cut,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"k_crit")==0) { /* Stress Intensity factor */
      getparam("k_crit",&kcrit,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"mue")==0) { /* Youngs modulus */
      getparam("mue",&mue,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"kel")==0) { /* Shear modulus */
      getparam("kel",&kel,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"crktip")==0) { /* Crack tip location */
      getparam("crktip",&tip,PARAM_REAL,DIM,DIM);
    }
#endif
#ifndef TWOD
    else if (strcasecmp(token,"view_pos")==0) { 
      /* view position */
      getparam("view_pos",&view_pos,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"view_dir")==0) {
      /* view direction */
      getparam("view_dir",&view_dir,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"projection")==0) {
      /* projection (0=orthogonal, 1=perspective) */
      getparam("projection",&projection,PARAM_INT,1,1);
    }
#endif
    else if (strcasecmp(token,"ecut_kin")==0) { 
      /* kinetic energy interval for pictures (min/max) */
      getparam("ecut_kin",&ecut_kin,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"ecut_pot")==0) { 
      /* potential energy interval for pictures (min/max) */
      getparam("ecut_pot",&ecut_pot,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"pic_ll")==0) { 
      /* lower left corner of picture */
      getparam("pic_ll", &pic_ll,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"pic_ur")==0) { 
      /* upper right corner of picture */
      getparam("pic_ur", &pic_ur,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"pic_res")==0) { 
      /* number of pixels in x/y direction */
      getparam("pic_res", &pic_res,PARAM_INT,1,2);
    }
    else if (strcasecmp(token,"numpix")==0) { 
      /* smearing radius in pixels */
      getparam("nsmear", &nsmear,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"pic_type")==0) { 
      /* number of pixels in x/y direction */
      getparam("pic_type", &pic_type,PARAM_INT,1,1);
    }
#ifdef HOMDEF
    else if (strcasecmp(token,"exp_interval")==0) {
      /* period of expansion intervals */
      getparam("exp_interval",&exp_interval,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"expansion")==0) {
      /* expansion */
      getparam("expansion",&expansion,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"hom_interval")==0) {
      /* period of homshear intervals */
      getparam("hom_interval",&hom_interval,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"shear_factor")==0) {
      /* maximum shear */
      getparam("shear_factor",&shear_factor,PARAM_REAL,1,1);
    }
#endif
#if defined(FRAC) || defined(DEFORM)
    else if (strcasecmp(token,"ekin_threshold")==0) {
      /* shear epsilon criterium, see imd_shear_new.c */
      getparam("ekin_threshold",&ekin_threshold,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"annealsteps")==0) {
      /* max nr of steps between shears */
      getparam("annealsteps",&annealsteps,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"max_deform_int")==0) {
      /* max nr of steps between shears */
      getparam("max_deform_int",&max_deform_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"strip_width")==0) { 
      /* strip width (in x dir.) */
      getparam("strip_width",&strip_width,PARAM_REAL,1,1);
    }
#endif
#ifdef SHOCK
    else if (strcasecmp(token,"shock_strip")==0) { 
      /* shock strip width (in x dir.) */
      getparam("shock_strip",&shock_strip,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"shock_speed")==0) { 
      /* shock speed (in x dir.) */
      getparam("shock_speed",&shock_speed,PARAM_REAL,1,1); 
    }
    else if (strcasecmp(token,"shock_mode")==0) { 
       /* shock type: plate or half */
       getparam("shock_mode",&shock_mode,PARAM_INT,1,1); 
       if (shock_mode > 1) shock_strip = 0;
       /* compatibility with old input files */
       if (shock_mode != 2 && shock_mode != 3) shock_mode = 1;
    }
#endif
#ifdef MPI
    else if (strcasecmp(token,"cpu_dim")==0) {
      /* CPU array dimension */
      getparam("cpu_dim",&cpu_dim,PARAM_INT,DIM,DIM);    
    }
    else if (strcasecmp(token,"parallel_output")==0) {
      /* parallel output flag */
      getparam("parallel_output",&parallel_output,PARAM_INT,1,1);    }
    else if (strcasecmp(token,"parallel_input")==0) {
      /* parallel input flag */
      getparam("parallel_input",&parallel_input,PARAM_INT,1,1);    }
#endif
#ifdef CORRELATE
    else if (strcasecmp(token,"correl_rmax")==0) {
      /* dimension of histogram in r domain */
      getparam("correl_rmax",&ncorr_rmax,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"correl_tmax")==0) {
      /* dimension of histogram in t domain */
      getparam("correl_tmax",&ncorr_tmax,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"correl_int")==0) {
      /* repeat interval for correlation */
      getparam("correl_int",&correl_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"correl_omode")==0) {
      /* repeat interval for correlation */
      getparam("correl_omode",&correl_omode,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"GS_rcut")==0) {
      /* cutoff radius for correlation data writes */
      getparam("GS_rcut",&GS_rcut,PARAM_REAL,1,1);
    }
#endif
#if defined(CORRELATE) || defined(MSQD)
    else if (strcasecmp(token,"correl_start")==0) {
      /* start time for correlation */
      getparam("correl_start",&correl_start,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"correl_end")==0) {
      /* end time for correlation */
      getparam("correl_end",&correl_end,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"correl_ts")==0) {
      /* sampling time interval for correlation */
      getparam("correl_ts",&correl_ts,PARAM_INT,1,1);
    }
#endif
#ifdef NVX
    else if (strcasecmp(token, "dTemp_start")==0){
      /* temperature asymmetry at start */
      getparam("dTemp_start", &dTemp_start, PARAM_REAL, 1,1);
    }
    else if (strcasecmp(token, "dTemp_end")==0){
      /* temperature asymmetry at end */
      getparam("dTemp_end", &dTemp_end, PARAM_REAL, 1,1);
    }
#endif
#ifdef RNEMD
    else if (strcasecmp(token, "exch_interval")==0){
      /* interval for particle exchange */
      getparam("exch_interval", &exch_interval, PARAM_INT, 1,1);
    }
#endif
#ifdef TRANSPORT
    else if (strcasecmp(token, "tran_nlayers")==0){
      /* number of layers */
      getparam("tran_nlayers", &tran_nlayers, PARAM_INT, 1,1);
    }
     else if (strcasecmp(token, "tran_interval")==0){
      /* number of steps between temp. writes  */
      getparam("tran_interval", &tran_interval, PARAM_INT, 1,1);
    }
#endif
#ifdef STRESS_TENS
    else if (strcasecmp(token, "press_dim")==0){
      /* pressure histogram dimension */
      getparam("press_dim", &press_dim, PARAM_INT, DIM,DIM);
    }
     else if (strcasecmp(token, "press_interval")==0){
      /*number of steps between press. writes  */
      getparam("press_interval", &press_interval, PARAM_INT, 1,1);
    }
#endif
#ifdef DISLOC
    else if (strcasecmp(token,"dem_int")==0) {
      /* number of steps between picture writes */
      getparam("dem_int",&dem_interval,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dsp_int")==0) {
      /* number of steps between picture writes */
      getparam("dsp_int",&dsp_interval,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"up_ort_ref")==0) {
      /* step number to compute ort_ref */
      getparam("update_ort_ref",&up_ort_ref,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"min_dpot")==0) {
      /* minimum Epot difference */
      getparam("min_dpot",&min_dpot,PARAM_REAL,1,1);
    }   
    else if (strcasecmp(token,"dpotsorte")==0) {
      /* atom type for Epot difference */
      getparam("dpotsorte",&dpotsorte,PARAM_INT,1,1);
    }   
    else if (strcasecmp(token,"ddelta")==0) {
      /* minimum distance for ddm */
      getparam("ddelta",&ddelta,PARAM_REAL,1,1);
    }   
    else if (strcasecmp(token,"reset_Epot_step")==0) {
      /* step at which to compute Epot_ref (if calc_Epot_ref==1) */
      getparam("reset_Epot_step",&reset_Epot_step,PARAM_INT,1,1);
    }   
    else if (strcasecmp(token,"calc_Epot_ref")==0) {
      /* read (0) or compute (1) reference potential energy */
      getparam("calc_Epot_ref",&calc_Epot_ref,PARAM_INT,1,1);
    }   
    else if (strcasecmp(token,"Epot_diff")==0) {
      /* write Epot (0) or Epot_diff (1) */
      getparam("Epot_diff",&Epot_diff,PARAM_INT,1,1);
    }   
#endif
#ifdef AVPOS
    else if (strcasecmp(token,"avpos_res")==0) {
      /* number of steps between coordinate addition */
      getparam("avpos_res",&avpos_res,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"avpos_int")==0) {
      /* number of steps between average position writes */
      getparam("avpos_int",&avpos_int,PARAM_INT,1,1);
    }
#endif
#ifdef ORDPAR
    else if (strcasecmp(token,"op_rcut")==0) {
      /* cutoff radius for order parameter */
      getparam("op_rcut",&op_r2_cut,PARAM_REAL,4,4);
      op_r2_cut[0][0] = SQR(op_r2_cut[0][0]);
      op_r2_cut[0][1] = SQR(op_r2_cut[0][1]);
      op_r2_cut[1][0] = SQR(op_r2_cut[1][0]);
      op_r2_cut[1][1] = SQR(op_r2_cut[1][1]);
    }   
    else if (strcasecmp(token,"op_weight")==0) {
      /* weights for order parameter */
      getparam("op_weight",&op_weight,PARAM_REAL,4,4);
    }   
#endif
#ifdef DEFORM
    else if (strcasecmp(token,"deform_shift")==0) {
      /* deform shift for virtual types */
      /* format: type shift.x shift.y (shift.z) read in a temp. vektor */
      getparam("deform_shift",&tempshift,PARAM_REAL,DIM+1,DIM+1);
      if (tempshift.x>vtypes-1)
       error("Shift defined for non existing virtual atom type\n");
      shift.x = tempshift.y;
      shift.y = tempshift.z;
#ifndef TWOD
      shift.z = tempshift.z2;
#endif
      *(deform_shift+(int)(tempshift.x)) = shift; 
    }
#endif
#ifdef USE_SOCKETS
    else if (strcasecmp(token,"socket_int")==0) {
      getparam("socket_int",&socket_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"display_host")==0) {
      getparam("display_host",display_host,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"baseport")==0) { int tmp;
      getparam("baseport",&tmp,PARAM_INT,1,1);
      baseport = tmp; /* conversion to unsigned short */
    }
    else if (strcasecmp(token,"use_socket_window")==0) {
      getparam("use_socket_window",&use_socket_window,PARAM_INT,1,1);
    }

#endif
#ifdef NPT
    else if (strcasecmp(token,"xi")==0) {
      /* xi variable for NPT thermostat */
      getparam("xi",&xi,PARAM_REAL,1,DIM);
    }
    else if (strcasecmp(token,"pressure_start")==0) {
      /* external starting pressure or stress for NPT */
      getparam("pressure_start",&pressure_ext,PARAM_REAL_COPY,1,DIM);
    }
    else if (strcasecmp(token,"use_current_pressure")==0) {
      /* set imposed pressure to current system pressure */
      use_curr_pressure = 1;
    }
    else if (strcasecmp(token,"pressure_end")==0) {
      /* external end pressure or stress for NPT */
      getparam("pressure_end",&pressure_end,PARAM_REAL_COPY,1,DIM);
    }
    else if (strcasecmp(token,"tau_xi")==0) {
      /* time constant tau_xi for NPT thermostat algorithm */
      getparam("tau_xi",&inv_tau_xi,PARAM_REAL,1,1);
      if (inv_tau_xi == (real)0) {
        error("tau_xi is zero.\n");
      }
      inv_tau_xi = 1.0/inv_tau_xi;
    }
    else if (strcasecmp(token,"isq_tau_xi")==0) {
      /* inverse of square of time constant tau_xi for NPT thermostat */
      getparam("isq_tau_xi",&inv_tau_xi,PARAM_REAL,1,1);
      inv_tau_xi = sqrt(inv_tau_xi);
    }
    else if (strcasecmp(token,"inv_tau_xi")==0) {
      /* inverse of time constant tau_xi for NPT thermostat */
      getparam("inv_tau_xi",&inv_tau_xi,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"cell_size_tol")==0) {
      /* rel. tolerance for volume rescaling during NPT sim. */
      getparam("cell_size_tol",&cell_size_tolerance,PARAM_REAL,1,1);
    }
#endif
#ifdef EAM
    else if (strcasecmp(token,"eam_len")==0) {
      /* EAM: number of neighbours */
      getparam("eam_len",&eam_len,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"eam_A")==0) {
      /* EAM: constant for cohesive function */
      getparam("eam_A",&eam_A,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"eam_r_cut")==0) {
      /* EAM: cutoff of cohesive function */
      getparam("eam_r_cut",&eam_r_cut,PARAM_REAL,1,1);
      eam_r2_cut = SQR(eam_r_cut);
      cellsz = MAX(cellsz,eam_r2_cut);
    }
    else if (strcasecmp(token,"eam_r_0")==0) {
      /* EAM: minimum distance of cohesive function */
      getparam("eam_r_0",&eam_r_0,PARAM_REAL,1,1);
    }
#endif
#ifdef EAM2
    else if (strcasecmp(token,"core_potential_file")==0) {
      /* EAM2:Filename for the tabulated Core-Core Potenial (r^2) */
      getparam("core_potential_file",eam2_core_pot_filename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"embedding_energy_file")==0) {
      /* EAM2:Filename for the tabulated Embedding Enery(rho_h) */
      getparam("embedding_energy_file",eam2_emb_E_filename,PARAM_STR,1,255);
    }
   else if (strcasecmp(token,"atomic_e-density_file")==0) {
      /* EAM2:Filename for the tabulated atomic electron density(r_ij^2) */
      getparam("atomic_e-density_file",eam2_at_rho_filename,PARAM_STR,1,255);
    }
#endif
#ifdef COVALENT
    else if (strcasecmp(token,"neigh_len")==0) {
      /* number of neighbors */
      getparam("neigh_len",&neigh_len,PARAM_INT,1,1);
    }
#endif
#ifdef TTBP
    else if (strcasecmp(token,"ttbp_constant")==0) {
      /* force constant (radians); type 0 */
      getparam("ttbp_constant",ttbp_constant,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"ttbp_sp")==0) {
      /* hybridization of the element type */
      getparam("ttbp_sp",ttbp_sp,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"ttbp_potfile")==0) {
      /* filename for ttbp potential data */
      getparam("ttbp_potfile",ttbp_potfilename,PARAM_STR,1,255);
    }
#endif
#ifdef TERSOFF
    /* Parameters for Tersoff potential */
    else if (strcasecmp(token,"ters_r_cut")==0) {     
      getparam("ters_r_cut",ters_r_cut,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"ters_r0")==0) {     
      getparam("ters_r0",ters_r0,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"ters_a")==0) {     
      getparam("ters_a",ters_a,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"ters_b")==0) {     
      getparam("ters_b",ters_b,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"ters_la")==0) {     
      getparam("ters_la",ters_la,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"ters_mu")==0) {     
      getparam("ters_mu",ters_mu,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"ters_chi")==0) {     
      getparam("ters_chi",ters_chi,PARAM_REAL,ntypes*(ntypes-1)/2,ntypes*(ntypes-1)/2);
    }
    else if (strcasecmp(token,"ters_om")==0) {     
      getparam("ters_om",ters_om,PARAM_REAL,ntypes*(ntypes-1)/2,ntypes*(ntypes-1)/2);
    }
    else if (strcasecmp(token,"ters_ga")==0) {     
      getparam("ters_ga",ters_ga,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"ters_n")==0) {     
      getparam("ters_n",ters_n,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"ters_c")==0) {     
      getparam("ters_c",ters_c,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"ters_d")==0) {     
      getparam("ters_d",ters_d,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"ters_h")==0) {     
      getparam("ters_h",ters_h,PARAM_REAL,ntypes,ntypes);
    }
#endif
#ifdef UNIAX
    else if (strcasecmp(token,"uniax_r_cut")==0) {
      /* UNIAX: cutoff radius of uniaxial molecules */
      getparam("uniax_r_cut",&uniax_r_cut,PARAM_REAL,1,1);
      uniax_r2_cut = SQR(uniax_r_cut);
      cellsz = MAX(cellsz,uniax_r2_cut);
    }
#endif 
    else {
      char msg[255];
      sprintf(msg,"****** Unknown TAG %s ignored ******",token);
      warning(msg);
    }
  } while (!feof(pf));
  if (feof(pf)) finished=1;
  fclose(pf);

} /* getparamfile */


/*****************************************************************
*
*   Check input for nonsense values
*
******************************************************************/

void check_parameters_complete()
{
  real tmp;

  if (ensemble == 0) {
    error("missing or unknown ensemble parameter.");
  }
  if (timestep == (real)0) {
    error("timestep is missing or zero.");
  }
  if (ntypes == 0) {
    error("ntypes is missing or zero.");
  }
#if defined(NPT) || defined(NVT) || defined(STM)
  if (temperature == 0) {
    error("starttemp is missing or zero.");
  }
#endif
#if defined(NPT) || defined(NVT) || defined(STM)
  if (end_temp == 0) {
    error("endtemp is missing or zero.");
  }
#endif
#if defined(CORRELATE) || defined(MSQD)
  if (correl_ts == 0) {
    if (eng_interval != 0) correl_ts = eng_interval;
    else {
      error("correl_ts is missing or zero.");
    }
  }
#endif
#ifdef CORRELATE
  if (ncorr_rmax == 0) {
    error("correl_rmax is missing or zero.");
  }
  if (ncorr_tmax == 0) {
    error("correl_tmax is zero.");
  }
#endif
#ifdef NVX
  	if (dTemp_start == 0){
		error ("dTemp_start is missing or zero.");
	}
	if (dTemp_end == 0){
		error ("dTemp_end is missing or zero.");
	}
	if (tran_interval == 0){
		error ("tran_interval is zero.");
	}
	if (tran_nlayers == 0){
		error ("tran_nlayers is zero.");
        }
#endif
#ifdef RNEMD
	if (tran_interval == 0){
		error ("tran_interval is zero.");
	}
	if (tran_nlayers == 0){
		error ("tran_nlayers is zero.");
        }
#endif
#ifdef STRESS_TENS
	if (press_interval == 0) {
		error ("press_interval is zero.");
	}
#endif
#ifdef MPI
#ifdef TWOD
        if ((cpu_dim.x==0) || (cpu_dim.y==0))
#else
        if ((cpu_dim.x==0) || (cpu_dim.y==0) || (cpu_dim.z==0))
#endif
	{
           error("cpu_dim is missing or zero.");
        }
#endif
#ifdef USE_SOCKETS
  if (display_host[0]=='\0') {
    error("display_host name or IP address missing.");
  }
#endif
#ifdef UNIAX
  if (uniax_r_cut == 0) {
    error("uniax_r_cut is missing or zero.");
  }
#endif
#ifdef MONOLJ
  if (monolj_r2_cut == 0) {
    error("monolj_r_cut is missing or zero.");
  }
  /* determine shift of potential at cutoff */
  pair_int_monolj(&monolj_shift,&tmp,monolj_r2_cut);
#endif
}

/*****************************************************************
*
*  read command line and first set of parameters
*
******************************************************************/

void read_parameters(int argc,char **argv)

{
  str255 fname;
  FILE *infile;
#if defined(__GNUC__) && defined(__STRICT_ANSI__)
  extern char *strdup(const char *);
#endif

  if ( 0 == myid ) { /* Read Parameters on Master Process */

    /* Check for Restart, process options */
    strcpy(progname,argv[0]);
    while ((argc > 1) && (argv[1][0] =='-')) {
      switch (argv[1][1]) {
        /* r - restart */
        case 'r':
          if (argv[1][2]=='\0') {
            if (NULL != argv[2]) {
              restart = atoi(argv[2]);
              --argc;
              ++argv;
            }
          }
          else restart = atoi(&argv[1][2]);
          break;
        case 'p':
          if (argv[1][2]=='\0') {
            if (NULL != argv[2]) {
              paramfilename = strdup(argv[2]);
              --argc;
              ++argv;
            }
          }
          else paramfilename = strdup(&argv[1][2]);
          break;
        default:
          printf("Illegal option %s \n",argv[1]);
          usage();
          exit(-1);
      }
      ++argv;
      --argc;
    }

    getparamfile(paramfilename,1);
    check_parameters_complete();

    /* Get restart parameters if restart */
    if (0 != restart) {
      sprintf(fname,"%s.%d.itr",outfilename,restart);
      sprintf(infilename,"%s.%d.%s",outfilename,restart,"chkpt");
      printf("Restarting from %s.\n",infilename);
      getparamfile(fname,1);
    } else {
      /* Delete energy file if not restart */
      sprintf(fname,"%s.eng",outfilename);
      unlink(fname);
      /* Delete distrib minmax file if not restart */
      sprintf(fname,"%s.minmax.dist",outfilename);
      unlink(fname);
      /* Delete tempdist file if not restart */
      sprintf(fname,"%s.tempdist",outfilename);
      unlink(fname);
      /* write header to minmax file, if we possibly need it */
      if ((dist_dim.x>0) && (dis_interval>0)) { 
        FILE *fl;
        fl = fopen(fname,"w");
#ifdef TWOD
        fprintf(fl,"# Dimension %d %d\n",dist_dim.x,dist_dim.y);
#else
        fprintf(fl,"# Dimension %d %d %d\n",dist_dim.x,dist_dim.y,dist_dim.z);
#endif
        fclose(fl);
      }
#ifdef MSQD
      /* Delete msqd file if not restart */
      sprintf(fname,"%s.msqd",outfilename);
      unlink(fname);
#endif
    }
  }
}


#ifdef MPI

/****************************************************************************
*
*  Broadcast all parameters to other CPU's (MPI only) 
*
*****************************************************************************/

void broadcast_params() {

  MPI_Bcast( &finished    , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &ensemble    , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &simulation  , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &seed        , 1, MPI_LONG, 0, MPI_COMM_WORLD); 

  MPI_Bcast( &steps_max   , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &steps_min   , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &restart     , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &rep_interval, 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &eng_interval, 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &dis_interval, 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &pic_interval, 1, MPI_INT,  0, MPI_COMM_WORLD); 

#ifdef TWOD
  MPI_Bcast( &pic_scale   , 2, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ecut_kin    , 2, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ecut_pot    , 2, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &pic_ll      , 2, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &pic_ur      , 2, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &pic_res     , 2, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &pic_type    , 1, MPI_INT, 0, MPI_COMM_WORLD); 
#endif

#ifdef EFILTER
  MPI_Bcast( &lower_e_pot,    1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &upper_e_pot,    1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &efrep_interval, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

  MPI_Bcast( &vtypes,         1, MPI_INT, 0, MPI_COMM_WORLD);

#ifdef FBC
  if (0!=myid) fbc_forces  = (vektor *) malloc(vtypes*DIM*sizeof(real));
  if (NULL==fbc_forces) 
    error("Cannot allocate memory for fbc_forces on client."); 
  MPI_Bcast( fbc_forces, vtypes*DIM, REAL, 0, MPI_COMM_WORLD);
 
  if (0!=myid) fbc_beginforces  = (vektor *) malloc(vtypes*DIM*sizeof(real));
  if (NULL==fbc_beginforces) 
    error("Cannot allocate memory for fbc_beginforces on client."); 
  MPI_Bcast( fbc_beginforces, vtypes*DIM, REAL, 0, MPI_COMM_WORLD); 
#ifdef SMIK
  MPI_Bcast( &relaxsteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
#ifdef MIK
  MPI_Bcast( &fbc_ekin_threshold , 1, REAL,    0, MPI_COMM_WORLD); 
  MPI_Bcast( &fbc_waitsteps      , 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &fbc_annealsteps    , 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (0!=myid) fbc_dforces  = (vektor *) malloc(vtypes*DIM*sizeof(real));
  if (NULL==fbc_dforces) 
    error("Cannot allocate memory for fbc_dforces on client."); 
  MPI_Bcast( fbc_dforces, vtypes*DIM, REAL, 0, MPI_COMM_WORLD); 
#else
  if (0!=myid) fbc_endforces  = (vektor *) malloc(vtypes*DIM*sizeof(real));
  if (NULL==fbc_endforces) 
    error("Cannot allocate memory for fbc_endforces on client."); 
  MPI_Bcast( fbc_endforces, vtypes*DIM, REAL, 0, MPI_COMM_WORLD); 
#endif
#endif
  if (0!=myid) restrictions  = (vektor *) malloc(vtypes*DIM*sizeof(real));
  if (NULL==restrictions) 
    error("Cannot allocate memory for restriction vectors on client."); 
  MPI_Bcast( restrictions, vtypes*DIM, REAL, 0, MPI_COMM_WORLD);  


  MPI_Bcast( &pbc_dirs    , DIM, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &box_x       , DIM, REAL,     0, MPI_COMM_WORLD); 
  MPI_Bcast( &box_y       , DIM, REAL,     0, MPI_COMM_WORLD);
#ifndef TWOD
  MPI_Bcast( &box_z       , DIM, REAL,     0, MPI_COMM_WORLD);
#endif 
  MPI_Bcast( &box_param   , DIM, MPI_INT,  0, MPI_COMM_WORLD); 

  MPI_Bcast( &timestep    ,   1, REAL,     0, MPI_COMM_WORLD); 
  MPI_Bcast( &ntypes      ,   1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &temperature ,   1, REAL,     0, MPI_COMM_WORLD); 
  MPI_Bcast( &cpu_dim     , DIM, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &dist_dim    , DIM, MPI_INT,  0, MPI_COMM_WORLD);

  MPI_Bcast( &parallel_output, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &parallel_input,  1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( outfilename , sizeof(outfilename), MPI_CHAR, 0, MPI_COMM_WORLD); 
  MPI_Bcast( infilename  , sizeof(infilename) , MPI_CHAR, 0, MPI_COMM_WORLD); 
  MPI_Bcast( reffilename , sizeof(reffilename), MPI_CHAR, 0, MPI_COMM_WORLD); 
  MPI_Bcast( potfilename , sizeof(potfilename), MPI_CHAR, 0, MPI_COMM_WORLD); 
#ifdef TTBP
  MPI_Bcast( ttbp_potfilename , sizeof(ttbp_potfilename), 
             MPI_CHAR, 0, MPI_COMM_WORLD); 
#endif

#if defined(AND) || defined(NVT) || defined(NPT) || defined(STM)
  MPI_Bcast( &end_temp,      1, REAL, 0, MPI_COMM_WORLD); 
#endif
#ifdef MONOLJ
  MPI_Bcast( &monolj_r2_cut, 1, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &monolj_shift,  1, REAL, 0, MPI_COMM_WORLD); 
#endif
  MPI_Bcast( &cellsz, 1, REAL,     0, MPI_COMM_WORLD); 
  MPI_Bcast( &initsz, 1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &incrsz, 1, MPI_INT,  0, MPI_COMM_WORLD);

#ifdef AND
  MPI_Bcast( &tmp_interval, 1, MPI_INT, 0, MPI_COMM_WORLD); 
#endif

#if defined(NVT) || defined(NPT) || defined(STM)
  MPI_Bcast( &eta ,         1 , REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &inv_tau_eta , 1 , REAL, 0, MPI_COMM_WORLD); 
#ifdef UNIAX
  MPI_Bcast( &eta_rot ,         1 , REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &inv_tau_eta_rot , 1 , REAL, 0, MPI_COMM_WORLD); 
#endif
#endif

#if defined(STM)
  MPI_Bcast( &stadium ,         2 , REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &center ,          2 , REAL, 0, MPI_COMM_WORLD); 
#endif

#ifdef NPT
  MPI_Bcast( &xi,                DIM, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &inv_tau_xi,          1, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &pressure_ext,      DIM, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &pressure_end,      DIM, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &cell_size_tolerance, 1, REAL, 0, MPI_COMM_WORLD); 
#endif

#ifdef MC
  MPI_Bcast( &mc_beta     , 1, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &mc_len      , 1, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &mc_seed     , 1, MPI_LONG, 0, MPI_COMM_WORLD); 
#endif

#if defined(CORRELATE)
  MPI_Bcast( &correl_omode, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &correl_int,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &correl_tmax,  1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &correl_rmax,  1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

#if defined(CORRELATE) || defined(MSQD)
  MPI_Bcast( &correl_start, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &correl_end,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &correl_ts,    1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

#ifdef NVX
  MPI_Bcast( &dTemp_start,   1, REAL,   0, MPI_COMM_WORLD); 
  MPI_Bcast( &dTemp_end,     1, REAL,   0, MPI_COMM_WORLD); 
#endif
#ifdef RNEMD
  MPI_Bcast( &exch_interval, 1, MPI_INT,  0, MPI_COMM_WORLD);
#endif
#ifdef TRANSPORT
  MPI_Bcast( &tran_nlayers,  1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &tran_interval, 1, MPI_INT,  0, MPI_COMM_WORLD);
#endif

#ifdef STRESS_TENS
  MPI_Bcast( &press_dim,    DIM, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &press_interval, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

#ifdef FRAC
  MPI_Bcast( &stadium , DIM, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &gamma_bar , 1, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &gamma_cut , 1, REAL, 0, MPI_COMM_WORLD);
#endif

#if defined(FRAC) || defined(DEFORM)
  MPI_Bcast( &strip_width, 1, REAL, 0, MPI_COMM_WORLD); 
  if (0!=myid) deform_shift  = (vektor *) malloc(vtypes*DIM*sizeof(real));
  if (NULL==deform_shift) 
    error("Cannot allocate memory for deform_shift on client."); 
  MPI_Bcast( deform_shift, vtypes*DIM, REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef SHOCK
  MPI_Bcast( &shock_strip, 1, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &shock_speed, 1, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &shock_mode,  1, MPI_INT, 0, MPI_COMM_WORLD); 
#endif

#ifdef FRAC
  MPI_Bcast( &kcrit , 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &mue   , 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &kel   , 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &tip.x , 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &tip.y , 1, REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef DISLOC
  MPI_Bcast( &min_dpot,        1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ddelta,          1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dem_interval,    1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dsp_interval,    1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &calc_Epot_ref,   1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &reset_Epot_step, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &Epot_diff,       1, MPI_INT, 0, MPI_COMM_WORLD); 
#endif

#ifdef AVPOS
  MPI_Bcast( &avp_int,         1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &avp_res,         1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

#ifdef ORDPAR
  MPI_Bcast( &op_r2_cut,       4, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &op_weight,       4, REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef HOMDEF
  MPI_Bcast( &hom_interval , 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &shear_factor , 1, REAL,    0, MPI_COMM_WORLD); 
  MPI_Bcast( &exp_interval , 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &expansion ,  DIM, REAL,    0, MPI_COMM_WORLD); 
#endif
#if defined(FRAC) || defined(DEFORM)
  MPI_Bcast( &ekin_threshold , 1, REAL,    0, MPI_COMM_WORLD); 
  MPI_Bcast( &annealsteps ,    1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &max_deform_int , 1, MPI_INT, 0, MPI_COMM_WORLD); 
#endif  

#ifdef USE_SOCKETS
  MPI_Bcast( &socket_int, 1, MPI_INT, 0, MPI_COMM_WORLD); 
#endif

#ifdef UNIAX
  MPI_Bcast( &uniax_r_cut,  1, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &uniax_r2_cut, 1, REAL, 0, MPI_COMM_WORLD); 
#endif

#ifdef EAM
  MPI_Bcast( &eam_len,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &eam_A,     1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &eam_r_cut, 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &eam_r2_cut,1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &eam_r_0,   1, REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef COVALENT
  MPI_Bcast( &neigh_len, 1, MPI_INT,  0, MPI_COMM_WORLD);
#endif

#ifdef TTBP
  MPI_Bcast( ttbp_constant,  ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ttbp_sp,        ntypes, REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef TERSOFF
  MPI_Bcast( ters_r_cut, ntypes*(ntypes+1)/2, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_r0,    ntypes*(ntypes+1)/2, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_a,     ntypes*(ntypes+1)/2, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_b,     ntypes*(ntypes+1)/2, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_la,    ntypes*(ntypes+1)/2, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_mu,    ntypes*(ntypes+1)/2, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_chi,   ntypes*(ntypes-1)/2, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_om,    ntypes*(ntypes-1)/2, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_ga,                 ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_n,                  ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_c,                  ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_d,                  ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_h,                  ntypes, REAL, 0, MPI_COMM_WORLD);
#endif
  
  /* broadcast integrator to other CPU's */

  switch (ensemble) {
    case ENS_NVE:       move_atoms = move_atoms_nve;       break;
    case ENS_MIK:       move_atoms = move_atoms_mik;       break;
    case ENS_NVT:       move_atoms = move_atoms_nvt;       break;
    case ENS_NPT_ISO:   move_atoms = move_atoms_npt_iso;   break;
    case ENS_NPT_AXIAL: move_atoms = move_atoms_npt_axial; break;
    case ENS_MC:        move_atoms = move_atoms_mc;        break;
    case ENS_FRAC:      move_atoms = move_atoms_frac;      break;
    case ENS_NVX:       move_atoms = move_atoms_nvx;       break;
    case ENS_STM:       move_atoms = move_atoms_stm;       break;  
    default: if (0==myid) error("unknown ensemble in broadcast"); break;
  }
  
}

#endif /* MPI */







