/*****************************************************************************
*
* read in parameter files (tag based)                               MH 260298
* 
* $RCSfile$
* $Revision$
* $Date$
*
******************************************************************************/

#include <stdio.h>
#include <string.h>

#include "imd.h"

#if defined(__GNUC__) && defined(__STRICT_ANSI__)
extern char *strdup(char *);
#endif

#ifdef __WATCOMC__
#define strcasecmp strcmpi
#endif

/* To do yet: improve checking of bad input files ! */
/* e.g. for forgotten starttemp/endtemp              */
/* clean up prototypes.h                            */

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
    str = strtok(NULL," \t\n");
    if (str == NULL) {
      fprintf(stderr,"parameter for %s missing in line %u\n",
              param_name,curline);
      error("string expected\n");
    }
    else strncpy((char *)param,str,pnum_max);
    numread++;
  }
  else if (ptype == PARAM_STRPTR) {
    str = strtok(NULL," \t\n");
    if (str == NULL) {
      fprintf(stderr,"parameter for %s missing in line %u\n",
              param_name,curline);
      error("string expected\n");
    }
    else *((char**)param) = strdup(str);
    numread++;
  }
  else if (ptype == PARAM_INT) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," \t\n");
      if (str == NULL) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"integer vector of dim %u expected\n",
                (unsigned)pnum_min);
        error(errmsg);
      }
      else ((int*)param)[i] = atoi(str);
      numread++;
    }
    for (i=pnum_min; i<pnum_max; i++) {
      if ((str = strtok(NULL," \t\n")) != NULL) {
        ((int*)param)[i] = atoi(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_INT_COPY) {
    int ival = 0;
    for (i=0; i<pnum_max; i++) {
      str = strtok(NULL," \t\n");
      if (str != NULL) {
        ival = atoi(str);
        numread++; /* return number of parameters actually read */
      }
      else if (i<pnum_min) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"integer vector of dim %u expected\n",
                (unsigned)pnum_min);
        error(errmsg);
      };
      ((int*)param)[i] = ival;
    }
  }
  else if (ptype == PARAM_INTEGER) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," \t\n");
      if (str == NULL) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"integer vector of dim %u expected\n",
                (unsigned)pnum_min);
        error(errmsg);
      }
      else ((integer*)param)[i] = atoi(str);
      numread++;
    }
    for (i=pnum_min; i<pnum_max; i++) {
      if ((str = strtok(NULL," \t\n")) != NULL) {
        ((integer*)param)[i] = atoi(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_INTEGER_COPY) {
    int ival = 0;
    for (i=0; i<pnum_max; i++) {
      str = strtok(NULL," \t\n");
      if (str != NULL) {
        ival = atoi(str);
        numread++; /* return number of parameters actually read */
      }
      else if (i<pnum_min) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"integer vector of dim %u expected\n",
                (unsigned)pnum_min);
        error(errmsg);
      };
      ((integer*)param)[i] = (integer)ival;
    }
  }
  else if (ptype == PARAM_REAL) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," \t\n");
      if (str == NULL) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"real vector of dim %u expected\n",
                (unsigned)pnum_min);
        error(errmsg);
      }
      else ((real*)param)[i] = atof(str);
      numread++;
    }
    for (i=pnum_min; i<pnum_max; i++) {
      if ((str = strtok(NULL," \t\n")) != NULL) {
        ((real*)param)[i] = atof(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_REAL_COPY) {
    real rval = 0;
    for (i=0; i<pnum_max; i++) {
      str = strtok(NULL," \t\n");
      if (str != NULL) {
        rval = atof(str);
        numread++; /* return number of parameters actually read */
      }
      else if (i<pnum_min) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"real vector of dim %u expected\n",
                (unsigned)pnum_min);
        error(errmsg);
      };
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
#ifdef MC
  int mc_seed_int;
#endif
  int i;

  curline = 0;
  pf = fopen(paramfname,"r");
  if (NULL == pf) {
    perror("getparam");
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    exit(10);
  };

  do {
    res=fgets(buffer,1024,pf);
    if (NULL == res) { finished=1; break; }; /* probably EOF reached */
    curline++;
    token = strtok(buffer," \t\n");
    if (NULL == token) continue; /* skip blank lines */
    if (token[0]=='#') continue; /* skip comments */

    if (strcasecmp(token,"simulation")==0) {
      /* file name for atom coordinate input data */
      getparam("simulation",&tmp,PARAM_INT,1,1);
      if (sim < tmp) break;
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
      else if (strcasecmp(tmpstr,"msd")==0) {
        ensemble = ENS_MSD;
        move_atoms = move_atoms_msd;
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
        ensemble = ENS_AND;
        move_atoms = move_atoms_and;
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
    }
    else if (strcasecmp(token,"starttemp")==0) {
      /* temperature at start of sim. */
      getparam("starttemp",&temperature,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"use_current_temp")==0) {
      /* set imposed temperature to current system temperature */
      use_curr_temp = 1;
    }
    else if (strcasecmp(token,"pn")==0) {
      /* z/y (3/2D)-coordinate of glideplane */
      getparam("pn",&pn,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"upperplane")==0) {
      /* z/y (3/2D)-coordinate of glideplane */
      getparam("upperplane",&upperplane,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"lowerplane")==0) {
      /* z/y (3/2D)-coordinate of glideplane */
      getparam("lowerplane",&lowerplane,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"burgersv")==0) {
      /* length of Burgers-vector */
      getparam("burgersv",&burgersv,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"width")==0) {
      /* width of dislocation */
      getparam("width",&width,PARAM_REAL,1,1);
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
#if MONOLJ
    else if (strcasecmp(token,"r_cut")==0) {
      /* cutoff radius */
      getparam("r_cut",&r2_cut,PARAM_REAL,1,1);
      r2_cut = SQR(r2_cut);
    }
    else if (strcasecmp(token,"cellsize")==0) {
      /* cell dimension */
      getparam("cellsize",&cellsz,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"initsize")==0) {
      /* initial cell size */
      getparam("initsize",&initsz,PARAM_INT,1,1);
    }
#endif
#if AND
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
      /* time constant for thermostat */
      getparam("tau_eta",&isq_tau_eta,PARAM_REAL,1,1);
      if (isq_tau_eta == (real)0) {
        error("tau_eta is zero.\n");
      }
      isq_tau_eta = (real)1/SQR(isq_tau_eta);
    }
    else if (strcasecmp(token,"isq_tau_eta")==0) {
      /* time constant for thermostat */
      getparam("isq_tau_eta",&isq_tau_eta,PARAM_REAL,1,1);
    }
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
    else if (strcasecmp(token,"pic_ll")==0) { 
      /* lower left front corner of configuration */
      getparam("pic_ll", &pic_ll,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"pic_ur")==0) { 
      /* upper right back corner of configuration */
      getparam("pic_ur", &pic_ur,PARAM_REAL,DIM,DIM);
    }
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
      /* number of pixels in x/y direction */
      getparam("numpix", &numpix,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"pic_type")==0) { 
      /* number of pixels in x/y direction */
      getparam("pic_type", &pic_type,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"pic_at_radius")==0) {
      /* radius (in A) for atoms in pictures */
      if (ntypes <= 0) {
        error("ntypes is less or equal zero or missing\n");
      }
      else {
        pic_at_radius = calloc(ntypes,sizeof(real));
        if (pic_at_radius == NULL) {
          error("Cannot allocate memory\n");
        };
      };
      getparam("pic_at_radius", pic_at_radius,PARAM_REAL_COPY,1,ntypes);
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
    else if (strcasecmp(token,"initial_shift")==0) {
      /* shall the whole sample be shifted before MD */
      getparam("initial_shift",&initial_shift,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"ins")==0) {
      /* shift vector (whole sample) */
      getparam("ins",&ins,PARAM_REAL,DIM,DIM);
    }   
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
    else if (strcasecmp(token,"strip_width")==0) { 
      /* strip width (in x dir.) */
      getparam("strip_width",&strip_width,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"shock_strip")==0) { 
      /* shock strip width (in x dir.) */
      getparam("shock_strip",&shock_strip,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"shock_speed")==0) { 
      /* shock speed (in x dir.) */
      getparam("shock_speed",&shock_speed,PARAM_REAL,1,1); 
    }
    else if (strcasecmp(token,"shock_elong")==0) { 
       /* shock elong (in x dir.) */
       getparam("shock_elong",&shock_elong,PARAM_REAL,1,1); 
    }
    else if (strcasecmp(token,"shock_mode")==0) { 
       /* shock elong (in x dir.) */
       getparam("shock_mode",&shock_mode,PARAM_INT,1,1); 
    }
#endif
#ifdef MPI
    else if (strcasecmp(token,"cpu_dim")==0) {
      /* CPU array dimension */
      getparam("cpu_dim",&cpu_dim,PARAM_INT,DIM,DIM);    }
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
      getparam("correl_rmax",&ncorr_rmax,PARAM_INTEGER,1,1);
    }
    else if (strcasecmp(token,"correl_tmax")==0) {
      /* dimension of histogram in t domain */
      getparam("correl_tmax",&ncorr_tmax,PARAM_INTEGER,1,1);
    }
    else if (strcasecmp(token,"correl_int")==0) {
      /* repeat interval for correlation */
      getparam("correl_int",&correl_int,PARAM_INTEGER,1,1);
    }
    else if (strcasecmp(token,"correl_omode")==0) {
      /* repeat interval for correlation */
      getparam("correl_omode",&correl_omode,PARAM_INTEGER,1,1);
    }
    else if (strcasecmp(token,"GS_rcut")==0) {
      /* cutoff radius for correlation data writes */
      getparam("GS_rcut",&GS_rcut,PARAM_REAL,1,1);
    }
#endif
#if defined(CORRELATE) || defined(MSQD)
    else if (strcasecmp(token,"correl_start")==0) {
      /* start time for correlation */
      getparam("correl_start",&correl_start,PARAM_INTEGER,1,1);
    }
    else if (strcasecmp(token,"correl_end")==0) {
      /* end time for correlation */
      getparam("correl_end",&correl_end,PARAM_INTEGER,1,1);
    }
    else if (strcasecmp(token,"correl_ts")==0) {
      /* sampling time interval for correlation */
      getparam("correl_ts",&correl_ts,PARAM_INTEGER,1,1);
    }
#endif
#ifdef TRANSPORT
    else if (strcasecmp(token, "dTemp_start")==0){
      /*  temperature ...*/
      getparam("dTemp_start", &dTemp_start, PARAM_REAL, 1,1);
    }
    else if (strcasecmp(token, "dTemp_end")==0){
      /* temperature ...*/
      getparam("dTemp_end", &dTemp_end, PARAM_REAL, 1,1);
    }
    else if (strcasecmp(token, "tran_nlayers")==0){
      /*number of layer  */
      getparam("tran_nlayers", &tran_nlayers, PARAM_INTEGER, 1,1);
    }
     else if (strcasecmp(token, "tran_interval")==0){
      /*number of steps between temp. writes  */
      getparam("tran_interval", &tran_interval, PARAM_INTEGER, 1,1);
    }
#endif
#ifdef STRESS_TENS
    else if (strcasecmp(token, "press_dim")==0){
      /* pressure histogram dimension */
      getparam("press_dim", &press_dim, PARAM_INTEGER, DIM,DIM);
    }
     else if (strcasecmp(token, "press_interval")==0){
      /*number of steps between press. writes  */
      getparam("press_interval", &press_interval, PARAM_INTEGER, 1,1);
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
      getparam("dpotsorte",&dpotsorte,PARAM_INTEGER,1,1);
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
    else if (strcasecmp(token,"strip_shift")==0) {
      /* strip move per timestep - this is a vector */
      getparam("strip_shift",&strip_shift,PARAM_REAL,DIM,DIM); 
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
      /* time constant tau for NPT thermostat algorithm */
      getparam("tau_xi",&isq_tau_xi,PARAM_REAL,1,1);
      if (isq_tau_xi == (real)0) {
        error("tau_xi is zero.\n");
      };
      isq_tau_xi = (real)1/SQR(isq_tau_xi);
    }
    else if (strcasecmp(token,"isq_tau_xi")==0) {
      /* inverse of square of time constant tau for NPT thermostat */
      getparam("isq_tau_xi",&isq_tau_xi,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"cell_size_tol")==0) {
      /* rel. tolerance for volume rescaling during NPT sim. */
      getparam("cell_size_tol",&cell_size_tolerance,PARAM_REAL,1,1);
    }
#endif
#ifdef EAM
    else if (strcasecmp(token,"eam_len")==0) {
      /* EAM: number of neighbours */
      getparam("eam_len",&eam_len,PARAM_INTEGER,1,1);
    }
    else if (strcasecmp(token,"eam_A")==0) {
      /* EAM: constant for cohesive function */
      getparam("eam_A",&eam_A,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"eam_r_cut")==0) {
      /* EAM: cutoff of cohesive function */
      getparam("eam_r_cut",&eam_r_cut,PARAM_REAL,1,1);
      eam_r2_cut = eam_r_cut * eam_r_cut;
    }
    else if (strcasecmp(token,"eam_r_0")==0) {
      /* EAM: minimum distance of cohesive function */
      getparam("eam_r_0",&eam_r_0,PARAM_REAL,1,1);
    }
#endif
#ifdef TTBP
    else if (strcasecmp(token,"ttbp_len")==0) {
      /* number of neighbors */
      getparam("ttbp_len",&ttbp_len,PARAM_INTEGER,1,1);
    }
    else if (strcasecmp(token,"ttbp_constant")==0) {
      /* force constant (radians); type 0 */
      getparam("ttbp_constant",&ttbp_constant,PARAM_REAL,ntypes,ntypes);
      for (i=0;i<ntypes;i++) {
        ttbp_constant[i] = ttbp_constant[i];
      };
    }
    else if (strcasecmp(token,"ttbp_theta")==0) {
      /* equilibrium angle in radians */
      getparam("ttbp_theta",&ttbp_theta,PARAM_REAL,ntypes,ntypes);
      for (i=0;i<ntypes;i++) {
        ttbp_theta[i] = ttbp_theta[i];
      };
    }
    else if (strcasecmp(token,"ttbp_potfile")==0) {
      /* filename for ttbp potential data */
      getparam("ttbp_potfile",ttbp_potfilename,PARAM_STR,1,255);
    }
#endif
    else {
      fprintf(stderr,"**** WARNING: Unknown TAG %s ignored ****\n",token);
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
  if (ensemble == 0) {
    error("missing or unknown ensemble parameter.\n");
  };
  if (timestep == (real)0) {
    error("timestep is missing or zero.\n");
  };
  if (ntypes == 0) {
    error("ntypes is missing or zero.\n");
  };
#if defined(NPT) || defined(NVT) || defined(STM)
  if (temperature == 0) {
    error("starttemp is missing or zero.\n");
  };
#endif
#if defined(NPT) || defined(NVT) || defined(STM)
  if (end_temp == 0) {
    error("endtemp is missing or zero.\n");
  };
#endif
#if defined(CORRELATE) || defined(MSQD)
  if (correl_ts == 0) {
    if (eng_interval != 0) correl_ts = eng_interval;
    else {
      error("correl_ts is missing or zero.\n");
    };
  };
#endif
#ifdef CORRELATE
  if (ncorr_rmax == 0) {
    error("correl_rmax is missing or zero.\n");
  }
  if (ncorr_tmax == 0) {
    error("correl_tmax is zero.\n");
  }
#endif
#ifdef TRANSPORT
  	if (dTemp_start == 0){
		error ("dTemp_start is missing or zero. \n ");
	}
	if (dTemp_end == 0){
		error ("dTemp_end is missing or zero. \n ");
	}
	if (tran_interval == 0){
		error ("tran_interval is zero. \n");
	}
	if (tran_nlayers == 0){
		error ("tran_nlayers is zero. \n");
        }
#endif
#ifdef STRESS_TENS
	if (press_interval == 0) {
		error ("press_interval is zero. \n");
	}
#endif
#ifdef MPI
#ifdef TWOD
        if ((cpu_dim.x==0) || (cpu_dim.y==0))
#else
        if ((cpu_dim.x==0) || (cpu_dim.y==0) || (cpu_dim.z==0))
#endif
	{
           error("cpu_dim is missing or zero\n");
        }
#endif
#ifdef USE_SOCKETS
  if (display_host[0]=='\0') {
    error("display_host name or IP address missing.\n");
  }
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

#ifdef MPI
  if ( 0 == myid ) { /* Read Parameters on Master Process */
#endif

/* Check for Restart, process options */

  strcpy(progname,argv[0]);
  while ((argc > 1) && (argv[1][0] =='-')) {
    switch (argv[1][1]) {
      /* r - restart */
    case 'r':
      restart = atoi(&argv[1][2]);
      break;
    case 'p':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          paramfilename = strdup(argv[2]);
          --argc;
          ++argv;
        };
      }
      else paramfilename = strdup(&argv[1][2]);
#ifdef DEBUG
      if (paramfilename) printf("%s\n",paramfilename);
#endif
      break;
    default:
      printf("Illegal option %s \n",argv[1]);
      usage();
      exit(-1);
    }
    ++argv;
    --argc;
  };

#ifdef DEBUG
  printf("getparamfile(%s);\n",paramfilename);
#endif
  getparamfile(paramfilename,1);
  check_parameters_complete();

  /* Get restart parameters if restart */
  if (0 != restart) {
    sprintf(fname,"%s.%d.itr",outfilename,restart);
    sprintf(infilename,"%s.%d",outfilename,restart);
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
  };

#ifdef MPI
  };
#endif

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

  MPI_Bcast( &steps_max   , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &steps_min   , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &restart     , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &rep_interval, 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &eng_interval, 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &dis_interval, 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &pic_interval, 1, MPI_INT,  0, MPI_COMM_WORLD); 

#ifdef TWOD
  MPI_Bcast( &pic_scale   , 2, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ecut_kin    , 2, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ecut_pot    , 2, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &pic_ll      , 2, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &pic_ur      , 2, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &pic_res     , 2, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &pic_type    , 1, MPI_INT,  0, MPI_COMM_WORLD); 
#endif

  MPI_Bcast( &box_x       , DIM, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &box_y       , DIM, MPI_REAL, 0, MPI_COMM_WORLD);
#ifndef TWOD
  MPI_Bcast( &box_z       , DIM, MPI_REAL, 0, MPI_COMM_WORLD);
#endif 
  MPI_Bcast( &box_param   , DIM, MPI_INT,  0, MPI_COMM_WORLD); 

  MPI_Bcast( &timestep    , 1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &ntypes      , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &temperature , 1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &cpu_dim     , DIM, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &dist_dim    , DIM, MPI_INT,  0, MPI_COMM_WORLD);

  MPI_Bcast( &parallel_output, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &parallel_input,  1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( outfilename , sizeof(outfilename), MPI_CHAR, 0, MPI_COMM_WORLD); 
  MPI_Bcast( infilename  , sizeof(infilename) , MPI_CHAR, 0, MPI_COMM_WORLD); 
  MPI_Bcast( reffilename , sizeof(reffilename), MPI_CHAR, 0, MPI_COMM_WORLD); 

#if defined(AND) || defined(NVT) || defined(NPT) || defined(STM)
  MPI_Bcast( &end_temp , 1 , MPI_REAL,  0, MPI_COMM_WORLD); 
#endif
#ifdef MONOLJ
  MPI_Bcast( &r2_cut     , 1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &cellsz      , 1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &initsz, 1, INTEGER, 0, MPI_COMM_WORLD);
#endif

#ifdef AND
  MPI_Bcast( &tmp_interval , 1 , MPI_INT, 0, MPI_COMM_WORLD); 
#endif

#if defined(NVT) || defined(NPT) || defined(STM)
  MPI_Bcast( &eta ,         1 , MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &isq_tau_eta , 1 , MPI_REAL, 0, MPI_COMM_WORLD); 
#endif

#if defined(STM)
  MPI_Bcast( &stadium ,         2 , MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &center ,          2 , MPI_REAL, 0, MPI_COMM_WORLD); 
#endif

#ifdef NPT
  MPI_Bcast( &xi,                DIM, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &isq_tau_xi,          1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &pressure_ext,      DIM, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &cell_size_tolerance, 1, MPI_REAL, 0, MPI_COMM_WORLD); 
#endif

#ifdef MC
  MPI_Bcast( &mc_beta     , 1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &mc_len      , 1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &mc_seed     , 1, MPI_LONG, 0, MPI_COMM_WORLD); 
#endif

#if defined(CORRELATE)
  MPI_Bcast( &correl_omode, 1, INTEGER, 0, MPI_COMM_WORLD);
  MPI_Bcast( &correl_int,   1, INTEGER, 0, MPI_COMM_WORLD);
  MPI_Bcast( &correl_tmax,  1, INTEGER, 0, MPI_COMM_WORLD);
  MPI_Bcast( &correl_rmax,  1, INTEGER, 0, MPI_COMM_WORLD);
#endif

#if defined(CORRELATE) || defined(MSQD)
  MPI_Bcast( &correl_start, 1, INTEGER, 0, MPI_COMM_WORLD);
  MPI_Bcast( &correl_end,   1, INTEGER, 0, MPI_COMM_WORLD);
  MPI_Bcast( &correl_ts,    1, INTEGER, 0, MPI_COMM_WORLD);
#endif

#ifdef TRANSPORT
  MPI_Bcast( &dTemp_start,   1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &dTemp_end,     1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &tran_nlayers,  1, INTEGER,  0, MPI_COMM_WORLD);
  MPI_Bcast( &tran_interval, 1, INTEGER,  0, MPI_COMM_WORLD);
#endif

#ifdef STRESS_TENS
  MPI_Bcast( &press_dim,      1, INTEGER,  0, MPI_COMM_WORLD);
  MPI_Bcast( &press_interval, 1, INTEGER,  0, MPI_COMM_WORLD);
#endif

#ifdef FRAC
  MPI_Bcast( &stadium   , DIM, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &gamma_bar , 1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &gamma_cut , 1, MPI_REAL, 0, MPI_COMM_WORLD);
#endif

#if defined(FRAC) || defined(DEFORM) || defined(SHOCK)
  MPI_Bcast( &strip_width, 1, MPI_REAL, 0, MPI_COMM_WORLD); 
#endif

#ifdef SHOCK
  MPI_Bcast( &shock_strip, 1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &shock_speed, 1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &shock_elong, 1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &shock_mode, 1, INTEGER, 0, MPI_COMM_WORLD); 
#endif

#ifdef FRAC
  MPI_Bcast( &kcrit , 1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &mue   , 1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &kel   , 1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &tip.x , 1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &tip.y , 1, MPI_REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef DISLOC
  MPI_Bcast( &min_dpot,        1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ddelta,          1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dem_interval,    1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &dsp_interval,    1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &calc_Epot_ref,   1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &reset_Epot_step, 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &Epot_diff,       1, MPI_INT,  0, MPI_COMM_WORLD); 
#endif

#ifdef ORDPAR
  MPI_Bcast( &op_r2_cut,       4, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &op_weight,       4, MPI_REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef HOMDEF
  MPI_Bcast( &exp_interval , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &expansion ,  DIM, MPI_REAL, 0, MPI_COMM_WORLD); 
#endif
#ifdef DEFORM
  MPI_Bcast( &hom_interval , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &shear_max ,    1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &strip_shift ,DIM, MPI_REAL, 0, MPI_COMM_WORLD); 
#endif
#if defined(FRAC) || defined(DEFORM)
  MPI_Bcast( &initial_shift ,  1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &ins ,          DIM, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &ekin_threshold , 1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &annealsteps ,    1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &max_deform_int , 1, MPI_INT,  0, MPI_COMM_WORLD); 
#endif  

#ifdef USE_SOCKETS
  MPI_Bcast( &socket_int, 1, MPI_INT, 0, MPI_COMM_WORLD); 
#endif

#ifdef EAM
  MPI_Bcast( &eam_len,   1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &eam_A,     1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &eam_r_cut, 1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &eam_r_0,   1, MPI_REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef TTBP
  MPI_Bcast( &ttbp_len,      1,      INTEGER,   0, MPI_COMM_WORLD);
  MPI_Bcast( &ttbp_constant, ntypes, MPI_REAL,  0, MPI_COMM_WORLD);
  MPI_Bcast( &ttbp_theta,    ntypes, MPI_REAL,  0, MPI_COMM_WORLD);
#endif

  /* broadcast integrator to other CPU's */

  switch (ensemble) {
    case ENS_NVE:       move_atoms = move_atoms_nve;       break;
    case ENS_MIK:       move_atoms = move_atoms_mik;       break;
    case ENS_NVT:       move_atoms = move_atoms_nvt;       break;
    case ENS_NPT_ISO:   move_atoms = move_atoms_npt_iso;   break;
    case ENS_NPT_AXIAL: move_atoms = move_atoms_npt_axial; break;
    case ENS_AND:       move_atoms = move_atoms_and;       break;
    case ENS_MC:        move_atoms = move_atoms_mc;        break;
    case ENS_FRAC:      move_atoms = move_atoms_frac;      break;
    case ENS_NVX:       move_atoms = move_atoms_nvx;       break;
    case ENS_MSD:       move_atoms = move_atoms_msd;       break;
    case ENS_STM:       move_atoms = move_atoms_stm;       break;  
    default: if (0==myid) error("unknown ensemble in broadcast"); break;
  }
  
}

#endif /* MPI */
