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
      else if (strcasecmp(tmpstr,"pull")==0) {
        ensemble = ENS_PULL;
        move_atoms = move_atoms_pull;
      }
      else if (strcasecmp(tmpstr,"mikshear")==0) {
        ensemble = ENS_MIKSHEAR;
        move_atoms = move_atoms_mik;
      }
      else {
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
#if defined(AND) || defined(NVT) || defined(NPT)
    else if (strcasecmp(token,"endtemp")==0) {
      /* temperature at end of sim. */
      getparam("endtemp",&end_temp,PARAM_REAL,1,1);
    }
#endif
#if AND
    else if (strcasecmp(token,"tempintv")==0) {
      /* temperature interval */
      getparam("tempintv",&tmp_interval,PARAM_INT,1,1);
    }
#endif
#if defined(NVT) || defined(NPT)
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
    else if (strcasecmp(token,"dampstadion")==0) { /* Damping stadion */
      getparam("dampstadion",&stadion,PARAM_REAL,DIM,DIM);
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
#ifdef TWOD
    else if (strcasecmp(token,"ecut_kin")==0) { 
      /* kinetic energy interval for pictures (min/max) */
      getparam("ecut_kin",&ecut_kin,PARAM_REAL,2,2);
    }
    else if (strcasecmp(token,"ecut_pot")==0) { 
      /* potential energy interval for pictures (min/max) */
      getparam("ecut_pot",&ecut_pot,PARAM_REAL,2,2);
    }
    else if (strcasecmp(token,"pic_scale")==0) { 
      /* picture scale (x y) */
      getparam("pic_scale", &pic_scale,PARAM_REAL,2,2);
    }
    else if (strcasecmp(token,"pic_ll")==0) { 
      /* lower left corner of picture */
      getparam("pic_ll", &pic_ll,PARAM_REAL,2,2);
    }
    else if (strcasecmp(token,"pic_ur")==0) { 
      /* upper right corner of picture */
      getparam("pic_ur", &pic_ur,PARAM_REAL,2,2);
    }
    else if (strcasecmp(token,"pic_res")==0) { 
      /* number of pixels in x/y direction */
      getparam("pic_res", &pic_res,PARAM_INT,1,2);
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
#endif
#if defined(FRAC) || defined(PULL) || defined(SHOCK) || defined(MIKSHEAR)
    else if (strcasecmp(token,"strip_width")==0) { 
      /* strip width (in x dir.) */
      getparam("strip_width",&strip,PARAM_REAL,1,1);
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
    else if (strcasecmp(token,"shock_elong")==0) { 
       /* shock elong (in x dir.) */
       getparam("shock_elong",&shock_elong,PARAM_REAL,1,1); 
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
#ifdef MIKSHEAR
    else if (strcasecmp(token,"shear_delta")==0) {
      /* shear delta per timestep */
      getparam("shear_delta",&shear_delta,PARAM_REAL,1,1);
    }   
    else if (strcasecmp(token,"shear_epsilon")==0) {
      /* shear epsilon criterium, see imd_shear_new.c */
      getparam("shear_epsilon",&shear_epsilon,PARAM_REAL,1,1);
    }   
    else if (strcasecmp(token,"glideplane")==0) {
      /* y component of glide plane coord. */
      getparam("glideplane",&glideplane,PARAM_REAL,1,1);
    }   
#endif
#ifdef DISLOC
    else if (strcasecmp(token,"dem_int")==0) {
      /* number of steps between picture writes */
      getparam("dem_int",&dem_interval,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"ddm_int")==0) {
      /* number of steps between picture writes */
      getparam("ddm_int",&ddm_interval,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dsp_int")==0) {
      /* number of steps between picture writes */
      getparam("dsp_int",&dsp_interval,PARAM_INT,1,1);
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
#ifdef PULL
    else if (strcasecmp(token,"strip_shift")==0) {
      /* strip move per timestep - this is a vector */
      getparam("strip_shift",&delta,PARAM_REAL,DIM,DIM); 
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
        ttbp_constant[i] = ttbp_constant[i] * 3.1415927 / 180.;
      };
    }
    else if (strcasecmp(token,"ttbp_theta")==0) {
      /* equilibrium angle in radians */
      getparam("ttbp_theta",&ttbp_theta,PARAM_REAL,ntypes,ntypes);
      for (i=0;i<ntypes;i++) {
        ttbp_theta[i] = ttbp_theta[i] * 3.1415927 / 180.;
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
#if defined(NPT) || defined(NVT)
  if (temperature == 0) {
    error("starttemp is missing or zero.\n");
  };
#endif
#if defined(NPT) || defined(NVT)
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
#ifdef MPI
  if ((cpu_dim.x==0) || (cpu_dim.y==0)
#ifndef TWOD
      || (cpu_dim.z==0)
#endif
     ) {
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
    printf("Restarting from %s.\n",fname);
    getparamfile(fname,1);
  } else {
    /* Delete energy file if not restart */
    sprintf(fname,"%s.eng",outfilename);
    unlink(fname);
    /* Delete distrib minmax file if not restart */
    sprintf(fname,"%s.minmax.dist",outfilename);
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

  MPI_Bcast( &timestep    , 1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &ntypes      , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &temperature , 1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &cpu_dim     , DIM, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &dist_dim    , DIM, MPI_INT,  0, MPI_COMM_WORLD);

  MPI_Bcast( &parallel_output, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &parallel_input,  1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( outfilename , sizeof(outfilename), MPI_CHAR, 0, MPI_COMM_WORLD); 
  MPI_Bcast( infilename  , sizeof(infilename) , MPI_CHAR, 0, MPI_COMM_WORLD); 

#if defined(AND) || defined(NVT) || defined(NPT)
  MPI_Bcast( &end_temp , 1 , MPI_REAL,  0, MPI_COMM_WORLD); 
#endif

#ifdef AND
  MPI_Bcast( &tmp_interval , 1 , MPI_INT, 0, MPI_COMM_WORLD); 
#endif

#if defined(NVT) || defined(NPT)
  MPI_Bcast( &isq_tau_eta , 1 , MPI_REAL, 0, MPI_COMM_WORLD); 
#endif

#ifdef NPT
  MPI_Bcast( &isq_tau_xi         , 1, MPI_REAL, 0, MPI_COMM_WORLD); 
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

#ifdef FRAC
  MPI_Bcast( &stadion   , DIM, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &gamma_bar , 1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &gamma_cut , 1, MPI_REAL, 0, MPI_COMM_WORLD);
#endif

#if defined(FRAC) || defined(PULL) || defined(SHOCK) || defined(MIKSHEAR)
  MPI_Bcast( &strip, 1, MPI_REAL, 0, MPI_COMM_WORLD); 
#endif

#ifdef SHOCK
  MPI_Bcast( &shock_strip, 1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &shock_speed, 1, MPI_REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &shock_elong, 1, MPI_REAL, 0, MPI_COMM_WORLD); 
#endif

#ifdef FRAC
  MPI_Bcast( &kcrit , 1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &mue   , 1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &kel   , 1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &tip.x , 1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &tip.y , 1, MPI_REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef MIKSHEAR
  MPI_Bcast( &shear_delta,   1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &shear_epsilon, 1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &glideplane,    1, MPI_REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef DISLOC
  MPI_Bcast( &min_dpot,        1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ddelta,          1, MPI_REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dem_interval,    1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &ddm_interval,    1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &dsp_interval,    1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &calc_Epot_ref,   1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &reset_Epot_step, 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &Epot_diff,       1, MPI_INT,  0, MPI_COMM_WORLD); 
#endif

#ifdef PULL
  MPI_Bcast( &delta , DIM, MPI_REAL, 0, MPI_COMM_WORLD); 
#endif  

#ifdef USE_SOCKETS
  MPI_Bcast( &socket_int, 1, MPI_INT, 0, MPI_COMM_WORLD); 
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
    case ENS_PULL:      move_atoms = move_atoms_pull;      break;
    case ENS_MIKSHEAR:  move_atoms = move_atoms_mik;       break;
    /* the following should never be executed */
    default: if (0==myid) error("unknown ensemble in broadcast"); break;
  };
  
}

#endif /* MPI */
