
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

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
  vektor shear, base;
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
      /* get number of the simulation phase */
      getparam("simulation",&tmp,PARAM_INT,1,1);
      if (sim < tmp) break;
    }
    else if (strcasecmp(token,"loop")==0) {
      /* looping for online visualisation */
      getparam(token,&loop,PARAM_INT,1,1);
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
    else if (strcasecmp(token,"itrame")==0) {
      /* file name for initial itr-file */
      getparam("itrname",itrfilename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"outfiles")==0) {
      /* output file basename */
      getparam("outfiles",outfilename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"potfile")==0) {
      /* filename for potential data */
      getparam("potfile",potfilename,PARAM_STR,1,255);
      have_potfile = 1;
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
      else if (strcasecmp(tmpstr,"frac")==0) {
        ensemble = ENS_FRAC;
        move_atoms = move_atoms_frac;
      }
      else if (strcasecmp(tmpstr,"ftg")==0) {
        ensemble = ENS_FTG;
        move_atoms = move_atoms_ftg;
      }
      else if (strcasecmp(tmpstr,"finnis")==0) {
        ensemble = ENS_FINNIS;
        move_atoms = move_atoms_finnis;
      }
      else if (strcasecmp(tmpstr,"sllod")==0) {
        ensemble = ENS_SLLOD;
        move_atoms = move_atoms_sllod;
      }
      else if (strcasecmp(tmpstr,"stm")==0) {
        ensemble = ENS_STM;
        move_atoms = move_atoms_stm;
      } 
      else if (strcasecmp(tmpstr,"cg")==0) {
        ensemble = ENS_CG;
        /* the following have different types; */
        /* move_atoms_cg is never used as move_atoms */
        /* move_atoms = move_atoms_cg; */
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
      getparam("checkpt_int",&checkpt_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"eng_int")==0) {
      /* energy data output interval */
      getparam("eng_int",&eng_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dist_int")==0) {
      /* number of steps between energy dist. writes */
      getparam(token,&dist_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dist_dim")==0) {
      /* dimension of distributions */
      getparam(token,&dist_dim,PARAM_INT,DIM,DIM);
    }
    else if (strcasecmp(token,"dist_ll")==0) {
      /* lower left corner of distributions */
      getparam(token,&dist_ll,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"dist_ur")==0) {
      /* upper right corner of distribution */
      getparam(token,&dist_ur,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"dist_Ekin_flag")==0) {
      /* write Ekin dist? */
      getparam(token,&dist_Ekin_flag,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dist_Epot_flag")==0) {
      /* write Epot dist? */
      getparam(token,&dist_Epot_flag,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dist_press_flag")==0) {
      /* write press dist? */
      getparam(token,&dist_press_flag,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dist_presstens_flag")==0) {
      /* write presstens dist? */
      getparam(token,&dist_presstens_flag,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dist_Ekin_long_flag")==0) {
      /* write longitudinal Ekin dist? */
      getparam(token,&dist_Ekin_long_flag,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dist_Ekin_trans_flag")==0) {
      /* write transversal Ekin dist? */
      getparam(token,&dist_Ekin_trans_flag,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dist_shock_shear_flag")==0) {
      /* write shock shear dist? */
      getparam(token,&dist_shock_shear_flag,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"pic_int")==0) {
      /* number of steps between picture writes */
      getparam("pic_int",&pic_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"pbc_dirs")==0) {
      /* directions with periodic boundary conditions */
      getparam("pbc_dirs",&pbc_dirs,PARAM_INT,DIM,DIM);
    }
#ifdef NBLIST
    else if (strcasecmp(token,"nbl_margin")==0) {
      /* margin of neighbor list */
      getparam(token,&nbl_margin,PARAM_REAL,1,1);
    }
#endif
#ifdef VEC
    else if (strcasecmp(token,"atoms_per_cpu")==0) {
      /* maximal number of atoms per CPU */
      getparam(token,&atoms_per_cpu,PARAM_INT,1,1);
    }
#endif
#ifdef EFILTER
    else if (strcasecmp(token,"ef_checkpt_int")==0) {
      /* number of steps between energy filtered checkpoints */
      getparam("ef_checkpt_int",&ef_checkpt_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"e_pot_lower")==0) {
      if (ntypes==0) error("specify parameter ntypes before e_pot_lower");
      getparam("e_pot_lower",lower_e_pot,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"e_pot_upper")==0) {
      if (ntypes==0) error("specify parameter ntypes before e_pot_upper");
      getparam("e_pot_upper",upper_e_pot,PARAM_REAL,ntypes,ntypes);
    }
#endif
#ifdef CLONE
    else if (strcasecmp(token,"nclones")==0) {
      /* number of clones to deal with*/
      getparam(token,&nclones,PARAM_INT,1,1);
    }
#endif
#ifdef NBFILTER
    else if (strcasecmp(token,"nb_checkpt_int")==0) {
      getparam("nb_checkpt_int",&nb_checkpt_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"nb_cut_lower")==0) {
      if (ntypes==0) error("specify parameter ntypes before nb_cut_lower");
      getparam("nb_cut_lower",lower_nb_cut,PARAM_INT,ntypes,ntypes);
    }
    else if (strcasecmp(token,"nb_cut_upper")==0) {
      if (ntypes==0) error("specify parameter ntypes before nb_cut_upper");
      getparam("nb_cut_upper",upper_nb_cut,PARAM_INT,ntypes,ntypes);
    }
#endif
#ifdef SNAPSHOT
    else if (strcasecmp(token,"sscount")==0) {
      /* actual snapshot nr., for restarting */
      getparam("sscount",&sscount,PARAM_INT,1,1);
    }
#endif
    else if (strcasecmp(token,"total_types")==0) {
      /* TOTAL nuber of atom types: ntypes + virtual types */
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
      deform_shift = (vektor *) malloc( vtypes * DIM * sizeof(real) );
      if (NULL==deform_shift)
	error("Cannot allocate memory for deform_shift\n");
      for(k=0; k<vtypes; k++)
       *(deform_shift+k) = nullv;

      /* Allocation & Initialisation of shear_def */
      shear_def = (int *) malloc( vtypes * sizeof(int) );
      if (NULL==shear_def)
	error("Cannot allocate memory for shear_def\n");
      for(k=0; k<vtypes; k++)
       *(shear_def+k) = 0;

      /* Allocation & Initialisation of deform_shear */
      deform_shear = (vektor *) malloc( vtypes * DIM * sizeof(real) );
      if (NULL==deform_shear)
	error("Cannot allocate memory for deform_shear\n");
      for(k=0; k<vtypes; k++)
       *(deform_shear+k) = nullv;

      /* Allocation & Initialisation of deform_base */
      deform_base = (vektor *) malloc( vtypes * DIM * sizeof(real) );
      if (NULL==deform_base)
	error("Cannot allocate memory for deform_base\n");
      for(k=0; k<vtypes; k++)
       *(deform_base+k) = nullv;
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
      /* box parameters for generated structures */
      getparam("box_param",&box_param,PARAM_INT,DIM,DIM);
    }
    else if (strcasecmp(token,"box_unit")==0) {
      /* lattice parameter for generated structures */
      getparam("box_unit",&box_unit,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"masses")==0) {
      /* masses for generated structures */
      if (ntypes==0) 
        error("specify parameter ntypes before parameter masses");
      getparam("masses",masses,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"timestep")==0) {
      /* size of timestep (in MD units) */
      getparam("timestep",&timestep,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"ntypes")==0) {
      /* number of atom types */
      getparam("ntypes",&ntypes,PARAM_INT,1,1);
#ifdef MONO
      if (ntypes!=1) error("this executable is for monoatomic systems only!");
#endif
#ifdef BINARY
      if (ntypes!=2) error("this executable is for binary systems only!");
#endif
      ntypepairs = ((ntypes+1)*ntypes)/2;
      ntypetriples = ntypes * ntypepairs;      
      /* if there are no virtual atom types */
      if (vtypes==0) vtypes=ntypes;
      restrictions=(vektor*)realloc(restrictions,vtypes*DIM*sizeof(real));
      if (NULL==restrictions)
	error("Cannot allocate memory for restriction vectors\n");
      for(k=0; k<ntypes; k++)
        *(restrictions+k) = einsv;
      /* array of masses for generated structures */
      masses=(real*)realloc(masses,ntypes*sizeof(real));
      if (NULL==masses)
	error("Cannot allocate memory for masses array\n");
      for(k=0; k<ntypes; k++)
        *(masses+k) = 1.0;
#ifdef EFILTER
      lower_e_pot = (real *) calloc(ntypes, sizeof(real));
      if (NULL==lower_e_pot)
	  error("Cannot allocate memory for lower_e_pot\n");
      upper_e_pot = (real *) calloc(ntypes, sizeof(real));
      if (NULL==upper_e_pot)
	  error("Cannot allocate memory for upper_e_pot\n");
#endif 
#ifdef NBFILTER
      lower_nb_cut = (int *) calloc(ntypes, sizeof(int));
      if (NULL==lower_nb_cut)
	  error("Cannot allocate memory for lower_nb_cut\n");
      upper_nb_cut = (int *) calloc(ntypes, sizeof(int));
      if (NULL==upper_nb_cut)
	  error("Cannot allocate memory for upper_nb_cut\n");
#endif 
    }
    else if (strcasecmp(token,"starttemp")==0) {
      /* temperature at start of sim. */
      getparam("starttemp",&temperature,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"use_current_temp")==0) {
      /* set imposed temperature to current system temperature */
      use_curr_temp = 1;
    }
#if defined(AND) || defined(NVT) || defined(NPT) || defined(STM) || defined(FRAC)|| defined(FINNIS)
    else if (strcasecmp(token,"endtemp")==0) {
      /* temperature at end of sim. */
      getparam("endtemp",&end_temp,PARAM_REAL,1,1);
    }
#endif
#if defined(STM) || defined(FRAC) || defined(FTG)
    else if (strcasecmp(token,"stadium")==0) {
      getparam("stadium",&stadium,PARAM_REAL,2,2);
    }
    else if (strcasecmp(token,"stadium2")==0) {
      getparam("stadium2",&stadium2,PARAM_REAL,2,2);
    }
    else if (strcasecmp(token,"center")==0) {
      getparam("center",&center,PARAM_REAL,2,2);
    }
#endif
    else if (strcasecmp(token,"cellsize")==0) {
      /* minimal cell diameter */
      getparam("cellsize",&rtmp,PARAM_REAL,1,1);
      cellsz = MAX(cellsz,SQR(rtmp));
    }
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
      getparam("tempintv",&tempintv,PARAM_INT,1,1);
    }
#endif
#if defined(NVT) || defined(NPT) || defined(STM)
    else if (strcasecmp(token,"eta")==0) {
      /* eta variable for NVT or NPT thermostat */
      getparam("eta",&eta,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"tau_eta")==0) {
      /* time constant tau_eta for thermostat */
      getparam("tau_eta",&tau_eta,PARAM_REAL,1,1);
      if (tau_eta == (real)0) {
        error("tau_eta is zero.\n");
      }
      isq_tau_eta = 1.0 / SQR(tau_eta);
    }
    else if (strcasecmp(token,"isq_tau_eta")==0) {
      /* inverse of square of time constant tau_eta for thermostat */
      getparam("isq_tau_eta",&isq_tau_eta,PARAM_REAL,1,1);
      if (isq_tau_eta == (real)0) tau_eta = 0.0;
      else tau_eta = 1.0 / sqrt(isq_tau_eta);
    }
    else if (strcasecmp(token,"inv_tau_eta")==0) {
      /* inverse of time constant tau_eta for thermostat */
      getparam("inv_tau_eta",&isq_tau_eta,PARAM_REAL,1,1);
      if (isq_tau_eta == (real)0) tau_eta = 0.0;
      else tau_eta = 1.0 / isq_tau_eta;
      isq_tau_eta = SQR(isq_tau_eta);
    }
#ifdef UNIAX
    else if (strcasecmp(token,"eta_rot")==0) {
      /* eta variable of rotational motion for NVT or NPT thermostat */
      getparam("eta_rot",&eta_rot,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"tau_eta_rot")==0) {
      /* time constant tau_eta for thermostat of rotational motion */
      getparam("tau_eta_rot",&isq_tau_eta_rot,PARAM_REAL,1,1);
      if (isq_tau_eta_rot == (real)0) {
        error("tau_eta_rot is zero.\n");
      }
      isq_tau_eta_rot = 1.0 / SQR(isq_tau_eta_rot);
    }
    else if (strcasecmp(token,"isq_tau_eta_rot")==0) {
      /* squared inverse of time constant for thermostat of rot. motion */
      getparam("isq_tau_eta_rot",&isq_tau_eta_rot,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"inv_tau_eta_rot")==0) {
      /* inverse of time constant for thermostat of rotational motion */
      getparam("inv_tau_eta_rot",&isq_tau_eta_rot,PARAM_REAL,1,1);
      isq_tau_eta_rot = SQR(isq_tau_eta_rot);
    }
#endif
#endif
#if defined(FRAC) || defined(FTG) 
    else if (strcasecmp(token,"strainrate")==0) {
	/* strain rate for crack loading */
	getparam("strainrate",&dotepsilon0,PARAM_REAL,1,1);
	dotepsilon = dotepsilon0;
    }
    else if (strcasecmp(token,"expansionmode")==0) {
	/* strain mode for crack loading */
	getparam("expansionmode",&expansionmode,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"gamma_bar")==0) { 
       /* Damping prefactor gamma_bar */
	getparam("gamma_bar",&gamma_bar,PARAM_REAL,1,1);
    }
  
    else if (strcasecmp(token,"gamma_damp")==0) { /* actual Damping factor */
	getparam("gamma_damp",&gamma_damp,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"dampingmode")==0) {
	/* damping mode for stadium geometry */
	getparam("dampingmode",&dampingmode,PARAM_INT,1,1);
    }
  
#endif
#ifdef FTG
  else if (strcasecmp(token,"delta_ftg")==0) {
      /* time constant delta for local temperature control  */
      getparam("delta_ftg",&delta_ftg,PARAM_REAL,1,1);
    } 
  else if (strcasecmp(token,"gamma_min")==0) { 
       /* minimal damping prefactor gamma_bar */
	getparam("gamma_min",&gamma_min,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"Tleft")==0) {
      /* damping mode for stadium geometry */
      getparam("Tleft",&Tleft,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"Tright")==0) {
      /* damping mode for stadium geometry */
      getparam("Tright",&Tright,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"nslices")==0) {
      /* nuber of slices*/
      getparam("nslices",&nslices,PARAM_INT,1,1);

      ninslice=(int*) malloc(nslices*sizeof(int));
      if (NULL==ninslice)
	error("Cannot allocate memory for ninslice vector\n");
      for(k=0; k<nslices; k++)
       *(ninslice+k) = 0;

      gamma_ftg=(real*) malloc(nslices*sizeof(real));
      if (NULL==gamma_ftg)
	error("Cannot allocate memory for gamma_ftg vector\n");
      for(k=0; k<nslices; k++)
	*(gamma_ftg+k) = 0.0;
      
      E_kin_ftg=(real*) malloc(nslices*sizeof(real));
      if (NULL==E_kin_ftg)
	error("Cannot allocate memory for E_kin_ftg vector\n");
      for(k=0; k<nslices; k++)
	*(E_kin_ftg+k) = 0.0;
    }
    else if (strcasecmp(token,"gamma_ftg")==0) {
      /* actual Damping factor for each slice */
      /* format: slice gamma_ftg  read in a temp. vektor */
      getparam("gamma_ftg",&tempvek,PARAM_REAL,2,2);
      if (tempvek.x>nslices-1)
	error("actual Damping factorfor non existing slice\n");
      *(gamma_ftg + (int)tempvek.x) = tempvek.y;
    }
    else if (strcasecmp(token,"nslices_Left")==0) {
      /* nuber of slices with Tleft*/
      getparam("nslices_Left",&nslices_Left,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"nslices_Right")==0) {
      /* nuber of slices with Right*/
      getparam("nslices_Right",&nslices_Right,PARAM_INT,1,1);
    }
#endif 
#ifdef FINNIS
  else if (strcasecmp(token,"delta_finnis")==0) {
      /* time constant delta for local temperature control  */
      getparam("delta_finnis",&delta_finnis,PARAM_REAL,1,1);
    } 
  else if (strcasecmp(token,"zeta_0")==0) {
      /* time constant delta for local temperature control  */
      getparam("zeta_0",&zeta_0,PARAM_REAL,1,1);
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
#ifdef SLLOD
    else if (strcasecmp(token,"shear_rate")==0) {
      /* shear strength, corresponds to xy-like entries in strain tensor */
      getparam("shear_rate",&shear_rate,PARAM_REAL,DIM,DIM);
    }
#ifndef TWOD
    else if (strcasecmp(token,"shear_rate2")==0) {
      /* shear strength, corresponds to yx-entry in strain tensor */
      getparam("shear_rate2",&shear_rate2,PARAM_REAL,DIM,DIM);
    }
#endif
#endif
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
      getparam("shear_factor",&shear_factor,PARAM_REAL,2,2);
    }
    else if (strcasecmp(token,"lindef_interval")==0) {
      /* period of linear deformation intervals */
      getparam("lindef_interval",&lindef_interval,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"deform_size")==0) { 
      /* scale factor for deformation */
      getparam("deform_size",&deform_size,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"lindef_x")==0) {
      /* first row of deformation matrix */
      getparam("lindef_x",&lindef_x,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"lindef_y")==0) {
      /* second row of deformation matrix */
      getparam("lindef_y",&lindef_y,PARAM_REAL,DIM,DIM);
    }
#ifndef TWOD
    else if (strcasecmp(token,"lindef_z")==0) {
      /* third row of deformation matrix */
      getparam("lindef_z",&lindef_z,PARAM_REAL,DIM,DIM);
    }
#endif
#endif
#if defined(DEFORM)
 else if (strcasecmp(token,"annealsteps")==0) {
      /* max nr of steps between shears */
      getparam("annealsteps",&annealsteps,PARAM_INT,1,1);
    }
 else if (strcasecmp(token,"ekin_threshold")==0) {
      /* shear epsilon criterium, see imd_shear_new.c */
      getparam("ekin_threshold",&ekin_threshold,PARAM_REAL,1,1);
    }
#endif
#if defined(GLOK)
 else if (strcasecmp(token,"glok_annealsteps")==0) {
      /* number of annealing steps */
      getparam(token,&glok_annealsteps,PARAM_INT,1,1);
    }
 else if (strcasecmp(token,"glok_ekin_threshold")==0) {
      /* threshold for ekin */
      getparam(token,&glok_ekin_threshold,PARAM_REAL,1,1);
    }
#endif
#ifdef DEFORM
    else if (strcasecmp(token,"max_deform_int")==0) {
      /* max nr of steps between shears */
      getparam("max_deform_int",&max_deform_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"fnorm_threshold")==0) {
      /* criterium to contiune deformation */
      getparam("fnorm_threshold",&fnorm_threshold,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"deform_size")==0) { 
      /* scale factor for deformation */
      getparam("deform_size",&deform_size,PARAM_REAL,1,1);
    }
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
    else if (strcasecmp(token,"deform_shear")==0) {
      /* deform shear for virtual types */
      /* format: type shear.x shear.y (shear.z) read in a temp. vektor */
      getparam("deform_shear",&tempshift,PARAM_REAL,DIM+1,DIM+1);
      if (tempshift.x>vtypes-1)
       error("Shear defined for non existing virtual atom type\n");
      shear.x = tempshift.y;
      shear.y = tempshift.z;
#ifndef TWOD
      shear.z = tempshift.z2;
#endif
      *(deform_shear+(int)(tempshift.x)) = shear; 
      *(shear_def+(int)(tempshift.x))    = 1;
    }
    else if (strcasecmp(token,"deform_base")==0) {
      /* deform base for virtual types */
      /* format: type shear.x shear.y (shear.z) read in a temp. vektor */
      getparam("deform_shear",&tempshift,PARAM_REAL,DIM+1,DIM+1);
      if (tempshift.x>vtypes-1)
       error("Shear base defined for non existing virtual atom type\n");
      base.x = tempshift.y;
      base.y = tempshift.z;
#ifndef TWOD
      base.z = tempshift.z2;
#endif
      *(deform_base+(int)(tempshift.x)) = base;
    }
#endif /* DEFORM */
#ifdef CG
    else if (strcasecmp(token,"cg_threshold")==0) {
      /* criterium to contiune relaxation */
      getparam("cg_threshold",&cg_threshold,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"cg_maxsteps")==0) {
      /* max steps of relaxation */
      getparam("cg_maxsteps",&cg_maxsteps,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"linmin_maxsteps")==0) {
      /* max steps to find min in one direction */
      getparam("linmin_maxsteps",&linmin_maxsteps,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"linmin_tol")==0) {
      /* tolerance to stop min search in one direction */
      getparam("linmin_tol",&linmin_tol,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"linmin_dmax")==0) {
      /* max. length of trial step in 1d minimum search */
      getparam("linmin_dmax",&linmin_dmax,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"annealsteps")==0) {
      /* max nr of steps between shears */
      getparam("annealsteps",&annealsteps,PARAM_INT,1,1);
    }
#endif /* CG */
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
    else if (strcasecmp(token,"msqd_ntypes")==0) {
      /* write msqd for real types */
      getparam("msqd_ntypes",&msqd_ntypes,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"msqd_vtypes")==0) {
      /* write msqd for virtual types */
      getparam("msqd_vtypes",&msqd_vtypes,PARAM_INT,1,1);
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
     else if (strcasecmp(token, "press_int")==0){
      /*number of steps between press. writes  */
      getparam("press_int", &press_int, PARAM_INT, 1,1);
    }
#endif
#ifdef DISLOC
    else if (strcasecmp(token,"dem_int")==0) {
      /* number of steps between picture writes */
      getparam("dem_int",&dem_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dsp_int")==0) {
      /* number of steps between picture writes */
      getparam("dsp_int",&dsp_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"update_ort_ref")==0) {
      /* step number to compute ort_ref */
      getparam("update_ort_ref",&up_ort_ref,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"min_dpot")==0) {
      /* minimum Epot difference */
      getparam("min_dpot",&min_dpot,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"min_dsp2")==0) {
      /* minimum square displacement in .dsp files */
      getparam(token,&min_dsp2,PARAM_REAL,1,1);
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
    else if (strcasecmp(token,"avpos_start")==0) {
      /* step at which coordinate addition begins */
      getparam("avpos_start",&avpos_start,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"avpos_end")==0) {
      /* step at which coordinate addition ends */
      getparam("avpos_end",&avpos_end,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"avpos_res")==0) {
      /* number of steps between coordinate addition */
      getparam("avpos_res",&avpos_res,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"avpos_int")==0) {
      /* number of steps between average position writes */
      getparam("avpos_int",&avpos_int,PARAM_INT,1,1);
    }
#endif
#ifdef FORCE
    else if (strcasecmp(token,"force_int")==0) {
      /* number of steps between average position writes */
      getparam("force_int",&force_int,PARAM_INT,1,1);
    }
#endif
#ifdef ATDIST
    else if (strcasecmp(token,"atdist_dim")==0) {
      /* dimension of atoms distribution array */
      getparam(token,&atdist_dim,PARAM_INT,DIM,DIM);
    }
    else if (strcasecmp(token,"atdist_int")==0) {
      /* interval between atoms distribution updates */
      getparam(token,&atdist_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"atdist_pos_int")==0) {
      /* interval between atom position writes */
      getparam(token,&atdist_pos_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"atdist_start")==0) {
      /* step when recording atoms distribution is started */
      getparam(token,&atdist_start,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"atdist_end")==0) {
      /* step when recording atoms distribution is stopped */
      getparam(token,&atdist_end,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"atdist_ll")==0) {
      /* lower left  corner of atoms distribution */
      getparam(token,&atdist_ll,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"atdist_ur")==0) {
      /* upper right corner of atoms distribution */
      getparam(token,&atdist_ur,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"atdist_per_ll")==0) {
      /* lower left of periodic extension */
      getparam(token,&atdist_per_ll,PARAM_INT,DIM,DIM);
    }
    else if (strcasecmp(token,"atdist_per_ur")==0) {
      /* upper right of periodic extension */
      getparam(token,&atdist_per_ur,PARAM_INT,DIM,DIM);
    }
    else if (strcasecmp(token,"atdist_phi")==0) {
      /* rotation angle around z-axis */
      getparam(token,&atdist_phi,PARAM_REAL,1,1);
      atdist_phi *= 8 * atan(1.0);
    }
#endif
#ifdef DIFFPAT
    else if (strcasecmp(token,"diffpat_dim")==0) {
      /* dimension of atoms distribution array */
      getparam(token,&diffpat_dim,PARAM_INT,DIM,DIM);
    }
    else if (strcasecmp(token,"diffpat_int")==0) {
      /* interval between diffraction pattern updates */
      getparam(token,&diffpat_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"diffpat_start")==0) {
      /* step when diffraction pattern recording is started */
      getparam(token,&diffpat_start,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"diffpat_end")==0) {
      /* step when diffraction pattern recording is stopped */
      getparam(token,&diffpat_end,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"diffpat_ur")==0) {
      /* upper right corner of atoms distribution */
      getparam(token,&diffpat_ur,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"diffpat_ll")==0) {
      /* lower left corner of atoms distribution */
      getparam(token,&diffpat_ll,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"diffpat_weight")==0) {
      /* scattering strength of different atom types */
      real w[10];
      if (ntypes==0) error("specify parameter ntypes before diffpat_weight");
      getparam(token,w,PARAM_REAL,ntypes,ntypes);
      for (i=0; i<ntypes; i++) diffpat_weight[i] = (float) w[i];
    }
#endif
#ifdef ORDPAR
    else if (strcasecmp(token,"op_rcut")==0) {
      /* cutoff radius for order parameter */
      getparam("op_rcut",op_r2_cut,PARAM_REAL,4,4);
      op_r2_cut[0][0] = SQR(op_r2_cut[0][0]);
      op_r2_cut[0][1] = SQR(op_r2_cut[0][1]);
      op_r2_cut[1][0] = SQR(op_r2_cut[1][0]);
      op_r2_cut[1][1] = SQR(op_r2_cut[1][1]);
    }   
    else if (strcasecmp(token,"op_weight")==0) {
      /* weights for order parameter */
      getparam("op_weight",op_weight,PARAM_REAL,4,4);
    }
#endif
#ifdef USE_SOCKETS
    else if (strcasecmp(token,"socket_mode")==0) {
      /* socket mode: client or server */
      getparam(token,tmpstr,PARAM_STR,1,255);
      if (strcasecmp(tmpstr,"client")==0) {
	server_socket = 0;
      }
      else if (strcasecmp(tmpstr,"server")==0) {
	server_socket = 1;
      }
      else {
        char msg[255];
        sprintf(msg,"****** Unknown socket mode %s ignored ******",tmpstr);
        warning(msg);
      }
    }
    else if (strcasecmp(token,"socket_int")==0) {
      getparam("socket_int",&socket_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"display_host")==0) {
      getparam("display_host",display_host,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"server_port")==0) { int tmp;
      getparam(token,&tmp,PARAM_INT,1,1);
      server_port = tmp; /* conversion to unsigned short */
    }
    else if (strcasecmp(token,"client_port")==0) { int tmp;
      getparam(token,&tmp,PARAM_INT,1,1);
      client_port = tmp; /* conversion to unsigned short */
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
      getparam("tau_xi",&isq_tau_xi,PARAM_REAL,1,1);
      if (isq_tau_xi == (real)0) {
        error("tau_xi is zero.\n");
      }
      isq_tau_xi = 1.0 / SQR(isq_tau_xi);
    }
    else if (strcasecmp(token,"isq_tau_xi")==0) {
      /* inverse of square of time constant tau_xi for NPT thermostat */
      getparam("isq_tau_xi",&isq_tau_xi,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"inv_tau_xi")==0) {
      /* inverse of time constant tau_xi for NPT thermostat */
      getparam("inv_tau_xi",&isq_tau_xi,PARAM_REAL,1,1);
      isq_tau_xi = SQR(isq_tau_xi);
    }
    else if (strcasecmp(token,"cell_size_tol")==0) {
      /* rel. tolerance for volume rescaling during NPT sim. */
      getparam("cell_size_tol",&cell_size_tolerance,PARAM_REAL,1,1);
    }
#endif
#ifdef EAM2
    else if (strcasecmp(token,"core_potential_file")==0) {
      /* EAM2:Filename for the tabulated Core-Core Potential (r^2) */
      getparam("core_potential_file",potfilename,PARAM_STR,1,255);
      have_potfile = 1;
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
#ifdef MEAM
    else if (strcasecmp(token,"core_potential_file")==0) {
      /* MEAM:Filename for the tabulated Core-Core Potential (r^2) */
      getparam("core_potential_file",potfilename,PARAM_STR,1,255);
      have_potfile = 1;
    }
    else if (strcasecmp(token,"embedding_energy_file")==0) {
      /* MEAM:Filename for the tabulated Embedding Enery(rho_h) */
      getparam("embedding_energy_file",meam_emb_E_filename,PARAM_STR,1,255);
      have_embed_potfile = 1;
    }
    else if (strcasecmp(token,"el_density_file")==0) {
      /* MEAM:Filename for the tabulated electron density */
      getparam("el_density_file",meam_eldensity_filename,PARAM_STR,1,255);
      have_eldensity_file = 1;
    }
    else if (strcasecmp(token,"meam_t_average")==0) {
      getparam(token, &meam_t_average, PARAM_INT, 1, 1);
    }
    else if (strcasecmp(token,"meam_t1")==0) {
      if (ntypes==0) error("specify parameter ntypes before meam_t1");
      getparam(token, meam_t1, PARAM_REAL, ntypes, ntypes);
    }
    else if (strcasecmp(token,"meam_t2")==0) {
      if (ntypes==0) error("specify parameter ntypes before meam_t2");
      getparam(token, meam_t2, PARAM_REAL, ntypes, ntypes);
    }
    else if (strcasecmp(token,"meam_t3")==0) {
      if (ntypes==0) error("specify parameter ntypes before meam_t3");
      getparam(token, meam_t3, PARAM_REAL, ntypes, ntypes);
    }
    else if (strcasecmp(token,"meam_f0")==0) {
      if (ntypes==0) error("specify parameter ntypes before meam_f0");
      getparam(token, meam_f0, PARAM_REAL, ntypes, ntypes);
    }
    else if (strcasecmp(token,"meam_r0")==0) {
      if (ntypes==0) error("specify parameter ntypes before meam_r0");
      getparam(token, meam_r0, PARAM_REAL, ntypes, ntypes);
    }
    else if (strcasecmp(token,"meam_beta0")==0) {
      if (ntypes==0) error("specify parameter ntypes before meam_beta0");
      getparam(token, meam_beta0, PARAM_REAL, ntypes, ntypes);
    }
    else if (strcasecmp(token,"meam_beta1")==0) {
      if (ntypes==0) error("specify parameter ntypes before meam_beta1");
      getparam(token, meam_beta1, PARAM_REAL, ntypes, ntypes);
    }
    else if (strcasecmp(token,"meam_beta2")==0) {
      if (ntypes==0) error("specify parameter ntypes before meam_beta2");
      getparam(token, meam_beta2, PARAM_REAL, ntypes, ntypes);
    }
    else if (strcasecmp(token,"meam_beta3")==0) {
      if (ntypes==0) error("specify parameter ntypes before meam_beta3");
      getparam(token, meam_beta3, PARAM_REAL, ntypes, ntypes);
    }
    else if (strcasecmp(token,"meam_rcut")==0) {
      if (ntypes==0) error("specify parameter ntypes before meam_rcut");
      getparam(token, meam_rcut_lin, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"meam_deltar")==0) {
      if (ntypes==0) error("specify parameter ntypes before meam_deltar");
      getparam(token, meam_deltar_lin, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"meam_cmin")==0) {
      if (ntypes==0) error("specify parameter ntypes before meam_cmin");
      getparam(token, meam_cmin_lin, PARAM_REAL, 1, ntypetriples);
    }
    else if (strcasecmp(token,"meam_cmax")==0) {
      if (ntypes==0) error("specify parameter ntypes before meam_cmax");
      getparam(token, meam_cmax_lin, PARAM_REAL, 1, ntypetriples);
    }
    else if (strcasecmp(token,"meam_a")==0) {
      if (ntypes==0) error("specify parameter ntypes before meam_a");
      getparam(token, meam_a, PARAM_REAL, 1, ntypes);
      have_pre_embed_pot = 1;
    }
    else if (strcasecmp(token,"meam_e")==0) {
      if (ntypes==0) error("specify parameter ntypes before meam_e");
      getparam(token, meam_e, PARAM_REAL, 1, ntypes);
    }
    else if (strcasecmp(token,"meam_rho0")==0) {
      if (ntypes==0) error("specify parameter ntypes before meam_rho0");
      getparam(token, meam_rho0, PARAM_REAL, 1, ntypes);
    }
#endif
    else if (strcasecmp(token,"debug_potential")==0) {
      /* write out interpolated potential */
      getparam(token, &debug_potential, PARAM_INT, 1, 1);
    }
    else if (strcasecmp(token,"debug_pot_res")==0) {
      /* resolution of test interpolation */
      getparam(token, &debug_pot_res, PARAM_INT, 1, 1);
    }
#ifdef PAIR
    /* analytically defined potentials */
    else if (strcasecmp(token,"r_cut")==0) {
      if (ntypes==0) error("specify parameter ntypes before r_cut");
      getparam(token, r_cut_lin, PARAM_REAL, ntypepairs, ntypepairs);
      have_pre_pot = 1;
    }
    else if (strcasecmp(token,"r_begin")==0) {
      if (ntypes==0) error("specify parameter ntypes before r_begin");
      getparam(token, r_begin, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"pot_res")==0) {     
      if (ntypes==0) error("specify parameter ntypes before pot_res");
      getparam(token, pot_res, PARAM_REAL, ntypepairs, ntypepairs);
    }
    /* Lennard-Jones */
    else if (strcasecmp(token,"lj_epsilon")==0) {
      if (ntypes==0) error("specify parameter ntypes before lj_epsilon");
      getparam(token ,lj_epsilon_lin, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"lj_sigma")==0) {
      if (ntypes==0) error("specify parameter ntypes before lj_sigma");
      getparam(token, lj_sigma_lin, PARAM_REAL, ntypepairs, ntypepairs);
    }
    /* Morse */
    else if (strcasecmp(token,"morse_epsilon")==0) {
      if (ntypes==0) error("specify parameter ntypes before morse_epsilon");
      getparam(token, morse_epsilon_lin, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"morse_sigma")==0) {
      if (ntypes==0) error("specify parameter ntypes before morse_sigma");
      getparam(token, morse_sigma_lin, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"morse_alpha")==0) {
      if (ntypes==0) error("specify parameter ntypes before morse_alpha");
      getparam(token, morse_alpha_lin, PARAM_REAL, ntypepairs, ntypepairs);
    }
    /* Buckingham */
    else if (strcasecmp(token,"buck_a")==0) {
      if (ntypes==0) error("specify parameter ntypes before buck_a");
      getparam(token, buck_a_lin, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"buck_c")==0) {
      if (ntypes==0) error("specify parameter ntypes before buck_c");
      getparam(token, buck_c_lin, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"buck_sigma")==0) {
      if (ntypes==0) error("specify parameter ntypes before buck_sigma");
      getparam(token, buck_sigma_lin, PARAM_REAL, ntypepairs, ntypepairs);
    }
    /* harmonic potential for shell model */
    else if (strcasecmp(token,"spring_const")==0) {
      if (ntypes==0) error("specify parameter ntypes before spring_const");
      getparam(token, spring_const, PARAM_REAL, 
               ntypepairs-ntypes, ntypepairs-ntypes);
    }
#endif
#ifdef COVALENT
    else if (strcasecmp(token,"neigh_len")==0) {
      /* number of neighbors */
      getparam(token, &neigh_len, PARAM_INT, 1, 1);
    }
#endif
#ifdef TTBP
    else if (strcasecmp(token,"ttbp_constant")==0) {
      /* force constant (radians); type 0 */
      if (ntypes==0) error("specify parameter ntypes before ttbp_constant");
      getparam(token, ttbp_constant, PARAM_REAL, ntypes, ntypes);
    }
    else if (strcasecmp(token,"ttbp_sp")==0) {
      /* hybridization of the element type */
      if (ntypes==0) error("specify parameter ntypes before ttbp_sp");
      getparam(token, ttbp_sp, PARAM_REAL, ntypes, ntypes);
    }
    else if (strcasecmp(token,"ttbp_potfile")==0) {
      /* filename for ttbp potential data */
      getparam(token, ttbp_potfilename, PARAM_STR, 1, 255);
      have_potfile = 1;
    }
#endif
#ifdef STIWEB
    else if (strcasecmp(token,"stiweb_a")==0) {
      if (ntypes==0) error("specify parameter ntypes before stiweb_a");
      getparam(token, stiweb_a, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"stiweb_b")==0) {
      if (ntypes==0) error("specify parameter ntypes before stiweb_b");
      getparam(token, stiweb_b, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"stiweb_p")==0) {
      if (ntypes==0) error("specify parameter ntypes before stiweb_p");
      getparam(token, stiweb_p, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"stiweb_q")==0) {
      if (ntypes==0) error("specify parameter ntypes before stiweb_q");
      getparam(token, stiweb_q, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"stiweb_a1")==0) {
      if (ntypes==0) error("specify parameter ntypes before stiweb_a1");
      getparam(token, stiweb_a1, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"stiweb_de")==0) {
      if (ntypes==0) error("specify parameter ntypes before stiweb_de");
      getparam(token, stiweb_de, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"stiweb_a2")==0) {
      if (ntypes==0) error("specify parameter ntypes before stiweb_a2");
      getparam(token, stiweb_a2, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"stiweb_ga")==0) {
      if (ntypes==0) error("specify parameter ntypes before stiweb_ga");
      getparam(token, stiweb_ga, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"stiweb_la")==0) {
      if (ntypes==0) error("specify parameter ntypes before stiweb_la");
      getparam(token, stiweb_la, PARAM_REAL, ntypepairs, ntypepairs);
    }
#endif
#ifdef TERSOFF
    /* Parameters for Tersoff potential */
    else if (strcasecmp(token,"ters_r_cut")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_r_cut");
      getparam(token, ters_r_cut, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"ters_r0")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_r0");
      getparam(token, ters_r0, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"ters_a")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_a");
      getparam(token, ters_a, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"ters_b")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_b");
      getparam(token, ters_b, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"ters_la")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_la");
      getparam(token, ters_la, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"ters_mu")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_mu");
      getparam(token, ters_mu, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"ters_chi")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_chi");
      getparam(token, ters_chi, PARAM_REAL, 
               ntypepairs-ntypes, ntypepairs-ntypes);
    }
    else if (strcasecmp(token,"ters_om")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_om");
      getparam(token, ters_om, PARAM_REAL,
               ntypepairs-ntypes, ntypepairs-ntypes);
    }
    else if (strcasecmp(token,"ters_ga")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_ga");
      getparam(token, ters_ga, PARAM_REAL, ntypes, ntypes);
    }
    else if (strcasecmp(token,"ters_n")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_n");
      getparam(token, ters_n, PARAM_REAL, ntypes, ntypes);
    }
    else if (strcasecmp(token,"ters_c")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_c");
      getparam(token, ters_c, PARAM_REAL, ntypes, ntypes);
    }
    else if (strcasecmp(token,"ters_d")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_d");
      getparam(token, ters_d, PARAM_REAL, ntypes, ntypes);
    }
    else if (strcasecmp(token,"ters_h")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_h");
      getparam(token, ters_h, PARAM_REAL, ntypes, ntypes);
    }
#endif
#ifdef KEATING
    /* Parameters for Keating potential */
    else if (strcasecmp(token,"keating_r_cut")==0) {
      if (ntypes==0) error("specify parameter ntypes before keating_r_cut");
      getparam(token, keating_r_cut, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"keating_alpha")==0) {
      if (ntypes==0) error("specify parameter ntypes before keating_alpha");
      getparam(token, keating_alpha, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"keating_d")==0) {
      if (ntypes==0) error("specify parameter ntypes before keating_d");
      getparam(token, keating_d, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"keating_beta")==0) {
      if (ntypes==0) error("specify parameter ntypes before keating_beta");
      getparam(token, keating_beta, PARAM_REAL, 
               ntypes*ntypepairs, ntypes*ntypepairs);
    }
#endif
#ifdef EWALD
    /* Parameters for Ewald summation */
    else if (strcasecmp(token,"charge")==0) {
      if (ntypes==0) error("specify parameter ntypes before charge");
      getparam("charge",charge,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"ew_kappa")==0) {
      getparam("ew_kappa",&ew_kappa,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"ew_kmax")==0) {
      getparam("ew_kmax",&ew_kmax,PARAM_INT,3,3);
    }
    else if (strcasecmp(token,"ew_nmax")==0) {
      getparam("ew_nmax",&ew_nmax,PARAM_INT,1,1);
    }
#endif
#ifdef EPITAX
    /* Parameters for option epitax */
    else if (strcasecmp(token,"epitax_rate")==0) {
      /* rate of creation of particles */
      if (ntypes==0) error("specify parameter ntypes before epitax_rate");
      getparam("epitax_rate",epitax_rate,PARAM_INT,ntypes,ntypes);
    }
    else if (strcasecmp(token,"epitax_type")==0) {
      /* type of particles to be created */
      if (ntypes==0) error("specify parameter ntypes before epitax_type");
      getparam("epitax_type",epitax_type,PARAM_INT,ntypes,ntypes);
    }
    else if (strcasecmp(token,"epitax_mass")==0) {
      /* mass of particles to be created */
      if (ntypes==0) error("specify parameter ntypes before epitax_mass");
      getparam("epitax_mass",epitax_mass,PARAM_REAL,ntypes,ntypes);
    } 
    else if (strcasecmp(token,"epitax_temp")==0) {
      /* temperature of particles to be created */
      if (ntypes==0) error("specify parameter ntypes before epitax_temp");
      getparam("epitax_temp",epitax_temp,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"epitax_cutoff")==0) {
      /* parameter for cutoff */
      getparam("epitax_cutoff",&epitax_cutoff,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"epitax_maxsteps")==0) {
      /* maximal steps in epitax simulation */
      getparam("epitax_maxsteps",&epitax_maxsteps,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"epitax_startstep")==0) {
      /* steps before atom creation starts */
      getparam("epitax_startstep",&epitax_startstep,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"epitax_ctrl")==0) {
      /* parameter for change of integrator  */
      getparam("epitax_ctrl",&epitax_ctrl,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"epitax_height")==0) {
      /* height of beam creation */
      getparam("epitax_height",&epitax_height,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"epitax_speed")==0) {
      /* height of beam creation */
      getparam("epitax_speed",&epitax_speed,PARAM_REAL,1,1);
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
    else if (strcasecmp(token,"use_header")==0) {
	/* shall a header be used */
      getparam("use_header",&use_header,PARAM_INT,1,1);
    }
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
#ifdef PAIR
  if ((have_potfile==0) && (have_pre_pot==0))
    error("You must specify a pair interaction!");
#endif
#if defined(NPT) || defined(NVT) || defined(FRAC) || defined(FINNIS)
  if (temperature == 0) {
    error("starttemp is missing or zero.");
  }
#endif
#if defined(NPT) || defined(NVT) || defined(FRAC) || defined(FINNIS)
  if (end_temp == 0) {
    error("endtemp is missing or zero.");
  }
#endif
#if defined(CORRELATE) || defined(MSQD)
  if (correl_ts == 0) {
    if (eng_int != 0) correl_ts = eng_int;
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
#ifdef FTG
	if (nslices < 2){
	  error ("nslices is missing or less than 2.");
	}
	if (Tleft == 0 ){
	  error ("Tleft is missing or zero.");
	}
	if (Tright == 0 ){
	  error ("Tright is missing or zero.");
	}
#endif
#ifdef MPI
        {
#ifdef TWOD
          int want_cpus = cpu_dim.x * cpu_dim.y;
#else
          int want_cpus = cpu_dim.x * cpu_dim.y * cpu_dim.z;
#endif
          if ( want_cpus != num_cpus) calc_cpu_dim();
          if ((want_cpus != num_cpus) && (want_cpus != 1)) 
            warning("cpu_dim incompatible with available CPUs, using default");
        }
#endif
#ifdef USE_SOCKETS
  if ((!server_socket) && (display_host[0]=='\0')) {
    error("display_host name or IP address missing.");
  }
#endif
#ifdef UNIAX
  if (uniax_r_cut == 0) {
    error("uniax_r_cut is missing or zero.");
  }
#endif

#if defined(FRAC) || defined(FTG) 
  if (stadium2.x==0 && stadium2.y==0 ){
    stadium2.x = box_x.x/2.0;
    stadium2.y = box_y.y/2.0;
  }
#endif
#ifdef AVPOS
  fprintf(stderr, "%d %d\n", avpos_start, imdrestart*checkpt_int);
  if (avpos_start <= imdrestart*checkpt_int)
    avpos_start = imdrestart*checkpt_int+1; /* do not ask me why +1 ;-) */
  /* Default initialisation of end time */ 
  if (0==avpos_end) avpos_end = steps_max;
#endif

#ifdef ATDIST
  if (0==atdist_end) atdist_end = steps_max;
#endif

#ifdef GLOKDEFORM
  if ((fnorm_threshold==0.0) && (max_deform_int == 0))
  {
       error("You have to define fnorm_threshold and max_deform_int ");
  }
#endif
#ifdef CG
  if ((cg_threshold==0.0) || (cg_maxsteps==0))
  {
       error("You have to define cg_threshold and cg_maxsteps  ");
  }
  if ((linmin_maxsteps==0) || (linmin_tol==0.0) || (linmin_dmax==0.0))
  {
       error("You have to parameters for the linmin search  ");
  }
#endif
}

/*****************************************************************
*
*  read command line on master process
*
******************************************************************/

void read_command_line(int argc,char **argv)
{
  if (0==myid) { 
    /* check for restart, process options */
    strcpy(progname,argv[0]);
    while ((argc > 1) && (argv[1][0] =='-')) {
      switch (argv[1][1]) {
        /* r - restart */
        case 'r':
          if (argv[1][2]=='\0') {
            if (NULL != argv[2]) {
              imdrestart = atoi(argv[2]);
              --argc;
              ++argv;
            }
          }
          else imdrestart = atoi(&argv[1][2]);
          break;
        case 'p':
          if (argv[1][2]=='\0') {
            if (NULL != argv[2]) {
              strcpy(paramfilename,argv[2]);
              --argc;
              ++argv;
            }
          }
          else strcpy(paramfilename,&argv[1][2]);
          break;
        default:
          printf("Illegal option %s \n",argv[1]);
          usage();
          exit(-1);
      }
      ++argv;
      --argc;
    }
  }
#ifdef MPI
  /* broadcast everything */
  MPI_Bcast( paramfilename, 255, MPI_CHAR, 0, MPI_COMM_WORLD); 
  MPI_Bcast( progname,      255, MPI_CHAR, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &imdrestart,     1, MPI_INT,  0, MPI_COMM_WORLD); 
#endif
}

/*****************************************************************
*
*  read parameters for the first simulation phase
*
******************************************************************/

void read_parameters(void)
{
  int i;
  str255 fname;
  FILE *testfile;

  if (0==myid) {
    getparamfile(paramfilename,1);
    /* read initial itr-file (if there is any), but keep steps_min value */
    if (0 < strlen(itrfilename)) {
      int tmp_steps = steps_min, tmp_finished = finished;
      getparamfile(itrfilename,1);
      steps_min = tmp_steps;
      finished  = tmp_finished;
    }
    check_parameters_complete();

    /* Get restart parameters if restart */
    if (0 != imdrestart) {
      int tmp = finished;
      sprintf(fname,"%s.%d.itr",outfilename,imdrestart);
      testfile = fopen(fname,"r");
      if (NULL==testfile) { 
        sprintf(fname,"%s.%05d.itr",outfilename,imdrestart);
      } else {
        fclose(testfile);
      }
      sprintf(infilename,"%s.%d.%s",outfilename,imdrestart,"chkpt");
      testfile = fopen(infilename,"r");
      if (NULL==testfile) { 
        sprintf(infilename,"%s.%05d.%s",outfilename,imdrestart,"chkpt");
      } else {
        fclose(testfile);
      }
      printf("Restarting from %s.\n",fname);
      getparamfile(fname,1);
      finished = tmp;
    } else {
      /* if not restart: delete files to which we append */
      sprintf(fname,"%s.eng",              outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.Ekin",      outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.Epot",      outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.press",     outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.presstens", outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.shock",     outfilename); unlink(fname);
      sprintf(fname,"%s.tempdist",         outfilename); unlink(fname);
      sprintf(fname,"%s.msqd",             outfilename); unlink(fname);
    }
  }
#ifdef MPI
  broadcast_params();
#endif

  /* make mapping of virtual types to real types */
  sort_tab = (int *) malloc( vtypes * sizeof(int) );
  if (NULL==sort_tab) error("Cannot allocate sort_tab");
  for (i=0; i<vtypes; i++) sort_tab[i] = i % ntypes;

}


#ifdef MPI

/****************************************************************************
*
*  Broadcast all parameters to other CPU's (MPI only) 
*
*****************************************************************************/

void broadcast_params() {

  int i;

  MPI_Bcast( &finished    , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &ensemble    , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &simulation  , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &loop        , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &seed        , 1, MPI_LONG, 0, MPI_COMM_WORLD); 

  MPI_Bcast( &steps_max   , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &steps_min   , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &checkpt_int , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &eng_int     , 1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &pic_int     , 1, MPI_INT,  0, MPI_COMM_WORLD); 

  MPI_Bcast( &dist_int,              1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &dist_dim,            DIM, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dist_ll,             DIM, REAL,    0, MPI_COMM_WORLD); 
  MPI_Bcast( &dist_ur,             DIM, REAL,    0, MPI_COMM_WORLD); 
  MPI_Bcast( &dist_Epot_flag,        1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &dist_Ekin_flag,        1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &dist_Ekin_long_flag,   1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &dist_Ekin_trans_flag,  1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &dist_press_flag,       1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &dist_presstens_flag,   1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &dist_shock_shear_flag, 1, MPI_INT, 0, MPI_COMM_WORLD); 

#ifdef TWOD
  MPI_Bcast( &pic_scale   , 2, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ecut_kin    , 2, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ecut_pot    , 2, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &pic_res     , 2, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &pic_type    , 1, MPI_INT, 0, MPI_COMM_WORLD); 
#endif
  MPI_Bcast( &pic_ll      , DIM, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &pic_ur      , DIM, REAL, 0, MPI_COMM_WORLD); 
#ifdef CLONE
  MPI_Bcast( &nclones, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
#ifdef SNAPSHOT
  MPI_Bcast( &sscount, 1, MPI_INT, 0, MPI_COMM_WORLD);
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
  MPI_Bcast( &box_unit    ,   1, REAL,     0, MPI_COMM_WORLD); 
  MPI_Bcast( &ntypes      ,   1, MPI_INT,  0, MPI_COMM_WORLD); 
  MPI_Bcast( &ntypepairs  ,   1, MPI_INT,  0, MPI_COMM_WORLD); 
  if (0!=myid) {
    masses=(real*)realloc(masses,ntypes*sizeof(real));
    if (NULL==masses) error("Cannot allocate memory for masses array\n");
  }

#ifdef NBLIST
  MPI_Bcast( &nbl_margin,    1, REAL, 0, MPI_COMM_WORLD);
#endif
#ifdef VEC
  MPI_Bcast( &atoms_per_cpu, 1,  INT, 0, MPI_COMM_WORLD);
#endif
#ifdef EFILTER
  if (0!=myid){
      lower_e_pot = (real *) calloc(ntypes, sizeof(real));
      if (NULL==lower_e_pot)
	  error("Cannot allocate memory for lower_e_pot\n");
      upper_e_pot = (real *) calloc(ntypes, sizeof(real));
      if (NULL==upper_e_pot)
	  error("Cannot allocate memory for upper_e_pot\n");
  }
  MPI_Bcast( lower_e_pot, ntypes,   REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( upper_e_pot, ntypes,   REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ef_checkpt_int, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

#ifdef NBFILTER
  if (0!=myid){
      lower_nb_cut = (int *) calloc(ntypes, sizeof(int));
      if (NULL==lower_nb_cut)
	  error("Cannot allocate memory for lower_nb_cut\n");
      upper_nb_cut = (int *) calloc(ntypes, sizeof(int));
      if (NULL==upper_nb_cut)
	  error("Cannot allocate memory for upper_nb_cut\n");
  }
  MPI_Bcast( lower_nb_cut,     ntypes, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( upper_nb_cut,     ntypes, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &nb_checkpt_int,  1,      MPI_INT, 0, MPI_COMM_WORLD);
#endif

  MPI_Bcast( masses,     ntypes, REAL,     0, MPI_COMM_WORLD); 

  MPI_Bcast( &timestep    ,   1, REAL,     0, MPI_COMM_WORLD); 
  MPI_Bcast( &temperature ,   1, REAL,     0, MPI_COMM_WORLD); 
  MPI_Bcast( &cpu_dim     , DIM, MPI_INT,  0, MPI_COMM_WORLD); 

  MPI_Bcast( &parallel_output, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &parallel_input,  1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( outfilename,            255, MPI_CHAR, 0, MPI_COMM_WORLD); 
  MPI_Bcast( infilename,             255, MPI_CHAR, 0, MPI_COMM_WORLD); 
  MPI_Bcast( reffilename,            255, MPI_CHAR, 0, MPI_COMM_WORLD); 
  MPI_Bcast( potfilename,            255, MPI_CHAR, 0, MPI_COMM_WORLD); 
#ifdef TTBP
  MPI_Bcast( ttbp_potfilename,       255, MPI_CHAR, 0, MPI_COMM_WORLD); 
#endif
#ifdef EAM2
  MPI_Bcast( eam2_emb_E_filename,    255, MPI_CHAR, 0, MPI_COMM_WORLD); 
  MPI_Bcast( eam2_at_rho_filename,   255, MPI_CHAR, 0, MPI_COMM_WORLD); 
#endif
#ifdef MEAM
  MPI_Bcast( meam_emb_E_filename,    255, MPI_CHAR, 0, MPI_COMM_WORLD); 
  MPI_Bcast( meam_eldensity_filename,255, MPI_CHAR, 0, MPI_COMM_WORLD); 
#endif

#if defined(AND) || defined(NVT) || defined(NPT) || defined(STM) || defined(FRAC)
  MPI_Bcast( &end_temp,      1, REAL, 0, MPI_COMM_WORLD); 
#endif
  MPI_Bcast( &cellsz, 1, REAL,     0, MPI_COMM_WORLD); 
  MPI_Bcast( &initsz, 1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &incrsz, 1, MPI_INT,  0, MPI_COMM_WORLD);

#ifdef AND
  MPI_Bcast( &tempintv, 1, MPI_INT, 0, MPI_COMM_WORLD); 
#endif

#if defined(NVT) || defined(NPT) || defined(STM)
  MPI_Bcast( &eta ,         1 , REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &isq_tau_eta , 1 , REAL, 0, MPI_COMM_WORLD); 
#ifdef UNIAX
  MPI_Bcast( &eta_rot ,         1 , REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &isq_tau_eta_rot , 1 , REAL, 0, MPI_COMM_WORLD); 
#endif
#endif

#if defined(STM) || defined(FRAC) || defined(FTG)
  MPI_Bcast( &stadium,          2 , REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &stadium2,         2 , REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &center,           2 , REAL, 0, MPI_COMM_WORLD); 
#endif

#ifdef NPT
  MPI_Bcast( &xi,                DIM, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &isq_tau_xi,          1, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &pressure_ext,      DIM, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &pressure_end,      DIM, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &cell_size_tolerance, 1, REAL, 0, MPI_COMM_WORLD); 
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
  MPI_Bcast( &msqd_ntypes,  1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &msqd_vtypes,  1, MPI_INT, 0, MPI_COMM_WORLD);
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
  MPI_Bcast( &press_int    , 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
#if defined(FRAC) || defined(FTG) 
  MPI_Bcast( &dotepsilon0   , 1, REAL   , 0, MPI_COMM_WORLD); 
  MPI_Bcast( &expansionmode , 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &gamma_bar     , 1, REAL   , 0, MPI_COMM_WORLD); 
  MPI_Bcast( &gamma_min     , 1, REAL   , 0, MPI_COMM_WORLD); 
  MPI_Bcast( &gamma_damp    , 1, REAL   , 0, MPI_COMM_WORLD); 
  MPI_Bcast( &dampingmode   , 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &delta_ftg     , 1, REAL   , 0, MPI_COMM_WORLD); 
#endif

#ifdef FTG
  MPI_Bcast( &Tleft,         1, REAL      , 0, MPI_COMM_WORLD);
  MPI_Bcast( &Tright,        1, REAL      , 0, MPI_COMM_WORLD);
  MPI_Bcast( &nslices,       1, MPI_INT   , 0, MPI_COMM_WORLD); 
  MPI_Bcast( &nslices_Left,  1, MPI_INT   , 0, MPI_COMM_WORLD); 
  MPI_Bcast( &nslices_Right, 1, MPI_INT   , 0, MPI_COMM_WORLD); 

  if (0!=myid){ 
    ninslice  = (int*) malloc(nslices*sizeof(int));
    if (NULL==ninslice)
      error("Cannot allocate memory for ninslice vector on client.\n");
  }                   
  if (0!=myid){ 
    E_kin_ftg = (real*) malloc(nslices*sizeof(real));
    if (NULL==E_kin_ftg) 
      error("Cannot allocate memory for E_kin_ftg vector on client.\n");
  }
  if (0!=myid){
    gamma_ftg = (real*) malloc(nslices*sizeof(real));
    if (NULL==gamma_ftg)
      error("Cannot allocate memory for gamma_ftg vector on client.\n");
    for(i=0;i<nslices;i++) *(gamma_ftg + i) = 0.0;
  }
  MPI_Bcast( gamma_ftg, nslices, REAL     , 0, MPI_COMM_WORLD);
#endif 

#ifdef FINNIS
 MPI_Bcast( &delta_finnis     , 1, REAL   , 0, MPI_COMM_WORLD); 
 MPI_Bcast( &zeta_0           , 1, REAL   , 0, MPI_COMM_WORLD); 
#endif
#if defined(GLOK)
 MPI_Bcast( &glok_annealsteps,    1, MPI_INT, 0, MPI_COMM_WORLD); 
 MPI_Bcast( &glok_ekin_threshold, 1, REAL,    0, MPI_COMM_WORLD); 
#endif
#ifdef DEFORM
 MPI_Bcast( &annealsteps,         1, MPI_INT, 0, MPI_COMM_WORLD); 
 MPI_Bcast( &ekin_threshold,      1, REAL,    0, MPI_COMM_WORLD); 
  MPI_Bcast( &max_deform_int,  1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &deform_size,     1, REAL,    0, MPI_COMM_WORLD); 
  MPI_Bcast( &fnorm_threshold, 1, REAL,    0, MPI_COMM_WORLD); 
  if (0!=myid) deform_shift = (vektor *) malloc( vtypes * DIM * sizeof(real) );
  if (NULL==deform_shift) 
    error("Cannot allocate memory for deform_shift on client."); 
  MPI_Bcast( deform_shift, vtypes * DIM, REAL, 0, MPI_COMM_WORLD);
  if (0!=myid) {
    shear_def = (int *) malloc( vtypes * sizeof(int) );
    if (NULL==shear_def) 
      error("Cannot allocate memory for shear_def on client."); 
    for(i=0; i<vtypes; i++) *(shear_def+i) = 0;
  }
  MPI_Bcast( shear_def, vtypes, MPI_INT, 0, MPI_COMM_WORLD);
  if (0!=myid) deform_shear = (vektor *) malloc( vtypes * DIM * sizeof(real) );
  if (NULL==deform_shear) 
    error("Cannot allocate memory for deform_shear on client."); 
  MPI_Bcast( deform_shear, vtypes * DIM, REAL, 0, MPI_COMM_WORLD);
  if (0!=myid) deform_base = (vektor *) malloc( vtypes * DIM * sizeof(real) );
  if (NULL==deform_base) 
    error("Cannot allocate memory for deform_base on client."); 
  MPI_Bcast( deform_base, vtypes * DIM, REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef HOMDEF
  MPI_Bcast( &hom_interval,    1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &shear_factor,    2, REAL,    0, MPI_COMM_WORLD); 
  MPI_Bcast( &exp_interval,    1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &expansion,     DIM, REAL,    0, MPI_COMM_WORLD); 

  MPI_Bcast( &deform_size,     1, REAL,    0, MPI_COMM_WORLD); 
  MPI_Bcast( &lindef_interval, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &lindef_x,      DIM, REAL,    0, MPI_COMM_WORLD); 
  MPI_Bcast( &lindef_y,      DIM, REAL,    0, MPI_COMM_WORLD); 
#ifndef TWOD
  MPI_Bcast( &lindef_z,      DIM, REAL,    0, MPI_COMM_WORLD); 
#endif
#endif

#ifdef SHOCK
  MPI_Bcast( &shock_strip, 1, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &shock_speed, 1, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &shock_mode,  1, MPI_INT, 0, MPI_COMM_WORLD); 
#endif

#ifdef DISLOC
  MPI_Bcast( &min_dpot,        1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &min_dsp2,        1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dem_int,         1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dsp_int,         1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &calc_Epot_ref,   1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &reset_Epot_step, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &Epot_diff,       1, MPI_INT, 0, MPI_COMM_WORLD); 
#endif

#ifdef AVPOS
  MPI_Bcast( &avpos_start,       1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &avpos_end,         1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &avpos_int,         1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &avpos_res,         1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

#ifdef ATDIST
  error("Option ATDIST is not supported under MPI");
#endif

#ifdef DIFFPAT
  error("Option DIFFPAT is not supported under MPI");
#ifdef TWOD
  error("Option DIFFPAT is not supported in 2D");
#endif
#endif

#ifdef ORDPAR
  MPI_Bcast( &op_r2_cut,       4, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &op_weight,       4, REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef CG
  MPI_Bcast( &cg_threshold,    1, REAL,    0, MPI_COMM_WORLD); 
  MPI_Bcast( &cg_maxsteps,     1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &linmin_maxsteps, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &linmin_tol,      1, REAL,    0, MPI_COMM_WORLD); 
  MPI_Bcast( &linmin_dmax,     1, REAL,    0, MPI_COMM_WORLD); 
#ifndef DEFORM
  MPI_Bcast( &annealsteps,     1, MPI_INT, 0, MPI_COMM_WORLD); 
#endif
#endif

#ifdef USE_SOCKETS
  MPI_Bcast( &socket_int, 1, MPI_INT, 0, MPI_COMM_WORLD); 
#endif

#ifdef UNIAX
  MPI_Bcast( &uniax_r_cut,  1, REAL, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &uniax_r2_cut, 1, REAL, 0, MPI_COMM_WORLD); 
#endif

#ifdef PAIR
  /* analytically defined potentials */
  MPI_Bcast( &have_pre_pot,   1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( &have_potfile,   1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast( r_cut_lin, ntypepairs, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( r_begin,   ntypepairs, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( pot_res,   ntypepairs, REAL, 0, MPI_COMM_WORLD);
  /* Lennard-Jones */
  MPI_Bcast( lj_epsilon_lin, ntypepairs, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( lj_sigma_lin,   ntypepairs, REAL, 0, MPI_COMM_WORLD);
  /* Morse */
  MPI_Bcast( morse_epsilon_lin, ntypepairs, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( morse_sigma_lin,   ntypepairs, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( morse_alpha_lin,   ntypepairs, REAL, 0, MPI_COMM_WORLD);
  /* Buckingham */
  MPI_Bcast( buck_a_lin,     ntypepairs, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( buck_c_lin,     ntypepairs, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( buck_sigma_lin, ntypepairs, REAL, 0, MPI_COMM_WORLD);
  /* harmonic potential for shell model */
  MPI_Bcast( spring_const, ntypepairs-ntypes, REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef COVALENT
  MPI_Bcast( &neigh_len, 1, MPI_INT,  0, MPI_COMM_WORLD);
#endif

#ifdef TTBP
  MPI_Bcast( ttbp_constant,  ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ttbp_sp,        ntypes, REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef SLLOD
  MPI_Bcast(&shear_rate, 2, REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef STIWEB
  MPI_Bcast( stiweb_a,         ntypepairs, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( stiweb_b,         ntypepairs, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( stiweb_p,         ntypepairs, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( stiweb_q,         ntypepairs, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( stiweb_a1,        ntypepairs, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( stiweb_de,        ntypepairs, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( stiweb_a2,        ntypepairs, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( stiweb_ga,        ntypepairs, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( stiweb_la, ntypes*ntypepairs, REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef TERSOFF
  MPI_Bcast( ters_r_cut, ntypepairs,        REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_r0,    ntypepairs,        REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_a,     ntypepairs,        REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_b,     ntypepairs,        REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_la,    ntypepairs,        REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_mu,    ntypepairs,        REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_chi,   ntypepairs-ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_om,    ntypepairs-ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_ga,               ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_n,                ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_c,                ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_d,                ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_h,                ntypes, REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef KEATING
  MPI_Bcast( keating_r_cut,       ntypepairs, REAL, 0,MPI_COMM_WORLD);
  MPI_Bcast( keating_alpha,       ntypepairs, REAL, 0,MPI_COMM_WORLD);
  MPI_Bcast( keating_d,           ntypepairs, REAL, 0,MPI_COMM_WORLD);
  MPI_Bcast( keating_beta, ntypes*ntypepairs, REAL, 0,MPI_COMM_WORLD);
#endif

#ifdef EPITAX
  MPI_Bcast( epitax_rate,     ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( epitax_type,     ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( epitax_mass,     ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( epitax_temp,     ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &epitax_maxsteps,     1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &epitax_startstep,    1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &epitax_ctrl,         1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &epitax_height,       1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &epitax_speed,        1, REAL, 0, MPI_COMM_WORLD);
#endif
  
  MPI_Bcast(&use_header,1, MPI_INT, 0, MPI_COMM_WORLD);

  /* broadcast integrator to other CPU's */

  switch (ensemble) {
    case ENS_NVE:       move_atoms = move_atoms_nve;       break;
    case ENS_MIK:       move_atoms = move_atoms_mik;       break;
    case ENS_NVT:       move_atoms = move_atoms_nvt;       break;
    case ENS_SLLOD:     move_atoms = move_atoms_sllod;     break;
    case ENS_NPT_ISO:   move_atoms = move_atoms_npt_iso;   break;
    case ENS_NPT_AXIAL: move_atoms = move_atoms_npt_axial; break;
    case ENS_FRAC:      move_atoms = move_atoms_frac;      break;
    case ENS_NVX:       move_atoms = move_atoms_nvx;       break;
    case ENS_STM:       move_atoms = move_atoms_stm;       break;  
    case ENS_FTG:       move_atoms = move_atoms_ftg;       break;  
    case ENS_FINNIS:    move_atoms = move_atoms_finnis;    break;  
    case ENS_CG:     /* move_atoms = move_atoms_cg; */     break;  
    default: if (0==myid) error("unknown ensemble in broadcast"); break;
  }
  
}

#endif /* MPI */
