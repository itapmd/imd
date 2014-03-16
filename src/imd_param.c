/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2012 Institute for Theoretical and Applied Physics,
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


/* the following is needed for gettimeofday */
#include <sys/time.h>

#include "imd.h"

#if defined(CBE)
#include "imd_cbe.h"
#endif


typedef enum ParamType {
  PARAM_STR, PARAM_STRPTR, PARAM_CHARPTR,
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
		  PARAM_CHARPTR : String, deklariert als char**
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

static int getparam(char *param_name, void *param, PARAMTYPE ptype,
		    int pnum_min, int pnum_max)
{
  static char errmsg[256];
  char *str;
  int i;
  int numread;

  numread = 0;
  if (ptype == PARAM_STR) {
    str = strtok(NULL," =\t\r\n");
    if (str == NULL) {
      sprintf(errmsg,"parameter for %s missing in line %u\nstring expected",
              param_name,curline);
      error(errmsg);
    }
    else strncpy((char *)param,str,pnum_max);
    numread++;
  }
  else if (ptype == PARAM_STRPTR) {
    str = strtok(NULL," =\t\r\n");
    if (str == NULL) {
      sprintf(errmsg,"parameter for %s missing in line %u\nstring expected",
              param_name,curline);
      error(errmsg);
    }
    else *((char**)param) = strdup(str);
    numread++;
  }
  else if (ptype == PARAM_CHARPTR) {
    i = 0;
    str = strtok(NULL," =\t\r\n");
    while (str != NULL && i < ntypes) {
      strncpy(((char **)param)[i],str,strlen(str));
      str = strtok(NULL," =\t\r\n");
      i++;
    };
    if (i < ntypes) {
      sprintf(errmsg, "parameter for %s missing in line %d\nlist of %d strings expected", param_name, curline, ntypes);
      error(errmsg);
    }
  }
  else if (ptype == PARAM_INT) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," =\t\r\n");
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
      if ((str = strtok(NULL," =\t\r\n")) != NULL) {
        ((int*)param)[i] = atoi(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_INT_COPY) {
    int ival = 0;
    for (i=0; i<pnum_max; i++) {
      str = strtok(NULL," =\t\r\n");
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
      str = strtok(NULL," =\t\r\n");
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
      if ((str = strtok(NULL," =\t\r\n")) != NULL) {
        ((integer*)param)[i] = atoi(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_INTEGER_COPY) {
    int ival = 0;
    for (i=0; i<pnum_max; i++) {
      str = strtok(NULL," =\t\r\n");
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
      str = strtok(NULL," =\t\r\n");
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
      if ((str = strtok(NULL," =\t\r\n")) != NULL) {
        ((real*)param)[i] = atof(str);
        numread++;
      }
      else break;
    }
  }
  else if (ptype == PARAM_REAL_COPY) {
    real rval = 0;
    for (i=0; i<pnum_max; i++) {
      str = strtok(NULL," =\t\r\n");
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

int getparamfile(char *paramfname, int phase)
{
  FILE *pf;
  char buffer[1024];
  char *token;
  char *res;
  str255 tmpstr;
  int tmp, finished = 0;
  real rtmp;
#ifdef EXTPOT
  real rtmp4[4];
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
  vektor shear, base;
  int k;
  int i;

  vektor3d tempv3d;

  real norm_bend_axis;
  ivektor2d tmp_ivec2d INIT(nullvektor2d);

  curline = 0;
  pf = fopen(paramfname,"r");
  if (NULL == pf) {
    error_str("Could not open parameter file %s", paramfname);
  }

  /* set the random number generator seed to the */
  /* negative of the current time in seconds */
  /* this will be superseded by a fixed value from the parameter file */
  {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    seed = (long) -tv.tv_sec;
  }

  do {
    res=fgets(buffer,1024,pf);
    if (NULL == res) { finished=1; break; }; /* probably EOF reached */
    curline++;
    /* delete comments */
    res = strchr(buffer, '#');
    if (res) *res = '\0';
    token = strtok(buffer," =\t\r\n");
    if (NULL == token) continue; /* skip blank lines */

    if (strcasecmp(token,"simulation")==0) {
      /* get number of the simulation phase */
      getparam(token,&tmp,PARAM_INT,1,1);
      if (phase < tmp) break;
    }
#ifdef DEBUG
    else if (strcasecmp(token,"force_celldim_divisor")==0) {
	getparam(token,&force_celldim_divisor,PARAM_INT,3,3);
	  }
#endif
    else if (strcasecmp(token,"maxwalltime")==0) {
      /* maximal walltime limit */
      getparam(token,&maxwalltime,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"watch_int")==0) {
      /* interval for checking write file */
      getparam(token,&watch_int,PARAM_INT,1,1);
      stop_int = watch_int;
    }
    else if (strcasecmp(token,"stop_int")==0) {
      /* interval for checking stop file */
      getparam(token,&stop_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"hyper_threads")==0) {
      /* number of hyperthreads per CPU */
      getparam(token,&hyper_threads,PARAM_INT,1,1);
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
    else if (strcasecmp(token,"do_maxwell")==0) {
      /* force temperature initialization */
      getparam(token,&do_maxwell,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"box_from_header")==0) {
      /* read box from config file */
      getparam(token,&box_from_header,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"coordname")==0) {
      /* file name for atom coordinate input data */
      getparam("coordname",infilename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"itrname")==0) {
      /* file name for initial itr-file */
      getparam(token,itrfilename,PARAM_STR,1,255);
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
    else if (strcasecmp(token,"ensemble")==0) {
      /* ensemble */
      getparam(token,tmpstr,PARAM_STR,1,255);
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
      else if (strcasecmp(tmpstr,"glok")==0) {
        ensemble = ENS_GLOK;
        move_atoms = move_atoms_nve;
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
#ifdef CG
      else if (strcasecmp(tmpstr,"cg")==0) {
        ensemble = ENS_CG;
        move_atoms = move_atoms_cg;
      }
#endif
      else if (strcasecmp(tmpstr,"ttm")==0) {
        ensemble = ENS_TTM;
	move_atoms = move_atoms_ttm;
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
    else if (strcasecmp(token,"flush_int")==0) {
      /* interval for flushing .eng file */
      getparam(token,&flush_int,PARAM_INT,1,1);
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
    else if (strcasecmp(token,"dist_pressoff_flag")==0) {
      /* write press dist? */
      getparam(token,&dist_pressoff_flag,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dist_presstens_flag")==0) {
      /* write pressoff dist? */
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
    else if (strcasecmp(token,"dist_Ekin_comp_flag")==0) {
      /* write difference Ekin dist? */
      getparam(token,&dist_Ekin_comp_flag,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dist_shock_shear_flag")==0) {
      /* write shock shear dist? */
      getparam(token,&dist_shock_shear_flag,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dist_shear_aniso_flag")==0) {
      /* write shear aniso dist? */
      getparam(token,&dist_shear_aniso_flag,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dist_dens_flag")==0) {
      /* write density dist? */
      getparam(token,&dist_dens_flag,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dist_vxavg_flag")==0) {
      /* write average sample velocity? */
      getparam(token,&dist_vxavg_flag,PARAM_INT,1,1);
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
    else if (strcasecmp(token,"nbl_size")==0) {
      /* size of neighbor list */
      getparam(token,&nbl_size,PARAM_REAL,1,1);
    }
#endif
#ifdef NEB
    else if (strcasecmp(token,"neb_nrep")==0) {
      /* number of NEB replicas */
      getparam(token,&neb_nrep,PARAM_INT,1,1);
      if (0==myrank)
	{
        if (num_cpus != neb_nrep)
          error("We need exactly neb_nrep MPI processes");
        if (neb_nrep>NEB_MAXNREP)
          error("Too many images for NEB");
	}
    }
    else if (strcasecmp(token,"neb_eng_int")==0) {
      /* interval of NEB energy writes */
      getparam(token,&neb_eng_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"neb_cineb_start")==0) {
      /* when to change to CINEB */
      getparam(token,&neb_cineb_start,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"neb_climbing_image")==0) {
      /* which image should be climbing */
      getparam(token,&neb_climbing_image,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"neb_vark_start")==0) {
      /* when to change to variable ks */
      getparam(token,&neb_vark_start,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"neb_k")==0) {
      /* spring constant of NEB */
      getparam(token,&neb_k,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"neb_maxmove")==0) {
      /* constrain relaxation steps */
      getparam(token,&neb_maxmove,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"neb_kmax")==0) {
      /* if >0 variable springs are used with  max. spring constant */
      getparam(token,&neb_kmax,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"neb_kmin")==0) {
      /* if >0 variable springs are used with  max. spring constant */
      getparam(token,&neb_kmin,PARAM_REAL,1,1);
    }
#endif
#ifdef BBOOST

    else if (strcasecmp(token,"bb_tot_bV")==0) {
      /* magnitude of boost potential, unit according to potential */
        getparam(token,&bb_tot_bV,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"bb_p1_2")==0) {
      /* curvature controller of the boost potential */
        getparam(token,&p1_2,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"bb_relaxsteps_max")==0) {
      /* max number of steps for the relaxation part of the bond boos method */
        getparam(token,&bb_relaxsteps_max,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"bb_shdn_max")==0) {
      /* max number of steps for the shutdown part of the bond boos method */
        getparam(token,&bb_shdn_max,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"bb_under_max")==0) {
      /* max number of steps for the under boosting part of the bond boos method */
        getparam(token,&bb_under_max,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"bb_epscrit")==0) {
      if (ntypes==0)
          error("specify parameter ntypes before bb_epscrit");
      /* critical fraction of bond length to consider a bond broken */
      /* format: type1 type2 epscrit */
      getparam(token,&tempv3d,PARAM_REAL,3,3);
      if (tempv3d.x>ntypes-1 || tempv3d.y>ntypes-1 )
          error("bb_epscrit defined for non existing type of bond\n");
      bb_epscrit[(int)(tempv3d.x)][(int)(tempv3d.y)] = tempv3d.z;
    }
    else if (strcasecmp(token,"bb_rcut")==0) {
      /* the cut off for the range of the bond boos method */
        getparam(token,&bb_rcut,PARAM_REAL,1,1);
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
#ifdef NNBR
    else if (strcasecmp(token,"nb_rcut")==0) {
      /* cutoff radius for coordination number */
      if (ntypes==0) error("specify parameter ntypes before nb_rcut");
      getparam(token,nb_r2_cut,PARAM_REAL,ntypes*ntypes,ntypes*ntypes);
      for (k=0; k<ntypes*ntypes; k++) nb_r2_cut[k] = SQR(nb_r2_cut[k]);
    }
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
    else if (strcasecmp(token,"total_types")==0) {
      /* TOTAL nuber of atom types: ntypes + virtual types */
      int vt, init = (vtypes==0);
      getparam(token,&vt,PARAM_INT,1,1);
      if (init) {
        vtypes = vt;
      }
      else if (vt != vtypes) {
        error("total_types must be constant during a simulation");
      }
      /* do some allocations and initialisations */
      if (init) {
        restrictions = (vektor *) malloc( vtypes * sizeof(vektor) );
        if (NULL==restrictions)
	  error("Cannot allocate memory for restriction vectors\n");
        for (k=0; k<vtypes; k++)
          restrictions[k] = einsv;
#ifdef FBC
        /* Allocation & Initialisation of fbc_forces */
        fbc_forces = (vektor *) malloc( vtypes * sizeof(vektor) );
        if (NULL==fbc_forces)
	  error("Cannot allocate memory for fbc_forces\n");
        for (k=0; k<vtypes; k++)
          fbc_forces[k] = nullv;
        fbc_beginforces = (vektor *) malloc( vtypes * sizeof(vektor) );
        if (NULL==fbc_beginforces)
          error("Cannot allocate memory for fbc_beginforces\n");
        for (k=0; k<vtypes; k++)
          fbc_beginforces[k] = nullv;
#ifdef RELAX
        fbc_dforces = (vektor *) malloc( vtypes * sizeof(vektor) );
        if (NULL==fbc_dforces)
	  error("Cannot allocate memory for fbc_dforces\n");
        for (k=0; k<vtypes; k++)
          fbc_dforces[k] = nullv;
#else
        fbc_endforces = (vektor *) malloc( vtypes * sizeof(vektor) );
        if (NULL==fbc_endforces)
	  error("Cannot allocate memory for fbc_endforces\n");
        for (k=0; k<vtypes; k++)
          fbc_endforces[k] = nullv;
#endif
         fbc_df = (vektor *) malloc( vtypes * sizeof(vektor) );
        if (NULL==fbc_df)
	  error("Cannot allocate memory for fbc_df\n");
        for (k=0; k<vtypes; k++)
          fbc_df[k] = nullv;
#endif /*FBC*/
#ifdef BEND
        /* Allocation & Initialisation of fbc_bforces */
        fbc_bforces = (vektor *) malloc( vtypes * sizeof(vektor) );
        if (NULL==fbc_bforces)
	  error("Cannot allocate memory for fbc_bforces\n");
        for (k=0; k<vtypes; k++)
            fbc_bforces[k] = nullv;
        /* Allocation & Initialisation of bend_forces */
        bend_forces = (vektor *) malloc( vtypes * sizeof(vektor) );
        if (NULL==bend_forces)
	  error("Cannot allocate memory for bend_forces\n");
        for (k=0; k<vtypes; k++)
          bend_forces[k] = nullv;
        fbc_beginbforces = (vektor *) malloc( vtypes * sizeof(vektor) );
        if (NULL==fbc_beginbforces)
          error("Cannot allocate memory for fbc_beginbforces\n");
        for (k=0; k<vtypes; k++)
          fbc_beginbforces[k] = nullv;
#ifdef RELAX
        fbc_bdforces = (vektor *) malloc( vtypes * sizeof(vektor) );
        if (NULL==fbc_bdforces)
	  error("Cannot allocate memory for fbc_bdforces\n");
        for (k=0; k<vtypes; k++)
          fbc_bdforces[k] = nullv;
#else
        fbc_endbforces = (vektor *) malloc( vtypes * sizeof(vektor) );
        if (NULL==fbc_endbforces)
	  error("Cannot allocate memory for fbc_endbforces\n");
        for (k=0; k<vtypes; k++)
          fbc_endbforces[k] = nullv;
#endif
         fbc_bdf = (vektor *) malloc( vtypes * sizeof(vektor) );
        if (NULL==fbc_bdf)
	  error("Cannot allocate memory for fbc_bdf\n");
        for (k=0; k<vtypes; k++)
          fbc_bdf[k] = nullv;
#endif /* BEND */
#ifdef RIGID
        /* Allocation & Initialization of superatom */
        superatom = (int *) malloc( vtypes * sizeof(int) );
        if (NULL==superatom)
	  error("Cannot allocate memory for superatom vector\n");
        for (k=0; k<vtypes; k++)
	  superatom[k] = -1;

        /* Allocation of superforce */
        superforce = (vektor *) malloc( vtypes * sizeof(vektor) );
        if (NULL==superforce)
          error("Cannot allocate memory for superforce vector\n");

        /* Allocation & Initialization of superrestrictions */
        superrestrictions = (vektor *) malloc( vtypes * sizeof(vektor) );
        if (NULL==superrestrictions)
	  error("Cannot allocate memory for superrestriction vectors\n");
        for (k=0; k<vtypes; k++)
          superrestrictions[k] = nullv;
#endif
#ifdef DEFORM
        /* Allocation & Initialisation of deform_shift */
        deform_shift = (vektor *) malloc( vtypes * sizeof(vektor) );
        if (NULL==deform_shift)
          error("Cannot allocate memory for deform_shift\n");
        for (k=0; k<vtypes; k++)
          deform_shift[k] = nullv;

        /* Allocation & Initialisation of shear_def */
        shear_def = (int *) malloc( vtypes * sizeof(int) );
        if (NULL==shear_def)
          error("Cannot allocate memory for shear_def\n");
        for (k=0; k<vtypes; k++)
          shear_def[k] = 0;

        /* Allocation & Initialisation of deform_shear */
        deform_shear = (vektor *) malloc( vtypes * sizeof(vektor) );
        if (NULL==deform_shear)
          error("Cannot allocate memory for deform_shear\n");
        for (k=0; k<vtypes; k++)
          deform_shear[k] = nullv;

        /* Allocation & Initialisation of deform_base */
        deform_base = (vektor *) malloc( vtypes * sizeof(vektor) );
        if (NULL==deform_base)
          error("Cannot allocate memory for deform_base\n");
        for (k=0; k<vtypes; k++)
          deform_base[k] = nullv;
#endif
      }
    }
#ifdef RIGID
    else if (strcasecmp(token,"rigid")==0) {
      int count, tmp, rigidv[15];
      if (vtypes==0)
        error("specify parameter total_types before rigid");
      /* virtual types forming superparticle */
      i = getparam(token,rigidv,PARAM_INT,1+DIM,vtypes+DIM);
      /* determine number of types in superparticle */
      count = i - DIM;
      /* construct superatom vector */
      tmp = superatom[rigidv[0]];
      for (i=0; i<count; i++) {
        if ( rigidv[i] > vtypes - 1 )
          error("Atom type in superparticle does not exist\n");
        if ( superatom[rigidv[i]] != tmp )
          error("Intersecting superparticles\n");
        if (tmp < 0)
          superatom[rigidv[i]] = nsuperatoms;
      }
      if (tmp < 0) {
        superrestrictions[nsuperatoms].x = rigidv[count  ];
        superrestrictions[nsuperatoms].y = rigidv[count+1];
#ifndef TWOD
        superrestrictions[nsuperatoms].z = rigidv[count+2];
#endif
        nsuperatoms++;
      }
    }
#endif

#ifdef RELAX
    else if (strcasecmp(token,"ekin_threshold")==0) {
      /* threshold for sufficient relaxation */
      getparam(token,&ekin_threshold,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"fnorm_threshold")==0) {
      /* threshold for sufficient relaxation */
      getparam(token,&fnorm_threshold,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"f_max_threshold")==0) {
      /* threshold for sufficient relaxation */
      getparam(token,&f_max_threshold,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"delta_epot_threshold")==0) {
      /* threshold for sufficient relaxation */
      getparam(token,&delta_epot_threshold,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"sscount")==0) {
      /* snapshot counter, for restarting */
      getparam(token,&sscount,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"nfc")==0) {
      /* nfc counter, for restart */
      getparam(token,&nfc,PARAM_INT,1,1);
    }
#endif

#ifdef FBC
    else if (strcasecmp(token,"extra_startforce")==0) {
      if (vtypes==0)
        error("specify parameter total_types before extra_startforce");
      /* extra force for virtual types */
      /* format: type force.x force.y (force.z) */
      getparam(token,&tempforce,PARAM_REAL,DIM+1,DIM+1);
      if (tempforce.x>vtypes-1)
        error("Force defined for non existing virtual atom type\n");
      force.x = tempforce.y;
      force.y = tempforce.z;
#ifndef TWOD
      force.z = tempforce.z2;
#endif
      fbc_beginforces[(int)(tempforce.x)] = force;
      fbc_forces     [(int)(tempforce.x)] = force;
    }
#ifdef RELAX
    else if (strcasecmp(token,"fbc_ekin_threshold")==0) {
      /* epsilon criterium to increment extra force*/
      getparam(token,&ekin_threshold,PARAM_REAL,1,1);
      warning("Parameter fbc_ekin_threshold replaced by ekin_threshold");
    }
    else if (strcasecmp(token,"max_fbc_int")==0) {
      /* max nr of steps between fbc increments */
      getparam(token,&max_fbc_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"fbc_waitsteps")==0) {
      /* max nr of steps between fbc increments */
      getparam(token,&max_fbc_int,PARAM_INT,1,1);
      warning("Parameter fbc_waitsteps replaced by max_fbc_int");
    }
    else if (strcasecmp(token,"extra_dforce")==0) {
      if (vtypes==0)
        error("specify parameter total_types before extra_dforce");
      /* extra force increment for virtual types */
      /* format: type force.x force.y (force.z)  */
      getparam(token,&tempforce,PARAM_REAL,DIM+1,DIM+1);
      if (tempforce.x>vtypes-1)
        error("Force increment defined for non existing virtual atom type\n");
      force.x = tempforce.y;
      force.y = tempforce.z;
#ifndef TWOD
      force.z = tempforce.z2;
#endif
      fbc_dforces[(int)(tempforce.x)] = force;
    }
#else
    else if (strcasecmp(token,"extra_endforce")==0) {
      if (vtypes==0)
        error("specify parameter total_types before extra_endforce");
      /* extra force for virtual types */
      /* format: type force.x force.y (force.z) */
      getparam(token,&tempforce,PARAM_REAL,DIM+1,DIM+1);
      if (tempforce.x>vtypes-1)
        error("Force defined for non existing virtual atom type\n");
      force.x = tempforce.y;
      force.y = tempforce.z;
#ifndef TWOD
      force.z = tempforce.z2;
#endif
      fbc_endforces[(int)(tempforce.x)] = force;
    }
#endif
#endif /* FBC */
#ifdef BEND
    else if (strcasecmp(token,"extra_startbforce")==0) {
      if (vtypes==0)
        error("specify parameter total_types before extra_bstartforce");
      /* extra force for virtual types */
      /* format: type force.x force.y (force.z) */
      getparam(token,&tempforce,PARAM_REAL,DIM+1,DIM+1);
      if (tempforce.x>vtypes-1)
        error("Force defined for non existing virtual atom type\n");
      force.x = tempforce.y;
      force.y = tempforce.z;
#ifndef TWOD
      force.z = tempforce.z2;
#endif
      fbc_beginbforces[(int)(tempforce.x)] = force;
      fbc_bforces     [(int)(tempforce.x)] = force;
    }
#ifdef RELAX
    else if (strcasecmp(token,"max_bfbc_int")==0) {
      /* max nr of steps between fbc increments */
      getparam(token,&max_bfbc_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"bfbc_waitsteps")==0) {
      /* max nr of steps between fbc increments */
      getparam(token,&max_bfbc_int,PARAM_INT,1,1);
      warning("Parameter bfbc_waitsteps replaced by max_fbc_int");
    }
    else if (strcasecmp(token,"extra_bdforce")==0) {
      if (vtypes==0)
        error("specify parameter total_types before extra_bdforce");
      /* extra force increment for virtual types */
      /* format: type force.x force.y (force.z)  */
      getparam(token,&tempforce,PARAM_REAL,DIM+1,DIM+1);
      if (tempforce.x>vtypes-1)
        error("Force increment defined for non existing virtual atom type\n");
      force.x = tempforce.y;
      force.y = tempforce.z;
#ifndef TWOD
      force.z = tempforce.z2;
#endif
      fbc_bdforces[(int)(tempforce.x)] = force;
    }
#else
    else if (strcasecmp(token,"extra_endbforce")==0) {
      if (vtypes==0)
        error("specify parameter total_types before extra_endforce");
      /* extra force for virtual types */
      /* format: type force.x force.y (force.z) */
      getparam(token,&tempforce,PARAM_REAL,DIM+1,DIM+1);
      if (tempforce.x>vtypes-1)
        error("Force defined for non existing virtual atom type\n");
      force.x = tempforce.y;
      force.y = tempforce.z;
#ifndef TWOD
      force.z = tempforce.z2;
#endif
      fbc_endbforces[(int)(tempforce.x)] = force;
    }
#endif
#endif /* BEND */


#ifdef FLAGEDATOMS
 else if (strcasecmp(token,"flagedatomstype")==0) {
        getparam(token,&flagedatomstype,PARAM_INT,1,1);
    }
#endif

#ifdef ZAPP
    else if (strcasecmp(token,"zapp_threshold")==0) {
        getparam(token,&zapp_threshold,PARAM_REAL,1,1);
    }
#endif
#ifdef BEND
    else if (strcasecmp(token,"bend_nmoments")==0) {
        /* nr of bending moments */
      getparam(token,&bend_nmoments,PARAM_INT,1,1);
      if (vtypes<bend_nmoments*2)
          error("need vtypes>=2*bend_nmoments") ;
      /* now do some allocations and initialisations */
      bend_axis = (vektor *) malloc( bend_nmoments * sizeof(vektor) );
      if (NULL==bend_axis)
          error("Cannot allocate memory for bend_axis vectors\n");
      for (k=0; k<bend_nmoments; k++)
      {bend_axis[k].x=0.0;bend_axis[k].y=0.0;bend_axis[k].z=0.0;}
      bend_origin = (vektor *) malloc( bend_nmoments * sizeof(vektor) );
      if (NULL==bend_origin)
          error("Cannot allocate memory for bend_origin vectors\n");
      for (k=0; k<bend_nmoments; k++)
      {bend_origin[k].x=0.0;bend_origin[k].y=0.0;bend_origin[k].z=0.0;}
      bend_cog = (vektor *) malloc( bend_nmoments * sizeof(vektor) );
      if (NULL==bend_cog)
          error("Cannot allocate memory for bend_cog vectors\n");
      for (k=0; k<bend_nmoments; k++)
      {bend_cog[k].x=0.0;bend_cog[k].y=0.0;bend_cog[k].z=0.0;}

      bend_vec = (vektor *) malloc( bend_nmoments * sizeof(vektor) );
      if (NULL==bend_vec)
          error("Cannot allocate memory for bend_vec vectors\n");
      for (k=0; k<bend_nmoments; k++)
          {bend_vec[k].x=0.0;bend_vec[k].y=0.0;bend_vec[k].z=0.0;}

    }
    else if (strcasecmp(token,"bend_axis")==0) {
      /* axis of all bending moments */
      getparam("bend_axis",&tempforce,PARAM_REAL,DIM+1,DIM+1);
      if (tempforce.x>bend_nmoments)
          error("Bend axis defined for non existing bending moment\n");
      force.x = tempforce.y;
      force.y = tempforce.z;
#ifndef TWOD
      force.z = tempforce.z2;
#endif
      if(SPROD(force,force)!=1)
      {
          warning("axis of bending moment not of unit length");
          norm_bend_axis = sqrt(SPROD(force,force));
          force.x /= norm_bend_axis;
          force.y /= norm_bend_axis;
          force.z /= norm_bend_axis;
      }
      bend_axis[(int)(tempforce.x)] = force;
    }
    else if (strcasecmp(token,"bend_vtype_of_origin")==0) {
        /* vtype of origin of bendingmoment i */
        /* format: bend_vtype_of_origin bendingmoment vtype */
      getparam("bend_vtype_of_origin",&tmp_ivec2d,PARAM_INT,1,2);
      bend_vtype_of_origin[tmp_ivec2d.x]=tmp_ivec2d.y;
    }
    else if (strcasecmp(token,"bend_vtype_of_force")==0) {
        /* vtype of atoms to apply FBC of bendingmoment i */
        /* format: bend_vtype_of_force bendingmoment vtype */
      getparam("bend_vtype_of_force",&tmp_ivec2d,PARAM_INT,1,2);
      bend_vtype_of_force[tmp_ivec2d.x]=tmp_ivec2d.y;
    }


#endif /* BEND */

    else if (strcasecmp(token,"restrictionvector")==0) {
      if (vtypes==0)
        error("specify parameter total_types before restrictionvector");
      /* restrictions for virtual types */
      /* format: type  1 1 (1) (=all directions ok) */
      getparam(token,&tempvek,PARAM_REAL,DIM+1,DIM+1);
      if (tempvek.x>vtypes-1)
        error("Restriction defined for non existing virtual atom type\n");
      vek.x = tempvek.y;
      vek.y = tempvek.z;
#ifndef TWOD
      vek.z = tempvek.z2;
#endif
      restrictions[(int)(tempvek.x)] = vek;
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
      getparam(token,&box_param,PARAM_INT,DIM,DIM);
    }
    else if (strcasecmp(token,"size_per_cpu")==0) {
      /* box parameters are given per CPU */
      getparam(token,&size_per_cpu,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"box_unit")==0) {
      /* lattice parameter for generated structures */
      getparam(token,&box_unit,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"masses")==0) {
      /* masses for generated structures */
      if (ntypes==0)
        error("specify parameter ntypes before parameter masses");
      getparam(token,masses,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"types")==0) {
      /* types for generated structures */
      if (ntypes==0)
        error("specify parameter ntypes before parameter types");
      getparam(token,gtypes,PARAM_INT,1,ntypes);
    }
    else if (strcasecmp(token,"timestep")==0) {
      /* size of timestep (in MD units) */
      getparam(token,&timestep,PARAM_REAL,1,1);
#ifdef ADAPTGLOK
      starttimestep = timestep;
#endif
    }
    else if (strcasecmp(token,"ntypes")==0) {
      /* number of atom types */
      int nt, init = (ntypes==0);
      getparam(token,&nt,PARAM_INT,1,1);
      if (init) {
        ntypes = nt;
      } else if (nt != ntypes) {
        error("ntypes must be constant during a simulation");
      }
      if (init) {
#ifdef MONO
        if (ntypes!=1)
        error("this executable is for monoatomic systems only!");
#endif
        ntypepairs = ((ntypes+1)*ntypes)/2;
        ntypetriples = ntypes * ntypepairs;
#if defined(TERSOFF2) || defined(TERSOFFMOD2) || defined(BRENNER)
        nvalues = ntypepairs;
#elif defined(TERSOFF) || defined(TERSOFFMOD) || defined(BRENNER)
        nvalues = ntypes;
#endif
        /* array of masses for generated structures */
        masses = (real *) malloc( ntypes * sizeof(real) );
        if (NULL==masses)
          error("Cannot allocate memory for masses array\n");
        for (k=0; k<ntypes; k++)
          masses[k] = 1.0;
        /* array of types for generated structures */
        gtypes = (int *) malloc( ntypes * sizeof(int) );
        if (NULL==gtypes)
          error("Cannot allocate memory for types array\n");
        for (k=0; k<ntypes; k++)
          gtypes[k] = k;
#ifdef EFILTER
        lower_e_pot = (real *) calloc(ntypes, sizeof(real));
        if (NULL==lower_e_pot)
          error("Cannot allocate memory for lower_e_pot\n");
        upper_e_pot = (real *) calloc(ntypes, sizeof(real));
        if (NULL==upper_e_pot)
          error("Cannot allocate memory for upper_e_pot\n");
#endif
#ifdef NNBR
        lower_nb_cut = (int *) calloc(ntypes, sizeof(int));
        if (NULL==lower_nb_cut)
          error("Cannot allocate memory for lower_nb_cut\n");
        upper_nb_cut = (int *) calloc(ntypes, sizeof(int));
        if (NULL==upper_nb_cut)
          error("Cannot allocate memory for upper_nb_cut\n");
        nb_r2_cut = (real *) calloc(ntypes*ntypes, sizeof(real));
        if (NULL==nb_r2_cut)
          error("Cannot allocate memory for nb_r2_cut");
#endif
#ifdef ORDPAR
        op_r2_cut = (real *) calloc(ntypes*ntypes, sizeof(real));
        if (NULL==op_r2_cut)
          error("Cannot allocate memory for op_r2_cut");
        op_weight = (real *) calloc(ntypes*ntypes, sizeof(real));
        if (NULL==op_weight)
          error("Cannot allocate memory for op_weight");
#endif
#ifdef KIM
	if (NULL == kim_el_names)
	  kim_el_names = (char **)calloc(ntypes, sizeof(char *));
	if (NULL == kim_el_names)
	  error("Cannot allocate memory for kim_el_names\n");
	for (i = 0;i<ntypes;i++) {
	  /* allocate 12 characters for each species - should be more than enough */
	  kim_el_names[i] = (char *)calloc(12, sizeof(char));
	  if (NULL == kim_el_names[i])
	    error("Cannot allocate memory for kim_el_names\n");
	  sprintf(kim_el_names[i],'\0');
	}
#endif
      }
    }
    else if (strcasecmp(token,"starttemp")==0) {
      /* temperature at start of sim. */
      getparam("starttemp",&temperature,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"use_current_temp")==0) {
      /* set imposed temperature to current system temperature */
      use_curr_temp = 1;
    }
#ifdef TEMPCONTROL
    else if (strcasecmp(token,"endtemp")==0) {
      /* temperature at end of sim. */
      getparam(token,&end_temp,PARAM_REAL,1,1);
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
#ifdef DAMP
    /* keep to old stadium convention  */
    else if (strcasecmp(token,"stadium")==0) {
      getparam("stadium",&stadium,PARAM_REAL,3,3);
    }
    else if (strcasecmp(token,"stadium2")==0) {
      getparam("stadium2",&stadium2,PARAM_REAL,3,3);
    }
    else if (strcasecmp(token,"center")==0) {
      getparam("center",&center,PARAM_REAL,3,3);
    }
    else if (strcasecmp(token,"damptemp")==0) { /* actual Damping factor */
        getparam("damptemp",&damptemp,PARAM_REAL,1,1);
    }
    /* Damping prefactor */
    else if (strcasecmp(token,"zeta_0")==0) {
        getparam("zeta_0",&zeta_0,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"delta_finnis")==0) { /* actual Damping factor */
        getparam("delta_finnis",&delta_finnis,PARAM_REAL,1,1);
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
    else if (strcasecmp(token,"outbuf_size")==0) {
      /* output buffer size in MB */
      getparam(token,&outbuf_size,PARAM_INT,1,1);
      outbuf_size *= 1048576;
    }
    else if (strcasecmp(token,"inbuf_size")==0) {
      /* total input buffer size in MB */
      getparam(token,&inbuf_size,PARAM_INT,1,1);
      inbuf_size *= 1048576;
    }
    else if (strcasecmp(token,"dist_chunk_size")==0) {
      /* size of MPI reduction in mega-floats */
      getparam(token,&dist_chunk_size,PARAM_INT,1,1);
      dist_chunk_size *= 1048576;
    }
#ifdef AND
    else if (strcasecmp(token,"tempintv")==0) {
      /* temperature interval */
      getparam("tempintv",&tempintv,PARAM_INT,1,1);
    }
#endif
#ifdef BER
    else if (strcasecmp(token,"tau_berendsen")==0) {
      /* temperature interval */
      getparam("tau_berendsen",&tauber,PARAM_REAL,1,1);
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
    else if (strcasecmp(token,"uniax_inert")==0) {
      /* moment of inertia */
      getparam(token,&uniax_inert,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"uniax_sig")==0) {
      /* nearest neighbor distances of potential in the three directions */
      getparam(token,&uniax_sig,PARAM_REAL,3,3);
      if (uniax_sig.x != uniax_sig.y)
        error("UNIAX molecules must be uniaxial!");
    }
    else if (strcasecmp(token,"uniax_eps")==0) {
      /* depth of potential in the three directions */
      getparam(token,&uniax_eps,PARAM_REAL,3,3);
      if (uniax_eps.x != uniax_eps.y)
         error("UNIAX molecules must be uniaxial!");
    }
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
      int ns, init = (nslices==0);
      getparam(token,&ns,PARAM_INT,1,1);
      if (init) {
        nslices = ns;
      } else if ( ns != nslices) {
        error("nslices must be constant during a simulation");
      }
      if (init) {
        ninslice = (int *) malloc(nslices*sizeof(int));
        if (NULL==ninslice)
          error("Cannot allocate memory for ninslice vector\n");
        for (k=0; k<nslices; k++)
          ninslice[k] = 0;
        gamma_ftg = (real *) malloc(nslices*sizeof(real));
        if (NULL==gamma_ftg)
          error("Cannot allocate memory for gamma_ftg vector\n");
        for (k=0; k<nslices; k++)
	  gamma_ftg[k] = 0.0;
        E_kin_ftg = (real*) malloc(nslices*sizeof(real));
        if (NULL==E_kin_ftg)
          error("Cannot allocate memory for E_kin_ftg vector\n");
        for (k=0; k<nslices; k++)
	  E_kin_ftg[k] = 0.0;
      }
    }
    else if (strcasecmp(token,"gamma_ftg")==0) {
      /* actual Damping factor for each slice */
      if (nslices==0)
        error("specify parameter nslices before gamma_ftg");
      /* format: slice gamma_ftg */
      getparam(token,&tempvek,PARAM_REAL,2,2);
      if (tempvek.x>nslices-1)
	error("actual Damping factorfor non existing slice\n");
      gamma_ftg[(int)(tempvek.x)] = tempvek.y;
    }
    else if (strcasecmp(token,"nslices_Left")==0) {
      /* nuber of slices with Tleft */
      getparam(token,&nslices_Left,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"nslices_Right")==0) {
      /* nuber of slices with Right */
      getparam(token,&nslices_Right,PARAM_INT,1,1);
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
    else if (strcasecmp(token,"nsmear")==0) {
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
#ifdef CYCLE
       else if (strcasecmp(token,"lindef_freq")==0) {
           /* frequency for deformation */
      getparam(token,&lindef_freq,PARAM_REAL,1,1);
    }
#endif
#ifdef HOMDEF
    else if (strcasecmp(token,"lindef_interval")==0) {
      /* period of linear deformation intervals */
      getparam(token,&lindef_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"lindef_size")==0) {
        /* scale factor for deformation */
        /* in case of CYCLE this is the strain amplitude */
      getparam(token,&lindef_size,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"lindef_x")==0) {
      /* first row of deformation matrix */
      getparam(token,&lindef_x,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"lindef_y")==0) {
      /* second row of deformation matrix */
      getparam(token,&lindef_y,PARAM_REAL,DIM,DIM);
    }
#ifndef TWOD
    else if (strcasecmp(token,"lindef_z")==0) {
      /* third row of deformation matrix */
      getparam(token,&lindef_z,PARAM_REAL,DIM,DIM);
    }
#endif
    else if (strcasecmp(token,"shear_module")==0) {
      /* estimate of shear module */
      getparam(token,&shear_module,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"bulk_module")==0) {
      /* estimate of bulk module */
      getparam(token,&bulk_module,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"relax_rate")==0) {
      /* pressure relaxation rate */
      getparam(token,&relax_rate,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"relax_mode")==0) {
      /* pressure relaxation mode */
      getparam(token,tmpstr,PARAM_STR,1,255);
      if      (strcasecmp(tmpstr,"full" )==0) relax_mode = RELAX_FULL;
      else if (strcasecmp(tmpstr,"axial")==0) relax_mode = RELAX_AXIAL;
      else if (strcasecmp(tmpstr,"iso"  )==0) relax_mode = RELAX_ISO;
      else    error_str("Unknown relax_mode %s", tmpstr);
    }
#endif
#if defined(HOMDEF) || defined(NPT_axial)
    else if (strcasecmp(token,"relax_dirs")==0) {
      /* box lengths which should be relaxed */
      getparam("relax_dirs",&relax_dirs,PARAM_INT,DIM,DIM);
    }
#endif
#ifdef GLOK
    else if (strcasecmp(token,"glok_ekin_threshold")==0) {
      /* threshold for ekin */
      getparam(token,&glok_ekin_threshold,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"fire_ekin_threshold")==0) {
      /* threshold for ekin */
      getparam(token,&glok_ekin_threshold,PARAM_REAL,1,1);
    }

#endif
#ifdef MIX
    else if (strcasecmp(token,"glok_mix")==0) {
      /* factor to turn velocities more parallel to forces */
      getparam(token,&glok_mix,PARAM_REAL,1,1);
    }
     else if (strcasecmp(token,"fire_mix")==0) {
      /* factor to turn velocities more parallel to forces */
      getparam(token,&glok_mix,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"glok_mixdec")==0) {
      /*decrease factor to turn velocities more parallel to forces */
      getparam(token,&glok_mixdec,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"fire_mixdec")==0) {
      /*decrease factor to turn velocities more parallel to forces */
      getparam(token,&glok_mixdec,PARAM_REAL,1,1);
    }
#endif
#ifdef ADAPTGLOK
    /* FIRE = mixedadaptglok, for cosmetic reasons & backward compatibility now all parameters
       can either be called glok_something or fire_something */
    else if (strcasecmp(token,"glok_minsteps")==0) {
      /* minimum of steps before increasing the timestep */
      getparam(token,&glok_minsteps,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"fire_minsteps")==0) {
      /* minimum of steps before increasing the timestep */
      getparam(token,&glok_minsteps,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"min_nPxF")==0) {
      /* minimum gloks before increasing the timestep */
      getparam(token,&min_nPxF,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"glok_fmaxcrit")==0) {
      /* critical max. force component  */
      getparam(token,&glok_fmaxcrit,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"fire_fmaxcrit")==0) {
      /* critical max. force component  */
      getparam(token,&glok_fmaxcrit,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"glok_incfac")==0) {
      /* factor to increase the timestep */
      getparam(token,&glok_incfac,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"fire_incfac")==0) {
      /* factor to increase the timestep */
      getparam(token,&glok_incfac,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"glok_decfac")==0) {
      /* factor to decrease the timestep */
      getparam(token,&glok_decfac,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"fire_decfac")==0) {
      /* factor to decrease the timestep */
      getparam(token,&glok_decfac,PARAM_REAL,1,1);
    }
   else if (strcasecmp(token,"glok_maxtimestep")==0) {
      /* max timestep */
      getparam(token,&glok_maxtimestep,PARAM_REAL,1,1);
    }
   else if (strcasecmp(token,"fire_maxtimestep")==0) {
      /* max timestep */
      getparam(token,&glok_maxtimestep,PARAM_REAL,1,1);
    }
   else if (strcasecmp(token,"glok_int")==0) {
      /* only needed for restarting */
      getparam(token,&glok_int,PARAM_INT,1,1);
    }
   else if (strcasecmp(token,"fire_int")==0) {
      /* only needed for restarting */
      getparam(token,&glok_int,PARAM_INT,1,1);
    }
#endif
#ifdef DEFORM
    else if (strcasecmp(token,"max_deform_int")==0) {
      /* max nr of steps between shears */
      getparam("max_deform_int",&max_deform_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"deform_size")==0) {
      /* scale factor for deformation */
      getparam("deform_size",&deform_size,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"deform_shift")==0) {
      /* deform shift for virtual types */
      /* format: type shift.x shift.y (shift.z) */
      if (vtypes==0)
        error("specify parameter total_types before deform_shift");
      getparam(token,&tempshift,PARAM_REAL,DIM+1,DIM+1);
      if (tempshift.x>vtypes-1)
        error("Shift defined for non existing virtual atom type\n");
      shift.x = tempshift.y;
      shift.y = tempshift.z;
#ifndef TWOD
      shift.z = tempshift.z2;
#endif
      deform_shift[(int)(tempshift.x)] = shift;
    }
    else if (strcasecmp(token,"deform_shear")==0) {
      /* deform shear for virtual types */
      /* format: type shear.x shear.y (shear.z) */
      if (vtypes==0)
        error("specify parameter total_types before deform_shear");
      getparam(token,&tempshift,PARAM_REAL,DIM+1,DIM+1);
      if (tempshift.x>vtypes-1)
        error("Shear defined for non existing virtual atom type\n");
      shear.x = tempshift.y;
      shear.y = tempshift.z;
#ifndef TWOD
      shear.z = tempshift.z2;
#endif
      deform_shear[(int)(tempshift.x)] = shear;
      shear_def   [(int)(tempshift.x)] = 1;
    }
    else if (strcasecmp(token,"deform_base")==0) {
      /* deform base for virtual types */
      /* format: type shear.x shear.y (shear.z) */
      if (vtypes==0)
        error("specify parameter total_types before deform_base");
      getparam(token,&tempshift,PARAM_REAL,DIM+1,DIM+1);
      if (tempshift.x>vtypes-1)
        error("Shear base defined for non existing virtual atom type\n");
      base.x = tempshift.y;
      base.y = tempshift.z;
#ifndef TWOD
      base.z = tempshift.z2;
#endif
      deform_base[(int)(tempshift.x)] = base;
    }
#endif /* DEFORM */
#ifdef CG
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
    else if (strcasecmp(token,"linmin_dmin")==0) {
      /* max. length of trial step in 1d minimum search */
      getparam("linmin_dmin",&linmin_dmin,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"cg_glimit")==0) {
      /* limit in mnbrak */
      getparam("cg_glimit",&cg_glimit,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"cg_zeps")==0) {
      /* in brent */
      getparam("cg_zeps",&cg_zeps,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"cg_fr")==0) {
      /* Fletcher-Reeves mode or not*/
      getparam(token,&cg_fr,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"cg_reset_int")==0) {
      /* interval for resetting cg */
      getparam(token,&cg_reset_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"cg_infolevel")==0) {
      /* cg_infolevel controls verbosity */
      getparam(token,&cg_infolevel,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"cg_mode")==0) {
      /* conjugate gradient mode - at present just the default one */
      getparam(token,tmpstr,PARAM_STR,1,255);
      if (strcasecmp(tmpstr,"cge")==0) {
        cg_mode = CGE;
      }
      /* not implemented yet
      else if (strcasecmp(tmpstr,"cgef")==0) {
        cg_mode = CGEF;
      }
      */
      else error_str("unknown CG mode %s",tmpstr);
    }
#endif /* CG */

#ifdef ACG
      else if (strcasecmp(token,"acg_alpha")==0) {
	/* starting alpha */
	getparam(token,&acg_init_alpha,PARAM_REAL,1,1);
      }
      else if (strcasecmp(token,"acg_incfac")==0) {
	/* increase alpha */
	getparam(token,&acg_incfac,PARAM_REAL,1,1);
      }
      else if (strcasecmp(token,"acg_decfac")==0) {
	/* decrease alpha */
	getparam(token,&acg_decfac,PARAM_REAL,1,1);
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
    else if (strcasecmp(token,"shock_speed_left")==0) {
      /* shock speed (in x dir.) */
      getparam("shock_speed_l",&shock_speed_l,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"shock_speed_right")==0) {
      /* shock speed (in x dir.) */
      getparam("shock_speed_r",&shock_speed_r,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"shock_incr")==0) {
      /* steps to achieve full velocity */
      getparam("shock_incr",&shock_incr,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"shock_mode")==0) {
       /* shock type: plate or half */
       getparam("shock_mode",&shock_mode,PARAM_INT,1,1);
       if (shock_mode > 1) shock_strip = 0;
       /* compatibility with old input files */
       if (shock_mode < 2 && shock_mode > 4) shock_mode = 1;
       /* */
       if (shock_mode == 4 && shock_speed_l ==0) shock_speed_l = shock_speed;
       if (shock_mode == 4 && shock_speed_r ==0) shock_speed_r = shock_speed;
    }
#endif
#ifdef MPI
    else if (strcasecmp(token,"cpu_dim")==0) {
      /* CPU array dimension */
      getparam(token,&cpu_dim,PARAM_INT,DIM,DIM);
    }
    else if (strcasecmp(token,"parallel_output")==0) {
      /* parallel output flag */
      getparam(token,&parallel_output,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"parallel_input")==0) {
      /* parallel input flag */
      getparam(token,&parallel_input,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"msgbuf_size")==0) {
      /* security factor of message buffer size */
      getparam(token,&msgbuf_size,PARAM_REAL,1,1);
    }
#endif
    else if (strcasecmp(token,"binary_output")==0) {
      /* binary output flag */
      getparam(token,&binary_output,PARAM_INT,1,1);
    }
#ifdef CORRELATE
    else if (strcasecmp(token,"correl_rmax")==0) {
      /* dimension of histogram in r domain */
      getparam("correl_rmax",&ncorr_rmax,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"correl_tmax")==0) {
      /* dimension of histogram in t domain */
      getparam("correl_tmax",&ncorr_tmax,PARAM_INT,1,1);
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
    else if (strcasecmp(token,"correl_int")==0) {
      /* repeat interval for correlation */
      getparam("correl_int",&correl_int,PARAM_INT,1,1);
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
#ifdef NMOLDYN
    else if (strcasecmp(token,"nmoldyn_int")==0) {
      /* interval for nmoldyn trajectory writes */
      getparam(token,&nmoldyn_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"nmoldyn_veloc")==0) {
      /* include velocities in nmoldyn trajectory? */
      getparam(token,&nmoldyn_veloc,PARAM_INT,1,1);
    }
#endif
#ifdef DSF
    else if (strcasecmp(token,"dsf_int")==0) {
      /* interval for dsf updates */
      getparam(token,&dsf_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dsf_weight")==0) {
      /* weights for dsf (usually coherent scattering length) */
      if (0==ntypes) error("specify parameter ntypes before dsf_weight");
      dsf_weight = (real *) malloc( ntypes * sizeof(real) );
      if (NULL==dsf_weight) error("cannot allocate dsf_weight");
      getparam(token,dsf_weight,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"dsf_nk")==0) {
      /* number of k-point series */
      getparam(token,&dsf_nkmax,PARAM_INT,1,1);
      dsf_k0   = (int *) malloc( dsf_nkmax * DIM * sizeof(int) );
      dsf_kdir = (int *) malloc( dsf_nkmax * DIM * sizeof(int) );
      dsf_kmax = (int *) malloc( dsf_nkmax       * sizeof(int) );
      if ((NULL==dsf_k0) || (NULL==dsf_kdir) || (NULL==dsf_kmax))
        error("cannot allocate dsf arrays");
    }
    else if (strcasecmp(token,"dsf_k")==0) {
      /* k-point series */
      int i=0, tmp[2*DIM+1];
      if (dsf_nk>=dsf_nkmax)
        error("number of k-point series exceeds dsf_nkmax");
      getparam(token,tmp,PARAM_INT,2*DIM+1,2*DIM+1);
      dsf_k0  [DIM*dsf_nk  ] = tmp[i++];
      dsf_k0  [DIM*dsf_nk+1] = tmp[i++];
#ifndef TWOD
      dsf_k0  [DIM*dsf_nk+2] = tmp[i++];
#endif
      dsf_kdir[DIM*dsf_nk  ] = tmp[i++];
      dsf_kdir[DIM*dsf_nk+1] = tmp[i++];
#ifndef TWOD
      dsf_kdir[DIM*dsf_nk+2] = tmp[i++];
#endif
      dsf_kmax[    dsf_nk  ] = tmp[i++];
      dsf_nk++;
    }
#endif
#if defined(HC) || defined(NVX)
    else if (strcasecmp(token, "hc_int")==0){
      /* number of steps between heat current or profile writes  */
      getparam(token, &hc_int, PARAM_INT, 1,1);
    }
    else if (strcasecmp(token, "hc_start")==0){
      /* start step for heat current or profile measurement  */
      getparam(token, &hc_start, PARAM_INT, 1,1);
    }
#endif
#ifdef HC
    else if (strcasecmp(token, "hc_av_start")==0){
      /* start step for energy averaging  */
      getparam(token, &hc_av_start, PARAM_INT, 1,1);
    }
#endif
#ifdef NVX
    else if (strcasecmp(token, "hc_nlayers")==0){
      /* number of layers */
      getparam(token, &hc_nlayers, PARAM_INT, 1,1);
    }
    else if (strcasecmp(token, "hc_count")==0){
      /* running index of temperature profile */
      getparam(token, &hc_count, PARAM_INT, 1,1);
    }
    else if (strcasecmp(token, "hc_heatcurr")==0){
      /* induced heat current density */
      getparam(token, &hc_heatcurr, PARAM_REAL, 1,1);
    }
#endif
#ifdef TTM
    else if (strcasecmp(token, "fd_g")==0){
      /* electron phonon coupling constant  */
      getparam("fd_g", &fd_g, PARAM_REAL, 1, 1);
    }
    else if (strcasecmp(token, "fd_update_steps")==0){
      /* how many steps before averaging over atoms to update FD cells  */
      getparam("fd_update_steps", &fd_update_steps, PARAM_INT, 1, 1);
    }
     else if (strcasecmp(token, "fd_ext")==0){
      /* how many MD cells in x,y,z-direction to one FD cell  */
      getparam("fd_ext", &fd_ext, PARAM_INT, DIM, DIM);
    }
    else if (strcasecmp(token, "fd_one_d")==0){
      /* FD lattice one dimensional in x or y or z if this is given  */
      getparam("fd_one_d", &fd_one_d_str, PARAM_STR, 1, 255);
    }
    else if (strcasecmp(token, "fd_k")==0){
      /* FD electronic heat conductivity  */
      getparam("fd_k", &fd_k, PARAM_REAL, 1, 1);
    }
    else if (strcasecmp(token, "fd_c")==0){
      /* FD electronic heat capacity  */
      getparam("fd_c", &fd_c, PARAM_REAL, 1, 1);
    }
    else if (strcasecmp(token, "fd_gamma")==0){
      /* FD electronic heat capacity / T_e (proport. const.) */
      getparam("fd_gamma", &fd_gamma, PARAM_REAL, 1, 1);
    }
    else if (strcasecmp(token, "fd_n_timesteps")==0){
      /* How many FD time steps to one MD time step?  */
      getparam("fd_n_timesteps", &fd_n_timesteps, PARAM_INT, 1, 1);
    }
    else if (strcasecmp(token, "ttm_int")==0){
      /* How many time steps between ttm writeouts?  */
      getparam("ttm_int", &ttm_int, PARAM_INT, 1, 1);
    }
    else if (strcasecmp(token, "init_t_el")==0){
      /* Initialize T_el to what temperature? */
      getparam("init_t_el", &init_t_el, PARAM_REAL, 1, 1);
    }
    else if (strcasecmp(token, "fix_t_el")==0){
      /* fix T_el to init_t_el? */
      getparam("fix_t_el", &fix_t_el, PARAM_INT, 1, 1);
    }
#endif
#ifdef PDECAY
else if (strcasecmp(token, "xipdecay")==0){
      /* value for the damping parameter added to the EQ's of motion */
      getparam("xipdecay", &xipdecay, PARAM_REAL, 1,1);
    }
else if (strcasecmp(token, "ramp_fraction")==0){
      /* fraction of the sample on which the damping ramp acts  */
      getparam("ramp_fraction", &ramp_fraction, PARAM_REAL, 1,1);
      if(ramp_fraction > 0.9 || ramp_fraction < 0.0)
	error("Nonsense value for ramp_fraction detected, please check prameter file!");
    }
else if (strcasecmp(token, "pdecay_mode")==0){
      /* mode for the damping function  */
      getparam("pdecay_mode", &pdecay_mode, PARAM_INT, 1,1);
    }
#endif
#ifdef LASER
else if (strcasecmp(token, "laser_delta_temp")==0){
      /* maximum heat added by laser (at the surface) (in maxwell routine) */
      getparam("laser_delta_temp", &laser_delta_temp, PARAM_REAL, 1,1);
    }
else if (strcasecmp(token, "laser_mu")==0){
      /* absorption coefficient (always needed)*/
      getparam("laser_mu", &laser_mu, PARAM_REAL, 1,1);
    }
else if (strcasecmp(token, "laser_offset")==0){
      /* offset of sample from origin */
      getparam("laser_offset", &laser_offset, PARAM_REAL, 1,1);
    }

else if (strcasecmp(token, "laser_dir")==0){
      /* direction of incidence of laser
       ( for now only along coordinate axes ) (always needed)*/
      getparam("laser_dir", &laser_dir, PARAM_INT, DIM, DIM);
    }
#ifdef LASERYZ
else if (strcasecmp(token, "laser_sigma_w_y")==0){
      /* y-center of gaussian laser-pulse  */
      getparam("laser_sigma_w_y", &laser_sigma_w_y, PARAM_REAL, 1,1);
    }
else if (strcasecmp(token, "laser_sigma_w_z")==0){
      /* z-center of gaussian laser-pulse  */
      getparam("laser_sigma_w_z", &laser_sigma_w_z, PARAM_REAL, 1,1);
    }
else if (strcasecmp(token, "laser_sigma_w0")==0){
      /* sigma_0 of spacial laser fluence  */
      getparam("laser_sigma_w0", &laser_sigma_w0, PARAM_REAL, 1,1);
    }
else if (strcasecmp(token,"laser_tem_mode")==0) {
      /* TEM_xy laser mode */
      getparam(token,&laser_tem_mode,PARAM_INT,3,3);
      switch ( laser_tem_mode.x ) /* Gauss Laguerre = 0, Gauss Hermite = 1 */
      {
	  case 0:
	  {
	    switch ( laser_tem_mode.y )
	    {
		case 0:
		{
		  switch ( laser_tem_mode.z )
		  {
		      case 0: laser_intensity_profile = laser_intensity_profile_laguerre_00; break;
		      case 1: laser_intensity_profile = laser_intensity_profile_laguerre_01; break;
		      case 2: laser_intensity_profile = laser_intensity_profile_laguerre_02; break;
		      case 3: laser_intensity_profile = laser_intensity_profile_laguerre_03; break;
		  }
		} break;

		case 1:
		{
		  switch ( laser_tem_mode.z )
		  {
		      case 0: laser_intensity_profile = laser_intensity_profile_laguerre_10; break;
		      case 1: laser_intensity_profile = laser_intensity_profile_laguerre_11; break;
		      case 2: laser_intensity_profile = laser_intensity_profile_laguerre_12; break;
		      case 3: laser_intensity_profile = laser_intensity_profile_laguerre_13; break;
		  }
		} break;

		case 2:
		{
		  switch ( laser_tem_mode.z )
		  {
		      case 0: laser_intensity_profile = laser_intensity_profile_laguerre_20; break;
		      case 1: laser_intensity_profile = laser_intensity_profile_laguerre_21; break;
		      case 2: laser_intensity_profile = laser_intensity_profile_laguerre_22; break;
		      case 3: laser_intensity_profile = laser_intensity_profile_laguerre_23; break;
		  }
		} break;

		case 3:
		{
		  switch ( laser_tem_mode.z )
		  {
		      case 0: laser_intensity_profile = laser_intensity_profile_laguerre_30; break;
		      case 1: laser_intensity_profile = laser_intensity_profile_laguerre_31; break;
		      case 2: laser_intensity_profile = laser_intensity_profile_laguerre_32; break;
		      case 3: laser_intensity_profile = laser_intensity_profile_laguerre_33; break;
		  }
		} break;

	    }
	  } break;
	  case 1:
	  {
	    switch ( laser_tem_mode.y )
	    {
		case 0:
		{
		  switch ( laser_tem_mode.z )
		  {
		      case 0: laser_intensity_profile = laser_intensity_profile_hermite_00; break;
		      case 1: laser_intensity_profile = laser_intensity_profile_hermite_01; break;
		      case 2: laser_intensity_profile = laser_intensity_profile_hermite_02; break;
		      case 3: laser_intensity_profile = laser_intensity_profile_hermite_03; break;
		  }
		} break;

		case 1:
		{
		  switch ( laser_tem_mode.z )
		  {
		      case 0: laser_intensity_profile = laser_intensity_profile_hermite_10; break;
		      case 1: laser_intensity_profile = laser_intensity_profile_hermite_11; break;
		      case 2: laser_intensity_profile = laser_intensity_profile_hermite_12; break;
		      case 3: laser_intensity_profile = laser_intensity_profile_hermite_13; break;
		  }
		} break;

		case 2:
		{
		  switch ( laser_tem_mode.z )
		  {
		      case 0: laser_intensity_profile = laser_intensity_profile_hermite_20; break;
		      case 1: laser_intensity_profile = laser_intensity_profile_hermite_21; break;
		      case 2: laser_intensity_profile = laser_intensity_profile_hermite_22; break;
		      case 3: laser_intensity_profile = laser_intensity_profile_hermite_23; break;
		  }
		} break;

		case 3:
		{
		  switch ( laser_tem_mode.z )
		  {
		      case 0: laser_intensity_profile = laser_intensity_profile_hermite_30; break;
		      case 1: laser_intensity_profile = laser_intensity_profile_hermite_31; break;
		      case 2: laser_intensity_profile = laser_intensity_profile_hermite_32; break;
		      case 3: laser_intensity_profile = laser_intensity_profile_hermite_33; break;
		  }
		} break;

	    }
	  } break;
      }
}

#endif
else if (strcasecmp(token, "laser_sigma_e")==0){
      /* area density of pulse energy (for rescaling method) */
      getparam("laser_sigma_e", &laser_sigma_e, PARAM_REAL, 1,1);
    }
else if (strcasecmp(token, "laser_sigma_t")==0){
      /* Pulse duration ( power is 1/e*P_max at t=t_0 +/- sigma_t ) */
      getparam("laser_sigma_t", &laser_sigma_t, PARAM_REAL, 1,1);
    }
else if (strcasecmp(token, "laser_t_0")==0){
      /* time of maximum pulse intensity */
      getparam("laser_t_0", &laser_t_0, PARAM_REAL, 1,1);
    }
else if (strcasecmp(token, "laser_sigma_e1")==0){
      /* area density of pulse energy (for rescaling method) */
      getparam("laser_sigma_e1", &laser_sigma_e1, PARAM_REAL, 1,1);
    }
else if (strcasecmp(token, "laser_sigma_t1")==0){
      /* Pulse duration ( power is 1/e*P_max at t=t_0 +/- sigma_t ) */
      getparam("laser_sigma_t1", &laser_sigma_t1, PARAM_REAL, 1,1);
    }
else if (strcasecmp(token, "laser_t_1")==0){
      /* time of maximum of second pulse intensity */
      getparam("laser_t_1", &laser_t_1, PARAM_REAL, 1,1);
    }
else if (strcasecmp(token, "laser_atom_vol")==0){
      /* Volume per particle (inverse density) */
      getparam("laser_atom_vol", &laser_atom_vol, PARAM_REAL, 1,1);
    }
else if (strcasecmp(token,"laser_rescale_mode")==0) {
      /* What rescale mode? */
      getparam("laser_rescale_mode",&laser_rescale_mode,PARAM_INT,1,1);
      switch ( laser_rescale_mode ) {
        case 0 :
          do_laser_rescale = laser_rescale_dummy;
          break;
        case 1 :
          do_laser_rescale = laser_rescale_1;
          break;
        case 2 :
          do_laser_rescale = laser_rescale_2;
          break;
        case 3 :
          do_laser_rescale = laser_rescale_3;
          break;
        case 4 :
#ifdef TTM
          do_laser_rescale = laser_rescale_ttm;
	  /* change electron temperature source terms, not atom velocities */
#else
	  error("Please compile with TTM if you want to use this laser rescale mode.\n");
#endif
          break;
        default :
          error("Illegal value for parameter laser_rescale_mode.\n");
          break;
      } /* switch */
    }
#endif

#ifdef STRESS_TENS
     else if (strcasecmp(token, "press_int")==0){
      /* number of steps between pressure writes  */
      getparam(token, &press_int, PARAM_INT, 1,1);
    }
     else if (strcasecmp(token, "presstens_ext")==0){
      /* external pressure tensor for relaxation */
      getparam(token, &presstens_ext, PARAM_REAL, DIM*(DIM+1)/2,DIM*(DIM+1)/2);
    }
#endif
#ifdef CNA
    else if (strcasecmp(token,"cna_start")==0) {
      /* step at which CNA begins */
      getparam("cna_start",&cna_start,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"cna_end")==0) {
      /* step at which CNA ends */
      getparam("cna_end",&cna_end,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"cna_int")==0) {
      /* number of steps between CNA */
      getparam("cna_int",&cna_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token, "cna_rcut")==0){
      /* cutoff */
      getparam("cna_rcut", &cna_rcut, PARAM_REAL, 1,1);
    }
    else if (strcasecmp(token,"cna_write")==0) {
      /* pair type to be written out */
      cna_write_n = getparam("cna_write",cna_writev,PARAM_INT,1,8);
    }
    else if (strcasecmp(token,"cna_crist")==0) {
      /* determine crystallinity of atoms */
      cna_crist_n = getparam("cna_crist",cna_cristv,PARAM_INT,1,4);
    }
    else if (strcasecmp(token,"cna_ll")==0) {
      /* lower left corner of partial box */
      getparam("cna_ll", &cna_ll,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"cna_ur")==0) {
      /* upper right corner of partial box */
      getparam("cna_ur", &cna_ur,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"cna_stat")==0) {
      /* write statistics */
      cna_write_statistics = 1;
    }
#endif
#ifdef ADA
   else if (strcasecmp(token, "ada_nbr_rcut") == 0) {
		/* cutoff radius for angle distribution analysis neighbor,*/
		/* squared value is needed*/
		getparam("ada_nbr_rcut", &ada_nbr_r2cut, PARAM_REAL, 1, 1);
		ada_nbr_r2cut = ada_nbr_r2cut * ada_nbr_r2cut;
	}
    else if (strcasecmp(token, "ada_write_int") == 0) {
   		/* write statistics interval*/
       	getparam("ada_write_int",&ada_write_int,PARAM_INT,1,1);
   	}
    else if (strcasecmp(token, "ada_crystal_structure") == 0) {
    	getparam(token,tmpstr,PARAM_STR,1,255);
    	if (strcasecmp(tmpstr,"fcc") == 0) {
    		ada_crystal_structure = ADA_FCC_CONFIG;
    	} else if (strcasecmp(tmpstr,"bcc") == 0) {
    		ada_crystal_structure = ADA_BCC_CONFIG;
    	} else if (strcasecmp(tmpstr,"ackland") == 0) {
    		ada_crystal_structure = ADA_ACKLAND_CONFIG;
    	}
    }
    else if (strcasecmp(token, "ada_latticeConst") == 0) {
      	/* lattice constant */
       	getparam("ada_latticeConst",&ada_latticeConst,PARAM_REAL,1,1);
    }
#endif
#ifdef NYETENSOR
    else if (strcasecmp(token,"nye_rotationAxis_x")==0) {
	  /* lattice orientation in x direction */
	  getparam("nye_rotationAxis_x",&nye_rotationAxis_x,PARAM_REAL,3,3);
	}
	else if (strcasecmp(token,"nye_rotationAxis_y")==0) {
	  /* lattice orientation in y direction */
	  getparam("nye_rotationAxis_y",&nye_rotationAxis_y,PARAM_REAL,3,3);
	}
	else if (strcasecmp(token,"nye_rotationAxis_z")==0) {
	  /* lattice orientation in z direction */
	  getparam("nye_rotationAxis_z",&nye_rotationAxis_z,PARAM_REAL,3,3);
	}
#endif

#ifdef DISLOC
     else if (strcasecmp(token,"reffile")==0) {
       /* filename for reference configuration */
       error(
       "Parameter reffile no longer supported - consult DISLOC documentation");
     }
    else if (strcasecmp(token,"dem_int")==0) {
      /* number of steps between picture writes */
      getparam(token,&dem_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"dsp_int")==0) {
      /* number of steps between picture writes */
      getparam(token,&dsp_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"update_ort_ref")==0) {
      /* step number to compute ort_ref */
      getparam(token,&up_ort_ref,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"min_dpot")==0) {
      /* minimum Epot difference */
      getparam(token,&min_dpot,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"min_dsp2")==0) {
      /* minimum square displacement in .dsp files */
      getparam(token,&min_dsp2,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"reset_Epot_step")==0) {
      /* step at which to compute Epot_ref (if calc_Epot_ref==1) */
      getparam(token,&reset_Epot_step,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"calc_Epot_ref")==0) {
      /* read (0) or compute (1) reference potential energy */
      getparam(token,&calc_Epot_ref,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"Epot_diff")==0) {
      /* write Epot (0) or Epot_diff (1) */
      getparam(token,&Epot_diff,PARAM_INT,1,1);
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
    else if (strcasecmp(token,"avpos_steps")==0) {
      /* number of steps to average over before position writes */
      getparam("avpos_steps",&avpos_steps,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"avpos_nwrites")==0) {
      /* number of position writes, only for processing the itr file */
      getparam("avpos_nwrites",&avpos_nwrites,PARAM_INT,1,1);
      printf("avpos_nwrites: %d\n",avpos_nwrites);fflush(stdout); 
    }
    else if (strcasecmp(token,"avpos_npwrites")==0) {
      /* number of pressure writes, only for processing the itr file */
      getparam("avpos_npwrites",&avpos_npwrites,PARAM_INT,1,1);
      printf("avpos_npwrites: %d\n",avpos_npwrites);fflush(stdout); 
    }
#endif
#if defined(FORCE) || defined(WRITEF)
    else if (strcasecmp(token,"force_int")==0) {
      /* number of steps between force writes */
      getparam(token,&force_int,PARAM_INT,1,1);
    }
#endif
#ifdef WRITEF
    else if (strcasecmp(token,"force_all")==0) {
      /* write all forces, or only those of atoms with virtual type */
      getparam(token,&force_all,PARAM_INT,1,1);
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
      if (ntypes==0) error("specify parameter ntypes before op_rcut");
      getparam(token,op_r2_cut,PARAM_REAL,ntypes*ntypes,ntypes*ntypes);
      for (k=0; k<ntypes*ntypes; k++) op_r2_cut[k] = SQR(op_r2_cut[k]);
    }
    else if (strcasecmp(token,"op_weight")==0) {
      /* weights for order parameter */
      if (ntypes==0) error("specify parameter ntypes before op_weight");
      getparam(token,op_weight,PARAM_REAL,ntypes*ntypes,ntypes*ntypes);
    }
#endif
#ifdef SOCKET_IO
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
      /* EAM2:Filename for the tabulated core-core potential (r^2) */
      getparam("core_potential_file",potfilename,PARAM_STR,1,255);
      have_potfile = 1;
    }
    else if (strcasecmp(token,"embedding_energy_file")==0) {
      /* EAM2:Filename for the tabulated embedding energy(rho_h) */
      getparam("embedding_energy_file",eam2_emb_E_filename,PARAM_STR,1,255);
    }
   else if (strcasecmp(token,"atomic_e-density_file")==0) {
      /* EAM2:Filename for the tabulated atomic electron density(r_ij^2) */
      getparam("atomic_e-density_file",eam2_at_rho_filename,PARAM_STR,1,255);
    }
#ifdef EEAM
    else if (strcasecmp(token,"eeam_energy_file")==0) {
      /* EEAM:Filename for the tabulated energy modification term(p_h) */
      getparam(token,eeam_mod_E_filename,PARAM_STR,1,255);
    }
#endif
#endif
#ifdef ADP
    else if (strcasecmp(token,"adp_upotfile")==0) {
      /* ADP dipole distortion potential */
      getparam(token,adp_upotfile,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"adp_wpotfile")==0) {
      /* ADP quadrupole distortion potential */
      getparam(token,adp_wpotfile,PARAM_STR,1,255);
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
    /* Gauss Part of Lennard-Jones-Gauss */
    else if (strcasecmp(token,"ljg_eps")==0) {
      if (ntypes==0) error("specify parameter ntypes before lj_sigma");
      getparam(token, ljg_eps_lin, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"ljg_r0")==0) {
      if (ntypes==0) error("specify parameter ntypes before lj_sigma");
      getparam(token, ljg_r0_lin, PARAM_REAL, ntypepairs, ntypepairs);
    }
    else if (strcasecmp(token,"ljg_sig")==0) {
      if (ntypes==0) error("specify parameter ntypes before lj_sigma");
      getparam(token, ljg_sig_lin, PARAM_REAL, ntypepairs, ntypepairs);
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
    else if (strcasecmp(token,"fix_bks")==0) {
      /* fix bks potential */
      getparam(token, &fix_bks, PARAM_INT, 1, 1);
    }
#endif
#ifdef SM
/* no charge update after each step */
    else if (strcasecmp(token,"charge_update_steps")==0) {
      getparam(token, &charge_update_steps, PARAM_INT, 1, 1);
    }
    /* keep charges fixed? */
    else if (strcasecmp(token,"sm_fixed_charges")==0) {
      getparam(token, &sm_fixed_charges, PARAM_INT, 1, 1);
    }
    /* Initial value of the electronegativity */
    else if (strcasecmp(token,"sm_chi_0")==0) {
      if (ntypes==0) error("specify parameter ntypes before sm_chi_0");
      getparam(token, sm_chi_0, PARAM_REAL, ntypes, ntypes);
    }
    /* Initial value of the effecitve core charge */
    else if (strcasecmp(token,"sm_Z")==0) {
      if (ntypes==0) error("specify parameter ntypes before sm_Z");
      getparam(token, sm_Z, PARAM_REAL, ntypes, ntypes);
    }
    /* SM zeta */
    else if (strcasecmp(token,"sm_zeta")==0) {
      if (ntypes==0) error("specify parameter ntypes before sm_zeta");
      getparam(token, sm_zeta, PARAM_REAL, ntypes, ntypes);
    }
    /* atomic hardness or self-Coulomb repulsion */
    else if (strcasecmp(token,"sm_J_0")==0) {
      if (ntypes==0) error("specify parameter ntypes before sm_J_0");
      getparam(token, sm_J_0, PARAM_REAL, ntypes, ntypes);
    }
    /* nuclear attraction potential */
    else if (strcasecmp(token,"na_pot_file")==0) {
      getparam(token, na_pot_filename,PARAM_STR,1,255);
    }
    /* coulomb repulsive potential */
    else if (strcasecmp(token,"cr_pot_file")==0) {
      getparam(token, cr_pot_filename,PARAM_STR,1,255);
    }
#ifndef NBLIST
    /* tabulated function erfc/r */
    else if (strcasecmp(token,"erfc_file")==0) {
      getparam(token, erfc_filename,PARAM_STR,1,255);
    }
#endif
#endif /* SM */
#ifdef FEFL
    /*harmonic potential for Einstein crystal */
    else if (strcasecmp(token,"spring_rate")==0) {
      if (ntypes==0) error("specify parameter ntypes before spring_rate");
      getparam(token, spring_rate, PARAM_REAL, ntypes, ntypes);
      have_pre_pot = 1;
    }
    else if (strcasecmp(token,"lambda")==0) {
      /* parameter for potential switching */
      getparam(token, &lambda, PARAM_REAL,1,1);
    }
#endif
#if defined(COVALENT) || defined(NNBR_TABLE)
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
    else if (strcasecmp(token,"ttbp_constant2")==0) {
      /* force constant for vashishta potential */
      getparam(token, ttbp_constant2, PARAM_REAL, 8, 8);
      ttbp_vas = 1;
      /* ntypes^3 values, vashishta-option only supported for ntypes = 2 */
    }
    else if (strcasecmp(token,"ttbp_sp")==0) {
      /* hybridization of the element type */
      if (ntypes==0) error("specify parameter ntypes before ttbp_sp");
      getparam(token, ttbp_sp, PARAM_REAL, ntypes, ntypes);
    }
    else if (strcasecmp(token,"ttbp_cut")==0) {
      /* cutoff for smoothing part of vashishta potential */
      getparam(token, ttbp_cut, PARAM_REAL, 1, 1);
    }
    else if (strcasecmp(token,"ttbp_potfile")==0) {
      /* filename for ttbp potential data */
      getparam(token, ttbp_potfilename, PARAM_STR, 1, 255);
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
      have_pre_pot = 1;
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
#if defined(TERSOFF) || defined(TERSOFFMOD) || defined(BRENNER)
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
#if defined(TERSOFF) || defined(BRENNER)
      have_pre_pot = 1;
#endif
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
#if defined(TERSOFF) || defined(BRENNER)
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
    /* nvalues is ntypes for TERSOFF or TERSOFFMOD
    and ntypepairs for TERSOFF2 or TERSOFFMOD2 */
    else if (strcasecmp(token,"ters_ga")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_ga");
      getparam(token, ters_ga, PARAM_REAL, nvalues, nvalues);
    }
    else if (strcasecmp(token,"ters_n")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_n");
      getparam(token, ters_n, PARAM_REAL, nvalues, nvalues);
    }
    else if (strcasecmp(token,"ters_c")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_c");
      getparam(token, ters_c, PARAM_REAL, nvalues, nvalues);
    }
    else if (strcasecmp(token,"ters_d")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_d");
      getparam(token, ters_d, PARAM_REAL, nvalues, nvalues);
    }
#else /* TERSOFFMOD */
    else if (strcasecmp(token,"ters_eta")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_eta");
      getparam(token, ters_eta, PARAM_REAL, nvalues, nvalues);
    }
    else if (strcasecmp(token,"ters_delta")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_delta");
      getparam(token, ters_delta, PARAM_REAL, nvalues, nvalues);
    }
    else if (strcasecmp(token,"ters_alpha")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_alpha");
      getparam(token, ters_alpha, PARAM_REAL, nvalues, nvalues);
    }
    else if (strcasecmp(token,"ters_beta")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_beta");
      getparam(token, ters_beta, PARAM_REAL, nvalues, nvalues);
    }
    else if (strcasecmp(token,"ters_c1")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_c1");
      getparam(token, ters_c1, PARAM_REAL, nvalues, nvalues);
    }
    else if (strcasecmp(token,"ters_c2")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_c2");
      getparam(token, ters_c2, PARAM_REAL, nvalues, nvalues);
    }
    else if (strcasecmp(token,"ters_c3")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_c3");
      getparam(token, ters_c3, PARAM_REAL, nvalues, nvalues);
    }
    else if (strcasecmp(token,"ters_c4")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_c4");
      getparam(token, ters_c4, PARAM_REAL, nvalues, nvalues);
    }
    else if (strcasecmp(token,"ters_c5")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_c5");
      getparam(token, ters_c5, PARAM_REAL, nvalues, nvalues);
    }
#endif
    else if (strcasecmp(token,"ters_h")==0) {
      if (ntypes==0) error("specify parameter ntypes before ters_h");
      getparam(token, ters_h, PARAM_REAL, nvalues, nvalues);
    }
#endif
#ifdef BRENNER
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
#if defined(EWALD) || defined(COULOMB) || defined(USEFCS)
    /* charge */
    else if (strcasecmp(token,"charge")==0) {
      if (ntypes==0) error("specify parameter ntypes before charge");
      getparam(token,charge,PARAM_REAL,ntypes,ntypes);
    }
    /* coul_eng */
    else if (strcasecmp(token,"coul_eng")==0) {
      getparam(token,&coul_eng,PARAM_REAL,1,1);
    }
#endif
#ifdef USEFCS
    else if (strcasecmp(token,"fcs_method")==0) {
      /* FCS method */
      getparam(token,tmpstr,PARAM_STR,1,255);
      if (strcasecmp(tmpstr,"direct")==0) {
        fcs_method = FCS_METH_DIRECT;
      }
      else if (strcasecmp(tmpstr,"pepc")==0) {
        fcs_method = FCS_METH_PEPC;
      }
      else if (strcasecmp(tmpstr,"fmm")==0) {
        fcs_method = FCS_METH_FMM;
      }
      else if (strcasecmp(tmpstr,"p3m")==0) {
        fcs_method = FCS_METH_P3M;
      }
      else if (strcasecmp(tmpstr,"p2nfft")==0) {
        fcs_method = FCS_METH_P2NFFT;
      }
      else if (strcasecmp(tmpstr,"vmg")==0) {
        fcs_method = FCS_METH_VMG;
      }
      else if (strcasecmp(tmpstr,"pp3mg")==0) {
        fcs_method = FCS_METH_PP3MG;
      }
    }
#ifdef PAIR
    /* delegate near-field to IMD? */
    else if (strcasecmp(token,"fcs_near_field_flag")==0) {
      getparam(token,&fcs_near_field_flag,PARAM_INT,1,1);
    }
    /* fcs_rcut for near-field delegation */
    else if (strcasecmp(token,"fcs_rcut")==0) {
      getparam(token,&fcs_rcut,PARAM_REAL,1,1);
      if (fcs_rcut > 0) have_pre_pot = 1;
    }
#endif
    /* fcs_tolerance */
    else if (strcasecmp(token,"fcs_tolerance")==0) {
      getparam(token,&fcs_tolerance,PARAM_REAL,1,1);
    }
    /* fcs_grid_dim */
    else if (strcasecmp(token,"fcs_grid_dim")==0) {
      getparam(token,&fcs_grid_dim,PARAM_INT,3,3);
    }
    /* fcs_max_iter */
    else if (strcasecmp(token,"fcs_max_iter")==0) {
      getparam(token,&fcs_max_iter,PARAM_INT,1,1);
    }
    /* fcs_iter_tolerance */
    else if (strcasecmp(token,"fcs_iter_tolerance")==0) {
      getparam(token,&fcs_iter_tolerance,PARAM_REAL,1,1);
    }
    /* fcs_pepc_eps */
    else if (strcasecmp(token,"fcs_pepc_eps")==0) {
      getparam(token,&fcs_pepc_eps,PARAM_REAL,1,1);
    }
    /* fcs_pepc_theta */
    else if (strcasecmp(token,"fcs_pepc_theta")==0) {
      getparam(token,&fcs_pepc_theta,PARAM_REAL,1,1);
    }
    /* fcs_pepc_nthreads */
    else if (strcasecmp(token,"fcs_pepc_nthreads")==0) {
      getparam(token,&fcs_pepc_nthreads,PARAM_INT,1,1);
    }
    /* fcs_fmm_absrel */
    else if (strcasecmp(token,"fcs_fmm_absrel")==0) {
      getparam(token,&fcs_fmm_absrel,PARAM_INT,1,1);
    }
    /* fcs_fmm_dcorr */
    else if (strcasecmp(token,"fcs_fmm_dcorr")==0) {
      getparam(token,&fcs_fmm_dcorr,PARAM_INT,1,1);
    }
    /* fcs_fmm_do_tune */
    else if (strcasecmp(token,"fcs_fmm_do_tune")==0) {
      getparam(token,&fcs_fmm_do_tune,PARAM_INT,1,1);
    }
    /* fcs_vmg_max_level */
    else if (strcasecmp(token,"fcs_vmg_max_level")==0) {
      getparam(token,&fcs_vmg_max_level,PARAM_INT,1,1);
    }
    /* fcs_vmg_smooth_steps */
    else if (strcasecmp(token,"fcs_vmg_smooth_steps")==0) {
      getparam(token,&fcs_vmg_smooth_steps,PARAM_INT,1,1);
    }
    /* fcs_vmg_gamma */
    else if (strcasecmp(token,"fcs_vmg_gamma")==0) {
      getparam(token,&fcs_vmg_gamma,PARAM_INT,1,1);
    }
    /* fcs_vmg_near_field_cells */
    else if (strcasecmp(token,"fcs_vmg_near_field_cells")==0) {
      getparam(token,&fcs_vmg_near_field_cells,PARAM_INT,1,1);
    }
    /* fcs_vmg_interpol_order */
    else if (strcasecmp(token,"fcs_vmg_interpol_order")==0) {
      getparam(token,&fcs_vmg_interpol_order,PARAM_INT,1,1);
    }
    /* fcs_vmg_discr_order */
    else if (strcasecmp(token,"fcs_vmg_discr_order")==0) {
      getparam(token,&fcs_vmg_discr_order,PARAM_INT,1,1);
    }
    /* fcs_pp3mg_ghosts */
    else if (strcasecmp(token,"fcs_pp3mg_ghosts")==0) {
      getparam(token,&fcs_pp3mg_ghosts,PARAM_INT,1,1);
    }
    /* fcs_pp3mg_degree */
    else if (strcasecmp(token,"fcs_pp3mg_degree")==0) {
      getparam(token,&fcs_pp3mg_degree,PARAM_INT,1,1);
    }
    /* fcs_pp3mg_max_part */
    else if (strcasecmp(token,"fcs_pp3mg_max_part")==0) {
      getparam(token,&fcs_pp3mg_max_part,PARAM_INT,1,1);
    }
    /* fcs_p2nfft_intpol_order */
    else if (strcasecmp(token,"fcs_p2nfft_intpol_order")==0) {
      getparam(token,&fcs_p2nfft_intpol_order,PARAM_INT,1,1);
    }
    /* fcs_p2nfft_epsI */
    else if (strcasecmp(token,"fcs_p2nfft_epsI")==0) {
      getparam(token,&fcs_p2nfft_epsI,PARAM_REAL,1,1);
    }
#endif /* USEFCS */
#if defined(EWALD) || defined(COULOMB)
    /* smoothing parameter */
    else if (strcasecmp(token,"ew_kappa")==0) {
      getparam(token,&ew_kappa,PARAM_REAL,1,1);
    }
    /* k-space cutoff */
    else if (strcasecmp(token,"ew_kcut")==0) {
      getparam(token,&ew_kcut,PARAM_REAL,1,1);
    }
    /* r-space cutoff */
    else if (strcasecmp(token,"ew_rcut")==0) {
      getparam(token,&rtmp,PARAM_REAL,1,1);
      ew_r2_cut = SQR(rtmp);
#ifdef DIPOLE
      dp_self=2./(3.*rtmp*ew_r2_cut*sqrt(2.*M_PI));
#endif /* DIPOLE */
#ifndef VARCHG
      have_pre_pot = 1;
#endif
    }
    /* number of image boxes */
    else if (strcasecmp(token,"ew_nmax")==0) {
      getparam(token,&ew_nmax,PARAM_INT,1,1);
    }
    /* test flag */
    else if (strcasecmp(token,"ew_test")==0) {
      getparam(token,&ew_test,PARAM_INT,1,1);
    }
    /* potential table resolution */
    else if (strcasecmp(token,"coul_res")==0) {
      getparam(token,&coul_res,PARAM_REAL,1,1);
    }
    /* potential table resolution - backwards compatibility */
    else if (strcasecmp(token,"dp_res")==0) {
      getparam(token,&coul_res,PARAM_REAL,1,1);
    }
    /* potential table resolution */
    else if (strcasecmp(token,"coul_begin")==0) {
      getparam(token,&coul_begin,PARAM_REAL,1,1);
    }
    /* potential table resolution - backwards compatibility */
    else if (strcasecmp(token,"dp_begin")==0) {
      getparam(token,&coul_begin,PARAM_REAL,1,1);
    }
#endif /* EWALD or COULOMB */
#ifdef DIPOLE
    /* dipole fixed? */
    else if (strcasecmp(token,"dp_fix")==0) {
      getparam(token,&dp_fix,PARAM_INT,1,1);
    }
    /* dipole field mixing param */
    else if (strcasecmp(token,"dp_mix")==0) {
      getparam(token,&dp_mix,PARAM_REAL,1,1);
    }
    /* dipole iteration precision */
    else if (strcasecmp(token,"dp_tol")==0) {
      getparam(token,&dp_tol,PARAM_REAL,1,1);
    }
    /* polarisability */
    else if (strcasecmp(token,"dp_alpha")==0) {
      if (ntypes==0) error("specify parameter ntypes before dp_alpha");
      dp_alpha = (real *) malloc( ntypes*sizeof(real));
      if (NULL==dp_alpha) error("cannot allocate dp_alpha");
      getparam(token,dp_alpha,PARAM_REAL,ntypes,ntypes);
    }
    /* short-range dipole parameter b */
    else if (strcasecmp(token,"dp_b")==0) {
      if (ntypepairs==0) error("specify parameter ntypes before dp_b");
      dp_b = (real *) malloc( ntypepairs*sizeof(real));
      if (NULL==dp_b) error("cannot allocate dp_b");
      getparam(token,dp_b,PARAM_REAL,ntypepairs,ntypepairs);
    }
    /* short-range dipole parameter c */
    else if (strcasecmp(token,"dp_c")==0) {
      if (ntypepairs==0) error("specify parameter ntypes before dp_c");
      dp_c = (real *) malloc( ntypepairs*sizeof(real));
      if (NULL==dp_c) error("cannot allocate dp_c");
      getparam(token,dp_c,PARAM_REAL,ntypepairs,ntypepairs);
    }
#endif /* DIPOLE */
#if ((defined(DIPOLE) || defined(MORSE)) && !defined(BUCK))
    /* Morse-Stretch parameter D */
    else if (strcasecmp(token,"ms_D")==0) {
      if (ntypes==0) error("specify parameter ntypes before ms_D");
      ms_D = (real *) malloc( ntypepairs*sizeof(real));
      if (NULL==ms_D) error("cannot allocate ms_D");
      getparam(token,ms_D,PARAM_REAL,ntypepairs,ntypepairs);
      have_pre_pot = 1;
    }
    /* Morse-Stretch parameter gamma */
    else if (strcasecmp(token,"ms_gamma")==0) {
      if (ntypes==0) error("specify parameter ntypes before ms_gamma");
      ms_gamma = (real *) malloc( ntypepairs*sizeof(real));
      if (NULL==ms_gamma) error("cannot allocate ms_gamma");
      getparam(token,ms_gamma,PARAM_REAL,ntypepairs,ntypepairs);
    }
    /* Morse-Stretch parameter r0 */
    else if (strcasecmp(token,"ms_r0")==0) {
      if (ntypes==0) error("specify parameter ntypes before ms_r0");
      ms_r0 = (real *) malloc( ntypepairs*sizeof(real));
      if (NULL==ms_r0) error("cannot allocate ms_r0");
      getparam(token,ms_r0,PARAM_REAL,ntypepairs,ntypepairs);
    }
    /* Morse-Stretch spring constant ms_harm_c */
    else if (strcasecmp(token,"ms_harm_c")==0) {
      if (ntypes==0) error("specify parameter ntypes before ms_harm_c");
      ms_harm_c = (real *) malloc( ntypepairs*sizeof(real));
      if (NULL==ms_harm_c) error("cannot allocate ms_harm_c");
      getparam(token,ms_harm_c,PARAM_REAL,ntypepairs,ntypepairs);
    }
    /* Morse-Stretch minimum distance ms_rmin */
    else if (strcasecmp(token,"ms_rmin")==0) {
      if (ntypes==0) error("specify parameter ntypes before ms_rmin");
      ms_r2_min = (real *) malloc( ntypepairs*sizeof(real));
      if (NULL==ms_r2_min) error("cannot allocate ms_r2_min");
      getparam(token,ms_r2_min,PARAM_REAL,ntypepairs,ntypepairs);
      for (i=0;i<ntypepairs;i++) {
	rtmp = SQR(ms_r2_min[i]);
	ms_r2_min[i] = rtmp;
      }
    }
#endif /* DIPOLE or MORSE */

#ifdef EXTF
    /* external homogeneous electrostatic field */
    else if (strcasecmp(token,"extf")==0) {
      getparam(token,&extf,PARAM_REAL,DIM,DIM);
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

#ifdef RELAX
    else if (strcasecmp(token,"max_sscount")==0) {
      /* max nr. of minimizations in quasistat sims */
        getparam(token,&max_sscount,PARAM_INT,1,1);
    }
#endif

#ifdef EXTPOT
    else if (strcasecmp(token,"ep_n")==0) {
      /* EXTPOT: number of external potentials */
      getparam(token,&ep_n,PARAM_INT,1,1);
      if (1 > ep_n) error("At least one extpot must be defined (ep_n < 1)");
      for (i=0; i<ep_n; i++) {
        ep_pos[i] = nullv;
        ep_vel[i] = nullv;
        ep_dir[i] = nullv;
      }
    }
    else if (strcasecmp(token,"ep_nind")==0) {
      /* EXTPOT: number of indentors (remaining extpots are walls) */
        getparam(token,&ep_nind,PARAM_INT,1,1);
        if (ep_nind > ep_n)
        error("Number of indeters exceeds ep_n");
    }
    else if (strcasecmp(token,"ep_key")==0) {
        /* EXTPOT:  potential key : which potential type to use*/
      getparam(token,&ep_key,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"ep_a")==0) {
      /* EXTPOT: strength of external potential */
      getparam(token,&ep_a,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"extpot_file")==0) {
      /* EXTPOT: Filename for the tabulated external potential (r^2) */
      getparam("expot_file",extpotfilename,PARAM_STR,1,255);
      if(ep_a!=0)
	 error("either use expotfile or define ep_a");
      have_extpotfile = 1;
    }
    else if (strcasecmp(token,"ep_max_int")==0) {
      /* EXTPOT: maximal wait steps during relaxation */
      getparam(token,&ep_max_int,PARAM_INT,1,1);
      if ((ep_max_int > 1) && (ep_max_int < 10)) ep_max_int = 10;
    }
    else if (strcasecmp(token,"ep_rcut")==0) {
      /* EXTPOT: cutoff radius of external potential */
      getparam(token,&ep_rcut,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"ep_pos")==0) {
      /* EXTPOT: position of external potential */
      getparam(token,&rtmp4,PARAM_REAL,4,4);
      i = (int) rtmp4[0];
      if (i >= ep_n) error("Number of external potential exceeds ep_n");
      ep_pos[i].x = rtmp4[1];
      ep_pos[i].y = rtmp4[2];
      ep_pos[i].z = rtmp4[3];
    }
    else if (strcasecmp(token,"ep_vel")==0) {
      /* EXTPOT: velocity of external potential */
      getparam(token,&rtmp4,PARAM_REAL,4,4);
      i = (int) rtmp4[0];
      if (i >= ep_n) error("Number of external potential exceeds ep_n");
      ep_vel[i].x = rtmp4[1];
      ep_vel[i].y = rtmp4[2];
      ep_vel[i].z = rtmp4[3];
    }
    else if (strcasecmp(token,"ep_dir")==0) {
      /* EXTPOT: direction of external potential */
      getparam(token,&rtmp4,PARAM_REAL,4,4);
      i = (int) rtmp4[0];
      if (i >= ep_n) error("Number of external potential exceeds ep_n");
      rtmp = SQRT( SQR(rtmp4[1]) + SQR(rtmp4[2]) + SQR(rtmp4[3]) );
      if (rtmp < 1e-6) error("parameter ep_dir requires non-zero direction");
      ep_dir[i].x = rtmp4[1] / rtmp;
      ep_dir[i].y = rtmp4[2] / rtmp;
      ep_dir[i].z = rtmp4[3] / rtmp;

    }
#endif /* EXTPOT */
#ifdef CBE
    else if (strcasecmp(token,"num_spus")==0) {
      /* number of SPUs to be used */
      getparam(token,&num_spus,PARAM_INT,1,1);
      /* Make sure parameter just read is in a valid range */
      if ( (num_spus<1) || (num_spus>N_SPU_THREADS_MAX) ) {
	 num_spus=N_SPU_THREADS_MAX;
      }
    }
    else if (strcasecmp(token, "num_bufs")==0) {
      /* Number of argument buffers per SPU */
      getparam(token,&num_bufs,PARAM_INT,1,1);
      if ( (num_bufs<1) || (num_bufs>N_ARGBUF) ) {
	 num_bufs=N_ARGBUF;
      }
    }
    else if (strcasecmp(token, "cbe_pot_steps")==0) {
      /* number of tabulation steps in potential table */
      getparam(token,&cbe_pot_steps,PARAM_INT,1,1);
    }
    else if (strcasecmp(token, "cbe_pot_max")==0) {
      /* maximum value in potential table */
      getparam(token,&cbe_pot_max,PARAM_REAL,1,1);
    }
#endif

#ifdef KIM
    else if (strcasecmp(token, "kim_model_name")==0) {
      /* name of the KIM model for force calculation */
      getparam(token, kim_model_name, PARAM_STR, 1, 255);
      have_potfile = 1;
    }
    else if (strcasecmp(token, "kim_el_names")==0) {
      /* element names, needed for the KIM matching process */
      if (ntypes == 0) error ("specify parameter ntypes before kim_elements\n");
      getparam(token, kim_el_names, PARAM_CHARPTR, ntypes, 4);
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

  return finished;

} /* getparamfile */


/*****************************************************************
*
*   Check input for nonsense values
*
******************************************************************/

void check_parameters_complete()
{
  real tmp;
  real norm_bend_axis;
  int  k;
#ifdef TWOD
  vektor einsv = {1.0,1.0};
#else
  vektor einsv = {1.0,1.0,1.0};
#endif
  vektor this_bend_axis;

  if (ensemble == 0) {
    error("missing or unknown ensemble parameter.");
  }

  if (timestep == (real)0) {
    error("timestep is missing or zero.");
  }

  if (ntypes == 0) {
    error("ntypes is missing or zero.");
  }

#ifdef BEND
  if(bend_nmoments >0)
  {
      if(bend_nmoments >6)
          error("currently only 6 bending moments are supported");

      for (k=0;k<bend_nmoments;k++)
      {
          this_bend_axis.x = (bend_axis + k)->x;
          this_bend_axis.y = (bend_axis + k)->y;
          this_bend_axis.z = (bend_axis + k)->z;
          if(SPROD(this_bend_axis,this_bend_axis)==0)
          error("definition of bending moment without axis");

#ifdef RELAX
          tmp=((fbc_bdforces + (bend_vtype_of_force[k]))->x)*
              ((fbc_bdforces + (bend_vtype_of_force[k]))->x) +
              ((fbc_bdforces + (bend_vtype_of_force[k]))->y)*
              ((fbc_bdforces + (bend_vtype_of_force[k]))->y) +
              ((fbc_bdforces + (bend_vtype_of_force[k]))->z)*
              ((fbc_bdforces + (bend_vtype_of_force[k]))->z);
#else
          tmp=((fbc_endbforces + (bend_vtype_of_force[k]))->x)*
              ((fbc_endbforces + (bend_vtype_of_force[k]))->x) +
              ((fbc_endbforces + (bend_vtype_of_force[k]))->y)*
              ((fbc_endbforces + (bend_vtype_of_force[k]))->y) +
              ((fbc_endbforces + (bend_vtype_of_force[k]))->z)*
              ((fbc_endbforces + (bend_vtype_of_force[k]))->z);
#endif
          if(tmp==0)
              error("definition of bending moment without force");
      }
  }
#endif

#ifdef EXTPOT
  if(ep_a !=0)
    printf("Usage of ep_a is depreciated, use extpot_file instead\n");
  if(ep_a !=0 && have_extpotfile==1)
    error("use either ep_a or extpotfile");
  if(ep_rcut <=0)
    error("need a value for ep_rcut");
#endif

#if defined(FBC) || defined(RIGID) || defined(DEFORM)
  if (vtypes == 0)
    error("FBC, RIGID, and DEFORM require parameter total_types to be set");
#endif
  if (vtypes == 0) {
    vtypes = ntypes;
    restrictions = (vektor *) malloc( vtypes * sizeof(vektor) );
    if (NULL==restrictions)
      error("Cannot allocate memory for restriction vectors\n");
    for (k=0; k<vtypes; k++)
      restrictions[k] = einsv;
  }
  if (vtypes < ntypes)
    error("total_types must not be smaller than ntypes");

#if defined(PAIR) && !defined(VARCHG) && !defined(TERSOFFMOD)
  if ((have_potfile==0) && (have_pre_pot==0))
    error("You must specify a pair interaction!");
#endif
#ifdef TEMPCONTROL
  if (temperature == 0) {
    error("starttemp is missing or zero.");
  }
  if (end_temp == 0) {
    end_temp = temperature;
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
  if (ensemble==ENS_NVX) {
    if (hc_int     == 0) error ("hc_int is zero.");
    if (hc_nlayers == 0) error ("hc_nlayers is zero.");
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

#ifdef LASER
  if (laser_dir.x!=0) {
    laser_dir.x=1;
    if (laser_dir.y!=0
#ifndef TWOD
      || laser_dir.z!=0
#endif
    ) error("Sorry: Laser incidence only along one coordinate axis.");
  }
  else if (laser_dir.y!=0) {
    laser_dir.y=1;
#ifndef TWOD
    if (laser_dir.z!=0){
      error("Sorry: Laser incidence only along one coordinate axis.");
    }
  }
  else if (laser_dir.z!=0) {
    laser_dir.z=1;
#endif
  }
  else {
    error("Parameter laser_dir (laser incidence direction) missing.");
  }
  if ( (laser_rescale_mode < 0) || (laser_rescale_mode > 4) ) {
    error("Parameter laser_rescale_mode must be a positive integer < 5 !");
  }
#endif /* LASER */
#ifdef LASERYZ
	if ( (laser_sigma_w_y == 0.0 ) || (laser_sigma_w_z == 0.0) ){
	  warning("laser_sigma_w_y and / or laser_sigma_w_z is set to 0.0 - this seems to be a nonsense value!\n");
	}
	if ( (laser_sigma_w_y < 0.0) || (laser_sigma_w_z < 0.0)){
	  error("laser_sigma_w_y and / or laser_sigma_w_z is smaller than zero - this is nonsense!\n");
	}
	if ( laser_sigma_w0 <= 0.0){
	  error("laser_sigma_w0 is equal or less than zero - which is nonsense or means zero energy density.");
	}

	if ( (laser_tem_mode.x < 0 ) || (laser_tem_mode.x > 1 ) )
	{
	  error("Laser TEM Mode has to be either Gauss-Laguerre (0) or Gauss-Hermite (1).");
	}

	if ( (laser_tem_mode.y < 0) || (laser_tem_mode.y > 3) || (laser_tem_mode.z < 0 ) || (laser_tem_mode.z > 3) )
	{
	  error("Only Laser TEM_xy modes for x=1...3 and y=1...3 are supported yet.");
	}
#endif
#ifdef TTM
  if (fd_update_steps <= 0) {
    warning("Ignoring illegal value of fd_update_steps, using 1\n");
    fd_update_steps=1;
  }
  if (init_t_el<0) {
    warning("Ignoring illegal value of init_t_el, using lattice temp\n");
    init_t_el=0.0;
  }
  if (fix_t_el!=0 && init_t_el==0.0)
    error("You need to specify init_t_el for enabled fix_t_el!\n");

  if (strcasecmp(fd_one_d_str,"x")==0 || strcasecmp(fd_one_d_str,"1")==0) {
    fd_one_d=1;
  }
  else if (strcasecmp(fd_one_d_str,"y")==0 || strcasecmp(fd_one_d_str,"2")==0){
    fd_one_d=2;
  }
  else if (strcasecmp(fd_one_d_str,"z")==0 || strcasecmp(fd_one_d_str,"3")==0){
    fd_one_d=3;
  }
  else if (strcasecmp(fd_one_d_str,"")!=0) {
    warning("Ignoring unknown value of fe_one_d\n");
  }
  if ((fd_gamma==0.0 && fd_c==0.0)||(fd_gamma!=0.0 && fd_c!=0.0)) {
    error ("You must specify either fd_gamma or fd_c for TTM simulations.");
  }
#endif /* TTM */
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
#ifdef SOCKET_IO
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
  fprintf(stdout, "avpos_start: %d imdrestart*checkpt_int: %d\n", avpos_start, imdrestart*checkpt_int);
  // if (avpos_start <= imdrestart*checkpt_int)
    //  avpos_start = imdrestart*checkpt_int+1; /* do not ask me why +1 ;-) */
    //  avpos_start = imdrestart*checkpt_int; 
  /* Default initialisation of end time */
  if (0==avpos_end) avpos_end = steps_max;
#ifdef STRESS_TENS
  if (press_int % avpos_int !=0 || avpos_res % eng_int !=0)
    error("For averaged pressure writes, press_int musst be an integer multiple of avpos_int and avpos_res must be an integer multiple of eng_int");
#endif
#endif

#ifdef ATDIST
  if (0==atdist_end) atdist_end = steps_max;
#endif

#ifdef CG
  if ((linmin_maxsteps==0) || (linmin_tol==0.0) )
    error("You have to set parameters for the linmin search");
#endif
#ifdef HOMDEF
  if (relax_rate > 0.0) {
#ifdef STRESS_TENS
    if (relax_mode == -1) relax_mode = RELAX_FULL;
#else
    if (relax_mode == -1) relax_mode = RELAX_ISO;
    if ((relax_mode == RELAX_FULL) || (relax_mode == RELAX_AXIAL))
      error("Pressure relaxation modes axial and full require option stress");
#endif
  }
#endif
#if defined(DIFFPAT) && defined(TWOD)
  error("Option DIFFPAT is not supported in 2D");
#endif

#ifdef KIM
  if (strcmp(kim_el_names[0],"\0")==0)
    error("kim_el_names is not properly set in parameter file");
#endif

#if defined(ADA) && defined(TWOD)
  error("Option ADA is not supported in 2D");
#endif

#ifdef ADA
  if (ada_nbr_r2cut == 0. && ada_latticeConst == 0.){
	  error("Nearest neighbor cutoff distance \"ada_nbr_rcut\" or lattice constant \"ada_latticeConst\" is missing or zero in the parameter file");
  }
#endif

#ifdef NYETENSOR
  if (ada_latticeConst == 0.){
  	  error("Lattice constant \"ada_latticeConst\" is missing or zero in the parameter file");
  }
  if (SPROD(nye_rotationAxis_x, nye_rotationAxis_x) == 0.)
	  error("Crystal orientation in x direction \"nye_rotationAxis_x\" is missing or zero in the parameter file");
  if (SPROD(nye_rotationAxis_y, nye_rotationAxis_y) == 0.)
  	  error("Crystal orientation in y direction \"nye_rotationAxis_y\" is missing or zero in the parameter file");
  if (SPROD(nye_rotationAxis_z, nye_rotationAxis_z) == 0.)
  	  error("Crystal orientation in x direction \"nye_rotationAxis_z\" is missing or zero in the parameter file");
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
*  read parameters
*
******************************************************************/

int read_parameters(char *paramfname, int phase)
{
  int i, finished = 0;
  str255 fname;
  FILE *testfile;

  if (0==myid) {

    /* write itr-file for the next phase */
    if (phase > 1) write_itr_file(-1, steps_max,"");
    finished = getparamfile(paramfname, phase);
    /* read initial itr-file (if there is any), but keep steps_min value */
    if ((phase == 1) && (0 < strlen(itrfilename))) {
      int tmp_steps = steps_min;
      getparamfile(itrfilename, 1);
      steps_min = tmp_steps;
    }
    /* read back itr-file for the next phase */
    if (phase > 1) {
#ifdef NEB
      sprintf(outfilename, "%s.%02d", neb_outfilename, myrank);
#endif
      sprintf( itrfilename,"%s-final.itr", outfilename );
      getparamfile(itrfilename, 1);
    }
    check_parameters_complete();

    /* Get restart parameters if restart */
    if (0 != imdrestart) {

      /* read itr-file */
      sprintf(fname,"%s.%d.itr",outfilename,imdrestart);
      testfile = fopen(fname,"r");
      if (NULL==testfile) {
        sprintf(fname,"%s.%05d.itr",outfilename,imdrestart);
        testfile = fopen(fname,"r");
        if (NULL==testfile) {
          error_str("file %s not found", fname);
        } else {
          fclose(testfile);
        }
      } else {
        fclose(testfile);
      }
      getparamfile(fname,1);

      /* get restart configuration */
      sprintf(infilename,"%s.%d.%s",outfilename,imdrestart,"chkpt");
      testfile = fopen(infilename,"r");
      if (NULL==testfile) {
        sprintf(infilename,"%s.%05d.%s",outfilename,imdrestart,"chkpt");
        testfile = fopen(infilename,"r");
        if (NULL==testfile) {
          error_str("file %s not found", infilename);
        } else {
          fclose(testfile);
        }
      } else {
        fclose(testfile);
      }
      printf("Restarting from %s.\n",infilename);

    } else if (phase == 1) {
      /* if not restart and phase 1: delete files to which we append */
      sprintf(fname,"%s.eng",                 outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.Ekin",         outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.Epot",         outfilename); unlink(fname);
#ifdef STRESS_TENS
      sprintf(fname,"%s.minmax.press",        outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.presstens",    outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.presstens_xx", outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.presstens_yy", outfilename); unlink(fname);
#ifndef TWOD
      sprintf(fname,"%s.minmax.presstens_zz", outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.presstens_yz", outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.presstens_zx", outfilename); unlink(fname);
#endif
      sprintf(fname,"%s.minmax.presstens_xy", outfilename); unlink(fname);
#endif /* STRESS_TENS */
#ifdef SHOCK
      sprintf(fname,"%s.minmax.vxavg",        outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.Ekin_long",    outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.Ekin_trans",   outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.Ekin_comp",    outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.shock_shear",  outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.shear_aniso",  outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.pressxy",      outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.pressyz",      outfilename); unlink(fname);
      sprintf(fname,"%s.minmax.presszx",      outfilename); unlink(fname);
#endif
      sprintf(fname,"%s.minmax.dens",         outfilename); unlink(fname);
      sprintf(fname,"%s.tempdist",            outfilename); unlink(fname);
      sprintf(fname,"%s.msqd",                outfilename); unlink(fname);
    }
  }
#ifdef MPI
  MPI_Bcast( &finished, 1, MPI_INT, 0, MPI_COMM_WORLD);
  broadcast_params();
#endif
  return finished;
}

#ifdef MPI

/****************************************************************************
*
*  Broadcast all parameters to other CPUs (MPI only)
*
*****************************************************************************/

void broadcast_params() {

  int i, k;
#ifdef TWOD
  vektor nullv = {0,0};
#else
  vektor nullv = {0,0,0};
#endif

  MPI_Bcast( &ensemble    , 1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &maxwalltime , 1, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( &hyper_threads,1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &watch_int   , 1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &stop_int    , 1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &loop        , 1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &seed        , 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast( &do_maxwell  , 1, MPI_INT,  0, MPI_COMM_WORLD);

  MPI_Bcast( &steps_max   , 1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &steps_min   , 1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &checkpt_int , 1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &eng_int     , 1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &flush_int   , 1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &pic_int     , 1, MPI_INT,  0, MPI_COMM_WORLD);
#if defined(FORCE) || defined(WRITEF)
  MPI_Bcast( &force_int   , 1, MPI_INT,  0, MPI_COMM_WORLD);
#endif
#ifdef WRITEF
  MPI_Bcast( &force_all   , 1, MPI_INT,  0, MPI_COMM_WORLD);
#endif

#ifdef DEBUG
  MPI_Bcast( &force_celldim_divisor, 3, MPI_INT, 0, MPI_COMM_WORLD);
#endif

  MPI_Bcast( &dist_int,              1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dist_dim,            DIM, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dist_ll,             DIM, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &dist_ur,             DIM, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &dist_Epot_flag,        1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dist_Ekin_flag,        1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dist_Ekin_long_flag,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dist_Ekin_trans_flag,  1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dist_Ekin_comp_flag,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dist_press_flag,       1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dist_pressoff_flag,    1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dist_presstens_flag,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dist_shock_shear_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dist_shear_aniso_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dist_dens_flag,        1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dist_vxavg_flag,       1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &box_from_header,       1, MPI_INT, 0, MPI_COMM_WORLD);

#ifdef TWOD
  /*  MPI_Bcast( &pic_scale   , 2, REAL, 0, MPI_COMM_WORLD); */
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

  MPI_Bcast( &vtypes, 1, MPI_INT, 0, MPI_COMM_WORLD);

#ifdef RELAX
  MPI_Bcast( &ekin_threshold,       1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &fnorm_threshold,      1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &f_max_threshold,      1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &delta_epot_threshold, 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &sscount,              1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &nfc,                  1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
#ifdef FBC
  if (NULL==fbc_forces) {
    fbc_forces = (vektor *) malloc( vtypes * sizeof(vektor) );
    if (NULL==fbc_forces)
      error("Cannot allocate memory for fbc_forces on client.");
  }
  MPI_Bcast( fbc_forces, vtypes * DIM, REAL, 0, MPI_COMM_WORLD);
  if (NULL==fbc_beginforces) {
    fbc_beginforces = (vektor *) malloc( vtypes * sizeof(vektor) );
    if (NULL==fbc_beginforces)
      error("Cannot allocate memory for fbc_beginforces on client.");
  }
  MPI_Bcast( fbc_beginforces, vtypes * DIM, REAL, 0, MPI_COMM_WORLD);
#ifdef RELAX
  MPI_Bcast( &max_fbc_int, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (NULL==fbc_dforces) {
    fbc_dforces = (vektor *) malloc( vtypes * sizeof(vektor) );
    if (NULL==fbc_dforces)
      error("Cannot allocate memory for fbc_dforces on client.");
  }
  MPI_Bcast( fbc_dforces, vtypes * DIM, REAL, 0, MPI_COMM_WORLD);
#else
  if (NULL==fbc_endforces) {
    fbc_endforces = (vektor *) malloc( vtypes * sizeof(vektor) );
    if (NULL==fbc_endforces)
      error("Cannot allocate memory for fbc_endforces on client.");
  }
  MPI_Bcast( fbc_endforces, vtypes * DIM, REAL, 0, MPI_COMM_WORLD);
#endif
 if (NULL==fbc_df) {
    fbc_df = (vektor *) malloc( vtypes * sizeof(vektor) );
    if (NULL==fbc_df)
      error("Cannot allocate memory for fbc_df on client.");
  }
#endif /*FBC*/
#ifdef BEND
   if (NULL==fbc_bforces) {
    fbc_bforces = (vektor *) malloc( vtypes * sizeof(vektor) );
    if (NULL==fbc_bforces)
      error("Cannot allocate memory for fbc_bforces on client.");
  }
  MPI_Bcast( fbc_bforces, vtypes * DIM, REAL, 0, MPI_COMM_WORLD);

  if (NULL==fbc_beginbforces) {
      fbc_beginbforces = (vektor *) malloc( vtypes * sizeof(vektor) );
    if (NULL==fbc_beginbforces)
      error("Cannot allocate memory for fbc_beginbforces on client.");
  }
  MPI_Bcast( fbc_beginbforces, vtypes * DIM, REAL, 0, MPI_COMM_WORLD);
#ifdef RELAX
  MPI_Bcast( &max_bfbc_int, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (NULL==fbc_bdforces) {
    fbc_bdforces = (vektor *) malloc( vtypes * sizeof(vektor) );
    if (NULL==fbc_bdforces)
      error("Cannot allocate memory for fbc_bdforces on client.");
  }
  MPI_Bcast( fbc_bdforces, vtypes * DIM, REAL, 0, MPI_COMM_WORLD);
#else
  if (NULL==fbc_endbforces) {
    fbc_endbforces = (vektor *) malloc( vtypes * sizeof(vektor) );
    if (NULL==fbc_endbforces)
      error("Cannot allocate memory for fbc_endbforces on client.");
  }
  MPI_Bcast( fbc_endbforces, vtypes * DIM, REAL, 0, MPI_COMM_WORLD);
#endif
 if (NULL==fbc_bdf) {
    fbc_bdf = (vektor *) malloc( vtypes * sizeof(vektor) );
    if (NULL==fbc_bdf)
      error("Cannot allocate memory for fbc_bdf on client.");
 }
 if (NULL==bend_forces) {
    bend_forces = (vektor *) malloc( vtypes * sizeof(vektor) );
    if (NULL==bend_forces)
      error("Cannot allocate memory for bend_forces on client.");
  }
#endif /* BEND */
#ifdef ZAPP
  MPI_Bcast( &zapp_threshold,         1, REAL,  0, MPI_COMM_WORLD);
#endif

#ifdef FLAGEDATOMS
 MPI_Bcast( &flagedatomstype,         1, MPI_INT,  0, MPI_COMM_WORLD);
#endif 
#ifdef BEND
    MPI_Bcast( &bend_nmoments,         1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast( &bend_vtype_of_origin,  6, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast( &bend_vtype_of_force,   6, MPI_INT,  0, MPI_COMM_WORLD);
     if (NULL==bend_axis) {
        bend_axis = (vektor *) malloc( bend_nmoments * sizeof(vektor) );
    if (NULL==bend_axis)
      error("Cannot allocate memory for bend_axis vector on client.\n");
    }
     MPI_Bcast( bend_axis, bend_nmoments * DIM, REAL, 0, MPI_COMM_WORLD);
     /* these variables will be calculated and do not need to be communicated here */
    if (NULL==bend_origin) {
        bend_origin = (vektor *) malloc( bend_nmoments * sizeof(vektor) );
    if (NULL==bend_origin)
      error("Cannot allocate memory for bend_origin vector on client.\n");
    }
    if (NULL==bend_cog) {
        bend_cog = (vektor *) malloc( bend_nmoments * sizeof(vektor) );
    if (NULL==bend_cog)
      error("Cannot allocate memory for bend_cog vector on client.\n");
    }
    if (NULL==bend_vec) {
        bend_vec = (vektor *) malloc( bend_nmoments * sizeof(vektor) );
    if (NULL==bend_vec)
      error("Cannot allocate memory for bend_vec vector on client.\n");
    }
#endif /* BEND */



  if (NULL==restrictions) {
    restrictions = (vektor *) malloc( vtypes * sizeof(vektor) );
    if (NULL==restrictions)
      error("Cannot allocate memory for restriction vectors on client.");
  }
  MPI_Bcast( restrictions, vtypes * DIM, REAL, 0, MPI_COMM_WORLD);

  MPI_Bcast( &pbc_dirs    , DIM, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &box_x       , DIM, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( &box_y       , DIM, REAL,     0, MPI_COMM_WORLD);
#ifndef TWOD
  MPI_Bcast( &box_z       , DIM, REAL,     0, MPI_COMM_WORLD);
#endif
  MPI_Bcast( &box_param,    DIM, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &size_per_cpu,   1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &box_unit,       1, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( &ntypes,         1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &ntypepairs,     1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &ntypetriples,   1, MPI_INT,  0, MPI_COMM_WORLD);
  if (NULL==masses) {
    masses = (real *) malloc( ntypes * sizeof(real) );
    if (NULL==masses)
      error("Cannot allocate memory for masses array\n");
  }
  MPI_Bcast( masses, ntypes, REAL,     0, MPI_COMM_WORLD);
  if (NULL==gtypes) {
    gtypes = (int *) malloc( ntypes * sizeof(int) );
    if (NULL==gtypes)
      error("Cannot allocate memory for types array\n");
  }
  MPI_Bcast( gtypes, ntypes, MPI_INT,  0, MPI_COMM_WORLD);
#ifdef NBLIST
  MPI_Bcast( &nbl_margin,    1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &nbl_size,      1, REAL, 0, MPI_COMM_WORLD);
#endif
#ifdef VEC
  MPI_Bcast( &atoms_per_cpu, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
#ifdef EFILTER
  if (NULL==lower_e_pot) {
    lower_e_pot = (real *) calloc(ntypes, sizeof(real));
    if (NULL==lower_e_pot)
      error("Cannot allocate memory for lower_e_pot\n");
  }
  if (NULL==upper_e_pot) {
    upper_e_pot = (real *) calloc(ntypes, sizeof(real));
    if (NULL==upper_e_pot)
      error("Cannot allocate memory for upper_e_pot\n");
  }
  MPI_Bcast( lower_e_pot, ntypes,   REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( upper_e_pot, ntypes,   REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ef_checkpt_int, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

#ifdef NNBR
  if (NULL==lower_nb_cut) {
    lower_nb_cut = (int *) calloc(ntypes, sizeof(int));
    if (NULL==lower_nb_cut)
      error("Cannot allocate memory for lower_nb_cut\n");
  }
  if (NULL==upper_nb_cut) {
    upper_nb_cut = (int *) calloc(ntypes, sizeof(int));
    if (NULL==upper_nb_cut)
      error("Cannot allocate memory for upper_nb_cut\n");
  }
  if (NULL==nb_r2_cut) {
    nb_r2_cut = (real *) calloc(ntypes*ntypes, sizeof(real));
    if (NULL==nb_r2_cut)
      error("Cannot allocate memory for nb_r2_cut\n");
  }
  MPI_Bcast( lower_nb_cut,     ntypes, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( upper_nb_cut,     ntypes, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( nb_r2_cut, ntypes*ntypes,    REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &nb_checkpt_int,       1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

  MPI_Bcast( &timestep    ,   1, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( &temperature ,   1, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( &use_curr_temp,  1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &cpu_dim     , DIM, MPI_INT,  0, MPI_COMM_WORLD);

  MPI_Bcast( &parallel_output, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &parallel_input,  1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &msgbuf_size,     1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &binary_output,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( outfilename,            255, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast( infilename,             255, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast( potfilename,            255, MPI_CHAR, 0, MPI_COMM_WORLD);
#ifdef TTBP
  MPI_Bcast( ttbp_potfilename,       255, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
#ifdef EAM2
  MPI_Bcast( eam2_emb_E_filename,    255, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast( eam2_at_rho_filename,   255, MPI_CHAR, 0, MPI_COMM_WORLD);
#ifdef EEAM
  MPI_Bcast( eeam_mod_E_filename,    255, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
#endif
#ifdef ADP
  MPI_Bcast( adp_upotfile,           255, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast( adp_wpotfile,           255, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
#ifdef MEAM
  MPI_Bcast( meam_emb_E_filename,    255, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast( meam_eldensity_filename,255, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast( meam_t1,             ntypes, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( meam_t2,             ntypes, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( meam_t3,             ntypes, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( meam_f0,             ntypes, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( meam_r0,             ntypes, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( meam_beta0,          ntypes, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( meam_beta1,          ntypes, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( meam_beta2,          ntypes, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( meam_beta3,          ntypes, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( meam_rcut_lin,   ntypepairs, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( meam_deltar_lin, ntypepairs, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( meam_cmin_lin, ntypetriples, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( meam_cmax_lin, ntypetriples, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( meam_e,              ntypes, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( meam_a,              ntypes, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( meam_rho0,           ntypes, REAL,     0, MPI_COMM_WORLD);
  MPI_Bcast( &meam_t_average,          1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &have_pre_embed_pot,      1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &have_potfile,            1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &have_eldensity_file,     1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &have_embed_potfile,      1, MPI_INT,  0, MPI_COMM_WORLD);
#endif

#ifdef TEMPCONTROL
  MPI_Bcast( &end_temp,        1, REAL,    0, MPI_COMM_WORLD);
#endif
  MPI_Bcast( &cellsz,          1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &initsz,          1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &incrsz,          1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &outbuf_size,     1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &inbuf_size,      1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dist_chunk_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

#ifdef AND
  MPI_Bcast( &tempintv, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
#ifdef BER
  MPI_Bcast( &tauber, 1, REAL, 0, MPI_COMM_WORLD);
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
#ifdef DAMP
  MPI_Bcast( &stadium,          3 , REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &stadium2,         3 , REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &center,           3 , REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &damptemp,         1 , REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &delta_finnis,     1 , REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &zeta_0,           1 , REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef NPT
  MPI_Bcast( &xi,                DIM, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &isq_tau_xi,          1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &pressure_ext,      DIM, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &use_curr_pressure,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &pressure_end,      DIM, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &cell_size_tolerance, 1, REAL, 0, MPI_COMM_WORLD);
#endif

#if defined(CORRELATE)
  MPI_Bcast( &correl_omode, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &correl_tmax,  1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &correl_rmax,  1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

#if defined(CORRELATE) || defined(MSQD)
  MPI_Bcast( &correl_start, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &correl_end,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &correl_ts,    1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &correl_int,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &msqd_ntypes,  1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &msqd_vtypes,  1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

#ifdef NMOLDYN
  MPI_Bcast( &nmoldyn_int,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &nmoldyn_veloc, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

#ifdef DSF
  MPI_Bcast( &dsf_int,       1, MPI_INT, 0, MPI_COMM_WORLD);
  if (NULL==dsf_weight) {
    dsf_weight = (real *) malloc( ntypes * sizeof(real) );
    if (NULL==dsf_weight) error("cannot allocate dsf_weight");
  }
  MPI_Bcast( dsf_weight,   ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dsf_nk,        1, MPI_INT, 0, MPI_COMM_WORLD);
  if ((myid>0) && (dsf_nk>0)) {
    dsf_k0   = (int *) malloc( DIM * dsf_nk * sizeof(int) );
    dsf_kdir = (int *) malloc( DIM * dsf_nk * sizeof(int) );
    dsf_kmax = (int *) malloc(       dsf_nk * sizeof(int) );
    if ((NULL==dsf_k0) || (NULL==dsf_kdir) || (NULL==dsf_kmax))
      error("cannot allocate dsf arrays");
  }
  MPI_Bcast( dsf_k0,   DIM*dsf_nk, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( dsf_kdir, DIM*dsf_nk, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( dsf_kmax,     dsf_nk, MPI_INT, 0, MPI_COMM_WORLD);
#endif

#if defined(HC) || defined(NVX)
  MPI_Bcast( &hc_int,        1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &hc_start,      1, MPI_INT,  0, MPI_COMM_WORLD);
#endif
#ifdef HC
  MPI_Bcast( &hc_av_start,   1, MPI_INT,  0, MPI_COMM_WORLD);
#endif
#ifdef NVX
  MPI_Bcast( &hc_nlayers,    1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &hc_count,      1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &hc_heatcurr,   1, REAL,     0, MPI_COMM_WORLD);
#endif
#ifdef LASER
  MPI_Bcast( &laser_rescale_mode,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast( &laser_dir,  DIM , MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &laser_mu,         1, REAL,  0, MPI_COMM_WORLD);
  MPI_Bcast( &laser_delta_temp, 1, REAL,  0, MPI_COMM_WORLD);
  MPI_Bcast( &laser_sigma_e,    1, REAL,  0, MPI_COMM_WORLD);
  MPI_Bcast( &laser_sigma_t,    1, REAL,  0, MPI_COMM_WORLD);
  MPI_Bcast( &laser_t_0,        1, REAL,  0, MPI_COMM_WORLD);
  MPI_Bcast( &laser_sigma_e1,    1, REAL,  0, MPI_COMM_WORLD);
  MPI_Bcast( &laser_sigma_t1,    1, REAL,  0, MPI_COMM_WORLD);
  MPI_Bcast( &laser_t_1,        1, REAL,  0, MPI_COMM_WORLD);
#ifdef LASERYZ
  MPI_Bcast( &laser_sigma_w_y,    1,    REAL,  0, MPI_COMM_WORLD);
  MPI_Bcast( &laser_sigma_w_z,    1,    REAL,  0, MPI_COMM_WORLD);
  MPI_Bcast( &laser_sigma_w0,     1,    REAL,  0, MPI_COMM_WORLD);
  MPI_Bcast( &laser_tem_mode,     3, MPI_INT,  0, MPI_COMM_WORLD);
#endif
#endif
#ifdef PDECAY
  MPI_Bcast( &xipdecay,     1, REAL,  0, MPI_COMM_WORLD);
  MPI_Bcast( &ramp_fraction,1, REAL,  0, MPI_COMM_WORLD);
  MPI_Bcast( &pdecay_mode,  1, REAL,  0, MPI_COMM_WORLD);
#endif
#ifdef TTM
  MPI_Bcast( &fd_g,     1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fd_update_steps,1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &fd_ext,       DIM, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &fd_one_d,       1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &fd_c,           1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fd_gamma,	      1, REAL,	  0, MPI_COMM_WORLD);
  MPI_Bcast( &fd_k,           1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fd_n_timesteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ttm_int,        1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &init_t_el,      1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fix_t_el,	      1, MPI_INT, 0, MPI_COMM_WORLD);
#endif /* TTM */
#ifdef STRESS_TENS
  MPI_Bcast( &press_int    , 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &presstens_ext, DIM*(DIM+1)/2, REAL, 0, MPI_COMM_WORLD);
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
  if (NULL==ninslice) {
    ninslice  = (int *) malloc(nslices*sizeof(int));
    if (NULL==ninslice)
      error("Cannot allocate memory for ninslice vector on client.\n");
  }
  if (NULL==E_kin_ftg) {
    E_kin_ftg = (real *) malloc(nslices*sizeof(real));
    if (NULL==E_kin_ftg)
      error("Cannot allocate memory for E_kin_ftg vector on client.\n");
  }
  if (NULL==gamma_ftg) {
    gamma_ftg = (real *) malloc(nslices*sizeof(real));
    if (NULL==gamma_ftg)
      error("Cannot allocate memory for gamma_ftg vector on client.\n");
    for (i=0;i<nslices;i++)
      gamma_ftg[i] = 0.0;
  }
  MPI_Bcast( gamma_ftg, nslices, REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef FINNIS
  MPI_Bcast( &delta_finnis     , 1, REAL   , 0, MPI_COMM_WORLD);
  MPI_Bcast( &zeta_0           , 1, REAL   , 0, MPI_COMM_WORLD);
#endif
#ifdef GLOK
  MPI_Bcast( &glok_ekin_threshold, 1, REAL, 0, MPI_COMM_WORLD);
#endif
#ifdef MIX
  MPI_Bcast( &glok_mix, 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &glok_mixdec, 1, REAL, 0, MPI_COMM_WORLD);
#endif
#ifdef ADAPTGLOK
  MPI_Bcast( &glok_fmaxcrit, 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &glok_incfac, 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &glok_decfac, 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &glok_maxtimestep, 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &starttimestep, 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &glok_minsteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &min_nPxF, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &glok_int, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
#ifdef RIGID
  MPI_Bcast( &nsuperatoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (NULL==superatom) {
    superatom = (int *) malloc( vtypes * sizeof(int) );
    if (NULL==superatom)
      error("Cannot allocate memory for superatom on client.");
    else
      for (k=0; k<vtypes; k++) superatom[k] = -1;
  }
  MPI_Bcast( superatom, vtypes, MPI_INT, 0, MPI_COMM_WORLD);
  if (NULL==superrestrictions) {
    superrestrictions = (vektor *) malloc( vtypes * sizeof(vektor) );
    if (NULL==superrestrictions)
      error("Cannot allocate memory for superrestrictions on client.");
    else
      for (k=0; k<vtypes; k++) superrestrictions[k] = nullv;
  }
  MPI_Bcast( superrestrictions, vtypes * DIM, REAL, 0, MPI_COMM_WORLD);
  if (NULL==superforce)
    superforce = (vektor *) malloc( vtypes * sizeof(vektor) );
  if (NULL==superforce)
    error("Cannot allocate memory for superforce on client.");
#endif
#ifdef DEFORM
  MPI_Bcast( &max_deform_int,  1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &deform_size,     1, REAL,    0, MPI_COMM_WORLD);
  if (NULL==deform_shift) {
    deform_shift = (vektor *) malloc( vtypes * sizeof(vektor) );
    if (NULL==deform_shift)
      error("Cannot allocate memory for deform_shift on client.");
  }
  MPI_Bcast( deform_shift, vtypes * DIM, REAL, 0, MPI_COMM_WORLD);
  if (NULL==shear_def) {
    shear_def = (int *) malloc( vtypes * sizeof(int) );
    if (NULL==shear_def)
      error("Cannot allocate memory for shear_def on client.");
    for (i=0; i<vtypes; i++) shear_def[i] = 0;
  }
  MPI_Bcast( shear_def, vtypes, MPI_INT, 0, MPI_COMM_WORLD);
  if (NULL==deform_shear) {
    deform_shear = (vektor *) malloc( vtypes * sizeof(vektor) );
    if (NULL==deform_shear)
      error("Cannot allocate memory for deform_shear on client.");
  }
  MPI_Bcast( deform_shear, vtypes * DIM, REAL, 0, MPI_COMM_WORLD);
  if (NULL==deform_base) {
    deform_base = (vektor *) malloc( vtypes * sizeof(vektor) );
    if (NULL==deform_base)
      error("Cannot allocate memory for deform_base on client.");
  }
  MPI_Bcast( deform_base, vtypes * DIM, REAL, 0, MPI_COMM_WORLD);
#endif
#ifdef CYCLE
   MPI_Bcast( &lindef_freq,     1, REAL,    0, MPI_COMM_WORLD);
#endif
#ifdef HOMDEF
  MPI_Bcast( &lindef_size,     1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &lindef_int     , 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &lindef_x,      DIM, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &lindef_y,      DIM, REAL,    0, MPI_COMM_WORLD);
#ifndef TWOD
  MPI_Bcast( &lindef_z,      DIM, REAL,    0, MPI_COMM_WORLD);
#endif
  MPI_Bcast( &shear_module,    1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &bulk_module,     1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &relax_rate,      1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &relax_mode,      1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
#if defined(HOMDEF) || defined(NPT_axial)
  MPI_Bcast( &relax_dirs,   DIM, MPI_INT,  0, MPI_COMM_WORLD);
#endif

#ifdef SHOCK
  MPI_Bcast( &shock_strip, 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &shock_speed, 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &shock_speed_l, 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &shock_speed_r, 1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &shock_incr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &shock_mode,  1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

#ifdef CNA
  MPI_Bcast( &cna_start,       1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &cna_end,         1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &cna_int,         1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &cna_rcut,        1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &cna_ll,          3, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &cna_ur,          3, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &cna_writev,      8, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &cna_write_n,     1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &cna_write_statistics, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &cna_cristv,      4, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &cna_crist_n,     1, MPI_INT, 0, MPI_COMM_WORLD);
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
  MPI_Bcast( &avpos_steps,       1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &avpos_nwrites,     1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &avpos_npwrites,    1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

#ifdef ATDIST
  error("Option ATDIST is not supported under MPI");
#endif

#ifdef DIFFPAT
  error("Option DIFFPAT is not supported under MPI");
#endif

#ifdef ORDPAR
  if (NULL==op_r2_cut) {
    op_r2_cut = (real *) calloc(ntypes*ntypes, sizeof(real));
    if (NULL==op_r2_cut)
      error("Cannot allocate memory for op_r2_cut\n");
  }
  MPI_Bcast( &op_r2_cut, ntypes*ntypes, REAL, 0, MPI_COMM_WORLD);
  if (NULL==op_weight) {
    op_weight = (real *) calloc(ntypes*ntypes, sizeof(real));
    if (NULL==op_weight)
      error("Cannot allocate memory for op_weight\n");
  }
  MPI_Bcast( &op_weight, ntypes*ntypes, REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef CG
  MPI_Bcast( &cg_fr,           1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &cg_reset_int,    1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &linmin_maxsteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &linmin_tol,      1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &linmin_dmax,     1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &linmin_dmin,     1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &cg_glimit,       1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &cg_zeps,         1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &cg_infolevel,    1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &cg_mode,         1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
#ifdef ACG
  MPI_Bcast( &acg_init_alpha,      1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &acg_decfac,     1, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &acg_incfac,     1, REAL,    0, MPI_COMM_WORLD);
#endif

#ifdef SOCKET_IO
  MPI_Bcast( &socket_int, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

#ifdef UNIAX
  MPI_Bcast( &uniax_inert,  1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &uniax_sig,    3, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &uniax_eps,    3, REAL, 0, MPI_COMM_WORLD);
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
  /* Gauss Part of Lennard-Jones-Gauss */
  MPI_Bcast( ljg_eps_lin, ntypepairs, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ljg_r0_lin,   ntypepairs, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ljg_sig_lin,   ntypepairs, REAL, 0, MPI_COMM_WORLD);
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
  MPI_Bcast( &fix_bks, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

#ifdef KIM
  MPI_Bcast( kim_model_name, 255, MPI_CHAR, 0, MPI_COMM_WORLD);
  if (NULL == kim_el_names)
    kim_el_names = (char **)calloc(ntypes, sizeof(char *));
  for (i=0;i<ntypes;i++) {
    if (NULL==kim_el_names[i])
      kim_el_names[i] = (char *)calloc(3, sizeof(char));
    MPI_Bcast(kim_el_names[i], 3, MPI_CHAR, 0, MPI_COMM_WORLD);
  }
#endif

#ifdef FEFL
  /* harmonic potential for Einstein crystal */
  MPI_Bcast( spring_rate, ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &lambda, 1, REAL, 0, MPI_COMM_WORLD);
#endif

#if defined(COVALENT) || defined(NNBR_TABLE)
  MPI_Bcast( &neigh_len, 1, MPI_INT,  0, MPI_COMM_WORLD);
#endif

#ifdef TTBP
  MPI_Bcast( ttbp_constant,  ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ttbp_constant2, 8,      REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ttbp_sp,        ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ttbp_cut,       1,      REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ttbp_vas,       1,      MPI_INT, 0, MPI_COMM_WORLD);
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

#if defined(TERSOFF) || defined(TERSOFFMOD) || defined(BRENNER)

#if defined(TERSOFF2) || defined(TERSOFFMOD2) || defined(BRENNER)
  nvalues = ntypepairs;
#else
  nvalues = ntypes;
#endif
  MPI_Bcast( ters_r_cut, ntypepairs,        REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_r0,    ntypepairs,        REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_a,     ntypepairs,        REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_b,     ntypepairs,        REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_la,    ntypepairs,        REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_mu,    ntypepairs,        REAL, 0, MPI_COMM_WORLD);
#if defined(TERSOFF) || defined(BRENNER)
  MPI_Bcast( ters_chi,   ntypepairs-ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_om,    ntypepairs-ntypes, REAL, 0, MPI_COMM_WORLD);
  /* nvalues is ntypes for TERSOFF and ntypepairs for TERSOFF2 */
  MPI_Bcast( ters_ga,        nvalues,       REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_n,         nvalues,       REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_c,         nvalues,       REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_d,         nvalues,       REAL, 0, MPI_COMM_WORLD);
#else /* TERSOFFMOD */
  MPI_Bcast( ters_eta,       nvalues,       REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_delta,     nvalues,       REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_alpha,     nvalues,       REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_beta,      nvalues,       REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_c1,        nvalues,       REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_c2,        nvalues,       REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_c3,        nvalues,       REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_c4,        nvalues,       REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ters_c5,        nvalues,       REAL, 0, MPI_COMM_WORLD);
#endif
  MPI_Bcast( ters_h,         nvalues,       REAL, 0, MPI_COMM_WORLD);
#endif
#ifdef BRENNER
#endif

#ifdef KEATING
  MPI_Bcast( keating_r_cut,       ntypepairs, REAL, 0,MPI_COMM_WORLD);
  MPI_Bcast( keating_alpha,       ntypepairs, REAL, 0,MPI_COMM_WORLD);
  MPI_Bcast( keating_d,           ntypepairs, REAL, 0,MPI_COMM_WORLD);
  MPI_Bcast( keating_beta, ntypes*ntypepairs, REAL, 0,MPI_COMM_WORLD);
#endif

#if defined(EWALD) || defined(COULOMB) || defined(USEFCS) || defined(SM)
  MPI_Bcast( charge,              ntypes, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &coul_eng,           1,      REAL,    0, MPI_COMM_WORLD);
#endif
#ifdef SM
  MPI_Bcast( &sm_fixed_charges,        1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( sm_chi_0,            ntypes, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( sm_J_0,              ntypes, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( sm_Z,                ntypes, REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( sm_zeta,             ntypes, REAL,    0, MPI_COMM_WORLD);
#endif
#ifdef USEFCS
  MPI_Bcast( &fcs_method,         1,   MPI_INT,    0, MPI_COMM_WORLD);
#ifdef PAIR
  MPI_Bcast( &fcs_near_field_flag,1,   MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_rcut,           1,      REAL,    0, MPI_COMM_WORLD);
#endif
  MPI_Bcast( &fcs_tolerance,      1,      REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_grid_dim,       3,   MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_max_iter,       1,   MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_iter_tolerance, 1,      REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_pepc_eps,       1,      REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_pepc_theta,     1,      REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_pepc_nthreads,  1,   MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_fmm_absrel,     1,   MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_fmm_dcorr,      1,   MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_fmm_do_tune ,   1,   MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_vmg_max_level,  1,   MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_vmg_smooth_steps,1,  MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_vmg_gamma,      1,   MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_vmg_near_field_cells,1,MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_vmg_interpol_order,1,MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_vmg_discr_order,1,   MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_pp3mg_ghosts,   1,   MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_pp3mg_degree,   1,   MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_pp3mg_max_part, 1,   MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_p2nfft_intpol_order, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &fcs_p2nfft_epsI,    1,      REAL,    0, MPI_COMM_WORLD);
#endif
#if defined(EWALD) || defined(COULOMB)
  MPI_Bcast( &ew_kappa,           1,      REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &ew_r2_cut,          1,      REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &ew_kcut,            1,      REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &ew_test,            1,      MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ew_nmax,            1,      MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &coul_res,           1,      REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &coul_begin,         1,      REAL,    0, MPI_COMM_WORLD);
#endif
#ifdef DIPOLE
  MPI_Bcast( &dp_fix,             1,      MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &dp_mix,             1,      REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &dp_tol,             1,      REAL,    0, MPI_COMM_WORLD);
  MPI_Bcast( &dp_self,            1,      REAL,    0, MPI_COMM_WORLD);
  if (NULL==dp_b) {
    dp_b = (real *) malloc( ntypepairs * sizeof(real) );
    if (NULL==dp_b)
      error("Cannot allocate memory for dp_b on client.");
  }
  MPI_Bcast( dp_b,            ntypepairs, REAL,    0, MPI_COMM_WORLD);
  if (NULL==dp_c) {
    dp_c = (real *) malloc( ntypepairs * sizeof(real) );
    if (NULL==dp_c)
      error("Cannot allocate memory for dp_c on client.");
  }
  MPI_Bcast( dp_c,            ntypepairs, REAL,    0, MPI_COMM_WORLD);
  if (NULL==dp_alpha) {
    dp_alpha = (real *) malloc( ntypes * sizeof(real) );
    if (NULL==dp_alpha)
      error("Cannot allocate memory for dp_alpha on client.");
  }
  MPI_Bcast( dp_alpha,            ntypes, REAL,    0, MPI_COMM_WORLD);
#endif
#if defined(DIPOLE) || defined(MORSE)
  if (NULL==ms_D) {
    ms_D = (real *) malloc( ntypepairs * sizeof(real) );
    if (NULL==ms_D)
      error("Cannot allocate memory for ms_D on client.");
  }
  MPI_Bcast( ms_D,            ntypepairs, REAL,    0, MPI_COMM_WORLD);
  if (NULL==ms_gamma) {
    ms_gamma = (real *) malloc( ntypepairs * sizeof(real) );
    if (NULL==ms_gamma)
      error("Cannot allocate memory for ms_gamma on client.");
  }
  MPI_Bcast( ms_gamma,        ntypepairs, REAL,    0, MPI_COMM_WORLD);
  if (NULL==ms_harm_c) {
    ms_harm_c = (real *) malloc( ntypepairs * sizeof(real) );
    if (NULL==ms_harm_c)
      error("Cannot allocate memory for ms_harm_c on client.");
  }
  MPI_Bcast( ms_harm_c,        ntypepairs, REAL,    0, MPI_COMM_WORLD);
  if (NULL==ms_r2_min) {
    ms_r2_min = (real *) malloc( ntypepairs * sizeof(real) );
    if (NULL==ms_r2_min)
      error("Cannot allocate memory for ms_r2_min on client.");
  }
  MPI_Bcast( ms_r2_min,      ntypepairs, REAL,    0, MPI_COMM_WORLD);
  if (NULL==ms_r0) {
    ms_r0 = (real *) malloc( ntypepairs * sizeof(real) );
    if (NULL==ms_r0)
      error("Cannot allocate memory for ms_r0 on client.");
  }
  MPI_Bcast( ms_r0,           ntypepairs, REAL,    0, MPI_COMM_WORLD);
#endif

#ifdef EPITAX
  MPI_Bcast( epitax_rate,     ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( epitax_type,     ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( epitax_mass,     ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( epitax_temp,     ntypes, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &epitax_cutoff,       1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &epitax_maxsteps,     1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &epitax_startstep,    1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &epitax_ctrl,         1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &epitax_height,       1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &epitax_speed,        1, REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef RELAX
  MPI_Bcast( &max_sscount,         1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
#ifdef EXTPOT
  MPI_Bcast( &have_extpotfile,            1, MPI_INT,  0, MPI_COMM_WORLD);
  MPI_Bcast( extpotfilename,            255, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ep_n,               1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ep_key,               1,    MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ep_nind,           1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ep_max_int,         1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( ep_pos,         3*ep_n,    REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ep_vel,         3*ep_n,    REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( ep_dir,         3*ep_n,    REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ep_a,               1,    REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ep_rcut,            1,    REAL, 0, MPI_COMM_WORLD);
#endif

#ifdef CBE
  MPI_Bcast( &num_spus,      1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &num_bufs,      1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &cbe_pot_steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &cbe_pot_max,   1, REAL,    0, MPI_COMM_WORLD);
#endif

  MPI_Bcast(&use_header,1, MPI_INT, 0, MPI_COMM_WORLD);

  /* broadcast integrator to other CPUs */
  switch (ensemble) {
    case ENS_NVE:       move_atoms = move_atoms_nve;       break;
    case ENS_TTM:       move_atoms = move_atoms_ttm;       break;
    case ENS_MIK:       move_atoms = move_atoms_mik;       break;
    case ENS_NVT:       move_atoms = move_atoms_nvt;       break;
    case ENS_NPT_ISO:   move_atoms = move_atoms_npt_iso;   break;
    case ENS_NPT_AXIAL: move_atoms = move_atoms_npt_axial; break;
    case ENS_GLOK:      move_atoms = move_atoms_nve;       break;
    case ENS_FRAC:      move_atoms = move_atoms_frac;      break;
    case ENS_SLLOD:     move_atoms = move_atoms_sllod;     break;
    case ENS_NVX:       move_atoms = move_atoms_nvx;       break;
    case ENS_STM:       move_atoms = move_atoms_stm;       break;
    case ENS_FTG:       move_atoms = move_atoms_ftg;       break;
    case ENS_FINNIS:    move_atoms = move_atoms_finnis;    break;
    case ENS_CG:                                           break;
    default: if (0==myid) error("unknown ensemble in broadcast"); break;
  }

#ifdef LASER
  /* broadcast laser rescaling routine to other CPUs */
  switch (laser_rescale_mode) {
    case 0:       do_laser_rescale = laser_rescale_dummy;  break;
    case 1:       do_laser_rescale = laser_rescale_1;   break;
    case 2:       do_laser_rescale = laser_rescale_2;   break;
    case 3:       do_laser_rescale = laser_rescale_3;   break;
#ifdef TTM
    case 4:       do_laser_rescale = laser_rescale_ttm; break;
#endif
    default: if (0==myid)
               error("unknown laser rescaling mode in broadcast"); break;
  }

#ifdef LASERYZ
  switch ( laser_tem_mode.x ) /* Gauss Laguerre = 0, Hermite = 1 */
      {
	  case 0:
	  {
	    switch ( laser_tem_mode.y )
	    {
		case 0:
		{
		  switch ( laser_tem_mode.z )
		  {
		      case 0: laser_intensity_profile = laser_intensity_profile_laguerre_00; break;
		      case 1: laser_intensity_profile = laser_intensity_profile_laguerre_01; break;
		      case 2: laser_intensity_profile = laser_intensity_profile_laguerre_02; break;
		      case 3: laser_intensity_profile = laser_intensity_profile_laguerre_03; break;
		  }
		} break;

		case 1:
		{
		  switch ( laser_tem_mode.z )
		  {
		      case 0: laser_intensity_profile = laser_intensity_profile_laguerre_10; break;
		      case 1: laser_intensity_profile = laser_intensity_profile_laguerre_11; break;
		      case 2: laser_intensity_profile = laser_intensity_profile_laguerre_12; break;
		      case 3: laser_intensity_profile = laser_intensity_profile_laguerre_13; break;
		  }
		} break;

		case 2:
		{
		  switch ( laser_tem_mode.z )
		  {
		      case 0: laser_intensity_profile = laser_intensity_profile_laguerre_20; break;
		      case 1: laser_intensity_profile = laser_intensity_profile_laguerre_21; break;
		      case 2: laser_intensity_profile = laser_intensity_profile_laguerre_22; break;
		      case 3: laser_intensity_profile = laser_intensity_profile_laguerre_23; break;
		  }
		} break;

		case 3:
		{
		  switch ( laser_tem_mode.z )
		  {
		      case 0: laser_intensity_profile = laser_intensity_profile_laguerre_30; break;
		      case 1: laser_intensity_profile = laser_intensity_profile_laguerre_31; break;
		      case 2: laser_intensity_profile = laser_intensity_profile_laguerre_32; break;
		      case 3: laser_intensity_profile = laser_intensity_profile_laguerre_33; break;
		  }
		} break;

	    }
	  } break;
	  case 1:
	  {
	    switch ( laser_tem_mode.y )
	    {
		case 0:
		{
		  switch ( laser_tem_mode.z )
		  {
		      case 0: laser_intensity_profile = laser_intensity_profile_hermite_00; break;
		      case 1: laser_intensity_profile = laser_intensity_profile_hermite_01; break;
		      case 2: laser_intensity_profile = laser_intensity_profile_hermite_02; break;
		      case 3: laser_intensity_profile = laser_intensity_profile_hermite_03; break;
		  }
		} break;

		case 1:
		{
		  switch ( laser_tem_mode.z )
		  {
		      case 0: laser_intensity_profile = laser_intensity_profile_hermite_10; break;
		      case 1: laser_intensity_profile = laser_intensity_profile_hermite_11; break;
		      case 2: laser_intensity_profile = laser_intensity_profile_hermite_12; break;
		      case 3: laser_intensity_profile = laser_intensity_profile_hermite_13; break;
		  }
		} break;

		case 2:
		{
		  switch ( laser_tem_mode.z )
		  {
		      case 0: laser_intensity_profile = laser_intensity_profile_hermite_20; break;
		      case 1: laser_intensity_profile = laser_intensity_profile_hermite_21; break;
		      case 2: laser_intensity_profile = laser_intensity_profile_hermite_22; break;
		      case 3: laser_intensity_profile = laser_intensity_profile_hermite_23; break;
		  }
		} break;

		case 3:
		{
		  switch ( laser_tem_mode.z )
		  {
		      case 0: laser_intensity_profile = laser_intensity_profile_hermite_30; break;
		      case 1: laser_intensity_profile = laser_intensity_profile_hermite_31; break;
		      case 2: laser_intensity_profile = laser_intensity_profile_hermite_32; break;
		      case 3: laser_intensity_profile = laser_intensity_profile_hermite_33; break;
		  }
		} break;

	    }
	  } break;
      }
#endif /* LASERYZ */

#endif /* LASER */

#ifdef ADA
  MPI_Bcast( &ada_nbr_r2cut,       1, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ada_write_int,  1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ada_crystal_structure,  1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ada_default_type,  1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &ada_latticeConst,  1, REAL, 0, MPI_COMM_WORLD);
#endif
#ifdef NYETENSOR
  MPI_Bcast( &nye_rotationAxis_x, 3, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &nye_rotationAxis_y, 3, REAL, 0, MPI_COMM_WORLD);
  MPI_Bcast( &nye_rotationAxis_z, 3, REAL, 0, MPI_COMM_WORLD);
#endif
}

#endif /* MPI */
