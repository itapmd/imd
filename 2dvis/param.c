#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

/*#include "imd.h"*/

#if defined(__GNUC__) && defined(__STRICT_ANSI__)
extern char *strdup(char *);
#endif

#ifdef __WATCOMC__
#define strcasecmp strcmpi
#endif

/* To do yet: improve checking of bad input files ! */
/* e.g. for forgotten starttemp/endtemp              */
/* clean up prototypes.h                            */

typedef short int integer; 
typedef float real;
typedef char str255[255];

typedef enum ParamType {
  PARAM_STR, PARAM_STRPTR,
  PARAM_INT, PARAM_INT_COPY,
  PARAM_INTEGER, PARAM_INTEGER_COPY,
  PARAM_REAL, PARAM_REAL_COPY
} PARAMTYPE;

int curline; /* number of current line */
int finished=0;

extern int scene_type,colmode,radectyp;
extern char *paramfilename;

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
  int i;

  curline = 0;
  pf = fopen(paramfname,"r");
  if (NULL == pf) {
    perror("getparam");
    exit(10);
  };

  /* set the random number generator seed to the current time in seconds */
  /* this will be superseeded by a fixed value from the parameter file */
  { 
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv,&tz);
  }

  do {
    res=fgets(buffer,1024,pf);
    if (NULL == res) { finished=1; break; }; /* probably EOF reached */
    curline++;
    token = strtok(buffer," \t\n");
    if (NULL == token) continue; /* skip blank lines */
    if (token[0]=='#') continue; /* skip comments */

    if (strcasecmp(token,"scene_type")==0) {
      /* type of scene*/
      getparam("scene_type",&scene_type,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"colmode")==0) {
      /* color mode */
      getparam("colmode",&colmode,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"radectyp")==0) {
      /* radius encodes type ? */
      getparam("radectyp",&radectyp,PARAM_INT,1,1);
    }
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


/* Check for Restart, process options */

  while ((argc > 1) && (argv[1][0] =='-')) {
    switch (argv[1][1]) {
    case 'p':
      if (argv[1][2]=='\0') {
        if (NULL != argv[2]) {
          paramfilename = strdup(argv[2]);
          --argc;
          ++argv;
        };
      }
      else paramfilename = NULL;
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
  if (paramfilename) {
    getparamfile(paramfilename,1);
    check_parameters_complete();
  }

}



