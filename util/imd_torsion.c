/******************************************************************************
*
* imd_torsion -- calculate torsion distribution functions
*
* - a descendant of imd.pair
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/


/******************************************************************************
*
*  some useful macros 
*
******************************************************************************/

/* number of slots in the histograms */
#ifndef SLOTS
#define SLOTS 1000 
#endif

#define SQR(a) (a)*(a)
#define TRUNC (int)

/* scalar product for vectors */
#define SPROD(a,b) (((a).x * (b).x) + ((a).y * (b).y) + ((a).z * (b).z))

/* Dynamically allocated 3D arrray -- half vector version */
#define PTR_3D_V(var,i,j,k,dim) (((var) + \
                                 ((i)*(dim.y)*(dim.z)) + \
                                 ((j)*(dim.z)) + \
                                 (k)))

/* Dynamically allocated 5D arrray -- half vector version */
#define PTR_5D_V(var,n,a,b,c,d,dim) (((var) + \
				      ((n)*(dim.i)*(dim.j)*(dim.k)*(dim.l)) + \
				      ((a)*(dim.j)*(dim.k)*(dim.l)) + \
				      ((b)*(dim.k)*(dim.l)) + \
				      ((c)*(dim.l)) + \
				      (d)))

/* Dynamically allocated 3D arrray -- full vector version */
#define PTR_VV(var,coord,dim) (((var) + \
				   ((coord.x)*(dim.y)*(dim.z)) + \
				   ((coord.y)*(dim.z)) + \
				   (coord.z)))

/* allocation increment */
#define CSTEP   10


/*****************************************************************************
*
*  include files
*
*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#ifdef sgi
#include <bstring.h>
#endif

/*****************************************************************************
*
*  typedefs
*
*****************************************************************************/

#ifdef SINGLE
typedef float real;
#else
typedef double real;
#endif

typedef struct {real x; real y; real z; } vektor;
typedef struct {int  x; int  y;  int z; } ivektor;

typedef struct {int  x; int  y; int z; } ivektor3d;

typedef struct {int n; int  i; int  j; int  k; int l; } ivektor5d;

/* Basic Data Type - The Cell */
typedef struct { vektor *ort;
		 int    *sorte;
		 int         n;
		 int     n_max;
	       } cell;

/* String used for Filenames etc. */
typedef char str255[255];

/* for parameter reading */
typedef enum ParamType {
  PARAM_STR, PARAM_STRPTR,
  PARAM_INT, PARAM_INT_COPY,
  PARAM_INTEGER, PARAM_INTEGER_COPY,
  PARAM_REAL, PARAM_REAL_COPY
} PARAMTYPE;


/*****************************************************************************
*
*  prototypes
*
*****************************************************************************/

int  main (int argc, char **argv);
void read_parameters(int argc, char **argv);
void read_atoms(str255 infilename);
void usage(void);
void getparamfile(char *infile);
int  getparam(char *param_name, void *param, PARAMTYPE ptype, 
              int pnum_min, int pnum_max);
void init_cells(void);
void error(char *mesg);
void move_atom(ivektor cellc, cell *from, int index);
void alloc_cell(cell *thecell, int count);
void do_angle(cell *p, cell *q, cell *r, cell *s, vektor pbc_q, vektor pbc_r, vektor pbc_s);
void calc_angles(void);
void write_histograms(int steps);
ivektor cell_coord(vektor pos);


/*****************************************************************************
*
*  global variables
*
*****************************************************************************/

cell *cell_array;    /* Array of Cells */
ivektor cell_dim;    /* Dimension of above */
int natoms;          /* Total number of atoms */
int ntypes;          /* Total number of different atom types */
str255 progname;     /* Name of current executable argv[0] */
int curline;         /* Number of current line for parameter reading */
str255 error_msg;    /* string for error message */

/* The simulation box and its inverse */
vektor box_x, ibox_x, tbox_x;
vektor box_y, ibox_y, tbox_y;
vektor box_z, ibox_z, tbox_z;

ivektor pbc_dirs = {1,1,1};  /* directions with pbc - default is PBC */

/* Filenames */
str255 infilename;    /* Input File */
str255 outfilename;   /* Output File */
char *paramfilename;  /* parameter file *(

/* The histograms */
real      *histogram;
ivektor5d hist_dim;

/* Bookkeeping */
int  restart;
int  nangles = 0;
real r_max = 1.0;     /* default value */
real r2_cut;

/*****************************************************************************
*
*  main
*
*****************************************************************************/

int main(int argc, char **argv)

{
  int tablesize;
  int n,i,j,k,l;

  /* Read Parameters from parameter file */
  read_parameters(argc,argv);

  tablesize = SLOTS*ntypes*ntypes*ntypes*ntypes*sizeof(real);
  histogram = (real *) malloc(tablesize);
  if (NULL==histogram) error("Cannot allocate memory for histograms.");
  hist_dim.n = SLOTS;
  hist_dim.i = ntypes;
  hist_dim.j = ntypes;
  hist_dim.k = ntypes;
  hist_dim.l = ntypes;

  for (n=0; n<SLOTS; ++n)
    for (i=0; i<ntypes; ++i)
      for (j=0; j<ntypes; ++j)
	for (k=0; k<ntypes; ++k)
	  for (l=0; l<ntypes; ++l)
	    *PTR_5D_V(histogram,n,i,j,k,l,hist_dim) = 0.0;

  r2_cut = SQR(r_max);

  /* Initialize cell data structures */
  init_cells();

  /* Read atoms */
  read_atoms(infilename);

  /* Calculate the torsion angles */
  calc_angles();

  /* Output results */
  write_histograms(restart);

  exit(0);

}


/*****************************************************************************
*
*  init_cell determines the minimal size of the cells such that their
*  diameter is less than the cutoff radius.
*
*  A suitable cell array os then allocated.
*
*****************************************************************************/ 

void init_cells(void)

{

  vektor hx,hy,hz; 
  vektor cell_scale;
  real s1,s2;
  cell *p;
  real det;
  int i,j,k;

  /* Calculate smallest possible cell (i.e. the height==cutoff) */

  /* Height x */
  s1 = SPROD(box_x,box_y)/SPROD(box_y,box_y);
  s2 = SPROD(box_x,box_z)/SPROD(box_z,box_z);

  hx.x = box_x.x - s1 * box_y.x - s2 * box_z.x;
  hx.y = box_x.y - s1 * box_y.y - s2 * box_z.y;
  hx.z = box_x.z - s1 * box_y.z - s2 * box_z.z;

  /* Height y */
  s1 = SPROD(box_y,box_x)/SPROD(box_x,box_x);
  s2 = SPROD(box_y,box_z)/SPROD(box_z,box_z);

  hy.x = box_y.x - s1 * box_x.x - s2 * box_z.x;
  hy.y = box_y.y - s1 * box_x.y - s2 * box_z.y;
  hy.z = box_y.z - s1 * box_x.z - s2 * box_z.z;

  /* Height z */
  s1 = SPROD(box_z,box_x)/SPROD(box_x,box_x);
  s2 = SPROD(box_z,box_y)/SPROD(box_y,box_y);

  hz.x = box_z.x - s1 * box_x.x - s2 * box_y.x;
  hz.y = box_z.y - s1 * box_x.y - s2 * box_y.y;
  hz.z = box_z.z - s1 * box_x.z - s2 * box_y.z;

  /* Scaling factors box/cell */
  cell_scale.x = sqrt( (double)(r2_cut / SPROD(hx,hx)) );
  cell_scale.y = sqrt( (double)(r2_cut / SPROD(hy,hy)) );
  cell_scale.z = sqrt( (double)(r2_cut / SPROD(hz,hz)) );

  printf("Minimal cell size: \n");
  printf("\t ( %f %f %f ) \n\t ( %f %f %f ) \n\t ( %f %f %f )\n",
     box_x.x * cell_scale.x, box_x.y * cell_scale.x, box_x.z * cell_scale.x,
     box_y.x * cell_scale.y, box_y.y * cell_scale.y, box_y.z * cell_scale.y,
     box_z.x * cell_scale.z, box_z.y * cell_scale.z, box_z.z * cell_scale.z);

  /* Set up Cell-Array */
  cell_dim.x = (int) ( 1.0 / cell_scale.x );
  cell_dim.y = (int) ( 1.0 / cell_scale.y );
  cell_dim.z = (int) ( 1.0 / cell_scale.z );
  
  /* If an integer number of cells does not fit exactly into the box, the
     cells are enlarged accordingly */
  cell_scale.x = 1.0 / cell_dim.x;
  cell_scale.y = 1.0 / cell_dim.y;
  cell_scale.z = 1.0 / cell_dim.z;

  printf("Actual cell size: \n");
  printf("\t ( %f %f %f ) \n\t ( %f %f %f ) \n\t ( %f %f %f )\n",
     box_x.x * cell_scale.x, box_x.y * cell_scale.x, box_x.z * cell_scale.x,
     box_y.x * cell_scale.y, box_y.y * cell_scale.y, box_y.z * cell_scale.y,
     box_z.x * cell_scale.z, box_z.y * cell_scale.z, box_z.z * cell_scale.z);

  printf("Cell array dimensions: %d %d %d\n",cell_dim.x,cell_dim.y,cell_dim.z);
  cell_array = (cell *) malloc(
		     cell_dim.x * cell_dim.y * cell_dim.z * sizeof(cell));
  if (NULL == cell_array) error("Cannot allocate memory for cells");

  /* Initialize cells */
  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
      for (k=0; k < cell_dim.z; ++k) 
      {
	p = PTR_3D_V(cell_array, i, j, k, cell_dim);

	p->n_max=0;
	alloc_cell(p, CSTEP);
  };

  /* To determine the cell into which a given particle belongs, we
     have to transform the cartesian coordinates of the particle into
     the coordinate system that is spanned by the vectors of the box
     edges. This yields coordinates in the interval [0..1] that are
     multiplied bye global_cell_dim to get the cells index.

     Here we calculate the transformation matrix. */

  /* Calculate inverse of coordinate matrix */

  /* Determinant first */
  det = box_x.y * box_y.z * box_z.x +
        box_x.z * box_y.x * box_z.y +
	box_x.x * box_y.y * box_z.z -
        box_x.z * box_y.y * box_z.x -
        box_x.x * box_y.z * box_z.y -
	box_x.y * box_y.x * box_z.z;
  if ( 0 == det) error("Box Edges are paralell.");

  /* Inverse second */
  ibox_x.x = ( box_y.y * box_z.z - box_y.z * box_z.y ) / det;
  ibox_x.y = ( box_x.z * box_z.y - box_x.y * box_z.z ) / det;
  ibox_x.z = ( box_x.y * box_y.z - box_x.z * box_y.y ) / det;

  ibox_y.x = ( box_y.z * box_z.x - box_y.x * box_z.z ) / det;
  ibox_y.y = ( box_x.x * box_z.z - box_x.z * box_z.x ) / det;
  ibox_y.z = ( box_x.z * box_y.x - box_x.x * box_y.z ) / det;

  ibox_z.x = ( box_y.x * box_z.y - box_y.y * box_z.x ) / det;
  ibox_z.y = ( box_x.y * box_z.x - box_x.x * box_z.y ) / det;
  ibox_z.z = ( box_x.x * box_y.y - box_x.y * box_y.x ) / det;

  /* Transpose */
  tbox_x.x = ibox_x.x;
  tbox_x.y = ibox_y.x;
  tbox_x.z = ibox_z.x;

  tbox_y.x = ibox_x.y;
  tbox_y.y = ibox_y.y;
  tbox_y.z = ibox_z.y;

  tbox_z.x = ibox_x.z;
  tbox_z.y = ibox_y.z;
  tbox_z.z = ibox_z.z;

}


/*****************************************************************************
*
*  cell_coord gives the (integral) cell_coordinates of a position
*
*****************************************************************************/

ivektor cell_coord(vektor ort)

{
  ivektor coord;

  /* Map positions to boxes */
  coord.x = (int) TRUNC( cell_dim.x * SPROD(ort,tbox_x) );
  coord.y = (int) TRUNC( cell_dim.y * SPROD(ort,tbox_y) );
  coord.z = (int) TRUNC( cell_dim.z * SPROD(ort,tbox_z) );

  /* Roundoff errors put atoms slightly out of the simulation cell */
  /* Great mess, needs more investigation */
  coord.x = coord.x <   0         ?             0 : coord.x;
  coord.x = coord.x >= cell_dim.x ? cell_dim.x -1 : coord.x;
  coord.y = coord.y <   0         ?             0 : coord.y;
  coord.y = coord.y >= cell_dim.y ? cell_dim.y -1 : coord.y;
  coord.z = coord.z <   0         ?             0 : coord.z;
  coord.z = coord.z >= cell_dim.z ? cell_dim.z -1 : coord.z;

  return(coord);
}


/******************************************************************************
*
*  move_atoms moves an atom from one cell to another
*
******************************************************************************/

void move_atom(ivektor cellc, cell *from, int index)

{
  cell *to;

  to = PTR_VV(cell_array,cellc,cell_dim);

  /* Check the parameters */
  if ((0 > index) || (index >= from->n)) 
    error("move_atom: index argument out of range.");
  
  /* See if we need some space */
  if (to->n >= to->n_max) alloc_cell(to,to->n_max+CSTEP);

  /* Got some space, move atom */
  to->ort[to->n].x = from->ort[index].x; 
  to->ort[to->n].y = from->ort[index].y; 
  to->ort[to->n].z = from->ort[index].z; 
  to->sorte[to->n] = from->sorte[index]; 
  ++to->n;

  /* Delete atom in original cell */
  --from->n;
  if (0 < from->n) {
    from->ort[index].x = from->ort[from->n].x; 
    from->ort[index].y = from->ort[from->n].y; 
    from->ort[index].z = from->ort[from->n].z; 
    from->sorte[index] = from->sorte[from->n]; 
  };

}


/******************************************************************************
*
*  alloc_cell allocates memory for a cell
*
******************************************************************************/

void alloc_cell(cell *thecell, int count)

{
  cell new;
  int i;

  new.ort    = (vektor *) malloc( count * sizeof(vektor));
  new.sorte  = (int    *) malloc( count * sizeof(int   ));

  if ((NULL==new.ort) || (NULL==new.sorte))
    error("Cannot allocate memory for cell.");

  if (0 == thecell->n_max) {
    /* cell is just initialized */
    thecell->n = 0;
  } else {
    if ( count >= thecell->n_max ) {
      /* cell is enlarged, copy data from old to new location */
      memcpy(new.ort   , thecell->ort,   thecell->n * sizeof(vektor));
      memcpy(new.sorte , thecell->sorte, thecell->n * sizeof(int)   );
     } else {
      /* cell is shrinking, data gets invalid */
      thecell->n = 0;
    };
    free(thecell->ort);
    free(thecell->sorte);
   };

  /* set pointers accordingly */
  thecell->ort    = new.ort;
  thecell->sorte  = new.sorte;
  thecell->n_max =  count;

}


/******************************************************************************
*
*   read_parameters reads the command line parameters, 
*   and the parameter file given on the command line.
*
******************************************************************************/

void read_parameters(int argc,char **argv)

{
  str255 tmpname;
  str255 fname;
  str255 parfilename;
  FILE *infile;
  FILE *tmpfile;

/* Check for Restart, process options */

  strcpy(progname,argv[0]);
  while ((argc > 1) && (argv[1][0] =='-')) {
    switch (argv[1][1]) {
      /* r - restart */
    case 'r':
      restart = atoi(&argv[1][2]);
      break;
      /* e - maximun radius */
    case 'e':
      r_max = atof(&argv[1][2]);
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
      break;
    default:
      printf("Illegal option %s \n",argv[1]);
      usage();
      exit(-1);
    }
    ++argv;
    --argc;
  };

  getparamfile(paramfilename);

  /* Get restart parameters if restart */
  if (0 != restart) {
    sprintf(fname,"%s.%d.itr",outfilename,restart);
    getparamfile(fname);
  };

  if (0!=restart) sprintf(infilename,"%s.%u",outfilename,restart);

} 


/*****************************************************************************
*
*  read one parameter (routine from imd_param.c)
*
*****************************************************************************/

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
      error("string expected");
    }
    else strncpy((char *)param,str,pnum_max);
    numread++;
  }
  else if (ptype == PARAM_STRPTR) {
    str = strtok(NULL," \t\n");
    if (str == NULL) {
      fprintf(stderr,"parameter for %s missing in line %u\n",
              param_name,curline);
      error("string expected");
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
        sprintf(errmsg,"integer vector of dim %u expected",
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
        sprintf(errmsg,"integer vector of dim %u expected",
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
        sprintf(errmsg,"integer vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      }
      else ((int *)param)[i] = atoi(str);
      numread++;
    }
    for (i=pnum_min; i<pnum_max; i++) {
      if ((str = strtok(NULL," \t\n")) != NULL) {
        ((int *)param)[i] = atoi(str);
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
        sprintf(errmsg,"integer vector of dim %u expected",
                (unsigned)pnum_min);
        error(errmsg);
      };
      ((int *)param)[i] = (int)ival;
    }
  }
  else if (ptype == PARAM_REAL) {
    for (i=0; i<pnum_min; i++) {
      str = strtok(NULL," \t\n");
      if (str == NULL) {
        fprintf(stderr,"parameter for %s missing in line %u\n",
                param_name,curline);
        sprintf(errmsg,"real vector of dim %u expected",
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
        sprintf(errmsg,"real vector of dim %u expected",
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
*  read tag-based parameter file
*  lines beginning with comment characters '#' or blank lines are skipped
*
*****************************************************************************/

void getparamfile(char *paramfname)
{
  FILE *pf;
  char buffer[1024];
  char *token;
  char *res;

  curline = 0;
  pf = fopen(paramfname,"r");
  if (NULL == pf) {
    sprintf(error_msg,"Cannot open parameter file %s",paramfname);
    error(error_msg);
  };

  do {
    res=fgets(buffer,1024,pf);
    if (NULL == res) { break; }; /* probably EOF reached */
    curline++;
    token = strtok(buffer," \t\n");
    if (NULL == token) continue; /* skip blank lines */
    if (token[0]=='#') continue; /* skip comments */

    if (strcasecmp(token,"coordname")==0) {
      /* file name for atom coordinate input data */
      getparam("coordname",infilename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"outfiles")==0) {
      /* output file basename */
      getparam("outfiles",outfilename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"box_x")==0) {
      /* 'x' or first vector for box */
      getparam("box_x",&box_x,PARAM_REAL,3,3);
    }
    else if (strcasecmp(token,"box_y")==0) {
      /* 'y' or second vector for box */
      getparam("box_y",&box_y,PARAM_REAL,3,3);
    }
    else if (strcasecmp(token,"box_z")==0) {
      /* 'z' or third vector for box */
      getparam("box_z",&box_z,PARAM_REAL,3,3);
    }
    else if (strcasecmp(token,"pbc_dirs")==0) {
      /* directions with periodic boundary conditions */
      getparam("pbc_dirs",&pbc_dirs,PARAM_INT,3,3);
    }
    else if (strcasecmp(token,"ntypes")==0) {
      /* number of atom types */
      getparam("ntypes",&ntypes,PARAM_INT,1,1);
    }
  } while (!feof(pf));
  fclose(pf);

} /* getparamfile */


/******************************************************************************
*
*  Usage -- educate users
*
******************************************************************************/

void usage(void)

{ 
  printf("%s [-r<nnn>] [-e<nnn>] [-p paramter-file]\n",progname); 
  exit(1); 
}

/******************************************************************************
*
*  error -- Complain and abort
*
******************************************************************************/

void error(char *mesg)

{
  printf("Error: %s\n",mesg);
  exit(2);
}


/******************************************************************************
*
*  read_atoms - reads atoms and velocities into the cell-array
*
*  The file format is flat ascii, one atom per line, lines beginning
*  with '#' denote comments. Each line consists of
*
*  number type mass x y z [rest]
*
*  where
*
*  number   is an arbitrary integer number assigned to each atom
*  type     is the atom's index to the potenital table
*  mass     is the mass of the atom
*  x,y,z    are the atom's coordinates
*  rest     is ignored until end of line
*
******************************************************************************/

void read_atoms(str255 infilename)

{
  cell *input;
  FILE *infile;
  char buf[512];
  int p;
  vektor pos;
  real m;
  int s,n;
  cell *to;
  ivektor cellc;

  /* we first try the old checkpoint name, then the new */
  infile = fopen(infilename,"r");
  if ((NULL==infile) && (restart!=0)) {
    infilename = strcat(infilename,".chkpt");
    infile = fopen(infilename,"r");
  }
  if (NULL==infile) {
    sprintf(error_msg,"Cannot open atoms file %s",infilename);
    error(error_msg);
  }

  /* Set up 1 atom input cell */

  input = (cell *) malloc(sizeof(cell));
  if (0==input) error("Cannot allocate input cell.");
  input->n_max=0;
  alloc_cell(input, 1);

  natoms=0;
  
  /* Read the input file line by line */
  while(!feof(infile)) {

    buf[0] = (char) NULL;
    fgets(buf,sizeof(buf),infile);
    while ('#'==buf[1]) fgets(buf,sizeof(buf),infile); /* eat comments */

#ifdef SINGLE
    p = sscanf(buf,"%d %d %f %f %f %f",&n,&s,&m,&pos.x,&pos.y,&pos.z);
#else
    p = sscanf(buf,"%d %d %lf %lf %lf %lf",&n,&s,&m,&pos.x,&pos.y,&pos.z);
#endif
    if (p>0) {
      ++natoms;
      input->n = 1;
      input->sorte[0] = s;
      input->ort[0] = pos;
      cellc = cell_coord(pos);
      move_atom(cellc, input, 0);
    };
  };
  fclose(infile);  

}


/******************************************************************************
*
*  write_histograms writes histogram to *.torsion file
*
******************************************************************************/

void write_histograms(int steps)

{
  FILE *out;
  str255 fname;
  int n,i,j,k,l;
  real phi;
  real f;

  if (0 == restart)
    sprintf(fname,"%s.torsion",infilename);
  else
    sprintf(fname,"%s.%u.torsion",outfilename,restart);

  out = fopen(fname,"w");
  if (NULL == out) error("Cannot open histograms file.");

  for (n = 0; n < SLOTS; ++n) {
    phi = ((float) n / SLOTS * 180);
    fprintf(out,"%f ", phi);
    for (i = 0; i < ntypes; ++i)
      for (j = i; j < ntypes; ++j)
	for (k = 0; k < ntypes; ++k)
	  for (l = (i != j ? 0 : k); l < ntypes; ++l)
	    if (nangles != 0)
	      fprintf(out,"%f ",*PTR_5D_V(histogram,n,i,j,k,l,hist_dim)/nangles);
	    else 
	      error("No neighbouring atoms in this distance range.");
    fprintf(out,"\n");
  };

  printf("%d torsion angles computed\n", nangles);

  fclose(out);

}


/******************************************************************************
*
*  do_angle calulates the torsion angles for atoms in four cells
*
******************************************************************************/

void do_angle(cell *p, cell *q, cell *r, cell *s, vektor pbc_q, vektor pbc_r, vektor pbc_s)

{
  int i, j, k, l;
  int ang;
  int temp;
  vektor d_ij, d_ik, d_jk, d_il, d_jl, d_kl;
  vektor v_k, v_l;
  real radius_ij, radius_ik, radius_jk, radius_il, radius_jl, radius_kl;
  real phi;
  real betrag_v_k, betrag_v_l;
  double sprod;
  int p_typ, q_typ, r_typ, s_typ;
  
  /*
                      l
                     /
                    /
             i-----j
            /
           /
          k

	  */

  /* For each atom in cell p */
  for (i = 0;i < p->n; ++i) 
    /* For each atom in neighbouring cell q */
    for (j = (p == q ? i+1 : 0); j < q->n; ++j) 
      /* For each atom in neighbouring cell r of p */
      for (k = 0; k < r->n; ++k) 
	/* for each atom in neighbouring cell s of q */
	for (l = 0; l < s->n; ++l)
	  {

	    /* Calculate distance vectors */
	    d_ij.x = q->ort[j].x - p->ort[i].x + pbc_q.x;
	    d_ij.y = q->ort[j].y - p->ort[i].y + pbc_q.y;
	    d_ij.z = q->ort[j].z - p->ort[i].z + pbc_q.z;

	    d_ik.x = r->ort[k].x - p->ort[i].x + pbc_r.x;
	    d_ik.y = r->ort[k].y - p->ort[i].y + pbc_r.y;
	    d_ik.z = r->ort[k].z - p->ort[i].z + pbc_r.z;

	    d_jl.x = s->ort[l].x - q->ort[j].x + pbc_s.x;
	    d_jl.y = s->ort[l].y - q->ort[j].y + pbc_s.y;
	    d_jl.z = s->ort[l].z - q->ort[j].z + pbc_s.z;

	    d_jk.x = d_ik.x - d_ij.x;
	    d_jk.y = d_ik.y - d_ij.y;
	    d_jk.z = d_ik.z - d_ij.z;

	    d_il.x = d_ij.x + d_jl.x;
	    d_il.y = d_ij.y + d_jl.y;
	    d_il.z = d_ij.z + d_jl.z;

	    d_kl.x = d_jl.x - d_jk.x;
	    d_kl.y = d_jl.y - d_jk.y;
	    d_kl.z = d_jl.z - d_jk.z;

	    radius_ij = sqrt( (double)(SPROD(d_ij,d_ij)) );
	    radius_ik = sqrt( (double)(SPROD(d_ik,d_ik)) );
	    radius_jl = sqrt( (double)(SPROD(d_jl,d_jl)) );
	    radius_jk = sqrt( (double)(SPROD(d_jk,d_jk)) );
	    radius_il = sqrt( (double)(SPROD(d_il,d_il)) );
	    radius_kl = sqrt( (double)(SPROD(d_kl,d_kl)) );

	    /* v_k = d_ik x d_jk */
	    v_k.x = d_ik.y * d_jk.z - d_ik.z * d_jk.y; 
	    v_k.y = d_ik.z * d_jk.x - d_ik.x * d_jk.z;
	    v_k.z = d_ik.x * d_jk.y - d_ik.y * d_jk.x;

	    betrag_v_k = sqrt( (double)(SPROD(v_k,v_k)) );
	    
	    /* v_l = d_il x d_jl */
	    v_l.x = d_il.y * d_jl.z - d_il.z * d_jl.y; 
	    v_l.y = d_il.z * d_jl.x - d_il.x * d_jl.z;
	    v_l.z = d_il.x * d_jl.y - d_il.y * d_jl.x;

	    betrag_v_l = sqrt( (double)(SPROD(v_l,v_l)) );

	    /* Calculate torsion angles */
	    if ( (radius_ij < r_max) && (radius_ik < r_max) && (radius_jl < r_max) && (radius_ij*radius_ik*radius_jl*radius_jk*radius_il*radius_kl > 0.0) && (betrag_v_k > 0.0) && (betrag_v_l > 0.0) ) 
	      {
		++nangles;
		sprod = (double)( SPROD(v_k,v_l)/( betrag_v_k * betrag_v_l ));
		if (sprod < -1.0) sprod = -1.0;
		  
		phi = (double) (acos(sprod));
  
		ang = (int) ( SLOTS * phi / 3.141592654 );

		p_typ = p->sorte[i];
		q_typ = q->sorte[j];
		r_typ = r->sorte[k];
		s_typ = s->sorte[l];


		if (p_typ == q_typ) { 
		  if ( r_typ > s_typ ) {
		    temp = s_typ; s_typ = r_typ; r_typ = temp;
		  }
		}
		else if (p_typ > q_typ) {
		  temp = q_typ; q_typ = p_typ; p_typ = temp;
		  temp = s_typ; s_typ = r_typ; r_typ = temp;
		}

		if ((ang >= 0) && (ang < SLOTS))
		  ++*PTR_5D_V(histogram, ang , p_typ, q_typ, r_typ, s_typ, hist_dim);
	      } 
	  }

}


/******************************************************************************
*
*  calc_angles calulates the torsion angles for all atoms
*
******************************************************************************/

void calc_angles(void)

{
  cell *p, *q,* r, *s;
  int i1, i2, i3, j1, j2, j3;
  int k1, k2, k3, l1, l2, l3;
  int a1, a2, a3, b1, b2, b3, c1, c2, c3;
  vektor pbc_q, pbc_r, pbc_s;

  /* for each cell p */
  for (i1 = 0; i1 < cell_dim.x; ++i1)
    for (i2 = 0; i2 < cell_dim.y; ++i2)
      for (i3 = 0; i3 < cell_dim.z; ++i3)

	/* For half of the neighbours q of this cell */
	for (j1 = 0; j1 <= 1; ++j1)
	  for (j2 = -j1; j2 <= 1; ++j2)
	    for (j3= (j1 == 0 ? -j2 : -j1); j3 <= 1; ++j3)

	      /* For the neighbours r of cell p */
	      for (k1 = -1; k1 <= 1; ++k1)
		for (k2 = -1; k2 <= 1; ++k2)
		  for (k3 = -1; k3 <= 1; ++k3)

		    /* for the neighbours s of cell q */
		    for (l1 = -1; l1 <= 1; ++l1)
		      for (l2 = -1; l2 <= 1; ++l2)
			for (l3 = -1; l3 <= 1; ++l3)

			  {
			    /* Given cell */

			    p = PTR_3D_V(cell_array,i1,i2,i3,cell_dim);

			    /* Calculate Indices of Neighbour */
			    a1 = i1 + j1;  pbc_q.x = 0;
			    a2 = i2 + j2;  pbc_q.y = 0;
			    a3 = i3 + j3;  pbc_q.z = 0;

			    /* deal with periodic boundary conditions if necessary */
			    if (a1 < 0) {
			      if (pbc_dirs.x == 1) {
				a1 = cell_dim.x-1; 
				pbc_q.x -= box_x.x;      
				pbc_q.y -= box_x.y;
				pbc_q.z -= box_x.z;

			      } else continue;
			    }
			    if (a2 < 0) {
			      if (pbc_dirs.y == 1) {
				a2 = cell_dim.y-1;
				pbc_q.x -= box_y.x;      
				pbc_q.y -= box_y.y;
				pbc_q.z -= box_y.z;

			      } else continue;
			    }
			    if (a3 < 0) {
			      if (pbc_dirs.z == 1) {
				a3 = cell_dim.z-1;
				pbc_q.x -= box_z.x;      
				pbc_q.y -= box_z.y;
				pbc_q.z -= box_z.z;
			      } else continue;
			    }
			    if (a1 > cell_dim.x-1) {
			      if (pbc_dirs.x == 1) {
				a1 = 0; 
				pbc_q.x += box_x.x;      
				pbc_q.y += box_x.y;
				pbc_q.z += box_x.z;
			      } else continue;
			    }
			    if (a2 > cell_dim.y-1) {
			      if (pbc_dirs.y == 1) {
				a2 = 0; 
				pbc_q.x += box_y.x;      
				pbc_q.y += box_y.y;
				pbc_q.z += box_y.z;
			      } else continue;
			    }
			    if (a3 > cell_dim.z-1) {
			      if (pbc_dirs.z == 1) {
				a3 = 0; 
				pbc_q.x += box_z.x;      
				pbc_q.y += box_z.y;
				pbc_q.z += box_z.z;
			      } else continue;
			    }

			    /* Neighbour cell q (note that p==q ist possible) */

			    q = PTR_3D_V(cell_array,a1,a2,a3,cell_dim);

			    /* Calculate Indices of neighbour r of p */
			    b1 = i1 + k1; pbc_r.x = 0;
			    b2 = i2 + k2; pbc_r.y = 0;
			    b3 = i3 + k3; pbc_r.z = 0;

			    /* deal with periodic boundary conditions if necessary */
			    if (b1 < 0) {
			      if (pbc_dirs.x == 1) {
				b1 = cell_dim.x-1; 
				pbc_r.x -= box_x.x;      
				pbc_r.y -= box_x.y;
				pbc_r.z -= box_x.z;
			      } else continue;
			    }
			    if (b2 < 0) {
			      if (pbc_dirs.y == 1) {
				b2 = cell_dim.y-1;
				pbc_r.x -= box_y.x;      
				pbc_r.y -= box_y.y;
				pbc_r.z -= box_y.z;
			      } else continue;
			    }
			    if (b3 < 0) {
			      if (pbc_dirs.z == 1) {
				b3 = cell_dim.z-1;
				pbc_r.x -= box_z.x;      
				pbc_r.y -= box_z.y;
				pbc_r.z -= box_z.z;
			      } else continue;
			    }
			    if (b1 > cell_dim.x-1) {
			      if (pbc_dirs.x == 1) {
				b1 = 0; 
				pbc_r.x += box_x.x;      
				pbc_r.y += box_x.y;
				pbc_r.z += box_x.z;
			      } else continue; 
			    }
			    if (b2 > cell_dim.y-1) {
			      if (pbc_dirs.y == 1) {
				b2 = 0; 
				pbc_r.x += box_y.x;      
				pbc_r.y += box_y.y;
				pbc_r.z += box_y.z;
			      } else continue; 
			    }
			    if (b3 > cell_dim.z-1) {
			      if (pbc_dirs.z == 1) {
				b3 = 0; 
				pbc_r.x += box_z.x;      
				pbc_r.y += box_z.y;
				pbc_r.z += box_z.z;
			      } else continue;
			    } 

			    /* Neighbour cell r */

			    r = PTR_3D_V(cell_array, b1, b2, b3, cell_dim);

			    /* Calculate Indices of neighbour s of q */
			    c1 = j1 + l1; pbc_s.x = 0;
			    c2 = j2 + l2; pbc_s.y = 0;
			    c3 = j3 + l3; pbc_s.z = 0;

			    /* deal with periodic boundary conditions if necessary */
			    if (c1 < 0) {
			      if (pbc_dirs.x == 1) {
				c1 = cell_dim.x-1; 
				pbc_s.x -= box_x.x;      
				pbc_s.y -= box_x.y;
				pbc_s.z -= box_x.z;
			      } else continue;
			    }
			    if (c2 < 0) {
			      if (pbc_dirs.y == 1) {
				c2 = cell_dim.y-1;
				pbc_s.x -= box_y.x;      
				pbc_s.y -= box_y.y;
				pbc_s.z -= box_y.z;
			      } else continue;
			    }
			    if (c3 < 0) {
			      if (pbc_dirs.z == 1) {
				c3 = cell_dim.z-1;
				pbc_s.x -= box_z.x;      
				pbc_s.y -= box_z.y;
				pbc_s.z -= box_z.z;
			      } else continue;
			    }
			    if (c1 > cell_dim.x-1) {
			      if (pbc_dirs.x == 1) {
				c1 = 0; 
				pbc_s.x += box_x.x;      
				pbc_s.y += box_x.y;
				pbc_s.z += box_x.z;
			      } else continue; 
			    }
			    if (c2 > cell_dim.y-1) {
			      if (pbc_dirs.y == 1) {
				c2 = 0; 
				pbc_s.x += box_y.x;      
				pbc_s.y += box_y.y;
				pbc_s.z += box_y.z;
			      } else continue; 
			    }
			    if (c3 > cell_dim.z-1) {
			      if (pbc_dirs.z == 1) {
				c3 = 0; 
				pbc_s.x += box_z.x;      
				pbc_s.y += box_z.y;
				pbc_s.z += box_z.z;
			      } else continue;
			    } 

			    /* Neighbour cell s */

			    s = PTR_3D_V(cell_array, c1, c2, c3, cell_dim);

			    /* Do the work */
			    do_angle(p,q,r,s,pbc_q,pbc_r,pbc_s);
			  }
  
}



































































