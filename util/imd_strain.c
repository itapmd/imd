/******************************************************************************
*
* imd_strain -- calculate strain tensor
*
* uses cell division routines of imd_pair
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/


/******************************************************************************
*
*  some usful macros 
*
******************************************************************************/

#ifdef TWOD
#define DIM 2
#else
#define DIM 3
#endif

#define SQR(a) (a)*(a)
#define TRUNC (int)

/* Scalar product for vectors */
#ifdef TWOD
#define SPROD(a,b) (((a).x * (b).x) + ((a).y * (b).y))
#else
#define SPROD(a,b) (((a).x * (b).x) + ((a).y * (b).y) + ((a).z * (b).z))
#endif


/* Dynamically allocated 3D arrray -- half vector version */
#define PTR_3D_V(var,i,j,k,dim) (((var) + \
                                 ((i)*(dim.y)*(dim.z)) + \
                                 ((j)*(dim.z)) + \
                                 (k)))

/* Dynamically allocated 3D arrray -- full vector version */
#define PTR_3D_VV(var,coord,dim) (((var) + \
                                 ((coord.x)*(dim.y)*(dim.z)) + \
                                 ((coord.y)*(dim.z)) + \
                                 (coord.z)))

/* Dynamically allocated 2D arrray -- half vector version */
#define PTR_2D_V(var,i,j,dim) (((var) + \
                               ((i)*(dim.y)) + \
                                (j)))
                                

/* Dynamically allocated 2D arrray -- full vector version */
#define PTR_2D_VV(var,coord,dim) (((var) + \
                                 ((coord.x)*(dim.y)) + \
                                  (coord.y)))

#ifdef TWOD
#define PTR_VV   PTR_2D_VV
#else
#define PTR_VV   PTR_3D_VV
#endif

/* Allocation increment */
#define CSTEP   2


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

#ifdef TWOD
typedef struct {real x; real y; } vektor;
typedef struct {int  x; int  y; } ivektor;
#else
typedef struct {real x; real y; real z; } vektor;
typedef struct {int  x; int  y;  int z; } ivektor;
#endif

typedef struct {int  x; int  y; int z; } ivektor3d;

/* Basic Data Type - The Cell */
typedef struct { vektor *ort;
                 vektor *dsp;
                 vektor *strain;
                 vektor *strain_offdia;
                 short  *empty;
		 int         n;
		 int     n_max;
	       } cell;

/* String used for Filenames etc. */
typedef char str255[255];

/* For parameter reading */
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
void read_displacement(str255 infilename);
void usage(void);
void getparamfile(char *infile);
int  getparam(char *param_name, void *param, PARAMTYPE ptype, 
              int pnum_min, int pnum_max);
void init_cells(void);
void error(char *mesg);
void move_atom(ivektor cellc, cell *from, int index);
void alloc_cell(cell *thecell, int count);
void calc_strain(void);
void write_strain(void);
ivektor cell_coord(vektor pos);


/*****************************************************************************
*
*  global variables
*
*****************************************************************************/

cell *cell_array;    /* Array of Cells */
ivektor cell_dim;    /* Dimension of above */
int natoms;          /* Total number of atoms */
str255 progname;     /* Name of current executable argv[0] */
int curline;         /* Number of current line for parameter reading */
str255 error_msg;    /* string for error message */

/* The simulation box and its inverse */
vektor box_x, tbox_x;
vektor box_y, tbox_y;
#ifndef TWOD
vektor box_z, tbox_z;
#endif
#ifdef TWOD
ivektor pbc_dirs = {1,1};    /* directions with pbc - default is PBC */
#else
ivektor pbc_dirs = {1,1,1};  /* directions with pbc - default is PBC */
#endif

/* Filenames */
str255 infilename;    /* Input File */
str255 outfilename;   /* Output File */
char *paramfilename;  /* parameter file */

int restart;
real r_max = 1.0;     /* default value */
real r_cell = 1.0;    /* default value */
real r2_cut;


/*****************************************************************************
*
*  main
*
*****************************************************************************/

int main(int argc, char **argv)

{

  /* Read Parameters from parameter file */
  read_parameters(argc,argv);

  /* Initialize cell data structures */
  init_cells();

  /* Read coordinates and displacement */
  read_displacement(infilename);

  /* Calculate strain tensor */
  calc_strain();

  /* Output strain tensor */
  write_strain();

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

#ifdef TWOD

 /* Height x */
  s1 = SPROD(box_x,box_x) - SPROD(box_x,box_y);
  
  /* Height y */
  s2 = SPROD(box_y,box_y) - SPROD(box_x,box_y);
  
  cell_scale.x = sqrt( (double)(r2_cut / s1) );
  cell_scale.y = sqrt( (double)(r2_cut / s2) );

  printf("Minimal cell size: \n\t ( %f %f ) \n\t ( %f %f ) \n",
	 box_x.x * cell_scale.x, box_x.y * cell_scale.x, 
	 box_y.x * cell_scale.y, box_y.y * cell_scale.y); 

  /* Set up Cell-Array */
  cell_dim.x = (int) ( 1.0 / cell_scale.x );
  cell_dim.y = (int) ( 1.0 / cell_scale.y );
  
  /* If an integer number of cells does not fit exactly into the box, the
     cells are enlarged accordingly */
  cell_scale.x = 1.0 / cell_dim.x;
  cell_scale.y = 1.0 / cell_dim.y;

  printf("Minimal cell size: \n\t ( %f %f ) \n\t ( %f %f ) \n",
	 box_x.x * cell_scale.x, box_x.y * cell_scale.x, 
	 box_y.x * cell_scale.y, box_y.y * cell_scale.y); 

  printf("Cell array dimensions: %d %d\n", cell_dim.x,cell_dim.y);
  cell_array = (cell *) malloc(cell_dim.x * cell_dim.y * sizeof(cell));
  if (NULL == cell_array) error("Cannot allocate memory for cells");

#else /* 3D */

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
  
  cell_array = (cell *) malloc(cell_dim.x * cell_dim.y * cell_dim.z * sizeof(cell)); 
  if (NULL == cell_array) error("Cannot allocate memory for cells");
 
#endif /* TWOD */

  /* Initialize cells */
  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
#ifndef TWOD
      for (k=0; k < cell_dim.z; ++k) 
#endif
      {
#ifdef TWOD
	p = PTR_2D_V(cell_array, i, j,    cell_dim);
#else
	p = PTR_3D_V(cell_array, i, j, k, cell_dim);
#endif
	p->n_max=0;
	alloc_cell(p, CSTEP);
  }

  /* To determine the cell into which a given particle belongs, we
     have to transform the cartesian coordinates of the particle into
     the coordinate system that is spanned by the vectors of the box
     edges. This yields coordinates in the interval [0..1] that are
     multiplied bye global_cell_dim to get the cells index.

     Here we calculate the transformation matrix. */

  /* Calculate inverse of coordinate matrix */

#ifdef TWOD

  /* determinant */
  det = box_x.x * box_y.y - box_x.y * box_y.x;
  if (0==det) error("Box Edges are parallel.");

  /* transposed inverse box */
  tbox_x.x =   box_y.y / det;
  tbox_x.y = - box_y.x / det;
  tbox_y.x = - box_x.y / det;
  tbox_y.y =   box_x.x / det;

#else /* 3D */

  /* determinant */
  det = box_x.y * box_y.z * box_z.x +
        box_x.z * box_y.x * box_z.y +
	box_x.x * box_y.y * box_z.z -
        box_x.z * box_y.y * box_z.x -
        box_x.x * box_y.z * box_z.y -
	box_x.y * box_y.x * box_z.z;
  if (0==det) error("Box Edges are parallel.");

  /* transposed inverse box */
  tbox_x.x = ( box_y.y * box_z.z - box_y.z * box_z.y ) / det;
  tbox_y.x = ( box_x.z * box_z.y - box_x.y * box_z.z ) / det;
  tbox_z.x = ( box_x.y * box_y.z - box_x.z * box_y.y ) / det;

  tbox_x.y = ( box_y.z * box_z.x - box_y.x * box_z.z ) / det;
  tbox_y.y = ( box_x.x * box_z.z - box_x.z * box_z.x ) / det;
  tbox_z.y = ( box_x.z * box_y.x - box_x.x * box_y.z ) / det;

  tbox_x.z = ( box_y.x * box_z.y - box_y.y * box_z.x ) / det;
  tbox_y.z = ( box_x.y * box_z.x - box_x.x * box_z.y ) / det;
  tbox_z.z = ( box_x.x * box_y.y - box_x.y * box_y.x ) / det;

#endif /* not TWOD */

}


/*****************************************************************************
*
*  cell_coord gives the (integral) cell_coorinates of a position
*
*****************************************************************************/

ivektor cell_coord(vektor ort)

{
  ivektor coord;

  /* Map positions to boxes */
  coord.x = (int) TRUNC( cell_dim.x * SPROD(ort,tbox_x) );
  coord.y = (int) TRUNC( cell_dim.y * SPROD(ort,tbox_y) );
#ifndef TWOD
  coord.z = (int) TRUNC( cell_dim.z * SPROD(ort,tbox_z) );
#endif

  /* Roundoff errors put atoms slightly out of the simulation cell */
  /* Great mess, needs more investigation */
  coord.x = coord.x <   0         ?             0 : coord.x;
  coord.x = coord.x >= cell_dim.x ? cell_dim.x -1 : coord.x;
  coord.y = coord.y <   0         ?             0 : coord.y;
  coord.y = coord.y >= cell_dim.y ? cell_dim.y -1 : coord.y;
#ifndef TWOD
  coord.z = coord.z <   0         ?             0 : coord.z;
  coord.z = coord.z >= cell_dim.z ? cell_dim.z -1 : coord.z;
#endif

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
  to->dsp[to->n].x = from->dsp[index].x;
  to->dsp[to->n].y = from->dsp[index].y; 
#ifndef TWOD
  to->ort[to->n].z = from->ort[index].z;
  to->dsp[to->n].z = from->dsp[index].z;
#endif 
  ++to->n;

  /* Delete atom in original cell */
  --from->n;
  if (0 < from->n) {
    from->ort[index].x = from->ort[from->n].x; 
    from->ort[index].y = from->ort[from->n].y;
    from->dsp[index].x = from->dsp[from->n].x; 
    from->dsp[index].y = from->dsp[from->n].y; 
#ifndef TWOD
    from->ort[index].z = from->ort[from->n].z; 
    from->dsp[index].z = from->dsp[from->n].z; 
#endif 
  }

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

  new.ort           = (vektor *) malloc( count * sizeof(vektor));
  new.dsp           = (vektor *) malloc( count * sizeof(vektor));
  new.strain        = (vektor *) malloc( count * sizeof(vektor));
  new.strain_offdia = (vektor *) malloc( count * sizeof(vektor));
  new.empty         = (short *) malloc( count * sizeof(short));

  if ((NULL==new.ort) || (NULL==new.dsp) || (NULL==new.strain) || (NULL==new.strain_offdia) || (NULL==new.empty))
    error("Cannot allocate memory for cell.");

  if (0 == thecell->n_max) {
    /* cell is just initialized */
    thecell->n = 0;
  } else {
    if ( count >= thecell->n_max ) {
      /* cell is enlarged, copy data from old to new location */
      memcpy(new.ort          , thecell->ort,           thecell->n * sizeof(vektor));
      memcpy(new.dsp          , thecell->dsp,           thecell->n * sizeof(vektor));
      memcpy(new.strain       , thecell->strain,        thecell->n * sizeof(vektor));
      memcpy(new.strain_offdia, thecell->strain_offdia, thecell->n * sizeof(vektor));
      memcpy(new.empty        , thecell->empty,         thecell->n * sizeof(short));

     } else {
      /* cell is shrinking, data gets invalid */
      thecell->n = 0;
    }
    free(thecell->ort);
    free(thecell->dsp);
    free(thecell->strain);
    free(thecell->strain_offdia);
    free(thecell->empty);
   }

  /* set pointers accordingly */
  thecell->ort           = new.ort;
  thecell->dsp           = new.dsp;
  thecell->strain        = new.strain; 
  thecell->strain_offdia = new.strain_offdia;
  thecell->empty         = new.empty;
  thecell->n_max         = count;

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

/* Check for process options */

  strcpy(progname,argv[0]);
  while ((argc > 1) && (argv[1][0] =='-')) {
    switch (argv[1][1]) {
      /* r - restart */
    case 'r':
      restart = atoi(&argv[1][2]);
      break;
      /* e - maximun radius */
    case 'e':
      r_max  = atof(&argv[1][2]);
      r_cell = r_max;
      break;           
      /* c - cell width */
    case 'c':
      r_cell = atof(&argv[1][2]);
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

  getparamfile(paramfilename);

  sprintf(infilename, "%s.dsp", outfilename);

  if (0!=restart) sprintf(infilename,"%s.%u.dsp",outfilename,restart);

  /* r_cell >= r_max is required */
  if (r_cell < r_max) error("Cell smaller than cutoff radius!");

  r2_cut = SQR(r_cell);
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
      }
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
      }
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
      }
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
  }

  do {
    res=fgets(buffer,1024,pf);
    if (NULL == res) break; /* probably EOF reached */
    curline++;
    token = strtok(buffer," \t\n");
    if (NULL == token) continue; /* skip blank lines */
    if (token[0]=='#') continue; /* skip comments */

    if (strcasecmp(token,"outfiles")==0) {
      /* output file basename */
      getparam("outfiles",outfilename,PARAM_STR,1,255);
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
    else if (strcasecmp(token,"pbc_dirs")==0) {
      /* directions with periodic boundary conditions */
      getparam("pbc_dirs",&pbc_dirs,PARAM_INT,DIM,DIM);
    }
   
  } while (!feof(pf));
  fclose(pf);

  if ( box_x.x < r_cell || box_y.y < r_cell 
#ifndef TWOD
       || box_z.z < r_cell 
#endif
       )
    error("Cell size greater than box!");

} /* getparamfile */


/******************************************************************************
*
*  Usage -- educate users
*
******************************************************************************/

void usage(void)

{ 
  printf("%s [-r<nnn>] [-e<nnn>] [-c<nnn>] [-p paramter-file]\n",progname); 
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
*  read_displacement - reads atom positions and displacements into the cell-array
*
*  The file format is flat ascii, one atom per line, lines beginning
*  with '#' denote comments. Each line consists of
*
*  x y [z] dx dy [dz] [rest]
*
*  where
*
*  x,y,z    are the atom's coordinates
*  dx,dy,dz are the atom's displacement vector
*  rest     is ignored until end of line
*
******************************************************************************/

void read_displacement(str255 infilename)

{
  cell *input;
  FILE *infile;
  char buf[512];
  int p;
  vektor pos;
  vektor u;
  real m;
  int s;
  cell *to;
  ivektor cellc;

  infile = fopen(infilename,"r");
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
#ifdef TWOD
    p = sscanf(buf,"%f %f %f %f",&pos.x,&pos.y,&u.x,&u.y);
#else
    p = sscanf(buf,"%f %f %f %f %f %f",&pos.x,&pos.y,&pos.z,&u.x,&u.y,&u.z);
#endif
#else
#ifdef TWOD
    p = sscanf(buf,"%lf %lf %lf %lf",&pos.x,&pos.y,&u.x,&u.y);
#else
    p = sscanf(buf,"%lf %lf %lf %lf %lf %lf",&pos.x,&pos.y,&pos.z,&u.x,&u.y,&u.z);
#endif
#endif
    if (p>0) {
      ++natoms;
      input->n      = 1;
      input->ort[0] = pos;
      input->dsp[0] = u;
      cellc         = cell_coord(pos);
      move_atom(cellc, input, 0);
    }
  }
  fclose(infile);  

}


/******************************************************************************
*
*  write_strain writes strain to *.strain file
*
******************************************************************************/

void write_strain(void)

{
  FILE *out;
  str255 fname;
  int i,j,k,l;
  int number=0;
  cell *p; 
 
  sprintf(fname,"%s.strain",infilename);
  

  out = fopen(fname,"w");
  if (NULL == out) error("Cannot open strain file.");

  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
#ifndef TWOD
      for (k=0; k < cell_dim.z; ++k)
#endif
	{
#ifdef TWOD
        p = PTR_2D_V(cell_array,i,j  ,cell_dim);
#else
        p = PTR_3D_V(cell_array,i,j,k,cell_dim);
#endif
	for (l=0;l<p->n; ++l)

	  if (p->empty[l]==0)
	    {
#ifdef TWOD
	  fprintf(out, "%d %f %f %f %f %f %f\n", ++number, p->ort[l].x, p->ort[l].y, p->strain[l].x, p->strain[l].y, p->strain_offdia[l].x);
#else
	fprintf(out, "%d %f %f %f %f %f %f %f %f %f\n", ++number, p->ort[l].x, p->ort[l].y, p->ort[l].z, p->strain[l].x, p->strain[l].y, p->strain[l].z, p->strain_offdia[l].x, p->strain_offdia[l].y, p->strain_offdia[l].z);
#endif
	    }
	}
  fclose(out);

}


/******************************************************************************
*
*  calc_strain calculates the strain tensor for all atoms
*
******************************************************************************/

void calc_strain(void)

{
  cell *p,*q;
  int i,j,k;
  int l,m,n;
  int r,s,t;
  int u,v,w;
  int g,h;
  int number;                          /* number of neighbours */ 
  int totalnumber = 0;                 /* number of neighbours of all atoms */
  int emptynumber = 0;                 /* number of atoms with less than 3 neighbours */ 
  int maxnumber = 0, minnumber = 1000; /* maximal and minimal number of neighbours occuring */
  int num = 4;
  vektor pbc;
  vektor *d, *du, *tmp;
  real a[3][3], b[3][3];
  real radius;
  real det;

  /* Allocate memory for temporary variables */
  d   = (vektor *) malloc( num * sizeof(vektor));
  du  = (vektor *) malloc( num * sizeof(vektor));
  tmp = (vektor *) malloc( num * sizeof(vektor));


  /* Initialization */

  /* for each cell */
  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
#ifndef TWOD
      for (k=0; k < cell_dim.z; ++k)
#endif
	{
#ifdef TWOD
	  p = PTR_2D_V(cell_array,i,j  ,cell_dim);
#else
	  p = PTR_3D_V(cell_array,i,j,k,cell_dim);
#endif
	  /* For each atom in this cell */
	  for (u=0; u < p->n; ++u)
	    { 
	      p->strain[u].x        = 0.0;
	      p->strain[u].y        = 0.0;
#ifndef TWOD
	      p->strain[u].z        = 0.0;
#endif
	      p->strain_offdia[u].x = 0.0;
#ifndef TWOD
	      p->strain_offdia[u].y = 0.0;
	      p->strain_offdia[u].z = 0.0;
#endif
	      p->empty[u]           = 0;
	    }
	}

  /* Compute strain tensor */

  /* For each cell */
  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
#ifndef TWOD
      for (k=0; k < cell_dim.z; ++k)
#endif
	{
#ifdef TWOD
	  p = PTR_2D_V(cell_array,i,j  ,cell_dim);
#else
	  p = PTR_3D_V(cell_array,i,j,k,cell_dim);
#endif

	  /* For each atom in this first cell */
	  for (u=0; u<p->n; ++u)
	    {

	      number = 0;
	      
	      /* For the neighbours of the first cell */
	      for (l=-1; l <= 1; ++l)
		for (m=-1; m <= 1; ++m)
#ifndef TWOD
		  for (n=-1; n <= 1; ++n)
#endif
		    {
	      /* Calculate Indices of Neighbour */
	      r = i+l;  pbc.x = 0;
	      s = j+m;  pbc.y = 0;
#ifndef TWOD
	      t = k+n;  pbc.z = 0;
#endif
              /* Deal with periodic boundary conditions if necessary */
              if (r<0) {
                if (pbc_dirs.x==1) {
                  r = cell_dim.x-1; 
                  pbc.x -= box_x.x;      
                  pbc.y -= box_x.y;
#ifndef TWOD
                  pbc.z -= box_x.z;
#endif
                } else continue;
              }
              if (s<0) {
                if (pbc_dirs.y==1) {
                  s = cell_dim.y-1;
                  pbc.x -= box_y.x;      
                  pbc.y -= box_y.y;
#ifndef TWOD
                  pbc.z -= box_y.z;
#endif
                } else continue;
              }
#ifndef TWOD
              if (t<0) {
                if (pbc_dirs.z==1) {
                  t = cell_dim.z-1;
                  pbc.x -= box_z.x;      
                  pbc.y -= box_z.y;
                  pbc.z -= box_z.z;
                } else continue;
              }
#endif
              if (r>cell_dim.x-1) {
                if (pbc_dirs.x==1) {
                  r = 0; 
                  pbc.x += box_x.x;      
                  pbc.y += box_x.y;
#ifndef TWOD
                  pbc.z += box_x.z;
#endif
                } else continue;
              }
              if (s>cell_dim.y-1) {
                if (pbc_dirs.y==1) {
                  s = 0; 
                  pbc.x += box_y.x;      
                  pbc.y += box_y.y;
#ifndef TWOD
                  pbc.z += box_y.z;
#endif
                } else continue;
              }
#ifndef TWOD
              if (t>cell_dim.z-1) {
                if (pbc_dirs.z==1) {
                  t = 0; 
                  pbc.x += box_z.x;      
                  pbc.y += box_z.y;
                  pbc.z += box_z.z;
                } else continue;
              }
#endif

	      /* Neighbour cell (note that p==q ist possible) */
#ifdef TWOD
              q = PTR_2D_V(cell_array,r,s,cell_dim);
#else
              q = PTR_3D_V(cell_array,r,s,t,cell_dim);
#endif

	      /* For each atom in the second cell */

	      for( v=0; v<q->n; ++v)
		{
		  /* Check whether there is enough memory for tmp variables */
		  if( number >= num) {
			++num;
			d    = (vektor *) realloc( d,   num * sizeof(vektor)); 
			du   = (vektor *) realloc( du,  num * sizeof(vektor));
			tmp  = (vektor *) realloc( tmp, num * sizeof(vektor)); 
		      }
		  
		  /* Calculate distance */
		  d[number].x = q->ort[v].x - q->dsp[v].x - p->ort[u].x + p->dsp[u].x + pbc.x;
		  d[number].y = q->ort[v].y - q->dsp[v].y - p->ort[u].y + p->dsp[u].y + pbc.y;
#ifndef TWOD
		  d[number].z = q->ort[v].z - q->dsp[v].z - p->ort[u].z + p->dsp[u].z + pbc.z;
#endif

		  radius = sqrt( (double)(SPROD(d[number],d[number])) );

		  if ( radius > 0.01 && radius < r_max ) 
		    {

		      /* Calculate differences of displacements */
		      du[number].x = q->dsp[v].x - p->dsp[u].x;
		      du[number].y = q->dsp[v].y - p->dsp[u].y;
#ifndef TWOD
		      du[number].z = q->dsp[v].z - p->dsp[u].z;
#endif

		      ++number; /* Number of neighbours*/

		    }
		} /* v */
	    } /* Neighbours of p */
	      
	      /* Calculate transformation matrix of du */

	      /* a = dx * dxT */

	      totalnumber += number;
	      maxnumber    = (number>maxnumber) ? number : maxnumber;
	      minnumber    = (number<minnumber) ? number : minnumber;

	      /* At least 3 neighbour atoms are required */
	      if (number < 3) { 
		p->empty[u] = 1;
		++emptynumber;
	      }
	      else
#ifdef TWOD
		{
	  for(g=0; g<2; ++g)
	    for(h=0; h<2; ++h)
	       a[g][h] = 0.0;

	  for (w=0; w<number; ++w) 
	    {
	      a[0][0] += d[w].x * d[w].x;  a[0][1] += d[w].x * d[w].y; 
	      a[1][0] += d[w].y * d[w].x;  a[1][1] += d[w].y * d[w].y;
	    }

	  /* b = Inverse of a */

	  det = a[0][0] * a[1][1] - a[0][1] * a[1][0];

	  if (det == 0.0) error("Transformation matrix zero.");

	  b[0][0] =  a[1][1] / det;
	  b[0][1] = -a[0][1] / det;
	  b[1][0] = -a[1][0] / det;
	  b[1][1] =  a[0][0] / det;

	  /* tmp = dx * b  */

	  for (w=0; w<number; ++w)
	    {
	      tmp[w].x = d[w].x * b[0][0] + d[w].y * b[1][0];
	      tmp[w].y = d[w].x * b[0][1] + d[w].y * b[1][1]; 
	    }

	  /* strain = (symmetrized) du * tmp  */

	  for (w=0; w<number; ++w)
	    {
	      p->strain[u].x        += du[w].x * tmp[w].x;
	      p->strain[u].y        += du[w].y * tmp[w].y;
	      
	      p->strain_offdia[u].x += ( du[w].y * tmp[w].x + du[w].x * tmp[w].y ) / 2;
	    }
#else /* 3D */
	{

	  for (g=0; g<3; ++g)
	    for (h=0; h<3; ++h)
	      a[g][h] = 0.0;

	  for (w=0; w<number; ++w) 
	    {
	      a[0][0] += d[w].x * d[w].x;  a[0][1] += d[w].x * d[w].y;  a[0][2] += d[w].x * d[w].z;
	      a[1][0] += d[w].y * d[w].x;  a[1][1] += d[w].y * d[w].y;  a[1][2] += d[w].y * d[w].z;
	      a[2][0] += d[w].z * d[w].x;  a[2][1] += d[w].z * d[w].y;  a[2][2] += d[w].z * d[w].z;
	    }

	  /* b = Inverse of a */

	  det = a[0][0] * a[1][1] * a[2][2] +
	        a[0][1] * a[1][2] * a[2][0] +
	        a[0][2] * a[1][0] * a[2][1] -
	        a[2][0] * a[1][1] * a[0][2] -
	        a[2][1] * a[1][2] * a[0][0] -
	        a[2][2] * a[1][0] * a[0][1];

	  if (det == 0.0) error("Transformation matrix singular.");

	  b[0][0] = ( a[1][1] * a[2][2] - a[1][2] * a[2][1] ) / det;
	  b[0][1] = ( a[2][1] * a[0][2] - a[0][1] * a[2][2] ) / det;
	  b[0][2] = ( a[0][1] * a[1][2] - a[1][1] * a[0][2] ) / det;

	  b[1][0] = ( a[1][2] * a[2][0] - a[1][0] * a[2][2] ) / det;
	  b[1][1] = ( a[0][0] * a[2][2] - a[2][0] * a[0][2] ) / det;
	  b[1][2] = ( a[0][2] * a[1][0] - a[0][0] * a[1][2] ) / det;
	  
	  b[2][0] = ( a[1][0] * a[2][1] - a[1][1] * a[2][0] ) / det;
	  b[2][1] = ( a[0][1] * a[2][0] - a[0][0] * a[2][1] ) / det;
	  b[2][2] = ( a[0][0] * a[1][1] - a[0][1] * a[1][0] ) / det;

	  /* tmp = dx * b  */

	  for (w=0; w<number; ++w)
	    {
	      tmp[w].x = d[w].x * b[0][0] + d[w].y * b[1][0] + d[w].z * b[2][0];
	      tmp[w].y = d[w].x * b[0][1] + d[w].y * b[1][1] + d[w].z * b[2][1];
	      tmp[w].z = d[w].x * b[0][2] + d[w].y * b[1][2] + d[w].z * b[2][2];
	    }

	  /* strain = (symmetrized) du * tmp  */

	  for (w=0; w<number; ++w)
	    {
	      p->strain[u].x += du[w].x * tmp[w].x;
	      p->strain[u].y += du[w].y * tmp[w].y;
	      p->strain[u].z += du[w].z * tmp[w].z;

	      p->strain_offdia[u].x += ( du[w].y * tmp[w].z + du[w].z * tmp[w].y ) / 2;
	      p->strain_offdia[u].y += ( du[w].x * tmp[w].z + du[w].z * tmp[w].x ) / 2;
	      p->strain_offdia[u].z += ( du[w].x * tmp[w].y + du[w].y * tmp[w].x ) / 2;
	    }
#endif	 /* Not TWOD */
	  
	} /* number > 2 */

		} /* u */

	    } /* First cell */

	  /* Statistics */
	  printf("Maximal number of neighbour atoms: %d\n", maxnumber);
	  printf("Minimal number of neighbour atoms: %d\n", minnumber);
	  printf("Average number of neighbour atoms: %f\n", (float) totalnumber/natoms );
	  if(emptynumber>0)
	    printf("Number of omitted atoms: %d (%.2f %%)\n", emptynumber, (float) emptynumber/natoms );

}









