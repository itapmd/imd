
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

typedef double real;
typedef struct { real x; real y; real z; } vektor;

typedef struct {
  int    typ;
  real   r2;
  vektor dist;
} neigh_t;

typedef struct {
  int    typ;
  int    n_neigh;
  vektor pos;
  neigh_t *neigh;
} atom_t;

typedef struct {
  real *begin;      /* first value in the table */
  real *end;        /* last value in the table */
  real *step;       /* table increment */
  real *invstep;    /* inverse of increment */
  int  *first;      /* index of first entry */
  int  *last;       /* index of last entry */
  int  len;         /* total length of the table */
  int  ncols;       /* number of columns */
  real *table;      /* the actual data */
} pot_table_t;


/* the global variables */
int    ntypes=1;       /* number of atom types - where from? */
int    natoms=0;       /* number of atoms - from read_config */
atom_t *atoms=NULL;    /* atoms array - from read_config */
real   *force_0=NULL;  /* the forces we aim at - from read_config */
pot_table_t pair_pot;  /* the potential table, from read_pot_table */


/******************************************************************************
*
* error -- Complain and abort
*
******************************************************************************/

void error(char *msg)
{
  fprintf(stderr,"Error: %s\n",msg);
  fflush(stderr);
  exit(2);
}


/******************************************************************************
*
* read potential table
*
******************************************************************************/

void read_pot_table( pot_table_t *pt, char *filename, int ncols )
{
  FILE *infile;
  char buffer[1024], msg[255], *res;
  int  have_format=0, end_header;
  int  format, size, i, j, *nvals;
  real *val;

  /* open file */
  infile = fopen(filename,"r");
  if (NULL == infile) {
    sprintf(msg,"Could not open file %s\n",filename);
    error(msg);
  }

  /* read the header */
  do {
    /* read one line */
    res=fgets(buffer,1024,infile);
    if (NULL == res) {
      sprintf(msg,"Unexpected end of file in %s",filename);
      error(msg);
    }
    /* check if it is a header line */
    if (buffer[0]!='#') {
      sprintf(msg,"Header corrupt in file %s",filename);
      error(msg);
    }
    /* stop after last header line */
    if (buffer[1]=='E') {
      end_header = 1;
    }
    /* see if it is the format line */
    else if (buffer[1]=='F') {
      /* format complete? */
      if (2!=sscanf( (const char*)(buffer+2), "%d %d", &format, &size )) {
        sprintf(msg,"Corrupt format header line in file %s",filename);
        error(msg);
      }
      /* right number of columns? */
      if (size!=ncols) {
        sprintf(msg,"Wrong number of data columns in file %s",filename);
        error(msg);
      }
      /* recognized format? */
      if (format!=3) {
        sprintf(msg,"Unrecognized format specified for file %s",filename);
        error(msg);
      }
      have_format=1;
    } 
    else { 
      /* header does not end properly */
      sprintf(msg,"Corrupted header in file %s",filename);
      error(msg);
    }
  } while (!end_header);

  /* did we have a format in the header? */
  if (!have_format) {
    sprintf(msg,"Format not specified in header of file %s",filename);
    error(msg);
  }

  /* allocate info block of function table */
  pt->len     = 0;
  pt->ncols   = size;
  pt->begin   = (real *) malloc(size*sizeof(real));
  pt->end     = (real *) malloc(size*sizeof(real));
  pt->step    = (real *) malloc(size*sizeof(real));
  pt->invstep = (real *) malloc(size*sizeof(real));
  pt->first   = (int  *) malloc(size*sizeof(int ));
  pt->last    = (int  *) malloc(size*sizeof(int ));
  nvals       = (int  *) malloc(size*sizeof(int ));
  if ((pt->begin   == NULL) || (pt->end   == NULL) || (pt->step == NULL) || 
      (pt->invstep == NULL) || (pt->first == NULL) || (pt->last == NULL) || 
      (nvals       == NULL)) {
    sprintf(msg,"Cannot allocate info block for potential table %s",filename);
    error(msg);
  }

  /* read the info block of the function table */
  for(i=0; i<ncols; i++) {
    if (3>fscanf(infile,"%lf %lf %d", &pt->begin[i], &pt->end[i], &nvals[i])) {
        sprintf(msg, "Premature end of potential file %s", filename);
        error(msg);
    }
    pt->step[i] = (pt->end[i] - pt->begin[i]) / (nvals[i]-1);
    pt->invstep[i] = 1.0 / pt->step[i];
    if (i==0) pt->first[i] = 0; else pt->first[i] = pt->last[i-1] + 1;
    pt->last[i] = pt->first[i] + nvals[i] - 1;
    pt->len = pt->first[i] + nvals[i];
  }

  /* allocate the function table */
  pt->table = (real *) malloc(pt->len * sizeof(real));
  if (NULL==pt->table) {
    error("Cannot allocate memory for potential table");
  }

  /* input loop */
  val = pt->table;
  for (i=0; i<ncols; i++) {
    for (j=0; j<nvals[i]; j++) {
      if (1>fscanf(infile, "%lf", val)) {
        sprintf(msg, "Premature end of potential file %s", filename);
        error(msg);
      } else val++;
    }
  }

  fclose(infile);
}


/*****************************************************************************
*
*  Evaluate derivative of potential with quadratic interpolation. 
*  We need (1/r)(dV/dr) = 2 * dV/dr^2, so we use table with equidistant r^2. 
*  col is typ1 * ntypes + typ2.
*
******************************************************************************/

real grad2(pot_table_t *pt, int col, real r2)
{
  real r2a, istep, chi, p0, p1, p2, dv, d2v;
  int  k;

  /* check for distances shorter than minimal distance in table */
  r2a = r2 - pt->begin[col];
  if (r2a < 0) {
    r2a   = 0;
  }

  /* indices into potential table */
  istep = pt->invstep[col];
  k     = (int) (r2a * istep);
  chi   = (r2a - k * pt->step[col]) * istep;
  k    += pt->first[col];

  /* intermediate values */
  p0  = (k<=pt->last[col]) ? pt->table[k++] : 0.0;
  p1  = (k<=pt->last[col]) ? pt->table[k++] : 0.0;
  p2  = (k<=pt->last[col]) ? pt->table[k++] : 0.0;
  dv  = p1 - p0;
  d2v = p2 - 2 * p1 + p0;

  /* return twice the derivative */
  return 2 * istep * (dv + (chi - 0.5) * d2v);
}


/*****************************************************************************
*
*  read the configurations
*
******************************************************************************/

void read_config(char *filename)
{
  FILE    *infile;
  char    msg[255];
  int     nconf=0, count, i, j, k;
  atom_t  *atom;
  neigh_t *neigh;

  /* open file */
  infile = fopen(filename,"r");
  if (NULL == infile) {
    sprintf(msg,"Could not open file %s\n",filename);
    error(msg);
  }

  /* read configurations until the end of the file */
  do {

    /* number of atoms in this configuration */
    if (1>fscanf(infile, "%d", &count)) error("Unexpected end of file");

    /* increase memory for this many additional atoms */
    atoms = (atom_t *) realloc(atoms, (natoms+count) * sizeof(atom_t));
    if (NULL==atoms)   error("Cannot allocate memory for atoms");
    force_0 = (real *) realloc(force_0, 3 * (natoms+count) * sizeof(real));
    if (NULL==force_0) error("Cannot allocate memory for forces");

    /* read the neighbors */
    for (i=0; i<count; i++) {

      k    = 3 * (natoms + i);
      atom = atoms + natoms + i;
      if (8>fscanf(infile,"%d %lf %lf %lf %lf %lf %lf %d\n", &(atom->typ), 
                   &(atom->pos.x), &(atom->pos.y), &(atom->pos.z), 
                   force_0+k, force_0+k+1, force_0+k+2, &(atom->n_neigh)))
        error("Corrupt configuration file");

      /* allocate memory for neighbors */
      atom->neigh = (neigh_t *) malloc(atom->n_neigh * sizeof(neigh_t));
      if (NULL==atom->neigh) error("Cannot allocate memory for neighbors");

      /* read the neighbors */
      for (j=0; j<atom->n_neigh; j++) {
	neigh = atom->neigh + j;
        if (5>fscanf(infile,"%d %lf %lf %lf %lf\n", 
                     &(neigh->typ), &(neigh->r2), 
                     &(neigh->dist.x), &(neigh->dist.y), &(neigh->dist.z)))
          error("Corrupt configuration file");
      }
    }

    /* increment natoms and configuration number */
    natoms += count;
    nconf++;

  } while (!feof(infile));

  /* print diagnostic message and close file */
  printf("Read %d configurations with a total of %d atoms\n", nconf, natoms);
  fclose(infile);
}


/*****************************************************************************
*
*  compute forces using pair potentials 
*
******************************************************************************/

void calc_forces_pair(real *forces)
{
  int     i, j, k, typ1, typ2, col;
  atom_t  *atom;
  neigh_t *neigh;
  real    grad;

  for (i=0; i<natoms; i++) {

    atom = atoms + i;
    typ1 = atom->typ;
    k    = 3*i;

    for (j=0; j<atom->n_neigh; j++) {

      neigh = atom->neigh+j;
      typ2  = neigh->typ;
      col   = (typ1 <= typ2) ? typ1 * ntypes + typ2 : typ2 * ntypes + typ1;

      if (neigh->r2 < pair_pot.end[col] + pair_pot.step[col]) {
        grad = grad2( &pair_pot, col, neigh->r2);
        forces[k  ] += neigh->dist.x * grad;
        forces[k+1] += neigh->dist.y * grad;
        forces[k+2] += neigh->dist.z * grad;
      }
    }
  }
}


int main(int argc, char **argv)
{
  read_config(argv[1]);
}
