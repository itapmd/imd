
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "param.h"

typedef struct {real    x; real    y; real    z; }  vektor;
typedef struct {integer x; integer y; integer z; } ivektor;

vektor  box_x, box_y, box_z;
vektor  *atoms, *atoms0, *atoms_old, *shift, *msqd, *sqd;
int     ntypes, natoms, *types, tt[4];

/* parameters with defaults */
int     msqd_int = 1;                      /* time interval for checkpoints */
char    outfiles[255] = "vasp";               /* base name for output files */
real    timestep = 0.002;               /* timestep of VASP simulation [ps] */
int     rav_num = 0;                     /* atom number for running average */
int     rav_int = 0;                        /* interval for running average */

/*****************************************************************************
*
* read parameter file
*
*****************************************************************************/

void readparamfile(char *paramfname)
{
  FILE *pf;
  char buffer[1024], *token, *res;

  /* open parameter file */
  pf = fopen(paramfname,"r");
  if (NULL == pf) error("Could not open parameter file");

  do {
    /* read one line */
    res=fgets(buffer,1024,pf);
    if (NULL == res) break;
    token = strtok(buffer," =\t\r\n");
    if (NULL == token) continue; /* skip blank lines */
    if (token[0]=='#') continue; /* skip comments */

    if (strcasecmp(token,"msqd_int")==0) {
      /* time interval for msqd writes */
      getparam(token,&msqd_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"timestep")==0) {
      /* base name of output files */
      getparam(token,&timestep,PARAM_REAL,1,1);
    }
    else if (strcasecmp(token,"rav_int")==0) {
      /* interval for running average */
      getparam(token,&rav_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"rav_num")==0) {
      /* atom number for running average */
      getparam(token,&rav_num,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"outfiles")==0) {
      /* base name of output files */
      getparam(token,outfiles,PARAM_STR,1,254);
    }
    else {
      /* warn if unknown parameter is used */
      fprintf(stderr, "****** Unknown TAG %s ignored ******", token);
    }
  } while (!feof(pf));
  fclose(pf);
}

/******************************************************************************
*
*  return error message and stop
*
******************************************************************************/

void error(char *msg)
{
  fprintf(stderr,"Error: %s\n",msg);
  exit(2);
}

/******************************************************************************
*
*  endian returns 1 if system is big endian, 0 if little endian
*
******************************************************************************/

int endian(void)
{
  unsigned short int word = 0x0001;
  unsigned char  *byte    = (unsigned char *) &word;
  return (byte[0] ? 0 : 1);
}

/******************************************************************************
*
*  read box an types from POSCAR
*
******************************************************************************/

void read_poscar()
{
  FILE *infile;
  char buf[255];
  real sc, x, y, z;
  int  i, j, k;

  /* open POSCAR file */
  infile = fopen("POSCAR","r");
  if (NULL==infile) error("Cannot open POSCAR file.");
  fgets(buf, 255, infile);

  /* get box vectors */
  fscanf(infile, "%lf", &sc);
  fscanf(infile, "%lf %lf %lf", &x, &y, &z);
  box_x.x = sc * x; 
  box_x.y = sc * y;
  box_x.z = sc * z; 
  fscanf(infile, "%lf %lf %lf", &x, &y, &z);
  box_y.x = sc * x;
  box_y.y = sc * y;
  box_y.z = sc * z; 
  fscanf(infile, "%lf %lf %lf", &x, &y, &z);
  box_z.x = sc * x; 
  box_z.y = sc * y; 
  box_z.z = sc * z; 
  printf("box_x: %f %f %f\n", box_x.x, box_x.y, box_x.z);
  printf("box_y: %f %f %f\n", box_y.x, box_y.y, box_y.z);
  printf("box_z: %f %f %f\n", box_z.x, box_z.y, box_z.z);

  /* get natoms and ntypes */
  ntypes = fscanf(infile, "%d %d %d %d", tt, tt+1, tt+2, tt+3);
  natoms = 0;
  for (i=0; i<ntypes; i++) natoms += tt[i];
  types = (int *) malloc( natoms * sizeof(int) );
  if (NULL==types) error("Cannot allocate types vector.");
  k = 0;
  for (i=0; i<ntypes; i++)
    for (j=0; j<tt[i]; j++) types[k++] = i;
  printf("ntypes: %d\n", ntypes);
  printf("natoms: %d\n", natoms);

  /* close file */
  fclose(infile);
}

/******************************************************************************
*
*  write msqd
*
******************************************************************************/

void write_msqd(int step, FILE *outfile)
{
  vektor *a, *b, *c, *d, *m;
  real x, y, z, xx, yy, zz;
  int  i, k;

  for (k=0; k<ntypes; k++) {
    (msqd+k)->x = 0.0;
    (msqd+k)->y = 0.0;
    (msqd+k)->z = 0.0;
  }

  for (i=0; i<natoms; ++i) {

    a = atoms + i;
    b = atoms_old + i;
    c = atoms0 + i;
    d = shift + i;
    m = msqd + types[i];

    /* apply PBC */
    if      ((a->x - b->x) < -0.5) d->x += 1.0;
    else if ((a->x - b->x) >  0.5) d->x -= 1.0;
    if      ((a->y - b->y) < -0.5) d->y += 1.0;
    else if ((a->y - b->y) >  0.5) d->y -= 1.0;
    if      ((a->z - b->z) < -0.5) d->z += 1.0;
    else if ((a->z - b->z) >  0.5) d->z -= 1.0;
    
    /* get atom displacement */
    x  = a->x + d->x - c->x;
    y  = a->y + d->y - c->y;
    z  = a->z + d->z - c->z;
    xx = x * box_x.x + y * box_y.x + z * box_z.x;
    yy = x * box_x.y + y * box_y.y + z * box_z.y;
    zz = x * box_x.z + y * box_y.z + z * box_z.z;
    m->x += xx*xx;
    m->y += yy*yy;
    m->z += zz*zz;

    /* copy new config to old config */
    b->x = a->x;
    b->y = a->y;
    b->z = a->z;
  }

  /* write msqd file */
  if (0 == step % msqd_int) {
    fprintf(outfile, "%e", timestep * step);
    for (k=0; k<ntypes; k++) fprintf(outfile, " %e %e %e", 
      (msqd+k)->x/tt[k], (msqd+k)->y/tt[k], (msqd+k)->z/tt[k]);
    fprintf(outfile, "\n");
  }

}

/******************************************************************************
*
*  write sqd
*
******************************************************************************/

void write_sqd()
{
  FILE *outfile;
  char fname[255];
  vektor *a, *c, *d;
  real x, y, z, xx, yy, zz;
  int  i;

  sprintf(fname, "%s.sqd", outfiles);
  outfile=fopen(fname, "w");
  if (NULL==outfile) error("Cannot open sqd file");

  for (i=0; i<natoms; ++i) {

    a = atoms + i;
    c = atoms0 + i;
    d = shift + i;

    /* get atom displacement */
    x  = a->x + d->x - c->x;
    y  = a->y + d->y - c->y;
    z  = a->z + d->z - c->z;
    xx = x * box_x.x + y * box_y.x + z * box_z.x;
    yy = x * box_x.y + y * box_y.y + z * box_z.y;
    zz = x * box_x.z + y * box_y.z + z * box_z.z;
    fprintf(outfile, "%d %e %e %e\n", types[i], xx*xx, yy*yy, zz*zz);
  }
  fclose(outfile);
}

/******************************************************************************
*
*  read configs from XDATCAR
*
******************************************************************************/

void read_xdatcar()
{
  FILE *infile, *msqdfile;
  char buf[255], fname[255];
  int i, cnt, cnt2;
  vektor *a, *rav, norm, sum = {0.0, 0.0, 0.0};

  /* open POSCAR file */
  infile = fopen("XDATCAR","r");
  if (NULL==infile) error("Cannot open XDATCAR file.");
  fgets(buf, 255, infile);
  fgets(buf, 255, infile);
  fgets(buf, 255, infile);
  fgets(buf, 255, infile);
  fgets(buf, 255, infile);
  fgets(buf, 255, infile);
  atoms     = (vektor *) malloc( natoms * sizeof(vektor) );
  atoms0    = (vektor *) malloc( natoms * sizeof(vektor) );
  atoms_old = (vektor *) malloc( natoms * sizeof(vektor) );
  shift     = (vektor *) calloc( natoms,  sizeof(vektor) );
  msqd      = (vektor *) calloc( ntypes,  sizeof(vektor) );  
  if ((NULL==atoms) || (NULL==atoms0) || (NULL==atoms_old) || 
      (NULL==shift) || (msqd==NULL)) error("Cannot allocate atoms memory."); 
  if (rav_int) {
    rav = (vektor *) malloc( rav_int * sizeof(vektor) );
    if (NULL==rav) error("Cannot allocate memory.");
  }
  if (rav_int) {
    rav = (vektor *) malloc( rav_int * sizeof(vektor) );
    if (NULL==rav) error("Cannot allocate memory.");
  }

  /* read first config */
  for (i=0; i<natoms; i++) {
    fgets(buf, 255, infile);
    a = atoms0 + i;
    if (3!=sscanf(buf, "%lf %lf %lf\n", &(a->x), &(a->y), &(a->z)))
      error("Atom with less than three components!");
    (atoms_old+i)->x = a->x;
    (atoms_old+i)->y = a->y;
    (atoms_old+i)->z = a->z;
  } 
  fgets(buf, 255, infile);

  /* open msqd file */
  sprintf(fname, "%s.msqd", outfiles);
  msqdfile = fopen(fname,"w");
  if (NULL==msqdfile) error("Cannot open msqd file.");

  /* process configs */
  cnt = 0;
  while (!feof(infile)) {
    for (i=0; i<natoms; i++) {
      fgets(buf, 255, infile);
      a = atoms + i;
      if (3!=sscanf(buf, "%lf %lf %lf\n", &(a->x), &(a->y), &(a->z)))
	error("Atom with less than three components!");
    } 
    /* subtract running average */
    if (rav_int) {
      if (cnt < rav_int) {
        rav[cnt].x = atoms[rav_num].x;
	rav[cnt].y = atoms[rav_num].y;
	rav[cnt].z = atoms[rav_num].z;
        sum.x += rav[cnt].x;
        sum.y += rav[cnt].y;
        sum.z += rav[cnt].z;
      } else {
        if (cnt==rav_int) norm = sum;
        cnt2 = cnt % rav_int;
        sum.x -= rav[cnt2].x;
        sum.y -= rav[cnt2].y;
        sum.z -= rav[cnt2].z;
        rav[cnt2].x = atoms[rav_num].x;
        rav[cnt2].y = atoms[rav_num].y;
        rav[cnt2].z = atoms[rav_num].z;
        sum.x += rav[cnt2].x;
        sum.y += rav[cnt2].y;
        sum.z += rav[cnt2].z;
        for (i=0; i<natoms; i++) {
          atoms[i].x -= (sum.x - norm.x) / rav_int;
          atoms[i].y -= (sum.y - norm.y) / rav_int;
          atoms[i].z -= (sum.z - norm.z) / rav_int;
	}
      }
    }
    write_msqd(cnt,msqdfile);
    cnt++;
    fgets(buf, 255, infile);
  }
  fclose(infile);
  fclose(msqdfile);
  printf("Number of configs: %d\n", cnt+1);
}

/******************************************************************************
*
*  main
*
******************************************************************************/

int main(int argc, char **argv) 
{
  readparamfile(argv[1]);
  read_poscar();
  read_xdatcar();
  write_sqd();
  return 0;
}
