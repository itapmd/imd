
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
int     ntypes, natoms, *types;

/* parameters with defaults */
vektor  pic_ll  = { -4.8, -8.4, 0.00 };               /* lower left  corner */
vektor  pic_ur  = { 11.2,  0.8, 4.08 };               /* upper right corner */
int     pic_int = 1;                       /* time interval for checkpoints */
real    phi     = 0.0;       /* rotation angle -- gets multiplied with 2*Pi */
char    outfiles[255] = "vasp";               /* base name for output files */
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
  int  ii[2];

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

    if (strcasecmp(token,"pic_ll")==0) {
      /* lower left corner of box */
      getparam(token,&pic_ll.x,PARAM_REAL,3,3);
    }
    else if (strcasecmp(token,"pic_ur")==0) {
      /* upper right corner of box */
      getparam(token,&pic_ur.x,PARAM_REAL,3,3);
    }
    else if (strcasecmp(token,"pic_int")==0) {
      /* time interval for checkpoints */
      getparam(token,&pic_int,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"spec_typ")==0) {
      /* special type for certain atoms */
      getparam(token,ii,PARAM_INT,2,2);
      types[ii[0]] = ii[1];      
    }
    else if (strcasecmp(token,"phi")==0) {
      /* rotation angle */
      getparam(token,&phi,PARAM_REAL,1,1);
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
  int tt[4], i, j, k;

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
*  rotate box around z-axis
*
******************************************************************************/

void rotate_box()
{
  real x, y;

  phi *= 8 * atan(1.0);
  printf("phi: %f\n", phi);

  x =  cos(phi) * box_x.x + sin(phi) * box_x.y; 
  y = -sin(phi) * box_x.x + cos(phi) * box_x.y;
  box_x.x = x;
  box_x.y = y;

  x =  cos(phi) * box_y.x + sin(phi) * box_y.y; 
  y = -sin(phi) * box_y.x + cos(phi) * box_y.y;
  box_y.x = x;
  box_y.y = y;

  printf("box_x: %f %f %f\n", box_x.x, box_x.y, box_x.z);
  printf("box_y: %f %f %f\n", box_y.x, box_y.y, box_y.z);
  printf("box_z: %f %f %f\n", box_z.x, box_z.y, box_z.z);
}

/******************************************************************************
*
*  write IMD checkpoint file
*
******************************************************************************/

void write_chkpt(int num, vektor *atoms)
{
  FILE *outfile;
  char fname[255];
  real x, y, z, xx, yy, zz;
  int  k, ix, iy, iz;

  /* open checkpoint file */
  sprintf(fname, "%s.%05d.cpt", outfiles, num);
  outfile = fopen(fname, "w");
  if (NULL==outfile) error("Cannot open checkpoint file.");

  /* write header */
  fprintf(outfile, "#F A 0 1 0 3 0 0\n");
  fprintf(outfile, "#C type x y z\n");
  fprintf(outfile, "#X %e 0.0 0.0\n", pic_ur.x - pic_ll.x);
  fprintf(outfile, "#Y 0.0 %e 0.0\n", pic_ur.y - pic_ll.y);
  fprintf(outfile, "#Z 0.0 0.0 %e\n", pic_ur.z - pic_ll.z);
  fprintf(outfile, "#E\n");

  /* write atoms */
  for (k=0; k<natoms; ++k) 
    for (ix=-2; ix<=2; ix++)
      for (iy=-2; iy<=2; iy++)
        for (iz=-2; iz<=2; iz++) {

          /* get atom coordinates */
	  x  = (atoms+k)->x + ix;
	  y  = (atoms+k)->y + iy;
          z  = (atoms+k)->z + iz;
          xx = x * box_x.x + y * box_y.x + z * box_z.x;
          yy = x * box_x.y + y * box_y.y + z * box_z.y;
          zz = x * box_x.z + y * box_y.z + z * box_z.z;

          /* continue if atom is not inside selected box */
          if ((xx < pic_ll.x) || (xx > pic_ur.x) ||
              (yy < pic_ll.y) || (yy > pic_ur.y) ||
	      (zz < pic_ll.z) || (zz > pic_ur.z)) continue; 

          /* write atom */
          fprintf(outfile, "%d %e %e %e\n", 
	    types[k], xx - pic_ll.x, yy - pic_ll.y, zz - pic_ll.z );
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
  FILE *infile;
  char buf[255];
  vektor *a, *atoms, *rav, norm, sum = {0.0, 0.0, 0.0};
  int i, cnt, cnt2;

  /* open POSCAR file */
  infile = fopen("XDATCAR","r");
  if (NULL==infile) error("Cannot open XDATCAR file.");
  fgets(buf, 255, infile);
  fgets(buf, 255, infile);
  fgets(buf, 255, infile);
  fgets(buf, 255, infile);
  fgets(buf, 255, infile);
  fgets(buf, 255, infile);
  atoms = (vektor *) malloc( natoms * sizeof(vektor) );
  if (NULL==atoms) error("Cannot allocate atoms list."); 
  if (rav_int) {
    rav = (vektor *) malloc( rav_int * sizeof(vektor) );
    if (NULL==rav) error("Cannot allocate memory.");
  }
  
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
    if (0==cnt%pic_int) write_chkpt(cnt/pic_int, atoms);
    cnt++;
    fgets(buf, 255, infile);
  }
  printf("Number of configs: %d\n", cnt);
}

/******************************************************************************
*
*  main
*
******************************************************************************/

int main(int argc, char **argv) 
{
  read_poscar();
  readparamfile(argv[1]);
  rotate_box();
  read_xdatcar();
  return 0;
}
