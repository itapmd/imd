/**********************************************************/
/*   procdata.c                                           */
/*   Weiterbearbeitung der Daten aus IMD und VVMD         */
/*   [ Konvertieren ins MOLPIC / RASMOL-PDB - Format,     */
/*   Berechnen der Radialverteilungsfunktion g(r) ]       */
/*   process data from IMD molecular dynamics program	  */
/*   perform various analysis procedures                  */
/*   1996 by Martin Hohl                                  */
/*   date of last modification: 30-Jun-97                 */
/**********************************************************/

/* format written by IMD:				*/
/* NUMBER SORT M X Y Z VX VY VZ				*/
/* NUMBER: number of atom		1	INT	*/
/* SORT: sort of atom			1	INT	*/
/* M: mass of atom			1	REAL	*/
/* X,Y,Z: position of atom		3	REAL	*/
/* VX,VY,VZ: velocity of atom		3	REAL	*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef AMIGA
#ifdef __SASC_60
#ifdef _IEEE
#include <mieeedoub.h>
#endif
#ifdef _M68881
#include <m68881.h>
#endif
#endif
#endif

#if defined(__i386__) && defined(linux) && defined(__i387)
#include <i386/__math.h>
#endif

#ifndef PI
#define PI 3.14159265
#endif

#ifdef USE_GETUTIME
/* Items from getutime.c */
#ifdef __STDC__
extern unsigned long getutime(void);
#else
extern unsigned long getutime();
#endif
extern unsigned long Ticks;
#endif

#define ALLOC_SIZE 1000

unsigned asize;

#define HISTSIZE 1000
#define NUMHIST 256

/* maximale Anzahl Atomsorten */
#define MAX_SORT 20

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE !FALSE
#endif

typedef double REAL;

struct ClipRange {
  REAL min_x, max_x;
  REAL min_y, max_y;
  REAL min_z, max_z;
};

struct RotateAngles {
  REAL phi_x,phi_y,phi_z;
};

typedef REAL VEC[3];
typedef REAL MAT[3][3];

void Usage(char *s)
{
  fprintf(stderr,"Usage: %s [-i{format}] [-f{format}] inputfile outputfile\n",s);
  fprintf(stderr,"-i0 reads IMD or VVMD data format (default)\n");
  fprintf(stderr,"-i1 reads modifti.f/tiunmodif.f output format\n");
  fprintf(stderr,"-f0 converts to MOLPIC format (for molgraph, default)\n");
  fprintf(stderr,"-f1 converts to PDB format (for rasmol)\n");
  fprintf(stderr,"-f2 converts to IMD format\n");
  fprintf(stderr,"-k keeps previous numbers (default: sort and begin from 1\n");
  fprintf(stderr,"-fRDF calculates radial distribution function\n");
  fprintf(stderr,"-Cxmin limit\tcuts off atoms with x<limit\n");
  fprintf(stderr,"-Cxmax limit\tcuts off atoms with x>limit\n");
  fprintf(stderr,"-Cymin limit\tcuts off atoms with y<limit\n");
  fprintf(stderr,"-Cymax limit\tcuts off atoms with y>limit\n");
  fprintf(stderr,"-Czmin limit\tcuts off atoms with z<limit\n");
  fprintf(stderr,"-Czmax limit\tcuts off atoms with z>limit\n");
  fprintf(stderr,"-Sz zincr n\tstacks the structure in z-direction n times with increment zincr\n");
  fprintf(stderr,"\t\t(used for generation of 2D-quasiperiodic structures)\n");
  fprintf(stderr,"-Rx phi   rotate by angle phi around x-axis\n");
  fprintf(stderr,"-Ry phi   rotate by angle phi around y-axis\n");
  fprintf(stderr,"-Rz phi   rotate by angle phi around z-axis\n");
  fprintf(stderr,"very special options are: (see source code)\n");
  fprintf(stderr,"-RAlCuCo -CAlCuCo  rotates AlCuCo approximant\n");
  fprintf(stderr,"and cuts out a rectangular piece (for approx2.inp)\n");
  fflush(stderr);
}

typedef struct {
  double x,y,z;
} Coord3D;

struct AtomDaten {
  double vx,vy,vz; /* velocity */
  double m; /* mass */
  unsigned short sorte; /* sort of atom */
  unsigned short nr; /* number of atom */
};

Coord3D *atome;
struct AtomDaten *atomdaten;

double *radius;

double x_min = 1E30,x_max = 0; /* limits for box, x-coordinate */
double y_min = 1E30,y_max = 0; /* limits for box, y-coordinate */
double z_min = 1E30,z_max = 0; /* limits for box, z-coordinate */

double x_width = 0, x_width2 = 0, xw_inv = 0;
double y_width = 0, y_width2 = 0, yw_inv = 0;
double z_width = 0, z_width2 = 0, zw_inv = 0;

unsigned minnr = 0,maxnr = 0;

int keep;

/* read in atom data from tiunmodif.f - JR type format */
unsigned lies_atome_jr(FILE *fp, Coord3D **atome, struct AtomDaten **atomdaten)
{
  unsigned nr;
  double wx,wy,wz; /* box width */
  double x,y,z; /* position */
  unsigned short sorte; /* atom sort */
  int retcode;
  Coord3D *node;
  struct AtomDaten *nodedata;
  char line[256];

  asize = ALLOC_SIZE;
  nr = 0;
  /* Zuerst einmal die Box-Dimensionen einlesen */
  fscanf(fp,"%lg %lg %lg\n",&wx,&wy,&wz);

  while (!feof(fp)) {
    line[0]='\0';
    fgets(line,256,fp);
    if (line[0]=='#') continue;
    if (line[0]=='\n') continue;
    if (line[0]=='\0') continue;
    retcode = sscanf(line,"%lg %lg %lg %hd",&x,&y,&z,&sorte);
    if ((retcode != 0) && (retcode != EOF)) {

      if (nr >= asize) { /* need more storage */
        asize += ALLOC_SIZE;
        *atome = realloc(*atome,asize*sizeof(Coord3D));
        *atomdaten = realloc(*atomdaten,asize*sizeof(struct AtomDaten));
        if ((*atome == NULL) || (*atomdaten == NULL)) {
          perror("realloc_atome");
          exit(10);
        }
        memset(&(*atome)[asize-ALLOC_SIZE],0,sizeof(Coord3D)*ALLOC_SIZE);
        memset(&(*atomdaten)[asize-ALLOC_SIZE],0,sizeof(struct AtomDaten)*ALLOC_SIZE);
      }
      node = &((*atome)[nr]);

      node->x=x;
      node->y=y;
      node->z=z;
      nodedata = &((*atomdaten)[nr]);
      nodedata->nr = nr;
      nodedata->sorte=sorte-1;
      nodedata->m=1.0;
      nodedata->vx=0.0;
      nodedata->vy=0.0;
      nodedata->vz=0.0;
      if (x_min>x) x_min=x;
      if (x_max<x) x_max=x;
      if (y_min>y) y_min=y;
      if (y_max<y) y_max=y;
      if (z_min>z) z_min=z;
      if (z_max<z) z_max=z;

      nr++;
    }
    else break;
  }
  printf("x_min=%g y_min=%g z_min=%g\n",x_min,y_min,z_min);
  printf("x_max=%g y_max=%g z_max=%g\n",x_max,y_max,z_max);

  return nr; /* number of atoms read in */
} /* lies_atome_jr */

/* read in atom data from IMD or VVMD */
unsigned lies_atome(FILE *fp, Coord3D **atome, struct AtomDaten **atomdaten)
{
  unsigned nread,nr;
  double x,y,z; /* position */
  double vx,vy,vz; /* velocity */
  double m; /* mass */
  unsigned short sorte; /* atom sort */
  int retcode;
  Coord3D *node;
  struct AtomDaten *nodedata;
  char line[256];

  asize = ALLOC_SIZE;
  nread = 0;
  minnr = 65535;
  maxnr = 0;
  nread = 0;
  while (!feof(fp)) {
    line[0]='\0';
    fgets(line,256,fp);
    if (line[0]=='#') continue;
    if (line[0]=='\n') continue;
    if (line[0]=='\0') continue;
    vx = vy = vz = 0;
    retcode = sscanf(line,"%u %hu %lg %lg %lg %lg %lg %lg %lg",&nr,&sorte,&m,&x,&y,&z,&vx,&vy,&vz);
    if ((retcode != 0) && (retcode != EOF)) {
      if (nr<minnr) minnr=nr;
      if (nr>maxnr) maxnr=nr;

      if (nr >= asize) { /* need more storage */
        unsigned osize = asize;
        asize = ((nr+ALLOC_SIZE)/ALLOC_SIZE)*ALLOC_SIZE;
        if (asize > osize) {
/*          printf("osize=%u asize=%u\n",osize,asize); */
          *atome = realloc(*atome,asize*sizeof(Coord3D));
          *atomdaten = realloc(*atomdaten,asize*sizeof(struct AtomDaten));
          if ((*atome == NULL) || (*atomdaten == NULL)) {
            perror("realloc_atome");
            exit(10);
          }
          memset(&(*atome)[osize],0,sizeof(Coord3D)*(asize-osize));
          memset(&(*atomdaten)[osize],0,sizeof(struct AtomDaten)*(asize-osize));
        }
      }
      node = &((*atome)[nr]);

      node->x=x;
      node->y=y;
      node->z=z;
      nodedata = &((*atomdaten)[nr]);
      nodedata->nr=nr;
      nodedata->sorte=sorte;
      nodedata->m=m;
      nodedata->vx=vx;
      nodedata->vy=vy;
      nodedata->vz=vz;
      if (x_min>x) x_min=x;
      if (x_max<x) x_max=x;
      if (y_min>y) y_min=y;
      if (y_max<y) y_max=y;
      if (z_min>z) z_min=z;
      if (z_max<z) z_max=z;

      nread++;
/*      node++; */
    }
    else break;
  }
  printf("x_min=%g y_min=%g z_min=%g\n",x_min,y_min,z_min);
  printf("x_max=%g y_max=%g z_max=%g\n",x_max,y_max,z_max);

  return nread; /* number of atoms read in */
} /* lies_atome */

/* Schreiben einer MOLPIC-Datei fuer das molgraph-Visualisierungsprogramm */
void writemgdat(char *out_name, Coord3D *atome, struct AtomDaten *atomdata, unsigned num_atoms)
{
  unsigned i;
  Coord3D *node;
  struct AtomDaten *nodedata;
  double sx_min,sx_max, sy_min,sy_max, sz_min,sz_max;
  double x_offset,y_offset,z_offset;
  FILE *out;

  node = &atome[minnr];
  nodedata = &atomdaten[minnr];
  out = fopen(out_name,"w");
  if (out == NULL) {
    perror("WriteMGdat");
    exit(10);
  };

  /* Eckpunkte der Box sollten fuer MOLGRAPH symmetrisch bzgl. Nullpunkt sein */
  sx_min = -(floor(x_min)+ceil(x_max))/2.0;
  sx_max = +(floor(x_min)+ceil(x_max))/2.0;
  x_offset = floor(x_min)-sx_min;
  sy_min = -(floor(y_min)+ceil(y_max))/2.0;
  sy_max = +(floor(y_min)+ceil(y_max))/2.0;
  y_offset = floor(y_min)-sy_min;
  sz_min = -(floor(z_min)+ceil(z_max))/2.0;
  sz_max = +(floor(z_min)+ceil(z_max))/2.0;
  z_offset = floor(z_min)-sz_min;

  fputs("GEOMETRY BOX\n",out);
  fprintf(out,"RANGE X %g %g\n",sx_min,sx_max);
  fprintf(out,"RANGE Y %g %g\n",sy_min,sy_max);
  fprintf(out,"RANGE Z %g %g\n",sz_min,sz_max);

  fprintf(out,"ATOMS %u\n",num_atoms);

  for (i=num_atoms; i!=0; i--) {
    fprintf(out,"%g %g %g %g %u\n",node->x-x_offset,node->y-y_offset,node->z-z_offset,radius[nodedata->sorte],nodedata->sorte+1);
    node++;
    nodedata++;
  }

  fputs("Box-Geometrie:\n",stderr);
  fprintf(stderr,"X-Bereich %g .. %g\n",sx_min,sx_max);
  fprintf(stderr,"Y-Bereich %g .. %g\n",sy_min,sy_max);
  fprintf(stderr,"Z-Bereich %g .. %g\n",sz_min,sz_max);
  fprintf(stderr,"Die Datei %s wurde erstellt.\n",out_name);
  fclose(out);
} /* writemgdat */

#if defined(__GNUC__) || defined(__SASC_60)
__inline static double square(double x)
{
  return x*x;
}
#else
#define square(x) ((x)*(x))
#endif

/*
char *names[] = {"AL","CO","PD","MN","CU","FE","NI","CR","V","TI","SI"};
int ordnum[11] = {13,27,46,25,29,26,28,24,23,22,14};
*/


/* Schreiben einer PDB-Datei */
void writepdb(char *out_name, Coord3D *atome, struct AtomDaten *atomdaten, unsigned num_atoms)
{
  unsigned i;
  Coord3D *node;
  struct AtomDaten *nodedata;
  char *names[2]={"AL","CO"};
  int ordnum[2]={13,27};
  FILE *out;

  out = fopen(out_name,"w");
  if (out == NULL) {
    perror("WritePDB");
    exit(10);
  };

  node = &atome[minnr];
  nodedata = &atomdaten[minnr];

  if (keep) {
    for (i=0; i<num_atoms; i++) {
      fprintf(out,"ATOM  %5d  %-3s      %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",nodedata->nr,
        names[nodedata->sorte],ordnum[nodedata->sorte],node->x,node->y,node->z,(double)1.0,nodedata->m*(square(nodedata->vx)+square(nodedata->vy)+square(nodedata->vz)));
      node++;
      nodedata++;
    }
  }
  else {
    for (i=0; i<num_atoms; i++) {
      fprintf(out,"ATOM  %5d  %-3s      %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",i+1,
        names[nodedata->sorte],ordnum[nodedata->sorte],node->x,node->y,node->z,(double)1.0,nodedata->m*(square(nodedata->vx)+square(nodedata->vy)+square(nodedata->vz)));
      node++;
      nodedata++;
    }
  }
  fclose(out);
} /* writepdb */


/* Schreiben einer CFG-Datei fuer VVMD oder IMD */
void writecfg(char *out_name, Coord3D *atome, struct AtomDaten *atomdaten, unsigned num_atoms)
{
  unsigned i;
  Coord3D *node;
  struct AtomDaten *nodedata;
  FILE *out;
  double x_offset,y_offset,z_offset;

  node = &atome[minnr];
  /* Offsets auf einen definierten Wert setzen */
  x_offset = node->x;
  y_offset = node->y;
  z_offset = node->z;
  for (i=0; i<num_atoms; i++) {
    if (node->x < x_offset) x_offset = node->x;
    if (node->y < y_offset) y_offset = node->y;
    if (node->z < z_offset) z_offset = node->z;
    node++;
  }

  out = fopen(out_name,"w");
  if (out == NULL) {
    perror("WriteCFG");
    exit(10);
  };

  node = &atome[minnr];
  nodedata = &atomdaten[minnr];
  if (keep) {
    for (i=0; i<num_atoms; i++) {
      fprintf(out,"%lu %u %g %g %g %g %g %g %g\n",(unsigned long)nodedata->nr,nodedata->sorte,
        nodedata->m,node->x-x_offset,node->y-y_offset,node->z-z_offset,
        nodedata->vx,nodedata->vy,nodedata->vz);
      node++;
      nodedata++;
    }
  }
  else {
    for (i=0; i<num_atoms; i++) {
      fprintf(out,"%lu %u %g %g %g %g %g %g %g\n",(unsigned long)i,nodedata->sorte,
        nodedata->m,node->x-x_offset,node->y-y_offset,node->z-z_offset,
        nodedata->vx,nodedata->vy,nodedata->vz);
      node++;
      nodedata++;
    }
  }
  fclose(out);
} /* writecfg */

/* Berechnung der (totalen) Radialverteilung */
void radialdist(char *out_name, Coord3D *atom, struct AtomDaten *atomdaten, int N)
{
  register double rijsq;
  double celldiag, V; /* Diagonale, Volumen */
  double dr;
  double xwidth,ywidth,zwidth;
#ifdef VECTORIZE
  long *radhist[NUMHIST];
  double *delta;
  long sum;
  int n;
#else
  long *radhist;
#endif
  int k;
  int i,j;
  FILE *histfile;
#ifdef USE_GETUTIME
  long time0,time1;
#endif
  char use_minimg;

#ifdef USE_GETUTIME
  time0 = getutime();
#endif

  if ((x_width == 0.0) || (y_width == 0.0) || (z_width == 0.0)) use_minimg = 0;
  else use_minimg = 1;

  xwidth = ceil(x_max - x_min);
  ywidth = ceil(y_max - y_min);
  zwidth = ceil(z_max - z_min);

  if (use_minimg) {
    celldiag = sqrt(x_width2*x_width2 + y_width2*y_width2 + z_width2*z_width2);
    V = x_width2*y_width2*z_width2;
  }
  else {
    celldiag = sqrt(xwidth*xwidth + ywidth*ywidth + zwidth*zwidth);
    V = xwidth*ywidth*zwidth;
  };

#ifdef VECTORIZE
  for (n=0; n<NUMHIST; n++) {
    radhist[n] = calloc(HISTSIZE,sizeof(long));
    if (radhist[n] == NULL) {
      perror("alloc_radhist");
      exit(10);
    };
  }
  delta = calloc(N,sizeof(Coord3D));
  if (delta == NULL) {
    perror("alloc_delta");
    exit(10);
  };
#else
  radhist = calloc(HISTSIZE,sizeof(long));
  if (radhist == NULL) {
    perror("alloc_radhist");
    exit(10);
  };
#endif

  dr = celldiag/(double)HISTSIZE;

#ifdef VECTORIZE
  if (use_minimg) {
    for (i=minnr; i<=maxnr; i++) {
      for (j=minnr; j<=maxnr; j++) { double rijx, rijy, rijz;
        rijx = atom[i].x-atom[j].x;
        rijy = atom[i].y-atom[j].y;
        rijz = atom[i].z-atom[j].z;

        /* minimal image convention */
        rijx = (rijx-x_width)-x_width*floor((rijx-x_width2)*xw_inv);
        rijy = (rijy-y_width)-y_width*floor((rijy-y_width2)*yw_inv);
        rijz = (rijz-z_width)-z_width*floor((rijz-z_width2)*zw_inv);

        rijsq = square(rijx)+square(rijy)+square(rijz);
        delta[j] = sqrt(rijsq)/dr;
      }
#pragma vdir nodep
      for (j=minnr; j<i; j++) {
        k = (int)delta[j];
        if (k<HISTSIZE) radhist[j % NUMHIST][k]++;
      }
#pragma vdir nodep
      for (j=i+1; j<=maxnr; j++) {
        k = (int)delta[j];
        if (k<HISTSIZE) radhist[j % NUMHIST][k]++;
      }
    }
  }
  else {
    for (i=minnr; i<=maxnr; i++) {
      for (j=minnr; j<=maxnr; j++) { double rijx, rijy, rijz;
        rijx = atom[i].x-atom[j].x;
        rijy = atom[i].y-atom[j].y;
        rijz = atom[i].z-atom[j].z;
        rijsq = square(rijx)+square(rijy)+square(rijz);
        delta[j] = sqrt(rijsq)/dr;
      }
#pragma vdir nodep
      for (j=minnr; j<i; j++) {
        k = (int)delta[j];
        if (k<HISTSIZE) radhist[j % NUMHIST][k]++;
      }
#pragma vdir nodep
      for (j=i+1; j<=maxnr; j++) {
        k = (int)delta[j];
        if (k<HISTSIZE) radhist[j % NUMHIST][k]++;
      }
    }
  }
#else /* not vectorizing */
  if (use_minimg) {
    for (i=minnr; i<=maxnr; i++) {
      for (j=i+1; j<=maxnr; j++) { double rijx, rijy, rijz;
        rijx = atom[i].x-atom[j].x;
        rijy = atom[i].y-atom[j].y;
        rijz = atom[i].z-atom[j].z;
  
        /* minimal image convention */
        rijx = (rijx-x_width)-x_width*floor((rijx-x_width2)*xw_inv);
        rijy = (rijy-y_width)-y_width*floor((rijy-y_width2)*yw_inv);
        rijz = (rijz-z_width)-z_width*floor((rijz-z_width2)*zw_inv);
  
        rijsq = square(rijx)+square(rijy)+square(rijz);
        if ((rijsq/dr) < 0.5) printf("i=%d, j=%d, rijsq=%f\n",i,j,(double)rijsq);
        k = (int)(sqrt(rijsq)/dr);
        if (k==0) printf("k=0 bei i=%d j=%d\n",i,j);
        if (k<HISTSIZE) {
          radhist[k]++;
        }
        else printf("k=%d bei i=%d j=%d\n",k,i,j);
      }
    }
  }
  else
  for (i=minnr; i<=maxnr; i++) {

    for (j=minnr; j<i; j++) { double rijx, rijy, rijz;
      rijx = atom[i].x-atom[j].x;
      rijy = atom[i].y-atom[j].y;
      rijz = atom[i].z-atom[j].z;
      rijsq = square(rijx)+square(rijy)+square(rijz);
      if ((rijsq/dr) < 0.5) printf("i=%d, j=%d, rijsq=%f\n",i,j,(double)rijsq);
      k = (int)(sqrt(rijsq)/dr);
      if (k==0) printf("k=0 bei i=%d j=%d\n",i,j);
      if (k<HISTSIZE) {
        radhist[k]++;
      }
    }

    for (j=i+1; j<=maxnr; j++) { double rijx, rijy, rijz;
      rijx = atom[i].x-atom[j].x;
      rijy = atom[i].y-atom[j].y;
      rijz = atom[i].z-atom[j].z;

      rijsq = square(rijx)+square(rijy)+square(rijz);
      if ((rijsq/dr) < 0.5) printf("i=%d, j=%d, rijsq=%f\n",i,j,(double)rijsq);
      k = (int)(sqrt(rijsq)/dr);
      if (k==0) printf("k=0 bei i=%d j=%d\n",i,j);
      if (k<HISTSIZE) {
        radhist[k]++;
      }
    }
  }
#endif
#ifdef VECTORIZE
  if ((histfile = fopen(out_name,"w")) == NULL) {
    perror("open_radialout");
    exit(10);
  };
  for (k=0; k<HISTSIZE; k++) { double r,r1,prefac;
    sum = 0;
    for (i=0; i<NUMHIST; i++) sum += radhist[i][k];
    r = (double)k*dr;
    r1 = r+dr;
    prefac = V/((double)N * (double)N)/(4/3*PI*((r1*r1*r1)-(r*r*r)));
    fprintf(histfile,"%15.7g %15.7g %15.7g\n",r,(double)sum,(double)sum*prefac);
  }
  fclose(histfile);
  for (k=0; k<NUMHIST; k++) free(radhist[k]);
  free(delta);
#else
  if ((histfile = fopen(out_name,"w")) == NULL) {
    perror("open_radialout");
    exit(10);
  };
  for (k=0; k<HISTSIZE; k++) { double r,r1,prefac,prefac2;
    r = (double)k*dr;
    r1 = r+dr;
    prefac = V/((double)N * (double)N)/(4/3*PI*((r1*r1*r1)-(r*r*r)));
    prefac2 = V/((double)N * (double)N)/(4*PI*(r*r));
    fprintf(histfile,"%15.7g %15.7g %15.7g %15.7g\n",r,(double)radhist[k],(double)radhist[k]*prefac,(double)radhist[k]*prefac2);
  }
  fclose(histfile);
  free(radhist);
#endif
#ifdef USE_GETUTIME
  time1 = getutime()-time0;
  printf("Time for radialdist() = %lu.%02lu\n",time1/Ticks,(time1 % Ticks)*(100L/Ticks));
#endif
} /* radialdist */

/* calculate the length of vector <a> */
REAL vbetrag(const VEC a)
{
  REAL sum;
  int i;
  sum = 0.0;
  for (i=0; i<3; i++) sum += a[i]*a[i];
  return sqrt(sum);
} /* vbetrag */

/* Skalarprodukt (inneres Produkt) zweier Vektoren <a> und <b> */
/* calculate the dot product (inner product) of <a> and <b> */
REAL vdot(const VEC a, const VEC b)
{
  REAL sum;
  sum = a[0]*b[0];
  sum += a[1]*b[1];
  sum += a[2]*b[2];
  return sum;
} /* vdot */

/* subtract 3D vector <b> from vector <a> and put the result in <c> */
void vecsub(VEC c, const VEC a, const VEC b)
{
  c[0] = a[0]-b[0];
  c[1] = a[1]-b[1];
  c[2] = a[2]-b[2];
}

/* add the vectors <a> and <b>, and put the result in <c> */
void vecadd(VEC c, const VEC a, const VEC b)
{
  c[0] = a[0]+b[0];
  c[1] = a[1]+b[1];
  c[2] = a[2]+b[2];
}

/* scale the vector <a> with factor <scale> and put the result in <c> */
void vecscale(VEC c, const VEC a, const REAL scale)
{
  c[0] = a[0]*scale;
  c[1] = a[1]*scale;
  c[2] = a[2]*scale;
}

/* copy a 3D vector from source to dest */
void vcopy(VEC dest, const VEC source)
{
  dest[0] = source[0];
  dest[1] = source[1];
  dest[2] = source[2];
}

/* Set up rotation matrix for rotation around x-axis with angle phi */
void setrot_x(MAT U, const REAL phi)
{
  U[0][0]=1;
  U[0][1]=0;
  U[0][2]=0;
  U[1][0]=0;
  U[1][1]=cos(phi);
  U[1][2]=-sin(phi);
  U[2][0]=0;
  U[2][1]=sin(phi);
  U[2][2]=cos(phi);
} /* setrot_x */

/* Set up rotation matrix for rotation around y-axis with angle phi */
void setrot_y(MAT U, const REAL phi)
{
  U[0][0]=cos(phi);
  U[0][1]=0;
  U[0][2]=sin(phi);
  U[1][0]=0;
  U[1][1]=1;
  U[1][2]=0;
  U[2][0]=-sin(phi);
  U[2][1]=0;
  U[2][2]=cos(phi);
} /* setrot_y */

/* Set up rotation matrix for rotation around z-axis with angle phi */
void setrot_z(MAT U, const REAL phi)
{
  U[0][0]=cos(phi);
  U[0][1]=sin(phi);
  U[0][2]=0;
  U[1][0]=-sin(phi);
  U[1][1]=cos(phi);
  U[1][2]=0;
  U[2][0]=0;
  U[2][1]=0;
  U[2][2]=1;
} /* setrot_z */

/* Matrix-Transform of a vector <a>:   x = U a; */
/* U is an orthogonal 3 x 3 - Matrix */
void vmatmul(const MAT U, const VEC a, VEC x)
{
  int i,j;
  for (i=0; i<3; i++) {
    x[i]=0;
    for (j=0; j<3; j++) {
      x[i] += U[i][j]*a[j];
    }
  }
} /* vmatmul */

/* C_ik = A_ij * B_jk */
void mmatmul(const MAT A, const MAT B, MAT C)
{
  double s;
  int i,j,k;
  for (i=0; i<3; i++) {
    for (k=0; k<3; k++) {
      s = 0;
      for (j=0; j<3; j++) s += A[i][j]*B[j][k];
      C[i][k] = s;
    }
  }
} /* mmatmul */

#if 0
/* Eine Drehung um einen beliebigen Punkt laesst sich zerlegen in eine
   Translation (die diesen Punkt zum neuen Nullpunkt macht), eine
   Drehung um den Nullpunkt, und die inverse Translation */

/* Translation */
void translate(Coord3D *atom, int N, Coord3D *displace)
{
  int i;
  for (i=0; i<N; i++) {
    atom[i].x += displace->x;
    atom[i].y += displace->y;
    atom[i].z += displace->z;
  }
}
#endif

/* Drehung um Nullpunkt */
void rotate(Coord3D *atom, struct AtomDaten *atomdaten, int N, struct RotateAngles *rot)
{
  int i;
  MAT U,U1,U2; /* Drehmatrix */
  MAT Unity = {{1,0,0},{0,1,0},{0,0,1}};
  VEC a,x;

  if (rot->phi_x != 0) setrot_x(U,rot->phi_x);
  else memmove(U,Unity,sizeof(MAT));
  if (rot->phi_y != 0) {
    setrot_y(U1,rot->phi_y);
    mmatmul(U,U1,U2);
    memmove(U,U2,sizeof(MAT));
  };
  if (rot->phi_z != 0) {
    setrot_z(U1,rot->phi_z);
    mmatmul(U,U1,U2);
    memmove(U,U2,sizeof(MAT));
  };
  for (i=minnr; i<=maxnr; i++) {
    a[0]=atom[i].x;
    a[1]=atom[i].y;
    a[2]=atom[i].z;
    
    vmatmul(U,a,x);
    
    atom[i].x = x[0];
    atom[i].y = x[1];
    atom[i].z = x[2];
  }
} /* rotate */

/* Clip a rectangular box */
void cliprbox(Coord3D **atom, struct AtomDaten **atomdaten, int *numatoms, struct ClipRange *clip, unsigned clipmask)
{
  Coord3D *atom_tmp;
  struct AtomDaten *atomdaten_tmp;  
  int N;
  int i,j;

  N = *numatoms;
  atom_tmp = calloc(N,sizeof(Coord3D));
  atomdaten_tmp = calloc(N,sizeof(struct AtomDaten));

  j=0;
  for (i=minnr; i<maxnr; i++) {
    if (clipmask & 1) if ((*atom)[i].x < clip->min_x) continue;
    if (clipmask & 2) if ((*atom)[i].x > clip->max_x) continue;
    if (clipmask & 4) if ((*atom)[i].y < clip->min_y) continue;
    if (clipmask & 8) if ((*atom)[i].y > clip->max_y) continue;
    if (clipmask & 0x10) if ((*atom)[i].z < clip->min_z) continue;
    if (clipmask & 0x20) if ((*atom)[i].z > clip->max_z) continue;

    atom_tmp[j]=(*atom)[i];
    atomdaten_tmp[j]=(*atomdaten)[i];
    j++;
  }
  free(*atom);
  free(*atomdaten);
  *atom = atom_tmp;
  *atomdaten = atomdaten_tmp;
  *numatoms = j;
  minnr = 0;
  maxnr = j-1;
} /* cliprbox */

/* Stapeln der eingelesenen Atomschicht(en) in z-Richtung */
void StackLayers(Coord3D **atome,struct AtomDaten **atomdaten, int *numatoms, double z_incr, int z_times)
{
  Coord3D *atom_tmp;
  struct AtomDaten *atomdaten_tmp;  
  int k,i;
  int N;

  N = *numatoms;
  atom_tmp = calloc(N*z_times,sizeof(Coord3D));
  atomdaten_tmp = calloc(N*z_times,sizeof(struct AtomDaten));

  for (k=0; k<z_times; k++) {
#if 1
    memcpy(&atom_tmp[k*N],&(*atome)[minnr],sizeof(Coord3D)*(maxnr-minnr+1));
    memcpy(&atomdaten_tmp[k*N],&(*atomdaten)[minnr],sizeof(struct AtomDaten)*(maxnr-minnr+1));
#endif
    for (i=0; i<maxnr-minnr+1; i++) {
#if 0
      atom_tmp[k*N+i]=(*atome)[i+minnr];
      atomdaten_tmp[k*N+i]=(*atomdaten)[i+minnr];
#endif
      atom_tmp[k*N+i].z += (double)k * z_incr;
    };
  }
  free(*atome);
  free(*atomdaten);
  *atome = atom_tmp;
  *atomdaten = atomdaten_tmp;
  *numatoms = N*z_times;
  minnr = 0;
  maxnr = N*z_times;
} /* StackLayers */

void AlCuCo_DiagAngle(struct RotateAngles *rota)
{
  const VEC A = {49.049999,13.800000,0}, B={0,51.299999,0};
  REAL x,y;
  REAL phi;

  /* Drehwinkel fuer Schnitt entlang langer Diagonale des rauten-
    foermigen AlCuCo-Approximanten ausrechnen (zum Ausschneiden
    eines rechteckigen Teilbereichs) */

/*
  A[0]=49.049999;
  A[1]=13.800000;
  A[2]=0;

  B[0]=0;
  B[1]=51.299999;
  B[2]=0;
*/

  x = vdot(A,B);
  y = vbetrag(A)*vbetrag(B);
  phi = PI/2.0-acos(x/y)/2.0;

  rota->phi_x = 0;
  rota->phi_y = 0;
  rota->phi_z = phi;
/*  printf("x=%g y=%g phi=%g\n",x,y,phi); */
} /* AlCuCo_DiagAngle */

void AlCuCo_ClipSize(struct ClipRange *clip, unsigned *clipmask)
{
  const VEC A = {49.049999,13.800000,0}, B={0,51.299999,0};
  VEC C,D,E,F,G;
  REAL l1,l2,l3;
  vecscale(C,A,0.5);
  vecscale(D,B,0.5);
  vecadd(E,C,D);
  l1=vbetrag(E);
  vecsub(F,D,C);
  l2=vbetrag(F);
  vecscale(C,C,0.5);
  vecscale(D,D,0.5);
  vecadd(G,C,D);
  l3=vbetrag(G);
  clip->min_x = l3;
  clip->max_x = l3+l1;
  clip->min_y = -l2/2.0;
  clip->max_y = l2/2.0;
  *clipmask = 0x0F;
} /* AlCuCo_ClipSize */

void AllocClip(struct ClipRange **clip)
{
  *clip = malloc(sizeof(struct ClipRange));
  if (*clip == NULL) {
    perror("AllocCLip");
    exit(10);
  }
} /* AllocClip */              

int main(int argc,char **argv)
{
  char *in_name = 0L, *out_name = 0L;
  FILE *in, *out;
  unsigned num_atoms;
  unsigned i,k;
#if 0
  unsigned wahl;
#endif
  char inputformat = 0; /* 0 --> read imd format */
  char format = 0; /* 0 --> write molgraph format */
  char noninteractive = FALSE;
#ifdef USE_GETUTIME
  long time1,time0;
#endif
  char *parameterfile;
  struct RotateAngles rota = {0,0,0};
  unsigned clipmaske = 0; /* was alles soll abgeschnitten werden */
  struct ClipRange *clip = 0L;
  double z_incr = 0;
  int z_times = 0;

  keep = FALSE;
  atome = (Coord3D*)calloc(sizeof(Coord3D),ALLOC_SIZE);
  if (atome == NULL) {
    perror("Allocate_Atoms");
    exit(10);
  }
  atomdaten = (struct AtomDaten*)calloc(sizeof(struct AtomDaten),ALLOC_SIZE);
  if (atomdaten == NULL) {
    perror("Allocate_AtomDaten");
    exit(10);
  }

  radius = (double*)calloc(sizeof(double),MAX_SORT);
  if (radius == NULL) {
    perror("AllocRadius");
    exit(10);
  };
  for (i=0; i<MAX_SORT; i++) {
    radius[i]=1.0;
  };

  k = 0;
  for (i=1; i<argc; i++) {
    char *argstr = argv[i];
    if (*argstr == '-') {
      switch(*(++argstr)) {
        case 'n' : noninteractive = TRUE; break;
        case 'k' : keep = TRUE; break;
        case 'C' :
          ++argstr;
          if (strcmp(argstr,"xmin")==0) {
            if (clip == 0L) AllocClip(&clip);
            clip->min_x = atof(argv[++i]);
            clipmaske |= 1;
          }
          else if (strcmp(argstr,"xmax")==0) {
            if (clip == 0L) AllocClip(&clip);
            clip->max_x = atof(argv[++i]);
            clipmaske |= 2;
          }
          else if (strcmp(argstr,"ymin")==0) {
            if (clip == 0L) AllocClip(&clip);
            clip->min_y = atof(argv[++i]);
            clipmaske |= 4;
          }
          else if (strcmp(argstr,"ymax")==0) {
            if (clip == 0L) AllocClip(&clip);
            clip->max_y = atof(argv[++i]);
            clipmaske |= 8;
          }
          else if (strcmp(argstr,"zmin")==0) {
            if (clip == 0L) AllocClip(&clip);
            clip->min_z = atof(argv[++i]);
            clipmaske |= 0x10;
          }
          else if (strcmp(argstr,"zmax")==0) {
            if (clip == 0L) AllocClip(&clip);
            clip->max_z = atof(argv[++i]);
            clipmaske |= 0x20;
          }
          else if (strcmp(argstr,"AlCuCo")==0) {
            if (clip == 0L) AllocClip(&clip);
            AlCuCo_ClipSize(clip,&clipmaske);
          };
          break;
        case 'f' :
          ++argstr;
          switch(*argstr) {
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
              format=atoi(argstr); break;
            case 'R':
              if (strcmp(argstr,"RDF")==0) {
                format = *argstr; break;
              };
            default:
              format = 0;
          };
          break;
        case 'i' : inputformat=atoi(++argstr); break;
        case 'P' : parameterfile=strdup(++argstr); break;
        case 'R' :
          ++argstr;
          if (strcmp(argstr,"x")==0) rota.phi_x = atof(argv[++i]);
          else if (strcmp(argstr,"y")==0) rota.phi_y = atof(argv[++i]);
          else if (strcmp(argstr,"z")==0) rota.phi_z = atof(argv[++i]);
          else if (strcmp(argstr,"AlCuCo")==0) AlCuCo_DiagAngle(&rota);
          else fputs("What?\n",stderr);
          break;
        case 'S' : z_incr = atof(argv[++i]); z_times=atoi(argv[++i]); break;
        case 'X' : x_width = atof(++argstr); x_width2=x_width/2.0; xw_inv=1.0/x_width; break;
        case 'Y' : y_width = atof(++argstr); y_width2=y_width/2.0; yw_inv=1.0/y_width; break;
        case 'Z' : z_width = atof(++argstr); z_width2=z_width/2.0; zw_inv=1.0/z_width; break;
      default:
        fprintf(stderr,"Unknown option %s\n",argstr);
      };
    }
    else {
      if (k==0) in_name = argstr;
      else if (k==1) out_name = argstr;
      else {
        fprintf(stderr,"Extraneous argument %s ignored\n",argstr);
        printf("k=%d\n",k);
      }
      k++;
    }
  }
  if (k==0) {
    Usage(argv[0]);
    exit(10);
  };

#ifdef USE_GETUTIME
  time0 = getutime();
#endif
  in = fopen(in_name,"r");
  if (in != NULL) {
    printf("Datei %s wird eingelesen\n",in_name);
    if (inputformat == 0) num_atoms=lies_atome(in,&atome,&atomdaten);
    else if (inputformat==1) num_atoms=lies_atome_jr(in,&atome,&atomdaten);
    fclose(in);
#ifdef USE_GETUTIME
    time1 = getutime()-time0;
#endif
  }
  else { /* Pech gehabt ! */
    printf("file %s couldn't be opened\n",in_name);
    perror("lies_atome");
    exit(10);
  }
#ifdef USE_GETUTIME
  printf("%u atoms were read in %lu.%02lu seconds\n",num_atoms,
         time1/Ticks,(time1 % Ticks)*(100L/Ticks));
#else
  printf("%u atoms were read\n",num_atoms);
#endif
#if 0 /* Historische Relikte ... ;-) */
  if (!noninteractive) {
    do {
      printf("1 - Write out data file for molgraph\n");
      printf("0 - Quit\n");
      do {
        printf("Your choice: ");
        scanf("%d",&wahl);
      } while (wahl>4);

      switch(wahl) {
        case 1 : writemgdat(out,out_name,atome,num_atoms); break;
        default: goto quit;
      }
    } while (wahl!=0);
  }
  else
#endif
  {
    if (z_times > 0) StackLayers(&atome,&atomdaten,&num_atoms,z_incr,z_times);
    if ((rota.phi_x != 0) || (rota.phi_y != 0) || (rota.phi_z != 0)) rotate(atome,atomdaten,num_atoms,&rota);
    if (clipmaske != 0) cliprbox(&atome,&atomdaten,&num_atoms,clip,clipmaske);
    if (format == 0) writemgdat(out_name,atome,atomdaten,num_atoms);
    else if (format == 1) writepdb(out_name,atome,atomdaten,num_atoms);
    else if (format == 2) writecfg(out_name,atome,atomdaten,num_atoms);
    else if (format == 'R') radialdist(out_name,atome,atomdaten,num_atoms);
  }
quit:
  return 0;
}
