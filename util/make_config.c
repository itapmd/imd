
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* make_config.c -- generate configurations
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MODE_HEX         1
#define MODE_FCC         2
#define MODE_BCC         3
#define MODE_B2          4
#define MODE_NACL        5
#define MODE_DIAMOND     6
#define MODE_ZINCBLENDE  7
#define MODE_LAV         8

typedef double real;
typedef struct { real x; real y; real z; }  vektor;
typedef struct { int  x; int  y; int  z; } ivektor;

ivektor box_param = {1, 1, 1};
real    box_unit = 1, masses[2] = {1.0, 1.0};
vektor  box_x, box_y, box_z;
int     mode = 0, natoms = 0;
char    *outfilename;
FILE    *outfile;

/* online help message */
void usage()
{
  printf("\n");
  printf("   Usage: make_config <type> <outfile> [options] \n\n");
  printf("   <type> must be one of:\n\n");
  printf("    - hex        2D hexagonal crystal\n");
  printf("    - fcc        fcc structure\n");
  printf("    - bcc        bcc structure\n");
  printf("    - b2         B2 structure\n");
  printf("    - nacl       NaCl structure\n");
  printf("    - diamond    cubic diamond structure\n");
  printf("    - zincblende zincblende structure\n");
  printf("    - lav        cubic C15 Laves structure (MgCu2)\n\n");
  printf("   Options:  -s <sx> <sy> [<sz>]  size in unit cells (default 1)\n");
  printf("             -a <a>               lattice constant (default 1.0)\n");
  printf("             -m <m0> [<m1>]       mass(es)         (default 1.0)\n");
  printf("\n");
  exit(1);
}

/* error and usage message */
void error(char *msg)
{
  printf("\n   %s\n",msg);
  usage();
}

/* initialize for hexagonal crystal */
void init_hex(void)
{
  if ((box_param.x==0) || (box_param.y==0)) error("box_param not set!");
  box_x.x = box_param.x * sqrt(3.0) * box_unit;
  box_x.y = 0.0;
  box_y.x = 0.0;
  box_y.y = box_param.y * box_unit;
  /* write header */
  fprintf(outfile, "#F A 1 1 1 2 0 0\n");
  fprintf(outfile, "#C number type mass x y\n");
  fprintf(outfile, "#X %f %f\n", box_x.x, box_x.y);
  fprintf(outfile, "#Y %f %f\n", box_y.x, box_y.y);
  fprintf(outfile, "#E\n");
}

/* generate 2D hexagonal crystal */
void generate_hex()
{
  vektor  min, max;
  int     i, j, typ;
  real    x, y;

  min.x = 0; max.x = 2 * box_param.x;
  min.y = 0; max.y = 2 * box_param.y;

  natoms  = 0;

  for (i=min.x ; i<max.x; i++)
    for (j=min.y; j<max.y; j++) {

      typ  = (i+j) % 2;
      if (typ > 0) continue;

      x = (i+0.5) * sqrt(3.0) * 0.5 * box_unit;
      y = (j+0.5) * 0.5 * box_unit;
      natoms++;

      fprintf(outfile,"%d %d %f %f %f\n", natoms, typ, masses[typ], x, y);
  }
}

/* initialize for cubic crystals */
void init_cubic(void)
{
  box_x.x = box_param.x * box_unit;  box_x.y = 0.0;  box_x.z = 0.0;
  box_y.x = 0.0;  box_y.y = box_param.y * box_unit;  box_y.z = 0.0;
  box_z.x = 0.0;  box_z.y = 0.0;  box_z.z = box_param.z * box_unit;
  if ((mode==MODE_FCC) || (mode==MODE_BCC)  || 
      (mode==MODE_B2)  || (mode==MODE_NACL)) {
    box_param.x *= 2;
    box_param.y *= 2;
    box_param.z *= 2;
    box_unit    /= 2;
  } else if ((mode==MODE_DIAMOND) || (mode==MODE_ZINCBLENDE)) {
    box_param.x *= 4;
    box_param.y *= 4;
    box_param.z *= 4;
    box_unit    /= 4;
  }
  /* write header */
  fprintf(outfile, "#F A 1 1 1 3 0 0\n");
  fprintf(outfile, "#C number type mass x y z\n");
  fprintf(outfile, "#X %f %f %f\n", box_x.x, box_x.y, box_x.z);
  fprintf(outfile, "#Y %f %f %f\n", box_y.x, box_y.y, box_y.z);
  fprintf(outfile, "#Z %f %f %f\n", box_z.x, box_z.y, box_z.z);
  fprintf(outfile, "#E\n");
}

/* generate cubic crystals (except Laves) */
void generate_fcc(int maxtyp)
{
  vektor  min, max;
  int     x, y, z, typ;
  real    xx, yy, zz;

  min.x = 0; max.x = box_param.x;
  min.y = 0; max.y = box_param.y;
  min.z = 0; max.z = box_param.z;

  natoms = 0;

  for (x=min.x ; x<max.x; x++)
    for (y=min.y; y<max.y; y++)
      for (z=min.z; z<max.z; z++) {
 
        typ  = (x+y+z) % 2;

	/* cubic diamond and zincblende case */
	if (maxtyp == 4 || maxtyp == 5) {
	    if ( ((x+y+z)%4==0) && 
		 (z%2==0) && (y%2==0) && (x%2==0)  ) {
		typ=0;
	    }
	    else if ( ((x+y+z)%4==3) && 
		      (z%2==1) && (y%2==1) && (x%2==1) ) {
		if (maxtyp == 4)
		    typ=0;
		else
		    typ=1;
	    }
	    else continue;
	}

        /* B2 == CsCl structure */
        if (maxtyp ==3) {
          if ((z%2==0) && (y%2==0) && (x%2==0)) {
             typ=0;
          }
          else if ((z%2==1) && (y%2==1) && (x%2==1)) {
             typ=1;
          }
          else continue;
        }

        /* bcc case - all atoms get typ 0 */
        if (maxtyp == 2) {
          if (((x%2) != (y%2)) || ((y%2) != (z%2))) continue;
          typ = 0;
        }

        /* maxtyp == 0: fcc;  maxtyp == 1: NaCl */
        if (typ > maxtyp) continue;  /* if fcc, only atoms of type 0 */

	natoms++;
        xx = (x + 0.5) * box_unit;
        yy = (y + 0.5) * box_unit;
        zz = (z + 0.5) * box_unit;
        fprintf(outfile, "%d %d %f %f %f %f\n", 
                natoms, typ, masses[typ], xx, yy, zz);
  }
}

/* generate a cubic Laves structure crystal */
void generate_lav()
{
  vektor  min, max;
  real    px[24], py[24], pz[24];
  real    x, y, z;
  int     i, j, k, l, typ, pa[24];

  px[ 0] = 0; py[ 0] = 0; pz[ 0] = 0; pa[ 0] = 0;
  px[ 1] = 2; py[ 1] = 2; pz[ 1] = 0; pa[ 1] = 0;
  px[ 2] = 2; py[ 2] = 0; pz[ 2] = 2; pa[ 2] = 0;
  px[ 3] = 0; py[ 3] = 2; pz[ 3] = 2; pa[ 3] = 0;
  px[ 4] = 3; py[ 4] = 3; pz[ 4] = 3; pa[ 4] = 1;
  px[ 5] = 5; py[ 5] = 5; pz[ 5] = 5; pa[ 5] = 1;

  px[ 6] = 4; py[ 6] = 4; pz[ 6] = 0; pa[ 6] = 0;
  px[ 7] = 6; py[ 7] = 6; pz[ 7] = 0; pa[ 7] = 0;
  px[ 8] = 6; py[ 8] = 4; pz[ 8] = 2; pa[ 8] = 0;
  px[ 9] = 4; py[ 9] = 6; pz[ 9] = 2; pa[ 9] = 0;
  px[10] = 7; py[10] = 7; pz[10] = 3; pa[10] = 1;
  px[11] = 1; py[11] = 1; pz[11] = 5; pa[11] = 1;

  px[12] = 4; py[12] = 0; pz[12] = 4; pa[12] = 0;
  px[13] = 6; py[13] = 2; pz[13] = 4; pa[13] = 0;
  px[14] = 6; py[14] = 0; pz[14] = 6; pa[14] = 0;
  px[15] = 4; py[15] = 2; pz[15] = 6; pa[15] = 0;
  px[16] = 7; py[16] = 3; pz[16] = 7; pa[16] = 1;
  px[17] = 1; py[17] = 5; pz[17] = 1; pa[17] = 1;

  px[18] = 0; py[18] = 4; pz[18] = 4; pa[18] = 0;
  px[19] = 2; py[19] = 6; pz[19] = 4; pa[19] = 0;
  px[20] = 2; py[20] = 4; pz[20] = 6; pa[20] = 0;
  px[21] = 0; py[21] = 6; pz[21] = 6; pa[21] = 0;
  px[22] = 3; py[22] = 7; pz[22] = 7; pa[22] = 1;
  px[23] = 5; py[23] = 1; pz[23] = 1; pa[23] = 1;

  min.x = 0; max.x = box_param.x;
  min.y = 0; max.y = box_param.y;
  min.z = 0; max.z = box_param.z;

  natoms  = 0;

  box_unit /= 8.0;
  for (i=min.x; i<max.x; i++)
    for (j=min.y; j<max.y; j++)
      for (k=min.z; k<max.z; k++) 
	for (l=0; l<24; l++) {
	  natoms++;
	  typ = pa[l];
	  x = (px[l] + 8*i + 0.5) * box_unit;
	  y = (py[l] + 8*j + 0.5) * box_unit;
	  z = (pz[l] + 8*k + 0.5) * box_unit;
          fprintf(outfile, "%d %d %f %f %f %f\n", 
                  natoms, typ, masses[typ], x, y, z);
        }
} 

int main( int argc, char **argv ) 
{
  char *modestr;

  if (argc<3) usage(); 
  modestr = strdup(argv[1]);
  outfilename = strdup(argv[2]);

  if      (0 == strcmp(modestr, "hex"))        mode = MODE_HEX; 
  else if (0 == strcmp(modestr, "fcc"))        mode = MODE_FCC;
  else if (0 == strcmp(modestr, "bcc"))        mode = MODE_BCC;
  else if (0 == strcmp(modestr, "b2"))         mode = MODE_B2;
  else if (0 == strcmp(modestr, "nacl"))       mode = MODE_NACL;
  else if (0 == strcmp(modestr, "diamond"))    mode = MODE_DIAMOND;
  else if (0 == strcmp(modestr, "zincblende")) mode = MODE_ZINCBLENDE;
  else if (0 == strcmp(modestr, "lav"))        mode = MODE_LAV;
  else error("Unknown structure type!");
  argc -= 2;
  argv += 2;

  /* parse options */
  while ((argc > 1) && (argv[1][0] =='-')) {
    if (argv[1][1]=='s') {
      if (mode == MODE_HEX) {
        if ((argc<4) || (argv[2][0] == '-') || (argv[3][0] == '-'))  
          error("Not enough parameters!");
        box_param.x = atoi(argv[2]);
        box_param.y = atoi(argv[3]);
        argc -= 3;
        argv += 3;
      } else {
        if ((argc<5) || (argv[2][0] == '-') || (argv[3][0] == '-')
                     || (argv[4][0] == '-')) error("Not enough parameters!");
        box_param.x = atoi(argv[2]);
        box_param.y = atoi(argv[3]);
        box_param.z = atoi(argv[4]);
        argc -= 4;
        argv += 4;
      }
    }
    else if (argv[1][1]=='a') {
      if ((argc<3) || (argv[2][0] == '-')) error("Not enough parameters!");
      box_unit = atof(argv[2]);
      argc -= 2;
      argv += 2;
    }
    else if (argv[1][1]=='m') {
      if ((mode==MODE_NACL)       || (mode==MODE_B2) || 
          (mode==MODE_ZINCBLENDE) || (mode==MODE_LAV) ) {
        if ((argc<4) || (argv[2][0] == '-') || (argv[3][0] == '-'))  
          error("Not enough parameters!");
        masses[0] = atof(argv[2]);
        masses[1] = atof(argv[3]);
        argc -= 3;
        argv += 3;
      } else {
        if ((argc<3) || (argv[2][0] == '-')) 
          error("Not enough parameters!");
        masses[0] = atof(argv[2]);
        argc -= 2;
        argv += 2;
      }
    }
    else error("Unknown option!\n");
  }

  /* open output file */
  outfile = fopen(outfilename,"w");
  if (NULL==outfile) error("Cannot open output file.");

  /* generate structure */
  switch (mode) {
    case MODE_HEX:
      /* 2D hexagonal */
      init_hex();
      generate_hex();
      break;
    case MODE_FCC:
       /* fcc */
      init_cubic();
      generate_fcc(0);
      break;
    case MODE_NACL:
      /* NaCl structure */
      init_cubic();
      generate_fcc(1);
       break;
   case MODE_BCC:
      /* bcc */
      init_cubic();
      generate_fcc(2);
       break;
   case MODE_B2:
      /* B2 (CsCl structure) */
      init_cubic();
      generate_fcc(3);
      break;
    case MODE_DIAMOND:
      /* cubic diamond */
      init_cubic();
      generate_fcc(4);
      break;
    case MODE_ZINCBLENDE:
      /* zincblende structure (binary) */
      init_cubic();
      generate_fcc(5);
      break;
    case MODE_LAV:
      /* C15 Laves phase (binary) */
      init_cubic();
      generate_lav();
      break;
  }
  fclose(outfile);

  printf("Generated %s structure with %d atoms.\n", modestr, natoms);
  return 0;
}

