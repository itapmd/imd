#include <stdio.h>
#include <string.h>
#define TRUE 1
#define FALSE 0

int main(int argc, char **argv) {

  FILE *fpin, *fpout;
  int i, twod;
  char intype[2],outtype[3], line[255];
  char argument[255], infname[255], outfname[255];
  float x,y,z, mass, vx, vy, vz, prop;
  float ix,iy,iz, imass, ivx, ivy, ivz, iprop;
  float ox,oy,oz, omass, ovx, ovy, ovz, oprop;
  float boxx,boxy,boxz,hboxx,hboxy,hboxz;
  int inr, itype;
  int onr, otype;
  int nr, type;
  int velocs, props, imdanz,in_binary,out_binary;
  char infpstr[3];
  char outfpstr[3];
  char usagestring[2000];

  sprintf(usagestring, "atcon <in_format>:<out_format> [-2d] [-bi] [-bo] [-help] [-i] <infile> [-o] <outfile>\n\nSupported input formats:\nTEE,IMD,MOP\nSupported output formats:\nTEE,IMD,MOP,POS,VRM,PDB\n\nTEE: < x y z t >, the first line contains the box\nIMD: < nr t m x y z (vx vy vz) (prop) >\nMOP:\nPOS: postscript\nVRM: vrml\nPDB: proteine database format (rasmol)\n");

  strcpy(infpstr,"r");
  strcpy(outfpstr,"w");

  if (argc < 2) {
    printf("Usage: %s",usagestring);
    exit(1);
  }

  for (i=1;i<argc;i++) {
    strcpy(argument, argv[i]);
    if (strncmp(argument,"-",1)==0) {
      if (strcmp(argument, "-2d")==0) {
	twod=TRUE;
	continue;
      }
      if (strcmp(argument, "-h")==0) {
	printf("Usage: %s",usagestring);
	exit(1);
      }
      if (strcmp(argument, "-help")==0) {
	printf("Usage: %s",usagestring);
	exit(1);
      }
      if (strcmp(argument, "-bi")==0) {
	strcpy(infpstr,"rb");
	continue;
      }
      if (strcmp(argument, "-bo")==0) {
	strcpy(outfpstr,"wb");
	continue;
      }
      if (strcmp(argument, "-o")==0) {
	strcpy(outfname, argv[++i]);
	continue;
      }
      if (strcmp(argument, "-i")==0) {
	strcpy(infname, argv[++i]);
	continue;
      }
      printf("Unknown option %s. Usage: %s",argument,usagestring);
      exit(-1);
    }
    else {
      if (argument[3]==':') {
	strncpy(intype,argument,3);
        strcpy(outtype, strchr(argument,':'));
      }
    }
  }
  
  if (infname==0) { printf("No input file specified\n"); exit(2); }
  if (outfname==0) { printf("No output file specified\n"); exit(3); }
  printf("%s", infpstr);
  fpin=fopen(infname, infpstr);
  if (fpin==NULL) { printf("Cannot open input file %s\n",infname); exit(4); }
  fpout=fopen(outfname, outfpstr);
  if (fpout==NULL) { printf("Cannot open output file %s\n", outfname); exit(5); }

  /* beginning of output files */
  if (strncmp(outtype,":POS",4)==0) {
    fprintf(fpout,"%!\n");
    fprintf(fpout,"72 2.54 div dup scale 0.01 setlinewidth 2 2 translate\n");
    fprintf(fpout,"/n{newpath}def /m{moveto}def /rls{rlineto stroke}def\n");
    fprintf(fpout,"/dot{n 0.15 0 360 arc closepath fill}def\n");
    fprintf(fpout,"0.5 0.5 scale\n");  
  }
  if (strncmp(outtype,":VRM",4)==0) {
    fprintf(fpout,"#VRML V1.0 ascii\n");
  }


  nr=0;
  while(fgets(line, 200, fpin)) {

    /* inits */
    ++nr;
    mass=0;
    if (strncmp(intype,"IMD",3)==0) {
      velocs=FALSE;
      props=FALSE;
    }

    /* read and convert to independent quantities */
    if (strncmp(intype,"TEE",3)==0) {
      /* box */
      if (nr==1) {
	sscanf(line, "%f %f %f", &boxx, &boxy, &boxz);
	hboxx=.5*boxx;
	hboxy=.5*boxy;
	hboxz=.5*boxz;
	fgets(line,200,fpin);
      }
      sscanf(line, "%f %f %f %d", &ix, &iy, &iz, &itype);
      type=itype-1;
      x=ix+hboxx;
      y=iy+hboxy;
      z=iz+hboxz;
    }
    if (strncmp(intype,"IMD",3)==0) {
      imdanz=sscanf(line, "%d %d %f %f %f %f", &inr, &itype, &imass, &ix, 
	     &iy, &iz,&ivx,&ivy,&ivz,&iprop);
      nr=inr;
      x=ix;
      y=iy;
      z=iz;
      vx=ivx;
      vy=ivy;
      vz=ivz;
      prop=iprop;
      type=itype;
      mass=imass;
    }

    /* clean up read */
    switch(imdanz) {
    case 9:  velocs=TRUE;break;
    case 10: velocs=TRUE;props=TRUE;break;
    }

    /* all over conversions */
    ox=x;
    oy=y;
    oz=z;
    if (velocs) {
      ovx=vx;
      ovy=vy;
      ovz=vz;
    }
    if (props) {
      prop=prop;
    }

    /* convert from independent quantities and print */
    if (strncmp(outtype,":TEE",4)==0) {
      otype=type-1;
      fprintf(fpout, "%f %f %f %d\n", ox, oy, oz, otype);
    }
    if (strncmp(outtype,":IMD",4)==0) {
      onr=nr;
      otype=type;
      omass=mass>0?mass:1.0;
      if (props)
	oprop=prop;
      fprintf(fpout, "%d %d %f %f %f %f", onr, otype, omass, ox, oy, oz);
      if (velocs) fprintf(fpout, " %f %f %f", ovx, ovy, ovz);
      if (props) fprintf(fpout, " %f", oprop);
      fprintf(fpout, "\n");
    }
    if (strncmp(outtype,":POS",4)==0) {
      fprintf(fpout, "n %f %f dot\n", ox, oy);
    }
    if (strncmp(outtype,":VRM",4)==0) {
      otype=type;
      fprintf(fpout, "  Separator {\n    Translation { translation %f %f %f }\n    Material { diffuseColor %d %d %d }\n    Sphere { radius %f }\n  }\n", ox, oy, oz, otype, otype, otype,0.1*(otype+1));
    }
    if (strncmp(outtype,":PDB",4)==0) {
      onr=nr;
      otype=type+1;
      oprop=prop>0?prop:1.0;
      fprintf(fpout,
	      "ATOM  %5d  AL       %4d     %7.3f %7.3f %7.3f  1.00%6.2f\n",
	      onr, otype, ox, oy, oz, oprop);
    }
  }

  fclose(fpin);

  /* end of outfiles */
  if (strncmp(outtype,":POS",4)==0) {
    fprintf(fpout, "print \"showpage\"\n");
  }
  
  fclose(fpout);
}









