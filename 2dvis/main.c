#include <stdio.h>
#include <string.h>
#include <vogle.h>
#include <math.h>

#define MAIN

#include "prototypes.h"
#include "globals.h"
#include "makros.h"


int main(int argc, char **argv)
{
  int i;
  char fname[255];

  /* inits */
  putenv("USEOWNCMAP=1");
  strcpy(uvfname,"unitvecs");

  read_parameters(argc,argv);
  window_main(argc,argv);
}

/* init_graph initializes the graphics window */
void init_graph(void) {

  char device[10];
  int nplanes,i;

  strcpy(device, "X11");
  prefsize(768,640);
  vinit(device);                 /* Grafikwindow */
  nplanes = getdepth();          /* Grafikplanes */
  
  color(BLACK);
  clear();
  polyfill(1);

}

/* drawing of text */
void draw_text(void) {

  char str1[200],str2[200],str3[300];

  color(CYAN);
  if (scene_type)
    sprintf(str1,"Distribution, size %dx%d\n",x_res,y_res);
  else
    sprintf(str1,"Configuration with %d atoms\n",natoms);
  if (col_mode) {
    if (eng_mode)
      sprintf(str2,"Kinetic energy is encoded");
    else
      sprintf(str2,"Potential energy is encoded");
  }
  else
    sprintf(str2,"Atom type is encoded");
  if (bond_mode==1)
      sprintf(str3,"Static bonds");
  if (bond_mode==1)
      sprintf(str3,"Dynamic bonds");
  else
    sprintf(str3,"");

  move(0.0,0.7);
  drawstr(str1);  
  move(0.0,0.6);
  drawstr(str2);  
  move(0.0,0.5);
  drawstr(str3);  
  printf("%s\n%s\n%s\n", str1,str2,str3);fflush(stdout);
}

void draw_scene_wrapper(int scene_type) {
  char ch,mkey;
  float xloc,yloc;

  do {
    if (movie_mode)
      connect_client(9);
    draw_scene(scene_type);
    if (ch = checkkey()) { movie_mode=0; break;}
    if (mkey = slocator(&xloc, &yloc)) { movie_mode=0;break;}
  }
  while(movie_mode);
}

/* draw_scene - what a name */
void draw_scene(int scene_type) {

  int i,j,k,l,cv,nb,vgl;
  int n;
  float points[8][3];
  double xx,yy,xxj,yyj;
  int ixx,iyy;
  float epsilon,dx,dy,dr;
#ifndef TWOD
  double zz,zzj,dz;
#endif

  backbuffer();
  polyfill(1);
  color(BLACK);
  clear();

  /* scene_type determines whether we deal with a distr. or a conf.
     distribution first */
  if (scene_type==1) {
    color(RED);
    for (i=0;i<x_res*y_res;i++) {
      iyy=i%x_res;
      ixx=(i-iyy) / y_res;
      xx=(float)(ixx)*scalex-1;
      yy=(float)(iyy)*scaley-1;
      if (col_mode==2)
	cv=(int)floor(scalekin*(kinarray[i]-offskin));
      if (col_mode==3)
	cv=(int)floor(scalepot*(potarray[i]-offspot));
      color(cv+10);
      rect(xx,yy,xx+scalex,yy+scaley);
    } /* for */
  } /* distribution (scene_type==1) */

  /* now configuration */ 
  if ((scene_type==0)&&(atom_mode==1)) {
    for (i=0;i<natoms;i++) {
      xx=(x[i]-minx)*scalex-1.0;
      yy=(y[i]-miny)*scaley-1.0;
#ifndef TWOD
      zz=(z[i]-minz)*scalez*-1.0;
      zz*=.2; /* looks better */
#endif
      switch(col_mode) {
      case 0:
	color(1);
	break;
      case 1: 
	color(sorte[i]+1);
	break;
      case 2: 
	cv=(int)(scalepot*(pot[i]-offspot));
	color(cv+10);
	break;
      case 3: 
	cv=(int)(scalekin*(kin[i]-offskin));
	color(cv+10);
	break;
      case 4:
	color(nb+1);
	break;
      case 5:
	color(nummer[i]<0?YELLOW:BLUE);
	break;
      }

      /* if we shall draw bonds afterwards we compute
	 the number of neighbours here */
      if (bond_mode==1) {
	nb=0;
	for(k=0;k<nunits;k++)
	  if (bcode[i]&(int)pow(2,k))
	    nb++;

	/*	color(nb+1);*/
      }

      if (size_mode)
	circle(xx,yy,.5*(sorte[i]+1)*radius*scalex);
      else {
#ifndef TWOD
	/* draw a cube */
	makepoly();
	move(xx-.01,yy-.01,zz-.01);
	draw(xx-.01,yy+.01,zz-.01);
	draw(xx+.01,yy+.01,zz-.01);
	draw(xx+.01,yy-.01,zz-.01);
	draw(xx+.01,yy-.01,zz+.01);
	draw(xx-.01,yy-.01,zz+.01);
	draw(xx-.01,yy-.01,zz-.01);
	draw(xx+.01,yy-.01,zz-.01);
	closepoly();
	makepoly();
	move(xx-.01,yy+.01,zz+.01);
	draw(xx+.01,yy+.01,zz+.01);
	draw(xx+.01,yy+.01,zz-.01);
	draw(xx-.01,yy+.01,zz-.01);
	draw(xx-.01,yy+.01,zz+.01);
	draw(xx-.01,yy-.01,zz+.01);
	draw(xx+.01,yy-.01,zz+.01);
	draw(xx+.01,yy+.01,zz+.01);
	closepoly();
#else
	circle(xx,yy,.01);
#endif
      } /* radius encoding */
    } /* for i= */
  } /* configuration */

  if (bond_mode) draw_bonds();
  if (text) draw_text();
  swapbuffers();
}

void draw_bonds() {
  int i,j,k;
  float xx,yy,xxj,yyj,dx,dy,dr,epsilon;
#ifndef TWOD
  float zz,zzj,dz;
#endif
  epsilon=.1;

  color(MAGENTA);
  for (i=0;i<natoms;i++) {
    xx=(x[i]-minx)*scalex-1;
    yy=(y[i]-miny)*scaley-1;
#ifndef TWOD
    zz=(z[i]-minz)*scalez-1;
    zz*=.2; /* looks better */
#endif
    for (j=0;j<natoms;j++) {
      xxj=(x[j]-minx)*scalex-1;
      yyj=(y[j]-miny)*scaley-1;
#ifndef TWOD
      zzj=(z[j]-minz)*scalez-1;
      zzj*=.2; /* looks better */
#endif
      if (i==j) continue;
      if (bond_mode == 1) {
	if (qp) {
#ifdef TWOD
	  if (sorte[i]!=sorte[j]) continue;
#else
	  if (sorte[i]!=0) continue;
	  if (sorte[j]!=0) continue;
#endif
	}
	dx=x[i]-x[j];
	if (ABS(dx)>2.4) continue;
	dy=y[i]-y[j];
	if (ABS(dy)>2.4) continue;
#ifndef TWOD
	dz=z[i]-z[j];
	if (ABS(dz)>2.4) continue;	    
#endif
	for (k=0;k<nunits;k++) {
	  if (ABS(dx-ux[k])<epsilon)
	    if (ABS(dy-uy[k])<epsilon) {
#ifndef TWOD
	      if (ABS(dz-uz[k])<epsilon) {
#endif
		bcode[i]+=pow(2,k);
#ifndef TWOD
		color(MAGENTA);
		if (k<2) color(YELLOW);
	        if (k>3) color(CYAN);
		move(xx,yy,zz);
		draw(xxj,yyj,zzj);
	      }
#else
	      move2(xx,yy);
	      draw2(xxj,yyj);
#endif
	    }
	}
      }
      if (bond_mode==2) {
	if (qp) {
#ifndef TWOD
	  if (sorte[i]!=0) continue;
	  if (sorte[j]!=0) continue;
#else
	  if (sorte[i]==sorte[j]) continue;
#endif
	  dx=x[i]-x[j];
	  dy=y[i]-y[j];
#ifndef TWOD
	  dz=z[i]-z[j];
	  dr=dx*dx+dy*dy+dz*dz;
#else
	  dr=dx*dx+dy*dy;
#endif
	  if (ABS(dr)<.9) continue;
	  if (ABS(dr)>1.1) continue;
#ifndef TWOD
	  move(xx,yy,zz);
	  draw(xxj,yyj,zzj);
#else
	  move2(xx,yy);
	  draw2(xxj,yyj);
#endif
	  bcode[i]++;
	}
	else {
	  dx=x[i]-x[j];
	  dy=y[i]-y[j];
#ifndef TWOD
	  dz=z[i]-z[j];
	  dr=dx*dx+dy*dy+dz*dz;
#else
	  dr=dx*dx+dy*dy;
#endif
	  if (ABS(dr)>1.1) continue;
	  if (ABS(dr)<.9) continue;
#ifndef TWOD
	  move(xx,yy,zz);
	  draw(xxj,yyj,zzj);
#else
	  move2(xx,yy);
	  draw2(xxj,yyj);
#endif
	  bcode[i]++;
	}
      }
    } /* for j... */
  } /* for i... */

}

/* reading of a configuration from file */
int read_distribution(char *fname) {
}

/* reading of a configuration from file */
int read_configuration(char *fname) {

  FILE *fp;
  char line[200],str[255];
  int n;

  fp=fopen(fname, "r");
  if (fp==NULL) {
    printf("Datei %s existiert nicht\n",fname);
    return -1;
  }
  n=0;

  while(fgets(line,200,fp)) {
    n++;
  }

  fclose(fp);
  fp=fopen(fname, "r");

  nummer=(int*)calloc(n,sizeof(int));
  sorte=(short int*)calloc(n,sizeof(int));
  masse=(double*)calloc(n,sizeof(double));
  x=(double*)calloc(n,sizeof(double));
  y=(double*)calloc(n,sizeof(double));
#ifndef TWOD
  z=(double*)calloc(n,sizeof(double));
  vz=(double*)calloc(n,sizeof(double));
#endif
  vx=(double*)calloc(n,sizeof(double));
  vy=(double*)calloc(n,sizeof(double));
  pot=(double*)calloc(n,sizeof(double));
  kin=(double*)calloc(n,sizeof(double));
  bcode=(int*)calloc(n,sizeof(int));

  n=0;
  while (fgets(line,200,fp)) {
#ifdef TWOD
    columns=sscanf(line,"%d %d %lf %lf %lf %lf %lf %lf",&nummer[n],&sorte[n],&masse[n],&x[n],&y[n],&vx[n],&vy[n],&pot[n]);
    if (columns==8) {
      kin[n] = vx[n]*vx[n]+vy[n]*vy[n];
#else
    columns=sscanf(line,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf",&nummer[n],&sorte[n],&masse[n],&x[n],&y[n],&z[n],&vx[n],&vy[n],&vz[n],&pot[n]);
    if (columns==10) {
      kin[n] = vx[n]*vx[n]+vy[n]*vy[n]+vz[n]*vz[n];
#endif
      if (maxp<pot[n]) maxp=pot[n];
      if (minp>pot[n]) minp=pot[n];
      if (maxk<kin[n]) maxk=kin[n];
      if (mink>kin[n]) mink=kin[n];
    }
    if (maxx<x[n]) maxx=x[n];
    if (maxy<y[n]) maxy=y[n];
    if (minx>x[n]) minx=x[n];
    if (miny>y[n]) miny=y[n];
#ifndef TWOD
    if (maxz<z[n]) maxz=z[n];
    if (minz>z[n]) minz=z[n];
#endif
    n++;
  }

  if (maxx==minx) {
    minx=0;
    scalex=1.0;
  }
  else
    scalex=2.0/(maxx-minx);
  if (maxy==miny) {
    minx=0; 
    scaley=1.0;
  }
  else
    scaley=2.0/(maxy-miny);
#ifndef TWOD
  if (maxz==minz) {
    minx=0;
    scalez=1.0;
  }
  else
    scalez=2.0/(maxz-minz);
#endif
  if (maxp==minp)
    scalepot=1.0;
  else
    scalepot=COLRES/(maxp-minp);
  if (maxk==mink)
    scalekin=1.0;
  else
    scalekin=COLRES/(maxk-mink);
  offspot=minp;
  offskin=mink;
  fclose(fp);

  return n;
}

  /* reading of a configuration from file */
int write_configuration(char *fname) {

  FILE *fp;
  char str[255];
  int n;

  fp=fopen(fname, "w");
  if (fp==NULL) {
    printf("Kann Datei %s nicht anlegen.\n",fname);
    return -1;
  }

  for (n=0;n<=natoms;n++)
#ifdef TWOD
    fprintf(fp,"%d %d %lf %lf %lf %lf %lf %lf\n",nummer[n],sorte[n],masse[n],x[n],y[n],vx[n],vy[n],pot[n]);
#else
    fprintf(fp,"%d %d %lf %lf %lf %lf %lf %lf\n",nummer[n],sorte[n],masse[n],x[n],y[n],z[n],vx[n],vy[n],vz[n],pot[n]);
#endif

  fclose(fp);
  return n;
}

int write_distribution(char *fname) {
  FILE *fp;
  int i;

  fp=fopen("2dvis.dist","w");
  fprintf(fp,"%d %d %s\n", x_res,y_res,(eng_mode)?"kin":"pot");
  for (i=0;i<x_res*y_res;i++)
    fprintf(fp,"%f %f\n",potarray[i],kinarray[i]);

  fclose(fp);
}

void read_unit_vectors(void) {
  FILE *fp;
  char line[255];
  int i;

  fp=fopen(uvfname,"r");

  fgets(line,200,fp);
  if (strncmp(line,"period",6))
    qp=0;
  if (strncmp(line,"quasip",6))
    qp=1;

  fgets(line,200,fp);
  sscanf(line,"%d",&nunits);
  if (ux==NULL)
    ux=(float*)calloc(nunits,sizeof(float));
  if (uy==NULL)
    uy=(float*)calloc(nunits,sizeof(float));
#ifndef TWOD
  if (uz==NULL)
    uz=(float*)calloc(nunits,sizeof(float));
#endif

  i=0;
  while(fgets(line,200,fp)) {
#ifdef TWOD
    sscanf(line,"%f %f\n",&ux[i],&uy[i]);
#else
    sscanf(line,"%f %f %f\n",&ux[i],&uy[i],&uz[i]);
#endif
    i++;
  }

  fclose(fp);
}

/* sos */
void display_help(void) {

  char sysstring[1000];

  sprintf(sysstring,"Use the following keys\n\nLMB:\tmove configuration\nMMB:\trotate configuration around x-axis\nR+MMBs\trotate configuration around y-axis\nRMB:\tscale configuration\n1:\tmove configuration to the lower left\n2:\tmove configuration down\n3:\tmove configuration to the lower right\n4:\tmove configuration to the left\n6:\tmove configuration to the right\n7:\tmove configuration to the upper left\n8:\tmove configuration up\n9:\tmove configuration to the upper right\n0:\tincrease radius\n.:\tdecrease radius\n+:\tscale configuration (larger)\n-:\tscale configuration (smaller)\n*:\trotate about x-axis\n/:\trotate about y-axis\na:\ttoggle drawing of atoms\nb:\ttoggle drawing of bonds\nc:\tget configuration via socket\nd:\tget distribution via socket\ne:\ttoggle type/energy encoding\nh:\tdisplay this help message\nk:\ttoggle potential/kinetic energy\nl:\tload IMD-configuration from file\nm:\tmovie mode\np\tmake picture\nq:\tquit program\nr:\ttoggle type encoding via radius\nt:\ttoggle text on/off\nw:\twrite to file\n");
  printf(sysstring);
  fflush(stdout);
} 

void usage(void) {
}

void error(char *mesg)

{
  
#ifdef MPI
  fprintf(stderr,"CPU %d: %s\n",myid,mesg);
  MPI_Abort(MPI_COMM_WORLD, 1);
#else
  fprintf(stderr,"Error: %s\n",mesg);
#endif
  exit(2);
}
