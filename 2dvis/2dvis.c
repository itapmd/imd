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
  float xloc,yloc,xlocold,ylocold,delta;
  char ch;
  int mkey;
  char fname[255];

  /* inits */
  putenv("USEOWNCMAP=1");
  strcpy(uvfname,"unitvecs");

  read_parameters(argc,argv);

  init_graph();

  /* main loop */
  do {
    /* Taste, oder Mauskey erlaubt */
    while (!((ch = checkkey()) || (mkey = slocator (&xloc, &yloc))));
    if (mkey==1) {
      xlocold=xloc;
      ylocold=yloc;
      while ((mkey=slocator(&xloc,&yloc))==1)
	;
      translate(xloc-xlocold,yloc-ylocold,.000001);
      draw_scene(scene_type);
    }
    if (mkey==2) {
      xlocold=xloc;
      ylocold=yloc;
      while ((mkey=slocator(&xloc,&yloc))==2)
	;
      rotate(100*(xloc-xlocold),'y');
      rotate(100*(ylocold-yloc),'x');
      draw_scene(scene_type);
    }
    if (mkey==3) { 
      xlocold=xloc;
      ylocold=yloc;
      while ((mkey=slocator(&xloc,&yloc))==3)
	;
      delta=sqrt((xloc-xlocold)*(xloc-xlocold)+(yloc-ylocold)*(yloc-ylocold));
      if (atan2(yloc-ylocold,xloc-xlocold)>0)
	scale(1+.1*delta,1+.1*delta,0);
      else
	scale(1-.1*delta,1-.1*delta,0);
      draw_scene(scene_type);
    }
    switch (ch) {
    case '0' : 
      radius*=1.5;
      draw_scene(scene_type);
      break;
    case '.' : 
      radius/=1.5;
      draw_scene(scene_type);
      break;
    case '1' : 
      translate(-.1,-.1,.000001);
      draw_scene(scene_type);
      break;
    case '2' : 
      translate(.000001,-.1,.000001);
      draw_scene(scene_type);
      break;
    case '3' : 
      translate(.1,-.1,.000001);
      draw_scene(scene_type);
      break;
    case '4' : 
      translate(-.1,.000001,.000001);
      draw_scene(scene_type);
      break;
    case '6' : 
      translate(.1,.000001,.000001);
      draw_scene(scene_type);
      break;
    case '7' : 
      translate(-.1,.1,.000001);
      draw_scene(scene_type);
      break;
    case '8' : 
      translate(.000001,.1,.000001);
      draw_scene(scene_type);
      break;
    case '9' : 
      translate(.1,.1,.000001);
      draw_scene(scene_type);
      break;
    case '+' : 
      scale(1.1,1.1,1.1);
      draw_scene(scene_type);
      break;
    case '-' : 
      scale(.9,.9,.9);
      draw_scene(scene_type);
      break;
    case '*' : 
      rotate(.1,'x');
      draw_scene(scene_type);
      break;
    case '/' : 
      rotate(.1,'y');
      draw_scene(scene_type);
      break;
    case 'a' :
      if (atom_mode)
	atom_mode=0;
      else {
	atom_mode=1;
      }
      draw_scene(scene_type);
      break;
    case 'b' : 
      if (bond_mode)
	bond_mode=0;
      else {
	bond_mode=1;
      }
      draw_scene(scene_type);
      break;
    case 'c' : 
      scene_type=0;
      if (connect_client(9)==0) {
	printf("No atoms!\n");
	exit(-1);
      }
      draw_scene(scene_type);
      break;
    case 'd' : 
      scene_type=1;
      if (connect_client(1)==0) {
	printf("No distribution!\n");
	exit(-1);
      }
      draw_scene(scene_type);
      break;
    case 'e' : 
      if (col_mode) 
	col_mode=0; 
      else
	col_mode=1;
      draw_scene(scene_type);
      break;
    case 'h' : 
      display_help();
      break;
    case 'k' : 
      if (eng_mode)
	eng_mode=0;
      else 
	eng_mode=1;
      draw_scene(scene_type);
      break;
    case 'l' : 
      do { 
	printf("Enter Filename: ");
	scanf("%s", fname);
      } 
      while ((natoms=read_atoms(fname))<0);
      if (natoms==0) {
	printf("Kein Atom gelesen\n");
	exit(-1);
      }
      draw_scene(scene_type);
      break;
    case 'm' : 
      while (!((ch = checkkey()) || (mkey = slocator (&xloc, &yloc)))) {
	if (connect_client(9)==0) {
	  printf("No atoms!\n");
	  exit(-1);
	}
	draw_scene(scene_type);
      }
      break;
    case 'p' : 
      make_picture();
      break;
    case 'q' : 
      vexit();
      exit(0);
      break;
    case 'r' : 
      if (radectyp) 
	radectyp=0; 
      else 
	radectyp=1;
      draw_scene(scene_type);
      break;
    case 's':
      if (stat_bond)
	stat_bond=0;
      else
	stat_bond=1;
      read_unit_vectors();
      break;
    case 't' : 
      if (text) 
	text=0; 
      else 
	text=1;
      draw_scene(scene_type);
      break;
    case 'w' : 
      write_to_file();
      break;
    default: break;
    }
  }
  while (1);
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

  for (i=0;i<COLRES;i++)
    mapcolor(i+10,i,0,245-i);



}

/* draw_scene - what a name */
void draw_scene(int scene_type) {

  int i,j,k,cv,nb;
  int n;
  float points[8][3];
  double xx,yy,xxj,yyj;
  int ixx,iyy;
  char str[200];
  float epsilon,dx,dy,dr;
#ifndef TWOD
  double zz,minz;
#endif

  epsilon=.1;
  backbuffer();
  color(BLACK);
  clear();
  polyfill(1);

  /* scene_type determines whether we deal with a distr. or a conf. */
  if (scene_type==1) {
    color(RED);
    for (i=0;i<x_res*y_res;i++) {
      iyy=i%x_res;
      ixx=(i-iyy) / y_res;
      xx=(float)(ixx)*scalex-1;
      yy=(float)(iyy)*scaley-1;
      if (eng_mode)
	cv=(int)floor(scalekin*(kinarray[i]-offskin));
      else
	cv=(int)floor(scalepot*(potarray[i]-offspot));
      color(cv+10);
      rect(xx,yy,xx+scalex,yy+scaley);
    } /* for */

    /* drawing of text */
    if (text) {
      color(CYAN);
      move(0.0,0.7);
      sprintf(str,"Distribution, size %dx%d\n",x_res,y_res);
      drawstr(str);
      if (eng_mode)
	sprintf(str,"Kinetic energy is encoded");
      else
	sprintf(str,"Potential energy is encoded");
      drawstr(str);
    } /* text */
  } /* distribution (scene_type==1) */



  if (scene_type==0) {
    if (bond_mode==1) {
      if (text) {
	color(CYAN);
	move(0.0,0.7);
	sprintf(str,"Drawing bonds");
	drawstr(str);
      }
      for (i=0;i<natoms;i++) {
	xx=(x[i]-minx)*scalex-1;
	yy=(y[i]-miny)*scaley-1;
#ifndef TWOD
	zz=z[i]*scalez-1;
#endif

	bcode[i]=0;
      }
      color(MAGENTA);
      for (i=0;i<natoms;i++) {
	xx=(x[i]-minx)*scalex-1;
	yy=(y[i]-miny)*scaley-1;
	for (j=0;j<natoms;j++) {
	  xxj=(x[j]-minx)*scalex-1;
	  yyj=(y[j]-miny)*scaley-1;
	  if (i==j) continue;
	  if (stat_bond == 1) {
	    dx=x[i]-x[j];
	    if (ABS(dx)>1.4) continue;
	    dy=y[i]-y[j];
	    if (ABS(dy)>1.4) continue;
	    for (k=0;k<nunits;k++) {
	      if (ABS(dx-ux[k])<epsilon)
		if (ABS(dy-uy[k])<epsilon) {
		  move2(xx,yy);
		  draw2(xxj,yyj);
		  bcode[i]+=pow(2,k);
		}
	    }
	  }
	  else { /* stat_bond==0 */
	    if (qp) {
	      if (sorte[i]==sorte[j]) continue;
	      dx=x[i]-x[j];
	      dy=y[i]-y[j];
	      dr=dx*dx+dy*dy;
	      if (ABS(dr)>1.1) continue;
	      if (ABS(dr)<0.9) continue;
	      move2(xx,yy);
	      draw2(xxj,yyj);
	      bcode[i]++;
	    }
	    else {
	      dx=x[i]-x[j];
	      dy=y[i]-y[j];
	      dr=dx*dx+dy*dy;
	      if (ABS(dr)>1.1) continue;
	      if (ABS(dr)<0.9) continue;
	      move2(xx,yy);
	      draw2(xxj,yyj);
	      bcode[i]++;
	    }
	  }
	}
      } /* for i... */
    }
     
    if (atom_mode==1) {
      for (i=0;i<natoms;i++) {
	xx=(x[i]-minx)*scalex-1;
	yy=(y[i]-miny)*scaley-1;
	nb=0;
	for(k=0;k<nunits;k++)
	  if (bcode[i]&(int)pow(2,k)) {
	    nb++;
	  }
	if (col_mode==0) {
	  if (sorte[i]==0) color(RED);
	  if (sorte[i]==1) color(GREEN);
	  if (sorte[i]==2) color(MAGENTA);
	  if (sorte[i]==3) color(WHITE);
	  if (sorte[i]==4) color(BLUE);
	  if (sorte[i]==5) color(YELLOW);
	  if (sorte[i]==6) color(CYAN);
	}
	else {
	  if (bond_mode==1) {
	    color(WHITE);
	    if (nb==0) color(RED);
	    if (nb==1) color(BLUE);
	    if (nb==2) color(GREEN);
	    if (nb==3) color(YELLOW);
	    if (nb==4) color(MAGENTA);
	    if (nb==5) color(WHITE);
	    if (nb==6) color(CYAN);
	  }
	  else { /* bond_mode == 0 */
	    if (eng_mode)
	      cv=(int)(scalekin*(kin[i]+offskin));
	    else
	      cv=(int)(scalepot*(pot[i]+offspot));
	    mapcolor(i+8,cv,cv,cv);
	    color(i+8);
	  }
        }
	
	if (radectyp)
	  circle(xx,yy,.5*(sorte[i]+1)*radius*scalex);
	else {
#ifndef TWOD
	  /*
	    move(xx,yy,zz);
	    draw(xx,yy,zz+.1);
	    draw(xx+.01,yy,zz+.1);
	    draw(xx+.005,yy+.00866,zz+.1);
	    draw(xx,yy,zz+.1);
	    draw(xx+.005,yy+.00297,zz+.10866);
	    draw(xx+.01,yy,zz+.1);
	    draw(xx+.005,yy+.00866,zz+.1);
	    draw(xx+.005,yy+.00297,zz+.10866);
	  */
	  move(xx,yy,zz);
	  rdraw(0,0,.1);
	  rdraw(.01,0,.1);
	  rdraw(.005,.00866,.1);
	  rdraw(0,0,.1);
	  rdraw(.005,.00297,.10866);
	  rdraw(.01,0,.1);
	  rdraw(.005,.00866,.1);
	  rdraw(.005,.00297,.10866);

#else
	  circle(xx,yy,radius*scalex);
#endif
	}
	if (text) {
	  color(CYAN);
	  move(0.0,0.7);
	  sprintf(str,"Configuration with %d atoms\n",natoms);
	  drawstr(str);
	  if (col_mode)
	    if (eng_mode)
	      sprintf(str,"Kinetic energy is encoded");
	    else
	      sprintf(str,"Potential energy is encoded");
	  else
	    sprintf(str,"Atom type is encoded");
	  drawstr(str);
	}
      } /* atoms? (char "a") */
    } /* bonds */
  } /* atoms (scene_type==0) */
  swapbuffers();
}

/* reading of a configuration from file */
int read_atoms(char *fname) {

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
  vx=(double*)calloc(n,sizeof(double));
  vy=(double*)calloc(n,sizeof(double));
  pot=(double*)calloc(n,sizeof(double));
  kin=(double*)calloc(n,sizeof(double));
  bcode=(int*)calloc(n,sizeof(int));

  n=0;
  while (fgets(line,200,fp)) {
#ifdef TWOD
    columns=sscanf(line,"%d %d %lf %lf %lf %lf %lf %lf",&nummer[n],&sorte[n],&masse[n],&x[n],&y[n],&vx[n],&vy[n],&pot[n]);
#else
    columns=sscanf(line,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf",&nummer[n],&sorte[n],&masse[n],&x[n],&y[n],&z[n],&vx[n],&vy[n],&vz[n],&pot[n]);
#endif
    kin[n] = vx[n]*vx[n]+vy[n]*vy[n];
    if (maxx<x[n]) maxx=x[n];
    if (maxy<y[n]) maxy=y[n];
    if (minx>x[n]) minx=x[n];
    if (miny>y[n]) miny=y[n];
    if (maxp<pot[n]) maxp=pot[n];
    if (minp>pot[n]) minp=pot[n];
    if (maxk<kin[n]) maxk=kin[n];
    if (mink>kin[n]) mink=kin[n];
    n++;
  }

  scalex=2.0/(maxx-minx);
  scaley=2.0/(maxy-miny);
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

void make_picture(void) {

  char str[255];

  sprintf(str,"import -window 0x4800002 2dvis.gif");
  system(str);
}

void write_to_file(void) {
  FILE *fp;
  int i;

  if (scene_type) {
    fp=fopen("2dvis.dist","w");
    fprintf(fp,"%d %d %s\n", x_res,y_res,(eng_mode)?"kin":"pot");
    for (i=0;i<x_res*y_res;i++)
      fprintf(fp,"%f %f\n",potarray[i],kinarray[i]);
  }
  else {
    fp=fopen("2dvis.id","w");
    for (i=0;i<natoms;i++) {
      if (columns==6)
	fprintf(fp,"%d %d %f %f %f\n",nummer[i],sorte[i],masse[i],x[i],y[i]);
      else
	fprintf(fp,"%d %d %f %f %f %f %f %f\n",nummer[i],sorte[i],masse[i],x[i],y[i],vx[i],vy[i],pot[i]);
    }
  }
  fclose(fp);
}

void read_unit_vectors(void) {
  FILE *fp;
  char line[255];
  int i;

  fp=fopen(uvfname,"r");

  fgets(line,200,fp);
  if (strncmp(line,"periodic",8))
    qp=0;
  if (strncmp(line,"quasiper",8))
    qp=1;

  fgets(line,200,fp);
  sscanf(line,"%d",&nunits);
  ux=(float*)calloc(nunits,sizeof(float));
  uy=(float*)calloc(nunits,sizeof(float));

  i=0;
  while(fgets(line,200,fp)) {
    sscanf(line,"%f %f\n",&ux[i],&uy[i]);
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
