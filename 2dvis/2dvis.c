#include <stdio.h>
#include <vogle.h>
#include <math.h>

/* Abs berechnet den Betrag einer Zahl */
#define ABS(a) ((a) >0 ? (a) : -(a))

#define TWOD
#define MAXNATOMS 12000
#define COLRES 245

int *nummer,*sorte,columns,*bcode;
double *masse,*x,*y,*z,*vx,*vy,*vz,*pot,*kin;
float *potarray,*kinarray;
int x_res,y_res;
int natoms,nunits;
unsigned short base_port;
float scalex,scaley,scalepot,scalekin,radius,offspot,offskin;
int colmode,scene_type,text,engmode,qp,radectyp;
float *ux,*uy;

void draw_scene(int scene_type);
void draw_bonds(void);
void init_graph(void);
void display_help(void);
int read_atoms(char *fname);
void read_unit_vectors(void);
void write_to_file(void);

main()
{
  float xloc,yloc,xlocold,ylocold,delta;
  char ch;
  int mkey;
  char fname[255];

  /* inits */
  base_port=31913;
  colmode=0;
  scene_type=0;
  text=0;
  radius=.3;
  radectyp=0;

  /* allocs */
  nummer=(int*)calloc(MAXNATOMS,sizeof(int));
  sorte=(int*)calloc(MAXNATOMS,sizeof(int));
  masse=(double*)calloc(MAXNATOMS,sizeof(float));
  x=(double*)calloc(MAXNATOMS,sizeof(float));
  y=(double*)calloc(MAXNATOMS,sizeof(float));
  vx=(double*)calloc(MAXNATOMS,sizeof(float));
  vy=(double*)calloc(MAXNATOMS,sizeof(float));
  pot=(double*)calloc(MAXNATOMS,sizeof(float));
  kin=(double*)calloc(MAXNATOMS,sizeof(float));
  bcode=(int*)calloc(MAXNATOMS,sizeof(int));

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
      delta=sqrt((xloc-xlocold)*(xloc-xlocold)+(yloc-ylocold)*(yloc-ylocold));
      rotate(100*delta,'x');
      draw_scene(scene_type);
    }
    if (mkey==4) { 
      xlocold=xloc;
      ylocold=yloc;
      while ((mkey=slocator(&xloc,&yloc))==4)
	;
      delta=sqrt((xloc-xlocold)*(xloc-xlocold)+(yloc-ylocold)*(yloc-ylocold));
      rotate(100*delta,'y');
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
    case '0' : radius*=1.5;draw_scene(scene_type);break;
    case '.' : radius/=1.5;draw_scene(scene_type);break;
    case '1' : translate(-.1,-.1,.000001);draw_scene(scene_type);break;
    case '2' : translate(.000001,-.1,.000001);draw_scene(scene_type);break;
    case '3' : translate(.1,-.1,.000001);draw_scene(scene_type);break;
    case '4' : translate(-.1,.000001,.000001);draw_scene(scene_type);break;
    case '6' : translate(.1,.000001,.000001);draw_scene(scene_type);break;
    case '7' : translate(-.1,.1,.000001);draw_scene(scene_type);break;
    case '8' : translate(.000001,.1,.000001);draw_scene(scene_type);break;
    case '9' : translate(.1,.1,.000001);draw_scene(scene_type);break;
    case '+' : scale(1.1,1.1,1.1);draw_scene(scene_type);break;
    case '-' : scale(.9,.9,.9);draw_scene(scene_type);break;
    case '*' : rotate(.1,'x');draw_scene(scene_type);break;
    case '/' : rotate(.1,'y');draw_scene(scene_type);break;
    case 'b' : scene_type=2;draw_bonds();draw_scene(scene_type);break;
    case 'c' : scene_type=0;if (connect_client(9)==0) {printf("No atoms!\n");exit(-1);};draw_scene(scene_type);break;
    case 'd' : scene_type=1;if (connect_client(1)==0) {printf("No distribution!\n");exit(-1);}draw_scene(scene_type);break;
    case 'e' : if (colmode) colmode=0; else colmode=1;draw_scene(scene_type);break;
    case 'h' : display_help();break;
    case 'k' : if (engmode) engmode=0; else engmode=1;draw_scene(scene_type);break;
    case 'l' : do { printf("Enter Filename: ");scanf("%s", fname);} while ((natoms=read_atoms(fname))<0);if (natoms==0) { printf("Kein Atom gelesen\n");exit(-1); };draw_scene(scene_type);break;
    case 'p' : make_picture();break;
    case 'r' : if (radectyp) radectyp=0; else radectyp=1;draw_scene(scene_type);break;
    case 'q' : vexit(); exit(0);break;
    case 't' : if (text) text=0; else text=1;draw_scene(scene_type);break;
    case 'u' : read_unit_vectors();break;
    case 'w' : write_to_file();break;
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

  /* fill the colormap */
  for (i=10;i<COLRES+10;i++)
    mapcolor(i,i,i,i);

}

/* draw_scene - what a name */
void draw_scene(int scene_type) {

  int i,cv,nb,k;
  double xx,yy;
  int ixx,iyy;
  char str[200];
  color(BLACK);clear();

  /* scene_type determines whether we deal with a distr. or a conf. */
  if (scene_type==1) {
    color(RED);
    for (i=0;i<x_res*y_res;i++) {
      iyy=i%x_res;
      ixx=(i-iyy) / y_res;
      xx=(float)(ixx)*scalex-1;
      yy=(float)(iyy)*scaley-1;
      if (engmode)
	cv=(int)floor(scalekin*(kinarray[i]-offskin));
      else
	cv=(int)floor(scalepot*(potarray[i]-offspot));
      color(cv+10);
      rect(xx,yy,xx+scalex,yy+scaley);
    }
    /* drawing of text */
    if (text) {
      color(CYAN);
      move(0.0,0.7);
      sprintf(str,"Distribution, size %dx%d\n",x_res,y_res);
      drawstr(str);
      if (engmode)
	sprintf(str,"Kinetic energy is encoded");
      else
	sprintf(str,"Potential energy is encoded");
      drawstr(str);
    }
  }
  if (scene_type==0) {
    for (i=0;i<natoms;i++) {
      if (colmode==0) {
	if (sorte[i]==0) color(RED);
	if (sorte[i]==1) color(GREEN);
	if (sorte[i]==2) color(MAGENTA);
	if (sorte[i]==3) color(WHITE);
	if (sorte[i]==4) color(BLUE);
	if (sorte[i]==5) color(YELLOW);
	if (sorte[i]==6) color(CYAN);
      }
      else {
	if (engmode)
	  cv=(int)(scalekin*(kin[i]+offskin));
	else
	  cv=(int)(scalepot*(pot[i]+offspot));
	mapcolor(i+8,cv,cv,cv);
	color(i+8);
      }
      xx=x[i]*scalex-1;
      yy=y[i]*scaley-1;
      if (sorte[i]>10) sorte[i]=0;
      printf("%d %d %d %f %f %f\n", i, nummer[i], sorte[i], masse[i], x[i], y[i]);fflush(stdout);
      if (radectyp)
	circle(xx,yy,radius*scalex);
      else
	circle(xx,yy,.5*(sorte[i]+1.0)*radius*scalex);
    }
    if (text) {
      color(CYAN);
      move(0.0,0.7);
      sprintf(str,"Configuration with %d atoms\n",natoms);
      drawstr(str);
      if (colmode)
	if (engmode)
	  sprintf(str,"Kinetic energy is encoded");
	else
	  sprintf(str,"Potential energy is encoded");
      else
	sprintf(str,"Atom type is encoded");
      drawstr(str);
    }
  }

  if (scene_type==2) {
    for (i=0;i<natoms;i++) {
      nb=0;
      xx=x[i]*scalex-1;
      yy=y[i]*scaley-1;
      color(CYAN);
      for(k=0;k<nunits;k++)
	if (bcode[i]&(int)pow(2,k)) {
	  move2(xx,yy);
	  draw2(xx-ux[k]*.5*scalex,yy-uy[k]*.5*scaley);
	  nb++;
	}
      if (colmode==0) {
	if (sorte[i]==0) color(RED);
	if (sorte[i]==1) color(GREEN);
	if (sorte[i]==2) color(MAGENTA);
	if (sorte[i]==3) color(WHITE);
	if (sorte[i]==4) color(BLUE);
	if (sorte[i]==5) color(YELLOW);
	if (sorte[i]==6) color(CYAN);
      }
      else {
	color(WHITE);
	if (nb==0) color(RED);
	if (nb==1) color(BLUE);
	if (nb==2) color(GREEN);
	if (nb==3) color(YELLOW);
	if (nb==4) color(MAGENTA);
	if (nb==5) color(WHITE);
	if (nb==6) color(CYAN);
      }
      if (sorte[i]>10) sorte[i]==0;
      printf("%d %f\n", i, sorte[i]);fflush(stdout);
      if (radectyp)
	circle(xx,yy,radius*scalex);
      else
	circle(xx,yy,.5*(sorte[i]+1)*radius*scalex);
    }
    if (text) {
      color(CYAN);
      move(0.0,0.7);
      sprintf(str,"Configuration with %d atoms\n",natoms);
      drawstr(str);
      if (colmode)
	sprintf(str,"Potential energy is encoded");
      else
	sprintf(str,"Atom type is encoded");
      drawstr(str);
    }
  }
}

/* reading of a configuration from file */
int read_atoms(char *fname) {

  FILE *fp;
  char line[200],str[255];
  int n;
  float maxx,minx,maxy,miny,maxp,minp,maxk,mink;

  maxx=-1000;
  maxy=-1000;
  minx=1000;
  miny=1000;
  maxp=-1000;
  maxk=-1000;
  minp=1000;
  mink=1000;

#ifndef TWOD
  z=(double*)calloc(MAXNATOMS,sizeof(float));
#endif
  pot=(float*)calloc(MAXNATOMS,sizeof(float));

  fp=fopen(fname, "r");
  if (fp==NULL) {
    printf("Datei %s existiert nicht\n",fname);
    return -1;
  }
  n=0;
  while (fgets(line,200,fp)) {
#ifdef TWOD
    columns=sscanf(line,"%d %d %lf %lf %lf %lf %lf %lf",&nummer[n],&sorte[n],&masse[n],&x[n],&y[n],&vx[n],&vy[n],&pot[n]);
#else
    columns=sscanf(line,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf",&nummer[n],&sorte[n],&masse[n],&x[n],&y[n],&z[n],&vx[n],&vy[n],&vz[n],&pot[n]);
#endif

    if (maxx<x[n]) maxx=x[n];
    if (maxy<y[n]) maxy=y[n];
    if (minx>x[n]) minx=x[n];
    if (miny>y[n]) miny=y[n];
    if (maxp<pot[n]) minp=pot[n];
    if (minp>pot[n]) minp=pot[n];
    n++;
  }

  scalex=2.0/(maxx-minx);
  scaley=2.0/(maxy-miny);
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
    fprintf(fp,"%d %d %s\n", x_res,y_res,(engmode)?"kin":"pot");
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

  fp=fopen("unitvecs","r");

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

void draw_bonds(void) {

  float epsilon=.1;
  int i,j,k;
  float dx,dy;

  for (i=0;i<natoms;i++)
    bcode[i]=0;

  for (i=0;i<natoms;i++) {
    for (j=0;j<natoms;j++) {
      if (i==j) continue;
      dx=x[i]-x[j];
      if (ABS(dx)>1.4) continue;
      dy=y[i]-y[j];
      if (ABS(dy)>1.4) continue;
      for (k=0;k<nunits;k++) {
	if (ABS(dx-ux[k])<epsilon)
	  if (ABS(dy-uy[k])<epsilon)
	    bcode[i]+=pow(2,k);
      }
    }
  }
}

/* sos */
void display_help(void) {

  char sysstring[1000];

  sprintf(sysstring,"Use the following keys\n\nLMB:\tmove configuration\nMMB:\trotate configuration around x-axis\nR+MMBs\trotate configuration around y-axis\nRMB:\tscale configuration\n1:\tmove configuration to the lower left\n2:\tmove configuration down\n3:\tmove configuration to the lower right\n4:\tmove configuration to the left\n6:\tmove configuration to the right\n7:\tmove configuration to the upper left\n8:\tmove configuration up\n9:\tmove configuration to the upper right\n0:\tincrease radius\n.:\tdecrease radius\n+:\tscale configuration (larger)\n-:\tscale configuration (smaller)\n*:\trotate about x-axis\n/:\trotate about y-axis\na:\tauto-scale\nb:\tdraw bonds\nc:\tget configuration via socket\nd:\tget distribution via socket\ne:\ttoggle type/energy encoding\nh:\tdisplay this help message\nk:\ttoggle potential/kinetic energy\nl:\tload IMD-configuration from file\np\tmake picture\nq:\tquit program\nr:\ttoggle type encoding via radius\nt:\ttoggle text on/off\nu:\tread_unit_vectors\nw:\twrite to file\n");
  printf(sysstring);
  fflush(stdout);
	 
}


