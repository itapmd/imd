
/* ---------------------------------------------------------

   COMPILATION:  

       gcc -O -o imd2ppm imd2ppm.c -lm

   USAGE:  

       imd2ppm [-p] [-a] input Emin Emax x0 y0 b0 h0  b1 h1

       where:

           -p           use potential energy for the coloring,
                        rather than kinetic energy
           -a           expect ASCII-formatted IMD configuration files 
                        as input, rather than the binary .pic files
           input        name of input file (output = input.ppm)
           Emin, Emax   bounderies of energy interval for coloring
           x0 y0        lower left corner of selected range in user coords
           b0 h0        width and height in original coordinates
           b1 h1        width and height of screen output (in pixels)


   INPUT FORMAT:

       Input files can be ASCII-formatted IMD configuration files
       (flag -a), or binary .pic files. The latter consist of 
       structs of the following form:

       struct { 
           float    pos_x, pos_y, E_kin, E_pot;
           integer  type;
       } picbuf;

   -----------------------------------------------------------
*/

#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef struct { char r,g,b; } Triple;

/* 
  procedures to draw a circle
  taken from: Grundlagen der Computergrafik,
              Foley, van Dam and others
              page 97
*/

int ok(int x, int y, int col, int row)
    {
      return ((x>=0) && (y>=0) && (x<col) && (y<row));
    }

void CirclePoints(int cx,int cy, int mx, int my, int col, int row, 
                  Triple value, Triple *matrix )
    { 
      int v;
      /* outer borders... */
      for (v=mx-cx;v<=mx+cx;v++){
	if (ok(v,my+cy,col,row))  matrix[(my+cy)*col + v]=value;}
      for (v=mx-cx;v<=mx+cx;v++){
	if (ok(v,my-cy,col,row)) matrix[(my-cy)*col + v]=value;}     
      for (v=mx-cy;v<=mx+cy;v++){
        if (ok(v,my+cx,col,row)) matrix[(my+cx)*col + v]=value;}
      for (v=mx-cy;v<=mx+cy;v++){
        if (ok(v,my-cx,col,row)) matrix[(my-cx)*col + v]=value;}

    }

void MidpointCircle(int radius, int mx, int my, int col, int row,
		    Triple value, Triple *matrix)
    {
      int xx,yy,dd;

      /* initialize */

      xx=0;
      yy=radius;
      dd=1-radius;

      CirclePoints(xx,yy,mx,my,col,row,value,matrix);
      while (yy>xx) {
	if (dd<0) {
	  dd+=xx*2 + 3;
	  xx++;
	}else{
	  dd+=(xx-yy)*2 + 5;
	  xx++;
	  yy--;
	}
	CirclePoints(xx,yy,mx,my,col,row,value,matrix);
      }
    }

void MidpointLine(int x0, int y0, int dx, int dy, int col, int row, 
                  Triple value, Triple *matrix)
{
  int x1,y1,incrE,incrNE,d,x,y;
  int StepX, StepY;

  x1=x0+dx;
  y1=y0+dy;

  if ((dx!=0) || (dy!=0)) { 

  if (dx>0) StepX=1;
  if (dx==0) StepX=0;
  if (dx<0) {StepX=-1; dx=-dx;}
	     
  if (dy>0) StepY=1;
  if (dy==0) StepY=0;
  if (dy<0) {StepY=-1; dy=-dy;}


  d=dy*2-dx;
  incrE=dy*2;
  incrNE=(dy-dx)*2;
  x=x0;
  y=y0;
  if (ok(x,y,col,row)) matrix[y*col + x] = value;
  while ((x!=x0+dx*StepX)&&(y!=y0+dy*StepY)) {
    if (d<=0){
      d+=incrE;
      x=x+StepX;
    }else{
      d+=incrNE;
      x=x+StepX;
      y=y+StepY;
    }
    if (ok(x,y,col,row)) matrix[y*col + x] = value;
  }
  }
}

/******************* main *******************/

int main (int argc, char **argv) 

{
  /* (n lines and m columns) */
  int      i,j,col,row,p,q;
  double   x0,y0,b0,h0;
  double   Emin, Emax, E, H;
  double   fl_t;

  int      ascii=0, epot=0;  /* default is .pic files, Ekin for coloring */

  int      dx,dy,pointok,vectorok,num;

  int      maxcolor = 255;
  int      maxsize  = 15;  /* not used yet */

  double   mass,x,y,vx,vy,color_r,color_g,color_b,size;
  double   bratio, hratio;

  double   tabred[2][5],tabgreen[2][5],tabblue[2][5];

  int      ix,iy,isize;
  int      ivx,ivy;

  char     rbuf[512];
  char     *sourcename, *targetname;
  Triple   *matrix, *TopMatrix;
  Triple   point, BackColor, LineColor;

#ifdef CRAY
typedef short int integer;
#else
typedef       int integer;
#endif

  struct { 
    float    pos_x, pos_y, E_kin, E_pot;
    integer  type;
  } picbuf;

  FILE    *source, *target;
  char    *buf; 
 
  double  dE,val,red,green,blue,dummy,sat_r,sat_g,sat_b;
  int     t,maxc,ind;

  while ((argc > 1) && (argv[1][0] =='-')) {
    if (argv[1][1]=='a') {
      ascii=1;
      argc--;
      argv++;
    }
    else if (argv[1][1]=='p') {
      epot=1;
      argc--;
      argv++;
    }
  }

  if (argc < 10) {
    printf("Usage: imd2ppm [-p] [-a] Infile Emin Emax x0 y0 b0 h0 b1 h1\n");
    exit(1);
  }

  sourcename = argv[1];
  targetname = malloc( strlen(sourcename)+9 );
  strcpy(targetname,sourcename);
  if (epot) {
    strcat(targetname,".pot.ppm");
  } else {
    strcat(targetname,".kin.ppm");
  }
  sscanf(argv[2], "%lf", &Emin);
  sscanf(argv[3], "%lf", &Emax); 
  sscanf(argv[4], "%lf", &x0);
  sscanf(argv[5], "%lf", &y0);
  sscanf(argv[6], "%lf", &b0);
  sscanf(argv[7], "%lf", &h0);
  sscanf(argv[8], "%d",  &col);
  sscanf(argv[9], "%d",  &row);

  bratio=(double)col/b0;
  hratio=(double)row/h0;

  matrix    = calloc( col*row, sizeof( Triple ) );

  /* color table token from IMD sourcecode, imd_io_2d.c */
 
  /* color typ 0 */
  tabred[0][0] = 1.00; tabgreen[0][0] = 1.00; tabblue[0][0] = 1.00;
  tabred[0][1] = 1.00; tabgreen[0][1] = 0.75; tabblue[0][1] = 0.75;
  tabred[0][2] = 1.00; tabgreen[0][2] = 0.50; tabblue[0][2] = 0.50;
  tabred[0][3] = 1.00; tabgreen[0][3] = 0.25; tabblue[0][3] = 0.25;
  tabred[0][4] = 1.00; tabgreen[0][4] = 0.02; tabblue[0][4] = 0.02; 

  /* color typ 1 */
  tabred[1][0] = 0.02; tabgreen[1][0] = 0.02; tabblue[1][0] = 1.00;
  tabred[1][1] = 0.25; tabgreen[1][1] = 0.02; tabblue[1][1] = 0.75;
  tabred[1][2] = 0.50; tabgreen[1][2] = 0.02; tabblue[1][2] = 0.50;
  tabred[1][3] = 0.75; tabgreen[1][3] = 0.02; tabblue[1][3] = 0.25;
  tabred[1][4] = 1.00; tabgreen[1][4] = 0.02; tabblue[1][4] = 0.02; 

  BackColor.r=0;
  BackColor.g=0;
  BackColor.b=0;

  LineColor.r=255;  
  LineColor.g=255;  
  LineColor.b=255;

  for (i=0;i<(col*row);i++) {
      matrix[i]    = BackColor;
  }

  /*****************************************/

  /* open source file */
  source = fopen( sourcename, "r" );
  if ( NULL == source ) {
     printf("input file not found\n");
     return 1;
  }

  /* eat header in .pic files */
  if (0==ascii) {
    do {
      fgets(rbuf,sizeof(rbuf),source);
    } while (rbuf[2]!='E');
  }

  do   /* read one atom at a time */
  {
    /* input in ASCII format */
    if (ascii) {
      rbuf[0] = (char) NULL;
      fgets(rbuf,sizeof(rbuf),source);
      while ('#'==rbuf[0]) fgets(rbuf,sizeof(rbuf),source); /* eat comments */
      p = sscanf(rbuf,"%d %d %lf %lf %lf %lf %lf %lf",
	              &num,&t,&mass,&x,&y,&vx,&vy,&E);
      if (epot) {
        /* keep the potential energy in variable E to compute the color */
        if (p<8) break; 
      }
      else {
        if (p<7) break;
        /* compute the color using the kinetic energy */
        E = 0.5 * mass * (vx*vx + vy*vy); 
      }
    } 
    /* not ASCII - binary .pic input format */
    else { 
      if (1 != fread( &picbuf, sizeof( picbuf ), 1, source ) ) break;
      t  = picbuf.type;
      x  = picbuf.pos_x;
      y  = picbuf.pos_y;
      if (epot) {
        E  = picbuf.E_pot;   /* compute the color using the kinetic energy */
      } else {
        E  = picbuf.E_kin;   /* compute the color using the potential energy */
      }
    }

    dE = Emax-Emin;

    if (t==1) {
      size=0.2; 
    } else {
      t=0; 
      size=0.1;
    }

    val=(E-Emin)/dE;
    if (val<0) val=0;
    if (val>1) val=1;    
    ind = (int) 3.999*val;

    /* interpolate linearly between tabulated values */
    red   = tabred  [t][ind]+(4*val-ind)*(tabred  [t][ind+1]-tabred  [t][ind]);
    green = tabgreen[t][ind]+(4*val-ind)*(tabgreen[t][ind+1]-tabgreen[t][ind]);
    blue  = tabblue [t][ind]+(4*val-ind)*(tabblue [t][ind+1]-tabblue [t][ind]);

    maxc  = 215;
    sat_r = 1;
    sat_g = 0.85;
    sat_b = 1; 

    color_r = (int)(red*sat_r*maxc);
    color_g = (int)(green*sat_g*maxc);
    color_b = (int)(blue*sat_b*maxc);

    if (color_r > maxcolor) color_r=maxcolor;
    if (color_g > maxcolor) color_g=maxcolor;
    if (color_b > maxcolor) color_b=maxcolor;    
    if (color_r < 0) color_r=0;    
    if (color_g < 0) color_g=0;
    if (color_b < 0) color_b=0;

    point.r = (int)floor(color_r);
    point.g = (int)floor(color_g);
    point.b = (int)floor(color_b);

    ix = (int)floor((x-x0)*bratio);  
    iy = (int)floor((y-y0)*hratio);  

    /* don't draw velocities
    ivx = (int)floor(vx*bratio);  
    ivy = (int)floor(vy*hratio);
    */

    isize = (int)ceil(size*(hratio + bratio)/2);
    if (isize <1) isize = 1; 

    /* draw if point can be seen on screen */
    if ((ix>-isize) && (iy>-isize) && (ix<col+isize) && (iy<row+isize)) {
        MidpointCircle(isize, ix, iy, col,row, point, matrix);
        /* don't draw velocities  
        MidpointLine(ix,iy,ivx,ivy,col,row, LineColor, TopMatrix); 
        */
    }
  }
  while (1);
  fclose(source);

  /* open output file */
  target = fopen( targetname, "w" );
  if ( NULL == target ) {
     printf("can not open output file\n");
     return;
  };

  /* write header of ppm file */
  fprintf( target, "P6\n%d %d\n255\n", col, row );

  /* write data */
  buf = malloc(3*col);
  for (i=0; i<row; i++) {
     for (j=0; j<col; j++) {
        point      = matrix[i*col+j];
        buf[3*j]   = point.r;
        buf[3*j+1] = point.g;
        buf[3*j+2] = point.b;
     };
     fwrite( buf, 1, 3*col, target ); 
  };

  /* write data with velocity lines */
  /*
  for (i=0; i<col*row; i++) {
    if ((TopMatrix[i].r == BackColor.r) &&
        (TopMatrix[i].g == BackColor.g) &&
	(TopMatrix[i].b == BackColor.b)) {
        fwrite( matrix[i], sizeof( Triple ), 1, target ); 
    }
    else {
        fwrite( TopMatrix[i], sizeof( Triple ), 1, target ); 
    }
  }
  */

  fclose( target );
  exit(0);

}


 













