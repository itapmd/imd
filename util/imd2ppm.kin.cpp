
/* the most necessary routines to read and write a matrix of triples
   the main purpose is only to show how this feature can be done in C++
   basics by Andreas Menzel on 09/04/1997
   written by Christof Horn
*/


/* ---------------------------------------------------------

   DATA FORMAT:

   <command line>  input output x0 y0 b0 h0  b1 h1

   with input : name of input file (see below)
                (output = input.ppm)
        Emin Emax : bounderies of kinetic energie
        x0 y0 : left upper edge of selected range in original coord.
        b0 h0 : width and height in original coordinates
        b1 h1 : width and height of screen output

   <input file> 

   every row describes one atom with 

   x  y vx vy  r g b  size

   with  x y   : coordinates
         r g b : color of atom
         size  : size of atom in original coordinates  

   


   -----------------------------------------------------------
*/


#include <math.h>
#include "strstream.h"
#include "fstream.h"
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>




class Triple 
 {
  public: 
  char r,g,b;
 };

 /* 
     procedures to draw a circle
     taken out of: Grundlagen der Computergrafik,
                   Foley, van Dam and others
                   page 97
  */

int ok(int x, int y, int col, int row)
    {
      return ((x>=0) && (y>=0) && (x<col) && (y<row));
    }

void CirclePoints(int cx,int cy, int mx, int my, int col, int row,  Triple value, Triple *matrix )
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



void MidpointCircle(int radius, int mx, int my, int col, int row, \
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


void MidpointLine(int x0, int y0, int dx, int dy, int col, int row,  Triple value, Triple *matrix)
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



void main (int argc, char **argv) 

{
		/* (n lines and m columns) */
  int      i,j,col,row,p,q;
  double   x0,y0,b0,h0;
  double   Emin, Emax, Ekin;
  double   fl_t;

  int      dx,dy,pointok,vectorok;

  int      maxcolor = 255;
  int      maxsize  = 15;  /* not used yet */

  double   x,y,vx,vy,color_r,color_g,color_b,size;
  double   bratio, hratio;

  double   tabred[2][5],tabgreen[2][5],tabblue[2][5];

  int      ix,iy,isize;
  int      ivx,ivy;

  char    *sourcename,
          *targetname,
           pixel[1];
  Triple  *matrix, *point;
  Triple  *TopMatrix;
  Triple  *BackColor;
  Triple  *LineColor;
 
  double  dE,val,red,green,blue,dummy,sat_r,sat_g,sat_b;
  int     t,maxc,ind;




 
  if (argc < 10) 
    {
    cout << "\n------> Zuwenig Parameter <----- \n"
         << "Infile Emin Emax x0 y0 b0 h0 b1 h1\n\n\n";
    exit(1);
    }

  sourcename = argv[1];
  targetname = new char[strlen(sourcename)+5];
  strcpy(targetname,sourcename);
  strcat(targetname,".ppm");

  FILE *source;
  ofstream target (targetname, ios::out);

  sscanf(argv[2], "%lf", &Emin);
  sscanf(argv[3], "%lf", &Emax); 
  sscanf(argv[4], "%lf", &x0);
  sscanf(argv[5], "%lf", &y0);
  sscanf(argv[6], "%lf", &b0);
  sscanf(argv[7], "%lf", &h0);
  sscanf(argv[8], "%d", &col);
  sscanf(argv[9], "%d", &row);

  bratio=(double)col/b0;
  hratio=(double)row/h0;

  matrix = new Triple[col*row];
  TopMatrix = new Triple[col*row];
  BackColor = new Triple[1];
  LineColor = new Triple[1];
  point = new Triple[1];

    // color table token from IMD sourcecode, imd_io_2d.c
 
    // color typ 0
    tabred[0][0] = 1.00; tabgreen[0][0] = 1.00; tabblue[0][0] = 1.00;
    tabred[0][1] = 1.00; tabgreen[0][1] = 0.75; tabblue[0][1] = 0.75;
    tabred[0][2] = 1.00; tabgreen[0][2] = 0.50; tabblue[0][2] = 0.50;
    tabred[0][3] = 1.00; tabgreen[0][3] = 0.25; tabblue[0][3] = 0.25;
    tabred[0][4] = 1.00; tabgreen[0][4] = 0.02; tabblue[0][4] = 0.02; 

    // color typ 1
    tabred[1][0] = 0.02; tabgreen[1][0] = 0.02; tabblue[1][0] = 1.00;
    tabred[1][1] = 0.25; tabgreen[1][1] = 0.02; tabblue[1][1] = 0.75;
    tabred[1][2] = 0.50; tabgreen[1][2] = 0.02; tabblue[1][2] = 0.50;
    tabred[1][3] = 0.75; tabgreen[1][3] = 0.02; tabblue[1][3] = 0.25;
    tabred[1][4] = 1.00; tabgreen[1][4] = 0.02; tabblue[1][4] = 0.02; 


  BackColor[1].r=0;
  BackColor[1].g=0;
  BackColor[1].b=0;

  LineColor[1].r=255;  
  LineColor[1].g=255;  
  LineColor[1].b=255;

  for (i=0;i<(col*row);i++) {
      matrix[i]=BackColor[1];
      TopMatrix[i]=BackColor[1];
  }

  /*****************************************/

    /* reading of sourcefile */
    source = fopen( sourcename, "r" );
    if ( NULL == source ) {
        printf("input file not found\n");
        return;
    };

#ifdef CRAY
typedef short int integer;
#else
typedef       int integer;
#endif

    struct { 
        float    pos_x, pos_y, E_kin, E_pot;
        integer  type ;
    } picbuf;

  do 
  {
    if (1 != fread( &picbuf, sizeof( picbuf ), 1, source ) ) break;

    t    = picbuf.type;
    x    = picbuf.pos_x;
    y    = picbuf.pos_y;
    Ekin = picbuf.E_kin;

    /* Compute color using kinetic energy */

    dE=Emax-Emin;
   

    if (t==1) 
      size=0.2; 
    else {
      t=0; 
      size=0.1;
	}

    val=(Ekin-Emin)/dE;
    if (val<0) val=0;
    if (val>1) val=1;    
    // Umrechnung auf [0..1] 
    ind = int(3.999*val);
    // Bestimmung des Abschnittes [0 .. 4]
    if (ind>4)   ind=4;
    if (ind<0)   ind=0;

  red   = -4.0*(tabred[t][ind] - tabred[t][ind+1])*val + (tabred[t][ind] * (ind+1) - ind*tabred[t][ind+1]);

  green = -4.0*(tabgreen[t][ind] - tabgreen[t][ind+1])*val +(tabgreen[t][ind]*(ind+1) - ind*tabgreen[t][ind+1]);

  blue  = -4.0*(tabblue[t][ind] - tabblue[t][ind+1])*val +(tabblue[t][ind]*(ind+1) - ind*tabblue[t][ind+1]);  

 maxc  = 215;
 sat_r = 1;
 sat_g = 0.85;
 sat_b = 1; 

 color_r = (int)(red*sat_r*maxc);
 color_g = (int)(green*sat_g*maxc);
 color_b = (int)(blue*sat_b*maxc);


    /* without velocity vectors  */

    vx=0;
    vy=0;


    ix = (int)floor((x-x0)*bratio);  
    iy = (int)floor((y-y0)*hratio);  

    ivx = (int)floor(vx*bratio);  
    ivy = (int)floor(vy*hratio);
    
    if (color_r > maxcolor) color_r=maxcolor;
    if (color_g > maxcolor) color_g=maxcolor;
    if (color_b > maxcolor) color_b=maxcolor;    
    if (color_r < 0) color_r=0;    
    if (color_g < 0) color_g=0;
    if (color_b < 0) color_b=0;

    point[1].r = (int)floor(color_r);
    point[1].g = (int)floor(color_g);
    point[1].b = (int)floor(color_b);

    isize = (int)ceil(size*(hratio + bratio)/2);

    if (isize <1) isize = 1; 
//  if (isize > maxsize)  isize=maxsize;

   /* check if point can be seen on screen
      x = 0..col-1 
      y = 0..row-1  */

    pointok=((ix>-isize) && (iy>-isize) && (ix<col+isize) && (iy<row+isize));

   /* write lines into TopMatrix */ 

    if (pointok) {
      MidpointLine(ix,iy,ivx,ivy,col,row, LineColor[1], TopMatrix);
    }
   /* write circle */

    if (pointok) {
      MidpointCircle(isize, ix, iy, col,row, point[1], matrix);
    }
  }
  while (source);


   /*****************************************/

  fclose(source);
  

   /*****************************************/



  /* write matrix to target */

  
  target << "P6 \n"                  /* header for PPM P6 format */
         << col << " " << row << "\n"  /* width height */
         << "255\n";                /* max. number of colors */

  for (i=0; i<col*row; i++) {
    if ((TopMatrix[i].r== BackColor[1].r)&&\
	(TopMatrix[i].r== BackColor[1].r)&&\
	(TopMatrix[i].r== BackColor[1].r)) {
      target << matrix[i].r << matrix[i].g << matrix[i].b;
    }else{
      target << TopMatrix[i].r << TopMatrix[i].g << TopMatrix[i].b;
    }
    }

  target.close();

}


 













