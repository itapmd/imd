
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/* ---------------------------------------------------------

   COMPILATION:  

       gcc -O -o imddist2ppm [-DASCII] imddist2ppm.c -lm

       OPTIONS:

           -DASCII:  expect ASCII-formatted IMD configuration files 
                     as input, rather than the binary .pic files

   USAGE:  

       imd2ppm input Emin Emax xres yres [color_table]

       where:

           input        name of input file (output = input.ppm)
	   Emin, Emax   bounderies of energy interval for coloring 
	   xres, yres   resolution of the distributen in x-, y-direction
	   color_table  name of color table file (optional)
	   
   COLOR TABLE FORMAT
  
       the color table consist of 10  lines 
       red_1  green_2   blue_3 
        .       .         .
	.       .         .
	
       red_10 green _10 blue_10 

       where red_i, green_i, blue_i are values between 0.0 an 1.0 



   INPUT FORMAT:

       Input files can be ASCII-formatted IMD configuration files
       (compilation switch -DASCII), or binary files. 

   -----------------------------------------------------------
*/

#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


typedef short int shortint;

/******************* main *******************/

int main (int argc, char **argv) 

{

  int      i, j;
  int      xres, yres;
  double   Emin, Emax;
  float    E, dE, val;
  int      p;

  int      maxcolor = 255;
  int      ind;

  float  tabred [10],tabgreen [10],tabblue [10];

  shortint *redbit;
  shortint *greenbit;
  shortint *bluebit;
  
  float red, green, blue;


  char     rbuf[512];
  char     tabbuf[512];
  char     *sourcename, *targetname, *tabfilename;

#ifdef CRAY
typedef short int integer;
#else
typedef       int integer;
#endif

  FILE    *source, *target, *tabfile;
  char    *buf; 
 

  if (argc < 5) {
     printf("\n------> not enough parameters <----- \n");
     printf("Infile Emin Emax xres yres\n\n");
     exit(1);
  }

  sourcename = argv[1];
  targetname = malloc( strlen(sourcename)+5 );
  strcpy(targetname,sourcename);
  strcat(targetname,".ppm");
  
  sscanf(argv[2], "%lf", &Emin);
  sscanf(argv[3], "%lf", &Emax);
  sscanf(argv[4], "%d", &xres);
  sscanf(argv[5], "%d", &yres); 
   
  /*  printf("%f\n", Emin);
  printf("%f\n", Emax);
  printf("%d\n", xres);
  printf("%d\n", yres);*/



  dE = Emax-Emin;

  redbit = (shortint*)calloc(xres*yres,sizeof(shortint));
  greenbit = (shortint*)calloc(xres*yres,sizeof(shortint));
  bluebit = (shortint*)calloc(xres*yres,sizeof(shortint));

  if (argc > 6){
    tabfilename = argv[6];
    
    /* open  tabfile */
    tabfile = fopen( tabfilename, "r" );
    if ( NULL == tabfile ) {
      printf("tabfile file not found\n");
      return 1;
    };

    /* read table */

    i=0;
    do 
      {
	tabbuf[0] = (char) NULL;
	fgets(tabbuf,sizeof(tabbuf),tabfile);
	while ('#'==tabbuf[0]) fgets(tabbuf,sizeof(tabbuf),tabfile); /* eat comments */
	p = sscanf(tabbuf,"%f %f %f",\
		   &red,&green,&blue);
       	

	if(p< 3) break;
	tabred  [i] = red;
	tabgreen[i] = green;
	tabblue [i] = blue;
	/*printf("%f %f %f %d\n",\
	       tabred[i],tabgreen[i],tabblue[i], i); */
	i++;
	
      }
    while (1);
    if (i<10){
      printf("color table too short\n\n"); exit(1);
    }



    fclose(tabfile);
    
  }else { 
 
  
  tabred [0] = 0.00; tabgreen [0] = 0.00; tabblue [0] = 1.00;
  tabred [1] = 0.00; tabgreen [1] = 0.50; tabblue [1] = 1.00;
  tabred [2] = 0.00; tabgreen [2] = 1.00; tabblue [2] = 1.00;
  tabred [3] = 0.00; tabgreen [3] = 1.00; tabblue [3] = 0.50;
  tabred [4] = 0.00; tabgreen [4] = 1.00; tabblue [4] = 0.00; 
  tabred [5] = 0.50; tabgreen [5] = 1.00; tabblue [5] = 0.00;
  tabred [6] = 1.00; tabgreen [6] = 1.00; tabblue [6] = 0.00;
  tabred [7] = 1.00; tabgreen [7] = 0.50; tabblue [7] = 0.00;
  tabred [8] = 1.00; tabgreen [8] = 0.00; tabblue [8] = 0.00;
  tabred [9] = 1.00; tabgreen [9] = 0.20; tabblue [9] = 0.20;  

  }
  
  
  /*****************************************/

  /* open source file */
  source = fopen( sourcename, "r" );
  if ( NULL == source ) {
     printf("input file not found\n");
     return 1;
  };


  i = 0;

  do   /* read one energie at a time */
  {

#ifdef ASCII /* input in ASCII format */

    rbuf[0] = (char) NULL;
    fgets(rbuf,sizeof(rbuf),source);
    while ('#'==rbuf[0]) fgets(rbuf,sizeof(rbuf),source); /* eat comments */
    p = sscanf(rbuf,"%f",&E);
    if (p<1) break;

#else /* not ASCII - binary  input format */

    if (1 != fread( &E, sizeof( float ), 1, source ) ) break;

#endif /* ASCII */




    val=(E-Emin)/dE;
    if (val<0) val=0;
    if (val>1) val=1;    
    ind = (int) 8.999*val;

    red   = tabred   [ind] + (8*val-ind)*(tabred   [ind+1]-tabred   [ind]);
    green = tabgreen [ind] + (8*val-ind)*(tabgreen [ind+1]-tabgreen [ind]);
    blue  = tabblue  [ind] + (8*val-ind)*(tabblue  [ind+1]-tabblue  [ind]);

    redbit  [i] = (shortint) (red   * maxcolor);
    greenbit[i] = (shortint) (green * maxcolor);
    bluebit [i] = (shortint) (blue  * maxcolor);

    /* ueberfluessige vorsichtsmassnahme */
    if (redbit  [i] > maxcolor) redbit  [i] = maxcolor;
    if (greenbit[i] > maxcolor) greenbit[i] = maxcolor;
    if (bluebit [i] > maxcolor) bluebit [i] = maxcolor;    
    if (redbit  [i]  < 0)  redbit  [i] = 0;    
    if (greenbit[i]  < 0)  greenbit[i] = 0;
    if (bluebit [i]  < 0)  bluebit [i] = 0;


    /* background */
    if (E == 0.0){
      redbit  [i] = 0;
      greenbit[i] = 0;
      bluebit [i] = 0;
    }
    
    i++;
  }
  
  while (1);
  fclose(source);

  /*  for(i=0;i<5000;i++)
    printf("%d\t%d\t%d\t%d\n", i, redbit[i], greenbit[i], bluebit[i]);
    */



  /* open output file */
  target = fopen( targetname, "w" );
  if ( NULL == target ) {
     printf("can not open output file\n");
     return -1;
  };

  /* write header of ppm file */
  fprintf( target, "P6\n%d %d\n255\n", xres, yres );

  /* write data */
  buf = malloc(3*xres);
  for (i=0; i<yres; i++) {
     for (j=0; j<xres; j++) {
        buf[3*j]   = redbit  [j*yres+i];
        buf[3*j+1] = greenbit[j*yres+i];
        buf[3*j+2] = bluebit [j*yres+i];
     };
     fwrite( buf, 1, 3*xres, target ); 
  };

  fclose( target );
  exit(0);

}


 













