/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

void write_pic_cell( cell *p, FILE *out ) 
{
    struct { 
        float    pos_x, pos_y, E_kin, E_pot;
        integer  type;
    } picbuf;

    real px, py;
    int i;

    for (i = 0;i < p->n; ++i) {
        picbuf.pos_x = (float) p->ort X(i);
        picbuf.pos_y = (float) p->ort Y(i);
        if ( pic_ur.x != (real)0 ) /*if pic_ur still 0, write everything */
        if ( (picbuf.pos_x < pic_ll.x) || (picbuf.pos_x > pic_ur.x) ||
             (picbuf.pos_y < pic_ll.y) || (picbuf.pos_y > pic_ur.y) ) continue;
        px = p->impuls X(i);
        py = p->impuls Y(i);
        picbuf.E_kin = (float) ( (px*px + py*py) / (2 * p->masse[i]) );
#ifdef DISLOC
        if (Epot_diff==1) {
          picbuf.E_pot = (float) p->pot_eng[i] - p->Epot_ref[i];
        } else
#else
        picbuf.E_pot = (float) p->pot_eng[i];
#endif
        picbuf.type  = (integer) p->sorte[i];
        fwrite( &picbuf, sizeof( picbuf ), 1, out ); 
    }
}


void write_pictures_raw(int steps)
{ 
  FILE *out;

  str255 fname;
  int fzhlr;
  cell *p,*q;
  int i,j,k,l,m,tag;

  /* Dateiname fuer Ausgabedatei erzeugen */
  fzhlr = steps / pic_interval;

#ifdef MPI  
  if (0<parallel_output)
    sprintf(fname,"%s.%u.%u.pic",outfilename,fzhlr,myid);
  else
#endif
    sprintf(fname,"%s.%u.pic",outfilename,fzhlr);


#ifdef MPI

  if (0<parallel_output) {

    /* Ausgabedatei oeffnen */
    out = fopen(fname,"w");
    if (NULL == out) error("Can't open output file for config.");


    for (j = 1; j < cell_dim.x-1; ++j )
      for (k = 1; k < cell_dim.y-1; ++k ) {
 	  p = PTR_2D_V(cell_array, j, k, cell_dim);
	  write_pic_cell( p, out );
	};
    
    fclose(out);

  } else { 

    if (0==myid) {

      /* Ausgabedatei oeffnen */
      out = fopen(fname,"w");
      if (NULL == out) error("Can't open output file for config.");

      /* Write data on CPU 0 */

      /* Write own data */
      for (j = 1; j < cell_dim.x-1; ++j )
	for (k = 1; k < cell_dim.y-1; ++k ) {
	  p = PTR_2D_V(cell_array, j, k, cell_dim);
	  write_pic_cell( p, out );
	};

      /* Receive data from other cpus and write that */
      p   = PTR_2D_V(cell_array, 0, 0, cell_dim);
      for ( m = 1; m < num_cpus; ++m)
	for (j = 1; j < cell_dim.x-1; ++j )
	  for (k = 1; k < cell_dim.y-1; ++k ) {
	    tag = PTR_2D_V(CELL_TAG, j, k, cell_dim);
	    recv_cell( p, m, tag );
	    write_pic_cell( p, out );
	  };

      fclose(out);      
    } else { 
      /* Send data to cpu 0 */
      for (j = 1; j < cell_dim.x-1; ++j )
	for (k = 1; k < cell_dim.y-1; ++k ) {
	  p   = PTR_2D_V(cell_array, j, k, cell_dim);
	  tag = PTR_2D_V(CELL_TAG, j, k, cell_dim);
	  send_cell( p, 0, tag );
	};
    };
  };
#else

  /* Ausgabedatei oeffnen */
  out = fopen(fname,"w");
  if (NULL == out) error("Can't open output file for config.");


  for (p = cell_array; 
       p <= PTR_2D_V(cell_array,
		     cell_dim.x-1,
		     cell_dim.y-1,
		     cell_dim);
       ++p ) {
    write_pic_cell( p , out );
  };

  fclose(out);  

#endif

}

/******************************************************************************
*
* write_pictures writes bitmap pictures of configuration
*
******************************************************************************/

void write_pictures_bins(int steps)

#define XRES 576
#define YRES 720
#define SFACTOR 1.0
#define NUMPIX  2

{

  vektor scale;
  real xshift, yshift;

  shortint redbit  [XRES][YRES];
  shortint greenbit[XRES][YRES];
  shortint bluebit [XRES][YRES];

#ifdef MPI
  shortint sum_red[XRES][YRES];
  shortint sum_green[XRES][YRES];
  shortint sum_blue[XRES][YRES];
  ivektor2d maxcoord,mincoord;
#endif

  char buf[3*YRES];

  str255 fname;
  int fzhlr;
  int i,j,k,r,s;
  real val;
  cell *p;
  real red,green,blue;
  FILE *out;
  ivektor2d coord,pixcoord;
  real tabred[5],tabgreen[5],tabblue[5];
  int ind;
  int pix;
  int np;

  /* the dist bins are orthogonal boxes in space */

  scale.x = XRES / box_x.x;
  scale.y = YRES / box_y.y;

  /* Make scaling same in both directions */
  if (scale.y<scale.x)  scale.x = scale.y;

  scale.x *= pic_scale.x;
  scale.y *= pic_scale.y;

  /* kinetic energy first */

  /* create filename */
  /* Dateiname fuer Ausgabedatei erzeugen */
  fzhlr = steps / pic_interval;

  sprintf(fname,"%s.%u.kin.ppm",outfilename,fzhlr);

  /* Center image */
  xshift = (XRES - box_x.x * scale.x) / 2.0;
  yshift = (YRES - box_y.y * scale.y) / 2.0; 


  /* Zero bitmap */
  for (i=0; i<XRES; i++ ) {
    for (j=0; j<YRES; j++ ) {
      bluebit[i][j]  = 0;
      greenbit[i][j] = 0;
      redbit[i][j]   = 0;
    };
  };

  /* Color lookup table */


/*  Original table from RUS - we use more saturation

    tabred[0] = 0.10; tabgreen[0] = 0.20; tabblue[0] = 0.50;
    tabred[1] = 0.05; tabgreen[1] = 0.75; tabblue[1] = 0.75;
    tabred[2] = 0.10; tabgreen[2] = 0.50; tabblue[2] = 0.25;
    tabred[3] = 0.75; tabgreen[3] = 0.75; tabblue[3] = 0.05;
    tabred[4] = 0.75; tabgreen[4] = 0.05; tabblue[4] = 0.05; */

    tabred[0] = 0.02; tabgreen[0] = 0.02; tabblue[0] = 0.45;
    tabred[1] = 0.03; tabgreen[1] = 0.23; tabblue[1] = 0.23;
    tabred[2] = 0.02; tabgreen[2] = 0.45; tabblue[2] = 0.02;
    tabred[3] = 0.23; tabgreen[3] = 0.23; tabblue[3] = 0.03;
    tabred[4] = 0.45; tabgreen[4] = 0.02; tabblue[4] = 0.02; 
  
  /* loop over all atoms */
    for ( r = cellmin.x; r < cellmax.x; ++r )
      for ( s = cellmin.y; s < cellmax.y; ++s ) {
      
	p = PTR_2D_V(cell_array, r, s, cell_dim);

	for (i = 0;i < p->n; ++i) {
	  coord.x = (int) (p->ort X(i) * scale.x) + xshift;
	  coord.y = (int) (p->ort Y(i) * scale.y) + yshift;
	  /* Check bounds */
	  if ((coord.x >= 0) && (coord.x < (XRES-NUMPIX)) &&
              (coord.y >= 0) && (coord.y < (YRES-NUMPIX))) { 
	  
             val = SPRODN(p->impuls,i,p->impuls,i) / (2*p->masse[i]);

             /* Scale Value to [0..1]   */
	     val = (val - ecut_kin.x) / (ecut_kin.y - ecut_kin.x);
             val = val > 1.0 ? 1.0 : val;
             val = val < 0.0 ? 0.0 : val;
             np  = NUMPIX;
             /* Get index into table */
	     ind = (int)(val * 3.9999);

             /* Get RBG values from linear interpolation */
	     red = -4.0 * (tabred[ind] - tabred[ind+1]) * val + 
	         (tabred[ind] * (ind+1) - ind * tabred[ind+1] );

	     green = -4.0 * (tabgreen[ind] - tabgreen[ind+1]) * val +
                 (tabgreen[ind] * (ind+1) - ind * tabgreen[ind+1] );

             blue = -4.0 * (tabblue[ind] - tabblue[ind+1]) * val +
                 (tabblue[ind] * (ind+1) - ind * tabblue[ind+1] );

             /* Set & Copy Pixel */
             for (j=0;j<np;++j)
               for (k=0;k<np;++k) {
                 pixcoord.x = coord.x + j;
                 pixcoord.y = coord.y + k;

                 pix =  redbit  [pixcoord.x][pixcoord.y] + (SFACTOR * 255 * red  );
                 redbit  [pixcoord.x][pixcoord.y] = (shortint) pix < 255 ? pix : 255;
	      
                 pix = bluebit [pixcoord.x][pixcoord.y] + (SFACTOR * 255 * blue );
                 bluebit [pixcoord.x][pixcoord.y] = (shortint) pix < 255 ? pix : 255;

                 pix = greenbit[pixcoord.x][pixcoord.y] + (SFACTOR * 255 * green);
                 greenbit[pixcoord.x][pixcoord.y] = (shortint) pix < 255 ? pix : 255;

             }; /* for k */
          }; /* if */
       }; /* for i */
    }; /* for r */

#ifdef MPI
/* Add the bitmaps */
   MPI_Reduce( redbit,   sum_red,   XRES * YRES, MPI_SHORT, MPI_SUM, 0, cpugrid);
   MPI_Reduce( greenbit, sum_green, XRES * YRES, MPI_SHORT, MPI_SUM, 0, cpugrid);
   MPI_Reduce( bluebit , sum_blue,  XRES * YRES, MPI_SHORT, MPI_SUM, 0, cpugrid);

if (0==myid) { 

/* Clip max value bitmap, create white background */
   for (i=0; i<XRES; i++ ) {
     for (j=0; j<YRES; j++ ) { 

      redbit[i][j]   = sum_red[i][j]   < 255 ? sum_red[i][j]   : 200;
      greenbit[i][j] = sum_green[i][j] < 255 ? sum_green[i][j] : 200;
      bluebit[i][j]  = sum_blue[i][j]  < 255 ? sum_blue[i][j]  : 200;
    };
  }; 

  
};
#endif

  for (i=0; i<XRES; i++ ) 
    for (j=0; j<YRES; j++ ) 
/* background white */
      if ((0==redbit[i][j]) && (0==greenbit[i][j]) && (0==bluebit[i][j])) {
         redbit[i][j]   = 250;
         greenbit[i][j] = 250;
         bluebit[i][j]  = 250;
      };


  /* write ppm file */

#ifdef MPI
  if (0==myid) {
#endif

  out = fopen(fname,"w");
  if (NULL == out) error("Can`t open bitmap file.");

  fprintf(out,"P6 %d %d 255\n", YRES, XRES);


  for (i=0; i<XRES; i++ ) {
    for (j=0; j<YRES; j++ ) {
      buf[3*j  ] = (char) redbit[i][j];
      buf[3*j+1] = (char) greenbit[i][j];
      buf[3*j+2] = (char) bluebit[i][j];
    };
    fwrite(buf, sizeof(char), 3*YRES, out);
  };

  fclose(out);

#ifdef MPI
 };
#endif

/* Potential energy second */


  /* create filename */

  sprintf(fname,"%s.%u.pot.ppm",outfilename,fzhlr);


  /* Zero bitmap */
  for (i=0; i<XRES; i++ ) {
    for (j=0; j<YRES; j++ ) {
      bluebit[i][j]  = 0;
      greenbit[i][j] = 0;
      redbit[i][j]   = 0;
    };
  };

  /* loop over all atoms */
  for ( r = cellmin.x; r < cellmax.x; ++r )
    for ( s = cellmin.y; s < cellmax.y; ++s )
      {

	p = PTR_2D_V(cell_array, r, s, cell_dim);

	for (i = 0;i < p->n; ++i) {
	  coord.x = (int) (p->ort X(i) * scale.x) + xshift;
	  coord.y = (int) (p->ort Y(i) * scale.y) + yshift;
	  /* Check bounds */
	  if ((coord.x>=0) && (coord.x<(XRES-NUMPIX)) &&
	      (coord.y>=0) && (coord.y<(YRES-NUMPIX))) {

	  val = p->pot_eng[i];

          /* Scale Value to [0..1]   */
	  val = (val - ecut_pot.x) / (ecut_pot.y - ecut_pot.x);
/* Values that are not in the interval are set to MINIMUM */
          val = val > 1.0 ? 0.0 : val;
          val = val < 0.0 ? 0.0 : val;
/* Defects in potential energy are rather point-like, so we enlarge all
Pixels not in the default interval */
          if ((val<1.0) && (val>0.0)) np = 3 * NUMPIX; else np = NUMPIX;
          /* Get index into table */
	  ind = (int)(val * 3.9999);

	  red = -4.0 * (tabred[ind] - tabred[ind+1]) * val + 
	    (tabred[ind] * (ind+1) - ind * tabred[ind+1] );

	  green = -4.0 * (tabgreen[ind] - tabgreen[ind+1]) * val + 
	    (tabgreen[ind] * (ind+1) - ind * tabgreen[ind+1] );

          blue = -4.0 * (tabblue[ind] - tabblue[ind+1]) * val + 
	    (tabblue[ind] * (ind+1) - ind * tabblue[ind+1] );


          /* Set & Copy Pixel */
          for (j=0;j<np;++j)
             for (k=0;k<np;++k) {
              pixcoord.x = coord.x + j;
              pixcoord.y = coord.y + k;

              if ((pixcoord.x<XRES) && (pixcoord.y<YRES)) {

                 pix =  redbit  [pixcoord.x][pixcoord.y] + (SFACTOR * 255 * red  );
                 redbit  [pixcoord.x][pixcoord.y] = (shortint) pix < 255 ? pix : 255;
	      
                 pix = bluebit [pixcoord.x][pixcoord.y] + (SFACTOR * 255 * blue );
                 bluebit [pixcoord.x][pixcoord.y] = (shortint) pix < 255 ? pix : 255;

                 pix = greenbit[pixcoord.x][pixcoord.y] + (SFACTOR * 255 * green);
                 greenbit[pixcoord.x][pixcoord.y] = (shortint) pix < 255 ? pix : 255;
              };
           }; 
	};
       };
    };

#ifdef MPI
/* Add the bitmaps */
   MPI_Reduce( redbit,   sum_red,   XRES * YRES, MPI_SHORT, MPI_SUM, 0, cpugrid);
   MPI_Reduce( greenbit, sum_green, XRES * YRES, MPI_SHORT, MPI_SUM, 0, cpugrid);
   MPI_Reduce( bluebit , sum_blue,  XRES * YRES, MPI_SHORT, MPI_SUM, 0, cpugrid);

if (0==myid) { 

/* Clip max value bitmap, create white background */
   for (i=0; i<XRES; i++ ) {
     for (j=0; j<YRES; j++ ) { 

      redbit[i][j]   = sum_red[i][j]   < 255 ? sum_red[i][j]   : 200;
      greenbit[i][j] = sum_green[i][j] < 255 ? sum_green[i][j] : 200;
      bluebit[i][j]  = sum_blue[i][j]  < 255 ? sum_blue[i][j]  : 200;
    };
  }; 

  for (i=0; i<XRES; i++ ) 
    for (j=0; j<YRES; j++ ) 
/* background white */
      if ((0==redbit[i][j]) && (0==greenbit[i][j]) && (0==bluebit[i][j])) {
         redbit[i][j]   = 250;
         greenbit[i][j] = 250;
         bluebit[i][j]  = 250;
      };
  
};
#endif

  /* write ppm file */

#ifdef MPI
  if (0==myid) {
#endif

  out = fopen(fname,"w");
  if (NULL == out) error("Can`t open bitmap file.");

  fprintf(out,"P6 %d %d 255\n", YRES, XRES);


  for (i=0; i<XRES; i++ ) {
    for (j=0; j<YRES; j++ ) {
      buf[3*j  ] = (char) redbit[i][j];
      buf[3*j+1] = (char) greenbit[i][j];
      buf[3*j+2] = (char) bluebit[i][j];
    };
    fwrite(buf, sizeof(char), 3*YRES, out);
  };

  fclose(out);

#ifdef MPI
 };
#endif
}


void write_pictures( steps )
{
  switch (pic_type) {
    case 0: write_pictures_raw ( steps ); break;
    case 1: break;
    case 2: write_pictures_bins( steps ); break;
  };
}

