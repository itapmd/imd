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

#define XRES 720
#define YRES 576
#define SFACTOR 1.0
#define NUMPIX  2

{

  vektor scale;
  real xshift, yshift;

  shortint redbit  [YRES][XRES];
  shortint greenbit[YRES][XRES];
  shortint bluebit [YRES][XRES];

#ifdef MPI
  shortint sum_red[YRES][XRES];
  shortint sum_green[YRES][XRES];
  shortint sum_blue[YRES][XRES];
  ivektor2d maxcoord,mincoord;
#endif

  char buf[3*XRES];

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

  /* compute pic_scale */
  pic_scale.x = box_x.x / (conf_ur.x - conf_ll.x);
  pic_scale.y = box_y.y / (conf_ur.y - conf_ll.y);

  /* Make scaling same in both directions */
  if (scale.y<scale.x)  scale.x = scale.y;

  scale.x *= pic_scale.x;
  scale.y *= pic_scale.y;

  /* kinetic energy first */

  /* create filename */
  /* Dateiname fuer Ausgabedatei erzeugen */
  fzhlr = steps / pic_interval;

  sprintf(fname,"%s.%u.kin.ppm",outfilename,fzhlr);

  /* Zero bitmap */
  for (j=0; j<YRES; j++ ) {
    for (i=0; i<XRES; i++ ) {
      bluebit[j][i]  = 0;
      greenbit[j][i] = 0;
      redbit[j][i]   = 0;
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

                 pix =  redbit  [pixcoord.y][pixcoord.x] + (SFACTOR * 255 * red  );
                 redbit  [pixcoord.y][pixcoord.x] = (shortint) pix < 255 ? pix : 255;
	      
                 pix = bluebit [pixcoord.y][pixcoord.x] + (SFACTOR * 255 * blue );
                 bluebit [pixcoord.y][pixcoord.x] = (shortint) pix < 255 ? pix : 255;

                 pix = greenbit[pixcoord.y][pixcoord.x] + (SFACTOR * 255 * green);
                 greenbit[pixcoord.y][pixcoord.x] = (shortint) pix < 255 ? pix : 255;

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
   for (j=0; j<YRES; j++ ) {
     for (i=0; i<XRES; i++ ) { 

      redbit[j][i]   = sum_red[j][i]   < 255 ? sum_red[j][i]   : 200;
      greenbit[j][i] = sum_green[j][i] < 255 ? sum_green[j][i] : 200;
      bluebit[j][i]  = sum_blue[j][i]  < 255 ? sum_blue[j][i]  : 200;
    };
  }; 

  
};
#endif

  for (j=0; j<YRES; j++ ) 
    for (i=0; i<XRES; i++ ) 
/* background white */
      if ((0==redbit[j][i]) && (0==greenbit[j][i]) && (0==bluebit[j][i])) {
         redbit[j][i]   = 250;
         greenbit[j][i] = 250;
         bluebit[j][i]  = 250;
      };


  /* write ppm file */

#ifdef MPI
  if (0==myid) {
#endif

  out = fopen(fname,"w");
  if (NULL == out) error("Can`t open bitmap file.");

  fprintf(out,"P6 %d %d 255\n", XRES, YRES);


  for (j=0; j<YRES; j++ ) {
    for (i=0; i<XRES; i++ ) {
      buf[3*i  ] = (char) redbit[j][i];
      buf[3*i+1] = (char) greenbit[j][i];
      buf[3*i+2] = (char) bluebit[j][i];
    };
    fwrite(buf, sizeof(char), 3*XRES, out);
  };

  fclose(out);

#ifdef MPI
 };
#endif

/* Potential energy second */


  /* create filename */

  sprintf(fname,"%s.%u.pot.ppm",outfilename,fzhlr);


  /* Zero bitmap */
  for (i=0; i<YRES; i++ ) {
    for (j=0; j<XRES; j++ ) {
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

                 pix =  redbit  [pixcoord.y][pixcoord.x] + (SFACTOR * 255 * red  );
                 redbit  [pixcoord.y][pixcoord.x] = (shortint) pix < 255 ? pix : 255;
	      
                 pix = bluebit [pixcoord.y][pixcoord.x] + (SFACTOR * 255 * blue );
                 bluebit [pixcoord.y][pixcoord.x] = (shortint) pix < 255 ? pix : 255;

                 pix = greenbit[pixcoord.y][pixcoord.x] + (SFACTOR * 255 * green);
                 greenbit[pixcoord.y][pixcoord.x] = (shortint) pix < 255 ? pix : 255;
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
   for (j=0; j<YRES; j++ ) {
     for (i=0; i<XRES; i++ ) { 

      redbit[j][i]   = sum_red[j][i]   < 255 ? sum_red[j][i]   : 200;
      greenbit[j][i] = sum_green[j][i] < 255 ? sum_green[j][i] : 200;
      bluebit[j][i]  = sum_blue[j][i]  < 255 ? sum_blue[j][i]  : 200;
    };
  }; 

  for (j=0; j<YRES; j++ ) 
    for (i=0; i<XRES; i++ ) 
/* background white */
      if ((0==redbit[j][i]) && (0==greenbit[j][i]) && (0==bluebit[j][i])) {
         redbit[j][i]   = 250;
         greenbit[j][i] = 250;
         bluebit[j][i]  = 250;
      };
  
};
#endif

  /* write ppm file */

#ifdef MPI
  if (0==myid) {
#endif

  out = fopen(fname,"w");
  if (NULL == out) error("Can`t open bitmap file.");

  fprintf(out,"P6 %d %d 255\n", XRES, YRES);


  for (j=0; j<XRES; j++ ) {
    for (i=0; i<YRES; i++ ) {
      buf[3*i  ] = (char) redbit[j][i];
      buf[3*i+1] = (char) greenbit[j][i];
      buf[3*i+2] = (char) bluebit[j][i];
    };
    fwrite(buf, sizeof(char), 3*XRES, out);
  };

  fclose(out);

#ifdef MPI
 };
#endif
}



/******************************************************************************
*
* write_pictures writes bitmap pictures of configuration
*
******************************************************************************/

void write_pictures_cooked(int steps)

#ifdef XRES
#undef XRES
#endif
#ifdef YRES
#undef YRES
#endif
#ifdef NUMPIX
#undef NUMPIX
#endif
#define XRES 720
#define YRES 576
#define SFACTOR 1.0
#define NUMPIX  8

{

  vektor scale;

  shortint redbit  [YRES][XRES];
  shortint greenbit[YRES][XRES];
  shortint bluebit [YRES][XRES];

#ifdef MPI
  shortint sum_red[XRES][YRES];
  shortint sum_green[XRES][YRES];
  shortint sum_blue[XRES][YRES];
  ivektor2d maxcoord,mincoord;
#endif

  char buf[3*XRES];

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

  /* compute pic_scale */
  pic_scale.x = box_x.x / (conf_ur.x - conf_ll.x);
  pic_scale.y = box_y.y / (conf_ur.y - conf_ll.y);

  /* Make scaling same in both directions */
  if (scale.y<scale.x)  scale.x = scale.y;

  scale.x *= pic_scale.x;
  scale.y *= pic_scale.y;

  /* kinetic energy first */

  /* create filename */
  /* Dateiname fuer Ausgabedatei erzeugen */
  fzhlr = steps / pic_interval;

  sprintf(fname,"%s.%u.kin.atm",outfilename,fzhlr);

  /* Zero bitmap */
  for (j=0; j<YRES; j++ ) {
    for (i=0; i<XRES; i++ ) {
      bluebit[j][i]  = 0;
      greenbit[j][i] = 0;
      redbit[j][i]   = 0;
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
          if ( (p->ort X(i) < conf_ll.x) || (p->ort X(i) > conf_ur.x) ||
             (p->ort Y(i) < conf_ll.y) || (p->ort Y(i) > conf_ur.y) )
            continue;
	  coord.x = (int) (p->ort X(i) * scale.x);
	  coord.y = (int) (p->ort Y(i) * scale.y);
	  /* Check bounds */
	  if ((coord.x >= NUMPIX) & (coord.x < (XRES-NUMPIX)) &&
              (coord.y >= NUMPIX) && (coord.y < (YRES-NUMPIX))) { 

             coord.y = YRES - coord.y; /* in pic: from top to bottom */
             val = SPRODN(p->impuls,i,p->impuls,i) / (2*p->masse[i]);

             /* Scale Value to [0..1]   */
	     val = (val - ecut_kin.x) / (ecut_kin.y - ecut_kin.x);
             val = val > 1.0 ? 1.0 : val;
             val = val < 0.0 ? 0.0 : val;
             np  = NUMPIX * (1 + p->sorte[i]);
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
             for (j=-np;j<np;++j)
               for (k=-np;k<np;++k) {
                 if (k*k + j*j > np*np) continue;
                 pixcoord.x = coord.x + j;
                 pixcoord.y = coord.y + k;

                 pix =  redbit  [pixcoord.y][pixcoord.x] + (SFACTOR * 255 * red  );
                 redbit  [pixcoord.y][pixcoord.x] = (shortint) pix < 255 ? pix : 255;
	      
                 pix = bluebit [pixcoord.y][pixcoord.x] + (SFACTOR * 255 * blue );
                 bluebit [pixcoord.y][pixcoord.x] = (shortint) pix < 255 ? pix : 255;

                 pix = greenbit[pixcoord.y][pixcoord.x] + (SFACTOR * 255 * green);
                 greenbit[pixcoord.y][pixcoord.x] = (shortint) pix < 255 ? pix : 255;

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
   for (j=0; j<YRES; j++ ) {
     for (i=0; i<XRES; i++ ) { 

      redbit[j][i]   = sum_red[j][i]   < 255 ? sum_red[j][i]   : 200;
      greenbit[j][i] = sum_green[j][i] < 255 ? sum_green[j][i] : 200;
      bluebit[j][i]  = sum_blue[j][i]  < 255 ? sum_blue[j][i]  : 200;
    };
  }; 

  
};
#endif

  for (j=0; j<YRES; j++ ) 
    for (i=0; i<XRES; i++ ) 
/* background white */
      if ((0==redbit[j][i]) && (0==greenbit[j][i]) && (0==bluebit[j][i])) {
         redbit[j][i]   = 250;
         greenbit[j][i] = 250;
         bluebit[j][i]  = 250;
      };


  /* write ppm file */

#ifdef MPI
  if (0==myid) {
#endif

  out = fopen(fname,"w");
  if (NULL == out) error("Can`t open bitmap file.");

  fprintf(out,"P6 %d %d 255\n", XRES, YRES);


  for (j=0; j<YRES; j++ ) {
    for (i=0; i<XRES; i++ ) {
      buf[3*i  ] = (char) redbit[j][i];
      buf[3*i+1] = (char) greenbit[j][i];
      buf[3*i+2] = (char) bluebit[j][i];
    };
    fwrite(buf, sizeof(char), 3*XRES, out);
  };

  fclose(out);

#ifdef MPI
 };
#endif


/* Potential energy second */


  /* create filename */

  sprintf(fname,"%s.%u.pot.atm",outfilename,fzhlr);

  /* Zero bitmap */
  for (j=0; j<YRES; j++ ) {
    for (i=0; i<XRES; i++ ) {
      bluebit[j][i]  = 0;
      greenbit[j][i] = 0;
      redbit[j][i]   = 0;
    };
  };

  /* loop over all atoms */
  for ( r = cellmin.x; r < cellmax.x; ++r )
    for ( s = cellmin.y; s < cellmax.y; ++s )
      {

	p = PTR_2D_V(cell_array, r, s, cell_dim);

	for (i = 0;i < p->n; ++i) {
          if ( (p->ort X(i) < conf_ll.x) || (p->ort X(i) > conf_ur.x) ||
             (p->ort Y(i) < conf_ll.y) || (p->ort Y(i) > conf_ur.y) )
            continue;
	  coord.x = (int) (p->ort X(i) * scale.x);
	  coord.y = (int) (p->ort Y(i) * scale.y);
	  /* Check bounds */
	  if ((coord.x>=NUMPIX) && (coord.x<(XRES-NUMPIX)) &&
	      (coord.y>=NUMPIX) && (coord.y<(YRES-NUMPIX))) {

          coord.y = YRES - coord.y;
#ifdef DISLOC
          if (Epot_diff==1) {
            val = p->pot_eng[i] - p->Epot_ref[i];
          } else
#else
            val = p->pot_eng[i];
#endif

          /* Scale Value to [0..1]   */
	  val = (val - ecut_pot.x) / (ecut_pot.y - ecut_pot.x);
/* Values that are not in the interval are set to MINIMUM */
          val = val > 1.0 ? 0.0 : val;
          val = val < 0.0 ? 0.0 : val;
          np  = NUMPIX * (1 + p->sorte[i]);
/* Defects in potential energy are rather point-like, so we enlarge all
Pixels not in the default interval
          if ((val<1.0) && (val>0.0)) np = 3 * NUMPIX; else np = NUMPIX;
*/
          /* Get index into table */
	  ind = (int)(val * 3.9999);

	  red = -4.0 * (tabred[ind] - tabred[ind+1]) * val + 
	    (tabred[ind] * (ind+1) - ind * tabred[ind+1] );

	  green = -4.0 * (tabgreen[ind] - tabgreen[ind+1]) * val + 
	    (tabgreen[ind] * (ind+1) - ind * tabgreen[ind+1] );

          blue = -4.0 * (tabblue[ind] - tabblue[ind+1]) * val + 
	    (tabblue[ind] * (ind+1) - ind * tabblue[ind+1] );


          /* Set & Copy Pixel */
          for (j=-np;j<np;++j)
             for (k=-np;k<np;++k) {
               if (k*k + j*j > np*np) continue;
               pixcoord.x = coord.x + j;
               pixcoord.y = coord.y + k;

/*               if ((pixcoord.x<XRES) && (pixcoord.y<YRES)) {
*/
                 pix =  redbit  [pixcoord.y][pixcoord.x] + (SFACTOR * 255 * red  );
                 redbit  [pixcoord.y][pixcoord.x] = (shortint) pix < 255 ? pix : 255;
	      
                 pix = bluebit [pixcoord.y][pixcoord.x] + (SFACTOR * 255 * blue );
                 bluebit [pixcoord.y][pixcoord.x] = (shortint) pix < 255 ? pix : 255;

                 pix = greenbit[pixcoord.y][pixcoord.x] + (SFACTOR * 255 * green);
                 greenbit[pixcoord.y][pixcoord.x] = (shortint) pix < 255 ? pix : 255;
/*              };
*/           }; 
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
   for (j=0; j<YRES; j++ ) {
     for (i=0; i<XRES; i++ ) { 

      redbit[j][i]   = sum_red[j][i]   < 255 ? sum_red[j][i]   : 200;
      greenbit[j][i] = sum_green[j][i] < 255 ? sum_green[j][i] : 200;
      bluebit[j][i]  = sum_blue[j][i]  < 255 ? sum_blue[j][i]  : 200;
    };
  }; 

  for (j=0; j<YRES; j++ ) 
    for (i=0; i<XRES; i++ ) 
/* background white */
      if ((0==redbit[j][i]) && (0==greenbit[j][i]) && (0==bluebit[j][i])) {
         redbit[j][i]   = 250;
         greenbit[j][i] = 250;
         bluebit[j][i]  = 250;
      };
  
};
#endif

  /* write ppm file */

#ifdef MPI
  if (0==myid) {
#endif

  out = fopen(fname,"w");
  if (NULL == out) error("Can`t open bitmap file.");

  fprintf(out,"P6 %d %d 255\n", XRES, YRES);


  for (j=0; j<YRES; j++ ) {
    for (i=0; i<XRES; i++ ) {
      buf[3*i  ] = (char) redbit[j][i];
      buf[3*i+1] = (char) greenbit[j][i];
      buf[3*i+2] = (char) bluebit[j][i];
    };
    fwrite(buf, sizeof(char), 3*XRES, out);
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
    case 1: write_pictures_cooked ( steps ); break;
    case 2: write_pictures_bins( steps ); break;
  };
}
