/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

void write_pic_cell( cell *p, FILE *out ) 
{
    struct { 
        float    pos_x, pos_y, pos_z, E_kin, E_pot;
        integer  type;
    } picbuf;

    real px, py, pz;
    int i;

    for (i = 0;i < p->n; ++i) {
        picbuf.pos_x = (float) p->ort X(i);
        picbuf.pos_y = (float) p->ort Y(i);
        picbuf.pos_z = (float) p->ort Z(i);
        if ( pic_ur.x != (real)0 ) /*if pic_ur still 0, write everything */
        if ( (picbuf.pos_x < pic_ll.x) || (picbuf.pos_x > pic_ur.x) ||
             (picbuf.pos_y < pic_ll.y) || (picbuf.pos_y > pic_ur.y) ||
             (picbuf.pos_z < pic_ll.z) || (picbuf.pos_z > pic_ur.z) ) continue;
        px = p->impuls X(i);
        py = p->impuls Y(i);
        pz = p->impuls Z(i);
        picbuf.E_kin = (float) ( (px*px + py*py + pz*pz) / (2 * p->masse[i]) );
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


void write_vrml_cell( cell *p, FILE *out ) 
{

    real px, py, pz, E_pot, E_kin;
    real red, green, blue;
    int i, ind;
    real tabred[5],tabgreen[5],tabblue[5];

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
  
    for (i = 0;i < p->n; ++i) {
        px = p->impuls X(i);
        py = p->impuls Y(i);
        pz = p->impuls Z(i);
        E_kin = (float) ( (px*px + py*py + pz*pz) / (2 * p->masse[i]) );
#ifdef DISLOC
        if (Epot_diff==1) {
          E_pot = (float) p->pot_eng[i] - p->Epot_ref[i];
        } else
#else
        E_pot = (float) (p->pot_eng[i]+35)/82;
#endif
	ind = (int)(E_pot * 3.9999);

        /* Get RBG values from linear interpolation */
	red = -4.0 * (tabred[ind] - tabred[ind+1]) * E_pot + (tabred[ind] * (ind+1) - ind * tabred[ind+1] );

	green = -4.0 * (tabgreen[ind] - tabgreen[ind+1]) * E_pot + (tabgreen[ind] * (ind+1) - ind * tabgreen[ind+1] );

        blue = -4.0 * (tabblue[ind] - tabblue[ind+1]) * E_pot + (tabblue[ind] * (ind+1) - ind * tabblue[ind+1] );

        fprintf(out, "  Separator {\n");
	fprintf(out, "    Material { diffuseColor %f %f %f }\n", red, green, blue);
        fprintf(out, "    Translation { translation %f %f %f }\n", p->ort X(i), p->ort Y(i), p->ort Z(i));
	fprintf(out, "    Sphere { radius %f }\n", (float)((p->sorte[i] + 1)*.3));
	fprintf(out, " }\n"); }
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


    for (i = 1; i < cell_dim.x-1; ++i )
      for (j = 1; j < cell_dim.y-1; ++j )
        for (k = 1; k < cell_dim.z-1; ++k ) {
 	  p = PTR_3D_V(cell_array, i, j, k, cell_dim);
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
    for (i = 1; i < cell_dim.x-1; ++i )
      for (j = 1; j < cell_dim.y-1; ++j )
        for (k = 1; k < cell_dim.z-1; ++k ) {
	  p = PTR_3D_V(cell_array, i, j, k, cell_dim);
	  write_pic_cell( p, out );
	};

      /* Receive data from other cpus and write that */
      p   = PTR_2D_V(cell_array, 0, 0, cell_dim);
      for ( m = 1; m < num_cpus; ++m)
        for (i = 1; i < cell_dim.x-1; ++i )
          for (j = 1; j < cell_dim.y-1; ++j )
            for (k = 1; k < cell_dim.z-1; ++k ) {
	      tag = PTR_3D_V(CELL_TAG, i, j, k, cell_dim);
	      recv_cell( p, m, tag );
	      write_pic_cell( p, out );
	    };

      fclose(out);      
    } else { 
      /* Send data to cpu 0 */
        for (i = 1; i < cell_dim.x-1; ++i )
          for (j = 1; j < cell_dim.y-1; ++j )
            for (k = 1; k < cell_dim.z-1; ++k ) {
      	      p   = PTR_3D_V(cell_array, i, j, k, cell_dim);
	      tag = PTR_3D_V(CELL_TAG, i, j, k, cell_dim);
	      send_cell( p, 0, tag );
	    };
    };
  };
#else

  /* Ausgabedatei oeffnen */
  out = fopen(fname,"w");
  if (NULL == out) error("Can't open output file for config.");


  for (p = cell_array; 
       p <= PTR_3D_V(cell_array,
		     cell_dim.x-1,
		     cell_dim.y-1,
		     cell_dim.z-1,
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
#define NUMPIX  5

{

  vektor scale;

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

  real xshift, yshift;

  str255 fname;
  int fzhlr;
  int i,j,k,r,s,t;
  real phi;
  real val;
  cell *p;
  real red,green,blue;
  real xmin, xmax, ymin, ymax;
  FILE *out;
  vektor3d a, b;
  ivektor3d coord,pixcoord;
  real tabred[5],tabgreen[5],tabblue[5];
  int ind;
  int pix;
  int np;
  int ia, ib;

  /* base vectors normal to view_dir */
  a.x = - view_dir.y;
  a.y = view_dir.x;
  a.z = 0;
  a.x /= sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
  a.y /= sqrt(a.x*a.x + a.y*a.y + a.z*a.z);

  b.x = -view_dir.x*view_dir.z;
  b.y = view_dir.y*view_dir.z;
  b.z = view_dir.x*view_dir.x + view_dir.y*view_dir.y;
  b.x /= sqrt(b.x*b.x + b.y*b.y + b.z*b.z);
  b.y /= sqrt(b.x*b.x + b.y*b.y + b.z*b.z);
  b.z /= sqrt(b.x*b.x + b.y*b.y + b.z*b.z);

  /* minimal and maximal values out of box vectors */
  xmin = 0;
  ymin = 0;
  xmax = box_x.x;
  ymax = box_y.y;
  ia = (int)floor(box_x.x*b.x + box_y.y*b.y + box_z.z*b.z);
  ib = (int)floor(box_x.x*b.x + box_y.y*b.y + box_z.z*b.z);
  if (ia < xmin) xmin = ia;
  if (ib < ymin) ymin = ib;
  if (ia > xmax) xmax = ia;
  if (ib > ymax) ymax = ib;
  ia = (int)floor(box_x.x*a.x + box_y.y*a.y);
  ib = (int)floor(box_x.x*b.x + box_y.y*b.y);
  if (ia < xmin) xmin = ia;
  if (ib < ymin) ymin = ib;
  if (ia > xmax) xmax = ia;
  if (ib > ymax) ymax = ib;
  ia = (int)floor(box_x.x*a.x + box_z.z*a.z);
  ib = (int)floor(box_x.x*b.x + box_z.z*b.z);
  if (ia < xmin) xmin = ia;
  if (ib < ymin) ymin = ib;
  if (ia > xmax) xmax = ia;
  if (ib > ymax) ymax = ib;
  ia = (int)floor(box_y.y*a.y + box_z.z*a.z);
  ib = (int)floor(box_y.y*b.y + box_z.z*b.z);
  if (ia < xmin) xmin = ia;
  if (ib < ymin) ymin = ib;
  if (ia > xmax) xmax = ia;
  if (ib > ymax) ymax = ib;
  ia = (int)floor(box_x.x*a.x);
  ib = (int)floor(box_x.x*b.x);
  if (ia < xmin) xmin = ia;
  if (ib < ymin) ymin = ib;
  if (ia > xmax) xmax = ia;
  if (ib > ymax) ymax = ib;
  ia = (int)floor(box_z.z*a.z);
  ib = (int)floor(box_z.z*b.z);
  if (ia < xmin) xmin = ia;
  if (ib < ymin) ymin = ib;
  if (ia > xmax) xmax = ia;
  if (ib > ymax) ymax = ib;
  
  /* the dist bins are orthogonal boxes in space */

  xshift = -xmin;
  yshift = -ymin;
  scale.x = XRES / (xmax-xmin);
  scale.y = YRES / (ymax-ymin);

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
      for ( s = cellmin.y; s < cellmax.y; ++s )
        for ( t = cellmin.z; t < cellmax.z; ++t ) {
      
   	  p = PTR_3D_V(cell_array, r, s, t, cell_dim);

	  for (i = 0;i < p->n; ++i) {
          if ( (p->ort X(i) < conf_llf.x) || (p->ort X(i) > conf_urb.x) ||
               (p->ort Y(i) < conf_llf.y) || (p->ort Y(i) > conf_urb.y) ||
               (p->ort Z(i) < conf_llf.z) || (p->ort Z(i) > conf_urb.z) )
            continue;
          coord.x = (int)floor((p->ort X(i)*a.x + p->ort Y(i)*a.y + p->ort Z(i)*a.z + xshift)*scale.x);
          coord.y = (int)floor((p->ort X(i)*b.x + p->ort Y(i)*b.y + p->ort Z(i)*b.z + yshift)*scale.y);
	  /* Check bounds */
	  if ((coord.x >= NUMPIX) && (coord.x < (XRES-NUMPIX)) &&
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

  sprintf(fname,"%s.%u.pot.ppm",outfilename,fzhlr);


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
      for ( t = cellmin.z; t < cellmax.z; ++t )
      {

	p = PTR_3D_V(cell_array, r, s, t, cell_dim);

	for (i = 0;i < p->n; ++i) {
          if ( (p->ort X(i) < conf_llf.x) || (p->ort X(i) > conf_urb.x) ||
               (p->ort Y(i) < conf_llf.y) || (p->ort Y(i) > conf_urb.y) ||
               (p->ort Z(i) < conf_llf.z) || (p->ort Z(i) > conf_urb.z) )
            continue;
          coord.x = (int)floor((p->ort X(i)*a.x + p->ort Y(i)*a.y + p->ort Z(i)*a.z + xshift)*scale.x);
          coord.y = (int)floor((p->ort X(i)*b.x + p->ort Y(i)*b.y + p->ort Z(i)*b.z + yshift)*scale.y);
	  /* Check bounds */
	  if ((coord.x>=NUMPIX) && (coord.x<(XRES-NUMPIX)) &&
	      (coord.y>=NUMPIX) && (coord.y<(YRES-NUMPIX))) {

	  val = p->pot_eng[i];

          /* Scale Value to [0..1]   */
	  val = (val - ecut_pot.x) / (ecut_pot.y - ecut_pot.x);
/* Values that are not in the interval are set to MINIMUM */
          val = val > 1.0 ? 0.0 : val;
          val = val < 0.0 ? 0.0 : val;
          np  = NUMPIX * (1 + p->sorte[i]);
/* Defects in potential energy are rather point-like, so we enlarge all
Pixels not in the default interval
          if ((val<1.0) && (val>0.0)) np = 3 * NUMPIX; else np = NUMPIX; */

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

              if ((pixcoord.y<YRES) && (pixcoord.x<XRES)) {

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

void write_vrmls(int steps)
{ 
  FILE *out;

  str255 fname;
  int fzhlr;
  cell *p,*q;
  int i,j,k,l,m,tag;
  real phi;

  /* Dateiname fuer Ausgabedatei erzeugen */
  fzhlr = steps / pic_interval;
  phi = acos(-view_dir.z / sqrt(view_dir.x*view_dir.x+view_dir.y*view_dir.y+view_dir.z*view_dir.z));
  
#ifdef MPI  
  if (0<parallel_output)
    sprintf(fname,"%s.%u.%u.wrl",outfilename,fzhlr,myid);
  else
#endif
    sprintf(fname,"%s.%u.wrl",outfilename,fzhlr);


#ifdef MPI

  if (0<parallel_output) {

    /* Ausgabedatei oeffnen */
    out = fopen(fname,"w");
    if (NULL == out) error("Can't open output file for config.");
    fprintf(out, "#VRML V1.0 ascii\n\n");
    fprintf(out, "DEF StartPers PerspectiveCamera {\n");
    fprintf(out, "position 0 0 0\n");
    fprintf(out, "orientation 1 1 1 0\n");
    fprintf(out, "focalDistance 10\n");
    fprintf(out, "}\n\n");
    fprintf(out, "Separator {\n");

    for (i = 1; i < cell_dim.x-1; ++i )
      for (j = 1; j < cell_dim.y-1; ++j )
        for (k = 1; k < cell_dim.z-1; ++k ) {
 	  p = PTR_3D_V(cell_array, i, j, k, cell_dim);
	  write_vrml_cell( p, out );
	};
    
    fprintf(out, "}\n");
    fclose(out);

  } else {

    if (0==myid) {

      /* Ausgabedatei oeffnen */
      out = fopen(fname,"w");
      if (NULL == out) error("Can't open output file for config.");
      fprintf(out, "#VRML V1.0 ascii\n\n");
      fprintf(out, "DEF StartPers PerspectiveCamera {\n");
      fprintf(out, "position 2.5 2.5 15\n");
      fprintf(out, "orientation 0 0 1 0\n");
      fprintf(out, "focalDistance 10\n");
      fprintf(out, "}\n\n");
      fprintf(out, "Separator {\n");


      /* Write data on CPU 0 */

      /* Write own data */
    for (i = 1; i < cell_dim.x-1; ++i )
      for (j = 1; j < cell_dim.y-1; ++j )
        for (k = 1; k < cell_dim.z-1; ++k ) {
	  p = PTR_3D_V(cell_array, i, j, k, cell_dim);
	  write_vrml_cell( p, out );
	};

      /* Receive data from other cpus and write that */
      p   = PTR_3D_V(cell_array, 0, 0, cell_dim);
      for ( m = 1; m < num_cpus; ++m)
        for (i = 1; i < cell_dim.x-1; ++i )
          for (j = 1; j < cell_dim.y-1; ++j )
            for (k = 1; k < cell_dim.z-1; ++k ) {
	      tag = PTR_3D_V(CELL_TAG, i, j, k, cell_dim);
	      recv_cell( p, m, tag );
	      write_vrml_cell( p, out );
	    };

      fprintf(out, "}\n");
      fclose(out);      
    } else { 
      /* Send data to cpu 0 */
        for (i = 1; i < cell_dim.x-1; ++i )
          for (j = 1; j < cell_dim.y-1; ++j )
            for (k = 1; k < cell_dim.z-1; ++k ) {
      	      p   = PTR_3D_V(cell_array, i, j, k, cell_dim);
	      tag = PTR_3D_V(CELL_TAG, i, j, k, cell_dim);
	      send_cell( p, 0, tag );
	    };
    };
  };
#else

  /* Ausgabedatei oeffnen */
  out = fopen(fname,"w");
  if (NULL == out) error("Can't open output file for config.");
  fprintf(out, "#VRML V1.0 ascii\n\n");
  fprintf(out, "Separator {\n");
  if (projection)
    fprintf(out, "  PerspectiveCamera {\n");
  else
    fprintf(out, "  OrthographicCamera {\n");
  fprintf(out, "  position %f %f %f\n", view_pos.x, view_pos.y, view_pos.z);
  fprintf(out, "  orientation %f %f 0 %f\n", view_dir.y, -view_dir.x, phi);
  fprintf(out, "  focalDistance 1\n");
  fprintf(out, "  }\n");
  fprintf(out, "  Separator {\n");
  fprintf(out, "    Material { diffuseColor 1 1 1 transparency 0.7 }\n");
  fprintf(out, "    Translation { translation %f %f %f }\n", box_x.x/2, box_y.y/2, box_z.z/2);
  fprintf(out, "    Cube { width %f depth %f height %f }\n", box_x.x, box_y.y, box_z.z);
  fprintf(out, "  }\n");
  fprintf(out, "  Separator {\n");
  fprintf(out, "    Material { diffuseColor 1 0 0 }\n");
  fprintf(out, "    Translation { translation 0 0 %f }\n", box_z.z/2);
  fprintf(out, "    Rotation { rotation 1 0 0 1.5707963 }\n");
  fprintf(out, "    Cylinder { radius .1 height %f}\n", box_z.z*1.1);
  fprintf(out, "  }\n");
  fprintf(out, "  Separator {\n");
  fprintf(out, "    Material { diffuseColor 0 1 0}\n");
  fprintf(out, "    Translation { translation %f 0 0 }\n", box_x.x/2);
  fprintf(out, "    Rotation { rotation 0 0 -1 1.5707963 }\n");
  fprintf(out, "    Cylinder { radius .1 height %f}\n", box_x.x*1.1);
  fprintf(out, "  }\n");
  fprintf(out, "  Separator {\n");
  fprintf(out, "    Material { diffuseColor 0 0 1 }\n");
  fprintf(out, "    Translation { translation 0 %f 0 }\n", box_y.y/2);
  fprintf(out, "    Cylinder { radius .1 height %f}\n", box_y.y*1.1);
  fprintf(out, "  }\n");
  fprintf(out, "  Separator {\n");
  fprintf(out, "    Material { diffuseColor 1 0 0 }\n");
  fprintf(out, "    Translation { translation 0 0 %f }\n", box_z.z*1.1);
  fprintf(out, "    Rotation { rotation 1 0 0 1.5707963 }\n");
  fprintf(out, "    Cone { radius .1 height 1}\n");
  fprintf(out, "  }\n");
  fprintf(out, "  Separator {\n");
  fprintf(out, "    Material { diffuseColor 0 1 0}\n");
  fprintf(out, "    Translation { translation %f 0 0 }\n", box_x.x*1.1);
  fprintf(out, "    Rotation { rotation 0 0 -1 1.5707963 }\n");
  fprintf(out, "    Cone { radius .1 height 1}\n");
  fprintf(out, "  }\n");
  fprintf(out, "  Separator {\n");
  fprintf(out, "    Material { diffuseColor 0 0 1 }\n");
  fprintf(out, "    Translation { translation 0 %f 0 }\n", box_y.y*1.1);
  fprintf(out, "    Cone { radius .1 height 1}\n");
  fprintf(out, "  }\n");
  fprintf(out, "  Separator {\n");
  fprintf(out, "    Material { diffuseColor 1 0 0 }\n");
  fprintf(out, "    Translation { translation 0 0 %f }\n", box_z.z*1.1);
  fprintf(out, "    Rotation { rotation 1 0 0 1.5707963 }\n");
  fprintf(out, "    AsciiText { string \"z\"}\n");
  fprintf(out, "  }\n");
  fprintf(out, "  Separator {\n");
  fprintf(out, "    Material { diffuseColor 0 1 0}\n");
  fprintf(out, "    Translation { translation %f 0 0 }\n", box_x.x*1.1);
  fprintf(out, "    Rotation { rotation 0 0 -1 1.5707963 }\n");
  fprintf(out, "    AsciiText { string \"x\"}\n");
  fprintf(out, "  }\n");
  fprintf(out, "  Separator {\n");
  fprintf(out, "    Material { diffuseColor 0 0 1 }\n");
  fprintf(out, "    Translation { translation 0 %f 0 }\n", box_y.y*1.1);
  fprintf(out, "    AsciiText { string \"y\"}\n");
  fprintf(out, "  }\n");




  for (p = cell_array; 
       p <= PTR_3D_V(cell_array,
		     cell_dim.x-1,
		     cell_dim.y-1,
		     cell_dim.z-1,
		     cell_dim);
       ++p ) {
    write_vrml_cell( p , out );
  };

  fprintf(out, "}\n");
  fclose(out);  

#endif

}


void write_pictures( steps )
{
  switch (pic_type) {
    case 0: write_pictures_raw ( steps ); break;
    case 1: break;
    case 2: write_pictures_bins( steps ); break;
    case 3: write_vrmls( steps ); break;
  };
}
