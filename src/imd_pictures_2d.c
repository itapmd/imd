
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2008 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"


/******************************************************************************
*
* write_pictures writes bitmap pictures of configuration
*
******************************************************************************/

void write_pictures_bitmap(int steps)

#define SFACTOR 1.0

{
  vektor scale;

  static shortint *redbit    = NULL;
  static shortint *greenbit  = NULL;
  static shortint *bluebit   = NULL;

#ifdef MPI
  static shortint *sum_red   = NULL;
  static shortint *sum_green = NULL;
  static shortint *sum_blue  = NULL;
  ivektor2d maxcoord,mincoord;
#endif

  static unsigned char *buf  = NULL;
  str255 fname;
  int fzhlr;
  int i,j,k,l,r,s;
  real val;
  cell *p;
  real red,green,blue;
  FILE *out;
  ivektor2d coord,pixcoord;
  real tabred[5],tabgreen[5],tabblue[5];
  int ind;
  int pix;
  int np = nsmear;
  int size, ii;

  size = pic_res.x * pic_res.y;
#ifdef MPI
  sum_red   = (shortint*) realloc(sum_red,   size*sizeof(shortint));
  sum_green = (shortint*) realloc(sum_green, size*sizeof(shortint));
  sum_blue  = (shortint*) realloc(sum_blue,  size*sizeof(shortint));
#endif
  redbit    = (shortint*) realloc(redbit,    size*sizeof(shortint));
  greenbit  = (shortint*) realloc(greenbit,  size*sizeof(shortint));
  bluebit   = (shortint*) realloc(bluebit,   size*sizeof(shortint));
  buf       = (unsigned char*) realloc(buf, 3*pic_res.x);

  /* the dist bins are orthogonal boxes in space */
  scale.x = pic_res.x / box_x.x;
  scale.y = pic_res.y / box_y.y;

  /* compute pic_scale */
  pic_scale.x = box_x.x / (pic_ur.x - pic_ll.x);
  pic_scale.y = box_y.y / (pic_ur.y - pic_ll.y);

  scale.x *= pic_scale.x;
  scale.y *= pic_scale.y;

  /* kinetic energy first */

  /* create filename */
  /* Dateiname fuer Ausgabedatei erzeugen */
  fzhlr = steps / pic_int;

  sprintf(fname,"%s.%u.kin.ppm",outfilename,fzhlr);

  /* Zero bitmap */
  for (i=0; i<size; i++ ) {
    bluebit [i] = 0;
    greenbit[i] = 0;
    redbit  [i] = 0;
  }

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
  for (k=0; k<ncells; k++) {

    p = cell_array + CELLS(k);

    for (i=0; i<p->n; ++i) {

      coord.x = (int) ((ORT(p,i,X)-pic_ll.x) * scale.x);
      coord.y = (int) ((ORT(p,i,Y)-pic_ll.y) * scale.y);

      /* Check bounds */
      if ((coord.x >= 0) && (coord.x < pic_res.x) &&
          (coord.y >= 0) && (coord.y < pic_res.y)) { 

        val = SPRODN(IMPULS,p,i,IMPULS,p,i) / (2*MASSE(p,i));

        /* Scale Value to [0..1]   */
        val = (val - ecut_kin.x) / (ecut_kin.y - ecut_kin.x);
        val = val > 1.0 ? 1.0 : val;
        val = val < 0.0 ? 0.0 : val;
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
        for (j=-np; j<np; ++j)
          for (l=-np; l<np; ++l) {
            if (l*l + j*j > np*np) continue;
            if ((j<0) || (l<0) || (j>=pic_res.x) || (l>=pic_res.y)) continue; 
            ii = (coord.y+l) * pic_res.x + coord.x+j;
            pix =  redbit [ii] + (SFACTOR * 255 * red  );
            redbit  [ii] = (shortint) pix < 255 ? pix : 255;
            pix = bluebit [ii] + (SFACTOR * 255 * blue );
            bluebit [ii] = (shortint) pix < 255 ? pix : 255;
            pix = greenbit[ii] + (SFACTOR * 255 * green);
            greenbit[ii] = (shortint) pix < 255 ? pix : 255;
	  }
      }
    }
  }

#ifdef MPI
  /* Add the bitmaps */
  MPI_Reduce( redbit,   sum_red,   size, MPI_SHORT, MPI_SUM, 0, cpugrid);
  MPI_Reduce( greenbit, sum_green, size, MPI_SHORT, MPI_SUM, 0, cpugrid);
  MPI_Reduce( bluebit , sum_blue,  size, MPI_SHORT, MPI_SUM, 0, cpugrid);

  /* Clip at max value of bitmap */
  if (0==myid)
    for (i=0; i<size; i++) {
      redbit  [i] = sum_red  [i] < 255 ? sum_red  [i] : 200;
      greenbit[i] = sum_green[i] < 255 ? sum_green[i] : 200;
      bluebit [i] = sum_blue [i] < 255 ? sum_blue [i] : 200;
    }
#endif

  /* background white */
  if (0==myid)
    for (i=0; i<size; i++)
      if ((0==redbit[i]) && (0==greenbit[i]) && (0==bluebit[i])) {
        redbit  [i] = 250;
        greenbit[i] = 250;
        bluebit [i] = 250;
      }

  /* write ppm file */
  if (0==myid) {

    out = fopen(fname,"w");
    if (NULL == out) error("Cannot open bitmap file.");

    fprintf(out,"P6 %d %d 255\n", pic_res.x, pic_res.y);

    for (j=pic_res.y-1; j>=0; j-- ) {
      for (i=0; i<pic_res.x; i++ ) {
        ii = j * pic_res.x + i;
        buf[3*i  ] = (unsigned char) redbit  [ii];
        buf[3*i+1] = (unsigned char) greenbit[ii];
        buf[3*i+2] = (unsigned char) bluebit [ii];
      }
      fwrite(buf, 1, 3*pic_res.x, out);
    }
    fclose(out);
  }

  /* Potential energy second */
  /* create filename */

  sprintf(fname,"%s.%u.pot.ppm",outfilename,fzhlr);

  /* Zero bitmap */
  for (i=0; i<size; i++) {
    bluebit [i] = 0;
    greenbit[i] = 0;
    redbit  [i] = 0;
  }

  /* loop over all atoms */
  for (k=0; k<ncells; k++) {

    p = cell_array + CELLS(k);

    for (i=0; i<p->n; ++i) {

      coord.x = (int) ((ORT(p,i,X)-pic_ll.x) * scale.x);
      coord.y = (int) ((ORT(p,i,Y)-pic_ll.y) * scale.y);

      /* Check bounds */
      if ((coord.x>=0) && (coord.x<pic_res.x) &&
          (coord.y>=0) && (coord.y<pic_res.y)) {

        val = POTENG(p,i);

        /* Scale Value to [0..1]   */
        val = (val - ecut_pot.x) / (ecut_pot.y - ecut_pot.x);
        /* Values that are not in the interval are set to MINIMUM */
        val = val > 1.0 ? 0.0 : val;
        val = val < 0.0 ? 0.0 : val;
        /* Get index into table */
	ind = (int)(val * 3.9999);

	red = -4.0 * (tabred[ind] - tabred[ind+1]) * val + 
	    (tabred[ind] * (ind+1) - ind * tabred[ind+1] );

	green = -4.0 * (tabgreen[ind] - tabgreen[ind+1]) * val + 
	    (tabgreen[ind] * (ind+1) - ind * tabgreen[ind+1] );

        blue = -4.0 * (tabblue[ind] - tabblue[ind+1]) * val + 
	    (tabblue[ind] * (ind+1) - ind * tabblue[ind+1] );

        /* Set & Copy Pixel */
        for (j=-np; j<np; ++j)
          for (l=-np; l<np; ++l) {
            if (l*l + j*j > np*np) continue;
            if ((j<0) || (l<0) || (j>=pic_res.x) || (l>=pic_res.y)) continue; 
            ii = (coord.y+l) * pic_res.x + coord.x+j;
            pix =  redbit[ii]  + (SFACTOR * 255 * red  );
            redbit  [ii] = (shortint) pix < 255 ? pix : 255;
            pix = bluebit[ii]  + (SFACTOR * 255 * blue );
            bluebit [ii] = (shortint) pix < 255 ? pix : 255;
            pix = greenbit[ii] + (SFACTOR * 255 * green);
            greenbit[ii] = (shortint) pix < 255 ? pix : 255;
	  }
      }
    }
  }

#ifdef MPI
  /* Add the bitmaps */
  MPI_Reduce( redbit,   sum_red,   size, MPI_SHORT, MPI_SUM, 0, cpugrid);
  MPI_Reduce( greenbit, sum_green, size, MPI_SHORT, MPI_SUM, 0, cpugrid);   
  MPI_Reduce( bluebit , sum_blue,  size, MPI_SHORT, MPI_SUM, 0, cpugrid);

  /* Clip at max value of bitmap */
  if (0==myid)
    for (i=0; i<size; i++ ) { 
      redbit  [i] = sum_red  [i] < 255 ? sum_red  [i] : 200;
      greenbit[i] = sum_green[i] < 255 ? sum_green[i] : 200;
      bluebit [i] = sum_blue [i] < 255 ? sum_blue [i] : 200;
    }
#endif

  /* background white */
  if (0==myid)
    for (i=0; i<size; i++)
      if ((0==redbit[i]) && (0==greenbit[i]) && (0==bluebit[i])) {
        redbit  [i] = 250;
        greenbit[i] = 250;
        bluebit [i] = 250;
      }

  /* write ppm file */
  if (0==myid) {

    out = fopen(fname,"w");
    if (NULL == out) error("Cannot open bitmap file.");

    fprintf(out,"P6 %d %d 255\n", pic_res.x, pic_res.y);

    for (j=pic_res.y-1; j>=0; j-- ) {
      for (i=0; i<pic_res.x; i++ ) {
        ii = j * pic_res.x + i;
        buf[3*i  ] = (unsigned char) redbit  [ii];
        buf[3*i+1] = (unsigned char) greenbit[ii];
        buf[3*i+2] = (unsigned char) bluebit [ii];
      }
      fwrite(buf, 1, 3*pic_res.x, out);
    }
    fclose(out);
  }
}

void write_pictures(int steps)
{
  switch (pic_type) {
    case 0: 
      write_config_select(steps/pic_int, "pic",
                          write_atoms_pic, write_header_pic); 
      break;
    case 1: 
      write_pictures_bitmap(steps); break;
  }
}
