
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_ps_main -- Main routines of imd_ps
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "util.h"
#include <stdlib.h>
#include <time.h>
#include <ctype.h>

#ifndef TWOD

/******************************************************************************
*
*  draw_pictures -- 3d version
*
******************************************************************************/

void draw_picture3d(void) {

  int l, k, j, i, type, neigh_num, neigh_type;
  cell *p, *neigh_cell;
  vektor pos, npos;
  neightab *neigh;
  FILE *out;
  str255 fname;
  int llx=30, lly=40, urx=570, ury=750;
  vektor2d nullvek = { 0.0, 0.0 };
  real ratio, boxratio;
  real bondlength;
  vektor deltapos, unit, ns;
  real ratom, tmp;
  vektor2d origin, center, ort;
  real sclfct, offset = 20;
  vektor2d width, cornerps[8], scornerps[8], rcornerps[8], tmpv;
  vektor pnt;
  vektor2d nort, deltaort, unitps, nsps, normalps, trans;
  real atomdistance, rbond, nrbond, nratom;
  real angle, deltaangle;
  real radscl, atom_vol;
  struct { vektor pos; int ngh[3], ord; real dist; } rcorner[8], scorner[8];
  int scrnorder[8], rcrnorder[8], tmpi;
  int tmpiv[3];
  vektor2d maxim={Max,Max}, minim={Min,Min};
  vektor2d smaxim={Max,Max}, sminim={Min,Min};
  vektor2d rmaxim={Max,Max}, rminim={Min,Min};
  vektor2d maximage={Max,Max}, minimage={Min,Min};
  real pts;
  real spect_height = 72, spect_offset = 0.0;
  vektor2d size;
  List_elmt *listelmt;
  int corner1 = 4, corner2 = 4;
  vektor2d boxdist, zentrum;
  real ellps;
  vektor inters;
  vektor2d intersps;
  real interspsdist;
  vektor ninters;
  vektor2d nintersps;
  real ninterspsdist;
  real   corr, hordist, spsdist, nspsdist, brbond, tmp2;
  vektor2d urspr;
  real distance, ndistance, nreldistance, nsdistance;
  real distancep, ndistancep;
  vektor ndist, nsdist;
  real asize, adelta, ax, arad, ashade, atmp;
  int N;
  real d_rbond, irbond, mrbond;
  real mrbonddist2;
  time_t now;
  char* revision = "$Revision$";
  char version[10] = "         ";

  /* Estimate atom and bond radii */
  if ( wireframe < 0.0 ) {
    if ( r_atom == -1.0 ) {
      atom_vol = ( maxl.x - minl.x ) * ( maxl.y - minl.y ) 
	* ( maxl.z - minl.z ) / natoms;
      r_atom = pow ( atom_vol / 100.0 , 1.0 / 3.0 );
      printf("Largest atom radius (R): %.2f (estimated)\n", r_atom);
    }
    if ( r_bond == -1.0 ) {
      r_bond = r_atom / 6.0;
      printf("Largest bond radius (B): %.2f (estimated)\n", r_bond);    
    }
  }

  /* Adjust bond radii, r_atom must be greater than r_bond */
  if ( r_bond > r_atom )
    r_atom = r_bond;

  /* Take atom radii into account */
  maxl.x += r_atom;
  minl.x -= r_atom;
  maxl.y += r_atom;
  minl.y -= r_atom;

  real_maxl.x += r_atom;
  real_minl.x -= r_atom;
  real_maxl.y += r_atom;
  real_minl.y -= r_atom;  
  real_maxl.z += r_atom;
  real_minl.z -= r_atom;  

  /* Corners of bounding box */
  tmpv.x = ( maxl.x - minl.x ) / 2.0;
  tmpv.y = ( maxl.y - minl.y ) / 2.0;
 
  pnt.x = -tmpv.x;
  pnt.y = -tmpv.y;
  pnt.z = 0.0;
  cornerps[0] = proj(pnt, 1.0, nullvek); 

  pnt.x = tmpv.x;
  cornerps[1] = proj(pnt, 1.0, nullvek);
 
  pnt.y = tmpv.y;
  cornerps[2] = proj(pnt, 1.0, nullvek);

  pnt.x = -tmpv.x;
  cornerps[3] = proj(pnt, 1.0, nullvek);

  pnt.y = -tmpv.y;
  pnt.z = - ( maxl.z - minl.z );
  cornerps[4] = proj(pnt, 1.0, nullvek); 

  pnt.x = tmpv.x;
  cornerps[5] = proj(pnt, 1.0, nullvek);
 
  pnt.y = tmpv.y;
  cornerps[6] = proj(pnt, 1.0, nullvek);

  pnt.x = -tmpv.x;
  cornerps[7] = proj(pnt, 1.0, nullvek);

  for ( i=0; i<8; i++ ) {
    maxim.x = MAX( maxim.x, cornerps[i].x );
    maxim.y = MAX( maxim.y, cornerps[i].y );
    minim.x = MIN( minim.x, cornerps[i].x );
    minim.y = MIN( minim.y, cornerps[i].y );
  }

  /* Corners of simulation box */
  if ( sframe >= 0 && sframe_linew != 0.0 ) {
    for ( i=0; i<8; i++ )
      scrnorder[i] = i;

    pnt.x = - ( maxl.x + minl.x ) / 2.0;
    pnt.y = - ( maxl.y + minl.y ) / 2.0;
    pnt.z = - maxl.z;
    scorner[0].pos = pnt;
    scorner[0].ngh[0] = 1;    
    scorner[0].ngh[1] = 3;    
    scorner[0].ngh[2] = 4;
    scorner[0].dist = fdistance(scorner[0].pos);

    pnt.x += box_x.x;
    pnt.y += box_x.y;
    pnt.z += box_x.z;
    scorner[1].pos = pnt;
    scorner[1].ngh[0] = 0;    
    scorner[1].ngh[1] = 2;    
    scorner[1].ngh[2] = 5;    
    scorner[1].dist = fdistance(scorner[1].pos);

    pnt.x += box_y.x;
    pnt.y += box_y.y;
    pnt.z += box_y.z;
    scorner[2].pos = pnt;
    scorner[2].ngh[0] = 1;    
    scorner[2].ngh[1] = 3;    
    scorner[2].ngh[2] = 6;    
    scorner[2].dist = fdistance(scorner[2].pos);

    pnt.x -= box_x.x;
    pnt.y -= box_x.y;
    pnt.z -= box_x.z;
    scorner[3].pos = pnt;
    scorner[3].ngh[0] = 0;    
    scorner[3].ngh[1] = 2;    
    scorner[3].ngh[2] = 7;    
    scorner[3].dist = fdistance(scorner[3].pos);

    pnt.x -= box_y.x;
    pnt.y -= box_y.y;
    pnt.z -= box_y.z;
    pnt.x += box_z.x;
    pnt.y += box_z.y;
    pnt.z += box_z.z;
    scorner[4].pos = pnt;
    scorner[4].ngh[0] = 0;    
    scorner[4].ngh[1] = 5;    
    scorner[4].ngh[2] = 7;    
    scorner[4].dist = fdistance(scorner[4].pos);

    pnt.x += box_x.x;
    pnt.y += box_x.y;
    pnt.z += box_x.z;
    scorner[5].pos = pnt;
    scorner[5].ngh[0] = 1;    
    scorner[5].ngh[1] = 4;    
    scorner[5].ngh[2] = 6;    
    scorner[5].dist = fdistance(scorner[5].pos);

    pnt.x += box_y.x;
    pnt.y += box_y.y;
    pnt.z += box_y.z;
    scorner[6].pos = pnt;
    scorner[6].ngh[0] = 2;    
    scorner[6].ngh[1] = 5;    
    scorner[6].ngh[2] = 7;    
    scorner[6].dist = fdistance(scorner[6].pos);

    pnt.x -= box_x.x;
    pnt.y -= box_x.y;
    pnt.z -= box_x.z;
    scorner[7].pos = pnt;
    scorner[7].ngh[0] = 3;    
    scorner[7].ngh[1] = 4;    
    scorner[7].ngh[2] = 6;    
    scorner[7].dist = fdistance(scorner[7].pos);

    for ( i=0; i<7; i++ )
      for ( j=0; j<7; j++ ) 
	if ( scorner[scrnorder[j]].dist > scorner[scrnorder[j+1]].dist ) {
	  tmp            = scrnorder[j+1];
	  scrnorder[j+1] = scrnorder[j]; 
	  scrnorder[j]   = tmp;
	}

    for (i=0; i<8; i++ ) 
      scornerps[i] = proj(scorner[i].pos, 1.0, nullvek);

    for ( i=0; i<8; i++ ) {
      smaxim.x = MAX( smaxim.x, scornerps[i].x );
      smaxim.y = MAX( smaxim.y, scornerps[i].y );
      sminim.x = MIN( sminim.x, scornerps[i].x );
      sminim.y = MIN( sminim.y, scornerps[i].y );
    }
  }

  /* Corners of real bounding box */
  if ( rframe >= 0 && rframe_linew != 0.0 ) {
    for ( i=0; i<8; i++ )
      rcrnorder[i] = i;

    pnt.x = ursprung.x 
      - ( unitv[0].x + unitv[1].x + unitv[2].x ) * r_atom;
    pnt.y = ursprung.y 
      - ( unitv[0].y + unitv[1].y + unitv[2].y ) * r_atom;
    pnt.z = ursprung.z 
      - ( unitv[0].z + unitv[1].z + unitv[2].z ) * r_atom;

    rcorner[0].pos = pnt;
    rcorner[0].ngh[0] = 1;    
    rcorner[0].ngh[1] = 3;    
    rcorner[0].ngh[2] = 4;
    rcorner[0].dist = fdistance(rcorner[0].pos);
    
    pnt.x += unitv[0].x * ( real_maxl.x - real_minl.x );
    pnt.y += unitv[0].y * ( real_maxl.x - real_minl.x );
    pnt.z += unitv[0].z * ( real_maxl.x - real_minl.x );
    rcorner[1].pos = pnt;
    rcorner[1].ngh[0] = 0;    
    rcorner[1].ngh[1] = 2;    
    rcorner[1].ngh[2] = 5;    
    rcorner[1].dist = fdistance(rcorner[1].pos);
    
    pnt.x += unitv[1].x * ( real_maxl.y - real_minl.y );
    pnt.y += unitv[1].y * ( real_maxl.y - real_minl.y );
    pnt.z += unitv[1].z * ( real_maxl.y - real_minl.y );
    rcorner[2].pos = pnt;
    rcorner[2].ngh[0] = 1;    
    rcorner[2].ngh[1] = 3;    
    rcorner[2].ngh[2] = 6;    
    rcorner[2].dist = fdistance(rcorner[2].pos);

    pnt.x -= unitv[0].x * ( real_maxl.x - real_minl.x );
    pnt.y -= unitv[0].y * ( real_maxl.x - real_minl.x );
    pnt.z -= unitv[0].z * ( real_maxl.x - real_minl.x );
    rcorner[3].pos = pnt;
    rcorner[3].ngh[0] = 0;    
    rcorner[3].ngh[1] = 2;    
    rcorner[3].ngh[2] = 7;    
    rcorner[3].dist = fdistance(rcorner[3].pos);

    pnt.x += unitv[2].x * ( real_maxl.z - real_minl.z );
    pnt.y += unitv[2].y * ( real_maxl.z - real_minl.z );
    pnt.z += unitv[2].z * ( real_maxl.z - real_minl.z );
    rcorner[7].pos = pnt;
    rcorner[7].ngh[0] = 3;    
    rcorner[7].ngh[1] = 4;    
    rcorner[7].ngh[2] = 6;    
    rcorner[7].dist = fdistance(rcorner[7].pos);
   
    pnt.x += unitv[0].x * ( real_maxl.x - real_minl.x );
    pnt.y += unitv[0].y * ( real_maxl.x - real_minl.x );
    pnt.z += unitv[0].z * ( real_maxl.x - real_minl.x );
    rcorner[6].pos = pnt;
    rcorner[6].ngh[0] = 2;    
    rcorner[6].ngh[1] = 5;    
    rcorner[6].ngh[2] = 7;    
    rcorner[6].dist = fdistance(rcorner[6].pos);

    pnt.x -= unitv[1].x * ( real_maxl.y - real_minl.y );
    pnt.y -= unitv[1].y * ( real_maxl.y - real_minl.y );
    pnt.z -= unitv[1].z * ( real_maxl.y - real_minl.y );
    rcorner[5].pos = pnt;
    rcorner[5].ngh[0] = 1;    
    rcorner[5].ngh[1] = 4;    
    rcorner[5].ngh[2] = 6;    
    rcorner[5].dist = fdistance(rcorner[5].pos);
    
    pnt.x -= unitv[0].x * ( real_maxl.x - real_minl.x );
    pnt.y -= unitv[0].y * ( real_maxl.x - real_minl.x );
    pnt.z -= unitv[0].z * ( real_maxl.x - real_minl.x );
    rcorner[4].pos = pnt;
    rcorner[4].ngh[0] = 0;    
    rcorner[4].ngh[1] = 5;    
    rcorner[4].ngh[2] = 7;    
    rcorner[4].dist = fdistance(rcorner[4].pos); 

    for ( i=0; i<7; i++ )
      for ( j=0; j<7; j++ ) 
	if ( rcorner[rcrnorder[j]].dist > rcorner[rcrnorder[j+1]].dist ) {
	  tmp            = rcrnorder[j+1];
	  rcrnorder[j+1] = rcrnorder[j]; 
	  rcrnorder[j]   = tmp;
	}

    for (i=0; i<8; i++ )
      rcornerps[i] = proj(rcorner[i].pos, 1.0, nullvek);

    for ( i=0; i<8; i++ ) {
      rmaxim.x = MAX( rmaxim.x, rcornerps[i].x );
      rmaxim.y = MAX( rmaxim.y, rcornerps[i].y );
      rminim.x = MIN( rminim.x, rcornerps[i].x );
      rminim.y = MIN( rminim.y, rcornerps[i].y );
    }
  }

  /* Determine boundaries of image */
  if ( ( sframe < 0.0 || sframe_linew == 0.0 ) && 
       ( rframe < 0.0 || rframe_linew == 0.0 ) ) {
    maximage = maxim;
    minimage = minim;
  }
  else if ( frame < 0.0 || frame_linew == 0) {
    maximage.x = MAX( smaxim.x, rmaxim.x );
    maximage.y = MAX( smaxim.y, rmaxim.y );
    minimage.x = MIN( sminim.x, rminim.x );
    minimage.y = MIN( sminim.y, rminim.y );
  }
  else {
      maximage.x = MAX( smaxim.x, rmaxim.x );
      maximage.y = MAX( smaxim.y, rmaxim.y );
      minimage.x = MIN( sminim.x, rminim.x );
      minimage.y = MIN( sminim.y, rminim.y );  
      
      maximage.x = MAX( maximage.x, maxim.x );
      maximage.y = MAX( maximage.y, maxim.y );
      minimage.x = MIN( minimage.x, minim.x );
      minimage.y = MIN( minimage.y, minim.y );  
  }

  if ( spect >= 0.0 ) {
      lly += spect_height;
      spect_offset = spect_height;
  }
 
  if ( maximage.x == minimage.x ) {
      maximage.x += 0.5;
      minimage.x -= 0.5;
  }
  if ( maximage.y == minimage.y ) {
      maximage.y += 0.5;
      minimage.y -= 0.5;
  }

  /* Ratio of page (ps)*/
  ratio = (real) ( ury - lly ) / ( urx - llx );
  /* Ratio of image */
  boxratio = ( maximage.y - minimage.y ) / ( maximage.x - minimage.x );
  
  /* Scale factor image -> ps */
  if ( boxratio <= ratio ) 
    sclfct = ( urx - llx ) / ( maximage.x - minimage.x ); 
  else
    sclfct = ( ury - lly ) / ( maximage.y - minimage.y ); 

  /* Width of box (ps) */
  width.x = sclfct * ( maximage.x - minimage.x );
  width.y = sclfct * ( maximage.y - minimage.y );

  /* Diplacement with respect to projected values */
  trans.x = ( urx + llx - sclfct * (  maximage.x + minimage.x ) ) / 2.0;
  trans.y = ( ury + lly - sclfct * (  maximage.y + minimage.y ) ) / 2.0;

  /* Transform bounding box (ps) */
  if ( frame >= 0.0 )
    for ( i=0; i<8; i++ ) 
      cornerps[i] = transform( cornerps[i], sclfct, trans);

  if ( sframe >= 0.0 )
    /* Transform simulation box (ps) */
    for ( i=0; i<8; i++ ) 
      scornerps[i] = transform( scornerps[i], sclfct, trans);

  if ( rframe >= 0.0 ) 
    /* Transform real bounding box (ps) */
    for ( i=0; i<8; i++ ) 
      rcornerps[i] = transform( rcornerps[i], sclfct, trans);

  /* Origin of box (ps) */
  origin.x = ( urx + llx - width.x ) / 2.0;
  origin.y = ( ury + lly - width.y ) / 2.0;

  /* Center of box (ps) */
  center.x = origin.x + width.x / 2.0;
  center.y = origin.y + width.y / 2.0;


  /* Open output file */
  sprintf(fname, "%s.eps", infilename);
  out = fopen(fname,"w");
  if (NULL == out) 
    error("Cannot open eps file.\n");

  /* Crop offset */ 
  if ( crop != -1.0 )
      offset = crop + frame_linew;

  /* Box width and height prescribed in cm */
  if ( ps_width != -1.0 || ps_height != -1.0 ) {

    if ( ps_width == -1.0 ) {
      pts = ps_height / PT - 2 * offset;
      scaling = pts / width.y;
    }
    else if ( ps_width > 0.0 ) {
      pts = ps_width / PT - 2.0 * offset;
      scaling = pts / width.x;
    } 
  } 

  /* Get the creation time */
  time(&now);

  /* Get the version number of imd_ps */
  j = 0;
  for ( i=1; revision[i]!='$'; i++)
      if ( isdigit( revision[i] ) || revision[i] == '.' )
	  version[j++] = revision[i];

  /* Print header*/
  fprintf(out, "%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(out, "%%%%Creator: imd_ps 3D Version %s\n", version);
  fprintf(out, "%%%%CreationDate: %s", ctime(&now));
  fprintf(out, "%%%%For: %s@%s\n", getenv("USER"), getenv("HOST")) ;
  fprintf(out, "%%%%Title: %s.eps\n", infilename); 
  fprintf(out, "%%%%BoundingBox: %d %d %d %d\n", 
	  (int) ( origin.x - offset ), 
	  (int) ( origin.y - scaling * spect_offset - offset ), 
	  (int) ( origin.x + scaling * width.x + offset ),  
	  (int) ( origin.y + scaling * width.y + offset ) ); 
  fprintf(out, "%%%%EndComments\n\n");

  fprintf(out, "gsave\n9 dict begin\n");
  fprintf(out, "\n/Y {grestore} bind def\n");

  /* Atom */
  fprintf(out, "\n/Shade {dup 5 index mul 1 index 5 index mul\n 2 index 5 index mul setrgbcolor} def\n");

  fprintf(out, "\n/Draw {mul 0 360 arc closepath fill} def\n");

  fprintf(out, "\n/A {ucache gsave translate gsave\n newpath 0 0 2 index 0 360 arc\n ");

  if ( shading < 0.0 ) 
    fprintf(out, "3 index 3 index 3 index setrgbcolor closepath fill\n ");
  else {
    fprintf(out, "clip gsave %.2f rotate %.6f 1 scale\n ", 
	    shading, sqrt(2.0/3.0)); 
    
    asize =  r_atom * sclfct;
    N = (int) ( asize / 2.0 * zoom * scaling );

    if ( N<3 )
      N = 3;

    adelta = 1.0 / N;
      
    for ( i=0; i<N; i++ ) {
	
      ax = i * adelta;
      
      atmp = 1.0 - 1.07 * ax * ax;
      if ( atmp < 0.0 )
	break;
      
      arad = sqrt( atmp );

      if ( ambient >= 0.0 )
	  ashade = 0.1 + pow( ax, ambient) * 0.9;
      else
	  ashade = 0.1 + ( 1.0 - pow( ABS(ax-1.0), -ambient)) * 0.9;

      fprintf(out, "%.6f Shade pop dup %.6f mul 0  2 index %.6f Draw\n ", 
	      ashade, ax, arad);
    }    
    fprintf(out, "grestore ");
  }

  /* Draw atom border */
  if ( alinew > 0.0 )
    fprintf(out, "grestore newpath %.6f setlinewidth 0 setgray\n 0 0 2 index 0 360 arc closepath stroke\n ", alinew);
  else
    fprintf(out, "grestore\n ");

  /* Draw lightshadow */
  if ( llinew > 0.0 ) {
    printcolor(out, 0, dist_min, 3);
    fprintf(out, "setrgbcolor newpath 0 0 2 index %.6f add 0 360 arc\n %.6f setlinewidth stroke\n ", (alinew + llinew ) / 2.0, llinew);
  }
  fprintf(out, "pop pop pop pop grestore } bind def\n");

  /* Lower part of bond */ 
  fprintf(out, "\n/B1 {ucache gsave newpath rotate 2 index 5 index neg moveto\n 2 copy neg lineto lineto 0 translate\n 1 scale 0 0 3 -1 roll 90 270 arc\n closepath gsave fill grestore ");
  if ( blinew > 0.0 ) 
    fprintf( out, "0 setgray 1 scale stroke ");
  fprintf( out, "grestore } bind def\n");
  
  /* Upper part of bond */
  fprintf(out, "\n/B2 {ucache gsave newpath rotate 2 copy moveto\n 5 1 roll 5 1 roll dup neg 0 4 1 roll arcn\n 2 copy neg lineto gsave 3 1 roll 0 translate\n 1 scale 0 0 3 -1 roll -90 90 arc fill grestore\n ");
  if ( blinew > 0.0 ) 
    fprintf( out, "0 setgray stroke ");
  fprintf( out, "grestore } bind def\n"); 

  /* Light shadow for bonds */
  if ( llinew > 0.0 ) {

    fprintf(out, "\n/BL {ucache gsave newpath rotate 2 copy moveto\n 3 index 3 index lineto neg moveto neg lineto %.6f setlinewidth\n ", llinew);
    printcolor(out,0,dist_min,3);
    fprintf(out, "setrgbcolor stroke grestore} bind def\n");
  } 

  /* Bond separator */
  if ( slinew > 0.0 ) {

    fprintf(out, "\n/BS1 {ucache gsave newpath rotate 0 translate 1 index 1 scale\n 0 0 3 -1 roll -90 90 arc 1 exch div 1 scale\n %.6f setlinewidth 0 setgray stroke grestore } bind def\n", slinew);

    fprintf(out, "\n/BS2 {ucache gsave newpath rotate 1 index moveto neg 0 exch 2 mul rlineto\n %.6f setlinewidth 0 setgray stroke grestore } bind def\n", 
	    slinew);
  }

  fprintf(out, "\n/B {ucache gsave setrgbcolor ");
  if ( blinew > 0.0 )
      fprintf(out, "setlinewidth ");
  fprintf(out, "translate} bind def\n\n");

  fprintf(out, "1 setlinejoin\n");

  /* Set background color */
  fprintf(out, "\n%%%% Backgroundcolor:\ngsave\n");
  printcolor(out, 0, 0, 3);
  fprintf(out, "setrgbcolor\nnewpath %d %d moveto\n", 
	  (int) ( origin.x - offset ), 
	  (int) ( origin.y - scaling * spect_offset - offset ) );
  fprintf(out,"%d %d lineto\n", 
	  (int) ( origin.x +  scaling * width.x + offset ),
	  (int) ( origin.y - scaling * spect_offset - offset ) );
  fprintf(out,"%d %d lineto\n",
	  (int) ( origin.x +  scaling * width.x + offset ),
	  (int) ( origin.y +  scaling * width.y + offset ));
  fprintf(out,"%d %d lineto\n",
	  (int) ( origin.x - offset ),
	  (int) ( origin.y +  scaling * width.y + offset ));
  fprintf(out,"closepath fill grestore\n\n");

  if ( spect >= 0.0 ) {

      fprintf(out, "gsave\n");

      fprintf(out, "newpath %d %d moveto\n", 
	      (int) ( origin.x - offset ), 
	      (int) ( origin.y - offset ) );
      fprintf(out, "%d %d lineto\n", 
	      (int) ( origin.x +  scaling * width.x + offset ),
	      (int) ( origin.y - offset ) );
      fprintf(out, "%d %d lineto\n",
	      (int) ( origin.x +  scaling * width.x + offset ),
	      (int) ( origin.y +  scaling * width.y + offset ));
      fprintf(out, "%d %d lineto\n",
	      (int) ( origin.x - offset ),
	      (int) ( origin.y +  scaling * width.y + offset ));
      fprintf(out, "closepath clip\n\n");
  }

  /* Scaling (changes also bounding box of eps file) */
  if ( scaling != 1.0 ) {
    fprintf(out, "%.6f %.6f scale\n", scaling, scaling); 
    fprintf(out, "%.6f %.6f translate\n", 
	    ( 1.0 / scaling - 1.0 ) * origin.x, 
	    ( 1.0 / scaling - 1.0 ) * origin.y );
  }

  /* Translate picture */
  if ( translation.x != 0.0 || translation.y != 0.0 )
    fprintf(out, "%.6f %.6f translate\n", 
	    translation.x * sclfct, translation.y * sclfct);

  /* Zoom */
  if ( zoom != 1.0 ) {
    fprintf(out, "%.6f %.6f scale\n", zoom, zoom);
    fprintf(out, "%.6f %.6f translate\n", 
	    ( 1.0 / zoom - 1.0 ) * center.x, 
	    ( 1.0 / zoom - 1.0 ) * center.y );
  }

  /* Draw background part of bounding box */
  if ( frame >= 0 && frame_linew != 0.0 ) {
    fprintf(out, "\n%%%% Bounding box:\n");
    fprintf(out, "gsave 1 setlinecap\n");
    setlinew(out, frame_linew);
    printcolor(out, 0, 0, 4);
    fprintf(out, "setrgbcolor\nnewpath %.2f %.2f moveto\n", 
	    cornerps[4].x, cornerps[4].y);
    fprintf(out, "%.2f %.2f lineto\n", 
	    cornerps[5].x, cornerps[5].y);
    fprintf(out, "%.2f %.2f lineto\n", 
	    cornerps[6].x, cornerps[6].y);
    fprintf(out, "%.2f %.2f lineto closepath\n", cornerps[7].x, cornerps[7].y);
    if ( frame_border_linew > 0.0 ) {
	fprintf(out, "gsave ");
	printcolor(out, 3, 0, 4);
	fprintf(out, "setrgbcolor %f setlinewidth stroke grestore\n"
		, frame_linew + frame_border_linew);
    }
    fprintf(out, "stroke\n");    

    boxdist.x = ( cornerps[1].x - cornerps[0].x 
		  - cornerps[5].x + cornerps[4].x ) / 2.0;
    boxdist.y = ( cornerps[3].y - cornerps[0].y 
		  - cornerps[7].y + cornerps[4].y ) / 2.0;

    zentrum.x = ( cornerps[5].x + cornerps[4].x 
		  - cornerps[0].x - cornerps[1].x ) / 2.0; 
    zentrum.y = ( cornerps[7].y + cornerps[4].y 
		  - cornerps[0].y - cornerps[3].y ) / 2.0; 

    if ( zentrum.x <= boxdist.x && zentrum.x >= -boxdist.x &&
	 zentrum.y <= boxdist.y && zentrum.y >= -boxdist.y ) 
      ;
    else if ( zentrum.x > TOL && zentrum.y > TOL ) 
      corner1 = 2;
    else if ( zentrum.x < -TOL && zentrum.y > TOL ) 
      corner1 = 3;
    else if ( zentrum.x < -TOL && zentrum.y < -TOL ) 
      corner1 = 0;
    else if ( zentrum.x > TOL && zentrum.y < -TOL ) 
      corner1 = 1;
    else if ( zentrum.x <= TOL &&  zentrum.x >= -TOL && 
	      zentrum.y > 0.0 ) {
      corner1 = 2;
      corner2 = 3;
    }
    else if ( zentrum.x < 0.0 && zentrum.y <= TOL &&  
	      zentrum.y >= -TOL ) {
      corner1 = 0;
      corner2 = 3;
    }
    else if ( zentrum.x <= TOL &&  zentrum.x >= -TOL && 
	      zentrum.y < 0.0 ) {
      corner1 = 0;
      corner2 = 1;
    }
    else if ( zentrum.x > 0.0 && zentrum.y <= TOL &&  
	      zentrum.y >= -TOL ) {
      corner1 = 1;
      corner2 = 2;
    }

    for ( i=0; i<4; i++ ) 
      if ( corner1 != i && corner2 != i ) { 
	fprintf(out, "newpath %.2f %.2f moveto\n", 
		cornerps[i].x, cornerps[i].y);
	fprintf(out, "%.2f %.2f lineto\n", 
		cornerps[i+4].x, cornerps[i+4].y);
	if ( frame_border_linew > 0.0 ) {
	    fprintf(out, "gsave ");
	    printcolor(out, 3, 0, 4);
	    fprintf(out, "setrgbcolor %f setlinewidth stroke grestore\n"
		    , frame_linew + frame_border_linew);
	}
	fprintf(out, "stroke\n");
      }
    fprintf(out, "grestore\n\n");
  }

  /* Draw background part of simulation box */
  if ( sframe >= 0.0 && sframe_linew != 0.0 ) {
    fprintf(out, "%%%% Simulation box:\n");
    fprintf(out, "gsave 1 setlinecap\n");
    setlinew(out, sframe_linew);
    printcolor(out, 1, 0, 4);
    fprintf(out, "setrgbcolor\n"); 
    for ( i=4; i<8; i++ ) {
      tmpi =  scrnorder[i];
      for ( j=0; j<3; j++ )  
	tmpiv[j] = scorner[tmpi].ngh[j];

      for ( j=0; j<3; j++ ) 
	if ( tmpiv[j] >= 0 ) {
	  fprintf(out, "gsave newpath %.2f %.2f moveto\n", 
		  scornerps[tmpi].x, scornerps[tmpi].y );
	  fprintf(out, "%.2f %.2f lineto\n", 
		  scornerps[tmpiv[j]].x, scornerps[tmpiv[j]].y );
	  if ( frame_border_linew > 0.0 ) {
	      fprintf(out, "gsave ");
	      printcolor(out, 3, 0, 4);
	      fprintf(out, "setrgbcolor %f setlinewidth stroke grestore\n"
		      , sframe_linew + frame_border_linew);
	  }
	  if ( axes == 1 ) {
	      if ( (tmpi == 0 && tmpiv[j] == 1) || (tmpi == 1 && tmpiv[j] == 0) )
		  fprintf(out, "1 0 0 setrgbcolor ");
	      else if ( (tmpi == 0 && tmpiv[j] == 3) || (tmpi == 3 && tmpiv[j] == 0) )
		  fprintf(out, "1 1 0 setrgbcolor ");
	      else if ( (tmpi == 0 && tmpiv[j] == 4) || (tmpi == 4 && tmpiv[j] == 0) )
		  fprintf(out, "0 1 0 setrgbcolor ");
	  }
	  fprintf(out, "stroke grestore\n");
	  for ( k=0; k<3; k++ )
	    if ( scorner[tmpiv[j]].ngh[k] == tmpi )
	      scorner[tmpiv[j]].ngh[k] = -1;
	}
    }
    fprintf(out, "grestore\n\n");
  }

  /* Draw background part of real bounding box */
  if ( rframe >= 0.0 && rframe_linew != 0.0 ) {
    fprintf(out, "%%%% System bounding box:\n");
    fprintf(out, "gsave 1 setlinecap\n");
    setlinew(out, rframe_linew);
    printcolor(out, 2, 0, 4);
    fprintf(out, "setrgbcolor\n");
    for ( i=4; i<8; i++ ) {
      tmpi =  rcrnorder[i];
      for (j=0; j<3; j++ )  
	tmpiv[j] = rcorner[tmpi].ngh[j];

      for ( j=0; j<3; j++ ) 
	if ( tmpiv[j] >= 0 ) {
	  fprintf(out, "gsave newpath %.2f %.2f moveto\n", 
		  rcornerps[tmpi].x, rcornerps[tmpi].y );
	  fprintf(out, "%.2f %.2f lineto\n", 
		  rcornerps[tmpiv[j]].x, rcornerps[tmpiv[j]].y );
	  if ( frame_border_linew > 0.0 ) {
	      fprintf(out, "gsave ");
	      printcolor(out, 3, 0, 4);
	      fprintf(out, "setrgbcolor %f setlinewidth stroke grestore\n", 
		      rframe_linew + frame_border_linew);
	  }
	  if ( axes == 1 ) {
	      if ( (tmpi == 0 && tmpiv[j] == 1) || (tmpi == 1 && tmpiv[j] == 0) )
		  fprintf(out, "1 0 0 setrgbcolor ");
	      else if ( (tmpi == 0 && tmpiv[j] == 3) || (tmpi == 3 && tmpiv[j] == 0) )
		  fprintf(out, "1 1 0 setrgbcolor ");
	      else if ( (tmpi == 0 && tmpiv[j] == 4) || (tmpi == 4 && tmpiv[j] == 0) )
		  fprintf(out, "0 1 0 setrgbcolor ");
	  }
	  fprintf(out, "stroke grestore\n");
	  for ( k=0; k<3; k++ )
	    if ( rcorner[tmpiv[j]].ngh[k] == tmpi )
	      rcorner[tmpiv[j]].ngh[k] = -1;
	}
    }
    fprintf(out, "grestore\n");
  }


  fprintf(out, "\n\n");  

  /* Draw atoms */ 
  for ( listelmt = atomlist; listelmt != NULL; listelmt = listelmt->next ) {

    p = listelmt->cl;
    i = listelmt->num;
    distancep = listelmt->d;

    /* Position of atom (real) */
    pos = p->ort[i];

    distance = fdistance(pos);

    type = p->sorte[i];

    /* Position of atom (ps) */
    ort = proj( pos, sclfct, trans);

    /* Radius of atom (ps) */
    ratom  = scale(pos, sclfct) * r_atom * setradius(type);

    if ( ratom > 0.0 ) {
	
	/* Color of atom */
	if ( colorencoding )
	    printcolor(out, p->enc[i], distance, 0);
	else
	    printcolor(out, type,      distance, 0);

      fprintf(out, "%.3f %.3f %.3f A\n", ratom, ort.x, ort.y);

      ++natoms_drawn;
    }

    /* Draw bonds (belonging to current atom) */

    neigh = p->neightab_array + i;
    
    if ( neigh->n > 0 && r_bond > 0.0 && ratom > 0.0 ) {

      /* Starting point of bond */
      fprintf(out, "%.2f %.2f ", ort.x, ort.y);

      /* Set linewidth of bond borders */
      if ( blinew > 0.0 )
	fprintf(out, "%.2f ", blinew);
      
      /* Set color of bonds */
      if ( bondcolor == 0 && colorencoding )
	  printcolor(out, p->enc[i], distance, 1);
      else
	  printcolor(out, type,      distance, 1);

      fprintf(out, "B\n");

      /* For all neighbour atoms */
      for ( k=0; k<neigh->n; k++ ) {

	neigh_cell = neigh->cl[k];
	neigh_num  = neigh->num[k];
	neigh_type = neigh->typ[k];

	/* Position of neighbour atom (real) */
	npos = neigh_cell->ort[neigh_num];

	/* Distance of neighbour atom from focal point */
	ndist.x = foc.x - npos.x;
	ndist.y = foc.y - npos.y;
	ndist.z = foc.z - npos.z;
	ndistance = sqrt( SPROD(ndist,ndist) );

	ndistancep = pdistance(npos);

	/* Differences of distances of atom and neighbour */
	nreldistance = distancep - ndistancep;

	/* Draw bond only when neighbour atom is visible */
	if ( setradius(neigh_type) > 0.0 ) {
  
	  /* Position of neighbour atom (ps) */
	  nort = proj( npos, sclfct, trans);

	  /* Relative position (real) */
	  deltapos.x = npos.x - pos.x;
	  deltapos.y = npos.y - pos.y;
	  deltapos.z = npos.z - pos.z;
	  
	  /* Relative position (ps) */
	  deltaort.x = nort.x - ort.x;
	  deltaort.y = nort.y - ort.y;
	  
	  /* bondlength (real) */
	  bondlength = sqrt(SPROD(deltapos,deltapos));
	  
	  /* bondlength (ps) */
	  atomdistance = sqrt( deltaort.x * deltaort.x 
			       + deltaort.y * deltaort.y ); 

	  /* Atoms should have different x,y positions , preliminary */
	  if ( atomdistance > 0.0 ) { 
	  
	    /* Unit vector in bond direction (real) */
	    unit.x = deltapos.x / bondlength;
	    unit.y = deltapos.y / bondlength;
	    unit.z = deltapos.z / bondlength;
	  
	    /* Unit vector in bond direction (ps) */
	    unitps.x = deltaort.x / atomdistance;
	    unitps.y = deltaort.y / atomdistance;
	  
	    /* Bond radius (ps) */
	    radscl = ( setradius(type) <= setradius(neigh_type) ) ? 
	      setradius(type) : setradius(neigh_type);
	    rbond = scale(pos, sclfct)  * r_bond * radscl;

	    /* Bond radius at end point of bond (ps) */
	    nrbond = scale( npos, sclfct) * r_bond * radscl;

	    /* Starting point of bond (ps) 
	     in coordinate system of atom */

	    if ( rbond > ratom )
	      rbond = ratom;

	    /* Intersection point of bond with sphere (real) */
	    inters.x = pos.x + unit.x * r_atom * setradius(type);
	    inters.y = pos.y + unit.y * r_atom * setradius(type);
	    inters.z = pos.z + unit.z * r_atom * setradius(type);

	    /* Projection of intersection point (ps) */
	    intersps = proj(inters, sclfct, trans);

	    /* Intersection point relative to atom (ps) */
	    intersps.x -= ort.x;
	    intersps.y -= ort.y;

	    interspsdist = sqrt( intersps.x * intersps.x + intersps.y * intersps.y );

	    /* Intersection point with neighbour atom */
	    if ( needle >= 0.0 && needle <= 1.0 && nreldistance > 0.0 ) {

		ninters.x =  npos.x - unit.x * r_atom * setradius(neigh_type);
		ninters.y =  npos.y - unit.y * r_atom * setradius(neigh_type);
		ninters.z =  npos.z - unit.z * r_atom * setradius(neigh_type);

		nintersps = proj(ninters, sclfct, trans);

		nintersps.x -= ort.x;
		nintersps.y -= ort.y;

		ninterspsdist = sqrt( nintersps.x * nintersps.x 
				      + nintersps.y * nintersps.y );
	    }

	    d_rbond = nrbond - rbond;
	    tmp2 = SQR(atomdistance) + SQR(d_rbond);
	    corr = ratom - atomdistance / tmp2 
		* ( sqrt( SQR(ratom) * tmp2 - SQR(rbond * atomdistance) ) 
		    - rbond * d_rbond );

	    hordist = ratom - corr;

	    brbond = rbond + d_rbond / atomdistance * ( ratom - corr );

	    if ( nreldistance <= 0.0 )
		spsdist = interspsdist * ( 1.0 - corr / ratom );
	    else
		spsdist = hordist;

	    /* If bonds are drawn as needles */
	    if ( needle >= 0.0 && needle <= 1.0 ) {
	      if ( nreldistance < 0.0 )
		brbond *= needle + ( 1.0 - needle ) 
		    * ( hordist - spsdist ) / ( atomdistance - spsdist );
	      else if ( nreldistance > 0.0 ) 
		  brbond *= 1.0 - ( 1.0 - needle ) * hordist / ninterspsdist;	    
	    } 

	    irbond = rbond + d_rbond / atomdistance * spsdist;

	    /* If bonds are drawn as needles */
	    if ( needle >= 0.0 && needle <= 1.0 ) 
	      if ( nreldistance < 0.0 )
		irbond *= needle;

	    /* End point of bond (real) (half the bondlength) */
	    ns.x = pos.x + unit.x * bondlength / 2.0;
	    ns.y = pos.y + unit.y * bondlength / 2.0;
	    ns.z = pos.z + unit.z * bondlength / 2.0;

	    /* Distance of end point of bond from focal point */
	    nsdist.x = foc.x - ns.x;
	    nsdist.y = foc.y - ns.y;
	    nsdist.z = foc.z - ns.z;
	    nsdistance = sqrt( SPROD(nsdist,nsdist) );
	    
	    /* End point of bond (ps) (half the bondlength) */
	    nsps = proj(ns, sclfct, trans);

	    nsps.x -= ort.x;
	    nsps.y -= ort.y;
	    
	    nspsdist = sqrt( nsps.x * nsps.x + nsps.y * nsps.y );

	    /* Make distance 1 percent longer, prov. */
	    nspsdist += 0.01 * nspsdist;

	    /* Bond radius at half of bond (ps) */
	    mrbond = scale( ns, sclfct) * r_bond * radscl;

	    /* If bonds are drawn as needles */
	    if ( needle >= 0.0 && needle <= 1.0 ) {
		if ( nreldistance < 0.0 ) 
		    mrbond *= needle + ( 1.0 - needle ) 
			* ( nspsdist - spsdist ) / ( atomdistance - spsdist );
		else if ( nreldistance > 0.0 ) 
		    mrbond *= 1.0 - ( 1.0 - needle ) * nspsdist / ninterspsdist; 
	    }

	    /* Radius of neighbour atom (ps) */
	    nratom = scale( npos, sclfct) * r_atom * setradius( neigh_type );
    
	    /* Normal vector to bond (ps) */
	    normalps.x = - unitps.y;
	    normalps.y =   unitps.x;

	    mrbonddist2 = ( nsps.x + normalps.x * mrbond ) 
		* ( nsps.x + normalps.x * mrbond ) 
		+ ( nsps.y + normalps.y * mrbond ) 
		* ( nsps.y + normalps.y * mrbond ); 

	    /* Neighbour has larger z coordinate */
	    if ( nreldistance < 0.0 ) {
	      
	      /* Axes of ellipse */
	      ellps = SPROD(unit, nsdist) / nsdistance;
	      if ( ellps < 0.0 ) 
		  ellps = -ellps;

	      /* Angle of bond synapse */
	      if ( unitps.y > 0 )
		angle = acos( unitps.x ) / PIN;
	      else 
		angle = acos( - unitps.x ) / PIN + 180;
	      
	      /* Draw bond */
	      fprintf(out, "%.2f %.2f %.2f %.2f %.2f %.2f %.2f B1\n", 
		      1.0/ellps, irbond, ellps, spsdist, nspsdist, mrbond, angle);

	      ++nbonds_drawn;

	      /* Draw lightshadow */
	      if ( llinew > 0.0 ) 
		fprintf(out, "%.2f %.2f %.2f %.2f %.2f BL\n",
			nspsdist, mrbond + ( blinew + llinew ) / 2.0, 
			hordist + alinew / 2.0, 
			brbond + ( blinew + llinew ) / 2.0, angle);
	    }
	    
	    /* Neighbour has lower z coordinate */
	    else if ( nreldistance >= 0.0 ) {
	      
	      /* Axis of ellipse */
	      ellps = - SPROD(unit, nsdist) / nsdistance;
	      if ( ellps < 0.0 ) 
		  ellps = -ellps;
	      
	      if ( unitps.y > 0 )
		angle = acos( unitps.x ) / PIN;
	      else 
		angle = acos( - unitps.x ) / PIN + 180;
      
	      corr = interspsdist - hordist;

	      if (corr < 0.0 )
		corr = 0.0;
	      
	      if ( corr == 0.0 )
		tmp = ratom;
	      else
		tmp = sqrt( SQR(hordist + corr * 50.0) + SQR(brbond) );
	      
	      deltaangle = asin( brbond / tmp ) / PIN;
	      
	      urspr.x = - unitps.x * corr * 50.0;
	      urspr.y = - unitps.y * corr * 50.0;
 
      
	      /* Draw only if bond separator is outside atom */
	      if( mrbonddist2 > ratom * ratom ) {
		
		/* Draw bond */
 		fprintf(out, "%.2f %.2f %.2f %.2f %.2f %.2f %.2f B2\n", 
			ellps, -corr*50.0, tmp, deltaangle, 
			nspsdist, mrbond, angle);
	
		/* Draw lightshadow */
		if ( llinew > 0.0 ) 
		  fprintf(out, "%.2f %.2f %.2f %.2f %.2f BL\n",
			  nspsdist, mrbond + ( blinew + llinew ) / 2.0, 
			  hordist + alinew / 2.0, 
			  brbond + ( blinew + llinew ) / 2.0, angle);

		/* Draw bond separator */
		if ( slinew > 0.0 ) {
		  if ( ellps > 0.01 ) 
		    fprintf(out, "%.2f %.2f %.2f %.2f BS1\n", 
			    ellps, mrbond, nspsdist, angle);
		  else 
		    fprintf(out, "%.2f %.2f %.2f BS2\n", 
			    mrbond, nspsdist, angle); 
		}
		
	      }
	      
	    }
	    
	  }
	  
	}

      } /* For all neighbour atoms */

      fprintf(out,"Y\n");       
    }    
  }

  /* Draw foreground part of bounding box */
  if ( frame >= 0.0 && frame_linew != 0.0 ) {

    fprintf(out, "\n%%%% Bounding box:\n");
    fprintf(out, "gsave 1 setlinecap\n");
    setlinew(out, 2*frame_linew);
    printcolor(out, 0, dist_min, 4);
    
    fprintf(out, "setrgbcolor\nnewpath %.2f %.2f moveto\n", 
	    cornerps[0].x, cornerps[0].y);
    fprintf(out, "%.2f %.2f lineto\n", 
	    cornerps[1].x, cornerps[1].y);
    fprintf(out, "%.2f %.2f lineto\n", 
	    cornerps[2].x, cornerps[2].y);
    fprintf(out, "%.2f %.2f lineto closepath\n", 
	    cornerps[3].x, cornerps[3].y);
    if ( frame_border_linew > 0.0 ) {
	fprintf(out, "gsave ");
	printcolor(out, 3, dist_min, 4);	    
	fprintf(out, "setrgbcolor %f setlinewidth stroke grestore\n"
		, 2*frame_linew + frame_border_linew);
    }
    fprintf(out, "stroke\n");

    for ( i=0; i<4; i++ ) 
      if ( corner1 == i || corner2 == i ) { 
	fprintf(out, "newpath %.2f %.2f moveto\n", 
		cornerps[i].x, cornerps[i].y);
	fprintf(out, "%.2f %.2f lineto\n", 
		cornerps[i+4].x, cornerps[i+4].y);
	if ( frame_border_linew > 0.0 ) {
	    fprintf(out, "gsave ");
	    printcolor(out, 3, dist_min, 4);	    
	    fprintf(out, "setrgbcolor %f setlinewidth stroke grestore\n"
		    , 2*frame_linew + frame_border_linew);
	}
	fprintf(out, "stroke\n");

      }
    fprintf(out, "grestore\n");
  }

  /* Foreground part of simulation box */
  if ( sframe >= 0.0 && sframe_linew != 0.0 ) {

    fprintf(out, "\n\n%%%% Simulation box:\n");
    fprintf(out, "gsave 1 setlinecap\n");
    setlinew(out, 1.5*sframe_linew);
    printcolor(out, 1, dist_min, 4);
    fprintf(out, "setrgbcolor\n");

    for ( i=0; i<4; i++ ) {
      tmpi = scrnorder[i];
      for ( j=0; j<3; j++ )
	tmpiv[j] = scorner[tmpi].ngh[j];

      for ( j=0; j<4; j++ ) 
	for ( k=0; k<3; k++ )
	  if ( tmpiv[k] == scrnorder[j] ) {
	    fprintf(out, "gsave newpath %.2f %.2f moveto\n", 
		    scornerps[tmpi].x, scornerps[tmpi].y );
	    fprintf(out, "%.2f %.2f lineto\n", 
		    scornerps[tmpiv[k]].x, scornerps[tmpiv[k]].y );
	    if ( frame_border_linew > 0.0 ) {
		fprintf(out, "gsave ");
		printcolor(out, 3, dist_min, 4);	
		fprintf(out, "setrgbcolor %f setlinewidth stroke grestore\n"
			, 1.5*sframe_linew + frame_border_linew);
	    }
	    if ( axes == 1 ) {
		if ( (tmpi == 0 && tmpiv[k] == 1) || (tmpi == 1 && tmpiv[k] == 0) )
		    fprintf(out, "1 0 0 setrgbcolor ");
		if ( (tmpi == 0 && tmpiv[k] == 3) || (tmpi == 3 && tmpiv[k] == 0) )
		    fprintf(out, "1 1 0 setrgbcolor ");
		if ( (tmpi == 0 && tmpiv[k] == 4) || (tmpi == 4 && tmpiv[k] == 0) )
		    fprintf(out, "0 1 0 setrgbcolor ");
	    }
	    fprintf(out, "stroke grestore\n");
	    for ( l=0; l<3; l++ )
	      if ( scorner[tmpiv[k]].ngh[l] == tmpi )
		scorner[tmpiv[k]].ngh[l] = -1;
	  }
    }
    fprintf(out, "grestore\n");
  }

  /* Foreground part of real bounding box */
  if ( rframe >= 0.0 && rframe_linew != 0.0 ) {

    fprintf(out, "\n\n%%%% System bounding box:\n");
    fprintf(out, "gsave 1 setlinecap\n");
    setlinew(out, 1.5*rframe_linew);
    printcolor(out, 2, dist_min, 4);
    fprintf(out, "setrgbcolor\n");
    for ( i=0; i<4; i++ ) {
      tmpi = rcrnorder[i];
      for ( j=0; j<3; j++ )
	tmpiv[j] = rcorner[tmpi].ngh[j];

      for ( j=0; j<4; j++ ) 
	for ( k=0; k<3; k++ )
	  if ( tmpiv[k] == rcrnorder[j] ) {
	    fprintf(out, "gsave newpath %.2f %.2f moveto\n", 
		    rcornerps[tmpi].x, rcornerps[tmpi].y );
	    fprintf(out, "%.2f %.2f lineto\n"
		    , rcornerps[tmpiv[k]].x, rcornerps[tmpiv[k]].y );
	    if ( frame_border_linew > 0.0 ) {
		fprintf(out, "gsave ");
		printcolor(out, 3, dist_min, 4);
		fprintf(out, "setrgbcolor %f setlinewidth stroke grestore\n"
			, 1.5*rframe_linew + frame_border_linew);
	    }
	    if ( axes == 1 ) {
		if ( (tmpi == 0 && tmpiv[k] == 1) || (tmpi == 1 && tmpiv[k] == 0) )
		    fprintf(out, "1 0 0 setrgbcolor ");
		if ( (tmpi == 0 && tmpiv[k] == 3) || (tmpi == 3 && tmpiv[k] == 0) )
		    fprintf(out, "1 1 0 setrgbcolor ");
		if ( (tmpi == 0 && tmpiv[k] == 4) || (tmpi == 4 && tmpiv[k] == 0) )
		    fprintf(out, "0 1 0 setrgbcolor ");
	    }
	    fprintf(out, "stroke grestore\n");
	    for ( l=0; l<3; l++ )
	      if ( rcorner[tmpiv[k]].ngh[l] == tmpi )
		rcorner[tmpiv[k]].ngh[l] = -1;
	  }
    }
    fprintf(out, "grestore\n");
  }

  fprintf(out, "\n\n");
 
  if ( spect >= 0.0 ) {
  fprintf(out, "grestore\n");
    /* Draw line */
    if (spect > 0.0 ) {
      fprintf(out, "%.2f setlinewidth\n", spect);
      fprintf(out, "newpath %d %d moveto\n", 
	      (int) ( origin.x - offset ), 
	      (int) ( origin.y - offset ) );
      fprintf(out, " %d %d lineto stroke\n", 
	      (int) ( origin.x +  scaling * width.x + offset ),
	      (int) ( origin.y - offset ) ); 
    }  
    tmpv.x = origin.x + scaling * width.x / 4.0;
    tmpv.y = origin.y - offset - scaling * spect_height * 0.75;
    size.x = scaling * width.x / 2.0;
    size.y = scaling * spect_height / 2.0;
    
    if ( eng != 0 || kineng != 0 || color_enc != 0 )
      draw_colorencoding(out, tmpv, size);

  }

  fprintf(out, "end grestore\n%%%%EOF\n");
  fclose(out);

}

#else

/******************************************************************************
*
*  draw_pictures -- 2d version
*
******************************************************************************/

void draw_picture2d(void) {

  int l, k, j, i, type, neigh_num, neigh_type;
  cell *p, *neigh_cell;
  vektor pos, npos;
  neightab *neigh;
  FILE *out;
  str255 fname;
  int llx=30, lly=40, urx=570, ury=750;
  real ratio, boxratio;
  real ratom;
  vektor origin, center;
  vektor2d ort;
  real sclfct, offset = 20;
  vektor2d width, cornerps[4], scornerps[4], rcornerps[4], tmpv;
  vektor maxim={Max,Max}, minim={Min,Min};
  vektor2d pnt;
  vektor smaxim={Max,Max}, sminim={Min,Min}; 
  vektor rmaxim={Max,Max}, rminim={Min,Min};
  vektor maximage={Max,Max}, minimage={Min,Min};
  real spect_height = 72.0, spect_offset = 0.0;
  vektor2d nort, deltaort, unitps, normalps, trans;
  real atomdistance, rbond, nratom;
  real angle;
  real radscl, atom_vol;
  real pts;
  real asize, adelta, ax, atmp, arad, ashade;
  int N;
  real spsdist, nspsdist;
  vektor2d size;
  time_t now;
  char *revision = "$Revision$";
  char version[10]  = "        ";

  /* Estimate atom and bond radii */
  if ( wireframe < 0.0 ) {
    if ( r_atom == -1.0 ) {
      atom_vol = ( maxl.x - minl.x ) * ( maxl.y - minl.y ) / natoms;
      r_atom = 2*sqrt ( atom_vol / 100.0 );
      printf("Largest atom radius (R): %.2f (estimated)\n", r_atom);
    }
    if ( r_bond == -1.0 ) {
      r_bond = r_atom / 6.0;
      printf("Largest bond radius (B): %.2f (estimated)\n", r_bond);    
    }
  }

  /* Adjust bond radii, r_atom must be greater than r_bond */
  if ( r_bond > r_atom )
    r_atom = r_bond;

  /* Take atom radii into account */
  maxl.x += r_atom;
  minl.x -= r_atom;
  maxl.y += r_atom;
  minl.y -= r_atom;

  real_maxl.x += r_atom;
  real_minl.x -= r_atom;
  real_maxl.y += r_atom;
  real_minl.y -= r_atom;  

  /* Corners of bounding box in plane */
  tmpv.x = ( maxl.x - minl.x ) / 2.0;
  tmpv.y = ( maxl.y - minl.y ) / 2.0;
 
  pnt.x = -tmpv.x;
  pnt.y = -tmpv.y;
  cornerps[0] = pnt; 

  pnt.x = tmpv.x;
  cornerps[1] = pnt;
 
  pnt.y = tmpv.y;
  cornerps[2] = pnt;

  pnt.x = -tmpv.x;
  cornerps[3] = pnt;

  for ( i=0; i<4; i++ ) {
      maxim.x = MAX( maxim.x, cornerps[i].x );
      maxim.y = MAX( maxim.y, cornerps[i].y );
      minim.x = MIN( minim.x, cornerps[i].x );
      minim.y = MIN( minim.y, cornerps[i].y );
  }  

  /* Corners of simulation box */
  if ( sframe >= 0 && sframe_linew != 0.0 ) {

    pnt.x = - ( maxl.x + minl.x ) / 2.0;
    pnt.y = - ( maxl.y + minl.y ) / 2.0;
    scornerps[0] = pnt;

    pnt.x += box_x.x;
    pnt.y += box_x.y;
    scornerps[1] = pnt;

    pnt.x += box_y.x;
    pnt.y += box_y.y;
    scornerps[2] = pnt;

    pnt.x -= box_x.x;
    pnt.y -= box_x.y;
    scornerps[3] = pnt;

    for ( i=0; i<4; i++ ) {
	smaxim.x = MAX( smaxim.x, scornerps[i].x );
	smaxim.y = MAX( smaxim.y, scornerps[i].y );
	sminim.x = MIN( sminim.x, scornerps[i].x );
	sminim.y = MIN( sminim.y, scornerps[i].y );
    }
  }

  /* Corners of real bounding box */
  if ( rframe >= 0 && rframe_linew != 0.0 ) {

    pnt.x = ursprung.x
      - ( unitv[0].x + unitv[1].x ) * r_atom;
    pnt.y = ursprung.y
      - ( unitv[0].y + unitv[1].y ) * r_atom;

    rcornerps[0] = pnt;

    pnt.x += unitv[0].x * ( real_maxl.x - real_minl.x );
    pnt.y += unitv[0].y * ( real_maxl.x - real_minl.x );
    rcornerps[1] = pnt;

    pnt.x += unitv[1].x * ( real_maxl.y - real_minl.y );
    pnt.y += unitv[1].y * ( real_maxl.y - real_minl.y );
    rcornerps[2] = pnt;

    pnt.x -= unitv[0].x * ( real_maxl.x - real_minl.x );
    pnt.y -= unitv[0].y * ( real_maxl.x - real_minl.x );
    rcornerps[3] = pnt;

    for ( i=0; i<4; i++ ) {
      rmaxim.x = MAX( rmaxim.x, rcornerps[i].x );
      rmaxim.y = MAX( rmaxim.y, rcornerps[i].y );
      rminim.x = MIN( rminim.x, rcornerps[i].x );
      rminim.y = MIN( rminim.y, rcornerps[i].y );
    }
  }

  /* Determine boundaries of image */
  if ( ( sframe < 0.0 || sframe_linew == 0.0 ) &&
       ( rframe < 0.0 || rframe_linew == 0.0 ) ) {
    maximage = maxim;
    minimage = minim;
  }
  else if ( frame < 0.0 || frame_linew == 0) {
    maximage.x = MAX( smaxim.x, rmaxim.x );
    maximage.y = MAX( smaxim.y, rmaxim.y );
    minimage.x = MIN( sminim.x, rminim.x );
    minimage.y = MIN( sminim.y, rminim.y );
  }
  else {
    maximage.x = MAX( smaxim.x, rmaxim.x );
    maximage.y = MAX( smaxim.y, rmaxim.y );
    minimage.x = MIN( sminim.x, rminim.x );
    minimage.y = MIN( sminim.y, rminim.y );

    maximage.x = MAX( maximage.x, maxim.x );
    maximage.y = MAX( maximage.y, maxim.y );
    minimage.x = MIN( minimage.x, minim.x );
    minimage.y = MIN( minimage.y, minim.y );
  }

  if ( spect >= 0.0 ) {
    lly += spect_height;
    spect_offset = spect_height;
  }

  if ( maximage.x == minimage.x ) {
    maximage.x += 0.5;
    minimage.x -= 0.5;
    }
  if ( maximage.y == minimage.y ) {
    maximage.y += 0.5;
    minimage.y -= 0.5;
    }

  /* Ratio of page (ps)*/
  ratio = (real) ( ury - lly ) / ( urx - llx );
  /* Ratio of image */
  boxratio = ( maximage.y - minimage.y ) / ( maximage.x - minimage.x );

  /* Scale factor image -> ps */
  if ( boxratio <= ratio ) 
    sclfct = ( urx - llx ) / ( maximage.x - minimage.x ); 
  else
    sclfct = ( ury - lly ) / ( maximage.y - minimage.y ); 

  /* Width of box (ps) */
  width.x = sclfct * ( maximage.x - minimage.x );
  width.y = sclfct * ( maximage.y - minimage.y );

  /* Diplacement with respect to projected values */
  trans.x = ( urx + llx - sclfct * (  maximage.x + minimage.x ) ) / 2.0;
  trans.y = ( ury + lly - sclfct * (  maximage.y + minimage.y ) ) / 2.0;

  /* Transform bounding box (ps) */
  if ( frame >= 0.0 )
    for ( i=0; i<4; i++ )
      cornerps[i] = transform( cornerps[i], sclfct, trans);

  if ( sframe >= 0.0 )
    /* Transform simulation box (ps) */
    for ( i=0; i<4; i++ )
      scornerps[i] = transform( scornerps[i], sclfct, trans);

  if ( rframe >= 0.0 )
    /* Transform real bounding box (ps) */
    for ( i=0; i<4; i++ )
      rcornerps[i] = transform( rcornerps[i], sclfct, trans);


  /* Origin of box (ps) */
  origin.x = ( urx + llx - width.x ) / 2.0;
  origin.y = ( ury + lly - width.y ) / 2.0;

  /* Center of box (ps) */
  center.x = origin.x + width.x / 2.0;
  center.y = origin.y + width.y / 2.0;


  /* Open output file */
  sprintf(fname,"%s.eps",infilename);
  out = fopen(fname,"w");
  if (NULL == out) 
      error("Cannot open eps file.\n");

  /* Crop offset */
  if ( crop != -1.0 )
      offset = crop + (int) frame_linew;

  /* Box width and height prescribed in cm */
  if ( ps_width != -1.0 || ps_height != -1.0 ) {

    if ( ps_width == -1.0 ) {
      pts = ps_height / PT - 2 * offset;
      scaling = pts / width.y;
    }
    else if ( ps_width > 0.0 ) {
      pts = ps_width / PT - 2.0 * offset;
      scaling = pts / width.x;
    }
  }
  /* Get the creation time */
  time(&now);

  /* Get the version number of imd_ps */
  j = 0;
  for ( i=1; revision[i]!='$'; i++)
      if ( isdigit( revision[i] ) || revision[i] == '.' )
	  version[j++] = revision[i];

  /* Print header*/
  fprintf(out, "%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(out, "%%%%Creator: imd_ps 2D Version %s\n", version);
  fprintf(out, "%%%%CreationDate: %s", ctime(&now));
  fprintf(out, "%%%%For: %s@%s\n", getenv("USER"), getenv("HOST")) ;
  fprintf(out, "%%%%Title: %s.eps\n", infilename);
  fprintf(out, "%%%%BoundingBox: %d %d %d %d\n",
          (int) ( origin.x - offset ),
          (int) ( origin.y - scaling * spect_offset - offset ),
          (int) ( origin.x + scaling * width.x + offset ),
          (int) ( origin.y + scaling * width.y + offset ) );
  fprintf(out, "%%%%EndComments\n\n");

  fprintf(out, "gsave\n10 dict begin\n");

  fprintf(out, "\n/X {gsave} bind def\n");
  fprintf(out, "\n/Y {grestore} bind def\n");

  /* Atom */
  fprintf(out, "\n/Shade {dup 5 index mul 1 index 5 index mul\n 2 index 5 index mul setrgbcolor} def\n");

  fprintf(out, "\n/Draw {mul 0 360 arc closepath fill} def\n");

  fprintf(out, "\n/A {ucache gsave translate gsave\n newpath 0 0 2 index 0 360 arc\n ");

  if ( shading < 0.0 )
    fprintf(out, "3 index 3 index 3 index setrgbcolor closepath fill\n ");
  else {
    fprintf(out, "clip gsave %.2f rotate %.6f 1 scale\n "
	    , shading, sqrt(2.0/3.0));

    asize =  r_atom * sclfct;
    N = 2 * (int) ( asize / 2.0 * zoom * scaling );

    if ( N<3 )
      N = 3;

    adelta = 1.0 / N;

    for ( i=0; i<N; i++ ) {

      ax = i * adelta;

      atmp = 1.0 - 1.07 * ax * ax;
      if ( atmp < 0.0 )
        break;

      arad = sqrt( atmp );

      if ( ambient >= 0.0 )
	  ashade = 0.1 + pow( ax, ambient) * 0.9;
      else
	  ashade = 0.1 + ( 1.0 - pow( ABS(ax - 1.0), -ambient) ) * 0.9;

      fprintf(out, "%.6f Shade pop dup %.6f mul 0  2 index %.6f Draw\n ",
              ashade, ax, arad);
    }

    fprintf(out, "grestore ");

  }
  /* Draw atom border */
  if ( alinew > 0.0 )
    fprintf(out, "grestore newpath %.6f setlinewidth 0 setgray\n 0 0 2 index 0 360 arc closepath stroke\n "
	    , alinew);
  else
    fprintf(out, "grestore\n ");

  /* Draw lightshadow */
  if ( llinew > 0.0 ) {
    printcolor(out, 0, dist_min, 3);
    fprintf(out, "setrgbcolor newpath 0 0 2 index %.6f add 0 360 arc\n %.6f setlinewidth stroke\n "
	    , (alinew + llinew ) / 2.0, llinew);
  }
  fprintf(out, "pop pop pop pop grestore } bind def\n");

  fprintf(out, "\n/B {ucache gsave setrgbcolor");
  if ( blinew > 0.0 )
      fprintf(out, " setlinewidth");
  fprintf(out, " translate} bind def\n");

  fprintf(out, "\n/B1 {ucache gsave newpath rotate 2 copy moveto\n 2 index 1 index lineto 0 1 index neg 2 mul rlineto neg lineto\n closepath gsave fill grestore ");
  if ( blinew > 0.0 )
    fprintf( out, "0 setgray stroke ");
  fprintf( out, "grestore } bind def\n");

  /* Bond separator */
  if ( slinew > 0.0 ) 
    fprintf(out, "\n/BS {ucache gsave newpath rotate 1 index moveto neg 0 exch 2 mul rlineto\n %.6f setlinewidth 0 setgray stroke grestore } bind def\n", 
	    slinew);

  fprintf(out, "1 setlinejoin\n");

  /* Set background color */
  fprintf(out, "\n%%%% Backgroundcolor:\ngsave\n");
  printcolor(out, 0, 0, 3);
  fprintf(out, "setrgbcolor\nnewpath %d %d moveto\n",
          (int) ( origin.x - offset ),
          (int) ( origin.y - scaling * spect_offset - offset ) );
  fprintf(out,"%d %d lineto\n",
          (int) ( origin.x +  scaling * width.x + offset ),
          (int) ( origin.y - scaling * spect_offset - offset ) );
  fprintf(out,"%d %d lineto\n",
          (int) ( origin.x +  scaling * width.x + offset ),
          (int) ( origin.y +  scaling * width.y + offset ));
  fprintf(out,"%d %d lineto\n",
          (int) ( origin.x - offset ),
          (int) ( origin.y +  scaling * width.y + offset ));
  fprintf(out,"closepath fill grestore\n\n");

  if ( spect >= 0.0 ) {

      fprintf(out, "gsave\n");

      fprintf(out, "newpath %d %d moveto\n",
              (int) ( origin.x - offset ),
              (int) ( origin.y - offset ) );
      fprintf(out, "%d %d lineto\n",
              (int) ( origin.x +  scaling * width.x + offset ),
              (int) ( origin.y - offset ) );
      fprintf(out, "%d %d lineto\n",
              (int) ( origin.x +  scaling * width.x + offset ),
              (int) ( origin.y +  scaling * width.y + offset ));
      fprintf(out, "%d %d lineto\n",
              (int) ( origin.x - offset ),
              (int) ( origin.y +  scaling * width.y + offset ));
      fprintf(out, "closepath clip\n\n");
  }

  /* Scaling (changes also bounding box of eps file) */
  if ( scaling != 1.0 ) {
    fprintf(out, "%.6f %.6f scale\n", scaling, scaling);
    fprintf(out, "%.6f %.6f translate\n",
            ( 1.0 / scaling - 1.0 ) * origin.x,
            ( 1.0 / scaling - 1.0 ) * origin.y );
  }

  /* Translate picture */
  if ( translation.x != 0.0 || translation.y != 0.0 )
    fprintf(out, "%.6f %.6f translate\n",
            translation.x * sclfct, translation.y * sclfct);

  /* Zoom */
  if ( zoom != 1.0 ) {
    fprintf(out, "%.6f %.6f scale\n", zoom, zoom);
    fprintf(out, "%.6f %.6f translate\n",
            ( 1.0 / zoom - 1.0 ) * center.x,
            ( 1.0 / zoom - 1.0 ) * center.y );
  }

  /* Draw bounding box */
  if ( frame >= 0 && frame_linew != 0.0 ) {
    fprintf(out, "\n%%%% Bounding box:\n");
    fprintf(out, "gsave\n");
    setlinew(out, frame_linew);
    printcolor(out, 0, 0, 4);
    fprintf(out, "setrgbcolor\nnewpath %.2f %.2f moveto\n",
            cornerps[0].x, cornerps[0].y);
    fprintf(out, "%.2f %.2f lineto\n",
            cornerps[1].x, cornerps[1].y);
    fprintf(out, "%.2f %.2f lineto\n",
            cornerps[2].x, cornerps[2].y);
    fprintf(out, "%.2f %.2f lineto closepath\n",
            cornerps[3].x, cornerps[3].y);
    if ( frame_border_linew > 0.0 ) {
	fprintf(out, "gsave ");
	printcolor(out, 3, 0, 4);   
	fprintf(out, "setrgbcolor %f setlinewidth stroke grestore\n", 
	    frame_linew + frame_border_linew);
    }
    fprintf(out, "stroke\n");    

    fprintf(out, "grestore\n\n");
  }

  /* Draw simulation box */
  if ( sframe >= 0.0 && sframe_linew != 0.0 ) {
    fprintf(out, "\n%%%% Simulation box:\n");
    fprintf(out, "gsave\n");
    setlinew(out, sframe_linew);
    printcolor(out, 1, 0, 4);
    fprintf(out, "setrgbcolor\n");
    fprintf(out, "newpath %.2f %.2f moveto\n",scornerps[0].x, scornerps[0].y);
    for ( i=1; i<4; i++ ) 
	fprintf(out, "%.2f %.2f lineto\n",
            scornerps[i].x, scornerps[i].y);

    fprintf(out, "closepath\n");
    if ( frame_border_linew > 0.0 ) {
	fprintf(out, "gsave ");
	printcolor(out, 3, dist_min, 4);
	fprintf(out, "setrgbcolor %f setlinewidth stroke grestore\n", 
		sframe_linew + frame_border_linew);
    }
    fprintf(out, "stroke\n");
    if ( axes == 1 ) {
	fprintf(out, "newpath %.2f %.2f moveto\n",scornerps[0].x, scornerps[0].y);
	fprintf(out, "%.2f %.2f lineto\n", scornerps[1].x, scornerps[1].y);
        fprintf(out, "1 0 0 setrgbcolor stroke\n");
	fprintf(out, "newpath %.2f %.2f moveto\n",scornerps[0].x, scornerps[0].y);
	fprintf(out, "%.2f %.2f lineto\n", scornerps[3].x, scornerps[3].y);
        fprintf(out, "1 1 0 setrgbcolor stroke\n");
    }
    fprintf(out, "grestore\n\n");
  }

  /* Draw real bounding box */
  if ( rframe >= 0.0 && rframe_linew != 0.0 ) {
    fprintf(out, "%%%% System bounding box:\n");
    fprintf(out, "gsave\n");
    setlinew(out, rframe_linew);
    printcolor(out, 2, 0, 4);
    fprintf(out, "setrgbcolor\n");
    fprintf(out, "newpath %.2f %.2f moveto\n",rcornerps[0].x, rcornerps[0].y);
    for ( i=1; i<4; i++ ) 
	fprintf(out, "%.2f %.2f lineto\n",
            rcornerps[i].x, rcornerps[i].y);

    fprintf(out, "closepath\n");
    if ( frame_border_linew > 0.0 ) {
	fprintf(out, "gsave ");
	printcolor(out, 3, 0, 4);
	fprintf(out, "setrgbcolor %f setlinewidth stroke grestore\n", 
		rframe_linew + frame_border_linew);
    }
    fprintf(out, "stroke\n");
    if ( axes == 1 ) {
	fprintf(out, "newpath %.2f %.2f moveto\n",rcornerps[0].x, rcornerps[0].y);
	fprintf(out, "%.2f %.2f lineto\n", rcornerps[1].x, rcornerps[1].y);
        fprintf(out, "1 0 0 setrgbcolor stroke\n");
	fprintf(out, "newpath %.2f %.2f moveto\n",rcornerps[0].x, rcornerps[0].y);
	fprintf(out, "%.2f %.2f lineto\n", rcornerps[3].x, rcornerps[3].y);
        fprintf(out, "1 1 0 setrgbcolor stroke\n");
    }
    fprintf(out, "grestore\n\n");
  }

  fprintf(out, "\n\n");

  /* Draw atoms */
  for ( i=0; i<cell_dim.x; ++i)
      for ( j=0; j<cell_dim.y; ++j) {
      
	  p = PTR_2D_V(cell_array,i,j,cell_dim);
      
	  for ( l=0; l<p->n; ++l ) {

	      /* Position of atom (real) */
	      pos = p->ort[l];

	      type = p->sorte[l];
    
	      /* Position of atom (ps) */
	      ort = scl(pos, sclfct, trans);

	      ratom  = sclfct * r_atom * setradius(type);

	      if ( ratom > 0.0 ) {

		  /* Color of atom */
		  if ( colorencoding )
		      printcolor(out, p->enc[i], 0, 0);
		  else 
		      printcolor(out, type, 0, 0);
		  
		  fprintf(out, "%.3f %.3f %.3f A\n", ratom, ort.x, ort.y);

		  ++ natoms_drawn;
	      }

	      /* Draw bonds (belonging to current atom) */
	      neigh = p->neightab_array + l;
	      
	      if ( neigh->n > 0 && r_bond > 0.0 && ratom > 0.0 ) {
		  
		  /* Starting point of bond */
		  fprintf(out, "%.2f %.2f ", ort.x, ort.y);
		  
		  /* Set linewidth of bond borders */
		  if ( blinew > 0.0 )
		      fprintf(out, "%.2f ", blinew);
		  
		  /* Set color of bonds */
		  if ( bondcolor == 0 && colorencoding )
		      printcolor(out, p->enc[l], 0, 1);
		  else
		      printcolor(out, type, 0, 1);
		  
		  fprintf(out, "B\n");

		  /* For all neighbour atoms */
		  for ( k=0; k<neigh->n; k++ ) {
		      
		      neigh_cell = neigh->cl[k];
		      neigh_num  = neigh->num[k];
		      neigh_type = neigh->typ[k];
		      
		      /* Position of neighbour atom (real) */
		      npos = neigh_cell->ort[neigh_num];

		      /* Draw bond only when neighbour atom is visible */
		      if ( setradius(neigh_type) > 0.0 ) {

			  /* Position of neighbour atom (ps) */
			  nort = scl(npos, sclfct, trans);

			  /* Relative position (ps) */
			  deltaort.x = nort.x - ort.x;
			  deltaort.y = nort.y - ort.y;
		  
			  /* bondlength (ps) */
			  atomdistance = sqrt( deltaort.x * deltaort.x
					       + deltaort.y * deltaort.y );

			  /* Unit vector in bond direction (ps) */
			  unitps.x = deltaort.x / atomdistance;
			  unitps.y = deltaort.y / atomdistance;

			  /* Bond radius (ps) */
			  radscl = setradius(type) <= setradius(neigh_type) ? 
			      setradius(type) : setradius(neigh_type);
			  rbond = sclfct * r_bond * radscl;

			  /* Starting point of bond (ps) */
			  if ( rbond > ratom )
			      rbond = ratom;

			  spsdist = sqrt( ratom * ratom - rbond * rbond );
	  
			  nspsdist = atomdistance / 2.0;
		  
			  /* Radius of neighbour atom (ps) */
			  nratom = sclfct * r_atom * setradius(neigh_type);

			  /* Normal vector to bond (ps) */
			  normalps.x = - unitps.y;
			  normalps.y =   unitps.x;

			  /* Angle of bond synapse */
			  if ( unitps.y > 0 )
			      angle = acos(   unitps.x ) / PIN;
			  else
			      angle = acos( - unitps.x ) / PIN + 180;
			  
			  /* Draw bond */
			  fprintf(out, "%.2f %.2f %.2f %.2f B1\n",
				  nspsdist, spsdist, rbond, angle);


			  /* Draw bond separator */
			  if ( slinew > 0.0 ) 
			      fprintf(out, "%.2f %.2f %.2f BS\n",
				      rbond, nspsdist, angle);

			  ++nbonds_drawn;
		      }
		      
		  }

		  fprintf(out, "Y\n");
	      }
	  }

	  if ( spect >= 0.0 ) {
	      fprintf(out, "grestore\n");
	      /* Draw line */
	      if (spect > 0.0 ) {
		  fprintf(out, "%.2f setlinewidth\n", spect);
		  fprintf(out, "newpath %d %d moveto\n",
			  (int) ( origin.x - offset ),
			  (int) ( origin.y - offset ) );
		  fprintf(out, " %d %d lineto stroke\n",
			  (int) ( origin.x +  scaling * width.x + offset ),
			  (int) ( origin.y - offset ) );
	      }
	      tmpv.x = origin.x + scaling * width.x / 4.0;
	      tmpv.y = origin.y - offset - scaling * spect_height * 0.75;
	      size.x = scaling * width.x / 2.0;
	      size.y = scaling * spect_height / 2.0;
	      
	      if ( eng != 0 || kineng != 0 || color_enc != 0 )
		  draw_colorencoding(out, tmpv, size);	      
	  }
  }

  fprintf(out, "end grestore\n%%%%EOF\n");
  fclose(out);

}

#endif



