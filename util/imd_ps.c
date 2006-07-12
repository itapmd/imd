
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2006 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_ps -- create a postscript picture of an IMD configuration
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#ifndef PS
#define PS
#endif

#define MAIN

#include "util.h"

/******************************************************************************
*
*  Usage -- educate users
*
******************************************************************************/

void usage(void)
{ 
  printf("%s [-r<nnn>] [-a<nnn>] [-e<nnn>] [-p paramter-file] [-OPTION]\n",
	 progname); 
  exit(1); 
}


/*****************************************************************************
*
*  main
*
*****************************************************************************/

int main(int argc, char **argv)
{
  /* Read command line arguments */
  read_command_line(argc,argv);

  /* Read Parameters from parameter file */
  read_parameters();

  /* read box from file header */
  if (box_from_header) read_box(infilename);

  /* Compute parameters, setup colors, rotate box vectors */
  init_ps();

  /* Initialize cell data structures */
  init_cells();

  /* Read atoms */
  read_atoms(infilename);

  /* Sort atoms */
  sort_atoms();

  /* Calculate neighbour tables */
  if ( r_bond != 0.0 )
    do_work(do_neighbour_tables);

  /* Draw the picture */
#ifndef TWOD
  draw_picture3d();
#else
  draw_picture2d();
#endif

  /* Write some data to standard output */
  write_profile();

  /* Write values of options */
  if ( settings == 1 )
      write_settings();

  return 0;
}

/******************************************************************************
*
*  init_ps -- Initialize cutoff radii and color encodings
*
******************************************************************************/

void init_ps(void) {

  int i, j, n = 0;
  real tmp;
  int  tmpcol, tmpcol2, maxtype, maxcol;
  long colcode = 0;

  /* Check parameter values */

  if ( axes != 0 && axes != 1 )
      axes = 0;

  if ( bondbrightness < -1.0 )
      bondbrightness = -1.0;
  else if ( bondbrightness > 1.0 ) 
      bondbrightness = 1.0;

  if ( crop < 0.0 )
      crop = 20.0;

  if ( radii < 0.0 )
      radii = 0.0;

  if ( backgrd < -1.0 )
      backgrd = 0.0;

  if ( ps_height >= 0.0 && ps_height < 1.5 )
      error("Height of picture too small\n");

  if ( ps_width >= 0.0 &&  ps_width < 1.5 )
      error("Width of picture too small\n");

  if ( zoom <= 0.1 )
      zoom = 1.0;

  if ( scaling <= 0.1 )
      scaling = 1.0;

  if ( atombrightness < -1.0 )
      atombrightness = -1.0;
  else if ( atombrightness > 1.0 ) 
      atombrightness = 1.0;

  colorencoding = !( eng == 0 && color_enc == 0 && kineng == 0 );

  /* Cutoff radii for contruction of neighbour tables */
  for (i=0; i<ntypes; i++) {
    for (j=i; j<ntypes; j++) {
      if ( r_cut_vec[0] != -1.0 )
	/* Cutoffs are given as parameter rcut */
	tmp = r_cut_vec[n];
      else if ( ters_r_cut[0] != -1.0 )
	/* Cutoffs are given as Tersoff parameters */
	tmp = ters_r_cut[n];
      else
	/* If cutoffs are not given in the parameter file,
	   take r_max specified by the option -e */
	tmp = r_max;

      r_cut[i][j]  = r_cut[j][i]  = tmp;
      r2cut[i][j] = r2cut[j][i] = r_cut[i][j] * r_cut[i][j];
      ++n;      
    }
  }

  /* Cutoff radius for cell decomposition */
  tmp = 0.0;
  for (i=0; i<ntypes; ++i)
    for (j=0; j<ntypes; ++j)
      tmp = MAX( tmp, r2cut[i][j] );
  r2_cut = MAX(r2_cut,tmp);

  /* Chose parameters in case of wireframes */
  if ( wireframe > 0.0 ) {
    if ( r_atom <= wireframe )
      atomcolor = -999999999;
    r_bond = wireframe;
    bondcolor = -9.0;
    alinew = 0.0;
    blinew = 0.0;
    slinew = 0.0;
  }

  /* Chose parameters for borders of atoms,
     bonds, and bond separators */
  if ( linew >= 0.0 ) {
    alinew = linew;
    blinew = linew;
    slinew = linew;
  }

  /* Linewidths should not change for zoom and scaling */
  if ( zoom != 0.0 && scaling != 0.0 ) {
    alinew       /= zoom * scaling;
    blinew       /= zoom * scaling;
    slinew       /= zoom * scaling;
    frame_linew  /= zoom * scaling;
    rframe_linew /= zoom * scaling;
    sframe_linew /= zoom * scaling;
  }

  /* Specify colors */

  /* User defined color association with types */
  if ( atomcolor != 0 ) {

    tmpcol = atomcolor;
    if ( tmpcol < 0 ) tmpcol = -tmpcol;
    for ( maxtype=0; tmpcol>0; tmpcol/=10, maxtype++ )
      ; 

    tmpcol = atomcolor;
    if ( tmpcol < 0 ) tmpcol = -tmpcol;
    for ( ; maxtype>0; maxtype-- ) {
      tmpcol2 = ( tmpcol % 10 ) - 1;
      atomcol[maxtype-1] = tmpcol2;
      tmpcol /= 10;
    }
  }

  if ( bondcolor != 0 ) {

      tmpcol = bondcolor;
      if ( tmpcol < 0 ) tmpcol = -tmpcol;
      for ( maxtype=0; tmpcol>0; tmpcol/=10, maxtype++ )
	  ; 
      
      tmpcol = bondcolor;
      if ( tmpcol < 0 ) tmpcol = -tmpcol;
      for ( ; maxtype>0; maxtype-- ) {
	  tmpcol2 = ( tmpcol % 10 ) - 1;
	  bondcol[maxtype-1] = tmpcol2;
	  tmpcol /= 10;
      }
  }

  /* Predefined colors */
  /* red */
  color_array[0].x = 1.0;
  color_array[0].y = 0.0;
  color_array[0].z = 0.0;
  /* yellow */
  color_array[1].x = 1.0;
  color_array[1].y = 1.0;
  color_array[1].z = 0.0;
  /* green */
  color_array[2].x = 0.0;
  color_array[2].y = 1.0;
  color_array[2].z = 0.0;
  /* blue */
  color_array[3].x = 0.0;
  color_array[3].y = 0.0;
  color_array[3].z = 1.0;
  /* magenta */
  color_array[4].x = 1.0;
  color_array[4].y = 0.0;
  color_array[4].z = 1.0;
  /* cyan */
  color_array[5].x = 0.0;
  color_array[5].y = 1.0;
  color_array[5].z = 1.0;
  /* orange */
  color_array[6].x = 1.0;
  color_array[6].y = 0.7;
  color_array[6].z = 0.2;
  /* grey */
  color_array[7].x = 0.7;
  color_array[7].y = 0.7;
  color_array[7].z = 0.7;
  /* white */
  color_array[8].x = 1.0;
  color_array[8].y = 1.0;
  color_array[8].z = 1.0;


  /* Predefined shades */
  shade_array[0] = 1.0;
  shade_array[1] = 0.5;
  shade_array[2] = 0.8;
  shade_array[3] = 0.3;
  shade_array[4] = 0.9;
  shade_array[5] = 0.4;
  shade_array[6] = 0.7;
  shade_array[7] = 0.2;
  shade_array[8] = 0.0;

  /* Specify atomic radii */

  /* User defined radius association with types */
  if ( radii > 0 ) {

    for ( i=0; i<9; i++ )
      rad[i] = i;

    tmpcol = radii;
    if ( tmpcol < 0 ) tmpcol = -tmpcol;
    for ( maxtype=0; tmpcol>0; tmpcol/=10, maxtype++ )
      ; 

    tmpcol = radii;
    if ( tmpcol < 0 ) tmpcol = -tmpcol;
    for ( ; maxtype>0; maxtype-- ) {
      tmpcol2 = ( tmpcol % 10 ) - 1;
      rad[maxtype-1] = tmpcol2;
      tmpcol /= 10;
    }

    /* Predefined radii */
    radii_array[0] = 1.0;
    radii_array[1] = 0.5;
    radii_array[2] = 0.75;
    radii_array[3] = 0.125;
    radii_array[4] = 0.625;
    radii_array[5] = 0.375;
    radii_array[6] = 0.875;
    radii_array[7] = 0.125;
    radii_array[8] = 0.0;
  }

  /* Specify color encoding of energy */
  if ( eng != 0 || kineng != 0 || color_enc != 0 ) {

    if ( eng != 0 )
      colcode = eng;
    else if ( kineng != 0 )
      colcode = kineng;
    else if ( color_enc != 0 )
      colcode = color_enc;     

    if ( colcode > 0 && colcode < 9 ) {
      enc_len = 1;
      enc[0] = colcode - 1;
      enc[1] = 7;
    }
    else if ( colcode > 10 ) {
      tmpcol = colcode;
      for ( maxcol=0; tmpcol>0; tmpcol/=10, maxcol++ )
	; 
      enc_len = maxcol-1;
      
      tmpcol = colcode;
      for ( ; maxcol>0; maxcol-- ) {
	tmpcol2 = ( tmpcol % 10 ) - 1;
	enc[maxcol-1] = tmpcol2;
	tmpcol /= 10;
      }
    }
    else if ( colcode < 0 ) {
      if ( colcode == -1 ) {
	enc_len = 5;
	enc[0] = 1;
	enc[1] = 6;
	enc[2] = 2;
	enc[3] = 4;
	enc[4] = 3;
	enc[5] = 5;
      }
      else if ( colcode == -2 ) {
	enc_len = 2;
	enc[0] = 1;
	enc[1] = 6;
	enc[2] = 7;	
      }
      else if ( colcode == -3 ) {
	  enc_len = 3;
	  enc[0] = 1;
	  enc[1] = 7;
	  enc[2] = 4;
	  enc[3] = 3;
      }
      else if ( colcode == -4 ) {
	  enc_len = 3;
	  enc[0] = 1;
	  enc[1] = 7;
	  enc[2] = 6;
	  enc[3] = 2;
      }
      else if ( colcode == -5 ) {
	  enc_len = 2;
	  enc[0] = 1;
	  enc[1] = 5;
	  enc[2] = 3;
      }
      else if ( colcode == -6 ) {
	  enc_len = 2;
	  enc[0] = 2;
	  enc[1] = 4;
	  enc[2] = 3;
      }
    }

    /* Predefined colors */
    /* dark grey */
    enc_array[0].x = 0.3;
    enc_array[0].y = 0.3;
    enc_array[0].z = 0.3;
    /* red */
    enc_array[1].x = 1.0;
    enc_array[1].y = 0.0;
    enc_array[1].z = 0.0;
    /* green */
    enc_array[2].x = 0.0;
    enc_array[2].y = 1.0;
    enc_array[2].z = 0.0;
    /* blue */
    enc_array[3].x = 0.0;
    enc_array[3].y = 0.0;
    enc_array[3].z = 1.0;
    /* cyan */
    enc_array[4].x = 0.0;
    enc_array[4].y = 1.0;
    enc_array[4].z = 1.0;
    /* magenta */
    enc_array[5].x = 1.0;
    enc_array[5].y = 0.0;
    enc_array[5].z = 1.0;
    /* yellow */
    enc_array[6].x = 1.0;
    enc_array[6].y = 1.0;
    enc_array[6].z = 0.0;
    /* white */
    enc_array[7].x = 1.0;
    enc_array[7].y = 1.0;
    enc_array[7].z = 1.0;

  }

      /* Rotations of box vectors */
#ifndef TWOD
  if ( angx != 0.0 ) {
    box_x = xrotate ( box_x, angx );
    box_y = xrotate ( box_y, angx );
    box_z = xrotate ( box_z, angx );
  }
  if ( angy != 0.0 ) {
    box_x = yrotate ( box_x, angy );
    box_y = yrotate ( box_y, angy );
    box_z = yrotate ( box_z, angy );
  }
#endif
  if ( angz != 0.0 ) {
    box_x = zrotate ( box_x, angz );    
    box_y = zrotate ( box_y, angz );    
#ifndef TWOD
    box_z = zrotate ( box_z, angz );
#endif
  }
  
  /* Rotation of unit vectors */
#ifndef TWOD
  if ( angx != 0.0 ) 
    for ( i=0; i<3; i++ )
      unitv[i] = xrotate( unitv[i], angx );
  
  if ( angy != 0.0 ) 
    for ( i=0; i<3; i++ )
      unitv[i] = yrotate( unitv[i], angy );
 
  if ( angz != 0.0 ) 
    for ( i=0; i<3; i++ )
      unitv[i] = zrotate( unitv[i], angz );
#else
  if ( angz != 0.0 ) 
    for ( i=0; i<2; i++ )
      unitv[i] = zrotate( unitv[i], angz );  
#endif
}  

/******************************************************************************
*
*  sort_atoms -- Sort atoms according to increasing distance from focal point
*
******************************************************************************/

void sort_atoms(void) {

  int i, j, l;
#ifndef TWOD
  int k;
#endif
  cell *p;
#ifndef TWOD
  List_elmt *newelmt, *lel; 
  vektor dist;
  real distance, prdistance;
#endif
  int atcount = 0;

  if ( maxl.x == minl.x ) {
    maxl.x += 1.0;
    minl.x -= 1.0;
  }
  if ( maxl.y == minl.y ) {
    maxl.y += 1.0;
    minl.y -= 1.0;
  }
#ifndef TWOD
  if ( maxl.z == minl.z ) {
    maxl.z += 1.0;
    minl.z -= 1.0;
  }
  
  /* Estimate projection parameters */
  if ( zeta == Max ) {
    zeta = ( maxl.z - minl.z );
    printf("Distance of projection plane (T): %.2f (estimated)\n", zeta);
  }
  if ( foc.z == Max ) {
    foc.z = 2 * ( maxl.z - minl.z );
    printf("Distance of center of projection (Z): %.2f (estimated)\n", foc.z); 
  }
 
  /* Check projection paramters */
  if ( foc.z > 0.0 && zeta >= foc.z )
    error("Z must be greater than T!");

  /* Unit vektor from center of system to focal point */
  if ( foc.z > 0.0 ) {
    dist.x = foc.x;
    dist.y = foc.y;
    dist.z = foc.z + (maxl.z - minl.z) / 2.0;
    
    distance = sqrt(SPROD(dist,dist));

    cunit.x = dist.x / distance;
    cunit.y = dist.y / distance;
    cunit.z = dist.z / distance;
  }
  else {
    cunit.x = 0.0;
    cunit.y = 0.0;
    cunit.z = 1.0;
  }
#endif

  /* Move atoms to default position */
  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
#ifndef TWOD
      for (k=0; k < cell_dim.z; ++k)
#endif
      {
#ifdef TWOD
        p = PTR_2D_V(cell_array,i,j  ,cell_dim);
#else
        p = PTR_3D_V(cell_array,i,j,k,cell_dim);
#endif
	for ( l=0; l<p->n; ++l ) {

	  /* Move atoms into default positions */
	  p->ort[l].x -= ( minl.x + maxl.x ) / 2.0;
	  p->ort[l].y -= ( minl.y + maxl.y ) / 2.0;
#ifndef TWOD
	  p->ort[l].z -= maxl.z;

	  /* Distance of atom from focal point */
	  distance = fdistance(p->ort[l]);

	  dist_min = MIN(dist_min, distance);
	  dist_max = MAX(dist_max, distance);

	  prdistance = pdistance(p->ort[l]);


	  /* Create ordered list of atoms (according to distance) */

	  newelmt = (List_elmt *) malloc( sizeof( List_elmt));
	  if ( newelmt == NULL )
	    error("Cannot allocate memory for atom list.");

	  newelmt->cl   = p;
	  newelmt->num  = l;
	  newelmt->d    = prdistance;
	  newelmt->next = NULL;

	  if ( atcount == 0 )
	    atomlist = newelmt;

	  else {

	    if ( atomlist->d > prdistance ) {

	      /* Add new atom at beginning of list */
	      newelmt->next = atomlist;
	      atomlist      = newelmt;

	    }
	    else {

	      for ( lel = atomlist; lel->next != NULL; lel = lel->next ) {

		if ( (lel->next)->d > prdistance )
		  break;
	      }
	      
	      /* Add new atom after listelement */
	      newelmt->next = lel->next;
	      lel->next     = newelmt;

	    }
	  }
#endif	  
	  atcount++;
 
	}
      }
#ifndef TWOD
  for ( lel = atomlist; lel != NULL; lel = lel->next ) {
      
      p = lel->cl;
      i = lel->num;
      distance = lel->d;
      
      /* Position of atom (real) */
      dist = p->ort[i];
  }
#endif

  /* Origin of real bounding box */
  ursprung = real_minl;
#ifndef TWOD
  if ( angx != 0.0 )
      ursprung = xrotate( ursprung, angx );
  if ( angy != 0.0 )
      ursprung = yrotate( ursprung, angy );
#endif
  if ( angz != 0.0 )
      ursprung = zrotate( ursprung, angz );
  
  ursprung.x -= ( minl.x + maxl.x ) / 2.0;
  ursprung.y -= ( minl.y + maxl.y ) / 2.0;
#ifndef TWOD
  ursprung.z -= maxl.z;
#endif
}

/******************************************************************************
*
*  do_neighbour_tables -- Calculate neighbour tables.
*
******************************************************************************/

void do_neighbour_tables(cell *p, cell *q, vektor pbc) {

  int i, j;
#ifndef TWOD
  int k, l;
#endif
  int jstart;
  int q_typ, p_typ;
  vektor d, tmp_d;
  real radius;
#ifndef TWOD
  cell *neigh_cell;
  real distance, ndistance, kdistance;
#endif

  /* For each atom in first cell */
  for (i=0; i<p->n; ++i) {

    tmp_d.x = p->ort[i].x - pbc.x;
    tmp_d.y = p->ort[i].y - pbc.y;
#ifndef TWOD
    tmp_d.z = p->ort[i].z - pbc.z;

    distance = pdistance(p->ort[i]);
#endif
    p_typ   = p->sorte[i];

    jstart = ( p==q ? i+1 : 0 );
    
    /* For each atom in neighbouring cell */
    for (j = jstart; j < q->n; ++j) {

      q_typ = q->sorte[j];
	  
      /* Calculate distance  */
      d.x = q->ort[j].x - tmp_d.x;
      d.y = q->ort[j].y - tmp_d.y;
#ifndef TWOD
      d.z = q->ort[j].z - tmp_d.z;

      ndistance = pdistance(q->ort[j]);

#endif
      radius = sqrt(SPROD(d,d));

      /* Make neighbour tables */
      if (radius <= r_cut[p_typ][q_typ]) { 
      
	neightab *neigh;
	int newsz;

	nbonds++;
	  
	/* Update neighbour table of particle i */
	neigh = p->neightab_array + i;

	if (neigh->n_max <= neigh->n ) {
	  newsz = neigh->n_max + NSTEP;

	  neigh->typ = (short *) realloc(neigh->typ, newsz * sizeof(short) );
	  neigh->cl  = (void **) realloc(neigh->cl,  newsz * sizeof(cellptr) );
	  neigh->num = (int *)   realloc(neigh->num, newsz * sizeof(int) );
	  
	  if ( neigh->typ==NULL || neigh->cl==NULL || neigh->num==NULL )
	    error("Cannot allocate memory for neighbour tables!");
	  
	  neigh->n_max += NSTEP;
	}
#ifndef TWOD
	/* Sort neighbours by distance */
	if (neigh->n == 0) {
	  neigh->typ[0] = q_typ;
	  neigh->cl [0] = q;
	  neigh->num[0] = j;
	}
	else {
	  for ( k=0; k<neigh->n; k++ ) {
	    neigh_cell = neigh->cl[k];

	    kdistance = pdistance(neigh_cell->ort[neigh->num[k]]);

	    if ( kdistance > ndistance )
	      break;
	  }
	  
	  for ( l=neigh->n; l>k; l-- ) {
	    neigh->typ[l] = neigh->typ[l-1];
	    neigh->cl [l] = neigh->cl [l-1];
	    neigh->num[l] = neigh->num[l-1];
	  }
	  
	  neigh->typ[k] = q_typ;
	  neigh->cl [k] = q;
	  neigh->num[k] = j;
	}
#else	  
	neigh->typ[neigh->n] = q_typ;
	neigh->cl [neigh->n] = q;
	neigh->num[neigh->n] = j;
#endif
	neigh->n++;
	
	/* Update neighbour table of particle j */
	neigh      = q->neightab_array + j;
	
	if (neigh->n_max <= neigh->n ) {
	  newsz = neigh->n_max + NSTEP;
	  
	  neigh->typ = (short *) realloc(neigh->typ, newsz * sizeof(short) );
	  neigh->cl  = (void **) realloc(neigh->cl,  newsz * sizeof(cellptr) );
	  neigh->num = (int *)   realloc(neigh->num, newsz * sizeof(int) );
	  
	  if ( neigh->typ==NULL || neigh->cl==NULL || neigh->num==NULL )
	    error("Cannot allocate memory for neighbour tables!");
	  
	  neigh->n_max      += NSTEP;
	}
#ifndef TWOD
	/* Sort neighbours by distance */
	if (neigh->n == 0) {
	  neigh->typ[0] = p_typ;
	  neigh->cl [0] = p;
	  neigh->num[0] = i;
	}
	else {
	  for ( k=0; k<neigh->n; k++ ) {
	    neigh_cell = neigh->cl[k];

	    kdistance = pdistance(neigh_cell->ort[neigh->num[k]]);

	    if ( kdistance > distance )
	      break;
	  }
	  
	  for ( l=neigh->n; l>k; l-- ) {
	    neigh->typ[l] = neigh->typ[l-1];
	    neigh->cl [l] = neigh->cl [l-1];
	    neigh->num[l] = neigh->num[l-1];
	  }
	  
	    neigh->typ[k] = p_typ;
	    neigh->cl [k] = p;
	    neigh->num[k] = i;
	}
#else	    
	neigh->typ[neigh->n] = p_typ;
	neigh->cl [neigh->n] = p;
	neigh->num[neigh->n] = i;
#endif
	neigh->n++;
      }
    } /* for j */
  } /* for i */

}

/******************************************************************************
*
*  write_profile -- Write some data to standard output
*
******************************************************************************/

void write_profile(void) {

  printf("\n");
  printf("Total number of atoms: %d\n", natoms);
  printf("Total number of bonds: %d\n", nbonds);
  printf("Number of displayed atoms: %d\n", natoms_drawn);
#ifdef TWOD
  printf("Number of displayed bonds: %d\n", nbonds_drawn / 2);
#else
  printf("Number of displayed bonds: %d\n", nbonds_drawn);
#endif
}

/******************************************************************************
*
*  write_settings -- Write settings of options to .opt file
*
******************************************************************************/

void write_settings(void) {

    FILE *opt;
    str255 fname;

    /* Open settigs file */
    sprintf(fname, "%s.ps.opt", infilename);
    opt = fopen(fname,"w");
    if (NULL == opt) 
	error("Cannot open settings file.\n");

    fprintf( opt, "imd_ps -p%s ", paramfilename);
    if ( restart != -1 )
	fprintf( opt, "-r%d ", restart);
    else if ( avpos != -1 )
	fprintf( opt, "-A%d ", avpos);
    if ( r_cut_vec[0] == -1 && ters_r_cut[0] == -1 )
	fprintf( opt, "-e%f ", r_max);
    if ( use_vtypes == 1 )
	fprintf( opt, "-v ");
    if ( settings ==1 )
	fprintf( opt, "-K ");
    fprintf( opt, "-a%d -b%f -B%f -c%f -C%ld -D%ld ", 
	     axes, bondbrightness, r_bond, crop, atomcolor, radii);
    fprintf( opt, "-E%ld -f%f -F%f -g%f -G%ld -h%f -H%f ",
	     eng, frame, frame_linew, backgrd, color_enc, ps_height, linew);
    fprintf( opt, "-i%f -I%f -j%f -J%f -k%f -K -l%f -L%f -m%f ",
	     translation.x, ambient, translation.y, atombrightness, spect, 
	     blinew, alinew, zoom);
    fprintf( opt, "-o%f -O%f -P%f -q%f -Q%f -R%f ", 
	     sframe, sframe_linew, slinew, frame_border, 
	     frame_border_linew, r_atom);
    fprintf( opt, "-s%f -S%f -t%ld -u%f -U%f -V%ld ",
	     scaling, shading, kineng, rframe, rframe_linew, bondcolor);
    fprintf( opt, "-w%f -W%f -z%f ",
	     ps_width, wireframe, angz);
#ifndef TWOD
    fprintf( opt, "-d%f -M%f -n%f -N%f -T%f ", 
	     depth, user_scale, llinew, needle, zeta);
    fprintf( opt, "-x%f -y%f -X%f -Y%f -Z%f ", 
	     angx, angy, foc.x, foc.y, foc.z);
#endif

    fclose(opt);
}

#ifndef TWOD

/******************************************************************************
*
*  fdistance -- Compute distance of atoms from focal point
*
******************************************************************************/

real fdistance(vektor pos) {

  pos.x -= foc.x;
  pos.y -= foc.y;
  pos.z -= foc.z;
  
  return sqrt( SPROD(pos,pos) );
}

/******************************************************************************
*
*  pdistance -- Compute length of projection of position vector onto
*               vector connecting the focal point and the system center
*
******************************************************************************/

real pdistance(vektor pos) {

  pos.z += ( maxl.z - minl.z ) / 2.0;

  return SPROD(pos,cunit);
}

/******************************************************************************
*
*  scale -- Compute scale factor with depth
*
******************************************************************************/

real scale(vektor pos, real sclfct) {

  real tmp, tmpx, tmpy;

  if ( user_scale < 0.0 ) user_scale = -user_scale;

  if ( foc.z > 0 ) {
    tmp  = ( pos.z * user_scale - foc.z );
    tmpx = SQR( ( pos.x - foc.x ) / tmp );
    tmpy = SQR( ( pos.y - foc.y ) / tmp );

    /* Provisorisch */
    if ( tmpx > 0.1 || tmpy > 0.1 ) {
      tmpx = 0.1;
      tmpy = 0.1;
    }
    tmpx = 0;
    tmpy = 0;

    return sclfct * ( zeta - foc.z ) / tmp 
      * ( 1.0 - tmpx - tmpy );
  } 
  else
    return sclfct;
}

/******************************************************************************
*
*  proj -- Compute perspective projection
*
******************************************************************************/

vektor2d proj(vektor pos, real sclfct, vektor2d trans) {

  vektor2d tmp;

  if ( foc.z > 0 ) {

      tmp.x = sclfct * ( foc.x + ( pos.x - foc.x ) 
	    * ( zeta - foc.z ) / ( pos.z - foc.z ) ) + trans.x;
      tmp.y = sclfct * ( foc.y + ( pos.y - foc.y ) 
	    * ( zeta - foc.z ) / ( pos.z - foc.z ) ) + trans.y;
  }
  else {

    tmp.x = sclfct * pos.x + trans.x;
    tmp.y = sclfct * pos.y + trans.y;
  }
  return tmp;
}

#endif

/******************************************************************************
*
*  transform -- Transformation of points in the projection plane
*
******************************************************************************/

vektor2d transform(vektor2d vekt, real scl, vektor2d trans) {

  vektor2d tmpv;

  tmpv.x = scl * vekt.x + trans.x;
  tmpv.y = scl * vekt.y + trans.y;

  return tmpv;
}

/******************************************************************************
*
*  printcolor -- Write rgb values of colors
*
******************************************************************************/

void printcolor(FILE *out, real type, real dist, int obj) {

    /* obj = 0 : atom          *
     * obj = 1 : bond          *
     * obj = 3 : background    *
     * obj = 4 : box           */

  real tmp = 1.0, tmpr, val;
  int i;
  real  shade = 1.0;
  vektor3d tmpv = {0, 0, 0};
  real dist_norm;
  int color = 0; 

  /* Atom or bond color */
  if ( obj == 0 || obj == 1 ) {

      if ( dist_max == dist_min ) 
	  dist_norm = 1.0;
      else
	  dist_norm = ( dist - dist_min ) / ( dist_max - dist_min );

#ifndef TWOD
      if ( depth >= 0.0 && depth <= 1.0 )
	  tmp = 1.0 - depth * dist_norm;
#endif

      if ( obj == 0 || bondcolor == 0 ) {

	  if ( colorencoding ) {
	      if ( ( maxeng - mineng ) > 0.0 )
		  val = ( type - mineng ) / ( maxeng - mineng ) * enc_len;
	      else
		  val = enc_len;
	      
	      for ( i=0; i<enc_len; i++ )
		  if ( val >= i && val < (i+1) ) {
		      tmpv.x = ( i+1 - val ) * enc_array[enc[i]].x 
			  - ( i - val ) * enc_array[enc[i+1]].x;
		      tmpv.y = ( i+1 - val ) * enc_array[enc[i]].y 
			  - ( i - val ) * enc_array[enc[i+1]].y;
		      tmpv.z = ( i+1 - val ) * enc_array[enc[i]].z 
			  - ( i - val ) * enc_array[enc[i+1]].z;
		  }
	      if ( val == enc_len ) {
		  tmpv.x = enc_array[enc[enc_len]].x;
		  tmpv.y = enc_array[enc[enc_len]].y;
		  tmpv.z = enc_array[enc[enc_len]].z;
	      }
	  }
	  else {

	      if ( atomcolor < 0 ) {
		  tmpr = shade_array[atomcol[(int)type]];
		  
		  tmpv.x = tmpr;
		  tmpv.y = tmpr;
		  tmpv.z = tmpr;
	      }
	      else if ( atomcolor > 0 )
		  tmpv = color_array[atomcol[(int)type]];
	  }
      }
      else if ( obj == 1 && bondcolor != 0 ) {
	  if ( bondcolor < 0 ) {
	      tmpr = shade_array[bondcol[(int)type]];
	      
	      tmpv.x = tmpr;
	      tmpv.y = tmpr;
	      tmpv.z = tmpr;
	  }
	  else if ( bondcolor > 0 )
	      tmpv = color_array[bondcol[(int)type]];
      }

      if ( obj == 0 ) {
	  if ( atombrightness >= 0.0 )
	      fprintf(out, "%.2f %.2f %.2f ",
		      ( 1.0 - atombrightness ) * tmpv.x * tmp + atombrightness, 
		      ( 1.0 - atombrightness ) * tmpv.y * tmp + atombrightness, 
		      ( 1.0 - atombrightness ) * tmpv.z * tmp + atombrightness );	
	  else
	      fprintf(out, "%.2f %.2f %.2f ",
		      ( 1.0 + atombrightness ) * tmpv.x * tmp, 
		      ( 1.0 + atombrightness ) * tmpv.y * tmp, 
		      ( 1.0 + atombrightness ) * tmpv.z * tmp );			  
      }    
      else if ( obj == 1 ) {
	  if ( bondbrightness >= 0.0 )
	      fprintf(out, "%.2f %.2f %.2f ",
		      ( 1.0 - bondbrightness ) * tmpv.x * tmp + bondbrightness, 
		      ( 1.0 - bondbrightness ) * tmpv.y * tmp + bondbrightness, 
		      ( 1.0 - bondbrightness ) * tmpv.z * tmp + bondbrightness );	
	  else
	      fprintf(out, "%.2f %.2f %.2f ",
		      ( 1.0 + bondbrightness ) * tmpv.x * tmp, 
		      ( 1.0 + bondbrightness ) * tmpv.y * tmp, 
		      ( 1.0 + bondbrightness ) * tmpv.z * tmp );		
      }
  }
    
  /* Background color */
  else if ( obj == 3 ) {

    if ( backgrd < 0.0 && backgrd >= -1.0 )
      fprintf(out, "%.2f %.2f %.2f ", 
	      1.0 + backgrd,
	      1.0 + backgrd,
	      1.0 + backgrd); 
    else {
      color    = (int)backgrd;
      shade = 1.0 - backgrd + color;
    
      if ( color == 0 ) {
	tmpv.x = 1;
	tmpv.y = 1;
	tmpv.z = 1;
      }
      else if ( color < 10 )
	tmpv = color_array[color-1];
    
      fprintf(out, "%.2f %.2f %.2f ", 
	      tmpv.x * shade, 
	      tmpv.y * shade, 
	      tmpv.z * shade); 
    }
  }

  /* Color of boxes */
  else if ( obj == 4 ) {

    if ( type == 0 ) {
      color = (int)frame;
      shade = 1.0 - frame + color;
    }
    else if ( type == 1 ) {
      color = (int)sframe;
      shade = 1.0 - sframe + color;
    } 
    else if ( type == 2 ) {
      color = (int)rframe;
      shade = 1.0 - rframe + color;
    } 
    else if ( type == 3 ) {
      color = (int)frame_border;
      shade = 1.0 - frame_border + color;
    } 	

    tmpv = color_array[color-1];

    fprintf(out, "%.2f %.2f %.2f ", 
	    tmpv.x * shade, 
	    tmpv.y * shade, 
	    tmpv.z * shade); 
  }
}

/******************************************************************************
*
*  setlinew -- Set the linewidth
*
******************************************************************************/

void setlinew(FILE *out, real linew) {

  if ( linew > 0.0 )
    fprintf(out, "[] 0 setdash %.2f setlinewidth ", linew); 
  else
    fprintf(out, "[4 4] 0 setdash %.2f setlinewidth ", -linew);       
}

/******************************************************************************
*
*  setradius -- Set the radius of the atoms
*
******************************************************************************/

real setradius(int type) {

  if ( radii > 0 )
    return radii_array[rad[type]];
  else
    return 1.0;
}

/******************************************************************************
*
*  draw_colorencoding -- Draw spectrum of colors of colorencoding
*
******************************************************************************/

void draw_colorencoding(FILE *out, vektor2d position, vektor2d size) {

  int res = (int)size.x;
  int i, j;
  int col[3];
  real val;

  fprintf(out, "gsave newpath\n");
  fprintf(out, "%.2f %.2f translate %.2f %.2f scale\n", 
	  position.x, position.y, size.x, size.y );

  fprintf(out, "%d 1 8 [%d 0 0 -1 0 1]\n", res+1, res+1);
  fprintf(out, "{<");

  for( j=0; j<=res; j++) {
    val = (real) j * enc_len / res;

    for ( i=0; i<enc_len; i++ )
      if ( val >= i && val < (i+1) ) {
	col[0] = (int) (255 * ( ( i+1 - val ) * enc_array[enc[i]].x 
				- ( i - val ) * enc_array[enc[i+1]].x ) );
	col[1] = (int) (255 * ( ( i+1 - val ) * enc_array[enc[i]].y 
				- ( i - val ) * enc_array[enc[i+1]].y ) );
	col[2] = (int) (255 * ( ( i+1 - val ) * enc_array[enc[i]].z 
				- ( i - val ) * enc_array[enc[i+1]].z ) );
      }
    if ( val == enc_len ) {
      col[0] = (int) 255 * enc_array[enc[enc_len]].x;
      col[1] = (int) 255 * enc_array[enc[enc_len]].y;
      col[2] = (int) 255 * enc_array[enc[enc_len]].z;
    }    

    for( i=0; i<3; i++) {
      if ( col[i] < 16 )
	fprintf(out, "0%x", col[i]);
      else
	fprintf(out, "%x", col[i]);
    }
    if ( ( 6 * j ) % 120 == 0)
	fprintf(out, "\n");
  }
  
  fprintf(out, "\n>}\nfalse 3 colorimage grestore\n");
  
  /* Draw box */
  fprintf(out, "gsave newpath %.2f %.2f moveto\n", 
	  position.x, position.y);
  fprintf(out, "%.2f %.2f lineto\n", 
	  position.x + size.x, position.y);  
  fprintf(out, "%.2f %.2f lineto\n", 
	  position.x + size.x, position.y + size.y);
  fprintf(out, "%.2f %.2f lineto\n", 
	  position.x, position.y + size.y);
  fprintf(out, "closepath 2 setlinewidth 0 setgray stroke grestore\n");
}

/******************************************************************************
*
*  Rotations of vectors
*
******************************************************************************/

#ifndef TWOD
vektor xrotate(vektor pos, real angle) {

  vektor tmp_pos;

  tmp_pos.x = pos.x;
  tmp_pos.y = cos( PIN * angx ) * pos.y - sin( PIN * angx ) * pos.z;
  tmp_pos.z = sin( PIN * angx ) * pos.y + cos( PIN * angx ) * pos.z;

  return tmp_pos;
}

vektor yrotate(vektor pos, real angle) {

  vektor tmp_pos;

  tmp_pos.x = sin( PIN * angy ) * pos.z + cos( PIN * angy ) * pos.x;
  tmp_pos.y = pos.y;
  tmp_pos.z = cos( PIN * angy ) * pos.z - sin( PIN * angy ) * pos.x;

  return tmp_pos;
}
#endif

vektor zrotate(vektor pos, real angle) {

  vektor tmp_pos;

	tmp_pos.x = cos( PIN * angz ) * pos.x - sin( PIN * angz ) * pos.y; 
	tmp_pos.y = sin( PIN * angz ) * pos.x + cos( PIN * angz ) * pos.y;
#ifndef TWOD 
	tmp_pos.z =  pos.z;
#endif

  return tmp_pos;
}

/******************************************************************************
*
*  Maximum and minimum of vectors
*
******************************************************************************/

vektor maxvektor(vektor v1, vektor v2) {

  vektor tmpv;

  tmpv.x = MAX( v1.x, v2.x );
  tmpv.y = MAX( v1.y, v2.y );
#ifndef TWOD
  tmpv.z = MAX( v1.z, v2.z );
#endif

  return tmpv;
}

vektor minvektor(vektor v1, vektor v2) {

  vektor tmpv;

  tmpv.x = MIN( v1.x, v2.x );
  tmpv.y = MIN( v1.y, v2.y );
#ifndef TWOD
  tmpv.z = MIN( v1.z, v2.z );
#endif

  return tmpv;
}

/******************************************************************************
*
*  scl -- Transformations of points
*
******************************************************************************/

#ifdef TWOD

vektor2d scl(vektor pos, real sclfct, vektor2d trans) {

  vektor2d tmp;
  
  tmp.x = pos.x * sclfct + trans.x;
  tmp.y = pos.y * sclfct + trans.y;

  return tmp;
}

#endif
