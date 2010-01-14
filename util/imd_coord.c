
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
* imd_coord -- calculate coordination numbers
*
* A descendant of imd_pair
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#ifndef COORD
#define COORD
#endif

#ifdef TERSOFF
#ifdef TWOD
#undef TWOD
#endif
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
  printf("%s [-r<nnn>] [-A<nnn>] [-e<nnn>] [-C<nnn>] [-l] [-g] [-u] [-p paramter-file]\n",progname); 
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

  /* Initializations */
  init_coord();

#ifdef TERSOFF
  init_tersoff();
#endif

  /* Initialize cell data structures */
  init_cells();

  /* Read atoms */
  read_atoms(infilename);

  /* Do the work */
  do_work(do_cell_pair);

  /* Output results */
  write_data();

  return 0;
}

/******************************************************************************
*
*  init_coord -- Initializations
*
******************************************************************************/

void init_coord()
{
  int tablesize;
  int i, j, n=0;
  real tmp;

  /* Allocate and initialize table for bond occurences */
  tablesize = ntypes*ntypes*sizeof(real);
  numbers = (real *) malloc(tablesize);
  if (NULL==numbers) error("Cannot allocate memory for numbers.");
  num_dim.x = ntypes;
  num_dim.y = ntypes;

  for (i=0; i<ntypes; ++i)
    for (j=0; j<ntypes; ++j)
	*PTR_2D_V(numbers,i,j,num_dim) = 0.0;

  /* Cutoff radii for contruction of neighbour tables */
  for (i=0; i<ntypes; i++) {
    for (j=i; j<ntypes; j++) {
      /* imd_coord should behave as written on web site, i.e take cut off */
      /* from command line, not parameter file */
/*       if ( r_cut_vec[0] != -1.0 ) */
/* 	/\* Cutoffs are given as parameter r_cut *\/ */
/* 	tmp = r_cut_vec[n]; */
/*       else */
/* 	/\* If cutoffs are not given in the parameter file, */
/* 	   take r_max specified by the option -e *\/ */
      tmp = r_max;

      r_cut[i][j]  = r_cut[j][i]  = tmp;
      ++n;      
    }
  }

  /* Cutoff radius for cell decomposition */
  tmp = 0.0;
  for (i=0; i<ntypes; ++i)
    for (j=0; j<ntypes; ++j)
      tmp = MAX( tmp, r_cut[i][j]*r_cut[i][j] );
  r2_cut = MAX(r2_cut,tmp);
}

/******************************************************************************
*
*  write_data writes numbers to *.coord file
*
******************************************************************************/

void write_data()
{
  FILE *out;
  str255 fname;
  int  i, j, l, m;
#ifndef TWOD
  int k;
#endif
  cell *p;
  real ttl_coord;
  int  *spec;

  spec = (int *) malloc( (c_max + 1) * sizeof(int));
  if ( spec == NULL ) error("Cannot allocate memory for coordination numbers");

  /* Write local coordination numbers in .coord file */
  if ( local == 1 ) {

    if (-1==restart)
      sprintf(fname,"%s.coord",infilename);
    else
      sprintf(fname,"%s.%05d.coord",outfilename,restart);

    out = fopen(fname,"w");
    if (NULL == out) error("Cannot open outfile.");
  }

  for ( i=0; i<= c_max; i++ ) 
    spec[i] = 0;

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
	  
	  if ( local == 1 ) {

#ifdef TWOD
	    fprintf(out,"%d %d %f %f ", p->nummer[l], p->sorte[l], 
		    p->ort[l].x, p->ort[l].y );
#else
	    fprintf(out,"%d %d %f %f %f ", p->nummer[l], p->sorte[l], 
		    p->ort[l].x, p->ort[l].y, p->ort[l].z );
#endif

	  }

	  /* Coordination number for all types */
	  ttl_coord = 0.0;
	  for ( m=0; m<ntypes; m++ )
	    ttl_coord += p->coord[ l * ntypes + m ];

	  if ( local == 1 ) fprintf(out, "%f ", ttl_coord);

	  /* Spectrum of coordination numbers */
	  if ( ttl_coord <= c_max ) 
	  ++spec[(int)(rint(ttl_coord))];

	  if ( local == 1 ) {

	    /* Coordination number for each type */
	    for ( m=0; m<ntypes; m++ )
	      fprintf(out, "%f ", p->coord[ l * ntypes + m ]);

	    if ( write_poteng == 1 )
	      fprintf(out, "%f", p->poteng[l]);

	    fprintf(out, "\n");
	  }
	}
      }

  if ( local ==1 ) 
    fclose(out);

  /* Write global coordination numbers in .gcoord file */
  if ( global == 1 ) {

    if (-1==restart)
      sprintf(fname,"%s.gcoord",infilename);
    else
      sprintf(fname,"%s.%05d.gcoord",outfilename,restart);

    out = fopen(fname,"w");
    if (NULL == out) error("Cannot open outfile.");

    for ( i=0; i<c_max; ++i )
      fprintf(out, "%d\t%f\n", i, (real) spec[i]/natoms);

    fclose(out);

  }

  /* Write global data to standard output */
  printf("\n");
  if ( bonds == 0 )
    printf("No bonds found.\n\n");
  printf("---------------------\n");
  printf(" Bond   Occurence\n\n");

  for (i=0; i<ntypes; ++i)
    for (j=i; j<ntypes; ++j) {
      if ( bonds > 0 )
	printf("  %d%d     %f %%\n",
	       i, j, *PTR_2D_V(numbers,i,j,num_dim)/bonds*100);
      else
	printf("  %d%d     0.00 %%\n", i, j);
    }

  printf("\n");
  printf("---------------------\n");
  printf(" Coordination    Occurence\n\n");
  
  for ( i=0; i<=c_max; ++i )
    printf("  %d               %2.2f %%\n", i, (real) spec[i]/natoms*100);

  printf("\n");
  printf("------------------------------\n");  

}


/******************************************************************************
*
*  do_cell_pair calulates the coordination numbers for atoms in two cells
*
******************************************************************************/

void do_cell_pair(cell *p, cell *q, vektor pbc)
{
  int i,j,k;
  int temp;
  vektor d;
  real radius;
  int p_typ,q_typ;
#ifdef TERSOFF
  real fc;
#endif

  /* For each atom in first cell */
  for (i = 0;i < p->n; ++i) 
    /* For each atom in neighbouring cell */
    /* Nasty little trick: If p==q, use only rest of atoms */
    for (j = ((p==q) ? i+1 : 0); j < q->n; ++j) {
      
      /* Calculate distance */
      d.x =q->ort[j].x - p->ort[i].x + pbc.x;
      d.y =q->ort[j].y - p->ort[i].y + pbc.y;
#ifndef TWOD
      d.z =q->ort[j].z - p->ort[i].z + pbc.z;
#endif

      radius = sqrt( (double)(SPROD(d,d)) );
      p_typ = p->sorte[i];
      q_typ = q->sorte[j];

#ifndef TERSOFF
      if ( radius <= r_cut[p_typ][q_typ] && radius > r_min ) {

	/* Compute local coordination numbers */
	for ( k=0; k<ntypes; k++ ) {
	  if ( q_typ == k ) 
	    ++(p->coord[i * ntypes + k]);
	  if ( p_typ == k )
	    ++(q->coord[j * ntypes + k]);
	}

	/* Compute average coordination number for each bond type */
	if (q_typ > p_typ) {temp = p_typ; p_typ = q_typ; q_typ = temp;}
      
	*PTR_2D_V(numbers, q_typ, p_typ, num_dim) += 1.0;

	++bonds;
      }
#else
      if ( radius <= ter_r_cut[p_typ][q_typ] ) {

	if ( use_unity == 1 )
	  fc = 1.0;
	else {

	  /* Use Cutoff function of Tersoff potential */
	  if ( radius <= ter_r0[p_typ][q_typ] ) 
	    fc = 1.0;
	  else 
	    fc = 0.5 * ( 1.0 + cos( 3.141592654 * ( radius - ter_r0[p_typ][q_typ] ) 
				    / ( ter_r_cut[p_typ][q_typ] - ter_r0[p_typ][q_typ] )  ) );
	}

	/* Compute local coordination numbers */
	for ( k=0; k<ntypes; k++ ) {
	  if ( q_typ == k ) 
	    p->coord[i * ntypes + k] += fc;
	  if ( p_typ == k )  
	    q->coord[j * ntypes + k] += fc;
	}

	/* Compute average coordination number for each bond type */
	if (q_typ > p_typ) {temp = p_typ; p_typ = q_typ; q_typ = temp;}
      
	*PTR_2D_V(numbers, q_typ, p_typ, num_dim) += fc;

	++bonds;
      }
#endif
    }    
}

#ifdef TERSOFF

/******************************************************************************
*
*  init_tersoff
*
******************************************************************************/

void init_tersoff(void) {

  int  i, j, n = 0;
  real tmp;

  /* parameters for more than one atom type */
  for (i=0; i<ntypes; i++) 
    for (j=i; j<ntypes; j++) {
      ter_r_cut[i][j]  = ter_r_cut[j][i]  = ters_r_cut[n];
      ter_r2_cut[i][j] = ter_r2_cut[j][i] = ter_r_cut[i][j] * ter_r_cut[i][j];
      ter_r0[i][j]     = ter_r0[j][i]     = ters_r0[n];
    }

  /* Determine cutoff for cell decomposition */
  tmp = 0.0;
  for (i=0; i<ntypes; ++i)
    for (j=0; j<ntypes; ++j)
      tmp = MAX( tmp, ter_r2_cut[i][j] );
  r2_cut = MAX(r2_cut,tmp);
  
}

#endif








