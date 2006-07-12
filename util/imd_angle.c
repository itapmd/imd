
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
* imd_angle -- calculate angular distribution functions
* A descendant of imd_pair
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#ifndef ANGLE
#define ANGLE
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
  printf("%s [-r<nnn>] [-A<nnn>] [-a<nnn>] [-e<nnn>] [-v] [-p paramter-file]\n",progname); 
  exit(1); 
}

/*****************************************************************************
*
*  main
*
*****************************************************************************/

int main(int argc, char **argv)
{
  int tablesize;
  int i,j,k,l;

  /* Read command line arguments */
  read_command_line(argc,argv);

  /* Read Parameters from parameter file */
  read_parameters();

  tablesize = slots*ntypes*ntypes*ntypes*sizeof(real);
  histogram = (real *) malloc(tablesize);
  if (NULL==histogram) error("Cannot allocate memory for histograms.");
  hist_dim.i = slots;
  hist_dim.x = ntypes;
  hist_dim.y = ntypes;
  hist_dim.z = ntypes;

  for (i=0; i<slots; ++i)
    for (j=0; j<ntypes; ++j)
      for (k=0; k<ntypes; ++k)
	for (l=0; l<ntypes; ++l)
	  *PTR_4D_V(histogram,i,j,k,l,hist_dim) = 0.0;

  r2_cut = SQR(r_max);

  /* read box from file header */
  if (box_from_header) read_box(infilename);

  /* Initialize cell data structures */
  init_cells();

  /* Read atoms */
  read_atoms(infilename);

  /* Calculate the angles */
  calc_angles();

  /* Output results */
  write_data();

  return 0;

}

/******************************************************************************
*
*  write_data writes histogram to *.angle file
*
******************************************************************************/

void write_data()
{
  FILE *out;
  str255 fname;
  int i,j,k,l;
  real phi;

  if (-1==restart)
    sprintf(fname,"%s.angle",infilename);
  else
    sprintf(fname,"%s.%05d.angle",outfilename,restart);

  out = fopen(fname,"w");
  if (NULL == out) error("Cannot open histograms file.");

  for (i=1; i<slots; ++i) {
    phi = ((float) i / slots * 180);
    fprintf(out,"%f ", phi);
    for (j=0; j<ntypes; ++j)
      for (k=0; k<ntypes; ++k)
	for (l=k; l<ntypes; ++l)
	  /* j is the middle atom */
	  if (nangles != 0)
	    fprintf(out,"%f ",*PTR_4D_V(histogram,i,j,k,l,hist_dim)/nangles);
          else 
	    error("No neighbouring atoms in this distance range.");
    fprintf(out,"\n");
  };

  printf("Minimal radius: %.6f\n", r_min);
  printf("Maximal radius: %.6f\n", r_max);
  printf("%d angles computed\n", nangles);
  printf("Histogram with %d slots\n", slots);

  fclose(out);

}

/******************************************************************************
*
*  do_angle calulates the angles for atoms in three cells
*
******************************************************************************/

void do_angle(cell *p, cell *q, cell *qq, vektor pbc, vektor ppbc)
{
  int i,j,k;
  int ang;
  int temp;
  vektor d, dd;
  real radius, rradius;
  real phi;
  int p_typ, q_typ, qq_typ;

  /* For each atom in first cell */
  for ( i=0; i<p->n; ++i) 
    /* For each atom in neighbouring cell */
    for ( j=0; j<q->n; ++j) 
      /* For each atom in second neighbouring cell */
      for ( k=((q==qq) ? j+1 : 0); k<qq->n; ++k) {
      
	/* Calculate distance vectors */
	d.x = q->ort[j].x - p->ort[i].x + pbc.x;
	d.y = q->ort[j].y - p->ort[i].y + pbc.y;
#ifndef TWOD
	d.z = q->ort[j].z - p->ort[i].z + pbc.z;
#endif

	dd.x = qq->ort[k].x - p->ort[i].x + ppbc.x;
	dd.y = qq->ort[k].y - p->ort[i].y + ppbc.y;
#ifndef TWOD
	dd.z = qq->ort[k].z - p->ort[i].z + ppbc.z;
#endif

	radius  = sqrt( (double)(SPROD(d,d)) );
	rradius = sqrt( (double)(SPROD(dd,dd)) );

	/* Calculate angles */
	if ( (radius < r_max) && (rradius < r_max) 
	  && (radius > r_min) && (rradius > r_min) ) {
	  ++nangles;
	  phi = (double) (acos( (double)(SPROD(d,dd) )/ (radius * rradius)) );
	  
	  ang = (int) ( slots * phi / 3.141592654 );
	  
	  p_typ  = p->sorte[i];
	  q_typ  = q->sorte[j];
	  qq_typ = qq->sorte[k];

	  if (q_typ > qq_typ) {temp = qq_typ; qq_typ = q_typ; q_typ = temp;}; 

	  if ((ang>0) && (ang<slots))
	    ++*PTR_4D_V(histogram, ang , p_typ, q_typ, qq_typ, hist_dim);
	} 
    };

}

/******************************************************************************
*
*  calc_angles calulates the angles for all atoms
*
******************************************************************************/

void calc_angles(void)
{
  cell *p,*q,*qq;
  int i,j,k;
  int l,m,n;
  int f,g,h;
  int r,s,t;
  int u,v,w;
  vektor pbc,ppbc;

  /* for each cell */
  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
#ifndef TWOD
      for (k=0; k < cell_dim.z; ++k)
#endif
	/* For the neighbours of this cell */
	for (l=-1; l <= 1; ++l)
	  for (m=-1; m <= 1; ++m)
#ifndef TWOD
	    for (n=-1; n <= 1; ++n)
#endif
	      /* For all the following neighbouring cells */
	      for (f=-1; f <= 1; ++f)
		for (g=-1; g <= 1; ++g)
#ifndef TWOD
		  for (h=-1; h <= 1; ++h)
#endif
		    /* Avoid selecting pairs of cells twice */ 
#ifdef TWOD
                    if ((h+3*g)>=(n+3*m))
#else
		    if ((h+3*g+9*f)>=(n+3*m+9*l))
#endif

            {
#ifdef TWOD
	      p = PTR_2D_V(cell_array,i,j  ,cell_dim);
#else
	      p = PTR_3D_V(cell_array,i,j,k,cell_dim);
#endif
	      /* Calculate Indicies of Neighbour */
              r = i+l;  pbc.x = 0;
              s = j+m;  pbc.y = 0;
#ifndef TWOD
              t = k+n;  pbc.z = 0;
#endif

             /* deal with periodic boundary conditions if necessary */
              if (r<0) {
                if (pbc_dirs.x==1) {
                  r = cell_dim.x-1; 
                  pbc.x -= box_x.x;      
                  pbc.y -= box_x.y;
#ifndef TWOD
                  pbc.z -= box_x.z;
#endif
                } else continue;
              }
              if (s<0) {
                if (pbc_dirs.y==1) {
                  s = cell_dim.y-1;
                  pbc.x -= box_y.x;      
                  pbc.y -= box_y.y;
#ifndef TWOD
                  pbc.z -= box_y.z;
#endif
                } else continue;
              }
#ifndef TWOD
              if (t<0) {
                if (pbc_dirs.z==1) {
                  t = cell_dim.z-1;
                  pbc.x -= box_z.x;      
                  pbc.y -= box_z.y;
                  pbc.z -= box_z.z;
                } else continue;
              }
#endif
              if (r>cell_dim.x-1) {
                if (pbc_dirs.x==1) {
                  r = 0; 
                  pbc.x += box_x.x;      
                  pbc.y += box_x.y;
#ifndef TWOD
                  pbc.z += box_x.z;
#endif
                } else continue;
              }
              if (s>cell_dim.y-1) {
                if (pbc_dirs.y==1) {
                  s = 0; 
                  pbc.x += box_y.x;      
                  pbc.y += box_y.y;
#ifndef TWOD
                  pbc.z += box_y.z;
#endif
                } else continue;
              }
#ifndef TWOD
              if (t>cell_dim.z-1) {
                if (pbc_dirs.z==1) {
                  t = 0; 
                  pbc.x += box_z.x;      
                  pbc.y += box_z.y;
                  pbc.z += box_z.z;
                } else continue;
              }
#endif

              /* Neighbour cell (note that p==q ist possible) */
#ifdef TWOD
	      q = PTR_2D_V(cell_array,r,s,cell_dim);
#else
	      q = PTR_3D_V(cell_array,r,s,t,cell_dim);
#endif

	      /* Calculate Indices of second neighbour */
	      u = i+f; ppbc.x = 0;
	      v = j+g; ppbc.y = 0;
#ifndef TWOD
	      w = k+h; ppbc.z = 0;
#endif

	      /* deal with periodic boundary conditions if necessary */
	    if (u<0) {
	      if (pbc_dirs.x==1) {
		u = cell_dim.x-1; 
		ppbc.x -= box_x.x;      
		ppbc.y -= box_x.y;
#ifndef TWOD
		ppbc.z -= box_x.z;
#endif
	      } else continue;
	    }
	      if (v<0) {
		if (pbc_dirs.y==1) {
		  v = cell_dim.y-1;
		  ppbc.x -= box_y.x;      
		  ppbc.y -= box_y.y;
#ifndef TWOD
		  ppbc.z -= box_y.z;
#endif
		} else continue;
	      }
#ifndef TWOD
	      if (w<0) {
		if (pbc_dirs.z==1) {
		  w = cell_dim.z-1;
		  ppbc.x -= box_z.x;      
		  ppbc.y -= box_z.y;
		  ppbc.z -= box_z.z;
		} else continue;
	      }
#endif
	      if (u>cell_dim.x-1) {
		if (pbc_dirs.x==1) {
		  u = 0; 
		  ppbc.x += box_x.x;      
		  ppbc.y += box_x.y;
#ifndef TWOD
		  ppbc.z += box_x.z;
#endif
		 } else continue; 
	      }
	      if (v>cell_dim.y-1) {
		if (pbc_dirs.y==1) {
		  v = 0; 
		  ppbc.x += box_y.x;      
		  ppbc.y += box_y.y;
#ifndef TWOD
		  ppbc.z += box_y.z;
#endif
		  } else continue; 
	      }
#ifndef TWOD
	      if (w>cell_dim.z-1) {
		if (pbc_dirs.z==1) {
		  w = 0; 
		  ppbc.x += box_z.x;      
		  ppbc.y += box_z.y;
		  ppbc.z += box_z.z;
		} else continue;
	      } 
#endif

	      /* Second neighbour cell */
#ifdef TWOD 
	      qq = PTR_2D_V(cell_array, u, v, cell_dim);
#else
	      qq = PTR_3D_V(cell_array, u, v, w, cell_dim);
#endif

	      /* Do the work */
	      do_angle(p,q,qq,pbc,ppbc);
      }
}
