
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
* imd_strain -- calculate strain tensor
*
* uses cell division routines of imd_pair
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#ifndef STRAIN
#define STRAIN
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
  printf("%s [-r<nnn>] [-e<nnn>] [-c<nnn>] [-v] [-p paramter-file]\n",progname); 
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

  /* Initialize cell data structures */
  init_cells();

  /* Read coordinates and displacement */
  read_displacement(infilename);

  /* Calculate strain tensor */
  calc_strain();

  /* Output strain tensor */
  write_data();

  return 0;

}

/******************************************************************************
*
*  read_displacement - reads atom positions and displacements into the cell-array
*
*  The file format is flat ascii, one atom per line, lines beginning
*  with '#' denote comments. Each line consists of
*
*  x y [z] dx dy [dz] [rest]
*
*  where
*
*  x,y,z    are the atom's coordinates
*  dx,dy,dz are the atom's displacement vector
*  rest     is ignored until end of line
*
******************************************************************************/

void read_displacement(str255 infilename)
{
  FILE *infile;
  char buf[512];
  int p;
  vektor pos;
  vektor u;
  cell *to;
  ivektor cellc;

  infile = fopen(infilename,"r");
  if (NULL==infile) {
    sprintf(error_msg,"Cannot open atoms file %s",infilename);
    error(error_msg);
  }

  natoms=0;
  
  /* Read the input file line by line */
  while (!feof(infile)) {

    buf[0] = (char) NULL;
    fgets(buf,sizeof(buf),infile);
    while ('#'==buf[1]) fgets(buf,sizeof(buf),infile); /* eat comments */

#ifdef TWOD
    p = sscanf(buf,"%lf %lf %lf %lf",&pos.x,&pos.y,&u.x,&u.y);
#else
    p = sscanf(buf,"%lf %lf %lf %lf %lf %lf",
               &pos.x,&pos.y,&pos.z,&u.x,&u.y,&u.z);
#endif

    if (p>0) {
      /* compute target cell */
      cellc = cell_coord(pos);
      to = PTR_VV(cell_array,cellc,cell_dim);
      /* enlarge it if necessary */
      if (to->n >= to->n_max) alloc_cell(to,to->n_max+CSTEP);
      /* put the data */
      to->ort[to->n] = pos;
      to->dsp[to->n] = u;
      to->n++;
      natoms++;
    }
  }
  fclose(infile);  

}


/******************************************************************************
*
*  write_data writes strain to *.strain file
*
******************************************************************************/

void write_data(void)
{
  FILE *out;
  str255 fname;
  int i,j,k,l;
  int number=0;
  cell *p; 
 
  sprintf(fname,"%s.strain",infilename);
  
  out = fopen(fname,"w");
  if (NULL == out) error("Cannot open strain file.");

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
	for (l=0;l<p->n; ++l)

	  if (p->empty[l]==0)
	    {
#ifdef TWOD
	      fprintf(out, "%d %f %f %f %f %f\n", ++number, 
		      p->ort[l].x, p->ort[l].y, 
		      p->strain[l].x, p->strain[l].y, p->strain_offdia[l].x);
#else
	      fprintf(out, "%d %f %f %f %f %f %f %f %f %f\n", ++number, 
		      p->ort[l].x, p->ort[l].y, p->ort[l].z, 
		      p->strain[l].x, p->strain[l].y, p->strain[l].z, 
		      p->strain_offdia[l].x, p->strain_offdia[l].y, 
		      p->strain_offdia[l].z);
#endif
	    }
      }

  fclose(out);
}


/******************************************************************************
*
*  calc_strain calculates the strain tensor for all atoms
*
******************************************************************************/

void calc_strain(void)
{
  cell *p,*q;
  int i,j,k;
  int l,m,n;
  int r,s,t;
  int u,v,w;
  int g,h;
  int number;                          /* number of neighbours */ 
  int totalnumber = 0;                 /* number of neighbours of all atoms */
  int emptynumber = 0;                 /* number of atoms with less than 3 neighbours */ 
  int maxnumber = 0, minnumber = 1000; /* maximal and minimal number of neighbours occuring */
  int num = 4;
  vektor pbc;
  vektor *d, *du, *tmp;
  real a[3][3], b[3][3];
  real radius;
  real det;

  /* Allocate memory for temporary variables */
  d   = (vektor *) malloc( num * sizeof(vektor));
  du  = (vektor *) malloc( num * sizeof(vektor));
  tmp = (vektor *) malloc( num * sizeof(vektor));


  /* Initialization */

  /* for each cell */
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
	  /* For each atom in this cell */
	  for (u=0; u < p->n; ++u)
	    { 
	      p->strain[u].x        = 0.0;
	      p->strain[u].y        = 0.0;
#ifndef TWOD
	      p->strain[u].z        = 0.0;
#endif
	      p->strain_offdia[u].x = 0.0;
#ifndef TWOD
	      p->strain_offdia[u].y = 0.0;
	      p->strain_offdia[u].z = 0.0;
#endif
	      p->empty[u]           = 0;
	    }
	}

  /* Compute strain tensor */

  /* For each cell */
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

	  /* For each atom in this first cell */
	  for (u=0; u<p->n; ++u)
	    {

	      number = 0;
	      
	      /* For the neighbours of the first cell */
	      for (l=-1; l <= 1; ++l)
		for (m=-1; m <= 1; ++m)
#ifndef TWOD
		  for (n=-1; n <= 1; ++n)
#endif
		    {
	      /* Calculate Indices of Neighbour */
	      r = i+l;  pbc.x = 0;
	      s = j+m;  pbc.y = 0;
#ifndef TWOD
	      t = k+n;  pbc.z = 0;
#endif
              /* Deal with periodic boundary conditions if necessary */
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

	      /* For each atom in the second cell */

	      for( v=0; v<q->n; ++v)
		{
		  /* Check whether there is enough memory for tmp variables */
		  if( number >= num) {
			++num;
			d    = (vektor *) realloc( d,   num * sizeof(vektor)); 
			du   = (vektor *) realloc( du,  num * sizeof(vektor));
			tmp  = (vektor *) realloc( tmp, num * sizeof(vektor)); 
		      }
		  
		  /* Calculate distance */
		  d[number].x = q->ort[v].x - q->dsp[v].x - p->ort[u].x + p->dsp[u].x + pbc.x;
		  d[number].y = q->ort[v].y - q->dsp[v].y - p->ort[u].y + p->dsp[u].y + pbc.y;
#ifndef TWOD
		  d[number].z = q->ort[v].z - q->dsp[v].z - p->ort[u].z + p->dsp[u].z + pbc.z;
#endif

		  radius = sqrt( (double)(SPROD(d[number],d[number])) );

		  if ( radius > 0.01 && radius < r_max ) 
		    {

		      /* Calculate differences of displacements */
		      du[number].x = q->dsp[v].x - p->dsp[u].x;
		      du[number].y = q->dsp[v].y - p->dsp[u].y;
#ifndef TWOD
		      du[number].z = q->dsp[v].z - p->dsp[u].z;
#endif

		      ++number; /* Number of neighbours*/

		    }
		} /* v */
	    } /* Neighbours of p */
	      
	      /* Calculate transformation matrix of du */

	      /* a = dx * dxT */

	      totalnumber += number;
	      maxnumber    = (number>maxnumber) ? number : maxnumber;
	      minnumber    = (number<minnumber) ? number : minnumber;

	      /* At least 3 neighbour atoms are required */
	      if (number < 3) { 
		p->empty[u] = 1;
		++emptynumber;
	      }
	      else
#ifdef TWOD
		{
	  for(g=0; g<2; ++g)
	    for(h=0; h<2; ++h)
	       a[g][h] = 0.0;

	  for (w=0; w<number; ++w) 
	    {
	      a[0][0] += d[w].x * d[w].x;  a[0][1] += d[w].x * d[w].y; 
	      a[1][0] += d[w].y * d[w].x;  a[1][1] += d[w].y * d[w].y;
	    }

	  /* b = Inverse of a */

	  det = a[0][0] * a[1][1] - a[0][1] * a[1][0];

	  if (det == 0.0) error("Transformation matrix zero.");

	  b[0][0] =  a[1][1] / det;
	  b[0][1] = -a[0][1] / det;
	  b[1][0] = -a[1][0] / det;
	  b[1][1] =  a[0][0] / det;

	  /* tmp = dx * b  */

	  for (w=0; w<number; ++w)
	    {
	      tmp[w].x = d[w].x * b[0][0] + d[w].y * b[1][0];
	      tmp[w].y = d[w].x * b[0][1] + d[w].y * b[1][1]; 
	    }

	  /* strain = (symmetrized) du * tmp  */

	  for (w=0; w<number; ++w)
	    {
	      p->strain[u].x        += du[w].x * tmp[w].x;
	      p->strain[u].y        += du[w].y * tmp[w].y;
	      
	      p->strain_offdia[u].x += ( du[w].y * tmp[w].x + du[w].x * tmp[w].y ) / 2;
	    }
#else /* 3D */
	{

	  for (g=0; g<3; ++g)
	    for (h=0; h<3; ++h)
	      a[g][h] = 0.0;

	  for (w=0; w<number; ++w) 
	    {
	      a[0][0] += d[w].x * d[w].x;  a[0][1] += d[w].x * d[w].y;  a[0][2] += d[w].x * d[w].z;
	      a[1][0] += d[w].y * d[w].x;  a[1][1] += d[w].y * d[w].y;  a[1][2] += d[w].y * d[w].z;
	      a[2][0] += d[w].z * d[w].x;  a[2][1] += d[w].z * d[w].y;  a[2][2] += d[w].z * d[w].z;
	    }

	  /* b = Inverse of a */

	  det = a[0][0] * a[1][1] * a[2][2] +
	        a[0][1] * a[1][2] * a[2][0] +
	        a[0][2] * a[1][0] * a[2][1] -
	        a[2][0] * a[1][1] * a[0][2] -
	        a[2][1] * a[1][2] * a[0][0] -
	        a[2][2] * a[1][0] * a[0][1];

	  if (det == 0.0) error("Transformation matrix singular.");

	  b[0][0] = ( a[1][1] * a[2][2] - a[1][2] * a[2][1] ) / det;
	  b[0][1] = ( a[2][1] * a[0][2] - a[0][1] * a[2][2] ) / det;
	  b[0][2] = ( a[0][1] * a[1][2] - a[1][1] * a[0][2] ) / det;

	  b[1][0] = ( a[1][2] * a[2][0] - a[1][0] * a[2][2] ) / det;
	  b[1][1] = ( a[0][0] * a[2][2] - a[2][0] * a[0][2] ) / det;
	  b[1][2] = ( a[0][2] * a[1][0] - a[0][0] * a[1][2] ) / det;
	  
	  b[2][0] = ( a[1][0] * a[2][1] - a[1][1] * a[2][0] ) / det;
	  b[2][1] = ( a[0][1] * a[2][0] - a[0][0] * a[2][1] ) / det;
	  b[2][2] = ( a[0][0] * a[1][1] - a[0][1] * a[1][0] ) / det;

	  /* tmp = dx * b  */

	  for (w=0; w<number; ++w)
	    {
	      tmp[w].x = d[w].x * b[0][0] + d[w].y * b[1][0] + d[w].z * b[2][0];
	      tmp[w].y = d[w].x * b[0][1] + d[w].y * b[1][1] + d[w].z * b[2][1];
	      tmp[w].z = d[w].x * b[0][2] + d[w].y * b[1][2] + d[w].z * b[2][2];
	    }

	  /* strain = (symmetrized) du * tmp  */

	  for (w=0; w<number; ++w)
	    {
	      p->strain[u].x += du[w].x * tmp[w].x;
	      p->strain[u].y += du[w].y * tmp[w].y;
	      p->strain[u].z += du[w].z * tmp[w].z;

	      p->strain_offdia[u].x += ( du[w].y * tmp[w].z + du[w].z * tmp[w].y ) / 2;
	      p->strain_offdia[u].y += ( du[w].x * tmp[w].z + du[w].z * tmp[w].x ) / 2;
	      p->strain_offdia[u].z += ( du[w].x * tmp[w].y + du[w].y * tmp[w].x ) / 2;
	    }
#endif	 /* Not TWOD */
	  
	} /* number > 2 */

		} /* u */

	    } /* First cell */

	  /* Statistics */
	  printf("Maximal number of neighbour atoms: %d\n", maxnumber);
	  printf("Minimal number of neighbour atoms: %d\n", minnumber);
	  printf("Average number of neighbour atoms: %f\n", (float) totalnumber/natoms );
	  if(emptynumber>0)
	    printf("Number of omitted atoms: %d (%.2f %%)\n", emptynumber, (float) emptynumber/natoms );

}









