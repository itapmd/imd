
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
* imd_pair -- calculate pair distibution functions
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#define PAIR
#include "util.h"

/******************************************************************************
*
*  Usage -- educate users
*
*  Compilation: gcc -O [-DTWOD] [-DSLOTS=<nnn>] [-DSINGLE] imd_pair.c -lm
*
******************************************************************************/

void usage(void)
{ 
  printf("%s [-r<nnn>] [-a<nnn>] [-e<nnn>] [-p paramter-file]\n",progname); 
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
  int i,j,k;

  /* Read Parameters from parameter file */
  read_parameters(argc,argv);

  tablesize = SLOTS*ntypes*ntypes*sizeof(real);
  histogram = (real *) malloc(tablesize);
  if (NULL==histogram) error("Cannot allocate memory for histograms.");
  hist_dim.x = SLOTS;
  hist_dim.y = ntypes;
  hist_dim.z = ntypes;

  for (i=0; i<SLOTS; ++i)
    for (j=0; j<ntypes; ++j)
      for (k=0; k<ntypes; ++k)
	*PTR_3D_V(histogram,i,j,k,hist_dim) = 0.0;

  r2_cut = SQR(r_max);

  /* Initialize cell data structures */
  init_cells();

  /* Read atoms */
  read_atoms(infilename);

  /* Calculate the distances */
  do_work(do_cell_pair);

  /* Output results */
  write_data();

  return 0;

}


/******************************************************************************
*
*  write_data writes histogram to *.pair file
*
******************************************************************************/

void write_data()
{
  FILE *out;
  str255 fname;
  int i,j,k;
  real r;
  real f;

  if (0==restart)
    sprintf(fname,"%s.pair",infilename);
  else
    sprintf(fname,"%s.%u.pair",outfilename,restart);

  out = fopen(fname,"w");
  if (NULL == out) error("Cannot open histograms file.");

  for (i=1; i<SLOTS; ++i) {
    r = ((float) i / SLOTS * (r_max - r_min)) + r_min;
    fprintf(out,"%f ", r);
    f = natoms;
    for (j=0; j<ntypes; ++j)
      for (k=j; k<ntypes; ++k)
	fprintf(out,"%f ",*PTR_3D_V(histogram,i,j,k,hist_dim)/f);
    f = natoms * 4 * 3.14159265 * SQR(r);
    for (j=0; j<ntypes; ++j)
      for (k=j; k<ntypes; ++k)
	fprintf(out,"%f ",*PTR_3D_V(histogram,i,j,k,hist_dim)/f);
    fprintf(out,"\n");
  }
  fclose(out);
}


/******************************************************************************
*
*  do_cell_pair calulates the distances for atoms in two cells
*
******************************************************************************/

void do_cell_pair(cell *p, cell *q, vektor pbc)
{
  int i,j,k;
  int temp;
  vektor d;
  real radius;
  int p_typ,q_typ;

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

      k     = (int) ( SLOTS * (radius - r_min) / (r_max - r_min));
      p_typ = p->sorte[i];
      q_typ = q->sorte[j];

      if (q_typ > p_typ) {temp = p_typ; p_typ = q_typ; q_typ = temp;}

      if ((k>0) && (k<SLOTS))
	*PTR_3D_V(histogram, k , q_typ, p_typ, hist_dim) += 2;
    }
}




