
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
* imd_coord -- calculate coordination numbers
*
* A descendant of imd_pair
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#define COORD
#include "util.h"

/******************************************************************************
*
*  Usage -- educate users
*
*  Compilation: gcc -O [-DTWOD] [-DSLOTS=<nnn>] [-DSINGLE] imd_coord.c -lm 
*
******************************************************************************/

void usage(void)
{ 
  printf("%s [-r<nnn>] [-e<nnn>] [-p paramter-file]\n",progname); 
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

  tablesize = ntypes*ntypes*sizeof(real);
  numbers = (real *) malloc(tablesize);
  if (NULL==numbers) error("Cannot allocate memory for numbers.");
  num_dim.x = ntypes;
  num_dim.y = ntypes;

  for (i=0; i<ntypes; ++i)
    for (j=0; j<ntypes; ++j)
	*PTR_2D_V(numbers,i,j,num_dim) = 0.0;

  r2_cut = SQR(r_max);

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
*  write_data writes numbers to *.coord file
*
******************************************************************************/

void write_data()
{
  FILE *out;
  str255 fname;
  int i,j;

  if (-1==restart)
    sprintf(fname,"%s.coord",infilename);
  else
    sprintf(fname,"%s.%u.coord",outfilename,restart);

  out = fopen(fname,"w");
  if (NULL == out) error("Cannot open outfile.");

  for (i=0; i<ntypes; ++i)
    for (j=i; j<ntypes; ++j)
      fprintf(out,"%d%d %f\n",i,j,*PTR_2D_V(numbers,i,j,num_dim)/natoms);

  fprintf(out,"\n");

  fclose(out);

  printf("\n");
  printf("---------------------\n");
  printf(" Bond   Coordination\n\n");

  for (i=0; i<ntypes; ++i)
    for (j=i; j<ntypes; ++j)
      printf("  %d%d     %f\n",i,j,*PTR_2D_V(numbers,i,j,num_dim)/natoms);

  printf("\n");
  printf("---------------------\n");

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

      if (q_typ > p_typ) {temp = p_typ; p_typ = q_typ; q_typ = temp;}

      if ( radius <= r_max )
	*PTR_2D_V(numbers, q_typ, p_typ, num_dim) += 2;
    }
}



