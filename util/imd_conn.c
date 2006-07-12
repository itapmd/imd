
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
*  imd_conn -- calculate connection matrix
*
*  A descendent of imd_pair
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#ifndef CONN
#define CONN
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
  printf("%s [-r<nnn>] [-A<nnn>] [-v] [-p paramter-file]\n",progname); 
  exit(1); 
}


/*****************************************************************************
*
*  main
*
*****************************************************************************/

int main(int argc, char **argv)
{
  int i;

  /* Read command line arguments */
  read_command_line(argc,argv);

  /* Read Parameters from parameter file */
  read_parameters();

  /* Calculate cutoff radius */
  r2_cut = 0.0;
  for (i=0; i<SQR(ntypes); i++) r2_cut = MAX( r2_cut, SQR(r_cut[i]));

  /* read box from file header */
  if (box_from_header) read_box(infilename);

  /* Initialize cell data structures */
  init_cells();

  /* Read atoms */
  read_atoms(infilename);

  /* Allocate connection matrix */
  nn  = (int *)    calloc(natoms, sizeof(int));
  ind = (int *)    malloc((n_max-n_min+1) * sizeof(int));
  num = (int *)    malloc(natoms * sizeof(int));
  tp  = (int *)    malloc(natoms * sizeof(int));
  cm  = (int *)    malloc(natoms * maxneigh * sizeof(int));
  pos = (vektor *) malloc(natoms * sizeof(vektor));
  if ((nn==NULL) || (tp==NULL) | (cm==NULL) | (pos==NULL))
    error("cannot allocate connection matrix");

  /* Calculate connection matrix */
  make_numbers();
  do_work(do_cell_pair);

  /* Write connection matrix */
  write_data();

  return 0;

}


/******************************************************************************
*
*  write connection matrix
*
******************************************************************************/

void write_data()
{
  FILE *out;
  str255 fname;
  int i,j;

  if (-1==restart)
    sprintf(fname,"%s.conn",infilename);
  else
    sprintf(fname,"%s.%u.conn",outfilename,restart);
  out = fopen(fname,"w");
  if (NULL == out) error("Cannot open connection matrix file.");

  for (i=0; i<natoms; ++i) {
#ifdef TWOD
    fprintf(out,"%d %d %f %f %d",   
            num[i],tp[i],pos[i].x,pos[i].y,         nn[i]);
#else
    fprintf(out,"%d %d %f %f %f %d",
            num[i],tp[i],pos[i].x,pos[i].y,pos[i].z,nn[i]);
#endif
    for (j=0; j<nn[i]; j++) fprintf(out," %d",*PTR_2D_V(cm,i,j,cm_dim));
    fprintf(out,"\n");
  }
  fclose(out);
}



/******************************************************************************
*
*  translation from number to index, and vice versa
*
******************************************************************************/

void make_numbers(void)
{
  cell *p;
  int i,j=0,k,l,m;

  /* for each cell */
  for (k=0; k < cell_dim.x; ++k)
    for (l=0; l < cell_dim.y; ++l)
#ifndef TWOD
      for (m=0; m < cell_dim.z; ++m)
#endif
      {
#ifdef TWOD
        p = PTR_2D_V(cell_array,k,l  ,cell_dim);
#else
        p = PTR_3D_V(cell_array,k,l,m,cell_dim);
#endif
        for (i=0; i<p->n; ++i) {
          tp [j] = p->sorte [i];
          pos[j] = p->ort   [i];
          ind[p->nummer[i]-n_min] = j;  /* index  as function of number */
          num[j] = p->nummer[i];        /* number as function of index  */
          j++;
	}
      }

}

/******************************************************************************
*
*  calulate the distances and connection matrix entries for atoms in two cells
*
******************************************************************************/

void do_cell_pair(cell *p, cell *q, vektor pbc)
{
  int i,j,k;
  vektor d;
  real radius;

  /* For each atom in first cell */
  for (i=0; i<p->n; ++i) {

    /* For each atom in neighbouring cell */
    /* If p==q, use only rest of atoms */
    for (j = ((p==q) ? i+1 : 0); j < q->n; ++j) {
      
      /* Calculate distance */
      d.x = q->ort[j].x - p->ort[i].x + pbc.x;
      d.y = q->ort[j].y - p->ort[i].y + pbc.y;
#ifndef TWOD
      d.z = q->ort[j].z - p->ort[i].z + pbc.z;
#endif

      radius = sqrt( (double)(SPROD(d,d)) );
      if (radius < r_cut[ (p->sorte[i])*ntypes + q->sorte[j] ]) {

        /* update neighbor table of p-particle */
        k = ind[p->nummer[i]-n_min];
	if (nn[k]<maxneigh) {
          *PTR_2D_V(cm,k,nn[k],cm_dim) = q->nummer[j];  nn[k]++;
        } else error("maximum number of neighbors exceeded");

        /* update neighbor table of q-particle */
        k = ind[q->nummer[j]-n_min];
        if (nn[k]<maxneigh) {
          *PTR_2D_V(cm,k,nn[k],cm_dim) = p->nummer[i];  nn[k]++;
	} else error("maximum number of neighbors exceeded");

      }
    }
  }
}
