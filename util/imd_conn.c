
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

#define CONN
#include "util.h"

/******************************************************************************
*
*  Usage -- educate users
*
*  Compilation: gcc -O [-DTWOD] [-DSINGLE] imd_conn.c -lm 
*
******************************************************************************/

void usage(void)
{ 
  printf("%s [-r<nnn>] [-p paramter-file]\n",progname); 
  exit(1); 
}


/*****************************************************************************
*
*  main
*
*****************************************************************************/

int main(int argc, char **argv)
{
  int i,j;

  /* Read Parameters from parameter file */
  read_parameters(argc,argv);

  /* Calculate cutoff radius */
  r2_cut = 0.0;
  for (i=0; i<ntypes; i++)
    for (j=0; j<ntypes; j++) r2_cut = MAX( r2_cut, SQR(r_max2d[i][j]));

  /* Initialize cell data structures */
  init_cells();

  /* Read atoms */
  read_atoms(infilename);

  /* Allocate connection matrix */
  nn  = (int *)    calloc(natoms, sizeof(int));
  tp  = (int *)    malloc(natoms * sizeof(int));
  cm  = (int *)    malloc(natoms * MAXNEIGH * sizeof(int));
  pos = (vektor *) malloc(natoms * sizeof(vektor));
  if ((nn==NULL) || (tp==NULL) | (cm==NULL) | (pos==NULL))
    error("cannot allocate connection matrix");

  /* Calculate connection matrix */
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

  if (0==restart)
    sprintf(fname,"%s.conn",infilename);
  else
    sprintf(fname,"%s.%u.conn",outfilename,restart);
  out = fopen(fname,"w");
  if (NULL == out) error("Cannot open connection matrix file.");

  for (i=0; i<natoms; ++i) {
#ifdef TWOD
    fprintf(out,"%d %d %f %f %d",   i,tp[i],pos[i].x,pos[i].y,         nn[i]);
#else
    fprintf(out,"%d %d %f %f %f %d",i,tp[i],pos[i].x,pos[i].y,pos[i].z,nn[i]);
#endif
    for (j=0; j<nn[i]; j++) fprintf(out," %d",*PTR_2D_V(cm,i,j,cm_dim));
    fprintf(out,"\n");
  }
  fclose(out);
}


/******************************************************************************
*
*  calulate the distances and connection matrix entries for atoms in two cells
*
******************************************************************************/

void do_cell_pair(cell *p, cell *q, vektor pbc)
{
  int i,j,ii,jj;
  vektor d;
  real radius;
  int p_typ,q_typ;

  /* For each atom in first cell */
  for (i=0; i<p->n; ++i) {

    ii      = p->nummer[i];
    tp[ii]  = p->sorte[i];
    pos[ii] = p->ort[i];
    
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
      if (radius < r_max2d[p_typ][q_typ]) {
        ii = p->nummer[i];
        jj = q->nummer[j];
        *PTR_2D_V(cm,ii,nn[ii],cm_dim) = jj;  nn[ii]++;
        *PTR_2D_V(cm,jj,nn[jj],cm_dim) = ii;  nn[jj]++;
      }
    }
  }
}
