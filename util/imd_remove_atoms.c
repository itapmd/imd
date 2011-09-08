
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
* imd_remove_atoms  -- iteratively removes too close atoms
*
******************************************************************************/


#ifndef REMAT
#define REMAT
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
  printf("%s [-a<nnn>] [-p paramter-file]\n",progname); 
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
  int n_nnn_crit;
  int round;

  ndeletedatoms=0;

  /* Read command line arguments */
  read_command_line(argc,argv);

  /* Read Parameters from parameter file */
  read_parameters();
 
  r2_cut = SQR(r_min);
  printf("Using r_crit=%lf (good estimate: 0.45 * lattice parameter)\n",r_crit);
  if(r_crit<=0)
    exit(20);

  /* read box from file header */
  if (box_from_header) read_box(infilename);

  /* Initialize cell data structures */
  init_cells();

  /* Read atoms */
  printf("reading atoms\n"); fflush(stdout);
  read_atoms(infilename);
  printf("read %d atoms\n",natoms); fflush(stdout);

  n_nnn_crit=6;
  round =0;
  /* Calculate the distances */
  do 
    {
      do_work(do_cell_pair);
      printf("determined too close neighbors, round %d n_nnn_crit %d\n",round,n_nnn_crit); fflush(stdout);
      

  
      calc_n_toclose_nn(n_nnn_crit);

      for(j=0;j<12;j++)
	{
	  printf("atoms with %d too close neighbors: %d\n",j+1,n_nnn[j]);
	}
      n_nnn_crit--;
      round++;
      printf("deleted %d atoms in total\n\n",ndeletedatoms);
    }while (n_nnn[1]>0);

  do_work(do_cell_pair);
  printf("final determined too close neighbors, round %d n_nnn_crit %d\n",round,n_nnn_crit); fflush(stdout);
  calc_n_toclose_nn(-1);
  printf("deleted %d atoms in total \n\n",ndeletedatoms);
  for(j=0;j<12;j++)
    {
      printf("atoms with %d too close neighbors: %d\n",j+1,n_nnn[j]);
    }

  write_atoms();
  return 0;

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
    //for (j =  0; j < q->n; ++j) {

      /* Calculate distance */
      d.x =q->ort[j].x - p->ort[i].x + pbc.x;
      d.y =q->ort[j].y - p->ort[i].y + pbc.y;
#ifndef TWOD
      d.z =q->ort[j].z - p->ort[i].z + pbc.z;
#endif

      radius = sqrt( (double)(SPROD(d,d)) );
      if(radius<=r_crit)
	{
	  if(q->flag[j]==0 && p->flag[i]==0)
	    p->nnn[i]++;
	}
     
    }
}

void calc_n_toclose_nn(int n_nnn_crit)
{
  cell *p,*q;
  int i,j,k;
  int l,m,n;
  int r,s,t;
  vektor pbc;
  
  ndeletedatoms=0;
  for(j=0;j<12;j++)
    n_nnn[j]=0;

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

	 /* For each atom in cell */
	for (l = 0;l < p->n; ++l) 
	  {
	    for(m=0;m<12;m++)
	      {
		if(p->nnn[l]==m+1)
		  {
		    n_nnn[m]++;
		  }
	      }
	    if(n_nnn_crit >0)
	      {
		if(n_nnn_crit >1)
		  {
		    if(p->nnn[l]>=n_nnn_crit)
		      {
			p->flag[l]=1;
		      }
		  }
		else if(n_nnn_crit ==1)
		  {
		    if(p->nnn[l]>=n_nnn_crit && drand48()<0.5)
		      {
			p->flag[l]=1;
		      }
		  }
	      }
	    p->nnn[l]=0;
	    if(p->flag[l]==1)
	      ndeletedatoms++;
	  }


      }
}


void write_atoms()
{
  int i, j, k, l;
  cell *p;
  FILE *out, *remout;
  str255 ofname,rfname;

  
  sprintf(ofname,"%s.cleaned",infilename);
  sprintf(rfname,"%s.removed",infilename);
  
  out = fopen(ofname,"w");
  if (NULL == out) {printf("Cannot open %s for writing.",ofname);exit(20);}
  remout = fopen(rfname,"w");
  if (NULL == out) {printf("Cannot open %s for writing.",rfname);exit(20);}


  fprintf(out, "#F A 1 1 1 3 0 0\n");
  fprintf(out, "#C number type mass x y z\n");
  fprintf(out, "#X %f %f %f\n", box_x.x, box_x.y, box_x.z);
  fprintf(out, "#Y %f %f %f\n", box_y.x, box_y.y, box_y.z);
  fprintf(out, "#Z %f %f %f\n", box_z.x, box_z.y, box_z.z);
  fprintf(out, "#E\n");
  fprintf(remout, "#F A 1 1 1 3 0 0\n");
  fprintf(remout, "#C number type mass x y z\n");
  fprintf(remout, "#X %f %f %f\n", box_x.x, box_x.y, box_x.z);
  fprintf(remout, "#Y %f %f %f\n", box_y.x, box_y.y, box_y.z);
  fprintf(remout, "#Z %f %f %f\n", box_z.x, box_z.y, box_z.z);
  fprintf(remout, "#E\n");


  /* For each cell */
  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
      for (k=0; k < cell_dim.z; ++k)
      {
        p = PTR_3D_V(cell_array,i,j,k,cell_dim);
	/* For each atom in cell */
	for (l=0; l<p->n; l++)
	  if (p->flag[l] == 0 )
	    fprintf(out, "%d %d %f %f %f %f\n", 
		    p->nummer[l], p->sorte[l], p->masse[l], 
		    p->ort[l].x, p->ort[l].y, p->ort[l].z );
	  else
	    fprintf(remout, "%d %d %f %f %f %f\n", 
		    p->nummer[l], p->sorte[l], p->masse[l], 
		    p->ort[l].x, p->ort[l].y, p->ort[l].z );
      }
  fclose(out);
  fclose(remout);
}
