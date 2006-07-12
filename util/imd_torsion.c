
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
* imd_torsion -- calculate torsion distribution functions
*
* A descendant of imd_pair
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#ifndef TORSION
#define TORSION
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
  printf("%s [-r<nnn>] [-A<nnn>] [-e<nnn>] [-v] [-p paramter-file]\n",progname); 
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
  int n,i,j,k,l;

  /* Read command line arguments */
  read_command_line(argc,argv);

  /* Read Parameters from parameter file */
  read_parameters();

  tablesize = slots*ntypes*ntypes*ntypes*ntypes*sizeof(real);
  histogram = (real *) malloc(tablesize);
  if (NULL==histogram) error("Cannot allocate memory for histograms.");
  hist_dim.n = slots;
  hist_dim.i = ntypes;
  hist_dim.j = ntypes;
  hist_dim.k = ntypes;
  hist_dim.l = ntypes;

  for (n=0; n<slots; ++n)
    for (i=0; i<ntypes; ++i)
      for (j=0; j<ntypes; ++j)
	for (k=0; k<ntypes; ++k)
	  for (l=0; l<ntypes; ++l)
	    *PTR_5D_V(histogram,n,i,j,k,l,hist_dim) = 0.0;

  r2_cut = SQR(r_max);

  /* read box from file header */
  if (box_from_header) read_box(infilename);

  /* Initialize cell data structures */
  init_cells();

  /* Read atoms */
  read_atoms(infilename);

  /* Calculate the torsion angles */
  calc_angles();

  /* Output results */
  write_data();

  return 0;

}


/******************************************************************************
*
*  write_data writes histogram to *.torsion file
*
******************************************************************************/

void write_data()
{
  FILE *out;
  str255 fname;
  int n,i,j,k,l;
  real phi;

  if (-1==restart)
    sprintf(fname,"%s.torsion",infilename);
  else
    sprintf(fname,"%s.%u.torsion",outfilename,restart);

  out = fopen(fname,"w");
  if (NULL == out) error("Cannot open histograms file.");

  for (n = 0; n < slots; ++n) {
    phi = ((float) n / slots * 180);
    fprintf(out,"%f ", phi);
    for (i = 0; i < ntypes; ++i)
      for (j = i; j < ntypes; ++j)
	for (k = 0; k < ntypes; ++k)
	  for (l = (i != j ? 0 : k); l < ntypes; ++l)
	    if (nangles != 0)
	      fprintf(out,"%f ",*PTR_5D_V(histogram,n,i,j,k,l,hist_dim)/nangles);
	    else 
	      error("No neighbouring atoms in this distance range.");
    fprintf(out,"\n");
  };

  printf("%d torsion angles computed\n", nangles);

  fclose(out);

}


/******************************************************************************
*
*  do_angle calulates the torsion angles for atoms in four cells
*
******************************************************************************/

void do_angle(cell *p, cell *q, cell *r, cell *s, 
              vektor pbc_q, vektor pbc_r, vektor pbc_s)
{
  int i, j, k, l;
  int ang;
  int temp;
  vektor d_ij, d_ik, d_jk, d_il, d_jl, d_kl;
  vektor v_k, v_l;
  real radius_ij, radius_ik, radius_jk, radius_il, radius_jl, radius_kl;
  real phi;
  real betrag_v_k, betrag_v_l;
  double sprod;
  int p_typ, q_typ, r_typ, s_typ;
  
  /*
                      l
                     /
                    /
             i-----j
            /
           /
          k

	  */

  /* For each atom in cell p */
  for (i = 0;i < p->n; ++i) 
    /* For each atom in neighbouring cell q */
    for (j = (p == q ? i+1 : 0); j < q->n; ++j) 
      /* For each atom in neighbouring cell r of p */
      for (k = 0; k < r->n; ++k) 
	/* for each atom in neighbouring cell s of q */
	for (l = 0; l < s->n; ++l)
	  {

	    /* Calculate distance vectors */
	    d_ij.x = q->ort[j].x - p->ort[i].x + pbc_q.x;
	    d_ij.y = q->ort[j].y - p->ort[i].y + pbc_q.y;
	    d_ij.z = q->ort[j].z - p->ort[i].z + pbc_q.z;

	    d_ik.x = r->ort[k].x - p->ort[i].x + pbc_r.x;
	    d_ik.y = r->ort[k].y - p->ort[i].y + pbc_r.y;
	    d_ik.z = r->ort[k].z - p->ort[i].z + pbc_r.z;

	    d_jl.x = s->ort[l].x - q->ort[j].x + pbc_s.x;
	    d_jl.y = s->ort[l].y - q->ort[j].y + pbc_s.y;
	    d_jl.z = s->ort[l].z - q->ort[j].z + pbc_s.z;

	    d_jk.x = d_ik.x - d_ij.x;
	    d_jk.y = d_ik.y - d_ij.y;
	    d_jk.z = d_ik.z - d_ij.z;

	    d_il.x = d_ij.x + d_jl.x;
	    d_il.y = d_ij.y + d_jl.y;
	    d_il.z = d_ij.z + d_jl.z;

	    d_kl.x = d_jl.x - d_jk.x;
	    d_kl.y = d_jl.y - d_jk.y;
	    d_kl.z = d_jl.z - d_jk.z;

	    radius_ij = sqrt( (double)(SPROD(d_ij,d_ij)) );
	    radius_ik = sqrt( (double)(SPROD(d_ik,d_ik)) );
	    radius_jl = sqrt( (double)(SPROD(d_jl,d_jl)) );
	    radius_jk = sqrt( (double)(SPROD(d_jk,d_jk)) );
	    radius_il = sqrt( (double)(SPROD(d_il,d_il)) );
	    radius_kl = sqrt( (double)(SPROD(d_kl,d_kl)) );

	    /* v_k = d_ik x d_jk */
	    v_k.x = d_ik.y * d_jk.z - d_ik.z * d_jk.y; 
	    v_k.y = d_ik.z * d_jk.x - d_ik.x * d_jk.z;
	    v_k.z = d_ik.x * d_jk.y - d_ik.y * d_jk.x;

	    betrag_v_k = sqrt( (double)(SPROD(v_k,v_k)) );
	    
	    /* v_l = d_il x d_jl */
	    v_l.x = d_il.y * d_jl.z - d_il.z * d_jl.y; 
	    v_l.y = d_il.z * d_jl.x - d_il.x * d_jl.z;
	    v_l.z = d_il.x * d_jl.y - d_il.y * d_jl.x;

	    betrag_v_l = sqrt( (double)(SPROD(v_l,v_l)) );

	    /* Calculate torsion angles */
	    if ( (radius_ij < r_max) && (radius_ik < r_max) 
              && (radius_jl < r_max) 
              && (radius_ij*radius_ik*radius_jl*radius_jk*radius_il*radius_kl > 0.0) && (betrag_v_k > 0.0) && (betrag_v_l > 0.0) ) 
	      {
		++nangles;
		sprod = (double)( SPROD(v_k,v_l)/( betrag_v_k * betrag_v_l ));
		if (sprod < -1.0) sprod = -1.0;
		  
		phi = (double) (acos(sprod));
  
		ang = (int) ( slots * phi / 3.141592654 );

		p_typ = p->sorte[i];
		q_typ = q->sorte[j];
		r_typ = r->sorte[k];
		s_typ = s->sorte[l];


		if (p_typ == q_typ) { 
		  if ( r_typ > s_typ ) {
		    temp = s_typ; s_typ = r_typ; r_typ = temp;
		  }
		}
		else if (p_typ > q_typ) {
		  temp = q_typ; q_typ = p_typ; p_typ = temp;
		  temp = s_typ; s_typ = r_typ; r_typ = temp;
		}

		if ((ang >= 0) && (ang < slots))
		  ++*PTR_5D_V(histogram, ang , p_typ, q_typ, r_typ, s_typ, hist_dim);
	      } 
	  }

}


/******************************************************************************
*
*  calc_angles calulates the torsion angles for all atoms
*
******************************************************************************/

void calc_angles(void)
{
  cell *p, *q,* r, *s;
  int i1, i2, i3, j1, j2, j3;
  int k1, k2, k3, l1, l2, l3;
  int a1, a2, a3, b1, b2, b3, c1, c2, c3;
  vektor pbc_q, pbc_r, pbc_s;

  /* for each cell p */
  for (i1 = 0; i1 < cell_dim.x; ++i1)
    for (i2 = 0; i2 < cell_dim.y; ++i2)
      for (i3 = 0; i3 < cell_dim.z; ++i3)

	/* For half of the neighbours q of this cell */
	for (j1 = 0; j1 <= 1; ++j1)
	  for (j2 = -j1; j2 <= 1; ++j2)
	    for (j3= (j1 == 0 ? -j2 : -j1); j3 <= 1; ++j3)

	      /* For the neighbours r of cell p */
	      for (k1 = -1; k1 <= 1; ++k1)
		for (k2 = -1; k2 <= 1; ++k2)
		  for (k3 = -1; k3 <= 1; ++k3)

		    /* for the neighbours s of cell q */
		    for (l1 = -1; l1 <= 1; ++l1)
		      for (l2 = -1; l2 <= 1; ++l2)
			for (l3 = -1; l3 <= 1; ++l3)

			  {
			    /* Given cell */

			    p = PTR_3D_V(cell_array,i1,i2,i3,cell_dim);

			    /* Calculate Indices of Neighbour */
			    a1 = i1 + j1;  pbc_q.x = 0;
			    a2 = i2 + j2;  pbc_q.y = 0;
			    a3 = i3 + j3;  pbc_q.z = 0;

			    /* deal with periodic boundary conditions if necessary */
			    if (a1 < 0) {
			      if (pbc_dirs.x == 1) {
				a1 = cell_dim.x-1; 
				pbc_q.x -= box_x.x;      
				pbc_q.y -= box_x.y;
				pbc_q.z -= box_x.z;

			      } else continue;
			    }
			    if (a2 < 0) {
			      if (pbc_dirs.y == 1) {
				a2 = cell_dim.y-1;
				pbc_q.x -= box_y.x;      
				pbc_q.y -= box_y.y;
				pbc_q.z -= box_y.z;

			      } else continue;
			    }
			    if (a3 < 0) {
			      if (pbc_dirs.z == 1) {
				a3 = cell_dim.z-1;
				pbc_q.x -= box_z.x;      
				pbc_q.y -= box_z.y;
				pbc_q.z -= box_z.z;
			      } else continue;
			    }
			    if (a1 > cell_dim.x-1) {
			      if (pbc_dirs.x == 1) {
				a1 = 0; 
				pbc_q.x += box_x.x;      
				pbc_q.y += box_x.y;
				pbc_q.z += box_x.z;
			      } else continue;
			    }
			    if (a2 > cell_dim.y-1) {
			      if (pbc_dirs.y == 1) {
				a2 = 0; 
				pbc_q.x += box_y.x;      
				pbc_q.y += box_y.y;
				pbc_q.z += box_y.z;
			      } else continue;
			    }
			    if (a3 > cell_dim.z-1) {
			      if (pbc_dirs.z == 1) {
				a3 = 0; 
				pbc_q.x += box_z.x;      
				pbc_q.y += box_z.y;
				pbc_q.z += box_z.z;
			      } else continue;
			    }

			    /* Neighbour cell q (note that p==q ist possible) */

			    q = PTR_3D_V(cell_array,a1,a2,a3,cell_dim);

			    /* Calculate Indices of neighbour r of p */
			    b1 = i1 + k1; pbc_r.x = 0;
			    b2 = i2 + k2; pbc_r.y = 0;
			    b3 = i3 + k3; pbc_r.z = 0;

			    /* deal with periodic boundary conditions if necessary */
			    if (b1 < 0) {
			      if (pbc_dirs.x == 1) {
				b1 = cell_dim.x-1; 
				pbc_r.x -= box_x.x;      
				pbc_r.y -= box_x.y;
				pbc_r.z -= box_x.z;
			      } else continue;
			    }
			    if (b2 < 0) {
			      if (pbc_dirs.y == 1) {
				b2 = cell_dim.y-1;
				pbc_r.x -= box_y.x;      
				pbc_r.y -= box_y.y;
				pbc_r.z -= box_y.z;
			      } else continue;
			    }
			    if (b3 < 0) {
			      if (pbc_dirs.z == 1) {
				b3 = cell_dim.z-1;
				pbc_r.x -= box_z.x;      
				pbc_r.y -= box_z.y;
				pbc_r.z -= box_z.z;
			      } else continue;
			    }
			    if (b1 > cell_dim.x-1) {
			      if (pbc_dirs.x == 1) {
				b1 = 0; 
				pbc_r.x += box_x.x;      
				pbc_r.y += box_x.y;
				pbc_r.z += box_x.z;
			      } else continue; 
			    }
			    if (b2 > cell_dim.y-1) {
			      if (pbc_dirs.y == 1) {
				b2 = 0; 
				pbc_r.x += box_y.x;      
				pbc_r.y += box_y.y;
				pbc_r.z += box_y.z;
			      } else continue; 
			    }
			    if (b3 > cell_dim.z-1) {
			      if (pbc_dirs.z == 1) {
				b3 = 0; 
				pbc_r.x += box_z.x;      
				pbc_r.y += box_z.y;
				pbc_r.z += box_z.z;
			      } else continue;
			    } 

			    /* Neighbour cell r */

			    r = PTR_3D_V(cell_array, b1, b2, b3, cell_dim);

			    /* Calculate Indices of neighbour s of q */
			    c1 = j1 + l1; pbc_s.x = 0;
			    c2 = j2 + l2; pbc_s.y = 0;
			    c3 = j3 + l3; pbc_s.z = 0;

			    /* deal with periodic boundary conditions if necessary */
			    if (c1 < 0) {
			      if (pbc_dirs.x == 1) {
				c1 = cell_dim.x-1; 
				pbc_s.x -= box_x.x;      
				pbc_s.y -= box_x.y;
				pbc_s.z -= box_x.z;
			      } else continue;
			    }
			    if (c2 < 0) {
			      if (pbc_dirs.y == 1) {
				c2 = cell_dim.y-1;
				pbc_s.x -= box_y.x;      
				pbc_s.y -= box_y.y;
				pbc_s.z -= box_y.z;
			      } else continue;
			    }
			    if (c3 < 0) {
			      if (pbc_dirs.z == 1) {
				c3 = cell_dim.z-1;
				pbc_s.x -= box_z.x;      
				pbc_s.y -= box_z.y;
				pbc_s.z -= box_z.z;
			      } else continue;
			    }
			    if (c1 > cell_dim.x-1) {
			      if (pbc_dirs.x == 1) {
				c1 = 0; 
				pbc_s.x += box_x.x;      
				pbc_s.y += box_x.y;
				pbc_s.z += box_x.z;
			      } else continue; 
			    }
			    if (c2 > cell_dim.y-1) {
			      if (pbc_dirs.y == 1) {
				c2 = 0; 
				pbc_s.x += box_y.x;      
				pbc_s.y += box_y.y;
				pbc_s.z += box_y.z;
			      } else continue; 
			    }
			    if (c3 > cell_dim.z-1) {
			      if (pbc_dirs.z == 1) {
				c3 = 0; 
				pbc_s.x += box_z.x;      
				pbc_s.y += box_z.y;
				pbc_s.z += box_z.z;
			      } else continue;
			    } 

			    /* Neighbour cell s */

			    s = PTR_3D_V(cell_array, c1, c2, c3, cell_dim);

			    /* Do the work */
			    do_angle(p,q,r,s,pbc_q,pbc_r,pbc_s);
			  }
  
}



































































