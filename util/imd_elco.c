
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
* imd_elco -- calculate elastic constants
*
* uses cell division routines of imd_pair
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#ifndef ELCO
#define ELCO
#endif

#if defined(TERSOFF) || defined(STIWEB) || defined(KEATING)
#ifdef TWOD
#undef TWOD
#endif
#endif

/* Default is PAIR_POT. PAIR_POT is also needed for EAM */
#ifndef COVALENT
#ifndef PAIR_POT
#define PAIR_POT
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
  printf("%s [-r<nnn>] [-A<nnn>] [-v] [-e<nnn>] [-c] [-m] [-M] [-s] -p paramter-file]\n",progname); 
  
  exit(1); 
}

/*****************************************************************************
*
*  main
*
*****************************************************************************/

int main(int argc, char **argv)
{
  /* Read Parameters from parameter file */
  read_parameters(argc,argv);

#ifdef PAIR_POT
  /* Read pair potential */
  init_pair();
#endif
#ifdef EAM
  /* Read EAM potential */
  init_eam();
#endif
#ifdef TERSOFF
  /* Compute Tersoff parameters */
  init_tersoff();
#elif STIWEB
  /* Compute Stilliner-Weber parameters */
  init_stiweb();
#elif KEATING
  /* Compute Keating parameters */
  init_keating();
#endif

  /* Initialize cell data structures */
  init_cells();

  /* Read atoms */
  read_atoms(infilename);

#ifdef COVALENT
  /* Create neighbour tables */
  do_work(do_cell_pair);
#endif

  /* Initializations */
  init_elco();

#ifdef PAIR_POT
  /* Compute stress and elastic constants for pair potential */
  do_work(do_elco_pair);
#endif
#ifdef EAM
  /* Compute stress and elastic constants for EAM potential */
  do_work(do_elco_eam);
#endif
#ifdef TERSOFF
  /* Compute stress and elastic constants for Tersoff potential */
  do_elco_tersoff();
#elif STIWEB
  /* Compute stress and elastic constants for Stillinger-Weber potential */
  do_elco_stiweb();
#elif KEATING
  /* Compute stress and elastic constants for Keating potential */
  do_elco_keating();
#endif

  /* Write global data */
  write_data();

  if ( stresstens || moduli || all_moduli ) {
    /* Compute volumes of Voronoi cells */
    voronoi();
    
    if ( stresstens )
      /* Output stress tensor */
      write_stress();

    if ( moduli || all_moduli )
      /* Output tensor of elastic moduli */
      write_elco();
  }

  return 0;
}

/******************************************************************************
*
*  init_elco -- Initializations
*
******************************************************************************/

void init_elco(void) 
{

  int   i, j, k, l;
  cell  *p;

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
	for (l=0; l<p->n; ++l) 
	{
	  p->stress[l].xx = 0.0;
	  p->stress[l].xy = 0.0;
	  p->stress[l].yx = 0.0;
	  p->stress[l].yy = 0.0;
#ifndef TWOD
	  p->stress[l].xz = 0.0;
	  p->stress[l].yz = 0.0;
	  p->stress[l].zx = 0.0;
	  p->stress[l].zy = 0.0;
	  p->stress[l].zz = 0.0;
#endif
	  p->elco[l].c11  = 0.0;
	  p->elco[l].c12  = 0.0;
	  p->elco[l].c16  = 0.0;
	  p->elco[l].c22  = 0.0;
	  p->elco[l].c26  = 0.0;
	  p->elco[l].c66  = 0.0;
#ifndef TWOD
	  p->elco[l].c13  = 0.0;
	  p->elco[l].c14  = 0.0;
	  p->elco[l].c15  = 0.0;
	  p->elco[l].c23  = 0.0;
	  p->elco[l].c24  = 0.0;
	  p->elco[l].c25  = 0.0;
	  p->elco[l].c33  = 0.0;
	  p->elco[l].c34  = 0.0;
	  p->elco[l].c35  = 0.0;
	  p->elco[l].c36  = 0.0;
	  p->elco[l].c44  = 0.0;
	  p->elco[l].c45  = 0.0;
	  p->elco[l].c46  = 0.0;
	  p->elco[l].c55  = 0.0;
	  p->elco[l].c56  = 0.0;
#endif

#ifdef EAM
	  EAM_RHO(p,l)        = 0.0;
	  p->eam_stress[l].xx = 0.0;
	  p->eam_stress[l].xy = 0.0;
	  p->eam_stress[l].yx = 0.0;
	  p->eam_stress[l].yy = 0.0;
#ifndef TWOD
	  p->eam_stress[l].xz = 0.0;
	  p->eam_stress[l].yz = 0.0;
	  p->eam_stress[l].zx = 0.0;
	  p->eam_stress[l].zy = 0.0;
	  p->eam_stress[l].zz = 0.0;
#endif
	  p->eam_press[l]     = 0.0;
	  p->eam_bulkm[l]     = 0.0;
	  p->eam_dbulkm[l]    = 0.0;
#endif
	}
      }
}

#ifdef COVALENT

/******************************************************************************
*
*  do_cell_pair -- calulates neighbour tables
*
******************************************************************************/

void do_cell_pair(cell *p, cell *q, vektor pbc)
{
  int i, j;
  int jstart;
  int q_typ, p_typ;
  vektor d, tmp_d;
  real radius;

  /* For each atom in first cell */
  for (i=0; i<p->n; ++i) {

    tmp_d.x = p->ort[i].x - pbc.x;
    tmp_d.y = p->ort[i].y - pbc.y;
    tmp_d.z = p->ort[i].z - pbc.z;

    p_typ   = p->sorte[i];

    jstart = (p==q ? i+1 : 0);
    
    /* For each atom in neighbouring cell */
    for (j = jstart; j < q->n; ++j) {

      q_typ = q->sorte[j];
      
      /* Calculate distance  */
      d.x = q->ort[j].x - tmp_d.x;
      d.y = q->ort[j].y - tmp_d.y;
      d.z = q->ort[j].z - tmp_d.z;

      radius = sqrt(SPROD(d,d));

      /* Make neighbour tables */
#ifdef TERSOFF
      if (radius <= ter_r_cut[p_typ][q_typ])
#elif STIWEB
      if (radius <= sw_r_cut[p_typ][q_typ])
#elif KEATING
      if (radius <= keat_r_cut[p_typ][q_typ])
#endif
      {        
        neightab *neigh;
        real  *tmp_ptr;

        /* Update neighbour table of particle i */
        neigh = &p->neightab_array[i];

        if (neigh->n_max <= neigh->n) {
	  error("Neighbour table too small, increase neigh_len");
        }
        neigh->typ[neigh->n] = q_typ;
        neigh->cl [neigh->n] = q;
        neigh->num[neigh->n] = j;

        tmp_ptr  = &neigh->dist[3*neigh->n];
        *tmp_ptr = d.x; ++tmp_ptr; 
        *tmp_ptr = d.y; ++tmp_ptr; 
        *tmp_ptr = d.z;
        neigh->n++;

        /* Update neighbour table of particle j */
        neigh = &q->neightab_array[j];
        if (neigh->n_max <= neigh->n) {
	  error("Neighbour table too small, increase neigh_len");
        }
        neigh->typ[neigh->n] = p_typ;
        neigh->cl [neigh->n] = p;
        neigh->num[neigh->n] = i;
        tmp_ptr  = &neigh->dist[3*neigh->n];
        *tmp_ptr = -d.x; ++tmp_ptr; 
        *tmp_ptr = -d.y; ++tmp_ptr; 
        *tmp_ptr = -d.z;
        neigh->n++;
      }
    } /* for j */
  } /* for i */

}

#endif /* COVALENT */

/******************************************************************************
*
*  write_data -- writes global data to standard output
*
******************************************************************************/

void write_data(void)
{
  printf("\n");

  /* Number of atoms */
  printf("N = %d\n", natoms);

  /* Volume */
  printf("V = %.15f A^3\n", vol);

  /* Volume per atom */
  printf("v = %.15f A^3\n", vol/natoms);

  /* Potential energy per atom */
  printf("E = %.15f eV\n", epot/natoms);
  
  /* Hydrostatic pressure */
#ifndef PAIR_POT
  press = sigma.xx + sigma.yy + sigma.zz;
#endif
  printf("p = %.10f eV/A^3 = %.10f GPa\n", 
	 -press / ( 3.0 * vol ), -press / ( 3.0 * vol ) * CONV); 

  /* Bulk modulus */
#ifndef PAIR_POT
  bulkm = c.c11 + c.c22 + c.c33 + 2.0 * ( c.c12 + c.c13 + c.c23 );
#endif
  printf("B = %.10f eV/A^3 = %.10f GPa\n", 
	 bulkm / ( 9.0 * vol ), bulkm / (9.0 * vol ) * CONV);

  /* Pressure derivative of the bulk modulus */
  if ( bulkm != 0.0 ) {
    dbulkm_dp = 1.0 / 3.0 - dbulkm_dp / ( 3.0 * bulkm );
    printf("B'= %.10f \n", dbulkm_dp);
  }

  /* Stress tensor*/
  printf("\n");
  printf("sigma_xx = %.10f eV/A^3 = %.10f GPa\n", 
	 sigma.xx / vol, sigma.xx / vol * CONV ); 
  printf("sigma_yy = %.10f eV/A^3 = %.10f GPa\n", 
	 sigma.yy / vol, sigma.yy / vol * CONV ); 
  printf("sigma_zz = %.10f eV/A^3 = %.10f GPa\n", 
	 sigma.zz / vol, sigma.zz / vol * CONV ); 
  printf("sigma_yz = %.10f eV/A^3 = %.10f GPa\n", 
	 sigma.yz / vol, sigma.yz / vol * CONV ); 
  printf("sigma_zx = %.10f eV/A^3 = %.10f GPa\n", 
	 sigma.zx / vol, sigma.zx / vol * CONV ); 
  printf("sigma_xy = %.10f eV/A^3 = %.10f GPa\n", 
	 sigma.xy / vol, sigma.xy / vol * CONV ); 

  /* Tensor of elastic moduli */
  printf("\n");
  printf("C_11 = %2.10f eV/A^3 = %2.10f GPa\n", 
	 c.c11 / vol, c.c11 / vol * CONV );
  printf("C_12 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c12 / vol, c.c12 / vol * CONV );
  printf("C_13 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c13 / vol, c.c13 / vol * CONV );
  printf("C_22 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c22 / vol, c.c22 / vol * CONV );
  printf("C_23 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c23 / vol, c.c23 / vol * CONV );
  printf("C_33 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c33 / vol, c.c33 / vol * CONV );
  printf("C_44 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c44 / vol, c.c44 / vol * CONV );
  printf("C_45 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c45 / vol, c.c45 / vol * CONV );
  printf("C_46 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c46 / vol, c.c46 / vol * CONV );
  printf("C_55 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c55 / vol, c.c55 / vol * CONV );
  printf("C_56 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c56 / vol, c.c56 / vol * CONV );
  printf("C_66 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c66 / vol, c.c66 / vol * CONV );
  printf("C_14 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c14 / vol, c.c14 / vol * CONV );
  printf("C_15 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c15 / vol, c.c15 / vol * CONV );
  printf("C_16 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c16 / vol, c.c16 / vol * CONV );
  printf("C_24 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c24 / vol, c.c24 / vol * CONV );
  printf("C_25 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c25 / vol, c.c25 / vol * CONV );
  printf("C_26 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c26 / vol, c.c26 / vol * CONV );
  printf("C_34 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c34 / vol, c.c34 / vol * CONV );
  printf("C_35 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c35 / vol, c.c35 / vol * CONV );
  printf("C_36 = %.10f eV/A^3 = %.10f GPa\n", 
	 c.c36 / vol, c.c36 / vol * CONV );
  printf("\n");
}
    
/******************************************************************************
*
*  write_stress -- writes data to *.stress file
*
******************************************************************************/

void write_stress(void)
{
  FILE *out;
  str255 fname;
  int i,j,k,l;
  cell *p; 
  real tmp;

  sprintf(fname,"%s.stress",infilename);
  
  out = fopen(fname,"w");
  if (NULL == out) error("Cannot open stress file.");

  fprintf(out, "#No type\t x        y        z        s_xx         s_yy         s_zz         s_yz         s_zx         s_xy         Vol_Voronoi\n");

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
	{
	  if ( p->vol[l] > 0.0 ) { 
	    tmp = 1.0 / p->vol[l];
#ifndef TWOD	  
	    fprintf(out, "%d %d %f %f %f %.10f %.10f %.10f %.10f %.10f %.10f %f\n", 
		    p->nummer[l], p->sorte[l], p->ort[l].x, p->ort[l].y, p->ort[l].z, 
		    p->stress[l].xx * tmp, p->stress[l].yy * tmp, p->stress[l].zz * tmp, 
		    p->stress[l].yz * tmp, p->stress[l].zx * tmp, p->stress[l].xy * tmp, 
		    p->vol[l]);
#else
	    fprintf(out, "%d %d %f %f %.10f %.10f %.10f %f\n", 
		    p->nummer[l], p->sorte[l], p->ort[l].x, p->ort[l].y, 
		    p->stress[l].xx * tmp, p->stress[l].yy * tmp, p->stress[l].zz * tmp,  
		    p->vol[l]);
#endif
	  }
	}
      }

  fclose(out);

  /* Write Statistics */
  if( atomcount > 0 )
    {
      printf("\n");
      printf("Maximal number of neighbour atoms:           %d\n", maxneigh );
      printf("Average number of neighbour atoms:           %.2f\n\n", 
	     (real)(sumneigh)/atomcount );
      printf("Maximal number of vertices of Voronoi cells: %d\n", maxvert );
      printf("Average number of vertices:                  %.2f\n\n", 
	     (real)(sumvert)/atomcount );
      printf("Maximal number of edges of Voronoi cells:    %d\n", maxedges );
      printf("Average number of edges:                     %.2f\n\n", 
	     (real)(sumedges)/atomcount );
#ifndef TWOD
      printf("Maximal number of faces of Voronoi cells:    %d\n", maxfaces );
      printf("Average number of faces:                     %.2f\n\n", 
	     (real)(sumfaces)/atomcount );
#endif
      printf("Total number of atoms:   %d\n", natoms);
      printf("Number of omitted atoms: %d (%.2f %%) \n\n", natoms - atomcount, 
	     (real)(natoms-atomcount)/natoms*100.0 );
    }
  else
    printf("No Voronoi cell found.\n");

}

/******************************************************************************
*
*  write_elco -- writes data *.elco file
*
******************************************************************************/

void write_elco(void)
{
  FILE *out;
  str255 fname;
  int i,j,k,l;
  cell *p; 
  real tmp, bulkmod;

  sprintf(fname,"%s.elco",infilename);
  
  out = fopen(fname,"w");
  if (NULL == out) error("Cannot open elco file.");

  if ( all_moduli == 1 )
#ifndef TWOD
    fprintf(out, "#No type x        y        z         B       c_11     c_12     c_13     c_22     c_23     c_33     c_44     c_45     c_46     c_55     c_56     c_66     c_14     c_15     c_16     c_24     c_25     c_26     c_34     c_35     c_36\n");
#else
    fprintf(out, "#No type x        y        B       c_11     c_12     c_22     c_66     c_16     c_26\n");
#endif         
  else
#ifndef TWOD
    fprintf(out, "#No type x        y        z        B           c_11         c_12         c_44\n");
#else
    fprintf(out, "#No type x        y        B           c_11         c_12\n");
#endif

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
	{
	  if ( p->vol[l] > 0.0 ) {
	    tmp = 1.0 / p->vol[l];
#ifndef TWOD
	    bulkmod = ( p->elco[l].c11 + p->elco[l].c22 + p->elco[l].c33 
		       + 2.0 * ( p->elco[l].c12 + p->elco[l].c13 + p->elco[l].c23 ) ) / 9.0;
#else
	    bulkmod = ( p->elco[l].c11 + p->elco[l].c22 
		       + 2.0 * p->elco[l].c12 ) / 9.0;  /* ?? */
#endif	    
	    if ( all_moduli == 1 ) 
#ifndef TWOD
	      fprintf(out, "%d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", 
		      p->nummer[l], p->sorte[l], p->ort[l].x, p->ort[l].y, p->ort[l].z, 
		      bulkmod * tmp, p->elco[l].c11 * tmp, p->elco[l].c12 * tmp, 
		      p->elco[l].c13 * tmp, p->elco[l].c22 * tmp, p->elco[l].c23 * tmp, 
		      p->elco[l].c33 * tmp, p->elco[l].c44 * tmp, p->elco[l].c45 * tmp, 
		      p->elco[l].c46 * tmp, p->elco[l].c55 * tmp, p->elco[l].c56 * tmp, 
		      p->elco[l].c66 * tmp, p->elco[l].c14 * tmp, p->elco[l].c15 * tmp, 
		      p->elco[l].c16 * tmp, p->elco[l].c24 * tmp, p->elco[l].c25 * tmp, 
		      p->elco[l].c26 * tmp, p->elco[l].c34 * tmp, p->elco[l].c35 * tmp, 
		      p->elco[l].c36 * tmp );
#else
	      fprintf(out, "%d %d %f %f %f %f %f %f %f %f %f\n", 
		      p->nummer[l], p->sorte[l], p->ort[l].x, p->ort[l].y, 
		      bulkmod * tmp, p->elco[l].c11 * tmp, p->elco[l].c12 * tmp, 
		      p->elco[l].c22 * tmp,  
		      p->elco[l].c66 * tmp, p->elco[l].c16 * tmp, p->elco[l].c26 * tmp );
#endif
	    else 
#ifndef TWOD
	      fprintf(out, "%d %d %f %f %f %.10f %.10f %.10f %.10f\n", 
		      p->nummer[l], p->sorte[l], p->ort[l].x, p->ort[l].y, p->ort[l].z, 
		      bulkmod * tmp, p->elco[l].c11 * tmp, p->elco[l].c12 * tmp, 
		      p->elco[l].c44 * tmp );
#else
	      fprintf(out, "%d %d %f %f %.10f %.10f %.10f\n", 
		      p->nummer[l], p->sorte[l], p->ort[l].x, p->ort[l].y, 
		      bulkmod * tmp, p->elco[l].c11 * tmp, p->elco[l].c12 * tmp );
#endif
	  }
	}
      }

  fclose(out);

}

