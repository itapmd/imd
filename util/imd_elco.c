
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

#define ELCO

#ifndef STIWEB
#undef TERSOFF
#define TERSOFF
#endif

#if defined(TERSOFF) || defined(STIWEB)
#undef TWOD
#endif

#include "util.h"
#include "imd_stress.c"

/******************************************************************************
*
*  Usage -- educate users
*
*  Compilation: gcc -O [-DTERSOFF] [-DSTIWEB] [-DSINGLE] imd_elco.c -lm 
*
******************************************************************************/

void usage(void)
{ 
  printf("%s [-r<nnn>] [-A<nnn>] [-e<nnn>] [-c] [-m] [-M] [-s] -p paramter-file]\n",progname); 
  
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

#ifdef TERSOFF
  /* Compute Tersoff parameters */
  init_tersoff();
#elif defined(STIWEB)
  /* Compute Stilliner-Weber parameters */
  init_stiweb();
#endif

  /* Initialize cell data structures */
  init_cells();

  /* Read atoms */
  read_atoms(infilename);

  /* Create neighbour tables */
  do_work(do_cell_pair);

#ifdef TERSOFF
  /* Compute stress and elastic constants for Tersoff potential*/
  do_elco_tersoff();
#elif defined(STIWEB)
  /* Compute stress and elastic constants for Stillinger-Weber potential*/
  do_elco_stiweb();
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

#ifdef TERSOFF

/******************************************************************************
*
*  init_tersoff
*
******************************************************************************/

void init_tersoff(void) {

  int i, j, n = 0;
  real tmp;

  /* parameters for more than one atom type */
  for (i=0; i<ntypes; i++) {
    ter_c2[i] = ters_c[i] * ters_c[i];
    ter_d2[i] = ters_d[i] * ters_d[i];
    for (j=i; j<ntypes; j++) {
      ter_r_cut[i][j]  = ter_r_cut[j][i]  = ters_r_cut[n];
      ter_r2_cut[i][j] = ter_r2_cut[j][i] = ter_r_cut[i][j] * ter_r_cut[i][j];
      ter_r0[i][j]     = ter_r0[j][i]     = ters_r0[n];
      ter_a[i][j]      = ter_a[j][i]      = ters_a[n];
      ter_b[i][j]      = ter_b[j][i]      = ters_b[n];
      ter_la[i][j]     = ter_la[j][i]     = ters_la[n];
      ter_mu[i][j]     = ter_mu[j][i]     = ters_mu[n];
      ++n;      
    }
  }

  for (i=0; i<ntypes; i++) ter_chi[i][i] = 1.0;
  if ( ntypes>1 ) {
    for (i=0; i<(ntypes-1); i++)
      for (j=(i+1); j<ntypes; j++) {
        ter_chi[i][j] = ter_chi[j][i] 
                      = ters_chi[i * ( 2 * ntypes - i - 3 ) / 2 + j - 1]; 
      }
  }

  for (i=0; i<ntypes; i++) ter_om[i][i] = 1.0;
  if ( ntypes>1 ) {
    for (i=0; i<(ntypes-1); i++)
      for (j=(i+1); j<ntypes; j++) {
        ter_om[i][j] = ter_om[j][i] 
                     = ters_om[i * ( 2 * ntypes - i - 3 ) / 2 + j - 1]; 
      }
  }
 
  if( r_cell != -1.0)
    r2_cut = SQR(r_cell);
  else { 
    tmp = 0.0;
    for (i=0; i<ntypes; ++i)
      for (j=0; j<ntypes; ++j)
	tmp = MAX( tmp, ter_r2_cut[i][j] );
    r2_cut = MAX(r2_cut,tmp);
  }

}

#elif defined(STIWEB)

/******************************************************************************
*
*  init_stiweb
*
******************************************************************************/

void init_stiweb(void) {

  int  i, j, k, n, m;
  real tmp;

  /* parameters for more than one atom type */
  n = 0; m = 0;
  for (i=0; i<ntypes; i++) 
    for (j=i; j<ntypes; j++) {
      sw_a[i][j]  = sw_a[j][i]  = stiweb_a[n];
      sw_b[i][j]  = sw_b[j][i]  = stiweb_b[n];
      sw_p[i][j]  = sw_p[j][i]  = stiweb_p[n];
      sw_q[i][j]  = sw_q[j][i]  = stiweb_q[n];
      sw_a1[i][j] = sw_a1[j][i] = stiweb_a1[n];
      sw_de[i][j] = sw_de[j][i] = stiweb_de[n];
      sw_ga[i][j] = sw_ga[j][i] = stiweb_ga[n];
      sw_a2[i][j] = sw_a2[j][i] = stiweb_a2[n];
      n++;
      for (k=0; k<ntypes; k++) {
	sw_la[k][i][j] = sw_la[k][j][i] = stiweb_la[m];
	m++;
      }
    }

  for (i=0; i<ntypes; ++i)
    for (j=0; j<ntypes; ++j) 
      sw_2_a1[i][j] = sw_a1[i][j] * sw_a1[i][j];

  for (i=0; i<ntypes; ++i)
    for (j=0; j<ntypes; ++j) 
      sw_r_cut[i][j] = MAX(sw_a1[i][j], sw_a2[i][j] );
 
  if( r_cell != -1.0)
    r2_cut = SQR(r_cell);
  else { 
    tmp = 0.0;
    for (i=0; i<ntypes; ++i)
      for (j=0; j<ntypes; ++j) {
	tmp = MAX( tmp, sw_a1[i][j] );
	tmp = MAX( tmp, sw_a2[i][j] );
      }
    r2_cut = MAX(r2_cut,tmp*tmp);
  }
}

#endif

/******************************************************************************
*
*  do_cell_pair calulates neighbor tables
*
******************************************************************************/

void do_cell_pair(cell *p, cell *q, vektor pbc)
{
  int i, j, k;
  int jstart, jend;
  int q_typ, p_typ, column;
  vektor d, tmp_d;
  real *qptr, radius;

  /* For each atom in first cell */
  for (i=0; i<p->n; ++i) {

    tmp_d.x = p->ort[i].x - pbc.x;
    tmp_d.y = p->ort[i].y - pbc.y;
    tmp_d.z = p->ort[i].z - pbc.z;
    p_typ   = p->sorte[i];

    jstart = (p==q ? i+1 : 0);
    
    /* For each atom in neighboring cell */
    for (j = jstart; j < q->n; ++j) {

      q_typ = q->sorte[j];
      
      /* Calculate distance  */
      d.x = q->ort[j].x - tmp_d.x;
      d.y = q->ort[j].y - tmp_d.y;
      d.z = q->ort[j].z - tmp_d.z;

      radius = sqrt(SPROD(d,d));

      /* Make neighbor tables */
#ifdef TERSOFF
      if (radius <= ter_r_cut[p_typ][q_typ])
#elif STIWEB
      if (radius <= sw_r_cut[p_typ][q_typ])
#endif
      {        
        neightab *neigh;
        real  *tmp_ptr;

        /* Update neighbor table of particle i */
        neigh = &p->neightab_array[i];

        if (neigh->n_max <= neigh->n) {
	  error("Neighbor table too small, increase neigh_len");
        }
        neigh->typ[neigh->n] = q_typ;
        neigh->cl [neigh->n] = q;
        neigh->num[neigh->n] = j;

        tmp_ptr  = &neigh->dist[3*neigh->n];
        *tmp_ptr = d.x; ++tmp_ptr; 
        *tmp_ptr = d.y; ++tmp_ptr; 
        *tmp_ptr = d.z;
        neigh->n++;

        /* Update neighbor table of particle j */
        neigh = &q->neightab_array[j];
        if (neigh->n_max <= neigh->n) {
	  error("Neighbor table too small, increase neigh_len");
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

/******************************************************************************
*
*  write_data writes global data to standard output
*
******************************************************************************/

void write_data(void)
{
  real p, bulkm;

  /* Potential energy */
  printf("\nE_pot = %.15f eV\n", epot/natoms);
  
  /* Hydrostatic pressure */
  p = - ( sigma.xx + sigma.yy + sigma.zz ) / 3.0 / vol;
  printf("p = %.10f eV/A^3 = %.10f GPa\n", p, p * 160.218); 

  /* Bulk modulus */
  bulkm = ( c.c11 + c.c22 + c.c33 + 2.0 * ( c.c12 + c.c13 + c.c23 ) ) / 9.0;
  printf("B = %.10f eV/A^3 = %.10f GPa\n", bulkm/vol, bulkm/vol*160.218);

  /* Pressure derivative of the bulk modulus */
  dbulkm_dp = 1.0 / 3.0 - dbulkm_dp / ( 27.0 * bulkm );
  printf("B' = %.10f \n", dbulkm_dp);

  /* Stress tensor*/
  printf("\n");
  printf("sigma_xx = %.10f eV/A^3 = %.10f GPa\n", sigma.xx / vol, sigma.xx / vol * 160.218); 
  printf("sigma_yy = %.10f eV/A^3 = %.10f GPa\n", sigma.yy / vol, sigma.yy / vol * 160.218); 
  printf("sigma_zz = %.10f eV/A^3 = %.10f GPa\n", sigma.zz / vol, sigma.zz / vol * 160.218); 
  printf("sigma_yz = %.10f eV/A^3 = %.10f GPa\n", sigma.yz / vol, sigma.yz / vol * 160.218); 
  printf("sigma_zx = %.10f eV/A^3 = %.10f GPa\n", sigma.zx / vol, sigma.zx / vol * 160.218); 
  printf("sigma_xy = %.10f eV/A^3 = %.10f GPa\n", sigma.xy / vol, sigma.xy / vol * 160.218); 
  printf("\n");

  /* Tensor of elastic moduli */
  printf("C_11 = %.10f eV/A^3 = %.10f GPa\n", c.c11 / vol, c.c11 / vol * 160.218 );
  printf("C_12 = %.10f eV/A^3 = %.10f GPa\n", c.c12 / vol, c.c12 / vol * 160.218 );
  printf("C_13 = %.10f eV/A^3 = %.10f GPa\n", c.c13 / vol, c.c13 / vol * 160.218 );
  printf("C_22 = %.10f eV/A^3 = %.10f GPa\n", c.c22 / vol, c.c22 / vol * 160.218 );
  printf("C_23 = %.10f eV/A^3 = %.10f GPa\n", c.c23 / vol, c.c23 / vol * 160.218 );
  printf("C_33 = %.10f eV/A^3 = %.10f GPa\n", c.c33 / vol, c.c33 / vol * 160.218 );
  printf("C_44 = %.10f eV/A^3 = %.10f GPa\n", c.c44 / vol, c.c44 / vol * 160.218 );
  printf("C_45 = %.10f eV/A^3 = %.10f GPa\n", c.c45 / vol, c.c45 / vol * 160.218 );
  printf("C_46 = %.10f eV/A^3 = %.10f GPa\n", c.c46 / vol, c.c46 / vol * 160.218 );
  printf("C_55 = %.10f eV/A^3 = %.10f GPa\n", c.c55 / vol, c.c55 / vol * 160.218 );
  printf("C_56 = %.10f eV/A^3 = %.10f GPa\n", c.c56 / vol, c.c56 / vol * 160.218 );
  printf("C_66 = %.10f eV/A^3 = %.10f GPa\n", c.c66 / vol, c.c66 / vol * 160.218 );
  printf("C_14 = %.10f eV/A^3 = %.10f GPa\n", c.c14 / vol, c.c14 / vol * 160.218 );
  printf("C_15 = %.10f eV/A^3 = %.10f GPa\n", c.c15 / vol, c.c15 / vol * 160.218 );
  printf("C_16 = %.10f eV/A^3 = %.10f GPa\n", c.c16 / vol, c.c16 / vol * 160.218 );
  printf("C_24 = %.10f eV/A^3 = %.10f GPa\n", c.c24 / vol, c.c24 / vol * 160.218 );
  printf("C_25 = %.10f eV/A^3 = %.10f GPa\n", c.c25 / vol, c.c25 / vol * 160.218 );
  printf("C_26 = %.10f eV/A^3 = %.10f GPa\n", c.c26 / vol, c.c26 / vol * 160.218 );
  printf("C_34 = %.10f eV/A^3 = %.10f GPa\n", c.c34 / vol, c.c34 / vol * 160.218 );
  printf("C_35 = %.10f eV/A^3 = %.10f GPa\n", c.c35 / vol, c.c35 / vol * 160.218 );
  printf("C_36 = %.10f eV/A^3 = %.10f GPa\n", c.c36 / vol, c.c36 / vol * 160.218 );
  printf("\n");

}
    
/******************************************************************************
*
*  write_stress writes data to *.stress file
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

  fprintf(out, "# x        y        z        s_xx         s_yy         s_zz         s_yz         s_zx         s_xy         Vol_Voronoi\n");

  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
      for (k=0; k < cell_dim.z; ++k)
      {
	p = PTR_3D_V(cell_array,i,j,k,cell_dim);

	for (l=0;l<p->n; ++l)
	{
	  if ( p->vol[l] > 0.0 ) { 
	    tmp = 1.0 / p->vol[l];
	  
	    fprintf(out, "%f %f %f %.10f %.10f %.10f %.10f %.10f %.10f %f\n", p->ort[l].x, p->ort[l].y, p->ort[l].z, p->stress[l].xx * tmp, p->stress[l].yy * tmp, p->stress[l].zz * tmp , p->stress[l].yz * tmp, p->stress[l].zx * tmp, p->stress[l].xy * tmp, p->vol[l]);
	  }
	}
      }

  fclose(out);

}

/******************************************************************************
*
*  write_elco writes data *.elco file
*
******************************************************************************/

void write_elco(void)
{
  FILE *out;
  str255 fname;
  int i,j,k,l;
  int number=0;
  cell *p; 
  real tmp, bulkm;

  sprintf(fname,"%s.elco",infilename);
  
  out = fopen(fname,"w");
  if (NULL == out) error("Cannot open elco file.");

  if ( all_moduli == 1 )
    fprintf(out, "#  x        y        z         B       c_11     c_12     c_13     c_22     c_23     c_33     c_44     c_45     c_46     c_55     c_56     c_66     c_14     c_15     c_16     c_24     c_25     c_26     c_34     c_35     c_36\n");         
  else
    fprintf(out, "#  x        y        z        B           c_11         c_12         c_44\n");

  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
      for (k=0; k < cell_dim.z; ++k)
      {
	p = PTR_3D_V(cell_array,i,j,k,cell_dim);

	for (l=0;l<p->n; ++l)
	{
	  if ( p->vol[l] > 0.0 ) {
	    tmp = 1.0 / p->vol[l];
	    bulkm = ( p->elco[l].c11 + p->elco[l].c22 + p->elco[l].c33 + 2.0 * ( p->elco[l].c12 + p->elco[l].c13 + p->elco[l].c23 ) ) / 9.0;
	    
	    if ( all_moduli == 1 ) 
	      fprintf(out, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", p->ort[l].x, p->ort[l].y, p->ort[l].z, bulkm * tmp, p->elco[l].c11 * tmp, p->elco[l].c12 * tmp, p->elco[l].c13 * tmp, p->elco[l].c14 * tmp, p->elco[l].c15 * tmp, p->elco[l].c16 * tmp, p->elco[l].c22 * tmp, p->elco[l].c23 * tmp, p->elco[l].c24 * tmp, p->elco[l].c25 * tmp, p->elco[l].c26 * tmp, p->elco[l].c33 * tmp, p->elco[l].c34 * tmp, p->elco[l].c35 * tmp, p->elco[l].c36 * tmp, p->elco[l].c44 * tmp, p->elco[l].c45 * tmp, p->elco[l].c46 * tmp, p->elco[l].c55 * tmp, p->elco[l].c56 * tmp, p->elco[l].c66 * tmp );
	    else 
	      fprintf(out, "%f %f %f %.10f %.10f %.10f %.10f\n", p->ort[l].x, p->ort[l].y, p->ort[l].z, bulkm * tmp, p->elco[l].c11 * tmp, p->elco[l].c12 * tmp, p->elco[l].c44 * tmp );
	  }
	}
      }

  fclose(out);

}

#ifdef TERSOFF

/******************************************************************************
*
*  do_elco_tersoff -- computes stress tensor and elastic constants using the 
*             Tersoff potential
*
******************************************************************************/

void do_elco_tersoff(void)
{
  cell   *p;
  int    ic, jc, kc;
  real   *r, *fc, *dfc, *ddfc, *dddfc;
  vektor *d;
  neightab *neigh;
  vektor dcos_j, dcos_k, dzeta_j, dphi_j;
  tensor ddcos_jj, ddcos_jk, ddcos_kk, ddzeta_jj, ddphi_jj, *ddzeta_kk;
  real   *dzeta_bk, *ddzeta_bkl, *dddzeta_bklm;
  real   *db_k, *ddb_lk, *dddb_lkm;
  vektor *dzeta_k;
  tensor *ddphi_jk, *ddphi_lk, *ddzeta_jk;
  cell   *jcell, *kcell;
  int    i, j, k, l, m, p_typ, j_typ, k_typ, jnum;
  real   *tmpptr;
  real   pot_zwi, tmp_grad;
  real   cos_theta, cut_tmp, cut_tmp_j;
  real   zeta, g_theta, b_ij;
  real   tmp_jj, tmp_jk, tmp_kk, tmp_j2, tmp_k2;
  real   phi_r, phi_a;
  real   dphi_r, dphi_a, ddphi_r, ddphi_a, dddphi_r, dddphi_a;
  real   tmp_1, tmp_2, tmp_3, tmp_4, tmp_5, tmp_6;
  real   tmp1_zeta, tmp2_zeta, tmp3_zeta, tmp4_zeta, tmp5_zeta;
  real   tmp1_phi, tmp2_phi, tmp3_phi, tmp4_phi, tmp5_phi;
  real   tmp6_phi, tmp7_phi, tmp8_phi, tmp9_phi;
  real   tmp_b1, tmp_b2, tmp_b3;
  real   tmp, bulkm;

  d            = (vektor *) malloc( neigh_len                         * sizeof(vektor) );
  r            = (real *)   malloc( neigh_len                         * sizeof(real)   );
  fc           = (real *)   malloc( neigh_len                         * sizeof(real)   );
  dfc          = (real *)   malloc( neigh_len                         * sizeof(real)   );
  ddfc         = (real *)   malloc( neigh_len                         * sizeof(real)   );
  dddfc        = (real *)   malloc( neigh_len                         * sizeof(real)   );
  dzeta_k      = (vektor *) malloc( neigh_len                         * sizeof(vektor) );
  dzeta_bk     = (real   *) malloc( neigh_len                         * sizeof(real) );
  ddzeta_bkl   = (real   *) malloc( neigh_len                         * sizeof(real) );
  dddzeta_bklm = (real   *) malloc( neigh_len                         * sizeof(real) );
  db_k         = (real   *) malloc( neigh_len                         * sizeof(real) );
  ddb_lk       = (real   *) malloc( neigh_len * neigh_len             * sizeof(real) );
  dddb_lkm     = (real   *) malloc( neigh_len * neigh_len * neigh_len * sizeof(real) );
  ddphi_jk     = (tensor *) malloc( neigh_len                         * sizeof(tensor) );
  ddphi_lk     = (tensor *) malloc( neigh_len * neigh_len             * sizeof(tensor) );
  ddzeta_jk    = (tensor *) malloc( neigh_len                         * sizeof(tensor) );
  ddzeta_kk    = (tensor *) malloc( neigh_len                         * sizeof(tensor) );

  if ((d==NULL) || (r==NULL) || (fc==NULL) || (dfc==NULL) 
      || (ddfc==NULL) || (dddfc==NULL)
      || (dzeta_k==NULL) || (dzeta_bk==NULL) || (ddzeta_bkl==NULL) 
      || (ddzeta_jk==NULL) || (ddzeta_kk==NULL) 
      || (db_k==NULL) || (ddb_lk==NULL) || (dddb_lkm==NULL) 
      || (ddphi_jk==NULL) || (ddphi_lk==NULL) )
    error("cannot allocate memory for temporary neighbor data");

  /*     k
          \
           \
	    i----j  */

  /* Initializations */
  for (i=0; i<neigh_len; i++)
    for (j=0; j<neigh_len; ++j) 
    {
      ddphi_lk I(i,j) .xx = 0.0;
      ddphi_lk I(i,j) .xy = 0.0;
      ddphi_lk I(i,j) .yy = 0.0;
      ddphi_lk I(i,j) .yz = 0.0;
      ddphi_lk I(i,j) .zx = 0.0;
      ddphi_lk I(i,j) .zz = 0.0;
    }

  for (ic=0; ic < cell_dim.x; ++ic)
    for (jc=0; jc < cell_dim.y; ++jc)
      for (kc=0; kc < cell_dim.z; ++kc)
      {
        p = PTR_3D_V(cell_array,ic,jc,kc,cell_dim);
	
	for (i=0; i<p->n; ++i) 
	{
	  p->elco[i].c11 = 0.0;
	  p->elco[i].c12 = 0.0;
	  p->elco[i].c13 = 0.0;
	  p->elco[i].c14 = 0.0;
	  p->elco[i].c15 = 0.0;
	  p->elco[i].c16 = 0.0;
	  p->elco[i].c22 = 0.0;
	  p->elco[i].c23 = 0.0;
	  p->elco[i].c24 = 0.0;
	  p->elco[i].c25 = 0.0;
	  p->elco[i].c26 = 0.0;
	  p->elco[i].c33 = 0.0;
	  p->elco[i].c34 = 0.0;
	  p->elco[i].c35 = 0.0;
	  p->elco[i].c36 = 0.0;
	  p->elco[i].c44 = 0.0;
	  p->elco[i].c45 = 0.0;
	  p->elco[i].c46 = 0.0;
	  p->elco[i].c55 = 0.0;
	  p->elco[i].c56 = 0.0;
	  p->elco[i].c66 = 0.0;

	  p->stress[i].xx = 0.0;
	  p->stress[i].xy = 0.0;
	  p->stress[i].xz = 0.0;
	  p->stress[i].yx = 0.0;
	  p->stress[i].yy = 0.0;
	  p->stress[i].yz = 0.0;
	  p->stress[i].zx = 0.0;
	  p->stress[i].zy = 0.0;
	  p->stress[i].zz = 0.0;

	}
      }

  /* For each cell */
  for (ic=0; ic < cell_dim.x; ++ic)
    for (jc=0; jc < cell_dim.y; ++jc)
      for (kc=0; kc < cell_dim.z; ++kc)
      {
        p = PTR_3D_V(cell_array,ic,jc,kc,cell_dim);

	/* For each atom in cell */
	for (i=0; i<p->n; ++i) 
	{

	  p_typ   = p->sorte[i];
	  neigh   = &p->neightab_array[i];

	  /* Construct some data for all neighbors */
	  tmpptr = neigh->dist;
	  for (j=0; j<neigh->n; ++j) 
	  {

	    /* Type, distance vector, radii */
	    j_typ   = neigh->typ[j];
	    d[j].x  = *tmpptr++;
	    d[j].y  = *tmpptr++;
	    d[j].z  = *tmpptr++;
	    r[j]    = sqrt(SPROD(d[j],d[j]));

	    /* Cutoff function and its derivatives */
	    cut_tmp   = M_PI / ( ter_r_cut[p_typ][j_typ] - ter_r0[p_typ][j_typ] );
	    cut_tmp_j = cut_tmp * ( r[j] - ter_r0[p_typ][j_typ] );
	    if ( r[j] < ter_r0[p_typ][j_typ] ) {
	      fc[j]    = 1.0; 
	      dfc[j]   = 0.0;
	      ddfc[j]  = 0.0;
	      dddfc[j] = 0.0;
	    }
	    else if ( r[j] > ter_r_cut[p_typ][j_typ] ) {
	      fc[j]    = 0.0;
	      dfc[j]   = 0.0;
	      ddfc[j]  = 0.0;
	      dddfc[j] = 0.0;
	    }
	    else {
	      fc[j]    =   0.5 * ( 1.0 + cos( cut_tmp_j ) );
	      dfc[j]   = - 0.5 * cut_tmp * sin( cut_tmp_j );
	      ddfc[j]  = - 0.5 * cut_tmp * cut_tmp * cos( cut_tmp_j );
	      dddfc[j] = - cut_tmp * cut_tmp * dfc[j];
	    }      
	  } /* j */

          /*********************************************************************/

	  /* For each neighbor of i */
	  for (j=0; j<neigh->n; ++j)
	  {

	    j_typ  = neigh->typ[j];
	    jcell  = (cell *) neigh->cl [j];
	    jnum   = neigh->num[j];
	    tmp_jj = 1 / ( r[j] * r[j] );

	    /* Initializations */
	    zeta         = 0.0;     
	    dzeta_j.x    = 0.0; dzeta_j.y    = 0.0; dzeta_j.z    = 0.0;
	    ddzeta_jj.xx = 0.0; ddzeta_jj.xy = 0.0; ddzeta_jj.yy = 0.0;
	    ddzeta_jj.yz = 0.0; ddzeta_jj.zx = 0.0; ddzeta_jj.zz = 0.0;

	    /* For each neighbor of i other than j */
	    for (k=0; k<neigh->n; ++k) if (k!=j) 
	    {

	      k_typ = neigh->typ[k];

	      /* Angular term */
	      tmp_jk    = 1 / ( r[j] * r[k] );  
	      cos_theta = SPROD(d[j],d[k]) * tmp_jk;
	      tmp_1     = ters_h[p_typ] - cos_theta;
	      tmp_2     = 1 / ( ter_d2[p_typ] + tmp_1 * tmp_1 );
	      g_theta   = 1 + ter_c2[p_typ]/ter_d2[p_typ] - ter_c2[p_typ] * tmp_2;
	      
	      /* zeta */
	      zeta  += fc[k] * ter_om[p_typ][k_typ] * g_theta; 

	      /* tmp variables */
	      tmp_j2 = cos_theta / ( r[j] * r[j] );
	      tmp_k2 = cos_theta / ( r[k] * r[k] );
	      tmp_kk = 1 /  ( r[k] * r[k] );
	      
	      /* Derivatives of cos(theta) */
	      dcos_j.x = tmp_jk * d[k].x - tmp_j2 * d[j].x;
	      dcos_j.y = tmp_jk * d[k].y - tmp_j2 * d[j].y;
	      dcos_j.z = tmp_jk * d[k].z - tmp_j2 * d[j].z;
	      
	      dcos_k.x = tmp_jk * d[j].x - tmp_k2 * d[k].x;
	      dcos_k.y = tmp_jk * d[j].y - tmp_k2 * d[k].y;
	      dcos_k.z = tmp_jk * d[j].z - tmp_k2 * d[k].z;
	      

	      ddcos_jj.xx = - 2 * tmp_jk * tmp_jj * d[j].x * d[k].x 
		            + 3 * tmp_j2 * tmp_jj * d[j].x * d[j].x - tmp_j2;

	      ddcos_jj.xy = - tmp_jk * tmp_jj * ( d[j].x * d[k].y + d[j].y * d[k].x )
		            + 3 * tmp_j2 * tmp_jj * d[j].x * d[j].y;

	      ddcos_jj.yy = - 2 * tmp_jk * tmp_jj * d[j].y * d[k].y 
		            + 3 * tmp_j2 * tmp_jj * d[j].y * d[j].y - tmp_j2;

	      ddcos_jj.yz = - tmp_jk * tmp_jj * ( d[j].y * d[k].z + d[j].z * d[k].y )
		            + 3 * tmp_j2 * tmp_jj * d[j].y * d[j].z;

	      ddcos_jj.zx = - tmp_jk * tmp_jj * ( d[j].z * d[k].x + d[j].x * d[k].z )
		            + 3 * tmp_j2 * tmp_jj * d[j].z * d[j].x;

	      ddcos_jj.zz = - 2 * tmp_jk * tmp_jj * d[j].z * d[k].z 
		            + 3 * tmp_j2 * tmp_jj * d[j].z * d[j].z - tmp_j2;


	      ddcos_jk.xx = - tmp_jj * tmp_jk * d[j].x * d[j].x 
		            + tmp_k2 * tmp_jj * d[j].x * d[k].x
		            - tmp_kk * tmp_jk * d[k].x * d[k].x + tmp_jk;

	      ddcos_jk.xy = - tmp_jj * tmp_jk * d[j].x * d[j].y 
		            + tmp_k2 * tmp_jj * d[j].x * d[k].y
		            - tmp_kk * tmp_jk * d[k].x * d[k].y;

	      ddcos_jk.xz = - tmp_jj * tmp_jk * d[j].x * d[j].z 
		            + tmp_k2 * tmp_jj * d[j].x * d[k].z
		            - tmp_kk * tmp_jk * d[k].x * d[k].z;

	      ddcos_jk.yx = - tmp_jj * tmp_jk * d[j].y * d[j].x 
		            + tmp_k2 * tmp_jj * d[j].y * d[k].x
		            - tmp_kk * tmp_jk * d[k].y * d[k].x;

	      ddcos_jk.yy = - tmp_jj * tmp_jk * d[j].y * d[j].y 
		            + tmp_k2 * tmp_jj * d[j].y * d[k].y
		            - tmp_kk * tmp_jk * d[k].y * d[k].y + tmp_jk;

	      ddcos_jk.yz = - tmp_jj * tmp_jk * d[j].y * d[j].z 
		              + tmp_k2 * tmp_jj * d[j].y * d[k].z
		              - tmp_kk * tmp_jk * d[k].y * d[k].z;

	      ddcos_jk.zx = - tmp_jj * tmp_jk * d[j].z * d[j].x 
		            + tmp_k2 * tmp_jj * d[j].z * d[k].x
		            - tmp_kk * tmp_jk * d[k].z * d[k].x;

	      ddcos_jk.zy = - tmp_jj * tmp_jk * d[j].z * d[j].y 
		            + tmp_k2 * tmp_jj * d[j].z * d[k].y
		            - tmp_kk * tmp_jk * d[k].z * d[k].y;

	      ddcos_jk.zz = - tmp_jj * tmp_jk * d[j].z * d[j].z 
		            + tmp_k2 * tmp_jj * d[j].z * d[k].z
		            - tmp_kk * tmp_jk * d[k].z * d[k].z + tmp_jk;


	      ddcos_kk.xx = - 2 * tmp_jk * tmp_kk * d[j].x * d[k].x 
		            + 3 * tmp_k2 * tmp_kk * d[k].x * d[k].x - tmp_k2;

	      ddcos_kk.xy = - tmp_jk * tmp_kk * ( d[j].x * d[k].y + d[k].x * d[j].y ) 
		            + 3 * tmp_k2 * tmp_kk * d[k].x * d[k].y;

	      ddcos_kk.yy = - 2 * tmp_jk * tmp_kk * d[j].y * d[k].y 
		            + 3 * tmp_k2 * tmp_kk * d[k].y * d[k].y - tmp_k2;

	      ddcos_kk.yz = - tmp_jk * tmp_kk * ( d[j].y * d[k].z + d[k].y * d[j].z ) 
		            + 3 * tmp_k2 * tmp_kk * d[k].y * d[k].z;

	      ddcos_kk.zx = - tmp_jk * tmp_kk * ( d[j].z * d[k].x + d[k].z * d[j].x ) 
		            + 3 * tmp_k2 * tmp_kk * d[k].z * d[k].x;

	      ddcos_kk.zz = - 2 * tmp_jk * tmp_kk * d[j].z * d[k].z 
		            + 3 * tmp_k2 * tmp_kk * d[k].z * d[k].z - tmp_k2;

	      /* tmp variables for derivatives of zeta */
	      tmp_3     = 2 * ter_c2[p_typ] * tmp_1 * tmp_2 * tmp_2 
		            * fc[k] * ter_om[p_typ][k_typ];
	      tmp_grad  = dfc[k] / r[k] * g_theta * ter_om[p_typ][k_typ];
	      tmp1_zeta = ter_om[p_typ][k_typ] * fc[k] * 2 * ter_c2[p_typ] * tmp_2 * tmp_2;
	      tmp2_zeta = tmp_2 * ( ter_d2[p_typ] - 3 * tmp_1 * tmp_1 );
	      tmp3_zeta = - ter_om[p_typ][k_typ] * dfc[k] / r[k] 
		          * 2 * ter_c2[p_typ] * tmp_1 * tmp_2 * tmp_2;
	      tmp4_zeta = ter_om[p_typ][k_typ] * tmp_kk * ( ddfc[k] - dfc[k] / r[k] ) * g_theta;
	      tmp5_zeta = ter_om[p_typ][k_typ] * dfc[k] / r[k] * g_theta;

	      /* First and second derivatives of zeta */
	      dzeta_bk[k]  = dfc[k] * g_theta;          /* For B' */

	      dzeta_j.x   -= tmp_3 * dcos_j.x;
	      dzeta_j.y   -= tmp_3 * dcos_j.y;
	      dzeta_j.z   -= tmp_3 * dcos_j.z;

	      dzeta_k[k].x = tmp_grad * d[k].x - tmp_3 * dcos_k.x;
	      dzeta_k[k].y = tmp_grad * d[k].y - tmp_3 * dcos_k.y;
	      dzeta_k[k].z = tmp_grad * d[k].z - tmp_3 * dcos_k.z;

	      ddzeta_bkl[k] = ddfc[k] * g_theta; /* For B' */

	      ddzeta_jj.xx += tmp1_zeta * ( tmp2_zeta * dcos_j.x * dcos_j.x
		                                - tmp_1 * ddcos_jj.xx );

	      ddzeta_jj.xy += tmp1_zeta * ( tmp2_zeta * dcos_j.x * dcos_j.y
		                                - tmp_1 * ddcos_jj.xy );

	      ddzeta_jj.yy += tmp1_zeta * ( tmp2_zeta * dcos_j.y * dcos_j.y
		                                - tmp_1 * ddcos_jj.yy );

	      ddzeta_jj.yz += tmp1_zeta * ( tmp2_zeta * dcos_j.y * dcos_j.z
		                                - tmp_1 * ddcos_jj.yz );

	      ddzeta_jj.zx += tmp1_zeta * ( tmp2_zeta * dcos_j.z * dcos_j.x
		                                - tmp_1 * ddcos_jj.zx );

	      ddzeta_jj.zz += tmp1_zeta * ( tmp2_zeta * dcos_j.z * dcos_j.z
		                                - tmp_1 * ddcos_jj.zz );

		
	      ddzeta_jk[k].xx = tmp1_zeta * ( tmp2_zeta * dcos_j.x * dcos_k.x
		                                  - tmp_1 * ddcos_jk.xx )
		                              + tmp3_zeta * d[k].x * dcos_j.x; 

	      ddzeta_jk[k].xy = tmp1_zeta * ( tmp2_zeta * dcos_j.x * dcos_k.y
		                                  - tmp_1 * ddcos_jk.xy )
		                              + tmp3_zeta * d[k].y * dcos_j.x; 

	      ddzeta_jk[k].xz = tmp1_zeta * ( tmp2_zeta * dcos_j.x * dcos_k.z
		                                  - tmp_1 * ddcos_jk.xz )
		                              + tmp3_zeta * d[k].z * dcos_j.x; 

	      ddzeta_jk[k].yx = tmp1_zeta * ( tmp2_zeta * dcos_j.y * dcos_k.x
		                                  - tmp_1 * ddcos_jk.yx )
		                              + tmp3_zeta * d[k].x * dcos_j.y; 

	      ddzeta_jk[k].yy = tmp1_zeta * ( tmp2_zeta * dcos_j.y * dcos_k.y
						  - tmp_1 * ddcos_jk.yy )
		                              + tmp3_zeta * d[k].y * dcos_j.y; 

	      ddzeta_jk[k].yz = tmp1_zeta * ( tmp2_zeta * dcos_j.y * dcos_k.z
		                                  - tmp_1 * ddcos_jk.yz )
		                              + tmp3_zeta * d[k].z * dcos_j.y; 

	      ddzeta_jk[k].zx = tmp1_zeta * ( tmp2_zeta * dcos_j.z * dcos_k.x
		                                  - tmp_1 * ddcos_jk.zx )
		                              + tmp3_zeta * d[k].x * dcos_j.z; 

	      ddzeta_jk[k].zy = tmp1_zeta * ( tmp2_zeta * dcos_j.z * dcos_k.y
		                                  - tmp_1 * ddcos_jk.zy )
		                              + tmp3_zeta * d[k].y * dcos_j.z; 

	      ddzeta_jk[k].zz = tmp1_zeta * ( tmp2_zeta * dcos_j.z * dcos_k.z
		                                  - tmp_1 * ddcos_jk.zz )
		                              + tmp3_zeta * d[k].z * dcos_j.z; 


	      ddzeta_kk[k].xx = + tmp4_zeta * d[k].x * d[k].x
		            + tmp3_zeta * 2 * d[k].x * dcos_k.x
		  + tmp1_zeta * ( tmp2_zeta * dcos_k.x * dcos_k.x - tmp_1 * ddcos_kk.xx )
		  + tmp5_zeta;

	      ddzeta_kk[k].xy = tmp4_zeta * d[k].x * d[k].y 
		            + tmp3_zeta * ( d[k].x * dcos_k.y + d[k].y * dcos_k.x )
		+ tmp1_zeta * ( tmp2_zeta * dcos_k.x * dcos_k.y - tmp_1 * ddcos_kk.xy ) ;

	      ddzeta_kk[k].yy = + tmp4_zeta * d[k].y * d[k].y 
		            + tmp3_zeta * 2 * d[k].y * dcos_k.y
		  + tmp1_zeta * ( tmp2_zeta * dcos_k.y * dcos_k.y - tmp_1 * ddcos_kk.yy )
		  + tmp5_zeta;

	      ddzeta_kk[k].yz = tmp4_zeta * d[k].y * d[k].z 
		            + tmp3_zeta * ( d[k].y * dcos_k.z + d[k].z * dcos_k.y )
		+ tmp1_zeta * ( tmp2_zeta * dcos_k.y * dcos_k.z - tmp_1 * ddcos_kk.yz ) ;

	      ddzeta_kk[k].zx = tmp4_zeta * d[k].z * d[k].x 
		            + tmp3_zeta * ( d[k].z * dcos_k.x + d[k].x * dcos_k.z )
		+ tmp1_zeta * ( tmp2_zeta * dcos_k.z * dcos_k.x - tmp_1 * ddcos_kk.zx ) ;

	      ddzeta_kk[k].zz = + tmp4_zeta * d[k].z * d[k].z 
		            + tmp3_zeta * 2 * d[k].z * dcos_k.z
		  + tmp1_zeta * ( tmp2_zeta * dcos_k.z * dcos_k.z - tmp_1 * ddcos_kk.zz )
		  + tmp5_zeta;

	      /* Third derivative of zeta, for B' */
	      dddzeta_bklm[k] = dddfc[k] * g_theta;

	    } /* k */

	    phi_r  = 0.5 * ter_a[p_typ][j_typ] * exp( - ter_la[p_typ][j_typ] * r[j] );
	    phi_a  = 0.5 * ter_b[p_typ][j_typ] * exp( - ter_mu[p_typ][j_typ] * r[j] );
	    tmp_4  = pow( ters_ga[p_typ] * zeta, ters_n[p_typ] );

	    b_ij  = ter_chi[p_typ][j_typ] * pow( 1 + tmp_4, - 1 / ( 2 * ters_n[p_typ] ) );

	    pot_zwi  = phi_r - b_ij * phi_a;

	    epot += fc[j] * pot_zwi;

	    if ( zeta == 0.0 )   /* only one neighbor of i */
	      tmp_5 = 0.0;
	    else
	      tmp_5 = - b_ij * fc[j] * phi_a * tmp_4 / ( 2 * zeta * ( 1 + tmp_4 ) );

	    tmp_6   = - ( fc[j] * ( - phi_r * ter_la[p_typ][j_typ] 
			+ phi_a * ter_mu[p_typ][j_typ] * b_ij )  + dfc[j] * pot_zwi ) / r[j];

	    /* First derivatives of phi */
	    dphi_r   = - ter_la[p_typ][j_typ] * phi_r;
	    dphi_a   = - ter_mu[p_typ][j_typ] * phi_a;

	    dphi_j.x = tmp_6 * d[j].x + tmp_5 * dzeta_j.x;
	    dphi_j.y = tmp_6 * d[j].y + tmp_5 * dzeta_j.y;
	    dphi_j.z = tmp_6 * d[j].z + tmp_5 * dzeta_j.z;

	    /* tmp variables for derivatives of phi */
	    tmp1_phi = tmp_jj * ( ddfc[j] - dfc[j] / r[j] ) * pot_zwi;
	    tmp2_phi = dfc[j] / r[j] * pot_zwi;
	    tmp3_phi = 2 * dfc[j] * tmp_jj * (- phi_r * ter_la[p_typ][j_typ] 
					      + phi_a * ter_mu[p_typ][j_typ] * b_ij);
	    tmp4_phi = fc[j] * ter_la[p_typ][j_typ] * phi_r / r[j];
	    tmp5_phi = - fc[j] * b_ij * ter_mu[p_typ][j_typ] * phi_a / r[j];
	    tmp6_phi = ter_la[p_typ][j_typ] / r[j] + tmp_jj;
	    tmp7_phi = ter_mu[p_typ][j_typ] / r[j] + tmp_jj;
	    if ( zeta == 0.0 ) {
	      tmp8_phi = 0.0;
	      tmp9_phi = 0.0;
	    }
	    else {
	      tmp8_phi = phi_a / r[j] * ( dfc[j] - fc[j] * ter_mu[p_typ][j_typ] ) 
	               * b_ij * tmp_4 / ( 2.0 * zeta * ( 1.0 + tmp_4 ) );
	      tmp9_phi = - tmp_5 * ( ters_n[p_typ] - 1.0 - 1.5 * tmp_4 ) / ( zeta * ( 1.0 + tmp_4 ) );
	    }

	    /* tmp variables for computation of B' */
	    if ( zeta == 0.0 ) {   /* only one neighbor of i */
	      tmp_b1 = 0.0;
	      tmp_b2 = 0.0;
	      tmp_b3 = 0.0;
	    }
	    else {
	      tmp_b1 = - b_ij * tmp_4 / ( 2.0 * zeta * ( 1.0 + tmp_4 ) );

	      tmp_b2 = tmp_b1 * ( ters_n[p_typ] - 1.0 - 1.5 * tmp_4 ) / ( zeta * ( 1.0 + tmp_4 ) );

	      tmp_b3 = tmp_b1 * ( - 2.0 - ters_n[p_typ] * ters_n[p_typ] + 3.0 * ters_n[p_typ]  
		 + tmp_4 * ( 9.0 * ters_n[p_typ] - 11 ) / 2.0 
		 + tmp_4 * tmp_4 * ( 2.0 * ters_n[p_typ] - 15.0 ) / 4.0 ) 
	      / (zeta * zeta * ( 1.0 + tmp_4 ) * ( 1.0 + tmp_4 ) );
	    }

	    /* Second derivatives of phi */
	    ddphi_r = - ter_la[p_typ][j_typ] * dphi_r;
	    ddphi_a = - ter_mu[p_typ][j_typ] * dphi_a;

	    ddphi_jj.xx = tmp1_phi * d[j].x * d[j].x + tmp2_phi
	                + tmp3_phi * d[j].x * d[j].x
	      + tmp4_phi * ( d[j].x * d[j].x * tmp6_phi - 1 )
	      + tmp5_phi * ( d[j].x * d[j].x * tmp7_phi - 1 ) 
	      + 2 * tmp8_phi * dzeta_j.x * d[j].x 
	      + tmp9_phi * dzeta_j.x * dzeta_j.x - tmp_5 * ddzeta_jj.xx;

	    ddphi_jj.xy = tmp1_phi* d[j].x * d[j].y 
	               + tmp3_phi * d[j].x * d[j].y 
	      + ( tmp4_phi * tmp6_phi + tmp5_phi * tmp7_phi ) * d[j].x * d[j].y 
	      + tmp8_phi * ( dzeta_j.x * d[j].y + dzeta_j.y * d[j].x ) 
	      + tmp9_phi * dzeta_j.x * dzeta_j.y - tmp_5 * ddzeta_jj.xy;

	    ddphi_jj.yy = tmp1_phi * d[j].y * d[j].y + tmp2_phi
	                + tmp3_phi * d[j].y * d[j].y
	      + tmp4_phi * ( d[j].y * d[j].y * tmp6_phi - 1 )
	      + tmp5_phi * ( d[j].y * d[j].y * tmp7_phi - 1 ) 
	      + 2 * tmp8_phi * dzeta_j.y * d[j].y 
	      + tmp9_phi * dzeta_j.y * dzeta_j.y - tmp_5 * ddzeta_jj.yy;

	    ddphi_jj.yz = tmp1_phi* d[j].y * d[j].z 
	               + tmp3_phi * d[j].y * d[j].z 
	      + ( tmp4_phi * tmp6_phi + tmp5_phi * tmp7_phi ) * d[j].y * d[j].z 
	      + tmp8_phi * ( dzeta_j.y * d[j].z + dzeta_j.z * d[j].y ) 
	      + tmp9_phi * dzeta_j.y * dzeta_j.z - tmp_5 * ddzeta_jj.yz;

	    ddphi_jj.zx = tmp1_phi* d[j].z * d[j].x 
	               + tmp3_phi * d[j].z * d[j].x 
	      + ( tmp4_phi * tmp6_phi + tmp5_phi * tmp7_phi ) * d[j].z * d[j].x 
	      + tmp8_phi * ( dzeta_j.z * d[j].x + dzeta_j.x * d[j].z ) 
	      + tmp9_phi * dzeta_j.z * dzeta_j.x - tmp_5 * ddzeta_jj.zx;

	    ddphi_jj.zz = tmp1_phi * d[j].z * d[j].z + tmp2_phi
	                + tmp3_phi * d[j].z * d[j].z
	      + tmp4_phi * ( d[j].z * d[j].z * tmp6_phi - 1 )
	      + tmp5_phi * ( d[j].z * d[j].z * tmp7_phi - 1 ) 
	      + 2 * tmp8_phi * dzeta_j.z * d[j].z 
	      + tmp9_phi * dzeta_j.z * dzeta_j.z - tmp_5 * ddzeta_jj.zz;

	    for (k=0; k<neigh->n; ++k) if (k!=j) 
	    {

	      ddphi_jk[k].xx = tmp8_phi * dzeta_k[k].x * d[j].x
		+ tmp9_phi * dzeta_j.x * dzeta_k[k].x - tmp_5 * ddzeta_jk[k].xx;

	      ddphi_jk[k].xy = tmp8_phi * dzeta_k[k].y * d[j].x
		+ tmp9_phi * dzeta_j.x * dzeta_k[k].y - tmp_5 * ddzeta_jk[k].xy;

	      ddphi_jk[k].xz = tmp8_phi * dzeta_k[k].z * d[j].x
		+ tmp9_phi * dzeta_j.x * dzeta_k[k].z - tmp_5 * ddzeta_jk[k].xz;

	      ddphi_jk[k].yx = tmp8_phi * dzeta_k[k].x * d[j].y
		+ tmp9_phi * dzeta_j.y * dzeta_k[k].x - tmp_5 * ddzeta_jk[k].yx;

	      ddphi_jk[k].yy = tmp8_phi * dzeta_k[k].y * d[j].y
		+ tmp9_phi * dzeta_j.y * dzeta_k[k].y - tmp_5 * ddzeta_jk[k].yy;

	      ddphi_jk[k].yz = tmp8_phi * dzeta_k[k].z * d[j].y
		+ tmp9_phi * dzeta_j.y * dzeta_k[k].z - tmp_5 * ddzeta_jk[k].yz;

	      ddphi_jk[k].zx = tmp8_phi * dzeta_k[k].x * d[j].z
		+ tmp9_phi * dzeta_j.z * dzeta_k[k].x - tmp_5 * ddzeta_jk[k].zx;

	      ddphi_jk[k].zy = tmp8_phi * dzeta_k[k].y * d[j].z
		+ tmp9_phi * dzeta_j.z * dzeta_k[k].y - tmp_5 * ddzeta_jk[k].zy;

	      ddphi_jk[k].zz = tmp8_phi * dzeta_k[k].z * d[j].z
		+ tmp9_phi * dzeta_j.z * dzeta_k[k].z - tmp_5 * ddzeta_jk[k].zz;

	      /* First derivative of b_ij , for B' */	      
	      db_k[k] = tmp_b1 * dzeta_bk[k];      

	      for(l=0; l<neigh->n; ++l) if (l!=j) 
	      {

		ddphi_lk I(l,k) .xx = tmp9_phi * dzeta_k[l].x * dzeta_k[k].x;
		ddphi_lk I(l,k) .xy = tmp9_phi * dzeta_k[l].x * dzeta_k[k].y;
		ddphi_lk I(l,k) .yy = tmp9_phi * dzeta_k[l].y * dzeta_k[k].y;
		ddphi_lk I(l,k) .yz = tmp9_phi * dzeta_k[l].y * dzeta_k[k].z;
		ddphi_lk I(l,k) .zx = tmp9_phi * dzeta_k[l].z * dzeta_k[k].x;
		ddphi_lk I(l,k) .zz = tmp9_phi * dzeta_k[l].z * dzeta_k[k].z;

		/* Second derivative of b_ij, for B' */
		ddb_lk I(l,k) = tmp_b2 * dzeta_bk[k] * dzeta_bk[l];
					       
		if ( l==k ) {
 
		  ddphi_lk I(l,k) .xx += - tmp_5 * ddzeta_kk[k].xx; 
		  ddphi_lk I(l,k) .xy += - tmp_5 * ddzeta_kk[k].xy; 
		  ddphi_lk I(l,k) .yy += - tmp_5 * ddzeta_kk[k].yy; 
		  ddphi_lk I(l,k) .yz += - tmp_5 * ddzeta_kk[k].yz; 
		  ddphi_lk I(l,k) .zx += - tmp_5 * ddzeta_kk[k].zx; 
		  ddphi_lk I(l,k) .zz += - tmp_5 * ddzeta_kk[k].zz;

		  /* Second derivative of b_ij, for B' */
		  ddb_lk I(l,k) += - tmp_b1 * ddzeta_bkl[k];
		  
		}

		for(m=0; m<neigh->n; ++m) if (m!=j) 
		{
		  dddb_lkm J(l,k,m) = tmp_b3 * dzeta_bk[l] * dzeta_bk[k] * dzeta_bk[m];

		  if ( k==l )
		    dddb_lkm J(l,k,m) += tmp_b2 * dzeta_bk[m] * ddzeta_bkl[k];
		  if ( k==m )
		    dddb_lkm J(l,k,m) += tmp_b2 * dzeta_bk[l] * ddzeta_bkl[k];
		  if ( l==m ) {
		    dddb_lkm J(l,k,m) += tmp_b2 * dzeta_bk[k] * ddzeta_bkl[l];

		    if ( k==l )
		      dddb_lkm J(l,k,m) += - tmp_b1 * dddzeta_bklm[l];
		  }
		}

	      } /* l */
	    } /* k */

	    /* Third derivatives of phi */
	    dddphi_r = - ter_la[p_typ][j_typ] * ddphi_r;
	    dddphi_a = - ter_mu[p_typ][j_typ] * ddphi_a;

	    /* Compute stress and elastic constants, contribution from j */
	    tmp = 0.5 * d[j].x * dphi_j.x;
	    p->stress[i].xx        -= tmp;
	    jcell->stress[jnum].xx -= tmp;
	    sigma.xx               -= 2 * tmp;
	    tmp = 0.5 * d[j].y * dphi_j.y;
	    p->stress[i].yy        -= tmp;
	    jcell->stress[jnum].yy -= tmp;
	    sigma.yy               -= 2 * tmp;
	    tmp = 0.5 * d[j].z * dphi_j.z;
	    p->stress[i].zz        -= tmp;
	    jcell->stress[jnum].zz -= tmp;
	    sigma.zz               -= 2 * tmp;
	    tmp = 0.25 * ( d[j].y * dphi_j.z + d[j].z * dphi_j.y );
	    p->stress[i].yz        -= tmp;
	    jcell->stress[jnum].yz -= tmp;
	    sigma.yz               -= 2 * tmp;
	    tmp = 0.25 * ( d[j].z * dphi_j.x + d[j].x * dphi_j.z );
	    p->stress[i].zx        -= tmp;
	    jcell->stress[jnum].zx -= tmp;
	    sigma.zx               -= 2 * tmp;
	    tmp = 0.25 * ( d[j].x * dphi_j.y + d[j].y * dphi_j.x );
	    p->stress[i].xy        -= tmp;
	    jcell->stress[jnum].xy -= tmp;
	    sigma.xy               -= 2 * tmp;

	    tmp = 0.5   * ddphi_jj.xx * d[j].x * d[j].x;
	    p->elco[i].c11        += tmp;
	    jcell->elco[jnum].c11 += tmp;
	    c.c11                 += 2 * tmp;
	    tmp = 0.5   * ddphi_jj.xy * d[j].x * d[j].y;
	    p->elco[i].c12        += tmp;
	    jcell->elco[jnum].c12 += tmp;
	    c.c12                 += 2 * tmp;
	    tmp = 0.5   * ddphi_jj.zx * d[j].x * d[j].z;
	    p->elco[i].c13        += tmp;
	    jcell->elco[jnum].c13 += tmp;
	    c.c13                 += 2 * tmp;
	    tmp = 0.5   * ddphi_jj.yy * d[j].y * d[j].y;
	    p->elco[i].c22        += tmp;
	    jcell->elco[jnum].c22 += tmp;
	    c.c22                 += 2 * tmp;
	    tmp = 0.5   * ddphi_jj.yz * d[j].y * d[j].z;
	    p->elco[i].c23        += tmp;
	    jcell->elco[jnum].c23 += tmp;
	    c.c23                 += 2 * tmp;
	    tmp = 0.5   * ddphi_jj.zz * d[j].z * d[j].z;
	    p->elco[i].c33        += tmp;
	    jcell->elco[jnum].c33 += tmp;
	    c.c33                 += 2 * tmp;
	    tmp = 0.125 * ( ddphi_jj.xx * d[j].y * d[j].y 
	                  + ddphi_jj.yy * d[j].x * d[j].x
	              + 2 * ddphi_jj.xy * d[j].x * d[j].y );
	    p->elco[i].c44        += tmp;
	    jcell->elco[jnum].c44 += tmp;
	    c.c44                 += 2 * tmp;
	    tmp = 0.125 * ( ddphi_jj.yz * d[j].z * d[j].x
	                  + ddphi_jj.zz * d[j].y * d[j].x
	                  + ddphi_jj.xy * d[j].z * d[j].z
	                  + ddphi_jj.zx * d[j].y * d[j].z );
	    p->elco[i].c45        += tmp;
	    jcell->elco[jnum].c45 += tmp;
	    c.c45                 += 2 * tmp;
	    tmp = 0.125 * ( ddphi_jj.xy * d[j].y * d[j].z
	                  + ddphi_jj.yy * d[j].x * d[j].z
	                  + ddphi_jj.zx * d[j].y * d[j].y
	                  + ddphi_jj.yz * d[j].x * d[j].y );
	    p->elco[i].c46        += tmp;
	    jcell->elco[jnum].c46 += tmp;
	    c.c46                 += 2 * tmp;
	    tmp = 0.125 * ( ddphi_jj.zz * d[j].x * d[j].x 
	                  + ddphi_jj.xx * d[j].z * d[j].z
	              + 2 * ddphi_jj.zx * d[j].z * d[j].x );
	    p->elco[i].c55        += tmp;
	    jcell->elco[jnum].c55 += tmp;
	    c.c55                 += 2 * tmp;
	    tmp = 0.125 * ( ddphi_jj.zx * d[j].x * d[j].y
	                    + ddphi_jj.xx * d[j].z * d[j].y
	                    + ddphi_jj.yz * d[j].x * d[j].x
	                    + ddphi_jj.xy * d[j].z * d[j].x );
	    p->elco[i].c56        += tmp;
	    jcell->elco[jnum].c56 += tmp;
	    c.c56                 += 2 * tmp;
	    tmp = 0.125 * ( ddphi_jj.xx * d[j].y * d[j].y 
	                  + ddphi_jj.yy * d[j].x * d[j].x
	              + 2 * ddphi_jj.xy * d[j].x * d[j].y );
	    p->elco[i].c66        += tmp;
	    jcell->elco[jnum].c66 += tmp;
	    c.c66                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.xy * d[j].x * d[j].z
	                 + ddphi_jj.zx * d[j].x * d[j].y );
	    p->elco[i].c14        += tmp;
	    jcell->elco[jnum].c14 += tmp;
	    c.c14                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.zx * d[j].x * d[j].x
	                 + ddphi_jj.xx * d[j].z * d[j].x );
	    p->elco[i].c15        += tmp;
	    jcell->elco[jnum].c15 += tmp;
	    c.c15                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.xx * d[j].x * d[j].y
	                 + ddphi_jj.xy * d[j].x * d[j].x );
	    p->elco[i].c16        += tmp;
	    jcell->elco[jnum].c16 += tmp;
	    c.c16                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.yy * d[j].y * d[j].z
	                 + ddphi_jj.yz * d[j].y * d[j].y );
	    p->elco[i].c24        += tmp;
	    jcell->elco[jnum].c24 += tmp;
	    c.c24                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.yz * d[j].y * d[j].x
	                 + ddphi_jj.xy * d[j].y * d[j].z );
	    p->elco[i].c25        += tmp;
	    jcell->elco[jnum].c25 += tmp;
	    c.c25                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.xy * d[j].y * d[j].y
	                 + ddphi_jj.yy * d[j].x * d[j].y );
	    p->elco[i].c26        += tmp;
	    jcell->elco[jnum].c26 += tmp;
	    c.c26                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.yz * d[j].z * d[j].z
	                 + ddphi_jj.zz * d[j].y * d[j].z );
	    p->elco[i].c34        += tmp;
	    jcell->elco[jnum].c34 += tmp;
	    c.c34                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.zz * d[j].z * d[j].x
	                 + ddphi_jj.zx * d[j].z * d[j].z );
	    p->elco[i].c35        += tmp;
	    jcell->elco[jnum].c35 += tmp;
	    c.c35                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.zx * d[j].z * d[j].y
	                    + ddphi_jj.yz * d[j].z * d[j].x );
	    p->elco[i].c36        += tmp;
	    jcell->elco[jnum].c36 += tmp;
	    c.c36                 += 2 * tmp;	    

	    /* B', contribution from j */
	    dbulkm_dp += ( dddfc[j] * ( phi_r - b_ij * phi_a ) 
	      + 3.0 * ddfc[j] * ( dphi_r - b_ij * dphi_a )
	      + 3.0 * dfc[j] * ( ddphi_r - b_ij * ddphi_a )
	      + fc[j] * ( dddphi_r - b_ij * dddphi_a ) ) 
	      * r[j] * r[j] * r[j];

	    /* Compute stress and elastic constants, contributions from k */
	    for (k=0; k<neigh->n; ++k) if (k!=j) 
	    {
	     tmp = 0.5 * d[k].x * tmp_5 * dzeta_k[k].x; 
	     p->stress[i].xx        -= tmp;
	     jcell->stress[jnum].xx -= tmp;
	     sigma.xx               -= 2 * tmp;
	     tmp = 0.5 * d[k].y * tmp_5 * dzeta_k[k].y;
	     p->stress[i].yy        -= tmp;
	     jcell->stress[jnum].yy -= tmp;
	     sigma.yy               -= 2 * tmp;
	     tmp = 0.5 * d[k].z * tmp_5 * dzeta_k[k].z;
	     p->stress[i].zz        -= tmp;
	     jcell->stress[jnum].zz -= tmp;
	     sigma.zz               -= 2 * tmp;
	     tmp = 0.25 * ( d[k].y * tmp_5 * dzeta_k[k].z  
			  + d[k].z * tmp_5 * dzeta_k[k].y );
	     p->stress[i].yz        -= tmp;
	     jcell->stress[jnum].yz -= tmp;
	     sigma.yz               -= 2 * tmp;
	     tmp = 0.25 * ( d[k].z * tmp_5 * dzeta_k[k].x  
			  + d[k].x * tmp_5 * dzeta_k[k].z );
	     p->stress[i].zx        -= tmp;
	     jcell->stress[jnum].zx -= tmp;
	     sigma.zx               -= 2 * tmp;
	     tmp = 0.25 * ( d[k].x * tmp_5 * dzeta_k[k].y  
		          + d[k].y * tmp_5 * dzeta_k[k].x );
	     p->stress[i].xy        -= tmp;
	     jcell->stress[jnum].xy -= tmp;
	     sigma.xy               -= 2 * tmp;

	     tmp = ddphi_jk[k].xx * d[j].x * d[k].x;
	     p->elco[i].c11        += tmp;
	     jcell->elco[jnum].c11 += tmp;
	     c.c11                 += 2 * tmp;
	     tmp = 0.5 * ( ddphi_jk[k].xy * d[j].x * d[k].y 
			   + ddphi_jk[k].yx * d[j].y * d[k].x );
	     p->elco[i].c12        += tmp;
	     jcell->elco[jnum].c12 += tmp;
	     c.c12                 += 2 * tmp;
	     tmp = 0.5 * ( ddphi_jk[k].xz * d[j].x * d[k].z 
			   + ddphi_jk[k].zx * d[j].z * d[k].x );
	     p->elco[i].c13        += tmp;
	     jcell->elco[jnum].c13 += tmp;
	     c.c13                 += 2 * tmp;
	     tmp = ddphi_jk[k].yy * d[j].y * d[k].y;
	     p->elco[i].c22        += tmp;
	     jcell->elco[jnum].c22 += tmp;
	     c.c22                 += 2 * tmp;
	     tmp = 0.5 * ( ddphi_jk[k].yz * d[j].y * d[k].z 
			   + ddphi_jk[k].zy * d[j].z * d[k].y );
	     p->elco[i].c23        += tmp;
	     jcell->elco[jnum].c23 += tmp;
	     c.c23                 += 2 * tmp;
	     tmp = ddphi_jk[k].zz * d[j].z * d[k].z;
	     p->elco[i].c33        += tmp;
	     jcell->elco[jnum].c33 += tmp;
	     c.c33                 += 2 * tmp;
	     tmp = 0.25 * ( ddphi_jk[k].yy * d[j].z * d[k].z
			    + ddphi_jk[k].zz * d[j].y * d[k].y
			    + ddphi_jk[k].yz * d[j].z * d[k].y
			    + ddphi_jk[k].zy * d[j].y * d[k].z );
	     p->elco[i].c44        += tmp;
	     jcell->elco[jnum].c44 += tmp;
	     c.c44                 += 2 * tmp;
	     tmp = 0.125 * ( ddphi_jk[k].yz * d[j].z * d[k].x
			     + ddphi_jk[k].zy * d[j].x * d[k].z
			     + ddphi_jk[k].zz * ( d[j].y * d[k].x + d[j].x * d[j].y )
			     + 2 * ddphi_jk[k].yx * d[j].z * d[k].z
			     + ddphi_jk[k].zx * d[j].y * d[k].z
			     + ddphi_jk[k].xz * d[j].z * d[k].y );
	     p->elco[i].c45        += tmp;
	     jcell->elco[jnum].c45 += tmp;
	     c.c45                 += 2 * tmp;	     
	     tmp = 0.125 * ( ddphi_jk[k].xy * d[j].y * d[k].z
			     + ddphi_jk[k].yx * d[j].z * d[k].y
			     + ddphi_jk[k].yy * ( d[j].x * d[k].z + d[j].z * d[j].x )
			     + 2 * ddphi_jk[k].xz * d[j].y * d[k].y
			     + ddphi_jk[k].yz * d[j].x * d[k].y
			     + ddphi_jk[k].zy * d[j].y * d[k].x );
	     p->elco[i].c46        += tmp;
	     jcell->elco[jnum].c46 += tmp;
	     c.c46                 += 2 * tmp;
	     tmp = 0.25 * ( ddphi_jk[k].zz * d[j].x * d[k].x 
			    + ddphi_jk[k].xx * d[j].z * d[k].z
			    + ddphi_jk[k].zx * d[j].x * d[k].z
			    + ddphi_jk[k].xz * d[j].z * d[k].x );
	     p->elco[i].c55        += tmp;
	     jcell->elco[jnum].c55 += tmp;
	     c.c55                 += 2 * tmp;
	     tmp = 0.125 * ( ddphi_jk[k].zx * d[j].x * d[k].y
			     + ddphi_jk[k].xz * d[j].y * d[k].x
			     + ddphi_jk[k].xx * ( d[j].z * d[k].y + d[j].y * d[j].z )
			     + 2 * ddphi_jk[k].zy * d[j].x * d[k].x
			     + ddphi_jk[k].xy * d[j].z * d[k].x
			     + ddphi_jk[k].yx * d[j].x * d[k].z );
	     p->elco[i].c56        += tmp;
	     jcell->elco[jnum].c56 += tmp;
	     c.c56                 += 2 * tmp;
	     tmp = 0.25 * ( ddphi_jk[k].xx * d[j].y * d[k].y 
			    + ddphi_jk[k].yy * d[j].x * d[k].x
			    + ddphi_jk[k].xy * d[j].y * d[k].x
			    + ddphi_jk[k].yx * d[j].x * d[k].y );
	     p->elco[i].c66        += tmp;
	     jcell->elco[jnum].c66 += tmp;
	     c.c66                 += 2 * tmp;
	     tmp = 0.25 * ( ddphi_jk[k].xy * d[j].x * d[k].z
			    + ddphi_jk[k].yx * d[j].z * d[k].x
			    + ddphi_jk[k].xz * d[j].x * d[k].y
			    + ddphi_jk[k].zx * d[j].y * d[k].x );
	     p->elco[i].c14        += tmp;
	     jcell->elco[jnum].c14 += tmp;
	     c.c14                 += 2 * tmp;
	     tmp = 0.25 * ( ( ddphi_jk[k].xz + ddphi_jk[k].zx ) * d[j].x * d[k].x
			    + ddphi_jk[k].xx * ( d[j].x * d[k].z + d[j].z * d[k].x ) );
	     p->elco[i].c15        += tmp;
	     jcell->elco[jnum].c15 += tmp;
	     c.c15                 += 2 * tmp;
	     tmp = 0.25 * ( ddphi_jk[k].xx * ( d[j].x * d[k].y + d[j].y * d[k].x )
			    + ( ddphi_jk[k].xy + ddphi_jk[k].yx ) * d[j].x * d[k].x );
	     p->elco[i].c16        += tmp;
	     jcell->elco[jnum].c16 += tmp;
	     c.c16                 += 2 * tmp;
	     tmp = 0.25 * ( ddphi_jk[k].yy * ( d[j].y * d[k].z + d[j].z * d[k].y )
			    + ( ddphi_jk[k].yz + ddphi_jk[k].zy ) * d[j].y * d[k].y );
	     p->elco[i].c24        += tmp;
	     jcell->elco[jnum].c24 += tmp;
	     c.c24                 += 2 * tmp;
	     tmp = 0.25 * ( ddphi_jk[k].yz * d[j].y * d[k].x
			    + ddphi_jk[k].zy * d[j].x * d[k].y
			    + ddphi_jk[k].yx * d[j].y * d[k].z
			    + ddphi_jk[k].xy * d[j].z * d[k].y );
	     p->elco[i].c25        += tmp;
	     jcell->elco[jnum].c25 += tmp;
	     c.c25                 += 2 * tmp;
	     tmp = 0.25 * ( ( ddphi_jk[k].yx + ddphi_jk[k].xy ) * d[j].y * d[k].y
			    + ddphi_jk[k].yy * ( d[j].y * d[k].x + d[j].x * d[k].y ) );
	     p->elco[i].c26        += tmp;
	     jcell->elco[jnum].c26 += tmp;
	     c.c26                 += 2 * tmp;
	     tmp = 0.25 * ( ( ddphi_jk[k].zy + ddphi_jk[k].yz ) * d[j].z * d[k].z
			    + ddphi_jk[k].zz * ( d[j].z * d[k].y + d[j].y * d[k].z ) );
	     p->elco[i].c34        += tmp;
	     jcell->elco[jnum].c34 += tmp;
	     c.c34                 += 2 * tmp;
	     tmp = 0.25 * ( ddphi_jk[k].zz * ( d[j].z * d[k].x + d[j].x * d[k].z )
			    + ( ddphi_jk[k].zx + ddphi_jk[k].xz ) * d[j].z * d[k].z );
	     p->elco[i].c35        += tmp;
	     jcell->elco[jnum].c35 += tmp;
	     c.c35                 += 2 * tmp;
	     tmp = 0.25 * ( ddphi_jk[k].zx * d[j].z * d[k].y
			    + ddphi_jk[k].xz * d[j].y * d[k].z
			    + ddphi_jk[k].zy * d[j].z * d[k].x
			    + ddphi_jk[k].yz * d[j].x * d[k].z );
	     p->elco[i].c36        += tmp;
	     jcell->elco[jnum].c36 += tmp;
	     c.c36                 += 2 * tmp;

	     /* B', contribution from k */
	     dbulkm_dp -= 3.0 * db_k[k] 
	       * ( ddfc[j] * phi_a  + 2.0 * dfc[j] * dphi_a + fc[j] * ddphi_a )
	       * r[j] * r[j] * r[k];

		for (l=0; l<neigh->n; ++l) if (l!=j) 
		{
		  tmp = 0.5 * ddphi_lk I(l,k) .xx * d[l].x * d[k].x;
		  p->elco[i].c11        += tmp;
		  jcell->elco[jnum].c11 += tmp;
		  c.c11                 += 2 * tmp;
		  tmp = 0.5 * ddphi_lk I(l,k) .xy * d[l].x * d[k].y;
		  p->elco[i].c12        += tmp;
		  jcell->elco[jnum].c12 += tmp;
		  c.c12                 += 2 * tmp;
		  tmp = 0.5 * ddphi_lk I(l,k) .zx * d[l].z * d[k].x;
		  p->elco[i].c13        += tmp;
		  jcell->elco[jnum].c13 += tmp;
		  c.c13                 += 2 * tmp;
		  tmp = 0.5 * ddphi_lk I(l,k) .yy * d[l].y * d[k].y;
		  p->elco[i].c22        += tmp;
		  jcell->elco[jnum].c22 += tmp;
		  c.c22                 += 2 * tmp;
		  tmp = 0.5 * ddphi_lk I(l,k) .yz * d[l].y * d[k].z;
 		  p->elco[i].c23        += tmp;
		  jcell->elco[jnum].c23 += tmp;
		  c.c23                 += 2 * tmp;
		  tmp = 0.5 * ddphi_lk I(l,k) .zz * d[l].z * d[k].z;
 		  p->elco[i].c33        += tmp;
		  jcell->elco[jnum].c33 += tmp;
		  c.c33                 += 2 * tmp;
		  tmp = 0.125 * ( ddphi_lk I(l,k) .xx * d[l].y * d[k].y
		                + ddphi_lk I(l,k) .yy * d[l].x * d[k].x
		            + 2 * ddphi_lk I(l,k) .xy * d[l].y * d[k].x );
 		  p->elco[i].c44        += tmp;
		  jcell->elco[jnum].c44 += tmp;
		  c.c44                 += 2 * tmp;
		  tmp = 0.125 * ( ddphi_lk I(l,k) .yz * d[l].z * d[k].x
		                + ddphi_lk I(l,k) .zz * d[l].y * d[k].x 
		                + ddphi_lk I(l,k) .yx * d[l].z * d[k].z
		                + ddphi_lk I(l,k) .zx * d[l].y * d[k].z );
 		  p->elco[i].c45        += tmp;
		  jcell->elco[jnum].c45 += tmp;
		  c.c45                 += 2 * tmp;		  
		  tmp = 0.125 * ( ddphi_lk I(l,k) .xy * d[l].y * d[k].z
		                + ddphi_lk I(l,k) .yy * d[l].x * d[k].z 
		                + ddphi_lk I(l,k) .xz * d[l].y * d[k].y
		                + ddphi_lk I(l,k) .yz * d[l].x * d[k].y );
 		  p->elco[i].c46        += tmp;
		  jcell->elco[jnum].c46 += tmp;
		  c.c46                 += 2 * tmp;	
		  tmp = 0.125 * ( ddphi_lk I(l,k) .zz * d[l].x * d[k].x
		                + ddphi_lk I(l,k) .xx * d[l].z * d[k].z
		            + 2 * ddphi_lk I(l,k) .zx * d[l].x * d[k].z );
 		  p->elco[i].c55        += tmp;
		  jcell->elco[jnum].c55 += tmp;
		  c.c55                 += 2 * tmp;	
		  tmp = 0.125 * ( ddphi_lk I(l,k) .zx * d[l].x * d[k].y
		                + ddphi_lk I(l,k) .xx * d[l].z * d[k].y 
		                + ddphi_lk I(l,k) .zy * d[l].x * d[k].x
		                + ddphi_lk I(l,k) .xy * d[l].z * d[k].x );
 		  p->elco[i].c56        += tmp;
		  jcell->elco[jnum].c56 += tmp;
		  c.c56                 += 2 * tmp;	
		  tmp = 0.125 * ( ddphi_lk I(l,k) .xx * d[l].y * d[k].y
		                + ddphi_lk I(l,k) .yy * d[l].x * d[k].x
		            + 2 * ddphi_lk I(l,k) .xy * d[l].y * d[k].x );
 		  p->elco[i].c66        += tmp;
		  jcell->elco[jnum].c66 += tmp;
		  c.c66                 += 2 * tmp;			  
		  tmp = 0.25  * ( ddphi_lk I(l,k) .xy * d[l].x * d[k].z
		                 + ddphi_lk I(l,k) .zx * d[l].y * d[k].x );
 		  p->elco[i].c14        += tmp;
		  jcell->elco[jnum].c14 += tmp;
		  c.c14                 += 2 * tmp;	
		  tmp = 0.25 * ( ddphi_lk I(l,k) .zx * d[l].x * d[k].x
		               + ddphi_lk I(l,k) .xx * d[l].x * d[k].z );
 		  p->elco[i].c15        += tmp;
		  jcell->elco[jnum].c15 += tmp;
		  c.c15                 += 2 * tmp;	
		  tmp = 0.25  * ( ddphi_lk I(l,k) .xx * d[l].x * d[k].y
		                + ddphi_lk I(l,k) .xy * d[l].x * d[k].x );
 		  p->elco[i].c16        += tmp;
		  jcell->elco[jnum].c16 += tmp;
		  c.c16                 += 2 * tmp;	
		  tmp = 0.25  * ( ddphi_lk I(l,k) .yy * d[l].y * d[k].z 
		                + ddphi_lk I(l,k) .yz * d[l].y * d[k].y );
		  p->elco[i].c24        += tmp;
		  jcell->elco[jnum].c24 += tmp;
		  c.c24                 += 2 * tmp;	
		  tmp = 0.25  * ( ddphi_lk I(l,k) .yz * d[l].y * d[k].x
		                + ddphi_lk I(l,k) .xy * d[l].z * d[k].y );
		  p->elco[i].c25        += tmp;
		  jcell->elco[jnum].c25 += tmp;
		  c.c25                 += 2 * tmp;	
		  tmp = 0.25  * ( ddphi_lk I(l,k) .xy * d[l].y * d[k].y
		                + ddphi_lk I(l,k) .yy * d[l].y * d[k].x );
		  p->elco[i].c26        += tmp;
		  jcell->elco[jnum].c26 += tmp;
		  c.c26                 += 2 * tmp;	
		  tmp = 0.25  * ( ddphi_lk I(l,k) .zx * d[l].z * d[k].z
		                + ddphi_lk I(l,k) .zz * d[l].z * d[k].y );
		  p->elco[i].c34        += tmp;
		  jcell->elco[jnum].c34 += tmp;
		  c.c34                 += 2 * tmp;	
		  tmp = 0.25  * ( ddphi_lk I(l,k) .zz * d[l].z * d[k].x
		                + ddphi_lk I(l,k) .zx * d[l].z * d[k].z );
		  p->elco[i].c35        += tmp;
		  jcell->elco[jnum].c35 += tmp;
		  c.c35                 += 2 * tmp;	
		  tmp = 0.25  * ( ddphi_lk I(l,k) .zx * d[l].z * d[k].y
		                + ddphi_lk I(l,k) .yz * d[l].x * d[k].z );
		  p->elco[i].c36        += tmp;
		  jcell->elco[jnum].c36 += tmp;
		  c.c36                 += 2 * tmp;	

		  /* B', contribution from l */
		  dbulkm_dp -= 3.0 * ddb_lk I(l,k) 
		    * ( dfc[j] * phi_a + fc[j] * dphi_a  )
		    * r[j] * r[k] * r[l];
		  
		  /* B' contribution from m */
		  for (m=0; m<neigh->n; ++m) if (m!=j) 
		    dbulkm_dp -= fc[j] * dddb_lkm J(l,k,m) * phi_a 
		      * r[k] * r[l] * r[m];

		} /* l */

	    } /* k */
            
	  } /* neighbor j */

	} /* i */

      } /* loop over cells */ 

}

#elif defined(STIWEB)

/******************************************************************************
*
*  do_elco_stiweb -- computes stress tensor and elastic constants using the 
*             Stillinger-Weber potential
*
******************************************************************************/

void do_elco_stiweb(void)
{
  cell   *p;
  int    ic, jc, kc;
  real   *r, *fc, *dfc;
  vektor *d;
  neightab *neigh;
  vektor dcos_j, dcos_k;
  tensor ddcos_jj, ddcos_jk, ddcos_kk;
  vektor dphi_j, dphi_k;
  tensor ddphi_jj, ddphi_jk, ddphi_kk;
  cell   *jcell, *kcell;
  int    i, j, k, l, p_typ, j_typ, k_typ, jnum, knum;
  real   *tmpptr;
  real   tmp_r, tmp_jj;
  real   inv_c, inv_r, c_cut, tmp_frc1, tmp_frc2, tmp_frc;
  real   pot, f_cut;
  real   tmp_phi1, tmp_jk, tmp_kk, tmp_j2, tmp_k2;
  real   tmp_cos, tmp_l, tmp_j, tmp_k, cos_theta;
  real   phi_r, phi_a;
  real   tmp, tmp2_j, tmp2_k, tmp2_l, bulkm = 0.0;
  real   tmp1_b, tmp2_b, ddphi_jj_rr;
  real   tmp3_b, tmp_b_j, tmp_b_k, ddphi_jj_b, ddphi_kk_b, ddphi_jk_b;
  real   tmp_p, tmp_q, dddphi_jjj_rrr;
  real   tmp_aj, tmp_ak, dddphi_jjj, dddphi_jjk, dddphi_kkj, dddphi_kkk;

  d    = (vektor *) malloc( neigh_len * sizeof(vektor) );
  r    = (real *)   malloc( neigh_len * sizeof(real)   );
  fc   = (real *)   malloc( neigh_len * sizeof(real)   );
  dfc  = (real *)   malloc( neigh_len * sizeof(real)   );

  if ((d==NULL) || (r==NULL) || (fc==NULL) || (dfc==NULL))
    error("cannot allocate memory for temporary neighbor data");

  /* Initializations */
  for (ic=0; ic < cell_dim.x; ++ic)
    for (jc=0; jc < cell_dim.y; ++jc)
      for (kc=0; kc < cell_dim.z; ++kc)
      {
        p = PTR_3D_V(cell_array,ic,jc,kc,cell_dim);
	
	for (i=0; i<p->n; ++i) 
	{
	  p->elco[i].c11 = 0.0;
	  p->elco[i].c12 = 0.0;
	  p->elco[i].c13 = 0.0;
	  p->elco[i].c14 = 0.0;
	  p->elco[i].c15 = 0.0;
	  p->elco[i].c16 = 0.0;
	  p->elco[i].c22 = 0.0;
	  p->elco[i].c23 = 0.0;
	  p->elco[i].c24 = 0.0;
	  p->elco[i].c25 = 0.0;
	  p->elco[i].c26 = 0.0;
	  p->elco[i].c33 = 0.0;
	  p->elco[i].c34 = 0.0;
	  p->elco[i].c35 = 0.0;
	  p->elco[i].c36 = 0.0;
	  p->elco[i].c44 = 0.0;
	  p->elco[i].c45 = 0.0;
	  p->elco[i].c46 = 0.0;
	  p->elco[i].c55 = 0.0;
	  p->elco[i].c56 = 0.0;
	  p->elco[i].c66 = 0.0;

	  p->stress[i].xx = 0.0;
	  p->stress[i].xy = 0.0;
	  p->stress[i].xz = 0.0;
	  p->stress[i].yx = 0.0;
	  p->stress[i].yy = 0.0;
	  p->stress[i].yz = 0.0;
	  p->stress[i].zx = 0.0;
	  p->stress[i].zy = 0.0;
	  p->stress[i].zz = 0.0;

	}
      }

  /* For each cell */
  for (ic=0; ic < cell_dim.x; ++ic)
    for (jc=0; jc < cell_dim.y; ++jc)
      for (kc=0; kc < cell_dim.z; ++kc)
      {
        p = PTR_3D_V(cell_array,ic,jc,kc,cell_dim);

	/* For each atom in cell */
	for (i=0; i<p->n; ++i) 
	{

	  p_typ   = p->sorte[i];
	  neigh   = &p->neightab_array[i];

	  /* Construct some data for all neighbors */
	  tmpptr = neigh->dist;
	  for (j=0; j<neigh->n; ++j) 
	  {

	    /* Type, distance vector, radii */
	    j_typ   = neigh->typ[j];
	    d[j].x  = *tmpptr++;
	    d[j].y  = *tmpptr++;
	    d[j].z  = *tmpptr++;
	    r[j]    = sqrt(SPROD(d[j],d[j]));

	    /* cutoff function for three-body term */
	    tmp_r  = 1.0 / ( r[j] - sw_a2[p_typ][j_typ] );
	    fc[j]  = exp( sw_ga[p_typ][j_typ] * tmp_r );
	    dfc[j] = - sw_ga[p_typ][j_typ] * tmp_r * tmp_r;
	     
	  } /* j */

	  /*************************************************************/

	  /* For each neighbor of i */
	  for (j=0; j<neigh->n; ++j)
	  {

	    j_typ = neigh->typ[j];
	    jcell = (cell *) neigh->cl [j];
	    jnum  = neigh->num[j];
	    tmp_jj = 1 / ( r[j] * r[j] );

	    /* Potential energy of pair potential part */
	    phi_r  =   0.5 * sw_a[p_typ][j_typ] * pow( r[j], - sw_p[p_typ][j_typ] );
	    phi_a  = - 0.5 * sw_b[p_typ][j_typ] * pow( r[j], - sw_q[p_typ][j_typ] );
	    inv_c  = 1.0 / ( r[j] - sw_a1[p_typ][j_typ] );
	    inv_r  = 1.0 / r[j];
	    f_cut  = exp( sw_de[p_typ][j_typ] * inv_c );

	    pot  = ( phi_r + phi_a ) * f_cut;

	    epot += pot;

	    /* First derivative of pair potential part */
	    tmp_frc1 = pot * sw_de[p_typ][j_typ] * inv_c * inv_c;
	    tmp_frc2 = f_cut * inv_r * ( sw_p[p_typ][j_typ] * phi_r 
					 + sw_q[p_typ][j_typ] * phi_a );
	    tmp_frc = - inv_r * ( tmp_frc1 + tmp_frc2 );

	    dphi_j.x = tmp_frc * d[j].x;
	    dphi_j.y = tmp_frc * d[j].y;
	    dphi_j.z = tmp_frc * d[j].z;

	    /* Compute stress for pair potential part */
	    tmp = 0.5 * dphi_j.x * d[j].x;
	    p->stress[i].xx        += tmp;
	    jcell->stress[jnum].xx += tmp;
	    sigma.xx               += 2 * tmp;
	    tmp = 0.5 * dphi_j.y * d[j].y;
	    p->stress[i].yy        += tmp;
	    jcell->stress[jnum].yy += tmp;
	    sigma.yy               += 2 * tmp;
	    tmp = 0.5 * dphi_j.z * d[j].z;
	    p->stress[i].zz        += tmp;
	    jcell->stress[jnum].zz += tmp;
	    sigma.zz               += 2 * tmp;
	    tmp = 0.5 * dphi_j.y * d[j].z;
	    p->stress[i].yz        += tmp;
	    jcell->stress[jnum].yz += tmp;
	    sigma.yz               += 2 * tmp;
	    tmp = 0.5 * dphi_j.z * d[j].x;
	    p->stress[i].zx        += tmp;
	    jcell->stress[jnum].zx += tmp;
	    sigma.zx               += 2 * tmp;
	    tmp = 0.5 * dphi_j.x * d[j].y;
	    p->stress[i].xy        += tmp;
	    jcell->stress[jnum].xy += tmp;
	    sigma.xy               += 2 * tmp;

	    /* Second derivatives of pair potential part */
	    tmp_phi1 = ( ( tmp_frc1 + 2 * tmp_frc2 ) 
	      * ( sw_de[p_typ][j_typ] * inv_c * inv_c + inv_r ) 
	        + 2 * tmp_frc1 * inv_c + f_cut * inv_r * inv_r 
	      * ( sw_p[p_typ][j_typ] * sw_p[p_typ][j_typ] * phi_r 
		+ sw_q[p_typ][j_typ] * sw_q[p_typ][j_typ] * phi_a ) ) * inv_r * inv_r;

	    ddphi_jj.xx = tmp_phi1 * d[j].x * d[j].x + tmp_frc;
	    ddphi_jj.xy = tmp_phi1 * d[j].x * d[j].y;
	    ddphi_jj.yy = tmp_phi1 * d[j].y * d[j].y + tmp_frc;
	    ddphi_jj.yz = tmp_phi1 * d[j].y * d[j].z;
	    ddphi_jj.zx = tmp_phi1 * d[j].z * d[j].x;
	    ddphi_jj.zz = tmp_phi1 * d[j].z * d[j].z + tmp_frc;

	    /* For computation of Bulk modulus */
	    tmp1_b = sw_de[p_typ][j_typ] * r[j] * inv_c * inv_c;
	    tmp2_b = 2.0 * tmp1_b * r[j] * inv_c;
	    tmp_p = sw_p[p_typ][j_typ];
	    tmp_q = sw_q[p_typ][j_typ];

	    ddphi_jj_rr = ( tmp2_b + tmp_p 
	      + ( tmp1_b + tmp_p ) * ( tmp1_b + tmp_p ) ) * phi_r * f_cut
	     + ( tmp2_b + tmp_q 
	      + ( tmp1_b + tmp_q ) * ( tmp1_b + tmp_q ) ) * phi_a * f_cut;

	    bulkm += ddphi_jj_rr;

	    /* For computation of B' */
	    tmp3_b = 3.0 * tmp2_b * r[j] * inv_c;
	    
	    dddphi_jjj_rrr = - ( tmp3_b + 2.0 * tmp_p 
	      + 3.0 * ( tmp1_b + tmp_p ) * ( tmp2_b + tmp_p )
	      + ( tmp1_b + tmp_p ) * ( tmp1_b + tmp_p ) * ( tmp1_b + tmp_p ) ) * phi_r * f_cut
	     - ( tmp3_b + 2.0 * tmp_q 
	      + 3.0 * ( tmp1_b + tmp_q ) * ( tmp2_b + tmp_q )
	      + ( tmp1_b + tmp_q ) * ( tmp1_b + tmp_q ) * ( tmp1_b + tmp_q ) ) * phi_a * f_cut;

	    dbulkm_dp += dddphi_jjj_rrr;

	    /* Compute elastic constants for pair potential part */
	    tmp = 0.5   * ddphi_jj.xx * d[j].x * d[j].x;
	    p->elco[i].c11        += tmp;
	    jcell->elco[jnum].c11 += tmp;
	    c.c11                 += 2 * tmp;
	    tmp = 0.5   * ddphi_jj.xy * d[j].x * d[j].y;
	    p->elco[i].c12        += tmp;
	    jcell->elco[jnum].c12 += tmp;
	    c.c12                 += 2 * tmp;
	    tmp = 0.5   * ddphi_jj.zx * d[j].x * d[j].z;
	    p->elco[i].c13        += tmp;
	    jcell->elco[jnum].c13 += tmp;
	    c.c13                 += 2 * tmp;
	    tmp = 0.5   * ddphi_jj.yy * d[j].y * d[j].y;
	    p->elco[i].c22        += tmp;
	    jcell->elco[jnum].c22 += tmp;
	    c.c22                 += 2 * tmp;
	    tmp = 0.5   * ddphi_jj.yz * d[j].y * d[j].z;
	    p->elco[i].c23        += tmp;
	    jcell->elco[jnum].c23 += tmp;
	    c.c23                 += 2 * tmp;
	    tmp = 0.5   * ddphi_jj.zz * d[j].z * d[j].z;
	    p->elco[i].c33        += tmp;
	    jcell->elco[jnum].c33 += tmp;
	    c.c33                 += 2 * tmp;
	    tmp = 0.125 * ( ddphi_jj.yy * d[j].z * d[j].z 
	                  + ddphi_jj.zz * d[j].y * d[j].y
	              + 2 * ddphi_jj.yz * d[j].y * d[j].z );
	    p->elco[i].c44        += tmp;
	    jcell->elco[jnum].c44 += tmp;
	    c.c44                 += 2 * tmp;
	    tmp = 0.125 * ( ddphi_jj.yz * d[j].z * d[j].x
	                  + ddphi_jj.zz * d[j].y * d[j].x
	                  + ddphi_jj.xy * d[j].z * d[j].z
	                  + ddphi_jj.zx * d[j].y * d[j].z );
	    p->elco[i].c45        += tmp;
	    jcell->elco[jnum].c45 += tmp;
	    c.c45                 += 2 * tmp;
	    tmp = 0.125 * ( ddphi_jj.xy * d[j].y * d[j].z
	                  + ddphi_jj.yy * d[j].x * d[j].z
	                  + ddphi_jj.zx * d[j].y * d[j].y
	                  + ddphi_jj.yz * d[j].x * d[j].y );
	    p->elco[i].c46        += tmp;
	    jcell->elco[jnum].c46 += tmp;
	    c.c46                 += 2 * tmp;
	    tmp = 0.125 * ( ddphi_jj.zz * d[j].x * d[j].x 
	                  + ddphi_jj.xx * d[j].z * d[j].z
	              + 2 * ddphi_jj.zx * d[j].z * d[j].x );
	    p->elco[i].c55        += tmp;
	    jcell->elco[jnum].c55 += tmp;
	    c.c55                 += 2 * tmp;
	    tmp = 0.125 * ( ddphi_jj.zx * d[j].x * d[j].y
	                    + ddphi_jj.xx * d[j].z * d[j].y
	                    + ddphi_jj.yz * d[j].x * d[j].x
	                    + ddphi_jj.xy * d[j].z * d[j].x );
	    p->elco[i].c56        += tmp;
	    jcell->elco[jnum].c56 += tmp;
	    c.c56                 += 2 * tmp;
	    tmp = 0.125 * ( ddphi_jj.xx * d[j].y * d[j].y 
	                  + ddphi_jj.yy * d[j].x * d[j].x
	              + 2 * ddphi_jj.xy * d[j].x * d[j].y );
	    p->elco[i].c66        += tmp;
	    jcell->elco[jnum].c66 += tmp;
	    c.c66                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.xy * d[j].x * d[j].z
	                 + ddphi_jj.zx * d[j].x * d[j].y );
	    p->elco[i].c14        += tmp;
	    jcell->elco[jnum].c14 += tmp;
	    c.c14                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.zx * d[j].x * d[j].x
	                 + ddphi_jj.xx * d[j].z * d[j].x );
	    p->elco[i].c15        += tmp;
	    jcell->elco[jnum].c15 += tmp;
	    c.c15                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.xx * d[j].x * d[j].y
	                 + ddphi_jj.xy * d[j].x * d[j].x );
	    p->elco[i].c16        += tmp;
	    jcell->elco[jnum].c16 += tmp;
	    c.c16                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.yy * d[j].y * d[j].z
	                 + ddphi_jj.yz * d[j].y * d[j].y );
	    p->elco[i].c24        += tmp;
	    jcell->elco[jnum].c24 += tmp;
	    c.c24                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.yz * d[j].y * d[j].x
	                 + ddphi_jj.xy * d[j].y * d[j].z );
	    p->elco[i].c25        += tmp;
	    jcell->elco[jnum].c25 += tmp;
	    c.c25                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.xy * d[j].y * d[j].y
	                 + ddphi_jj.yy * d[j].x * d[j].y );
	    p->elco[i].c26        += tmp;
	    jcell->elco[jnum].c26 += tmp;
	    c.c26                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.yz * d[j].z * d[j].z
	                 + ddphi_jj.zz * d[j].y * d[j].z );
	    p->elco[i].c34        += tmp;
	    jcell->elco[jnum].c34 += tmp;
	    c.c34                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.zz * d[j].z * d[j].x
	                 + ddphi_jj.zx * d[j].z * d[j].z );
	    p->elco[i].c35        += tmp;
	    jcell->elco[jnum].c35 += tmp;
	    c.c35                 += 2 * tmp;
	    tmp = 0.25 * ( ddphi_jj.zx * d[j].z * d[j].y
	                    + ddphi_jj.yz * d[j].z * d[j].x );
	    p->elco[i].c36        += tmp;
	    jcell->elco[jnum].c36 += tmp;
	    c.c36                 += 2 * tmp;	    

	    for (k=j+1; k<neigh->n; ++k) {

	      k_typ = neigh->typ[k];
	      kcell = (cell *) neigh->cl [k];
	      knum  = neigh->num[k];

	      tmp_jk    = 1 / ( r[j] * r[k] );  
	      cos_theta = SPROD(d[j],d[k]) * tmp_jk;
	      tmp_cos = cos_theta + 1.0 / 3.0;
	      tmp_j2 = cos_theta * tmp_jj;
	      tmp_kk = 1 /  ( r[k] * r[k] );
	      tmp_k2 = cos_theta * tmp_kk;

	      /* Potential energy of three-body part */
	      epot += sw_la[p_typ][j_typ][k_typ] * fc[j] * fc[k] * tmp_cos * tmp_cos;

	      /* First derivatives of cos(theta) */
	      dcos_j.x = tmp_jk * d[k].x - tmp_j2 * d[j].x;
	      dcos_j.y = tmp_jk * d[k].y - tmp_j2 * d[j].y;
	      dcos_j.z = tmp_jk * d[k].z - tmp_j2 * d[j].z;
	      
	      dcos_k.x = tmp_jk * d[j].x - tmp_k2 * d[k].x;
	      dcos_k.y = tmp_jk * d[j].y - tmp_k2 * d[k].y;
	      dcos_k.z = tmp_jk * d[j].z - tmp_k2 * d[k].z;
	
	      /* First derivatives of three-body potential part */	      
	      tmp_l = sw_la[p_typ][j_typ][k_typ] * fc[j] * fc[k] * tmp_cos;
	      tmp_j = tmp_cos * dfc[j] / r[j];
	      tmp_k = tmp_cos * dfc[k] / r[k];

	      dphi_j.x = tmp_l * ( tmp_j * d[j].x + 2.0 * dcos_j.x );
	      dphi_j.y = tmp_l * ( tmp_j * d[j].y + 2.0 * dcos_j.y );
	      dphi_j.z = tmp_l * ( tmp_j * d[j].z + 2.0 * dcos_j.z );

	      dphi_k.x = tmp_l * ( tmp_k * d[k].x + 2.0 * dcos_k.x );
	      dphi_k.y = tmp_l * ( tmp_k * d[k].y + 2.0 * dcos_k.y );
	      dphi_k.z = tmp_l * ( tmp_k * d[k].z + 2.0 * dcos_k.z );     

	      /* Compute stress for three-body potential part */	      
	      tmp = 0.5 * ( dphi_j.x * d[j].x + dphi_k.x * d[k].x );
	      p->stress[i].xx        += tmp;
	      jcell->stress[jnum].xx += 0.5 * tmp;
	      kcell->stress[knum].xx += 0.5 * tmp;
	      sigma.xx               += 2 * tmp;
	      tmp = 0.5 * ( dphi_j.y * d[j].y + dphi_k.y * d[k].y );
	      p->stress[i].yy        += tmp;
	      jcell->stress[jnum].yy += 0.5 * tmp;
	      kcell->stress[knum].yy += 0.5 * tmp;
	      sigma.yy               += 2 * tmp;
	      tmp = 0.5 * ( dphi_j.z * d[j].z + dphi_k.z * d[k].z );
	      p->stress[i].zz        += tmp;
	      jcell->stress[jnum].zz += 0.5 * tmp;
	      kcell->stress[knum].zz += 0.5 * tmp;
	      sigma.zz               += 2 * tmp;
	      tmp = 0.25 * ( dphi_j.y * d[j].z + dphi_k.y * d[k].z
			   + dphi_j.z * d[j].y + dphi_k.z * d[k].y );
	      p->stress[i].yz        += tmp;
	      jcell->stress[jnum].yz += 0.5 * tmp;
	      kcell->stress[knum].yz += 0.5 * tmp;
	      sigma.yz               += 2 * tmp;
	      tmp = 0.25 * ( dphi_j.z * d[j].x + dphi_k.z * d[k].x
			   + dphi_j.x * d[j].z + dphi_k.x * d[k].z );
	      p->stress[i].zx        += tmp;
	      jcell->stress[jnum].zx += 0.5 * tmp;
	      kcell->stress[knum].zx += 0.5 * tmp;
	      sigma.zx               += 2 * tmp;
	      tmp = 0.25 * ( dphi_j.x * d[j].y + dphi_k.x * d[k].y
			   + dphi_j.y * d[j].x + dphi_k.y * d[k].x );
	      p->stress[i].xy        += tmp;
	      jcell->stress[jnum].xy += 0.5 * tmp;
	      kcell->stress[knum].xy += 0.5 * tmp;
	      sigma.xy               += 2 * tmp;
	      
	      /* Second derivatives of cos(theta) */ 
	      ddcos_jj.xx = - 2 * tmp_jk * tmp_jj * d[j].x * d[k].x 
		            + 3 * tmp_j2 * tmp_jj * d[j].x * d[j].x - tmp_j2;

	      ddcos_jj.xy = - tmp_jk * tmp_jj * ( d[j].x * d[k].y + d[j].y * d[k].x )
		            + 3 * tmp_j2 * tmp_jj * d[j].x * d[j].y;

	      ddcos_jj.yy = - 2 * tmp_jk * tmp_jj * d[j].y * d[k].y 
		            + 3 * tmp_j2 * tmp_jj * d[j].y * d[j].y - tmp_j2;

	      ddcos_jj.yz = - tmp_jk * tmp_jj * ( d[j].y * d[k].z + d[j].z * d[k].y )
		            + 3 * tmp_j2 * tmp_jj * d[j].y * d[j].z;

	      ddcos_jj.zx = - tmp_jk * tmp_jj * ( d[j].z * d[k].x + d[j].x * d[k].z )
		            + 3 * tmp_j2 * tmp_jj * d[j].z * d[j].x;

	      ddcos_jj.zz = - 2 * tmp_jk * tmp_jj * d[j].z * d[k].z 
		            + 3 * tmp_j2 * tmp_jj * d[j].z * d[j].z - tmp_j2;


	      ddcos_jk.xx = - tmp_jj * tmp_jk * d[j].x * d[j].x 
		            + tmp_k2 * tmp_jj * d[j].x * d[k].x
		            - tmp_kk * tmp_jk * d[k].x * d[k].x + tmp_jk;

	      ddcos_jk.xy = - tmp_jj * tmp_jk * d[j].x * d[j].y 
		            + tmp_k2 * tmp_jj * d[j].x * d[k].y
		            - tmp_kk * tmp_jk * d[k].x * d[k].y;

	      ddcos_jk.xz = - tmp_jj * tmp_jk * d[j].x * d[j].z 
		            + tmp_k2 * tmp_jj * d[j].x * d[k].z
		            - tmp_kk * tmp_jk * d[k].x * d[k].z;

	      ddcos_jk.yx = - tmp_jj * tmp_jk * d[j].y * d[j].x 
		            + tmp_k2 * tmp_jj * d[j].y * d[k].x
		            - tmp_kk * tmp_jk * d[k].y * d[k].x;

	      ddcos_jk.yy = - tmp_jj * tmp_jk * d[j].y * d[j].y 
		            + tmp_k2 * tmp_jj * d[j].y * d[k].y
		            - tmp_kk * tmp_jk * d[k].y * d[k].y + tmp_jk;

	      ddcos_jk.yz = - tmp_jj * tmp_jk * d[j].y * d[j].z 
		            + tmp_k2 * tmp_jj * d[j].y * d[k].z
		            - tmp_kk * tmp_jk * d[k].y * d[k].z;

	      ddcos_jk.zx = - tmp_jj * tmp_jk * d[j].z * d[j].x 
		            + tmp_k2 * tmp_jj * d[j].z * d[k].x
		            - tmp_kk * tmp_jk * d[k].z * d[k].x;

	      ddcos_jk.zy = - tmp_jj * tmp_jk * d[j].z * d[j].y 
		            + tmp_k2 * tmp_jj * d[j].z * d[k].y
		            - tmp_kk * tmp_jk * d[k].z * d[k].y;

	      ddcos_jk.zz = - tmp_jj * tmp_jk * d[j].z * d[j].z 
		            + tmp_k2 * tmp_jj * d[j].z * d[k].z
		            - tmp_kk * tmp_jk * d[k].z * d[k].z + tmp_jk;


	      ddcos_kk.xx = - 2 * tmp_jk * tmp_kk * d[j].x * d[k].x 
		            + 3 * tmp_k2 * tmp_kk * d[k].x * d[k].x - tmp_k2;

	      ddcos_kk.xy = - tmp_jk * tmp_kk * ( d[j].x * d[k].y + d[k].x * d[j].y ) 
		            + 3 * tmp_k2 * tmp_kk * d[k].x * d[k].y;

	      ddcos_kk.yy = - 2 * tmp_jk * tmp_kk * d[j].y * d[k].y 
		            + 3 * tmp_k2 * tmp_kk * d[k].y * d[k].y - tmp_k2;

	      ddcos_kk.yz = - tmp_jk * tmp_kk * ( d[j].y * d[k].z + d[k].y * d[j].z ) 
		            + 3 * tmp_k2 * tmp_kk * d[k].y * d[k].z;

	      ddcos_kk.zx = - tmp_jk * tmp_kk * ( d[j].z * d[k].x + d[k].z * d[j].x ) 
		            + 3 * tmp_k2 * tmp_kk * d[k].z * d[k].x;

	      ddcos_kk.zz = - 2 * tmp_jk * tmp_kk * d[j].z * d[k].z 
		            + 3 * tmp_k2 * tmp_kk * d[k].z * d[k].z - tmp_k2;

	      /* Second derivatives of three-body potential part */
	      tmp_j = dfc[j] / r[j];
	      tmp2_j = - ( 3 * r[j] - sw_a2[p_typ][j_typ] ) * tmp_jj / ( r[j] - sw_a2[p_typ][j_typ] );
	      tmp_k = dfc[k] / r[k];
	      tmp2_k = - ( 3 * r[k] - sw_a2[p_typ][k_typ] ) * tmp_kk / ( r[k] - sw_a2[p_typ][k_typ] );
	      tmp2_l = 2 * sw_la[p_typ][j_typ][k_typ] * fc[j] * fc[k];

	      ddphi_jj.xx = tmp_j * ( dphi_j.x * d[j].x   
	        + tmp_l * ( tmp2_j * d[j].x * d[j].x + 4 * d[j].x * dcos_j.x + 1.0 ) )
	        + 2 * tmp_l * ddcos_jj.xx
	       	+ tmp2_l * dcos_j.x * dcos_j.x;
 
	      ddphi_jj.xy = tmp_j * ( dphi_j.y * d[j].x  
	        + tmp_l * ( tmp2_j * d[j].x * d[j].y + 2 * ( d[j].x * dcos_j.y + d[j].y * dcos_j.x ) ) )
	        + 2 * tmp_l * ddcos_jj.xy
	       	+ tmp2_l * dcos_j.x * dcos_j.y;

	      ddphi_jj.yy = tmp_j * ( dphi_j.y * d[j].y   
	        + tmp_l * ( tmp2_j * d[j].y * d[j].y + 4 * d[j].y * dcos_j.y + 1.0 ) )
	        + 2 * tmp_l * ddcos_jj.yy
	       	+ tmp2_l * dcos_j.y * dcos_j.y;

	      ddphi_jj.yz = tmp_j * ( dphi_j.z * d[j].y  
	        + tmp_l * ( tmp2_j * d[j].y * d[j].z + 2 * ( d[j].y * dcos_j.z + d[j].z * dcos_j.y ) ) )
	        + 2 * tmp_l * ddcos_jj.yz
	       	+ tmp2_l * dcos_j.y * dcos_j.z; 

	      ddphi_jj.zx = tmp_j * ( dphi_j.x * d[j].z  
	        + tmp_l * ( tmp2_j * d[j].z * d[j].x + 2 * ( d[j].z * dcos_j.x + d[j].x * dcos_j.z ) ) )
	        + 2 * tmp_l * ddcos_jj.zx
	       	+ tmp2_l * dcos_j.z * dcos_j.x; 

	      ddphi_jj.zz = tmp_j * ( dphi_j.z * d[j].z   
	        + tmp_l * ( tmp2_j * d[j].z * d[j].z + 4 * d[j].z * dcos_j.z + 1.0 ) )
	        + 2 * tmp_l * ddcos_jj.zz
	       	+ tmp2_l * dcos_j.z * dcos_j.z; 

	      ddphi_jk.xx = tmp_j * d[j].x * dphi_k.x + 2 * tmp_l * ( tmp_k * dcos_j.x * d[k].x + ddcos_jk.xx )
		+ tmp2_l * dcos_j.x * dcos_k.x;

	      ddphi_jk.xy = tmp_j * d[j].x * dphi_k.y + 2 * tmp_l * ( tmp_k * dcos_j.x * d[k].y + ddcos_jk.xy )
		+ tmp2_l * dcos_j.x * dcos_k.y;

	      ddphi_jk.xz = tmp_j * d[j].x * dphi_k.z + 2 * tmp_l * ( tmp_k * dcos_j.x * d[k].z + ddcos_jk.xz )
		+ tmp2_l * dcos_j.x * dcos_k.z;

	      ddphi_jk.yx = tmp_j * d[j].y * dphi_k.x + 2 * tmp_l * ( tmp_k * dcos_j.y * d[k].x + ddcos_jk.yx )
		+ tmp2_l * dcos_j.y * dcos_k.x;

	      ddphi_jk.yy = tmp_j * d[j].y * dphi_k.y + 2 * tmp_l * ( tmp_k * dcos_j.y * d[k].y + ddcos_jk.yy )
		+ tmp2_l * dcos_j.y * dcos_k.y;

	      ddphi_jk.yz = tmp_j * d[j].y * dphi_k.z + 2 * tmp_l * ( tmp_k * dcos_j.y * d[k].z + ddcos_jk.yz )
		+ tmp2_l * dcos_j.y * dcos_k.z;

	      ddphi_jk.zx = tmp_j * d[j].z * dphi_k.x + 2 * tmp_l * ( tmp_k * dcos_j.z * d[k].x + ddcos_jk.zx )
		+ tmp2_l * dcos_j.z * dcos_k.x;

	      ddphi_jk.zy = tmp_j * d[j].z * dphi_k.y + 2 * tmp_l * ( tmp_k * dcos_j.z * d[k].y + ddcos_jk.zy )
		+ tmp2_l * dcos_j.z * dcos_k.y;

	      ddphi_jk.zz = tmp_j * d[j].z * dphi_k.z + 2 * tmp_l * ( tmp_k * dcos_j.z * d[k].z + ddcos_jk.zz )
		+ tmp2_l * dcos_j.z * dcos_k.z;

	      ddphi_kk.xx = tmp_k * ( d[k].x * dphi_k.x   
	        + tmp_l * ( tmp2_k * d[k].x * d[k].x + 4 * d[k].x * dcos_k.x + 1.0 ) )
	        + 2 * tmp_l * ddcos_kk.xx
	       	+ tmp2_l * dcos_k.x * dcos_k.x;

	      ddphi_kk.xy = tmp_k * ( d[k].x * dphi_k.y   
	        + tmp_l * ( tmp2_k * d[k].x * d[k].y + 2 * ( d[k].x * dcos_k.y + d[k].y * dcos_k.x ) ) )
	        + 2 * tmp_l * ddcos_kk.xy
	       	+ tmp2_l * dcos_k.x * dcos_k.y;

	      ddphi_kk.yy = tmp_k * ( d[k].y * dphi_k.y   
	        + tmp_l * ( tmp2_k * d[k].y * d[k].y + 4 * d[k].y * dcos_k.y + 1.0 ) )
	        + 2 * tmp_l * ddcos_kk.yy
	       	+ tmp2_l * dcos_k.y * dcos_k.y; 

	      ddphi_kk.yz = tmp_k * ( d[k].y * dphi_k.z   
	        + tmp_l * ( tmp2_k * d[k].y * d[k].z + 2 * ( d[k].y * dcos_k.z + d[k].z * dcos_k.y ) ) )
	        + 2 * tmp_l * ddcos_kk.yz
	       	+ tmp2_l * dcos_k.y * dcos_k.z;

	      ddphi_kk.zx = tmp_k * ( d[k].z * dphi_k.x   
	        + tmp_l * ( tmp2_k * d[k].z * d[k].x + 2 * ( d[k].z * dcos_k.x + d[k].x * dcos_k.z ) ) )
	        + 2 * tmp_l * ddcos_kk.zx
	       	+ tmp2_l * dcos_k.z * dcos_k.x; 

	      ddphi_kk.zz = tmp_k * ( d[k].z * dphi_k.z   
	        + tmp_l * ( tmp2_k * d[k].z * d[k].z + 4 * d[k].z * dcos_k.z + 1.0 ) )
	        + 2 * tmp_l * ddcos_kk.zz
	       	+ tmp2_l * dcos_k.z * dcos_k.z; 

	      /* For computation of bulk modulus */
	      tmp3_b = sw_la[p_typ][j_typ][k_typ] * fc[j] * fc[k] * tmp_cos * tmp_cos;
	      tmp_b_j = sw_ga[p_typ][j_typ] / ( r[j] - sw_a2[p_typ][j_typ] );
	      tmp_b_k = sw_ga[p_typ][k_typ] / ( r[k] - sw_a2[p_typ][k_typ] );
	      tmp_aj = sw_a2[p_typ][j_typ];
	      tmp_ak = sw_a2[p_typ][k_typ];

	      ddphi_jj_b = tmp_b_j / ( ( r[j] - tmp_aj ) * ( r[j] - tmp_aj ) )
		* ( 2.0 + tmp_b_j ) * tmp3_b;
	      ddphi_kk_b = tmp_b_k / ( r[k] - tmp_ak ) / ( r[k] - tmp_ak )
		* ( 2.0 + tmp_b_k ) * tmp3_b;
	      ddphi_jk_b = tmp_b_j / ( r[j] - tmp_aj )
		* tmp_b_k / ( r[k] - tmp_ak ) * tmp3_b;

	      bulkm += ddphi_jj_b * r[j] * r[j] 
		+ 2 * ddphi_jk_b * r[j] * r[k] 
		+ ddphi_kk_b * r[k] * r[k];

	      /* For computation of B' */	      
	      dddphi_jjj = - tmp_b_j / ( ( r[j] - tmp_aj ) * ( r[j] - tmp_aj ) * ( r[j] - tmp_aj ) )
		* ( 6.0 + tmp_b_j * ( 6.0 + tmp_b_j ) ) * tmp3_b;

	      dddphi_kkk = - tmp_b_k / ( ( r[k] - tmp_ak ) * ( r[k] - tmp_ak ) * ( r[k] - tmp_ak ) )
		* ( 6.0 + tmp_b_k * ( 6.0 + tmp_b_k ) ) * tmp3_b;

	      dddphi_jjk = - tmp_b_j * tmp_b_k / ( ( r[j] - tmp_aj ) * ( r[j] - tmp_aj ) * ( r[k] - tmp_ak ) )
		* ( tmp_b_j + 2.0 ) * tmp3_b;

	      dddphi_kkj = - tmp_b_j * tmp_b_k / ( ( r[j] - tmp_aj ) * ( r[k] - tmp_ak ) * ( r[k] - tmp_ak ) )
		* ( tmp_b_k + 2.0 ) * tmp3_b;

	      dbulkm_dp += dddphi_jjj * r[j] * r[j] * r[j]
                + 3.0 * dddphi_jjk * r[j] * r[j] * r[k]
                + 3.0 * dddphi_kkj * r[k] * r[k] * r[j]
		+  dddphi_kkk * r[k] * r[k] * r[k];

	      /* Compute elastic constants of the three-body potential part */
	      tmp = 0.5 * 
		      ( ddphi_jj.xx * d[j].x * d[j].x 
		  + 2 * ddphi_jk.xx * d[j].x * d[k].x 
		      + ddphi_kk.xx * d[k].x * d[k].x );
	      p->elco[i].c11        += tmp;
	      jcell->elco[jnum].c11 += 0.5 * tmp;
	      kcell->elco[knum].c11 += 0.5 * tmp;	      
	      c.c11                 += 2 * tmp;
	      tmp = 0.5 * 
		      ( ddphi_jj.xy * d[j].x * d[j].y 
		      + ddphi_jk.xy * d[j].x * d[k].y
		      + ddphi_jk.yx * d[k].x * d[j].y
		      + ddphi_kk.xy * d[k].x * d[k].y ); 
	      p->elco[i].c12        += tmp;
	      jcell->elco[jnum].c12 += 0.5 * tmp;
	      kcell->elco[knum].c12 += 0.5 * tmp;	      
	      c.c12                 += 2 * tmp;
	      tmp = 0.5 * 
		      ( ddphi_jj.zx * d[j].x * d[j].z 
		      + ddphi_jk.xz * d[j].x * d[k].z
		      + ddphi_jk.zx * d[k].x * d[j].z
		      + ddphi_kk.zx * d[k].x * d[k].z );   
	      p->elco[i].c13        += tmp;
	      jcell->elco[jnum].c13 += 0.5 * tmp;
	      kcell->elco[knum].c13 += 0.5 * tmp;	      
	      c.c13                 += 2 * tmp;
	      tmp = 0.5 * 
		      ( ddphi_jj.yy * d[j].y * d[j].y 
		  + 2 * ddphi_jk.yy * d[j].y * d[k].y 
		      + ddphi_kk.yy * d[k].y * d[k].y );  
	      p->elco[i].c22        += tmp;
	      jcell->elco[jnum].c22 += 0.5 * tmp;
	      kcell->elco[knum].c22 += 0.5 * tmp;	      
	      c.c22                 += 2 * tmp;
	      tmp = 0.5 * 
		      ( ddphi_jj.yz * d[j].y * d[j].z 
		      + ddphi_jk.yz * d[j].y * d[k].z
		      + ddphi_jk.zy * d[k].y * d[j].z
		      + ddphi_kk.yz * d[k].y * d[k].z );
	      p->elco[i].c23        += tmp;
	      jcell->elco[jnum].c23 += 0.5 * tmp;
	      kcell->elco[knum].c23 += 0.5 * tmp;	      
	      c.c23                 += 2 * tmp;
	      tmp = 0.5 * 
		      ( ddphi_jj.zz * d[j].z * d[j].z 
		  + 2 * ddphi_jk.zz * d[j].z * d[k].z 
		      + ddphi_kk.zz * d[k].z * d[k].z ); 
	      p->elco[i].c33        += tmp;
	      jcell->elco[jnum].c33 += 0.5 * tmp;
	      kcell->elco[knum].c33 += 0.5 * tmp;	      
	      c.c33                 += 2 * tmp;
	      tmp = 0.125 * 
		      ( ddphi_jj.yy * d[j].z * d[j].z 
		  + 2 * ddphi_jk.yy * d[j].z * d[k].z 
		      + ddphi_kk.yy * d[k].z * d[k].z 
		  + 2 * ddphi_jj.yz * d[j].y * d[j].z 
		   + 2 * ddphi_jk.yz * d[j].z * d[k].y + 2 * ddphi_jk.zy * d[j].y * d[k].z 
		   + 2 * ddphi_kk.yz * d[k].y * d[k].z 
		  + ddphi_jj.zz * d[j].y * d[j].y 
		   + 2 * ddphi_jk.zz * d[j].y * d[k].y 
		   + ddphi_kk.zz * d[k].y * d[k].y ); 
	      p->elco[i].c44        += tmp;
	      jcell->elco[jnum].c44 += 0.5 * tmp;
	      kcell->elco[knum].c44 += 0.5 * tmp;	      
	      c.c44                 += 2 * tmp;
	      tmp = 0.125 *
		( ddphi_jj.yz * d[j].z * d[j].x 
		  + ddphi_jk.yz * d[j].z * d[k].x
		   + ddphi_jk.zy * d[j].x * d[k].z
		   + ddphi_kk.yz * d[k].z * d[k].x 
		  + ddphi_jj.zz * d[j].y * d[j].x 
		  + ddphi_jk.zz * ( d[j].y * d[k].x + d[j].x * d[k].y ) 
		   + ddphi_kk.zz * d[k].y * d[k].x
		  + ddphi_jj.xy * d[j].z * d[j].z 
		  + ( ddphi_jk.yx + ddphi_jk.xy ) * d[j].z * d[k].z 
		   + ddphi_kk.xy * d[k].z * d[k].z
		  + ddphi_jj.zx * d[j].y * d[j].z 
		   + ddphi_jk.zx * d[j].y * d[k].z
		   + ddphi_jk.xz * d[j].z * d[k].y 
		   + ddphi_kk.zx * d[k].y * d[k].z ); 
	      p->elco[i].c45        += tmp;
	      jcell->elco[jnum].c45 += 0.5 * tmp;
	      kcell->elco[knum].c45 += 0.5 * tmp;	      
	      c.c45                 += 2 * tmp;              
	      tmp = 0.125 *
		( ddphi_jj.xy * d[j].y * d[j].z 
		  + ddphi_jk.xy * d[j].y * d[k].z
		   + ddphi_jk.yx * d[j].z * d[k].y
		   + ddphi_kk.xy * d[k].y * d[k].z 
		  + ddphi_jj.yy * d[j].x * d[j].z 
		  + ddphi_jk.yy * ( d[j].x * d[k].z + d[j].z * d[k].x ) 
		   + ddphi_kk.yy * d[k].x * d[k].z
		  + ddphi_jj.zx * d[j].y * d[j].y 
		  + ( ddphi_jk.xz + ddphi_jk.zx ) * d[j].y * d[k].y 
		   + ddphi_kk.zx * d[k].y * d[k].y
		  + ddphi_jj.yz * d[j].x * d[j].y 
		   + ddphi_jk.yz * d[j].x * d[k].y
		   + ddphi_jk.zy * d[j].y * d[k].x 
		   + ddphi_kk.yz * d[k].x * d[k].y ); 
	      p->elco[i].c46        += tmp;
	      jcell->elco[jnum].c46 += 0.5 * tmp;
	      kcell->elco[knum].c46 += 0.5 * tmp;	      
	      c.c46                 += 2 * tmp;              
	      tmp = 0.125 * 
		      ( ddphi_jj.zz * d[j].x * d[j].x 
		  + 2 * ddphi_jk.zz * d[j].x * d[k].x 
		      + ddphi_kk.zz * d[k].x * d[k].x 
		  + 2 * ddphi_jj.zx * d[j].z * d[j].x 
		   + 2 * ddphi_jk.zx * d[j].x * d[k].z + 2 * ddphi_jk.xz * d[j].z * d[k].x 
		   + 2 * ddphi_kk.zx * d[k].z * d[k].x 
		  + ddphi_jj.xx * d[j].z * d[j].z 
		   + 2 * ddphi_jk.xx * d[j].z * d[k].z 
		   + ddphi_kk.xx * d[k].z * d[k].z ); 
	      p->elco[i].c55        += tmp;
	      jcell->elco[jnum].c55 += 0.5 * tmp;
	      kcell->elco[knum].c55 += 0.5 * tmp;	      
	      c.c55                 += 2 * tmp;
	      tmp = 0.125 *
		( ddphi_jj.zx * d[j].x * d[j].y 
		  + ddphi_jk.zx * d[j].x * d[k].y
		   + ddphi_jk.xz * d[j].y * d[k].x
		   + ddphi_kk.zx * d[k].x * d[k].y 
		  + ddphi_jj.xx * d[j].z * d[j].y 
		  + ddphi_jk.xx * ( d[j].z * d[k].y + d[j].y * d[k].z ) 
		   + ddphi_kk.xx * d[k].z * d[k].y
		  + ddphi_jj.yz * d[j].x * d[j].x 
		  + ( ddphi_jk.zy + ddphi_jk.yz ) * d[j].x * d[k].x 
		   + ddphi_kk.yz * d[k].x * d[k].x
		  + ddphi_jj.xy * d[j].z * d[j].x 
		   + ddphi_jk.xy * d[j].z * d[k].x
		   + ddphi_jk.yx * d[j].x * d[k].z 
		   + ddphi_kk.xy * d[k].z * d[k].x ); 
	      p->elco[i].c56        += tmp;
	      jcell->elco[jnum].c56 += 0.5 * tmp;
	      kcell->elco[knum].c56 += 0.5 * tmp;	      
	      c.c56                 += 2 * tmp;              
	      tmp = 0.125 * 
		      ( ddphi_jj.xx * d[j].y * d[j].y 
		  + 2 * ddphi_jk.xx * d[j].y * d[k].y 
		      + ddphi_kk.xx * d[k].y * d[k].y 
		  + 2 * ddphi_jj.xy * d[j].x * d[j].y 
		   + 2 * ddphi_jk.xy * d[j].y * d[k].x + 2 * ddphi_jk.yx * d[j].x * d[k].y 
		   + 2 * ddphi_kk.xy * d[k].x * d[k].y 
		  + ddphi_jj.yy * d[j].x * d[j].x 
		   + 2 * ddphi_jk.yy * d[j].x * d[k].x 
		   + ddphi_kk.yy * d[k].x * d[k].x ); 
	      p->elco[i].c66        += tmp;
	      jcell->elco[jnum].c66 += 0.5 * tmp;
	      kcell->elco[knum].c66 += 0.5 * tmp;	      
	      c.c66                 += 2 * tmp;
	      tmp = 0.25 * 
		 ( ddphi_jj.xy * d[j].x * d[j].z 
		  + ddphi_jk.xy * d[j].x * d[k].z
 		  + ddphi_jk.yx * d[j].z * d[k].x
		  + ddphi_kk.xy * d[k].x * d[k].z 
		 + ddphi_jj.zx * d[j].x * d[j].y 
		  + ddphi_jk.xz * d[j].x * d[k].y
		  + ddphi_jk.zx * d[j].y * d[k].x 
		  + ddphi_kk.zx * d[k].x * d[k].y ); 
	      p->elco[i].c14        += tmp;
	      jcell->elco[jnum].c14 += 0.5 * tmp;
	      kcell->elco[knum].c14 += 0.5 * tmp;	      
	      c.c14                 += 2 * tmp;
	      tmp = 0.25 * 
		 ( ddphi_jj.zx * d[j].x * d[j].x 
		  + ( ddphi_jk.xz + ddphi_jk.zx ) * d[j].x * d[k].x
		  + ddphi_kk.zx * d[k].x * d[k].x 
		 + ddphi_jj.xx * d[j].x * d[j].z 
	          + ddphi_jk.xx * ( d[j].x * d[k].z + d[j].z * d[k].x )
		  + ddphi_kk.xx * d[k].x * d[k].z ); 
	      p->elco[i].c15        += tmp;
	      jcell->elco[jnum].c15 += 0.5 * tmp;
	      kcell->elco[knum].c15 += 0.5 * tmp;	      
	      c.c15                 += 2 * tmp;
	      tmp = 0.25 * 
		 ( ddphi_jj.xx * d[j].x * d[j].y 
	          + ddphi_jk.xx * ( d[j].x * d[k].y + d[j].y * d[k].x )
		  + ddphi_kk.xx * d[k].x * d[k].y 
		 + ddphi_jj.xy * d[j].x * d[j].x 
	          + ( ddphi_jk.xy + ddphi_jk.yx ) * d[j].x * d[k].x 
		  + ddphi_kk.xy * d[k].x * d[k].x ); 
	      p->elco[i].c16        += tmp;
	      jcell->elco[jnum].c16 += 0.5 * tmp;
	      kcell->elco[knum].c16 += 0.5 * tmp;	      
	      c.c16                 += 2 * tmp;
	      tmp = 0.25 * 
		 ( ddphi_jj.yy * d[j].y * d[j].z 
	          + ddphi_jk.yy * ( d[j].y * d[k].z + d[j].z * d[k].y )
		  + ddphi_kk.yy * d[k].y * d[k].z 
		 + ddphi_jj.yz * d[j].y * d[j].y 
	          + ( ddphi_jk.yz + ddphi_jk.zy ) * d[j].y * d[k].y 
		  + ddphi_kk.yz * d[k].y * d[k].y ); 
	      p->elco[i].c24        += tmp;
	      jcell->elco[jnum].c24 += 0.5 * tmp;
	      kcell->elco[knum].c24 += 0.5 * tmp;	      
	      c.c24                 += 2 * tmp;
	      tmp = 0.25 * 
		 ( ddphi_jj.yz * d[j].y * d[j].x 
		  + ddphi_jk.yz * d[j].y * d[k].x
 		  + ddphi_jk.zy * d[j].x * d[k].y
		  + ddphi_kk.yz * d[k].y * d[k].x 
		 + ddphi_jj.xy * d[j].y * d[j].z 
		  + ddphi_jk.yx * d[j].y * d[k].z
		  + ddphi_jk.xy * d[j].z * d[k].y 
		  + ddphi_kk.xy * d[k].y * d[k].z ); 
	      p->elco[i].c25        += tmp;
	      jcell->elco[jnum].c25 += 0.5 * tmp;
	      kcell->elco[knum].c25 += 0.5 * tmp;	      
	      c.c25                 += 2 * tmp;
	      tmp = 0.25 * 
		 ( ddphi_jj.xy * d[j].y * d[j].y 
		  + ( ddphi_jk.yx + ddphi_jk.xy ) * d[j].y * d[k].y
		  + ddphi_kk.xy * d[k].y * d[k].y
		 + ddphi_jj.yy * d[j].y * d[j].x 
	          + ddphi_jk.yy * ( d[j].y * d[k].x + d[j].x * d[k].y )
		  + ddphi_kk.yy * d[k].y * d[k].x ); 
	      p->elco[i].c26        += tmp;
	      jcell->elco[jnum].c26 += 0.5 * tmp;
	      kcell->elco[knum].c26 += 0.5 * tmp;	      
	      c.c26                 += 2 * tmp;
	      tmp = 0.25 * 
		 ( ddphi_jj.yz * d[j].z * d[j].z 
		  + ( ddphi_jk.zy + ddphi_jk.yz ) * d[j].z * d[k].z
		  + ddphi_kk.yz * d[k].z * d[k].z
		 + ddphi_jj.zz * d[j].z * d[j].y 
	          + ddphi_jk.zz * ( d[j].z * d[k].y + d[j].y * d[k].z )
		  + ddphi_kk.zz * d[k].z * d[k].y ); 
	      p->elco[i].c34        += tmp;
	      jcell->elco[jnum].c34 += 0.5 * tmp;
	      kcell->elco[knum].c34 += 0.5 * tmp;
	      c.c34                 += 2 * tmp;
	      tmp = 0.25 * 
		 ( ddphi_jj.zz * d[j].z * d[j].x 
	          + ddphi_jk.zz * ( d[j].z * d[k].x + d[j].x * d[k].z )
		  + ddphi_kk.zz * d[k].z * d[k].x 
		 + ddphi_jj.zx * d[j].z * d[j].z 
	          + ( ddphi_jk.zx + ddphi_jk.xz ) * d[j].z * d[k].z 
		  + ddphi_kk.zx * d[k].z * d[k].z ); 
	      p->elco[i].c35        += tmp;
	      jcell->elco[jnum].c35 += 0.5 * tmp;
	      kcell->elco[knum].c35 += 0.5 * tmp;	      
	      c.c35                 += 2 * tmp;
	      tmp = 0.25 * 
		 ( ddphi_jj.zx * d[j].z * d[j].y 
		  + ddphi_jk.zx * d[j].z * d[k].y
 		  + ddphi_jk.xz * d[j].y * d[k].z
		  + ddphi_kk.zx * d[k].z * d[k].y 
		 + ddphi_jj.yz * d[j].z * d[j].x 
		  + ddphi_jk.zy * d[j].z * d[k].x
		  + ddphi_jk.yz * d[j].x * d[k].z 
		  + ddphi_kk.yz * d[k].z * d[k].x ); 
	      p->elco[i].c36        += tmp;
	      jcell->elco[jnum].c36 += 0.5 * tmp;
	      kcell->elco[knum].c36 += 0.5 * tmp;	      
	      c.c36                 += 2 * tmp;

	    } /* k */

	  } /* j */

	}

      }

}

#endif
