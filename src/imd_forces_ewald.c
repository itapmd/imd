
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2005 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/*****************************************************************************
*
* imd_forces_ewald -- force loops for Ewald summation
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

void clear_forces(void)
{
  int k;
  for (k=0; k<nallcells; k++) {
    cell *p = cell_array + k;
    int i;
    for (i=0; i<p->n; i++) {
      KRAFT(p,i,X) = 0.0;
      KRAFT(p,i,Y) = 0.0;
      KRAFT(p,i,Z) = 0.0;
    }
  }
}

real force_norm(void)
{
  int k;
  double norm = 0.0;
  for (k=0; k<nallcells; k++) {
    cell *p = cell_array + k;
    int i;
    for (i=0; i<p->n; i++) {
      norm += SQR( KRAFT(p,i,X) );
      norm += SQR( KRAFT(p,i,Y) );
      norm += SQR( KRAFT(p,i,Z) );
    }
  }
  return norm / (DIM*natoms);
}

/******************************************************************************
*
*  do_forces_ewald
*
******************************************************************************/

void do_forces_ewald(int steps)
{
  int n;
  real old1, old2;

  if ((steps==0) && (ew_test)) {
    old1 = tot_pot_energy / natoms;
    printf( "Ewald pair int: Epot = %e, force = %e time = %e\n", 
            old1, force_norm(), ewald_time.total);
    ewald_time.total = 0.0;
    imd_start_timer( &ewald_time );
    clear_forces();
  }

  /* real space part; if ew_nmax < 0, we do it with the other pair forces */
  if (ew_nmax >= 0) do_forces_ewald_real();

  if ((steps==0) && (ew_test)) {
    imd_stop_timer( &ewald_time );
    old2 = tot_pot_energy / natoms;
    printf( "Ewald real space: Epot = %e, force = %e, time = %e\n", 
            old2 - old1, force_norm(), ewald_time.total);
    ewald_time.total = 0.0;
    imd_start_timer( &ewald_time );
    clear_forces();
  }

  /* Fourier space part */
  if (ew_kcut > 0) do_forces_ewald_fourier();

  if ((steps==0) && (ew_test)) {
    imd_stop_timer( &ewald_time );
    old1 = tot_pot_energy / natoms;
    printf( "Ewald Fourier space: Epot = %e, force = %e, time = %e\n", 
            old1 - old2, force_norm(), ewald_time.total);
    ewald_time.total = 0.0;
    imd_start_timer( &ewald_time );
    clear_forces();
  }

  /* self energy part */
  for (n=0; n<ntypes; n++)
    tot_pot_energy -= ew_vorf * num_sort[n] * SQR( charge[n] ) * ew_eps;

  if ((steps==0) && (ew_test)) {
    imd_stop_timer( &ewald_time );
    old2 = tot_pot_energy / natoms;
    printf( "Ewald self-energy: Epot = %e, force = %e, time = %e\n", 
            old2 - old1, force_norm(), ewald_time.total);
    printf( "Ewald total: Epot = %f\n", old2);
    clear_forces();
  }
}

/******************************************************************************
*
*  do_forces_ewald_fourier
*
*  computes the fourier part of the Ewald sum
*
******************************************************************************/

void do_forces_ewald_fourier(void)
{
  int    i, j, k, l, m, n, c;
  int    count, typ;
  real   tmp, tmp_virial=0.0, sum_cos, sum_sin;
  real   kforce, kpot;

  /* Compute exp(ikr) recursively */
  count = 0;
  for (c=0; c<ncells; c++) {

    cell *p = CELLPTR(c);

    for (i=0; i<p->n; i++) {

      tmp = twopi * SPRODX(&ORT(p,i,X),tbox_x);
      coskx I(ew_nx,count,ew_nx+1) =  cos(tmp);
      coskx I(ew_nx,count,ew_nx-1) =  coskx I(ew_nx,count,ew_nx+1);
      sinkx I(ew_nx,count,ew_nx+1) =  sin(tmp);
      sinkx I(ew_nx,count,ew_nx-1) = -sinkx I(ew_nx,count,ew_nx+1); 

      tmp = twopi * SPRODX(&ORT(p,i,X),tbox_y);
      cosky I(ew_ny,count,ew_ny+1) =  cos(tmp);
      cosky I(ew_ny,count,ew_ny-1) =  cosky I(ew_ny,count,ew_ny+1);
      sinky I(ew_ny,count,ew_ny+1) =  sin(tmp);
      sinky I(ew_ny,count,ew_ny-1) = -sinky I(ew_ny,count,ew_ny+1);

      tmp = twopi * SPRODX(&ORT(p,i,X),tbox_z);
      coskz I(ew_nz,count,ew_nz+1) = cos(tmp);
      sinkz I(ew_nz,count,ew_nz+1) = sin(tmp);

      for (j=2; j<=ew_nx; j++) {
        coskx I(ew_nx,count,ew_nx+j) =   coskx I(ew_nx,count,ew_nx+j-1) * coskx I(ew_nx,count,ew_nx+1) 
                                       - sinkx I(ew_nx,count,ew_nx+j-1) * sinkx I(ew_nx,count,ew_nx+1);
        coskx I(ew_nx,count,ew_nx-j) =   coskx I(ew_nx,count,ew_nx+j);
        sinkx I(ew_nz,count,ew_nx+j) =   coskx I(ew_nx,count,ew_nx+j-1) * sinkx I(ew_nx,count,ew_nx+1) 
                                       - sinkx I(ew_nx,count,ew_nx+j-1) * coskx I(ew_nx,count,ew_nx+1);
        sinkx I(ew_nx,count,ew_nx-j) = - sinkx I(ew_nx,count,ew_nx+j);
      }

      for (j=2; j<=ew_ny; j++) {
        cosky I(ew_ny,count,ew_ny+j) =   cosky I(ew_ny,count,ew_ny+j-1) * cosky I(ew_ny,count,ew_ny+1) 
                                       - sinky I(ew_ny,count,ew_ny+j-1) * sinky I(ew_ny,count,ew_ny+1);
        cosky I(ew_ny,count,ew_ny-j) =   cosky I(ew_ny,count,ew_ny+j);
        sinky I(ew_ny,count,ew_ny+j) =   cosky I(ew_ny,count,ew_ny+j-1) * sinky I(ew_ny,count,ew_ny+1) 
                                       - sinky I(ew_ny,count,ew_ny+j-1) * cosky I(ew_ny,count,ew_ny+1);
        sinky I(ew_ny,count,ew_ny-j) = - sinky I(ew_ny,count,ew_ny+j);
      }

      for (j=2; j<=ew_nz; j++) {
        coskz I(ew_nz,count,ew_nz+j) =   coskz I(ew_nz,count,ew_nz+j-1) * coskz I(ew_nz,count,ew_nz+1) 
                                       - sinkz I(ew_nz,count,ew_nz+j-1) * sinkz I(ew_nz,count,ew_nz+1);
        sinkz I(ew_nz,count,ew_nz+j) =   coskz I(ew_nz,count,ew_nz+j-1) * sinkz I(ew_nz,count,ew_nz+1) 
                                       - sinkz I(ew_nz,count,ew_nz+j-1) * coskz I(ew_nz,count,ew_nz+1);
      }
      count++;
    }
  }

  /* Loop over all reciprocal vectors */
  for (k=0; k<ew_totk; k++) {

    l = ew_ivek[k].x;
    m = ew_ivek[k].y;
    n = ew_ivek[k].z;

    sum_cos = 0.0;
    sum_sin = 0.0;

    /* Compute exp(ikr) and sums thereof */
    count = 0;
    for (c=0; c<ncells; c++) {

      cell *p = CELLPTR(c);

      for (i=0; i<p->n; i++) {

        coskr[count] =   coskx I(ew_nx,count,ew_nx+l) * cosky I(ew_ny,count,ew_ny+m) * coskz I(ew_nz,count,ew_nz+n)
                       - sinkx I(ew_nx,count,ew_nx+l) * sinky I(ew_ny,count,ew_ny+m) * coskz I(ew_nz,count,ew_nz+n)
                       - sinkx I(ew_nx,count,ew_nx+l) * cosky I(ew_ny,count,ew_ny+m) * sinkz I(ew_nz,count,ew_nz+n) 
                       - coskx I(ew_nx,count,ew_nx+l) * sinky I(ew_ny,count,ew_ny+m) * sinkz I(ew_nz,count,ew_nz+n);

        sinkr[count] = - sinkx I(ew_nx,count,ew_nx+l) * sinky I(ew_ny,count,ew_ny+m) * sinkz I(ew_nz,count,ew_nz+n)
                       + sinkx I(ew_nx,count,ew_nx+l) * cosky I(ew_ny,count,ew_ny+m) * coskz I(ew_nz,count,ew_nz+n)
                       + coskx I(ew_nx,count,ew_nx+l) * sinky I(ew_ny,count,ew_ny+m) * coskz I(ew_nz,count,ew_nz+n)
                       + coskx I(ew_nx,count,ew_nx+l) * cosky I(ew_ny,count,ew_ny+m) * sinkz I(ew_nz,count,ew_nz+n);

        typ      = SORTE(p,i);
        sum_cos += charge[typ] * coskr[count];
        sum_sin += charge[typ] * sinkr[count]; 

        count++;
      }
    }

    /* update energies and forces */
    count = 0;
    for (c=0; c<ncells; c++) {

      cell *p = CELLPTR(c);

      for (i=0; i<p->n; i++) {

        typ    = SORTE(p,i);

        /* update potential energy */
        kpot   = charge[typ] * ew_expk[k]
	         * (sinkr[count] * sum_sin + coskr[count] * sum_cos);
        POTENG(p,i)    += kpot;
        tot_pot_energy += kpot;

        /* update force */
        kforce = charge[typ] * ew_expk[k] 
                 * (sinkr[count] * sum_cos - coskr[count] * sum_sin);
        KRAFT(p,i,X) += ew_kvek[k].x * kforce;
        KRAFT(p,i,Y) += ew_kvek[k].y * kforce;
        KRAFT(p,i,Z) += ew_kvek[k].z * kforce;
        tmp_virial   += kforce * SPRODX( &ORT(p,i,X), ew_kvek[k] );
	  
        count++;
      }
    }

  }  /* k */

  virial += tmp_virial;

}

/******************************************************************************
*
*  do_forces_ewald_real
*
*  computes the real part of the Ewald sum
*
******************************************************************************/

void do_forces_ewald_real(void)
{
  int i, j, k, l, m, r, s, jstart; 
  int p_typ, q_typ;
  cell *p, *q;
  vektor d, tmp_d, rforce;
  real tmp_sprod, tmp_boxl, radius2, radius, rpot, charge2, tmp_rforce;
  real tmp_rpot, tmp_virial=0.0;

  /* Loop over all pairs of cells */
  for(r=0; r<ncells; ++r)
    for(s=r; s<ncells; ++s)
      {
	p = CELLPTR(r);
	q = CELLPTR(s);

	/* For each atom in first cell */
	for (i=0; i<p->n; ++i) 
	  {
	    p_typ = SORTE(p,i);

	    /* For each atom in second cell */
            jstart = (p==q ? i : 0 );

	    for (j=jstart; j<q->n; ++j)
	      {
	  
		q_typ = SORTE(q,j);

		/* Distance vector */

		d.x = ORT(p,i,X) - ORT(q,j,X);
		d.y = ORT(p,i,Y) - ORT(q,j,Y);
		d.z = ORT(p,i,Z) - ORT(q,j,Z);

		/* Apply the minimum image convention */
		tmp_sprod = SPROD(d,box_x);
		tmp_boxl  = 0.5 * SPROD(box_x,box_x);
		if ( tmp_sprod > tmp_boxl )
		  {
		    d.x -= box_x.x;
		    d.y -= box_x.y;
		    d.z -= box_x.z;
		  }
		if ( tmp_sprod < -tmp_boxl ) 
		  {
		    d.x += box_x.x;
		    d.y += box_x.y;
		    d.z += box_x.z;
		  }
	      
		tmp_sprod = SPROD(d,box_y);
		tmp_boxl  = 0.5 * SPROD(box_y,box_y);
		if ( tmp_sprod > tmp_boxl ) 
		  {
		    d.x -= box_y.x;
		    d.y -= box_y.y;
		    d.z -= box_y.z;
		  }
		if ( tmp_sprod < -tmp_boxl ) 
		  {
		    d.x += box_y.x;
		    d.y += box_y.y;
		    d.z += box_y.z;
		  }
	      
		tmp_sprod = SPROD(d,box_z);
		tmp_boxl  = 0.5 * SPROD(box_z,box_z);
		if ( tmp_sprod > tmp_boxl ) 
		  {
		    d.x -= box_z.x;
		    d.y -= box_z.y;
		    d.z -= box_z.z;
		  }
		if ( tmp_sprod < -tmp_boxl ) 
		  {
		    d.x += box_z.x;
		    d.y += box_z.y;
		    d.z += box_z.z;
		  }
	      
		rpot = 0.0; rforce.x = 0.0; rforce.y = 0.0; rforce.z = 0.0;
		
		charge2 = charge[p_typ] * charge[q_typ] * ew_eps;

		/* Include image boxes */
		for( k=ew_nmax; k>=-ew_nmax; k--)
		  for( l=ew_nmax; l>=-ew_nmax; l--)
		    for( m=ew_nmax; m>=-ew_nmax; m--)
		      {
		
			tmp_d.x = d.x + k*box_x.x + l*box_y.x + m*box_z.x;
			tmp_d.y = d.y + k*box_x.y + l*box_y.y + m*box_z.y;
			tmp_d.z = d.z + k*box_x.z + l*box_y.z + m*box_z.z;

			radius2 = SPROD(tmp_d,tmp_d);

			if( radius2 > 1.0e-6 )
			  {
			    radius  = sqrt(radius2);
			    
			    /* real space part of potential and force */
			    tmp_rpot = charge2 * erfc1(ew_kappa*radius)/radius;
			    
			    rpot  += tmp_rpot;
			  
			    tmp_rforce = ( tmp_rpot + charge2 * 2.0 * ew_vorf 
                                       * exp(-SQR(ew_kappa*radius))) / radius2;

			    rforce.x += tmp_d.x * tmp_rforce;
			    rforce.y += tmp_d.y * tmp_rforce;
			    rforce.z += tmp_d.z * tmp_rforce;

			    tmp_virial += radius2 * tmp_rforce;
			}
		      }

		/* update potential energy and forces */

		POTENG(p,i) += rpot;
		POTENG(q,j) += rpot;

                /*db_rpot += rpot;*/

                tot_pot_energy += rpot;

		KRAFT(p,i,X) += rforce.x;
		KRAFT(q,j,X) -= rforce.x;
	
		KRAFT(p,i,Y) += rforce.y;
		KRAFT(q,j,Y) -= rforce.y;
	
		KRAFT(p,i,Z) += rforce.z;
		KRAFT(q,j,Z) -= rforce.z;

	      } /* for j */
	    
	  } /* for i */

	virial += tmp_virial;
     
      }

}


/******************************************************************************
*
*  init_ewald
*
*  setup reciprocal lattice vectors
*
******************************************************************************/

void init_ewald(void)
{
  int    i, j, k;
  real   kvek2, vorf1;
  int    count, num;

  /* we implicitly assume a system of units, in which lengths are
     measured in Angstrom [A], energies in electron volts [eV], 
     and charges in units of the elementary charge e */

  /* ew_eps is already initialized to 14.38 in globals.h    */
  /* this is e^2 / (2*pi*epsilon_0) in [eV A]               */
  /* we need it earlier for analytically defined potentials */
  /* ew_eps   = 14.38; */
  twopi    = 2.0 * M_PI;
  ew_vorf  = ew_kappa / sqrt( M_PI );
  vorf1    = 2.0 * twopi * ew_eps / volume;

  ew_nx = (int) (ew_kcut * SQRT( SPROD(box_x,box_x) ) / twopi) + 1;
  ew_ny = (int) (ew_kcut * SQRT( SPROD(box_y,box_y) ) / twopi) + 1;
  ew_nz = (int) (ew_kcut * SQRT( SPROD(box_z,box_z) ) / twopi) + 1;
  printf("EWALD: nx = %d, ny = %d, nz = %d\n", ew_nx, ew_ny, ew_nz);

  /* Allocate memory for k-vectors */
  num     = (2*ew_nx+1) * (2*ew_ny+1) * (ew_nz+1);
  ew_kvek = (vektor  *) malloc( num * sizeof(vektor) );
  ew_ivek = (ivektor *) malloc( num * sizeof(ivektor) );
  ew_expk = (real    *) malloc( num * sizeof(real) );

  if( (ew_kvek == NULL) || (ew_ivek == NULL) || (ew_expk == NULL) )
    error("EWALD: Cannot allocate memory for k-vectors");

  /* Determine k-vectors */
  ew_totk = 0;
  for (i = -ew_nx; i <= ew_nx; i++)
    for (j = -ew_ny; j <= ew_ny; j++)
      for (k = 0; k <= ew_nz; k++) {

        ew_kvek[ew_totk].x = twopi * (i*tbox_x.x + j*tbox_y.x + k*tbox_z.x);
        ew_kvek[ew_totk].y = twopi * (i*tbox_x.y + j*tbox_y.y + k*tbox_z.y);
        ew_kvek[ew_totk].z = twopi * (i*tbox_x.z + j*tbox_y.z + k*tbox_z.z);

        /* skip, if outside cutoff radius */
        kvek2 = SPROD( ew_kvek[ew_totk], ew_kvek[ew_totk] );
        if (kvek2 > SQR(ew_kcut)) continue;

        /* only half of the k-vectors are needed */
        if ( (k==0) && ( (j<0) || ( (j==0) && (i<1) ) ) ) continue;

        ew_ivek[ew_totk].x = i;
        ew_ivek[ew_totk].y = j;
        ew_ivek[ew_totk].z = k;

        ew_expk[ew_totk] = vorf1 * exp( -kvek2 / (4.0*SQR(ew_kappa)) ) / kvek2;

        ew_totk++;

      }

  printf("EWALD: %d k-vectors\n", ew_totk);

  /* Allocate memory for exp(ikr) */
  coskx = (real *) malloc( natoms * (2*ew_nx+1) * sizeof(real));
  sinkx = (real *) malloc( natoms * (2*ew_nx+1) * sizeof(real));
  cosky = (real *) malloc( natoms * (2*ew_ny+1) * sizeof(real));
  sinky = (real *) malloc( natoms * (2*ew_ny+1) * sizeof(real));
  coskz = (real *) malloc( natoms * (2*ew_nz+1) * sizeof(real));
  sinkz = (real *) malloc( natoms * (2*ew_nz+1) * sizeof(real));
  coskr = (real *) malloc( natoms               * sizeof(real));
  sinkr = (real *) malloc( natoms               * sizeof(real));

  if( coskx == NULL || sinkx == NULL || cosky == NULL || sinky == NULL
      || coskz == NULL || sinkz == NULL || coskr == NULL || sinkr == NULL )
    error("EWALD: Cannot allocate memory for exp(ikr)");

  /* Position independent initializations */
  count = 0;
  for (k=0; k<ncells; k++) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; i++) {
      coskx I(ew_nx,count,ew_nx) = 1.0;
      sinkx I(ew_nx,count,ew_nx) = 0.0;
      cosky I(ew_ny,count,ew_ny) = 1.0;
      sinky I(ew_ny,count,ew_ny) = 0.0;
      coskz I(ew_nz,count,ew_nz) = 1.0;
      sinkz I(ew_nz,count,ew_nz) = 0.0;
      count++;
    }
  }
}


/******************************************************************************
*
*  erfc1
*
*  Approximation of erfc() 
*
******************************************************************************/

real erfc1(real x)
{

#define A_1 0.254829592
#define A_2 -0.284496736
#define A_3 1.421413741
#define A_4 -1.453152027
#define A_5 1.061405429
#define P   0.3275911

  real t, xsq, tp;

  t = 1.0 / ( 1.0 + P * x );

  xsq = x * x;

  tp = t * ( A_1 + t * ( A_2 + t * ( A_3 + t * ( A_4 + t * A_5 ) ) ) );
  
  return ( tp * exp( -xsq ) );

} 
