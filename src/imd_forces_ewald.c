
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

/******************************************************************************
*
*  do_forces_ewald
*
******************************************************************************/

void do_forces_ewald(int steps)
{
  int n;
  real old1, old2;

  if (steps==0) {
    old1 = tot_pot_energy / natoms;
    printf( "Ewald pair int: Epot = %f, time = %f\n", 
            old1, ewald_time.total);
    ewald_time.total = 0.0;
    imd_start_timer( &ewald_time );
  }

  /* real space part; if ew_nmax < 0, we do it with the other pair forces */
  if (ew_nmax >= 0) do_forces_ewald_real();

  if (steps==0) {
    imd_stop_timer( &ewald_time );
    old2 = tot_pot_energy / natoms;
    printf( "Ewald real space: Epot = %f, time = %f\n", 
            old2 - old1, ewald_time.total);
    ewald_time.total = 0.0;
    imd_start_timer( &ewald_time );
  }

  /* Fourier space part */
  do_forces_ewald_fourier();

  if (steps==0) {
    imd_stop_timer( &ewald_time );
    old1 = tot_pot_energy / natoms;
    printf( "Ewald Fourier space: Epot = %f, time = %f\n", 
            old1 - old2, ewald_time.total);
    ewald_time.total = 0.0;
    imd_start_timer( &ewald_time );
  }

  /* self energy part */
  for (n=0; n<ntypes; n++)
    tot_pot_energy -= ew_vorf * num_sort[n] * SQR( charge[n] ) * ew_eps;

  if (steps==0) {
    imd_stop_timer( &ewald_time );
    old2 = tot_pot_energy / natoms;
    printf( "Ewald self-energy: Epot = %f, time = %f\n", 
            old2 - old1, ewald_time.total);
    printf( "Ewald total: Epot = %f\n", old2);
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

  int    i, j, k, l, m, n;
  cell   *p;
  int    count, totk, typ;
  vektor coord;
  real   tmp, tmp_virial=0.0, sum_cos, sum_sin;
  real   kforce, kpot;

  count = 0;

  /* Compute exp(ikr) recursively */
  for( k=0; k<ncells; ++k)
    {
      p = CELLPTR(k);
      for( i=0; i<p->n; ++i)
	{
	  coord.x = ORT(p,i,X);
	  coord.y = ORT(p,i,Y);
	  coord.z = ORT(p,i,Z);

	  tmp = twopi * SPROD(tbox_x, coord);
	  coskx I(ew_nx,count,ew_nx+1) =  cos( tmp );
	  coskx I(ew_nx,count,ew_nx-1) =  coskx I(ew_nx,count,ew_nx+1);
	  sinkx I(ew_nx,count,ew_nx+1) =  sin( tmp );
	  sinkx I(ew_nx,count,ew_nx-1) = -sinkx I(ew_nx,count,ew_nx+1); 

	  tmp = twopi * SPROD(tbox_y, coord);
	  cosky I(ew_ny,count,ew_ny+1) =  cos( tmp );
	  cosky I(ew_ny,count,ew_ny-1) =  cosky I(ew_ny,count,ew_ny+1);
	  sinky I(ew_ny,count,ew_ny+1) =  sin( tmp);
	  sinky I(ew_ny,count,ew_ny-1) = -sinky I(ew_ny,count,ew_ny+1);

	  tmp = twopi * SPROD(tbox_z, coord);
	  coskz I(ew_nz,count,ew_nz+1) = cos( tmp );
	  sinkz I(ew_nz,count,ew_nz+1) = sin( tmp );

	  for( j=2; j<=ew_nx; ++j)
	    {
	      coskx I(ew_nx,count,ew_nx+j) =   coskx I(ew_nx,count,ew_nx+j-1) * coskx I(ew_nx,count,ew_nx+1) 
		                             - sinkx I(ew_nx,count,ew_nx+j-1) * sinkx I(ew_nx,count,ew_nx+1);
	      coskx I(ew_nx,count,ew_nx-j) =   coskx I(ew_nx,count,ew_nx+j);
	      sinkx I(ew_nz,count,ew_nx+j) =   coskx I(ew_nx,count,ew_nx+j-1) * sinkx I(ew_nx,count,ew_nx+1) 
		                             - sinkx I(ew_nx,count,ew_nx+j-1) * coskx I(ew_nx,count,ew_nx+1);
	      sinkx I(ew_nx,count,ew_nx-j) = - sinkx I(ew_nx,count,ew_nx+j);
	    }

	  for( j=2; j<=ew_ny; ++j)
	    {
	      cosky I(ew_ny,count,ew_ny+j) =   cosky I(ew_ny,count,ew_ny+j-1) * cosky I(ew_ny,count,ew_ny+1) 
		                             - sinky I(ew_ny,count,ew_ny+j-1) * sinky I(ew_ny,count,ew_ny+1);
	      cosky I(ew_ny,count,ew_ny-j) =   cosky I(ew_ny,count,ew_ny+j);
	      sinky I(ew_ny,count,ew_ny+j) =   cosky I(ew_ny,count,ew_ny+j-1) * sinky I(ew_ny,count,ew_ny+1) 
		                             - sinky I(ew_ny,count,ew_ny+j-1) * cosky I(ew_ny,count,ew_ny+1);
	      sinky I(ew_ny,count,ew_ny-j) = - sinky I(ew_ny,count,ew_ny+j);
	    }

	  for( j=2; j<=ew_nz; ++j)
	    {
	      coskz I(ew_nz,count,ew_nz+j) =   coskz I(ew_nz,count,ew_nz+j-1) * coskz I(ew_nz,count,ew_nz+1) 
		                             - sinkz I(ew_nz,count,ew_nz+j-1) * sinkz I(ew_nz,count,ew_nz+1);
	      sinkz I(ew_nz,count,ew_nz+j) =   coskz I(ew_nz,count,ew_nz+j-1) * sinkz I(ew_nz,count,ew_nz+1) 
		                             - sinkz I(ew_nz,count,ew_nz+j-1) * coskz I(ew_nz,count,ew_nz+1);
	    }

	  ++count;

	}
    }

  /* Loop over all reciprocal vectors */
  for (l=-ew_nx; l<=ew_nx; ++l)
    for (m=-ew_ny; m<=ew_ny; ++m)
      for (n=0; n<=ew_nz; ++n) {

        totk = n + (m+ew_ny)*(ew_nz+1) + (l+ew_nx)*(2*ew_ny+1)*(ew_nz+1);

        sum_cos = 0.0;
        sum_sin = 0.0;

        count = 0;

        /* Compute exp(ikr) and sums thereof */
        for (k=0; k<ncells; ++k) {

          cell *p = CELLPTR(k);

	  for( i=0; i<p->n; ++i) {

            typ = SORTE(p,i);

            coskr[count] =   coskx I(ew_nx,count,ew_nx+l) * cosky I(ew_ny,count,ew_ny+m) * coskz I(ew_nz,count,ew_nz+n)
                           - sinkx I(ew_nx,count,ew_nx+l) * sinky I(ew_ny,count,ew_ny+m) * coskz I(ew_nz,count,ew_nz+n)
                           - sinkx I(ew_nx,count,ew_nx+l) * cosky I(ew_ny,count,ew_ny+m) * sinkz I(ew_nz,count,ew_nz+n) 
                           - coskx I(ew_nx,count,ew_nx+l) * sinky I(ew_ny,count,ew_ny+m) * sinkz I(ew_nz,count,ew_nz+n);

            sinkr[count] = - sinkx I(ew_nx,count,ew_nx+l) * sinky I(ew_ny,count,ew_ny+m) * sinkz I(ew_nz,count,ew_nz+n)
                           + sinkx I(ew_nx,count,ew_nx+l) * cosky I(ew_ny,count,ew_ny+m) * coskz I(ew_nz,count,ew_nz+n)
                           + coskx I(ew_nx,count,ew_nx+l) * sinky I(ew_ny,count,ew_ny+m) * coskz I(ew_nz,count,ew_nz+n)
                           + coskx I(ew_nx,count,ew_nx+l) * cosky I(ew_ny,count,ew_ny+m) * sinkz I(ew_nz,count,ew_nz+n);

            sum_cos += charge[typ] * coskr[count];
            sum_sin += charge[typ] * sinkr[count]; 

            ++count;
	  }
	}

        /* update potential energy */
        kpot = ew_vorf1 * ew_expk[totk] * ( SQR(sum_sin) + SQR(sum_cos) );		  
        POTENG(p,i)    += kpot;
        tot_pot_energy += kpot;

        count = 0;

        /* update forces */
        for (k=0; k<ncells; ++k) {

          cell *p = CELLPTR(k);

          for (i=0; i<p->n; ++i) {

            typ = SORTE(p,i);
            kforce = ew_vorf1 * charge[typ] * ew_expk[totk] 
                     * (sinkr[count] * sum_cos - coskr[count] * sum_sin);

            KRAFT(p,i,X) += ew_kvek[totk].x * kforce;
            KRAFT(p,i,Y) += ew_kvek[totk].y * kforce;
            KRAFT(p,i,Z) += ew_kvek[totk].z * kforce;

            tmp_virial += kforce * 
                          ( ew_kvek[totk].x * ORT(p,i,X) + 
                            ew_kvek[totk].y * ORT(p,i,Y) + 
                            ew_kvek[totk].z * ORT(p,i,Z) );  
	  
            ++count;

	  }
	}

      }  /* l, m, n */

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
  real   kvek2;
  int    count, num;
  cell   *p;

  /* Some constants and prefactors */
  /* first calculate volume of box */
  volume = box_x.x * ( box_y.y * box_z.z - box_y.z * box_z.y)
           - box_x.y * ( box_y.x * box_z.z - box_y.z * box_z.x)
           + box_x.z * ( box_y.x * box_z.y - box_y.y * box_z.x); 

  /* we implicitly assume a system of units, in which lengths are
     measured in Angstrom [A], energies in electron volts [eV], 
     and charges in units of the elementary charge e */

  /* ew_eps is already initialized to 14.38 in globals.h    */
  /* this is e^2 / (2*pi*epsilon_0) in [eV A]               */
  /* we need it earlier for analytically defined potentials */
  /* ew_eps   = 14.38; */
  twopi    = 2.0 * M_PI;
  ew_vorf  = ew_kappa / sqrt( M_PI );
  ew_vorf1 = 2.0 * twopi * ew_eps / volume;

  ew_nx = ew_kmax.x;
  ew_ny = ew_kmax.y;
  ew_nz = ew_kmax.z;

  /* Allocate memory for k-vectors */
  num     = (2*ew_nx+1) * (2*ew_ny+1) * (ew_nz+1);
  ew_kvek = (vektor *) malloc( num * sizeof(vektor) );
  ew_expk = (real   *) malloc( num * sizeof(real) );

  if( (ew_kvek == NULL) || (ew_expk == NULL) )
    error("EWALD: Cannot allocate memory for k-vectors");

  ew_totk  = 0;
	  
  /* Determine k-vectors */
  for( i = -ew_nx; i <= ew_nx; ++i )
    for( j = -ew_ny; j <= ew_ny; ++j )
      for( k = 0; k <= ew_nz; ++k )
	{
	  
	  ew_kvek[ew_totk].x = twopi * (i*tbox_x.x + j*tbox_y.x + k*tbox_z.x );
	  ew_kvek[ew_totk].y = twopi * (i*tbox_x.y + j*tbox_y.y + k*tbox_z.y );
	  ew_kvek[ew_totk].z = twopi * (i*tbox_x.z + j*tbox_y.z + k*tbox_z.z );

	  kvek2 = SPROD( ew_kvek[ew_totk], ew_kvek[ew_totk] );

	  /* Only half of the k-vectors will be needed */
	  if ( k==0 && ( j<0 || ( j==0 && i<1 ) ) )
	    {
	      ew_kvek[ew_totk].x = 0.0;
	      ew_kvek[ew_totk].y = 0.0;
	      ew_kvek[ew_totk].z = 0.0;
	      
	      ew_expk[ew_totk]   = 0.0;
	    }
	  else
	    ew_expk[ew_totk] = exp( -kvek2 / (4.0*ew_kappa*ew_kappa) ) / kvek2;

	  ++ew_totk;

	}

  printf("EWALD: %d k-vectors\n", (2*ew_nx+1)*(2*ew_ny+1)*(2*ew_nz+1));

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

  count = 0;

  /* Position independent initializations */
  for( k=0; k<ncells; ++k)
    {
      p = CELLPTR(k);
      for( i=0; i<p->n; ++i)
	{
	  coskx I(ew_nx,count,ew_nx) = 1.0;
	  sinkx I(ew_nx,count,ew_nx) = 0.0;

	  cosky I(ew_ny,count,ew_ny) = 1.0;
	  sinky I(ew_ny,count,ew_ny) = 0.0;

	  coskz I(ew_nz,count,ew_nz) = 1.0;
	  sinkz I(ew_nz,count,ew_nz) = 0.0;
      
	  ++count;

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
