
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
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
      p = cell_array + k;
      for( i=0; i<p->n; ++i)
	{
	  coord.x = p->ort X(i);
	  coord.y = p->ort Y(i);
	  coord.z = p->ort Z(i);

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
  for( l=-ew_nx; l<=ew_nx; ++l)
    for( m=-ew_ny; m<=ew_ny; ++m)
      for( n=0; n<=ew_nz; ++n)
	{
	  totk = n + (m+ew_ny)*(ew_nz+1) + (l+ew_nx)*(2*ew_ny+1)*(ew_nz+1);

	  sum_cos = 0.0;
	  sum_sin = 0.0;

	  count = 0;

	  /* Compute exp(ikr) and sums thereof */
	  for( k=0; k<ncells; ++k)
	    {
	      p = cell_array + k;
	      for( i=0; i<p->n; ++i)
		{


		  typ = SORTE(p,i);
	  
		  coskr J(ew_totk,count,totk) =   coskx I(ew_nx,count,ew_nx+l) * cosky I(ew_ny,count,ew_ny+m) * coskz I(ew_nz,count,ew_nz+n)
		                       - sinkx I(ew_nx,count,ew_nx+l) * sinky I(ew_ny,count,ew_ny+m) * coskz I(ew_nz,count,ew_nz+n)
		                       - sinkx I(ew_nx,count,ew_nx+l) * cosky I(ew_ny,count,ew_ny+m) * sinkz I(ew_nz,count,ew_nz+n) 
		                       - coskx I(ew_nx,count,ew_nx+l) * sinky I(ew_ny,count,ew_ny+m) * sinkz I(ew_nz,count,ew_nz+n);

		  sinkr J(ew_totk,count,totk) = - sinkx I(ew_nx,count,ew_nx+l) * sinky I(ew_ny,count,ew_ny+m) * sinkz I(ew_nz,count,ew_nz+n)
		                       + sinkx I(ew_nx,count,ew_nx+l) * cosky I(ew_ny,count,ew_ny+m) * coskz I(ew_nz,count,ew_nz+n)
		                       + coskx I(ew_nx,count,ew_nx+l) * sinky I(ew_ny,count,ew_ny+m) * coskz I(ew_nz,count,ew_nz+n)
		                       + coskx I(ew_nx,count,ew_nx+l) * cosky I(ew_ny,count,ew_ny+m) * sinkz I(ew_nz,count,ew_nz+n);

		  sum_cos += charge[typ] * coskr J(ew_totk,count,totk);
		  sum_sin += charge[typ] * sinkr J(ew_totk,count,totk); 

		  ++count;

		}
	    }

	  count = 0;

	  /* Update potential and forces */
	  for( k=0; k<ncells; ++k)
	    {
	      p = cell_array + k;
	      for( i=0; i<p->n; ++i)
		{

		  typ = SORTE(p,i);
	
		  kpot = ew_vorf1 * ew_expk[totk] * charge[typ] 
		             * ( sinkr J(ew_totk,count,totk) * sum_sin + coskr J(ew_totk,count,totk) * sum_cos );
/*db_kpot += kpot;*/		  
		  p->pot_eng[i] += kpot;

		  tot_pot_energy += kpot;

		  kforce = ew_vorf1 * charge[typ] * ew_expk[totk] 
		               * ( sinkr J(ew_totk,count,totk) * sum_cos - coskr J(ew_totk,count,totk) * sum_sin );
	  		  
		  p->kraft X(i) += ew_kvek[totk].x * kforce;
		  p->kraft Y(i) += ew_kvek[totk].y * kforce;
		  p->kraft Z(i) += ew_kvek[totk].z * kforce;

		  tmp_virial += kforce * 
		                ( ew_kvek[totk].x * p->ort X(i) + ew_kvek[totk].y * p->ort Y(i) + ew_kvek[totk].z * p->ort Z(i) );  
		  
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
  real tmp_sprod, tmp_boxl, radius2, radius, rpot, spot, charge2, tmp_rforce;
  real tmp_rpot, tmp_virial=0.0;

  spot = 0.0;

  /* Loop over all pairs of cells */
  for(r=0; r<ncells; ++r)
    for(s=r; s<ncells; ++s)
      {
	p = cell_array + r;
	q = cell_array + s;

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

		d.x = p->ort X(i) - q->ort X(j);
		d.y = p->ort Y(i) - q->ort Y(j);
		d.z = p->ort Z(i) - q->ort Z(j);

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
			    tmp_rpot = charge2 * erfc1(ew_kappa * radius) / radius;
			    
			    rpot  += tmp_rpot;
			  
			    tmp_rforce = ( tmp_rpot + charge2 * 2.0 * ew_vorf * exp(-SQR(ew_kappa*radius))) / radius2; 

			    rforce.x += tmp_d.x * tmp_rforce;
			    rforce.y += tmp_d.y * tmp_rforce;
			    rforce.z += tmp_d.z * tmp_rforce;

			    tmp_virial += radius2 * tmp_rforce;
			}
		      }

		/* update potential energy and forces */

		p->pot_eng[i]  += rpot;
		q->pot_eng[j]  += rpot;

/*db_rpot += rpot;*/
	  
                tot_pot_energy += rpot;

		p->kraft X(i) += rforce.x;
		q->kraft X(j) -= rforce.x;
	
		p->kraft Y(i) += rforce.y;
		q->kraft Y(j) -= rforce.y;
	
		p->kraft Z(i) += rforce.z;
		q->kraft Z(j) -= rforce.z;

	      } /* for j */
	    
	    /* self energy part */
	    if( p == q )
	      spot += charge[p_typ] * charge[p_typ] * ew_eps;
  
	  } /* for i */

	virial += tmp_virial;
     
      }

/*db_spot -= ew_vorf * spot;*/ 
     tot_pot_energy -= ew_vorf * spot;

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
  int    count, num=100;
  cell   *p;

  /* Some constants and prefactors */
  /* first calculate volume of box */
  volume = box_x.x * ( box_y.y * box_z.z - box_y.z * box_z.y)
           - box_x.y * ( box_y.x * box_z.z - box_y.z * box_z.x)
           + box_x.z * ( box_y.x * box_z.y - box_y.y * box_z.x); 

  ew_eps   = 14.38;
  twopi    = 2.0 * M_PI;
  ew_vorf  = ew_kappa / sqrt( M_PI );
  ew_vorf1 = 2.0 * twopi * ew_eps / volume;

  ew_nx = ew_kmax.x;
  ew_ny = ew_kmax.y;
  ew_nz = ew_kmax.z;

  /* Allocate memory for k-vectors */
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

	  if( ew_totk >= num )
	    {
	      num += 10;
	      ew_kvek = (vektor *) realloc( ew_kvek, num * sizeof(vektor));
	      ew_expk = (real   *) realloc( ew_expk, num * sizeof(real));
	    }  
	  
	  ew_kvek[ew_totk].x = twopi * ( i * tbox_x.x + j * tbox_y.x + k * tbox_z.x );
	  ew_kvek[ew_totk].y = twopi * ( i * tbox_x.y + j * tbox_y.y + k * tbox_z.y );
	  ew_kvek[ew_totk].z = twopi * ( i * tbox_x.z + j * tbox_y.z + k * tbox_z.z );

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
	    ew_expk[ew_totk] = exp( - kvek2 / ( 4.0 * ew_kappa * ew_kappa ) ) / kvek2;

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
  coskr = (real *) malloc( natoms * ew_totk     * sizeof(real));
  sinkr = (real *) malloc( natoms * ew_totk     * sizeof(real));

  if( coskx == NULL || sinkx == NULL || cosky == NULL || sinky == NULL
      || coskz == NULL || sinkz == NULL || coskr == NULL || sinkr == NULL )
    error("EWALD: Cannot allocate memory for exp(ikr)");

  count = 0;

  /* Position independent initializations */
  for( k=0; k<ncells; ++k)
    {
      p = cell_array + k;
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
	      







