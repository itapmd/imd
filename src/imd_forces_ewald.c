
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2012 Institute for Theoretical and Applied Physics,
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
  tot_pot_energy = 0.0;
}

real force_norm(void)
{
#ifdef DEBUG
  printf("force_norm\n");
#endif

  int k;
  double norm = 0.0;

#ifdef BUFCELLS
  /* collect forces from buffer cells */
  send_forces(add_forces,pack_forces,unpack_forces);
#endif

  for (k=0; k<ncells; k++) {
    cell *p = CELLPTR(k);
    int i;
    for (i=0; i<p->n; i++) {
      norm += SQR( KRAFT(p,i,X) );
      norm += SQR( KRAFT(p,i,Y) );
      norm += SQR( KRAFT(p,i,Z) );
    }
  }
  return sqrt(norm / (DIM*natoms));
}

/******************************************************************************
*
*  do_forces_ewald
*
******************************************************************************/

void do_forces_ewald(int steps)
{

  int n, k, i, typ;
  real tot=0.0, pot;

  /* real space part; if ew_nmax < 0, we do it with the other pair forces */
  if (ew_nmax >= 0) do_forces_ewald_real();

  if ((steps==0) && (ew_test)) {
    imd_stop_timer( &ewald_time );
    tot += tot_pot_energy / natoms;
    printf( "Ewald pair:   Epot = %e, force = %e, time = %e\n", 
            tot_pot_energy / natoms, force_norm(), ewald_time.total);
    ewald_time.total = 0.0;
    imd_start_timer( &ewald_time );
    clear_forces();
  }

  /* Fourier space part */
  if (ew_kcut > 0) do_forces_ewald_fourier();

  if ((steps==0) && (ew_test) && (ew_kcut>0)) {
    imd_stop_timer( &ewald_time );
    tot += tot_pot_energy / natoms;
    printf( "Ewald smooth: Epot = %e, force = %e, time = %e\n", 
            tot_pot_energy / natoms, force_norm(), ewald_time.total);
    ewald_time.total = 0.0;
    imd_start_timer( &ewald_time );
    clear_forces();
  }

  /* self energy part */
  for (k=0; k<ncells; k++) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; i++) {
#ifdef VARCHG
      pot  = ew_vorf * SQR( CHARGE(p,i) ) * coul_eng;
#else
      typ  = SORTE(p,i);
      pot  = ew_vorf * SQR( charge[typ] ) * coul_eng;
#endif
      tot_pot_energy -= pot;
      POTENG(p,i)    -= pot;
    }
  }

  if ((steps==0) && (ew_test)) {
    tot += tot_pot_energy / natoms;
    printf( "Ewald self:   Epot = %e\n", tot_pot_energy / natoms );
    printf( "Ewald total:  Epot = %f\n", tot);
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
  int    cnt, typ, offx, offy, offz;
  int    px, py, pz, mx, my, mz;
  real   tmp, tmp_virial=0.0, sum_cos, sum_sin;
  real   kforce, kpot;

  /* Compute exp(ikr) recursively */
  px = (ew_nx+1) * natoms;  mx = (ew_nx-1) * natoms;
  py = (ew_ny+1) * natoms;  my = (ew_ny-1) * natoms;
  pz = (ew_nz+1) * natoms;  mz = (ew_nz-1) * natoms;
  cnt = 0;
  for (c=0; c<ncells; c++) {
    cell *p = CELLPTR(c);
    for (i=0; i<p->n; i++) {
      tmp = twopi * SPRODX(ORT,p,i,tbox_x);
      coskx[px+cnt] =  cos(tmp);
      coskx[mx+cnt] =  coskx[px+cnt];
      sinkx[px+cnt] =  sin(tmp);
      sinkx[mx+cnt] = -sinkx[px+cnt]; 
      tmp = twopi * SPRODX(ORT,p,i,tbox_y);
      cosky[py+cnt] =  cos(tmp);
      cosky[my+cnt] =  cosky[py+cnt];
      sinky[py+cnt] =  sin(tmp);
      sinky[my+cnt] = -sinky[py+cnt];
      tmp = twopi * SPRODX(ORT,p,i,tbox_z);
      coskz[pz+cnt] = cos(tmp);
      sinkz[pz+cnt] = sin(tmp);
      cnt++;
    }
  }

  for (j=2; j<=ew_nx; j++) {
    int pp, qq, mm, ee, i;
    pp  = (ew_nx+j  ) * natoms;
    qq  = (ew_nx+j-1) * natoms;
    mm  = (ew_nx-j  ) * natoms;
    ee  = (ew_nx  +1) * natoms;
    for (i=0; i<natoms; i++) {
      coskx[pp+i] =   coskx[qq+i] * coskx[ee+i] - sinkx[qq+i] * sinkx[ee+i];
      coskx[mm+i] =   coskx[pp+i];
      sinkx[pp+i] =   coskx[qq+i] * sinkx[ee+i] - sinkx[qq+i] * coskx[ee+i];
      sinkx[mm+i] = - sinkx[pp+i];
    }
  }

  for (j=2; j<=ew_ny; j++) {
    int pp, qq, mm, ee, i;
    pp  = (ew_ny+j  ) * natoms;
    qq  = (ew_ny+j-1) * natoms;
    mm  = (ew_ny-j  ) * natoms;
    ee  = (ew_ny  +1) * natoms;
    for (i=0; i<natoms; i++) {
      cosky[pp+i] =   cosky[qq+i] * cosky[ee+i] - sinky[qq+i] * sinky[ee+i];
      cosky[mm+i] =   cosky[pp+i];
      sinky[pp+i] =   cosky[qq+i] * sinky[ee+i] - sinky[qq+i] * cosky[ee+i];
      sinky[mm+i] = - sinky[pp+i];
    }
  }

  for (j=2; j<=ew_nz; j++) {
    int pp, qq, ee, i;
    pp  = (ew_nz+j  ) * natoms;
    qq  = (ew_nz+j-1) * natoms;
    ee  = (ew_nz  +1) * natoms;
    for (i=0; i<natoms; i++) {
      coskz[pp+i] =   coskz[qq+i] * coskz[ee+i] - sinkz[qq+i] * sinkz[ee+i];
      sinkz[pp+i] =   coskz[qq+i] * sinkz[ee+i] - sinkz[qq+i] * coskz[ee+i];
    }
  }

  /* Loop over all reciprocal vectors */
  for (k=0; k<ew_totk; k++) {

    cnt     = 0; 
    offx    = ew_ivek[k].x * natoms;
    offy    = ew_ivek[k].y * natoms;
    offz    = ew_ivek[k].z * natoms;
    sum_cos = 0.0;
    sum_sin = 0.0;

    /* Compute exp(ikr) and sums thereof */
    for (c=0; c<ncells; c++) {

      cell *p = CELLPTR(c);

      for (i=0; i<p->n; i++) {

        coskr[cnt] =   coskx[offx+cnt] * cosky[offy+cnt] * coskz[offz+cnt]
                     - sinkx[offx+cnt] * sinky[offy+cnt] * coskz[offz+cnt]
                     - sinkx[offx+cnt] * cosky[offy+cnt] * sinkz[offz+cnt] 
                     - coskx[offx+cnt] * sinky[offy+cnt] * sinkz[offz+cnt];

        sinkr[cnt] = - sinkx[offx+cnt] * sinky[offy+cnt] * sinkz[offz+cnt]
                     + sinkx[offx+cnt] * cosky[offy+cnt] * coskz[offz+cnt]
                     + coskx[offx+cnt] * sinky[offy+cnt] * coskz[offz+cnt]
                     + coskx[offx+cnt] * cosky[offy+cnt] * sinkz[offz+cnt];

#ifdef VARCHG
        sum_cos += CHARGE(p,i) * coskr[cnt];
        sum_sin += CHARGE(p,i) * sinkr[cnt];
#else
        typ      = SORTE(p,i);
        sum_cos += charge[typ] * coskr[cnt];
        sum_sin += charge[typ] * sinkr[cnt];
#endif
        cnt++;
      }
    }

    /* update total potential energy */
    tot_pot_energy += ew_expk[k] * (SQR(sum_sin) + SQR(sum_cos));

    /* updates */
    cnt = 0;
    for (c=0; c<ncells; c++) {

      cell *p = CELLPTR(c);

      for (i=0; i<p->n; i++) {

        /* update potential energy */
#ifdef VARCHG
        kpot   = CHARGE(p,i) * ew_expk[k]
	         * (sinkr[cnt] * sum_sin + coskr[cnt] * sum_cos);
#else
        typ    = SORTE(p,i);
        kpot   = charge[typ] * ew_expk[k]
	         * (sinkr[cnt] * sum_sin + coskr[cnt] * sum_cos);
#endif
        POTENG(p,i)  += kpot;

        /* update force */
#ifdef VARCHG
        kforce = CHARGE(p,i) * ew_expk[k] 
                 * (sinkr[cnt] * sum_cos - coskr[cnt] * sum_sin);
#else
        kforce = charge[typ] * ew_expk[k] 
                 * (sinkr[cnt] * sum_cos - coskr[cnt] * sum_sin);
#endif
        KRAFT(p,i,X) += ew_kvek[k].x * kforce;
        KRAFT(p,i,Y) += ew_kvek[k].y * kforce;
        KRAFT(p,i,Z) += ew_kvek[k].z * kforce;
        tmp_virial   += kforce * SPRODX(ORT,p,i,ew_kvek[k]);

        cnt++;
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

#ifdef VARCHG
		charge2 = CHARGE(p,i) * CHARGE(q,i) * coul_eng;
#else
		charge2 = charge[p_typ] * charge[q_typ] * coul_eng;
#endif
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

  int    i, j, k, offx, offy, offz, count, num;
  real   kvek2, vorf1;

  /* we implicitly assume a system of units, in which lengths are
     measured in Angstrom [A], energies in electron volts [eV], 
     and charges in units of the elementary charge e */

  /* coul_eng is already initialized to 14.40 in globals.h  */
  /* this is e^2 / (4*pi*epsilon_0) in [eV A]               */
  /* we need it earlier for analytically defined potentials */

  twopi    = 2.0 * M_PI;

  ew_vorf  = ew_kappa / SQRT( M_PI );
  vorf1    = twopi * coul_eng / volume;

  /* SM ?
  ew_vorf  = 2.0 * ew_kappa / SQRT( M_PI );
  vorf1    = 2.0 * twopi / volume; */

  if (!(ew_kcut > 0)) return;

  ew_nx = (int) (ew_kcut * SQRT( SPROD(box_x,box_x) ) / twopi) + 1;
  ew_ny = (int) (ew_kcut * SQRT( SPROD(box_y,box_y) ) / twopi) + 1;
  ew_nz = (int) (ew_kcut * SQRT( SPROD(box_z,box_z) ) / twopi) + 1;

  if (0==myid)
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

        ew_ivek[ew_totk].x = i + ew_nx;
        ew_ivek[ew_totk].y = j + ew_ny;
        ew_ivek[ew_totk].z = k + ew_nz;

        ew_expk[ew_totk] = vorf1 * exp( -kvek2 / (4.0*SQR(ew_kappa)) ) / kvek2;

        ew_totk++;

      }

  if (0==myid)
    printf("EWALD: %d k-vectors\n", ew_totk); 

  /* Allocate memory for exp(ikr) */
  ew_dx = 2 * ew_nx + 1;
  ew_dy = 2 * ew_ny + 1;
  ew_dz = 2 * ew_nz + 1;
  coskx = (real *) malloc( natoms * ew_dx * sizeof(real));
  sinkx = (real *) malloc( natoms * ew_dx * sizeof(real));
  cosky = (real *) malloc( natoms * ew_dy * sizeof(real));
  sinky = (real *) malloc( natoms * ew_dy * sizeof(real));
  coskz = (real *) malloc( natoms * ew_dz * sizeof(real));
  sinkz = (real *) malloc( natoms * ew_dz * sizeof(real));
  coskr = (real *) malloc( natoms         * sizeof(real));
  sinkr = (real *) malloc( natoms         * sizeof(real));

  if( coskx == NULL || sinkx == NULL || cosky == NULL || sinky == NULL
      || coskz == NULL || sinkz == NULL || coskr == NULL || sinkr == NULL )
    error("EWALD: Cannot allocate memory for exp(ikr)");

  /* Position independent initializations */
  offx = ew_nx * natoms; 
  offy = ew_ny * natoms; 
  offz = ew_nz * natoms;
  for (i=0; i<natoms; i++) {
    coskx[offx+i] = 1.0;
    sinkx[offx+i] = 0.0;
    cosky[offy+i] = 1.0;
    sinky[offy+i] = 0.0;
    coskz[offz+i] = 1.0;
    sinkz[offz+i] = 0.0;
  }
}
