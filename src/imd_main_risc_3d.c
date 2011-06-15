
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2011 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_main_risc_3d.c -- main loop, risc specific part, three dimensions
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"


/*****************************************************************************
*
*  calc_forces
*
*****************************************************************************/

void calc_forces(int steps)
{
  int n, k;

  /* clear global accumulation variables */
  tot_pot_energy = 0.0;
  virial = 0.0;
  vir_xx = 0.0;
  vir_yy = 0.0;
  vir_zz = 0.0;
  vir_yz = 0.0;
  vir_zx = 0.0;
  vir_xy = 0.0;
  nfc++;

  /* clear per atom accumulation variables */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (k=0; k<ncells; ++k) {
    int  i;
    cell *p;
    p = cell_array + k;
    for (i=0; i<p->n; ++i) {
      KRAFT(p,i,X) = 0.0;
      KRAFT(p,i,Y) = 0.0;
      KRAFT(p,i,Z) = 0.0;
#ifdef UNIAX
      DREH_MOMENT(p,i,X) = 0.0;
      DREH_MOMENT(p,i,Y) = 0.0;
      DREH_MOMENT(p,i,Z) = 0.0;
#endif
#if defined(STRESS_TENS)
      PRESSTENS(p,i,xx) = 0.0;
      PRESSTENS(p,i,yy) = 0.0;
      PRESSTENS(p,i,zz) = 0.0;
      PRESSTENS(p,i,yz) = 0.0;
      PRESSTENS(p,i,zx) = 0.0;
      PRESSTENS(p,i,xy) = 0.0;
#endif     
#ifndef MONOLJ
      POTENG(p,i) = 0.0;
#endif
#ifdef CNA
      if (cna)
	MARK(p,i) = 0;
#endif
#ifdef NNBR
      NBANZ(p,i) = 0;
#endif
#ifdef COVALENT
      NEIGH(p,i)->n = 0;
#endif
#ifdef EAM2
      EAM_RHO(p,i) = 0.0; /* zero host electron density at atom site */
#ifdef EEAM
      EAM_P(p,i) = 0.0; /* zero host electron density at atom site */
#endif
#endif
    }
  }
#ifdef RIGID
  /* clear total forces */
  if ( nsuperatoms>0 ) 
    for(k=0; k<nsuperatoms; k++) {
      superforce[k].x = 0.0;
      superforce[k].y = 0.0;
      superforce[k].z = 0.0;
    }
#endif

#ifdef EWALD
  if (steps==0) {
    ewald_time.total = 0.0;
    imd_start_timer( &ewald_time );
  }
#endif

  /* compute forces for all pairs of cells */
  for (n=0; n<nlists; ++n) {
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime) \
  reduction(+:tot_pot_energy,virial,vir_xx,vir_yy,vir_zz,vir_yz,vir_zx,vir_xy)
#endif
    for (k=0; k<npairs[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n]+k;
      pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
      pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
      pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
      do_forces(cell_array + P->np, cell_array + P->nq, pbc,
                &tot_pot_energy, &virial, &vir_xx, &vir_yy, &vir_zz,
                                          &vir_yz, &vir_zx, &vir_xy);
    }
  }

#ifdef EWALD
  if (steps==0) {
    imd_stop_timer( &ewald_time );
  }
#endif

#ifdef EAM2
  /* compute embedding energy and its derivative */
  do_embedding_energy();

  for (n=0; n<nlists; ++n) {
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime) \
  reduction(+:virial,vir_xx,vir_yy,vir_zz,vir_yz,vir_zx,vir_xy)
#endif
    for (k=0; k<npairs[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n]+k;
      pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
      pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
      pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
      do_forces_eam2(cell_array + P->np, cell_array + P->nq, pbc,
        &virial, &vir_xx, &vir_yy, &vir_zz, &vir_yz, &vir_zx, &vir_xy);
    }
  }
#endif

#if defined(COVALENT) && !defined(CNA)
/* does not work correctly - different threads may write to same variables 
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime) \
  reduction(+:tot_pot_energy,virial,vir_xx,vir_yy,vir_zz,vir_yz,vir_zx,vir_xy)
#endif
*/
  for (k=0; k<ncells; ++k) {
    do_forces2(cell_array+k, &tot_pot_energy, &virial, 
               &vir_xx, &vir_yy, &vir_zz, &vir_yz, &vir_zx, &vir_xy);
  }
#endif

#ifdef EWALD 
  do_forces_ewald(steps);
#endif 

}

/******************************************************************************
*
*  fix_cells
*
*  check if each atom is in the correct cell; 
*  move atoms that have left their cells
*
******************************************************************************/

void fix_cells(void)
{
  int i,j,k,l,clone;
  cell *p, *q;
  ivektor coord, lcoord;

  /* apply periodic boundary conditions */
  do_boundaries();

  /* for each cell in bulk */
  for (i=cellmin.x; i < cellmax.x; ++i)
    for (j=cellmin.y; j < cellmax.y; ++j)
      for (k=cellmin.z; k < cellmax.z; ++k) {

	p = PTR_3D_V(cell_array, i, j, k, cell_dim);
	/*printf(" cell %d %d %d \n",i,j,k);fflush(stdout);*/
	/* loop over atoms in cell */
	l=0;
	while( l<p->n ) {
          coord = cell_coord( ORT(p,l,X), ORT(p,l,Y), ORT(p,l,Z) );
          q = PTR_3D_VV(cell_array,coord,cell_dim);
          /* if it's in the wrong cell, move it to the right cell */
          if (p != q) {
            MOVE_ATOM(q,p,l); 
#ifdef CLONE
            if (l < p->n-nclones)
              for (clone=1; clone<nclones; clone++)
                MOVE_ATOM(q, p, l+clone);
            else /* we are dealing with the last in the stack */
              for (clone=1; clone<nclones; clone++)
                MOVE_ATOM(q, p, l);
#endif
	  }
          else ++l;
	}
      }
}
