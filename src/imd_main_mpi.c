
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
* imd_main_mpi_3d.c -- main loop, mpi specific part, three dimensions
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/******************************************************************************
*
* calc_forces 
*
* The forces of the atoms are calulated here. To achive this, atoms on
* the surface of a cpu are exchanged with the neigbours.
*
* If AR is defined, we use actio=reactio even across CPUs, otherwise we don't
*
* The force calculation is split into those steps:
*
* i)   send atoms positions of cells on surface neighbours, 
*      receive atom positions from neigbours
* ii)  zero forces on all cells (local and buffer)
* iii) calculate forces in local cells, use lower half of neigbours 
*      for each cell and use actio==reactio
* iv)  calculate forces also for upper half of neighbours for all cells
*      that are on the upper surface
* iv)  or send forces back and add them
*
******************************************************************************/

void calc_forces(int steps)
{
  int n, k;
  real tmpvec1[8], tmpvec2[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  /* fill the buffer cells */
  if ((steps == steps_min) || (0 == steps % BUFSTEP)) setup_buffers();
  send_cells(copy_cell,pack_cell,unpack_cell);

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
  for (k=0; k<nallcells; ++k) {
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
#ifdef NNBR
      NBANZ(p,i) = 0;
#endif
#ifdef CNA
      if (cna)
	MARK(p,i) = 0;
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

  /* What follows is the standard one-cpu force 
     loop acting on our local data cells */

  /* compute forces for all pairs of cells */
  for (n=0; n<nlists; ++n) {
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime) \
  reduction(+:tot_pot_energy,virial,vir_xx,vir_yy,vir_zz,vir_yz,vir_zx,vir_xy)
#endif
    for (k=0; k<npairs[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n] + k;
      pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
      pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
      pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
      do_forces(cell_array + P->np, cell_array + P->nq, pbc,
                &tot_pot_energy, &virial, &vir_xx, &vir_yy, &vir_zz,
                                          &vir_yz, &vir_zx, &vir_xy);
    }
  }

#ifdef COVALENT
  /* complete neighbor tables for remaining pairs of cells */
  for (n=0; n<nlists; ++n) {
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
#endif
    for (k=npairs[n]; k<npairs2[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n] + k;
      pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
      pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
      pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
      do_neightab(cell_array + P->np, cell_array + P->nq, pbc);
    }
  }

#ifndef CNA
  /* second force loop for covalent systems */
/* does not work correctly - different threads may write to same variables 
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime) \
  reduction(+:tot_pot_energy,virial,vir_xx,vir_yy,vir_zz,vir_yz,vir_zx,vir_xy)
#endif
*/
  for (k=0; k<ncells; ++k) {
    do_forces2(cell_array + CELLS(k),
               &tot_pot_energy, &virial, &vir_xx, &vir_yy, &vir_zz,
                                         &vir_yz, &vir_zx, &vir_xy);
  }
#endif
#endif /* COVALENT */

#ifndef AR
  /* If we don't use actio=reactio accross the cpus, we have do do
     the force loop also on the other half of the neighbours for the 
     cells on the surface of the CPU */

  /* compute forces for remaining pairs of cells */
  for (n=0; n<nlists; ++n) {
/* does not work correctly - different threads may write to same variables 
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
#endif
*/
    for (k=npairs[n]; k<npairs2[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n] + k;
      pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
      pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
      pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
      /* potential energy and virial are already complete;          */
      /* to avoid double counting, we update only the dummy tmpvec2 */
      do_forces(cell_array + P->np, cell_array + P->nq, pbc,
                tmpvec2, tmpvec2+1, tmpvec2+2, tmpvec2+3, tmpvec2+4,
                                    tmpvec2+5, tmpvec2+6, tmpvec2+7);
    }
  }
#endif  /* not AR */

#ifdef EAM2

#ifdef AR
  /* collect host electron density */
  send_forces(add_rho,pack_rho,unpack_add_rho);
#endif
  /* compute embedding energy and its derivative */
  do_embedding_energy();
  /* distribute derivative of embedding energy */
  send_cells(copy_dF,pack_dF,unpack_dF);

  /* second EAM2 loop over all cells pairs */
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

#ifndef AR
  /* If we don't use actio=reactio accross the cpus, we have do do
     the force loop also on the other half of the neighbours for the 
     cells on the surface of the CPU */

  /* compute forces for remaining pairs of cells */
  for (n=0; n<nlists; ++n) {
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
#endif
    for (k=npairs[n]; k<npairs2[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n]+k;
      pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
      pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
      pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
      /* potential energy and virial are already complete;          */
      /* to avoid double counting, we update only the dummy tmpvec2 */
      do_forces_eam2(cell_array + P->np, cell_array + P->nq, pbc,
                     tmpvec2, tmpvec2+1, tmpvec2+2, tmpvec2+3, tmpvec2+4,
                                         tmpvec2+5, tmpvec2+6, tmpvec2+7);
    }
  }
#endif /* not AR */

#endif /* EAM2 */

  /* sum up results of different CPUs */
  tmpvec1[0] = tot_pot_energy;
  tmpvec1[1] = virial;
  tmpvec1[2] = vir_xx;
  tmpvec1[3] = vir_yy;
  tmpvec1[4] = vir_zz;
  tmpvec1[5] = vir_yz;
  tmpvec1[6] = vir_zx;
  tmpvec1[7] = vir_xy;

  MPI_Allreduce( tmpvec1, tmpvec2, 8, REAL, MPI_SUM, cpugrid); 

  tot_pot_energy = tmpvec2[0];
  virial         = tmpvec2[1];
  vir_xx         = tmpvec2[2];
  vir_yy         = tmpvec2[3];
  vir_zz         = tmpvec2[4];
  vir_yz         = tmpvec2[5];
  vir_zx         = tmpvec2[6];
  vir_xy         = tmpvec2[7];

#ifdef AR
  send_forces(add_forces,pack_forces,unpack_forces);
#endif

}

