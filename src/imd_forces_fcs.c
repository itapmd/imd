
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
* imd_forces_fcs.c -- interface to ScaFaCoS library
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"
#ifdef PEPC
#include <fcs_pepc.h>
#endif
#ifdef PP3MG
#include <fcs_pp3mg.h>
#endif

double *fcs_pos=NULL, *fcs_chg=NULL, *fcs_force=NULL, *fcs_pot=NULL; 
int loc_atoms_max=0;

void fcsAssert(FCSError err) {
  if(fcsError_getErrorCode(err) != FCS_SUCCESS) {
    fcsError_printError(err);
    MPI_Finalize();
    exit(-1);
  }
  fcsError_destroy(&err);
}

void clear_forces(void) {
  int k,i;
  tot_pot_energy = 0.0;
  virial = vir_xx = vir_yy = vir_xy = vir_zz = vir_yz = vir_zx = 0.0;
  for (k=0; k<nallcells; k++) {
    cell *p = cell_array + k;
    for (i=0; i<p->n; i++) {
      KRAFT(p,i,X) = 0.0;
      KRAFT(p,i,Y) = 0.0;
      KRAFT(p,i,Z) = 0.0;
#if defined (STRESS_TENS)
      PRESSTENS(p,i,xx) = 0.0;
      PRESSTENS(p,i,yy) = 0.0;
      PRESSTENS(p,i,xy) = 0.0;
      PRESSTENS(p,i,zz) = 0.0;
      PRESSTENS(p,i,yz) = 0.0;
      PRESSTENS(p,i,zx) = 0.0;
#endif
      POTENG(p,i) = 0.0;
    }
  }
}


/******************************************************************************
*
* calc_forces_pepc
*
******************************************************************************/

#ifdef PEPC

void calc_forces_pepc(void) {

  FCSError err;
  FCSInput input_handle;
  FCSOutput output_handle;
  FCSParameter param_handle = NULL;
  FCSMethodProperties props_handle = NULL;
  str255 paramstring;
  double box_length[3] = { box_x.x,    box_y.y,    box_z.z    };
  char   pbc_flags [3] = { pbc_dirs.x, pbc_dirs.y, pbc_dirs.z };
  double e, pot1, pot2;
  int loc_atoms, k, n, m, i;
  cell *p;
  FILE *out;
  str255 fname;

  sprintf(paramstring, "%s,%s", PEPC_THETA, PEPC_EPS);

  /* (re-)allocate data structures if necessary */ 
  loc_atoms = 0;
  for (k=0; k<NCELLS; ++k) loc_atoms += CELLPTR(k)->n;
  if (loc_atoms_max < loc_atoms) {
    loc_atoms_max = (int) (1.1 * loc_atoms);
    free(fcs_pos); 
    fcs_pos   = (double *) malloc( DIM * loc_atoms_max * sizeof(double) );
    free(fcs_chg); 
    fcs_chg   = (double *) malloc(       loc_atoms_max * sizeof(double) );
  }

  /* if we do only fcs, we have to clear forces and energies */
#ifdef PAIR
  if (!have_potfile && !have_pre_pot) 
#endif
    clear_forces();

  /* collect data from cell array */
  n=0; m=0;
  for (k=0; k<NCELLS; ++k) {
    p = CELLPTR(k);
    for (i=0; i<p->n; ++i) { 
      fcs_pos[n++] = ORT(p,i,X); 
      fcs_pos[n++] = ORT(p,i,Y); 
      fcs_pos[n++] = ORT(p,i,Z); 
      fcs_chg[m++] = CHARGE(p,i);
    }
  }

  /* create parameter handle - move to init function */
  /* we currently can't use loc_atoms_max */
  fcsAssert(fcsParameter_create(&param_handle, natoms, loc_atoms, loc_atoms,
                                DIM, box_length, pbc_flags, cpugrid, 1));

  /* initialize pepc method - move to init function */
  fcsAssert(fcsinit_pepc(param_handle, &props_handle));


  /* create fcs input structure */
  fcsAssert(fcsInput_create( &input_handle, natoms, loc_atoms,
                             fcs_pos, DIM, fcs_chg));

  /* create fcs output structure */
  fcsAssert(fcsOutput_create(&output_handle, loc_atoms, DIM));

  /* run pepc via fcs interface */
  fcsAssert(fcsrun_pepc_with_opt_param(props_handle,input_handle,output_handle, 
                                   paramstring, fcs_pepc_theta, fcs_pepc_eps));
            
  /* extract output and distribute it to cell array */
  fcs_force = fcsOutput_getForces(output_handle);
  fcs_pot   = fcsOutput_getPotentials(output_handle);
  n=0; m=0; pot1=0.0;
  for (k=0; k<NCELLS; ++k) {
    p = CELLPTR(k);
    for (i=0; i<p->n; ++i) { 
      KRAFT(p,i,X) += fcs_force[n++] * coul_eng; 
      KRAFT(p,i,Y) += fcs_force[n++] * coul_eng; 
      KRAFT(p,i,Z) += fcs_force[n++] * coul_eng;
      e = CHARGE(p,i) * fcs_pot[m++] * coul_eng * 0.5;
      POTENG(p,i)  += e;
      pot1         += e;
    }
  }
  fcsAssert(fcsInput_destroy(&input_handle));
  fcsAssert(fcsOutput_destroy(&output_handle));

  /* deallocate stuff */
  fcsAssert(fcsfree_pepc());
  fcsMethodProperties_destroy(&props_handle);
  fcsParameter_destroy(&param_handle);

#ifdef MPI
  MPI_Allreduce( &pot1, &pot2, 1, MPI_DOUBLE, MPI_SUM, cpugrid);
  tot_pot_energy += pot2;
#else
  tot_pot_energy += pot1;
#endif

}

#endif

/******************************************************************************
*
* calc_forces_pp3mg
*
******************************************************************************/

#ifdef PP3MG

void calc_forces_pp3mg(void) {

  FCSError err;
  FCSInput input_handle;
  FCSOutput output_handle;
  FCSParameter param_handle = NULL;
  FCSMethodProperties props_handle = NULL;
  str255 pp3mg_paramstring;
  int pp3mg_ghosts  =   4;
  int pp3mg_cells_x = 128;
  int pp3mg_cells_y = 128;
  int pp3mg_cells_z = 128;
  int    mpi_dims  [3] = { cpu_dim.x,  cpu_dim.y,  cpu_dim.z  };
  double box_length[3] = { box_x.x,    box_y.y,    box_z.z    };
  char   pbc_flags [3] = { pbc_dirs.x, pbc_dirs.y, pbc_dirs.z };
  double e, pot1, pot2;
  int loc_atoms, k, n, m, i, max_part = 70000;
  cell *p;
  FILE *out;
  str255 fname;

  /* further parameters:
  PP3MG_PERIODIC periodic boundary conditions? (TRUE|FALSE)
  PP3MG_DEGREE  Degree of the interpolation spline 
  PP3MG_MAX_ITERATIONS Maximum number of iterations 
  PP3MG_ERR_BOUND required accuracy (relative) */

  /* parameters not passed are max_iterations, nu1, nu2, periodic, 
     pol_degree, err_bound */
  sprintf(pp3mg_paramstring, "%s,%s,%s,%s,%s,%s", 
          PP3MG_MPI_DIMS, PP3MG_GHOST_CELLS, PP3MG_MAX_PARTICLES,
          PP3MG_CELLS_X, PP3MG_CELLS_Y, PP3MG_CELLS_Z);

  /* (re-)allocate data structures if necessary */ 
  loc_atoms = 0;
  for (k=0; k<NCELLS; ++k) loc_atoms += CELLPTR(k)->n;
  if (loc_atoms_max < loc_atoms) {
    loc_atoms_max = (int) (1.1 * loc_atoms);
    free(fcs_pos); 
    fcs_pos   = (double *) malloc( DIM * loc_atoms_max * sizeof(double) );
    free(fcs_chg); 
    fcs_chg   = (double *) malloc(       loc_atoms_max * sizeof(double) );
  }

  /* if we do only fcs, we have to clear forces and energies */
#ifdef PAIR
  if (!have_potfile && !have_pre_pot) 
#endif
    clear_forces();

  /* collect data from cell array */
  n=0; m=0;
  for (k=0; k<NCELLS; ++k) {
    p = CELLPTR(k);
    for (i=0; i<p->n; ++i) { 
      fcs_pos[n++] = ORT(p,i,X); 
      fcs_pos[n++] = ORT(p,i,Y); 
      fcs_pos[n++] = ORT(p,i,Z); 
      fcs_chg[m++] = CHARGE(p,i);
    }
  }

  /* create parameter handle - move to init function */
  fcsAssert(fcsParameter_create(&param_handle, natoms, loc_atoms, loc_atoms,
                                DIM, box_length, pbc_flags, cpugrid, 1));

  /* create fcs input structure */
  fcsAssert(fcsInput_create(&input_handle,natoms,loc_atoms,fcs_pos,DIM,fcs_chg));

  /* create fcs output structure */
  fcsAssert(fcsOutput_create(&output_handle, loc_atoms, DIM));

  /* initialize pp3mg method - move to init function? */
  fcsAssert(fcsinit_pp3mg_with_opt_param(param_handle, &props_handle,
                                         "PP3MG_MAX_PARTICLES", max_part));
  /*
    pp3mg_paramstring, mpi_dims, pp3mg_ghosts, max_part, 
    pp3mg_cells_x, pp3mg_cells_y, pp3mg_cells_z);
  */

  /* run ppp3mg via fcs interface */
  fcsAssert(err = fcsrun_pp3mg(props_handle, input_handle, output_handle));

  /* extract output and distribute it to cell array */
  fcs_force = fcsOutput_getForces(output_handle);
  fcs_pot   = fcsOutput_getPotentials(output_handle);
  n=0; m=0; pot1=0.0;
  for (k=0; k<NCELLS; ++k) {
    p = CELLPTR(k);
    for (i=0; i<p->n; ++i) { 
      KRAFT(p,i,X) += fcs_force[n++] * coul_eng; 
      KRAFT(p,i,Y) += fcs_force[n++] * coul_eng; 
      KRAFT(p,i,Z) += fcs_force[n++] * coul_eng;
      e = CHARGE(p,i) * fcs_pot[m++] * coul_eng;
      POTENG(p,i)  += e;
      pot1         += e;
    }
  }
  fcsAssert(fcsInput_destroy(&input_handle));
  fcsAssert(fcsOutput_destroy(&output_handle));

  /* deallocate stuff */
  fcsAssert(fcsfree_pp3mg());
  fcsAssert(fcsMethodProperties_destroy(&props_handle));
  fcsAssert(fcsParameter_destroy(&param_handle));

#ifdef MPI
  MPI_Allreduce( &pot1, &pot2, 1, MPI_DOUBLE, MPI_SUM, cpugrid);
  tot_pot_energy += pot2;
#else
  tot_pot_energy += pot1;
#endif

}

#endif

