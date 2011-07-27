
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
#include <fcs.h>

fcs_float *pos=NULL, *chg=NULL, *field=NULL, *pot=NULL; 
int       nloc, nloc_max=0;
FCS       handle=NULL;
FCSOutput output=NULL;
 
#define ASSERT_FCS(err) \
  do { \
    if(err) { \
      fcsResult_printResult(err); MPI_Finalize(); exit(-1); \
    } \
  } while (0)

/******************************************************************************
*
* clear_forces
*
******************************************************************************/

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
* pack_fcs
*
******************************************************************************/

void pack_fcs(void) {

  int k, n, m, i;

  /* (re-)allocate data structures if necessary */ 
  nloc = 0;
  for (k=0; k<NCELLS; ++k) nloc += CELLPTR(k)->n;
  if (nloc_max < nloc) {
    nloc_max = (int) (1.1 * nloc);
    free(pos);   pos   = (fcs_float*) malloc(DIM * nloc_max*sizeof(fcs_float));
    free(chg);   chg   = (fcs_float*) malloc(      nloc_max*sizeof(fcs_float));
    free(field); field = (fcs_float*) malloc(DIM * nloc_max*sizeof(fcs_float));
    free(pot);   pot   = (fcs_float*) malloc(      nloc_max*sizeof(fcs_float));
    if ((NULL==pos) || (NULL==chg) || (NULL==field) || (NULL==pot)) 
      error("Could not allocate fcs data");
  }

  /* if we do only fcs, we have to clear forces and energies */
#ifdef PAIR
  if (!have_potfile && !have_pre_pot) 
#endif
    clear_forces();

  /* collect data from cell array */
  n=0; m=0;
  for (k=0; k<NCELLS; ++k) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; ++i) { 
      pos[n++] = ORT(p,i,X); 
      pos[n++] = ORT(p,i,Y); 
      pos[n++] = ORT(p,i,Z); 
      chg[m++] = CHARGE(p,i);
    }
  }
}

/******************************************************************************
*
* unpack_fcs
*
******************************************************************************/

void unpack_fcs(FCSOutput output) {

  fcs_float *vir;
  real pot1, pot2, e, c, sum=0.0;
  int n, m, k, i;

  /* extract output and distribute it to cell array */
  vir = fcsOutput_getVirial(output);
  n=0; m=0; pot1=0.0;
  for (k=0; k<NCELLS; ++k) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; ++i) { 
      c = CHARGE(p,i) * coul_eng;
      KRAFT(p,i,X) += field[n++] * c; 
      KRAFT(p,i,Y) += field[n++] * c; 
      KRAFT(p,i,Z) += field[n++] * c;
      e = pot[m++] * c * 0.5;
      POTENG(p,i)  += e;
      pot1         += e;
    }
  }

  /* unpack virial */
#ifdef P_AXIAL
  vir_xx += vir[0];
  vir_yy += vir[4];
  vir_zz += vir[8];
#else
  virial += vir[0] + vir[4] + vir[8];
#endif
#ifdef STRESS_TENS
  if (do_press_calc) {
    /* distribute virial tensor evenly on all atoms */
    sym_tensor pp;
    pp.xx = vir[0] / natoms;
    pp.yy = vir[4] / natoms;
    pp.zz = vir[8] / natoms;
    pp.yz = (vir[5]+vir[7]) / (2*natoms);
    pp.zx = (vir[2]+vir[6]) / (2*natoms);
    pp.xy = (vir[1]+vir[3]) / (2*natoms);
    for (k=0; k<NCELLS; ++k) {
      cell *p = CELLPTR(k);
      for (i=0; i<p->n; ++i) { 
        PRESSTENS(p,i,xx) += pp.xx;
        PRESSTENS(p,i,yy) += pp.yy;
        PRESSTENS(p,i,zz) += pp.zz;
        PRESSTENS(p,i,yz) += pp.yz;
        PRESSTENS(p,i,zx) += pp.zx;
        PRESSTENS(p,i,xy) += pp.xy;
      }
    }
  }
#endif

  /* sum up potential energy */
#ifdef MPI
  MPI_Allreduce( &pot1, &pot2, 1, MPI_DOUBLE, MPI_SUM, cpugrid);
  tot_pot_energy += pot2;
#else
  tot_pot_energy += pot1;
#endif
}

/******************************************************************************
*
* init_fcs
*
******************************************************************************/

void init_fcs(void) {

  FCSResult result;
  fcs_int srf = 1, cflag = 0;
  char *method;

  fcs_int   pbc [3] = { pbc_dirs.x, pbc_dirs.y, pbc_dirs.z };
  fcs_float BoxX[3] = { box_x.x, box_x.y, box_x.z };
  fcs_float BoxY[3] = { box_y.x, box_y.y, box_y.z };
  fcs_float BoxZ[3] = { box_z.x, box_z.y, box_z.z };

  switch (fcs_method) {
    case FCS_METH_PEPC:   method = "PEPC";   break;
    case FCS_METH_FMM:    method = "FMM";    break;
    case FCS_METH_PP3MG:  method = "PP3MG";  break;
    case FCS_METH_VMG:    method = "VMG";    break;
    case FCS_METH_P3M:    method = "P3M";    srf = fcs_rcut > 0 ? 0 : 1; break;
    case FCS_METH_MEMD:   method = "MEMD";   break;
    case FCS_METH_NFFT:   method = "NFFT";   break;
    case FCS_METH_DIRECT: method = "DIRECT"; break;
  }

  /* initialize handle and set common parameters */
  result = fcs_init(&handle, method, cpugrid); 
  ASSERT_FCS(result);
  result = fcs_common_set(handle, srf, BoxX, BoxY, BoxZ, pbc, natoms, natoms);
  ASSERT_FCS(result);
  result = fcsOutput_create(&output);
  ASSERT_FCS(result);

  /* set method specific parameters */
  switch (fcs_method) {
#ifdef FCS_ENABLE_PEPC
    case FCS_METH_PEPC:
      result = fcs_PEPC_setup(handle, (fcs_float)fcs_pepc_eps, 
           (fcs_float)fcs_pepc_theta, (fcs_int)fcs_debug_level );
      ASSERT_FCS(result);
      break;
#endif
#ifdef FCS_ENABLE_FMM
    case FCS_METH_FMM:
      result = fcs_FMM_setup(handle, (fcs_int)fcs_fmm_absrel, 
           (fcs_float)fcs_fmm_deltaE, (fcs_int)fcs_fmm_dcorr);
      ASSERT_FCS(result);
      break;
#endif
#ifdef FCS_ENABLE_PP3MG_PMG
      /* PP3MG */
#endif
#ifdef FCS_ENABLE_VMG
      /* VMG */
#endif
#ifdef FCS_ENABLE_P3M
    case FCS_METH_P3M:
      fcs_P3M_set_r_cut(handle, (fcs_float)fcs_rcut);
      fcs_P3M_set_required_accuracy(handle, (fcs_float)fcs_p3m_accuracy);
      //fcs_P3M_set_mesh(handle, 64);
      //fcs_P3M_set_cao(handle, 7);
      //fcs_P3M_set_alpha(handle, 2.5);
      break;
#endif
#ifdef FCS_ENABLE_MEMD
      /* MEMD */
#endif
#ifdef FCS_ENABLE_NFFT
      /* NFFT */
#endif
#ifdef FCS_ENABLE_DIRECT
    case FCS_METH_DIRECT:
      /* no need to do anything */
      break;
#endif
    default: 
      error("FCS method unknown or not implemented"); 
      break;
  }
  pack_fcs();
  result = fcs_tune(handle, nloc, nloc_max, pos, chg, cflag);
  ASSERT_FCS(result);

  /* add near-field potential, after fcs_tune */
  if (0==srf) fcs_update_pottab();

}

/******************************************************************************
*
* calc_forces_fcs
*
******************************************************************************/

void calc_forces_fcs(void) {
  FCSResult result;
  fcs_int flag=0;  /* component flag - what for? */
  pack_fcs();
  result = fcs_run(handle, nloc, nloc_max, pos, chg, field, pot, flag, output);
  ASSERT_FCS(result);
  unpack_fcs(output);
}


/******************************************************************************
*
* pair_int_fcs (near field)
*
******************************************************************************/

void fcs_pair_int(real *pot, real *grad, real r2)
{
  fcs_float r, fcs_pot, fcs_grad;
  FCSResult result;
  r = SQRT(r2);
  result = fcs_compute_near(handle, r, &fcs_pot, &fcs_grad);
  ASSERT_FCS(result);
  *pot  = fcs_pot;
  /* return (1/r)*dV/dr as derivative */
  *grad = fcs_grad / r; 
}

/******************************************************************************
*
* add FCS near-field to pair_pot
*
******************************************************************************/

void fcs_update_pottab(void) {

  int i, j, k, col;
  real fcs_shift, fcs_fshift, r2, pot, grad, fcs_r2cut = SQR(fcs_rcut);
  pot_table_t *pt = &pair_pot;

  /* FCS near-field shifts */
  fcs_pair_int(&fcs_shift, &fcs_fshift, fcs_r2cut);
  if (0==myid) {
    printf("FCS near-field shifted by: %g , %g\n", fcs_shift, fcs_fshift);
  }

  /* add near-field to potential table */
  for (i=0; i<ntypes; i++)
    for (j=0; j<ntypes; j++) {
      col = i * ntypes + j;
      for (k=0; k < pt->len[col]; k++) {
        r2 = pt->begin[col] + k * pt->step[col];
        if (r2 > fcs_r2cut) continue;
        fcs_pair_int(&pot, &grad, r2);
        pot -= fcs_shift;
        pot -= 0.5 * fcs_fshift * (r2 - fcs_r2cut);
        pot *= coul_eng * charge[i] * charge[j];
        *PTR_2D(pt->table, k, col, pt->maxsteps, pt->ncols) = pot;
      }
    }
}
