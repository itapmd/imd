
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
* imd_main_2d.c -- main loop in two dimensions
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"


/******************************************************************************
*
*  main_loop
*
******************************************************************************/

void main_loop(void)
{
  int i,j,k;
  cell *p;
  real tmp_pot_energy;
  real tmp_kin_energy;
  real dtemp;
  vektor d_pressure;
#if defined(CORRELATE) || defined(MSQD)
  int ref_step = correl_start;
#endif
#ifdef FBC 
#ifdef MIK
  int nofbcsteps=0;
#endif
  int l;
  vektor nullv={0.0,0.0};
  vektor temp_df;
  vektor *fbc_df;
  fbc_df = (vektor *) malloc(vtypes*DIM*sizeof(real));
  if (NULL==fbc_df)
    error("Can't allocate memory for fbc_df\n");
#endif 

  /* initialize temperature if necessary */
  if (0==restart) { 
    if (do_maxwell) maxwell(temperature);
    do_maxwell=0;
  }

#if defined(FRAC) || defined(FTG) 
  if (0==myid) {
      printf( "Strain rate is  %1.10f\n", dotepsilon0 );
      printf( "Damping mode  %d\n", dampingmode );
      printf( "Damping prefactor is %1.10f\n\n", gamma_bar );
  }
#endif

  if (0==myid) printf( "Starting simulation %d\n", simulation );

#if defined(AND) || defined(NVT) || defined(NPT) || defined(STM) || defined(FRAC)
  dtemp = (end_temp - temperature) / (steps_max - steps_min);
#endif

#ifdef FBC
#ifndef MIK
/* dynamic loading, increment linearly each timestep */
  for (l=0;l<vtypes;l++){
    temp_df.x = (((fbc_endforces+l)->x) - ((fbc_beginforces+l)->x))/(steps_max - steps_min);
    temp_df.y = (((fbc_endforces+l)->y) - ((fbc_beginforces+l)->y))/(steps_max - steps_min);
    
    *(fbc_df+l) = temp_df;
  }
#endif
#endif

#ifdef NVX
  dtemp = (dTemp_end - dTemp_start)/(steps_max - steps_min);
  tran_Tleft  = temperature + dTemp_start;
  tran_Tright = temperature - dTemp_start;
#endif

#ifdef NPT
  d_pressure.x = (pressure_end.x - pressure_ext.x) / (steps_max - steps_min);
  d_pressure.y = (pressure_end.y - pressure_ext.y) / (steps_max - steps_min);
  calc_dyn_pressure();
  if (isq_tau_xi==0.0) {
    xi.x = 0.0;
    xi.y = 0.0;
  }
#endif

#if defined(CORRELATE) || defined(MSQD)
  init_correl(ncorr_rmax,ncorr_tmax);
#endif

#ifdef ATDIST
  init_atoms_dist();
#endif

#ifdef DEFORM
  deform_int = 0;
#endif

  for (steps=steps_min; steps <= steps_max; ++steps) { 

#ifdef STRESS_TENS
    do_press_calc = (((eng_interval  >0) && (0==steps%eng_interval)) ||
                     ((press_interval>0) && (0==steps%eng_interval)));
#endif

#ifdef EPITAX
    for (i=0; i<ntypes; ++i ) {
      if ( (steps >= epitax_startstep) && (steps <= epitax_maxsteps) ) {  
	if ( (epitax_rate[i] != 0)&&(steps-steps_min>0) && ((steps-steps_min)%epitax_rate[i]) == 0 ) {
	  delete_atoms(); 
	  create_atom(i, epitax_mass[i], epitax_temp[i]);
	}
	else if ( (epitax_rate[i] != 0) && (steps > epitax_maxsteps) && ((steps-steps_min)%epitax_rate[i]) == 0 ) 
	  delete_atoms();
      }
    }
#endif

#ifdef FBC
#ifdef MIK  
    /* just increment the force if under threshold of e_kin or after 
       waitsteps and after annelasteps */
    temp_df.x = 0.0;
    temp_df.y = 0.0;
    for (l=0;l<vtypes;l++) *(fbc_df+l) = temp_df;

    if (steps > fbc_annealsteps) {
      nofbcsteps++; 
      if ((2.0*tot_kin_energy/nactive < fbc_ekin_threshold) ||
          (nofbcsteps==fbc_waitsteps)) {
        nofbcsteps=0;
        for (l=0;l<vtypes;l++) *(fbc_df+l) = *(fbc_dforces+l);
      }
    }
#endif /* MIK */
#endif /* FBC */

#ifdef HOMDEF
    if ((exp_interval > 0) && (0 == steps%exp_interval)) expand_sample();
    if ((hom_interval > 0) && (0 == steps%hom_interval)) shear_sample();
    if ((lindef_interval > 0) && (0 == steps%lindef_interval)) lin_deform();
#endif
#ifdef DEFORM
    if (steps > annealsteps) {
      deform_int++;
      if ((2.0*tot_kin_energy/nactive < ekin_threshold) || 
          (deform_int==max_deform_int)) {
        deform_sample();
        deform_int=0;
      }
    }
#endif

#ifdef FBC
    /* fbc_forces is already initialised with beginforces */
    for (l=0;l<vtypes;l++) {
      /* fbc_df=0 if MIK && ekin> ekin_threshold */ 
      (fbc_forces+l)->x += (fbc_df+l)->x; 
      (fbc_forces+l)->y += (fbc_df+l)->y;
    } 
#endif

    calc_forces(steps); 

#ifdef EPITAX
    if (steps == steps_min) {
      calc_poteng_min();
      if (0 == myid) 
	printf("EPITAX: Minimal potential energy = %lf\n", epitax_poteng_min);
    }
#endif

#ifdef DISLOC
    if ((steps==reset_Epot_step) && (calc_Epot_ref==1)) reset_Epot_ref();
#endif
#ifdef AVPOS
    if ((steps == steps_min) || (steps == avpos_start)) {
       update_avpos();
    }
#endif
#ifdef ATDIST
    if ((atoms_dist_int>0) && (0==steps % atoms_dist_int) &&
        (steps>=atoms_dist_start) && (steps<=atoms_dist_end)) 
      update_atoms_dist();
#endif

#if defined(CORRELATE) || defined(MSQD)
    if ((steps >= correl_start) && ((steps < correl_end) || (correl_end==0))) {
      int istep = steps - correl_start;
      if (istep % correl_ts == 0) correlate(steps,ref_step,istep/correl_ts);
      if ((correl_int != 0) && (steps-ref_step+1 >= correl_int)) 
        ref_step += correl_int;
    }
#endif

#ifdef GLOK 
    /* "globale konvergenz": set impulses=0 if P*F <0 (global vectors) */
    if (PxF<0.0)
      for (k=0; k<ncells; ++k) {
	p = cell_array + CELLS(k);
	for (i=0; i<p->n; ++i) {
	  p->impuls X(i) = 0.0;
	  p->impuls Y(i) = 0.0;
	}
      }
#endif

    move_atoms(); 

#if defined(EPITAX) && !defined(NVE)
    /* beam atoms are always integrated by NVE */
    move_atoms_nve();
#endif 

#if defined(AND) || defined(NVT) || defined(NPT) || defined(STM) || defined(FRAC)
    if ((steps==steps_min) && (use_curr_temp==1)) {
#ifdef STM
      temperature = 2 * tot_kin_energy / (nactive-n_stadium);
#else
      temperature = 2 * tot_kin_energy / nactive;
#endif
      dtemp = (end_temp - temperature) / (steps_max - steps_min);
      use_curr_temp = 0;
    }
#endif

#ifdef NPT_iso
    if ((steps==steps_min) && (ensemble==ENS_NPT_ISO) && 
        (use_curr_pressure==1)) {
      pressure_ext.x = pressure;
      d_pressure.x = (pressure_end.x-pressure_ext.x) / (steps_max-steps_min);
      d_pressure.y = (pressure_end.y-pressure_ext.y) / (steps_max-steps_min);
      use_curr_pressure = 0;
    }
#endif

#ifdef NPT_axial
    if ((steps==steps_min) && (ensemble==ENS_NPT_AXIAL) && 
        (use_curr_pressure==1)) {
      pressure_ext.x = stress_x;
      pressure_ext.y = stress_y;
      d_pressure.x = (pressure_end.x-pressure_ext.x) / (steps_max-steps_min);
      d_pressure.y = (pressure_end.y-pressure_ext.y) / (steps_max-steps_min);
      use_curr_pressure = 0;
    }
#endif

#if defined(AND) || defined(NVT) || defined(NPT) || defined(STM) || defined(FRAC) 
    temperature += dtemp;
#endif

#ifdef NVX
    tran_Tleft   += dtemp;
    tran_Tright  -= dtemp;
#endif

#ifdef NPT
    pressure_ext.x += d_pressure.x;
    pressure_ext.y += d_pressure.y;
#endif

    /* Periodic I/O */
#ifdef TIMING
    imd_start_timer(&time_io);
#endif
    if ((rep_interval > 0) && (0 == steps%rep_interval)) write_config(steps);
    if ((eng_interval > 0) && (0 == steps%eng_interval) && (0==myid)) 
       write_eng_file(steps);
    if ((dis_interval > 0) && (0 == steps%dis_interval)) write_distrib(steps);
    if ((pic_interval > 0) && (0 == steps%pic_interval)) write_pictures(steps);

#ifdef EFILTER  /* just print atoms in an energy-window */ 
    if ((efrep_interval > 0) && (0 == steps%efrep_interval)) 
       write_config_select(steps/efrep_interval, "ef",
                           write_atoms_ef, write_header_ef);
#endif
#ifdef DISLOC
    if (steps == up_ort_ref) update_ort_ref();
    if ((dem_interval > 0) && (0 == steps%dem_interval)) 
       write_config_select(steps, "dem", write_atoms_dem, write_header_dem);
    if ((dsp_interval > up_ort_ref) && (0 == steps%dsp_interval)) 
       write_config_select(steps, "dsp", write_atoms_dsp, write_header_dsp);
#endif
#ifdef AVPOS
    if ( steps <= avpos_end ) {
	if ((avpos_res > 0) && (0 == (steps - avpos_start) % avpos_res) && steps > avpos_start)
	    add_position();
	if ((avpos_int > 0) && (0 == (steps - avpos_start) % avpos_int) && steps > avpos_start) {
	    write_config_select((steps-avpos_start)/avpos_int,"avp",
				write_atoms_avp,write_header_avp);
	    update_avpos();
	}
    }
#endif

#ifdef TRANSPORT
    if ((tran_interval > 0) && (steps > 0) && (0 == steps%tran_interval)) 
       write_temp_dist(steps);
#endif
#ifdef RNEMD
    if ((exch_interval > 0) && (0 == steps%exch_interval)) 
       rnemd_heat_exchange();
#endif

#ifdef STRESS_TENS
#ifdef SHOCK
    if ((press_interval > 0) && (0 == steps%press_interval)) 
       write_press_dist_shock(steps);
#else
    if ((press_interval > 0) && (0 == steps%press_interval)) {
      if (0==press_dim.x) 
        write_config_select(steps/press_interval, "press",
                            write_atoms_press, write_header_press);
      else 
        write_press_dist(steps);
    }
#endif
#endif

#ifdef USE_SOCKETS
    if ((socket_int>0) && (0==steps%socket_int)) check_socket(steps);
#endif

#ifdef TIMING
    imd_stop_timer(&time_io);
#endif

    do_boundaries();    
    fix_cells();

#ifdef ATDIST
    if (steps==atoms_dist_end) write_atoms_dist();
#endif
#ifdef MSQD
    if ((correl_end >0) && (steps==correl_end) || 
        (correl_end==0) && (steps==steps_max))   
      write_config_select(0, "sqd", write_atoms_sqd, write_header_sqd);
#endif
  }

  /* clean up the current phase, and clear restart flag */
  restart=0;
  if (0==myid) {
    write_itr_file(-1, steps_max);
    printf( "End of simulation %d\n", simulation );
  }

}

/******************************************************************************
*
*  do_boundaries
*
*  Apply periodic boundaries to all atoms
*  Could change so that only cells on surface to some work
*
******************************************************************************/

void do_boundaries(void)
{
  int k;

  /* for each cell in bulk */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (k=0; k<ncells; ++k) {

    int l,i;
    cell *p;

    p = cell_array + CELLS(k);

    /* PBC in x direction */
    if (pbc_dirs.x==1)
    for (l=0; l<p->n; ++l) {
      i = -FLOOR(SPRODX(p->ort,l,tbox_x));
      p->ort X(l)    += i * box_x.x;
      p->ort Y(l)    += i * box_x.y;
#ifdef MSQD
      p->refpos X(l) += i * box_x.x;
      p->refpos Y(l) += i * box_x.y;
#ifdef AVPOS
      p->sheet X(l)  -= i * box_x.x;
      p->sheet Y(l)  -= i * box_x.y;
#endif
#endif
    }

    /* PBC in y direction */
    if (pbc_dirs.y==1)
    for (l=0; l<p->n; ++l) {
      i = -FLOOR(SPRODX(p->ort,l,tbox_y));
      p->ort X(l)    += i * box_y.x;
      p->ort Y(l)    += i * box_y.y;
#ifdef MSQD
      p->refpos X(l) += i * box_y.x;
      p->refpos Y(l) += i * box_y.y;
#endif
#ifdef AVPOS
      p->sheet X(l)  -= i * box_y.x;
      p->sheet Y(l)  -= i * box_y.y;
#endif
    }

  }
}

#ifdef STRESS_TENS

/******************************************************************************
*
*  calc_tot_presstens
*
******************************************************************************/

void calc_tot_presstens(void)
{
  int i;

  real tmp_presstens1[3], tmp_presstens2[3];

  tot_presstens.xx = 0.0; 
  tot_presstens.yy = 0.0; 
  tot_presstens.xy = 0.0;

  /*sum up total pressure tensor */
  for (i=0; i<ncells; ++i) {
    int j;
    cell *p;
    p = cell_array + CELLS(i);
    for (j=0; j<p->n; ++j) {
      tot_presstens.xx += p->presstens[j].xx;
      tot_presstens.yy += p->presstens[j].yy;
      tot_presstens.xy += p->presstens[j].xy;
    }
  }

#ifdef MPI
  /* add up results from different CPUs */
  tmp_presstens1[0] = tot_presstens.xx; 
  tmp_presstens1[1] = tot_presstens.yy; 
  tmp_presstens1[2] = tot_presstens.xy;

  MPI_Allreduce( tmp_presstens1, tmp_presstens2, 3, REAL, MPI_SUM, cpugrid);

  tot_presstens.xx = tmp_presstens2[0];
  tot_presstens.yy = tmp_presstens2[1]; 
  tot_presstens.xy = tmp_presstens2[2]; 
#endif /* MPI */

}

#endif/* STRESS_TENS */

