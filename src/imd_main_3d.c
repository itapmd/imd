
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2008 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_main_3d.c -- main loop, used for both two and three dimensions
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/*****************************************************************************
*
*  main_loop
*
*****************************************************************************/

int main_loop(int simulation)
{
  int  finished = 0;
  int  i, j, k, l;
  int  steps_diff = steps_max - steps_min;
  int  deform_int = 0, do_fbc_incr = 0;
  int  have_fbc_incr = 0;
  real dtemp, dshock_speed;
  vektor d_pressure, *fbc_df;
  real tmpvec1[DIM], tmpvec2[DIM];
  char tmp_str[9];
#ifdef GLOK
  int glok_start = steps_min; 
#endif
#ifdef ACG
  acg_alpha = acg_init_alpha;
#endif
#ifdef RELAX
  is_relaxed = 0;
#endif

#ifdef SHOCK
  /* compute speed increase */
  if (shock_mode == 3) { 
    dshock_speed=0;
    if (shock_incr>0) {
      dshock_speed=shock_speed/(real)shock_incr;
      shock_speed=0.0;
    }
  }
#endif

  if (0==imdrestart) {
    /* initialize temperature, if necessary */
    if (do_maxwell) maxwell(temperature);
    do_maxwell=0;
#ifdef TTM
    /* let the electron system know about the new temperatures */
    update_fd();
    ttm_overwrite();
    ttm_fill_ghost_layers();
#endif
  }

#ifdef LASER
  if (0==myid) {
    printf( "Parameter laser_rescale_mode is %d\n", laser_rescale_mode );
    printf( "Parameter laser_delta_temp is  %1.10f\n", laser_delta_temp );
#ifndef TWOD
    printf( "Laser irradiates from direction (%d, %d, %d)\n", laser_dir.x,
	    laser_dir.y, laser_dir.z);
#else
    printf( "Laser irradiates from direction (%d, %d)\n", laser_dir.x,
	    laser_dir.y);
#endif /*TWOD*/

    if (laser_mu==0.0) {
      printf( "Absorption length is infinite.\n" );
    } else {
      printf( "Absorption length is %1.10f\n", 1.0/laser_mu );
    }
    printf( "Laser energy density is %1.10f\n", laser_sigma_e);
    printf( "Laser pulse duration (sigma) is %1.10f\n", laser_sigma_t);
    printf( "Time t_0 of laser pulse is %1.10f (%1.10f time steps after start of simulation)\n\n", laser_t_0, laser_t_0/timestep);
#ifdef TTM
    printf( "Using Two Temperature Model TTM\n");
#endif /*TTM*/
  }
#endif /*LASER*/

#if defined(FRAC) || defined(FTG) 
  if (0==myid) {
      printf( "Strain rate is  %1.10f\n", dotepsilon0 );
      printf( "Damping mode  %d\n", dampingmode );
      printf( "Damping prefactor is %1.10f\n\n", gamma_bar );
  }
#endif

#ifdef FBC
  init_fbc();
#endif

#if defined(CORRELATE) || defined(MSQD)
  init_correl(ncorr_rmax,ncorr_tmax);
#endif

#ifdef NMOLDYN
  if (nmoldyn_int > 0) init_nmoldyn();
#endif

#ifdef ATDIST
  if (atdist_int > 0) init_atdist();
#endif

#ifdef DIFFPAT
  if (diffpat_int > 0) init_diffpat();
#endif

#ifdef CG
  if (ensemble == ENS_CG) reset_cg();
#endif

#ifdef DEFORM
  deform_int = 0; 
#endif

  /* simulation loop */
  for (steps=steps_min; steps <= steps_max; ++steps) {

#ifdef SHOCK
    /* accelerate blocks */
    if (shock_mode == 3 && shock_incr>0) { 
      if (steps<=shock_incr){  
	  shock_speed+=dshock_speed;
	  for (k=0; k<ncells; ++k) {
	      cell *p;
	      p = cell_array + CELLS(k);
	      for (i=0; i<p->n; ++i) {
		  IMPULS(p,i,X) += dshock_speed * MASSE(p,i);
	      }
	  }
      }
    }
#endif

#ifdef STRESS_TENS
    do_press_calc = (((eng_int   > 0) && (0 == steps % eng_int  )) ||
                     ((press_int > 0) && (0 == steps % press_int)) ||
#ifdef FORCE
		     ((force_int  > 0) && (0 == steps % force_int )) ||
#endif /* FORCE */
                     ((dist_int  > 0) && (0 == steps % dist_int )) ||
                     (relax_rate > 0.0) );
#endif

#ifdef EPITAX
    for (i=0; i<ntypes; ++i ) {
      if ( (steps >= epitax_startstep) && (steps <= epitax_maxsteps) ) {  
	if ( (epitax_rate[i] != 0) && ((steps-steps_min)%epitax_rate[i])==0 ) {
	  delete_atoms(); 
	  create_atom(i, epitax_mass[i], epitax_temp[i]);
	}
	else if ( (epitax_rate[i] != 0) && (steps > epitax_maxsteps) && 
                  ( (steps-steps_min)%epitax_rate[i]) == 0 ) 
	  delete_atoms();
      }
    }
#endif

#ifdef HOMDEF
    if ((lindef_int > 0) && (0 == steps % lindef_int)) 
#ifdef TWOD
      lin_deform(lindef_x, lindef_y,           lindef_size);
#else
      lin_deform(lindef_x, lindef_y, lindef_z, lindef_size);
#endif
#endif

#ifdef DEFORM
    if (max_deform_int > 0) {
#ifdef RELAX
      if ((is_relaxed) || (deform_int == max_deform_int))
#else
      if (deform_int == max_deform_int)
#endif
      {
        deform_sample();
        deform_int=0;
#ifdef CG
        if (ensemble == ENS_CG) reset_cg();
#endif
      }
      deform_int++;
    }
#endif

#ifdef CNA
    if (steps <= cna_end) {
      if (0 == (steps - cna_start)%(cna_int)) {
	/* activate computation of neighbour tables */ 
	for (i=0; i<ntypes*ntypes; i++)
	  neightab_r2cut[i] = cna_rcut * cna_rcut;
	/* activate CNA */
	cna = 1;
	cna_pairs = 0;
      }
    }
#endif

#ifdef AVPOS
    if ((steps == steps_min) || (steps == avpos_start)) {
       update_avpos();
    }
#endif
#ifdef ATDIST
    if ((atdist_int > 0) && (steps >= atdist_start) && 
        (steps <= atdist_end)) update_atdist();
#endif
#ifdef DIFFPAT
    if ((diffpat_int > 0) && (steps >= diffpat_start) && 
        (steps <= diffpat_end)) update_diffpat(steps);
#endif

#ifdef FBC
    update_fbc();
#endif

#ifdef TIMING
    imd_start_timer(&time_forces);
#endif
#if defined (CG) && !defined(ACG)
    if (ensemble == ENS_CG) cg_step(steps);
    else
#elif defined(ACG)
    if (ensemble == ENS_CG) acg_step(steps);
    else
#endif
    calc_forces(steps);
#ifdef RIGID
    calc_superforces();
#endif
#ifdef NEB
    calc_forces_neb();
#endif
#ifdef TIMING
    imd_stop_timer(&time_forces);
#endif

#ifdef FORCE
    /* we have to write the forces *before* the atoms are moved */
    if ((force_int > 0) && (0 == steps % force_int)) 
       write_config_select( steps/force_int, "force",
                            write_atoms_force, write_header_force);
#endif

#ifdef WRITEF /* can be used as tool for postprocessing */
    if ((force_int > 0) && (0 == steps % force_int)) 
      write_config_select( steps/force_int, "wf", 
                           write_atoms_wf, write_header_wf);
#endif

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

#ifdef CNA
    if (cna) {
      do_cna();
      if (0==myid && cna_write_statistics) {
	/* works not correctly in parallel version */
	sort_pair_types();
	write_statistics();
      }
      for (k=0; k<MAX_TYPES; k++)
	type_count[k] = 0;
      /* deactivate computation of neighbour tables */
      for (i=0; i<ntypes*ntypes; i++)
	neightab_r2cut[i] = -1.0;
      /* deactivate CNA */
      cna = 0;
      /* write CNA atoms */
      cna_writec = 1;
      for (k=0; k<cna_write_n; k++) {
	sprintf(tmp_str, "%d.cna", cna_writev[k]);
	write_config_select((steps-1)/cna_int, tmp_str,
			    write_atoms_cna, write_header_cna);
	cna_writec *= 2;
      }
    }
#endif

#if defined(CORRELATE) || defined(MSQD)
    if ((steps >= correl_start) && ((steps < correl_end) || (correl_end==0))) {
      int istep = steps - correl_start;
      if (istep % correl_ts == 0) 
        correlate(steps, correl_refstep, istep/correl_ts);
      if ((correl_int != 0) && (steps-correl_refstep+1 >= correl_int)) 
        correl_refstep += correl_int;
    }
#endif

#ifdef GLOK 
    /* "global convergence": set momenta to 0 if P*F < 0 (global vectors) */
    if (ensemble == ENS_GLOK) {
      update_glok();
    }
#endif

#ifdef LASER
    /* do rescaling of atom velocities/electron temperature source terms */
    /* to simulate absorption of laser pulse */
    do_laser_rescale();
#endif
#ifdef TTM
    calc_ttm();
#endif

#ifdef TIMING
    imd_start_timer(&time_integrate);
#endif
#if !defined(CBE) || !defined(SPU_INT)
    if (ensemble != ENS_CG) move_atoms(); /* here PxF is recalculated */
#endif
#ifdef TIMING
    imd_stop_timer(&time_integrate);
#endif

#ifdef EPITAX
    /* beam atoms are always integrated by NVE */
    if (ensemble != ENS_NVE) move_atoms_nve();
#endif

#ifdef TEMPCONTROL 
    increment_temperature();
#endif

    /* Periodic I/O */
#ifdef TIMING
    imd_start_timer(&time_output);
#endif
    if ((checkpt_int > 0) && (0 == steps % checkpt_int)) 
       write_config( steps/checkpt_int, steps);
    if ((eng_int  > 0) && (0 == steps % eng_int )) write_eng_file(steps);
    if ((dist_int > 0) && (0 == steps % dist_int)) write_distrib(steps);
    if ((pic_int  > 0) && (0 == steps % pic_int )) write_pictures(steps);

#ifdef TTM
    if ((ttm_int > 0) && (0 == steps % ttm_int)) ttm_writeout(steps/ttm_int);
#endif
#ifdef EFILTER  /* just print atoms if in an energy-window */ 
    if ((ef_checkpt_int > 0) && (0 == steps % ef_checkpt_int)) 
       write_config_select( steps/ef_checkpt_int, "ef",
                            write_atoms_ef, write_header_ef);
#endif
#ifdef NNBR  /* just print atoms by neighbour condition */ 
    if ((nb_checkpt_int > 0) && (0 == steps % nb_checkpt_int)) 
       write_config_select( steps/nb_checkpt_int, "nb",
                            write_atoms_nb, write_header_nb);
#endif
#ifdef ATDIST
    if ((atdist_pos_int > 0) && (0 == steps % atdist_pos_int))
       write_config_select( steps / atdist_pos_int, "cpt", 
                            write_atoms_atdist_pos, write_header_atdist_pos);
#endif
#ifdef DISLOC
    if (steps == up_ort_ref) update_ort_ref();
    if ((dem_int > 0) && (0 == steps % dem_int)) 
       write_config_select(steps, "dem", write_atoms_dem, write_header_dem);
    if ((dsp_int > 0) && (steps > up_ort_ref) && (0 == steps % dsp_int)) 
       write_config_select(steps, "dsp", write_atoms_dsp, write_header_dsp);
#endif
#ifdef AVPOS
    if ( steps <= avpos_end ){
      if ((avpos_res > 0) && (0 == (steps - avpos_start) % avpos_res) && 
          (steps > avpos_start)) add_positions();
      if ((avpos_int > 0) && (0 == (steps - avpos_start) % avpos_int) && 
          (steps > avpos_start)) {
        write_config_select((steps - avpos_start) / avpos_int,"avp",
                            write_atoms_avp,write_header_avp);
        write_avpos_itr_file((steps - avpos_start) / avpos_int, steps);
        update_avpos();
      }
    }
#endif
#ifdef TRANSPORT 
    if ((tran_int > 0) && (steps > 0) && (0 == steps%tran_int)) 
       write_temp_dist(steps);
#endif
#ifdef RNEMD
    if ((exch_int > 0) && (0 == steps%exch_int)) 
       rnemd_heat_exchange();
#endif

#ifdef STRESS_TENS
    if ((press_int > 0) && (0 == steps % press_int)) {
      if (!do_press_calc) error("pressure tensor incomplete");
       write_config_select( steps/press_int, "press",
                            write_atoms_press, write_header_press);
    }
#endif
#ifdef NMOLDYN
    if ((nmoldyn_int > 0) && (0 == steps % nmoldyn_int)) write_nmoldyn(steps);
#endif
#ifdef DSF
    if ((dsf_int > 0) && (0 == steps % dsf_int)) write_dsf();
#endif

#ifdef SOCKET_IO
    if ((socket_int > 0) && (0 == steps % socket_int)) check_socket();
#endif

#ifdef TIMING
    imd_stop_timer(&time_output);
#endif

#ifdef HOMDEF
    if (relax_rate > 0.0) relax_pressure();
#endif

#ifdef RELAX
    check_relaxed();
#endif

#ifdef NEB
    if ((0==myrank) && (neb_eng_int > 0) && (0 == steps % neb_eng_int ))
      write_neb_eng_file(steps);
#endif

#ifdef NBLIST
    check_nblist();
#else
    fix_cells();  
#endif

#ifdef ATDIST
    if ((atdist_int > 0) && (steps == atdist_end)) write_atdist();
#endif
#ifdef DIFFPAT
    if ((diffpat_int > 0) && (steps == diffpat_end)) write_diffpat();
#endif
#ifdef MSQD
    if ((correl_end >0) && (steps==correl_end) || 
        (correl_end==0) && (steps==steps_max))   
      write_config_select(0, "sqd", write_atoms_sqd, write_header_sqd);
#endif

    /* write checkpoint, if empty write file is found */
    if ((watch_int > 0) && (0==steps%watch_int)) check_write();

    /* finish, if stop file is found */
    if ((stop_int > 0) && (0==steps%stop_int)) {
      if ((finished = check_stop())) break;
    }

    /* finish, if maxwalltime is reached */
    if (maxwalltime > 0) {
      if ((finished = check_walltime())) break;
    }
  }

  /* clean up the current phase, and clear restart flag */
  imdrestart=0;
  if (0==myid) {
    write_itr_file(-1, steps_max,"");
    if (0==myrank) printf( "End of simulation %d\n", simulation );
  }  
  return finished;
}

/*****************************************************************************
*
*  close .eng file and Co.
*
*****************************************************************************/

void close_files(void)
{
  if (0 == myid) {
    if (NULL!= eng_file) { fclose( eng_file); eng_file  = NULL; }
    if (NULL!=msqd_file) { fclose(msqd_file); msqd_file = NULL; }
  }
}

#ifdef GLOK

/*****************************************************************************
*
*  update state of (adaptive) glok integrator
*
*****************************************************************************/

void update_glok(void)
{
  int i, k;

  if (steps == steps_min) {
    glok_start = steps_min; 
#ifdef MIX
    mix = glok_mix;
#endif
  }
  glok_int = steps - glok_start;

  /* always start glok with new dynamics, not with old velocities */
  if (glok_int == 0) {
    for (k=0; k<NCELLS; ++k) {
      cell *p = CELLPTR(k);
      for (i=0; i<p->n; ++i) {
        IMPULS(p,i,X) = 0.0;
        IMPULS(p,i,Y) = 0.0;
#ifndef TWOD
        IMPULS(p,i,Z) = 0.0;
#endif
      }
    }
  }

#ifdef ADAPTGLOK
  /* increase the timestep, but not immediately after P*F was < 0 */ 
  if ( (nPxF>= min_nPxF)  && (glok_int > glok_minsteps)) {
    timestep = (timestep * glok_incfac > glok_maxtimestep) ? 
               glok_maxtimestep : timestep * glok_incfac;
#ifdef MIX
    mix *= glok_mixdec;
#endif
  }
#endif

#ifdef MIX
  if (fnorm >=1e-20) mixforcescalefac = SQRT(pnorm/fnorm);
#endif
  real ekin = 2 * tot_kin_energy / nactive;

  if ((PxF < 0.0) || (ekin > glok_ekin_threshold)  || 
      (sqrt(f_max2) >= glok_fmaxcrit) && steps >5) {
#ifdef ADAPTGLOK
    if (PxF < 0.0) nPxF++;
    /* decrease timestep, but only when it has been increased before */
    if (glok_int > glok_minsteps ) {
      if (timestep > glok_maxtimestep/50.0) timestep *=glok_decfac;
    }
#endif

#ifdef MIX
    mix = glok_mix;
#endif
    for (k=0; k<NCELLS; ++k) {
      cell *p = CELLPTR(k);
      for (i=0; i<p->n; ++i) {
        ORT(p,i,X) -= 0.5* timestep / MASSE(p,i) * IMPULS(p,i,X);
        ORT(p,i,Y) -= 0.5* timestep / MASSE(p,i) * IMPULS(p,i,Y);
#ifndef TWOD
        ORT(p,i,Z) -= 0.5* timestep / MASSE(p,i) * IMPULS(p,i,Z);
#endif
        IMPULS(p,i,X) = 0.0;
        IMPULS(p,i,Y) = 0.0;
#ifndef TWOD
        IMPULS(p,i,Z) = 0.0;
#endif
      }
    }
    glok_start = steps;
  }
}

#endif

#ifdef TEMPCONTROL

/*****************************************************************************
*
*  increment temperatue
*
*****************************************************************************/

void increment_temperature(void)
{
  static real dtemp = 0.0;

  if (steps==steps_min) {
    if (use_curr_temp==1) {
#ifdef UNIAX
      temperature = 2.0 * tot_kin_energy / (nactive + nactive_rot);
#else
      temperature = 2.0 * tot_kin_energy / nactive;
#endif
      use_curr_temp = 0;
    }
    dtemp = (end_temp - temperature) / (steps_max - steps_min);
  }
  temperature += dtemp;
}

#endif

#ifdef FBC

/*****************************************************************************
*
*  initialize FBC
*
*****************************************************************************/

void init_fbc(void)
{
  int l, steps_diff = steps_max - steps_min;

  have_fbc_incr = 0;
  do_fbc_incr = 0;

#ifdef RELAX
  fbc_int = 0;
  for (l=0; l<vtypes; l++) {
    if ((fbc_dforces+l)->x != 0.0) have_fbc_incr = 1;
    if ((fbc_dforces+l)->y != 0.0) have_fbc_incr = 1;
#ifndef TWOD
    if ((fbc_dforces+l)->z != 0.0) have_fbc_incr = 1;
#endif
  }
#else
  /* dynamic loading, increment linearly at each timestep */
  if ((ensemble!=ENS_MIK) && (ensemble!=ENS_GLOK) && (ensemble!=ENS_CG)) {
    for (l=0;l<vtypes;l++){
      (fbc_df+l)->x = ((fbc_endforces+l)->x-(fbc_beginforces+l)->x)/steps_diff;
      (fbc_df+l)->y = ((fbc_endforces+l)->y-(fbc_beginforces+l)->y)/steps_diff;
#ifndef TWOD
      (fbc_df+l)->z = ((fbc_endforces+l)->z-(fbc_beginforces+l)->z)/steps_diff;
#endif
    }
    do_fbc_incr = 1;
  }
#endif /* RELAX */

}

/*****************************************************************************
*
*  update FBC
*
*****************************************************************************/

void update_fbc()
{
  int l;
#ifdef TWOD
  vektor nullv={0.0,0.0};
#else
  vektor nullv={0.0,0.0,0.0};
#endif

#ifdef RELAX
  /* set fbc increment if necessary */
  if ((ensemble==ENS_MIK) || (ensemble==ENS_GLOK) || (ensemble==ENS_CG)) {
    if ((is_relaxed) || (fbc_int > max_fbc_int)) { 
      for (l=0; l<vtypes; l++) fbc_df[l] = fbc_dforces[l];
      fbc_int = 0;
      do_fbc_incr = 1;
    }
    else {
      for (l=0; l<vtypes; l++) *(fbc_df+l) = nullv;
      do_fbc_incr = 0;
      fbc_int++;
    }
  }
#endif

  /* apply fbc increment if necessary */
  if (do_fbc_incr == 1) {
    for (l=0; l<vtypes; l++){ 
      (fbc_forces+l)->x += (fbc_df+l)->x;  
      (fbc_forces+l)->y += (fbc_df+l)->y;
#ifndef TWOD
      (fbc_forces+l)->z += (fbc_df+l)->z;
#endif
    } 
#ifdef CG
    if (ensemble == ENS_CG) reset_cg();
#endif
  }
}

#endif /* FBC */

/*****************************************************************************
*
*  check if sample is relaxed
*
*****************************************************************************/

#ifdef RELAX

void check_relaxed(void)
{
  is_relaxed = 0;

  if ((ensemble==ENS_MIK) || (ensemble==ENS_GLOK) || (ensemble==ENS_CG)) {

    int stop = 0;
    real fnorm2, ekin, epot, delta_epot;
#ifdef NEB
    MPI_Allreduce( &fnorm, &neb_fnorm, 1, REAL, MPI_SUM, MPI_COMM_WORLD);
    neb_fnorm = SQRT( neb_fnorm / (nactive * (neb_nrep-2)) );
    if (neb_fnorm < fnorm_threshold) is_relaxed = 1;
    else is_relaxed = 0;
#else
    fnorm2 = SQRT( fnorm / nactive );
    ekin   = 2 * tot_kin_energy / nactive;
    epot   = tot_pot_energy / natoms;
    delta_epot = old_epot - epot;
    if (delta_epot < 0) delta_epot = -delta_epot;

    if ((ekin  <  ekin_threshold) || (fnorm2 < fnorm_threshold) || 
        (delta_epot < delta_epot_threshold)) is_relaxed = 1;
    else is_relaxed = 0;

    old_epot = epot;
#endif

    if (is_relaxed) {
      stop = 1;
      write_eng_file(steps);
      write_ssconfig(steps);
#ifdef NEB
      if (0==myrank) write_neb_eng_file(steps);
#else
      if (0==myid) {
        printf("nfc = %d epot = %22.16f\n", nfc, epot );
        printf("ekin = %e fnorm = %e f_max = %e delta_epot = %e\n", 
               ekin, fnorm2, f_max, delta_epot);
      }
#endif
    }

#ifdef DEFORM
    if (max_deform_int > 0) stop=0;
#endif
#ifdef FBC
    if (have_fbc_incr) stop=0;
#endif
    if (stop) steps_max = steps;
  }

}

#endif /* RELAX */


#ifdef RIGID

/*****************************************************************************
*
*  calculate force on superparticles
*
*****************************************************************************/

void calc_superforces(void)
{
  int i, k;
  real tmpvec1[DIM], tmpvec2[DIM];

  /* total force on superparticles (for each cpu) */
  for(k=0; k<ncells; k++) {
    cell *p;
    int sorte;
    p = CELLPTR(k);
    for (i=0; i<p->n; i++) {
      sorte = VSORTE(p,i);
      if ( superatom[sorte] > -1 ) {
        superforce[superatom[sorte]].x += KRAFT(p,i,X);
	superforce[superatom[sorte]].y += KRAFT(p,i,Y);
#ifndef TWOD
	superforce[superatom[sorte]].z += KRAFT(p,i,Z);
#endif

#ifdef FBC
	superforce[superatom[sorte]].x += (fbc_forces+sorte)->x;
	superforce[superatom[sorte]].y += (fbc_forces+sorte)->y;
#ifndef TWOD
	superforce[superatom[sorte]].z += (fbc_forces+sorte)->z;
#endif
#endif
      }
    }
  }

#ifdef MPI
  /* total force on superparticles */
  for (i=0; i<nsuperatoms; i++) {
    tmpvec1[0] = superforce[i].x;
    tmpvec1[1] = superforce[i].y;
#ifndef TWOD
    tmpvec1[2] = superforce[i].z;
#endif
    MPI_Allreduce( tmpvec1, tmpvec2, DIM, REAL, MPI_SUM, cpugrid); 
    superforce[i].x = tmpvec2[0];
    superforce[i].y = tmpvec2[1];
#ifndef TWOD
    superforce[i].z = tmpvec2[2];
#endif
  }
#endif

}

#endif /* RIGID */

/*****************************************************************************
*
*  check write file, and stop if present
*
*****************************************************************************/

void check_write(void)
{
  int write = 0;
  if (0 == myid) {
    FILE *testfile = fopen("write","r");
    if (NULL!=testfile && fgetc(testfile)==EOF) { 
      fclose(testfile);
      remove("write");
      write = 1;
    }
  }
#ifdef MPI
  MPI_Bcast( &write, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  if (write) {
    if (myid == 0) printf("Write file found after %d steps\n", steps);
    write_config(-2,steps);
  }
}

/*****************************************************************************
*
*  check stop file, and stop if present
*
*****************************************************************************/

int check_stop(void)
{
  int stop = 0;
  if (0 == myid) {
    FILE *testfile = fopen("stop","r");
    if (NULL!=testfile) { 
      fclose(testfile);
      unlink("stop");
      stop = 1;
    }
  }
#ifdef MPI
  MPI_Bcast( &stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  if (stop) {
    if (myid == 0) printf("Stop file found after %d steps\n", steps);
    write_config(-1,steps);
    steps_max = steps;
  }
  return stop;
}

/*****************************************************************************
*
*  check walltime, and stop if necessary
*
*****************************************************************************/

int check_walltime(void)
{
  int stop = 0;
  double tdiff = difftime(time(&tend), tstart);
#ifdef MPI
  MPI_Bcast( &tdiff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
  if (tdiff > maxwalltime) {
    if (myid == 0) 
      printf("Maximal allowed walltime reached after %d steps\n", steps);
    write_config(-1,steps);
    steps_max = steps;
    stop = 1;
  }
  return stop;
}

/******************************************************************************
*
* do_boundaries
*
* Apply periodic boundaries to all atoms
* Could change so that only cells on surface to some work
*
******************************************************************************/

void do_boundaries(void)
{
  int k;

  /* for each cell in bulk */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (k=0; k<NCELLS; ++k) {

    int  l;
    real i;  /* FLOOR returns a real */
    cell *p;

    p = CELLPTR(k);

    /* PBC in x direction */
    if (pbc_dirs.x==1)
    for (l=0; l<p->n; ++l) {
      i = -FLOOR( SPRODX(ORT,p,l,tbox_x) );
      ORT(p,l,X)     += i * box_x.x;
      ORT(p,l,Y)     += i * box_x.y;
#ifndef TWOD
      ORT(p,l,Z)     += i * box_x.z;
#endif
#if defined(MSQD) || defined(NMOLDYN)
      REF_POS(p,l,X) += i * box_x.x;
      REF_POS(p,l,Y) += i * box_x.y;
#ifndef TWOD
      REF_POS(p,l,Z) += i * box_x.z;
#endif
#endif
#ifdef AVPOS
      SHEET(p,l,X)   -= i * box_x.x;
      SHEET(p,l,Y)   -= i * box_x.y;
#ifndef TWOD
      SHEET(p,l,Z)   -= i * box_x.z;
#endif
#endif
    }

    /* PBC in y direction */
    if (pbc_dirs.y==1)
    for (l=0; l<p->n; ++l) {
      i = -FLOOR( SPRODX(ORT,p,l,tbox_y) );
      ORT(p,l,X)     += i * box_y.x;
      ORT(p,l,Y)     += i * box_y.y;
#ifndef TWOD
      ORT(p,l,Z)     += i * box_y.z;
#endif
#if defined(MSQD) || defined(NMOLDYN)
      REF_POS(p,l,X) += i * box_y.x;
      REF_POS(p,l,Y) += i * box_y.y;
#ifndef TWOD
      REF_POS(p,l,Z) += i * box_y.z;
#endif
#endif
#ifdef AVPOS
      SHEET(p,l,X)   -= i * box_y.x;
      SHEET(p,l,Y)   -= i * box_y.y;
#ifndef TWOD
      SHEET(p,l,Z)   -= i * box_y.z;
#endif
#endif
    }

#ifndef TWOD
    /* PBC in z direction */
    if (pbc_dirs.z==1)
    for (l=0; l<p->n; ++l) {
      i = -FLOOR( SPRODX(ORT,p,l,tbox_z) );
      ORT(p,l,X)     += i * box_z.x;
      ORT(p,l,Y)     += i * box_z.y;
      ORT(p,l,Z)     += i * box_z.z;
#if defined(MSQD) || defined(NMOLDYN)
      REF_POS(p,l,X) += i * box_z.x;
      REF_POS(p,l,Y) += i * box_z.y;
      REF_POS(p,l,Z) += i * box_z.z;
#endif
#ifdef AVPOS
      SHEET(p,l,X)   -= i * box_z.x;
      SHEET(p,l,Y)   -= i * box_z.y;
      SHEET(p,l,Z)   -= i * box_z.z;
#endif
    }
#endif
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

  real tmp_presstens1[6], tmp_presstens2[6];

  if (!do_press_calc) error("pressure tensor incomplete");
  tot_presstens.xx = 0.0; 
  tot_presstens.yy = 0.0; 
  tot_presstens.xy = 0.0;
#ifndef TWOD
  tot_presstens.zz = 0.0; 
  tot_presstens.yz = 0.0;
  tot_presstens.zx = 0.0;
#endif

  /* sum up total pressure tensor */
  for (i=0; i<NCELLS; ++i) {
    int j;
    cell *p;
    p = CELLPTR(i);
    for (j=0; j<p->n; ++j) {
      tot_presstens.xx += PRESSTENS(p,j,xx);
      tot_presstens.yy += PRESSTENS(p,j,yy);
      tot_presstens.xy += PRESSTENS(p,j,xy);  
#ifndef TWOD
      tot_presstens.zz += PRESSTENS(p,j,zz);
      tot_presstens.yz += PRESSTENS(p,j,yz);  
      tot_presstens.zx += PRESSTENS(p,j,zx);  
#endif
    }
  }

#ifdef MPI

  tmp_presstens1[0] = tot_presstens.xx; 
  tmp_presstens1[1] = tot_presstens.yy; 
  tmp_presstens1[2] = tot_presstens.xy;
#ifndef TWOD
  tmp_presstens1[3] = tot_presstens.zz; 
  tmp_presstens1[4] = tot_presstens.yz;
  tmp_presstens1[5] = tot_presstens.zx;
#endif

#ifdef TWOD
  MPI_Allreduce( tmp_presstens1, tmp_presstens2, 3, REAL, MPI_SUM, cpugrid);
#else
  MPI_Allreduce( tmp_presstens1, tmp_presstens2, 6, REAL, MPI_SUM, cpugrid);
#endif

  tot_presstens.xx  = tmp_presstens2[0];
  tot_presstens.yy  = tmp_presstens2[1]; 
  tot_presstens.xy  = tmp_presstens2[2]; 
#ifndef TWOD
  tot_presstens.zz  = tmp_presstens2[3];
  tot_presstens.yz  = tmp_presstens2[4]; 
  tot_presstens.zx  = tmp_presstens2[5]; 
#endif

#endif /* MPI */

}

#endif/* STRESS_TENS */
