
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

/* if bond boost method is used, this main loop is the main loop for the
   MD part. The relaxation is done in a copy of this main loop in imd_bboost.c
   function bb_minimize. if you change something here, please also change in 
   bb_minimize */

int main_loop(int simulation)
{
  int  finished = 0;
  int  i, j, k, l;
  int  steps_diff = steps_max - steps_min;
  int  deform_int = 0, do_fbc_incr = 0;
  int  have_fbc_incr = 0;
  real dtemp, dshock_speed;
  real ri;
  vektor d_pressure;
  real tmpvec1[DIM], tmpvec2[DIM];
  char tmp_str[9];

  int tmpsteps;

  real fnorm2,ekin,epot,delta_epot;
#ifdef ADA
  int adaDone;
#endif
#ifdef NYETENSOR
int nyeDone;
#endif
#ifdef debugLo
    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("passing main loop from checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("    ************************* \n");fflush(stdout);
#endif
  
#ifdef GLOK
  if(glok_start <=steps_min)
      glok_start=steps_min; 
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

#if defined(FRAC) || defined(FTG) 
  if (0==myid) {
      printf( "Strain rate is  %1.10f\n", dotepsilon0 );
      printf( "Damping mode  %d\n", dampingmode );
      printf( "Damping prefactor is %1.10f\n\n", gamma_bar );
  }
#endif

#ifndef BBOOST /* in bond boost, these initializations are done in init_bboost */
#ifdef EXTPOT
  init_extpot();
#endif  

#ifdef USEFCS
  init_fcs();
#endif
  
#ifdef FBC
  init_fbc();
#ifdef BEND
  init_bfbc();
#endif
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
#if defined(HOMDEF) && defined(RELAX)
    deform_int = 0; 
#endif

#endif /* not BBOOST */ 

  imd_start_timer(&time_main);
  
  /* simulation loop */
  for (steps=steps_min; steps <= steps_max; ++steps) {

#ifdef NNBR_TABLE
  nnbr_done = 0;
#endif
#ifdef ADA
  adaDone = 0;
#endif
#ifdef NYETENSOR
  nyeDone = 0;
#endif

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
#ifdef HC
                     ((hc_int > 0) && (steps >= hc_av_start) && 
                      ((steps < hc_start) || ((steps-hc_start)%hc_int==0))) ||
#endif
                     (relax_rate > 0.0) );
#endif /* STRESS_TENS */

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

#ifdef EXTPOT
    /* update extpot position if necessary */
#ifdef RELAX
    if ( ((ensemble==ENS_MIK) || (ensemble==ENS_GLOK) || (ensemble==ENS_CG)) &&
         (ep_max_int > 0) ) {
        if ((is_relaxed) || (ep_int > ep_max_int)) {
            write_ssdef(steps);    /* write info for quasistat simulations */
            write_fext(steps);     /* update .ind file */
            write_ssconfig(steps); /* write config, even when not fully relaxed */
            move_extpot(1.0);
            is_relaxed=0;
            ep_int = 0;
            is_relaxed = 0;
#ifdef GLOK
            if (ensemble==ENS_GLOK)
            {
                reset_glok();
            }
#endif
        }
        ep_int++;
    }
    else 
#endif
    move_extpot(timestep);
#endif /* EXTPOT */

#if defined(HOMDEF) && defined(CYCLE)
    if ((lindef_int > 0) && (0 == steps % lindef_int))
    {
        ri =( (lindef_size -1.0) * sin(2.0*3.141592653589793238*lindef_freq*timestep*(steps+1))+1.0 ) /
            ( (lindef_size -1.0) * sin(2.0*3.141592653589793238*lindef_freq*timestep*steps)+1.0 ) - 1.0 ;
        //     printf(" ri = %f box_x= %f\n",ri,box_x.x*(ri+1.0)); fflush(stdout);
#ifdef TWOD
    lin_deform(lindef_x, lindef_y,       ri);
#else
    lin_deform(lindef_x, lindef_y, lindef_z, ri);
#endif
    }
#endif

#if defined(HOMDEF) && defined(RELAX)
    if ( (ensemble==ENS_MIK) || (ensemble==ENS_GLOK) || (ensemble==ENS_CG) )
    {
        if(lindef_int >0)
        {
            //lindef_int now plays the role of  max_deform_int
            // have to find something which plays the fole of deformint
#ifdef GLOK
            if ( (is_relaxed && (steps-glok_start>=10)) ||  (0 == deform_int % lindef_int))
#else  /* GLOK */
            if ( (is_relaxed ) ||  (0 == deform_int % lindef_int))
#endif
            {
                write_ssdef(steps);
                write_ssconfig(steps); /* write config, even when not fully relaxed */
#ifdef TWOD
                lin_deform(lindef_x, lindef_y,           lindef_size);
#else
                lin_deform(lindef_x, lindef_y, lindef_z, lindef_size);
#endif
                deform_int=0;
                is_relaxed=0;

                make_box();
#ifdef NBLIST
                have_valid_nbl = 0;
#endif
                fix_cells();  

#ifdef GLOK
               if (ensemble==ENS_GLOK)
               {
                   reset_glok();
               }
#endif
           }
            deform_int++;
        }
    }
#endif
    
#if defined(HOMDEF) && !defined(RELAX)
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
      {
          write_ssdef(steps);
          write_ssconfig(steps); /* write config, even when not fully relaxed */
#ifdef GLOK
          if (ensemble==ENS_GLOK)
          {
              reset_glok();
          }
#endif
#else
      if (deform_int == max_deform_int)
      {
#endif
      
        deform_sample();
#ifdef RELAX
        is_relaxed=0;
#endif
        deform_int=0;

#ifdef NBLIST
    check_nblist();
#else
    fix_cells();  
#endif

#ifdef CG
        if (ensemble == ENS_CG) reset_cg();
#endif
      }
      deform_int++;
    }
#endif

#ifdef CNA
#ifdef debugLo
    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("passing "" input cna "" checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("    ************************* \n");fflush(stdout);
#endif 
    if (steps <= cna_end) {
      if (0 == (steps - cna_start)%(cna_int)
	  || ((cna_crist>0) && (checkpt_int > 0) 
	      && (0 == steps % checkpt_int))) {
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
#ifdef BEND
    update_bfbc();
#endif
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
#ifdef USEFCS
#ifdef PAIR
      calc_forces(steps);
#endif
      calc_forces_fcs(steps);
#else
      calc_forces(steps);
#endif
/* #ifdef BBOOST
    calc_bondboost(steps);
#endif */
#ifdef EXTPOT
    calc_extpot();
#endif
#ifdef RIGID
    calc_superforces();
#endif    
#ifdef NEB
    calc_forces_neb();
#endif
#ifdef FEFL
    calc_fefl();
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
#ifdef debugLo
    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("passing "" go into do_cna "" checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("    ************************* \n");fflush(stdout);
#endif 
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
      if( cna_crist > 0 )
	write_config_select(steps/cna_int, "crist",
			    write_atoms_crist, write_header_crist);
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


#ifdef GETSADDLE
  if(last_PxF<0.0 && PxF>=0)
  {
      write_saddleconfig(steps);
  }
  last_PxF=PxF;
#endif

#ifdef GETMIN
  if(last_PxF>=0.0 && PxF<0)
  {
      write_minconfig(steps);
  }
  last_PxF=PxF;
#endif



#ifdef LASER
    /* do rescaling of atom velocities/electron temperature source terms */
    /* to simulate absorption of laser pulse */
    do_laser_rescale();
#endif
#ifdef TTM
    calc_ttm();
#endif

#ifdef SM
#ifdef NBLIST
    if ((!sm_fixed_charges) && ((charge_update_steps > 0) && steps % charge_update_steps == 0)){
      charge_update_sm();
	}
#else
    if ((!sm_fixed_charges) && ((charge_update_steps > 0) && steps % charge_update_steps == 0)){
      do_charge_update();
       }
#endif
#endif

#ifdef HC
    if ((hc_int > 0) && (steps >= hc_av_start) && 
        ((steps < hc_start) || ((steps - hc_start) % hc_int == 0)))
      do_heat_cond(steps);  
#endif

#ifdef TIMING
    imd_start_timer(&time_integrate);
#endif
#if !defined(CBE) || !defined(SPU_INT)
#ifdef NEB
    if(myrank != 0 && myrank != neb_nrep-1)
    {
#endif
    /* move atoms */
    if (ensemble != ENS_CG) move_atoms(); /* here PxF is recalculated */
#ifdef NEB
    }
#endif
    
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

#ifdef ADA
		if ( ((checkpt_int > 0) && (0 == steps % checkpt_int)) ||
				((ada_write_int > 0) && (0 == steps % ada_write_int)) ) {
			if (adaDone==0){
				do_ada();
				adaDone = 1;
			}
#ifdef NYETENSOR
			buildHopsToDefect();
			calculateNyeTensorData();
			nyeDone = 1;
#endif
		}

		if ((ada_write_int > 0) && (0 == steps % ada_write_int)) {
			write_config_select(steps/ada_write_int, "ada",
					write_atoms_ada, write_header_ada);
		}
#endif

#ifdef NEB
    if(myrank != 0 && myrank != neb_nrep-1)
    {
        constrain_move();
    }
#endif

#ifdef ZAPP
    zapp();
#endif
    
#ifdef BEND
    update_bend();
#endif 

/* #ifdef BBOOST
    postpro_boost(steps);
#endif  */

    /* Periodic I/O */
#ifdef TIMING
    imd_start_timer(&time_output);
#endif

#ifdef AVPOS
    if ( steps <= avpos_end && steps > avpos_start ){
      if( avpos_steps == 0)    /* default, for backwards compatibility */
	{ 
	  if ((avpos_res > 0) && (0 == (steps - avpos_start) % avpos_res) )
	    add_positions();
	  if ((avpos_int > 0) && (0 == (steps - avpos_start) % avpos_int) ) 
	    {
	      write_config_select((steps - avpos_start) / avpos_int,"avp",
				  write_atoms_avp,write_header_avp);
	      write_avpos_itr_file((steps - avpos_start) / avpos_int, steps);
	      update_avpos();
	    }
	}
      else
	{
	  if ( steps >= avpos_start + (avpos_nwrites+1)*avpos_int - avpos_res*avpos_steps)
	    {
	      tmpsteps = avpos_start + (avpos_nwrites+1)*avpos_int - steps;
	      //  printf("steps = %d tmpsteps =%d\n",steps,tmpsteps);fflush(stdout); 
	      if ( tmpsteps == avpos_res*avpos_steps)
		{ 
		  //  	  printf("avpos_start = %d avpos_nwrites=%d avpos_int=%d avpos_steps=%d \n",avpos_start,avpos_nwrites, avpos_int,avpos_steps );fflush(stdout); 
		  //	  printf("updating positions\n");fflush(stdout); 
		  update_avpos();
		}
	      else if ( tmpsteps % (avpos_res) == 0)
		{
		  add_positions(); 
		  //	  printf("adding positions\n");fflush(stdout); 
		}
	      if (tmpsteps == 0)
		{
		  avpos_nwrites++;
		  write_config_select(avpos_nwrites,"avp",write_atoms_avp,write_header_avp);
		  write_avpos_itr_file(avpos_nwrites, steps); /* I don't think that this is needed? */
		  //	  printf("writing out: %d\n",avpos_nwrites);fflush(stdout); 
		}

	    }
#ifdef STRESS_TENS
	  if ( (press_int > 0) && (steps >= (avpos_npwrites+1)*press_int - avpos_res*avpos_steps))
	    {
	      tmpsteps = avpos_start + (avpos_npwrites+1)*press_int - steps;
	      if ( tmpsteps == avpos_res*avpos_steps)
		{ 		  
		  update_avpress();
		}
	      else if ( tmpsteps % (avpos_res) == 0)
		{
		  add_presstensors();
		}
	      if (tmpsteps == 0)
		{
		  avpos_npwrites++;
		  write_config_select(avpos_npwrites , "press",
				       write_atoms_press, write_header_press);
		  //	  printf("writing out: %d\n",avpos_npwrites);fflush(stdout); 
		}

	    }
#endif /* STRESS_TENS */

	}
    
    }
#endif /* AVPOS, has to be before write_config so that the correct avpos_nwrites is in the itr file */

    if ((checkpt_int > 0) && (0 == steps % checkpt_int)) 
       write_config( steps/checkpt_int, steps);
    if ((eng_int  > 0) && (0 == steps % eng_int )) write_eng_file(steps);
    if ((dist_int > 0) && (0 == steps % dist_int)) write_distrib(steps);
    if ((pic_int  > 0) && (0 == steps % pic_int )) write_pictures(steps);
#ifdef EXTPOT
#ifdef RELAX
    if ( ((ensemble==ENS_MIK) || (ensemble==ENS_GLOK) || (ensemble==ENS_CG)) &&
         (ep_max_int <= 0) )
#endif
    if ((eng_int > 0) && (0 == steps % eng_int )) write_fext(steps);
#endif
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

#ifdef NVX 
    if ((ensemble == ENS_NVX) && (hc_int > 0) && (steps > hc_start)) 
       write_temp_dist(steps - hc_start);
#endif

#ifdef STRESS_TENS
    if ((press_int > 0) && (0 == steps % press_int)) {
      if (!do_press_calc) error("pressure tensor incomplete for writing .press file");
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

#ifdef NYETENSOR
    if (nyeDone == 1)
    	removeNyeTensorData();
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

      /* finish, if max deformation steps in quasistatic simulation are done */
#ifdef RELAX
#if defined (DEFORM) || defined (HOMDEF) || defined (EXTPOT) || defined (FBC)
    if ( (ensemble==ENS_MIK) || (ensemble==ENS_GLOK) || (ensemble==ENS_CG) ) {
        if ( (max_sscount>0) && (sscount>max_sscount) ) {
                 finished = 1 ;
                 steps_max = steps;
                 break;
            }
    }
#endif
#endif
    
    /* finish, if maxwalltime is reached */
    if (maxwalltime > 0) {
      if ((finished = check_walltime())) break;
    }
  }

  imd_stop_timer(&time_main);
#ifdef USEFCS
  fcs_cleanup();
#endif

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

#ifdef NEB
/*****************************************************************************
*
*  globally constrain the max step size
*
*****************************************************************************/

void constrain_move(void)
{
  int i, k;
  real normp;
  real newxnorm=0;
  real tmp_x_max2=0.0;
  
      if(neb_maxmove !=0.0)
      {
          if (SQRT(xnorm) >neb_maxmove)
          {
              normp=sqrt(pnorm);
              printf("step %d myrank:%d xnorm = %lf maxmove = %lf  normp =%lf x_max =%lf \n",steps,myrank,SQRT(xnorm),neb_maxmove,normp,SQRT(x_max2));
	 
              for (k=0; k<NCELLS; ++k) {
                  cell *p = CELLPTR(k);
                  for (i=0; i<p->n; ++i) {
                      ORT(p,i,X) -= timestep / MASSE(p,i) * IMPULS(p,i,X);
                      ORT(p,i,Y) -= timestep / MASSE(p,i) * IMPULS(p,i,Y);
#ifndef TWOD
                      ORT(p,i,Z) -= timestep / MASSE(p,i) * IMPULS(p,i,Z);
#endif
                      if(normp>0.0)
                      {
                          ORT(p,i,X) += neb_maxmove * IMPULS(p,i,X)/normp;
                          ORT(p,i,Y) += neb_maxmove * IMPULS(p,i,Y)/normp;
#ifndef TWOD
                          ORT(p,i,Z) += neb_maxmove * IMPULS(p,i,Z)/normp;
#endif
                          tmp_x_max2 =  MAX(SQR(neb_maxmove*IMPULS(p,i,X)/normp),tmp_x_max2);
                          tmp_x_max2 =  MAX(SQR(neb_maxmove*IMPULS(p,i,Y)/normp),tmp_x_max2);
                          tmp_x_max2 =  MAX(SQR(neb_maxmove*IMPULS(p,i,Z)/normp),tmp_x_max2);
                          
                          newxnorm += neb_maxmove * neb_maxmove / normp /normp * SPRODN(IMPULS,p,i,IMPULS,p,i);
                      }
                      

                      //IMPULS(p,i,X) = 0.0;
                      // IMPULS(p,i,Y) = 0.0;
#ifndef TWOD
                      // IMPULS(p,i,Z) = 0.0;
#endif
                  }
              }
              printf("myrank:%d newxnorm = %lf new xmax %lf\n",myrank,SQRT(newxnorm),SQRT(tmp_x_max2));
              fflush(stdout);
              xnorm=newxnorm;
              x_max2 = tmp_x_max2;
              
          }
      }
  }
  

#endif


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
      (sqrt(f_max2) >= glok_fmaxcrit))
  {
#ifdef ADAPTGLOK
      if(steps >glok_minsteps) {
          if (PxF < 0.0)
          {
              nPxF++;
              /* decrease timestep, but only when it has been increased before */
        //    if (glok_int > glok_minsteps ) {
              if (timestep > glok_maxtimestep/50.0) timestep *=glok_decfac;
              // }
          }
      }
#endif

#ifdef ADAPTGLOK
     if (glok_int > glok_minsteps)
        mix = glok_mix;
    else
        //   mix=0.5;
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



/*****************************************************************************
*
*  reset state of (adaptive) glok integrator, e.g. after deformation step
*
*****************************************************************************/

void reset_glok(void)
{

    fnorm=9.99e99;
    glok_start = steps;
    // reset velocities
    int i,k;
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
    
    
#ifdef ADAPTGLOK
    timestep = starttimestep;
    nPxF = 0;
    glok_int = 0;
#endif

#ifdef MIX
    mix = 0.0;
#endif
    maxwell(temperature);
    // move_atoms();
#ifdef MIX
    mix = glok_mix; 
#endif  
    
    if (0 == myid)
    {
        printf("Reseting GLOK: T= %f, glok_start=%d, timestep =%f, mix=%f \n\n",temperature, glok_start,timestep,mix);fflush(stdout);
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
    int l, steps_diff;

#ifdef TWOD
  vektor nullv={0.0,0.0};
#else
  vektor nullv={0.0,0.0,0.0};
#endif
    
#if !defined(RELAX)
    steps_diff = steps_max - steps_min;
#endif
  have_fbc_incr = 0;
  do_fbc_incr = 0;
#ifdef RELAX
  fbc_int = 0;
  for (l=0; l<vtypes; l++) {
      /* initialize fbc_df with 0*/
      /* will be updated in fbc_update*/
      *(fbc_df+l) = nullv;
    if ((fbc_dforces+l)->x != 0.0) have_fbc_incr = 1;
    if ((fbc_dforces+l)->y != 0.0) have_fbc_incr = 1;
#ifndef TWOD
    if ((fbc_dforces+l)->z != 0.0) have_fbc_incr = 1;
#endif   
  }
#else
  /* dynamic loading, increment linearly at each timestep */
  if (0 == myid) printf("FBC: vtype  fbc_df.x fbc_df.y fbc_df.z\n");
  if ((ensemble!=ENS_MIK) && (ensemble!=ENS_GLOK) && (ensemble!=ENS_CG)) {
    for (l=0;l<vtypes;l++){
      (fbc_df+l)->x = ((fbc_endforces+l)->x-(fbc_beginforces+l)->x)/steps_diff;
      (fbc_df+l)->y = ((fbc_endforces+l)->y-(fbc_beginforces+l)->y)/steps_diff;
#ifndef TWOD
      (fbc_df+l)->z = ((fbc_endforces+l)->z-(fbc_beginforces+l)->z)/steps_diff;
      if (0 == myid) printf("     %d   %e   %e   %e\n",l,(fbc_df+l)->x,(fbc_df+l)->y,(fbc_df+l)->z);
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
        write_ssdef(steps);
        write_ssconfig(steps); /* write config, even when not fully relaxed */
        for (l=0; l<vtypes; l++) fbc_df[l] = fbc_dforces[l];
        fbc_int = 0;
        do_fbc_incr = 1;
#ifdef CG
    if (ensemble == ENS_CG) reset_cg();
#endif
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
#ifdef GLOK
    if (ensemble==ENS_GLOK)
    {
        reset_glok();
    }
#endif
  }
}

#endif /* FBC */


#ifdef BEND

/*****************************************************************************
*
*  initialize BFBC
*
*****************************************************************************/

void init_bfbc(void)
{
    int l, steps_diff;

#ifdef TWOD
  vektor nullv={0.0,0.0};
#else
  vektor nullv={0.0,0.0,0.0};
#endif
    
#if !defined(RELAX)
    steps_diff = steps_max - steps_min;
#endif
  have_bfbc_incr = 0;
  do_bfbc_incr = 0;
#ifdef RELAX
  bfbc_int = 0;
  for (l=0; l<vtypes; l++) {
      /* initialize fbc_df with 0*/
      /* will be updated in fbc_update*/
      *(fbc_bdf+l) = nullv;
    if ((fbc_bdforces+l)->x != 0.0) have_bfbc_incr = 1;
    if ((fbc_bdforces+l)->y != 0.0) have_bfbc_incr = 1;
#ifndef TWOD
    if ((fbc_bdforces+l)->z != 0.0) have_bfbc_incr = 1;
#endif   
  }
#else
  /* dynamic loading, increment linearly at each timestep */
  if ((ensemble!=ENS_MIK) && (ensemble!=ENS_GLOK) && (ensemble!=ENS_CG)) {
    for (l=0;l<vtypes;l++){
      (fbc_bdf+l)->x = ((fbc_endbforces+l)->x-(fbc_beginbforces+l)->x)/steps_diff;
      (fbc_bdf+l)->y = ((fbc_endbforces+l)->y-(fbc_beginbforces+l)->y)/steps_diff;
#ifndef TWOD
      (fbc_bdf+l)->z = ((fbc_endbforces+l)->z-(fbc_beginbforces+l)->z)/steps_diff;
#endif
    }
    do_bfbc_incr = 1;
  }
#endif /* RELAX */

}

/*****************************************************************************
*
*  update BFBC
*
*****************************************************************************/

void update_bfbc()
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
    if ((is_relaxed) || (bfbc_int > max_fbc_int)) {
        write_ssdef(steps);
        write_ssconfig(steps); /* write config, even when not fully relaxed */
        for (l=0; l<vtypes; l++) fbc_dbf[l] = fbc_dbforces[l];
        bfbc_int = 0;
        do_bfbc_incr = 1;
#ifdef CG
    if (ensemble == ENS_CG) reset_cg();
#endif
    }
    else {
      for (l=0; l<vtypes; l++) *(fbc_bdf+l) = nullv;
      do_bfbc_incr = 0;
      bfbc_int++;
    }
  }
#endif
  /* apply fbc increment if necessary */
  if (do_fbc_incr == 1) {
    for (l=0; l<vtypes; l++){     
      (fbc_bforces+l)->x += (fbc_bdf+l)->x;  
      (fbc_bforces+l)->y += (fbc_bdf+l)->y;
#ifndef TWOD
      (fbc_bforces+l)->z += (fbc_bdf+l)->z;
#endif
    } 
#ifdef GLOK
    if (ensemble==ENS_GLOK)
    {
        reset_glok();
    }
#endif
  }
}

#endif /* BEND */



#ifdef ZAPP
void init_zapp(void)
{
    int         k;
    vektor vectmp;
    total_impuls.x = 0.0;   nactive_vect.x = 0.0;
    total_impuls.y = 0.0;   nactive_vect.y = 0.0;
#ifndef TWOD
    total_impuls.z = 0.0;   nactive_vect.z = 0.0;
#endif

    /* calc total impuls */
   for (k=0; k<NCELLS; ++k) {
      int i;
      cell *p;
      vektor *rest;
      p = CELLPTR(k);
      for (i=0; i<p->n; ++i) {
          rest = restrictions + VSORTE(p,i);
          nactive_vect.x += (int) rest->x;
          nactive_vect.y += (int) rest->y;
#ifndef TWOD
          nactive_vect.z += (int) rest->z;
#endif
          total_impuls.x += IMPULS(p,i,X);
          total_impuls.y += IMPULS(p,i,Y);
#ifndef TWOD
          total_impuls.z += IMPULS(p,i,Z);
#endif
      }
          
   }
#ifdef MPI
   MPI_Allreduce( &total_impuls,&vectmp, DIM, REAL, MPI_SUM, cpugrid);
   total_impuls = vectmp;
   MPI_Allreduce( &nactive_vect,&vectmp, DIM, REAL, MPI_SUM, cpugrid);
   nactive_vect = vectmp;
#endif   
    total_impuls.x = nactive_vect.x == 0 ? 0.0 : total_impuls.x / nactive_vect.x;
    total_impuls.y = nactive_vect.y == 0 ? 0.0 : total_impuls.y / nactive_vect.y;
#ifndef TWOD
    total_impuls.z = nactive_vect.z == 0 ? 0.0 : total_impuls.z / nactive_vect.z;
#endif

    
    if(SPROD(total_impuls,total_impuls)>=zapp_threshold*zapp_threshold)
    {
        for (k=0; k<NCELLS; ++k) {
            int i;
            cell *p;
            vektor *rest;
            p = CELLPTR(k);
            for (i=0; i<p->n; ++i) {
                rest = restrictions + VSORTE(p,i);
                IMPULS(p,i,X) -= total_impuls.x * rest->x;
                IMPULS(p,i,Y) -= total_impuls.y * rest->y;
#ifndef TWOD
                IMPULS(p,i,Z) -= total_impuls.z * rest->z;
#endif
            }
        }
    }
}


void zapp(void)
{
    int         k;
    vektor vectmp;
    total_impuls.x = 0.0;   
    total_impuls.y = 0.0;   
#ifndef TWOD
    total_impuls.z = 0.0;   
#endif

    /* calc total impuls */
   for (k=0; k<NCELLS; ++k) {
      int i;
      cell *p;
      vektor *rest;
      p = CELLPTR(k);
      for (i=0; i<p->n; ++i) {
          total_impuls.x += IMPULS(p,i,X);
          total_impuls.y += IMPULS(p,i,Y);
#ifndef TWOD
          total_impuls.z += IMPULS(p,i,Z);
#endif
      }
          
   }
#ifdef MPI
   MPI_Allreduce( &total_impuls,&vectmp, DIM, REAL, MPI_SUM, cpugrid);
   total_impuls = vectmp;
#endif   
    total_impuls.x = nactive_vect.x == 0 ? 0.0 : total_impuls.x / nactive_vect.x;
    total_impuls.y = nactive_vect.y == 0 ? 0.0 : total_impuls.y / nactive_vect.y;
#ifndef TWOD
    total_impuls.z = nactive_vect.z == 0 ? 0.0 : total_impuls.z / nactive_vect.z;
#endif
    if(SPROD(total_impuls,total_impuls)>=zapp_threshold*zapp_threshold)
    {
        for (k=0; k<NCELLS; ++k) {
            int i;
            cell *p;
            vektor *rest;
            p = CELLPTR(k);
            for (i=0; i<p->n; ++i) {
                rest = restrictions + VSORTE(p,i);
                IMPULS(p,i,X) -= total_impuls.x * rest->x;
                IMPULS(p,i,Y) -= total_impuls.y * rest->y;
#ifndef TWOD
                IMPULS(p,i,Z) -= total_impuls.z * rest->z;
#endif
            }
        }
    }
}
#endif


#ifdef BEND
/******************************************************************************
*
* initialize centers of gravity etc of bending moments
*
******************************************************************************/

void init_bend(void)
{
    int k,j,i;
    real norm_bend_vec, norm_fbc_bforce;
    int tmpivec1[12], tmpivec2[12];
    real tmpvec1[6],tmpvec2[6];
    vektor tmp_force,fbc_bforce, tmp_bend_vec,this_bend_axis;
    int sort;
    
    for (j=0; j<bend_nmoments; j++){
        (bend_origin + j)->x = 0.0;
        (bend_origin + j)->y = 0.0;
        (bend_origin + j)->z = 0.0;
        bend_natomsvtype_origin[j]=0;

        (bend_cog + j)->x = 0.0;
        (bend_cog + j)->y = 0.0;
        (bend_cog + j)->z = 0.0;
        bend_natomsvtype_force[j]=0;
    }
    for (k=0; k<NCELLS; k++) {
        int i;
        cell* p;
        p = CELLPTR(k);
        for (i=0; i<p->n; i++) {
            for (j=0; j<bend_nmoments; j++){
                if(VSORTE(p,i)==bend_vtype_of_origin[j])
                {
                    (bend_origin + j)->x += ORT(p,i,X);
                    (bend_origin + j)->y += ORT(p,i,Y);
                    (bend_origin + j)->z += ORT(p,i,Z);
                    bend_natomsvtype_origin[j]++;
                }
                else if(VSORTE(p,i)==bend_vtype_of_force[j])
                {
                    (bend_cog + j)->x += ORT(p,i,X);
                    (bend_cog + j)->y += ORT(p,i,Y);
                    (bend_cog + j)->z += ORT(p,i,Z);
                    bend_natomsvtype_force[j]++;
                }
            }
        }
    }
   
#ifdef MPI
    for(i=0;i<6;i++)
    {
        tmpivec1[i]= bend_natomsvtype_origin[i];
    }
    for(i=0;i<6;i++)
    {
        tmpivec1[i+6]= bend_natomsvtype_force[i];
    }
    MPI_Allreduce( tmpivec1, tmpivec2,12, MPI_INT, MPI_SUM, cpugrid);
    for(i=0;i<6;i++)
    {
        bend_natomsvtype_origin[i]=tmpivec2[i];
        if(i<bend_nmoments && bend_natomsvtype_origin[i]==0)
            error("bending moment defined without atoms at origin");
    }
    for(i=0;i<6;i++)
    {
        bend_natomsvtype_force[i]=tmpivec2[i+6];
        if(i<bend_nmoments && bend_natomsvtype_origin[i]==0)
            error("bending moment defined without atoms to apply force to");
    }

    for (j=0; j<bend_nmoments; j++){
        tmpvec1[0]=(bend_origin + j)->x;
        tmpvec1[1]=(bend_origin + j)->y;
        tmpvec1[2]=(bend_origin + j)->z;
        tmpvec1[3]=(bend_cog + j)->x;
        tmpvec1[4]=(bend_cog + j)->y;
        tmpvec1[5]=(bend_cog + j)->z;
        MPI_Allreduce(tmpvec1 , tmpvec2, 6, REAL, MPI_SUM, cpugrid);
        (bend_origin + j)->x = tmpvec2[0];
        (bend_origin + j)->y = tmpvec2[1];
        (bend_origin + j)->z = tmpvec2[2];
        (bend_cog + j)->x    = tmpvec2[3];
        (bend_cog + j)->y    = tmpvec2[4];
        (bend_cog + j)->z    = tmpvec2[5];
    }
#endif /* MPI */

    for (j=0; j<bend_nmoments; j++){
        (bend_origin + j)->x /= bend_natomsvtype_origin[j];
        (bend_origin + j)->y /= bend_natomsvtype_origin[j];
        (bend_origin + j)->z /= bend_natomsvtype_origin[j];

        (bend_cog + j)->x /= bend_natomsvtype_force[j];
        (bend_cog + j)->y /= bend_natomsvtype_force[j];
        (bend_cog + j)->z /= bend_natomsvtype_force[j];

        (bend_vec + j)->x  = (bend_cog + j)->x - (bend_origin + j)->x;
        (bend_vec + j)->y  = (bend_cog + j)->y - (bend_origin + j)->y;
        (bend_vec + j)->z  = (bend_cog + j)->z - (bend_origin + j)->z;

        sort = bend_vtype_of_force[j];
        fbc_bforce.x = (fbc_bforces+sort)->x;
        fbc_bforce.y = (fbc_bforces+sort)->y;
        fbc_bforce.z = (fbc_bforces+sort)->z;

        tmp_bend_vec.x = (bend_vec+j)->x;
        tmp_bend_vec.y = (bend_vec+j)->y;
        tmp_bend_vec.z = (bend_vec+j)->z;
        
        norm_bend_vec  = sqrt(SPROD(tmp_bend_vec,tmp_bend_vec));
       
        norm_fbc_bforce = sqrt(SPROD(fbc_bforce,fbc_bforce));

        tmp_bend_vec.x /= norm_bend_vec;
        tmp_bend_vec.y /= norm_bend_vec;
        tmp_bend_vec.z /= norm_bend_vec;

        this_bend_axis.x = (bend_axis + j)->x;
        this_bend_axis.y = (bend_axis + j)->y;
        this_bend_axis.z = (bend_axis + j)->z;
        
        CROSS3D(tmp_bend_vec,this_bend_axis,tmp_force);

        (bend_forces + sort)->x = tmp_force.x * norm_fbc_bforce;
        (bend_forces + sort)->y = tmp_force.y * norm_fbc_bforce;
        (bend_forces + sort)->z = tmp_force.z * norm_fbc_bforce;

        if (myid==0)
        {
            printf("Bending moment %d: bend_forces %d %lf %lf %lf compare to fbc_bforces  %d %lf %lf %lf\n",
                   j,sort,(bend_forces + sort)->x,(bend_forces + sort)->y,(bend_forces + sort)->z,
                   sort,(fbc_bforces + sort)->x,(fbc_bforces + sort)->y,(fbc_bforces + sort)->z);
        }
        
    }
    
    
}
#endif /* BEND */

#ifdef BEND
/******************************************************************************
*
* update direction of bending forces
*
******************************************************************************/

void update_bend(void)
{ 
    int k,j;
    real norm_bend_vec, norm_fbc_bforce;
    int tmpivec1[12], tmpivec2[12];
    real tmpvec1[6],tmpvec2[6];
    vektor tmp_force,fbc_bforce, tmp_bend_vec,this_bend_axis;
    int sort;
    
    for (j=0; j<bend_nmoments; j++){
        (bend_origin + j)->x = 0.0;
        (bend_origin + j)->y = 0.0;
        (bend_origin + j)->z = 0.0;
        (bend_cog + j)->x = 0.0;
        (bend_cog + j)->y = 0.0;
        (bend_cog + j)->z = 0.0;
     }
    for (k=0; k<NCELLS; k++) {
        int i;
        cell* p;
        p = CELLPTR(k);
        for (i=0; i<p->n; i++) {
            for (j=0; j<bend_nmoments; j++){
                if(VSORTE(p,i)==bend_vtype_of_origin[j])
                {
                    (bend_origin + j)->x += ORT(p,i,X);
                    (bend_origin + j)->y += ORT(p,i,Y);
                    (bend_origin + j)->z += ORT(p,i,Z);
                }
                else if(VSORTE(p,i)==bend_vtype_of_force[j])
                {
                    (bend_cog + j)->x += ORT(p,i,X);
                    (bend_cog + j)->y += ORT(p,i,Y);
                    (bend_cog + j)->z += ORT(p,i,Z);
                }
            }
        }
    }
  
#ifdef MPI
     for (j=0; j<bend_nmoments; j++){
    
        tmpvec1[0]=(bend_origin + j)->x;
        tmpvec1[1]=(bend_origin + j)->y;
        tmpvec1[2]=(bend_origin + j)->z;
        tmpvec1[3]=(bend_cog + j)->x;
        tmpvec1[4]=(bend_cog + j)->y;
        tmpvec1[5]=(bend_cog + j)->z;
        MPI_Allreduce(tmpvec1 , tmpvec2, 6, REAL, MPI_SUM, cpugrid);
        (bend_origin + j)->x = tmpvec2[0];
        (bend_origin + j)->y = tmpvec2[1];
        (bend_origin + j)->z = tmpvec2[2];
        (bend_cog + j)->x    = tmpvec2[3];
        (bend_cog + j)->y    = tmpvec2[4];
        (bend_cog + j)->z    = tmpvec2[5];
        }
#endif /* MPI */
    for (j=0; j<bend_nmoments; j++){
        (bend_origin + j)->x /= bend_natomsvtype_origin[j];
        (bend_origin + j)->y /= bend_natomsvtype_origin[j];
        (bend_origin + j)->z /= bend_natomsvtype_origin[j];

        (bend_cog + j)->x /= bend_natomsvtype_force[j];
        (bend_cog + j)->y /= bend_natomsvtype_force[j];
        (bend_cog + j)->z /= bend_natomsvtype_force[j];

     
        
        (bend_vec + j)->x  = (bend_cog + j)->x - (bend_origin + j)->x;
        (bend_vec + j)->y  = (bend_cog + j)->y - (bend_origin + j)->y;
        (bend_vec + j)->z  = (bend_cog + j)->z - (bend_origin + j)->z;

        sort = bend_vtype_of_force[j];
        fbc_bforce.x = (fbc_bforces+sort)->x;
        fbc_bforce.y = (fbc_bforces+sort)->y;
        fbc_bforce.z = (fbc_bforces+sort)->z;

        tmp_bend_vec.x = (bend_vec+j)->x;
        tmp_bend_vec.y = (bend_vec+j)->y;
        tmp_bend_vec.z = (bend_vec+j)->z;
        
        norm_bend_vec  = sqrt(SPROD(tmp_bend_vec,tmp_bend_vec));
       
        norm_fbc_bforce = sqrt(SPROD(fbc_bforce,fbc_bforce));

        tmp_bend_vec.x /= norm_bend_vec;
        tmp_bend_vec.y /= norm_bend_vec;
        tmp_bend_vec.z /= norm_bend_vec;

        this_bend_axis.x = (bend_axis + j)->x;
        this_bend_axis.y = (bend_axis + j)->y;
        this_bend_axis.z = (bend_axis + j)->z;
        
        CROSS3D(tmp_bend_vec,this_bend_axis,tmp_force);

        (bend_forces + sort)->x = tmp_force.x * norm_fbc_bforce;
        (bend_forces + sort)->y = tmp_force.y * norm_fbc_bforce;
        (bend_forces + sort)->z = tmp_force.z * norm_fbc_bforce;

  
    }
    
    
}
#endif /* BEND */


/*****************************************************************************
*
*  check if sample is relaxed
*
*****************************************************************************/

#ifdef RELAX

void check_relaxed(void)
{
    int write_ss=1;
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

 

            
            /* if we are doing quasistatic simulations, we write out     *
             * in case the sample is relaxed or the max. nr of relaxation*
             * steps has been reached. Here we are only writing out if   *
             * the sample is truely relaxed */
#ifdef DEFORM
            if (max_deform_int > 0) write_ss=0; 
#endif
#ifdef HOMDEF
            if (lindef_int > 0) write_ss=0;
#endif
#ifdef FBC
            if (have_fbc_incr) write_ss=0;
#endif
#ifdef EXTPOT
            if (ep_max_int > 0) write_ss=0;
#endif
            if(write_ss==1)
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
#ifdef HOMDEF
    if (lindef_int > 0) stop=0;
#endif
#ifdef FBC
    if (have_fbc_incr) stop=0;
#endif
#ifdef EXTPOT
    if (ep_max_int > 0) stop = 0;
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
#if defined(MSQD) || defined(NMOLDYN) || defined(FEFL)
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
#if defined(MSQD) || defined(NMOLDYN) || defined(FEFL)
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
#if defined(MSQD) || defined(NMOLDYN) || defined(FEFL)
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

  if (!do_press_calc) error("pressure tensor incomplete in calculation of total presstens");
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
