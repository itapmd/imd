/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2007 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_bboost.c -- Routines for the Bond Boost Method, see ....
*
******************************************************************************/

#include "imd.h"


void  init_bboost(void)
{
    int n, k, i;
#ifdef debugLo
    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("print "" check init_bboost start "" checking by Lo! \n");fflush(stdout);
#endif
    /* initializations usually done in main_loop */   
#ifdef EXTPOT
    init_extpot();
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
#if defined(HOMDEF) && defined(RELAX)
    deform_int = 0; 
#endif

    printf("In init_bboost: determining reference position \n");fflush(stdout);
    
    /* find reference positions by energy minimization */
    bb_minimize(0);
    
    
    /* make relaxed position reference position */
    for (k=0; k<NCELLS; ++k) { /* loop over all cells */
        int  i,j, sort;
        cell *p;
        p = CELLPTR(k);
        for (i=0; i<p->n; ++i) { /* loop over all atoms in the cell */
            REFPOSONE(p,i,X) = ORT(p,i,X);
            REFPOSONE(p,i,Y) = ORT(p,i,Y);
            REFPOSONE(p,i,Z) = ORT(p,i,Z);
           
        }
        
    }


  /* update bb_neighbor table cutoff */
  if (NULL==bb_neightab_r2cut) {
     bb_neightab_r2cut = (real *) calloc( ntypes * ntypes, sizeof(real) );
    if (NULL==bb_neightab_r2cut) 
       error("cannot allocate memory for bb_neightab_r2cut");
  }

  for (i=0; i<ntypes*ntypes; i++) {
      bb_neightab_r2cut[i] = bb_rcut * bb_rcut;
  }

  /* loop over all pairs of cells */
  for (n=0; n<nlists; ++n) {
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
#endif

        for (k=npairs[n]; k<npairs[n]; ++k) { //???
            vektor pbc;
            pair *P;
            P = pairs[n]+k;
            pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
            pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
            pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
            do_bb_neightab(cell_array + P->np, cell_array + P->nq, pbc);
        }
    }


} /* end of init_bboost */




/******************************************************************************
*
*  do_bb_neightab - compute neighbor table
*
******************************************************************************/

void do_bb_neightab(cell *p, cell *q, vektor pbc)
{
  int i, j, k;
  int jstart, jend;
  int q_typ, p_typ, column;
  vektor d, tmp_d;
  real *qptr, radius2;

  /* For each atom in first cell */
  for (i=0; i<p->n; ++i) {

      tmp_d.x = ORT(p,i,X) - pbc.x;
      tmp_d.y = ORT(p,i,Y) - pbc.y;
      tmp_d.z = ORT(p,i,Z) - pbc.z;
      p_typ   = SORTE(p,i);

#ifdef TWOD
      jstart  = (((p==q) && (pbc.x==0) && (pbc.y==0))               ? i+1 : 0);
#else
      jstart  = (((p==q) && (pbc.x==0) && (pbc.y==0) && (pbc.z==0)) ? i+1 : 0);
#endif
      qptr    = &ORT(q,jstart,X);
    
    /* For each atom in neighbouring cell */
    for (j = jstart; j < q->n; ++j) {

        q_typ = SORTE(q,j);
      
      /* Calculate distance  */
        d.x = *qptr++ - tmp_d.x;
        d.y = *qptr++ - tmp_d.y;
        d.z = *qptr++ - tmp_d.z;

        column  = p_typ * ntypes + q_typ;
        radius2 = SPROD(d,d);

        if (0==radius2) { char msgbuf[256];
           sprintf(msgbuf,
                "Distance is zero: nrs=%d %d\norte: %f %f %f, %f %f %f\n",
                NUMMER(p,i),NUMMER(q,i),
                ORT(p,i,X), ORT(p,i,Y), ORT(p,i,Z),
                ORT(q,j,X), ORT(q,j,Y), ORT(q,j,Z) );
           error(msgbuf);
        }

      /* make neighbor tables for boost methods */
        if (radius2 <= bb_neightab_r2cut[column]) {
           bb_neightab *bb_neigh;
           real  *tmp_ptr;
           bb_neigh = BBNEIGH(p,i);
           bb_neigh->numref1[bb_neigh->nbondsref1] = NUMMER(q,i);
           tmp_ptr  = &bb_neigh->distref1[bb_neigh->nbondsref1];
           *tmp_ptr = radius2;
           bb_neigh->nbondsref1++;
        }
    } /* for j */
  } /* for i */

}  /* the end of do_bb_neightab */






/* this is the minimization part of the bond boost method */
/* minimization is carried out according to which option you compiled
   and the parameters in the param file
   this is a copy of main_loop 
 */
int bb_minimize(int simulation)
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

  real fnorm2,ekin,epot,delta_epot;

  
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



  
  /* simulation loop */
  for (steps=steps_min; steps <= bb_relaxsteps_max; ++steps) {

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



      

#ifdef NBLIST
    check_nblist();
#else
    fix_cells();  
#endif

#ifdef CG
        if (ensemble == ENS_CG) reset_cg();
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

        /* calculation of forces */
        calc_forces(steps);

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


#ifdef GLOK
    /* "global convergence": set momenta to 0 if P*F < 0 (global vectors) */
    if (ensemble == ENS_GLOK) {
      update_glok();     
    }
#endif

#ifdef TIMING
    imd_start_timer(&time_integrate);
#endif
#if !defined(CBE) || !defined(SPU_INT)
    /* move atoms */
    if (ensemble != ENS_CG) move_atoms(); /* here PxF is recalculated */
#endif
#ifdef TIMING
    imd_stop_timer(&time_integrate);
#endif


#ifdef BEND
    update_bend();
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

#ifdef NBLIST
    check_nblist();
#else
    fix_cells();  
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

  /* clean up the current phase, and clear restart flag */
  imdrestart=0;
  if (0==myid) {
    write_itr_file(-1, steps_max,"");
    if (0==myrank) printf( "End of simulation %d\n", simulation );
  }  
  return finished;
}
    
