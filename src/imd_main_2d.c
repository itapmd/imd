
/******************************************************************************
*
* imd_main_2d.c -- main loop in two dimensions
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
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
  int steps;
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

  if (0==myid) printf( "Starting simulation %d\n", simulation );

#if defined(AND) || defined(NVT) || defined(NPT) || defined(STM)
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
#endif

#if defined(CORRELATE) || defined(MSQD)
  init_correl(ncorr_rmax,ncorr_tmax);
#endif

  calc_properties();

#ifdef MPI
  mpi_addtime(&time_setup);
#endif

  /* initializations for the current simulation phase, if not yet done */
  if (0==restart) init();

#ifdef DEFORM
  deform_int = 0;
#endif

  for (steps=steps_min; steps <= steps_max; ++steps) { 

#ifdef FBC
#ifdef MIK  
/* just increment the force if under threshold of e_kin or after waitsteps
  and after annelasteps */
  temp_df.x = 0.0;
  temp_df.y = 0.0;
   for (l=0;l<vtypes;l++)
      *(fbc_df+l) = temp_df;

 if (steps > fbc_annealsteps)
 {
   nofbcsteps++; 
   if((tot_kin_energy/nactive < fbc_ekin_threshold) ||
        (nofbcsteps==fbc_waitsteps)) 
     {
      nofbcsteps=0;
      for (l=0;l<vtypes;l++)
         *(fbc_df+l) = *(fbc_dforces+l) ;
      /* MIK affects the total impuls, especially in inhomogenous samples,
	 so we set the velocities to 0 befor each force increment 
      for (k=0; k<ncells; ++k) {
	p = cell_array + CELLS(k);
	for (i=0; i<p->n; ++i) {
	  p->impuls X(i) = 0.0;
	  p->impuls Y(i) = 0.0;
	}
      }
      */
     }
 }

#endif
#endif

#ifdef HOMDEF
    if ((exp_interval > 0) && (0 == steps%exp_interval)) expand_sample();
    if ((hom_interval > 0) && (0 == steps%hom_interval)) shear_sample();
#endif
#ifdef DEFORM
    if (steps > annealsteps) {
      deform_int++;
      if ((tot_kin_energy/nactive < ekin_threshold) || 
          (deform_int==max_deform_int)) {
        deform_sample();
        deform_int=0;
      }
    }
#endif
#ifdef FBC
 for (l=0;l<vtypes;l++){               /* fbc_forces is already initialised with beginforces */
   (fbc_forces+l)->x += (fbc_df+l)->x; /* fbc_df=0 if MIK && ekin> ekin_threshold */ 
   (fbc_forces+l)->y += (fbc_df+l)->y;
  } 
#endif
#ifdef MPI
    /* we should do this in more appropriate intervals */
    if (((eng_interval != 0) && (0 == steps%eng_interval)) || 
        (steps == steps_min)) setup_buffers();
    send_atoms_force();
    mpi_addtime(&time_comm);
#endif

#ifndef MC
    calc_forces(); 
#endif

#ifdef MPI
    mpi_addtime(&time_calc);
#endif

#ifdef DISLOC
    if ((steps==reset_Epot_step) && (calc_Epot_ref==1)) reset_Epot_ref();
#endif

#if defined(CORRELATE) || defined(MSQD)
    if ((steps >= correl_start) && ((steps < correl_end) || (correl_end==0))) {
      int istep = steps - correl_start;
      if (istep % correl_ts == 0) correlate(steps,ref_step,istep/correl_ts);
      if ((correl_int != 0) && (steps-ref_step+1 >= correl_int)) 
        ref_step += correl_int;
    }
#endif
#ifdef GLOK /* "globale konvergenz": set impulses=0 if P*F <0 (global vectors) */
    if (PxF<0.0)
      for (k=0; k<ncells; ++k) {
	p = cell_array + CELLS(k);
	for (i=0; i<p->n; ++i) {
	  p->impuls X(i) = 0.0;
	  p->impuls Y(i) = 0.0;
	  p->impuls Z(i) = 0.0;
	}
      }
#endif
    move_atoms(); 

#if defined(AND) || defined(NVT) || defined(NPT) || defined(STM)
    if ((steps==steps_min) && (use_curr_temp==1)) {
      temperature = 2 * tot_kin_energy / (DIM * natoms);
      dtemp = (end_temp - temperature) / (steps_max - steps_min);
      use_curr_temp = 0;
    };
#endif

#ifdef NPT_iso
    if ((steps==steps_min) && (ensemble==ENS_NPT_ISO) && 
        (use_curr_pressure==1)) {
      pressure_ext.x = pressure;
      d_pressure.x = (pressure_end.x-pressure_ext.x) / (steps_max-steps_min);
      d_pressure.y = (pressure_end.y-pressure_ext.y) / (steps_max-steps_min);
      use_curr_pressure = 0;
    };
#endif

#ifdef NPT_axial
    if ((steps==steps_min) && (ensemble==ENS_NPT_AXIAL) && 
        (use_curr_pressure==1)) {
      pressure_ext = stress;
      d_pressure.x = (pressure_end.x-pressure_ext.x) / (steps_max-steps_min);
      d_pressure.y = (pressure_end.y-pressure_ext.y) / (steps_max-steps_min);
      use_curr_pressure = 0;
    };
#endif

#if defined(AND) || defined(NVT) || defined(NPT) || defined(STM)
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

#ifdef MPI
    mpi_addtime(&time_calc);
#endif

    calc_properties(); 

#ifdef MPI
    mpi_addtime(&time_calc);
#endif

    /* Periodic I/O */
    if ((rep_interval > 0) && (0 == steps%rep_interval)) write_config(steps);
    if ((eng_interval > 0) && (0 == steps%eng_interval) && (0==myid)) 
       write_properties(steps);
    if ((dis_interval > 0) && (0 == steps%dis_interval)) write_distrib(steps);
    if ((pic_interval > 0) && (0 == steps%pic_interval)) write_pictures(steps);

#ifdef EFILTER  /* just print atoms in an energy-window */ 
    if ((efrep_interval > 0) && (0 == steps%efrep_interval)) 
       write_config_select(steps/efrep_interval,"ef",write_cell_ef);
#endif
#ifdef DISLOC
    if (steps == up_ort_ref) update_ort_ref();
    if ((dem_interval > 0) && (0 == steps%dem_interval)) 
       write_config_select(steps,"dem",write_cell_dem);
    if ((dsp_interval > up_ort_ref) && (0 == steps%dsp_interval)) 
       write_config_select(steps,"dsp",write_cell_dsp);
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
      if (0==press_dim.x) write_press_atoms(steps);
      else write_press_dist(steps);
    }
#endif
#endif

#ifdef USE_SOCKETS
    if ((socket_int>0) && (0==steps%socket_int)) check_socket(steps);
#endif

#ifdef MPI
    mpi_addtime(&time_io);
#endif

    do_boundaries();    

    fix_cells();

#ifdef NPT
#ifdef DYNAMIC_CELLS
    /* revise cell division if necessary */
    if (revise_cell_division==1) {
      init_cells();
      fix_cells();
    }  
#else
    /* check if cells have become too small */
    if (myid == 0)
    if (cells_too_small) {
       /* write config if not yet written */
       if ((rep_interval == 0) || (0 != steps%rep_interval)){
          write_config(steps);
          write_properties(steps);
       };
       printf("cells too small -- start afresh with larger cells!\n");
       break;
    }
#endif
#endif

  }

  /* clean up the current phase, and clear restart flag */
  restart=0;

  if (0==myid) printf( "End of simulation %d\n", simulation );

}


/******************************************************************************
*
*  fix_cells
*
*  check if each atom is in the correct cell and on the correct CPU 
*  move atoms that have left their cells in the last timestep
*
******************************************************************************/

void fix_cells(void)
{
  int i,j,l;
  cell *p, *q;
  ivektor coord, dcpu, to_coord;

#ifdef MPI
  empty_mpi_buffers();
#endif

  /* for each cell in bulk */
  for (i=cellmin.x; i < cellmax.x; ++i)
    for (j=cellmin.y; j < cellmax.y; ++j) {

      p = PTR_2D_V(cell_array, i, j, cell_dim);

      /* loop over atoms in cell */
      l=0;
      while( l < p->n ) {

#ifndef MPI
        coord = cell_coord(p->ort X(l),p->ort Y(l));
        q = PTR_2D_VV(cell_array,coord,cell_dim);
        /* if it's in the wrong cell, move it to the right cell */
        if (p != q) 
          move_atom(coord,p,l); 
        else 
          ++l;
#else
        coord = local_cell_coord(p->ort X(l),p->ort Y(l));
	/* see if atom is in wrong cell */
        if ((coord.x == i) && (coord.y == j)) {
          l++;
        } else {

          /* Calculate distance on CPU grid */
          to_coord = cpu_coord_v( cell_coord( p->ort X(l),p->ort Y(l) ));
          dcpu.x = to_coord.x - my_coord.x;
          dcpu.y = to_coord.y - my_coord.y;

          /* Consider PBC */
          if (pbc_dirs.x == 1) {
            if (cpu_dim.x == 1) dcpu.x = 0; 
            else dcpu.x -= ((int) (dcpu.x / (cpu_dim.x/2)) * cpu_dim.x);
          }
          if (pbc_dirs.y == 1) {
            if (cpu_dim.y == 1) dcpu.y = 0;
            else dcpu.y -= ((int) (dcpu.y / (cpu_dim.y/2)) * cpu_dim.y);
          }

          /* Check, if atom is on my cpu */
          /* If not, copy into send buffer else move to correct cell */
          if      ( 0 < dcpu.x ) copy_one_atom( &send_buf_west,  p, l);
          else if ( 0 > dcpu.x ) copy_one_atom( &send_buf_east,  p, l);
          else if ( 0 < dcpu.y ) copy_one_atom( &send_buf_south, p, l); 
          else if ( 0 > dcpu.y ) copy_one_atom( &send_buf_north, p, l); 
          else { /* atom is on my cpu */
            move_atom(coord, p, l);
          }
        }
#endif /* MPI */
      }
    }

#ifdef MPI
  /* send border cells to neighbbours */
  send_atoms(ATOMS);
#endif

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
      p->ort X(l) += i * box_x.x;
      p->ort Y(l) += i * box_x.y;
    }

    /* PBC in y direction */
    if (pbc_dirs.y==1)
    for (l=0; l<p->n; ++l) {
      i = -FLOOR(SPRODX(p->ort,l,tbox_y));
      p->ort X(l) += i * box_y.x;
      p->ort Y(l) += i * box_y.y;
    }

  }
}


/******************************************************************************
*
*  calc_properties
*
******************************************************************************/

void calc_properties(void)
{
  /* volume and pressure */
  volume   = box_x.x * box_y.y - box_x.y * box_y.x;
  pressure = ( 2.0 * tot_kin_energy + virial ) / ( 2 * volume );
}


/*****************************************************************************
*
*  ensemble specific initializations
*
*****************************************************************************/

void init(void)
{
  /* Set Up Initial Temperature */
  if (do_maxwell) maxwell(temperature);
  do_maxwell=0;

#ifdef FRAC
  if(ensemble==ENS_FRAC) load_sample(); 
#endif
  
}







