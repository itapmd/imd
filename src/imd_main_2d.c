
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
  real tmp_pot_energy;
  real tmp_kin_energy;
  real dtemp;
  vektor d_pressure;
#if defined(CORRELATE) || defined(MSQD)
  int ref_step = correl_start;
#endif

  if (0==myid) printf( "Starting simulation %d\n", simulation );

#if defined(AND) || defined(NVT) || defined(NPT)
  dtemp = (end_temp - temperature) / (steps_max - steps_min);
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

#ifdef MIKSHEAR
  stepssincelastshear = 0;
#endif

  for (steps=steps_min; steps <= steps_max; ++steps) { 

#ifdef MIKSHEAR
  stepssincelastshear++;
  if ((ensemble==ENS_MIKSHEAR) && (steps > annealsteps) && ((tot_kin_energy/natoms < shear_epsilon) || ((stepssincelastshear % maxshearrelaxsteps) == 0))) {
    shear1step(steps);
    stepssincelastshear = 0;
  }
#endif

#ifdef MPI
    /* we should do this in more appropriate intervals */
    if (((eng_interval != 0) && (0 == steps%eng_interval)) || 
        (steps == steps_min)) setup_buffers();
    send_atoms(FORCE);
    mpi_addtime(&time_comm);
#endif

#ifndef MC
    calc_forces(); 
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

    move_atoms(); 

    if ((steps==steps_min) && (use_curr_temp==1)) {
      temperature = 2 * tot_kin_energy / (DIM * natoms);
      dtemp = (end_temp - temperature) / (steps_max - steps_min);
      use_curr_temp = 0;
    };

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

#ifndef NPPBC
    do_boundaries();    
#endif

#ifdef MPI
    mpi_addtime(&time_calc);
    MPI_Allreduce(&tot_pot_energy,&tmp_pot_energy,1,MPI_REAL,MPI_SUM,cpugrid); 
    tot_pot_energy = tmp_pot_energy; 
    mpi_addtime(&time_comm);
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

#ifdef DISLOC
    if (steps == up_ort_ref) update_ort_ref();
    if ((dem_interval > 0) && (0 == steps%dem_interval)) write_demmaps(steps);
    if ((dsp_interval > up_ort_ref) && (0 == steps%dsp_interval)) write_dspmaps(steps);
#endif

#ifdef USE_SOCKETS
    if ((socket_int>0) && (0==steps%socket_int)) check_socket(steps);
#endif

#if defined(AND) || defined(NVT) || defined(NPT)
    temperature += dtemp;
#endif

#ifdef NPT
    pressure_ext.x += d_pressure.x;
    pressure_ext.y += d_pressure.y;
#endif

#ifdef MPI
    mpi_addtime(&time_io);
#endif

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
  epilogue();
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

#ifndef NOPBC
          /* Consider PBC */
          dcpu.x -= ((int) (dcpu.x / (cpu_dim.x/2)) * cpu_dim.x);
          dcpu.y -= ((int) (dcpu.y / (cpu_dim.y/2)) * cpu_dim.y);
#endif
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
  int i,j,k,l;
  cell *p;
  vektor d;

  /* for each cell in bulk */
  for (i=cellmin.x; i < cellmax.x; ++i)
    for (j=cellmin.y; j < cellmax.y; ++j) {
      p = PTR_2D_V(cell_array, i, j, cell_dim);
        for (l=0; l<p->n; ++l) {

          /* Apply periodic boundaries */
          d.x = -FLOOR(SPRODX(p->ort,l,tbox_x)) * box_x.x;
          d.y = -FLOOR(SPRODX(p->ort,l,tbox_x)) * box_x.y;

          d.x -= FLOOR(SPRODX(p->ort,l,tbox_y)) * box_y.x;
          d.y -= FLOOR(SPRODX(p->ort,l,tbox_y)) * box_y.y;

#ifndef SHOCK
	  p->ort X(l) += d.x;
#endif
	  p->ort Y(l) += d.y;
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
  int r, s, i;
  cell *p;

  /* Set Up Initial Temperature */
  if (do_maxwell) maxwell(temperature);
  do_maxwell=0;

#ifdef FRAC
  if(ensemble==ENS_FRAC) load_sample(); 
#endif

#if defined(NVT) || defined(NPT)
  eta = 0.0;
#endif

#ifdef NPT
  xi.x = 0.0;
  xi.y = 0.0;
#endif
  
#if defined(PULL) || defined(FRAC) || defined(MIKSHEAR)
  /* Atoms with negative numbers are not moved */
  /* loop over all atoms */
  /* for each cell in bulk */
  if ((ensemble==ENS_PULL) || (ensemble==ENS_FRAC) ||
      (ensemble==ENS_MIKSHEAR)) 
  for (r=cellmin.x; r < cellmax.x; ++r)
    for (s=cellmin.x; s < cellmax.y; ++s) {
      p = PTR_2D_V(cell_array, r, s, cell_dim);
      for (i = 0;i < p->n; ++i) {
        /* Make Atom numbers in strip negative */
        if ((p->ort X(i) < strip) || (p->ort X(i) > (box_x.x - strip))) {
          if (0<p->nummer[i]) p->nummer[i] = - p->nummer[i];
        };
      };
    };
#endif
  
}


/*****************************************************************************
*
*  ensemble specific epilogue
*
*****************************************************************************/

void epilogue(void)
{
  int r, s, i;
  cell *p;
  
#if defined(PULL) || defined(FRAC) || defined(MIKSHEAR)
  /* make all atoms mobile again */
  if ((ensemble==ENS_PULL) || (ensemble==ENS_FRAC) || 
      (ensemble==ENS_MIKSHEAR)) 
  for (r=1; r < cell_dim.x-1; ++r)
    for (s=1; s < cell_dim.y-1; ++s) {
      p = PTR_2D_V(cell_array, r, s, cell_dim);
      for (i = 0;i < p->n; ++i) {
        if (0>p->nummer[i]) p->nummer[i] = - p->nummer[i];
      };
    };
#endif

}






