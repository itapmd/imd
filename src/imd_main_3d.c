/******************************************************************************
*
* imd_main_3d.c -- main loop, three dimensions
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"


/*****************************************************************************
*
*  main_loop
*
*****************************************************************************/

void main_loop(void)
{
  int steps;
  real tmp_pot_energy;
  real tmp_kin_energy;
  int i,j,k;
  cell *p;
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
  d_pressure.z = (pressure_end.z - pressure_ext.z) / (steps_max - steps_min);
#endif

#ifdef EAM
  /* memory allocation */
  eam_rho   =  calloc(natoms,sizeof(real));
  if(NULL   == eam_rho) error("Cannot allocate memory for eam_rho");
  eam_ij    =  calloc(natoms*eam_len,sizeof(integer));
  if(NULL   == eam_ij)  error("Cannot allocate memory for eam_ij");
  eam_dij_x =  calloc(natoms*eam_len,sizeof(real));
  if(NULL   == eam_dij_x) error("Cannot allocate memory for eam_dij_x");
  eam_dij_y =  calloc(natoms*eam_len,sizeof(real));
  if(NULL   == eam_dij_y) error("Cannot allocate memory for eam_dij_y");
  eam_dij_z =  calloc(natoms*eam_len,sizeof(real));
  if(NULL   == eam_dij_z) error("Cannot allocate memory for eam_dij_z");
#endif /* EAM */

#ifdef TTBP
  /* memory allocation */
  ttbp_ij   =  calloc(natoms*ttbp_len*2,sizeof(integer));
  if(NULL   == ttbp_ij)  error("Can't allocate memory for ttbp_ij");
  ttbp_j    =  calloc((natoms+1)*ttbp_len,sizeof(real));
  if(NULL   == ttbp_j) error("Can't allocate memory for ttbp_j");
  ttbp_force=  calloc((natoms+1)*3,sizeof(real));
  if(NULL   == ttbp_force)  error("Can't allocate memory for ttbp_force");
#endif /* TTBP */

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
  if (maxshearrelaxsteps == 0) maxshearrelaxsteps = 32767;
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
#ifdef SAVEMEM
    send_atoms_by_cell();
#else
    /* we should do this in more appropriate intervals */
    if (((eng_interval != 0) && (0 == steps%eng_interval)) || 
        (steps == steps_min)) setup_buffers();
    mpi_addtime(&time_io);
#ifdef AR
    send_atoms_ar();
#else
    send_atoms(FORCE);
#endif
#endif
    mpi_addtime(&time_comm_force);
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

#ifdef MPI
#ifdef AR    
    send_forces();
    mpi_addtime(&time_comm_ar);
#endif
#endif

    move_atoms();

#if defined(AND) || defined(NVT) || defined(NPT)
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
      d_pressure.z = (pressure_end.z-pressure_ext.z) / (steps_max-steps_min);
      use_curr_pressure = 0;
    };
#endif

#ifdef NPT_axial
    if ((steps==steps_min) && (ensemble==ENS_NPT_AXIAL) && 
        (use_curr_pressure==1)) {
      pressure_ext = stress;
      d_pressure.x = (pressure_end.x-pressure_ext.x) / (steps_max-steps_min);
      d_pressure.y = (pressure_end.y-pressure_ext.y) / (steps_max-steps_min);
      d_pressure.z = (pressure_end.z-pressure_ext.z) / (steps_max-steps_min);
      use_curr_pressure = 0;
    };
#endif

#ifndef NOPBC
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
    pressure_ext.z += d_pressure.z;
#endif

#ifdef MPI
    mpi_addtime(&time_io);
#endif

    /* fix_cells redistributes atoms across the cpus. Putting at the bottom 
       of the force loop enables us to calculate/write properties locally */
#ifdef SAVEMEM
#ifdef MPI
    /* Deallocate buffer cells each timestep to save memory */
    for (i=0; i < cell_dim.x; ++i)
        for (j=0; j < cell_dim.y; ++j)
            for (k=0; k < cell_dim.z; ++k) {
              p = PTR_3D_V(cell_array, i, j, k, cell_dim); 
            /* Check for buffer cell */
              if ((0==i) || (0==j) || (0==k) ||
                  (i == cell_dim.x-1) ||
                  (j == cell_dim.y-1) ||
                  (k == cell_dim.z-1))
		/* and deallocate */
		alloc_cell(p, 0);
            };
    fix_cells_by_cell();
#endif
#else
    fix_cells();  
#endif

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
       printf("cells too small - start afresh with larger cells!\n");
       break;
    }
#endif
#endif

  }

  /* clean up the current phase, and clear restart flag */
  epilogue();
  restart=0;

#ifdef MPI
  mpi_addtime(&time_comm);
#endif

#ifdef EAM
  free(eam_rho);
  free(eam_ij);
  free(eam_dij_x);
  free(eam_dij_y);
  free(eam_dij_z);
#endif

#ifdef TTBP
  free(ttbp_ij);
  free(ttbp_j);
  free(ttbp_force);
#endif

  if (0==myid) printf( "End of simulation %d\n", simulation );
  
}


/******************************************************************************
*
*  fix_cells
*
*  check if each atom is in the correct cell and on the correct CPU 
*  move atoms the have left their cells in the last timestep
*
*  this also uses Plimpton's comm scheme
*
******************************************************************************/

void fix_cells(void)

{
  int i,j,k,l;
  cell *p, *q;
  ivektor coord, to_coord;
  int to_cpu;

#ifdef MPI
  empty_mpi_buffers();
  empty_buffer_cells();
#endif

  /* for each cell in bulk */
  for (i=cellmin.x; i < cellmax.x; ++i)
    for (j=cellmin.y; j < cellmax.y; ++j)
      for (k=cellmin.z; k < cellmax.z; ++k) {

	p = PTR_3D_V(cell_array, i, j, k, cell_dim);

	/* loop over atoms in cell */
	l=0;
	while( l<p->n ) {

#ifndef MPI
          coord = cell_coord(p->ort X(l),p->ort Y(l),p->ort Z(l));
          q = PTR_3D_VV(cell_array,coord,cell_dim);
          /* if it's in the wrong cell, move it to the right cell */
          if  (p != q) 
            move_atom(coord,p,l); 
          else
            ++l;
#else
	  coord = local_cell_coord(p->ort X(l),p->ort Y(l),p->ort Z(l));
 	  /* see if atom is in wrong cell */
	  if ((coord.x == i) && (coord.y == j) && (coord.z == k)) {
            l++;
          } 
          else {

            to_cpu = cpu_coord( global_cell_coord( coord ));
            /* west */
            if ((to_cpu==nbwest) || (to_cpu==nbnw)  || (to_cpu==nbws) ||
                (to_cpu==nbuw  ) || (to_cpu==nbunw) || (to_cpu==nbuws)||
                (to_cpu==nbdw  ) || (to_cpu==nbdwn) || (to_cpu==nbdsw))
                copy_one_atom( &send_buf_west,  p, l);
            
            /* east */
            else if
                ((to_cpu==nbeast) || (to_cpu==nbse)  || (to_cpu==nben) ||
                 (to_cpu==nbue  ) || (to_cpu==nbuse) || (to_cpu==nbuen)||
                 (to_cpu==nbde  ) || (to_cpu==nbdes) || (to_cpu==nbdne))
                copy_one_atom( &send_buf_east,  p, l);
                        
            /* south  */
            else if ((to_cpu==nbsouth) || (to_cpu==nbus)  || (to_cpu==nbds))
                copy_one_atom( &send_buf_south,  p, l);
                        
            /* north  */
            else if ((to_cpu==nbnorth) || (to_cpu==nbun)  || (to_cpu==nbdn))
                copy_one_atom( &send_buf_north,  p, l);
            
            /* down  */
            else if (to_cpu==nbdown)
                copy_one_atom( &send_buf_down,  p, l);
            
            /* up  */
            else if (to_cpu==nbup)
                copy_one_atom( &send_buf_up  ,  p, l);

            /* atom is on my cpu */
            else if (to_cpu==myid)
                move_atom(coord, p, l);
            
            else error("Atom jumped multiple CPUs");
                        
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
* do_boundaries
*
* Apply periodic boundaries to all atoms
* Could change so that only cells on surface to some work
*
******************************************************************************/

void do_boundaries(void)
{
  int i,j,k,l;
  cell *p;
  vektor d;

  /* for each cell in bulk */
  for (i=cellmin.x; i < cellmax.x; ++i)
    for (j=cellmin.y; j < cellmax.y; ++j)
      for (k=cellmin.z; k < cellmax.z; ++k) {
	p = PTR_3D_V(cell_array, i, j, k, cell_dim);
	for (l=0; l<p->n; ++l) {

	  /* Apply periodic boundaries */
          d.x = -FLOOR(SPRODX(p->ort,l,tbox_x)) * box_x.x;
          d.y = -FLOOR(SPRODX(p->ort,l,tbox_x)) * box_x.y;
          d.z = -FLOOR(SPRODX(p->ort,l,tbox_x)) * box_x.z;

          d.x -= FLOOR(SPRODX(p->ort,l,tbox_y)) * box_y.x;
          d.y -= FLOOR(SPRODX(p->ort,l,tbox_y)) * box_y.y;
          d.z -= FLOOR(SPRODX(p->ort,l,tbox_y)) * box_y.z;

          d.x -= FLOOR(SPRODX(p->ort,l,tbox_z)) * box_z.x;
          d.y -= FLOOR(SPRODX(p->ort,l,tbox_z)) * box_z.y;
          d.z -= FLOOR(SPRODX(p->ort,l,tbox_z)) * box_z.z;

#ifndef SHOCK
          p->ort X(l) += d.x;
#endif
          p->ort Y(l) += d.y;
          p->ort Z(l) += d.z;
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
  volume   = box_x.x * ( box_y.y * box_z.z - box_y.z * box_z.y)
           - box_x.y * ( box_y.x * box_z.z - box_y.z * box_z.x)
           + box_x.z * ( box_y.x * box_z.y - box_y.x * box_z.y);
  pressure = ( 2.0 * tot_kin_energy + virial ) / ( 3 * volume );
}


/*****************************************************************************
*
*  ensemble specific initializations
*
*****************************************************************************/

void init(void)
{
  int i, r, s, t;
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
  xi.z = 0.0;
#endif
  
#if defined(PULL) || defined(FRAC) || defined(MIKSHEAR)
  /* Atoms with negative numbers are immobile */
  if ((ensemble==ENS_PULL) || (ensemble==ENS_FRAC) ||
      (ensemble==ENS_MIKSHEAR)) 
  for (r=cellmin.x; r < cellmax.x; ++r)
    for (s=cellmin.y; s < cellmax.y; ++s)
      for (t=cellmin.z; t < cellmax.z; ++t) {
	p = PTR_3D_V(cell_array, r, s, t, cell_dim);
        for (i = 0; i < p->n; ++i) {
          /* Make Atom numbers in strip negative */
	  if ((p->ort X(i) < strip) || (p->ort X(i) > (box_x.x - strip))) {
            if (0<p->nummer[i]) p->nummer[i] = - p->nummer[i];
          }
        }
  }
#endif
  
}


/*****************************************************************************
*
*  ensemble specific epilogue
*
*****************************************************************************/

void epilogue(void)
{
  int r, s, t, i;
  cell *p;
  
#if defined(PULL) || defined(FRAC) || defined(MIKSHEAR)
  /* make all atoms mobile again */
  if ((ensemble==ENS_PULL) || (ensemble==ENS_FRAC)) 
  for (r=cellmin.x; r < cellmax.x; ++r)
    for (s=cellmin.y; s < cellmax.y; ++s)
      for (t=cellmin.z; t < cellmax.z; ++t) {
	p = PTR_3D_V(cell_array, r, s, t, cell_dim);
        for (i = 0;i < p->n; ++i) {
          if (0>p->nummer[i]) p->nummer[i] = - p->nummer[i];
        }
  }
#endif

}
