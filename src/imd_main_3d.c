
/******************************************************************************
*
* imd_main_3d.c -- main loop, three dimensions
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/*#define ATNR 2826   easier debuging... */

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
#ifdef GLOK
  real old_PxF=1.0;
#endif
#if defined(CORRELATE) || defined(MSQD)
  int ref_step = correl_start;
#endif

#ifdef FBC
#ifdef MIK
  int nofbcsteps=0;
#endif 
  int l;
  vektor nullv={0.0,0.0,0.0};
  vektor temp_df;
  vektor *fbc_df;
  fbc_df = (vektor *) malloc(vtypes*DIM*sizeof(real));
  if (NULL==fbc_df)
    error("Can't allocate memory for fbc_df\n");
#endif

  if (0==myid) printf( "Starting simulation %d\n", simulation );

#if defined(AND) || defined(NVT) || defined(NPT)
  dtemp = (end_temp - temperature) / (steps_max - steps_min);
#endif



#ifdef FBC
#ifndef MIK
/* dynamic loading, increment linearly each timestep */
  for (l=0;l<vtypes;l++){
    temp_df.x = (((fbc_endforces+l)->x) - ((fbc_beginforces+l)->x))/(steps_max - steps_min);
    temp_df.y = (((fbc_endforces+l)->y) - ((fbc_beginforces+l)->y))/(steps_max - steps_min);
    temp_df.z = (((fbc_endforces+l)->z) - ((fbc_beginforces+l)->z))/(steps_max - steps_min);
    
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
  d_pressure.z = (pressure_end.z - pressure_ext.z) / (steps_max - steps_min);
#endif

#ifdef EAM
  /* memory allocation */
  eam_rho   =  calloc((natoms+1),sizeof(real));
  if(NULL   == eam_rho) error("Cannot allocate memory for eam_rho");
  eam_ij    =  calloc((natoms+1)*eam_len,sizeof(real));
  if(NULL   == eam_ij) error("Cannot allocate memory for eam_ij");
  eam_dij_x =  calloc((natoms+1)*eam_len,sizeof(real));
  if(NULL   == eam_dij_x) error("Cannot allocate memory for eam_dij_x");
  eam_dij_y =  calloc((natoms+1)*eam_len,sizeof(real));
  if(NULL   == eam_dij_y) error("Cannot allocate memory for eam_dij_y");
  eam_dij_z =  calloc((natoms+1)*eam_len,sizeof(real));
  if(NULL   == eam_dij_z) error("Cannot allocate memory for eam_dij_z");
#endif /* EAM */

#if defined(CORRELATE) || defined(MSQD)
  init_correl(ncorr_rmax,ncorr_tmax);
#endif

  calc_properties();

  /* initializations for the current simulation phase, if not yet done */
  if (0==restart) init();

#ifdef DEFORM
  deform_int = 0;
#endif

  /* simulation loop */
  for (steps=steps_min; steps <= steps_max; ++steps) {

#ifdef FBC
#ifdef MIK  
    /* just increment the force if under threshold of e_kin or after 
       waitsteps and after annelasteps */
    temp_df.x = 0.0;
    temp_df.y = 0.0;
    temp_df.z = 0.0;  
    for (l=0; l<vtypes; l++) *(fbc_df+l) = temp_df;

    if (steps > fbc_annealsteps) {
      nofbcsteps++; 

      if ((2.0*tot_kin_energy/nactive < fbc_ekin_threshold) ||
          (nofbcsteps==fbc_waitsteps)) 
      {
        nofbcsteps=0;
        for (l=0;l<vtypes;l++) {
          (fbc_df+l)->x = (fbc_dforces+l)->x;
          (fbc_df+l)->y = (fbc_dforces+l)->y;
          (fbc_df+l)->z = (fbc_dforces+l)->z;

          /* MIK affects the total impuls, especially in inhomogenous samples,
             so we set the velocities to 0 befor each force increment 

             for (k=0; k<ncells; ++k) {
               p = cell_array + CELLS(k);
               for (i=0; i<p->n; ++i) {
                 p->impuls X(i) = 0.0;
                 p->impuls Y(i) = 0.0;
                 p->impuls Z(i) = 0.0;
               }
             }
          */
        }
        /* to avoid an imense total impuls and to see more of the phase space,
	   we do some AND steps here.... */
#ifdef MIKAND
        if ((nofbcsteps!=0) && (0==nofbcsteps%tmp_interval)) 
          maxwell(temperature);
#endif
      }
    }
#endif /* MIK */
#endif /* FBC */

#ifdef HOMDEF
    if ((exp_interval > 0) && (0 == steps%exp_interval)) expand_sample();
    if ((hom_interval > 0) && (0 == steps%hom_interval)) shear_sample();
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
    for (l=0; l<vtypes; l++){ 
      /* fbc_df=0 if MIK && ekin> ekin_threshold */
      (fbc_forces+l)->x += (fbc_df+l)->x;  
      (fbc_forces+l)->y += (fbc_df+l)->y;
      (fbc_forces+l)->z += (fbc_df+l)->z;
    } 
#ifdef ATNR
    printf(" step: %d \n",steps); 
    fflush(stdout);
#endif
#endif

#ifdef MPI
#ifdef TIMING
    imd_start_timer(&time_force_comm);
#endif
#ifdef SAVEMEM
    send_atoms_by_cell();
#else
    /* we should do this in more appropriate intervals */
    if (((eng_interval != 0) && (0 == steps%eng_interval)) || 
        (steps == steps_min)) setup_buffers();
    send_cells(copy_cell,pack_cell,unpack_cell);
#endif
#ifdef TIMING
    imd_stop_timer(&time_force_comm);
#endif
#endif

#ifndef MC
#ifdef TIMING
    imd_start_timer(&time_force_calc);
#endif
    calc_forces();
#ifdef TIMING
    imd_stop_timer(&time_force_calc);
#endif
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
#ifdef TIMING
    imd_start_timer(&time_force_comm);
#endif
    send_forces(add_forces,pack_forces,unpack_forces);
#ifdef TIMING
    imd_stop_timer(&time_force_comm);
#endif
#endif
#endif

#ifdef SMIK /* have relaxsteps nvesteps before mik */
    if (steps<relaxsteps)
      move_atoms_nve();
    else
      move_atoms();
#else 

#ifdef GLOK 
    /* "globale konvergenz": set impulses=0 if P*F <0 (global vectors) */
    if (PxF<0.0) {
      for (k=0; k<ncells; ++k) {
        p = cell_array + CELLS(k);
        for (i=0; i<p->n; ++i) {
          p->impuls X(i) = 0.0;
          p->impuls Y(i) = 0.0;
          p->impuls Z(i) = 0.0;
        }
      }
      if (myid==0) write_properties(steps); 
    }
    /* properties as they were after setting p=0 and 
       calculating ekin_ and p. p should then = f*dt
       therefore PxF >0 and f*f = 2*M*ekin
    */
    if ((myid==0) && (old_PxF<0.0)) write_properties(steps); 
    old_PxF = PxF;
#endif

    move_atoms(); /* here PxF is recalculated */

#endif /* not SMIK */

#if defined(AND) || defined(NVT) || defined(NPT)
    if ((steps==steps_min) && (use_curr_temp==1)) {
#ifdef UNIAX
      temperature = 2.0 * tot_kin_energy / (5.0 * natoms);
#else
      temperature = 2.0 * tot_kin_energy / (DIM * natoms);
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
      d_pressure.z = (pressure_end.z-pressure_ext.z) / (steps_max-steps_min);
      use_curr_pressure = 0;
    }
#endif

#ifdef NPT_axial
    if ((steps==steps_min) && (ensemble==ENS_NPT_AXIAL) && 
        (use_curr_pressure==1)) {
      pressure_ext = stress;
      d_pressure.x = (pressure_end.x-pressure_ext.x) / (steps_max-steps_min);
      d_pressure.y = (pressure_end.y-pressure_ext.y) / (steps_max-steps_min);
      d_pressure.z = (pressure_end.z-pressure_ext.z) / (steps_max-steps_min);
      use_curr_pressure = 0;
    }
#endif

#if defined(AND) || defined(NVT) || defined(NPT)
    temperature += dtemp;
#endif

#ifdef NVX
    tran_Tleft   += dtemp;
    tran_Tright  -= dtemp;
#endif


#ifdef NPT
    pressure_ext.x += d_pressure.x;
    pressure_ext.y += d_pressure.y;
    pressure_ext.z += d_pressure.z;
#endif

    calc_properties(); 

    /* Periodic I/O */
#ifdef TIMING
    imd_start_timer(&time_io);
#endif
    if ((rep_interval > 0) && (0 == steps%rep_interval)) write_config(steps);
    if ((eng_interval > 0) && (0 == steps%eng_interval) && (0==myid)) 
       write_properties(steps);
    if ((dis_interval > 0) && (0 == steps%dis_interval)) write_distrib(steps);
    if ((pic_interval > 0) && (0 == steps%pic_interval)) write_pictures(steps);

#ifdef EFILTER  /* just print atoms if in an energy-window */ 
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
      if (0==press_dim.x) 
        write_config_select(steps/press_interval,"press",write_cell_press);
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

    /* fix_cells redistributes atoms across the cpus. Putting at the bottom 
       of the force loop enables us to calculate/write properties locally */
#if (defined(MPI) && defined(SAVEMEM))
    /* Deallocate buffer cells each timestep to save memory */
    for (i=0; i < cell_dim.x; ++i)
        for (j=0; j < cell_dim.y; ++j)
            for (k=0; k < cell_dim.z; ++k) {
              p = PTR_3D_V(cell_array, i, j, k, cell_dim); 
              /* check for buffer cell and deallocate */
              if ((0==i) || (0==j) || (0==k) ||
                  (i == cell_dim.x-1) ||
                  (j == cell_dim.y-1) ||
                  (k == cell_dim.z-1))
		alloc_cell(p, 0);
            }
    fix_cells_by_cell();
#else
    fix_cells();  
#endif

#ifdef NPT
    /* revise cell division if necessary */
    if (revise_cell_division==1) {
       init_cells();
       fix_cells();
    }
#endif

  }

  /* clean up the current phase, and clear restart flag */
  restart=0;

#ifdef EAM
  free(eam_rho);
  free(eam_ij);
  free(eam_dij_x);
  free(eam_dij_y);
  free(eam_dij_z);
#endif

  if (0==myid) printf( "End of simulation %d\n", simulation );
  
}


/******************************************************************************
*
*  fix_cells
*
*  check if each atom is in the correct cell and on the correct CPU 
*  move atoms that have left their cells in the last timestep
*
*  this also uses Plimpton's comm scheme
*
******************************************************************************/

void fix_cells(void)

{
  int i,j,k,l;
  cell *p, *q;
  ivektor coord, lcoord;
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
	  lcoord = local_cell_coord(p->ort X(l),p->ort Y(l),p->ort Z(l));
 	  /* see if atom is in wrong cell */
	  if ((lcoord.x == i) && (lcoord.y == j) && (lcoord.z == k)) {
            l++;
          } 
          else {

            /* global cell coord and CPU */
            coord  = cell_coord(p->ort X(l),p->ort Y(l),p->ort Z(l));
            to_cpu = cpu_coord(coord);

            /* atom is on my cpu */
            if (to_cpu==myid)
                move_atom(lcoord, p, l);
            
            /* west */
            else if 
               ((to_cpu==nbwest) || (to_cpu==nbnw)  || (to_cpu==nbws) ||
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
      p->ort Z(l) += i * box_x.z;
    }

    /* PBC in y direction */
    if (pbc_dirs.y==1)
    for (l=0; l<p->n; ++l) {
      i = -FLOOR(SPRODX(p->ort,l,tbox_y));
      p->ort X(l) += i * box_y.x;
      p->ort Y(l) += i * box_y.y;
      p->ort Z(l) += i * box_y.z;
    }

    /* PBC in z direction */
    if (pbc_dirs.z==1)
    for (l=0; l<p->n; ++l) {
      i = -FLOOR(SPRODX(p->ort,l,tbox_z));
      p->ort X(l) += i * box_z.x;
      p->ort Y(l) += i * box_z.y;
      p->ort Z(l) += i * box_z.z;
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
#ifdef UNIAX
  pressure = ( 2.0 / 5.0 * tot_kin_energy + virial / 3.0 ) / volume ;
#else
  pressure = ( 2.0 * tot_kin_energy + virial ) / ( 3.0 * volume );
#endif
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

