
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
* globals.h -- Global Variables for the imd Package
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/


/* MAIN is defined only once in the main module */
#ifdef MAIN 
#define EXTERN        /* define Variables in main */
#define INIT(data) =data

#ifdef TWOD
#define nullvektor  { 0.0, 0.0 }
#define nullivektor { 0, 0 }
#define einsvektor  { 1.0, 1.0 }
#define einsivektor { 1, 1 }
#define parteinsivektor { 1, 0 }
#define nullsymtensor { 0.0, 0.0, 0.0 }
#else
#define nullvektor  { 0.0, 0.0, 0.0 }
#define nullivektor { 0, 0, 0 }
#define einsvektor  { 1.0, 1.0, 1.0 }
#define einsivektor { 1, 1, 1 }
#define parteinsivektor { 1, 1, 0 }
#define nullsymtensor { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
#endif

#define nullvektor2d { 0.0, 0.0 }
#define nullivektor2d { 0, 0 }
#define nullvektor3d { 0.0, 0.0, 0.0 }
#define nullivektor3d { 0, 0, 0 }

#define nullbuffer  {NULL, 0, 0 }

#define zero10       {0,0,0,0,0,0,0,0,0,0}

#define zero45       {0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0}

#define zero55       {0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0}

#define zero550      {0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0}

#else
#define EXTERN extern /* declare them extern otherwise */
#define INIT(data)
#endif

/* Where all the data lives */
EXTERN cell *cell_array INIT(NULL);      /* array of cells */
EXTERN integer *cells INIT(NULL);        /* list if inner cell indices */
#ifdef TWOD
EXTERN pair *pairs[9];                   /* arrays of cell pairs */  
EXTERN int  npairs[9], npairs2[9];       /* number of cell pairs */
#else
EXTERN pair *pairs[27];                  /* arrays of cell pairs */  
EXTERN int  npairs[27], npairs2[27];     /* number of cell pairs */
#endif
EXTERN int ncells, nallcells INIT(0);    /* number of cells */
EXTERN int nlists;                       /* number of cell pair lists */
EXTERN ivektor cell_dim;                 /* dimension of cell array (per cpu)*/
EXTERN ivektor global_cell_dim;          /* dimension of cell array */

/* Boundary Conditions */
#ifdef EPITAX
EXTERN ivektor pbc_dirs INIT(parteinsivektor);
#else
EXTERN ivektor pbc_dirs INIT(einsivektor); /* directions with pbc */
#endif

EXTERN int vtypes INIT(0);   /* vtypes = ntypes + ntypes*(types of forces) */
EXTERN vektor *restrictions INIT(NULL);  /* directions the atom is allowed to move in */

#ifdef FBC    /* FBC uses the virtual atom types, */
EXTERN vektor *fbc_forces;       /* each vtype has its force */
EXTERN vektor *fbc_beginforces;  /* begin & endvalues for linear interpolation */ 
#if defined(MIK) ||defined (CG)                       /* or ... */
EXTERN vektor *fbc_dforces;      /* each vtype has its force increment */
EXTERN real   fbc_ekin_threshold INIT(0.0);/* threshold for ekin */ 
EXTERN int    fbc_waitsteps INIT(1); /* between increase of forces */
EXTERN int    fbc_annealsteps INIT(1); /* times before 1. + df */
EXTERN vektor *fbc_endforces;   
#else
EXTERN vektor *fbc_endforces;   
#endif

#endif

/* Global bookkeeping */
EXTERN int is_big_endian;        /* 1 if big endian, 0 if little endian */
EXTERN int natoms  INIT(0);      /* Total number of atoms */
EXTERN int nactive INIT(0);      /* number of transl. degrees of freedom */
EXTERN int nactive_rot INIT(0);  /* number of rot. degrees of freedom */
EXTERN int ntypes INIT(0);       /* Total number of different atom types */
EXTERN int ntypepairs INIT(0);   /* Total number of different types pairs */
EXTERN int ntypetriples INIT(0); /* Total number of different types triples 
				    symmetric in the last two indices */
EXTERN int *num_sort INIT(NULL); /* number of atoms for each real type */
EXTERN int *num_vsort INIT(NULL); /* number of atoms for each virtual type */
EXTERN int steps INIT(0);        /* number of current MD step */
EXTERN int steps_max INIT(0);    /* Maximum number of MD steps */
EXTERN int steps_min INIT(0);    /* starting step nr for current phase */
EXTERN int restart INIT(0);      /* file number for restart */
EXTERN int checkpt_int INIT(0);  /* Period of checkpoints */
EXTERN int eng_int INIT(0);      /* Period of data output */
EXTERN int force_int INIT(0);    /* Period of force file writing */
EXTERN str255 progname; /* Name of current executable argv[0] */
EXTERN ivektor cellmin; /* Minimum index of local cells (1 with MPI, 0 otherwise) */
EXTERN ivektor cellmax; /* Maximum index of local cells  */
EXTERN int use_header INIT(1);   /* shall a header be written */

/* controlling distribution output */
EXTERN int dist_Ekin_flag        INIT(0); /* write Ekin dists? */
EXTERN int dist_Epot_flag        INIT(0); /* write Epot dists? */
EXTERN int dist_Ekin_long_flag   INIT(0); /* write Ekin_long dists? */
EXTERN int dist_Ekin_trans_flag  INIT(0); /* write Ekin_trans dists? */
EXTERN int dist_shock_shear_flag INIT(0); /* write shock shear dists? */
EXTERN int dist_press_flag       INIT(0); /* write press dists? */
EXTERN int dist_presstens_flag   INIT(0); /* write presstens dists? */
EXTERN int dist_int              INIT(0); /* Period of distribution writes */
EXTERN ivektor dist_dim          INIT(einsivektor); /* resolution of dist */
EXTERN vektor  dist_ur           INIT(nullvektor);  /* lower left  corner */
EXTERN vektor  dist_ll           INIT(nullvektor);  /* upper right corner */

#ifdef EFILTER
EXTERN int  ef_checkpt_int INIT(0);  /* Period of ef writes */
EXTERN real lower_e_pot INIT(0.0);   /* lower end of energy window */
EXTERN real upper_e_pot INIT(0.0);   /* upper end of energy window */
#endif
#ifdef CLONE
EXTERN int  nclones INIT(0);         /* number of periodic clones */
#endif

#ifdef NBFILTER
EXTERN int nb_checkpt_int INIT(0);   /* Period of nb writes */
EXTERN int *lower_nb_cut INIT(NULL); /* lower number of neighbours  */
EXTERN int *upper_nb_cut INIT(NULL); /* upper number of neighbours  */
#endif
/* data for generated structures */
EXTERN ivektor box_param INIT(nullivektor);  /* box parameters */
EXTERN real box_unit INIT(1.0);              /* lattice parameter */
EXTERN real *masses INIT(NULL);              /* masses */

/* The simulation box and its inverse */
EXTERN vektor box_x  INIT(nullvektor);
EXTERN vektor box_y  INIT(nullvektor);
#ifndef TWOD
EXTERN vektor box_z  INIT(nullvektor);
#endif
EXTERN vektor tbox_x INIT(nullvektor);
EXTERN vektor tbox_y INIT(nullvektor);
#ifndef TWOD
EXTERN vektor tbox_z INIT(nullvektor);
#endif
EXTERN vektor     height INIT(nullvektor);
EXTERN vektor max_height INIT(nullvektor);
EXTERN vektor min_height INIT(nullvektor);

/* Filenames */
EXTERN char outbuf[OUTPUT_BUF_SIZE] INIT("\0");  /* output buffer */
EXTERN str255 infilename INIT("\0");    /* Input File */
EXTERN str255 itrfilename INIT("\0");   /* initial itr-file */
EXTERN str255 outfilename INIT("\0");   /* Output File */
EXTERN str255 potfilename INIT("\0");   /* Potential */
EXTERN str255 paramfilename INIT("\0"); /* Parameter File */

/* Parameters for displacement map */
EXTERN str255 reffilename INIT("\0");   /* Parameter File */

/* pointers to files that are kept open */
EXTERN FILE  *eng_file INIT(NULL);      /* pointer to .eng file  */
EXTERN FILE *msqd_file INIT(NULL);      /* pointer to .msqd file */
EXTERN int flush_int INIT(50);          /* flush .eng and .msqd files
                                           every flush_int writes */

/* Parameters for pictures */
EXTERN int pic_int INIT(0);                 /* Period for picture output */
EXTERN vektor2d ecut_kin INIT(nullvektor2d);/* Kin. Energy interval for pictures */
EXTERN vektor2d ecut_pot INIT(nullvektor2d);/* Pot. Energy interval for pictures */   
EXTERN vektor pic_scale INIT(nullvektor);   /* Scale factor x/y for pictures     */
EXTERN vektor pic_ll  INIT(nullvektor);     /* lower left (front) corner */
EXTERN vektor pic_ur  INIT(nullvektor);     /* upper right (back) corner */
EXTERN ivektor   pic_res INIT(nullivektor); /* number of pixels in x/y dir.*/
EXTERN int       pic_type INIT(0);          /* picture type 0/1 */
EXTERN int       nsmear   INIT(5);          /* smearing radius in pixels */
#ifndef TWOD
EXTERN vektor3d view_dir INIT(nullvektor);  /* view direction */
EXTERN vektor3d view_pos INIT(nullvektor);  /* view position */
EXTERN int      projection INIT(0);         /* projection type 0/1 */
#endif

/* MD Stuff */
EXTERN real timestep INIT(0.0);
EXTERN real tot_pot_energy INIT(0.0);
EXTERN real tot_kin_energy INIT(0.0);
EXTERN real Ekin_old INIT(0.0);
EXTERN real Erot_old INIT(0.0);
EXTERN real pressure INIT(0.0);
EXTERN real volume INIT(0.0);
EXTERN real volume_init INIT(0.0);
EXTERN real virial INIT(0.0);
EXTERN real temperature INIT(0.0);
EXTERN int  use_curr_temp INIT(0);  /* which starting temp to use (flag) */
EXTERN int  do_maxwell INIT(0.0);
EXTERN long seed INIT(0);           /* seed for random number generator */

/* square of global force vector f=(f1.x, f1.y,...,fn.z) */
EXTERN real fnorm INIT(0.0);  
/* scalar product of global force and momentum vectors */ 
EXTERN real PxF INIT(0.0);
/* Einstein frequency is similar as fnorm, but divided by the masses */ 
EXTERN real omega_E INIT(0.0);

/* Potential Table */
EXTERN pot_table_t pair_pot;         /* potential data structure */
EXTERN real cellsz INIT(0);          /* minimal cell diameter */
EXTERN int  initsz INIT(10);         /* initial number of atoms in cell */
EXTERN int  incrsz INIT(10);         /* increment of number of atoms in cell */


/* MPI housekeeping */
EXTERN int myid INIT(0);                  /* Who am I? (0 if serial) */
EXTERN int num_cpus INIT(1);              /* How many cpus are there */
EXTERN int parallel_output INIT(0);       /* Flag for parallel output */
EXTERN int parallel_input  INIT(1);       /* Flag for parallel input */
#ifdef MPI
EXTERN ivektor cpu_dim INIT(nullivektor); /* Dimensions of CPU-Array */
#else
EXTERN ivektor cpu_dim INIT(einsivektor); /* Dimensions of CPU-Array */
#endif
#ifdef MPI
EXTERN int binc INIT(0);                  /* buffer size per atom */
EXTERN int *cpu_ranks INIT(0);            /* Mapping of coords to ranks */
EXTERN cell buf_one_atom;                 /* Buffer that holds one Atom */
EXTERN MPI_Comm cpugrid;                  /* Cartesian MPI communicator */
EXTERN ivektor my_coord;                  /* Cartesian coordinates of cpu */

/* Send and Receive buffers */
EXTERN msgbuf send_buf_east  INIT(nullbuffer);
EXTERN msgbuf send_buf_west  INIT(nullbuffer);
EXTERN msgbuf send_buf_north INIT(nullbuffer);
EXTERN msgbuf send_buf_south INIT(nullbuffer);
EXTERN msgbuf send_buf_up    INIT(nullbuffer);
EXTERN msgbuf send_buf_down  INIT(nullbuffer);
EXTERN msgbuf recv_buf_east  INIT(nullbuffer);
EXTERN msgbuf recv_buf_west  INIT(nullbuffer);
EXTERN msgbuf recv_buf_north INIT(nullbuffer);
EXTERN msgbuf recv_buf_south INIT(nullbuffer);
EXTERN msgbuf recv_buf_up    INIT(nullbuffer);
EXTERN msgbuf recv_buf_down  INIT(nullbuffer);

/* Neighbours */
EXTERN int nbwest, nbeast, nbnorth, nbsouth, nbup, nbdown; /* Faces */
EXTERN int nbnw, nbws, nbse, nben,                         /* Edges */
           nbun, nbuw, nbus, nbue,
           nbdn, nbde, nbds, nbdw;
EXTERN int nbunw, nbuws, nbuse, nbuen,                     /* Corners */
           nbdne, nbdes, nbdsw, nbdwn;
#endif

/* Timers */
EXTERN imd_timer time_total;
EXTERN imd_timer time_setup;
EXTERN imd_timer time_main;
EXTERN imd_timer time_io;
EXTERN imd_timer time_fft;
EXTERN imd_timer time_fft_plan;

/* Parameters for the various ensembles */

#ifdef AND
EXTERN int tempintv INIT(0);        /* Interval in which the thermostat */
                                    /* kicks (in timesteps) */
#endif

#if defined(NVT) || defined(NPT) || defined(STM)
EXTERN real eta INIT(0.0);          /* Nose-Hoover heat bath variable */
EXTERN real tau_eta INIT(0.0);      /* Nose-Hoover heat bath relax. time */
EXTERN real isq_tau_eta INIT(0.0);  /* isq_tau_eta : inverse of tau_eta^2 */
#ifdef UNIAX
/* Nose-Hoover heat bath variable for rotational motion: */
EXTERN real eta_rot INIT(0.0);      
/* inverse of tau_eta_rot^2, with tau_eta_rot the Nose-Hoover */
/* heat bath relaxation time for rotational motion: */
EXTERN real isq_tau_eta_rot INIT(0.0);  
#endif
#endif

/* virial tensor */
EXTERN real vir_xx INIT(0.0), vir_yy INIT(0.0), vir_zz INIT(0.0);
EXTERN real vir_yz INIT(0.0), vir_zx INIT(0.0), vir_xy INIT(0.0);
/* diagonal of stress tensor */
EXTERN real stress_x INIT(0.0), stress_y INIT(0.0), stress_z INIT(0.0);
EXTERN real dyn_stress_x INIT(0.0), dyn_stress_y INIT(0.0), dyn_stress_z INIT(0.0);

#ifdef NPT
EXTERN real   isq_tau_xi INIT(0.0); /* inverse of tau_xi^2 */
EXTERN real   cell_size_tolerance INIT(0.05);
EXTERN vektor pressure_ext INIT(nullvektor);
EXTERN vektor pressure_end INIT(nullvektor);
EXTERN int    use_curr_pressure INIT(0);  /* which starting pressure to use */
EXTERN vektor xi INIT(nullvektor), xi_old INIT(nullvektor);
#endif

EXTERN real end_temp INIT(0.0);        /* Temperature and at of simulation */
#if defined(DEFORM)
EXTERN real   ekin_threshold INIT(0.0); /* threshold for ekin */    
EXTERN int    annealsteps INIT(0);      /* number of annealing steps */    
#endif

#if defined(GLOK)
EXTERN real   glok_ekin_threshold INIT(100.0); /* threshold for ekin */  
EXTERN int    glok_annealsteps INIT(0);      /* number of annealing steps */
#endif
#ifdef DEFORM
EXTERN int    deform_int INIT(0);       /* counting steps between 2 shears */
EXTERN int    max_deform_int INIT(0);   /* max. steps between 2 shear steps */
EXTERN real   fnorm_threshold INIT(0.0);/* threshold for fnorm */    
EXTERN vektor *deform_shift;            /* shift for each vtype */
EXTERN vektor *deform_shear;            /* shear for each vtype */
EXTERN vektor *deform_base;             /* base point for shear deformation */
EXTERN int    *shear_def;               /* shear flag for each vtype */
#endif

#if defined(DEFORM) || defined(HOMDEF)
EXTERN real   deform_size INIT(1.0);    /* scale factor for deformation */
#endif

#ifdef HOMDEF
EXTERN int    exp_interval INIT(0);       /* period of expansion steps */
EXTERN vektor expansion INIT(einsvektor); /* expansion factors in x/y/z-dir */
EXTERN int    hom_interval INIT(0);       /* period of homshear steps */
EXTERN vektor2d shear_factor INIT(nullvektor2d);/* shear factor in x,y-dir */
EXTERN int    lindef_interval INIT(0);    /* period of linear deform. steps */
EXTERN vektor lindef_x INIT(nullvektor);  /* \               */
EXTERN vektor lindef_y INIT(nullvektor);  /*  |  linear      */
#ifndef TWOD                              /*   > deformation */
EXTERN vektor lindef_z INIT(nullvektor);  /*  |  matrix      */
#endif                                    /* /               */
#endif

#ifdef SLLOD
EXTERN vektor shear_rate   INIT(nullvektor); /* shear rate as a vector */
#ifndef TWOD
EXTERN vektor shear_rate2  INIT(nullvektor); /* shear rate as a vector */
#endif
#endif

#if defined(FRAC) || defined(STM) || defined(FTG) 
EXTERN vektor2d   center  INIT(nullvektor2d); /* center of stadium */
EXTERN vektor2d  stadium  INIT(nullvektor2d); /* half axes of Damping stadium */
EXTERN vektor2d  stadium2 INIT(nullvektor2d); /* half axes where the max. 
						damping factor  is reached  */
#endif

#if defined(FRAC) || defined(STM) 
EXTERN real   E_kin_stadium INIT(0.0);      /* kin energy of the stadium */
EXTERN int    n_stadium INIT(0);            /* number of transl. degrees 
					       of freedom in the stadium */
#endif

#if defined(FRAC) || defined(FTG)  
EXTERN real gamma_damp  INIT(0.0);        /* Damping factor */
EXTERN real gamma_min  INIT(0.0);         /* minimal Damping factor */
EXTERN real gamma_bar   INIT(0.0);        /* Damping prefactor */
EXTERN int  dampingmode INIT(0);          /* damping mode  */
                                          /* 0: ramped viscous damping       */
                                          /* 1: Nose-Hoover           */
EXTERN real delta_ftg   INIT(10.0);        /* free parameter in calculation 
					    of local temperature */

EXTERN real dotepsilon INIT(0.0);         /* strain rate for crack loading */
EXTERN real dotepsilon0 INIT(0.0);        /* initial strain rate */
EXTERN int  expansionmode INIT(1);        /* mode for loading */
#endif

#ifdef FINNIS
EXTERN real delta_finnis INIT(10.0);       /* time constant in finnis */ 
EXTERN real zeta_0       INIT(0.0);        /* prefactor for finnis*/
#endif


#ifdef FRAC
EXTERN real E_kin_damp INIT(0.0);         /* weighted !!  kin dampenergy  */
EXTERN real sum_f INIT(0.0);              /* Sum of stadium function */
#endif

#ifdef FTG
EXTERN real  Tleft   INIT(0.0);          /* Temperature of the left  wall*/
EXTERN real  Tright  INIT(0.0);          /* Temperature of the right wall*/
EXTERN int   nslices INIT(1);            /* Number of slices*/
EXTERN int   nslices_Left  INIT(1);      /* Number of slices with Tleft */
EXTERN int   nslices_Right INIT(1);      /* Number of slices with Tright*/
EXTERN int  *ninslice;                   /* Number of atoms in slice*/
EXTERN real *gamma_ftg;                  /* Damping prefactor for slices*/
EXTERN real *E_kin_ftg;                  /* kin energy of the slices */
#endif

#ifdef CG
/* Parameters used by CG */
EXTERN real   cg_threshold INIT(0.0);   /* threshold for cg */    
EXTERN int    cg_maxsteps  INIT(0);     /* max number of cg steps */    
EXTERN int    linmin_maxsteps  INIT(0); /* max number of linmin steps */    
EXTERN real   linmin_tol  INIT(0.0);    /* tolerance between 2 linmin steps */
EXTERN real   linmin_dmax  INIT(0.0);   /* max. search steps in linmin  */ 
#ifndef DEFORM                          /* no double definition */
EXTERN int    annealsteps INIT(0);      /* number of annealing steps */    
#endif
/* Variables needed by CG */
EXTERN int    cgsteps   INIT(0);        /* current nr of CG step */
EXTERN real   fmax2     INIT(0.0);      /* max. force comp. ^2  */    
EXTERN real   old_cgval  INIT(0.0);      /* old value to minimize */ 
EXTERN real   gg         INIT(0.0);      /* see Num. Rec. p.320 */       
EXTERN real   dgg        INIT(0.0);      /* see Num. Rec. p.320 */    
EXTERN real   cg_gamma   INIT(0.0);      /* see Num. Rec. p.320 */    
#endif

#ifdef SNAPSHOT
EXTERN int sscount INIT(0);
#endif

#ifdef DISLOC
EXTERN int  dem_int INIT(0);          /* Period of dem output */
EXTERN int  dsp_int INIT(0);          /* Period of dsp output */
EXTERN int  up_ort_ref INIT(0);       /* time to update ort_ref ? */
EXTERN real min_dpot INIT(1.0);       /* difference for dem */
EXTERN real min_dsp2 INIT(0.0);       /* minimal square displacement in .dsp */
EXTERN int  reset_Epot_step INIT(0);  /* step at which Epot_ref is computed */
EXTERN int  calc_Epot_ref INIT(0);    /* flag whether to compute Epot_ref */
EXTERN int  Epot_diff INIT(1);        /* flag whether to write Epot_diff */
#endif

#ifdef AVPOS
EXTERN int avpos_start INIT(0);       /* Start time for avpos */
EXTERN int avpos_end INIT(0);         /* End time for avpos */
EXTERN int avpos_int INIT(0);         /* Period of avp output ==0 */
EXTERN int avpos_res INIT(0);         /* Period of coordinate addition */
#ifdef NPT
EXTERN vektor av_box_x  INIT(nullvektor); /* Average of box vectors */
EXTERN vektor av_box_y  INIT(nullvektor); 
#ifndef TWOD
EXTERN vektor av_box_z  INIT(nullvektor);
#endif
#endif
#endif

#ifdef ATDIST
EXTERN float  *atdist INIT(NULL);                /* atoms distribution */
EXTERN ivektor atdist_dim    INIT(nullivektor);  /* dimension of atoms_dist */
EXTERN ivektor atdist_per_ll INIT(nullivektor);  /* ll of periodic ext. */
EXTERN ivektor atdist_per_ur INIT(einsivektor);  /* ur of periodic ext. */
EXTERN vektor  atdist_scale  INIT(nullvektor );  /* scale of atoms dist bins */
EXTERN vektor  atdist_ll     INIT(nullvektor);   /* lower left of dist */
EXTERN vektor  atdist_ur     INIT(nullvektor);   /* upper right of dist */
EXTERN int atdist_int INIT(0);        /* interval between atoms dist updates */
EXTERN int atdist_start INIT(0);      /* start step of atoms distribution */
EXTERN int atdist_end INIT(0);        /* stop step of atoms distribution */
EXTERN int atdist_size INIT(0);       /* size of atoms distribution */
EXTERN int atdist_pos_int INIT(0);    /* period of position writes */
EXTERN real atdist_phi INIT(0.0);     /* rotation angle around z-axis */
#endif

#ifdef DIFFPAT
EXTERN float *diffdist INIT(NULL);              /* atoms distribution */
EXTERN float *diffpat  INIT(NULL);              /* diffraction pattern */
EXTERN float diffpat_weight[10] INIT(zero10);   /* diffraction strengths */
EXTERN ivektor diffpat_dim   INIT(nullivektor);  /* dimension of diffdist */
EXTERN vektor  diffpat_scale INIT(nullvektor );  /* scale of atoms dist bins */
EXTERN vektor  diffpat_ur INIT(nullvektor );   /* upper right corner of dist */
EXTERN vektor  diffpat_ll INIT(nullvektor );   /* lower left  corner of dist */
EXTERN int diffpat_int INIT(0);     /* interval between atoms dist updates */
EXTERN int diffpat_start INIT(0);   /* start step of atoms distribution */
EXTERN int diffpat_end INIT(0);     /* stop step of atoms distribution */
EXTERN int diffpat_size INIT(0);    /* size of atoms distribution */
EXTERN fftwf_plan diffpat_plan;     /* plan for FFT */
#endif

#ifdef ORDPAR
#define nullvektor4d { 0.0, 0.0, 0.0, 0.0 }
EXTERN real op_weight[2][2] INIT(nullvektor4d);
EXTERN real op_r2_cut[2][2] INIT(nullvektor4d);
#endif

/* Global data for MSQD & correlation */
#if defined(MSQD) || defined(CORRELATE)
EXTERN real *msqd INIT(NULL);        /* (local) array of mean square disp. */
EXTERN real *msqd_global INIT(NULL); /* global array of mean square disp. */
EXTERN real *msqdv INIT(NULL);       /* the same for virtual types */
EXTERN real *msqdv_global INIT(NULL);/*         - '' -             */
EXTERN int  msqd_ntypes INIT(1);     /* write msqd for real types */
EXTERN int  msqd_vtypes INIT(0);     /* write msqd for virtual types */
EXTERN int  correl_int INIT(0);      /* repeat interval for correlation */
EXTERN int  correl_start INIT(0);    /* start time for correlation */
EXTERN int  correl_end INIT(0);      /* end time for correlation */
EXTERN int  correl_ts INIT(0);  /* sampling time interval for correlation */
EXTERN int  ncorr_rmax INIT(0); /* dimension of histogram in r domain */
EXTERN int  ncorr_tmax INIT(1); /* dimension of histogram in t domain */
EXTERN real GS_rcut INIT(0);    /* cutoff radius for correlation data writes */
#endif

/* Global data for correlation */
#ifdef CORRELATE
EXTERN real inv_dr;                /* inverse of step size delta r */
EXTERN int correl_omode INIT(1);   /* output mode for histogram */
EXTERN integer ***GS INIT(NULL);   /* histogram array for self correlation */
#endif

#ifdef NVX
EXTERN real dTemp_start   INIT(0.0);   /* deviation of starting heat bath 
                                          temperature from mean temp. */
EXTERN real dTemp_end 	  INIT(0.0);   /* deviation of final heat bath 
                                          temperature from mean temp. */ 
EXTERN real tran_Tleft 	  INIT(0.0);   /* temperature of the left=hot wall */
EXTERN real tran_Tright   INIT(0.0);   /* temperature of the right=cold  wall*/
EXTERN real heat_cond     INIT(0.0);   /* heat conductivity */
#endif
#ifdef RNEMD
EXTERN real heat_transfer INIT(0.0);   /* total (integrated) heat transfer */
EXTERN int  exch_interval INIT(0);     /* interval between particle exchange */
#endif
#ifdef TRANSPORT
EXTERN int  tran_interval INIT(0);     /* Intervalle der Temperaturaufz.*/
EXTERN int  tran_nlayers  INIT(0);     /* number of layers*/
#endif

#ifdef STRESS_TENS
EXTERN sym_tensor tot_presstens INIT(nullsymtensor);/* global pressure tens. */
EXTERN int press_int INIT(0);    /* interval for writing the pressure tensor */
EXTERN int do_press_calc INIT(0);           /* flag whether to do press calc */
#endif

/* I/O via sockets */
#ifdef USE_SOCKETS
EXTERN int server_socket INIT(1);          /* socket mode: client or server */
EXTERN int socket_int INIT(1);               /* interval for reading socket */
EXTERN unsigned short client_port INIT(0);        /* client port for socket */
EXTERN unsigned short server_port INIT(31050);    /* server port for socket */
EXTERN char display_host[256] INIT("");      /* name of controlling machine */
EXTERN int  use_socket_window INIT(0);  /* flag for using a window to write */
EXTERN vektor socketwin_ll  INIT(nullvektor);  /* lower left (front) corner */
EXTERN vektor socketwin_ur  INIT(nullvektor);  /* upper right (back) corner */
#endif

EXTERN int  have_potfile INIT(0); 
#ifdef PAIR
/* analytically defined pair potentials */
EXTERN real r_cut_lin[55] INIT(zero55); 
EXTERN real r_cut [10][10]; 
EXTERN real r2_cut[10][10]; 
EXTERN real r_begin[55] INIT(zero55);
EXTERN real pot_res[55] INIT(zero55);
EXTERN int  have_pre_pot INIT(0); 
/* Lennard-Jones */
EXTERN real lj_epsilon_lin[55] INIT(zero55);
EXTERN real lj_epsilon[10][10];
EXTERN real lj_sigma_lin[55] INIT(zero55);
EXTERN real lj_sigma[10][10];
EXTERN real lj_shift[10][10];
/* Morse */
EXTERN real morse_epsilon_lin[55] INIT(zero55);
EXTERN real morse_epsilon[10][10];
EXTERN real morse_sigma_lin[55] INIT(zero55);
EXTERN real morse_sigma[10][10];
EXTERN real morse_alpha_lin[55] INIT(zero55);
EXTERN real morse_alpha[10][10];
EXTERN real morse_shift[10][10];
/* Buckingham */
EXTERN real buck_a_lin[55] INIT(zero55);
EXTERN real buck_a[10][10];
EXTERN real buck_c_lin[55] INIT(zero55);
EXTERN real buck_c[10][10];
EXTERN real buck_sigma_lin[55] INIT(zero55);
EXTERN real buck_sigma[10][10];
EXTERN real buck_shift[10][10];
/* harmonic potential for shell model */
EXTERN real spring_const[45] INIT(zero45);
EXTERN real spring_cst[10][10];
#endif

#ifdef EAM2
EXTERN pot_table_t embed_pot;                     /* embedding energy table  */
EXTERN pot_table_t rho_h_tab;                     /* electron transfer table */
EXTERN str255 eam2_emb_E_filename INIT("\0");     /* embedding energy file   */
EXTERN str255 eam2_at_rho_filename INIT("\0");    /* electron transfer file  */
#endif

#ifdef MEAM
EXTERN int have_embed_potfile INIT(0);
EXTERN int have_pre_embed_pot INIT(0);
EXTERN pot_table_t embed_pot;                     /* embedding energy table  */
EXTERN str255 meam_emb_E_filename INIT("\0");     /* embedding energy file   */
EXTERN int have_eldensity_file INIT(0);
EXTERN pot_table_t el_density;
EXTERN str255 meam_eldensity_filename INIT("\0");
EXTERN real meam_r2_cut[10][10];
EXTERN int  meam_t_average INIT(0);
EXTERN real meam_t1[10] INIT(zero10);
EXTERN real meam_t2[10] INIT(zero10);
EXTERN real meam_t3[10] INIT(zero10);
EXTERN real meam_f0[10] INIT(zero10);
EXTERN real meam_r0[10] INIT(zero10);
EXTERN real invmeam_r0[10];
EXTERN real meam_beta0[10] INIT(zero10);
EXTERN real meam_beta1[10] INIT(zero10);
EXTERN real meam_beta2[10] INIT(zero10);
EXTERN real meam_beta3[10] INIT(zero10);
EXTERN real meam_rcut_lin[55] INIT(zero55);
EXTERN real meam_rcut[10][10];
EXTERN real meam_deltar_lin[55] INIT(zero55);
EXTERN real meam_deltar[10][10];
EXTERN real meam_cmax_lin[550] INIT(zero550);
EXTERN real meam_cmax[10][10][10];
EXTERN real meam_cmin_lin[550] INIT(zero550);
EXTERN real meam_cmin[10][10][10];
EXTERN real meam_a[10] INIT(zero10);
EXTERN real meam_e[10] INIT(zero10);
EXTERN real meam_rho0[10] INIT(zero10);
EXTERN real invmeam_rho0[10];
EXTERN vektor nullvek INIT(nullvektor);
#endif

#ifdef TTBP
EXTERN str255 ttbp_potfilename INIT("\0"); /* TTBP smoothing potential file */
EXTERN pot_table_t smooth_pot;    /* TTBP smoothing potential */
EXTERN real ttbp_constant[10];	  /* constants; less than 11 atom types! */
EXTERN real ttbp_sp[10];          /* constants; less than 11 atom types! */
#endif

#ifdef STIWEB
/* Parameters for Stillinger-Weber potential */
/* Not more than 10 atom types!              */
EXTERN real stiweb_a[55] INIT(zero55); 
EXTERN real sw_a[10][10];
EXTERN real stiweb_b[55] INIT(zero55);
EXTERN real sw_b[10][10];
EXTERN real stiweb_p[55] INIT(zero55);
EXTERN real sw_p[10][10];
EXTERN real stiweb_q[55] INIT(zero55);
EXTERN real sw_q[10][10];
EXTERN real stiweb_a1[55] INIT(zero55);
EXTERN real sw_a1[10][10];
EXTERN real stiweb_de[55] INIT(zero55);
EXTERN real sw_de[10][10];
EXTERN real sw_2_a1[10][10];
EXTERN real stiweb_a2[55] INIT(zero55);
EXTERN real sw_a2[10][10];
EXTERN real stiweb_la[550] INIT(zero550);
EXTERN real sw_la[10][10][10];
EXTERN real stiweb_ga[55] INIT(zero55);
EXTERN real sw_ga[10][10];
#endif

#ifdef TERSOFF
/* Parameters for Tersoff potential  */
/* Not more than 10 atom types!      */
EXTERN real ters_r_cut[55] INIT(zero55);  /* cutoff^2 */
EXTERN real ter_r_cut[10][10];
EXTERN real ter_r2_cut[10][10];
EXTERN real ters_r0[55] INIT(zero55);        
EXTERN real ter_r0[10][10];
EXTERN real ters_a[55] INIT(zero55);       
EXTERN real ter_a[10][10];
EXTERN real ters_b[55] INIT(zero55);
EXTERN real ter_b[10][10];
EXTERN real ters_la[55] INIT(zero55);
EXTERN real ter_la[10][10];
EXTERN real ters_mu[55] INIT(zero55);
EXTERN real ter_mu[10][10];
EXTERN real ters_chi[45] INIT(zero45);
EXTERN real ter_chi[10][10];
EXTERN real ters_om[45] INIT(zero45);
EXTERN real ter_om[10][10];
EXTERN real ters_ga[10];
EXTERN real ters_n[10] INIT(zero10);
EXTERN real ters_c[10] INIT(zero10);
EXTERN real ter_c2[10] INIT(zero10);
EXTERN real ters_d[10] INIT(zero10);
EXTERN real ter_d2[10] INIT(zero10);
EXTERN real ters_h[10] INIT(zero10);
#endif

#ifdef KEATING
EXTERN real keating_alpha[55] INIT(zero55);
EXTERN real keat_alpha[10][10];
EXTERN real keating_d[55] INIT(zero55);
EXTERN real keat_d[10][10];
EXTERN real keating_r_cut[55] INIT(zero55);
EXTERN real keat_r_cut[10][10];
EXTERN real keat_r2_cut[10][10];
EXTERN real keating_beta[550] INIT(zero550);
EXTERN real keat_beta[10][10][10];
#endif

#ifdef COVALENT  
EXTERN int neigh_len INIT(NEIGH_LEN_INIT); /* initial neighbor table length */
#endif 

#ifdef EPITAX
EXTERN int  epitax_rate[10] INIT(zero10);  /* creation rate of atoms */
EXTERN int  epitax_type[10] INIT(zero10);  /* type of atom to be created */
EXTERN real epitax_mass[10] INIT(zero10);  /* mass of atom to be created */
EXTERN real epitax_temp[10] INIT(zero10);  /* temp. of atom to be created */
EXTERN real epitax_cutoff INIT(0.0);  
EXTERN int epitax_maxsteps;        /* number of timesteps with atom creation */
EXTERN int epitax_startstep INIT(0);  /* timestep where cretion begins */
EXTERN real epitax_ctrl INIT(1.0);    /* parameter for change NVE -> NVT */ 
EXTERN real epitax_height;
EXTERN real epitax_level;
EXTERN int epitax_number;
EXTERN int nepitax INIT(0);
EXTERN int epitax_sub_n INIT(0);
EXTERN real epitax_poteng_min INIT(0);
EXTERN real epitax_speed INIT(1.0);
#endif

#ifdef EWALD
EXTERN real     charge[10] INIT(zero10); /* Charge of atoms */
EXTERN real     ew_kappa;                /* Parameter kappa */
EXTERN ivektor  ew_kmax;                 /* Number of k-vectors in Ewald sum */
EXTERN int      ew_nmax;                 /* Number of image boxes */
EXTERN int      ew_nx;
EXTERN int      ew_ny;
EXTERN int      ew_nz;
EXTERN vektor   *ew_kvek;
EXTERN real     *ew_expk;
EXTERN int      ew_totk  INIT(0);
EXTERN real     *coskx;
EXTERN real     *sinkx;
EXTERN real     *cosky;
EXTERN real     *sinky;
EXTERN real     *coskz;
EXTERN real     *sinkz;
EXTERN real     *coskr;
EXTERN real     *sinkr;
EXTERN real     ew_vorf;
EXTERN real     ew_vorf1;
EXTERN real     ew_eps;
EXTERN real     twopi;
#endif

/* generate quasicrystal */
#ifdef QUASI
EXTERN int appr[3];                                /* approximant order */
EXTERN int k1min[6],k1max[6];                      /* grid boundaries  */
EXTERN real gx[6],gy[6],gz[6];                     /* grid vectors */
EXTERN real gam[6];                                /* grid offset */
EXTERN real tx[6],ty[6],tz[6];                     /* tiling vectors */

EXTERN vektor lmin,lmax;                           /* cell borders per pe */
EXTERN vektor gmin,gmax;                           /* global borders */

EXTERN real perm[3],perp[3];                       /* cell borders per pe */
EXTERN real iper[3];                               /* inverse cell size */
EXTERN int kf[6];                                  /* index borders, index */
EXTERN cell *input;                                /* save data */
#endif

/* simulate shock wave */
#ifdef SHOCK
EXTERN real shock_strip;                          /* width of shock strip */
EXTERN real shock_speed INIT(0.0);                /* velocity in shock strip */
EXTERN int  shock_mode INIT(2);                   /* type of shock */
#endif

EXTERN int ensemble INIT(ENS_EMPTY);    /* active ensemble type */
EXTERN int simulation INIT(1);          /* number of current simulation */
EXTERN int finished INIT(0);            /* last phase of simulation? */
EXTERN int loop INIT(0);                /* loop for online visualisation */
EXTERN void (*move_atoms)(void);        /* active integrator routine */

/* global parameters for UNIAX */
#ifdef UNIAX
EXTERN real uniax_r_cut  INIT(0.0);/* cutoff radius for uniaxial molecules */
EXTERN real uniax_r2_cut INIT(0.0);/* cutoff radius^2 for uniaxial molecules */
#endif

