
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

#else
#define nullvektor  { 0.0, 0.0, 0.0 }
#define nullivektor { 0, 0, 0 }
#define einsvektor  { 1.0, 1.0, 1.0 }
#define einsivektor { 1, 1, 1 }
#endif

#define nullvektor2d { 0.0, 0.0 }
#define nullivektor2d { 0, 0 }
#define nullvektor3d { 0.0, 0.0, 0.0 }
#define nullivektor3d { 0, 0, 0 }

#define nullbuffer  {NULL, 0, 0 }

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
EXTERN ivektor pbc_dirs INIT(einsivektor); /* directions with pbc */

EXTERN int vtypes INIT(0);   /* vtypes = ntypes + ntypes*(types of forces) */
EXTERN vektor *restrictions INIT(NULL);  /* directions the atom is allowed to move in */

#ifdef FBC    /* FBC uses the virtual atom types, */
EXTERN vektor *fbc_forces;       /* each vtype has its force */
EXTERN vektor *fbc_beginforces;  /* begin & endvalues for linear interpolation */ 
#ifdef MIK                       /* or ... */
EXTERN vektor *fbc_dforces;      /* each vtype has its force increment */
EXTERN real   fbc_ekin_threshold INIT(0.0);/* threshold for ekin */ 
EXTERN int    fbc_waitsteps INIT(1); /* between increase of forces */
EXTERN int    fbc_annealsteps INIT(1); /* times before 1. + df */
#else
EXTERN vektor *fbc_endforces;   
#endif

#endif

/* Global bookkeeping */

EXTERN int natoms  INIT(0);      /* Total number of atoms */
EXTERN int nactive INIT(0);      /* number of transl. degrees of freedom */
EXTERN int nactive_rot INIT(0);  /* number of rot. degrees of freedom */
EXTERN int ntypes INIT(0);       /* Total number of different atom types */
EXTERN int *num_sort INIT(NULL); /* number of atoms for each type */
EXTERN int steps_max INIT(0);    /* Maximum number of MD steps */
EXTERN int steps_min INIT(0);    /* starting step nr for current phase */
EXTERN int restart INIT(0);      /* file number for restart */
EXTERN int rep_interval INIT(0); /* Period of checkpoints ==0 for no checkpoints */
EXTERN int eng_interval INIT(0); /* Period of data output ==0 for no energy data */
EXTERN int dis_interval INIT(0); /* Period of spatial eng. distrib output ==0 for no data */
EXTERN int pic_interval INIT(0); /* Period of data output ==0 for no energy data */
EXTERN int onl_interval INIT(0); /* Period of online visualization */
EXTERN int dist_binary_io INIT(0); /* Flag for binary I/O */
EXTERN ivektor dist_dim INIT(einsivektor); /* Resolution of spatial distrib */
EXTERN str255 rundesc; /* Description */
EXTERN str255 progname; /* Name of current executable argv[0] */
EXTERN ivektor cellmin; /* Minimum index of local cells (1 with MPI, 0 otherwise) */
EXTERN ivektor cellmax; /* Maximum index of local cells  */

#ifdef EFILTER
EXTERN int efrep_interval INIT(0); /* Period of checkpoints ==0 for no checkpoints */
EXTERN real  lower_e_pot INIT(0.0); /* lower end of energy window */
EXTERN real  upper_e_pot INIT(0.0); /* upper end of energy window */
#endif

/* box parameters for generated structures */
EXTERN ivektor box_param INIT(nullivektor);

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
EXTERN vektor height INIT(nullvektor);

/* Filenames */
EXTERN char outbuf[OUTPUT_BUF_SIZE] INIT("\0");  /* output buffer */
EXTERN str255 infilename INIT("\0");    /* Input File */
EXTERN str255 outfilename INIT("\0");   /* Output File */
EXTERN str255 potfilename INIT("\0");   /* Potential */
EXTERN char *paramfilename INIT(0L);    /* Parameter File */

/* Parameters for displacement map */
EXTERN str255 reffilename INIT("\0");   /* Parameter File */

/* Parameters for pictures */
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
EXTERN real pressure INIT(0.0);
EXTERN real volume INIT(0.0);
EXTERN real virial INIT(0.0);
EXTERN real temperature INIT(0.0);
EXTERN int  use_curr_temp INIT(0);  /* which starting temp to use (flag) */
EXTERN int  do_maxwell INIT(0.0);
EXTERN long seed INIT(0);           /* seed for random number generator */

/* scalar product of global force vector f=(f1.x, f1.y,...,fn.z) */
EXTERN real fnorm INIT(0.0);  
/* scalar product of global force and momentum vectors */ 
EXTERN real PxF INIT(0.0);

/* Potential Table */
EXTERN pot_table_t pair_pot;         /* potential data structure */
EXTERN real monolj_r2_cut INIT(0.0); /* cutoff^2 for MONOLJ */   
EXTERN real monolj_shift INIT(0.0);  /* shift of monolj potential */   
EXTERN real cellsz INIT(0);          /* minimal cell diameter */
EXTERN int  initsz INIT(10);         /* initial number of atoms in cell */
EXTERN int  incrsz INIT(10);         /* increment of number of atoms in cell */


/* MPI housekeeping */
EXTERN int myid INIT(0);                  /* Who am I? (0 if RISC) */
EXTERN int parallel_output INIT(0);       /* Flag for parallel output */
EXTERN int parallel_input  INIT(1);       /* Flag for parallel input */
#ifdef MPI
EXTERN int binc INIT(0);                  /* buffer size per atom */
EXTERN int num_cpus INIT(0);              /* How many cpus are there */
EXTERN int *cpu_ranks INIT(0);            /* Mapping of coords to ranks */
EXTERN ivektor cpu_dim INIT(nullivektor); /* Dimensions of CPU-Array */
EXTERN cell buf_one_atom;                 /* Buffer that holds one Atom */
EXTERN MPI_Comm cpugrid;
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
EXTERN imd_timer time_force_comm;
EXTERN imd_timer time_force_calc;

/* Parameters for the various ensembles */

#ifdef AND
EXTERN int tmp_interval INIT(0);     /* Interval in which the thermostat */
                                     /* kicks (in timesteps) */
#endif

#if defined(NVT) || defined(NPT) || defined(STM)
EXTERN real eta INIT(0.0);          /* Nose-Hoover heat bath variable */
EXTERN real inv_tau_eta INIT(0.0);  /* tau_eta: Nose-Hoover heat bath 'mass' */
                                    /* inv_tau_eta : inverse of tau_eta */
#ifdef UNIAX
/* Nose-Hoover heat bath variable for rotational motion */
EXTERN real eta_rot INIT(0.0);      
/* tau_eta_rot: Nose-Hoover heat bath 'mass' for rotational motion */
/* inv_tau_eta_rot : inverse of tau_eta_rot */
EXTERN real inv_tau_eta_rot INIT(0.0);  
#endif
#endif

/* diagonal of virial tensor */
EXTERN real vir_x INIT(0.0), vir_y INIT(0.0), vir_z INIT(0.0);
/* diagonal of stress tensor */
EXTERN real stress_x INIT(0.0), stress_y INIT(0.0), stress_z INIT(0.0);

EXTERN int    revise_cell_division INIT(0);

#ifdef NPT
EXTERN real   inv_tau_xi INIT(0.0); /* inverse of tau_xi */
EXTERN real   cell_size_tolerance INIT(0.05);
EXTERN vektor pressure_ext INIT(nullvektor);
EXTERN vektor pressure_end INIT(nullvektor);
EXTERN int    use_curr_pressure INIT(0);  /* which starting pressure to use */
EXTERN vektor xi INIT(nullvektor), xi_old INIT(nullvektor);
EXTERN vektor box_size INIT(einsvektor); 
EXTERN vektor actual_shrink INIT(einsvektor), limit_shrink INIT(einsvektor);
EXTERN vektor limit_growth INIT(einsvektor);
#endif

#if defined(AND) || defined(NVT) || defined(NPT) || defined(STM)
EXTERN real end_temp INIT(0.0);        /* Temperature and at of simulation */
#endif

#if defined(FRAC) || defined(PULL)
EXTERN vektor stadion INIT(nullvektor);      /* Damping stadion */
EXTERN real   gamma_bar INIT(0.0);    /* Damping prefactor */
EXTERN real   gamma_cut INIT(0.0);    /* Damping cutoff */
#endif

#ifdef HOMDEF
EXTERN int    exp_interval INIT(0);       /* period of expansion steps */
EXTERN vektor expansion INIT(einsvektor); /* expansion factors in x/y/z-dir */
EXTERN int    hom_interval INIT(0);       /* period of homshear steps */
EXTERN vektor2d shear_factor INIT(nullvektor2d);/* shear factor in x,y-direction */
#endif

#ifdef SLLOD
EXTERN real   epsilon  INIT(0.0);         /* shear factor in x-direction */
#endif

#if defined(FRAC) || defined(STM)
EXTERN vektor stadium INIT(nullvektor); /* Damping stadium */
EXTERN vektor center  INIT(nullvektor); /* center of stadium */
EXTERN real   gamma_bar INIT(0.0);      /* Damping prefactor */
EXTERN real   gamma_cut INIT(0.0);      /* Damping cutoff */
#endif

#if defined(FRAC) || defined(DEFORM)
EXTERN int    deform_int INIT(0);      /* counting steps between 2 shears */
EXTERN real   strip_width INIT(0.0);   /* Strip width */    
EXTERN real   ekin_threshold INIT(0.0);/* threshold for ekin */    
EXTERN int    annealsteps INIT(0);     /* number of annealing steps */    
EXTERN int    max_deform_int INIT(0);  /* max. steps between 2 shear steps */  
#endif

#ifdef FRAC
EXTERN real kcrit INIT(0.0);          /* Stress Intensity Factor */
EXTERN real mue INIT(0.0);            /* Youngs Modulus */
EXTERN real kel INIT(0.0);            /* Shear Modulus */
EXTERN vektor2d tip INIT(nullvektor2d); /* Location of crack Tip */
#endif

#ifdef DEFORM
EXTERN vektor *deform_shift;       /* shift for each vtype */
#endif

#ifdef DISLOC
EXTERN int  dem_interval INIT(0);     /* Period of dem output ==0 */
EXTERN int  dsp_interval INIT(0);     /* Period of dsp output ==0 */
EXTERN int  up_ort_ref INIT(0);       /* time to update ort_ref ? */
EXTERN real min_dpot INIT(1.0);       /* difference for dem */
EXTERN int  dpotsorte INIT(0);        /* type to compute dem map */
EXTERN real ddelta INIT(1.0);         /* distance for ddm */
EXTERN int  reset_Epot_step INIT(0);  /* step at which Epot_ref is computed */
EXTERN int  calc_Epot_ref INIT(0);    /* flag whether to compute Epot_ref */
EXTERN int  Epot_diff INIT(1);        /* flag whether to write Epot_diff */
#endif

#ifdef AVPOS
EXTERN int avpos_int INIT(0);         /* Period of avp output ==0 */
EXTERN int avpos_res INIT(0);         /* Period of coordinate addition */
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
EXTERN int press_interval INIT(0);          /* Intervalle der Aufzeichnung */ 
                                            /* des Drucktensors */
EXTERN ivektor press_dim INIT(nullivektor); /* pressure histogram dimension */
#endif

/* I/O via sockets */
#ifdef USE_SOCKETS
EXTERN int socket_int INIT(1); /* interval for reading socket */
EXTERN struct sigaction act;   /* sigaction structure for signal handling */
EXTERN struct sigaction oact;  /* sigaction structure for signal handling */
EXTERN int socket_id;          /* initial socket parameter */
EXTERN int soc;                /* socket parameter after data request */
EXTERN unsigned short baseport INIT(31913);  /* base port for VolIMD socket */
EXTERN char display_host[256] INIT("");      /* name of controlling machine */
EXTERN unsigned long varIP INIT(0);
EXTERN int  use_socket_window INIT(0);        /* flag for using a window to write*/
EXTERN vektor socketwin_ll  INIT(nullvektor);     /* lower left (front) corner */
EXTERN vektor socketwin_ur  INIT(nullvektor);     /* upper right (back) corner */
EXTERN int  socket_atoms INIT(0);        /* counter for atoms to send window */
#endif

/* EAM stuff */
#ifdef EAM
/* neighborlist for EAM forces in 2D array */
EXTERN real *eam_ij    INIT(NULL);  /* neighbor array */
EXTERN int  eam_len    INIT(50);    /* max neighbors */
EXTERN real *eam_rho   INIT(NULL);  /* cohesive function density */ 
				    /* (Finnis/Sinclair) */
EXTERN real *eam_dij_x INIT(NULL);  /* distance in x direction */
EXTERN real *eam_dij_y INIT(NULL);  /* distance in y direction */
EXTERN real *eam_dij_z INIT(NULL);  /* distance in z direction */
EXTERN real eam_r_cut  INIT(0.0);   /* EAM cutoff radius */
EXTERN real eam_r_0    INIT(1.0);   /* EAM minimum distance */
EXTERN real eam_r2_cut;             /* EAM cutoff radius ^2 */
EXTERN real eam_A      INIT(0.0);   /* EAM cohesive function constant A */
#endif

/* EAM2 stuff */
#ifdef EAM2
EXTERN str255 eam2_emb_E_filename INIT("\0"); /* filenames for the tabulated functions */
EXTERN str255 eam2_at_rho_filename INIT("\0");
EXTERN str255 eam2_core_pot_filename INIT("\0");
/* function tables */
EXTERN real *eam2_f_i;                 /* table for the Embedding Energy as function of (rho_h) */
EXTERN real *eam2_rho_at;              /* table for the electron density  as function of the distance */
EXTERN real *eam2_phi;                 /* table for the core-core 2body Potential */
/* layout of the tables: rho_at[atom_i][atom_j][mapping index]
   i.e.                         0        1     int (r-r_end)/rstep 
   the same for phi
   and  for f_i[atom_i][mapping index]
*/
EXTERN real *eam2_rho_begin; /* for the most general case: the tables have  */
EXTERN real *eam2_rho_end;   /* different stepsizes, etc. for each atomtype */
EXTERN real *eam2_rho_step;  /* layout of the fields:                       */
EXTERN real *eam2_r_begin;   /* rho_begin[atom_i][atom_j][rho], etc.        */
EXTERN real *eam2_r_end;
EXTERN real *eam2_r_step;
EXTERN real *eam2_phi_r_begin;
EXTERN real *eam2_phi_r_end;
EXTERN real *eam2_phi_r_step;
EXTERN int eam2_max_r_steps INIT(0);  /* info needed to access the tables */
EXTERN int eam2_max_rho_steps INIT(0);
EXTERN int eam2_max_phi_r_steps INIT(0);
#endif

#ifdef TTBP
EXTERN str255 ttbp_potfilename INIT("\0"); /* TTBP smoothing potential file */
EXTERN pot_table_t smooth_pot;    /* TTBP smoothing potential */
EXTERN real ttbp_constant[10];	  /* constants; less than 11 atom types! */
EXTERN real ttbp_sp[10];          /* constants; less than 11 atom types! */
#endif

#ifdef TERSOFF
EXTERN real ters_r_cut[55];  /* cutoff^2;  less than 11 atom types! */
EXTERN real ter_r_cut[10][10];
EXTERN real ter_r2_cut[10][10];
EXTERN real ters_r0[55];        
EXTERN real ter_r0[10][10];
EXTERN real ters_a[55];        /* Parameters for Tersoff potential  */
EXTERN real ter_a[10][10];
EXTERN real ters_b[55];
EXTERN real ter_b[10][10];
EXTERN real ters_la[55];
EXTERN real ter_la[10][10];
EXTERN real ters_mu[55];
EXTERN real ter_mu[10][10];
EXTERN real ters_chi[45];
EXTERN real ter_chi[10][10];
EXTERN real ters_om[45];
EXTERN real ter_om[10][10];
EXTERN real ters_ga[10];
EXTERN real ters_n[10];
EXTERN real ters_c[10];
EXTERN real ter_c2[10];
EXTERN real ters_d[10];
EXTERN real ter_d2[10];
EXTERN real ters_h[10];
#endif

/* for TTBP and TERSOFF */
#ifdef COVALENT  
EXTERN int neigh_len INIT(50);     /* max neighbors */
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
EXTERN void (*move_atoms)(void);        /* active integrator routine */

/* global parameters for UNIAX */
#ifdef UNIAX
EXTERN real uniax_r_cut  INIT(0.0);/* cutoff radius for uniaxial molecules */
EXTERN real uniax_r2_cut INIT(0.0);/* cutoff radius^2 for uniaxial molecules */
#endif











