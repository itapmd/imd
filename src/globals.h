
/******************************************************************************
*
* globals.h -- Global Variables for the imd Package
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
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
#else
#define EXTERN extern /* declare them extern otherwise */
#define INIT(data)
#define nullvektor
#define nullivektor
#define einsvektor
#define einsivektor
#endif

/* Where all the data lives */
EXTERN cell *cell_array INIT(NULL);    /* 3d Array of Cells */
EXTERN ivektor cell_dim;    /* Dimension of above (per cpu)*/
EXTERN ivektor global_cell_dim;    /* Dimension of above */

/* Global bookkeeping */

EXTERN int natoms INIT(0);       /* Total number of atoms */
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
EXTERN ivektor dist_dim INIT(nullivektor); /* Resolution of spatial distrib */
EXTERN str255 rundesc; /* Description */
EXTERN str255 progname; /* Name of current executable argv[0] */
EXTERN ivektor cellmin; /* Minimum index of local cells (1 with MPI, 0 otherwise) */
EXTERN ivektor cellmax; /* Maximum index of local cells  */

/* The simulation box and its inverse */
EXTERN vektor box_x INIT(nullvektor);
EXTERN vektor box_y INIT(nullvektor);
#ifndef TWOD
EXTERN vektor box_z INIT(nullvektor);
#endif
/* box parameters for generated structures */
EXTERN ivektor box_param INIT(nullivektor);
/* if initialization will be forgotten by user, then this could lead to division by zero, so we may have to check box_{x,y,z} later */
EXTERN vektor ibox_x;
EXTERN vektor ibox_y;
#ifndef TWOD
EXTERN vektor ibox_z;
#endif
EXTERN vektor tbox_x;
EXTERN vektor tbox_y;
#ifndef TWOD
EXTERN vektor tbox_z;
#endif

/* Filenames */
EXTERN str255 infilename INIT("\0");    /* Input File */
EXTERN str255 outfilename INIT("\0");   /* Output File */
EXTERN str255 potfilename INIT("\0");   /* Potential */
EXTERN char *paramfilename INIT(0L);    /* Parameter File */

/* construction of PN-dislocations */
EXTERN int pn INIT(0);  /* construct a Peierls-Nabarro-dislocation? (0==no) */
EXTERN real burgersv INIT(0);/* Burgers-vector */
EXTERN real width INIT(1);              /* width of Burgers-vector density */
EXTERN real upperplane INIT(0);         /* y/z (2/3D)-coordinate of glidep. */
EXTERN real lowerplane INIT(0);         /* y/z (2/3D)-coordinate of glidep. */

/* Parameters for displacement map */
EXTERN str255 reffilename INIT("\0");   /* Parameter File */

/* Parameters for pictures */
EXTERN vektor2d ecut_kin INIT(nullvektor2d);/* Kin. Energy interval for pictures */
EXTERN vektor2d ecut_pot INIT(nullvektor2d);/* Pot. Energy interval for pictures */   
EXTERN vektor pic_scale INIT(nullvektor);   /* Scale factor x/y for pictures     */
EXTERN vektor pic_ll  INIT(nullvektor);     /* lower left (front) corner */
EXTERN vektor pic_ur  INIT(nullvektor);     /* upper right (back) corner */
EXTERN ivektor   pic_res INIT(nullivektor); /* number of pixels in x/y dir.*/
EXTERN int       numpix  INIT(1);           /* number of pixels in x/y dir.*/
EXTERN int       pic_type INIT(0);          /* picture type 0/1/2 */
EXTERN real      *pic_at_radius INIT(NULL); /* atom radius for pictures */
#ifndef TWOD
EXTERN vektor3d view_dir INIT(nullvektor);  /* view direction */
EXTERN vektor3d view_pos INIT(nullvektor);  /* view position */
EXTERN int      projection INIT(0);         /* projection type 0/1 */
#endif

/* Monte Carlo stuff */
#ifdef MC
EXTERN real mc_beta   INIT(0.0);
EXTERN real mc_len    INIT(1.0);
EXTERN real mc_accept INIT(0.0);
EXTERN long mc_seed   INIT(0);
EXTERN int  mc_count  INIT(0);
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

#ifndef MONOLJ
/* Potential Table */
#ifdef STATIC_POT
/* Dynamic array inhibits vectorization on SX4 */
/* We have several copies of the pot table */
/* this reduces bank conflicts on vector computers */
#define MAXPOTLEN 4096
#define MAXATOMTYPES 2
EXTERN real potential[MAXATOMTYPES][MAXATOMTYPES][MAXPOTLEN];
#else
EXTERN real *potential;       /* Potential array */
#endif
EXTERN ivektor3d pot_dim INIT(nullivektor3d);       /* its dimensions */
EXTERN real r2_step;          /* delta between potential samples */
EXTERN real inv_r2_step;      /* inverse thereof */
EXTERN real r2_0;             /* minimum r^2 */
EXTERN real r2_end;           /* maximum r^2 */
#endif /* MONOLJ */
EXTERN real r2_cut;           /* cutoff^2 */   
EXTERN real cellsz INIT(0);   /* cell size */
EXTERN int  initsz INIT(0);   /* cell size */
/* MPI housekeeping */
EXTERN int myid INIT(0);                  /* Who am I? (0 if RISC) */
#ifdef MPI
EXTERN int num_cpus INIT(0);              /* How many cpus are there */
EXTERN int *cpu_ranks INIT(0);            /* Mapping of coords to ranks */
EXTERN ivektor cpu_dim INIT(nullivektor); /* Dimensions of CPU-Array */
EXTERN cell buf_one_atom;                 /* Buffer that holds one Atom */
EXTERN MPI_Comm cpugrid;
EXTERN ivektor my_coord;                  /* Cartesian coordinates of cpu */
EXTERN int parallel_output INIT(0);       /* Flag for parallel output */
EXTERN int parallel_input  INIT(1);       /* Flag for parallel input */

/* Send and Receive buffers */
EXTERN msgbuf send_buf_east;
EXTERN msgbuf send_buf_west;
EXTERN msgbuf send_buf_north;
EXTERN msgbuf send_buf_south;
EXTERN msgbuf send_buf_up;
EXTERN msgbuf send_buf_down;
EXTERN msgbuf recv_buf_east;
EXTERN msgbuf recv_buf_west;
EXTERN msgbuf recv_buf_north;
EXTERN msgbuf recv_buf_south;
EXTERN msgbuf recv_buf_up;
EXTERN msgbuf recv_buf_down;

/* Neighbours */
EXTERN int nbwest, nbeast, nbnorth, nbsouth, nbup, nbdown; /* Faces */
EXTERN int nbnw, nbws, nbse, nben,                         /* Edges */
           nbun, nbuw, nbus, nbue,
           nbdn, nbde, nbds, nbdw;
EXTERN int nbunw, nbuws, nbuse, nbuen,                     /* Corners */
           nbdne, nbdes, nbdsw, nbdwn;

/* MPI timing */
EXTERN double time_start;
EXTERN double time_last;
EXTERN double time_now;
EXTERN double time_io;
EXTERN double time_comm;
EXTERN double time_calc;
EXTERN double time_setup;
EXTERN double time_stop;
EXTERN double time_comm_force;
EXTERN double time_comm_ar;
EXTERN double time_calc_local;
EXTERN double time_calc_nonlocal;
#endif

/* Parameters for the various ensembles */

#ifdef AND
EXTERN int tmp_interval INIT(0);     /* Interval in which the thermostat */
                                     /* kicks (in timesteps) */
#endif

#if defined(NVT) || defined(NPT)
EXTERN real eta INIT(0.0);          /* Nose-Hoover heat bath variable */
EXTERN real isq_tau_eta INIT(0.0);  /* tau_eta: Nose-Hoover heat bath 'mass' */
                               /* isq_tau_eta : inverse of square of tau_eta */
#endif

#ifdef P_AXIAL
EXTERN vektor vir_vect INIT(nullvektor);
EXTERN vektor stress INIT(nullvektor);
#endif

EXTERN int    cells_too_small INIT(0);
EXTERN int    revise_cell_division INIT(0);

#ifdef NPT
EXTERN real   isq_tau_xi INIT(0.0); /* inverse of square of tau_xi */
EXTERN real   cell_size_tolerance INIT(0.05);
EXTERN vektor pressure_ext INIT(nullvektor);
EXTERN vektor pressure_end INIT(nullvektor);
EXTERN int    use_curr_pressure INIT(0);  /* which starting pressure to use */
EXTERN vektor xi INIT(nullvektor), xi_old INIT(nullvektor);
EXTERN vektor box_size INIT(einsvektor); 
EXTERN vektor actual_shrink INIT(einsvektor), limit_shrink INIT(einsvektor);
EXTERN vektor limit_growth INIT(einsvektor);
#endif

#if defined(AND) || defined(NVT) || defined(NPT)
EXTERN real end_temp INIT(0.0);        /* Temperature and at of simulation */
#endif

#if defined HOM
EXTERN int    hom_interval INIT(0);    /* period of homshear steps */
EXTERN real   shear_max INIT(0.0);     /* max shear in y-direction */
EXTERN int    exp_interval INIT(0);    /* period of expansion steps */
EXTERN real   expansion INIT(1.0);     /* max expansion in y-direction */
#endif

#if defined(FRAC) || defined(PULL)
EXTERN vektor stadion INIT(nullvektor);      /* Damping stadion */
EXTERN real   gamma_bar INIT(0.0);    /* Damping prefactor */
EXTERN real   gamma_cut INIT(0.0);    /* Damping cutoff */
#endif

#if defined(FRAC) || defined(PULL) || defined(SHOCK)
EXTERN int    dnoshsteps INIT(0);     /* counting steps between 2 shears */
EXTERN real   strip INIT(0.0);        /* Strip width */    
EXTERN real   ekin_threshold INIT(1.0e+20); /* threshold for ekin */    
EXTERN int    annealsteps INIT(0);    /* number of annealing steps */    
EXTERN int    maxdnoshsteps INIT(0);  /* max. steps between 2 shear steps */  
EXTERN int initial_shift INIT(0);  /* flag whether the sample is shifted */
EXTERN vektor ins INIT(nullvektor);/* initial shift */
EXTERN int shear_steps INIT(0);   /* number of shear_steps */
#endif

/*  #ifdef FRAC */
EXTERN real kcrit INIT(0.0);          /* Stress Intensity Factor */
EXTERN real mue INIT(0.0);            /* Youngs Modulus */
EXTERN real kel INIT(0.0);            /* Shear Modulus */
EXTERN vektor2d tip INIT(nullvektor2d);          /* Location of crack Tip */
/* #endif */

#ifdef PULL
EXTERN vektor   delta;      /* atoms in strip move by this amount */
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

#ifdef ORDPAR
#define nullvektor4d { 0.0, 0.0, 0.0, 0.0 }
EXTERN real op_weight[2][2] INIT(nullvektor4d);
EXTERN real op_r2_cut[2][2] INIT(nullvektor4d);
#endif

/* Global data for MSQD & correlation */
#if defined(MSQD) || defined(CORRELATE)
EXTERN real *msqd INIT(NULL);        /* (local) array of mean square disp. */
EXTERN real *msqd_global INIT(NULL); /* global array of mean square disp. */
EXTERN integer correl_int INIT(0);   /* repeat interval for correlation */
EXTERN integer correl_start INIT(0); /* start time for correlation */
EXTERN integer correl_end INIT(0);   /* end time for correlation */
EXTERN integer correl_ts INIT(0); /* sampling time interval for correlation */
EXTERN integer ncorr_rmax INIT(0); /* dimension of histogram in r domain */
EXTERN integer ncorr_tmax INIT(1); /* dimension of histogram in t domain */
EXTERN real GS_rcut INIT(0); /* cutoff radius for correlation data writes */
#endif

/* Global data for correlation */
#ifdef CORRELATE
EXTERN real inv_dr;                /* inverse of step size delta r */
EXTERN integer ***GS INIT(NULL);    /* histogram array for self correlation */
EXTERN integer correl_omode INIT(1); /* output mode for histogram */
#endif

#if defined(TRANSPORT)
EXTERN real dTemp_start   INIT(0.0);   /* deviation of starting heat bath 
                                          temperature from mean temp. */
EXTERN real dTemp_end 	  INIT(0.0);   /* deviation of final heat bath 
                                          temperature from mean temp. */ 
EXTERN real tran_Tleft 	  INIT(0.0);   /* temperature of the left=hot wall */
EXTERN real tran_Tright   INIT(0.0);   /* temperature of the right=cold  wall*/
EXTERN integer tran_interval INIT(0);  /* Intervalle der Temperaturaufz.*/
EXTERN integer tran_nlayers  INIT(0);  /* number of layers*/
EXTERN real heat_cond     INIT(0.0);   /* heat conductivity */
#endif

#ifdef STRESS_TENS
EXTERN integer press_interval INIT(0);  /* Intervalle der Aufzeichnung */ 
                                          /* des Drucktensors */
EXTERN ivektor press_nlayers;     /* Zahl der Schichten */
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
#endif

/* EAM stuff */
#ifdef EAM
/* neighborlist for EAM forces in 2D array */
EXTERN integer *eam_ij INIT(NULL);  	/* neighbor array */
EXTERN integer eam_len INIT(50);  	/* max neighbors */
EXTERN real *eam_rho   INIT(NULL);	/* cohesive function density (Finnis/Sinclair) */
EXTERN real *eam_dij_x INIT(NULL);	/* distance in x direction */
EXTERN real *eam_dij_y INIT(NULL);	/* distance in y direction */
EXTERN real *eam_dij_z INIT(NULL);	/* distance in z direction */
EXTERN real eam_r_cut  INIT(0.0);	/* EAM cutoff radius */
EXTERN real eam_r_0    INIT(1.0);	/* EAM minimum distance */
EXTERN real eam_r2_cut;  		/* EAM cutoff radius ^2 */
EXTERN real eam_A      INIT(0.0);	/* EAM cohesive function constant A */
#endif

/* TTBP */
#ifdef TTBP
/* neighborlist for TTBP */
EXTERN integer *ttbp_ij     INIT(NULL); /* neighbor array */
EXTERN integer ttbp_len     INIT(50);	/* max neighbors */
EXTERN real *ttbp_j 	    INIT(NULL);	/* position array */
EXTERN real *ttbp_force     INIT(NULL); /* force array */
EXTERN str255 ttbp_potfilename INIT("\0");   /* TTBP Potential */
EXTERN real *ttbp_potential;            /* TTBP Potential array */
EXTERN real ttbp_r2_cut[10][10];        /* cutoff^2;  less than 10 atom types! */   
EXTERN real ttbp_constant[10];		/* constants; less than 10 atom types! */ 
EXTERN real ttbp_theta[10];		/* theta;     less than 10 atom types! */ 
EXTERN real ttbp_r2_0;             	/* TTBP minimum r^2 */
EXTERN real ttbp_r2_end;           	/* TTBP maximum r^2 */
EXTERN real ttbp_r2_step;          	/* delta between potential samples */
EXTERN real ttbp_inv_r2_step;      	/* inverse thereof */
EXTERN ivektor3d ttbp_pot_dim INIT(nullivektor3d); /* pot dimensions */
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
EXTERN real shock_speed;                          /* velocity in shock strip */
EXTERN real shock_elong;                          /* atom elongation */
EXTERN int shock_mode;                           /* type of shock */
#endif

EXTERN int ensemble INIT(ENS_EMPTY);    /* active ensemble type */
EXTERN int simulation INIT(1);          /* number of current simulation */
EXTERN int finished INIT(0);            /* last phase of simulation? */
EXTERN void (*move_atoms)(void);        /* active integrator routine */

