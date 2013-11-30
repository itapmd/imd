/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2012 Institute for Theoretical and Applied Physics,
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
#define xivektor    { 1, 0 }
#define nullsymtensor { 0.0, 0.0, 0.0 }
#else
#define nullvektor  { 0.0, 0.0, 0.0 }
#define nullivektor { 0, 0, 0 }
#define einsvektor  { 1.0, 1.0, 1.0 }
#define einsivektor { 1, 1, 1 }
#define parteinsivektor { 1, 1, 0 }
#define xivektor    { 1, 0, 0 }
#define nullsymtensor { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
#endif

#define nullvektor2d { 0.0, 0.0 }
#define nullivektor2d { 0, 0 }
#define nullvektor3d { 0.0, 0.0, 0.0 }
#define nullivektor3d { 0, 0, 0 }
#define nullvektor4d { 0.0, 0.0, 0.0, 0.0 }

#define nullbuffer  {NULL, 0, 0 }

#define zero10       {0,0,0,0,0,0,0,0,0,0}

#define zero45       {0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0}

#define one45        {1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1,\
 1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1}

#define zero55       {0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0}
#define zero100       {0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,\
}

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

#define minusone50            {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,\
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,\
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,-1,-1}

#else
#define EXTERN extern /* declare them extern otherwise */
#define INIT(data)
#endif

/* Where all the data lives */
EXTERN minicell *cell_array INIT(NULL);  /* array of cells */
#ifdef VEC
EXTERN cell atoms;                       /* big cell in vector mode */
EXTERN int  atoms_per_cpu INIT(150000);  /* estimated size of big cell */
#endif
EXTERN integer *cells INIT(NULL);        /* list if inner cell indices */
#ifdef TWOD
EXTERN pair *pairs[9];                   /* arrays of cell pairs */
EXTERN int  npairs[9], npairs2[9];       /* number of cell pairs */
#else
EXTERN pair *pairs[27];                  /* arrays of cell pairs */
EXTERN int  npairs[27], npairs2[27];     /* number of cell pairs */
#if defined(VEC) || defined(NBLIST)
EXTERN cell_nbrs_t *cnbrs INIT(NULL);    /* neighbors of each cell */
#endif
#endif
EXTERN int ncells, nallcells INIT(0);    /* number of cells */
EXTERN int ncells2;                      /* cells on lower bondary (for nbl) */
EXTERN int nlists;                       /* number of cell pair lists */
EXTERN ivektor cell_dim;                 /* dimension of cell array (per cpu)*/
EXTERN ivektor global_cell_dim;          /* dimension of cell array */

/* Boundary Conditions */
#ifdef EPITAX
EXTERN ivektor pbc_dirs INIT(parteinsivektor);
#else
EXTERN ivektor pbc_dirs INIT(einsivektor); /* directions with pbc */
#endif

EXTERN vektor *restrictions INIT(NULL);  /* directions the atom is allowed to move in */

#ifdef RELAX
EXTERN real ekin_threshold       INIT(0.0); /* threshold for Ekin */
EXTERN real fnorm_threshold      INIT(0.0); /* threshold for force norm */
EXTERN real f_max_threshold      INIT(-1.0);/* threshold for force maximum */
EXTERN real delta_epot_threshold INIT(0.0); /* threshold for delta Epot */
EXTERN int  is_relaxed           INIT(0);   /* sample is relaxed? */
#endif

#ifdef FBC                                 /* FBC uses virtual atom types, */
EXTERN vektor *fbc_forces      INIT(NULL); /* each vtype has its force */
EXTERN vektor *fbc_beginforces INIT(NULL); /* begin values for interpolation */
EXTERN vektor *fbc_endforces   INIT(NULL); /* end values for interpolation */
EXTERN vektor *fbc_df          INIT(NULL); /* computed force increment */
#ifdef RELAX
EXTERN vektor *fbc_dforces     INIT(NULL); /* force increment */
EXTERN int    max_fbc_int      INIT(1);    /* max int for force increment */
EXTERN int    fbc_int          INIT(0);    /* time since last FBC increase */
#endif
EXTERN int    have_fbc_incr    INIT(0);    /* fbc increment given */
EXTERN int    do_fbc_incr      INIT(0);    /* whether to apply FBC increment */
#endif /* FBC */

#ifdef BEND
EXTERN vektor *fbc_bforces      INIT(NULL); /* each vtype has its force */
EXTERN vektor *fbc_beginbforces INIT(NULL); /* begin values for interpolation */
EXTERN vektor *fbc_endbforces   INIT(NULL); /* end values for interpolation */
EXTERN vektor *fbc_bdf          INIT(NULL); /* computed force increment */
#ifdef RELAX
EXTERN vektor *fbc_bdforces     INIT(NULL); /* force increment */
EXTERN int    max_bfbc_int      INIT(1);    /* max int for force increment */
EXTERN int    bfbc_int          INIT(0);    /* time since last FBC increase */
#endif
EXTERN int    have_bfbc_incr    INIT(0);    /* fbc increment given */
EXTERN int    do_bfbc_incr      INIT(0);    /* whether to apply FBC increment */


EXTERN int    bend_nmoments    INIT(0);    /* number of bending moments,
                                              one moment for each vtype, e.g. in alloys  */
EXTERN int    bend_natomsvtype_origin[6];         /* number of atoms for each vtype which is
                                                     part of bending, to calc center of gravity */
EXTERN int    bend_natomsvtype_force[6];         /* number of atoms for each vtype which is
                                                     part of bending, to calc center of gravity */
EXTERN int    bend_vtype_of_origin[6];     /* match origin of bendmoment to a vtype      */
EXTERN int    bend_vtype_of_force[6];      /* match atoms with added force to a vtyp     */
EXTERN vektor *bend_axis       INIT(NULL); /* axis of bending moment*/
EXTERN vektor *bend_origin     INIT(NULL); /* origin of each moment */
EXTERN vektor *bend_cog        INIT(NULL); /* center of grav of atoms w. added force     */
EXTERN vektor *bend_vec        INIT(NULL); /* bend_cog - bend_origin    */
EXTERN vektor *bend_forces      INIT(NULL);/* each vtype has its force according to fbc_forces
                                              and the direction according to the bending moment*/
#endif

#ifdef ZAPP
EXTERN real   zapp_threshold INIT(0);
EXTERN vektor nactive_vect   INIT(nullvektor);
EXTERN vektor total_impuls   INIT(nullvektor);
#endif

#ifdef FLAGEDATOMS
EXTERN int flagedatomstype INIT(666);
#endif

/* Global bookkeeping */
EXTERN time_t tstart, tend;
EXTERN real maxwalltime INIT(0);  /* maximal allowed walltime */
EXTERN int watch_int INIT(0);     /* interval for checking write file */
EXTERN int stop_int INIT(0);      /* interval for checking stop file */
EXTERN int is_big_endian;         /* 1 if big endian, 0 if little endian */
EXTERN long natoms  INIT(0);      /* Total number of atoms */
EXTERN long nactive INIT(0);      /* number of transl. degrees of freedom */
EXTERN long nactive_rot INIT(0);  /* number of rot. degrees of freedom */
EXTERN int ntypes INIT(0);        /* number of real atom types */
EXTERN int vtypes INIT(0);        /* number of virtual atom types */
EXTERN int ntypepairs INIT(0);    /* Total number of different types pairs */
EXTERN int ntypetriples INIT(0);  /* Total number of different types triples
				     symmetric in the last two indices */
EXTERN int nvalues INIT(0);       /* either ntypes or ntypepairs */
EXTERN long *num_sort INIT(NULL); /* number of atoms for each real type */
EXTERN long *num_vsort INIT(NULL);/* number of atoms for each virtual type */
EXTERN int steps INIT(0);        /* number of current MD step */
EXTERN int steps_max INIT(0);    /* Maximum number of MD steps */
EXTERN int steps_min INIT(0);    /* starting step nr for current phase */
EXTERN int imdrestart INIT(0);   /* file number for restart */
EXTERN int checkpt_int INIT(0);  /* Period of checkpoints */
EXTERN int eng_int INIT(0);      /* Period of data output */
EXTERN int force_int INIT(0);    /* Period of force file writing */
EXTERN str255 progname; /* Name of current executable argv[0] */
EXTERN char *version_str;
EXTERN ivektor cellmin; /* Minimum index of local cells (1 with BUFCELLS,
                                                         0 otherwise) */
EXTERN ivektor cellmax; /* Maximum index of local cells  */
EXTERN int use_header INIT(1);   /* shall a header be written */
EXTERN int hyper_threads INIT(1); /* number of hyperthreads per CPU */

/* controlling distribution output */
EXTERN int dist_Ekin_flag        INIT(0); /* write Ekin dists? */
EXTERN int dist_Epot_flag        INIT(0); /* write Epot dists? */
EXTERN int dist_Ekin_long_flag   INIT(0); /* write Ekin_long dists? */
EXTERN int dist_Ekin_trans_flag  INIT(0); /* write Ekin_trans dists? */
EXTERN int dist_Ekin_comp_flag   INIT(0); /* write Ekin_comp dists? */
EXTERN int dist_shock_shear_flag INIT(0); /* write shock shear dists? */
EXTERN int dist_shear_aniso_flag INIT(0); /* write shear aniso dists? */
EXTERN int dist_press_flag       INIT(0); /* write press dists? */
EXTERN int dist_pressoff_flag    INIT(0); /* write press off diag dists? */
EXTERN int dist_presstens_flag   INIT(0); /* write presstens dists? */
EXTERN int dist_dens_flag   INIT(0); /* write density dists? */
EXTERN int dist_vxavg_flag  INIT(0); /* write average sample velocity dists? */
EXTERN int dist_int              INIT(0); /* Period of distribution writes */
EXTERN int dist_chunk_size       INIT(2*1024*1024); /* size of dist reduct. */
EXTERN ivektor dist_dim          INIT(einsivektor); /* resolution of dist */
EXTERN vektor  dist_ur           INIT(nullvektor);  /* lower left  corner */
EXTERN vektor  dist_ll           INIT(nullvektor);  /* upper right corner */

EXTERN int binary_output INIT(0);  /* write binary atoms data? */

#ifdef WRITEF
EXTERN int force_all INIT(0); /* write all forces, or only of virtual types */
#endif
#ifdef EFILTER
EXTERN int  ef_checkpt_int INIT(0);  /* Period of ef writes */
EXTERN real *lower_e_pot INIT(NULL); /* lower end of energy window */
EXTERN real *upper_e_pot INIT(NULL); /* upper end of energy window */
#endif
EXTERN int  nclones INIT(1);         /* number of clones */

#ifdef NNBR
EXTERN int  nb_checkpt_int INIT(0);    /* Period of nb writes */
EXTERN int  *lower_nb_cut  INIT(NULL); /* lower bound for neighbours  */
EXTERN int  *upper_nb_cut  INIT(NULL); /* upper bound for neighbours  */
EXTERN real *nb_r2_cut     INIT(NULL); /* cutoff for neighbor determination */
#endif

/* data for generated structures */
EXTERN ivektor box_param INIT(nullivektor);  /* box parameters */
EXTERN int  size_per_cpu INIT(0);            /* box_param is given per cpu */
EXTERN real box_unit INIT(1.0);              /* lattice parameter */
EXTERN real *masses INIT(NULL);              /* masses */
EXTERN int  *gtypes INIT(NULL);              /* types */

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
EXTERN char *outbuf INIT(NULL);         /* output buffer */
EXTERN int outbuf_size INIT(OUTPUT_BUF_SIZE * 1048576);
EXTERN int inbuf_size  INIT(INPUT_BUF_SIZE  * 1048576);
EXTERN str255 infilename INIT("\0");    /* Input File */
EXTERN str255 itrfilename INIT("\0");   /* initial itr-file */
EXTERN str255 outfilename INIT("\0");   /* Output File */
EXTERN str255 potfilename INIT("\0");   /* Potential */
EXTERN str255 paramfilename INIT("\0"); /* Parameter File */

/* pointers to files that are kept open */
EXTERN FILE  *eng_file INIT(NULL);      /* pointer to .eng file  */
EXTERN FILE *msqd_file INIT(NULL);      /* pointer to .msqd file */
EXTERN int flush_int INIT(50);          /* flush .eng and .msqd files
                                           every flush_int writes */
#ifdef EXTPOT
EXTERN FILE *ind_file INIT(NULL);       /* pointer to .ind file */
EXTERN str255 extpotfilename INIT("\0");   /* Potential */
#endif

#ifdef RELAX
EXTERN FILE *ssdef_file INIT(NULL);     /* pointer to .ssdef file */
EXTERN int max_sscount INIT(0);         /* max  quasistat. simulation steps */
#endif

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
EXTERN int  use_curr_temp INIT(0);   /* which starting temp to use (flag) */
EXTERN int  do_maxwell INIT(0);
EXTERN long seed INIT(0);            /* seed for random number generator */
EXTERN int  box_from_header INIT(0); /* read box from config file */
#ifdef NBLIST
EXTERN real nbl_margin INIT(0.4);    /* neighbor list margin */
EXTERN real nbl_size   INIT(1.1);    /* neighbor list size */
EXTERN int  nbl_count  INIT(0);      /* counting neighbor list rebuild */
EXTERN int  have_valid_nbl INIT(0);
EXTERN int  last_nbl_len   INIT(0);
#endif

/* quantities relevant for checking the relaxation process */
/* square of global force vector f=(f1.x, f1.y,...,fn.z) */
EXTERN real old_epot INIT(0.0);
EXTERN real xnorm INIT(0.0);
EXTERN real x_max2 INIT(0.0);
EXTERN real fnorm INIT(0.0);
EXTERN real pnorm INIT(0.0);
EXTERN real fnorm_old INIT(0.0);
EXTERN real f_max INIT(0.0);
EXTERN real f_max2 INIT(0.0);
/* scalar product of global force and momentum vectors */
EXTERN real PxF INIT(0.0);
EXTERN real last_PxF INIT(0.0);
/* Einstein frequency is similar as fnorm, but divided by the masses */
EXTERN real omega_E INIT(0.0);

/* Potential Table */
EXTERN pot_table_t pair_pot;         /* potential data structure */
#ifdef EXTPOT
EXTERN pot_table_t ext_pot;         /* potential data structure */
#endif
#ifdef MULTIPOT
EXTERN pot_table_t pair_pot_ar[N_POT_TAB]; /* array of potential tables */
#endif
#ifdef LINPOT
EXTERN lin_pot_table_t pair_pot_lin; /* potential data structure */
#endif
EXTERN real cellsz INIT(0);          /* minimal cell diameter */
EXTERN int  initsz INIT(10);         /* initial number of atoms in cell */
EXTERN int  incrsz INIT(10);         /* increment of number of atoms in cell */
EXTERN int  debug_potential INIT(0);   /* write out interpolated potential */
EXTERN int  debug_pot_res INIT(10000); /* resolution of the above */

/* MPI housekeeping */
EXTERN int myid INIT(0);                  /* Who am I? (0 if serial) */
EXTERN int num_cpus INIT(1);              /* How many cpus are there */
EXTERN int parallel_output INIT(0);       /* Flag for parallel output */
EXTERN int parallel_input  INIT(1);       /* Flag for parallel input */
EXTERN ivektor my_coord INIT(nullivektor);/* Cartesian coordinates of cpu */
EXTERN ivektor cpu_dim INIT(einsivektor); /* Dimensions of CPU-Array */
EXTERN int binc INIT(0);                  /* buffer size per atom */
EXTERN int *cpu_ranks  INIT(NULL);        /* Mapping of coords to ranks */
EXTERN int *io_grps    INIT(NULL);  /* mapping of ranks to IO groups */
EXTERN int my_inp_grp   INIT(0);    /* my input group */
EXTERN int my_inp_id    INIT(0);    /* input rank for my input group */
EXTERN int inp_grp_size INIT(1);    /* size of my input group */
EXTERN int n_inp_grps   INIT(1);    /* number of input groups */
EXTERN int my_out_grp   INIT(0);    /* my output group */
EXTERN int my_out_id    INIT(0);    /* output rank for my output group */
EXTERN int out_grp_size INIT(1);    /* size of my output group */
EXTERN int n_out_grps   INIT(1);    /* number of output groups */
#ifdef MPI
EXTERN MPI_Comm cpugrid;                  /* Cartesian MPI communicator */

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
EXTERN msgbuf dump_buf       INIT(nullbuffer);
EXTERN real   msgbuf_size    INIT(1.2);
EXTERN int    atom_size      INIT(0);

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
EXTERN imd_timer time_output;
EXTERN imd_timer time_fft;
EXTERN imd_timer time_fft_plan;
EXTERN imd_timer time_input;
EXTERN imd_timer time_integrate;
EXTERN imd_timer time_forces;

/* Parameters for the various ensembles */

#ifdef AND
EXTERN int tempintv INIT(0);        /* Interval in which the thermostat */
                                    /* kicks (in timesteps) */
#endif

#ifdef BER
EXTERN real tauber INIT(0);           /* Relaxation constant Berendsen Thermo */
                                     /* (in timesteps) */
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
#ifdef EXTPOT
EXTERN int  have_extpotfile INIT(0);
EXTERN int ep_key INIT(0);
EXTERN int ep_n INIT(0);
EXTERN int ep_nind INIT(1);
EXTERN int ep_int INIT(0), ep_max_int INIT(0);
EXTERN real ep_rcut INIT(0.0), ep_a INIT(0.0);
EXTERN real ep_fext[10];
EXTERN real ep_xmin[10], ep_xmax[10], ep_ymin[10], ep_ymax[10];
EXTERN vektor ep_pos[10];
EXTERN vektor ep_vel[10];
EXTERN vektor ep_dir[10];
EXTERN long nactive_vect[3];
EXTERN int ep_atomsincontact[10];
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

#ifdef GLOK
EXTERN real   glok_ekin_threshold INIT(100.0); /* threshold for ekin */
EXTERN int    glok_int       INIT(0);
EXTERN int    glok_start     INIT(0);
#endif
#ifdef MIX
EXTERN real mix              INIT(0.0);
EXTERN real glok_mix         INIT(0.0);
EXTERN real glok_mixdec      INIT(1.0);
EXTERN real mixforcescalefac  INIT(0.0);
#endif
#ifdef ADAPTGLOK
EXTERN real glok_incfac      INIT(1.02);
EXTERN real glok_decfac      INIT(0.5);
EXTERN real glok_maxtimestep INIT(0.0);
EXTERN real starttimestep    INIT(0.0);
EXTERN int  glok_minsteps    INIT(5); /* threshold for minsteps */
EXTERN int  nPxF             INIT(0);
EXTERN int  min_nPxF         INIT(0);
#endif

EXTERN real glok_fmaxcrit    INIT(10000);


#ifdef DEFORM
EXTERN int    max_deform_int INIT(0);   /* max. steps between 2 shear steps */
EXTERN int    deform_int     INIT(0);   /* curr. step between 2 shear steps */
EXTERN real   deform_size INIT(1.0);    /* scale factor for deformation */
EXTERN vektor *deform_shift;            /* shift for each vtype */
EXTERN vektor *deform_shear;            /* shear for each vtype */
EXTERN vektor *deform_base;             /* base point for shear deformation */
EXTERN int    *shear_def;               /* shear flag for each vtype */
#endif
#ifdef CYCLE
EXTERN real   lindef_freq INIT(0.0);      /* frequency for deformation */
#endif
#ifdef HOMDEF
EXTERN int    lindef_int INIT(0);         /* period of linear deform. steps */
EXTERN real   lindef_size INIT(1.0);      /* scale factor of deformation */
EXTERN vektor lindef_x INIT(nullvektor);  /* \               */
EXTERN vektor lindef_y INIT(nullvektor);  /*  |  linear      */
#ifndef TWOD                              /*   > deformation */
EXTERN vektor lindef_z INIT(nullvektor);  /*  |  matrix      */
#endif                                    /* /               */
EXTERN real shear_module INIT(1.0);       /* estimate of the shear module */
EXTERN real bulk_module  INIT(1.0);       /* estimate of the bulk module */
EXTERN int  relax_mode   INIT(-1);        /* pressure relaxation mode */
#endif
EXTERN ivektor relax_dirs INIT(einsivektor); /* directions in which to relax pressure */
EXTERN real relax_rate   INIT(0.0);       /* pressure relaxation rate */

#ifdef RIGID
EXTERN int  nsuperatoms INIT(0);          /* number of superatoms */
EXTERN int  *superatom INIT(NULL);   /* maps virtual types to superatoms */
EXTERN vektor *superrestrictions INIT(NULL); /* restricted rigidity */
EXTERN vektor *superforce INIT(NULL);     /* total force on superatoms */
EXTERN real *supermass INIT(NULL);        /* masses of superatoms */
EXTERN int  *num_ssort INIT(NULL);   /* number of atoms for each superatom */
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

#ifdef DAMP
EXTERN vektor   center  INIT(nullvektor); /* center of stadium */
EXTERN vektor  stadium  INIT(nullvektor); /* half axes of Damping stadium */
EXTERN vektor  stadium2 INIT(nullvektor); /* half axes where the max.
                                                damping factor  is reached  */
EXTERN real delta_finnis INIT(0.05);       /* time constant in finnis */
EXTERN real zeta_0       INIT(0.0);        /* prefactor for finnis*/

EXTERN real damptemp INIT(0.0);
EXTERN real tot_kin_energy_damp INIT(0.0); /* kin energy of damping layers  */
EXTERN int n_damp INIT(0);                 /* number of transl. degrees of  */
                                           /* freedom in the damping layers */
#endif

#ifdef FINNIS
EXTERN real delta_finnis INIT(10.0);       /* time constant in finnis */
EXTERN real zeta_0       INIT(0.0);        /* prefactor for finnis */
#endif

#ifdef FRAC
EXTERN real E_kin_damp INIT(0.0);         /* weighted !!  kin dampenergy  */
EXTERN real sum_f INIT(0.0);              /* Sum of stadium function */
#endif

#ifdef FTG
EXTERN real  Tleft   INIT(0.0);          /* Temperature of the left  wall*/
EXTERN real  Tright  INIT(0.0);          /* Temperature of the right wall*/
EXTERN int   nslices INIT(0);            /* Number of slices*/
EXTERN int   nslices_Left  INIT(1);      /* Number of slices with Tleft */
EXTERN int   nslices_Right INIT(1);      /* Number of slices with Tright*/
EXTERN int  *ninslice;                   /* Number of atoms in slice*/
EXTERN real *gamma_ftg;                  /* Damping prefactor for slices*/
EXTERN real *E_kin_ftg;                  /* kin energy of the slices */
#endif

#ifdef CG
/* Parameters used by CG */
EXTERN int    linmin_maxsteps  INIT(100);/* max number of linmin steps */
EXTERN real   linmin_tol  INIT(0.0002);  /* tolerance between 2 linmin steps */
EXTERN real   linmin_dmax  INIT(0.01);   /* max. search steps in linmin  */
EXTERN real   linmin_dmin  INIT(0.001);  /* min. search steps in linmin  */
EXTERN real   cg_glimit  INIT(100);
EXTERN real   cg_zeps  INIT(1e-10);
EXTERN int    cg_infolevel INIT(0);     /* cg_infolevel controls verbosity */
EXTERN int    cg_mode INIT(0);          /* CG mode */
EXTERN int    cg_fr INIT(0);            /* Fletcher-Reeves mode or not */
EXTERN int    cg_reset_int INIT(0);     /* interval between cg resetting */

/* Variables needed by CG */
EXTERN real   cg_poteng     INIT(0.0);      /* potential energy per atom */
EXTERN real   old_cg_poteng INIT(0.0);      /* old poteng value */
EXTERN real   gg            INIT(0.0);      /* see Num. Rec. p.320 */
EXTERN real   dgg           INIT(0.0);      /* see Num. Rec. p.320 */
EXTERN real   cg_gamma      INIT(0.0);      /* see Num. Rec. p.320 */
#endif

#ifdef ACG
EXTERN real   acg_alpha         INIT(0.005);  /* Kai Nordlunds adaptive CG */
EXTERN real   acg_init_alpha    INIT(0.005);  /* Kai Nordlunds adaptive CG */
EXTERN real   acg_incfac        INIT(1.05);   /* Kai Nordlunds adaptive CG */
EXTERN real   acg_decfac        INIT(0.5);    /* Kai Nordlunds adaptive CG */
#endif

#ifdef RELAX
EXTERN int sscount INIT(0);           /* snapshot counter */
#endif
EXTERN int nfc INIT(0);               /* counts force computations */

#ifdef CNA
EXTERN int cna_start INIT(0);
EXTERN int cna_end INIT(0);
EXTERN int cna_int INIT(0);
EXTERN real cna_rcut INIT(1.0);
EXTERN real cna_r2cut INIT(1.0);
EXTERN vektor cna_ll  INIT(nullvektor);     /* lower left (front) corner */
EXTERN vektor cna_ur  INIT(nullvektor);     /* upper right (back) corner */
EXTERN int cna INIT(0);
EXTERN int cna_pairs INIT(0);
EXTERN int cna_writev[8];
EXTERN int cna_writec INIT(0);
EXTERN int cna_write_n INIT(0);
EXTERN int cna_write_statistics INIT(0);
EXTERN char cna_crist INIT(0);
EXTERN int cna_cristv[4];
EXTERN int cna_crist_n INIT(0);
EXTERN int bondlist[MAX_BONDS][3];
EXTERN int type_list[MAX_TYPES];
EXTERN int type_count[MAX_TYPES];
EXTERN int type_sort[MAX_TYPES];
EXTERN int type_list_length INIT(0);
#endif

#ifdef ADA
EXTERN real ada_nbr_r2cut INIT(0.0);		/* Squared nearest neighbor cutoff radius */
EXTERN shortint ada_default_type INIT(127); /* Atom with matching ada-type are ignored in ada-files */
EXTERN int ada_write_int INIT(0);		/* Interval to write ada-files*/
EXTERN int ada_crystal_structure INIT(ADA_FCC_CONFIG);
EXTERN real ada_latticeConst INIT(0.);
#ifdef NYETENSOR
EXTERN vektor nye_rotationAxis_x  INIT(nullvektor); /* Crystal orientation in x-direction e.g. 1 -1 0*/
EXTERN vektor nye_rotationAxis_y  INIT(nullvektor); /* Crystal orientation in y-direction e.g. 0 0 1*/
EXTERN vektor nye_rotationAxis_z  INIT(nullvektor); /* Crystal orientation in z-direction e.g. 1 1 0 */
#endif
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
EXTERN int avpos_steps INIT(0);       /* Number of steps to average over before position writes*/
EXTERN int avpos_nwrites INIT(0);     /* Number of position writes performed */
EXTERN int avpos_npwrites INIT(0);    /* Number of pressure writes performed */
EXTERN int avpos_cnt INIT(0);         /* Number of positions added */
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
EXTERN real *op_weight INIT(NULL);
EXTERN real *op_r2_cut INIT(NULL);
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
EXTERN int  correl_refstep INIT(0);  /* reference step for correlation */
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

/* data for heat conductivity measurements */
#if defined(HC) || defined(NVX)
EXTERN int hc_start INIT(2000);        /* heat current starting time */
EXTERN int hc_int   INIT(0);           /* heat current writing interval */
#endif
#ifdef HC
EXTERN int hc_av_start INIT(1000);     /* energy average starting time  */
EXTERN vektor hc;                      /* heat current */
EXTERN FILE *hc_file   INIT(NULL);     /* heat current file */
EXTERN real heat_cond  INIT(0.0);      /* heat conductivity */
#endif
#ifdef NVX
EXTERN int  hc_nlayers  INIT(0);       /* number of layers */
EXTERN int  hc_count    INIT(0);       /* running index of temp. profile */
EXTERN real hc_heatcurr INIT(0.0);     /* induced heat current density */
#endif

#ifdef DEBUG
EXTERN ivektor force_celldim_divisor INIT(einsivektor); /* if you want cell dimensions to be
							   divisible by certain numbers, this is your
							   parameter. Useful for comparing
							   serial and parallel simulations if
							   cell dimensions make a difference */
#endif /* DEBUG*/

#ifdef TTM /* two temperature model */
  /* These will point to the calculation net in a nice 3D array fashion */
EXTERN  ttm_Element *** l1, *** l2, *** l3;
  /* These will be used to allocate and free the nets in bulk */
EXTERN  ttm_Element * lattice1, * lattice2;
EXTERN real fd_k INIT(1.0); /* electronic thermal conductivity */
EXTERN real fd_c INIT(0.0); /* electronic thermal capacity */
EXTERN real fd_gamma INIT(0.0); /* fd_c / T_e, proport. const. */
EXTERN real fd_g INIT(1.0);        /* electron-phonon coupling constant */
EXTERN int fd_n_timesteps INIT(1); /* how many FD steps to a MD timestep? */
EXTERN int fd_update_steps INIT(1);/* how often are FD cells updated
				      by averaging over atoms ? */
EXTERN int fd_min_atoms INIT(3);   /* minimum number of atoms needed in a
				      FD cell for it to be active */
EXTERN int ttm_int INIT(0); /* How many steps before ttm writeouts? */
EXTERN real ttm_eng INIT(0.0); /* Electronic heat energy per atom for .eng file*/
EXTERN str255 fd_one_d_str INIT("\0"); /* string for param file input */
EXTERN int fd_one_d	  INIT(0); /* FD lattice is one dimensional in direction 1,2,3 (x,y,z), else 3D */
EXTERN ivektor fd_ext	  INIT(einsivektor); /* how many MD cells in x,y,z per FE cell
                            (two of the values will be ignored if fd_one_d != 0)*/
EXTERN int n_md_cells     INIT(1); /* number of MD cells in one FD cell */
EXTERN int natoms_local INIT(0); /* number of atoms in this process */
EXTERN real init_t_el INIT(0.0); /* temperature to initialize electron system
				    to, for simple relaxation simulations.
				    Use lattice temp if ==0 */
EXTERN int fix_t_el INIT(0); /* fix electron temperature to init_t_el? */
EXTERN vektor fd_h	  INIT(nullvektor); /* lattice constants in coordinate directions of FE part */
EXTERN ivektor global_fd_dim INIT(nullivektor); /* global FD dimensions, w/o BC cells */
EXTERN ivektor local_fd_dim INIT(nullivektor); /* local FD dimensions incl. 2 ghost layers in every direction */
#ifdef MPI
/* MPI Datatypes for parts of the 3d Array.
 *  mpi_*_block are the ones to be used,
 * the others are intermediate and used to build other types */
EXTERN MPI_Datatype  mpi_element, mpi_element2,
	             mpi_zrow, mpi_zrow_block,
		     mpi_xplane, mpi_xplane_block,
		     mpi_yplane, mpi_yplane_block,
		     mpi_yrow, mpi_yrow_block,
		     mpi_zplane, mpi_zplane_block;
EXTERN MPI_Status stati[6];
EXTERN MPI_Request reque[6];
#endif /*MPI*/
EXTERN double E_new INIT(0.0); /* Energy of newly created FD cells */
EXTERN double E_new_local INIT(0.0);
#ifdef DEBUG /* for debugging (duh) */
EXTERN double E_el_ab INIT(0.0);  /* Energy taken from electrons */
EXTERN double E_el_ab_local INIT(0.0);
EXTERN double E_ph_auf INIT(0.0); /* Energy put into lattice */
EXTERN double E_ph_auf_local INIT(0.0);
#endif /*DEBUG*/
#endif /*TTM*/

#ifdef LASER /* Laser parameters */

EXTERN real laser_delta_temp INIT(0.0);   /* maximum Temperature added
				             (at surface of the sample)      */
EXTERN real laser_mu	  INIT(0.0);   /* absorption coefficient             */
EXTERN ivektor laser_dir  INIT(xivektor); /* direction of laser incidence    */
EXTERN real laser_offset INIT(0.0);     /* how much "space" between 0 and the surface of sample */
EXTERN real laser_sigma_e INIT(0.0);   /* area density of pulse energy (for rescaling method) */
EXTERN real laser_sigma_t INIT(0.5);   /* half pulse duration (sigma of gaussian pulse) (rescaling) */
EXTERN real laser_sigma_t_squared INIT(0.25); /* same, squared */
EXTERN real laser_t_0	  INIT(1.0);   /* time of maximum intensity of pulse (rescaling) */
EXTERN real laser_sigma_e1 INIT(0.0);   /* area density of second pulse energy (for rescaling method) */
EXTERN real laser_sigma_t1 INIT(0.5);   /* half pulse duration of the second pulse (sigma of gaussian pulse) (rescaling) */
EXTERN real laser_sigma_t1_squared INIT(0.25); /* same, squared */
EXTERN real laser_t_1	  INIT(1.0);   /* time of maximum intensity of the second pulse (rescaling) */
EXTERN real laser_p_peak  INIT(0.0);   /* Peak power density (calculated in imd.c from previous parameters)*/
EXTERN real laser_p_peak1  INIT(0.0);   /* Peak power density (calculated in imd.c from previous parameters)*/
EXTERN real laser_atom_vol INIT(16.6);  /* Volume per particle (inverse density) ATTENTION: THIS VALUE IS VOR ALUMINUM ONLY*/
EXTERN int  laser_rescale_mode INIT(1); /* Mode for laser velocity rescaling */

#ifdef LASERYZ
EXTERN real laser_sigma_w_y INIT(0.0); /* y-center of gaussian laser-pulse  */
EXTERN real laser_sigma_w_z INIT(0.0); /* z-center of gaussian laser-pulse  */
EXTERN real laser_sigma_w0 INIT(10.0); /* diameter of the laser beam at 1/e */
EXTERN ivektor laser_tem_mode INIT(nullivektor); /* Defines the laser TEM_xy mode */
#endif

EXTERN void (*do_laser_rescale)(void);  /* Function pointer for rescaling routine */
EXTERN double (*laser_intensity_profile)(double, double, double);

#endif

#ifdef PDECAY
EXTERN real xipdecay INIT (0.0);        /* damping parameter */
EXTERN real ramp_start INIT (0.9);      /* start of damping ramp in % of box_size */
EXTERN real ramp_end INIT (1.0);        /* end of damping ramp */
EXTERN real ramp_fraction INIT (0.2);   /* fraction of the sample (from right side) on which the damping ramp acts */
EXTERN int pdecay_mode INIT (1);        /* which form of the damping function is used */
#endif

#ifdef SM
EXTERN int  charge_update_steps INIT(0); /* number of steps between charge updates */
EXTERN int  sm_fixed_charges INIT(0);    /* if 1, keep charges fixed */
EXTERN real sm_chi_0[2]; /* Initial value of the electronegativity */
EXTERN real sm_Z[2];     /* Initial value of the effecitve core charge */
EXTERN real sm_J_0[2];   /* atomic hardness or self-Coulomb repulsion */
EXTERN real sm_zeta[2];
EXTERN str255 na_pot_filename INIT("\0");     /* nuclear attraction potential file   */
EXTERN str255 cr_pot_filename INIT("\0");     /* coulomb repulsive potential file   */
EXTERN str255 erfc_filename INIT("\0");     /* tabulated function erfc/r file   */
#endif


#ifdef STRESS_TENS
EXTERN sym_tensor tot_presstens INIT(nullsymtensor);/* global pressure tens. */
EXTERN sym_tensor presstens_ext INIT(nullsymtensor);  /* ext. pressure tens. */
EXTERN int press_int INIT(0);    /* interval for writing the pressure tensor */
EXTERN int do_press_calc INIT(0);   /* flag whether to do press calc */
#endif

/* I/O via sockets */
#ifdef SOCKET_IO
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
EXTERN real lj_aaa  [10][10];
#ifdef VEC
EXTERN real lj_epsilon_vec[55] INIT(zero55);
EXTERN real lj_sigma2_vec[55] INIT(zero55);
#endif
/* Lennard-Jones-Gauss */
EXTERN real ljg_eps_lin[55] INIT(zero55);
EXTERN real ljg_eps[10][10];
EXTERN real ljg_sig_lin[55] INIT(zero55);
EXTERN real ljg_sig[10][10];
EXTERN real ljg_r0_lin[55] INIT(zero55);
EXTERN real ljg_r0[10][10];
/* Morse */
EXTERN real morse_epsilon_lin[55] INIT(zero55);
EXTERN real morse_epsilon[10][10];
EXTERN real morse_sigma_lin[55] INIT(zero55);
EXTERN real morse_sigma[10][10];
EXTERN real morse_alpha_lin[55] INIT(zero55);
EXTERN real morse_alpha[10][10];
EXTERN real morse_shift[10][10];
EXTERN real morse_aaa  [10][10];
/* Buckingham */
EXTERN real buck_a_lin[55] INIT(zero55);
EXTERN real buck_a[10][10];
EXTERN real buck_c_lin[55] INIT(zero55);
EXTERN real buck_c[10][10];
EXTERN real buck_sigma_lin[55] INIT(zero55);
EXTERN real buck_sigma[10][10];
EXTERN real buck_shift[10][10];
EXTERN real buck_aaa  [10][10];
/* harmonic potential for shell model */
EXTERN real spring_const[45] INIT(zero45);
EXTERN real spring_cst[10][10];
#ifdef EWALD
EXTERN real ew_shift [10][10];
EXTERN real ew_fshift[10][10];
#endif
EXTERN int  fix_bks INIT(0);
#endif

#ifdef FEFL
EXTERN real spring_rate[10] INIT(zero10);
/* EXTERN real spring_rt[10]; */
EXTERN real lambda INIT(1.0);
EXTERN real tot_harm_energy INIT(0.0);
#endif

#ifdef EAM2
EXTERN pot_table_t embed_pot;                     /* embedding energy table  */
EXTERN pot_table_t rho_h_tab;                     /* electron transfer table */
EXTERN str255 eam2_emb_E_filename INIT("\0");     /* embedding energy file   */
EXTERN str255 eam2_at_rho_filename INIT("\0");    /* electron transfer file  */
#ifdef EEAM
EXTERN pot_table_t emod_pot;                      /* energy mod. term table  */
EXTERN str255 eeam_mod_E_filename INIT("\0");     /* energy mod. term file   */
#endif
#endif

#ifdef ADP
EXTERN pot_table_t adp_upot;      /* dipole     distortion potential      */
EXTERN pot_table_t adp_wpot;      /* quadrupole distortion potential      */
EXTERN str255      adp_upotfile;  /* dipole     distortion potential name */
EXTERN str255      adp_wpotfile;  /* quadrupole distortion potential name */
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
EXTERN real ttbp_constant2[8];	  /* constants for Vashishta potential, only for two atom types! */
EXTERN real ttbp_sp[10];          /* constants; less than 11 atom types! */
EXTERN real ttbp_cut[1];          /* cutoff for smoothing part of vashishta potential */
EXTERN int ttbp_vas INIT(0);      /* switch for vashishta potential */
EXTERN real B[2][2][2];           /* factor B_ijk for vashishta potential */
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
EXTERN real stiweb_a2[55] INIT(zero55);
EXTERN real sw_a2[10][10];
EXTERN real stiweb_la[550] INIT(zero550);
EXTERN real sw_la[10][10][10];
EXTERN real stiweb_ga[55] INIT(zero55);
EXTERN real sw_ga[10][10];
#endif

#if defined(TERSOFF) || defined(TERSOFFMOD) || defined(BRENNER)
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
#ifdef TERSOFF
EXTERN real ters_chi[45] INIT(one45);
EXTERN real ters_om[45] INIT(one45);
EXTERN real ter_om[10][10];
#ifdef TERSOFF2
EXTERN real ters_ga[55] INIT(zero55);
EXTERN real ter_ga[10][10];
EXTERN real ters_n[55] INIT(zero55);
EXTERN real ter_n[10][10];
EXTERN real ters_c[55] INIT(zero55);
EXTERN real ter_c[10][10];
EXTERN real ter_c2[10][10];
EXTERN real ters_d[55] INIT(zero55);
EXTERN real ter_d[10][10];
EXTERN real ter_d2[10][10];
EXTERN real ters_h[55] INIT(zero55);
EXTERN real ter_h[10][10];
#else
EXTERN real ters_ga[10];
EXTERN real ters_n[10] INIT(zero10);
EXTERN real ters_c[10] INIT(zero10);
EXTERN real ter_c2[10] INIT(zero10);
EXTERN real ters_d[10] INIT(zero10);
EXTERN real ter_d2[10] INIT(zero10);
EXTERN real ters_h[10] INIT(zero10);
#endif
#else  /* TERSOFFMOD */
#ifdef TERSOFFMOD2
EXTERN real ters_eta[55] INIT(zero55);
EXTERN real ter_eta[10][10];
EXTERN real ters_delta[55] INIT(zero55);
EXTERN real ter_delta[10][10];
EXTERN real ters_alpha[55] INIT(zero55);
EXTERN real ter_alpha[10][10];
EXTERN real ters_beta[55] INIT(zero55);
EXTERN real ter_beta[10][10];
EXTERN real ters_c1[55] INIT(zero55);
EXTERN real ter_c1[10][10];
EXTERN real ters_c2[55] INIT(zero55);
EXTERN real ter_c2[10][10];
EXTERN real ters_c3[55] INIT(zero55);
EXTERN real ter_c3[10][10];
EXTERN real ters_c4[55] INIT(zero55);
EXTERN real ter_c4[10][10];
EXTERN real ters_c5[55] INIT(zero55);
EXTERN real ter_c5[10][10];
EXTERN real ters_h[55] INIT(zero55);
EXTERN real ter_h[10][10];
#else
EXTERN real ters_eta[10] INIT(zero10);
EXTERN real ters_delta[10] INIT(zero10);
EXTERN real ters_alpha[10] INIT(zero10);
EXTERN real ters_beta[10] INIT(zero10);
EXTERN real ters_c1[10] INIT(zero10);
EXTERN real ters_c2[10] INIT(zero10);
EXTERN real ters_c3[10] INIT(zero10);
EXTERN real ters_c4[10] INIT(zero10);
EXTERN real ters_c5[10] INIT(zero10);
EXTERN real ters_h[10] INIT(zero10);
#endif
#endif
#endif

#ifdef BRENNER
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

#if defined(COVALENT) || defined(NNBR_TABLE)
EXTERN int neigh_len INIT(NEIGH_LEN_INIT); /* initial neighbor table length */
EXTERN real *neightab_r2cut INIT(NULL);    /* cutoff of neighbor table */
#endif

#ifdef NNBR_TABLE
EXTERN int nnbr_done INIT(0);              /* Flag indicating if nearest neighbors are computed during this time step */
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

EXTERN real     ew_r2_cut INIT(0.0);     /* EWALD r-space cutoff */
#ifdef SM
EXTERN pot_table_t erfc_r_tab;    /* tabulated function erfc/r */
EXTERN pot_table_t cr_pot_tab;    /* tabulated coulomb repulsive potential*/
EXTERN pot_table_t na_pot_tab;    /* tabulated nuclear attraction potential */
EXTERN real tot_sm_es_energy INIT(0.0); /* electrostatic energy per atom for .eng file */
#endif

#if defined(EWALD) || defined(COULOMB) || defined(USEFCS) || defined(SM)
EXTERN real     charge[10] INIT(zero10); /* Charge of atoms */
EXTERN real     coul_eng INIT(14.40);    /* this is e^2/(4*pi*epsilon_0) in eV A */
#endif
#ifdef USEFCS
EXTERN int      fcs_method               INIT(0);
EXTERN real     fcs_tolerance            INIT(1e-3);
EXTERN int      fcs_near_field_flag      INIT(1);
EXTERN real     fcs_rcut                 INIT(0.0);
EXTERN ivektor  fcs_grid_dim             INIT(nullivektor);
EXTERN int      fcs_max_iter             INIT(0);
EXTERN real     fcs_iter_tolerance       INIT(1e-4);
EXTERN real     fcs_pepc_eps             INIT(0.0);
EXTERN real     fcs_pepc_theta           INIT(0.3);
EXTERN int      fcs_pepc_nthreads        INIT(1);
EXTERN int      fcs_fmm_absrel           INIT(2);
EXTERN int      fcs_fmm_dcorr            INIT(0);
EXTERN int      fcs_fmm_do_tune          INIT(0);
EXTERN int      fcs_vmg_max_level        INIT(6);
EXTERN int      fcs_vmg_smooth_steps     INIT(3);
EXTERN int      fcs_vmg_gamma            INIT(2);
EXTERN int      fcs_vmg_near_field_cells INIT(6);
EXTERN int      fcs_vmg_interpol_order   INIT(5);
EXTERN int      fcs_vmg_discr_order      INIT(4);
EXTERN int      fcs_pp3mg_ghosts         INIT(0);
EXTERN int      fcs_pp3mg_degree         INIT(0);
EXTERN int      fcs_pp3mg_max_part       INIT(0);
EXTERN int      fcs_p2nfft_intpol_order  INIT(0);
EXTERN real     fcs_p2nfft_epsI          INIT(0.0);
#endif
#if defined(EWALD) || defined(COULOMB)
EXTERN imd_timer ewald_time;
EXTERN real     ew_kappa;                /* Parameter kappa */
EXTERN real     ew_kcut;                 /* k-space cutoff */
#ifndef SM
EXTERN int      ew_nmax INIT(-1);        /* Number of image boxes */
#else
EXTERN int      ew_nmax INIT(0);        /* Number of image boxes */
#endif
EXTERN int      ew_totk;
EXTERN int      ew_nx;
EXTERN int      ew_ny;
EXTERN int      ew_nz;
EXTERN int      ew_dx;
EXTERN int      ew_dy;
EXTERN int      ew_dz;
EXTERN int      ew_test INIT(0);
EXTERN vektor   *ew_kvek;
EXTERN ivektor  *ew_ivek;
EXTERN real     *ew_expk;
EXTERN real     *coskx;
EXTERN real     *sinkx;
EXTERN real     *cosky;
EXTERN real     *sinky;
EXTERN real     *coskz;
EXTERN real     *sinkz;
EXTERN real     *coskr;
EXTERN real     *sinkr;
EXTERN real     ew_vorf;
EXTERN real     twopi;
EXTERN pot_table_t coul_table; /* one table to hold all coul. and dipole fn */
EXTERN real     coul_res   INIT(0); /* function table resolution */
EXTERN real     coul_begin INIT(0.2); /* start of function table */
EXTERN real     coul_shift;
EXTERN real     coul_fshift;
#endif /* EWALD or COULOMB */
#ifdef DIPOLE
EXTERN int      dp_fix     INIT(0); /* Keep dipoles fixed? */
EXTERN real     dp_mix     INIT(0.8); /* dipole field mixing parameter */
EXTERN real     dp_tol     INIT(1.e-7); /* dipole iteration precision */
EXTERN real     dp_self;     	        /* dipole self field factor */
EXTERN real     *dp_alpha  INIT(NULL); /* in e^2 A^2 / eV^2 */
EXTERN real     *dp_b      INIT(NULL);		/* in eV A / e^2 */
EXTERN real     *dp_c      INIT(NULL);		/* in 1/A */
#endif
#if defined(DIPOLE)|| defined(MORSE)
EXTERN real     *ms_D      INIT(NULL); /* in eV */
EXTERN real     *ms_gamma  INIT(NULL);
EXTERN real     *ms_harm_a  INIT(NULL);
EXTERN real     *ms_harm_b  INIT(NULL);
EXTERN real     *ms_harm_c  INIT(NULL);
EXTERN real     *ms_r2_min INIT(NULL);
EXTERN real     *ms_r0     INIT(NULL); /* in A */
EXTERN real     *ms_shift  INIT(NULL);
EXTERN real     *ms_fshift INIT(NULL);
#endif /* DIPOLE or MORSE */
#if defined(BUCK)
EXTERN real     *bk_shift  INIT(NULL);
EXTERN real     *bk_fshift INIT(NULL);
#endif /* BUCK */
#ifdef EXTF
/* external homogeneous electrostatic field */
EXTERN vektor extf INIT(nullvektor);    /* in e eV */
#endif /* EXTF */

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
EXTERN real shock_speed_l INIT(0.0);              /* velocity of left mirror */
EXTERN real shock_speed_r INIT(0.0);              /* velocity of right mirror */
EXTERN int shock_incr INIT(0);      /* steps to accelerate to full velocity */
EXTERN int shock_mode INIT(2);                   /* type of shock */
#endif

EXTERN int ensemble INIT(ENS_EMPTY);    /* active ensemble type */
EXTERN int loop INIT(0);                /* loop for online visualisation */
EXTERN void (*move_atoms)(void);        /* active integrator routine */

/* global parameters for UNIAX */
#ifdef UNIAX
EXTERN real uniax_inert  INIT(1.0);        /* moment of intertia */
EXTERN vektor uniax_sig  INIT(nullvektor); /* shape of potential */
EXTERN vektor uniax_eps  INIT(nullvektor); /* depth of potential */
EXTERN real uniax_r_cut  INIT(0.0);        /* cutoff radius */
EXTERN real uniax_r2_cut INIT(0.0);        /* square of cutoff radius */
#endif

#ifdef NMOLDYN
EXTERN int nmoldyn_int   INIT(0);
EXTERN int nmoldyn_veloc INIT(1);
#endif

#ifdef DSF
EXTERN int  dsf_int     INIT(0);
EXTERN int  dsf_nk      INIT(0);
EXTERN int  dsf_nkmax   INIT(0);
EXTERN int  *dsf_k0     INIT(NULL);
EXTERN int  *dsf_kdir   INIT(NULL);
EXTERN int  *dsf_kmax   INIT(NULL);
EXTERN real *dsf_weight INIT(NULL);
#endif

#ifdef CBE
EXTERN int num_spus INIT(6);
EXTERN int num_bufs INIT(2);
EXTERN int cbe_pot_steps INIT(100);
EXTERN real cbe_pot_max INIT(20.0);
#endif

EXTERN int myrank   INIT(0);
#ifdef NEB
#define NEB_MAXNREP 100
EXTERN int  neb_nrep INIT(0);
EXTERN int  neb_cineb_start INIT(99999999);
EXTERN int  neb_climbing_image INIT(-1);
EXTERN int  neb_vark_start INIT(99999999);
EXTERN int  neb_eng_int INIT(0);
EXTERN real neb_k INIT(1.0);                /* spring constant */
EXTERN real neb_kmax INIT(0.0);             /* max spring constant for variable spring method*/
EXTERN real neb_kmin INIT(0.0);             /* min constant for variable spring method*/
EXTERN real neb_fnorm INIT(0.0);            /* total force norm */
EXTERN char *neb_outfilename INIT(NULL);    /* name of NEB .eng file */
EXTERN FILE *neb_eng_file INIT(NULL);       /* pointer to NEB .eng file */
EXTERN real neb_image_energies[NEB_MAXNREP] INIT(zero100); /* energies of the individual images */
EXTERN real neb_epot_im[NEB_MAXNREP] INIT(zero100);
EXTERN real neb_ks[NEB_MAXNREP] INIT(zero100);
EXTERN real phi_dl INIT(0.0);
EXTERN real phi_dr INIT(0.0);
EXTERN real phi_lr INIT(0.0);
EXTERN real neb_maxmove INIT(0.0);
#endif


#ifdef BBOOST
#define BB_MAXNEIGH     30
#define BB_MAXCOMPONENTS 5                 /* maximal number of alloying components  = ntypes */
EXTERN real bb_btime INIT(0.0);            /* number of timestep * boostfactor during boost */
EXTERN real bb_tot_bV INIT(0.0);           /* magnitude of boost potential */
EXTERN real sum_bfcr INIT(0.0);           /* summation of boost factor */
EXTERN real C_x INIT(0.0);		/* section of shut down function */
EXTERN real C_x2 INIT(0.0);		/* section of shut down function */
EXTERN real B_x INIT(0.0);		/* section of shut down function */
EXTERN real B_x2 INIT(0.0);		/* section of shut down function */
EXTERN real A_e_max INIT(0.0);		/* shut down function */
EXTERN real p1_2 INIT(0.98);           	/* curvature controller of the boost potential */
EXTERN real bb_rcut INIT(1.0);		/* the cut_off for the bb_neight */
EXTERN real *bb_neightab_r2cut INIT(NULL);  /* the cut_off for the bb_neight */
EXTERN real bb_epscrit[BB_MAXCOMPONENTS][BB_MAXCOMPONENTS]; /* largest fraction of bondlength to determine whether a bond is broken */
EXTERN real bb_eps[BB_MAXCOMPONENTS][BB_MAXCOMPONENTS];
EXTERN real bb_eps2[BB_MAXCOMPONENTS][BB_MAXCOMPONENTS];
EXTERN real bb_epsold[BB_MAXCOMPONENTS][BB_MAXCOMPONENTS];
EXTERN int bb_relaxsteps_max INIT(1000); /* maximal relaxed time for the beginning a boost activity */
EXTERN int bb_relaxsteps INIT(0);	/* count the relaxed time */
EXTERN int bb_shdn_max INIT(200);	/* maximal shutdown time that bond fraction/strain exceed bb_epscrit */
EXTERN int bb_shdn INIT(0);		/* count for the bb_shdn_max */
EXTERN int bb_under_max INIT(200);	/* maximal time under boosting that bond fraction/strain not-exceed the bb_epscrit */
EXTERN int bb_under INIT(0);		/* count for the bb_under_max */
EXTERN int bflag1 INIT(0);		/* flag to switch the boost method ( = 0 is non-boosting, = 1 is to launch the boosting) */
EXTERN int bflag2 INIT(0);		/* flag to switch the testing mode or safe mode for a chosen bb_epscrit */
EXTERN int bflag3 INIT(0);		/* flag for time windows used to apply minimization to boost MD */
#endif

#ifdef KIM
EXTERN str255 kim_model_name;
EXTERN char **kim_el_names INIT(NULL);
EXTERN imd_kim_t kim;
#endif
