
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
* prototypes.h -- Function prototypes for IMD
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/


/* the main routines - files imd.c, imd_main_2/3d.c */
int  main(int argc, char **argv);
void main_loop(void);
void init(void);

/* read parameters - file imd_param.c */
void read_parameters(int argc, char **argv);
void getparamfile(char *paramfname, int sim);
void check_parameters_complete(void);
void broadcast_params(void);

/* read and access potential tables - file imd_potential.c */
void read_pot_table1( pot_table_t *pt, char *filename );
void read_pot_table2( pot_table_t *pt, char *filename, int cols );
void pair_int_monolj(real *pot, real *grad, real r2);
void pair_int2  (real*, real*, int*, pot_table_t*, int, int, real);
void pair_int3  (real*, real*, int*, pot_table_t*, int, int, real);
void   val_func2(real*,        int*, pot_table_t*, int, int, real);
void   val_func3(real*,        int*, pot_table_t*, int, int, real);
void deriv_func2(       real*, int*, pot_table_t*, int, int, real);
void deriv_func3(       real*, int*, pot_table_t*, int, int, real);

/* read configuration - files imd_io_2/3d.c */
void   read_atoms(str255 infilename);
#ifdef MPI
void   recv_atoms(void);
#endif

/* generate configuration - file imd_generate.c */
void generate_atoms(str255 infilename);
void init_cubic(void);
void generate_fcc(int maxtyp);
void generate_lav(void);
void init_hex(void);
void generate_hex(void);
#ifdef QUASI /* generate quasicrystal - file imd_qc.c */
void init_qc(void);
void generate_qc(void);
#endif

/* miscellaneous routines - files imd_misc.c, imd_time.c, imd_maxwell.c */
void usage(void);
void error(char *msg);
void warning(char *msg);
void imd_start_timer(imd_timer *timer);
void imd_stop_timer(imd_timer *timer);
void maxwell(real TEMP);
int  endian(void);

/* start and stop MPI - files imd_mpi_util.c, imd_geom_mpi_*.c */
#ifdef MPI
void init_mpi(int *argc_pointer, char **argv);
void setup_mpi_topology(void);
void shutdown_mpi(void);
#endif

/* manage MPI buffers and buffer cells - file imd_mpi_util.c */
#ifdef MPI
void setup_buffers(void);
void empty_mpi_buffers(void);
#ifdef SR
int  sendrecv_buf(msgbuf *send, int to_cpu, 
                  msgbuf *recv, int from_cpu, MPI_Status *status);
#else
int  irecv_buf(msgbuf *b, int from_cpu, MPI_Request *req);
int  isend_buf(msgbuf *b, int to_cpu,   MPI_Request *req);
#endif
void empty_buffer_cells(void);
#ifdef SAVEMEM  /* imd_savemem_3d.c */
void dealloc_buffer_cells(void);
#endif
#endif

/* make and maintain cells and their geometry - files imd_geom_*.c */
ivektor maximal_cell_dim(void);
vektor  back_into_box(vektor pos);
vektor  vec_prod(vektor u, vektor v);
void make_box(void);
void init_cells(void);
void make_cell_lists(void);
void check_pairs(void);
void move_atom(cell *to, cell *from, int index);
void alloc_cell(cell *thecell, int count);
#ifdef TWOD
ivektor       cell_coord(real x, real y);
ivektor local_cell_coord(real x, real y);
#else
ivektor       cell_coord(real x, real y, real z);
ivektor local_cell_coord(real x, real y, real z);
#endif
#ifdef MPI
int     cpu_coord(ivektor cellc);
ivektor cpu_coord_v(ivektor cellc);
int     cpu_grid_coord(ivektor cellc);
#endif

/* force computation - files imd_main_*.c, imd_forces_*.c */
void calc_forces(void);
void do_forces(cell*, cell*, vektor, real*, real*, real*, real*, real*);
#ifdef COVALENT
void do_forces2(cell*, real*, real*, real*, real*, real*);
#endif
#ifdef EAM2
void do_forces_eam2(cell*, cell*, vektor, real*, real*, real*, real*, real*);
#endif
#ifdef TERSOFF
void init_tersoff(void);
#endif
#ifdef UNIAX
void gay_berne ( vektor r12, vektor e1, vektor e2, 
		 real rsqr, vektor s1, vektor w1, 
		 real *pot12, vektor *force12, 
		 vektor *torque12, vektor *torque21 );
#endif

/* communication for force computation - files imd_comm_force_2/3d.c */
#ifdef MPI
#ifdef TWOD
void send_cells (void (*copy_func)  (int, int, int, int),
                 void (*pack_func)  (msgbuf*, int, int),
                 void (*unpack_func)(msgbuf*, int, int));
void send_forces(void (*add_func)   (int, int, int, int),
                 void (*pack_func)  (msgbuf*, int, int),
                 void (*unpack_func)(msgbuf*, int, int));
void copy_cell    ( int j, int k, int l, int m );
void pack_cell    ( msgbuf *b, int j, int k );
void unpack_cell  ( msgbuf *b, int j, int k );
void add_forces   ( int j, int k, int l, int m );
void pack_forces  ( msgbuf *b, int j, int k );
void unpack_forces( msgbuf *b, int j, int k );
#else  /* 3D */
void send_cells (void (*copy_func)  (int, int, int, int, int, int),
                 void (*pack_func)  (msgbuf*, int, int, int),
                 void (*unpack_func)(msgbuf*, int, int, int));
void send_forces(void (*copy_func)  (int, int, int, int, int, int),
                 void (*pack_func)  (msgbuf*, int, int, int),
                 void (*unpack_func)(msgbuf*, int, int, int));
void copy_cell    ( int k, int l, int m, int r, int s, int t );
void pack_cell    ( msgbuf *b, int k, int l, int m );
void unpack_cell  ( msgbuf *b, int k, int l, int m );
void add_forces   ( int k, int l, int m, int r, int s, int t );
void pack_forces  ( msgbuf *b, int k, int l, int m);
void unpack_forces( msgbuf *b, int k, int l, int m );
#ifdef EAM2
void copy_rho_h  ( int k, int l, int m, int r, int s, int t );
void pack_rho_h  ( msgbuf *b, int k, int l, int m );
void unpack_rho_h( msgbuf *b, int k, int l, int m );
#endif
#ifdef MONOLJ   /* imd_main_mpi_3d.c */
vektor global_pbc(int i,int j, int k);
#endif
#ifdef SAVEMEM  /* imd_savemem_3d.c */
void send_cells_by_cell(void);
void send_recv_cell(int i, int j, int k, int l, int m, int n);
void send_cell_force(cell *p, int to_cpu, int tag);
void recv_cell_force(cell *p, int from_cpu,int tag);
#endif
#endif /* 3D  */
#endif /* MPI */

/* integrators - file imd_integrate.c */
void move_atoms_nve(void);
void move_atoms_mik(void);
void move_atoms_nvt(void);
void calc_dyn_pressure(void);
void move_atoms_npt_iso(void);
void move_atoms_npt_axial(void);
void move_atoms_mc(void);
void move_atoms_frac(void);
void move_atoms_nvx(void);
void move_atoms_msd(void);
void move_atoms_stm(void);

/* fix distribution on cells - files imd_main_*.c, imd_mpi_util.c */
void do_boundaries(void);
void fix_cells(void);
#ifdef MPI
void copy_atoms_buf(msgbuf *to, msgbuf *from);
void copy_one_atom(msgbuf *to, cell *from, int index, int delete);
void process_buffer(msgbuf *b, cell *p);
void send_atoms(void);
#ifdef SAVEMEM   /* imd_savemem_3d.c */
void fix_cells_by_cell(void);
#endif
#endif
/* write properties - file imd_io_*.c */
void write_end_file(int steps);
void write_eng_file_header(void);

/* write configurations - files imd_io.c, imd_io_*.c */
void write_itr_file(int fzhlr, int steps);
void write_config(int steps);
void write_config_select(int fzhlr, char *suffix,
  void (*write_atoms_fun)(FILE *out), void (*write_header_fun)(FILE *out));
void write_atoms_config(FILE *out);
void write_header_config(FILE *out);
void write_atoms_pic(FILE *out); 
void write_header_pic(FILE *out); 
#ifdef DISLOC
void write_atoms_dem(FILE *out);
void write_header_dem(FILE *out);
void write_atoms_dsp(FILE *out);
void write_header_dsp(FILE *out);
#endif
#ifdef EFILTER
void write_atoms_ef(FILE *out);
void write_header_ef(FILE *out);
#endif
#ifdef STRESS_TENS
void write_atoms_press(FILE *out);
void write_header_press(FILE *out);
void calc_tot_presstens(void);
#endif
#ifdef AVPOS
void write_atoms_avp(FILE *out);
void write_header_avp(FILE *out);
#endif

void reduce_displacement(vektor *d);
#ifdef MPI
void recv_cell_old(cell *p, int from_cpu, int tag);
void send_cell_old(cell *p, int to_cpu, int tag);
void recv_cell(cell *p, int from_cpu, int tag);
void send_cell(cell *p, int to_cpu, int tag);
#endif
#ifdef USE_SOCKETS
void write_conf_using_sockets(void);
void write_ras_using_sockets(void);
#endif

/* write distributions - files imd_histogram.c, socket_io.c */
void make_histograms(hist_t *hist);
void write_distrib(int steps);
void write_distrib_header(FILE *out, char *type);
#ifdef STRESS_TENS
void write_press_dist(int steps);
void write_press_dist_header(FILE *out);
#endif
#ifdef SHOCK
void write_press_dist_shock(int steps);
void write_press_dist_shock_header(FILE *out);
#endif
#ifdef USE_SOCKETS
void write_distrib_using_sockets(void);
#endif

/* write pictures - files imd_io_*.c, socket_io.c */
void write_pictures(int steps);
void write_pictures_bitmap(int steps);
#ifdef USE_SOCKETS
#ifdef TWOD
void write_rgb_picture_to_socket(void);
#endif
#endif

/* shear, deform of load sample - files imd_deform.c, imd_load.c */
#ifdef HOMDEF
void shear_sample(void);
void expand_sample(void);
#endif
#ifdef DEFORM
void deform_sample(void);
#endif
#ifdef FRAC
void load_sample(void);
#endif

/* support for neighbor tables - files imd_alloc.c, imd_forces_covalent.c */
#ifdef COVALENT
void do_neightab(cell *p, cell *q, vektor pbc);
neightab *alloc_neightab(neightab *neigh, int count);
#endif

#ifdef EWALD
/* support for computation of Coulomb forces */
void do_forces_ewald_real(void);
void do_forces_ewald_fourier(void);
void init_ewald(void);
real erfc1(real x);
#endif

/* support for dislocations - file imd_io.c */
#if defined(DISLOC) || defined(AVPOS)
void reset_Epot_ref(void);
void update_ort_ref(void);
#endif

/* support for average over positions */
#ifdef AVPOS
void add_positions(void);
#endif

/* support for correlation functions - file imd_correl.c */
#if defined(CORRELATE) || defined(MSQD)
void alloc_correl(int, int);
void correlate(int istep, int refstep, unsigned seqnum);
#endif
#ifdef CORRELATE
void write_add_corr(int it, int steps, unsigned seqnum);
#endif

/* support for heat transport - file imd_transport.c */
#ifdef TRANSPORT
void write_temp_dist(int steps);
#endif
#ifdef RNEMD
void rnemd_heat_exchange();
#endif

/* support for socket communication - file socket_io.c */
#ifdef USE_SOCKETS
void check_socket(int steps);
#endif
