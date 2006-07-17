
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2006 Institute for Theoretical and Applied Physics,
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
void read_command_line(int argc, char **argv);
void read_parameters(void);
void getparamfile(char *paramfname, int sim);
void check_parameters_complete(void);
void broadcast_params(void);

/* read and access potential tables - file imd_potential.c */
void read_pot_table ( pot_table_t*, char*, int, int );
void read_pot_table1( pot_table_t*, int, char*, FILE*, int );
void read_pot_table2( pot_table_t*, int, char*, FILE*, int );
void free_pot_table ( pot_table_t*);
#ifdef LINPOT
void make_lin_pot_table( pot_table_t, lin_pot_table_t* );
#endif
#ifdef PAIR
void pair_int_lj(real *pot, real *grad, int p_typ, int q_typ, real r2);
void pair_int_morse(real *pot, real *grad, int p_typ, int q_typ, real r2);
void pair_int_buck(real *pot, real *grad, int p_typ, int q_typ, real r2);
#ifdef EWALD
void pair_int_ewald(real *pot, real *grad, int p_typ, int q_typ, real r2);
#endif
#endif
void pair_int2  (real*, real*, int*, pot_table_t*, int, int, real);
void pair_int3  (real*, real*, int*, pot_table_t*, int, int, real);
void   val_func2(real*,        int*, pot_table_t*, int, int, real);
void   val_func3(real*,        int*, pot_table_t*, int, int, real);
void deriv_func2(       real*, int*, pot_table_t*, int, int, real);
void deriv_func3(       real*, int*, pot_table_t*, int, int, real);
void init_pre_pot(void);
void create_pot_table(pot_table_t *pt);
void test_potential(pot_table_t, char*, int);
#if   defined(FOURPOINT)
void init_fourpoint(pot_table_t*, int);
#elif defined(SPLINE)
void init_spline(pot_table_t*, int, int);
#else
void init_threepoint(pot_table_t*, int);
#endif

/* read configuration - files imd_io_2/3d.c */
void read_atoms(str255 infilename);
void read_atoms_cleanup(void);  
#ifdef MPI
void recv_atoms(void);
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
void error_str(char *msg, char *str);
void error_str_str(char *msg, char *str1, char *str2);
void warning(char *);
void warning_str(char *, char *);
void warning_str_str(char *, char *, char *);
void imd_start_timer(imd_timer *timer);
void imd_stop_timer(imd_timer *timer);
void maxwell(real TEMP);
int  endian(void);
integer SwappedInteger(integer);
float   SwappedFloat  (float  );
double  SwappedDouble (double );

/* start and stop MPI - files imd_mpi_util.c, imd_geom_mpi_*.c */
#ifdef MPI
void init_mpi(int *argc_pointer, char **argv);
void shutdown_mpi(void);
void alloc_msgbuf(msgbuf*, int);
void realloc_msgbuf(msgbuf*, int);
void free_msgbuf(msgbuf*);
#endif
void setup_mpi_topology(void);
void calc_cpu_dim(void);

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
void init_io(void);
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
#ifdef VEC
void move_atom_mini(minicell *, minicell *, int);
void insert_atom(minicell *, cell *, int);
void alloc_minicell(minicell *, int);
#endif
void alloc_cell(cell *thecell, int count);
#ifdef TWOD
ivektor cell_coord(real x, real y);
#else
ivektor cell_coord(real x, real y, real z);
#endif
ivektor local_cell_coord(ivektor cellc);
int     cpu_coord(ivektor cellc);
ivektor cpu_coord_v(ivektor cellc);
int     cpu_grid_coord(ivektor cellc);

/* force computation - files imd_main_*.c, imd_forces_*.c */
void calc_forces(int steps);
void do_forces(cell*, cell*, vektor, real*, real*, real*, real*, real*, real*, real*, real*);
#ifdef COVALENT
void do_forces2(cell*, real*, real*, real*, real*, real*, real*, real*, real*);
#endif
#ifdef EAM2
void do_forces_eam2(cell*, cell*, vektor, real*, real*, real*, real*, real*, real*, real*);
void do_embedding_energy(void);
#endif
#ifdef NBLIST
int  estimate_nblist_size(void);
void make_nblist(void);
void check_nblist(void);
#endif
#ifdef MEAM
void init_meam(void);
#endif
#ifdef KEATING
void init_keating(void);
#endif
#ifdef STIWEB
void init_stiweb(void);
void pair_int_stiweb(real *pot, real *grad, int p_typ, int q_typ, real r2);
#endif
#ifdef TTBP
void init_ttbp(void);
#endif
#ifdef TERSOFF
void init_tersoff(void);
void pair_int_tersoff(real *pot, int p_typ, int q_typ, real r2);
#endif
#ifdef UNIAX
void gay_berne ( vektor r12, vektor e1, vektor e2, 
		 real rsqr, vektor s1, vektor w1, 
		 real *pot12, vektor *force12, 
		 vektor *torque12, vektor *torque21 );
#endif

/* communication for force computation - files imd_comm_force_2/3d.c */
#if defined(MPI) || defined(NBLIST)
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
void send_cells (void (*copy_func)  (int, int, int, int, int, int, vektor),
                 void (*pack_func)  (msgbuf*, int, int, int, vektor),
                 void (*unpack_func)(msgbuf*, int, int, int));
void send_forces(void (*copy_func)  (int, int, int, int, int, int),
                 void (*pack_func)  (msgbuf*, int, int, int),
                 void (*unpack_func)(msgbuf*, int, int, int));
void copy_cell    ( int k, int l, int m, int r, int s, int t, vektor v );
void pack_cell    ( msgbuf *b, int k, int l, int m, vektor v );
void unpack_cell  ( msgbuf *b, int k, int l, int m );
void add_forces   ( int k, int l, int m, int r, int s, int t );
void pack_forces  ( msgbuf *b, int k, int l, int m);
void unpack_forces( msgbuf *b, int k, int l, int m );
#ifdef EAM2
void copy_dF       ( int k, int l, int m, int r, int s, int t, vektor );
void add_rho       ( int k, int l, int m, int r, int s, int t );
void pack_dF       ( msgbuf *b, int k, int l, int m, vektor );
void pack_rho      ( msgbuf *b, int k, int l, int m );
void unpack_dF     ( msgbuf *b, int k, int l, int m );
void unpack_add_rho( msgbuf *b, int k, int l, int m );
#endif
#endif /* 3D  */
#endif /* MPI or NBLIST */

/* integrators - file imd_integrate.c */
void move_atoms_nve(void);
void move_atoms_mik(void);
void move_atoms_nvt(void);
void calc_dyn_pressure(void);
void move_atoms_npt_iso(void);
void move_atoms_npt_axial(void);
void move_atoms_mc(void);
void move_atoms_frac(void);
void move_atoms_sllod(void);
void move_atoms_nvx(void);
void move_atoms_msd(void);
void move_atoms_stm(void);
void move_atoms_ftg(void);
void move_atoms_finnis(void);
#ifdef SHOCK
void calc_pxavg(void);
#endif

/* fix distribution on cells - files imd_main_*.c, imd_mpi_util.c */
void do_boundaries(void);
void fix_cells(void);
#ifdef MPI
void copy_atoms_buf(msgbuf *to, msgbuf *from);
void copy_one_atom (msgbuf *to, int to_cpu, minicell *from, int index,int del);
void copy_atom     (msgbuf *to, int to_cpu, cell *from, int index );
void process_buffer(msgbuf *b,  cell *p);
void send_atoms(void);
#endif
/* write properties - file imd_io_*.c */
void write_eng_file(int steps);
void write_eng_file_header(void);

/* write configurations - files imd_io.c, imd_io_*.c */
void read_box(str255);
int  read_header(header_info_t *, str255);
#ifdef MPI
void broadcast_header(header_info_t *);
#endif
void flush_outbuf(FILE *out, int *len, int tag);
void write_itr_file(int fzhlr, int steps,char *suffix);
void write_config(int fzhlr, int steps);
#ifdef RELAX
void write_ssconfig(int steps);
#endif
void write_config_select(int fzhlr, char *suffix,
  void (*write_atoms_fun)(FILE *out), void (*write_header_fun)(FILE *out));
void write_atoms_config(FILE *out);
void write_header_config(FILE *out);
void write_atoms_pic(FILE *out); 
void write_header_pic(FILE *out); 
void write_atoms_pos(FILE *out); 
void write_header_pos(FILE *out); 
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
#ifdef NNBR
void write_atoms_nb(FILE *out);
void write_header_nb(FILE *out);
#endif
#ifdef WRITEF
void write_atoms_wf(FILE *out);
void write_header_wf(FILE *out);
#endif
#ifdef STRESS_TENS
void write_atoms_press(FILE *out);
void write_header_press(FILE *out);
void calc_tot_presstens(void);
#endif
#ifdef REFPOS
void init_refpos(void);
#endif
#ifdef AVPOS
void write_atoms_avp(FILE *out);
void write_header_avp(FILE *out);
void write_avpos_itr_file(int fzhlr, int steps);
#endif
#ifdef FORCE
void write_atoms_force(FILE *out);
void write_header_force(FILE *out);
#endif
#ifdef MSQD
void write_atoms_sqd(FILE *out);
void write_header_sqd(FILE *out);
#endif
void reduce_displacement(vektor *d);
#ifdef MPI
void recv_cell(minicell *p, int from_cpu, int tag);
void send_cell(minicell *p, int to_cpu, int tag);
#endif

#ifdef USE_SOCKETS
void init_socket(void);
int  connect_visualization(void);
void close_socket(void);
void check_socket(void);
void vis_init(void);
void vis_write_config_quit(void);
void vis_send_msg(char *msg);
void vis_check_atoms_flags(void);
void vis_write_atoms_buf(int *len, int tag);
void vis_write_atoms_fun(void);
void vis_init_atoms(void);
void vis_write_atoms(void);
void vis_change_params(void);
void vis_change_params_deform(integer flag);
void vis_restart_simulation(void);
void write_ras_using_sockets(void);
void write_conf_using_sockets(void);
void write_distrib_using_sockets(void);
void write_rgb_picture_to_socket(void);
#endif

/* write distributions - file imd_distrib.c */
void make_distrib_select(dist_t*, int, int*, 
                         void (*fun)(float*, cell*, int));
void write_distrib_select(dist_t*, int, int, int, char*, char*);
void write_distrib_header(FILE*, dist_t*, int, int, char*);
void write_distrib(int);
void dist_Ekin_fun       (float*, cell*, int);
void dist_Epot_fun       (float*, cell*, int);
void dist_Ekin_long_fun  (float*, cell*, int);
void dist_Ekin_trans_fun (float*, cell*, int);
void dist_Ekin_comp_fun (float*, cell*, int);
void dist_shock_shear_fun(float*, cell*, int);
void dist_shear_aniso_fun(float*, cell*, int);
void dist_press_fun      (float*, cell*, int);
void dist_presstens_fun  (float*, cell*, int);
void dist_presstens_xx_fun  (float*, cell*, int);
void dist_presstens_yy_fun  (float*, cell*, int);
#ifndef TWOD
void dist_presstens_zz_fun  (float*, cell*, int);
void dist_presstens_yz_fun  (float*, cell*, int);
void dist_presstens_zx_fun  (float*, cell*, int);
#endif
void dist_presstens_xy_fun  (float*, cell*, int);
void dist_pressoff_fun   (float*, cell*, int);
void dist_pressxy_fun   (float*, cell*, int);
void dist_pressyz_fun   (float*, cell*, int);
void dist_presszx_fun   (float*, cell*, int);
void dist_vxavg_fun   (float*, cell*, int);

#ifdef ATDIST
void   init_atdist(void);
void update_atdist(void);
void  write_atdist(void);
void write_atoms_atdist_pos(FILE*);
void write_header_atdist_pos(FILE*);
#endif

#ifdef DIFFPAT
void   init_diffpat(void);
void update_diffpat(int );
void  write_diffpat(void);
#endif

/* write pictures - files imd_io_*.c */
void write_pictures(int steps);
void write_pictures_bitmap(int steps);

/* shear, deform of load sample - files imd_deform.c, imd_load.c */
#ifdef HOMDEF
void shear_sample(void);
void expand_sample(void);
#ifdef TWOD
void lin_deform(vektor,vektor,real);
#else
void lin_deform(vektor,vektor,vektor,real);
#endif
void relax_pressure(void);
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
void increase_neightab(neightab *neigh, int count);
#endif

#ifdef EWALD
/* support for computation of Coulomb forces */
void do_forces_ewald(int);
void do_forces_ewald_real(void);
void do_forces_ewald_fourier(void);
void init_ewald(void);
real erfc1(real x);
#endif

#ifdef EPITAX
/* support for molecular beam epitaxy - file imd_epitax.c */
void create_atom(int type, real mass, real temp);
void delete_atoms(void);
real substrate_level(void);
void calc_poteng_min(void);
void check_boxheight(void);
ivektor cell_map(ivektor cellc);
#endif

/* support for dislocations - file imd_io.c */
#ifdef DISLOC
void reset_Epot_ref(void);
void update_ort_ref(void);
#endif

/* support for average over positions */
#ifdef AVPOS
void update_avpos(void);
void add_positions(void);
#endif

/* support for correlation functions - file imd_correl.c */
#if defined(CORRELATE) || defined(MSQD)
void init_correl(int, int);
void alloc_correl(int, int);
void correlate(int istep, int refstep, unsigned seqnum);
void write_msqd(int);
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

#ifdef CG
void reset_cg(void);
void cg_step(int steps);
void write_cgconfig(int steps);
real fonedim (real) ;
real fonedim_sd (real) ;
void cg_calcgamma(void);
void set_hg(void);
void calc_fnorm(void);
void calc_fnorm_g_h(void);
void move_atoms_cg(real);
void move_atoms_sd(real);
int linmin();
int mnbrak(real *ax, real *bx, real *cx, real *fa, real *fb, real *fc);
int brent(real ax, real bx, real cx, real fa, real fb, real fc,real *alphamin);
#endif

#ifdef ACG
void acg_step(int steps);
int findalpha();

#endif

#ifdef NMOLDYN
void init_nmoldyn(void);
void write_nmoldyn(int);
#endif

