/******************************************************************************
*
* prototypes.h -- Function prototypes for the imd Package
*
*
******************************************************************************/


int  main (int argc, char **argv);
void main_loop(void);
void read_parameters(int argc, char **argv);
void broadcast_params(void);
void check_parameters_complete(void);
void read_potential(str255 potfilename);
void read_atoms(str255 infilename);
void generate_atoms(str255 infilename);
void init_cubic(void);
void generate_fcc(int maxtyp);
void generate_lav(void);
void init_hex(void);
void generate_hex(void);
void init(void);
void usage(void);
void getparamfile(char *paramfname, int sim);
void make_box(void);
ivektor maximal_cell_dim(void);
void init_cells(void);
void make_cell_lists(void);
void error(char *mesg);
#ifdef TWOD
ivektor cell_coord(real x, real y);
#else
ivektor cell_coord(real x, real y, real z);
#endif
vektor back_into_box(vektor pos);
void move_atom(ivektor cellc, cell *from, int index);
void alloc_cell(cell *thecell, int count);
void calc_forces(void);
void move_atoms_nve(void);
void move_atoms_mik(void);
void move_atoms_nvt(void);
void move_atoms_npt_iso(void);
void move_atoms_npt_axial(void);
void move_atoms_and(void);
void move_atoms_mc(void);
void move_atoms_frac(void);
void move_atoms_nvx(void);
void move_atoms_msd(void);
void move_atoms_stm(void);
void do_boundaries(void);
void fix_cells(void);
void calc_properties(void);
void write_cell(FILE *out, cell *p);
void write_config(int steps);
void write_properties(int steps);
void write_surroundings(int steps);
void write_distrib(int steps);
void write_pic_cell( cell *p, FILE *out ); 
void write_vrml_cell( cell *p, FILE *out ); 
void write_pictures(int steps);
void write_pictures_raw(int steps);
void write_pictures_atoms(int steps);
void write_pictures_bins(int steps);
void write_vrmls(int steps);
void check_socket(int steps);
void do_forces(cell *p, cell *q, vektor pbc);
void maxwell(real TEMP);
float gasdev(long *idum);
float ran1(long *idum);

void write_conf_using_sockets(void);
void write_ras_using_sockets(void);
void write_distrib_using_sockets(void);
#ifdef TWOD
void write_rgb_picture_to_socket(void);
#endif

#ifdef MC
real mc_epot_diff( vektor old_pos, vektor new_pos, 
                   int p_num, int p_typ, ivektor cellc );
real mc_epot_atom( vektor pos, int p_num, int p_typ, ivektor cellc );
real mc_epot_part( void );
void one_mc_step();
#endif

#ifdef MPI
void setup_mpi_topology( void );
void init_mpi(int argc,char **argv);
void shutdown_mpi(void);
int  cpu_coord(ivektor cellc);
int  cpu_grid_coord(ivektor cellc);
#ifdef TWOD
ivektor local_cell_coord(real x, real y);
#else
ivektor local_cell_coord(real x, real y, real z);
#endif
void recv_atoms(void);
void send_atoms(int mode);
void send_atoms_force(void);
void send_atoms_ar(void);
void send_forces(void);
void send_forces_full(void);
#ifdef TWOD
void copy_cell_force(int i, int j, int k, int l);
void copy_atoms_force(msgbuf *b, int k, int l);
void move_atoms_force(msgbuf *b, int k, int l);
#else
vektor global_pbc(int i,int j, int k);
void copy_cell_force(int i, int j, int k, int l, int m, int n);
void add_cell_force(int i, int j, int k, int l, int m, int n);
void move_atoms_force(msgbuf *b, int k, int l, int m);
void copy_atoms_force(msgbuf *b, int k, int l, int m);
void add_forces(msgbuf *b, int k, int l, int m);
void copy_forces(msgbuf *b, int k, int l, int m);
#endif
void recv_cell(cell *p, int from_cpu, int tag);
void send_cell(cell *p, int to_cpu, int tag);
int  irecv_buf( msgbuf *b, int from_cpu, MPI_Request *req);
int  isend_buf( msgbuf *b, int to_cpu,   MPI_Request *req);
void mpi_addtime(double *timer);
void copy_atoms_buf(msgbuf *to, msgbuf *from);
void process_buffer(msgbuf *b, int mode);
void setup_buffers(void);
void empty_buffer_cells(void);
void copy_one_atom(msgbuf *to, cell *from, int index);
void empty_mpi_buffers(void);
ivektor cpu_coord_v(ivektor cellc);
#ifdef SAVEMEM
ivektor global_cell_coord(ivektor coords);
void fix_cells_by_cell(void);
void send_atoms_by_cell(void);
void send_recv_cell(int i,int j,int k,int l,int m,int n);
#endif
#endif

#ifdef PACX
ivektor my_cart_coords(int myid);
void my_cart_rank(ivektor my_coord);
#endif

#ifdef HOMDEF
void shear_sample(void);
void expand_sample(void);
#endif
#ifdef DEFORM
void deform_sample(void);
#endif

#ifdef DISLOC
void reset_Epot_ref();
void write_demmaps(int steps);
void write_dspmaps(int steps);
void update_ort_ref(void);
#endif

#ifdef EAM
void do_forces_eam_1(cell *p, cell *q, vektor pbc);
void do_forces_eam_2(cell *p);
#endif

#ifdef EAM2
void eam2_read_core_pot(str255 core_pot_filename);
void eam2_read_embedding_energy(str255 emb_E_filename);
void eam2_read_atomic_rho(str255 at_rho_filename);
void eam2_do_forces1(cell *p, cell *q, vektor pbc);
void eam2_do_forces2(cell *p, cell *q, vektor pbc);
void copy_cell_eam2_rho_h(int i, int j, int k, int l, int m, int n);
void move_eam2_rho_h(msgbuf *b, int k, int l, int m);
void copy_eam2_rho_h(msgbuf *b, int k, int l, int m);
void send_eam2_rho_h(void);
#endif

#ifdef COVALENT
void do_neightab(cell *p, cell *q, vektor pbc);
neightab *alloc_neightab(neightab *neigh, int count);
#endif

#ifdef TTBP
void do_forces_ttbp(cell *p);
void read_ttbp_potential(str255 ttbp_potfilename);
#endif

#ifdef TERSOFF
void do_forces_tersoff(cell *p);
#endif

#ifdef UNIAX
void do_forces_uniax(cell *p, cell *q, vektor pbc);
void gay_berne ( vektor r12, vektor e1, vektor e2, 
		 real rsqr, vektor s1, vektor w1, 
		 real *pot12, vektor *force12, 
		 vektor *torque12, vektor *torque21 );
#endif

#if defined(CORRELATE) || defined(MSQD)
void alloc_correl(int, int);
void correlate(int istep, int refstep, unsigned seqnum);
#endif

#ifdef CORRELATE
void write_add_corr(int it, int steps, unsigned seqnum);
#endif

#ifdef TRANSPORT
void write_temp_dist(int steps);
#endif
#ifdef RNEMD
void rnemd_heat_exchange();
#endif

#ifdef STRESS_TENS
void write_press_dist(int steps);
void write_press_dist_shock(int steps);
void write_press_atoms(int steps);
#endif

/* generate quasicrystal */
#ifdef QUASI
real r2 (real ai, real aj, real bi, real bj, real b, real ci, real cj,
	 real c, int ks, real gi, real gj);
real r3 (real ai, real aj, real ak, real bi, real bj, real bk, real b,
	 real ci, real cj, real ck, int ks1, int ks2, real gi, real gj,
	 real gk);
real det(real ai, real aj, real ak, real bi, real bj, real bk, real ci, 
	 real cj, real ck);
void sortin (int ifeld[]);
void adjust(void);
void decorate(int i, int j, int k);
void locate(real x, real y, real z, int i, int j, int k);
void borders(void);
#endif
void display_conf(int steps);
#ifdef EFILTER
void efwrite_config(int steps);
#endif






