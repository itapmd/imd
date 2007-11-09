%module IMD
%{
#include "imd.h"
%}

%include "constants.h"
%include "version.h"
%include "config.h"

typedef double real;
typedef struct {real x; real y; real z; } vektor;
typedef struct {int  x; int  y; int  z; } ivektor;
typedef char str255[255];

int  read_parameters(char *, int);
void setup_potentials(void);
void generate_atoms(char *);
void read_atoms(char *);
void maxwell(real);

long natoms;
int  num_cpus;
ivektor box_param;
real box_unit;
real *masses;
vektor box_x;
vektor box_y;
vektor box_z;

int  imdrestart;
real temperature;
int  do_maxwell;
int  myid;
int  eng_int, dist_int, checkpt_int, pic_int;
int  steps, steps_min, steps_max;
int  ensemble;
real relax_rate;

str255 infilename;
str255 outfilename;

real RealGetElm( real *p, int i);
void RealSetElm( real *p, int i, real val); 
int  IntGetElm( int *p, int i);
void IntSetElm( int *p, int i, int val);

void make_box(void);
void calc_forces(int);
void move_atoms(void);
void write_config(int,int);
void write_eng_file(int steps);
void write_eng_file_header(void);
void write_distrib(int);
void write_pictures(int);

void fix_cells(void);
void close_files(void);
void check_write(void);
int  check_stop(void);
int  check_walltime(void);
int  watch_int, stop_int;
real maxwalltime;

#ifdef MPI
void init_mpi(void);
void shutdown_mpi(void);
#endif

#ifdef NBLIST
void check_nblist(void);
#endif

#ifdef TEMPCONTROL
void increment_temperature(void);
#endif

#ifdef GLOK
void update_glok(void);
#endif

#ifdef EFILTER
void write_config_ef(int nr);
int  ef_checkpt_int;
#endif

#ifdef NNBR
void write_config_nb(int nr);
int  nb_checkpt_int;
#endif

#ifdef DISLOC
void write_config_dem(int nr);
void write_config_dsp(int nr);
void reset_Epot_ref(void);
void update_ort_ref(void);
int  reset_Epot_step, calc_Epot_ref, up_ort_ref;
int  dem_int, dsp_int;
#endif

#ifdef ATDIST
void init_atdist(void);
void update_atdist(void);
void write_atdist(void);
void write_config_atdist_pos(int nr);
int  atdist_int, atdist_start, atdist_end, atdist_pos_int;
#endif

#ifdef DIFFPAT
void init_diffpat(void);
void update_diffpat(int );
void write_diffpat(void);
int  diffpat_int, diffpat_start, diffpat_end;
#endif

#ifdef AVPOS
void write_config_avpos(int nr);
void update_avpos(void);
void add_positions(void);
int  avpos_int, avpos_start, avpos_end, avpos_res;
#endif

#ifdef STRESS_TENS
void calc_tot_presstens(void);
void write_config_press(int nr);
int  do_press_calc, press_int;
#endif

#if defined(MSQD) || defined(CORRELATE)
void write_config_sqd(int nr);
void init_correl(int, int);
void correlate(int istep, int refstep, unsigned seqnum);
void write_msqd(int);
int  correl_int;
int  correl_start;
int  correl_refstep;
int  correl_end;
int  correl_ts;
int  ncorr_rmax;
int  ncorr_tmax;
#endif

#ifdef FORCE
void write_config_force(int nr);
int  force_int;
#endif

#ifdef WRITEF
void write_config_wf(int nr);
#endif

#ifdef CG
void reset_cg(void);
void cg_step(int steps);
void write_cgconfig(int steps);
#endif

#ifdef ACG
void acg_step(int steps);
real acg_alpha;
real acg_init_alpha;
#endif

#ifdef NMOLDYN
void init_nmoldyn(void);
void write_nmoldyn(int);
int  nmoldyn_int;
#endif

#ifdef HOMDEF
#ifdef TWOD
void lin_deform(vektor,vektor,real);
#else
void lin_deform(vektor,vektor,vektor,real);
#endif
void relax_pressure(void);
int  lindef_int;
real lindef_size;
vektor lindef_x;
vektor lindef_y;
#ifndef TWOD
vektor lindef_z;
#endif
#endif

#ifdef DEFORM
void deform_sample(void);
int  deform_int;
int  max_deform_int;
#endif

#ifdef FBC
void init_fbc(void);
void update_fbc();
#endif

#ifdef SOCKET_IO
void init_socket(void);
void check_socket(void);
int  socket_int;
#endif

#ifdef EWALD
void init_ewald(void);
#endif

#ifdef REFPOS
void init_refpos(void);
#endif

#ifdef RIGID
void calc_superforces(void);
#endif

#ifdef RELAX
void check_relaxed(void);
int  is_relaxed;
#endif
