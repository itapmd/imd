#ifdef MAIN
#define EXTERN
#define INIT(a) =a
#else
#define EXTERN extern
#define INIT(a)
#endif

EXTERN char sim_host[256] INIT("");      /* name of sim host */
EXTERN unsigned long varIP INIT(0);

EXTERN int *nummer;
EXTERN int columns;
EXTERN int *bcode;
EXTERN int allocated INIT(0);
EXTERN int movie_mode INIT(0);
EXTERN int bond_type INIT(1);
EXTERN int bond_mode INIT(0);
EXTERN int atom_mode INIT(1);
EXTERN int savimg_mode INIT(0);
EXTERN int col_mode INIT(1);
EXTERN int size_mode INIT(0);
EXTERN int stat_bond INIT(1);
EXTERN int scene_type INIT(0);
EXTERN int endian_byte_swap INIT(0);
EXTERN int text INIT(0);
EXTERN int eng_mode;
EXTERN int eng_minmax INIT(0);
EXTERN int qp INIT(1);
EXTERN int x_res;
EXTERN int y_res;
EXTERN int natoms;
EXTERN int nunits;
EXTERN int socket_id_int;      /* initial socket parameter for iact */

EXTERN short int *sorte;

EXTERN double *masse;
EXTERN double *x;
EXTERN double *y;
#ifndef TWOD
EXTERN double *z;
#endif
EXTERN double *vx;
EXTERN double *vy;
#ifndef TWOD
EXTERN double *vz;
#endif
EXTERN double *pot;
EXTERN double *kin;

EXTERN unsigned short base_port INIT(31913);
EXTERN unsigned short base_port_int INIT(31914);

EXTERN float temperature;
EXTERN float maxx INIT(-1000);
EXTERN float minx INIT(1000);
EXTERN float maxy INIT(-1000);
EXTERN float miny INIT(1000);
#ifndef TWOD
EXTERN float minz INIT(1000);
EXTERN float maxz INIT(-1000);
#endif
EXTERN float maxp INIT(-1000);
EXTERN float minp INIT(1000);
EXTERN float maxk INIT(-1000);
EXTERN float mink INIT(1000);

EXTERN float maxbl INIT(.9);
EXTERN float minbl INIT(1.1);

EXTERN float scalex;
EXTERN float scaley;
#ifndef TWOD
EXTERN float scalez;
#endif
EXTERN float scalepot;
EXTERN float scalekin;
EXTERN float radius INIT(.3);
EXTERN float offspot;
EXTERN float offskin;
EXTERN float *potarray;
EXTERN float *kinarray;
EXTERN float *ux;
EXTERN float *uy;
#ifndef TWOD
EXTERN float *uz;
#endif

EXTERN char *paramfilename;
EXTERN char uvfname[255];
EXTERN float lbondx[10000],lbondy[10000],lbondz[10000];
EXTERN float rbondx[10000],rbondy[10000],rbondz[10000];
EXTERN int banz;
