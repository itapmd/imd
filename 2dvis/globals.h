#ifdef MAIN
#define EXTERN
#define INIT(a) =a
#else
#define EXTERN extern
#define INIT(a)
#endif

EXTERN int *nummer;
EXTERN int columns;
EXTERN int *bcode;
EXTERN int bond_mode;
EXTERN int atom_mode INIT(1);
EXTERN int col_mode INIT(0);
EXTERN int scene_type INIT(0);
EXTERN int text INIT(0);
EXTERN int eng_mode;
EXTERN int qp;
EXTERN int radectyp INIT(0);
EXTERN int x_res;
EXTERN int y_res;
EXTERN int natoms;
EXTERN int nunits;

EXTERN short int *sorte;

EXTERN double *masse;
EXTERN double *x;
EXTERN double *y;
EXTERN double *z;
EXTERN double *vx;
EXTERN double *vy;
EXTERN double *vz;
EXTERN double *pot;
EXTERN double *kin;

EXTERN unsigned short base_port INIT(31913);

EXTERN float maxx;
EXTERN float minx;
EXTERN float maxy;
EXTERN float miny;
EXTERN float maxp;
EXTERN float minp;
EXTERN float maxk;
EXTERN float mink;
EXTERN float scalex;
EXTERN float scaley;
EXTERN float scalepot;
EXTERN float scalekin;
EXTERN float radius INIT(.3);
EXTERN float offspot;
EXTERN float offskin;
EXTERN float *potarray;
EXTERN float *kinarray;
EXTERN float *ux;
EXTERN float *uy;

EXTERN char *paramfilename;
EXTERN char uvfname[255];
