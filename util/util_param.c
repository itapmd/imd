
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
*  Routines for reading of parameters from command line and parameter file
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "util.h"

/******************************************************************************
*
*  Read the command line parameters
*
******************************************************************************/

void read_command_line(int argc, char **argv)
{
  strcpy(progname,argv[0]);

  while ((argc > 1)) {

    /* All parameters start with '-' */
    if ( argv[1][0] != '-' ) {
      printf("Illegal option %s\n", argv[1]);
      usage();
      exit(-1);
    }

    switch (argv[1][1]) {

      /* r - restart */
    case 'r':
      read_arg_int(&argc, &argv, &restart);
      break;

      /* A - use avpos outfiles */
    case 'A':
      read_arg_int(&argc, &argv, &avpos);
      break;

    case 'p':
      read_arg_string(&argc, &argv, &paramfilename);
      break;
      
    case 'v':
      read_arg_bool(&argc, &argv, &use_vtypes);
      break;

#if defined(ANGLE) || defined(PAIR) || defined(COORD) || defined(CNA) ||defined(REMAT)
      /* a - minimum radius */
    case 'a':
      read_arg_real(&argc, &argv, &r_min);
      break;
#endif

#if defined(ANGLE) || defined(PAIR) || defined(TORSION) || defined(CNA)
      /* n - number of slots in histogram */
    case 'n':
      read_arg_int(&argc, &argv, &slots);
      break;
#endif

#ifdef CNA
      /* f - consider only nearest neighbour atoms */
    case 'f':
      read_arg_bool(&argc, &argv, &nearestneighbour);
      break;
      /* g - write radial distribution function */
    case 'g':
      read_arg_bool(&argc, &argv, &rdf);
      break;
     /* l - lower left coordinates of partial box */
    case 'l':
      read_arg_vektor(&argc, &argv, &ll);
      use_ll = 1;
      break;
      /* u - upper right coordinates of partial box */
    case 'u':
      read_arg_vektor(&argc, &argv, &ur);
      use_ur = 1;
      break;
    case 'w':
      read_arg_int(&argc, &argv, &pair_type_short);
      writeatoms = 1;
      break;
      /* W - pair type to be written out */
    case 'W':
      read_arg_ivektor4d(&argc, &argv, &pair_type);
      writeatoms = 1;
      break;
#endif

#ifdef CONN
      /* n - number of maximal neighbours */
    case 'n':
      read_arg_int(&argc, &argv, &maxneigh);
      cm_dim.x = 0;
      cm_dim.y = maxneigh;
      break;
#endif

#ifndef CONN
      /* e - maximum radius */
    case 'e':
      read_arg_real(&argc, &argv, &r_max);
#ifdef STRAIN
      r_cell = r_max;
#endif
      break;
#endif

#ifdef COORD
      /* u - use cutoff-function equal to 1 */
    case 'u':
      read_arg_bool(&argc, &argv, &use_unity);
      break;

      /* l - write local coordination numbers */
    case 'l':
      read_arg_bool(&argc, &argv, &local);
      break;

      /* g - write global coordination numbers */
    case 'g':
      read_arg_bool(&argc, &argv, &global);
      break; 

      /* C - largest coordination number */
    case 'C':
      read_arg_int(&argc, &argv, &c_max);
      break;

      /* E - write potential energy */
    case 'E':
      read_arg_bool(&argc, &argv, &write_poteng);
      break; 
#endif
 
#if defined(STRAIN) || defined(ELCO)
      /* c - cell width */
    case 'c':
      read_arg_real(&argc, &argv, &r_cell);
      break;
#endif

#ifdef ELCO
      /* s - Ouput local stress tensor */
    case 's':
      read_arg_bool(&argc, &argv, &stresstens);
      break;

      /* m - Output local elastic moduli c11, c12, c44 */
    case 'm':
      read_arg_bool(&argc, &argv, &moduli);
      break;

      /* M - Output all local elastic moduli */
    case 'M':
      read_arg_bool(&argc, &argv, &all_moduli);
      break;
      /* w - Output local elastic constant C_ij */
    case 'w':
      read_arg_int(&argc, &argv, &cindex);
      select_moduli = 1;
      break;
#endif

#ifdef RING
      /* l - Maximal length of rings */
    case 'l':
      read_arg_int(&argc, &argv, &max_length);
      break;

      /* n - Size of neighbour tables */
    case 'n':
      read_arg_int(&argc, &argv, &neigh_len);
      break;
#endif

#ifdef PS
      /* a - draw coordinate axes */
    case 'a':
      read_arg_int(&argc, &argv, &axes);
      break;

      /* b - brightness of bonds */
    case 'b':
      read_arg_real(&argc, &argv, &bondbrightness);
      break;

      /* B - bond radius */
    case 'B':
      read_arg_real(&argc, &argv, &r_bond);
      break;

      /* c - crop boundary strip */
    case 'c':
      read_arg_real(&argc, &argv, &crop);
      break;

      /* C - atom colors */
    case 'C':
      read_arg_long(&argc, &argv, &atomcolor);
      break;

#ifndef TWOD
      /* d - shading with depth */
    case 'd':
      read_arg_real(&argc, &argv, &depth);
      break;
#endif

      /* D - atom radii */
    case 'D':
      read_arg_long(&argc, &argv, &radii);
      break; 

      /* E - Use color encoding for energy */
    case 'E':
      read_arg_long(&argc, &argv, &eng);
      break;

      /* f - draw boxes */
    case 'f':
      read_arg_real(&argc, &argv, &frame);
      break;

      /* F - linewidth of bounding box */
    case 'F':
      read_arg_real(&argc, &argv, &frame_linew);
      break; 

      /* g - background color */
    case 'g':
      read_arg_real(&argc, &argv, &backgrd);
      break;

       /* G - colorencoding of column G */ 
    case 'G':
      read_arg_long(&argc, &argv, &color_enc);
      break;

      /* h - height of picture */
    case 'h':
      read_arg_real(&argc, &argv, &ps_height);
      break;

      /* H - linewidth of atoms, bonds, and bond separators */
    case 'H':
      read_arg_real(&argc, &argv, &linew);
      break;

      /* i - x translation */
    case 'i':
      read_arg_real(&argc, &argv, &translation.x);
      break;

      /* I - Brigthness of atom shading */
    case 'I':
      read_arg_real(&argc, &argv, &ambient);
      break;

      /* j - y translation */
    case 'j':
      read_arg_real(&argc, &argv, &translation.y);
      break;

      /* J - brightness of atoms */
    case 'J':
      read_arg_real(&argc, &argv, &atombrightness);
      break;

     /* k - show color encoding of energy */
    case 'k':
      read_arg_real(&argc, &argv, &spect);
      break;

     /* K write file with parameter settings */
    case 'K':
      read_arg_bool(&argc, &argv, &settings);
      break;

     /* l linewidth of bond borders */
    case 'l':
      read_arg_real(&argc, &argv, &blinew);
      break;

     /* L - line width of atom borders */
    case 'L':
      read_arg_real(&argc, &argv, &alinew);
      break;

      /* m - zoom */
    case 'm':
      read_arg_real(&argc, &argv, &zoom);
      break;

      /* M - Scale factor with depth */
    case 'M':
      read_arg_real(&argc, &argv, &user_scale);
      break;

#ifndef TWOD
      /* n - draw lightshadow */
    case 'n':
      read_arg_real(&argc, &argv, &llinew);
      break;

     /* N - bonds like needles */
    case 'N':
      read_arg_real(&argc, &argv, &needle);
      break;
#endif

     /* o - draw simulation box */
    case 'o':
      read_arg_real(&argc, &argv, &sframe);
      break;

     /* O - linewidth of simulation box */
    case 'O':
      read_arg_real(&argc, &argv, &sframe_linew);
      break;

     /* P - linewidth of bond separators */
    case 'P':
      read_arg_real(&argc, &argv, &slinew);
      break;

     /* q - color of border of frames */
    case 'q':
      read_arg_real(&argc, &argv, &frame_border);
      break;

     /* Q - linewidth of border of frames */
    case 'Q':
      read_arg_real(&argc, &argv, &frame_border_linew);
      break;

      /* R - atom radius */
    case 'R':
      read_arg_real(&argc, &argv, &r_atom);
      break;

      /* s - scale factor */
    case 's':
      read_arg_real(&argc, &argv, &scaling);
      break; 

      /* S - shading of atoms */
    case 'S':
      read_arg_real(&argc, &argv, &shading);
      break; 

      /* t - kinetic energy */
    case 't':
      read_arg_long(&argc, &argv, &kineng);
      break; 

#ifndef TWOD
      /* T - Location of projection plane */
    case 'T':
      read_arg_real(&argc, &argv, &zeta);
      break;
#endif

      /* u - draw real bounding box */
    case 'u':
      read_arg_real(&argc, &argv, &frame);
      break;

      /* U - draw real bounding box */
    case 'U':
      read_arg_real(&argc, &argv, &rframe_linew);
      break;

      /* V - color of bonds */
    case 'V':
      read_arg_long(&argc, &argv, &bondcolor);
      break;

      /* w - width of picture */
    case 'w':
      read_arg_real(&argc, &argv, &ps_width);
      break;

     /* W - draw wireframe */
    case 'W':
      read_arg_real(&argc, &argv, &wireframe);
      break;

#ifndef TWOD 
      /* x - Rotation */
    case 'x':
      read_arg_real(&argc, &argv, &angx);
      break;

      /* y - Rotation */
    case 'y':
      read_arg_real(&argc, &argv, &angy);
      break;
#endif

      /* z - Rotation */
    case 'z':
      read_arg_real(&argc, &argv, &angz);
      break;

#ifndef TWOD
      /* X - Center of projection */
    case 'X':
      read_arg_real(&argc, &argv, &foc.x);
      break;

     /* Y - Center of projection */
    case 'Y':
      read_arg_real(&argc, &argv, &foc.y);
      break;

      /* Z - Center of projection */
    case 'Z':
      read_arg_real(&argc, &argv, &foc.z);
      break;
#endif
#endif

    default:
      printf("Illegal option %s \n",argv[1]);
      usage();
      exit(-1);
    }
  }
}

/*****************************************************************************
*
*  Subroutines for reading of command line arguments
*
*****************************************************************************/

void read_arg_bool(int *argcptr, char ***argvptr, int *parptr)
{
  int n;

  *parptr = 1;

  if ( (*argvptr)[1][2] != '\0' ) {
    /* remove leading boolean parameter */
    for( n=1; ( (*argvptr)[1][n] = (*argvptr)[1][n+1] ) != '\0'; n++)
      ;
  }
  else {
    *argcptr -= 1;
    *argvptr += 1;
  }
}

void read_arg_int(int *argcptr, char ***argvptr, int *parptr) 
{   
  if ( (*argvptr)[1][2]=='\0' && NULL != (*argvptr)[2] ) {
    *parptr = atoi( (*argvptr)[2] );
    *argcptr -= 2;
    *argvptr += 2;
  }
  else { 
    *parptr = atoi( &(*argvptr)[1][2] );  
    *argcptr -= 1;
    *argvptr += 1;
  }
}

void read_arg_long(int *argcptr, char ***argvptr, long *parptr) 
{   
  if ( (*argvptr)[1][2]=='\0' && NULL != (*argvptr)[2] ) {
    *parptr = atol( (*argvptr)[2] );
    *argcptr -= 2;
    *argvptr += 2;
  }
  else { 
    *parptr = atol( &(*argvptr)[1][2] );  
    *argcptr -= 1;
    *argvptr += 1;
  }
}

void read_arg_real(int *argcptr, char ***argvptr, real *parptr)
{
  if ( (*argvptr)[1][2]=='\0' && NULL != (*argvptr)[2] ) {
    *parptr = atof( (*argvptr)[2] );
    *argcptr -= 2;
    *argvptr += 2;
  }
  else { 
    *parptr = atof( &(*argvptr)[1][2] );
    *argcptr -= 1;
    *argvptr += 1;
  }
}

void read_arg_string(int *argcptr, char ***argvptr, char **parptr)
{
  if ( (*argvptr)[1][2]=='\0' && NULL != (*argvptr)[2] ) {
    *parptr = strdup( (*argvptr)[2] );
    *argcptr -= 2;
    *argvptr += 2;
  }
  else { 
    *parptr = strdup( &(*argvptr)[1][2] );
    *argcptr -= 1;
    *argvptr += 1;
  }
}

void read_arg_vektor(int *argcptr, char ***argvptr, vektor *parptr)
{
  if ( (*argvptr)[1][2]=='\0' && NULL != (*argvptr)[2] && NULL != (*argvptr)[3]
#ifndef TWOD
      && NULL != (*argvptr)[4]
#endif
      ) {
    parptr->x = atof( (*argvptr)[2] );
    parptr->y = atof( (*argvptr)[3] );
#ifndef TWOD
    parptr->z = atof( (*argvptr)[4] );
#endif
    *argcptr -= DIM + 1;
    *argvptr += DIM + 1;
  }
  else if ( NULL != (*argvptr)[2]
#ifndef TWOD
	   && NULL != (*argvptr)[3]
#endif
	   ) {
    parptr->x = atof( &(*argvptr)[1][2] );
    parptr->y = atof( (*argvptr)[2] );
#ifndef TWOD
    parptr->z = atof( (*argvptr)[3] );
#endif
    *argcptr -= DIM;
    *argvptr += DIM;
  }
}

#ifdef CNA
void read_arg_ivektor4d(int *argcptr, char ***argvptr, ivektor4d *parptr)
{
  if ( (*argvptr)[1][2]=='\0' && NULL != (*argvptr)[2] && NULL != (*argvptr)[3]
       && NULL != (*argvptr)[4]  && NULL != (*argvptr)[5] ) {

    parptr->i = atoi( (*argvptr)[2] );
    parptr->x = atoi( (*argvptr)[3] );
    parptr->y = atoi( (*argvptr)[4] );
    parptr->z = atoi( (*argvptr)[5] );

    *argcptr -= 5;
    *argvptr += 5;
  }
  else if ( NULL != (*argvptr)[2] && NULL != (*argvptr)[3] 
	    && NULL != (*argvptr)[4] ) {

    parptr->i = atof( &(*argvptr)[1][2] );
    parptr->x = atof( (*argvptr)[2] );
    parptr->y = atof( (*argvptr)[3] );
    parptr->z = atof( (*argvptr)[4] );

    *argcptr -= 4;
    *argvptr += 4;
  }
}
#endif

/*****************************************************************************
*
*  Read parameters from parameter file
*
*****************************************************************************/

void read_parameters(void)
{
  str255 fname;
  FILE *testfile;

  getparamfile(paramfilename);
  if (use_vtypes==1) ntypes = vtypes;

  /* Get restart parameters if restart */
  if (-1 != restart) {
    sprintf(fname,"%s.%d.itr",outfilename,restart);
    testfile = fopen(fname,"r");
    if (NULL==testfile) { 
      sprintf(fname,"%s.%05d.itr",outfilename,restart);
    } else {
      fclose(testfile);
    }
    getparamfile(fname);
  }
  else if  (-1 != avpos) {
    sprintf(fname,"%s.%d.avp.itr",outfilename,avpos);
    testfile = fopen(fname,"r");
    if (NULL==testfile) { 
      sprintf(fname,"%s.%05d.avp.itr",outfilename,avpos);
    } else {
      fclose(testfile);
    }
    getparamfile(fname);
  }
#ifdef STRAIN
  if (-1 == restart) {
    sprintf(infilename, "%s.dsp", outfilename);
  }
  else {
    sprintf(infilename,"%s.%u.dsp",outfilename,restart);
    testfile = fopen(infilename,"r");
    if (NULL==testfile) { 
      sprintf(infilename,"%s.%05d.dsp",outfilename,restart);
    } else {
      fclose(testfile);
    }
  }
#elif defined(STRESS)
  if (-1 == restart) {
    sprintf(infilename, "%s.press",outfilename);
  }
  else {
    sprintf(infilename,"%s.%u.press",outfilename,restart);
    testfile = fopen(infilename,"r");
    if (NULL==testfile) { 
      sprintf(infilename,"%s.%05d.press",outfilename,restart);
    } else {
      fclose(testfile);
    }
  }
#else
  if (-1 != restart) {
    sprintf(infilename,"%s.%u.chkpt",outfilename,restart);
    testfile = fopen(infilename,"r");
    if (NULL==testfile) { 
      sprintf(infilename,"%s.%05d.chkpt",outfilename,restart);
    } else {
      fclose(testfile);
    }
  }
  else if (-1 != avpos) {
    sprintf(infilename,"%s.%u.avp",outfilename,avpos);
    testfile = fopen(infilename,"r");
    if (NULL==testfile) { 
      sprintf(infilename,"%s.%05d.avp",outfilename,avpos);
    } else {
      fclose(testfile);
    }
  }
#endif

#ifdef STRAIN
  /* r_cell >= r_max is required */
  if (r_cell < r_max) error("Cell smaller than cutoff radius!");
  r2_cut = SQR(r_cell);
#endif
#if defined(STRESS) || defined (RING)
  r2_cut = SQR(r_max);  
#endif

} 

/*****************************************************************************
*
*  Read tag-based parameter file
*  Lines beginning with comment characters '#' or blank lines are skipped
*
*****************************************************************************/

void getparamfile(char *paramfname)
{
  FILE *pf;
  char buffer[1024];
  char *token;
  char *res;

  pf = fopen(paramfname,"r");
  if (NULL == pf) error_str("Cannot open parameter file %s", paramfname);

  do {
    res=fgets(buffer,1024,pf);
    if (NULL == res) break; /* probably EOF reached */
    token = strtok(buffer," =\t\r\n");
    if (NULL == token) continue; /* skip blank lines */
    if (token[0]=='#') continue; /* skip comments */

    if (strcasecmp(token,"coordname")==0) {
      /* file name for atom coordinate input data */
      getparam("coordname",infilename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"outfiles")==0) {
      /* output file basename */
      getparam("outfiles",outfilename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"box_x")==0) {
      /* 'x' or first vector for box */
      getparam("box_x",&box_x,PARAM_REAL,DIM,DIM);
    }
    else if (strcasecmp(token,"box_from_header")==0) {
      /* read box from config file header */
      getparam(token,&box_from_header,PARAM_INT,1,1);
    }
    else if (strcasecmp(token,"box_y")==0) {
      /* 'y' or second vector for box */
      getparam("box_y",&box_y,PARAM_REAL,DIM,DIM);
    }
#ifndef TWOD
    else if (strcasecmp(token,"box_z")==0) {
      /* 'z' or third vector for box */
      getparam("box_z",&box_z,PARAM_REAL,DIM,DIM);
    }
#endif
#ifndef PS
    else if (strcasecmp(token,"pbc_dirs")==0) {
      /* directions with periodic boundary conditions */
      getparam("pbc_dirs",&pbc_dirs,PARAM_INT,DIM,DIM);
    }
#endif
    else if (strcasecmp(token,"ntypes")==0) {
      /* number of atom types */
      getparam("ntypes",&ntypes,PARAM_INT,1,1);
      vtypes = MAX(vtypes,ntypes);
    }
    else if (strcasecmp(token,"total_types")==0) {
      /* number of virtual atom types */
      getparam("total_types",&vtypes,PARAM_INT,1,1);
    }
#ifdef CONN
    else if (strcasecmp(token,"r_cut")==0) {
      /* cutoff radii */
      int nn;
      if (use_vtypes) nn = SQR(vtypes);
      else            nn = SQR(ntypes);
      r_cut = (real *) calloc(nn,sizeof(real));
      if (NULL == r_cut) error("cannot allocate r_cut");
      getparam("r_cut",r_cut,PARAM_REAL,nn,nn);
    }
#endif
#ifdef REMAT
 else if (strcasecmp(token,"r_crit")==0) {
      getparam("r_crit",&r_crit,PARAM_REAL,1,1);
    }
#endif
#ifdef PAIR_POT
    else if (strcasecmp(token,"potfile")==0) {
      getparam("potfile",potfilename,PARAM_STR,1,255);
    }
#endif
#ifdef COVALENT
    else if (strcasecmp(token,"neigh_len")==0) {
      /* number of neighbors */
      getparam("neigh_len",&neigh_len,PARAM_INT,1,1);
    }
#endif
#if defined(TERSOFF) || defined(PS)
    else if (strcasecmp(token,"ters_r_cut")==0) {     
      getparam("ters_r_cut",ters_r_cut,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
#endif
#ifdef TERSOFF
    else if (strcasecmp(token,"ters_r0")==0) {     
      getparam("ters_r0",ters_r0,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"ters_a")==0) {     
      getparam("ters_a",ters_a,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"ters_b")==0) {     
      getparam("ters_b",ters_b,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"ters_la")==0) {     
      getparam("ters_la",ters_la,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"ters_mu")==0) {     
      getparam("ters_mu",ters_mu,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"ters_chi")==0) {     
      getparam("ters_chi",ters_chi,PARAM_REAL,ntypes*(ntypes-1)/2,ntypes*(ntypes-1)/2);
    }
    else if (strcasecmp(token,"ters_om")==0) {     
      getparam("ters_om",ters_om,PARAM_REAL,ntypes*(ntypes-1)/2,ntypes*(ntypes-1)/2);
    }
    else if (strcasecmp(token,"ters_ga")==0) {     
      getparam("ters_ga",ters_ga,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"ters_n")==0) {     
      getparam("ters_n",ters_n,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"ters_c")==0) {     
      getparam("ters_c",ters_c,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"ters_d")==0) {     
      getparam("ters_d",ters_d,PARAM_REAL,ntypes,ntypes);
    }
    else if (strcasecmp(token,"ters_h")==0) {     
      getparam("ters_h",ters_h,PARAM_REAL,ntypes,ntypes);
    }
#elif STIWEB
    else if (strcasecmp(token,"stiweb_a")==0) {     
      getparam("stiweb_a",stiweb_a,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"stiweb_b")==0) {     
      getparam("stiweb_b",stiweb_b,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"stiweb_p")==0) {     
      getparam("stiweb_p",stiweb_p,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"stiweb_q")==0) {     
      getparam("stiweb_q",stiweb_q,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"stiweb_a1")==0) {     
      getparam("stiweb_a1",stiweb_a1,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"stiweb_de")==0) {     
      getparam("stiweb_de",stiweb_de,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"stiweb_a2")==0) {     
      getparam("stiweb_a2",stiweb_a2,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"stiweb_ga")==0) {     
      getparam("stiweb_ga",stiweb_ga,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"stiweb_la")==0) {     
      getparam("stiweb_la",stiweb_la,PARAM_REAL,ntypes*ntypes*(ntypes+1)/2,ntypes*ntypes*(ntypes+1)/2);
    }
#elif KEATING
    else if (strcasecmp(token,"keating_alpha")==0) {     
      getparam("keating_alpha",keating_alpha,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"keating_d")==0) {     
      getparam("keating_d",keating_d,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"keating_r_cut")==0) {     
      getparam("keating_r_cut",keating_r_cut,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
    else if (strcasecmp(token,"keating_beta")==0) {     
      getparam("keating_beta",keating_beta,PARAM_REAL,ntypes*ntypes*(ntypes+1)/2,ntypes*ntypes*(ntypes+1)/2);
    }
#endif
#ifdef EAM
    else if (strcasecmp(token,"core_potential_file")==0) {
      /* EAM:Filename for the tabulated Core-Core Potential (r^2) */
      getparam("core_potential_file",potfilename,PARAM_STR,1,255);
    }
    else if (strcasecmp(token,"embedding_energy_file")==0) {
      /* EAM:Filename for the tabulated Embedding Enery(rho_h) */
      getparam("embedding_energy_file",eam_emb_E_filename,PARAM_STR,1,255);
    }
   else if (strcasecmp(token,"atomic_e-density_file")==0) {
      /* EAM:Filename for the tabulated atomic electron density(r_ij^2) */
      getparam("atomic_e-density_file",eam_at_rho_filename,PARAM_STR,1,255);
    }
#endif
#if defined(COORD) || defined(PS) || defined(CNA)
    else if (strcasecmp(token,"r_cut")==0) {
      /* Cutoff parameter */
      getparam("r_cut",&r_cut_vec,PARAM_REAL,ntypes*(ntypes+1)/2,ntypes*(ntypes+1)/2);
    }
#endif
  } while (!feof(pf));
  fclose(pf);

#ifdef STRAIN
  if ( box_x.x < r_cell || box_y.y < r_cell 
#ifndef TWOD
       || box_z.z < r_cell 
#endif
       )
    error("Cell size greater than box!");
#endif

} /* getparamfile */

