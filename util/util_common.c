
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2007 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
*  Common routines for IMD configuration utility programs
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "util.h"
#include "conf_tools.h"

/******************************************************************************
*
*  To determine the cell into which a given particle belongs, we
*  have to transform the cartesian coordinates of the particle into
*  the coordinate system spanned by the vectors of the box edges. 
*  This yields coordinates in the interval [0..1] that are
*  multiplied by cell_dim to get the cell's index.
*
******************************************************************************/

#ifdef TWOD

/* compute transformation matrix */
void make_box(void)
{
  real volume;

  /* volume */
  volume = box_x.x * box_y.y - box_x.y * box_y.x;
  if (0==volume) error("Box Edges are parallel.");

  /* compute tbox_j such that SPROD(box_i,tbox_j) == delta_ij */
  tbox_x.x =   box_y.y / volume;
  tbox_x.y = - box_y.x / volume;
  tbox_y.x = - box_x.y / volume;
  tbox_y.y =   box_x.x / volume;

  /* squares of the box heights perpendicular to the faces */
  height.x = 1.0 / SPROD(tbox_x,tbox_x);
  height.y = 1.0 / SPROD(tbox_y,tbox_y);
}

#else

/* vector product */ 
vektor vec_prod(vektor u, vektor v)
{
  vektor w;
  w.x = u.y * v.z - u.z * v.y;
  w.y = u.z * v.x - u.x * v.z;
  w.z = u.x * v.y - u.y * v.x;
  return w;
}

/* compute transformation matrix */
void make_box(void)
{
  real volume;

  /* compute tbox_j such that SPROD(box_i,tbox_j) == delta_ij */
  /* first unnormalized */
  tbox_x = vec_prod( box_y, box_z );
  tbox_y = vec_prod( box_z, box_x );
  tbox_z = vec_prod( box_x, box_y );

  /* volume */
  volume = SPROD( box_x, tbox_x );
  if (0==volume) error("Box Edges are parallel.");

#ifdef ELCO
  vol = volume;
#endif

  /* normalization */
  tbox_x.x /= volume;  tbox_x.y /= volume;  tbox_x.z /= volume;
  tbox_y.x /= volume;  tbox_y.y /= volume;  tbox_y.z /= volume;
  tbox_z.x /= volume;  tbox_z.y /= volume;  tbox_z.z /= volume;

  /* squares of the box heights perpendicular to the faces */
  height.x = 1.0 / SPROD(tbox_x,tbox_x);
  height.y = 1.0 / SPROD(tbox_y,tbox_y);
  height.z = 1.0 / SPROD(tbox_z,tbox_z);
}

#endif


/******************************************************************************
*
*  Compute the size of the cells, and initialize the cell array.
*  There must be at least two cells in each direction.
*
******************************************************************************/

void init_cells(void)
{
  vektor cell_scale;
  cell   *p;
  int    nmax, i;

  make_box();

  /* scaling factors box/cell */
  cell_scale.x = sqrt( (double)(r2_cut / height.x) );
  cell_scale.y = sqrt( (double)(r2_cut / height.y) );
#ifndef TWOD
  cell_scale.z = sqrt( (double)(r2_cut / height.z) );
#endif

#ifdef TWOD
  printf("Minimal cell size: \n\t ( %f %f ) \n\t ( %f %f ) \n",
	 box_x.x * cell_scale.x, box_x.y * cell_scale.x, 
	 box_y.x * cell_scale.y, box_y.y * cell_scale.y); 
#else
  printf("Minimal cell size: \n");
  printf("\t ( %f %f %f ) \n\t ( %f %f %f ) \n\t ( %f %f %f )\n",
     box_x.x * cell_scale.x, box_x.y * cell_scale.x, box_x.z * cell_scale.x,
     box_y.x * cell_scale.y, box_y.y * cell_scale.y, box_y.z * cell_scale.y,
     box_z.x * cell_scale.z, box_z.y * cell_scale.z, box_z.z * cell_scale.z);
#endif

  /* set up cell array */
  cell_dim.x = (int) ( 1.0 / cell_scale.x );
  cell_dim.y = (int) ( 1.0 / cell_scale.y );
#ifndef TWOD
  cell_dim.z = (int) ( 1.0 / cell_scale.z );
#endif

#ifdef TWOD
  if ((cell_dim.x==0) || (cell_dim.y==0))
#else
  if ((cell_dim.x==0) || (cell_dim.y==0) || (cell_dim.z==0))
#endif
    error("cutoff radius must be smaller than box diameter");
  
  /* If an integer number of cells does not fit exactly into the box, the
     cells are enlarged accordingly */
  cell_scale.x = 1.0 / cell_dim.x;
  cell_scale.y = 1.0 / cell_dim.y;
#ifndef TWOD
  cell_scale.z = 1.0 / cell_dim.z;
#endif

#ifdef TWOD
  printf("Actual cell size: \n\t ( %f %f ) \n\t ( %f %f ) \n",
	 box_x.x * cell_scale.x, box_x.y * cell_scale.x, 
	 box_y.x * cell_scale.y, box_y.y * cell_scale.y); 
  printf("Cell array dimensions: %d %d\n", cell_dim.x,cell_dim.y);
  cell_array = (cell *) malloc(cell_dim.x * cell_dim.y * sizeof(cell));
#else
  printf("Actual cell size: \n");
  printf("\t ( %f %f %f ) \n\t ( %f %f %f ) \n\t ( %f %f %f )\n",
     box_x.x * cell_scale.x, box_x.y * cell_scale.x, box_x.z * cell_scale.x,
     box_y.x * cell_scale.y, box_y.y * cell_scale.y, box_y.z * cell_scale.y,
     box_z.x * cell_scale.z, box_z.y * cell_scale.z, box_z.z * cell_scale.z);
  printf("Cell array dimensions: %d %d %d\n",cell_dim.x,cell_dim.y,cell_dim.z);
  cell_array = (cell *) malloc(cell_dim.x*cell_dim.y*cell_dim.z*sizeof(cell));
#endif
  if (NULL == cell_array) error("Cannot allocate memory for cells");

  nmax = cell_dim.x * cell_dim.y;
#ifndef TWOD
  nmax = nmax * cell_dim.z;
#endif

  /* initialize cells */
  for (i=0; i<nmax; ++i) {
    p = cell_array + i;
    p->n_max = 0;
    alloc_cell(p, CSTEP);
  }

}


/*****************************************************************************
*
*  cell_coord gives the (integral) cell_coorinates of a position
*
*****************************************************************************/

ivektor cell_coord(vektor ort)
{
  ivektor coord;

  /* Map positions to boxes */
  coord.x = (int) TRUNC( cell_dim.x * SPROD(ort,tbox_x) );
  coord.y = (int) TRUNC( cell_dim.y * SPROD(ort,tbox_y) );
#ifndef TWOD
  coord.z = (int) TRUNC( cell_dim.z * SPROD(ort,tbox_z) );
#endif

  /* Roundoff errors put atoms slightly out of the simulation cell */
  coord.x = coord.x <   0         ?             0 : coord.x;
  coord.x = coord.x >= cell_dim.x ? cell_dim.x -1 : coord.x;
  coord.y = coord.y <   0         ?             0 : coord.y;
  coord.y = coord.y >= cell_dim.y ? cell_dim.y -1 : coord.y;
#ifndef TWOD
  coord.z = coord.z <   0         ?             0 : coord.z;
  coord.z = coord.z >= cell_dim.z ? cell_dim.z -1 : coord.z;
#endif
  return coord;
}


/******************************************************************************
*
*  allocate memory for a cell
*
******************************************************************************/

void alloc_cell(cell *cl, int count)
{

#ifdef COVALENT
  neightab *neigh;
  int i;
#endif

  /* initialize if it is the first call */
  if (cl->n_max == 0) {
    cl->ort           = NULL;
#ifdef STRAIN
    cl->dsp           = NULL;
    cl->strain        = NULL;
    cl->strain_offdia = NULL;
    cl->empty         = NULL;
#elif defined(STRESS)
    cl->stress        = NULL;
    cl->stress_offdia = NULL;
    cl->vol           = NULL;
#else
    cl->sorte         = NULL;
#endif
#if defined(CONN) || defined(ELCO) || defined(COORD) || defined(CNA) || defined(REMAT)
    cl->nummer        = NULL;
#endif
#ifdef REMAT
    cl->masse         = NULL;
    cl->nnn           = NULL;
    cl->flag          = NULL;
#endif
#ifdef CNA
    cl->mark          = NULL;
    cl->masse         = NULL;
#endif
#ifdef COORD
    cl->coord         = NULL;
    cl->poteng        = NULL;
#endif
#ifdef COVALENT
    cl->neightab_array= NULL;
#endif
#ifdef ELCO
    cl->stress        = NULL;
    cl->elco          = NULL;
    cl->vol           = NULL;
#ifdef EAM
    cl->eam_stress    = NULL;
    cl->eam_press     = NULL;
    cl->eam_bulkm     = NULL;
    cl->eam_dbulkm    = NULL;
#endif
#endif
#ifdef EAM
    cl->eam_rho_h     = NULL;
#endif
#ifdef RING
    cl->del            = NULL;
    cl->hops           = NULL;
    cl->sp_hops        = NULL;
    cl->color          = NULL;
    cl->sp_color       = NULL;
    cl->perm_neightab_array = NULL;
#endif
#ifdef PS
    cl->enc            = NULL;
#endif
    cl->n             = 0;
  }

  /* (re)allocate */
  cl->ort    = (vektor *) realloc(cl->ort,    count * sizeof(vektor));
#ifdef STRAIN
  cl->dsp    = (vektor *) realloc(cl->dsp,    count * sizeof(vektor));
  cl->strain = (vektor *) realloc(cl->strain, count * sizeof(vektor));
  cl->strain_offdia = (vektor *) realloc(cl->strain_offdia,
                                              count * sizeof(vektor));
  cl->empty  = (short  *) realloc(cl->empty,  count * sizeof(short));
#elif defined(STRESS)
  cl->stress = (vektor *) realloc(cl->stress, count * sizeof(vektor));
  cl->stress_offdia = (vektor *) realloc(cl->stress_offdia,
                                              count * sizeof(vektor));
  cl->vol    = (real   *) realloc(cl->vol,    count * sizeof(real));
#else
  cl->sorte  = (int    *) realloc(cl->sorte,  count * sizeof(int));
#endif
#if defined(CONN) || defined(ELCO) || defined(COORD) || defined(CNA)
  cl->nummer = (int    *) realloc(cl->nummer, count * sizeof(int));
#endif
#ifdef REMAT
  cl->nummer = (int    *) realloc(cl->nummer, count * sizeof(int));
  cl->nnn = (int    *) realloc(cl->nnn, count * sizeof(int));
  cl->flag = (int    *) realloc(cl->flag, count * sizeof(int));
  cl->masse  = (real   *) realloc(cl->masse,  count * sizeof(real));
#endif
#ifdef CNA
  cl->mark   = (short  *) realloc(cl->mark,   count * sizeof(short));
  cl->masse  = (real   *) realloc(cl->masse,  count * sizeof(real));
#endif
#ifdef COORD
  cl->coord  = (real   *) realloc(cl->coord,  ntypes * count * sizeof(real));
  cl->poteng = (real   *) realloc(cl->poteng, ntypes * count * sizeof(real));
#endif
#ifdef ELCO
  cl->stress = (tensor *) realloc(cl->stress, count * sizeof(tensor));
  cl->elco   = (tensor6d *) realloc(cl->elco, count * sizeof(tensor6d));
  cl->vol    = (real   *) realloc(cl->vol,    count * sizeof(real));
#ifdef EAM
  cl->eam_stress = (tensor *) realloc(cl->eam_stress, count * sizeof(tensor));
  cl->eam_press  = (real *)   realloc(cl->eam_press, count * sizeof(real));
  cl->eam_bulkm  = (real *)   realloc(cl->eam_bulkm, count * sizeof(real));
  cl->eam_dbulkm = (real *)   realloc(cl->eam_dbulkm, count * sizeof(real));
#endif
#endif
#ifdef EAM
  cl->eam_rho_h = (real *) realloc(cl->eam_rho_h, count * sizeof(real));
#endif
#ifdef RING
  cl->del      = (int    *) realloc(cl->del,      count * sizeof(int));
  cl->hops     = (int    *) realloc(cl->hops,     count * sizeof(int));
  cl->sp_hops  = (int    *) realloc(cl->sp_hops,  count * sizeof(int));
  cl->color    = (int    *) realloc(cl->color,    count * sizeof(int));
  cl->sp_color = (int    *) realloc(cl->sp_color, count * sizeof(int));
#endif
#ifdef PS
  cl->enc     = (real   *) realloc(cl->enc,     count * sizeof(real) ); 
#endif
#ifdef COVALENT
  /* if cell shrinks, deallocate neighbor tables before shrinking the array */
  if (count < cl->n_max) {
    for (i = count; i < cl->n_max; ++i) {
      neigh = cl->neightab_array + i;
#if !defined(RING) && !defined(CNA) && !defined(PS)
      free(neigh->dist);
#endif
      free(neigh->typ);
      free(neigh->cl );
      free(neigh->num);
    }
  }
  cl->neightab_array = (neightab *) realloc( cl->neightab_array, 
					      count * sizeof(neightab));
  if (NULL == cl->neightab_array)
    error("Cannot allocate neighbor tables");

  /* Allocate memory for neighbour tables */
  for (i = cl->n_max; i < count; ++i) {
    neigh = cl->neightab_array + i;

    neigh->n     = 0;
    neigh->n_max = neigh_len;
#if !defined(RING) && !defined(CNA) && !defined(PS)
    neigh->dist  = (real *)  malloc( neigh_len * 3 * sizeof(real) );
#endif
    neigh->typ   = (short *) malloc( neigh_len * sizeof(short) );
    neigh->cl    = (void **) malloc( neigh_len * sizeof(cellptr) );
    neigh->num   = (int *)   malloc( neigh_len * sizeof(int) );

    if (
#if !defined(RING) && !defined(CNA) && !defined(PS)
      (neigh->dist==NULL) || 
#endif
      (neigh->typ ==NULL) || 
      (neigh->cl  ==NULL) || (neigh->num==NULL) )
      error("Cannot allocate memory for neighbor table");
  }
#endif
#ifdef RING
  /* if cell shrinks, deallocate neighbor tables before shrinking the array */
  if (count < cl->n_max) {
    for (i = count; i < cl->n_max; ++i) {
      neigh = cl->perm_neightab_array + i;
      free(neigh->typ);
      free(neigh->cl );
      free(neigh->num);
    }
  }
  cl->perm_neightab_array = (neightab *) realloc( cl->perm_neightab_array, 
					      count * sizeof(neightab));
  if (NULL == cl->perm_neightab_array) 
    error("Cannot allocate permanent neighbor tables");

  /* Allocate memory for permanent neighbour tables */
  for (i = cl->n_max; i < count; ++i) {
    neigh = cl->perm_neightab_array + i;

    neigh->n     = 0;
    neigh->n_max = neigh_len;
    neigh->typ   = (short *) malloc( neigh_len * sizeof(short) );
    neigh->cl    = (void **) malloc( neigh_len * sizeof(cellptr) );
    neigh->num   = (int *)   malloc( neigh_len * sizeof(int) );

    if ((neigh->typ==NULL) || (neigh->cl==NULL) || (neigh->num==NULL) )
      error("Cannot allocate memory for permanent neighbor table");
  }
#endif

  /* check if it worked */
  if ( (NULL==cl->ort)
#ifdef STRAIN
    || (NULL==cl->dsp)
    || (NULL==cl->strain)
    || (NULL==cl->strain_offdia)
    || (NULL==cl->empty)
#elif defined(STRESS)
    || (NULL==cl->stress)
    || (NULL==cl->stress_offdia)
    || (NULL==cl->vol)
#else   
    || (NULL==cl->sorte)
#endif
#if defined(CONN) || defined(ELCO) || defined(COORD) || defined(CNA)
    || (NULL==cl->nummer)
#endif
#ifdef REMAT
       || (NULL==cl->nummer)
       || (NULL==cl->nnn)
       || (NULL==cl->masse)
       || (NULL==cl->flag)
#endif
#ifdef CNA
    || (NULL==cl->mark)
    || (NULL==cl->masse)
#endif
#ifdef COORD
    || (NULL==cl->coord)
    || (NULL==cl->poteng)
#endif
#ifdef ELCO
    || (NULL==cl->stress)
    || (NULL==cl->elco)
    || (NULL==cl->vol)
#ifdef EAM
    || (NULL==cl->eam_stress)
    || (NULL==cl->eam_press)
    || (NULL==cl->eam_bulkm)
    || (NULL==cl->eam_dbulkm)
#endif
#endif
#ifdef EAM
    || (NULL==cl->eam_rho_h)
#endif
#ifdef RING
    || (NULL==cl->del)
    || (NULL==cl->hops)
    || (NULL==cl->sp_hops)
    || (NULL==cl->color)
    || (NULL==cl->sp_color)
#endif
#ifdef PS
   || (NULL==cl->enc)
#endif
  ) error("Cannot allocate memory for cell.");
  cl->n_max  = count;
}

/******************************************************************************
*
* read_box reads the box from the config file
*
******************************************************************************/

void read_box(str255 infilename)
{
  FILE   *infile;
  str255 line, fname;

  infile = fopen(infilename,"r");
  if (NULL==infile) error_str("cannot open input file %s", infilename);
  fgets(line, 255, infile);
  while ((line[0]=='#') && (line[1]!='E')) {
#ifdef TWOD
    if      (line[1]=='X') 
      sscanf(line+2, "%lf %lf", &box_x.x, &box_x.y);
    else if (line[1]=='Y') 
      sscanf(line+2, "%lf %lf", &box_y.x, &box_y.y);
#else
    if      (line[1]=='X') 
      sscanf(line+2, "%lf %lf %lf", &box_x.x, &box_x.y, &box_x.z);
    else if (line[1]=='Y') 
      sscanf(line+2, "%lf %lf %lf", &box_y.x, &box_y.y, &box_y.z);
    else if (line[1]=='Z') 
      sscanf(line+2, "%lf %lf %lf", &box_z.x, &box_z.y, &box_z.z);
#endif
    fgets(line, 255, infile);
    if (feof(infile)) break;
  }
  fclose(infile);
}

/******************************************************************************
*
*  read_atoms - reads atoms into the cell-array
*
******************************************************************************/

void read_atoms(str255 infilename)
{
  FILE *infile;
  char buf[512];
  int  p, s, n, i, have_header;
  real m, e, c, tmp;
  cell *to;
  ivektor cellc;
  header_info_t info;
  atom_t atom;

  /* we first try the old checkpoint name, then the new */
  infile = fopen(infilename,"r");
  if ((NULL == infile) && (restart != -1)) {
    infilename = strcat(infilename,".chkpt");
    infile = fopen(infilename,"r");
  }
  if (NULL==infile) error_str("Cannot open atoms file %s", infilename);

  /* read file header */
  have_header = read_header(&info, infilename);
  natoms=0;

  /* eat header, if there is one */
  if (have_header) {
    fgets(buf,sizeof(buf),infile);
    while (('#'==buf[0]) && ('E'!=buf[1]) && !feof(infile)) {
      fgets(buf,sizeof(buf),infile);
    }
  }
 
  /* read the config file */
  while (!feof(infile)) {

    p = read_atom(&info, infile, &atom);
    if (info.n_number == 0) atom.number = natoms;

    if (p>0) {

   
#ifdef PS
      /* Determine bounding box before rotation */
      real_minl = minvektor(real_minl, atom.pos);
      real_maxl = maxvektor(real_maxl, atom.pos);

      /* Rotations of atom positions */
#ifndef TWOD
      if ( angx != 0.0 ) 
	atom.pos = xrotate(atom.pos, angx);
      if ( angy != 0.0 ) 
	atom.pos = yrotate(atom.pos, angy);
#endif
      if ( angz != 0.0 ) 
	atom.pos = zrotate(atom.pos, angz);
#endif

      /* compute target cell */
      cellc = cell_coord(atom.pos);
      to = PTR_VV(cell_array,cellc,cell_dim);
      /* enlarge it if necessary */
      if (to->n >= to->n_max) alloc_cell(to,to->n_max+CSTEP);
   
      /* put the data */

      to->ort   [to->n] = atom.pos;
     

#if (!defined(STRAIN) && !defined(STRESS))
      to->sorte [to->n] = MOD(atom.type,ntypes);
#endif

#if defined(CONN) || defined(ELCO) || defined(COORD) || defined(CNA)
      to->nummer[to->n] = atom.number;
#endif

#ifdef REMAT
   
      to->nummer[to->n] = atom.number;
      to->masse[to->n] = atom.mass;
      to->nnn[to->n] = 0;
      to->flag[to->n] = 0;
   
#endif
#ifdef CNA
      to->mark[to->n]  = 0;
      to->masse[to->n] = atom.mass;
#endif
#ifdef CONN
      n_min = MIN(n_min,atom.number);
      n_max = MAX(n_max,atom.number);
#endif
#ifdef COORD
      /* Initialization */
      for ( i=0; i<ntypes; i++) 
	to->coord[to->n * ntypes + i] = 0.0;
      /* store potential energy */
      to->poteng[to->n] = atom.data[0];
#endif
#ifdef RING
      to->del[to->n] = 0;
#endif

#ifdef PS
      /* Determine bounding box after rotation */
      minl = minvektor( minl, atom.pos);
      maxl = maxvektor( maxl, atom.pos);
      if ( p > 6 && ( eng != 0 || kineng != 0 || color_enc != 0 ) ) {
        if ( eng != 0 )
          tmp = atom.data[0];
        else if ( kineng != 0 ) 
          tmp = 0.5 * atom.mass * SPROD(atom.vel,atom.vel);
        else if ( color_enc != 0 )
          tmp = atom.data[1];
        maxeng = MAX(maxeng, tmp);
        mineng = MIN(mineng, tmp);
        to->enc[to->n] = tmp;
      }
#endif

      to->n++;
      natoms++;
    }
  }
  fclose(infile);  
}

/******************************************************************************
*
*  Loop over all cell pairs
*
******************************************************************************/

void do_work(void (*do_cell_pair)(cell *p, cell *q, vektor pbc))
{
  cell *p,*q;
  int i,j,k;
  int l,m,n;
  int r,s,t;
  vektor pbc;

  /* for each cell */
  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
#ifndef TWOD
      for (k=0; k < cell_dim.z; ++k)
#endif
      {
#ifdef TWOD
        p = PTR_2D_V(cell_array,i,j  ,cell_dim);
#else
        p = PTR_3D_V(cell_array,i,j,k,cell_dim);
#endif
        /* For half of the neighbours of this cell */
        for (l=0; l <= 1; ++l)
          for (m=-l; m <= 1; ++m)
#ifndef TWOD
            for (n=(l==0 ? -m  : -l ); n <= 1; ++n)
#endif
            {
              /* Calculate Indicies of Neighbour */
              r = i+l;  pbc.x = 0;
              s = j+m;  pbc.y = 0;
#ifndef TWOD
              t = k+n;  pbc.z = 0;
#endif

              /* deal with periodic boundary conditions if necessary */
              if (r<0) {
                if (pbc_dirs.x==1) {
                  r = cell_dim.x-1; 
                  pbc.x -= box_x.x;      
                  pbc.y -= box_x.y;
#ifndef TWOD
                  pbc.z -= box_x.z;
#endif
                } else continue;
              }
              if (s<0) {
                if (pbc_dirs.y==1) {
                  s = cell_dim.y-1;
                  pbc.x -= box_y.x;      
                  pbc.y -= box_y.y;
#ifndef TWOD
                  pbc.z -= box_y.z;
#endif
                } else continue;
              }
#ifndef TWOD
              if (t<0) {
                if (pbc_dirs.z==1) {
                  t = cell_dim.z-1;
                  pbc.x -= box_z.x;      
                  pbc.y -= box_z.y;
                  pbc.z -= box_z.z;
                } else continue;
              }
#endif
              if (r>cell_dim.x-1) {
                if (pbc_dirs.x==1) {
                  r = 0; 
                  pbc.x += box_x.x;      
                  pbc.y += box_x.y;
#ifndef TWOD
                  pbc.z += box_x.z;
#endif
                } else continue;
              }
              if (s>cell_dim.y-1) {
                if (pbc_dirs.y==1) {
                  s = 0; 
                  pbc.x += box_y.x;      
                  pbc.y += box_y.y;
#ifndef TWOD
                  pbc.z += box_y.z;
#endif
                } else continue;
              }
#ifndef TWOD
              if (t>cell_dim.z-1) {
                if (pbc_dirs.z==1) {
                  t = 0; 
                  pbc.x += box_z.x;      
                  pbc.y += box_z.y;
                  pbc.z += box_z.z;
                } else continue;
              }
#endif

	      /* Neighbour cell (note that p==q ist possible) */
#ifdef TWOD
              q = PTR_2D_V(cell_array,r,s,cell_dim);
#else
              q = PTR_3D_V(cell_array,r,s,t,cell_dim);
#endif
              /* Do the work */
              do_cell_pair(p,q,pbc);
            }
      }
}

#ifdef CNA

/******************************************************************************
*
*  atom_in_pbox --  returns 1 if atom is in partial box defined by ll
*                   and ur, 0 otherwise
*
******************************************************************************/

int atom_in_pbox(vektor pos)
{
  if ( pos.x > ll.x && pos.y > ll.y && pos.x < ur.x && pos.y < ur.y
#ifndef TWOD
      && pos.z > ll.z && pos.z < ur.z 
#endif
      )
    return 1;
  else
    return 0;
}

#endif
