
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
* imd_geom_3d.c -- domain decomposition routines, 3d version
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#define INDEXED_ACCESS
#include "imd.h"

/******************************************************************************
*
*  To determine the cell into which a given particle belongs, we
*  have to transform the cartesian coordinates of the particle into
*  the coordinate system spanned by the vectors of the box edges. 
*  This yields coordinates in the interval [0..1] that are
*  multiplied by global_cell_dim to get the cell's index.
*
******************************************************************************/

/* vector product */ 
vektor vec_prod(vektor u, vektor v)
{
  vektor w;
  w.x = u.y * v.z - u.z * v.y;
  w.y = u.z * v.x - u.x * v.z;
  w.z = u.x * v.y - u.y * v.x;
  return w;
}

/******************************************************************************
*
*  compute box transformation matrix;
*  initialize or revise the cell array if nessessary
*
******************************************************************************/

void make_box( void )
{

#ifdef debugLo
    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("print "" make_box start "" checking by Lo! \n");fflush(stdout);
#endif
  /* compute tbox_j such that SPROD(box_i,tbox_j) == delta_ij */
  /* first unnormalized */
  tbox_x = vec_prod( box_y, box_z );
  tbox_y = vec_prod( box_z, box_x );
  tbox_z = vec_prod( box_x, box_y );

  /* volume */
  volume = SPROD( box_x, tbox_x );
  if ((0==myid) && (0==volume)) error("Box Edges are parallel.");

  /* normalization */
  tbox_x.x /= volume;  tbox_x.y /= volume;  tbox_x.z /= volume;
  tbox_y.x /= volume;  tbox_y.y /= volume;  tbox_y.z /= volume;
  tbox_z.x /= volume;  tbox_z.y /= volume;  tbox_z.z /= volume;

  //  printf("tbox_x: %lf %lf %lf tbox_y: %lf %lf %lf tbox_z: %lf %lf %lf \n",tbox_x.x,tbox_x.y,tbox_x.z,tbox_y.x,tbox_y.y,tbox_y.z,tbox_z.x,tbox_z.y,tbox_z.z); fflush(stdout);

  /* squares of the box heights perpendicular to the faces */
  height.x = 1.0 / SPROD(tbox_x,tbox_x);
  height.y = 1.0 / SPROD(tbox_y,tbox_y);
  height.z = 1.0 / SPROD(tbox_z,tbox_z);

  /* initialize or revise cell division if necessary */
  if ( (height.x < min_height.x) || (height.x > max_height.x)
    || (height.y < min_height.y) || (height.y > max_height.y)
    || (height.z < min_height.z) || (height.z > max_height.z)
  ) init_cells();

  /* do some sanity checks */
  if (0 > volume) {
    volume = -volume;
    warning("System of box vectors is left-handed!");
  }
  if (volume_init==0) {
    volume_init = volume;
  } else {
    if ((myid==0) && (volume>8*volume_init)) error("system seems to explode!");
  }

#ifdef debugLo
    printf("print "" make_box end "" checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("    ************************* \n");fflush(stdout);
#endif
}


/******************************************************************************
*
*  Compute the size of the cells, and initialize the cell array.
*
******************************************************************************/

void init_cells( void )
{
  int i, j, k, l;
  real tmp, tol=1.0;
  vektor cell_scale;
  ivektor next_cell_dim, cell_dim_old, cd, cellc;
  minicell *p, *cell_array_old, *to;
  str255 msg;

#ifdef NBLIST
  /* add neighbor list margin (only the first time) */
  if (NULL == cell_array)
    cellsz = SQR( sqrt((double) cellsz) + nbl_margin );
#endif

#ifdef NPT
  /* if NPT, we need some tolerance */
  if ((ensemble == ENS_NPT_ISO) || (ensemble == ENS_NPT_AXIAL)) {
    tol = SQR(1.0 + cell_size_tolerance);
  }
#endif

  /* compute scaling factors */
  cell_scale.x = sqrt( tol * cellsz / height.x );
  cell_scale.y = sqrt( tol * cellsz / height.y );
  cell_scale.z = sqrt( tol * cellsz / height.z );

  /* set up cell array dimensions */
  global_cell_dim.x = (int) ( 1.0 / cell_scale.x );
  global_cell_dim.y = (int) ( 1.0 / cell_scale.y );
  global_cell_dim.z = (int) ( 1.0 / cell_scale.z );

  if ((0 == myid ) && (0 == myrank))
  printf("Minimal cell size: \n\t ( %f %f %f ) \n\t ( %f %f %f ) \n\t ( %f %f %f )\n",
    box_x.x * cell_scale.x, box_x.y * cell_scale.x, box_x.z * cell_scale.x,
    box_y.x * cell_scale.y, box_y.y * cell_scale.y, box_y.z * cell_scale.y,
    box_z.x * cell_scale.z, box_z.y * cell_scale.z, box_z.z * cell_scale.z);

  /* global_cell_dim must be a multiple of cd */
#ifdef OMP
  cd.x = cpu_dim.x * 2;
  cd.y = cpu_dim.y * 2;
  cd.z = cpu_dim.z * 2;
#else
#ifdef TTM
  cd.x = cpu_dim.x * fd_ext.x;
  cd.y = cpu_dim.y * fd_ext.y;
  cd.z = cpu_dim.z * fd_ext.z;
#else
  cd.x = cpu_dim.x;
  cd.y = cpu_dim.y;
  cd.z = cpu_dim.z;
#endif /*TTM*/
#ifdef DEBUG
  if (cd.x % force_celldim_divisor.x != 0) cd.x *= force_celldim_divisor.x;
  if (cd.y % force_celldim_divisor.y != 0) cd.y *= force_celldim_divisor.y;
  if (cd.z % force_celldim_divisor.z != 0) cd.z *= force_celldim_divisor.z;
#endif /*DEBUG*/
#endif /*OMP*/
  if (0 != (global_cell_dim.x % cd.x))
     global_cell_dim.x = ((int)(global_cell_dim.x / cd.x)) * cd.x;
  if (0 != (global_cell_dim.y % cd.y))
     global_cell_dim.y = ((int)(global_cell_dim.y / cd.y)) * cd.y;
  if (0 != (global_cell_dim.z % cd.z))
     global_cell_dim.z = ((int)(global_cell_dim.z / cd.z)) * cd.z;

  /* Check if cell array is large enough */
  if ( 0 == myid ) {
    if (global_cell_dim.x < cd.x) {
      sprintf(msg,"global_cell_dim.x too small, need at least %d",cd.x);
      error(msg);
    }
    if (global_cell_dim.y < cd.y) {
      sprintf(msg,"global_cell_dim.y too small, need at least %d",cd.y);
      error(msg);
    }
    if (global_cell_dim.z < cd.z) {
      sprintf(msg,"global_cell_dim.z too small, need at least %d",cd.z);
      error(msg);
    }
  }

  /* if system grows, the next cell division should have more cells */
  next_cell_dim.x = global_cell_dim.x + cd.x;
  next_cell_dim.y = global_cell_dim.y + cd.y;
  next_cell_dim.z = global_cell_dim.z + cd.z;

  /* maximal and minimal heights before a new cell division is needed */
  min_height.x = cellsz * SQR(global_cell_dim.x);
  min_height.y = cellsz * SQR(global_cell_dim.y);
  min_height.z = cellsz * SQR(global_cell_dim.z);
  max_height.x = cellsz * SQR(  next_cell_dim.x) * tol;
  max_height.y = cellsz * SQR(  next_cell_dim.y) * tol;
  max_height.z = cellsz * SQR(  next_cell_dim.z) * tol;

  /* If an integer number of cells does not fit exactly into the box, the
     cells are enlarged accordingly */
  cell_scale.x = 1.0 / global_cell_dim.x;
  cell_scale.y = 1.0 / global_cell_dim.y;
  cell_scale.z = 1.0 / global_cell_dim.z;

  if ((0 == myid ) && (0 == myrank))
    printf("Actual cell size: \n\t ( %f %f %f ) \n\t ( %f %f %f ) \n\t ( %f %f %f )\n",
      box_x.x * cell_scale.x, box_x.y * cell_scale.x, box_x.z * cell_scale.x,
      box_y.x * cell_scale.y, box_y.y * cell_scale.y, box_y.z * cell_scale.y,
      box_z.x * cell_scale.z, box_z.y * cell_scale.z, box_z.z * cell_scale.z);

  if ((0 == myid ) && (0 == myrank))
    printf("Global cell array dimensions: %d %d %d\n",
      global_cell_dim.x,global_cell_dim.y,global_cell_dim.z);

#ifdef TTM
  /* set distances of FD lattice points in x,y,z-direction */
  /*      (fd_ext.?) * (MD cell size in ? direction) */
  fd_h.x = fd_ext.x * box_x.x * cell_scale.x;
  fd_h.y = fd_ext.y * box_y.y * cell_scale.y;
  fd_h.z = fd_ext.z * box_z.z * cell_scale.z;
#endif

  /* keep a copy of cell_dim, so that we can redistribute the atoms */
  cell_dim_old = cell_dim;

#ifdef BUFCELLS
  cell_dim.x = global_cell_dim.x / cpu_dim.x + 2;  
  cell_dim.y = global_cell_dim.y / cpu_dim.y + 2;
  cell_dim.z = global_cell_dim.z / cpu_dim.z + 2;

  cellmin.x = 1;   cellmax.x = cell_dim.x - 1;
  cellmin.y = 1;   cellmax.y = cell_dim.y - 1;
  cellmin.z = 1;   cellmax.z = cell_dim.z - 1;

  if ((0 == myid ) && (0 == myrank))
    printf("Local cell array dimensions (incl buffer): %d %d %d\n",
	   cell_dim.x,cell_dim.y,cell_dim.z);
#else
  cell_dim.x = global_cell_dim.x;
  cell_dim.y = global_cell_dim.y;
  cell_dim.z = global_cell_dim.z;

  cellmin.x = 0;   cellmax.x = cell_dim.x;
  cellmin.y = 0;   cellmax.y = cell_dim.y;
  cellmin.z = 0;   cellmax.z = cell_dim.z;

  printf("Local cell array dimensions: %d %d %d\n",
	 cell_dim.x,cell_dim.y,cell_dim.z);
#endif

  /* save old cell_array (if any), and allocate new one */
  cell_array_old = cell_array;
  cell_array = (minicell *) malloc(
               cell_dim.x * cell_dim.y * cell_dim.z * sizeof(minicell));

  if (0==myid)
    if (NULL == cell_array) error("Cannot allocate memory for cells");

  /* Initialize cells */
  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
      for (k=0; k < cell_dim.z; ++k) {

	p = PTR_3D_V(cell_array, i, j, k, cell_dim);
	p->n_max=0;
        p->n=0;
#ifdef BUFCELLS
        /* don't alloc data space for buffer cells */
        if ((0 != i) && (0 != j) && (0 != k) &&
            (i != cell_dim.x-1) &&
            (j != cell_dim.y-1) &&
            (k != cell_dim.z-1))
#endif

#ifdef debugLo
    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("print "" check init_cells start 111"" checking by Lo! \n");fflush(stdout);
#endif 

            ALLOC_MINICELL(p, initsz);

#ifdef debugLo
            printf("print "" check init_cells end "" checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("    ************************* \n");fflush(stdout);
#endif
  }

  /* on the first invocation we have to set up the MPI process topology */
#ifdef BUFCELLS
  if (cell_array_old == NULL) setup_mpi_topology();
#endif
  /* this is also the moment to inform about the number of threads */
#ifdef OMP
  if ((cell_array_old == NULL) && (myid == 0))
    printf("Computing with %d thread(s) per process.\n",omp_get_max_threads());
#endif

  /* redistribute atoms */
  if (cell_array_old != NULL) {
    for (j=0; j < cell_dim_old.x; j++)
      for (k=0; k < cell_dim_old.y; k++)
        for (l=0; l < cell_dim_old.z; l++) {
          p = PTR_3D_V(cell_array_old, j, k, l, cell_dim_old);

#ifdef BUFCELLS
          /* redistribute only contents of real cells */
          if ((0 != j) && (0 != k) && (0 != l) &&
              (j != cell_dim_old.x-1) &&
              (k != cell_dim_old.y-1) &&
              (l != cell_dim_old.z-1))
#endif
            for (i = p->n - 1; i >= 0; i--) {
              cellc = cell_coord( ORT(p,i,X), ORT(p,i,Y), ORT(p,i,Z) );
#ifdef BUFCELLS
              cellc = local_cell_coord( cellc );
              /* make sure atoms don't end up in buffer cells */
              if      (cellc.x <  cellmin.x) cellc.x = cellmin.x; 
              else if (cellc.x >= cellmax.x) cellc.x = cellmax.x-1;
              if      (cellc.y <  cellmin.y) cellc.y = cellmin.y; 
              else if (cellc.y >= cellmax.y) cellc.y = cellmax.y-1;
              if      (cellc.z <  cellmin.z) cellc.z = cellmin.z; 
              else if (cellc.z >= cellmax.z) cellc.z = cellmax.z-1;
#endif
              to = PTR_VV(cell_array,cellc,cell_dim);
              MOVE_ATOM( to, p, i );
            }

          ALLOC_MINICELL( p, 0 );  /* free old cell */
    }
    free(cell_array_old);

#ifdef debugLo

    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("passing "" go into make_cell_lists A ! "" checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);

#endif

    make_cell_lists();
    fix_cells();
#ifdef MPI
    setup_buffers();
#endif
  } else {

#ifdef debugLo
    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("passing "" go into make_cell_lists B ! "" checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("    ************************* \n");fflush(stdout);
#endif
    
    make_cell_lists();
  }
}


#ifndef NBLIST

/******************************************************************************
*
*  make_cell_lists creates a list of indices of all inner cells
*  (only if BUFCELLS), and a list of all pairs of interacting cells.
*  These lists make it easy to loop over these cells and pairs.
*
******************************************************************************/

void make_cell_lists(void)
{
  int i,j,k,l,m,n,r,s,t,nn,nnx,nny,nnz;
  ivektor ipbc, neigh;
  pair *P;

#ifdef debugLo
    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("passing "" make_cell_lists 1 start!, with npairs "" checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
#endif

#ifdef OMP
  nlists = 27;
#else
  nlists = 1;
#endif

  /* initialize pairs before creating the first pair lists */
  if (nallcells==0) for (i=0; i<nlists; ++i) pairs[i] = NULL;

  nallcells = cell_dim.x * cell_dim.y * cell_dim.z;
#ifdef BUFCELLS
  ncells = (cell_dim.x-2) * (cell_dim.y-2) * (cell_dim.z-2);
  /* make list of inner cell indices */
  cells  = (integer*) realloc( cells, ncells * sizeof(integer) );
  l = 0;
  for (i=cellmin.x; i<cellmax.x; ++i)
    for (j=cellmin.y; j<cellmax.y; ++j)
      for (k=cellmin.z; k<cellmax.z; ++k)
        cells[l++] = i * cell_dim.y * cell_dim.z + j * cell_dim.z + k;
#else
  ncells = cell_dim.x * cell_dim.y * cell_dim.z;
#endif

  /* Make lists with pairs of interacting cells, taking account of 
     the boundary conditions. We distribute pairs on several lists 
     such that among the pairs in any list there is no cell that 
     occurs twice. This allows to update forces independently
     for the pairs from the same list.
  */

#ifdef OMP
  nn = sizeof(pair) * (cell_dim.x * cell_dim.y * cell_dim.z);
  pairs[0] = (pair *) realloc( pairs[0], nn );  npairs[0] = 0;

  nn = sizeof(pair) * (cell_dim.x * cell_dim.y * ((cell_dim.z+1)/2));
  pairs[1] = (pair *) realloc( pairs[1], nn );  npairs[1] = 0;
  pairs[2] = (pair *) realloc( pairs[2], nn );  npairs[2] = 0;

  nn = sizeof(pair) * (cell_dim.x * ((cell_dim.y+1)/2) * cell_dim.z);
  for (i=3; i<9; i++) {
    pairs[i] = (pair *) realloc( pairs[i], nn );  npairs[i] = 0;
  }
  nn = sizeof(pair) * (((cell_dim.x+1)/2) * cell_dim.y * cell_dim.z);
  for (i=9; i<27; i++) {
    pairs[i] = (pair *) realloc( pairs[i], nn );  npairs[i] = 0;
  }

  if ((pairs[ 0]==NULL) || (pairs[ 1]==NULL) || (pairs[ 2]==NULL) ||
      (pairs[ 3]==NULL) || (pairs[ 4]==NULL) || (pairs[ 5]==NULL) ||
      (pairs[ 6]==NULL) || (pairs[ 7]==NULL) || (pairs[ 8]==NULL) ||
      (pairs[ 9]==NULL) || (pairs[10]==NULL) || (pairs[11]==NULL) ||
      (pairs[12]==NULL) || (pairs[13]==NULL) || (pairs[14]==NULL) ||
      (pairs[15]==NULL) || (pairs[16]==NULL) || (pairs[17]==NULL) ||
      (pairs[18]==NULL) || (pairs[19]==NULL) || (pairs[20]==NULL) ||
      (pairs[21]==NULL) || (pairs[22]==NULL) || (pairs[23]==NULL) ||
      (pairs[24]==NULL) || (pairs[25]==NULL) || (pairs[26]==NULL)) 
    error("cannot allocate pair lists");
#else
  nn = sizeof(pair) * cell_dim.x * cell_dim.y * cell_dim.z * 14;
  pairs[0] = (pair *) realloc( pairs[0], nn );  npairs[0] = 0;
  if (pairs[0]==NULL) error("cannot allocate pair list");
#endif

  /* for each cell */
  for (i=cellmin.x; i<cellmax.x; ++i)
    for (j=cellmin.y; j<cellmax.y; ++j)
      for (k=cellmin.z; k<cellmax.z; ++k) {

#ifdef OMP
        if (i % 2 == 0) nnx =  9; else nnx = 10;
        if (j % 2 == 0) nny =  3; else nny =  4;
        if (k % 2 == 0) nnz =  1; else nnz =  2;
#endif

	/* For half of the neighbours of this cell */
	for (l=0; l <= 1; ++l)
	  for (m=-l; m <= 1; ++m)
	    for (n=(l==0 ? -m  : -l ); n <= 1; ++n) { 

#ifdef OMP
              /* array where to put the pairs */
              if (l==0) {
                if (m==0) { 
                  if (n==0) nn = 0;
                  else    { nn = nnz; nnz += 2; }
		} else    { nn = nny; nny += 2; }
	      } else      { nn = nnx; nnx += 2; }
#else
              nn = 0;
#endif

#ifdef BUFCELLS
              r = i+l - 1 + my_coord.x * (cell_dim.x - 2);
              s = j+m - 1 + my_coord.y * (cell_dim.y - 2);
              t = k+n - 1 + my_coord.z * (cell_dim.z - 2);
#else
              r = i+l;
              s = j+m;
              t = k+n;
#endif

              /* Apply periodic boundaries */
              ipbc.x = 0;
              if (r<0) ipbc.x--; else if (r>global_cell_dim.x-1) ipbc.x++;

              ipbc.y = 0;
              if (s<0) ipbc.y--; else if (s>global_cell_dim.y-1) ipbc.y++;

              ipbc.z = 0;
              if (t<0) ipbc.z--; else if (t>global_cell_dim.z-1) ipbc.z++;

#ifdef BUFCELLS
              r = i+l;
              s = j+m;
              t = k+n;
#else
              if (r<0) r=cell_dim.x-1; 
              else if (r>cell_dim.x-1) r=0;

              if (s<0) s=cell_dim.y-1; 
              else if (s>cell_dim.y-1) s=0;

              if (t<0) t=cell_dim.z-1; 
              else if (t>cell_dim.z-1) t=0;
#endif

              if (((pbc_dirs.x==1) || (pbc_dirs.x==ipbc.x)) &&
                  ((pbc_dirs.y==1) || (pbc_dirs.y==ipbc.y)) &&
                  ((pbc_dirs.z==1) || (pbc_dirs.z==ipbc.z)))
              {
                /* add pair to list */
                P = pairs[nn] + npairs[nn];
                P->np = i*cell_dim.y*cell_dim.z + j*cell_dim.z + k;
                P->nq = r*cell_dim.y*cell_dim.z + s*cell_dim.z + t;
                P->ipbc[0] = ipbc.x;
                P->ipbc[1] = ipbc.y;
                P->ipbc[2] = ipbc.z;
                npairs[nn]++;
	      }
	    }
      }

#ifdef BUFCELLS

  /* If we don't use actio=reactio accross cpus, we have to do
     the force loop also on the other half of the neighbours for the 
     cells on the surface of the CPU */

  for (i=0; i<nlists; ++i) npairs2[i] = npairs[i];

  /* for each cell */
  for (i=cellmin.x; i<cellmax.x; ++i)
    for (j=cellmin.y; j<cellmax.y; ++j)
      for (k=cellmin.z; k<cellmax.z; ++k) {

#ifdef OMP
        if (i % 2 == 0) nnx =  9; else nnx = 10;
        if (j % 2 == 0) nny =  3; else nny =  4;
        if (k % 2 == 0) nnz =  1; else nnz =  2;
#endif

	/* for the other half of the neighbours of this cell */
	for (l=0; l <= 1; ++l)
	  for (m=-l; m <= 1; ++m)
	    for (n=(l==0 ? -m  : -l ); n <= 1; ++n) { 

              neigh.x = i-l;
              neigh.y = j-m;
              neigh.z = k-n;

#ifdef OMP
              /* array where to put the pairs */
              if (l==0) {
                if (m==0) { 
                  if (n==0) nn = 0;
                  else    { nn = nnz; nnz += 2; }
		} else    { nn = nny; nny += 2; }
	      } else      { nn = nnx; nnx += 2; }
#else
              nn = 0;
#endif

              /* if second cell is a buffer cell */
              if ((neigh.x == 0) || (neigh.x == cell_dim.x-1) || 
                  (neigh.y == 0) || (neigh.y == cell_dim.y-1) ||
                  (neigh.z == 0) || (neigh.z == cell_dim.z-1)) 
              {
                /* Apply periodic boundaries */
                ipbc.x = 0; r = neigh.x - 1 + my_coord.x * (cell_dim.x - 2);
                if (r<0) ipbc.x--; else if (r>global_cell_dim.x-1) ipbc.x++;
                r = neigh.x;

                ipbc.y = 0; s = neigh.y - 1 + my_coord.y * (cell_dim.y - 2);
                if (s<0) ipbc.y--; else if (s>global_cell_dim.y-1) ipbc.y++;
                s = neigh.y;

                ipbc.z = 0; t = neigh.z - 1 + my_coord.z * (cell_dim.z - 2);
                if (t<0) ipbc.z--; else if (t>global_cell_dim.z-1) ipbc.z++;
                t = neigh.z;

                if (((pbc_dirs.x==1) || (pbc_dirs.x==ipbc.x)) &&
                    ((pbc_dirs.y==1) || (pbc_dirs.y==ipbc.y)) &&
                    ((pbc_dirs.z==1) || (pbc_dirs.z==ipbc.z)))
                {
                  /* add pair to list */
                  P = pairs[nn] + npairs2[nn];
                  P->np = i*cell_dim.y*cell_dim.z + j*cell_dim.z + k;
                  P->nq = r*cell_dim.y*cell_dim.z + s*cell_dim.z + t;
                  P->ipbc[0] = ipbc.x;
                  P->ipbc[1] = ipbc.y;
                  P->ipbc[2] = ipbc.z;
                  npairs2[nn]++;
	        }
	      }
	    }
      }

#endif /* BUFCELLS */

#ifdef OMP
    check_pairs();
#endif


#ifdef debugLo
    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("passing "" make_cell_lists 1 end!, with npairs"" checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("    ************************* \n");fflush(stdout);
#endif

}

/******************************************************************************
*
*  check pair lists in OMP mode
*
******************************************************************************/

void check_pairs()
{
  int i, j, max, *lst;
  pair *p;

  max = cell_dim.x * cell_dim.y * cell_dim.z;
  lst = (int *) malloc(max*sizeof(int));

  for (i=0; i<nlists; i++) {
    for (j=0; j<max; j++) lst[j]=0;
    for (j=0; j<npairs[i]; j++) {
      p = pairs[i]+j;
      if (lst[p->np]>0)                     error("pair list corruption!"); 
      if (p->np>=max)                       error("pair overflow!");
      lst[p->np]=1;
      if ((lst[p->nq]>0) && (p->np!=p->nq)) error("pair list corruption!"); 
      if (p->nq>=max)                       error("pair overflow!");
      lst[p->nq]=1;
    }
  }
  free(lst);
}

#endif /* not NBLIST */

#ifdef NBLIST

/******************************************************************************
*
*  In the neighbor list version, make_cell_lists creates for each cell 
*  a list of neighbor cells
*
******************************************************************************/

void make_cell_lists(void)
{
  int i,j,k, l,m,n, r,s,t, nn, qq, flag, ncnbrs=0;
  cell_nbrs_t *CN;
  ivektor ipbc;

#ifdef debugLo

    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("passing "" make_cell_lists 2 start!, without npairs "" checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("    ************************* \n");fflush(stdout);
#endif

  nallcells = cell_dim.x * cell_dim.y * cell_dim.z;
  ncells = (cell_dim.x-2) * (cell_dim.y-2) * (cell_dim.z-2);
  /* make list of inner cell indices */
  cells  = (integer*) realloc( cells, ncells * sizeof(integer) );
  l = 0;
  for (i=cellmin.x; i<cellmax.x; ++i)
    for (j=cellmin.y; j<cellmax.y; ++j)
      for (k=cellmin.z; k<cellmax.z; ++k)
        cells[l++] = i * cell_dim.y * cell_dim.z + j * cell_dim.z + k;

  ncells2 = ncells;
  cnbrs = (cell_nbrs_t *) realloc( cnbrs, nallcells * sizeof(cell_nbrs_t) );
  if (cnbrs==NULL) error("cannot allocate cell neighbor list");
  CN = cnbrs;

  /* for each cell */
  for (i=cellmin.x; i<cellmax.x; ++i)
    for (j=cellmin.y; j<cellmax.y; ++j)
      for (k=cellmin.z; k<cellmax.z; ++k) {

        CN->np = i * cell_dim.y * cell_dim.z + j * cell_dim.z + k;
        CN->nq[0] = CN->np;
        nn = 1;

#ifdef AR
	/* For half of the neighbours of this cell */
	for (l=0; l <= 1; ++l)
	  for (m=-l; m <= 1; ++m)
	    for (n=(l==0 ? -m  : -l ); n <= 1; ++n)
#else
	/* For all neighbours of this cell */
	for (l=-1; l <= 1; ++l)
	  for (m=-1; m <= 1; ++m)
	    for (n=-1; n <= 1; ++n) 
#endif
	    {
              r = i+l;
              s = j+m;
              t = k+n;

              qq = r * cell_dim.y * cell_dim.z + s * cell_dim.z + t;
              if (qq==CN->np) continue;

              /* apply periodic boundaries */
              r += -1 + my_coord.x * (cell_dim.x - 2);
              s += -1 + my_coord.y * (cell_dim.y - 2);
              t += -1 + my_coord.z * (cell_dim.z - 2);

              ipbc.x = 0;
              if (r<0) ipbc.x--; else if (r>global_cell_dim.x-1) ipbc.x++;

              ipbc.y = 0;
              if (s<0) ipbc.y--; else if (s>global_cell_dim.y-1) ipbc.y++;

              ipbc.z = 0;
              if (t<0) ipbc.z--; else if (t>global_cell_dim.z-1) ipbc.z++;

              /* add cell to the list */
              if ( ((pbc_dirs.x==1) || (ipbc.x==0)) &&
                   ((pbc_dirs.y==1) || (ipbc.y==0)) &&
                   ((pbc_dirs.z==1) || (ipbc.z==0)) )
                CN->nq[nn] = qq;
              else
                CN->nq[nn] = -1;
              nn++;
	    }
        CN++;
      }

#if defined(COVALENT) || defined(NNBR_TABLE)

  /* for each cell */
  for (i=cellmin.x; i<cellmax.x; ++i)
    for (j=cellmin.y; j<cellmax.y; ++j)
      for (k=cellmin.z; k<cellmax.z; ++k) {

        CN->np = i * cell_dim.y * cell_dim.z + j * cell_dim.z + k;
        for (nn=0; nn<14; nn++) CN->nq[nn] = -1;
        flag = 0;
        nn   = 0;

	/* for the other half of the neighbours of this cell */
	for (l=0; l <= 1; ++l)
	  for (m=-l; m <= 1; ++m)
	    for (n=(l==0 ? -m  : -l ); n <= 1; ++n) {
              r = i-l;
              s = j-m;
              t = k-n;

              /* if second cell is a buffer cell */
              if ((r == 0) || (r == cell_dim.x-1) || 
                  (s == 0) || (s == cell_dim.y-1) ||
                  (t == 0) || (t == cell_dim.z-1)) {

                qq = r * cell_dim.y * cell_dim.z + s * cell_dim.z + t;
              
                /* apply periodic boundaries */
                r += -1 + my_coord.x * (cell_dim.x - 2);
                s += -1 + my_coord.y * (cell_dim.y - 2);
                t += -1 + my_coord.z * (cell_dim.z - 2);

                ipbc.x = 0;
                if (r<0) ipbc.x--; else if (r>global_cell_dim.x-1) ipbc.x++;

                ipbc.y = 0;
                if (s<0) ipbc.y--; else if (s>global_cell_dim.y-1) ipbc.y++;

                ipbc.z = 0;
                if (t<0) ipbc.z--; else if (t>global_cell_dim.z-1) ipbc.z++;

                /* add cell to the list */
                if ( ((pbc_dirs.x==1) || (ipbc.x==0)) &&
                     ((pbc_dirs.y==1) || (ipbc.y==0)) &&
                     ((pbc_dirs.z==1) || (ipbc.z==0)) ) {
                  CN->nq[nn] = qq;
                  flag=1;
		}
	      }
              nn++;
	    }

        /* increase pointer only if there were some neighbor BUFFER cells */
        if (flag) {
          CN++;
          ncells2++;
        }

      }

#endif /* COVALENT */

#ifdef debugLo
    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("passing "" make_cell_lists 2 end!, without npairs "" checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("    ************************* \n");fflush(stdout);
#endif

}

#endif

/******************************************************************************
*
*  cell_coord computes the (global) cell coordinates of a position
*
******************************************************************************/

ivektor cell_coord(real x, real y, real z)
{
  ivektor coord;

  /* Map positions to boxes */
  coord.x = (int)(global_cell_dim.x * (x*tbox_x.x + y*tbox_x.y + z*tbox_x.z));
  coord.y = (int)(global_cell_dim.y * (x*tbox_y.x + y*tbox_y.y + z*tbox_y.z));
  coord.z = (int)(global_cell_dim.z * (x*tbox_z.x + y*tbox_z.y + z*tbox_z.z));

  /* rounding errors may put atoms slightly outside the simulation cell */
  /* in the case of no pbc they may even be far outside */
  if      (coord.x >= global_cell_dim.x) coord.x = global_cell_dim.x - 1;
  else if (coord.x < 0)                  coord.x = 0;
  if      (coord.y >= global_cell_dim.y) coord.y = global_cell_dim.y - 1;
  else if (coord.y < 0)                  coord.y = 0;
  if      (coord.z >= global_cell_dim.z) coord.z = global_cell_dim.z - 1;
  else if (coord.z < 0)                  coord.z = 0;

  return coord;

}


/******************************************************************************
*
*  map vektor back into simulation box
*
******************************************************************************/

vektor back_into_box(vektor pos)
{
  real i;

  if (pbc_dirs.x==1) {
    i = FLOOR(SPROD(pos,tbox_x));
    pos.x  -= i *  box_x.x;
    pos.y  -= i *  box_x.y;
    pos.z  -= i *  box_x.z;
  }

  if (pbc_dirs.y==1) {
    i = FLOOR(SPROD(pos,tbox_y));
    pos.x  -= i *  box_y.x;
    pos.y  -= i *  box_y.y;
    pos.z  -= i *  box_y.z;
  }

  if (pbc_dirs.z==1) {
    i = FLOOR(SPROD(pos,tbox_z));
    pos.x  -= i *  box_z.x;
    pos.y  -= i *  box_z.y;
    pos.z  -= i *  box_z.z;
  }

  return pos;

}

