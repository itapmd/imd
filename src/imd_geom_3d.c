
/******************************************************************************
*
* imd_geom_3d.c -- domain decomposition routines, 3d version
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

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

/* compute transformation matrix */
void make_box( void )
{
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

  /* squares of the box heights perpendicular to the faces */
  height.x = 1.0 / SPROD(tbox_x,tbox_x);
  height.y = 1.0 / SPROD(tbox_y,tbox_y);
  height.z = 1.0 / SPROD(tbox_z,tbox_z);
}


/******************************************************************************
*
*  Compute the size of the cells, and initialize the cell array.
*  There must be at least two cells in each direction.
*
******************************************************************************/

void init_cells( void )
{
  int i, j, k, l;
  real tmp;
  vektor cell_scale;
  ivektor next_cell_dim, cell_dim_old;
  ivektor cellmin_old, cellmax_old, cellc;
  cell *p, *cell_array_old, *to;

  /* compute scaling factors */
  make_box();
  cell_scale.x = sqrt( cellsz / height.x );
  cell_scale.y = sqrt( cellsz / height.y );
  cell_scale.z = sqrt( cellsz / height.z );

#ifdef NPT
  /* if NPT, we need some tolerance */
  if ((ensemble == ENS_NPT_ISO) || (ensemble == ENS_NPT_AXIAL)) {
    cell_scale.x *= (1.0 + cell_size_tolerance);
    cell_scale.y *= (1.0 + cell_size_tolerance);
    cell_scale.z *= (1.0 + cell_size_tolerance);
  }
#endif

  /* set up cell array dimensions */
  global_cell_dim.x = (int) ( 1.0 / cell_scale.x );
  global_cell_dim.y = (int) ( 1.0 / cell_scale.y );
  global_cell_dim.z = (int) ( 1.0 / cell_scale.z );

  if (0 == myid )
  printf("Minimal cell size: \n\t ( %f %f %f ) \n\t ( %f %f %f ) \n\t ( %f %f %f )\n",
    box_x.x * cell_scale.x, box_x.y * cell_scale.x, box_x.z * cell_scale.x,
    box_y.x * cell_scale.y, box_y.y * cell_scale.y, box_y.z * cell_scale.y,
    box_z.x * cell_scale.z, box_z.y * cell_scale.z, box_z.z * cell_scale.z);

#ifdef MPI
  /* cpu_dim must be a divisor of global_cell_dim */
  if (0 != (global_cell_dim.x % cpu_dim.x))
     global_cell_dim.x = ((int)(global_cell_dim.x/cpu_dim.x))*cpu_dim.x;
  if (0 != (global_cell_dim.y % cpu_dim.y))
     global_cell_dim.y = ((int)(global_cell_dim.y/cpu_dim.y))*cpu_dim.y;
  if (0 != (global_cell_dim.z % cpu_dim.z))
     global_cell_dim.z = ((int)(global_cell_dim.z/cpu_dim.z))*cpu_dim.z;
#elif defined(MPI) && defined(OMP) 
  /* cpu_dim must be a divisor of 2 * global_cell_dim */
  if (0 != (global_cell_dim.x % (2*cpu_dim.x)))
     global_cell_dim.x = ((int)(global_cell_dim.x/(2*cpu_dim.x)))*2*cpu_dim.x;
  if (0 != (global_cell_dim.y % cpu_dim.y))
     global_cell_dim.y = ((int)(global_cell_dim.y/(2*cpu_dim.y)))*2*cpu_dim.y;
  if (0 != (global_cell_dim.z % cpu_dim.z))
     global_cell_dim.z = ((int)(global_cell_dim.z/(2*cpu_dim.z)))*2*cpu_dim.z;
#elif defined(OMP)
  /* global_cell_dim must be even */
  if (0 != (global_cell_dim.x % 2)) global_cell_dim.x -= 1;
  if (0 != (global_cell_dim.y % 2)) global_cell_dim.y -= 1;
  if (0 != (global_cell_dim.z % 2)) global_cell_dim.z -= 1;
#endif

  /* Check if cell array is large enough */
  if ( 0 == myid ) {
#ifdef MPI
    if (global_cell_dim.x < cpu_dim.x) error("global_cell_dim.x < cpu_dim.x");
    if (global_cell_dim.y < cpu_dim.y) error("global_cell_dim.y < cpu_dim.y");
    if (global_cell_dim.z < cpu_dim.z) error("global_cell_dim.z < cpu_dim.z");
#endif
    if (global_cell_dim.x < 2) error("global_cell_dim.x < 2");
    if (global_cell_dim.y < 2) error("global_cell_dim.y < 2");
    if (global_cell_dim.z < 2) error("global_cell_dim.z < 2");
  }

#ifdef MPI
  next_cell_dim.x = global_cell_dim.x + cpu_dim.x;
  next_cell_dim.y = global_cell_dim.y + cpu_dim.y;
  next_cell_dim.z = global_cell_dim.z + cpu_dim.z;
#elif defined(MPI) && defined(OMP)
  next_cell_dim.x = global_cell_dim.x + 2 * cpu_dim.x;
  next_cell_dim.y = global_cell_dim.y + 2 * cpu_dim.y;
  next_cell_dim.z = global_cell_dim.z + 2 * cpu_dim.z;
#elif defined(OMP)
  next_cell_dim.x = global_cell_dim.x + 2;
  next_cell_dim.y = global_cell_dim.y + 2;
  next_cell_dim.z = global_cell_dim.z + 2;
#else
  next_cell_dim.x = global_cell_dim.x + 1;
  next_cell_dim.y = global_cell_dim.y + 1;
  next_cell_dim.z = global_cell_dim.z + 1;
#endif

  /* maximal and minimal heights before a new cell division is needed */
  min_height.x = cellsz * SQR(global_cell_dim.x);
  min_height.y = cellsz * SQR(global_cell_dim.y);
  min_height.z = cellsz * SQR(global_cell_dim.z);
  max_height.x = cellsz * SQR(  next_cell_dim.x);
  max_height.y = cellsz * SQR(  next_cell_dim.y);
  max_height.z = cellsz * SQR(  next_cell_dim.z);

#ifdef NPT
  /* if NPT, we should let it grow a a little more */
  if ((ensemble == ENS_NPT_ISO) || (ensemble == ENS_NPT_AXIAL)) {
    max_height.x *= SQR(1 + cell_size_tolerance);
    max_height.y *= SQR(1 + cell_size_tolerance);
    max_height.z *= SQR(1 + cell_size_tolerance);
  }
#endif

  /* If an integer number of cells does not fit exactly into the box, the
     cells are enlarged accordingly */
  cell_scale.x = 1.0 / global_cell_dim.x;
  cell_scale.y = 1.0 / global_cell_dim.y;
  cell_scale.z = 1.0 / global_cell_dim.z;

  if (0 == myid )
  printf("Actual cell size: \n\t ( %f %f %f ) \n\t ( %f %f %f ) \n\t ( %f %f %f )\n",
    box_x.x * cell_scale.x, box_x.y * cell_scale.x, box_x.z * cell_scale.x,
    box_y.x * cell_scale.y, box_y.y * cell_scale.y, box_y.z * cell_scale.y,
    box_z.x * cell_scale.z, box_z.y * cell_scale.z, box_z.z * cell_scale.z);

  if (0==myid)
  printf("Global cell array dimensions: %d %d %d\n",
          global_cell_dim.x,global_cell_dim.y,global_cell_dim.z);

  /* keep a copy of cell_dim & Co., so that we can redistribute the atoms */
  cell_dim_old = cell_dim;
  cellmin_old  = cellmin;
  cellmax_old  = cellmax;

#ifdef MPI
  cell_dim.x = global_cell_dim.x / cpu_dim.x + 2;  
  cell_dim.y = global_cell_dim.y / cpu_dim.y + 2;
  cell_dim.z = global_cell_dim.z / cpu_dim.z + 2;

  cellmin.x = 1;   cellmax.x = cell_dim.x - 1;
  cellmin.y = 1;   cellmax.y = cell_dim.y - 1;
  cellmin.z = 1;   cellmax.z = cell_dim.z - 1;

  if (0==myid) 
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
  cell_array = (cell *) malloc(
		     cell_dim.x * cell_dim.y * cell_dim.z * sizeof(cell));
  if ( 0 == myid )
    if (NULL == cell_array) error("Cannot allocate memory for cells");

  /* Initialize cells */
  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
      for (k=0; k < cell_dim.z; ++k) {

	p = PTR_3D_V(cell_array, i, j, k, cell_dim);
	p->n_max=0;
        p->n=0;
#ifdef MPI
        /* don't alloc data space for buffer cells */
        if ((0 != i) && (0 != j) && (0 != k) &&
            (i != cell_dim.x-1) &&
            (j != cell_dim.y-1) &&
            (k != cell_dim.z-1))
#endif
            alloc_cell(p, initsz);
  }

  /* on the first invocation we have to set up the MPI process topology */
#ifdef MPI
  if (cell_array_old == NULL) setup_mpi_topology();
#endif

  /* redistribute atoms */
  if (cell_array_old != NULL) {
    for (j=cellmin_old.x; j < cellmax_old.x; j++)
      for (k=cellmin_old.y; k < cellmax_old.y; k++)
        for (l=cellmin_old.z; l < cellmax_old.z; l++) {
          p = PTR_3D_V(cell_array_old, j, k, l, cell_dim_old);
          for (i = p->n - 1; i >= 0; i--) {
#ifdef MPI
            cellc = local_cell_coord(p->ort X(i),p->ort Y(i),p->ort Z(i));
            /* strangly, some atoms get into buffer cells; 
               we push them back into the real cells, 
               so that we don't lose them  */
            if (cellc.x == 0) cellc.x++;
            if (cellc.y == 0) cellc.y++;
            if (cellc.z == 0) cellc.z++;
            if (cellc.x == cellmax.x) cellc.x--;
            if (cellc.y == cellmax.y) cellc.y--;
            if (cellc.z == cellmax.z) cellc.z--;
#else
            cellc = cell_coord(p->ort X(i),p->ort Y(i),p->ort Z(i));
#endif
            to = PTR_VV(cell_array,cellc,cell_dim);
            move_atom( to, p, i );
          }
          alloc_cell( p, 0 );  /* free old cell */
    }
    free(cell_array_old);
  }

  make_cell_lists();

}


/******************************************************************************
*
*  make_cell_lists creates a list of indices of all inner cells
*  (only if MPI), and a list of all pairs of interacting cells.
*  These lists make it easy to loop over these cells and pairs.
*
******************************************************************************/

void make_cell_lists(void)
{
  int i,j,k,l,m,n,r,s,t,nn,nnx,nny,nnz;
  ivektor ipbc, neigh;
  pair *P;

#ifdef OMP
  nlists = 27;
#else
  nlists = 1;
#endif

  /* initialize pairs before creating the first pair lists */
  if (nallcells==0) for (i=0; i<nlists; ++i) pairs[i] = NULL;

  nallcells = cell_dim.x * cell_dim.y * cell_dim.z;
#ifdef MPI
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

#if !(defined(MPI) && defined(MONOLJ)) /* i.e., not world record */ 

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

#ifdef MPI
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

#ifdef MPI
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

#ifdef MPI

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

#endif /* MPI */

#endif /* !(defined(MPI) && definded(MONOLJ)), i.e., not world record */

#ifdef OMP
    check_pairs();
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

  max = global_cell_dim.x * global_cell_dim.y * global_cell_dim.y;
  lst = (int *) malloc(max*sizeof(int));

  for (i=0; i<nlists; i++) {
    for (j=0; j<max; j++) lst[j]=0;
    for (j=0; j<npairs[i]; j++) {
      p = pairs[i]+j;
      if (lst[p->np]>0)                     error("pair list corruption!"); 
      lst[p->np]=1;
      if ((lst[p->nq]>0) && (p->np!=p->nq)) error("pair list corruption!"); 
      lst[p->nq]=1;
    }
  }
  free(lst);
}


/******************************************************************************
*
*  cell_coord computes the (global) cell coorinates of a position
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

  return(coord);

}


/******************************************************************************
*
*  map vektor back into simulation box
*
******************************************************************************/

vektor back_into_box(vektor pos)
{
  int i;

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

