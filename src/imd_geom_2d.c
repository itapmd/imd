
/******************************************************************************
*
* imd_geom_2d.c -- domain decomposition routines, 2d version
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

/* compute transformation matrix */
void make_box( void )
{
  /* volume */
  volume = box_x.x * box_y.y - box_x.y * box_y.x;
  if ((0==myid) && (0==volume)) error("Box Edges are parallel.");

  /* compute tbox_j such that SPROD(box_i,tbox_j) == delta_ij */
  tbox_x.x =   box_y.y / volume;
  tbox_x.y = - box_y.x / volume;
  tbox_y.x = - box_x.y / volume;
  tbox_y.y =   box_x.x / volume;

  /* squares of the box heights perpendicular to the faces */
  height.x = 1.0 / SPROD(tbox_x,tbox_x);
  height.y = 1.0 / SPROD(tbox_y,tbox_y);
}

/* compute maximal dimension of cell array */
ivektor maximal_cell_dim( void )
{
  ivektor max_cell_dim;
  max_cell_dim.x = (int) (1.0 / sqrt( cellsz / height.x ) );
  max_cell_dim.y = (int) (1.0 / sqrt( cellsz / height.y ) );
  return max_cell_dim;
}


/******************************************************************************
*
*  Compute the size of the cells, and initialize the cell array.
*  There must be at least two cells in each direction.
*
******************************************************************************/

void init_cells( void )
{
  int i, j, k;
  real tmp; 
  vektor cell_scale;
  ivektor next_cell_dim, cell_dim_old;
  ivektor cellmin_old, cellmax_old, cellc;
  cell *p, *cell_array_old, *to; 

#ifdef NPT
  if (0 == myid) {
    if (ensemble == ENS_NPT_ISO) {
      printf("actual_shrink=%f limit_shrink=%f limit_growth=%f\n",
              actual_shrink.x, limit_shrink.x, limit_growth.x );
    }
    else if (ensemble == ENS_NPT_AXIAL) {
      printf("actual_shrink.x=%f limit_shrink.x=%f limit_growth.x=%f\n",
              actual_shrink.x,   limit_shrink.x,   limit_growth.x );
      printf("actual_shrink.y=%f limit_shrink.y=%f limit_growth.y=%f\n",
              actual_shrink.y,   limit_shrink.y,   limit_growth.y );
    }
  }
#endif

  /* compute scaling factors */
  make_box();
  cell_scale.x = sqrt( cellsz / height.x );
  cell_scale.y = sqrt( cellsz / height.y );

#ifdef NPT
  if ((ensemble == ENS_NPT_ISO) || (ensemble == ENS_NPT_AXIAL)) {
    /* the NEXT cell array for a GROWING system; 
       needed to determine when to recompute the cell division */
    next_cell_dim.x = (int) ( (1.0 + cell_size_tolerance) / cell_scale.x );
    next_cell_dim.y = (int) ( (1.0 + cell_size_tolerance) / cell_scale.y );
    cell_scale.x   /= (1.0 - cell_size_tolerance);
    cell_scale.y   /= (1.0 - cell_size_tolerance);
  }
#endif
  /* set up the CURRENT cell array dimensions */
  global_cell_dim.x = (int) ( 1.0 / cell_scale.x );
  global_cell_dim.y = (int) ( 1.0 / cell_scale.y );

  if (0 == myid )
  printf("Minimal cell size: \n\t ( %f %f ) \n\t ( %f %f ) \n",
	 box_x.x * cell_scale.x, box_x.y * cell_scale.x, 
	 box_y.x * cell_scale.y, box_y.y * cell_scale.y); 

#ifdef MPI
  /* cpu_dim must be a divisor of global_cell_dim */
  if (0 != (global_cell_dim.x % cpu_dim.x))
     global_cell_dim.x = ((int)(global_cell_dim.x/cpu_dim.x))*cpu_dim.x;
  if (0 != (global_cell_dim.y % cpu_dim.y))
     global_cell_dim.y = ((int)(global_cell_dim.y/cpu_dim.y))*cpu_dim.y;
#ifdef NPT
  /* cpu_dim must be a divisor of next_cell_dim */
  if ((ensemble == ENS_NPT_ISO) || (ensemble == ENS_NPT_AXIAL)) {
    if (0 != (next_cell_dim.x % cpu_dim.x))
       next_cell_dim.x = ((int)(next_cell_dim.x/cpu_dim.x))*cpu_dim.x;
    if (0 != (next_cell_dim.y % cpu_dim.y))
       next_cell_dim.y = ((int)(next_cell_dim.y/cpu_dim.y))*cpu_dim.y;
  }
#endif
#elif defined(OMP)
  /* global_cell_dim must be even */
  if (0 != (global_cell_dim.x % 2)) global_cell_dim.x -= 1;
  if (0 != (global_cell_dim.y % 2)) global_cell_dim.y -= 1;
#ifdef NPT
  /* next_cell_dim must be even */
  if ((ensemble == ENS_NPT_ISO) || (ensemble == ENS_NPT_AXIAL)) {
    if (0 != (next_cell_dim.x % 2)) next_cell_dim.x -= 1;
    if (0 != (next_cell_dim.y % 2)) next_cell_dim.y -= 1;
  }
#endif
#endif
  
  /* Check if cell array is large enough */
  if (0 == myid ) {
#ifdef MPI
    if (global_cell_dim.x < cpu_dim.x) error("global_cell_dim.x < cpu_dim.x");
    if (global_cell_dim.y < cpu_dim.y) error("global_cell_dim.y < cpu_dim.y");
#endif
    if (global_cell_dim.x < 2) error("global_cell_dim.x < 2");
    if (global_cell_dim.y < 2) error("global_cell_dim.y < 2");
  }

#ifdef NPT

  /* if system grows, the next cell division should have more cells */
  if ((ensemble == ENS_NPT_ISO) || (ensemble == ENS_NPT_AXIAL)) {
#ifdef MPI
    if (next_cell_dim.x == global_cell_dim.x) next_cell_dim.x += cpu_dim.x;
    if (next_cell_dim.y == global_cell_dim.y) next_cell_dim.y += cpu_dim.y;
#elif defined(OMP)
    if (next_cell_dim.x == global_cell_dim.x) next_cell_dim.x += 2;
    if (next_cell_dim.y == global_cell_dim.y) next_cell_dim.y += 2;
#else
    if (next_cell_dim.x == global_cell_dim.x) next_cell_dim.x += 1;
    if (next_cell_dim.y == global_cell_dim.y) next_cell_dim.y += 1;
#endif
  }

  if (ensemble == ENS_NPT_ISO) {
    /* factor by which a cell can grow before a change of
       the cell division becomes worthwhile */
    limit_growth.x = cell_scale.x * next_cell_dim.x;
    tmp            = cell_scale.y * next_cell_dim.y;
    /* getting more cells in at least one direction is enough */
    limit_growth.x = MIN( limit_growth.x, tmp );
    /* factor by which a cell can safely shrink */
    limit_shrink.x = cell_scale.x * global_cell_dim.x;
    tmp            = cell_scale.y * global_cell_dim.y;
    limit_shrink.x = MAX( limit_shrink.x, tmp );
    limit_shrink.x = limit_shrink.x * (1.0 - cell_size_tolerance);
  }

  if (ensemble == ENS_NPT_AXIAL) {
    /* factors by which a cell can grow before a change of
       the cell division becomes worthwhile */
    limit_growth.x = cell_scale.x * next_cell_dim.x;
    limit_growth.y = cell_scale.y * next_cell_dim.y;
    /* factor by which a cell can safely shrink */
    limit_shrink.x = cell_scale.x*global_cell_dim.x*(1.0-cell_size_tolerance);
    limit_shrink.y = cell_scale.y*global_cell_dim.y*(1.0-cell_size_tolerance);
  }

#endif /* NPT */

  /* If an integer number of cells does not fit exactly into the box, the
     cells are enlarged accordingly */
  cell_scale.x = 1.0 / global_cell_dim.x;
  cell_scale.y = 1.0 / global_cell_dim.y;

  if (0 == myid) {
    printf("Actual cell size: \n\t ( %f %f ) \n\t ( %f %f ) \n",
           box_x.x * cell_scale.x, box_x.y * cell_scale.x,
           box_y.x * cell_scale.y, box_y.y * cell_scale.y);
    printf("Global cell array dimensions: %d %d\n",
           global_cell_dim.x,global_cell_dim.y);
  }

  /* keep a copy of cell_dim & Co., so that we can redistribute the atoms */
  cell_dim_old = cell_dim;
  cellmin_old  = cellmin;
  cellmax_old  = cellmax;

#ifdef MPI
  cell_dim.x = global_cell_dim.x / cpu_dim.x + 2;  
  cell_dim.y = global_cell_dim.y / cpu_dim.y + 2;
  
  cellmin.x = 1;   cellmax.x = cell_dim.x - 1;
  cellmin.y = 1;   cellmax.y = cell_dim.y - 1;

  if (0==myid) 
    printf("Local cell array dimensions (incl buffer): %d %d\n",
	   cell_dim.x,cell_dim.y);
#else
  cell_dim.x = global_cell_dim.x;
  cell_dim.y = global_cell_dim.y;

  cellmin.x = 0;   cellmax.x = cell_dim.x;
  cellmin.y = 0;   cellmax.y = cell_dim.y;

  printf("Local cell array dimensions: %d %d\n", cell_dim.x, cell_dim.y);
#endif

  /* save old cell_array (if any), and allocate new one */
  cell_array_old = cell_array;
  cell_array = (cell *) malloc( cell_dim.x * cell_dim.y * sizeof(cell) );

  if (0 == myid)
  if (NULL == cell_array) error("Cannot allocate memory for cells");

  /* initialize cells */
  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j) {
      p = PTR_2D_V(cell_array, i, j, cell_dim);
      p->n_max=0;
      p->n=0;
#ifdef MPI
      /* don't alloc data space for buffer cells */
      if ((0 != i) && (0 != j) && (i != cell_dim.x-1) && (j != cell_dim.y-1))
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
      for (k=cellmin_old.y; k < cellmax_old.y; k++) {
	p = PTR_2D_V(cell_array_old, j, k, cell_dim_old);
        for (i = p->n - 1; i >= 0; i--) {
#ifdef MPI
          cellc = local_cell_coord(p->ort X(i),p->ort Y(i));
          /* strangly, some atoms get into buffer cells; 
             we push them back into the real cells, 
             so that we do not lose them  */
          if (cellc.x == 0) cellc.x++;
          if (cellc.y == 0) cellc.y++;
          if (cellc.x == cellmax.x) cellc.x--;
          if (cellc.y == cellmax.y) cellc.y--;
#else
          cellc = cell_coord(p->ort X(i),p->ort Y(i));
#endif
          to = PTR_VV(cell_array,cellc,cell_dim);
          move_atom( to, p, i );
        }
        alloc_cell( p, 0 );  /* free old cell */
    }
    free(cell_array_old);
  }

#ifdef NPT

  if ((ensemble == ENS_NPT_ISO) || (ensemble == ENS_NPT_AXIAL)) {
    revise_cell_division = 0;
    actual_shrink.x      = 1.0;
    actual_shrink.y      = 1.0;
  }

  if (0==myid) {
    if (ensemble == ENS_NPT_ISO) {
      printf("actual_shrink=%f limit_shrink=%f limit_growth=%f\n",
              actual_shrink.x, limit_shrink.x, limit_growth.x );
    } 
    else if (ensemble == ENS_NPT_AXIAL) { 
      printf("actual_shrink.x=%f limit_shrink.x=%f limit_growth.x=%f\n",
              actual_shrink.x,   limit_shrink.x,   limit_growth.x );
      printf("actual_shrink.y=%f limit_shrink.y=%f limit_growth.y=%f\n",
              actual_shrink.y,   limit_shrink.y,   limit_growth.y );
    } 
  }

#endif /* NPT */

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
  int i,j,k,l,m,r,s,nn,nnx,nny;
  ivektor ipbc, neigh;
  pair *P;

#ifdef OMP
  nlists = 9;
#else
  nlists = 1;
#endif

  /* initialize pairs before creating the first pair lists */
  if (nallcells==0) for (i=0; i<nlists; ++i) pairs[i] = NULL;

  nallcells = cell_dim.x * cell_dim.y;
#ifdef MPI
  ncells = (cell_dim.x-2) * (cell_dim.y-2);
  /* make list of inner cell indices */
  cells  = (integer*) realloc( cells, ncells * sizeof(integer) );
  k = 0;
  for (i=cellmin.x; i<cellmax.x; ++i)
    for (j=cellmin.y; j<cellmax.y; ++j)
      cells[k++] = i * cell_dim.y + j;
#else
  ncells = cell_dim.x * cell_dim.y;
#endif

  /* Make lists with pairs of interacting cells, taking account of 
     the boundary conditions. We distribute pairs on several lists 
     such that among the pairs in any list there is no cell that 
     occurs twice. This allows to update forces independently
     for the pairs from the same list.
  */

#ifdef OMP
  nn = sizeof(pair) * (cell_dim.x * cell_dim.y);
  pairs[0] = (pair *) realloc( pairs[0], nn );  npairs[0] = 0;
  nn = sizeof(pair) * (cell_dim.x * ((cell_dim.y+1)/2));
  pairs[1] = (pair *) realloc( pairs[1], nn );  npairs[1] = 0;
  pairs[2] = (pair *) realloc( pairs[2], nn );  npairs[2] = 0;
  nn = sizeof(pair) * (((cell_dim.x+1)/2) * cell_dim.y);
  pairs[3] = (pair *) realloc( pairs[3], nn );  npairs[3] = 0;
  pairs[4] = (pair *) realloc( pairs[4], nn );  npairs[4] = 0;
  pairs[5] = (pair *) realloc( pairs[5], nn );  npairs[5] = 0;
  pairs[6] = (pair *) realloc( pairs[6], nn );  npairs[6] = 0;
  pairs[7] = (pair *) realloc( pairs[7], nn );  npairs[7] = 0;
  pairs[8] = (pair *) realloc( pairs[8], nn );  npairs[8] = 0;
  if ((pairs[0]==NULL) || (pairs[1]==NULL) || (pairs[2]==NULL) ||
      (pairs[3]==NULL) || (pairs[4]==NULL) || (pairs[5]==NULL) ||
      (pairs[6]==NULL) || (pairs[7]==NULL) || (pairs[8]==NULL)) 
    error("cannot allocate pair lists");
#else
  nn = sizeof(pair) * cell_dim.x * cell_dim.y * 5;
  pairs[0] = (pair *) realloc( pairs[0], nn );  npairs[0] = 0;
  if (pairs[0]==NULL) error("cannot allocate pair list");
#endif

  /* for each cell */
  for (i=cellmin.x; i<cellmax.x; ++i)
    for (j=cellmin.y; j<cellmax.y; ++j) {

#ifdef OMP
      if (i % 2 == 0) nnx = 3; else nnx = 4;
      if (j % 2 == 0) nny = 1; else nny = 2;
#endif

      /* for half of the neighbours of this cell */
      for (l=0; l<=1; ++l)
        for (m=-l; m<=1; ++m) {

#ifdef OMP
          /* array where to put the pairs */
          if (l==0) {
            if (m==0) nn = 0; 
            else      nn = nny;
          } else {    nn = nnx; nnx +=2; }
#else
          nn = 0;
#endif

#ifdef MPI
          r = i+l - 1 + my_coord.x * (cell_dim.x - 2);
          s = j+m - 1 + my_coord.y * (cell_dim.y - 2);
#else
          r = i+l;
          s = j+m;
#endif

          /* Apply periodic boundaries */
          ipbc.x = 0;
          if (r<0) ipbc.x--; else if (r>global_cell_dim.x-1) ipbc.x++;

          ipbc.y = 0;
          if (s<0) ipbc.y--; else if (s>global_cell_dim.y-1) ipbc.y++;

#ifdef MPI
          r = i+l;
          s = j+m;
#else
          if (r<0) r=cell_dim.x-1; 
          else if (r>cell_dim.x-1) r=0;

          if (s<0) s=cell_dim.y-1; 
          else if (s>cell_dim.y-1) s=0;
#endif

          if (((pbc_dirs.x==1) || (pbc_dirs.x==ipbc.x)) &&
              ((pbc_dirs.y==1) || (pbc_dirs.y==ipbc.y)))
          {
            /* add pair to list */
            P = pairs[nn] + npairs[nn];
            P->np = i * cell_dim.y + j;
            P->nq = r * cell_dim.y + s;
            P->ipbc[0] = ipbc.x;
            P->ipbc[1] = ipbc.y;
            npairs[nn]++;
          }
	}
    }

#ifdef MPI

  /* If we do not use actio=reactio accross cpus, we have to do
     the force loop also on the other half of the neighbours for the 
     cells on the surface of the CPU */

  for (i=0; i<nlists; ++i) npairs2[i] = npairs[i];

  /* for each cell */
  for (i=cellmin.x; i<cellmax.x; ++i)
    for (j=cellmin.y; j<cellmax.y; ++j) {

#ifdef OMP
     if (i % 2 == 0) nnx = 3; else nnx = 4;
     if (j % 2 == 0) nny = 1; else nny = 2;
#endif

     /* for the other half of the neighbours of this cell */
      for (l=0; l<=1; ++l)
        for (m=-l; m<=1; ++m) {

          neigh.x = i-l;
          neigh.y = j-m;

#ifdef OMP
          /* array where to put the pairs */
          if (l==0) {
            if (m==0) nn = 0; 
            else      nn = nny;
          } else {    nn = nnx; nnx +=2; }
#else
          nn = 0;
#endif

          /* if second cell is a buffer cell */
          if ((neigh.x == 0) || (neigh.x == cell_dim.x-1) || 
              (neigh.y == 0) || (neigh.y == cell_dim.y-1)) 
          {

            /* Apply periodic boundaries */
            ipbc.x = 0; r = neigh.x - 1 + my_coord.x * (cell_dim.x - 2);
            if (r<0) ipbc.x--; else if (r>global_cell_dim.x-1) ipbc.x++;
            r = neigh.x;

            ipbc.y = 0; s = neigh.y - 1 + my_coord.y * (cell_dim.y - 2);
            if (s<0) ipbc.y--; else if (s>global_cell_dim.y-1) ipbc.y++;
            s = neigh.y;

            if (((pbc_dirs.x==1) || (pbc_dirs.x==ipbc.x)) &&
                ((pbc_dirs.y==1) || (pbc_dirs.y==ipbc.y)))
            {
              /* add pair to list */
              P = pairs[nn] + npairs2[nn];
              P->np = i * cell_dim.y + j;
              P->nq = r * cell_dim.y + s;
              P->ipbc[0] = ipbc.x;
              P->ipbc[1] = ipbc.y;
              npairs2[nn]++;
	    }
	  }
        }
    }

#endif /* MPI */

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

  max = global_cell_dim.x * global_cell_dim.y;
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

ivektor2d cell_coord(real x,real y)
{
  ivektor2d coord;

  /* Map positions to boxes */
  coord.x = (int) ( global_cell_dim.x * (x * tbox_x.x + y * tbox_x.y) );
  coord.y = (int) ( global_cell_dim.y * (x * tbox_y.x + y * tbox_y.y) );

  /* rounding errors may put atoms slightly outside the simulation cell */
  /* in the case of no pbc they may even be far outside */
  if      (coord.x >= global_cell_dim.x) coord.x = global_cell_dim.x - 1;
  else if (coord.x < 0)                  coord.x = 0;
  if      (coord.y >= global_cell_dim.y) coord.y = global_cell_dim.y - 1;
  else if (coord.y < 0)                  coord.y = 0;

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
  }

  if (pbc_dirs.y==1) {
    i = FLOOR(SPROD(pos,tbox_y));
    pos.x -= i *  box_y.x;
    pos.y -= i *  box_y.y;
  }

  return pos;

}


