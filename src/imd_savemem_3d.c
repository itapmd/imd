
/******************************************************************************
*
* imd_savemem_3d.c -- routines for SAVEMEM option, three dimensions
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/******************************************************************************
*
*  global_cell_coord computes the global coordinates of a cell from the
*  local ones. PBC are always applied. Is used only for SAVEMEM.
*
******************************************************************************/

ivektor global_cell_coord(ivektor local_coord)
{
  ivektor global_coord;

  global_coord.x = local_coord.x - 1 + my_coord.x * (cell_dim.x - 2);
  global_coord.y = local_coord.y - 1 + my_coord.y * (cell_dim.y - 2);
  global_coord.z = local_coord.z - 1 + my_coord.z * (cell_dim.z - 2);

  if (global_coord.x < 0)                  global_coord.x += global_cell_dim.x;
  if (global_coord.x >= global_cell_dim.x) global_coord.x -= global_cell_dim.x;
  if (global_coord.y < 0)                  global_coord.y += global_cell_dim.y;
  if (global_coord.y >= global_cell_dim.y) global_coord.y -= global_cell_dim.y;
  if (global_coord.z < 0)                  global_coord.z += global_cell_dim.z;
  if (global_coord.z >= global_cell_dim.z) global_coord.z -= global_cell_dim.z;

  return(global_coord);
}


/******************************************************************************
*
* fix_cells_by_cell
*
* check if each atom is in the correct cell and on the correct CPU 
* move atoms the have left their cells in the last timestep
*
* This does not use Plimptons Comm Scheme
* Takes more time but saves on memory
*
******************************************************************************/

void fix_cells_by_cell(void)

{
  int i,j,k,l;
  int r,s,t;
  int odd;
  int tag;
  cell *p, b;
  ivektor coord,opposite,buf_coord;
  int to_cpu,from_cpu;
  int sum;
  int count;

  /* empty border cells and deallocate */
  empty_buffer_cells();
  
      /* move atoms, use border cells when necessary */

  b.n_max = 0;
  b.n     = 0;

  /* for each cell in bulk */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.y-1; ++j)
      for (k=1; k < cell_dim.z-1; ++k) {

	p = PTR_3D_V(cell_array, i, j, k, cell_dim);
	/* loop over atoms in cell */

	l=0;
	while( l<p->n ) {

	  coord = local_cell_coord(p->ort X(l),p->ort Y(l), p->ort Z(l));

	  /* see if atom is in wrong cell */
	  if ((coord.x != i) || (coord.y != j) || (coord.z != k)) {

	    /* check if atom is on correct CPU */
	    to_cpu = cpu_coord( global_cell_coord( coord ));
	    
	    if (myid!=to_cpu) {
	      /* atom is moved to different cpu */
	      tag = 0;
	      /* Search through buffer cells for correct cpu, 
                 set flag if found */
	      for (r = 0; r <= 2; ++r)
		for (s = 0; s <= 2; ++s)
		  for (t = 0; t <= 2; ++t) {
		    buf_coord.x = r;
		    buf_coord.y = s;
		    buf_coord.z = t;
		    if (2 == buf_coord.x) buf_coord.x = cell_dim.x - 1;
		    if (2 == buf_coord.y) buf_coord.y = cell_dim.y - 1;
		    if (2 == buf_coord.z) buf_coord.z = cell_dim.z - 1;
		    if (to_cpu == cpu_coord( global_cell_coord( buf_coord ))) {
		      tag = 1;
		      coord.x = buf_coord.x;
		      coord.y = buf_coord.y;
		      coord.z = buf_coord.z;
		    };
		  };
	      if (0==tag)  {
		error("No matching buffer cell.");
	      };
	    };

	    move_atom(coord, p, l);

	  } else {
	    /* move_atom shuffles the atom in cell p around, 
	       so we have to recheck atom number l after a move. */
	    ++l;
	  };

	};
      };

  /* send border cells to neighbbours */
  for (i = 0; i <= 2; ++i)
    for (j = 0; j <= 2; ++j)
      for (k = 0; k <= 2; ++k) {

	if      (1!=i) odd = my_coord.x % 2;
	else if (1!=j) odd = my_coord.y % 2;
	else           odd = my_coord.z % 2;

	coord.x = i;
	coord.y = j;
	coord.z = k;

	if (2 == coord.x) coord.x = cell_dim.x - 1;
	if (2 == coord.y) coord.y = cell_dim.y - 1;
	if (2 == coord.z) coord.z = cell_dim.z - 1;

	opposite.x = i;
	opposite.y = j;
	opposite.z = k;

	if (0==i)
	  opposite.x = cell_dim.x-1;
	if (2==i)
	  opposite.x = 0;
	if (0==j)
	  opposite.y = cell_dim.y-1;
	if (2==j)
	  opposite.y = 0;
	if (0==k)
	  opposite.z = cell_dim.z-1;
	if (2==k)
	  opposite.z = 0;

	to_cpu   = cpu_coord( global_cell_coord( coord ));
	from_cpu = cpu_coord( global_cell_coord( opposite ));

#ifndef MONOLJ
	tag = PTR_3D_VV(CELL_TAG, coord, cell_dim);
#else
	tag = ORT_TAG;
#endif
	p   = PTR_3D_VV(cell_array, coord, cell_dim);

	if (to_cpu != myid) {
	  if (0 == odd) {
	    send_cell( p, to_cpu, tag );
	    recv_cell( &b, from_cpu, tag );
	  } else {
	    recv_cell( &b, from_cpu, tag );
	    send_cell( p, to_cpu, tag );
	  };

	  /* move_atom decrements b.n internally */
	  while(0 < b.n) {
	    coord = local_cell_coord(b.ort X(0), b.ort Y(0), b.ort Z(0));
	    if ( myid !=  cpu_coord( global_cell_coord( coord ))) {
	      error("Atom moved to wrong CPU.");
	    }
	    move_atom(coord, &b, 0);
	  };
	  if (0<b.n) error("Atom left in buffer.");
	};
      };

  /* Deallocate temorary workspace */
  alloc_cell(&b,0);
    
}


/******************************************************************************
*
* send_atoms_by_cell
*
* This sends atoms to the neighbours on a cell by cell basis
*
* It generates more MPI Calls that send_atoms, but saves on memory
*
******************************************************************************/

void send_atoms_by_cell()
{
  int i,j,k,l,m,n,u,v,w;
  ivektor neighbour;
  ivektor nbcoord;
  int nbrank_cell;
  int nbrank_dir;
  
  /* send postitions to neighbours */

  /* for each cell with local data, we generate its neigbour cells.
  *  for each neighbour, we check if its located on our processor if not,
  *  the neighbour is in a buffer cell and we calculate the direction to
  *  the neighbours cpu we receive data from the neighbours cpu into its
  *  buffer cell we send the local data cell that is opposite to the
  *  neighbour to the cpu in the opposite direction */

  /* for each cell in bulk */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.y-1; ++j)
      for (k=1; k < cell_dim.z-1; ++k)

	/* For all of the neighbours of this cell */

	for (l=-1; l <= 1; ++l)
	  for (m=-1; m <= 1; ++m)
	    for (n=-1; n <= 1; ++n) {

	      /* Calculate Indicies of Neighbour */
	      neighbour.x = i + l;
	      neighbour.y = j + m;
	      neighbour.z = k + n;

	      /* Calculate Indices of Neighour CPU */
	      nbcoord.x = my_coord.x + l;
	      nbcoord.y = my_coord.y + m;
	      nbcoord.z = my_coord.z + n;

	      /* Calculate neighbours CPU via buffer cell and direction */
	      nbrank_cell = cpu_coord(global_cell_coord( neighbour ));
	      nbrank_dir = cpu_grid_coord( nbcoord );

	      /* Do send if neighbour is on different CPU */
	      if ((nbrank_cell != myid) && (nbrank_dir == nbrank_cell)) {

		/* calculate opposite of neighbour i.e. the cell to send */
		if      (nbcoord.x>my_coord.x) u=1;
		else if (nbcoord.x<my_coord.x) u=cell_dim.x-2;
		else                           u=neighbour.x;
                if      (nbcoord.y>my_coord.y) v=1;
		else if (nbcoord.y<my_coord.y) v=cell_dim.y-2;
		else                           v=neighbour.y;
		if      (nbcoord.z>my_coord.z) w=1;
		else if (nbcoord.z<my_coord.z) w=cell_dim.z-2;
		else                           w=neighbour.z;
                /* Exchange data */
		send_recv_cell(u, v, w, neighbour.x, neighbour.y, neighbour.z);
	      }
  }
}


/******************************************************************************
*
* send_cell_force lean version of send_cell for force_loop - unused
*
******************************************************************************/

void send_cell_force(cell *p, int to_cpu, int tag)
{
#if defined(PACX) || defined(MONOLJ)
  MPI_Send( &(p->n), 1,         MPI_INT, to_cpu, tag + SIZE_TAG,  cpugrid);
#endif
  MPI_Send( p->ort, DIM * p->n, REAL,    to_cpu, tag + ORT_TAG,   cpugrid);
#ifndef MONOLJ
  MPI_Send( p->sorte,     p->n, SHORT,   to_cpu, tag + SORTE_TAG, cpugrid);
#endif

}


/******************************************************************************
*
* recv_cell_force lean version of recv_cell for force_loop - unused
*
******************************************************************************/

void recv_cell_force(cell *p, int from_cpu,int tag)
{
  int size;
  int newsize;
  MPI_Status status;

#if defined(PACX) || defined(MONOLJ)
  MPI_Recv( &size, 1, MPI_INT, from_cpu, tag + SIZE_TAG , cpugrid, &status );
#else
  MPI_Probe( from_cpu, tag + ORT_TAG, cpugrid, &status );
  MPI_Get_count( &status, REAL, &size );
  size /= DIM;
#endif
    
  /* realloc cell if necessary */
  newsize = p->n_max; 
  while ( newsize < size ) newsize += incrsz;
  if (newsize > p->n_max) {
    /* Inihibit superfluous copy operation */
    p->n = 0;
    alloc_cell(p,newsize);
  }
  p->n = size;

  MPI_Recv(p->ort, DIM*size, REAL,from_cpu, tag + ORT_TAG,   cpugrid, &status);
#ifndef MONOLJ
  MPI_Recv(p->sorte, size, SHORT, from_cpu, tag + SORTE_TAG, cpugrid, &status);
#endif

}


/******************************************************************************
*
* Exchange the data of two cells
*
* This sends cell (i,j,k) into direction (i-l,j-m,k-n)
* and receives cell (l,m,n) from the opposite direction
*
******************************************************************************/

void send_recv_cell(int i, int j, int k, int l, int m, int n)
{
  int odd;
  int to_cpu,from_cpu;
  int tag;

  ivektor local;
  ivektor neighbour;
  ivektor opposite;

  cell *p,*q;

  if      (i!=l) odd = my_coord.x % 2;
  else if (j!=m) odd = my_coord.y % 2;
  else           odd = my_coord.z % 2;

  local.x = i;
  local.y = j;
  local.z = k;

  neighbour.x = l;
  neighbour.y = m;
  neighbour.z = n;

  if (neighbour.x < local.x)
    opposite.x = cell_dim.x - 1;
  else if (neighbour.x > local.x)
    opposite.x = 0;
  else
    opposite.x = local.x;

  if (neighbour.y < local.y)
    opposite.y = cell_dim.y - 1;
  else if (neighbour.y > local.y)
    opposite.y = 0;
  else
    opposite.y = local.y;

  if (neighbour.z < local.z)
    opposite.z = cell_dim.z - 1;
  else if (neighbour.z > local.z)
    opposite.z = 0;
  else
    opposite.z = local.z;

  /* Get neigbour CPU */
  to_cpu   = cpu_coord(global_cell_coord(opposite ));
  from_cpu = cpu_coord(global_cell_coord(neighbour));

  /* Calculate tag */

#ifndef MONOLJ
  tag = PTR_3D_VV(CELL_TAG,neighbour,cell_dim);
#else
  tag = ORT_TAG;
#endif

  /* DEBUG */
  if (myid == to_cpu) 
    error("Neighbour is on same cpu while sending.");

  if (myid == from_cpu) 
    error("Neighbour is on same cpu while receiving");

  /* Do ordered send/receive */
  p = PTR_3D_VV(cell_array, local    , cell_dim);
  q = PTR_3D_VV(cell_array, neighbour, cell_dim);

  if ( 0==odd ) {
    send_cell_force(p , to_cpu  , tag);
    recv_cell_force(q , from_cpu, tag);
  } else {
    recv_cell_force(q , from_cpu, tag);
    send_cell_force(p , to_cpu  , tag);
  };

}



