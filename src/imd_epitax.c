/******************************************************************************
*
* epitax.c - Routines for molecular beam epitaxy simulations in 2d and 3d
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

#ifndef TWOD

/******************************************************************************
*
* create_atom(type,mass,temp), 3d version
*
* creates a beam atom 
*
******************************************************************************/

void create_atom(int type, real mass, real temp)
{
  cell *p;
  ivektor cellc, neigh_cellc;
  vektor pos, dist;
  real ran1, ran2;
  int pbc_x, pbc_y;
  real min_dist2, dist2;
  int i, l, m, n, cpu;
  int j = 0;

  /* data for new atom */
  do {
    if (j==0)
      epitax_height += epitax_level * epitax_speed/(natoms - nepitax);
    else if (j==10000) error("EPITAX: 10000 search loops");
    ++j; 

    /* trial position of new particle (on master processor) */
    /* x,y coordinates at random, z coordinate shifted continuously */
    if (0 == myid) {
      ran1 = drand48();
      ran2 = drand48();
      pos.x = ran1 * box_x.x + ran2 * box_y.x;
      pos.y = ran1 * box_x.y + ran2 * box_y.y;
      pos.z = epitax_height;
    }

   /* broadcast position */ 
#ifdef MPI
    MPI_Bcast( &pos, DIM, REAL, 0, cpugrid);
#endif

    dist2  = box_z.z * box_z.z;
    dist.x = box_z.z;
    dist.y = box_z.z;
    dist.z = box_z.z;

    cellc = cell_coord(pos.x,pos.y,pos.z);

    /* minimal distance from atoms in neighboring cells */
    for ( l=-1; l<2; ++l)
      for ( m=-1; m<2; ++m)
	for ( n=-1; n<2; ++n) {

	  /* Neighbor cells */
	  neigh_cellc.x = cellc.x + l;
	  neigh_cellc.y = cellc.y + m;
	  neigh_cellc.z = cellc.z + n;

	  /* Apply periodic boundary conditions */
	  pbc_x = 0; pbc_y = 0;
	  if ( neigh_cellc.x < 0 ) {
	    neigh_cellc.x += global_cell_dim.x;    pbc_x = 1;
	  }
	  else if ( neigh_cellc.x >= global_cell_dim.x ) {
	    neigh_cellc.x -= global_cell_dim.x;    pbc_x = -1;
	  }
	  if ( neigh_cellc.y < 0 ) {
	    neigh_cellc.y += global_cell_dim.y;    pbc_y = 1;
	  }
	  else if ( neigh_cellc.y >= global_cell_dim.y ) {
	    neigh_cellc.y -= global_cell_dim.y;    pbc_y = -1;
	  }

#ifdef MPI
	  cpu = cpu_coord(neigh_cellc);
	  if (cpu == myid) {
	    neigh_cellc = cell_map(neigh_cellc);
#endif
	    p = PTR_3D_VV( cell_array, neigh_cellc, cell_dim);
	  
	    for ( i=0; i<p->n; ++i) {
	      
	      dist.x = pos.x - ORT(p,i,X) + pbc_x * box_x.x;
	      dist.y = pos.y - ORT(p,i,Y) + pbc_y * box_y.y;
	      dist.z = pos.z - ORT(p,i,Z);

	      dist2  = MIN(dist2,SPROD(dist,dist));
	    }
#ifdef MPI
	  }
#endif
	} /* for l,m,n */
#ifdef MPI
    MPI_Allreduce( &dist2, &min_dist2, 1, REAL, MPI_MIN, cpugrid); 
#else
    min_dist2 = dist2;
#endif

    /* check if distance is smaller than cutoff radius */
  } while (min_dist2 < SQR(epitax_cutoff));


  /* create new atom */

#ifdef MPI
  cpu = cpu_coord(cellc);
  if (cpu == myid) {
    cellc = cell_map(cellc);
#endif

  p = PTR_3D_VV(cell_array, cellc, cell_dim);

  if ( p->n >= p->n_max) alloc_cell(p, p->n_max + CSTEP);

  ORT(p,p->n,X) = pos.x;
  ORT(p,p->n,Y) = pos.y;
  ORT(p,p->n,Z) = pos.z;

  IMPULS(p,p->n,X) = 0.0;
  IMPULS(p,p->n,Y) = 0.0;
  IMPULS(p,p->n,Z) = - sqrt( 3 * temp * mass );

  KRAFT(p,p->n,X) = 0.0;
  KRAFT(p,p->n,Y) = 0.0;
  KRAFT(p,p->n,Z) = 0.0;
  
  NUMMER(p,p->n)  = epitax_number + 1;
  VSORTE(p,p->n)  = type;
  MASSE(p,p->n)   = mass;
  POTENG(p,p->n)  = 0.0;

  ++p->n;

#ifndef MPI
  printf("Atom %d created, coordinates x=%lf, y=%lf, z=%lf\n", epitax_number + 1, pos.x, pos.y, pos.z);
#endif
#ifdef MPI
  }
#endif
 
  ++natoms;  
  ++nactive;
  ++nepitax;
  ++epitax_number;

}

/******************************************************************************
*
* delete_atoms, 3d version
*
* deletes atoms that leave the simulation box in z-direction
*
******************************************************************************/

void delete_atoms(void)
{
  cell *p;
  ivektor cellc;
  int i, l, m, cpu;
  int n_glob, n_loc = 0; 
  int nepitax_glob, nepitax_loc = 0;

#ifdef MPI
  if ( (cpu_dim.z == 1) || ((myid+1)%cpu_dim.z == 0) ) {
    cellc.z = cell_dim.z - 2;
    for ( l = 0; l < (cell_dim.x - 2); ++l)
      for ( m = 0; m < (cell_dim.y - 2); ++m) { 
	cellc.x = l + 1;
	cellc.y = m + 1;
#else
  cellc.z = cell_dim.z - 1;  
  for ( l = 0; l < cell_dim.x; ++l)
    for ( m = 0; m < cell_dim.y; ++m) { 
      cellc.x = l;
      cellc.y = m;
#endif

	/* atoms in top cells are deleted */
	p = PTR_3D_VV( cell_array, cellc, cell_dim);
     
	if ( p->n > 0) {

	  for ( i=0; i<p->n; ++i) {
#ifndef MPI
	    printf("Atom %d deleted. z-coordinate = %lf\n", 
                   NUMMER(p,i), ORT(p,i,Z) );
#endif
	    if ( NUMMER(p,i) > epitax_sub_n)  ++nepitax_loc;
	    ++n_loc;
	  }

	  p->n = 0;
	
	} /* if */
      }  /* for */

#ifdef MPI
  }
  MPI_Allreduce( &nepitax_loc, &nepitax_glob, 1, MPI_INT, MPI_SUM, cpugrid);
  nepitax  -= nepitax_glob;
  MPI_Allreduce( &n_loc, &n_glob, 1, MPI_INT, MPI_SUM, cpugrid);
  natoms   -= n_glob;
  nactive  -= n_glob;
#else
  nepitax  -= nepitax_loc;
  natoms   -= n_loc;
  nactive  -= n_loc;
#endif

}

/******************************************************************************
*
*  substrate_level, 3d version 
*
*  computes z-coordinate of substrate
*
******************************************************************************/

real substrate_level(void)
{
  cell *p;
  int i, k;
  real level, tmp_level = 0.0;

  for (k=0; k<ncells; ++k ) {
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) tmp_level = MAX(tmp_level, ORT(p,i,Z) );
  }

#ifdef MPI
  MPI_Allreduce( &tmp_level, &level, 1, REAL, MPI_MAX, cpugrid);
#else
  level = tmp_level;
#endif

  return level;
}

/******************************************************************************
*
*  calc_poteng_min, 3d version 
*
*  computes minimal potential energy per atom
*
******************************************************************************/

void calc_poteng_min(void)
{
  cell *p;
  int i, k;
  real tmp_epot = 0.0;

  for (k=0; k<ncells; ++k ) {
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i)  
      tmp_epot = MIN(tmp_epot, POTENG(p,i) );
  }

#ifdef MPI
  MPI_Allreduce( &tmp_epot, &epitax_poteng_min, 1, REAL, MPI_MIN, cpugrid);
#else
  epitax_poteng_min = tmp_epot;
#endif
  
}

/******************************************************************************
*
*  check_boxheight, 3d version 
*
*  checks whether boxheight is large enough
*
******************************************************************************/

void check_boxheight(void)
{
  real epitax_nmax = 0.0;
  int i;

  for ( i=0; i<ntypes; ++i)
    if ( epitax_rate[i] != 0)
      epitax_nmax += ( epitax_maxsteps - epitax_startstep ) / epitax_rate[i];

  if ( ( epitax_height + epitax_level * epitax_nmax / natoms 
       + box_z.z / ( cell_dim.z + 1 ) ) > box_z.z)
    error("EPITAX: Boxheight too small. Increase box_z.z!");
    
}


/******************************************************************************
*
*  cell_map(cellc), 3d version (MPI) 
*
*  computes local cell coordinates from global cell coordinates 
*
******************************************************************************/

ivektor cell_map(ivektor cellc)
{
  cellc.x = cellc.x % ( cell_dim.x - 2 ) + 1;
  cellc.y = cellc.y % ( cell_dim.y - 2 ) + 1;
  cellc.z = cellc.z % ( cell_dim.z - 2 ) + 1;

  return cellc;
}


#else 

/******************************************************************************
*
* create_atom(type,mass,temp), 2d version
*
* creates a beam atom
*
******************************************************************************/

void create_atom(int type, real mass, real temp)
{
  cell *p;
  ivektor cellc, neigh_cellc;
  vektor pos, dist;
  int pbc_x;
  real min_dist2, dist2;
  int i, l, m, cpu; 
  int j = 0;

  /* data for new atom */
  do {
    if (j==0)
      epitax_height += epitax_level * epitax_speed/(natoms - nepitax);
    else if (j==10000) error("EPITAX: 10000 search loops");
    ++j; 

    /* trial position of new particle (on master processor) */
    /* x coordinate at random, y coordinate shifted continuously */
    if (0 == myid) {
    pos.x = drand48() * box_x.x;
    pos.y = epitax_height;
    }

   /* broadcast position */ 
#ifdef MPI
    MPI_Bcast( &pos, DIM, REAL, 0, cpugrid);
#endif

    dist2  = box_y.y * box_y.y;
    dist.x = box_y.y;
    dist.y = box_y.y;

    cellc = cell_coord(pos.x,pos.y);

    /* minimal distance from atoms in neighboring cells */
    for ( l=-1; l<2; ++l)
      for ( m=-1; m<2; ++m) {

	  /* Neighbor cells */
	  neigh_cellc.x = cellc.x + l;
	  neigh_cellc.y = cellc.y + m;

	  /* Apply periodic boundary conditions */
	  pbc_x = 0;
	  if ( neigh_cellc.x < 0 ) {
	    neigh_cellc.x += global_cell_dim.x;    pbc_x = 1;
	  }
	  else if ( neigh_cellc.x >= global_cell_dim.x ) {
	    neigh_cellc.x -= global_cell_dim.x;    pbc_x = -1;
	  }

#ifdef MPI
	  cpu = cpu_coord(neigh_cellc);
	  if (cpu == myid) {
	    neigh_cellc = cell_map(neigh_cellc);
#endif
	    p = PTR_2D_VV( cell_array, neigh_cellc, cell_dim);
	  
	    for ( i=0; i<p->n; ++i) {
	      
	      dist.x = pos.x - ORT(p,i,X) + pbc_x * box_x.x;
	      dist.y = pos.y - ORT(p,i,Y);

	      dist2  = MIN(dist2,SPROD(dist,dist));
	    }
#ifdef MPI
	  }
#endif
	} /* for l,m,n */
#ifdef MPI
    MPI_Allreduce( &dist2, &min_dist2, 1, REAL, MPI_MIN, cpugrid); 
#else
    min_dist2 = dist2;
#endif

    /* check if distance is smaller than cutoff radius */
  } while (min_dist2 < r2_cut);


  /* create new atom */

#ifdef MPI
  cpu = cpu_coord(cellc);
  if (cpu == myid) {
    cellc = cell_map(cellc);
#endif

  p = PTR_2D_VV(cell_array, cellc, cell_dim);

  if ( p->n >= p->n_max) alloc_cell(p, p->n_max + CSTEP);

  ORT(p,p->n,X)    = pos.x;
  ORT(p,p->n,Y)    = pos.y;

  IMPULS(p,p->n,X) = 0.0;
  IMPULS(p,p->n,Y) = - sqrt( 3 * temp * mass );

  KRAFT(p,p->n,X)  = 0.0;
  KRAFT(p,p->n,Y)  = 0.0;
  
  NUMMER(p,p->n)   = epitax_number + 1;
  VSORTE(p,p->n)   = type;
  MASSE(p,p->n)    = mass;
  POTENG(p,p->n)   = 0.0;

  ++p->n;

#ifndef MPI
  printf("Atom %d created, coordinates x=%lf, y=%lf\n", 
         epitax_number + 1, pos.x, pos.y);
#endif
#ifdef MPI
  }
#endif
 
  ++natoms;  
  ++nactive;
  ++nepitax;
  ++epitax_number;

}

/******************************************************************************
*
* delete_atoms, 2d version 
*
* deletes atoms that leave the simulation box in z-direction
*
******************************************************************************/

void delete_atoms(void)
{
  cell *p;
  ivektor cellc;
  int i, l, cpu;
  int n_glob, n_loc = 0; 
  int nepitax_glob, nepitax_loc = 0;

#ifdef MPI
  if ( (cpu_dim.y == 1) || ((myid+1)%cpu_dim.y == 0) ) {
    cellc.y = cell_dim.y - 2;
    for ( l = 0; l < (cell_dim.x - 2); ++l) {
      cellc.x = l + 1;
#else
  cellc.y = cell_dim.y - 1;    
  for ( l = 0; l < cell_dim.x; ++l) {  
    cellc.x = l;
#endif

      /* atoms in top cells are deleted */
      p = PTR_2D_VV( cell_array, cellc, cell_dim);
      
      if ( p->n > 0) {

	for ( i=0; i<p->n; ++i) {
#ifndef MPI
	  printf("Atom %d deleted. y-coordinate = %lf\n", 
                 NUMMER(p,i), ORT(p,i,Y) );
#endif
	  if ( NUMMER(p,i) > epitax_sub_n)  ++nepitax_loc;
	  ++n_loc;
	}

	p->n = 0;
	
      }  /* if */
    }   /* for */

#ifdef MPI
  }
  MPI_Allreduce( &nepitax_loc, &nepitax_glob, 1, MPI_INT, MPI_SUM, cpugrid);
  nepitax  -= nepitax_glob;
  MPI_Allreduce( &n_loc, &n_glob, 1, MPI_INT, MPI_SUM, cpugrid);
  natoms   -= n_glob;
  nactive  -= n_glob;
#else
  nepitax  -= nepitax_loc;
  natoms   -= n_loc;
  nactive  -= n_loc;
#endif

}

/******************************************************************************
*
*  substrate_level, 2d version 
*
*  computes y-coordinate of substrate
*
******************************************************************************/

real substrate_level(void)
{
  cell *p;
  int i, k;
  real level, tmp_level = 0.0;

  for (k=0; k<ncells; ++k ) {
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) tmp_level = MAX(tmp_level, ORT(p,i,Y) );
  }

#ifdef MPI
  MPI_Allreduce( &tmp_level, &level, 1, REAL, MPI_MAX, cpugrid);
#else
  level = tmp_level;
#endif

  return level;
}

/******************************************************************************
*
*  calc_poteng_min, 2d version 
*
*  computes minimal potential energy per atom
*
******************************************************************************/

void calc_poteng_min(void)
{
  cell *p;
  int i, k;
  real tmp_epot = 0.0;

  for (k=0; k<ncells; ++k ) {
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i)  
      tmp_epot = MIN(tmp_epot, POTENG(p,i) );
  }

#ifdef MPI
  MPI_Allreduce( &tmp_epot, &epitax_poteng_min, 1, REAL, MPI_MIN, cpugrid);
#else
  epitax_poteng_min = tmp_epot;
#endif
  
}

/******************************************************************************
*
*  check_boxheight, 2d version 
*
*  checks whether boxheight is large enough
*
******************************************************************************/

void check_boxheight(void)
{
  real epitax_nmax = 0.0;
  int i;

  for ( i=0; i<ntypes; ++i)
    if ( epitax_rate[i] != 0)
      epitax_nmax += ( epitax_maxsteps - epitax_startstep ) / epitax_rate[i];

  if ( ( epitax_height + epitax_level * epitax_nmax / natoms 
       + box_y.y / ( cell_dim.y + 1 ) ) > box_y.y)
    error("EPITAX: Boxheight too small. Increase box_y.y!");
    
}


/******************************************************************************
*
*  cell_map(cellc), 2d version (MPI) 
*
*  computes local cell coordinates from global cell coordinates
*
******************************************************************************/

ivektor cell_map(ivektor cellc)
{
  cellc.x = cellc.x % ( cell_dim.x - 2 ) + 1;
  cellc.y = cellc.y % ( cell_dim.y - 2 ) + 1;

  return cellc;
}

#endif



