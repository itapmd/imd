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
* imd_ttm.c -- Code for two-temperature model calculations
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"
#ifdef DEBUG
#include <assert.h>
#endif

#ifdef BUFCELLS
#define NBUFFC 2
#else
#define NBUFFC 0
#endif /*BUFCELLS*/

/* Macro to allow electr. heat capacity to vary
 * with electronic temperature (C_e=gamma*T_e) */
#define FD_C ((fd_c==0)?(fd_gamma*l1[i][j][k].temp):(fd_c))

/* update_fd(): update natoms_local, fd_min_atoms and natoms,
 * md_temp and v_com in FD lattice cells 
 * watch out for activated or deactivated cells */
void update_fd()
{
  int i,j,k;

  natoms_local=0;
  for (i=1; i<local_fd_dim.x-1; ++i)
  { 
    for (j=1; j<local_fd_dim.y-1; ++j)
    {
      for (k=1; k<local_fd_dim.z-1; ++k)
      {
	int loopvar;
	int natoms_previous;
	real tot_mass=0; /* mass of all atoms in the FD cell */
	int l;
	cellptr p;

	l1[i][j][k].v_com.x=0.0;
	l1[i][j][k].v_com.y=0.0;
	l1[i][j][k].v_com.z=0.0;
	natoms_previous=l1[i][j][k].natoms;
	l1[i][j][k].natoms=0;

	/* loop over encompassed MD cells */
	for (loopvar=0; loopvar<n_md_cells; loopvar++)
	{
	  p=l1[i][j][k].md_cellptrs[loopvar];

	  /* add number of atoms of this MD cell*/
	  l1[i][j][k].natoms += p->n;
	  natoms_local += p->n;

	  /* next we need the velocity of the center of mass of the FD cell.
	   * For this, sum up momenta and masses of particles seperately.
	   */

	  /* loop over atoms in MD cell */
	  for (l=0; l<p->n; l++)
	  {
	    l1[i][j][k].v_com.x += IMPULS(p,l,X);
	    l1[i][j][k].v_com.y += IMPULS(p,l,Y);
	    l1[i][j][k].v_com.z += IMPULS(p,l,Z);
	    tot_mass += MASSE(p,l);
	  } 

	}


	l2[i][j][k].natoms  = l1[i][j][k].natoms;

	if (tot_mass!=0) /* tot_mass==0 can happen when no atoms in cell */
	{
	  l1[i][j][k].v_com.x/=tot_mass;
	  l1[i][j][k].v_com.y/=tot_mass;
	  l1[i][j][k].v_com.z/=tot_mass;
	  l2[i][j][k].v_com=l1[i][j][k].v_com;
	}

	/* next, MD temperature of cell, from kinetic energies of particles
	 * in the reference frame of the center of mass */
	l1[i][j][k].md_temp=0;

	/* loop over MD cells and atoms */
	for (loopvar=0; loopvar<n_md_cells; loopvar++)
	{
	  p=l1[i][j][k].md_cellptrs[loopvar];
	  for (l=0; l<p->n; l++)
	  {
	    l1[i][j][k].md_temp += MASSE(p,l)*SQR(IMPULS(p,l,X)/MASSE(p,l)
		                   - l1[i][j][k].v_com.x);
	    l1[i][j][k].md_temp += MASSE(p,l)*SQR(IMPULS(p,l,Y)/MASSE(p,l)
		                   - l1[i][j][k].v_com.y);
	    l1[i][j][k].md_temp += MASSE(p,l)*SQR(IMPULS(p,l,Z)/MASSE(p,l)
		                   - l1[i][j][k].v_com.z);
	  }
	}
	if (l1[i][j][k].natoms != 0)
	{
	  l1[i][j][k].md_temp /= 3*l1[i][j][k].natoms;
	  /* E_kin per particle=3/2*T (assuming here all particles
	   * are unrestricted 3D)                                    */
	}
	l2[i][j][k].md_temp = l1[i][j][k].md_temp;

        if (natoms_previous>=fd_min_atoms && l1[i][j][k].natoms<fd_min_atoms)
	{
#ifdef DEBUG
	  warning("FD cell deactivated\n");
#endif
          /* Cell deactivated. Deduce its electronic energy from E_new_local */
	  l1[i][j][k].xi=0.0;
	  E_new_local -= (fd_c==0)?
	                 (0.5*fd_gamma*SQR(l1[i][j][k].temp)):
			 (fd_c*l1[i][j][k].temp);
	}

	if (natoms_previous<fd_min_atoms && l1[i][j][k].natoms>=fd_min_atoms)
	{
#ifdef DEBUG
	  warning("New FD cell activated\n");
#endif
	  /* Freshly activated cell. Gets avg. electron energy of active 
             neighbor cells, the created energy is added to E_new_local */
	  int n_neighbors=0;
	  double E_el_neighbors=0.0;
	  /* 6 indices: -x,x,-y,y,-z,z */

	  if (l1[i+1][j][k].natoms >= fd_min_atoms)
	  {
	    E_el_neighbors += (fd_c==0)?
	                      (SQR(l1[i+1][j][k].temp)):
			      (l1[i+1][j][k].temp);
	    n_neighbors++;
	  }
	  if (l1[i-1][j][k].natoms >= fd_min_atoms)
	  {
	    E_el_neighbors += (fd_c==0)?
	                      (SQR(l1[i-1][j][k].temp)):
			      (l1[i-1][j][k].temp);
	    n_neighbors++;
	  }
	  if (l1[i][j+1][k].natoms >= fd_min_atoms)
	  {
	    E_el_neighbors += (fd_c==0)?
	                      (SQR(l1[i][j+1][k].temp)):
			      (l1[i][j+1][k].temp);
	    n_neighbors++;
	  }
	  if (l1[i][j-1][k].natoms >= fd_min_atoms)
	  {
	    E_el_neighbors += (fd_c==0)?
	                      (SQR(l1[i][j-1][k].temp)):
			      (l1[i][j-1][k].temp);
	    n_neighbors++;
	  }
	  if (l1[i][j][k+1].natoms >= fd_min_atoms)
	  {
	    E_el_neighbors += (fd_c==0)?
	                      (SQR(l1[i][j][k+1].temp)):
			      (l1[i][j][k+1].temp);
	    n_neighbors++;
	  }
	  if (l1[i][j][k-1].natoms >= fd_min_atoms)
	  {
	    E_el_neighbors +=(fd_c==0)?
	                     (SQR(l1[i][j][k-1].temp)):
			     (l1[i][j][k-1].temp);
	    n_neighbors++;
	  }

	  E_new_local += (fd_c==0)?
	                 (0.5*fd_gamma*SQR(l1[i][j][k].temp)):
			 (fd_c*l1[i][j][k].temp);

	  if (n_neighbors != 0)
	  {
	    l1[i][j][k].temp = (fd_c==0)?
	                       (sqrt(E_el_neighbors/n_neighbors)):
			       (E_el_neighbors/n_neighbors);
	  }
	  else{
	    l1[i][j][k].temp=l1[i][j][k].md_temp;
	  }
	  E_new_local += (fd_c==0)?
	                 (0.5*fd_gamma*SQR(l1[i][j][k].temp)):
			 (fd_c*l1[i][j][k].temp);

	}
      }
    }
  }

  fd_min_atoms = natoms_local /
                 ((local_fd_dim.x-2)*(local_fd_dim.y-2)*(local_fd_dim.z-2)*10);
  fd_min_atoms = MAX(fd_min_atoms,10);
#ifdef MPI
  {
    double E_new_reduced;
    MPI_Reduce(&E_new_local,&E_new_reduced,1,MPI_DOUBLE,MPI_SUM,0,cpugrid);
    if (myid==0) E_new += E_new_reduced * fd_h.x*fd_h.y*fd_h.z / natoms;
  }
#else
  E_new += E_new_local * fd_h.x*fd_h.y*fd_h.z / natoms;
#endif /*MPI*/

  E_new_local=0.0;

}


#define TARGET_TEMP (init_t_el==0?l1[i][j][k].md_temp:init_t_el)

/** Function to set electron temperature **/
void ttm_overwrite()
{
  int i,j,k;

  for (i=1; i<local_fd_dim.x-1; ++i)
  {
    for (j=1; j<local_fd_dim.y-1; ++j)
    {
      for (k=1; k<local_fd_dim.z-1; ++k)
	l1[i][j][k].temp = l2[i][j][k].temp = TARGET_TEMP;
    }
  }
}

/* init_ttm(): initialize FD stuff */
void init_ttm()
{
  int i,j,k;

  natoms_local=0;
  n_md_cells=fd_ext.x*fd_ext.y*fd_ext.z;

  /* Check if cell_dim and fd_ext are commensurate */
  if ( fd_one_d == 1 || fd_one_d == 0 )
    if ( (cell_dim.x-NBUFFC) % fd_ext.x != 0 )
    {
      char buf[255];
      sprintf(buf,
     "cell_dim and fd_ext are not commensurate:\ncell_dim.x=%d; fd_ext.x=%d\n",
              cell_dim.x, fd_ext.x );
      error (buf);
    }
  if ( fd_one_d == 2 || fd_one_d == 0 )
    if ( (cell_dim.y-NBUFFC) % fd_ext.y != 0 )
    {
      char buf[255];
      sprintf(buf,
     "cell_dim and fd_ext are not commensurate:\ncell_dim.y=%d; fd_ext.y=%d\n",
              cell_dim.y, fd_ext.y );
      error (buf);
    }
  if ( fd_one_d == 3 || fd_one_d == 0 )
    if ( (cell_dim.z-NBUFFC) % fd_ext.z != 0 )
    {
      char buf[255];
      sprintf(buf,
     "cell_dim and fd_ext are not commensurate:\ncell_dim.z=%d; fd_ext.z=%d\n",
              cell_dim.z, fd_ext.z );
      error (buf);
    }

  /* local size of FD lattice, will add 2 layers later for boundaries */
  local_fd_dim.x=(cell_dim.x-NBUFFC)/fd_ext.x;
  local_fd_dim.y=(cell_dim.y-NBUFFC)/fd_ext.y;
  local_fd_dim.z=(cell_dim.z-NBUFFC)/fd_ext.z;

  /* check if we want a 1D TTM simulation (this probably needs some work)*/
  if (fd_one_d ==1) 
  {
    if (cpu_dim.y!=1 || cpu_dim.z!=1)
      error("cpu_dim.y and cpu_dim.z must be 1 for fd_one_d==x\n");
    local_fd_dim.y=local_fd_dim.z = 1;
  }
  if (fd_one_d ==2)
  {
    local_fd_dim.x=local_fd_dim.z = 1;
    if (cpu_dim.x!=1 || cpu_dim.z!=1)
      error("cpu_dim.x and cpu_dim.z must be 1 for fd_one_d==y\n");
  }
  if (fd_one_d ==3)
  {
    local_fd_dim.x=local_fd_dim.y = 1;
    if (cpu_dim.x!=1 || cpu_dim.y!=1)
      error("cpu_dim.x and cpu_dim.y must be 1 for fd_one_d==z\n");
  }

  /* get these to the right sizes */
  global_fd_dim.x = local_fd_dim.x * cpu_dim.x; 
  global_fd_dim.y = local_fd_dim.y * cpu_dim.y;
  global_fd_dim.z = local_fd_dim.z * cpu_dim.z;
  local_fd_dim.x += 2; /* 2 for ghost layers */
  local_fd_dim.y += 2; /* 2 for ghost layers */
  local_fd_dim.z += 2; /* 2 for ghost layers */

  /* Time to initialize our FD lattice... */

  /* Allocate memory for two lattices */
  lattice1=(ttm_Element*) malloc(
           (local_fd_dim.x)*(local_fd_dim.y)*(local_fd_dim.z)*sizeof(ttm_Element));
  lattice2=(ttm_Element*) malloc(
           (local_fd_dim.x)*(local_fd_dim.y)*(local_fd_dim.z)*sizeof(ttm_Element));
  l1 = (ttm_Element***) malloc( local_fd_dim.x * sizeof(ttm_Element**) );
  l2 = (ttm_Element***) malloc( local_fd_dim.x * sizeof(ttm_Element**) );
  for (i=0; i<local_fd_dim.x; i++)
  {
    l1[i] = (ttm_Element**) malloc( local_fd_dim.y * sizeof(ttm_Element*) );
    l2[i] = (ttm_Element**) malloc( local_fd_dim.y * sizeof(ttm_Element*) );

    for (j=0; j<local_fd_dim.y; j++)
    {
      l1[i][j] = lattice1 + i*local_fd_dim.y*local_fd_dim.z + j*local_fd_dim.z;
      l2[i][j] = lattice2 + i*local_fd_dim.y*local_fd_dim.z + j*local_fd_dim.z;

      if (i!=0 && j!=0 && i!=local_fd_dim.x-1 && j!=local_fd_dim.y-1)
      {	/* we are not in a ghost layer*/
	for (k=1; k<local_fd_dim.z-1; k++)
	{ 
	  int tmpindex=0;
	  int xc,yc,zc;
	  int loopvar;

	  /* Initialize md_cellptrs, temp,
	   * and source of this FE Cell.
	   * Also set fd_cell_idx in encompassed MD cells.
	   **********/

	  /* allocate MD-cell pointer arrays */
	  l1[i][j][k].md_cellptrs=(cellptr*)malloc(n_md_cells*sizeof(cellptr));
	  l2[i][j][k].md_cellptrs=(cellptr*)malloc(n_md_cells*sizeof(cellptr));

	  /* loop over encompassed MD cells */
	  for (xc=(i-1)*fd_ext.x+NBUFFC/2; xc<i*fd_ext.x+NBUFFC/2; xc++)
	  {
	    for (yc=(j-1)*fd_ext.y+NBUFFC/2; yc<j*fd_ext.y+NBUFFC/2; yc++)
	    {
	      for (zc=(k-1)*fd_ext.z+NBUFFC/2; zc<k*fd_ext.z+NBUFFC/2; zc++)
	      {
		cellptr p;
		
		/* pointer to this MD cell */
		p = l1[i][j][k].md_cellptrs[tmpindex] = 
		    l2[i][j][k].md_cellptrs[tmpindex] =
		    PTR_3D_V(cell_array,xc,yc,zc,cell_dim);

		/* write array indices of our FD cell to this MD cell*/
		p->fd_cell_idx.x = i;
		p->fd_cell_idx.y = j;
		p->fd_cell_idx.z = k;

		natoms_local += p->n;
		tmpindex++;
	      }
	    }
	  }

#ifdef DEBUG
	  assert(tmpindex==n_md_cells);
#endif

	  /* no incoming thermal power per default */
	  l2[i][j][k].source = l1[i][j][k].source = 0.0;

	}
      }
    }
  }

  if (myid==0)
  {
    printf("Global FD cell array dimensions: %d %d %d\n",
	global_fd_dim.x, global_fd_dim.y, global_fd_dim.z);
    printf("Local FD cell array dimensions: %d %d %d\n",
	local_fd_dim.x, local_fd_dim.y, local_fd_dim.z);

    printf("Volume of one FD cell: %e [cubic Angstroms]\n",
           fd_h.x*fd_h.y*fd_h.z );
    printf("Volume of whole sample: %e [cubic Angstroms]\n",
           fd_h.x*fd_h.y*fd_h.z *
           global_fd_dim.x*global_fd_dim.y*global_fd_dim.z );
  }
#ifdef DEBUG
    printf(
      "Found %d atoms initializing FD lattice in process number %d.\n",
           natoms_local, myid );
#endif

  update_fd(); /* get md_temp and v_com etc. */

  ttm_overwrite(); /* electron temperature is initialized */


#ifdef MPI
  /* create MPI datatypes */
  ttm_create_mpi_datatypes();

  /* time to contact neighbors and fill ghost layers */
#endif /*MPI*/
  ttm_fill_ghost_layers();

  if (fix_t_el!=0) /* T_el will remain fixed from now on.
		      Calculate value ttm_eng,
		      which will remain fixed, too */
  {
    /* temp * capacity * volume / natoms */
    ttm_eng = init_t_el * 
              global_fd_dim.x * global_fd_dim.y * global_fd_dim.z *
              ((fd_c==0)?(0.5*fd_gamma*init_t_el):(fd_c)) * 
              fd_h.x * fd_h.y * fd_h.z / natoms;
  }

  if (myid==0) printf( "Using Two Temperature Model TTM\n");

}

/* solve heat diffusion equation for electronic system */
void calc_ttm()
{
  int i,j,k;
  int fd_timestep;
  int xmin, xmax, /* these are neighboring indices           */
      ymin, ymax, /* (to account for bc & deactivated cells) */
      zmin, zmax;

  if(fix_t_el==0) /* T_el is not fixed, otherwise no big calculations needed */
  {
    if (steps%fd_update_steps==0)
    { /* we need new lattice temperature and number of atoms etc. */
      update_fd();
    }

    /* set all xi to zero */
    for (i=1; i<local_fd_dim.x-1; ++i)
    {
      for (j=1; j<local_fd_dim.y-1; ++j)
      {
	for (k=1; k<local_fd_dim.z-1; ++k)
	  l1[i][j][k].xi = l2[i][j][k].xi = 0.0;
      }
    }

#ifdef DEBUG
    E_el_ab_local = 0.0;
#endif

    for (fd_timestep=1; fd_timestep<=fd_n_timesteps; ++fd_timestep)
    {

      for (i=1; i<local_fd_dim.x-1; ++i)
      {
	for (j=1; j<local_fd_dim.y-1; ++j)
	{
	  for (k=1; k<local_fd_dim.z-1; ++k)
	  {
	    /* only do calculation if cell is not deactivated */
	    if (l1[i][j][k].natoms < fd_min_atoms) 
	    {
	      continue;
	    }

	    if (l1[i-1][j][k].natoms < fd_min_atoms)
	      xmin=i;
	    else
	      xmin=i-1;

	    if (l1[i+1][j][k].natoms < fd_min_atoms)
	      xmax=i;
	    else
	      xmax=i+1;

	    if (l1[i][j-1][k].natoms < fd_min_atoms)
	      ymin=j;
	    else
	      ymin=j-1;

	    if (l1[i][j+1][k].natoms < fd_min_atoms)
	      ymax=j;
	    else
	      ymax=j+1;

	    if (l1[i][j][k-1].natoms < fd_min_atoms)
	      zmin=k;
	    else
	      zmin=k-1;

	    if (l1[i][j][k+1].natoms < fd_min_atoms)
	      zmax=k;
	    else
	      zmax=k+1;

#ifdef DEBUG
	    E_el_ab_local += l1[i][j][k].temp - l1[i][j][k].md_temp ;
#endif

	    /* NOW calculate */

	    l2[i][j][k].temp = timestep/fd_n_timesteps *
	      ( fd_k/FD_C *
		(   1.0/(fd_h.x * fd_h.x) * ( l1[xmin][j][k].temp + l1[xmax][j][k].temp - 2*l1[i][j][k].temp )
		    + 1.0/(fd_h.y * fd_h.y) * ( l1[i][ymin][k].temp + l1[i][ymax][k].temp - 2*l1[i][j][k].temp )
		    + 1.0/(fd_h.z * fd_h.z) * ( l1[i][j][zmin].temp + l1[i][j][zmax].temp - 2*l1[i][j][k].temp ) )
		- 1.0/FD_C * fd_g * ( l1[i][j][k].temp - l1[i][j][k].md_temp )
		+ 1.0/FD_C * l1[i][j][k].source )
	      + l1[i][j][k].temp;


	    l2[i][j][k].xi += (l2[i][j][k].temp-l2[i][j][k].md_temp);
	    l1[i][j][k].xi = l2[i][j][k].xi;
	  }
	}
      }


      /* take care - l1 must always be the updated lattice */
      l3=l1;
      l1=l2;
      l2=l3;

      /* MPI communication / pbc / reflecting bc */
      ttm_fill_ghost_layers();

    }

    ttm_eng=0.0;

    /* summed xi still need a factor,
     * and we update ttm_eng */
    for (i=1; i<local_fd_dim.x-1; ++i)
    {
      for (j=1; j<local_fd_dim.y-1; ++j)
      {
	for (k=1; k<local_fd_dim.z-1; ++k)
	{
	  if(l1[i][j][k].natoms>=fd_min_atoms)
	  {
	    l1[i][j][k].xi *= fd_g * fd_h.x*fd_h.y*fd_h.z / 
	      (fd_n_timesteps * l1[i][j][k].md_temp * 3 * l1[i][j][k].natoms);
	  } else 
	  {
	    l1[i][j][k].xi = 0.0;
	  }
	  l2[i][j][k].xi = l1[i][j][k].xi;
	  /* E=\gamma/2*T^2 or E=c_e*T?*/
	  ttm_eng += (fd_c==0)?
	             (0.5*fd_gamma*SQR(l1[i][j][k].temp)):
		     (fd_c*l1[i][j][k].temp);
	}
      }
    }

#ifdef DEBUG
#ifdef MPI
    {
      double E_el_ab_reduced;
      MPI_Reduce(&E_el_ab_local,&E_el_ab_reduced, 1, MPI_DOUBLE, MPI_SUM, 0, cpugrid);
      if (myid==0)
      {
	E_el_ab += E_el_ab_reduced *  fd_h.x*fd_h.y*fd_h.z *
	           fd_g * timestep / fd_n_timesteps ;
      }
    }
#else
    E_el_ab += E_el_ab_local * fd_h.x*fd_h.y*fd_h.z *
               fd_g * timestep / fd_n_timesteps ;
#endif /*MPI*/
#endif /*DEBUG*/

    ttm_eng *= fd_h.x*fd_h.y*fd_h.z;
#ifdef MPI
    {
      double ttm_eng_reduced;
      MPI_Reduce(&ttm_eng, &ttm_eng_reduced, 1, MPI_DOUBLE, MPI_SUM, 0, cpugrid);
      if(myid==0) ttm_eng = ttm_eng_reduced / natoms;
    }
#else
    ttm_eng /= natoms;
#endif /*MPI*/
  } else
  { /* T_el is fixed. But we need xi */
    for (i=1; i<local_fd_dim.x-1; ++i)
    {
      for (j=1; j<local_fd_dim.y-1; ++j)
      {
	for (k=1; k<local_fd_dim.z-1; ++k)
	{
	  if(l1[i][j][k].natoms != 0)
	  {
	    l1[i][j][k].xi = fd_g * fd_h.x*fd_h.y*fd_h.z *
	      (l1[i][j][k].temp-l1[i][j][k].md_temp) / 
	      (l1[i][j][k].md_temp * 3 * l1[i][j][k].natoms);
	  } else {
	    l1[i][j][k].xi = 0.0;
	  }
	  l2[i][j][k].xi = l1[i][j][k].xi;
	}
      }
    }
  }
}


void ttm_writeout(int number)
{
  int n,nlocal;
  int i,j,k;
  ttm_Element * lglobal;
  ttm_Element * llocal;

  n = global_fd_dim.x*global_fd_dim.y*global_fd_dim.z;
  nlocal=(local_fd_dim.x-2)*(local_fd_dim.y-2)*(local_fd_dim.z-2);
#ifdef DEBUG
  assert(nlocal == n/num_cpus);
#endif

#ifdef MPI2
  MPI_Alloc_mem ( nlocal * sizeof(ttm_Element), MPI_INFO_NULL, &llocal );
#else
  llocal=malloc(nlocal * sizeof(ttm_Element));
#endif

  for (i=1; i<local_fd_dim.x-1; ++i)
  {
    for (j=1; j<local_fd_dim.y-1; ++j)
    {
      for (k=1; k<local_fd_dim.z-1; ++k)
      { /* all the "-1" and "-2" because we don't store ghost layers
           and don't want to waste the space in llocal */
	llocal[ (i-1)*(local_fd_dim.y-2)*(local_fd_dim.z-2)
	       +(j-1)*(local_fd_dim.z-2)
	       + k-1 ] = l1[i][j][k];
      }
    }
  }

#ifdef MPI
  if (myid==0)
  {
#ifdef MPI2
    MPI_Alloc_mem ( n * sizeof(ttm_Element), MPI_INFO_NULL, &lglobal );
#else
    lglobal = malloc (n*sizeof(ttm_Element));
#endif /*MPI2*/
  }

  /* Note: This only works because mpi_element2 includes the last element
   * of ttm_Element (v_com.z). If this is changed, you need to work with
   * displacements or something (e.g. MPI standard, Ex. 4.5) */
  MPI_Gather( llocal, nlocal, mpi_element2,
      lglobal, nlocal, mpi_element2, 0, cpugrid );

#else /* no MPI */
  lglobal = llocal;
#endif /* MPI */

  if (myid==0)
  {
    FILE *outfile;
    char fname[255];
    sprintf(fname,"%s.%d.ttm",outfilename,number);
    outfile=fopen(fname, "w");
    if (NULL==outfile) error ("Cannot open ttm file for writing.\n");
    fprintf(outfile,
       	"#x y z natoms temp md_temp xi source v_com.x v_com.y v_com.z\n");
    for (i=0;i<global_fd_dim.x;++i)
    {
      for (j=0;j<global_fd_dim.y;++j)
      {
	for (k=0;k<global_fd_dim.z;++k)
	{
	  int index;

#ifdef MPI
	  ivektor from_process; /* we need to look in the data from
				   the process with these grid coords */
	  /* data from the processes is sorted by rank. */
	  from_process.x=i/(local_fd_dim.x-2);
	  from_process.y=j/(local_fd_dim.y-2);
	  from_process.z=k/(local_fd_dim.z-2);
	  index=cpu_grid_coord(from_process)*nlocal + 
	    (i%(local_fd_dim.x-2))*(local_fd_dim.y-2)*(local_fd_dim.z-2) +
	    (j%(local_fd_dim.y-2))*(local_fd_dim.z-2) +
	    (k%(local_fd_dim.z-2));
#else /* no MPI */
	  index=
	    i*(local_fd_dim.y-2)*(local_fd_dim.z-2) +
	    j*(local_fd_dim.z-2) +
	    k;
#endif /* MPI*/

	  fprintf(outfile, "%d %d %d %d %e %e %e %e %e %e %e\n",
	      i,j,k,lglobal[index].natoms,lglobal[index].temp,
	      lglobal[index].md_temp, lglobal[index].xi,
	      lglobal[index].source,
	      lglobal[index].v_com.x, lglobal[index].v_com.y,
	      lglobal[index].v_com.z);
	}
      }
    }
    fclose(outfile);
  }

#ifdef MPI
  if (myid==0)
  {
#ifdef MPI2
    MPI_Free_mem (lglobal);
#else
    free (lglobal);
#endif
  }
#endif

#ifdef MPI2
  MPI_Free_mem (llocal);
#else
  free(llocal);
#endif

}


#ifdef MPI
void ttm_create_mpi_datatypes(void)
{
  { /* type for our basic struct */

    /* we don't send unneeded elements of struct, i.e. 
     *                             md_cellptrs, xi, md_temp, v_com, source */

    /* elements to be sent:        natoms (to determine if cell is active)
     *                             temp (electron temperature).
     *                             (MPI_UB to set upper bound and skip the rest) */

    ttm_Element tmpelement;
    ttm_Element * tmpelement_pointer;
    tmpelement_pointer=&tmpelement;
    MPI_Aint tmpaddr;
    int blockcounts[3]={1,1,1};
    MPI_Datatype types[3]={MPI_INT, MPI_DOUBLE, MPI_UB};
    MPI_Aint displs[3];  

    MPI_Address(tmpelement_pointer, &tmpaddr);
    MPI_Address(&tmpelement.natoms, &displs[0]);
    MPI_Address(&tmpelement.temp, &displs[1]);
    tmpelement_pointer++;
    MPI_Address(tmpelement_pointer, &displs[2]);

    displs[2]-=tmpaddr;
    displs[1]-=tmpaddr;
    displs[0]-=tmpaddr;

    MPI_Type_struct(3,blockcounts,displs,types,&mpi_element);
    MPI_Type_commit(&mpi_element);
  }
  { /* type for our basic struct, used for ttm file output */

    /* we don't send unneeded elements of struct, i.e. 
     *                             md_cellptrs */

    /* elements to be sent:        natoms, temp,
     *                             xi, md_temp, v_com, source. */

    int i;
    ttm_Element tmpelement;
    ttm_Element * tmpelement_pointer;
    tmpelement_pointer=&tmpelement;
    MPI_Aint tmpaddr;
    int blockcounts[9]={1,1,1,1,1,1,1,1,1};
    MPI_Datatype types[9]={MPI_INT,
                           MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
                           MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_UB};
    MPI_Aint displs[9];  

    MPI_Address(&tmpelement, &tmpaddr);
    MPI_Address(&tmpelement.natoms, &displs[0]);
    MPI_Address(&tmpelement.temp, &displs[1]);
    MPI_Address(&tmpelement.xi, &displs[2]);
    MPI_Address(&tmpelement.md_temp, &displs[3]);
    MPI_Address(&tmpelement.source, &displs[4]);
    MPI_Address(&tmpelement.v_com.x, &displs[5]);
    MPI_Address(&tmpelement.v_com.y, &displs[6]);
    MPI_Address(&tmpelement.v_com.z, &displs[7]);
    tmpelement_pointer++;
    MPI_Address(tmpelement_pointer, &displs[8]);

    for (i=0; i<9; ++i)
    {
      displs[i]-=tmpaddr;
    }

    MPI_Type_struct(9,blockcounts,displs,types,&mpi_element2);
    MPI_Type_commit(&mpi_element2);
  }

  /* datatype for one string of elements along z (short of 2 lattice points) */
  MPI_Type_contiguous(local_fd_dim.z-2, mpi_element, &mpi_zrow);
  MPI_Type_commit(&mpi_zrow);

  { /* add displacements to skip ghost layers (mpi_zrow_block) */
    ttm_Element tmpel;
    ttm_Element * tmpstore;
    tmpstore=&tmpel;

    int blockcounts[3]={1,1,1};
    MPI_Datatype types[3]={MPI_LB, mpi_zrow, MPI_UB};
    MPI_Aint displs[3];

    MPI_Address(tmpstore,&displs[0]);
    tmpstore += 1;
    MPI_Address(tmpstore,&displs[1]);
    tmpstore += 1 + (local_fd_dim.z-2);
    MPI_Address(tmpstore,&displs[2]);

    displs[2]-=displs[0];	
    displs[1]-=displs[0];
    displs[0]-=displs[0];

    MPI_Type_struct(3, blockcounts, displs, types, &mpi_zrow_block);
    MPI_Type_commit(&mpi_zrow_block);
  }


  /* datatype for one layer of elements perpendicular to x
   * (short of 2 strings along z-axis)                    */
  MPI_Type_contiguous(local_fd_dim.y-2, mpi_zrow_block, &mpi_xplane);
  MPI_Type_commit(&mpi_xplane);

  { /* add displacements to skip ghost layers (mpi_xplane_block) */
    ttm_Element tmpel;
    ttm_Element * tmpstore;
    tmpstore=&tmpel;

    int blockcounts[3]={1,1,1};
    MPI_Datatype types[3]={MPI_LB, mpi_xplane, MPI_UB};
    MPI_Aint displs[3];

    MPI_Address(tmpstore,&displs[0]);
    tmpstore += 2 + (local_fd_dim.z-2);
    MPI_Address(tmpstore,&displs[1]);
    tmpstore += (2 + (local_fd_dim.z-2))*((local_fd_dim.y-2) + 1);
    MPI_Address(tmpstore,&displs[2]);

    displs[2]-=displs[0];	
    displs[1]-=displs[0];
    displs[0]-=displs[0];

    MPI_Type_struct(3, blockcounts, displs, types, &mpi_xplane_block);
    MPI_Type_commit(&mpi_xplane_block);
  }

  /* datatype for one layer of elements perpendicular to y (short of 2) */
  MPI_Type_vector((local_fd_dim.x-2), 1, 2+(local_fd_dim.y-2),
                  mpi_zrow_block, &mpi_yplane);
  MPI_Type_commit(&mpi_yplane);

  { /* add displacements to skip ghost layers (mpi_yplane_block) */
    ttm_Element tmpel;
    ttm_Element * tmpstore;
    tmpstore=&tmpel;

    int blockcounts[3]={1,1,1};
    MPI_Datatype types[3]={MPI_LB, mpi_yplane, MPI_UB};
    MPI_Aint displs[3];

    MPI_Address(tmpstore,&displs[0]);
    tmpstore += (2 + (local_fd_dim.z-2))*
                (2 + (local_fd_dim.y-2));
    MPI_Address(tmpstore,&displs[1]);
    tmpstore += (1 + (local_fd_dim.x-2))*
                (2 + (local_fd_dim.z-2))*
		(2 + (local_fd_dim.y-2));
    MPI_Address(tmpstore,&displs[2]);

    displs[2]-=displs[0];	
    displs[1]-=displs[0];
    displs[0]-=displs[0];

    MPI_Type_struct(3, blockcounts, displs, types, &mpi_yplane_block);
    MPI_Type_commit(&mpi_yplane_block);
  }


  /* datatype for one string of elements along y (short of 2) */
  MPI_Type_vector((local_fd_dim.y-2), 1, 2+(local_fd_dim.z-2),
                  mpi_element, &mpi_yrow);
  MPI_Type_commit(&mpi_yrow);  

  /* add displacements to create datatype from which mpi_zplane will be built
   * (again, skipping ghost layers) */
  { /* mpi_yrow_block */ 
    ttm_Element tmpel;
    ttm_Element * tmpstore;
    tmpstore=&tmpel;

    int blockcounts[3] = {1,1,1};
    MPI_Datatype types[3]={MPI_LB, mpi_yrow, MPI_UB};
    MPI_Aint displs[3];

    MPI_Address(tmpstore,&displs[0]);
    tmpstore += 2 + (local_fd_dim.z-2);
    MPI_Address(tmpstore,&displs[1]);
    tmpstore += (1 + (local_fd_dim.y-2)) *
                (2 + (local_fd_dim.z-2));
    MPI_Address(tmpstore,&displs[2]);

    displs[2]-=displs[0];
    displs[1]-=displs[0];
    displs[0]-=displs[0];

    MPI_Type_struct(3, blockcounts, displs, types, &mpi_yrow_block);
    MPI_Type_commit(&mpi_yrow_block);
  }

#ifdef DEBUG
  if(myid==0)
  { /* output size/extent comparisons */
    int size;
    MPI_Aint extent;

    MPI_Type_size(mpi_zrow, &size);
    MPI_Type_extent(mpi_zrow, &extent);
    printf("Size / Extent of mpi_zrow: %d / %ld\n", size, (long)extent);

    MPI_Type_size(mpi_yrow, &size);
    MPI_Type_extent(mpi_yrow, &extent);
    printf("Size / Extent of mpi_yrow: %d / %ld\n", size, (long)extent);

    MPI_Type_size(mpi_yrow_block, &size);
    MPI_Type_extent(mpi_yrow_block, &extent);
    printf("Size / Extent of mpi_yrow_block: %d / %ld\n", size, (long)extent);

    MPI_Type_size(mpi_element, &size);
    MPI_Type_extent(mpi_element, &extent);
    printf("Size / Extent of mpi_element: %d / %ld\n", size, (long)extent);

  }
#endif /*DEBUG*/

  /* datatype for one layer of elements perpendicular to z
   * (short of 2 strings) */
  MPI_Type_contiguous((local_fd_dim.x-2), mpi_yrow_block, &mpi_zplane);
  MPI_Type_commit(&mpi_zplane); 

  { /* add displacements to skip ghost layers (mpi_zplane_block) */
    ttm_Element tmpel;
    ttm_Element * tmpstore;
    tmpstore=&tmpel;

    int blockcounts[3]={1,1,1};
    MPI_Datatype types[3]={MPI_LB, mpi_zplane, MPI_UB};
    MPI_Aint displs[3];

    MPI_Address(tmpstore,&displs[0]);
    tmpstore += (2 + (local_fd_dim.z-2)) *
                (2 + (local_fd_dim.y-2));
    MPI_Address(tmpstore,&displs[1]);
    tmpstore += (1 + (local_fd_dim.x-2)) *
                (2 + (local_fd_dim.z-2)) *
		(2 + (local_fd_dim.y-2));
    MPI_Address(tmpstore,&displs[2]);

    displs[2]-=displs[0];	
    displs[1]-=displs[0];
    displs[0]-=displs[0];

    MPI_Type_struct(3, blockcounts, displs, types, &mpi_zplane_block);
    MPI_Type_commit(&mpi_zplane_block);
  }
}

void ttm_fill_ghost_layers(void)
{
  /** MPI communication (can occur before and/or during MD calculations?) */
  /* Remember:
   * east -> -x
   * west -> +x
   * north-> -y
   * south-> +y
   * up   -> -z
   * down -> +z
   * *************/
  int i,j,k;

  /* x direction */
  if(pbc_dirs.x==1 || (my_coord.x != 0 && my_coord.x != cpu_dim.x-1) )
  {
    /* send left slice to left neighbor. */
    /* Simultaneously receive slice from right neighbor */
    MPI_Sendrecv(&l1[1][0][0],1,mpi_xplane_block,nbeast,7100,
	&l1[(local_fd_dim.x-2)+1][0][0],1,mpi_xplane_block,nbwest,7100,
	cpugrid,&stati[0]);


    /* send right slice to right neighbor. */
    /* Simultaneously receive slice from left neighbor */
    MPI_Sendrecv(&l1[(local_fd_dim.x-2)][0][0],1,mpi_xplane_block,nbwest,7200,
	&l1[0][0][0],1,mpi_xplane_block,nbeast,7200,
	cpugrid,&stati[1]); 
  }
  else /* no pbc and we are at the surface */
    if (my_coord.x==0 && my_coord.x!=cpu_dim.x-1) /* left surface */
    { 


      /* only receive from right */
      MPI_Recv(&l1[(local_fd_dim.x-2)+1][0][0],1,mpi_xplane_block,nbwest,7100,
	  cpugrid,&stati[0]);

      /* only send to right */
      MPI_Send(&l1[(local_fd_dim.x-2)][0][0],1,mpi_xplane_block,nbwest,7200,
	  cpugrid);

      /* left ghost layer receives reflecting bc */
      for (j=1;j<=(local_fd_dim.y-2);j++)
      {
	for (k=1;k<=(local_fd_dim.z-2);k++)
	{
	  /* no atoms -> no conduction. */
	  l1[0][j][k].natoms=0;
	}
      }


    }
    else if (my_coord.x==cpu_dim.x-1 && my_coord.x!=0) /* right surface */
    { 

      /* only send to left */
      MPI_Send(&l1[1][0][0],1,mpi_xplane_block,nbeast,7100,
	  cpugrid);

      /* only receive from left */
      MPI_Recv(&l1[0][0][0],1,mpi_xplane_block,nbeast,7200,
	  cpugrid,&stati[1]);

      /* right ghost layer receives reflecting bc */
      for (j=1;j<=(local_fd_dim.y-2);j++)
      {
	for (k=1;k<=(local_fd_dim.z-2);k++)
	{
	  /* no atoms -> no conduction. */
	  l1[local_fd_dim.x-1][j][k].natoms=0;
	}
      }


    } else
    { /* two surfaces, just apply reflecting bc */
      for (j=1;j<=(local_fd_dim.y-2);j++)
      {
	for (k=1;k<=(local_fd_dim.z-2);k++)
	{
	  l1[0][j][k].natoms=l1[local_fd_dim.x-1][j][k].natoms=0;
	}
      }

    }



    /* y direction */
    if(pbc_dirs.y==1 || (my_coord.y != 0 && my_coord.y != cpu_dim.y-1) )
    {
      /* send left slice to left neighbor. */
      /* Simultaneously receive slice from right neighbor */
      MPI_Sendrecv(&l1[0][1][0],1,mpi_yplane_block,nbnorth,710,
	  &l1[0][(local_fd_dim.y-2)+1][0],1,mpi_yplane_block,nbsouth,710,
	  cpugrid,&stati[2]);
      /* send right slice to right neighbor. */
      /* Simultaneously receive slice from left neighbor */
      MPI_Sendrecv(&l1[0][(local_fd_dim.y-2)][0],1,mpi_yplane_block,nbsouth,720,
	  &l1[0][0][0],1,mpi_yplane_block,nbnorth,720,
	  cpugrid,&stati[3]);
    }
    else /* no pbc and we are at the surface */
      if (my_coord.y==0) /* left surface */
      { 
	/* only receive from right */
	MPI_Recv(&l1[0][(local_fd_dim.y-2)+1][0],1,mpi_yplane_block,nbsouth,710,
	    cpugrid,&stati[2]);
	/* only send to right */
	MPI_Send(&l1[0][(local_fd_dim.y-2)][0],1,mpi_yplane_block,nbsouth,720,
	    cpugrid);

	/* left ghost layer receives reflecting bc */
	for (i=1;i<=(local_fd_dim.x-2);i++)
	{
	  for (k=1;k<=(local_fd_dim.z-2);k++)
	  {
	    /* no atoms -> no conduction. */
	    l1[i][0][k].natoms=0;
	  }
	}
      }
      else if (my_coord.y==cpu_dim.y-1) /* right surface */
      { 
	/* only send to left */
	MPI_Send(&l1[0][1][0],1,mpi_yplane_block,nbnorth,710,
	    cpugrid);
	/* only receive from left */
	MPI_Recv(&l1[0][0][0],1,mpi_yplane_block,nbnorth,720,
	    cpugrid,&stati[3]);

	/* right ghost layer receives reflecting bc */
	for (i=1;i<=(local_fd_dim.x-2);i++)
	{
	  for (k=1;k<=(local_fd_dim.z-2);k++)
	  {
	    /* no atoms -> no conduction. */
	    l1[i][local_fd_dim.y-1][k].natoms=0;
	  }
	}
      } else { error("This should be logically impossible.\n");}


    /* z direction */
    if(pbc_dirs.z==1 || (my_coord.z != 0 && my_coord.z != cpu_dim.z-1) )
    {
      /* send left slice to left neighbor. */
      /* Simultaneously receive slice from right neighbor */
      MPI_Sendrecv(&l1[0][0][1],1,mpi_zplane_block,nbup,71,
	  &l1[0][0][(local_fd_dim.z-2)+1],1,mpi_zplane_block,nbdown,71,
	  cpugrid,&stati[4]);
      /* send right slice to right neighbor. */
      /* Simultaneously receive slice from left neighbor */
      MPI_Sendrecv(&l1[0][0][(local_fd_dim.z-2)],1,mpi_zplane_block,nbdown,72,
	  &l1[0][0][0],1,mpi_zplane_block,nbup,72,
	  cpugrid,&stati[5]);
    }
    else /* no pbc and we are at the surface */
      if (my_coord.z==0) /* left surface */
      { 
	/* only receive from right */
	MPI_Recv(&l1[0][0][(local_fd_dim.z-2)+1],1,mpi_zplane_block,nbdown,71,
	    cpugrid,&stati[4]);
	/* only send to right */
	MPI_Send(&l1[0][0][(local_fd_dim.x-2)],1,mpi_zplane_block,nbdown,72,
	    cpugrid);

	/* left ghost layer receives reflecting bc */
	for (i=1;i<=(local_fd_dim.x-2);i++)
	{
	  for (j=1;j<=(local_fd_dim.y-2);j++)
	  {
	    /* no atoms -> no conduction. */
	    l1[i][j][0].natoms=0;
	  }
	}
      }
      else if (my_coord.z==cpu_dim.z-1) /* right surface */
      { 
	/* only send to left */
	MPI_Send(&l1[0][0][1],1,mpi_zplane_block,nbup,71,
	    cpugrid);
	/* only receive from left */
	MPI_Recv(&l1[0][0][0],1,mpi_zplane_block,nbup,72,
	    cpugrid,&stati[5]);

	/* right ghost layer receives reflecting bc */
	for (i=1;i<=(local_fd_dim.x-2);i++)
	{
	  for (j=1;j<=(local_fd_dim.y-2);j++)
	  {
	    /* no atoms -> no conduction. */
	    l1[i][j][local_fd_dim.z-1].natoms=0;
	  }
	}
      } else { error("This should be logically impossible.\n");}

}

#else

/* Serial version */
void ttm_fill_ghost_layers(void)
{
  /** Serial computation.
   * Copy into ghost layers on opposite sides of the sample
   ***/
  int i,j,k;

  /* x direction */
  for (j=1;j<=(local_fd_dim.y-2);j++)
  {
    for (k=1;k<=(local_fd_dim.z-2);k++)
    {
      if (pbc_dirs.x==1)
      {
	l1[0][j][k]=l1[local_fd_dim.x-2][j][k];
	l1[local_fd_dim.x-1][j][k]=l1[1][j][k];
      } else
      {
	/* no atoms -> no conduction. */
	l1[0][j][k].natoms=l1[local_fd_dim.x-1][j][k].natoms=0;
      }
    }
  }
  /* y direction */
  for (i=1;i<=(local_fd_dim.x-2);i++)
  {
    for (k=1;k<=(local_fd_dim.z-2);k++)
    {	
      if (pbc_dirs.y==1)
      {
	l1[i][0][k]=l1[i][local_fd_dim.y-2][k];
	l1[i][local_fd_dim.y-1][k]=l1[i][1][k];
      } else
      {
	/* no atoms -> no conduction. */
	l1[i][0][k].natoms=l1[i][local_fd_dim.y-1][k].natoms=0;
      }

    }
  }
  /* z direction */
  for (i=1;i<=(local_fd_dim.x-2);i++)
  {
    for (j=1;j<=(local_fd_dim.y-2);j++)
    {
      if (pbc_dirs.z==1)
      {
	l1[i][j][0]=l1[i][j][local_fd_dim.z-2];
	l1[i][j][local_fd_dim.z-1]=l1[i][j][1];
      } else
      {
	/* no atoms -> no conduction. */
	l1[i][j][0].natoms=l1[i][j][local_fd_dim.z-1].natoms=0;
      }
    }
  }


}
#endif /*MPI*/

