
/******************************************************************************
*
* imd_main_mpi_3d.c -- main loop, mpi specific part, three dimensions
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"


/******************************************************************************
*
* calc_forces 
*
* The forces of the atoms are calulated here. To achive this, atoms on
* the surface of a cpu are exchanged with the neigbours.
*
* If AR is defined, we use actio=reactio even across CPUs, otherwise we don't
*
* The force calculation is split into those steps:
*
* i)   send atoms positions of cells on surface neighbours, 
*      receive atom positions from neigbours
* ii)  zero forces on all cells (local and buffer)
* iii) calculate forces in local cells, use lower half of neigbours 
*      for each cell and use actio==reactio
* iv)  calculate forces also for upper half of neighbours for all cells
*      that are on the upper surface
* iv)  or send forces back and add them
*
******************************************************************************/

#ifndef MONOLJ

void calc_forces(void)
{
  int n, k;
  real tmp, tmpvir;
  vektor tmpvec;

  /* clear global accumulation variables */
  tot_pot_energy = 0.0;
  virial         = 0.0;
  vir_vect.x     = 0.0;
  vir_vect.y     = 0.0;
  vir_vect.z     = 0.0;

#ifdef EAM
#ifdef AR
  error("EAM force routine not defined in case of actio = reactio.");
#endif
  memset(eam_rho,   0, (natoms+1)*        sizeof(real));
  memset(eam_ij,    0, (natoms+1)*eam_len*sizeof(real));
  memset(eam_dij_x, 0, (natoms+1)*eam_len*sizeof(real));
  memset(eam_dij_y, 0, (natoms+1)*eam_len*sizeof(real));
  memset(eam_dij_z, 0, (natoms+1)*eam_len*sizeof(real));
#endif /* EAM */

  /* clear per atom accumulation variables */
  for (k=0; k<nallcells; ++k) {
    int  i,j;
    cell *p;
    p = cell_array + k;
    for (i = 0;i < p->n; ++i) {
      p->kraft X(i) = 0.0;
      p->kraft Y(i) = 0.0;
      p->kraft Z(i) = 0.0;
#ifdef UNIAX
      p->dreh_moment X(i) = 0.0;
      p->dreh_moment Y(i) = 0.0;
      p->dreh_moment Z(i) = 0.0;
#endif
#ifdef TRANSPORT
      p->heatcond[i] = 0.0;
#endif      
#ifdef STRESS_TENS
      p->presstens X(i) = 0.0;
      p->presstens Y(i) = 0.0;
      p->presstens Z(i) = 0.0;
      p->presstens_offdia X(i) = 0.0;
      p->presstens_offdia Y(i) = 0.0;
      p->presstens_offdia Z(i) = 0.0;
#endif      
#ifndef MONOLJ
      p->pot_eng[i] = 0.0;
#endif
#ifdef COVALENT
      for (j=0; j<p->n; ++j) {
        p->neigh[j]->n = 0;
      }
#endif
#ifdef EAM2
      p->eam2_rho_h[i] = 0.0; /* zero host electron density at atom site */
#endif
    }
  }

  /* What follows is the standard one-cpu force 
     loop acting on our local data cells */

  /* compute forces for all pairs of cells */
  for (n=0; n<6; ++n) {
    for (k=0; k<npairs[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n] + k;
      pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
      pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
      pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
#ifdef EAM
      /* first EAM call */
      do_forces_eam_1(cell_array + P->np, cell_array + P->nq, pbc);
#elif defined(EAM2)
      eam2_do_forces1(cell_array + P->np, cell_array + P->nq, pbc);
#elif defined(UNIAX)
      do_forces_uniax(cell_array + P->np, cell_array + P->nq, pbc);
#else
      do_forces(cell_array + P->np, cell_array + P->nq, pbc);
#endif
    }
  }

#ifdef EAM2
  /* if EAM2, we have to loop a second time over pairs of distinct cells */
  for (n=0; n<6; ++n) {
    for (k=0; k<npairs[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n]+k;
      pbc.x = -(P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x);
      pbc.y = -(P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y);
      pbc.z = -(P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z);
      if (P->np != P->nq)
        eam2_do_forces1(cell_array + P->nq, cell_array + P->np, pbc);
    }
  }
#endif

#ifdef COVALENT
  /* complete neighbor tables for remaining pairs of cells */
  for (n=0; n<6; ++n) {
    for (k=npairs[n]; k<npairs2[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n] + k;
      pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
      pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
      pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
      do_neightab(cell_array + P->np, cell_array + P->nq, pbc);
    }
  }
#endif

  mpi_addtime(&time_calc_local);
  
#ifndef AR  

  /* If we don't use actio=reactio accross the cpus, we have do do
     the force loop also on the other half of the neighbours for the 
     cells on the surface of the CPU */

  /* potential energy and virial are already complete; to avoid double
     counting, we keep a copy of the current value, which we use later */

  tmp      = tot_pot_energy;
  tmpvir   = virial;
  tmpvec.x = vir_vect.x;
  tmpvec.y = vir_vect.y;
  tmpvec.z = vir_vect.z;

  /* compute forces for remaining pairs of cells */
  for (n=0; n<6; ++n) {
    for (k=npairs[n]; k<npairs2[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n] + k;
      pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
      pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
      pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
#ifdef EAM
      /* first EAM call */
      do_forces_eam_1(cell_array + P->np, cell_array + P->nq, pbc);
#elif defined(EAM2)
      eam2_do_forces1(cell_array + P->np, cell_array + P->nq, pbc);
#elif defined(UNIAX)
      do_forces_uniax(cell_array + P->np, cell_array + P->nq, pbc);
#else
      do_forces(cell_array + P->np, cell_array + P->nq, pbc);
#endif
    }
  }

  /* use the previously saved values of potential energy and virial */
  tot_pot_energy = tmp;
  virial     = tmpvir;
  vir_vect.x = tmpvec.x;
  vir_vect.y = tmpvec.y;
  vir_vect.z = tmpvec.z;

  mpi_addtime(&time_calc_nonlocal);

#endif  /* ... ifndef AR */

#ifdef EAM
  /* EAM cohesive function potential */
  for (k=0; k<ncells; ++k) {
    do_forces_eam_2(cell_array + CELLS(k));
  }
#endif /* EAM */

#ifdef TTBP
  /* TTBP: three body potential */
  for (k=0; k<ncells; ++k) {
    do_forces_ttbp(cell_array + CELLS(k));
  }
#endif /* TTBP */

#ifdef TERSOFF
  /* Tersoff potential */
  for (k=0; k<ncells; ++k) {
    do_forces_tersoff(cell_array + CELLS(k));
  }
#endif /* TERSOFF */

#ifdef EAM2

  send_eam2_rho_h();

  /* second EAM2 loop over all cells pairs */
  for (n=0; n<6; ++n) {
    for (k=0; k<npairs[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n]+k;
      pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
      pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
      pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
      eam2_do_forces2(cell_array + P->np, cell_array + P->nq, pbc);
    }
  }

  /* if EAM2, we have to loop a second time over pairs of distinct cells */
  for (n=0; n<6; ++n) {
    for (k=0; k<npairs[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n]+k;
      pbc.x = -(P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x);
      pbc.y = -(P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y);
      pbc.z = -(P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z);
      if (P->np != P->nq)
        eam2_do_forces2(cell_array + P->nq, cell_array + P->np, pbc);
    }
  }

#ifndef AR

  /* If we don't use actio=reactio accross the cpus, we have do do
     the force loop also on the other half of the neighbours for the 
     cells on the surface of the CPU */

  /* potential energy and virial are already complete; to avoid double
     counting, we keep a copy of the current value, which we use later */

  tmp      = tot_pot_energy;
  tmpvir   = virial;
  tmpvec.x = vir_vect.x;
  tmpvec.y = vir_vect.y;
  tmpvec.z = vir_vect.z;

  /* compute forces for remaining pairs of cells */
  for (n=0; n<6; ++n) {
    for (k=npairs[n]; k<npairs2[n]; ++k) {
      vektor pbc;
      pair *P;
      P = pairs[n]+k;
      pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
      pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
      pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
      eam2_do_forces2(cell_array + P->np, cell_array + P->nq, pbc);
    }
  }

  /* use the previously saved values of potential energy and virial */
  tot_pot_energy = tmp;
  virial     = tmpvir;
  vir_vect.x = tmpvec.x;
  vir_vect.y = tmpvec.y;
  vir_vect.z = tmpvec.z;

#endif /* not AR */

#endif /* EAM2 */

  /* sum up results of different CPUs */
  MPI_Allreduce( &tot_pot_energy, &tmp, 1, MPI_REAL, MPI_SUM, cpugrid); 
  tot_pot_energy = tmp; 

  MPI_Allreduce( &virial,   &tmpvir,    1, MPI_REAL, MPI_SUM, cpugrid);
  virial = tmpvir;

  MPI_Allreduce( &vir_vect, &tmpvec,  DIM, MPI_REAL, MPI_SUM, cpugrid);
  vir_vect.x = tmpvec.x;
  vir_vect.y = tmpvec.y;
  vir_vect.z = tmpvec.z;

}

#else /* MONOLJ */

/******************************************************************************
*
*  this version is kept for world records only - it uses less memory
*
******************************************************************************/

void calc_forces(void)
{
  cell *p,*q;
  int i,j,k;
  int l,m,n;
  int r,s,t;
  int u,v,w;
  vektor pbc = {0.0,0.0,0.0};
  real tmp, tmpvir;
#ifdef P_AXIAL
  vektor tmpvec;
#endif

  tot_pot_energy = 0.0;
  virial         = 0.0;
#ifdef P_AXIAL
  vir_vect.x     = 0.0;
  vir_vect.y     = 0.0;
  vir_vect.z     = 0.0;
#endif

  /* Zero Forces */
  for (k=0; k<nallcells; ++k) {
    int  i;
    cell *p;
    p = cell_array + k;
    for (i=0; i<p->n; ++i) {
      p->kraft X(i) = 0.0;
      p->kraft Y(i) = 0.0;
      p->kraft Z(i) = 0.0;
    }
  }

  /* What follows is the standard one-cpu force 
     loop acting on our local data cells */

  /* for each cell in bulk */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.y-1; ++j)
      for (k=1; k < cell_dim.z-1; ++k) {

        p = PTR_3D_V(cell_array,i,j,k,cell_dim);

	/* For half of the neighbours of this cell */
	for (l=0; l <= 1; ++l)
	  for (m=-l; m <= 1; ++m)
	    for (n=(l==0 ? -m  : -l ); n <= 1; ++n) {

	      /* Calculate Indicies of Neighbour */
	      r = i + l;
	      s = j + m;
	      t = k + n; 

	      /* Neighbour (note that p==q ist possible) */
	      q = PTR_3D_V(cell_array,r,s,t,cell_dim);

	      /* Apply periodic boundaries */
	      pbc = global_pbc(r,s,t);
	      do_forces(p,q,pbc);
	    }
      }

  mpi_addtime(&time_calc_local);
  
#ifndef AR  
  /* Calculate forces on boundary half of cell */
  /* potential energy and virial are already complete; to avoid double
     counting, we keep a copy of the current value, which we use later */

  tmp      = tot_pot_energy;
  tmpvir   = virial;
#ifdef P_AXIAL
  tmpvec.x = vir_vect.x;
  tmpvec.y = vir_vect.y;
  tmpvec.z = vir_vect.z;
#endif

  /* for each cell in bulk */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.y-1; ++j)
      for (k=1; k < cell_dim.z-1; ++k) {

        p = PTR_3D_V(cell_array,i,j,k,cell_dim);

	/* For half of the neighbours of this cell */
	for (l=0; l <= 1; ++l)
	  for (m=-l; m <= 1; ++m)
	    for (n=(l==0 ? -m  : -l ); n <= 1; ++n) {

	      /* Calculate Indicies of Neighbour */
	      r = i - l;
	      s = j - m;
	      t = k - n;

              /* if second cell is a buffer cell */
              if ((r == 0) || (r == cell_dim.x-1) || 
                  (s == 0) || (s == cell_dim.y-1) ||
                  (t == 0) || (t == cell_dim.z-1)) 
              {
		q = PTR_3D_V(cell_array,r,s,t,cell_dim);
		/* Apply periodic boundaries */
		pbc = global_pbc(r,s,t);
		do_forces(p,q,pbc);
	      }
	    }
      }

  /* use the previously saved values of potential energy and virial */
  tot_pot_energy = tmp;
  virial     = tmpvir;
#ifdef P_AXIAL
  vir_vect.x = tmpvec.x;
  vir_vect.y = tmpvec.y;
  vir_vect.z = tmpvec.z;
#endif

  mpi_addtime(&time_calc_nonlocal);
#endif  /* ... ifndef AR */

  /* sum up results of different CPUs */
  MPI_Allreduce( &tot_pot_energy, &tmp, 1, MPI_REAL, MPI_SUM, cpugrid); 
  tot_pot_energy = tmp; 

  MPI_Allreduce( &virial, &tmp,         1, MPI_REAL, MPI_SUM, cpugrid);
  virial = tmp;

#ifdef P_AXIAL
  MPI_Allreduce( &vir_vect, &tmpvec,  DIM, MPI_REAL, MPI_SUM, cpugrid);
  vir_vect.x = tmpvec.x;
  vir_vect.y = tmpvec.y;
  vir_vect.z = tmpvec.z;
#endif

}


/******************************************************************************
*
*  global_pbc tells if a local buffer cell is across the boundaries
*
*  only used for MONOLJ
*
******************************************************************************/

vektor global_pbc(int i, int j, int k)    
{
  ivektor global_coord;
  ivektor local_coord;
  vektor pbc = { 0.0, 0.0, 0.0};

  local_coord.x = i;
  local_coord.y = j;
  local_coord.z = k;

  global_coord.x = local_coord.x - 1 + my_coord.x * (cell_dim.x - 2);
  global_coord.y = local_coord.y - 1 + my_coord.y * (cell_dim.y - 2);
  global_coord.z = local_coord.z - 1 + my_coord.z * (cell_dim.z - 2);

  if (global_coord.x < 0) {
    pbc.x -= box_x.x;      
    pbc.y -= box_x.y;
    pbc.z -= box_x.z;
  }

  if (global_coord.x >= global_cell_dim.x) {
    pbc.x += box_x.x;      
    pbc.y += box_x.y;
    pbc.z += box_x.z;
  }

  if (global_coord.y < 0) {
    pbc.x -= box_y.x;      
    pbc.y -= box_y.y;
    pbc.z -= box_y.z;
  }

  if (global_coord.y >= global_cell_dim.y) {
    pbc.x += box_y.x;      
    pbc.y += box_y.y;
    pbc.z += box_y.z;
  }

  if (global_coord.z < 0) {
    pbc.x -= box_z.x;      
    pbc.y -= box_z.y;
    pbc.z -= box_z.z;
  }

  if (global_coord.z >= global_cell_dim.z) {
    pbc.x += box_z.x;      
    pbc.y += box_z.y;
    pbc.z += box_z.z;
  }

  return pbc;

}

#endif /* MONOLJ */


/******************************************************************************
*
* set up mpi
*
******************************************************************************/

void init_mpi(int argc,char *argv[])
{
  ivektor nbcoord;

  /* Initialize MPI */
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&num_cpus);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  if (0 == myid) { 
    printf("%s\n", argv[0]);
    printf("Starting up MPI on %d nodes.\n",num_cpus);
  };

  /* Setup send/receive buffers */
  send_buf_north.n_max = 0;
  send_buf_south.n_max = 0;
  send_buf_east.n_max  = 0;
  send_buf_west.n_max  = 0;
  send_buf_up.n_max    = 0;
  send_buf_down.n_max  = 0;
  recv_buf_north.n_max = 0;
  recv_buf_south.n_max = 0;
  recv_buf_east.n_max  = 0;
  recv_buf_west.n_max  = 0;
  recv_buf_up.n_max    = 0;
  recv_buf_down.n_max  = 0;

}


/******************************************************************************
*
* shut down mpi
*
******************************************************************************/

void shutdown_mpi(void)
{
  MPI_Barrier(MPI_COMM_WORLD);  /* Wait for all processes to arrive */
  MPI_Finalize();               /* Shutdown */
}


/******************************************************************************
*
* send_atoms
*
* This sends the atom positions to the neighbouring cpus using
* Steve Plimptons comm scheme
*
* Atoms are sent only over the faces of the processors box in 
* the order east-west, north-south, up-down
*
* This is only used for atom redistribuion (mode==ATOMS)
* where all of the data are sent.
*
* The drawback of the Plimpton scheme are the huge comm-buffers,
* that use lots of memory. It is fast, however!!
*
******************************************************************************/

void send_atoms(int mode)
{
  MPI_Status  stateast[2],  statwest[2];
  MPI_Status statnorth[2], statsouth[2];
  MPI_Status    statup[2],  statdown[2];

  MPI_Request  reqeast[2],   reqwest[2];
  MPI_Request reqnorth[2],  reqsouth[2];
  MPI_Request    requp[2],   reqdown[2];

  /* send east */
  irecv_buf( &recv_buf_west, nbwest, &reqwest[1] );
  isend_buf( &send_buf_east, nbeast, &reqwest[0] );

  /* send west */
  irecv_buf( &recv_buf_east, nbeast, &reqeast[1] );
  isend_buf( &send_buf_west, nbwest, &reqeast[0] );

  /* Wait for atoms from west, move them to cells */
  MPI_Waitall(2, reqwest, statwest);
  MPI_Get_count( &statwest[1], MPI_REAL, &recv_buf_west.n );
  process_buffer( &recv_buf_west, mode);

  /* Wait for atoms from east, move them to cells */
  MPI_Waitall(2, reqeast, stateast);
  MPI_Get_count( &stateast[1], MPI_REAL, &recv_buf_east.n );
  process_buffer( &recv_buf_east, mode);
  
  /* Append atoms from west into north send buffer */
  copy_atoms_buf( &send_buf_north, &recv_buf_west );
  /* Append atoms from east into north send buffer */
  copy_atoms_buf( &send_buf_north, &recv_buf_east );
  /* check special case cpu_dim.y==2 */ 
  if (nbsouth!=nbnorth) {
    /* append atoms from east & west to south send buffer */
    copy_atoms_buf( &send_buf_south, &recv_buf_east );
    copy_atoms_buf( &send_buf_south, &recv_buf_west );
  }

  /* Send atoms north */
  irecv_buf( &recv_buf_south, nbsouth, &reqsouth[1] );
  isend_buf( &send_buf_north, nbnorth, &reqsouth[0] );
  
  /* Send atoms south */
  irecv_buf( &recv_buf_north, nbnorth, &reqnorth[1] );
  isend_buf( &send_buf_south, nbsouth, &reqnorth[0] );

  /* Wait for atoms from south, move them to cells */
  MPI_Waitall(2, reqsouth, statsouth);
  MPI_Get_count( &statsouth[1], MPI_REAL, &recv_buf_south.n );
  process_buffer( &recv_buf_south, mode);

  /* Wait for atoms from north, move them to cells */
  MPI_Waitall(2, reqnorth, statnorth);
  MPI_Get_count( &statnorth[1], MPI_REAL, &recv_buf_north.n );
  process_buffer( &recv_buf_north, mode);

  /* Append atoms from north to up send buffer */
  copy_atoms_buf( &send_buf_up, &recv_buf_north );
  /* Append atoms from south to up send buffer */
  copy_atoms_buf( &send_buf_up, &recv_buf_south );
  /* Append atoms from east to up send buffer */
  copy_atoms_buf( &send_buf_up, &recv_buf_east );
  /* Append atoms from west to up send buffer */
  copy_atoms_buf( &send_buf_up, &recv_buf_west );
  /* check special case cpu_dim.z==2 */ 
  if (nbdown!=nbup) {
    /* append atoms from north,south,east,west to down send buffer */
    copy_atoms_buf( &send_buf_down, &recv_buf_north );
    copy_atoms_buf( &send_buf_down, &recv_buf_south );
    copy_atoms_buf( &send_buf_down, &recv_buf_east  );
    copy_atoms_buf( &send_buf_down, &recv_buf_west  );
  }
  
  /* Send atoms up */
  irecv_buf( &recv_buf_down , nbdown, &reqdown[1]);
  isend_buf( &send_buf_up   , nbup  , &reqdown[0]);
  
  /* Send atoms down */
  irecv_buf( &recv_buf_up  , nbup  , &requp[1] );
  isend_buf( &send_buf_down, nbdown, &requp[0] );

  /* Wait for atoms from down, move them to cells */
  MPI_Waitall(2, reqdown, statdown);
  MPI_Get_count( &statdown[1], MPI_REAL, &recv_buf_down.n );
  process_buffer( &recv_buf_down, mode);   

  /* Wait for atoms from up, move them to cells */
  MPI_Waitall(2, requp, statup);
  MPI_Get_count( &statup[1], MPI_REAL, &recv_buf_up.n );
  process_buffer( &recv_buf_up, mode);  

}


/******************************************************************************
*
*  copy_cell_force  -  copy contents of one cell to another cell
*
******************************************************************************/

void copy_cell_force( int k, int l, int m, int r, int s, int t)
{
  int i, tmp_n;
  cell *from, *to;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  to   = PTR_3D_V(cell_array, r, s, t, cell_dim);

  tmp_n = from->n;
  if (tmp_n >= to->n_max) {
    to->n = 0;
    alloc_cell(to, tmp_n);
  }
  
  to->n = tmp_n;
  for (i=0; i<to->n; ++i) {
    to->ort X(i) = from->ort X(i);
    to->ort Y(i) = from->ort Y(i);
    to->ort Z(i) = from->ort Z(i);
#ifndef MONOLJ
    to->sorte[i] = from->sorte[i];
#ifdef UNIAX
    to->achse X(i) = from->achse X(i);
    to->achse Y(i) = from->achse Y(i);
    to->achse Z(i) = from->achse Z(i);
    to->shape X(i) = from->shape X(i);
    to->shape Y(i) = from->shape Y(i);
    to->shape Z(i) = from->shape Z(i);
    to->pot_well X(i) = from->pot_well X(i);
    to->pot_well Y(i) = from->pot_well Y(i);
    to->pot_well Z(i) = from->pot_well Z(i);
#endif
#endif 
#ifdef EAM
    to->nummer[i] = from->nummer[i];
#endif
  }
}


/******************************************************************************
*
*  move_atoms_force  -  move atoms from MPI buffer to cells
*
******************************************************************************/

void move_atoms_force( msgbuf *b, int k, int l, int m )
{
  int i;
  int tmp_n;
  cell *to;

  to = PTR_3D_V(cell_array, k, l, m, cell_dim);

  tmp_n = (int) b->data[ b->n++ ];
  
  if (tmp_n >= to->n_max) {
    to->n = 0;
    alloc_cell(to, tmp_n);
  }
  
  to->n = tmp_n;
  for (i=0; i<to->n; ++i) {
    to->ort X(i) = b->data[ b->n++ ];
    to->ort Y(i) = b->data[ b->n++ ];
    to->ort Z(i) = b->data[ b->n++ ];
#ifndef MONOLJ
    to->sorte[i] = (shortint) b->data[ b->n++ ];
#ifdef UNIAX
    to->achse X(i) = b->data[ b->n++ ];
    to->achse Y(i) = b->data[ b->n++ ];
    to->achse Z(i) = b->data[ b->n++ ];
    to->shape X(i) = b->data[ b->n++ ];
    to->shape Y(i) = b->data[ b->n++ ];
    to->shape Z(i) = b->data[ b->n++ ];
    to->pot_well X(i) = b->data[ b->n++ ];
    to->pot_well Y(i) = b->data[ b->n++ ];
    to->pot_well Z(i) = b->data[ b->n++ ];
#endif
#endif 
#ifdef EAM
    to->nummer[i] = (integer) b->data[ b->n++ ];
#endif
  }
  if (b->n_max <= b->n) error("Buffer overflow in move_atoms_force");
}


/******************************************************************************
*
*  copy_atoms_force  -  copy atoms to MPI buffer
*
******************************************************************************/

void copy_atoms_force( msgbuf *b, int k, int l, int m)
{
  int i;
  cell *from;
    
  from = PTR_3D_V(cell_array, k, l, m, cell_dim);

  b->data[ b->n++ ] = (real) from->n;
    
  for (i=0; i<from->n; ++i) {
    b->data[ b->n++ ] = from->ort X(i);
    b->data[ b->n++ ] = from->ort Y(i);
    b->data[ b->n++ ] = from->ort Z(i);
#ifndef MONOLJ
    b->data[ b->n++ ] = (real) from->sorte[i];
#ifdef UNIAX
    b->data[ b->n++ ] = from->achse X(i);
    b->data[ b->n++ ] = from->achse Y(i);
    b->data[ b->n++ ] = from->achse Z(i);
    b->data[ b->n++ ] = from->shape X(i);
    b->data[ b->n++ ] = from->shape Y(i);
    b->data[ b->n++ ] = from->shape Z(i);
    b->data[ b->n++ ] = from->pot_well X(i);
    b->data[ b->n++ ] = from->pot_well Y(i);
    b->data[ b->n++ ] = from->pot_well Z(i);
#endif
#endif
#ifdef EAM
    b->data[ b->n++ ] = (real) from->nummer[i];
#endif
  }
  if (b->n_max <= b->n)  error("Buffer overflow in copy_atoms_force");
}


/******************************************************************************
*
* send_atoms_force
*
* This sends the atom positions to the neighbouring cpus using
* Steve Plimptons comm scheme
*
* Atoms are sent only over the faces of the processors box in 
* the order east-west, north-south, up-down
*
* For the force calculation (mode==FORCE), only atom
* type and position are sent.
*
* The drawback of the Plimpton scheme are the huge comm-buffers,
* that use lots of memory. It is fast, however!!
*
******************************************************************************/

void send_atoms_force()
{
  int i,j;

  MPI_Status  stateast[2],  statwest[2];
  MPI_Status statnorth[2], statsouth[2];
  MPI_Status    statup[2],  statdown[2];

  MPI_Request  reqeast[2],   reqwest[2];
  MPI_Request reqnorth[2],  reqsouth[2];
  MPI_Request    requp[2],   reqdown[2];

  empty_mpi_buffers();
  empty_buffer_cells();

  /* exchange east/west */
  if (cpu_dim.x==1) {
    /* simply copy east/west atoms to buffer cells */
    for (i=1; i < cell_dim.y-1; ++i)
      for (j=1; j < cell_dim.z-1; ++j) { 
        copy_cell_force( 1, i, j, cell_dim.x-1, i, j );
        copy_cell_force( cell_dim.x-2, i, j, 0, i, j );
      }
  } else {
    /* copy east atoms into send buffer, send east */
    for (i=1; i < cell_dim.y-1; ++i)
      for (j=1; j < cell_dim.z-1; ++j) 
	copy_atoms_force( &send_buf_east, 1, i, j);
    irecv_buf( &recv_buf_west, nbwest, &reqwest[1] );
    isend_buf( &send_buf_east, nbeast, &reqwest[0] );

    /* copy west atoms into send buffer, send west */
    for (i=1; i < cell_dim.y-1; ++i) 
      for (j=1; j < cell_dim.z-1; ++j) 
	copy_atoms_force( &send_buf_west, cell_dim.x-2, i, j);
    irecv_buf( &recv_buf_east, nbeast, &reqeast[1] );
    isend_buf( &send_buf_west, nbwest, &reqeast[0] );

    /* wait for atoms from west, move them to buffer cells */
    MPI_Waitall(2, reqwest, statwest);
    recv_buf_west.n = 0;
    for (i=1; i < cell_dim.y-1; ++i)
      for (j=1; j < cell_dim.z-1; ++j) 
	move_atoms_force( &recv_buf_west, cell_dim.x-1, i, j);

    /* wait for atoms from east, move them to buffer cells */
    MPI_Waitall(2, reqeast, stateast);
    recv_buf_east.n = 0;
    for (i=1; i < cell_dim.y-1; ++i) 
      for (j=1; j < cell_dim.z-1; ++j) 
	move_atoms_force( &recv_buf_east, 0, i, j);
  }
  
  /* exchange north/south */
  if (cpu_dim.y==1) {
    /* simply copy north/south atoms to buffer cells */
    for (i=0; i < cell_dim.x; ++i)
      for (j=1; j < cell_dim.z-1; ++j) {
        copy_cell_force( i, 1, j, i, cell_dim.y-1, j );
        copy_cell_force( i, cell_dim.y-2, j, i, 0, j );
      }
  } else {
    /* copy north atoms into send buffer, send north */
    for (i=0; i < cell_dim.x; ++i) 
      for (j=1; j < cell_dim.z-1; ++j)
	copy_atoms_force( &send_buf_north, i, 1, j);
    irecv_buf( &recv_buf_south, nbsouth, &reqsouth[1] );
    isend_buf( &send_buf_north, nbnorth, &reqsouth[0] );

    /* copy south atoms into send buffer, send south */
    for (i=0; i < cell_dim.x; ++i) 
      for (j=1; j < cell_dim.z-1; ++j)
        copy_atoms_force( &send_buf_south, i, cell_dim.y-2, j);
    irecv_buf( &recv_buf_north, nbnorth, &reqnorth[1] );
    isend_buf( &send_buf_south, nbsouth, &reqnorth[0] );

    /* wait for atoms from south, move them to buffer cells */
    MPI_Waitall(2, reqsouth, statsouth);
    recv_buf_south.n = 0;
    for (i=0; i < cell_dim.x; ++i) 
      for (j=1; j < cell_dim.z-1; ++j)
	move_atoms_force( &recv_buf_south, i, cell_dim.y-1, j);

    /* wait for atoms from north, move them to buffer cells */
    MPI_Waitall(2, reqnorth, statnorth);
    recv_buf_north.n = 0;
    for (i=0; i < cell_dim.x; ++i) 
      for (j=1; j < cell_dim.z-1; ++j)
        move_atoms_force( &recv_buf_north, i, 0, j);
  }

  /* exchange up/down */
  if (cpu_dim.z==1) {
    /* simply copy up/down atoms to buffer cells */
    for (i=0; i < cell_dim.x; ++i) 
      for (j=0; j < cell_dim.y; ++j) {
	copy_cell_force( i, j, 1, i, j, cell_dim.z-1 );
        copy_cell_force( i, j, cell_dim.z-2, i, j, 0 );
      }
  } else {
    /* copy up atoms into send buffer, send up */
    for (i=0; i < cell_dim.x; ++i) 
      for (j=0; j < cell_dim.y; ++j)
	copy_atoms_force( &send_buf_up, i, j, 1);
    irecv_buf( &recv_buf_down , nbdown, &reqdown[1]);
    isend_buf( &send_buf_up   , nbup  , &reqdown[0]);

    /* copy down atoms into send buffer, send down */
    for (i=0; i < cell_dim.x; ++i) 
      for (j=0; j < cell_dim.y; ++j)
        copy_atoms_force( &send_buf_down, i, j, cell_dim.z-2);
    irecv_buf( &recv_buf_up  , nbup  , &requp[1] );
    isend_buf( &send_buf_down, nbdown, &requp[0] );

    /* wait for atoms from down, move them to buffer cells */
    MPI_Waitall(2, reqdown, statdown);
    recv_buf_down.n = 0;
    for (i=0; i < cell_dim.x; ++i) 
      for (j=0; j < cell_dim.y; ++j)
	move_atoms_force( &recv_buf_down, i, j, cell_dim.z-1);

    /* wait for atoms from up, move them to buffer cells */
    MPI_Waitall(2, requp, statup);
    recv_buf_up.n = 0;
    for (i=0; i < cell_dim.x; ++i) 
      for (j=0; j < cell_dim.y; ++j)
        move_atoms_force( &recv_buf_up, i, j, 0);
  }
}


/******************************************************************************
*
* send_atoms_ar
*
* This sends the atoms in surface cells
* back to the neighbouring cpus and copies them to
* the buffer cells
*
* The data sent is:
*
* Faces:   up, south, down, west
* Edges:   up-south, down-south, up-west, north-west, down-west, west-south
* Corners: down-north-west, down-south-west, up-north-west, up-west-south
*
* The data is received on the canonical 'opposite side' of the data sent.
*
* A processor's area:
*    
*  from top    from botton
*
*      N           N
*    |---|       |---|
*  W | U | E    E| D |W    the coordinates origin is in 
*    |---|       |---|     the upper, north, west corner (unw)
*      S           S
*
*  Labeling:
*
*  Faces:    north, west, south, east, up, down
*  Edges:    north-west, west-south, south-east, east-north
*            up-north, up-west, up-south, up-east
*            down-north, down-east, down-south, down-west
*  Corners:  up-north-west, up-west-south, up-south-east, up-east-north
*            down-north-east, down-east-south, down-south-west, down-west-north
*   
*  Abbreviatations are e.g. dne for down-north-east etc.
*
******************************************************************************/

void send_atoms_ar(void)
{
  int i,j,k,l;

  MPI_Status  stateast[2],  statwest[2];
  MPI_Status statnorth[2], statsouth[2];
  MPI_Status    statup[2],  statdown[2];

  MPI_Request  reqeast[2],   reqwest[2];
  MPI_Request reqnorth[2],  reqsouth[2];
  MPI_Request    requp[2],   reqdown[2];

  /* Clear buffers */
  empty_mpi_buffers();
  empty_buffer_cells();
  
  /* Exchange faces */
  /* Fill send buffers */
  /* up, down */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.y-1; ++j) {
      copy_atoms_force( &send_buf_up  , i, j,            1 );
      copy_atoms_force( &send_buf_down, i ,j, cell_dim.z-2 );
  }
  /* south, north */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.z-1; ++j) {
      copy_atoms_force( &send_buf_south, i, cell_dim.y-2, j );
      copy_atoms_force( &send_buf_north, i,            1, j );
  }
  /* west, east */
  for (i=1; i < cell_dim.y-1; ++i)
    for (j=1; j < cell_dim.z-1; ++j) {
      copy_atoms_force( &send_buf_west, cell_dim.x-2, i, j );
      copy_atoms_force( &send_buf_east,            1, i, j );  
  }

  /* Exchange faces */
  irecv_buf( &recv_buf_west , nbeast , &reqwest[1] );
  irecv_buf( &recv_buf_south, nbnorth, &reqsouth[1]);
  irecv_buf( &recv_buf_down , nbup   , &reqdown[1] );
  irecv_buf( &recv_buf_up   , nbdown , &requp[1]   );
  irecv_buf( &recv_buf_north, nbsouth, &reqnorth[1]);
  irecv_buf( &recv_buf_east , nbwest , &reqeast[1]);

  isend_buf( &send_buf_west , nbwest  , &reqwest[0] );
  isend_buf( &send_buf_south, nbsouth , &reqsouth[0]);
  isend_buf( &send_buf_down , nbdown  , &reqdown[0] );
  isend_buf( &send_buf_up   , nbup    , &requp[0]   );
  isend_buf( &send_buf_north, nbnorth , &reqnorth[0]);
  isend_buf( &send_buf_east , nbeast  , &reqeast[0] );

  MPI_Waitall(2, reqwest, statwest);
  recv_buf_west.n = 0;
  MPI_Waitall(2, reqsouth, statsouth);
  recv_buf_south.n = 0;
  MPI_Waitall(2, reqdown, statdown);
  recv_buf_down.n = 0;
  MPI_Waitall(2, requp, statup);
  recv_buf_up.n = 0;
  MPI_Waitall(2, reqnorth, statnorth);
  recv_buf_north.n = 0;
  MPI_Waitall(2, reqeast, stateast); 
  recv_buf_east.n = 0;

  /* Add forces from faces */
  /* up, down */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.y-1; ++j) {
      move_atoms_force( &recv_buf_up  , i, j, cell_dim.z-1 );
      move_atoms_force( &recv_buf_down, i, j,            0 );
  }
  /* south, north */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.z-1; ++j) {
      move_atoms_force( &recv_buf_south, i,            0, j);
      move_atoms_force( &recv_buf_north, i, cell_dim.y-1, j);
  }
  /* west, east */
  for (i=1; i < cell_dim.y-1; ++i)
    for (j=1; j < cell_dim.z-1; ++j) {
      move_atoms_force( &recv_buf_west,            0, i, j );      
      move_atoms_force( &recv_buf_east, cell_dim.x-1, i, j ); 
  }

  /* send edges Part 1 */
  /* mapping buffers to edges 

     buffer     send        recv
     west    -  up-south    down-north
     east    -  down-south  up-north
     north   -  up-west     down-east
     south   -  north-west  south-east
     up      -  down-west   up-east
     down    -  west-south  east-north
  */

  empty_mpi_buffers();
  /* north-west, west-south */
  for (i=1; i < cell_dim.z-1; ++i) {
    copy_atoms_force( &send_buf_south, cell_dim.x-2,            1, i);
    copy_atoms_force( &send_buf_down,  cell_dim.x-2, cell_dim.y-2, i);
  }
  /* up-west, down-west */
  for (i=1; i < cell_dim.y-1; ++i) {
    copy_atoms_force( &send_buf_north, cell_dim.x-2, i,            1);
    copy_atoms_force( &send_buf_up   , cell_dim.x-2, i, cell_dim.z-2);
  }
  /* up-south, down-south */
  for (i=1; i < cell_dim.x-1; ++i) {
    copy_atoms_force( &send_buf_west , i, cell_dim.y-2,            1);
    copy_atoms_force( &send_buf_east , i, cell_dim.y-2, cell_dim.z-2);
  }

  /* Exchange edges */
  irecv_buf( &recv_buf_west , nbdn , &reqwest[1] );
  irecv_buf( &recv_buf_east , nbun , &reqeast[1] );
  irecv_buf( &recv_buf_north, nbde , &reqnorth[1]);
  irecv_buf( &recv_buf_south, nbse , &reqsouth[1]);
  irecv_buf( &recv_buf_up   , nbue , &requp[1]   );
  irecv_buf( &recv_buf_down , nben , &reqdown[1] );

  isend_buf( &send_buf_west , nbus , &reqwest[0] );
  isend_buf( &send_buf_east , nbds , &reqeast[0] );
  isend_buf( &send_buf_north, nbuw , &reqnorth[0]);
  isend_buf( &send_buf_south, nbnw , &reqsouth[0]);
  isend_buf( &send_buf_up   , nbdw , &requp[0]   );
  isend_buf( &send_buf_down , nbws , &reqdown[0] );

  MPI_Waitall(2, reqwest , statwest);
  MPI_Waitall(2, reqeast , stateast);
  MPI_Waitall(2, reqnorth, statnorth);
  MPI_Waitall(2, reqsouth, statsouth);
  MPI_Waitall(2, requp   , statup);
  MPI_Waitall(2, reqdown , statdown);

  recv_buf_west.n  = 0;
  recv_buf_east.n  = 0;
  recv_buf_north.n = 0;
  recv_buf_south.n = 0;
  recv_buf_up.n    = 0;
  recv_buf_down.n  = 0;

  /* Add edges */
  /* south-east, east-north  */
  for (i=1; i < cell_dim.z-1; ++i) {
    move_atoms_force( &recv_buf_south, 0, cell_dim.y-1, i);
    move_atoms_force( &recv_buf_down , 0,            0, i);
  }
  /* down-east, up-east */
  for (i=1; i < cell_dim.y-1; ++i) {
    move_atoms_force( &recv_buf_north, 0, i, cell_dim.z-1);
    move_atoms_force( &recv_buf_up   , 0, i,            0);
  }
  /* down-north, up-north */
  for (i=1; i < cell_dim.x-1; ++i) {
    move_atoms_force( &recv_buf_west, i, 0, cell_dim.z-1);
    move_atoms_force( &recv_buf_east, i, 0,            0);
  }

  /* send edges Part 2  */
  /* mapping buffers to edges 

     buffer     send        recv
     west    -  down-north  up-south
     east    -  up-north    down-south
     north   -  south-east  north-west  
     south   -  east-north  west-south 
     up      -  up-east     down-west
     down    -  down-east   up-west
  */

  empty_mpi_buffers();
  /* down-north, up-north */
  for (i=1; i < cell_dim.x-1; ++i) {
    copy_atoms_force( &send_buf_west, i, 1, cell_dim.z-2);
    copy_atoms_force( &send_buf_east, i, 1,            1);
  }
  /* south-east, east-north */
  for (i=1; i < cell_dim.z-1; ++i) { 
    copy_atoms_force( &send_buf_north, 1, cell_dim.y-2, i); 
    copy_atoms_force( &send_buf_south, 1,            1, i); 
  }
  /* up-east, down-east */ 
  for (i=1; i < cell_dim.y-1; ++i) { 
    copy_atoms_force( &send_buf_down, 1,  i,  cell_dim.z-2); 
    copy_atoms_force( &send_buf_up  , 1,  i,             1); 
  }

  /* Exchange edges */
  irecv_buf( &recv_buf_west , nbus , &reqwest[1] );
  irecv_buf( &recv_buf_east , nbds , &reqeast[1] );
  irecv_buf( &recv_buf_north, nbnw , &reqnorth[1]); 
  irecv_buf( &recv_buf_south, nbws , &reqsouth[1]); 
  irecv_buf( &recv_buf_up   , nbdw , &requp[1]   ); 
  irecv_buf( &recv_buf_down , nbuw , &reqdown[1] ); 

  isend_buf( &send_buf_west , nbdn , &reqwest[0] );
  isend_buf( &send_buf_east , nbun , &reqeast[0] );
  isend_buf( &send_buf_north, nbse , &reqnorth[0]); 
  isend_buf( &send_buf_south, nben , &reqsouth[0]); 
  isend_buf( &send_buf_up   , nbue , &requp[0]   ); 
  isend_buf( &send_buf_down , nbde , &reqdown[0] ); 

  MPI_Waitall(2, reqwest , statwest);
  MPI_Waitall(2, reqeast , stateast);
  MPI_Waitall(2, reqnorth, statnorth);
  MPI_Waitall(2, reqsouth, statsouth);
  MPI_Waitall(2, requp   , statup);
  MPI_Waitall(2, reqdown , statdown);

  recv_buf_west.n  = 0;
  recv_buf_east.n  = 0;
  recv_buf_north.n = 0;
  recv_buf_south.n = 0;
  recv_buf_up.n    = 0;
  recv_buf_down.n  = 0;

  /* Add edges */
  /* down-north, up-north */
  for (i=1; i < cell_dim.x-1; ++i) {
    move_atoms_force( &recv_buf_west, i, cell_dim.y-1,            0);
    move_atoms_force( &recv_buf_east, i, cell_dim.y-1, cell_dim.z-1);
  }
  /* south-east, east-north */
  for (i=1; i < cell_dim.z-1; ++i) {
    move_atoms_force( &recv_buf_north, cell_dim.x-1,            0, i);
    move_atoms_force( &recv_buf_south, cell_dim.x-1, cell_dim.y-1, i);
  }
  /* up-east, down-east */
  for (i=1; i < cell_dim.y-1; ++i) {
    move_atoms_force( &recv_buf_down, cell_dim.x-1, i,            0); 
    move_atoms_force( &recv_buf_up  , cell_dim.x-1, i, cell_dim.z-1); 
  }

  /* send corners */
  /* mapping buffers - to corners

     buffer    send             recv 
     west   -  up-north-west    down-east-south
     east   -  up-west-south    down-north-east
     south  -  down-south-west  up-east-north
     north  -  down-west-north  up-south-east

  */

  empty_mpi_buffers();
  /* unw, uws, dsw, dwn */

  copy_atoms_force( &send_buf_west , 1,            1,            1 );
  copy_atoms_force( &send_buf_east , 1, cell_dim.y-2,            1 );
  copy_atoms_force( &send_buf_south, 1, cell_dim.y-2, cell_dim.z-2 );
  copy_atoms_force( &send_buf_north, 1,            1, cell_dim.z-2 );

  /* Exchange corners */
  irecv_buf( &recv_buf_west ,  nbdsw , &reqwest[1] );
  irecv_buf( &recv_buf_east ,  nbdwn , &reqeast[1] );
  irecv_buf( &recv_buf_south,  nbunw , &reqsouth[1]);
  irecv_buf( &recv_buf_north,  nbuws , &reqnorth[1]);

  isend_buf( &send_buf_west ,  nbuen , &reqwest[0] );
  isend_buf( &send_buf_east ,  nbuse , &reqeast[0] );
  isend_buf( &send_buf_south,  nbdes , &reqsouth[0]);
  isend_buf( &send_buf_north,  nbdne , &reqnorth[0]);

  MPI_Waitall(2, reqwest, statwest);
  recv_buf_west.n = 0;
  MPI_Waitall(2, reqeast, stateast);
  recv_buf_east.n = 0;
  MPI_Waitall(2, reqnorth, statnorth);
  recv_buf_north.n = 0;
  MPI_Waitall(2, reqsouth, statsouth);
  recv_buf_south.n = 0;

  /* add corners */
  move_atoms_force( &recv_buf_west ,cell_dim.x-1, cell_dim.y-1, cell_dim.z-1);
  move_atoms_force( &recv_buf_east ,cell_dim.x-1,            0, cell_dim.z-1);
  move_atoms_force( &recv_buf_south,cell_dim.x-1,            0, 0);
  move_atoms_force( &recv_buf_north,cell_dim.x-1, cell_dim.y-1, 0);

}


/******************************************************************************
*
*  add_forces
*
******************************************************************************/

void add_forces( msgbuf *b, int k, int l, int m )
 {
  int i;
  cell *from;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  for (i=0; i<from->n; ++i) {
    from->kraft X(i) += b->data[ b->n++ ];
    from->kraft Y(i) += b->data[ b->n++ ];
    from->kraft Z(i) += b->data[ b->n++ ];
#ifndef MONOLJ
    from->pot_eng[i] += b->data[ b->n++ ];
#endif
  }
  if (b->n_max <= b->n) error("Buffer overflow in add_forces.");
}


/******************************************************************************
*
*  copy_forces
*
******************************************************************************/

void copy_forces( msgbuf *b, int k, int l, int m)
 {
  int i;
  cell *from;
    
  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  for (i=0; i<from->n; ++i) {
    b->data[ b->n++ ] = from->kraft X(i);
    b->data[ b->n++ ] = from->kraft Y(i);
    b->data[ b->n++ ] = from->kraft Z(i);
#ifndef MONOLJ
    b->data[ b->n++ ] = from->pot_eng[i];
#endif
  }
  if (b->n_max <= b->n) error("Buffer overflow in copy_forces.");
}


/******************************************************************************
*
*  send_forces
*
* This sends the forces accumulated by atoms in buffer cells
* back to the original cpus of these atoms and adds this 
* forces to each atoms total.
*
* The data sent is:
*
* Faces:   up, south, down, west
* Edges:   up-south, down-south, up-west, north-west, down-west, west-south
* Corners: down-north-west, down-south-west, up-north-west, up-west-south
*
* The data is received on the canonical 'opposite side' of the data sent.
*
* A processor's area:
*    
*  from top    from botton
*
*      N           N
*    |---|       |---|
*  W | U | E    E| D |W    the coordinates origin is in 
*    |---|       |---|     the upper, north, west corner (unw)
*      S           S
*
*  Labeling:
*
*  Faces:    north, west, south, east, up, down
*  Edges:    north-west, west-south, south-east, east-north
*            up-north, up-west, up-south, up-east
*            down-north, down-east, down-south, down-west
*  Corners:  up-north-west, up-west-south, up-south-east, up-east-north
*            down-north-east, down-east-south, down-south-west, down-west-north
*   
*  Abbreviatations are e.g. dne for down-north-east etc.
*
******************************************************************************/

void send_forces(void)
{
  int i,j,k,l;

  MPI_Status  stateast[2],  statwest[2];
  MPI_Status statnorth[2], statsouth[2];
  MPI_Status    statup[2],  statdown[2];

  MPI_Request  reqeast[2],   reqwest[2];
  MPI_Request reqnorth[2],  reqsouth[2];
  MPI_Request    requp[2],   reqdown[2];

  /* Clear buffers */
  empty_mpi_buffers();

  /* Fill send buffers */
  /* up, down */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.y-1; ++j) {
      copy_forces( &send_buf_up  , i, j,            0 );
      copy_forces( &send_buf_down, i ,j, cell_dim.z-1 );
  }
  /* south, north */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.z-1; ++j) {
      copy_forces( &send_buf_south, i, cell_dim.y-1, j );
      copy_forces( &send_buf_north, i,            0, j );
  }
  /* west */
  for (i=1; i < cell_dim.y-1; ++i)
    for (j=1; j < cell_dim.z-1; ++j) {
      copy_forces( &send_buf_west, cell_dim.x-1, i, j );
      /* copy_forces( &send_buf_east, 0, i, j );  */
  }

  /* Exchange faces */
  irecv_buf( &recv_buf_west , nbeast , &reqwest[1] );
  irecv_buf( &recv_buf_south, nbnorth, &reqsouth[1]);
  irecv_buf( &recv_buf_down , nbup   , &reqdown[1] );
  irecv_buf( &recv_buf_up   , nbdown , &requp[1]   );
  irecv_buf( &recv_buf_north, nbsouth, &reqnorth[1]);
  /* irecv_buf( &recv_buf_east, nbwest, &reqeast[1] ); */

  isend_buf( &send_buf_west , nbwest  , &reqwest[0] );
  isend_buf( &send_buf_south, nbsouth , &reqsouth[0]);
  isend_buf( &send_buf_down , nbdown  , &reqdown[0] );
  isend_buf( &send_buf_up   , nbup    , &requp[0]   );
  isend_buf( &send_buf_north, nbnorth , &reqnorth[0]);
  /* isend_buf( &send_buf_east, nbeast, &reqeast[0]); */

  MPI_Waitall(2, reqwest, statwest);
  recv_buf_west.n = 0;
  MPI_Waitall(2, reqsouth, statsouth);
  recv_buf_south.n = 0;
  MPI_Waitall(2, reqdown, statdown);
  recv_buf_down.n = 0;
  MPI_Waitall(2, requp, statup);
  recv_buf_up.n = 0;
  MPI_Waitall(2, reqnorth, statnorth);
  recv_buf_north.n = 0;
  /* MPI_Waitall(2, reqeast, stateast); */
  recv_buf_east.n = 0;

  /* Add forces from faces */
  /* up, down */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.y-1; ++j) {
      add_forces( &recv_buf_up  , i, j, cell_dim.z-2 );
      add_forces( &recv_buf_down, i, j,            1 );
  }
  /* south, north */
  for (i=1; i < cell_dim.x-1; ++i)
    for (j=1; j < cell_dim.z-1; ++j) {
      add_forces( &recv_buf_south, i,            1, j);
      add_forces( &recv_buf_north, i, cell_dim.y-2, j);
  }
  /* west, east */
  for (i=1; i < cell_dim.y-1; ++i)
    for (j=1; j < cell_dim.z-1; ++j) {
      add_forces( &recv_buf_west,            1, i, j );      
     /* add_forces( &recv_buf_east, cell_dim.x-2, i, j ); */
  }

  /* send edges Part 1 */
  /* mapping buffers to edges 

     buffer     send        recv
     west    -  up-south    down-north
     east    -  down-south  up-north
     north   -  up-west     down-east
     south   -  north-west  south-east
     up      -  down-west   up-east
     down    -  west-south  east-north
  */

  /* fill buffers */
  empty_mpi_buffers();
  /* north-west, west-south */
  for (i=1; i < cell_dim.z-1; ++i) {
    copy_forces( &send_buf_south, cell_dim.x-1,            0, i);
    copy_forces( &send_buf_down,  cell_dim.x-1, cell_dim.y-1, i);
  }
  /* up-west, down-west */
  for (i=1; i < cell_dim.y-1; ++i) {
    copy_forces( &send_buf_north, cell_dim.x-1, i,            0);
    copy_forces( &send_buf_up   , cell_dim.x-1, i, cell_dim.z-1);
  }
  /* up-south, down-south */
  for (i=1; i < cell_dim.x-1; ++i) {
    copy_forces( &send_buf_west , i, cell_dim.y-1,            0);
    copy_forces( &send_buf_east , i, cell_dim.y-1, cell_dim.z-1);
  }

  /* Exchange edges */
  irecv_buf( &recv_buf_west , nbdn , &reqwest[1] );
  irecv_buf( &recv_buf_east , nbun , &reqeast[1] );
  irecv_buf( &recv_buf_north, nbde , &reqnorth[1]);
  irecv_buf( &recv_buf_south, nbse , &reqsouth[1]);
  irecv_buf( &recv_buf_up   , nbue , &requp[1]   );
  irecv_buf( &recv_buf_down , nben , &reqdown[1] );

  isend_buf( &send_buf_west , nbus , &reqwest[0] );
  isend_buf( &send_buf_east , nbds , &reqeast[0] );
  isend_buf( &send_buf_north, nbuw , &reqnorth[0]);
  isend_buf( &send_buf_south, nbnw , &reqsouth[0]);
  isend_buf( &send_buf_up   , nbdw , &requp[0]   );
  isend_buf( &send_buf_down , nbws , &reqdown[0] );

  MPI_Waitall(2, reqwest , statwest );
  MPI_Waitall(2, reqeast , stateast );
  MPI_Waitall(2, reqnorth, statnorth);
  MPI_Waitall(2, reqsouth, statsouth);
  MPI_Waitall(2, requp   , statup   );
  MPI_Waitall(2, reqdown , statdown );

  recv_buf_west.n  = 0;
  recv_buf_east.n  = 0;
  recv_buf_north.n = 0;
  recv_buf_south.n = 0;
  recv_buf_up.n    = 0;
  recv_buf_down.n  = 0;

  /* Add edges */
  /* south-east, east-north  */
  for (i=1; i < cell_dim.z-1; ++i) {
    add_forces( &recv_buf_south, 1, cell_dim.y-2, i);
    add_forces( &recv_buf_down , 1,            1, i);
  };
  /* down-east, up-east */
  for (i=1; i < cell_dim.y-1; ++i) {
    add_forces( &recv_buf_north, 1, i, cell_dim.z-2);
    add_forces( &recv_buf_up   , 1, i,            1);
  };
  /* down-north, up-north */
  for (i=1; i < cell_dim.x-1; ++i) {
    add_forces( &recv_buf_west, i, 1, cell_dim.z-2);
    add_forces( &recv_buf_east, i, 1,            1);
  };

  /* send edges Part 2  */
  /* mapping buffers to edges 

     buffer     send        recv
     west    -  down-north  up-south
     east    -  up-north    down-south
     north   -  south-east  north-west  
     south   -  east-north  west-south 
     up      -  up-east     down-west
     down    -  down-east   up-west
  */

  empty_mpi_buffers();
  /* down-north, up-north */
  for (i=1; i < cell_dim.x-1; ++i) {
    copy_forces( &send_buf_west, i, 0, cell_dim.z-1);
    copy_forces( &send_buf_east, i, 0,            0);
  };

  /* Exchange edges */
  irecv_buf( &recv_buf_west , nbus , &reqwest[1]);
  irecv_buf( &recv_buf_east , nbds , &reqeast[1]);
  isend_buf( &send_buf_west , nbdn , &reqwest[0]);
  isend_buf( &send_buf_east , nbun , &reqeast[0]);

  MPI_Waitall(2, reqwest , statwest);
  MPI_Waitall(2, reqeast , stateast);

  recv_buf_west.n  = 0;
  recv_buf_east.n  = 0;

  /* Add edges */
  /* down-north, up-north */
  for (i=1; i < cell_dim.x-1; ++i) {
    add_forces( &recv_buf_west, i, cell_dim.y-2,            1);
    add_forces( &recv_buf_east, i, cell_dim.y-2, cell_dim.z-2);
  }

  /* send corners */
  /* mapping buffers - to corners

     buffer    send             recv 
     west   -  up-north-west    down-east-south
     east   -  up-west-south    down-north-east
     south  -  down-south-west  up-east-north
     north  -  down-west-north  up-south-east

  */

  empty_mpi_buffers();

  /* unw, uws, dsw, dwn */
  copy_forces( &send_buf_west , cell_dim.x-1,            0,            0 );
  copy_forces( &send_buf_east , cell_dim.x-1, cell_dim.y-1,            0 );
  copy_forces( &send_buf_south, cell_dim.x-1, cell_dim.y-1, cell_dim.z-1 );
  copy_forces( &send_buf_north, cell_dim.x-1,            0, cell_dim.z-1 );

  /* Exchange corners */
  irecv_buf( &recv_buf_west , nbdes, &reqwest[1] );
  irecv_buf( &recv_buf_east , nbdne, &reqeast[1] );
  irecv_buf( &recv_buf_south, nbuen, &reqsouth[1]);
  irecv_buf( &recv_buf_north, nbuse, &reqnorth[1]);

  isend_buf( &send_buf_west , nbunw, &reqwest[0] );
  isend_buf( &send_buf_east , nbuws, &reqeast[0] );
  isend_buf( &send_buf_south, nbdsw, &reqsouth[0]);
  isend_buf( &send_buf_north, nbdwn, &reqnorth[0]);  

  MPI_Waitall(2, reqwest, statwest );
  recv_buf_west.n = 0;
  MPI_Waitall(2, reqeast, stateast );
  recv_buf_east.n = 0;
  MPI_Waitall(2, reqnorth, statnorth );
  recv_buf_north.n = 0;
  MPI_Waitall(2, reqsouth, statsouth );
  recv_buf_south.n = 0;

  /* add corners */
  add_forces( &recv_buf_west , 1, cell_dim.y-2, cell_dim.z-2);
  add_forces( &recv_buf_east , 1,            1, cell_dim.z-2);
  add_forces( &recv_buf_south, 1,            1, 1);
  add_forces( &recv_buf_north, 1, cell_dim.y-2, 1);

}


/******************************************************************************
*
*  add_cell_force  -  add forces of one cell to those of another cell
*
******************************************************************************/

void add_cell_force( int k, int l, int m, int r, int s, int t)
{
  int i;
  cell *from, *to;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  to   = PTR_3D_V(cell_array, r, s, t, cell_dim);

  for (i=0; i<to->n; ++i) {
    to->kraft X(i) += from->kraft X(i);
    to->kraft Y(i) += from->kraft Y(i);
    to->kraft Z(i) += from->kraft Z(i);
#ifndef MONOLJ
    to->pot_eng[i] += from->pot_eng[i];
#endif
  }
}


/******************************************************************************
*
*  send_forces_full  -  send forces of ALL buffer cells
*
******************************************************************************/

void send_forces_full(void)
{
  int i,j;

  MPI_Status  stateast[2],  statwest[2];
  MPI_Status statnorth[2], statsouth[2];
  MPI_Status    statup[2],  statdown[2];

  MPI_Request  reqeast[2],   reqwest[2];
  MPI_Request reqnorth[2],  reqsouth[2];
  MPI_Request    requp[2],   reqdown[2];

  empty_mpi_buffers();

  /* send forces up/down */
  if (cpu_dim.z==1) {
    /* simply add up/down forces to original cells */
    for (i=0; i < cell_dim.x; ++i) 
      for (j=0; j < cell_dim.y; ++j) {
	add_cell_force( i, j, 0, i, j, cell_dim.z-2 );
        add_cell_force( i, j, cell_dim.z-1, i, j, 1 );
      }
  } else {
    /* copy up forces into send buffer, send up */
    for (i=0; i < cell_dim.x; ++i) 
      for (j=0; j < cell_dim.y; ++j)
	copy_forces( &send_buf_up, i, j, 0);
    irecv_buf( &recv_buf_down , nbdown, &reqdown[1]);
    isend_buf( &send_buf_up   , nbup  , &reqdown[0]);

    /* copy down forces into send buffer, send down */
    for (i=0; i < cell_dim.x; ++i) 
      for (j=0; j < cell_dim.y; ++j)
        copy_forces( &send_buf_down, i, j, cell_dim.z-1);
    irecv_buf( &recv_buf_up  , nbup  , &requp[1] );
    isend_buf( &send_buf_down, nbdown, &requp[0] );

    /* wait for forces from down, add them to original cells */
    MPI_Waitall(2, reqdown, statdown);
    recv_buf_down.n = 0;
    for (i=0; i < cell_dim.x; ++i) 
      for (j=0; j < cell_dim.y; ++j)
	add_forces( &recv_buf_down, i, j, cell_dim.z-2);

    /* wait for forces from up, add them to original cells */
    MPI_Waitall(2, requp, statup);
    recv_buf_up.n = 0;
    for (i=0; i < cell_dim.x; ++i) 
      for (j=0; j < cell_dim.y; ++j)
        add_forces( &recv_buf_up, i, j, 1);
  }

  /* send forces north/south */
  if (cpu_dim.y==1) {
    /* simply add north/south forces to original cells */
    for (i=0; i < cell_dim.x; ++i)
      for (j=1; j < cell_dim.z-1; ++j) {
        add_cell_force( i, 0, j, i, cell_dim.y-2, j );
        add_cell_force( i, cell_dim.y-1, j, i, 1, j );
      }
  } else {
    /* copy north forces into send buffer, send north */
    for (i=0; i < cell_dim.x; ++i) 
      for (j=1; j < cell_dim.z-1; ++j)
	copy_forces( &send_buf_north, i, 0, j);
    irecv_buf( &recv_buf_south, nbsouth, &reqsouth[1] );
    isend_buf( &send_buf_north, nbnorth, &reqsouth[0] );

    /* copy south forces into send buffer, send south */
    for (i=0; i < cell_dim.x; ++i) 
      for (j=1; j < cell_dim.z-1; ++j)
        copy_forces( &send_buf_south, i, cell_dim.y-1, j);
    irecv_buf( &recv_buf_north, nbnorth, &reqnorth[1] );
    isend_buf( &send_buf_south, nbsouth, &reqnorth[0] );

    /* wait for forces from south, add them to original cells */
    MPI_Waitall(2, reqsouth, statsouth);
    recv_buf_south.n = 0;
    for (i=0; i < cell_dim.x; ++i) 
      for (j=1; j < cell_dim.z-1; ++j)
	add_forces( &recv_buf_south, i, cell_dim.y-2, j);

    /* wait for forces from north, add them to original cells */
    MPI_Waitall(2, reqnorth, statnorth);
    recv_buf_north.n = 0;
    for (i=0; i < cell_dim.x; ++i) 
      for (j=1; j < cell_dim.z-1; ++j)
        add_forces( &recv_buf_north, i, 1, j);
  }

  /* send forces east/west */
  if (cpu_dim.x==1) {
    /* simply add east/west forces to original cells */
    for (i=1; i < cell_dim.y-1; ++i)
      for (j=1; j < cell_dim.z-1; ++j) { 
        add_cell_force( 0, i, j, cell_dim.x-2, i, j );
        add_cell_force( cell_dim.x-1, i, j, 1, i, j );
      }
  } else {
    /* copy east forces into send buffer, send east */
    for (i=1; i < cell_dim.y-1; ++i)
      for (j=1; j < cell_dim.z-1; ++j) 
	copy_forces( &send_buf_east, 0, i, j);
    irecv_buf( &recv_buf_west, nbwest, &reqwest[1] );
    isend_buf( &send_buf_east, nbeast, &reqwest[0] );

    /* copy west forces into send buffer, send west */
    for (i=1; i < cell_dim.y-1; ++i) 
      for (j=1; j < cell_dim.z-1; ++j) 
	copy_forces( &send_buf_west, cell_dim.x-1, i, j);
    irecv_buf( &recv_buf_east, nbeast, &reqeast[1] );
    isend_buf( &send_buf_west, nbwest, &reqeast[0] );

    /* wait for forces from west, add them to original cells */
    MPI_Waitall(2, reqwest, statwest);
    recv_buf_west.n = 0;
    for (i=1; i < cell_dim.y-1; ++i)
      for (j=1; j < cell_dim.z-1; ++j) 
	add_forces( &recv_buf_west, cell_dim.x-2, i, j);

    /* wait for forces from east, add them to original cells */
    MPI_Waitall(2, reqeast, stateast);
    recv_buf_east.n = 0;
    for (i=1; i < cell_dim.y-1; ++i) 
      for (j=1; j < cell_dim.z-1; ++j) 
	add_forces( &recv_buf_east, 1, i, j);
  }
}

#ifdef EAM2

/******************************************************************************
*
*  copy_cell_eam2_rho_h  -  copy eam2_rho_h of one cell to another cell
*
******************************************************************************/

void copy_cell_eam2_rho_h( int k, int l, int m, int r, int s, int t)
{
  int i;
  cell *from, *to;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);
  to   = PTR_3D_V(cell_array, r, s, t, cell_dim);

  for (i=0; i<to->n; ++i) {
    to->eam2_rho_h[i] = from->eam2_rho_h[i];
  }
}


/******************************************************************************
*
*  move_eam2_rho_h  -  move eam2_rho_h from MPI buffer to cells
*
******************************************************************************/

void move_eam2_rho_h( msgbuf *b, int k, int l, int m )
{
  int i;
  cell *to;

  to = PTR_3D_V(cell_array, k, l, m, cell_dim);

  for (i=0; i<to->n; ++i) {
    to->eam2_rho_h[i] = b->data[ b->n++ ];
  }
}


/******************************************************************************
*
*  copy_eam2_rho_h  -  copy eam2_rho_h to MPI buffer
*
******************************************************************************/

void copy_eam2_rho_h( msgbuf *b, int k, int l, int m)
{
  int i;
  cell *from;
    
  from = PTR_3D_V(cell_array, k, l, m, cell_dim);

  for (i=0; i<from->n; ++i) {
    b->data[ b->n++ ] = from->eam2_rho_h[i];
  }
}


/******************************************************************************
*
* send_eam2_rho_h
*
* This sends eam2_rho_h to the neighbouring cpus using
* Steve Plimptons comm scheme
*
* Data is sent only over the faces of the processors box in 
* the order east-west, north-south, up-down
*
******************************************************************************/

void send_eam2_rho_h()
{
  int i,j;

  MPI_Status  stateast[2],  statwest[2];
  MPI_Status statnorth[2], statsouth[2];
  MPI_Status    statup[2],  statdown[2];

  MPI_Request  reqeast[2],   reqwest[2];
  MPI_Request reqnorth[2],  reqsouth[2];
  MPI_Request    requp[2],   reqdown[2];

  empty_mpi_buffers();

  /* exchange east/west */
  if (cpu_dim.x==1) {
    /* simply copy east/west atoms to buffer cells */
    for (i=1; i < cell_dim.y-1; ++i)
      for (j=1; j < cell_dim.z-1; ++j) { 
        copy_cell_eam2_rho_h( 1, i, j, cell_dim.x-1, i, j );
        copy_cell_eam2_rho_h( cell_dim.x-2, i, j, 0, i, j );
      }
  } else {
    /* copy east atoms into send buffer, send east */
    for (i=1; i < cell_dim.y-1; ++i)
      for (j=1; j < cell_dim.z-1; ++j) 
	copy_eam2_rho_h( &send_buf_east, 1, i, j);
    irecv_buf( &recv_buf_west, nbwest, &reqwest[1] );
    isend_buf( &send_buf_east, nbeast, &reqwest[0] );

    /* copy west atoms into send buffer, send west */
    for (i=1; i < cell_dim.y-1; ++i) 
      for (j=1; j < cell_dim.z-1; ++j) 
	copy_eam2_rho_h( &send_buf_west, cell_dim.x-2, i, j);
    irecv_buf( &recv_buf_east, nbeast, &reqeast[1] );
    isend_buf( &send_buf_west, nbwest, &reqeast[0] );

    /* wait for atoms from west, move them to buffer cells */
    MPI_Waitall(2, reqwest, statwest);
    recv_buf_west.n = 0;
    for (i=1; i < cell_dim.y-1; ++i)
      for (j=1; j < cell_dim.z-1; ++j) 
	move_eam2_rho_h( &recv_buf_west, cell_dim.x-1, i, j);

    /* wait for atoms from east, move them to buffer cells */
    MPI_Waitall(2, reqeast, stateast);
    recv_buf_east.n = 0;
    for (i=1; i < cell_dim.y-1; ++i) 
      for (j=1; j < cell_dim.z-1; ++j) 
	move_eam2_rho_h( &recv_buf_east, 0, i, j);
  }
  
  /* exchange north/south */
  if (cpu_dim.y==1) {
    /* simply copy north/south atoms to buffer cells */
    for (i=0; i < cell_dim.x; ++i)
      for (j=1; j < cell_dim.z-1; ++j) {
        copy_cell_eam2_rho_h( i, 1, j, i, cell_dim.y-1, j );
        copy_cell_eam2_rho_h( i, cell_dim.y-2, j, i, 0, j );
      }
  } else {
    /* copy north atoms into send buffer, send north */
    for (i=0; i < cell_dim.x; ++i) 
      for (j=1; j < cell_dim.z-1; ++j)
	copy_eam2_rho_h( &send_buf_north, i, 1, j);
    irecv_buf( &recv_buf_south, nbsouth, &reqsouth[1] );
    isend_buf( &send_buf_north, nbnorth, &reqsouth[0] );

    /* copy south atoms into send buffer, send south */
    for (i=0; i < cell_dim.x; ++i) 
      for (j=1; j < cell_dim.z-1; ++j)
        copy_eam2_rho_h( &send_buf_south, i, cell_dim.y-2, j);
    irecv_buf( &recv_buf_north, nbnorth, &reqnorth[1] );
    isend_buf( &send_buf_south, nbsouth, &reqnorth[0] );

    /* wait for atoms from south, move them to buffer cells */
    MPI_Waitall(2, reqsouth, statsouth);
    recv_buf_south.n = 0;
    for (i=0; i < cell_dim.x; ++i) 
      for (j=1; j < cell_dim.z-1; ++j)
	move_eam2_rho_h( &recv_buf_south, i, cell_dim.y-1, j);

    /* wait for atoms from north, move them to buffer cells */
    MPI_Waitall(2, reqnorth, statnorth);
    recv_buf_north.n = 0;
    for (i=0; i < cell_dim.x; ++i) 
      for (j=1; j < cell_dim.z-1; ++j)
        move_eam2_rho_h( &recv_buf_north, i, 0, j);
  }

  /* exchange up/down */
  if (cpu_dim.z==1) {
    /* simply copy up/down atoms to buffer cells */
    for (i=0; i < cell_dim.x; ++i) 
      for (j=0; j < cell_dim.y; ++j) {
	copy_cell_eam2_rho_h( i, j, 1, i, j, cell_dim.z-1 );
        copy_cell_eam2_rho_h( i, j, cell_dim.z-2, i, j, 0 );
      }
  } else {
    /* copy up atoms into send buffer, send up */
    for (i=0; i < cell_dim.x; ++i) 
      for (j=0; j < cell_dim.y; ++j)
	copy_eam2_rho_h( &send_buf_up, i, j, 1);
    irecv_buf( &recv_buf_down , nbdown, &reqdown[1]);
    isend_buf( &send_buf_up   , nbup  , &reqdown[0]);

    /* copy down atoms into send buffer, send down */
    for (i=0; i < cell_dim.x; ++i) 
      for (j=0; j < cell_dim.y; ++j)
        copy_eam2_rho_h( &send_buf_down, i, j, cell_dim.z-2);
    irecv_buf( &recv_buf_up  , nbup  , &requp[1] );
    isend_buf( &send_buf_down, nbdown, &requp[0] );

    /* wait for atoms from down, move them to buffer cells */
    MPI_Waitall(2, reqdown, statdown);
    recv_buf_down.n = 0;
    for (i=0; i < cell_dim.x; ++i) 
      for (j=0; j < cell_dim.y; ++j)
	move_eam2_rho_h( &recv_buf_down, i, j, cell_dim.z-1);

    /* wait for atoms from up, move them to buffer cells */
    MPI_Waitall(2, requp, statup);
    recv_buf_up.n = 0;
    for (i=0; i < cell_dim.x; ++i) 
      for (j=0; j < cell_dim.y; ++j)
        move_eam2_rho_h( &recv_buf_up, i, j, 0);
  }
}

#endif /* EAM2 */
