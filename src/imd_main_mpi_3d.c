
/******************************************************************************
*
* imd_main_mpi_3d.c -- main loop, mpi specific part, three dimensions
*
******************************************************************************/

/******************************************************************************
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
    for (i=0; i<p->n; ++i) {
      p->kraft X(i) = 0.0;
      p->kraft Y(i) = 0.0;
      p->kraft Z(i) = 0.0;
#ifdef UNIAX
      p->dreh_moment X(i) = 0.0;
      p->dreh_moment Y(i) = 0.0;
      p->dreh_moment Z(i) = 0.0;
#endif
#ifdef NVX
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
#ifdef ORDPAR
      p->nbanz[i] = 0;
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

  send_cells(copy_rho_h,pack_rho_h,unpack_rho_h);

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
  MPI_Allreduce( &tot_pot_energy, &tmp, 1, REAL, MPI_SUM, cpugrid); 
  tot_pot_energy = tmp; 

  MPI_Allreduce( &virial,   &tmpvir,    1, REAL, MPI_SUM, cpugrid);
  virial = tmpvir;

  MPI_Allreduce( &vir_vect, &tmpvec,  DIM, REAL, MPI_SUM, cpugrid);
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

#endif  /* ... ifndef AR */

  /* sum up results of different CPUs */
  MPI_Allreduce( &tot_pot_energy, &tmp, 1, REAL, MPI_SUM, cpugrid); 
  tot_pot_energy = tmp; 

  MPI_Allreduce( &virial, &tmp,         1, REAL, MPI_SUM, cpugrid);
  virial = tmp;

#ifdef P_AXIAL
  MPI_Allreduce( &vir_vect, &tmpvec,  DIM, REAL, MPI_SUM, cpugrid);
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
* send_atoms
*
* This sends the atom positions to the neighbouring cpus using
* Steve Plimptons comm scheme
*
* Atoms are sent only over the faces of the processors box in 
* the order east-west, north-south, up-down
*
* This is only used for atom redistribuion
* where all of the data are sent.
*
* The drawback of the Plimpton scheme are the huge comm-buffers,
* that use lots of memory. It is fast, however!!
*
******************************************************************************/

void send_atoms()
{
  MPI_Status  stateast[2],  statwest[2];
  MPI_Status statnorth[2], statsouth[2];
  MPI_Status    statup[2],  statdown[2];

  MPI_Request  reqeast[2],   reqwest[2];
  MPI_Request reqnorth[2],  reqsouth[2];
  MPI_Request    requp[2],   reqdown[2];

  if (cpu_dim.x > 1) {
    /* send east */
    irecv_buf( &recv_buf_west, nbwest, &reqwest[1] );
    isend_buf( &send_buf_east, nbeast, &reqwest[0] );

    /* send west */
    irecv_buf( &recv_buf_east, nbeast, &reqeast[1] );
    isend_buf( &send_buf_west, nbwest, &reqeast[0] );

    /* Wait for atoms from west, move them to cells */
    MPI_Waitall(2, reqwest, statwest);
    MPI_Get_count( &statwest[1], REAL, &recv_buf_west.n );
    process_buffer( &recv_buf_west, (cell *) NULL );

    /* Wait for atoms from east, move them to cells */
    MPI_Waitall(2, reqeast, stateast);
    MPI_Get_count( &stateast[1], REAL, &recv_buf_east.n );
    process_buffer( &recv_buf_east, (cell *) NULL );

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
  }

  if (cpu_dim.y > 1) {
    /* Send atoms north */
    irecv_buf( &recv_buf_south, nbsouth, &reqsouth[1] );
    isend_buf( &send_buf_north, nbnorth, &reqsouth[0] );
  
    /* Send atoms south */
    irecv_buf( &recv_buf_north, nbnorth, &reqnorth[1] );
    isend_buf( &send_buf_south, nbsouth, &reqnorth[0] );

    /* Wait for atoms from south, move them to cells */
    MPI_Waitall(2, reqsouth, statsouth);
    MPI_Get_count( &statsouth[1], REAL, &recv_buf_south.n );
    process_buffer( &recv_buf_south, (cell *) NULL );

    /* Wait for atoms from north, move them to cells */
    MPI_Waitall(2, reqnorth, statnorth);
    MPI_Get_count( &statnorth[1], REAL, &recv_buf_north.n );
    process_buffer( &recv_buf_north, (cell *) NULL );

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
  }

  if (cpu_dim.z > 1) {
    /* Send atoms up */
    irecv_buf( &recv_buf_down , nbdown, &reqdown[1]);
    isend_buf( &send_buf_up   , nbup  , &reqdown[0]);
  
    /* Send atoms down */
    irecv_buf( &recv_buf_up  , nbup  , &requp[1] );
    isend_buf( &send_buf_down, nbdown, &requp[0] );

    /* Wait for atoms from down, move them to cells */
    MPI_Waitall(2, reqdown, statdown);
    MPI_Get_count( &statdown[1], REAL, &recv_buf_down.n );
    process_buffer( &recv_buf_down, (cell *) NULL );

    /* Wait for atoms from up, move them to cells */
    MPI_Waitall(2, requp, statup);
    MPI_Get_count( &statup[1], REAL, &recv_buf_up.n );
    process_buffer( &recv_buf_up, (cell *) NULL );
  }

}



