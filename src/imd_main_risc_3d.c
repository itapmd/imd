/******************************************************************************
*
* imd_main_risc_3d.c -- main loop, risc specific part, three dimensions
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"


/*****************************************************************************
*
* calc_forces()
*
*****************************************************************************/

void calc_forces(void)
{
  cell *p,*q;
  int i,j,k;
  int l,m,n;
  int r,s,t;
  vektor pbc;

  tot_pot_energy = 0.0;
  virial         = 0.0;
#ifdef P_AXIAL
  vir_vect.x     = 0.0;
  vir_vect.y     = 0.0;
  vir_vect.z     = 0.0;
#endif
  
#ifdef EAM
  memset(eam_rho,   0, natoms*        sizeof(real));
  memset(eam_ij,    0, natoms*eam_len*sizeof(integer));
  memset(eam_dij_x, 0, natoms*eam_len*sizeof(real));
  memset(eam_dij_y, 0, natoms*eam_len*sizeof(real));
  memset(eam_dij_z, 0, natoms*eam_len*sizeof(real));
#endif /* EAM */

#ifdef TTBP
  memset(ttbp_ij,   0, natoms    *ttbp_len*2*sizeof(integer));
  memset(ttbp_j,    0, (natoms+1)*ttbp_len  *sizeof(real));
  memset(ttbp_force,0, (natoms+1)         *3*sizeof(real));
#endif /* TTBP */

  /* Zero Forces and potential energy */
  for (p = cell_array; 
       p <= PTR_3D_V(cell_array,
		     cell_dim.x-1,
		     cell_dim.y-1,
		     cell_dim.z-1,
		     cell_dim);
       ++p ) 
    for (i = 0;i < p->n; ++i) {
      p->kraft X(i) = 0.0;
      p->kraft Y(i) = 0.0;
      p->kraft Z(i) = 0.0;
#ifdef TRANSPORT
      p->heatcond[i] = 0.0;
#endif     
#ifndef MONOLJ
      p->pot_eng[i] = 0.0;
#endif
    };

  /* for each cell */
  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
      for (k=0; k < cell_dim.z; ++k)

	/* For half of the neighbours of this cell */
	for (l=0; l <= 1; ++l)
	  for (m=-l; m <= 1; ++m)
	    for (n=(l==0 ? -m  : -l ); n <= 1; ++n) { 

	      /* Given cell */
              /* Hey optimizer, this is invariant to the last three loops!! */
	      p = PTR_3D_V(cell_array,i,j,k,cell_dim);
	      /* Calculate Indicies of Neighbour */
	      r = i+l;
	      s = j+m;
	      t = k+n;
	      /* Apply periodic boundaries */
	      pbc.x = 0;
	      pbc.y = 0;
	      pbc.z = 0;

	      if (r<0) {
		r = cell_dim.x-1; 
		pbc.x -= box_x.x;      
		pbc.y -= box_x.y;
		pbc.z -= box_x.z;
	      };
	      if (s<0) {
		s = cell_dim.y-1;
		pbc.x -= box_y.x;      
		pbc.y -= box_y.y;
		pbc.z -= box_y.z;
	      };
	      if (t<0) {
		t = cell_dim.z-1;
		pbc.x -= box_z.x;      
		pbc.y -= box_z.y;
		pbc.z -= box_z.z;
	      };
	      if (r>cell_dim.x-1) {
		r = 0; 
		pbc.x += box_x.x;      
		pbc.y += box_x.y;
		pbc.z += box_x.z;
	      };
	      if (s>cell_dim.y-1) {
		s = 0; 
		pbc.x += box_y.x;      
		pbc.y += box_y.y;
		pbc.z += box_y.z;
	      };
	      if (t>cell_dim.z-1) {
		t = 0; 
		pbc.x += box_z.x;      
		pbc.y += box_z.y;
		pbc.z += box_z.z;
	      };

	      /* Neighbour (note that p==q ist possible) */
	      q = PTR_3D_V(cell_array,r,s,t,cell_dim);
	      /* Do the work */

#ifdef SHOCK
              if (0 == pbc.x)
#endif
#ifdef NOPBC
              if ((0 == pbc.x) && (0 == pbc.y) && (0 == pbc.z))
#endif
#ifdef EAM
	      do_forces_eam_1(p,q,pbc);
#elif TTBP
	      do_forces_ttbp_1(p,q,pbc);
#else
	      do_forces(p,q,pbc);
#endif /* EAM TTBP classical */

      };

#ifdef EAM
  /* EAM cohesive function potential: for each cell */
  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
      for (k=0; k < cell_dim.z; ++k) {
	      /* Given cell */
              /* Hey optimizer, this is invariant to the last three loops!! */
	      p = PTR_3D_V(cell_array,i,j,k,cell_dim);
	      pbc.x = 0;
	      pbc.y = 0;
	      pbc.z = 0;
	      /* Neighbour (dummy; p==q) */
              q = p;
	      /* Do the work */
#ifdef SHOCK
              if (0 == pbc.x)
#endif
#ifdef NOPBC
              if ((0 == pbc.x) && (0 == pbc.y) && (0 == pbc.z))
#endif
	      do_forces_eam_2(p,q,pbc);
      };
#endif /* EAM */

#ifdef TTBP
  /* Part II. TTBP: three body potential */
  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
      for (k=0; k < cell_dim.z; ++k) {
	      /* Given cell */
              /* Hey optimizer, this is invariant to the last three loops!! */
	      p = PTR_3D_V(cell_array,i,j,k,cell_dim);
	      pbc.x = 0;
	      pbc.y = 0;
	      pbc.z = 0;
	      /* Neighbour (dummy; p==q) */
              q = p;
	      /* Do the work */
#ifdef SHOCK
              if (0 == pbc.x)
#endif
#ifdef NOPBC
              if ((0 == pbc.x) && (0 == pbc.y) && (0 == pbc.z))
#endif
	      do_forces_ttbp_2(p,q,pbc);
      };

  /* Part III. TTBP: three body potential */
  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
      for (k=0; k < cell_dim.z; ++k) {
	      /* Given cell */
              /* Hey optimizer, this is invariant to the last three loops!! */
	      p = PTR_3D_V(cell_array,i,j,k,cell_dim);
	      pbc.x = 0;
	      pbc.y = 0;
	      pbc.z = 0;
	      /* Neighbour (dummy; p==q) */
              q = p;
	      /* Do the work */
#ifdef SHOCK
              if (0 == pbc.x)
#endif
#ifdef NOPBC
              if ((0 == pbc.x) && (0 == pbc.y) && (0 == pbc.z))
#endif
	      do_forces_ttbp_3(p,q,pbc);
      };
#endif /* TTBP */

}


