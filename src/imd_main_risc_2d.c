/******************************************************************************
*
* imd_main_risc.c -- main loop, risc specific part, two dimensions
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
*  calc_forces()
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
#endif

  /* Zero Forces */
  for (p = cell_array; 
       p <= PTR_2D_V(cell_array,
		     cell_dim.x-1,
		     cell_dim.y-1,
		     cell_dim);
       ++p ) 
    for (i = 0;i < p->n; ++i) {
      p->kraft X(i) = 0.0;
      p->kraft Y(i) = 0.0;
      p->pot_eng[i] = 0.0;
      
    };

  /* for each cell */
  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)

	/* For half of the neighbours of this cell */
	for (l=0; l <= 1; ++l)
	  for (m=-l; m <= 1; ++m) {

	      /* Given cell */
              /* Hey optimizer, this is invariant to the last three loops!! */
	      p = PTR_2D_V(cell_array,i,j,cell_dim);
	      /* Calculate Indicies of Neighbour */
	      r = i+l;
	      s = j+m;
	      /* Apply periodic boundaries */
	      pbc.x = 0;
	      pbc.y = 0;
#ifndef NOPBC
	      if (r<0) {
		r = cell_dim.x-1; 
		pbc.x -= box_x.x;      
		pbc.y -= box_x.y;
	      };
	      if (s<0) {
		s = cell_dim.y-1;
		pbc.x -= box_y.x;      
		pbc.y -= box_y.y;
	      };
	      if (r>cell_dim.x-1) {
		r = 0; 
		pbc.x += box_x.x;      
		pbc.y += box_x.y;
	      };
	      if (s>cell_dim.y-1) {
		s = 0; 
		pbc.x += box_y.x;      
		pbc.y += box_y.y;
	      };
#endif

#ifdef SHOCK
              if (0 == pbc.x)
#endif
#ifdef NOPBC
              if ((0 == pbc.x) && (0 == pbc.y))
#endif
	      {
                /* Neighbour (note that p==q ist possible) */
                q = PTR_2D_V(cell_array,r,s,cell_dim);

                /* Do the work */
                do_forces(p,q,pbc);
              }
      }
}



