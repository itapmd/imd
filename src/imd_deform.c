/******************************************************************************
*
* imd_deform.c -- deform sample
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

#ifdef HOMDEF   /* homogeneous deformation with pbc */

/*****************************************************************************
*
* expand_sample()
*
*****************************************************************************/

void expand_sample(void)
{
  int i,r,s,t;
  ivektor max_cell_dim;
  cell *p;
  
  /* Apply expansion */
  for ( r = cellmin.x; r < cellmax.x; ++r )
    for ( s = cellmin.y; s < cellmax.y; ++s )
#ifndef TWOD
      for ( t = cellmin.z; t < cellmax.z; ++t )
#endif
      {      
#ifdef TWOD
        p = PTR_2D_V(cell_array, r, s, cell_dim);
#else
        p = PTR_3D_V(cell_array, r, s, t, cell_dim);
#endif
        for (i = 0;i < p->n; ++i) {
          p->ort X(i) *= expansion.x;
          p->ort Y(i) *= expansion.y;
#ifndef TWOD
          p->ort Z(i) *= expansion.z;
#endif
        }
      }
  /* new box size */
#ifdef TWOD
  box_x.x *= expansion.x;  box_y.x *= expansion.x;
  box_x.y *= expansion.y;  box_y.y *= expansion.y;
#else
  box_x.x *= expansion.x;  box_x.y *= expansion.y;  box_x.z *= expansion.z;
  box_y.x *= expansion.x;  box_y.y *= expansion.y;  box_y.z *= expansion.z;
  box_z.x *= expansion.x;  box_z.y *= expansion.y;  box_z.z *= expansion.z;
#endif
  make_box();

  /* revise cell decomposition if necessary */
  max_cell_dim = maximal_cell_dim();
  if ((max_cell_dim.x<global_cell_dim.x) || (max_cell_dim.y<global_cell_dim.y)
#ifndef TWOD
      || (max_cell_dim.z<global_cell_dim.z)
#endif
  ) {
    init_cells();
    fix_cells();
  }

} /* expand sample */


/*****************************************************************************
*
* shear_sample()
*
*****************************************************************************/

void shear_sample(void)
{
  int i,r,s,t;
  ivektor max_cell_dim;
  cell *p;

  /* Apply shear */
  for ( r = cellmin.x; r < cellmax.x; ++r )
    for ( s = cellmin.y; s < cellmax.y; ++s )
#ifndef TWOD
      for ( t = cellmin.z; t < cellmax.z; ++t )
#endif
      {
#ifdef TWOD
        p = PTR_2D_V(cell_array, r, s, cell_dim);
#else
        p = PTR_3D_V(cell_array, r, s, t, cell_dim);
#endif
        for (i = 0; i < p->n; ++i)
          p->ort Y(i) += shear_factor * p->ort X(i);
      }

  /* new box size */
  box_x.y += shear_factor * box_x.x;
  make_box();

  /* revise cell decomposition if necessary */
  max_cell_dim = maximal_cell_dim();
  if ((max_cell_dim.x<global_cell_dim.x) || (max_cell_dim.y<global_cell_dim.y)
#ifndef TWOD
      || (max_cell_dim.z<global_cell_dim.z)
#endif
  ) {
    init_cells();
    fix_cells();
  }

} /* shear sample */

#endif /* HOMDEF */


#ifdef DEFORM

/*****************************************************************************
*
* deform_sample()
*
*****************************************************************************/

void deform_sample(void) {

  cell *p;
  int i;
  int r,s,t;
  real box_x_half;

  box_x_half = 0.5 * box_x.x;

  /* loop over all atoms */
    for ( r = cellmin.x; r < cellmax.x; ++r )
      for ( s = cellmin.y; s < cellmax.y; ++s )
#ifndef TWOD
	for ( t = cellmin.z; t < cellmax.z; ++t )
#endif
	{

#ifndef TWOD
	  p = PTR_3D_V(cell_array, r, s, t, cell_dim);
#else
	  p = PTR_2D_V(cell_array, r, s,    cell_dim);
#endif

	  for (i = 0;i < p->n; ++i) {
            /* move only atoms with negative number */
            if (NUMMER(p,i) > 0) continue;
	    /* which direction of pulling? */
	    if (p->ort X(i) <= box_x_half) {
	      p->ort X(i) += strip_shift.x;
	      p->ort Y(i) += strip_shift.y;
#ifndef TWOD
	      p->ort Z(i) += strip_shift.z;
#endif
	    } else {
	      p->ort X(i) -= strip_shift.x;
	      p->ort Y(i) -= strip_shift.y;
#ifndef TWOD
	      p->ort Z(i) -= strip_shift.z;
#endif        
	    }
	  } /* i - loop */
	} /* cell loop */

} /* deform_atoms */

#endif /* DEFORM */

