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
  int k;
  ivektor max_cell_dim;
  
  /* Apply expansion */
#pragma omp parallel for
  for (k=0; k<ncells; ++k) {
    int i;
    cell *p;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) {
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
  int k;
  ivektor max_cell_dim;

  /* Apply shear */
#pragma omp parallel for
  for (k=0; k<ncells; ++k) {
    int i;
    cell *p;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) {
      p->ort Y(i) += shear_factor * p->ort X(i);
    }
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

void deform_sample(void) 
{
  int k;
  real box_x_half = 0.5 * box_x.x;

  /* loop over all atoms */
#pragma omp parallel for
  for (k=0; k<ncells; ++k) {
    int i;
    cell *p;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) {
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
    }
  }
}

#endif /* DEFORM */

