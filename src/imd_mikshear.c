/******************************************************************************
*
* imd_mikshear.c -- New shear routines by Gunther Schaaf
*
* $RCSfile$
* $Revision$
* $Date$
*
******************************************************************************/

#include "imd.h"

#ifdef MIKSHEAR
/******************************************************************************
* apply_shear: probe is being sheared into x-direction by shear_delta         *
*              this is done if shear_flag == TRUE                             *
*              that is when the difference in pot eng is less than epsilon    *
******************************************************************************/

void apply_shear()
{
  cell *p;
  int i,j,k,l,m,tag;

#ifdef MPI
  if (1==parallel_output) {
    for (j = 1; j < cell_dim.x-1; ++j )
      for (k = 1; k < cell_dim.y-1; ++k )
#ifndef TWOD
	for (l = 1; l < cell_dim.z-1; ++l ) {
 	  p = PTR_3D_V(cell_array, j, k, l, cell_dim);
#else
 	p = PTR_2D_V(cell_array, j, k, cell_dim);
#endif
          for (i = 0;i < p->n; ++i)
            if (0 > p->nummer[i]) /* immoveable atoms: negative number */
              if (p->ort X(i) < strip) /* upper edge: forward */
                p->ort Y(i) -= shear_delta;
              else
                p->ort Y(i) += shear_delta;
 	};
    
  } else { 

    if (0==myid) {

      /* Manipulate data on CPU 0 */

      /* Manipulate own data */
      for (j = 1; j < cell_dim.x-1; ++j )
	for (k = 1; k < cell_dim.y-1; ++k )
#ifndef TWOD
	  for (l = 1; l < cell_dim.z-1; ++l ) {
	    p = PTR_3D_V(cell_array, j, k, l, cell_dim);
#else
 	p = PTR_2D_V(cell_array, j, k, cell_dim);
#endif
            for (i = 0;i < p->n; ++i)
            if (0 > p->nummer[i]) /* immoveable atoms: negative number */
              if (p->ort X(i) < strip) /* upper edge: forward */
                p->ort Y(i) -= shear_delta;
              else
                p->ort Y(i) += shear_delta;
	  };

      /* Receive data from other cpus and manipulate that */
      p   = PTR_3D_V(cell_array, 0, 0, 0, cell_dim);
      for ( m = 1; m < num_cpus; ++m)
	for (j = 1; j < cell_dim.x-1; ++j )
	  for (k = 1; k < cell_dim.y-1; ++k )
#ifndef TWOD
	    for (l = 1; l < cell_dim.z-1; ++l ) {
#else
 	p = PTR_2D_V(cell_array, j, k, cell_dim);
#endif
	      tag = PTR_3D_V(CELL_TAG, j, k, l, cell_dim);
	      recv_cell( p, m, tag );
              for (i = 0;i < p->n; ++i)
            if (0 > p->nummer[i]) /* immoveable atoms: negative number */
              if (p->ort X(i) < strip) /* upper edge: forward */
                p->ort Y(i) -= shear_delta;
              else
                p->ort Y(i) += shear_delta;
	    };

    } else { 
      /* Send data to cpu 0 */
      for (j = 1; j < cell_dim.x-1; ++j )
	for (k = 1; k < cell_dim.y-1; ++k )
#ifndef TWOD
	  for (l = 1; l < cell_dim.z-1; ++l ) {
#else
 	p = PTR_2D_V(cell_array, j, k, cell_dim);
#endif
	    p   = PTR_3D_V(cell_array, j, k, l, cell_dim);
	    tag = PTR_3D_V(CELL_TAG, j, k, l, cell_dim);
	    send_cell( p, 0, tag );
	  };
    };
  };
#else

#ifdef TWOD
  for (p = cell_array; /* loop over all atoms */
       p <= PTR_2D_V(cell_array,
		     cell_dim.x-1,
		     cell_dim.y-1,
		     cell_dim);
       ++p ) 
#else
  for (p = cell_array; /* loop over all atoms */
       p <= PTR_3D_V(cell_array,
		     cell_dim.x-1,
		     cell_dim.y-1,
		     cell_dim.z-1,
		     cell_dim);
       ++p ) 
#endif

     for (i = 0;i < p->n; ++i)
       if (0 > p->nummer[i]) /* immoveable atoms: negative number */
	 if (p->ort X(i) < strip) /* upper edge: forward */
	   p->ort Y(i) -= shear_delta;
	 else
	   p->ort Y(i) += shear_delta;

#endif
} /* apply_shear */


/******************************************************************************
*
* write_shear_energy writes selected properties to *sheng-file
*
******************************************************************************/

void write_shear_energy(int steps, int shear_steps)

{
  FILE *out;
  str255 fname;
  real  vol; 
  real part_kin_energy;
  real part_pot_energy; 

  /* Energiefile schreiben */

  sprintf(fname,"%s.sheng",outfilename);

  /* Groessen pro Tln ausgeben */
  part_pot_energy =                tot_pot_energy / natoms;
  part_kin_energy = ( 2.0 / 3.0 )* tot_kin_energy / natoms;
  vol = volume / natoms;

  out = fopen(fname,"a");
  if (NULL == out) error("Can't open shear-properties file.");

  fprintf(out,"%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
          shear_steps * shear_delta,
	  steps * timestep,
	  part_pot_energy,
	  part_kin_energy,
	  pressure,
	  vol);

  fclose(out);


} /* write_shear_energy */


/******************************************************************************
*                                                                             *
* do one shear step in main loop if ekin-value less than threshold            *
*                                                                             *
******************************************************************************/


void shear1step(int steps)
{
#ifdef MPI
  if (myid == 0)
#endif 
  write_shear_energy(steps, shear_steps);
  write_config(shear_steps);
  apply_shear();
  maxwell(temperature);
  shear_steps++;
}

#endif /* MIKSHEAR */
