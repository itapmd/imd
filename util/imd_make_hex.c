
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

#include <math.h>
#include <stdio.h>

int main()
{ 
  FILE *out;
  int imax,jmax, typ, natoms, i, j;
  float x, y, dimx, dimy, strip;
  char fname[30];

  /* ask for parameters */
  printf("\nOutput file name: "); scanf("%s",fname );
  printf("\nDimension_x (n*sqrt(3), n integer > 0): "); scanf("%d", &imax );
  printf("\nDimension_y (integer > 0): "); scanf("%d", &jmax );
  printf("\nBoundary strip width: "); scanf("%f", &strip );

  /* open output file */
  if ((out = fopen(fname,"w")) == NULL) return 2;

  natoms = 0;
  dimx = (float) imax * sqrt(3.0);
  dimy = (float) jmax;

  for (i=0; i<2*imax; i++)
    for (j=0; j<2*jmax; j++) {

      typ = (i+j) % 2;
      if (typ>0) continue;

      x = (i+0.5) * sqrt(3.0) * 0.5;
      y = (j+0.5) * 0.5;

      if ((x<strip) || (x>dimx-strip) || (y<strip) || (y>dimy-strip))
        continue;

      fprintf(out,"%d %d 1.0 %f %f\n", natoms++, typ, x, y);
    }  

  fclose(out);
  return 0;

}


