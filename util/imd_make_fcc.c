
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

#include <stdio.h>

int main()
{ 
  FILE *out;
  int x,y,z,typ,dimx,dimy,dimz,s,n;
  char fname[30];

  /* ask for parameters */
  printf("\nOutput file name: "); scanf("%s",fname );
  printf("\nDimension_x (even integer > 0): "); scanf("%d", &dimx );
  printf("\nDimension_y (even integer > 0): "); scanf("%d", &dimy );
  printf("\nDimension_z (even integer > 0): "); scanf("%d", &dimz );
  printf("\nBoundary strip with (integer >= 0): "); scanf("%d", &s );

  /* open output file */
  if ((out = fopen(fname,"w")) == NULL) return 2;
  n=0;
  for(z=s;z<dimz-s;z++)
    for(x=s;x<dimx-s;x++)
      for(y=s;y<dimy-s;y++) {
	typ = (x+y+z) % 2;
        if (typ==0) 
          fprintf(out,"%d %d 1.0 %f %f %f \n",n++,typ,x+.5,y+.5,z+.5);
      }
  
  fclose(out);
  
  return 0;
}


