
/******************************************************************************
*
* imd_misc.c -- Some Misc. Routines for the imd package
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/******************************************************************************
*
*  Usage -- educate users
*
******************************************************************************/

void usage(void)
{ 
  if (myid==0) {
    fprintf(stderr,"%s [-r<nnn>] [-p paramter-file]\n",progname); 
    fflush(stderr);
  }
#ifdef MPI
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  exit(1); 
}

/******************************************************************************
*
*  print a warning
*
******************************************************************************/

void warning(char *msg)
{ 
  if (myid==0) {
    fprintf(stderr,"WARNING: %s\n",msg);
    fflush(stderr);
  }
}

/******************************************************************************
*
* error -- Complain and abort
*
******************************************************************************/

void error(char *msg)
{
#ifdef MPI
  fprintf(stderr,"Error on CPU %d: %s\n",myid,msg);
  fflush(stderr);
  MPI_Abort(MPI_COMM_WORLD, 1);
#else
  fprintf(stderr,"Error: %s\n",msg);
  fflush(stderr);
#endif
  exit(2);
}

#ifdef AVPOS

/******************************************************************************
*
* add_position adds coordinates and potential energy for computation of average
* 
* position and potential energy
*
******************************************************************************/

void add_position(void)
{
  int k;
  for (k=0; k<ncells; k++) {
    int i;
    cell* p;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; i++) {
        p->ort_ref X(i) += p->ort X(i) + p->sheet X(i) * box_x.x;
        p->ort_ref Y(i) += p->ort Y(i) + p->sheet Y(i) * box_y.y;
#ifndef TWOD
        p->ort_ref Z(i) += p->ort Z(i) + p->sheet Z(i) * box_z.z;
#endif
        p->Epot_ref[i]  += p->pot_eng[i];
    }
  }
}
  
#endif
