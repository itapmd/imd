
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

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
#else
  fprintf(stderr,"Error: %s\n",msg);
#endif
  /* try to flush and close whatever we can */
  fflush(stderr);
  fflush(stdout);
  if ((myid==0) && ( eng_file!=NULL)) fclose( eng_file);
  if ((myid==0) && (msqd_file!=NULL)) fclose(msqd_file);
#ifdef MPI
  MPI_Abort(MPI_COMM_WORLD, 1);
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
        p->avpos X(i) += p->ort X(i) + p->sheet X(i);
        p->avpos Y(i) += p->ort Y(i) + p->sheet Y(i);
#ifndef TWOD
        p->avpos Z(i) += p->ort Z(i) + p->sheet Z(i);
#endif
        p->av_epot[i] += p->pot_eng[i];
    }
  }
}
  
#endif

/******************************************************************************
*
*  endian returns 1 if system is big endian, 0 if little endian
*
******************************************************************************/

int endian(void)
{
  unsigned short int word = 0x0001;
  unsigned char  *byte    = (unsigned char *) &word;
  return (byte[0] ? 0 : 1);
}

