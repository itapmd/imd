
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
#ifdef USE_SOCKETS
  close_socket();
#endif
  exit(2);
}

/* error message built from two strings */
void error_str(char *msg, char *str)
{
  char buf[255];
  sprintf(buf, msg, str);
  error(buf);
}

/* error message built from three strings */
void error_str_str(char *msg, char *str1, char *str2)
{
  char buf[255];
  sprintf(buf, msg, str1, str2);
  error(buf);
}

#ifdef AVPOS

/******************************************************************************
*
* add_positions adds coordinates and potential energy for computation of average
* 
* position and potential energy
*
******************************************************************************/

void add_positions(void)
{
  int k;
  for (k=0; k<ncells; k++) {
    int i;
    cell* p;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; i++) {
        AV_POS(p,i,X) += ORT(p,i,X) + SHEET(p,i,X);
        AV_POS(p,i,Y) += ORT(p,i,Y) + SHEET(p,i,Y);
#ifndef TWOD
        AV_POS(p,i,Z) += ORT(p,i,Z) + SHEET(p,i,Z);
#endif
        AV_EPOT(p,i)  += POTENG(p,i);
    }
  }

#ifdef NPT
  av_box_x.x += box_x.x;
  av_box_x.y += box_x.y;
  av_box_y.x += box_y.x;
  av_box_y.y += box_y.y;
#ifndef TWOD
  av_box_x.z += box_x.z;  
  av_box_y.z += box_y.z;
  av_box_z.x += box_z.x;
  av_box_z.y += box_z.y;  
  av_box_z.z += box_z.z;  
#endif
#endif

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

