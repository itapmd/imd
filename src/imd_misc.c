
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2006 Institute for Theoretical and Applied Physics,
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
#if defined (MPI) || defined(NEB)
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
  fprintf(stderr,"WARNING: %s \n",msg);
    fflush(stderr);
  }
}

/* warning message built from two strings */
void warning_str(char *msg, char *str)
{ 
  char buf[255];
  sprintf(buf, msg, str);
  warning(buf);
}

/* warning message built from three strings */
void warning_str_str(char *msg, char *str1, char *str2)
{ 
  char buf[255];
  sprintf(buf, msg, str1, str2);
  warning(buf);
}

/******************************************************************************
*
* error -- Complain and abort
*
******************************************************************************/

void imderror(char *msg)
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
#if defined (MPI) || defined(NEB)
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
#ifdef SOCKET_IO
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
* add_positions adds coordinates and potential energy for computation of 
* average position and potential energy
*
******************************************************************************/

void add_positions(void)
{
  int k;
  for (k=0; k<NCELLS; k++) {
    int i;
    cell* p;
    p = CELLPTR(k);
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
  avpos_cnt++;

}

#ifdef STRESS_TENS
/******************************************************************************
*
* add_stresstensors for computation of average stress tensor
*
******************************************************************************/

void add_presstensors(void)
{
  int k;
  for (k=0; k<NCELLS; k++) {
    int i;
    cell* p;
    p = CELLPTR(k);
    for (i=0; i<p->n; i++) {
      AVPRESSTENS(p,i,xx) += PRESSTENS(p,i,xx);
      AVPRESSTENS(p,i,yy) += PRESSTENS(p,i,yy);
      AVPRESSTENS(p,i,xy) += PRESSTENS(p,i,xy);
#ifndef TWOD
      AVPRESSTENS(p,i,zz) += PRESSTENS(p,i,zz);
      AVPRESSTENS(p,i,zx) += PRESSTENS(p,i,zx);
      AVPRESSTENS(p,i,yz) += PRESSTENS(p,i,yz);
#endif
    }
  }
}
#endif /*stress_tens */  
#endif /*AVPOS*/

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

integer SwappedInteger( integer i )
{
  union {
    integer i;
    unsigned char b[4];
  } dat1, dat2;
  dat1.i = i;
  dat2.b[0] = dat1.b[3];
  dat2.b[1] = dat1.b[2];
  dat2.b[2] = dat1.b[1];
  dat2.b[3] = dat1.b[0];
  return dat2.i;
}

float SwappedFloat( float f )
{
  union {
    float f;
    unsigned char b[4];
  } dat1, dat2;
  dat1.f = f;
  dat2.b[0] = dat1.b[3];
  dat2.b[1] = dat1.b[2];
  dat2.b[2] = dat1.b[1];
  dat2.b[3] = dat1.b[0];
  return dat2.f;
}

double SwappedDouble( double d )
{
  union {
    double d;
    unsigned char b[8];
  } dat1, dat2;
  dat1.d = d;
  dat2.b[0] = dat1.b[7];
  dat2.b[1] = dat1.b[6];
  dat2.b[2] = dat1.b[5];
  dat2.b[3] = dat1.b[4];
  dat2.b[4] = dat1.b[3];
  dat2.b[5] = dat1.b[2];
  dat2.b[6] = dat1.b[1];
  dat2.b[7] = dat1.b[0];
  return dat2.d;
}

#ifdef BGL

int get_free_mem(void)
{
  unsigned long tmp[1024], mem;
  tmp[0] = 12345;
  mem = (unsigned long) &tmp[1023] - (unsigned long) sbrk(0);
  return (int) (mem / 1048576);
}

#endif
