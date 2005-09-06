
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2005 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#ifndef _SOCKUTIL_H_
#define _SOCKUTIL_H_
#include <sys/types.h>
int WriteFull(int filedes, const void *buffer, int nbytes);
int ReadFull(int filedes, const void *buffer, int nbytes);
unsigned long GetIP(const char *name);
int OpenServerSocket(u_short MyPort);
int OpenClientSocket(u_long toIP, u_short toPort, u_short locPort);
void WriteSync(int fd);
void ReadSync(int fd);
#endif
#ifdef CRAY
#define NO_FASYNC
#endif

