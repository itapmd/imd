
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

typedef double real;
typedef int integer;

#include "socket_io.h"

void error(char *msg)
{
  fprintf(stderr,"Error: %s\n",msg);
  exit(2);
}

void init_socket(unsigned short locPort, char *host, unsigned short toPort)
{
  if (host==NULL) {
    /* open server socket */
    fprintf(stderr, "Opening server socket on port %d ...", locPort);
    soc = OpenServerSocket(htons(locPort));
    if (soc < 1) error("Cannot open server socket");
    fprintf(stderr, " connected\n");
  }
  else {
    /* determine remote host */
    varIP = GetIP(host);
    if (varIP==0) {
      char msg[255];
      sprintf(msg, "Cannot find host %s", host);
      error(msg);
    }

    /* info message */
    if (locPort>0)
      fprintf(stderr, "Connecting from port %d to %s port %d ...",
              locPort, host, toPort);
    else
      fprintf(stderr, "Connecting to %s port %d ...", host, toPort);

    /* open client socket */
    soc=OpenClientSocket(varIP,htons(toPort),htons(locPort));
    if (soc <=0) error("Cannot open client socket");
    else fprintf(stderr," done\n");
  }
}

void vis_init()
{
  unsigned char b = VIS_INIT, res[4];
  WriteFull(soc,&b,1);
  ReadFull(soc,&res,4);
  printf("Protocol = %u.%u, endian = %u, dimension = %u\n",
         res[0], res[1], res[2], res[3] );
}

void vis_quit()
{
  unsigned char b = VIS_QUIT;
  WriteFull(soc,&b,1);
}

void vis_write_config_quit()
{
  unsigned char b = VIS_WRITE_QUIT;
  WriteFull(soc,&b,1);
}

void vis_restart_simulation()
{
  unsigned char b = VIS_RESTART;
  WriteFull(soc,&b,1);
}

void vis_init_atoms()
{
  unsigned char b = VIS_INIT_ATOMS;
  atoms_flag_t send_flags = {1,1,1,1,1,0}; 
  atoms_flag_t filt_flags = {0,0,0,0,0,0}; 
  atoms_filt_t filt_max   = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  atoms_filt_t filt_min   = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  WriteFull(soc,&b,1);

  ReadFull(soc,&send_flags,4*ATOMS_FLAG_SIZE);
  printf("Received send flags:\n");
  printf("  sorte=%d, ort=%d, impuls=%d, Ekin=%d, Epot=%d, nbanz=%d\n",
         send_flags.sorte, send_flags.ort, send_flags.impuls, 
         send_flags.Ekin, send_flags.Epot, send_flags.nbanz );

  ReadFull(soc,&filt_min,  4*ATOMS_FILT_SIZE);
  printf("Received filter min:\n");
  printf("  sorte=%e, x=%e, y=%e, z=%e, Ekin=%e, Epot=%e, nbanz=%e\n",
         filt_min.sorte, filt_min.x, filt_min.y, filt_min.z,
         filt_min.Ekin, filt_min.Epot, filt_min.nbanz );

  ReadFull(soc,&filt_max,  4*ATOMS_FILT_SIZE);
  printf("Received filter max:\n");
  printf("  sorte=%e, x=%e, y=%e, z=%e, Ekin=%e, Epot=%e, nbanz=%e\n",
         filt_max.sorte, filt_max.x, filt_max.y, filt_max.z,
         filt_max.Ekin, filt_max.Epot, filt_max.nbanz );
}

void vis_write_atoms()
{
  unsigned char b = VIS_WRITE_ATOMS;
  atoms_flag_t send_flags = {1,1,1,1,1,0}; 
  atoms_flag_t filt_flags = {0,0,0,0,0,0}; 
  atoms_filt_t filt_max   = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  atoms_filt_t filt_min   = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  float *a=NULL;
  integer step, atlen, num=0;
  int i, j;

  WriteFull(soc,&b,1);
  WriteFull(soc,&send_flags,4*ATOMS_FLAG_SIZE);
  WriteFull(soc,&filt_flags,4*ATOMS_FLAG_SIZE);
  WriteFull(soc,&filt_max,  4*ATOMS_FILT_SIZE);
  WriteFull(soc,&filt_min,  4*ATOMS_FILT_SIZE);

  ReadFull(soc,&step,4);
  ReadFull(soc,&atlen,4);
  printf("Step %d, atlen = %d\n", step, atlen);
  a = (float *) realloc(a,4*atlen);

  do {
    ReadFull(soc,&num,4);
    printf("Receiving block of %d atoms... ", num);
    for (i=0; i<num; i++) {
      ReadFull(soc,a,4*atlen);
#ifdef DEBUG
      for (j=0; j<atlen; j++) printf("%e ",a[j]);
      printf("\n");
#endif
    }
    printf("done\n");
  } while (num); 
  printf("Receiving atoms finished\n");
}

void vis_change_params()
{
  unsigned char b = VIS_CHANGE_PARAMS;
  integer pargrp=VIS_PARAM_DEFORM, flag0=0, flag1=1, step;
  float   tmp, newparam=1.1;

  WriteFull(soc,&b,1);
  WriteFull(soc,&pargrp,4);
  WriteFull(soc,&flag0,4);
  ReadFull(soc,&step,4);
  ReadFull(soc,&tmp,4);
  printf("step = %d, old deform_size = %f\n", step, tmp);

  WriteFull(soc,&b,1);
  WriteFull(soc,&pargrp,4);
  WriteFull(soc,&flag1,4);
  WriteFull(soc,&newparam,4);
  ReadFull(soc,&step,4);
  ReadFull(soc,&tmp,4);
  printf("step = %d, new deform_size = %f\n", step, tmp);
}

int main(int argc, char **argv)
{
  int  n;
  char *host=NULL;
  unsigned short locPort=0, toPort=0;

  if (argc==2) {      /* server socket */
    locPort = (unsigned short) atoi(argv[1]);
  }
  else if (argc==4) { /* client socket */
    locPort = (unsigned short) atoi(argv[1]);
    host    = strdup(argv[2]);
    toPort  = (unsigned short) atoi(argv[3]);
  }
  else {
    printf("Usage: socktest <port>                    # open server socket\n");
    printf("       socktest <locPort> <host> <toPort> # open client socket\n");
    return 0;
  }

  init_socket(locPort, host, toPort);

  do {

    scanf("%d", &n);
    printf("Read flag %d\n", n);

    switch (n) {
      case VIS_INIT:
        vis_init();
        break;
      case VIS_QUIT:
        vis_quit();
	break;
      case VIS_WRITE_QUIT:
        vis_write_config_quit();
	break;
      case VIS_INIT_ATOMS:
        vis_init_atoms();
        break;
      case VIS_WRITE_ATOMS:
        vis_write_atoms();
        break;
      case VIS_CHANGE_PARAMS:
        vis_change_params();
        break;
      case VIS_RESTART:
        vis_restart_simulation();
        break;
      default:
        if (n>0) fprintf(stderr, "Unknown command %d\n", n);
        break;
    }
  } while (n>0);

  shutdown(soc,2);
  close(soc);
  return 0;

}
