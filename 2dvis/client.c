#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include "sockutil.h"
#include "tokens.h"
#include "client.h"
#include "globals.h"
#include "prototypes.h"
#include "makros.h"


void iendian(int *number);
void sendian(short int *number);
void dendian(double *number);

char server_name[256]="visrn";
/* array pointers for received data */
static int socket_id=0;

/*-----------------------------------------------------------------*/
int receive_conf()
{
  int i,j,k,size,itmp;
  short int stmp;
  char ctmp;
  double tmp;
  int anz;
  float f;

  maxx=-1000;
  maxy=-1000;
#ifndef TWOD
  maxz=-1000;
#endif
  minx=1000;
  miny=1000;
#ifndef TWOD
  minz=1000;
#endif
  maxp=-1000;
  maxk=-1000;
  minp=1000;
  mink=1000;

  ReadFull(socket_id, (void *)&f, sizeof(float));
  if (f == 1)
    endian_byte_swap = 0;
  else
    endian_byte_swap = 1;
  endian_byte_swap = 0;

  ReadFull(socket_id, (void *)&anz, sizeof(int));

  if (endian_byte_swap) {
    iendian(&anz);
  }

  nummer = (int *)calloc(anz, sizeof(int));
  sorte  = (short int *)calloc(anz, sizeof(int));
  masse  = (double *)calloc(anz, sizeof(double));
  x      = (double *)calloc(anz, sizeof(double));
  y      = (double *)calloc(anz, sizeof(double));
#ifndef TWOD
  z      = (double *)calloc(anz, sizeof(double));
#endif
  vx     = (double *)calloc(anz, sizeof(double));
  vy     = (double *)calloc(anz, sizeof(double));
#ifndef TWOD
  vz     = (double *)calloc(anz, sizeof(double));
#endif
  pot    = (double *)calloc(anz, sizeof(double));
  kin    = (double *)calloc(anz, sizeof(double));
  bcode  = (int    *)calloc(anz, sizeof(double));

  size=anz*sizeof(int);
  ReadFull(socket_id,(void *) &nummer[0], size);
  size=anz*sizeof(short int);
  ReadFull(socket_id,(void *) &sorte[0], size);
  size=anz*sizeof(double);
  ReadFull(socket_id,(void *) &masse[0], size);
  ReadFull(socket_id,(void *) &x[0], size);
  ReadFull(socket_id,(void *) &y[0], size);
#ifndef TWOD
  ReadFull(socket_id,(void *) &z[0], size);
#endif
  ReadFull(socket_id,(void *) &vx[0], size);
  ReadFull(socket_id,(void *) &vy[0], size);
#ifndef TWOD
  ReadFull(socket_id,(void *) &vz[0], size);
#endif
  ReadFull(socket_id,(void *) &pot[0], size);

  for (i=0;i<anz;i++) {
    if (endian_byte_swap) {
      iendian(&nummer[i]);
      sendian(&sorte[i]);
      dendian(&masse[i]);
      dendian(&x[i]);
      dendian(&y[i]);
#ifndef TWOD
      dendian(&z[i]);
#endif
      dendian(&vx[i]);
      dendian(&vy[i]);
#ifndef TWOD
      dendian(&vz[i]);
#endif
      dendian(&pot[i]);
    }
    kin[i] = vx[i]*vx[i]+vy[i]*vy[i];
    if (maxx<x[i]) maxx=x[i];
    if (minx>x[i]) minx=x[i];
    if (maxy<y[i]) maxy=y[i];
    if (miny>y[i]) miny=y[i];
#ifndef TWOD
    if (maxz<z[i]) maxz=z[i];
    if (minz>z[i]) minz=z[i];
#endif
    if (maxp<pot[i]) maxp=pot[i];
    if (minp>pot[i]) minp=pot[i];
    if (maxk<kin[i]) maxk=kin[i];
    if (mink>kin[i]) mink=kin[i];
  }

  scalex=2.0/(maxx-minx);
  scaley=2.0/(maxy-miny);
#ifndef TWOD
  scalez=2.0/(maxz-minz);
#endif
  if (maxp==minp)
    scalepot=1.0;
  else
    scalepot=COLRES/(maxp-minp);
  if (maxk==mink)
    scalekin=1.0;
  else
    scalekin=COLRES/(maxk-mink);
  offspot=minp;
  offskin=mink;

  return anz;
}

/*-----------------------------------------------------------------*/
int receive_dist()
{
  int size, i;
  unsigned char distribution_resolution[6];
  float maxx,maxy,minx,miny,maxp,minp,maxk,mink,dummy;
  float patarray[400];
  extern float *potarray;
  extern float *kinarray;

  maxx=-1000;
  maxy=-1000;
  minx=1000;
  miny=1000;
  maxp=-1000;
  maxk=-1000;
  minp=1000;
  mink=1000;
  size=0;

  ReadFull(socket_id,(void *) &dummy, sizeof(float));
  /* assuming for instance that compute machine has same sizeof(float) */
  /* and that sizes of arrays are less than 65536 in each direction */
  /* length of integer may be different on different machine types ! */
  ReadFull(socket_id,(void *) distribution_resolution, 4);
  x_res = 256 * distribution_resolution[0] + distribution_resolution[1];
  y_res = 256 * distribution_resolution[2] + distribution_resolution[3];
  size = x_res * y_res * sizeof(float);
  printf("%d\n", size);
  /*
  if (size > 0) {
    if (potarray==NULL) potarray = (float *) malloc(size);
    if (kinarray==NULL) kinarray = (float *) malloc(size); 
    ReadFull(socket_id, (void *) &potarray[0], size);
    printf("socket #%d: received size is %d\n",socket_id, size);
    ReadFull(socket_id, (void *) &kinarray[0], size);
    printf("socket #%d: received size is %d\n",socket_id, size);
  }
  */
  potarray = (float *) calloc(x_res*y_res,sizeof(float));
  kinarray = (float *) calloc(x_res*y_res,sizeof(float));

  if (potarray==NULL) {printf("No memory for pot\n");return 0;}
  if (kinarray==NULL) {printf("No memory for kin\n");return 0;}

  ReadFull(socket_id,(void *) &potarray[0], size);
  ReadFull(socket_id,(void *) &kinarray[0], size);
  /*
  for (i=0;i<x_res*y_res;i++) { 
    ReadFull(socket_id, (void *) &dummy, sizeof(float));
    potarray[i]=dummy;
  }
  for (i=0;i<x_res*y_res;i++) { 
    ReadFull(socket_id, (void *) &dummy, sizeof(float));
    kinarray[i]=dummy;
  }
  */
  for (i=0;i<x_res*y_res;i++) {
    if (maxp<potarray[i]) maxp=potarray[i];
    if (minp>potarray[i]) minp=potarray[i];
    if (maxk<kinarray[i]) maxk=kinarray[i];
    if (mink>kinarray[i]) mink=kinarray[i];
  }
  scalex=2.0/(float)x_res;
  scaley=2.0/(float)y_res;
  offspot=minp;
  offskin=mink;
  scalepot=COLRES*1.0/(maxp-minp);
  scalekin=COLRES*1.0/(maxk-mink);

  return x_res*y_res;
}

/*-----------------------------------------------------------------*/
int initServer(int base_port){
/******************************************************************/
/* initializes server part for a socket-socket connection         */
/******************************************************************/
#ifdef DEBUG
  printf("Start to establish connection ...\n");
#endif
  socket_id=OpenServerSocket(htons(base_port));
  if (socket_id < 1) return -1;
  return socket_id;
}

/*****************************************************************************
*
* init_client, if the visualization acts as a client
*
*****************************************************************************/

void init_client() {

  fprintf(stderr, "baseport is %d\n", base_port_int);
  base_port_int = htons(base_port_int); /* we need it in network byte order */
  varIP = GetIP(sim_host);
  if (varIP==0) {
    error("gethostbyname() failed, check sim_host or specify IP number\n");
  } 
  else {
    printf("display_host is %s\n", sim_host);
  }
}


/*-----------------------------------------------------------------*/

int connect_client(char token)
/************************************************************************/
/* Connection to compute server and data handling	         	*/
/************************************************************************/
{
  int size;

  natoms=0;
  size=0;
  /* establish client-server connection */
  if (initServer(base_port)==-1) {
    printf("Client-Server connection not established",
	   "check port adress %d\n",base_port);
  }
  else {
#ifdef DEBUG
    printf("using filedescriptor #%d as socket\n",socket_id);
#endif
    WriteFull(socket_id,(void *) &token, 1);
    WriteFull(socket_id,(void *) &temperature, 4);
    switch (token) {
      case T_CONF:
	natoms=receive_conf();
	break;
      case T_DISTRIBUTION:
	size=receive_dist();
        break;
      case T_PICTURE:
        break;
      case T_QUIT:
        break;
      default:
        fprintf(stderr,"Unimplemented token case!\n");
        break;
    }
    close(socket_id);
  }
  return size+natoms;
}

void iendian(int *anz) {
  int i,s;
  char *p,h;

  s=sizeof(int);
  p = (char *)calloc(s,sizeof(char));
  p = (char *)&anz[0];

  for (i=0;i<s/2;i++) {
    h = p[s-1-i];
    p[s-1-i] = p[i];
    p[i] = h; 
  }
}

void sendian(short int *anz) {
  int i,s;
  char *p,h;

  s=sizeof(short int);
  p = (char *)calloc(s,sizeof(char));
  p = (char *)&anz[0];

  for (i=0;i<s/2;i++) {
    h = p[s-1-i];
    p[s-1-i] = p[i];
    p[i] = h; 
  }
}

void dendian(double *anz) {
  int i,s;
  char *p,h;

  s=sizeof(double);
  p = (char *)calloc(s,sizeof(char));
  p = (char *)&anz[0];

  for (i=0;i<s/2;i++) {
    h = p[s-1-i];
    p[s-1-i] = p[i];
    p[i] = h; 
  }
}

