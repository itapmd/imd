#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "sockutil.h"
#include "tokens.h"
#include "client.h"

/* definition of function for byte swapping (little endian <--> big endian) */
#define swap(x,y)  x ^= y;  y ^= x;  x ^= y
#define T_CONF 9
#define T_DIST 1
#define DIM 2
#define COLRES 245.0

extern unsigned short base_port;
char server_name[256]="visrn";
/* array pointers for received data */
static int socket_id;
extern int x_res,y_res;
extern int natoms;
extern float scalex,scaley,scalepot,scalekin,offspot,offskin;

/*-----------------------------------------------------------------*/
int receive_conf()
{
  int i;
  double tmp;
  int itmp;
  int anz;
  extern double *x, *y, *z, *vx, *vy, *vz;
  extern double *pot, *kin, *masse;
  extern int *nummer, *bcode;
  extern short int *sorte;
  float maxx,maxy,minx,miny,maxp,minp,maxk,mink;

  maxx=-1000;
  maxy=-1000;
  minx=1000;
  miny=1000;
  maxp=-1000;
  maxk=-1000;
  minp=1000;
  mink=1000;

  ReadFull(socket_id, (void *)&anz, sizeof(int));
  printf("%d\n", anz);
  nummer = (int *)calloc(anz, sizeof(int));
  sorte  = (short int *)calloc(anz, sizeof(int));
  masse  = (double *)calloc(anz, sizeof(double));
  x      = (double *)calloc(anz, sizeof(double));
  y      = (double *)calloc(anz, sizeof(double));
#ifdef TWOD
  z      = (double *)calloc(anz, sizeof(double));
#endif
  vx     = (double *)calloc(anz, sizeof(double));
  vy     = (double *)calloc(anz, sizeof(double));
#ifdef TWOD
  vz     = (double *)calloc(anz, sizeof(double));
#endif
  pot    = (double *)calloc(anz, sizeof(double));
  kin    = (double *)calloc(anz, sizeof(double));
  bcode  = (int    *)calloc(anz, sizeof(double));

  for (i=0;i<anz;i++) {
    kin[i]=0;
    ReadFull(socket_id, (void *)&itmp, sizeof(int));
    printf("%d %d ", i, itmp);
    sorte[i] = itmp;
    ReadFull(socket_id, (void *)&itmp, sizeof(int));
    nummer[i] = itmp;
    printf("%d\n",itmp);
    ReadFull(socket_id, (void *)&tmp, sizeof(double));
    masse[i] = tmp;
    ReadFull(socket_id, (void *)&tmp, sizeof(double));
    if (maxx<tmp) maxx=tmp;
    if (minx>tmp) minx=tmp;
    x[i] = tmp;
    ReadFull(socket_id, (void *)&tmp, sizeof(double));
    if (maxy<tmp) maxy=tmp;
    if (miny>tmp) miny=tmp;
    y[i] = tmp;
#ifndef TWOD
    ReadFull(socket_id, (void *)&tmp, sizeof(double));
    z[i] = tmp;
#endif
    ReadFull(socket_id, (void *)&tmp, sizeof(double));
    vx[i] = tmp;
    kin[i]+=tmp*tmp;
    ReadFull(socket_id, (void *)&tmp, sizeof(double));
    vy[i] = tmp;
    kin[i]+=tmp*tmp;
#ifndef TWOD
    ReadFull(socket_id, (void *)&tmp, sizeof(double));
    vz[i] = tmp;
    kin[i]+=tmp*tmp;
#endif
    ReadFull(socket_id, (void *)&tmp, sizeof(double));
    if (maxp<tmp) maxp=tmp;
    if (minp>tmp) minp=tmp;
    pot[i] = tmp;
    if (maxk<kin[i]) maxk=kin[i];
    if (mink>kin[i]) mink=kin[i];
  }

  scalex=2.0/(maxx-minx);
  scaley=2.0/(maxy-miny);
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

  for (i=0;i<x_res*y_res;i++) { 
    ReadFull(socket_id, (void *) &dummy, sizeof(float));
    potarray[i]=dummy;
  }
  for (i=0;i<x_res*y_res;i++) { 
    ReadFull(socket_id, (void *) &dummy, sizeof(float));
    kinarray[i]=dummy;
  }
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
  printf("hello: %f %f %f %f\n", maxk,mink,scalekin);

  return x_res*y_res;
}
/*-----------------------------------------------------------------*/

/*-----------------------------------------------------------------*/

void connect_server(char token)
/************************************************************************/
/* Connection to compute server and data handling	         	*/
/************************************************************************/
{
  /* check if server is the paragon, where little endian is used 
  if ( !strcmp(server_name,"paragon") || !strcmp(server_name,"paragon-fd") 
        || !strcmp(server_name,"paragon-fd") ){
    little_endian = 1;
  }
  else{
    printf("server is not the paragon --> little endian = 0\n");
    little_endian = 0;
  }
  */
  /* establish client-server connection */
  token = 10;
  if (initClient(&(socket_id),0,server_name)==-1){
    printf("Client-Server connection not established",
                 "check server name","");
  }
  else{  
    printf("using filedescriptor #%d as socket\n",socket_id);
    WriteFull(socket_id,(void *) &token, 1);
    switch (token) {
      case T_DISTRIBUTION:
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
}

/*-----------------------------------------------------------------*/
int initServer(){
/******************************************************************/
/* initializes server part for a socket-socket connection         */
/******************************************************************/
  printf("Start to establish connection ...\n");
  socket_id=OpenServerSocket(htons(base_port));
  if (socket_id < 1) return -1;
  return socket_id;
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
  if (initServer()==-1) {
    printf("Client-Server connection not established",
                 "check port adress","");
  }
  else {  
    printf("using filedescriptor #%d as socket\n",socket_id);
    WriteFull(socket_id,(void *) &token, 1);
    switch (token) {
      case T_CONF:
	natoms=receive_conf();
	printf("Got configuration with %d atoms.\n",natoms);
	break;
      case T_DISTRIBUTION:
	size=receive_dist();
	printf("Got distribution of size %d.\n",size);
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


int initClient(int *soc, int i, char *server_name){
/******************************************************************/
/* makes the socket-socket connection to server server_name       */
/******************************************************************/
  unsigned long var1;
  printf("Open client socket ...\n");
  var1 = GetIP(server_name);
  printf("GetIP: %u\n",var1);
  (*soc)=OpenClientSocket(GetIP(server_name),base_port+i);
  if (*soc < 1){
    return -1;
  }
  else{
    return 0;
  }
}

