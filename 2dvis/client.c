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

/* definition of variables for local (routines in this file) and extern usage */
/* boolean for determining if kinetic or potential energy data 
   should be displayed */
extern unsigned char KineticEnergyFile=0;
/* boolean for 2D Simulation */
extern unsigned char TwoDSimulation;
/* default for compute server */
extern char server_name[256] = "visrn";
extern int endianByteSwap=0;
/* for hints in the histogram display the maxima and minima of a received
   data set is stored in last_received... */
extern unsigned char usePreviousLimits;
extern double last_received_kinetic_minimum = 0; 
extern double last_received_kinetic_maximum = 0;
extern double last_received_potential_minimum = 0;
extern double last_received_potential_maximum = 0;
/* number of first port that will be connected */
/*unsigned short base_port = 12955;*/
extern unsigned short base_port;
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
  extern double *pot, *masse;
  extern int *nummer, *sorte;
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
  nummer = (int *)calloc(anz, sizeof(int));
  sorte = (int *)calloc(anz, sizeof(int));
  masse = (double *)calloc(anz, sizeof(double));
  x = (double *)calloc(anz, sizeof(double));
  y = (double *)calloc(anz, sizeof(double));
  z = (double *)calloc(anz, sizeof(double));
  vx = (double *)calloc(anz, sizeof(double));
  vy = (double *)calloc(anz, sizeof(double));
  vz = (double *)calloc(anz, sizeof(double));
  pot = (double *)calloc(anz, sizeof(double));

  /*  ReadFull(socket_id, (void *)&boxx, sizeof(double));
  ReadFull(socket_id, (void *)&boxy, sizeof(double));*/
  for (i=0;i<anz;i++) {
    ReadFull(socket_id, (void *)&itmp, sizeof(int));
    nummer[i] = itmp;
    ReadFull(socket_id, (void *)&itmp, sizeof(int));
    sorte[i] = itmp;
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
    ReadFull(socket_id, (void *)&tmp, sizeof(double));
    vy[i] = tmp;
#ifndef TWOD
    ReadFull(socket_id, (void *)&tmp, sizeof(double));
    vz[i] = tmp;
#endif
    ReadFull(socket_id, (void *)&tmp, sizeof(double));
    if (maxp<tmp) maxp=tmp;
    if (minp>tmp) minp=tmp;
    pot[i] = tmp;
  }

  scalex=2.0/(maxx-minx);
  scaley=2.0/(maxy-miny);
  scalepot=COLRES/(maxp-minp);
  offspot=minp;

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

double get_maximum(double *array, int size)
/******************************************************************/
/* returns the maximum of array[size]                             */
/******************************************************************/
{
  double *float_tmp, max;
  int i;
  float_tmp = array;
  max = *array++;
  for(i=1; i<size; i++){
    max = (max > *array ? max : *array); array++;
  }
  array = float_tmp;
  printf("maximum is %f\n",max);
  return max;
}

double get_minimum(double *array, int size)
/******************************************************************/
/* returns the minimum of array[size]                             */
/******************************************************************/
{
  double *float_tmp, min;
  int i;
  float_tmp = array;
  min = *array++;
  for(i=1; i<size; i++){
    min = (min < *array ? min : *array); array++;
  }
  array = float_tmp;
  printf("minimum is %f\n",min);
  return min;
}

void scale_array(double *array, int size, double min, double max, double *buffer)
/******************************************************************/
/* scale an array with values between min and max                 */
/******************************************************************/
{
  double grayvalue, width, *float_tmp;
  int i; 
  width = max - min;
  printf("scale_array: width of data is %f\n",width);
  float_tmp = array;
  if (width == 0){
    printf("cannot scale due to width 0");
  }
  else{
    for(i=0; i<size; i++){
      grayvalue =((*array++)-min)/width;
      grayvalue = (grayvalue > 1.0) ? grayvalue = 1.0 : grayvalue;
      grayvalue = (grayvalue < 0.0) ? grayvalue = 0.0 : grayvalue;
      *buffer++ = grayvalue;
    }
  }
  array = float_tmp;
}

void write_bytes_to_file(char filename[256],char *byte_array, int size)
/******************************************************************/
/* write size bytes to file filename                              */
/******************************************************************/
{
  int fd;
  if( ( fd=creat(filename,0640)) == -1 ) {
    printf("Cannot create file named %s for output\n",filename);
  }
  if ((fd = open(filename, O_CREAT | O_RDWR)) < 0) {
    perror("open");
    fprintf(stderr, "could not open %s\n", filename);
  }
  else{
    printf("writing data to file ...");
    write(fd,byte_array,size);
    close(fd);
    printf("finished\n");
  }
}
