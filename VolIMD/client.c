/* program name ..............  VolIMD  
   file name .................  client.c
   author ....................  Roland Niemeier
   version ...................  
   date of last change .......  11/06/96
   description : interface for interactive visualization of
                 simulation on (parallel) computers. Derived
                 from VolVR
   purpose of client.c: connecting server and socket-socket exchange
*/
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
// for starting rendering functions
#include "ImgRenderAction.h"
#include "global.h"
extern GlobalInfoTyp global; /* viewer pointer, etc. */
unsigned  char *get_image_pointer(void);

/* declaration of local functions */
int initClient(int *soc, int i, char  *server_name);
void connect_server(char token);
float get_maximum(float *array, int size);
float get_minimum(float *array, int size);
void scale_array(float *array, int size, float min, float max, char *char_buffer);
void write_bytes_to_file(char filename[256],char *byte_array, int size);
extern multi_rotate_image(int axis, float angle, int number_of_rotations);
extern void new_classification(void);
extern void load_classified_volume(char *clVOLUME_FILE);
void automatic_load_and_classify_volume(char *filename, int xdim, int ydim, int zdim);
void automatic_load_image(char *filename, int xdim, int ydim);
void receive_picture();
void receive_distribution();
void show_message(const char *s1, const char *s2, const char *s3);

/* definition of function for byte swapping (little endian <--> big endian) */
#define swap(x,y)  x ^= y;  y ^= x;  x ^= y

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
extern float last_received_kinetic_minimum = 0; 
extern float last_received_kinetic_maximum = 0;
extern float last_received_potential_minimum = 0;
extern float last_received_potential_maximum = 0;
/* number of first port that will be connected */
unsigned short base_port = 31913;
/* array pointers for received data */
float *array;
char *byte_array;
static int socket_id;

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
  if (initClient(&(socket_id),0,server_name)==-1){
    show_message("Client-Server connection not established",
                 "check server name","");
  }
  else{  
    printf("using filedescriptor #%d as socket\n",socket_id);
    WriteFull(socket_id,(void *) &token, 1);
    switch (token) {
      case T_DISTRIBUTION:
        receive_distribution();
        break;
      case T_PICTURE:
        receive_picture();
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
  return 0;
}
/*-----------------------------------------------------------------*/

void connect_client(char token)
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
  if (initServer()==-1){
    show_message("Client-Server connection not established",
                 "check port adress","");
  }
  else{  
    printf("using filedescriptor #%d as socket\n",socket_id);
    WriteFull(socket_id,(void *) &token, 1);
    switch (token) {
      case T_DISTRIBUTION:
        receive_distribution();
        break;
      case T_PICTURE:
        receive_picture();
        break;
      default:
        fprintf(stderr,"Unimplemented token case!\n");
        break;
    }
  close(socket_id);
  }
}

/*-----------------------------------------------------------------*/

void receive_picture()     
/**************************************************************************/
/* Receives a picture of the kinetic energy distribution		  */
/**************************************************************************/
{
  unsigned char *this_image, image_resolution[4];
  int i=0,offset=0, size=0, x_res, y_res;
  global.viewer->getImgRenderAction()->getImgSize(&x_res,&y_res);
  this_image = get_image_pointer();
  /* length of integer may be different on different machine types ! */
  image_resolution[0] = x_res / 256; /* high byte */
  image_resolution[1] = x_res % 256; /* low byte */
  image_resolution[2] = y_res / 256;
  image_resolution[3] = y_res % 256;
  WriteFull(socket_id,(void *) &image_resolution, 4);
  size=3*x_res;  
  for (i = 0; i < y_res; i++){
    ReadFull(socket_id,(void *) &(this_image[offset]), size);
    offset += size;
  }
  printf("Calling DrawImage2D\n");
  global.viewer->getImgRenderAction()->DrawImage2D();
}

/*-----------------------------------------------------------------*/

void receive_distribution(){      
/************************************************************************/
/* Receives the kinetic and potential energy distribution	        */
/************************************************************************/
  int size=0, i, k, size_of_float, size_of_real;
  unsigned char distribution_resolution[6];
  int x_res, y_res, z_res;
  char load_filename[256] = "kin.dat";
  char *swap_buf;
  float min, max, f;
  /* assuming for instance that compute machine has same sizeof(float) */
  /* and that sizes of arrays are less than 65536 in each direction */
  /* length of integer may be different on different machine types ! */
  size_of_float = sizeof(float);
  size_of_real = size_of_float;

  /* FG - read test float, and check if endian swap is needed */
  ReadFull(socket_id, (void *) &f, sizeof(float));
  endianByteSwap = (f==(float)1.0) ? 0 : 1;
  if (endianByteSwap) {
    printf("endianByteSwap: yes\n");
  } else {
    printf("endianByteSwap: no\n");
  }

  if (TwoDSimulation){
    ReadFull(socket_id,(void *) distribution_resolution, 4);
    x_res = 256 * distribution_resolution[0] + distribution_resolution[1];
    y_res = 256 * distribution_resolution[2] + distribution_resolution[3];
    size = x_res * y_res * size_of_real;
  }
  else{
    ReadFull(socket_id,(void *) distribution_resolution, 6);
    x_res = 256 * distribution_resolution[0] + distribution_resolution[1];
    y_res = 256 * distribution_resolution[2] + distribution_resolution[3];
    z_res = 256 * distribution_resolution[4] + distribution_resolution[5];
    size = x_res * y_res * z_res * size_of_real;
  }
  printf("socket #%d: received size is %d\n",socket_id, size);
  if (size > 0){
    /* currently floats are assumed */
    array = (float *) malloc(size); 
    byte_array = (char *) malloc(size); 
    ReadFull(socket_id,(void *) &array[0], size);
    if (endianByteSwap){
      swap_buf = (char *)&(array[0]);
      for (k=0; k < size / size_of_real; k++){
        for (i=0; i < size_of_real/2; i++){
          swap(*(swap_buf+i), *(swap_buf+size_of_real-1-i));
        }
        swap_buf = swap_buf + size_of_real;
      } 
    }
    if (usePreviousLimits && 
        (last_received_potential_maximum != last_received_potential_minimum)
       ){
      min = last_received_potential_minimum;
      max = last_received_potential_maximum;
    }
    else{
      if (usePreviousLimits){
        printf("option use previous limits is not used (last minimum = last maximum)\n");
      }
      min = get_minimum(array,size/size_of_real);
      max = get_maximum(array,size/size_of_real);
      last_received_potential_minimum = min;
      last_received_potential_maximum = max;
    }
    scale_array(array, size/size_of_real, min, max, byte_array);
    write_bytes_to_file("pot.dat",byte_array,size/size_of_real);
    printf("socket #%d: received size is %d\n",socket_id, size);
    ReadFull(socket_id, (void *) &array[0], size);
    if (endianByteSwap){
      swap_buf = (char *)&(array[0]);
      for (k=0; k < size / size_of_real; k++){
        for (i=0; i < size_of_real/2; i++){
          swap(*(swap_buf+i), *(swap_buf+size_of_real-1-i));
        }
        swap_buf = swap_buf + size_of_real;
      } 
    }
    if (usePreviousLimits && 
        (last_received_potential_maximum != last_received_potential_minimum)
       ){
      min = last_received_potential_minimum;
      max = last_received_potential_maximum;
    }
    else{
      if (usePreviousLimits){
        printf("option use previous limits is not used (last minimum = last maximum)\n");
      }
      min = get_minimum(array,size/size_of_real);
      max = get_maximum(array,size/size_of_real);
      last_received_potential_minimum = min;
      last_received_potential_maximum = max;
    }
    scale_array(array, size/size_of_real, min, max, byte_array);
    write_bytes_to_file("kin.dat",byte_array,size/size_of_real);
    free(array); free(byte_array);
    /* loading either potential or kinetic energy volume file */
    if (!KineticEnergyFile) sprintf(load_filename,"pot.dat");
    if (TwoDSimulation){
      printf("x_res, y_res is %d %d\n",x_res,y_res);
      automatic_load_image(load_filename, x_res, y_res);   
    }
    else{
      automatic_load_and_classify_volume(load_filename, x_res, y_res, z_res);
    }
  }
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

float get_maximum(float *array, int size)
/******************************************************************/
/* returns the maximum of array[size]                             */
/******************************************************************/
{
  float *float_tmp, max;
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

float get_minimum(float *array, int size)
/******************************************************************/
/* returns the minimum of array[size]                             */
/******************************************************************/
{
  float *float_tmp, min;
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

void scale_array(float *array, int size, float min, float max, char *char_buffer)
/******************************************************************/
/* scale an array with values between min and max                 */
/******************************************************************/
{
  float grayvalue, width, *float_tmp;
  char *char_tmp;
  int i; 
  width = max - min;
  printf("scale_array: width of data is %f\n",width);
  float_tmp = array;
  char_tmp = char_buffer;
  if (width == 0){
    printf("cannot scale due to width 0");
  }
  else{
    for(i=0; i<size; i++){
      grayvalue =(255.0 * ((*array++)-min))/width;
      grayvalue = (grayvalue > 255.0) ? grayvalue = 255 : grayvalue;
      grayvalue = (grayvalue < 0.0) ? grayvalue = 0 : grayvalue;
      *char_buffer++ = (unsigned char) grayvalue;
    }
  }
  char_buffer = char_tmp;
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
