#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <netdb.h>
#include <fcntl.h>
#include "string.h"
#include "sockutil.h"

#define min(a,b) ( ( (a) < (b) ) ? (a) : (b) )

/* ######################################
   ### Host-Name --> IP-Nummer        ###
   ### Ergebnis in Host Byte Order    ###
   ###################################### */

unsigned long GetIP(const char *name)
{
  struct hostent *host;
  struct in_addr inaddr;
  unsigned long IP=0;

  host = gethostbyname(name);
  if (host) {
    bcopy(host->h_addr, (char *)&inaddr, host->h_length);
    IP = inaddr.s_addr;
  }
  return IP;
}



/* ####################################################
   ### Schreibe <bytes> Byte aus <buffer> nach <fd> ###
   #################################################### */


int WriteFull(int fd, const void *buffer, int bytes)
{
  register char *bptr=(char *)buffer;
  register int written;
  register int nbytes=bytes;
  
  while (nbytes>0) {
    written = write(fd,(void *) bptr,nbytes);
    printf("Sent %d Bytes package\n",written);
    if (written < 0) return written; /* ERROR */
    nbytes-=written;
    bptr+=written;
    if (written==0) 
      fprintf(stderr,"Network write() returned 0: retrying...");
  }
  return 0;
}

/* #################################################
   ### Liest <bytes> Byte aus <fd> nach <buffer> ###
   ################################################# */

int ReadFull(int filedes, const void *buffer, int bytes)
{
  register char *bptr=(char *)buffer;
  register int nread;
  register int nbytes=bytes;
  
  while (nbytes>0) {
    nread = read(filedes,(void *) bptr,nbytes);
    printf("Received %d Byte package\n",nread);
    if (nread < 0) return nread; /* ERROR */
    nbytes-=nread;
    bptr+=nread;
    /*if (nread==0) break;*/
  }
  return bytes-nbytes;
}

/* #########################################################
   ### Open Socket for SERVER: Port in Network-ByteOrder ###
   ######################################################### */

int OpenServerSocket(u_short MyPort)
{
  struct sockaddr_in ServAddr,CliAddr;
  int soc,soc1,CliAddrLen;

  soc = socket(PF_INET,SOCK_STREAM,0);
  if (soc == -1) {
    perror("Socket creation failed:");
    return -1;
  }
  printf("server socket #%d\n",soc);
  
  bzero((char *) &ServAddr, sizeof(ServAddr));
  ServAddr.sin_family        = PF_INET;
  ServAddr.sin_addr.s_addr   = INADDR_ANY;
  ServAddr.sin_port          = MyPort;
  
  printf("calling bind \n");
  if ( bind(soc,(struct sockaddr *)&ServAddr,sizeof(ServAddr)) ) {
    perror("Bind failed");
    return -1;
  }
  
  printf("calling listening \n");
  if ( listen(soc,5) ) {
    perror("Listen failed");
    return -1;
  }
  
  bzero((char *) &CliAddr, sizeof(CliAddr));
  CliAddr.sin_family        = PF_INET;
  CliAddrLen=sizeof(CliAddr);
  
  printf("calling accept \n");
  soc1=accept(soc,(struct sockaddr *)  &CliAddr,&CliAddrLen);
  if (soc1 == -1) {
    perror("Accept failed");
    return -1;
  }
  printf("returning\n");
  close(soc);
  return soc1;
}


/* #####################################################################
   ### Open Socket for SERVER without accept: Port in Network-ByteOrder###
   ##################################################################### */

int OpenNBServerSocket(u_short MyPort)
{
  struct sockaddr_in ServAddr,CliAddr;
  int soc,soc1,CliAddrLen;

  soc = socket(AF_INET,SOCK_STREAM,0);
  if (soc == -1) {
    perror("Socket creation failed:");
    return -1;
  }
  printf("server socket #%d\n",soc);
  
  bzero((char *) &ServAddr, sizeof(ServAddr));
  ServAddr.sin_family        = AF_INET;
  ServAddr.sin_addr.s_addr   = INADDR_ANY;
  ServAddr.sin_port          = MyPort;
  
  if(fcntl(soc, F_SETFL, O_NDELAY)<0){
    perror("fcntl F_SETOWN");
    exit(1);
  }
  printf("calling bind \n");
  if ( bind(soc,(struct sockaddr *)&ServAddr,sizeof(ServAddr)) ) {
    perror("Bind failed");
    return -1;
  }
  
  printf("calling listening \n");
  if ( listen(soc,5) ) {
    perror("Listen failed");
    return -1;
  }
  if(fcntl(soc, F_SETOWN, getpid())<0){
    perror("fcntl F_SETOWN");
    exit(1);
  }
#ifdef NO_FASYNC
  if (fcntl(soc, F_SETFL, O_ASYNC) < 0) {
    perror("fcntl F_SETFL, O_ASYNC");
    exit(1);
  }
#else
  if (fcntl(soc, F_SETFL, FASYNC) < 0) {
    perror("fcntl F_SETFL, O_ASYNC");
    exit(1);
  }
#endif
  
  printf("returning \n");
  return soc;
}


/* ##########################################################
   ### Open Socket for CLIENT: Port, IP in Host-ByteOrder ###
   ########################################################## */

int OpenClientSocket(u_long toIP, u_short toPort)
{
  extern int errno;
  int soc, con_return;
  struct sockaddr_in ServAddr;

  soc = socket(PF_INET,SOCK_STREAM,0);
 /* printf("client socket #%d\n",soc); */
  if (soc == -1) {
    perror("Socket creation failed:");
    return -1;
  }
  
  bzero((char *) &ServAddr, sizeof(ServAddr));
  ServAddr.sin_family        = PF_INET;
  ServAddr.sin_port          = toPort;
  ServAddr.sin_addr.s_addr   = toIP;
  
  if ((con_return = connect(soc,(struct sockaddr *)&ServAddr,sizeof(ServAddr))) < 0) {
    /*perror("Connect failed");
    printf("Error bei connect()\n");
    printf("error return is %d\n",errno);*/
    shutdown(soc,2);
	close(soc);
	return -1;
  }
 printf("return socket #%d\n",soc);
 return soc;
}

/* #####################################
   ### Schreibe Sync-Byte nach <soc> ###
   ##################################### */


void WriteSync(int fd)
{
  char sync=0x2f;
  int res;
  while (1!=(res=write(fd,&sync,1)))
    if (res==-1) {
      perror("Write Sync:");
      return;
    }
}

/* ################################
   ### Lies Sync-Byte von <soc> ###
   ################################ */

void ReadSync(int fd)
{
  char sync=0x2f;
  int res;
  while (1!=(res=read(fd,&sync,1)))
    if (res==-1) {
      perror("Read Sync:");
      return;
    }
}

