#include <stdio.h>
#include <strings.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <fcntl.h>
/*#include "bstring.h"*/
#include "sockutil.h"

#define min(a,b) ( ( (a) < (b) ) ? (a) : (b) )

/* ######################################
   ### Host-Name --> IP-Nummer        ###
   ### Ergebnis in Network Byte Order ###
   ###################################### */

unsigned long GetIP(const char *name)
{
  struct hostent *hostinfo;
  hostinfo=gethostbyname(name);
  if (hostinfo)
    return *((unsigned long *)*hostinfo->h_addr_list);
  else
    return 0;
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
#ifdef DEBUG
    printf("Sent %d Bytes package\n",written);
#endif
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
#ifdef DEBUG
    printf("nbytes=%d, Received %d Byte package\n",nbytes,nread);*/
#endif
    if (nread < 0) return nread; /* ERROR */
    nbytes-=nread;
    bptr+=nread;
    if (nread==0) break;
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
#ifdef DEBUG
  printf("server socket #%d\n",soc);
#endif

  bzero((char *) &ServAddr, sizeof(ServAddr));
  ServAddr.sin_family        = PF_INET;
  ServAddr.sin_addr.s_addr   = INADDR_ANY;
  ServAddr.sin_port          = MyPort;
  
#ifdef DEBUG
  printf("calling bind \n");
#endif
  if ( bind(soc,&ServAddr,sizeof(ServAddr)) ) {
    perror("Bind failed");
    return -1;
  }
  
#ifdef DEBUG
  printf("calling listening \n");
#endif
  if ( listen(soc,5) ) {
    perror("Listen failed");
    return -1;
  }
  
  bzero((char *) &CliAddr, sizeof(CliAddr));
  CliAddr.sin_family        = PF_INET;
  CliAddrLen=sizeof(CliAddr);
  
#ifdef DEBUG
  printf("calling accept \n");
#endif
  soc1=accept(soc,&CliAddr,&CliAddrLen);
  if (soc1 == -1) {
    perror("Accept failed");
    return -1;
  }
#ifdef DEBUG
  printf("returning\n");
#endif
  close(soc);
  return soc1;
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
#ifdef DEBUG
  printf("client socket #%d\n",soc);
#endif
  if (soc == -1) {
    perror("Socket creation failed:");
    return -1;
  }
  
  bzero((char *) &ServAddr, sizeof(ServAddr));
  ServAddr.sin_family        = PF_INET;
  ServAddr.sin_port          = toPort;
  ServAddr.sin_addr.s_addr   = toIP;
  
  if ((con_return = connect(soc,&ServAddr,sizeof(ServAddr))) < 0) {
    perror("Connect failed");
    printf("Error bei connect()\n");
    printf("error return is %d\n",errno);
    return -1;
  }
#ifdef DEBUG
 printf("return socket #%d\n",soc);
#endif
 return soc;
}











