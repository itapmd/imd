
#ifndef USE_SOCKETS
#define USE_SOCKETS
#endif

#include "imd.h" 

/*****************************************************************************
*
*  write float array to socket (unused - not endianness proof!)
*
*****************************************************************************/

void write_to_socket(float *dist, int size)
{
  printf("write_to_socket size: %d\n",size);
  WriteFull(soc, (void *)&size, sizeof(int));
  if (size > 0) WriteFull(soc, (void *)dist, size );
}


/*****************************************************************************
*
*  write quit to socket (unused - needed only if IMD is server)
*
*****************************************************************************/

void write_quit_to_socket()
{
  int quit = T_QUIT;
#ifdef MPI
  if (0==myid)
#endif
  {
    printf("write_to_socket quit: %d\n",quit);
    WriteFull(soc, (void *)&quit, sizeof(int));
  }
}


void init_server(){
/***********************************************************/
/* initializes nonblocking socket and SIGIO handler        */
/* initializes socket_flag and baseport  		   */
/***********************************************************/
  /* initialize socket(s) */
    socket_flag = 0; /* flag for the main_loop */
    soc = 0; /* soc is the filedescriptor if a connection is established */

    /* instead of using a socket for every process like */
    /* socket_id = OpenNBServerSocket(baseport+myid);*/
    socket_id = OpenNBServerSocket(baseport);

    /* setup signal_handler check_io_signal */
    /* better use sigaction */
    /* act.sa_handler = (void *)check_io_signal_using_select; */
    act.sa_handler = check_io_signal_using_select; 
    sigemptyset(&act.sa_mask);
    act.sa_flags = 0;
    act.sa_flags |= SA_RESTART;
    if (sigaction(SIGIO, &act, &oact) != 0){
      printf("sigaction for SIGIO failed\n");
    }
    else{
      printf("sigaction installed handler for SIGIO\n");
    }
    /* than the older signal function 
    signal(SIGIO,(void *)check_io_signal); */
}

/* ---------------------------------------------------------------- */

void check_io_signal_using_poll(){
/****************************************************************/
/* checks the SIGIO signal if the request is from  the         	*/
/* socket. If request is from the socket the socket-socket     	*/
/* connection is established. Then reads socket-token from     	*/ 
/* socket and sets socket_flag         				*/
/****************************************************************/
/*  int i=0,  quit = T_QUIT;
  char token=0;
  struct sockaddr_in from;
  int len = sizeof (from); 
  struct pollfd socket_info;
  struct sigaction act, oact; 

  printf("enter check_io_signal_using_poll\n");


  act.sa_handler = SIG_IGN;
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;
  act.sa_flags |= SA_RESTART;
  sigaction(SIGIO,&act, &oact);

  socket_info.fd = socket_id;

#ifdef MPI
  socket_info.events = POLLNORM;
#else
  socket_info.events = POLLRDNORM;
#endif
  poll(&socket_info,1,0);
#ifdef MPI
  printf("node %d: socket_info.revents is %d\n",myid,socket_info.revents);
#endif
#ifdef MPI
  if (socket_info.revents == 1){
#else
  if (socket_info.revents == 64){
#endif
    soc = accept(socket_id,(struct sockaddr *)&from,&len);
#ifdef MPI
    printf("Node %d: return value of accept is %d\n",myid, soc);
#endif
    if(soc <= 0){
      printf("Connection not established, return code of accept is %d\n",soc);
    }
    else{
      ReadFull(soc,&socket_flag,1);
    }
  }
  act.sa_handler = (void *)check_io_signal_using_poll;
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;
  act.sa_flags |= SA_RESTART;
  sigaction(SIGIO,&act, &oact); */
  /* signal(SIGIO,(void *)check_io_signal_using_poll);  */
  printf("Returning from check_io_signal_using_poll\n");
}

/* ---------------------------------------------------------------- */

void check_io_signal_using_select(){
/****************************************************************/
/* checks the SIGIO signal if the request is from  the         	*/
/* socket. If request is from the socket the socket-socket     	*/
/* connection is established. Then reads socket-token from     	*/ 
/* socket and sets socket_flag         				*/
/****************************************************************/
  int i=0,  quit = T_QUIT;
  char token=0;
  struct sockaddr_in from;
  int len = sizeof (from); 
/*  struct pollfd socket_info;
  struct sigaction act, oact; 
*/
  struct timeval dontwait;
  fd_set readset, writeset, exceptset;
  printf("enter check_io_signal_using_select\n");
  dontwait.tv_sec = 0; dontwait.tv_usec = 0;
  FD_ZERO(&readset);
  FD_ZERO(&writeset);
  FD_ZERO(&exceptset);
  FD_SET(socket_id, &readset);
  FD_SET(socket_id, &writeset);
  FD_SET(socket_id, &exceptset);

  /* instead of signal(SIGIO,SIG_IGN);  */
  act.sa_handler = SIG_IGN;
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;
  act.sa_flags |= SA_RESTART;
  sigaction(SIGIO,&act, &oact);

  if (select(socket_id+1,&readset,&writeset,&exceptset,&dontwait) > 0){
    printf("Something is ringing the socket!\n");
    if (FD_ISSET(socket_id, &exceptset)){
      printf("The exceptset case!\n");
    }    
    if (FD_ISSET(socket_id, &readset) || FD_ISSET(socket_id, &writeset)){   
      printf("The readset or writeset case!\n");
      soc = accept(socket_id,(struct sockaddr *)&from,&len);
#ifdef MPI
      printf("Node %d: return value of accept is %d\n",myid, soc);
#endif /* MPI */
      if(soc <= 0){
        printf("Connection not established, return code of accept is %d\n",soc);
      }
      else{
        ReadFull(soc,&socket_flag,1);
      }
    }
  }
  /* act.sa_handler = (void *)check_io_signal_using_select; */
  act.sa_handler = check_io_signal_using_select;
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;
  act.sa_flags |= SA_RESTART;
  sigaction(SIGIO,&act, &oact);
  /* signal(SIGIO,(void *)check_io_signal_using_select);  */
  printf("Returning from check_io_signal_using_select\n");
}

