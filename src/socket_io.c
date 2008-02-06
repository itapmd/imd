
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2008 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include <sys/ioctl.h>
#include "imd.h" 
#include "sockets.h"
#include "socket_io.h"

/*****************************************************************************
*
*  init_socket (executed only on MPI rank 0)
*
*****************************************************************************/

void init_socket() 
{
  static int done = 0;
  str255 hostname, msg;

  /* determine simulation host (the one with MPI rank 0 */
  if (gethostname(hostname,255)) {
    error("Cannot determine simulation host name");
  } else {
    fprintf(stderr, "Simulation host is %s\n", hostname); 
  }
 
  /* IMD acts as server */
  if ((server_socket) && (!done)) {
    /* open server socket */
    fprintf(stderr, "Opening server socket on port %d ...", server_port);
    soc = OpenServerSocket(htons(server_port));
    if (soc < 1) error("Cannot open server socket");
    else fprintf(stderr," connected\n");
  }

  /* IMD acts as client */
  if ((!server_socket) && (!done)) {

    /* determine display host */
    varIP = GetIP(display_host);
    if (varIP==0) {
      sprintf(msg, "Cannot find display_host %s", display_host);
      error(msg);
    }

    /* info message */
    if (client_port > 0)
      fprintf(stderr, "Connecting from port %d to %s port %d ...",
              client_port, display_host, server_port);
    else
      fprintf(stderr, "Connecting to %s port %d ...", 
              display_host, server_port);
    fflush(stderr);
   
    /* open client socket */
    soc=OpenClientSocket(varIP, htons(server_port), htons(client_port));
    if (soc < 1) error("Cannot open client socket");
    else fprintf(stderr," done\n");
  }

  /* do initialization only once */
  done = 1;
}


/*****************************************************************************
*
*  Connect to visualization (executed only on MPI rank 0)
*
*****************************************************************************/

int connect_visualization()
{
  unsigned char byte;
  int flag = 0;
  size_t nbytes;

  /* switch temporarily to non-blocking I/O on socket */
  /*
  fcntl(soc, F_SETFL, O_NONBLOCK);
  if (1==read(soc, &byte, 1)) {
    flag = byte; 
    fprintf(stderr, "received socket_flag: %d\n", flag);
  }
  fcntl(soc, F_SETFL, O_SYNC);
  */

  /* check whether there is something to read */
  ioctl(soc,FIONREAD,&nbytes);
  if (nbytes>0) {
    ReadFull(soc,&byte,1);
    flag=byte;
    fprintf(stderr, "received socket_flag: %d\n", flag);
  }

  return flag;
}


/*****************************************************************************
*
*  close socket
*
*****************************************************************************/

void close_socket() {
  shutdown(soc,2);
  close(soc);
}


/*****************************************************************************
*
*  check if there is a request on socket
*
*****************************************************************************/

void check_socket() {

  int socket_flag;

  if (0==myid) socket_flag = connect_visualization();

#ifdef MPI
  MPI_Bcast(&socket_flag, 1, MPI_INT, 0, cpugrid);
#endif  

  if (socket_flag) {
    switch (socket_flag) {
      case VIS_INIT:
        vis_init();
        break;
      case VIS_QUIT:
        if (0==myid) {
	  close_socket();
          error("Termination request received.");
	}
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
        if (0==myid) printf("unknown token: %d\n", socket_flag);
        break;
    }
  }
#ifdef MPI
  MPI_Barrier(cpugrid);
#endif
}

/******************************************************************************
*
*  deliver initialization data
*
******************************************************************************/

void vis_init()
{
  unsigned char data[4];
  if (0==myid) {
    data[0] = PROTOCOL_VERSION_MAJOR;
    data[1] = PROTOCOL_VERSION_MINOR;
    data[2] = endian();
    data[3] = DIM;
    WriteFull( soc, (void *) data, 4 );
  }
}

/******************************************************************************
*
*  write a final checkpoint and quit 
*
******************************************************************************/

void vis_write_config_quit()
{
  write_config(-1, steps);
#ifdef MPI
  MPI_Barrier(cpugrid);
#endif
  if (0==myid) {
    close_socket();
    error("Termination request received.");
  }
}

/******************************************************************************
*
*  inform on available data and the range of its values
*
******************************************************************************/

void vis_init_atoms()
{
  atoms_flag_t flags = {1,1,1,1,1,0};
  atoms_filt_t min, max, tot_min, tot_max;
  int first=1, k, i;
  float Ekin;

#ifdef NNBR
  flags.nbanz = 1;
#endif

  /* loop over all atoms */
  for (k=0; k<ncells; k++) {
    cell *p;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; i++) {
      Ekin = SPRODN(IMPULS,p,i,IMPULS,p,i) / (2 * MASSE(p,i));
      /* the first atom on this CPU; we hope there is at least one... */
      if (first) {
        min.sorte = VSORTE(p,i); 
	max.sorte = VSORTE(p,i); 
        min.x     = ORT(p,i,X); 
	max.x     = ORT(p,i,X); 
        min.y     = ORT(p,i,Y); 
	max.y     = ORT(p,i,Y); 
#ifndef TWOD
        min.z     = ORT(p,i,Z); 
	max.z     = ORT(p,i,Z); 
#else
        min.z     = 0.0; 
	max.z     = 0.0; 
#endif
        min.Ekin  = Ekin; 
	max.Ekin  = Ekin; 
        min.Epot  = POTENG(p,i); 
	max.Epot  = POTENG(p,i); 
#ifdef NNBR
        min.nbanz = NBANZ(p,i); 
	max.nbanz = NBANZ(p,i); 
#else
        min.nbanz = 0.0; 
	max.nbanz = 0.0; 
#endif
        first = 0;
      } 
      /* the remaining atoms on this CPU */
      else {
        min.sorte = MIN( min.sorte, VSORTE(p,i) ); 
	max.sorte = MAX( max.sorte, VSORTE(p,i) ); 
        min.x     = MIN( min.x,     ORT(p,i,X)  ); 
	max.x     = MAX( max.x,     ORT(p,i,X)  ); 
        min.y     = MIN( min.y,     ORT(p,i,Y)  ); 
	max.y     = MAX( max.y,     ORT(p,i,Y)  ); 
#ifndef TWOD
        min.z     = MIN( min.z,     ORT(p,i,Z)  ); 
	max.z     = MAX( max.z,     ORT(p,i,Z)  ); 
#endif
        min.Ekin  = MIN( min.Ekin,  Ekin        ); 
	max.Ekin  = MAX( max.Ekin,  Ekin        ); 
        min.Epot  = MIN( min.Epot,  POTENG(p,i) ); 
	max.Epot  = MAX( max.Epot,  POTENG(p,i) ); 
#ifdef NNBR
        min.nbanz = MIN( min.nbanz, NBANZ(p,i)  ); 
	max.nbanz = MAX( max.nbanz, NBANZ(p,i)  ); 
#endif
      }
    }
  }

  /* add up results of different CPUs */
#ifdef MPI
  MPI_Reduce( &min, &tot_min, ATOMS_FILT_SIZE, MPI_FLOAT, MPI_MIN, 0, cpugrid);
  MPI_Reduce( &max, &tot_max, ATOMS_FILT_SIZE, MPI_FLOAT, MPI_MAX, 0, cpugrid);
#else
  tot_min = min;
  tot_max = max;
#endif

  /* write the result to socket */
  if (0==myid) {
    WriteFull( soc, &flags,   ATOMS_FLAG_SIZE * sizeof(integer) );
    WriteFull( soc, &tot_min, ATOMS_FILT_SIZE * sizeof(float)   );
    WriteFull( soc, &tot_max, ATOMS_FILT_SIZE * sizeof(float)   );
#ifdef DEBUG
    printf("Sent atom flags:\n");
    printf("  sorte=%d ort=%d impuls=%d Ekin=%d Epot=%d nbanz=%d\n",
      flags.sorte, flags.ort,  flags.impuls, 
      flags.Ekin,  flags.Epot, flags.nbanz ); 
    printf("Sent atom min/max values:\n");
    printf("  Min: sorte=%f x=%f y=%f z=%f Ekin=%f Epot=%f nbanz=%f\n",
      min.sorte, min.x, min.y, min.z, min.Ekin, min.Epot, min.nbanz ); 
    printf("  Max: sorte=%f x=%f y=%f z=%f Ekin=%f Epot=%f nbanz=%f\n",
      max.sorte, max.x, max.y, max.z, max.Ekin, max.Epot, max.nbanz );
#endif
  }

}


/******************************************************************************
*
*  read and check flags for writing atoms to socket
*
******************************************************************************/

void vis_check_atoms_flags()
{
  int stop = 0;

  if (myid==0) {

    /* read send flags */
    ReadFull( soc, &at_send_flags, ATOMS_FLAG_SIZE * sizeof(integer) );

    /* read filter flags */
    ReadFull( soc, &at_filt_flags, ATOMS_FLAG_SIZE * sizeof(integer) );

    /* read filter values */
    ReadFull( soc, &at_filt_min,   ATOMS_FILT_SIZE * sizeof(float)   );
    ReadFull( soc, &at_filt_max,   ATOMS_FILT_SIZE * sizeof(float)   );

    /* check send_flags */
    atlen = 0; 
    if (at_send_flags.sorte ) atlen += 1;
    if (at_send_flags.ort   ) atlen += 3;
    if (at_send_flags.impuls) atlen += 3;
    if (at_send_flags.Ekin  ) atlen += 1;
    if (at_send_flags.Epot  ) atlen += 1;
    if (at_send_flags.nbanz ) {
#ifdef NNBR
      atlen += 1;
#else
      printf("Coordination number not available");
      stop = 1;
#endif
    }
    if (stop) atlen = -1;
#ifdef DEBUG
    printf("atlen is %d\n", atlen);
    printf("Received atom send flags:\n");
    printf("  sorte=%d ort=%d impuls=%d Ekin=%d Epot=%d nbanz=%d\n",
      at_send_flags.sorte, at_send_flags.ort,  at_send_flags.impuls, 
      at_send_flags.Ekin,  at_send_flags.Epot, at_send_flags.nbanz ); 
    printf("Received atom filter flags:\n");
    printf("  sorte=%d ort=%d impuls=%d Ekin=%d Epot=%d nbanz=%d\n",
      at_filt_flags.sorte, at_filt_flags.ort,  at_filt_flags.impuls, 
      at_filt_flags.Ekin,  at_filt_flags.Epot, at_filt_flags.nbanz ); 
    printf("Received atom filter values:\n");
    printf("  Min: sorte=%f x=%f y=%f z=%f Ekin=%f Epot=%f nbanz=%f\n",
      at_filt_min.sorte, at_filt_min.x,    at_filt_min.y, 
      at_filt_min.z,     at_filt_min.Ekin, at_filt_min.Epot, 
      at_filt_min.nbanz ); 
    printf("  Max: sorte=%f x=%f y=%f z=%f Ekin=%f Epot=%f nbanz=%f\n",
      at_filt_max.sorte, at_filt_max.x,    at_filt_max.y, 
      at_filt_max.z,     at_filt_max.Ekin, at_filt_max.Epot, 
      at_filt_max.nbanz );
#endif
  }

#ifdef MPI
  /* all CPUs need to know when to stop, so we distribute atlen */
  MPI_Bcast( &atlen, 1, INTEGER, 0, cpugrid );
  /* broadcast flags and filters */
  if (atlen > 0) {
    MPI_Bcast(&at_send_flags, ATOMS_FLAG_SIZE, INTEGER,   0, cpugrid);
    MPI_Bcast(&at_filt_flags, ATOMS_FLAG_SIZE, INTEGER,   0, cpugrid);
    MPI_Bcast(&at_filt_min ,  ATOMS_FILT_SIZE, MPI_FLOAT, 0, cpugrid);
    MPI_Bcast(&at_filt_max,   ATOMS_FILT_SIZE, MPI_FLOAT, 0, cpugrid);
  }
#endif

}

/******************************************************************************
*
*  if on CPU 0, write atoms to socket; otherwise, send buffer to CPU 0;
*  the last buffer is sent with a different tag
*
******************************************************************************/

void vis_write_atoms_buf(int *len, int tag)
{
  integer num;
  if (myid==0) {
    /* write only non-empty blocks */
    if (*len>0) {
      num = *len / atlen;
      WriteFull( soc, (void *) &num, sizeof(integer) );  
      WriteFull( soc, (void *) sock_buf_at, *len * sizeof(float) );
#ifdef DEBUG
      printf("Sent block of %d atoms\n", num);
#endif
    }
  }
#ifdef MPI
  else {
    /* we increase the length so that even an empty message is sent */
    sock_buf_at[*len] = 0.0;
    MPI_Send( sock_buf_at, *len+1, MPI_FLOAT, 0, tag, cpugrid );
  }
#endif
  *len=0;
}

/******************************************************************************
*
*  write atoms on a CPU to at buffer, and send the buffer
*
******************************************************************************/

void vis_write_atoms_fun()
{
  int i, k, len=0;
  cell *p;
  float Ekin;
  atoms_filt_t *min, *max;

  min = &at_filt_min;
  max = &at_filt_max;

  for (k=0; k<ncells; k++) {
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; i++) {

      Ekin = SPRODN(IMPULS,p,i,IMPULS,p,i) / (2 * MASSE(p,i));

      /* skip atom if it does not satisfy all filters */
      if (at_filt_flags.sorte) {
        if ((min->sorte > VSORTE(p,i)) || (max->sorte < VSORTE(p,i))) continue;
      }
      if (at_filt_flags.ort) {
        if ((min->x > ORT(p,i,X)) || (max->x < ORT(p,i,X))) continue;
        if ((min->y > ORT(p,i,Y)) || (max->y < ORT(p,i,Y))) continue;
#ifndef TWOD
        if ((min->z > ORT(p,i,Z)) || (max->z < ORT(p,i,Z))) continue;
#endif
      }
      if (at_filt_flags.Ekin) {
        if ((min->Ekin > Ekin) || (max->Ekin < Ekin)) continue;
      }
      if (at_filt_flags.Epot) {
        if ((min->Epot > POTENG(p,i)) || (max->Epot < POTENG(p,i))) continue;
      }
#ifdef NNBR
      if (at_filt_flags.nbanz) {
        if ((min->nbanz > NBANZ(p,i)) || (max->nbanz < NBANZ(p,i))) continue;
      }
#endif

      /* pack the requested data in buffer */
      if (at_send_flags.sorte) {
        sock_buf_at[len++] = VSORTE(p,i);
      }
      if (at_send_flags.ort) {
        sock_buf_at[len++] = ORT(p,i,X);
        sock_buf_at[len++] = ORT(p,i,Y);
#ifndef TWOD
        sock_buf_at[len++] = ORT(p,i,Z);
#else
        sock_buf_at[len++] = 0.0;
#endif
      }
      if (at_send_flags.impuls) {
        sock_buf_at[len++] = IMPULS(p,i,X);
        sock_buf_at[len++] = IMPULS(p,i,Y);
#ifndef TWOD
        sock_buf_at[len++] = IMPULS(p,i,Z);
#else
        sock_buf_at[len++] = 0.0;
#endif
      }
      if (at_send_flags.Ekin) {
        sock_buf_at[len++] = Ekin;
      }
      if (at_send_flags.Epot) {
        sock_buf_at[len++] = POTENG(p,i);
      }
#ifdef NNBR
      if (at_send_flags.nbanz) {
        sock_buf_at[len++] = NBANZ(p,i);
      }
#endif
    }

    /* send buffer if it is full */
    if (len > SOCK_BUF_AT_SIZE - 256) vis_write_atoms_buf( &len, AT_BUF_TAG );
  }

  /* send the last buffer */
  vis_write_atoms_buf( &len, AT_BUF_TAG+1 );
}

/******************************************************************************
*
*  write atoms to socket for visualization
*
******************************************************************************/

void vis_write_atoms()
{
  integer stp;
  integer zero = 0;

  /* get and distribute flags and filters */
  vis_check_atoms_flags();

  /* write step number and atlen */
  if (0==myid) {
    stp = steps;
    WriteFull( soc, &stp,   sizeof(integer) );
    WriteFull( soc, &atlen, sizeof(integer) );
  }

  /* return zero atoms if request makes no sense */
  if (atlen < 0) {
    if (0==myid) WriteFull( soc, &zero, sizeof(integer) );
#ifdef DEBUG
    if (0==myid) printf("Atom send request cannot be satisfied\n");
#endif
    return;
  }

  /* write or send own data */
  vis_write_atoms_fun();

#ifdef MPI
  /* receive and write foreign data */
  if (0==myid) {
    MPI_Status status;
    int m=1, len;
    while (m < num_cpus) {
      MPI_Recv( sock_buf_at, SOCK_BUF_AT_SIZE, MPI_FLOAT, MPI_ANY_SOURCE, 
                MPI_ANY_TAG, cpugrid, &status );
      MPI_Get_count( &status, MPI_FLOAT, &len ); 
      len--; /* correct increased length */
      if ((status.MPI_TAG != AT_BUF_TAG+1) && (status.MPI_TAG != AT_BUF_TAG))
        error("messages mixed up");
      if (status.MPI_TAG == AT_BUF_TAG+1) m++;
      if (len>0) vis_write_atoms_buf( &len, AT_BUF_TAG );
    }
  }
  /* do not send other messages before we are finished */
  MPI_Barrier(cpugrid);
#endif /* MPI */

  /* the last block with zero atoms */
  if (0==myid) WriteFull( soc, (void *) &zero, sizeof(integer) );
#ifdef DEBUG
  if (0==myid) printf("Sent last block with 0 atoms\n");
#endif
}

/*****************************************************************************
*
*  change or report parameters 
*
*****************************************************************************/

void vis_change_params()
{
  integer par_group, flag;

  /* read parameter group and flag */
  if (0==myid) {
    ReadFull( soc, &par_group, sizeof(integer) );
    ReadFull( soc, &flag,      sizeof(integer) );
#ifdef DEBUG
    printf("par_group=%d, change_flag=%d\n", par_group, flag);
#endif
  }
#ifdef MPI
  MPI_Bcast(&par_group, 1, INTEGER, 0, cpugrid);  
  MPI_Bcast(&flag     , 1, INTEGER, 0, cpugrid);  
#endif

  /* change parameters */
  switch (par_group) {
#if defined(DEFORM) || defined(HOMDEF)
    case VIS_PARAM_DEFORM:
      vis_change_params_deform(flag);
      break;
#endif
    default:
      if (0==myid) 
        printf("Unknown or unsupported parameter group: %d\n", par_group);
      break;
  }
}

/*****************************************************************************
*
*  change or report deform parameters 
*
*****************************************************************************/

#if defined(DEFORM) || defined(HOMDEF)
void vis_change_params_deform(integer flag)
{
  float   dsz;
  integer stp;

  /* if flag is set, we change the parameters */
  if (0==myid) {
    if (flag) {
      ReadFull( soc, &dsz, sizeof(float) );
      deform_size = dsz;
#ifdef DEBUG
      printf("Received new deform_size: %f\n", deform_size);
#endif
    }
  }
#ifdef MPI
  if (flag) MPI_Bcast(&deform_size, 1, REAL, 0, cpugrid);
#endif

  /* send back the changed parameters and other info */
  if (0==myid) {
    stp = steps;
    dsz = deform_size;
    WriteFull( soc, &stp, sizeof(integer) );
    WriteFull( soc, &dsz, sizeof(float)   );
#ifdef DEBUG
    printf("Sent changed parameters: deform_size = %f, current step = %d\n",
           deform_size, steps); 
#endif
  }
}
#endif

/*****************************************************************************
*
*  restart the simulation 
*
*****************************************************************************/

void vis_restart_simulation()
{
  /* decrease steps_max, so that we exit the main loop on the next iteration */
  steps_max=steps;
}
