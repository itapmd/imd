
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#ifndef USE_SOCKETS
#define USE_SOCKETS
#endif

#include "imd.h" 
#include "socket_io.h"


/*****************************************************************************
*
* init_client, if the simulation acts as a client
*
*****************************************************************************/

void init_client() {

  fprintf(stderr, "baseport is %d\n", baseport);fflush(stderr);
  baseport = htons(baseport); /* we need it in network byte order */
  varIP = GetIP(display_host);
  if (0==myid) {
    if (varIP==0) {
      error("gethostbyname() failed, check display_host or specify IP number\n");
    } 
    else {
      printf("display_host is %s\n", display_host);
    }
  }
}


/*****************************************************************************
*
*  Try to connect to the visualization server
*
*****************************************************************************/

int connect_server() {

  unsigned char byte;
  int flag = 0;
  
  /* try to connect */
  soc=OpenClientSocket(varIP,baseport);
  if (soc > 0) {
    ReadFull(soc,&byte,1);
    flag = byte; /* convert to int */
    fprintf(stderr, "received socket_flag: %d\n", flag);fflush(stderr);
  }
  else {
    close_socket();
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

  if (0==myid) socket_flag = connect_server();

#ifdef MPI
  MPI_Bcast(&socket_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif  

  if (socket_flag) {
    switch (socket_flag) {
      case VIS_INIT:
        vis_init();
        break;
      case VIS_QUIT:
        if (0==myid) error("Termination request received.");
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
      case T_WRITE_RAS:
	write_ras_using_sockets();
        break;
      case T_WRITE_DISTRIB:
        write_distrib_using_sockets();
        break;
      case T_WRITE_CONF_SOCKET:
	write_conf_using_sockets();
	break;
#ifdef TWOD
      case T_WRITE_PICTURE:
        write_rgb_picture_to_socket();
        break;
#endif
      default:
        if (0==myid) printf("unknown token: %d\n", socket_flag);
        break;
    }
    if (0==myid) close_socket();
  }

}

/******************************************************************************
*
*  deliver initialization data
*
******************************************************************************/

void vis_init()
{
  unsigned char data[4];
  data[0] = PROTOCOL_VERSION_MAJOR;
  data[1] = PROTOCOL_VERSION_MINOR;
  data[2] = endian();
  data[3] = DIM;
  WriteFull( soc, (void *) data, 4 );
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
  if (0==myid) error("Termination request received.");
}

/******************************************************************************
*
*  write an error message to socket; to announce it, 
*  it is preceded by the negative of its length
*
******************************************************************************/

void vis_send_msg(char *msg)
{
  integer len;
  len = -(strlen(msg) + 1);
  WriteFull( soc, (void *) &len, sizeof(integer) );
  WriteFull( soc, (void *) msg, len );
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

#ifdef ORDPAR
  flags.nbanz = 1;
#endif

  /* loop over all atoms */
  for (k=0; k<ncells; k++) {
    cell *p;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; i++) {
      Ekin = SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / (2 * MASSE(p,i));
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
#ifdef ORDPAR
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
#ifdef ORDPAR
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
  WriteFull( soc, &flags,   ATOMS_FLAG_SIZE * sizeof(integer) );
  WriteFull( soc, &tot_min, ATOMS_FILT_SIZE * sizeof(float)   );
  WriteFull( soc, &tot_max, ATOMS_FILT_SIZE * sizeof(float)   );
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
#ifdef ORDPAR
      atlen += 1;
#else
      printf("Coordination number not available");
      stop = 1;
#endif
    }
    if (stop) atlen = -1;
  }

#ifdef MPI
  /* all CPUs need to know when to stop, so we distribute atlen */
  MPI_Bcast( &atlen, 1, INTEGER, 0, MPI_COMM_WORLD );
  /* broadcast flags and filters */
  if (atlen > 0) {
    MPI_Bcast(&at_send_flags, ATOMS_FLAG_SIZE, INTEGER,   0, MPI_COMM_WORLD);
    MPI_Bcast(&at_filt_flags, ATOMS_FLAG_SIZE, INTEGER,   0, MPI_COMM_WORLD);
    MPI_Bcast(&at_filt_min ,  ATOMS_FILT_SIZE, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&at_filt_max,   ATOMS_FILT_SIZE, MPI_FLOAT, 0, MPI_COMM_WORLD);
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

      Ekin = SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / (2 * MASSE(p,i));

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
#ifdef ORDPAR
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
#ifdef ORDPAR
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
  float zero = 0.0;

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
    if (0==myid) WriteFull( soc, &zero, sizeof(float) );
    return;
  }

  /* write or send own data */
  vis_write_atoms_fun();

#ifdef MPI
  /* receive and write foreign data */
  if (0==myid) {
    MPI_Status status;
    int m=1, len;
    do {
      MPI_Recv( sock_buf_at, SOCK_BUF_AT_SIZE, MPI_FLOAT, MPI_ANY_SOURCE, 
                MPI_ANY_TAG, cpugrid, &status );
      MPI_Get_count( &status, MPI_FLOAT, &len ); 
      len--; /* correct increased length */
      if ((status.MPI_TAG != AT_BUF_TAG+1) && (status.MPI_TAG != AT_BUF_TAG))
        error("messages mixed up");
      if (status.MPI_TAG == AT_BUF_TAG+1) m++;
      if (len>0) write_sock_buf_at( &len, AT_BUF_TAG );
    } while (m < num_cpus);
  }
  /* do not send other messages before we are finished */
  MPI_Barrier(cpugrid);
#endif /* MPI */

  /* the last block with zero atoms */
  if (0==myid) WriteFull( soc, (void *) &zero, sizeof(float) );

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
  }
#ifdef MPI
  MPI_Bcast(&par_group, 1, INTEGER, 0, MPI_COMM_WORLD);  
  MPI_Bcast(&flag     , 1, INTEGER, 0, MPI_COMM_WORLD);  
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
    }
  }
#ifdef MPI
  if (flag) MPI_Bcast(&deform_size, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif

  /* send back the changed parameters and other info */
  if (0==myid) {
    stp = steps;
    dsz = deform_size;
    WriteFull( soc, &stp, sizeof(integer) );
    WriteFull( soc, &dsz, sizeof(float)   );
  }
}
#endif

/*****************************************************************************
*
*  sends atomic configuration to socket 
*
*****************************************************************************/

void write_ras_using_sockets() 
{
   int r,s,t,i,k, ready; 
   int a; 
   cell *p; 
   real tmp;

   WriteFull(soc,&natoms,sizeof(int)); 
   /*  loop over all atoms */
   for (k=0; k<ncells; k++) {

     p = cell_array + CELLS(k);

     for (i=0; i<p->n; ++i) { 
       WriteFull(soc,&NUMMER(p,i),sizeof(integer)); 
       WriteFull(soc,&VSORTE(p,i),sizeof(shortint)); 
       WriteFull(soc,&MASSE(p,i),sizeof(real)); 
       WriteFull(soc,&ORT(p,i,X),sizeof(real)); 
       WriteFull(soc,&ORT(p,i,Y),sizeof(real)); 
#ifndef TWOD 
       WriteFull(soc,&ORT(p,i,Z),sizeof(real)); 
#endif 
       WriteFull(soc,&IMPULS(p,i,X),sizeof(real)); 
       WriteFull(soc,&IMPULS(p,i,Y),sizeof(real)); 
#ifndef TWOD
       WriteFull(soc,&IMPULS(p,i,Z),sizeof(real));
#endif
       WriteFull(soc,&POTENG(p,i),sizeof(real));
#ifdef ORDPAR
       WriteFull(soc,&NBANZ(p,i),sizeof(shortint)); 
#endif
     }
   }
}


/*****************************************************************************
*
*  sends atomic configuration to socket 
*
*****************************************************************************/

void write_conf_using_sockets() 
{
   int i, j, ready; 
   int a; 
   cell *p;
   int k=0, m;
   int *nummer;
   short int *sorte;
   real *masse, *x, *y, *vx, *vy, *pot;
   float f;
#ifndef TWOD
   real *z, *vz;
#endif

   /* allocs */
   nummer = (int *)       calloc(natoms, sizeof(int));
   sorte  = (short int *) calloc(natoms, sizeof(short int));
   masse  = (real *)    calloc(natoms, sizeof(real));
   x      = (real *)    calloc(natoms, sizeof(real));
   y      = (real *)    calloc(natoms, sizeof(real));
#ifndef TWOD
   z      = (real *)    calloc(natoms, sizeof(real));
#endif
   vx     = (real *)    calloc(natoms, sizeof(real));
   vy     = (real *)    calloc(natoms, sizeof(real));
#ifndef TWOD
   vz     = (real *)    calloc(natoms, sizeof(real));
#endif
   pot    = (real *)    calloc(natoms, sizeof(real));


   if (use_socket_window) { /* write all atoms in the box */    
     socketwin_ll.x = pic_ll.x;
     socketwin_ll.y = pic_ll.y;
     socketwin_ur.x = pic_ur.x;
     socketwin_ur.y = pic_ur.y;
   } else {              /* write only atoms in the window */    
     socketwin_ll.x = box_y.x;
     socketwin_ll.y = box_x.y;
     socketwin_ur.x = box_x.x;
     socketwin_ur.y = box_y.y;
   }

   /*  loop over all atoms */

   /* write own data */
   if (0 == myid) { 
      k=0;
      for (j=0; j<ncells; ++j) {
        p = cell_array + CELLS(j);
        for (i=0; i < p->n; ++i) {
          if ( (ORT(p,i,X) >= socketwin_ll.x) &&
               (ORT(p,i,X) <= socketwin_ur.x) &&
               (ORT(p,i,Y) >= socketwin_ll.y) &&
               (ORT(p,i,Y) <= socketwin_ur.y) ) {
            nummer[k] = NUMMER(p,i);
            sorte[k]  = VSORTE(p,i);
            masse[k]  = MASSE(p,i);
            x[k]      = ORT(p,i,X);
            y[k]      = ORT(p,i,Y);
#ifndef TWOD
            z[k]      = ORT(p,i,Z);
#endif
            vx[k]     = IMPULS(p,i,X);
            vy[k]     = IMPULS(p,i,Y);
#ifndef TWOD
            vz[k]     = IMPULS(p,i,Z);
#endif
            pot[k]    = POTENG(p,i);
            k++;
          }
        }
      }
   }

#ifdef MPI
   /* data of other CPUs */
   if (0 == myid) {
     /* Receive data from other cpus and write that */
      p = cell_array;  /* this is a pointer to the first (buffer) cell */
      for (m=1; m<num_cpus; ++m)
      for (j=0; j<ncells; ++j) {
        recv_cell(p,MPI_ANY_SOURCE,CELL_TAG);
        for (i=0; i < p->n; ++i) {
          if ( (ORT(p,i,X) >= socketwin_ll.x) &&
               (ORT(p,i,X) <= socketwin_ur.x) &&
               (ORT(p,i,Y) >= socketwin_ll.y) &&
               (ORT(p,i,Y) <= socketwin_ur.y) ) {
            nummer[k] = NUMMER(p,i);
            sorte[k]  = VSORTE(p,i);
            masse[k]  = MASSE(p,i);
            x[k]      = ORT(p,i,X);
            y[k]      = ORT(p,i,Y);
#ifndef TWOD
            z[k]      = ORT(p,i,Z);
#endif
            vx[k]     = IMPULS(p,i,X);
            vy[k]     = IMPULS(p,i,Y);
#ifndef TWOD
            vz[k]     = IMPULS(p,i,Z);
#endif
            pot[k]    = POTENG(p,i);
            k++;
          }
	}
        p->n=0;
      }
   } else { /* myid != 0 */
     /* Send data to cpu 0 */
     for (j=0; j<ncells; j++) {
       p = cell_array + CELLS(j);
       send_cell(p,0,CELL_TAG);
     }
   }
#endif

   /* there are k atoms to write */
   if (0==myid) {
     f=(float)1.0;
     WriteFull(soc,&f, sizeof(float));
     WriteFull(soc,&k,sizeof(int)); 
     WriteFull(soc,(void *) nummer, k*sizeof(int)); 
     WriteFull(soc,(void *) sorte, k*sizeof(short int)); 
     WriteFull(soc,(void *) masse, k*sizeof(real)); 
     WriteFull(soc,(void *) x, k*sizeof(real)); 
     WriteFull(soc,(void *) y, k*sizeof(real)); 
#ifndef TWOD
     WriteFull(soc,(void *) z, k*sizeof(real)); 
#endif
     WriteFull(soc,(void *) vx, k*sizeof(real)); 
     WriteFull(soc,(void *) vy, k*sizeof(real)); 
#ifndef TWOD
     WriteFull(soc,(void *) vz, k*sizeof(real)); 
#endif
     WriteFull(soc,(void *) pot, k*sizeof(real));
   }
}


/******************************************************************************
*
* write_distrib_using_sockets (only energy distributions)
*
******************************************************************************/

void write_distrib_using_sockets()
{
  hist_t Ekin_dist, Epot_dist;
  unsigned char resolution_buffer[6];
  int *flag;
  float f;

  Ekin_dist.dim = dist_dim; Epot_dist.dim = dist_dim;
  Ekin_dist.ur  = hist_ur;  Epot_dist.ur  = hist_ur;
  Ekin_dist.ll  = hist_ll;  Epot_dist.ll  = hist_ll;

  *flag = 0;
  make_distrib_select(&Epot_dist, 1, flag, dist_Epot_fun); 
  make_distrib_select(&Ekin_dist, 1, flag, dist_Ekin_fun); 

  if (0==myid) {

    resolution_buffer[0] = dist_dim.x / 256; /* high byte */
    resolution_buffer[1] = dist_dim.x % 256; /* low byte */
    resolution_buffer[2] = dist_dim.y / 256;
    resolution_buffer[3] = dist_dim.y % 256;
#ifndef TWOD
    resolution_buffer[4] = dist_dim.z / 256;
    resolution_buffer[5] = dist_dim.z % 256;
#endif

#ifdef TWOD
    if ((dist_dim.x > 65535) || (dist_dim.y > 65535))
#else
    if ((dist_dim.x > 65535) || (dist_dim.y > 65535) || (dist_dim.z > 65535))
#endif
      printf("Warning: sizes are larger than maximum size 65535!!!\n");

    f = (float)1.0;
    WriteFull(soc,(void *) &f, sizeof(float));
    WriteFull(soc,(void *) resolution_buffer, 2*DIM);
    WriteFull(soc,(void *) Epot_dist.dat, Epot_dist.size * sizeof(float));
    WriteFull(soc,(void *) Ekin_dist.dat, Ekin_dist.size * sizeof(float));
  }
}

#ifdef TWOD

/******************************************************************************
*
*  write_rgb_picture_to_socket writes pictures of configuration to sockets
*
******************************************************************************/

void write_rgb_picture_to_socket()
{
  int XRES, YRES;
  unsigned char image_resolution[4];
  vektor scale;
  real xshift, yshift;

  static shortint *redbyte   = NULL;
  static shortint *greenbyte = NULL;
  static shortint *bluebyte  = NULL;

#ifdef MPI
  static shortint *sum_red   = NULL;
  static shortint *sum_green = NULL;
  static shortint *sum_blue  = NULL;
#endif

  static unsigned char *buf  = NULL;

  int i,j,k,r,s;
  real val, delta, one_minus_delta;
  cell *p;
  real red,green,blue;
  FILE *out;
  /* real ecut; */
  ivektor2d coord;
  real tabred[5],tabgreen[5],tabblue[5];
  int ind, size, ii;

  /* get the current image size from socket */
  if (myid == 0) {
    ReadFull(soc, (void *) image_resolution, 4);
    XRES = (int) image_resolution[0] * 256 + (int) image_resolution[1];
    YRES = (int) image_resolution[2] * 256 + (int) image_resolution[3];
  }
#ifdef MPI
  MPI_Bcast( &XRES, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &YRES, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

  /* the dist bins are orthogonal boxes in space */
  size = XRES * YRES;
  scale.x = XRES / box_x.x;
  scale.y = YRES / box_y.y;

  /* Center image */
  xshift = (XRES - box_x.x * scale.x) / 2.0;
  yshift = (YRES - box_y.y * scale.y) / 2.0;

#ifdef MPI
  sum_red   = (shortint*) realloc(sum_red,   size*sizeof(shortint));
  sum_green = (shortint*) realloc(sum_green, size*sizeof(shortint));
  sum_blue  = (shortint*) realloc(sum_blue,  size*sizeof(shortint));
#endif
  redbyte   = (shortint*) realloc(redbyte,   size*sizeof(shortint));
  greenbyte = (shortint*) realloc(greenbyte, size*sizeof(shortint));
  bluebyte  = (shortint*) realloc(bluebyte,  size*sizeof(shortint));
  buf       = (unsigned char*) realloc(buf, 3*XRES);

  /* Zero bitmaps */
  for (i=0; i<size; i++ ) {
    bluebyte [i] = 0;
    greenbyte[i] = 0;
    redbyte  [i] = 0;
  }

  /* Color lookup table */
  tabred[0] = 0.10; tabgreen[0] = 0.20; tabblue[0] = 0.50;
  tabred[1] = 0.05; tabgreen[1] = 0.75; tabblue[1] = 0.75;
  tabred[2] = 0.10; tabgreen[2] = 0.50; tabblue[2] = 0.25;
  tabred[3] = 0.75; tabgreen[3] = 0.75; tabblue[3] = 0.05;
  tabred[4] = 0.75; tabgreen[4] = 0.05; tabblue[4] = 0.05;
  
  /* loop over all atoms */
  for (k=0; k<ncells; k++)
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) {
      coord.x = (int) (ORT(p,i,X) * scale.x) + xshift;
      coord.y = (int) (ORT(p,i,Y) * scale.y) + yshift;
      /* Check bounds */
      if (coord.x <  0   ) coord.x = 0;
      if (coord.x >= XRES) coord.x = XRES-1;
      if (coord.y <  0   ) coord.y = 0;
      if (coord.y >= YRES) coord.y = YRES-1;

      val = SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / (2*MASSE(p,i));

      val /= ecut_kin.y;
      if (1.0<val) val=0.9999;
      ind = (int)(val * 4.0);
      delta = (val * 4.0 - ind)/4.0; one_minus_delta = 1.0 - delta;

      red   = tabred  [ind+1] * delta + one_minus_delta * tabred  [ind];
      green = tabgreen[ind+1] * delta + one_minus_delta * tabgreen[ind];
      blue  = tabblue [ind+1] * delta + one_minus_delta * tabblue [ind];

      ii = coord.x * YRES + coord.y; 
      redbyte  [ii] = (shortint) (255 * red  );
      bluebyte [ii] = (shortint) (255 * blue );
      greenbyte[ii] = (shortint) (255 * green);
    }

#ifdef MPI
  /* Add the bytemaps */
  MPI_Reduce( redbyte,   sum_red,   size, SHORT, MPI_SUM, 0, cpugrid);
  MPI_Reduce( greenbyte, sum_green, size, SHORT, MPI_SUM, 0, cpugrid);
  MPI_Reduce( bluebyte , sum_blue,  size, SHORT, MPI_SUM, 0, cpugrid);

  /* clip bitmap to maximum */
  if (0==myid) 
  for (i=0; i<size; i++ ) {
    redbyte  [i] = sum_red  [i] < 255 ? sum_red  [i] : 255;
    greenbyte[i] = sum_green[i] < 255 ? sum_green[i] : 255;
    bluebyte [i] = sum_blue [i] < 255 ? sum_blue [i] : 255;
  }
#endif
  
  /* write bytes to socket */
  if (0==myid)
  for (i=0; i<YRES; i++ ) {
    for (j=0; j<XRES; j++ ) {
      ii = i * YRES + j;
      buf[3*j  ] = (unsigned char) redbyte  [ii];
      buf[3*j+1] = (unsigned char) greenbyte[ii];
      buf[3*j+2] = (unsigned char) bluebyte [ii];
    }
    WriteFull(soc,(void *)buf, 3*XRES);
  }
}

#endif  /* TWOD */
