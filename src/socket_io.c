
/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#ifndef USE_SOCKETS
#define USE_SOCKETS
#endif

#include "imd.h" 


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

void check_socket(int steps) {

  int socket_flag;

  if (0==myid)
    socket_flag = connect_server();

#ifdef MPI
  MPI_Bcast(&socket_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif  

  if (socket_flag) {
    switch (socket_flag) {
      case T_QUIT:
        if (0==myid) error("Termination request received.");
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
        if (0==myid) printf("unknown token: %d\n",socket_flag);
        break;
    }
    if (0==myid)
      close_socket();
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
      coord.x = (int) (p->ort X(i) * scale.x) + xshift;
      coord.y = (int) (p->ort Y(i) * scale.y) + yshift;
      /* Check bounds */
      if (coord.x <  0   ) coord.x = 0;
      if (coord.x >= XRES) coord.x = XRES-1;
      if (coord.y <  0   ) coord.y = 0;
      if (coord.y >= YRES) coord.y = YRES-1;

      val = SPRODN(p->impuls,i,p->impuls,i) / (2*MASSE(p,i));

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
       WriteFull(soc,&p->nummer[i],sizeof(int)); 
       WriteFull(soc,&p->sorte[i],sizeof(int)); 
       WriteFull(soc,&p->masse[i],sizeof(real)); 
       WriteFull(soc,&p->ort X(i),sizeof(real)); 
       WriteFull(soc,&p->ort Y(i),sizeof(real)); 
#ifndef TWOD 
       WriteFull(soc,&p->ort Z(i),sizeof(real)); 
#endif 
       WriteFull(soc,&p->impuls X(i),sizeof(real)); 
       WriteFull(soc,&p->impuls Y(i),sizeof(real)); 
#ifndef TWOD
       WriteFull(soc,&p->impuls Z(i),sizeof(real));
#endif
#ifdef ORDPAR
#ifndef TWOD
       tmp = p->pot_eng[i];
       if (p->nbanz[i]>0) 
         tmp /= p->nbanz[i];
       else
         tmp=0.0;
       WriteFull(soc,&tmp,sizeof(real));
#endif
#endif
       WriteFull(soc,&p->pot_eng[i],sizeof(real));
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

#ifdef MPI
   if (0 == myid) { /* own data */
      k=0;
      for (j=0; j<ncells; ++j) {
        p = cell_array + CELLS(j);
        for (i=0; i < p->n; ++i) {
          if ( (p->ort X(i) >= socketwin_ll.x) &&
               (p->ort X(i) <= socketwin_ur.x) &&
               (p->ort Y(i) >= socketwin_ll.y) &&
               (p->ort Y(i) <= socketwin_ur.y) ) {
            nummer[k] = p->nummer[i];
            sorte[k]  = p->sorte[i];
            masse[k]  = p->masse[i];
            x[k]      = p->ort X(i);
            y[k]      = p->ort Y(i);
#ifndef TWOD
            z[k]      = p->ort Z(i);
#endif
            vx[k]     = p->impuls X(i);
            vy[k]     = p->impuls Y(i);
#ifndef TWOD
            vz[k]     = p->impuls Z(i);
#endif
            pot[k]    = p->pot_eng[i];
            k++;
          }
        }
      }

      /* Receive data from other cpus and write that */
      p = cell_array;  /* this is a pointer to the first (buffer) cell */
      for (m=1; m<num_cpus; ++m)
      for (j=0; j<ncells; ++j) {
        recv_cell(p,MPI_ANY_SOURCE,CELL_TAG);
        for (i=0; i < p->n; ++i) {
          if ( (p->ort X(i) >= socketwin_ll.x) &&
               (p->ort X(i) <= socketwin_ur.x) &&
               (p->ort Y(i) >= socketwin_ll.y) &&
               (p->ort Y(i) <= socketwin_ur.y) ) {
            nummer[k] = p->nummer[i];
            sorte[k]  = p->sorte[i];
            masse[k]  = p->masse[i];
            x[k]      = p->ort X(i);
            y[k]      = p->ort Y(i);
#ifndef TWOD
            z[k]      = p->ort Z(i);
#endif
            vx[k]     = p->impuls X(i);
            vy[k]     = p->impuls Y(i);
#ifndef TWOD
            vz[k]     = p->impuls Z(i);
#endif
            pot[k]    = p->pot_eng[i];
            k++;
	  }
	}
      }

   } else { /* myid != 0 */
     /* Send data to cpu 0 */
     for (j=0; j<ncells; j++) {
       p = cell_array + CELLS(j);
       send_cell(p,0,CELL_TAG);
     }
   }

#else /* #ifdef MPI */

   k=0;
   for (j=0; j<ncells; j++) {
     p = cell_array + CELLS(j);
     for (i = 0;i < p->n; ++i) {
       if (p->ort X(i) > socketwin_ll.x &&
           p->ort X(i) < socketwin_ur.x &&
           p->ort Y(i) > socketwin_ll.y &&
           p->ort Y(i) < socketwin_ur.y) {
         nummer[k] = p->nummer[i];
         sorte[k]  = p->sorte[i];
         masse[k]  = p->masse[i];
         x[k]      = p->ort X(i);
         y[k]      = p->ort Y(i);
#ifndef TWOD
         z[k]      = p->ort Z(i);
#endif
         vx[k]     = p->impuls X(i);
         vy[k]     = p->impuls Y(i);
#ifndef TWOD
         vz[k]     = p->impuls Z(i);
#endif
         pot[k]    = p->pot_eng[i];
         k++;
       }
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
  hist_t hist;
  unsigned char resolution_buffer[6];
  float f;

  hist.dim = dist_dim;
  make_histograms(&hist);

  if (0==myid) {

    resolution_buffer[0] = hist.dim.x / 256; /* high byte */
    resolution_buffer[1] = hist.dim.x % 256; /* low byte */
    resolution_buffer[2] = hist.dim.y / 256;
    resolution_buffer[3] = hist.dim.y % 256;
#ifndef TWOD
    resolution_buffer[4] = hist.dim.z / 256;
    resolution_buffer[5] = hist.dim.z % 256;
#endif

#ifdef TWOD
    if ((hist.dim.x > 65535) || (hist.dim.y > 65535))
#else
    if ((hist.dim.x > 65535) || (hist.dim.y > 65535) || (hist.dim.z > 65535))
#endif
      printf("Warning: sizes are larger than maximum size 65535!!!\n");

    f = (float)1.0;
    WriteFull(soc,(void *) &f, sizeof(float));
    WriteFull(soc,(void *) resolution_buffer, 2*DIM);
    WriteFull(soc,(void *) hist.pot_hist, hist.size * sizeof(float));
    WriteFull(soc,(void *) hist.kin_hist, hist.size * sizeof(float));
  }
}
