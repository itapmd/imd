/******************************************************************************
* $RCSfile:
* $Revision:
* $Date:
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

  fprintf(stderr, "baseport is %d\n", baseport);
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
    fprintf(stderr, "received socket_flag: %d\n", flag);
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

#define XMAXRES 400
#define YMAXRES 400
  int XRES, YRES;
  unsigned char image_resolution[4];
  vektor scale;
  real xshift, yshift;

  short int redbyte[XMAXRES][YMAXRES];
  short int greenbyte[XMAXRES][YMAXRES];
  short int bluebyte[XMAXRES][YMAXRES];

#ifdef MPI
  short int sum_red[XMAXRES][YMAXRES];
  short int sum_green[XMAXRES][YMAXRES];
  short int sum_blue[XMAXRES][YMAXRES];
#endif

  unsigned char buf[3*YMAXRES];

  int i,j,r,s;
  real val, delta, one_minus_delta;
  cell *p;
  real red,green,blue;
  FILE *out;
  /* real ecut; */
  ivektor2d coord;
  real tabred[5],tabgreen[5],tabblue[5];
  int ind;

  /* get the current image size on the Display form socket */
#ifdef MPI
  if (myid == 0) {
#endif

    ReadFull(soc, (void *) image_resolution, 4);
    XRES = (int) image_resolution[0] * 256 + (int) image_resolution[1];
    YRES = (int) image_resolution[2] * 256 + (int) image_resolution[3];

#ifdef MPI
  }
  MPI_Bcast( &XRES, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast( &YRES, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

  /* the dist bins are orthogonal boxes in space */
  scale.x = XRES / box_x.x;
  scale.y = YRES / box_y.y;

  /* Center image */
  xshift = (XRES - box_x.x * scale.x) / 2.0;
  yshift = (YRES - box_y.y * scale.y) / 2.0;

  /* Zero bytemaps */
  for (i=0; i<XRES; i++ ) {
    for (j=0; j<YRES; j++ ) {
      bluebyte[i][j]  = 0;
      greenbyte[i][j] = 0;
      redbyte[i][j]   = 0;
    }
  }

  /* Color lookup table */
  tabred[0] = 0.10; tabgreen[0] = 0.20; tabblue[0] = 0.50;
  tabred[1] = 0.05; tabgreen[1] = 0.75; tabblue[1] = 0.75;
  tabred[2] = 0.10; tabgreen[2] = 0.50; tabblue[2] = 0.25;
  tabred[3] = 0.75; tabgreen[3] = 0.75; tabblue[3] = 0.05;
  tabred[4] = 0.75; tabgreen[4] = 0.05; tabblue[4] = 0.05;
  
  /* loop over all atoms */
  for ( r = cellmin.x; r < cellmax.x; ++r )
    for ( s = cellmin.y; s < cellmax.y; ++s ) {
      p = PTR_2D_V(cell_array, r, s,    cell_dim);
      for (i = 0;i < p->n; ++i) {
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
	
	red = tabred[ind+1] * delta + one_minus_delta * tabred[ind];
	green = tabgreen[ind+1] * delta + one_minus_delta * tabgreen[ind];
	blue = tabblue[ind+1] * delta + one_minus_delta * tabblue[ind];
	    
	redbyte[coord.x][coord.y]   = (short int) (255 * red  );
	bluebyte[coord.x][coord.y]  = (short int) (255 * blue );
	greenbyte[coord.x][coord.y] = (short int) (255 * green);
      }
    }

#ifdef MPI
  /* Add the bytemaps */
  MPI_Reduce( redbyte,   sum_red,   XRES * YRES, SHORT, MPI_SUM, 0, cpugrid);
  MPI_Reduce( greenbyte, sum_green, XRES * YRES, SHORT, MPI_SUM, 0, cpugrid);
  MPI_Reduce( bluebyte , sum_blue,  XRES * YRES, SHORT, MPI_SUM, 0, cpugrid);

  /* Zero bitmap */
  if (0==myid) 
    for (i=0; i<XRES; i++ ) {
      for (j=0; j<YRES; j++ ) {
	redbyte[i][j]   = sum_red[i][j]   < 255 ? sum_red[i][j]   : 255;
	greenbyte[i][j] = sum_green[i][j] < 255 ? sum_green[i][j] : 255;
	bluebyte[i][j]  = sum_blue[i][j]  < 255 ? sum_blue[i][j]  : 255;
      }
    }
#endif
  
  /* write bytes to socket */
#ifdef MPI
  if (0==myid) {
#endif
    for (i=0; i<YRES; i++ ) {
      for (j=0; j<XRES; j++ ) {
        buf[3*j  ] = (unsigned char) redbyte[i][j];
        buf[3*j+1] = (unsigned char) greenbyte[i][j];
        buf[3*j+2] = (unsigned char) bluebyte[i][j];
      }
      WriteFull(soc,(void *)buf, 3*XRES);
    }
#ifdef MPI
  }
#endif
}
#endif  /* TWOD */

/*****************************************************************************
*
*  sends atomic configuration to socket 
*
*****************************************************************************/

void write_ras_using_sockets() {
   int r,s,t,i, ready; 
   int a; 
   cell *p; 
   real tmp;

   WriteFull(soc,&natoms,sizeof(int)); 
   /*   WriteFull(soc,&box_x.x,sizeof(real)); 
   WriteFull(soc,&box_y.y,sizeof(real)); */
   /*  loop over all atoms */
   for ( r = cellmin.x; r < cellmax.x; ++r ) 
     for ( s = cellmin.y; s < cellmax.y; ++s ) 
#ifdef TWOD
       {
       p = PTR_2D_V(cell_array, r, s, cell_dim); 
#else
       for ( t = cellmin.z; t < cellmax.z; ++t ) { 
	   p = PTR_3D_V(cell_array, r, s, t, cell_dim); 
#endif
	   for (i = 0;i < p->n; ++i) { 
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
   int k=0;
   int *nummer;
   short int *sorte;
   real *masse, *x, *y, *vx, *vy, *pot;
   float f;
#ifndef TWOD
   real *z, *vz;
#endif
#ifdef MPI
   int m,tag;
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


   if (use_socket_window){ /* write all atoms in the box */    
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
        p = cell_array + CELLS(j);
        tag = CELL_TAG + j;
        recv_cell(p,m,tag);
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
       tag = CELL_TAG + j;
       send_cell(p,0,tag);
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


/*****************************************************************************
*
*  sends 2D or 3D kinetic and potential energy distribution to the socket 
*
*****************************************************************************/

void write_distrib_using_sockets()
{
  static size_t size;
  vektor scale;
  ivektor coord;
  cell *p;
  float *pot, *kin, f;
  shortint *num;
  int i,r,s,t;
  unsigned char resolution_buffer[6];
  static float    *pot_hist_local=NULL;
  static float    *kin_hist_local=NULL;
  static shortint *num_hist_local=NULL;
#ifdef MPI
  static float    *pot_hist_global=NULL;
  static float    *kin_hist_global=NULL;
  static shortint *num_hist_global=NULL;
#endif
  float *pot_hist, *kin_hist;
  shortint *num_hist;

#ifdef TWOD
  size = dist_dim.x * dist_dim.y;
#else
  size = dist_dim.x * dist_dim.y * dist_dim.z;
#endif

  /* allocate histogram arrays */
  if (NULL==pot_hist_local) {
    pot_hist_local = (float *) malloc(size*sizeof(float));
    if (NULL==pot_hist_local) 
      error("Cannot allocate distrib array.");
  }
  if (NULL==kin_hist_local) {
    kin_hist_local = (float *) malloc(size*sizeof(float));
    if (NULL==kin_hist_local) 
      error("Cannot allocate distrib array.");
  }
  if (NULL==num_hist_local) {
    num_hist_local = (shortint *) malloc(size*sizeof(shortint));
    if (NULL==num_hist_local) 
      error("Cannot allocate distrib array.");
  }
#ifdef MPI
  if (NULL==pot_hist_global) {
    pot_hist_global = (float *) malloc(size*sizeof(float));
    if (NULL==pot_hist_global) 
      error("Cannot allocate distrib array.");
  }
  if (NULL==kin_hist_global) {
    kin_hist_global = (float *) malloc(size*sizeof(float));
    if (NULL==kin_hist_global) 
      error("Cannot allocate distrib array.");
  }
  if (NULL==num_hist_global) {
    num_hist_global = (shortint *) malloc(size*sizeof(shortint));
    if (NULL==num_hist_global) 
      error("Cannot allocate distrib array.");
  }
#endif

  for (i=0; i<size; i++) {
    pot_hist_local[i]=0.0;
    kin_hist_local[i]=0.0;
    num_hist_local[i]=0;
  }

  /* the dist bins are orthogonal boxes in space */
  scale = box_x; 
  if (scale.x < box_y.x) scale.x = box_y.x; 
  if (scale.y < box_y.y) scale.y = box_y.y; 
#ifndef TWOD
  if (scale.z < box_y.z) scale.z = box_y.z; 
  
  if (scale.x < box_z.x) scale.x = box_z.x; 
  if (scale.y < box_z.y) scale.y = box_z.y; 
  if (scale.z < box_z.z) scale.z = box_z.z; 
#endif

  scale.x = dist_dim.x / scale.x;
  scale.y = dist_dim.y / scale.y;
#ifndef TWOD
  scale.z = dist_dim.z / scale.z;
#endif

  /* loop over all atoms */
  for ( r = cellmin.x; r < cellmax.x; ++r )
    for ( s = cellmin.y; s < cellmax.y; ++s )
#ifndef TWOD
      for ( t = cellmin.z; t < cellmax.z; ++t ) 
#endif
      {	
#ifdef TWOD
        p = PTR_2D_V(cell_array, r, s,    cell_dim);
#else
        p = PTR_3D_V(cell_array, r, s, t, cell_dim);
#endif
	for (i = 0;i < p->n; ++i) {
          coord.x = (int) (p->ort X(i) * scale.x);
          coord.y = (int) (p->ort Y(i) * scale.y);
#ifndef TWOD
          coord.z = (int) (p->ort Z(i) * scale.z);
#endif
          /* Check bounds */
          if (coord.x<0          ) coord.x = 0;
          if (coord.x>=dist_dim.x) coord.x = dist_dim.x-1;
          if (coord.y<0          ) coord.y = 0;
          if (coord.y>=dist_dim.y) coord.y = dist_dim.y-1;
#ifndef TWOD
          if (coord.z<0          ) coord.z = 0;
          if (coord.z>=dist_dim.z) coord.z = dist_dim.z-1;
#endif
	  pot = PTR_VV(pot_hist_local, coord, dist_dim);
	  kin = PTR_VV(kin_hist_local, coord, dist_dim);
          num = PTR_VV(num_hist_local, coord, dist_dim);
          (*num)++;
#ifdef DISLOC
          if (Epot_diff==1) { 
            *pot += p->pot_eng[i] - p->Epot_ref[i];
          } else
#endif
	  *pot += p->pot_eng[i];
	  *kin += SPRODN(p->impuls,i,p->impuls,i) / (2*MASSE(p,i));
        }
      }

#ifdef MPI
  MPI_Reduce(pot_hist_local,pot_hist_global,size,MPI_FLOAT,MPI_SUM,0,cpugrid);
  MPI_Reduce(kin_hist_local,kin_hist_global,size,MPI_FLOAT,MPI_SUM,0,cpugrid);
  MPI_Reduce(num_hist_local,num_hist_global,size,    SHORT,MPI_SUM,0,cpugrid);
  pot_hist=pot_hist_global;
  kin_hist=kin_hist_global;
  num_hist=num_hist_global;
#else
  pot_hist=pot_hist_local;
  kin_hist=kin_hist_local;
  num_hist=num_hist_local;
#endif

#ifdef MPI
  if (0==myid)
#endif
  {
    for (i=0; i<size; i++) {
      if (num_hist[i]>0) {
         pot_hist[i] /= num_hist[i];
         kin_hist[i] /= num_hist[i];
      }
    }

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
    f=(float)1.0;
    WriteFull(soc,(void *) &f, sizeof(float));
    WriteFull(soc,(void *) resolution_buffer, 2*DIM);
    WriteFull(soc,(void *) pot_hist, size*sizeof(float));
    WriteFull(soc,(void *) kin_hist, size*sizeof(float));
  }
}





