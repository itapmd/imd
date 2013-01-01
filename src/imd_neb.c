
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2007 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_neb -- functions for the NEB method
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

#ifdef TWOD
#define nebSPRODN(x,y) ( (x)[0]*(y)[0] + (x)[1]*(y)[1] )
#else
#define nebSPRODN(x,y) ( (x)[0]*(y)[0] + (x)[1]*(y)[1] + (x)[2]*(y)[2] )
#endif

/* auxiliary arrays */
real *pos=NULL, *pos_l=NULL, *pos_r=NULL, *f=NULL, *tau=NULL, *dRleft=NULL, *dRright=NULL;

/******************************************************************************
*
*  initialize MPI (NEB version)
*
******************************************************************************/

void init_mpi(void)
{
  /* Initialize MPI */
  MPI_Comm_size(MPI_COMM_WORLD,&num_cpus);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  if (0 == myrank) { 
    printf("NEB: Starting up MPI with %d processes.\n", num_cpus);
  }
}

/******************************************************************************
*
*  shutdown MPI (NEB version)
*
******************************************************************************/

void shutdown_mpi(void)
{
  MPI_Barrier(MPI_COMM_WORLD);   /* Wait for all processes to arrive */
#ifdef MPELOG
  MPE_Log_sync_clocks();
#ifdef NO_LLMPE
  MPE_Finish_log( progname );
#endif
#endif
  MPI_Finalize();                /* Shutdown */
}

/******************************************************************************
*
*  allocate auxiliary arrays
*
******************************************************************************/

void alloc_pos(void) 
{
  pos   = (real *) malloc( DIM * natoms * sizeof(real ) );
  pos_l = (real *) malloc( DIM * natoms * sizeof(real ) );
  pos_r = (real *) malloc( DIM * natoms * sizeof(real ) );
  f     = (real *) malloc( DIM * natoms * sizeof(real ) );
  tau   = (real *) malloc( DIM * natoms * sizeof(real ) );
  dRleft= (real *) malloc( DIM * natoms * sizeof(real ) );
  dRright= (real *) malloc( DIM * natoms * sizeof(real ) );
  if ((NULL==pos) || (NULL==pos_l) || (NULL==pos_r) || (NULL==f)|| (NULL==tau)|| (NULL==dRleft) || (NULL==dRright))
    error("cannot allocate NEB position arrays");
}

/******************************************************************************
*
*  read all configurations (including initial and final)
*
******************************************************************************/

void read_atoms_neb(str255 infilename)
{
  str255 fname;
  int i, k, n;

  /* keep a copy of the outfile name without replica suffix */
  neb_outfilename = strdup(outfilename);

  /* read positions of initial configuration */
  if (0==myrank) {
    sprintf(fname, "%s.%02d", infilename, 0);
    myrank = 1;  /* avoid double info messages */
    read_atoms(fname);
    myrank = 0;
    alloc_pos();

    /* compute and write energy of initial configuration */
    calc_forces(0);
    neb_image_energies[0]=tot_pot_energy;
    sprintf(outfilename, "%s.%02d", neb_outfilename, 0);
    write_eng_file_header();
    write_eng_file(0);
    fclose(eng_file);
    eng_file = NULL;
  }

  /* read positions of final configuration */
  else if (neb_nrep-1==myrank) {
    sprintf(fname, "%s.%02d", infilename, neb_nrep-1);
    read_atoms(fname);
    if (NULL==pos) alloc_pos();

    /* compute and write energy of initial configuration */
    calc_forces(0);
    neb_image_energies[ neb_nrep-1]=tot_pot_energy;
    sprintf(outfilename, "%s.%02d", neb_outfilename, neb_nrep-1);
    write_eng_file_header();
    write_eng_file(0);
    fclose(eng_file);
    eng_file = NULL;
  }

  else
  {
      /* read positions of my configuration */
      sprintf(fname, "%s.%02d", infilename, myrank);
      printf("rank: %d reading  %s.%02d\n",myrank, infilename, myrank);fflush(stdout);
      read_atoms(fname);
      if (NULL==pos) alloc_pos();
      sprintf(outfilename, "%s.%02d", neb_outfilename, myrank);
  }
}

/******************************************************************************
*
*  exchange positions with neighbor replicas
*
******************************************************************************/

void neb_sendrecv_pos(void)
{
  int i, k, n, cpu_l, cpu_r;
  MPI_Status status;

  /* fill pos array */
  for (k=0; k<NCELLS; k++) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; i++) { 
      n = NUMMER(p,i);
      pos X(n) = ORT(p,i,X);
      pos Y(n) = ORT(p,i,Y);
      pos Z(n) = ORT(p,i,Z);
    }
  }

  /* ranks of left/right cpus */
  cpu_l = (0            == myrank) ? MPI_PROC_NULL : myrank - 1;
  cpu_r = (neb_nrep - 1 == myrank) ? MPI_PROC_NULL : myrank + 1;

  /* send positions to right, receive from left */
  MPI_Sendrecv(pos,   DIM*natoms, REAL, cpu_r, BUFFER_TAG,
	       pos_l, DIM*natoms, REAL, cpu_l, BUFFER_TAG,
	       MPI_COMM_WORLD, &status );

  /* send positions to left, receive from right */
  MPI_Sendrecv(pos,   DIM*natoms, REAL, cpu_l, BUFFER_TAG,
	       pos_r, DIM*natoms, REAL, cpu_r, BUFFER_TAG,
	       MPI_COMM_WORLD, &status );
}

/******************************************************************************
*
*  modify forces according to NEB
*
******************************************************************************/

void calc_forces_neb(void)
{
  real dl2=0.0, dr2=0.0, drl=0.0, d2=0.0, f2=0.0, f2max=0.0, drlmax=0.0,df=0.0;
  real tmp,tmp1,tmp2, cosphi, fphi, src[3], dest[3], *d=pos;
  real kr,kl;
  int k, i;
  int var_k=0;
 
  int myimage,maximage;
  real V_previous, V_actual, V_next;
  real deltaVmin,deltaVmax;
  real normdr,normdl,inormd;
  real Eref,Emax,Emin,delta_E;
  real ratio_plus,ratio_minus,abs_next,abs_previous ;
  real k_sum, k_diff,tmpl,tmpr;
  real tmp_neb_ks[NEB_MAXNREP] INIT(zero100);
  real felastfact=0.0;

  myimage = myrank;

  /* get info about the energies of the different images */
  neb_image_energies[ myimage]=tot_pot_energy;
  MPI_Allreduce(neb_image_energies , neb_epot_im, NEB_MAXNREP, REAL, MPI_SUM, MPI_COMM_WORLD);
  Emax=-999999999999999;
  Emin=999999999999999;
  for(i=0;i<neb_nrep;i++)
    {
      if(neb_epot_im[i]>=Emax)
	{
	  Emax=neb_epot_im[i];
	  maximage=i;
	}
      if(neb_epot_im[i]<=Emin)
	{
	  Emin=neb_epot_im[i];
	}
    }
  if(steps == neb_cineb_start)
    {
      if(neb_climbing_image > 0)
	{
	  if(myrank==0)
	    {
	      if( neb_climbing_image == maximage)
		printf("Starting climbing image = %d (= max_Epot = %lf)\n",neb_climbing_image, Emax);
	      else
		printf("Starting climbing image = %d \n WARNING: %d != %d with max_Epot = %lf)\n", \
		       neb_climbing_image,neb_climbing_image,maximage, Emax);
	    }
	}
      else
	{
	  neb_climbing_image = maximage;
	  if(myrank==0)
	    {
	      printf("Starting climbing image, image set to %d (= max_Epot = %lf)\n",maximage, Emax);
	    }
	}

    }

  /* determine variable spring constants (jcp113 p. 9901) */

  tmp_neb_ks[myimage]=0;

  if(myrank != 0 && myrank != neb_nrep-1)
  { 

    V_previous = neb_epot_im[myimage-1];
    V_actual   = neb_epot_im[myimage];
    V_next     = neb_epot_im[myimage+1];	

    if ( neb_kmax > 0 & neb_kmin >0 &&  steps > neb_vark_start)
      {
	var_k=1;
	k_sum  = neb_kmax + neb_kmin;
	k_diff = neb_kmax - neb_kmin;    
	delta_E = Emax - Emin;
	if (delta_E > 1.0e-12)
	  {
	    tmp_neb_ks[myimage] = 0.5 *(k_sum - k_diff * cos(3.141592653589793238*( neb_epot_im[myimage] - Emin )/delta_E ));
	  }
      }
    else
      {
	tmp_neb_ks[myimage] = neb_k;
      }
  }
  MPI_Allreduce(tmp_neb_ks , neb_ks, NEB_MAXNREP, REAL, MPI_SUM, MPI_COMM_WORLD); 

  /* exchange positions with neighbor replicas */
  neb_sendrecv_pos();
   
  /* determine tangent vector and the elastic spring force */
  if(myrank != 0 && myrank != neb_nrep-1)
  { 
      dl2=0.0;d2=0.0;dr2=0.0;
      
      kr = 0.5 * (neb_ks[myimage]+neb_ks[myimage+1]);
      kl = 0.5 * (neb_ks[myimage]+neb_ks[myimage-1]);

      /* preparation: calculate distance to left and right immage */
       for (i=0; i<DIM*natoms; i+=DIM) {
	 vektor dr,dl;	
	 real x;
	 dl.x = pos  [i  ] - pos_l[i  ];
	 dl.y = pos  [i+1] - pos_l[i+1];
	 dl.z = pos  [i+2] - pos_l[i+2];
	 dr.x = pos_r[i  ] - pos  [i  ];
	 dr.y = pos_r[i+1] - pos  [i+1];
	 dr.z = pos_r[i+2] - pos  [i+2];
	 
	    /* apply periodic boundary conditions */
	 if (1==pbc_dirs.x) {
	   x = - round( SPROD(dl,tbox_x) );
	   dl.x += x * box_x.x;
	   dl.y += x * box_x.y;
	   dl.z += x * box_x.z;
	   x = - round( SPROD(dr,tbox_x) );
	   dr.x += x * box_x.x;
	   dr.y += x * box_x.y;
	   dr.z += x * box_x.z;
	 }
	 if (1==pbc_dirs.y) {
	   x = - round( SPROD(dl,tbox_y) );
	   dl.x += x * box_y.x;
	   dl.y += x * box_y.y;
	   dl.z += x * box_y.z;
	   x = - round( SPROD(dr,tbox_y) );
	   dr.x += x * box_y.x;
	   dr.y += x * box_y.y;
	   dr.z += x * box_y.z;
	 }
	 if (1==pbc_dirs.z) {
	   x = - round( SPROD(dl,tbox_z) );
	   dl.x += x * box_z.x;
	   dl.y += x * box_z.y;
	   dl.z += x * box_z.z;
	   x = - round( SPROD(dr,tbox_z) );
	   dr.x += x * box_z.x;
	   dr.y += x * box_z.y;
	   dr.z += x * box_z.z;
	 }
	 dRleft[i  ] = dl.x; 
	 dRleft[i+1] = dl.y; 
	 dRleft[i+2] = dl.z; 
	 dRright[i  ] = dr.x; 
	 dRright[i+1] = dr.y; 
	 dRright[i+2] = dr.z; 
 	    
       }

      /* computation of the tangent requires 2 steps: determination of the direction and then normalization */
      /* here we use only the improved tangent method */

      if ( ( V_next > V_actual ) && ( V_actual > V_previous ) )
	{

	  for (i=0; i<DIM*natoms; i+=DIM) {
	    tau[i  ] = dRright[i  ];
	    tau[i+1] = dRright[i+1];
	    tau[i+2] = dRright[i+2];
	    d2  += dRright[i  ]*dRright[i  ];
	    d2  += dRright[i+1]*dRright[i+1];
	    d2  += dRright[i+2]*dRright[i+2];
	  }
	  
	  tmp=1.0/sqrt(d2);
	  for (i=0; i<DIM*natoms; i+=DIM) {
	    tau[i  ] *= tmp;
	    tau[i+1] *= tmp;
	    tau[i+2] *= tmp;
	    	 
	    if (var_k==1)
	      {
		felastfact += tau[i  ] * ( -kr *dRright[i  ] + kl * dRleft[i  ]);
		felastfact += tau[i+1] * ( -kr *dRright[i+1] + kl * dRleft[i+1]);
		felastfact += tau[i+2] * ( -kr *dRright[i+2] + kl * dRleft[i+2]);
	      }
	    else 
	      {
		felastfact +=  tau[i  ] * ( - dRright[i  ] + dRleft[i  ]);
		felastfact +=  tau[i+1] * ( - dRright[i+1] + dRleft[i+1]);
		felastfact +=  tau[i+2] * ( - dRright[i+2] + dRleft[i+2]);
	      }
	  }
	}
      else if ( ( V_next < V_actual ) && ( V_actual < V_previous ) ) 
	{
	  for (i=0; i<DIM*natoms; i+=DIM) {
	    tau[i  ] = dRleft[i  ];
	    tau[i+1] = dRleft[i+1];
	    tau[i+2] = dRleft[i+2];
	    d2  += dRleft[i  ]*dRleft[i  ];
	    d2  += dRleft[i+1]*dRleft[i+1];
	    d2  += dRleft[i+2]*dRleft[i+2];
	  }
	  
	  tmp=1.0/sqrt(d2);
	  for (i=0; i<DIM*natoms; i+=DIM) {
	    tau[i  ] *= tmp;
	    tau[i+1] *= tmp;
	    tau[i+2] *= tmp;
	    	  

	    if (var_k==1)
	      {
		felastfact += tau[i  ] * ( -kr *dRright[i  ] + kl * dRleft[i  ]);
		felastfact += tau[i+1] * ( -kr *dRright[i+1] + kl * dRleft[i+1]);
		felastfact += tau[i+2] * ( -kr *dRright[i+2] + kl * dRleft[i+2]);
	      }
	    else 
	      {
		felastfact +=  tau[i  ] * ( - dRright[i  ] + dRleft[i  ]);
		felastfact +=  tau[i+1] * ( - dRright[i+1] + dRleft[i+1]);
		felastfact +=  tau[i+2] * ( - dRright[i+2] + dRleft[i+2]);

	
	      }
	  }
	}      
      else
	{
	  abs_next     = FABS( V_next     - V_actual );
	  abs_previous = FABS( V_previous - V_actual );
	  deltaVmax    = MAX( abs_next, abs_previous );
	  deltaVmin    = MIN( abs_next, abs_previous );

	  for (i=0; i<DIM*natoms; i+=DIM) {
	    dr2  += dRright[i  ]*dRright[i  ];
	    dr2  += dRright[i+1]*dRright[i+1];
	    dr2  += dRright[i+2]*dRright[i+2];
	    dl2  += dRleft[i  ]*dRleft[i  ];
	    dl2  += dRleft[i+1]*dRleft[i+1];
	    dl2  += dRleft[i+2]*dRleft[i+2];
	  }
	  tmpl=1.0/sqrt(dl2);
	  tmpr=1.0/sqrt(dr2);

	  for (i=0; i<DIM*natoms; i+=DIM) {
	    vektor dl, dr;
	    dr.x =  dRright[i  ]*tmpr;
	    dr.y =  dRright[i+1]*tmpr;
	    dr.z =  dRright[i+2]*tmpr;

	    dl.x =  dRleft[i  ]*tmpl;
	    dl.y =  dRleft[i+1]*tmpl;
	    dl.z =  dRleft[i+2]*tmpl;

	    if (V_next > V_previous ) 
	      {
		tau[i  ] = dr.x * deltaVmax + dl.x * deltaVmin ;
		tau[i+1] = dr.y * deltaVmax + dl.y * deltaVmin;
		tau[i+2] = dr.z * deltaVmax + dl.z * deltaVmin;
	      }
	    else if ( V_next < V_previous ) 
	      {
		tau[i  ] = dr.x * deltaVmin + dl.x * deltaVmax ;
		tau[i+1] = dr.y * deltaVmin + dl.y * deltaVmax;
		tau[i+2] = dr.z * deltaVmin + dl.z * deltaVmax;
	      }
	    else
	      {
		tau[i  ] = dr.x + dl.x;
		tau[i+1] = dr.y + dl.y;
		tau[i+2] = dr.z + dl.z;
	      }
	    d2  += tau[i  ]*tau[i  ];
	    d2  += tau[i+1]*tau[i+1];
	    d2  += tau[i+2]*tau[i+2];
	  }
	  tmp=1.0/sqrt(d2);
	  for (i=0; i<DIM*natoms; i+=DIM) {	 
	    tau[i  ] *= tmp;
	    tau[i+1] *= tmp;
	    tau[i+2] *= tmp;
	    	   
	    if (var_k==1)
	      {
		felastfact += tau[i  ] * ( -kr *dRright[i  ] + kl * dRleft[i  ]);
		felastfact += tau[i+1] * ( -kr *dRright[i+1] + kl * dRleft[i+1]);
		felastfact += tau[i+2] * ( -kr *dRright[i+2] + kl * dRleft[i+2]);
	      }
	    else 
	      {
		felastfact +=  tau[i  ] * ( - dRright[i  ] + dRleft[i  ]);
		felastfact +=  tau[i+1] * ( - dRright[i+1] + dRleft[i+1]);
		felastfact +=  tau[i+2] * ( - dRright[i+2] + dRleft[i+2]);	
	      }
	    

	  }
	}

      /* finally construct the spring force */
      for (i=0; i<DIM*natoms; i+=DIM) {
	if (var_k==1)
	  {
	    f[i  ] = - tau[i  ] *felastfact;
	    f[i+1] = - tau[i+1] *felastfact;
	    f[i+2] = - tau[i+2] *felastfact;
	  }
	else
	  {
	    f[i  ] = - neb_k * tau[i  ] *felastfact;
	    f[i+1] = - neb_k * tau[i+1] *felastfact;
	    f[i+2] = - neb_k * tau[i+2] *felastfact;
	  }
      }
  }// end  if(myrank != 0 && mrank != neb_nrep-1)


  /* calculate the neb-force */
 if(myrank != 0 && myrank != neb_nrep-1)
  {

    // first scalar product of -force and tangent vector 
      tmp = 0.0;
      for (k=0; k<NCELLS; k++) {
	cell *p = CELLPTR(k);
	for (i=0; i<p->n; i++) { 
	  int n = NUMMER(p,i);
	  tmp -= tau X(n) * KRAFT(p,i,X);
	  tmp -= tau Y(n) * KRAFT(p,i,Y);
	  tmp -= tau Z(n) * KRAFT(p,i,Z);
	}
      }
     
      // add tmp times the tangent vector
      // and the spring force
      for (k=0; k<NCELLS; k++) {
	cell *p = CELLPTR(k);
	for (i=0; i<p->n; i++) { 
	  int n = NUMMER(p,i);
	  if(myimage == neb_climbing_image && (steps >= neb_cineb_start))
	    {
	      KRAFT(p,i,X) += 2.0*tmp * tau X(n);
	      KRAFT(p,i,Y) += 2.0*tmp * tau Y(n);
	      KRAFT(p,i,Z) += 2.0*tmp * tau Z(n);
	    }
	  else
	    {
	      KRAFT(p,i,X) += tmp * tau X(n) + f X(n);
	      KRAFT(p,i,Y) += tmp * tau Y(n) + f Y(n);
	      KRAFT(p,i,Z) += tmp * tau Z(n) + f Z(n);
	    }

	  
	}
      }
      
  } //  if(myrank != 0 && mrank != neb_nrep-1)
  
}

/******************************************************************************
*
*  write file with total fnorm, for monitoring convergence
*
******************************************************************************/
void write_neb_eng_file(int steps)
{
  static int flush_count=0;
  str255 fname;
  int i;

  /* write header */
  if (steps==0) {
    sprintf(fname, "%s.eng", neb_outfilename);
    neb_eng_file = fopen(fname,"a");
    if (NULL == neb_eng_file) 
      error_str("Cannot open properties file %s", fname);
    fprintf(neb_eng_file, "# nfc fnorm neb_k Epot_0 Epot_1 ... Epot_nrep\n");
  }

  /* open .eng file if not yet open */
  if (NULL == neb_eng_file) {
    sprintf(fname, "%s.eng", neb_outfilename);
    neb_eng_file = fopen(fname,"a");
    if (NULL == neb_eng_file) 
      error_str("Cannot open properties file %s.eng", outfilename);
  }

  fprintf(neb_eng_file, "%d %e %e   ", nfc, neb_fnorm, neb_k);
  for(i=0;i<neb_nrep;i++)
	{
	  fprintf(neb_eng_file,"%lf ", neb_epot_im[i]);
	}
  fprintf(neb_eng_file,"\n ");

  /* flush .eng file every flush_int writes */
  if (flush_count++ > flush_int) {
    fflush(neb_eng_file);
    flush_count=0;
  }
}

