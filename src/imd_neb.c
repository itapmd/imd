
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

/* auxiliary arrays */
real *pos=NULL, *pos_l=NULL, *pos_r=NULL, *f=NULL;

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
  if ((NULL==pos) || (NULL==pos_l) || (NULL==pos_r) || (NULL==f))
    error("cannot allocate NEB position array");
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
    for (k=0; k<NCELLS; k++) {
      cell *p = CELLPTR(k);
      for (i=0; i<p->n; i++) { 
        n = NUMMER(p,i);
        pos_l X(n) = ORT(p,i,X);
        pos_l Y(n) = ORT(p,i,Y);
        pos_l Z(n) = ORT(p,i,Z);
      }
    }
    /* compute and write energy of initial configuration */
    calc_forces(0);
    sprintf(outfilename, "%s.%02d", neb_outfilename, 0);
    write_eng_file_header();
    write_eng_file(0);
    fclose(eng_file);
    eng_file = NULL;
  }

  /* read positions of final configuration */
  if (neb_nrep-3==myrank) {
    sprintf(fname, "%s.%02d", infilename, neb_nrep-1);
    read_atoms(fname);
    if (NULL==pos) alloc_pos();
    for (k=0; k<NCELLS; k++) {
      cell *p = CELLPTR(k);
      for (i=0; i<p->n; i++) { 
        n = NUMMER(p,i);
        pos_r X(n) = ORT(p,i,X);
        pos_r Y(n) = ORT(p,i,Y);
        pos_r Z(n) = ORT(p,i,Z);
      }
    }
    /* compute and write energy of initial configuration */
    calc_forces(0);
    sprintf(outfilename, "%s.%02d", neb_outfilename, neb_nrep-1);
    write_eng_file_header();
    write_eng_file(0);
    fclose(eng_file);
    eng_file = NULL;
  }

  /* read positions of my configuration */
  sprintf(fname, "%s.%02d", infilename, myrank+1);
  read_atoms(fname);
  if (NULL==pos) alloc_pos();
  sprintf(outfilename, "%s.%02d", neb_outfilename, myrank+1);

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
  cpu_r = (neb_nrep - 3 == myrank) ? MPI_PROC_NULL : myrank + 1;

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
  real dl2=0.0, dr2=0.0, drl=0.0, d2=0.0, f2=0.0;
  real tmp, cosphi, fphi, src[2], dest[2], *d=pos;
  int k, i;

  /* exchange positions with neighbor replicas */
  neb_sendrecv_pos();

  /* compute tangent of current NEB path */
  for (i=0; i<DIM*natoms; i+=DIM) {

    vektor dl, dr;
    real x;

    /* distance to left and right replica */
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

    /* unnormalized tangent vector */
    d[i  ] = dr.x + dl.x;
    d[i+1] = dr.y + dl.y;
    d[i+2] = dr.z + dl.z;

    /* unmodified spring force */
    f[i  ] = dr.x - dl.x; 
    f[i+1] = dr.y - dl.y; 
    f[i+2] = dr.z - dl.z; 

    /* add up norms */
    dl2 += SPROD(dl,dl); 
    dr2 += SPROD(dr,dr); 
    drl += SPROD(dr,dl); 
    d2  += SPRODN(d+i,d+i); 
  }

  /* project internal force onto perpendicular direction */
  tmp = 0.0;
  for (k=0; k<NCELLS; k++) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; i++) { 
      int n = NUMMER(p,i);
      tmp += d X(n) * KRAFT(p,i,X);
      tmp += d Y(n) * KRAFT(p,i,Y);
      tmp += d Z(n) * KRAFT(p,i,Z);
    }
  }
  tmp /= d2;
  for (k=0; k<NCELLS; k++) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; i++) { 
      int n = NUMMER(p,i);
      KRAFT(p,i,X) -= tmp * d X(n);
      KRAFT(p,i,Y) -= tmp * d Y(n);
      KRAFT(p,i,Z) -= tmp * d Z(n);
      f2 += SPRODN( &KRAFT(p,i,X), &KRAFT(p,i,X) );
    }
  }

  /* estimate spring constant */
  src[0] = sqrt(dr2);
  src[1] = sqrt(f2);
  if (0==myrank) src[0] += sqrt(dl2);
  MPI_Allreduce( src, dest, 2, REAL, MPI_SUM, MPI_COMM_WORLD);
  neb_k = dest[1] / dest[0];

  /* project spring force onto parallel direction */
  cosphi = drl / sqrt(dl2 * dr2);
  if (cosphi > 0.0) fphi = 0.5 * (1.0 + cos(M_PI * cosphi));
  else fphi = 1.0;
  tmp   = (1 - fphi) * (dr2 - dl2) * neb_k / d2;
  fphi *= neb_k;
  for (k=0; k<NCELLS; k++) {
    cell *p = CELLPTR(k);
    for (i=0; i<p->n; i++) { 
      int n = NUMMER(p,i);
      KRAFT(p,i,X) += tmp * d X(n) + fphi * f X(n);
      KRAFT(p,i,Y) += tmp * d Y(n) + fphi * f Y(n);
      KRAFT(p,i,Z) += tmp * d Z(n) + fphi * f Z(n);
    }
  }
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

  /* write header */
  if (steps==0) {
    sprintf(fname, "%s.eng", neb_outfilename);
    neb_eng_file = fopen(fname,"a");
    if (NULL == neb_eng_file) 
      error_str("Cannot open properties file %s", fname);
    fprintf(neb_eng_file, "# nfc fnorm\n");
  }

  /* open .eng file if not yet open */
  if (NULL == neb_eng_file) {
    sprintf(fname, "%s.eng", neb_outfilename);
    neb_eng_file = fopen(fname,"a");
    if (NULL == neb_eng_file) 
      error_str("Cannot open properties file %s.eng", outfilename);
  }

  fprintf(neb_eng_file, "%d %e\n", nfc, neb_fnorm);

  /* flush .eng file every flush_int writes */
  if (flush_count++ > flush_int) {
    fflush(neb_eng_file);
    flush_count=0;
  }
}
