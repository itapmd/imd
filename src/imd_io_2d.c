/******************************************************************************
*
* imd_io_2d.c -- IO routines for the imd package 2D version
*
******************************************************************************/

/******************************************************************************
* $RCSfile$
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"


/******************************************************************************
*
* read_atoms - reads atoms and velocities into the cell-array
*
* The file format is flat ascii, one atom per line, lines beginning
* with '#' denote comments. Each line consists of
*
* number type mass x y z vx vy vz rest
*
* where
*
* number   is an arbitrary integer number assigned to each atom
* type     is the atom's index to the potenital table
* mass     is the mass of the atom
* x,y,z    are the atom's coordinates
* vx,vy,vz are the atom's velocity
* rest     is ignored until end of line
*
******************************************************************************/

void read_atoms(str255 infilename)
{
  cell *input;
  FILE *infile;
  char buf[512];
  int p;
  vektor2d pos;
  vektor2d vau;
#ifdef DISLOC
  FILE *reffile;
  int pref;
  char refbuf[512];
  real refeng, fdummy;
  int refn, idummy;
  vektor2d refpos;
#endif
  real m;
  int i,s,n;
  cell *to;
  ivektor2d cellc;
  int to_cpu;
  int addnumber = 0;

  /* allocate num_sort on all CPUs */
  if ((num_sort=calloc(ntypes,sizeof(int)))==NULL) {
      error("cannot allocate memory for num_sort\n");
  };

#ifdef MPI
  /* Try opening a per cpu file first when parallel_input is active */
  if (1==parallel_input) {
    sprintf(buf,"%s.%u",infilename,myid); 
    infile = fopen(buf,"r");
#ifdef DISLOC
  if (calc_Epot_ref == 0) {
    sprintf(buf,"%s.%u",reffilename,myid); 
    reffile = fopen(reffilename,"r");
  }
#endif
    /* When each cpu reads only part of the atoms, we have to add the
       number of atoms together to get the correct natoms. We set a
       flag here */
    if (NULL!=infile) addnumber=1;
  } else if (0!=myid) {
    recv_atoms(); 
    /* If CPU 0 found velocities in its data, no initialisation is done */
    MPI_Bcast( &natoms,       1, MPI_INT, 0, MPI_COMM_WORLD);  
    MPI_Bcast( num_sort, ntypes, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( &do_maxwell,   1, MPI_INT, 0, MPI_COMM_WORLD);
    return;
  };

  if ((1!=parallel_input) || (NULL==infile))
    infile = fopen(infilename,"r");
#ifdef DISLOC
  if ((1!=parallel_input) || (NULL==reffile))
    if (calc_Epot_ref == 0)
      reffile = fopen(reffilename,"r");
#endif

#else
  infile = fopen(infilename,"r");
#ifdef DISLOC
  if (calc_Epot_ref == 0)
    reffile = fopen(reffilename,"r");
#endif
#endif
  if (NULL==infile) error("Cannot open atoms file.");

#ifdef DISLOC
  if (calc_Epot_ref == 0)
    if (NULL==reffile) error("Cannot open reference file.");
#endif

  /* Set up 1 atom input cell */
  input = (cell *) malloc(sizeof(cell));
  if (0==input) error("Cannot allocate input cell.");
  input->n_max=0;
  alloc_cell(input, 1);

  natoms=0;
#ifdef SHOCK
  do_maxwell=1;
#else
  do_maxwell=0;
#endif

  /* Read the input file line by line */
  while(!feof(infile)) {

    buf[0] = (char) NULL;
    buf[0] = (char) NULL;
    fgets(buf,sizeof(buf),infile);
    while ('#'==buf[1]) fgets(buf,sizeof(buf),infile); /* eat comments */

#ifdef DISLOC
    refbuf[0] = (char) NULL;
    fgets(refbuf,sizeof(refbuf),reffile);
    while ('#'==refbuf[1]) fgets(refbuf,sizeof(refbuf),reffile); /* eat comments */
#endif

    /* Should use temporary variable */
#ifdef DISLOC
#ifdef DOUBLE
    p = sscanf(buf,"%d %d %lf %lf %lf %lf %lf %lf",
	      &n,&s,&m,&pos.x,&pos.y,&vau.x,&vau.y, &refeng);
    if (calc_Epot_ref == 0)
      pref = sscanf(refbuf,"%d %d %lf %lf %lf %lf",
		    &refn,&idummy,&fdummy,&refpos.x,&refpos.y,&refeng);
#else
    p = sscanf(buf,"%d %d %f %f %f %f %f %f",
	      &n,&s,&m,&pos.x,&pos.y,&vau.x,&vau.y, &refeng);  
    if (calc_Epot_ref == 0)
      pref = sscanf(refbuf,"%d %d %f %f %f %f",
		    &refn,&idummy,&fdummy,&refpos.x,&refpos.y,&refeng);
#endif
    if (calc_Epot_ref == 0)
      if (ABS(refn) != ABS(n)) error("Numbers in infile and reffile are different.\n");

#else
#ifdef DOUBLE
    p = sscanf(buf,"%d %d %lf %lf %lf %lf %lf",
	      &n,&s,&m,&pos.x,&pos.y,&vau.x,&vau.y);
#else
    p = sscanf(buf,"%d %d %f %f %f %f %f",
	      &n,&s,&m,&pos.x,&pos.y,&vau.x,&vau.y);  
#endif
#endif

#ifndef NOPBC
    pos = back_into_box(pos);
#endif

    if (0>=m) error("Mass zero or negative.\n");
    if (p>0) {
      switch( p ) {
      case(5):  /* n, m, s, ort */
	do_maxwell=1;
	input->n = 1;
	input->nummer[0] = n;
	input->sorte[0] = s;
	input->masse[0] = m;
	input->ort    X(0) = pos.x;
	input->ort    Y(0) = pos.y;
	input->impuls X(0) = 0;
	input->impuls Y(0) = 0;
	input->kraft  X(0) = 0;
	input->kraft  Y(0) = 0;
#ifdef DISLOC
	input->ort_ref X(0) = refpos.x;
	input->ort_ref Y(0) = refpos.y;
	input->Epot_ref[0] = refeng;
#endif
	break;
      case(7):  /* n, m, s, ort, vau */
	input->n = 1;
	input->nummer[0] = n;
	input->sorte[0] = s;
	input->masse[0] = m;
        input->ort    X(0) = pos.x;
	input->ort    Y(0) = pos.y;
	input->impuls X(0) = vau.x * m;
	input->impuls Y(0) = vau.y * m;
	input->kraft  X(0) = 0;
	input->kraft  Y(0) = 0;
#ifdef DISLOC
	input->ort_ref X(0) = refpos.x;
	input->ort_ref Y(0) = refpos.y;
	input->Epot_ref[0] = refeng;
#endif
break;
      default:
        printf("%d is where we find an...\n", natoms+1);
	error("Incomplete line in atoms file at line.\n");
      };

      cellc = cell_coord(pos.x,pos.y);

#ifdef MPI

      to_cpu = cpu_coord(cellc);
      
      if ((myid != to_cpu) && (1!=parallel_input)) {
        natoms++;
        num_sort[input->sorte[0]]++;
	MPI_Send( input->ort,     DIM, MPI_REAL, to_cpu, ORT_TAG,    cpugrid);
	MPI_Send( input->sorte,    1,  SHORT,    to_cpu, SORTE_TAG , cpugrid);
	MPI_Send( input->masse,    1,  MPI_REAL, to_cpu, MASSE_TAG , cpugrid);
	MPI_Send( input->nummer,   1,  INTEGER,  to_cpu, NUMMER_TAG, cpugrid);
	MPI_Send( input->impuls,  DIM, MPI_REAL, to_cpu, IMPULS_TAG, cpugrid);
#ifdef DISLOC
	MPI_Send( input->ort_ref, DIM, MPI_REAL, to_cpu, ORT_REF_TAG,cpugrid);
	MPI_Send( input->Epot_ref, 1,  MPI_REAL, to_cpu, POT_REF_TAG,cpugrid);
#endif
      } else if (to_cpu==myid) {  
        natoms++;
        num_sort[input->sorte[0]]++;
	cellc = local_cell_coord(pos.x,pos.y);
	move_atom(cellc, input, 0);
      };
#else
      natoms++;
      num_sort[input->sorte[0]]++;
      move_atom(cellc, input, 0);
#endif
    };
  };

  fclose(infile);  
#ifdef DISLOC
  if (calc_Epot_ref == 0)
    fclose(reffile);
#endif

#ifdef MPI

  /* Tell other CPUs that reading atoms is finished. (tag==0) */
  if (1!=parallel_input)
    for (s=1; s<num_cpus; s++)
      MPI_Ssend( input->ort, 2, MPI_REAL, s, 0, cpugrid);

  /* Add the number of atoms read (and kept) by each CPU */
  if (1==parallel_input) {
    MPI_Allreduce( &natoms, &addnumber, 1, MPI_INT, MPI_SUM, cpugrid);
    natoms = addnumber;
    for (i=0; i<ntypes; i++) {
      MPI_Allreduce(&num_sort[i], &addnumber, 1, MPI_INT, MPI_SUM, cpugrid);
      num_sort[i]=addnumber;
    };
  } else { /* broadcast */
    MPI_Bcast( &natoms ,      1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( num_sort, ntypes, MPI_INT, 0, MPI_COMM_WORLD);
  };

  /* If CPU 0 found velocities in its data, no initialisation is done */
  MPI_Bcast( &do_maxwell, 1, MPI_INT, 0, MPI_COMM_WORLD);

#endif /* MPI */

  /* print number of atoms */
  if (0==myid) {
    printf("Read structure with %d atoms.\n",natoms);
    addnumber=num_sort[0];
    printf("num_sort = [ %u",num_sort[0]);
    for (i=1; i<ntypes; i++) {
      printf(", %u",num_sort[i]);
      addnumber+=num_sort[i];
    };
    printf(" ],  total = %u\n",addnumber);
  };

}


#ifdef MPI

/******************************************************************************
*
*  recv_atoms
*
*  recveive atoms one at a time from CPU 0
*
*  this is only used when parallel_input==0
*
******************************************************************************/

void recv_atoms(void)
{
  cell *input, *target;
  MPI_Status status;
  ivektor2d cellc;
  ivektor2d local_cellc;

  printf("Node %d listening.\n",myid);

  input = (cell *) malloc(sizeof(cell));
  if (0==input) error("Cannot allocate input cell.");
  input->n_max = 0;
  alloc_cell(input,1);
  
  while ( 1 ) {

    MPI_Recv(input->ort, DIM, MPI_REAL, 0, MPI_ANY_TAG   , cpugrid, &status );

    if ((0 != status.MPI_TAG) && (ORT_TAG != status.MPI_TAG)) 
       error("Messages mixed up.");

    if ( 0 == status.MPI_TAG ) break;

    MPI_Recv(input->sorte,  1, SHORT,     0, SORTE_TAG , cpugrid, &status );
    MPI_Recv(input->masse,  1, MPI_REAL,  0, MASSE_TAG , cpugrid, &status );
    MPI_Recv(input->nummer, 1, INTEGER,   0, NUMMER_TAG, cpugrid, &status );
    MPI_Recv(input->impuls,DIM, MPI_REAL, 0, IMPULS_TAG, cpugrid, &status );
#ifdef DISLOC
    MPI_Recv(input->Epot_ref,1, MPI_REAL, 0, POT_REF_TAG, cpugrid, &status);
    MPI_Recv(input->ort_ref, 2, MPI_REAL, 0, ORT_REF_TAG, cpugrid, &status);
#endif

    local_cellc = local_cell_coord(input->ort X(0),input->ort Y(0));
    target = PTR_2D_VV(cell_array,local_cellc,cell_dim);

    /* See if we need some space */
    if (target->n >= target->n_max) alloc_cell(target,target->n_max+CSTEP);

    target->ort X(target->n)    = input->ort X(0);
    target->ort Y(target->n)    = input->ort Y(0);
    target->impuls X(target->n) = input->impuls X(0);
    target->impuls Y(target->n) = input->impuls Y(0);
    target->masse[target->n]  = input->masse[0];
    target->sorte[target->n]  = input->sorte[0];
    target->nummer[target->n] = input->nummer[0];
#ifdef DISLOC
    MPI_Recv(input->Epot_ref,1, MPI_REAL, 0, POT_REF_TAG, cpugrid, &status);
    MPI_Recv(input->ort_ref,    2, MPI_REAL, 0, ORT_REF_TAG, cpugrid, &status);
#endif

    ++target->n;

  };
  printf("Node %d leaves listen.\n",myid);
}

#endif 


/******************************************************************************
*
* write_properties writes selected properties to *.eng file
*
******************************************************************************/

void write_properties(int steps)
{
  FILE *out;
  str255 fname;
  int i;
  real  vol; 
  real part_kin_energy;
  real part_pot_energy; 

  /* Energiefile schreiben */
  sprintf(fname,"%s.eng",outfilename);

  /* Groessen pro Tln ausgeben */
#ifdef MC
  part_pot_energy = mc_epot_part();
#else
  part_pot_energy = tot_pot_energy / natoms;
  part_kin_energy = tot_kin_energy / natoms;
#endif
  vol = volume / natoms;

  out = fopen(fname,"a");
  if (NULL == out) error("Cannot open properties file.");

  fprintf(out,"%10.4e", (double)(steps * timestep));
  fprintf(out," %10.4e", (double)part_pot_energy);
#ifndef MC
  fprintf(out," %10.4e", (double)part_kin_energy);
  fprintf(out," %10.4e", (double)pressure);
  /* heat_cond is now written to the file tempdist
     #ifdef TRANSPORT
     fprintf(out," %10.4e", (double)heat_cond); 
     #endif */
#else
  fprintf(out," %10.4e", (double)(mc_accept/(real)mc_count));
  mc_accept = (real)0;
  mc_count  = 0;
#endif
  fprintf(out," %10.4e", (double)vol);

#ifdef PAXTEST
  if (ensemble==ENS_NPT_AXIAL) {
    fprintf(out," %10.4e %10.4e", 
                  (double) stress.x, (double) stress.y );
    fprintf(out," %10.4e %10.4e", 
                  (double) box_x.x,  (double) box_y.y  );
  };
#endif

  putc('\n',out);

  fclose(out);

}


#ifdef MSQD

/******************************************************************************
*
* write_msqd writes mean square displacement to *.msqd file
*
******************************************************************************/

void write_msqd(int steps)
{
  FILE *out;
  str255 fname;
  int i;

  sprintf(fname,"%s.msqd",outfilename);
  out = fopen(fname,"a");
  if (NULL == out) error("Cannot open msqd file.");

  fprintf(out, "%10.4e", (double)(steps * timestep));
  for (i=0; i<ntypes; i++) {
    fprintf(out," %10.4e", (double)(msqd_global[i] / num_sort[i]));
  };
  putc('\n',out);

  fclose(out);

}

#endif


/******************************************************************************
*
* write_config writes a configuration to a numbered file
* also creates a checkpoint
*
******************************************************************************/

void write_config(int steps)

/* Makro to write data of cell p to file out */
#ifdef ZOOM
/* write only atoms inside the picture frame */
#define WRITE_CELL     for (i = 0;i < p->n; ++i) \
            if( (p->ort X(i) >= pic_ll.x) && (p->ort X(i) <= pic_ur.x) && \
	        (p->ort Y(i) >= pic_ll.y) && (p->ort Y(i) <= pic_ur.y) ) \
             fprintf(out,"%d %d %12f %12f %12f %12f %12f %12f\n",\
	     p->nummer[i],\
	     p->sorte[i],\
	     p->masse[i],\
	     p->ort X(i),\
	     p->ort Y(i),\
	     p->impuls X(i) / p->masse[i],\
	     p->impuls Y(i) / p->masse[i],\
             p->pot_eng[i])
#else
#define WRITE_CELL     for (i = 0;i < p->n; ++i) \
             fprintf(out,"%d %d %12f %12f %12f %12f %12f %12f\n",\
	     p->nummer[i],\
	     p->sorte[i],\
	     p->masse[i],\
	     p->ort X(i),\
	     p->ort Y(i),\
	     p->impuls X(i) / p->masse[i],\
	     p->impuls Y(i) / p->masse[i],\
             p->pot_eng[i])
#endif

{ 
  FILE *out;
  str255 fname;
  int fzhlr;
  cell *p,*q;
  int i,j,k,l,m,tag;

  /* Dateiname fuer Ausgabedatei erzeugen */
#ifdef SHEAR
  fzhlr = steps;
#else
  fzhlr = steps / rep_interval;
#endif

#ifdef MPI  
  if (1==parallel_output)
    sprintf(fname,"%s.%u.%u",outfilename,fzhlr,myid);
  else
#endif
    sprintf(fname,"%s.%u",outfilename,fzhlr);

#ifdef MPI

  if (1==parallel_output) {

    /* Ausgabedatei oeffnen */
    out = fopen(fname,"w");
    if (NULL == out) error("Cannot open output file for config.");

    for (j = 1; j < cell_dim.x-1; ++j )
      for (k = 1; k < cell_dim.y-1; ++k ) {
 	  p = PTR_2D_V(cell_array, j, k, cell_dim);
	  WRITE_CELL;
	};
    
    fclose(out);

  } else { 

    if (0==myid) {

      /* Ausgabedatei oeffnen */
      out = fopen(fname,"w");
      if (NULL == out) error("Cannot open output file for config.");

      /* Write data on CPU 0 */

      /* Write own data */
      for (j = 1; j < cell_dim.x-1; ++j )
	for (k = 1; k < cell_dim.y-1; ++k ) {
	  p = PTR_2D_V(cell_array, j, k, cell_dim);
	  WRITE_CELL;
	};

      /* Receive data from other cpus and write that */
      p   = PTR_2D_V(cell_array, 0, 0, cell_dim);
      for ( m = 1; m < num_cpus; ++m)
	for (j = 1; j < cell_dim.x-1; ++j )
	  for (k = 1; k < cell_dim.y-1; ++k ) {
	    tag = PTR_2D_V(CELL_TAG, j, k, cell_dim);
	    recv_cell( p, m, tag );
	    WRITE_CELL;
	  };

      fclose(out);      
    } else { 
      /* Send data to cpu 0 */
      for (j = 1; j < cell_dim.x-1; ++j )
	for (k = 1; k < cell_dim.y-1; ++k ) {
	  p   = PTR_2D_V(cell_array, j, k, cell_dim);
	  tag = PTR_2D_V(CELL_TAG, j, k, cell_dim);
	  send_cell( p, 0, tag );
	}
    }
  }

#else

  /* Ausgabedatei oeffnen */
  out = fopen(fname,"w");
  if (NULL == out) error("Cannot open output file for config.");

  for (p = cell_array; 
       p <= PTR_2D_V(cell_array,
		     cell_dim.x-1,
		     cell_dim.y-1,
		     cell_dim);
       ++p ) 

    WRITE_CELL;

  fclose(out);  

#endif

  /* write iteration file */
#ifdef MPI
  if (myid == 0) {
#endif
  sprintf(fname,"%s.%u.itr",outfilename,fzhlr);

  out = fopen(fname,"w");
  if (NULL == out) error("Cannot write checkpoint file.");

  fprintf(out,"startstep \t%d\n",steps);
  fprintf(out,"box_x \t%f %f\n",box_x.x,box_x.y);
  fprintf(out,"box_y \t%f %f\n",box_y.x,box_y.y);
  fprintf(out,"starttemp \t%f\n",temperature);
#ifdef NPT
  if (ensemble==ENS_NPT_ISO) {
    fprintf(out,"pressure_ext \t%f\n",pressure_ext.x);
  };
  if (ensemble==ENS_NPT_AXIAL) {
    fprintf(out,"pressure_ext \t%f %f\n",pressure_ext.x,pressure_ext.y);
  };
#endif

  fclose(out);

#ifdef MPI
  }
#endif

}


#ifdef DISLOC

/******************************************************************************
*
* write_demmaps writes a differential energy map to file *.dem.x
*
******************************************************************************/

void write_demmaps(int steps)

#define WRITE_CELL_DEM     for (i = 0;i < p->n; ++i) {\
             if (p->sorte[i] == dpotsorte) {\
	       dpot = ABS(p->pot_eng[i] - p->Epot_ref[i]);\
               if (dpot > min_dpot)\
                 fprintf(demout,"%12f %12f %12f\n",\
	         p->ort X(i),\
	         p->ort Y(i),\
                 dpot);\
            }\
          }
{
  FILE *demout;
  str255 demfname;
  int fzhlr;
  cell *p,*q;
  int i,j,k,m,tag;
  real dpot;

  /* Dateiname fuer Ausgabedatei erzeugen */
  fzhlr = steps;
  sprintf(demfname,"%s.dem.%u",outfilename,fzhlr);
#ifdef MPI

  if (0==myid) {

    /* Ausgabedatei oeffnen */
    demout = fopen(demfname,"w");
    if (NULL == demout) error("Cannot open output file for dem.");

    /* Write data on CPU 0 */

    /* Write own data */
    for (j = 1; j < cell_dim.x-1; ++j )
      for (k = 1; k < cell_dim.y-1; ++k ) { 
	p = PTR_2D_V(cell_array, j, k, cell_dim);
	WRITE_CELL_DEM;
      };

    /* Receive data from other cpus and write that */
    p   = PTR_2D_V(cell_array, 0, 0, cell_dim);
    for ( m = 1; m < num_cpus; ++m)
      for (j = 1; j < cell_dim.x-1; ++j )
	for (k = 1; k < cell_dim.y-1; ++k ) {
	  tag = PTR_2D_V(CELL_TAG, j, k, cell_dim);
	  recv_cell( p, m, tag );
	  WRITE_CELL_DEM;
	};

    fclose(demout);      

  } else { 
    /* Send data to cpu 0 */
    for (j = 1; j < cell_dim.x-1; ++j )
      for (k = 1; k < cell_dim.y-1; ++k ) {
        p   = PTR_2D_V(cell_array, j, k, cell_dim);
        tag = PTR_2D_V(CELL_TAG, j, k, cell_dim);
        send_cell( p, 0, tag );
      };
  };

#else

  /* Ausgabedatei oeffnen */
  demout = fopen(demfname,"w");
  if (NULL == demout) error("Cannot open output file for dem.");

  for (p = cell_array; 
       p <= PTR_2D_V(cell_array,
		     cell_dim.x-1,
		     cell_dim.y-1,
		     cell_dim);
       ++p ) {
    WRITE_CELL_DEM;
    }

  fclose(demout);

#endif
  
}

/******************************************************************************
*
* write_dspmaps writes a differential displacement map to file *.ddm.x
*
******************************************************************************/

void write_dspmaps(int steps)

#ifdef MONOLJ
#define WRITE_CELL     for (i = 0;i < p->n; ++i) \
             fprintf(out,"%12f %12f %12f %12f\n",\
             p->ort X(i),\
             p->ort Y(i),\
	     p->impuls X(i),\
	     p->impuls Y(i));
		     
#else

#define WRITE_CELL_DSP     for (i = 0;i < p->n; ++i) {\
             dx = p->ort X(i) - p->ort_ref X(i);\
             dy = p->ort Y(i) - p->ort_ref Y(i);\
             fprintf(dspout,"%12f %12f %12f %12f\n",\
	     p->ort X(i),\
	     p->ort Y(i),\
	     dx,\
             dy);\
          }
#endif

{
  FILE *dspout;
  str255 dspfname;
  int fzhlr;
  cell *p,*q;
  int i,j,k,m,tag;
  real dx, dy, boxx, boxy;

  /* Dateiname fuer Ausgabedatei erzeugen */
  fzhlr = steps;
  sprintf(dspfname,"%s.dsp.%u",outfilename,fzhlr);

  /* 1/2 of the boxlength, this is not correct for non-cartesian boxes ! */
  boxx = box_x.x/2;
  boxy = box_y.y/2;
#ifdef MPI

  if (0==myid) {

  /* Ausgabedatei oeffnen */
    dspout = fopen(dspfname,"w");
    if (NULL == dspout) error("Cannot open output file for dsp.");

    /* Write data on CPU 0 */

    /* Write own data */
    for (j = 1; j < cell_dim.x-1; ++j )
      for (k = 1; k < cell_dim.y-1; ++k ) {
	p = PTR_2D_V(cell_array, j, k, cell_dim);
	WRITE_CELL_DSP;
      }

    /* Receive data from other cpus and write that */
    p   = PTR_2D_V(cell_array, 0, 0, cell_dim);
    for ( m = 1; m < num_cpus; ++m)
      for (j = 1; j < cell_dim.x-1; ++j )
        for (k = 1; k < cell_dim.y-1; ++k ) {
          tag = PTR_2D_V(CELL_TAG, j, k, cell_dim);
          recv_cell( p, m, tag );
          WRITE_CELL_DSP;
	}

    fclose(dspout);      
    
  } else { 
  /* Send data to cpu 0 */
  for (j = 1; j < cell_dim.x-1; ++j )
    for (k = 1; k < cell_dim.y-1; ++k ) {
      p   = PTR_2D_V(cell_array, j, k, cell_dim);
      tag = PTR_2D_V(CELL_TAG, j, k, cell_dim);
      send_cell( p, 0, tag );
    };
  };

#else

  /* Ausgabedatei oeffnen */
  dspout = fopen(dspfname,"w");
  if (NULL == dspout) error("Cannot open output file for dsp.");

  for (p = cell_array; 
       p <= PTR_2D_V(cell_array,
		     cell_dim.x-1,
		     cell_dim.y-1,
		     cell_dim);
       ++p ) {
    WRITE_CELL_DSP;

  }

  fclose(dspout);

#endif
}


/******************************************************************************
*
* update_ort_ref updates ort_ref
*
******************************************************************************/

void update_ort_ref(void)

{ 
  cell *p,*q;
  int i,j,k,m,tag;
  real dx, dy, boxx, boxy;

  printf("hallo");fflush(stdout);
#ifdef MPI

  if (0==myid) {

  /* Update own data */
  for (j = 1; j < cell_dim.x-1; ++j )
    for (k = 1; k < cell_dim.y-1; ++k ) 
      p = PTR_2D_V(cell_array, j, k, cell_dim);
      for (i = 0;i < p->n; ++i) {
	p->ort_ref X(i) = p->ort X(i);
	p->ort_ref Y(i) = p->ort Y(i);
      }

  /* Receive data from other cpus and Update that */
  p = PTR_2D_V(cell_array, 0, 0, cell_dim);
  for ( m = 1; m < num_cpus; ++m)
    for (j = 1; j < cell_dim.x-1; ++j )
      for (k = 1; k < cell_dim.y-1; ++k ) {
	tag = PTR_2D_V(CELL_TAG, j, k, cell_dim);
        recv_cell( p, m, tag );
	  for (i = 0;i < p->n; ++i) {
	    p->ort_ref X(i) = p->ort X(i);
	    p->ort_ref Y(i) = p->ort Y(i);
	  }
	};

  } else { 
  /* Send data to cpu 0 */
  for (j = 1; j < cell_dim.x-1; ++j )
    for (k = 1; k < cell_dim.y-1; ++k ) {
      p   = PTR_2D_V(cell_array, j, k, cell_dim);
      tag = PTR_2D_V(CELL_TAG, j, k, cell_dim);
      send_cell( p, 0, tag );
    };
  };

#else

  for (p = cell_array; 
       p <= PTR_2D_V(cell_array,
		     cell_dim.x-1,
		     cell_dim.y-1,
		     cell_dim);
       ++p ) {
    for (i = 0;i < p->n; ++i) {
      p->ort_ref X(i) = p->ort X(i);
      p->ort_ref Y(i) = p->ort Y(i);
    }
  }

#endif

}


#endif

/******************************************************************************
*
* write_distrib write spatial distribution of potential and kinetic energy
*
******************************************************************************/

void write_distrib(int steps)
{
  FILE *outpot, *outkin, *outminmax;
  str255 fnamepot, fnamekin, fnameminmax;
  size_t size, count_pot, count_kin;
  vektor scale;
  ivektor coord;
  cell *p;
  float *pot;
  float *kin;
  shortint *num;
  float minpot, maxpot, minkin, maxkin;
  int fzhlr,i,j,r,s,t;
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

  size = dist_dim.x * dist_dim.y;
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

  /* create filename */
  /* Dateiname fuer Ausgabedatei erzeugen */
  fzhlr = steps / dis_interval;

  sprintf(fnamepot,"%s.%u.dist",outfilename,fzhlr);
  sprintf(fnamekin,"%s.%u.kin.dist",outfilename,fzhlr);
  sprintf(fnameminmax,"%s.minmax.dist",outfilename);

  /* the dist bins are orthogonal boxes in space */
  scale = box_x; 
  if (scale.x < box_y.x) scale.x = box_y.x; 
  if (scale.y < box_y.y) scale.y = box_y.y; 

  scale.x = dist_dim.x / scale.x;
  scale.y = dist_dim.y / scale.y;

  /* loop over all atoms */
  for ( r = cellmin.x; r < cellmax.x; ++r )
    for ( s = cellmin.y; s < cellmax.y; ++s ) {
      p = PTR_2D_V(cell_array, r, s, cell_dim);
      for (i = 0;i < p->n; ++i) {
        coord.x = (int) (p->ort X(i) * scale.x);
        coord.y = (int) (p->ort Y(i) * scale.y);
        /* Check bounds */
        if (coord.x<0          ) coord.x = 0;
        if (coord.x>=dist_dim.x) coord.x = dist_dim.x-1;
        if (coord.y<0          ) coord.y = 0;
        if (coord.y>=dist_dim.y) coord.y = dist_dim.y-1;
        /* Add up distribution */
        pot = PTR_2D_VV(pot_hist_local, coord, dist_dim);
        kin = PTR_2D_VV(kin_hist_local, coord, dist_dim);
        num = PTR_2D_VV(num_hist_local, coord, dist_dim);
        (*num)++;
#ifdef DISLOC
        if (Epot_diff==1) {
          *pot += p->pot_eng[i] - p->Epot_ref[i];
        } else
#endif
        *pot += p->pot_eng[i];
        *kin += SPRODN(p->impuls,i,p->impuls,i) / (2*p->masse[i]);
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
    outpot = fopen(fnamepot,"w");
    if (NULL == outpot) error("Cannot open pot distrib file.");
    outkin = fopen(fnamekin,"w");
    if (NULL == outkin) error("Cannot open kin distrib file.");

    for (i=0; i<size; i++) {
      if (num_hist[i]>0) {
         pot_hist[i] /= num_hist[i];
         kin_hist[i] /= num_hist[i];
      }
    }

    j=0;
    while (num_hist[j]==0) j++;
    minpot = pot_hist[j];
    maxpot = pot_hist[j];
    minkin = kin_hist[j];
    maxkin = kin_hist[j];
    for (i=j+1; i<size; i++) {
      if (num_hist[i]>0) {
        if (maxpot<pot_hist[i]) maxpot=pot_hist[i];
        if (minpot>pot_hist[i]) minpot=pot_hist[i];
        if (maxkin<kin_hist[i]) maxkin=kin_hist[i];
        if (minkin>kin_hist[i]) minkin=kin_hist[i];
      }
    }

    outminmax = fopen(fnameminmax, "a");
    fprintf(outminmax, "%d %f %f %f %f\n", 
            fzhlr, minpot, maxpot, minkin, maxkin);
    fclose(outminmax);

    if (dist_binary_io) {
      count_pot=fwrite(pot_hist, sizeof(float), size, outpot);
      count_kin=fwrite(kin_hist, sizeof(float), size, outkin);
      if ((count_pot!=size) || (count_kin!=size)) {
        fprintf(stderr,"dist write incomplete - cnt_pot = %d, cnt_kin = %d\n", 
                count_pot, count_kin );
      }
    } else {
      for ( r = 0; r < dist_dim.x; ++r )
        for ( s = 0; s < dist_dim.y; ++s ) {
          pot = PTR_2D_V(pot_hist, r, s, dist_dim);
          kin = PTR_2D_V(kin_hist, r, s, dist_dim);
          fprintf(outpot,"%f\n", *pot);
          fprintf(outkin,"%f\n", *kin);
        };
    }

    fclose(outpot);
    fclose(outkin);
  }
}
