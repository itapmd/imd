/******************************************************************************
*
* imd_io.c -- IO routines for the imd package
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
  vektor pos;
  vektor vau;
#ifdef DISLOC
  FILE *reffile;
  int pref;
  char refbuf[512];
  real refeng=0;
  real fdummy;
  int refn, idummy;
  vektor3d refpos;
#endif
  real m;
#ifdef UNIAX
  real the;
  vektor axe;
  vektor sig;
  vektor eps;
  vektor ome;
#endif
  int i,s,n;
  cell *to;
  ivektor cellc;
  int to_cpu;
  int addnumber = 0;

  /* allocate num_sort on all CPUs */
  if ((num_sort=calloc(ntypes,sizeof(int)))==NULL)
    error("cannot allocate memory for num_sort\n");

#ifdef MPI

  /* Try opening a per cpu file first when parallel_input is active */
  if (1==parallel_input) {
    sprintf(buf,"%s.%u",infilename,myid); 
    infile = fopen(buf,"r");
    /* When each cpu reads only part of the atoms, we have to add the
       number of atoms together to get the correct natoms. We set a
       flag here */
    if (NULL!=infile) addnumber=1;
  } else if (0!=myid) {
    recv_atoms();
    /* If CPU 0 found velocities in its data, no initialisation is done */
    MPI_Bcast( &natoms,       1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( &nactive,      1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( num_sort, ntypes, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( &do_maxwell,   1, MPI_INT, 0, MPI_COMM_WORLD);
    return;
  }

  if ((1!=parallel_input) || (NULL==infile)) 
    infile = fopen(infilename,"r");
#ifdef DISLOC
  if ((calc_Epot_ref == 0) || (NULL==reffile))
    reffile = fopen(reffilename,"r");
#endif

#else /* not MPI */

  infile = fopen(infilename,"r");
#ifdef DISLOC
  if (calc_Epot_ref == 0) 
    reffile = fopen(reffilename,"r");
#endif

#endif /* MPI or not MPI */

  if (NULL==infile) error("Cannot open atoms file.");
#ifdef DISLOC
  if ((calc_Epot_ref == 0) && (NULL==reffile)) 
    error("Cannot open reference file.");
#endif

 /* Set up 1 atom input cell */
  input = (cell *) malloc(sizeof(cell));
  if (0==input) error("Cannot allocate input cell.");
  input->n_max=0;
  alloc_cell(input, 1);

  natoms=0;
  nactive=0;
#ifdef SHOCK
  do_maxwell=1;
#else
  do_maxwell=0;
#endif

  /* Read the input file line by line */
  while(!feof(infile)) {

    buf[0] = (char) NULL;
    fgets(buf,sizeof(buf),infile);
    /* eat comments */
    while ('#'==buf[1]) fgets(buf,sizeof(buf),infile); 

#ifdef DISLOC
    if (calc_Epot_ref == 0) {
      refbuf[0] = (char) NULL;
      fgets(refbuf,sizeof(refbuf),reffile);
      /* eat comments */
      while ('#'==refbuf[1]) fgets(refbuf,sizeof(refbuf),reffile);
    }
#endif

#ifdef UNIAX 

#ifdef DOUBLE
    p = sscanf(buf,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	      &n,&s,&m,&the,&pos.x,&pos.y,&pos.z,&axe.x,&axe.y,&axe.z,&sig.x,&sig.y,&sig.z,&eps.x,&eps.y,&eps.z,&vau.x,&vau.y,&vau.z,&ome.x,&ome.y,&ome.z);  
    if ( axe.x * axe.x + axe.y * axe.y + axe.z * axe.z 
	 - 1.0 > 1.0e-04 )
      error("Molecular axis not a unit vector!");
    if ( sig.x != sig.y )
      error("This is not a uniaxial molecule!");
    if ( eps.x != eps.y )
      error("This is not a uniaxial molecule!");
#else
    p = sscanf(buf,"%d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
	      &n,&s,&m,&the,&pos.x,&pos.y,&pos.z,&axe.x,&axe.y,&axe.z,&sig.x,&sig.y,&sig.z,&eps.x,&eps.y,&eps.z,&vau.x,&vau.y,&vau.z,&ome.x,&ome.y,&ome.z);  
    if ( axe.x * axe.x + axe.y * axe.y + axe.z * axe.z 
	 - 1.0 > 1.0e-04 )
      error("Molecular axis not a unit vector!");
    if ( sig.x != sig.y )
      error("This is not a uniaxial molecule!");
    if ( eps.x != eps.y )
      error("This is not a uniaxial molecule!");
#endif

#else /* not UNIAX */

#ifdef DOUBLE
    p = sscanf(buf,"%d %d %lf %lf %lf %lf %lf %lf %lf",
	      &n,&s,&m,&pos.x,&pos.y,&pos.z,&vau.x,&vau.y,&vau.z);  
#else
    p = sscanf(buf,"%d %d %f %f %f %f %f %f %f",
	      &n,&s,&m,&pos.x,&pos.y,&pos.z,&vau.x,&vau.y,&vau.z);
#endif

#endif /* UNIAX or not UNIAX */

    if (p>0) {  /* skip empty lines at end of file */

#ifdef DISLOC
      if (calc_Epot_ref == 0) {
#ifdef DOUBLE
        pref = sscanf(refbuf,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf",
	  	    &refn,&idummy,&fdummy,&refpos.x,&refpos.y,&refpos.z,
                    &fdummy,&fdummy,&fdummy,&refeng);
#else
        pref = sscanf(refbuf,"%d %d %f %f %f %f %f %f %f %f",
		    &refn,&idummy,&fdummy,&refpos.x,&refpos.y,&refpos.z,
                    &fdummy,&fdummy,&fdummy,&refeng);
#endif
        if ((ABS(refn)) != (ABS(n))) {
          printf("%d %d\n", n, refn); 
          error("Numbers in infile and reffile are different.\n");
        }
      }
#endif /* DISLOC */

      pos = back_into_box(pos);

      if (0>=m) error("Mass zero or negative.\n");
#ifdef UNIAX
      if ((p!=16) && (p<22)) error("incorrect line in configuration file.");
#else
      if ((p!=6) && (p<9)) error("incorrect line in configuration file.");
#endif

      input->n = 1;
#ifndef MONOLJ
      input->nummer[0] = n;
      input->sorte[0]  = s;
      input->masse[0]  = m;
#endif
      input->ort X(0) = pos.x;
      input->ort Y(0) = pos.y;
      input->ort Z(0) = pos.z;
#ifdef UNIAX
      input->traeg_moment[0] = the;
      input->achse X(0) = axe.x ;
      input->achse Y(0) = axe.y ;
      input->achse Z(0) = axe.z ;
      input->shape X(0) = sig.x ;
      input->shape Y(0) = sig.y ;
      input->shape Z(0) = sig.z ;
      input->pot_well X(0) = eps.x ;
      input->pot_well Y(0) = eps.y ;
      input->pot_well Z(0) = eps.z ;
      if (p==16) {
	do_maxwell=1;
	input->impuls X(0) = 0 ;
	input->impuls Y(0) = 0 ;
	input->impuls Z(0) = 0 ;
	input->dreh_impuls X(0) = 0 ;
	input->dreh_impuls Y(0) = 0 ;
	input->dreh_impuls Z(0) = 0 ;
      } else {
	input->impuls X(0) = vau.x * m;
	input->impuls Y(0) = vau.y * m;
	input->impuls Z(0) = vau.z * m;
	input->dreh_impuls X(0) = ome.x * the;
	input->dreh_impuls Y(0) = ome.y * the;
	input->dreh_impuls Z(0) = ome.z * the;
      }
      input->dreh_moment X(0) = 0 ;
      input->dreh_moment Y(0) = 0 ;
      input->dreh_moment Z(0) = 0 ;
#else /* not UNIAX */
      if (p==6) {
        do_maxwell=1;
        input->impuls X(0) = 0;
        input->impuls Y(0) = 0;
        input->impuls Z(0) = 0;
      } else {
        input->impuls X(0) = vau.x * m;
        input->impuls Y(0) = vau.y * m;
        input->impuls Z(0) = vau.z * m;
      }
#endif /* UNIAX or not UNIAX */
      input->kraft  X(0) = 0;
      input->kraft  Y(0) = 0;
      input->kraft  Z(0) = 0;
#ifdef DISLOC
      input->ort_ref X(0) = refpos.x;
      input->ort_ref Y(0) = refpos.y;
      input->ort_ref Z(0) = refpos.z;
      input->Epot_ref[0]  = refeng;
#endif

      cellc = cell_coord(pos.x,pos.y,pos.z);

#ifdef MPI

      to_cpu = cpu_coord(cellc);
      if ((myid != to_cpu) && (1!=parallel_input)) {
        natoms++;
        if (NUMMER(input,0)>=0) nactive++;
        num_sort[SORTE(input,0)]++;
	MPI_Ssend( input->ort,     DIM, MPI_REAL,to_cpu,ORT_TAG,    cpugrid);
#ifdef DISLOC
	MPI_Ssend( input->ort_ref, DIM, MPI_REAL,to_cpu,ORT_REF_TAG,cpugrid);
	MPI_Ssend( input->Epot_ref,1,MPI_REAL,to_cpu,POT_REF_TAG,cpugrid);
#endif
#ifndef MONOLJ
	MPI_Ssend( input->sorte,  1, SHORT    , to_cpu, SORTE_TAG , cpugrid);
	MPI_Ssend( input->masse,  1, MPI_REAL , to_cpu, MASSE_TAG , cpugrid);
        MPI_Ssend( input->nummer, 1, INTEGER  , to_cpu, NUMMER_TAG, cpugrid);
#ifdef UNIAX
	MPI_Ssend( input->traeg_moment,  1, MPI_REAL , to_cpu, TRAEG_MOMENT_TAG , cpugrid);
	MPI_Ssend( input->achse,   DIM, MPI_REAL,to_cpu, ACHSE_TAG, cpugrid);
	MPI_Ssend( input->shape, DIM, MPI_REAL,to_cpu, SHAPE_TAG, cpugrid);
	MPI_Ssend( input->pot_well, DIM, MPI_REAL,to_cpu, POT_WELL_TAG, cpugrid);
	MPI_Ssend( input->dreh_impuls, DIM, MPI_REAL,to_cpu, DREH_IMPULS_TAG, cpugrid);
#endif
#endif
	MPI_Ssend( input->impuls,DIM, MPI_REAL, to_cpu, IMPULS_TAG, cpugrid);
      } else if (to_cpu==myid) {
        natoms++;  
        if (NUMMER(input,0)>=0) nactive++;
        num_sort[SORTE(input,0)]++;
	cellc = local_cell_coord(pos.x,pos.y,pos.z);
	move_atom(cellc, input, 0);
      }

#else /* not MPI */

      natoms++;  
      if (NUMMER(input,0)>=0) nactive++;
      num_sort[SORTE(input,0)]++;
      move_atom(cellc, input, 0);

#endif /* MPI or not MPI */

    } /* (p>0) */
  } /* !feof(infile) */

  fclose(infile);  
#ifdef DISLOC
  if (calc_Epot_ref==0) fclose(reffile);
#endif

#ifdef MPI

/* Tell other CPUs that reading atoms is finished. (tag==0) */
  if (1!=parallel_input)
    for (s=1; s<num_cpus; s++)
      MPI_Ssend( input->ort, 3, MPI_REAL, s, 0, cpugrid);

  /* Add the number of atoms read (and kept) by each CPU */
  if (1==parallel_input) {
    MPI_Allreduce( &natoms,  &addnumber, 1, MPI_INT, MPI_SUM, cpugrid);
    natoms = addnumber;
    MPI_Allreduce( &nactive, &addnumber, 1, MPI_INT, MPI_SUM, cpugrid);
    nactive = addnumber;
    for (i=0; i<ntypes; i++) {
      MPI_Allreduce( &num_sort[i], &addnumber, 1, MPI_INT, MPI_SUM, cpugrid);
      num_sort[i] = addnumber;
    };
  } else { /* broadcast if serial io */
    MPI_Bcast( &natoms ,      1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( &nactive,      1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( num_sort, ntypes, MPI_INT, 0, MPI_COMM_WORLD);
  };

  /* If CPU 0 found velocities in its data, no initialisation is done */
  MPI_Bcast( &do_maxwell , 1, MPI_INT, 0, MPI_COMM_WORLD);

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
  }
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
  ivektor cellc;
  ivektor local_cellc;

  printf("Node %d listening.\n",myid);

  input = (cell *) malloc(sizeof(cell));
  if (0==input) error("Cannot allocate input cell.");
  input->n_max = 0;
  alloc_cell(input,1);

  while ( 1 ) {

    MPI_Recv(input->ort,  DIM, MPI_REAL, 0, MPI_ANY_TAG   , cpugrid, &status );

    if ((0 != status.MPI_TAG) && (ORT_TAG != status.MPI_TAG)) 
       error("Messages mixed up.");

    if ( 0 == status.MPI_TAG ) break;

#ifndef MONOLJ
    MPI_Recv(input->sorte,  1, SHORT,     0, SORTE_TAG , cpugrid, &status );
    MPI_Recv(input->masse,  1, MPI_REAL,  0, MASSE_TAG , cpugrid, &status );
    MPI_Recv(input->nummer, 1, INTEGER,   0, NUMMER_TAG, cpugrid, &status );
#ifdef UNIAX
    MPI_Recv(input->traeg_moment,  1, MPI_REAL,  0, TRAEG_MOMENT_TAG , cpugrid, &status );
    MPI_Recv(input->achse, DIM, MPI_REAL, 0, ACHSE_TAG , cpugrid, &status );
    MPI_Recv(input->shape, DIM, MPI_REAL, 0, SHAPE_TAG , cpugrid, &status );
    MPI_Recv(input->pot_well, DIM, MPI_REAL, 0, POT_WELL_TAG , cpugrid, &status );
    MPI_Recv(input->dreh_impuls,DIM, MPI_REAL, 0, DREH_IMPULS_TAG, cpugrid, &status );
#endif
#endif
    MPI_Recv(input->impuls,DIM, MPI_REAL, 0, IMPULS_TAG, cpugrid, &status );
#ifdef DISLOC
    MPI_Recv(input->Epot_ref,1, MPI_REAL, 0, POT_REF_TAG, cpugrid, &status);
    MPI_Recv(input->ort_ref, DIM, MPI_REAL, 0, ORT_REF_TAG,  cpugrid, &status);
#endif

    local_cellc = local_cell_coord(input->ort X(0),input->ort Y(0),
                                   input->ort Z(0));
    target = PTR_3D_VV(cell_array,local_cellc,cell_dim);

    /* See if we need some space */
    if (target->n >= target->n_max) alloc_cell(target,target->n_max+incrsz);

    target->ort X(target->n)    = input->ort X(0);
    target->ort Y(target->n)    = input->ort Y(0);
    target->ort Z(target->n)    = input->ort Z(0);
#ifndef MONOLJ
    target->sorte[target->n]  = input->sorte[0];
    target->masse[target->n]  = input->masse[0];
    target->nummer[target->n] = input->nummer[0];
#ifdef UNIAX
    target->traeg_moment[target->n]  = input->traeg_moment[0];
    target->achse X(target->n)  = input->achse X(0);
    target->achse Y(target->n)  = input->achse Y(0);
    target->achse Z(target->n)  = input->achse Z(0);
    target->shape X(target->n)  = input->shape X(0);
    target->shape Y(target->n)  = input->shape Y(0);
    target->shape Z(target->n)  = input->shape Z(0);
    target->pot_well X(target->n)  = input->pot_well X(0);
    target->pot_well Y(target->n)  = input->pot_well Y(0);
    target->pot_well Z(target->n)  = input->pot_well Z(0);
    target->dreh_impuls X(target->n) = input->dreh_impuls X(0);
    target->dreh_impuls Y(target->n) = input->dreh_impuls Y(0);
    target->dreh_impuls Z(target->n) = input->dreh_impuls Z(0);
#endif
    target->impuls X(target->n) = input->impuls X(0);
    target->impuls Y(target->n) = input->impuls Y(0);
    target->impuls Z(target->n) = input->impuls Z(0);
#endif /* not MONOLJ */
#ifdef DISLOC
    target->Epot_ref[target->n]  = input->Epot_ref[0];
    target->ort_ref X(target->n) = input->ort_ref X(0);
    target->ort_ref Y(target->n) = input->ort_ref Y(0);
    target->ort_ref Z(target->n) = input->ort_ref Z(0);
#endif
    ++target->n;

  }
  printf("Node %d leaves listen.\n",myid);
}

#endif /* MPI */ 


/******************************************************************************
*
* write_properties writes selected properties to *.eng file
*
******************************************************************************/


void write_properties(int steps)
{
  FILE *out;
  str255 fname;
#ifdef EAM2
  char *format=" %.16e";
#else
  char *format=" %e";
#endif
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
  part_pot_energy =                tot_pot_energy / nactive;
#ifdef UNIAX
  part_kin_energy = ( 2.0 / 5.0 )* tot_kin_energy / nactive;
#else
  part_kin_energy = ( 2.0 / 3.0 )* tot_kin_energy / nactive;
#endif
#endif
  vol = volume / natoms;

  out = fopen(fname,"a");
  if (NULL == out) error("Cannot open properties file.");

  fprintf(out, "%e",   (double)(steps * timestep));
  fprintf(out, format, (double)part_pot_energy);
 
#ifndef MC
  fprintf(out, format, (double)part_kin_energy);
 
  fprintf(out," %e",   (double)pressure);
#else
  fprintf(out," %e",   (double)(mc_accept/(real)mc_count));
  mc_accept = (real)0;
  mc_count  = 0;
#endif
  fprintf(out," %e",   (double)vol);

#ifdef PAXTEST
  if (ensemble==ENS_NPT_AXIAL) {
    fprintf(out," %e %e %e", 
                  (double) stress.x, (double) stress.y, (double) stress.z );
    fprintf(out," %e %e %e", 
                  (double) box_x.x,  (double) box_y.y,  (double) box_z.z );
  }
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
#ifdef UNIAX
#define WRITE_CELL     for (i = 0;i < p->n; ++i) \
             fprintf(out,"%d %d %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f\n",\
	     p->nummer[i],\
	     p->sorte[i],\
	     p->masse[i],\
	     p->traeg_moment[i],\
	     p->ort X(i),\
	     p->ort Y(i),\
	     p->ort Z(i),\
	     p->achse X(i),\
	     p->achse Y(i),\
	     p->achse Z(i),\
	     p->shape X(i),\
	     p->shape Y(i),\
	     p->shape Z(i),\
	     p->pot_well X(i),\
	     p->pot_well Y(i),\
	     p->pot_well Z(i),\
	     p->impuls X(i) / p->masse[i],\
	     p->impuls Y(i) / p->masse[i],\
	     p->impuls Z(i) / p->masse[i],\
	     p->dreh_impuls X(i) / p->traeg_moment[i],\
	     p->dreh_impuls Y(i) / p->traeg_moment[i],\
	     p->dreh_impuls Z(i) / p->traeg_moment[i]) 

#else /* not UNIAX */
#define WRITE_CELL     for (i = 0;i < p->n; ++i) \
             fprintf(out,"%d %d %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f %12.16f\n",\
	     NUMMER(p,i),\
	     p->sorte[i],\
	     MASSE(p,i),\
	     p->ort X(i),\
	     p->ort Y(i),\
	     p->ort Z(i),\
	     p->impuls X(i) / MASSE(p,i),\
	     p->impuls Y(i) / MASSE(p,i),\
	     p->impuls Z(i) / MASSE(p,i),\
             POTENG(p,i)) 
#endif /* UNIAX or not UNIAX */

{ 
  FILE *out;
  str255 fname;
  int fzhlr;
  cell *p,*q;
  int i,j,k,l,m,n,tag;

  /* Dateiname fuer Ausgabedatei erzeugen */
  fzhlr = steps / rep_interval;

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
      for (k = 1; k < cell_dim.y-1; ++k )
	for (l = 1; l < cell_dim.z-1; ++l ) {
 	  p = PTR_3D_V(cell_array, j, k, l, cell_dim);
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
	for (k = 1; k < cell_dim.y-1; ++k )
	  for (l = 1; l < cell_dim.z-1; ++l ) {
	    p = PTR_3D_V(cell_array, j, k, l, cell_dim);
	    WRITE_CELL;
	  };

      /* Receive data from other cpus and write that */
      p   = PTR_3D_V(cell_array, 0, 0, 0, cell_dim);
      for ( m = 1; m < num_cpus; ++m)
	for (j = 1; j < cell_dim.x-1; ++j )
	  for (k = 1; k < cell_dim.y-1; ++k )
	    for (l = 1; l < cell_dim.z-1; ++l ) {
#ifndef MONOLJ
	      tag = PTR_3D_V(CELL_TAG, j, k, l, cell_dim);
	      recv_cell( p, m, tag );
#else
	      recv_cell( p, m, ORT_TAG );
#endif
	      WRITE_CELL;
	    };

      fclose(out);      
    } else { 
      /* Send data to cpu 0 */
      for (j = 1; j < cell_dim.x-1; ++j )
	for (k = 1; k < cell_dim.y-1; ++k )
	  for (l = 1; l < cell_dim.z-1; ++l ) {
	    p   = PTR_3D_V(cell_array, j, k, l, cell_dim);
#ifndef MONOLJ
	    tag = PTR_3D_V(CELL_TAG, j, k, l, cell_dim);
	    send_cell( p, 0, tag );
#else
	    send_cell( p, 0, ORT_TAG );
#endif
	  };
    };
  };
#else

  /* Ausgabedatei oeffnen */
  out = fopen(fname,"w");
  if (NULL == out) error("Cannot open output file for config.");

  for (p = cell_array; 
       p <= PTR_3D_V(cell_array,
                     cell_dim.x-1,
                     cell_dim.y-1,
                     cell_dim.z-1,
                     cell_dim);
       ++p ) 
    
      WRITE_CELL;

  fclose(out);  

#endif

  /* write iteration file */
  if (myid == 0) {
    sprintf(fname,"%s.%u.itr",outfilename,fzhlr);
    out = fopen(fname,"w");
    if (NULL == out) error("Cannot write iteration file.");

    fprintf(out,"startstep \t%d\n",steps+1);
    fprintf(out,"box_x \t%.16f %.16f %.16f\n",box_x.x,box_x.y,box_x.z);
    fprintf(out,"box_y \t%.16f %.16f %.16f\n",box_y.x,box_y.y,box_y.z);
    fprintf(out,"box_z \t%.16f %.16f %.16f\n",box_z.x,box_z.y,box_z.z);
    fprintf(out,"starttemp \t%f\n",temperature);
#if defined(NVT) || defined(NPT)
    fprintf(out,"eta \t%f\n",eta);
#ifdef UNIAX
    fprintf(out,"eta_rot \t%f\n",eta_rot);
#endif
#endif

#ifdef FBC
    for(n=0; n<vtypes;n++)
      fprintf(out,"extra_startforce %d %.21g %.21g %.21g \n",
	      n,(fbc_forces+n)->x,(fbc_forces+n)->y,(fbc_forces+n)->z);
#endif

#ifdef NPT
    if (ensemble==ENS_NPT_ISO) {
      fprintf(out,"pressure_ext \t%f\n",pressure_ext.x);
      fprintf(out,"xi \t%f\n",xi.x);
    }
    if (ensemble==ENS_NPT_AXIAL) {
      fprintf(out,"pressure_ext \t%f %f %f\n",
              pressure_ext.x,pressure_ext.y,pressure_ext.z);
      fprintf(out,"xi \t%f %f %f\n", xi.x,xi.y,xi.z);
    }
#endif /* NPT */
    fclose(out);
  }
}


#ifdef DISLOC

/******************************************************************************
*
* write_demmaps writes a differential energy map to file *.dem.x
*
******************************************************************************/

void write_demmaps(int steps)

#ifdef MONOLJ

#define WRITE_CELL     for (i = 0;i < p->n; ++i) \
             fprintf(out,"%d %d %12f %12f %12f %12f %12f %12f %12f %12f\n",\
	     NUMMER(p,i),\
	     SORTE(p,i),\
	     MASSE(p,i),\
	     p->ort X(i),\
	     p->ort Y(i),\
	     p->ort Z(i),\
	     p->impuls X(i) / MASSE(p,i),\
	     p->impuls Y(i) / MASSE(p,i),\
	     p->impuls Z(i) / MASSE(p,i),\
             POTENG(p,i)) 

#else

#define WRITE_CELL_DEM     for (i = 0;i < p->n; ++i) {\
             if (p->sorte[i] == dpotsorte) {\
               dpot = p->pot_eng[i] - p->Epot_ref[i];\
               if (ABS(dpot) > min_dpot)\
                 fprintf(demout,"%12f %12f %12f %12f\n",\
	         p->ort X(i),\
	         p->ort Y(i),\
	         p->ort Z(i),\
                 dpot);\
	     }\
          }

#endif

{ 
  FILE *demout;
  str255 demfname;
  int fzhlr;
  cell *p,*q;
  int i,j,k,l,m,tag;
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
    for (k = 1; k < cell_dim.y-1; ++k )
      for (l = 1; l < cell_dim.z-1; ++l ) {
	p = PTR_3D_V(cell_array, j, k, l, cell_dim);
	WRITE_CELL_DEM;
      };

  /* Receive data from other cpus and write that */
  p   = PTR_3D_V(cell_array, 0, 0, 0, cell_dim);
    for ( m = 1; m < num_cpus; ++m)
      for (j = 1; j < cell_dim.x-1; ++j )
        for (k = 1; k < cell_dim.y-1; ++k )
          for (l = 1; l < cell_dim.z-1; ++l ) {
#ifndef MONOLJ
            tag = PTR_3D_V(CELL_TAG, j, k, l, cell_dim);
            recv_cell( p, m, tag );
#else
            recv_cell( p, m, POT_REF_TAG );
#endif
	    WRITE_CELL_DEM;
	  };

  fclose(demout);      

  } else { 
  /* Send data to cpu 0 */
  for (j = 1; j < cell_dim.x-1; ++j )
    for (k = 1; k < cell_dim.y-1; ++k )
      for (l = 1; l < cell_dim.z-1; ++l ) {
        p   = PTR_3D_V(cell_array, j, k, l, cell_dim);
#ifndef MONOLJ
        tag = PTR_3D_V(CELL_TAG, j, k, l, cell_dim);
        send_cell( p, 0, tag );
#else
        send_cell( p, 0, POT_REF_TAG );
#endif
      };
  };

#else

  /* Ausgabedatei oeffnen */
  demout = fopen(demfname,"w");
  if (NULL == demout) error("Cannot open output file for dem.");

  for (p = cell_array; 
       p <= PTR_3D_V(cell_array,
		     cell_dim.x-1,
		     cell_dim.y-1,
		     cell_dim.z-1,
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
             fprintf(out,"%d %d %12f %12f %12f %12f %12f %12f %12f %12f\n",\
	     NUMMER(p,i),\
	     SORTE(p,i),\
	     MASSE(p,i),\
	     p->ort X(i),\
	     p->ort Y(i),\
	     p->ort Z(i),\
	     p->impuls X(i) / MASSE(p,i),\
	     p->impuls Y(i) / MASSE(p,i),\
	     p->impuls Z(i) / MASSE(p,i),\
             POTENG(p,i)) 

#else

#define WRITE_CELL_DDM     for (i = 0;i < p->n; ++i) {\
             dx = p->ort X(i) - p->ort_ref X(i);\
             dy = p->ort Y(i) - p->ort_ref Y(i);\
             dz = p->ort Z(i) - p->ort_ref Z(i);\
             if ((dx < boxx)&&(dy < boxy)&&(dz < boxz))\
               fprintf(ddmout,"%12f %12f %12f %12f %12f %12f\n",\
	       p->ort X(i),\
	       p->ort Y(i),\
	       p->ort Z(i),\
               dx,\
               dy,\
               dz);\
          }

#endif

{ 
  FILE *ddmout;
  str255 ddmfname;
  int fzhlr;
  cell *p,*q;
  int i,j,k,l,m,tag;
  real dx, dy, dz, boxx, boxy, boxz;

  /* Dateiname fuer Ausgabedatei erzeugen */
  fzhlr = steps;
  sprintf(ddmfname,"%s.dsp.%u",outfilename,fzhlr);

  /* 1/2 of the boxlength, this is not correct for non-cartesian boxes ! */
  boxx = box_x.x/2;
  boxy = box_y.y/2;
  boxz = box_z.z/2;

#ifdef MPI

  if (0==myid) {

  /* Ausgabedatei oeffnen */
  ddmout = fopen(ddmfname,"w");
  if (NULL == ddmout) error("Cannot open output file for ddm.");

  /* Write data on CPU 0 */

  /* Write own data */
  for (j = 1; j < cell_dim.x-1; ++j )
    for (k = 1; k < cell_dim.y-1; ++k )
      for (l = 1; l < cell_dim.z-1; ++l ) {
	p = PTR_3D_V(cell_array, j, k, l, cell_dim);
	WRITE_CELL_DDM;
      }

  /* Receive data from other cpus and write that */
  p   = PTR_3D_V(cell_array, 0, 0, 0, cell_dim);
    for ( m = 1; m < num_cpus; ++m)
      for (j = 1; j < cell_dim.x-1; ++j )
        for (k = 1; k < cell_dim.y-1; ++k )
          for (l = 1; l < cell_dim.z-1; ++l ) {
#ifndef MONOLJ
            tag = PTR_3D_V(CELL_TAG, j, k, l, cell_dim);
            recv_cell( p, m, tag );
#else
            recv_cell( p, m, ORT_REF_TAG );
#endif
	    WRITE_CELL_DDM;
	  }
  fclose(ddmout);      

  } else { 
  /* Send data to cpu 0 */
  for (j = 1; j < cell_dim.x-1; ++j )
    for (k = 1; k < cell_dim.y-1; ++k )
      for (l = 1; l < cell_dim.z-1; ++l ) {
        p   = PTR_3D_V(cell_array, j, k, l, cell_dim);
#ifndef MONOLJ
        tag = PTR_3D_V(CELL_TAG, j, k, l, cell_dim);
        send_cell( p, 0, tag );
#else
        send_cell( p, 0, ORT_REF_TAG );
#endif
      };
  };

#else

  /* Ausgabedatei oeffnen */
  ddmout = fopen(ddmfname,"w");
  if (NULL == ddmout) error("Cannot open output file for ddm.");

  for (p = cell_array; 
       p <= PTR_3D_V(cell_array,
		     cell_dim.x-1,
		     cell_dim.y-1,
		     cell_dim.z-1,
		     cell_dim);
       ++p ) {
    WRITE_CELL_DDM;

    }

  fclose(ddmout);

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
  int i,j,k,l,m,tag;
  real dx, dy, dz, boxx, boxy, boxz;

#ifdef MPI

  if (0==myid) {

  /* Update own data */
  for (j = 1; j < cell_dim.x-1; ++j )
    for (k = 1; k < cell_dim.y-1; ++k )
      for (l = 1; l < cell_dim.z-1; ++l ) {
	p = PTR_3D_V(cell_array, j, k, l, cell_dim);
	for (i = 0;i < p->n; ++i) {
	  p->ort_ref X(i) = p->ort X(i);
	  p->ort_ref Y(i) = p->ort Y(i);
	  p->ort_ref Z(i) = p->ort Z(i);
	}
      };

  /* Receive data from other cpus and Update that */
  p   = PTR_3D_V(cell_array, 0, 0, 0, cell_dim);
    for ( m = 1; m < num_cpus; ++m)
      for (j = 1; j < cell_dim.x-1; ++j )
        for (k = 1; k < cell_dim.y-1; ++k )
          for (l = 1; l < cell_dim.z-1; ++l ) {
#ifndef MONOLJ
            tag = PTR_3D_V(CELL_TAG, j, k, l, cell_dim);
            recv_cell( p, m, tag );
#else
            recv_cell( p, m, ORT_REF_TAG );
#endif
	    for (i = 0;i < p->n; ++i) {
	      p->ort_ref X(i) = p->ort X(i);
	      p->ort_ref Y(i) = p->ort Y(i);
	      p->ort_ref Z(i) = p->ort Z(i);
	    }
	  };

  } else { 
  /* Send data to cpu 0 */
  for (j = 1; j < cell_dim.x-1; ++j )
    for (k = 1; k < cell_dim.y-1; ++k )
      for (l = 1; l < cell_dim.z-1; ++l ) {
        p   = PTR_3D_V(cell_array, j, k, l, cell_dim);
#ifndef MONOLJ
        tag = PTR_3D_V(CELL_TAG, j, k, l, cell_dim);
        send_cell( p, 0, tag );
#else
        send_cell( p, 0, ORT_REF_TAG );
#endif
      };
  };

#else

  for (p = cell_array; 
       p <= PTR_3D_V(cell_array,
		     cell_dim.x-1,
		     cell_dim.y-1,
		     cell_dim.z-1,
		     cell_dim);
       ++p ) {
    for (i = 0;i < p->n; ++i) {
      p->ort_ref X(i) = p->ort X(i);
      p->ort_ref Y(i) = p->ort Y(i);
      p->ort_ref Z(i) = p->ort Z(i);
    }
  }

#endif

}

#endif /* DISLOC */

/******************************************************************************
*
* write_distrib write spatial distribution of potential and kinetice energy
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

  size = dist_dim.x * dist_dim.y * dist_dim.z;
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

  sprintf(fnamepot,"%s.%u.pot.dist",outfilename,fzhlr);
  sprintf(fnamekin,"%s.%u.kin.dist",outfilename,fzhlr);
  sprintf(fnameminmax,"%s.minmax.dist",outfilename);

  /* the dist bins are orthogonal boxes in space */
  scale = box_x; 
  if (scale.x < box_y.x) scale.x = box_y.x; 
  if (scale.y < box_y.y) scale.y = box_y.y; 
  if (scale.z < box_y.z) scale.z = box_y.z; 
  
  if (scale.x < box_z.x) scale.x = box_z.x; 
  if (scale.y < box_z.y) scale.y = box_z.y; 
  if (scale.z < box_z.z) scale.z = box_z.z; 

  scale.x = dist_dim.x / scale.x;
  scale.y = dist_dim.y / scale.y;
  scale.z = dist_dim.z / scale.z;

  /* loop over all atoms */
  for ( r = cellmin.x; r < cellmax.x; ++r )
    for ( s = cellmin.y; s < cellmax.y; ++s )
      for ( t = cellmin.z; t < cellmax.z; ++t ) {
	p = PTR_3D_V(cell_array, r, s, t, cell_dim);
        for (i = 0;i < p->n; ++i) {
          coord.x = (int) (p->ort X(i) * scale.x);
          coord.y = (int) (p->ort Y(i) * scale.y);
          coord.z = (int) (p->ort Z(i) * scale.z);
	  /* Check bounds */
	  if (coord.x<0          ) coord.x = 0;
          if (coord.x>=dist_dim.x) coord.x = dist_dim.x-1;
          if (coord.y<0          ) coord.y = 0;
          if (coord.y>=dist_dim.y) coord.y = dist_dim.y-1;
          if (coord.z<0          ) coord.z = 0;
          if (coord.z>=dist_dim.z) coord.z = dist_dim.z-1;
          /* Add up distribution */
          pot = PTR_3D_VV(pot_hist_local, coord, dist_dim);
          kin = PTR_3D_VV(kin_hist_local, coord, dist_dim);
          num = PTR_3D_VV(num_hist_local, coord, dist_dim);
          (*num)++;
#ifndef MONOLJ
#ifdef DISLOC
          if (Epot_diff==1) {
            *pot += p->pot_eng[i] - p->Epot_ref[i];
          } else
#endif
          *pot += p->pot_eng[i];
#endif
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
        for ( s = 0; s < dist_dim.y; ++s )
	  for ( t = 0; t < dist_dim.z; ++t ) {
            pot = PTR_3D_V(pot_hist, r, s, t, dist_dim);
            kin = PTR_3D_V(kin_hist, r, s, t, dist_dim);
            fprintf(outpot,"%f\n", *pot);
            fprintf(outkin,"%f\n", *kin);
          };
    }

    fclose(outpot);
    fclose(outkin);
  }
}
