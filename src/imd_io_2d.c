
/******************************************************************************
*
* imd_io_2d.c -- 2D-specific IO routines
*
******************************************************************************/

/******************************************************************************
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
  real refeng=0;
  real fdummy;
  int refn, idummy;
  vektor2d refpos;
#endif
  real m;
  int i,s,n;
  cell *to;
  ivektor2d cellc;
  int to_cpu;
  int addnumber = 0;
#ifdef MPI
  real a[MAX_ATOM_SIZE];
  msgbuf mpibuf;
#endif

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
  if (calc_Epot_ref == 0) reffile = fopen(reffilename,"r");
#endif

  mpibuf.data = a;

#else /* not MPI */

  infile = fopen(infilename,"r");
#ifdef DISLOC
  if (calc_Epot_ref == 0) reffile = fopen(reffilename,"r");
#endif

#endif /* MPI */

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

  natoms  = 0;
  nactive = 0;
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

#ifdef DOUBLE
    p = sscanf(buf,"%d %d %lf %lf %lf %lf %lf",
	      &n,&s,&m,&pos.x,&pos.y,&vau.x,&vau.y);
#else
    p = sscanf(buf,"%d %d %f %f %f %f %f",
	      &n,&s,&m,&pos.x,&pos.y,&vau.x,&vau.y);  
#endif

    if (p>0) {  /* skip empty lines at end of file */

#ifdef DISLOC
      if (calc_Epot_ref == 0) {
#ifdef DOUBLE
        pref = sscanf(refbuf,"%d %d %lf %lf %lf %lf %lf %lf",
	 	      &refn,&idummy,&fdummy,&refpos.x,&refpos.y,
                      &fdummy,&fdummy,&refeng);
#else
        pref = sscanf(refbuf,"%d %d %f %f %f %f %f %f",
		      &refn,&idummy,&fdummy,&refpos.x,&refpos.y,
                      &fdummy,&fdummy,&refeng);
#endif
        if (ABS(refn) != ABS(n)) 
          error("Numbers in infile and reffile are different.\n");
      }
#endif /* DISLOC */

      pos = back_into_box(pos);

      if (0>=m) error("Mass zero or negative.\n");
      if ((p!=5) && (p<7)) error("incorrect line in configuration file.");

      input->n = 1;
      input->nummer[0] = n;
      input->sorte[0]  = s;
      input->masse[0]  = m;
      input->ort X(0)  = pos.x;
      input->ort Y(0)  = pos.y;
      if (p==5) {
        do_maxwell=1;
        input->impuls X(0) = 0;
        input->impuls Y(0) = 0;
      } else {
	input->impuls X(0) = vau.x * m * (restrictions+s)->x;
	input->impuls Y(0) = vau.y * m * (restrictions+s)->y;
      }
      input->kraft X(0) = 0;
      input->kraft Y(0) = 0;
#ifdef DISLOC
      input->ort_ref X(0) = refpos.x;
      input->ort_ref Y(0) = refpos.y;
      input->Epot_ref[0]  = refeng;
#endif

      cellc = cell_coord(pos.x,pos.y);

#ifdef MPI

      to_cpu = cpu_coord(cellc);
      if ((myid != to_cpu) && (1!=parallel_input)) {
        natoms++;
        /* we still have s == input->sorte[0] */
        if (s < ntypes)
          nactive += DIM;
        else {
          nactive += (int) (restrictions+s)->x;
          nactive += (int) (restrictions+s)->y;
        }
        num_sort[SORTE(input,0)]++;
        mpibuf.n = 0;
        copy_one_atom(&mpibuf, input, 0, 0);
        MPI_Send(mpibuf.data, mpibuf.n, REAL, to_cpu, CELL_TAG, cpugrid);
      } else if (to_cpu==myid) {  
        natoms++;
        /* we still have s == input->sorte[0] */
        if (s < ntypes)
          nactive += DIM;
        else {
          nactive += (int) (restrictions+s)->x;
          nactive += (int) (restrictions+s)->y;
        }
        num_sort[SORTE(input,0)]++;
	cellc = local_cell_coord(pos.x,pos.y);
        to = PTR_VV(cell_array,cellc,cell_dim);
	move_atom(to, input, 0);
      }

#else /* not MPI */

      natoms++;
      /* we still have s == input->sorte[0] */
      if (s < ntypes)
        nactive += DIM;
      else {
        nactive += (int) (restrictions+s)->x;
        nactive += (int) (restrictions+s)->y;
      }
      num_sort[SORTE(input,0)]++;
      to = PTR_VV(cell_array,cellc,cell_dim);
      move_atom(to, input, 0);

#endif /* MPI */

    } /* (p>0) */
  } /* !feof(infile) */

  fclose(infile);  
#ifdef DISLOC
  if (calc_Epot_ref == 0) fclose(reffile);
#endif

#ifdef MPI

  /* Tell other CPUs that reading atoms is finished. (tag==0) */
  if (1!=parallel_input)
    for (s=1; s<num_cpus; s++)
      MPI_Send(mpibuf.data, 1, REAL, s, 0, cpugrid);

  /* Add the number of atoms read (and kept) by each CPU */
  if (1==parallel_input) {
    MPI_Allreduce( &natoms,  &addnumber, 1, MPI_INT, MPI_SUM, cpugrid);
    natoms = addnumber;
    MPI_Allreduce( &nactive, &addnumber, 1, MPI_INT, MPI_SUM, cpugrid);
    nactive = addnumber;
    for (i=0; i<ntypes; i++) {
      MPI_Allreduce(&num_sort[i], &addnumber, 1, MPI_INT, MPI_SUM, cpugrid);
      num_sort[i]=addnumber;
    }
  } else { /* broadcast */
    MPI_Bcast( &natoms ,      1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( &nactive,      1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( num_sort, ntypes, MPI_INT, 0, MPI_COMM_WORLD);
  }

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
    }
    printf(" ],  total = %u\n",addnumber);
  }
}


#ifdef MPI

/******************************************************************************
*
*  recveive atoms one at a time from CPU 0
*  this is only used when parallel_input==0
*
******************************************************************************/

void recv_atoms(void)
{
  MPI_Status status;
  real a[MAX_ATOM_SIZE];
  msgbuf mpibuf;   
  mpibuf.data = a;

  printf("Node %d listening.\n",myid);
  while ( 1 ) {
    MPI_Recv(mpibuf.data, 100, REAL, 0, MPI_ANY_TAG, cpugrid, &status);
    if (0 == status.MPI_TAG) break;
    process_buffer( &mpibuf, (cell *) NULL );
  }
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
  part_pot_energy =       tot_pot_energy / natoms;
  part_kin_energy = 2.0 * tot_kin_energy / nactive;
#endif
  vol = volume / natoms;

  out = fopen(fname,"a");
  if (NULL == out) error("Cannot open properties file.");

  fprintf(out,"%e",  (double)(steps * timestep));
  fprintf(out," %e", (double)part_pot_energy);
#ifndef MC
  fprintf(out," %e", (double)part_kin_energy);
#ifdef FNORM
  fprintf(out, "%e", (double)fnorm);
#endif
  fprintf(out," %e", (double)pressure);
#else
  fprintf(out," %e", (double)(mc_accept/(real)mc_count));
  mc_accept = (real)0;
  mc_count  = 0;
#endif
  fprintf(out," %e", (double)vol);
  if (ensemble==ENS_NPT_AXIAL) {
    fprintf(out," %e %e", 
                  (double) stress.x, (double) stress.y );
    fprintf(out," %e %e", 
                  (double) box_x.x,  (double) box_y.y  );
  }
#if defined(NVT) || defined(NPT) || defined(STM)
  fprintf(out," %e", eta );
#endif

  putc('\n',out);

  fclose(out);

}


/******************************************************************************
*
*  write_cell - utility routine for write_config
*
******************************************************************************/

void write_cell(FILE *out, cell *p)
{
  int i;
  double h;

  for (i=0; i<p->n; i++) {
#ifdef NVX
      h  = SPRODN(p->impuls,i,p->impuls,i)/(2*p->masse[i])+p->heatcond[i];
      h *=  p->impuls X(i) /  p->masse[i];
      fprintf(out,"%d %d %12f %12f %12f %12f %12f %12f %12f\n",
#else
      fprintf(out,"%d %d %12f %12f %12f %12f %12f %12f\n",
#endif
        p->nummer[i],
        p->sorte[i],
        p->masse[i],
        p->ort X(i),
        p->ort Y(i),
        p->impuls X(i) / p->masse[i],
        p->impuls Y(i) / p->masse[i],
        p->pot_eng[i]
#ifdef NVX
        ,h
#endif
      );
  }
}

/******************************************************************************
*
* write_config writes a configuration to a numbered file,
* which can serve as a checkpoint; uses write_cell
*
******************************************************************************/

void write_config(int steps)
{ 
  FILE *out;
  str255 fname;
  int fzhlr;

  fzhlr = steps / rep_interval;

  /* write checkpoint */
  write_config_select(fzhlr,"chkpt",write_cell);

  /* write iteration file */
  if (myid == 0) {
    sprintf(fname,"%s.%u.itr",outfilename,fzhlr);

    out = fopen(fname,"w");
    if (NULL == out) error("Cannot write iteration file.");

    fprintf(out,"# checkpoint %d\n",fzhlr);
    fprintf(out,"startstep \t%d\n",steps+1);
    fprintf(out,"box_x \t%f %f\n",box_x.x,box_x.y);
    fprintf(out,"box_y \t%f %f\n",box_y.x,box_y.y);
    fprintf(out,"starttemp \t%f\n",temperature);
#if defined(NVT) || defined(NPT) || defined(STM) 
    fprintf(out,"eta \t%f\n",eta);
#endif

#ifdef FBC
    for(m=0; m<vtypes;m++)
      fprintf(out,"extra_startforce %d %.21g %.21g \n",
	      m,(fbc_forces+m)->x,(fbc_forces+m)->y);
#endif

#ifdef NPT
    if (ensemble==ENS_NPT_ISO) {
      fprintf(out,"pressure_ext \t%f\n",pressure_ext.x);
      fprintf(out,"xi \t%f\n",xi.x);
    }
    if (ensemble==ENS_NPT_AXIAL) {
      fprintf(out,"pressure_ext \t%f %f\n",pressure_ext.x,pressure_ext.y);
      fprintf(out,"xi \t%f %f\n",xi.x,xi.y);
    }
#endif
    fclose(out);
  }
}
