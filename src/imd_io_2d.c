
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2004 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

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
  str255 msg;
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
  long addnumber = 0;
#ifdef MPI
  msgbuf *input_buf, *b;
#endif

  /* allocate num_sort and num_vsort on all CPUs */
  if ((num_sort = (long *) calloc(ntypes,sizeof(long)))==NULL)
    error("cannot allocate memory for num_sort\n");
  if ((num_vsort = (long *) calloc(vtypes,sizeof(long)))==NULL)
    error("cannot allocate memory for num_vsort\n");

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
    MPI_Bcast( &natoms,        1, MPI_LONG, 0, MPI_COMM_WORLD);  
    MPI_Bcast( &nactive,       1, MPI_LONG, 0, MPI_COMM_WORLD);  
    MPI_Bcast( num_sort,  ntypes, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast( num_vsort, vtypes, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast( &do_maxwell,    1, MPI_INT,  0, MPI_COMM_WORLD);
    return;
  }

  if ((0==parallel_input) || (NULL==infile)) infile = fopen(infilename,"r");
#ifdef DISLOC
  if (calc_Epot_ref == 0) reffile = fopen(reffilename,"r");
#endif

  /* allocate temporary input buffer */
  if (0==parallel_input) {
    input_buf = (msgbuf *) malloc( num_cpus * sizeof(msgbuf) );
    if (NULL==input_buf) error("cannot allocate input buffers");
    input_buf[0].data  = (real *) NULL;
    input_buf[0].n_max = 0;
    input_buf[0].n     = 0;
    for (i=1; i<num_cpus; i++) {
      input_buf[i].data  = (real *) malloc( INPUT_BUF_SIZE * sizeof(real) );
      if (NULL==input_buf[i].data) error("cannot allocate input buffer");
      input_buf[i].n_max = INPUT_BUF_SIZE;
      input_buf[i].n     = 0;
    }
  }

#else /* not MPI */

  infile = fopen(infilename,"r");
#ifdef DISLOC
  if (calc_Epot_ref == 0) reffile = fopen(reffilename,"r");
#endif

#endif /* MPI */

  if (NULL==infile) {
    sprintf(msg,"File %s not found",infilename);
    error(msg);
  }
#ifdef DISLOC
  if ((calc_Epot_ref == 0) && (NULL==reffile)) {
    sprintf(msg,"File %s not found",reffilename);
    error(msg);
  }
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

    buf[0] = '\0';
    fgets(buf,sizeof(buf),infile);
    /* eat comments */
    while ('#'==buf[1]) fgets(buf,sizeof(buf),infile); 

#ifdef DISLOC
    if (calc_Epot_ref == 0) {
      refbuf[0] = '\0';
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
      NUMMER(input,0) = n;
      VSORTE(input,0) = s;
      MASSE(input,0)  = m;
      ORT(input,0,X)  = pos.x;
      ORT(input,0,Y)  = pos.y;
      if (p==5) {
        do_maxwell=1;
        IMPULS(input,0,X) = 0;
        IMPULS(input,0,Y) = 0;
      } else {
	IMPULS(input,0,X) = vau.x * m * (restrictions+s)->x;
	IMPULS(input,0,Y) = vau.y * m * (restrictions+s)->y;
      }
      KRAFT(input,0,X) = 0;
      KRAFT(input,0,Y) = 0;
#ifdef DISLOC
      ORT_REF(input,0,X) = refpos.x;
      ORT_REF(input,0,Y) = refpos.y;
      EPOT_REF(input,0)  = refeng;
#endif

#ifdef EPITAX
      /* determine largest atom number of substrate atoms */
      epitax_sub_n = MAX( n, epitax_sub_n );
#endif

      cellc = cell_coord(pos.x,pos.y);

#ifdef MPI

      to_cpu = cpu_coord(cellc);
      if ((myid != to_cpu) && (0==parallel_input)) {
        natoms++;
        /* we still have s == input->sorte[0] */
        if (s < ntypes)
          nactive += DIM;
        else {
          nactive += (long) (restrictions+s)->x;
          nactive += (long) (restrictions+s)->y;
        }
        num_sort[SORTE(input,0)]++;
        num_vsort[VSORTE(input,0)]++;
        b = input_buf + to_cpu;
        copy_atom(b, to_cpu, input, 0);
        if (b->n_max - b->n < MAX_ATOM_SIZE) {
          MPI_Send(b->data, b->n, REAL, to_cpu, INBUF_TAG, cpugrid);
          b->n = 0;
        }
      } else if (to_cpu==myid) {  
        natoms++;
        /* we still have s == input->sorte[0] */
        if (s < ntypes)
          nactive += DIM;
        else {
          nactive += (long) (restrictions+s)->x;
          nactive += (long) (restrictions+s)->y;
        }
        num_sort[SORTE(input,0)]++;
        num_vsort[VSORTE(input,0)]++;
	cellc = local_cell_coord(cellc);
        to = PTR_VV(cell_array,cellc,cell_dim);
	INSERT_ATOM(to, input, 0);
      }

#else /* not MPI */

      natoms++;
      /* we still have s == input->sorte[0] */
      if (s < ntypes)
        nactive += DIM;
      else {
        nactive += (long) (restrictions+s)->x;
        nactive += (long) (restrictions+s)->y;
      }
      num_sort[SORTE(input,0)]++;
      num_vsort[VSORTE(input,0)]++;
      to = PTR_VV(cell_array,cellc,cell_dim);
      INSERT_ATOM(to, input, 0);

#endif /* MPI */

    } /* (p>0) */
  } /* !feof(infile) */

  fclose(infile);  
#ifdef DISLOC
  if (calc_Epot_ref == 0) fclose(reffile);
#endif

#ifdef MPI

  /* The last buffer is sent with a different tag, which tells the
     target CPU that reading is finished; we increase the size by
     one, so that the buffer is sent even if it is empty */
  if (0==parallel_input)
    for (s=1; s<num_cpus; s++) {
      b = input_buf + s;
      b->data[b->n++] = 0.0;
      MPI_Send(b->data, b->n, REAL, s, INBUF_TAG+1, cpugrid);
      free(b->data);
    }

  /* Add the number of atoms read (and kept) by each CPU */
  if (1==parallel_input) {
    MPI_Allreduce( &natoms,  &addnumber, 1, MPI_LONG, MPI_SUM, cpugrid);
    natoms = addnumber;
    MPI_Allreduce( &nactive, &addnumber, 1, MPI_LONG, MPI_SUM, cpugrid);
    nactive = addnumber;
    for (i=0; i<ntypes; i++) {
      MPI_Allreduce(&num_sort[i], &addnumber, 1, MPI_LONG, MPI_SUM, cpugrid);
      num_sort[i]=addnumber;
    }
    for (i=0; i<vtypes; i++) {
      MPI_Allreduce(&num_vsort[i], &addnumber, 1, MPI_LONG, MPI_SUM, cpugrid);
      num_vsort[i]=addnumber;
    }
  } else { /* broadcast */
    MPI_Bcast( &natoms ,       1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast( &nactive,       1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast( num_sort,  ntypes, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast( num_vsort, vtypes, MPI_LONG, 0, MPI_COMM_WORLD);
  }

  /* If CPU 0 found velocities in its data, no initialisation is done */
  MPI_Bcast( &do_maxwell, 1, MPI_INT, 0, MPI_COMM_WORLD);

#endif /* MPI */

  /* print number of atoms */
  if (0==myid) {
    printf("Read structure with %ld atoms.\n",natoms);
    addnumber=num_sort[0];
    printf("num_sort = [ %ld",num_sort[0]);
    for (i=1; i<ntypes; i++) {
      printf(", %ld",num_sort[i]);
      addnumber+=num_sort[i];
    }
    printf(" ],  total = %ld\n",addnumber);
    addnumber=num_vsort[0];
    printf("num_vsort = [ %ld",num_vsort[0]);
    for (i=1; i<vtypes; i++) {
      printf(", %ld",num_vsort[i]);
      addnumber+=num_vsort[i];
    }
    printf(" ],  total = %ld\n",addnumber);
  }
}


#ifdef MPI

/******************************************************************************
*
*  recveive atoms in several chunks from CPU 0
*  this is only used when parallel_input==0
*
******************************************************************************/

void recv_atoms(void)
{
  MPI_Status status;
  int finished=0;
  msgbuf b;   

  b.data  = (real *) malloc( INPUT_BUF_SIZE * sizeof(real) );
  b.n_max = INPUT_BUF_SIZE;
  if (NULL == b.data) error("cannot allocate input receive buffer");

  printf("Node %d listening.\n",myid);
  do {
    MPI_Recv(b.data, INPUT_BUF_SIZE, REAL, 0, MPI_ANY_TAG, cpugrid, &status);
    MPI_Get_count(&status, REAL, &b.n);
    if (status.MPI_TAG==INBUF_TAG+1) { b.n--; finished=1; } /* last buffer */
    process_buffer( &b, (cell *) NULL );
  } while (0==finished);
  printf("Node %d leaves listen.\n",myid);
  free(b.data);
}

#endif 


/******************************************************************************
*
*  filter function for write_config_select
*  writes data of all atoms for checkpoints
*
******************************************************************************/

void write_atoms_config(FILE *out)
{
  int i, k, len=0;
  cell *p;
  double h;

#ifdef HPO
#define RESOL1 " %12.16f"
#define RESOL2 " %12.16f %12.16f"
#else
#define RESOL1 " %f"
#define RESOL2 " %f %f"
#endif

  for (k=0; k<ncells; k++) {
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; i++) {
#ifdef NVX
      h  = SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) /
                                      (2*MASSE(p,i)) + HEATCOND(p,i);
      h *= IMPULS(p,i,X) / MASSE(p,i);
#endif
      len += sprintf(outbuf+len, "%d %d", NUMMER(p,i), VSORTE(p,i) );
      len += sprintf(outbuf+len, RESOL1, MASSE(p,i) );
      len += sprintf(outbuf+len, RESOL2, ORT(p,i,X), ORT(p,i,Y) );
      len += sprintf(outbuf+len, RESOL2, IMPULS(p,i,X) / MASSE(p,i),
                                         IMPULS(p,i,Y) / MASSE(p,i) );
      len += sprintf(outbuf+len, RESOL1, POTENG(p,i) );
#ifdef NVX
      len += sprintf(outbuf+len, RESOL1, h);
#endif
#ifdef ORDPAR
      len += sprintf(outbuf+len," %d", NBANZ(p,i));
#endif
      len += sprintf(outbuf+len,"\n");
      /* flush or send outbuf if it is full */
      if (len > OUTPUT_BUF_SIZE - 256) flush_outbuf(out,&len,OUTBUF_TAG);
    }
  }
  flush_outbuf(out,&len,OUTBUF_TAG+1);
}

/******************************************************************************
*
*  write iteration file
*
******************************************************************************/

void write_itr_file(int fzhlr, int steps,char *suffix)
{
  FILE *out;
  str255 fname;
  int m;

#ifdef SNAPSHOT
  if (fzhlr>=0) sprintf(fname,"%s.%05d.%sitr",outfilename,fzhlr,suffix);
  else          sprintf(fname,"%s-final.%sitr",outfilename,suffix);
#else
  if (fzhlr>=0) sprintf(fname,"%s.%05d.itr",outfilename,fzhlr);
  else          sprintf(fname,"%s-final.itr",outfilename);
#endif

  out = fopen(fname,"w");
  if (NULL == out) error("Cannot write iteration file.");

  fprintf(out,"# checkpoint %d\n",fzhlr);
  fprintf(out,"startstep \t%d\n",steps+1);
#ifdef SNAPSHOT
  fprintf(out,"sscounr \t%d\n",sscount);
#endif
  fprintf(out,"box_x \t%f %f\n",box_x.x,box_x.y);
  fprintf(out,"box_y \t%f %f\n",box_y.x,box_y.y);

#if defined(NVT) || defined(NPT) || defined(STM)
  /* if we have temperature control, write external temperature and eta */
  if ( (ensemble==ENS_NVT)     || (ensemble==ENS_NPT_AXIAL) || 
       (ensemble==ENS_NPT_ISO) || (ensemble==ENS_STM) ) {
    fprintf(out,"starttemp \t%f\n",temperature);
    fprintf(out,"eta \t%f\n",eta);
  }
#endif

#ifdef FRAC 
  /* with FRAC ensemble, write actual damping factor and strainrate*/
  fprintf(out,"gamma_damp \t%f\n",gamma_damp);
  fprintf(out,"strainrate \t%f\n",dotepsilon);
#endif

#ifdef FTG
  for(m=0; m<nslices;m++)
  fprintf(out,"gamma_ftg %d\t%f\n",m , *(gamma_ftg + m));
#endif

#ifdef AND
  /* with Anderson thermostat, write external temperature */
  if (tempintv>0) fprintf(out,"starttemp \t%f\n",temperature);
#endif

#ifdef FBC
  for(m=0; m<vtypes;m++)
    fprintf(out,"extra_startforce %d %.21g %.21g \n",
            m,(fbc_forces+m)->x,(fbc_forces+m)->y);
#endif

#ifdef NPT
  /* if we have pressure control, write external pressure and xi */
  if ((ensemble==ENS_NPT_ISO) && (isq_tau_xi>0)) {
    fprintf(out,"pressure_start \t%f\n",pressure_ext.x);
    fprintf(out,"xi \t%f\n",xi.x);
  }
  if ((ensemble==ENS_NPT_AXIAL) && (isq_tau_xi>0)) {
    fprintf(out,"pressure_start \t%f %f\n",pressure_ext.x,pressure_ext.y);
    fprintf(out,"xi \t%f %f\n",xi.x,xi.y);
  }
#endif

  fclose(out);
}

#ifdef AVPOS

/******************************************************************************
*
*  write avpos iteration file
*
******************************************************************************/

void write_avpos_itr_file(int fzhlr, int steps)
{
  FILE *out;
  str255 fname;
  real tmp;

  if (fzhlr>=0) sprintf(fname,"%s.%05d.avp.itr",outfilename,fzhlr);
  else          sprintf(fname,"%s-final.avp.itr",outfilename);

  out = fopen(fname,"w");
  if (NULL == out) error("Cannot write avpos iteration file.");

  fprintf(out,"# checkpoint %d\n",fzhlr);
  fprintf(out,"startstep \t%d\n",steps+1);
#ifdef NPT
  /* Take average of box vectors */
  tmp = (real) avpos_res / ( avpos_int + avpos_res );
  fprintf(out,"box_x \t%f %f\n",av_box_x.x * tmp,av_box_x.y * tmp);
  fprintf(out,"box_y \t%f %f\n",av_box_y.x * tmp,av_box_y.y * tmp);
#else
  fprintf(out,"box_x \t%f %f\n",box_x.x,box_x.y);
  fprintf(out,"box_y \t%f %f\n",box_y.x,box_y.y);
#endif

  fclose(out);
}

#endif



