
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2012 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_io_3d.c -- 3D-specific IO routines (now also used for 2D)
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
  header_info_t info;
  cell *input;
  FILE *infile=NULL;
  char buf[1024];
  long addnumber = 0;
  int  p=0, k, maxc1=0, maxc2, count_atom;
  int  i, s, n, to_cpu, have_header=0, count;
  vektor   pos, axe;
  ivektor  cellc;
  real     m, d[MAX_ITEMS_CONFIG];
  minicell *to;
#ifdef MPI
  msgbuf   *input_buf=NULL, *b;
#endif

   if ((0 == myid) && (0==myrank)) {
    printf("Reading atoms from %s.\n", infilename); 
    fflush(stdout);
  }

  /* setup box and cell array; empty old cells, if there are any */
  if (box_from_header == 1) read_box(infilename);
  if (natoms > 0) {
    for (k=0; k<nallcells; k++) cell_array[k].n = 0;
  }
#ifdef BBOOST
    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("print "" read_atoms start "" checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("    ************************* \n");fflush(stdout);
#endif
  make_box();
#ifdef BBOOST
    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("print "" read_atoms end "" checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("    ************************* \n");fflush(stdout);
#endif
  /* initialize random number generator */
  srand48(seed);
  is_big_endian = endian();

  /* allocate num_sort and num_vsort on all CPUs */
  if ((num_sort = (long *) calloc(ntypes,sizeof(long)))==NULL)
    error("cannot allocate memory for num_sort\n");
  if ((num_vsort = (long *) calloc(vtypes,sizeof(long)))==NULL)
    error("cannot allocate memory for num_vsort\n");
#ifdef RIGID
  /* allocate num_ssort and supermass on all CPUs */
  if ( nsuperatoms>0 ) {
    if ((num_ssort = (int *) calloc(nsuperatoms,sizeof(int)))==NULL)
      error("cannot allocate memory for num_ssort\n");
    if((supermass = (real *) calloc(nsuperatoms,sizeof(real)))==NULL)
      error("cannot allocate memory for supermass\n"); 
    for (i=0; i<nsuperatoms; i++)
      supermass[i] = 0.0;
  }
#endif

#ifdef VEC
  /* allocate the space for all atoms in one step */
  atoms.n = 0;
  atoms.n_max = 0;
  atoms.n_buf = 0;
  alloc_cell(&atoms, atoms_per_cpu);
#endif

#ifdef MPI

  /* size of temporary input buffers */
  if (inp_grp_size > 1) inbuf_size /= (sizeof(real) * (inp_grp_size-1)); 

#ifndef BG
  /* Try opening first a per cpu file - not supported on BlueGene/L */
  if (1==parallel_input) {
    sprintf(buf,"%s.%u",infilename,myid); 
    infile = fopen(buf,"r");
    if (NULL!=infile) {
      if (0==myid) {
        sprintf(buf,"%s.head",infilename); 
        have_header = read_header(&info, buf);
        /* have_header==2 indicates header in separate file */
        if (have_header) have_header++;
      }
      MPI_Bcast( &have_header, 1, MPI_INT, 0, MPI_COMM_WORLD); 
      /* When each cpu reads only part of the atoms, we have to add up the
         number of atoms to get the correct natoms. We set a flag here. */
      addnumber=1;
      if (have_header) broadcast_header(&info);
    }
  } else 
#endif
  if (myid != my_inp_id) {
    recv_atoms();
    read_atoms_cleanup();
    return;
  }

  if (NULL==infile) {
    infile = fopen(infilename,"r");
    if (NULL==infile) error_str("File %s not found", infilename);
    have_header = read_header( &info, infilename );
  }

  /* allocate temporary input buffer */
  if (inp_grp_size > 1) {
    input_buf = (msgbuf *) malloc( num_cpus * sizeof(msgbuf) );
    if (NULL==input_buf) error("cannot allocate input buffers");
    input_buf[0].data  = NULL;
    input_buf[0].n_max = 0;
    input_buf[0].n     = 0;
    for (i=0; i<num_cpus; i++) {
      input_buf[i].data  = NULL;
      input_buf[i].n     = 0;
      input_buf[i].n_max = 0;
#ifdef BG
      if ((i != myid) && ((parallel_input==0) || (my_inp_grp == io_grps[i])))
#else
      if (i != myid)
#endif
        alloc_msgbuf(input_buf+i, inbuf_size);
    }
  }

#else /* not MPI */

  infile = fopen(infilename, "r");
  if (NULL==infile) error_str("File %s not found", infilename);
  have_header = read_header(&info, infilename);

#endif /* MPI or not MPI */

  /* limited backwards compatibility */
  if (have_header==0) {
    info.format   = 'A';
    info.n_number = 1;
    info.n_type   = 1;
    info.n_mass   = 1;
#ifdef UNIAX
    info.n_pos    = 6;
    info.n_vel    = 6;
    info.n_data   = 0;
    info.n_items  = 3+2*6;
#else
    info.n_pos    = DIM;
    info.n_vel    = DIM;
    info.n_data   = 0;
    info.n_items  = 3+2*DIM;
#endif
#ifdef REFPOS
    info.n_refpos_x = -1;
    if (imdrestart)
      error("Restart with ref. positions requires config file with header");
#endif
#ifdef DISLOC
    info.n_x_ref = -1;
    if ((up_ort_ref < 0) && (info.n_x_ref < 0))
      error("Reference positions require configuration file with header");
    info.n_Epot_ref = -1;
    if ((calc_Epot_ref == 0) && (info.n_Epot_ref < 0))
      error("Reference energy requires configuration file with header");
#endif
#ifdef VARCHG
    error("Option VARCHG requires configuration file with header");
#endif
  }

  /* Set up 1 atom input cell */
  input = (cell *) malloc(sizeof(cell));
  if (0==input) error("Cannot allocate input cell.");
  input->n_max=0;
  alloc_cell(input, 1);

  natoms  = 0;
  nactive = 0;
#ifdef SHOCK
  if ( steps_min ==0 ) do_maxwell=1;
#endif

  /* read away header; if have_header==2, header is in separate file */
  if (have_header==1) {
    do {
      char *s;
      s=fgets(buf,sizeof(buf),infile);
    } while (('#'!=buf[0]) || ('E'!=buf[1])); 
  }

  /* Read the input file line by line */
  while(!feof(infile)) {

    /* ASCII input */
    if (info.format == 'A') {
      p = 1;
      if (NULL==fgets(buf,sizeof(buf),infile)) p=0;
      /* eat comments */
      while (('#'==buf[0]) && !feof(infile))
        if (NULL==fgets(buf,sizeof(buf),infile)) p=0;
      if (p==0) break;
#ifdef DOUBLE
      p = sscanf( buf,
       "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
        &n,&s,&d[0],&d[1],&d[2],&d[3],&d[4],&d[5],&d[6],&d[7],&d[8],
        &d[9],&d[10],&d[11],&d[12],&d[13],&d[14],&d[15]);
#else
      p = sscanf( buf,
       "%d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
        &n,&s,&d[0],&d[1],&d[2],&d[3],&d[4],&d[5],&d[6],&d[7],&d[8],
        &d[9],&d[10],&d[11],&d[12],&d[13],&d[14],&d[15]);
#endif
    }
    /* double precision input */
    else if ((info.format=='B') || (info.format=='L')) {
      i_or_d *data = (i_or_d *) buf;
      p = fread(buf, sizeof(i_or_d), info.n_items-1, infile);
      if (p>0) p++; /* first value contains two items */
      if (info.endian == is_big_endian) {
        n = data[0].i[0];
        s = data[0].i[1];
        for (k=0; k < info.n_items-2; k++) { 
          d[k] = (real) data[k+1].d;
	}
      }
      else {
        n = SwappedInteger(data[0].i[0]);
        s = SwappedInteger(data[0].i[1]);
        for (k=0; k < info.n_items-2; k++) 
          d[k] = (real) SwappedDouble(data[k+1].d);
      }
    }
    /* single precision input */
    else if ((info.format=='b') || (info.format=='l')) {
      i_or_f *data = (i_or_f *) buf;
      p = fread(buf, sizeof(i_or_f), info.n_items, infile);
      if (info.endian == is_big_endian) {
        n = data[0].i;
        s = data[1].i;
        for (k=0; k < info.n_items-2; k++) 
          d[k] = (real) data[k+2].f;
      }
      else {
        n = SwappedInteger(data[0].i);
        s = SwappedInteger(data[1].i);
        for (k=0; k < info.n_items-2; k++) 
          d[k] = (real) SwappedFloat(data[k+2].f);
      }
    }
    /* skip empty lines at end of file */
    if (p>0) {

      if ((p < info.n_items) && (p != info.n_items - info.n_vel)) {
        error("incorrect line in configuration file.");
      }

      input->n = 1;
      count = 0;
#ifndef MONOLJ
      NUMMER(input,0) = n;
#ifndef MONO
      SORTE (input,0) = MOD(s,ntypes);
#endif
      if (s>=vtypes) error("atom type must not exceed total_types");
      VSORTE(input,0) = s;
      if (info.n_mass) m = d[count++];
      else             m = masses[SORTE(input,0)];
      if (0>=m) error("Mass zero or negative.\n");
      MASSE (input,0) = m;
#ifdef RIGID
      if ( nsuperatoms > 0 )
	if ( superatom[s] > -1 ) 
	  supermass[superatom[s]] += m;
#endif
#endif
      pos.x = d[count++];
      pos.y = d[count++];
#ifndef TWOD
      pos.z = d[count++];
#endif
      /* with per CPU input files, avoid back_into_box */
      if (0==addnumber) pos = back_into_box(pos);
      ORT(input,0,X) = pos.x;
      ORT(input,0,Y) = pos.y;
#ifndef TWOD
      ORT(input,0,Z) = pos.z;
#endif
#ifdef UNIAX
      ACHSE(input,0,X) = axe.x = d[count++];
      ACHSE(input,0,Y) = axe.y = d[count++];
      ACHSE(input,0,Z) = axe.z = d[count++];
      if (ABS( SPROD(axe,axe) - 1.0 ) > 1.0e-04)
        error("Molecular axis not a unit vector!");
      if ((0 == info.n_vel) || (p < info.n_items)) {
	do_maxwell=1;
	IMPULS(input,0,X) = 0 ;
	IMPULS(input,0,Y) = 0 ;
	IMPULS(input,0,Z) = 0 ;
	DREH_IMPULS(input,0,X) = 0 ;
	DREH_IMPULS(input,0,Y) = 0 ;
	DREH_IMPULS(input,0,Z) = 0 ;
      } else {
	IMPULS(input,0,X) = d[count++] * m;
	IMPULS(input,0,Y) = d[count++] * m;
	IMPULS(input,0,Z) = d[count++] * m;
	DREH_IMPULS(input,0,X) = d[count++] * uniax_inert;
	DREH_IMPULS(input,0,Y) = d[count++] * uniax_inert;
	DREH_IMPULS(input,0,Z) = d[count++] * uniax_inert;
      }
      DREH_MOMENT(input,0,X) = 0 ;
      DREH_MOMENT(input,0,Y) = 0 ;
      DREH_MOMENT(input,0,Z) = 0 ;
#else /* not UNIAX */
      if ((0 == info.n_vel) || (p < info.n_items)) {
        do_maxwell=1;
        IMPULS(input,0,X) = 0;
        IMPULS(input,0,Y) = 0;
#ifndef TWOD
        IMPULS(input,0,Z) = 0;
#endif
      } else {
        IMPULS(input,0,X) = d[count++] * m * (restrictions+s)->x;
        IMPULS(input,0,Y) = d[count++] * m * (restrictions+s)->y;
#ifndef TWOD
        IMPULS(input,0,Z) = d[count++] * m * (restrictions+s)->z;
#endif
      }
#endif /* UNIAX or not UNIAX */
      KRAFT(input,0,X) = 0;
      KRAFT(input,0,Y) = 0;
#ifndef TWOD
      KRAFT(input,0,Z) = 0;
#endif
#ifdef DIPOLE 			/* only 3D */
      DP_P_IND(input,0,X) = d[count++];
      DP_P_IND(input,0,Y) = d[count++];
      DP_P_IND(input,0,Z) = d[count++];
#endif /* DIPOLE */
#ifdef REFPOS
      if (info.n_refpos_x > 0) {
        REF_POS(input,0,X) = d[info.n_refpos_x-2];
        REF_POS(input,0,Y) = d[info.n_refpos_x-1];
#ifndef TWOD
        REF_POS(input,0,Z) = d[info.n_refpos_x  ];
#endif
      }
#endif
#ifdef DISLOC
      if (info.n_x_ref > 0) {
        ORT_REF(input,0,X) = d[info.n_x_ref-2];
        ORT_REF(input,0,Y) = d[info.n_x_ref-1];
#ifndef TWOD
        ORT_REF(input,0,Z) = d[info.n_x_ref  ];
#endif
      }
      if (info.n_Epot_ref > 0) {
        EPOT_REF(input,0)  = d[info.n_Epot_ref-2];
      }
#endif
#ifdef VARCHG
      if (info.n_charge > 0) {
        CHARGE(input,0)    = d[info.n_charge-2];
      } 
      else {
        CHARGE(input,0)    = charge[ SORTE(input,0) ];
      }
#endif

#ifdef EPITAX
      /* determine largest atom number of substrate atoms */
      epitax_sub_n = MAX( n, epitax_sub_n );
#endif

#ifdef TWOD
      cellc = cell_coord(pos.x,pos.y);
#else
      cellc = cell_coord(pos.x,pos.y,pos.z);
#endif

#ifdef BUFCELLS

      to_cpu = cpu_coord(cellc);
      count_atom = 0;

#ifdef MPI

      /* to_cpu is in my input group, but not myself */
      if ((inp_grp_size > 1) && (myid != to_cpu)) {
        b = input_buf + to_cpu;
        if (b->data != NULL) {
          copy_atom_cell_buf(b, to_cpu, input, 0);
          if (b->n_max - b->n < atom_size) {
            MPI_Send(b->data, b->n, REAL, to_cpu, INBUF_TAG, cpugrid);
            b->n = 0;
          }
          count_atom = 1;
	}
      } else

#endif /* MPI */

      /* with per CPU input files, make sure we keep all atoms */
      if ((to_cpu==myid) || ((1==parallel_input) && (1==addnumber))) {
	cellc = local_cell_coord(cellc);
        to = PTR_VV(cell_array,cellc,cell_dim);
	INSERT_ATOM(to, input, 0);
        count_atom = 1;
      }

#else /* not BUFCELLS */

      to = PTR_VV(cell_array,cellc,cell_dim);
      INSERT_ATOM(to, input, 0);
      count_atom = 1;

#endif /* BUFCELLS or not BUFCELLS */

      if (count_atom) {
        natoms++;  
        /* we still have s == input->vsorte[0] */
        if (s < ntypes) {
          nactive += DIM;
#ifdef UNIAX
          nactive_rot += 2;
#endif
        } else {
          nactive += (long) (restrictions+s)->x;
          nactive += (long) (restrictions+s)->y;
#ifndef TWOD
          nactive += (long) (restrictions+s)->z;
#endif
        }
        num_sort [ SORTE(input,0)]++;
        num_vsort[VSORTE(input,0)]++;
      }
    } /* (p>0) */
  } /* !feof(infile) */
  fclose(infile);  

#ifdef MPI
  if (inp_grp_size > 1) {
    /* The last buffer is sent with a different tag, which tells the
       target CPU that reading is finished; we increase the size by
       one, so that the buffer is sent even if it is empty */
    for (s=0; s<num_cpus; s++) {
      b = input_buf + s;
      if (b->data) {
        b->data[b->n++] = 0.0;
        MPI_Send(b->data, b->n, REAL, s, INBUF_TAG+1, cpugrid);
        free_msgbuf(b);
      }
    }
    free(input_buf);
  }
#endif

  /* the following routine contains stuff that must be executed
     on all CPUs, also if parallel_input = 0 */
  read_atoms_cleanup();

}

/******************************************************************************
*
*  read_atoms_cleanup contains stuff that must be executed
*  on all CPUs, also if parallel_input = 0
*
******************************************************************************/

void read_atoms_cleanup(void)
{
  int  k, i, s, maxc1=0, maxc2;
  long tmp;

#ifdef MPI

  /* Add the number of atoms read (and kept) by each CPU */
  if (1==parallel_input) {
    MPI_Allreduce( &natoms,  &tmp, 1, MPI_LONG, MPI_SUM, cpugrid);
    natoms = tmp;
    MPI_Allreduce( &nactive, &tmp, 1, MPI_LONG, MPI_SUM, cpugrid);
    nactive = tmp;
#ifdef UNIAX
    MPI_Allreduce( &nactive_rot, &tmp, 1, MPI_LONG, MPI_SUM, cpugrid);
    nactive_rot = tmp;
#endif
    for (i=0; i<ntypes; i++) {
      MPI_Allreduce( &num_sort[i], &tmp, 1, MPI_LONG, MPI_SUM, cpugrid);
      num_sort[i] = tmp;
    }
    for (i=0; i<vtypes; i++) {
      MPI_Allreduce( &num_vsort[i], &tmp, 1, MPI_LONG, MPI_SUM, cpugrid);
      num_vsort[i] = tmp;
    }
  } else { /* broadcast if serial io */
    MPI_Bcast( &natoms ,      1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast( &nactive,      1, MPI_LONG, 0, MPI_COMM_WORLD);
#ifdef UNIAX
    MPI_Bcast( &nactive_rot,  1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif
    MPI_Bcast( num_sort, ntypes, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast( num_vsort,vtypes, MPI_LONG, 0, MPI_COMM_WORLD);
  }

  /* If CPU 0 found velocities in its data, no initialisation is done */
  MPI_Bcast( &do_maxwell, 1, MPI_INT, 0, MPI_COMM_WORLD);

#endif /* MPI */

#ifdef RIGID
  if ( nsuperatoms > 0 ) {
    if (0==myid) {
      /* number of atoms belonging to superatoms */
      for (i=0; i<vtypes; i++) {
	if ( superatom[i] > -1 )
	  num_ssort[superatom[i]] += num_vsort[i];
      }
    }
#ifdef MPI
    MPI_Bcast( num_ssort, nsuperatoms, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

    if (0==myid) {
      /* recompute nactive */
      nactive = 0;
      for (i=0; i<vtypes; i++) {
	/* translational degrees of freedom of nonrigid vtypes */
	if (superatom[i]==-1 || ((superrestrictions+superatom[i])->x)==0 )
	  nactive += num_vsort[i] * (long) (restrictions+i)->x;
	if (superatom[i]==-1 || ((superrestrictions+superatom[i])->y)==0 )
	  nactive += num_vsort[i] * (long) (restrictions+i)->y;
#ifndef TWOD
	if (superatom[i]==-1 || ((superrestrictions+superatom[i])->z)==0 )
	  nactive += num_vsort[i] * (long) (restrictions+i)->z;
#endif
      }
      for (s=0; s<nsuperatoms; s++) {
	/* check whether superatom s is mobile */
        vektor3d mobile = {1, 1, 1};
	for(i=0; i<vtypes; i++)
	  if (superatom[i]==s ) {
	    if ((restrictions+i)->x == 0) 
	      mobile.x = 0;
	    if ((restrictions+i)->y == 0)
	      mobile.y = 0;
#ifndef TWOD
	    if ((restrictions+i)->z == 0)
	      mobile.z = 0;
#endif
	  }
	/* translational degrees of freedom of rigid vtypes */
	nactive += (superrestrictions+s)->x * mobile.x
                +  (superrestrictions+s)->y * mobile.y;
#ifndef TWOD
        nactive += (superrestrictions+s)->z * mobile.z;
#endif
      }
    }
  }
#ifdef MPI
  MPI_Bcast( &nactive, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

#endif /* RIGID */

  /* print number of atoms */
  if ((0==myid) && (0==myrank)) {
    printf("Read structure with %ld atoms.\n",natoms);
    tmp = num_sort[0];
    printf("num_sort = [ %ld",num_sort[0]);
    for (i=1; i<ntypes; i++) {
      printf(", %ld",num_sort[i]);
      tmp += num_sort[i];
    }
    printf(" ],  total = %ld\n",tmp);
    tmp = num_vsort[0];
    printf("num_vsort = [ %ld",num_vsort[0]);
    for (i=1; i<vtypes; i++) {
      printf(", %ld",num_vsort[i]);
      tmp += num_vsort[i];
    }
    printf(" ],  total = %ld\n",tmp);
#ifdef RIGID
    if ( nsuperatoms>0 ) {
      printf("num_ssort = [ %u",num_ssort[0]);
      for (i=1; i<nsuperatoms; i++)
	printf(", %u",num_ssort[i]);
      printf(" ]\n");
    }
#endif
  }

  /* determine maximal cell occupancy */
  for (k=0; k<ncells; k++) maxc1 = MAX( maxc1, (cell_array+CELLS(k))->n );
#ifdef MPI
  MPI_Reduce( &maxc1, &maxc2, 1, MPI_INT, MPI_MAX, 0, cpugrid);
#else
  maxc2 = maxc1;
#endif
  if ((0 == myid) && (0==myrank)) { 
    printf("maximal cell occupancy: %d\n", maxc2);
    fflush(stdout);  
  }
}

#ifdef MPI

/******************************************************************************
*
*  recveive atoms in several chunks from CPU 0
*  this is only used when parallel_input==0, or on BlueGene/L
*
******************************************************************************/

void recv_atoms(void)
{
  MPI_Status status;
  int finished = 0, src = my_inp_id;
  msgbuf b = {NULL,0,0};

  alloc_msgbuf(&b, inbuf_size);
  do {
    MPI_Recv(b.data, inbuf_size, REAL, src, MPI_ANY_TAG, cpugrid, &status);
    MPI_Get_count(&status, REAL, &b.n);
    if (status.MPI_TAG==INBUF_TAG+1) { b.n--; finished=1; } /* last buffer */
    process_buffer( &b );
  } while (0==finished);
  free_msgbuf(&b);
}

#endif /* MPI */ 


/******************************************************************************
*
*  filter function for write_config_select
*  writes data of all atoms for checkpoints
*
******************************************************************************/

void write_atoms_config(FILE *out)
{
  int i, k, n, len=0;
  i_or_r *data;
#ifdef CNA
  int crist, nn, nn_other, nn_1421, nn_1422;
#endif

#ifdef HPO
#define RESOL1 " %12.16f"
#define RESOL2 " %12.16f %12.16f"
#define RESOL3 " %12.16f %12.16f %12.16f"
#else
#define RESOL1 " %f"
#define RESOL2 " %f %f"
#define RESOL3 " %f %f %f"
#endif

  for (k=0; k<NCELLS; k++) {
    cell *p;
    p = CELLPTR(k);
    for (i=0; i<p->n; i++) {

#ifdef CNA
      if (cna_crist>0) {
	nn_other = MARK(p,i) % 100; 
	nn_1422  = ( MARK(p,i) / 100 ) % 100;
	nn_1421  = ( MARK(p,i) / 10000 ) % 100;
	nn       = nn_1421 + nn_1422 + nn_other;
	  
	/* fcc */
	if ( nn == 12 && nn_1421 == 12 )
	  crist = 0;
	/* hcp */
	else if ( nn == 12 && nn_1421 == 6 && nn_1422 == 6 )
	  crist = 1;
	/* other 12 */
	else if ( nn == 12 )
	  crist = 2;
	/* other */
	else
	  crist = 3;
      }
#endif

      if (binary_output) {
        n = 0;
        data = (i_or_r *) (outbuf+len);
#ifdef DOUBLE
        data[n  ].i[0] = NUMMER(p,i);
        data[n++].i[1] = VSORTE(p,i);
#else
        data[n++].i    = NUMMER(p,i);
        data[n++].i    = VSORTE(p,i);
#endif
        data[n++].r = MASSE(p,i);
        data[n++].r = ORT(p,i,X);
        data[n++].r = ORT(p,i,Y);
#ifndef TWOD
        data[n++].r = ORT(p,i,Z);
#endif
#ifdef UNIAX
        data[n++].r = ACHSE(p,i,X);
        data[n++].r = ACHSE(p,i,Y);
        data[n++].r = ACHSE(p,i,Z);
#endif
        if (ensemble != ENS_CG) {
          data[n++].r = (IMPULS(p,i,X) / MASSE(p,i));
          data[n++].r = (IMPULS(p,i,Y) / MASSE(p,i));
#ifndef TWOD
          data[n++].r = (IMPULS(p,i,Z) / MASSE(p,i));
#endif
	}
#ifdef UNIAX
        data[n++].r = DREH_IMPULS(p,i,X) / uniax_inert;
        data[n++].r = DREH_IMPULS(p,i,Y) / uniax_inert;
        data[n++].r = DREH_IMPULS(p,i,Z) / uniax_inert; 
#endif
        data[n++].r = POTENG(p,i);
#if defined(VARCHG) || defined(EWALD) || defined(USEFCS)
        data[n++].r = CHARGE(p,i);
#endif
#ifdef NNBR
        data[n++].r = (real) NBANZ(p,i);
#endif
#ifdef REFPOS
        data[n++].r = REF_POS(p,i,X);
        data[n++].r = REF_POS(p,i,Y);
#ifndef TWOD
        data[n++].r = REF_POS(p,i,Z);
#endif
#endif
#ifdef DISLOC
        data[n++].r = ORT_REF(p,i,X);
        data[n++].r = ORT_REF(p,i,Y);
#ifndef TWOD
        data[n++].r = ORT_REF(p,i,Z);
#endif
        data[n++].r = EPOT_REF(p,i);
#endif
#if defined(EAM2) && !defined(NORHOH)
        data[n++].r = EAM_RHO(p,i);
#ifdef EEAM
        data[n++].r = EAM_P(p,i);
#endif
#endif
#ifdef DAMP
        data[n++].r = DAMPF(p,i);
#endif
#ifdef ADA
#ifdef DOUBLE
        data[n++].i[0] = ADATYPE(p,i);
#else
        data[n++].i = ADATYPE(p,i);
#endif
#endif
#ifdef NYETENSOR
		nyeTensorInfo *info = NYE(p,i);
		if (info != NULL){
			data[n++].r = info->ls[0];
			data[n++].r = info->ls[1];
			data[n++].r = info->ls[2];
			data[n++].r = info->bv[0];
			data[n++].r = info->bv[1];
			data[n++].r = info->bv[2];
		} else {
			data[n++].r	= 0.; data[n++].r = 0.; data[n++].r = 0.;
			data[n++].r	= 0.; data[n++].r = 0.; data[n++].r = 0.;
		}
#endif
#ifdef CNA
	if (cna_crist>0) {	  
	  data[n++].r = (real) crist;
	}
#endif
        len += n * sizeof(i_or_r);
      }
      else {
        len += sprintf(outbuf+len, "%d %d", NUMMER(p,i), VSORTE(p,i));
        len += sprintf(outbuf+len, RESOL1, MASSE(p,i));
#ifdef TWOD
        len += sprintf(outbuf+len, 
          RESOL2, ORT(p,i,X), ORT(p,i,Y) );
#else
        len += sprintf(outbuf+len, 
          RESOL3, ORT(p,i,X), ORT(p,i,Y), ORT(p,i,Z) );
#endif
#ifdef UNIAX
        len += sprintf(outbuf+len, 
          RESOL3, ACHSE(p,i,X), ACHSE(p,i,Y), ACHSE(p,i,Z) );
#endif
#ifdef TWOD
        if (ensemble != ENS_CG)
          len += sprintf(outbuf+len, RESOL2,
            IMPULS(p,i,X) / MASSE(p,i), 
            IMPULS(p,i,Y) / MASSE(p,i) );
#else
        if (ensemble != ENS_CG)
          len += sprintf(outbuf+len, RESOL3,
            IMPULS(p,i,X) / MASSE(p,i), 
            IMPULS(p,i,Y) / MASSE(p,i), 
            IMPULS(p,i,Z) / MASSE(p,i) );
#endif
#ifdef UNIAX
        len += sprintf(outbuf+len, RESOL3,
          DREH_IMPULS(p,i,X) / uniax_inert,
          DREH_IMPULS(p,i,Y) / uniax_inert,
          DREH_IMPULS(p,i,Z) / uniax_inert ); 
#endif
        len += sprintf(outbuf+len, RESOL1, POTENG(p,i));
#if defined(VARCHG) || defined(EWALD) || defined(USEFCS)
        len += sprintf(outbuf+len, RESOL1, CHARGE(p,i));
#endif
#ifdef NNBR
        len += sprintf(outbuf+len, " %d",  NBANZ(p,i));
#endif
#ifdef REFPOS
#ifdef TWOD
        len += sprintf(outbuf+len, 
          RESOL2, REF_POS(p,i,X), REF_POS(p,i,Y));
#else
        len += sprintf(outbuf+len, 
          RESOL3, REF_POS(p,i,X), REF_POS(p,i,Y), REF_POS(p,i,Z));
#endif
#endif
#ifdef DISLOC
#ifdef TWOD
        len += sprintf(outbuf+len, 
          RESOL2, ORT_REF(p,i,X), ORT_REF(p,i,Y));
#else
        len += sprintf(outbuf+len, 
          RESOL3, ORT_REF(p,i,X), ORT_REF(p,i,Y), ORT_REF(p,i,Z));
#endif
        len += sprintf(outbuf+len, RESOL1, EPOT_REF(p,i));
#endif
#if defined(EAM2) && !defined(NORHOH)
        len += sprintf(outbuf+len, RESOL1, EAM_RHO(p,i));
#ifdef EEAM
        len += sprintf(outbuf+len, RESOL1, EAM_P(p,i));
#endif
#endif
#ifdef DIPOLE
	len += sprintf(outbuf+len, RESOL3, 
		       DP_P_IND(p,i,X), DP_P_IND(p,i,Y), DP_P_IND(p,i,Z)); 
#endif	/* DIPOLE */
#ifdef DAMP
        len += sprintf(outbuf+len, RESOL1, DAMPF(p,i));
#endif
#ifdef ADA
        len += sprintf(outbuf+len, " %d", ADATYPE(p,i));
#endif
#ifdef NYETENSOR
        nyeTensorInfo *info = NYE(p,i);
       	if (info != NULL) {
       		len += sprintf(outbuf+len, RESOL3, info->ls[0], info->ls[1], info->ls[2] );
       		len += sprintf(outbuf+len, RESOL3, info->bv[0], info->bv[1], info->bv[2] );
       	} else {
       		len += sprintf(outbuf + len," 0. 0. 0. 0. 0. 0.");
       	}
#endif
#ifdef CNA
	if (cna_crist>0) 
	  len += sprintf(outbuf+len, " %d", crist);
#endif
        len += sprintf(outbuf+len,"\n");
      }
      /* flush or send outbuf if it is full */
      if (len > outbuf_size - 256) flush_outbuf(out,&len,OUTBUF_TAG);
    }
  }
  flush_outbuf(out,&len,OUTBUF_TAG+1);
}

/******************************************************************************
*
*  write iteration file
*
******************************************************************************/

void write_itr_file(int fzhlr, int steps, char *suffix)
{
  FILE *out;
  str255 fname;
  int m, n;

  if (strcasecmp(suffix,"ss")==0) {
    if (fzhlr>=0) sprintf(fname,"%s.%05d.%sitr",outfilename,fzhlr,suffix);
    else          sprintf(fname,"%s-final.%sitr",outfilename,suffix);
  }
  else {
    if (fzhlr>=0)       sprintf(fname,"%s.%05d.itr",outfilename,fzhlr);
    else if (fzhlr==-1) sprintf(fname,"%s-final.itr",outfilename);
    else                sprintf(fname,"%s-interm.itr",outfilename);
  }

  out = fopen(fname,"w");
  if (NULL == out) error("Cannot write iteration file.");

  fprintf(out, "# checkpoint %d\n", fzhlr);
  fprintf(out, "startstep \t%d\n", steps+1);
#ifdef RELAX
  fprintf(out, "sscount \t%d\n", sscount);
  fprintf(out, "nfc \t%d\n", nfc);
#endif

#ifdef TWOD
  fprintf(out, "box_x \t%.16f %.16f\n", box_x.x, box_x.y);
  fprintf(out, "box_y \t%.16f %.16f\n", box_y.x, box_y.y);
#else
  fprintf(out, "box_x \t%.16f %.16f %.16f\n", box_x.x, box_x.y, box_x.z);
  fprintf(out, "box_y \t%.16f %.16f %.16f\n", box_y.x, box_y.y, box_y.z);
  fprintf(out, "box_z \t%.16f %.16f %.16f\n", box_z.x, box_z.y, box_z.z);
#endif

#if defined(NVT) || defined(NPT) 
  /* if we have temperature control, write external temperature and eta */
  if (((ensemble==ENS_NVT) || (ensemble==ENS_NPT_AXIAL) || 
       (ensemble==ENS_NPT_ISO)) && (isq_tau_eta>0)) {
    fprintf(out, "starttemp \t%f\n", temperature);
    fprintf(out, "eta \t%f\n", eta);
#ifdef UNIAX
    fprintf(out, "eta_rot \t%f\n", eta_rot);
#endif
  }
#endif

#ifdef FRAC 
 /* with FRAC ensemble, write actual damping factor and strainrate*/
  fprintf(out, "gamma_damp \t%f\n",gamma_damp);
  fprintf(out, "strainrate \t%f\n",dotepsilon);
#endif

#ifdef DAMP
#ifdef TWOD
  fprintf(out, "center   \t%f %f\n", center.x,   center.y  );
  fprintf(out, "stadium  \t%f %f\n", stadium.x,  stadium.y );
#else
  fprintf(out, "center   \t%f %f %f\n", center.x,   center.y,   center.z  );
  fprintf(out, "stadium  \t%f %f %f\n", stadium.x,  stadium.y,  stadium.z );
  fprintf(out, "stadium2 \t%f %f %f\n", stadium2.x, stadium2.y, stadium2.z);
#endif
#endif

#ifdef FTG
  for(m=0; m<nslices;m++)
    fprintf(out, "gamma_ftg %d\t%f\n", m, *(gamma_ftg + m));
#endif

#ifdef AND
  /* with Anderson thermostat, write external temperature */
  if (tempintv>0) fprintf(out, "starttemp \t%f\n", temperature);
#endif

#ifdef BER
  /* with Berendsen thermostat, write external temperature */
  if (tauber>0) fprintf(out, "starttemp \t%f\n", temperature);
#endif

#ifdef FBC
  for(n=0; n<vtypes;n++)
#ifdef TWOD
    fprintf(out, "extra_startforce %d %.21g %.21g\n",
            n, (fbc_forces+n)->x, (fbc_forces+n)->y );
#else
    fprintf(out, "extra_startforce %d %.21g %.21g %.21g \n",
            n, (fbc_forces+n)->x, (fbc_forces+n)->y, (fbc_forces+n)->z);
#endif
#endif
    
#ifdef BEND
  for(n=0; n<vtypes;n++)
#ifdef TWOD
    fprintf(out, "extra_startbforce %d %.21g %.21g\n",
            n, (fbc_bforces+n)->x, (fbc_bforces+n)->y );
#else
    fprintf(out, "extra_startbforce %d %.21g %.21g %.21g \n",
            n, (fbc_bforces+n)->x, (fbc_bforces+n)->y, (fbc_bforces+n)->z);
#endif
#endif

#ifdef NPT
  /* if we have pressure control, write external pressure and xi */
  if ((ensemble==ENS_NPT_ISO) && (isq_tau_xi>0)) {
    fprintf(out, "pressure_start \t%f\n", pressure_ext.x);
    fprintf(out, "xi \t%f\n", xi.x);
  }
  if ((ensemble==ENS_NPT_AXIAL) && (isq_tau_xi>0)) {
#ifdef TWOD
    fprintf(out, "pressure_start \t%f %f\n",
            pressure_ext.x, pressure_ext.y);
    fprintf(out,"xi \t%f %f\n", xi.x, xi.y);
#else
    fprintf(out, "pressure_start \t%f %f %f\n",
            pressure_ext.x, pressure_ext.y, pressure_ext.z);
    fprintf(out,"xi \t%f %f %f\n", xi.x, xi.y, xi.z);
#endif
  }
#endif /* NPT */

#ifdef NVX
  fprintf(out, "hc_count \t%d\n", hc_count);
#endif

#ifdef EXTPOT
  /* write positions of external potential */
  for (n=0; n<ep_n; ++n) {
    fprintf(out,"ep_pos    %d %.21g %.21g %.21g \n",
            n, ep_pos[n].x, ep_pos[n].y, ep_pos[n].z);
  }
#endif

#ifdef AVPOS
    fprintf(out,"avpos_nwrites \t%d\n",avpos_nwrites);
    fprintf(out,"avpos_npwrites \t%d\n",avpos_npwrites);
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

  fprintf(out, "# avpos checkpoint %d\n", fzhlr);
  fprintf(out, "startstep \t%d\n", steps+1);

#ifdef NPT

  /* Take average of box vectors */
  tmp = (real) avpos_res / ( avpos_int + avpos_res );
#ifdef TWOD
  fprintf(out,"box_x \t%.16f %.16f\n", av_box_x.x * tmp, av_box_x.y * tmp);
  fprintf(out,"box_y \t%.16f %.16f\n", av_box_y.x * tmp, av_box_y.y * tmp);
#else
  fprintf(out, "box_x \t%.16f %.16f %.16f\n",
          av_box_x.x * tmp, av_box_x.y * tmp, av_box_x.z * tmp);
  fprintf(out, "box_y \t%.16f %.16f %.16f\n",
          av_box_y.x * tmp, av_box_y.y * tmp, av_box_y.z * tmp);
  fprintf(out, "box_z \t%.16f %.16f %.16f\n",
          av_box_z.x * tmp, av_box_z.y * tmp, av_box_z.z * tmp);
#endif

#else 

#ifdef TWOD
  fprintf(out, "box_x \t%.16f %.16f\n", box_x.x, box_x.y);
  fprintf(out, "box_y \t%.16f %.16f\n", box_y.x, box_y.y);
#else
  fprintf(out, "box_x \t%.16f %.16f %.16f\n", box_x.x, box_x.y, box_x.z);
  fprintf(out, "box_y \t%.16f %.16f %.16f\n", box_y.x, box_y.y, box_y.z);
  fprintf(out, "box_z \t%.16f %.16f %.16f\n", box_z.x, box_z.y, box_z.z);
#endif

#endif

  fclose(out);
}

#endif
