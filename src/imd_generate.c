
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
* imd_generate.c -- generate configurations online
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/******************************************************************************
*
* generate_atoms - generate atoms for initial configuration
*
* filenames starting with an _ don't specify a file to read from
* but a crystal structure to generate as an initial configuration:
*
* _hex        -- generates 2D hexagonal crystal
* _fcc        -- generates FCC structure
* _bcc        -- generates BCC structure
* _b2         -- generates B2 structure
* _l12        -- generates L1_2 structure
* _nacl       -- generates NaCl structure
* _diamond    -- generates cubic diamond structure
* _zincblende -- generates zincblende structure
* _lav        -- generates a cubic Laves structure C15 (MgCu2)
* _tiqc       -- generates a truncated icosahedra quasicrystal
*
* The lattice constant of the conventional unit cell of the crystal 
* structures is box_unit.
*
******************************************************************************/

void generate_atoms(str255 mode)
{
  cell *to;
  long addnumber, tmp, ninc, i, k;
#ifdef MPI
  MPI_Status status;
#endif

  if ((num_sort = (long *) calloc(ntypes,sizeof(long)))==NULL) {
      error("cannot allocate memory for num_sort\n");
  }
  if ((num_vsort = (long *) calloc(vtypes,sizeof(long)))==NULL) {
      error("cannot allocate memory for num_vsort\n");
  }

#ifdef DISLOC
  calc_Epot_ref=1;
#endif

  /* We always have to initialize the velocities on generated structures */
  do_maxwell=1;

#ifdef TWOD
  if (0 == strcmp(mode,"_hex")) {          /* 2D hexagonal crystal */
    init_hex();
    generate_hex();
#else /* 3D */
  if (0 == strcmp(mode,"_fcc")) {          /* FCC */
    init_cubic();
    generate_fcc(0);
  } else if (0 == strcmp(mode,"_nacl")) {  /* NaCl */
    init_cubic();
    generate_fcc(1);
  } else if (0 == strcmp(mode,"_bcc")) {   /* BCC */
    init_cubic();
    generate_fcc(2);
  } else if (0 == strcmp(mode,"_b2")) {    /* B2, CsCl */
    init_cubic();
    generate_fcc(3);
  } else if (0 == strcmp(mode,"_diamond")) {  /* cubic diamond */
    init_cubic();
    generate_fcc(4);
  } else if (0 == strcmp(mode,"_zincblende")) { /* zincblende */
    init_cubic();
    generate_fcc(5);
  } else if (0 == strcmp(mode,"_cu3au")) { /* L1_2, Cu3Au */
    init_cubic();
    generate_fcc(6);
  } else if (0 == strcmp(mode,"_cuau3")) { /* L1_2, CuAu3 */
    init_cubic();
    generate_fcc(7);
  } else if (0 == strcmp(mode,"_lav")) {   /* C15 Laves */
    init_cubic();
    generate_lav();
#ifdef QUASI
  } else if (0 == strcmp(mode,"_tiqc")) {  /* quasicrystal */
    init_qc();
    generate_qc();
#endif
#endif /* 3D */
  } else if (0==myid) error("Filename with _ specifies unknown structure.");

#ifdef MPI
#ifndef MONOLJ
  /* Get numbering of atoms right even across CPUs */
  ninc = 0;
  if ((0==myid) && (num_cpus>1))
     MPI_Send( &natoms, 1, MPI_LONG, myid + 1, 0, cpugrid );
  else if ((0<myid) && (myid<(num_cpus-1))) {
     MPI_Recv( &tmp,    1, MPI_LONG, myid - 1, 0, cpugrid, &status );
     ninc = tmp;
     tmp += natoms;
     MPI_Send( &tmp,    1, MPI_LONG, myid + 1, 0, cpugrid );
  }
  else if (myid==(num_cpus-1)) {
     MPI_Recv( &tmp,    1, MPI_LONG, myid - 1, 0, cpugrid, &status );
     ninc = tmp;
  }
  /* loop over all atoms, fix numbers */
  for (k=0; k<NCELLS; ++k) {
    int i;
    cell *p;
    p = CELLPTR(k);
    for (i=0; i<p->n; ++i) NUMMER(p,i) += ninc;
  }
#endif /* MONOLJ */

  MPI_Allreduce( &natoms,  &addnumber, 1, MPI_LONG, MPI_SUM, cpugrid);
  natoms = addnumber;
  MPI_Allreduce( &nactive, &addnumber, 1, MPI_LONG, MPI_SUM, cpugrid);
  nactive = addnumber;
  for (i=0; i<ntypes; i++) {
    MPI_Allreduce(&num_sort[i], &addnumber, 1, MPI_LONG, MPI_SUM, cpugrid);
    num_sort[i] = addnumber;
  }
#endif /* MPI */

  for (i=0; i<ntypes; i++) num_vsort[i] = num_sort[i];
  if (0 == myid) {
    printf("Generated %s structure with %ld atoms.\n",mode,natoms);
    addnumber=num_sort[0];
    printf("num_sort = [ %ld",num_sort[0]);
    for (i=1; i<ntypes; i++) {
      printf(", %ld",num_sort[i]);
      addnumber += num_sort[i];
    }
    printf(" ],  total = %ld\n",addnumber);
  }
}


#ifdef TWOD

/* initialize for hexagonal crystal */
void init_hex(void)
{
  if ((box_param.x==0) || (box_param.y==0)) error("box_param not set!");
  box_x.x = box_param.x * sqrt(3.0) * box_unit;
  box_x.y = 0.0;
  box_y.x = 0.0;
  box_y.y = box_param.y * box_unit;
  make_box();
}

/* generate hexagonal crystal */
void generate_hex()
{
  cell    *input, *to;
  ivektor min, max, cellc, lcellc;
  int     to_cpu;
  int     i, j, typ;
  real    x, y;

#ifdef MPI
  if (myid==0)
    if ((box_param.x % cpu_dim.x) || (box_param.y % cpu_dim.y))
      error("box_param must be commensurate with cpu_dim");
#endif

#ifdef MPI
  min.x =  my_coord.x      * 2 * box_param.x / cpu_dim.x;
  max.x = (my_coord.x + 1) * 2 * box_param.x / cpu_dim.x;
  min.y =  my_coord.y      * 2 * box_param.y / cpu_dim.y;
  max.y = (my_coord.y + 1) * 2 * box_param.y / cpu_dim.y;
#else
  min.x = 0; max.x = 2 * box_param.x;
  min.y = 0; max.y = 2 * box_param.y;
#endif

  /* Set up 1 atom input cell */
  input = (cell *) malloc(sizeof(cell));
  if (0==input) error("Cannot allocate input cell");
  input->n_max = 0;
  alloc_cell(input,1);

  natoms  = 0;
  nactive = 0;

  for (i=min.x ; i<max.x; i++)
    for (j=min.y; j<max.y; j++) {

      typ  = (i+j) % 2;
      if (typ > 0) continue;

      x = (i+0.5) * sqrt(3.0) * 0.5 * box_unit;
      y = (j+0.5) * 0.5 * box_unit;

      natoms++;
      nactive += 2;

      input->n = 1;
      ORT(input,0,X)  = x;
      ORT(input,0,Y)  = y;
      NUMMER(input,0) = natoms;
#ifndef MONO
      SORTE (input,0) = typ;
#endif
      VSORTE(input,0) = typ;
      MASSE(input,0)  = masses[typ];
      num_sort[typ]++;
      cellc = cell_coord(x,y);
#ifdef MPI
      to_cpu = cpu_coord(cellc);
      if (to_cpu==myid) {
        lcellc = local_cell_coord(cellc);
        to = PTR_VV(cell_array,lcellc,cell_dim);
        INSERT_ATOM(to, input, 0);
      }
      else error("imd_generate: atom on wrong CPU");
#else
      to = PTR_VV(cell_array,cellc,cell_dim);
      INSERT_ATOM(to, input, 0);
#endif
  }
} 

#endif /* TWOD */


#ifndef TWOD

/* initialize for cubic crystals */
void init_cubic(void)
{
  if ((box_param.x==0) || (box_param.y==0) || (box_param.z==0))
    error("box_param not set!");
  box_x.x = box_param.x * box_unit;  box_x.y = 0.0;  box_x.z = 0.0;
  box_y.x = 0.0;  box_y.y = box_param.y * box_unit;  box_y.z = 0.0;
  box_z.x = 0.0;  box_z.y = 0.0;  box_z.z = box_param.z * box_unit;
  make_box();
}

/* generate cubic crystal structures (not just fcc...) */
void generate_fcc(int maxtyp)
{
  ivektor  min, max, cellc, lcellc;
  minicell *to;
  cell     *input;
  real     xx, yy, zz;
  int      to_cpu;
  int      x, y, z, typ;
  int      count;

#ifdef MPI
  if (myid==0)
    if ((box_param.x % cpu_dim.x) || (box_param.y % cpu_dim.y) || 
        (box_param.z % cpu_dim.z))
      error("box_param must be commensurate with cpu_dim");
#endif

  /* conventional unit cell has size box_unit */
  /* FCC, BCC, B2, NaCl, L1_2 */
  if ((maxtyp < 4) || (maxtyp==6) || (maxtyp==7)) {
    box_param.x *= 2; 
    box_param.y *= 2; 
    box_param.z *= 2; 
    box_unit    /= 2;
  }
  /* diamond and zincblende */
  else {
    box_param.x *= 4; 
    box_param.y *= 4; 
    box_param.z *= 4; 
    box_unit    /= 4;
  }

#ifdef BUFCELLS
  min.x =  my_coord.x      * box_param.x / cpu_dim.x;
  max.x = (my_coord.x + 1) * box_param.x / cpu_dim.x;
  min.y =  my_coord.y      * box_param.y / cpu_dim.y;
  max.y = (my_coord.y + 1) * box_param.y / cpu_dim.y;
  min.z =  my_coord.z      * box_param.z / cpu_dim.z;
  max.z = (my_coord.z + 1) * box_param.z / cpu_dim.z;
#else
  min.x = 0; max.x = box_param.x;
  min.y = 0; max.y = box_param.y;
  min.z = 0; max.z = box_param.z;
#endif

#ifdef VEC
  /* estimate number of atoms per CPU, and allocate cell */
  count = (max.x-min.x) * (max.y-min.y) * (max.z-min.z);
  if ((maxtyp==0) || (maxtyp==6) || (maxtyp==7)) count /=2; /* fcc */ 
  else if ((maxtyp==2) || (maxtyp==3))           count /=4; /* bcc */
  else if ((maxtyp==4) || (maxtyp==5))           count /=8; /* diamond */
  count = (int) (count * (1.25 * nallcells / ncells));
  atoms.n = 0;
  atoms.n_max = 0;
  atoms.n_buf = 0;
  alloc_cell(&atoms, count);
#endif

  /* Set up 1 atom input cell */
  input = (cell *) malloc(sizeof(cell));
  if (0==input) error("Cannot allocate input cell.");
  input->n_max = 0;
  alloc_cell(input,1);

  natoms  = 0;
  nactive = 0;

  for (x=min.x ; x<max.x; x++)
    for (y=min.y; y<max.y; y++)
      for (z=min.z; z<max.z; z++) {
 
        typ  = (x+y+z) % 2;

	/* L1_2 structure, cuau3 */
	if (maxtyp == 7) {
	  if (typ==1) 
	    continue;
	  else if (x%2==1 || y%2==1 || z%2==1)
	    typ = 1;
	}

	/* L1_2 structure, cu3au */
	if (maxtyp == 6) {
	  if (typ==1) 
	    continue;
	  else if (x%2==1 || y%2==1 || z%2==1)
	    typ = 0;
	  else
	    typ = 1;
	}

	/* cubic diamond and zincblende case */
	if (maxtyp == 4 || maxtyp == 5) {
	    if ( ((x+y+z)%4==0) && 
		 (z%2==0) && (y%2==0) && (x%2==0)  ) {
		typ=0;
	    }
	    else if ( ((x+y+z)%4==3) && 
		      (z%2==1) && (y%2==1) && (x%2==1) ) {
		if (maxtyp == 4)
		    typ=0;
		else
		    typ=1;
	    }
	    else continue;
	}

        /* B2 == CsCl structure */
        if (maxtyp == 3) {
          if ((z%2==0) && (y%2==0) && (x%2==0)) {
             typ=0;
          }
          else if ((z%2==1) && (y%2==1) && (x%2==1)) {
             typ=1;
          }
          else continue;
        }

        /* bcc case - all atoms get typ 0 */
        if (maxtyp == 2) {
          if (((x%2) != (y%2)) || ((y%2) != (z%2))) continue;
          typ = 0;
        }

        /* maxtyp == 0: fcc;  maxtyp == 1: NaCl */
        if (typ > maxtyp) continue;  /* if fcc, only atoms of type 0 */

	natoms++;
        nactive += 3;

        input->n = 1;
        xx = (x + 0.5) * box_unit;
        yy = (y + 0.5) * box_unit;
        zz = (z + 0.5) * box_unit;
        cellc = cell_coord( xx, yy, zz );
        ORT(input,0,X) = xx;
        ORT(input,0,Y) = yy;
        ORT(input,0,Z) = zz;
#ifndef MONOLJ
        NUMMER(input,0) = natoms;
#ifndef MONO
        SORTE (input,0) = typ;
#endif
        VSORTE(input,0) = typ;
        MASSE (input,0) = masses[typ];
#endif
        num_sort[typ]++;

#ifdef BUFCELLS
	to_cpu = cpu_coord(cellc);
        if (to_cpu==myid) {
	  lcellc = local_cell_coord( cellc );
          to = PTR_VV(cell_array,lcellc,cell_dim);
	  INSERT_ATOM(to, input, 0);
	}
        else error("atom on wrong CPU");
#else
	to = PTR_VV(cell_array,cellc,cell_dim);
        INSERT_ATOM(to, input, 0);
#endif
  }
}

/* generate a cubic Laves structure crystal */

void generate_lav()
{
  minicell *to;
  cell     *input;
  ivektor  min, max, cellc, lcellc;
  int      to_cpu;
  real     px[24],py[24],pz[24];
  real     xx, yy, zz;
  int      i,j,k,l,typ,pa[24];
  int      count;

#ifdef MPI
  if (myid==0)
    if ((box_param.x % cpu_dim.x) || (box_param.y % cpu_dim.y) || 
        (box_param.z % cpu_dim.z))
      error("box_param must be commensurate with cpu_dim");
#endif

  px[ 0] = 0; py[ 0] = 0; pz[ 0] = 0; pa[ 0] = 0;
  px[ 1] = 2; py[ 1] = 2; pz[ 1] = 0; pa[ 1] = 0;
  px[ 2] = 2; py[ 2] = 0; pz[ 2] = 2; pa[ 2] = 0;
  px[ 3] = 0; py[ 3] = 2; pz[ 3] = 2; pa[ 3] = 0;
  px[ 4] = 3; py[ 4] = 3; pz[ 4] = 3; pa[ 4] = 1;
  px[ 5] = 5; py[ 5] = 5; pz[ 5] = 5; pa[ 5] = 1;

  px[ 6] = 4; py[ 6] = 4; pz[ 6] = 0; pa[ 6] = 0;
  px[ 7] = 6; py[ 7] = 6; pz[ 7] = 0; pa[ 7] = 0;
  px[ 8] = 6; py[ 8] = 4; pz[ 8] = 2; pa[ 8] = 0;
  px[ 9] = 4; py[ 9] = 6; pz[ 9] = 2; pa[ 9] = 0;
  px[10] = 7; py[10] = 7; pz[10] = 3; pa[10] = 1;
  px[11] = 1; py[11] = 1; pz[11] = 5; pa[11] = 1;

  px[12] = 4; py[12] = 0; pz[12] = 4; pa[12] = 0;
  px[13] = 6; py[13] = 2; pz[13] = 4; pa[13] = 0;
  px[14] = 6; py[14] = 0; pz[14] = 6; pa[14] = 0;
  px[15] = 4; py[15] = 2; pz[15] = 6; pa[15] = 0;
  px[16] = 7; py[16] = 3; pz[16] = 7; pa[16] = 1;
  px[17] = 1; py[17] = 5; pz[17] = 1; pa[17] = 1;

  px[18] = 0; py[18] = 4; pz[18] = 4; pa[18] = 0;
  px[19] = 2; py[19] = 6; pz[19] = 4; pa[19] = 0;
  px[20] = 2; py[20] = 4; pz[20] = 6; pa[20] = 0;
  px[21] = 0; py[21] = 6; pz[21] = 6; pa[21] = 0;
  px[22] = 3; py[22] = 7; pz[22] = 7; pa[22] = 1;
  px[23] = 5; py[23] = 1; pz[23] = 1; pa[23] = 1;

#ifdef BUFCELLS
  min.x =  my_coord.x      * box_param.x / cpu_dim.x;
  max.x = (my_coord.x + 1) * box_param.x / cpu_dim.x;
  min.y =  my_coord.y      * box_param.y / cpu_dim.y;
  max.y = (my_coord.y + 1) * box_param.y / cpu_dim.y;
  min.z =  my_coord.z      * box_param.z / cpu_dim.z;
  max.z = (my_coord.z + 1) * box_param.z / cpu_dim.z;
#else
  min.x = 0; max.x = box_param.x;
  min.y = 0; max.y = box_param.y;
  min.z = 0; max.z = box_param.z;
#endif

#ifdef VEC
  /* estimate number of atoms per CPU, and allocate cell */
  count = (max.x-min.x) * (max.y-min.y) * (max.z-min.z) * 24;
  count = (int) (count * (1.25 * nallcells / ncells));
  atoms.n = 0;
  atoms.n_max = 0;
  atoms.n_buf = 0;
  alloc_cell(&atoms, count);
#endif

  /* Set up 1 atom input cell */
  input = (cell *) malloc(sizeof(cell));
  if (0==input) error("Cannot allocate input cell.");
  input->n_max = 0;
  alloc_cell(input,1);

  natoms  = 0;
  nactive = 0;

  box_unit /= 8.0;
  for (i=min.x; i<max.x; i++)
    for (j=min.y; j<max.y; j++)
      for (k=min.z; k<max.z; k++) 
	for (l=0; l<24; l++) {

	  natoms++;
          nactive += 3;
	  typ = pa[l];

	  input->n = 1;
	  xx  = (px[l] + 8*i + 0.5) * box_unit;
	  yy  = (py[l] + 8*j + 0.5) * box_unit;
	  zz  = (pz[l] + 8*k + 0.5) * box_unit;
          ORT(input,0,X)  = xx;
	  ORT(input,0,Y)  = yy;
	  ORT(input,0,Z)  = zz;
#ifndef MONOLJ
	  NUMMER(input,0) = natoms;
#ifndef MONO
	  SORTE (input,0) = typ;
#endif
	  VSORTE(input,0) = typ;
          MASSE (input,0) = masses[typ];
#endif
	  cellc = cell_coord(xx,yy,zz);
	  num_sort[typ]++;

#ifdef BUFCELLS
	  to_cpu = cpu_coord(cellc);
	  if (to_cpu==myid) {
	    lcellc = local_cell_coord(cellc);
            to = PTR_VV(cell_array,lcellc,cell_dim);
            INSERT_ATOM(to, input, 0);
          }
          else error("imd_generate: atom on wrong CPU");
#else
          to = PTR_VV(cell_array,cellc,cell_dim);
          INSERT_ATOM(to, input, 0);
#endif
        }
} 

#endif /* not TWOD */

