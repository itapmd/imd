
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
* _fcc      -- generates fcc structure
* _bcc      -- generates fcc structure
* _b2       -- generates B2 structure
* _nacl     -- generates binary nacl structure (atom type 0 and 1)
* _tiqc     -- generates a truncated icosahedra quasicrystal
* _hex      -- generates 2D hexagonal crystal
* _lav      -- generates a cubic Laves structure A15 (MgCu2)
*
* The lattice constant of the crystal structures (fcc and nacl) is 2.0.
*
******************************************************************************/

void generate_atoms(str255 mode)
{
  cell *to;
  int addnumber, tmp, ninc;
  int r, s, t, i, k;
#ifdef MPI
  MPI_Status status;
#endif

  if ((num_sort=calloc(ntypes,sizeof(int)))==NULL) {
      error("cannot allocate memory for num_sort\n");
  }

#ifdef DISLOC
  calc_Epot_ref=1;
#endif

  /* We always have to initialize the velocities on generated structures */
  do_maxwell=1;

#ifdef TWOD
  if (0 == strcmp(mode,"_hex")) {          /* hexagonal crystal */
    init_hex();
    init_cells();
    generate_hex();
#else /* 3D */
  if (0 == strcmp(mode,"_fcc")) {          /* fcc */
    init_cubic();
    init_cells();
    generate_fcc(0);
  } else if (0 == strcmp(mode,"_nacl")) {  /* NaCl */
    init_cubic();
    init_cells();
    generate_fcc(1);
  } else if (0 == strcmp(mode,"_bcc")) {  /* bcc */
    init_cubic();
    init_cells();
    generate_fcc(2);
  } else if (0 == strcmp(mode,"_b2")) {  /* B2 */
    init_cubic();
    init_cells();
    generate_fcc(3);
  } else if (0 == strcmp(mode,"_lav")) {   /* Laves */
    init_cubic();
    init_cells();
    generate_lav();
#ifdef QUASI
  } else if (0 == strcmp(mode,"_tiqc")) {  /* quasicrystal */
    init_qc();
    init_cells();
    generate_qc();
#endif
#endif /* 3D */
  } else if (0==myid) error("Filename with _ specifies unknown structure.");

#ifdef MPI
#ifndef MONOLJ
  /* Get numbering of atoms right even across CPUs */
  ninc = 0;
  if (0==myid)
     MPI_Send( &natoms, 1, MPI_INT, myid + 1, 0, cpugrid );
  else if ((0<myid) && (myid<(num_cpus-1))) {
     MPI_Recv( &tmp,    1, MPI_INT, myid - 1, 0, cpugrid, &status );
     ninc = tmp;
     tmp += natoms;
     MPI_Send( &tmp,    1, MPI_INT, myid + 1, 0, cpugrid );
  }
  else if (myid==(num_cpus-1)) {
     MPI_Recv( &tmp,    1, MPI_INT, myid - 1, 0, cpugrid, &status );
     ninc = tmp;
  }
  /* loop over all atoms, fix numbers */
  for (k=0; k<ncells; ++k) {
    int i;
    cell *p;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) p->nummer[i] += ninc;
  }
#endif /* MONOLJ */

  MPI_Allreduce( &natoms,  &addnumber, 1, MPI_INT, MPI_SUM, cpugrid);
  natoms = addnumber;
  MPI_Allreduce( &nactive, &addnumber, 1, MPI_INT, MPI_SUM, cpugrid);
  nactive = addnumber;
  for (i=0; i<ntypes; i++) {
    MPI_Allreduce(&num_sort[i], &addnumber, 1, MPI_INT, MPI_SUM, cpugrid);
    num_sort[i] = addnumber;
  }
#endif /* MPI */

  if (0== myid) {
    printf("Generated %s structure with %d atoms.\n",mode,natoms);
    addnumber=num_sort[0];
    printf("num_sort = [ %u",num_sort[0]);
    for (i=1; i<ntypes; i++) {
      printf(", %u",num_sort[i]);
      addnumber += num_sort[i];
    }
    printf(" ],  total = %u\n",addnumber);
  }
}


#ifdef TWOD

/* initialize for hexagonal crystal */
void init_hex(void)
{
  if ((box_param.x==0) || (box_param.y==0)) error("box_param not set!");
  box_x.x = box_param.x * sqrt(3.0);
  box_x.y = 0.0;
  box_y.x = 0.0;
  box_y.y = box_param.y;
}

/* generate hexagonal crystal */
void generate_hex()
{
  cell    *input, *to;
  vektor  min, max;
  ivektor cellc;
  int     to_cpu;
  int     i, j, typ;
  real    x, y;

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
  if (0==input) error("Can't allocate input cell.");
  input->n_max = 0;
  alloc_cell(input,1);
  input->masse[0] = 1.0;

  natoms  = 0;
  nactive = 0;

  for (i=min.x ; i<max.x; i++)
    for (j=min.y; j<max.y; j++) {

      typ  = (i+j) % 2;
      if (typ > 0) continue;

      x = (i+0.5) * sqrt(3.0) * 0.5;
      y = (j+0.5) * 0.5;

      natoms++;
      nactive += 2;

      input->n = 1;
      input->ort X(0) = x;
      input->ort Y(0) = y;
      input->nummer[0] = natoms;
      input->sorte[0] = typ;
      num_sort[typ]++;
      cellc = cell_coord(x,y);
#ifdef MPI
      to_cpu = cpu_coord(cellc);
      if (to_cpu==myid) {
        cellc = local_cell_coord(input->ort X(0), input->ort Y(0));
        to = PTR_VV(cell_array,cellc,cell_dim);
        move_atom(to, input, 0);
      }
#else
      to = PTR_VV(cell_array,cellc,cell_dim);
      move_atom(to, input, 0);
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
  box_x.x = box_param.x; box_x.y = 0.0;         box_x.z = 0.0;
  box_y.x = 0.0;         box_y.y = box_param.y; box_y.z = 0.0;
  box_z.x = 0.0;         box_z.y = 0.0;         box_z.z = box_param.z;
}

/* generate fcc or NaCl crystal */
void generate_fcc(int maxtyp)
{
  cell    *input, *to;
  vektor  min, max;
  ivektor cellc;
  int     to_cpu;
  int     x, y, z, typ;

#ifdef MPI
  min.x =  my_coord.x      * box_x.x / cpu_dim.x;
  max.x = (my_coord.x + 1) * box_x.x / cpu_dim.x;
  min.y =  my_coord.y      * box_y.y / cpu_dim.y;
  max.y = (my_coord.y + 1) * box_y.y / cpu_dim.y;
  min.z =  my_coord.z      * box_z.z / cpu_dim.z;
  max.z = (my_coord.z + 1) * box_z.z / cpu_dim.z;
#else
  min.x = 0; max.x = box_x.x;
  min.y = 0; max.y = box_y.y;
  min.z = 0; max.z = box_z.z;
#endif

  /* Set up 1 atom input cell */
  input = (cell *) malloc(sizeof(cell));
  if (0==input) error("Can't allocate input cell.");
  input->n_max = 0;
  alloc_cell(input,1);
#ifndef MONOLJ
  input->masse[0] = 1.0;
#endif

  natoms  = 0;
  nactive = 0;

  for (x=min.x ; x<max.x; x++)
    for (y=min.y; y<max.y; y++)
      for (z=min.z; z<max.z; z++) {
 
        typ  = (x+y+z) % 2;

        if (maxtyp ==3) {
          if ((z%2==0)&&(y%2==0)&&(x%2==0)) {
             typ=0;
          }
          else if ((z%2==1)&&(y%2==1)&&(x%2==1)) {
             typ=1;
          }
          else
            continue;
        }

        /* bcc case - all atoms get typ 0 */
        if (maxtyp == 2) {
          if (((x%2) != (y%2)) || ((y%2) != (z%2))) continue;
          typ = 0;
        }

        if (typ > maxtyp) continue;  /* if fcc, only atoms of type 0 */

	natoms++;
        nactive += 3;

        input->n = 1;
        input->ort X(0) = x + 0.5;
        input->ort Y(0) = y + 0.5;
        input->ort Z(0) = z + 0.5;
        cellc = cell_coord(input->ort X(0), input->ort Y(0), input->ort Z(0));
#ifndef MONOLJ
        input->nummer[0] = natoms;
        input->sorte[0] = typ;
#endif
        num_sort[typ]++;

#ifdef MPI
	to_cpu = cpu_coord(cellc);
        if (to_cpu==myid) {
	  cellc = local_cell_coord(input->ort X(0), input->ort Y(0), 
                                   input->ort Z(0));
          to = PTR_VV(cell_array,cellc,cell_dim);
	  move_atom(to, input, 0);
	}
#else
	to = PTR_VV(cell_array,cellc,cell_dim);
        move_atom(to, input, 0);
#endif
  }
} 

/* generate a cubic Laves structure crystal */

void generate_lav()
{
  cell    *input, *to;
  vektor  min, max, lmin, lmax, rmax, rmin;
  ivektor cellc;
  int     to_cpu;
  real     x, y, z, co;
  real     px[24],py[24],pz[24];
  int     i,j,k,l,typ,pa[24];

  co = sqrt(2.)/8.;

  px[0] = 0; py[0] = 0; pz[0] = 0; pa[0] = 0;
  px[1] = 2; py[1] = 2; pz[1] = 0; pa[1] = 0;
  px[2] = 2; py[2] = 0; pz[2] = 2; pa[2] = 0;
  px[3] = 0; py[3] = 2; pz[3] = 2; pa[3] = 0;
  px[4] = 3; py[4] = 3; pz[4] = 3; pa[4] = 1;
  px[5] = 5; py[5] = 5; pz[5] = 5; pa[5] = 1;

  px[6] = 4; py[6] = 4; pz[6] = 0; pa[6] = 0;
  px[7] = 6; py[7] = 6; pz[7] = 0; pa[7] = 0;
  px[8] = 6; py[8] = 4; pz[8] = 2; pa[8] = 0;
  px[9] = 4; py[9] = 6; pz[9] = 2; pa[9] = 0;
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

#ifdef MPI

  max.x = (int)(box_x.x/sqrt(8.)+0.5);
  max.y = (int)(box_y.y/sqrt(8.)+0.5);
  max.z = (int)(box_z.z/sqrt(8.)+0.5);

  if (myid==0) printf("Zahl der Atome: %d\n",(int)(max.x*max.y*max.z)*24);

  rmin.x = my_coord.x      * box_x.x / cpu_dim.x;
  rmax.x = (my_coord.x + 1) * box_x.x / cpu_dim.x;
  rmin.y = my_coord.y      * box_y.y / cpu_dim.y;
  rmax.y = (my_coord.y + 1) * box_y.y / cpu_dim.y;
  rmin.z = my_coord.z      * box_z.z / cpu_dim.z;
  rmax.z = (my_coord.z + 1) * box_z.z / cpu_dim.z;

  lmin.x = (my_coord.x -1)     * box_x.x / cpu_dim.x;
  lmax.x = (my_coord.x + 2) * box_x.x / cpu_dim.x;
  lmin.y = (my_coord.y -1)     * box_y.y / cpu_dim.y;
  lmax.y = (my_coord.y + 2) * box_y.y / cpu_dim.y;
  lmin.z = (my_coord.z -1)     * box_z.z / cpu_dim.z;
  lmax.z = (my_coord.z + 2) * box_z.z / cpu_dim.z;

  min.x = (int)(lmin.x/sqrt(8.)+0.5);
  max.x = (int)(lmax.x/sqrt(8.)+0.5);
  min.y = (int)(lmin.y/sqrt(8.)+0.5);
  max.y = (int)(lmax.y/sqrt(8.)+0.5);
  min.z = (int)(lmin.z/sqrt(8.)+0.5);
  max.z = (int)(lmax.z/sqrt(8.)+0.5);

  /*
  printf("%d %f %f %f %f %f %f\n",myid,min.x,max.x,min.y,max.y,min.z,max.z);
  */

#else

  min.x=0;
  min.y=0;
  min.z=0;
  max.x=(int)(box_x.x/sqrt(8.)+0.5);
  max.y=(int)(box_y.y/sqrt(8.)+0.5);
  max.z=(int)(box_z.z/sqrt(8.)+0.5);

  printf("Zahl der Atome: %d\n",(int)(max.x*max.y*max.z)*24);

#endif

  /* Set up 1 atom input cell */
  input = (cell *) malloc(sizeof(cell));
  if (0==input) error("Can't allocate input cell.");
  input->n_max = 0;
  alloc_cell(input,1);
#ifndef MONOLJ
  input->masse[0] = 1.0;
#endif

  natoms  = 0;
  nactive = 0;

  for (i=min.x ; i<max.x; i++)
    for (j=min.y; j<max.y; j++)
      for (k=min.z; k<max.z; k++) 

	for (l=0;l<24;l++) {
	  x=(px[l]+8.*i)/sqrt(8.);
	  y=(py[l]+8.*j)/sqrt(8.);
	  z=(pz[l]+8.*k)/sqrt(8.);
	  typ=pa[l];

/*	  if ((x+co < rmin.x) || (x+co > rmax.x) ||
	      (y+co < rmin.y) || (y+co > rmax.y) ||
	      (z+co < rmin.z) || (z+co > rmax.z)) continue; */

	  natoms++;
          nactive += 3;

	  input->n = 1;
	  input->ort X(0) = x + co;
	  input->ort Y(0) = y + co;
	  input->ort Z(0) = z + co;
	  cellc = cell_coord(input->ort X(0),input->ort Y(0),input->ort Z(0));
#ifndef MONOLJ
	  input->nummer[0] = natoms;
	  input->sorte[0] = typ;
#endif
	  num_sort[typ]++;
#ifdef MPI
	  to_cpu = cpu_coord(cellc);
	  if (to_cpu==myid) {
	    cellc = local_cell_coord(input->ort X(0), input->ort Y(0), 
				     input->ort Z(0));
            to = PTR_VV(cell_array,cellc,cell_dim);
            move_atom(to, input, 0);
          }
#else
          to = PTR_VV(cell_array,cellc,cell_dim);
          move_atom(to, input, 0);
#endif
        }
} 


#endif /* not TWOD */


