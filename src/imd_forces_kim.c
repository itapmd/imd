/****************************************************************
 *
 * IMD -- The ITAP Molecular Dynamics Program
 *
 * Copyright 1996-2012 Institute for Theoretical and Applied Physics,
 * University of Stuttgart, D-70550 Stuttgart
 *
 ****************************************************************/

/****************************************************************
 *
 * imd_forces_kim.c -- call force loop from openKIM.org project
 *
 ****************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#ifdef KIM

#include "imd.h"

#include "KIM_API_C.h"
#include "KIM_API_status.h"

#define NBLMINLEN 100000
#define KIM_NBCELLS 27

/****************************************************************
 *
 * global variables
 *
 ****************************************************************/

/* pointers for IMD neighbor lists */
int  *tl = NULL, *tb = NULL, *cl_off = NULL, *cl_num = NULL, nb_max = 0;
/* pointer for KIM neighbor list */
int  *kim_nbl = NULL;
/* pointer to relative position vector */
real *Rij = NULL;

/****************************************************************
 *
 * test_buffer
 *
 ****************************************************************/

/* test buffer */
typedef struct {
  /* the cell pointer for the current cell */
  cell *p;

  /* contiguous memory for the KIM API data structures */
  int   max_len;		/* length of following arrays in units of atoms */
  real *coords;			/* coordinates = 3*max_len */
  real *forces;			/* forces = 3*max_len */
  real *penergy;		/* particleEnergy = max_len */
  int  *ptypes;			/* particleTypes = max_len */
} test_buffer_t;

/****************************************************************
 *
 * function declarations
 *
 ****************************************************************/

/* function declarations */
void  init_kim_object(imd_kim_t *);
void  init_kim_info();
void  init_kim();
void  destroy_kim();
void  write_descriptor(char **);
void  model_has_flags();
void  deallocate_nblist();
int   estimate_nblist_size();
void  make_nblist();
void  calc_forces(int);
void  check_nblist();
int   get_neigh(void **, int *, int *, int *, int *, int **, double **);
void  extend_test_buffer(int);
void  kim_warning(char *);
void  kim_error(char *);
int imd_process_dEdr(void **, double *, double *, double **, int *, int *);

/****************************************************************
 *
 *  void init_kim_object(imd_kim_t *);
 *
 *    properly initialize all variables
 *
 ****************************************************************/

void init_kim_object(imd_kim_t *kim)
{
  kim->pkim = NULL;

  kim->model_using_half = 0;
  kim->model_using_cluster = 0;
  kim->model_using_Rij = 0;

  kim->model_has_forces = 0;
  kim->model_has_energy = 0;
  kim->model_has_particleEnergy = 0;
  kim->model_has_process_dEdr = 0;

  kim->ind_coordinates = 0;
  kim->ind_numberOfParticles = 0;
  kim->ind_numberContributingParticles = 0;
  kim->ind_numberParticleTypes = 0;
  kim->ind_particleTypes = 0;
  kim->ind_get_neigh = 0;
  kim->ind_neighObject = 0;
  kim->ind_cutoff = 0;
  kim->ind_forces = 0;
  kim->ind_energy = 0;
  kim->ind_particleEnergy = 0;
  kim->ind_process_dEdr = 0;

  kim->cell_ind = 0;
  kim->cell_list = NULL;
  kim->cell_offset = NULL;
  kim->cell_atom_ind = 0;
  kim->cell_index_atom = NULL;

  kim->kim_particle_codes = NULL;
  kim->iterator_position = 0;
}

/****************************************************************
 *
 *  void destroy_kim();
 *
 *    this is the routine for deallocating the kim object
 *    it is directly called from the main imd.c routine
 *
 ****************************************************************/

void destroy_kim()
{
  int   kimerror = 0;
  test_buffer_t *test_buffer = NULL;

  test_buffer = (test_buffer_t *) KIM_API_get_test_buffer(kim.pkim, &kimerror);
  free(test_buffer->coords);
  free(test_buffer->forces);
  free(test_buffer->penergy);
  free(test_buffer->ptypes);
  free(test_buffer);

  kimerror = KIM_API_model_destroy(kim.pkim);

  KIM_API_free(&kim.pkim, &kimerror);

  free(kim.kim_particle_codes);
  free(kim.cell_index_atom);
  free(kim.cell_offset);
  free(kim.cell_list);

  free(tl);
  free(tb);
  free(cl_off);
  free(cl_num);
  free(kim_nbl);
  if (kim.model_using_Rij)
    free(Rij);
}

/****************************************************************
 *
 *  void init_kim_info();
 *    this function is called from the setup_potentials routine
 *    its purpose is to check for the flags of the given model
 *    also read the cutoff, which is required for reading atoms into cells
 *
 ****************************************************************/

void init_kim_info()
{
  char *test_descriptor_string = NULL;
  void *pkim;
  int   kimerror;
  real  rc = 0.0;

  /* initialize the variables in the kim struct */
  init_kim_object(&kim);

  /* print info about KIM interface */
  if ((0 == myid) && (0 == myrank)) {
    printf("Setting up KIM interface with model %s.\n", kim_model_name);
    fflush(stdout);
    kim_warning("THIS INTERFACE IS STILL IN DEVELOPMENT!!!! DO NOT USE !!!!");
  }

  /* set up a bogus kim descriptor file */
  write_descriptor(&test_descriptor_string);

  /* get KIM API object representing the KIM Model only */
  kimerror = KIM_API_string_init(&pkim, test_descriptor_string, kim_model_name);
  if (KIM_STATUS_OK > kimerror) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM initialization failed", kimerror);
    exit(EXIT_FAILURE);
  } else {
    free(test_descriptor_string);
    test_descriptor_string = NULL;
  }

  /* set the cutoff pointer to the correct location */
  kimerror = KIM_API_set_data(pkim, "cutoff", 1, (void *)&rc);
  if (KIM_STATUS_OK > kimerror) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_set_data failed", kimerror);
    exit(EXIT_FAILURE);
  }

  /* determine if we have half neighbor lists */
  kim.model_using_half = KIM_API_is_half_neighbors(pkim, &kimerror);
  if (KIM_STATUS_OK > kimerror) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_is_half_neighbors failed", kimerror);
    exit(EXIT_FAILURE);
  }

  /* ask about the NBC method and write cluster and Rij flags */
  char *NBC_method = (char *)KIM_API_get_NBC_method(pkim, &kimerror);
  if (KIM_STATUS_OK > kimerror) {
    KIM_API_report_error(__LINE__, __FILE__, "NBC method not set", kimerror);
    exit(EXIT_FAILURE);
  }
  /* check for CLUSTER mode */
  kim.model_using_cluster = (strcmp(NBC_method, "CLUSTER") == 0);
  /* check if Rij needed for get_neigh */
  kim.model_using_Rij = ((strcmp(NBC_method, "NEIGH_RVEC_F") == 0) || (strcmp(NBC_method, "NEIGH_RVEC_H") == 0));
  free((void *)NBC_method);

  /* initialize the model */
  kimerror = KIM_API_model_init(pkim);
  if (KIM_STATUS_OK > kimerror) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_init failed", kimerror);
    exit(EXIT_FAILURE);
  }

  /* set the cutoff */
  cellsz = SQR(rc);

  kimerror = KIM_API_model_destroy(pkim);
  if (KIM_STATUS_OK > kimerror) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_destroy failed", kimerror);
    exit(EXIT_FAILURE);
  }

  /* deallocate the memory */
  KIM_API_free(&pkim, &kimerror);
  if (KIM_STATUS_OK > kimerror) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_free failed", kimerror);
    exit(EXIT_FAILURE);
  }
}

/****************************************************************
 *
 *  void model_has_flags();
 *
 *    this function determines the capabilities of the requested KIM model
 *    currently we use the LAMMPS philosophy:
 *      >> If we can run anything, then we do it. <<
 *    a warning is printed at the beginning if any key capabilities are missing
 *    and the simulation will probably not work properly
 *    in some special cases this however might be desired
 *
 ****************************************************************/

void model_has_flags()
{
  void *pkim;
  int   kimerror;

  /* get KIM API object representing the KIM Model only */
  kimerror = KIM_API_model_info(&pkim, kim_model_name);
  if (KIM_STATUS_OK > kimerror) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM initialization failed.", kimerror);
    exit(EXIT_FAILURE);
  }

  /* determine if the KIM Model can compute the total energy */
  KIM_API_get_index(pkim, (char *)"energy", &kimerror);
  kim.model_has_energy = (kimerror == KIM_STATUS_OK);
  if (!kim.model_has_energy)
    kim_warning("KIM Model does not provide 'energy' - the potential energy will be zero!");

  /* determine if the KIM Model can compute the forces */
  KIM_API_get_index(pkim, (char *)"forces", &kimerror);
  kim.model_has_forces = (kimerror == KIM_STATUS_OK);
  if (!kim.model_has_forces)
    kim_warning("KIM Model does not provide 'forces' -  the forces will be zero!");

  /* determine if the KIM Model can compute the particleEnergy */
  KIM_API_get_index(pkim, (char *)"particleEnergy", &kimerror);
  kim.model_has_particleEnergy = (kimerror == KIM_STATUS_OK);
  if (!kim.model_has_particleEnergy)
    kim_warning("KIM Model does not provide 'particleEnergy' - the energy per atom will be zero!");

  /* determine if the KIM Model uses the process_dEdr routine */
  KIM_API_get_index(pkim, (char *)"process_dEdr", &kimerror);
  kim.model_has_process_dEdr = (kimerror == KIM_STATUS_OK);
  if (!kim.model_has_process_dEdr)
    kim_warning("KIM Model does not use the 'process_dEdr' routine - the virial and pressure will be zero!");

  /* tear down KIM API object */
  KIM_API_free(&pkim, &kimerror);
  if (KIM_STATUS_OK > kimerror) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_free failed.", kimerror);
    exit(EXIT_FAILURE);
  }
}

/****************************************************************
 *
 *  void init_kim();
 *
 *    this is the actual initialization of the kim object
 *    it is directly called from the main imd.c routine
 *
 ****************************************************************/

void init_kim()
{
  char *test_descriptor_string = NULL;
  int   i, kimerror = 0;
  real  rc = 0.0;

  /* determine the capabilities of the requested model */
  model_has_flags();

  /* generate KIM file in a string */
  write_descriptor(&test_descriptor_string);

  /* initialize the KIM object from the string */
  kimerror = KIM_API_string_init(&kim.pkim, test_descriptor_string, kim_model_name);
  if (KIM_STATUS_OK > kimerror)
    KIM_API_report_error(__LINE__, __FILE__, "KIM initialization failed", kimerror);
  else {
    free(test_descriptor_string);
    test_descriptor_string = NULL;
  }

  /* set the pointer for the cutoff */
  kimerror = KIM_API_set_data(kim.pkim, "cutoff", 1, (void *)&rc);
  if (KIM_STATUS_OK > kimerror) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_set_data failed", kimerror);
    exit(EXIT_FAILURE);
  }

  /* get correct index of each variable in kim_api object */
  /* *INDENT-OFF* */
  KIM_API_getm_index(kim.pkim, &kimerror, 3 * 12,
    "coordinates", 			&kim.ind_coordinates, 			1,
    "cutoff", 				&kim.ind_cutoff, 			1,
    "numberOfParticles", 		&kim.ind_numberOfParticles, 		1,
    "numberParticleTypes", 		&kim.ind_numberParticleTypes, 		1,
    "particleTypes", 			&kim.ind_particleTypes, 		1,
    "numberContributingParticles", 	&kim.ind_numberContributingParticles, 	kim.model_using_half,
    "forces", 				&kim.ind_forces, 			kim.model_has_forces,
    "energy", 				&kim.ind_energy, 			kim.model_has_energy,
    "particleEnergy", 			&kim.ind_particleEnergy, 		kim.model_has_particleEnergy,
    "neighObject", 			&kim.ind_neighObject, 			!kim.model_using_cluster,
    "get_neigh", 			&kim.ind_get_neigh, 			!kim.model_using_cluster,
    "process_dEdr",			&kim.ind_process_dEdr, 			kim.model_has_process_dEdr);
  /* *INDENT-ON* */
  if (KIM_STATUS_OK > kimerror) {
    KIM_API_report_error(__LINE__, __FILE__, "getm_index", kimerror);
    exit(EXIT_FAILURE);
  }

  /* initialize the KIM API */
  kimerror = KIM_API_model_init(kim.pkim);
  if (KIM_STATUS_OK > kimerror) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_model_init", kimerror);
    exit(EXIT_FAILURE);
  }

  /* initialize the test buffer */
#define KIM_MAX_LEN 10000
  test_buffer_t *test_buffer = NULL;
  test_buffer = malloc(1 * sizeof(test_buffer_t));
  test_buffer->coords = malloc(3 * KIM_MAX_LEN * sizeof(real));
  test_buffer->forces = malloc(3 * KIM_MAX_LEN * sizeof(real));
  test_buffer->penergy = malloc(KIM_MAX_LEN * sizeof(real));
  test_buffer->ptypes = malloc(KIM_MAX_LEN * sizeof(int));
  if (NULL == test_buffer || NULL == test_buffer->coords || NULL == test_buffer->forces ||
	NULL == test_buffer->penergy || NULL == test_buffer->ptypes) {
    kim_error("Could not allocate memory for the test_buffer");
  }
  for (i = 0; i < KIM_MAX_LEN; i++) {
    test_buffer->coords[3 * i + 0] = 0.0;
    test_buffer->coords[3 * i + 1] = 0.0;
    test_buffer->coords[3 * i + 2] = 0.0;
    test_buffer->forces[3 * i + 0] = 0.0;
    test_buffer->forces[3 * i + 1] = 0.0;
    test_buffer->forces[3 * i + 2] = 0.0;
    test_buffer->penergy[i] = 0.0;
    test_buffer->ptypes[i] = -1;
  }
  test_buffer->max_len = KIM_MAX_LEN;
  kim.cell_index_atom = (int *)malloc(KIM_MAX_LEN * sizeof(int));
#undef KIM_MAX_LEN

  /* register the test buffer with the API */
  KIM_API_set_test_buffer(kim.pkim, test_buffer, &kimerror);
  if (KIM_STATUS_OK > kimerror) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_set_test_buffer", kimerror);
    exit(EXIT_FAILURE);
  }

  /* setup mapping between IMD and KIM particle type codes */
  kim.kim_particle_codes = malloc(ntypes * sizeof(int));
  for (i = 0; i < ntypes; i++) {
    kim.kim_particle_codes[i] = KIM_API_get_partcl_type_code(kim.pkim, kim_el_names[i], &kimerror);
    if (KIM_STATUS_OK > kimerror) {
      KIM_API_report_error(__LINE__, __FILE__, "create_kim_particle_codes: symbol not found ", kimerror);
      exit(EXIT_FAILURE);
    }
  }

  /* initialize the kim.cell_list and kim.cell_offset */
  kim.cell_list = malloc(nallcells * sizeof(int));
  kim.cell_offset = malloc(nallcells * sizeof(int));
  if (NULL == kim.cell_list || NULL == kim.cell_offset)
    kim_error("Could not allocate memory for kim.cells");
  for (i = 0; i < nallcells; i++) {
    kim.cell_list[i] = -1;
    kim.cell_offset[i] = -1;
  }
}

/****************************************************************
 *
 *  void write_descriptor(char *);
 *
 *    create a string that can be used as a virtual KIM descriptor file
 *
 ****************************************************************/

void write_descriptor(char **string)
{
  int   i = 0;
  char  tmp[100];

  if (*string != NULL)
    error("test_descriptor_string already allocated.");
  *string = (char *)malloc(100 * 75 * sizeof(char));
  if (NULL == *string)
    kim_error("Could not allocate memory for kim descriptor string");
  strcpy(*string, "");

  /* *INDENT-OFF* */
  strcat(*string,
    "# This file is automatically generated from IMD with openKIM.org support\n"
    "TEST_NAME        := test_IMD\n\n");

  strcat(*string,
    "# Base units\n"
    "Unit_length      := A\n"
    "Unit_energy      := eV\n"
    "Unit_charge      := e\n"
    "Unit_temperature := K\n"
    "Unit_time        := ps\n\n");

  strcat(*string,
    "SUPPORTED_ATOM/PARTICLES_TYPES:\n"
    "# Symbol/name           Type            code");

  for (i = 0; i < vtypes; i++) {
    sprintf(tmp, "\n%-24s%-16s%-3i", kim_el_names[i], "spec", i);
    strcat(*string, tmp);
  }

  strcat(*string,
    "\n\nCONVENTIONS:\n"
    "# Name 			Type\n"
    "ZeroBasedLists 		flag\n"
    "Neigh_BothAccess 		flag\n\n");

  strcat(*string,
    "NEIGH_RVEC_H 	     flag\n"
    "NEIGH_PURE_H            flag\n"
    "NEIGH_RVEC_F            flag\n\n"
    "NEIGH_PURE_F            flag\n\n");

  /* Write input section */
  strcat(*string,
    "MODEL_INPUT:\n"
    "# Name                         Type         Unit    Shape\n"
    "numberOfParticles              integer      none    []\n"
    "numberContributingParticles    integer      none    []\n"
    "numberParticleTypes            integer      none    []\n"
    "particleTypes                  integer      none    [numberOfParticles]\n"
    "coordinates                    double       length  [numberOfParticles,3]\n"
    "neighObject                    pointer      none    []\n"
    "process_dEdr 		    method 	 none    []\n"
    "get_neigh                      method       none    []\n\n");

  /* Write output section */
  strcat(*string,
    "MODEL_OUPUT:\n"
    "# Name                         Type         Unit    Shape\n"
    "compute                        method       none    []\n"
    "destroy                        method       none    []\n"
    "cutoff                         double       length  []\n");
  if (kim.model_has_energy)
    strcat(*string, "energy                         double       energy  []\n\n");
  if (kim.model_has_forces)
    strcat(*string, "forces                         double       force   [numberOfParticles,3]\n\n");
  if (kim.model_has_particleEnergy)
    strcat(*string, "particleEnergy                 double       energy  [numberOfParticles]\n\n");
  /* *INDENT-ON* */
}


/****************************************************************
 *
 *  void deallocate_nblist();
 *    deallocate (largest part of) neighbor list
 *
 ****************************************************************/

void deallocate_nblist()
{
#if defined(DEBUG) || defined(TIMING)
  if (myid == 0)
    printf("Size of neighbor table: %d MB\n", (int)(nb_max * sizeof(int) / SQR(1024)));
#endif
  if (tb)
    free(tb);
  if (kim_nbl)
    free(kim_nbl);
  if (Rij)
    free(Rij);
  tb = NULL;
  kim_nbl = NULL;
  Rij = NULL;
  have_valid_nbl = 0;
}

/****************************************************************
 *
 *  int estimate_nblist_size();
 *    estimate_nblist_size
 *
 ****************************************************************/

int estimate_nblist_size()
{
  int   c, tn = 1;

  /* for all cells */
  for (c = 0; c < ncells2; c++) {

    int   i, c1 = cnbrs[c].np;
    cell *p = cell_array + c1;

    /* for each atom in cell */
    for (i = 0; i < p->n; i++) {

      int   m;
      vektor d1;

      d1.x = ORT(p, i, X);
      d1.y = ORT(p, i, Y);
#ifndef TWOD
      d1.z = ORT(p, i, Z);
#endif

      /* for each neighboring atom */
      for (m = 0; m < KIM_NBCELLS; m++) {
	int   c2, jstart, j;
	real  r2;
	cell *q;
	c2 = cnbrs[c].nq[m];
	if (c2 < 0)
	  continue;
	if (c2 == c1)
	  jstart = i + 1;
	else
	  jstart = 0;
	q = cell_array + c2;
#ifdef ia64
#pragma ivdep
#endif
	for (j = jstart; j < q->n; j++) {
	  vektor d;
	  d.x = ORT(q, j, X) - d1.x;
	  d.y = ORT(q, j, Y) - d1.y;
#ifndef TWOD
	  d.z = ORT(q, j, Z) - d1.z;
#endif
	  r2 = SPROD(d, d);
	  if (r2 < cellsz)
	    tn++;
	}
      }
    }
  }
  return tn;
}

/****************************************************************
 *
 *  void make_nblist();
 *    create the neighbor list for all cells
 *
 ****************************************************************/

void make_nblist()
{
  static int at_max = 0, pa_max = 0, ncell_max = 0;
  int   c, i, k, n, tn, at, cc;

  /* update reference positions */
  for (k = 0; k < ncells; k++) {
    cell *p = cell_array + cnbrs[k].np;
#ifdef ia64
#pragma ivdep,swp
#endif
    for (i = 0; i < p->n; i++) {
      NBL_POS(p, i, X) = ORT(p, i, X);
      NBL_POS(p, i, Y) = ORT(p, i, Y);
#ifndef TWOD
      NBL_POS(p, i, Z) = ORT(p, i, Z);
#endif
    }
  }

  /* (re)allocate cl_off */
  if (nallcells > ncell_max) {
    cl_off = (int *)realloc(cl_off, nallcells * sizeof(int));
    if (cl_off == NULL)
      error("cannot allocate neighbor table");
    ncell_max = nallcells;
  }

  /* count atom numbers (including buffer atoms) */
  at = 0;
  for (k = 0; k < nallcells; k++) {
    cell *p = cell_array + k;
    cl_off[k] = at;
    at += p->n;
  }

  /* (re-)allocate neighbor table */
  if (at >= at_max) {
    free(tl);
    free(cl_num);
    at_max = (int)(1.1 * at);
    tl = (int *)malloc(at_max * sizeof(int));
    cl_num = (int *)malloc(at_max * sizeof(int));
  }
  if (NULL == tb && NULL == kim_nbl) {
    if (0 == last_nbl_len)
      nb_max = MAX(((int)(nbl_size * estimate_nblist_size())), (NBLMINLEN));
    else
      nb_max = (int)(nbl_size * last_nbl_len);
    tb = (int *)malloc(nb_max * sizeof(int));
    kim_nbl = (int *)malloc(nb_max * sizeof(int));
    if (kim.model_using_Rij == 1 && NULL == Rij)
      Rij = (real *)malloc(3 * nb_max * sizeof(real));
  } else if (last_nbl_len * sqrt(nbl_size) > nb_max) {
    free(tb);
    free(kim_nbl);
    nb_max = (int)(nbl_size * last_nbl_len);
    tb = (int *)malloc(nb_max * sizeof(int));
    kim_nbl = (int *)malloc(nb_max * sizeof(int));
    if (kim.model_using_Rij == 1) {
      free(Rij);
      Rij = (real *)malloc(3 * nb_max * sizeof(real));
    }
  }
  if ((tl == NULL) || (tb == NULL) || (cl_num == NULL) || (kim_nbl == NULL))
    error("cannot allocate neighbor table");

  /* set cl_num */
  n = 0;
  for (k = 0; k < nallcells; k++) {
    cell *p = cell_array + k;
    for (i = 0; i < p->n; i++)
      cl_num[n++] = k;
  }

  /* for all cells */
  n = 0;
  tn = 0;
  tl[0] = 0;
  for (c = 0; c < ncells2; c++) {

    int   c1 = cnbrs[c].np;
    cell *p = cell_array + c1;

    /* for each atom in cell */
    for (i = 0; i < p->n; i++) {

      int   m;
      vektor d1;

      d1.x = ORT(p, i, X);
      d1.y = ORT(p, i, Y);
#ifndef TWOD
      d1.z = ORT(p, i, Z);
#endif

      /* for each neighboring atom */
      for (m = 0; m < KIM_NBCELLS; m++) {
	int   c2, jstart, j;
	cell *q;
	c2 = cnbrs[c].nq[m];
	if (c2 < 0)
	  continue;
	if (c2 == c1)
	  if (kim.model_using_half)
	    jstart = i + 1;
	  else
	    jstart = 0;
	else
	  jstart = 0;
	q = cell_array + c2;
#ifdef ia64
#pragma ivdep
#endif
	for (j = jstart; j < q->n; j++) {
	  if (c2 != c1 || j != i) {
	    vektor d;
	    real  r2;
	    d.x = ORT(q, j, X) - d1.x;
	    d.y = ORT(q, j, Y) - d1.y;
#ifndef TWOD
	    d.z = ORT(q, j, Z) - d1.z;
#endif
	    r2 = SPROD(d, d);
	    if (r2 < cellsz) {
	      tb[tn++] = cl_off[c2] + j;
	    }
	  }
	}
      }
      pa_max = MAX(pa_max, tn - tl[n]);
      tl[++n] = tn;
      if (tn > nb_max - 2 * pa_max) {
	error("neighbor table full - increase nbl_size");
      }
    }
  }
  last_nbl_len = tn;
  have_valid_nbl = 1;
  nbl_count++;
}

/****************************************************************
 *
 *  void calc_forces(int);
 *    this is the routine that calculates the forces on all atoms
 *    we prepare the data in the KIM format and call the API
 *    to do the actual force computation
 *
 ****************************************************************/

void calc_forces(int steps)
{
  void *pkim = kim.pkim;
  int   i, k, n = 0;
  int   kimerror = 0;
  test_buffer_t *buf;

  /* static arrays for data exchange with KIM */
  buf = (test_buffer_t *) KIM_API_get_test_buffer(pkim, &kimerror);
  if (KIM_STATUS_OK > kimerror) {
    KIM_API_report_error(__LINE__, __FILE__, "Could not get test buffer", kimerror);
    exit(EXIT_FAILURE);
  }

  if (0 == have_valid_nbl) {
#ifdef MPI
    /* check message buffer size */
    if (0 == nbl_count % BUFSTEP)
      setup_buffers();
#endif
    /* update cell decomposition */
    fix_cells();
  }

  /* fill the buffer cells */
  send_cells(copy_cell, pack_cell, unpack_cell);

  /* make new neighbor lists */
  if (0 == have_valid_nbl)
    make_nblist();

  /* clear global accumulation variables */
  tot_pot_energy = 0.0;
  virial = 0.0;
  vir_xx = 0.0;
  vir_yy = 0.0;
  vir_xy = 0.0;
  vir_zz = 0.0;
  vir_yz = 0.0;
  vir_zx = 0.0;
  nfc++;

  /* clear per atom accumulation variables, also in buffer cells */
  n = 0;
  for (k = 0; k < nallcells; k++) {
    cell *p = cell_array + k;
    n += p->n;
#ifdef ia64
#pragma ivdep,swp
#endif
    for (i = 0; i < p->n; i++) {
      KRAFT(p, i, X) = 0.0;
      KRAFT(p, i, Y) = 0.0;
      KRAFT(p, i, Z) = 0.0;
#ifdef STRESS_TENS
      PRESSTENS(p, i, xx) = 0.0;
      PRESSTENS(p, i, yy) = 0.0;
      PRESSTENS(p, i, xy) = 0.0;
      PRESSTENS(p, i, zz) = 0.0;
      PRESSTENS(p, i, yz) = 0.0;
      PRESSTENS(p, i, zx) = 0.0;
#endif
      POTENG(p, i) = 0.0;
#ifdef CNA
      if (cna)
	MARK(p, i) = 0;
#endif
    }
  }

  if (n * 1.1 > (*buf).max_len)
    extend_test_buffer((int)(n * KIM_NBCELLS * 1.1));

  /* interactions for all atoms */
  n = 0;
  kim.cell_atom_ind = 0;
  for (k = 0; k < ncells; k++) {
    int   c2, ind = 0;
    cell *p = cell_array + cnbrs[k].np;
    cell *q;

    /* convert data into KIM format */
    for (n = 0; n < KIM_NBCELLS; n++) {
      c2 = cnbrs[k].nq[n];
      kim.cell_list[c2] = n;
      kim.cell_offset[n] = ind;
      q = cell_array + c2;
      memcpy(buf->coords + 3 * ind, q->ort, 3 * q->n * sizeof(real));
      for (i = 0; i < q->n; i++) {
	(buf->ptypes + ind)[i] = kim.kim_particle_codes[SORTE(q, i)];
	kim.cell_index_atom[ind + i] = c2;
      }
      ind += q->n;
    }

    /* set pointers to correct location */
    real  pot = 0.0;

    /* *INDENT-OFF* */
    KIM_API_setm_data_by_index(pkim, &kimerror, 4 * 10,
      kim.ind_coordinates, 			3 * ind, 	(void *)buf->coords, 		1,
      kim.ind_particleTypes, 			ind, 		(void *)buf->ptypes, 		1,
      kim.ind_forces, 				3 * ind, 	(void *)buf->forces, 		kim.model_has_forces,
      kim.ind_energy, 				1, 		(void *)&pot, 			kim.model_has_energy,
      kim.ind_particleEnergy, 			ind, 		(void *)buf->penergy, 		kim.model_has_particleEnergy,
      kim.ind_numberOfParticles, 		1, 		(void *)&ind,			1,
      kim.ind_numberContributingParticles, 	1, 		(void *)&p->n,			1,
      kim.ind_numberParticleTypes, 		1, 		(void *)&ntypes, 		1,
      kim.ind_process_dEdr, 			1, 		(void *)imd_process_dEdr, 	kim.model_has_process_dEdr,
      kim.ind_get_neigh, 			1, 		(void *)get_neigh, 		1);
    /* *INDENT-ON* */
    if (KIM_STATUS_OK > kimerror) {
      KIM_API_report_error(__LINE__, __FILE__, "Error setting data by index", kimerror);
      exit(EXIT_FAILURE);
    }

    kim.cell_ind = cnbrs[k].nq[0];

    KIM_API_set_compute_by_index(pkim, kim.ind_particleEnergy, 1, &kimerror);
    if (KIM_STATUS_OK > kimerror) {
      KIM_API_report_error(__LINE__, __FILE__, "Error setting compute by index", kimerror);
      exit(EXIT_FAILURE);
    }

    /* call the KIM API to do the actual force calculation */
    kimerror = KIM_API_model_compute(pkim);
    if (KIM_STATUS_OK > kimerror) {
      KIM_API_report_error(__LINE__, __FILE__, "Error in KIM_API_model_compute", kimerror);
      exit(EXIT_FAILURE);
    }

    /* read all data from kim_arrays */
    /* all data has to be added and not only copied ! */
    ind = 0;
    for (n = 0; n < KIM_NBCELLS; n++) {
      c2 = cnbrs[k].nq[n];
      kim.cell_list[c2] = -1;
      kim.cell_offset[c2] = -1;
      q = cell_array + c2;
      for (i = 0; i < q->n; i++) {
	KRAFT(q, i, X) += (buf->forces + 3 * ind)[3 * i + 0];
	KRAFT(q, i, Y) += (buf->forces + 3 * ind)[3 * i + 1];
	KRAFT(q, i, Z) += (buf->forces + 3 * ind)[3 * i + 2];
	if (c2 == cnbrs[k].np)
	  POTENG(q, i) += (buf->penergy + ind)[i];
      }
      ind += q->n;
    }
    tot_pot_energy += pot;
    kim.cell_atom_ind += p->n;
  }

  /* add forces back to original cells/cpus */
  send_forces(add_forces, pack_forces, unpack_forces);
}

/****************************************************************
 *
 *  check_nblist
 *
 ****************************************************************/

void check_nblist()
{
  real  r2, max1 = 0.0, max2;
  vektor d;
  int   k;

  /* compare with reference positions */
  for (k = 0; k < NCELLS; k++) {
    int   i;
    cell *p = CELLPTR(k);
#ifdef ia64
#pragma ivdep,swp
#endif
    for (i = 0; i < p->n; i++) {
      d.x = ORT(p, i, X) - NBL_POS(p, i, X);
      d.y = ORT(p, i, Y) - NBL_POS(p, i, Y);
#ifndef TWOD
      d.z = ORT(p, i, Z) - NBL_POS(p, i, Z);
#endif
      r2 = SPROD(d, d);
      if (r2 > max1)
	max1 = r2;
    }
  }

#ifdef MPI
  MPI_Allreduce(&max1, &max2, 1, REAL, MPI_MAX, cpugrid);
#else
  max2 = max1;
#endif
  if (max2 > SQR(0.5 * nbl_margin))
    have_valid_nbl = 0;
}

/****************************************************************
 *
 *  get_neigh
 *    This is the neighbor function for the KIM API.
 *    It translates the IMD neighbors into KIM neighbors.
 *
 ****************************************************************/

int get_neigh(void **kimmdl, int *mode, int *request, int *atom, int *numnei, int **nei1atom, double **pRij)
{
  int   i, imd_nr;
  int   ii, in_cell;
  int   offset, atom_offset;

  /* determine on which cell we are working */
  cell *p = cell_array + kim.cell_ind;
  cell *q;

  vektor d1;

  /* iterator mode */
  if (*mode == 0) {
    /* increment iterator */
    if (*request == 1) {
      if (kim.iterator_position < p->n) {
	/* get global atom index */
	imd_nr = kim.cell_atom_ind + kim.iterator_position;
	/* get number of neighbors */
	*numnei = tl[imd_nr + 1] - tl[imd_nr];
	/* store position of the current atom */
	if (kim.model_using_Rij) {
	  d1.x = ORT(p, kim.iterator_position, X);
	  d1.y = ORT(p, kim.iterator_position, Y);
	  d1.z = ORT(p, kim.iterator_position, Z);
	}
	/* create neighbor list */
	for (i = 0; i < *numnei; i++) {
	  /* get atom index */
	  ii = tb[tl[imd_nr] + i];
	  /* get cell index */
	  in_cell = kim.cell_list[cl_num[ii]];
	  /* check if cell is in the test buffer */
	  if (in_cell != -1) {
	    /* offset of the cell in the test buffer */
	    offset = kim.cell_offset[in_cell];
	    /* offset of the atom in the cell */
	    atom_offset = ii - cl_off[cl_num[ii]];
	    kim_nbl[i] = offset + atom_offset;
	    /* calculate relative postition vector if required */
	    if (kim.model_using_Rij) {
	      q = cell_array + cl_num[ii];
	      Rij[3 * i + 0] = ORT(q, atom_offset, X) - d1.x;
	      Rij[3 * i + 1] = ORT(q, atom_offset, Y) - d1.y;
	      Rij[3 * i + 2] = ORT(q, atom_offset, Z) - d1.z;
	    }
	  }
	  /* if this is not the case, then something went wrong */
	  else {
	    kim_error("Cell index error");
	  }
	}
	*atom = kim.iterator_position++;
	*nei1atom = kim_nbl;
	if (kim.model_using_Rij)
	  *pRij = Rij;
	return KIM_STATUS_OK;
      } else if (kim.iterator_position == p->n) {
	*numnei = 0;
	return KIM_STATUS_NEIGH_ITER_PAST_END;
      } else if (kim.iterator_position > p->n) {
	error("KIM neighbor iterator exceeded range");
      }
    } else
      /* reset iterator */
    if (*request == 0) {
      kim.iterator_position = 0;
      *numnei = 0;
      return KIM_STATUS_NEIGH_ITER_INIT_OK;
    }
  } else
    /* locator mode */
  if (*mode == 1) {
    if (*request < p->n && *request >= 0) {
      /* get global atom index */
      imd_nr = kim.cell_atom_ind + *request;
      /* get number of neighbors */
      *numnei = tl[imd_nr + 1] - tl[imd_nr];
      /* store position of the current atom */
      if (kim.model_using_Rij) {
	d1.x = ORT(p, *request, X);
	d1.y = ORT(p, *request, Y);
	d1.z = ORT(p, *request, Z);
      }
      /* create neighbor list */
      for (i = 0; i < *numnei; i++) {
	/* get atom index */
	ii = tb[tl[imd_nr] + i];
	/* get cell index */
	in_cell = kim.cell_list[cl_num[ii]];
	/* check if cell is in the test buffer */
	if (in_cell != -1) {
	  /* offset of the cell in the test buffer */
	  offset = kim.cell_offset[in_cell];
	  /* offset of the atom in the cell */
	  atom_offset = ii - cl_off[cl_num[ii]];
	  kim_nbl[i] = offset + atom_offset;
	  /* calculate relative postition vector if required */
	  if (kim.model_using_Rij) {
	    q = cell_array + cl_num[ii];
	    Rij[3 * i + 0] = ORT(q, atom_offset, X) - d1.x;
	    Rij[3 * i + 1] = ORT(q, atom_offset, Y) - d1.y;
	    Rij[3 * i + 2] = ORT(q, atom_offset, Z) - d1.z;
	  }
	}
	/* if this is not the case, then something went wrong */
	else {
	  kim_error("Cell index error");
	}
      }
      *nei1atom = kim_nbl;
      if (kim.model_using_Rij)
	*pRij = Rij;
      return KIM_STATUS_OK;
    } else if (*request >= p->n) {
      *atom = *request;
      *numnei = 0;
      return KIM_STATUS_OK;
    } else if (*request < 0) {
      return KIM_STATUS_NEIGH_INVALID_REQUEST;
    }
  } else
    /* there is only 0 or 1, so return INVALID MODE */
    return KIM_STATUS_NEIGH_INVALID_MODE;

  return 0;
}

/****************************************************************
 *
 *  void extend_test_buffer(int)
 *    allocate more memory for the data in the test buffer
 *    should not be needed very often
 *
 ****************************************************************/

void extend_test_buffer(int n)
{
  int   kimerror = KIM_STATUS_OK;
  test_buffer_t *buf;

  /* static arrays for data exchange with KIM */
  buf = (test_buffer_t *) KIM_API_get_test_buffer(kim.pkim, &kimerror);
  if (KIM_STATUS_OK > kimerror) {
    KIM_API_report_error(__LINE__, __FILE__, "Could not get the test buffer", kimerror);
    exit(EXIT_FAILURE);
  }

  buf->coords = realloc(buf->coords, 3 * n * sizeof(real));
  buf->forces = realloc(buf->forces, 3 * n * sizeof(real));
  buf->penergy = realloc(buf->penergy, n * sizeof(real));
  buf->ptypes = realloc(buf->ptypes, n * sizeof(int));
  kim.cell_index_atom = realloc(kim.cell_index_atom, n * sizeof(int));
  if (NULL == buf->coords || NULL == buf->forces || NULL == buf->penergy ||
	NULL == buf->ptypes || NULL == kim.cell_index_atom)
    kim_error("Could not increase test buffer size");
  buf->max_len = n;
}

/****************************************************************
 *
 *  void kim_warning(char *);
 *    print a fancy warning
 *
 ****************************************************************/

void kim_warning(char *warn)
{
  printf("\n*****************************************************************\n");
  printf("* WARNING: %s\n", warn);
  printf("*****************************************************************\n\n");
  fflush(stdout);
}

/****************************************************************
 *
 *  void kim_error(char *);
 *    print a fancy error message
 *
 ****************************************************************/

void kim_error(char *error)
{
  fprintf(stderr, "\n*****************************************************************\n");
  fprintf(stderr, "* ERROR: %s\n", error);
  fprintf(stderr, "*****************************************************************\n\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

/****************************************************************
 *
 *  int imd_process_dEdr(void **km, double* dEdr, double *r, double **dx,
 *  	int *i, int *j);
 *  our own process_dEdr function
 *
 ****************************************************************/

int imd_process_dEdr(void **km, double *grad, double *r, double **Rij, int *i, int *j)
{
  if (0.0 == *r)
    return KIM_STATUS_FAIL;
#if defined P_AXIAL || defined STRESS_TENS
  double fx = (*Rij)[0] * *grad / *r;
  double fy = (*Rij)[1] * *grad / *r;
  double fz = (*Rij)[2] * *grad / *r;

  /* determine on which cell we are working */
  cell *p = cell_array + kim.cell_ind;
  cell *q = cell_array + kim.cell_index_atom[*j];
#endif

#ifdef P_AXIAL
  vir_xx -= (*Rij)[0] * fx;
  vir_yy -= (*Rij)[1] * fy;
#ifndef TWOD
  vir_zz -= (*Rij)[2] * fz;
#endif
#else
  virial -= *r * *grad;
#endif

#ifdef STRESS_TENS
  if (do_press_calc) {
    int jj = *j - kim.cell_offset[kim.cell_list[kim.cell_index_atom[*j]]];
    fx *= 0.5;
    fy *= 0.5;
    fz *= 0.5;

    PRESSTENS(p, *i, xx) -= (*Rij)[0] * fx;
    PRESSTENS(p, *i, yy) -= (*Rij)[1] * fy;
    PRESSTENS(p, *i, zz) -= (*Rij)[2] * fz;
    PRESSTENS(p, *i, xy) -= (*Rij)[0] * fy;
    PRESSTENS(p, *i, yz) -= (*Rij)[1] * fz;
    PRESSTENS(p, *i, zx) -= (*Rij)[2] * fx;

    PRESSTENS(q, jj, xx) -= (*Rij)[0] * fx;
    PRESSTENS(q, jj, yy) -= (*Rij)[1] * fy;
    PRESSTENS(q, jj, zz) -= (*Rij)[2] * fz;
    PRESSTENS(q, jj, xy) -= (*Rij)[0] * fy;
    PRESSTENS(q, jj, yz) -= (*Rij)[1] * fz;
    PRESSTENS(q, jj, zx) -= (*Rij)[2] * fx;
  }
#endif

  return KIM_STATUS_OK;
}

#endif /* KIM */
