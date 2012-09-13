
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2011 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* Allocation and moving of atom data
*
*   Contains routines move_atom() and alloc_cell(), which were in
*   imd_geom.c and imd_geom_2d.c, but are dimension independent.  
*
* $Revision$
* $Date$
*
******************************************************************************/

#define INDEXED_ACCESS
#include "imd.h"

#ifdef VEC

/*****************************************************************************
*
*  Move atom from minicell to minicell. Atom stays in the same cell. 
*
******************************************************************************/

void move_atom_mini(minicell *to, minicell *from, int index)
{
  /* check the parameters */
  if ((0 > index) || (index >= from->n)) 
    error("move_atom: index argument out of range.");

  /* see if we need some space */
  if (to->n >= to->n_max) alloc_minicell(to,to->n_max+incrsz);

#ifdef MPI
  /* pointer from atom to index in its minicell */
  atoms.ind[from->ind[index]] = to->n;
#endif

  /* move the indices */
  to->ind[to->n++] = from->ind[index];
  if (index < from->n-1) from->ind[index] = from->ind[from->n-1];  
  from->n--;
}

/*****************************************************************************
*
*  Move atom from temporary cell to cell_array. 
*
******************************************************************************/

void insert_atom(minicell *to, cell *from, int index)
{
  /* see if we need some space */
  if (to->n >= to->n_max) alloc_minicell(to,to->n_max+incrsz);

  /* put index of the inserted atom */
  to->ind[to->n++] = atoms.n;

  /* put atom at the end of the cell; increase size if necessary */
  if (atoms.n >= atoms.n_max) alloc_cell(&atoms, 2*atoms.n_max);
  move_atom(&atoms, from, index);

#ifdef MPI
  /* pointer from atom to index in its minicell */
  atoms.ind[atoms.n-1] = to->n-1;
#endif
}

/******************************************************************************
*
*  Allocate a minicell
*
******************************************************************************/

void alloc_minicell(minicell *p, int count)
{
  if (0==count) {  /* cell is freed */
    free(p->ind);
  }
  else {  /* cell is allocated or enlarged */
    if (p->n_max==0) p->ind = NULL;
    p->ind = (int *) realloc( p->ind, count * sizeof(int) );
    if (NULL==p->ind) error("minicell allocation failed");
    p->n_max = count;
  }
}

#endif

/*****************************************************************************
*
*  Moves an atom from one cell to another. 
*  Does not move atoms between CPUs!
*
******************************************************************************/

void move_atom(cell *to, cell *from, int index)
{
  /* check the parameters */
  if ((0 > index) || (index >= from->n)) 
    error("move_atom: index argument out of range.");
 
  /* see if we need more space */
  if (to->n >= to->n_max) alloc_cell(to,to->n_max+incrsz);
  
  /* append atom to target cell */
  copy_atom_cell_cell(to, to->n, from, index);
  ++to->n;

  /* delete atom in original cell, by moving the last one to the empty slot */
  --from->n;
  if (index < from->n) copy_atom_cell_cell(from, index, from, from->n);

}

/*****************************************************************************
*
*  Low-level copying of atom data from one cell to another 
*
******************************************************************************/

void copy_atom_cell_cell(cell *to, int i, cell *from, int j)
{
    int k;
  to->ort X(i) = from->ort X(j); 
  to->ort Y(i) = from->ort Y(j); 
#ifndef TWOD
  to->ort Z(i) = from->ort Z(j); 
#endif

#ifdef BBOOST
  to->bb_refposone X(i) = from->bb_refposone X(j);
  to->bb_refposone Y(i) = from->bb_refposone Y(j); 
#ifndef TWOD
  to->bb_refposone Z(i) = from->bb_refposone Z(j); 
#endif
  to->bb_refpostwo X(i) = from->bb_refpostwo X(j);
  to->bb_refpostwo Y(i) = from->bb_refpostwo Y(j); 
#ifndef TWOD
  to->bb_refpostwo Z(i) = from->bb_refpostwo Z(j);
#endif
  to->bb_oldpos X(i) = from->bb_oldpos X(j);
  to->bb_oldpos Y(i) = from->bb_oldpos Y(j); 
#ifndef TWOD
  to->bb_oldpos Z(i) = from->bb_oldpos Z(j); 
#endif
  to->bb_neigh[i]->nbondsref1 = from->bb_neigh[j]->nbondsref1;
  to->bb_neigh[i]->nbondsref2 = from->bb_neigh[j]->nbondsref2;

  for (k=0; k< BB_MAXNEIGH; ++k) {
      to->bb_neigh[i]->distref1[k] =  from->bb_neigh[j]->distref1[k];
      to->bb_neigh[i]->distref2[k] =  from->bb_neigh[j]->distref2[k];
      to->bb_neigh[i]->numref1[k]  =  from->bb_neigh[j]->numref1[k];
      to->bb_neigh[i]->numref2[k]  =  from->bb_neigh[j]->numref2[k];      
    }
  
#endif /*BBOOST*/
  
#ifndef MONOLJ
  to->nummer [i] = from->nummer [j]; 
  to->sorte  [i] = from->sorte  [j]; 
  to->vsorte [i] = from->vsorte [j]; 
  to->masse  [i] = from->masse  [j]; 
#ifndef CBE
  to->pot_eng[i] = from->pot_eng[j];
#endif
#endif
#ifdef EAM2
  to->eam_rho[i] = from->eam_rho[j]; 
  to->eam_dF [i] = from->eam_dF [j]; 
#endif
#ifdef EEAM
  to->eeam_p_h[i] = from->eeam_p_h[j]; 
  to->eeam_dM [i] = from->eeam_dM [j]; 
#endif
#ifdef DAMP
  to->damp_f[i] = from->damp_f[j];
#endif
#ifdef ADP
  to->adp_mu   X(i)    = from->adp_mu   X(j); 
  to->adp_mu   Y(i)    = from->adp_mu   Y(j); 
  to->adp_mu   Z(i)    = from->adp_mu   Z(j); 
  to->adp_lambda[i].xx = from->adp_lambda[j].xx;   
  to->adp_lambda[i].yy = from->adp_lambda[j].yy;   
  to->adp_lambda[i].zz = from->adp_lambda[j].zz;   
  to->adp_lambda[i].yz = from->adp_lambda[j].yz;   
  to->adp_lambda[i].zx = from->adp_lambda[j].zx;   
  to->adp_lambda[i].xy = from->adp_lambda[j].xy;   
#endif
#ifdef VARCHG
  to->charge[i] = from->charge[j];
#endif
#ifdef DIPOLE
  to->dp_p_ind  X(i)   = from->dp_p_ind X(j);
  to->dp_p_ind  Y(i)   = from->dp_p_ind Y(j);
  to->dp_p_ind  Z(i)   = from->dp_p_ind Z(j);
  to->dp_E_old_1 X(i)  = from->dp_E_old_1 X(j);
  to->dp_E_old_1 Y(i)  = from->dp_E_old_1 Y(j);
  to->dp_E_old_1 Z(i)  = from->dp_E_old_1 Z(j);
  to->dp_E_old_2 X(i)  = from->dp_E_old_2 X(j);
  to->dp_E_old_2 Y(i)  = from->dp_E_old_2 Y(j);
  to->dp_E_old_2 Z(i)  = from->dp_E_old_2 Z(j);
  to->dp_E_old_3 X(i)  = from->dp_E_old_3 X(j);
  to->dp_E_old_3 Y(i)  = from->dp_E_old_3 Y(j);
  to->dp_E_old_3 Z(i)  = from->dp_E_old_3 Z(j);
#endif /* Dipole */
#ifdef CG
  to->h  X(i) = from->h X(j); 
  to->h  Y(i) = from->h Y(j); 
#ifndef TWOD
  to->h  Z(i) = from->h Z(j);
#endif 
  to->g  X(i) = from->g X(j); 
  to->g  Y(i) = from->g Y(j); 
#ifndef TWOD
  to->g  Z(i) = from->g Z(j);
#endif  
  to->old_ort  X(i) = from->old_ort X(j); 
  to->old_ort  Y(i) = from->old_ort Y(j); 
#ifndef TWOD
  to->old_ort  Z(i) = from->old_ort Z(j); 
#endif
#endif
#ifdef DISLOC
  to->Epot_ref  [i] = from->Epot_ref[j];
  to->ort_ref X (i) = from->ort_ref X(j);
  to->ort_ref Y (i) = from->ort_ref Y(j);
#ifndef TWOD
  to->ort_ref Z (i) = from->ort_ref Z(j);
#endif
#endif /* DISLOC */
#ifdef CNA
  to->mark[i] = from->mark[j];
#endif
#ifdef ADA
  to->adaType[i] = from->adaType[j];
  to->hopsToDefect[i] = from->hopsToDefect[j];
#endif
#ifdef AVPOS
  to->av_epot[i] = from->av_epot[j];
  to->avpos X(i) = from->avpos X(j);
  to->avpos Y(i) = from->avpos Y(j);
  to->sheet X(i) = from->sheet X(j);
  to->sheet Y(i) = from->sheet Y(j);
#ifndef TWOD
  to->avpos Z(i) = from->avpos Z(j);
  to->sheet Z(i) = from->sheet Z(j);
#endif
#endif /* AVPOS */
#ifdef NNBR
  to->nbanz[i] = from->nbanz[j]; 
#endif
#ifdef REFPOS
  to->refpos X(i) = from->refpos X(j);
  to->refpos Y(i) = from->refpos Y(j);
#ifndef TWOD
  to->refpos Z(i) = from->refpos Z(j);
#endif
#endif /* REFPOS */
#ifdef HC
  to->hcaveng[i]  = from->hcaveng[j];   
#endif
#ifdef STRESS_TENS
  to->presstens[i].xx = from->presstens[j].xx;   
  to->presstens[i].yy = from->presstens[j].yy;   
  to->presstens[i].xy = from->presstens[j].xy;   
#ifndef TWOD
  to->presstens[i].zz = from->presstens[j].zz;   
  to->presstens[i].yz = from->presstens[j].yz;   
  to->presstens[i].zx = from->presstens[j].zx;   
#endif
#ifdef AVPOS
  to->avpresstens[i].xx = from->avpresstens[j].xx;   
  to->avpresstens[i].yy = from->avpresstens[j].yy;   
  to->avpresstens[i].xy = from->avpresstens[j].xy;   
#ifndef TWOD
  to->avpresstens[i].zz = from->avpresstens[j].zz;   
  to->avpresstens[i].yz = from->avpresstens[j].yz;   
  to->avpresstens[i].zx = from->avpresstens[j].zx;  
#endif 
#endif /* AVPOS */
#endif /* STRESS_TENS */
#ifdef SHOCK
  to->pxavg[i] = from->pxavg[j];   
#endif
  to->impuls X(i) = from->impuls X(j); 
  to->impuls Y(i) = from->impuls Y(j); 
#ifndef TWOD
  to->impuls Z(i) = from->impuls Z(j); 
#endif
  to->kraft X(i) = from->kraft X(j); 
  to->kraft Y(i) = from->kraft Y(j); 
#ifndef TWOD
  to->kraft Z(i) = from->kraft Z(j); 
#endif
#ifdef CBE
  to->kraft W(i) = from->kraft W(j); 
#endif
#ifdef COVALENT
  /* neighbor table is not copied */
#endif
#ifdef NBLIST
  /* reference positions of neighbor list are not copied */
#endif
#ifdef UNIAX
  to->achse X(i) = from->achse X(j); 
  to->achse Y(i) = from->achse Y(j); 
  to->achse Z(i) = from->achse Z(j); 
  to->dreh_impuls X(i) = from->dreh_impuls X(j); 
  to->dreh_impuls Y(i) = from->dreh_impuls Y(j); 
  to->dreh_impuls Z(i) = from->dreh_impuls Z(j); 
  to->dreh_moment X(i) = from->dreh_moment X(j); 
  to->dreh_moment Y(i) = from->dreh_moment Y(j); 
  to->dreh_moment Z(i) = from->dreh_moment Z(j); 
#endif

}

#if defined(COVALENT) || defined(NNBR_TABLE)
/******************************************************************************
*
*  allocate neighbor table for one particle
*
******************************************************************************/

neightab *alloc_neightab(neightab *neigh, int count)
{
  if (0 == count) { /* deallocate */
    free(neigh->dist);
    free(neigh->typ);
    free(neigh->cl);
    free(neigh->num);
    free(neigh);
  } else { /* allocate */
    neigh = (neightab *) malloc(sizeof(neightab));
    if (neigh==NULL) {
      error("COVALENT: cannot allocate memory for neighbor table\n");
    }
    neigh->n     = 0;
    neigh->n_max = count;
    neigh->dist  = (real *)     malloc( count * SDIM * sizeof(real) );
    neigh->typ   = (shortint *) malloc( count * sizeof(shortint) );
    neigh->cl    = (void **)    malloc( count * sizeof(cellptr) );
    neigh->num   = (integer *)  malloc( count * sizeof(integer) );
    if ((neigh->dist==NULL) || (neigh->typ==0) ||
        (neigh->cl  ==NULL) || (neigh->num==0)) {
      error("COVALENT: cannot allocate memory for neighbor table");
    }
  }
  return(neigh);
}
#endif

#ifdef BBOOST
/******************************************************************************
*
*  allocate neighbor table for one particle
*
******************************************************************************/

bb_neightab *alloc_bb_neightab(bb_neightab *bb_neigh, int count)
{
  if (0 == count) { /* deallocate */
    free(bb_neigh->distref1);
    free(bb_neigh->distref2);
    free(bb_neigh->numref1);
    free(bb_neigh->numref2);
    free(bb_neigh); /* should I add this line ? */
  } else { /* allocate */
    bb_neigh = (bb_neightab *) malloc(sizeof(bb_neightab));
    if (bb_neigh==NULL) {
      error("BBOOST: cannot allocate memory for neighbor table\n");
    }
    bb_neigh->nbondsref1     = 0;
    bb_neigh->nbondsref2     = 0;
    /* careful: COVALENT uses dist to store in consecutive entries the vector between two atoms, we want only the distance */
    bb_neigh->distref1  = (real *)     malloc( count * sizeof(real) );
    bb_neigh->distref2  = (real *)     malloc( count * sizeof(real) );
    bb_neigh->numref1   = (integer *)  malloc( count * sizeof(integer) );
    bb_neigh->numref2   = (integer *)  malloc( count * sizeof(integer) );
    if ((bb_neigh->distref1==NULL) || (bb_neigh->distref2==NULL) ||
        (bb_neigh->numref1==NULL)  || (bb_neigh->numref2==NULL) ) {
      error("BBOOST: cannot allocate memory for neighbor table");
    }
  }
  return(bb_neigh);
}

#endif /*BBOOST*/

#if defined(COVALENT) || defined(NNBR_TABLE)
/******************************************************************************
*
*  increase already existing neighbor table for one particle
*
******************************************************************************/

void increase_neightab(neightab *neigh, int count)
{
  neigh->dist = (real *) realloc( neigh->dist, count * SDIM * sizeof(real));
  neigh->typ  = (shortint *) realloc( neigh->typ,  count * sizeof(shortint) );
  neigh->cl   = (void **)    realloc( neigh->cl,   count * sizeof(cellptr) );
  neigh->num  = (integer *)  realloc( neigh->num,  count * sizeof(integer) );
  if ((neigh->dist==NULL) || (neigh->typ==0) ||
      (neigh->cl  ==NULL) || (neigh->num==0)) {
    error("COVALENT: cannot extend memory for neighbor table");
  }
  neigh->n_max = count;
  /* update maximal neighbor table length on *this* CPU */
  neigh_len = MAX( neigh_len, count );
}
#endif



/******************************************************************************
*
*  Allocate memory with prescribed alignment 
*
*  This routine can replace malloc, calloc, realloc, free.
*
*  The parameters are:
*
*    p:      on entry, address of pointer to old memory block
*            on exit, address of pointer to new memory block
*    count:  size of new memory (number of items)
*            if zero, old memory is just deallocated
*    size:   size of one item (in bytes)
*    align:  memory is aligned at <align> byte boundaries
*            must be a multiple of the pointer size
*    ncopy:  number of items to copy from old to new memory
*            if non-zero, array at old location is deallocated afterwards
*            if negative, nothing is copied, but old memory is deallocated
*            if zero, there is no old memory to deallocate or copy
*    clear:  if non-zero, new memory is zeroed
*    name:   name of memory block (for error messages)
*
******************************************************************************/

void memalloc(void *p, int count, int size, int align, int ncopy, int clear,
              char *name)
{
  void *new, **old = (void **)p;
  int  ret, len, a = align - 1;

// #ifdef debugLo
//     printf("print "" count = %d, size = %d, align = %d, ncopy = %d, clear = %d, name = %s "" checking by Lo! \n", count,size,align,ncopy,clear,name);fflush(stdout);
// #endif
  if (count>0) {  /* allocate count * size bytes */
    len = (count * size + a) & (~a);  /* enlarge to multiple of align */
#ifdef MEMALIGN
    ret = posix_memalign(&new, align, len);
    if (ret==EINVAL) { /* align must be a multiple of the pointer size */
      error("invalid alignment request in memory allocation");
    }
    else if (ret==ENOMEM) { /* out of memory */
      error_str("Cannot allocate memory for %s", name);
    }
#else

    new = malloc( len );

    if (NULL==new) { /* out of memory */
      error_str("Cannot allocate memory for %s", name);
    }
#endif
    else {  /* allocation succeded */
      if (clear  ) memset(new, 0, len);              /* zero new memory */

#ifdef debugLo
  //   printf("    ************************* \n");fflush(stdout);
//     printf("********************************* \n");fflush(stdout);
//     printf("print "" check memcpy start "" checking by Lo! \n");fflush(stdout);
//     printf("print "" ncopy = %d "" checking by Lo! \n", ncopy);fflush(stdout); 
#endif
      if (ncopy>0) memcpy(new, *old, ncopy * size);  /* copy old data */
#ifdef debugLo
 //   printf("print "" check memcpy end "" checking by Lo! \n");fflush(stdout);
//     printf("********************************* \n");fflush(stdout);
//     printf("    ************************* \n");fflush(stdout); 
#endif
      if (ncopy  ) free(*old); 
                      /* deallocate old data */
      *old = new;
    }
  }
  else {  /* deallocate */
    if (ncopy) free(*old);
    *old = NULL;
  }

}

/******************************************************************************
*
*  Allocate memory for a cell
*
******************************************************************************/

void alloc_cell(cell *p, int n)
{
  int i, ncopy;
#ifdef CBE
  int al=128;
#else
  int al=8;
#endif

  /* cells are either deallocated or increased; they never shrink */

  if ((n>0) && (n < p->n_max)) error("cells cannot shrink");

  /* cell is to be deallocated or has just been initialized -> no valid data */
  if ((0==n) || (0==p->n_max)) {
    p->n = 0;
#ifdef VEC
    p->n_buf = 0;
#endif
  }

#if (defined(COVALENT) || defined(NNBR_TABLE)) && !defined(TWOD)
  /* if cell is to be deallocated, begin with neighbor tables */
  if (0==n) {
    for (i=0; i<p->n_max; ++i) {
      p->neigh[i] = alloc_neightab(p->neigh[i], 0);
    }
  }
#endif

#if (defined(BBOOST) && !defined(TWOD))
  /* if cell is to be deallocated, begin with neighbor tables */
  if (0==n) {
    for (i=0; i< BB_MAXNEIGH; ++i) {
      p->bb_neigh[i] = alloc_bb_neightab(p->bb_neigh[i], 0);
    }
  }
#endif /*BBOOST*/

  /* if there are valid particles in cell, copy them to new storage */
  if (0==p->n_max) {
    ncopy = 0;
  }
  else {
#ifdef VEC
      ncopy = ( (0==n) || (0==MAX(p->n,p->n_buf) )) ? -1 : MAX(p->n,p->n_buf);
#else
      ncopy = ( (0==n) || (0==p->n)              ) ? -1 : p->n;
#endif
  }

  /* allocate memory */
  memalloc( &p->ort,      n*SDIM, sizeof(real), al, ncopy*SDIM, 0, "ort" );

  
#ifdef BBOOST
  memalloc( &p->bb_refposone,      n*SDIM, sizeof(real), al, ncopy*SDIM, 0, "bb_refposone" );
  memalloc( &p->bb_refpostwo,      n*SDIM, sizeof(real), al, ncopy*SDIM, 0, "bb_refpostwo" );
  memalloc( &p->bb_oldpos,         n*SDIM, sizeof(real), al, ncopy*SDIM, 0, "bb_oldpos" );
#endif /*BBOOST*/
  
  memalloc( &p->impuls,   n*SDIM, sizeof(real), al, ncopy*SDIM, 0, "impuls" );
  memalloc( &p->kraft,    n*SDIM, sizeof(real), al, ncopy*SDIM, 0, "kraft" );
#ifndef MONOLJ
  memalloc( &p->nummer,   n, sizeof(integer),  al, ncopy, 0, "nummer" );
  memalloc( &p->sorte,    n, sizeof(shortint), al, ncopy, 0, "sorte" );
  memalloc( &p->vsorte,   n, sizeof(shortint), al, ncopy, 0, "vsorte" );
  memalloc( &p->masse,    n, sizeof(real),     al, ncopy, 0, "masse" );
#ifndef CBE
  memalloc( &p->pot_eng,  n, sizeof(real),     al, ncopy, 0, "pot_eng" );
#endif
#endif
#ifdef EAM2
  memalloc( &p->eam_rho,  n, sizeof(real),     al, ncopy, 0, "eam_rho" );
  memalloc( &p->eam_dF,   n, sizeof(real),     al, ncopy, 0, "eam_dF" );
#endif
#ifdef EEAM
  memalloc( &p->eeam_p_h, n, sizeof(real),     al, ncopy, 0, "eeam_p_h" );
  memalloc( &p->eeam_dM,  n, sizeof(real),     al, ncopy, 0, "eeam_dM" );
#endif
#ifdef DAMP
  memalloc( &p->damp_f,   n, sizeof(real),     al, ncopy, 0, "damp_f" );
#endif
#ifdef ADP
  memalloc( &p->adp_mu,   n*SDIM, sizeof(real),   al, ncopy*SDIM, 0, "adp_mu");
  memalloc( &p->adp_lambda, n, sizeof(sym_tensor), al, ncopy, 0, "adp_lambda");
#endif
#ifdef VARCHG
  memalloc( &p->charge,    n,     sizeof(real), al, ncopy, 0, "charge" );
#endif
#ifdef SM
  memalloc(&p->chi_sm,   n, sizeof(real),   al, ncopy, 0, "chi_sm");
  memalloc(&p->z_sm,   n, sizeof(real),   al, ncopy, 0, "z_sm");
  memalloc(&p->j_sm,   n, sizeof(real),   al, ncopy, 0, "j_sm");
  memalloc(&p->v_sm,   n, sizeof(real),   al, ncopy, 0, "v_sm");

  memalloc(&p->b_sm,   n, sizeof(real),   al, ncopy, 0, "b_sm");
  memalloc(&p->x_sm,   n, sizeof(real),   al, ncopy, 0, "x_sm");
  memalloc(&p->r_sm,   n, sizeof(real),   al, ncopy, 0, "r_sm");
  memalloc(&p->d_sm,   n, sizeof(real),   al, ncopy, 0, "d_sm");
  memalloc(&p->s_sm,   n, sizeof(real),   al, ncopy, 0, "s_sm");
  memalloc(&p->q_sm,   n, sizeof(real),   al, ncopy, 0, "q_sm");
#endif
#ifdef DIPOLE
  memalloc( &p->dp_E_stat, n*DIM, sizeof(real), al, ncopy*DIM, 0, "dp_E_stat");
  memalloc( &p->dp_E_ind,  n*DIM, sizeof(real), al, ncopy*DIM, 0, "dp_E_ind");
  memalloc( &p->dp_E_old_1,n*DIM, sizeof(real), al, ncopy*DIM, 1, "dp_E_old_1");
  memalloc( &p->dp_E_old_2,n*DIM, sizeof(real), al, ncopy*DIM, 1, "dp_E_old_2");
  memalloc( &p->dp_E_old_3,n*DIM, sizeof(real), al, ncopy*DIM, 1, "dp_E_old_3");
  memalloc( &p->dp_p_stat, n*DIM, sizeof(real), al, ncopy*DIM, 0, "dp_p_stat");
  memalloc( &p->dp_p_ind , n*DIM, sizeof(real), al, ncopy*DIM, 1, "dp_p_ind");
#endif
#ifdef CG
  memalloc( &p->h,        n*SDIM, sizeof(real), al, ncopy*SDIM, 0, "h" );
  memalloc( &p->g,        n*SDIM, sizeof(real), al, ncopy*SDIM, 0, "g" );
  memalloc( &p->old_ort,  n*SDIM, sizeof(real), al, ncopy*SDIM, 0, "old_ort" );
#endif
#ifdef NNBR
  memalloc( &p->nbanz,    n, sizeof(shortint), al, ncopy, 0, "nbanz" );
#endif
#ifdef DISLOC
  memalloc( &p->Epot_ref, n,      sizeof(real), al, ncopy,      1,"Epot_ref" );
  memalloc( &p->ort_ref,  n*SDIM, sizeof(real), al, ncopy*SDIM, 1, "ort_ref" );
#endif
#ifdef CNA
  memalloc( &p->mark,     n, sizeof(long), al, ncopy, 1, "mark" );
#endif
#ifdef ADA
  memalloc( &p->adaType,  n, sizeof(shortint), al, ncopy, 127, "adaType" );
  memalloc( &p->hopsToDefect,  n, sizeof(shortint), al, ncopy, 127, "hopsToDefect" );
#endif
#ifdef NYETENSOR
  memalloc( &p->nyeTens, n, sizeof(nyeTensorInfo*), al, ncopy, 1, "nyeTensorInfo" );
#endif
#ifdef AVPOS
  memalloc( &p->av_epot,  n,      sizeof(real), al, ncopy,      1, "av_epot" );
  memalloc( &p->avpos,    n*SDIM, sizeof(real), al, ncopy*SDIM, 1, "avpos"  );
  memalloc( &p->sheet,    n*SDIM, sizeof(real), al, ncopy*SDIM, 0, "sheet"  );
#endif
#ifdef REFPOS
  memalloc( &p->refpos,   n*SDIM, sizeof(real), al, ncopy*SDIM, 1, "refpos" );
#endif
#ifdef HC
  memalloc( &p->hcaveng,  n,      sizeof(real), al, ncopy,      0, "hcaveng");
#endif
#ifdef STRESS_TENS
  memalloc( &p->presstens, n, sizeof(sym_tensor), al, ncopy, 0, "presstens" );
#ifdef AVPOS
  memalloc( &p->avpresstens, n, sizeof(sym_tensor), al, ncopy, 0, "avpresstens" );
#endif
#endif
#ifdef SHOCK
  memalloc( &p->pxavg, n, sizeof(real), al, ncopy, 1, "pxavg" );
#endif
#if defined(COVALENT) || defined(NNBR_TABLE)
  memalloc( &p->neigh, n, sizeof(neighptr), al, p->n_max, 0, "neigh" );
  for (i=p->n_max; i<n; ++i) {
    p->neigh[i] = alloc_neightab(p->neigh[i], neigh_len);
  }
#endif
#ifdef BBOOST
// #ifdef debugLo
//     printf("    ************************* \n");fflush(stdout);
//     printf("********************************* \n");fflush(stdout);
//     printf("print "" check alloc_cell start "" checking by Lo! \n");fflush(stdout);
// #endif  

    memalloc( &p->bb_neigh, n, sizeof(bb_neighptr), al, p->n_max, 0, "bb_neigh" );

// #ifdef debugLo
//     printf("print "" check alloc_cell end "" checking by Lo! \n");fflush(stdout);
//     printf("********************************* \n");fflush(stdout);
//     printf("    ************************* \n");fflush(stdout);
// #endif
    // for (i=BB_MAXNEIGH; i<n; ++i) {
     for (i=p->n_max; i<n; ++i) {
    p->bb_neigh[i] = alloc_bb_neightab(p->bb_neigh[i], BB_MAXNEIGH);
  }
#endif /* BBOOST */
  
#ifdef NBLIST
  memalloc( &p->nbl_pos,  n*SDIM, sizeof(real), al, ncopy*SDIM, 0, "nbl_pos" );
#endif
#ifdef UNIAX
  memalloc( &p->achse,       n*SDIM, sizeof(real), al, ncopy*SDIM, 0, "achse");
  memalloc( &p->dreh_impuls, n*SDIM, sizeof(real), al, ncopy*SDIM, 0, "dreh_impuls");
  memalloc( &p->dreh_moment, n*SDIM, sizeof(real), al, ncopy*SDIM, 0, "dreh_moment");
#endif
#if defined(VEC) && defined(MPI)
  memalloc( &p->ind, n, sizeof(integer), al, ncopy, 0, "ind" );
#endif
 
  p->n_max = n;

}
