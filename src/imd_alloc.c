/******************************************************************************
*
* Allocation and moving of atom data
*
*   Contains routines move_atom() and alloc_cell(), which were in
*   imd_geom.c and imd_geom_2d.c, but are dimension independent.  
*
* $RCSfile$
* $Revision$
* $Date$
*
******************************************************************************/

#include "imd.h"


/******************************************************************************
*
*  Moves an atom from one cell to another. 
*  Does not move atoms between CPUs!
*  Does not really belong to file imd_alloc.c either!
*
******************************************************************************/

void move_atom(ivektor cellc, cell *from, int index)

{
  cell *to;

#ifdef TWOD
  to = PTR_2D_VV(cell_array,cellc,cell_dim);
#else
  to = PTR_3D_VV(cell_array,cellc,cell_dim);
#endif

  /* Check the parameters */
  if ((0 > index) || (index >= from->n)) 
    error("move_atom: index argument out of range.");
  
  /* See if we need some space */
  if (to->n >= to->n_max) alloc_cell(to,to->n_max+incrsz);

  /* Got some space, move atom */
  to->ort X(to->n) = from->ort X(index); 
  to->ort Y(to->n) = from->ort Y(index); 
#ifndef TWOD
  to->ort Z(to->n) = from->ort Z(index); 
#endif
  
  to->kraft X(to->n) = from->kraft X(index); 
  to->kraft Y(to->n) = from->kraft Y(index); 
#ifndef TWOD
  to->kraft Z(to->n) = from->kraft Z(index); 
#endif
  
  to->impuls X(to->n) = from->impuls X(index); 
  to->impuls Y(to->n) = from->impuls Y(index); 
#ifndef TWOD
  to->impuls Z(to->n) = from->impuls Z(index); 
#endif
#ifndef MONOLJ    
  to->masse[to->n] = from->masse[index]; 
  to->sorte[to->n] = from->sorte[index]; 
  to->nummer[to->n] = from->nummer[index]; 
  to->pot_eng[to->n] = from->pot_eng[index];
#ifdef DISLOC
  to->Epot_ref[to->n] = from->Epot_ref[index];
  to->ort_ref X (to->n) = from->ort_ref X(index);
  to->ort_ref Y (to->n) = from->ort_ref Y(index);
#ifndef TWOD
  to->ort_ref Z (to->n) = from->ort_ref Z(index);
#endif
#endif
#ifdef REFPOS
  to->refpos X(to->n) = from->refpos X(index);
  to->refpos Y(to->n) = from->refpos Y(index);
#ifndef TWOD
  to->refpos Z(to->n) = from->refpos Z(index);
#endif
#endif /* REFPOS */
#endif /* not MONOLJ */
  ++to->n;

  /* Delete atom in original cell */

  --from->n;

  if (0 < from->n) {

    from->ort X(index) = from->ort X(from->n); 
    from->ort Y(index) = from->ort Y(from->n); 
#ifndef TWOD
    from->ort Z(index) = from->ort Z(from->n); 
#endif

    from->kraft X(index) = from->kraft X(from->n); 
    from->kraft Y(index) = from->kraft Y(from->n); 
#ifndef TWOD
    from->kraft Z(index) = from->kraft Z(from->n); 
#endif

    from->impuls X(index) = from->impuls X(from->n); 
    from->impuls Y(index) = from->impuls Y(from->n); 
#ifndef TWOD
    from->impuls Z(index) = from->impuls Z(from->n); 
#endif
#ifndef MONOLJ
    from->masse[index] = from->masse[from->n]; 
    from->sorte[index] = from->sorte[from->n]; 
    from->nummer[index] = from->nummer[from->n]; 
    from->pot_eng[index] = from->pot_eng[from->n];
#ifdef DISLOC
    from->Epot_ref[index]  = from->Epot_ref[from->n];
    from->ort_ref X(index) = from->ort_ref X(from->n);
    from->ort_ref Y(index) = from->ort_ref Y(from->n);
#ifndef TWOD
    from->ort_ref Z(index) = from->ort_ref Z(from->n);
#endif
#endif
#ifdef REFPOS
    from->refpos X(index) = from->refpos X(from->n);
    from->refpos Y(index) = from->refpos Y(from->n);
#ifndef TWOD
    from->refpos Z(index) = from->refpos Z(from->n);
#endif
#endif /* REFPOS */
#endif /* not MONOLJ */
  }
}


#ifdef DISLOC
void reset_Epot_ref(void)
{
  int  k;
#pragma omp parallel for
  for (k=0; k<ncells; ++k) {
    int  i;
    cell *p;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) {
      p->Epot_ref[i] = p->pot_eng[i];
    }
  }
}
#endif


#ifdef TTBP
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
    free(neigh);
  } else { /* allocate */
    neigh = (neightab *) malloc(sizeof(neightab));
    if (neigh==NULL) {
      error("TTBP: cannot allocate memory for neighbor table\n");
    }
    neigh->n     = 0;
    neigh->n_max = count;
    neigh->dist  = (real *)     malloc( count * DIM * sizeof(real) );
    neigh->typ   = (shortint *) malloc( count *  sizeof(shortint) );
    if ((neigh->dist==NULL) || (neigh->typ==0)) {
      error("TTBP: cannot allocate memory for neighbor table");
    }
  }
  return(neigh);
}
#endif


/******************************************************************************
*
*   Allocates memory for a cell. Space is allocated in a single
*  chunk, so that we can send a cell in one MPI operation.
*
******************************************************************************/

void alloc_cell(cell *thecell, int count)

{
  void *space;
  cell newcell;
  int i;
  real *tmp;
  int newcellsize;

/* Cell is freed */
  if (0==count) {
    thecell->n = 0;
    
    newcell.ort =    NULL;
    newcell.impuls = NULL;
    newcell.kraft  = NULL;
#ifndef MONOLJ
    newcell.nummer = NULL;
    newcell.pot_eng= NULL;
#ifdef DISLOC
    newcell.Epot_ref = NULL;
    newcell.ort_ref = NULL;
#endif
#ifdef REFPOS
    newcell.refpos = NULL;
#endif
#ifdef TRANSPORT
    newcell.heatcond = NULL;
#endif
#ifdef STRESS_TENS
    newcell.presstens = NULL;
    newcell.presstens_offdia = NULL;
#endif
    newcell.sorte  = NULL;
    newcell.masse  = NULL;
#endif
#ifdef TTBP
    newcell.neigh  = NULL;
#endif
  }
  else {
        /* Cell IS allocated */
    
    /* Calulate newcell size */
    newcellsize = count *
        ( DIM * sizeof(real) + /* ort */
          DIM * sizeof(real) + /* impuls */
          DIM * sizeof(real) ); /* kraft */
  
        /* Get some space */
    space = malloc(newcellsize);


        /* Calculate Pointers */
  
    tmp = (real *) space;
    newcell.ort    = tmp; tmp += DIM * count;
    newcell.impuls = tmp; tmp += DIM * count;
    newcell.kraft  = tmp; tmp += DIM * count;
#ifndef MONOLJ
        /* Allocate rest of variables */
    newcell.nummer = (integer * ) malloc(count * sizeof(integer) );
    newcell.pot_eng= (real    * ) malloc(count * sizeof(real)    );
#ifdef DISLOC
    newcell.Epot_ref = (real *) calloc(count,sizeof(real));
    newcell.ort_ref = (real *) calloc(count*DIM, sizeof(real));
#endif
#ifdef REFPOS
    newcell.refpos = (real *) malloc(count*DIM*sizeof(real));
#endif
#ifdef TRANSPORT
    newcell.heatcond = (real *) malloc(count*sizeof(real));
#endif
#ifdef STRESS_TENS
    newcell.presstens = (real *) malloc(count*DIM*sizeof(real));
    newcell.presstens_offdia = (real *) malloc(count*DIM*sizeof(real));
#endif
#if (defined(TTBP) && !defined(TWOD))
    newcell.neigh = (neightab **) malloc( count * sizeof(neighptr) );
    if (NULL == newcell.neigh) {
      error("TTBP: cannot allocate neighbor tables");
    }
    for (i=0; i<thecell->n_max; ++i) {
      newcell.neigh[i] = thecell->neigh[i];
    }
    for (i=thecell->n_max; i<count; ++i) {
      newcell.neigh[i] = alloc_neightab(newcell.neigh[i], ttbp_len);
    }
#endif
    newcell.sorte  = (shortint* ) malloc(count * sizeof(shortint));
    newcell.masse  = (real    * ) malloc(count * sizeof(real)    );
#endif

    if ((NULL == space)
#ifndef MONOLJ
        || (NULL == newcell.nummer)
        || (NULL == newcell.pot_eng)
#ifdef DISLOC
        || (NULL == newcell.Epot_ref)
        || (NULL == newcell.ort_ref)
#endif
#ifdef REFPOS
        || (NULL == newcell.refpos)
#endif
#ifdef TRANSPORT
        || (NULL == newcell.heatcond)
#endif
#ifdef STRESS_TENS
        || (NULL == newcell.presstens)
        || (NULL == newcell.presstens_offdia)
#endif
        || (NULL == newcell.sorte)
        || (NULL == newcell.masse)
#endif
        ) {
      printf("Want %d bytes.\n",newcellsize);
#ifdef CRAY
      malloc_stats(0);
#endif
      printf("Have %d atoms.\n",natoms);
      error("Cannot allocate memory for cell.");
    }
  }
  
  if (0 == thecell->n_max) { /* cell is just initialized */
    thecell->n = 0;
  } else {

    if (count < thecell->n_max) { /* cell shrinks, data is invalidated */
      thecell->n = 0;
#if (defined(TTBP) && !defined(TWOD))
      /* deallocate all neighbor tables */
      for (i=0; i<thecell->n_max; ++i) {
        thecell->neigh[i] = alloc_neightab(thecell->neigh[i],0);
      }
      free(thecell->neigh);
#endif
    }

    /* if there are valid particles in cell, copy them to new cell */
    if (thecell->n > 0) {
      /* cell is enlarged, copy data from old to newcell location */
      memcpy(newcell.ort   , thecell->ort,    thecell->n * DIM * sizeof(real));
      memcpy(newcell.impuls, thecell->impuls, thecell->n * DIM * sizeof(real));
      memcpy(newcell.kraft , thecell->kraft,  thecell->n * DIM * sizeof(real));
#ifndef MONOLJ
      memcpy(newcell.nummer,  thecell->nummer,  thecell->n * sizeof(integer));
      memcpy(newcell.pot_eng, thecell->pot_eng, thecell->n * sizeof(real));
#ifdef DISLOC
      memcpy(newcell.Epot_ref, thecell->Epot_ref, 
                               thecell->n * sizeof(real));
      memcpy(newcell.ort_ref,  thecell->ort_ref, 
                               thecell->n * DIM * sizeof(real));
#endif
#ifdef REFPOS
      memcpy(newcell.refpos, thecell->refpos, thecell->n * DIM * sizeof(real));
#endif
#ifdef TRANSPORT
      memcpy(newcell.heatcond,  thecell->heatcond,  thecell->n * sizeof(real));
#endif
#ifdef STRESS_TENS
      memcpy(newcell.presstens, thecell->presstens, 
             thecell->n * DIM * sizeof(real));
      memcpy(newcell.presstens_offdia, thecell->presstens_offdia, 
             thecell->n * DIM * sizeof(real));
#endif
      memcpy(newcell.sorte ,  thecell->sorte,   thecell->n * sizeof(shortint));
      memcpy(newcell.masse ,  thecell->masse,   thecell->n * sizeof(real));
#endif
    }

    /* deallocate old cell */
    free(thecell->ort);
#ifndef MONOLJ
    free(thecell->nummer);
    free(thecell->pot_eng);
#ifdef DISLOC
    free(thecell->Epot_ref);
    free(thecell->ort_ref);
#endif
#ifdef REFPOS
    free(thecell->refpos);
#endif
#ifdef TRANSPORT
    free(thecell->heatcond);
#endif
#ifdef STRESS_TENS
    free(thecell->presstens);
    free(thecell->presstens_offdia);
#endif
    free(thecell->sorte);
    free(thecell->masse);
#endif
#ifdef TTBP
    free(thecell->neigh);
#endif
  }

  /* set pointers to contents of new cell */
  thecell->ort      = newcell.ort;
  thecell->impuls   = newcell.impuls;
  thecell->kraft    = newcell.kraft;
#ifndef MONOLJ
  thecell->nummer   = newcell.nummer;
  thecell->pot_eng  = newcell.pot_eng;
#ifdef DISLOC
  thecell->Epot_ref = newcell.Epot_ref;
  thecell->ort_ref  = newcell.ort_ref;
#endif
#ifdef REFPOS
  thecell->refpos   = newcell.refpos;
#endif
#ifdef TRANSPORT
  thecell->heatcond = newcell.heatcond;
#endif
#ifdef STRESS_TENS
  thecell->presstens        = newcell.presstens;
  thecell->presstens_offdia = newcell.presstens_offdia;
#endif
#ifdef TTBP
  thecell->neigh    = newcell.neigh;
#endif
  thecell->sorte    = newcell.sorte;
  thecell->masse    = newcell.masse;
#endif

  thecell->n_max    = count;

}

