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

/*****************************************************************************
*
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
#ifdef EAM2
  to->eam2_rho_h[to->n] = from->eam2_rho_h[index]; 
#endif
#ifdef ORDPAR
#ifndef TWOD
  to->nbanz[to->n] = from->nbanz[index]; 
#endif
#endif
#ifdef DISLOC
  to->Epot_ref[to->n] = from->Epot_ref[index];
  to->ort_ref X (to->n) = from->ort_ref X(index);
  to->ort_ref Y (to->n) = from->ort_ref Y(index);
#ifndef TWOD
  to->ort_ref Z (to->n) = from->ort_ref Z(index);
#endif
#endif /* DISLOC */
#ifdef REFPOS
  to->refpos X(to->n) = from->refpos X(index);
  to->refpos Y(to->n) = from->refpos Y(index);
#ifndef TWOD
  to->refpos Z(to->n) = from->refpos Z(index);
#endif
#endif /* REFPOS */
#ifdef UNIAX
  to->traeg_moment[to->n] = from->traeg_moment[index]; 
  to->achse X(to->n) = from->achse X(index); 
  to->achse Y(to->n) = from->achse Y(index); 
  to->achse Z(to->n) = from->achse Z(index); 
  to->shape X(to->n) = from->shape X(index); 
  to->shape Y(to->n) = from->shape Y(index); 
  to->shape Z(to->n) = from->shape Z(index); 
  to->pot_well X(to->n) = from->pot_well X(index); 
  to->pot_well Y(to->n) = from->pot_well Y(index); 
  to->pot_well Z(to->n) = from->pot_well Z(index); 
  to->dreh_moment X(to->n) = from->dreh_moment X(index); 
  to->dreh_moment Y(to->n) = from->dreh_moment Y(index); 
  to->dreh_moment Z(to->n) = from->dreh_moment Z(index); 
  to->dreh_impuls X(to->n) = from->dreh_impuls X(index); 
  to->dreh_impuls Y(to->n) = from->dreh_impuls Y(index); 
  to->dreh_impuls Z(to->n) = from->dreh_impuls Z(index); 
#endif
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
#ifdef EAM2
    from->eam2_rho_h[index] = from->eam2_rho_h[from->n];
#endif
#ifdef ORDPAR
#ifndef TWOD
    from->masse[index] = from->masse[from->n];
#endif
#endif
#ifdef DISLOC
    from->Epot_ref[index]  = from->Epot_ref[from->n];
    from->ort_ref X(index) = from->ort_ref X(from->n);
    from->ort_ref Y(index) = from->ort_ref Y(from->n);
#ifndef TWOD
    from->ort_ref Z(index) = from->ort_ref Z(from->n);
#endif
#endif /* DISLOC */
#ifdef REFPOS
    from->refpos X(index) = from->refpos X(from->n);
    from->refpos Y(index) = from->refpos Y(from->n);
#ifndef TWOD
    from->refpos Z(index) = from->refpos Z(from->n);
#endif
#endif /* REFPOS */
#ifdef UNIAX
    from->traeg_moment[index] = from->traeg_moment[from->n]; 
    from->achse X(index) = from->achse X(from->n); 
    from->achse Y(index) = from->achse Y(from->n); 
    from->achse Z(index) = from->achse Z(from->n); 
    from->shape X(index) = from->shape X(from->n); 
    from->shape Y(index) = from->shape Y(from->n); 
    from->shape Z(index) = from->shape Z(from->n); 
    from->pot_well X(index) = from->pot_well X(from->n); 
    from->pot_well Y(index) = from->pot_well Y(from->n); 
    from->pot_well Z(index) = from->pot_well Z(from->n); 
    from->dreh_moment X(index) = from->dreh_moment X(from->n); 
    from->dreh_moment Y(index) = from->dreh_moment Y(from->n); 
    from->dreh_moment Z(index) = from->dreh_moment Z(from->n); 
    from->dreh_impuls X(index) = from->dreh_impuls X(from->n); 
    from->dreh_impuls Y(index) = from->dreh_impuls Y(from->n); 
    from->dreh_impuls Z(index) = from->dreh_impuls Z(from->n); 
#endif
#endif /* not MONOLJ */
  }
}


#ifdef DISLOC
void reset_Epot_ref(void)
{
  int  k;
#ifdef _OPENMP
#pragma omp parallel for
#endif
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


#ifdef COVALENT
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
    neigh->dist  = (real *)     malloc( count * DIM * sizeof(real) );
    neigh->typ   = (shortint *) malloc( count * sizeof(shortint) );
    neigh->cl    = (void *)     malloc( count * sizeof(cellptr) );
    neigh->num   = (integer *)  malloc( count * sizeof(integer) );
    if ((neigh->dist==NULL) || (neigh->typ==0) ||
        (neigh->cl  ==NULL) || (neigh->num==0)) {
      error("COVALENT: cannot allocate memory for neighbor table");
    }
  }
  return(neigh);
}
#endif


/******************************************************************************
*
*  Allocates memory for a cell. Space is allocated in a single
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
    
    newcell.ort = NULL;
    newcell.impuls = NULL;
    newcell.kraft = NULL;
#ifndef MONOLJ
    newcell.nummer = NULL;
    newcell.sorte  = NULL;
    newcell.masse  = NULL;
    newcell.pot_eng= NULL;
#ifdef EAM2
    newcell.eam2_rho_h =NULL;
#endif
#ifdef ORDPAR
#ifndef TWOD
    newcell.nbanz = NULL;
#endif
#endif
#ifdef DISLOC
    newcell.Epot_ref = NULL;
    newcell.ort_ref = NULL;
#endif
#ifdef REFPOS
    newcell.refpos = NULL;
#endif
#ifdef NVX
    newcell.heatcond = NULL;
#endif
#ifdef STRESS_TENS
    newcell.presstens = NULL;
    newcell.presstens_offdia = NULL;
#endif
#ifdef COVALENT
    newcell.neigh  = NULL;
#endif
#ifdef UNIAX
    newcell.traeg_moment = NULL;
    newcell.achse = NULL;
    newcell.shape = NULL;
    newcell.pot_well = NULL;
    newcell.dreh_impuls = NULL;
    newcell.dreh_moment = NULL;
#endif
#endif /* not MONOLJ */
  }
  else {
        /* Cell IS allocated */
    
    /* Calulate newcell size */
#ifdef UNIAX
    newcellsize = count *
        ( DIM * sizeof(real) + /* ort */
          DIM * sizeof(real) + /* achse */
          DIM * sizeof(real) + /* impuls */
          DIM * sizeof(real) + /* drehimpuls */
          DIM * sizeof(real) + /* kraft */
          DIM * sizeof(real) ); /* drehmoment */
#else
    newcellsize = count *
        ( DIM * sizeof(real) + /* ort */
          DIM * sizeof(real) + /* impuls */
          DIM * sizeof(real) ); /* kraft */
#endif  
        /* Get some space */
    space = malloc(newcellsize);


        /* Calculate Pointers */
  
    tmp = (real *) space;
    newcell.ort = tmp; tmp += DIM * count;
    newcell.impuls = tmp; tmp += DIM * count;
    newcell.kraft  = tmp; tmp += DIM * count;
#ifdef UNIAX
    newcell.achse = tmp; tmp += DIM * count;
    newcell.dreh_impuls = tmp; tmp += DIM * count;
    newcell.dreh_moment = tmp; tmp += DIM * count;
#endif
#ifndef MONOLJ
    /* Allocate rest of variables */
    newcell.nummer = (integer * ) malloc(count * sizeof(integer) );
    newcell.sorte  = (shortint* ) malloc(count * sizeof(shortint));
    newcell.masse  = (real    * ) malloc(count * sizeof(real)    );
    newcell.pot_eng= (real    * ) malloc(count * sizeof(real)    );
#ifdef EAM2
    newcell.eam2_rho_h = (real *) malloc(count * sizeof(real)    );
#endif
#ifdef ORDPAR
#ifndef TWOD
    newcell.nbanz = (shortint *) malloc(count * sizeof(shortint));
#endif
#endif
#ifdef DISLOC
    newcell.Epot_ref = (real *) calloc(count,sizeof(real));
    newcell.ort_ref = (real *) calloc(count*DIM, sizeof(real));
#endif
#ifdef REFPOS
    newcell.refpos = (real *) malloc(count*DIM*sizeof(real));
#endif
#ifdef NVX
    newcell.heatcond = (real *) malloc(count*sizeof(real));
#endif
#ifdef STRESS_TENS
    newcell.presstens = (real *) malloc(count*DIM*sizeof(real));
    newcell.presstens_offdia = (real *) malloc(count*DIM*sizeof(real));
#endif
#if (defined(COVALENT) && !defined(TWOD))
    newcell.neigh = (neightab **) malloc( count * sizeof(neighptr) );
    if (NULL == newcell.neigh) {
      error("COVALENT: cannot allocate neighbor tables");
    }
    for (i=0; i<thecell->n_max; ++i) {
      newcell.neigh[i] = thecell->neigh[i];
    }
    for (i=thecell->n_max; i<count; ++i) {
      newcell.neigh[i] = alloc_neightab(newcell.neigh[i], neigh_len);
    }
#endif
#ifdef UNIAX
    newcell.traeg_moment  = (real    * ) malloc(count * sizeof(real)    );
    newcell.shape = (real *) malloc(count*DIM*sizeof(real));
    newcell.pot_well = (real *) malloc(count*DIM*sizeof(real));
#endif
#endif /* not MONOLJ */

    if ((NULL == space)
#ifndef MONOLJ
        || (NULL == newcell.nummer)
        || (NULL == newcell.sorte)
        || (NULL == newcell.masse)
        || (NULL == newcell.pot_eng)
#ifdef EAM2
	|| (NULL == newcell.eam2_rho_h)
#endif
#ifdef ORDPAR
#ifndef TWOD
	|| (NULL == newcell.nbanz)
#endif
#endif
#ifdef DISLOC
        || (NULL == newcell.Epot_ref)
        || (NULL == newcell.ort_ref)
#endif
#ifdef REFPOS
        || (NULL == newcell.refpos)
#endif
#ifdef NVX
        || (NULL == newcell.heatcond)
#endif
#ifdef STRESS_TENS
        || (NULL == newcell.presstens)
        || (NULL == newcell.presstens_offdia)
#endif
#ifdef UNIAX
        || (NULL == newcell.traeg_moment)
	|| (NULL == newcell.shape)
	|| (NULL == newcell.pot_well)
#endif
#endif /* not MONOLJ */
        ) {
      printf("Want %d bytes.\n",newcellsize);
#ifdef CRAY
      malloc_stats(0);
#endif
      printf("Have %d atoms.\n",natoms);
      error("Cannot allocate memory for cell.");
    }
  }
  
  if (0 == thecell->n_max) {
    /* cell is just initialized */
    thecell->n = 0;
  } else {

    if (count < thecell->n_max) { /* cell shrinks, data is invalidated */
      thecell->n = 0;
#if (defined(COVALENT) && !defined(TWOD))
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
      memcpy(newcell.sorte ,  thecell->sorte,   thecell->n * sizeof(shortint));
      memcpy(newcell.masse ,  thecell->masse,   thecell->n * sizeof(real));
      memcpy(newcell.pot_eng, thecell->pot_eng, thecell->n * sizeof(real));
#ifdef EAM2
      memcpy(newcell.eam2_rho_h, thecell->eam2_rho_h, 
             thecell->n * sizeof(real));
#endif
#ifdef ORDPAR
#ifndef TWOD
      memcpy(newcell.nbanz, thecell->nbanz,
             thecell->n * sizeof(shortint));
#endif
#endif
#ifdef DISLOC
      memcpy(newcell.Epot_ref, thecell->Epot_ref, 
                               thecell->n * sizeof(real));
      memcpy(newcell.ort_ref,  thecell->ort_ref, 
                               thecell->n * DIM * sizeof(real));
#endif
#ifdef REFPOS
      memcpy(newcell.refpos, thecell->refpos, thecell->n * DIM * sizeof(real));
#endif
#ifdef NVX
      memcpy(newcell.heatcond,  thecell->heatcond,  thecell->n * sizeof(real));
#endif
#ifdef STRESS_TENS
      memcpy(newcell.presstens,  thecell->presstens,  
             thecell->n * DIM * sizeof(real));
      memcpy(newcell.presstens_offdia,  thecell->presstens_offdia,  
             thecell->n * DIM * sizeof(real));
#endif
#ifdef UNIAX
      memcpy(newcell.traeg_moment, thecell->traeg_moment, 
                               thecell->n * sizeof(real));
      memcpy(newcell.achse , thecell->achse,  thecell->n * DIM * sizeof(real));
      memcpy(newcell.shape, thecell->shape, 
                               thecell->n * DIM * sizeof(real));
      memcpy(newcell.pot_well, thecell->pot_well, 
                               thecell->n * DIM * sizeof(real));
      memcpy(newcell.dreh_impuls, thecell->dreh_impuls, 
                               thecell->n * DIM * sizeof(real));
      memcpy(newcell.dreh_moment, thecell->dreh_moment, 
                               thecell->n * DIM * sizeof(real));
#endif
#endif /* not MONOLJ */
    }

    /* deallocate old cell */
    free(thecell->ort);
#ifndef MONOLJ
    free(thecell->nummer);
    free(thecell->sorte);
    free(thecell->masse);
    free(thecell->pot_eng);
#ifdef EAM2
    free(thecell->eam2_rho_h);
#endif
#ifdef ORDPAR
#ifndef TWOD
    free(thecell->nbanz);
#endif
#endif
#ifdef DISLOC
    free(thecell->Epot_ref);
    free(thecell->ort_ref);
#endif
#ifdef REFPOS
    free(thecell->refpos);
#endif
#ifdef NVX
    free(thecell->heatcond);
#endif
#ifdef STRESS_TENS
    free(thecell->presstens);
    free(thecell->presstens_offdia);
#endif
#ifdef COVALENT
    free(thecell->neigh);
#endif
#ifdef UNIAX
    free(thecell->traeg_moment);
    free(thecell->shape);
    free(thecell->pot_well);
#endif
#endif /* not MONOLJ */
  }

  /* set pointers accordingly */
  thecell->ort    = newcell.ort;
  thecell->impuls = newcell.impuls;
  thecell->kraft  = newcell.kraft;
#ifndef MONOLJ
  thecell->nummer   = newcell.nummer;
  thecell->sorte    = newcell.sorte;
  thecell->masse    = newcell.masse;
  thecell->pot_eng  = newcell.pot_eng;
#ifdef EAM2
  thecell->eam2_rho_h = newcell.eam2_rho_h;
#endif
#ifdef ORDPAR
  thecell->nbanz = newcell.nbanz;
#endif
#ifdef DISLOC
  thecell->Epot_ref = newcell.Epot_ref;
  thecell->ort_ref = newcell.ort_ref;
#endif
#ifdef REFPOS
  thecell->refpos = newcell.refpos;
#endif
#ifdef NVX
  thecell->heatcond = newcell.heatcond;
#endif
#ifdef STRESS_TENS
  thecell->presstens = newcell.presstens;
  thecell->presstens_offdia = newcell.presstens_offdia;
#endif
#ifdef COVALENT
  thecell->neigh = newcell.neigh;
#endif
#ifdef UNIAX
  thecell->traeg_moment = newcell.traeg_moment;
  thecell->achse  = newcell.achse;
  thecell->dreh_impuls = newcell.dreh_impuls;
  thecell->dreh_moment = newcell.dreh_moment;
  thecell->shape = newcell.shape;
  thecell->pot_well = newcell.pot_well;
#endif
#endif /* not MONOLJ */

  thecell->n_max = count;

}
