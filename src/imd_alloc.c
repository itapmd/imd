
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
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

#include "imd.h"

/*****************************************************************************
*
*
*  Moves an atom from one cell to another. 
*  Does not move atoms between CPUs!
*  Does not really belong to file imd_alloc.c either!
*
******************************************************************************/

void move_atom(cell *to, cell *from, int index)
{
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
#ifndef MONOLJ
  to->nummer [to->n] = from->nummer[index]; 
  to->sorte  [to->n] = from->sorte[index]; 
  to->masse  [to->n] = from->masse[index]; 
  to->pot_eng[to->n] = from->pot_eng[index];
#endif
#ifdef EAM2
  to->eam2_rho_h[to->n] = from->eam2_rho_h[index]; 
#endif
#ifdef CG
  to->h  X(to->n) = from->h X(index); 
  to->h  Y(to->n) = from->h Y(index); 
#ifndef TWOD
  to->h  Z(to->n) = from->h Z(index);
#endif 
  to->g  X(to->n) = from->g X(index); 
  to->g  Y(to->n) = from->g Y(index); 
#ifndef TWOD
  to->g  Z(to->n) = from->g Z(index);
#endif  
  to->old_ort  X(to->n) = from->old_ort X(index); 
  to->old_ort  Y(to->n) = from->old_ort Y(index); 
#ifndef TWOD
  to->old_ort  Z(to->n) = from->old_ort Z(index); 
#endif
#endif
#ifdef DISLOC
  to->Epot_ref  [to->n] = from->Epot_ref[index];
  to->ort_ref X (to->n) = from->ort_ref X(index);
  to->ort_ref Y (to->n) = from->ort_ref Y(index);
#ifndef TWOD
  to->ort_ref Z (to->n) = from->ort_ref Z(index);
#endif
#endif /* DISLOC */
#ifdef AVPOS
  to->av_epot[to->n] = from->av_epot[index];
  to->avpos X(to->n) = from->avpos X(index);
  to->avpos Y(to->n) = from->avpos Y(index);
  to->sheet X(to->n) = from->sheet X(index);
  to->sheet Y(to->n) = from->sheet Y(index);
#ifndef TWOD
  to->avpos Z(to->n) = from->avpos Z(index);
  to->sheet Z(to->n) = from->sheet Z(index);
#endif
#endif /* AVPOS */
#ifdef ORDPAR
  to->nbanz[to->n] = from->nbanz[index]; 
#endif
#ifdef REFPOS
  to->refpos X(to->n) = from->refpos X(index);
  to->refpos Y(to->n) = from->refpos Y(index);
#ifndef TWOD
  to->refpos Z(to->n) = from->refpos Z(index);
#endif
#endif /* REFPOS */
#ifdef NVX
  to->head_cond[to->n] = from->head_cond[index];   
#endif
#ifdef STRESS_TENS
  to->presstens[to->n].xx = from->presstens[index].xx;   
  to->presstens[to->n].yy = from->presstens[index].yy;   
  to->presstens[to->n].xy = from->presstens[index].xy;   
#ifndef TWOD
  to->presstens[to->n].zz = from->presstens[index].zz;   
  to->presstens[to->n].yz = from->presstens[index].yz;   
  to->presstens[to->n].zx = from->presstens[index].zx;   
#endif
#endif /* STRESS_TENS */
  to->impuls X(to->n) = from->impuls X(index); 
  to->impuls Y(to->n) = from->impuls Y(index); 
#ifndef TWOD
  to->impuls Z(to->n) = from->impuls Z(index); 
#endif
  to->kraft X(to->n) = from->kraft X(index); 
  to->kraft Y(to->n) = from->kraft Y(index); 
#ifndef TWOD
  to->kraft Z(to->n) = from->kraft Z(index); 
#endif
#ifdef COVALENT
  /* neighbor table is not copied */
#endif
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
  to->dreh_impuls X(to->n) = from->dreh_impuls X(index); 
  to->dreh_impuls Y(to->n) = from->dreh_impuls Y(index); 
  to->dreh_impuls Z(to->n) = from->dreh_impuls Z(index); 
  to->dreh_moment X(to->n) = from->dreh_moment X(index); 
  to->dreh_moment Y(to->n) = from->dreh_moment Y(index); 
  to->dreh_moment Z(to->n) = from->dreh_moment Z(index); 
#endif

  ++to->n;

  /* Delete atom in original cell */

  --from->n;

  if (0 < from->n) {

    from->ort X(index) = from->ort X(from->n); 
    from->ort Y(index) = from->ort Y(from->n); 
#ifndef TWOD
    from->ort Z(index) = from->ort Z(from->n); 
#endif
#ifndef MONOLJ
    from->nummer [index] = from->nummer[from->n]; 
    from->sorte  [index] = from->sorte[from->n]; 
    from->masse  [index] = from->masse[from->n]; 
    from->pot_eng[index] = from->pot_eng[from->n];
#endif
#ifdef EAM2
    from->eam2_rho_h[index] = from->eam2_rho_h[from->n];
#endif
#ifdef CG
  from->h X(index) = from->h X(from->n); 
  from->h Y(index) = from->h Y(from->n); 
#ifndef TWOD
  from->h Z(index) = from->h Z(from->n); 
#endif 
  from->g X(index) = from->g X(from->n); 
  from->g Y(index) = from->g Y(from->n); 
#ifndef TWOD
  from->g Z(index) = from->g Z(from->n); 
#endif  
  from->old_ort X(index) = from->old_ort X(from->n); 
  from->old_ort Y(index) = from->old_ort Y(from->n); 
#ifndef TWOD
  from->old_ort Z(index) = from->old_ort Z(from->n); 
#endif
#endif
#ifdef DISLOC
    from->Epot_ref [index] = from->Epot_ref[from->n];
    from->ort_ref X(index) = from->ort_ref X(from->n);
    from->ort_ref Y(index) = from->ort_ref Y(from->n);
#ifndef TWOD
    from->ort_ref Z(index) = from->ort_ref Z(from->n);
#endif
#endif /* DISLOC */
#ifdef AVPOS
    from->av_epot[index] = from->av_epot[from->n];
    from->avpos X(index) = from->avpos X(from->n);
    from->avpos Y(index) = from->avpos Y(from->n);
    from->sheet X(index) = from->sheet X(from->n);
    from->sheet Y(index) = from->sheet Y(from->n);
#ifndef TWOD
    from->avpos Z(index) = from->avpos Z(from->n);
    from->sheet Z(index) = from->sheet Z(from->n);
#endif
#endif /* AVPOS */
#ifdef ORDPAR
    from->nbanz[index] = from->nbanz[from->n];
#endif
#ifdef REFPOS
    from->refpos X(index) = from->refpos X(from->n);
    from->refpos Y(index) = from->refpos Y(from->n);
#ifndef TWOD
    from->refpos Z(index) = from->refpos Z(from->n);
#endif
#endif /* REFPOS */
#ifdef NVX
    from->head_cond[index] = from->head_cond[from->n];
#endif
#ifdef STRESS_TENS
    from->presstens[index].xx = from->presstens[from->n].xx;   
    from->presstens[index].yy = from->presstens[from->n].yy;   
    from->presstens[index].xy = from->presstens[from->n].xy;   
#ifndef TWOD
    from->presstens[index].zz = from->presstens[from->n].zz;   
    from->presstens[index].yz = from->presstens[from->n].yz;   
    from->presstens[index].zx = from->presstens[from->n].zx;   
#endif
#endif /* STRESS_TENS */
    from->impuls X(index) = from->impuls X(from->n); 
    from->impuls Y(index) = from->impuls Y(from->n); 
#ifndef TWOD
    from->impuls Z(index) = from->impuls Z(from->n); 
#endif
    from->kraft X(index) = from->kraft X(from->n); 
    from->kraft Y(index) = from->kraft Y(from->n); 
#ifndef TWOD
    from->kraft Z(index) = from->kraft Z(from->n); 
#endif
#ifdef COVALENT
    /* neighbor table is not copied */
#endif
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
    from->dreh_impuls X(index) = from->dreh_impuls X(from->n); 
    from->dreh_impuls Y(index) = from->dreh_impuls Y(from->n); 
    from->dreh_impuls Z(index) = from->dreh_impuls Z(from->n); 
    from->dreh_moment X(index) = from->dreh_moment X(from->n); 
    from->dreh_moment Y(index) = from->dreh_moment Y(from->n); 
    from->dreh_moment Z(index) = from->dreh_moment Z(from->n); 
#endif
  }
}

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
#ifdef CG
    newcell.h = NULL;
    newcell.g = NULL;
    newcell.old_ort = NULL;
#endif
#ifdef ORDPAR
    newcell.nbanz = NULL;
#endif
#ifdef DISLOC
    newcell.Epot_ref = NULL;
    newcell.ort_ref = NULL;
#endif
#ifdef AVPOS
    newcell.av_epot = NULL;
    newcell.avpos   = NULL;
    newcell.sheet   = NULL;
#endif
#ifdef REFPOS
    newcell.refpos = NULL;
#endif
#ifdef NVX
    newcell.heatcond = NULL;
#endif
#ifdef STRESS_TENS
    newcell.presstens = NULL;
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
          DIM * sizeof(real)   /* kraft */ 
#ifdef CG
	  + DIM * sizeof(real) /* h */
	  + DIM * sizeof(real) /* g */
	  + DIM * sizeof(real) /* old_ort */
#endif
	    ); 
#endif  

        /* Get some space */
    space = malloc(newcellsize);


        /* Calculate Pointers */
  
    tmp = (real *) space;
    newcell.ort = tmp; tmp += DIM * count;
    newcell.impuls = tmp; tmp += DIM * count;
    newcell.kraft  = tmp; tmp += DIM * count;
#ifdef CG
    newcell.h = tmp; tmp += DIM * count;
    newcell.g = tmp; tmp += DIM * count;
    newcell.old_ort = tmp; tmp += DIM * count;
#endif
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
    newcell.nbanz = (shortint *) malloc(count * sizeof(shortint));
#endif
#ifdef DISLOC
    newcell.Epot_ref = (real *) calloc(count,sizeof(real));
    newcell.ort_ref = (real *) calloc(count*DIM, sizeof(real));
#endif
#ifdef AVPOS
    newcell.av_epot = (real *) calloc(count,     sizeof(real));
    newcell.avpos   = (real *) calloc(count*DIM, sizeof(real));
    newcell.sheet   = (real *) malloc(count * DIM * sizeof(real));
#endif
#ifdef REFPOS
    newcell.refpos = (real *) malloc(count*DIM*sizeof(real));
#endif
#ifdef NVX
    newcell.heatcond = (real *) malloc(count*sizeof(real));
#endif
#ifdef STRESS_TENS
    newcell.presstens = (sym_tensor *) malloc(count*sizeof(sym_tensor));
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
	|| (NULL == newcell.nbanz)
#endif
#ifdef DISLOC
        || (NULL == newcell.Epot_ref)
        || (NULL == newcell.ort_ref)
#endif
#ifdef AVPOS
	|| (NULL == newcell.av_epot)
        || (NULL == newcell.avpos)
        || (NULL == newcell.sheet)
#endif
#ifdef REFPOS
        || (NULL == newcell.refpos)
#endif
#ifdef NVX
        || (NULL == newcell.heatcond)
#endif
#ifdef STRESS_TENS
        || (NULL == newcell.presstens)
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
#ifdef CG
      memcpy(newcell.h   , thecell->h,    thecell->n * DIM * sizeof(real));
      memcpy(newcell.g   , thecell->g,    thecell->n * DIM * sizeof(real));
      memcpy(newcell.old_ort   , thecell->old_ort,    thecell->n * DIM * sizeof(real));
#endif
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
      memcpy(newcell.nbanz, thecell->nbanz,
             thecell->n * sizeof(shortint));
#endif
#ifdef DISLOC
      memcpy(newcell.Epot_ref, thecell->Epot_ref, 
                               thecell->n * sizeof(real));
      memcpy(newcell.ort_ref,  thecell->ort_ref, 
                               thecell->n * DIM * sizeof(real));
#endif
#ifdef AVPOS
      memcpy(newcell.av_epot, thecell->av_epot, 
                              thecell->n * sizeof(real));
      memcpy(newcell.avpos,   thecell->avpos, 
                              thecell->n * DIM * sizeof(real));
      memcpy(newcell.sheet,   thecell->sheet, thecell->n * DIM * sizeof(real));
#endif
#ifdef REFPOS
      memcpy(newcell.refpos, thecell->refpos, thecell->n * DIM * sizeof(real));
#endif
#ifdef NVX
      memcpy(newcell.heatcond,  thecell->heatcond,  thecell->n * sizeof(real));
#endif
#ifdef STRESS_TENS
      memcpy(newcell.presstens,  thecell->presstens,  
             thecell->n * sizeof(sym_tensor));
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
    free(thecell->nbanz);
#endif
#ifdef DISLOC
    free(thecell->Epot_ref);
    free(thecell->ort_ref);
#endif
#ifdef AVPOS
    free(thecell->av_epot);
    free(thecell->avpos);
    free(thecell->sheet);
#endif
#ifdef REFPOS
    free(thecell->refpos);
#endif
#ifdef NVX
    free(thecell->heatcond);
#endif
#ifdef STRESS_TENS
    free(thecell->presstens);
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
#ifdef CG
  thecell->h    = newcell.h;
  thecell->g    = newcell.g;
  thecell->old_ort    = newcell.old_ort;
#endif
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
  thecell->ort_ref  = newcell.ort_ref;
#endif
#ifdef AVPOS
  thecell->av_epot = newcell.av_epot;
  thecell->avpos   = newcell.avpos;
  thecell->sheet   = newcell.sheet;
#endif
#ifdef REFPOS
  thecell->refpos = newcell.refpos;
#endif
#ifdef NVX
  thecell->heatcond = newcell.heatcond;
#endif
#ifdef STRESS_TENS
  thecell->presstens = newcell.presstens;
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

























