
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2006 Institute for Theoretical and Applied Physics,
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
  to->ort X(i) = from->ort X(j); 
  to->ort Y(i) = from->ort Y(j); 
#ifndef TWOD
  to->ort Z(i) = from->ort Z(j); 
#endif
#ifndef MONOLJ
  to->nummer [i] = from->nummer [j]; 
  to->sorte  [i] = from->sorte  [j]; 
  to->vsorte [i] = from->vsorte [j]; 
  to->masse  [i] = from->masse  [j]; 
  to->pot_eng[i] = from->pot_eng[j];
#endif
#ifdef EAM2
  to->eam_rho[i] = from->eam_rho[j]; 
  to->eam_dF [i] = from->eam_dF [j]; 
#ifdef EEAM
  to->eeam_p_h[i] = from->eeam_p_h[j]; 
  to->eeam_dM [i] = from->eeam_dM [j]; 
#endif
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
#ifdef NVX
  to->heatcond[i] = from->heatcond[j];   
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

/******************************************************************************
*
*  increase already existing neighbor table for one particle
*
******************************************************************************/

void increase_neightab(neightab *neigh, int count)
{
  neigh->dist = (real *)     realloc( neigh->dist, count * DIM * sizeof(real));
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
*  Allocate memory for a cell. 
*
******************************************************************************/

void alloc_cell(cell *thecell, int count)
{
  void *space;
  cell newcell;
  int i, ncopy;
  real *tmp;
  int newcellsize;

  /* cell is freed */
  if (0==count) {
    thecell->n = 0;
#ifdef VEC
    thecell->n_buf = 0;
#endif
    newcell.ort = NULL;
    newcell.impuls = NULL;
    newcell.kraft = NULL;
#ifndef MONOLJ
    newcell.nummer = NULL;
    newcell.sorte  = NULL;
    newcell.vsorte = NULL;
    newcell.masse  = NULL;
    newcell.pot_eng= NULL;
#ifdef EAM2
    newcell.eam_rho = NULL;
    newcell.eam_dF  = NULL;
#ifdef EEAM
    newcell.eeam_p_h = NULL;
    newcell.eeam_dM  = NULL;
#endif
#endif
#ifdef DAMP
    newcell.damp_f = NULL;
#endif
#ifdef ADP
    newcell.adp_mu     = NULL;
    newcell.adp_lambda = NULL;
#endif
#ifdef CG
    newcell.h = NULL;
    newcell.g = NULL;
    newcell.old_ort = NULL;
#endif
#ifdef NNBR
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
#ifdef SHOCK
    newcell.pxavg = NULL;
#endif
#ifdef COVALENT
    newcell.neigh  = NULL;
#endif
#ifdef NBLIST
    newcell.nbl_pos  = NULL;
#endif
#ifdef UNIAX
    newcell.achse = NULL;
    newcell.dreh_impuls = NULL;
    newcell.dreh_moment = NULL;
#endif
#if defined(VEC) && defined(MPI)
    newcell.ind = NULL;
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
    newcell.vsorte = (shortint* ) malloc(count * sizeof(shortint));
    newcell.masse  = (real    * ) malloc(count * sizeof(real)    );
    newcell.pot_eng= (real    * ) malloc(count * sizeof(real)    );
#ifdef EAM2
    newcell.eam_rho = (real *) malloc(count * sizeof(real));
    newcell.eam_dF  = (real *) malloc(count * sizeof(real));
#ifdef EEAM
    newcell.eeam_p_h = (real *) malloc(count * sizeof(real));
    newcell.eeam_dM  = (real *) malloc(count * sizeof(real));
#endif
#endif
#ifdef DAMP
    newcell.damp_f = (real *) malloc(count * sizeof(real));
#endif
#ifdef ADP
    newcell.adp_mu     = (real       *) malloc(count * DIM * sizeof(real));
    newcell.adp_lambda = (sym_tensor *) malloc(count * sizeof(sym_tensor));
#endif
#ifdef NNBR
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
#ifdef SHOCK
    newcell.pxavg = (real *) malloc(count*sizeof(real));
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
#ifdef NBLIST
    newcell.nbl_pos = (real *) malloc(count*DIM*sizeof(real));
#endif
#if defined(VEC) && defined(MPI)
    newcell.ind = (integer *) malloc(count * sizeof(integer));
#endif
#endif /* not MONOLJ */

    if ((NULL == space)
#ifndef MONOLJ
        || (NULL == newcell.nummer)
        || (NULL == newcell.sorte)
        || (NULL == newcell.vsorte)
        || (NULL == newcell.masse)
        || (NULL == newcell.pot_eng)
#ifdef EAM2
	|| (NULL == newcell.eam_rho)
	|| (NULL == newcell.eam_dF)
#ifdef EEAM
	|| (NULL == newcell.eeam_p_h)
	|| (NULL == newcell.eeam_dM)
#endif
#endif
#ifdef DAMP
        || (NULL == newcell.damp_f)
#endif
#ifdef ADP
	|| (NULL == newcell.adp_mu)
	|| (NULL == newcell.adp_lambda)
#endif
#ifdef NNBR
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
#ifdef SHOCK
        || (NULL == newcell.pxavg)
#endif
#ifdef NBLIST
        || (NULL == newcell.nbl_pos)
#endif
#if defined(VEC) && defined(MPI)
        || (NULL == newcell.ind)
#endif
#endif /* not MONOLJ */
        ) {
      printf("Want %d bytes.\n",newcellsize);
#ifdef CRAY
      malloc_stats(0);
#endif
      printf("Have %ld atoms.\n",natoms);
      error("Cannot allocate memory for cell.");
    }
  }
  
  if (0 == thecell->n_max) {
    /* cell is just initialized */
    thecell->n = 0;
#ifdef VEC
    thecell->n_buf = 0;
#endif
  } else {

    if (count < thecell->n_max) { /* cell shrinks, data is invalidated */
      thecell->n = 0;
#ifdef VEC
      thecell->n_buf = 0;
#endif
#if (defined(COVALENT) && !defined(TWOD))
      /* deallocate all neighbor tables */
      for (i=0; i<thecell->n_max; ++i) {
        thecell->neigh[i] = alloc_neightab(thecell->neigh[i],0);
      }
      /* free(thecell->neigh); is freed later again */
#endif
    }

    /* if there are valid particles in cell, copy them to new cell */
#ifdef VEC
    ncopy = MAX(thecell->n,thecell->n_buf);
#else
    ncopy = thecell->n;
#endif
    if (ncopy > 0) {
      /* cell is enlarged, copy data from old to newcell location */
      memcpy(newcell.ort,     thecell->ort,     ncopy * DIM * sizeof(real));
      memcpy(newcell.impuls,  thecell->impuls,  ncopy * DIM * sizeof(real));
      memcpy(newcell.kraft,   thecell->kraft,   ncopy * DIM * sizeof(real));
#ifdef CG
      memcpy(newcell.h,       thecell->h,       ncopy * DIM * sizeof(real));
      memcpy(newcell.g,       thecell->g,       ncopy * DIM * sizeof(real));
      memcpy(newcell.old_ort, thecell->old_ort, ncopy * DIM * sizeof(real));
#endif
#ifndef MONOLJ
      memcpy(newcell.nummer,  thecell->nummer,  ncopy * sizeof(integer));
      memcpy(newcell.sorte,   thecell->sorte,   ncopy * sizeof(shortint));
      memcpy(newcell.vsorte,  thecell->vsorte,  ncopy * sizeof(shortint));
      memcpy(newcell.masse,   thecell->masse,   ncopy * sizeof(real));
      memcpy(newcell.pot_eng, thecell->pot_eng, ncopy * sizeof(real));
#ifdef EAM2
      memcpy(newcell.eam_rho, thecell->eam_rho, ncopy * sizeof(real));
      memcpy(newcell.eam_dF,  thecell->eam_dF,  ncopy * sizeof(real));
#ifdef EEAM
      memcpy(newcell.eeam_p_h, thecell->eeam_p_h, ncopy * sizeof(real));
      memcpy(newcell.eeam_dM,  thecell->eeam_dM,  ncopy * sizeof(real));
#endif
#endif
#ifdef DAMP
      memcpy(newcell.damp_f, thecell->damp_f, ncopy * sizeof(real));
#endif
#ifdef ADP
      memcpy(newcell.adp_mu,     thecell->adp_mu,  ncopy * DIM * sizeof(real));
      memcpy(newcell.adp_lambda, thecell->adp_lambda,ncopy*sizeof(sym_tensor));
#endif
#ifdef NNBR
      memcpy(newcell.nbanz, thecell->nbanz, ncopy * sizeof(shortint));
#endif
#ifdef DISLOC
      memcpy(newcell.Epot_ref, thecell->Epot_ref, ncopy * sizeof(real));
      memcpy(newcell.ort_ref,  thecell->ort_ref,  ncopy * DIM * sizeof(real));
#endif
#ifdef AVPOS
      memcpy(newcell.av_epot, thecell->av_epot, ncopy * sizeof(real));
      memcpy(newcell.avpos,   thecell->avpos, ncopy * DIM * sizeof(real));
      memcpy(newcell.sheet,   thecell->sheet, ncopy * DIM * sizeof(real));
#endif
#ifdef REFPOS
      memcpy(newcell.refpos, thecell->refpos, ncopy * DIM * sizeof(real));
#endif
#ifdef NVX
      memcpy(newcell.heatcond,  thecell->heatcond,  ncopy * sizeof(real));
#endif
#ifdef STRESS_TENS
      memcpy(newcell.presstens, thecell->presstens, ncopy*sizeof(sym_tensor));
#endif
#ifdef SHOCK
      memcpy(newcell.pxavg,  thecell->pxavg,  ncopy * sizeof(real));
#endif
#ifdef NBLIST
      memcpy(newcell.nbl_pos,  thecell->nbl_pos,  ncopy * sizeof(real));
#endif
#ifdef UNIAX
      memcpy(newcell.achse ,   thecell->achse,    ncopy * DIM * sizeof(real));
      memcpy(newcell.dreh_impuls, thecell->dreh_impuls,ncopy*DIM*sizeof(real));
      memcpy(newcell.dreh_moment, thecell->dreh_moment,ncopy*DIM*sizeof(real));
#endif
#if defined(VEC) && defined(MPI)
      memcpy(newcell.ind, thecell->ind, ncopy * sizeof(integer));
#endif
#endif /* not MONOLJ */
    }

    /* deallocate old cell */
    free(thecell->ort);
#ifndef MONOLJ
    free(thecell->nummer);
    free(thecell->sorte);
    free(thecell->vsorte);
    free(thecell->masse);
    free(thecell->pot_eng);
#ifdef EAM2
    free(thecell->eam_rho);
    free(thecell->eam_dF);
#ifdef EEAM
    free(thecell->eeam_p_h);
    free(thecell->eeam_dM);
#endif
#endif
#ifdef DAMP
    free(thecell->damp_f);
#endif
#ifdef ADP
    free(thecell->adp_mu);
    free(thecell->adp_lambda);
#endif
#ifdef NNBR
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
#ifdef SHOCK
    free(thecell->pxavg);
#endif
#ifdef COVALENT
    free(thecell->neigh);
#endif
#ifdef NBLIST
    free(thecell->nbl_pos);
#endif
#if defined(VEC) && defined(MPI)
    free(thecell->ind);
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
  thecell->vsorte   = newcell.vsorte;
  thecell->masse    = newcell.masse;
  thecell->pot_eng  = newcell.pot_eng;
#ifdef EAM2
  thecell->eam_rho = newcell.eam_rho;
  thecell->eam_dF  = newcell.eam_dF;
#ifdef EEAM
  thecell->eeam_p_h = newcell.eeam_p_h;
  thecell->eeam_dM  = newcell.eeam_dM;
#endif
#endif
#ifdef DAMP
  thecell->damp_f = newcell.damp_f;
#endif
#ifdef ADP
  thecell->adp_mu     = newcell.adp_mu;
  thecell->adp_lambda = newcell.adp_lambda;
#endif
#ifdef NNBR
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
#ifdef SHOCK
  thecell->pxavg = newcell.pxavg;
#endif
#ifdef COVALENT
  thecell->neigh = newcell.neigh;
#endif
#ifdef NBLIST
  thecell->nbl_pos = newcell.nbl_pos;
#endif
#ifdef UNIAX
  thecell->achse  = newcell.achse;
  thecell->dreh_impuls = newcell.dreh_impuls;
  thecell->dreh_moment = newcell.dreh_moment;
#endif
#if defined(VEC) && defined(MPI)
  thecell->ind = newcell.ind;
#endif
#endif /* not MONOLJ */

  thecell->n_max = count;

}
