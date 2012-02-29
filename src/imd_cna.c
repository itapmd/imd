
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2007 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_cna.c -- Routines for Common-Neighbour Analysis
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/******************************************************************************
*
*  do_cna -- Perform Common-Neighbour Analysis
*
******************************************************************************/

#ifdef NBLIST

void do_cna(void)
{
  vektor pbc = {0.0, 0.0, 0.0}; /* atoms in buffer cells have pbc applied */
  int c, c1, c2, m;
#ifdef BBOOST
    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("passing "" do_cna 1 start, without ipbc"" checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("    ************************* \n");fflush(stdout);
#endif
  for (c=0; c<ncells; c++) {
    c1 = cnbrs[c].np;
    for (m=0; m<14; m++) {
      c2 = cnbrs[c].nq[m];
      if (c2<0) continue;
      do_cna_func(cell_array + c1, cell_array + c2, pbc);
    }
  }
  /* collect mark variables */
  send_forces(add_mark,pack_mark,unpack_add_mark);
#ifdef BBOOST
    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("passing "" do_cna 1 end, without ipbc"" checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("    ************************* \n");fflush(stdout);
#endif

}

#else

void do_cna(void)
{
  vektor pbc;
  int n, k;
    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("passing "" do_cna 2 start, with ipbc"" checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("    ************************* \n");fflush(stdout); 
  for (n=0; n<nlists; ++n) {
    for (k=0; k<npairs[n]; ++k) {
      pair *P = pairs[n] + k;
      pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
      pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
      pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
      do_cna_func(cell_array + P->np, cell_array + P->nq, pbc);
    }
  }
#ifdef MPI
  /* collect mark variables */
  send_forces(add_mark,pack_mark,unpack_add_mark);
#endif
    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("passing "" do_cna 2 end, with ipbc"" checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("    ************************* \n");fflush(stdout); 
}

#endif

void do_cna_func(cell *p, cell *q, vektor pbc) {

  int       i, j, jstart, k, l, m;
  neightab  *ineigh, *neigh;
  vektor    tmp_d;
  real      *qptr;
  int       cna_neigh, cna_atoms, cna_bonds, cna_chain, tmp_cna_chain;
  cell      *cna_cell[MAX_NEIGH];
  int       cna_num[MAX_NEIGH];
  vektor    cna_d[MAX_NEIGH];
  cell      *cptr;
  int       start, end;
  int       type;
  int       pair_index;
  vektor    d, dj;
  static vektor *di = NULL;
  static int    curr_len = 0;
  real      *tmpptr, dj2;
  real      r2, distance;
  int       index;
  vektor    dk, dlk;
  real      dlk2;
  real      EPS = 0.001;

  if (curr_len < neigh_len) {
    di = (vektor *) realloc( di, neigh_len * sizeof(vektor) );
   if ( di==NULL )
      error("cannot allocate memory for temporary neighbor data");
    curr_len = neigh_len;
  }

  /* For each atom in first cell */
  for (i=0; i<p->n; ++i) {

    tmp_d.x = ORT(p,i,X) - pbc.x;
    tmp_d.y = ORT(p,i,Y) - pbc.y;
    tmp_d.z = ORT(p,i,Z) - pbc.z;

    ineigh  = NEIGH(p,i);

    jstart = (((p==q) && (pbc.x==0) && (pbc.y==0) && (pbc.z==0)) ? i+1 : 0);

    qptr   = &ORT(q,jstart,X);

    /* For each atom in neighbouring cell */
    for (j=jstart; j<q->n; ++j) {

      /* calculate distance */
      d.x = *qptr - tmp_d.x; ++qptr;
      d.y = *qptr - tmp_d.y; ++qptr;
      d.z = *qptr - tmp_d.z; ++qptr;
	
      cna_neigh = 2;
      cna_atoms = 0;
      cna_bonds = 0;
      cna_chain = 0; 
	
      /* distance between neighbours */
      r2 = SPROD(d,d);

      /* consider only first neighbours */
      if (r2 < cna_r2cut) {

	tmpptr = ineigh->dist;

	/* For all neighbours of atom i */
	for (k=0; k<ineigh->n; k++) {
	    
	  /* Check whether j is neighbour of i */
	  /* ist nur notwendig, wenn auch zweite Nachbarn betrachtet werden */
	  if ( ineigh->cl[k] == q && ineigh->num[k] == j )
	    cna_neigh = 1;

	  /* distance between i and k */
	  di[k].x    = *tmpptr++;
	  di[k].y    = *tmpptr++;
	  di[k].z    = *tmpptr++;	  

	  /* Check whether k is neighbour of j */
	  /* distance between j and k */
	  dj.x = di[k].x - d.x;
	  dj.y = di[k].y - d.y;
	  dj.z = di[k].z - d.z;

	  dj2 = SPROD(dj,dj);

	  /* Count number of common neighbours and store them */
	  if ( dj2 < cna_r2cut && dj2 > EPS) {

	    if ( cna_atoms >= MAX_NEIGH )
	      error("Too many common neighbours");
	    cna_cell[cna_atoms]  = ineigh->cl[k];
	    cna_num[cna_atoms]   = ineigh->num[k];
	    cna_d[cna_atoms].x   = di[k].x;
	    cna_d[cna_atoms].y   = di[k].y;
	    cna_d[cna_atoms].z   = di[k].z;
	    ++cna_atoms;
	  }
	} /* k */

	if( cna_atoms > 0 ) {

	  /* Count total number of pairs */
	  ++cna_pairs;
	    
	  /* Count number of bonds and store them */
	  for (k=0; k<cna_atoms; k++) {
	    dk.x = cna_d[k].x;
	    dk.y = cna_d[k].y;
	    dk.z = cna_d[k].z;

	    for (l=k+1; l<cna_atoms; l++) {

	      dlk.x = dk.x - cna_d[l].x; 
	      dlk.y = dk.y - cna_d[l].y; 
	      dlk.z = dk.z - cna_d[l].z; 
	      dlk2  = SPROD(dlk,dlk);

	      if ( dlk2 < cna_r2cut ) {

		if (cna_bonds>=MAX_BONDS)
		  error("Too many bonds");
		bondlist[cna_bonds][0] = k;
		bondlist[cna_bonds][1] = l;
		bondlist[cna_bonds][2] = 0;
		cna_bonds++;
	      }
	    }
	  }
	      
	  /* Count longest continuous chain of bonds */
	  for (k=0; k<cna_bonds; k++) {
	      
	    /* Initialize bond data */
	    start          = bondlist[k][0];
	    end            = bondlist[k][1];
	    for (l=0; l<cna_bonds; l++)
	      bondlist[l][2] = 0;
	    bondlist[k][2] = 1;
	      
	    tmp_cna_chain = 1;
	    cna_chain = MAX(cna_chain,tmp_cna_chain);
	    if (cna_chain==cna_bonds) 
	      break;
	      
	    /* Add further bonds to start bond recursively */
	    domino(start, end, cna_bonds, &cna_chain, &tmp_cna_chain);
	      
	    if (cna_chain==cna_bonds) 
	      break;		
	      
	  }

	  /* convert pair_type into integer form */
	  if (cna_atoms<10 && cna_bonds<10)
	    pair_index = ((cna_neigh*10+cna_atoms)*10+cna_bonds)*10+cna_chain;
	  else
	    pair_index = -1;

	  /* count pair types according to crystallinity */
	  if ( cna_crist > 0 ) {
	    if ( pair_index == 1421 ) {
	      MARK(p,i) += 10000;
	      MARK(q,j) += 10000;
	    }
	    else if ( pair_index == 1422 ) {
	      MARK(p,i) += 100;
	      MARK(q,j) += 100;
	    }
	    else {
	      ++MARK(p,i);
	      ++MARK(q,j);
	    }
	  }

	  /* Mark atoms to be written out */
	  if ( cna_write_n > 0 ) {
	    l = 1;
	    for(k=0; k<cna_write_n; k++) {
	      if ( pair_index == cna_writev[k] ) {
		if ( MARK(p,i)%(2*l) < l ) MARK(p,i) += l;
		if ( MARK(q,j)%(2*l) < l ) MARK(q,j) += l;
	      }
	      l *= 2;
	    }
	  }

	  /* Count number of pairs of specific type */
	  if (cna_write_statistics) {

	    if ( type_list_length == 0 ) {
	      type_list[type_list_length++] = pair_index;
	      type = 0;
	      type_count[type] = 0;
	    }
	    else {
	      type = -1;
	      for (k=0; k<type_list_length; k++)
		if ( type_list[k] == pair_index ) {
		  type = k;
		  break;
		}
	    }
	    
	    if ( type == -1 ) {
	      if (type_list_length>=MAX_TYPES)
		error("Too many pair types");
	      type_list[type_list_length++] = pair_index;
	      type = type_list_length - 1;
	      type_count[type]    = 0;
	    }
	    
	    ++type_count[type];
	  } 

	} /* cna_atoms > 0 && ... */
      } /* radius < r_cut */
      
    } /* j */
  } /* i */
}

/******************************************************************************
*
*  domino -- Adds a bond to the chain of bonds
*
******************************************************************************/

void domino(int start, int end, int listlength, int *max_chain, int *chain) 
{
  int i, start_old, end_old;

  /* Check all unused bonds */
  for(i=0; i<listlength; i++)
    if ( bondlist[i][2]==0 ) { 

      start_old = start;
      end_old   = end;
      
      if (      bondlist[i][0] == start )
	start = bondlist[i][1];
      else if ( bondlist[i][0] == end )
	end   = bondlist[i][1];
      else if ( bondlist[i][1] == start )
	start = bondlist[i][0];
      else if ( bondlist[i][1] == end )
	end   = bondlist[i][0];
      else
	continue;
      
      /* If a bond is found, remove it from the list of bonds */
      /* and invoke domino recursively */
      
      /* Update bond data */
      bondlist[i][2] = 1;
      ++(*chain);

      *max_chain = MAX(*max_chain,*chain);      
      if (*max_chain==listlength)
	break;
      
      domino(start, end, listlength, max_chain, chain);
      
      /* Reset bond data */
      --(*chain);
      start = start_old;
      end   = end_old;
      bondlist[i][2] = 0;
    }
}

/******************************************************************************
*
*  init_cna
*
******************************************************************************/

void init_cna(void)
{
  int k;
    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("passing "" init_cna start"" checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("    ************************* \n");fflush(stdout); 
  /* default value of cna_end */
  if (0==cna_end) cna_end = steps_max;
  cna_r2cut = cna_rcut * cna_rcut;
  /* update neighbor table cutoff */
  if (NULL==neightab_r2cut) {
    neightab_r2cut = (real *) calloc( ntypes * ntypes, sizeof(real) );
    if (NULL==neightab_r2cut) 
      error("cannot allocate memory for neightab_r2cut");
  }
  /* deactivate computation of neighbour tables */
  for(k=0; k<ntypes*ntypes; k++)
    neightab_r2cut[k] = -1.0;

  /* binary encoding of types of crystallinity to be written out */
  for(k=0; k<cna_crist_n; k++)
    cna_crist += ( 1 << cna_cristv[k] );

  /* if crystallinity is studied do not write out cna pairs */
  if ( cna_crist > 0 )
    cna_write_n = 0;
    printf("    ************************* \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("passing "" init_cna end"" checking by Lo! \n");fflush(stdout);
    printf("********************************* \n");fflush(stdout);
    printf("    ************************* \n");fflush(stdout); 
}

/******************************************************************************
*
*  write_statistics -- Writes statistics to standard output
*
******************************************************************************/

void write_statistics(void) 
{
  int i;

  printf("\nCNA: Found %d pairs.\n\n", cna_pairs);

  if (cna_pairs>0) {
    printf("  pair type   occurrence     relative occurrence\n");
    printf("--------------------------------------------------------\n");
    for (i=0; i<type_list_length; i++) 
      if(1)
	printf("   %d          %10d         %3.2f %%\n",
	       type_list[type_sort[i]],   
	       type_count[type_sort[i]], 
	       100.0*(real)type_count[type_sort[i]]/cna_pairs);
    printf("\n");  
  }
}

/******************************************************************************
*
*  sort_pair_types -- Sorts pair types according to number of occurrence
*
******************************************************************************/
 
void sort_pair_types(void)
{
  int i, j, tmp;

  for(i=0; i<type_list_length; i++)
    type_sort[i] = i;

  for(i=0; i<type_list_length; i++)
    for(j=0; j<type_list_length-1; j++)
      if(type_count[type_sort[j+1]]>type_count[type_sort[j]]) {
	tmp            = type_sort[j+1];
	type_sort[j+1] = type_sort[j];
	type_sort[j]   = tmp;
      }
}
