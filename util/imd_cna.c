
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
* imd_cna -- perform Common-Neighbour Analysis
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#ifndef CNA
#define CNA
#endif

#define MAIN

#include "util.h"

/******************************************************************************
*
*  Usage
*
******************************************************************************/

void usage(void)
{ 
  printf("%s [-r <nnn>] [-A <nnn>] [-a <nnn>] [-e <nnn>] [-v] [-f] [-g] [-l <nnn> <nnn> <nnn>] [-n <nnn>] [-u <nnn> <nnn> <nnn>] [-w <nnn>] [-W <n> <n> <n> <n>] [-p paramter-file]\n", progname); 
  exit(1); 
}

/*****************************************************************************
*
*  main
*
*****************************************************************************/

int main(int argc, char **argv)
{
  /* Read command line arguments */
  read_command_line(argc,argv);

  /* Read Parameters from parameter file */
  read_parameters();

  /* read box from file header */
  if (box_from_header) read_box(infilename);

  /* Initializations */
  init_cna();

  /* Initialize cell data structures */
  init_cells();

  /* Read atoms */
  read_atoms(infilename);

  /* Calculate neighbour tables */
  do_work(do_neighbour_tables);

  /* Perform Common-Neighbour Analysis */
  do_work(do_cna);

  /* Sort pair types */
  sort_pair_types();

  /* Write statistics */
  write_statistics();

  /* Write atoms */
  if (writeatoms) write_atoms();

  /* Write decomposition of RDF */
  if (rdf) write_rdf();

  return 0;
}

/******************************************************************************
*
*  init_cna -- Initializations
*
******************************************************************************/

void init_cna()
{
  int  i, j, n = 0;
  real tmp;

  /* Cutoff radii for definition of nearest neighbour atoms */
  for (i=0; i<ntypes; i++) {
    for (j=i; j<ntypes; j++) {
      if ( r_max == -1.0 && r_cut_vec[0] != -1.0 )
	/* Cutoffs are given by parameter r_cut */
	tmp = r_cut_vec[n];
      else if ( r_max == -1.0 )
	/* Default value for r_cut */
	tmp = 1.0;
      else
	/* Take r_max given by the option -e */
	tmp = r_max;

      r_cut[i][j]  = r_cut[j][i]  = tmp;
      ++n;      
    }
  }

  /* Cutoff radius for cell decomposition */
  tmp = 0.0;
  for (i=0; i<ntypes; ++i)
    for (j=0; j<ntypes; ++j)
      tmp = MAX( tmp, r_cut[i][j]*r_cut[i][j] );
  r2_cut = MAX(r2_cut,tmp);
  if(!nearestneighbour)
    r2_cut *= 4.0;  /* r_cut is twice the cutoff length */
  rcut = sqrt(r2_cut);

  /* Pair type to be written out */
  if( pair_type_short<0 || pair_type_short>2999 )
    error("invalid pair_type!");
  else {
    pair_type.z = pair_type_short % 10;
    pair_type_short /= 10;
    pair_type.y = pair_type_short % 10;
    pair_type_short /= 10;
    pair_type.x = pair_type_short % 10;
    pair_type_short /= 10;
    pair_type.i = pair_type_short % 10;
  }
}

/******************************************************************************
*
*  do_neighbour_tables -- Calculates neighbour tables
*
******************************************************************************/

void do_neighbour_tables(cell *p, cell *q, vektor pbc) 
{
  int    i, j;
  int    jstart;
  int    q_typ, p_typ;
  vektor d, tmp_d;
  real   radius;

  /* For each atom in first cell */
  for (i=0; i<p->n; ++i) {

    tmp_d.x = p->ort[i].x - pbc.x;
    tmp_d.y = p->ort[i].y - pbc.y;
#ifndef TWOD
    tmp_d.z = p->ort[i].z - pbc.z;
#endif

    p_typ   = p->sorte[i];
    jstart = ( p==q ? i+1 : 0 );

    /* For each atom in neighbouring cell */
    for (j=jstart; j<q->n; ++j) {

      q_typ = q->sorte[j];
	  
      /* Calculate distance  */
      d.x = q->ort[j].x - tmp_d.x;
      d.y = q->ort[j].y - tmp_d.y;
#ifndef TWOD
      d.z = q->ort[j].z - tmp_d.z;
#endif

      radius = sqrt(SPROD(d,d));

      /* Make neighbour tables */
      if (radius <= r_cut[p_typ][q_typ]) { 

	neightab *neigh;
	int newsz;

	/* Update neighbour table of particle i */
	neigh = p->neightab_array + i;

	if (neigh->n_max <= neigh->n ) {
	  newsz = neigh->n_max + NSTEP;

	  neigh->typ = (short *) realloc(neigh->typ, newsz * sizeof(short) );
	  neigh->cl  = (void **) realloc(neigh->cl,  newsz * sizeof(cellptr) );
	  neigh->num = (int *)   realloc(neigh->num, newsz * sizeof(int) );

	  if ( neigh->typ==NULL || neigh->cl==NULL || neigh->num==NULL )
	    error("Cannot allocate memory for neighbour tables!");
	  
	  neigh->n_max += NSTEP;
	}	  
	neigh->typ[neigh->n] = q_typ;
	neigh->cl [neigh->n] = q;
	neigh->num[neigh->n] = j;
	neigh->n++;

	/* Update neighbour table of particle j */
	neigh      = q->neightab_array + j;
	
	if (neigh->n_max <= neigh->n ) {
	  newsz = neigh->n_max + NSTEP;
	  
	  neigh->typ = (short *) realloc(neigh->typ, newsz * sizeof(short) );
	  neigh->cl  = (void **) realloc(neigh->cl,  newsz * sizeof(cellptr) );
	  neigh->num = (int *)   realloc(neigh->num, newsz * sizeof(int) );

	  if ( neigh->typ==NULL || neigh->cl==NULL || neigh->num==NULL )
	    error("Cannot allocate memory for neighbour tables!");
	  
	  neigh->n_max      += NSTEP;
	}
	    
	neigh->typ[neigh->n] = p_typ;
	neigh->cl [neigh->n] = p;
	neigh->num[neigh->n] = i;
	neigh->n++;
      }
    } /* for j */
  } /* for i */
}

/******************************************************************************
*
*  do_cna -- Perform Common-Neighbour Analysis
*
******************************************************************************/

void do_cna(cell *p, cell *q, vektor pbc) {

  int       i, j, jstart, k, l, m;
  int       p_typ, q_typ;
  real      p_mass, q_mass;
  neightab  *ineigh, *jneigh, *neigh;
  int       cna_neigh, cna_atoms, cna_bonds, cna_chain, tmp_cna_chain;
  cell      *cna_cell[MAX_NEIGH];
  int       cna_num[MAX_NEIGH];
  cell      *cptr;
  int       start, end;
  int       type;
  ivektor4d pair_index;
  vektor    d;
  real      radius2, distance;
  int       index;

  /* For each atom in first cell */
  for (i=0; i<p->n; ++i) {

    p_typ   = p->sorte[i];
    p_mass  = p->masse[i];
    ineigh  = p->neightab_array + i;
    jstart  = ( p==q ? i+1 : 0 );

    /* For each atom in neighbouring cell */
    for (j=jstart; j<q->n; ++j) {


      if ( atom_in_pbox(p->ort[i]) || atom_in_pbox(q->ort[j]) ) {

	q_typ  = q->sorte[j];
	q_mass = q->masse[j];
	jneigh = q->neightab_array + j;
	
	cna_neigh = 2;
	cna_atoms = 0;
	cna_bonds = 0;
	cna_chain = 0; 
	
	d.x = q->ort[j].x - p->ort[i].x + pbc.x;
	d.y = q->ort[j].y - p->ort[i].y + pbc.y;
#ifndef TWOD
	d.z = q->ort[j].z - p->ort[i].z + pbc.z;
#endif
	radius2 = (double)(SPROD(d,d));

	if (radius2 < r2_cut) {
 
	  /* For all neighbours of atom i */
	  for (k=0; k<ineigh->n; k++) {
	    
	    /* Check whether j is neighbour of i */
	    if ( ineigh->cl[k] == q && ineigh->num[k] == j )
	      cna_neigh = 1;
	    
	    /* For all neighbours of atom j */
	    for (l=0; l<jneigh->n; l++) {
	      
	      /* Count number of common neighbours and store them */	  
	      if ( ineigh->cl[k] == jneigh->cl[l]
		   && ineigh->num[k] == jneigh->num[l] ) {
		if ( cna_atoms >= MAX_NEIGH )
		  error("Too many common neighbours");
		cna_cell[cna_atoms]  = ineigh->cl[k];
		cna_num[cna_atoms]   = ineigh->num[k];
		++cna_atoms;
	      }
	    } /* l */
	  } /* k */
	  
	  if( cna_atoms > 0 && ( !nearestneighbour || cna_neigh==1 ) ) {

	    /* Count total number of pairs */
	    ++cna_pairs;
	    
	    /* Count number of bonds and store them */
	    for (k=0; k<cna_atoms; k++) {
	      neigh = cna_cell[k]->neightab_array + cna_num[k];
	      for (l=k+1; l<cna_atoms; l++)
		for (m=0; m<neigh->n; m++) {
		  if ( cna_cell[l] == neigh->cl[m]
		       && cna_num[l] == neigh->num[m] ) {
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

	    pair_index.i = cna_neigh;
	    pair_index.x = cna_atoms;
	    pair_index.y = cna_bonds;
	    pair_index.z = cna_chain;

	    /* Mark atoms to be written out */
	    if ( writeatoms ) {	    
	      if ( pair_index.i    == pair_type.i
		   && pair_index.x == pair_type.x
		   && pair_index.y == pair_type.y
		   && pair_index.z == pair_type.z ) {
		if ( p->mark[i] == 0 )
		  p->mark[i] = 1;
		if ( q->mark[j] == 0 )
		  q->mark[j] = 1;
	      }
	    }

	    /* Count number of pairs of specific type */
	    if ( rdf ) {
	      distance = sqrt(radius2);
	      if ( distance > r_min )
		index = (int) ( slots * ( distance - r_min ) 
				/ ( rcut - r_min ) );
	    }
	    
	    if ( type_list_length == 0 ) {
	      type_list[type_list_length++] = pair_index;
	      type = 0;
	      if (rdf) {
		rdf_tab[type] = (int *) malloc(slots*sizeof(int));
		for(k=0; k<slots; k++)
		  rdf_tab[type][k] = 0;
	      }
	      type_count[type] = 0;
	      }
	    else {
	      type = -1;
	      for (k=0; k<type_list_length; k++)
		if ( type_list[k].i == pair_index.i
		     && type_list[k].x == pair_index.x
		     && type_list[k].y == pair_index.y
		     && type_list[k].z == pair_index.z ) {
		  type = k;
		  break;
		}
	    }

	    if ( type == -1 ) {
	      if (type_list_length>=MAX_TYPES)
		error("Too many pair types");
	      type_list[type_list_length++] = pair_index;
	      type = type_list_length - 1;
	      if (rdf) {
		rdf_tab[type] = (int *) malloc(slots*sizeof(int));
		for(k=0; k<slots; k++)
		  rdf_tab[type][k] = 0;
	      }
	      type_count[type]    = 0;
	    }
	    
	    if (rdf) 
	      ++rdf_tab[type][index];
	    ++type_count[type];
	    
	  } /* cna_atoms > 0 && ... */
	} /* radius < r_cut */
      }
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

/******************************************************************************
*
*  write_atoms -- Writes marked atoms to file .ijkl.cna
*
******************************************************************************/

void write_atoms(void)
{
  int i, j, k, l;
  cell *p;
  FILE *out;
  str255 fname;

  if (-1==restart)
    sprintf(fname,"%s.%d%d%d%d.cna",infilename,
	    pair_type.i, pair_type.x, pair_type.y, pair_type.z);
  else
    sprintf(fname,"%s.%05d.%d%d%d%d.cna",outfilename,restart,
	    pair_type.i, pair_type.x, pair_type.y, pair_type.z);
  
  out = fopen(fname,"w");
  if (NULL == out) error("Cannot open cna file.");

  /* For each cell */
  for (i=0; i < cell_dim.x; ++i)
    for (j=0; j < cell_dim.y; ++j)
#ifndef TWOD
      for (k=0; k < cell_dim.z; ++k)
#endif
      {
#ifdef TWOD
        p = PTR_2D_V(cell_array,i,j  ,cell_dim);
#else
        p = PTR_3D_V(cell_array,i,j,k,cell_dim);
#endif
	/* For each atom in cell */
	for (l=0; l<p->n; l++)
	  if (p->mark[l] == 1 )
	    fprintf(out, "%d %d %f %f %f %f\n", 
		    p->nummer[l], p->sorte[l], p->masse[l], 
		    p->ort[l].x, p->ort[l].y
#ifndef TWOD
		    , p->ort[l].z
#endif
		    );
      }
  fclose(out);
}

/******************************************************************************
*
*  write_statistics -- Writes statistics to standard output
*
******************************************************************************/

void write_statistics(void) 
{
  int i;

  printf("\nFound %d pairs.\n\n", cna_pairs);

  if (cna_pairs>0) {
    printf("  pair type   occurrence     relative occurrence\n");
    printf("--------------------------------------------------------\n");
    for (i=0; i<type_list_length; i++) 
      if(!nearestneighbour || type_list[i].i==1)
	printf("   %d%d%d%d    %10d         %3.2f %%\n",
	       type_list[type_sort[i]].i, type_list[type_sort[i]].x, 
	       type_list[type_sort[i]].y, type_list[type_sort[i]].z,   
	       type_count[type_sort[i]], 
	       100.0*(real)type_count[type_sort[i]]/cna_pairs);
    printf("\n");  
  }
}

/******************************************************************************
*
*  write_rdf -- Writes decomposition of RDF to file .rdf
*
******************************************************************************/

void write_rdf(void) 
{
  int i, j, total;
  static FILE *out;
  str255 fname;

  if (-1==restart)
    sprintf(fname,"%s.rdf",infilename);
  else
    sprintf(fname,"%s.%05d.rdf",outfilename,restart);

  out = fopen(fname,"w");
  if (NULL == out) error("Cannot open rdf file.");

  /* Header */
  fprintf(out, "#r total");
  for (i=0; i<type_list_length; i++)
    fprintf(out, " %d%d%d%d", 
	    type_list[type_sort[i]].i, type_list[type_sort[i]].x, 
	    type_list[type_sort[i]].y, type_list[type_sort[i]].z);
  fprintf(out, "\n");

  /* RDF */
  for (i=0; i<slots; i++) {
    fprintf(out, "%f", r_min + i * ( rcut - r_min ) / slots );
    total = 0;
    for (j=0; j<type_list_length; j++)
      total += rdf_tab[j][i];
    fprintf(out, " %f", (real)total/cna_pairs);
    for (j=0; j<type_list_length; j++)
      fprintf(out, " %f", (real) rdf_tab[type_sort[j]][i] / cna_pairs);
    fprintf(out, "\n");
  }

  fclose(out);
}
