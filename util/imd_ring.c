
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
* imd_ring -- calculate ring statistics
*
* uses cell division routines of imd_pair
*
* Algorithm taken from D.S.Franzblau, Phys.Rev. B44 (1991) 4925
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#ifndef RING
#define RING
#endif

#ifdef TERSOFF
#undef TWOD
#endif

#define MAIN

#include "util.h"

/******************************************************************************
*
*  Usage -- educate users
*
******************************************************************************/

void usage(void)
{ 
  printf("%s [-r<nnn>] [-A<nnn>] [-e<nnn>] [-l<nnn>] [-v] [-p paramter-file]\n",progname); 
  
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

#ifdef TERSOFF
  /* Compute Tersoff parameters */
  init_tersoff();
#endif

  /* Initialize cell data structures */
  init_cells();

  /* Read atoms */
  read_atoms(infilename);

  /* Compute neighbour tables */
  do_work(do_neighbour_tables); 
  first = 0;
  do_work(do_neighbour_tables); 

  /* Search for rings */
  search_rings();

  /* Output ring statistics */
  write_data();

  return 0;

}

#ifdef TERSOFF

/******************************************************************************
*
*  init_tersoff -- initialize cutoff radii
*
******************************************************************************/

void init_tersoff(void) {

  int i, j, n = 0;
  real tmp;

  /* Cutoff radii for more than one atom type */
  for (i=0; i<ntypes; i++) {
    for (j=i; j<ntypes; j++) {
      ter_r_cut[i][j]  = ter_r_cut[j][i]  = ters_r_cut[n];
      ter_r2_cut[i][j] = ter_r2_cut[j][i] = ter_r_cut[i][j] * ter_r_cut[i][j];
      ++n;      
    }
  }

  /* Cutoff radius for cell decomposition */
  tmp = 0.0;
  for (i=0; i<ntypes; ++i)
    for (j=0; j<ntypes; ++j)
      tmp = MAX( tmp, ter_r2_cut[i][j] );
  r2_cut = MAX(r2_cut,tmp);

}  

#endif

/******************************************************************************
*
*  do_neighbour_tables -- Calculates neighbour tables
*
******************************************************************************/

void do_neighbour_tables(cell *p, cell *q, vektor pbc)
{
  int i, j;
  int jstart;
  int q_typ, p_typ;
  vektor d, tmp_d;
  real radius;

  /* For each atom in first cell */
  for (i=0; i<p->n; ++i) {

    /* Compute only neighbour tables for existing atoms */
    if ( p->del[i] == 0 ) {

      /* All atoms are white */ 
      p->color[i]    = -1;

      tmp_d.x = p->ort[i].x - pbc.x;
      tmp_d.y = p->ort[i].y - pbc.y;
#ifndef TWOD
      tmp_d.z = p->ort[i].z - pbc.z;
#endif
      p_typ   = p->sorte[i];

      jstart = (p==q ? i+1 : 0);
    
      /* For each atom in neighbouring cell */
      for (j = jstart; j < q->n; ++j) {

	if ( q->del[j] == 0 ) {
	  q_typ = q->sorte[j];
	  
	  /* Calculate distance  */
	  d.x = q->ort[j].x - tmp_d.x;
	  d.y = q->ort[j].y - tmp_d.y;
#ifndef TWOD
	  d.z = q->ort[j].z - tmp_d.z;
#endif

	  radius = sqrt(SPROD(d,d));

	  /* Make neighbour tables */
#ifdef TERSOFF
	  if (radius <= ter_r_cut[p_typ][q_typ])
#else
	    if (radius <= r_max)
#endif 
	    {        
	      neightab *neigh;

	      /* Update neighbour table of particle i */
	      if ( first == 0 )
		neigh = &p->neightab_array[i];
	      else 
		neigh = &p->perm_neightab_array[i];

	      if (neigh->n_max <= neigh->n ) {
		error("Neighbour table too small, increase neigh_len");
	      }

	      neigh->typ[neigh->n] = q_typ;
	      neigh->cl [neigh->n] = q;
	      neigh->num[neigh->n] = j;

	      neigh->n++;

	      /* Update neighbour table of particle j */
	      if ( first == 0 )
		neigh = &q->neightab_array[j];
	      else
		neigh = &q->perm_neightab_array[j];

	      if (neigh->n_max <= neigh->n ) {
		error("Neighbour table too small, increase neigh_len");
	      }

	      neigh->typ[neigh->n] = p_typ;
	      neigh->cl [neigh->n] = p;
	      neigh->num[neigh->n] = i;

	      neigh->n++;
	    }
	}
      } /* for j */
    }
  } /* for i */

}

/******************************************************************************
*
*  update_neighbour_tables -- Calculates new neighbour tables
*
******************************************************************************/
void update_neighbour_tables(cell *p, int i) {

  neightab *neigh, *next_neigh;
  int j, k, l;
  cell *q;

  /* Neighbour table of p, i */
  neigh = &p->neightab_array[i];

  /* For all neighbours of p, i */
  for ( k=0; k<neigh->n; k++) {
    q = neigh->cl[k];
    j = neigh->num[k];
    if ( q->del[j] == 0 ) {
      /* Neighbour tables of neighbours */
      next_neigh = &q->neightab_array[j];

      for ( l=0; l<next_neigh->n; l++ )
	if ( next_neigh->cl[l] == p && next_neigh->num[l] == i ) {
	  /* Remove atom p, i from neighbour tables */
	  next_neigh->cl[l]  = next_neigh->cl[next_neigh->n - 1];
	  next_neigh->num[l] = next_neigh->num[next_neigh->n - 1];
	  next_neigh->typ[l] = next_neigh->typ[next_neigh->n - 1];
	  next_neigh->n--;
	}
    }
  }
}

/******************************************************************************
*
*  write_data -- output ring statistics
*
******************************************************************************/
void write_data(void) {

  FILE *out;
  str255 fname;
  int i;

  if ( total_rings > 0 ) {

    /* Output on stdout */
    printf("--------------------------------------------------\n");
    printf(" Ring length  Total Occurence  Relative Occurence\n\n");
    for ( i=3; i<=max_length; i++ )
      printf("      %d            %d                %.3f %%\n", i, histogram[i]/2, (real)100*histogram[i]/(total_rings));
    printf("\n--------------------------------------------------\n");

    /* Output in outfile */
    sprintf(fname,"%s.ring",infilename);
  
    out = fopen(fname,"w");
    if (NULL == out) 
      error("Cannot open ring file.");

    fprintf(out, "#Ring length     Relative Occurence\n");
    for ( i=3; i<=max_length; i++ )
      fprintf(out, "%d\t%.6f\n", i, (real)histogram[i]/(total_rings));

    fclose(out);
  }

}

/******************************************************************************
*
*  search_rings -- searches for rings
*
******************************************************************************/

void search_rings(void) {

  cell *p, *act_cell;
  int  ic, jc, kc, intv, stars = 0, count = 0;
  int  i, act_num, j;
  neightab *pre_neigh, *act_neigh;
  int    icc, jcc, kcc, m;
  cell *q;

  /* Status line */
  intv = natoms / 50;
  printf("\nSearching rings.\n");
  printf("|--------------------------------------------------|\n|");
  fflush(stdout);

  /* Initialize histogram */
  if ( (histogram = (int *)malloc((max_length+1)*sizeof(int))) == NULL )
    error("Cannot allocate memory for histogram!\n");
  for( i=0; i<=max_length; i++)
    histogram[i] = 0;

  /* Allocate memory for stack */
  stack = (atom * ) malloc( (max_length+1) * sizeof(atom) );
  if ( stack == NULL )
    error("Cannot allocate memory for stack!");

  /* For each cell */
  for (ic=0; ic < cell_dim.x; ++ic)
    for (jc=0; jc < cell_dim.y; ++jc)
#ifndef TWOD
      for (kc=0; kc < cell_dim.z; ++kc)
#endif
      {
#ifdef TWOD
        p = PTR_2D_V(cell_array,ic,jc  ,cell_dim);
#else
        p = PTR_3D_V(cell_array,ic,jc,kc,cell_dim);
#endif

	/* For each atom in cell */
	  for ( i=0; i<p->n; ++i ) {

	    /* Initialize neighbour tables */
	    for (icc=0; icc < cell_dim.x; ++icc)
	      for (jcc=0; jcc < cell_dim.y; ++jcc)
#ifndef TWOD
		for (kcc=0; kcc < cell_dim.z; ++kcc)
#endif
		{
#ifdef TWOD
		  q = PTR_2D_V(cell_array,icc,jcc  ,cell_dim);
#else
		  q = PTR_3D_V(cell_array,icc,jcc,kcc,cell_dim);
#endif
		  for ( m=0; m<q->n; ++m )
		    if ( q->del[m] == 0 ) 
		      q->color[m] = -1;
		}

	    /* Compute distances of neighbours */

	    /* First atom has distance 0 and gets grey */
	    p->hops[i] = 0;
	    p->color[i] = 0;

	    /* Put base atom into queue (for computation of distances) */
	    queue = queue_create( p, i);
	    queue_length = 1;

	    /* Compute distances from atom p,i */
	    compute_hops();


	    /* Push base atom onto stack */	    
	    stack[0].cl  = p;
	    stack[0].num = i;
	    stack[0].status = 1;
	    stack_end = 0;

	    /* Neighbour table of base atom */
	    pre_neigh = &p->neightab_array [i];
	    
	    /* For all neighbour atoms */
	    for ( j=0; j<pre_neigh->n; j++ ) {
 
	      /* Cell and number of first neighbour atom */
	      act_cell = (cell *) pre_neigh->cl [j];
	      act_num  =          pre_neigh->num[j];

	      /* Push first neighbour atom onto stack */
	      stack[1].cl = act_cell;
	      stack[1].num = act_num;
	      stack[1].status = 1;
	      stack_end = 1;

	      /* status =  1: distance from base atom increases,
		 status =  0: distance from base atom does not change,
		 status = -1: distance from base atom decreases */

	      /* Neighbour table of first neighbour atom */
	      act_neigh = &act_cell->neightab_array [act_num];

	      /* If first neighbour has more than one neighbour, go on */
	      if ( act_neigh->n > 1 ) 
		go_forward(); 

	    }

	    /* Remove atom p,i from system (i.e., mark it)
	     and update neighbour tables */
	    p->del[i] = 1;
	    update_neighbour_tables( p, i );

	    /* Status line */
	    count++;
	    if ( count > intv ) {
	      printf("*");
	      fflush(stdout);
	      count = 0;
	      stars++;
	    }

	  }
      }

  /* Status line */
  if ( stars <= 50 ) 
    for ( i=0; i<(50-stars); i++)
      printf("*");
  printf("|\n");

  if ( total_rings > 0 )
    printf("Total number of rings found: %d\n\n",total_rings/2);
  else
    printf("No rings found.\n\n");
}

/******************************************************************************
*
*  go_forward -- next step in ring
*
******************************************************************************/

void go_forward(void) {

  int   j, k, cand, visited = 0;
  cell  *p, *act_cell, *pre_cell, *neigh_cell;
  int   i, act_num, pre_num, neigh_num;
  neightab *act_neigh, *next_neigh;
  int   delta;

  act_cell  = stack[stack_end].cl;
  act_num   = stack[stack_end].num;
  act_neigh = &act_cell->neightab_array [act_num];

  pre_cell  = stack[stack_end - 1].cl;
  pre_num   = stack[stack_end - 1].num;

  p         = stack[0].cl;
  i         = stack[0].num;

  /* Test whether neighbour atoms of actual atom
     are vertices of possible ring */
  for ( j=0; j<act_neigh->n; j++ ) {

    /* cell and number of neighbour atom */
    neigh_cell = (cell *) act_neigh->cl [j];
    neigh_num  =          act_neigh->num[j];

    /* Check whether neighbour has already been visited */
    visited = 0;
    for ( k=1; k<=stack_end; k++)
      if ( (neigh_cell == stack[k].cl) && (neigh_num == stack[k].num) )
	visited = 1;

    /* Neighbour has not been visited. Disregard rings of length 2 */
    if ( visited == 0 && (neigh_cell != p || neigh_num != i || stack_end > 1 ) ) {

      /* Is path a closed ring ? */
      if ( neigh_cell == p && neigh_num == i ) {

	/* Test whether ring is shortest path ring */
	if ( sp_ring() ) {
	  ++histogram[stack_end+1];
	  ++total_rings;
	}
      }
      /* Path is not closed */
      else {

	/* Push new vertex on stack */
	++stack_end;
	stack[stack_end].cl = neigh_cell;
	stack[stack_end].num = neigh_num;

	/* Check path for unimodularity.
	   delta = 1:  path length increases,
	   delta = 0:  path length does not change,
	   delta = -1: path length decreases */
	delta = neigh_cell->hops[neigh_num] - act_cell->hops[act_num];

        if ( delta == 1 ) {
	  if ( stack[stack_end-1].status == 1 ) { 
	    cand = 1;
	    stack[stack_end].status = 1;
	  }
	  else if ( stack[stack_end-1].status == 0 ) { 
	    cand = 0;
	    stack[stack_end].status = 0;
	  }
	  else if ( stack[stack_end-1].status == -1 ) {
	    cand = 0;
	    stack[stack_end].status = -1;
	  }
	}
	else if (delta == 0 ) {
	  if ( stack[stack_end-1].status == 1 ) {
	    cand = 1;
	    stack[stack_end].status = 0;
	  }
	  else if ( stack[stack_end-1].status == 0 ) {
	    cand = 0;
	    stack[stack_end].status = 0;
	  }
	  else if (stack[stack_end-1].status == -1 ) {
	    cand = 0;
	    stack[stack_end].status = -1;
	  }
	}
	else if ( delta == -1 ) {
	  if ( stack[stack_end-1].status == 1 ) {
	    cand = 1;
	    stack[stack_end].status = -1;
	  }
	  else if (stack[stack_end-1].status == 0 ) {
	    cand = 1;
	    stack[stack_end].status = -1;
	  }
	  else if (stack[stack_end-1].status == -1 ) {
	    cand = 1;
	    stack[stack_end].status = -1;
	  }
	}
	/* Next vertex is a candidate for a path vertex */
	if ( cand == 1 ) {
	  next_neigh = &neigh_cell->neightab_array[neigh_num];

	  /* If there are further neighbours and the distance is not greater
	     than max_length/2, iterate ring search */ 
	  if ( next_neigh->n > 1 && stack_end <= max_length && neigh_cell->hops[neigh_num] <= max_length/2 ) {
	    go_forward();
	  }
	  else {
	    /* Only one neighbour or path longer than max_length/2 */
	    stack_end--;
	  }
	} 
	else { 
	  /* Path not unimodular */
	  stack_end--;
	}
      } /* no ring found */
      
    } /* neighbour is different from previous atoms. Not ring of length 2 */  
    
  }
  stack_end--;
}
	    
/******************************************************************************
*
*  compute_hops -- compute distances from base atom, breadth-first search
*
******************************************************************************/

void compute_hops() {

  neightab *neigh;
  cell *p, *neigh_cell;
  int  i, j, neigh_num;

  while ( queue_length > 0 ) {
   
    /* Take atom from queue */
    queue_dequeue(&queue, &p, &i);
    queue_length--;

    /* Atom gets black */
    p->color[i] = 1;

    /* Neighbour table of atom */
    neigh = &p->neightab_array[i];
    
    /* For all neighbour atoms */
    for ( j=0; j<neigh->n; ++j ) {

      neigh_cell = (cell *) neigh->cl [j];
      neigh_num  = neigh->num [j];
      
      /* Neighbour atom exists and is white */ 
      if ( neigh_cell->del[neigh_num] == 0 && neigh_cell->color[neigh_num] == -1 ) {
	
	/* Put neighbour atom into queue */
	queue = queue_enqueue(queue, neigh_cell, neigh_num);
	queue_length++;

	/* Compute distance of neighbour atom and color it grey */
	neigh_cell->hops[neigh_num] = p->hops[i] + 1;
	neigh_cell->color[neigh_num] = 0;
      }
    }

  } /* while */
}

/******************************************************************************
*
*  sp_ring -- checks whether ring is shortest path ring
*
******************************************************************************/

int sp_ring(void) {

  int i, j, k, l, ic, jc, kc, lc;
  cell *p, *q, *pc, *act_cell, *neigh_cell;
  int act_num, neigh_num;
  int r_dist, even = 0;
  neightab *neigh;

  /* Check whether number of vertices in ring is even */
  if ( ( stack_end + 1 ) % 2 == 0 ) 
    even = 1;

  /* Diameter of ring */
  r_dist = ( stack_end + 1 ) / 2;

  /* For half the number of vertices in ring */
  for ( l=0; l<(even==1?r_dist:(r_dist+1)); l++ ) {

    p = stack[l].cl;
    i = stack[l].num;

    /* Compute distances of neighbours of p,i. 
       Breadth-first search, for all atoms */

    /* Initializations */
    for (ic=0; ic < cell_dim.x; ++ic)
      for (jc=0; jc < cell_dim.y; ++jc)
#ifndef TWOD
	for (kc=0; kc < cell_dim.z; ++kc)
#endif
	{
#ifdef TWOD
	  pc = PTR_2D_V(cell_array,ic,jc  ,cell_dim);
#else
	  pc = PTR_3D_V(cell_array,ic,jc,kc,cell_dim);
#endif
	  for ( lc=0; lc<pc->n; lc++ ) 
	    pc->sp_color[lc] = -1;
	}

    p->sp_hops[i]  = 0;
    p->sp_color[i] = 0;

    sp_queue = queue_create(p,i);
    sp_queue_length = 1;

    while ( sp_queue_length > 0 ) {

      queue_dequeue( &sp_queue, &act_cell, &act_num );
      sp_queue_length--;

      /* Actual atom gets black */
      act_cell->sp_color[act_num] = 1;

      /* Neighbour table of actual atom */
      neigh = &act_cell->perm_neightab_array[act_num];
    
      /* For all neighbours of actual atom */
      for ( k=0; k<neigh->n; k++ ) {

	neigh_cell = (cell *) neigh->cl [k];
	neigh_num  = neigh->num [k];
      
	/* If neighbour is white and not farther away than max_length, 
	   enqueue it */
	if ( neigh_cell->sp_color[neigh_num] == -1 
	     && act_cell->sp_hops[act_num] <= (max_length/2 + 1) ) {
	  
	  sp_queue = queue_enqueue( sp_queue, neigh_cell, neigh_num);
	  sp_queue_length++;

	  /* Distance of neighbour atom, which now gets grey */
	  neigh_cell->sp_hops[neigh_num] = act_cell->sp_hops[act_num] + 1;
	  neigh_cell->sp_color[neigh_num] = 0;
	}
      }

    } /* while */

    /* Antipodal atom */
    q = stack[l+r_dist].cl;
    j = stack[l+r_dist].num;
    
    /* Ring is not shortest path */
    if ( q->sp_hops[j] < r_dist ) {
      return 0;
    }
    else if ( even == 0 && l+r_dist+1<=stack_end ) {
      q = stack[l+r_dist+1].cl;
      j = stack[l+r_dist+1].num;
      if  ( q->sp_hops[j] < r_dist )
	return 0;
    }
    
  } /* For antipodal vertices */

  return 1;

}

/******************************************************************************
*
*  Functions for queues
*
******************************************************************************/

Queue_elmt *queue_create(cell *cl, int num) {

  Queue_elmt *newelmt;

  newelmt = (Queue_elmt *) malloc( sizeof( Queue_elmt));
  if ( newelmt == NULL )
    error("Cannot allocate memory for queue.");

  newelmt->cl   = cl;
  newelmt->num  = num;
  newelmt->next = NULL;

  return newelmt;
}

Queue_elmt *queue_enqueue(Queue_elmt *queue, cell *cl, int num) {

  Queue_elmt *newelmt;

  if ((newelmt = (Queue_elmt *)malloc(sizeof(Queue_elmt))) == NULL )
    error("Cannot allocate memory for queue.");

  newelmt->cl   = cl;
  newelmt->num  = num;
  newelmt->next = queue;

  return newelmt;

}

void queue_dequeue(Queue_elmt **qptr, cell **clptr, int *numptr) {

  Queue_elmt *q;
  
  /* Only one element in queue */
  if ( (*qptr)->next == NULL ) {

    *clptr  = (*qptr)->cl;
    *numptr = (*qptr)->num;

    free(*qptr);
    *qptr = NULL;


  }

  else {

    q = *qptr;

    while ( q->next->next!=NULL ) {

      q=q->next;

    }

    *clptr = q->next->cl;
    *numptr = q->next->num;

    free(q->next);
    q->next = NULL;

  }

  return;

}


    





