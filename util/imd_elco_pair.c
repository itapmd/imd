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
* imd_elco_pair
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#ifndef ELCO
#define ELCO
#endif

#ifndef PAIR_POT
#define PAIR_POT
#endif

#include "util.h"
#include "potaccess.h"

/******************************************************************************
*
*  init_pair -- Initializations
*
******************************************************************************/

void init_pair(void) 
{
  /* read pair potential file */
  read_pot_table(&pair_pot, potfilename, ntypes*ntypes);
}

void read_pot_table( pot_table_t *pt, char *filename, int ncols )
{
  FILE *infile;
  char buffer[1024], msg[255];
  char *res;
  int  have_header=0, have_format=0, end_header;
  int  format=2;
  int  size=ncols;
  int i;

  /* open file */
  infile = fopen(filename,"r");
  if (NULL == infile) {
    sprintf(msg,"Could not open file %s\n",filename);
    error(msg);
  }

  /* read the header */
  do {
    /* read one line */
    res=fgets(buffer,1024,infile);
    if (NULL == res) {
      sprintf(msg,"Unexpected end of file in %s",filename);
      error(msg);
    }
    
    /* see if it is a header line */
    if (buffer[0]=='#') {
      have_header = 1;
      /* stop after last header line */
      end_header = (buffer[1]=='E');
      /* see if it is the format line */
      if (buffer[1]=='F') {
	/* format complete? */
	if (2!=sscanf( (const char*)(buffer+2), "%d%d", &format, &size )) {
	  sprintf(msg,"Corrupted format header line in file %s",filename);
	  error(msg);
	}

	/* right number of columns? */
	if (size!=ncols) {
	  sprintf(msg,"Wrong number of data columns in file %s",filename);
	  error(msg);
	}
	/* recognized format? */
	if ((format!=1) && (format!=2)) {
	  sprintf(msg,"Unrecognized format specified for file %s",filename);
	  error(msg);
	}
	have_format=1;
      }
    } else if (have_header) { 
      /* header does not end properly */
      sprintf(msg,"Corrupted header in file %s",filename);
      error(msg);
    } else {
      /* we have no header, stop reading further */
      end_header=1;
    }
  } while (!end_header);
  
  /* did we have a format in the header */
  if ((have_header) && (!have_format)) {
    sprintf(msg,"Format not specified in header of file %s",filename);
    error(msg);
  }
  
  /* rewind if there was no header */
  if (!have_header) rewind(infile);
  
  /* warn if we have no header */
  if (!have_header) {
    sprintf(msg,"File %s has no header",filename);
    printf(msg);
  }
  
  /* allocate info block of function table */
  pt->maxsteps = 0;
  pt->begin    = (real *) malloc(size*sizeof(real));
  pt->end      = (real *) malloc(size*sizeof(real));
  pt->step     = (real *) malloc(size*sizeof(real));
  pt->invstep  = (real *) malloc(size*sizeof(real));
  if ((pt->begin   == NULL) || (pt->end == NULL) || (pt->step == NULL) || 
      (pt->invstep == NULL)) {
    sprintf(msg,"Cannot allocate info block for function table %s.",filename);
    error(msg);
  }

  /* catch the case where potential is identically zero */
  for (i=0; i<size; ++i) {
    pt->end[i] = 0.0;
  }

  if (format==1) read_pot_table1(pt, size, filename, infile);
  if (format==2) read_pot_table2(pt, size, filename, infile);

  fclose(infile);

}

void read_pot_table1(pot_table_t *pt, int ncols, char *filename, FILE *infile)
{
  int i, k;
  int tablesize, npot=0;
  real val, delta;
  real r2, r2_start, r2_step;
  str255 msg;
  real cellsz = 0.0;

  /* allocate the function table */
  pt->maxsteps = PSTEP;
  tablesize = ncols * pt->maxsteps;
  pt->table = (real *) malloc(tablesize*sizeof(real));
  if (NULL==pt->table) {
    sprintf(msg,"Cannot allocate memory for function table %s.",filename);
    error(msg);
  }

  /* input loop */
  while (!feof(infile)) {

    /* still some space left? */ 
    if (((npot%PSTEP) == 0) && (npot>0)) {
      pt->maxsteps += PSTEP;
      tablesize = ncols * pt->maxsteps;
      pt->table = (real *) realloc(pt->table, tablesize*sizeof(real));
      if (NULL==pt->table) {
        sprintf(msg,"Cannot extend memory for function table %s.",filename);
        error(msg);
      }
    }
    /*  read in potential */
    if ( 1 != fscanf(infile,"%lf",&r2) ) break;
    if (npot==0) r2_start = r2;  /* catch first value */
    for (i=0; i<ncols; ++i) {
      if (( 1 != fscanf(infile,"%lf", &val))) 
        error("Line incomplete in potential file.");
      *PTR_2D(pt->table,npot,i,ncols) = val;
      if (val!=0.0) pt->end[i] = r2; /* catch last non-zero value */
    }
    ++npot;
  }

  r2_step = (r2 - r2_start) / (npot-1);

  if (ncols==ntypes) {
    printf("Read tabulated function %s for %d atoms types.\n",
	   filename,ncols);
  } else {
    printf("Read tabulated function %s for %d pairs of atoms types.\n",
	   filename,ncols);
  }


  /* fill info block, and shift potential to zero */
  for (i=0; i<ncols; ++i) {
    pt->begin[i] = r2_start;
    pt->step[i] = r2_step;
    pt->invstep[i] = 1.0 / r2_step;
    delta = *PTR_2D(pt->table,(npot-1),i,ncols);
    /* do not shift embedding energy of EAM */
    if (ncols==ntypes*ntypes) {
      if (delta!=0.0) {
	printf("Potential %1d%1d shifted by %f\n",
	       (i/ntypes),(i%ntypes),delta);
        for (k=0; k<npot; ++k) *PTR_2D(pt->table,k,i,ncols) -= delta;
      } else {
        pt->end[i] += r2_step;
      }
    }
    if (ncols==ntypes*ntypes) cellsz = MAX(cellsz,pt->end[i]);
  }
  printf("\n");
  /* The interpolation uses k+1 and k+2, so we make a few copies 
     of the last value at the end of the table */
  for (k=1; k<=5; ++k) {
    /* still some space left? */ 
    if (((npot%PSTEP) == 0) && (npot>0)) {
      pt->maxsteps += PSTEP;
      tablesize = ncols * pt->maxsteps;
      pt->table = (real *) realloc(pt->table, tablesize*sizeof(real));
      if (NULL==pt->table) {
        sprintf(msg,"Cannot extend memory for function table %s.",filename);
        error(msg);
      }
    }
    for (i=0; i<ncols; ++i)
      *PTR_2D(pt->table,npot,i,ncols) 
          = *PTR_2D(pt->table,npot-1,i,ncols);
    ++npot;
  }
  r2_cut = MAX(r2_cut,cellsz);
}

void read_pot_table2(pot_table_t *pt, int ncols, char *filename, FILE *infile)
{
  int i, k, *len;
  int tablesize;
  real val, numstep;
  str255 msg;
  real cellsz = 0.0;

  len = (int  *) malloc(ncols * sizeof(real));
  if (len==NULL) error("allocation failed in read_pot_table");

  /* read the info block of the function table */
  for(i=0; i<ncols; i++) {
    if (3 != fscanf(infile, "%lf %lf %lf",
		    &pt->begin[i], &pt->end[i], &pt->step[i])) { 
      sprintf(msg, "Info line in %s corrupt.", filename);
      error(msg);
    }
    if (ncols==ntypes*ntypes) cellsz = MAX(cellsz,pt->end[i]);

    pt->invstep[i] = 1.0 / pt->step[i];
    numstep        = 1 + (pt->end[i] - pt->begin[i]) / pt->step[i];
    len[i]         = (int) (numstep+0.49);  
    pt->maxsteps   = MAX(pt->maxsteps, len[i]);

    /* some security against rounding errors */
    if ((fabs(len[i] - numstep) >= 0.1)) {
      char msg[255];
      sprintf(msg,"numstep = %f rounded to %d in file %s.",
              numstep, len[i], filename);
      printf(msg);
    }
  }

  /* allocate the function table */
  /* allow some extra values at the end for interpolation */
  tablesize = ncols * (pt->maxsteps+3);
  pt->table = (real *) malloc(tablesize * sizeof(real));
  if (NULL==pt->table) {
    sprintf(msg,"Cannot allocate memory for function table %s.",filename);
    error(msg);
  }

  /* input loop */
  for (i=0; i<ncols; i++) {
    for (k=0; k<len[i]; k++) {
      if (1 != fscanf(infile,"%lf", &val)) {
	sprintf(msg, "wrong format in file %s.", filename);
	error(msg);
      }
      *PTR_2D(pt->table,k,i,ncols) = val;
    }
    /* make some copies of the last value for interpolation */
    for (k=len[i]; k<len[i]+3; k++)
      *PTR_2D(pt->table,k,i,ncols) = val;
  }
  if (ncols==ntypes) {
    printf("Read tabulated function %s for %d atoms types.\n",
	   filename,ncols);
  } else {
    printf("Read tabulated function %s for %d pairs of atoms types.\n",
	   filename,ncols);
  }
    printf("Maximal length of table is %d.\n",pt->maxsteps);

    r2_cut = MAX(cellsz,r2_cut);
}

/******************************************************************************
*
*  do_elco_pair -- computes stress tensor and elastic constants using  
*                  tabulated pair potentials
*
******************************************************************************/

void do_elco_pair(cell *p, cell *q, vektor pbc)
{
  int i, j;
  vektor d;
  int  p_typ, q_typ, col1, inc = ntypes * ntypes, is_short = 0;
  real r2, tmp, phi, dphi, ddphi, dddphi;
#ifdef EAM
  real rho_h;
  int col2;
  real rho_i_strich, rho_i_zweistrich, rho_i_dreistrich;
  real rho_j_strich, rho_j_zweistrich, rho_j_dreistrich;
#endif

  /* For each atom in first cell */
  for (i=0; i<p->n; ++i) 
    /* For each atom in neighbouring cell */
    for (j=((p==q) ? i+1 : 0); j<q->n; ++j) {
      
      /* Calculate distance */
      d.x = q->ort[j].x - p->ort[i].x + pbc.x;
      d.y = q->ort[j].y - p->ort[i].y + pbc.y;
#ifndef TWOD
      d.z = q->ort[j].z - p->ort[i].z + pbc.z;
#endif
      r2    = SPROD(d,d);

      p_typ = p->sorte[i];
      q_typ = q->sorte[j];
      col1  = p_typ * ntypes + q_typ;

      if ( r2 <= pair_pot.end[col1] ) {

	PAIR_INT4(phi, dphi, ddphi, dddphi, pair_pot, col1, inc, r2, is_short)

	/* Compute potential energy */
	epot += phi;

	/* Compute stress and elastic constants */
	tmp = d.x * d.x * dphi; 
	p->stress[i].xx += tmp;
	q->stress[j].xx += tmp;
	sigma.xx        += 2.0 * tmp;
	tmp = d.x * d.y * dphi; 
	p->stress[i].xy += tmp;
	q->stress[j].xy += tmp;
	sigma.xy        += 2.0 * tmp;
	tmp = d.y * d.y * dphi; 
	p->stress[i].yy += tmp;
	q->stress[j].yy += tmp;
	sigma.yy        += 2.0 * tmp;
#ifndef TWOD
	tmp = d.y * d.z * dphi; 
	p->stress[i].yz += tmp;
	q->stress[j].yz += tmp;
	sigma.yz        += 2.0 * tmp;
	tmp = d.z * d.z * dphi; 
	p->stress[i].zz += tmp;
	q->stress[j].zz += tmp;
	sigma.zz        += 2.0 * tmp;
	tmp = d.z * d.x * dphi; 
	p->stress[i].zx += tmp;
	q->stress[j].zx += tmp;
	sigma.zx        += 2.0 * tmp;
#endif

	tmp = 2.0 * d.x * d.x * d.x * d.x * ddphi + d.x * d.x * dphi;
	p->elco[i].c11 += tmp;
	q->elco[j].c11 += tmp;
	c.c11          += 2.0 * tmp;
	tmp = 2.0 * d.x * d.x * d.y * d.y * ddphi;
	p->elco[i].c12 += tmp;
	q->elco[j].c12 += tmp;
	c.c12          += 2.0 * tmp;
	tmp += 0.5 * ( d.x * d.x + d.y * d.y ) * dphi; 
	p->elco[i].c66 += tmp;
	q->elco[j].c66 += tmp;
	c.c66          += 2.0 * tmp;
	tmp = 2.0 * d.y * d.y * d.y * d.y * ddphi + d.y * d.y * dphi;
	p->elco[i].c22 += tmp;
	q->elco[j].c22 += tmp;
	c.c22          += 2.0 * tmp;
#ifndef TWOD
	tmp = 2.0 * d.x * d.x * d.z * d.z * ddphi;
	p->elco[i].c13 += tmp;
	q->elco[j].c13 += tmp;
	c.c13          += 2.0 * tmp;
	tmp += 0.5 * ( d.x * d.x + d.z * d.z ) * dphi; 
	p->elco[i].c55 += tmp;
	q->elco[j].c55 += tmp;
	c.c55          += 2.0 * tmp;
	tmp = 2.0 * d.y * d.y * d.z * d.z * ddphi;
	p->elco[i].c23 += tmp;
	q->elco[j].c23 += tmp;
	c.c23          += 2.0 * tmp;
	tmp += 0.5 * ( d.y * d.y + d.z * d.z ) * dphi; 
	p->elco[i].c44 += tmp;
	q->elco[j].c44 += tmp;
	c.c44          += 2.0 * tmp;
	tmp = 2.0 * d.z * d.z * d.z * d.z * ddphi + d.z * d.z * dphi;
	p->elco[i].c33 += tmp;
	q->elco[j].c33 += tmp;
	c.c33          += 2.0 * tmp;
	tmp = 2.0 * d.x * d.y * d.z * d.z * ddphi + 0.25 * d.x * d.y * dphi;
	p->elco[i].c45 += tmp;
	q->elco[j].c45 += tmp;
	c.c45          += 2.0 * tmp;
	tmp = 2.0 * d.x * d.y * d.y * d.z * ddphi + 0.25 * d.x * d.z * dphi;
	p->elco[i].c46 += tmp;
	q->elco[j].c46 += tmp;
	c.c46          += 2.0 * tmp;
	tmp -= 0.25 * d.x * d.z * dphi;
	p->elco[i].c25 += tmp;
	q->elco[j].c25 += tmp;
	c.c25          += 2.0 * tmp;
	tmp = 2.0 * d.x * d.x * d.y * d.z * ddphi + 0.25 * d.y * d.z * dphi;
	p->elco[i].c56 += tmp;
	q->elco[j].c56 += tmp;
	c.c56          += 2.0 * tmp;
	tmp -= 0.25 * d.y * d.z * dphi;
	p->elco[i].c14 += tmp;
	q->elco[j].c14 += tmp;
	c.c14          += 2.0 * tmp;
	tmp = 2.0 * d.x * d.x * d.x * d.z * ddphi + 0.5 * d.x * d.z * dphi;
	p->elco[i].c15 += tmp;
	q->elco[j].c15 += tmp;
	c.c15          += 2.0 * tmp;
#endif
	tmp = 2.0 * d.x * d.x * d.x * d.y * ddphi + 0.5 * d.x * d.y * dphi;
	p->elco[i].c16 += tmp;
	q->elco[j].c16 += tmp;
	c.c16          += 2.0 * tmp;
	tmp = 2.0 * d.x * d.y * d.y * d.y * ddphi + 0.5 * d.x * d.y * dphi;
	p->elco[i].c26 += tmp;
	q->elco[j].c26 += tmp;
	c.c26          += 2.0 * tmp;
#ifndef TWOD
	tmp = 2.0 * d.y * d.y * d.y * d.z * ddphi + 0.5 * d.y * d.z * dphi;
	p->elco[i].c24 += tmp;
	q->elco[j].c24 += tmp;
	c.c24          += 2.0 * tmp;
	tmp = 2.0 * d.y * d.z * d.z * d.z * ddphi + 0.5 * d.y * d.z * dphi;
	p->elco[i].c34 += tmp;
	q->elco[j].c34 += tmp;
	c.c34          += 2.0 * tmp;
	tmp = 2.0 * d.x * d.z * d.z * d.z * ddphi + 0.5 * d.x * d.z * dphi;
	p->elco[i].c35 += tmp;
	q->elco[j].c35 += tmp;
	c.c35          += 2.0 * tmp;
	tmp = 2.0 * d.x * d.y * d.z * d.z * ddphi;
	p->elco[i].c36 += tmp;
	q->elco[j].c36 += tmp;
	c.c36          += 2.0 * tmp;
#endif
	press          +=   2.0 * dphi * r2;

	bulkm          += ( 2.0 * ddphi * r2 + dphi ) * 2.0 * r2;

	dbulkm_dp      += ( 2.0 * dddphi * r2 + 3.0 * ddphi ) * 4.0 * r2 * r2;

      }

#ifdef EAM
      col2  = q_typ * ntypes + p_typ;

      /* compute host electron density */
      if ( r2 < rho_h_tab.end[col1] )  {
        VAL_FUNC(rho_h, rho_h_tab, col1,  inc, r2, is_short);
        EAM_RHO(p,i) += rho_h; 
      }
      if ( p_typ == q_typ ) {
        if ( r2 < rho_h_tab.end[col1] ) 
	  EAM_RHO(q,j) += rho_h; 
      } else {
        if ( r2 < rho_h_tab.end[col2] ) {
          VAL_FUNC(rho_h, rho_h_tab, col2, inc, r2, is_short);
          EAM_RHO(q,j) += rho_h; 
        }
      }

      /* Compute stress for EAM potential */
      if ( (r2 < rho_h_tab.end[col1]) || (r2 < rho_h_tab.end[col2]) ) {

        DERIV_FUNC(rho_i_strich, rho_i_zweistrich, rho_i_dreistrich, 
		   rho_h_tab, col2, inc, r2, is_short);

	q->eam_stress[j].xx += 2.0 * rho_i_strich * d.x * d.x;
	q->eam_stress[j].xy += 2.0 * rho_i_strich * d.x * d.y;
	q->eam_stress[j].yy += 2.0 * rho_i_strich * d.y * d.y;
#ifndef TWOD
	q->eam_stress[j].yz += 2.0 * rho_i_strich * d.y * d.z;
	q->eam_stress[j].zz += 2.0 * rho_i_strich * d.z * d.z;
	q->eam_stress[j].zx += 2.0 * rho_i_strich * d.z * d.x;
#endif
	q->eam_press[j]     +=   2.0 * rho_i_strich * r2; 
	q->eam_bulkm[j]     += ( 2.0 * rho_i_zweistrich * r2 + rho_i_strich )
	                         * 2.0 * r2; 
	q->eam_dbulkm[j]    += ( 2.0 * rho_i_dreistrich * r2 
				+ 3.0 * rho_i_zweistrich ) 
	                         * 4.0 * r2 * r2;

        if ( col1 == col2 ) {
          rho_j_strich     = rho_i_strich;
          rho_j_zweistrich = rho_i_zweistrich;
	  rho_j_dreistrich = rho_i_dreistrich;
	} 
	else {

          DERIV_FUNC(rho_j_strich, rho_j_zweistrich, rho_j_dreistrich, 
		     rho_h_tab, col1, inc, r2, is_short);
	}

	p->eam_stress[i].xx += 2.0 * rho_j_strich * d.x * d.x;
	p->eam_stress[i].xy += 2.0 * rho_j_strich * d.x * d.y;
	p->eam_stress[i].yy += 2.0 * rho_j_strich * d.y * d.y;
#ifndef TWOD
	p->eam_stress[i].yz += 2.0 * rho_j_strich * d.y * d.z;
	p->eam_stress[i].zz += 2.0 * rho_j_strich * d.z * d.z;
	p->eam_stress[i].zx += 2.0 * rho_j_strich * d.z * d.x;
#endif
	p->eam_press[i]     +=   2.0 * rho_j_strich * r2; 
	p->eam_bulkm[i]     += ( 2.0 * rho_j_zweistrich * r2 + rho_j_strich )
	                         * 2.0 * r2; 
	p->eam_dbulkm[i]    += ( 2.0 * rho_j_dreistrich * r2 
				+ 3.0 * rho_j_zweistrich )
	                         * 4.0 * r2 * r2;
      }
#endif /* EAM */
    }
}
