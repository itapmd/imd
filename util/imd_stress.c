
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
* imd_stress -- calculates stress field 
*         
* Voronoi analysis according to Allen, Tildesley
* Parameter reading and cell decomposition from imd_pair
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#ifndef STRESS
#include "util.h"

#else
#define MAIN

#include "util.h"

/******************************************************************************
*
*  Usage -- educate users
*
******************************************************************************/

void usage(void)

{ 
  printf("%s [-r<nnn>] [-e<nnn>] [-v] [-p paramter-file]\n", progname); 
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

  /* Initialize cell data structures */
  init_cells();

  /* Read atoms */
  read_stress(infilename);

  /* Calculate volume of Voronoi cells */
  voronoi();

  /* Output results */
  write_stress(restart);

  exit(0);

}


/******************************************************************************
*
*  read_stress -- reads coordinates and stress tensor into the cell-array
*
*  The file format is flat ascii, one atom per line, lines beginning
*  with '#' denote comments. Each line consists of
*
*  x y [z] s_xx s_yy [s_zz] [s_yz] [s_zx] s_xy [rest]
*
*  where
*
*  x,y,z    are the atom's coordinates
*  s_ij     are the components of the stress tensor
*  rest     is ignored until end of line
*
******************************************************************************/

void read_stress(str255 infilename)

{
  cell *to;
  FILE *infile;
  char buf[512];
  int p,nr,typ;
  double mass;
  vektor pos, sigma, sigma_offdia;
  ivektor cellc;

  infile = fopen(infilename,"r");

  if (NULL==infile) {
    sprintf(error_msg,"Cannot open atoms file %s",infilename);
    error(error_msg);
  }

  natoms=0;
  
  /* Read the input file line by line */
  while(!feof(infile)) {

    buf[0] = (char) NULL;
    fgets(buf,sizeof(buf),infile);
    while ('#'==buf[1]) fgets(buf,sizeof(buf),infile); /* eat comments */

#ifdef TWOD
    p = sscanf(buf,"%lf %lf %lf %lf %lf",
               &pos.x,&pos.y,&sigma.x,&sigma.y,&sigma_offdia.x);
#else
    p = sscanf(buf,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &nr,&typ,&mass,&pos.x,&pos.y,&pos.z,&sigma.x,&sigma.y,&sigma.z,
               &sigma_offdia.x,&sigma_offdia.y,&sigma_offdia.z);
#endif

    if (p>0) {
      /* compute target cell */
      cellc = cell_coord(pos);
      to = PTR_VV(cell_array,cellc,cell_dim);
      /* enlarge it if necessary */
      if (to->n >= to->n_max) alloc_cell(to,to->n_max+CSTEP);
      /* put the data */
      to->nummer[to->n] = nr;
      to->sorte[to->n] = typ;
      to->masse[to->n] = mass;
      to->ort[to->n] = pos;
      to->stress[to->n] = sigma;
      to->stress_offdia[to->n] = sigma_offdia;
      to->n++;
      natoms++;
    }
  }
  fclose(infile);  

}


/******************************************************************************
*
*  write_stress -- writes stress to *.stress file
*
******************************************************************************/

void write_stress(int restart)

{
  FILE *out;
  str255 fname;
  int i,j,k,l;
  cell *p; 
  real tmp;

  sprintf(fname,"%s.%u.stress",outfilename, restart);
  

  out = fopen(fname,"w");
  if (NULL == out) error("Cannot open stress file.");

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
	for (l=0; l<p->n; ++l)
	  {
	    if ( p->vol[l] > 0.0 )
	      {
		tmp = 1.0 / p->vol[l];
#ifdef TWOD
		fprintf(out, "%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n", p->ort[l].x, p->ort[l].y, p->stress[l].x*tmp, p->stress[l].y*tmp, p->stress_offdia[l].x*tmp, p->vol[l]);
#else
		fprintf(out, "%d %d %f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n",p->nummer[l],p->sorte[l],p->masse[l], p->ort[l].x, p->ort[l].y, p->ort[l].z, p->stress[l].x*tmp, p->stress[l].y*tmp, p->stress[l].z*tmp, p->stress_offdia[l].x*tmp, p->stress_offdia[l].y*tmp, p->stress_offdia[l].z*tmp, p->vol[l]);
#endif
	      }
	  }
	}

  fclose(out);

  /* Write Statistics */
  if( atomcount > 0 )
    {
      printf("\n");
      printf("Maximal number of neighbour atoms:           %d\n", maxneigh );
      printf("Average number of neighbour atoms:           %.2f\n\n", (real)(sumneigh)/atomcount );
      printf("Maximal number of vertices of Voronoi cells: %d\n", maxvert );
      printf("Average number of vertices:                  %.2f\n\n", (real)(sumvert)/atomcount );
      printf("Maximal number of edges of Voronoi cells:    %d\n", maxedges );
      printf("Average number of edges:                     %.2f\n\n", (real)(sumedges)/atomcount );
#ifndef TWOD
      printf("Maximal number of faces of Voronoi cells:    %d\n", maxfaces );
      printf("Average number of faces:                     %.2f\n\n", (real)(sumfaces)/atomcount );
#endif

      printf("Total number of atoms:   %d\n", natoms);
      printf("Number of omitted atoms: %d (%.2f %%) \n\n", natoms - atomcount, (real)(natoms-atomcount)/natoms*100.0 );
    }
  else
    printf("No Voronoi cell found.\n");
}

#endif /* STRESS */

/******************************************************************************
*
*  voronoi -- calculates neighbouring points for Voronoi construction
*             and calls calculation of volume/area
*
******************************************************************************/

void voronoi(void)

{
  cell *p, *q;
  int i, j, k, l, m, n, r, s, t;
  int v, w, neighcount;
  vektor tmp;
  vektor pbc;
  real tmpdist2;
  int num = NUM;

  candcoord = (vektor *) malloc( num * sizeof(vektor));
  canddist2 = (real   *) malloc( num * sizeof(real));

  /* for each cell */
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

	/* for each atom in first cell */
	for (v=0; v<p->n; ++v)
	  {
	    neighcount = 0;

	    /* For each neighbour of this cell */
	    for (l=-1; l <= 1; ++l)
	      for (m=-1; m <= 1; ++m)
#ifndef TWOD
		for (n=-1; n <= 1; ++n)
#endif
		  {
		    /* Calculate Indicies of Neighbour */
		    r = i+l;  pbc.x = 0;
		    s = j+m;  pbc.y = 0;
#ifndef TWOD
		    t = k+n;  pbc.z = 0;
#endif

		    /* deal with periodic boundary conditions if necessary */
		    if (r<0) {
		      if (pbc_dirs.x==1) {
			r = cell_dim.x-1; 
			pbc.x -= box_x.x;      
			pbc.y -= box_x.y;
#ifndef TWOD
			pbc.z -= box_x.z;
#endif
		      } else continue;
		    }
		    if (s<0) {
		      if (pbc_dirs.y==1) {
			s = cell_dim.y-1;
			pbc.x -= box_y.x;      
			pbc.y -= box_y.y;
#ifndef TWOD
			pbc.z -= box_y.z;
#endif
		      } else continue;
		    }
#ifndef TWOD
		    if (t<0) {
		      if (pbc_dirs.z==1) {
			t = cell_dim.z-1;
			pbc.x -= box_z.x;      
			pbc.y -= box_z.y;
			pbc.z -= box_z.z;
		      } else continue;
		    }
#endif
		    if (r>cell_dim.x-1) {
		      if (pbc_dirs.x==1) {
			r = 0; 
			pbc.x += box_x.x;      
			pbc.y += box_x.y;
#ifndef TWOD
			pbc.z += box_x.z;
#endif
		      } else continue;
		    }
		    if (s>cell_dim.y-1) {
		      if (pbc_dirs.y==1) {
			s = 0; 
			pbc.x += box_y.x;      
			pbc.y += box_y.y;
#ifndef TWOD
			pbc.z += box_y.z;
#endif
		      } else continue;
		    }
#ifndef TWOD
		    if (t>cell_dim.z-1) {
		      if (pbc_dirs.z==1) {
			t = 0; 
			pbc.x += box_z.x;      
			pbc.y += box_z.y;
			pbc.z += box_z.z;
		      } else continue;
		    }
#endif
	    
		    /* Neighbour cell (note that p==q ist possible) */
#ifdef TWOD
		    q = PTR_2D_V(cell_array,r,s,cell_dim);
#else
		    q = PTR_3D_V(cell_array,r,s,t,cell_dim);
#endif
		    /* for each particle in second cell */
		    for (w=0; w<q->n; ++w)
		      {
			tmp.x    = q->ort[w].x - p->ort[v].x + pbc.x;
			tmp.y    = q->ort[w].y - p->ort[v].y + pbc.y;
#ifndef TWOD
			tmp.z    = q->ort[w].z - p->ort[v].z + pbc.z;
#endif
			tmpdist2 = SPROD( tmp, tmp );
		       
			/* Candidates for Voronoi cells */
#ifdef STRESS
			if( (tmpdist2 <= r2_cut) && (tmpdist2 > TOL2))
#else
			if( (tmpdist2 <= SQR(r_max)) && (tmpdist2 > TOL2))
#endif
			  {
			    if( neighcount > num-1 )
			      {
				num += 10;
				candcoord = (vektor *) realloc(candcoord, num * sizeof(vektor));
				canddist2 = (real   *) realloc(canddist2, num * sizeof(real));
			      }

			    candcoord[neighcount].x    = tmp.x;
			    candcoord[neighcount].y    = tmp.y;
#ifndef TWOD
			    candcoord[neighcount].z    = tmp.z;
#endif
			    canddist2[neighcount]      = tmpdist2;

			    ++ neighcount;

			  }

		      } /* w */

		  } /* lmn */

	    /* If there are less than four (three) points, a polyhedron (polygon) cannot
	       be constructed */
#ifndef TWOD
	    if( (neighnum = neighcount) < 4 )  volume = 0.0;
#else
	    if( (neighnum = neighcount) < 3 )  area   = 0.0;
#endif
	    else
	      {
		/* Sort candidates in ascending order of distance */
		sort();

		/* Perform Voronoi analysis */
#ifdef TWOD
		do_voronoi_2d();
#else
		do_voronoi_3d(); 
#endif
	      }
 
#ifndef TWOD
	    p->vol[v] = volume;
#else
	    p->vol[v] = area;
#endif

	  } /* v */

      } /* ijk */

}


/******************************************************************************
*
*  sort -- Sorts candidates for neighbour atoms in increasing order of distance
*
******************************************************************************/

void sort(void)

{
  int i, j;
  vektor tmp;
  real tmpdist2;

  for (i=(neighnum-1); i>0; --i)
    for (j=0; j<i; ++j)
      if( canddist2[j] > canddist2[j+1] )
	{
	  tmp.x            = candcoord[j].x;
	  tmp.y            = candcoord[j].y;
#ifndef TWOD
	  tmp.z            = candcoord[j].z;
#endif
	  tmpdist2         = canddist2[j];

	  candcoord[j].x   = candcoord[j+1].x;
	  candcoord[j].y   = candcoord[j+1].y;
#ifndef TWOD
	  candcoord[j].z   = candcoord[j+1].z;
#endif
	  canddist2[j]     = canddist2[j+1];

	  candcoord[j+1].x = tmp.x;
	  candcoord[j+1].y = tmp.y;
#ifndef TWOD
	  candcoord[j+1].z = tmp.z;
#endif
	  canddist2[j+1]   = tmpdist2;
	}

}


#ifndef TWOD

/******************************************************************************
*
*  do_voronoi_3d -- Calculates Voronoi cells and volume, 3d version
*
******************************************************************************/

void do_voronoi_3d(void)

{
  int       i, j, k, l, n;
  real      ab, bc, ca, da, db, dc, det, detinv, tmp;
  vektor    icoord, jcoord, kcoord;
  real      idist2, jdist2, kdist2;
  vektor    tmpvek, tmpvertex, vertex[NUM];
  int       ok, vertexcount, vertexnum, facesnum, edgesnum;
  int       *vertexnumi;
  real      area_i, height;
  real      sin, cos, maxcos;
  int       mink, ord[NUM], index[NUM], surfind[NUM];
  vektor    *coord, center;
  vektorstr *vertexloc;


  /* Allocate memory for vertices */
  vertexnumi  = (int *) malloc( neighnum * sizeof(int));
  coord       = (vektor *) malloc( neighnum * sizeof(vektor));
  vertexloc   = (vektorstr *) malloc( neighnum * sizeof(vektorstr));

  if( vertexloc == NULL || coord == NULL || vertexnumi == NULL ) 
    error("Cannot allocate memory for vertices!\n");

  vertexcount = 0;
  volume      = 0.0;

  /* Each possible vertex defined by the intersection of 3 planes is examined */
  for (i=0; i<neighnum-2; ++i)
    {
      icoord.x = candcoord[i].x;
      icoord.y = candcoord[i].y;
      icoord.z = candcoord[i].z;
      idist2   = -canddist2[i];

      for (j=i+1; j<neighnum-1; ++j)
	{
	  jcoord.x = candcoord[j].x;
	  jcoord.y = candcoord[j].y;
	  jcoord.z = candcoord[j].z;
	  jdist2   = -canddist2[j];

	  ab = icoord.x * jcoord.y - jcoord.x * icoord.y;
	  bc = icoord.y * jcoord.z - jcoord.y * icoord.z;
	  ca = icoord.z * jcoord.x - jcoord.z * icoord.x;
	  da = idist2   * jcoord.x - jdist2   * icoord.x;
	  db = idist2   * jcoord.y - jdist2   * icoord.y;
	  dc = idist2   * jcoord.z - jdist2   * icoord.z;

	  for (k=j+1; k<neighnum; ++k)
	    {
	      kcoord.x = candcoord[k].x;
	      kcoord.y = candcoord[k].y;
	      kcoord.z = candcoord[k].z;
	      kdist2   = -canddist2[k];

	      det = kcoord.x * bc + kcoord.y * ca + kcoord.z * ab;

	      /* Check whether planes intersect */
	      if( SQR(det) > TOL2 )
		{
		  detinv = 1.0 / det;

		  tmpvertex.x = ( -kdist2 * bc + kcoord.y * dc - kcoord.z * db ) * detinv;
		  tmpvertex.y = ( -kcoord.x * dc - kdist2 * ca + kcoord.z * da ) * detinv;
		  tmpvertex.z = (  kcoord.x * db - kcoord.y * da - kdist2 * ab ) * detinv;

		  /* Check whether vertex belongs to the Voronoi cell */
		  l  = 0;
		  ok = 1;

		  do {
		    if( l!=i && l!=j && l!=k)
		      ok = ( SPROD( candcoord[l], tmpvertex ) <= canddist2[l] + TOL_VERT2 )    ;
	      
		    ++l;

		  } while( ok && (l<neighnum));
		    
		    if( ok )
		      {
			vertex[vertexcount].x  = 0.5 * tmpvertex.x;
			vertex[vertexcount].y  = 0.5 * tmpvertex.y;
			vertex[vertexcount].z  = 0.5 * tmpvertex.z;
	
			++vertexcount;
		      }

		} /* Planes intersect */

	    } /* k */
	} /* j */
    } /* i */

  vertexnum = vertexcount;

  /* Check whether some vertices coincide */
  for ( i=0; i<vertexnum; ++i )
    {
      index[i] = 1;
      for ( j=i+1; j<vertexnum; ++j )
	{
	  tmpvek.x = vertex[j].x - vertex[i].x;
	  tmpvek.y = vertex[j].y - vertex[i].y;
	  tmpvek.z = vertex[j].z - vertex[i].z;

	  if ( SPROD( tmpvek, tmpvek) < TOL2 )
	    index[i] = 0;
	}  
    }

  /* Remove coincident vertices */
  j = 0;
  for ( i=0; i<vertexnum; ++i )
    if ( index[i] != 0 )
      {
	vertex[j].x  = vertex[i].x;
	vertex[j].y  = vertex[i].y;
	vertex[j].z  = vertex[i].z;

	++j;
      }

  vertexnum = j;

  /* Number of vertices of Voronoi cell must be greater than 3 */
  if(vertexnum > 3 )  
    {
      /* Check whether faces exist */
      facesnum = 0;

      /* Each neighbour atom i corresponds to at most one surface * 
       * Sum over all surfaces */
      for (i=0; i<neighnum; ++i)
	{
	  /* Coordinates of center of surface i */
	  coord[i].x = 0.5 * candcoord[i].x;
	  coord[i].y = 0.5 * candcoord[i].y;
	  coord[i].z = 0.5 * candcoord[i].z;

          /* Look for vertices that belong to surface i */
          vertexnumi[i] = 0;
	  for (j=0; j<vertexnum; ++j)
	    {
	      surfind[j] = 0;

	      vertexloc[i][j].x = vertex[j].x - coord[i].x;
	      vertexloc[i][j].y = vertex[j].y - coord[i].y;
	      vertexloc[i][j].z = vertex[j].z - coord[i].z;

	      tmp = SPROD(coord[i],vertexloc[i][j]);

	      if( SQR(tmp) < TOL_DIST2 )
		{
		  /* vertex j belongs to surface i */
		  surfind[j] = 1; 
		  ++vertexnumi[i];
		}
	    }

	  /* Surface i exists */
	  if (vertexnumi[i] > 2)
	    {
	      ++facesnum;

	      /* Compute coordinates of vertices belonging to surface i */
	      k = 0;
	      for (j=0; j<vertexnum; ++j)
		if( surfind[j] == 1)
		  {
		    vertexloc[i][k].x = vertexloc[i][j].x;
		    vertexloc[i][k].y = vertexloc[i][j].y;
		    vertexloc[i][k].z = vertexloc[i][j].z;

		    ++k;
		  }
	    }
	  /* Transform into center of mass system */
	  center.x = 0.0;
	  center.y = 0.0; 
	  center.z = 0.0;
	  
	  if( vertexnumi[i] > 2)
	    {
	      for ( j=0; j<vertexnumi[i]; ++j)
		{
		  center.x += vertexloc[i][j].x;
		  center.y += vertexloc[i][j].y;
		  center.z += vertexloc[i][j].z;
		}
	      
	      tmp       = 1.0 / vertexnumi[i];
	      center.x *= tmp;
	      center.y *= tmp;
	      center.z *= tmp;
	      
	      for ( j=0; j<vertexnumi[i]; ++j)
		{
		  vertexloc[i][j].x -= center.x;
		  vertexloc[i][j].y -= center.y;
		  vertexloc[i][j].z -= center.z;
		}

	    }

	} /* i */

      /* Number of edges of Voronoi cell */
      edgesnum = 0;

      for ( n=0; n<neighnum; ++n)
	if( vertexnumi[n] > 2)
	  edgesnum += vertexnumi[n];
      
      edgesnum /= 2;

      /* Check whether Euler relation holds */
      if ( (vertexnum - edgesnum + facesnum) == 2 )
	{      

	  /* Statistics */
	  if( neighnum  > maxneigh  ) maxneigh = neighnum;
	  sumneigh += neighnum;
	  if( vertexnum > maxvert  ) maxvert  = vertexnum;
	  sumvert  += vertexnum;
	  if( edgesnum  > maxedges ) maxedges = edgesnum;
	  sumedges += edgesnum;
	  if( facesnum  > maxfaces ) maxfaces = facesnum;
	  sumfaces += facesnum;

	  ++ atomcount;

	  /* Compute volume of Voronoi cell */

	  /* For all potential faces */
	  for (i=0; i<neighnum; ++i)
	    /* Surface i exists */
	    if(vertexnumi[i] > 2)
	      {
		/* Sort vertices of face i */
		ord[0] = 0;
		for (j=0; j<vertexnumi[i]-1; ++j)
		  {
		    maxcos = -1.0;
		    for (k=0; k<vertexnumi[i]; ++k)
		      {
			tmpvek.x = vertexloc[i][k].y * vertexloc[i][ord[j]].z - vertexloc[i][k].z * vertexloc[i][ord[j]].y;
			tmpvek.y = vertexloc[i][k].z * vertexloc[i][ord[j]].x - vertexloc[i][k].x * vertexloc[i][ord[j]].z;
			tmpvek.z = vertexloc[i][k].x * vertexloc[i][ord[j]].y - vertexloc[i][k].y * vertexloc[i][ord[j]].x; 
 
			sin = SPROD( tmpvek, coord[i]);
		      
			if( sin > TOL )
			  {
			    cos = SPROD( vertexloc[i][k], vertexloc[i][ord[j]] )/ sqrt(SPROD(vertexloc[i][k],vertexloc[i][k]))/ sqrt(SPROD(vertexloc[i][ord[j]],vertexloc[i][ord[j]]));
			    if( cos > maxcos )
			      {
				maxcos = cos;
				mink   = k;
			      }
			  }
		      }

		    ord[j+1] = mink;
		  
		  }

		/* Compute area of surface i */
		area_i = 0.0;
		height = sqrt(SPROD( coord[i], coord[i] ));
		tmp    = 1.0 / height;

		for (j=0; j<vertexnumi[i]-1; ++j)
		    {
		      tmpvek.x = vertexloc[i][ord[j+1]].y * vertexloc[i][ord[j]].z - vertexloc[i][ord[j+1]].z * vertexloc[i][ord[j]].y;
		      tmpvek.y = vertexloc[i][ord[j+1]].z * vertexloc[i][ord[j]].x - vertexloc[i][ord[j+1]].x * vertexloc[i][ord[j]].z;
		      tmpvek.z = vertexloc[i][ord[j+1]].x * vertexloc[i][ord[j]].y - vertexloc[i][ord[j+1]].y * vertexloc[i][ord[j]].x; 
			      
		      area_i += 0.5 * SPROD( tmpvek, coord[i] ) * tmp;
	
		    }
		tmpvek.x = vertexloc[i][ord[0]].y * vertexloc[i][ord[vertexnumi[i]-1]].z - vertexloc[i][ord[0]].z * vertexloc[i][ord[vertexnumi[i]-1]].y;
		tmpvek.y = vertexloc[i][ord[0]].z * vertexloc[i][ord[vertexnumi[i]-1]].x - vertexloc[i][ord[0]].x * vertexloc[i][ord[vertexnumi[i]-1]].z;
		tmpvek.z = vertexloc[i][ord[0]].x * vertexloc[i][ord[vertexnumi[i]-1]].y - vertexloc[i][ord[0]].y * vertexloc[i][ord[vertexnumi[i]-1]].x; 
		  
		area_i +=  0.5 * SPROD( tmpvek, coord[i] ) * tmp;
	
		/* Volume of Voronoi cell */	  
		volume += area_i * height / 3.0;

	      } /* vertexnum[i] > 2 */

	} /* Euler relation holds */

    } /* Number of vertices > 3 */


  free(vertexloc);
  free(coord);
  free(vertexnumi);

}


#else

/******************************************************************************
*
*  do_voronoi_2d -- Calculates Voronoi cells and volume, 2d version
*
******************************************************************************/

void do_voronoi_2d(void)

{
  int     i, j, l, n;
  real    det, detinv;
  vektor  icoord, jcoord, tmpvertex, vertex[NUM];
  int     ok, index[NUM];
  int     vertexcount, vertexnum, edgesnum;
  int     *edges;
  ivektor ivertex[NUM];
  real    idist2, jdist2;
  real    sin, cos, maxcos;
  int     minj, ord[NUM];

  /* Allocate memory for data of edges */
  edges = (int *) malloc( neighnum * sizeof(int) );

  if( edges == NULL )
    error("Cannot allocate memory for vertices");

  vertexcount = 0;
  area        = 0.0;

  /* Each possible vertex defined by the intersection of 2 lines is examined */  
  for (i=0; i<(neighnum-1); ++i)
    {
      icoord.x = candcoord[i].x;
      icoord.y = candcoord[i].y;
      idist2   = -canddist2[i];

      for (j=i+1; j<neighnum; ++j)
	{
	  jcoord.x = candcoord[j].x;
	  jcoord.y = candcoord[j].y;
	  jdist2   = -canddist2[j];

	  det = icoord.x * jcoord.y - icoord.y * jcoord.x;

	  /* check whether edges intersect */
	  if ( SQR(det) > TOL2)
	    {
	      detinv = 1.0 / det;

	      tmpvertex.x = ( icoord.y * jdist2 - jcoord.y * idist2 ) * detinv;
	      tmpvertex.y = ( jcoord.x * idist2 - icoord.x * jdist2 ) * detinv;

	      /* Check whether vertex belongs to voronoi cell */
	      l = 0;
	      ok = 1;

	      do {
		if( l!=i && l!=j )
		  ok = ( SPROD(candcoord[l] , tmpvertex ) <= canddist2[l] );

		++l;

	      } while ( (ok==1) && (l<neighnum) );

	      if( ok==1 )
		{ 
		  ivertex[vertexcount].x = i;
		  ivertex[vertexcount].y = j;
  
		  vertex[vertexcount].x  = 0.5 * tmpvertex.x;
		  vertex[vertexcount].y  = 0.5 * tmpvertex.y;

		  ++vertexcount;
		}
	    }
	}
    }

  vertexnum = vertexcount;

  /* Check whether some vertices coincide */
  for ( i=0; i<vertexnum; ++i )
    {
      index[i] = 1;
      for (j=i+1; j<vertexnum; ++j )
	if ( (SQR(vertex[j].x-vertex[i].x)<TOL2) && (SQR(vertex[j].y-vertex[i].y)<TOL2) )
	  index[i] = 0;
    }

  /* Remove coincident vertices */
  j = 0;
  for ( i=0; i<vertexnum; ++i )
    if ( index[i] != 0 )
      {
	ivertex[j].x = ivertex[i].x;
	ivertex[j].y = ivertex[i].y;

	vertex[j].x  = vertex[i].x;
	vertex[j].y  = vertex[i].y;
	
        ++j;
      }

  vertexnum = j;    

  /* Number of vertices of Voronoi cell must be greater than 2 */
  if ( vertexnum < 3 )  area = 0.0;
  else
    {
      /* Initialization */
      for (n=0; n<neighnum; ++n)
	edges[n] = 0;

      for (n=0; n<vertexnum; ++n)
	{
	  ++edges[ivertex[n].x];
	  ++edges[ivertex[n].y];
	}
     
      /* Number of edges of Voronoi cell */
      edgesnum = 0;
      for (n=0; n<neighnum; ++n)
	{
	  edgesnum += edges[n];
	}

      edgesnum /= 2;

      /* Check whether number of vertices equals number of edges */
      if ( edgesnum == vertexnum )
	{
	  
	  /* Statistics */
	  if( neighnum  > maxneigh ) maxneigh = neighnum;
	  sumneigh += neighnum;
	  if( vertexnum > maxvert  ) maxvert  = vertexnum;
	  sumvert  += vertexnum;
	  if( edgesnum  > maxedges ) maxedges = edgesnum;
	  sumedges += edgesnum;
	
	  ++atomcount;

	  /* Order vertices */
	  ord[0] = 0;
	  for ( i=0; i<vertexnum-1; ++i)
	    {
	      maxcos = -1.0;
	      for (j=0; j<vertexnum; ++j)
		{
		  sin = vertex[j].x * vertex[ord[i]].y - vertex[j].y * vertex[ord[i]].x;
		  if ( sin > TOL )
		    {
		      cos = SPROD( vertex[j], vertex[ord[i]] )/ sqrt(SPROD(vertex[j],vertex[j]))/ sqrt(SPROD(vertex[ord[i]],vertex[ord[i]]));
		      if ( cos > maxcos ) 
			{
			  maxcos = cos;
			  minj   = j;
			}
		    }
		}

	      ord[i+1] = minj;

	    }

	  /* Compute area of voronoi cell */
	  for (i=0; i<vertexnum-1; ++i)
	    area += 0.5 * (vertex[ord[i+1]].x * vertex[ord[i]].y - vertex[ord[i+1]].y * vertex[ord[i]].x);
	  
	  area += 0.5 * (vertex[ord[0]].x * vertex[ord[vertexnum-1]].y - vertex[ord[0]].y * vertex[ord[vertexnum-1]].x);

	} /* number of edges == number of vertices */

      else area = 0.0;

    } /* vertesnum < 3 */

  free(edges);
								  
}

#endif






