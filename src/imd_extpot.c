
/******************************************************************************
 *
 * IMD -- The ITAP Molecular Dynamics Program
 *
 * Copyright 1996-2008 Institute for Theoretical and Applied Physics,
 * University of Stuttgart, D-70550 Stuttgart
 *
 ******************************************************************************/

/******************************************************************************
 *
 * imd_extpot.c -- external potentials
 *
 ******************************************************************************/

/******************************************************************************
 * $Revision$
 * $Date$
 ******************************************************************************/

#include "imd.h"

/******************************************************************************
 *
 * compute forces due to external potentials
 *
 ******************************************************************************/

#define  UPPER_EXP (20.0)
#define  LOWER_EXP (0.03)

void calc_extpot(void)
{
  int k, i, n;
  int isinx,isiny,isinz;
  real tmpvec1[4], tmpvec2[4];
  
  for (k=0; k<ep_n; k++) {
    ep_fext[k] = 0.0;
    ep_xmax[k] = 0.0;
    ep_ymax[k] = 0.0;
    ep_xmin[k] = 1.e8;
    ep_ymin[k] = 1.e8;
  }
  
  if(ep_key == 0) {   /* default: original harmonic potential */
    for (k=0; k<NCELLS; ++k) {
      cell *p = CELLPTR(k);
      for (i=0; i<p->n; ++i) {
	for (n=0; n<ep_n; ++n) {
	  
	  isinx= ep_dir[n].x;                  
	  isiny= ep_dir[n].y;                  
	  isinz= ep_dir[n].x;                  
	  
	  vektor d;
	  real   dn;
	  d.x = ep_pos[n].x - ORT(p,i,X); 
	  d.y = ep_pos[n].y - ORT(p,i,Y); 
	  d.z = ep_pos[n].z - ORT(p,i,Z); 
	  dn  = SPROD(d,ep_dir[n]);
	  
	  /* spherical indentor*/
	  if (n<ep_nind) {
	    if (dn > -ep_rcut) {
	      real d2 = SPROD(d,d);
	      real d1 = SQRT(d2);
	      real dd = ep_rcut - d1;
	      if (dd > 0.0) {
		real f = ep_a * dd * dd / d1;   /* force on atoms and indentor */
		KRAFT(p,i,X) -= f * d.x; 
		KRAFT(p,i,Y) -= f * d.y; 
		KRAFT(p,i,Z) -= f * d.z;
		ep_fext[n]   += f * ABS(dn); /* normal force on indentor */
		
		/* for determination of contact area */
		if(isinz)
		  {				  
		    ep_xmax[n] = MAX(ep_xmax[n], ORT(p,i,X) );    
		    ep_ymax[n] = MAX(ep_ymax[n], ORT(p,i,Y) );    
		    ep_xmin[n] = MIN(ep_xmin[n], ORT(p,i,X) );    
		    ep_ymin[n] = MIN(ep_ymin[n], ORT(p,i,Y) );    
		  }
		else if(isiny)
		  {
		    ep_xmax[n] = MAX(ep_xmax[n], ORT(p,i,X) );    
		    ep_ymax[n] = MAX(ep_ymax[n], ORT(p,i,Z) );    
		    ep_xmin[n] = MIN(ep_xmin[n], ORT(p,i,X) );    
		    ep_ymin[n] = MIN(ep_ymin[n], ORT(p,i,Z) );    
		  }
		else
		  {
		    ep_xmax[n] = MAX(ep_xmax[n], ORT(p,i,Y) );    
		    ep_ymax[n] = MAX(ep_ymax[n], ORT(p,i,Z) );    
		    ep_xmin[n] = MIN(ep_xmin[n], ORT(p,i,Y) );    
		    ep_ymin[n] = MIN(ep_ymin[n], ORT(p,i,Z) );    
		  }
	      }
	    }
	  }          
	  /*  potential wall */
	  else {  
	    if (dn*dn < ep_rcut*ep_rcut) {
	      real d1 = (dn>0) ? dn : -1*dn ;
	      real dd = ep_rcut - d1;
	      if (dd > 0.0) {
		real f = ep_a * dd * dd / d1;   /* force on atoms and indentor */
		KRAFT(p,i,X) += f * ep_dir[n].x;
		KRAFT(p,i,Y) += f * ep_dir[n].y;
		KRAFT(p,i,Z) += f * ep_dir[n].z;
		ep_fext[n]   += f;  /* magnitude of force on wall */
	      }
	    } 
	  }
	}
      }
    }
  }
  else if(ep_key == 1) /* Ju Li's spherical indenter, see PRB 67, 104105 */
    {
      vektor d,addforce,totaddforce;
      real   dd,cc;
      real   dn,ddn,ee;
      
      totaddforce.x=0.0;
      totaddforce.y=0.0;
      totaddforce.z=0.0;
      
      for (k=0; k<NCELLS; ++k) {
	cell *p = CELLPTR(k);
	for (i=0; i<p->n; ++i) {
	  for (n=0; n<ep_n; ++n) {
	    
	    isinx= ep_dir[n].x;                  
	    isiny= ep_dir[n].y;                  
	    isinz= ep_dir[n].x;      
	    
	    if (n<ep_nind) {
	      vektor d;
	      real   dn;
	      d.x =   ORT(p,i,X)-ep_pos[n].x; 
	      d.y =   ORT(p,i,Y)-ep_pos[n].y; 
	      d.z =   ORT(p,i,Z)-ep_pos[n].z; 
	      dn  = SPROD(d,ep_dir[n]);
	      dd  = SPROD(d,d);
	      
	      if ( dd < ep_rcut*ep_rcut)
		{
		  /* for the determination of the contact area */
		  if(isinz)
		    {				  
		      ep_xmax[n] = MAX(ep_xmax[n], ORT(p,i,X) );    
		      ep_ymax[n] = MAX(ep_ymax[n], ORT(p,i,Y) );    
		      ep_xmin[n] = MIN(ep_xmin[n], ORT(p,i,X) );    
		      ep_ymin[n] = MIN(ep_ymin[n], ORT(p,i,Y) );    
		    }
		  else if(isiny)
		    {
		      ep_xmax[n] = MAX(ep_xmax[n], ORT(p,i,X) );    
		      ep_ymax[n] = MAX(ep_ymax[n], ORT(p,i,Z) );    
		      ep_xmin[n] = MIN(ep_xmin[n], ORT(p,i,X) );    
		      ep_ymin[n] = MIN(ep_ymin[n], ORT(p,i,Z) );    
		    }
		  else
		    {
		      ep_xmax[n] = MAX(ep_xmax[n], ORT(p,i,Y) );    
		      ep_ymax[n] = MAX(ep_ymax[n], ORT(p,i,Z) );    
		      ep_xmin[n] = MIN(ep_xmin[n], ORT(p,i,Y) );    
		      ep_ymin[n] = MIN(ep_ymin[n], ORT(p,i,Z) );    
		  }
		  
		  ddn= sqrt(dd);
		  cc = (ep_rcut - ddn)/ep_a;
		  if (cc > UPPER_EXP) cc = UPPER_EXP;
		  if (cc < LOWER_EXP) cc = LOWER_EXP;
		  ee = exp(cc - 1.0/cc);
			
		  tot_pot_energy += ee;
		  POTENG(p,i) += ee;
			
		  ee = ee / ep_a / ddn * (1.0 + 1.0 /(cc*cc));
                        
		  KRAFT(p,i,X) += ee * d.x; 
		  KRAFT(p,i,Y) += ee * d.y; 
		  KRAFT(p,i,Z) += ee * d.z;
			
		  totaddforce.x += ee * d.x;
		  totaddforce.y += ee * d.y;
		  totaddforce.z += ee * d.z; 
                  
		  ep_fext[n]   += ee * ABS(dn); /* normal force on indentor */
		  
		}
	    } 
	  }
	}
      }
      
#ifdef MPI
      tmpvec1[0] =  totaddforce.x ;
      tmpvec1[1] =  totaddforce.y ;
      tmpvec1[2] =  totaddforce.z ;
      //    printf("before totaddforcereduce allreduce\n");fflush(stdout);
      MPI_Allreduce( tmpvec1, tmpvec2, 4, REAL, MPI_SUM, cpugrid);
      // printf("after totaddforce allreduce\n");fflush(stdout);
      totaddforce.x = tmpvec2[0];
      totaddforce.y = tmpvec2[1];
      totaddforce.z = tmpvec2[2];      
#endif
      /* no need for a wall as the total additional impuls is substracted */
      
      totaddforce.x *= 1.0/nactive_vect[0];
      totaddforce.y *= 1.0/nactive_vect[1];
      totaddforce.z *= 1.0/nactive_vect[2];
      
      
      for (k=0; k<NCELLS; ++k) {
	cell *p = CELLPTR(k);
	for (i=0; i<p->n; ++i) {
	  KRAFT(p,i,X) -= totaddforce.x;
	  KRAFT(p,i,Y) -= totaddforce.y;
	  KRAFT(p,i,Z) -= totaddforce.z;
	}
      }
    }
  
  else if(ep_key == 2) 
    /* Ju Li's spherical indenter made flat, see PRB 67, 104105
       with subtraction of total additional impulse
       works only with indentation directions parallel to box vectors*/ 
    {
      vektor d,addforce,totaddforce;
      real   dd,cc;
      real   dn,ddn,ee;
      
      totaddforce.x=0.0;
      totaddforce.y=0.0;
      totaddforce.z=0.0;
      
      
      for (k=0; k<NCELLS; ++k) {
	cell *p = CELLPTR(k);
	for (i=0; i<p->n; ++i) {
	  for (n=0; n<ep_n; ++n) {
	    
	    isinx= ep_dir[n].x;                  
	    isiny= ep_dir[n].y;                  
	    isinz= ep_dir[n].x;      
	    
	    vektor d;
	    real   dn;
	    d.x =  (ep_dir[n].x==0)  ? 0 : (ORT(p,i,X)-ep_pos[n].x) ; 
	    d.y =   (ep_dir[n].y==0)  ? 0 : (ORT(p,i,Y)-ep_pos[n].y); 
	    d.z =   (ep_dir[n].z==0)  ? 0 : (ORT(p,i,Z)-ep_pos[n].z); 
	    dn  = SPROD(d,ep_dir[n]);
	    dd  = SPROD(d,d);
		
	    if ( dd < ep_rcut*ep_rcut)
	      {
		/* for the determination of the contact area */
		if(isinz)
		    {				  
		      ep_xmax[n] = MAX(ep_xmax[n], ORT(p,i,X) );    
		      ep_ymax[n] = MAX(ep_ymax[n], ORT(p,i,Y) );    
		      ep_xmin[n] = MIN(ep_xmin[n], ORT(p,i,X) );    
		      ep_ymin[n] = MIN(ep_ymin[n], ORT(p,i,Y) );    
		    }
		  else if(isiny)
		    {
		      ep_xmax[n] = MAX(ep_xmax[n], ORT(p,i,X) );    
		      ep_ymax[n] = MAX(ep_ymax[n], ORT(p,i,Z) );    
		      ep_xmin[n] = MIN(ep_xmin[n], ORT(p,i,X) );    
		      ep_ymin[n] = MIN(ep_ymin[n], ORT(p,i,Z) );    
		    }
		  else
		    {
		      ep_xmax[n] = MAX(ep_xmax[n], ORT(p,i,Y) );    
		      ep_ymax[n] = MAX(ep_ymax[n], ORT(p,i,Z) );    
		      ep_xmin[n] = MIN(ep_xmin[n], ORT(p,i,Y) );    
		      ep_ymin[n] = MIN(ep_ymin[n], ORT(p,i,Z) );    
		    }

		ddn= sqrt(dd);
		cc = (ep_rcut - ddn)/ep_a;
		if (cc > UPPER_EXP) cc = UPPER_EXP;
		if (cc < LOWER_EXP) cc = LOWER_EXP;
		ee = exp(cc - 1.0/cc);
		    
		tot_pot_energy += ee;
		POTENG(p,i) += ee;
		    
		ee = ee / ep_a / ddn * (1.0 + 1.0 /(cc*cc));
                    
		KRAFT(p,i,X) += ee * d.x; 
		KRAFT(p,i,Y) += ee * d.y; 
		KRAFT(p,i,Z) += ee * d.z;
		    
		totaddforce.x += ee * d.x;
		totaddforce.y += ee * d.y;
		totaddforce.z += ee * d.z; 
                    
		ep_fext[n]   += ee * ABS(dn); /* normal force on indentor */
                    
	      }
		
	  } 
	      
	}
      }
      
      
#ifdef MPI
      tmpvec1[0] =  totaddforce.x ;
      tmpvec1[1] =  totaddforce.y ;
      tmpvec1[2] =  totaddforce.z ;
      //    printf("before totaddforcereduce allreduce\n");fflush(stdout);
      MPI_Allreduce( tmpvec1, tmpvec2, 4, REAL, MPI_SUM, cpugrid);
      // printf("after totaddforce allreduce\n");fflush(stdout);
      totaddforce.x = tmpvec2[0];
      totaddforce.y = tmpvec2[1];
      totaddforce.z = tmpvec2[2];      
#endif
      /* no need for a wall as the total additional impuls is substracted */
      
      totaddforce.x *= 1.0/nactive_vect[0];
      totaddforce.y *= 1.0/nactive_vect[1];
      totaddforce.z *= 1.0/nactive_vect[2];
      
      
      for (k=0; k<NCELLS; ++k) {
	cell *p = CELLPTR(k);
	for (i=0; i<p->n; ++i) {
	  KRAFT(p,i,X) -= totaddforce.x;
	  KRAFT(p,i,Y) -= totaddforce.y;
	  KRAFT(p,i,Z) -= totaddforce.z;
	}
      }
    }

  else if(ep_key == 3) 
    /* Ju Li's spherical indenter made flat, see PRB 67, 104105
       without subtraction of the total additional impulse
       works only with indentation directions parallel to box vectors*/ 
    {
      vektor d,addforce,totaddforce;
      real   dd,cc;
      real   dn,ddn,ee;
      
      totaddforce.x=0.0;
      totaddforce.y=0.0;
      totaddforce.z=0.0;
      
      
      for (k=0; k<NCELLS; ++k) {
	cell *p = CELLPTR(k);
	for (i=0; i<p->n; ++i) {
	  for (n=0; n<ep_n; ++n) {
	    
	    isinx= ep_dir[n].x;                  
	    isiny= ep_dir[n].y;                  
	    isinz= ep_dir[n].x;      
	    
	    vektor d;
	    real   dn;
	    d.x =  (ep_dir[n].x==0)  ? 0 : (ORT(p,i,X)-ep_pos[n].x) ; 
	    d.y =   (ep_dir[n].y==0)  ? 0 : (ORT(p,i,Y)-ep_pos[n].y); 
	    d.z =   (ep_dir[n].z==0)  ? 0 : (ORT(p,i,Z)-ep_pos[n].z); 
	    dn  = SPROD(d,ep_dir[n]);
	    dd  = SPROD(d,d);
		
	    if ( dd < ep_rcut*ep_rcut)
	      {
		/* for the determination of the contact area */
		if(isinz)
		    {				  
		      ep_xmax[n] = MAX(ep_xmax[n], ORT(p,i,X) );    
		      ep_ymax[n] = MAX(ep_ymax[n], ORT(p,i,Y) );    
		      ep_xmin[n] = MIN(ep_xmin[n], ORT(p,i,X) );    
		      ep_ymin[n] = MIN(ep_ymin[n], ORT(p,i,Y) );    
		    }
		  else if(isiny)
		    {
		      ep_xmax[n] = MAX(ep_xmax[n], ORT(p,i,X) );    
		      ep_ymax[n] = MAX(ep_ymax[n], ORT(p,i,Z) );    
		      ep_xmin[n] = MIN(ep_xmin[n], ORT(p,i,X) );    
		      ep_ymin[n] = MIN(ep_ymin[n], ORT(p,i,Z) );    
		    }
		  else
		    {
		      ep_xmax[n] = MAX(ep_xmax[n], ORT(p,i,Y) );    
		      ep_ymax[n] = MAX(ep_ymax[n], ORT(p,i,Z) );    
		      ep_xmin[n] = MIN(ep_xmin[n], ORT(p,i,Y) );    
		      ep_ymin[n] = MIN(ep_ymin[n], ORT(p,i,Z) );    
		    }

		ddn= sqrt(dd);
		cc = (ep_rcut - ddn)/ep_a;
		if (cc > UPPER_EXP) cc = UPPER_EXP;
		if (cc < LOWER_EXP) cc = LOWER_EXP;
		ee = exp(cc - 1.0/cc);
		    
		tot_pot_energy += ee;
		POTENG(p,i) += ee;
		    
		ee = ee / ep_a / ddn * (1.0 + 1.0 /(cc*cc));
		
		KRAFT(p,i,X) += ee * d.x; 
		KRAFT(p,i,Y) += ee * d.y; 
		KRAFT(p,i,Z) += ee * d.z;
		ep_fext[n]   += ee * ABS(dn); /* normal force on indentor */
	      }	    
	  } 
	}
      }     
    }
  
  else
    {
      error("Error: external potential ep_key not defined.\n");
    }
}

/******************************************************************************
 *
 * move external potentials
 *
 ******************************************************************************/

void move_extpot(real factor)
{
  int n;
  for (n=0; n<ep_n; ++n) {
    ep_pos[n].x += factor * ep_vel[n].x;
    /*
      if (ep_pos[n].x < 0.0    ) ep_pos[n].x += box_x.x;
      if (ep_pos[n].x > box_x.x) ep_pos[n].x -= box_x.x;
    */
    ep_pos[n].y += factor * ep_vel[n].y;
    /*
      if (ep_pos[n].y < 0.0    ) ep_pos[n].y += box_y.y;
      if (ep_pos[n].y > box_y.y) ep_pos[n].y -= box_y.y;
    */
    ep_pos[n].z += factor * ep_vel[n].z;
  }
}

void init_extpot(void)
{
  long tmpvec1[3], tmpvec2[3];
  int i,k,sort;
  nactive_vect[0]=0;
  nactive_vect[1]=0;
  nactive_vect[2]=0;
  if (0==myid)
    {
      printf( "EXTPOT: choosen potential ep_key = %d\n", ep_key);
      printf( "EXTPOT: number of indenters = %d\n", ep_n);
      printf( "EXTPOT: number of walls = %d\n", ep_n-ep_nind);
      printf( "EXTPOT: external potential constant = %f\n",ep_a );
      printf( "EXTPOT: cutoff radius = %f\n", ep_rcut );
      for (i=0; i<ep_n; i++)
        {
	  printf("EXTPOT: ep_pos #%d   %.10g %.10g %.10g \n",
		 i, ep_pos[i].x, ep_pos[i].y, ep_pos[i].z);
	  printf("EXTPOT: ep_vel #%d   %.10g %.10g %.10g \n",
		 i, ep_vel[i].x, ep_vel[i].y, ep_vel[i].z);
	  printf("EXTPOT: ep_dir #%d   %.10g %.10g %.10g \n",
		 i, ep_dir[i].x, ep_dir[i].y, ep_dir[i].z);
        }
    }
  /* needed if the indenter impulse should be balanced */
  /* might not work with 2D, epitax, clone, superatom  */
  if(ep_key==1 || ep_key==2  ){
    for (k=0; k<NCELLS; ++k) {
      cell *p = CELLPTR(k);
      for (i=0; i<p->n; ++i) {
	sort = VSORTE(p,i);
	nactive_vect[0] += (restrictions + sort)->x;
	nactive_vect[1] += (restrictions + sort)->y;
	nactive_vect[2] += (restrictions + sort)->z;
      }
    }
#ifdef MPI
    tmpvec1[0] =  nactive_vect[0] ;
    tmpvec1[1] =  nactive_vect[1] ;
    tmpvec1[2] =  nactive_vect[2] ;
    MPI_Allreduce( tmpvec1, tmpvec2, 3, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    nactive_vect[0] = tmpvec2[0];
    nactive_vect[1] = tmpvec2[1];
    nactive_vect[2] = tmpvec2[2];
#endif
    if (0==myid)
      {
	if (nactive_vect[0] == 0 ||  nactive_vect[1] == 0 || nactive_vect[2] == 0)
	  error ("ep_key=1 requires atoms free to move in all directions\n");
	printf("EXTPOT: active degrees of freedom: x %ld y %ld z %ld\n",nactive_vect[0],nactive_vect[1],nactive_vect[2]);
      }
#ifdef RELAX
    printf( "EXTPOT: max number of relaxation steps = %d\n", ep_max_int);
    printf( "EXTPOT: ekin_threshold = %f\n", glok_ekin_threshold);
    printf( "EXTPOT: fnorm_threshold = %f\n", fnorm_threshold);
#endif
        
      
  }
}
