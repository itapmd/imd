
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2004 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_qc.c -- generates a quasicrystal
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/* prototypes needed only here */
real r2 (real ai, real aj, real bi, real bj, real b, real ci, real cj,
	 real c, int ks, real gi, real gj);
real r3 (real ai, real aj, real ak, real bi, real bj, real bk, real b,
	 real ci, real cj, real ck, int ks1, int ks2, real gi, real gj,
	 real gk);
real det(real ai, real aj, real ak, real bi, real bj, real bk, real ci, 
	 real cj, real ck);
void sortin (int ifeld[]);
void adjust(void);
void decorate(int i, int j, int k);
void locate(real x, real y, real z, int i, int j, int k);
void borders(void);


/******************************************************************************
* 
* generate quasicrystal approximant box size from approximant order
*
******************************************************************************/

void init_qc(void)    
{
  integer p[3],q[3],hv,i,j;
  real    tau[3],period[3],tautl;
  integer np,no,na,nb,nc,nt;

  if (size_per_cpu) {
    box_param.x *= cpu_dim.x;
    box_param.y *= cpu_dim.y;
    box_param.z *= cpu_dim.z;
  }
  if ((box_param.x==0) || (box_param.y==0) || (box_param.z==0))
    error("box_param not set!");

  appr[0] = box_param.x;
  appr[1] = box_param.y; 
  appr[2] = box_param.z;

  /* defining constants */
  tautl = (sqrt(5.0)+1.0)*0.5;

  if (0==myid) printf("Approximant %d %d %d\n\n",appr[0],appr[1],appr[2]);

  /* generate fibonacci numbers */

  for (j=0;j<3;j++) {
      p[j] =1;q[j]=0;hv=0;
      for (i=0;i<appr[j];i++) {
	  hv=q[j];
          q[j]=p[j];
          p[j]=p[j]+hv;
      }
      tau[j]=((real) p[j])/((real) q[j]);
      period[j]=4*(tautl*p[j]+q[j])/sqrt(tautl+2.0);

      if (0==myid) {
         printf("p(%1d)= %3d,q(%1d)= %3d,tau(%1d)= %3f\n",
	        j,p[j],j,q[j],j,tau[j]);
         printf("period in atomic units %10.5f\n",period[j]);
      }
  }

  /* analytical computation of size */

  no=4*(p[0]*p[1]*p[2]+p[2]*q[0]*q[1]+p[1]*q[0]*q[2]+
	p[0]*q[1]*q[2]-q[0]*q[1]*q[2]);
  np=4*(p[0]*p[1]*p[2]+p[1]*p[2]*q[0]+p[0]*p[2]*q[1]+
	p[0]*p[1]*q[2]+q[0]*q[1]*q[2]);
  na=np+no;nb=3*na;nc=2*np;nt=na+nb+nc;

  if (0==myid) {
     printf("\nnumber of prolate rhombohedra (p) %d\n",np);
     printf("number of oblate rhombohedra (o) %d\n",no);
     printf("number of vertex atoms (a) %d\n",na);
     printf("number of edge atoms (b) %d\n",nb);
     printf("number of large atoms (c) %d\n",nc);
     printf("total number %d\n\n ",nt);
  }
  
  box_x.x=period[0];
  box_y.y=period[1];
  box_z.z=period[2];

  box_x.y=0;
  box_x.z=0;
  box_y.x=0;
  box_y.z=0;
  box_z.x=0;
  box_z.y=0;
  make_box();
}

/******************************************************************************
* 
* generate quasicrystal approximant box size from approximant order
*
******************************************************************************/

void generate_qc( void )
{
  real c=0.80,d=0.5;                        /* constants */
  int a1,a2,a3,a4,a5,a6,a7,a8;              /* grid boundaries */
  int hv,i,j;                               /* auxiliary variables */
  int p[3],q[3];                            /* approximant */
  int no,np,na,nb,nc,nt;                    /* geometry numbers */
  real tau[3],betrag[3],tau0[3],tau1[3];    /* grid vectors */
  real perkah[3];                           /* period in grid space */
  real tautl;                               /* golden number */
  real tautl0,tautl1,betrtl;                /* tiling vectors (lt) */
  real box[3];
  vektor cen;                               /* cell center */
  vektor perpro;                            /* per pe */

#ifndef BUFCELLS

  ivektor cpu_dim;
  int num_cpus = 0;
  int myid = 0;
  ivektor my_coord;

  /* this can't possibly work... */
  cpu_dim.x=1;
  cpu_dim.y=1;
  cpu_dim.z=1;

  /* this can't possibly work... */
  my_coord.x=0;
  my_coord.y=0;
  my_coord.z=0;

#endif  

  tautl = (sqrt(5.0)+1.0)*0.5;
  gam[0]=0.14;gam[1]=-0.25;gam[2]=0.33;gam[3]=-0.41;gam[4]=0.52;gam[5]=-0.33; 
  num_sort[0]=0;num_sort[1]=0;num_sort[2]=0;

 /* generating grid vectors */

 for (j=0;j<3;j++)
   {
      p[j] =1;q[j]=0;hv=0;
      for (i=0;i<appr[j];i++)
	{
	  hv=q[j];
          q[j]=p[j];
          p[j]=p[j]+hv;
	}
      tau[j]=((real) p[j])/((real) q[j]);
      perkah[j]=(tautl*p[j]+q[j])/sqrt(tautl+2.0);
#ifdef MPI
 if (0==myid)
   {
      printf("p(%1d)= %3d,q(%1d)= %3d,tau(%1d)= %3f\n",
	     j,p[j],j,q[j],j,tau[j]);
      printf("period in atomic units %10.5f\n",perkah[j]);
   };
#endif
   }

  for (i=0;i<3;i++)
    {
      betrag[i]=sqrt(tau[i]*tau[i]+1.0);
      tau0[i]=tau[i]/betrag[i];
      tau1[i]=1.0/betrag[i];
    } 
  gx[0]= tau0[0];gy[0]= 0;      gz[0]=-tau1[2];
  gx[1]= tau1[0];gy[1]= tau0[1];gz[1]= 0;
  gx[2]= 0;      gy[2]= tau1[1];gz[2]= tau0[2];
  gx[3]= 0;      gy[3]=-tau1[1];gz[3]= tau0[2];
  gx[4]= tau1[0];gy[4]=-tau0[1];gz[4]= 0;
  gx[5]= tau0[0];gy[5]= 0;      gz[5]= tau1[2];

  /* generating sort tiling vectors */
  
  betrtl=sqrt(tautl+2.0);
  tautl0=tautl/betrtl;tautl1=1.0/betrtl;
  tx[0]= tautl0;ty[0]= 0;     tz[0]=-tautl1;
  tx[1]= tautl1;ty[1]= tautl0;tz[1]= 0;
  tx[2]= 0;     ty[2]= tautl1;tz[2]= tautl0;
  tx[3]= 0;     ty[3]=-tautl1;tz[3]= tautl0;
  tx[4]= tautl1;ty[4]=-tautl0;tz[4]= 0;
  tx[5]= tautl0;ty[5]= 0;     tz[5]= tautl1;

  /* distributing computation */

  gmin.x=-perkah[0];gmax.x=perkah[0];
  gmin.y=-perkah[1];gmax.y=perkah[1];
  gmin.z=-perkah[2];gmax.z=perkah[2];
  
  perpro.x=(gmax.x-gmin.x)/((real) cpu_dim.x);
  perpro.y=(gmax.y-gmin.y)/((real) cpu_dim.y);
  perpro.z=(gmax.z-gmin.z)/((real) cpu_dim.z);
  
  for (i=0;i<3;i++) 
    {
      box[i]=4.0*perkah[i];
      iper[i]=1.0/box[i];
    }

      /* Set up 1 atom input cell */

      input = (cell *) malloc(sizeof(cell));
      if (0==input) error("Can't allocate input cell.") ;
      input->n_max=0;
      alloc_cell(input, 1);

      /* Set up cpu parts */

      cen.x=(2.*my_coord.x-cpu_dim.x+1.)*perpro.x*0.5*d;
      lmin.x=cen.x-perpro.x*0.5-c;
      lmax.x=cen.x+perpro.x*0.5+c;
      
      lmin.x=MAX(lmin.x,gmin.x);
      lmax.x=MIN(lmax.x,gmax.x);
      
      perm[0]=2.0*(my_coord.x*perpro.x);
      perp[0]=perm[0]+2.0*perpro.x;
      
      cen.y=(2.*my_coord.y-cpu_dim.y+1.)*perpro.y*0.5*d;
      lmin.y=cen.y-perpro.y*0.5-c;
      lmax.y=cen.y+perpro.y*0.5+c;

      lmin.y=MAX(lmin.y,gmin.y);
      lmax.y=MIN(lmax.y,gmax.y);
      
      perm[1]=2.0*(my_coord.y*perpro.y);
      perp[1]=perm[1]+2.0*perpro.y;
      
      cen.z=(2.*my_coord.z-cpu_dim.z+1.)*perpro.z*0.5*d;
      lmin.z=cen.z-perpro.z*0.5-c;
      lmax.z=cen.z+perpro.z*0.5+c;
      
      lmin.z=MAX(lmin.z,gmin.z);
      lmax.z=MIN(lmax.z,gmax.z);
      
      perm[2]=2.0*(my_coord.z*perpro.z);
      perp[2]=perm[2]+2.0*perpro.z;
      
      /* computing basic grid border parameters */
      
      for (j=0;j<6;j++)
	{
	  a1=floor(gx[j]*lmin.x+gy[j]*lmin.y+gz[j]*lmin.z-gam[j]+0.5);
	  a2=floor(gx[j]*lmax.x+gy[j]*lmax.y+gz[j]*lmax.z-gam[j]+0.5);
	  a3=floor(gx[j]*lmin.x+gy[j]*lmax.y+gz[j]*lmin.z-gam[j]+0.5);
	  a4=floor(gx[j]*lmax.x+gy[j]*lmin.y+gz[j]*lmax.z-gam[j]+0.5);
	  a5=floor(gx[j]*lmin.x+gy[j]*lmin.y+gz[j]*lmax.z-gam[j]+0.5);
	  a6=floor(gx[j]*lmax.x+gy[j]*lmax.y+gz[j]*lmin.z-gam[j]+0.5);
	  a7=floor(gx[j]*lmin.x+gy[j]*lmax.y+gz[j]*lmax.z-gam[j]+0.5);
	  a8=floor(gx[j]*lmax.x+gy[j]*lmin.y+gz[j]*lmin.z-gam[j]+0.5);
	  k1min[j]=MIN(MIN(MIN(a1,a2),MIN(a3,a4)),MIN(MIN(a5,a6),MIN(a7,a8)));
	  k1max[j]=MAX(MAX(MAX(a1,a2),MAX(a3,a4)),MAX(MAX(a5,a6),MAX(a7,a8)));
	}      
      natoms=0;
      nactive=0;
      
      borders();
      
      adjust();

      printf("Number of PE, number of atoms: %d %d\n",myid,natoms);
}


/******************************************************************************
*
* computation of the grid boundaries
*
******************************************************************************/

void borders(void)
{
  int i,j,k,ks1,ks2,ks3;
  int k2min[6],k2max[6],k3min[6],k3max[6];
  real r11,r12,r13,r14,r21,r22,r23,r24,r31,r32,r33,r34,x,y,z;
  real rmin1,rmax1,rmin2,rmax2,rmin3,rmax3,rmin,rmax;
  
  /* through all grid directions */

  for (i=0;i<4;i++)
    for (j=i+1;j<5;j++)
      for (k=j+1;k<6;k++)

	{
	  /* first grid */
	  
	  for (ks1=k1min[i];ks1<k1max[i]+1;ks1++)
	    {
	      kf[i]=ks1;
	      if (gx[i] != 0)
		{
		  r11=r2(gx[i],gx[j],gy[i],gy[j],lmin.y,gz[i],gz[j],lmin.z,
			 ks1,gam[i],gam[j]);
		  r12=r2(gx[i],gx[j],gy[i],gy[j],lmin.y,gz[i],gz[j],lmax.z,
			 ks1,gam[i],gam[j]);	
		  r13=r2(gx[i],gx[j],gy[i],gy[j],lmax.y,gz[i],gz[j],lmin.z,
			 ks1,gam[i],gam[j]);
		  r14=r2(gx[i],gx[j],gy[i],gy[j],lmax.y,gz[i],gz[j],lmax.z,
			 ks1,gam[i],gam[j]);
		  rmin1=MIN(MIN(r11,r12),MIN(r13,r14));
		  rmax1=MAX(MAX(r11,r12),MAX(r13,r14));
		}
	      else
		{
		  rmin1=k1min[j];
		  rmax1=k1max[j];
		}
	      if (gy[i] != 0)
		{
		  r21=r2(gy[i],gy[j],gz[i],gz[j],lmin.z,gx[i],gx[j],lmin.x,
			 ks1,gam[i],gam[j]);
		  r22=r2(gy[i],gy[j],gz[i],gz[j],lmin.z,gx[i],gx[j],lmax.x,
			 ks1,gam[i],gam[j]);
		  r23=r2(gy[i],gy[j],gz[i],gz[j],lmax.z,gx[i],gx[j],lmin.x,
		       ks1,gam[i],gam[j]);
		  r24=r2(gy[i],gy[j],gz[i],gz[j],lmax.z,gx[i],gx[j],lmax.x,
			 ks1,gam[i],gam[j]);
		  rmin2=MIN(MIN(r21,r22),MIN(r23,r24));
		  rmax2=MAX(MAX(r21,r22),MAX(r23,r24));
		}
	      else
		{
		  rmin2=k1min[j];
		  rmax2=k1max[j];
		}
	      if (gz[i] != 0) 
		{
		  r31=r2(gz[i],gz[j],gx[i],gx[j],lmin.x,gy[i],gy[j],lmin.y,
			 ks1,gam[i],gam[j]);
		  r32=r2(gz[i],gz[j],gx[i],gx[j],lmin.x,gy[i],gy[j],lmax.y,
			 ks1,gam[i],gam[j]);
		  r33=r2(gz[i],gz[j],gx[i],gx[j],lmax.x,gy[i],gy[j],lmin.y,
			 ks1,gam[i],gam[j]);
		  r34=r2(gz[i],gz[j],gx[i],gx[j],lmax.x,gy[i],gy[j],lmax.y,
			 ks1,gam[i],gam[j]);
		  rmin3=MIN(MIN(r31,r32),MIN(r33,r34));
		  rmax3=MAX(MAX(r31,r32),MAX(r33,r34));
		}
	      else
		{
		  rmin3=k1min[j];
		  rmax3=k1max[j];
		}

	      rmin=MAX(MAX(rmin1,rmin2),rmin3);
	      rmax=MIN(MIN(rmax1,rmax2),rmax3);

	      k2min[j]=floor(rmin+0.5);
	      k2max[j]=floor(rmax+0.5);
	      
	      /* second grid */
	      
	      for (ks2=k2min[j];ks2<k2max[j]+1;ks2++)
		{
		  kf[j]=ks2;
		  if (gx[i]*gz[j]-gx[j]*gz[i] != 0) 
		    {
		      r11=r3(gx[i],gx[j],gx[k],gy[i],gy[j],gy[k],lmin.y,
			     gz[i],gz[j],gz[k],ks1,ks2,gam[i],gam[j],
			     gam[k]);
		      r12=r3(gx[i],gx[j],gx[k],gy[i],gy[j],gy[k],lmax.y,
			     gz[i],gz[j],gz[k],ks1,ks2,gam[i],gam[j],
			     gam[k]);
		      rmin1=MIN(r11,r12);
		      rmax1=MAX(r11,r12);
		    }
		  else
		    {
		      rmin1=k1min[k];
		      rmax1=k1max[k];
		    }                   
		  if (gx[i]*gy[j]-gx[j]*gy[i] != 0) 
		    {
		      r21=r3(gy[i],gy[j],gy[k],gz[i],gz[j],gz[k],lmin.z,
			     gx[i],gx[j],gx[k],ks1,ks2,gam[i],gam[j],
			     gam[k]);
		      r22=r3(gy[i],gy[j],gy[k],gz[i],gz[j],gz[k],lmax.z,
			     gx[i],gx[j],gx[k],ks1,ks2,gam[i],gam[j],
			     gam[k]);
		      rmin2=MIN(r21,r22);
		      rmax2=MAX(r21,r22);
		    }
		  else
		    {
		      rmin2=k1min[k];
		      rmax2=k1max[k];
		    }
		  if (gy[i]*gz[j]-gy[j]*gz[i] != 0) 
		    {
		      r31=r3(gz[i],gz[j],gz[k],gx[i],gx[j],gx[k],lmin.x,
			     gy[i],gy[j],gy[k],ks1,ks2,gam[i],gam[j],
			     gam[k]);
		      r32=r3(gz[i],gz[j],gz[k],gx[i],gx[j],gx[k],lmax.x,
			     gy[i],gy[j],gy[k],ks1,ks2,gam[i],gam[j],
			     gam[k]);
		      rmin3=MIN(r31,r32);
		      rmax3=MAX(r31,r32);
		    }  
		  else 
		    {
		      
		      rmin3=k1min[k];
		      rmax3=k1max[k];
		    }
		  rmin=MAX(MAX(rmin1,rmin2),rmin3);
		  rmax=MIN(MIN(rmax1,rmax2),rmax3);
		  k3min[k]=floor(rmin+0.5);
		  k3max[k]=floor(rmax+0.5);
		  
		  /* third grid */
		  
		  for (ks3=k3min[k];ks3<k3max[k]+1;ks3++)
		    {
		      kf[k]=ks3;
		      
		      /* compute the intersection point */
		      
		      z=det(gx[i],gx[j],gx[k],gy[i],gy[j],gy[k],
			    ks1+gam[i],ks2+gam[j],ks3+gam[k])/
			det(gx[i],gx[j],gx[k],gy[i],gy[j],gy[k],
			    gz[i],gz[j],gz[k]);
		      y=det(gx[i],gx[j],gx[k],ks1+gam[i],ks2+gam[j],
			    ks3+gam[k],gz[i],gz[j],gz[k])/
			det(gx[i],gx[j],gx[k],gy[i],gy[j],gy[k],
			  gz[i],gz[j],gz[k]);
		      x=det(ks1+gam[i],ks2+gam[j],ks3+gam[k],
			    gy[i],gy[j],gy[k],gz[i],gz[j],gz[k])/
			det(gx[i],gx[j],gx[k],gy[i],gy[j],gy[k],
			    gz[i],gz[j],gz[k]);
		      
		      /* compute the six-dimensional coordinates */

		      locate(x,y,z,i,j,k);
		      
		      /* generate the quasilattice */
		      
		      decorate(i,j,k);
		    }
		}
	    }
	}
}

/******************************************************************************
*
* auxiliary functions
*
******************************************************************************/

/* two grids and an edge */

real r2 (real ai, real aj, real bi, real bj, real b, real ci, real cj, 
	 real c, int ks, real gi, real gj) 
{ 
  real ra2;
  ra2=(bj-aj*bi/ai)*b+(cj-aj*ci/ai)*c+aj/ai*(ks+gi)-gj;
  return ra2;
}

/* one grid and a plane */

real r3 (real ai, real aj, real ak, real bi, real bj, real bk, real b,
	 real ci, real cj, real ck, int ks1, int ks2, real gi, real gj,
	 real gk) 
{
  real ra3;
  ra3=-(((ai*bj-aj*bi)*ck+(ak*bi-ai*bk)*cj+(aj*bk-ak*bj)*ci)*b+
       (ak*ci-ai*ck)*ks2+(aj*ck-ak*cj)*ks1+(gi*aj-gj*ai)*ck+
       (gk*ai-gi*ak)*cj+(gj*ak-gk*aj)*ci)/(ai*cj-aj*ci);
  return ra3;
}  

/* determinant */

real det(real ai, real aj, real ak, real bi, real bj, real bk, real ci, 
	 real cj, real ck)
{
  real deter;
  deter=ai*(bj*ck-bk*cj)-aj*(bi*ck-bk*ci)+ak*(bi*cj-ci*bj);
  return deter;
}

/******************************************************************************
*
* computation of the 6d grid index
*
******************************************************************************/


void locate(real x, real y, real z, int i, int j, int k)
{
  int kp1,kp2,kp3,l;
  real zw;
    
  kp1=kf[i];
  kp2=kf[j];
  kp3=kf[k];
  for (l=0;l<6;l++)
    {
      zw=gx[l]*x+gy[l]*y+gz[l]*z-gam[l];
      kf[l]=ceil(zw);
      /*      if (zw<0.0) kf[l]++; */
    }
  kf[i]=kp1;
  kf[j]=kp2;
  kf[k]=kp3;
}

/******************************************************************************
*
* decoration of the tiling 
*
******************************************************************************/

void decorate(int i, int j, int k)
{
  int l,n,kfeld[22][7],ifeld[7],m;

  for (l=0;l<22;l++)
    for (n=0;n<7;n++)
      kfeld[l][n]=0;

  /* vertices */
      
      kfeld[0][i]=2*kf[i]+2;kfeld[0][j]=2*kf[j];  kfeld[0][k]=2*kf[k];
      kfeld[1][i]=2*kf[i]+2;kfeld[1][j]=2*kf[j]+2;kfeld[1][k]=2*kf[k];
      kfeld[2][i]=2*kf[i];  kfeld[2][j]=2*kf[j]+2;kfeld[2][k]=2*kf[k];
      kfeld[3][i]=2*kf[i];  kfeld[3][j]=2*kf[j];  kfeld[3][k]=2*kf[k]+2;
      kfeld[4][i]=2*kf[i]+2;kfeld[4][j]=2*kf[j];  kfeld[4][k]=2*kf[k]+2;
      kfeld[5][i]=2*kf[i];  kfeld[5][j]=2*kf[j]+2;kfeld[5][k]=2*kf[k]+2;
      kfeld[6][i]=2*kf[i];  kfeld[6][j]=2*kf[j];  kfeld[6][k]=2*kf[k];
      kfeld[7][i]=2*kf[i]+2;kfeld[7][j]=2*kf[j]+2;kfeld[7][k]=2*kf[k]+2;
      
      /* edge atoms */

      kfeld[8][i]=2*kf[i]+1; kfeld[8][j]=2*kf[j];   kfeld[8][k]=2*kf[k];
      kfeld[9][i]=2*kf[i];   kfeld[9][j]=2*kf[j]+1; kfeld[9][k]=2*kf[k];
      kfeld[10][i]=2*kf[i];  kfeld[10][j]=2*kf[j];  kfeld[10][k]=2*kf[k]+1;
      kfeld[11][i]=2*kf[i]+2;kfeld[11][j]=2*kf[j]+2;kfeld[11][k]=2*kf[k]+1;
      kfeld[12][i]=2*kf[i]+2;kfeld[12][j]=2*kf[j]+1;kfeld[12][k]=2*kf[k]+2;
      kfeld[13][i]=2*kf[i]+1;kfeld[13][j]=2*kf[j]+2;kfeld[13][k]=2*kf[k]+2;
      kfeld[14][i]=2*kf[i]+2;kfeld[14][j]=2*kf[j]+1;kfeld[14][k]=2*kf[k];
      kfeld[15][i]=2*kf[i]+1;kfeld[15][j]=2*kf[j]+2;kfeld[15][k]=2*kf[k];
      kfeld[16][i]=2*kf[i]+2;kfeld[16][j]=2*kf[j];  kfeld[16][k]=2*kf[k]+1;
      kfeld[17][i]=2*kf[i];  kfeld[17][j]=2*kf[j]+2;kfeld[17][k]=2*kf[k]+1;
      kfeld[18][i]=2*kf[i];  kfeld[18][j]=2*kf[j]+1;kfeld[18][k]=2*kf[k]+2;
      kfeld[19][i]=2*kf[i]+1;kfeld[19][j]=2*kf[j];  kfeld[19][k]=2*kf[k]+2;
      
      /*large atoms in prolate rhombohedra */
      
      if (i == 0 && j == 1 && k == 3 ) 
	{
	  kfeld[20][0]=2*kf[0]+1;kfeld[20][1]=2*kf[1]+1;kfeld[20][3]=2*kf[3]+1;
	  kfeld[20][2]=2*kf[2]-1;kfeld[20][4]=2*kf[4]-1;kfeld[20][5]=2*kf[5]+1;
	  kfeld[21][0]=2*kf[0]+1;kfeld[21][1]=2*kf[1]+1;kfeld[21][3]=2*kf[3]+1;
	  kfeld[21][2]=2*kf[2]+1;kfeld[21][4]=2*kf[4]+1;kfeld[21][5]=2*kf[5]-1;
	}
      if (i == 0 && j == 1 && k == 5 ) 
	{
	  kfeld[20][0]=2*kf[0]+1;kfeld[20][1]=2*kf[1]+1;kfeld[20][5]=2*kf[5]+1;
	  kfeld[20][2]=2*kf[2]+1;kfeld[20][3]=2*kf[3]-1;kfeld[20][4]=2*kf[4]+1;
	  kfeld[21][0]=2*kf[0]+1;kfeld[21][1]=2*kf[1]+1;kfeld[21][5]=2*kf[5]+1;
	  kfeld[21][2]=2*kf[2]-1;kfeld[21][3]=2*kf[3]+1;kfeld[21][4]=2*kf[4]-1;
	}
      if (i == 0 && j == 2 && k == 3 ) 
	{
	  kfeld[20][0]=2*kf[0]+1;kfeld[20][2]=2*kf[2]+1;kfeld[20][3]=2*kf[3]+1;
	  kfeld[20][1]=2*kf[1]-1;kfeld[20][4]=2*kf[4]-1;kfeld[20][5]=2*kf[5]+1;
	  kfeld[21][0]=2*kf[0]+1;kfeld[21][2]=2*kf[2]+1;kfeld[21][3]=2*kf[3]+1;
	  kfeld[21][1]=2*kf[1]+1;kfeld[21][4]=2*kf[4]+1;kfeld[21][5]=2*kf[5]-1;
	}
      if (i == 0 && j == 2 && k == 4 )
	{
	  kfeld[20][0]=2*kf[0]+1;kfeld[20][2]=2*kf[2]+1;kfeld[20][4]=2*kf[4]+1;
	  kfeld[20][1]=2*kf[1]-1;kfeld[20][3]=2*kf[3]-1;kfeld[20][5]=2*kf[5]+1;
	  kfeld[21][0]=2*kf[0]+1;kfeld[21][2]=2*kf[2]+1;kfeld[21][4]=2*kf[4]+1;
	  kfeld[21][1]=2*kf[1]+1;kfeld[21][3]=2*kf[3]+1;kfeld[21][5]=2*kf[5]-1;
	}
      if (i == 0 && j == 4 && k == 5 )
	{
	  kfeld[20][0]=2*kf[0]+1;kfeld[20][4]=2*kf[4]+1;kfeld[20][5]=2*kf[5]+1;
	  kfeld[20][1]=2*kf[1]+1;kfeld[20][2]=2*kf[2]-1;kfeld[20][3]=2*kf[3]+1;
	  kfeld[21][0]=2*kf[0]+1;kfeld[21][4]=2*kf[4]+1;kfeld[21][5]=2*kf[5]+1;
	  kfeld[21][1]=2*kf[1]-1;kfeld[21][2]=2*kf[2]+1;kfeld[21][3]=2*kf[3]-1;
	}
      if (i == 1 && j == 2 && k == 4 )
	{
	  kfeld[20][1]=2*kf[1]+1;kfeld[20][2]=2*kf[2]+1;kfeld[20][4]=2*kf[4]+1;
	  kfeld[20][0]=2*kf[0]-1;kfeld[20][3]=2*kf[3]-1;kfeld[20][5]=2*kf[5]+1;
	  kfeld[21][1]=2*kf[1]+1;kfeld[21][2]=2*kf[2]+1;kfeld[21][4]=2*kf[4]+1;
	  kfeld[21][0]=2*kf[0]+1;kfeld[21][3]=2*kf[3]+1;kfeld[21][5]=2*kf[5]-1;
	}
      if (i == 1 && j == 2 && k == 5 )
	{
	  kfeld[20][1]=2*kf[1]+1;kfeld[20][2]=2*kf[2]+1;kfeld[20][5]=2*kf[5]+1;
	  kfeld[20][0]=2*kf[0]+1;kfeld[20][3]=2*kf[3]+1;kfeld[20][4]=2*kf[4]-1;
	  kfeld[21][1]=2*kf[1]+1;kfeld[21][2]=2*kf[2]+1;kfeld[21][5]=2*kf[5]+1;
	  kfeld[21][0]=2*kf[0]-1;kfeld[21][3]=2*kf[3]-1;kfeld[21][4]=2*kf[4]+1;
	}
      if (i == 1 && j == 3 && k == 4 )
	{
	  kfeld[20][1]=2*kf[1]+1;kfeld[20][3]=2*kf[3]+1;kfeld[20][4]=2*kf[4]+1;
	  kfeld[20][0]=2*kf[0]-1;kfeld[20][2]=2*kf[2]-1;kfeld[20][5]=2*kf[5]+1;
	  kfeld[21][1]=2*kf[1]+1;kfeld[21][3]=2*kf[3]+1;kfeld[21][4]=2*kf[4]+1;
	  kfeld[21][0]=2*kf[0]+1;kfeld[21][2]=2*kf[2]+1;kfeld[21][5]=2*kf[5]-1;
	}
      if (i == 2 && j == 3 && k == 5 ) 
	{
	  kfeld[20][2]=2*kf[2]+1;kfeld[20][3]=2*kf[3]+1;kfeld[20][5]=2*kf[5]+1;
	  kfeld[20][0]=2*kf[0]-1;kfeld[20][1]=2*kf[1]+1;kfeld[20][4]=2*kf[4]+1;
	  kfeld[21][2]=2*kf[2]+1;kfeld[21][3]=2*kf[3]+1;kfeld[21][5]=2*kf[5]+1;
	  kfeld[21][0]=2*kf[0]+1;kfeld[21][1]=2*kf[1]-1;kfeld[21][4]=2*kf[4]-1;
	}
      if (i == 3 && j == 4 && k == 5 ) 
	{
	  kfeld[20][3]=2*kf[3]+1;kfeld[20][4]=2*kf[4]+1;kfeld[20][5]=2*kf[5]+1;
	  kfeld[20][0]=2*kf[0]+1;kfeld[20][1]=2*kf[1]-1;kfeld[20][2]=2*kf[2]+1;
	  kfeld[21][3]=2*kf[3]+1;kfeld[21][4]=2*kf[4]+1;kfeld[21][5]=2*kf[5]+1;
	  kfeld[21][0]=2*kf[0]-1;kfeld[21][1]=2*kf[1]+1;kfeld[21][2]=2*kf[2]-1;
	}
      
      /* setting the atom type */

      for (l=0;l<8;l++) kfeld[l][6]=1;
      for (l=8;l<20;l++) kfeld[l][6]=2;
      kfeld[20][6]=3;kfeld[21][6]=3;

      /* coordinates for the other grid directions */

      for (l=0;l<6;l++)
	if (l != i && l != j && l != k)
	    for (n=0;n<20;n++) kfeld[n][l]=2*kf[l];
      
      /* collecting the atoms */

      for (l=0;l<22;l++)
	{
	  for (m=0;m<7;m++) ifeld[m]=kfeld[l][m];

	  sortin (ifeld);
	}
}

/******************************************************************************
*
* link cell collecting of atoms 
*
******************************************************************************/
      
void sortin (int ifeld[])
{
  int typ,sign,icell,i,hv,it,to_cpu;
  ivektor cellc;
  real x,y,z,dx,dy,dz,dist;
  cell *p, *q;
  
  x=tx[0]*ifeld[0]+tx[1]*ifeld[1]+tx[2]*ifeld[2]+
    tx[3]*ifeld[3]+tx[4]*ifeld[4]+tx[5]*ifeld[5]+0.1-2.*gmin.x;
  y=ty[0]*ifeld[0]+ty[1]*ifeld[1]+ty[2]*ifeld[2]+
    ty[3]*ifeld[3]+ty[4]*ifeld[4]+ty[5]*ifeld[5]+0.1-2.*gmin.y;
  z=tz[0]*ifeld[0]+tz[1]*ifeld[1]+tz[2]*ifeld[2]+
    tz[3]*ifeld[3]+tz[4]*ifeld[4]+tz[5]*ifeld[5]+0.1-2.*gmin.z;

  if (x < perp[0] && y < perp[1] && z < perp[2] && 
      x > perm[0] && y > perm[1] && z > perm[2]) 
    {   
      cellc = cell_coord(x, y, z);
#ifdef BUFCELLS
      cellc = local_cell_coord(cellc);
#endif
      p = PTR_3D_VV(cell_array,cellc,cell_dim);
      
      hv=1;

      for (i = 0 ; i < p->n; i++) {
	 dx = x - ORT(p,i,X);
	 dy = y - ORT(p,i,Y);
	 dz = z - ORT(p,i,Z);
	 dist = dx*dx+dy*dy+dz*dz;
	
	 if (dist < 0.01) {
	     hv=0;
	     break;
	 }
      }

      if (hv == 1) {

	      natoms++;
              nactive +=3;

	      input->n = 1;
	      ORT(input,0,X)  = x;
	      ORT(input,0,Y)  = y;
	      ORT(input,0,Z)  = z;
	      NUMMER(input,0) = natoms;
	      typ=ifeld[6]-1;

	      if (FABS(x+2.*gmin.x) < 0.0001 && FABS(y+2.*gmin.y) < 0.0001 && 
		  FABS(z+2.*gmin.z) < 0.0001) typ=0;

	      if (typ == 1) typ=0; 
	      if (typ == 2) typ=1;

	      SORTE (input,0) = typ;
	      VSORTE(input,0) = typ;
              MASSE(input,0)  = masses[typ];
              INSERT_ATOM(p, input, 0);
      }
    }
}

/*****************************************************************************
*
* what do we do here?
*
******************************************************************************/

void adjust()
{
  int k;

  for (k=0; k<ncells; ++k) {

    cell *p;
    int  i,typ;
    real x,y,z;

    p = cell_array + CELLS(k);

    for (i=0; i<p->n; ++i) {
	
	x = ORT(p,i,X)-0.1;
	y = ORT(p,i,Y)-0.1;
	z = ORT(p,i,Z)-0.1;
	typ = VSORTE(p,i);
	
	/* fix the type of the first atom */
	if (FABS(x+2.*gmin.x) < 0.0001 && FABS(y+2.*gmin.y) < 0.0001 && 
	    FABS(z+2.*gmin.z) < 0.0001) 
        {
           typ=0;
           SORTE (p,i) = typ;
           VSORTE(p,i) = typ;
        }
        num_sort[typ]++;
    }
  }
}



      
