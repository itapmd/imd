/******************************************************************************
*
* imd_nyeTensorAnalysis_3d.c -- Routines for Nye tensor analysis to identify Burgers vectors
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/
#include "imd.h"

/******************************************************************************
*
* Implements the scheme to approximate Burgers vectors for atoms
* as published in Begau et al. JMPS 60 (2012) 711–722
*
* Supports both single crystalline FCC and BCC materials
*
******************************************************************************/

real neighPerf[14][3];
real neighPerfDistance[14];
int neighPerfLength;
real integrationRadius;

real icoVertices[12][3];
real icoNormals[12][3];

/*
 * Vertices of icosaeder triangles
 * Vertices itself are defined in init_NyeTensor()
 */
int icoFaces[20][3] = {
	{0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},
	{8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},
	{7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6},
	{6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11}
};

nyeTensorInfo* alloc_nyeTensorInfo(void) {
	nyeTensorInfo *info = (nyeTensorInfo*) malloc(sizeof(nyeTensorInfo));

	info->lcm[0][0] = 1.;  info->lcm[1][0] = 0.; info->lcm[2][0] = 0.;
	info->lcm[0][1] = 0.;  info->lcm[1][1] = 1.; info->lcm[2][1] = 0.;
	info->lcm[0][2] = 0.;  info->lcm[1][2] = 0.; info->lcm[2][2] = 1.;

	info->bv[0] = 0.; info->bv[1] = 0.; info->bv[2] = 0.;
	info->ls[0] = 0.; info->ls[1] = 0.; info->ls[2] = 0.;
	return info;
}

/*
 * Free all memory allocated for nyeTensorInfo including the buffer cells
 */
void removeNyeTensorData(){
	int i,j,k,l;
	for (i = 0; i < cell_dim.x; ++i){
		for (j = 0; j < cell_dim.y; ++j){
			for (k = 0; k < cell_dim.z; ++k){
				cell *p = PTR_3D_V(cell_array, i, j, k, cell_dim);
				for (l = 0; l < p->n; l++) {
					if (p->nyeTens[l] != NULL) {
						free(p->nyeTens[l]);
						p->nyeTens[l] = NULL;
					}
				}
			}
		}
	}
}

void init_NyeTensor(){
	int i;

	if (ada_nbr_r2cut == 0.) {
		if (ada_crystal_structure == ADA_FCC_CONFIG) {
			ada_nbr_r2cut = SQR(0.862*ada_latticeConst);
		} else if (ada_crystal_structure == ADA_BCC_CONFIG) {
			ada_nbr_r2cut = SQR(1.22*ada_latticeConst);
		} else if (ada_crystal_structure == ADA_ACKLAND_CONFIG) {
			ada_nbr_r2cut = SQR(0.862*ada_latticeConst);
		}
	}

	if (ada_crystal_structure == ADA_FCC_CONFIG || ada_crystal_structure == ADA_ACKLAND_CONFIG){
		integrationRadius = 0.707107 * ada_latticeConst; /*Nearest neighbor distance*/
		neighPerfLength = 12;
		neighPerf[0][0] = 0.; neighPerf[0][1] = 0.5; neighPerf[0][2] = 0.5;
		neighPerf[1][0] = 0.; neighPerf[1][1] =-0.5; neighPerf[1][2] =-0.5;
		neighPerf[2][0] = 0.; neighPerf[2][1] =-0.5; neighPerf[2][2] = 0.5;
		neighPerf[3][0] = 0.; neighPerf[3][1] = 0.5; neighPerf[3][2] =-0.5;

		neighPerf[4][0] = 0.5; neighPerf[4][1] = 0.; neighPerf[4][2] = 0.5;
		neighPerf[5][0] =-0.5; neighPerf[5][1] = 0.; neighPerf[5][2] =-0.5;
		neighPerf[6][0] =-0.5; neighPerf[6][1] = 0.; neighPerf[6][2] = 0.5;
		neighPerf[7][0] = 0.5; neighPerf[7][1] = 0.; neighPerf[7][2] =-0.5;

		neighPerf[8][0] = 0.5; neighPerf[8][1] = 0.5; neighPerf[8][2] = 0.;
		neighPerf[9][0] =-0.5; neighPerf[9][1] =-0.5; neighPerf[9][2] = 0.;
		neighPerf[10][0]=-0.5; neighPerf[10][1]= 0.5; neighPerf[10][2]= 0.;
		neighPerf[11][0]= 0.5; neighPerf[11][1]=-0.5; neighPerf[11][2]= 0.;
	} else if (ada_crystal_structure == ADA_BCC_CONFIG){
		integrationRadius = 0.866025*ada_latticeConst; /*First nearest neighbor distance*/
		neighPerfLength = 14;

		neighPerf[0][0] = 0.5; 	neighPerf[0][1] = 0.5;	neighPerf[0][2] = 0.5;
		neighPerf[1][0] = 0.5; 	neighPerf[1][1] = 0.5;	neighPerf[1][2] =-0.5;
		neighPerf[2][0] = 0.5; 	neighPerf[2][1] =-0.5;	neighPerf[2][2] = 0.5;
		neighPerf[3][0] = 0.5; 	neighPerf[3][1] =-0.5; 	neighPerf[3][2] =-0.5;
		neighPerf[4][0] =-0.5; 	neighPerf[4][1] = 0.5; 	neighPerf[4][2] = 0.5;
		neighPerf[5][0] =-0.5; 	neighPerf[5][1] = 0.5;	neighPerf[5][2] =-0.5;
		neighPerf[6][0] =-0.5; 	neighPerf[6][1] =-0.5; 	neighPerf[6][2] = 0.5;
		neighPerf[7][0] =-0.5; 	neighPerf[7][1] =-0.5; 	neighPerf[7][2] =-0.5;

		neighPerf[8][0] = 1.; 	neighPerf[8][1] = 0.; 	neighPerf[8][2] = 0.;
		neighPerf[9][0] =-1; 	neighPerf[9][1] = 0.; 	neighPerf[9][2] = 0.;
		neighPerf[10][0]= 0.;	neighPerf[10][1]= 1.; 	neighPerf[10][2]= 0.;
		neighPerf[11][0]= 0.; 	neighPerf[11][1]=-1.;	neighPerf[11][2]= 0.;
		neighPerf[12][0]= 0.;	neighPerf[12][1]= 0.; 	neighPerf[12][2]= 1.;
		neighPerf[13][0]= 0.; 	neighPerf[13][1]= 0.;	neighPerf[13][2]=-1.;
	} else {
		error("Crystal structure not supported by NYETENSOR");
	}

	real rot[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

	rot[0][0] = nye_rotationAxis_x.x/SQRT(SPROD3D(nye_rotationAxis_x, nye_rotationAxis_x));
	rot[1][0] = nye_rotationAxis_x.y/SQRT(SPROD3D(nye_rotationAxis_x, nye_rotationAxis_x));
	rot[2][0] = nye_rotationAxis_x.z/SQRT(SPROD3D(nye_rotationAxis_x, nye_rotationAxis_x));

	rot[0][1] = nye_rotationAxis_y.x/SQRT(SPROD3D(nye_rotationAxis_y, nye_rotationAxis_y));
	rot[1][1] = nye_rotationAxis_y.y/SQRT(SPROD3D(nye_rotationAxis_y, nye_rotationAxis_y));
	rot[2][1] = nye_rotationAxis_y.z/SQRT(SPROD3D(nye_rotationAxis_y, nye_rotationAxis_y));

	rot[0][2] = nye_rotationAxis_z.x/SQRT(SPROD3D(nye_rotationAxis_z, nye_rotationAxis_z));
	rot[1][2] = nye_rotationAxis_z.y/SQRT(SPROD3D(nye_rotationAxis_z, nye_rotationAxis_z));
	rot[2][2] = nye_rotationAxis_z.z/SQRT(SPROD3D(nye_rotationAxis_z, nye_rotationAxis_z));

	real l;

	real icoConstX  = 1.;
	real icoConstZ  = 1.6180339887;

	for (i=0; i<neighPerfLength;i++){
		real n[3];
		n[0] = (neighPerf[i][0] * rot[0][0] + neighPerf[i][1] * rot[1][0] + neighPerf[i][2] * rot[2][0]) * ada_latticeConst;
		n[1] = (neighPerf[i][0] * rot[0][1] + neighPerf[i][1] * rot[1][1] + neighPerf[i][2] * rot[2][1]) * ada_latticeConst;
		n[2] = (neighPerf[i][0] * rot[0][2] + neighPerf[i][1] * rot[1][2] + neighPerf[i][2] * rot[2][2]) * ada_latticeConst;
		neighPerf[i][0] = n[0];
		neighPerf[i][1] = n[1];
		neighPerf[i][2] = n[2];
	}

	for (i=0; i<neighPerfLength; ++i){
		neighPerfDistance[i] = SQRT(SPRODA3D(neighPerf[i],neighPerf[i]));
	}

	//Create icoseader vertices to approximate a sphere with a radius of the perfect neighbor length for FCC
	real r = (2 * integrationRadius) / SQRT(10 + 2 * SQRT(5));

	icoVertices[0][0] = -icoConstX*r; icoVertices[0][1] = 0.; icoVertices[0][2] = icoConstZ*r;
	icoVertices[1][0] = icoConstX*r; icoVertices[1][1] = 0.; icoVertices[1][2] = icoConstZ*r;
	icoVertices[2][0] = -icoConstX*r; icoVertices[2][1] = 0.; icoVertices[2][2] = -icoConstZ*r;
	icoVertices[3][0] = icoConstX*r; icoVertices[3][1] = 0.; icoVertices[3][2] = -icoConstZ*r;
	icoVertices[4][0] = 0.; icoVertices[4][1] = icoConstZ*r; icoVertices[4][2] = icoConstX*r;
	icoVertices[5][0] = 0.; icoVertices[5][1] = icoConstZ*r; icoVertices[5][2] = -icoConstX*r;
	icoVertices[6][0] = 0.; icoVertices[6][1] = -icoConstZ*r; icoVertices[6][2] = icoConstX*r;
	icoVertices[7][0] = 0.; icoVertices[7][1] = -icoConstZ*r; icoVertices[7][2] = -icoConstX*r;
	icoVertices[8][0] = icoConstZ*r; icoVertices[8][1] = icoConstX*r; icoVertices[8][2] = 0.;
	icoVertices[9][0] = -icoConstZ*r; icoVertices[9][1] = icoConstX*r; icoVertices[9][2] = 0.;
	icoVertices[10][0] = icoConstZ*r; icoVertices[10][1] = -icoConstX*r; icoVertices[10][2] = 0.;
	icoVertices[11][0] = -icoConstZ*r; icoVertices[11][1] = -icoConstX*r; icoVertices[11][2] = 0.;

	for (i = 0; i < 12; i++) {
		l = SQRT(icoVertices[i][0]*icoVertices[i][0]
		        +icoVertices[i][1]*icoVertices[i][1]
		        +icoVertices[i][2]*icoVertices[i][2]);
		icoNormals[i][0] = icoVertices[i][0] / l;
		icoNormals[i][1] = icoVertices[i][1] / l;
		icoNormals[i][2] = icoVertices[i][2] / l;
	}

	if (myid == 0){
		if (ada_crystal_structure == ADA_FCC_CONFIG) {
			printf("NYETENSOR: FCC analysis scheme selected \n");
		} else if (ada_crystal_structure == ADA_BCC_CONFIG) {
			printf("NYETENSOR: BCC analysis scheme selected \n");
		} else if (ada_crystal_structure == ADA_ACKLAND_CONFIG) {
			printf("NYETENSOR: ACKLAND analysis scheme selected \n");
		}
		printf("NYETENSOR: lattice constant: %f \n", ada_latticeConst);
		printf("NYETENSOR: given crystal orientation x : %f %f %f\n", nye_rotationAxis_x.x, nye_rotationAxis_x.y, nye_rotationAxis_x.z);
		printf("NYETENSOR: given crystal orientation y : %f %f %f\n", nye_rotationAxis_y.x, nye_rotationAxis_y.y, nye_rotationAxis_y.z);
		printf("NYETENSOR: given crystal orientation z : %f %f %f\n", nye_rotationAxis_z.x, nye_rotationAxis_z.y, nye_rotationAxis_z.z);
	}
}

/*
 * Subroutines to compute the integral of surface in a vector field
 */
real tetraederVolume(real p1[3], real p2[3], real p3[3], real p4[3]){
	return ((p2[0]-p1[0])*((p3[1]-p1[1])*(p4[2]-p1[2])-(p3[2]-p1[2])*(p4[1]-p1[1]))
	          +(p2[1]-p1[1])*((p3[2]-p1[2])*(p4[0]-p1[0])-(p3[0]-p1[0])*(p4[2]-p1[2]))
	          +(p2[2]-p1[2])*((p3[0]-p1[0])*(p4[1]-p1[1])-(p3[1]-p1[1])*(p4[0]-p1[0])))/6.;
}

real integral(real p0[3], real p1[3], real p2[3],
		real p3[3], real p4[3], real p5[3]){
	real volume;
	volume = tetraederVolume(p3, p5, p1, p4);
	volume += tetraederVolume(p3, p5, p2, p1);
	volume += tetraederVolume(p3, p2, p0, p1);

	return volume;
}

int matrixInverse(real m[3][3]){
	real a = m[0][0]; real b = m[0][1]; real c = m[0][2];
	real d = m[1][0]; real e = m[1][1]; real f = m[1][2];
	real g = m[2][0]; real h = m[2][1]; real i = m[2][2];

	real det = (a*e*i + b*f*g + c*d*h - c*e*g - a*f*h - b*d*i);
	if (ABS(det)<0.001) return 0;
	real deti = 1./det;

	m[0][0] = deti * (e*i-h*f); m[0][1] = deti * (c*h-b*i); m[0][2] = deti * (b*f-c*e);
	m[1][0] = deti * (f*g-d*i); m[1][1] = deti * (a*i-c*g); m[1][2] = deti * (c*d-a*f);
	m[2][0] = deti * (d*h-e*g); m[2][1] = deti * (b*g-a*h); m[2][2] = deti * (a*e-b*d);
	return 1;
}

/*
 * Compute the 3x3 lattice correspondence matrix as the least square solution to
 * match the local configuration of nearest neighbor bonds to a reference structure
 */
void calculateLcm(cell *p, int n) {
	int i, j, ii, invertible, best;
	real l, bestAngle, angle;
	int nneigh = NEIGH(p, n)->n;
	real neigh[nneigh][3];
	nyeTensorInfo *nti;
	real *nei;

	for (i = 0; i < nneigh; i++) {
		cell *q = NEIGH(p, n)->cl[i];
		ii = NEIGH(p, n)->num[i];

		neigh[i][0] = ORT(q,ii,X) - ORT(p,n,X);
		neigh[i][1] = ORT(q,ii,Y) - ORT(p,n,Y);
		neigh[i][2] = ORT(q,ii,Z) - ORT(p,n,Z);

		if (pbc_dirs.x) {
			if (neigh[i][0] < -box_x.x * 0.5)
				neigh[i][0] += box_x.x;
			else if (neigh[i][0] > box_x.x * 0.5)
				neigh[i][0] -= box_x.x;
		}
		if (pbc_dirs.y) {
			if (neigh[i][1] < -box_y.y * 0.5)
				neigh[i][1] += box_y.y;
			else if (neigh[i][1] > box_y.y * 0.5)
				neigh[i][1] -= box_y.y;
		}
		if (pbc_dirs.z) {
			if (neigh[i][2] < -box_z.z * 0.5)
				neigh[i][2] += box_z.z;
			else if (neigh[i][2] > box_z.z * 0.5)
				neigh[i][2] -= box_z.z;
		}
	}

	real a[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
	real b[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};

	for (i = 0; i < nneigh; i++) {
		bestAngle = -1.;
		best = 0;
		nei = neigh[i];
		l = SQRT(SPRODA3D(nei,nei));
		for (j = 0; j < neighPerfLength; j++) {
			angle = SPRODA3D(nei, neighPerf[j]) / (l*neighPerfDistance[j]);
			if (angle > bestAngle) {
				best = j;
				bestAngle = angle;
			}
		}
		if (bestAngle > 0.9396926207859084279 ) { /* Angle < 20°, 0.939=cos(20*PI/180.);*/
			a[0][0] += nei[0] * neighPerf[best][0];
			a[0][1] += nei[0] * neighPerf[best][1];
			a[0][2] += nei[0] * neighPerf[best][2];
			a[1][0] += nei[1] * neighPerf[best][0];
			a[1][1] += nei[1] * neighPerf[best][1];
			a[1][2] += nei[1] * neighPerf[best][2];
			a[2][0] += nei[2] * neighPerf[best][0];
			a[2][1] += nei[2] * neighPerf[best][1];
			a[2][2] += nei[2] * neighPerf[best][2];

			b[0][0] += nei[0] * nei[0];
			b[0][1] += nei[0] * nei[1];
			b[0][2] += nei[0] * nei[2];
			b[1][0] += nei[1] * nei[0];
			b[1][1] += nei[1] * nei[1];
			b[1][2] += nei[1] * nei[2];
			b[2][0] += nei[2] * nei[0];
			b[2][1] += nei[2] * nei[1];
			b[2][2] += nei[2] * nei[2];
		}
	}

	invertible = matrixInverse(a);

	nti = NYE(p,n);
	if (invertible) {
		nti->lcm[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0] +a [0][2] * b[2][0];
		nti->lcm[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1] +a [0][2] * b[2][1];
		nti->lcm[0][2] = a[0][0] * b[0][2] + a[0][1] * b[1][2] +a [0][2] * b[2][2];

		nti->lcm[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0] +a [1][2] * b[2][0];
		nti->lcm[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1] +a [1][2] * b[2][1];
		nti->lcm[1][2] = a[1][0] * b[0][2] + a[1][1] * b[1][2] +a [1][2] * b[2][2];

		nti->lcm[2][0] = a[2][0] * b[0][0] + a[2][1] * b[1][0] +a [2][2] * b[2][0];
		nti->lcm[2][1] = a[2][0] * b[0][1] + a[2][1] * b[1][1] +a [2][2] * b[2][1];
		nti->lcm[2][2] = a[2][0] * b[0][2] + a[2][1] * b[1][2] +a [2][2] * b[2][2];
	}
}

void calculateNye(cell *p, int n){
	int i, j, ii, k;
	int nneigh = NEIGH(p, n)->n;
	real neigh[nneigh][3];
	real neighT[nneigh][3][3];
	real grd[3][3][3];
	nyeTensorInfo *nti, *ntiNeigh;
	real e0, de;

	nti = NYE(p,n);

	for (i = 0; i < nneigh; i++) {
		cell *q = NEIGH(p, n)->cl[i];
		ii = NEIGH(p, n)->num[i];

		neigh[i][0] = ORT(q,ii,X) - ORT(p,n,X);
		neigh[i][1] = ORT(q,ii,Y) - ORT(p,n,Y);
		neigh[i][2] = ORT(q,ii,Z) - ORT(p,n,Z);

		if (pbc_dirs.x) {
			if (neigh[i][0] < -box_x.x * 0.5)
				neigh[i][0] += box_x.x;
			else if (neigh[i][0] > box_x.x * 0.5)
				neigh[i][0] -= box_x.x;
		}
		if (pbc_dirs.y) {
			if (neigh[i][1] < -box_y.y * 0.5)
				neigh[i][1] += box_y.y;
			else if (neigh[i][1] > box_y.y * 0.5)
				neigh[i][1] -= box_y.y;
		}
		if (pbc_dirs.z) {
			if (neigh[i][2] < -box_z.z * 0.5)
				neigh[i][2] += box_z.z;
			else if (neigh[i][2] > box_z.z * 0.5)
				neigh[i][2] -= box_z.z;
		}
	}

	for (i = 0; i < nneigh; i++) {
		neighT[i][0][0] = neigh[i][0] * neigh[i][0];
		neighT[i][0][1] = neigh[i][0] * neigh[i][1];
		neighT[i][0][2] = neigh[i][0] * neigh[i][2];

		neighT[i][1][0] = neigh[i][1] * neigh[i][0];
		neighT[i][1][1] = neigh[i][1] * neigh[i][1];
		neighT[i][1][2] = neigh[i][1] * neigh[i][2];

		neighT[i][2][0] = neigh[i][2] * neigh[i][0];
		neighT[i][2][1] = neigh[i][2] * neigh[i][1];
		neighT[i][2][2] = neigh[i][2] * neigh[i][2];
	}

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			e0 = nti->lcm[i][j];
			real a[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
			real c[3] = {0.,0.,0.};
			for (k = 0; k < nneigh; k++) {
				a[0][0] += neighT[k][0][0];
				a[0][1] += neighT[k][0][1];
				a[0][2] += neighT[k][0][2];

				a[1][0] += neighT[k][1][0];
				a[1][1] += neighT[k][1][1];
				a[1][2] += neighT[k][1][2];

				a[2][0] += neighT[k][2][0];
				a[2][1] += neighT[k][2][1];
				a[2][2] += neighT[k][2][2];

				cell *q = NEIGH(p, n)->cl[k];
				ii = NEIGH(p, n)->num[k];
				ntiNeigh = NYE(q,ii);

				de = ntiNeigh->lcm[i][j] - e0;

				c[0] += neigh[k][0] * de;
				c[1] += neigh[k][1] * de;
				c[2] += neigh[k][2] * de;
			}

			matrixInverse(a);
			grd[i][j][0] = c[0] * a[0][0] + c[1] * a[0][1] + c[2] * a[0][2];
			grd[i][j][1] = c[0] * a[1][0] + c[1] * a[1][1] + c[2] * a[1][2];
			grd[i][j][2] = c[0] * a[2][0] + c[1] * a[2][1] + c[2] * a[2][2];
		}
	}


	for (i = 0; i < 3; i++) {
		nti->nyeTensor[0][i] = -grd[2][i][1] + grd[1][i][2];
		nti->nyeTensor[1][i] = -grd[0][i][2] + grd[2][i][0];
		nti->nyeTensor[2][i] = -grd[1][i][0] + grd[0][i][1];
	}
}


int interpolateLineSense(real x, real y, real z, real distSqr,
		cell **neiCells, int *neiIndices, int n, real interpolated[3][3]) {

	real distances[n];	/*Distances from atom to interpolation point*/
	real sum = 0.;
	int i;
	real d_x, d_y, d_z, dis, d;

	interpolated[0][0] = 0.; interpolated[0][1] = 0.; interpolated[0][2] = 0.;
	interpolated[1][0] = 0.; interpolated[1][1] = 0.; interpolated[1][2] = 0.;
	interpolated[2][0] = 0.; interpolated[2][1] = 0.; interpolated[2][2] = 0.;

	for (i=0; i<n;i++) {
		d_x = x - ORT(neiCells[i],neiIndices[i],X);
		d_y = y - ORT(neiCells[i],neiIndices[i],Y);
		d_z = z - ORT(neiCells[i],neiIndices[i],Z);

		if (pbc_dirs.x) {
			if (d_x < -box_x.x * 0.5)
			d_x += box_x.x;
			else if (d_x > box_x.x * 0.5)
			d_x -= box_x.x;
		}
		if (pbc_dirs.y) {
			if (d_y < -box_y.y * 0.5)
			d_y += box_y.y;
			else if (d_y > box_y.y * 0.5)
			d_y -= box_y.y;
		}
		if (pbc_dirs.z) {
			if (d_z < -box_z.z * 0.5)
			d_z += box_z.z;
			else if (d_z > box_z.z * 0.5)
			d_z -= box_z.z;
		}

		dis = d_x*d_x + d_y*d_y + d_z*d_z;
		//Store the distance if it is below the threshold, otherwise save 0.
		if ( dis < distSqr) {
			distances[i] = 1./dis;
			sum += distances[i];
		} else {
			distances[i] = 0.;
		}
	}

	if (sum>0) {
		for (i=0; i<n; i++) {
			nyeTensorInfo *nyeInfo = NYE(neiCells[i],neiIndices[i]);
			interpolated[0][0] += (nyeInfo->nyeTensor[0][0]) * distances[i];
			interpolated[0][1] += (nyeInfo->nyeTensor[1][0]) * distances[i];
			interpolated[0][2] += (nyeInfo->nyeTensor[2][0]) * distances[i];

			interpolated[1][0] += (nyeInfo->nyeTensor[0][1]) * distances[i];
			interpolated[1][1] += (nyeInfo->nyeTensor[1][1]) * distances[i];
			interpolated[1][2] += (nyeInfo->nyeTensor[2][1]) * distances[i];

			interpolated[2][0] += (nyeInfo->nyeTensor[0][2]) * distances[i];
			interpolated[2][1] += (nyeInfo->nyeTensor[1][2]) * distances[i];
			interpolated[2][2] += (nyeInfo->nyeTensor[2][2]) * distances[i];
		}
		d = 1./sum;
		interpolated[0][0] *= d; interpolated[1][0] *= d; interpolated[2][0] *= d;
		interpolated[0][1] *= d; interpolated[1][1] *= d; interpolated[2][1] *= d;
		interpolated[0][2] *= d; interpolated[1][2] *= d; interpolated[2][2] *= d;
	} else return 0;
	return 1;
}

void interpolateBurgersVector(real x, real y, real z, real distSqr,
		cell **neiCells, int *neiIndices, int n, real norm[3], real nyeInter[3]) {

	real distances[n];	/*Distances from atom to interpolation point*/
	real sum = 0.;
	int i;
	real d_x, d_y, d_z, dis, d;

	nyeInter[0] = 0.; nyeInter[1] = 0.; nyeInter[2] = 0.;

	for (i=0; i<n;i++) {
		d_x = x - ORT(neiCells[i],neiIndices[i],X);
		d_y = y - ORT(neiCells[i],neiIndices[i],Y);
		d_z = z - ORT(neiCells[i],neiIndices[i],Z);

		if (pbc_dirs.x) {
			if (d_x < -box_x.x * 0.5)
			d_x += box_x.x;
			else if (d_x > box_x.x * 0.5)
			d_x -= box_x.x;
		}
		if (pbc_dirs.y) {
			if (d_y < -box_y.y * 0.5)
			d_y += box_y.y;
			else if (d_y > box_y.y * 0.5)
			d_y -= box_y.y;
		}
		if (pbc_dirs.z) {
			if (d_z < -box_z.z * 0.5)
			d_z += box_z.z;
			else if (d_z > box_z.z * 0.5)
			d_z -= box_z.z;
		}

		dis = d_x*d_x + d_y*d_y + d_z*d_z;
		//Store the distance if it is below the threshold, otherwise save 0.
		if ( dis < distSqr) {
			distances[i] = 1./dis;
			sum += distances[i];
		} else {
			distances[i] = 0.;
		}
	}

	if (sum>0) {
		for (i=0; i<n; i++) {
			nyeTensorInfo *nyeInfo = NYE(neiCells[i],neiIndices[i]);
			nyeInter[0] += (norm[0]*nyeInfo->nyeTensor[0][0]
			            + norm[1]*nyeInfo->nyeTensor[1][0] + norm[2]*nyeInfo->nyeTensor[2][0]) * distances[i];
			nyeInter[1] += (norm[0]*nyeInfo->nyeTensor[0][1]
			            + norm[1]*nyeInfo->nyeTensor[1][1] + norm[2]*nyeInfo->nyeTensor[2][1]) * distances[i];
			nyeInter[2] += (norm[0]*nyeInfo->nyeTensor[0][2]
			            + norm[1]*nyeInfo->nyeTensor[1][2] + norm[2]*nyeInfo->nyeTensor[2][2]) * distances[i];
		}
		d = 1./sum;
		nyeInter[0] *= d;
		nyeInter[1] *= d;
		nyeInter[2] *= d;
	}
}

int calculateLineSense(real nnDist, cell *central, int centralInd, cell **neiCells, int *neiIndices, int n, real ls[3], int dir[3]) {
	real sqr = nnDist*nnDist;

	real lineSenseInter[12][3][3];
	real x = ORT(central, centralInd, X);
	real y = ORT(central, centralInd, Y);
	real z = ORT(central, centralInd, Z);
	ls[0] = 0.; ls[1] = 0.; ls[2] = 0.;
	int ok = 0;
	int i;

	for (i = 0; i<12; i++) {
		ok = interpolateLineSense(x+icoVertices[i][0], y+icoVertices[i][1], z+icoVertices[i][2], sqr-0.001, neiCells, neiIndices, n, lineSenseInter[i]);
		if (!ok) return 0;
	}

	/*Integrate the Nye tensor field on the surface of a icosahedron, as an approximation of a sphere, to estimate the line direction*/
	for (i = 0; i<20; i++){
		int v1 = icoFaces[i][0];
		int v2 = icoFaces[i][1];
		int v3 = icoFaces[i][2];

		real i1[3];
		real i2[3];
		real i3[3];

		real m = dir[0]*lineSenseInter[v1][0][0]+dir[1]*lineSenseInter[v1][1][0]+dir[2]*lineSenseInter[v1][2][0];
		i1[0] = icoVertices[v1][0]+icoNormals[v1][0]*m;
		i1[1] = icoVertices[v1][1]+icoNormals[v1][1]*m;
		i1[2] = icoVertices[v1][2]+icoNormals[v1][2]*m;

		i2[0] = icoVertices[v2][0]+icoNormals[v2][0]*m;
		i2[1] = icoVertices[v2][1]+icoNormals[v2][1]*m;
		i2[2] = icoVertices[v2][2]+icoNormals[v2][2]*m;

		i3[0] = icoVertices[v3][0]+icoNormals[v3][0]*m;
		i3[1] = icoVertices[v3][1]+icoNormals[v3][1]*m;
		i3[2] = icoVertices[v3][2]+icoNormals[v3][2]*m;

		ls[0] += integral(icoVertices[v1], icoVertices[v2], icoVertices[v3],i1, i2, i3);


		m = dir[0]*lineSenseInter[v1][0][1]+dir[1]*lineSenseInter[v1][1][1]+dir[2]*lineSenseInter[v1][2][1];
		i1[0] = icoVertices[v1][0]+icoNormals[v1][0]*m;
		i1[1] = icoVertices[v1][1]+icoNormals[v1][1]*m;
		i1[2] = icoVertices[v1][2]+icoNormals[v1][2]*m;

		i2[0] = icoVertices[v2][0]+icoNormals[v2][0]*m;
		i2[1] = icoVertices[v2][1]+icoNormals[v2][1]*m;
		i2[2] = icoVertices[v2][2]+icoNormals[v2][2]*m;

		i3[0] = icoVertices[v3][0]+icoNormals[v3][0]*m;
		i3[1] = icoVertices[v3][1]+icoNormals[v3][1]*m;
		i3[2] = icoVertices[v3][2]+icoNormals[v3][2]*m;

		ls[1] += integral(icoVertices[v1], icoVertices[v2], icoVertices[v3],i1, i2, i3);


		m = dir[0]*lineSenseInter[v1][0][2]+dir[1]*lineSenseInter[v1][1][2]+dir[2]*lineSenseInter[v1][2][2];
		i1[0] = icoVertices[v1][0]+icoNormals[v1][0]*m;
		i1[1] = icoVertices[v1][1]+icoNormals[v1][1]*m;
		i1[2] = icoVertices[v1][2]+icoNormals[v1][2]*m;

		i2[0] = icoVertices[v2][0]+icoNormals[v2][0]*m;
		i2[1] = icoVertices[v2][1]+icoNormals[v2][1]*m;
		i2[2] = icoVertices[v2][2]+icoNormals[v2][2]*m;

		i3[0] = icoVertices[v3][0]+icoNormals[v3][0]*m;
		i3[1] = icoVertices[v3][1]+icoNormals[v3][1]*m;
		i3[2] = icoVertices[v3][2]+icoNormals[v3][2]*m;

		ls[2] += integral(icoVertices[v1], icoVertices[v2], icoVertices[v3],i1, i2, i3);
	}

	return 1;
}

void calculateBurgersVector(real nnDist, cell *central, int centralInd, cell **neiCells, int *neiIndices, int n, real lsNorm[3], real bv[3]) {
	real sqr = nnDist*nnDist;

	real bvInter[12][3];
	real x = ORT(central, centralInd, X);
	real y = ORT(central, centralInd, Y);
	real z = ORT(central, centralInd, Z);
	int i;

	bv[0] = 0.; bv[1] = 0.; bv[2] = 0.;

	for (i = 0; i<12; i++) {
		interpolateBurgersVector(x+icoVertices[i][0], y+icoVertices[i][1], z+icoVertices[i][2], sqr, neiCells, neiIndices, n, lsNorm, bvInter[i]);
	}
	/*Integrate the Nye tensor field on the surface of a icosahedron, as an approximation of a sphere, to get the Burgers vector*/
	for (i = 0; i<20; i++) {
		int v1 = icoFaces[i][0];
		int v2 = icoFaces[i][1];
		int v3 = icoFaces[i][2];

		real i1[3] = {icoVertices[v1][0]+icoNormals[v1][0]*bvInter[v1][0],
					  icoVertices[v1][1]+icoNormals[v1][1]*bvInter[v1][0],
					  icoVertices[v1][2]+icoNormals[v1][2]*bvInter[v1][0]};
		real i2[3] = {icoVertices[v2][0]+icoNormals[v2][0]*bvInter[v2][0],
					  icoVertices[v2][1]+icoNormals[v2][1]*bvInter[v2][0],
					  icoVertices[v2][2]+icoNormals[v2][2]*bvInter[v2][0]};
		real i3[3] = {icoVertices[v3][0]+icoNormals[v3][0]*bvInter[v3][0],
					  icoVertices[v3][1]+icoNormals[v3][1]*bvInter[v3][0],
					  icoVertices[v3][2]+icoNormals[v3][2]*bvInter[v3][0]};

		bv[0] += integral(icoVertices[v1], icoVertices[v2], icoVertices[v3], i1, i2, i3);

		i1[0] = icoVertices[v1][0]+icoNormals[v1][0]*bvInter[v1][1];
		i1[1] = icoVertices[v1][1]+icoNormals[v1][1]*bvInter[v1][1];
		i1[2] = icoVertices[v1][2]+icoNormals[v1][2]*bvInter[v1][1];

		i2[0] = icoVertices[v2][0]+icoNormals[v2][0]*bvInter[v2][1];
		i2[1] = icoVertices[v2][1]+icoNormals[v2][1]*bvInter[v2][1];
		i2[2] = icoVertices[v2][2]+icoNormals[v2][2]*bvInter[v2][1];

		i3[0] = icoVertices[v3][0]+icoNormals[v3][0]*bvInter[v3][1];
		i3[1] = icoVertices[v3][1]+icoNormals[v3][1]*bvInter[v3][1];
		i3[2] = icoVertices[v3][2]+icoNormals[v3][2]*bvInter[v3][1];

		bv[1] += integral(icoVertices[v1], icoVertices[v2], icoVertices[v3], i1, i2, i3);

		i1[0] = icoVertices[v1][0]+icoNormals[v1][0]*bvInter[v1][2];
		i1[1] = icoVertices[v1][1]+icoNormals[v1][1]*bvInter[v1][2];
		i1[2] = icoVertices[v1][2]+icoNormals[v1][2]*bvInter[v1][2];

		i2[0] = icoVertices[v2][0]+icoNormals[v2][0]*bvInter[v2][2];
		i2[1] = icoVertices[v2][1]+icoNormals[v2][1]*bvInter[v2][2];
		i2[2] = icoVertices[v2][2]+icoNormals[v2][2]*bvInter[v2][2];

		i3[0] = icoVertices[v3][0]+icoNormals[v3][0]*bvInter[v3][2];
		i3[1] = icoVertices[v3][1]+icoNormals[v3][1]*bvInter[v3][2];
		i3[2] = icoVertices[v3][2]+icoNormals[v3][2]*bvInter[v3][2];

		bv[2] += integral(icoVertices[v1], icoVertices[v2], icoVertices[v3], i1, i2, i3);
	}
}

/******************************************************************************
*  pack nyeTensorInfo from buffer cell into MPI buffer
******************************************************************************/
void pack_nyeTensorInfo(msgbuf *b, int k, int l, int m) {
	int i, j = b->n;
	minicell *from;

	from = PTR_3D_V(cell_array, k, l, m, cell_dim);

	for (i = 0; i < from->n; ++i) {
		if (NYE(from, i) == NULL){
			b->data[j++] = 0;
		} else {
			nyeTensorInfo *info = NYE(from, i);
			b->data[j++] = 1;
			b->data[j++] = info->lcm[0][0];
			b->data[j++] = info->lcm[0][1];
			b->data[j++] = info->lcm[0][2];
			b->data[j++] = info->lcm[1][0];
			b->data[j++] = info->lcm[1][1];
			b->data[j++] = info->lcm[1][2];
			b->data[j++] = info->lcm[2][0];
			b->data[j++] = info->lcm[2][1];
			b->data[j++] = info->lcm[2][2];

			b->data[j++] = info->nyeTensor[0][0];
			b->data[j++] = info->nyeTensor[0][1];
			b->data[j++] = info->nyeTensor[0][2];
			b->data[j++] = info->nyeTensor[1][0];
			b->data[j++] = info->nyeTensor[1][1];
			b->data[j++] = info->nyeTensor[1][2];
			b->data[j++] = info->nyeTensor[2][0];
			b->data[j++] = info->nyeTensor[2][1];
			b->data[j++] = info->nyeTensor[2][2];
		}
		//LineSense and resultant burgers vector do not need to be copied
	}

	b->n = j;
	if (b->n_max < b->n)
		error("Buffer overflow in pack_nyeTensorInfo - increase msgbuf_size");
}

/******************************************************************************
*  unpack nyeTensorInfo from MPI buffer, and add them to those of the original cell
******************************************************************************/

void unpack_nyeTensorInfo(msgbuf *b, int k, int l, int m) {
	int i, isNull, j = b->n;
	minicell *to;

	to = PTR_3D_V(cell_array, k, l, m, cell_dim);

	for (i = 0; i < to->n; ++i) {
		isNull = b->data[j++];
		if (isNull){
			/*
			 * Buffer cells are only allocated on demand
			 */
			if (NYE(to,i) == NULL){
				NYE(to,i) = alloc_nyeTensorInfo();
			}
			nyeTensorInfo *info = NYE(to, i);

			info->lcm[0][0] = b->data[j++];
			info->lcm[0][1] = b->data[j++];
			info->lcm[0][2] = b->data[j++];
			info->lcm[1][0] = b->data[j++];
			info->lcm[1][1] = b->data[j++];
			info->lcm[1][2] = b->data[j++];
			info->lcm[2][0] = b->data[j++];
			info->lcm[2][1] = b->data[j++];
			info->lcm[2][2] = b->data[j++];

			info->nyeTensor[0][0] = b->data[j++];
			info->nyeTensor[0][1] = b->data[j++];
			info->nyeTensor[0][2] = b->data[j++];
			info->nyeTensor[1][0] = b->data[j++];
			info->nyeTensor[1][1] = b->data[j++];
			info->nyeTensor[1][2] = b->data[j++];
			info->nyeTensor[2][0] = b->data[j++];
			info->nyeTensor[2][1] = b->data[j++];
			info->nyeTensor[2][2] = b->data[j++];
			//LineSense and resultant burgers vector do not need to be copied
		}
	}

	b->n = j;
	if (b->n_max < b->n)
		error("Buffer overflow in unpack_nyeTensorInfo - increase msgbuf_size");
}

void copy_nyeTensorInfo(int k, int l, int m, int r, int s, int t) {
	int i;
	minicell *from, *to;

	from = PTR_3D_V(cell_array, k, l, m, cell_dim);
	to = PTR_3D_V(cell_array, r, s, t, cell_dim);

	for (i = 0; i < to->n; ++i) {
		if (from->nyeTens[i] != NULL) {
			if (to->nyeTens[i] == NULL) {
				to->nyeTens[i] = alloc_nyeTensorInfo();
			}

			to->nyeTens[i]->lcm[0][0] = from->nyeTens[i]->lcm[0][0];
			to->nyeTens[i]->lcm[0][1] = from->nyeTens[i]->lcm[0][1];
			to->nyeTens[i]->lcm[0][2] = from->nyeTens[i]->lcm[0][2];
			to->nyeTens[i]->lcm[1][0] = from->nyeTens[i]->lcm[1][0];
			to->nyeTens[i]->lcm[1][1] = from->nyeTens[i]->lcm[1][1];
			to->nyeTens[i]->lcm[1][2] = from->nyeTens[i]->lcm[1][2];
			to->nyeTens[i]->lcm[2][0] = from->nyeTens[i]->lcm[2][0];
			to->nyeTens[i]->lcm[2][1] = from->nyeTens[i]->lcm[2][1];
			to->nyeTens[i]->lcm[2][2] = from->nyeTens[i]->lcm[2][2];

			to->nyeTens[i]->nyeTensor[0][0] = from->nyeTens[i]->nyeTensor[0][0];
			to->nyeTens[i]->nyeTensor[0][1] = from->nyeTens[i]->nyeTensor[0][1];
			to->nyeTens[i]->nyeTensor[0][2] = from->nyeTens[i]->nyeTensor[0][2];
			to->nyeTens[i]->nyeTensor[1][0] = from->nyeTens[i]->nyeTensor[1][0];
			to->nyeTens[i]->nyeTensor[1][1] = from->nyeTens[i]->nyeTensor[1][1];
			to->nyeTens[i]->nyeTensor[1][2] = from->nyeTens[i]->nyeTensor[1][2];
			to->nyeTens[i]->nyeTensor[2][0] = from->nyeTens[i]->nyeTensor[2][0];
			to->nyeTens[i]->nyeTensor[2][1] = from->nyeTens[i]->nyeTensor[2][1];
			to->nyeTens[i]->nyeTensor[2][2] = from->nyeTens[i]->nyeTensor[2][2];
		}
	}
}

void calculateNyeTensorData(){
	int k,i,j,ii,l,iii,m;

	cell *neiCells[100];
	int neighIndices[100];
	int numNeighs;
	int ok, nneigh1, nneigh2;
	int dir1[3];
	int dir2[3];
	int dir3[3];

	/* Three orthogonal assumed burgers vectors */
	dir1[0] = 1; dir1[1] = 0; dir1[2] = 0;
	dir2[0] = 0; dir2[1] = 1; dir2[2] = 0;
	dir3[0] = 0; dir3[1] = 0; dir3[2] = 1;

	/* Compute neighbor table (if already computed, nothing is done)*/
	do_neightab_complete();

	/*
	 * Allocate memory and calculate lattice correspondence matrices (LCM)
	 */
	for (k = 0; k < ncells; k++) {
		cell *p = CELLPTR(k);

		for (i = 0; i < p->n; i++) {
			if (HOPSTODEFECT(p,i) <= 3)
			{
				p->nyeTens[i] = alloc_nyeTensorInfo();
				calculateLcm(p,i);
			}
		}
	}

	/*Exchange LCM between cells*/
	send_fromCellsToBuffer(copy_nyeTensorInfo,pack_nyeTensorInfo,unpack_nyeTensorInfo);

	/**
	 * Calculate the nye tensor
	 */
	for (k = 0; k < ncells; k++) {
		cell *p = CELLPTR(k);
		for (i = 0; i < p->n; i++) {
			if (HOPSTODEFECT(p, i) <= 2)
			{
				calculateNye(p, i);
			}
		}
	}

	/*Exchange nye tensor between cells*/
	send_fromCellsToBuffer(copy_nyeTensorInfo,pack_nyeTensorInfo,unpack_nyeTensorInfo);

	/*Calculate lineSense and burgers vectors for defects*/
	for (k = 0; k < ncells; k++) {
		cell *p = CELLPTR(k);
		for (i = 0; i < p->n; i++) {
			if ( ADATYPE(p, i) != ada_default_type &&
					((ada_crystal_structure == ADA_FCC_CONFIG && ADATYPE(p,i)>=3 && ADATYPE(p,i)<=5) ||
					 (ada_crystal_structure == ADA_BCC_CONFIG && ADATYPE(p,i)!=6)) )
				{
				nyeTensorInfo *info = NYE(p,i);
				/*Copy the central atom, its nearest neighbors and their nearest neighbors into a list*/
				numNeighs = 0;
				/*Add the central atom*/
				neiCells[numNeighs] = p;
				neighIndices[numNeighs] = i;
				numNeighs++;

				/*Add its nearest neighbors*/
				nneigh1 = NEIGH(p, i)->n;
				for (j = 0; j < nneigh1; j++) {
					cell *q = NEIGH(p, i)->cl[j];
					ii = NEIGH(p, i)->num[j];
					neiCells[numNeighs] = q;
					neighIndices[numNeighs] = ii;
					numNeighs++;
				}

				/*Add its second nearest neighbors*/
				nneigh1 = NEIGH(p, i)->n;
				for (j = 0; j < nneigh1; j++) {
					cell *q = NEIGH(p, i)->cl[j];
					ii = NEIGH(p, i)->num[j];
					nneigh2 = NEIGH(q, ii)->n;
					for (l = 0; l < nneigh2; l++) {
						cell *r = NEIGH(q, ii)->cl[l];
						iii = NEIGH(q, ii)->num[l];
						/*Test if this atom is already inserted in the list*/
						ok = 1;
						for (m=0; m<numNeighs;m++){
							if (neiCells[m] == r && neighIndices[m] == iii){
								ok = 0;
								break;
							}
						}
						/*If not insert it too*/
						if (ok){
							neiCells[numNeighs] = r;
							neighIndices[numNeighs] = iii;
							numNeighs++;
						}
					}
				}

				ok = calculateLineSense(integrationRadius,p,i,neiCells,neighIndices, numNeighs, info->ls, dir1);
				if (ok){
					real l = SQRT(info->ls[0]*info->ls[0] + info->ls[1]*info->ls[1] + info->ls[2]*info->ls[2]);
					/*If the line sense vector is very small, test another assumed burgers vector*/
					if (l<0.316) {
						calculateLineSense(integrationRadius,p,i,neiCells,neighIndices, numNeighs, info->ls, dir2);
						l = SQRT(info->ls[0]*info->ls[0] + info->ls[1]*info->ls[1] + info->ls[2]*info->ls[2]);
						/*And finally the last of the three*/
						if (l<0.316) {
							calculateLineSense(integrationRadius,p,i,neiCells,neighIndices, numNeighs, info->ls, dir3);
							l = SQRT(info->ls[0]*info->ls[0] + info->ls[1]*info->ls[1] + info->ls[2]*info->ls[2]);
						}
					}

					if (l>0.){
						info->ls[0] /= l; info->ls[1] /= l; info->ls[2] /= l;
						/*Calculate burgers vector*/
						calculateBurgersVector(integrationRadius,p,i,neiCells,neighIndices,numNeighs,info->ls,info->bv);
					}
				}
			}
		}
	}

	/*Remove rbv with an insignificant burgers vector*/
	for (k = 0; k < ncells; k++) {
			cell *p = CELLPTR(k);

			for (i = 0; i < p->n; i++) {
				nyeTensorInfo *info = NYE(p,i);
				if (info != NULL){
					real l = info->bv[0]*info->bv[0] + info->bv[1]*info->bv[1] + info->bv[2]*info->bv[2];
					if (l < 0.02*integrationRadius*integrationRadius){
						free(p->nyeTens[i]);
						p->nyeTens[i] = NULL;
					}
				}
			}
		}
}
