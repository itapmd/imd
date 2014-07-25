/**
 * imd_loadBalance_direct -- primary routines for dynamic load-balancing in IMD
 * Communication routines for load balancing using direct communication between CPUs
 */
#include "imd.h"

void lb_initDirect() {
	int i, x, y, z;

	if (myid == 0){
		if (lb_balancingType == 1)
			printf("LOAD BALANCING: Communication with any CPU enabled by \"lb_balancingType 1\"\n");
		else if (lb_balancingType == 2)
			printf("LOAD BALANCING: Balancing using orthogonal domains \"lb_balancingType 2\"\n");
		else printf("LOAD BALANCING: Communication limited to direct neighbors by \"lb_balancingType 0\"\n");
	}

	lb_halfspaceLUT = malloc(num_cpus * sizeof *lb_halfspaceLUT);
	lb_pbcFlag = malloc(num_cpus * sizeof *lb_pbcFlag);
	if (lb_pbcFlag == NULL || lb_halfspaceLUT == NULL)
		error("cannot allocate lookup tables in lb_initDirect");

	for (i = 0; i < num_cpus; i++) {
		lb_halfspaceLUT[i] = -1;
		lb_pbcFlag[i].x = 0;
		lb_pbcFlag[i].y = 0;
		lb_pbcFlag[i].z = 0;
	}

	/* Divide the geometry of cpus into halfspaces, into one halfspace forces will be send */
	/* The other will send cells, required to properly setup send buffers and in make_cell_lists*/
	lb_halfspaceLUT[nbuen  ] = LB_SEND_CELL;    lb_halfspaceLUT[nbdown ] = LB_SEND_FORCE;
	lb_halfspaceLUT[nben   ] = LB_SEND_CELL;    lb_halfspaceLUT[nbus   ] = LB_SEND_FORCE;
	lb_halfspaceLUT[nbdne  ] = LB_SEND_CELL;    lb_halfspaceLUT[nbsouth] = LB_SEND_FORCE;
	lb_halfspaceLUT[nbue   ] = LB_SEND_CELL;    lb_halfspaceLUT[nbds   ] = LB_SEND_FORCE;
	lb_halfspaceLUT[nbeast ] = LB_SEND_CELL;    lb_halfspaceLUT[nbunw  ] = LB_SEND_FORCE;
	lb_halfspaceLUT[nbde   ] = LB_SEND_CELL;    lb_halfspaceLUT[nbnw   ] = LB_SEND_FORCE;
    lb_halfspaceLUT[nbuse  ] = LB_SEND_CELL;    lb_halfspaceLUT[nbdwn  ] = LB_SEND_FORCE;
    lb_halfspaceLUT[nbse   ] = LB_SEND_CELL;    lb_halfspaceLUT[nbuw   ] = LB_SEND_FORCE;
	lb_halfspaceLUT[nbdes  ] = LB_SEND_CELL;    lb_halfspaceLUT[nbwest ] = LB_SEND_FORCE;
	lb_halfspaceLUT[nbun   ] = LB_SEND_CELL;    lb_halfspaceLUT[nbdw   ] = LB_SEND_FORCE;
	lb_halfspaceLUT[nbnorth] = LB_SEND_CELL;    lb_halfspaceLUT[nbuws  ] = LB_SEND_FORCE;
	lb_halfspaceLUT[nbdn   ] = LB_SEND_CELL;    lb_halfspaceLUT[nbws   ] = LB_SEND_FORCE;
	lb_halfspaceLUT[nbup   ] = LB_SEND_CELL;    lb_halfspaceLUT[nbdsw  ] = LB_SEND_FORCE;
	lb_halfspaceLUT[myid   ] = LB_DONT_SEND;

	/*Override the halfspaces if only 2 cpus exist in a direction*/
	if (cpu_dim.x == 2){
		/* Eastern neighbors                   Western neighbors */
		lb_halfspaceLUT[nbuen  ] = -1;   lb_halfspaceLUT[nbunw  ] = -1;
		lb_halfspaceLUT[nben   ] = -1;   lb_halfspaceLUT[nbnw   ] = -1;
		lb_halfspaceLUT[nbdne  ] = -1;   lb_halfspaceLUT[nbdwn  ] = -1;
		lb_halfspaceLUT[nbue   ] = -1;   lb_halfspaceLUT[nbuw   ] = -1;
		lb_halfspaceLUT[nbeast ] = -1;   lb_halfspaceLUT[nbwest ] = -1;
		lb_halfspaceLUT[nbde   ] = -1;   lb_halfspaceLUT[nbdw   ] = -1;
		lb_halfspaceLUT[nbuse  ] = -1;   lb_halfspaceLUT[nbuws  ] = -1;
		lb_halfspaceLUT[nbse   ] = -1;   lb_halfspaceLUT[nbws   ] = -1;
		lb_halfspaceLUT[nbdes  ] = -1;   lb_halfspaceLUT[nbdsw  ] = -1;
	}

	if (cpu_dim.y == 2){
		/* Northern neighbors               Southern neighbors */
		lb_halfspaceLUT[nbuen  ] = -1;	lb_halfspaceLUT[nbuse  ] = -1;
		lb_halfspaceLUT[nben   ] = -1;	lb_halfspaceLUT[nbse   ] = -1;
		lb_halfspaceLUT[nbdne  ] = -1;	lb_halfspaceLUT[nbdes  ] = -1;
		lb_halfspaceLUT[nbun   ] = -1;	lb_halfspaceLUT[nbus   ] = -1;
		lb_halfspaceLUT[nbnorth] = -1;	lb_halfspaceLUT[nbsouth] = -1;
		lb_halfspaceLUT[nbdn   ] = -1;	lb_halfspaceLUT[nbds   ] = -1;
		lb_halfspaceLUT[nbunw  ] = -1;	lb_halfspaceLUT[nbuws  ] = -1;
		lb_halfspaceLUT[nbnw   ] = -1;	lb_halfspaceLUT[nbws   ] = -1;
		lb_halfspaceLUT[nbdwn  ] = -1;	lb_halfspaceLUT[nbdsw  ] = -1;
	}

	if (cpu_dim.z == 2){
		/*    neighbors up                neighbors down*/
		lb_halfspaceLUT[nbuen  ] = -1;	lb_halfspaceLUT[nbdne  ] = -1;
		lb_halfspaceLUT[nbue   ] = -1;	lb_halfspaceLUT[nbde   ] = -1;
		lb_halfspaceLUT[nbuse  ] = -1;	lb_halfspaceLUT[nbdes  ] = -1;
		lb_halfspaceLUT[nbun   ] = -1;	lb_halfspaceLUT[nbdn   ] = -1;
		lb_halfspaceLUT[nbup   ] = -1;	lb_halfspaceLUT[nbdown ] = -1;
		lb_halfspaceLUT[nbus   ] = -1;	lb_halfspaceLUT[nbds   ] = -1;
		lb_halfspaceLUT[nbunw  ] = -1;	lb_halfspaceLUT[nbdwn  ] = -1;
		lb_halfspaceLUT[nbuw   ] = -1;	lb_halfspaceLUT[nbdw   ] = -1;
		lb_halfspaceLUT[nbuws  ] = -1;	lb_halfspaceLUT[nbdsw  ] = -1;
	}

	/* send cells to one half of non-neighbor domains */
	/* send forces to the other half, receiving domain is inverse*/
	for (i = 0; i < num_cpus; i++) {
		if (lb_halfspaceLUT[i] == -1){
			if (myid % 2 == i % 2) {
				if (myid<i)
					lb_halfspaceLUT[i] = LB_SEND_FORCE;
				else lb_halfspaceLUT[i] = LB_SEND_CELL;
			} else {
				if (myid<i)
					lb_halfspaceLUT[i] = LB_SEND_CELL;
				else lb_halfspaceLUT[i] = LB_SEND_FORCE;
			}
		}
	}

	/* setup the pbc-correction vectors for domains at boundaries*/
	ivektor v;
	if (my_coord.x == 0 && pbc_dirs.x == 1) {
		for (y = 0; y < cpu_dim.y; y++)
			for (z = 0; z < cpu_dim.z; z++){
				v.x = cpu_dim.x-1; v.y = y; v.z = z;
				lb_pbcFlag[cpu_grid_coord(v)].x = 1;
			}
	}
	if (my_coord.x == cpu_dim.x - 1 && pbc_dirs.x == 1) {
		for (y = 0; y < cpu_dim.y; y++)
			for (z = 0; z < cpu_dim.z; z++){
				v.x = 0; v.y = y; v.z = z;
				lb_pbcFlag[cpu_grid_coord(v)].x = 1;
			}
	}

	if (my_coord.y == 0 && pbc_dirs.y == 1) {
		for (x = 0; x < cpu_dim.x; x++)
			for (z = 0; z < cpu_dim.z; z++){
				v.x = x; v.y = cpu_dim.y-1; v.z = z;
				lb_pbcFlag[cpu_grid_coord(v)].y = 1;
			}
	}
	if (my_coord.y == cpu_dim.y - 1 && pbc_dirs.y == 1) {
		for (x = 0; x < cpu_dim.x; x++)
			for (z = 0; z < cpu_dim.z; z++){
				v.x = x; v.y = 0; v.z = z;
				lb_pbcFlag[cpu_grid_coord(v)].y = 1;
			}
	}

	if (my_coord.z == 0 && pbc_dirs.z == 1) {
		for (x = 0; x < cpu_dim.x; x++)
			for (y = 0; y < cpu_dim.y; y++){
				v.x = x; v.y = y; v.z = cpu_dim.z - 1;
				lb_pbcFlag[cpu_grid_coord(v)].z = 1;
			}
	}
	if (my_coord.z == cpu_dim.z - 1 && pbc_dirs.z == 1) {
		for (x = 0; x < cpu_dim.x; x++)
			for (y = 0; y < cpu_dim.y; y++){
				v.x = x; v.y = y; v.z = 0;
				lb_pbcFlag[cpu_grid_coord(v)].z = 1;
			}
	}
}

int lb_syncBufferCellAffinity() {
	int i,j, valid;
	int x,y,z;
	int x2,y2,z2;
	cell *cell, *cell2;
	lb_domainInfo neighborDomains[26];

	real *allCorners;
	real *corners;
	lb_domainInfo *allDomains;

	valid = 1;

	for (i = 0; i < lb_nTotalComms; i++){
		/* free old send lists */
		free(lb_sendForces[i]);
		free(lb_sendCells[i]);
	}

	/* The number of cells for which information is to be exchanged is counted in this array*/
	int *numBufferCellsToCPU = malloc(num_cpus * sizeof *numBufferCellsToCPU);
	if (numBufferCellsToCPU == NULL)
			error("Cannot allocate array for LookUpTables in syncBufferCellAffinity");
	/* Look-up table to get the index which buffer has be be used to communicate */
	/* with another CPU. Since the number of partners can vary, the number of buffers */
	/* changes as well. */
	int *lutCommIndex = malloc(num_cpus * sizeof *lutCommIndex);

	if (lutCommIndex == NULL)
			error("Cannot allocate array for LookUpTables in syncBufferCellAffinity");

	for (i = 0; i < num_cpus; i++) {
		numBufferCellsToCPU[i] = 0;
		lutCommIndex[i] = -1;
	}

	int listNeighbors[26];
	listNeighbors[0] = nbuen;    listNeighbors[1] = nben;
	listNeighbors[2] = nbdne;    listNeighbors[3] = nbue;
	listNeighbors[4] = nbeast;   listNeighbors[5] = nbde;
    listNeighbors[6] = nbuse;    listNeighbors[7] = nbse;
	listNeighbors[8] = nbdes;    listNeighbors[9] = nbun;
	listNeighbors[10] = nbnorth; listNeighbors[11] = nbdn;
	listNeighbors[12] = nbup;    listNeighbors[13] = nbdown;
	listNeighbors[14] = nbus;    listNeighbors[15] = nbsouth;
	listNeighbors[16] = nbds;    listNeighbors[17] = nbunw;
	listNeighbors[18] = nbnw;    listNeighbors[19] = nbdwn;
	listNeighbors[20] = nbuw;    listNeighbors[21] = nbwest;
	listNeighbors[22] = nbdw;    listNeighbors[23] = nbuws;
	listNeighbors[24] = nbws;    listNeighbors[25] = nbdsw;

	allDomains = malloc(num_cpus * sizeof *allDomains);

#ifdef MPI2
	MPI_Alloc_mem(8*3*sizeof *corners, MPI_INFO_NULL, &corners);
	MPI_Alloc_mem(8*3*num_cpus * sizeof *allCorners, MPI_INFO_NULL, &allCorners);
#else
	corners = malloc(8*3*sizeof *corners);
	allCallCorners = malloc(8*3*num_cpus * sizeof *allCorners);
#endif

	for (i = 0; i<8; i++){
		corners[i*3+0] = lb_domain.corners[i].p.x;
		corners[i*3+1] = lb_domain.corners[i].p.y;
		corners[i*3+2] = lb_domain.corners[i].p.z;
	}

	/*
	 * Distribute and gather all domains
	 */
	MPI_Allgather(corners, 24, REAL, allCorners, 24, REAL, MPI_COMM_WORLD);
	//Create neighbor domains
	for (i=0; i<26;i++){
		for (j = 0; j < 8; j++) {
			neighborDomains[i].corners[j].p.x = allCorners[j*3+0 + listNeighbors[i]*24];
			neighborDomains[i].corners[j].p.y = allCorners[j*3+1 + listNeighbors[i]*24];
			neighborDomains[i].corners[j].p.z = allCorners[j*3+2 + listNeighbors[i]*24];
		}
		lb_updateDomain(&neighborDomains[i]);
	}
	//Create neighbor domains of all CPUs
	for (i = 0; i < num_cpus; i++) {
		for (j = 0; j<8; j++){
			allDomains[i].corners[j].p.x = allCorners[j*3+0 + i*24];
			allDomains[i].corners[j].p.y = allCorners[j*3+1 + i*24];
			allDomains[i].corners[j].p.z = allCorners[j*3+2 + i*24];
		}
		lb_updateDomain(&allDomains[i]);
	}

	/* Assign buffer cells affinity to the domains holding the corresponding real cell*/
	for (x = 0; x < cell_dim.x; ++x) {
		for (y = 0; y < cell_dim.y; ++y) {
			for (z = 0; z < cell_dim.z; ++z) {
				cell = PTR_3D_V(cell_array, x, y, z, cell_dim);
				if (cell->lb_cell_type == LB_BUFFER_CELL || cell->lb_cell_type == LB_EMPTY_CELL) {
					int indexX = x + lb_cell_offset.x;
					int indexY = y + lb_cell_offset.y;
					int indexZ = z + lb_cell_offset.z;

					/* Wrap around PBCs*/
					ivektor3d pbcWrap = {0, 0, 0};
					if (indexX < 0) {indexX += global_cell_dim.x; pbcWrap.x = 1;}
					if (indexY < 0) {indexY += global_cell_dim.y; pbcWrap.y = 1;}
					if (indexZ < 0) {indexZ += global_cell_dim.z; pbcWrap.z = 1;}
					if (indexX >= global_cell_dim.x) {indexX -= global_cell_dim.x; pbcWrap.x = 1;}
					if (indexY >= global_cell_dim.y) {indexY -= global_cell_dim.y; pbcWrap.y = 1;}
					if (indexZ >= global_cell_dim.z) {indexZ -= global_cell_dim.z; pbcWrap.z = 1;}

					vektor center = lb_getCellCenter(indexX, indexY, indexZ);
					int cellAssigned = 0;
					for (j = 0; j < 26; j++) { /*Test first the 26 nearest neighbors */
						if (lb_isPointInDomain(center, &neighborDomains[j])) {
							cell->lb_cpu_affinity = listNeighbors[j];
							cellAssigned = 1;
							break;
						}
					}
					/* None of the 26 accepted, test all cpus if they have the cell*/
					/* If communication is restricted to 26 neighbors, reject step here*/
					if (!cellAssigned) {
						if (lb_balancingType!=1){
							if(cell->lb_cell_type == LB_EMPTY_CELL) cellAssigned = 1;
							if(cell->lb_cell_type == LB_BUFFER_CELL) valid = 0;
						}
						else {
							for (j = 0; j < num_cpus; j++) {
								if (lb_isPointInDomain(center, &allDomains[j])) {
									cell->lb_cpu_affinity = j;
									cellAssigned = 1;
									lb_pbcFlag[j] = pbcWrap;
									break;
								}
							}
						}
					}
					if (!cellAssigned){
#ifdef DEBUG
						printf("LB: Cell could not be assigned to any cpu on %i\n", myid);
#endif
						valid = 0; /* causing the loadbalancing step to be reverted*/
						/* some values need to be set to avoid crashes */
						cell->lb_cpu_affinity = myid;
						cell->lb_cell_type = LB_REAL_CELL;
					}
				}
			}
		}
	}

#ifdef MPI2
	MPI_Free_mem(corners);
	MPI_Free_mem(allCorners);
#else
	free(corners);
	free(allCorners);
#endif

	free(allDomains);

	/* Creation of communication lists */

	lb_nTotalComms = 0;
	/* Count how many cells / forces have to be send to different CPUs*/

	/* Identify to which cpus communication is required, look which buffer cells exist */
	/* First identify which CPUs must receive cells */
	for (x = 0; x < cell_dim.x; ++x) {
		for (y = 0; y < cell_dim.y; ++y) {
			for (z = 0; z < cell_dim.z; ++z) {
				cell = PTR_3D_V(cell_array, x, y, z, cell_dim);
				int index = cell->lb_cpu_affinity;
				if (cell->lb_cell_type == LB_BUFFER_CELL && lb_halfspaceLUT[index] == LB_SEND_CELL) {
					numBufferCellsToCPU[index]++;
					if (lutCommIndex[index] == -1)	/* Keep the index to communicate*/
						lutCommIndex[index] = lb_nTotalComms++;
					cell->lb_neighbor_index = lutCommIndex[index];
				}
			}
		}
	}

	lb_nForceComms = 0;
	/* Then identify which CPUs receive forces*/
	for (x = 0; x < cell_dim.x; ++x) {
		for (y = 0; y < cell_dim.y; ++y) {
			for (z = 0; z < cell_dim.z; ++z) {
				cell = PTR_3D_V(cell_array, x, y, z, cell_dim);
				int index = cell->lb_cpu_affinity;
				if (cell->lb_cell_type == LB_BUFFER_CELL && lb_halfspaceLUT[index] == LB_SEND_FORCE) {
					numBufferCellsToCPU[index]++;
					if (lutCommIndex[index] == -1){	/* Keep the index to communicate*/
						lutCommIndex[index] = lb_nTotalComms++;
						lb_nForceComms++;	/* Count the number of unique partners to send forces */
					}
					cell->lb_neighbor_index = lutCommIndex[index];
				}
			}
		}
	}

	/* lb_nTotalComms now contains the number of communications required to communicate with all partners */
	/* the set is split into two parts for half-space communication*/
	/* the first half of neighbors must receive cells, the second half forces*/
	/* Additionally, the number of buffer cell send to different neighbors is stored in lutNumSendForcesToCPU[]*/

	/* Allocate the index buffers for forces*/
	lb_commIndexToCpu = realloc(lb_commIndexToCpu, lb_nTotalComms * sizeof(int));
	lb_nSendForces = realloc(lb_nSendForces, lb_nTotalComms * sizeof(int));
	for (i = 0; i < num_cpus; i++) {
		if (lutCommIndex[i] != -1){
			lb_nSendForces[lutCommIndex[i]] = numBufferCellsToCPU[i];
			//The lists are send to the CPUs whose ranks are stored in lb_commIndexToCPU
			lb_commIndexToCpu[lutCommIndex[i]] = i;
		}
	}

	int *elementsInList = malloc(lb_nTotalComms * sizeof *elementsInList);
	lb_sendForces = realloc(lb_sendForces, lb_nTotalComms * sizeof *lb_sendForces);

	for (i = 0; i < lb_nTotalComms; i++){
		elementsInList[i] = 0;
		//Allocate the send list with proper sizes
		lb_sendForces[i] = malloc(lb_nSendForces[i] * sizeof *lb_sendForces[i]);
	}

	/*Now fill the send lists with the indices of cells that must be send for communication of forces*/
	for (x = 0; x < cell_dim.x; ++x) {
		for (y = 0; y < cell_dim.y; ++y) {
			for (z = 0; z < cell_dim.z; ++z) {
				cell = PTR_3D_V(cell_array, x, y, z, cell_dim);
				if (cell->lb_cell_type == LB_BUFFER_CELL) {
					int index = cell->lb_neighbor_index;
					i = elementsInList[index];

					lb_sendForces[index][i].x = x;
					lb_sendForces[index][i].y = y;
					lb_sendForces[index][i].z = z;
					elementsInList[index]++;
				}
			}
		}
	}

	free(numBufferCellsToCPU);

	/* Communication for forces is now ready, create lists for communications of cells as well */
	int *targeted = malloc(lb_nTotalComms * sizeof *targeted);
	lb_nSendCells = realloc(lb_nSendCells, lb_nTotalComms * sizeof *lb_nSendCells);

	for (i = 0; i < lb_nTotalComms; i++){
		targeted[i] = 0;
		lb_nSendCells[i] = 0;
	}

	for (x = 1; x < cell_dim.x-1; ++x) {
		for (y = 1; y < cell_dim.y-1; ++y) {
			for (z = 1; z < cell_dim.z-1; ++z) {
				cell = PTR_3D_V(cell_array, x, y, z, cell_dim);
				if (cell->lb_cell_type == LB_REAL_CELL) {
					/* Check all neighbors of the cell if they are next to a buffer */
					/* At first clear buffer indicating if a neighbor has been found */
					for (i = 0; i < lb_nTotalComms; i++)
						targeted[i] = 0;

					for (x2 = -1; x2 <= 1; ++x2) {
						for (y2 = -1; y2 <= 1; ++y2) {
							for (z2 = -1; z2 <= 1; ++z2) {
								cell2 = PTR_3D_V(cell_array, x+x2, y+y2, z+z2, cell_dim);
								if (cell2->lb_cell_type == LB_BUFFER_CELL){
									int index = cell2->lb_neighbor_index;
									/* Ensure that cell is only included only once in a list*/
									if (!targeted[index]){
										targeted[index] = 1;
										lb_nSendCells[index]++;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	/*The size of each cell list is known, allocate the lists */
	lb_sendCells = realloc(lb_sendCells, lb_nTotalComms * sizeof *lb_sendCells);
	for (i = 0; i < lb_nTotalComms; i++){
		elementsInList[i] = 0;
		//Allocate the send list with proper sizes
		lb_sendCells[i] = malloc(lb_nSendCells[i] * sizeof *lb_sendCells[i]);
	}

	/*Fill the send lists */
	for (x = 1; x < cell_dim.x-1; ++x) {
		for (y = 1; y < cell_dim.y-1; ++y) {
			for (z = 1; z < cell_dim.z-1; ++z) {
				cell = PTR_3D_V(cell_array, x, y, z, cell_dim);
				if (cell->lb_cell_type == LB_REAL_CELL) {
					/* Check all neighbors of the cell if they are next to a buffer */
					/* At first clear buffer indicating if a neighbor has been found */
					for (i = 0; i < lb_nTotalComms; i++)
						targeted[i] = 0;

					for (x2 = -1; x2 <= 1; ++x2) {
						for (y2 = -1; y2 <= 1; ++y2) {
							for (z2 = -1; z2 <= 1; ++z2) {
								cell2 = PTR_3D_V(cell_array, x+x2, y+y2, z+z2, cell_dim);
								if (cell2->lb_cell_type == LB_BUFFER_CELL){
									int index = cell2->lb_neighbor_index;
									if (!targeted[index]){
										targeted[index] = 1;
										i = elementsInList[index];
										lb_sendCells[index][i].x = x;
										lb_sendCells[index][i].y = y;
										lb_sendCells[index][i].z = z;
										elementsInList[index]++;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	free(elementsInList);
	free(targeted);

#ifdef DEBUG
	for (x=0; x<cell_dim.x; ++x){
		for (y=0; y<cell_dim.y; ++y){
			for (z=0; z<cell_dim.z; ++z){
				cell = PTR_3D_V(cell_array, x, y, z, cell_dim);
				if (cell->lb_cell_type == LB_BUFFER_CELL && cell->lb_cpu_affinity<0)
					error("Broken lb config");
				if (cell->lb_cell_type == LB_BUFFER_CELL && cell->lb_neighbor_index<0)
					error("Broken lb config");
			}
		}
	}
#endif

	return valid;
}

/******************************************************************************
 *
 *  Send cells into buffer cells on neighboring CPUs for the force computation.
 *  What exactly is sent is determined by the parameter functions.
 *
 ******************************************************************************/

void send_cells(void (*copy_func)(int, int, int, int, int, int, vektor),
		void (*pack_func)(msgbuf*, int, int, int, vektor),
		void (*unpack_func)(msgbuf*, int, int, int)) {
	sync_cells_direct(*copy_func, *pack_func, *unpack_func, 0);
}

/******************************************************************************
 *
 *  Add forces in buffer cells on neighboring CPUs back to the original cells.
 *  What exactly is sent is determined by the parameter functions.
 *
 ******************************************************************************/

void send_forces(void (*add_func)   (int, int, int, int, int, int),
                 void (*pack_func)  (msgbuf*, int, int, int),
                 void (*unpack_func)(msgbuf*, int, int, int)) {
	int i,k;

	int sendForces = lb_nForceComms;
	int recvForces = lb_nTotalComms - lb_nForceComms;

	MPI_Status stat;

	empty_mpi_buffers();
	int offset = lb_nTotalComms - lb_nForceComms;

	for (i = 0; i<sendForces;++i){
		/*Send data away*/
		k= i+offset;
		lb_copyForcesDataToSend(&lb_send_buf[k], lb_sendForces[k], lb_nSendForces[k], pack_func);
		isend_buf(&lb_send_buf[k], lb_commIndexToCpu[k], &lb_req_send[k]);
		lb_requests[k] = lb_req_send[k];
		lb_request_indices[k] = -1; /* Indicates no processing required */
	}

	for (i = 0; i<recvForces;++i){
		/*Start receiving data*/
		irecv_buf(&lb_recv_buf[i], lb_commIndexToCpu[i], &lb_req_recv[i]);
		lb_requests[i] = lb_req_recv[i];
		lb_request_indices[i] = i;
	}

	/*Receive and process data as soon as something is available*/
	for (i = lb_nTotalComms; i>0; i--){
		int finished;
		MPI_Waitany(i, lb_requests, &finished, &stat);
		int ind = lb_request_indices[finished];
		if (ind != -1){
			MPI_Get_count(&stat, REAL, &lb_recv_buf[ind].n);
			lb_unpackForcesDataFromBuffer(&lb_recv_buf[ind], (*unpack_func));
		}
		lb_requests[finished] = lb_requests[i-1];
		lb_request_indices[finished] = lb_request_indices[i-1];
	}
}

void sync_cells(void (*copy_func)(int, int, int, int, int, int, vektor),
				void (*pack_func)(msgbuf*, int, int, int, vektor),
				void (*unpack_func)(msgbuf*, int, int, int)) {
	sync_cells_direct((*copy_func), (*pack_func), (*unpack_func), 1);
}

/******************************************************************************
 *
 *  Synchronize cells into buffer cells on neighboring CPUs.
 *  What exactly is sent is determined by the parameter functions.
 *  Uses direct communication.
 *  all==0 Send cells only in the halfspace opposite to force communication
 *  all!=0 Send cells in all directions
 *
 ******************************************************************************/

void sync_cells_direct(void (*copy_func)(int, int, int, int, int, int, vektor),
				void (*pack_func)(msgbuf*, int, int, int, vektor),
				void (*unpack_func)(msgbuf*, int, int, int), int all) {
	int i,k;

	int sendCells;
	int recvCells;
	int totalOperations;

	if (all){
		sendCells = lb_nTotalComms;
		recvCells = lb_nTotalComms;
	} else {
		sendCells = lb_nTotalComms-lb_nForceComms;
		recvCells = lb_nForceComms;
	}
	totalOperations = sendCells + recvCells;

	MPI_Status stat;

	empty_mpi_buffers();

	for (i = 0; i<sendCells;++i){
		/*Send data away*/
		lb_copyCellDataToSend(&lb_send_buf[i], lb_sendCells[i], lb_nSendCells[i], pack_func, lb_commIndexToCpu[i]);
		isend_buf(&lb_send_buf[i], lb_commIndexToCpu[i], &lb_req_send[i]);
		lb_requests[i] = lb_req_send[i];
		lb_request_indices[i] = -1; /* Indicates no processing required */
	}

	for (i = 0; i<recvCells;++i){
		/*Start receiving data*/
		k = (lb_nTotalComms-1)-i;
		irecv_buf(&lb_recv_buf[k], lb_commIndexToCpu[k], &lb_req_recv[k]);
		lb_requests[i+sendCells] = lb_req_recv[k];
		lb_request_indices[i+sendCells] = k;
	}

	/*Receive and process data as soon as something is available*/
	for (i = totalOperations; i>0; i--){
		int finished;
		MPI_Waitany(i, lb_requests, &finished, &stat);
		int ind = lb_request_indices[finished];
		if (ind != -1){
			MPI_Get_count(&stat, REAL, &lb_recv_buf[ind].n);

			lb_unpackCellDataFromBuffer(&lb_recv_buf[ind], lb_commIndexToCpu[ind], (*unpack_func));
		}
		lb_requests[finished] = lb_requests[i-1];
		lb_request_indices[finished] = lb_request_indices[i-1];
	}
}

/**
 * Fills the message buffer with the information to be send to other CPUs
 * using the given pack function
 * Additionally to the packed information, a header with cell index, is created
 */
void lb_copyCellDataToSend(msgbuf* buf, ivektor* cellList, int cellListSize,
		void (*pack_func)(msgbuf*, int, int, int, vektor), int cpu){
	int k,j;
	for (k = 0; k < cellListSize; ++k) {
		j = buf->n;
		buf->data[j++] = (real)(cellList[k].x + lb_cell_offset.x);
		buf->data[j++] = (real)(cellList[k].y + lb_cell_offset.y);
		buf->data[j++] = (real)(cellList[k].z + lb_cell_offset.z);
		buf->n += 3;

		vektor vec = {0., 0., 0.};
#ifdef NBLIST
		if (pbc_dirs.x && cellList[k].x + lb_cell_offset.x == 0 && lb_pbcFlag[cpu].x)
			vec.x = (box_x.x + box_y.x + box_z.x);
		if (pbc_dirs.y && cellList[k].y + lb_cell_offset.y == 0 && lb_pbcFlag[cpu].y)
			vec.y = (box_x.y + box_y.y + box_z.y);
		if (pbc_dirs.z && cellList[k].z + lb_cell_offset.z == 0 && lb_pbcFlag[cpu].z)
			vec.z = (box_x.z + box_y.z + box_z.z);

		if (pbc_dirs.x && cellList[k].x + lb_cell_offset.x == global_cell_dim.x-1 && lb_pbcFlag[cpu].x)
			vec.x = -(box_x.x + box_y.x + box_z.x);
		if (pbc_dirs.y && cellList[k].y + lb_cell_offset.y == global_cell_dim.y-1 && lb_pbcFlag[cpu].y)
			vec.y = -(box_x.y + box_y.y + box_z.y);
		if (pbc_dirs.z && cellList[k].z + lb_cell_offset.z == global_cell_dim.z-1 && lb_pbcFlag[cpu].z)
			vec.z = -(box_x.z + box_y.z + box_z.z);
#endif
		(*pack_func)(buf, cellList[k].x, cellList[k].y, cellList[k].z, vec);
	}
}

void lb_copyForcesDataToSend(msgbuf* buf, ivektor* cellList, int cellListSize,
		void (*pack_func)(msgbuf*, int, int, int)){
	int i ,x,y,z;
	for (i = 0; i < cellListSize; ++i) {
		x = cellList[i].x + lb_cell_offset.x;
		y = cellList[i].y + lb_cell_offset.y;
		z = cellList[i].z + lb_cell_offset.z;

		/* Correct index for PBC, always send correct cell indices*/
		if (x<0) x += global_cell_dim.x;
		if (x>=global_cell_dim.x) x -= global_cell_dim.x;

		if (y<0) y += global_cell_dim.y;
		if (y>=global_cell_dim.y) y -= global_cell_dim.y;

		if (z<0) z += global_cell_dim.z;
		if (z>=global_cell_dim.z) z -= global_cell_dim.z;

		buf->data[buf->n++] = (real)x;
		buf->data[buf->n++] = (real)y;
		buf->data[buf->n++] = (real)z;
		(*pack_func)(buf, cellList[i].x, cellList[i].y, cellList[i].z);
	}
}

/**
 * Reads data from the message buffer and process this with the unpack_func and copy_func
 * If the returned value of forward is not zero, the cell must travel over edges
 * The value of dataChunkSize indicates the size of data for one cell excluding the values in
 * the header. The value "buf->n" is not changed.
 */
void lb_unpackCellDataFromBuffer(msgbuf* buf, int originCPU,
		void (*unpack_func)(msgbuf*, int, int, int)) {
	int j;
	ivektor pos;
	j = buf->n;
	buf->n=0;
	while (buf->n < j){
		/*read entries of the header*/
		pos.x = (int)buf->data[buf->n++];
		pos.y = (int)buf->data[buf->n++];
		pos.z = (int)buf->data[buf->n++];
#ifdef DEBUG
		if (pos.x<0 || pos.y<0 || pos.z<0 ||
				pos.x>=global_cell_dim.x || pos.y>=global_cell_dim.y || pos.z>=global_cell_dim.z)
			error("LB: received impossible cell index");
#endif
		/*Correct positions over periodic images*/
		if (pbc_dirs.x == 1 && lb_pbcFlag[originCPU].x){
			if (pos.x <= 0) pos.x += global_cell_dim.x;
			else if (pos.x >= global_cell_dim.x - 1) pos.x -= global_cell_dim.x;
		}

		if (pbc_dirs.y == 1 && lb_pbcFlag[originCPU].y){
			if (pos.y <= 0) pos.y += global_cell_dim.y;
			else if (pos.y == global_cell_dim.y - 1) pos.y -= global_cell_dim.y;
		}

		if (pbc_dirs.z == 1 && lb_pbcFlag[originCPU].z){
			if (pos.z <= 0) pos.z += global_cell_dim.z;
			else if (pos.z == global_cell_dim.z - 1) pos.z -= global_cell_dim.z;
		}

		pos.x -= lb_cell_offset.x;
		pos.y -= lb_cell_offset.y;
		pos.z -= lb_cell_offset.z;
#ifdef DEBUG
		if (PTR_3D_VV(cell_array, pos, cell_dim)->lb_cell_type != LB_BUFFER_CELL){
			error("LB: received data for non-buffer cell");
		}
#endif

		(*unpack_func)( buf, pos.x, pos.y, pos.z );
	}
}

void lb_unpackForcesDataFromBuffer(msgbuf* buf,
		void (*unpack_func)(msgbuf*, int, int, int)) {
	int j;
	ivektor pos;
	j = buf->n;
	buf->n = 0;

	while (buf->n < j){
		/*read entries of the header*/
		pos.x = (int)buf->data[buf->n++];
		pos.y = (int)buf->data[buf->n++];
		pos.z = (int)buf->data[buf->n++];
#ifdef DEBUG
		if (pos.x<0 || pos.y<0 || pos.z<0 ||
				pos.x>=global_cell_dim.x || pos.y>=global_cell_dim.y || pos.z>=global_cell_dim.z){
			error("LB: received impossible cell index for forces");
		}
#endif

		pos.x -= lb_cell_offset.x;
		pos.y -= lb_cell_offset.y;
		pos.z -= lb_cell_offset.z;

#ifdef DEBUG
		if (PTR_3D_VV(cell_array, pos, cell_dim)->lb_cell_type != LB_REAL_CELL)
			error("LB: received data for non-real cell");
#endif
		(*unpack_func)( buf, pos.x, pos.y, pos.z );
	}
}

void lb_relocateParticles(ivektor oldOffset, ivektor oldSize, cell* oldCells){
	int x,y,z, x2, y2, z2, i, to_cpu;
	cell *cell_new, *cell_old;
	msgbuf *buf;
	int *numCellsToSend = malloc(num_cpus*sizeof *numCellsToSend);
	int *numCellsToReceive = malloc(num_cpus*sizeof *numCellsToReceive);
	if (numCellsToSend == NULL || numCellsToReceive == NULL)
		error("Cannot allocate send/recv Buffer in lb_relocateParticles");

	for (x = 0; x < num_cpus; x++){
		numCellsToReceive[x] = 0;
		numCellsToSend[x] = 0;
	}

	/*Count how many cells have been lost to which domain*/
	for (x=0; x<cell_dim.x; ++x){
		for (y=0; y<cell_dim.y; ++y){
			for (z=0; z<cell_dim.z; ++z){
				cell_new = PTR_3D_V(cell_array, x, y, z, cell_dim);

				x2 = x+lb_cell_offset.x-oldOffset.x;
				y2 = y+lb_cell_offset.y-oldOffset.y;
				z2 = z+lb_cell_offset.z-oldOffset.z;
				cell_old = lb_accessCell(oldCells, x2, y2, z2, oldSize);

				/*Cell changed from real to buffer*/
				if (cell_old != NULL &&
						cell_old->lb_cell_type == LB_REAL_CELL &&
						cell_new->lb_cell_type != LB_REAL_CELL){
					numCellsToSend[cell_new->lb_cpu_affinity]++;
				}
			}
		}
	}

	/* Count how many cells have been gained from which domain*/
	for (x=0; x<oldSize.x; ++x){
		for (y=0; y<oldSize.y; ++y){
			for (z=0; z<oldSize.z; ++z){
				cell_old = PTR_3D_V(oldCells, x, y, z, oldSize);

				x2 = x-lb_cell_offset.x+oldOffset.x;
				y2 = y-lb_cell_offset.y+oldOffset.y;
				z2 = z-lb_cell_offset.z+oldOffset.z;
				cell_new = lb_accessCell(cell_array,x2,y2,z2,cell_dim);

				if (cell_old->lb_cell_type == LB_REAL_CELL && cell_new == NULL){
					error("Cannot relocate cell, new cell is null\n");
				}

				/*Cell changed from real to buffer*/
				if (cell_new != NULL &&
						cell_old->lb_cell_type != LB_REAL_CELL &&
						cell_new->lb_cell_type == LB_REAL_CELL){
					numCellsToReceive[cell_old->lb_cpu_affinity]++;
				}
			}
		}
	}

	/*Allocate buffer for send & receive*/
	msgbuf *sendBuf = NULL;
	memalloc(&sendBuf, num_cpus, sizeof(msgbuf), sizeof(void*), 0, 1, "sendBuf");
	msgbuf *recvBuf = NULL;
	memalloc(&recvBuf, num_cpus, sizeof(msgbuf), sizeof(void*), 0, 1, "recvBuf");

	int numSendTo = 0;
	int numReceiveFrom = 0;

	if (sendBuf == NULL || recvBuf == NULL)
		error("Cannot allocate send/recv Buffer in lb_relocateParticles");

	for (x = 0; x < num_cpus; x++){
		if (numCellsToReceive[x] != 0) {
			alloc_msgbuf(&recvBuf[x], atom_size*numCellsToReceive[x]*lb_largest_cell);
			numReceiveFrom++;
		}
		if (numCellsToSend[x] != 0){
			alloc_msgbuf(&sendBuf[x], atom_size*numCellsToSend[x]*lb_largest_cell);
			numSendTo++;
		}
	}

	/*Put everything in the buffer*/
	/*new geometry is valid, now copy the content from the old to the new domains*/
	for (x=1; x<cell_dim.x-1; ++x){
		for (y=1; y<cell_dim.y-1; ++y){
			for (z=1; z<cell_dim.z-1; ++z){
				cell_new = PTR_3D_V(cell_array, x, y, z, cell_dim);
				if(cell_new->lb_cell_type == LB_REAL_CELL){
					x2 = x+lb_cell_offset.x-oldOffset.x;
					y2 = y+lb_cell_offset.y-oldOffset.y;
					z2 = z+lb_cell_offset.z-oldOffset.z;
					cell_old = lb_accessCell(oldCells,x2,y2,z2,oldSize);

					if (cell_old->lb_cell_type == LB_REAL_CELL){
						/*Alloc new cell and copy content from old*/
						alloc_cell(cell_new, cell_old->n_max);
						for (i=0; i<cell_old->n; ++i){
							copy_atom_cell_cell(cell_new, cell_new->n, cell_old, i);
							++cell_new->n;
						}
						alloc_cell(cell_old, 0);
					}
				}
			}
		}
	}

	/*Send data from cells that transformed from real cells to buffer cells*/
	/*to the CPUs that now hold this data as a real cell*/
	for (x=0; x<cell_dim.x; ++x){
		for (y=0; y<cell_dim.y; ++y){
			for (z=0; z<cell_dim.z; ++z){
				cell_new = PTR_3D_V(cell_array, x, y, z, cell_dim);

				x2 = x+lb_cell_offset.x-oldOffset.x;
				y2 = y+lb_cell_offset.y-oldOffset.y;
				z2 = z+lb_cell_offset.z-oldOffset.z;
				cell_old = lb_accessCell(oldCells,x2,y2,z2,oldSize);

				/*Cell changed from real to buffer*/
				if (cell_old != NULL &&
						cell_old->lb_cell_type == LB_REAL_CELL &&
						cell_new->lb_cell_type != LB_REAL_CELL){

					to_cpu = cell_new->lb_cpu_affinity;
					buf = &sendBuf[to_cpu];
					for (i=0; i<cell_old->n; i++) {
						copy_one_atom(buf, to_cpu, cell_old, i, 0);
#ifdef CLONE
						int clone;
						if (l < cell_old->n-nclones)
						for (clone=1; clone<nclones; clone++)
						copy_one_atom( buf, to_cpu, cell_old, i+clone, 0);
						else /* we are dealing with the last in the stack */
						for (clone=1; clone<nclones; clone++)
						copy_one_atom( buf, to_cpu, cell_old, i, 0);
#endif
					}

				}
			}
		}
	}

	/*Send/receive buffers*/
	int totalOperations = numReceiveFrom + numSendTo;
	MPI_Request *requests = malloc(totalOperations * sizeof *requests);
	int *indices = malloc(totalOperations * sizeof *indices);
	MPI_Status stat;

	x = 0;
	for (i = 0; i < num_cpus; ++i) {
		/*Send data away*/
		if (numCellsToSend[i] != 0){
			isend_buf(&sendBuf[i], i, &requests[x]);
			indices[x++] = -1;
		}
		/*Start receiving data*/
		if (numCellsToReceive[i] != 0){
			irecv_buf(&recvBuf[i], i, &requests[x]);
			indices[x++] = i;
		}
	}

	/*Receive and process data as soon as something is available*/
	for (i = totalOperations; i>0; i--){
		int finished;
		MPI_Waitany(i, requests, &finished, &stat);
		int ind = indices[finished];
		if (ind != -1){
			MPI_Get_count(&stat, REAL, &recvBuf[ind].n);
			process_buffer( &recvBuf[ind]);
		}
		requests[finished] = requests[i-1];
		indices[finished] = indices[i-1];
	}

	free(requests);
	free(indices);

	/*Clean up*/
	for (x = 0; x < num_cpus; x++){
		if (numCellsToReceive[x] != 0) free_msgbuf(&recvBuf[x]);
		if (numCellsToSend[x] != 0) free_msgbuf(&sendBuf[x]);
	}
	free(numCellsToSend);
	free(numCellsToReceive);
}


int lb_isGeometryChangeValid(){
	int i,j,k,x,y,z;
	int i2,j2,k2;
	int minCell[3];
	int maxCell[3];
	ivektor index, newDomainOffset, newDomainSize;
	cell *cell;

	/*Basic check of the geometry*/
	if (!lb_isGeometryValid(&lb_domain)) return 0;

	/* Further checks are not based on the geometry of the domain, but on the level
	 * of cells directly to identify special cases if the geometry has changed to fast*/
	maxCell[0] = 0; maxCell[1] = 0; maxCell[2] = 0;
	minCell[0] = global_cell_dim.x;  minCell[1] = global_cell_dim.y; minCell[2] = global_cell_dim.z;

	/*Compute the boundaries of the modified domain*/
	for (i=0; i<8; i++){
		index = cell_coord(lb_domain.corners[i].discretizedP.x,
						   lb_domain.corners[i].discretizedP.y,
						   lb_domain.corners[i].discretizedP.z);

		if (index.x<minCell[0]) minCell[0] = index.x;
		if (index.y<minCell[1]) minCell[1] = index.y;
		if (index.z<minCell[2]) minCell[2] = index.z;
		if (index.x>maxCell[0]) maxCell[0] = index.x;
		if (index.y>maxCell[1]) maxCell[1] = index.y;
		if (index.z>maxCell[2]) maxCell[2] = index.z;
	}

	newDomainOffset.x = minCell[0]-1;
	newDomainOffset.y = minCell[1]-1;
	newDomainOffset.z = minCell[2]-1;

	newDomainSize.x = maxCell[0]-minCell[0]+3;
	newDomainSize.y = maxCell[1]-minCell[1]+3;
	newDomainSize.z = maxCell[2]-minCell[2]+3;
	/*Check if the domain size has been modified more than one layer of cells*/
	if (abs(lb_cell_offset.x - newDomainOffset.x) > 1 ||
		abs(lb_cell_offset.y - newDomainOffset.y) > 1 ||
		abs(lb_cell_offset.z - newDomainOffset.z) > 1 )
		return 0;

	if (abs(lb_cell_offset.x-newDomainOffset.x)-abs(cell_dim.x - newDomainSize.x) > 1 ||
		abs(lb_cell_offset.y-newDomainOffset.y)-abs(cell_dim.y - newDomainSize.y) > 1 ||
		abs(lb_cell_offset.z-newDomainOffset.z)-abs(cell_dim.z - newDomainSize.z) > 1 )
			return 0;

	/* Check if an empty cell has been turned into a real cell in the new domain.
	 * Reject geometry if variable communication is disabled*/
	if (lb_balancingType==0){
		for (i = 0; i<newDomainSize.x; i++){
			for (j = 0; j<newDomainSize.y; j++){
				for (k = 0; k<newDomainSize.z; k++){
					x = i+newDomainOffset.x - lb_cell_offset.x;
					y = j+newDomainOffset.y - lb_cell_offset.y;
					z = k+newDomainOffset.z - lb_cell_offset.z;
					cell = lb_accessCell(cell_array,x,y,z, cell_dim);
					if (cell == NULL || cell->lb_cell_type == LB_EMPTY_CELL){
						if (lb_isPointInDomain(lb_getCellCenter(x+lb_cell_offset.x, y+lb_cell_offset.y, z+lb_cell_offset.z),&lb_domain)){
							return 0;
						}
					}
				}
			}
		}
		/* Check if a real cell has been turned into an empty cell in the new domain
		 * This is the case, if not at least the cell or one of its 26 neighbors was inside*/
		for (i = 1; i<cell_dim.x-1; i++){
			for (j = 1; j<cell_dim.y-1; j++){
				for (k = 1; k<cell_dim.z-1; k++){
					x = i+lb_cell_offset.x;
					y = j+lb_cell_offset.y;
					z = k+lb_cell_offset.z;
					cell = lb_accessCell(cell_array,i,j,k, cell_dim);
					if (cell->lb_cell_type == LB_REAL_CELL){
						int oneInside = 0;
						for (i2=-1; i2<=1; i2++){
							for (j2=-1; j2<=1; j2++){
								for (k2=-1; k2<=1; k2++){
									if (oneInside==0 && lb_isPointInDomain(lb_getCellCenter(x+i2, y+j2, z+k2),&lb_domain)) oneInside++;
								}
							}
						}
						if (oneInside==0) return 0;
					}
				}
			}
		}
	}

	/*All test successful, new geometry seems valid*/
	return 1;
}


int lb_isGeometryValid(lb_domainInfo *dom){

	if (lb_balancingType == 1){
	/* In comparision to the set of rules in the plimpton scheme, these values are less strict
	 * Except for the tetrahedron volumes, the checks are not necessary, but prevent the domain
	 * from degenerating to severely, which can cause problems in the convergence behavior.
	 * The movement of corners is otherwise likely to get stuck in local minima.
	 */
		/* Tetrahedral subvolumes in the domain must be positively oriented */
		/* Self-intersecting cubes are bad, very bad*/
		/* Test all four permutations of how the cube can be split into tetrahedrons*/
		if (lb_getTetraederVolumeIndexed(0, 5, 4, 7, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(0, 3, 1, 7, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(0, 1, 5, 7, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(0, 4, 6, 7, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(0, 6, 2, 7, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(0, 2, 3, 7, &lb_domain) <= 0) return 0;

		if (lb_getTetraederVolumeIndexed(1, 7, 5, 6, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(1, 2, 3, 6, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(1, 3, 7, 6, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(1, 5, 4, 6, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(1, 4, 0, 6, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(1, 0, 2, 6, &lb_domain) <= 0) return 0;

		if (lb_getTetraederVolumeIndexed(2, 4, 6, 5, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(2, 1, 0, 5, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(2, 0, 4, 5, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(2, 6, 7, 5, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(2, 7, 3, 5, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(2, 3, 1, 5, &lb_domain) <= 0) return 0;

		if (lb_getTetraederVolumeIndexed(3, 6, 7, 4, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(3, 0, 2, 4, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(3, 2, 6, 4, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(3, 7, 5, 4, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(3, 5, 1, 4, &lb_domain) <= 0) return 0;
		if (lb_getTetraederVolumeIndexed(3, 1, 0, 4, &lb_domain) <= 0) return 0;

		//Additionally prevent the collapse of the corners in the domain
		//This would yield a topological different domain geometry
		if (lb_getTetraederVolumeIndexed(0, 1, 4, 2, &lb_domain) < 0) return 0;
		if (lb_getTetraederVolumeIndexed(1, 5, 4, 7, &lb_domain) < 0) return 0;
		if (lb_getTetraederVolumeIndexed(2, 4, 6, 7, &lb_domain) < 0) return 0;
		if (lb_getTetraederVolumeIndexed(2, 7, 3, 1, &lb_domain) < 0) return 0;
		if (lb_getTetraederVolumeIndexed(0, 4, 6, 5, &lb_domain) < 0) return 0;
		if (lb_getTetraederVolumeIndexed(0, 1, 5, 3, &lb_domain) < 0) return 0;
		if (lb_getTetraederVolumeIndexed(0, 6, 2, 3, &lb_domain) < 0) return 0;
		if (lb_getTetraederVolumeIndexed(5, 3, 7, 6, &lb_domain) < 0) return 0;

	} else if (lb_balancingType==0) {
		/* rules that enforce communication with nearest neigbors*/
		const real minCellSize = MIN(MIN(lb_cell_size.x, lb_cell_size.y), lb_cell_size.z);
		const real minDistanceFromDividingPlane = minCellSize * 1.733;

		//Prevent the collapse of the corners in the domain
		//This would yield a topological different domain geometry
		if (lb_getTetraederVolumeIndexed(0, 1, 4, 2, &lb_domain) < 0) return 0;
		if (lb_getTetraederVolumeIndexed(1, 5, 4, 7, &lb_domain) < 0) return 0;
		if (lb_getTetraederVolumeIndexed(2, 4, 6, 7, &lb_domain) < 0) return 0;
		if (lb_getTetraederVolumeIndexed(2, 7, 3, 1, &lb_domain) < 0) return 0;
		if (lb_getTetraederVolumeIndexed(0, 4, 6, 5, &lb_domain) < 0) return 0;
		if (lb_getTetraederVolumeIndexed(0, 1, 5, 3, &lb_domain) < 0) return 0;
		if (lb_getTetraederVolumeIndexed(0, 6, 2, 3, &lb_domain) < 0) return 0;
		if (lb_getTetraederVolumeIndexed(5, 3, 7, 6, &lb_domain) < 0) return 0;

		/*From each corner, six planes (the three faces on the original cube,
		  eachsplit into two triangles) must have a minimum distance
		  that guarantees that at least one cell is in between.*/
		if (lb_getDistanceToPlane(0, 4, 6, 7, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(0, 4, 5, 7, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(0, 2, 3, 7, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(0, 2, 6, 7, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(0, 1, 3, 7, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(0, 1, 5, 7, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(1, 2, 7, 3, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(1, 7, 2, 6, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(1, 0, 6, 2, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(1, 0, 4, 6, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(1, 4, 6, 7, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(1, 4, 7, 5, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(2, 1, 3, 7, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(2, 1, 7, 5, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(2, 0, 1, 5, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(2, 0, 5, 4, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(2, 4, 5, 7, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(2, 4, 7, 6, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(3, 2, 0, 6, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(3, 0, 6, 4, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(3, 4, 7, 5, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(3, 4, 7, 6, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(3, 0, 1, 5, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(3, 0, 4, 5, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(4, 0, 1, 3, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(4, 0, 3, 2, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(4, 7, 6, 2, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(4, 7, 3, 2, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(4, 3, 1, 7, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(4, 1, 7, 5, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(5, 2, 3, 7, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(5, 2, 7, 6, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(5, 0, 6, 4, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(5, 0, 2, 6, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(5, 1, 0, 3, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(5, 3, 2, 0, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(6, 4, 5, 0, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(6, 0, 1, 5, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(6, 1, 7, 3, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(6, 1, 7, 5, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(6, 0, 2, 3, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(6, 3, 1, 0, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(7, 0, 3, 1, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(7, 0, 3, 2, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(7, 0, 6, 4, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(7, 0, 2, 6, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(7, 0, 1, 5, dom) < minDistanceFromDividingPlane) return 0;
		if (lb_getDistanceToPlane(7, 0, 4, 5, dom) < minDistanceFromDividingPlane) return 0;


//		const real minAngle = cos(22.5 * 3.141592654/180.);
//		const real warpMax =  sin(15 * 3.141592654/180.);
//		real a;
//
//		/*Check minimum distance from a corner to opposite faces*/
//		if (lb_getDistanceToPlane(0, 4, 1, 3, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(0, 1, 2, 6, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(0, 2, 4, 5, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(1, 0, 5, 7, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(1, 5, 3, 2, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(1, 3, 0, 4, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(2, 0, 3, 7, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(2, 3, 6, 4, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(2, 6, 0, 1, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(3, 1, 7, 6, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(3, 7, 2, 0, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(3, 2, 1, 5, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(4, 6, 5, 1, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(4, 5, 0, 2, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(4, 0, 6, 7, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(5, 1, 4, 6, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(5, 4, 7, 3, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(5, 7, 1, 0, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(6, 4, 2, 3, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(6, 2, 7, 5, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(6, 7, 4, 0, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(7, 5, 6, 2, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(7, 6, 3, 1, dom) < minDistanceFromDividingPlane) return 0;
//		if (lb_getDistanceToPlane(7, 3, 5, 4, dom) < minDistanceFromDividingPlane) return 0;
//		/* Check all 24 angles at all six faces, must not be smaller than the threshold to avoid
//		 * illegal geometries */
//		a = lb_angleCos(4, 0, 2, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(4, 0, 1, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(2, 0, 1, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(5, 1, 3, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(3, 1, 0, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(5, 1, 0, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(1, 3, 7, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(7, 3, 2, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(2, 3, 1, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(6, 2, 0, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(0, 2, 3, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(3, 2, 6, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(6, 4, 0, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(0, 4, 5, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(5, 4, 6, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(4, 5, 7, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(7, 5, 1, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(1, 5 ,4, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(6, 7, 5, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(5, 7, 3, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(3, 7, 6, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(7, 6, 4, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(2, 6, 7, dom); if (a>minAngle) return 0;
//		a = lb_angleCos(4, 6, 2, dom); if (a>minAngle) return 0;
//
//		/*Check the warpage on the six faces*/
//		if (lb_computeWarpage(1, 0, 5, 4, dom) > warpMax) return 0;
//		if (lb_computeWarpage(3, 1, 7, 5, dom) > warpMax) return 0;
//		if (lb_computeWarpage(3, 2, 7, 6, dom) > warpMax) return 0;
//		if (lb_computeWarpage(2, 6, 0, 4, dom) > warpMax) return 0;
//		if (lb_computeWarpage(6, 4, 7, 5, dom) > warpMax) return 0;
//		if (lb_computeWarpage(2, 0, 3, 1, dom) > warpMax) return 0;
	}
	return 1;
}
