#include "imd.h"

/******************************************************************************
*
* imd_loadBalance -- primary routines for dynamic load-balancing in IMD
* between MPI processes
*
* The load balancer requires certain modification in the communication
* systems, the required functions (send_cells, send_forces plus helper methods)
* are stored in separate files. Two different communication schemes are available
*
* imd_loadBalance_plimpton: Communication using the plimpton scheme
* imd_loadBalance_direct: Communication using direct communication between CPUs
*
* The plimpton scheme is faster on computing cluster higher latency communication.
* Fewer messages are exchanged, but the overhead is larger to pass information of
* cells several times.
* In very low latency systems, the use of direct communication can be beneficial.
* Here each cpu communicates with all neighbors immediately, reducing overhead,
* but more message are needed to send cells and forces
******************************************************************************/

void init_loadBalance() {
	int i,j;
	int x, y, z;

	if (lb_init) return;
	lb_init = 1;

	if (lb_contractionRate == -1)
		lb_contractionRate = SQRT(cellsz);

#ifndef AR
	error("actio=reactio condition must be fulfilled for Load Balancing");
#endif
#ifdef TWOD
	error("2D not supported by Load Balancing");
#endif
#if defined(NPT_AXIAL) || defined(NPT_ISO) || defined(NPT_FULL) || defined(DEFORM) || defined(HOMDEF)
	error("options causing the volume/geometry of the box to be changed are currently not supported using load balancing");
#endif
	if ( (cpu_dim.x % 2 == 1 && pbc_dirs.x == 1) ||
			(cpu_dim.y % 2 == 1 && pbc_dirs.y == 1) ||
			(cpu_dim.z % 2 == 1 && pbc_dirs.z == 1))
			error("Load balancing requires that cpu_dim must be a multiple of two in periodic directions");
	if (box_x.y != 0. || box_x.z != 0. || box_y.x != 0. || box_y.z != 0. || box_z.x != 0. || box_z.y != 0.)
		warning("WARNING!!! Non-orthogonal boxes are not tested together with load balancing.");

	if (myid==0){
		printf("LOAD BALANCING: load balancing steps before simulation %i\n", lb_preRuns);
		printf("LOAD BALANCING: contraction rate %f\n", lb_contractionRate);
		printf("LOAD BALANCING: balancing every %i steps\n", lb_frequency);
	}

	lb_initDirect();

	for (i=0; i<26;i++){
		lb_pbcCorrection[i].x = 0;
		lb_pbcCorrection[i].y = 0;
		lb_pbcCorrection[i].z = 0;
	}

	/* setup the pbc-correction vectors, the vectors are always set even in non-periodic systems*/
	if (my_coord.x == 0){
		lb_pbcCorrection[0].x  = 1; lb_pbcCorrection[1].x  = 1;
		lb_pbcCorrection[2].x  = 1; lb_pbcCorrection[3].x  = 1;
		lb_pbcCorrection[4].x  = 1; lb_pbcCorrection[5].x  = 1;
		lb_pbcCorrection[6].x  = 1; lb_pbcCorrection[7].x  = 1;
		lb_pbcCorrection[8].x  = 1;
	}
	if (my_coord.x == cpu_dim.x-1){
		lb_pbcCorrection[17].x = -1; lb_pbcCorrection[18].x = -1;
		lb_pbcCorrection[19].x = -1; lb_pbcCorrection[20].x = -1;
		lb_pbcCorrection[21].x = -1; lb_pbcCorrection[22].x = -1;
		lb_pbcCorrection[23].x = -1; lb_pbcCorrection[24].x = -1;
		lb_pbcCorrection[25].x = -1;
	}
	if (my_coord.y == 0){
		lb_pbcCorrection[0].y  = 1; lb_pbcCorrection[1].y  = 1;
		lb_pbcCorrection[2].y  = 1; lb_pbcCorrection[9].y  = 1;
		lb_pbcCorrection[10].y = 1; lb_pbcCorrection[11].y = 1;
		lb_pbcCorrection[17].y = 1; lb_pbcCorrection[18].y = 1;
		lb_pbcCorrection[19].y = 1;
	}
	if (my_coord.y == cpu_dim.y-1){
		lb_pbcCorrection[6].y  = -1; lb_pbcCorrection[7].y  = -1;
		lb_pbcCorrection[8].y  = -1; lb_pbcCorrection[14].y = -1;
		lb_pbcCorrection[15].y = -1; lb_pbcCorrection[16].y = -1;
		lb_pbcCorrection[23].y = -1; lb_pbcCorrection[24].y = -1;
		lb_pbcCorrection[25].y = -1;
	}
	if (my_coord.z == 0){
		lb_pbcCorrection[0].z  = 1; lb_pbcCorrection[3].z  = 1;
		lb_pbcCorrection[6].z  = 1; lb_pbcCorrection[9].z  = 1;
		lb_pbcCorrection[12].z = 1; lb_pbcCorrection[14].z = 1;
		lb_pbcCorrection[17].z = 1; lb_pbcCorrection[20].z = 1;
		lb_pbcCorrection[23].z = 1;
	}
	if (my_coord.z == cpu_dim.z-1){
		lb_pbcCorrection[2].z  = -1; lb_pbcCorrection[5].z  = -1;
		lb_pbcCorrection[8].z  = -1; lb_pbcCorrection[11].z = -1;
		lb_pbcCorrection[13].z = -1; lb_pbcCorrection[16].z = -1;
		lb_pbcCorrection[19].z = -1; lb_pbcCorrection[22].z = -1;
		lb_pbcCorrection[25].z = -1;
	}

	/* Initializing corners of CPU domain boundaries
	 * Initially the domains have the same dimension as without load balancing.
	 * The eight corners are moved independently and cells inside the domain
	 * are handled by the CPU assigned to this domain.
	 * Corners placed directly at the periodic boundary are partially restricted in their
	 * degrees of freedom to move, since the shape of the total simulation box is not changed
	 * (property "fixed" is set to 1).
	 */
	for (i = 0; i<8; i++){
		lb_domain.corners[i].fixed.x = 0;
		lb_domain.corners[i].fixed.y = 0;
		lb_domain.corners[i].fixed.z = 0;
		for (j=0; j<11; j++)
			lb_localCommPartners[i][j] = 0;
	}

	/*
	 *	The domain corner positions
	 *	If corners are moved in their positions, the faces are not necessarily planar
	 *	By definition, the faces are splitted into two triangles along the diagonals
	 *	0-5, 0-6, 0-3, 4-7, 1-7, 2-7
	 *	The geometry
	 *
	 *   6----7
	 *  /|   /|   z
	 * 4----5 |   ^  y
	 * | 2--|-3   | ^
	 * |/   |/    |/
	 * 0----1     --->x
	 */
	for (i = 0; i<8; i++){
		z = (i&4)>>2;
		y = (i&2)>>1;
		x = (i&1);
		/* Current position of corner */
		lb_domain.corners[i].p.x = ((box_x.x + box_y.x + box_z.x) / cpu_dim.x) * (my_coord.x + x);
		lb_domain.corners[i].p.y = ((box_x.y + box_y.y + box_z.y) / cpu_dim.y) * (my_coord.y + y);
		lb_domain.corners[i].p.z = ((box_x.z + box_y.z + box_z.z) / cpu_dim.z) * (my_coord.z + z);

		/* Reference position of corner, will never be changed */
		lb_domain.corners[i].ref = lb_domain.corners[i].p;

		/* Fixing corners at the simulation box boundaries */
		if (x == 0 && my_coord.x == 0) lb_domain.corners[i].fixed.x = 1;
		if (y == 0 && my_coord.y == 0) lb_domain.corners[i].fixed.y = 1;
		if (z == 0 && my_coord.z == 0) lb_domain.corners[i].fixed.z = 1;
		if (x == 1 && my_coord.x == cpu_dim.x - 1) lb_domain.corners[i].fixed.x = 1;
		if (y == 1 && my_coord.y == cpu_dim.y - 1) lb_domain.corners[i].fixed.y = 1;
		if (z == 1 && my_coord.z == cpu_dim.z - 1) lb_domain.corners[i].fixed.z = 1;
	}

	lb_updateDomain(&lb_domain);


	/* Create local groups
	 * To move a corner, the load of eight adjacent domains is needed
	 * Thus each CPU is part of eight different communication groups.
	 * Since some MPI implementations do support only a limited number of
	 * MPI_Communicators (e.g. 2048 in MPICH2), using MPI communications cannot
	 * be used if there are several thousand CPUs to be used in any case.
	 * Therefore, collective operations in the group of up to eight CPUs
	 * are performed by conventional MPI_Send/MPI_Recv point-to-point communication
	 */
	int procs[8];
	ivektor tmp;
	int my_index_inGroup;

	int cornersX = pbc_dirs.x == 1 ? cpu_dim.x : cpu_dim.x + 1;
	int cornersY = pbc_dirs.y == 1 ? cpu_dim.y : cpu_dim.y + 1;
	int cornersZ = pbc_dirs.z == 1 ? cpu_dim.z : cpu_dim.z + 1;

	/* Create a communication group for each corner
	 * Each of the 1-8 CPUs that share the corner are included in the group.
	 * The CPU with rank 0 is the master for the node. In fully periodic systems, each domain
	 * is master of corner 0. In open systems, CPUs at the upper boundaries
	 * are masters of multiple corners
	 */
	for (x=0; x<cornersX; ++x){
		for (y=0; y<cornersY; ++y){
			for (z=0; z<cornersZ; ++z){
				my_index_inGroup = -1;
				int numDomains = 1;
				int domIndex = 0;

				if (pbc_dirs.x == 1 || (x>0 && x<cpu_dim.x)) numDomains*=2;
				if (pbc_dirs.y == 1 || (y>0 && y<cpu_dim.y)) numDomains*=2;
				if (pbc_dirs.z == 1 || (z>0 && z<cpu_dim.z)) numDomains*=2;

				for (i=0; i<8; i++){
					int assign = 1;
					tmp.x = x - (i & 1);
					tmp.y = y - ((i & 2) >> 1);
					tmp.z = z - ((i & 4) >> 2);

					if (pbc_dirs.x == 1 && tmp.x == -1) tmp.x = cpu_dim.x-1;
					if (pbc_dirs.y == 1 && tmp.y == -1) tmp.y = cpu_dim.y-1;
					if (pbc_dirs.z == 1 && tmp.z == -1) tmp.z = cpu_dim.z-1;

					if (tmp.x<0 || tmp.y<0 || tmp.z<0) assign = 0;
					if (x == cpu_dim.x && tmp.x != cpu_dim.x) assign = 0;
					if (y == cpu_dim.y && tmp.y != cpu_dim.y) assign = 0;
					if (z == cpu_dim.z && tmp.z != cpu_dim.z) assign = 0;

					if (tmp.x == cpu_dim.x) tmp.x--;
					if (tmp.y == cpu_dim.y) tmp.y--;
					if (tmp.z == cpu_dim.z) tmp.z--;

					int id = *PTR_3D_VV(cpu_ranks, tmp, cpu_dim);;

					if (assign == 1){
						if (id == myid) my_index_inGroup = domIndex;
						procs[domIndex++] = id;
					}
				}

				/* Each of the up to eight domains must move the same corner at once.
				 * Depending on the coordinate of the CPU (odd or even in each direction),
 	 	 	 	 * the order of moving its local corners differs, but is synchronized with its neighbors*/
				if (my_index_inGroup != -1){
					int pos;
					tmp.x = x - my_coord.x;
					tmp.y = y - my_coord.y;
					tmp.z = z - my_coord.z;
					if (tmp.x<0) tmp.x = 1;
					if (tmp.y<0) tmp.y = 1;
					if (tmp.z<0) tmp.z = 1;
					pos = tmp.x + tmp.y*2 + tmp.z*4;

					lb_localCommPartners[pos][0] = numDomains;
					lb_localCommPartners[pos][1] = my_index_inGroup;
					for (i = 0; i<numDomains; i++)
						lb_localCommPartners[pos][i+2] = procs[i];
				}
			}
		}
	}

	if (lb_balancingType == 2){
		x_bounds = malloc((cpu_dim.x+1) * sizeof *x_bounds);
		y_bounds = malloc((cpu_dim.y+1) * sizeof *y_bounds);
		z_bounds = malloc((cpu_dim.z+1) * sizeof *z_bounds);

		for (i=0; i<=cpu_dim.x;i++)
			x_bounds[i] = global_cell_dim.x/cpu_dim.x*i;
		for (i=0; i<=cpu_dim.y;i++)
			y_bounds[i] = global_cell_dim.y/cpu_dim.y*i;
		for (i=0; i<=cpu_dim.z;i++)
			z_bounds[i] = global_cell_dim.z/cpu_dim.z*i;
	}
}

/*
 * Perform a load balancing steps.
 * Make sure to call fix_cells() before this method to ensure all atoms are located on the correct CPUs.
 * Parameters:
 * (*getLoad) - function defining the load of a CPUs, should be normalized, that the average load is equal to 1
 * (*getCog) - function for defining the center of gravity, i.e. the point at which the load is concentrated
 * reset - if not 0, a partial reset is performed in which corners are not moved according to the load,
 * but instead are moving towards their reference coordinate. Useful to get of situation where the loadbalancer has
 * converged into a local minimum
 * cycle - if called externally, always set to 0, internally used for handling errors which are solved recursively
 * Return:
 * 1: balance step successful
 * 0: balancing step has been rejected, illegal geometries have been created. Geometry rolled backed to last valid state.
 */
int balanceLoad(real (*getLoad)(void), vektor (*getCog)(void), int reset, int cycle) {
	int x,y,z, x2,y2,z2, i;
	ivektor lb_cell_offset_old, cell_dim_old, index;
	vektor p;
	minicell *cell_array_old;
	minicell *c, *cell_old, *cell_new;
	int minCell[3];
	int maxCell[3];
	real domainInfo[8];

	vektor oldPositions[8];

	for (i=0; i<8;i++){
		oldPositions[i].x = lb_domain.corners[i].p.x;
		oldPositions[i].y = lb_domain.corners[i].p.y;
		oldPositions[i].z = lb_domain.corners[i].p.z;
	}

#ifdef DEBUG
	int atomsBefore, atomsAfter;
	int totalAtomsBefore, totalsAtomsAfter;
	atomsBefore = lb_countAtoms();
#endif

	lb_cell_offset_old.x = lb_cell_offset.x;
	lb_cell_offset_old.y = lb_cell_offset.y;
	lb_cell_offset_old.z = lb_cell_offset.z;

	cell_dim_old.x = cell_dim.x;
	cell_dim_old.y = cell_dim.y;
	cell_dim_old.z = cell_dim.z;


	if(lb_balancingType == 2){
		balanceOrtho();
	} else {
		/*The actual load balancing, moving the boundaries according to load on the CPUs*/
		if (reset){
			domainInfo[0] = 0.;
			domainInfo[1] = 0.;
		} else {
			domainInfo[0] = (*getLoad)();
			domainInfo[1] = lb_getVolume();
		}
		p = (*getCog)();
		domainInfo[2] = p.x;
		domainInfo[3] = p.y;
		domainInfo[4] = p.z;

		domainInfo[5] = my_coord.x;
		domainInfo[6] = my_coord.y;
		domainInfo[7] = my_coord.z;


		/* Instead of trying to move the corners each time in a fixed order, */
		/* use every time a new randomly shuffled permutation*/
		int indizes[8];
		for (i=0; i<8;i++)
			indizes[i] = i;
		srand(lb_randomNumberGeneratorState);

		for (i=0; i<8;i++){
			int tmp, tmp2;
			tmp = rand()&7;		/*Pick a random position between 0-7*/
			tmp2 = indizes[i];	/*Swap the indices */
			indizes[i] = indizes[tmp];
			indizes[tmp] = tmp2;
		}
		lb_randomNumberGeneratorState = rand();


		/**************************/
		/* Move all corners  	  */
		/**************************/
		for (i=0; i<8;i++){
			if (reset)
				lb_moveCornersReset(domainInfo, indizes[i]);
			else
				lb_moveAllCorners(domainInfo, indizes[i]);
			lb_updateDomain(&lb_domain);
		}

	}

	/* All corners have been moved to new positions, create new geometry and communication lists*/
	maxCell[0] = 0; maxCell[1] = 0; maxCell[2] = 0;
	minCell[0] = global_cell_dim.x;  minCell[1] = global_cell_dim.y; minCell[2] = global_cell_dim.z;
	/*Find min/max of bounding box for domain*/
	for (i=0; i<8; i++){
		index.x = lb_domain.corners[i].index.x / LB_CELL_SUBLEVELS;
		index.y = lb_domain.corners[i].index.y / LB_CELL_SUBLEVELS;
		index.z = lb_domain.corners[i].index.z / LB_CELL_SUBLEVELS;

		if (index.x<minCell[0]) minCell[0] = index.x;
		if (index.y<minCell[1]) minCell[1] = index.y;
		if (index.z<minCell[2]) minCell[2] = index.z;
		if (index.x>maxCell[0]) maxCell[0] = index.x;
		if (index.y>maxCell[1]) maxCell[1] = index.y;
		if (index.z>maxCell[2]) maxCell[2] = index.z;
	}
	/*Update variables defining the domain siez*/
	lb_cell_offset.x = minCell[0]-1;
	lb_cell_offset.y = minCell[1]-1;
	lb_cell_offset.z = minCell[2]-1;

	cell_dim.x = maxCell[0]-minCell[0]+3;
	cell_dim.y = maxCell[1]-minCell[1]+3;
	cell_dim.z = maxCell[2]-minCell[2]+3;

	cellmin.x = 1; cellmin.y = 1; cellmin.z = 1;
	cellmax.x = cell_dim.x-1; cellmax.y = cell_dim.y-1; cellmax.z = cell_dim.z-1;

#ifdef DEBUG
	/*Plasusibility checks*/
	ivektor diffSize;
	ivektor diffOffset;
	diffSize.x = cell_dim_old.x - cell_dim.x;
	diffSize.y = cell_dim_old.y - cell_dim.y;
	diffSize.z = cell_dim_old.z - cell_dim.z;

	diffOffset.x = lb_cell_offset_old.x - lb_cell_offset.x;
	diffOffset.y = lb_cell_offset_old.y - lb_cell_offset.y;
	diffOffset.z = lb_cell_offset_old.z - lb_cell_offset.z;

	if (ABS(diffOffset.x)>1 || ABS(diffOffset.y)>1 || ABS(diffOffset.z)>1){
		error("Load Balance: cell_offset change too large");
	}

	if (ABS(diffOffset.x)-ABS(diffSize.x)>1 ||
			ABS(diffOffset.y)-ABS(diffSize.y)>1 || ABS(diffOffset.z)-ABS(diffSize.z)>1){
		error("Load Balance: cell_dim change too large");
	}
#endif

	cell_array_old = cell_array;

	/*Allocate and clear memory for new cell_array*/
	cell_array = malloc(cell_dim.x * cell_dim.y * cell_dim.z * sizeof *cell_array);
	if (NULL == cell_array)
		error("Cannot allocate memory for cells");
	memset(cell_array, 0, cell_dim.x * cell_dim.y * cell_dim.z * sizeof *cell_array);

	int valid = 1;

	/*Assign real cells*/
	for (x=1; x<cell_dim.x-1; ++x){
		for (y=1; y<cell_dim.y-1; ++y){
			for (z=1; z<cell_dim.z-1; ++z){
				//Identify new cell type
				cell_new = PTR_3D_V(cell_array, x, y, z, cell_dim);

				/* Flag all cells within the domain as real cells, reals cells are handled by this CPU*/
				if(lb_isPointInDomain(lb_getCellCenter(x+lb_cell_offset.x,y+lb_cell_offset.y,z+lb_cell_offset.z),&lb_domain)){
					cell_new->lb_cell_type = LB_REAL_CELL;	/* Update info in the nulled cell */
					cell_new->lb_cpu_affinity = myid;
					cell_new->lb_neighbor_index = -LB_REAL_CELL;
					/*Check state of the corresponding cell in the old domain*/
					x2 = x+lb_cell_offset.x-lb_cell_offset_old.x;
					y2 = y+lb_cell_offset.y-lb_cell_offset_old.y;
					z2 = z+lb_cell_offset.z-lb_cell_offset_old.z;
					cell_old = lb_accessCell(cell_array_old,x2,y2,z2,cell_dim_old);

					/*Illegal geometry, set valid to 0 to ensure rollback*/
					if (cell_old == NULL){
						valid = 0;
#ifdef DEBUG
						printf("Load Balance: Attempted to change a non-existing cell to real cell %i.\n", myid);
#endif
					} else if (cell_old->lb_cell_type == LB_EMPTY_CELL){
#ifdef DEBUG
						/*valid = 0;*/
						printf("Load Balance: Changed an empty cell to real cell %i.\n", myid);
#endif
					}

				}
			}
		}
	}

	/* Assign state of empty and buffer cells*/
	for (x=0; x<cell_dim.x; ++x){
		for (y=0; y<cell_dim.y; ++y){
			for (z=0; z<cell_dim.z; ++z){
				c = PTR_3D_V(cell_array, x, y, z, cell_dim);
				if (c->lb_cell_type != LB_REAL_CELL){
					c->lb_cell_type = lb_identifyCellType(x,y,z);
					c->lb_neighbor_index = -c->lb_cell_type;
					c->lb_cpu_affinity = -1;
				}
			}
		}
	}

	/*Create the new communication lists*/
	/*Includes further tests if new geometry is requiring impossible communication*/
	if (valid)
		valid = lb_syncBufferCellAffinity();


	/****************************************************************/
	/* Test if all geometries are valid, if not perform rollback.   */
	/* Furthermore, disable CPUs causing problems and try once more */
	/****************************************************************/
	int allValid;
	MPI_Allreduce(&valid, &allValid, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);
	if (!allValid){
		if (!valid) /* Disable movement of corners in the next cycle*/
			lb_disabledAtThisCPU = 1;

		/*Perform rollpack*/
		for (x=0; x<cell_dim.x; ++x){
			for (y=0; y<cell_dim.y; ++y){
				for (z=0; z<cell_dim.z; ++z){
					c = PTR_3D_V(cell_array, x, y, z, cell_dim);
					alloc_cell( c, 0 );
				}
			}
		}
		free(cell_array);

		cell_array = cell_array_old;

		lb_cell_offset.x = lb_cell_offset_old.x;
		lb_cell_offset.y = lb_cell_offset_old.y;
		lb_cell_offset.z = lb_cell_offset_old.z;
		cell_dim.x = cell_dim_old.x;
		cell_dim.y = cell_dim_old.y;
		cell_dim.z = cell_dim_old.z;
		cellmin.x = 1; cellmin.y = 1; cellmin.z = 1;
		cellmax.x = cell_dim.x-1; cellmax.y = cell_dim.y-1; cellmax.z = cell_dim.z-1;

		for (i=0; i<8;i++){
			lb_domain.corners[i].p.x = oldPositions[i].x;
			lb_domain.corners[i].p.y = oldPositions[i].y;
			lb_domain.corners[i].p.z = oldPositions[i].z;
		}
		lb_updateDomain(&lb_domain);
#ifdef DEBUG
		if (valid == 0)
			printf("CPU %i caused load balancing to be reverted on step %i\n",myid, steps);
#endif
		/* Try again with disabled movement of corners at CPUs that rejected the step */
		if (cycle == 0){
			int success = balanceLoad(lb_getLoad, lb_getCenterOfGravity, reset, 1);
			if (success) return 1;
			else if (myid == 0)
				printf("Load Balance step was globally rejected at step %i\n", steps);
		}

		lb_syncBufferCellAffinity();
		make_cell_lists();
		setup_buffers();
#ifdef NBLIST
		lb_need_nbl_update = 1;
		have_valid_nbl = 0;
#endif

		return 0;	/* Balancing was rejected twice in a row, give up and continue with the last valid configuration*/
	} else {
		/* The step was globally accepted, all CPUs are allowed to participate again*/
		lb_disabledAtThisCPU = 0;
	}


	/*Exchange particles between cpus*/
	lb_relocateParticles(lb_cell_offset_old, cell_dim_old, cell_array_old);

	/*Dealloc old cells*/
	for (x=0; x<cell_dim_old.x; ++x){
		for (y=0; y<cell_dim_old.y; ++y){
			for (z=0; z<cell_dim_old.z; ++z){
				c = PTR_3D_V(cell_array_old, x, y, z, cell_dim_old);
				alloc_cell( c, 0 );
			}
		}
	}
	free(cell_array_old);

	make_cell_lists();


#ifdef NBLIST
	lb_need_nbl_update = 1;
	have_valid_nbl = 0;
#endif

	/*Check*/
#ifdef DEBUG
	atomsAfter = lb_countAtoms();
	MPI_Reduce(&atomsBefore, &totalAtomsBefore, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&atomsAfter, &totalsAtomsAfter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (myid==0){
		if (totalAtomsBefore != totalsAtomsAfter)
			error("LB: Load balance corrupt: Number of atoms has changed");
	}
#endif

	setup_buffers();

	return 1;
}

void lb_moveAllCorners(real *domainInfo, int iteration){
	int i,j,x,y,z, rank, partners, master;
	vektor cog, p, cpuCoord;
	real localDens[8], force[3];
	vektor u[8];
	real tmp;
	real domainInfos[64];
	real averageDens;
	MPI_Status status;

	x = iteration&1;
	y = (iteration&2)>>1;
	z = (iteration&4)>>2;

	/* In order to avoid biasing if two corners of a domain move at the same time,
	 * only every second node (or 1/8 of the total nodes) is moved at a time.
	 * Depending on the index of the current cpu, identify the index of the corner on
	 * the cpu that is the handler of the corner (the one where the moving corner is on position 0)
	 */
	int com = ((my_coord.x+x) & 1) + ((my_coord.y+y) & 1)*2 + ((my_coord.z+z) & 1)*4;

	if (lb_localCommPartners[com][0] == 0) return;

	rank = lb_localCommPartners[com][1];
	partners = lb_localCommPartners[com][0];
	master = lb_localCommPartners[com][2];

	/* Receive the domain of all processes in the local group*/
	if (myid == master){
		for (i=0;i<8;i++)
			domainInfos[i] = domainInfo[i];
		for (i=1; i<partners;i++){
			int partner = lb_localCommPartners[com][i+2];
			MPI_Recv(&domainInfos[i*8], 8, REAL, partner, partner, MPI_COMM_WORLD, &status);
		}
	} else {
		MPI_Send(domainInfo, 8, REAL, master, myid, MPI_COMM_WORLD);
	}

	if (rank == 0){
		//Compute force
		averageDens = 0.;
		for (j=0;j<partners;j++){
			cog.x = domainInfos[j*8+2];
			cog.y = domainInfos[j*8+3];
			cog.z = domainInfos[j*8+4];

			cpuCoord.x = domainInfos[j*8+5];
			cpuCoord.y = domainInfos[j*8+6];
			cpuCoord.z = domainInfos[j*8+7];
			/*Checking periodicity and wrap around if necessary*/
			if (cpuCoord.x > my_coord.x) cog.x -= (box_x.x + box_y.x + box_z.x);
			if (cpuCoord.y > my_coord.y) cog.y -= (box_x.y + box_y.y + box_z.y);
			if (cpuCoord.z > my_coord.z) cog.z -= (box_x.z + box_y.z + box_z.z);

			p.x = lb_domain.corners[0].p.x;
			p.y = lb_domain.corners[0].p.y;
			p.z = lb_domain.corners[0].p.z;
			tmp = 1./SQRT( (cog.x-p.x)*(cog.x-p.x)+(cog.y-p.y)*(cog.y-p.y)+(cog.z-p.z)*(cog.z-p.z)  );
			u[j].x = (cog.x-p.x)*tmp;
			u[j].y = (cog.y-p.y)*tmp;
			u[j].z = (cog.z-p.z)*tmp;

			localDens[j] = domainInfos[j*8];
			averageDens += localDens[j];
		}
		if (averageDens == 0.){
			real totalBoxVolume = (box_x.x + box_y.x + box_z.x) *
					(box_x.y + box_y.y + box_z.y) *
					(box_x.z + box_y.z + box_z.z);
			for (j=0; j<partners; j++){
				localDens[j] = domainInfos[j*8+1]/totalBoxVolume;
				localDens[j] *= cpu_dim.x*cpu_dim.y*cpu_dim.z;
				averageDens += localDens[j];
			}
		}

		averageDens /= partners;
		force[0] = 0.; force[1] = 0.; force[2] = 0.;

		for (j=0; j<partners; j++){
			force[0] += u[j].x*(localDens[j]-averageDens);
			force[1] += u[j].y*(localDens[j]-averageDens);
			force[2] += u[j].z*(localDens[j]-averageDens);
		}

		real var = 0.;
		for (j=0; j<partners; j++){
			var += (localDens[j]-averageDens)*(localDens[j]-averageDens);
		}
		var = SQRT(var/partners);

		force[0] *= lb_contractionRate*(1+var);
		force[1] *= lb_contractionRate*(1+var);
		force[2] *= lb_contractionRate*(1+var);
		force[0] = lb_domain.corners[com].fixed.x ? 0.: force[0];
		force[1] = lb_domain.corners[com].fixed.y ? 0.: force[1];
		force[2] = lb_domain.corners[com].fixed.z ? 0.: force[2];

		real l = SQRT(force[0]*force[0] + force[1]*force[1] + force[2]*force[2]);
		real maxMove = MIN(MIN(lb_cell_size.x, lb_cell_size.y), lb_cell_size.z);
		if (l>maxMove){
			force[0] = force[0] * maxMove/l;
			force[1] = force[1] * maxMove/l;
			force[2] = force[2] * maxMove/l;
		}
	}

	/* broadcast force to all processes in the local group*/
	if (myid == master){
		for (i=1; i<partners;i++){
			int partner = lb_localCommPartners[com][i+2];
			MPI_Send(force, 3, REAL, partner, myid, MPI_COMM_WORLD);
		}
	} else {
		MPI_Recv(force, 3, REAL, master, master, MPI_COMM_WORLD, &status);
	}
	
	/*Try to move node  */
	if (force[0]*force[0] + force[1]*force[1] + force[2]*force[2] > 0.001){
		/*Test first what is the smallest fraction of the move that satisfies the geometrical constraints*/
		int smallestValidStep = 0;

		if (!lb_disabledAtThisCPU){
			vektor defaultPos;
			defaultPos.x = lb_domain.corners[com].p.x;
			defaultPos.y = lb_domain.corners[com].p.y;
			defaultPos.z = lb_domain.corners[com].p.z;
			for (i=1; i<4; i++){
				lb_moveCorner(com, force[0]/(1<<(3-i)), force[1]/(1<<(3-i)), force[2]/(1<<(3-i)));
				lb_updateDomain(&lb_domain);
				int valid = lb_isGeometryChangeValid();
				/* Reset corner */
				lb_domain.corners[com].p.x = defaultPos.x;
				lb_domain.corners[com].p.y = defaultPos.y;
				lb_domain.corners[com].p.z = defaultPos.z;

				if (valid)
					smallestValidStep=i;
				else break;
			}
		}

		int smallestCommonValid;

		/* Identify the smallest step size, valid in all domains */
		if (myid == master){
			smallestCommonValid = smallestValidStep;
			for (i=1; i<partners;i++){
				int partner = lb_localCommPartners[com][i+2];
				MPI_Recv(&smallestValidStep, 1, MPI_INT, partner, partner, MPI_COMM_WORLD, &status);
				smallestCommonValid = MIN(smallestCommonValid, smallestValidStep);
			}
			for (i=1; i<partners;i++){
				int partner = lb_localCommPartners[com][i+2];
				MPI_Send(&smallestCommonValid, 1, MPI_INT, partner, master, MPI_COMM_WORLD);
			}
		} else {
			MPI_Send(&smallestValidStep, 1, MPI_INT, master, myid, MPI_COMM_WORLD);
			MPI_Recv(&smallestCommonValid, 1, MPI_INT, master, master, MPI_COMM_WORLD, &status);
		}

		if (smallestCommonValid != 0){
			real fraction = 1<<(3-smallestCommonValid);
			lb_moveCorner(com, force[0]/fraction, force[1]/fraction, force[2]/fraction);
		}
		lb_updateDomain(&lb_domain);
	}
}

void lb_moveCornersReset(real *domainInfo, int iteration){
	int i,x,y,z, rank, partners, master;
	real force[3];
	MPI_Status status;

	x = iteration&1;
	y = (iteration&2)>>1;
	z = (iteration&4)>>2;

	/* In order to avoid biasing if two corners of a domain move at the same time,
	 * only every second node (or 1/8 of the total nodes) is moved at a time.
	 * Depending on the index of the current cpu, identify the index of the corner on
	 * the cpu that is the handler of the corner (the one where the moving corner is on position 0)
	 */
	int com = ((my_coord.x+x) & 1) + ((my_coord.y+y) & 1)*2 + ((my_coord.z+z) & 1)*4;

	if (lb_localCommPartners[com][0] == 0) return;

	rank = lb_localCommPartners[com][1];
	partners = lb_localCommPartners[com][0];
	master = lb_localCommPartners[com][2];

	if (rank == 0){
		/*Move node in the direction of its position on an regular grid*/
		force[0] = lb_domain.corners[com].ref.x - lb_domain.corners[com].p.x ;
		force[1] = lb_domain.corners[com].ref.y - lb_domain.corners[com].p.y ;
		force[2] = lb_domain.corners[com].ref.z - lb_domain.corners[com].p.z ;

		force[0] = lb_domain.corners[com].fixed.x ? 0.: force[0];
		force[1] = lb_domain.corners[com].fixed.y ? 0.: force[1];
		force[2] = lb_domain.corners[com].fixed.z ? 0.: force[2];

		if (ABS(force[0]) > lb_cell_size.x) force[0] = lb_cell_size.x*SIGNUM(force[0]);
		if (ABS(force[1]) > lb_cell_size.y) force[1] = lb_cell_size.y*SIGNUM(force[1]);
		if (ABS(force[2]) > lb_cell_size.z) force[2] = lb_cell_size.z*SIGNUM(force[2]);

		/* Limit length of vector*/
		real l = SQRT(force[0]*force[0] + force[1]*force[1] + force[2]*force[2]);
		real maxMove = MIN(MIN(lb_cell_size.x, lb_cell_size.y), lb_cell_size.z)*0.9;
		if (l>maxMove){
			force[0] = force[0] * maxMove/l;
			force[1] = force[1] * maxMove/l;
			force[2] = force[2] * maxMove/l;
		}
	}

	/* broadcast force to all processes in the local group*/
	if (myid == master){
		for (i=1; i<partners;i++){
			int partner = lb_localCommPartners[com][i+2];
			MPI_Send(force, 3, REAL, partner, myid, MPI_COMM_WORLD);
		}
	} else {
		MPI_Recv(force, 8, REAL, master, master, MPI_COMM_WORLD, &status);
	}

	/*Try to move node  */
	if (force[0]*force[0] + force[1]*force[1] + force[2]*force[2] > 0.001){
		/*Test first what is the smallest fraction of the move that satisfies the geometrical constraints*/
		int smallestValidStep = 0;

		if (!lb_disabledAtThisCPU){
			vektor defaultPos;
			defaultPos.x = lb_domain.corners[com].p.x;
			defaultPos.y = lb_domain.corners[com].p.y;
			defaultPos.z = lb_domain.corners[com].p.z;
			for (i=1; i<4; i++){
				lb_moveCorner(com, force[0]/(1<<(3-i)), force[1]/(1<<(3-i)), force[2]/(1<<(3-i)));
				lb_updateDomain(&lb_domain);
				int valid = lb_isGeometryChangeValid();
				/* Reset corner */
				lb_domain.corners[com].p.x = defaultPos.x;
				lb_domain.corners[com].p.y = defaultPos.y;
				lb_domain.corners[com].p.z = defaultPos.z;

				if (valid)
					smallestValidStep=i;
				else break;
			}
		}

		int smallestCommonValid;

		/* Identify the smallest step size, valid in all domains */
		if (myid == master){
			smallestCommonValid = smallestValidStep;
			for (i=1; i<partners;i++){
				int partner = lb_localCommPartners[com][i+2];
				MPI_Recv(&smallestValidStep, 1, MPI_INT, partner, partner, MPI_COMM_WORLD, &status);
				smallestCommonValid = MIN(smallestCommonValid, smallestValidStep);
			}
			for (i=1; i<partners;i++){
				int partner = lb_localCommPartners[com][i+2];
				MPI_Send(&smallestCommonValid, 1, MPI_INT, partner, master, MPI_COMM_WORLD);
			}
		} else {
			MPI_Send(&smallestValidStep, 1, MPI_INT, master, myid, MPI_COMM_WORLD);
			MPI_Recv(&smallestCommonValid, 1, MPI_INT, master, master, MPI_COMM_WORLD, &status);
		}

		if (smallestCommonValid != 0){
			real fraction = 1<<(3-smallestCommonValid);
			lb_moveCorner(com, force[0]/fraction, force[1]/fraction, force[2]/fraction);
		}
		lb_updateDomain(&lb_domain);
	}
}

int lb_identifyCellType(int x, int y, int z){
	cell *c;
	c = PTR_3D_V(cell_array, x, y, z, cell_dim);
	int i,j,k;
	/* Assign cell next to a real cell as buffer, otherwise as empty */

	if ( (pbc_dirs.x == 0 && x+lb_cell_offset.x < 0) ||
			(pbc_dirs.x == 0 && x+lb_cell_offset.x > global_cell_dim.x - 1) ||
			(pbc_dirs.y == 0 && y+lb_cell_offset.y < 0) ||
			(pbc_dirs.y == 0 && y+lb_cell_offset.y > global_cell_dim.y - 1) ||
			(pbc_dirs.z == 0 && z+lb_cell_offset.z < 0) ||
			(pbc_dirs.z == 0 && z+lb_cell_offset.z > global_cell_dim.z - 1)){
		return LB_NON_PBC_BUFFER_CELL;
	}

	for (i = -1; i <= 1; i++)
		for (j = -1; j <= 1; j++)
			for (k = -1; k <= 1; k++)
				if ( (c = lb_accessCell(cell_array,x-i, y-j, z-k, cell_dim)) != NULL
						&& c->lb_cell_type == LB_REAL_CELL)
					return LB_BUFFER_CELL;

	return LB_EMPTY_CELL;
}

cell* lb_accessCell(cell* cellArray, int x, int y, int z, ivektor cellDim){
	if (x<0 || y<0 || z<0) return NULL;
	if (x>=cellDim.x || y>=cellDim.y || z>=cellDim.z) return NULL;
	cell *c = PTR_3D_V(cellArray,x, y, z, cellDim);
	return c;
}

real lb_getLoad(){
	real load = lb_countAtoms();
	/*Scale load, in case of an evenly distribution, each cpu has a load of exactly 1.*/
	load *= (cpu_dim.x*cpu_dim.y*cpu_dim.z) / (real)natoms;
	return load;
}

void lb_makeNormals(lb_domainInfo *dom){
	lb_makeNormal(0 , 0, 1, 5, dom);
	lb_makeNormal(1 , 0, 5, 4, dom);
	lb_makeNormal(2 , 7, 3, 2, dom);
	lb_makeNormal(3 , 6, 7, 2, dom);
	lb_makeNormal(4 , 1, 3, 7, dom);
	lb_makeNormal(5 , 1, 7, 5, dom);
	lb_makeNormal(6 , 6, 2, 0, dom);
	lb_makeNormal(7 , 4, 6, 0, dom);
	lb_makeNormal(8 , 4, 5, 7, dom);
	lb_makeNormal(9 , 4, 7, 6, dom);
	lb_makeNormal(10, 3, 1, 0, dom);
	lb_makeNormal(11, 2, 3, 0, dom);

	dom->faceConvexity[0] = lb_getTetraederVolumeIndexed(0, 1, 5, 4, dom);
	dom->faceConvexity[1] = lb_getTetraederVolumeIndexed(3, 2, 6, 7, dom);
	dom->faceConvexity[2] = lb_getTetraederVolumeIndexed(1, 3, 7, 5, dom);
	dom->faceConvexity[3] = lb_getTetraederVolumeIndexed(2, 0, 4, 6, dom);
	dom->faceConvexity[4] = lb_getTetraederVolumeIndexed(4, 5, 7, 6, dom);
	dom->faceConvexity[5] = lb_getTetraederVolumeIndexed(2, 3, 1, 0, dom);
}

void lb_makeNormal(int n, int p0, int p1, int p2, lb_domainInfo *dom){
	dom->normals[n].x = (dom->corners[p1].index.y - dom->corners[p0].index.y) *
	       (dom->corners[p2].index.z - dom->corners[p0].index.z) -
	       (dom->corners[p1].index.z - dom->corners[p0].index.z) *
	       (dom->corners[p2].index.y - dom->corners[p0].index.y);
	dom->normals[n].y = (dom->corners[p1].index.z - dom->corners[p0].index.z) *
	       (dom->corners[p2].index.x - dom->corners[p0].index.x) -
	       (dom->corners[p1].index.x - dom->corners[p0].index.x) *
	       (dom->corners[p2].index.z - dom->corners[p0].index.z);
	dom->normals[n].z = (dom->corners[p1].index.x - dom->corners[p0].index.x) *
	  	   (dom->corners[p2].index.y - dom->corners[p0].index.y) -
	  	   (dom->corners[p1].index.y - dom->corners[p0].index.y) *
	  	   (dom->corners[p2].index.x - dom->corners[p0].index.x);
}

int lb_getTetraederVolumeIndexed(int c1, int c2, int c3, int c4, lb_domainInfo *dom) {
	int dir1_0, dir1_1, dir1_2;
	int dir2_0, dir2_1, dir2_2;
	int dir3_0, dir3_1, dir3_2;

	dir1_0 = dom->corners[c2].index.x - dom->corners[c1].index.x;
	dir1_1 = dom->corners[c2].index.y - dom->corners[c1].index.y;
	dir1_2 = dom->corners[c2].index.z - dom->corners[c1].index.z;

	dir2_0 = dom->corners[c3].index.x - dom->corners[c1].index.x;
	dir2_1 = dom->corners[c3].index.y - dom->corners[c1].index.y;
	dir2_2 = dom->corners[c3].index.z - dom->corners[c1].index.z;

	dir3_0 = dom->corners[c4].index.x - dom->corners[c1].index.x;
	dir3_1 = dom->corners[c4].index.y - dom->corners[c1].index.y;
	dir3_2 = dom->corners[c4].index.z - dom->corners[c1].index.z;

	return (dir1_0 * (dir3_1 * dir2_2 - dir3_2 * dir2_1) +
			dir1_1 * (dir3_2 * dir2_0 - dir3_0 * dir2_2) +
			dir1_2 * (dir3_0 * dir2_1 - dir3_1 * dir2_0));
}

real lb_getTetrahedronVolume(int c1, int c2, int c3, int c4, lb_domainInfo *dom){
	real dir1_0, dir1_1, dir1_2;
	real dir2_0, dir2_1, dir2_2;
	real dir3_0, dir3_1, dir3_2;

	dir1_0 = dom->corners[c2].discretizedP.x - dom->corners[c1].discretizedP.x;
	dir1_1 = dom->corners[c2].discretizedP.y - dom->corners[c1].discretizedP.y;
	dir1_2 = dom->corners[c2].discretizedP.z - dom->corners[c1].discretizedP.z;

	dir2_0 = dom->corners[c3].discretizedP.x - dom->corners[c1].discretizedP.x;
	dir2_1 = dom->corners[c3].discretizedP.y - dom->corners[c1].discretizedP.y;
	dir2_2 = dom->corners[c3].discretizedP.z - dom->corners[c1].discretizedP.z;

	dir3_0 = dom->corners[c4].discretizedP.x - dom->corners[c1].discretizedP.x;
	dir3_1 = dom->corners[c4].discretizedP.y - dom->corners[c1].discretizedP.y;
	dir3_2 = dom->corners[c4].discretizedP.z - dom->corners[c1].discretizedP.z;

	return (dir1_0 * (dir3_1 * dir2_2 - dir3_2 * dir2_1) +
			dir1_1 * (dir3_2 * dir2_0 - dir3_0 * dir2_2) +
			dir1_2 * (dir3_0 * dir2_1 - dir3_1 * dir2_0))/6.;
}

/**
 * Returns the center of a cell on coordinate x,y,z in the global cell grid
 */
vektor lb_getCellCenter(int x, int y, int z){
	vektor v;
	v.x = (x+0.5)*lb_cell_size.x;
	v.y = (y+0.5)*lb_cell_size.y;
	v.z = (z+0.5)*lb_cell_size.z;

	return v;
}



real lb_angleCos(int u1, int c, int v1, lb_domainInfo *dom){
	real u_0 = dom->corners[u1].discretizedP.x - dom->corners[c].discretizedP.x;
	real u_1 = dom->corners[u1].discretizedP.y - dom->corners[c].discretizedP.y;
	real u_2 = dom->corners[u1].discretizedP.z - dom->corners[c].discretizedP.z;

	real v_0 = dom->corners[v1].discretizedP.x - dom->corners[c].discretizedP.x;
	real v_1 = dom->corners[v1].discretizedP.y - dom->corners[c].discretizedP.y;
	real v_2 = dom->corners[v1].discretizedP.z - dom->corners[c].discretizedP.z;

	real a = (u_0*v_0+u_1*v_1+u_2*v_2) /(SQRT(u_0*u_0+u_1*u_1+u_2*u_2) * SQRT(v_0*v_0+v_1*v_1+v_2*v_2));
	if (a>1) a = 1;		//Rounding errors produce in some cases something like 1.00000001
	if (a<-1) a =- 1;
	return a;
}

real lb_computeWarpage(int p0, int p1, int p2, int p, lb_domainInfo *dom){
	real distFromPlane, sinus;
	vektor v;

	distFromPlane = lb_getDistanceToPlane(p, p0, p1, p2, dom);

	v.x = dom->corners[p0].discretizedP.x - dom->corners[p].discretizedP.x;
	v.y = dom->corners[p0].discretizedP.y - dom->corners[p].discretizedP.y;
	v.z = dom->corners[p0].discretizedP.z - dom->corners[p].discretizedP.z;

	sinus = distFromPlane / SQRT(v.x*v.x + v.y*v.y + v.z*v.z);
	return sinus;
}

/*
 * Returns the squared distance from the domain corner with index p
 * to the plane formed by the domain corners with index p0,p1,p2
 */
real lb_getDistanceToPlane(int p, int p0, int p1, int p2, lb_domainInfo *dom){
	real d, d2, dist;
	vektor v;
	/* Normal of the plane*/
	v.x = (dom->corners[p1].discretizedP.y - dom->corners[p0].discretizedP.y) *
		  (dom->corners[p2].discretizedP.z - dom->corners[p0].discretizedP.z) -
		  (dom->corners[p1].discretizedP.z - dom->corners[p0].discretizedP.z) *
		  (dom->corners[p2].discretizedP.y - dom->corners[p0].discretizedP.y);
	v.y = (dom->corners[p1].discretizedP.z - dom->corners[p0].discretizedP.z) *
		  (dom->corners[p2].discretizedP.x - dom->corners[p0].discretizedP.x) -
		  (dom->corners[p1].discretizedP.x - dom->corners[p0].discretizedP.x) *
		  (dom->corners[p2].discretizedP.z - dom->corners[p0].discretizedP.z);
	v.z = (dom->corners[p1].discretizedP.x - dom->corners[p0].discretizedP.x) *
		  (dom->corners[p2].discretizedP.y - dom->corners[p0].discretizedP.y) -
		  (dom->corners[p1].discretizedP.y - dom->corners[p0].discretizedP.y) *
		  (dom->corners[p2].discretizedP.x - dom->corners[p0].discretizedP.x);
	/* setup equations of two parallel planes (v*p0=d, v*p=d2) */
	d  = dom->corners[p0].discretizedP.x * v.x +
		 dom->corners[p0].discretizedP.y * v.y +
		 dom->corners[p0].discretizedP.z * v.z;
	d2 = dom->corners[p].discretizedP.x * v.x +
		 dom->corners[p].discretizedP.y * v.y +
		 dom->corners[p].discretizedP.z * v.z;
	/*distance between the planes*/
	dist = abs(d2-d) / SQRT(v.x*v.x + v.y*v.y + v.z*v.z);
	return dist;
}

/**
 * Tests if a given point p is inside the domain or not
 * Returns: 1 - point is inside, 0 otherwise
 */
int lb_isPointInDomain(vektor p, lb_domainInfo *dom){
	ivektor pi;
	ivektor *v;

	/* Discretize the point onto the sublevel grid
	 * Integer arithmetic makes sure that a point is not accepted by two domains if the point
	 * is placed directly on the boundary. In that case one domain will reject and the other will accept
	 * by construction
	 */
	pi.x = (int)((p.x/ lb_cell_size.x )*LB_CELL_SUBLEVELS);
	pi.y = (int)((p.y/ lb_cell_size.y )*LB_CELL_SUBLEVELS);
	pi.z = (int)((p.z/ lb_cell_size.z )*LB_CELL_SUBLEVELS);

	/* Test the six faces of the initial cube*/
	/* Each initial face is splitted into two triangles*/
	/* Compute the product of the triangles normal vector and the vector to pi*/
	v = &dom->normals[0];
	int d1 = (pi.x-dom->corners[0].index.x)*v->x+(pi.y-dom->corners[0].index.y)*v->y+(pi.z-dom->corners[0].index.z)*v->z;
	v = &dom->normals[1];
	int d2 = (pi.x-dom->corners[0].index.x)*v->x+(pi.y-dom->corners[0].index.y)*v->y+(pi.z-dom->corners[0].index.z)*v->z;

	/* There are two cases. The two triangles on the face can be either convex or concave.
	 * A positive volume of the domain tetrahedron indicates convexity. In this case on normal product must be positive
	 * and the point p is not in the domain. In case of a concave face, both normals must be larger than zero.
	 * If the products are both zero, the point is not accepted as well, but the domain on the other side of the face will accept.
	 * (There the normals are tested for >=0)
	 */
	if (dom->faceConvexity[0]<= 0){
		if (d1>0 && d2>0) return 0;
	} else if (d1>0 || d2>0) return 0;

	/*The other faces are tested similarly*/
	v = &dom->normals[2];
	d1 = (pi.x-dom->corners[2].index.x)*v->x+(pi.y-dom->corners[2].index.y)*v->y+(pi.z-dom->corners[2].index.z)*v->z;
	v = &dom->normals[3];
	d2 = (pi.x-dom->corners[2].index.x)*v->x+(pi.y-dom->corners[2].index.y)*v->y+(pi.z-dom->corners[2].index.z)*v->z;
	if (dom->faceConvexity[1]>= 0){
		if (d1>=0 && d2>=0) return 0;
	} else if (d1>=0 || d2>=0) return 0;

	v = &dom->normals[4];
	d1 = (pi.x-dom->corners[1].index.x)*v->x+(pi.y-dom->corners[1].index.y)*v->y+(pi.z-dom->corners[1].index.z)*v->z;
	v = &dom->normals[5];
	d2 = (pi.x-dom->corners[1].index.x)*v->x+(pi.y-dom->corners[1].index.y)*v->y+(pi.z-dom->corners[1].index.z)*v->z;
	if (dom->faceConvexity[2]<= 0){
		if (d1>=0 && d2>=0) return 0;
	} else if (d1>=0 || d2>=0) return 0;

	v = &dom->normals[6];
	d1 = (pi.x-dom->corners[0].index.x)*v->x+(pi.y-dom->corners[0].index.y)*v->y+(pi.z-dom->corners[0].index.z)*v->z;
	v = &dom->normals[7];
	d2 = (pi.x-dom->corners[0].index.x)*v->x+(pi.y-dom->corners[0].index.y)*v->y+(pi.z-dom->corners[0].index.z)*v->z;
	if (dom->faceConvexity[3]>= 0){
		if (d1>0 && d2>0) return 0;
	} else if (d1>0 || d2>0) return 0;


	v = &dom->normals[8];
	d1 = (pi.x-dom->corners[4].index.x)*v->x+(pi.y-dom->corners[4].index.y)*v->y+(pi.z-dom->corners[4].index.z)*v->z;
	v = &dom->normals[9];
	d2 = (pi.x-dom->corners[4].index.x)*v->x+(pi.y-dom->corners[4].index.y)*v->y+(pi.z-dom->corners[4].index.z)*v->z;
	if (dom->faceConvexity[4]<= 0){
		if (d1>=0 && d2>=0) return 0;
	} else if (d1>=0 || d2>=0) return 0;

	v = &dom->normals[10];
	d1 = (pi.x-dom->corners[0].index.x)*v->x+(pi.y-dom->corners[0].index.y)*v->y+(pi.z-dom->corners[0].index.z)*v->z;
	v = &dom->normals[11];
	d2 = (pi.x-dom->corners[0].index.x)*v->x+(pi.y-dom->corners[0].index.y)*v->y+(pi.z-dom->corners[0].index.z)*v->z;
	if (dom->faceConvexity[5]>= 0){
		if (d1>0 && d2>0) return 0;
	} else if (d1>0 || d2>0) return 0;

	return 1;
}

/**
 * Return the center of gravity of this domain
 * If there are particles inside the domain, the cog of all particles is returned
 * If there are no particles, the cog of all real cells is returned instead
 */
vektor lb_getCenterOfGravity(){
	vektor cog, center;
	int x,y,z,i, n = 0;
	cell *c;
	cog.x = 0.; cog.y = 0.; cog.z = 0.;

	/*Compute center of gravity of all particles in real cells on this CPU*/
	for (x=0; x<nallcells;++x){
		c = cell_array+x;
		if (c->lb_cell_type == LB_REAL_CELL){
			for (i=0; i<c->n; i++){
				cog.x += ORT(c,i,X);
				cog.y += ORT(c,i,Y);
				cog.z += ORT(c,i,Z);
			}
			n += c->n;
		}
	}

	if (n == 0){
		/*There are no particles on this CPU, compute the center of gravity of the domain's volume*/
		for (x=1; x<cell_dim.x-1; ++x){
			for (y=1; y<cell_dim.y-1; ++y){
				for (z=1; z<cell_dim.z-1; ++z){
					c = PTR_3D_V(cell_array, x, y, z, cell_dim);
					if (c->lb_cell_type == LB_REAL_CELL){
						center = lb_getCellCenter(x+lb_cell_offset.x,y+lb_cell_offset.y,z+lb_cell_offset.z);
						cog.x += center.x;
						cog.y += center.y;
						cog.z += center.z;
						n++;
					}
				}
			}
		}
	}

	cog.x /= n;
	cog.y /= n;
	cog.z /= n;

	return cog;
}

/*
 * Compute and cache additional values needed later to check which cells belong to this CPU,
 * based on the position of the eight domain corners
 */
void lb_updateDomain(lb_domainInfo *dom){
	int i;
	for (i = 0; i < 8; i++) {
		/* Transform position into an integer index, required for integer arithmetic in order to identify if a
		 * cell is inside the domain without the possibility of rounding errors.
		 * Otherwise there is the chance that if the center of a cell is very close to the domain interface
		 * that either both or no CPU will handle the cell. Both cases are fatal and must be avoided.
		 * Each cell is subdivided into a certain number of levels. The index is a discretization onto one level.
		 */
		dom->corners[i].index.x = (int) ((dom->corners[i].p.x / lb_cell_size.x) * LB_CELL_SUBLEVELS);
		dom->corners[i].index.y = (int) ((dom->corners[i].p.y / lb_cell_size.y) * LB_CELL_SUBLEVELS);
		dom->corners[i].index.z = (int) ((dom->corners[i].p.z / lb_cell_size.z) * LB_CELL_SUBLEVELS);

		/*
		 * Retransform the integer coordinate back into the world coordinates
		 */
		dom->corners[i].discretizedP.x = dom->corners[i].index.x * lb_cell_size.x / LB_CELL_SUBLEVELS;
		dom->corners[i].discretizedP.y = dom->corners[i].index.y * lb_cell_size.y / LB_CELL_SUBLEVELS;
		dom->corners[i].discretizedP.z = dom->corners[i].index.z * lb_cell_size.z / LB_CELL_SUBLEVELS;
	}
	/*Caching the normals of the 12 faces, defining the domain boundary*/
	lb_makeNormals(dom);
}

/**
 * Move the corner i (=[0..7]) by the given values dx,dy,dz
 * If necessary, the movement is limited to the size of one cell
 * If the corner is fixed in a direction, the movement in that direction will be 0.
 */
void lb_moveCorner(int i, real dx, real dy, real dz) {
	/*Limit movement velocity to cell size*/
	if (ABS(dx) > lb_cell_size.x) dx = lb_cell_size.x * SIGNUM(dx);
	if (ABS(dy) > lb_cell_size.y) dy = lb_cell_size.y * SIGNUM(dy);
	if (ABS(dz) > lb_cell_size.z) dz = lb_cell_size.z * SIGNUM(dz);
	/*Moving the point, if not restricted*/
	lb_domain.corners[i].p.x += lb_domain.corners[i].fixed.x ? 0. : dx;
	lb_domain.corners[i].p.y += lb_domain.corners[i].fixed.y ? 0. : dy;
	lb_domain.corners[i].p.z += lb_domain.corners[i].fixed.z ? 0. : dz;
}

/*
 * Returns the volume of the domain
 */
real lb_getVolume(){
	real volume = 0.;
	/* Split the domain into six tetrahedrons */
	volume += lb_getTetrahedronVolume(0, 5, 4, 7, &lb_domain);
	volume += lb_getTetrahedronVolume(0, 3, 1, 7, &lb_domain);
	volume += lb_getTetrahedronVolume(0, 1, 5, 7, &lb_domain);
	volume += lb_getTetrahedronVolume(0, 4, 6, 7, &lb_domain);
	volume += lb_getTetrahedronVolume(0, 6, 2, 7, &lb_domain);
	volume += lb_getTetrahedronVolume(0, 2, 3, 7, &lb_domain);
	return volume;
}

/**
 * Return the total number of atoms inside real cells in this domain
 */
int lb_countAtoms(){
	cell *c;
	int i, n = 0;
	for (i=0; i<nallcells;++i){
		c = cell_array+i;
		if (c->lb_cell_type == LB_REAL_CELL)
			n += c->n;
	}
	return n;
}

void lb_computeVariance(){
	int i;
	real load = lb_getLoad();
	real *allLoads;
	allLoads = malloc(cpu_dim.x*cpu_dim.y*cpu_dim.z * sizeof *allLoads);

	MPI_Gather(&load, 1, REAL, allLoads, 1, REAL, 0, MPI_COMM_WORLD );
	if (myid == 0){
		lb_maxLoad = load;
		lb_maxLoadOnCPU = myid;
		lb_minLoad = load;
		lb_loadVariance = 0.;
		for (i = 0; i < cpu_dim.x * cpu_dim.y * cpu_dim.z; i++) {
			if (allLoads[i] > lb_maxLoad){
				lb_maxLoad = allLoads[i];
				lb_maxLoadOnCPU = i;
			}
			if (allLoads[i] < lb_minLoad) lb_minLoad = allLoads[i];
			lb_loadVariance += (allLoads[i] - 1.) * (allLoads[i] - 1.);
		}
		lb_loadVariance /= cpu_dim.x * cpu_dim.y * cpu_dim.z;
	}
	free(allLoads);

	MPI_Bcast(&lb_loadVariance, 1, REAL, 0, MPI_COMM_WORLD);
	MPI_Bcast(&lb_minLoad, 1, REAL, 0, MPI_COMM_WORLD);
	MPI_Bcast(&lb_maxLoad, 1, REAL, 0, MPI_COMM_WORLD);
}

/* Load balancing with orthogonal domains */
void balanceOrtho(){
	int i,x,y,z;

	int *x_count_local = malloc(global_cell_dim.x * sizeof *x_count_local);
	int *y_count_local = malloc(global_cell_dim.y * sizeof *y_count_local);
	int *z_count_local = malloc(global_cell_dim.z * sizeof *z_count_local);
	for (i = 0; i<global_cell_dim.x; i++)
		x_count_local[i] = 0;
	for (i = 0; i<global_cell_dim.y; i++)
		y_count_local[i] = 0;
	for (i = 0; i<global_cell_dim.z; i++)
		z_count_local[i] = 0;

	//count & reduce number of atoms in each cell layer
	for (x=1; x<cell_dim.x-1; ++x){
		for (y=1; y<cell_dim.y-1; ++y){
			for (z=1; z<cell_dim.z-1; ++z){
				cell *c = PTR_3D_V(cell_array, x, y, z, cell_dim);
				if (c->lb_cell_type == LB_REAL_CELL){
					x_count_local[x+lb_cell_offset.x]+=c->n;
					y_count_local[y+lb_cell_offset.y]+=c->n;
					z_count_local[z+lb_cell_offset.z]+=c->n;
				}
			}
		}
	}

	int *x_count = NULL ,*y_count = NULL, *z_count = NULL;
	if (myid==0){
		x_count = malloc(global_cell_dim.x * sizeof *x_count);
		y_count = malloc(global_cell_dim.y * sizeof *y_count);
		z_count = malloc(global_cell_dim.z * sizeof *z_count);
	}


	MPI_Reduce(x_count_local, x_count, global_cell_dim.x, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(y_count_local, y_count, global_cell_dim.y, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(z_count_local, z_count, global_cell_dim.z, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	free(x_count_local);
	free(y_count_local);
	free(z_count_local);

	if (myid==0){
		lb_balanceOneAxisOrthogonal(global_cell_dim.z, cpu_dim.z, z_count, z_bounds);
		lb_balanceOneAxisOrthogonal(global_cell_dim.y, cpu_dim.y, y_count, y_bounds);
		lb_balanceOneAxisOrthogonal(global_cell_dim.x, cpu_dim.x, x_count, x_bounds);

		free(x_count);
		free(y_count);
		free(z_count);
	}

	MPI_Bcast(x_bounds, cpu_dim.x+1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(y_bounds, cpu_dim.y+1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(z_bounds, cpu_dim.z+1, MPI_INT, 0, MPI_COMM_WORLD);
	/* Set new domain coordinates*/
	lb_domain.corners[0].p.x = lb_cell_size.x * (x_bounds[my_coord.x]);
	lb_domain.corners[0].p.y = lb_cell_size.y * (y_bounds[my_coord.y]);
	lb_domain.corners[0].p.z = lb_cell_size.z * (z_bounds[my_coord.z]);
	lb_domain.corners[1].p.x = lb_cell_size.x * (x_bounds[my_coord.x+1]);
	lb_domain.corners[1].p.y = lb_cell_size.y * (y_bounds[my_coord.y]);
	lb_domain.corners[1].p.z = lb_cell_size.z * (z_bounds[my_coord.z]);
	lb_domain.corners[2].p.x = lb_cell_size.x * (x_bounds[my_coord.x]);
	lb_domain.corners[2].p.y = lb_cell_size.y * (y_bounds[my_coord.y+1]);
	lb_domain.corners[2].p.z = lb_cell_size.z * (z_bounds[my_coord.z]);
	lb_domain.corners[3].p.x = lb_cell_size.x * (x_bounds[my_coord.x+1]);
	lb_domain.corners[3].p.y = lb_cell_size.y * (y_bounds[my_coord.y+1]);
	lb_domain.corners[3].p.z = lb_cell_size.z * (z_bounds[my_coord.z]);
	lb_domain.corners[4].p.x = lb_cell_size.x * (x_bounds[my_coord.x]);
	lb_domain.corners[4].p.y = lb_cell_size.y * (y_bounds[my_coord.y]);
	lb_domain.corners[4].p.z = lb_cell_size.z * (z_bounds[my_coord.z+1]);
	lb_domain.corners[5].p.x = lb_cell_size.x * (x_bounds[my_coord.x+1]);
	lb_domain.corners[5].p.y = lb_cell_size.y * (y_bounds[my_coord.y]);
	lb_domain.corners[5].p.z = lb_cell_size.z * (z_bounds[my_coord.z+1]);
	lb_domain.corners[6].p.x = lb_cell_size.x * (x_bounds[my_coord.x]);
	lb_domain.corners[6].p.y = lb_cell_size.y * (y_bounds[my_coord.y+1]);
	lb_domain.corners[6].p.z = lb_cell_size.z * (z_bounds[my_coord.z+1]);
	lb_domain.corners[7].p.x = lb_cell_size.x * (x_bounds[my_coord.x+1]);
	lb_domain.corners[7].p.y = lb_cell_size.y * (y_bounds[my_coord.y+1]);
	lb_domain.corners[7].p.z = lb_cell_size.z * (z_bounds[my_coord.z+1]);

	lb_updateDomain(&lb_domain);
}

void lb_balanceOneAxisOrthogonal(int numCellsInDirection, int numProcessorLayer, int* atomsPerCellLayer, int* bounds){
	int i;
	int remainingAtoms = 0;

	int *tmp_bounds = malloc((numProcessorLayer+1)*sizeof *tmp_bounds);

	for (i=0; i<=numProcessorLayer; i++)
		tmp_bounds[i] = bounds[i];

	for (i=0; i<numCellsInDirection; i++)
		remainingAtoms+=atomsPerCellLayer[i];

	for (i = 0; i<numProcessorLayer-1;i++){
		int currentNumberOfAtoms = 0;
		int targetNumberOfAtoms = remainingAtoms/(numProcessorLayer-i);
		int j;
		for (j=tmp_bounds[i]; j<tmp_bounds[i+1];j++)
			currentNumberOfAtoms += atomsPerCellLayer[j];
		/* Test if it better to move the boundary up, down or leave it at its current position*/
		int newLayerPosition = tmp_bounds[i+1];
		if (currentNumberOfAtoms < targetNumberOfAtoms){
			if (numCellsInDirection-tmp_bounds[i+1]>numProcessorLayer-i){	//Leave enough cells for the following layers
				int numberOfAtomsMovedUp = currentNumberOfAtoms+atomsPerCellLayer[tmp_bounds[i+1]+1];
				if (abs(targetNumberOfAtoms-numberOfAtomsMovedUp) < abs(targetNumberOfAtoms-currentNumberOfAtoms)){

					newLayerPosition++;
					remainingAtoms-=numberOfAtomsMovedUp;
				} else remainingAtoms -= currentNumberOfAtoms;
			} else remainingAtoms -= currentNumberOfAtoms;
		} else {
			if (tmp_bounds[i] != tmp_bounds[i+1]-1){	//Do not collapse the layer to zero
				int numberOfAtomsMovedDown = currentNumberOfAtoms-atomsPerCellLayer[tmp_bounds[i+1]-1];
				if (abs(targetNumberOfAtoms-numberOfAtomsMovedDown) < abs(targetNumberOfAtoms-currentNumberOfAtoms)){
					newLayerPosition--;
					remainingAtoms-=numberOfAtomsMovedDown;
				} else remainingAtoms -= currentNumberOfAtoms;
			} else remainingAtoms -= currentNumberOfAtoms;
		}

		tmp_bounds[i+1] = newLayerPosition;
	}

	for (i=0; i<=numProcessorLayer; i++){
		bounds[i] = tmp_bounds[i];
	}
	free(tmp_bounds);
}

