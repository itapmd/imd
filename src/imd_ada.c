/******************************************************************************
*
* imd_ada.c -- Routines for Bond angular distribution analysis
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/******************************************************************************
*
* Implements the scheme to distinguish atoms in different crystal structures
* published in Ackland & Jones PhysRevB 73 054104 (2006)
*
* Does not compute the the nearest neighbor distance dynamically, but with a
* fixed radius. For pure FCC and BCC materials, additional patterns to classify
* atoms are implemented.
*
******************************************************************************/

void init_ada(void) {
	/* if cutoff distance is not explicitly given, calculate it from the lattice constant*/
	if (ada_nbr_r2cut == 0.){
		if (ada_crystal_structure == ADA_FCC_CONFIG) {
			/*cutoff = ~1.22*nearest neighbor dist*/
			ada_nbr_r2cut = SQR(0.862*ada_latticeConst);
		} else if (ada_crystal_structure == ADA_BCC_CONFIG) {
			/*cutoff = ~1.22*second nearest neighbor dist*/
			ada_nbr_r2cut = SQR(1.22*ada_latticeConst);
		} else if (ada_crystal_structure == ADA_ACKLAND_CONFIG) {
			/*cutoff = ~1.22*nearest neighbor dist*/
			ada_nbr_r2cut = SQR(0.862*ada_latticeConst);
		}
	}

	if (ada_crystal_structure == ADA_FCC_CONFIG) {
		ada_default_type = 1;
	} else if (ada_crystal_structure == ADA_BCC_CONFIG) {
		ada_default_type = 0;
	} else if (ada_crystal_structure == ADA_ACKLAND_CONFIG) {
		ada_default_type = 1;
	} else {
		error("Crystal structure not supported by ADA");
	}

	if (myid == 0){
		if (ada_crystal_structure == ADA_FCC_CONFIG) {
			printf( "ADA: FCC analysis scheme selected \n");
		} else if (ada_crystal_structure == ADA_BCC_CONFIG) {
			printf( "ADA: BCC analysis scheme selected \n");
		} else if (ada_crystal_structure == ADA_ACKLAND_CONFIG) {
			printf( "ADA: ACKLAND analysis scheme selected \n");
		} else {
			error("Crystal structure not supported by ADA");
		}

		printf( "ADA: nearest neighbor cut-off radius: %f \n", sqrt(ada_nbr_r2cut));
		printf( "ADA: write interval: %d \n", ada_write_int);
	}
}

void do_ada(void) {
	real v1_length, v2_length, angle;
	int co_x0, co_x1, co_x2, co_x3, co_x4, co_x5, co_x6, co_x7;
	real delta_bcc, delta_cp, delta_fcc, delta_hcp;
	int nneigh, k, i, j, ii, jj;
	shortint type = 0;
	vektor v1, v2;

	do_neightab_complete();
	if (ada_crystal_structure == ADA_FCC_CONFIG) {
		/*
		 * typ=0: bcc
		 * typ=1: fcc
		 * typ=2: hcp
		 * typ=3: 12 neighbors, not fcc/hcp
		 * typ=4: less than 12 neighbors
		 * typ=5: more than 12 neighbors
		 * typ=6: less than 10 neighbors
		 * typ=7: more than 14 neighbors
		 */

		for (k = 0; k < ncells; k++) {
			cell *p = CELLPTR(k);
			for (i = 0; i < p->n; i++) {
				nneigh = NEIGH(p, i)->n;

				if (nneigh < 10)
					type = 6;
				else if (nneigh < 12)
					type = 4;
				else if (nneigh > 14)
					type = 7;
				else {
					co_x0 = 0;
					co_x1 = 0;
					co_x2 = 0;

					for (j = 0; j < nneigh; j++) {
						real *tmp_ptr = &(NEIGH(p,i)->dist[3*j]);
						v1.x = *tmp_ptr; tmp_ptr++;
						v1.y = *tmp_ptr; tmp_ptr++;
						v1.z = *tmp_ptr;
						v1_length = sqrt(SPROD(v1,v1));
						for (jj = 0; jj < j; jj++) {
							real *tmp_ptr  = &NEIGH(p,i)->dist[3*jj];
							v2.x = *tmp_ptr; tmp_ptr++;
							v2.y = *tmp_ptr; tmp_ptr++;
							v2.z = *tmp_ptr;
							v2_length = sqrt(SPROD(v2,v2));
							angle = SPROD(v1,v2);
							angle /= (v1_length * v2_length);
							if (angle < -.945)
								co_x0++;
							else if (angle < -.915)
								co_x1++;
							else if (angle < -.775)
								co_x2++;
						}
					}

					if (co_x0 == 7 && nneigh == 14)
						type = 0;
					else if (co_x0 == 6 && nneigh == 12)
						type = 1;
					else if (co_x0 == 3 && co_x1 <= 1 && co_x2 > 2 && nneigh == 12)
						type = 2;
					else if (nneigh > 12)
						type = 5;
					else if (nneigh == 12)
						type = 3;
					else
						type = 4;
				}
				ADATYPE(p, i) = type;
			}
		}
	} else if (ada_crystal_structure == ADA_BCC_CONFIG) {
		/*
		 * type=0: bcc
		 * type=1: fcc
		 * type=2: hcp
		 * type=3: 14 neighbors, no configuration
		 * type=4: 12-13 neighbors
		 * type=5: 15 neighbors
		 * type=6: less than 11 neighbors
		 * type=7: unknown
		 */
		for (k = 0; k < ncells; k++) {
			cell *p = CELLPTR(k);
			for (i = 0; i < p->n; i++) {
				nneigh = NEIGH(p, i)->n;

				if (nneigh < 11) type = 6;
				else if (nneigh == 13) type = 4;
				else if (nneigh == 15) type = 5;
				else if (nneigh > 15) type = 7;
				else {
					co_x0 = 0;
					co_x1 = 0;
					co_x2 = 0;
					co_x3 = 0;
					co_x4 = 0;

					for (j = 0; j < nneigh; j++) {
						real *tmp_ptr = &(NEIGH(p,i)->dist[3*j]);
						v1.x = *tmp_ptr; tmp_ptr++;
						v1.y = *tmp_ptr; tmp_ptr++;
						v1.z = *tmp_ptr;
						v1_length = sqrt(SPROD(v1,v1));
						for (jj = 0; jj < j; jj++) {
							real *tmp_ptr  = &NEIGH(p,i)->dist[3*jj];
							v2.x = *tmp_ptr; tmp_ptr++;
							v2.y = *tmp_ptr; tmp_ptr++;
							v2.z = *tmp_ptr;
							v2_length = sqrt(SPROD(v2,v2));
							angle = SPROD(v1,v2);
							angle /= (v1_length * v2_length);

							if (angle < -.945)
								co_x0++;
							else if (angle < -.915)
								co_x1++;
							else if (angle < -.755)
								co_x2++;
							else if (angle < -.705)
								co_x3++;
							else if ((angle > -.195) & (angle < .195))
								co_x4++;
						}
					}

					if (co_x0 > 5 && co_x0+co_x1==7 && nneigh==14)
						type = 0;
					else if (co_x0 == 6 && nneigh == 12)
						type = 1;
					else if (co_x0 == 3 && nneigh == 12)
						type = 2;
					else if (nneigh == 12)
						type = 4;
					else type = 3;
				}
				ADATYPE(p, i) = type;
			}
		}
	} else if (ada_crystal_structure == ADA_ACKLAND_CONFIG) {
		/*
		 * type=0: bcc
		 * type=1: fcc
		 * type=2: hcp
		 * type=3: unassigned
		 * type=4: unknown (usually surface)
		 */
		for (k = 0; k < ncells; k++) {
			cell *p = CELLPTR(k);
			for (i = 0; i < p->n; i++) {
				nneigh = NEIGH(p, i)->n;

				co_x0 = 0;
				co_x1 = 0;
				co_x2 = 0;
				co_x3 = 0;
				co_x4 = 0;
				co_x5 = 0;
				co_x6 = 0;
				co_x7 = 0;

				for (j = 0; j < nneigh; j++) {
					real *tmp_ptr = &(NEIGH(p,i)->dist[3*j]);
					v1.x = *tmp_ptr; tmp_ptr++;
					v1.y = *tmp_ptr; tmp_ptr++;
					v1.z = *tmp_ptr;
					v1_length = sqrt(SPROD(v1,v1));
					for (jj = 0; jj < j; jj++) {
						real *tmp_ptr  = &NEIGH(p,i)->dist[3*jj];
						v2.x = *tmp_ptr; tmp_ptr++;
						v2.y = *tmp_ptr; tmp_ptr++;
						v2.z = *tmp_ptr;
						v2_length = sqrt(SPROD(v2,v2));
						angle = SPROD(v1,v2);
						angle /= (v1_length * v2_length);

						if (angle < -.945)
							co_x0++;
						else if (angle < -.915)
							co_x1++;
						else if (angle < -.755)
							co_x2++;
						else if (angle < -.705)
							co_x3++;
						else if ((angle > -.195) && (angle < 195))
							co_x4++;
						else if (angle < -.245)
							co_x5++;
						else if (angle < -.795)
							co_x6++;
						else co_x7++;
					}
					delta_bcc = (0.35 * co_x4)
							/ (co_x5 + co_x6 + co_x7 - co_x4);
					delta_cp = 0.61 * ABS(1.-co_x6/24.);
					delta_fcc = 0.61 * (ABS(co_x0+co_x1-6) + co_x2) / 6.;
					delta_hcp = (ABS(co_x0-3) + ABS(co_x0+co_x1+co_x2+co_x3-9))
							/ 12.;

					if (nneigh < 11 || co_x7 > 0)
						type = 4;
					else if (co_x0 == 7)
						type = 0;
					else if (co_x0 == 6)
						type = 1;
					else if (co_x0 == 3)
						type = 2;
					else if (delta_bcc < 0.1 && delta_fcc < 0.1
							&& delta_hcp < 0.1 && delta_cp < 0.1)
						type = 3;
					else if (delta_bcc < delta_cp && 10 < nneigh && nneigh < 13)
						type = 0;
					else if (nneigh > 12)
						type = 3;
					else if (delta_hcp < delta_fcc)
						type = 2;
					else
						type = 1;
				}
				ADATYPE(p, i) = type;
			}
		}
	} else {
		error("Crystal structure not supported by ADA");
	}
}

void buildHopsToDefect(){
	int k,i,j, nneigh ,ii, l, m;

	/* Init all cells including buffer cells */
	for (k = 0; k < cell_dim.x; k++) {
		for (l = 0; l < cell_dim.y; l++) {
			for (m = 0; m < cell_dim.z; m++) {
				minicell *p = PTR_3D_V(cell_array, k, l, m, cell_dim);
				for	(i = 0; i < p->n; i++)
					HOPSTODEFECT(p,i) = 127;
			}
		}
	}

	/*mark all defects*/
	for (k = 0; k < ncells; k++) {
		cell *p = CELLPTR(k);
		for (i = 0; i < p->n; i++) {
			if (ADATYPE(p,i)!=ada_default_type){
				HOPSTODEFECT(p,i) = 0;
			}
		}
	}

	/* Mark all neighbors to a defect, which are not a defect itself */
	for (k = 0; k < ncells; k++) {
		cell *p = CELLPTR(k);
		for (i = 0; i < p->n; i++) {
			if (HOPSTODEFECT(p,i) == 0) {
				nneigh = NEIGH(p, i)->n;
				for (j = 0; j < nneigh; j++) {
					cell *q = NEIGH(p, i)->cl[j];
					ii = NEIGH(p, i)->num[j];
					if (HOPSTODEFECT(q,ii) == 127)
						HOPSTODEFECT(q,ii) = 1;
				}
			}
		}
	}

	send_fromBufferToCells(copy_hopsToDefect,pack_hopsToDefect,unpack_hopsToDefect);

	/* Mark second layer */
	for (k = 0; k < ncells; k++) {
		cell *p = CELLPTR(k);
		for (i = 0; i < p->n; i++) {
			if (HOPSTODEFECT(p,i) == 1) {
				nneigh = NEIGH(p, i)->n;
				for (j = 0; j < nneigh; j++) {
					cell *q = NEIGH(p, i)->cl[j];
					ii = NEIGH(p, i)->num[j];
					if (HOPSTODEFECT(q,ii) == 127)
						HOPSTODEFECT(q,ii) = 2;
				}
			}
		}
	}

	send_fromBufferToCells(copy_hopsToDefect,pack_hopsToDefect,unpack_hopsToDefect);

	/* Mark third layer */
	for (k = 0; k < ncells; k++) {
		cell *p = CELLPTR(k);
		for (i = 0; i < p->n; i++) {
			if (HOPSTODEFECT(p,i) == 2) {
				nneigh = NEIGH(p, i)->n;
				for (j = 0; j < nneigh; j++) {
					cell *q = NEIGH(p, i)->cl[j];
					ii = NEIGH(p, i)->num[j];
					if (HOPSTODEFECT(q,ii) == 127)
						HOPSTODEFECT(q,ii) = 3;
				}
			}
		}
	}

	send_fromBufferToCells(copy_hopsToDefect,pack_hopsToDefect,unpack_hopsToDefect);
}

/******************************************************************************
*
*  pack hopsToDefect from buffer cell into MPI buffer
*
******************************************************************************/
void pack_hopsToDefect( msgbuf *b, int k, int l, int m)
{
  int i, j = b->n;
  minicell *from;

  from = PTR_3D_V(cell_array, k, l, m, cell_dim);

  for (i=0; i<from->n; ++i) {
    b->data[ j++ ] = HOPSTODEFECT(from,i);
  }

  b->n = j;
  if (b->n_max < b->n)
    error("Buffer overflow in pack_hopsToDefect - increase msgbuf_size");
}

/******************************************************************************
*
*  unpack hopsToDefect from MPI buffer, and add them to those of the original cell
*
******************************************************************************/

void unpack_hopsToDefect( msgbuf *b, int k, int l, int m )
 {
	 int i, j = b->n;
	  minicell *to;
	  char hops;

	  to = PTR_3D_V(cell_array, k, l, m, cell_dim);

	  for (i=0; i<to->n; ++i) {
		hops = (char)(b->data[ j++ ]);
		if (hops<HOPSTODEFECT(to,i))
			HOPSTODEFECT(to,i) = hops;
	  }

	  b->n = j;
	  if (b->n_max < b->n)
	    error("Buffer overflow in unpack_hopsToDefect - increase msgbuf_size");
}

void copy_AdaAndNumber(int k, int l, int m, int r, int s, int t) {
	int i;
	minicell *from, *to;

	from = PTR_3D_V(cell_array, k, l, m, cell_dim);
	to = PTR_3D_V(cell_array, r, s, t, cell_dim);

	for (i = 0; i < to->n; ++i) {
		ADATYPE(to, i) = ADATYPE(from, i);
		NUMMER(to, i) = NUMMER(from, i);
	}
}

void pack_AdaAndNumber(msgbuf *b, int k, int l, int m) {
	int i, j = b->n;
	minicell *from;

	from = PTR_3D_V(cell_array, k, l, m, cell_dim);

	for (i = 0; i < from->n; ++i) {
		b->data[j++] = ADATYPE(from, i);
		b->data[j++] = NUMMER(from, i);
	}

	b->n = j;
	if (b->n_max < b->n)
		error("Buffer overflow in pack_ADA - increase msgbuf_size");
}

void unpack_AdaAndNumber(msgbuf *b, int k, int l, int m) {
	int i, j = b->n;
	minicell *to;

	to = PTR_3D_V(cell_array, k, l, m, cell_dim);

	for (i = 0; i < to->n; ++i) {
		ADATYPE(to, i) = b->data[j++];
		NUMMER(to, i) = b->data[j++];
	}

	b->n = j;
	if (b->n_max < b->n)
		error("Buffer overflow in unpack_ADA - increase msgbuf_size");
}

/******************************************************************************
 *
 *  copy contents of one cell to another (buffer) cell
 *
 ******************************************************************************/
void copy_hopsToDefect(int k, int l, int m, int r, int s, int t) {
	int i;
	minicell *from, *to;

	from = PTR_3D_V(cell_array, k, l, m, cell_dim);
	to = PTR_3D_V(cell_array, r, s, t, cell_dim);

	for (i = 0; i < to->n; ++i) {
		if (HOPSTODEFECT(from, i) < HOPSTODEFECT(to, i))
			HOPSTODEFECT(to, i) = HOPSTODEFECT(from, i);
	}
}

void send_fromBufferToCells(void(*add_func)(int, int, int, int, int, int),
		void(*pack_func)(msgbuf*, int, int, int), void(*unpack_func)(msgbuf*, int, int, int)) {
	int i, j;

#ifdef MPI
	MPI_Status stateast[2], statwest[2];
	MPI_Status statnorth[2], statsouth[2];
	MPI_Status statup[2], statdown[2];

	MPI_Request reqeast[2], reqwest[2];
	MPI_Request reqnorth[2], reqsouth[2];
	MPI_Request requp[2], reqdown[2];

	empty_mpi_buffers();
#endif
	/* send hopsToDefect east/west */
	if (cpu_dim.x == 1) {
		for (i = 0; i < cell_dim.y; ++i)
			for (j = 0; j < cell_dim.z; ++j) {
				(*add_func)( 0, i, j, cell_dim.x - 2, i, j);
				(*add_func)( cell_dim.x - 1, i, j, 1, i, j);
			}
	} else {
#ifdef MPI
		/* copy east hopsToDefect into send buffer, send east */
		for (i = 0; i < cell_dim.y; ++i)
			for (j = 0; j < cell_dim.z; ++j)
				(*pack_func)(&send_buf_east, 0, i, j);
		irecv_buf(&recv_buf_west, nbwest, &reqwest[1]);
		isend_buf(&send_buf_east, nbeast, &reqwest[0]);

		/* copy west hopsToDefect into send buffer, send west */
		for (i = 0; i < cell_dim.y; ++i)
			for (j = 0; j < cell_dim.z; ++j)
				(*pack_func)(&send_buf_west, cell_dim.x - 1, i, j);
		irecv_buf(&recv_buf_east, nbeast, &reqeast[1]);
		isend_buf(&send_buf_west, nbwest, &reqeast[0]);

		/* wait for hopsToDefect from west, add them to original cells */
		MPI_Waitall(2, reqwest, statwest);
		recv_buf_west.n = 0;
		for (i = 0; i < cell_dim.y; ++i)
			for (j = 0; j < cell_dim.z; ++j)
				(*unpack_func)(&recv_buf_west, cell_dim.x - 2, i, j);

		/* wait for hopsToDefect from east, add them to original cells */
		MPI_Waitall(2, reqeast, stateast);
		recv_buf_east.n = 0;
		for (i = 0; i < cell_dim.y; ++i)
			for (j = 0; j < cell_dim.z; ++j)
				(*unpack_func)(&recv_buf_east, 1, i, j);
#endif
	}

	/* send hopsToDefect north/south */
	if (cpu_dim.y == 1) {
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 0; j < cell_dim.z; ++j) {
				(*add_func)( i, 0, j, i, cell_dim.y - 2, j);
				(*add_func)( i, cell_dim.y - 1, j, i, 1, j);
			}
	} else {
#ifdef MPI
		/* copy north hopsToDefect into send buffer, send north */
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 0; j < cell_dim.z; ++j)
				(*pack_func)(&send_buf_north, i, 0, j);
		irecv_buf(&recv_buf_south, nbsouth, &reqsouth[1]);
		isend_buf(&send_buf_north, nbnorth, &reqsouth[0]);

		/* copy south hopsToDefect into send buffer, send south */
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 0; j < cell_dim.z; ++j)
				(*pack_func)(&send_buf_south, i, cell_dim.y - 1, j);
		irecv_buf(&recv_buf_north, nbnorth, &reqnorth[1]);
		isend_buf(&send_buf_south, nbsouth, &reqnorth[0]);

		/* wait for hopsToDefect from south, add them to original cells */
		MPI_Waitall(2, reqsouth, statsouth);
		recv_buf_south.n = 0;
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 0; j < cell_dim.z; ++j)
				(*unpack_func)(&recv_buf_south, i, cell_dim.y - 2, j);

		/* wait for hopsToDefect from north, add them to original cells */
		MPI_Waitall(2, reqnorth, statnorth);
		recv_buf_north.n = 0;
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 0; j < cell_dim.z; ++j)
				(*unpack_func)(&recv_buf_north, i, 1, j);
#endif
	}

	/* send hopsToDefect up/down */
	if (cpu_dim.z == 1) {
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 1; j < cell_dim.y - 1; ++j) {
				(*add_func)(i, j, 0, i, j, cell_dim.z - 2);
				(*add_func)(i, j, cell_dim.z - 1, i, j, 1);
			}
	} else {
#ifdef MPI
		/* copy up hopsToDefect into send buffer, send up */
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 1; j < cell_dim.y - 1; ++j)
				(*pack_func)(&send_buf_up, i, j, 0);
		irecv_buf(&recv_buf_down, nbdown, &reqdown[1]);
		isend_buf(&send_buf_up, nbup, &reqdown[0]);

		/* copy down hopsToDefect into send buffer, send down */
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 1; j < cell_dim.y - 1; ++j)
				(*pack_func)(&send_buf_down, i, j, cell_dim.z - 1);
		irecv_buf(&recv_buf_up, nbup, &requp[1]);
		isend_buf(&send_buf_down, nbdown, &requp[0]);

		/* wait for hopsToDefect from down, add them to original cells */
		MPI_Waitall(2, reqdown, statdown);
		recv_buf_down.n = 0;
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 1; j < cell_dim.y - 1; ++j)
				(*unpack_func)(&recv_buf_down, i, j, cell_dim.z - 2);

		/* wait for hopsToDefect from up, add them to original cells */
		MPI_Waitall(2, requp, statup);
		recv_buf_up.n = 0;
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 1; j < cell_dim.y - 1; ++j)
				(*unpack_func)(&recv_buf_up, i, j, 1);
#endif
	}
}

/******************************************************************************
 *copy data of one cell into buffer cells
 ******************************************************************************/
void send_fromCellsToBuffer(void(*copy_func)(int, int, int, int, int, int),
		void(*pack_func)(msgbuf*, int, int, int), void(*unpack_func)(msgbuf*, int, int, int)) {
	int i, j;

#ifdef MPI
	MPI_Status stateast[2], statwest[2];
	MPI_Status statnorth[2], statsouth[2];
	MPI_Status statup[2], statdown[2];

	MPI_Request reqeast[2], reqwest[2];
	MPI_Request reqnorth[2], reqsouth[2];
	MPI_Request requp[2], reqdown[2];

	empty_mpi_buffers();
#endif
	/* exchange up/down */
	if (cpu_dim.z == 1) {
		/* simply copy up/down atoms to buffer cells*/
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 1; j < cell_dim.y - 1; ++j) {
				(*copy_func)(i, j, 1, i, j, cell_dim.z - 1);
				(*copy_func)(i, j, cell_dim.z - 2, i, j, 0);
			}
	}

	else {
#ifdef MPI
		/* copy up atoms into send buffer, send up */
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 1; j < cell_dim.y - 1; ++j)
				(*pack_func)(&send_buf_up, i, j, 1);
		irecv_buf(&recv_buf_down, nbdown, &reqdown[1]);
		isend_buf(&send_buf_up, nbup, &reqdown[0]);

		/* copy down atoms into send buffer, send down */
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 1; j < cell_dim.y - 1; ++j)
				(*pack_func)(&send_buf_down, i, j, cell_dim.z - 2);
		irecv_buf(&recv_buf_up, nbup, &requp[1]);
		isend_buf(&send_buf_down, nbdown, &requp[0]);

		/* wait for atoms from down, move them to buffer cells */
		MPI_Waitall(2, reqdown, statdown);
		recv_buf_down.n = 0;
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 1; j < cell_dim.y - 1; ++j)
				(*unpack_func)(&recv_buf_down, i, j, cell_dim.z - 1);

		/* wait for atoms from up, move them to buffer cells*/
		MPI_Waitall(2, requp, statup);
		recv_buf_up.n = 0;
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 1; j < cell_dim.y - 1; ++j)
				(*unpack_func)(&recv_buf_up, i, j, 0);
#endif
	}

	/* exchange north/south */
	if (cpu_dim.y == 1) {
		/* simply copy north/south atoms to buffer cells */
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 0; j < cell_dim.z; ++j) {
				(*copy_func)(i, 1, j, i, cell_dim.y - 1, j);
				(*copy_func)(i, cell_dim.y - 2, j, i, 0, j);
			}
	}

	else {
#ifdef MPI
		/* copy north atoms into send buffer, send north */
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 0; j < cell_dim.z; ++j)
				(*pack_func)(&send_buf_north, i, 1, j);
		irecv_buf(&recv_buf_south, nbsouth, &reqsouth[1]);
		isend_buf(&send_buf_north, nbnorth, &reqsouth[0]);

		/* copy south atoms into send buffer, send south*/
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 0; j < cell_dim.z; ++j)
				(*pack_func)(&send_buf_south, i, cell_dim.y - 2, j);
		irecv_buf(&recv_buf_north, nbnorth, &reqnorth[1]);
		isend_buf(&send_buf_south, nbsouth, &reqnorth[0]);

		/* wait for atoms from south, move them to buffer cells */
		MPI_Waitall(2, reqsouth, statsouth);
		recv_buf_south.n = 0;
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 0; j < cell_dim.z; ++j)
				(*unpack_func)(&recv_buf_south, i, cell_dim.y - 1, j);

		/* wait for atoms from north, move them to buffer cells*/
		MPI_Waitall(2, reqnorth, statnorth);
		recv_buf_north.n = 0;
		for (i = 1; i < cell_dim.x - 1; ++i)
			for (j = 0; j < cell_dim.z; ++j)
				(*unpack_func)(&recv_buf_north, i, 0, j);
#endif
	}

	/* exchange east/west*/
	if (cpu_dim.x == 1) {
		/* simply copy east/west atoms to buffer cells*/
		for (i = 0; i < cell_dim.y; ++i)
			for (j = 0; j < cell_dim.z; ++j) {
				(*copy_func)(1, i, j, cell_dim.x - 1, i, j);

				(*copy_func)(cell_dim.x - 2, i, j, 0, i, j);

			}
	}

	else {
#ifdef MPI
		/* copy east atoms into send buffer, send east*/
		for (i = 0; i < cell_dim.y; ++i)
			for (j = 0; j < cell_dim.z; ++j)
				(*pack_func)(&send_buf_east, 1, i, j);
		irecv_buf(&recv_buf_west, nbwest, &reqwest[1]);
		isend_buf(&send_buf_east, nbeast, &reqwest[0]);

		/* copy west atoms into send buffer, send west*/
		for (i = 0; i < cell_dim.y; ++i)
			for (j = 0; j < cell_dim.z; ++j)
				(*pack_func)(&send_buf_west, cell_dim.x - 2, i, j);
		irecv_buf(&recv_buf_east, nbeast, &reqeast[1]);
		isend_buf(&send_buf_west, nbwest, &reqeast[0]);

		/* wait for atoms from west, move them to buffer cells*/
		MPI_Waitall(2, reqwest, statwest);
		recv_buf_west.n = 0;
		for (i = 0; i < cell_dim.y; ++i)
			for (j = 0; j < cell_dim.z; ++j)
				(*unpack_func)(&recv_buf_west, cell_dim.x - 1, i, j);

		/* wait for atoms from east, move them to buffer cells*/
		MPI_Waitall(2, reqeast, stateast);
		recv_buf_east.n = 0;
		for (i = 0; i < cell_dim.y; ++i)
			for (j = 0; j < cell_dim.z; ++j)
				(*unpack_func)(&recv_buf_east, 0, i, j);
#endif
	}
}
