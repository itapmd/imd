
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
* imd_misc.c -- Some Misc. Routines for the imd package
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/******************************************************************************
*
*  Usage -- educate users
*
******************************************************************************/

void usage(void)
{ 
  if (myid==0) {
    fprintf(stderr,"%s [-r<nnn>] [-p paramter-file]\n",progname); 
    fflush(stderr);
  }
#ifdef MPI
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  exit(1); 
}

/******************************************************************************
*
*  print a warning
*
******************************************************************************/

void warning(char *msg)
{ 
  if (myid==0) {
    fprintf(stderr,"WARNING: %s\n",msg);
    fflush(stderr);
  }
}

/******************************************************************************
*
* error -- Complain and abort
*
******************************************************************************/

void error(char *msg)
{
#ifdef MPI
  fprintf(stderr,"Error on CPU %d: %s\n",myid,msg);
#else
  fprintf(stderr,"Error: %s\n",msg);
#endif
  /* try to flush and close whatever we can */
  fflush(stderr);
  fflush(stdout);
  if ((myid==0) && ( eng_file!=NULL)) fclose( eng_file);
  if ((myid==0) && (msqd_file!=NULL)) fclose(msqd_file);
#ifdef MPI
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  exit(2);
}

#ifdef AVPOS

/******************************************************************************
*
* add_positions adds coordinates and potential energy for computation of average
* 
* position and potential energy
*
******************************************************************************/

void add_positions(void)
{
  int k;
  for (k=0; k<ncells; k++) {
    int i;
    cell* p;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; i++) {
        p->avpos X(i) += p->ort X(i) + p->sheet X(i);
        p->avpos Y(i) += p->ort Y(i) + p->sheet Y(i);
#ifndef TWOD
        p->avpos Z(i) += p->ort Z(i) + p->sheet Z(i);
#endif
        p->av_epot[i] += p->pot_eng[i];
    }
  }

#ifdef NPT
  av_box_x.x += box_x.x;
  av_box_x.y += box_x.y;
  av_box_y.x += box_y.x;
  av_box_y.y += box_y.y;
#ifndef TWOD
  av_box_x.z += box_x.z;  
  av_box_y.z += box_y.z;
  av_box_z.x += box_z.x;
  av_box_z.y += box_z.y;  
  av_box_z.z += box_z.z;  
#endif
#endif

}
  
#endif

/******************************************************************************
*
*  endian returns 1 if system is big endian, 0 if little endian
*
******************************************************************************/

int endian(void)
{
  unsigned short int word = 0x0001;
  unsigned char  *byte    = (unsigned char *) &word;
  return (byte[0] ? 0 : 1);
}

#ifdef PAIR_PRE

/******************************************************************************
*
*  init_pot_par -- initialize parameters for predefinded pair potentials
*
******************************************************************************/

void init_pot_par(void) {

  int  i, j, n = 0;
  real tmp = 0.0;

  for (i=0; i<ntypes; i++) 
    for (j=i; j<ntypes; j++) {

      r_cut[i][j]  = r_cut[j][i]  = r_cut_lin[n];
      r2_cut[i][j] = r2_cut[j][i] = SQR(r_cut[i][j]);

#ifdef LJ
      lj_epsilon[i][j] = lj_epsilon[j][i] = lj_epsilon_lin[n]; 
      lj_sigma[i][j]   = lj_sigma[j][i]   = lj_sigma_lin[n]; 
#endif
#ifdef MORSE
      morse_epsilon[i][j] = morse_epsilon[j][i] = morse_epsilon_lin[n]; 
      morse_sigma[i][j]   = morse_sigma[j][i]   = morse_sigma_lin[n]; 
      morse_alpha[i][j]   = morse_alpha[j][i]   = morse_alpha_lin[n]; 
#endif
#ifdef BUCK
      buck_a[i][j]      = buck_a[j][i]      = buck_a_lin[n];
      buck_c[i][j]      = buck_c[j][i]      = buck_c_lin[n];
      buck_sigma[i][j]  = buck_sigma[j][i]  = buck_sigma_lin[n];
#endif      
      n++;
    }

  for (i=0; i<ntypes; ++i)
    for (j=0; j<ntypes; ++j)
      tmp = MAX( tmp, r2_cut[i][j] );
  
  cellsz = MAX(cellsz,tmp);

  /* Shift of potentials */
  for (i=0; i<ntypes; i++) 
      for (j=0; j<ntypes; j++) {
#ifdef LJ
	  tmp = pow( lj_sigma[i][j] / r_cut[i][j], 6 );
	  lj_shift[i][j] = lj_epsilon[i][j] * ( tmp * tmp - 2 * tmp );
	  printf("Lennard-Jones potential %1d %1d shifted by %f\n", 
		 i, j, -lj_shift[i][j]);
#endif
#ifdef MORSE
	  tmp = 1.0 - exp( - morse_alpha[i][j] * 
			   ( r_cut[i][j] - morse_sigma[i][j] ) );
	  morse_shift[i][j] = morse_epsilon[i][j] *  ( tmp * tmp - 1.0 );
	  printf("Morse potential %1d %1d shifted by %f\n", 
		 i, j, -morse_shift[i][j]);
#endif
#ifdef BUCK
	  buck_shift[i][j] = buck_a[i][j] 
	      * exp( - r_cut[i][j] / buck_sigma[i][j] )
	      - buck_c[i][j] * pow( buck_sigma[i][j] / r_cut[i][j] , 6 );
	  printf("Buckingham potential %1d %1d shifted by %f\n", 
		 i, j, -buck_shift[i][j]); 
#endif
    }
}

#endif




