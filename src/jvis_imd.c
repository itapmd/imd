#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define MAIN
#include "imd.h"

void   output(void);
double calc_diff(void);

double box[2][3];
double Ekin_av, Epot_av, press_av;
int    sample=100, run=1; 

int main(int argc, char **argv) {

   int i, j;
   char input[256];
   cell *incell;

   // read parameter file
   strcpy(paramfilename, argv[1]);
   read_parameters();

   // setup potential
   read_pot_table(&pair_pot, potfilename, ntypes*ntypes, 1);

   // read box size
   for (i = 0; i < 2; ++i) {
      for (j = 0; j < 3; ++j) {
         scanf("%lf", &(box[i][j]));
      }
   }
#ifdef TWOD
   box_x.x = box[1][0] - box[0][0]; box_x.y = 0.0;
   box_y.x = 0.0; box_y.y = box[1][1] - box[0][1];
#else
   box_x.x = box[1][0] - box[0][0]; box_x.y = 0.0; box_x.z = 0.0;
   box_y.x = 0.0; box_y.y = box[1][1] - box[0][1]; box_y.z = 0.0;
   box_z.x = 0.0; box_z.y = 0.0; box_y.y = box[1][2] - box[0][2];
#endif
   make_box();

   // allocate input cell for one atom
   incell = (cell *) malloc(sizeof(cell));
   if (0==incell) error("Cannot allocate input cell.");
   incell->n_max=0;
   alloc_cell(incell, 1);

   // scan atoms
   scanf("%d", &natoms);
   nactive = DIM * natoms;
   for (i = 0; i < natoms; ++i) {
      double tmp[3];
      ivektor cellc;
      cell* to;
      scanf("%s %lf %lf %lf", input, tmp, tmp+1, tmp+2);
      for (j = 0; j < DIM; ++j) {
	 tmp[j]           -= box[0][j];
	 incell->ort[j]    = tmp[j];
	 incell->impuls[j] = 0.0;
	 incell->kraft[j]  = 0.0;
      }
      incell->nummer[0] = i;
      incell->sorte[0]  = 0;
      incell->vsorte[0] = 0;
      incell->masse[0]  = 1.0;
#ifdef TWOD
      cellc = cell_coord(tmp[0],tmp[1]);
#else
      cellc = cell_coord(tmp[0],tmp[1],tmp[2]);
#endif
#ifdef BUFCELLS
      cellc = local_cell_coord(cellc);
#endif
      to = PTR_VV(cell_array,cellc,cell_dim);
      incell->n=1;
      INSERT_ATOM(to, incell, 0);
   }

   maxwell(temperature);
   init_refpos();
   output();

   do {
      do {
         double val, val2;
         scanf("%s", input);

         if (strcmp(input, "quit") == 0) {
            return 0;
         } 
         else if (strcmp(input, "next") == 0) {
            break;
         } 
         else if (strcmp(input, "ekin") == 0) {
            scanf("%lf", &val);
            temperature = (2.0 / DIM) * val;
            eta  = 0.0;
            xi.x = 0.0;
         } 
         else if (strcmp(input, "vel") == 0) {
            scanf("%lf", &val);
            temperature = (2.0 / DIM) * val;
            maxwell(temperature);
         } 
/*
         else if (strcmp(input, "a") == 0) {
            scanf("%lf", &val);
            A_par = val;
            pot_init();
         } 
         else if (strcmp(input, "b") == 0) {
            scanf("%lf", &val);
            B_par = val;
            pot_init();
         } 
         else if (strcmp(input, "c") == 0) {
            scanf("%lf", &val);
            C_par = val;
            pot_init();
         }
         else if (strcmp(input, "size") == 0) {
            scanf("%lf", &val);
            scanf("%lf", &val2);
#ifdef TWOD
            resize(val, val2);
#else
            double val3;
            scanf("%lf", &val3);
            resize(val, val2, val3);
#endif
         }
*/ 
         else if (strcmp(input, "NPT") == 0) {
            ensemble   = ENS_NPT_ISO;
            move_atoms = move_atoms_npt_iso;
            eta  = 0.0;
            xi.x = 0.0;
         } 
         else if (strcmp(input, "NVT") == 0) {
            ensemble   = ENS_NVT;
            move_atoms = move_atoms_nvt;
            eta = 0.0;
         } 
         else if (strcmp(input, "NVE") == 0) {
	    ensemble   = ENS_NVE;
            move_atoms = move_atoms_nve;
         }
/* 
         else if (strcmp(input, "perbound") == 0) {
            scanf("%s", input);
            if (strcmp(input, "true") == 0) perbound = 1;
            if (strcmp(input, "false") == 0) perbound = 0;
         } 
         else if (strcmp(input, "run") == 0) {
            scanf("%s", input);
            if (strcmp(input, "true") == 0) run = 1;
            if (strcmp(input, "false") == 0) run = 0;
         } 
*/
         else if (strcmp(input, "sample") == 0) {
            scanf("%d", &sample);
         } 
         else if (strcmp(input, "press") == 0) {
            scanf("%lf", &val);
            pressure_ext.x = val;
         }
/*
         else if (strcmp(input, "alpha") == 0) {
            scanf("%lf", &val);
            alpha = val;
         } 
         else if (strcmp(input, "min") == 0) {
            minimize();
         } 
         else if (strcmp(input, "beta") == 0) {
            scanf("%lf", &val);
            beta = val;
         } 
         else if (strcmp(input, "start") == 0) {
            start_dsf();
            doDSF = 1;
         }
*/
      } while (1);

      Ekin_av  = 0.0;
      Epot_av  = 0.0;
      press_av = 0.0;
      init_refpos();
//      if (run) { 
         for (i = 0; i < sample; i++) {
            calc_forces(i);
            move_atoms();
            fix_cells();  
            // if (doDSF) loop_dsf();
            Ekin_av  += tot_kin_energy / natoms;
            Epot_av  += tot_pot_energy / natoms;
            press_av += (2.0 * tot_kin_energy + virial) / (DIM * volume); 
	 }
//      } else {
//         fill_cells();
//         calc_en();
//      }
      Ekin_av  /= sample;
      Epot_av  /= sample;
      press_av /= sample;
      output();
   } while (1);
}

double calc_diff() {
   double diff = 0.0;
   int i, k;
   for (k=0; k<NCELLS; ++k) {
      cell *p = CELLPTR(k);
      for (i = 0; i < p->n; ++i) {
         diff += SQR( ORT(p,i,X) - REF_POS(p,i,X) );
         diff += SQR( ORT(p,i,Y) - REF_POS(p,i,Y) );
#ifndef TWOD
         diff += SQR( ORT(p,i,Z) - REF_POS(p,i,Z) );
#endif
      }
   }
   return diff / (natoms * sample * timestep);
}

void output() {

   int i, k;

   // mean square displacement
   double d, diff = calc_diff();

   // energies, pressure, ...
   printf("start_output_here\n");
   printf("Ekin=%.4lf, Epot=%.4lf, press=%.4lf, diff=%.4lf, vol=%.0lf\n", 
          Ekin_av, Epot_av, press_av, diff, volume);
   printf("%.4lf\n", Ekin_av);
   printf("%.4lf\n", Epot_av);
   printf("%.4lf\n", press_av);
   printf("%.4lf\n", diff);

   // write box; keep it centered
   d = box[1][0] - box[0][0] - box_x.x;
   box[0][0] += d/2; box[1][0] -= d/2;
   d = box[1][1] - box[0][1] - box_y.y;
   box[0][1] += d/2; box[1][1] -= d/2;
#ifndef TWOD
   d = box[1][2] - box[0][2] - box_z.z;
   box[0][2] += d/2; box[1][2] -= d/2;
#endif
   for (i = 0; i < 2; ++i) {
#ifdef TWOD
      printf("%f\t %f\t %f\n", box[i][0], box[i][1], 0.0);
#else
      printf("%f\t %f\t %f\n", box[i][0], box[i][1], box[i][2]);
#endif
   }

   // write atom coordinates
   for (k=0; k<NCELLS; ++k) {
      cell *p = CELLPTR(k);
      for (i = 0; i < p->n; ++i) {
          printf("%f\t", ORT(p,i,X) + box[0][0] );
          printf("%f\t", ORT(p,i,Y) + box[0][1] );
#ifdef TWOD
          printf("0.0\n");
#else
          printf("%f\n", ORT(p,i,Z) + box[0][2] );
#endif
      }
   }
   fflush(stdout);
}

