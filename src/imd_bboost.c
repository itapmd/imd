/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2007 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
*
* imd_bboost.c -- Routines for the Bond Boost Method, see ....
*
******************************************************************************/

#include "imd.h"


void  init_bboost(void)
{
    int n, k, i, j;
    //M3 *V=NULL, *W=NULL;
    /* initializations usually done in main_loop */
#ifndef SBOOST
    printf(" begin bond boost MD !!! \n");fflush(stdout);      
#else
    printf(" begin strain boost MD !!! \n");fflush(stdout);  
#endif    

#ifdef EXTPOT
    init_extpot();
#endif  
  
#ifdef FBC
    init_fbc();
#endif

#if defined(CORRELATE) || defined(MSQD)
    init_correl(ncorr_rmax,ncorr_tmax);
#endif

#ifdef NMOLDYN
    if (nmoldyn_int > 0) init_nmoldyn();
#endif

#ifdef ATDIST
    if (atdist_int > 0) init_atdist();
#endif

#ifdef DIFFPAT
    if (diffpat_int > 0) init_diffpat();
#endif

#ifdef CG
    if (ensemble == ENS_CG) reset_cg();
#endif

#ifdef DEFORM
    deform_int = 0; 
#endif
#if defined(HOMDEF) && defined(RELAX)
    deform_int = 0; 
#endif
 
  /* allocate NOWPOS and REPOST on all CPUs */
    if ( natoms>0 ) {
        if((NOWPOS = (real *) calloc((natoms+1)*SDIM,sizeof(real)))==NULL)
         error("cannot allocate memory for NOWPOS\n"); 
        if((REPOST = (real *) calloc((natoms+1)*SDIM,sizeof(real)))==NULL)
         error("cannot allocate memory for REPOST\n");
        if((REACFJ = (real *) calloc((natoms+1)*SDIM,sizeof(real)))==NULL)
         error("cannot allocate memory for REPOST\n");        
       for (i=0; i<=SDIM*natoms; i++) {
           NOWPOS[i] = 0.0;
           REPOST[i] = 0.0;
           REACFJ[i] = 0.0;
           //     printf(" !!!!! myid: %d  myrank: %d !!!!! \n", myid,myrank);
       }
    }
    // timll = 0;
#ifndef SBOOST
    if ( natoms>0 ) {

        if((bb_vtypej  = (int *) calloc((natoms+1),sizeof(int)))==NULL)
         error("cannot allocate memory for bb_vtypej\n");                
       for (i=0; i<=natoms; i++) {
           bb_vtypej[i]  = 0;
           //           printf(" !!!!! VInv[%d] = %f !!!!! \n", i,VInv[i]);
       }       
    }
#else
    if ( natoms>0 ) {
        if((VInv  = (real *) calloc((natoms+1)*SDIM*SDIM,sizeof(real)))==NULL)
         error("cannot allocate memory for VInv\n");
        if((M3eta = (real *) calloc((natoms+1)*SDIM*SDIM,sizeof(real)))==NULL)
         error("cannot allocate memory for M3eta\n");
        if((M3J   = (real *) calloc((natoms+1)*SDIM*SDIM,sizeof(real)))==NULL)
         error("cannot allocate memory for M3J\n");
        if((Mises  = (real *) calloc((natoms+1),sizeof(real)))==NULL)
         error("cannot allocate memory for Mises\n");                
       for (i=0; i<=SDIM*SDIM*natoms; i++) {
           VInv[i] = 0.0;
           M3eta[i]= 0.0;
           M3J[i]  = 0.0;
           //           printf(" !!!!! VInv[%d] = %f !!!!! \n", i,VInv[i]);
       }
       for (i=0; i<=natoms; i++) {
           Mises[i]  = 0.0;
           //           printf(" !!!!! VInv[%d] = %f !!!!! \n", i,VInv[i]);
       }       
    }
    
#endif    
    
     printf(" bb_tot_bV = %f, bbp1_2 = %f, bb_rcut = %f: \n", bb_tot_bV,bbp1_2,bb_rcut);fflush(stdout);
     printf(" bb_epscrit = %f, delta_bb_epscrit = %f, bb_bound_div = %f: \n", bb_epscrit,delta_bb_epscrit,bb_bound_div);fflush(stdout);
     bb_box_x.x = box_x.x;
     bb_box_y.y = box_y.y;
     bb_box_z.z = box_z.z;
    /* store the positions and velocities before doing minimization */
    for (k=0; k<NCELLS; ++k) { /* loop over all cells */
        int  i,j, sort;
        cell *p;
        p = CELLPTR(k);
        for (i=0; i<p->n; ++i) { /* loop over all atoms in the cell */
            OLDPOS(p,i,X)    = ORT(p,i,X);
            OLDPOS(p,i,Y)    = ORT(p,i,Y);
            OLDPOS(p,i,Z)    = ORT(p,i,Z);
            OLDIMPULS(p,i,X) = IMPULS(p,i,X);
            OLDIMPULS(p,i,Y) = IMPULS(p,i,Y);
            OLDIMPULS(p,i,Z) = IMPULS(p,i,Z);
#ifndef SBOOST
            bb_vtypej[NUMMER(p,i)] = VSORTE(p,i);
            // printf(" bb_vtypej[%d] = %d \n",bb_vtypej[NUMMER(p,i)],NUMMER(p,i),VSORTE(p,i));fflush(stdout);
#else
#endif            
	  /*  BBNEIGH(p,i)->nbondsref1 = 0; */
        }
        //  if(myid==1)

    }
    //     error("check the code before here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    /* find reference positions by energy minimization */
    bb_minimize(0);
    max_epot   = tot_pot_energy / natoms;
    min_epot   = bb_epot;
    /* make relaxed position reference position */
    for (k=0; k<NCELLS; ++k) { /* loop over all cells */
        int  i,j, sort;
        cell *p;
        p = CELLPTR(k);
        for (i=0; i<p->n; ++i) { /* loop over all atoms in the cell */
            REFPOSONE(p,i,X) = REFPOSTWO(p,i,X);
            REFPOSONE(p,i,Y) = REFPOSTWO(p,i,Y);
            REFPOSONE(p,i,Z) = REFPOSTWO(p,i,Z);
        }
        
    }


    if ((bb_checkpt_int > 0) && (0 == 0 % bb_checkpt_int)) 
       write_config_select( 0/bb_checkpt_int, "bb_old",
                            write_atoms_bb, write_header_bb);

  /* update bb_neighbor table cutoff */
  if (NULL==bb_neightab_r2cut) {
     bb_neightab_r2cut = (real *) calloc( ntypes * ntypes, sizeof(real) );
    if (NULL==bb_neightab_r2cut) 
       error("cannot allocate memory for bb_neightab_r2cut");
  }
      
  for (i=0; i<ntypes*ntypes; i++) {
      bb_neightab_r2cut[i] = bb_rcut * bb_rcut;
      // printf("myid=%d myrank: %d bb_neightab_r2cut[%d] =%f bb_rcut = %f: \n",  myid,myrank,i,bb_neightab_r2cut[i],bb_rcut);fflush(stdout);
  }

  /* loop over all pairs of cells */
  for (n=0; n<nlists; ++n) {
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
#endif
             /* if not define AR, using k=npairs[n]; k<npairs2[n]*/
 /*       for (k=npairs[n]; k<npairs2[n]; ++k) { //??? */
          for (k=0; k<npairs[n]; ++k) {
              // printf("myid=%d : k = %d, steps = %d \n",myid,k,steps);fflush(stdout);           
            vektor pbc;
            pair *P;
            P = pairs[n]+k;
            pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
            pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
            pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
            do_bb_neightab(cell_array + P->np, cell_array + P->nq, pbc);
        }
    }

 
} /* end of init_bboost */




/******************************************************************************
*
*  do_bb_neightab - compute neighbor table
*
******************************************************************************/

void do_bb_neightab(cell *p, cell *q, vektor pbc)
{
  int i, j, k;
  int jstart, jend;
  int q_typ, p_typ, column;
  vektor d, tmp_d;
  real *qptr, radius2;

  /* For each atom in first cell */
  for (i=0; i<p->n; ++i) {
       tmp_d.x = REFPOSONE(p,i,X) - pbc.x;
       tmp_d.y = REFPOSONE(p,i,Y) - pbc.y;
#ifndef TWOD
       tmp_d.z = REFPOSONE(p,i,Z) - pbc.z;
#endif
       p_typ   = SORTE(p,i);

#ifdef TWOD
       jstart  = (((p==q) && (pbc.x==0) && (pbc.y==0))               ? i+1 : 0);
#else
       jstart  = (((p==q) && (pbc.x==0) && (pbc.y==0) && (pbc.z==0)) ? i+1 : 0);
#endif
       qptr    = &REFPOSONE(q,jstart,X);

    /* For each atom in neighbouring cell */
      for (j = jstart; j < q->n; ++j) {

           q_typ = SORTE(q,j);
  
      /* Calculate distance dxji[]=x_j[]-x_i[] and |dxji|^2 */
           //           d.x = *qptr - tmp_d.x; ++qptr;
           // d.y = *qptr - tmp_d.y; ++qptr;
           d.x = *qptr++ - tmp_d.x;
           d.y = *qptr++ - tmp_d.y;           
#ifndef TWOD
           //           d.z = *qptr - tmp_d.z; ++qptr;
           d.z = *qptr++ - tmp_d.z;
#endif  
 
           column  = p_typ * ntypes + q_typ;
           radius2 = SPROD(d,d);
           if (0==radius2) { char msgbuf[256];
               sprintf(msgbuf,
               "Distance is zero: nrs=%d %d\norte: %f %f %f, %f %f %f\n",
               NUMMER(p,i),NUMMER(q,j),
               REFPOSONE(p,i,X), REFPOSONE(p,i,Y), REFPOSONE(p,i,Z),
               REFPOSONE(q,j,X), REFPOSONE(q,j,Y), REFPOSONE(q,j,Z) );
               error(msgbuf);
           }
           
      /* make neighbor tables for boost methods */

           if (radius2 <= bb_neightab_r2cut[column]) {
               bb_neightab *bb_neigh;
               real  *tmp_ptr;
               bb_neigh = BBNEIGH(p,i);
               bb_neigh->numref1[bb_neigh->nbondsref1] = NUMMER(q,j);
               bb_neigh->distref1[bb_neigh->nbondsref1] = radius2;
               tmp_ptr  = &bb_neigh->vectref1[SDIM*bb_neigh->nbondsref1];
               *tmp_ptr = d.x; ++tmp_ptr; 
               *tmp_ptr = d.y; ++tmp_ptr; 
               *tmp_ptr = d.z;               
               //*tmp_ptr = radius2;
               bb_neigh->nbondsref1++;
               //if(NUMMER(p,i) == 769)
               //printf("nbondsref1 = %d\n",bb_neigh->nbondsref1);fflush(stdout);
	/* update bb_neightab of j for the actio = reactio */
#ifdef SBOOST
               bb_neigh = BBNEIGH(q,j);
               bb_neigh->numref1[bb_neigh->nbondsref1] = NUMMER(p,i);
               bb_neigh->distref1[bb_neigh->nbondsref1] = radius2;
               tmp_ptr  = &bb_neigh->vectref1[SDIM*bb_neigh->nbondsref1];
               *tmp_ptr = -d.x; ++tmp_ptr; 
               *tmp_ptr = -d.y; ++tmp_ptr; 
               *tmp_ptr = -d.z;                
               //*tmp_ptr = radius2;
               bb_neigh->nbondsref1++;
#endif
               //if(NUMMER(q,j) == 769)
               //printf("nbondsref1 = %d\n",bb_neigh->nbondsref1);fflush(stdout);
           }
      } /* for j */
  } /* for i */

}  /* the end of do_bb_neightab */






/* this is the minimization part of the bond boost method */
/* minimization is carried out according to which option you compiled
   and the parameters in the param file
   this is a copy of main_loop 
 */
int bb_minimize(int simulation)
{
  int  finished = 0;
  int  i, j, k, l;
  int  steps_diff = steps_max - steps_min;
  int  deform_int = 0, do_fbc_incr = 0;
  int  have_fbc_incr = 0;
  real dtemp, dshock_speed;
  vektor d_pressure;
  real tmpvec1[DIM], tmpvec2[DIM];
  char tmp_str[9];
  real fnorm2,ekin,epot,delta_epot;
  real bb_tot_pot_energy;

#ifdef GLOK
  if(glok_start <=bbminsteps_min)
      glok_start=bbminsteps_min; 
#endif
#ifdef ACG
  acg_alpha = acg_init_alpha;
#endif
#ifdef RELAX
  is_relaxed = 0;
#endif
    /* record bounadry imformation before doing relaxation, this is for some minimize methods that boundary can change */
  tmp_box_x.x = box_x.x;
  tmp_box_y.y = box_y.y;
  tmp_box_z.z = box_z.z;
    /* record atoms imformation before doing relaxation, we needn't record velocity, velocity will not be changed before and after cg */
  for (k=0; k<NCELLS; ++k) { /* loop over all cells */
       int  i,j, sort;
       cell *p;
       p = CELLPTR(k);
      for (i=0; i<p->n; ++i) { /* loop over all atoms in the cell */
           TMPPOS(p,i,X) = ORT(p,i,X);
           TMPPOS(p,i,Y) = ORT(p,i,Y);
           TMPPOS(p,i,Z) = ORT(p,i,Z);
  //      if(NUMMER(p,i) == 769)
  //      printf(" IMPULS_X = %f, IMPULS_Y = %f, IMPULS_Z = %f\n", IMPULS(p,i,X),IMPULS(p,i,Y),IMPULS(p,i,Z));fflush(stdout);
           TMPKRAFT(p,i,X) = KRAFT(p,i,X);
           TMPKRAFT(p,i,Y) = KRAFT(p,i,Y);
           TMPKRAFT(p,i,Z) = KRAFT(p,i,Z);
//        if(NUMMER(p,i) == 769)
//        printf(" KRAFT_X = %f, KRAFT_Y = %f, KRAFT_Z = %f\n", KRAFT(p,i,X),KRAFT(p,i,Y),KRAFT(p,i,Z));fflush(stdout);
      }
        
  }
 
  bbnfc=nfc;
  bb_tot_pot_energy = tot_pot_energy;
  /*  nfc = 0; */
  bbminsteps_tmp = bbminsteps_max;
  //printf("myid=%d myrank: %d bbnfc=%d, nfc = %d, bb_tot_pot_energy = %f, tot_pot_energy = %f: \n",  myid,myrank,bbnfc, nfc,bb_tot_pot_energy,tot_pot_energy);fflush(stdout); 
  /* minimization loop */
  for (bbminsteps= 0; bbminsteps <= bbminsteps_tmp; ++bbminsteps) {

#ifdef NBLIST
    check_nblist();
#else
    fix_cells();  
#endif

#ifdef CG
        if (bb_ensemble == ENS_CG) reset_cg();
#endif
 
#ifdef FBC
    update_fbc();
#endif

#ifdef TIMING
    imd_start_timer(&time_forces);
#endif
#if defined (CG) && !defined(ACG)
    if (bb_ensemble == ENS_CG) cg_step(bbminsteps);
    else
#elif defined(ACG)
    if (bb_ensemble == ENS_CG) acg_step(bbminsteps);
    else
#endif

        /* calculation of forces */

        calc_forces(bbminsteps);

#ifdef EXTPOT
    calc_extpot();
#endif
#ifdef RIGID
    calc_superforces();
#endif    
#ifdef NEB
    calc_forces_neb();
#endif
#ifdef FEFL
    calc_fefl();
#endif
#ifdef TIMING
    imd_stop_timer(&time_forces);
#endif


#ifdef GLOK
    /* "global convergence": set momenta to 0 if P*F < 0 (global vectors) */
    if (bb_ensemble == ENS_GLOK) {
      update_glok();     
    }
#endif
 
#ifdef TIMING
    imd_start_timer(&time_integrate);
#endif
#if !defined(CBE) || !defined(SPU_INT)
    /* move atoms */    
    if (bb_ensemble != ENS_CG) bb_move_atoms(); /* here PxF is recalculated */
#endif
#ifdef TIMING
    imd_stop_timer(&time_integrate);
#endif

    /* Periodic I/O */
#ifdef TIMING
    imd_start_timer(&time_output);
#endif
    //printf(" bbminsteps = %d \n",bbminsteps);fflush(stdout);
    //if ((checkpt_int > 0) && (0 == bbminsteps % checkpt_int)) 
       //write_bb_config( bbminsteps/checkpt_int, bbminsteps);
    if ((eng_int  > 0) && (0 == bbminsteps % eng_int )) {
        write_bbeng_file(bbminsteps);
    }

#ifdef TIMING
    imd_stop_timer(&time_output);
#endif


#ifdef RELAX
    check_bbrelaxed();
#endif

#ifdef NBLIST
    check_nblist();
#else
    fix_cells();  
#endif



    /* write checkpoint, if empty write file is found */
    if ((watch_int > 0) && (0==bbminsteps%watch_int)) check_write();

    /* finish, if stop file is found */
    if ((stop_int > 0) && (0==bbminsteps%stop_int)) {
      if ((finished = check_stop())) break;
    }

      /* finish, if max deformation bbminsteps in quasistatic simulation are done */
#ifdef RELAX
#if defined (DEFORM) || defined (HOMDEF) || defined (EXTPOT) || defined (FBC)
    if ( (bb_ensemble==ENS_MIK) || (bb_ensemble==ENS_GLOK) || (bb_ensemble==ENS_CG) ) {
        if ( (max_sscount>0) && (sscount>max_sscount) ) {
                 finished = 1 ;
                 bbminsteps_tmp = bbminsteps;
                 break;
            }
    }
#endif
#endif
    

    
    /* finish, if maxwalltime is reached */
    if (maxwalltime > 0) {
      if ((finished = check_walltime())) break;
    }
  } /* end of md loop*/
  tot_pot_energy = bb_tot_pot_energy;
  /* clean up the current phase, and clear restart flag */
  imdrestart=0;
  if (0==myid) {
    write_itr_file(-1, bbminsteps_tmp,"");
    if (0==myrank) {
       printf("bbminsteps_tmp = %d \n", bbminsteps_tmp);
       printf( "End of bb_simulation %d\n", simulation );
       if (bbminsteps_tmp >=  bbminsteps_max) {
           error("bbminsteps_tmp >  bbminsteps_max");
       }
    }
  }  
  return finished;
}
    



/*****************************************************************************
*
*  check if sample is relaxed
*
*****************************************************************************/

#ifdef RELAX

void check_bbrelaxed(void)
{
    int k;
    int write_ss=1;
    is_relaxed = 0;

    if ((bb_ensemble==ENS_MIK) || (bb_ensemble==ENS_GLOK) || (bb_ensemble==ENS_CG)) {
        
        int stop = 0;
        real fnorm2, ekin, epot, delta_epot;
#ifdef NEB
        MPI_Allreduce( &fnorm, &neb_fnorm, 1, REAL, MPI_SUM, MPI_COMM_WORLD);
        neb_fnorm = SQRT( neb_fnorm / (nactive * (neb_nrep-2)) );
        if (neb_fnorm < fnorm_threshold) is_relaxed = 1;
        else is_relaxed = 0;
#else
        fnorm2 = SQRT( fnorm / nactive );
        ekin   = 2 * tot_kin_energy / nactive;
        epot   = tot_pot_energy / natoms;
        delta_epot = old_epot - epot;
        if (delta_epot < 0) delta_epot = -delta_epot;
        
        if ((ekin  <  ekin_threshold) || (fnorm2 < fnorm_threshold) || 
            (delta_epot < delta_epot_threshold)) is_relaxed = 1;
        else is_relaxed = 0;
        
        old_epot = epot;
#endif


        if (is_relaxed) {
            nfc=bbnfc;
            stop = 1;
            //write_bbeng_file(bbminsteps);
       //   if(write_ss==1)
            //write_ssconfig(steps);
            ref2box_x.x = box_x.x;
            ref2box_y.y = box_y.y;
            ref2box_z.z = box_z.z;
            box_x.x = tmp_box_x.x;
            box_y.y = tmp_box_y.y;
            box_z.z = tmp_box_z.z;            
            /* return recorded atoms imformation after doing relaxation */
           for (k=0; k<NCELLS; ++k) { /* loop over all cells */
                int  i,j,sort,num;
                cell *p;
                p = CELLPTR(k);
               for (i=0; i<p->n; ++i) { /* loop over all atoms in the cell */
                    num = NUMMER(p,i);
                    REFPOSTWO(p,i,X)   = ORT(p,i,X);
                    REFPOSTWO(p,i,Y)   = ORT(p,i,Y);
                    REPOST[num*SDIM]   = REFPOSTWO(p,i,X);
                    REPOST[num*SDIM+1] = REFPOSTWO(p,i,Y);
                    ORT(p,i,X)         = TMPPOS(p,i,X);
                    ORT(p,i,Y)         = TMPPOS(p,i,Y);
                    KRAFT(p,i,X)       = TMPKRAFT(p,i,X);
                    KRAFT(p,i,Y)       = TMPKRAFT(p,i,Y);
#ifndef TWOD
                    REFPOSTWO(p,i,Z)   = ORT(p,i,Z);
                    REPOST[num*SDIM+2] = REFPOSTWO(p,i,Z);
                    ORT(p,i,Z)         = TMPPOS(p,i,Z);
                    KRAFT(p,i,Z)       = TMPKRAFT(p,i,Z);
#endif
               }
           } 

            /* if we are doing quasistatic simulations, we write out     *
             * in case the sample is relaxed or the max. nr of relaxation*
             * steps has been reached. Here we are only writing out if   *
             * the sample is truely relaxed */
#ifdef DEFORM
            if (max_deform_int > 0) write_ss=0; 
#endif
#ifdef HOMDEF
            if (lindef_int > 0) write_ss=0;
#endif
#ifdef FBC
            if (have_fbc_incr) write_ss=0;
#endif
#ifdef EXTPOT
            if (ep_max_int > 0) write_ss=0;
#endif
      //      if(write_ss==1)
      //          write_ssconfig(bbminsteps);
            
#ifdef NEB
            if (0==myrank) write_neb_eng_file(bbminsteps);
#else
            
            if (0==myid) {
                bb_epot = epot;
                printf("nfc = %d epot = %22.16f\n", nfc, epot );
                printf("ekin = %e fnorm = %e f_max = %e delta_epot = %e\n", 
                       ekin, fnorm2, f_max, delta_epot);
            }
#endif
                       
        }

        
#ifdef DEFORM
    if (max_deform_int > 0) stop=0;
#endif
#ifdef HOMDEF
    if (lindef_int > 0) stop=0;
#endif
#ifdef FBC
    if (have_fbc_incr) stop=0;
#endif
#ifdef EXTPOT
    if (ep_max_int > 0) stop = 0;
#endif
    if (stop) bbminsteps_tmp = bbminsteps;

  }

}

#endif /* RELAX */


/******************************************************************************
*
*  main program of bond boost
*
******************************************************************************/

void calc_bondboost(int steps)
{
  int  i, b, k, n;
  real ad_bV, bst_fcr, Temp, epot;
  real *pptr, *rptr;

#ifdef UNIAX
  Temp = 2.0 * tot_kin_energy / (nactive + nactive_rot);
#elif defined(DAMP)
  Temp = 2.0 * tot_kin_energy / (nactive - n_damp);
#else
  Temp = 2.0 * tot_kin_energy / nactive;
#endif
  epot = tot_pot_energy / natoms;
  //bb_ad_bV     = 0.0;
  ad_bV        = 0.0;
  Nb           = 0.0;
  A_e_max      = 0.0;
  bb_sum_bV    = 0.0;
  tmp_maxbond1 = 0.0;
  tmp_maxbond2 = 0.0; 
  //printf(" steps = %d \n",steps);fflush(stdout);
  /* update total required atoms positions for the use of bb_neigh */
  for (k=0; k<NCELLS; ++k) { 
      int  i,j,sort,num;
      cell *p;
      p = CELLPTR(k);
      for (i=0; i<p->n; ++i) {    
          num = NUMMER(p,i);
          NOWPOS[(num)*SDIM]   = ORT(p,i,X);
          NOWPOS[(num)*SDIM+1] = ORT(p,i,Y);
          REACFJ[(num)*SDIM]   = 0.0;
          REACFJ[(num)*SDIM+1] = 0.0;
#ifndef TWOD
          NOWPOS[(num)*SDIM+2] = ORT(p,i,Z);
          REACFJ[(num)*SDIM+2] = 0.0;
#endif
      } 
  } 

  /* first loop for the calculation of bias(boost) potential and the count of tagged bond number */
  /* AR is not used in this version */
#ifndef SBOOST
  for (k=0; k<ncells; ++k) {
      do_bb_pot(cell_array+k);
  }
#else
  for (k=0; k<ncells; ++k) {
      do_sbb_pot(cell_array+k);
  }
#endif

#ifdef MPI
  MPI_Allreduce( &tmp_maxbond1, &eci_ref1max, 1, REAL, MPI_MAX, cpugrid); 
#else
  eci_ref1max = tmp_maxbond1;
#endif
  bb_ref1max  = eci_ref1max;
  //printf(" bb_ref1max = %f, eci_ref1max = %f \n",bb_ref1max,eci_ref1max);fflush(stdout); 
  //max_epot    = MAX(epot,max_epot);
  //if (epot > max_epot) {  /* for capturing the strain at saddle point by watching largest potential energy */
    // max_epot = epot;
    // saddleta = bb_ref1max;
  //}
     cons_V   = 0.0;
  if (Nb != 0.0) 
     cons_V   = bb_tot_bV/Nb;
  bb_sum_bV   = cons_V*bb_sum_bV;
  bb_eps2     = bb_epscrit * bb_epscrit;
#ifdef ADAPTIVE
  if (bflag1 == 0) {
     bb_relaxsteps++;
     bflag3  = 0;
     if (bb_relaxsteps >= bb_relaxsteps_max) {
	bflag3  = 1;
	bb_minimize(1);
	bb_relaxsteps = 0;
     }
  }
  else if(bflag1 == 1) {
     if (bflag2 == 0) {
         //bb_maxwellrelax++;         
        bflag3     = 0;
        A_e_max    = 0.0;
        if (bb_ref1max <= bb_epscrit) {
            //printf(" bb_ref1max = %f, and bb_epscrit = %f \n", bb_ref1max,bb_epscrit);fflush(stdout);
	   C_x     = 1.0 - (bb_ref1max/bb_epscrit)*(bb_ref1max/bb_epscrit);
	   C_x2    = C_x * C_x;
	   B_x     = 1.0 - bbp1_2*(bb_ref1max/bb_epscrit)*(bb_ref1max/bb_epscrit);
	   B_x2    = B_x * B_x;
	   A_e_max = C_x2/B_x;
	   bb_under++;
	   bb_shdn = 0;
	   if (bb_under >= bb_under_max) {
	      bflag3   = 1;
	      bb_under = 0;
	   }
        }
        else if (bb_ref1max > bb_epscrit) {
            //printf(" bb_ref1max = %f, and bb_epscrit = %f \n", bb_ref1max,bb_epscrit);fflush(stdout);
	   bb_shdn++;
	   if (bb_shdn >= bb_shdn_max) {
	      bflag3   = 1;
 	      bb_minimize(2);
	      bb_shdn  = 0;
	   }
        }
        ad_bV     = A_e_max * bb_sum_bV;
        //printf(" ad_bV = %f \n", ad_bV);fflush(stdout);
        //bb_ad_bV  = epot + ad_bV;
     }
     else if (bflag2 == 1) {
        bb_btime++;
        bflag3     = 0;
        A_e_max    = 0.0;
        if (bb_ref1max <= bb_epscrit) {
	   C_x     = 1.0 - (bb_ref1max/bb_epscrit)*(bb_ref1max/bb_epscrit);
	   C_x2    = C_x * C_x;
	   B_x     = 1.0 - bbp1_2*(bb_ref1max/bb_epscrit)*(bb_ref1max/bb_epscrit);
	   B_x2    = B_x * B_x;
	   A_e_max = C_x2/B_x;
	   bb_shdn = 0;
        }
        else if (bb_ref1max > bb_epscrit) {
	   bb_shdn++;
	   if (bb_shdn >= bb_shdn_max) {
	      bflag3   = 1;
 	      bb_minimize(3);
	      bb_shdn  = 0;
	   }
        }
        ad_bV     = A_e_max * bb_sum_bV;
        bst_fcr   = exp(ad_bV/Temp);
        sum_bfcr += bst_fcr;
     }
  }
#else
  if (bflag1 == 0) {
     bb_relaxsteps++;
     bflag3  = 0;
     if (bb_relaxsteps >= bb_relaxsteps_max) {
	bflag3  = 1;
	bb_minimize(1);
	bb_relaxsteps = 0;
     }
  }
  else if(bflag1 == 1) {
     bb_btime++;
     bflag3     = 0;
     A_e_max    = 0.0;
     if (bb_ref1max <= bb_epscrit) {
        C_x     = 1.0 - (bb_ref1max/bb_epscrit)*(bb_ref1max/bb_epscrit);
        C_x2    = C_x * C_x;
        B_x     = 1.0 - bbp1_2*(bb_ref1max/bb_epscrit)*(bb_ref1max/bb_epscrit);
        B_x2    = B_x * B_x;
        A_e_max = C_x2/B_x;
        bb_shdn = 0;
     }
     else if (bb_ref1max > bb_epscrit) {
        bb_shdn++;
	if (bb_shdn >= bb_shdn_max) {
	   bflag3   = 1;
 	   bb_minimize(3);
	   bb_shdn  = 0;
        }
     }
     ad_bV     = A_e_max * bb_sum_bV;
     bst_fcr   = exp(ad_bV/Temp);
     sum_bfcr += bst_fcr;
     //printf(" bb_btime = %f, bst_fcr = %f, sum_bfcr = %f\n",bb_btime,bst_fcr,sum_bfcr);fflush(stdout);
  }  
#endif
  
  bb_ad_bV  = epot + ad_bV/natoms;
#ifndef SBOOST 
  for (k=0; k<ncells; ++k) {
      do_bb_forc(cell_array+k);
  }
#else
  for (k=0; k<ncells; ++k) {
      do_sbb_forc(cell_array+k);
  }
#endif
  
#ifdef MPI
  MPI_Allreduce( &tmp_maxbond2, &eci_ref2max, 1, REAL, MPI_MAX, cpugrid);
#else
  eci_ref2max = tmp_maxbond2;
#endif

  bb_ref2max = eci_ref2max;

  /* output of the information of bond ratio, strain */
  if ((eng_int  > 0) && (0 == steps % eng_int )) write_bbtran_file(steps);

}



/******************************************************************************
*
*  the post processes are finished in this function 
*
******************************************************************************/
#ifdef ADAPTIVE
void postpro_boost(int steps)
{
  int k, n;
  if (bflag1 == 0) {         /* annealing in regular MD */
     if ((bflag3 == 1) && (bb_ref2max > bb_bound_div)) { /* transition even occurs in regular MD, seeking next well */
        for (k=0; k<NCELLS; ++k) { /* loop over all cells */
            int  i,j, sort;
            cell *p;
            p = CELLPTR(k);
            for (i=0; i<p->n; ++i) { /* loop over all atoms in the cell */
                REFPOSONE(p,i,X) = REFPOSTWO(p,i,X);
                REFPOSONE(p,i,Y) = REFPOSTWO(p,i,Y);
                REFPOSONE(p,i,Z) = REFPOSTWO(p,i,Z);
	        BBNEIGH(p,i)->nbondsref1 = 0;            
            }
        }
        for (n=0; n<nlists; ++n) {
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
#endif
            for (k=0; k<npairs[n]; ++k) { //???
                vektor pbc;
                pair *P;
                P = pairs[n]+k;
                pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
                pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
                pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
                do_bb_neightab(cell_array + P->np, cell_array + P->nq, pbc);
            }
        }
     }  
     else if ((bflag3 == 1) && (bb_ref2max <= bb_bound_div)) { /* fall into a well, need to launch a boost MD */
        bb_epsold    = bb_epscrit;
        bb_btimeold  = bb_btime;
        sum_bfcrold  = sum_bfcr;
        bb_tot_bVold = bb_tot_bV;
        //min_saddleta = 100.0;
        bb_overbs    = 0;
        bb_index     = 0;
        bflag1       = 1;
        bflag2       = 0;
        bflag3       = 0;
        bb_box_x.x   = box_x.x;
        bb_box_y.y   = box_y.y;
        bb_box_z.z   = box_z.z;        
        for (k=0; k<NCELLS; ++k) { /* loop over all cells */
            int  i,j, sort;
            cell *p;
            p = CELLPTR(k);
            for (i=0; i<p->n; ++i) { /* loop over all atoms in the cell */
                REFPOSONE(p,i,X) = REFPOSTWO(p,i,X);
                REFPOSONE(p,i,Y) = REFPOSTWO(p,i,Y);
                REFPOSONE(p,i,Z) = REFPOSTWO(p,i,Z);
                OLDPOS(p,i,X)    = ORT(p,i,X);
                OLDPOS(p,i,Y)    = ORT(p,i,Y);
                OLDPOS(p,i,Z)    = ORT(p,i,Z);
                OLDIMPULS(p,i,X) = IMPULS(p,i,X);
                OLDIMPULS(p,i,Y) = IMPULS(p,i,Y);
                OLDIMPULS(p,i,Z) = IMPULS(p,i,Z);
	        BBNEIGH(p,i)->nbondsref1 = 0;            
            }
        }
        for (n=0; n<nlists; ++n) {
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
#endif
            for (k=0; k<npairs[n]; ++k) { //???
                vektor pbc;
                pair *P;
                P = pairs[n]+k;
                pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
                pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
                pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
                do_bb_neightab(cell_array + P->np, cell_array + P->nq, pbc);
            }
        }
        if (bb_evens == 0) {
           write_config_select( bb_evens/bb_checkpt_int, "bb_evens", write_atoms_bb, write_header_bb);
           write_bbcheck_file(steps);
        }
     }	
  }
  else if (bflag1 == 1) {    /* under boost MD */
     if (bflag2 == 0) {      /* mode 1: searching largest bb_epscrit */
	if ((bflag3 == 1) && (bb_ref1max <= bb_epscrit)) { /* under boosting for a long time, we need to enlarge bias voltage */
           bflag1 = 0;
           bb_epscrit = bb_epsold;
           box_x.x    = bb_box_x.x;
           box_y.y    = bb_box_y.y;
           box_z.z    = bb_box_z.z;
           for (k=0; k<NCELLS; ++k) { /* let system be boosted from ref1 */
               int  i,j,sort;
               cell *p;
               p = CELLPTR(k);
               for (i=0; i<p->n; ++i) {                  
                   ORT(p,i,X) = REFPOSONE(p,i,X);
                   ORT(p,i,Y) = REFPOSONE(p,i,Y);
                   ORT(p,i,Z) = REFPOSONE(p,i,Z);
               } 
           } 
           maxwell(temperature); /* Initialize the velocities */
            //bb_epscrit += 0.05;
	}
	else if ((bflag3 == 1) && (bb_ref1max > bb_epscrit)) { /* stoping boosting for a long time, checking transition even */
           if (bb_ref2max <= bb_bound_div) {      /* transition has not happened yet, enlarge bb_epscrit */
              bb_epscrit += delta_bb_epscrit;
           }
           else if (bb_ref2max > bb_bound_div) {  /* transition has happened */
              bb_qsafety++;
              bb_evens++;
              bflag1 = 0;
              //bb_maxwellrelax = 0;
              if (bb_epot < min_epot) {
                 bflag4   = 1;
                 min_epot = bb_epot;
                 tmp2_box_x.x = ref2box_x.x;
                 tmp2_box_y.y = ref2box_y.y;
                 tmp2_box_z.z = ref2box_z.z;
                 for (k=0; k<NCELLS; ++k) { /* let system be boosted from ref1 */
                     int  i,j,sort;
                     cell *p;
                     p = CELLPTR(k);
                     for (i=0; i<p->n; ++i) {                  
                         TMP2POS(p,i,X) = REFPOSTWO(p,i,X);
                         TMP2POS(p,i,Y) = REFPOSTWO(p,i,Y);
                         TMP2POS(p,i,Z) = REFPOSTWO(p,i,Z);
                     } 
                 }
              }
                 write_bbcheck_file(steps);
              if (bb_qsafety < bb_qsafety_max) {
                  //bb_epscrit = bb_epscrit_old;
                 bb_epscrit = bb_epsold;
                 box_x.x    = bb_box_x.x;
                 box_y.y    = bb_box_y.y;
                 box_z.z    = bb_box_z.z;                 
                 for (k=0; k<NCELLS; ++k) { /* let system be boosted from ref1 */
                     int  i,j,sort;
                     cell *p;
                     p = CELLPTR(k);
                     for (i=0; i<p->n; ++i) {                  
                         ORT(p,i,X) = REFPOSONE(p,i,X);
                         ORT(p,i,Y) = REFPOSONE(p,i,Y);
                         ORT(p,i,Z) = REFPOSONE(p,i,Z);
                     } 
                 } 
                 maxwell(temperature); /* Initialize the velocities */
                 write_config_select( bb_evens/bb_checkpt_int, "bb_evens", write_atoms_bb, write_header_bb);
              }
              else if (bb_qsafety >= bb_qsafety_max) {
                 //bb_epscrit = min_saddleta*0.9; /* fix a largest bb_epscrit */
		 //bflag2     = 1;
                 //bb_maxwellrelax = 0;
                 //bb_tot_bV *= 0.5;
                 //box_x.x    = bb_box_x.x;
                 //box_y.y    = bb_box_y.y;
                 //box_z.z    = bb_box_z.z;
                 //for (k=0; k<NCELLS; ++k) {
                   //  int  i,j, sort;
                     //cell *p;
                     //p = CELLPTR(k);
                     //for (i=0; i<p->n; ++i) { 
                       //  ORT(p,i,X)    = OLDPOS(p,i,X);
                       //  ORT(p,i,Y)    = OLDPOS(p,i,Y);
                       //  ORT(p,i,Z)    = OLDPOS(p,i,Z);
                       //  IMPULS(p,i,X) = OLDIMPULS(p,i,X);
                       //  IMPULS(p,i,Y) = OLDIMPULS(p,i,Y);
                       //  IMPULS(p,i,Z) = OLDIMPULS(p,i,Z);           
                    // } /* loop over all atoms in the cell */
                 //} /* loop over all cells */
                 write_config_select( bb_evens/bb_checkpt_int, "bb_evens", write_atoms_bb, write_header_bb);
                 if (bflag4 == 1) {
                    bb_epscrit = bb_epsold;
		    bb_qsafety = 0;
                    bflag4     = 0;
                    box_x.x    = tmp2_box_x.x;
                    box_y.y    = tmp2_box_y.y;
                    box_z.z    = tmp2_box_z.z;
                    for (k=0; k<NCELLS; ++k) {
                        int  i,j, sort;
                        cell *p;
                        p = CELLPTR(k);
                        for (i=0; i<p->n; ++i) { 
                            ORT(p,i,X)    = TMP2POS(p,i,X);
                            ORT(p,i,Y)    = TMP2POS(p,i,Y);
                            ORT(p,i,Z)    = TMP2POS(p,i,Z);
                        } /* loop over all atoms in the cell */
                    } /* loop over all cells */
                    maxwell(temperature); /* Initialize the velocities */
                 }
                 else if (bflag4 == 0) {
                    error(" transition occur bflag4 == 0 !!! \n");
                 }
              }
              //write_config_select( bb_evens/bb_checkpt_int, "bb_evens", write_atoms_bb, write_header_bb);
           }
	}
     }
     else if (bflag2 == 1) {  /* mode 2: running boost MD with a safe largest bb_epscrit, and calculate boost factor */
        if ((bflag3 == 1) && (bb_ref2max > bb_bound_div)) {        /* transition occurs! */
           bb_evens++;
           boost_fcr  = sum_bfcr;
           bb_tot_bV  = bb_tot_bVold;
           bb_epscrit = bb_epsold;
           bb_btime   = 0;
           sum_bfcr   = 0.0;
           bflag1     = 0;
           bflag2     = 0;
           for (k=0; k<NCELLS; ++k) { /* loop over all cells */
               int  i,j, sort;
               cell *p;
               p = CELLPTR(k);
               for (i=0; i<p->n; ++i) { /* loop over all atoms in the cell */
                   REFPOSONE(p,i,X) = REFPOSTWO(p,i,X);
                   REFPOSONE(p,i,Y) = REFPOSTWO(p,i,Y);
                   REFPOSONE(p,i,Z) = REFPOSTWO(p,i,Z);
	           BBNEIGH(p,i)->nbondsref1 = 0;
               }
           }
           for (n=0; n<nlists; ++n) {
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
#endif
               for (k=0; k<npairs[n]; ++k) { //???
                   vektor pbc;
                   pair *P;
                   P = pairs[n]+k;
                   pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
                   pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
                   pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
                   do_bb_neightab(cell_array + P->np, cell_array + P->nq, pbc);
               }
           }
           //write_config_select( bb_evens/bb_checkpt_int, "bb_evens", write_atoms_bb, write_header_bb);
           error(" transition occur !!! \n");
        }
	else if ((bflag3 == 1) && (bb_ref2max <= bb_bound_div)) {  /* watching window is too short for transition */
           error(" watching window is too short for transition !!! \n");
	}
     }
  }
}
#else
void postpro_boost(int steps)
{
  int k, n;
  if (bflag1 == 0) {         /* annealing in regular MD */
     if ((bflag3 == 1) && (bb_ref2max > bb_bound_div)) { /* transition even occurs in regular MD, seeking next well */
        for (k=0; k<NCELLS; ++k) { /* loop over all cells */
            int  i,j, sort;
            cell *p;
            p = CELLPTR(k);
            for (i=0; i<p->n; ++i) { /* loop over all atoms in the cell */
                REFPOSONE(p,i,X) = REFPOSTWO(p,i,X);
                REFPOSONE(p,i,Y) = REFPOSTWO(p,i,Y);
                REFPOSONE(p,i,Z) = REFPOSTWO(p,i,Z);
	        BBNEIGH(p,i)->nbondsref1 = 0;            
            }
        }
        for (n=0; n<nlists; ++n) {
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
#endif
            for (k=0; k<npairs[n]; ++k) { //???
                vektor pbc;
                pair *P;
                P = pairs[n]+k;
                pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
                pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
                pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
                do_bb_neightab(cell_array + P->np, cell_array + P->nq, pbc);
            }
        }
     }  
     else if ((bflag3 == 1) && (bb_ref2max <= bb_bound_div)) { /* fall into a well, need to launch a boost MD */
        bb_epsold    = bb_epscrit;
        bb_btimeold  = bb_btime;
        sum_bfcrold  = sum_bfcr;
        bb_tot_bVold = bb_tot_bV;
        //min_saddleta = 100.0;
        bb_overbs    = 0;
        bb_index     = 0;
        bflag1       = 1;
        bflag2       = 0;
        bflag3       = 0;
        bb_box_x.x   = box_x.x;
        bb_box_y.y   = box_y.y;
        bb_box_z.z   = box_z.z;        
        for (k=0; k<NCELLS; ++k) { /* loop over all cells */
            int  i,j, sort;
            cell *p;
            p = CELLPTR(k);
            for (i=0; i<p->n; ++i) { /* loop over all atoms in the cell */
                REFPOSONE(p,i,X) = REFPOSTWO(p,i,X);
                REFPOSONE(p,i,Y) = REFPOSTWO(p,i,Y);
                REFPOSONE(p,i,Z) = REFPOSTWO(p,i,Z);
                OLDPOS(p,i,X)    = ORT(p,i,X);
                OLDPOS(p,i,Y)    = ORT(p,i,Y);
                OLDPOS(p,i,Z)    = ORT(p,i,Z);
                OLDIMPULS(p,i,X) = IMPULS(p,i,X);
                OLDIMPULS(p,i,Y) = IMPULS(p,i,Y);
                OLDIMPULS(p,i,Z) = IMPULS(p,i,Z);
	        BBNEIGH(p,i)->nbondsref1 = 0;            
            }
        }
        for (n=0; n<nlists; ++n) {
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
#endif
            for (k=0; k<npairs[n]; ++k) { //???
                vektor pbc;
                pair *P;
                P = pairs[n]+k;
                pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
                pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
                pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
                do_bb_neightab(cell_array + P->np, cell_array + P->nq, pbc);
            }
        }
        if (bb_evens == 0) {
           write_config_select( bb_evens/bb_checkpt_int, "bb_evens", write_atoms_bb, write_header_bb);
           write_bbcheck_file(steps);
        }
     }	
  }
  else if (bflag1 == 1) {    /* under boost MD */
        if ((bflag3 == 1) && (bb_ref2max > bb_bound_div)) {        /* transition occurs! */
           bb_evens++;
           boost_fcr  = sum_bfcr/bb_btime;
           hypertime  = sum_bfcr*timestep;
           bb_tot_bV  = bb_tot_bVold;
           bb_epscrit = bb_epsold;
           printf(" bb_btime = %f, sum_bfcr = %f\n",bb_btime,sum_bfcr);fflush(stdout);
           printf(" boost_fcr = %f, hypertime = %f\n",boost_fcr,hypertime);fflush(stdout);
           bb_btime   = 0.0;
           sum_bfcr   = 0.0;
           bflag1     = 0;
           bflag2     = 0;
           for (k=0; k<NCELLS; ++k) { /* loop over all cells */
               int  i,j, sort;
               cell *p;
               p = CELLPTR(k);
               for (i=0; i<p->n; ++i) { /* loop over all atoms in the cell */
                   REFPOSONE(p,i,X) = REFPOSTWO(p,i,X);
                   REFPOSONE(p,i,Y) = REFPOSTWO(p,i,Y);
                   REFPOSONE(p,i,Z) = REFPOSTWO(p,i,Z);
	           BBNEIGH(p,i)->nbondsref1 = 0;
               }
           }
           for (n=0; n<nlists; ++n) {
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
#endif
               for (k=0; k<npairs[n]; ++k) { //???
                   vektor pbc;
                   pair *P;
                   P = pairs[n]+k;
                   pbc.x = P->ipbc[0]*box_x.x + P->ipbc[1]*box_y.x + P->ipbc[2]*box_z.x;
                   pbc.y = P->ipbc[0]*box_x.y + P->ipbc[1]*box_y.y + P->ipbc[2]*box_z.y;
                   pbc.z = P->ipbc[0]*box_x.z + P->ipbc[1]*box_y.z + P->ipbc[2]*box_z.z;
                   do_bb_neightab(cell_array + P->np, cell_array + P->nq, pbc);
               }
           }
           write_config_select( bb_evens/bb_checkpt_int, "bb_evens", write_atoms_bb, write_header_bb);
           write_bbcheck_file(steps);
           error(" transition occur !!! \n");
        }
	else if ((bflag3 == 1) && (bb_ref2max <= bb_bound_div)) {  /* watching window is too short for transition */
           error(" watching window is too short for transition or epscrit is too small !!! \n");
	}
  }
}
#endif

/******************************************************************************
*
*  cmoputer the boost potential for the satrin boost method
*
******************************************************************************/ 
#ifdef SBOOST
void do_sbb_pot(cell *p)
{
 /* static vektor *d  = NULL; */
  bb_neightab *bb_neigh;
  vektor pos;
  int    i, j, k, p_typ, j_typ, jnum, ii, inum;
  int    doboost=0;
  real   r, r2, req1, req1_2, req2, req2_2;
  real   ecij_ref1, *tmpptr;
  M3 V, W, J, eta, tmp;
  V3 dxji, dxji0;  
 
  /* For each atom in cell */
  for (i=0; i<p->n; ++i) {
      if ((VSORTE(p,i)== 0) || (VSORTE(p,i)== 1)) {
         M3ZERO(V);
         M3ZERO(W);  
         bb_neigh = BBNEIGH(p,i);  
     /* get the position of i atom at the present time step */
         pos.x = ORT(p,i,X);  
         pos.y = ORT(p,i,Y);
#ifndef TWOD
         pos.z = ORT(p,i,Z);
#endif
         inum  = NUMMER(p,i);
     /* for each atom in neighbouring cell according to bb_neigh */
         tmpptr = bb_neigh->vectref1;
         for (j=0; j<bb_neigh->nbondsref1; ++j) {
		
             jnum   = bb_neigh->numref1[j];
    //  printf(" VSORTE(p,i) = %d , NUMMER = %d, jnum = %d",VSORTE(p,i),NUMMER(p,i), jnum);fflush(stdout);
      /* calculate dxji[]=x_j[]-x_i[] and |dxji|^2 at present state */
             dxji[0] = -pos.x + NOWPOS[jnum*3]; 
             dxji[1] = -pos.y + NOWPOS[jnum*3+1]; 
#ifndef TWOD
             dxji[2] = -pos.z + NOWPOS[jnum*3+2];
#endif
      /* correction periodic boundary conditions */

	     if(dxji[0] >= box_x.x/2.0)
 	       { dxji[0] = dxji[0] - box_x.x;}
 	     else if(dxji[0] <= -box_x.x/2.0)
 	       { dxji[0] = dxji[0] + box_x.x;}
 	     if(dxji[1] >= box_y.y/2.0)
 	       { dxji[1] = dxji[1] - box_y.y;}
 	     else if(dxji[1] <= -box_y.y/2.0)
 	       { dxji[1] = dxji[1] + box_y.y;}
#ifndef TWOD
 	     if(dxji[2] >= box_z.z/2.0)
 	       { dxji[2] = dxji[2] - box_z.z;}
 	     else if(dxji[2] <= -box_z.z/2.0)
 	       { dxji[2] = dxji[2] + box_z.z;}
#endif
                   
             /* calculate dxji[] at ref_1 state */
 
             dxji0[0]  = *tmpptr++;
             dxji0[1]  = *tmpptr++;
             dxji0[2]  = *tmpptr++;
	     /* C[][] := a[]' * b[] */
             M3ASSIGNV3V3 (dxji0, dxji0, tmp);
	     /* B[][] := A[][] + B[][] */
             M3AdD (tmp, V);
             M3ASSIGNV3V3 (dxji0, dxji, tmp);
             M3AdD (tmp, W);

#ifdef DEBUG
           if (0==r2) { char msgbuf[256];
              sprintf(msgbuf, "Distance is zero between particles %d and %d!\n",
              NUMMER(p,i), jnum);
              error(msgbuf);
           }
#endif

         } /* for j */
	     /* M3VOLUME(A) fabs(M3DETERMINANT(A)) */
	     /* A[][] := A[][]^-1; return original det(A) */
	     /* M3MUL(A,B,C) => C[][] :=A[][]*B[][] */
             if (M3VOLUME(V) != 0) {
                 M3Inv (V);
                 M3MUL (V, W, J);
             }
	     else error(" V matrix is zero !!! \n");
             //else M3IDENTITY(J);
             // 2D matrix to 1D
             M3ToV9(inum,J,M3J);
             M3ToV9(inum,V,VInv);
	     /* C[][] := A[][] * A'[][] */
             M3MULT (J, eta);
	     /* A[][] := A[][] - a * I */
             M3SubdiaG (eta, 1);
	     /* A[][] := A[][] / divisor */
             M3DividE (eta, 2);
             Mises[inum]     = SymmetricM3MisesInvariant(eta);
             ecij_ref1       = Mises[inum];
             M3ToV9(inum,eta,M3eta);
             //printf("  inum = %d, SymmetricM3MisesInvariant = %f , Mises = %f\n",inum,SymmetricM3MisesInvariant(eta),Mises[inum]);fflush(stdout); 

             if (ecij_ref1 <= bb_epscrit) {
                 Nb++;
                 bb_sum_bV += (1.0 - (ecij_ref1/bb_epscrit)*(ecij_ref1/bb_epscrit));
             }
             tmp_maxbond1 = MAX(ecij_ref1,tmp_maxbond1);
             
      } /* for bb_vtype on inum of bb_vtype */
      else {
             inum  = NUMMER(p,i);
             M3ZERO(J);
             M3ZERO(V);
             M3ZERO(eta);
             M3ToV9(inum,J,M3J);
             M3ToV9(inum,V,VInv);
             M3ToV9(inum,eta,M3eta);
	     Mises[inum] = 0.0;
      } /* for bb_vtype on inum of not bb_vtype */
  } /* for i */
#ifdef DEBUG
  if (is_short==1) printf("\n Short distance in BB!\n");
#endif

} /* do_bb_pot */

/******************************************************************************
*
*  cmoputer the boost force for the bond boost method
*
******************************************************************************/
/*void do_bb_forc(cell *p, cell *q, vektor pbc, real *BB_ref1max, real *Cons_V)
{ */
void do_sbb_forc(cell *p)
{
  bb_neightab *bb_neigh;
  int i, j, jnum, inum;
  vektor pos, pos2, d, d1, force;
  real *tmpptr;
  real ecij_ref2;
  real eci, feci, grad_bf, d_Cx2, d_Bx, d_A_ema;
  V3 ddxji, dxji0, dMises, d_Vec, d_Vecj, dMisesjmax;
  V3 ddxij, dxij0, dxji2;
  M3 V, J, dW, dJ, dJx, dJy, dJz, eta, tmp, dWj, W2, tmp2;
  M3 detax, detay, detaz, dJj, dJjx, dJjy, dJjz, J2, eta2;

  /* for each atom in first cell */
  for (i=0; i<p->n; ++i) {

         V3ZERO(d_Vecj);
	 V3ZERO(dMisesjmax);
         M3ZERO(dW);
         M3ZERO(W2);
         bb_neigh = BBNEIGH(p,i);
         pos2.x = REFPOSTWO(p,i,X);  
         pos2.y = REFPOSTWO(p,i,Y);
#ifndef TWOD
         pos2.z = REFPOSTWO(p,i,Z);
#endif
         inum  = NUMMER(p,i);
         tmpptr = bb_neigh->vectref1;
         if (inum>natoms)
             error("BBOOST requires serial labeling of atoms from 0...natoms");
         for (j=0; j<bb_neigh->nbondsref1; ++j) {
	      jnum = bb_neigh->numref1[j];
             if (jnum>natoms)
                 error("BBOOST requires serial labeling of atoms from 0...natoms");

             dxji0[0] = *tmpptr++;
             dxji0[1] = *tmpptr++;
             dxji0[2] = *tmpptr++;
             ddxji[0] = -1;
             ddxji[1] = -1;    
             ddxji[2] = -1;
             M3ASSIGNV3V3 (dxji0, ddxji, tmp);
             M3AdD (tmp, dW);
	     /* the force of dMises[j]/dxi */
	     dxij0[0] = -dxji0[0];
	     dxij0[1] = -dxji0[1];
	     dxij0[2] = -dxji0[2];
             ddxij[0] = 1;
             ddxij[1] = 1;    
             ddxij[2] = 1;
	     M3ASSIGNV3V3 (dxij0, ddxij, dWj);
             V9ToM3(jnum,M3J,J);
             V9ToM3(jnum,VInv,V);
             V9ToM3(jnum,M3eta,eta);
	     /* M3MUL(A,B,C) => C[][] :=A[][]*B[][] */
	     M3MUL (V, dWj, dJj);
             M3ZERO(dJjx);
             M3ZERO(dJjy);
             M3ZERO(dJjz);
	     M3decomposing(dJj, dJjx, dJjy, dJjz);
	     /* C[][] := A[][] * B'[][] */
	     M3MULT1 (dJjx, J, detax);
	     M3MULT1 (dJjy, J, detay);
             M3MULT1 (dJjz, J, detaz);
	     /* A[][] := A[][] + A'[][] */
	     M3AdDT (detax);
	     M3AdDT (detay);
	     M3AdDT (detaz);
	     /* A[][] := A[][] / divisor */
	     M3DividE(detax,2);
	     M3DividE(detay,2);
	     M3DividE(detaz,2);
             dMises[0] = DecompoM3MisesInvariant(eta,detax);
             dMises[1] = DecompoM3MisesInvariant(eta,detay);
             dMises[2] = DecompoM3MisesInvariant(eta,detaz);
	     d_Vec[0] = -2.0*cons_V*dMises[0]/bb_eps2;
	     d_Vec[1] = -2.0*cons_V*dMises[1]/bb_eps2;
	     d_Vec[2] = -2.0*cons_V*dMises[2]/bb_eps2;
	     V3AdD(d_Vec,d_Vecj);
	     if (Mises[jnum] == bb_ref1max) {
                 dMisesjmax[0] = dMises[0];
                 dMisesjmax[1] = dMises[1];
                 dMisesjmax[2] = dMises[2];
	     }
 
	     if ((VSORTE(p,i)== 0) || (VSORTE(p,i)== 1)) {
                 dxji2[0] = -pos2.x + REPOST[jnum*3]; 
                 dxji2[1] = -pos2.y + REPOST[jnum*3+1]; 
#ifndef TWOD
                 dxji2[2] = -pos2.z + REPOST[jnum*3+2];
#endif
	         if(dxji2[0] >= ref2box_x.x/2.0)
 	          { dxji2[0] = dxji2[0] - ref2box_x.x;}
 	         else if(dxji2[0] <= -ref2box_x.x/2.0)
 	          { dxji2[0] = dxji2[0] + ref2box_x.x;}
 	         if(dxji2[1] >= ref2box_y.y/2.0)
 	          { dxji2[1] = dxji2[1] - ref2box_y.y;}
 	         else if(dxji2[1] <= -ref2box_y.y/2.0)
 	          { dxji2[1] = dxji2[1] + ref2box_y.y;}
#ifndef TWOD
 	         if(dxji2[2] >= ref2box_z.z/2.0)
 	          { dxji2[2] = dxji2[2] - ref2box_z.z;}
 	         else if(dxji2[2] <= -ref2box_z.z/2.0)
 	          { dxji2[2] = dxji2[2] + ref2box_z.z;}
#endif
                 M3ASSIGNV3V3 (dxji0, dxji2, tmp2);
                 M3AdD (tmp2, W2);
	     } /* for bb_vtype */             
	   
         } /* for j */
         V9ToM3(inum,M3J,J);
         V9ToM3(inum,VInv,V);
	 V9ToM3(inum,M3eta,eta);
         M3MUL (V, dW, dJ);
         M3ZERO(dJx);
         M3ZERO(dJy);
         M3ZERO(dJz);
         M3decomposing(dJ, dJx, dJy, dJz);
	 /* C[][] := A[][] * B'[][] */
	 M3MULT1 (dJx, J, detax);
	 M3MULT1 (dJy, J, detay);
         M3MULT1 (dJz, J, detaz);
	 /* A[][] := A[][] + A'[][] */
	 M3AdDT (detax);
	 M3AdDT (detay);
	 M3AdDT (detaz);
	 /* A[][] := A[][] / divisor */
	 M3DividE(detax,2);
	 M3DividE(detay,2);
	 M3DividE(detaz,2);
         dMises[0] = DecompoM3MisesInvariant(eta,detax);
         dMises[1] = DecompoM3MisesInvariant(eta,detay);
         dMises[2] = DecompoM3MisesInvariant(eta,detaz);

	 if((Mises[inum] == bb_ref1max)&&(bb_ref1max <= bb_epscrit) && (bflag1 == 1)) {

	 d_Vec[0] = -2.0*cons_V*dMises[0]/bb_eps2;
	 d_Vec[1] = -2.0*cons_V*dMises[1]/bb_eps2;
	 d_Vec[2] = -2.0*cons_V*dMises[2]/bb_eps2;

	 d_Cx2    = -4.0*C_x/bb_eps2;
	 d_Bx     = -2.0*bbp1_2/bb_eps2;
	 d_A_ema  = (d_Cx2*B_x-C_x2*d_Bx)/B_x2;

	 force.x  = -A_e_max * (d_Vec[0]+d_Vecj[0]) - d_A_ema*dMises[0]*bb_sum_bV;
         force.y  = -A_e_max * (d_Vec[1]+d_Vecj[1]) - d_A_ema*dMises[1]*bb_sum_bV;
         force.z  = -A_e_max * (d_Vec[2]+d_Vecj[2]) - d_A_ema*dMises[2]*bb_sum_bV;
	                
         KRAFT(p,i,X) += force.x; 
         KRAFT(p,i,Y) += force.y; 
#ifndef TWOD
         KRAFT(p,i,Z) += force.z;
#endif
	 }
	 else if ((Mises[inum] < bb_ref1max)&&(bb_ref1max <= bb_epscrit) && (bflag1 == 1)){
	 
	 d_Vec[0] = -2.0*cons_V*dMises[0]/bb_eps2;
	 d_Vec[1] = -2.0*cons_V*dMises[1]/bb_eps2;
	 d_Vec[2] = -2.0*cons_V*dMises[2]/bb_eps2;

	 d_Cx2    = -4.0*C_x/bb_eps2;
	 d_Bx     = -2.0*bbp1_2/bb_eps2;
	 d_A_ema  = (d_Cx2*B_x-C_x2*d_Bx)/B_x2;
         
	 force.x  = -A_e_max * (d_Vec[0]+d_Vecj[0]) - d_A_ema*dMisesjmax[0]*bb_sum_bV;
         force.y  = -A_e_max * (d_Vec[1]+d_Vecj[1]) - d_A_ema*dMisesjmax[1]*bb_sum_bV;
         force.z  = -A_e_max * (d_Vec[2]+d_Vecj[2]) - d_A_ema*dMisesjmax[2]*bb_sum_bV;
	                 
         KRAFT(p,i,X) += force.x; 
         KRAFT(p,i,Y) += force.y; 
#ifndef TWOD
         KRAFT(p,i,Z) += force.z;
#endif
	 } 
	     


	     if ((VSORTE(p,i)== 0) || (VSORTE(p,i)== 1)) {

             M3MUL (V, W2, J2);
             M3MULT (J2, eta2);
             M3SubdiaG (eta2, 1);
             M3DividE (eta2, 2);
             ecij_ref2   = SymmetricM3MisesInvariant(eta2);

             //printf("  inum = %d, SymmetricM3MisesInvariant2 = %f \n",inum,SymmetricM3MisesInvariant(eta2));fflush(stdout);

             tmp_maxbond2 = MAX(ecij_ref2,tmp_maxbond2);
	     } /* for bb_vtype */

  } /* for i */
#ifdef DEBUG
  if (is_short==1) printf("\n Short distance!\n");
#endif

}  /*  do_bb_forc  */


/* A[][] := A[][]^-1; return original det(A) */
double M3Inv (double A[3][3])
{
    double determinant, D11, D22, D33, D12, D23, D31, D13, D21, D32;
    D11 = A[1][1]*A[2][2]-A[1][2]*A[2][1];
    D22 = A[2][2]*A[0][0]-A[2][0]*A[0][2];
    D33 = A[0][0]*A[1][1]-A[0][1]*A[1][0];
    D12 = A[1][2]*A[2][0]-A[1][0]*A[2][2];
    D23 = A[2][0]*A[0][1]-A[2][1]*A[0][0];
    D31 = A[0][1]*A[1][2]-A[0][2]*A[1][1];
    D13 = A[1][0]*A[2][1]-A[2][0]*A[1][1];
    D21 = A[2][1]*A[0][2]-A[0][1]*A[2][2];
    D32 = A[0][2]*A[1][0]-A[1][2]*A[0][0];
    determinant = A[0][0]*D11+A[0][1]*D12+A[0][2]*D13;
    A[0][0] = D11/determinant;
    A[1][1] = D22/determinant;
    A[2][2] = D33/determinant;
    A[0][1] = D21/determinant;
    A[1][2] = D32/determinant;
    A[2][0] = D13/determinant;
    A[1][0] = D12/determinant;
    A[2][1] = D23/determinant;
    A[0][2] = D31/determinant;
    return (determinant);
} /* end M3Inv() */
#else

/******************************************************************************
*
*  cmoputer the boost potential for the bond boost method
*
******************************************************************************/ 

void do_bb_pot(cell *p)
{
 /* static vektor *d  = NULL; */
  bb_neightab *bb_neigh;
  vektor pos, d;
  int    i, j, k, p_typ, j_typ, jnum, ii, inum;
  int    doboost=0;
  real   r, r2, req1, req1_2, req2, req2_2;
  real   ecij_ref1, tmp;

  /* For each atom in cell */
  for (i=0; i<p->n; ++i) {
         bb_neigh = BBNEIGH(p,i);  
     /* get the position of i atom at the present time step */
         pos.x = ORT(p,i,X);  
         pos.y = ORT(p,i,Y);
#ifndef TWOD
         pos.z = ORT(p,i,Z);
#endif
         inum  = NUMMER(p,i);
     /* for each atom in neighbouring cell according to bb_neigh */

         for (j=0; j<bb_neigh->nbondsref1; ++j) {
		
             jnum   = bb_neigh->numref1[j];

             if ( (bb_boosttypes[VSORTE(p,i)]!= 0) ||  (bb_boosttypes[bb_vtypej[jnum]]!= 0) ) {             
      /* calculate distance at present state */
             d.x = pos.x - NOWPOS[jnum*3]; 
             d.y = pos.y - NOWPOS[jnum*3+1]; 
#ifndef TWOD
             d.z = pos.z - NOWPOS[jnum*3+2];
#endif
      /* correction periodic boundary conditions */

	     if(d.x >= box_x.x/2.0)
 	       { d.x =  d.x - box_x.x;}
 	     else if(d.x <= -box_x.x/2.0)
 	       { d.x =  d.x + box_x.x;}
 	     if(d.y >= box_y.y/2.0)
 	       { d.y =  d.y - box_y.y;}
 	     else if(d.y <= -box_y.y/2.0)
 	       { d.y =  d.y + box_y.y;}
#ifndef TWOD
 	     if(d.z >= box_z.z/2.0)
 	       { d.z =  d.z - box_z.z;}
 	     else if(d.z <= -box_z.z/2.0)
 	       { d.z =  d.z + box_z.z;}
#endif
             

             r2     = SPROD(d,d);		
             /* calculate distance at ref_1 state */
             req1_2 = bb_neigh->distref1[j];
          

#ifdef DEBUG
           if (0==r2) { char msgbuf[256];
              sprintf(msgbuf, "Distance is zero between particles %d and %d!\n",
              NUMMER(p,i), jnum);
              error(msgbuf);
           }
#endif

	     r    = sqrt(r2);
	     req1 = sqrt(req1_2);
 	     ecij_ref1 = fabs((r-req1)/req1);
	     if (ecij_ref1 <= bb_epscrit) {
	        Nb++;
	        bb_sum_bV += (1.0 - (ecij_ref1/bb_epscrit)*(ecij_ref1/bb_epscrit));
	     }
             tmp_maxbond1 = MAX(ecij_ref1,tmp_maxbond1);
           }
                  } /* for j */
  } /* for i */
#ifdef DEBUG
  if (is_short==1) printf("\n Short distance in BB!\n");
#endif

} /* do_bb_pot */


/******************************************************************************
*
*  cmoputer the boost force for the bond boost method
*
******************************************************************************/
void do_bb_forc(cell *p)
{
  bb_neightab *bb_neigh;
  int i, j, jnum, inum;
  vektor pos, pos2, d, d1, force, rij;
  real ecij_ref2, r, r2, req1, req1_2, req2, req2_2;
  real eci, feci, d_Vec, d_Cx2, d_Bx, d_A_ema, grad_bf, bb_tor;

  /* for each atom in first cell */
  for (i=0; i<p->n; ++i) {

         bb_neigh = BBNEIGH(p,i);  /* should we build a new BBNEIGH(p,i)? */
         pos.x  = ORT(p,i,X);
         pos2.x = REFPOSTWO(p,i,X);  
         pos.y  = ORT(p,i,Y);
         pos2.y = REFPOSTWO(p,i,Y);
#ifndef TWOD
         pos.z  = ORT(p,i,Z);
         pos2.z = REFPOSTWO(p,i,Z);
#endif
         inum  = NUMMER(p,i);
         if (inum>natoms)
             error("BBOOST requires serial labeling of atoms from 0...natoms");
         for (j=0; j<bb_neigh->nbondsref1; ++j) {
		
             jnum = bb_neigh->numref1[j];
             if (jnum>natoms)
                 error("BBOOST requires serial labeling of atoms from 0...natoms");

             if ( (bb_boosttypes[VSORTE(p,i)]!= 0) ||  (bb_boosttypes[bb_vtypej[jnum]]!= 0) ) {             

         /* calculate distance at present state */
             d.x  = pos.x - NOWPOS[jnum*3]; 
             d1.x = pos2.x - REPOST[jnum*3];

             d.y  = pos.y - NOWPOS[jnum*3+1]; 
             d1.y = pos2.y - REPOST[jnum*3+1];

#ifndef TWOD
             d.z  = pos.z - NOWPOS[jnum*3+2];
             d1.z = pos2.z - REPOST[jnum*3+2];
#endif
             /* correction periodic boundary conditions */

	     if(d.x >= box_x.x/2.0)
 	       { d.x =  d.x - box_x.x;}
 	     else if(d.x <= -box_x.x/2.0)
 	       { d.x =  d.x + box_x.x;}
 	     if(d.y >= box_y.y/2.0)
 	       { d.y =  d.y - box_y.y;}
 	     else if(d.y <= -box_y.y/2.0)
 	       { d.y =  d.y + box_y.y;}
#ifndef TWOD
 	     if(d.z >= box_z.z/2.0)
 	       { d.z =  d.z - box_z.z;}
 	     else if(d.z <= -box_z.z/2.0)
 	       { d.z =  d.z + box_z.z;}
#endif
             
             r2     = SPROD(d,d);
	     r      = sqrt(r2);
         /* calculate unit vector of rij */
             rij.x  = d.x/r;
             rij.y  = d.y/r;
#ifndef TWOD
             rij.z  = d.z/r;
#endif
 //      /* correction periodic boundary conditions */

       	     if(d1.x >= box_x.x/2.0)
 	       { d1.x =  d1.x - box_x.x;}
 	     else if(d1.x <= -box_x.x/2.0)
 	       { d1.x =  d1.x + box_x.x;}
 	     if(d1.y >= box_y.y/2.0)
 	       { d1.y =  d1.y - box_y.y;}
 	     else if(d1.y <= -box_y.y/2.0)
 	       { d1.y =  d1.y + box_y.y;}
#ifndef TWOD
 	     if(d1.z >= box_z.z/2.0)
 	       { d1.z =  d1.z - box_z.z;}
 	     else if(d1.z <= -box_z.z/2.0)
 	       { d1.z =  d1.z + box_z.z;}
#endif

             req2_2 = SPROD(d1,d1);

#ifdef DEBUG
        if (0==r2) { char msgbuf[256];
           sprintf(msgbuf, "Distance is zero between particles %d and %d!\n",
                NUMMER(p,i), jnum);
           error(msgbuf);
        }
#endif
 
             req1_2 = bb_neigh->distref1[j];
	     req1   = sqrt(req1_2);
	     req2   = sqrt(req2_2);

	     ecij_ref2    = fabs((req2-req1)/req1);
             tmp_maxbond2 = MAX(ecij_ref2,tmp_maxbond2);
	     eci          = (r-req1)/req1;
	     feci         = fabs(eci);
             bb_tor       = req1*bb_eps2;

	     if ((feci == bb_ref1max) && (bb_ref1max <= bb_epscrit) && (bflag1 == 1)) {

	        d_Vec   = -2.0*cons_V*eci/(bb_tor); 
	        if (eci > 0.0) {
		d_Cx2   = -4.0*C_x*bb_ref1max/(bb_tor);
	        d_Bx    = -2.0*bbp1_2*bb_ref1max/(bb_tor);
                }
	        else if (eci < 0.0) {
		d_Cx2   =  4.0*C_x*bb_ref1max/(bb_tor);
		d_Bx    =  2.0*bbp1_2*bb_ref1max/(bb_tor);
	        }

	        d_A_ema = (d_Cx2*B_x-C_x2*d_Bx)/B_x2;
	        grad_bf = - A_e_max * d_Vec - d_A_ema*bb_sum_bV;

                /* store force in temporary variable */
                force.x = rij.x * grad_bf;
                force.y = rij.y * grad_bf;
#ifndef TWOD
                force.z = rij.z * grad_bf;
#endif
                /* accumulate forces */
                KRAFT(p,i,X) += force.x; 
                KRAFT(p,i,Y) += force.y; 
#ifndef TWOD
                KRAFT(p,i,Z) += force.z;
#endif
                REACFJ[jnum*3]   -= force.x;
                REACFJ[jnum*3+1] -= force.y; 
#ifndef TWOD
                REACFJ[jnum*3+2] -= force.z;
#endif
             }
	     else if ((feci < bb_ref1max) && (bb_ref1max <= bb_epscrit) && (bflag1 == 1)) {

	        d_Vec   = -2.0*cons_V*eci/(bb_tor);
	        grad_bf = - A_e_max * d_Vec;

                /* store force in temporary variable */
                force.x = rij.x * grad_bf;
                force.y = rij.y * grad_bf;
#ifndef TWOD
                force.z = rij.z * grad_bf;
#endif
                /* accumulate forces */
                KRAFT(p,i,X) += force.x; 
                KRAFT(p,i,Y) += force.y; 
#ifndef TWOD
                KRAFT(p,i,Z) += force.z;
#endif
                REACFJ[jnum*3]   -= force.x;
                REACFJ[jnum*3+1] -= force.y; 
#ifndef TWOD
                REACFJ[jnum*3+2] -= force.z;
#endif
	     }
           }
         } /* for j */
         KRAFT(p,i,X) += REACFJ[inum*3]; 
         KRAFT(p,i,Y) += REACFJ[inum*3+1]; 
#ifndef TWOD
         KRAFT(p,i,Z) += REACFJ[inum*3+2];
#endif

  } /* for i */
#ifdef DEBUG
  if (is_short==1) printf("\n Short distance!\n");
#endif

}  /*  do_bb_forc  */
#endif
