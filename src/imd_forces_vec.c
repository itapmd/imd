

#define MAXCELL 500

/******************************************************************************
*
*  do_forces, version for vector processors
*
*  computes the forces between atoms in two given cells
*
******************************************************************************/

void do_forces(cell *ap, cell *aq, vektor apbc)
{

  int nmax;
  int jstart,jend;
  int i,j;
  int n;
  int c,k;
  int np,nq;
  vektor force[MAXCELL][MAXCELL];
  vektor tmp_force;
  real radius2;
  vektor d;
#ifdef MONOLJ
  real sig_d_rad2,sig_d_rad6,sig_d_rad12;
#else
  real pot_k0,pot_k1,pot_k2;
  real dv, d2v;
  int q_typ,p_typ;
  real *potptr;
  real chi;
#endif
  real pot_zwi,pot_grad;
  cell *p,*q;
  vektor pbc;
  real tmp_energy, tmp_virial;
#ifdef P_AXIAL
  vektor tmp_vir_vect;
#endif
  real inv_nq;
  int  inc;
  vektor p_force,q_force;

  /* Init force array and pot_eng array */
  for (i=0; i<MAXCELL; ++i) 
    for (j=0; j<MAXCELL; ++j) {
      force[i][j].x = 0.0;
      force[i][j].y = 0.0;
      force[i][j].z = 0.0;
    };


/*
   The double loop over the atoms of both cells in the risc version
   is converted to one (longer) single loop. The particle indexes i,j
   have to be calculated in the loop. The algorithm we use assumes that
   q->n > p->n. It generates more evenly distributed memory accesses (
   which is bad on machines with a cache, but fast on bank-interleaved
   systems).
*/

  /* Sx4 doesn't like writes to a function's arguments,
     so we copy them here to local variables. */

  if (aq->n > ap->n) {
    p = ap;
    q = aq;
    pbc.x = apbc.x;
    pbc.y = apbc.y;
    pbc.z = apbc.z; 
  } else {
    q = ap;
    p = aq;
    pbc.x = - apbc.x;
    pbc.y = - apbc.y;
    pbc.z = - apbc.z; 
  };


  np=p->n;
  nq=q->n;
  inv_nq = 1.0 / (real) nq;

  if ((MAXCELL<=np) || (MAXCELL<=nq)) error("Cell too large. Boost MAXCELL.");

  tmp_energy = 0.0;
  tmp_virial = 0.0;
#ifdef P_AXIAL
  tmp_vir_vect.x = 0.0;
  tmp_vir_vect.y = 0.0;
#ifndef TWOD
  tmp_vir_vect.z = 0.0;
#endif
#endif

#ifdef PVPCRAY
#pragma vdir ivdep
#endif
#ifdef SX4
#pragma vdir vector,nodep
#endif
  for (n=0; n<np*nq; ++n) {

    j   = n % nq;
    inc = (int) (n * inv_nq);
    i   = (j + inc) % np; 

    /* Calculate distance */
    d.x = q->ort X(j) - p->ort X(i) + pbc.x;
    d.y = q->ort Y(j) - p->ort Y(i) + pbc.y;
#ifndef TWOD
    d.z = q->ort Z(j) - p->ort Z(i) + pbc.z;
#endif        
    /* Square of distance */
    radius2 = SPROD(d,d);
    
    /* Calculate force, if distance smaller than cutoff */ 

   if ((radius2 <= r2_cut) &&  ((p!=q) || (i<j))) { 
#ifdef MONOLJ

      sig_d_rad2  = 2.0 / radius2;
      sig_d_rad6  = sig_d_rad2 * sig_d_rad2 * sig_d_rad2;
      sig_d_rad12 = sig_d_rad6 * sig_d_rad6;

      pot_grad = - 6 * sig_d_rad2 * ( sig_d_rad12 - sig_d_rad6 );
      pot_zwi  = sig_d_rad12 - 2.0 * sig_d_rad6;

#else
  
      /* Check for distances, shorter than minimal distance in pot. table */
      radius2 = (radius2 <= r2_0) ? r2_0 : radius2;

      /* Indices into potential table */
      k     = (int) ((radius2 - r2_0) * inv_r2_step);
      q_typ = q->sorte[j];
      p_typ = p->sorte[i];
      chi = (radius2 - r2_0 - k * r2_step) * inv_r2_step;

#ifdef STATIC_POT
      pot_k0 = potential[p_typ][q_typ][k    ];
      pot_k1 = potential[p_typ][q_typ][k + 1];
      pot_k2 = potential[p_typ][q_typ][k + 2];  
#else
      potptr = PTR_3D_V(potential, k, p_typ, q_typ , pot_dim);
      pot_k0 = *potptr; potptr += pot_dim.y * pot_dim.z;
      pot_k1 = *potptr; potptr += pot_dim.y * pot_dim.z;
      pot_k2 = *potptr;
#endif /* STATIC_POT */
      dv  = pot_k1 - pot_k0;
      d2v = pot_k2 - 2 * pot_k1 + pot_k0;
      
      /* Norm of Gradient */
      pot_grad = 2 * inv_r2_step * ( dv + (chi - 0.5) * d2v );
	
      /* Potential energy of atom */
      pot_zwi =  pot_k0 + chi * dv + 0.5 * chi * (chi - 1) * d2v;
#endif /* MONOLJ */
      
      /* Store forces in temp */
      force[i][j].x = d.x * pot_grad;
      force[i][j].y = d.y * pot_grad;
#ifndef TWOD
      force[i][j].z = d.z * pot_grad;
#endif
      if ((NUMMER(p,i)>=0) || (NUMMER(q,j)>=0))
        tmp_energy += pot_zwi;

#ifdef P_AXIAL
      tmp_vir_vect.x -= d.x * d.x * pot_grad;
      tmp_vir_vect.y -= d.y * d.y * pot_grad;
#ifndef TWOD
      tmp_vir_vect.z -= d.z * d.z * pot_grad;
#endif
#else
      tmp_virial     -= radius2 * pot_grad;  
#endif
	/* negativ, da pot_grad gleich force !! */
#ifdef STRESS_TENS
      p->presstens X(i) -= d.x * d.x * pot_grad;
      p->presstens Y(i) -= d.y * d.y * pot_grad;
      q->presstens X(i) -= d.x * d.x * pot_grad;
      q->presstens Y(i) -= d.y * d.y * pot_grad;
#ifndef TWOD
      p->presstens Z(i) -= d.z * d.z * pot_grad;
      q->presstens Z(i) -= d.z * d.z * pot_grad;
#endif
#endif
#ifdef NVX
        p->heatcond[i] += pot_zwi - radius2 * pot_grad;
        q->heatcond[j] += pot_zwi - radius2 * pot_grad;
#endif


    };
  };
  tot_pot_energy += tmp_energy;

#ifdef P_AXIAL
  vir_vect.x += tmp_vir_vect.x;
  vir_vect.y += tmp_vir_vect.y;
  virial     += tmp_vir_vect.x;
  virial     += tmp_vir_vect.y;
#ifndef TWOD
  vir_vect.z += tmp_vir_vect.z;
  virial     += tmp_vir_vect.z;
#endif
#else
  virial     += tmp_virial;
#endif


#ifdef SX4
#pragma vdir novector
#endif
  for (i=0; i<np; ++i) {
    tmp_force.x = 0.0;
    tmp_force.y = 0.0;
    tmp_force.z = 0.0;
#ifdef PVPCRAY
#pragma vdir shortloop
#endif
#ifdef SX4
#pragma vdir loopcnt=MAXCELL
#endif
    for (j=0; j<nq; ++j) {
      tmp_force.x += force[i][j].x;
      tmp_force.y += force[i][j].y;
      tmp_force.z += force[i][j].z;
    };
    p->kraft X(i) += tmp_force.x;
    p->kraft Y(i) += tmp_force.y;
    p->kraft Z(i) += tmp_force.z;
  };

#ifdef SX4
#pragma vdir novector
#endif
  for (j=0; j<nq; ++j) {
    tmp_force.x = 0.0;
    tmp_force.y = 0.0;
    tmp_force.z = 0.0;
#ifdef PVPCRAY
#pragma vdir shortloop
#endif
#ifdef SX4
#pragma vdir loopcnt=MAXCELL
#endif
    for (i=0; i<np; ++i) {
      tmp_force.x += force[i][j].x;
      tmp_force.y += force[i][j].y;
      tmp_force.z += force[i][j].z;
    };
    q->kraft X(j) -= tmp_force.x;
    q->kraft Y(j) -= tmp_force.y;
    q->kraft Z(j) -= tmp_force.z;
  };

  return;

}


