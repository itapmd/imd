
/******************************************************************************
*
* Monte-Carlo-Routinen von Franz Gaehler
*
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"
#include "potaccess.h"

/* potential energy difference between old and new position of an atom */

real mc_epot_diff( vektor old_pos, vektor new_pos, 
                   int p_num, int p_typ, ivektor cellc )
{
  int    q_typ, col, inc=ntypes*ntypes, is_short=0, i, j, k, l, m, n, r, s, t;
  vektor pbc, tmp_new, tmp_old, d_new, d_old;
  real   *qptr, pot_tmp, grad_tmp, radius2, result = 0.0;
  cell   *q;
  real   t_old=0.0,t_new=0.0,r2_short=0.0;
  ivektor ipbc;

  /* For all neighbours of this cell */
  for (l=-1; l <= 1; l++)
  for (m=-1; m <= 1; m++)
#ifndef TWOD
  for (n=-1; n <= 1; n++) 
#endif
  {
    r = cellc.x + l;  ipbc.x = 0;
    s = cellc.y + m;  ipbc.y = 0;
#ifndef TWOD
    t = cellc.z + n;  ipbc.z = 0;
#endif

    /* Apply periodic boundaries */

    if (r<0) {
       r = cell_dim.x-1; 
       ipbc.x--;
    }
    if (s<0) {
      s = cell_dim.y-1;
      ipbc.y--;
    }
#ifndef TWOD
    if (t<0) {
      t = cell_dim.z-1;
      ipbc.z--;
    }
#endif
    if (r>cell_dim.x-1) {
      r = 0; 
      ipbc.x++;
    }
    if (s>cell_dim.y-1) {
      s = 0; 
      ipbc.y++;
    }
#ifndef TWOD
    if (t>cell_dim.z-1) {
      t = 0; 
      ipbc.z++;
    }
#endif

    if (   ((pbc_dirs.x==1) || (pbc_dirs.x==ipbc.x))
#ifndef TWOD
        && ((pbc_dirs.z==1) || (pbc_dirs.z==ipbc.z)) 
#endif
        && ((pbc_dirs.y==1) || (pbc_dirs.y==ipbc.y)))
    {

#ifdef TWOD
      q = PTR_2D_V(cell_array,r,s,cell_dim);
      pbc.x = ipbc.x * box_x.x + ipbc.y * box_y.x; 
      pbc.y = ipbc.x * box_x.y + ipbc.y * box_y.y; 
#else
      q = PTR_3D_V(cell_array,r,s,t,cell_dim);
      pbc.x = ipbc.x * box_x.x + ipbc.y * box_y.x + ipbc.z * box_z.x;
      pbc.y = ipbc.x * box_x.y + ipbc.y * box_y.y + ipbc.z * box_z.y;
      pbc.z = ipbc.x * box_x.z + ipbc.y * box_y.z + ipbc.z * box_z.z;
#endif
      tmp_old.x = old_pos.x - pbc.x;
      tmp_new.x = new_pos.x - pbc.x;
      tmp_old.y = old_pos.y - pbc.y;
      tmp_new.y = new_pos.y - pbc.y;
#ifndef TWOD
      tmp_old.z = old_pos.z - pbc.z;
      tmp_new.z = new_pos.z - pbc.z;
#endif

      qptr   = q->ort;
      for (i = 0; i < q->n; i++) {
        d_new.x = *qptr - tmp_new.x;
        d_old.x = *qptr - tmp_old.x; ++qptr;
        d_new.y = *qptr - tmp_new.y;
        d_old.y = *qptr - tmp_old.y; ++qptr;
#ifndef TWOD
        d_new.z = *qptr - tmp_new.z;
        d_old.z = *qptr - tmp_old.z; ++qptr;
#endif
        if (p_num != q->nummer[i]) {

          q_typ  = q->sorte[i];
          col    = p_typ * ntypes + q_typ;

          /* contribution of new position */
          radius2 = SPROD(d_new,d_new);
          if (radius2 <= pair_pot.end[col]) {
            PAIR_INT2(pot_tmp,grad_tmp,pair_pot,col,inc,radius2,r2_short,is_short)
            t_new  += pot_tmp;
            result += pot_tmp;
          }

          /* contribution of old position */
          radius2 = SPROD(d_old,d_old);
          if (radius2 <= pair_pot.end[col]) {
            PAIR_INT2(pot_tmp,grad_tmp,pair_pot,col,inc,radius2,r2_short,is_short)
            t_old  += pot_tmp;
            result -= pot_tmp;
          }
        }
      }
    }
  }
  /*  printf("new %10.4e old %10.4e\n",t_new,t_old);  */
  return result;
}


/* potential energy of one atom */

real mc_epot_atom( vektor pos, int p_num, int p_typ, ivektor cellc )
{
  int    q_typ, col, inc=ntypes*ntypes, is_short=0, i, j, k, l, m, n, r, s, t;
  vektor pbc, tmp, d;
  real   *qptr, pot_tmp, grad_tmp, radius2, result = 0.0, r2_short=0;
  cell   *q;
  ivektor ipbc;

  /* For all neighbours of this cell */
  for (l=-1; l <= 1; l++)
  for (m=-1; m <= 1; m++)
#ifndef TWOD
  for (n=-1; n <= 1; n++) 
#endif
  {
    r = cellc.x + l;  pbc.x = 0;
    s = cellc.y + m;  pbc.y = 0;
#ifndef TWOD
    t = cellc.z + n;  pbc.z = 0;
#endif

    /* Apply periodic boundaries */

    if (r<0) {
       r = cell_dim.x-1; 
       ipbc.x--;
    }
    if (s<0) {
      s = cell_dim.y-1;
      ipbc.y--;
    }
#ifndef TWOD
    if (t<0) {
      t = cell_dim.z-1;
      ipbc.z--;
    }
#endif
    if (r>cell_dim.x-1) {
      r = 0; 
      ipbc.x++;
    }
    if (s>cell_dim.y-1) {
      s = 0; 
      ipbc.y++;
    }
#ifndef TWOD
    if (t>cell_dim.z-1) {
      t = 0; 
      ipbc.z++;
    }
#endif

    if (   ((pbc_dirs.x==1) || (pbc_dirs.x==ipbc.x))
#ifndef TWOD
        && ((pbc_dirs.z==1) || (pbc_dirs.z==ipbc.z)) 
#endif
        && ((pbc_dirs.y==1) || (pbc_dirs.y==ipbc.y))) 
      {

#ifdef TWOD
      q = PTR_2D_V(cell_array,r,s,cell_dim);
      pbc.x = ipbc.x * box_x.x + ipbc.y * box_y.x; 
      pbc.y = ipbc.x * box_x.y + ipbc.y * box_y.y; 
#else
      q = PTR_3D_V(cell_array,r,s,t,cell_dim);
      pbc.x = ipbc.x * box_x.x + ipbc.y * box_y.x + ipbc.z * box_z.x;
      pbc.y = ipbc.x * box_x.y + ipbc.y * box_y.y + ipbc.z * box_z.y;
      pbc.z = ipbc.x * box_x.z + ipbc.y * box_y.z + ipbc.z * box_z.z;
#endif
      tmp.x = pos.x - pbc.x;
      tmp.y = pos.y - pbc.y;
#ifndef TWOD
      tmp.z = pos.z - pbc.z;
#endif

      qptr   = q->ort;
      for (i = 0; i < q->n; i++) {
        d.x = *qptr - tmp.x; ++qptr;
        d.y = *qptr - tmp.y; ++qptr;
#ifndef TWOD
        d.z = *qptr - tmp.z; ++qptr;
#endif
        if (p_num != q->nummer[i]) {

          q_typ  = q->sorte[i];
          col    = p_typ * ntypes + q_typ;

          /* potential energy */
          radius2 = SPROD(d,d);
          if (radius2 <= pair_pot.end[col]) {
            PAIR_INT2(pot_tmp,grad_tmp,pair_pot,col,inc,radius2,r2_short,is_short)
            result += pot_tmp;
          }
	}
      }
    }
  }
  return result;
}


/* one MC sweep - this is the move_atoms routine for MC */

void one_mc_step( void )
{
  int     p_num, p_typ;
  int     i, k, l, m;
  int     total=0, accepted=0;
  ivektor cellc, newcellc;
  vektor  old_pos, new_pos;
  real    ediff;
  cell    *p, *q;

  /* for each cell */
  for (k=0; k < cell_dim.x; ++k)
  for (l=0; l < cell_dim.y; ++l)
#ifndef TWOD
  for (m=0; m < cell_dim.z; ++m)
#endif
  {
#ifndef TWOD
    p = PTR_3D_V(cell_array,k,l,m,cell_dim);
#else
    p = PTR_2D_V(cell_array,k,l,cell_dim);
#endif   

    cellc.x = k;
    cellc.y = l;
#ifndef TWOD
    cellc.z = m;
#endif   

    for (i=0; i < p->n; ++i) {
       p_typ     = p->sorte[i];
       p_num     = p->nummer[i];
       old_pos.x = p->ort X(i);
       old_pos.y = p->ort Y(i);
#ifndef TWOD
       old_pos.z = p->ort Z(i);
#endif
       new_pos.x = old_pos.x + mc_len * (real)gasdev( &mc_seed );
       new_pos.y = old_pos.y + mc_len * (real)gasdev( &mc_seed );
#ifndef TWOD
       new_pos.z = old_pos.z + mc_len * (real)gasdev( &mc_seed );
#endif

       ediff = mc_epot_diff( old_pos, new_pos, p_num, p_typ, cellc);
       /*       printf("%10.4e\n",ediff);    */
       total++;

       if ((ediff <= 0.0) || (exp(-mc_beta*ediff) >= (real)ran1(&mc_seed))) {
         accepted++;
         p->ort X(i) = new_pos.x;
         p->ort Y(i) = new_pos.y;
#ifndef TWOD
         p->ort Z(i) = new_pos.z;
#endif

#ifndef TWOD
         newcellc = cell_coord( new_pos.x, new_pos.y, new_pos.z );
#else
         newcellc = cell_coord( new_pos.x, new_pos.y );
#endif
         q = PTR_VV(cell_array,newcellc,cell_dim);
         move_atom( q, p, i );
       }
    }
  }
  mc_accept += (real)accepted/(real)total;
  mc_count++;
}


/* total potential energy per particle */

real mc_epot_part( void )
{
  int     num, typ, i, k, l, m;
  ivektor cellc;
  vektor  pos;
  real    result=0.0;
  cell    *p;

  /* for each cell */
  for (k=0; k < cell_dim.x; ++k)
  for (l=0; l < cell_dim.y; ++l)
#ifndef TWOD
  for (m=0; m < cell_dim.z; ++m)
#endif
  {
#ifndef TWOD
    p = PTR_3D_V(cell_array,k,l,m,cell_dim);
#else
    p = PTR_2D_V(cell_array,k,l,cell_dim);
#endif   

    cellc.x = k;
    cellc.y = l;
#ifndef TWOD
    cellc.z = m;
#endif   

    for (i=0; i < p->n; ++i) {
       typ   = p->sorte[i];
       num   = p->nummer[i];
       pos.x = p->ort X(i);
       pos.y = p->ort Y(i);
#ifndef TWOD
       pos.z = p->ort Z(i);
#endif
       p->pot_eng[i] = mc_epot_atom( pos, num, typ, cellc ) * 0.5;
       result       += p->pot_eng[i];
    };
  };

  result /= (real)natoms;

  return(result);

}




