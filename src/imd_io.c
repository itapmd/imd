
/******************************************************************************
*
* imd_io.c -- dimension independent IO routines 
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include "imd.h"

/******************************************************************************
*
*  write_config_select writes selected data of a configuration to a 
*  file with number fzhlr and specified suffix. The data is written
*  (and selected) by the function *write_cell_fun, which is supposed 
*  to write the data of one cell.
*
******************************************************************************/

void write_config_select(int fzhlr, char *suffix, 
                         void (*write_cell_fun)(FILE *out, cell *p))
{
  FILE *out;
  str255 fname,msg;
  cell *p;
  int k,m;

  /* create file name, and open file */
#ifdef MPI
  if (1==parallel_output) {
    sprintf(fname,"%s.%u.%s.%u",outfilename,fzhlr,suffix,myid);
    out = fopen(fname,"w");
    if (NULL == out) { 
       sprintf(msg,"Cannot open output file %s",fname);
       error(msg);
    }
  } else
#endif
  if (0==myid) {
    sprintf(fname,"%s.%u.%s",outfilename,fzhlr,suffix);
    out = fopen(fname,"w");
    if (NULL == out) {
       sprintf(msg,"Cannot open output file %s",fname);
       error(msg);
    }
  }

  /* write data */
#ifndef MPI
  for (k=0; k<ncells; k++) {
    p = cell_array + CELLS(k);
    (*write_cell_fun)(out,p);
  }
  fclose(out);
#else
  if (1==parallel_output) {
    for (k=0; k<ncells; k++) {
      p = cell_array + CELLS(k);
      (*write_cell_fun)(out,p);
    }
    fclose(out);
  } else if (0==myid) {
    /* write own data */
    for (k=0; k<ncells; k++) {
      p = cell_array + CELLS(k);
      (*write_cell_fun)(out,p);
    }
    /* write foreign data */
    p = cell_array;  /* this is a pointer to the first (buffer) cell */
    for (m=1; m<num_cpus; ++m)
      for (k=0; k<ncells; k++) {
        recv_cell(p,MPI_ANY_SOURCE,CELL_TAG); /* accept cells in any order */
        (*write_cell_fun)(out,p);
      }
    fclose(out);      
  } else { 
    /* send data to cpu 0 */
    for (k=0; k<ncells; k++) {
      p = cell_array + CELLS(k);
      send_cell(p,0,CELL_TAG);
    }
  }
#endif /* MPI */

}

#ifdef EFILTER

/******************************************************************************
*
*  filter function for write_config_select
*  writes an 'energy filtered' configuration
*
******************************************************************************/

void write_cell_ef(FILE *out, cell *p)
{
  int i;
  double h;

  for (i=0; i<p->n; i++)
#ifdef ZOOM
    if ((p->ort X(i) >= pic_ll.x) && (p->ort X(i) <= pic_ur.x) &&
#ifndef TWOD
        (p->ort Z(i) >= pic_ll.z) && (p->ort Z(i) <= pic_ur.z) && 
#endif
        (p->ort Y(i) >= pic_ll.y) && (p->ort Y(i) <= pic_ur.y))
#endif
      if ( (POTENG(p,i)>=lower_e_pot) && (POTENG(p,i)<=upper_e_pot) ) {
#ifdef TWOD
        fprintf(out,"%d %d %12f %12f %12f %12f %12f %12f\n",
#else
        fprintf(out,"%d %d %12f %12f %12f %12f %12f %12f %12f %12f\n",
#endif
          NUMMER(p,i),
          p->sorte[i],
          MASSE(p,i),
          p->ort X(i),
          p->ort Y(i),
#ifndef TWOD
          p->ort Z(i),
#endif
          p->impuls X(i) / MASSE(p,i),
          p->impuls Y(i) / MASSE(p,i),
#ifndef TWOD
          p->impuls Z(i) / MASSE(p,i),
#endif
          POTENG(p,i)
        );
      }
}
#endif /* EFILTER */


#ifdef STRESS_TENS

/******************************************************************************
*
*  filter function for write_config_select
*  writes pressure tensor for each atom to files *.nr.press
*
******************************************************************************/

void write_cell_press(FILE *out, cell *p)
{
  int i;
  for (i=0; i<p->n; ++i) {
#ifdef TWOD
    fprintf(out,"%10.4e %10.4e %10.4e %10.4e %10.4e\n", 
      p->ort X(i),p->ort Y(i),
      p->presstens X(i),p->presstens Y(i),
      p->presstens_offdia[i]);
#else
    fprintf(out,
      "%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n", 
      p->ort X(i),p->ort Y(i),p->ort Z(i),
      p->presstens X(i),p->presstens Y(i),p->presstens Z(i),
      p->presstens_offdia X(i),p->presstens_offdia Y(i),
      p->presstens_offdia Z(i));
#endif
  }
}
#endif /* STRESS_TENS */


/******************************************************************************
*
*  filter function for write_config_select
*  writes .pic files (raw data for pictures)
*
******************************************************************************/

void write_cell_pic(FILE *out, cell *p) 
{
  struct { 
    float   pos_x, pos_y;
#ifndef TWOD
    float   pos_z; 
#endif
    float   E_kin, E_pot;
    integer type;
  } picbuf;

  int i;

  for (i=0; i<p->n; ++i) {
    picbuf.pos_x = (float) p->ort X(i);
    picbuf.pos_y = (float) p->ort Y(i);
#ifndef TWOD
    picbuf.pos_z = (float) p->ort Z(i);
#endif
    if ( pic_ur.x != (real)0 ) /*if pic_ur still 0, write everything */
    if ( (picbuf.pos_x < pic_ll.x) || (picbuf.pos_x > pic_ur.x) ||
#ifndef TWOD
         (picbuf.pos_z < pic_ll.z) || (picbuf.pos_z > pic_ur.z) ||
#endif
         (picbuf.pos_y < pic_ll.y) || (picbuf.pos_y > pic_ur.y) ) continue;

    picbuf.E_kin = (float) SPRODN(p->impuls,i,p->impuls,i)/(2*MASSE(p,i));
#ifdef DISLOC
    if (Epot_diff==1)
      picbuf.E_pot = (float) POTENG(p,i) - p->Epot_ref[i];
    else
#endif
#if defined(ORDPAR) && !defined(TWOD)
    picbuf.E_pot = (p->nbanz[i]==0) ? 0 : p->pot_eng[i]/p->nbanz[i];
#else
    picbuf.E_pot = POTENG(p,i);
#endif
    picbuf.type  = (integer) VSORTE(p,i);
    fwrite(&picbuf, sizeof(picbuf), 1, out); 
  }
}


#ifdef DISLOC

/******************************************************************************
*
*  filter function for write_config_select
*  writes a differential energy map to files *.dem.x
*
******************************************************************************/

void write_cell_dem(FILE *out, cell *p)
{
  int i;
  real dpot;
  for (i=0; i<p->n; ++i) {
    if (p->sorte[i] == dpotsorte) {
      dpot = ABS(p->pot_eng[i] - p->Epot_ref[i]);
      if (dpot > min_dpot)
#ifdef TWOD
        fprintf(out,"%12f %12f %12f\n",
                p->ort X(i),p->ort Y(i),dpot);
#else
        fprintf(out,"%12f %12f %12f %12f\n",
                p->ort X(i),p->ort Y(i),p->ort Z(i),dpot);
#endif
    }
  }
}

/******************************************************************************
*
*  filter function for write_config_select
*  writes a differential displacement map to files *.dsp.x
*
******************************************************************************/

void write_cell_dsp(FILE *out, cell *p)
{
  int i;
  vektor d;
  for (i=0; i<p->n; ++i) {
    d.x = p->ort X(i) - p->ort_ref X(i);
    d.y = p->ort Y(i) - p->ort_ref Y(i);
#ifndef TWOD 
    d.z = p->ort Z(i) - p->ort_ref Z(i);
#endif
    reduce_displacement(&d);
#ifdef TWOD
    fprintf(out,"%12f %12f %12f %12f\n",
            p->ort X(i),p->ort Y(i),d.x,d.y);
#else
    fprintf(out,"%12f %12f %12f %12f %12f %12f\n",
            p->ort X(i),p->ort Y(i),p->ort Z(i),d.x,d.y,d.z);
#endif
  }
}


/******************************************************************************
*
*  reset_Epot_ref
*
******************************************************************************/

void reset_Epot_ref(void)
{
  int  k;
  for (k=0; k<ncells; ++k) {
    int  i;
    cell *p;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) {
      p->Epot_ref[i] = p->pot_eng[i];
    }
  }
}


/******************************************************************************
*
* update_ort_ref updates ort_ref
*
******************************************************************************/

void update_ort_ref(void)
{ 
  int k;
  for (k=0; k<ncells; k++) {
    int i;
    cell* p;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; i++) {
      p->ort_ref X(i) = p->ort X(i);
      p->ort_ref Y(i) = p->ort Y(i);
#ifndef TWOD
      p->ort_ref Z(i) = p->ort Z(i);
#endif
    }
  }
}

#endif /* DISLOC */


#ifdef MSQD

/******************************************************************************
*
* write_msqd writes mean square displacement to *.msqd file
*
******************************************************************************/

void write_msqd(int steps)
{
  FILE *out;
  str255 fname;
  int i;

  sprintf(fname,"%s.msqd",outfilename);
  out = fopen(fname,"a");
  if (NULL == out) error("Cannot open msqd file.");

  fprintf(out, "%10.4e", (double)(steps * timestep));
  for (i=0; i<ntypes; i++) {
    fprintf(out," %10.4e", (double)(msqd_global[i] / num_sort[i]));
  }
  putc('\n',out);
  fclose(out);
}

#endif

/******************************************************************************
*
* reduce displacement vector modulo box vectors, so that it is in a
* box centered at the origin
*
******************************************************************************/

void reduce_displacement(vektor *dist)
{
  vektor d;
  d = *dist;
  while (SPROD(d,tbox_x) > 0.5) {
    d.x -= box_x.x;
    d.y -= box_x.y;
#ifndef TWOD
    d.z -= box_x.z;
#endif
  }
  while (SPROD(d,tbox_x) < -0.5) {
    d.x += box_x.x;
    d.y += box_x.y;
#ifndef TWOD
    d.z += box_x.z;
#endif
  }
  while (SPROD(d,tbox_y) > 0.5) {
    d.x -= box_y.x;
    d.y -= box_y.y;
#ifndef TWOD
    d.z -= box_y.z;
#endif
  }
  while (SPROD(d,tbox_y) < -0.5) {
    d.x += box_y.x;
    d.y += box_y.y;
#ifndef TWOD
    d.z += box_y.z;
#endif
  }
#ifndef TWOD
  while (SPROD(d,tbox_z) > 0.5) {
    d.x -= box_z.x;
    d.y -= box_z.y;
    d.z -= box_z.z;
  }
  while (SPROD(d,tbox_z) < -0.5) {
    d.x += box_z.x;
    d.y += box_z.y;
    d.z += box_z.z;
  }
#endif
  *dist = d;
}
