
/******************************************************************************
*
* imd_io.c -- dimension independent IO  routines 
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
  int k,m,tag;

  if (0==myid) {

    /* open output file */
    sprintf(fname,"%s.%s.%u",outfilename,suffix,fzhlr);
    out = fopen(fname,"w");
    if (NULL == out) {
       sprintf(msg,"Cannot open output file %s",fname);
       error(msg);
    }

    /* write own data */
    for (k=0; k<ncells; k++) {
      p = cell_array + CELLS(k);
      (*write_cell_fun)(out,p);
    }

#ifdef MPI
    p = cell_array;  /* this is a pointer to the first (buffer) cell */
    for (m=1; m<num_cpus; ++m)
      for (k=0; k<ncells; k++) {
        tag = CELL_TAG + CELLS(k);
        recv_cell(p,m,tag);
        (*write_cell_fun)(out,p);
      }
#endif

    fclose(out);      
  } 

#ifdef MPI
  else { 
    /* send data to cpu 0 */
    for (k=0; k<ncells; k++) {
      p = cell_array + CELLS(k);
      tag = CELL_TAG + CELLS(k);
      send_cell(p,0,tag);
    }
  }
#endif

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
          SORTE(p,i),
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
*  writes pressure tensor for each atom to files *.press.x
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
  real dx,dy,dz;
  for (i=0; i<p->n; ++i) {
    dx = p->ort X(i) - p->ort_ref X(i);
    dy = p->ort Y(i) - p->ort_ref Y(i);
#ifndef TWOD 
    dz = p->ort Z(i) - p->ort_ref Z(i);
    if ((ABS(dx)<box_x.x/2) && (ABS(dy)<box_y.y/2) && (ABS(dz)<box_z.z/2))
    fprintf(out,"%12f %12f %12f %12f %12f %12f\n",
            p->ort X(i),p->ort Y(i),p->ort Z(i),dx,dy,dz);
#else
    if ((ABS(dx)<box_x.x/2) && (ABS(dy)<box_y.y/2))
    fprintf(out,"%12f %12f %12f %12f\n",
            p->ort X(i),p->ort Y(i),dx,dy);
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


