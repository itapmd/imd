
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
*  flush outbuf to disk, or send it to CPU 0
*  the last buffer is sent with a different tag
*
******************************************************************************/

void flush_outbuf(FILE *out, int *len, int tag)
{
  if ((myid==0) || (parallel_output==1)) {
    if (*len>0) fwrite(outbuf, 1, *len, out);
  }
#ifdef MPI
  else /* we add 1 to the length so that even an empty message is sent */
    MPI_Send( outbuf, *len+1, MPI_CHAR, 0, tag, cpugrid );
#endif
  *len=0;
}

/******************************************************************************
*
*  write_config_select writes selected data of a configuration to a 
*  file with number fzhlr and specified suffix. The data is written
*  (and selected) by the function *write_atoms_fun, which is supposed 
*  to write the data of one cell.
*
******************************************************************************/

void write_config_select(int fzhlr, char *suffix, 
                         void (*write_atoms_fun)(FILE *out))
{
  FILE *out=NULL;
  str255 fname, msg;

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

  /* write or send own data */
  (*write_atoms_fun)(out);
#ifdef MPI
  /* if serial output, receive and write foreign data */
  if ((0==myid) && (parallel_output==0)) {
    MPI_Status status;
    int m=1, len;
    do {
      MPI_Recv(outbuf, OUTPUT_BUF_SIZE, MPI_CHAR, MPI_ANY_SOURCE, 
               MPI_ANY_TAG, cpugrid, &status);
      MPI_Get_count(&status, MPI_CHAR, &len);
      if (status.MPI_TAG==OUTBUF_TAG+1) m++;
      if (len>1) fwrite(outbuf, 1, len-1, out);
    } while (m < num_cpus);
  }
#endif /* MPI */
  if ((0==myid) || (1==parallel_output)) fclose(out);
}

/******************************************************************************
*
*  write_config writes a configuration to a numbered file,
*  which can serve as a checkpoint; uses write_atoms
*
******************************************************************************/

void write_config(int steps)
{ 
  int fzhlr = steps / rep_interval;

  /* first make sure that every atom is inside the box and on the right CPU */
  if (1==parallel_output) {
    do_boundaries();
    fix_cells();
  }

  /* write checkpoint */
  write_config_select(fzhlr, "chkpt", write_atoms);

  /* write iteration file */
  if (myid == 0) write_itr_file(fzhlr, steps);
}

#ifdef EFILTER

/******************************************************************************
*
*  filter function for write_config_select
*  writes an 'energy filtered' configuration
*
******************************************************************************/

void write_atoms_ef(FILE *out)
{
  int i, k, len=0;
  cell *p;
  double h;

  for (k=0; k<ncells; k++) {
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; i++) {
      if ( pic_ur.x != (real)0 ) /*if pic_ur.x still 0, write everything */
        if ((p->ort X(i) < pic_ll.x) && (p->ort X(i) > pic_ur.x) &&
#ifndef TWOD
            (p->ort Z(i) < pic_ll.z) && (p->ort Z(i) > pic_ur.z) && 
#endif
            (p->ort Y(i) < pic_ll.y) && (p->ort Y(i) > pic_ur.y)) continue;

      if ( (POTENG(p,i)>=lower_e_pot) && (POTENG(p,i)<=upper_e_pot) ) {
        len += sprintf( outbuf+len,
#ifdef TWOD
          "%d %d %12f %12f %12f %12f %12f %12f\n",
#else
          "%d %d %12f %12f %12f %12f %12f %12f %12f %12f\n",
#endif
          NUMMER(p,i), VSORTE(p,i), MASSE(p,i),
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
        /* flush or send outbuf if it is full */
        if (len > OUTPUT_BUF_SIZE - 256) flush_outbuf(out,&len,OUTBUF_TAG);
      }
    }
  }
  flush_outbuf(out,&len,OUTBUF_TAG+1);
}
#endif /* EFILTER */


#ifdef STRESS_TENS

/******************************************************************************
*
*  filter function for write_config_select
*  writes pressure tensor for each atom to files *.nr.press
*
******************************************************************************/

void write_atoms_press(FILE *out)
{
  int i, k, len=0;
  cell *p;

  for (k=0; k<ncells; k++) {
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) {
#ifdef TWOD
      len += sprintf( outbuf+len, 
        "%10.4e %10.4e %10.4e %10.4e %10.4e\n", 
        p->ort X(i),p->ort Y(i),
        p->presstens X(i),p->presstens Y(i),
        p->presstens_offdia[i]);
#else
      len += sprintf( outbuf+len,
        "%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n", 
        p->ort X(i),p->ort Y(i),p->ort Z(i),
        p->presstens X(i),p->presstens Y(i),p->presstens Z(i),
        p->presstens_offdia X(i),p->presstens_offdia Y(i),
        p->presstens_offdia Z(i));
#endif
      /* flush or send outbuf if it is full */
      if (len > OUTPUT_BUF_SIZE - 256) flush_outbuf(out,&len,OUTBUF_TAG);
    }
  }
  flush_outbuf(out,&len,OUTBUF_TAG+1);
}
#endif /* STRESS_TENS */


/******************************************************************************
*
*  filter function for write_config_select
*  writes .pic files (raw data for pictures)
*
******************************************************************************/

void write_atoms_pic(FILE *out) 
{
  typedef struct { 
    float   pos_x, pos_y;
#ifndef TWOD
    float   pos_z; 
#endif
    float   E_kin, E_pot;
    integer type;
  } picbuf_t;

  int i, k, len=0, sz=sizeof(picbuf_t);
  picbuf_t *picbuf;
  cell *p;

  for (k=0; k<ncells; k++) {
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) {
      picbuf = (picbuf_t *) (outbuf+len);
      picbuf->pos_x = (float) p->ort X(i);
      picbuf->pos_y = (float) p->ort Y(i);
#ifndef TWOD
      picbuf->pos_z = (float) p->ort Z(i);
#endif
      if ( pic_ur.x != (real)0 ) /*if pic_ur still 0, write everything */
      if ( (picbuf->pos_x < pic_ll.x) || (picbuf->pos_x > pic_ur.x) ||
#ifndef TWOD
           (picbuf->pos_z < pic_ll.z) || (picbuf->pos_z > pic_ur.z) ||
#endif
           (picbuf->pos_y < pic_ll.y) || (picbuf->pos_y > pic_ur.y) ) continue;

      picbuf->E_kin = (float) SPRODN(p->impuls,i,p->impuls,i)/(2*MASSE(p,i));
#ifdef DISLOC
      if (Epot_diff==1)
        picbuf->E_pot = (float) POTENG(p,i) - p->Epot_ref[i];
      else
#endif
#if defined(ORDPAR) && !defined(TWOD)
      picbuf->E_pot = (p->nbanz[i]==0) ? 0 : p->pot_eng[i]/p->nbanz[i];
#else
      picbuf->E_pot = POTENG(p,i);
#endif
      picbuf->type  = (integer) VSORTE(p,i);
      len += sz;
      /* flush or send outbuf if it is full */
      if (len > OUTPUT_BUF_SIZE - 256) flush_outbuf(out,&len,OUTBUF_TAG);
    }
  }
  flush_outbuf(out,&len,OUTBUF_TAG+1);
}


#ifdef DISLOC

/******************************************************************************
*
*  filter function for write_config_select
*  writes a differential energy map to files *.dem.x
*
******************************************************************************/

void write_atoms_dem(FILE *out)
{
  int i, k, len=0;
  cell *p;
  real dpot;

  for (k=0; k<ncells; k++) {
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) {
      if (p->sorte[i] == dpotsorte) {
        dpot = ABS(p->pot_eng[i] - p->Epot_ref[i]);
        if (dpot > min_dpot) {
#ifdef TWOD
          len += sprintf( outbuf+len, "%12f %12f %12f\n",
            p->ort X(i), p->ort Y(i), dpot);
#else
          len += sprintf( outbuf+len, "%12f %12f %12f %12f\n",
            p->ort X(i), p->ort Y(i), p->ort Z(i), dpot);
#endif
          /* flush or send outbuf if it is full */
          if (len > OUTPUT_BUF_SIZE - 256) flush_outbuf(out,&len,OUTBUF_TAG);
	}
      }
    }
  }
  flush_outbuf(out,&len,OUTBUF_TAG+1);
}

/******************************************************************************
*
*  filter function for write_config_select
*  writes a differential displacement map to files *.dsp.x
*
******************************************************************************/

void write_atoms_dsp(FILE *out)
{
  int i, k, len=0;
  cell *p;
  vektor d;

  for (k=0; k<ncells; k++) {
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) {
      d.x = p->ort X(i) - p->ort_ref X(i);
      d.y = p->ort Y(i) - p->ort_ref Y(i);
#ifndef TWOD 
      d.z = p->ort Z(i) - p->ort_ref Z(i);
#endif
      reduce_displacement(&d);
#ifdef TWOD
      len += sprintf( outbuf+len, "%12f %12f %12f %12f\n",
        p->ort X(i), p->ort Y(i), d.x, d.y);
#else
      len += sprintf( outbuf+len, "%12f %12f %12f %12f %12f %12f\n",
        p->ort X(i), p->ort Y(i), p->ort Z(i), d.x, d.y, d.z);
#endif
      /* flush or send outbuf if it is full */
      if (len > OUTPUT_BUF_SIZE - 256) flush_outbuf(out,&len,OUTBUF_TAG);
    }
  }
  flush_outbuf(out,&len,OUTBUF_TAG+1);
}

#endif

#if defined(DISLOC) || defined(AVPOS)

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

#endif /* DISLOC, AVPOS */

#ifdef AVPOS

/******************************************************************************
*
*  filter function for write_config_select
*  writes average position to files *.nr.avp
*
******************************************************************************/

void write_atoms_avp(FILE *out)
{
  int i, k, len=0;
  cell *p;
  real x, y, z;

  for (k=0; k<ncells; k++) {
    cell* p;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; i++) {
      x = p->ort_ref X(i) * avpos_res / avpos_int;
      if ( pbc_dirs.x == 1 && x > box_x.x)  x -= box_x.x;
      else if ( pbc_dirs.x == 1 && x < 0.0) x += box_x.x;
      y = p->ort_ref Y(i) * avpos_res / avpos_int;
      if ( pbc_dirs.y == 1 && y > box_y.y)  y -= box_y.y;
      else if ( pbc_dirs.y == 1 && y < 0.0) y += box_y.y;
#ifndef TWOD
      z = p->ort_ref Z(i) * avpos_res / avpos_int;
      if ( pbc_dirs.z == 1 && z > box_z.z)  z -= box_z.z;
      else if ( pbc_dirs.z == 1 && z < 0.0) z += box_z.z;

      len += sprintf( outbuf+len,
        "%d %d %12.16f %12.16f %12.16f %12.16f %12.16f\n ",
        NUMMER(p,i), VSORTE(p,i), MASSE(p,i), 
        x, y, z, p->Epot_ref[i] * avpos_res / avpos_int);
#else
      len += sprintf( outbuf+len,
        "%d %d %12.16f %12.16f %12.16f %12.16f\n ",
        NUMMER(p,i), VSORTE(p,i), MASSE(p,i), 
        x, y, p->Epot_ref[i] * avpos_res / avpos_int);
#endif
      /* flush or send outbuf if it is full */
      if (len > OUTPUT_BUF_SIZE - 256) flush_outbuf(out,&len,OUTBUF_TAG);
    }
  }
  flush_outbuf(out,&len,OUTBUF_TAG+1);
}

#endif

/******************************************************************************
*
*  write header of properties file  - keep in sync with write_properties
*
******************************************************************************/

void write_eng_file_header()
{
  str255 fname;
  FILE *fl;

  if (myid == 0) {

    sprintf(fname,"%s.eng",outfilename);
    fl = fopen(fname,"w");
    if (NULL == fl) error("Cannot open properties file.");

    fprintf(fl, "# time Epot ");
    fprintf(fl, "temperature ");
#ifdef FNORM
    fprintf(fl, "fnorm ");
#endif
#ifdef GLOK
    fprintf(fl, "PxF ");
#endif
    fprintf(fl, "pressure ");
    fprintf(fl, "volume ");
#if defined(NVT) || defined(NPT)
    fprintf(fl, "eta ");
#endif
#ifdef NPT_axial
#ifdef TWOD
    fprintf(fl, "stress_x stress_y ");
    fprintf(fl, "box_x.x box_y.y ");
#else
    fprintf(fl, "stress_x stress_y stress_z ");
    fprintf(fl, "box_x.x box_y.y box_z.z ");
#endif
#endif
    putc('\n',fl);

    fclose(fl);
  }
}

/******************************************************************************
*
*  write selected properties to *.eng file
*  keep in sync with write_eng_file_header
*
******************************************************************************/

void write_properties(int steps)
{
  FILE *out;
  str255 fname;
#ifdef HPO
  char *format=" %.16e";
#else
  char *format=" %e";
#endif
  real Epot, Temp, vol;

  Epot =       tot_pot_energy / natoms;
#ifdef UNIAX
  Temp = 2.0 * tot_kin_energy / (nactive + nactive_rot);
#else
  Temp = 2.0 * tot_kin_energy / nactive;
#endif
  vol  = volume / natoms;

  sprintf(fname,"%s.eng",outfilename);
  out = fopen(fname,"a");
  if (NULL == out) error("Cannot open properties file.");

  fprintf(out, "%e",     (double) (steps * timestep));
  fprintf(out, " %.16e ", (double) Epot);
 
  fprintf(out, format,   (double) Temp);
#ifdef FNORM
  fprintf(out, format,   (double) fnorm);
#endif
#ifdef GLOK
  fprintf(out, format,   (double) PxF);
#endif
  fprintf(out," %e",     (double) pressure);
  fprintf(out," %e",     (double) vol);
#if defined(NVT) || defined(NPT)
  fprintf(out," %e",     (double) eta );
#endif
  if (ensemble==ENS_NPT_AXIAL) {
#ifdef TWOD
    fprintf(out," %e %e", (double) stress_x, (double) stress_y );
    fprintf(out," %e %e", (double)  box_x.x, (double)  box_y.y );
#else
    fprintf(out," %e %e %e", 
                  (double) stress_x, (double) stress_y, (double) stress_z );
    fprintf(out," %e %e %e", 
                  (double)  box_x.x, (double)  box_y.y, (double)  box_z.z );
#endif
  }
  putc('\n',out);
  fclose(out);
}

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
