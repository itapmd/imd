
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

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
  if (*len+1>OUTPUT_BUF_SIZE) error("outbuf overflow");
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
  void (*write_atoms_fun)(FILE *out), void (*write_header_fun)(FILE *out))
{
  FILE *out=NULL;
  str255 fname, msg;

#ifdef MPI
  if (1==parallel_output) {
    /* write header to separate file */
    if ((myid==0) && (use_header)) {
      if (fzhlr>=0) sprintf(fname,"%s.%05d.%s.head",outfilename,fzhlr,suffix);
      else          sprintf(fname,"%s-final.%s.head",outfilename,suffix);
      out = fopen(fname, "w");
      if (NULL == out) { 
        sprintf(msg,"Cannot open output file %s",fname);
	error(msg);
      }
      (*write_header_fun)(out);
      fclose(out);
    }
    /* open output file */
    if (fzhlr>=0) sprintf(fname,"%s.%05d.%s.%u",outfilename,fzhlr,suffix,myid);
    else          sprintf(fname,"%s-final.%s.%u", outfilename,suffix,myid);
    out = fopen(fname,"w");
    if (NULL == out) { 
       sprintf(msg,"Cannot open output file %s",fname);
       error(msg);
    }
  } else
#endif
  if (0==myid) {
    /* open output file */
    if (fzhlr >= 0) sprintf(fname,"%s.%05d.%s", outfilename,fzhlr,suffix);
    else            sprintf(fname,"%s-final.%s", outfilename,suffix);
    out = fopen(fname,"w");
    if (NULL == out) {
       sprintf(msg,"Cannot open output file %s",fname);
       error(msg);
    }
    /* write header */
    if (use_header)
      (*write_header_fun)(out);
  }

  /* write or send own data */
  (*write_atoms_fun)(out);
#ifdef MPI
  /* if serial output, receive and write foreign data */
  if ((0==myid) && (parallel_output==0)) {
    MPI_Status status;
    int m=1, len;
    while (m < num_cpus) {
      MPI_Recv(outbuf, OUTPUT_BUF_SIZE, MPI_CHAR, MPI_ANY_SOURCE, 
               MPI_ANY_TAG, cpugrid, &status);
      MPI_Get_count(&status, MPI_CHAR, &len);
      if ((status.MPI_TAG!=OUTBUF_TAG+1) && (status.MPI_TAG!=OUTBUF_TAG))
        error("messages mixed up");
      if (status.MPI_TAG==OUTBUF_TAG+1) m++;
      if (len>1) fwrite(outbuf, 1, len-1, out);
    }
  }
  /* don't send non-io messages before we are finished */
  MPI_Barrier(cpugrid);
#endif /* MPI */
  if ((0==myid) || (1==parallel_output)) fclose(out);
}

/******************************************************************************
*
*  write_config writes a configuration to a numbered file,
*  which can serve as a checkpoint; uses write_atoms
*
******************************************************************************/

void write_config(int fzhlr, int steps)
{ 
  /* first make sure that every atom is inside the box and on the right CPU */
  if (1==parallel_output) {
    do_boundaries();
    fix_cells();
  }

  /* write checkpoint */
  write_config_select(fzhlr, "chkpt", write_atoms_config, write_header_config);

  /* write iteration file */
  if (myid == 0) write_itr_file(fzhlr, steps,"");
}

#ifdef SNAPSHOT
/******************************************************************************
*
*  write_ssconfig writes a configuration to a numbered file,
*  which can serve as a checkpoint; uses write_atoms
*
******************************************************************************/

void write_ssconfig(int steps)
{ 
  
  /* first make sure that every atom is inside the box and on the right CPU */
  if (1==parallel_output) {
    do_boundaries();
    fix_cells();
  }

  /* write checkpoint */
  write_config_select(sscount, "ss", write_atoms_config, write_header_config);

  /* write iteration file */
  if (myid == 0) write_itr_file(sscount, steps,"ss");
}
#endif

#ifdef CG
/******************************************************************************
*
*  write_cgconfig writes a configuration to a numbered file,
*  which can serve as a checkpoint; uses write_atoms
* 
******************************************************************************/

void write_cgconfig(int steps)
{ 
  
  /* first make sure that every atom is inside the box and on the right CPU */
  if (1==parallel_output) {
    do_boundaries();
    fix_cells();
  }

  /* write checkpoint */
  write_config_select(steps, "cgchkpt", write_atoms_config, write_header_config);

  /* write iteration file */
  if (myid == 0) write_itr_file(steps, steps,"cg");
}
#endif


#ifdef EFILTER

/******************************************************************************
*
*  writes header for ef-file
*
******************************************************************************/

void write_header_ef(FILE *out)
{
  char c;
  time_t now;

  /* format line */
  /* binary_io for checkpoints not implemented
   if (binary_io)
     c = is_big_endian ? 'B' : 'L';
   else
    */
    c = 'A';
  fprintf(out, "#F %c 1 1 1 %d %d 1\n", c, DIM, DIM);
  
  /* contents line */
#ifdef TWOD
  fprintf(out, "#C number type mass x y vx vy Epot\n");
#else
  fprintf(out, "#C number type mass x y z vx vy vz Epot\n");
#endif

  /* box lines */
#ifdef TWOD
  fprintf(out, "#X \t%.16e %.16e\n", box_x.x , box_x.y);
  fprintf(out, "#Y \t%.16e %.16e\n", box_y.x , box_y.y);
#else
  fprintf(out, "#X \t%.16e %.16e %.16e\n", box_x.x , box_x.y , box_x.z);
  fprintf(out, "#Y \t%.16e %.16e %.16e\n", box_y.x , box_y.y , box_y.z);
  fprintf(out, "#Z \t%.16e %.16e %.16e\n", box_z.x , box_z.y , box_z.z);
#endif

  /* endheader line */
  time(&now);
  fprintf(out, "## Generated by %s on %s", progname, ctime(&now));
  fprintf(out, "#E\n");

}

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
      if ( pic_ur.x != (real)0 ) /* if pic_ur.x still 0, write everything */
        if ((ORT(p,i,X) < pic_ll.x) || (ORT(p,i,X) > pic_ur.x) ||
#ifndef TWOD
            (ORT(p,i,Z) < pic_ll.z) || (ORT(p,i,Z) > pic_ur.z) || 
#endif
            (ORT(p,i,Y) < pic_ll.y) || (ORT(p,i,Y) > pic_ur.y)) continue;

      if ( (POTENG(p,i)>=lower_e_pot) && (POTENG(p,i)<=upper_e_pot) ) {
        len += sprintf( outbuf+len,
#ifdef TWOD
          "%d %d %12f %12f %12f %12f %12f %12f\n",
#else
          "%d %d %12f %12f %12f %12f %12f %12f %12f %12f\n",
#endif
          NUMMER(p,i), VSORTE(p,i), MASSE(p,i),
          ORT(p,i,X),
          ORT(p,i,Y),
#ifndef TWOD
          ORT(p,i,Z),
#endif
          IMPULS(p,i,X) / MASSE(p,i),
          IMPULS(p,i,Y) / MASSE(p,i),
#ifndef TWOD
          IMPULS(p,i,Z) / MASSE(p,i),
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


#ifdef NBFILTER

/******************************************************************************
*
*  writes header for ef-file
*
******************************************************************************/

void write_header_nb(FILE *out)
{
  char c;
  time_t now;

  /* format line */
  /* binary_io for checkpoints not implemented
   if (binary_io)
     c = is_big_endian ? 'B' : 'L';
   else
    */
    c = 'A';
  fprintf(out, "#F %c 1 1 1 %d %d 1\n", c, DIM, DIM);
  
  /* contents line */
#ifdef TWOD
  fprintf(out, "#C number type mass x y vx vy Epot\n");
#else
  fprintf(out, "#C number type mass x y z vx vy vz Epot\n");
#endif

  /* box lines */
#ifdef TWOD
  fprintf(out, "#X \t%.16e %.16e\n", box_x.x , box_x.y);
  fprintf(out, "#Y \t%.16e %.16e\n", box_y.x , box_y.y);
#else
  fprintf(out, "#X \t%.16e %.16e %.16e\n", box_x.x , box_x.y , box_x.z);
  fprintf(out, "#Y \t%.16e %.16e %.16e\n", box_y.x , box_y.y , box_y.z);
  fprintf(out, "#Z \t%.16e %.16e %.16e\n", box_z.x , box_z.y , box_z.z);
#endif

  /* endheader line */
  time(&now);
  fprintf(out, "## Generated by %s on %s", progname, ctime(&now));
  fprintf(out, "#E\n");

}

/******************************************************************************
*
*  filter function for write_config_select
*  writes an 'neighbours filtered' configuration
*
******************************************************************************/

void write_atoms_nb(FILE *out)
{
  int i, k, len=0;
  cell *p;
  double h;

  for (k=0; k<ncells; k++) {
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; i++) {
      if ( pic_ur.x != (real)0 ) /* if pic_ur.x still 0, write everything */
        if ((ORT(p,i,X) < pic_ll.x) || (ORT(p,i,X) > pic_ur.x) ||
#ifndef TWOD
            (ORT(p,i,Z) < pic_ll.z) || (ORT(p,i,Z) > pic_ur.z) || 
#endif
            (ORT(p,i,Y) < pic_ll.y) || (ORT(p,i,Y) > pic_ur.y)) continue;

      if ( (SORTE(p,i) == VSORTE(p,i)) && 
	   ((NBANZ(p,i)< lower_nb_cut[SORTE(p,i)]) || 
	    (NBANZ(p,i)> upper_nb_cut[SORTE(p,i)]))){
	  len += sprintf( outbuf+len,
#ifdef TWOD
          "%d %d %12f %12f %12f %12f %12f %12f\n",
#else
          "%d %d %12f %12f %12f %12f %12f %12f %12f %12f\n",
#endif
          NUMMER(p,i), VSORTE(p,i), MASSE(p,i),
          ORT(p,i,X),
          ORT(p,i,Y),
#ifndef TWOD
          ORT(p,i,Z),
#endif
          IMPULS(p,i,X) / MASSE(p,i),
          IMPULS(p,i,X) / MASSE(p,i),
#ifndef TWOD
          IMPULS(p,i,Z) / MASSE(p,i),
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
#endif /* NBFILTER */

#ifdef WRITEF

/******************************************************************************
*
*  writes header for wf-file
*
******************************************************************************/

void write_header_wf(FILE *out)
{
  char c;
  time_t now;

  /* format line */
  /* binary_io for checkpoints not implemented
   if (binary_io)
     c = is_big_endian ? 'B' : 'L';
   else
    */
    c = 'A';
  fprintf(out, "#F %c 1 1 1 %d %d 1\n", c, DIM, DIM);
  
  /* contents line */
#ifdef TWOD
  fprintf(out, "#C number type mass x y fx fy Epot\n");
#else
  fprintf(out, "#C number type mass x y z fx fy fz Epot\n");
#endif

  /* box lines */
#ifdef TWOD
  fprintf(out, "#X \t%.16e %.16e\n", box_x.x , box_x.y);
  fprintf(out, "#Y \t%.16e %.16e\n", box_y.x , box_y.y);
#else
  fprintf(out, "#X \t%.16e %.16e %.16e\n", box_x.x , box_x.y , box_x.z);
  fprintf(out, "#Y \t%.16e %.16e %.16e\n", box_y.x , box_y.y , box_y.z);
  fprintf(out, "#Z \t%.16e %.16e %.16e\n", box_z.x , box_z.y , box_z.z);
#endif

  /* endheader line */
  time(&now);
  fprintf(out, "## Generated by %s on %s", progname, ctime(&now));
  fprintf(out, "#E\n");

}

/******************************************************************************
*
*  filter function for write_config_select
*  writes an 'energy filtered' configuration
*
******************************************************************************/

void write_atoms_wf(FILE *out)
{
  int i, k, len=0;
  cell *p;
  double h;

  for (k=0; k<ncells; k++) {
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; i++) {
      if(SORTE(p,i) != VSORTE(p,i)){
        len += sprintf( outbuf+len,
#ifdef TWOD
          "%d %d %12f %12f %12f %12g %12g %12f\n",
#else
          "%d %d %12f %12f %12f %12f %12e %12e %12e %12f\n",
#endif
          NUMMER(p,i), VSORTE(p,i), MASSE(p,i),
          ORT(p,i,X),
          ORT(p,i,Y),
#ifndef TWOD
          ORT(p,i,Z),
#endif
          KRAFT(p,i,X),
          KRAFT(p,i,Y),
#ifndef TWOD
          KRAFT(p,i,Z),
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
#endif /* WRITEF */


#ifdef STRESS_TENS

/******************************************************************************
*
*  writes header for press-file
*
******************************************************************************/

void write_header_press(FILE *out)
{
  char c;
  time_t now;

  /* format line */
  /* binary_io for checkpoints not implemented
   if (binary_io)
     c = is_big_endian ? 'B' : 'L';
   else
    */
    c = 'A';
  
  /* format and contents lines */
#ifdef TWOD
  fprintf(out, "#F %c 0 0 0 2 0 3\n", c);
  fprintf(out, "#C x y P_xx P_yy P_xy\n");
#else
  fprintf(out, "#F %c 0 0 0 3 0 6\n", c);
  fprintf(out, "#C x y z P_xx P_yy P_zz P_yz P_zx P_xy\n");
#endif

  /* box lines */
#ifdef TWOD
  fprintf(out, "#X \t%.16e %.16e\n", box_x.x , box_x.y);
  fprintf(out, "#Y \t%.16e %.16e\n", box_y.x , box_y.y);
#else
  fprintf(out, "#X \t%.16e %.16e %.16e\n", box_x.x , box_x.y , box_x.z);
  fprintf(out, "#Y \t%.16e %.16e %.16e\n", box_y.x , box_y.y , box_y.z);
  fprintf(out, "#Z \t%.16e %.16e %.16e\n", box_z.x , box_z.y , box_z.z);
#endif

  /* endheader line */
  time(&now);
  fprintf(out, "## Generated by %s on %s", progname, ctime(&now));
  fprintf(out, "#E\n");

}

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
        ORT(p,i,X),ORT(p,i,Y),
        PRESSTENS(p,i,xx), PRESSTENS(p,i,yy), PRESSTENS(p,i,xy) );
#else
      len += sprintf( outbuf+len,
        "%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n", 
        ORT(p,i,X),ORT(p,i,Y),ORT(p,i,Z),
        PRESSTENS(p,i,xx), PRESSTENS(p,i,yy), PRESSTENS(p,i,zz),
        PRESSTENS(p,i,yz), PRESSTENS(p,i,zx), PRESSTENS(p,i,xy) );
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
*  writes header for pic-file
*
******************************************************************************/

void write_header_pic(FILE *out)
{
  char c;
  time_t now;

  /* format line; format is always binary */
  c = is_big_endian ? 'B' : 'L';
  fprintf(out, "#F %c 0 0 0 %d 0 3\n", c, DIM);

  /* contents line */
#ifdef TWOD
  fprintf(out, "#C x y Ekin");
#else
  fprintf(out, "#C x y z Ekin");
#endif
#ifdef DISLOC
  if (Epot_diff==1)
    fprintf(out, " Epot-Epot_ref type\n");
  else
#endif
#ifdef ORDPAR 
    fprintf(out, " ordpar type\n");
#else
    fprintf(out, " Epot type\n");
#endif

  /* box lines */
#ifdef TWOD
  fprintf(out, "#X \t%.16e %.16e\n", box_x.x , box_x.y);
  fprintf(out, "#Y \t%.16e %.16e\n", box_y.x , box_y.y);
#else
  fprintf(out, "#X \t%.16e %.16e %.16e\n", box_x.x , box_x.y , box_x.z);
  fprintf(out, "#Y \t%.16e %.16e %.16e\n", box_y.x , box_y.y , box_y.z);
  fprintf(out, "#Z \t%.16e %.16e %.16e\n", box_z.x , box_z.y , box_z.z);
#endif

  /* endheader line */
  time(&now);
  fprintf(out, "## Generated by %s on %s", progname, ctime(&now));
  fprintf(out, "#E\n");

}

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
    integer end;   /* this must be the last component of picbuf_t! */
  } picbuf_t;

  int i, k, len=0, sz;
  picbuf_t *picbuf, dummy_buf;
  cell *p;

  /* get the real size of picbuf_t in bytes */
  sz = ((char *) &dummy_buf.end) - ((char *) &dummy_buf);

  for (k=0; k<ncells; k++) {
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; ++i) {
      picbuf = (picbuf_t *) (outbuf+len);
      picbuf->pos_x = (float) ORT(p,i,X);
      picbuf->pos_y = (float) ORT(p,i,Y);
#ifndef TWOD
      picbuf->pos_z = (float) ORT(p,i,Z);
#endif
      if ( pic_ur.x != (real)0 ) /* if pic_ur still 0, write everything */
      if ( (picbuf->pos_x < pic_ll.x) || (picbuf->pos_x > pic_ur.x) ||
#ifndef TWOD
           (picbuf->pos_z < pic_ll.z) || (picbuf->pos_z > pic_ur.z) ||
#endif
           (picbuf->pos_y < pic_ll.y) || (picbuf->pos_y > pic_ur.y) ) continue;

      picbuf->E_kin = (float) SPRODN( &IMPULS(p,i,X), &IMPULS(p,i,X) ) / 
                                                            (2 * MASSE(p,i));
#ifdef DISLOC
      if (Epot_diff==1)
        picbuf->E_pot = (float) (POTENG(p,i) - EPOT_REF(p,i));
      else
#endif
      picbuf->E_pot = POTENG(p,i);
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
*  write header for dem-file
*
******************************************************************************/

void write_header_dem(FILE *out)
{
  char c;
  time_t now;

  /* format line */
  /* binary_io for checkpoints not implemented
   if (binary_io)
     c = is_big_endian ? 'B' : 'L';
   else
    */
    c = 'A';
  fprintf(out, "#F %c 0 0 0 %d 0 1\n", c, DIM);
  
  /* contents line */
#ifdef TWOD
  fprintf(out, "#C x y dpot\n");
#else
  fprintf(out, "#C x y z dpot\n");
#endif

  /* box lines */
#ifdef TWOD
  fprintf(out, "#X \t%.16e %.16e\n", box_x.x , box_x.y);
  fprintf(out, "#Y \t%.16e %.16e\n", box_y.x , box_y.y);
#else
  fprintf(out, "#X \t%.16e %.16e %.16e\n", box_x.x , box_x.y , box_x.z);
  fprintf(out, "#Y \t%.16e %.16e %.16e\n", box_y.x , box_y.y , box_y.z);
  fprintf(out, "#Z \t%.16e %.16e %.16e\n", box_z.x , box_z.y , box_z.z);
#endif

  /* endheader line */
  time(&now);
  fprintf(out, "## Generated by %s on %s", progname, ctime(&now));
  fprintf(out, "#E\n");

}

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
      if (VSORTE(p,i) == dpotsorte) {
        dpot = ABS(POTENG(p,i) - EPOT_REF(p,i));
        if (dpot > min_dpot) {
#ifdef TWOD
          len += sprintf( outbuf+len, "%12f %12f %12f\n",
            ORT(p,i,X), ORT(p,i,Y), dpot);
#else
          len += sprintf( outbuf+len, "%12f %12f %12f %12f\n",
            ORT(p,i,X), ORT(p,i,Y), ORT(p,i,Z), dpot);
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
*  writes header for dsp-file
*
******************************************************************************/

void write_header_dsp(FILE *out)
{
  char c;
  time_t now;

  /* format line */
  /* binary_io for checkpoints not implemented
   if (binary_io)
     c = is_big_endian ? 'B' : 'L';
   else
    */
    c = 'A';
  fprintf(out, "#F %c 0 0 0 %d 0 %d\n", c, DIM, DIM);
  
  /* contents line */
#ifdef TWOD
  fprintf(out, "#C number type mass x y dx dy\n");
#else
  fprintf(out, "#C number type mass x y z dx dy dz\n");
#endif

  /* box lines */
#ifdef TWOD
  fprintf(out, "#X \t%.16e %.16e\n", box_x.x , box_x.y);
  fprintf(out, "#Y \t%.16e %.16e\n", box_y.x , box_y.y);
#else
  fprintf(out, "#X \t%.16e %.16e %.16e\n", box_x.x , box_x.y , box_x.z);
  fprintf(out, "#Y \t%.16e %.16e %.16e\n", box_y.x , box_y.y , box_y.z);
  fprintf(out, "#Z \t%.16e %.16e %.16e\n", box_z.x , box_z.y , box_z.z);
#endif

  /* endheader line */
  time(&now);
  fprintf(out, "## Generated by %s on %s", progname, ctime(&now));
  fprintf(out, "#E\n");

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
      d.x = ORT(p,i,X) - ORT_REF(p,i,X);
      d.y = ORT(p,i,Y) - ORT_REF(p,i,Y);
#ifndef TWOD 
      d.z = ORT(p,i,Z) - ORT_REF(p,i,Z);
#endif
      reduce_displacement(&d);
#ifdef TWOD
      len += sprintf( outbuf+len, "%12f %12f %12f %12f\n",
        ORT(p,i,X), ORT(p,i,Y), d.x, d.y);
#else
      len += sprintf( outbuf+len, "%12f %12f %12f %12f %12f %12f\n",
        ORT(p,i,X), ORT(p,i,Y), ORT(p,i,Z), d.x, d.y, d.z);
#endif
      /* flush or send outbuf if it is full */
      if (len > OUTPUT_BUF_SIZE - 256) flush_outbuf(out,&len,OUTBUF_TAG);
    }
  }
  flush_outbuf(out,&len,OUTBUF_TAG+1);
}

#endif

#ifdef DISLOC

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
      EPOT_REF(p,i) = POTENG(p,i);
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
      ORT_REF(p,i,X) = ORT(p,i,X);
      ORT_REF(p,i,Y) = ORT(p,i,Y);
#ifndef TWOD
      ORT_REF(p,i,Z) = ORT(p,i,Z);
#endif
    }
  }
}

#endif /* DISLOC */

#ifdef AVPOS

/******************************************************************************
*
* update_avpos updates avpos and resets av_epot
*
******************************************************************************/

void update_avpos(void)
{ 
  int k;
  for (k=0; k<ncells; k++) {
    int i;
    cell* p;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; i++) {
      AV_EPOT(p,i)  = POTENG(p,i);
      AV_POS(p,i,X) = ORT(p,i,X);
      AV_POS(p,i,Y) = ORT(p,i,Y);
      SHEET(p,i,X)  = 0.0;
      SHEET(p,i,Y)  = 0.0;
#ifndef TWOD
      AV_POS(p,i,Z) = ORT(p,i,Z);
      SHEET(p,i,Z)  = 0.0;
#endif
    }
  }
#ifdef NPT
  av_box_x.x = box_x.x;
  av_box_x.y = box_x.y;
  av_box_y.x = box_y.x;
  av_box_y.y = box_y.y;
#ifndef TWOD
  av_box_x.z = box_x.z;
  av_box_y.z = box_y.z;
  av_box_z.x = box_z.x;
  av_box_z.y = box_z.y;
  av_box_z.z = box_z.z;
#endif
#endif
}

/******************************************************************************
*
*  writes header for avp-file
*
******************************************************************************/

void write_header_avp(FILE *out)
{
  char c;
  time_t now;

  /* format line */
  /* binary_io for checkpoints not implemented
   if (binary_io)
     c = is_big_endian ? 'B' : 'L';
   else
    */
    c = 'A';
  fprintf(out, "#F %c 0 0 0 %d 0 1\n", c, DIM);
  
  /* contents line */
#ifdef TWOD
  fprintf(out, "#C number type mass x_av y_av Epot_av\n");
#else
  fprintf(out, "#C number type mass x_av y_av z_av Epot_av\n");
#endif

  /* box lines */
#ifdef TWOD
  fprintf(out, "#X \t%.16e %.16e\n", box_x.x , box_x.y);
  fprintf(out, "#Y \t%.16e %.16e\n", box_y.x , box_y.y);
#else
  fprintf(out, "#X \t%.16e %.16e %.16e\n", box_x.x , box_x.y , box_x.z);
  fprintf(out, "#Y \t%.16e %.16e %.16e\n", box_y.x , box_y.y , box_y.z);
  fprintf(out, "#Z \t%.16e %.16e %.16e\n", box_z.x , box_z.y , box_z.z);
#endif

  /* endheader line */
  time(&now);
  fprintf(out, "## Generated by %s on %s", progname, ctime(&now));
  fprintf(out, "#E\n");

}

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
  vektor avp_pos, coeff;

  for (k=0; k<ncells; k++) {
    cell* p;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; i++) {

	/* Averaged coordinates of atoms */
	avp_pos.x = AV_POS(p,i,X) * avpos_res / ( avpos_int + avpos_res );
	avp_pos.y = AV_POS(p,i,Y) * avpos_res / ( avpos_int + avpos_res );
#ifndef TWOD
	avp_pos.z = AV_POS(p,i,Z) * avpos_res / ( avpos_int + avpos_res );
#endif
	/* Coefficients of coordinates with respect to box vectors */
	coeff.x = SPROD( avp_pos, tbox_x );
	coeff.y = SPROD( avp_pos, tbox_y );
#ifndef TWOD
	coeff.z = SPROD( avp_pos, tbox_z );
#endif
	/* For periodic boundary conditions map coordinates into box */
	if( pbc_dirs.x == 1 ) 
	    coeff.x -= floor(coeff.x);
	if( pbc_dirs.y == 1 )
	    coeff.y -= floor(coeff.y);
#ifndef TWOD
	if( pbc_dirs.z == 1 ) 
	    coeff.z -= floor(coeff.z);
	
	x = coeff.x * box_x.x + coeff.y * box_y.x + coeff.z * box_z.x;
	y = coeff.x * box_x.y + coeff.y * box_y.y + coeff.z * box_z.y;
	z = coeff.x * box_x.z + coeff.y * box_y.z + coeff.z * box_z.z;
		
	len += sprintf( outbuf+len,
			"%d %d %12.16f %12.16f %12.16f %12.16f %12.16f\n",
			NUMMER(p,i), VSORTE(p,i), MASSE(p,i), 
			x, y, z, AV_EPOT(p,i) * avpos_res / avpos_int);
#else
	x = coeff.x * box_x.x + coeff.y * box_y.x;
	y = coeff.x * box_x.y + coeff.y * box_y.y;

	len += sprintf( outbuf+len,
			"%d %d %12.16f %12.16f %12.16f %12.16f\n",
			NUMMER(p,i), VSORTE(p,i), MASSE(p,i), x, y, 
			AV_EPOT(p,i) * avpos_res / (avpos_int + avpos_res ));
#endif
      /* flush or send outbuf if it is full */
      if (len > OUTPUT_BUF_SIZE - 256) flush_outbuf(out,&len,OUTBUF_TAG);
    }
  }
  flush_outbuf(out,&len,OUTBUF_TAG+1);
}

#endif /* AVPOS */

#ifdef FORCE

/******************************************************************************
*
*  writes header for force file
*  this is a special purpose header for potfit
*
******************************************************************************/

void write_header_force(FILE *out)
{
  /* number of atoms */
  fprintf(out, "%d\n", natoms);
  
  /* cohesive energy */
  fprintf(out, "%.16e\n",tot_pot_energy / natoms);

  /* box lines */
#ifdef TWOD
  fprintf(out, "%.16e %.16e\n", box_x.x, box_x.y );
  fprintf(out, "%.16e %.16e\n", box_y.x, box_y.y );
#else
  fprintf(out, "%.16e %.16e %.16e\n", box_x.x, box_x.y, box_x.z );
  fprintf(out, "%.16e %.16e %.16e\n", box_y.x, box_y.y, box_y.z );
  fprintf(out, "%.16e %.16e %.16e\n", box_z.x, box_z.y, box_z.z );
#endif

}

/******************************************************************************
*
*  filter function for write_config_select
*  writes forces to files *.nr.force
*
******************************************************************************/

void write_atoms_force(FILE *out)
{
  int k, len=0;

  for (k=0; k<ncells; k++) {
    cell* p;
    int i;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; i++) {
#ifdef TWOD
      len += sprintf( outbuf+len, "%d %.16e %.16e %.16e %.16e\n",
                      SORTE(p,i), ORT(p,i,X), ORT(p,i,Y),
                      KRAFT(p,i,X), KRAFT(p,i,Y) );
#else
      len += sprintf( outbuf+len, "%d %.16e %.16e %.16e %.16e %.16e %.16e\n",
                      SORTE(p,i), ORT(p,i,X), ORT(p,i,Y), ORT(p,i,Z),
                      KRAFT(p,i,X), KRAFT(p,i,Y), KRAFT(p,i,Z) );
#endif
      /* flush or send outbuf if it is full */
      if (len > OUTPUT_BUF_SIZE - 256) flush_outbuf(out,&len,OUTBUF_TAG);
    }
  }
  flush_outbuf(out,&len,OUTBUF_TAG+1);
}

#endif /* FORCE */

#ifdef ATDIST

/******************************************************************************
*
*  writes header for position file
*
******************************************************************************/

void write_header_atdist_pos(FILE *out)
{
  char c;
  time_t now;

  /* format line */
  /* binary_io for checkpoints not implemented
   if (binary_io)
     c = is_big_endian ? 'B' : 'L';
   else
    */
    c = 'A';
  fprintf(out, "#F %c 0 1 0 %d 0 0\n", c, DIM);
  
  /* contents line */
#ifdef TWOD
  fprintf(out, "#C type x y\n");
#else
  fprintf(out, "#C type x y z\n");
#endif

  /* box lines */
#ifdef TWOD
  fprintf(out, "#X %e 0.0\n", pic_ur.x - pic_ll.x);
  fprintf(out, "#Y 0.0 %e\n", pic_ur.y - pic_ll.y);
#else
  fprintf(out, "#X %e 0.0 0.0\n", pic_ur.x - pic_ll.x);
  fprintf(out, "#Y 0.0 %e 0.0\n", pic_ur.y - pic_ll.y);
  fprintf(out, "#Z 0.0 0.0 %e\n", pic_ur.z - pic_ll.z);
#endif

  /* endheader line */
  time(&now);
  fprintf(out, "## Generated by %s on %s", progname, ctime(&now));
  fprintf(out, "#E\n");
}

/******************************************************************************
*
*  filter function for write_config_select
*  writes types and positions, with periodic extension
*
******************************************************************************/

void write_atoms_atdist_pos(FILE *out)
{
  int k, len=0, ix, iy, iz;
  real x, y, z, t, co, si;

  co = cos(atdist_phi);
  si = sin(atdist_phi);

  for (k=0; k<ncells; k++) {
    cell* p;
    int i;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; i++)

      /* periodic continuation */
      for (ix=atdist_per_ll.x; ix<=atdist_per_ur.x; ix++)
#ifndef TWOD
        for (iz=atdist_per_ll.z; iz<=atdist_per_ur.z; iz++)
#endif
          for (iy=atdist_per_ll.y; iy<=atdist_per_ur.y; iy++) {
#ifdef TWOD
            x = ORT(p,i,X) + ix * box_x.x + iy * box_y.x;
            y = ORT(p,i,Y) + ix * box_x.y + iy * box_y.y;
#else
            x = ORT(p,i,X) + ix * box_x.x + iy * box_y.x + iz * box_z.x;
            y = ORT(p,i,Y) + ix * box_x.y + iy * box_y.y + iz * box_z.y;
            z = ORT(p,i,Z) + ix * box_x.z + iy * box_y.z + iz * box_z.y;
#endif
            t =  co * x + si * y;
            y = -si * x + co * y;
            x = t;

            /* continue if atom is not inside selected box */
            if ((x < pic_ll.x) || (x > pic_ur.x) ||
#ifndef TWOD
                (z < pic_ll.z) || (z > pic_ur.z) || 
#endif
                (y < pic_ll.y) || (y > pic_ur.y)) continue;

#ifdef TWOD
            len += sprintf( outbuf+len, "%d %e %e\n", SORTE(p,i), 
                            x - pic_ll.x, y - pic_ll.y );
#else
            len += sprintf( outbuf+len, "%d %e %e %e\n", SORTE(p,i), 
                            x - pic_ll.x, y - pic_ll.y, z - pic_ll.z );
#endif
            /* flush or send outbuf if it is full */
            if (len > OUTPUT_BUF_SIZE - 256) flush_outbuf(out,&len,OUTBUF_TAG);
	  }
  }
  flush_outbuf(out,&len,OUTBUF_TAG+1);
}

#endif /* ATDIST */

/******************************************************************************
*
*  write eng file header  - keep in sync with write_eng_file
******************************************************************************/

void write_eng_file_header()
{
  str255 fname;
  FILE *fl;

  int i;

  if (myid == 0) {

    sprintf(fname,"%s.eng",outfilename);
    fl = fopen(fname,"w");
    if (NULL == fl) error("Cannot open properties file.");

#ifdef CG
    fprintf(fl, "# ctf Epot fnorm");
#else
    fprintf(fl, "# time Epot ");
    fprintf(fl, "temperature ");

#if defined(STM) || defined(FRAC) 
    fprintf(fl, "stadiontemp ");
#endif

#ifdef FRAC
    fprintf(fl, "dampingtemp ");
#endif

#ifdef FTG
    for(i=0;i<nslices;i++)
      fprintf(fl, "temp_%d ", i);
#endif

#ifdef FNORM
    fprintf(fl, "fnorm ");
#endif
#ifdef GLOK
    fprintf(fl, "PxF ");
#endif
#ifdef EINSTEIN
    fprintf(fl, "omega_E ");
#endif
    fprintf(fl, "pressure ");
    fprintf(fl, "volume ");
#if defined(NVT) || defined(NPT) || defined(STM) 
    fprintf(fl, "eta * tau_eta ");
#endif
#ifdef FRAC  
    fprintf(fl, "gamma_damp ");
    fprintf(fl, "strainrate ");
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
#ifdef STRESS_TENS
    fprintf(fl, "Press_xx Press_yy ");
#ifdef TWOD
    fprintf(fl, "Press_xy ");
#else 
    fprintf(fl, "Press_zz ");
    fprintf(fl, "Press_yz Press_xz Press_xy");
#endif    
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

void write_eng_file(int steps)
{
  str255 fname;
  int i;
  static int flush_count=0;
  real tmp;  

#ifdef HPO
  char *format=" %.16e";
#else
  char *format=" %e";
#endif

  real Epot, Temp, vol;

#if defined(STM) || defined(FRAC)
  real Temp_damp, Temp_stadium = 0.0;
#endif

#ifdef STRESS_TENS
  real Press_xx,Press_yy, Press_xy;
#ifndef TWOD
  real Press_zz,Press_yz, Press_zx;
#endif
  calc_tot_presstens();
#endif

  /* write only on CPU 0; 
     calc_tot_presstensor() above must be executed on all CPUs */
  if (myid>0) return;

#ifdef STRESS_TENS
  Press_xx = tot_presstens.xx / volume; 
  Press_yy = tot_presstens.yy / volume; 
#ifndef TWOD
  Press_zz = tot_presstens.zz / volume;
  Press_yz = tot_presstens.yz / volume; 
  Press_zx = tot_presstens.zx / volume; 
#endif
  Press_xy = tot_presstens.xy / volume; 
#endif

  Epot =       tot_pot_energy / natoms;

#ifndef CG
#ifdef UNIAX
  Temp = 2.0 * tot_kin_energy / (nactive + nactive_rot);
#else
  Temp = 2.0 * tot_kin_energy / nactive;
#endif
#endif

#ifdef STM 
  Temp     = 2.0 * tot_kin_energy / (nactive - n_stadium);
#endif 

#if defined(STM) || defined(FRAC)
  if(n_stadium != 0)  Temp_stadium = 2.0 * E_kin_stadium / n_stadium;
#endif

#if defined(FRAC)
  if(sum_f !=0){
      Temp_damp = 2.0 * E_kin_damp / (sum_f * DIM);
  } else {
      Temp_damp = 0.0;
  }
#endif

  vol  = volume / natoms;
#ifndef CG
  pressure = Temp / vol + virial / (DIM * volume);
#endif

  /* open .eng file if it is not yet open */
  if (NULL == eng_file) {
    sprintf(fname,"%s.eng",outfilename);
    eng_file = fopen(fname,"a");
    if (NULL == eng_file) error("Cannot open properties file.");
  }

#ifdef CG
  fprintf(eng_file, "%d",     steps);
  fprintf(eng_file, " %.16e", (double) Epot);
  fprintf(eng_file, format,   (double) sqrt( (double) fnorm / nactive ) );
#else
  fprintf(eng_file, "%e",     (double) (steps * timestep));
  fprintf(eng_file, " %.16e", (double) Epot);
  fprintf(eng_file, format,   (double) Temp);
#if defined(STM) || defined(FRAC)
  fprintf(eng_file, format,   (double) Temp_stadium);
#endif

#ifdef FRAC
  fprintf(eng_file, format,   (double) Temp_damp);
#endif

#ifdef FTG
  for(i=0;i<nslices;i++){
    Temp =  0.0;
    if (0 !=  *(ninslice + i))
      Temp =  2.0* *(E_kin_ftg+i)/ *(ninslice + i);
    fprintf(eng_file, format, Temp); 
  }
#endif


#ifdef FNORM
  fprintf(eng_file, format,   (double) sqrt( (double) fnorm / nactive ) );
#endif
#ifdef GLOK
  fprintf(eng_file, format,   (double) PxF);
#endif
#ifdef EINSTEIN
  fprintf(eng_file, format,   sqrt((double) omega_E / (nactive * Temp)));
#endif
  fprintf(eng_file," %e",     (double) pressure);
  fprintf(eng_file," %e",     (double) vol);
#if defined(NVT) || defined(NPT) || defined(STM)
  fprintf(eng_file," %e",     (double) (eta * tau_eta) );
#endif
#ifdef FRAC 
  fprintf(eng_file," %e",     (double) gamma_damp );
  fprintf(eng_file," %e",     (double) dotepsilon );
#endif

  if (ensemble==ENS_NPT_AXIAL) {
#ifdef TWOD
    fprintf(eng_file," %e %e", (double) stress_x, (double) stress_y );
    fprintf(eng_file," %e %e", (double)  box_x.x, (double)  box_y.y );
#else
    fprintf(eng_file," %e %e %e", 
	    (double) stress_x, (double) stress_y, (double) stress_z );
    fprintf(eng_file," %e %e %e", 
	    (double)  box_x.x, (double)  box_y.y, (double)  box_z.z );
#endif
  }
#ifdef STRESS_TENS
  fprintf(eng_file," %e %e", (double) Press_xx, (double) Press_yy);
#ifdef TWOD
  fprintf(eng_file," %e", (double) Press_xy);
#else 
  fprintf(eng_file," %e", (double) Press_zz);
  fprintf(eng_file," %e %e %e",
	  (double) Press_yz, (double) Press_zx, (double) Press_xy);
#endif    
#endif
#endif /* CG */
  putc('\n',eng_file);
  flush_count++;

  /* flush .eng file every flush_int writes */
  if (flush_count > flush_int) {
    fflush(eng_file);
    flush_count=0;
  }

}

#ifdef MSQD

/******************************************************************************
*
* write_msqd writes mean square displacement to *.msqd file
*
******************************************************************************/

void write_msqd(int steps)
{
  str255 fname;
  int i,j;
  static int flush_count=0;

  /* open .msqd file if not yet open */
  if (NULL == msqd_file) {
    sprintf(fname,"%s.msqd",outfilename);
    msqd_file = fopen(fname,"a");
    if (NULL == msqd_file) error("Cannot open msqd file.");
  }

  /* write the mean square displacements */
  fprintf(msqd_file, "%10.4e", (double)(steps * timestep));
  for (i=0; i<ntypes; i++) 
    for (j=0; j<DIM; j++) {
    fprintf(msqd_file," %10.4e", (double)(msqd_global[i*DIM+j] / num_sort[i]));
  }
  putc('\n',msqd_file);
  flush_count++;

  /* flush .msqd file every flush_int writes */
  if (flush_count > flush_int){
    fflush(msqd_file);
    flush_count=0;
  }

}

/******************************************************************************
*
*  writes header for sqd-file
*
******************************************************************************/

void write_header_sqd(FILE *out)
{
  char c;
  time_t now;

  /* format line */
  /* binary_io for checkpoints not implemented
   if (binary_io)
     c = is_big_endian ? 'B' : 'L';
   else
    */
  c = 'A';
  fprintf(out, "#F %c 0 1 0 0 0 %d\n", c, DIM);
  
  /* contents line */
#ifdef TWOD
  fprintf(out, "#C type d^2_x d^2_y\n");
#else
  fprintf(out, "#C type d^2_x d^2_y d^2_z\n");
#endif

  /* endheader line */
  time(&now);
  fprintf(out, "## Generated by %s on %s", progname, ctime(&now));
  fprintf(out, "#E\n");
}

/******************************************************************************
*
*  filter function for write_config_select
*  writes average position to files *.sqd
*
******************************************************************************/

void write_atoms_sqd(FILE *out)
{
  int i, k, len=0;
  vektor d;

  for (k=0; k<ncells; k++) {
    cell* p;
    p = cell_array + CELLS(k);
    for (i=0; i<p->n; i++) {
      d.x = SQR( ORT(p,i,X) - REF_POS(p,i,X) );
      d.y = SQR( ORT(p,i,Y) - REF_POS(p,i,Y) );
#ifndef TWOD
      d.z = SQR( ORT(p,i,Z) - REF_POS(p,i,Z) );
#endif
#ifdef TWOD
      len += sprintf( outbuf+len, "%d %e %e\n",    SORTE(p,i), d.x, d.y     );
#else
      len += sprintf( outbuf+len, "%d %e %e %e\n", SORTE(p,i), d.x, d.y, d.z);
#endif
      /* flush or send outbuf if it is full */
      if (len > OUTPUT_BUF_SIZE - 256) flush_outbuf(out,&len,OUTBUF_TAG);
    }
  }
  flush_outbuf(out,&len,OUTBUF_TAG+1);
}

#endif /* MSQD */

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

/******************************************************************************
*
* write_config_header writes a header to a config file
*
******************************************************************************/

void write_header_config(FILE *out)
{
  char c;
  time_t now;

  /* format line */
  /* binary_io for checkpoints not implemented
   if (binary_io)
     c = is_big_endian ? 'B' : 'L';
   else
    */
    c = 'A';
#ifdef UNIAX
  fprintf(out, "#F %c 1 1 2 %d %d 1\n", c, 4*DIM, 2*DIM);
#else
  fprintf(out, "#F %c 1 1 1 %d %d ", c, DIM, DIM);
#ifdef ORDPAR
  fprintf(out, "2\n" );
#else
  fprintf(out, "1\n" );
#endif
#endif

  /* contents line */
#ifdef CG
  fprintf(out, "#C ask Erik\n");
#elif defined(UNIAX)
  fprintf(out, "#C number type mass inertia pos(3) axis(3) shape(3)");
  fprintf(out, " pot_depth(3) velocities(3) ang_veloc(3) Epot\n");
#else
#ifdef TWOD
  fprintf(out, "#C number type mass x y vx vy");
#else
  fprintf(out, "#C number type mass x y z vx vy vz");
#endif
#ifdef ORDPAR
  fprintf(out, " ordpar nbanz" );
#else
  fprintf(out, " Epot" );
#endif
  fprintf(out, "\n" );
#endif


  /* box lines */
#ifdef TWOD
  fprintf(out, "#X \t%.16e %.16e\n", box_x.x , box_x.y);
  fprintf(out, "#Y \t%.16e %.16e\n", box_y.x , box_y.y);
#else
  fprintf(out, "#X \t%.16e %.16e %.16e\n", box_x.x , box_x.y , box_x.z);
  fprintf(out, "#Y \t%.16e %.16e %.16e\n", box_y.x , box_y.y , box_y.z);
  fprintf(out, "#Z \t%.16e %.16e %.16e\n", box_z.x , box_z.y , box_z.z);
#endif

  /* endheader line */
  time(&now);
  fprintf(out, "## Generated by %s on %s", progname, ctime(&now));
  fprintf(out, "#E\n");

}

