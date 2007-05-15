
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
*  Routines for dealing with IMD configurations
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "misc.h"
#include "util.h"
#include "conf_tools.h"

/******************************************************************************
*
* read the header of a config file
*
******************************************************************************/

int read_header(header_info_t *info, str255 infilename)
{
  FILE   *infile;
  str255 line, fname, str;
  int    have_format=0, have_header=0;
  int    p, n, t, m, np, nv, nd, i;

  infile = fopen(infilename,"r");
  if (NULL==infile) error_str("cannot open input file %s", infilename);

  info->n_lines = 0;
  fgets(line, 255, infile);
  while (line[0]=='#') {
    /* format line */
    if      (line[1]=='F') {
      p = sscanf(line+2, "%s %d %d %d %d %d %d", str,&n,&t,&m,&np,&nv,&nd );
      if (p<7) error_str("Format line of file %s corrupt", infilename);
      have_format    = 1;
      info->format   = str[0];
      info->n_number = n;
      info->n_type   = t;
      info->n_mass   = m;
      info->n_pos    = np;
      info->n_vel    = nv;
      info->n_data   = nd;
      info->n_items  = n + t + m + np + nv + nd;
      if ((info->format == 'B') || (info->format == 'b'))
        info->endian = 1;
      else
        info->endian = 0;
    }
    /* contents line */
    else if (line[1]=='C') {
      char *token = strtok(line+2, " \t\r\n");
      for (i=1; i < n+t+m+np+nv; i++) {
        token = strtok(NULL, " \t\r\n");
      }
      for (i=0; i < nd; i++) {
        token = strtok(NULL, " \t\r\n");
        info->contents[i] = strdup(token);
      }
    }
    /* endheader line */
    else if (line[1]=='E') {
      if (have_format) have_header = 1;
    }
    info->lines[info->n_lines++] = strdup(line);
    fgets(line, 255, infile);
    if (feof(infile)) break;
  }

  /* limited backwards compatibility, it there is no header */
  if (have_header==0) {
    int    p, n, s;
    double d[MAX_ITEMS_CONFIG];
    info->format   = 'A';
    info->n_number = 1;
    info->n_type   = 1;
    info->n_mass   = 1;
    info->n_pos    = DIM;
    p = sscanf( line,
      "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
       &n,&s,&d[0],&d[1],&d[2],&d[3],&d[4],&d[5],&d[6],&d[7],&d[8],
       &d[9],&d[10],&d[11],&d[12],&d[13],&d[14],&d[15]);
    if (p < 3 + DIM)
      error("corrupt first line in config file");
    if (p >= 3+2*DIM) {
      info->n_vel  = DIM;
      info->n_data = p - (3+2*DIM);
      if (info->n_data > 0) {
        info->contents[0] = "Epot";
        for (i=1; i<info->n_data; i++) info->contents[i] = "unknown";
      }
    }
    else {
      info->n_vel  = 0;
      info->n_data = 0;
    }
    info->n_items  = 3 + DIM + info->n_vel + info->n_data;
  }
  fclose(infile);

  /* check whether file contains what we need */
  if (have_header) {
    if (info->n_items > MAX_ITEMS_CONFIG + 2)
      error("too many items per atom in config file");
  } 
  return have_header;
}

/******************************************************************************
*
*  map vektor back into simulation box
*
******************************************************************************/

#ifdef C2C

#define back_into_box(pos) (pos)

#else

vektor back_into_box(vektor pos)
{
  double i;

  if (pbc_dirs.x==1) {
    i = floor(SPROD(pos,tbox_x));
    pos.x  -= i *  box_x.x;
    pos.y  -= i *  box_x.y;
#ifndef TWOD
    pos.z  -= i *  box_x.z;
#endif
  }

  if (pbc_dirs.y==1) {
    i = floor(SPROD(pos,tbox_y));
    pos.x  -= i *  box_y.x;
    pos.y  -= i *  box_y.y;
#ifndef TWOD
    pos.z  -= i *  box_y.z;
#endif
  }

#ifndef TWOD
  if (pbc_dirs.z==1) {
    i = floor(SPROD(pos,tbox_z));
    pos.x  -= i *  box_z.x;
    pos.y  -= i *  box_z.y;
    pos.z  -= i *  box_z.z;
  }
#endif
  return pos;

}

#endif

/******************************************************************************
*
* read one atom from a config file
*
******************************************************************************/

int read_atom(header_info_t *info, FILE *infile, atom_t *atom)
{
  double m, d[MAX_ITEMS_CONFIG];
  vektor pos;
  char buf[1024];
  int  p=0, n, s, i, k, is_big_endian = endian();

  /* ASCII input */
  if (info->format == 'A') {
    p = 1;
    if (NULL==fgets(buf,sizeof(buf),infile)) p=0;
    /* eat comments */
    while (('#'==buf[0]) && !feof(infile))
      if (NULL==fgets(buf,sizeof(buf),infile)) p=0;
    if (p==0) return 0;
    p = sscanf( buf,
      "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
       &n,&s,&d[0],&d[1],&d[2],&d[3],&d[4],&d[5],&d[6],&d[7],&d[8],
       &d[9],&d[10],&d[11],&d[12],&d[13],&d[14],&d[15]);
  }
  /* double precision input */
  else if ((info->format=='B') || (info->format=='L')) {
    i_or_d *data = (i_or_d *) buf;
    p = fread(buf, sizeof(i_or_d), info->n_items-1, infile);
    if (p>0) p++; /* first value contains two items */
    if (info->endian == is_big_endian) {
      n = data[0].i[0];
      s = data[0].i[1];
      for (k=0; k < info->n_items-2; k++) { 
        d[k] = (double) data[k+1].d;
      }
    }
    else {
      n = SwappedInteger(data[0].i[0]);
      s = SwappedInteger(data[0].i[1]);
      for (k=0; k < info->n_items-2; k++) 
        d[k] = (double) SwappedDouble(data[k+1].d);
    }
  }
  /* single precision input */
  else if ((info->format=='b') || (info->format=='l')) {
    i_or_f *data = (i_or_f *) buf;
    p = fread(buf, sizeof(i_or_f), info->n_items, infile);
    if (info->endian == is_big_endian) {
      n = data[0].i;
      s = data[1].i;
      for (k=0; k < info->n_items-2; k++) 
        d[k] = (double) data[k+2].f;
    }
    else {
      n = SwappedInteger(data[0].i);
      s = SwappedInteger(data[1].i);
      for (k=0; k < info->n_items-2; k++) 
        d[k] = (double) SwappedFloat(data[k+2].f);
    }
  }

  atom->number = n;
  atom->type = s;
  k = 0;
  if (info->n_mass == 1) atom->mass = d[k++];
  if (info->n_pos == 2) {
    pos.x = d[k++];
    pos.y = d[k++];
    atom->pos = back_into_box(pos);
  }
#ifndef TWOD
  else if (info->n_pos == 3) {
    pos.x = d[k++];
    pos.y = d[k++];
    pos.z = d[k++];
    atom->pos = back_into_box(pos);
  }
#endif
  if (info->n_vel == 2) {
    atom->vel.x = d[k++];
    atom->vel.y = d[k++];
  }
#ifndef TWOD
  else if (info->n_vel == 3) {
    atom->vel.x = d[k++];
    atom->vel.y = d[k++];
    atom->vel.z = d[k++];
  }
#endif
  for (i=0; i<info->n_data; i++)
    atom->data[i] = d[k++];

  return p;

}

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

int SwappedInteger( int i )
{
  union {
    int i;
    unsigned char b[4];
  } dat1, dat2;
  dat1.i = i;
  dat2.b[0] = dat1.b[3];
  dat2.b[1] = dat1.b[2];
  dat2.b[2] = dat1.b[1];
  dat2.b[3] = dat1.b[0];
  return dat2.i;
}

float SwappedFloat( float f )
{
  union {
    float f;
    unsigned char b[4];
  } dat1, dat2;
  dat1.f = f;
  dat2.b[0] = dat1.b[3];
  dat2.b[1] = dat1.b[2];
  dat2.b[2] = dat1.b[1];
  dat2.b[3] = dat1.b[0];
  return dat2.f;
}

double SwappedDouble( double d )
{
  union {
    double d;
    unsigned char b[8];
  } dat1, dat2;
  dat1.d = d;
  dat2.b[0] = dat1.b[7];
  dat2.b[1] = dat1.b[6];
  dat2.b[2] = dat1.b[5];
  dat2.b[3] = dat1.b[4];
  dat2.b[4] = dat1.b[3];
  dat2.b[5] = dat1.b[2];
  dat2.b[6] = dat1.b[1];
  dat2.b[7] = dat1.b[0];
  return dat2.d;
}

