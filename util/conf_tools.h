
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

#define MAX_ITEMS_CONFIG 16
#define MAX_LINES_HEADER 100
#define OUTPUT_BUF_SIZE  1000000

#ifndef UTIL_H

/* if util.h is not read, do some of the declarations here */

#ifdef TWOD
#define DIM 2
typedef struct {
  double x, y;
} vektor;
#else
#define DIM 3
typedef struct {
  double x, y, z;
} vektor;
#endif

typedef char str255[255];

#endif

typedef struct {
  char format;
  int  endian;
  int  n_number;
  int  n_type;
  int  n_mass;
  int  n_pos;
  int  n_vel;
  int  n_data;
  int  n_items;
  int  n_lines;
  char *contents[MAX_ITEMS_CONFIG];
  char *lines[MAX_LINES_HEADER];
} header_info_t;

typedef struct {
  int    number, type;
  double mass, Ekin;
  vektor pos, vel;
  double data[MAX_ITEMS_CONFIG];
} atom_t;

typedef union {
  int   i;
  float f;
} i_or_f;

typedef union {
  int    i[2];
  double d;
} i_or_d;


int    read_header(header_info_t *, str255);
int    read_atom(header_info_t *, FILE *, atom_t *);
int    SwappedInteger(int i);
float  SwappedFloat(float f);
double SwappedDouble(double d);
int    endian(void);
