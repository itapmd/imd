
/******************************************************************************
*
*  IMD -- The ITAP Molecular Dynamics Program
*
*  Copyright 1996-2006 Institute for Theoretical and Applied Physics,
*  University of Stuttgart, D-70550 Stuttgart
*
*  $Revision$
*  $Date$
*
******************************************************************************/

#define SPROD(v,w) (v.x*w.x + v.y*w.y + v.z*w.z)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "misc.h"
#include "conf_tools.h"

void usage(char *progname) 
{
  printf("\n");
  printf("   Usage: %s [options] <infile> <outfile>\n\n", progname);
  printf("   Converts IMD configuration files between different formats,\n");
  printf("   and possibly strips off unwanted data.\n\n");
  printf("   Options:  -f <F>     Output format; <F> can be:\n\n");
  printf("                          - A for ASCII (default)\n");
  printf("                          - l for single precision binary (little endian)\n");
  printf("                          - b for single precision binary (big endian)\n");
  printf("                          - L for double precision binary (little endian)\n");
  printf("                          - B for double precision binary (big endian)\n\n");
  printf("             -v         Include velocities in output\n\n");
  printf("             -m         Include mass in output\n\n");
  printf("             -k         Compute Ekin from mass and velocities, and append it\n\n");
  printf("             -d <cols>  Comma separated indices of data columns to be\n");
  printf("                        included in output, like 0,2\n\n");
  printf("             -h         This help\n\n");
  exit(1);
}

int main(int argc, char **argv) 
{
  char *progname, *infilename, *outfilename, *str, *token, format='A';
  char outbuf[OUTPUT_BUF_SIZE], line[1024];
  int  with_vel=0, with_Ekin=0, with_mass=0;
  int  n_data_out=0, data_out[MAX_ITEMS_CONFIG];
  int  i, p, n, my_endian, out_endian, len=0, natoms=0, have_header;
  FILE *infile, *outfile;
  header_info_t info;
  atom_t atom;

  /* parse command line options */
  progname = strdup(argv[0]);
  while ((argc > 1) && (argv[1][0] =='-')) {
    if (argv[1][1]=='f') {
      format = argv[2][0];
      if ((format != 'A') && (format != 'B') && (format != 'b')  
                          && (format != 'L') && (format != 'l') )
        error("unkown output format");
      argc -= 2;
      argv += 2;
    }
    else if (argv[1][1]=='v') {
      with_vel = 1;
      argc -= 1;
      argv += 1;
    }
    else if (argv[1][1]=='k') {
      with_Ekin = 1;
      argc -= 1;
      argv += 1;
    }
    else if (argv[1][1]=='m') {
      with_mass = 1;
      argc -= 1;
      argv += 1;
    }
    else if (argv[1][1]=='d') {
      str = argv[2];
      token = strtok(str, ",");
      while (token!=NULL) {
        data_out[n_data_out++] = atoi(token);
        token = strtok(NULL, ",");
      }
      argc -= 2;
      argv += 2;
    }
    else if (argv[1][1]=='h') {
      usage(progname);
    }
  }
  if (argc<3) {
    usage(progname);
  }

  infilename  = strdup(argv[1]);
  outfilename = strdup(argv[2]);
  my_endian = endian();
  if ((format=='B') || (format=='b')) out_endian=1;
  else                                out_endian=0;

  /* open input file */
  infile = fopen(infilename,"r");
  if (NULL==infile) error_str("Cannot open atoms file %s", infilename);

  /* read file header */
  have_header = read_header(&info, infilename);
  /* eat header, if there is one */
  if (have_header) {
    fgets(line,sizeof(line),infile);
    while (('#'==line[0]) && ('E'!=line[1]) && !feof(infile)) {
      fgets(line,sizeof(line),infile);
    }
  }
  if ((with_Ekin) && ((info.n_mass==0) || (info.n_vel==0)))
    error("Cannot compute kinetic energy");
  if ((with_vel) && (info.n_vel==0))
    error("No velocities available");
  if ((with_mass) && (info.n_mass==0))
    error("No mass available");

  /* open output file */
  outfile = fopen(outfilename,"w");
  if (NULL==outfile) error_str("Cannot open output file %s", outfilename);

  /* write output header */
  /* format line */
  fprintf(outfile, "#F %c %d %d", format, info.n_number, info.n_type);
  if (with_mass)
    fprintf(outfile, " 1 %d", info.n_pos);
  else
    fprintf(outfile, " 0 %d", info.n_pos);
  if (with_vel)
    fprintf(outfile, " %d %d\n", info.n_vel, n_data_out + with_Ekin);
  else
    fprintf(outfile, " 0 %d\n", n_data_out + with_Ekin);
  /* contents line */
  fprintf(outfile, "#C");
  if (info.n_number) fprintf(outfile, " number");
  if (info.n_type) fprintf(outfile, " type");
  if (with_mass) fprintf(outfile, " mass");
  if      (info.n_pos==2) fprintf(outfile, " x y");
  else if (info.n_pos==3) fprintf(outfile, " x y z");
  if (with_vel) {
    if      (info.n_vel==2) fprintf(outfile, " vx vy");
    else if (info.n_vel==3) fprintf(outfile, " vx vy vz");
  }
  for (i=0; i<n_data_out; i++) {
    if (data_out[i] >= info.n_data)
      error("non-existing data column");
    fprintf(outfile, " %s", info.contents[data_out[i]]);
  }
  if (with_Ekin) fprintf(outfile, " Ekin");
  fprintf(outfile, "\n");
  /* copy all other lines */
  for (i=0; i<info.n_lines; i++) 
    if ((info.lines[i][1] != 'F') && (info.lines[i][1] != 'C')
      && (info.lines[i][1] != 'E')) fprintf(outfile, "%s", info.lines[i]);
  fprintf(outfile, "#E\n");

  /* process the config file */
  while (!feof(infile)) {

    p = read_atom(&info, infile, &atom);
    if (p==0) break;
    if (info.n_number == 0) atom.number = natoms++;
    if (info.n_vel == 2) atom.vel.z = 0.0;
    if (with_Ekin) atom.Ekin = 0.5 * atom.mass * SPROD(atom.vel,atom.vel);

    /* 32-bit binary output */
    if ((format=='l') || (format=='b')) {
      i_or_f *data = (i_or_f *) (outbuf+len);
      n = 0;
      if (info.n_number) data[n++].i = (int)   atom.number;
      if (info.n_type)   data[n++].i = (int)   atom.type;
      if (with_mass)     data[n++].f = (float) atom.mass;
      if (info.n_pos==2) {
        data[n++].f = (float) atom.pos.x;
        data[n++].f = (float) atom.pos.y;
      }
      else if (info.n_pos==3) {
        data[n++].f = (float) atom.pos.x;
        data[n++].f = (float) atom.pos.y;
        data[n++].f = (float) atom.pos.z;
      }
      if (with_vel) { 
        if (info.n_pos==2) {
          data[n++].f = (float) atom.vel.x;
          data[n++].f = (float) atom.vel.y;
        }
        else if (info.n_pos==3) {
          data[n++].f = (float) atom.vel.x;
          data[n++].f = (float) atom.vel.y;
          data[n++].f = (float) atom.vel.z;
        }
      }
      for (i=0; i<n_data_out; i++)
        data[n++].f = (float) atom.data[data_out[i]];
      if (with_Ekin)
        data[n++].f = (float) atom.Ekin;
      if (my_endian != out_endian)
        for (i=0; i<n; i++) data[i].f = SwappedFloat( data[i].f );
      len += n * sizeof(i_or_f); 
    }

    /* 64-bit binary output */
    if ((format=='L') || (format=='B')) {
      i_or_d *data = (i_or_d *) (outbuf+len);
      n = 0;
      if (info.n_type) {
        data[n].i[0] = (int) atom.number;
        data[n].i[1] = (int) atom.type;
        n++;
      }
      if (with_mass) data[n++].d = (double) atom.mass;
      if (info.n_pos==2) {
        data[n++].d = (double) atom.pos.x;
        data[n++].d = (double) atom.pos.y;
      }
      else if (info.n_pos==3) {
        data[n++].d = (double) atom.pos.x;
        data[n++].d = (double) atom.pos.y;
        data[n++].d = (double) atom.pos.z;
      }
      if (with_vel) { 
        if (info.n_pos==2) {
          data[n++].d = (double) atom.vel.x;
          data[n++].d = (double) atom.vel.y;
        }
        else if (info.n_pos==3) {
          data[n++].d = (double) atom.vel.x;
          data[n++].d = (double) atom.vel.y;
          data[n++].d = (double) atom.vel.z;
        }
      }
      for (i=0; i<n_data_out; i++)
        data[n++].d = (double) atom.data[data_out[i]];
      if (with_Ekin)
        data[n++].d = (double) atom.Ekin;
      if (my_endian != out_endian) {
        if (info.n_type) {
          data[0].i[0] = SwappedInteger( data[0].i[0] );
          data[0].i[1] = SwappedInteger( data[0].i[1] );
          for (i=1; i<n; i++) data[i].d = SwappedDouble( data[i].d );
	} 
        else {
          for (i=0; i<n; i++) data[i].d = SwappedDouble( data[i].d );
	}
      }
      len += n * sizeof(i_or_d); 
    }

    /* ASCII output */
    else if (format=='A') {
      if (info.n_number) 
        len += sprintf( outbuf+len, " %d", atom.number);
      if (info.n_type)
        len += sprintf( outbuf+len, " %d", atom.type);
      if (with_mass)
        len += sprintf( outbuf+len, " %e", atom.mass);
      if (info.n_pos==2) {
        len += sprintf( outbuf+len, " %e %e", atom.pos.x, atom.pos.y);
      }
      else if (info.n_pos==3) {
        len += sprintf( outbuf+len, " %e %e %e", 
                        atom.pos.x, atom.pos.y, atom.pos.z);
      }
      if (with_vel) { 
        if (info.n_vel==2) {
          len += sprintf( outbuf+len, " %e %e", atom.vel.x, atom.vel.y);
        }
        else if (info.n_vel==3) {
          len += sprintf( outbuf+len, " %e %e %e", 
                          atom.vel.x, atom.vel.y, atom.vel.z);
        }
      }
      for (i=0; i<n_data_out; i++)
        len += sprintf( outbuf+len, " %e", atom.data[data_out[i]]);
      if (with_Ekin)
        len += sprintf( outbuf+len, " %e", atom.Ekin);
      len += sprintf( outbuf+len, "\n");
    }

    /* flush output */
    if (len > OUTPUT_BUF_SIZE - 256) {
      fwrite(outbuf, sizeof(char), len, outfile);
      len = 0;
    }
    natoms++;
  }
  fwrite(outbuf, sizeof(char), len, outfile);
  fclose(outfile);
  fclose(infile);
 
  return 0;
}
