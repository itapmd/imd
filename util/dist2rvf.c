/******************************************************************************
*
*  dist2rvf is a utility program to convert distribution files to 
*  virvo .rvf format.
*
*  Compilation:
*
*    gcc -o dist2rvf dist2rvf.c
*
*  Usage:
*
*     dist2rvf <distfile>
*
*  It writes to the file dist2rvf. The distribution file must contain the
*  standard header.
*
******************************************************************************/

#include <stdio.h>
#include <limits.h>

void error(char *errstr);

void error(char *errstr) {
  fprintf(stderr, "%s\n", errstr);
  exit(-1);
}

int main(int argc, char *argv[]) {
  char line[255], distname[255], rvfname[255], sdum[255];
  FILE *fpin, *fpout;
  int dimx, dimy, dimz;
  float scale, val, min, max;
  unsigned char c;
  int i=0, n;

  if (argc != 2)
    error("Usage: dist2rvf <filename>");

  sprintf(distname, "%s", argv[1]);
  sprintf(rvfname,  "%s.rvf", argv[1]);
  fpin = fopen(distname, "r");
  if (fpin==NULL) error("File not found!\n");

  while ((fgets(line, 255, fpin))&&(line[0]=='#')) {
      if (line[1]=='D') {
	sscanf(line, "%s %d %d %d", sdum, &dimx, &dimy, &dimz);
	continue;
      }	else continue;
  }

  min = 100000.0;
  max = -100000.0;

  do {
    sscanf(line, "%f", &val);
    if (min>val) min=val;
    if (max<val) max=val;
  } while (fgets(line, 255, fpin));
  scale = 256.0/(max-min);

  fclose(fpin);

  fpin = fopen(distname, "r");
  fpout = fopen(rvfname, "w");

  c = (unsigned char)((dimx & 65280)>>8);
  fputc(c, fpout);
  c = (unsigned char)(dimx & 255); 
  fputc(c, fpout);
  c = (unsigned char)((dimy & 65280)>>8);
  fputc(c, fpout);
  c = (unsigned char)(dimy & 255); 
  fputc(c, fpout);
  c = (unsigned char)((dimz & 65280)>>8);
  fputc(c, fpout);
  c = (unsigned char)(dimz & 255); 
  fputc(c, fpout);
  while (fgets(line, 255, fpin)) {
    if (line[0]=='#') continue;
    n = sscanf(line, "%f", &val);
    if (n==-1) continue;
    val -= min;
    val *= scale;
    c = (unsigned char)(val);
    fputc(c, fpout);
    i++;
  }

  printf("%d %d %d %d\n", i, dimx, dimy, dimz);
  fclose(fpin);
  fclose(fpout);
}

