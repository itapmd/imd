
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define BUFSIZE 1024

int main (int argc, char **argv) 
{
  FILE  *infile, *outfile;
  char  *infilename, *outfilename;
  float Emin, Emax, scale, tmp;
  int   i,readcount, writecount;
  float floatbuf[BUFSIZE];
  unsigned char charbuf[BUFSIZE];

  if (argc < 3) {
    printf("Usage: %s Infile Emin Emax\n", argv[0]);
    exit(1);
  }

  infilename  = argv[1];
  outfilename = (char *)malloc(strlen(infilename)+5);
  strcpy(outfilename,infilename);
  strcat(outfilename,".dat");
  sscanf(argv[2], "%f", &Emin);
  sscanf(argv[3], "%f", &Emax); 

  infile  = fopen(infilename,"r");
  if (infile==NULL) {
    printf("Cannot open input file %s\n", infilename );
    exit(1);
  }

  outfile = fopen(outfilename,"w");
  if (outfile==NULL) {
    printf("Cannot open output file %s\n", outfilename );
    exit(1);
  }

  scale = 255.9999/(Emax-Emin);
  do {
    readcount = fread(floatbuf,sizeof(float),BUFSIZE,infile);
    if ((readcount<BUFSIZE) && (!feof(infile))) {
      printf("reading error\n");
      exit(1);
    }
    if (readcount>0) {
      for (i=0; i<readcount; i++) {
        tmp = (floatbuf[i]<Emin) ? Emin : floatbuf[i];
        tmp = (tmp>Emax)         ? Emax : tmp; 
        charbuf[i] = (unsigned char) scale*(tmp-Emin);
      }
      writecount = fwrite(charbuf,sizeof(unsigned char),readcount,outfile);
      if (writecount<readcount) {
        printf("writing error\n");
        exit(1);
      }
    }
  } while (!feof(infile));

  close(infile);
  close(outfile);

}
  
