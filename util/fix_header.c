#include <stdio.h>
#include <string.h>

/* replace the first header line (in place) by the value of REPLACEMENT */

#define PRESS_DIST_REP "#F A 3 0 9"

#define SHOCK_DIST_REP "#F A 3 0 13"

#define REPLACEMENT PRESS_DIST_REP

void error(char *msg)
{
  fprintf(stderr,"Error: %s\n",msg);
  exit(2);
}

int main( int argc, char **argv ) 
{
  FILE *fp;
  char line[255];
  int ln1, ln2, i;

  fp = fopen(argv[1],"r+");
  fgets(line, 254,fp);
  ln1 = strlen(line);
  ln2 = strlen(REPLACEMENT);
  if (ln2>ln1-1) error("replacement header too long");
  sprintf(line,"%s",REPLACEMENT);
  for (i=ln2; i<ln1-1; i++) line[i]=' ';
  line[ln1-1]='\n';
  rewind(fp);
  fwrite(line,sizeof(char),ln1,fp);
  fclose(fp);
  return 0;
}
