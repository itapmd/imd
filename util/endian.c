/* Routine zum Testen der Hardware-Architektur */

#ifdef TEST
int endian(void);

char *lbtext[2] = {"little", "big"};

int main()
{
  printf("This machine is %s endian\n",lbtext[endian()]);
  exit(0);
}

#endif

/* Testen der Hardware-Architektur, evtl. fuer die Socket-Kommunikation */
/* Diese Routine liefert 0 fuer little endian und 1 fuer big endian zurueck */
int endian(void)
{
  unsigned short a;
  unsigned char *b;

  a = 0x4501;
  b = (unsigned char *) (&a);
  if (b[0]==1) return 0; /* little */
  else return 1;
}
