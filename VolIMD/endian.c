#ifdef TEST
int endian(void);

char *lbtext[2] = {"little", "big"};

main()
{
  printf("This machine is %s endian\n",lbtext[endian()]);
}

#endif


int endian(void)
{
  unsigned short a;
  unsigned char *b;

  a = 0x4501;
  b = (unsigned char *) (&a);
  if (b[0]==1) return 0; /* little */
  else return 1;
}
