
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2001 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/****************************************************************************
*
* mklj makes an binary lj-parameterfile
*
******************************************************************************/

/* Quick and dirty -- needs lots of work */

#include <stdio.h>
#include <math.h>
#define SQR(a) ((a)*(a))


int main(void)

{ FILE *out;

  float r_0,r_end,r_step,r,r2,p1,p2,p3;
  float energie[2][2];
  float sigma[2][2];
  float p,q;
  int n_step;
  char fname[30];

/* Benutzer abfragen, erstmal Dateinamen */

  printf("\nAusgabedatei: "); scanf("%s",fname );
  printf("\nRadius Anf. : "); scanf("%f",&r_0 );
  printf("\nRadius Ende : "); scanf("%f",&r_end);
  printf("\nZahl Schritte: "); scanf("%d",&n_step);
  printf("\nExponent   p: "); scanf("%f",&p);
  printf("\nExponent   q: "); scanf("%f",&q);
  printf("\nEnergie11   : "); scanf("%f",&energie[0][0]);
  printf("\nEnergie22   : "); scanf("%f",&energie[1][1]);
  printf("\nEnergie12   : "); scanf("%f",&energie[0][1]);
  printf("\nSigma11     : "); scanf("%f",&sigma[0][0]);
  printf("\nSigma22     : "); scanf("%f",&sigma[1][1]);
  printf("\nSigma12     : "); scanf("%f",&sigma[0][1]);

  r_step = (SQR(r_end) - SQR(r_0)) / n_step;

/* Ausgabedatei oeffnen */
  out = fopen(fname,"w");
  if (NULL == out) {
    printf("Can't open %s for output.\n",fname);
    exit(1);
  };
  r2     = SQR(r_0);
  do { r = sqrt(r2);
       p1 = energie[0][0] * (pow((sigma[0][0]/r),p) - p/q * pow((sigma[0][0]/r),q));
       p2 = energie[1][1] * (pow((sigma[1][1]/r),p) - p/q * pow((sigma[1][1]/r),q));
       p3 = energie[0][1] * (pow((sigma[0][1]/r),p) - p/q * pow((sigma[0][1]/r),q));
       fprintf(out,"%f %f %f %f %f\n",r2,p1,p3,p3,p2);
       r2 += r_step;
  } while (r2<SQR(r_end));
  fclose(out);

  exit(0);

}






