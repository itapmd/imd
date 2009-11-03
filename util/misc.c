
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2006 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/

/******************************************************************************
* $Revision$
* $Date$
******************************************************************************/

/******************************************************************************
*
*  error -- Complain and abort
*
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

void error(char *msg)
{
  printf("Error: %s\n",msg);
  exit(2);
}

/* error message built from two strings */
void error_str(char *msg, char *str)
{
  char buf[255];
  sprintf(buf, msg, str);
  error(buf);
}

/* error message built from three strings */
void error_str_str(char *msg, char *str1, char *str2)
{
  char buf[255];
  sprintf(buf, msg, str1, str2);
  error(buf);
}

