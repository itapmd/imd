#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include "ImgExaminerViewer.h"

typedef struct
{
  ImgExaminerViewer *viewer;
  int soc;
  unsigned char Camera_Type;
} 
   GlobalInfoTyp;

#endif
