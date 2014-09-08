#ifndef MDR_H
#define MDR_H
//------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
//------------------------------------------------------------------------------
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define MAX(x,y) (x>y?x:y)
#define MIN(x,y) (x<y?x:y)
//------------------------------------------------------------------------------
void sran1(long value);
double ran1();
//------------------------------------------------------------------------------
#endif // MDR_H
