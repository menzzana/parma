#include "global.h"
//---------------------------------------------------------------------------
double global::ran1(long option) {
  int j;
  long k;
  static long rseed=-123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (option!=0)
    rseed=option;
  if (rseed<=0 || !iy) {
    if (-rseed<1)
      rseed=1;
    else
      rseed=-rseed;
    for (j=NTAB+7; j>=0; j--) {
      k=rseed/IQ;
      rseed=IA*(rseed-k*IQ)-IR*k;
      if (rseed<0)
        rseed+=IM;
      if (j<NTAB)
        iv[j]=rseed;
      }
    iy=iv[0];
    }
  k=(rseed)/IQ;
  rseed=IA*(rseed-k*IQ)-IR*k;
  if (rseed<0)
    rseed+=IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j]=rseed;
  temp=AM*iy;
  if (temp>RNMX)
    return RNMX;
  return temp;
  }
//---------------------------------------------------------------------------
