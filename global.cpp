#include "global.h"
//---------------------------------------------------------------------------
void CALC::sran1(long seedvalue) {
  // Initialize with negative number
  if (rseed<0)
    rseed=seedvalue;
  }
//---------------------------------------------------------------------------
double CALC::ran1() {
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  double temp;

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
unsigned long long CALC::C(unsigned long long n, unsigned long long k) {
  unsigned long long d,r=1;

  if (k>n)
    return 0;
  for (d=1; d<=k; d++) {
    r*=n--;
    r/=d;
    }
  return r;
  }
//---------------------------------------------------------------------------
