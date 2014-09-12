#include "global.h"
#include "mdr.h"
//---------------------------------------------------------------------------
void init() {
  rseed=12345678;
  }
//---------------------------------------------------------------------------
void sran1(long value) {
  rseed=value;
  }
//---------------------------------------------------------------------------
double ran1() {
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
void initialMask(unsigned char *mask, int nanalysedmarkers, int nmarkers) {
  int i1;

  for (i1=0; i1<nmarkers; i1++)
    mask[i1]=i1<nanalysedmarkers?1:0;
  }
//---------------------------------------------------------------------------
bool incMask(unsigned char *mask, int nmarkers) {
  int i1,i2,nset;

  for (i1=nmarkers-1,nset=0; i1>=0; i1--)
    if (mask[i1]==1) {
      nset++;
      if (mask[i1+1]==1 || (i1+1)==nmarkers)
        continue;
      if ((nset+i1)>=nmarkers)
        return false;
      for (i2=i1; i2<nmarkers; i2++)
        mask[i2]=i2<(nset+i1+1) && i2!=i1?1:0;
      return true;
      }
  return false;
  }
//---------------------------------------------------------------------------
void andMask(char *oldmask,char *mask, int nmarkers) {
  int i1;

  for (i1=0; i1<nmarkers; i1++)
    oldmask[i1]=oldmask[i1]==mask[i1]?oldmask[i1]:0;
  }
//---------------------------------------------------------------------------
bool TestMask(char *oldmask,char *mask, int nmarkers) {
  int i1;

  for (i1=0; i1<nmarkers; i1++)
    if (oldmask[i1]==1 && mask[i1]==0)
      return false;
  return true;
  }
//---------------------------------------------------------------------------
void randomshuffle(unsigned char *data, int length) {
  int i1,i2;
  unsigned char c1;

  for (i1=0; i1<length; i1++) {
    i2=(int)ran1()*length;
    c1=data[i1];
    data[i1]=data[i2];
    data[i2]=c1;
    }
  }
//---------------------------------------------------------------------------
void mdr(unsigned char **data, unsigned char *phenotype, int nindividuals, int nmarkers) {
  unsigned char *mdrparts;
  int i1,i2;

  mdrparts=malloc(nindividuals);
  for (i1=i2=0; i1<nindividuals; i1++,i2++) {
    if (i2==N_MDR_PARTS)
      i2=0;
    mdrparts[i1]=i2;
    }
  randomshuffle(mdrparts,nindividuals);



  free(mdrparts);
  }
//---------------------------------------------------------------------------

