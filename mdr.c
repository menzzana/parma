#include "global.h"
#include "mdr.h"
//---------------------------------------------------------------------------
void init() {
  rseed=-12345678;
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
void setInitialCombination(int *markercombo, int markerno, int nmarkers, int ncombo) {
  int i1;

  markercombo[0]=markerno;
  if (nmarkers>1)
    for (i1=1; i1<ncombo; i1++)
      markercombo[i1]=markerno==0?1:0;
  }
//---------------------------------------------------------------------------
bool increaseCombination(int *markercombo, int markerno, int nmarkers, int ncombo) {
  int i1;

  for (i1=ncombo-1; i1>0; i1++) {
    markercombo[i1]++;
    if (markercombo[i1]==markerno)
      markercombo[i1]++;
    if (markercombo[i1]<nmarkers)
      return false;
    markercombo[i1]=markerno==0?1:0;
    }
  return true;
  }
//---------------------------------------------------------------------------
void populateMDRParts(unsigned char *data, int length) {
  int i1,i2;

  for (i1=i2=0; i1<length; i1++,i2++) {
    if (i2==N_MDR_PARTS)
      i2=0;
    data[i1]=i2;
    }
  }
//---------------------------------------------------------------------------
void randomShuffle(unsigned char *data, int length) {
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
void mdr(unsigned char **gendata, unsigned char *phenotype, int *marker,
         int nmarkers, int nindividuals, int nselectmarkers, int permutations) {
  int ncombo,i1;
  int *markercombo;
  unsigned char **permpheno, *mdrparts;

  mdrparts=malloc(nindividuals);
  populateMDRParts(mdrparts,nindividuals);
  randomShuffle(mdrparts,nindividuals);
  permpheno=malloc(permutations);
  for (i1=0; i1<permutations; i1++) {
    permpheno[i1]=malloc(nindividuals);
    memcpy(permpheno[i1],phenotype,nindividuals);
    randomShuffle(permpheno[i1],nindividuals);
    }
  markercombo=malloc(MAX_MARKER_COMBINATIONS);
  for (ncombo=1; ; ncombo++) {
    for (i1=0; i1<nselectmarkers; i1++) {
      setInitialCombination(markercombo,marker[i1],nmarkers,ncombo);
      do {

        } while (increaseCombination(markercombo,marker[i1],nmarkers,ncombo));

      }
    }
  for (i1=0; i1<permutations; i1++)
    free(permpheno[i1]);
  free(permpheno);
  free(markercombo);
  free(mdrparts);
  }
//---------------------------------------------------------------------------

