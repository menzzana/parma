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
void setInitialCombination(int markercombo[MAX_MARKER_COMBINATIONS], int markerno,
                           int nmarkers, int combinations) {
  int i1;

  markercombo[0]=markerno;
  if (nmarkers>1)
    for (i1=1; i1<combinations; i1++)
      markercombo[i1]=markerno==0?1:0;
  }
//---------------------------------------------------------------------------
bool increaseCombination(int markercombo[MAX_MARKER_COMBINATIONS], int markerno,
                         int nmarkers, int combinations) {
  int i1;

  for (i1=combinations-1; i1>0; i1++) {
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
void swap(unsigned char a, unsigned char b) {
  a^=b;
  b^=a;
  a^=b;
  }
//---------------------------------------------------------------------------
void randomShuffle(unsigned char *data, int length) {
  int i1;

  for (i1=0; i1<length; i1++)
    swap(data[i1],data[(int)ran1()*length]);
  }
//---------------------------------------------------------------------------
void clearMDRResults(int mdrpartres[N_MDR_PARTS][PHENOTYPE_COMBINATIONS][LIST_ALLELE_MARKER_COMBINATIONS],
                     int combinations) {
  int i1,i2;

  i2=pow(combinations,ALLELE_COMBINATIONS);
  for (i1=0; i1<N_MDR_PARTS; i1++) {
    memset(mdrpartres[i1][CONTROL],0,i2);
    memset(mdrpartres[i1][CASE],0,i2);
    }
  }
//---------------------------------------------------------------------------
float balancedAccuracy(struct MDRData res) {
  float sens,spec;

  sens=res.tp/(res.tp+res.fn);
  spec=res.tn/(res.fp+res.tn);
  return (sens+spec)/2;
  }
//---------------------------------------------------------------------------
struct MDRAccuracy analyseAlleles(unsigned char **gendata, int markercombo[MAX_MARKER_COMBINATIONS],
                                 unsigned char *phenotype,unsigned char *mdrparts, int ncombo, int nindividuals) {
  struct MDRAccuracy accres;
  struct MDRData test,train;
  int mdrres[N_MDR_PARTS][PHENOTYPE_COMBINATIONS][LIST_ALLELE_MARKER_COMBINATIONS];
  int mdrsumres[PHENOTYPE_COMBINATIONS][LIST_ALLELE_MARKER_COMBINATIONS];
  int idxmark,idxind,idxres,idxparts;

  clearMDRResults(mdrres,ncombo);
  for (idxind=0; idxind<nindividuals; idxind++) {
    idxres=0;
    for(idxmark=0; idxmark<ncombo; idxmark++)
      idxres+=pow(ALLELE_COMBINATIONS,idxmark)*gendata[idxind][markercombo[idxmark]];
    mdrres[mdrparts[idxind]][phenotype[idxind]][idxres]++;
    mdrsumres[phenotype[idxind]][idxres]++;
    }
  test.tp=test.tn=test.fp=test.fn=0;
  train.tp=train.tn=train.fp=train.fn=0;
  for (idxparts=0; idxparts<N_MDR_PARTS; idxparts++)
    for(idxres=0; idxres<pow(ncombo,ALLELE_COMBINATIONS); idxres++) {
      if (mdrres[idxparts][CONTROL][idxres]==0 && mdrres[idxparts][CASE][idxres]==0)
        continue;
      if ((mdrres[CONTROL][idxres]-mdrres[idxparts][CONTROL][idxres])>
          (mdrres[CASE][idxres]-mdrres[idxparts][CASE][idxres])) {
        train.tn+=mdrres[idxparts][CONTROL][idxres];
        train.fp+=mdrres[idxparts][CASE][idxres];
        if (mdrres[idxparts][CONTROL][idxres]>mdrres[idxparts][CASE][idxres]) {
          test.fp+=mdrres[idxparts][CASE][idxres];
          test.tn+=mdrres[idxparts][CONTROL][idxres];
          }
        else {
          test.tp+=mdrres[idxparts][CASE][idxres];
          test.fn+=mdrres[idxparts][CONTROL][idxres];
          }
        }
      else {
        train.fn+=mdrres[idxparts][CONTROL][idxres];
        train.tp+=mdrres[idxparts][CASE][idxres];
        if (mdrres[idxparts][CONTROL][idxres]>mdrres[idxparts][CASE][idxres]) {
          test.fp+=mdrres[idxparts][CASE][idxres];
          test.tn+=mdrres[idxparts][CONTROL][idxres];
          }
        else {
          test.tp+=mdrres[idxparts][CASE][idxres];
          test.fn+=mdrres[idxparts][CONTROL][idxres];
          }
        }
      }
  memcpy(accres.markercombo,markercombo,ncombo);
  accres.pvaluetrain=accres.pvaluetest=0;
  accres.train=balancedAccuracy(train);
  accres.test=balancedAccuracy(test);
  return accres;
  }
//---------------------------------------------------------------------------
void mdr(unsigned char **gendata, unsigned char *phenotype, int *marker,
         int nmarkers, int nindividuals, int nselectmarkers, int permutations) {
  int ncombo,selmarkno,i1;
  struct MDRAccuracy mdracc,mdrpermacc;
  int markercombo[MAX_MARKER_COMBINATIONS];
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
  for (ncombo=1; ncombo<=MAX_MARKER_COMBINATIONS; ncombo++) {
    for (selmarkno=0; selmarkno<nselectmarkers; selmarkno++) {
      setInitialCombination(markercombo,marker[selmarkno],nmarkers,ncombo);
      do {
        mdracc=analyseAlleles(gendata,markercombo,phenotype,mdrparts,ncombo,nindividuals);
        //Permutations should be added here
        printf("Markers:");
        for (i1=0; i1<ncombo; i1++)
          printf(" %d",mdracc.markercombo[i1]);
        printf("\tTrain: %f\tTest: %f\n",mdracc.train,mdracc.test);
        /*
        for (i1=0; i1<permutation; i1++) {
          mdrpermacc=analyseAlleles(gendata,markercombo,permpheno[i1],mdrparts,ncombo,nindividuals);
          if (mdrpermacc.train>mdracc.train)
            mdracc.pvaluetrain++;
          if (mdrpermacc.test>mdracc.test)
            mdracc.pvaluetest++;
          }

        */
        } while (increaseCombination(markercombo,marker[selmarkno],nmarkers,ncombo));
      }
    }
  for (i1=0; i1<permutations; i1++)
    free(permpheno[i1]);
  free(permpheno);
  free(mdrparts);
  }
//---------------------------------------------------------------------------

