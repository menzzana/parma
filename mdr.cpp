#include "global.h"
#include "mdr.h"
//---------------------------------------------------------------------------
using namespace MDR;
//---------------------------------------------------------------------------
Analysis::Analysis() {
  permutations=nindividuals=nmarkers=frommarker=tomarker=0;
  gendata=NULL;
  phenotype=NULL;
  marker=NULL;
  }
//---------------------------------------------------------------------------
void Analysis::setInitialArrays() {
  parts=new unsigned char[nindividuals];
  populateMDRParts();
  randomShuffle(parts);
  permpheno=new unsigned char*[permutations];
  for (int i1=0; i1<permutations; i1++) {
    permpheno[i1]=new unsigned char[nindividuals];
    memcpy(permpheno[i1],phenotype,nindividuals);
    randomShuffle(permpheno[i1]);
    }
  }
//---------------------------------------------------------------------------
void Analysis::populateMDRParts() {
  int i1,i2;

  for (i1=i2=0; i1<nindividuals; i1++,i2++) {
    if (i2==N_MDR_PARTS)
      i2=0;
    parts[i1]=i2;
    }
  }
//---------------------------------------------------------------------------
void Analysis::randomShuffle(unsigned char *data) {
  for (int i1=0; i1<nindividuals; i1++)
    std::swap(data[i1],data[(int)ran1(0)*nindividuals]);
  }
//---------------------------------------------------------------------------
void Analysis::setInitialCombination(int markercombo[MAX_MARKER_COMBINATIONS], int idxmark, int combinations) {
  markercombo[0]=idxmark;
  if (nmarkers>1)
    for (int i1=1; i1<combinations; i1++)
      markercombo[i1]=idxmark==0?1:0;
  }
//---------------------------------------------------------------------------
bool Analysis::increaseCombination(int markercombo[MAX_MARKER_COMBINATIONS], int idxmark, int combinations) {
  for (int i1=combinations-1; i1>0; i1++) {
    markercombo[i1]++;
    if (markercombo[i1]==idxmark)
      markercombo[i1]++;
    if (markercombo[i1]<nmarkers)
      return false;
    markercombo[i1]=idxmark==0?1:0;
    }
  return true;
  }
//---------------------------------------------------------------------------
void Analysis::clearMDRResults(int mdrpartres[N_MDR_PARTS][PHENOTYPE_COMBINATIONS][LIST_ALLELE_MARKER_COMBINATIONS],
                     int mdrsumres[PHENOTYPE_COMBINATIONS][LIST_ALLELE_MARKER_COMBINATIONS],
                     int combinations) {
  int vlength;

  vlength=pow(combinations,ALLELE_COMBINATIONS);
  memset(mdrsumres[CONTROL],0,vlength);
  memset(mdrsumres[CASE],0,vlength);
  for (int i1=0; i1<N_MDR_PARTS; i1++) {
    memset(mdrpartres[i1][CONTROL],0,vlength);
    memset(mdrpartres[i1][CASE],0,vlength);
    }
  }
//---------------------------------------------------------------------------
void SummedData::clear() {
  tp=tn=fp=fn=0;
  accuracy=npospermutations=0;
  }
//---------------------------------------------------------------------------
void SummedData::setAccuracy() {
  float sens,spec;

  sens=tp/(tp+fn);
  spec=tn/(fp+tn);
  accuracy=(sens+spec)/2;
  }
//---------------------------------------------------------------------------
Result Analysis::analyseAlleles(int markercombo[MAX_MARKER_COMBINATIONS],unsigned char *vpheno,int ncombo) {
  Result accres;
  int mdrres[N_MDR_PARTS][PHENOTYPE_COMBINATIONS][LIST_ALLELE_MARKER_COMBINATIONS];
  int mdrsumres[PHENOTYPE_COMBINATIONS][LIST_ALLELE_MARKER_COMBINATIONS];
  int idxmark,idxind,idxres,idxparts;

  clearMDRResults(mdrres,mdrsumres,ncombo);
  for (idxind=0; idxind<nindividuals; idxind++) {
    idxres=0;
    for(idxmark=0; idxmark<ncombo; idxmark++)
      idxres+=pow(ALLELE_COMBINATIONS,idxmark)*gendata[idxind][markercombo[idxmark]];
    mdrres[parts[idxind]][vpheno[idxind]][idxres]++;
    mdrsumres[vpheno[idxind]][idxres]++;
    }
  accres.test.clear();
  accres.train.clear();
  for (idxparts=0; idxparts<N_MDR_PARTS; idxparts++)
    for(idxres=0; idxres<pow(ncombo,ALLELE_COMBINATIONS); idxres++) {
      if (mdrres[idxparts][CONTROL][idxres]==0 && mdrres[idxparts][CASE][idxres]==0)
        continue;
      if ((mdrres[CONTROL][idxres]-mdrres[idxparts][CONTROL][idxres])>
          (mdrres[CASE][idxres]-mdrres[idxparts][CASE][idxres])) {
        accres.train.tn+=mdrres[idxparts][CONTROL][idxres];
        accres.train.fp+=mdrres[idxparts][CASE][idxres];
        if (mdrres[idxparts][CONTROL][idxres]>mdrres[idxparts][CASE][idxres]) {
          accres.test.fp+=mdrres[idxparts][CASE][idxres];
          accres.test.tn+=mdrres[idxparts][CONTROL][idxres];
          }
        else {
          accres.test.tp+=mdrres[idxparts][CASE][idxres];
          accres.test.fn+=mdrres[idxparts][CONTROL][idxres];
          }
        }
      else {
        accres.train.fn+=mdrres[idxparts][CONTROL][idxres];
        accres.train.tp+=mdrres[idxparts][CASE][idxres];
        if (mdrres[idxparts][CONTROL][idxres]>mdrres[idxparts][CASE][idxres]) {
          accres.test.fp+=mdrres[idxparts][CASE][idxres];
          accres.test.tn+=mdrres[idxparts][CONTROL][idxres];
          }
        else {
          accres.test.tp+=mdrres[idxparts][CASE][idxres];
          accres.test.fn+=mdrres[idxparts][CONTROL][idxres];
          }
        }
      }
  memcpy(accres.markercombo,markercombo,ncombo);
  accres.train.setAccuracy();
  accres.test.setAccuracy();
  return accres;
  }
//---------------------------------------------------------------------------
void Analysis::Run() {
  int ncombo,idxmark;
  Result origaccuracy,permaccuracy;
  int markercombo[MAX_MARKER_COMBINATIONS];

  for (ncombo=1; ncombo<=MAX_MARKER_COMBINATIONS; ncombo++) {
    for (idxmark=frommarker; idxmark<tomarker; idxmark++) {
      setInitialCombination(markercombo,idxmark,ncombo);
      do {
        origaccuracy=analyseAlleles(markercombo,phenotype,ncombo);
        //Permutations should be added here
        printf("Markers:");
        for (int i1=0; i1<ncombo; i1++)
          printf(" %d",marker[origaccuracy.markercombo[i1]]);
        printf("\tTrain: %f\tTest: %f\n",origaccuracy.train.accuracy,origaccuracy.test.accuracy);
        /*
        for (i1=0; i1<permutation; i1++) {
          permaccuracy=analyseAlleles(markercombo,permpheno[i1],ncombo);
          if (permaccuracy.train>origaccuracy.train)
            origaccuracy.train.npospermutations++;
          if (permaccuracy.test>origaccuracy.test)
            origaccuracy.test.npospermutations++;
          }

        */
        } while (increaseCombination(markercombo,idxmark,ncombo));
      }
    }
  }
//---------------------------------------------------------------------------
Analysis::~Analysis() {
  for (int i1=0; i1<permutations; i1++)
    delete permpheno[i1];
  delete permpheno;
  delete parts;
  }
//---------------------------------------------------------------------------

