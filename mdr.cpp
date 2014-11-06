#include "mdr.h"
//---------------------------------------------------------------------------
using namespace MDR;
//---------------------------------------------------------------------------
SummedData::SummedData() {
  tp=tn=fp=fn=0;
  accuracy=nnegpermutations=0;
  }
//---------------------------------------------------------------------------
void SummedData::copy(SummedData summeddata) {
  tp=summeddata.tp;
  tn=summeddata.tn;
  fp=summeddata.fp;
  fn=summeddata.fn;
  accuracy=summeddata.accuracy;
  nnegpermutations=summeddata.nnegpermutations;
  }
//---------------------------------------------------------------------------
void SummedData::setAccuracy() {
  float sens,spec;

  sens=(tp+fn)==0?0:tp/(tp+fn);
  spec=(fp+tn)==0?0:tn/(fp+tn);
  accuracy=(sens+spec)/2;
  }
//---------------------------------------------------------------------------
double SummedData::getPvaluePermutations(int npermutations) {
  return (npermutations-nnegpermutations)/(double)(npermutations==0?1:npermutations);
  }
//---------------------------------------------------------------------------
Result::Result() {
  combinations=0;
  train=SummedData();
  test=SummedData();
  }
//---------------------------------------------------------------------------
void Result::copy(Result result) {
  combinations=result.combinations;
  memcpy(markercombo,result.markercombo,sizeof(int)*combinations);
  train.copy(result.train);
  test.copy(result.test);
  }
//---------------------------------------------------------------------------
bool Result::testBestCombination(double nnegpermutations1, double nnegpermutations2,
                                 double accuracy1, double accuracy2) {
  if (nnegpermutations1!=nnegpermutations2)
    return nnegpermutations1<nnegpermutations2;
  return accuracy1<accuracy2;
  }
//---------------------------------------------------------------------------
void Result::testBestCombination(Result result) {
  if (testBestCombination(test.nnegpermutations,result.test.nnegpermutations,
                          test.accuracy,result.test.accuracy))
    copy(result);
  }
//---------------------------------------------------------------------------
void Result::printHeader(bool ispermutation) {
  cout << "Markers"<< delimiter << "Train"<< delimiter << "Test"<< delimiter;
  if (ispermutation)
    cout << "Train p-value"<< delimiter << "Test p-value"<< delimiter;
  cout << "Best value" << endl;
  }
//---------------------------------------------------------------------------
void Result::print(char **marker,int npermutations, bool highest) {
  for (int i1=0; i1<combinations; i1++)
    cout << (i1==0?"":",")<<marker[markercombo[i1]];
  cout  << delimiter << train.accuracy << delimiter << test.accuracy;
  if (npermutations>0) {
    cout << delimiter << train.getPvaluePermutations(npermutations);
    cout << delimiter << test.getPvaluePermutations(npermutations);
    }
  if (highest)
    cout << delimiter << "*";
  cout << endl;
  }
//---------------------------------------------------------------------------
Analysis::Analysis() {
  npermutations=nindividuals=nmarkers=0;
  cutpvalue=NO_CUTOFF;
  gendata=NULL;
  phenotype=NULL;
  marker=NULL;
  permpheno=NULL;
  parts=NULL;
  }
//---------------------------------------------------------------------------
void Analysis::createDataBuffers(bool initthisrank) {
  if (!initthisrank)
    return;
  phenotype=new unsigned char[nindividuals];
  gendata=new unsigned char*[nindividuals];
  gendata[0]=new unsigned char[nindividuals*nmarkers];
  for (int i1=1; i1<nindividuals; i1++)
    gendata[i1]=&gendata[0][i1*nmarkers];
  marker=new char*[nmarkers];
  marker[0]=new char[nmarkers*global::MAX_LENGTH_MARKER_NAME];
  for (int i1=1; i1<nmarkers; i1++)
    marker[i1]=&marker[0][i1*global::MAX_LENGTH_MARKER_NAME];
  }
//---------------------------------------------------------------------------
void Analysis::setInitialArrays() {
  parts=new unsigned char[nindividuals];
  populateMDRParts();
  randomShuffle(parts);
  if (npermutations==0)
    return;
  permpheno=new unsigned char*[npermutations];
  for (int i1=0; i1<npermutations; i1++) {
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
    swap(data[i1],data[(int)(RND::ran1()*nindividuals)]);
  }
//---------------------------------------------------------------------------
int Analysis::getAlleleCombinations(int combinations) {
  return pow(ALLELE_COMBINATIONS,combinations);
  }
//---------------------------------------------------------------------------
bool Analysis::setInitialCombination(int idxmark, int combinations) {
  if (nmarkers<(combinations+idxmark))
    return false;
  markercombo[0]=idxmark;
  for (int i1=1; i1<combinations; i1++)
    markercombo[i1]=idxmark+combinations-i1;
  return true;
  }
//---------------------------------------------------------------------------
bool Analysis::increaseCombination(int idxmarkcombo, int combinations) {
  if (idxmarkcombo==combinations)
    return false;
  if (increaseMarker(idxmarkcombo))
    return true;
  if (!increaseCombination(idxmarkcombo+1,combinations))
    return false;
  markercombo[idxmarkcombo]=markercombo[idxmarkcombo+1];
  return increaseMarker(idxmarkcombo);
  }
//---------------------------------------------------------------------------
bool Analysis::increaseMarker(int idxmarkcombo) {
  markercombo[idxmarkcombo]+=(markercombo[0]==(markercombo[idxmarkcombo]+1)?2:1);
  return markercombo[idxmarkcombo]<nmarkers;
  }
//---------------------------------------------------------------------------
void Analysis::clearMDRResults(int combinations) {
  int vlength;

  vlength=getAlleleCombinations(combinations);
  memset(mdrsumres[CONTROL],0,sizeof(int)*vlength);
  memset(mdrsumres[CASE],0,sizeof(int)*vlength);
  for (int i1=0; i1<N_MDR_PARTS; i1++) {
    memset(mdrpartres[i1][CONTROL],0,sizeof(int)*vlength);
    memset(mdrpartres[i1][CASE],0,sizeof(int)*vlength);
    }
  }
//---------------------------------------------------------------------------
Result Analysis::analyseAlleles(unsigned char *vpheno, int combinations) {
  Result accres;
  int idxmark,idxind,idxres,idxparts;

  clearMDRResults(combinations);
  for (idxind=0; idxind<nindividuals; idxind++) {
    idxres=0;
    for(idxmark=0; idxmark<combinations; idxmark++)
      idxres+=getAlleleCombinations(idxmark)*gendata[idxind][markercombo[idxmark]];
    mdrpartres[(int)parts[idxind]][(int)vpheno[idxind]][idxres]++;
    mdrsumres[(int)vpheno[idxind]][idxres]++;
    }
  accres=Result();
  for (idxparts=0; idxparts<N_MDR_PARTS; idxparts++)
    for(idxres=0; idxres<getAlleleCombinations(combinations); idxres++) {
      if (mdrpartres[idxparts][CONTROL][idxres]==0 && mdrpartres[idxparts][CASE][idxres]==0)
        continue;
      bool controlratio;
      controlratio=mdrpartres[idxparts][CONTROL][idxres]>mdrpartres[idxparts][CASE][idxres];
      (controlratio?accres.train.fp:accres.train.tp)+=mdrpartres[idxparts][CASE][idxres];
      (controlratio?accres.train.tn:accres.train.fn)+=mdrpartres[idxparts][CONTROL][idxres];
      controlratio=(mdrsumres[CONTROL][idxres]-mdrpartres[idxparts][CONTROL][idxres])>
        (mdrsumres[CASE][idxres]-mdrpartres[idxparts][CASE][idxres]) &&
        mdrpartres[idxparts][CONTROL][idxres]>mdrpartres[idxparts][CASE][idxres];
      (controlratio?accres.test.fp:accres.test.tp)+=mdrpartres[idxparts][CASE][idxres];
      (controlratio?accres.test.tn:accres.test.fn)+=mdrpartres[idxparts][CONTROL][idxres];
      }
  memcpy(accres.markercombo,markercombo,sizeof(int)*combinations);
  accres.combinations=combinations;
  accres.train.setAccuracy();
  accres.test.setAccuracy();
  return accres;
  }
//---------------------------------------------------------------------------
bool Analysis::Run(int rank, int blocksize, int combination) {
  int idxmark;
  Result origaccuracy,permaccuracy;

  try {
    maxaccuracy=Result();
    for (idxmark=rank; idxmark<nmarkers; idxmark+=blocksize) {
      if (!setInitialCombination(idxmark,combination))
        continue;
      do {
        origaccuracy=analyseAlleles(phenotype,combination);
        for (int i1=0; i1<npermutations; i1++) {
          permaccuracy=analyseAlleles(permpheno[i1],combination);
          if (permaccuracy.train.accuracy<origaccuracy.train.accuracy)
            origaccuracy.train.nnegpermutations++;
          if (permaccuracy.test.accuracy<origaccuracy.test.accuracy)
            origaccuracy.test.nnegpermutations++;
          }
        if (origaccuracy.test.getPvaluePermutations(npermutations)<=cutpvalue)
          origaccuracy.print(marker,npermutations,false);
        maxaccuracy.testBestCombination(origaccuracy);
        } while (increaseCombination(1,combination));
      }

    return true;
    }
  catch(exception &e) {
    cerr << e.what() << endl;
    return false;
    }
  }
//---------------------------------------------------------------------------
void Analysis::printBestResult() {
  maxaccuracy.print(marker,npermutations,true);
  }
//---------------------------------------------------------------------------
Analysis::~Analysis() {
  for (int i1=0; i1<npermutations; i1++)
    delete permpheno[i1];
  delete permpheno;
  delete parts;
  delete phenotype;
  delete[] marker;
  delete[] gendata;
  }
//---------------------------------------------------------------------------
