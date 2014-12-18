#include "mdr.h"
//---------------------------------------------------------------------------
using namespace MDR;
//---------------------------------------------------------------------------
SummedData::SummedData() {
  calc.accuracy=calc.nnegpermutations=0;
  tp=fp=tn=fn=0;
  memset(partaccuracy,0,sizeof(double)*N_MDR_PARTS);
  }
//---------------------------------------------------------------------------
void SummedData::clearPartData() {
  tp=fp=tn=fn=0;
  }
//---------------------------------------------------------------------------
void SummedData::addAccuracy(int idxpart) {
  float sens,spec;

  sens=(tp+fn)==0?0:tp/(tp+fn);
  spec=(fp+tn)==0?0:tn/(fp+tn);
  partaccuracy[idxpart]=(sens+spec)/2;
  calc.accuracy+=partaccuracy[idxpart];
  }
//---------------------------------------------------------------------------
double SummedData::getPvaluePermutations(int npermutations) {
  return (N_MDR_PARTS*npermutations-calc.nnegpermutations)/
      (double)(npermutations==0?1:npermutations)/
      (double)N_MDR_PARTS;
  }
//---------------------------------------------------------------------------
bool SummedData::testBestCombination(Calculated calc1, Calculated calc2) {
  if (calc1.nnegpermutations!=calc2.nnegpermutations)
    return calc1.nnegpermutations<calc2.nnegpermutations;
  return calc1.accuracy<calc2.accuracy;
  }
//---------------------------------------------------------------------------
#ifndef SERIAL
  void SummedData::procTestBestCombination(Calculated *in, Calculated *inout, int *len, MPI_Datatype *type) {
    if (testBestCombination(*inout,*in)) {
      inout->nnegpermutations=in->nnegpermutations;
      inout->accuracy=in->accuracy;
      inout->rank=in->rank;
      }
    }
#endif
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
  memcpy(&train,&result.train,sizeof(SummedData));
  memcpy(&test,&result.test,sizeof(SummedData));
  }
//---------------------------------------------------------------------------
void Result::setBestCombination(Result result) {
  if (test.testBestCombination(test.calc,result.test.calc))
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
  for (int i1=0; i1<combinations; i1++) {
    cout <<marker[markercombo[i1]];
    if (i1==0) {
      cout  << delimiter << train.calc.accuracy << delimiter << test.calc.accuracy;
      if (npermutations>0) {
        cout << delimiter << train.getPvaluePermutations(npermutations);
        cout << delimiter << test.getPvaluePermutations(npermutations);
        }
      if (highest)
        cout << delimiter << "*";
      }
    cout << endl;
    }
  }
//---------------------------------------------------------------------------
Analysis::Analysis() {
  cutpvalue=NO_CUTOFF;
  gendata=NULL;
  phenotype=NULL;
  marker=NULL;
  permpheno=NULL;
  parts=NULL;
  param.npermutations=0;
  param.maxcombinations=MAX_MARKER_COMBINATIONS;
  param.nmarkers=0;
  param.nindividuals=0;
  param.randomseed=0;
  param.cutpvalue=NO_CUTOFF;
  }
//---------------------------------------------------------------------------
void Analysis::createMasterDataBuffers(int nmarker, int nindividual) {
  param.nmarkers=nmarker;
  param.nindividuals=nindividual;
  createDataBuffers(true);
  }
//---------------------------------------------------------------------------
void Analysis::createDataBuffers(bool initthisrank) {
  if (!initthisrank)
    return;
  phenotype=new unsigned char[param.nindividuals];
  gendata=new unsigned char*[param.nindividuals];
  gendata[0]=new unsigned char[param.nindividuals*param.nmarkers];
  for (int i1=1; i1<param.nindividuals; i1++)
    gendata[i1]=&gendata[0][i1*param.nmarkers];
  memset(&gendata[0][0],255,param.nindividuals*param.nmarkers);
  marker=new char*[param.nmarkers];
  marker[0]=new char[param.nmarkers*global::MAX_LENGTH_MARKER_NAME];
  for (int i1=1; i1<param.nmarkers; i1++)
    marker[i1]=&marker[0][i1*global::MAX_LENGTH_MARKER_NAME];
  }
//---------------------------------------------------------------------------
void Analysis::setInitialArrays() {
  parts=new unsigned char[param.nindividuals];
  populateMDRParts();
  randomShuffle(parts);
  if (param.npermutations==0)
    return;
  permpheno=new unsigned char*[param.npermutations];
  for (int i1=0; i1<param.npermutations; i1++) {
    permpheno[i1]=new unsigned char[param.nindividuals];
    memcpy(permpheno[i1],phenotype,param.nindividuals);
    randomShuffle(permpheno[i1]);
    }
  }
//---------------------------------------------------------------------------
void Analysis::populateMDRParts() {
  int i1,i2;

  for (i1=i2=0; i1<param.nindividuals; i1++,i2++) {
    if (i2==N_MDR_PARTS)
      i2=0;
    parts[i1]=i2;
    }
  }
//---------------------------------------------------------------------------
void Analysis::randomShuffle(unsigned char *data) {
  for (int i1=0; i1<param.nindividuals; i1++)
    swap(data[i1],data[(int)(CALC::ran1()*param.nindividuals)]);
  }
//---------------------------------------------------------------------------
int Analysis::getAlleleCombinations(int combinations) {
  return pow(ALLELE_COMBINATIONS,combinations);
  }
//---------------------------------------------------------------------------
void Analysis::setMarkerCombination(unsigned long long cidx, int combinations) {
  int i1,i2;
  unsigned long long cidx1;

  cidx1=cidx;
  for (i1=0; i1<combinations; i1++) {
    for (i2=0; CALC::C(i2+1,combinations-i1)<=cidx1; i2++);
    cidx1-=CALC::C(i2,combinations-i1);
    markercombo[i1]=i2;
    }
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
void Analysis::removeNonGenotypeIndividuals() {
  int idxindividual,idxnewindividual,idxmarker;

  for (idxindividual=idxnewindividual=0; idxindividual<param.nindividuals; idxindividual++) {
    phenotype[idxnewindividual]=phenotype[idxindividual];
    for (idxmarker=0; idxmarker<param.nmarkers; idxmarker++) {
      if (gendata[idxindividual][idxmarker]==255)
        break;
      gendata[idxnewindividual][idxmarker]=gendata[idxindividual][idxmarker];
      }
    if (idxmarker==param.nmarkers)
      idxnewindividual++;
    }
  cerr << "Removed " << param.nindividuals-idxnewindividual << " individuals" << endl;
  param.nindividuals=idxnewindividual;
  }
//---------------------------------------------------------------------------
Result Analysis::analyseAlleles(unsigned char *vpheno, int combinations, Result *original) {
  Result accres;
  int idxmark,idxind,idxres,idxpart;
  int traincase,traincontrol;

  clearMDRResults(combinations);
  for (idxind=0; idxind<param.nindividuals; idxind++) {
    idxres=0;
    for(idxmark=0; idxmark<combinations; idxmark++)
      idxres+=getAlleleCombinations(idxmark)*gendata[idxind][markercombo[idxmark]];
    mdrpartres[(int)parts[idxind]][(int)vpheno[idxind]][idxres]++;
    mdrsumres[(int)vpheno[idxind]][idxres]++;
    }
  accres=Result();
  for (idxpart=0; idxpart<N_MDR_PARTS; idxpart++) {
    accres.train.clearPartData();
    accres.test.clearPartData();
    for(idxres=0; idxres<getAlleleCombinations(combinations); idxres++) {
      traincase=mdrsumres[CASE][idxres]-mdrpartres[idxpart][CASE][idxres];
      traincontrol=mdrsumres[CONTROL][idxres]-mdrpartres[idxpart][CONTROL][idxres];
      if (traincase<traincontrol) {
        accres.train.fp+=traincase;
        accres.train.tn+=traincontrol;
        accres.test.fp+=mdrpartres[idxpart][CASE][idxres];
        accres.test.tn+=mdrpartres[idxpart][CONTROL][idxres];
        }
      else {
        accres.train.tp+=traincase;
        accres.train.fn+=traincontrol;
        accres.test.tp+=mdrpartres[idxpart][CASE][idxres];
        accres.test.fn+=mdrpartres[idxpart][CONTROL][idxres];
        }
      }
    accres.train.addAccuracy(idxpart);
    accres.test.addAccuracy(idxpart);
    if (original!=NULL) {
      if (accres.train.partaccuracy[idxpart]<original->train.partaccuracy[idxpart])
        original->train.calc.nnegpermutations++;
      if (accres.test.partaccuracy[idxpart]<original->test.partaccuracy[idxpart])
        original->test.calc.nnegpermutations++;
      }
    }
  memcpy(accres.markercombo,markercombo,sizeof(int)*combinations);
  accres.combinations=combinations;
  accres.train.calc.accuracy/=(double)N_MDR_PARTS;
  accres.test.calc.accuracy/=(double)N_MDR_PARTS;
  return accres;
  }
//---------------------------------------------------------------------------
bool Analysis::Run(int rank, int mpisize, int combination) {
  unsigned long long cidx,maxmarkercombos;
  Result origaccuracy;

  try {
    maxaccuracy=Result();
    maxmarkercombos=CALC::C(param.nmarkers,combination);
    for (cidx=rank; cidx<maxmarkercombos; cidx+=mpisize) {
      setMarkerCombination(cidx,combination);
      origaccuracy=analyseAlleles(phenotype,combination,NULL);
      for (int i1=0; i1<param.npermutations; i1++)
        analyseAlleles(permpheno[i1],combination,&origaccuracy);
      if (origaccuracy.test.getPvaluePermutations(param.npermutations)<=param.cutpvalue)
        origaccuracy.print(marker,param.npermutations,false);
      maxaccuracy.setBestCombination(origaccuracy);
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
  maxaccuracy.print(marker,param.npermutations,true);
  }
//---------------------------------------------------------------------------
void Analysis::printParameters() {
  clog << "Markers: " << param.nmarkers << endl;
  clog << "Individuals: " << param.nindividuals << endl;
  if (param.npermutations==0)
    clog << "Permutations: None" << endl;
  else
    clog << "Permutations: " << param.npermutations << endl;
  clog << "Max combinations: " << param.maxcombinations << endl;
  if (param.cutpvalue==NO_CUTOFF)
    clog << "Max p-value: Best value only" << endl;
  else
    clog << "Max p-value: " << param.cutpvalue << endl;
  if (param.randomseed==0)
    clog << "Seed: Std" << endl;
  else
    clog << "Seed: " << param.randomseed << endl;
  }
//---------------------------------------------------------------------------
Analysis::~Analysis() {
  for (int i1=0; i1<param.npermutations; i1++)
    delete permpheno[i1];
  delete permpheno;
  delete parts;
  delete phenotype;
  delete[] marker;
  delete[] gendata;
  }
//---------------------------------------------------------------------------
