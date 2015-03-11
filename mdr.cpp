#include "mdr.h"
//---------------------------------------------------------------------------
using namespace MDR;
//---------------------------------------------------------------------------
SummedData::SummedData() {
  calc.error=0;
  calc.pospermutations=0;
  tp=fp=tn=fn=0;
  parterror=0;
  originalparterror=0;
  }
//---------------------------------------------------------------------------
void SummedData::clearPartData() {
  tp=fp=tn=fn=0;
  parterror=0;
  }
//---------------------------------------------------------------------------
void SummedData::calculateError(int idxperm) {
  double caseerror,controlerror;

  caseerror=(tp+fp)==0?0:fp/(tp+fp);
  controlerror=(fn+tn)==0?0:fn/(fn+tn);
  parterror=(caseerror+controlerror)/2;
  if (idxperm==0) {
    calc.error+=parterror;
    originalparterror=parterror;
    }
  else
    if (parterror<=originalparterror)
      calc.pospermutations++;
  }
//---------------------------------------------------------------------------
double SummedData::getPvaluePermutations(int npermutations) {
  return calc.pospermutations/
      (double)(npermutations==0?1:npermutations)/
      (double)N_MDR_PARTS;
  }
//---------------------------------------------------------------------------
bool SummedData::testBestCombination(Calculated calc1, Calculated calc2) {
  if (calc1.pospermutations!=calc2.pospermutations)
    return calc1.pospermutations>calc2.pospermutations;
  return calc1.error>calc2.error;
  }
//---------------------------------------------------------------------------
#ifndef SERIAL
  void SummedData::procTestBestCombination(Calculated *in, Calculated *inout, int *len, MPI_Datatype *type) {
    if (testBestCombination(*inout,*in)) {
      inout->pospermutations=in->pospermutations;
      inout->error=in->error;
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
  for (int i1=0; i1<HEADER_LENGTH; i1++) {
    if (!ispermutation && (i1==3 || i1==4))
      continue;
    cout << HEADER[i1]<<delimiter;
    }
  cout << endl;
  }
//---------------------------------------------------------------------------
void Result::print(char **marker,int npermutations, bool lowest) {
  for (int i1=0; i1<combinations; i1++) {
    cout <<marker[markercombo[i1]];
    if (i1==0) {
      cout  << delimiter << train.calc.error << delimiter << test.calc.error;
      if (npermutations>0) {
        cout << delimiter << train.getPvaluePermutations(npermutations);
        cout << delimiter << test.getPvaluePermutations(npermutations);
        }
      if (lowest)
        cout << delimiter << "*";
      }
    cout << endl;
    }
  }
//---------------------------------------------------------------------------
Analysis::Analysis() {
  gendata=NULL;
  phenotype=NULL;
  marker=NULL;
  parts=NULL;
  allelegroup=NULL;
  param.onlypermuteone=false;
  param.npermutations=0;
  param.maxcombinations=MAX_MARKER_COMBINATIONS;
  param.mincombinations=MIN_MARKER_COMBINATIONS;
  param.nmarkers=0;
  param.nindividuals=0;
  param.randomseed=0;
  param.cutpvalue=NO_CUTOFF;
  }
//---------------------------------------------------------------------------
void Analysis::createDataBuffers() {
  phenotype=new unsigned char*[param.npermutations+1];
  phenotype[0]=new unsigned char[param.nindividuals*(param.npermutations+1)];
  for (int i1=1; i1<=param.npermutations; i1++)
    phenotype[i1]=&phenotype[0][i1*param.nindividuals];
  parts=new unsigned char[param.nindividuals];
  allelegroup=new short[param.nindividuals];
  mdrsumres[CONTROL]=new short*[param.npermutations+1];
  mdrsumres[CASE]=new short*[param.npermutations+1];
  mdrsumres[CONTROL][0]=new short[param.nindividuals*(param.npermutations+1)];
  memset(mdrsumres[CONTROL][0],0,param.nindividuals*(param.npermutations+1)*sizeof(short));
  mdrsumres[CASE][0]=new short[param.nindividuals*(param.npermutations+1)];
  memset(mdrsumres[CASE][0],0,param.nindividuals*(param.npermutations+1)*sizeof(short));
  for (int i1=1; i1<=param.npermutations; i1++) {
    mdrsumres[CONTROL][i1]=&mdrsumres[CONTROL][0][i1*param.nindividuals];
    mdrsumres[CASE][i1]=&mdrsumres[CASE][0][i1*param.nindividuals];
    }
  for (int i1=0; i1<N_MDR_PARTS; i1++) {
    mdrpartres[i1][CONTROL]=new short*[param.npermutations+1];
    mdrpartres[i1][CASE]=new short*[param.npermutations+1];
    mdrpartres[i1][CONTROL][0]=new short[param.nindividuals*(param.npermutations+1)];
    memset(mdrpartres[i1][CONTROL][0],0,param.nindividuals*(param.npermutations+1)*sizeof(short));
    mdrpartres[i1][CASE][0]=new short[param.nindividuals*(param.npermutations+1)];
    memset(mdrpartres[i1][CASE][0],0,param.nindividuals*(param.npermutations+1)*sizeof(short));
    for (int i2=1; i2<=param.npermutations; i2++) {
      mdrpartres[i1][CONTROL][i2]=&mdrpartres[i1][CONTROL][0][i2*param.nindividuals];
      mdrpartres[i1][CASE][i2]=&mdrpartres[i1][CASE][0][i2*param.nindividuals];
      }
    }
  gendata=new unsigned char*[param.nindividuals];
  gendata[0]=new unsigned char[param.nindividuals*param.nmarkers];
  for (int i1=1; i1<param.nindividuals; i1++)
    gendata[i1]=&gendata[0][i1*param.nmarkers];
  marker=new char*[param.nmarkers];
  marker[0]=new char[param.nmarkers*global::MAX_LENGTH_MARKER_NAME];
  for (int i1=1; i1<param.nmarkers; i1++)
    marker[i1]=&marker[0][i1*global::MAX_LENGTH_MARKER_NAME];
  }
//---------------------------------------------------------------------------
void Analysis::initializePartPermutationArrays() {
  populateMDRParts();
  randomShuffle(parts);
  for (int i1=1; i1<=param.npermutations; i1++) {
    memcpy(phenotype[i1],phenotype[0],param.nindividuals);
    randomShuffle(phenotype[i1]);
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
void Analysis::checkParameters() {
  param.maxcombinations=min(param.maxcombinations,param.nmarkers);
  param.onlypermuteone=param.onlypermuteone && param.cutpvalue==NO_CUTOFF;
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
void Analysis::clearMDRResults(int groupn) {
  int idxperm,idxpart;

  if (groupn>0)
    for (idxperm=0; idxperm<=param.npermutations; idxperm++) {
      memset(mdrsumres[CASE][idxperm],0,groupn*sizeof(short));
      memset(mdrsumres[CONTROL][idxperm],0,groupn*sizeof(short));
      for (idxpart=0; idxpart<N_MDR_PARTS; idxpart++) {
        memset(mdrpartres[idxpart][CASE][idxperm],0,groupn*sizeof(short));
        memset(mdrpartres[idxpart][CONTROL][idxperm],0,groupn*sizeof(short));
        }
      }
  }
//---------------------------------------------------------------------------
Result Analysis::analyseAlleles(int combinations,bool nopermutation) {
  Result accres;
  int idxind,idxgroup,idxmark,idxpart,idxperm;
  int traincase,traincontrol,permutations;
  static int groupn=0;

  permutations=nopermutation?0:param.npermutations;
  clearMDRResults(groupn);
  groupn=0;
  for (idxind=0; idxind<param.nindividuals; idxind++) {
    for (idxgroup=0; idxgroup<groupn; idxgroup++) {
      for (idxmark=0; idxmark<combinations; idxmark++)
        if (gendata[idxind][markercombo[idxmark]]!=gendata[allelegroup[idxgroup]][markercombo[idxmark]])
          break;
      if (idxmark==combinations)
       break;
      }
    if (idxgroup==groupn) {
      groupn++;
      allelegroup[idxgroup]=idxind;
      }
    for (idxperm=0; idxperm<=permutations; idxperm++) {
      mdrpartres[(int)parts[idxind]][(int)phenotype[idxperm][idxind]][idxperm][idxgroup]++;
      mdrsumres[(int)phenotype[idxperm][idxind]][idxperm][idxgroup]++;
      }
    }
  accres=Result();
  for (idxpart=0; idxpart<N_MDR_PARTS; idxpart++) {
    for (idxperm=0; idxperm<=permutations; idxperm++) {
      accres.train.clearPartData();
      accres.test.clearPartData();
      for(idxgroup=0; idxgroup<groupn; idxgroup++) {
        traincase=mdrsumres[CASE][idxperm][idxgroup]-mdrpartres[idxpart][CASE][idxperm][idxgroup];
        traincontrol=mdrsumres[CONTROL][idxperm][idxgroup]-mdrpartres[idxpart][CONTROL][idxperm][idxgroup];
        if (traincase<traincontrol) {
          accres.train.fp+=traincase;
          accres.train.tn+=traincontrol;
          accres.test.fp+=mdrpartres[idxpart][CASE][idxperm][idxgroup];
          accres.test.tn+=mdrpartres[idxpart][CONTROL][idxperm][idxgroup];
          }
        else {
          accres.train.tp+=traincase;
          accres.train.fn+=traincontrol;
          accres.test.tp+=mdrpartres[idxpart][CASE][idxperm][idxgroup];
          accres.test.fn+=mdrpartres[idxpart][CONTROL][idxperm][idxgroup];
          }
        }
      accres.train.calculateError(idxperm);
      accres.test.calculateError(idxperm);
      }
    }
  memcpy(accres.markercombo,markercombo,sizeof(int)*combinations);
  accres.combinations=combinations;
  accres.train.calc.error/=(double)N_MDR_PARTS;
  accres.test.calc.error/=(double)N_MDR_PARTS;
  return accres;
  }
//---------------------------------------------------------------------------
bool Analysis::Run(int rank, int mpisize, int combination) {
  unsigned long long cidx,maxmarkercombos;
  Result origerror;

  try {
    minerror=Result();
    minerror.test.calc.error=MAX_ERROR;
    minerror.test.calc.pospermutations=param.npermutations*N_MDR_PARTS;
    maxmarkercombos=CALC::C(param.nmarkers,combination);
    if (maxmarkercombos==0)
      THROW_ERROR("Max Combination overflow");
    for (cidx=rank; cidx<maxmarkercombos; cidx+=mpisize) {
      setMarkerCombination(cidx,combination);
      origerror=analyseAlleles(combination,param.onlypermuteone);
      if (origerror.test.getPvaluePermutations(param.npermutations)<=param.cutpvalue)
        origerror.print(marker,param.npermutations,false);
      minerror.setBestCombination(origerror);
      }
    if (param.onlypermuteone) {
      memcpy(markercombo,minerror.markercombo,sizeof(int)*combination);
      minerror=analyseAlleles(combination,false);
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
  minerror.print(marker,param.npermutations,true);
  }
//---------------------------------------------------------------------------
void Analysis::printParameters() {
  clog << PARAMETERS[0] << param.nmarkers << endl;
  clog << PARAMETERS[1] << param.nindividuals << endl;
  clog << PARAMETERS[2];
  if (param.onlypermuteone)
    clog << PARAMETERS[9];
  else
    clog << param.npermutations;
  clog << endl;
  clog << PARAMETERS[3] << param.mincombinations << endl;
  clog << PARAMETERS[4] << param.maxcombinations << endl;
  if (param.cutpvalue==NO_CUTOFF)
    clog << PARAMETERS[6] << endl;
  else
    clog << PARAMETERS[5] << param.cutpvalue << endl;
  if (param.randomseed==0)
    clog << PARAMETERS[8] << endl;
  else
    clog << PARAMETERS[7] << param.randomseed << endl;
  }
//---------------------------------------------------------------------------
Analysis::~Analysis() {
  delete parts;
  delete allelegroup;
  delete[] phenotype;
  delete[] marker;
  delete[] gendata;
  delete[] mdrsumres[CONTROL];
  delete[] mdrsumres[CASE];
  for (int i1=0; i1<N_MDR_PARTS; i1++) {
    delete[] mdrpartres[i1][CONTROL];
    delete[] mdrpartres[i1][CASE];
    }
  }
//---------------------------------------------------------------------------
