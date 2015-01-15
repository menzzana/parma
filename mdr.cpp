#include "mdr.h"
//---------------------------------------------------------------------------
using namespace MDR;
//---------------------------------------------------------------------------
SummedData::SummedData() {
  calc.error=1;
  calc.nnegpermutations=0;
  tp=fp=tn=fn=0;
  memset(parterror,0,sizeof(double)*N_MDR_PARTS);
  }
//---------------------------------------------------------------------------
void SummedData::clearPartData() {
  tp=fp=tn=fn=0;
  }
//---------------------------------------------------------------------------
void SummedData::addError(int idxpart) {
  float caseerror,controlerror;

  caseerror=(tp+fp)==0?0:fp/(tp+fp);
  controlerror=(fn+tn)==0?0:fn/(fn+tn);
  parterror[idxpart]=(caseerror+controlerror)/2;
  calc.error+=parterror[idxpart];
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
  return calc1.error>calc2.error;
  }
//---------------------------------------------------------------------------
#ifndef SERIAL
  void SummedData::procTestBestCombination(Calculated *in, Calculated *inout, int *len, MPI_Datatype *type) {
    if (testBestCombination(*inout,*in)) {
      inout->nnegpermutations=in->nnegpermutations;
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
  cutpvalue=NO_CUTOFF;
  gendata=NULL;
  phenotype=NULL;
  marker=NULL;
  permpheno=NULL;
  parts=NULL;
  allelegroup=NULL;
  param.npermutations=0;
  param.maxcombinations=MAX_MARKER_COMBINATIONS;
  param.mincombinations=MIN_MARKER_COMBINATIONS;
  groupn=0;
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
  parts=new unsigned char[param.nindividuals];
  allelegroup=new short[param.nindividuals];
  mdrsumres[CONTROL]=new short[param.nindividuals];
  mdrsumres[CASE]=new short[param.nindividuals];
  for (int i1=0; i1<N_MDR_PARTS; i1++) {
    mdrpartres[CONTROL][i1]=new short[param.nindividuals];
    mdrpartres[CASE][i1]=new short[param.nindividuals];
    }
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
void Analysis::initializePartPermutationArrays() {
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
void Analysis::checkMaxCombination() {
  param.maxcombinations=min(param.maxcombinations,param.nmarkers);
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
void Analysis::clearMDRResults() {
  memset(mdrsumres[CONTROL],0,sizeof(short)*groupn);
  memset(mdrsumres[CASE],0,sizeof(short)*groupn);
  for (int i1=0; i1<N_MDR_PARTS; i1++) {
    memset(mdrpartres[CONTROL][i1],0,sizeof(short)*groupn);
    memset(mdrpartres[CASE][i1],0,sizeof(short)*groupn);
    }
  }
//---------------------------------------------------------------------------
void Analysis::setAlleleGroups(int combinations) {
  int idxind1,idxind2,idxmark;

  memset(allelegroup,0,sizeof(short)*param.nindividuals);
  groupn=1;
  for (idxind1=0; idxind1<param.nindividuals; idxind1++) {
    if (allelegroup[idxind1]>0)
      continue;
    allelegroup[idxind1]=groupn;
    groupn++;
    for (idxind2=idxind1+1; idxind2<param.nindividuals; idxind2++) {
      for (idxmark=0; idxmark<combinations; idxmark++)
        if (gendata[idxind1][markercombo[idxmark]]!=gendata[idxind2][markercombo[idxmark]])
          break;
      if (idxmark==combinations)
        allelegroup[idxind2]=allelegroup[idxind1];
      }
    }
  }
//---------------------------------------------------------------------------
Result Analysis::analyseAlleles(unsigned char *vpheno, int combinations, Result *original) {
  Result accres;
   int idxind,idxgroup,idxpart;
   int traincase,traincontrol;

   clearMDRResults();
   for (idxind=0; idxind<param.nindividuals; idxind++) {
     mdrpartres[(int)vpheno[idxind]][(int)parts[idxind]][(int)allelegroup[idxind]]++;
     mdrsumres[(int)vpheno[idxind]][(int)allelegroup[idxind]]++;
     }
   accres=Result();
   for (idxpart=0; idxpart<N_MDR_PARTS; idxpart++) {
     accres.train.clearPartData();
     accres.test.clearPartData();
     for(idxgroup=0; idxgroup<groupn; idxgroup++) {
       traincase=mdrsumres[CASE][idxgroup]-mdrpartres[CASE][idxpart][idxgroup];
       traincontrol=mdrsumres[CONTROL][idxgroup]-mdrpartres[CONTROL][idxpart][idxgroup];
       if (traincase<traincontrol) {
         accres.train.fp+=traincase;
         accres.train.tn+=traincontrol;
         accres.test.fp+=mdrpartres[CASE][idxpart][idxgroup];
         accres.test.tn+=mdrpartres[CONTROL][idxpart][idxgroup];
         }
       else {
         accres.train.tp+=traincase;
         accres.train.fn+=traincontrol;
         accres.test.tp+=mdrpartres[CASE][idxpart][idxgroup];
         accres.test.fn+=mdrpartres[CONTROL][idxpart][idxgroup];
         }
       }
     accres.train.addError(idxpart);
     accres.test.addError(idxpart);
     if (original!=NULL) {
       if (accres.train.parterror[idxpart]>original->train.parterror[idxpart])
         original->train.calc.nnegpermutations++;
       if (accres.test.parterror[idxpart]>original->test.parterror[idxpart])
         original->test.calc.nnegpermutations++;
       }
     }
  if (original==NULL) {
    memcpy(accres.markercombo,markercombo,sizeof(int)*combinations);
    accres.combinations=combinations;
    accres.train.calc.error/=(double)N_MDR_PARTS;
    accres.test.calc.error/=(double)N_MDR_PARTS;
    }
  return accres;
  }
//---------------------------------------------------------------------------
bool Analysis::Run(int rank, int mpisize, int combination) {
  unsigned long long cidx,maxmarkercombos;
  Result origerror;

  try {
    minerror=Result();
    maxmarkercombos=CALC::C(param.nmarkers,combination);
    for (cidx=rank; cidx<maxmarkercombos; cidx+=mpisize) {
      setMarkerCombination(cidx,combination);
      setAlleleGroups(combination);
      origerror=analyseAlleles(phenotype,combination,NULL);
      for (int i1=0; i1<param.npermutations; i1++)
        analyseAlleles(permpheno[i1],combination,&origerror);
      if (origerror.test.getPvaluePermutations(param.npermutations)<=param.cutpvalue)
        origerror.print(marker,param.npermutations,false);
      minerror.setBestCombination(origerror);
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
  clog << PARAMETERS[2] << param.npermutations << endl;
  clog << PARAMETERS[3] << param.mincombinations << endl;
  clog << PARAMETERS[4] << param.maxcombinations << endl;
  if (param.cutpvalue==NO_CUTOFF)
    clog << PARAMETERS[6] << endl;
  else
    clog << PARAMETERS[5] << endl;
  if (param.randomseed==0)
    clog << PARAMETERS[8] << endl;
  else
    clog << PARAMETERS[7] << param.randomseed << endl;
  }
//---------------------------------------------------------------------------
Analysis::~Analysis() {
  for (int i1=0; i1<param.npermutations; i1++)
    delete permpheno[i1];
  delete permpheno;
  delete parts;
  delete allelegroup;
  delete phenotype;
  delete[] marker;
  delete[] gendata;
  delete mdrsumres[CONTROL];
  delete mdrsumres[CASE];
  for (int i1=0; i1<N_MDR_PARTS; i1++) {
    delete mdrpartres[CONTROL][i1];
    delete mdrpartres[CASE][i1];
    }
  }
//---------------------------------------------------------------------------
