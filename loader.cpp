#include "loader.h"
//---------------------------------------------------------------------------
Loader::Loader() {
  selmarker.clear();
  }
//---------------------------------------------------------------------------
bool Loader::loadSelectedMarkers(string filename) {
  ifstream fpr;
  string fstr;
  vector<string> selmarker;

  try {
    fpr.open(filename.c_str());
    while (!getline(fpr,fstr))
      selmarker.push_back(fstr);
    fpr.close();
    return true;
    }
  catch(exception &e) {
    cerr << e.what() << endl;
    return false;
    }
  }
//---------------------------------------------------------------------------
int Loader::getIndex(vector <string> myvec, string searchstring) {
  vector <string>::iterator i1;
  int i2;

  for (i1=myvec.begin(),i2=0; i1!=myvec.end(); i1++,i2++)
    if (searchstring.compare(*i1)==0)
      return i2;
  return -1;
  //i1=find(myvec.begin(),myvec.end(),searchstring);
  //return (i1==myvec.end()?-1:i1-myvec.begin());
  }
//---------------------------------------------------------------------------
void Loader::splitDataString(string fstr, string *data) {
  int i1,i2;

  for (i1=i2=0; i2<MAX_DATA_COLUMNS; i2++) {
    data[i2]=fstr.substr(i1,fstr.find(delimiter,i1)-i1);
    i1=fstr.find(delimiter,i1)+1;
    if (i1==0)
      break;
    }
  }
//---------------------------------------------------------------------------
Loader::~Loader() {
  selmarker.clear();
  }
//---------------------------------------------------------------------------
bool ExampleLoader::loadFile(string filename, MDR::Analysis *analysis) {
  ifstream fpr;
  string fstr;
  int idxind,idxmark,nmarker,nindividual,i1;

  try {
    fpr.open(filename.c_str());
    if (!fpr.good())
      return false;
    getline(fpr,fstr);
    for (idxmark=nmarker=0; idxmark<fstr.length(); idxmark++)
      if (fstr[idxmark]==delimiter)
        nmarker++;
    for (nindividual=0; getline(fpr,fstr); nindividual++);
    analysis->createMasterDataBuffers(nmarker,nindividual);
    fpr.clear();
    fpr.seekg(0,ios_base::beg);
    getline(fpr,fstr);
    for (idxmark=i1=0; idxmark<analysis->param.nmarkers; idxmark++) {
      strncpy(analysis->marker[idxmark],fstr.substr(i1,fstr.find(delimiter,i1)-i1).c_str(),
              global::MAX_LENGTH_MARKER_NAME-1);
      i1=fstr.find(delimiter,i1)+1;
      }
    for (idxind=0; idxind<analysis->param.nindividuals; idxind++) {
      if (!getline(fpr,fstr))
        throw runtime_error("Unexpected End of file");
      for (idxmark=i1=0; idxmark<analysis->param.nmarkers; idxmark++) {
        analysis->gendata[idxind][idxmark]=fstr[i1]-ASCII0;
        i1=fstr.find(delimiter,i1)+1;
        }
      analysis->phenotype[idxind]=fstr[i1]-ASCII0;
      }
    fpr.close();
    return true;
    }
  catch(exception &e) {
    cerr << e.what() << endl;
    return false;
    }
  }
//------------------------------------------------------------------------------
template<typename T> T *getEntry(T *first,string name) {
  T *tl1;

  for (tl1=first; tl1!=NULL; tl1=tl1->Next)
    if (name.compare(tl1->name)==0)
      return tl1;
  return NULL;
  }
//------------------------------------------------------------------------------
template<typename T> T *addEntry(T *first,string name) {
  T *tl1,*tl2;

  for (tl1=tl2=first; tl1!=NULL; tl2=tl1,tl1=tl1->Next)
    if (name.compare(tl1->name)==0)
      return tl1;
  tl1=new T();
  tl1->name=name;
  if (tl2!=NULL) {
    tl2->Next=tl1;
    tl1->index=first->size;
    first->size++;
    }
  return tl1;
  }
//------------------------------------------------------------------------------
bool SPLoader::loadFile(string filename, MDR::Analysis *analysis) {
  const char BASES[]="ACGTacgt";
  ifstream fpr;
  string fstr,data[MAX_DATA_COLUMNS];
  MarkerList *marker,*pmarker,*ml1;
  IndividualList *individual,*il1,*pindividual;
  int nmarkers,i1;

  try {
    fpr.open(filename.c_str());
    if (!fpr.good())
      return false;
    getline(fpr,fstr);
    marker=NULL;
    individual=NULL;
    while (getline(fpr,fstr)) {
      splitDataString(fstr,data);
      if (strchr(BASES,data[SPALLELE1][0])==NULL && data[SPALLELE1].length()>1)
        continue;
      if (selmarker.size()!=0 && getIndex(selmarker,data[SPMARKER])<0)
        continue;
      ml1=marker->add(data[SPMARKER],min(toupper(data[SPALLELE1][0]),toupper(data[SPALLELE2][0])));
      if (marker==NULL)
        marker=ml1;
      il1=individual->add(data[SPINDIVIDUAL]);
      if (individual==NULL)
        individual=il1;
      }
    nmarkers=marker->size;
    for (ml1=marker,i1=0; ml1!=NULL; ml1=ml1->Next)
      if (ml1->count<(SPPROCCOMPLETEANALYSIS*(float)individual->size)) {
        marker->size--;
        ml1->complete_analysis=false;
        }
      else
        ml1->index=i1++;
    cerr << "Removed " << nmarkers-marker->size << " markers" << endl;
    analysis->createMasterDataBuffers(marker->size,individual->size);
    for (ml1=marker; ml1!=NULL; ml1=ml1->Next)
      if (ml1->complete_analysis)
        strcpy(analysis->marker[ml1->index],ml1->name.c_str());
    fpr.clear();
    fpr.seekg(0,ios_base::beg);
    getline(fpr,fstr);
    while (getline(fpr,fstr)) {
      splitDataString(fstr,data);
      pmarker=getEntry(marker,data[SPMARKER]);
      pindividual=getEntry(individual,data[SPINDIVIDUAL]);
      if (pindividual==NULL || pmarker==NULL)
        continue;
      if (!pmarker->complete_analysis)
        continue;
      analysis->gendata[pindividual->index][pmarker->index]=
          data[SPALLELE1][0]!=data[SPALLELE2][0]?1:data[SPALLELE1][0]==pmarker->alleletype?0:2;
      analysis->phenotype[pindividual->index]=(data[SPPHENOTYPE][0]-ASCII0)>0?1:0;
      }
    fpr.close();
    analysis->removeNonGenotypeIndividuals();
    delete marker;
    delete individual;
    return true;
    }
  catch(exception &e) {
    cerr << e.what() << endl;
    return false;
    }
  }
//------------------------------------------------------------------------------
MarkerList::MarkerList() {
  name="";
  Next=NULL;
  size=1;
  index=0;
  alleletype=255;
  complete_analysis=true;
  count=0;
  }
//------------------------------------------------------------------------------
MarkerList::~MarkerList() {
  delete Next;
  }
//------------------------------------------------------------------------------
MarkerList *MarkerList::add(string name,unsigned char alleletype) {
  MarkerList *tl1;

  tl1=addEntry(this,name);
  tl1->count++;
  tl1->alleletype=min(tl1->alleletype,alleletype);
  return tl1;
  }
//------------------------------------------------------------------------------
IndividualList::IndividualList() {
  name="";
  Next=NULL;
  size=1;
  index=0;
  }
//------------------------------------------------------------------------------
IndividualList::~IndividualList() {
  delete Next;
  }
//------------------------------------------------------------------------------
IndividualList *IndividualList::add(string name) {
  return addEntry(this,name);
  }
//------------------------------------------------------------------------------
