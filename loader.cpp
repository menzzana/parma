#include "loader.h"
//---------------------------------------------------------------------------
Loader::Loader() {
  selmarker.clear();
  nmarkers=nindividuals=0;
  phenotype=NULL;
  gendata=NULL;
  marker=NULL;
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
  delete phenotype;
  delete[] gendata;
  delete[] marker;
  }
//---------------------------------------------------------------------------
void Loader::createDataBuffers() {
  phenotype=new unsigned char[nindividuals];
  gendata=new unsigned char*[nindividuals];
  gendata[0]=new unsigned char[nindividuals*nmarkers];
  for (int i1=1; i1<nindividuals; i1++)
    gendata[i1]=&gendata[0][i1*nmarkers];
  memset(&gendata[0][0],255,nindividuals*nmarkers);
  marker=new char*[nmarkers];
  marker[0]=new char[nmarkers*global::MAX_LENGTH_MARKER_NAME];
  for (int i1=1; i1<nmarkers; i1++)
    marker[i1]=&marker[0][i1*global::MAX_LENGTH_MARKER_NAME];
  }
//---------------------------------------------------------------------------
int Loader::removeNonGenotypeIndividuals() {
  int idxindividual,idxnewindividual,idxmarker;

  for (idxindividual=idxnewindividual=0; idxindividual<nindividuals; idxindividual++) {
    phenotype[idxnewindividual]=phenotype[idxindividual];
    for (idxmarker=0; idxmarker<nmarkers; idxmarker++) {
      if (gendata[idxindividual][idxmarker]==255)
        break;
      gendata[idxnewindividual][idxmarker]=gendata[idxindividual][idxmarker];
      }
    if (idxmarker==nmarkers)
      idxnewindividual++;
    }
  cerr << "Removed " << nindividuals-idxnewindividual << " individuals" << endl;
  return idxnewindividual;
  }
//---------------------------------------------------------------------------
void Loader::copy(MDR::Analysis *analysis) {
  memcpy(analysis->phenotype[0],phenotype,nindividuals);
  memcpy(analysis->gendata[0],gendata[0],nindividuals*nmarkers);
  memcpy(analysis->marker[0],marker[0],nmarkers*global::MAX_LENGTH_MARKER_NAME);
  }
//---------------------------------------------------------------------------
bool ExampleLoader::loadFile(string filename) {
  ifstream fpr;
  string fstr,markerstr;
  int idxind,idxmark,i1;

  try {
    fpr.open(filename.c_str());
    if (!fpr.good())
      return false;
    getline(fpr,fstr);
    for (idxmark=nmarkers=0; idxmark<fstr.length(); idxmark++)
      if (fstr[idxmark]==delimiter)
        nmarkers++;
    for (nindividuals=0; getline(fpr,fstr); nindividuals++);
    createDataBuffers();
    fpr.clear();
    fpr.seekg(0,ios_base::beg);
    getline(fpr,markerstr);
    for (idxind=0; idxind<nindividuals; idxind++) {
      if (!getline(fpr,fstr))
        throw runtime_error("Unexpected End of file");
      for (idxmark=i1=0; idxmark<nmarkers; idxmark++) {
        gendata[idxind][idxmark]=fstr[i1]-ASCII0;
        i1=fstr.find(delimiter,i1)+1;
        }
      phenotype[idxind]=fstr[i1]-ASCII0;
      }
    fpr.close();
    for (idxmark=i1=0; idxmark<nmarkers; idxmark++) {
      strncpy(marker[idxmark],markerstr.substr(i1,markerstr.find(delimiter,i1)-i1).c_str(),
              global::MAX_LENGTH_MARKER_NAME-1);
      i1=markerstr.find(delimiter,i1)+1;
      }
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
bool SPLoader::loadFile(string filename) {
  const char BASES[]="ACGTacgt";
  ifstream fpr;
  string fstr,data[MAX_DATA_COLUMNS];
  MarkerList *marker1,*pmarker,*ml1;
  IndividualList *individual1,*il1,*pindividual;
  int i1;

  try {
    fpr.open(filename.c_str());
    if (!fpr.good())
      return false;
    getline(fpr,fstr);
    marker1=NULL;
    individual1=NULL;
    while (getline(fpr,fstr)) {
      splitDataString(fstr,data);
      if (strchr(BASES,data[SPALLELE1][0])==NULL && data[SPALLELE1].length()>1)
        continue;
      if (selmarker.size()!=0 && getIndex(selmarker,data[SPMARKER])<0)
        continue;
      ml1=marker1->add(data[SPMARKER],min(toupper(data[SPALLELE1][0]),toupper(data[SPALLELE2][0])));
      if (marker1==NULL)
        marker1=ml1;
      il1=individual1->add(data[SPINDIVIDUAL]);
      if (individual1==NULL)
        individual1=il1;
      }
    nmarkers=marker1->size;
    for (ml1=marker1,i1=0; ml1!=NULL; ml1=ml1->Next)
      if (ml1->count<(SPPROCCOMPLETEANALYSIS*(float)individual1->size)) {
        marker1->size--;
        ml1->complete_analysis=false;
        }
      else
        ml1->index=i1++;
    cerr << "Removed " << nmarkers-marker1->size << " markers" << endl;
    nmarkers=marker1->size;
    nindividuals=individual1->size;
    createDataBuffers();
    fpr.clear();
    fpr.seekg(0,ios_base::beg);
    getline(fpr,fstr);
    while (getline(fpr,fstr)) {
      splitDataString(fstr,data);
      pmarker=getEntry(marker1,data[SPMARKER]);
      pindividual=getEntry(individual1,data[SPINDIVIDUAL]);
      if (pindividual==NULL || pmarker==NULL)
        continue;
      if (!pmarker->complete_analysis)
        continue;
      gendata[pindividual->index][pmarker->index]=
          data[SPALLELE1][0]!=data[SPALLELE2][0]?MDR::HETEROZYGOTE:
          data[SPALLELE1][0]==pmarker->alleletype?MDR::HOMOZYGOTE1:MDR::HOMOZYGOTE2;
      phenotype[pindividual->index]=(data[SPPHENOTYPE][0]-ASCII0)>0?MDR::CASE:MDR::CONTROL;
      }
    fpr.close();
    nindividuals=removeNonGenotypeIndividuals();
    for (ml1=marker1; ml1!=NULL; ml1=ml1->Next)
      if (ml1->complete_analysis)
        strcpy(marker[ml1->index],ml1->name.c_str());
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
