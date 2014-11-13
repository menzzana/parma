#include "loader.h"
//---------------------------------------------------------------------------
Loader::Loader() {
  nselmarkers=0;
  selmarker=NULL;
  }
//---------------------------------------------------------------------------
bool Loader::loadSelectedMarkers(string filename) {
  ifstream fpr;
  string fstr;
  vector<string> tmpmarker;

  try {
    fpr.open(filename.c_str());
    while (!getline(fpr,fstr))
      tmpmarker.push_back(fstr);
    fpr.close();
    nselmarkers=tmpmarker.size();
    selmarker=new char*[nselmarkers];
    selmarker[0]=new char[nselmarkers*global::MAX_LENGTH_MARKER_NAME];
    strncpy(selmarker[0],tmpmarker[0].c_str(),global::MAX_LENGTH_MARKER_NAME-1);
    for (int i1=1; i1<nselmarkers; i1++) {
      selmarker[i1]=&selmarker[0][i1*global::MAX_LENGTH_MARKER_NAME];
      strncpy(selmarker[i1],tmpmarker[0].c_str(),global::MAX_LENGTH_MARKER_NAME-1);
      }
    return true;
    }
  catch(exception &e) {
    cerr << e.what() << endl;
    return false;
    }
  }
//---------------------------------------------------------------------------
Loader::~Loader() {
  delete[] selmarker;
  }
//---------------------------------------------------------------------------
bool ExampleLoader::loadFile(string filename, string phenoname, MDR::Analysis *analysis) {
  ifstream fpr;
  string fstr;
  int idxind,idxmark,i1;

  try {
    fpr.open(filename.c_str());
    if (!fpr.good())
      return false;
    getline(fpr,fstr);
    for (idxmark=0; idxmark<fstr.length(); idxmark++)
      if (fstr[idxmark]==delimiter)
        analysis->param.nmarkers++;
    for (analysis->param.nindividuals=0; getline(fpr,fstr); analysis->param.nindividuals++);
    analysis->createDataBuffers(true);
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
      analysis->phenotype[idxind]=(phenoname[0]==fstr[i1]?1:0);
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
