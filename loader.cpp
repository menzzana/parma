#include "loader.h"
//------------------------------------------------------------------------------
Loader::Loader() {
  nindividuals=nmarkers=0;
  marker=NULL;
  gendata=NULL;
  phenotype=NULL;
  }
//------------------------------------------------------------------------------
Loader::~Loader() {
  delete[] marker;
  delete phenotype;
  for (int i1=0; i1<nindividuals; i1++)
    delete gendata[i1];
  delete gendata;
  }
//------------------------------------------------------------------------------
void Loader::setData(MDR::Analysis data) {
  data.nmarkers=nmarkers;
  data.gendata=gendata;
  data.phenotype=phenotype;
  }
//---------------------------------------------------------------------------
bool ExampleLoader::loadFile(string filename, string phenoname) {
  ifstream fpr;
  streampos fpos;
  string fstr;
  int idxind,idxmark,i1;

  try {
    fpr.open(filename.c_str());
    getline(fpr,fstr);
    for (idxmark=0; idxmark<fstr.length(); idxmark++)
      if (fstr[idxmark]==delimiter)
        nmarkers++;
    marker=new string[nmarkers];
    for (idxmark=i1=0; idxmark<nmarkers; idxmark++) {
      marker[idxmark]=fstr.substr(i1,fstr.find(delimiter,i1)-i1+1);
      i1=fstr.find(delimiter,i1)+1;
      }
    fpos=fpr.tellg();
    for (nindividuals=0; getline(fpr,fstr); nindividuals++);
    fpr.clear();
    fpr.seekg(fpos,ios_base::beg);
    gendata=new unsigned char*[nindividuals];
    for (idxind=0; idxind<nindividuals; idxind++)
      gendata[idxind]=new unsigned char[nmarkers];
    phenotype=new unsigned char[nindividuals];
    for (idxind=0; idxind<nindividuals; idxind++) {
      if (!getline(fpr,fstr))
        THROW_ERROR("Unexpected End of file");
      for (idxmark=i1=0; idxmark<nmarkers; idxmark++) {
        gendata[idxind][idxmark]=fstr[i1];
        i1=fstr.find(delimiter,i1)+1;
        }
      phenotype[idxind]=fstr[i1];
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
