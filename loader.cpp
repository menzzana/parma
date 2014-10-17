#include "loader.h"
//---------------------------------------------------------------------------
bool ExampleLoader::loadFile(string filename, string phenoname, MDR::Analysis *analysis) {
  ifstream fpr;
  string fstr;
  int idxind,idxmark,i1;

  try {
    fpr.open(filename.c_str());
    getline(fpr,fstr);
    for (idxmark=0; idxmark<fstr.length(); idxmark++)
      if (fstr[idxmark]==delimiter)
        analysis->nmarkers++;
    for (analysis->nindividuals=0; getline(fpr,fstr); analysis->nindividuals++);
    analysis->createDataBuffers(true);
    fpr.clear();
    fpr.seekg(0,ios_base::beg);
    getline(fpr,fstr);
    for (idxmark=i1=0; idxmark<analysis->nmarkers; idxmark++) {
      strcpy(analysis->marker[idxmark],fstr.substr(i1,fstr.find(delimiter,i1)-i1).c_str());
      i1=fstr.find(delimiter,i1)+1;
      }
    for (idxind=0; idxind<analysis->nindividuals; idxind++) {
      if (!getline(fpr,fstr))
        throw runtime_error("Unexpected End of file");
      for (idxmark=i1=0; idxmark<analysis->nmarkers; idxmark++) {
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
