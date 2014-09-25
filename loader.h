#ifndef LOADER_H
#define LOADER_H
//------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <ctype.h>
#include "global.h"
#include "mdr.h"
//------------------------------------------------------------------------------
const char delimiter='\t';
//------------------------------------------------------------------------------
using namespace std;
//------------------------------------------------------------------------------
class Loader {
  public:
    int nindividuals,nmarkers;
    string *marker;
    unsigned char **gendata,*phenotype;

    Loader();
    virtual bool loadFile(string filename, string phenoname)=0;
    void setData(MDR::Analysis data);
    ~Loader();
  };
//------------------------------------------------------------------------------
class ExampleLoader : public Loader {
  public:
    bool loadFile(string filename, string phenoname);
  };
//------------------------------------------------------------------------------
#endif // LOADER_H
