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
using namespace std;
//------------------------------------------------------------------------------
class Loader {
  public:
    static const char delimiter='\t';
    static const unsigned char ASCII0=48;
    static const int STD=1;
    static const int DB=2;
    int nindividuals,nmarkers;
    string *marker;
    unsigned char **gendata,*phenotype;
    int *selectedmarkers;

    Loader();
    void setSelectedMarkers();
    virtual bool loadFile(string filename, string phenoname)=0;
    ~Loader();
  };
//------------------------------------------------------------------------------
class ExampleLoader : public Loader {
  public:
    bool loadFile(string filename, string phenoname);
  };
//------------------------------------------------------------------------------
#endif // LOADER_H
