#ifndef LOADER_H
#define LOADER_H
//------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <ctype.h>
#include <vector>
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
    int nselmarkers;
    char **selmarker;

    Loader();
    bool loadSelectedMarkers(string filename);
    virtual bool loadFile(string filename, string phenoname, MDR::Analysis *analysis)=0;
    ~Loader();
  };
//------------------------------------------------------------------------------
class ExampleLoader : public Loader {
  public:
    bool loadFile(string filename, string phenoname, MDR::Analysis *analysis);
  };
//------------------------------------------------------------------------------
#endif // LOADER_H
