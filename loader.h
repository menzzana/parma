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
class Loader {
  public:
    static const char delimiter='\t';
    static const unsigned char ASCII0=48;
    static const int STD=1;               // Data format Example Loader class
    static const int SPDB=2;              // Data format DB schizophrenia class
    static const int SPMARKER=0;
    static const int SPINDIVIDUAL=1;
    static const int SPPHENOTYPE=2;
    static const int SPALLELE1=3;
    static const int SPALLELE2=4;
    static const int MAX_DATA_COLUMNS=5;
    static const int BED=3;               // Data format for BED plink files
    vector <string> selmarker;

    Loader();
    bool loadSelectedMarkers(string filename);
    int getIndex(vector <string> myvec, string searchstring);
    int getIndex(vector <char[30]> myvec, string searchstring);
    void splitDataString(string fstr,string *data);
    virtual bool loadFile(string filename, string phenoname, MDR::Analysis *analysis)=0;
    ~Loader();
  };
//------------------------------------------------------------------------------
class ExampleLoader : public Loader {
  public:
    bool loadFile(string filename, string phenoname, MDR::Analysis *analysis);
  };
//------------------------------------------------------------------------------
class SPLoader : public Loader {
  public:
    bool loadFile(string filename, string phenoname, MDR::Analysis *analysis);
  };
//------------------------------------------------------------------------------
class MarkerList {
  public:
    string name;
    int size,index,count;
    unsigned char alleletype;
    bool complete_analysis;
    MarkerList *Next;

    MarkerList();
    MarkerList *add(string name, unsigned char alleletype);
    ~MarkerList();
  };
//------------------------------------------------------------------------------
class IndividualList {
  public:
    string name;
    int size,index;
    IndividualList *Next;

    IndividualList();
    IndividualList *add(string name);
    ~IndividualList();
  };
//------------------------------------------------------------------------------
#endif // LOADER_H
