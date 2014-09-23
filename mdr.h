#ifndef MDR_H
#define MDR_H
//------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <iostream>
#include <algorithm>
//------------------------------------------------------------------------------
using namespace std;
//------------------------------------------------------------------------------
//Namespace MDR for specific MDR related functionslity
// Prior of running the variables...
// permutations nindividuals nmarkers frommarker tomarker
// .. and the arrays...
// gendata phenotype marker
// ...must be set
// Then call setInitialArrays() & Run() for calculations
//------------------------------------------------------------------------------
namespace MDR {
  static const unsigned char HOMOZYGOTE1=0;
  static const unsigned char HETEROZYGOTE=1;
  static const unsigned char HOMOZYGOTE2=2;
  static const unsigned char CONTROL=0;
  static const unsigned char CASE=1;
  static const int MAX_MARKER_COMBINATIONS=100;
  static const int N_MDR_PARTS=10;
  static const int PHENOTYPE_COMBINATIONS=2;
  static const int ALLELE_COMBINATIONS=3;
  // LIST_ALLELE_MARKER_COMBINATIONS=MAX_MARKER_COMBINATIONS^ALLELE_COMBINATIONS
  static const int LIST_ALLELE_MARKER_COMBINATIONS=10E6;
//------------------------------------------------------------------------------
  class SummedData {
    public:
      float tp,fp,tn,fn;
      float accuracy,npospermutations;

      void clear();
      void setAccuracy();
    };
//------------------------------------------------------------------------------
  class Result {
    public:
      int markercombo[MAX_MARKER_COMBINATIONS];
      SummedData train,test;
    };
//------------------------------------------------------------------------------
  class Analysis {
    private:
      unsigned char **permpheno,*parts;
      int markercombo[MAX_MARKER_COMBINATIONS];
      int mdrpartres[N_MDR_PARTS][PHENOTYPE_COMBINATIONS][LIST_ALLELE_MARKER_COMBINATIONS];
      int mdrsumres[PHENOTYPE_COMBINATIONS][LIST_ALLELE_MARKER_COMBINATIONS];

      void populateMDRParts();
      void randomShuffle(unsigned char *data);
      void setInitialCombination(int *markercombo, int idxmark, int combinations);
      bool increaseCombination(int *markercombo, int idxmark, int combinations);
      void clearMDRResults(int combinations);
      Result analyseAlleles(int *markercombo, unsigned char *vpheno, int combinations);

    public:
      int permutations,nindividuals,nmarkers,frommarker,tomarker;
      int *marker;
      unsigned char **gendata,*phenotype;

      Analysis();
      void setInitialArrays();
      void Run();
      ~Analysis();
    };
  }
//------------------------------------------------------------------------------
#endif // MDR_H
