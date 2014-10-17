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
#include "global.h"
//------------------------------------------------------------------------------
using namespace std;
//------------------------------------------------------------------------------
//Namespace MDR for specific MDR related functionslity
// Prior of running the variables...
// permutations nindividuals nmarkers frommarker tomarker
// .. and the arrays...
// gendata phenotype marker
// ...must be set
// Then call Run(rank,blocksize) for calculations
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
  static const double NO_CUTOFF=-1;
//------------------------------------------------------------------------------
  class SummedData {
    public:
      float tp,fp,tn,fn;
      float accuracy,nnegpermutations;

      SummedData();
      void copy(SummedData summeddata);
      void setAccuracy();
      double getPvaluePermutations(int npermutations);
    };
//------------------------------------------------------------------------------
  class Result {
    public:
      Result();
      void copy(Result result);
      void testBestCombination(Result result, int npermutations);
      void print(char **marker, int npermutations);
      int markercombo[MAX_MARKER_COMBINATIONS],combinations;
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
      void setInitialArrays();
      void randomShuffle(unsigned char *data);
      int getAlleleCombinations(int combinations);
      bool setInitialCombination(int idxmark, int combinations);
      bool increaseCombination(int idxmarkcombo, int combinations);
      bool increaseMarker(int idxmarkcombo);
      void clearMDRResults(int combinations);
      Result analyseAlleles(unsigned char *vpheno, int combinations);

    public:
      int nindividuals,nmarkers,npermutations,maxcombinations;
      unsigned char **gendata,*phenotype;
      double cutpvalue;
      char **marker;

      Analysis();
      void createDataBuffers(bool initthisrank);
      bool Run(int rank, int blocksize);
      ~Analysis();
    };
  }
//------------------------------------------------------------------------------
#endif // MDR_H
