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
/*
  Namespace MDR for specific MDR related functionality
  Prior of running the variables...
  permutations nindividuals nmarkers
  .. and the arrays...
  gendata phenotype marker
  ...must be set
  Do use function createDataBuffers() and setInitialArrays() on each process
  Then call Run(rank,blocksize) for calculations
*/
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
      double tp,fp,tn,fn;
      struct Calculated {
        double nnegpermutations;
        double accuracy;
        int rank;
        } calc;

      SummedData();
      void copy(SummedData summeddata);
      void setAccuracy();
      double getPvaluePermutations(int npermutations);
      static bool testBestCombination(Calculated calc1, Calculated calc2);
      static void procTestBestCombination(Calculated *in, Calculated *inout, int *len, MPI_Datatype *type);
    };
//------------------------------------------------------------------------------
  class Result {
    public:
      static const char delimiter='\t';
      int markercombo[MAX_MARKER_COMBINATIONS],combinations;
      SummedData train,test;

      Result();
      void copy(Result result);
      void testBestCombination(Result result);
      static void printHeader(bool ispermutation);
      void print(char **marker, int npermutations, bool highest);
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
      int getAlleleCombinations(int combinations);
      bool setInitialCombination(int idxmark, int combinations);
      bool increaseCombination(int idxmarkcombo, int combinations);
      bool increaseMarker(int idxmarkcombo);
      void clearMDRResults(int combinations);
      Result analyseAlleles(unsigned char *vpheno, int combinations);

    public:
      struct Param {
        int npermutations,maxcombinations,nmarkers,nindividuals;
        long randomseed;
        double cutpvalue;
        } param;
      unsigned char **gendata,*phenotype;
      double cutpvalue;
      char **marker;
      Result maxaccuracy;

      Analysis();
      void setInitialArrays();
      void createDataBuffers(bool initthisrank);
      bool Run(int rank, int blocksize, int combination);
      void printBestResult();
      ~Analysis();
    };
  }
//------------------------------------------------------------------------------
#endif // MDR_H
