#ifndef MDR_H
#define MDR_H
//------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
//------------------------------------------------------------------------------
// Constants
//------------------------------------------------------------------------------
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define MAX(x,y) (x>y?x:y)
#define MIN(x,y) (x<y?x:y)

#define MAX_MARKER_COMBINATIONS 100
static const int N_MDR_PARTS=10;
static const int PHENOTYPE_COMBINATIONS=2;
static const int ALLELE_COMBINATIONS=3;
// LIST_ALLELE_MARKER_COMBINATIONS=MAX_MARKER_COMBINATIONS^ALLELE_COMBINATIONS
static const int LIST_ALLELE_MARKER_COMBINATIONS=10E6;
//------------------------------------------------------------------------------
// Definitions
//------------------------------------------------------------------------------
typedef enum {
  false, true
  } bool;
struct MDRData {
  float tp,fp,tn,fn;
  };
struct MDRAccuracy {
  int markercombo[MAX_MARKER_COMBINATIONS];
  float train,test,pvaluetrain,pvaluetest;
  };
//------------------------------------------------------------------------------
// Variables
//------------------------------------------------------------------------------
static long rseed;
//------------------------------------------------------------------------------
// Functions
//------------------------------------------------------------------------------
void init();
void sran1(long value);
double ran1();
void setInitialCombination(int markercombo[MAX_MARKER_COMBINATIONS], int markerno,
                           int nmarkers, int combinations);
bool increaseCombination(int markercombo[MAX_MARKER_COMBINATIONS], int markerno,
                         int nmarkers, int ncombo);
void populateMDRParts(unsigned char *data, int length);
void swap(unsigned char a, unsigned char b);
void randomShuffle(unsigned char *data, int length);
void clearMDRResults(int mdrpartres[N_MDR_PARTS][PHENOTYPE_COMBINATIONS][LIST_ALLELE_MARKER_COMBINATIONS],
                     int combinations);
float balancedAccuracy(struct MDRData res);
struct MDRAccuracy analyseAlleles(unsigned char **gendata, int markercombo[MAX_MARKER_COMBINATIONS],
                                 unsigned char *phenotype,unsigned char *mdrparts, int ncombo, int nindividuals);
void mdr(unsigned char **gendata, unsigned char *phenotype, int *marker,
         int nmarkers, int nindividuals, int nselectmarkers, int permutations);
//------------------------------------------------------------------------------
#endif // MDR_H
