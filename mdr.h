#ifndef MDR_H
#define MDR_H
//------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
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

static const int N_MDR_PARTS=10;
static const int MAX_MARKER_COMBINATIONS=100;

typedef enum {
  false, true
  } bool;
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
void setInitialCombination(int *markercombo, int markerno, int nmarkers, int ncombo);
bool increaseCombination(int *markercombo, int markerno, int nmarkers, int ncombo);
void populateMDRParts(unsigned char *data, int length);
void randomShuffle(unsigned char *data, int length);
void mdr(unsigned char **gendata, unsigned char *phenotype, int *marker,
         int nmarkers, int nindividuals, int nselectmarkers, int permutations);
//------------------------------------------------------------------------------
#endif // MDR_H
