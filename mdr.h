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
typedef enum { false, true } bool;

static const int N_MDR_PARTS=10;
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
void initialMask(unsigned char *mask, int nanalysedmarkers, int nmarkers);
bool incMask(unsigned char *mask, int nmarkers);
void andMask(char *oldmask,char *mask, int nmarkers);
bool TestMask(char *oldmask,char *mask, int nmarkers);
void randomshuffle(unsigned char *data, int length);
void mdr(unsigned char **data, unsigned char *phenotype, int nindividuals, int nphenotype);
//------------------------------------------------------------------------------
#endif // MDR_H
