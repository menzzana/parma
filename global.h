#ifndef GLOBAL_H
#define GLOBAL_H

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdexcept>
#include "mpi.h"
//------------------------------------------------------------------------------
// global constants
//------------------------------------------------------------------------------
using namespace std;
//------------------------------------------------------------------------------
namespace global {
  static const int END_OF_OPTIONS=-1;
  static const int MAX_LENGTH_MARKER_NAME=30;
  static const int MAX_LENGTH_STRING=100;
  }
//------------------------------------------------------------------------------
// global functions
//------------------------------------------------------------------------------
#define THROW_ERROR(text) throw runtime_error(text)
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long option);
//------------------------------------------------------------------------------
#endif // GLOBAL_H
