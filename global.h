#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdexcept>
#include <iostream>
#include <sstream>

#ifndef SERIAL
  #include "mpi.h"
#endif

#define THROW_ERROR(text) throw runtime_error(text)
#define THROW_ERROR_VALUE(text,value) throw runtime_error(string(text)+value)
//------------------------------------------------------------------------------
// global constants
//------------------------------------------------------------------------------
namespace global {
  static const int END_OF_OPTIONS=-1;
  static const int MAX_LENGTH_MARKER_NAME=30;
  static const int MPIROOT=0;
  #ifndef SERIAL
    static const int LENGTH_2DOUBLE_INT=3;
    static MPI_Datatype TYPE_2DOUBLE_INT[LENGTH_2DOUBLE_INT]={
      MPI_DOUBLE,MPI_DOUBLE,MPI_INT
      };
    static int BLOCK_2DOUBLE_INT[LENGTH_2DOUBLE_INT]={
      1,1,1
      };
    static MPI_Aint DISP_2DOUBLE_INT[LENGTH_2DOUBLE_INT]={
      0,sizeof(double),2*sizeof(double)
      };
    static const int LENGTH_5INT_LONG_DOUBLE=7;
    static MPI_Datatype TYPE_5INT_LONG_DOUBLE[LENGTH_5INT_LONG_DOUBLE]={
      MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_LONG,MPI_DOUBLE
      };
    static int BLOCK_5INT_LONG_DOUBLE[LENGTH_5INT_LONG_DOUBLE]={
      1,1,1,1,1,1,1
      };
    static MPI_Aint DISP_5INT_LONG_DOUBLE[LENGTH_5INT_LONG_DOUBLE]={
      0,sizeof(int),2*sizeof(int),3*sizeof(int),4*sizeof(int),
      5*sizeof(int),5*sizeof(int)+sizeof(long)
      };
  #endif
  }
//------------------------------------------------------------------------------
namespace CALC {
  #define IA 16807
  #define IM 2147483647
  #define AM (1.0/IM)
  #define IQ 127773
  #define IR 2836
  #define NTAB 32
  #define NDIV (1+(IM-1)/NTAB)
  #define EPS 1.2e-7
  #define RNMX (1.0-EPS)
  static long rseed=-123456789;

  void sran1(long rseed);
  double ran1();
  unsigned long long C(unsigned long long n, unsigned long long k);
  }
//------------------------------------------------------------------------------
namespace ERRORTEXT {
  const char NO_MPI[]="Cannot init MPI";
  const char TYPE_README[]="See README for a list of options.";
  const char NO_GENOTYPE_DATA_FORMAT[]="Genotype data format not set";
  const char NO_MARKERS[]="Cannot load markers from file: ";
  const char NO_FILE_LOAD[]="Cannot load data file: ";
  const char NO_PARAMETER_SEND[]="Cannot broadcast parameters";
  const char NO_PHENOTYPE_SEND[]="Cannot broadcast phenotype data";
  const char NO_GENETIC_DATA_SEND[]="Cannot broadcast genetic data";
  const char NO_MARKER_SEND[]="Cannot broadcast marker names";
  const char NO_DATA_ANALYSED[]="Cannot analyse data";
  const char NO_REDUCTION[]="Cannot reduce max results from all processes";
  }
//------------------------------------------------------------------------------
#endif // GLOBAL_H
