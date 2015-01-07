/********************************************************************************

Parma

Authors: Henric Zazzi, Hannes Leskelä
Copyright (C) PDC Center for High Performance Computing, KTH 2014  

Licensed under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Copy of the GNU General Public License can be onbtained from
see <http://www.gnu.org/licenses/>.
*******************************************************************************/
#include "version_config.h"
#include "global.h"
#include "mdr.h"
#include "loader.h"
#include <getopt.h>
//------------------------------------------------------------------------------
static struct option long_options[]={
  {"file",required_argument,0,'f'},
  {"permutations",required_argument,0,'p'},
  {"markerfile",required_argument,0,'m'},
  {"seed",required_argument,0,'s'},
  {"datatype",required_argument,0,'d'},
  {"maxcombinations",required_argument,0,'a'},
  {"mincombinations",required_argument,0,'i'},
  {"cutoffpvalue",required_argument,0,'u'},
  {0,0,0,0}
  };
Loader *mydata;
MDR::Analysis *myanalysis;
//------------------------------------------------------------------------------
void cleanUp() {
  delete mydata;
  delete myanalysis;
  }
//------------------------------------------------------------------------------
int main(int argc, char **argv) {
  string filename,markerfilename;
  int optionvalue,mpirank,mpisize,optionindex,exitvalue;
  MDR::SummedData::Calculated maxresult;
  #ifndef SERIAL
    MPI_Datatype MPI_2DOUBLE_INT,MPI_5INT_LONG_DOUBLE;
    MPI_Op MPI_BESTCOMBINATION;
  #endif

  try {
    #ifndef SERIAL
      if (MPI_Init(&argc,&argv)!=MPI_SUCCESS)
        throw runtime_error("Cannot init MPI");
      MPI_Comm_rank(MPI_COMM_WORLD,&mpirank);
      MPI_Comm_size(MPI_COMM_WORLD,&mpisize);
      MPI_Type_create_struct(global::LENGTH_2DOUBLE_INT,global::BLOCK_2DOUBLE_INT,global::DISP_2DOUBLE_INT,
                             global::TYPE_2DOUBLE_INT,&MPI_2DOUBLE_INT);
      MPI_Type_commit(&MPI_2DOUBLE_INT);
      MPI_Type_create_struct(global::LENGTH_5INT_LONG_DOUBLE,global::BLOCK_5INT_LONG_DOUBLE,
                             global::DISP_5INT_LONG_DOUBLE,global::TYPE_5INT_LONG_DOUBLE,&MPI_5INT_LONG_DOUBLE);
      MPI_Type_commit(&MPI_5INT_LONG_DOUBLE);
      MPI_Op_create((MPI_User_function *)MDR::SummedData::procTestBestCombination, true, &MPI_BESTCOMBINATION);
    #else
      mpirank=global::MPIROOT;
      mpisize=1;
      maxresult.rank=global::MPIROOT;
    #endif
    myanalysis=new MDR::Analysis();
    markerfilename="";
    if (mpirank==global::MPIROOT) {
      printVersion();
      while ((optionvalue=getopt_long_only(argc,argv,"f:p:m:s:t:d:c:",long_options,&optionindex))!=global::END_OF_OPTIONS)
        switch (optionvalue) {
          case 'f':
            filename=optarg;
            break;
          case 'p':
            myanalysis->param.npermutations=atoi(optarg);
            break;
          case 'u':
            myanalysis->param.cutpvalue=atof(optarg);
            break;
          case 'a':
            myanalysis->param.maxcombinations=atoi(optarg);
            break;
          case 'i':
            myanalysis->param.mincombinations=atoi(optarg);
            break;
          case 'm':
            markerfilename=optarg;
            break;
          case 's':
            myanalysis->param.randomseed=-atol(optarg);
            break;
          case 'd':
            switch(atoi(optarg)) {
              case Loader::STD:
                mydata=new ExampleLoader();
                break;
              case Loader::SPDB:
                mydata=new SPLoader();
                break;
              }
            break;
          default:
            throw runtime_error("See README for a list of options.");
          }
      if (mydata==NULL)
        throw runtime_error("Genotype data format not set");
      if (markerfilename.length()>0)
        if (!mydata->loadSelectedMarkers(markerfilename))
          throw runtime_error("Cannot load markers from file: "+markerfilename);
      if (!mydata->loadFile(filename, myanalysis))
        throw runtime_error("Cannot load data file: "+filename);
      myanalysis->checkMaxCombination();
      myanalysis->printParameters();
      MDR::Result::printHeader(myanalysis->param.npermutations>0);
      }
    #ifndef SERIAL
      if (MPI_Bcast(&myanalysis->param,1,MPI_5INT_LONG_DOUBLE,global::MPIROOT,MPI_COMM_WORLD)!=MPI_SUCCESS)
        throw runtime_error("Cannot broadcast parameters");
    #endif
    CALC::sran1(myanalysis->param.randomseed);
    myanalysis->createDataBuffers(mpirank!=global::MPIROOT);
    #ifndef SERIAL
      if (MPI_Bcast(&myanalysis->phenotype[0],myanalysis->param.nindividuals,
                    MPI_UNSIGNED_CHAR,global::MPIROOT,MPI_COMM_WORLD)!=MPI_SUCCESS)
        throw runtime_error("Cannot broadcast phenotype vector");
      if (MPI_Bcast(&myanalysis->gendata[0][0],myanalysis->param.nmarkers*myanalysis->param.nindividuals,
                    MPI_UNSIGNED_CHAR,global::MPIROOT,MPI_COMM_WORLD)!=MPI_SUCCESS)
        throw runtime_error("Cannot broadcast genetic data");
      if (MPI_Bcast(&myanalysis->marker[0][0],myanalysis->param.nmarkers*global::MAX_LENGTH_MARKER_NAME,
                    MPI_CHAR,global::MPIROOT,MPI_COMM_WORLD)!=MPI_SUCCESS)
        throw runtime_error("Cannot broadcast marker names");
    #endif
    myanalysis->setInitialArrays();
    for (int ncombo=myanalysis->param.mincombinations; ncombo<=myanalysis->param.maxcombinations; ncombo++) {
      if (!myanalysis->Run(mpirank,mpisize,ncombo))
        throw runtime_error("Cannot analyse data");
      myanalysis->maxaccuracy.test.calc.rank=mpirank;
      #ifndef SERIAL
        if (MPI_Allreduce(&myanalysis->maxaccuracy.test.calc,&maxresult,1,
                          MPI_2DOUBLE_INT,MPI_BESTCOMBINATION,MPI_COMM_WORLD)!=MPI_SUCCESS)
          throw runtime_error("Cannot reduce max results from all processes");
      #endif
      if (maxresult.rank==mpirank)
        myanalysis->printBestResult();
      }
    exitvalue=EXIT_SUCCESS;
    }
  catch(exception &e) {
    cerr << e.what() << endl;
    exitvalue=EXIT_FAILURE;
    }
  cleanUp();
  #ifndef SERIAL
    MPI_Finalize();
  #endif
  exit(exitvalue);
  }
//------------------------------------------------------------------------------
