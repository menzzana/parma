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
  {"perm",required_argument,0,'p'},
  {"mfile",required_argument,0,'m'},
  {"seed",required_argument,0,'s'},
  {"pheno",required_argument,0,'t'},
  {"dtype",required_argument,0,'d'},
  {"comb",required_argument,0,'c'},
  {"cpvalue",required_argument,0,'u'},
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
  string filename,phenoname, markerfilename;
  int optionvalue,mpirank,mpisize,maxcombinations,optionindex,exitvalue;
  long randomseed;
  MPI_Datatype MPI_2DOUBLE_INT;
  MPI_Op MPI_BESTPVALUE;
  MDR::SummedData::Calculated procresult,maxresult;

  try {
    if (MPI_Init(&argc,&argv)!=MPI_SUCCESS)
      throw runtime_error("Cannot init MPI");
    MPI_Comm_rank(MPI_COMM_WORLD,&mpirank);
    MPI_Comm_size(MPI_COMM_WORLD,&mpisize);
    MPI_Type_create_struct(global::LENGTH_2DOUBLE_INT,global::BLOCK_2DOUBLE_INT,global::DISP_2DOUBLE_INT,
                           global::TYPE_2DOUBLE_INT,&MPI_2DOUBLE_INT);
    MPI_Type_commit(&MPI_2DOUBLE_INT);
    MPI_Op_create((MPI_User_function *)MDR::SummedData::procTestBestCombination, true, &MPI_BESTPVALUE);
    myanalysis=new MDR::Analysis();
    markerfilename="";
    maxcombinations=MDR::MAX_MARKER_COMBINATIONS;
    randomseed=0;
    if (mpirank==global::MPIROOT) {
      printVersion();
      while ((optionvalue=getopt_long_only(argc,argv,"f:p:m:s:t:d:c:",long_options,&optionindex))!=global::END_OF_OPTIONS)
        switch (optionvalue) {
          case 'f':
            filename=optarg;
            break;
          case 'p':
            myanalysis->npermutations=atoi(optarg);
            break;
          case 'c':
            maxcombinations=atoi(optarg);
            break;
          case 'm':
            markerfilename=optarg;
            break;
          case 's':
            randomseed=-atol(optarg);
            break;
          case 't':
            phenoname=optarg;
            break;
          case 'd':
            switch(atoi(optarg)) {
              case Loader::STD:
                mydata=new ExampleLoader();
                break;
              case Loader::HZDB:
                break;
              }
            break;
          case 'u':
            myanalysis->cutpvalue=atof(optarg);
            break;
          default:
            throw runtime_error("See README for a list of options.");
          }
      if (mydata==NULL)
        throw runtime_error("Genotype data format not set");
      if (markerfilename.length()>0)
        if (!mydata->loadSelectedMarkers(markerfilename))
          throw runtime_error("Cannot load markers from file: "+markerfilename);
      if (!mydata->loadFile(filename, phenoname, myanalysis))
        throw runtime_error("Cannot load data file: "+filename);
      MDR::Result::printHeader(myanalysis->npermutations>0);
      }
    if (MPI_Bcast(&myanalysis->npermutations,1,MPI_INT,global::MPIROOT,MPI_COMM_WORLD)!=MPI_SUCCESS)
      throw runtime_error("Cannot broadcast permutations");
    if (MPI_Bcast(&maxcombinations,1,MPI_INT,global::MPIROOT,MPI_COMM_WORLD)!=MPI_SUCCESS)
      throw runtime_error("Cannot broadcast combinations");
    if (MPI_Bcast(&myanalysis->nmarkers,1,MPI_INT,global::MPIROOT,MPI_COMM_WORLD)!=MPI_SUCCESS)
      throw runtime_error("Cannot broadcast number of markers");
    if (MPI_Bcast(&myanalysis->nindividuals,1,MPI_INT,global::MPIROOT,MPI_COMM_WORLD)!=MPI_SUCCESS)
      throw runtime_error("Cannot broadcast individuals");
    if (MPI_Bcast(&randomseed,1,MPI_LONG,global::MPIROOT,MPI_COMM_WORLD)!=MPI_SUCCESS)
      throw runtime_error("Cannot broadcast random seed");
    if (MPI_Bcast(&myanalysis->cutpvalue,1,MPI_DOUBLE,global::MPIROOT,MPI_COMM_WORLD)!=MPI_SUCCESS)
      throw runtime_error("Cannot broadcast cutoff p-value");
    myanalysis->createDataBuffers(mpirank!=global::MPIROOT);
    if (MPI_Bcast(&myanalysis->phenotype[0],myanalysis->nindividuals,
                  MPI_UNSIGNED_CHAR,global::MPIROOT,MPI_COMM_WORLD)!=MPI_SUCCESS)
      throw runtime_error("Cannot broadcast phenotype vector");
    if (MPI_Bcast(&myanalysis->gendata[0][0],myanalysis->nmarkers*myanalysis->nindividuals,
                  MPI_UNSIGNED_CHAR,global::MPIROOT,MPI_COMM_WORLD)!=MPI_SUCCESS)
      throw runtime_error("Cannot broadcast genetic data");
    if (MPI_Bcast(&myanalysis->marker[0][0],myanalysis->nmarkers*global::MAX_LENGTH_MARKER_NAME,
                  MPI_CHAR,global::MPIROOT,MPI_COMM_WORLD)!=MPI_SUCCESS)
      throw runtime_error("Cannot broadcast marker names");
    RND::sran1(randomseed);
    myanalysis->setInitialArrays();
    for (int ncombo=1; ncombo<=maxcombinations; ncombo++) {
      if (!myanalysis->Run(mpirank,mpisize,ncombo))
        throw runtime_error("Cannot analyse data");

      myanalysis->maxaccuracy.test.calc.rank=mpirank;

      procresult.nnegpermutations=myanalysis->maxaccuracy.test.nnegpermutations;
      procresult.accuracy=myanalysis->maxaccuracy.test.accuracy;
      procresult.rank=mpirank;


      if (MPI_Allreduce(&procresult,&maxresult,1,MPI_2DOUBLE_INT,MPI_BESTPVALUE,MPI_COMM_WORLD)!=MPI_SUCCESS)
        throw runtime_error("Cannot reduce max results from all processes");
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
  MPI_Finalize();
  exit(exitvalue);
  }
//------------------------------------------------------------------------------
