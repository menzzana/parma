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
using namespace std;
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

  MPI::Init(argc,argv);
  MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
  try {
    myanalysis=new MDR::Analysis();
    markerfilename="";
    mpisize=MPI::COMM_WORLD.Get_size();
    mpirank=MPI::COMM_WORLD.Get_rank();
    if (mpirank==global::MPIROOT) {
      printVersion();
      maxcombinations=MDR::MAX_MARKER_COMBINATIONS;
      while ((optionvalue=getopt_long_only(argc,argv,"f:p:m:s:t:d:c:",long_options,&optionindex))!=global::END_OF_OPTIONS)
        switch (optionvalue) {
          case 'f':
            filename=optarg;
            break;
          case 'p':
            myanalysis->npermutations=atoi(optarg);
            break;
          case 'c':
            myanalysis->maxcombinations=atoi(optarg);
            break;
          case 'm':
            markerfilename=optarg;
            break;
          case 's':
            RND::sran1(-atol(optarg));
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
      }
    MPI::COMM_WORLD.Bcast(&myanalysis->nmarkers,1,MPI_INT,global::MPIROOT);
    MPI::COMM_WORLD.Bcast(&myanalysis->nindividuals,1,MPI_INT,global::MPIROOT);
    myanalysis->createDataBuffers(mpirank!=global::MPIROOT);
    MPI::COMM_WORLD.Bcast(&myanalysis->phenotype[0],myanalysis->nindividuals,MPI_UNSIGNED_CHAR,global::MPIROOT);
    MPI::COMM_WORLD.Bcast(&myanalysis->gendata[0][0],myanalysis->nmarkers*myanalysis->nindividuals,
                          MPI_UNSIGNED_CHAR,global::MPIROOT);
    MPI::COMM_WORLD.Bcast(&myanalysis->marker[0][0],myanalysis->nmarkers*global::MAX_LENGTH_MARKER_NAME,
                          MPI_CHAR,global::MPIROOT);
    if (!myanalysis->Run(mpirank,mpisize))
      throw runtime_error("Cannot analyse data");
    exitvalue=EXIT_SUCCESS;
    }
  catch (MPI::Exception e) {
    cerr << "MPI ERROR: " << e.Get_error_code() << " - " << e.Get_error_string() << endl;
    exitvalue=EXIT_FAILURE;
    }
  catch(exception &e) {
    cerr << e.what() << endl;
    exitvalue=EXIT_FAILURE;
    }
  cleanUp();
  MPI::Finalize();
  exit(exitvalue);
  }
//------------------------------------------------------------------------------
