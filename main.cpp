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
//------------------------------------------------------------------------------
using namespace std;
//------------------------------------------------------------------------------
Loader *mydata;
MDR::Analysis myanalysis;
//------------------------------------------------------------------------------
void cleanUp(int exitvalue) {
  delete mydata;
  exit(exitvalue);
  }
//------------------------------------------------------------------------------
int main(int argc, char **argv) {
  string filename,phenoname, filenamemarkers;
  int optionvalue,mpie,rank,numtasks,maxcombinations;

  printVersion();
  maxcombinations=MDR::MAX_MARKER_COMBINATIONS;
  try {
    while ((optionvalue=getopt(argc,argv,"f:p:m:s:t:d:c:"))!=global::END_OF_OPTIONS)
      switch (optionvalue) {
        case 'f':
          filename=optarg;
          break;
        case 'p':
          myanalysis.npermutations=atoi(optarg);
          break;
        case 'c':
          maxcombinations=atoi(optarg);
          break;
        case 'm':
          filenamemarkers=optarg;
          break;
        case 's':
          ran1(atol(optarg));
          break;
        case 't':
          phenoname=optarg;
          break;
        case 'd':
          switch(atoi(optarg)) {
            case Loader::STD: // Data format Example Loader class
              mydata=new ExampleLoader();
              break;
            case Loader::DB: // Data format DB schizophrenia class
              break;
            }
          break;
        default:
          throw runtime_error("Unknown option: -"+optopt);
        }
    if (mydata==NULL)
      throw runtime_error("Genotype data format not set");
    if (!mydata->loadFile(filename, phenoname))
      throw runtime_error("Cannot load data file: "+filename);
    mydata->setSelectedMarkers();
    myanalysis.setParameters(mydata->nmarkers,mydata->nindividuals,mydata->gendata,mydata->phenotype,mydata->selectedmarkers);
    myanalysis.setInitialArrays();
    /*
    mpie=MPI_Init(&argc,&argv);
    if (mpie!=MPI_SUCCESS)
      throw runtime_error("Error starting MPI program. Terminating.");

    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    printf ("Hello, World from rank %d out of %d\n", rank, numtasks);
    */
    for (int i1=1; i1<=maxcombinations; i1++)
      if (!myanalysis.Run(0,mydata->nmarkers,i1))
        throw runtime_error("Cannot analyse data");
    //MPI_Finalize();
    cleanUp(EXIT_SUCCESS);
    }
  catch(exception &e) {
    cerr << e.what() << endl;
    //MPI_Abort(MPI_COMM_WORLD, mpie);
    cleanUp(EXIT_FAILURE);
    }
  }
//------------------------------------------------------------------------------
