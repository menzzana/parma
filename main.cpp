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
  int optionvalue;

  printVersion();
  try {
    while ((optionvalue=getopt(argc,argv,"f:p:m:s:t:d:"))!=global::END_OF_OPTIONS)
      switch (optionvalue) {
        case 'f':
          filename=optarg;
          break;
        case 'p':
          myanalysis.permutations=atoi(optarg);
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
            case 1: // Data format Example Loader class
              mydata=new ExampleLoader();
              break;
            case 2: // Data format DB schizophrenia class
              break;
            }
          break;
        default:
          THROW_ERROR("Unknown option: -"+optopt);
        }
    if (mydata==NULL)
      THROW_ERROR("Genotype data format not set");
    if (!mydata->loadFile(filename, phenoname))
      THROW_ERROR("Cannot load data file: "+filename);
    mydata->setData(myanalysis);
    myanalysis.frommarker=0;
    myanalysis.tomarker=mydata->nmarkers;
    myanalysis.setInitialArrays();
    myanalysis.Run();
    cleanUp(EXIT_SUCCESS);
    }
  catch(exception &e) {
    cerr << e.what() << endl;
    cleanUp(EXIT_FAILURE);
    }
  }
//------------------------------------------------------------------------------
