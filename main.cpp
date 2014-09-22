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
//------------------------------------------------------------------------------
int main(int argc, char **argv) {
  int optionvalue, nindividuals, nmarkers, permutations, combinations;
  unsigned char **gendata, **markername, *phenotype;
  char filename[global::MAX_LENGTH_STRING],phenoname[global::MAX_LENGTH_STRING];
  char filenamemarkers[global::MAX_LENGTH_STRING];

  printVersion();
  while ((optionvalue=getopt(argc,argv,"f:p:m:s:t:c:"))!=global::END_OF_OPTIONS)
    switch (optionvalue) {
      case 'f':
        strncpy(filename,optarg,global::MAX_LENGTH_STRING-1);
        break;
      case 'p':
        permutations=atoi(optarg);
        break;
      case 'm':
        strncpy(filenamemarkers,optarg,global::MAX_LENGTH_STRING-1);
        break;
      case 's':
        ran1(atol(optarg));
        break;
      case 't':
        strncpy(phenoname,optarg,global::MAX_LENGTH_STRING-1);
        break;
      case 'c':
        combinations=atoi(optarg);
        break;
      default:

        fprintf(stderr,"Unknown option `-%c'.\n", optopt);
        exit(EXIT_FAILURE);
      }

  exit(EXIT_SUCCESS);
  }
//------------------------------------------------------------------------------
