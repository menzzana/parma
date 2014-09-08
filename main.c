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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>     
#include "mpi.h"
#include "version_config.h"
#include "global.h"
#include "mdr.h"
//------------------------------------------------------------------------------
int main(int argc, char **argv) {
  int optionvalue;

  printVersion();
  while ((optionvalue=getopt(argc,argv,"f:p:m:s:"))!=END_OF_OPTIONS)
    switch (optionvalue) {
      case 'f':
        // filename=optarg;
        break;
      case 'p':
        // permutations=optarg
        break;
      case 'm':
        // number of markers=optarg
        break;
      case 's':
        sran1(atol(optarg));
        break;
      default:
        exit(EXIT_FAILURE);
      }
  exit(EXIT_SUCCESS);
  }
//------------------------------------------------------------------------------
