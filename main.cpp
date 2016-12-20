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
#include "messages.h"
#include "boost/program_options.hpp"
//------------------------------------------------------------------------------
namespace prgm_opt=boost::program_options;
//------------------------------------------------------------------------------
Loader *mydata;
MDR::Analysis *myanalysis;
//------------------------------------------------------------------------------
void cleanUp(bool exitvalue) {
  delete mydata;
  delete myanalysis;
  exit(exitvalue);
  }
//------------------------------------------------------------------------------
int main(int argc, char **argv) {
  string filename,markerfilename;
  int mpirank,mpisize,exitvalue;
  prgm_opt::variables_map option_map;
  prgm_opt::options_description options("Options");
  MDR::SummedData::Calculated maxresult;
  #ifndef SERIAL
    MPI_Datatype MPI_2DOUBLE_INT,MPI_5INT_LONG_DOUBLE_BOOL;
    MPI_Op MPI_BESTCOMBINATION;
  #endif

  try {
    #ifndef SERIAL
      if (MPI_Init(&argc,&argv)!=MPI_SUCCESS)
        THROW_ERROR(ERRORTEXT::NO_MPI);
      MPI_Comm_rank(MPI_COMM_WORLD,&mpirank);
      MPI_Comm_size(MPI_COMM_WORLD,&mpisize);
      MPI_Type_create_struct(global::LENGTH_2DOUBLE_INT,global::BLOCK_2DOUBLE_INT,global::DISP_2DOUBLE_INT,
                             global::TYPE_2DOUBLE_INT,&MPI_2DOUBLE_INT);
      MPI_Type_commit(&MPI_2DOUBLE_INT);
      MPI_Type_create_struct(global::LENGTH_5INT_LONG_DOUBLE_BOOL,global::BLOCK_5INT_LONG_DOUBLE_BOOL,
                             global::DISP_5INT_LONG_DOUBLE_BOOL,global::TYPE_5INT_LONG_DOUBLE_BOOL,&MPI_5INT_LONG_DOUBLE_BOOL);
      MPI_Type_commit(&MPI_5INT_LONG_DOUBLE_BOOL);
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
      // Program options
      prgm_opt::arg="[Value]";
      options.add_options()
        (CMDOPTIONS::HELP_OPTION[0],CMDOPTIONS::HELP_OPTION[2])
        (CMDOPTIONS::FILE_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::FILE_OPTION[2])
        (CMDOPTIONS::MARKERFILE_OPTION[0],prgm_opt::value<string>()->required(),CMDOPTIONS::MARKERFILE_OPTION[2])
        (CMDOPTIONS::PERMUTATION_OPTION[0],prgm_opt::value<int>()->required(),CMDOPTIONS::PERMUTATION_OPTION[2])
        (CMDOPTIONS::MAXCOMBO_OPTION[0],prgm_opt::value<int>()->required(),CMDOPTIONS::MAXCOMBO_OPTION[2])
        (CMDOPTIONS::MINCOMBO_OPTION[0],prgm_opt::value<int>()->required(),CMDOPTIONS::MINCOMBO_OPTION[2])
        (CMDOPTIONS::DATATYPE_OPTION[0],prgm_opt::value<int>()->required(),CMDOPTIONS::DATATYPE_OPTION[2])
        (CMDOPTIONS::CUTOFF_OPTION[0],prgm_opt::value<double>()->required(),CMDOPTIONS::CUTOFF_OPTION[2])
        (CMDOPTIONS::PERMUTEONE_OPTION[0],CMDOPTIONS::PERMUTEONE_OPTION[2])
        (CMDOPTIONS::GLOBALCOMBO_OPTION[0],prgm_opt::value<int>()->required(),CMDOPTIONS::GLOBALCOMBO_OPTION[2]);
      if (option_map.count(CMDOPTIONS::HELP_OPTION[1])) {
        cout << options;
        cleanUp(EXIT_SUCCESS);
        }
      if (option_map.count(CMDOPTIONS::PERMUTATION_OPTION[1]))
        myanalysis->param.npermutations=option_map[CMDOPTIONS::PERMUTATION_OPTION[1]].as<int>();
      if (option_map.count(CMDOPTIONS::FILE_OPTION[1]))
        filename=option_map[CMDOPTIONS::FILE_OPTION[1]].as<string>();
      if (option_map.count(CMDOPTIONS::CUTOFF_OPTION[1]))
        myanalysis->param.cutpvalue=option_map[CMDOPTIONS::CUTOFF_OPTION[1]].as<double>();
      if (option_map.count(CMDOPTIONS::MAXCOMBO_OPTION[1]))
        myanalysis->param.maxcombinations=option_map[CMDOPTIONS::MAXCOMBO_OPTION[1]].as<int>();
      if (option_map.count(CMDOPTIONS::MINCOMBO_OPTION[1]))
        myanalysis->param.mincombinations=option_map[CMDOPTIONS::MINCOMBO_OPTION[1]].as<int>();
      if (option_map.count(CMDOPTIONS::MARKERFILE_OPTION[1]))
        markerfilename=option_map[CMDOPTIONS::MARKERFILE_OPTION[1]].as<string>();
      if (option_map.count(CMDOPTIONS::SEED_OPTION[1]))
        myanalysis->param.randomseed=-option_map[CMDOPTIONS::SEED_OPTION[1]].as<int>();
      if (option_map.count(CMDOPTIONS::DATATYPE_OPTION[1]))
        switch(option_map[CMDOPTIONS::DATATYPE_OPTION[1]].as<int>()) {
          case Loader::STD:
            mydata=new ExampleLoader();
            break;
          case Loader::SPDB:
            mydata=new SPLoader();
            break;
          case Loader::BED:
            break;
          }
      if (option_map.count(CMDOPTIONS::PERMUTEONE_OPTION[1]))
        myanalysis->param.onlypermuteone=true;
      if (mydata==NULL)
        THROW_ERROR(ERRORTEXT::NO_GENOTYPE_DATA_FORMAT);
      if (markerfilename.length()>0)
        if (!mydata->loadSelectedMarkers(markerfilename))
          THROW_ERROR_VALUE(ERRORTEXT::NO_MARKERS,markerfilename);
      if (!mydata->loadFile(filename))
        THROW_ERROR_VALUE(ERRORTEXT::NO_FILE_LOAD,filename);
      myanalysis->param.nindividuals=mydata->nindividuals;
      myanalysis->param.nmarkers=mydata->nmarkers;
      myanalysis->checkParameters();
      myanalysis->printParameters();
      MDR::Result::printHeader(myanalysis->param.npermutations>0);
      }
    #ifndef SERIAL
      if (MPI_Bcast(&myanalysis->param,1,MPI_5INT_LONG_DOUBLE_BOOL,global::MPIROOT,MPI_COMM_WORLD)!=MPI_SUCCESS)
        THROW_ERROR(ERRORTEXT::NO_PARAMETER_SEND);
    #endif
    CALC::sran1(myanalysis->param.randomseed);
    myanalysis->createDataBuffers();
    if (mpirank==global::MPIROOT)
      mydata->copy(myanalysis);
      #ifndef SERIAL
      if (MPI_Bcast(&myanalysis->phenotype[0][0],myanalysis->param.nindividuals,
                    MPI_UNSIGNED_CHAR,global::MPIROOT,MPI_COMM_WORLD)!=MPI_SUCCESS)
        THROW_ERROR(ERRORTEXT::NO_PHENOTYPE_SEND);
      if (MPI_Bcast(&myanalysis->gendata[0][0],myanalysis->param.nmarkers*myanalysis->param.nindividuals,
                    MPI_UNSIGNED_CHAR,global::MPIROOT,MPI_COMM_WORLD)!=MPI_SUCCESS)
        THROW_ERROR(ERRORTEXT::NO_GENETIC_DATA_SEND);
      if (MPI_Bcast(&myanalysis->marker[0][0],myanalysis->param.nmarkers*global::MAX_LENGTH_MARKER_NAME,
                    MPI_CHAR,global::MPIROOT,MPI_COMM_WORLD)!=MPI_SUCCESS)
        THROW_ERROR(ERRORTEXT::NO_MARKER_SEND);
      #endif
    myanalysis->initializePartPermutationArrays();
    for (int ncombo=myanalysis->param.mincombinations; ncombo<=myanalysis->param.maxcombinations; ncombo++) {
      if (!myanalysis->Run(mpirank,mpisize,ncombo))
        THROW_ERROR(ERRORTEXT::NO_DATA_ANALYSED);
      myanalysis->minerror.test.calc.rank=mpirank;
      #ifndef SERIAL
        if (MPI_Allreduce(&myanalysis->minerror.test.calc,&maxresult,1,
                          MPI_2DOUBLE_INT,MPI_BESTCOMBINATION,MPI_COMM_WORLD)!=MPI_SUCCESS)
          THROW_ERROR(ERRORTEXT::NO_REDUCTION);
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
  #ifndef SERIAL
    MPI_Finalize();
  #endif
    cleanUp(exitvalue);
  }
//------------------------------------------------------------------------------
