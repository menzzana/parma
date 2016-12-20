#ifndef MESSAGES_H
#define MESSAGES_H
//------------------------------------------------------------------------------
// strings used for options
//------------------------------------------------------------------------------
namespace CMDOPTIONS {
  const char *const HELP_OPTION[]={"help,h","help","Displays available commands\n"};
  const char *const PERMUTATION_OPTION[]={"permutations,p","permutations","Sets the number of permutations. Default: 0\n"};
  const char *const FILE_OPTION[]={"file,f","file","The file(s) that should contain the genetic data\n"};
  const char *const MARKERFILE_OPTION[]={"markerfile,m","markerfile","The file containing the markers to be analyzed. Default: All\n"};
  const char *const MAXCOMBO_OPTION[]={"maxcombinations,a","maxcombinations","Max number of marker combinations test\n"};
  const char *const MINCOMBO_OPTION[]={"mincombinations,i","mincombinations","Min number of marker combinations tested\n"};
  const char *const DATATYPE_OPTION[]={"datatype,d","datatype","Input data format. 1) MDR format 2) Genetic DB format 3) plink format Default: 1\n"};
  const char *const CUTOFF_OPTION[]={"cutoffpvalue,u","cutoffpvalue","Show only test pvalues below this limit. Default: Best combination\n"};
  const char *const PERMUTEONE_OPTION[]={"onlypermuteone,o","onlypermuteone","Only performs permutation analysis with the loci combination with the lowest prediction error\n"};
  const char *const GLOBALCOMBO_OPTION[]={"globalcombinations,c","globalcombinations","Maximum number of analyzed combinations per marker combinations\n"};
  const char *const SEED_OPTION[]={"seed,s","seed","Specifies the random seed used by the analysis [Default: 123456789]\n"};
  }
//------------------------------------------------------------------------------
// Error Messages
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
#endif // MESSAGES_H
