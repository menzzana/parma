Parma
=====

Genes are the building blocks of heredity. They are passed from parent 
to child. They hold DNA, the instructions for making proteins, which in 
turn, produce energy, move molecules from one place to another, build 
structures, break down toxins, and do many other maintenance jobs.
Sometimes there is a mutation, a change in a gene or genes, which changes 
the gene's instructions for making a protein, so the protein does not work 
properly or is missing entirely. This can cause a medical condition called 
a genetic disorder. 
Identifying monogenic diseases, aka a genetic disorder where only one gene 
is involved, is straightforward since you only need to find one locus in 
close proximity to the gene in question where the distribution of the alleles 
differs between sick and healthy individuals.
In polygenic diseases the statistics can be very complex since there can be 
any number of genes involved and not all contributing genes have the same 
impact on the disorder.
In order to discover which combination of loci has the greatest impact on 
a disease, this software analyzes all loci combinations, and to further ascertain that 
the best combination of loci has been found, there is an option to permute the dataset 
a number of times, at the phenotype level, and directly compare the permuted 
results to the results obtained from the non-permuted dataset.
This software performs this analysis by using an algorithm called  
Multi‑Dimensionality Reduction (MDR) developed by Hahn et al.

L. W. Hahn, M. D. Ritchie, J. H. 2003. Multifactor dimensionality reduction software for
detecting gene–gene and gene–environment interactions. Bioinformatics 19:3. 376-382

BUILDING (from source bundle)
-----------------------------

Dependencies
^^^^^^^^^^^^

==================== ===============================================================
C++                  Works well with the GNU compiler
cmake                version 2.8+
Boost                Version 1.36+ http://www.boost.org/
MPI                  If not found a serial version will be compiled
==================== ===============================================================

How to build
^^^^^^^^^^^^

In order to build just run::

  mkdir [BUILD DIR]
  cd [BUILD DIR]
  cmake [SOURCE DIR]
  make

In case you want to run parma serially, uncomment find_package(MPI)
in *CMakeLists.txt*

Files
^^^^^

::

  CMakeLists.txt
  main.cpp
  mdr.cpp/.h
  global.cpp/.h
  loader.cpp/.h
  messages.h
  version_config.h.in

OPTIONS
-------

==============================  =================================================================================
-s,--seed [value]               Set [value] as new seed for random calculations
                                [value] should be a positive integer
-p,--permutations [value]       Set [value] as the number of permutations
-f,--file [file]                The [file] that should contain the genetic data
-m,--markerfile [file]          The [file] containing the markers to be analyzed
                                If not included all markers will be analyzed.
                                Markers in this file should be ordered in one marker per line with no header
-a,--maxcombinations [value]    Max number of marker combinations tested
-i,--mincombinations [value]    Min number of marker combinations tested
-d,--datatype [value]           Input data format

                                1) Example MDR data format from MDR
                                2) Genetic Database file format
                                3) plink binary file format
                                
-u,--cutoffpvalue [value]       Show only test pvalues below this limit.
                                Without this option only the best combinations will be shown
-o,--onlypermuteone             Only performs permutation analysis with the loci combination
                                with the lowest prediction error.
-c,--globalcombinations [value] Maximum number of analyzed combinations per marker combinations
==============================  =================================================================================

OUTPUT
------

The software outputs results in a TAB delimited format.

=========================== ==================================================
Markers                     The combination of markers analyzed
ClassificationError         The frequency of missclassified individuals in the
                            training dataset
PredictionError             The frequency of misspredicted individuals in the
                            testing dataset
ClassificationError p-value p-value of ClassificationError in the
                            permuted dataset
PredictionError p-value     p-value of PredictionError in the 
                            permuted dataset
MinError                    The best value in this specific marker combination
                            highlighted with a *
                            MinError is evaluated using PredictionError or
                            PredictionError p-value
=========================== ==================================================

USING
-----

Parma Can analyze a maxium of 3E18 combinations, but that will take a huge amount
of time, even if using a HPC cluster.
One process is capable of analyzing 3E6 combinations in one hour. Using parallel
execution you can multiple that number with the number of processes.
In order to restrict the number of combinations use the maxcombinations flag.

FILE FORMAT
-----------

PARMA is capable of reading different input files

Markerfile
^^^^^^^^^^

This should be a textbased file with no header and one marker per line

MDR Example file
^^^^^^^^^^^^^^^^

This is the original MDR file. More information about
the file can be found at http://www.multifactordimensionalityreduction.org/

Genetic Database file
^^^^^^^^^^^^^^^^^^^^^

This consist of a TAB delimited file sorted according to::

  [Marker]<TAB>[Individual Identity]<TAB>[Phenotype]<TAB>Allele 1<TAB>Allele 2
  
The first line should be the header of this file and phenotype is either 0 
for Healthy or 1 for Afflicted.

PLINK Binary file
^^^^^^^^^^^^^^^^^

Genetic data can also be in PLINK Binary file format (.bim/.fam/.bed)
