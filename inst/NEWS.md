Changes in version 0.99.15 (2021-10-20)
+ Submitted to Bioconductor

Changes in version 0.99.16 (2021-10-20)
+ Fixed a bug that probabilities can be 0 for all genotypes when mismap rate is 1.

Changes in version 0.99.26 (2021-11-24)
+ Reformatted the scripts and the vignette to meet Bioconductor's instruction.

Changes in version 0.99.27 (2021-12-21)
+ Added detailed instruction for building a scheme object in the vignette.

Changes in version 0.99.28 (2021-12-21)
+ Added a code to close the file connection to GDS after opening file via 
+ openGDS() called in the exmple section of the openGDS function.

Changes in version 0.99.29 (2021-12-23)
+ Fixed a bug in getRead() to get read counts of specified chromosome only.
+ Removed duplicated entries in Methods-GbsrGenotypeData.R.

Changes in version 0.99.30 (2022-1-15)
+ Fixed a bug in gbsrGDS2VCF().
+ Reorganized cpp script.

Changes in version 0.99.31 (2022-2-14)
+ Reorganize scripts.

Changes in version 0.99.32 (2022-2-14)
+ Fix typo in gbsrGDS2VCF().
+ Removed the test script for unused function.
