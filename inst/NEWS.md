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

Changes in version 0.99.32 (2022-2-16)
+ Fix typo in gbsrGDS2VCF().
+ Removed the test script for unused function.

Changes in version 0.99.33 (2022-2-16)
+ Reduced the example codes for the Scheme object.

Changes in version 0.99.34 (2022-2-24)
+ Fix a bug in getGenotype not to flip corrected genotypes.

Changes in version 0.99.35 (2022-3-16)
+ Fix a bug in the test script test_subset_n_gds2vcf.R.

Changes in version 0.99.36 (2022-3-18)
+ Add closeGDS(gds) at the end of sample script for gbsrGDS2CSV().

Changes in version 0.99.37 (2022-3-20)
+ No change, but just for version bump.

Changes in version 0.99.38 (2022-3-20)
+ Change the man page for gbsrGDS2CSV().

Changes in version 0.99.39 (2022-3-24)
+ Reorganized NAMESPACE.

Changes in version 0.99.40 (2022-4-25)
+ Fix a bug in plotReadRatio().
