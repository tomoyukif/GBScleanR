Changes in version 1.9.12 (2024-5-22)
+ Fix a bug in setParents() with setting `bi = TRUE`.

Changes in version 1.9.10 (2024-5-21)
+ Fix minor bugs in closeGDS(), estGeno(), and the documents.

Changes in version 1.9.3 (2024-5-16)
+ Modify the estGeno() function to make it fix the parental genotypes at the 
+ first half of markers, for the reverse-direction genotype estimation step, as 
+ those were estimated in the forward-direction genotype estimation step. This 
+ modification make the function to avoid estimating the different haplotype 
+ combinations as the border of the first half and last half of markers after 
+ the concatenation of the first half and last half of estimated haplotypes in 
+ the .halfJoint() function.

Changes in version 1.9.2 (2024-5-16)
+ Update the vignette.

Changes in version 1.9.1 (2024-5-15)
+ Change .minimizeRecombination() to make it ignore the cases that are not listed
+ on the possible haplotype combinations.
+ Change the coding of .solveBorderConflict() to make it faster.
+ Add gbsrGDS2CSV() to output a CSV file.
+ Add setReplicates() to allow replicates in the data and modified estGeno() to
+ sum up read counts of replicates for genotype estimation.
+ Add makeScheme() to allow a simpler assignment of the scheme information.

Changes in version 1.5.13 (2023-10-17)
+ Update manual pages for plotDosage() and plotReadRatio().

Changes in version 1.5.11 (2023-10-2)
+ Fix a bug in setParents() that produces NA if a user set mono = TRUE and/or bi = TRUE.

Changes in version 1.5.8 (2023-8-14)
+ A user of GBScleanR reported the case in which there are NA values in read
+ count data in the AD node after conversion from a VCF file and results in a 
+ fatal error during the genotype estimation process.
+ Since the VCF file that causes NA values in the AD node has not provided and 
+ we have no clue to fix the bug in the data conversion, we tentatively added 
+ the codes to set 0 to NA values in the read count data.

Changes in version 1.5.7 (2023-8-2)
+ Fix a bug in plotDosage() and plotReadRatio() where the selection of sample
+ ID had been wrongly conducted when sample filter have been set.

Changes in version 1.5.5 (2023-8-1)
+ Add a code to make EDS node in the given GDS after esetGeno() to store dosage information
+ when the given population is inbred biparental. 

Changes in version 1.5.4 (2023-6-24)
+ Add an argument `node` in `gbsrGDS2VCF()` to switch genotype data to be output
+ into a VCF file. The setting `node == cor` replace data in the `genotype` node
+ of the given GDS file with the corrected genotype data stored in the
+ `annotation/format/CGT` node.
+ Add arguments `info.import` and `fmt.import` in `gbsrVCF2GDS`, and 
+ `info.export`, `fmt.export` in `gbsrGDS2VCF` to specify the data to be 
+ imported from a given VCF file or exported to a VCF file.

Changes in version 1.5.3 (2023-5-11)
+ Fix a bug in setInfoFilter() in which markers having no value for specified INFO
+ had been coded to be filtered out but now those can be retained.

Changes in version 1.4.2 (2023-4-26)
+ Fix a typo in a man page.

Changes in version 1.4.1 (2023-4-26)
+ Fix a bug in the vignette.

Changes in version 1.3.18 (2023-4-20)
+ Fix a bug in an example code in a manual page.

Changes in version 1.3.17 (2023-4-20)
+ Edit some codes to fit Bioconductor's requirements.

Changes in version 1.3.16 (2023-4-18)
+ Add the parentless mode in which estGeno() assign the number of dummy reads 
`dummy_reads` to dummy parents and estimate genotypes.

Changes in version 1.2.15 (2023-4-18)
+ Fix error found in build reports of Bioconductor

Changes in version 1.2.14 (2023-4-6)
+ Upon the request for multiple pedigree in the scheme, estGeno() supports a 
population consisting of samples that belongs to different pedigrees. 
+ Deprecated the crosstype argument in initScheme().
+ Add assignScheme() to let user specify which sample belongs to which pedigree.

Changes in version 1.2.9 (2023-2-14)
+ Minor change in genotype probability calculation algorithm which offset 0.5% 
when the minimum genotype probability is zero.

Changes in version 1.2.7 (2023-2-10)
+ Minor change in genotype probability calculation algorithm which set 
heterozygous calls impossible when the founders are inbred lines.

Changes in version 1.2.6 (2023-2-1)
+ Bug fix in HMM when a genotype is completely impossible to occur.
+ Bug fix in genotype pattern list for HMM when the pop is F1.

Changes in version 1.2.4 (2022-11-14)
+ Minor update in gbsrGDS2CSV().

Changes in version 1.2.3 (2022-11-14)
+ Bug fix in GBSR_HMM.cpp that caused a compile error in installation via bioconda.

Changes in version 1.2.1 (2022-11-11)
+ Minor modification in reference manual.

Changes in version 1.1.8 (2022-11-05)
+ Minor modification in gbsrGDS2CSV().

Changes in version 1.1.7 (2022-11-01)
+ Minor bug fix in gbsrGDS2CSV().

Changes in version 1.1.6 (2022-11-01)
+ Minor bug fix in gbsrGDS2CSV().

Changes in version 1.1.4 (2022-10-19)
+ Update in the refernece manual

Changes in version 1.1.4 (2022-10-19)
+ Minor bug fix in gbsrGDS2VCF().

Changes in version 1.1.3 (2022-10-19)
+ Minor bug fix in getRead().

Changes in version 1.1.2 (2022-10-19)
+ Fix a typo in the vignette

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

Changes in version 1.0.2 (2022-9-7)
+ Fix a bug in GBSR_HMM.cpp.

Changes in version 1.0.5 (2022-9-7)
+ Rewrote some scripts and replaced some function names for better user experience and code readability

Changes in version 1.0.6 (2022-10-17)
+ Rewrote vignette and reference manuals

