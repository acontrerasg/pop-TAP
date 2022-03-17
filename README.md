# popTAP


Rationale: <br>
Most of available tools developed to detect genome-wide copy number variations (CNV) make use complex algorithms because they are challenged by two issues:
the length and location of the CNV is unknown and also the type of CNV to detect is unkwnown. A review of the topic can be found here: https://doi.org/10.3389/fgene.2015.00138. Transposon absence polymorhpisms, however, are a type of CNV in which location, length and type of CNV to detect are known before hand and, therefore, more simple algorithms are sufficient to reliable detect this type of CNV. <br>

This pipeline asses transposon absence polymorphisms (TAPs) using pair-end short read illumina data. It first calculates GC corrected read depth coverage (RD) per bin  across whole genome and then compares a test sample against a user defined control sample to detected significant changes of RD for every TE of at least 300 bp. A  non-parametric approach comparing a sample versus control TE related bin populations (with a minimun of 30 bins per population) is used to test significance. If the tested TE has a significant sample vs control differnce in RD and also its median RD shows at least a 10 fold reduction than median sample genome-wide RD, it is considered a TE absence polymorphism. 

To run this pipeline several packages and tools must be installed and be accessible through  the $PATH : 
 
 - bwa
 - mosdepth
 - samtools
 - bedtools
 - awk
 - bioawk
 - samblaster
 - vroom (R library)
 - dplyr (R library)

It also requeires a GFF file with the TE annotation of the reference genome. 
The reference genome must have been indexed with bwa index. 

<br>

### Example Usage:

**First Step:**

Run get_cov_sample.sh for each one of the samples (including controls) as folllows:

> bash ./main/get_cov_sample.sh \
> -g genome_file \
> -w 10 \
> -f first_pair.fq.gz \
> -t 8

This will create a file with the extension "depthGCcorrected.per10_bp.bed.gz" to the base name of the read pair.
Important: This script expects the pair end raw reads to be adaptor free and to have the following extensions:
 - First mate of the pair: "_1.fastq.gz"
 - Second mate of the pair: "_2.fastq.gz"

**Second Step:** 

Run the main script TE_CNV.sh for each sample to test (you can use the same control  everytime):

> bash ./main/TE_CNV.sh \
> -i sample_depthGCcorrected.per10_bp.bed.gz \
> -r control_depthGCcorrected.per10_bp.bed.gz \
> -a TE_annotation_gff \
> -t 300  \
> -S ./main/Wilcoxon_test.R

Note: Because TE_CNV.sh only uses one thread, one can use GNU parallel to run all your samples faster:

> parallel --jobs 8 "bash ./main/TE_CNV.sh -i {} -r control_depthGCcorrected.per10_bp.bed.gz -a TE_annotation.gff \
> -t 300  -S ./main/Wilcoxon_test.R" ::: *.depthGCcorrected.per10_bp.bed.gz

In this example we run 8 samples at the same time. 

**Final output**

For each sample file run with TE_CNV.sh. Two final tables must have beein created:
 -  summarized_output.tsv file  with  the  postion of all tested reference TE plus Pvalues and median calculated coverages. 
 -  TAPs.tsv  same format as previous file  but only the tested reference TEs considered absent in the  tested sample are present. 
