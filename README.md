# pop-TAP

This repository contains custom scripts used to asses transposon absence polymorphisms (TAPs) using pair-end short read illumina data.
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

> bash  ./main/get_cov_sample.sh \
> -g ${genome_file} \
> -w 10 \
> -f ${first_pair.fq.gz} \
> -t 8

This will create a file with the  extension  "depthGCcorrected.per10_bp.bed.gz" to the base name of the  read pair.
Important: This script expects the pair end raw reads to be adaptor free and to have the following extensions:
 - First mate of the pair: "_1.fastq.gz"
 - Second mate of the pair: "_2.fastq.gz"

**Second Step:** 

Run the main script TE_CNV.sh for each sample to test ( you can use the same control  everytime):

> bash ./main/TE_CNV.sh \
> -i sample_depthGCcorrected.per10_bp.bed.gz \
> -r control_depthGCcorrected.per10_bp.bed.gz \
> -a TE_annotation_gff \
> -t 300  \
> -S ./main/Wilcoxon_test.R

Note: Because  TE_CNV.sh  only uses one thread, one can use GNU parallel to run all your samples faster:

> parallel --jobs 8 "bash ./main/TE_CNV.sh -i {} -r control_depthGCcorrected.per10_bp.bed.gz -a TE_annotation.gff \
> -t 300  -S ./main/Wilcoxon_test.R" ::: *.depthGCcorrected.per10_bp.bed.gz

In this example we run 8 samples at the same time. 

**Final output**

For each sample file run with TE_CNV.sh. Two final tables must have beein created:
 -  summarized_output.tsv file  with  the  postion of all tested reference TE plus Pvalues and median calculated coverages. 
 -  TAPs.tsv  same format as previous file  but only the tested reference TEs considered absent in the  tested sample are present. 
