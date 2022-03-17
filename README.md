# pop-TAP

This  repository contains custom scripts used to asses transposon absence polymorphisms (TAPs) using pair-end short read illumina data.


### Example Usage:

**First Step:**

Run get_cov_sample.sh for each one of the samples (including controls) as folllows:

> bash  ./main/get_cov_sample.sh \
> -g ${genome_file} \
> -w 10 \
> -f ${first_pair.fq.gz} \
> -t 8

**Second Step:** 

Run the main script TE_CNV.sh for each sample to test ( you can use the same control  everytime):

> bash ./main/TE_CNV.sh \
> -i sample_depthGCcorrected.per10_bp.bed.gz \
> -r control_depthGCcorrected.per10_bp.bed.gz \
> -a TE_annotation_gff \
> -t 300  \
> -S ./main/Wilcoxon_test.R

Note: Because  TE_CNV.sh  only uses one thread, one can use GNU parallel to run all your samples faster:

> parallel --jobs 8 "bash ./main/TE_CNV.sh -i {} -r control_depthGCcorrected.per10_bp.bed.gz -a TE_annotation_gff \
> -t 300  -S ./main/Wilcoxon_test.R" ::: *.depthGCcorrected.per10_bp.bed.gz

In this example we run 8 samples at the same time. 
