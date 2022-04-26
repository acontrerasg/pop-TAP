# popTAP


**Rationale** <br>
Most available tools developed to detect genome-wide copy number variations (CNV) make use of complex algorithms because they are challenged by two issues: the length and location of the CNV is unknown and also the type of CNV to detect is unknown. A review of the topic can be found here: https://doi.org/10.3389/fgene.2015.00138.  Transposon absence polymorphisms, however, are a type of CNV in which location, length and variant signature are known beforehand and, therefore, more straightforward algorithms are sufficient to reliably detect this type of CNV. <br>

**Pipeline description** <br>
This pipeline evaluates transposon absence polymorphisms (TAPs) using pair-end short read illumina data. It first calculates GC corrected read depth coverage (RD) per user-defined bin (recommended 10bp) across the whole genome of an array of samples and controls. Then, for every TE present in the provided gff annotation file (recommended to exclude TEs smaller than 300 bp), a non-parametric approach is used to do sample versus control comparisons of bin populations (with a minimum of 30 bins per population if TE < 300 bp were excluded). If the tested TE shows a significant difference and also its median RD shows at least a 10 fold reduction than the median sample genome-wide RD, then that sample TE is considered an absence polymorphism.

**Requeriments** <br>
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

It also requeires a GFF file with only the TE annotation of the reference genome. 
The reference genome must have been indexed with bwa index. 

<br>

### Example Usage

**First Step**

Run get_cov_sample.sh for each one of the samples (including controls) as folllows:

> bash ./main/get_cov_sample.sh \
> -g genome_file \
> -w 10 \
> -f first_pair.fq.gz \
> -t 8

This will create a file with the extension "depthGCcorrected.per10_bp.bed.gz"  with the base name of the provided read pair.
Important: This script expects the pair end raw reads to be adaptor free and to have the same name just differing at the following extensions:
 - First mate of the pair: "_1.fastq.gz"
 - Second mate of the pair: "_2.fastq.gz"

**Second Step** 

Run the main script TE_CNV.sh for each sample to test (you can use the same control  everytime):

> bash ./main/TE_CNV.sh \
> -i sample_depthGCcorrected.per10_bp.bed.gz \
> -r control_depthGCcorrected.per10_bp.bed.gz \
> -a TE_annotation_gff \
> -t 300  \
> -S ./main/Wilcoxon_test.R

Note: Because TE_CNV.sh only uses one thread, one can use GNU parallel to speed up this step:

> parallel --jobs 8 "bash ./main/TE_CNV.sh -i {} -r control_depthGCcorrected.per10_bp.bed.gz -a TE_annotation.gff \
> -t 300  -S ./main/Wilcoxon_test.R" ::: *.depthGCcorrected.per10_bp.bed.gz

In this example we run 8 samples at the same time. 

**Final output**

For each sample file run with TE_CNV.sh. Two final tables must have beein created:
 -  summarized_output.tsv file  with  the  postion of all tested reference TE plus Pvalues and median calculated coverages. 
 -  TAPs.tsv  same format as previous file  but only the reference TEs considered absent in the tested sample are present. 
