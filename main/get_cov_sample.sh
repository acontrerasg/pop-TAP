#!/usr/bin/env bash

#This script first maps samples to the reference genome and then uses mosdepth to calculate the coverage in windows.
#uses bioawk, bwa mem, mosdepth and bedtools.

usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }

[ $# -eq 0 ] && usage

while getopts ":h:w:f:g:b:t:" arg; do
  case $arg in
    g) # path to genome in fasta format
      genome=${OPTARG}
      ;;
    w) # Window size (in integers) to calculate coverage
         window_size=${OPTARG}
      ;;
    f) # path to the first read of the pair. Importat! Both reads shoudl be located in the same folder.
      fq=${OPTARG}
      ;;
    t) # threads (in integers)
      threads=${OPTARG}
      ;;
    h | *) # Display help.
      usage
      exit 0
      ;;
  esac
done



echo  "Start mapping with bwa mem..."
#Get read names and bam names:
read_name=`echo $fq | sed 's/_1.fastq.gz//'` ;
out_name=`basename $read_name` ;

#Map usig bwa men and 8 threads.

bwa mem -t $threads -c 1 ${genome} ${read_name}_1.fastq.gz ${read_name}_2.fastq.gz  |
samblaster |
samtools sort -@ ${threads} - |
samtools view -Sb -@ $threads > ${out_name}.bam
samtools index ${out_name}.bam


echo "Mapping done, starting the read depth coverage GC correction."
#Get names:
bam_name=${out_name}
win_size_str=$(echo $window_size"_bp")

#Temp file names:
temp_chrom_sizes=$(mktemp --suffix=.txt temp_chromsizes.XXXXXXXX)
temp_win100_gen=$(mktemp --suffix=.bed temp_100win_genome.XXXXXXXX)
temp_customwin_gen=$(mktemp --suffix=.bed temp_${win_size_str}win_genome.XXXXXXXX)
temp_win_gen_GC=$(mktemp --suffix=.bed temp_win_GC_genome.XXXXXXX)
temp_GC_calc=$(mktemp --suffix=.bed temp_windowed_genome.XXXXXXXX)
temp_GC_perwin=$(mktemp --suffix=.bed temp_windowed_GCarray.XXXXXXXX)
temp_median=$(mktemp --suffix=.txt temp_median.XXXXXXXX)
temp_mCG=$(mktemp --suffix=.txt temp_mCG.XXXXXXXX)
temp_median_per_bp_GCcont=$(mktemp --suffix=.bed temp_median_per_bp_GCcont.XXXXXXXX)
temp_cov_GC_100bp=$(mktemp --suffix=.bed temp_cov_GC_100bp.XXXXX)

echo Uses bioawk, bedtools and mosdepth

#####Main    ######

#calcualte length of each chromosome
echo calculating the length of each chromosome for $genome
bioawk -c fastx '{ print $name, length($seq) }' ${genome} >  ${temp_chrom_sizes}
#Slice the chromosome in non overalpping windows of 100 bp to calculate the GC content per 100bp windows.
bedtools makewindows -g ${temp_chrom_sizes} -w 100 > ${temp_win100_gen}
bedtools  nuc -fi ${genome} -bed ${temp_win100_gen} | tail -n +2 |
awk '{printf "%.2f\n", $5}' >  ${temp_GC_calc}
paste ${temp_win100_gen} ${temp_GC_calc} > ${temp_win_gen_GC}

#Calculate median read coverage per window:
echo Using mosdepth to calculate median read coverage per 100bp win in ${bam_name}.bam
echo Mosdepth skips as default reads marked as 1796 UNMAP,SECONDARY,QCFAIL,DUP.
mosdepth -t ${threads} -m -b  ${temp_win100_gen} ${bam_name}.100GC ${bam_name}.bam
#calculate overall median from the  same 100 bp windows used for GC correction.
echo calculating overall median
zcat ${bam_name}.100GC.per-base.bed.gz   | cut -f4 |
awk ' { a[i++]=$1; } END { print a[int(i/2)]; }' > ${temp_median}
echo overall median for ${bam_name} is  ${temp_median}
#GetGCperCentwinds
echo retrieving GC percentage per windowns of 100bp size....
zcat ${bam_name}.100GC.regions.bed.gz |
paste - ${temp_win_gen_GC} | cut -f1-4,8 > ${temp_cov_GC_100bp}


#calculate median reads per GC fraction value:
echo calculate median reads per GC fraction value... slow step....
cat ${temp_win_gen_GC} | cut -f4 | sort | uniq > ${temp_GC_perwin}
for value in $(cat ${temp_GC_perwin}) ;
do
cat ${temp_cov_GC_100bp} |
awk -v value=$value  '{if($5 == value) print ($4+1) }' |
sort -n |
awk  -v GCvalue=$value '{ a[i++]=$1; }END { x=int((i+1)/2);
if (x < (i+1)/2) print GCvalue"\t"(a[x-1]+a[x])/2;
else print GCvalue"\t"a[x-1]; }' ;
done > ${temp_mCG}

#Get custom window size coverage
echo  Get custom window size coverage
echo Size selected  ${window_size} bp
echo The smaller the size, the slower the process and the bigger the file.
bedtools makewindows -g ${temp_chrom_sizes} -w ${window_size} > ${temp_customwin_gen}
echo calculating the coverage of the sampel in ${win_size_str} bins. Ignoring the 1796 flag.
mosdepth -t ${threads} -m -n -b ${temp_customwin_gen} ${bam_name}.custom_coverage ${bam_name}.bam

####### CALCULATIONS BACKGROUND   ########
##########################################
#We sought to adjust the 100-bp window read counts based on the observed deviation incoverage for a given G+C percentage.
#For G+C percentages of 0, 1,2, 3,..., 100%, we determined the deviation of coverage from thegenome average.
#Then a simple adjustment was made accordingto the equation
#      ~rcorr =ri* m/mGC
#where ri are read counts of the ith window, mGC is the median read counts of all windows that have the same G+C percentage as theith window,
#and m is theoverall median of all the windows.


#From: Sensitive and accurate detection of copy numbervariants using read depth of coverage. Genome res. Yoon et al 2009.
##########################################


#Get GC corrected  coverage:
echo Get GC corrected  coverage for ${bam_name}... slow step...
#First step, assign a value of GC content ot each median read value and assign  its overall mCG median
zcat ${bam_name}.custom_coverage.regions.bed.gz |
bedtools  intersect -wb -a - -b  ${temp_win_gen_GC} |
cut -f1-4,8 |
awk 'NR==FNR {h[$1] = $2; next} {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" h[$5]}' ${temp_mCG} -  > ${temp_median_per_bp_GCcont}
#  The clever awk from above is from:
#  https://unix.stackexchange.com/questions/113898/how-to-merge-two-files-based-on-the-matching-of-two-columns
#It maps the given mGC of a corresponding RD median to each  file entry.

#Second step, perform the GC correction calculation.
echo perform the GC correction calculation
cat  ${temp_median_per_bp_GCcont} |
awk -v total_median=$(cat ${temp_median}) 'BEGIN{OFS="\t" ; FS=" "}
{RCGCcorr=(total_median/$6) ;  print $1,$2,$3,$4,($4*RCGCcorr)}' |
sed '1 i #Chrom\tStart\tEnd\tRC\tRCcorr' |
gzip >  ${bam_name}.depthGCcorrected.per${win_size_str}.bed.gz

##remove used temp files:
rm ${temp_chrom_sizes}
rm ${temp_win_gen}
rm ${temp_win_gen_GC}
rm ${temp_GC_calc}
rm ${temp_GC_perwin}
rm ${temp_median}
rm ${temp_mCG}
rm ${temp_median_per_bp_GCcont}
rm ${temp_custom_cov}
rm ${temp_customwin_gen}
rm ${temp_win100_gen}
rm ${temp_cov_GC_100bp}

rm ${bam_name}.100GC.mosdepth.summary.txt
rm ${bam_name}.100GC.mosdepth.global.dist.txt
rm ${bam_name}.100GC.regions.bed.gz
rm ${bam_name}.100GC.regions.bed.gz.csi
rm ${bam_name}.100GC.per-base.bed.gz
rm ${bam_name}.100GC.per-base.bed.gz.csi
rm ${bam_name}.100GC.mosdepth.region.dist.txt


rm ${bam_name}.custom_coverage.mosdepth.summary.txt
rm ${bam_name}.custom_coverage.mosdepth.global.dist.txt
rm ${bam_name}.custom_coverage.regions.bed.gz
rm ${bam_name}.custom_coverage.regions.bed.gz.csi
rm ${bam_name}.custom_coverage.mosdepth.region.dist.txt
