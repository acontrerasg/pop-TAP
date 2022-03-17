#!/usr/bin/env Rscript

suppressMessages(library(vroom))
suppressMessages(library(dplyr))

#Lazy way of parsing args.... maybe  use library("optparse")?
# for now use:
#${Rscritp}  ${temp_Rinput} ${temp_mean_input}  ${temp_mean_ref} ${TE_anno}

#load up files and variables:
##############
#Main file
args = commandArgs(trailingOnly=TRUE)
input_file <- vroom(args[1], col_names = TRUE, progress = FALSE, col_types = cols())
#means
input_sum=scan(args[2], quiet=TRUE)
ref_sum=scan(args[3], quiet=TRUE)
#TEanno
annotation <- vroom(args[4], col_names = FALSE, progress = FALSE)
colnames(annotation) <- c("#chrom", "source" , "feature", "start", "end", ".", "strand", "score", "TE_ID")

#Outrput names:
#Get path and basename of inpit file 
name <- as.character(args[1])
path <- dirname(name)
basename <- basename(name)
name <-  gsub(pattern = ".depthGCcorrected.per10_bp.bed.output.tsv.gz$", "", basename)

summarized_table=paste0(basename,"_summarized_output.tsv")
filtered_table=paste0(basename,"_filtered_TAPS.tsv")

##############

#calculate Zscore to normalize Ref and input samples:
input_file <- input_file  %>%  mutate(CPM_ref = (input_file$ref / ref_sum)*10^6)
input_file <- input_file  %>%  mutate(CPM_input = (input_file$input / input_sum)*10^6)

##### MAIN TEST ######

test_outoput <- input_file %>%  group_by(TE_ID) %>%
    do(
      w = wilcox.test(x=.$CPM_input, y=.$CPM_ref,
                      paired=FALSE, exact = FALSE, conf.int = TRUE, alternative = "less")
      ) %>%
      summarise(TE_ID, Wpval = w$p.value, West = w$estimate, Wconf = w$conf.int)


final_output <- as.data.frame(dplyr::full_join(input_file, test_outoput, by="TE_ID"))

#summarize for each individual and TE.

input_file_summarized <- final_output %>%  group_by(TE_ID)  %>%
  summarise(Wpval = unique(Wpval),
            median_coverage = median(input),
            median_CPM = median(CPM_input),
            median_CPM_ref=median(CPM_ref)) %>%
  mutate(padjust_FDR = p.adjust(Wpval, method = "fdr", n = length(Wpval)))

#MERGE TE annotation and   summarized file 
input_file)summarized <- dplyr::left_join(annotation, input_file_summarized, by="TE_ID")

#Write out output without filtering:
write.table(input_file_summarized , paste0("./",summarized_table), quote=F, sep='\t', row.names=F, col.names=T)

#Filter putative absences by pvale and  10 times less than mean coverage 
input_file_filtered <- input_file  %>%  
  filter(padjust_FDR < 0.05) %>%  
  filter(mean_coverage < median(input_sum)/10)  

#Write out filterined: 
write.table(input_file_filtered, paste0(path,"/",filtered_table), quote=F, sep='\t', row.names=F, col.names=T)
