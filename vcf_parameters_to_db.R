## Script for loading vcf files into sqlite db
library(VariantAnnotation)
library(data.table)
library(plyr)
library(stringr)
library(dplyr)
library(tidyr)
library(reshape2)

library(doMC)
## Specify the number of cores
registerDoMC(8)

ref = "../../data/RM8375/ref/CFSAN008157.HGAP.fasta"

#per base purity
calc_purity <- function(I16){
  return( sum(I16[1:2]) / sum(I16[1:4]) )
}

#per base purity probabilities
calc_pure_prob <- function(I16,p){
  return(pbinom(q=sum(I16[1:2]),size = sum(I16[1:4]),prob = p,))
}

## calculating position probabilites and purity
process_vcf_purity <- function (vcf_file, plat, vial, rep, vcf_db){
  #read vcf
  vcf <- readVcf(vcf_file, geno=ref)
  
  #get I16 info into a data.table  
  I16_names <- c("R_Q13_F","R_Q13_R","NR_Q13_F","NR_Q13_R","RS_BQ",
                 "R_BQ_SSq","NR_BQ_S","NR_BQ_SSq","R_MQ_S","R_MQ_SSq",
                 "NR_MQ_S","NR_MQ_SSq","R_TD_S","R_TD_SSq","NR_TD_S","NR_TD_SSq")
  I16 <- ldply(info(vcf)$I16,.parallel = TRUE) %>% tbl_df() %>% setnames(I16_names)
  
  # calculate purity
  PUR <- sapply(info(vcf)$I16,FUN = calc_purity)
  PUR_prob97 <- sapply(info(vcf)$I16,FUN=calc_pure_prob, p = 0.97)
  
  #generate datatable
  vcf_tbl <- data.table(PLAT=plat, VIAL=vial,REP=rep, CHROM = str_sub(string = rownames(info(vcf)),start = 1,end = 8), 
                        POS = ranges(vcf)@start, DP = info(vcf)$DP, QUAL = vcf@fixed$QUAL, RPB = info(vcf)$RPB, MQB = info(vcf)$MQB, 
                            BQB = info(vcf)$BQB, MBSQ = info(vcf)$MQSB, MQ0F = info(vcf)$MQ0F, PUR, PUR_prob97)
  I16$POS <- vcf_tbl$POS
  vcf_join <- join(vcf_tbl,I16)
  
  #move to database
  ##%#%#%# NEED TO CHANGE THE NEW TABLE NAME
  copy_to(vcf_db, vcf_join, name = "S1", temporary = FALSE, indexes = list("PLAT","VIAL","REP","CHROM","POS"))
  rm(vcf,vcf_join,vcf_tbl,I16,PUR,PUR_prob97)
}



#initiate sqlite database
vcf_db <- src_sqlite("../../data/RM8375/vcf_db.sqlite3", create = T)

# extrac desired parameters from vcf files and load into sqlite db
# file:///home/nolson/R/x86_64-pc-linux-gnu-library/3.1/dplyr/doc/databases.html - dplyr databases vinette
process_vcf_purity(vcf_file = "../../data/RM8375//MiSeq/mpileup/mpileup_vcf/S0h-1_S1_L001_R1_001.bwa.dedup.vcf",
                             plat = "MiSeq",vial = 0,rep = 1, vcf_db)

#write code for processing all data sets
#need to parse dataset names
#alternatively using abreviated name and generate lookup table for 

#3. summarizing data
#pairwise and "multi"wise plots

