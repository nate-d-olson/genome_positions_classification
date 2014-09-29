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

#global variables
ref = "../../data/RM8375/ref/CFSAN008157.HGAP.fasta"

#per base purity
calc_purity <- function(I16){
  return( sum(I16[1:2]) / sum(I16[1:4]) )
}

#per base purity probabilities
calc_pure_prob <- function(I16,p){
  return(pbinom(q=sum(I16[1:2]),size = sum(I16[1:4]),prob = p,))
}

#per base purity quantiles
calc_pure_prob <- function(I16,p){
  return(pbinom(q=sum(I16[1:2]),size = sum(I16[1:4]),prob = p,))
}

#calculate purity quantiles 
calc_quantile <- function(I16,p){
  return(qbinom(p,size = sum(I16[1:4]), sum(I16[1:2]) / sum(I16[1:4]))/sum(I16[1:04]))
}

## Get file info
parse_vcf_filename <- function(vcf_filename){
  full_split <- str_split(vcf_filename,pattern = "mpileup_vcf/")[[1]]
  if(grepl("MiSeq",full_split[1])){
    sample_name <- str_split(full_split[2],"_")[[1]]
    return(c("name" = sample_name[1],"PLAT"= "MiSeq","VIAL" = str_sub(sample_name[1],2,2) , "REP" = str_sub(sample_name[1],5)))
  } else {
    vial = str_sub(full_split[2], 14,14)
    return(c("name" = str_c("PGM",vial,sep = "-"),"PLAT"= "PGM","VIAL" = vial, "REP" = 1))
  }  
}

## calculating position probabilites and purity
process_vcf_purity <- function (vcf_file, vcf_db){
  #read vcf
  vcf <- readVcf(vcf_file, geno=ref)
  
  # get metadata
  vcf_meta <- parse_vcf_filename(vcf_file)
  #get I16 info into a data.table  
  I16_names <- c("R_Q13_F","R_Q13_R","NR_Q13_F","NR_Q13_R","RS_BQ",
                 "R_BQ_SSq","NR_BQ_S","NR_BQ_SSq","R_MQ_S","R_MQ_SSq",
                 "NR_MQ_S","NR_MQ_SSq","R_TD_S","R_TD_SSq","NR_TD_S","NR_TD_SSq")
  I16 <- ldply(info(vcf)$I16,.parallel = TRUE) %>% tbl_df() %>% setnames(I16_names)
  
  # calculate purity
  PUR <- sapply(info(vcf)$I16,FUN = calc_purity)
  PUR_prob97 <- sapply(info(vcf)$I16,FUN=calc_pure_prob, p = 0.97)
  PUR_Q2.5 <- sapply(info(vcf)$I16,FUN=calc_quantile, p = 0.025)
  PUR_Q50 <- sapply(info(vcf)$I16,FUN=calc_quantile, p = 0.5)
  PUR_Q97.5 <- sapply(info(vcf)$I16,FUN=calc_quantile, p = 0.975)
  
  #generate datatable
  vcf_tbl <- data.table(PLAT=vcf_meta["PLAT"], VIAL=vcf_meta["VIAL"],REP=vcf_meta["REP"], CHROM = str_sub(string = rownames(info(vcf)),start = 1,end = 8), 
                        POS = ranges(vcf)@start, DP = info(vcf)$DP, QUAL = vcf@fixed$QUAL, RPB = info(vcf)$RPB, MQB = info(vcf)$MQB, 
                            BQB = info(vcf)$BQB, MBSQ = info(vcf)$MQSB, MQ0F = info(vcf)$MQ0F, PUR, PUR_prob97, str_c("PUR_Q", c(2.5,50,97.5),sep = ""))
  I16$POS <- vcf_tbl$POS
  vcf_join <- join(vcf_tbl,I16)
  
  #move to database
  copy_to(vcf_db, vcf_join, name = vcf_meta["name"], temporary = FALSE, indexes = list("PLAT","VIAL","REP","CHROM","POS"))
  rm(vcf,vcf_join,vcf_tbl,I16,PUR,PUR_prob97)
}



#initiate sqlite database
vcf_db <- src_sqlite("../../data/RM8375/vcf_db.sqlite3", create = T)

#processing all mpileup vcf files
vcf_files <- list.files(path = "../../data/RM375/*/mipleup/mpileup_vcf/", pattern = "*.vcf", full.names = TRUE)

for(vcf in vcf_files){
  process_vcf_purity(vcf_file = vcf, vcf_db)
}
