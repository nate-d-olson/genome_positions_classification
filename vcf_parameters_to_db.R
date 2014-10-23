## Script for loading vcf files into sqlite db - changed to postgres
library(VariantAnnotation) #
library(data.table) 
library(plyr) 
library(stringr) 
library(dplyr) 
#library(tidyr)
#library(reshape2)

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
  full_split <- str_split(vcf_filename,pattern = "mpileup_vcf//")[[1]]
  if(grepl("MiSeq",full_split[1])){
    sample_name <- str_split(full_split[2],"_")[[1]]
    return(c("name" = sample_name[1],"PLAT"= "MiSeq","VIAL" = str_sub(sample_name[1],2,2) , "REP" = str_sub(sample_name[1],5)))
  } else {
    vial = str_sub(full_split[2], 13,13)
    return(c("name" = str_c("PGM",vial,sep = "-"),"PLAT"= "PGM","VIAL" = vial, "REP" = 1))
  }  
}

## calculating position probabilites and purity
process_vcf_purity <- function (vcf_file, vcf_db){
  # get metadata
  vcf_meta <- parse_vcf_filename(vcf_file)
  tbl_name <- str_replace(string = unname(vcf_meta["name"]),pattern = "-",replacement = "_")
  #if(tbl_name %in% dbListTables(vcf_db$con)){ # can potentially replace dbListTables with src_tbls
   # return("Next Dataset")
  #}
  
  #read vcf
  vcf <- readVcf(vcf_file, geno=ref)
  

#   #get I16 info into a data.table  
    I16_names <- c("R_Q13_F","R_Q13_R","NR_Q13_F","NR_Q13_R","RS_BQ",
                   "R_BQ_SSq","NR_BQ_S","NR_BQ_SSq","R_MQ_S","R_MQ_SSq",
                   "NR_MQ_S","NR_MQ_SSq","R_TD_S","R_TD_SSq","NR_TD_S","NR_TD_SSq")
    I16 <- ldply(info(vcf)$I16) %>% tbl_df() %>% setnames(I16_names)
#   
#   # calculate purity
   PUR <- sapply(info(vcf)$I16,FUN = calc_purity)
   PUR_prob97 <- sapply(info(vcf)$I16,FUN=calc_pure_prob, p = 0.97)
   PUR_Q2.5 <- sapply(info(vcf)$I16,FUN=calc_quantile, p = 0.025)
   PUR_Q50 <- sapply(info(vcf)$I16,FUN=calc_quantile, p = 0.5)
   PUR_Q97.5 <- sapply(info(vcf)$I16,FUN=calc_quantile, p = 0.975)
#   
  #generate datatable
  vcf_tbl <- data.table(CHROM = str_sub(string = rownames(info(vcf)),start = 1,end = 8), 
                        POS = ranges(vcf)@start, WIDTH = ranges(vcf)@width, 
                        DP = info(vcf)$DP, QUAL = vcf@fixed$QUAL, PUR, PUR_prob97, PUR_Q2.5, PUR_Q50, PUR_Q97.5)

   I16$POS <- vcf_tbl$POS
   I16$WIDTH <- vcf_tbl$WIDTH
   I16$CHROM <- vcf_tbl$CHROM
#   vcf_join <- join(vcf_tbl,I16)
  copy_to(vcf_db, vcf_tbl, name = tbl_name, temporary = FALSE, indexes = list("CHROM","POS","WIDTH")) 
  copy_to(vcf_db, I16, name = str_c("I16",tbl_name, sep = "_"), temporary = FALSE, indexes = list("CHROM","POS","WIDTH")) 
#  copy_to(vcf_db, vcf_join, name = tbl_name, temporary = FALSE, indexes = list("CHROM","POS", "END_POS"))
  #dbWriteTable(conn = vcf_db$con, value =vcf_join, name = "pruity_test", append =TRUE,)
  rm(vcf,vcf_tbl,I16,PUR,PUR_prob97)
}

#initiate sqlite database
vcf_db <- src_sqlite("../../data/RM8375/vcf_db_split_2.sqlite", create = TRUE)

#processing all mpileup vcf files
vcf_dir_list <- list.files(c("../../data//RM8375//PGM//mpileup/mpileup_vcf/","../../data//RM8375//MiSeq//mpileup/mpileup_vcf/"),full.names = TRUE)
vcf_dir_list <- grep("Undetermined",vcf_dir_list,invert =  TRUE,value = TRUE)
vcf_dir_list <- grep("nomatch",vcf_dir_list,invert = TRUE, TRUE,value = TRUE)

for(vcf in vcf_dir_list){
  process_vcf_purity(vcf_file = vcf, vcf_db)
}
