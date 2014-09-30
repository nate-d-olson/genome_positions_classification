library(data.table)
library(plyr)
library(stringr)
library(dplyr)
library(reshape2)

## get vcf purity parameters --------------------------------------------------
parse_vcf_filename <- function(vcf_filename){
  full_split <- str_split(vcf_filename,pattern = "mpileup_vcf/")[[1]]
  if(grepl("MiSeq",full_split[1])){
    sample_name <- str_split(full_split[2],"_")[[1]]
    return(c("name" = sample_name[1],"PLAT"= "MiSeq",
             "VIAL" = str_sub(sample_name[1],2,2) ,
             "REP" = str_sub(sample_name[1],5)))
  } else {
    vial = str_sub(full_split[2], 14,14)
    return(c("name" = str_c("PGM",vial,sep = "-"),
             "PLAT"= "PGM","VIAL" = vial, "REP" = 1))
  }  
}


## db table metadata ----------------------------------------------------------
vcf_dir_list <- list("../../data//RM8375//PGM//mpileup/mpileup_vcf/","../../data//RM8375//MiSeq//mpileup/mpileup_vcf/")
vcf_dir_list <- grep("nomatch\|Undetermined",vcf_dir_list,inverse = TRUE,value = TRUE)

vcf_meta_df <- dlply(vcf_dir_list,parse_vcf_filename)
vcf_meta_df$tbl_name <- str_replace(string = vcf_meta_df$name ,pattern = "-",replacement = "_")


## create pure_db -------------------------------------------------------------
vcf_db <- src_sqlite("../../data/RM8375/vcf_db.sqlite3", create = T)

get_pure_quants <- function(db_tbl){
  #selects and retrieves the purity quantile values for each PGM datatable
  return(data.table(tbl(vcf_db(), db_tbl) %>% select(POS, CHROM, VIAL,REP, PLAT, PUR_Q2.5, PUR_Q50, PUR_Q97.5))
}

#will need to exclude undetermined
purity_tbl <- alply(vcf_meta_df$tbl_name,.fun = get_pure_quants) %>% rbind_all() 
copy_to(vcf_db, pruity_tbl, name = purity, temporary = FALSE, indexes = list("PLAT","VIAL","REP","CHROM","POS"),)


