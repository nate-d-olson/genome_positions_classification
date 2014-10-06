## Script to create a single sqlite data base with purity values for all datasets.

library(plyr)
library(stringr)
library(dplyr)

## get vcf purity parameters --------------------------------------------------
parse_vcf_filename <- function(vcf_filename){
  full_split <- str_split(vcf_filename,pattern = "mpileup_vcf//")[[1]]
  if(grepl("MiSeq",full_split[1])){
    sample_name <- str_split(full_split[2],"_")[[1]]
    return(c("name" = sample_name[1],"PLAT"= "MiSeq",
             "VIAL" = str_sub(sample_name[1],2,2) ,
             "REP" = str_sub(sample_name[1],5,5)))
  } else {
    vial = str_sub(full_split[2], 13,13)
    return(c("name" = str_c("PGM",vial,sep = "-"),
             "PLAT"= "PGM","VIAL" = vial, "REP" = 1))
  }  
}


## db table metadata ----------------------------------------------------------
vcf_dir_list <- list.files(c("../../data//RM8375//PGM//mpileup/mpileup_vcf/","../../data//RM8375//MiSeq//mpileup/mpileup_vcf/"),full.names = TRUE)
vcf_dir_list <- grep("Undetermined",vcf_dir_list,invert =  TRUE,value = TRUE)
vcf_dir_list <- grep("nomatch",vcf_dir_list,invert = TRUE, TRUE,value = TRUE)

vcf_meta_df <- ldply(vcf_dir_list,parse_vcf_filename)
vcf_meta_df$tbl_name <- str_replace(string = vcf_meta_df$name ,pattern = "-",replacement = "_")


## create pure_db -------------------------------------------------------------
vcf_db <- src_sqlite("../../data/RM8375/vcf_db.sqlite3", create = T)

get_pure_quants <- function(db_tbl, src = vcf_db){
  #selects and retrieves the purity quantile values for each PGM datatable
  return(tbl(src = vcf_db, db_tbl) %>% select(POS, CHROM, VIAL,REP, PLAT, PUR_Q50))
}

#will need to exclude undetermined
#purity_tbl <- alply(vcf_meta_df$tbl_name,.fun = get_pure_quants, src = vcf_db,.margins = 1) %>% rbind_all() 
for(vcf_tbl in vcf_meta_df$tbl_name[1:3]){
  print(vcf_tbl)
  copy_to(
          dest = vcf_db, name = "purity_test", 
          temporary = FALSE, 
          indexes = list("PLAT","VIAL","REP","CHROM","POS"),
          append=TRUE)
}
alply(vcf_meta_df$tbl_name[1:3],.fun=get_pure_quants) 
%>% copy_to(dest = vcf_db, name = "purity_test", 
                                                             temporary = FALSE, 
                                                             indexes = list("PLAT","VIAL","REP","CHROM","POS"),
                                                             append=TRUE)

for(vcf_tbl in vcf_meta_df$tbl_name){
  print(vcf_tbl)
  print(tbl(vcf_db,from = vcf_tbl))
}

tbl1 <- vcf_meta_df$tbl_name[1]
tbl2 <- vcf_meta_df$tbl_name[2]
tbl(vcf_db,from = tbl1) %>%
  full_join(tbl2)


get_pure_quants <- function(db_tbl, src = vcf_db){
  #selects and retrieves the purity quantile values for each PGM datatable
  return(tbl(src = vcf_db, db_tbl) %>% select(PUR_Q50))
}
purity <- laply(vcf_meta_df$tbl_name[1:3],.fun=get_pure_quants) %>% rowwise() %>% summarise(median)
