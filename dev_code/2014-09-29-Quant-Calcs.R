## 9/29 read vcf into db methods development


#initiate sqlite database
vcf_db <- src_sqlite("../../data/RM8375/vcf_db.sqlite3", create = T)

vcf_file = "../../data/RM8375//MiSeq/mpileup/mpileup_vcf/S0h-1_S1_L001_R1_001.bwa.dedup.vcf"
ref = "../../data/RM8375/ref/CFSAN008157.HGAP.fasta"

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

#comparison of sapply with dplyr functions

#sapply on CompressedNumericList
system.time(sapply(info(vcf)$I16,FUN = calc_purity)) #28 s
system.time(sapply(info(vcf)$I16,FUN=calc_pure_prob, p = 0.97)) #38 s

#dply on data.table
# was going to take estimates 40 h

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


calc_pur <- function(I16_table){
  I16_table %>% rowwise() %>% do(i = (.$R_Q13_F + .$R_Q13_R)/(.$R_Q13_F + .$R_Q13_R + .$NR_Q13_F + .$NR_Q13_R))
  return(cbind(I16_table,.Last.value %>% summarise(PUR = length(i))))
}
#system.time(sapply(calc_pur(I16)))
#> system.time(sapply(calc_pur(I16)))
#|                                                                                        |  0% ~40 h remaining 


# will use sapply method
#checking multiple values for sapply
calc_quantiles <- function(I16,p){
  quants = qbinom(p,size = sum(I16[1:4]), sum(I16[1:2]) / sum(I16[1:4]))/sum(I16[1:04])
  #names(quants) = c("Q2.5","Q50","Q97.5")
  return(quants)
}

system.time(seqapply(info(vcf)$I16,FUN = calc_quantiles, p = 0.025))
#user  system elapsed 
#63.228   0.164  63.381 

system.time(sapply(info(vcf)$I16,FUN = calc_quantiles, p = 0.025))
#user  system elapsed 
#56.143   0.073  56.198 