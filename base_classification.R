## Base classification algorithm
library(VariantAnnotation);library(data.table);library(reshape2);library(plyr);library(dplyr);library(knitr);library(stringr);library(ggplot2)

#read dataset
#change to accept command line arguments
setwd("/media//nolson/second//mirror/Micro_RM/data/RM8375/MiSeq/mpileup/")
ref = "/media//nolson/second//mirror/Micro_RM/data/RM8375/ref/CFSAN008157.HGAP.fasta"

#Needed to decompress the vcf file before using readVcf
vcf <- readVcf("S0h-1-test.vcf", geno=ref)


# Calc prob and purity values ------------------------------------

calc_purity <- function(I16){
  return( sum(I16[1:2]) / sum(I16[1:4]) )
}

calc_pure_prob <- function(I16,p){
  return(pbinom(q=sum(I16[1:2]),size = sum(I16[1:4]),prob = p,))
}

calc_poly_prob <- function(I16,p){
  return(pbinom(q=sum(I16[3:4]),size = sum(I16[1:4]),prob = p))
}

## calculating position probabilites and purity
PUR <- sapply(vcf@info@listData$I16,FUN = calc_purity)
HC_prob90 <- sapply(vcf@info@listData$I16,FUN=calc_pure_prob, p = 0.9)
HC_prob99 <- sapply(vcf@info@listData$I16,FUN=calc_pure_prob, p = 0.99)
CP_prob <- sapply(vcf@info@listData$I16,FUN=calc_poly_prob, p = 0.1)


# Variant Annotations Data Table ===============================
annotations <- data.table(CHROM = str_sub(string = rownames(info(vcf)),start = 1,end = 8), POS = ranges(vcf)@start,
                          DP = info(vcf)$DP, QUAL = vcf@fixed$QUAL, RPB = info(vcf)$RPB, MQB = info(vcf)$MQB, 
                          BQB = info(vcf)$BQB, MBSQ = info(vcf)$MQSB, MQ0F = info(vcf)$MQ0F, PUR, HC_prob90, HC_prob99, CP_prob)


## Splitting data to make more manageable size
library(caret)
set.seed(3456)
trainIndex <- createDataPartition(annotations$HC_prob90, p = .25,
                                  list = FALSE,
                                  times = 1)
annoTrain <- annotations[as.vector(trainIndex),]

# Variant Annotations ==========================================
# need to work out methods for defining functions for use with dplyr
base_classification <- function(Bias, HC_prob90, HP_prob99, CPprob){
  if(bias > 3) return(IB)
  if(DP < 10) return("LC")
  if(HC_prop99 == 1) return("3C")
  if(HC_prob90 == 1) return("2C")
  if(HC_prob90 > 0.95) return("1C")
  if(CPprob == 1) return("3P")
  if(CPprob > 0.95 & Bias > 1) return("1P")
  if(CPprob > 0.95) return("2P")
}
  
bias_characterization <- function(RPB,MBQ,BQB,MBSQ){
  bias_count = 0
  for(bias in c(RPB,MBQ,BQB,MBSQ)){
    if(!is.na(bias) & bias > 0.95)
      bias_count = bias_count + 1
  }
  return(bias_count)  
}

annoTrain2 <- annoTrain
#annotations2 <- summarize(annotations, BIAS = bias_characterization, CLASS = base_classification)
#note for trainning dataset that is 1/4 original size BIAS takes ~10 minutes
#does not work return a row of single row dataframe with <fnct[1]>
#annoTrain2 %>% rowwise() %>% do(BIAS = bias_characterization)




## Classification analysis

