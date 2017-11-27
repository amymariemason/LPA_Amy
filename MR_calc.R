######
# Author: Amy Mason ; based heavily on code by James Staley  
# Date: Oct 2017
# Goal: Perform Mendelian randomization  
# Inputs: list of variants, data to search for those variants
# Outputs: log file reporting missing variants, dataframe containing the subset of larger data matching the variants given 
######

#####################################################
##### Set-up #####
#####################################################
#Run for
# 1) Platelet count

rm(list=ls())
setwd("C://Users/am2609/Dropbox (Personal)/lpamaster/")
library(ggplot2)
library(MendelianRandomization)
library(plotly)
library(htmlwidgets)
library(assertthat)
library(stringr)

# outcome to test
outcomespec<-"plt"
outcomename<-"platelet count"

# create log file 
logfile<-paste ("\\\\me-filer1/home$/am2609/My Documents/Programs/Amy 1/Logs/mr_working_", outcomespec, ".log", sep = "")

sink(logfile, append=FALSE, split=TRUE)
#add comments with cat()
cat("\n This log file is showing working on the merge and mr calculations with ", outcomespec, " file \n")

# load data
loadplace<- paste ("//me-filer1/home$/am2609/My Documents/Blood Cell Traits Data/Output/", outcomespec, ".Rda", sep = "")
load(loadplace)
outcome<-output

data_rho_master <- read.table("lpa/data/LPA_master_dataset_pcs_EUwinsor_withoutM.txt", header=T, sep="\t", colClasses="character")

ivw_corr <- function(bx, by, byse, rho){
  Omega <- byse%o%byse*rho
  beta <- solve(t(bx)%*%solve(Omega)%*%bx)*t(bx)%*%solve(Omega)%*%by
  se <- sqrt(solve(t(bx)%*%solve(Omega)%*%bx))
  df <- length(bx) - 1
  lci <- beta - qt(df=df, 0.975)*se
  uci <- beta + qt(df=df, 0.975)*se
  p <- 2*(1-pt(abs(beta/se),df))
  results <- data.frame(beta,se,lci,uci,p)
  return(results)
}

## External

lpa <- read.table("lpa/data/LPA_Variants_EUwinsor_withoutM.txt", sep="\t", header=T, colClasses="character")
# lpa <- read.table("lpa/data/LPA_Variants_EUnowinsor_withoutM.txt", sep="\t", header=T, colClasses="character")
# lpa <- read.table("lpa/data/LPA_Variants_EUquantile_withoutM.txt", sep="\t", header=T, colClasses="character")
# lpa <- read.table("lpa/data/LPA_Variants_EUrobust_withoutM.txt", sep="\t", header=T, colClasses="character")
fsteps <- read.table("lpa/data/fstep_snps_0.4_EUwinsor.txt", sep="\t", header=T, colClasses="character")
# fsteps = list(snp = c("rs10455872_g", "rs3798220_t"))
lpa <- lpa[(lpa$variantID %in% fsteps$snp),]
lpa$chr.pos <- paste0("chr", lpa$chr, ":", lpa$pos)
lpa$snp <- lpa$chr.pos; lpa$a1 <- lpa$allele1; lpa$a2 <- lpa$allele2
lpa <- lpa[, c("variantID", "snp", "chr.pos", "chr", "pos" , "a1", "a2", "beta", "se")]
#Steve says: A1 is the effect allele

# CHECK: I am assuming these are from positive strand - is this true?



# list of fields 
# : ref allele: alternate allele
#VARIANT (hg19) [CHROM:POS:REF:ALT]	Variant position as [Chromosome : hg19 Position : reference allele : alternate allele]
#rsid	SNP rsID as provided by UK Biobank
#nCompleteSamples	Number of samples analyzed with non-missing phenotypes
#AC	Dosage allele count from all samples with non-missing phenotypes
#ytx	Dosage alternate allele count (ytx means "y * x" where y=phenotype and x=alternate allele dosage) 
#beta	linear regression beta coefficient
#se	linear regression standard error
#pval	linear regression p-value

#steve says: reference allele is the effect allele


# change so everything is positive strand

temp<-outcome
temp$REF2<-temp$REF


temp[temp$DIRECTION=="-"&temp$REF=="G",]$REF2<-str_replace(temp[temp$DIRECTION=="-"&temp$REF=="G",]$REF, "G", "C" )
temp[temp$DIRECTION=="-"&temp$REF=="C",]$REF2<-str_replace(temp[temp$DIRECTION=="-"&temp$REF=="C",]$REF, "C", "G" )
temp[temp$DIRECTION=="-"&temp$REF=="A",]$REF2<-str_replace(temp[temp$DIRECTION=="-"&temp$REF=="A",]$REF, "A", "T" )
temp[temp$DIRECTION=="-"&temp$REF=="T",]$REF2<-str_replace(temp[temp$DIRECTION=="-"&temp$REF=="T",]$REF, "T", "A" )
is.flag(nrow(temp[temp$DIRECTION=="-"& nchar(temp$REF)!=1,])==0)

# all variables only single value, so replace
temp[temp$DIRECTION=="-"]$REF<-temp[temp$DIRECTION=="-"]$REF2

# now do the ALT value column
temp$ALT2<-temp$ALT
temp[temp$DIRECTION=="-"&temp$ALT=="G",]$ALT2<-str_replace(temp[temp$DIRECTION=="-"&temp$ALT=="G",]$ALT, "G", "C" )
temp[temp$DIRECTION=="-"&temp$ALT=="C",]$ALT2<-str_replace(temp[temp$DIRECTION=="-"&temp$ALT=="C",]$ALT, "C", "G" )
temp[temp$DIRECTION=="-"&temp$ALT=="A",]$ALT2<-str_replace(temp[temp$DIRECTION=="-"&temp$ALT=="A",]$ALT, "A", "T" )
temp[temp$DIRECTION=="-"&temp$ALT=="T",]$ALT2<-str_replace(temp[temp$DIRECTION=="-"&temp$ALT=="T",]$ALT, "T", "A" )
is.flag(nrow(temp[temp$DIRECTION=="-"& nchar(temp$ALT)!=1,])==0)

# all variables only single value, so replace
temp[temp$DIRECTION=="-"]$ALT<-temp[temp$DIRECTION=="-"]$ALT2


# all the negative strands have been replaced now so
temp[temp$DIRECTION=="-"]$DIRECTION<-"+"

# I'm happy this is correct, so use as outcome
#outcome2<-temp[,c("CHR", "BP", "REF", "ALT", "EFFECT", "SE")]

# this added because not matching strand pos/neg with steven's dataset
outcome2<-outcome[,c("CHR", "BP", "REF", "ALT", "EFFECT", "SE")]

names(outcome2)<-c("chr","pos","A1","A2", "outcomeBeta", "outcomeSE")
outcome2$chr<-as.character(outcome2$chr)
outcome2$pos<-as.character(outcome2$pos)

# change to D/I format for insertions to match lpa data set
outcome2[outcome2$pos=="160899049" & outcome2$A2=="CA","A1"] <- "D" 
outcome2[outcome2$pos=="160899049" & outcome2$A2=="CA","A2"] <- "I" 


# attempt to merge with lpa dataset to highlight any alleles in wrong direction

# subset lpa to those variants in the outcome dataset
outcome3<-merge(outcome2, lpa, by=c("chr", "pos"))

#WARNING - alignment issues ????

# check not lost any variants from outcome
cat("\nThis shows an error if variants missing compared to outcome \n")
tryCatch({stopifnot(nrow(outcome3) == nrow(outcome2))
}, error = function(err.msg){
  # Add error message to the error log file
  cat("/n", toString(err.msg), "/n")
}
)

cat("\nThis shows an error if variants missing compared to lpa \n")
tryCatch({stopifnot(nrow(outcome3) == nrow(lpa))
}, error = function(err.msg){
  # Add error message to the error log file
  cat("\n", toString(err.msg), "\n", nrow(lpa)-nrow(outcome3), "variants missing from outcome data\n")
}
)


# create correlation matrix

data_rho <- data_rho_master[(data_rho_master$prev_chd==0), ]
data_rho <- data_rho[, match(outcome3$variantID, names(data_rho))]
data_rho <- as.matrix(data_rho)
class(data_rho) <- "numeric"
rho <- cor(data_rho, use="complete.obs")


# Set-up for mr
bx_effect <- outcome3$a1
by_effect <- outcome3$A1
bx <- as.numeric(outcome3$beta)
by <- as.numeric(outcome3$outcomeBeta)
by = ifelse(bx_effect == by_effect, by, -by)
byse <- as.numeric(outcome3$outcomeSE)
bxse <- as.numeric(outcome3$se)
rho <- as.matrix(rho)

# Analysis
MRdata_input <- mr_input(bx, bxse, by, byse, corr=rho, outcome=outcomename, exposure="Lp(a)", snps=outcome3$snp)

mr = mr_ivw(mr_input(bx, bxse, by, byse, corr=rho))

mr

# attempt at a graph
mr_plot(MRdata_input, interactive=TRUE)

sink()
