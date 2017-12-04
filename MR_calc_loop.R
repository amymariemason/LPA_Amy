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


rm(list=ls())


#install.packages('devtools',lib='C://Users/am2609/R/library')
library(digest,lib.loc='C://Users/am2609/R/library')
library(mime,lib.loc='C://Users/am2609/R/library')
library(ggplot2, lib.loc='C://Users/am2609/R/library')
library(MendelianRandomization, lib.loc='C://Users/am2609/R/library')
library(plotly, lib.loc='C://Users/am2609/R/library')
library(htmlwidgets, lib.loc='C://Users/am2609/R/library')
library(assertthat, lib.loc='C://Users/am2609/R/library')
library(stringr, lib.loc='C://Users/am2609/R/library')

#ggplot2 in devtools is struggling to find local packages; load explicity
library(labeling, lib.loc='C://Users/am2609/R/library')
library(tidyselect, lib.loc='C://Users/am2609/R/library')
library(crosstalk, lib.loc='C://Users/am2609/R/library')
library(yaml, lib.loc='C://Users/am2609/R/library')

#in order to run Rmarkdown
library(tools)
library(utils)
library(methods)
library(knitr, lib.loc='C://Users/am2609/R/library')
library(evaluate, lib.loc='C://Users/am2609/R/library')
library(base64enc, lib.loc='C://Users/am2609/R/library')
library(rprojroot, lib.loc='C://Users/am2609/R/library')
library(highr, lib.loc='C://Users/am2609/R/library')
library(markdown, lib.loc='C://Users/am2609/R/library')
library(caTools, lib.loc='C://Users/am2609/R/library')
library(htmltools, lib.loc='C://Users/am2609/R/library')

#set working directory
#setwd("C://Users/am2609/Dropbox (Personal)/lpamaster/")
setwd("C://Users/am2609/Blood Cell Traits Data/")

# load LPA data
data_rho_master <- read.table("LPA/LPA_master_dataset_pcs_EUwinsor_withoutM.txt", header=T, sep="\t", colClasses="character")

# load lpa effects

lpa <- read.table("LPA/LPA_Variants_EUwinsor_withoutM.txt", sep="\t", header=T, colClasses="character")
fsteps <- read.table("LPA/fstep_snps_0.4_EUwinsor.txt", sep="\t", header=T, colClasses="character")
lpa <- lpa[(lpa$variantID %in% fsteps$snp),]
lpa$chr.pos <- paste0("chr", lpa$chr, ":", lpa$pos)
lpa$snp <- lpa$chr.pos; lpa$a1 <- lpa$allele1; lpa$a2 <- lpa$allele2
lpa <- lpa[, c("variantID", "snp", "chr.pos", "chr", "pos" , "a1", "a2", "beta", "se")]
#Steve says: A1 is the effect allele


# setup testing loop
#1 import data
#2 setup mr comparision
#3 output report

# import list of traits to include
TL<-scan("./Traitlist.txt", what="", sep="\n")
TLframe<-as.data.frame(TL)
names(TLframe)<-"trait"
TLframe$traitmin<- substring(TLframe$trait,1,10)

# load names for traits 
traitnames<-read.csv("./Traitnames.csv", header=FALSE)
traitnames<-traitnames[,1:3]
names(traitnames)<-c("trait", "name", "units")
TLall<-merge(TLframe, traitnames, by="trait", all.x=TRUE)

# create loop

for (i in 1:length(TL)){
  # print outcome name
  outcomespec<-as.character(TLall[i,"trait"])
  outcomename<-as.character(TLall[i,"name"])
  outcomeunit<-as.character(TLall[i,"units"])
  print(as.character(outcomespec))
  #add report to new logfile
  # create log file 
  logfile<-paste ("./Logs/", outcomespec,"_mr.log", sep="")
  sink(logfile, append=FALSE, split=TRUE)
  #add comments with cat()
  cat("\n This log file is showing working on the merge and mr calculations with ", outcomespec, " file \n")
# loop working starts here
    # load data as file called output (see bottom for file details)
      load(paste ("./Output/", outcomespec, ".Rda", sep = ""))
      outcome<-output[,c("CHR", "BP", "REF", "ALT", "EFFECT", "SE")]
      names(outcome)<-c("chr","pos","A1","A2", "outcomeBeta", "outcomeSE")
      outcome$chr<-as.character(outcome$chr)
      outcome$pos<-as.character(outcome$pos)
    # change to D/I format for insertions to match lpa data set
      outcome[outcome$pos=="160899049" & outcome$A2=="CA","A1"] <- "D" 
      outcome[outcome$pos=="160899049" & outcome$A2=="CA","A2"] <- "I" 
    
    # subset lpa to those variants in the outcome dataset
      outcome3<-merge(outcome, lpa, by=c("chr", "pos"))
    # check not lost any variants from outcome/lpa; report to log if so
      cat("\nThis shows an error if variants missing compared to outcome \n")
      tryCatch({stopifnot(nrow(outcome3) == nrow(outcome2))
      }, error = function(err.msg){
      # Add error message to the error log file
      cat("ERROR",  "\n", nrow(outcome3)- nrow(lpa), "variants missing from outcome data\n")
      write.table(outcome3[!(outcome3$pos%in%lpa$pos),c("chr", "pos")], row.names=FALSE, quote=FALSE)
      }
      )
      #
      cat("\nThis shows an error if variants missing compared to lpa \n")
      tryCatch({stopifnot(nrow(outcome3) == nrow(lpa))
      }, error = function(err.msg){
      # Add error message to the error log file
      cat("ERROR", "\n", nrow(lpa)-nrow(outcome3), "variants missing from outcome data\n")
      write.table(lpa[!(lpa$pos%in%outcome3$pos),c("chr", "pos")], row.names=FALSE, quote=FALSE)
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
      
   
     tryCatch({mr = mr_ivw(mr_input(bx, bxse, by, byse, corr=rho))
      }, 
       error = function(error_message){
               cat("ERROR", "\n", outcomespec, " did not run correctly in mr_ivw\n")
              message(error_message)
                },
       warning=function(warn){
              cat("WARNING\n", "outcomespec gave warning in mr_ivw: \n")
              message(warn)
                }, 
       finally ={
              mr = mr_ivw(mr_input(bx, bxse, by, byse, corr=rho))
                }
            )
      # if ran correctly, add to report 
        mr
        # attempt at a graph
        mr_plot(MRdata_input, interactive=TRUE)
     
      
  
  
# end log file
  sink()
  sink.number()==0
}
  


## External





# list of fields in output file
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