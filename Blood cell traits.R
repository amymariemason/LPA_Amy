###### Blood Cell Traits
# Author: Amy Mason
# Date: Oct 2017
# Goal: Search large data sets of blood cell data for relevant variants, create subset if present
# Inputs: list of variants, data to search for those variants
# Outputs: log file reporting missing variants, dataframe containing the subset of larger data matching results for the variants given 
######

###




# clean workshpace
rm(list = ls())

# libraries
library(data.table)

# input file (large data files)
inputfile<- '//me-filer1/home$/am2609/My Documents/Blood Cell Traits Data/plt_build37_166066_20161212.tsv'

# name output file
outputfile<- '//me-filer1/home$/am2609/My Documents/Blood Cell Traits Data/Output/plt.Rda'

# logfile 
logfile<-'//me-filer1/home$/am2609/My Documents/Blood Cell Traits Data/Logs/plt.log'


# import 43 variants as table
var43<-read.table("//me-filer1/home$/am2609/My Documents/Programs/Amy 1/Data/43var.txt")

# no chr notation in database, so strip those off the front (I'm assuming the chr is implict in the tsv file)
# split into chromosome and position
var43$V2<-substring(var43$V1, 4, 100)
var43$Missing<-rep(0,43)
var43 <- cbind(var43, do.call("rbind", strsplit(var43[, "V2"], ":")))
names(var43)<-c("var", "V2", "missing", "chr", "pos")


# create empty data frame for variant data
output<-fread(file=inputfile,nrows=1)
output$duplicate<-0
output<- output[0,]

# loop to read in data from outside file; assigned as missing added to output file
num_er =0 # error counter
for (j in var43$V2){
  # add error handling
  #   print(j)
  if (exists("test")) rm(test)
  invisible(test <- try(fread(file=inputfile,nrows=5, skip=j)))
  # if error, report variant as missing  
  if("try-error" %in% class(test)){
    var43[var43$V2==j,]$missing<-1;
    num_er <- num_er +1;
  } 
  
  if(!("try-error" %in% class(test))) {
    test$duplicate<-rep(0, nrow(test))
    names(test)<-names(output);
    testSubset <- test[grep(j, test$VARIANT), ];
    #check for multiple matching rows			
    if(nrow(testSubset)!=1) testSubset$duplicate<-rep(1, nrow(testSubset));
    output <- rbind(output, testSubset) 
  }
}



# save output file
save(output, file=outputfile)

# create log file 
sink(logfile, append=FALSE, split=TRUE)
#add comments with cat()
writeLines("This log file is showing which variants are missing from the I73 tsv file")

# create record of missing variants
cat("There are a total of  ", num_er, " variants are missing from the input data \n")
write.table(var43[var43$Missing==1,]$V1, row.names=FALSE, quote=FALSE, col.names=FALSE)

cat("More than one record found for these variants below \n")
write.table(unique(output[output$duplicate==1,]$V1), row.names=FALSE, quote=FALSE, col.names=FALSE)


# end log file
sink()
sink.number()==0

