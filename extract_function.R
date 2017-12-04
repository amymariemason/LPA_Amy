###### Extract function
# Author: Amy Mason
# Date: Oct 2017
# Goal: Search large data sets of blood cell data for relevant variants, create subset if present
# Inputs: list of variants, data to search for those variants
# Outputs: log file reporting missing variants, dataframe containing the subset of larger data matching results for the variants given
######



# clean workshpace
rm(list = ls())

# libraries
#install.packages('data.table',lib='C://Users/am2609/R/library')
library(data.table,lib.loc='C://Users/am2609/R/library')
#library(data.table)

# input file (large data files)
#inputfile<-'Blood Cell Traits Data/plt_build37_166066_20161212.tsv'



####################################
# Set folder for outputs
#####################################
#setwd("//me-filer1/home$/am2609/My Documents/")
# setwd("Z://My Documents/") 
setwd("C://Users/am2609/")

outputfolder<- './Blood Cell Traits Data/'

##################################
#data input
##################################

# import list of traits to include
TL<-scan("./Blood Cell Traits Data/Traitlist.txt", what="", sep="\n")
#end of trait file
fileend<-'_20161212.tsv'
# set file to output to
outputfolder<- './Blood Cell Traits Data/'


####################################
#variant input
####################################
#varfile<-"./Programs/Amy 1/Data/43var.txt"
varfile<-"./Blood Cell Traits Data/43var.txt"
# varfile should have format: single column "chr" + chromosome number + ":" + chromosome position
# import 43 variants as table
var43<-read.table(varfile)

# no chr notation in database, so strip those off the front (I'm assuming the chr is implict in the tsv file)
# split into chromosome and position

var43$V2<-substring(var43$V1, 4, 100)
var43$Missing<-rep(0,43)
var43 <- cbind(var43, do.call("rbind", strsplit(var43[, "V2"], ":")))
names(var43)<-c("var", "V2", "missing", "chr", "pos")
var43_refresh<-var43


##########################################
# Functions
#########################################

# create function to loop round variables


create_blank_file<- function(trait, inputfile, outputfolder){

# create empty data frame for variant data
output<-fread(file=inputfile,nrows=1)
output$duplicate<-0
output<- output[0,]

return(output)
}

import_data<-function(inputfile, output){

# loop to read in data from outside file; assigned as missing added to output file
num_er =0 # error counter
for (j in var43[var43$missing!=1,]$V2){
#for (j in c("6:160998199","6:161162290")){
  # add error handling
     print(j)
  if (exists("test")) rm(test)
  test <- tryCatch(
    # this is what I want it to do:  
    fread(file=inputfile,nrows=5, skip=j)
    ,
    # if error occurs
    error=function(error_message) {
      message("Error")
      message(error_message)
      return("MISSING")
    }
  )
  # if error, report variant as missing
  if(is.character(test)){
    var43[var43$V2==j,]$missing<-1;
    num_er <- num_er +1;
  }
  if(!(is.character(test))) {
    test$duplicate<-rep(0, nrow(test))
    names(test)<-names(output);
    testSubset <- test[grep(j, test$VARIANT), ];
    #check for multiple matching rows			
    if(nrow(testSubset)!=1) testSubset$duplicate<-rep(1, nrow(testSubset));
    output <- rbind(output, testSubset)
  }
}  
outlist<-list(num_er,var43, output)
return(outlist)
}

save_output_log<-function(output, outputfile, logfile){

# save output file
save(output, file=outputfile)

# create log file
sink(logfile, append=FALSE, split=TRUE)
#add comments with cat()
writeLines("This log file is showing which variants are missing from the I73 tsv file")

# create record of missing variants
cat("There are a total of  ", num_error, " variants are missing from the input data \n")
write.table(varerror[varerror$missing==1, "var"], row.names=FALSE, quote=FALSE, col.names=FALSE)

cat("More than one record found for these variants below \n")
write.table(unique(output[output$duplicate==1,"VARIANT"]), row.names=FALSE, quote=FALSE, col.names=FALSE)


# end log file
sink()
sink.number()==0

}

############################################
 # create loop round data
###########################################


var43<-var43_refresh
missinglog<-paste (outputfolder, "Logs/missingtraits.log", sep="")
 for (i in 1:length(TL)){
 trait<-TL[i]
 traitmin<- substring(trait,1,10)
 print(traitmin)
#var43<-var43_refresh
   # name input file
 inputfile<-paste(outputfolder, trait, '_20161212.tsv', sep="")
 if (file.exists(inputfile)){
 # name output file
 outputfile<- paste (outputfolder, "Output/", trait,".Rda", sep="")
 # name logfile
 logfile<-paste (outputfolder, "Logs/", trait,".log", sep="")
 #create blank output
 output<-create_blank_file(trait,inputfile = inputfile, outputfolder=outputfolder)
 #fill output
 keep<-import_data(inputfile, output)
 num_error<-keep[[1]]
 varerror<-keep[[2]]
 output<-keep[[3]]
 # save output
 save_output_log(output, outputfile, logfile)
 }
 if (!(file.exists(inputfile))){
#add to missing file if file cannot be found
   sink(missinglog, append=TRUE, split=TRUE)  
cat(trait, "\n")
 # end log file
 sink()
 sink.number()==0
  }
}
 

######################################################## 
 