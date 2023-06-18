#!/usr/bin/env Rscript
rm(list=ls())
library(plyr)
library(dplyr)

##################################  ####################################
reportIdentity <- function (inFile, splID) {

idnData <- read.delim(inFile, header= TRUE, sep="\t")

idnSum <- ddply(idnData, .(ref), summarize, 
				aveIdn=mean(identity),
				avRefseq=mean(refseq))
	
#idn <- idnSum %>% filter(aveIdn >0)
idn <- format(idnSum, digits =1, nsmall =1)

write.table(idn, paste(splID, '.identSummary', sep=""), sep="\t", row.names=FALSE, quote=FALSE)

}

################################## run function ###################################
args = commandArgs(trailingOnly=TRUE)
inFile <- args[1]
splID <- args[2]

reportIdentity(inFile, splID)

