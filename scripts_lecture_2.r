#
# Lecture_2.r
#
# Created by David Ruau on 2011-04-26.  
# Copyright (c) 2011 Dept. of Pediatrics and Anesthesia. All rights reserved.

## ACCESS R-CLOUD
# OPEN THE TUTORIAL PROJECT.

#####################################################
## PREREQUISITE PACKAGE TO INSTALL
#####################################################
source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")

## load the required library
library(GEOquery)
library(affy)

#####################################################
## OBTAIN RAW DATA FROM GEO
#####################################################
# GSE12499 raw data set with 3 cell types
getGEOSuppFiles("GSE12499")

# uncompress the archive and clean up to keep only the essential
system('tar -xf GSE12499/GSE12499_RAW.tar -C GSE12499/')
system('rm GSE12499/*.CHP*')
system('rm GSE12499/*.tar')

#####################################################
## IMPORT THE AFFYMETRIX DATA IN R
#####################################################
da <- ReadAffy(celfile.path="./GSE12499/", compress=TRUE)

#####################################################
## ADDING SOME PHENODATA
#####################################################
# what are the phenoData by default
pData(da)
# 
sampleNames(da)

# upload the new phenodata with the one from the text file
# upload the treatment.txt file in R
URL <- "http://www.stanford.edu/~druau/treatment.txt"
download.file(URL, "./treatment.txt")
pd <- read.table("treatment.txt", sep='\t', header=TRUE)
pd
# comma separated files
URL <- "http://www.stanford.edu/~druau/treatment.csv"
download.file(URL, "./treatment.csv")
pd <- read.table("treatment.csv", sep=',', header=TRUE)
pd

## update teh phenoData of your data set
pData(da) <- pd
pData(da)

#####################################################
## BACKGROUND CORRECT + NORMALIZE YOUR DATA USING RMA
#####################################################
library(gcrma)
da.n <- rma(da)
dim(da.n)

#####################################################
## EXTRACT THE NORMALIZED EXPRESSION VALUES
#####################################################
da.eset <- exprs(da.n)

# look what is inside this normalized expression matrix
dim(da.eset)
colnames(da.eset)
# [1] "NSC_1F_iPS_1" "NSC_1F_iPS_2" "NSC_1F_iPS_3" "1F_iPS_NSC_1" "1F_iPS_NSC_2"
#  [6] "1F_iPS_NSC_3" "NSC_1"        "NSC_2"        "NSC_3"        "NSC_4"

#####################################################
## EXTRACT THE GENE SIGNIFICANTLY REGULATED USING RANK PRODUCT
#####################################################
# loading the Rank Product library
library(RankProd)
# Loading the annotation for the microarray
library(mouse4302.db)
gnames <- as.vector(unlist(as.list(mouse4302SYMBOL)))

# build a vector of the condition you want to compare
# here we will compare NSC to 1F iPS cells
cl <- c(rep(0,3), rep(1,4))
da.rp <- RP(da.eset[,c(7:10, 4:6)], cl=cl, logged=TRUE, num.perm=100, plot=FALSE, rand=5432)

# now we extract the genes with corrected p-value above 0.05
r.nsc.1fipsc <- topGene(da.rp, cutoff=0.05, method="pfp",logged=TRUE,logbase=2, gene.names=gnames)

#####################################################
## THE LISTS
#####################################################
# The genes significantly up-regulated
head(r.nsc.1fipsc$Table1, 20)
# how many genes are in table 2
dim(r.nsc.1fipsc$Table1)

# The genes significantly down-regulated
head(r.nsc.1fipsc$Table2, 20)
# how many genes are in table 2
dim(r.nsc.1fipsc$Table2)

#####################################################
## EXPORT THE LIST TO A TABULATED FILE
#####################################################
# up reg
x <- r.nsc.1fipsc$Table1
# replace the fold change value by their log2 counterpart
x[,3] <- log2(1/x[,3])
write.table(x, file = "NSC_vs_1F_iPSC_up.txt", sep = "\t", quote = FALSE, row.names=TRUE)

# down reg
x <- r.nsc.1fipsc$Table2
# replace the fold change value by their log2 counterpart
x[,3] <- log2(1/x[,3])
write.table(x, file = "NSC_vs_1F_iPSC_down.txt", sep = "\t", quote = FALSE, row.names=TRUE)

#####################################################
## SAVE YOUR WORKSPACE
#####################################################
save.image()

