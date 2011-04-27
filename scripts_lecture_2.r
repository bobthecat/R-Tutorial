#
# Lecture_2.r
#
# Created by David Ruau on 2011-04-26.  
# Copyright (c) 2011 Department of Pediatrics/Cancer Biology Stanford University. All rights reserved.

## access R-cloud
# open the tutorial project.

## Prerequisite package to install
source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")

## load package library
library(GEOquery)

## obtain raw data from GEO
# GSE12499 raw data set with 3 cell types
getGEOSuppFiles("GSE12499")

# uncompress the archive and clean up to keep only the essential
system('tar -xf GSE12499/GSE12499_RAW.tar -C GSE12499/')
system('rm GSE12499/*.CHP*')
system('rm GSE12499/*.tar')

## Import the Affymetrix data in R
da <- ReadAffy(celfile.path="./GSE12499/", compress=TRUE)

## Add some phenoData
sampleNames(da)


