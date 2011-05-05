#
# scripts_lecture_3.r
#
# Created by David Ruau on 2011-05-04.  
# Copyright (c) 2011 Dept. of Pediatrics and Anesthesia. All rights reserved.

## ACCESS R-CLOUD
# OPEN THE TUTORIAL PROJECT.

## LOAD THE LIBRARY
library(RankProd)
library(affy)
colnames(da.eset)

## COMPARE NSC TO 1F_NSC
cl <- c(rep(0,4), rep(1,3))
da.rp <- RP(da.eset[,c(7:10, 1:3)], cl=cl, logged=TRUE, num.perm=100, plot=FALSE, rand=5432)

## ANNOTATION PACKAGE
library(mouse4302.db)
gnames <- as.vector(unlist(as.list(mouse4302SYMBOL)))

#####################################################
## EXTRACTING THE SIGNIFICANTLY REGULATED GENES
#####################################################
r.nsc.nsc_1f <- topGene(da.rp, cutoff = 0.05, method = "pfp", logged = TRUE, 
	logbase = 2, gene.names=gnames)
# The genes significantly up-regulated
head(r.nsc.nsc_1f$Table1, 20)
# how many genes are in table 2
dim(r.nsc.nsc_1f$Table1)

# The genes significantly down-regulated
head(r.nsc.nsc_1f$Table2, 20)
# how many genes are in table 2
dim(r.nsc.nsc_1f$Table2)

#####################################################
## COMPUTING THE CORRELATION MATRIX
#####################################################
## We are interested in only a subset of the genes
## dissimilarity
d = as.dist(1 - cor(da.eset, use = "complete.obs", method ="pearson"))
## dendrogram
hc = hclust(d, method = "average")
# pdf(file='dendogram_half_quantile[average].pdf', width=8, height=8)
plot(hc, labels = colnames(da.eset), main = "Hier. clust. Pearson", hang=-1)
# dev.off()

#####################################################
## HEAMAP AND HIERARCHICAL CLUSTERING
#####################################################
# Now we build the dissimilarity matrix
d = as.dist(1 - cor(t(da.eset[r.nsc.nsc_1f$Table1[,1],]), use = "complete.obs", method ="pearson"))
# use revC = TRUE for inverting the heatmap
# Build the color palette
hmcol = colorRampPalette(brewer.pal(10, "RdBu"))(256)
# pdf(file='heatmap_all_array.pdf', height=8, width=8)
heatmap.2(da.eset[r.nsc.nsc_1f$Table1[,1],], 
  distfun=function(x){as.dist(1 - cor(t(x), use = "complete.obs", method ="pearson"))}, 
  hclustfun=function(m){hclust(m, method="average")},
  symm=F, col=hmcol, trace='none', notecol='black', 
  denscol='black', notecex=0.8, dendrogram="column")
# dev.off()






