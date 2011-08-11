###################################################
### chunk number 1: 
###################################################
#line 71 "R_for_gene_exp.Rnw"
options(width=80, continue=" ")


###################################################
### chunk number 2: loading the libraries
###################################################
#line 83 "R_for_gene_exp.Rnw"
## load package library
library(GEOquery)
library(affy)
library(RColorBrewer)
library(affyPLM)
library(mouse4302cdf)


###################################################
### chunk number 3: importAffy
###################################################
#line 107 "R_for_gene_exp.Rnw"
library(affy)
## Import the Affymetrix data in R
da <- ReadAffy(celfile.path="./GSE12499/", compress=TRUE)


###################################################
### chunk number 4: pdata
###################################################
#line 114 "R_for_gene_exp.Rnw"
## Information about your data
da
## ADDING SOME PHENODATA
# what are the phenoData by default
pData(da)
# 
sampleNames(da)


###################################################
### chunk number 5: import
###################################################
#line 126 "R_for_gene_exp.Rnw"
URL <- "http://www.stanford.edu/~druau/treatment.txt"
pd <- read.table(URL, sep='\t', header=TRUE)
pd


###################################################
### chunk number 6: pdata2
###################################################
#line 133 "R_for_gene_exp.Rnw"
pData(da) <- pd
pData(da)
sampleNames(da) <- pd[,1]


###################################################
### chunk number 7: pset
###################################################
#line 148 "R_for_gene_exp.Rnw"
pset <- fitPLM(da)


###################################################
### chunk number 8: img.test
###################################################
#line 152 "R_for_gene_exp.Rnw"
img.Test <- function(batch,pset,x) {
	par(mfrow = c(2,2))
	image(batch[,x])
	image(pset, type = "weights", which = x)
	image(pset, type = "resids", which = x)
	image(pset, type = "sign.resids", which = x)
}


###################################################
### chunk number 9: fig2plot
###################################################
#line 201 "R_for_gene_exp.Rnw"
cols <- brewer.pal(12, "Set3")
Mbox(pset, col = cols, main ="RLE (Relative Log Expression)", 
	xlab="Assuming that the majority of the gene are not changing\n Ideally these boxes would have small spread and be centered at M=0")


###################################################
### chunk number 10: fig2
###################################################
#line 208 "R_for_gene_exp.Rnw"
#line 201 "R_for_gene_exp.Rnw"
cols <- brewer.pal(12, "Set3")
Mbox(pset, col = cols, main ="RLE (Relative Log Expression)", 
	xlab="Assuming that the majority of the gene are not changing\n Ideally these boxes would have small spread and be centered at M=0")
#line 209 "R_for_gene_exp.Rnw"


###################################################
### chunk number 11: fig3plot
###################################################
#line 219 "R_for_gene_exp.Rnw"
boxplot(pset, col=cols, main= "NUSE (Normalized Unscaled Standard Error)", 
	xlab="High values of median NUSE are indicative of a problematic array")


###################################################
### chunk number 12: fig3
###################################################
#line 225 "R_for_gene_exp.Rnw"
#line 219 "R_for_gene_exp.Rnw"
boxplot(pset, col=cols, main= "NUSE (Normalized Unscaled Standard Error)", 
	xlab="High values of median NUSE are indicative of a problematic array")
#line 226 "R_for_gene_exp.Rnw"


###################################################
### chunk number 13: fig4plot
###################################################
#line 236 "R_for_gene_exp.Rnw"
RNAdeg <- AffyRNAdeg(da)
plotAffyRNAdeg(RNAdeg, cols=cols)
legend("topleft", sampleNames(da), lty=1, col=cols)
box()


###################################################
### chunk number 14: fig4
###################################################
#line 244 "R_for_gene_exp.Rnw"
#line 236 "R_for_gene_exp.Rnw"
RNAdeg <- AffyRNAdeg(da)
plotAffyRNAdeg(RNAdeg, cols=cols)
legend("topleft", sampleNames(da), lty=1, col=cols)
box()
#line 245 "R_for_gene_exp.Rnw"


###################################################
### chunk number 15: loading agilent
###################################################
#line 266 "R_for_gene_exp.Rnw"
library(Agi4x44PreProcess)
targets <- read.targets(infile='mRNA_labelling.txt')


###################################################
### chunk number 16: read Agilent data eval=FALSE
###################################################
## #line 271 "R_for_gene_exp.Rnw"
## ## read the data in
## da <- read.AgilentFE(targets, makePLOT=FALSE)
## 
## ## background correction and normalization
## da.n <- BGandNorm(da, BGmethod = "normexp", NORMmethod = "quantile", foreground = "ProcessedSignal", background = "BGUsed", offset = 50, makePLOTpre = FALSE, makePLOTpost = FALSE)
## ## Data are log-2 transformed see ?BGandNorm


###################################################
### chunk number 17: unspecific filtering eval=FALSE
###################################################
## #line 320 "R_for_gene_exp.Rnw"
## ## Unspecific filtering
## da.f <- filter.probes(da.n,
##                 control=TRUE,
##                 wellaboveBG=TRUE,
##                 isfound=TRUE,
##                 wellaboveNEG=TRUE,
##                 sat=TRUE,
##                 PopnOL=TRUE,
##                 NonUnifOL=T,
##                 nas=TRUE,
##                 limWellAbove=75,
##                 limISF=75,
##                 limNEG=75,
##                 limSAT=75,
##                 limPopnOL=75,
##                 limNonUnifOL=75,
##                 limNAS=100,
##                 makePLOT=F,annotation.package="HsAgilentDesign026652.db",flag.counts=T,targets)
## 
## ## Summarize signal
## da.s <- summarize.probe(da.f, makePLOT=FALSE, targets)
## 
## rownames(targets) <- colnames(da.s)
## 
## ## building a expression set
## da.eset <- build.eset(da.s, targets, makePLOT = FALSE, annotation.package = "HsAgilentDesign026652.db")


###################################################
### chunk number 18: RMA
###################################################
#line 400 "R_for_gene_exp.Rnw"
da.rma <- rma(da)


###################################################
### chunk number 19: exprs
###################################################
#line 404 "R_for_gene_exp.Rnw"
da.eset <- exprs(da.rma)
dim(da.eset)
colnames(da.eset)


###################################################
### chunk number 20: RP
###################################################
#line 432 "R_for_gene_exp.Rnw"
library(RankProd)
cl <- c(rep(0,3), rep(1,4))
options(width=100)
head(da.eset[,c(4:6, 7:10)])
options(width=80)
da.rp <- RP(da.eset[,c(4:6, 7:10)], cl=cl, logged=TRUE, num.perm=100, plot=FALSE, rand=5432)


###################################################
### chunk number 21: extracting_gene
###################################################
#line 443 "R_for_gene_exp.Rnw"
library(mouse4302.db)
gnames <- as.vector(unlist(as.list(mouse4302SYMBOL)))
r.nsc.1fipsc <- topGene(da.rp, cutoff = 0.05, method = "pfp", logged = TRUE, 
	logbase = 2, gene.names=gnames)
# The genes significantly up-regulated
head(r.nsc.1fipsc$Table1, 20)
# How many genes are in table 2
dim(r.nsc.1fipsc$Table1)

# The genes significantly down-regulated
head(r.nsc.1fipsc$Table2, 20)
# how many genes are in table 2
dim(r.nsc.1fipsc$Table2)


###################################################
### chunk number 22: export
###################################################
#line 461 "R_for_gene_exp.Rnw"
# up reg
x <- r.nsc.1fipsc$Table1
# replace the fold change value by their log2 conterpart
x[,3] <- log2(1/x[,3])
write.table(x, file = "NSC_vs_1F_iPSC_up.txt", sep = "\t", quote = FALSE)

# down reg
x <- r.nsc.1fipsc$Table2
# replace the fold change value by their log2 conterpart
x[,3] <- log2(1/x[,3])
write.table(x, file = "NSC_vs_1F_iPSC_down.txt", sep = "\t", quote = FALSE)


###################################################
### chunk number 23: gene list 2
###################################################
#line 500 "R_for_gene_exp.Rnw"
library(RankProd)
library(affy)
cl <- c(rep(0,4), rep(1,3))
da.rp <- RP(da.eset[,c(7:10, 1:3)], cl=cl, logged=TRUE, num.perm=100, plot=FALSE, rand=5432)


###################################################
### chunk number 24: annot_package
###################################################
#line 508 "R_for_gene_exp.Rnw"
## ANNOTATION PACKAGE
library(mouse4302.db)
gnames <- as.vector(unlist(as.list(mouse4302SYMBOL)))


###################################################
### chunk number 25: genes
###################################################
#line 515 "R_for_gene_exp.Rnw"
r.nsc.nsc_1f <- topGene(da.rp, cutoff = 0.05, method = "pfp", logged = TRUE, logbase = 2, gene.names=gnames)
# The genes significantly up-regulated
head(r.nsc.nsc_1f$Table1, 20)
# how many genes are in table 2
dim(r.nsc.nsc_1f$Table1)

# The genes significantly down-regulated
head(r.nsc.nsc_1f$Table2, 20)
# how many genes are in table 2
dim(r.nsc.nsc_1f$Table2)


###################################################
### chunk number 26: fig5plot
###################################################
#line 546 "R_for_gene_exp.Rnw"
library(bioDist)
## Pearson correlation dissimilarity
d <- cor.dist(t(da.eset)) # careful here don't forget to transpose your matrix
dim(as.matrix(d))

#####################################################
## HIERARCHICAL CLUSTERING
#####################################################
## dendrogram
hc = hclust(d, method = "average")
plot(hc, labels = colnames(da.eset), main = "Hier. clust. Pearson", hang=-1)
# effect of the hang = -1
plot(hc, labels = colnames(da.eset), main = "Hier. clust. Pearson")


###################################################
### chunk number 27: fig5
###################################################
#line 563 "R_for_gene_exp.Rnw"
#line 546 "R_for_gene_exp.Rnw"
library(bioDist)
## Pearson correlation dissimilarity
d <- cor.dist(t(da.eset)) # careful here don't forget to transpose your matrix
dim(as.matrix(d))

#####################################################
## HIERARCHICAL CLUSTERING
#####################################################
## dendrogram
hc = hclust(d, method = "average")
plot(hc, labels = colnames(da.eset), main = "Hier. clust. Pearson", hang=-1)
# effect of the hang = -1
plot(hc, labels = colnames(da.eset), main = "Hier. clust. Pearson")
#line 564 "R_for_gene_exp.Rnw"


###################################################
### chunk number 28: fig6_1
###################################################
#line 573 "R_for_gene_exp.Rnw"
library(graphics)
d <- cor.dist(as.matrix(USArrests))

# Let's compare the linkage methods
par(mfrow = c(2,1))
hc <- hclust(d, "average")
plot(hc, main = "AVERAGE", hang = -1)

hc = hclust(d, method = "single")
plot(hc, main = "SINGLE", hang=-1)


###################################################
### chunk number 29: fig6_2
###################################################
#line 585 "R_for_gene_exp.Rnw"
hc = hclust(d, method = "complete")
plot(hc, main = "COMPLETE", hang=-1)

hc = hclust(d, method = "ward")
plot(hc, main = "WARD", hang=-1)

# ressetting the graphic device layout to default
par(mfrow = c(1,1))


###################################################
### chunk number 30: fig6_1
###################################################
#line 597 "R_for_gene_exp.Rnw"
#line 573 "R_for_gene_exp.Rnw"
library(graphics)
d <- cor.dist(as.matrix(USArrests))

# Let's compare the linkage methods
par(mfrow = c(2,1))
hc <- hclust(d, "average")
plot(hc, main = "AVERAGE", hang = -1)

hc = hclust(d, method = "single")
plot(hc, main = "SINGLE", hang=-1)
#line 598 "R_for_gene_exp.Rnw"


###################################################
### chunk number 31: fig6_2
###################################################
#line 607 "R_for_gene_exp.Rnw"
#line 585 "R_for_gene_exp.Rnw"
hc = hclust(d, method = "complete")
plot(hc, main = "COMPLETE", hang=-1)

hc = hclust(d, method = "ward")
plot(hc, main = "WARD", hang=-1)

# ressetting the graphic device layout to default
par(mfrow = c(1,1))
#line 608 "R_for_gene_exp.Rnw"


###################################################
### chunk number 32: fig7plot
###################################################
#line 616 "R_for_gene_exp.Rnw"
library(cluster)
hc.d <- diana(d)
plot(hc.d, which.plots=2, main = "DIVISIVE HIERARCHICAL CLUSTERING")


###################################################
### chunk number 33: fig7
###################################################
#line 624 "R_for_gene_exp.Rnw"
#line 616 "R_for_gene_exp.Rnw"
library(cluster)
hc.d <- diana(d)
plot(hc.d, which.plots=2, main = "DIVISIVE HIERARCHICAL CLUSTERING")
#line 625 "R_for_gene_exp.Rnw"


###################################################
### chunk number 34: fig8plot
###################################################
#line 635 "R_for_gene_exp.Rnw"
library(gplots)
library(RColorBrewer)
library(bioDist)
# First we build a color palette
hmcol = colorRampPalette(brewer.pal(10, "RdBu"))(256)

# Recomputing the dissimilarity matrix for the gene expression values
d <- cor.dist(t(da.eset))
## Heatmap for the samples
heatmap.2(as.matrix(d), 
  distfun=function(x){as.dist(x)}, 
  hclustfun=function(m){hclust(m, method="average")},
  symm=F, col=hmcol, trace='none', notecol='black', 
  denscol='black', notecex=0.8, dendrogram="column")


###################################################
### chunk number 35: fig8
###################################################
#line 653 "R_for_gene_exp.Rnw"
#line 635 "R_for_gene_exp.Rnw"
library(gplots)
library(RColorBrewer)
library(bioDist)
# First we build a color palette
hmcol = colorRampPalette(brewer.pal(10, "RdBu"))(256)

# Recomputing the dissimilarity matrix for the gene expression values
d <- cor.dist(t(da.eset))
## Heatmap for the samples
heatmap.2(as.matrix(d), 
  distfun=function(x){as.dist(x)}, 
  hclustfun=function(m){hclust(m, method="average")},
  symm=F, col=hmcol, trace='none', notecol='black', 
  denscol='black', notecex=0.8, dendrogram="column")
#line 654 "R_for_gene_exp.Rnw"


###################################################
### chunk number 36: fig9plot
###################################################
#line 662 "R_for_gene_exp.Rnw"
r.nsc.nsc_1f <- topGene(da.rp, cutoff = 0.001, method = "pfp", logged = TRUE, logbase = 2, gene.names=gnames)
heatmap.2(da.eset[r.nsc.nsc_1f$Table1[,1],], 
  distfun=function(x){as.dist(1 - cor(t(x), use = "complete.obs", method ="pearson"))}, 
  hclustfun=function(m){hclust(m, method="average")},
  symm=F, col=hmcol, trace='none', notecol='black', 
  denscol='black', notecex=0.8, dendrogram="column")


###################################################
### chunk number 37: fig9
###################################################
#line 672 "R_for_gene_exp.Rnw"
#line 662 "R_for_gene_exp.Rnw"
r.nsc.nsc_1f <- topGene(da.rp, cutoff = 0.001, method = "pfp", logged = TRUE, logbase = 2, gene.names=gnames)
heatmap.2(da.eset[r.nsc.nsc_1f$Table1[,1],], 
  distfun=function(x){as.dist(1 - cor(t(x), use = "complete.obs", method ="pearson"))}, 
  hclustfun=function(m){hclust(m, method="average")},
  symm=F, col=hmcol, trace='none', notecol='black', 
  denscol='black', notecex=0.8, dendrogram="column")
#line 673 "R_for_gene_exp.Rnw"


###################################################
### chunk number 38: over-representation
###################################################
#line 757 "R_for_gene_exp.Rnw"
library(GOstats)
library(mouse4302.db)
library(RankProd)

r.nsc.nsc_1f <- topGene(da.rp, cutoff = 0.001, method = "pfp", 
logged = TRUE, logbase = 2, gene.names = rownames(da.rp$AveFC))

uniqueId <- mouse4302ENTREZID
entrezUniverse <- unique(as.character(uniqueId))

# Convert our list of Affymetrix probe IDs to
# NCBI Entrez gene ID
g.list <- c(rownames(r.nsc.nsc_1f$Table1), rownames(r.nsc.nsc_1f$Table2))
ourlist <- mouse4302ENTREZID[g.list]
ourlist <- unique(as.character(ourlist))

# creating the GOHyperGParams object
params = new("GOHyperGParams", geneIds=ourlist, 
universeGeneIds=entrezUniverse, annotation='mouse4302.db', 
ontology="BP", pvalueCutoff=0.001, conditional=FALSE, testDirection="over")

# running the test
mfhyper = hyperGTest(params)


###################################################
### chunk number 39: hypergeo
###################################################
#line 783 "R_for_gene_exp.Rnw"
mfhyper

head(summary(mfhyper))

# grabbing detail of a GO category
GOTERM[["GO:0032502"]]

# Information on the Directed Acyclic Graph (DAG)
goDag(mfhyper)

# how many gene were mapped in the end?
geneMappedCount(mfhyper)

# how many genes are in the universe
universeMappedCount(mfhyper)

# html output
htmlReport(mfhyper, file="BP_list_significant.html")


###################################################
### chunk number 40: fig13plot
###################################################
#line 805 "R_for_gene_exp.Rnw"
library(Rgraphviz)
g1 <- GOGraph(head(summary(mfhyper))$GOBPID, GOBPPARENTS)
plot(g1)

# display the label in the nodes
my.labels <- vector()
for(i in 1:length(slot(g1, "nodes"))){
  my.labels[i] <- Term(get(slot(g1, "nodes")[[i]], GOTERM))
}
my.labels

nodattr <- makeNodeAttrs(g1, label=my.labels, 
  shape = "ellipse", fillcolor = "#f2f2f2", fixedsize = FALSE)

plot(g1, nodeAttrs = nodattr)


###################################################
### chunk number 41: fig13
###################################################
#line 824 "R_for_gene_exp.Rnw"
#line 805 "R_for_gene_exp.Rnw"
library(Rgraphviz)
g1 <- GOGraph(head(summary(mfhyper))$GOBPID, GOBPPARENTS)
plot(g1)

# display the label in the nodes
my.labels <- vector()
for(i in 1:length(slot(g1, "nodes"))){
  my.labels[i] <- Term(get(slot(g1, "nodes")[[i]], GOTERM))
}
my.labels

nodattr <- makeNodeAttrs(g1, label=my.labels, 
  shape = "ellipse", fillcolor = "#f2f2f2", fixedsize = FALSE)

plot(g1, nodeAttrs = nodattr)
#line 825 "R_for_gene_exp.Rnw"


###################################################
### chunk number 42: PAM
###################################################
#line 835 "R_for_gene_exp.Rnw"
library(cluster)
library(lattice)
library(bioDist)
# K-means
# we want to find the different type of gene expression pattern among 
# the genes significantly regulated 
sub.da <- da.eset[g.list,]
# compute the correlation matrix between the genes
d <- cor.dist(sub.da)
# Cluster the genes for k = 3
g.pam <- pam(d, k = 3)


###################################################
### chunk number 43: fig14plot
###################################################
#line 850 "R_for_gene_exp.Rnw"
x <- vector()
for(i in 3:20){
  g.pam <- pam(d, k = i)
  x[i] <- g.pam$silinfo$avg.width
}

plot(x, xlab = 'K (cluster number requested)', ylab = 'average silhouette width')


###################################################
### chunk number 44: fig14
###################################################
#line 861 "R_for_gene_exp.Rnw"
#line 850 "R_for_gene_exp.Rnw"
x <- vector()
for(i in 3:20){
  g.pam <- pam(d, k = i)
  x[i] <- g.pam$silinfo$avg.width
}

plot(x, xlab = 'K (cluster number requested)', ylab = 'average silhouette width')
#line 862 "R_for_gene_exp.Rnw"


###################################################
### chunk number 45: fig15plot
###################################################
#line 870 "R_for_gene_exp.Rnw"
g.pam <- pam(d, k = 8)

# transform sub.da
x <- data.frame()
for(i in 1:ncol(sub.da)){
  x <- rbind(x, data.frame(gene.expression = sub.da[,i], cell.type = colnames(sub.da)[i]))
}
# plot the data
df <- data.frame(cluster = as.vector(g.pam$clustering), x)
bwplot(gene.expression ~ cell.type | cluster, data = df)


