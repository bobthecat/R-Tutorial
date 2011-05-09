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

## COMPARE NSC TO 1F_NSC
cl <- c(rep(0,4), rep(1,3))
da.rp <- RP(da.eset[,c(7:10, 1:3)], cl=cl, logged=TRUE, num.perm=100, plot=FALSE, rand=5432)

## ANNOTATION PACKAGE
library(mouse4302.db)
gnames <- as.vector(unlist(as.list(mouse4302SYMBOL)))

#####################################################
## EXTRACTING THE SIGNIFICANTLY REGULATED GENES
#####################################################
r.nsc.nsc_1f <- topGene(da.rp, cutoff = 0.05, method = "pfp", logged = TRUE, logbase = 2, gene.names=gnames)
# The genes significantly up-regulated
head(r.nsc.nsc_1f$Table1, 20)
# how many genes are in table 2
dim(r.nsc.nsc_1f$Table1)

# The genes significantly down-regulated
head(r.nsc.nsc_1f$Table2, 20)
# how many genes are in table 2
dim(r.nsc.nsc_1f$Table2)

#####################################################
## CORRELATION MATRIX
#####################################################
## required library
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

# Alternative hierarchical clustering method exist
library(graphics)
# This data set contains statistics, in arrests per 100,000 residents for assault, murder, and rape in each of the 50 US states in 1973. Also given is the percent of the population living in urban areas.
d <- cor.dist(as.matrix(USArrests))

# Let's compare the linkage methods
par(mfrow = c(4,1))
hc <- hclust(d, "average")
plot(hc, main = "AVERAGE", hang = -1)

hc = hclust(d, method = "single")
plot(hc, main = "SINGLE", hang=-1)

hc = hclust(d, method = "complete")
plot(hc, main = "COMPLETE", hang=-1)

hc = hclust(d, method = "ward")
plot(hc, main = "WARD", hang=-1)

# ressetting the graphic device layout to default
par(mfrow = c(1,1))

# DIVISIVE HIERARCHICAL CLUSTERING
library(cluster)
hc.d <- diana(d)
plot(hc.d, which.plots=2, main = "DIVISIVE HIERARCHICAL CLUSTERING")

# to save your plot into a PDF
pdf(file='dendogram.pdf', width=8, height=8)
plot(hc, labels = colnames(da.eset), main = "Hier. clust. Pearson", hang=-1)
dev.off()

#####################################################
## HEAMAP AND HIERARCHICAL CLUSTERING
#####################################################
# LIBRARY REQUIRED
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

## Heatmap for the genes
# pdf(file='heatmap_all_array.pdf', height=8, width=8)
r.nsc.nsc_1f <- topGene(da.rp, cutoff = 0.001, method = "pfp", logged = TRUE, logbase = 2, gene.names=gnames)
heatmap.2(da.eset[r.nsc.nsc_1f$Table1[,1],], 
  distfun=function(x){as.dist(1 - cor(t(x), use = "complete.obs", method ="pearson"))}, 
  hclustfun=function(m){hclust(m, method="average")},
  symm=F, col=hmcol, trace='none', notecol='black', 
  denscol='black', notecex=0.8, dendrogram="column")
# dev.off()

#####################################################
## VOLCANO PLOTS
#####################################################
library(annotate)
library(mouse4302.db)
library(simpleaffy)
# VOLCANO PLOT: "fold change" x "p-values"
# with rank product p-values
pval <- apply(da.rp$pval, 1, min)
volcano <- data.frame(FC = -{da.rp$AveFC}, pval = -log(pval))
plot(-{da.rp$AveFC}, -log(pval))

# using simpleaffy to extract the fold change p-values
results <- pairwise.comparison(da.n, "cell_type", c("NSC","NSC_1F"), spots=da)
plot(-{da.rp$AveFC}, -log(results@tt))

# displaying the gene significantly regulated in the plot
# the gene list
g.list <- c(rownames(r.nsc.nsc_1f$Table1), rownames(r.nsc.nsc_1f$Table2))

# determining the coordinates of the genes
x <- -{da.rp$AveFC[g.list,]}
y <- -{log(results@tt[g.list])}
z <- getSYMBOL(g.list, "mouse4302")
points(x, y, col = "green", pch = 19)
# to display the gene names
# text(x, y, label = z, pos = 4, offset = 0.5)


#####################################################
## SCATTER PLOTS
#####################################################
plot(results@means, 
  xlab = "NSC", ylab = "NSC_1F", 
  main = "Gene differentially expressed between NSC and NSC_1F")

# Highliting the genes found significantly regulated
x <- results@means[g.list, "NSC"]
y <- results@means[g.list, "NSC_1F"]
z <- getSYMBOL(g.list, "mouse4302")
points(x, y, col = "green", pch = 19)
# to display the gene names
# text(x, y, label = z, pos = 4, offset = 0.5)


#####################################################
## GENE ONTOLOGY OVER-REPRESENTATION ANALYSIS
#####################################################
library(GOstats)
library(mouse4302.db)

uniqueId <- mouse4302ENTREZID
entrezUniverse <- unique(as.character(uniqueId))

# Convert our list of Affymetrix probe IDs to
# NCBI Entrez gene ID
ourlist <- mouse4302ENTREZID[g.list]
ourlist <- unique(as.character(ourlist))

# creating the GOHyperGParams object
params = new("GOHyperGParams", geneIds=ourlist, 
universeGeneIds=entrezUniverse, annotation='mouse4302.db', 
ontology="BP", pvalueCutoff=0.001, conditional=FALSE, testDirection="over")

# running the test
mfhyper = hyperGTest(params)

# Looking at the output of the test
hist(pvalues(mfhyper), breaks=50, col="mistyrose")

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

# graph
library(Rgraphviz)
g1 <- GOGraph(head(summary(mfhyper))$GOBPID, GOBPPARENTS)
plot(g1)

# display the label in the nodes
my.labels <- vector()
for(i in 1:length(g1@nodes)){
  my.labels[i] <- Term(get(g1@nodes[[i]], GOTERM))
}
my.labels

nodattr <- makeNodeAttrs(g1, label=my.labels, 
  shape = "ellipse", fillcolor = "#f2f2f2", fixedsize = FALSE)

plot(g1, nodeAttrs = nodattr)

#####################################################
## FINDING SIMILAR PATTERN OF GENE EXPRESSION
#####################################################
library(cluster)
library(lattice)
library(bioDist)
# K-means
# we want to find the different type of gene expression pattern among 
# the genes significantly regulated 
sub.da <- da.eset[g.list,]
# compute the correlation matrix between the genes
d <- cor.dist(sub.da)
# Cluster the genes for k = 6
g.pam <- pam(d, k = 6)

# SILHOUETTE
x <- vector()
for(i in 3:20){
  g.pam <- pam(d, k = i)
  x[i] <- g.pam$silinfo$avg.width
}

plot(x, xlab = 'K (cluster number requested)', ylab = 'average silhouette width')

## redo the clsutering for the optimum number of cluster
g.pam <- pam(d, k = 8)

# transform sub.da
x <- data.frame()
for(i in 1:ncol(sub.da)){
  x <- rbind(x, data.frame(gene.expression = sub.da[,i], cell.type = colnames(sub.da)[i]))
}
# plot the data
df <- data.frame(cluster = as.vector(g.pam$clustering), x)
bwplot(gene.expression ~ cell.type | cluster, data = df)




