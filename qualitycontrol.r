# R script for automatic quality control
# argument: the affybatch

qa <- function(abatch) {
	
	require(affy)
	require(simpleaffy)
	require(RColorBrewer)
	require(affyPLM)
	
	if (class(abatch)!= 'AffyBatch') {
		stop("argument must be AffyBatch!")
		}
	
	# colors
	cols <- brewer.pal(12, "Set3")
	
    # Boxplot
    pdf(file='boxplot.pdf', height=8, width=10)
    boxplot(abatch, col=cols, main="Unprocessed log scale probe-level data", xlab="If discrepancy, they are not conclusive\n Difference can be reduce by normalization")
    dev.off()
    
    # Histogram
    pdf(file='histogram.pdf', height=8, width=8)
    hist(abatch, col=cols, xlab="Log(base2) intensities; Bimodal distribution indicate spatial artifact\n Second mode is the result of array(s) having abnormally high value")
    legend("topright", sampleNames(abatch), lty=1,col=cols)
    dev.off()
    
	#RNA degradation
	pdf(file="RNAdeg.pdf", height=8, width=8)
	RNAdeg <- AffyRNAdeg(abatch)
	plotAffyRNAdeg(RNAdeg, cols=cols)
    legend("topleft", sampleNames(abatch), lty=1,col=cols)
    box()
	dev.off()
	
	# simpleaffy graph
	abatch.qc <- qc(abatch)
	pdf(file="QC-simpleaffy.pdf", height=8, width=10)
	plot(abatch.qc)
	dev.off()
	
	source("/Users/druau/Desktop/Arbeit/scripts/my_img_Test.r")
	pset <- fitPLM(abatch)
	
	# false color image control
	for (n in 1:length(abatch)) {
		filename <- paste("QC",as.vector(sampleNames(abatch))[n],".png")
		png(file=filename, height=900, width=800)
		img.Test(abatch,pset,n)
		dev.off()
	}
	
	# RLE plot
	pdf(file="RLE.pdf", height=8, width= 8)
	Mbox(pset, col = cols, main ="RLE (Relative Log Expression)", xlab="Assuming that the majority of the gene are not changing\n Ideally these boxes would have small spread and be centered at M=0")
	dev.off()
	
	# NUSE plot
	pdf(file="NUSE.pdf", height=8, width= 8)
	boxplot(pset, col=cols, main= "NUSE (Normalized Unscaled Standard Error)", xlab="High values of median NUSE are indicative of a problematic array")
	dev.off()    
}
