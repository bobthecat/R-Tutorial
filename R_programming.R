###################################################
### chunk number 1: 
###################################################
#line 121 "R_programming.Rnw"
options(width=60, continue=" ")


###################################################
### chunk number 2: vectors1
###################################################
#line 170 "R_programming.Rnw"
# characters
x <- c("Lincoln", "Roosevelt", "Jackson")
# access the second elment of the vector
x[2]
# numeric vectors
x <- 1:20
x <- rep("A", 20)


###################################################
### chunk number 3: vector2
###################################################
#line 180 "R_programming.Rnw"
# operation on vectors
# an operation apply to all the elements
x <- c(1, 1, 2, 3, 5, 8)
x
x - 1
x * 2
x^2


###################################################
### chunk number 4: vector3
###################################################
#line 190 "R_programming.Rnw"
# basic functions
mean(x)
median(x)
sd(x)
sum(x)
summary(x)


###################################################
### chunk number 5: matrix
###################################################
#line 201 "R_programming.Rnw"
# Build a matrix with 2 columns
matrix(1:20, ncol=2)
# Similarly you can have a matrix with 2 rows
matrix(1:20, nrow=2)
# fill up by rows
matrix(1:20, nrow=2, byrow=TRUE)
# build a matrix object
mat <- matrix(x, ncol=2)

# dimension of the matrix
dim(mat)
# accessing a row or a column
mat[2,]
mat[,2]

# operation on a matrix are vectorized
mat * 2
# useful functions
colMeans(mat)
rowMeans(mat)


###################################################
### chunk number 6: data.frame
###################################################
#line 226 "R_programming.Rnw"
df <- data.frame(
  fruits = c("oranges", "bananas", "apples", "strawberries"), 
  color = c("orange", "yellow", "red", "red")
)

# information on the data frame
attributes(df)
names(df)
# accessing data
df$fruits
df$color
# subsetting a dataframe
subset(df, color=='red')


###################################################
### chunk number 7: list
###################################################
#line 244 "R_programming.Rnw"
my.list <- list(
  fruits = c("oranges", "bananas", "apples", "strawberries"), 
  color = c("orange", "yellow", "red", "red")
)
my.list

my.list$fruits

my.list[1]

# list can store any other object even list
my.list <- list(
  fruits = df,
  mat = mat
)
my.list

# access the second row of the matrix we stored in the list.
my.list[[2]][2,]



###################################################
### chunk number 8: reshaping
###################################################
#line 271 "R_programming.Rnw"
library(reshape)
# pivot file example from Digithead's Lab Notebook
URL <- "http://www.stanford.edu/~druau/pivot_table.csv"
pivot <- read.table(URL, sep=',', header=TRUE)
head(pivot, 20)
# 
cast(pivot, gene ~ condition)


###################################################
### chunk number 9: install from console eval=FALSE
###################################################
## #line 287 "R_programming.Rnw"
## install.package("mypackage")


###################################################
### chunk number 10: install from terminal eval=FALSE
###################################################
## #line 294 "R_programming.Rnw"
## # R CMD INSTALL mypackage.tar.gz


###################################################
### chunk number 11: importing text files
###################################################
#line 305 "R_programming.Rnw"
URL <- "http://www.stanford.edu/~druau/treatment.txt"
# download.file() download a file from the internet.
download.file(URL, "./treatment.txt")
pd <- read.table("treatment.txt", sep='\t', header=TRUE)
pd
# comma separated files
URL <- "http://www.stanford.edu/~druau/treatment.csv"
download.file(URL, "./treatment.csv")
pd <- read.table("treatment.csv", sep=',', header=TRUE)
pd


###################################################
### chunk number 12: Excel files eval=FALSE
###################################################
## #line 323 "R_programming.Rnw"
## library(xlsx)
## # there are many function in this package. A help page can be found through:
## vignette('xlsx')
## # to read Excel files use
## ?read.xlsx


###################################################
### chunk number 13: importing other format eval=FALSE
###################################################
## #line 334 "R_programming.Rnw"
## library(foreign)
## # to list all the function in a package
## ls('package:foreign')
## # to import STATA file
## read.dta()
## # Import SAS files
## read.xport()
## # Import SPSS
## read.spss()


###################################################
### chunk number 14: FOR example
###################################################
#line 376 "R_programming.Rnw"
for(i in 1:10){
	# do something
	print(i)
}


###################################################
### chunk number 15: IF example
###################################################
#line 387 "R_programming.Rnw"
i <- 1
if(i == 1){
  print("i is equal 1")
}


###################################################
### chunk number 16: foreach example
###################################################
#line 425 "R_programming.Rnw"
library(foreach)
library(doMC)
library(multicore)

# Here we detect the number of cores available on your machine
ncore = multicore:::detectCores()
# And here we define the the number of cores we want to use.
registerDoMC(cores = ncore)

results <- foreach(i = 1:10, .combine=c) %dopar% {
	i+i
}
# displaying the results
results


###################################################
### chunk number 17: apply
###################################################
#line 446 "R_programming.Rnw"
# Let's build a matrix with 5 rows and 5 columns.
mat <- matrix(1:25, nrow=5, byrow=T)
mat
# applying apply on the rows
apply(mat, 1, sum)
# on the column
apply(mat, 2, sum)


###################################################
### chunk number 18: tapply
###################################################
#line 457 "R_programming.Rnw"
n <- 17
fac <- factor(rep(1:3, length = n), levels = 1:5)
fac
# summary
table(fac)
# now we apply the function in a class specific manner.
tapply(1:n, fac, sum)


###################################################
### chunk number 19: simple function
###################################################
#line 473 "R_programming.Rnw"
Mr.euclide <- function(x, y){
	dist <- sqrt(sum((x - y)^2))
	return(dist)
}
x <- c(1, 1)
y <- c(2, 2)
Mr.euclide(x, y)


###################################################
### chunk number 20: cars
###################################################
#line 498 "R_programming.Rnw"
library(datasets)
data(cars)
head(cars)


###################################################
### chunk number 21: fig1plot
###################################################
#line 505 "R_programming.Rnw"
attach(cars)
plot(speed, dist)
detach(cars)


###################################################
### chunk number 22: fig1
###################################################
#line 512 "R_programming.Rnw"
#line 505 "R_programming.Rnw"
attach(cars)
plot(speed, dist)
detach(cars)
#line 513 "R_programming.Rnw"


###################################################
### chunk number 23: fig2plot
###################################################
#line 522 "R_programming.Rnw"
normal <- rnorm(1000, 1)
par(mfrow=c(1,3))
hist(normal, main = "HISTOGRAM")
barplot(normal, main = "BARPLOT")
plot(density(normal), main = "DENSITY")


###################################################
### chunk number 24: fig2
###################################################
#line 531 "R_programming.Rnw"
#line 522 "R_programming.Rnw"
normal <- rnorm(1000, 1)
par(mfrow=c(1,3))
hist(normal, main = "HISTOGRAM")
barplot(normal, main = "BARPLOT")
plot(density(normal), main = "DENSITY")
#line 532 "R_programming.Rnw"


###################################################
### chunk number 25: fig3plot
###################################################
#line 541 "R_programming.Rnw"
# very handy function here to "attach" a data frame to the
# environment. Allowing to access the variable much more
# easily 
attach(cars)
par(mfrow=c(2,1))
plot(speed, dist)
boxplot(dist ~ speed)
# finally we "detach" the data frame
detach(cars)


###################################################
### chunk number 26: fig3
###################################################
#line 554 "R_programming.Rnw"
#line 541 "R_programming.Rnw"
# very handy function here to "attach" a data frame to the
# environment. Allowing to access the variable much more
# easily 
attach(cars)
par(mfrow=c(2,1))
plot(speed, dist)
boxplot(dist ~ speed)
# finally we "detach" the data frame
detach(cars)
#line 555 "R_programming.Rnw"


###################################################
### chunk number 27: fig4plot
###################################################
#line 564 "R_programming.Rnw"
myc.ctl <- rnorm(20, mean = 6, sd = 2)
myc.drug <- rnorm(20, mean = 10, sd = 2)
condition <- c(rep("ctl", 10),rep("drug", 10))
experiment <- data.frame(MYC = c(myc.ctl, myc.drug), condition)
par(mfrow=c(1,2))
boxplot(MYC ~ condition, data = experiment)
boxplot(MYC ~ condition, data = experiment, col=c("tomato", "dodgerblue"))


###################################################
### chunk number 28: fig4
###################################################
#line 575 "R_programming.Rnw"
#line 564 "R_programming.Rnw"
myc.ctl <- rnorm(20, mean = 6, sd = 2)
myc.drug <- rnorm(20, mean = 10, sd = 2)
condition <- c(rep("ctl", 10),rep("drug", 10))
experiment <- data.frame(MYC = c(myc.ctl, myc.drug), condition)
par(mfrow=c(1,2))
boxplot(MYC ~ condition, data = experiment)
boxplot(MYC ~ condition, data = experiment, col=c("tomato", "dodgerblue"))
#line 576 "R_programming.Rnw"


###################################################
### chunk number 29: venn diagram install eval=FALSE
###################################################
## #line 589 "R_programming.Rnw"
## install.packages("Vennerable", repos="http://R-Forge.R-project.org")


###################################################
### chunk number 30: venn_diagram
###################################################
#line 595 "R_programming.Rnw"
library(Vennerable)
listA <- LETTERS[1:10]
listB <- LETTERS[5:20]
vv <- list(drug=listA, no.drug=listB)
vv <- Venn(vv)
plot(vv)


###################################################
### chunk number 31: venn_diagram
###################################################
#line 605 "R_programming.Rnw"
#line 595 "R_programming.Rnw"
library(Vennerable)
listA <- LETTERS[1:10]
listB <- LETTERS[5:20]
vv <- list(drug=listA, no.drug=listB)
vv <- Venn(vv)
plot(vv)
#line 606 "R_programming.Rnw"


###################################################
### chunk number 32: googleVis_earthquakes eval=FALSE
###################################################
## #line 615 "R_programming.Rnw"
## library(googleVis)
## library(XML)
## # Get earthquake data of the last 30 days
## eq=readHTMLTable("http://www.iris.edu/seismon/last30.html")
## eq <- eq[[2]] ## extract the earthquake table
## # Filter for Japan earthquakes
## jp <- grep("*japan*", eq$REGION, ignore.case=T)
## eq <- eq[jp,]
## eq$loc=paste(eq$LAT, eq$LON, sep=":") ## create a lat:long location variable
## M <- gvisGeoMap(eq, locationvar="loc", numvar="MAG", hovervar="MAG",  options=list(region="JP"))
## plot(M)


###################################################
### chunk number 33: slopegraph
###################################################
#line 641 "R_programming.Rnw"
table.graph <- function(df, line.col=c("grey", "black"), label.cex=1, title.cex=1, width = 6, digits = 2, rounding.method = NULL, ...) {
  xmin <- min(df)
  xmax <- max(df)
  X1 <- as.numeric(as.vector(df[,1]))
  X2 <- as.numeric(as.vector(df[,2]))
  # original settings
  old.par <- par(no.readonly = TRUE)
  # par settings usually margins
  par(...)
  # rounding
  fmt <- .rd.method(rounding.method, width, digits)
  # left
  plot(rep(0, nrow(df)), X1, xlim=c(0,1), ylim=c(xmin, xmax), 
    axes=FALSE, xlab='', ylab='', type='n')
  mtext(text=paste(rownames(df), sprintf(fmt, X1), sep='  '), side=2, at=X1, las=1, cex=label.cex)
  par(new=TRUE)
  # right
  plot(rep(1, nrow(df)), X2, xlim=c(0,1), ylim=c(xmin, xmax), 
    axes=FALSE, xlab='', ylab='', type='n')
  mtext(text=paste(sprintf(fmt, X2), rownames(df), sep='  '), side=4, at=X2, las=1, cex=label.cex)
  # class label
  mtext(colnames(df)[1], side=3, at=0, cex=title.cex)
  mtext(colnames(df)[2], side=3, at=1, cex=title.cex)
  # lines
  segments(x0 = rep(0, nrow(df)), y0 = X1, x1 = rep(1, nrow(df)), y1 = X2,
   col=ifelse({X1 - X2} < 0, line.col[1], line.col[2]))
  # restore original settings
  par(old.par)
}

.rd.method <- function(rounding.method, width, digits){
  if(is.null(rounding.method)){
    fmt = "%s"
  } 
  else{
    rounding.character <- switch(match(rounding.method, c("round", "signif")), "f", "g")
    fmt = paste("%", width, ".", digits, rounding.character, sep = "")
  }
  return(fmt)
}


###################################################
### chunk number 34: fig_slopegraph
###################################################
#line 686 "R_programming.Rnw"
table.graph(WorldPhones[,1:2], label.cex=0.7, title.cex=1.2, mar=c(5, 5, 1, 5))


###################################################
### chunk number 35: welch
###################################################
#line 701 "R_programming.Rnw"
t.test(MYC ~ condition, data = experiment, alternative = "two.sided")


###################################################
### chunk number 36: vartest
###################################################
#line 706 "R_programming.Rnw"
var.test(MYC ~ condition, data = experiment)


###################################################
### chunk number 37: stu.test
###################################################
#line 711 "R_programming.Rnw"
t.test(MYC ~ condition, data = experiment, alternative = "two.sided", var.equal = TRUE)


###################################################
### chunk number 38: wilcox.text
###################################################
#line 717 "R_programming.Rnw"
wilcox.test(MYC ~ condition, data = experiment, alternative = "two.sided")


###################################################
### chunk number 39: intake1
###################################################
#line 723 "R_programming.Rnw"
library(ISwR)
attach(intake)
intake


###################################################
### chunk number 40: intake2
###################################################
#line 730 "R_programming.Rnw"
t.test(pre, post, paired = T)


###################################################
### chunk number 41: intake3
###################################################
#line 734 "R_programming.Rnw"
t.test(pre, post)


###################################################
### chunk number 42: fig5plot
###################################################
#line 741 "R_programming.Rnw"
attach(faithful)
plot(faithful, main="Eruption of Old Faithful",
  xlab = "Eruption time (min)",
  ylab = "Waiting time to next eruption (min)")

# fit a linear regression
fit <- lm(waiting ~ eruptions)
# plot it
abline(fit)
# compute the confidence intervals
pc <- predict(fit, int="c")
# plot the confidence intervals
matlines(eruptions, pc, lty = c(0, 2, 2), col = "blue")


###################################################
### chunk number 43: fig5
###################################################
#line 758 "R_programming.Rnw"
#line 741 "R_programming.Rnw"
attach(faithful)
plot(faithful, main="Eruption of Old Faithful",
  xlab = "Eruption time (min)",
  ylab = "Waiting time to next eruption (min)")

# fit a linear regression
fit <- lm(waiting ~ eruptions)
# plot it
abline(fit)
# compute the confidence intervals
pc <- predict(fit, int="c")
# plot the confidence intervals
matlines(eruptions, pc, lty = c(0, 2, 2), col = "blue")
#line 759 "R_programming.Rnw"


###################################################
### chunk number 44: faithful
###################################################
#line 767 "R_programming.Rnw"
fit <- lm(waiting ~ eruptions)
summary(fit)


###################################################
### chunk number 45: fig6plot
###################################################
#line 776 "R_programming.Rnw"
x <- matrix(data=c(c(2,3,7,5,4), c(1,2,6,4,2), c(4,6,14,10,8)) ,nrow=3, ncol=5, byrow=T, dimnames=list(c('A','B','C'), c('A','B','C','D','E')))
par(cex=1.5)
plot(x[1,], type='l', ylim=c(0,15), axes=F, xlab=NA, ylab=NA)
points(x[1,], pch=21, bg='black')
par(new=T)
plot(x[2,], type='l', ylim=c(0,15), axes=F, xlab=NA, ylab=NA)
points(x[2,], pch=22, bg='black')
par(new=T)
plot(x[3,], type='l', ylim=c(0,15), xlab="Experiments", ylab="Gene expression")
points(x[3,], pch=23, bg='black')
legend('topright', legend=c("A","B","C"), pch=c(21,22,23), pt.bg='black')


###################################################
### chunk number 46: fig6
###################################################
#line 791 "R_programming.Rnw"
#line 776 "R_programming.Rnw"
x <- matrix(data=c(c(2,3,7,5,4), c(1,2,6,4,2), c(4,6,14,10,8)) ,nrow=3, ncol=5, byrow=T, dimnames=list(c('A','B','C'), c('A','B','C','D','E')))
par(cex=1.5)
plot(x[1,], type='l', ylim=c(0,15), axes=F, xlab=NA, ylab=NA)
points(x[1,], pch=21, bg='black')
par(new=T)
plot(x[2,], type='l', ylim=c(0,15), axes=F, xlab=NA, ylab=NA)
points(x[2,], pch=22, bg='black')
par(new=T)
plot(x[3,], type='l', ylim=c(0,15), xlab="Experiments", ylab="Gene expression")
points(x[3,], pch=23, bg='black')
legend('topright', legend=c("A","B","C"), pch=c(21,22,23), pt.bg='black')
#line 792 "R_programming.Rnw"


###################################################
### chunk number 47: pearson
###################################################
#line 800 "R_programming.Rnw"
cor.test(x[1,], x[2,])
cor.test(x[1,], x[3,])


###################################################
### chunk number 48: euclidean
###################################################
#line 806 "R_programming.Rnw"
dist(x[1:2,], method='euc')
dist(x[c(1,3),], method='euc')


