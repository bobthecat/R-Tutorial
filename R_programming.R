###################################################
### chunk number 1: vectors1
###################################################
#line 183 "R_programming.Rnw"
# characters
x <- c("Lincoln", "Roosevelt", "Jackson")
# access the second elment of the vector
x[2]
# numeric vectors
x <- 1:20
x <- rep("A", 20)


###################################################
### chunk number 2: vector2
###################################################
#line 193 "R_programming.Rnw"
# operation on vectors
# an operation apply to all the elements
x <- c(1, 1, 2, 3, 5, 8)
x
x - 1
x * 2
x^2


###################################################
### chunk number 3: vector3
###################################################
#line 203 "R_programming.Rnw"
# basic functions
mean(x)
median(x)
sd(x)
sum(x)
summary(x)


###################################################
### chunk number 4: matrix
###################################################
#line 214 "R_programming.Rnw"
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
### chunk number 5: data.frame
###################################################
#line 239 "R_programming.Rnw"
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
### chunk number 6: list
###################################################
#line 257 "R_programming.Rnw"
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


###################################################
### chunk number 7: reshaping
###################################################
#line 280 "R_programming.Rnw"
library(reshape)
# pivot file example from Digithead's Lab Notebook
URL <- "http://www.stanford.edu/~druau/pivot_table.csv"
pivot <- read.table(URL, sep=',', header=TRUE)
head(pivot, 20)
# 
cast(pivot, gene ~ condition)


###################################################
### chunk number 8: install from console eval=FALSE
###################################################
## #line 296 "R_programming.Rnw"
## install.package("mypackage")


###################################################
### chunk number 9: install from terminal eval=FALSE
###################################################
## #line 305 "R_programming.Rnw"
## # R CMD INSTALL mypackage.tar.gz


###################################################
### chunk number 10: importing text files
###################################################
#line 318 "R_programming.Rnw"
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
### chunk number 11: Excel files eval=FALSE
###################################################
## #line 336 "R_programming.Rnw"
## library(xlsx)
## # there are many function in this package. A help page can be found through:
## vignette('xlsx')
## # to read Excel files use
## ?read.xlsx


###################################################
### chunk number 12: importing other format
###################################################
#line 347 "R_programming.Rnw"
library(foreign)
# to list all the function in a package
ls('package:foreign')


###################################################
### chunk number 13: foreach eval=FALSE
###################################################
## #line 430 "R_programming.Rnw"
## library(foreach)
## library(doMC)
## library(multicore)
## 
## # Here we detect the number of cores available on your machine
## ncore = multicore:::detectCores()
## # And here we define the 
## registerDoMC(cores = ncore)
## 
## fib <- foreach(i = 1:50, .combine=c) %dopar% {
## 	# The first 2 number are predefined
##   	if(i <= 2){
##     	fib[i] <- 1
##   	}
##   	else{
##     	fib[i] <- fib[i-1] + fib[i-2]
##   	}
## }
## # displaying the sequence
## fib


###################################################
### chunk number 14: apply
###################################################
#line 458 "R_programming.Rnw"
# let's build a matrix with 5 rows and 5 columns.
mat <- matrix(1:25, nrow=5, byrow=T)
mat
# applying apply on the rows
apply(mat, 1, sum)
# on the column
apply(mat, 2, sum)


###################################################
### chunk number 15: tapply
###################################################
#line 469 "R_programming.Rnw"
n <- 17
fac <- factor(rep(1:3, length = n), levels = 1:5)
fac
# summary
table(fac)
# now we apply the function in a class specific manner.
tapply(1:n, fac, sum)


###################################################
### chunk number 16: simple function
###################################################
#line 485 "R_programming.Rnw"
Mr.euclide <- function(x, y){
	dist <- sqrt(sum((x - y)^2))
	return(dist)
}
x <- c(1, 1)
y <- c(2, 2)
Mr.euclide(x, y)


###################################################
### chunk number 17: cars
###################################################
#line 505 "R_programming.Rnw"
library(datasets)
data(cars)
head(cars)


###################################################
### chunk number 18: fig1plot
###################################################
#line 512 "R_programming.Rnw"
attach(cars)
plot(speed, dist)
detach(cars)


###################################################
### chunk number 19: fig1
###################################################
#line 519 "R_programming.Rnw"
#line 512 "R_programming.Rnw"
attach(cars)
plot(speed, dist)
detach(cars)
#line 520 "R_programming.Rnw"


###################################################
### chunk number 20: fig2plot
###################################################
#line 529 "R_programming.Rnw"
normal <- rnorm(1000, 1)
par(mfrow=c(1,3))
hist(normal, main = "HISTOGRAM")
barplot(normal, main = "BARPLOT")
plot(density(normal), main = "DENSITY")


###################################################
### chunk number 21: fig2
###################################################
#line 538 "R_programming.Rnw"
#line 529 "R_programming.Rnw"
normal <- rnorm(1000, 1)
par(mfrow=c(1,3))
hist(normal, main = "HISTOGRAM")
barplot(normal, main = "BARPLOT")
plot(density(normal), main = "DENSITY")
#line 539 "R_programming.Rnw"


###################################################
### chunk number 22: fig3plot
###################################################
#line 548 "R_programming.Rnw"
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
### chunk number 23: fig3
###################################################
#line 561 "R_programming.Rnw"
#line 548 "R_programming.Rnw"
# very handy function here to "attach" a data frame to the
# environment. Allowing to access the variable much more
# easily 
attach(cars)
par(mfrow=c(2,1))
plot(speed, dist)
boxplot(dist ~ speed)
# finally we "detach" the data frame
detach(cars)
#line 562 "R_programming.Rnw"


###################################################
### chunk number 24: fig4plot
###################################################
#line 571 "R_programming.Rnw"
myc.ctl <- rnorm(20, mean = 6, sd = 2)
myc.drug <- rnorm(20, mean = 10, sd = 2)
condition <- c(rep("ctl", 10),rep("drug", 10))
experiment <- data.frame(MYC = c(myc.ctl, myc.drug), condition)
par(mfrow=c(1,2))
boxplot(MYC ~ condition, data = experiment)
boxplot(MYC ~ condition, data = experiment, col=c("tomato", "dodgerblue"))


###################################################
### chunk number 25: fig4
###################################################
#line 582 "R_programming.Rnw"
#line 571 "R_programming.Rnw"
myc.ctl <- rnorm(20, mean = 6, sd = 2)
myc.drug <- rnorm(20, mean = 10, sd = 2)
condition <- c(rep("ctl", 10),rep("drug", 10))
experiment <- data.frame(MYC = c(myc.ctl, myc.drug), condition)
par(mfrow=c(1,2))
boxplot(MYC ~ condition, data = experiment)
boxplot(MYC ~ condition, data = experiment, col=c("tomato", "dodgerblue"))
#line 583 "R_programming.Rnw"


###################################################
### chunk number 26: venn diagram install eval=FALSE
###################################################
## #line 596 "R_programming.Rnw"
## install.packages("Vennerable", repos="http://R-Forge.R-project.org")


###################################################
### chunk number 27: venn_diagram
###################################################
#line 602 "R_programming.Rnw"
library(Vennerable)
listA <- LETTERS[1:10]
listB <- LETTERS[5:20]
vv <- list(drug=listA, no.drug=listB)
vv <- Venn(vv)
plot(vv)


###################################################
### chunk number 28: venn_diagram
###################################################
#line 612 "R_programming.Rnw"
#line 602 "R_programming.Rnw"
library(Vennerable)
listA <- LETTERS[1:10]
listB <- LETTERS[5:20]
vv <- list(drug=listA, no.drug=listB)
vv <- Venn(vv)
plot(vv)
#line 613 "R_programming.Rnw"


###################################################
### chunk number 29: googleVis_earthquakes eval=FALSE
###################################################
## #line 622 "R_programming.Rnw"
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
### chunk number 30: slopegraph
###################################################
#line 648 "R_programming.Rnw"
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
### chunk number 31: fig_slopegraph
###################################################
#line 693 "R_programming.Rnw"
table.graph(WorldPhones[,1:2], label.cex=0.7, title.cex=1.2, mar=c(5, 5, 1, 5))


###################################################
### chunk number 32: welch
###################################################
#line 708 "R_programming.Rnw"
t.test(MYC ~ condition, data = experiment, alternative = "two.sided")


###################################################
### chunk number 33: vartest
###################################################
#line 713 "R_programming.Rnw"
var.test(MYC ~ condition, data = experiment)


###################################################
### chunk number 34: stu.test
###################################################
#line 718 "R_programming.Rnw"
t.test(MYC ~ condition, data = experiment, alternative = "two.sided", var.equal = TRUE)


###################################################
### chunk number 35: wilcox.text
###################################################
#line 724 "R_programming.Rnw"
wilcox.test(MYC ~ condition, data = experiment, alternative = "two.sided")


###################################################
### chunk number 36: intake1
###################################################
#line 730 "R_programming.Rnw"
library(ISwR)
attach(intake)
intake


###################################################
### chunk number 37: intake2
###################################################
#line 737 "R_programming.Rnw"
t.test(pre, post, paired = T)


###################################################
### chunk number 38: intake3
###################################################
#line 741 "R_programming.Rnw"
t.test(pre, post)


###################################################
### chunk number 39: fig5plot
###################################################
#line 748 "R_programming.Rnw"
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
### chunk number 40: fig5
###################################################
#line 765 "R_programming.Rnw"
#line 748 "R_programming.Rnw"
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
#line 766 "R_programming.Rnw"


###################################################
### chunk number 41: faithful
###################################################
#line 774 "R_programming.Rnw"
fit <- lm(waiting ~ eruptions)
summary(fit)


###################################################
### chunk number 42: fig6plot
###################################################
#line 783 "R_programming.Rnw"
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
### chunk number 43: fig6
###################################################
#line 798 "R_programming.Rnw"
#line 783 "R_programming.Rnw"
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
#line 799 "R_programming.Rnw"


###################################################
### chunk number 44: pearson
###################################################
#line 807 "R_programming.Rnw"
cor.test(x[1,], x[2,])
cor.test(x[1,], x[3,])


###################################################
### chunk number 45: euclidean
###################################################
#line 813 "R_programming.Rnw"
dist(x[1:2,], method='euc')
dist(x[c(1,3),], method='euc')


