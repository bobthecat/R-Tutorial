##
## -- Lecture 1 --
##  R for beginner
##
##  Created by David Ruau on 2011-04-22.
##  Copyright (c) 2011 Dept. of Pediatrics and Anesthesia. All rights reserved.
##

##########################################
##########################################
## FIRST STEPS
##########################################
##########################################
# R is a fancy calculator
# numerical operation on integers
# square
3 * 3

x <- 3
x

x * 2
# square
x^2

y <- x * 2
1:5
x <- 1:10
# accessing what's inside x
x[2]
x[10]

# R workspace
# to list the content of the workspace
ls()
# to erase one or more object fro the workspace
rm(x)
# to quit R
q()

# to obtain help and documentation on a function
?ls
?mean

##########################################
##########################################
# DATA STRUCTURES AND DATA MANIPULATION
##########################################
##########################################
## VECTORS
# different types data can be store in vectors
# characters
x <- c("Lincoln", "Roosevelt", "Jackson")
# numeric
x <- 1:20
x <- rep("A", 20)

# operation on vectors
# an operation apply to all the elements
x <- c(1, 1, 2, 3, 5, 8)
x - 1
x * 2
x^2

# basic functions
mean(x)
median(x)
sd(x)
sum(x)
summary(x)

## MATRIX
# the really basic way to build a matrix
# lets build a example data set
x <- c(rep(1, 10), rep(3, 10))
x
# Build a matrix with 2 columns
matrix(x, ncol=2)
# Similarly you can have a matrix with 2 rows
matrix(x, nrow=2)
# fill up by rows
matrix(x, nrow=2, byrow=TRUE)


mat <- matrix(x, ncol=2)

# dimension of the matrix
dim(mat)
# accessing one row or column
mat[2,]
mat[,2]

# operation on a matrix
mat * 2

colMeans(mat)
rowMeans(mat)

## DATA FRAME
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

## LIST
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

##########################################
##########################################
# IMPORT EXPORT
##########################################
##########################################
# tab delimited files
# download the file
download.file("https://sites.google.com/site/appliedbioinfo/home/introduction-to-r/treatment.txt", "./treatment.txt")
# import it
pd <- read.table("treatment.txt", sep='\t', header=TRUE)
# check it
pd
# comma separated files
download.file("https://sites.google.com/site/appliedbioinfo/home/introduction-to-r/treatment.txt", "./treatment.txt")
pd <- read.table("treatment.csv", sep=',', header=TRUE)
pd

##########################################
##########################################
# FUNCTIONS
##########################################
##########################################
# "FOR" LOOP
for(i in 1:10){
  print(i)
}
# "IF" CONDITION
i <- 1
if(i == 1){
  print("i is equal 1")
}
# conbine it in a loop !
for(i in 1:10){
  if(i == 2){
    print("hello")
  }
  else{
    print(i)
  }
}

# WRITE YOUR FIRST OWN FOR LOOP
# Goal: generate the sequence of the first 50 Fibonacci numbers
# 1, 1, 2, 3, 5, 8, 13, 21...
# Fi = F[i-1] + F[i-2]

#
# your code here
#

##########################################
##########################################
# PLOT
##########################################
##########################################
## FIRST OVERVIEW OF THE PLOT FUNCTIONALITY
# SCATTER PLOTS
# first lets use a classic data set: cars
data(cars)

# have a look at your data
head(cars)
# how big is it
dim(cars)
# attach te data set to the environment
attach(cars)
# do a scatter plot of the 
plot(speed, dist)
detach(cars)

# HISTOGRAMS AND BARPLOTS
# let's generate a normal distribution with mean 1 and standard devistion 1
# This is the assumption for microarray data expression fold change
normal <- rnorm(1000, 1)
# histogram
hist(normal, main = "HISTOGRAM")
# barplots
barplot(normal, main = "BARPLOT")
# density
plot(density(normal), main = "DENSITY")

# combine graphics
par(mfrow=c(1,3))
hist(normal, main = "HISTOGRAM")
barplot(normal, main = "BARPLOT")
plot(density(normal), main = "DENSITY")
# reset the graphic device
par(mfrow=c(1,1))

# BOXPLOTS
# using the cars data set on the distance needed to stop at a certain speed.
data(cars)
attach(cars)
plot(speed, dist)
# boxplot accept formulas
boxplot(dist ~ speed)
detach(cars)

# another example for the gene expression of MYC in an hypothetical experiment
myc.ctl <- rnorm(10, mean = 6, sd = 2)
myc.drug <- rnorm(10, mean = 10, sd = 2)
condition <- c(rep("ctl", 10),rep("drug", 10))
experiment <- data.frame(MYC = c(myc.ctl, myc.drug), condition)
par(mfrow=c(1,2))
boxplot(MYC ~ condition, data = experiment)
boxplot(MYC ~ condition, data = experiment, col=c("tomato", "dodgerblue"))

##########################################
##########################################
# IMPORT - EXPORT DATA
##########################################
##########################################
# tab delimited file
GSE12499 <- read.table('treatment.txt', sep='\t', header=TRUE)
GSE12499
# data are imported into data frame
class(GSE12499)

# csv file
GSE12499 <- read.table('treatment.csv', sep=',', header=TRUE)
GSE12499

##########################################
##########################################
# BASIC STATISTICS
##########################################
##########################################
## calculate t-values and p-values
## -1- TWO SAMPLES T-TEST NORMALLY DISTRIBUTED
# Welch t-test
t.test(MYC ~ condition, data = experiment, alternative = "two.sided")
# Student's t-test for equal variance
t.test(MYC ~ condition, data = experiment, alternative = "two.sided", var.equal = TRUE)

# Equal variance test
var.test(MYC ~ condition, data = experiment)

## -2- TWO SAMPLES T-TEST NOT NORMALLY DISTRIBUTED
# Wilcoxon test
wilcox.test(MYC ~ condition, data = experiment, alternative = "two.sided")

## -3- PAIRED SAMPLE MEASUREMENT
library(ISwR)
attach(intake)
intake
# 11 women are measured twice
# We spcify that the sample are paired
t.test(pre, post, paired = T)
# if we consider the sample not paired
t.test(pre, post)

## LINEAR REGRESSION

attach(faithful)
plot(faithful, main="Eruption of Old Faithful",
  xlab = "Eruption time (min)",
  ylab = "Waiting time to next eruption (min)")

# fit a linear regression
fit <- lm(waiting ~ eruptions)
# plot it
abline(fit)
# Obtain a summary
summary(fit)

# compute the confidence intervals
pc <- predict(fit, int="c")

# plot the confidence intervals
matlines(eruptions, pc, lty = c(0, 2, 2), col = "blue")

## CORRELATION
# let's consider this fake gene expression matrix
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

# Pearson correlation
cor.test(x[1,], x[2,])
cor.test(x[1,], x[3,])

# Euclidean distance
dist(x[1:2,], method='euc')
dist(x[c(1,3),], method='euc')

############# END #############




