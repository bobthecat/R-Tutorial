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
3 * 3
x <- 3
x
x * 2
y <- x * 2

# R workspace
# to list the content of the workspace
ls()
# to erase one or more object fro the workspace
rm(x)

# to obtain help and documentation on a function
# Functions are organized into packages like files into folders
?ls

##########################################
##########################################
# DATA STRUCTURES
##########################################
##########################################
# INTEGERS

# VECTORS

# MATRIX

# DATA FRAME

# LIST

##########################################
##########################################
# DATA MANIPULATION
##########################################
##########################################
# OPERATION ON VECTORS

# OPERATION ON MATRICES

# OPERATION ON DATA FRAME

# OPERATION ON LIST


##########################################
##########################################
# FUNCTIONS
##########################################
##########################################

?t.test # belong to the "stats" package

# to search the help
??heatmap


## basic plot functionality
# scatter plots



# histogram
x <- rnorm(1000, 1)
y <- rnorm(1000, 2)

hist(x)
hist(y)

plot(density(x))
plot(density(y))
fivenum(x)
fivenum(y)

plot(density(x), xlim = c(-2,6), ylim = c(0, 0.4))
par(new = T)
plot(density(y), xlim = c(-2,6), ylim = c(0, 0.4), col = 'red')
legend("topleft", legend = c("x", "y"), fill = c("black", "red"))

t.test(x, y)


