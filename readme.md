# Introduction for R package kurtICM
We proposed an efficient estimation and statsitical inference on kurtosis of independent component model.
# Getting Started
These instructions will give you a toy example for implementing the package.
## Prerequisites
```
install.packages("devtools")
```
## Install KurtosisICM

```
library("devtools")
devtools::install_github("Bowen-stat/KurtosisICM")
```
## Toy example 
```
rm(list = ls())
set.seed(123)
library(KurtosisICM)
n=200
p=100
X = matrix(rnorm(n*p),nrow = n)
#estimation
ka = kurtICM(X)
# test if the kurtosis is equal to k0
k0=0
obj = Ktest(X,k0)
obj
```