# Collection of packages for data wrangling and visualization
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# NorTA (Normal To Anything): simulates correlated count data via Gaussian copula
# https://github.com/nome-repo/ToyModel
library(ToyModel)

# Quantile functions for zero-inflated and extended distributions (qzinegbin)
# https://cran.r-project.org/web/packages/VGAM/index.html
library(VGAM)


# -------- SIMULATION PARAMETERS --------

# Sample size
n <- 10^3

# Correlation matrix: D variables with a strong negative correlation between var 2 and var 4
D <- 5
R <- diag(D)
R[4,2] <- R[2,4] <- -.9


# -------- NorTA PROCEDURE --------

# generate a multivariate Normal distribution with the desired correlation structure
rnorm <- mvtnorm::rmvnorm(n=n, sigma=R)

# choose margins to Uniform [0,1] via the Normal CDF 
unif <- stats::pnorm(rnorm)

# apply inverse CDF of a zero-inflated negative binomial to each margin (NorTA step)
NorTA <- matrix(0, nrow=n, ncol=D)
NorTA[,1] <- qzinegbin(p=unif[,1], munb=20, size=30, pstr0=.25)
NorTA[,2] <- qzinegbin(p=unif[,2], munb=20, size=30, pstr0=.25)
NorTA[,3] <- qzinegbin(p=unif[,3], munb=20, size=30, pstr0=.25)
NorTA[,4] <- qzinegbin(p=unif[,4], munb=20, size=30, pstr0=.25)
NorTA[,5] <- qzinegbin(p=unif[,5], munb=20, size=30, pstr0=.25)

# verify the empirical correlation between var 2 and var 4 (expected ≈ −0.9)
cor(NorTA[,2], NorTA[,4])