#### Models proposed in the  paper titled: ####
# A multivariate statistical approach to predict 
# COVID-19 count data with epidemiological interpretation and uncertainty quantification
# Authors: Bartolucci, F., Pennoni, F. and Mira, A. #
#######################################################
rm(list = ls())
require("splines")
require("extraDistr")
source("dirmultAR_mcmcFB.R")
source("multinomAR_mcmcFB.R")

#### Load Italian data ####

Y <- read.csv("Italy.csv")[,-1]

#### Set LAmax ####

ORmax <- matrix(c(NA, 10^-7, 0.001, 0.0001, 10^-6, 10^-7,
                 NA, NA,    0.001, 0.0001, 10^-6, 10^-7,
                 NA, 0.1,      NA,    0.1, 10^-5, 10^-6,
                 NA, 0.1,     0.1,     NA,   0.1,  0.01,
                 NA, 10^-7, 10^-7,   0.25,    NA,  0.25,
                 NA, NA,       NA,     NA,    NA,    NA),
               6,byrow=TRUE)

#### Define burn in and iterations  ####
burnin <- 100*10^3; R <- 500*10^3

#### Estimate Models form 1 to 4 of the paper ####
est1 <- multinomAR_mcmc(Y[1:61,], tahead=50, tint=c(7,20), burnin=burnin, R=R, mra=50, degree=2)
est2 <- multinomAR_mcmc(Y[1:61,], tahead=50, tint=c(7,20), burnin=burnin, R=R, mra=50, degree=2, LAmax=ORmax)
est3 <- multinomAR_mcmc(Y[1:61,], tahead=50, tint=c(7,20), burnin=burnin, R=R, mra=50, degree=3)
est4 <- multinomAR_mcmc(Y[1:61,], tahead=50, tint=c(7,20), burnin=burnin, R=R, mra=50, degree=3, LAmax=ORmax)


#### Estimate Models form 5 to 8 of the paper ####
est5 <- dirmultAR_mcmc(Y[1:61,], tahead=50, tint=c(7,20), burnin=burnin, R=R, mra=50, degree=2,Y_new=Y[62:71,])
est6 <- dirmultAR_mcmc(Y[1:61,], tahead=50, tint=c(7,20), burnin=burnin, R=R, mra=50, degree=2,ORmax=ORmax,Y_new=Y[62:71,])
est7 <- dirmultAR_mcmc(Y[1:61,], tahead=50, tint=c(7,20), burnin=burnin, R=R, mra=50, degree=3,Y_new=Y[62:71,])
est8 <- dirmultAR_mcmc(Y[1:61,], tahead=50, tint=c(7,20), burnin=burnin, R=R, mra=50, degree=3,ORmax=ORmax,Y_new=Y[62:71,])

