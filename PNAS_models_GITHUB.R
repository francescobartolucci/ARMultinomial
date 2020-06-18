#### Models of the PNAS paper titled: ####
# A multivariate statistical model to predict 
# COVID-19 count data with epidemiological interpretation and uncertainty quantification
# Authors: Bartolucci, F., Pennoni, F. and Mira, A. #
#######################################################
# load packages and source functions
rm(list = ls())
require(splines)
source("multinomAR_mcmcFB.R")
source("print.multinomAR.R")
source("plot.multinomAR.R")

# load data for Italy
Y = read.csv("Italy.csv")[,-1]

# Set LAmax
LAmax = matrix(c(NA, 10^-7, 0.001, 0.0001, 10^-6, 10^-7,
                 NA, NA,    0.001, 0.0001, 10^-6, 10^-7,
                 NA, 0.1,      NA,    0.1, 10^-5, 10^-6,
                 NA, 0.1,     0.1,     NA,   0.1,  0.01,
                 NA, 10^-7, 10^-7,   0.25,    NA,  0.25,
                 NA, NA,       NA,     NA,    NA,    NA),
               6,byrow=TRUE)

#
# Estimate the model for period 1 and  prediction 5 days ahead
#
est1 = multinomAR_mcmc(Y[1:11,], tahead=5, 
                        burnin=100*10^3, 
                        R=500*10^3, 
                        mra=50, 
                        LAmax=LAmax)

#
# Estimate the model for period 2 and  prediction 5 days ahead
#
est2 = multinomAR_mcmc(Y[1:31,], 
                        tahead=50, 
                        burnin=100*10^3, 
                        R=500*10^3, 
                        mra=50, LAmax=LAmax)

#
# Estimate the model for period 3 with the effects of 2 interventions at day 7 and 20
# and  prediction 5 days ahead
#
est3 = multinomAR_mcmc(Y[1:61,], 
                        tahead=50, 
                        tint=c(7,20), 
                        burnin=100*10^3, 
                        R=500*10^3, mra=50, 
                        LAmax=LAmax)
#
#save.image("example_Italy.RData")
#
plot(est1)   # plot_pred1
plot(est1,type="Rt")   # plot_Rt1

