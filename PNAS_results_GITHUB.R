#### RESULTS of the PNAS paper titled: ####
# A multivariate statistical model to predict 
# COVID-19 count data with epidemiological interpretation and uncertainty quantification
# Authors: Bartolucci, F., Pennoni, F. and Mira, A. #
#######################################################
#### Results for Italy  periods 1,2,3 ####
rm(list = ls())
load("example_Italy.Rdata")
#### LAmax####
LAmax
#
####  TABLE 1 #######
#
# Realized values of the Dist_k measure
chisq1<-(Y[12:16,]-est1$mTot_new[1:5,])^2/(est1$mTot_new[1:5,])
TOTchisq1<-apply(chisq1,2,sum); round(TOTchisq1,0)
chi1 <-c(TOTchisq1,sum(TOTchisq1)); chi1
#
chisq2<-(Y[32:36,]-est2$mTot_new[1:5,])^2/(est2$mTot_new[1:5,])
TOTchisq2<-apply(chisq2,2,sum); round(TOTchisq2,0)
chi2 <-c(TOTchisq2,sum(TOTchisq2))
#
chisq3<-(Y[62:66,]-est3$mTot_new[1:5,])^2/(est3$mTot_new[1:5,])
TOTchisq3<-apply(chisq3,2,sum); round(TOTchisq2,0)
chi3 <-c(TOTchisq3,sum(TOTchisq3))
Tab1<-rbind(chi1, chi2, chi3); round(Tab1,0)
#
#### TABLE 2 #####
#
# Estimated posterior predicted transitions
tab <- est1$mTAB_new[,,1]; round(tab,0)
#
#### TABLE 3  ####
# Estimated posterior 95\% CI
i<-5
citot1<-cbind(est1$lwTAB_new[i,1:6,1],
              est1$upTAB_new[i,1:6,1]);round(citot1,0)
#
#### TABLE 4 #### 
# Estimated predicted increase
#
pred<-est1$mDiff_new[1:5,4:5];
ci1_pred1<-est1$lwDiff_new[1:5,4:5]; 
ci2_pred1<-est1$upDiff_new[1:5,4:5]; 
#
H<-cbind(pred[,1], ci1_pred1[,1],ci2_pred1[,1],
         pred[,2], ci1_pred1[,2],ci2_pred1[,2])
round(H,0)
#
#### TABLE 5 #####
#
# Estimated posterior predicted transitions
tab <- est2$mTAB_new[,,1]; round(tab,0)

#
#### TABLE 6 #### 
#
#predictions of the increase  for H and ICU
pred2<-est2$mDiff_new[1:5,4:5]
ci1_pred2<-est2$lwDiff_new[1:5,4:5]
ci2_pred2<-est2$upDiff_new[1:5,4:5]
#
H2<-cbind(pred2[,1], ci1_pred2[,1],ci2_pred2[,1],
          pred2[,2], ci1_pred2[,2],ci2_pred2[,2])
round(H2,0)
#
#### TABLE 7 #####
#
# Estimated posterior predicted transitions
tab <- est3$mTAB_new[,,1]; round(tab,0)

#
#### TABLE 8 #### 
# predictions of the increase  for H and ICU
#
pred3<-est3$mDiff_new[1:5,4:5];round(pred3,0)
ci1_pred3<-est3$lwDiff_new[1:5,4:5]; 
ci2_pred3<-est3$upDiff_new[1:5,4:5];
#
H3<-cbind(pred3[,1], ci1_pred3[,1],ci2_pred3[,1],
          pred3[,2], ci1_pred3[,2],ci2_pred3[,2]);
round(H3,0)

