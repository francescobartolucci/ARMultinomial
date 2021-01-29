#### Results of the models proposed in the  paper titled: ####
# A multivariate statistical approach to predict 
# COVID-19 count data with epidemiological interpretation and uncertainty quantification
# Authors: Bartolucci, F., Pennoni, F. and Mira, A. #
#######################################################

####  Table 3  ####
A <- rbind(
  cbind(mean(est5$CHI2),# hat(Dist)
        mean(est5$CHI2_SIM), # tilde(Dist)
        est5$ppvalue), 
  cbind(mean(est6$CHI2),
        mean(est6$CHI2_SIM),
        est6$ppvalue), 
  cbind(mean(est7$CHI2),
        mean(est7$CHI2_SIM),
        est7$ppvalue),
  cbind(mean(est8$CHI2),
        mean(est8$CHI2_SIM),
        est8$ppvalue)
); A


#### Table 4 realized and predicted discrepancy measures####
P <-cbind(est8$CHI2[1:10], #hat Dist_t
          est8$CHI2_new[1:10],# tilde Dist_t
          est8$ppvalue_new[1:10])
rownames(P)<-c("25th April",
                "26th April",
                "27th April",
                "28th April",
                "29th April",
                "30th April",
                "1st May",
                "2nd May",
                "3rd May",
                "4th May"); P


#### Table 4 posterior p-values ####
PP <-cbind(est8$CHI2[1:10], #hat Dist_t
           est8$CHI2_new[1:10],# tilde Dist_t
           est8$ppvalue_new[1:10])
rownames(PP)<-c("25th April",
                 "26th April",
                 "27th April",
                 "28th April",
                 "29th April",
                 "30th April",
                 "1st May",
                 "2nd May",
                 "3rd May",
                 "4th May"); PP
PP

#### Table 6 Estimated posterior predicted transitions  ####
tab6 <- est8$mTAB_new[,,1]; tab6
#
#### Table 7 Estimated posterior 95% prediction upper and lower bounds for the transitions ####
#
cbind(est8$lwTAB_new[1:6,2,1],
      est8$upTAB_new[1:6,2,1]);
cbind(est8$lwTAB_new[1:6,3,1],
      est8$upTAB_new[1:6,3,1]);
cbind(est8$lwTAB_new[1:6,4,1],
      est8$upTAB_new[1:6,4,1]);
cbind(est8$lwTAB_new[1:6,5,1],
      est8$upTAB_new[1:6,5,1]);
cbind(est8$lwTAB_new[1:6,6,1],
      est8$upTAB_new[1:6,6,1]);

#### Table 8 Estimated posterior mean increase ####
pred8 <- est8$mDiff_new[1:10,4:5]
ci1_pred8<-est8$lwDiff_new[1:10,4:5] 
ci2_pred8<-est8$upDiff_new[1:10,4:5] 
H<-cbind(pred8[,1], ci1_pred8[,1],ci2_pred8[,1],
         pred8[,2], ci1_pred8[,2],ci2_pred8[,2])


#### Table 8 Estimated posterior increase 95% prediction interval ####

pred8<-est8$mDiff_new[1:10,4:5]
ci1_pred8<-est8$lwDiff_new[1:10,4:5] 
ci2_pred8<-est8$upDiff_new[1:10,4:5] 
H<-cbind(pred8[,1], ci1_pred8[,1],ci2_pred8[,1],
         pred8[,2], ci1_pred8[,2],ci2_pred8[,2])

#### FIGURE 1 Observed  and Predicted frequencies ####

ytick<-c(0,22000,45000,65500,90000, 108300)
xtick<-c(1,11,21,31,41,51,61,71)
#
categories<-c(
  "R",
  "Q",
  "H", 
  "ICU",
  "D", 
  "NP")
Tmp <- rbind(cbind(est8$Y[,-1],
                  now_pos = rowSums(est8$Y[,3:5])),
            est8$mTot_new[1:10,-c(1,8)])

ylim = c(0,max(max(Tmp),max(est8$upTot_new[1:10,-c(1,8)])))
col=c(14:19,1) 
col = c("royalblue2",
        "red",
        "olivedrab4",
        "gold2",
        "gray19",
        "wheat4")
op <- par(mar = c(5,6,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1
matplot(Tmp,
        ylab=expression("Predicted frequency ("*italic(hat(y)[tk]^(s))*")"),
        xlab=expression("Day ("*italic(t)*")"), 
        col = col,
        type="l",
        ylim=ylim, 
        lwd=2, 
        lty=1,
        cex.axis = 1, 
        cex.lab = 0.85,
        yaxt="n",
        xaxt="n")
#
abline(v=62,lty = 3)# Fig3
#

axis(side=2, at=ytick, labels = FALSE,
     cex.lab = 0.5, tck=-0.009)
text(par("usr")[1], 
     ytick, 
     # line=1,
     labels = ytick, 
     srt = 45, pos = 2, xpd = TRUE)
#
axis(side=1, at=xtick, labels = FALSE, 
     cex=0.5, 
     padj = 5, tck=-0.009)
text(x=xtick,  
     # line = 2,
     par("usr")[3], 
     labels = xtick, 
     pos = 1, xpd = TRUE)


legend("topleft",
       legend = categories, 
       col = col, 
       pch = 20, bty = "n", 
       x.intersp = 0.8,
       cex= 1,  pt.cex = 1,
       xpd = TRUE, text.width = 0.003)
par(op)
TT = nrow(est8$Y); TT
#
tahead = 10
ind = c((TT+1):(TT+tahead),(TT+tahead):(TT+1)); ind
for(j in 1:6) {
  polygon(ind,c(est8$lwTot_new[1:tahead,j+1],
                rev(est8$upTot_new[1:tahead,j+1])),
          col="gray80",border = NA, lwd=5)
}
for(j in 1:6) lines(Tmp[,j],col=col[j],lwd = 2)

cbind(est8$lwTot_new[1:10,5],est8$upTot_new[1:10,5])

#### F2gure 2 Plot Rt #####

m <- c(est8$mRt,est8$mRt_new)
m<-m[1:71]
TT = nrow(est8$Y); TT
lw = c(est8$lwRt,est8$lwRt_new)
lw<-lw[1:71]
up = c(est8$upRt,est8$upRt_new)
up<-up[1:71];
le = length(m)
ylim = c(0,max(up+0.5,na.rm=TRUE))
xlim = c(1,73)
xtick<-c(1,11,21,31,41,51,61,71)
op <- par(mar = c(5,6,4,2) + 0.1) 
plot(m,
     type="l",
     xaxt="n",
     xlim = xlim,
     xlab=expression("Day ("*italic(t)*")"),
     ylab=expression("Estimated reproduction number ("*italic(hat(R)[t]^(s))*")"),
     ylim=ylim)
axis(side=1, at=xtick, labels = FALSE, 
     cex=0.5, 
     padj = 5, tck=-0.009)
text(x=xtick,  
     # line = 2,
     par("usr")[3], 
     labels = xtick, 
     pos = 1, xpd = TRUE)

polygon(c(1:le,le:1),c(lw,rev(up)),col="gray80",border = NA)
lines(m,type="l")
le1 = length(est4$mRt)
abline(v=61, lty = 3)
