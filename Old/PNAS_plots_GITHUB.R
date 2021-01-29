#### PLOTS of the PNAS paper titled: ####
# A multivariate statistical model to predict 
# COVID-19 count data with epidemiological interpretation and uncertainty quantification
# Authors: Bartolucci, F., Pennoni, F. and Mira, A. #
#######################################################
#### Plot RESULTS for Italy first period ####
#### Author: Fulvia Pennoni ####
rm(list = ls())
load("example_Italy.Rdata")
#### define plots ####
ytick<-c(0,2000,5000,8000)
xtick<-c(1,6,11,16)
#
categories<-c(
  "R",
  "Q",
  "H", 
  "ICU",
  "D", 
  "NP")
Tmp = rbind(cbind(est1$Y[,-1],
                  now_pos = rowSums(est1$Y[,3:5])),
            est1$mTot_new[1:5,-c(1,8)])
ylim = c(0,max(max(Tmp),max(est1$upTot_new[1:5,-c(1,8)])))
col=c(14:19,1) 
col = c("royalblue2",
        "red",
        "olivedrab4",
        "gold2",
        "gray19",
        "wheat4")
matplot(Tmp,
        ylab=expression("Predicted frequency ("*italic(hat(y)[tk])*")"),
        xlab=expression("Day ("*italic(t)*")"), 
        col = col,
        type="l",
        ylim=ylim, 
        lwd=2, 
        lty=1,
        cex.axis = 1, 
        cex.lab = 1,
        yaxt="n",
        xaxt="n")
#
abline(v=11, lty = 3)
#
axis(side=2, at=ytick, labels = FALSE)
text(par("usr")[1], 
     ytick,  
     labels = ytick, srt = 45, pos = 2, xpd = TRUE)
#
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  
     par("usr")[3], 
     labels = xtick, 
     pos = 1, xpd = TRUE)
legend("topleft",
       legend = categories, 
       col = col, 
       pch = 20, bty = "n", 
       x.intersp = 0.5,
       cex= 0.7,  pt.cex = 1,
       xpd = TRUE, text.width = 0.0001)
TT = nrow(est1$Y); TT
tahead = 5
ind = c((TT+1):(TT+tahead),(TT+tahead):(TT+1)); ind
for(j in 1:6) polygon(ind,c(est1$lwTot_new[1:tahead,j+1],
                            rev(est1$upTot_new[1:tahead,j+1])),
                      col="gray80",border = NA, lwd=5)
for(j in 1:6) lines(Tmp[,j],col=col[j],lwd = 2)

#### Plot Rt  period 3 #####
m = c(est3$mRt,est3$mRt_new)
m<-m[1:71]
lw = c(est3$lwRt,est3$lwRt_new)
lw<-lw[1:71]
up = c(est3$upRt,est3$upRt_new)
up<-up[1:71];
le = length(m)
ylim = c(0,max(up,na.rm=TRUE))
plot(m,type="l",xlab=expression("Day ("*italic(t)*")"),
     ylab=expression("Reproduction number ("*italic(R[t])*")"),ylim=ylim)
polygon(c(1:le,le:1),c(lw,rev(up)),col="gray80",border = NA)
lines(m,type="l")
le1 = length(est1$mRt)
lines(c(TT+0.5,TT+0.5),ylim,lty=2)

