plot.multinomAR <- function(out,ylim=NULL,type=c("tot","diff","Rt"),...){

# plot object of class mulinomAR
#
# INPUT:
# out  = output from multinomAR_mcmc
# ylim = limits of the y-axis
# type = tot to represent total cases and diff to represent the increase in the total cases
  
#preliminaries
  type = match.arg(type)
  TT = nrow(out$Y)
  tahead = nrow(out$mTot_new)
  categories = colnames(out$mTot_new)

# plot
  if(type=="tot"){
    Tmp = rbind(cbind(out$Y[,-1],
                      now_pos = rowSums(out$Y[,3:5]),
                      tot_pos = rowSums(out$Y[,-1])),
                      out$mTot_new[,-1])
    if(is.null(ylim)) ylim = c(0,max(max(Tmp),max(out$upTot_new[,-1])))
    plot(Tmp[,6],xlab=expression("day ("*italic(t)*")"),
         ylab=expression("predicted frequency ("*italic(hat(y)[tk])*")"),type="l",ylim=ylim,...)
    ind = c((TT+1):(TT+tahead),(TT+tahead):(TT+1))
    col = c(2:7,1)
    legend("topleft", legend = categories[-1], col = col, pch = 19, bty = "n")
    for(j in 1:7) polygon(ind,c(out$lwTot_new[,j+1],rev(out$upTot_new[,j+1])),col="gray80",border = NA)
    for(j in 1:7) lines(Tmp[,j],col=col[j])
    lines(c(TT+0.5,TT+0.5),ylim,lty=2)
  }else if(type=="diff"){
    Tmp = rbind(diff(cbind(out$Y[,-1],now_pos = rowSums(out$Y[,3:5]),tot_pos = rowSums(out$Y[,-1]))),
                out$mDiff_new[,-1])
    if(is.null(ylim)) ylim = c(min(min(Tmp),min(out$lwDiff_new[,-1])),
                               max(max(Tmp),max(out$upDiff_new[,-1])))
    plot(2:(nrow(Tmp)+1),Tmp[,6],xlab="day",ylab="frequency",type="l",ylim=ylim,...)
    ind = c((TT+1):(TT+tahead),(TT+tahead):(TT+1))
    col = c(2:7,1)
    legend("topleft", legend = categories[-1], col = col, pch = 19, bty = "n")
    for(j in 1:7) polygon(ind,c(out$lwDiff_new[,j+1],rev(out$upDiff_new[,j+1])),col="gray80",border = NA)
    for(j in 1:7) lines(2:(nrow(Tmp)+1),Tmp[,j],col=col[j])
  }else if(type=="Rt"){
    m = c(out$mRt,out$mRt_new)
    lw = c(out$lwRt,out$lwRt_new)
    up = c(out$upRt,out$upRt_new)
    le = length(m)
    ylim = c(0,max(up,na.rm=TRUE))
    plot(m,type="l",xlab=expression("day ("*italic(t)*")"),
    ylab=expression("reproduction number ("*italic(R[t])*")"),ylim=ylim)
    polygon(c(1:le,le:1),c(lw,rev(up)),col="gray80",border = NA)
    lines(m,type="l")
    le1 = length(out$mRt)
    lines(c(TT+0.5,TT+0.5),ylim,lty=2)
  }
}