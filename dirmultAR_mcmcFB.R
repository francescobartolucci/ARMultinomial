dirmultAR_mcmc <- function(Y,R=10000,burnin=1000,mBE=NULL,VBE=NULL,Tau=matrix(0.1,6,6),
                             mra=10,tint=NULL,tahead=1,output=FALSE,disp=TRUE,tin=10,
                             ORmin=matrix(NA,6,6),ORmax=matrix(NA,6,6),conbe=FALSE,
                             degree=2,Y_new=NULL){

  
# WITH DIRICHLET-MULTINOMIAL MODEL with dispersion parameters

# Estimate multinomial AR model for COVID-19 data by an MCMC algorithm
# Categories must be ordered as:
# "susceptible", "recovered", "quarantene", "hospitalized", "intensive_care", "deceased"
#
# INPUT:
# Y      = matrix (TT x 6) of observed daily frequencies for the 6 categories ordered from suscptibles do died
# burnin = number of iterations of burnin
# mBE    = array of means of priors on regression parameters
# VBE    = array of variance-covarince matrices of priors on regression parameters
# Tau    = matrix of standard deviations for each proposal
# mra    = maximum variation of each cell in updating the tables
# tint   = times of interventions (NA for no intervention)
# thaed  = numer of time ahead for prediction
# output = TRUE for providing all estimated matrices at each MCMC iteration
# disp   = to disply prediction during estimation
# tin    = store results each tin iterations
# ORmin  = matrix of lower limits of OR
# ORmax = matrix of upper limits of OR
# conbe  = contraint 2° order coefficients in the polinomials to be non-positive
#
# OUTPUT:
# mTAB       = posterior mean of transition tables
# lwTAB      = posterior 95% CI lwb of transition tables
# upTAB      = posterior 95% CI upb of transition tables
# mBE        = posterior mean of regression coefficients
# VBE        = posterior variance of regression coefficients
# lwBE       = posterior 95% CI lwb of regression coefficients
# upBE       = posterior 95% CI upb of regression coefficients
# mTAB_new   = posterior mean of predicted transition tables
# lwTAB_new  = posterior 95% CI lwb of predicted transition tables
# upTAB_new  = posterior 95% CI upb of predicted transition tables
# mTot_new   = posterior mean of predicted totals
# lwTot_new  = posterior 95% CI lwb  of predicted totals
# upTot_new  = posterior 95% CI upb  of predicted totals
# mDiff_new  = posterior mean of predicted increase in totals
# lwDiff_new = posterior 95% CI lwb  of predicted increase in totals
# upTot_new  = posterior 95% CI upb  of predicted increase in totals
# acctab     = overal acceptance rate for tables
# accbe      = overal acceptance rate for regression parameters
# Accbe      = acceptance rate for each regression parameter vector
# Y          = matrix of observed daily frequencies
# TTAB       = array of transition tables for each MCMC iteration
# BBE        = array of regression coefficients for each MCMC iteration
# TTAB_new   = array of predicted transition matrices for each MCMC iteration
# TOT_new    = array of predicted totals for each MCMC iteration
# DIFF_new   = array of predicted increase of totals for each MCMC iteration
  
# prelinaries
  Y = as.matrix(Y)
  TT = nrow(Y)
  categories = colnames(Y)
  if(is.null(mBE)){
    ldnorm1 <- function(x1,x2,si2=1) -sum(x1^2-x2^2)/(2*si2)
  }else{
  	iVBE = array(0,dim(VBE))
  	for(i in 1:6) for(j in 1:6) iVBE[i,j,,] = solve(VBE[i,j,,])
  	ldmvnorm1 <- function(x1,x2,mu=rep(0,length(x1)),iSi=diag(length(x1))) 
  	               -c((x1-mu)%*%iSi%*%(x1-mu)-(x2-mu)%*%iSi%*%(x2-mu))/2	
  }

# initial tables
  TAB = array(0,c(6,6,TT))
  TAB[,,1] = NA
  for(t in 2:TT){
    mi = round(pmin(Y[t,],Y[t-1,])*c(1,0.9,0.9,0.9,0.9,1))
    rtot = Y[t-1,]-mi
    ctot = Y[t,]-mi
    TAB[,,t] = diag(mi)+r2dtable(1,rtot,ctot)[[1]]
  }

# design matrices
  if(is.null(tint)){
    XX = cbind(1,bs(1:(TT-1+tahead),degree=degree))
  }else{
    XX = cbind(1,bs(1:(TT-1+tahead),degree=degree,knots = tint-1))
  }
  XXC = XX
  XX_new = XX[-(1:(TT-1)),]
  XX = XX[1:(TT-1),]

# initial values of the parameters and multinomial probabilities
  mTAB = apply(TAB[,,2:TT],c(1,2),mean)
  nbe = degree+1+length(tint)
  BE = array(0,c(6,6,nbe))
  Ind = rbind(cbind(1,1:6),cbind(2,2:6),cbind(3,2:6),cbind(4,2:6),cbind(5,2:6),c(6,6))
  BE[2:6,1,] = BE[6,2:5,] = NA
  for(h in 1:27){
    i = Ind[h,1]; j = Ind[h,2]
    if(j!=i){
      if(!is.null(mBE)) BE[i,j,] = mBE[i,j,]
      if(!is.na(ORmin[i,j]) & is.na(ORmax[i,j])) BE[i,j,1] = max(BE[i,j,1],log(ORmin[i,j]))
      if(is.na(ORmin[i,j]) & !is.na(ORmax[i,j])) BE[i,j,1] = min(BE[i,j,1],log(ORmax[i,j]))
      if(!is.na(ORmin[i,j]) & !is.na(ORmax[i,j])){
        if(BE[i,j,1]<log(ORmin[i,j]) | BE[i,j,1]>log(ORmax[i,j])) BE[i,j,1] = log((ORmin[i,j]+ORmax[i,j])/2)
      }
    }
  }
  LA = PP = array(0,c(6,6,TT))
  LA[,,1] = PP[,,1] = NA
  for(h in 1:27) LA[Ind[h,1],Ind[h,2],2:TT] = exp(XX%*%BE[Ind[h,1],Ind[h,2],])
  for(t in 2:TT) PP[,,t] = (1/rowSums(LA[,,t]))*LA[,,t]
  LAC = array(0,c(6,6,TT+tahead))
  LAC[,,1] = NA
  for(h in 1:27) LAC[Ind[h,1],Ind[h,2],2:(TT+tahead)] = exp(XXC%*%BE[Ind[h,1],Ind[h,2],])
  ORC = LAC
  if(is.null(mBE)) si2 = 100

# iterate
  acctab = accbe = accnu = 0
  Accbe = matrix(0,5,6)
  Accbe[2:5,1] = NA
  TTAB = array(as.integer(0),c(6,6,TT,R/tin))
  TTAB[,,,1] = NA
  BBE = array(0,c(6,6,nbe,R/tin))
  BBE[2:6,1,,] = BBE[6,2:5,,] = NA
  RT = matrix(0,TT,R/tin)
  TTAB_new = array(as.integer(0),c(6,6,tahead,R/tin))
  TOT_new = DIFF_new = array(as.integer(0),c(tahead,8,R/tin))
  RT_new = matrix(0,tahead,R/tin)
  CHI2 = CHI2_SIM = rep(0,R/tin)
  if(!is.null(Y_new)){
    tahead1 = nrow(Y_new)
    CHI2_new = CHI2_SIM_new = matrix(0,R/tin,tahead1)
    Y_new = as.matrix(Y_new)
  }
  it = 0
  t0 = proc.time()[3]; names(t0) = NULL
  seqra = c(-mra:-1,1:mra)
  for(r in -(burnin-1):R){
    it = it+1

# update table
    for(t in 2:TT){
      La = LA[,,t]; P = PP[,,t]
      for(it1 in 1){
      	tmp = sample(1:5,2); i1 = min(tmp); i2 = max(tmp)
      	tmp = sample(2:6,2); j1 = min(tmp); j2 = max(tmp)
        ra = sample(seqra,1)
        Tab = Tabs = TAB[,,t]
        Tabs[i1,j1] = Tab[i1,j1]+ra
        Tabs[i1,j2] = Tab[i1,j2]-ra
        Tabs[i2,j1] = Tab[i2,j1]-ra
        Tabs[i2,j2] = Tab[i2,j2]+ra
        # if(t==TT) browser()
        if(all(Tabs>=0)){
          lnum = sum(lgamma(Tabs[i1,c(j1,j2)]+La[i1,c(j1,j2)])-lgamma(Tabs[i1,c(j1,j2)]+1)+
                     lgamma(Tabs[i2,c(j1,j2)]+La[i2,c(j1,j2)])-lgamma(Tabs[i2,c(j1,j2)]+1))
          lden = sum(lgamma(Tab[i1,c(j1,j2)]+La[i1,c(j1,j2)])-lgamma(Tab[i1,c(j1,j2)]+1)+
                     lgamma(Tab[i2,c(j1,j2)]+La[i2,c(j1,j2)])-lgamma(Tab[i2,c(j1,j2)]+1))
          al = exp(lnum-lden)
          # al1 = exp(lnum-lden)
          # print(c(al,al1,al/al1-1))
          if(is.nan(al)){
            print("NA in updating tables")
            browser()
          }
          if(runif(1)<=al){
            TAB[,,t] = Tabs
            acctab = acctab + 1/(1*(TT-1))
          }
        }
      }
    }

# update beta
    for(h in 1:27){
      i = Ind[h,1]; j = Ind[h,2]
      if(i<6){
        ind = Ind[Ind[,1]==i,2]
        BES = BE; LAS = LA; LACS = LAC; ORCS = ORC; PPS = PP
        BES[i,j,] = BE[i,j,]+rnorm(nbe,0,Tau[i,j])
        LAS[i,j,2:TT] = exp(XX%*%BES[i,j,])
        LACS[i,j,2:(TT+tahead)] = exp(XXC%*%BES[i,j,])
        for(t in 2:TT) PPS[i,,t] = LAS[i,,t]/sum(LAS[i,,t])
        for(t in 2:(TT+tahead)) ORCS[i,,t] = LACS[i,,t]/LACS[i,i,t]
        check = TRUE
        if(conbe) if(BES[i,j,3]>0) check=FALSE
        for(j1 in ind) if(j1!=i){
          if(!is.na(ORmin[i,j1])) if(any(ORCS[i,j1,2:(TT+tahead)]<ORmin[i,j1])) check = FALSE
          if(!is.na(ORmax[i,j1])) if(any(ORCS[i,j1,2:(TT+tahead)]>ORmax[i,j1])) check = FALSE
        }
        if(check){
          if(is.null(mBE)){
            tmp = ldnorm1(BES[i,j,],BE[i,j,],si2)
          }else{
            tmp = ldmvnorm1(BES[i,j,],BE[i,j,],mBE[i,j,],iVBE[i,j,,])
          }
          for(t in 2:TT){
            La = LA[,,t]; Las = LAS[,,t]
            tla = sum(La[i,ind]); tlas = sum(Las[i,ind])
            ntab = sum(TAB[i,ind,t])
            tmp = tmp + lgamma(tlas)-lgamma(ntab+tlas)+sum(lgamma(TAB[i,ind,t]+Las[i,ind])-lgamma(Las[i,ind]))-
                        lgamma(tla)+lgamma(ntab+tla)-sum(lgamma(TAB[i,ind,t]+La[i,ind])-lgamma(La[i,ind]))
          }
          al = exp(tmp)
          if(is.nan(al)){
            print("NA in updating Beta")
            browser()
          }
          if(runif(1)<=al){
            BE = BES; LA = LAS; LAC = LACS; ORC = ORCS; PP = PPS
            accbe = accbe + 1/26
            Accbe[i,j] = Accbe[i,j]+1
          }
        }
      }
    }

# Chi2 for prediction error
if((r>0 & r%%tin==0) | r%%250==0){
    chi2 = chi2_sim = rep(0,6)
    for(t in 2:TT){
      tmp = colSums(Y[t-1,]*PP[,,t])
      chi2 = chi2+(Y[t,]-tmp)^2/tmp
      Tmp = matrix(0,6,6)
      for(i in 1:5){
        ind = Ind[Ind[,1]==i,2]
        Tmp[i,ind] = rdirmnom(1,Y[t-1,i],LA[i,ind,t])
      }
      Tmp[6,6] = Y[t-1,6]
      tmp2 = colSums(Tmp)
      chi2_sim = chi2_sim+(tmp2-tmp)^2/tmp
    }
    Chi2 = sum(chi2); Chi2_sim = sum(chi2_sim)
}

# Estimate Rt
    if((r>0 & r%%tin==0) | r%%250==0){
      newill = rep(0,TT+tahead); Rt = rep(0,TT)
      newill[1] = Rt[1] = NA
      for(t in 2:TT){
        newill[t] = sum(Y[t-1,1]*PP[1,-1,t])
        phi = dgamma(1:(t-2),shape=1.87,rate=0.28)
        phi = phi/sum(phi)
        Rt[t] = newill[t]/sum(phi*newill[(t-1):2])
      }
    }

# Prediction
    if((r>0 & r%%tin==0) | r%%250==0){
      LA_new = PP_new = TAB_new = array(0,c(6,6,tahead))
      Tot_new = Diff_new = matrix(0,tahead,8)
      Rt_new = rep(0,tahead)
      for(h in 1:27) LA_new[Ind[h,1],Ind[h,2],] = exp(XX_new%*%BE[Ind[h,1],Ind[h,2],])
      for(t in 1:tahead) PP_new[,,t] = (1/rowSums(LA_new[,,t]))*LA_new[,,t]
      for(i in 1:5){
        ind = Ind[Ind[,1]==i,2]
        TAB_new[i,ind,1] = rdirmnom(1,Y[TT,i],LA_new[i,ind,1])
      }
      TAB_new[6,6,1] = Y[TT,6]
      Tot_new[1,1:6] = colSums(TAB_new[,,1])
      Tot_new[1,7] = sum(Tot_new[1,3:5])
      Tot_new[1,8] = sum(Tot_new[1,2:6])
      if(!is.null(Y_new)){
        Chi2_new = Chi2_sim_new = rep(0,tahead1)
        pred_new = colSums(Y[TT,]*PP_new[,,1])
        Chi2_new[1] = sum((Y_new[1,]-pred_new)^2/pred_new)
        Chi2_sim_new[1] = sum((Tot_new[1,1:6]-pred_new)^2/pred_new)
      }
      Diff_new[1,] = Tot_new[1,]-c(Y[TT,],sum(Y[TT,3:5]),sum(Y[TT,-1]))
      t1 = TT+1
      newill[t1] = sum(Y[TT,1]*PP_new[1,-1,1])
      phi = dgamma(1:(t1-2),shape=1.87,rate=0.28)
      phi = phi/sum(phi)
      Rt_new[1] = newill[t1]/sum(phi*newill[(t1-1):2])
      if(tahead>1) for(t in 2:tahead){
        for(i in 1:5){
          ind = Ind[Ind[,1]==i,2]
          TAB_new[i,ind,t] = rdirmnom(1,Tot_new[t-1,i],LA_new[i,ind,t])
        }
        TAB_new[6,6,t] = Tot_new[t-1,6]
        Tot_new[t,1:6] = colSums(TAB_new[,,t])
        Tot_new[t,7] = sum(Tot_new[t,3:5])
        Tot_new[t,8] = sum(Tot_new[t,2:6])
        if(!is.null(Y_new)) if(t<=tahead1){
          pred_new = colSums(pred_new*PP_new[,,t])
          Chi2_new[t] = sum((Y_new[t,]-pred_new)^2/pred_new)
          Chi2_sim_new[t] = sum((Tot_new[t,1:6]-pred_new)^2/pred_new)
        }
        Diff_new[t,] = Tot_new[t,]-Tot_new[t-1,]
        t1 = TT+t
        newill[t1] = sum(Tot_new[t-1,1]*PP_new[1,-1,t])
        phi = dgamma(1:(t1-2),shape=1.87,rate=0.28)
        phi = phi/sum(phi)
        Rt_new[t] = newill[t1]/sum(phi*newill[(t1-1):2])
      }
    }

# Store values
    if(r>0 & r%%tin==0){
      r1 = round(r/tin)
      TTAB[,,,r1] = TAB
      BBE[,,,r1] = BE
      RT[,r1] = Rt
      TTAB_new[,,,r1] = TAB_new
      TOT_new[,,r1] = Tot_new
      DIFF_new[,,r1] = Diff_new
      RT_new[,r1] = Rt_new
      CHI2[r1] = Chi2
      CHI2_SIM[r1] = Chi2_sim
      if(!is.null(Y_new)){
        CHI2_new[r1,] = Chi2_new
        CHI2_SIM_new[r1,] = Chi2_sim_new
      }
    }

# display output
    if(r%%250==0){
      tt = proc.time()[3]; names(tt) = NULL
      print(c(iteration=r,acc_tab=acctab/it,acc_be=accbe/it,time100=100*(tt-t0)/it))
      if(disp){
        print(Accbe)
        cat("\n")
      	print("last matrix of OR")
        print(round(ORC[,,TT],6))
      	cat("\n")
      	print("last table")
        print(cbind(TAB[,,TT],NA,round(Y[TT-1,]*PP[,,TT])))
      	cat("\n")
      	print("Predicted totals last occasion")
      	print(rbind(Y[TT,],tmp))
      	cat("\n")
      	print("chi2 for prediction error")
      	print(round(chi2,2))
      	print("Chi2 for observed and predicted data")
      	print(round(c(Chi2,Chi2_sim),2))
      	if(!is.null(Y_new)){
        	print("Chi2 for observed and predicted data (forecast)")
        	print(round(cbind(Chi2_new,Chi2_sim_new),2))
      	}
      	if(r>=tin){
          print("mean of chi2 and posterior p-value")
          print(c(mean(CHI2[1:r1]),mean(CHI2_SIM[1:r1])))
      	  print(mean(CHI2_SIM[1:r1]>CHI2[1:r1]))
      	  if(!is.null(Y_new)){
        	  print("mean of chi2 and posterior p-value (forecast)")
            Tmp = NULL
            for(t in 1:tahead1) Tmp = rbind(Tmp,c(mean(CHI2_new[1:r1,t]),mean(CHI2_SIM_new[1:r1,t]),
                                                  mean(CHI2_SIM_new[1:r1,t]>CHI2_new[1:r1,t])))
            print(Tmp)
      	  }
      	}
      	cat("\n")
        cat("Predicted Rt\n")
        print(head(Rt))
        cat("\n")
      }
      if(r>tin & disp){
        mTot_new = apply(TOT_new[,,1:round(r/tin)],1:2,mean)
        colnames(mTot_new) = c(categories,"now_pos","tot_pos")
        cat("Predicted totals\n")
        print(head(mTot_new))
        cat("\n")
        cat("Predicted_new Rt\n")
        print(head(Rt_new))
        cat("\n")
      }
    }
  }

# output
  mTAB = apply(TTAB,1:3,mean)
  lwTAB = apply(TTAB,1:3,quantile,probs=0.025,na.rm=TRUE)
  upTAB = apply(TTAB,1:3,quantile,probs=0.975,na.rm=TRUE)
  dimnames(mTAB) = dimnames(lwTAB) = dimnames(upTAB) = list(categories,categories,NULL)
  mBE = apply(BBE,1:3,mean)
  VBE = array(0,c(6,6,nbe,nbe))
  for(i in 1:6) for(j in 1:6) VBE[i,j,,] = cov(t(BBE[i,j,,]))
  lwBE = apply(BBE,1:3,quantile,probs=0.025,na.rm=TRUE)
  upBE = apply(BBE,1:3,quantile,probs=0.975,na.rm=TRUE)
  mTAB_new = apply(TTAB_new,1:3,mean)
  lwTAB_new = apply(TTAB_new,1:3,quantile,probs=0.025,na.rm=TRUE)
  upTAB_new = apply(TTAB_new,1:3,quantile,probs=0.975,na.rm=TRUE)
  dimnames(mTAB_new) = dimnames(lwTAB_new) = dimnames(upTAB_new) = list(categories,categories,NULL)
  mRt = rowMeans(RT)
  lwRt = apply(RT,1,quantile,probs=0.025,na.rm=TRUE)
  upRt = apply(RT,1,quantile,probs=0.975,na.rm=TRUE)
  mTot_new = apply(TOT_new,1:2,mean)
  lwTot_new = apply(TOT_new,1:2,quantile,probs=0.025,na.rm=TRUE)
  upTot_new = apply(TOT_new,1:2,quantile,probs=0.975,na.rm=TRUE)
  colnames(mTot_new) = colnames(lwTot_new) = colnames(upTot_new) = c(categories,"now_pos","tot_pos")
  mDiff_new = apply(DIFF_new,1:2,mean)
  lwDiff_new = apply(DIFF_new,1:2,quantile,probs=0.025,na.rm=TRUE)
  upDiff_new = apply(DIFF_new,1:2,quantile,probs=0.975,na.rm=TRUE)
  colnames(mDiff_new) = colnames(lwDiff_new) = colnames(upDiff_new) = c(categories,"now_pos","tot_pos")
  mRt_new = rowMeans(RT_new)
  lwRt_new = apply(RT_new,1,quantile,probs=0.025,na.rm=TRUE)
  upRt_new = apply(RT_new,1,quantile,probs=0.975,na.rm=TRUE)
  out = list(mTAB=mTAB,lwTAB=lwTAB,upTAB=upTAB,mBE=mBE,lwBE=lwBE,upBE=upBE,VBE=VBE,
             mRt=mRt,lwRt=lwRt,upRt=upRt,
             mTAB_new=mTAB_new,lwTAB_new=lwTAB_new,upTAB_new=upTAB_new,
             mTot_new=mTot_new,lwTot_new=lwTot_new,upTot_new=upTot_new,
             mDiff_new=mDiff_new,lwDiff_new=lwDiff_new,upDiff_new=upDiff_new,
             mRt_new=mRt_new,lwRt_new=lwRt_new,upRt_new=upRt_new,
             acctab=acctab/(R+burnin),accbe=accbe/(R+burnin),Accbe=Accbe/(R+burnin),
             CHI2=CHI2,CHI2_SIM=CHI2_SIM,
             CHI2_new = CHI2_new,
             ppvalue_new =  colMeans(CHI2_SIM_new>CHI2_new),
             ppvalue = mean(CHI2_SIM>CHI2),
             Y=Y,call = match.call())
  if(is.null(Y_new)){
    out$CHI2_new=CHI2_new; out$CHI2_SIM_new=CHI2_SIM_new
    out$ppvalue_new = colMeans(CHI2_SIM_new>CHI2_new)
  }
  if(output){
    out$TTAB = TTAB; out$BBE = BBE; out$RT = RT; out$TTAB_new = TTAB_new
    out$TOT_new = TOT_new; out$RT_new = RT_new; out$DIFF_new = DIFF_new
  }
  class(out) = "dirmultAR"
  return(out)

}