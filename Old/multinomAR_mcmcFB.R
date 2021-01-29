multinomAR_mcmc <- function(Y,R=10000,burnin=1000,mBE=NULL,VBE=NULL,Tau=matrix(0.1,6,6),
                            mra=10,tint=NULL,tahead=1,output=FALSE,disp=TRUE,tin=10,
                            LAmin=matrix(NA,6,6),LAmax=matrix(NA,6,6),conbe=FALSE){

  
# ADD TOT E DIFF

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
# LAmin  = matrix of lower limits of OR
# LAmax  = matrix of upper limits of OR
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
    XX = cbind(1,bs(1:(TT-1+tahead),degree=2))
  }else{
    XX = cbind(1,bs(1:(TT-1+tahead),degree=2,knots = tint-1))
  }
  XXC = XX
  XX_new = XX[-(1:(TT-1)),]
  XX = XX[1:(TT-1),]

# initial values of the parameters and multinomial probabilities
  mTAB = apply(TAB[,,2:TT],c(1,2),mean)
  nbe = 3+length(tint)
  BE = array(0,c(6,6,nbe))
  Ind = rbind(cbind(1,1:6),cbind(2,2:6),cbind(3,2:6),cbind(4,2:6),cbind(5,2:6),c(6,6))
  BE[2:6,1,] = BE[6,2:5,] = NA
  for(h in 1:27){
    i = Ind[h,1]; j = Ind[h,2]
    if(j!=i){
      if(!is.null(mBE)) BE[i,j,] = mBE[i,j,]
      if(!is.na(LAmin[i,j]) & is.na(LAmax[i,j])) BE[i,j,1] = max(BE[i,j,1],log(LAmin[i,j]))
      if(is.na(LAmin[i,j]) & !is.na(LAmax[i,j])) BE[i,j,1] = min(BE[i,j,1],log(LAmax[i,j]))
      if(!is.na(LAmin[i,j]) & !is.na(LAmax[i,j])){
        if(BE[i,j,1]<log(LAmin[i,j]) | BE[i,j,1]>log(LAmax[i,j])) BE[i,j,1] = log((LAmin[i,j]+LAmax[i,j])/2)
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
  if(is.null(mBE)) si2 = 100

# iterate
  acctab = accbe = 0
  Accbe = matrix(0,6,6)
  Accbe[2:6,1] = Accbe[6,2:5] = NA
  diag(Accbe) = NA
  TTAB = array(as.integer(0),c(6,6,TT,R/tin))
  TTAB[,,,1] = NA
  BBE = array(0,c(6,6,nbe,R/tin))
  BBE[2:6,1,,] = BBE[6,2:5,,] = NA
  RT = matrix(0,TT,R/tin)
  TTAB_new = array(as.integer(0),c(6,6,tahead,R/tin))
  TOT_new = DIFF_new = array(as.integer(0),c(tahead,8,R/tin))
  RT_new = matrix(0,tahead,R/tin)
  it = 0
  t0 = proc.time()[3]; names(t0) = NULL
  seqra = c(-mra:-1,1:mra)
  for(r in -(burnin-1):R){
    it = it+1

# update table
    for(t in 2:TT){
      P = PP[,,t]
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
          P1 = c(P[i1,j1],P[i1,j2],P[i2,j1],P[i2,j2])
          Tab1 = c(Tab[i1,j1],Tab[i1,j2],Tab[i2,j1],Tab[i2,j2]) 
          Tabs1 = c(Tabs[i1,j1],Tabs[i1,j2],Tabs[i2,j1],Tabs[i2,j2])
          al = exp(sum(log(P1)*(Tabs1-Tab1)+lgamma(Tab1+1)-lgamma(Tabs1+1))) 
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
      if(j!=i){
        BES = BE; LAS = LA; LACS = LAC; PPS = PP
        BES[i,j,] = BE[i,j,]+rnorm(nbe,0,Tau[i,j])
        LAS[i,j,2:TT] = exp(XX%*%BES[i,j,])
        LACS[i,j,2:(TT+tahead)] = exp(XXC%*%BES[i,j,])
        for(t in 2:TT) PPS[i,,t] = LAS[i,,t]/sum(LAS[i,,t])
        check = TRUE
        if(conbe) if(BES[i,j,3]>0) check=FALSE
        if(!is.na(LAmin[i,j])) if(any(LACS[i,j,2:(TT+tahead)]<LAmin[i,j])) check = FALSE
        if(!is.na(LAmax[i,j])) if(any(LACS[i,j,2:(TT+tahead)]>LAmax[i,j])) check = FALSE
        if(check){
          ind = Ind[Ind[,1]==i,2]
          if(is.null(mBE)){
            tmp = ldnorm1(BES[i,j,],BE[i,j,],si2)
          }else{
            tmp = ldmvnorm1(BES[i,j,],BE[i,j,],mBE[i,j,],iVBE[i,j,,])
          }
          for(t in 2:TT){
            P = PP[,,t]; PS = PPS[,,t]
            tmp = tmp+c(TAB[i,ind,t]%*%(log(PS[i,ind])-log(P[i,ind])))
          }
          al = exp(tmp)
          if(is.nan(al)){
            print("NA in updating Beta")
            browser()
          }
          if(runif(1)<=al){
            BE = BES; LA = LAS; LAC = LACS; PP = PPS
            accbe = accbe + 1/21
            Accbe[i,j] = Accbe[i,j]+1
          }
        }
      }
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
      for(i in 1:6) TAB_new[i,,1] = rmultinom(1,Y[TT,i],PP_new[i,,1])
      Tot_new[1,1:6] = colSums(TAB_new[,,1])
      Tot_new[1,7] = sum(Tot_new[1,3:5])
      Tot_new[1,8] = sum(Tot_new[1,2:6])
      Diff_new[1,] = Tot_new[1,]-c(Y[TT,],sum(Y[TT,3:5]),sum(Y[TT,-1]))
      t1= TT+1
      newill[t1] = sum(Y[TT,1]*PP_new[1,-1,1])
      phi = dgamma(1:(t1-2),shape=1.87,rate=0.28)
      phi = phi/sum(phi)
      Rt_new[1] = newill[t1]/sum(phi*newill[(t1-1):2])
      if(tahead>1) for(t in 2:tahead){
        for(i in 1:6) TAB_new[i,,t] = rmultinom(1,Tot_new[t-1,i],PP_new[i,,t])
        Tot_new[t,1:6] = colSums(TAB_new[,,t])
        Tot_new[t,7] = sum(Tot_new[t,3:5])
        Tot_new[t,8] = sum(Tot_new[t,2:6])
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
    }

# display output
    if(r%%250==0){
      tt = proc.time()[3]; names(tt) = NULL
      print(c(iteration=r,acc_tab=acctab/it,acc_be=accbe/it,time100=100*(tt-t0)/it))
      if(disp){
      	cat("\n")
      	print("last matrix of OR")
        print(round(LA[,,TT],6))
      	cat("\n")
      	print("last table")
        print(cbind(TAB[,,TT],NA,round(Y[TT-1,]*PP[,,TT])))
      	cat("\n")
        print("chi2 for prediction error")
        chi2 = 0
        for(t in 2:TT){
          tmp = colSums(Y[t-1,]*PP[,,t])
          chi2 = chi2+(Y[t,]-tmp)^2/tmp
        }
        print(round(chi2,2))
        cat("\n")
        print("Predicted totals last occasion")
        print(rbind(Y[t,],tmp))
        cat("\n")
        cat("Predicted Rt\n")
        print(Rt)
        cat("\n")
      }      
      if(r>tin & disp){
        mTot_new = apply(TOT_new[,,1:round(r/tin)],1:2,mean)
        colnames(mTot_new) = c(categories,"now_pos","tot_pos")
        cat("Predicted totals\n")
        print(head(mTot_new))
        cat("\n")
        cat("Predicted_new Rt\n")
        print(Rt_new)
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
             Y=Y,call = match.call())
  if(output){
    out$TTAB = TTAB; out$BBE = BBE; out$RT = RT; out$TTAB_new = TTAB_new
    out$TOT_new = TOT_new; out$RT_new = RT_new; out$DIFF_new = DIFF_new
  }
  class(out) = "multinomAR"
  return(out)

}