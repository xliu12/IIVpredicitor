
library(lavaan)

IVAR.TPB <- function(
  compx, # an N*T (N: number of individuals; T: number of measurement occasions) matrix containing wide format repeated measures of a time-varying independent variable X
  y # a column vector of length N containing a single measurement of the outcome variable y
  , nB=1000 # number of bootstrap samples
  # ,numpar=3 # number of parcels
){
  # Time parceling with bootstrap ---------------------
  numpar=3
  T=ncol(compx)
  indexC=as.numeric(cut(1:T, numpar))

  beta_mu.tpB<-array(0,nB)
  beta_IIV.tpB<-array(0,nB)
  conv.tpB<-array(0,nB)
  for(B in 1:nB){
    indexB=sample(1:N,replace = TRUE)
    compxtpB <- compx[indexB,]
    ytpB=y[indexB]

    mean <- IIV <- NULL
    for(i in 1:numpar){
      x <- compxtpB[, indexC==i]
      IIV[[i]] <- apply(x,1,var)
      mean[[i]] <- apply(x,1,mean)
    }
    IIVx1 <- IIV[[1]]
    IIVx2 <- IIV[[2]]
    IIVx3 <- IIV[[3]]
    meanx1 <- mean[[1]]
    meanx2 <- mean[[2]]
    meanx3 <- mean[[3]]
    fitdat <- data.frame(ytpB,IIVx1,IIVx2,IIVx3,meanx1,meanx2,meanx3)

    model<-'
    mu =~ meanx1+meanx2+meanx3
    IIV =~ IIVx1+IIVx2+IIVx3
    ytpB ~ mu+IIV
    mu ~~ 0*IIV
    meanx1~0
    meanx2~0
    meanx3~0
    IIVx1~0
    IIVx2~0
    IIVx3~0
    mu~1
    IIV~1
    '
    fittpB <- sem(model,data=fitdat,mimic="Mplus")
    # summary(fit)
    est.coef.tpB<- parameterEstimates(fittpB)
    beta_mu.tpB[B] <- est.coef.tpB$est[7]
    beta_IIV.tpB[B] <- est.coef.tpB$est[8]
    conv.tpB[B] <- fittpB@optim$converged

    if(!fittpB@optim$converged){
      beta_mu.tpB[B]=NA
      beta_IIV.tpB[B]=NA
    }
  }

  est.beta_mu.tpB=mean(beta_mu.tpB,na.rm=T)
  se.beta_mu.tpB<- sd(beta_mu.tpB,na.rm=T)
  l.beta_mu.tpB<- quantile(beta_mu.tpB,0.025,na.rm=T)
  u.beta_mu.tpB <- quantile(beta_mu.tpB,0.975,na.rm=T)

  est.beta_IIV.tpB=mean(beta_IIV.tpB,na.rm=T)
  se.beta_IIV.tpB<- sd(beta_IIV.tpB,na.rm=T)
  l.beta_IIV.tpB<- quantile(beta_IIV.tpB,0.025,na.rm=T)
  u.beta_IIV.tpB <- quantile(beta_IIV.tpB,0.975,na.rm=T)

  conv.tpB <- mean(conv.tpB)


  resv=c(
    "IIVtpB"=est.beta_IIV.tpB, "se.IIVtpB"=se.beta_IIV.tpB,
    "l.IIVtpB"=l.beta_IIV.tpB, "u.IIVtpB"=u.beta_IIV.tpB,
    "mutpB"=est.beta_mu.tpB, "se.mutpB"=se.beta_mu.tpB,
    "l.mutpB"=l.beta_mu.tpB, "u.mutpB"=u.beta_mu.tpB
  )
  res = data.frame(matrix(resv, nrow = 2, byrow = T))
  colnames(res) = c('Estimate','Std.err', 'CI.lower', 'CI.upper')
  rownames(res)=c('beta_IIV','beta_IM')

  reslist=list(res=res, convergence=conv.tpB)

  return(reslist)
}





















