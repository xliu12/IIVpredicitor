
library(R2jags)



IVAR.BV = function(
  compx, # an N*T (N: number of individuals; T: number of measurement occasions) matrix containing wide format repeated measures of a time-varying independent variable X
  y, # a column vector of length N containing a single measurement of the outcome variable y
  ReMeasure=0.88 # scale reliability
  , nchains=2, niter=8000
){

  compxl=matrix(c(t(compx)),ncol = 1)
  yl=matrix(rep(y,each=T),ncol = 1)
  IDl=matrix(rep(1:N,each=T),ncol = 1)
  datal=cbind(yl,compxl,IDl)
  colnames(datal)=c("Y","X","ID")

  # jags data
  jagsdat=list(
    N=nrow(compx)
    ,J=ncol(compx)
    ,X=compx
    ,y=y
    ,ReMeasure=ReMeasure
  )

  # jags model
  write("
model {

            for(i in 1:N){
            for(j in 1:J){
            #x_it = mu_i + v_it + e_it
            X[i,j] ~ dnorm(mux[i],taux[i])
            }
            mux[i] ~ dnorm(mu_mux,tau_mux)

            sigmasq_v[i] ~ dgamma(shape_v,rate_v)
            sigma_v[i] <- sqrt(sigmasq_v[i])
            taux[i] <- 1/(sigmasq_v[i]+sigmasq_e)

            y[i] ~ dnorm(muy[i],tauy)
            muy[i] <- beta0+beta1*mux[i]+beta2*sigmasq_v[i]
            }
            mu_mux ~ dnorm(0,1e-6)
            sigma_mux ~ dt(0, pow(100,-2), 1) T(0,)
            sigmasq_mux <- sigma_mux^2
            tau_mux <- 1/sigmasq_mux

            shape_v ~ dt(0, pow(100,-2), 1) T(0,)
            #scale_v ~ dt(0, pow(100,-2), 1) T(0,)
            #rate_v <- 1/scale_v
            rate_v ~ dt(0, pow(100,-2), 1) T(0,)
            scale_v <- 1/rate_v

            sigmasq_e <- (1-ReMeasure)/ReMeasure * (1/tau_mux+shape_v/rate_v)

            beta0 ~ dnorm(0,1e-6)
            beta1 ~ dnorm(0,1e-6)
            beta2 ~ dnorm(0,1e-6)
            sigma_y ~ dt(0, pow(100,-2), 1) T(0,)
            sigmasq_y <- sigma_y^2
            tauy <- 1/sigmasq_y

            }
        ","jagsmodel_ISD2.txt")
  # run jags
  jagsout=jags(data = jagsdat
               , parameters.to.save = c( "beta1", "beta2","beta0", "shape_v", "scale_v", "sigmasq_y", "sigmasq_mux")
               # , inits = inits_tr
               ,model.file = "jagsmodel_ISD2.txt",n.chains = nchains,
               n.iter = niter)
  jagsres=jagsout$BUGSoutput$summary

  est.beta_mu.BV <- jagsres[1,1]
  est.beta_IIV.BV <- jagsres[2,1]

  se.beta_mu.BV <- jagsres[1,2]
  se.beta_IIV.BV <- jagsres[2,2]

  l.beta_IIV.BV <- jagsres[2,3]
  u.beta_IIV.BV <- jagsres[2,7]
  l.beta_mu.BV <- jagsres[1,3]
  u.beta_mu.BV <- jagsres[1,7]

  Rhat.beta_mu <- jagsres[1,c("Rhat")]
  Rhat.beta_IIV <- jagsres[2,c("Rhat")]

  resv=c(
    IIVBV=est.beta_IIV.BV, se.IIVBV=se.beta_IIV.BV,
    l.IIVBV=l.beta_IIV.BV, u.IIVBV=u.beta_IIV.BV,
    muBV=est.beta_mu.BV, se.muBV=se.beta_mu.BV,
    l.muBV=l.beta_mu.BV, u.muBV=u.beta_mu.BV,
    Rhat.beta_mu=Rhat.beta_mu, Rhat.beta_IIV=Rhat.beta_IIV )

  res = data.frame(matrix(resv, nrow = 2, byrow = T))
  colnames(res) = c('Estimate','Std.err', '2.5% quantile', '97.5% quantile', 'Rhat')
  rownames(res)=c('beta_IIV','beta_IM')

  return(res)
}


