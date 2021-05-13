\name{ISD.BV}
\alias{ISD.BV}
%  Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian variability modeling method for modeling intraindividual standard deviation as a predictor
%%  ::function to do ... ::
}
\description{
Intraindividual mean and standard deviation are modeled as predictors of an outcome variable using Bayesian approach, where hyper priors are specified for the true (i.e., latent) intraindividual mean and standard deviation.
%%  :: A concise (1 5 lines) description of what the function does. ::
}
\usage{
ISD.BV(
  compx,
  y,
  ReMeasure
  , nchains=2, niter=8000
)
}
%  maybe also 'usage' for other objects documented here.
\arguments{
\item{compx}{an N*T (N: number of individuals; T: number of measurement occasions) matrix containing wide format repeated measures of a time varying independent variable X}
\item{y}{a column vector of length N containing a single measurement of the outcome variable y
}
\item{ReMeasure}{measurement scale reliability of the time varying independent variable X}
\item{nchains}{number of Markov chains (default: 2)}
\item{niter}{number of total iterations per chain (including burn in; default: 8000)}
}
\details{
For each individual, the predictor score at occasion t is assumed to follow a normal distribution:

\emph{x_it : N( mux_i, (sigmav_i)^2 + sigmaex^2 ),}

where the intraindividual variance (sigmav_i)^2 is modeled as \emph{sigmav_i^2 : Gamma(location, scale)} with unknown location and scale parameters; and the intraindividual mean (mux_i) as \emph{mux_i : N(mu_mux, sigma_mux^2)}. The measurement scale reliability of X ("ReMeasure") is used to help identify the model, where we set \emph{sigmaex^2 = (1--
ReMeasure )*(sigma_mux^2+location*scale)/ReMeasure}.

For predicting outcome Y, we assume
\emph{y_i : N( beta_0+beta_IM*mux_i+beta_IIV*sigmav_i, sigmay^2 )}.

Half Cauchy priors are used as a weakly informative prior for location, 1/scale, sigma_mux, and sigmay. Noninformative normal priors,N(0,0.001^2) are used for the regression coefficients. R JAGS is used to implement the MCMC sampling.

}
\value{
%%  :Describe the value returned
IVAR.BV returns a data.frame named "res" containing the mean("Estimate"), standard deviation ("Std.err"), 0.025 and 0.975 quantiles of the posterior distributions of
the effects of the intraindividual standard deviation ("beta_IIV") and intraindividual mean ("beta_IM") of X (the time varying indenpendent variable) on Y (the outcome variable),
, and their corresponding Rhat (the potential scale reduction factor of the Gelman Rubin convergence diagnositic) values.
%%  If it is a LIST, use list(gammDiff,se_sw, z.wald)
%%  \item{comp1 }{Description of 'comp1'}
}
%%  \references{
%% :put references to the literature/web site here : }
%%  \author{
%%  ::who you are::
%%  }
%%  \note{
%%  ::further notes:: }

%% :Make other sections like Warning with \section{Warning }{....} :
%%  \seealso{
%% ::objects to See Also as \code{\link{help}}, :::}

\examples{%\dontrun{
\donttest{
%##   ==>  Define data, use random,
# Don't run
#  Population Values
# IM: mux_i : N(0, var=9)
sigma_mu=sqrt(9);
# IIV: sigmasqx_i : Gamma(shape=0.25, scale=17)
shape=0.25;
scale=17;
# measurement error
sigma_e=sqrt(1.8) # such that the measurement scale reliability is .88
# y_i = beta_mu*IM_i + beta_IIV*ISD_i + epsilon_i
# epsilon_i : N(0,1)
beta_mu=1;
sigma_ep=1

# beta_IIV is beta2 for ISD
beta_IIV=c(0,0.3,0.8,1.50)[2]
N=c(200,500)[1]
T=c(6,9,12)[2]

#  Data Generation
set.seed(12)
mu< rnorm(N,0,sigma_mu)
IIV< rgamma(N,shape=shape,scale=scale) #sigmasq_vi
isd< sqrt(IIV)
# population regression model
# y< 0+beta_mu*mu+beta_IIV*IIV+rnorm(N,0,sigma_ep)
# use ISD as predictor
y< 0+beta_mu*mu+beta_IIV*isd+rnorm(N,0,sigma_ep)

x=v=compx=array(0,c(N,T))
for (i in 1:N) {
  v[i,]=mu[i]+rnorm(T,0,isd[i])
  x[i,]=v[i,]+rnorm(T,0,sigma_e)
  compx[i,]=x[i,]
}


ISD.BV(
  compx,
  y,
  ReMeasure=0.88
  , nchains=2, niter=8000
)

%}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ :design }% use one of  RShowDoc("KEYWORDS")
%\keyword{ :models }% __ONLY ONE__ keyword per line
