\name{IVAR.TPB}
\alias{IVAR.TPB}
%  Also NEED an '\alias' for EACH other topic documented here.
\title{
Time parceling (with bootstrap) method for modeling intraindividual variance as a predictor
%%  ::function to do ... ::
}
\description{
Use observed intraindividual means and variances of multiple (e.g., 3) time parcels as manifest indicators of true (latent) intraindividual mean and variance respectively to predict an outcome variable in a structural equation model (SEM).
}
\usage{
IVAR.TPB(
  compx,
  y,
  nB=1000
)
}
%  maybe also 'usage' for other objects documented here.
\arguments{
\item{compx}{an N*T (N: number of individuals; T: number of measurement occasions) matrix containing wide format repeated measures of a time varying independent variable X.}
\item{y}{a column vector of length N containing the outcome variable Y.
}
\item{nB}{number of bootstrap samples (default: 1000).}
}
\details{
A fixed number (numpar = 3) of time parcels are formed for each individual using the individual's repeated measures of the independent variable X at T occasions. The time parcels are created by cutting the sequence 1,2,...,T into three intervals of (approximately) equal length (about T/numpar).
Then for each individual, sample intraindividual means and variances of the created time parcels are computed and used as manifest indicators of the individual's true (latent) intraindividual mean and variance respectively to predict the individual's outcome score in a structural equation model (SEM).
Sampling (i.e., sampling individuals from population) variablility of the parameter estimates is accounted for using the Bootstrap approach. #nB bootstrap samples are drawn to estimate the standard errors and 0.95 confidence intervals.
%%  :: If necessary, more details than the description above ::
}
\value{
%%  :Describe the value returned
IVAR.TPB returns a list "reslist" containing two components.
%%  If it is a LIST, use list(gammDiff,se_sw, z.wald)
\item{res }{a data.frame named "res" containing the estimates ("Estimate"), standard error estimates ("Std.err"), and lower and upper bounds of 0.95 confidence intervals of
the effects of the intraindividual variance ("beta_IIV") and intraindividual mean ("beta_IM") of X (the time varying indenpendent variable) on Y (the outcome variable).}
\item{convergence }{convergence rate of the bootstrap samples.}
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
# Don't run
#  Population Values
# IM: mux_i  N(0, var=9)
sigma_mu=sqrt(9);
# IIV: sigmasqx_i  Gamma(shape=0.25, scale=17)
shape=0.25;
scale=17;
# measurement error
sigma_e=sqrt(1.8) # such that the measurement scale reliability is .88
# y_i = beta_mu*IM_i + beta_IIV*ISD_i + epsilon_i
# epsilon_i  N(0,1)
beta_mu=1;
sigma_ep=1

# beta_IIV is beta2 for ISD
beta_IIV=c(0,0.3,0.8,1.50)[2]
N=c(200,500)[1]
T=c(6,9,12)[2]

#  Data Generation
set.seed(12)
mu=rnorm(N,0,sigma_mu)
IIV=rgamma(N,shape=shape,scale=scale) #sigmasq_vi
isd=sqrt(IIV)
# population regression model
# y=0+beta_mu*mu+beta_IIV*IIV+rnorm(N,0,sigma_ep)
# use ISD as predictor
y=0+beta_mu*mu+beta_IIV*isd+rnorm(N,0,sigma_ep)

x=v=compx=array(0,c(N,T))
for (i in 1:N) {
  v[i,]=mu[i]+rnorm(T,0,isd[i])
  x[i,]=v[i,]+rnorm(T,0,sigma_e)
  compx[i,]=x[i,]
}

IVAR.TPB(
  compx,
  y,
  nB=10
)

%}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ :design }% use one of  RShowDoc("KEYWORDS")
%\keyword{ :models }% __ONLY ONE__ keyword per line
