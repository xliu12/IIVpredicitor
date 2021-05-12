\name{TPB.latentIVAR}
\alias{TPB.latentIVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
TPB.latentIVAR
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TPB.latentIVAR(
  compx,
  y,
  nB
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{compx}{an N*T (N: number of individuals; T: number of measurement occasions) matrix containing wide format repeated measures of a time-varying independent variable X}
\item{y}{a column vector of length N containing a single measurement of the outcome variable y
}
\item{nB}{number of bootstrap samples}
}
\details{
The sandwisch type standard error estimation accounts for the clustering in the treatment arm.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
IPW.MLMPN returns a list "ipw10.mlmpn" containing the following components:
%%  If it is a LIST, use list(gammDiff,se_sw, z.wald)
%%  \item{comp1 }{Description of 'comp1'}
\item{gammDiff }{the treatment effect estimate, i.e., the mean difference between the treatment and control arm. }
\item{se_sw }{the sandwich-type standard error estimate of the treatment effect estimate. }
\item{z.wald }{the Wald statistic (i.e., gammDiff/se_sw). }
}
%%  \references{
%% ~put references to the literature/web site here ~ }
%%  \author{
%%  ~~who you are~~
%%  }
%%  \note{
%%  ~~further notes~~ }

%% ~Make other sections like Warning with \section{Warning }{....} ~
%%  \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~}

\examples{
%\dontrun{
\donttest{
%##-- ==>  Define data, use random,
#  Population Values  ----------------
# IM: mux_i ~ N(0, var=9)
sigma_mu=sqrt(9);
# IIV: sigmasqx_i ~ Gamma(shape=0.25, scale=17)
shape=0.25;
scale=17;
# measurement error
sigma_e=sqrt(1.8) # such that the measurement scale reliability is .88
# y_i = beta_mu*IM_i + beta_IIV*ISD_i + epsilon_i
# epsilon_i ~ N(0,1)
beta_mu=1;
sigma_ep=1

# beta_IIV is beta2 for ISD
beta_IIV=c(0,0.3,0.8,1.50)[2]
N=c(200,500)[1]
T=c(6,9,12)[2]

#  Data Generation -----
set.seed(12)
mu<-rnorm(N,0,sigma_mu)
IIV<-rgamma(N,shape=shape,scale=scale) #sigmasq_vi
isd<-sqrt(IIV)
# population regression model
# y<-0+beta_mu*mu+beta_IIV*IIV+rnorm(N,0,sigma_ep)
# use ISD as predictor
y<-0+beta_mu*mu+beta_IIV*isd+rnorm(N,0,sigma_ep)

x=v=compx=array(0,c(N,T))
for (i in 1:N) {
  v[i,]=mu[i]+rnorm(T,0,isd[i])
  x[i,]=v[i,]+rnorm(T,0,sigma_e)
  compx[i,]=x[i,]
}

TPB.latentIVAR(
  compx,
  y,
  nB
)

%}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~design }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~models }% __ONLY ONE__ keyword per line
