\name{HDBRR}
\alias{HDBRR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High Dimensional Bayesian Ridge Regression without MCMC.
}
\description{
Ridge regression provide biased estimators of the regression parameters with lower variance. The HDBRR ("High Dimensional Bayesian Ridge Regression") function fits Bayesian Ridge regression without MCMC, this one uses the SVD or QR decomposition for the posterior computation.
}
\usage{
HDBRR(y, X, a1 = 1, a2 = 1, b1 = 1, b2 = 1/ncol(X),
    intercept = TRUE, npts = NULL, c = NULL, corpred = NULL,
    method = c("svd", "qr"),bigmat = TRUE, ncores = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
The data vector (numeric, n) NAs allowed.
}
  \item{X}{
Design Matrix of dimension \code{n x p}.
}
  \item{a1}{
\code{n} a priori. \code{a1/2} is the shape parameter of the Gamma Inverse prior assigned to the residual variance. The default value for the \code{a1} parameter is 1.
}
  \item{a2}{
\code{p} a priori. \code{a2/2} is the shape parameter of the Gamma Inverse prior assigned to the Beta's variance. The default value for the \code{a2} parameter is 1.
}
  \item{b1}{
\code{b1/2} is the scale parameter of the Gamma Inverse prior assigned to the residual variance. The default value for the \code{b1} parameter is 1.
}
  \item{b2}{
\code{b2/2} is the scale parameter of the Gamma Inverse prior assigned to the Beta's variance. The default value for the \code{b2} parameter is \code{1/ncol(X)}, since the beta's variance is small, in this case \code{ncol(X) = p}, if \code{p} is large then \code{b2} is small.
}
  \item{intercept}{
Logic value. The default value for the intercept is TRUE.
}
  \item{npts}{
Number (integer) of points used to evaluate the u's density for the numeric aprroach. The default value for the npts parameter is 200.
}
  \item{c}{
Auxiliary constant for the prior density for Beta in the variable selection, in this case the prior density is a mixture of a normal distribution.
}
  \item{corpred}{
The method for the compute of the correlation, there are two methods, Empirical Basis (\code{"eb"}) and Bayesian (\code{"b"}) method. The default value for the parameter corpred is NULL. If the values is NULL then the corr and edf values will be NULL.
}
  \item{method}{
Options for the posterior computation. There are two methods available: \code{"qr"} and \code{"svd"} decomposition. The default value for the method is SVD decomposition.
}
 \item{bigmat}{
Use of the bigstatsr package. The default value for bigmat is \code{TRUE}.
}
  \item{ncores}{
Number of the cores for the auxiliary functions. The default value for the ncores is 2, you can detect your number of cores with detectCores() and use it.
}
}
\details{
Ridge regression is a useful tool to deal with colinerity in the homocesastic linear regression model providing biased estimators of the regression parameters with lower variance than the least square estimators.
The model \deqn{y = X\beta + \epsilon} where \eqn{\epsilon} vector is assumed Normal with mean vector 0 and covariance matrix \eqn{\sigma^2 I_n}. See more \url{https://drive.google.com/drive/folders/1xJw7gM5_XiJipQ3grTZkfc6q4K0hzuCx?usp=sharing}
}
\value{
List containing the following components:
\item{betahat}{
Vector (numeric, \code{p})  with the betas estimates.
}
\item{yhat}{
Vector (numeric, \code{n}) with the y's estimates.
}
\item{varb}{
Vector (numeric, \code{p}) with the beta's variance.
}
\item{sigsqhat}{
Value (numeric) of the residual variance estimate.
}
\item{sigbsqhat}{
Value (numeric) of the Beta's variance estimate.
}
\item{u}{
Vector (numeric, \code{npts}) with the u's values.
}
\item{postu}{
Vector (numeric, \code{npts}) with the values of the u posterior.
}
\item{uhat}{
Value (numeric) of u estimated.
}
\item{umode}{
Value (numeric) of the posterior mode of u.
}
\item{whichNa}{
Value (integer) of NAs in the y vector.
}
\item{phat}{
Vector (numeric, \code{p}), selection probability of x_i.
}
\item{delta}{
Used in the variable selection.
}
\item{edf}{
Value (numeric) of the effective degrees of freedom for regression.
}
\item{corr}{
Vector (numeric, \code{n}) of the correlation between \code{y_i} estimates and \code{y_i}.
}
}

\author{
Sergio Perez-Elizalde, Blanca E. Monroy-Castillo, Paulino Perez-Rodriguez.
}

\examples{
n <- 30
p <- 100
X <- matrix(rnorm(n*(p-1),1,1/p),nrow = n,ncol = p-1)
Beta <- sample(1:p,p-1,rep = FALSE)
Beta <- c(1,Beta)
y <- cbind(rep(1,n),X) \%*\% Beta+rnorm(n,0,1)
fit <- HDBRR(y,X/sqrt(p), intercept = TRUE, c = 100, cor = "b")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BRR.}% use one of  RShowDoc("KEYWORDS")
\keyword{ BRR without MCMC.}% __ONLY ONE__ keyword per line
