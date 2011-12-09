\name{MCMC.parallel}
\alias{MCMC.parallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parallel computation of MCMC()
}
\description{
Several independent Markov chains are generated in parallel by using the
\pkg{snowfall} package.
}
\usage{
MCMC.parallel(p, n, init, n.chain = 4, n.cpu = 2, packages = NULL,
    scale = rep(1, length(init)), log = TRUE, adapt = !is.null(acc.rate),
    acc.rate = NULL, gamma = 0.55, list = FALSE, parallel = TRUE,
    verbose = 100, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
function that returns the (log) probability density to sample from. Must have two or more dimensions.
}
  \item{n}{
number of samples.
}
  \item{init}{
vector with initial values.
}
  \item{n.chain}{
number of independent chains.
}
  \item{n.cpu}{
number of CPUs that can be used in parallel
}
  \item{packages}{
    vector with name of packages to load into each instance. (Typically,
    all packages on which \code{p} depends.)
}
  \item{scale}{
vector with the variances \emph{or} covariance matrix of the jump distribution.
}
  \item{log}{
logical. If \code{TRUE}, a log density is expected from \code{p} (strongly recommended).
}
  \item{adapt}{
 if \code{TRUE}, adaptive sampling is used, if \code{FALSE} classic metropolis
 sampling, if a positive integer the adaption stops after \code{adapt} iterations.
}
  \item{acc.rate}{
desired acceptance rate (ignored if \code{adapt=FALSE})
}
  \item{gamma}{
 controls the speed of adaption. Should be between 0.5 and 1. A lower
 gamma leads to faster adaption.
}
  \item{list}{
logical. If \code{TRUE} a list of lits is returned otherwise a list of matrices with the samples.
}
\item{parallel}{
  logical. Should the chains be calculated in parallel?
}
  \item{verbose}{
every \code{verbose}th iteration is printed.
}
  \item{\dots}{
    further arguments passed to \code{p}
}

}

\details{
This function is just a wrapper to use \code{MCMC} in parallel. It is
based on \pkg{snowfall}.
}
\value{
A list with a list or matrix for each chain. See \code{\link{MCMC}} for details.
}

\author{Andreas Scheidegger, \email{andreas.scheidegger@eawag.ch} or
  \email{scheidegger.a@gmail.com}
}

\seealso{
\code{\link{MCMC}}
}
\examples{

## ----------------------
## Banana shaped distribution

## log-pdf to sample from
p.log <- function(x) {
  B <- 0.03                              # controls 'bananacity'
  -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2
}

## ----------------------
## generate samples
## compute 4 independent chains on 2 CPU's in parallel

samp <- MCMC.parallel(p.log, n=200, init=c(x1=0, x2=1),
    n.chain=4, n.cpu=2, scale=c(1, 0.1),
    adapt=TRUE, acc.rate=0.234)

str(samp)

}


