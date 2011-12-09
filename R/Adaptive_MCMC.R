## =======================================================
## (Adaptive) Metropolis Sampler
##
## Implementation of the RAM (robust adaptive Metropolis)
## sampler of
## Vihola, M. (2011) Robust adaptive Metropolis algorithm with
## coerced acceptance rate. Statistics and Computing.
## [online] http://www.springerlink.com/content/672270222w79h431/
## (Accessed December 8, 2011).

##
## December 9, 2011 -- Andreas Scheidegger
## =======================================================


MCMC <- function(p, n, init, scale=rep(1, length(init)), log=TRUE,
                 adapt=!is.null(acc.rate), acc.rate=NULL, gamma=0.55, list=TRUE, verbose=100, n.start=0, ...) {

  ## checks
  if(adapt & !is.numeric(acc.rate)) stop('Argument "acc.rate" is missing!')
  if(gamma<0.5 | gamma>1) stop('Argument "gamma" must be between 0.5 and 1!')

  ## number of adaption steps
  if(is.numeric(adapt)) n.adapt <- adapt
  if(adapt==TRUE) n.adapt <- Inf
  if(adapt==FALSE) n.adapt <- 0

  ## number of parameter
  d <- length(init)

  ## matrix to store MC chain
  X <- matrix(NA, ncol=d, nrow=n)
  colnames(X) <- names(init)
  X[1,] <- init

  ## vector to store densities p(x)
  p.val <- rep(NA, n)
  p.val[1] <- p(X[1,], ...)

  ## initial S
  if(length(scale)==d) {
    M <- diag(scale)
  } else {
    M <- scale
  }
  ## check
  if(ncol(M) != length(init)) stop("Length or dimension of 'init' and 'scale' do not match!")
  S <-  t(chol(M))

  k <- 0
  for(i in 2:n){
    if(i %% verbose == 0) {
      cat('iteration', i, 'of', n, '\n')
    }

    ## proposal value
    U <- rnorm(d)
    X.prop <- c( X[i-1,] + S %*% U )

    ## calculate density at X.prop
    p.val.prop <- p(X.prop, ...)

    ## acceptance probability
    if(log) {
      alpha <- min(1, exp( p.val.prop - p.val[i-1] )) # for log density
    } else {
      alpha <- min(1, p.val.prop/p.val[i-1] )
    }
    if(!is.finite(alpha)) alpha <- 0    # if zero divided by zero

    ## accept with P=alpha
    if(runif(1)<alpha) {
      X[i,] <- X.prop                   # accept
      p.val[i] <- p.val.prop
      k <- k+1
    } else {
      X[i,] <- X[i-1,]                  # or not
      p.val[i] <- p.val[i-1]
    }

    ## compute new S
    ii <- i+n.start
    if(ii < n.adapt) {
      adapt.rate <-  min(1, d*ii^(-gamma))
      M <- S %*% (diag(d) + adapt.rate*(alpha - acc.rate) * U%*%t(U)/sqrt(sum(U^2))) %*% t(S)

      ## check if M is positive definite. If not, use nearPD().
      eig <- eigen(M, only.values = TRUE)$values
      tol <- ncol(M)*max(abs(eig))*.Machine$double.eps

      if( !isSymmetric(M) | is.complex(eig) | !all(Re(eig)>tol) ){
        ## nearPD() computes the 'nearest' positive definite matrix
        M <- as.matrix(Matrix::nearPD(M)$mat)
      }

      S <- t(chol(M))

    }
  }

  ## calculate accpetance rate
  acceptance.rate <- round(k/(n-1), 3)

  if(list) {
      return(list(samples=X,
                  cov.jump=M,
                  n.sample=n,
                  acceptance.rate=acceptance.rate,
                  adaption=adapt,
                  sampling.parameters=list(sample.density=p,
                                          log=log,
                                          acc.rate=acc.rate,
                                          gamma=gamma)
                  ) )
  } else {
    cat("Acceptance rate:", acceptance.rate, "\n")
    return(X)
  }
}



## ----------------------
## add more samples to an existing chain
## does work with objet generated with MCMC.parallel()
## but updating is not computed  parallel.

MCMC.add.samples <- function(MCMC.object, n.update, verbose=100, ...) {

    ## if single chain
    if(!is.null(names(MCMC.object))) {
        if(is.matrix(MCMC.object)) stop("Only MCMC objects generated with option 'list=TRUE' can be updated!")

        ## get values from last sampling
        p <- MCMC.object$sampling.parameters$sample.density
        init <- MCMC.object$samples[nrow(MCMC.object$samples),]
        scale <- MCMC.object$cov.jump
        log <- MCMC.object$sampling.parameters$log
        acc.rate <- MCMC.object$sampling.parameters$acc.rate
        gamma <- MCMC.object$sampling.parameters$gamma
        n.before <- MCMC.object$n.sample    # number of existing samples
        adapt <- MCMC.object$adaption

        ## generate new samples
        samp.update <- MCMC(p=p, n=n.update, init=init, scale=scale, log=log,  adapt=adapt, acc.rate=acc.rate,
                            gamma=gamma, list=TRUE, verbose=verbose, n.start=n.before, ...)

        ## update old sampling object
        MCMC.object$cov.jump <- samp.update$cov.jump
        m <- c(MCMC.object$n.sample, samp.update$n.sample)
        MCMC.object$acceptance.rate <-  1/sum(m)*(m[1]*MCMC.object$acceptance.rate + m[2]*samp.update$acceptance.rate)
        MCMC.object$n.sample <- MCMC.object$n.sample + n.update

        MCMC.object$samples <- rbind(MCMC.object$samples, samp.update$samples)

        ## return the updated list
        return(MCMC.object)
    }

    ## if list of chains
    if(is.null(names(MCMC.object))) {
        ## recursive call of MCMC.add.samples() to update single chains
        MCMC.object <- lapply(MCMC.object, function(x) MCMC.add.samples(x, n.update=n.update, ...))
        return(MCMC.object)
    }
}


## ----------------------
## parallel version
MCMC.parallel <- function(p, n, init, n.chain=4, n.cpu=2, packages=NULL, scale=rep(1, length(init)), log=TRUE,
                 adapt=!is.null(acc.rate), acc.rate=NULL, gamma=0.55, list=FALSE, parallel=TRUE, verbose=100, ...) {

  require(snowfall)

  ## initialisation of (local) parallel computing
  sfInit(parallel=parallel, cpus=min(n.chain, n.cpu), type="SOCK")

  ## stop parallel computing on exit
  on.exit(sfStop())

  ## export complete work space of master
  sfExportAll()

  ## exports library
  sfLibrary("Matrix", character.only=TRUE)

  ## load additional packages
  if(!is.null(packages)) {
      sapply(packages, function(x) sfLibrary(x, character.only=TRUE) )
  }

  ## init random generators
  sfClusterSetupRNG()

  ## wrapper
  MCMC.wrap <- function(x, ...) {
      MCMC(...)
  }

  ## sample n.chin in parallel
  result <- sfLapply(1:n.chain, MCMC.wrap, p=p, n=n, init=init,
                        scale=scale, log=log, adapt=adapt, acc.rate=acc.rate,
                        gamma=gamma, list=list, verbose=verbose, ...)


  return(result)
}



## ----------------------
## converts a sample into coda object
convert.to.coda <- function(sample) {
    ## if single chain
    if(!is.null(names(sample))) {
        if(is.matrix(sample)) {
            obj <- coda::mcmc(sample)
        }
        if(is.list(sample)) {
            obj <- coda::mcmc(sample$samples)
        }
        return(obj)
    } else {

        ## if sample is a list of chains
        if(is.matrix(sample[[1]])) {
            obj <- as.mcmc.list(lapply(sample, coda::mcmc))
        }
        if(is.list(sample[[1]])) {
            obj <- as.mcmc.list(lapply(sample, function(x) {coda::mcmc(x$samples)}))
        }
        return(obj)
    }
}

