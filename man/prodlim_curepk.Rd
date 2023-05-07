\name{prodlim_curepk
}
\alias{prodlim_curepk
}
\title{Compute Product-Limit Estimator of Conditional Survival Function when
    Cure Status is Partially Known
}
\description{This function computes the nonparametric estimator of the
    conditional survival function when cure status is partially known
    proposed by Safari \emph{et al} (2021).
}
\usage{
    prodlim_curepk(x, t, d, xinu, dataset, x0, h, local = TRUE,
                   bootpars = if (!missing(h)) NULL else controlpars())
}
\arguments{
    \item{x}{If \code{dataset} is missing, a numeric object giving the
        covariate values. If \code{dataset} is a data frame, it is interpreted
        as the name of the variable corresponding to the covariate in the data
        frame.}
    \item{t}{If \code{dataset} is missing, a numeric object giving the observed
        times. If \code{dataset} is a data frame, it is interpreted as the name
        of the variable corresponding to the observed times in the data frame.}
    \item{d}{If \code{dataset} is missing, an integer object giving the values
        of the uncensoring indicator. Censored observations must be coded as 0,
        uncensored ones as 1. If dataset is a data frame, it is interpreted as
        the name of the variable corresponding to the uncensoring indicator in
        the data frame.}
    \item{xinu}{If \code{dataset} is missing, an integer object giving the
        values of the cure status indicator. Uncensored and unknown censored
        observations must be coded as 0, known to be cured censored ones as 1.
        If \code{dataset} is a data frame, it is interpreted as the name of the
        variable corresponding to the cure status indicator in the data frame.}
    \item{dataset}{An optional data frame in which the variables named in
        \code{x, t, d} and \code{xinu} are interpreted. If it is missing,
        \code{x, t, d} and \code{xinu} must be objects of the workspace.}
    \item{x0}{A numeric vector of covariate values where the estimates of the
        conditional survival function will be computed.}
    \item{h}{A numeric vector of bandwidths.}
    \item{local}{A logical value, TRUE by default, specifying whether
        \code{local} or \code{global} bandwidths are used.}
    \item{bootpars}{A list of parameters controlling the bootstrap when
        computing the bootstrap bandwidths of the product-limit estimator.
        \code{B}, the number of bootstrap resamples, and \code{nnfrac},
        the fraction of the sample size that determines the order of the nearest
        neighbor used for choosing a pilot bandwidth. If \code{h} is missing
        the list of parameters is extended to be the same used for computing the
        bootstrap bandwidth. The default is the value returned by the
        \code{controlpars} function called without arguments.}
}
\value{A list of components:
    \item{h}{The numeric vector of bandwidths used in the estimation. If
        \code{h} argument is missing, the bootstrap bandwidth computed with the
        control parameters in argument \code{bootpars}.}
    \item{x0}{The numeric vector of covariate values where the estimate of the
        conditional survival function is computed.}
    \item{t}{The observed time values, where the conditional survival function
        is estimated.}
    \item{surv}{Estimates of the survival function for each one of the covariate
        values specified by the \code{x0} argument and the bandwidths in
        \code{h}. It is a matrix of dimension \eqn{n\times length(x0)}{(n, length(x0))}
        if local bandwidths or bootstrap bandwidths are used, or an array of
        dimension \eqn{n\times length(x0)\times length(h)}{(n, length(x0), length(h))}
        for global bandwidths instead.}
}
\details{Mixture cure model writes the conditional survival function
    \eqn{S(t\mid x)=P(Y>t\mid X=x)}{S(t|x)=P(Y>t|X=x)} as
    \eqn{S(t\mid x)=1-p(x)+p(x)S_0(t\mid x)}{S(t|x)=1-p(x)+p(x)S_0(t|x)}. 
    This function computes the generalized product-limit estimator of the
    conditional survival function \eqn{S(t \mid x)}{S(t|x)}, using the
    Nadaraya-Watson weights, when the cure status is partially known,
    introduced in Safari \emph{et al} (2021). If there are not individuals
    known to be cured (\code{xinu=0}), then the usual generalized product-limit
    estimator in Beran (1981) is computed.
    
    The Epanechnikov kernel is used. If the smoothing parameter \code{h} is not
    provided, then the bootstrap bandwidth selector in Safari \emph{et al}
    (2021) is used. The function is available only for one continuous covariate
    \eqn{X}.
}
\examples{
library(npcurePK)

## Data-generating function
## n: sample size
## x_cov_range: range of covariate values
## p_knowncure: probability of known cure
data_gen <- function(n, x_cov_range, p_knowncure) {
  ## probability of being susceptible
  p0 <- function(x) exp(2*x)/(1 + exp(2*x))
  ## covariate values
  x <- runif(n, x_cov_range[1], x_cov_range[2])
  ## censoring times
  c <- rexp(n)
  u <- runif(n)
  v <- runif(n)
  data <- data.frame(matrix(0, nrow = n, ncol = 4L,
                            dimnames = list(NULL, c("x", "t", "d", "xinu"))))
  data[, "x"] <- x
  for (i in 1:n) {
    if (u[i] > p0(x[i])) {
      ## Cured individuals (all of them are censored: Yi = infty,
      ## Ti = Ci, delta = 0, nu = 1)
      data[i, "t"] <- c[i]
      if (v[i] < p_knowncure)
        data[i, "xinu"]  <- 1 
    } else {
      ## Uncured individual (Yi < infty, Ti = min(Yi, Ci),
      ## delta = 1(Yi < Ci), nu = 0)
      ## Uncensored individual (d = 1): cure status is
      ## observed (xi = 1), i.e., xinu = 0
      ## Censored individual (d = 0): cure status is
      ## unknown (xi = 0), i.e., xi.nu = 0
      y <- rweibull(1, shape = 0.5 * (x[i] + 4))
      data[i, "t"]  <- ifelse(v[i] < p_knowncure, y, min(y, c[i]))
      if (data[i, "t"] == y) data[i, "d"] <- 1
    }
  }
  return(data)
}

set.seed(123)
data <- data_gen(n = 100, x_cov_range = c(-2, 2), p_knowncure = 0.8)

## Covariate values where the conditional survival function is estimated
x0 <- c(0, 0.5)

## Survival estimates for covariate values x0 = c(0, 0.5)
## ... (a) with 3 global bandwidths (0.5, 1, 1.25)
## The survival function S(t|x) is estimated for each value of x0 with the three
## bandwidths (local == FALSE).
## The estimates are saved in an array (n x length(x0) x length(h))
S1 <- prodlim_curepk(x, t, d, xinu, data, x0 = x0, h = c(0.5, 1, 1.25), local = FALSE)

## Plot predicted survival curve for covariate value x0 = 0.5 and bandwidth
## h = 0.5
x0 <- 0.5
plot(S1$t, S1$surv[, 2, 1], type = "s", xlab = "Time",
     ylab = "Survival probability", ylim = c(0, 1))
## The true survival curve is included as reference
lines(S1$t, 1 - exp(2*x0)/(1 + exp(2*x0)) + exp(2*x0)/(1 + exp(2*x0))*
            (1 - pweibull(S1$t, shape = 0.5 * (x0 + 4))), lwd = 2)     
       
## Plot predicted survival curve for covariate value x0 = 0.5 and all
## bandwidths 
plot(S1$t, S1$surv[, 2, 1], type = "s", xlab = "Time",
     ylab = "Survival probability", ylim = c(0, 1))
lines(S1$t, S1$surv[, 2, 2], type = "s", lwd = 2)
lines(S1$t, S1$surv[, 2, 3],  type = "s", lwd = 3)
# The true survival curve is included as reference
lines(S1$t, 1 - exp(2*x0)/(1 + exp(2*x0)) + exp(2*x0)/(1 + exp(2*x0))*
            (1 - pweibull(S1$t, shape = 0.5 * (x0 + 4))), lwd = 2)     
     
## ... (b) with local bandwidths h = (3, 1)
## The survival function S(t|x) is estimated for each value of x0 with the
## corresponding assigned bandwidth (local == TRUE).
## Note that the length of the vector x0 and the bandwidth h must be the same.
## The estimates are saved in a matrix of dimension (n, length(x0))
x0 <- c(0, 0.5)
h  <- c(3, 1)
S3 <- prodlim_curepk(x, t, d, xinu, data, x0 = x0, h = h, local = TRUE)
## Plot predicted survival curve for covariate value x = 0 and its bandwidth
## (h = 3)
plot(S3$t, S3$surv[, 1], type = "s", xlab = "Time",
     ylab = "Survival probability", ylim = c(0, 1))
## The true survival curve is included as reference     
x0 <- 0 
lines(S3$t, 1 - exp(2*x0)/(1 + exp(2*x0)) + exp(2*x0)/(1 + exp(2*x0))*
            (1 - pweibull(S3$t, shape = 0.5 * (x0 + 4))), lwd = 2)
\donttest{
## ... (c) with the bootstrap bandwidth selector in x0 = 0 (the default
## when the bandwidth argument h is not provided).
## The bootstrap bandwidth is searched in a grid of 10 bandwidths (hl = 10)
## between 0.2 and 2 times the standardized interquartile range of the
## covariate values (hbound = c(0.1, 2)).
x0 <- 0
(S4 <- prodlim_curepk(x, t, d, xinu, data, x0 = x0))
## Equivalently
(S4 <- prodlim_curepk(x, t, d, xinu, data, x0 = x0,
                      bootpars = controlpars(hl = 10, hbound = c(0.1, 2))))
## Plot predicted survival curve for covariate value x = 0 and the bootstrap
## bandwidth
plot(S4$t, S4$surv[, 1], type = "s", xlab = "Time",
     ylab = "Survival probability", ylim = c(0, 1))
## The true survival curve is included as reference     
lines(S4$t, 1 - exp(2*x0)/(1 + exp(2*x0)) + exp(2*x0)/(1 + exp(2*x0))*
            (1 - pweibull(S4$t, shape = 0.5 * (x0 + 4))), lwd = 2)

## ... (d) with parallel computation (The bootstrap bandwidth is searched with
## b = 100 bootstrap resamples and 2 cores)
library(doParallel)
(S5 <- prodlim_curepk(x, t, d, xinu, data, x0 = x0,
                      bootpars = controlpars(b = 100, ncores = 2)))
## Plot predicted survival curve for covariate value x = 0 and the bootstrap
## bandwidth
plot(S5$t, S5$surv[, 1], type = "s", xlab = "Time",
     ylab = "Survival probability", ylim = c(0, 1))
## The true survival curve is included as reference     
lines(S5$t, 1 - exp(2*x0)/(1 + exp(2*x0)) + exp(2*x0)/(1 + exp(2*x0))*
            (1 - pweibull(S5$t, shape = 0.5 * (x0 + 4))), lwd = 2)
    }
}
\references{
    Beran, R. (1981). Nonparametric regression with randomly censored survival
    data. Technical Report. Berkeley, University of California.
    
    Safari, W. C., López-de-Ullibarri I., Jácome, M. A. (2021). A product-limit
    estimator of the conditional survival function when cure status is partially
    known. Biometrical Journal, 63(5): 984-1005. \doi{10.1002/bimj.202000173}.
}
\seealso{\code{\link{controlpars}}
}
