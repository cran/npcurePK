\name{prob_curepk
}
\alias{prob_curepk
}
\title{Compute Estimator of Cure Probability when Cure Status is Partially Known
}
\description{This function computes the nonparametric estimator of the cure
    probability when cure status is partially known proposed by Safari
    \emph{et al} (2022).
}
\usage{
    prob_curepk(x, t, d, xinu, dataset, x0, h, local = TRUE,
                bootpars = if (!missing(h)) NULL else controlpars())
}
\arguments{
    \item{x}{If \code{dataset} is missing, a numeric object giving the covariate
        values. If \code{dataset} is a data frame, it is interpreted as the name
        of the variable corresponding to the covariate in the data frame.}
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
        cure probability will be computed.}
    \item{h}{A numeric vector of bandwidths.}
    \item{local}{A logical value, \code{TRUE} by default, specifying whether
        \code{local} or \code{global} bandwidths are used.}
    \item{bootpars}{A list of parameters controlling the bootstrap when
        computing the bootstrap bandwidths of the cure probability estimator.
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
        cure probability is computed.}
    \item{prob_cure}{The estimate of the cure probability \code{1-p(x0)} with
        bandwidth \code{h}. It is a vector of the same length as \code{x0}.}
}
\details{Mixture cure model writes the conditional survival function
    \eqn{S(t\mid x)=P(Y>t\mid X=x)}{S(t|x)=P(Y>t|X=x)} as
    \deqn{S(t\mid x)=1-p(x)+p(x)S_0(t\mid x)}{S(t|x)=1-p(x)+p(x)S_0(t|x)} where
    \eqn{1-p(x)=P(Y=\infty\mid X=x)}{1-p(x)=P(Y=infinity|X=x)} is the
    probability of cure.

    This function computes the kernel estimator of the probability of cure
    \eqn{1-p(x)} in Safari \emph{et al} (2022). It is based on the previous
    relationship and the generalized product-limit estimator of the  conditional
    survival function \eqn{S(t\mid x)}{S(t|x)} in Safari \emph{et al} (2021),
    using the Nadaraya-Watson weights, when the cure status is partially known.
    If there are not individuals known to be cured (\code{xinu=0}), then the
    nonparametric estimator of the cure rate in López-Cheda \emph{et al} (2017)
    is computed.
    
    The Epanechnikov kernel is used. If the smoothing parameter \code{h} is not
    provided, then the bootstrap bandwidth selector in Safari \emph{et al}
    (2022) is used. The function is available only for one continuous covariate
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

## Cure rate estimates for one single covariate value x0 = 0 and using ...
## ... (a) one single fixed bandwidth h = 0.5 
p1 <- prob_curepk(x, t, d, xinu, data, x0 = 0, 
                  h = 0.5, local = TRUE)

## ... (b) a vector of bandwidths h = c(0.25, 0.5, 0.75, 1) 
p2 <- prob_curepk(x, t, d, xinu, data, x0 = c(0, 0, 0, 0), 
                  h = c(0.25, 0.5, 0.75, 1), local = TRUE)

## ... (c) a bootstrap bandwidth (the default when the bandwidths
## argument h is not provided).
## The bootstrap bandwidth is searched in a grid of 10 bandwidths (hl = 10)
## between 0.2 and 2 times the standardized interquartile range of the
## covariate values (hbound = c(0.1, 3)).
(p3 <- prob_curepk(x, t, d, xinu, data, x0 = 0))
## Equivalently
\donttest{    
    (p3 <- prob_curepk(x, t, d, xinu, data, x0 = 0,
                       bootpars = controlpars(hl = 10, hbound = c(0.1, 3))))
}

## Cure rate estimates for a vector of 20 covariate values and using ...
x0 = seq(from = min(data$x), to = max(data$x), length.out = 15)
## ... (a) one single fixed bandwidth h = 0.5 
p4 <- prob_curepk(x, t, d, xinu, data, x0 = x0, h = 0.5, local = FALSE)
## Plot predicted cure probabilities for covariate values x0 and bandwidths
## h = 0.5
plot(p4$x0, p4$prob_cure, xlab = "Covariate X", type = "l",
     ylab = "Probability of cure", ylim = c(0, 1))
## The true cure rate is included as reference
lines(p4$x0, 1 - exp(2*x0)/(1 + exp(2*x0)), lwd = 2)

## ... (b) a vector of bandwidths h = c(0.5, 0.75, 1) 
p5 <- prob_curepk(x, t, d, xinu, data, x0 = x0, h = c(0.5, 0.75, 1),
                  local = FALSE)
## Plot predicted cure probabilities for covariate values x0 and bandwidths
## h = 0.5
plot(p5$x0, p5$prob_cure[1, ], xlab = "Covariate X", type = "l",
     ylab = "Probability of cure", ylim = c(0, 1))
## The estimates with bandwidth h = 0.75 and h = 1 are added
lines(p5$x0, p5$prob_cure[2, ])
lines(p5$x0, p5$prob_cure[3, ])
## The true cure rate is included as reference
lines(p5$x0, 1 - exp(2*x0)/(1 + exp(2*x0)), lwd = 2)

\donttest{
    ## ... (c) the bootstrap bandwidth
    (p6 <- prob_curepk(x, t, d, xinu, data, x0 = x0,
                       bootpars = controlpars(b = 50, ncores = 2, seed = 123)))
    ## Plot predicted cure probabilities for covariate values x0 and bootstrap
    ## bandwidths
    plot(p6$x0, p6$prob_cure, xlab = "Covariate X", type = "l",
         ylab = "Probability of cure", ylim = c(0, 1))
    ## The true cure rate is included as reference
    lines(p6$x0, 1 - exp(2*x0)/(1 + exp(2*x0)), lwd = 2)
}
}
\references{
    Beran, R. (1981). Nonparametric regression with randomly censored survival
    data. Technical Report. Berkeley, University of California.
    
    López-Cheda, A. Cao, R., Jácome, M.A., Van Keilegom, I. (2017).
    Nonparametric incidence estimation and bootstrap bandwidth selection in
    mixture cure models. Computational Statistics and Data Analysis 105:144-165.
    \doi{10.1016/j.csda.2016.08.002}.

    Safari, W. C., López-de-Ullibarri I., Jácome, M. A. (2021). A product-limit
    estimator of the conditional survival function when cure status is
    partially known. Biometrical Journal, 63(5): 984-1005.
    \doi{10.1002/bimj.202000173}.

    Safari, W. C., López-de-Ullibarri I., Jácome, M. A. (2022). Nonparametric
    kernel estimation of the probability of cure in a mixture cure model when
    the cure status is partially observed. Statistical Methods in Medical
    Research, 31(11):2164-2188. 
    \doi{10.1177/09622802221115880}.
}
\seealso{\code{\link{controlpars}}
}
