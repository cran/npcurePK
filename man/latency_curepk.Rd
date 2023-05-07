\name{latency_curepk
}
\alias{latency_curepk
}
\title{Compute Estimator of Latency Function  when Cure Status is Partially
    Known
}
\description{This function computes the nonparametric estimator of the latency
    function when cure status is partially known proposed by Safari \emph{et al}
    (2023).
}
\usage{
    latency_curepk(x, t, d, xinu, dataset, x0, h, local = TRUE,
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
        latency function will be computed.}
    \item{h}{A numeric matrix of bandwidths.}
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

    \item{h}{The numeric matrix \code{(2 x length(x0))} of bandwidths used in
        the estimation. One bandwidth \code{h[1, ]} is used for the estimation
        of \eqn{1-p(x)} and another bandwidth \code{h[2, ]} is used for the
        estimation of \eqn{S(t \mid x)}{S(t|x)}. If \code{h} argument is
        missing, the bootstrap bandwidth computed with the control parameters
        in argument \code{bootpars}.}
    \item{x0}{The numeric vector of covariate values where the estimate of the
        latency function is computed.}
    \item{prob_cure}{The estimate of the cure probability \code{1 - p(x0)}
        with bandwidth \code{h[1, ]}. It is a vector of the same length as
        \code{x0}.}
    \item{t}{The observed time values, where the latency function is estimated.}
    \item{surv}{Estimates of the survival function for each one of the
        covariate values specified by the \code{x0} argument and the bandwidths
        in \code{h[2, ]}. It is a matrix of dimension \eqn{n\times length(x0)}{(n, length(x0))}
        if local bandwidths or bootstrap bandwidths are used, or an array for
        global bandwidths instead.}
    \item{latency}{Estimates of the latency for each one of the covariate
        values specified by the \code{x0} argument and the bandwidths in
        \code{h}. It is a matrix of dimension \eqn{n\times length(x0)}{(n, length(x0))}
        if local bandwidths or bootstrap bandwidths are used, or an array for
        global bandwidths instead.}
}
\details{This function computes an estimator of the latency function
    \eqn{S_0(t\mid x)=P(Y>t\mid Y<\infty, X=x)}{S_0(t|x)=P(Y>t|Y<infinity, X=x)}
    when the cure status is partially known, introduced in Safari \emph{et al}
    (2023). It is based on the relationship \deqn{S(t\mid x)=1-p(x)+p(x)S_0(t\mid x)}{S(t|x)=1-p(x)+p(x)S_0(t| x)},
    using the kernel estimator of the cure rate \eqn{1-p(x)} in Safari
    \emph{et al} (2022) and the survival function \eqn{S(t\mid x)}{S(t|x)}
    in Safari \emph{et al} (2021), with Nadaraya-Watson weights and bandwidth
    \eqn{h_1} for the cure rate and \eqn{h_2} for the survival function. If
    there are not individuals known to be cured (\code{xinu=0}), then the
    kernel estimator of the cure rate in López-Cheda \emph{et al} (2017) is
    computed.

    The latency estimator is computed with the pair of bandwidths in \code{h}.
    One bandwidth \code{h[1, ]} is used for the estimation of \eqn{1-p(x)} and
    another bandwidth \code{h[2, ]} is used for the estimation of
    \eqn{S(t\mid x)}{S(t|x)}. If the smoothing parameter \code{h} is not
    provided, then the bootstrap bandwidth selector in Safari \emph{et al}
    (2023) is used. The kernel considered is Epanechnikov kernel. The function
    is available only for one continuous covariate \eqn{X}.
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

## Latency estimates for one single covariate value x0 = 0 and using...
x0 <- 0

## ... (a) one single fixed bandwidth h = [1.1, 1] 
## h[1,] = 1.1 is used for estimating p(x) at x0 = 0
## h[2,] = 1 is used for estimating S(t|x) at x0 = 0
## The latency estimates are saved in an array (n × 1)
S0_1 <- latency_curepk(x, t, d, xinu, data, x0 = 0, 
                       h = matrix(c(1.1, 1), nrow = 2, ncol = 1, byrow = TRUE),
                       local = TRUE)
## Plot predicted latency curve for covariate value x0 = 0 and bandwidths
## h = [1.1, 1]
plot(S0_1$t, S0_1$latency, type = "s", xlab = "Time",
     ylab = "Latency function", ylim = c(0, 1))
## The true latency function is included as reference     
lines(S0_1$t, 1 - pweibull(S0_1$t, shape = 0.5 * (x0 + 4)))

## ... (b) two fixed bandwidths h = [1.1, 1] and h = [1.5, 2]
## One estimate of the latency S0(t|x0 = 0) is obtained using h[1, 1] = 1.1
## for estimating p(x) and h[2,1] = 1 for estimating S(t|x)
## Second estimate of the latency S0(t|x0 = 0) is obtained using h[1, 2] = 1.5
## using h[1,2] = 1.5 for estimating p(x) and h[2,2] = 2 for estimating S(t|x)
## The estimates are saved in an array (n × 2)
S0_2 <- latency_curepk(x, t, d, xinu, data, x0 = c(0, 0), 
                       h = matrix(c(1.1, 1, 1.5, 2), nrow = 2, ncol = 2,
                                  byrow = FALSE), local = TRUE)
## Plot predicted latency curve for covariate value x0 = 0 and bandwidths
## h = [1.1, 1] and and h = [1.5, 2]
plot(S0_2$t, S0_2$latency[, 1], type = "s", xlab = "Time",
     ylab = "Latency function", ylim = c(0, 1))
lines(S0_2$t, S0_2$latency[, 2], type = "s", lwd = 2)
## The true latency function is included as reference     
lines(S0_2$t, 1 - pweibull(S0_2$t, shape = 0.5 * (x0 + 4)))

\donttest{
    ## ... (c) with the bootstrap bandwidth selector (the default when the 
    ## bandwidth argument h is not provided).
    ## The bootstrap bandwidth is searched with parallel computation 
    ## (ncores = 2) in a grid of 9 bandwidths (hl = 9) between 0.2 and 2 times
    ## the standardized interquartile range of the covariate values
    ## (hbound = c(0.1, 2)). The latency estimates are saved in an array of
    ## dimension (n, 1)
    library(doParallel)
    (S0_3 <- latency_curepk(x, t, d, xinu, data, x0 = 0,
                            bootpars = controlpars(b = 50, hl = 9, 
                            hbound = c(0.1, 2), ncores = 2)))
    plot(S0_3$t, S0_3$latency[, 1], type = "s", xlab = "Time",
         ylab = "Latency function", ylim = c(0, 1))                      
    ## The true latency function is included as reference     
    lines(S0_3$t, 1 - pweibull(S0_3$t, shape = 0.5 * (x0 + 4)))
}
}
\references{
    López-Cheda, A., Jácome, M.A., Cao, R. (2017). Nonparametric latency
    estimation for mixture cure models. TEST 26:353-376.
    \doi{10.1007/s11749-016-0515-1}.

    Safari, W. C., López-de-Ullibarri I., Jácome, M. A. (2021). A product-limit
    estimator of the conditional survival function when cure status is
    partially known. Biometrical Journal, 63(5): 984-1005.
    \doi{10.1002/bimj.202000173}.

    Safari, W. C., López-de-Ullibarri I., Jácome, M. A. (2022). Nonparametric
    kernel estimation of the probability of cure in a mixture cure model when
    the cure status is partially observed. Statistical Methods in Medical
    Research, 31(11):2164-2188.
    \doi{10.1177/09622802221115880}.

    Safari, W. C., López-de-Ullibarri I., Jácome, M. A. (2023). Latency
    function estimation under the mixture cure model when the cure status is
    available. Lifetime Data Analysis.
    \doi{10.1007/s10985-023-09591-x}.
}
\seealso{\code{\link{controlpars}}
}
