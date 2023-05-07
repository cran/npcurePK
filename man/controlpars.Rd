\name{controlpars
}
\alias{controlpars
}
\title{Control Values for the Bootstrap
}
\description{This function returns a list of values for the control parameters
    of the functions that will be used for the bootstrap bandwidth
    selector.
}
\usage{
    controlpars(b = 100L, hbound = c(0.1, 3), hl = 30L, hgrid_save = FALSE,
                nnfrac = 0.25, fpilot = NULL, qt = 0.9, ncores = 1L,
                seed = NULL, \dots)
}
\arguments{
    \item{b}{An integer giving the number of bootstrap resamples, \code{100} by
        default.
    }   
    \item{hbound}{A numeric vector of length 2 specifying the minimum
        (default, \code{0.1}) and maximum (default, \code{3}), respectively,
        of the initial grid of bandwidths as a multiple of the standardized
        interquartile range of the covariate values.
    }
    \item{hl}{A numeric value giving the length of the initial grid of
        bandwidths. The default is \code{10}.
    }
    \item{hgrid_save}{A logical value specifying if the grids of bandwidths
        must be saved as a component of the list returned by the
        \code{prodlim_curepk_boot} function. The default is \code{FALSE}.
    }
    \item{nnfrac}{A numeric value giving the fraction of the sample size that
        determines the order of the nearest neighbor used when choosing the
        pilot bandwidth. The default is \code{0.25}.
    }
    \item{fpilot}{A function name or \code{NULL}. If \code{NULL}, the default,
        the \code{hpilot} function is used for computing a pilot bandwidth in
        case that one is needed. If not \code{NULL}, it must be the name of a
        user-defined function, given as a function name or as a character
        string. This function must necessarily have an argument \code{x0},
        playing the same role than in \code{hpilot}, and must return a value of
        the same length than \code{x0}. If \code{fpilot} has more arguments,
        they are passed through the \dots argument (see below).
    }
    \item{qt}{In bandwidth selection for the product-limit estimator, a numeric
        value specifying the order of a quantile of the observed times. It
        determines the right boundary of the integration interval in the
        computation of the ISE (the lower limit is \eqn{0}). The default is
        \code{0.9} (90th quantile).
    }
    \item{ncores}{The number of cores used in parallel computations.
    }
    \item{seed}{An optional integer passed to \code{set.seed()} to set the
        randomization seed.
    }
    \item{\dots}{Arguments of \code{fpilot}, if \code{fpilot} is not
        \code{NULL}.
    }
}
\details{The output of \code{controlpars} is a list of control parameters
    required by the functions which use the bootstrap. This is mainly the case
    of the \code{prodlimcurePKhboot} function, which compute the bootstrap
    bandwidth selectors of the estimators of the survival, latency and the
    probability of cure. Since these functions are indirectly called by
    \code{prodlimcurePKhboot} function when their \code{h} argument is missing,
    the output of \code{controlpars} is also the expected (and default) way of
    passing to them the parameters for bandwidth selection.
}
