\docType{data}
\name{sarcoma
}
\alias{sarcoma
}
\title{Sarcoma Dataset
}
\description{Sarcoma is a rare type of cancer that represents 1\% of all adult solid malignancies (Choy, 2014). If a tumor can be surgically removed to render the patient with sarcoma free of detectable disease, 5 years is the survival time at which sarcoma oncologists assume long-term remissions.
\code{sarcoma} dataset contains the observed survival time of 232 patients until death from sarcoma, and covariates such as the age at diagnosis. Patients tumor free for more than 5 years were assumed to be long-term survivors (known to be cured, \code{xinu = 1}).
}
\format{A data frame with 232 rows and 4 variables:
    \describe{
        \item{\code{x}}{Age (years) of patients at diagnosis.}
        \item{\code{t}}{Observed time until death from sarcoma.}
        \item{\code{d}}{Censoring status (0 = censored, 1 = death from sarcoma).
        }
        \item{\code{xinu}}{Cure status (0 = dead or unknown, 1 = tumor free for
                                        more than 5 years).
        }
    }
}
\source{Provided by the authors to serve as an example.
}
\usage{
    sarcoma
    data(sarcoma, package = "npcurePK")
}
\keyword{datasets
}
\references{
Choy, E. (2014). Sarcoma after 5 years of progression-free survival: Lessons from the French sarcoma group. Cancer, 120(19), 2942-2943.
}
\examples{
    \donttest{
## Survival estimates of patients aged 40 and 90 years old
## computed with bootstrap bandwidths
## (seed is used for bootstrap resampling)
(S1 <- prodlim_curepk(x, t, d, xinu, sarcoma, x0 = c(40, 90),
                      bootpars = controlpars(b = 100, ncores = 2, seed = 123)))
plot(S1$t, S1$surv[, 1], type = "s", xlab = "Time",
     ylab = "Survival probability", ylim = c(0, 1))
lines(S1$t, S1$surv[, 2], type = "s")
    }
}
