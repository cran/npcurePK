## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, fig.width = 4.8, fig.height = 4.8)

## -----------------------------------------------------------------------------
library(npcurePK)

## -----------------------------------------------------------------------------
S <- prodlim_curepk(x, t, d, xinu, sarcoma, x0 = c(40, 90))
plot(S$t, S$surv[, 1], type = "s", xlab = "Time", ylab = "Survival probability", ylim = c(0, 1))
lines(S$t, S$surv[, 2], type = "s")  

## -----------------------------------------------------------------------------
S0 <- latency_curepk(x, t, d, xinu, sarcoma, x0 = 60,
                      bootpars = controlpars(b = 50, hl = 10, ncores = 2, hbound = c(0.1, 2)))
plot(S0$t, S0$latency[, 1], type = "s", xlab = "Time", ylab = "Latency function", ylim = c(0, 1))

## -----------------------------------------------------------------------------
x0 <-  seq(from = min(sarcoma$x), to = max(sarcoma$x), length.out = 50)

## -----------------------------------------------------------------------------
p <- prob_curepk(x, t, d, xinu, sarcoma, x0 = x0, h = c(20, 25, 30), local = FALSE)

## -----------------------------------------------------------------------------
plot(p$x0, p$prob_cure[1, ], xlab = "Age", type = "l", ylab = "Probability of cure", ylim = c(0, 1)) 
lines(p$x0, p$prob_cure[2, ], lwd = 1.5)
lines(p$x0, p$prob_cure[3, ], lwd = 3)

## -----------------------------------------------------------------------------
library(doParallel)
p2 <- prob_curepk(x, t, d, xinu, sarcoma, x0 = x0,
                   bootpars = controlpars(b = 50, ncores = 2, seed = 123))
plot(p2$x0, p2$prob_cure, xlab = "Age", type = "l", ylab = "Probability of cure", ylim = c(0, 1))

