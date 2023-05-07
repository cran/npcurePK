controlpars <- function(b = 100L, hbound = c(0.1, 3), hl = 30L, hgrid_save = FALSE,
                        nnfrac = 0.25, fpilot = NULL, qt = 0.9, ncores = 1L,
                        seed = NULL, ...) {
    stopifnot("Incorrect 'hbound' parameter" = {
        length(hbound) == 2
        any(hbound > 0)
    }
    )
    stopifnot("Incorrect 'hl' parameter" = hl > 1)
    stopifnot("Incorrect 'nnfrac' parameter" = {
        nnfrac > 0
        nnfrac < 1
    }
    )
    stopifnot("Incorrect 'qt' parameter" = {
        qt > 0
        nnfrac < 1
    }
    )
    return(list(b = as.integer(b), hbound = as.numeric(hbound),
                hl = as.integer(hl), hgrid_save = hgrid_save, nnfrac = nnfrac,
                fpilot = fpilot, qt = qt, ncores = ncores, seed = seed,
                dots = list(...)))
}

prodlim_curepk_boot <- function(dfr, x0, bootpars) {
    prodlim_curepk_inner <- function(dfr, x0, g) {
        n_row <- nrow(dfr)
        surv <- numeric(length = n_row)
        res <- 1
        xg <- (x0 - dfr$x) / g
        k <- (abs(xg) < 1) * 0.75 * (1 - xg^2)
        bcx <- sum(k * (dfr$xinu == 1))
        sumbx <- rev(cumsum(rev(k * (dfr$xinu == 0))))
        for (i in 1:n_row) {
            denom <- sumbx[i] + bcx
            if (!is.na(denom) && denom > 0 && (dfr$d[i] == 1))
                res <- res * (1 - k[i] / denom)
            surv[i] <- res
        }
        return(surv)
    }
    prodlim_curepk_boot_mise <- function(dfr, x0, g, bootpars, hgrid) {
        n_row <- nrow(dfr)
        survgt <- prodlim_curepk_inner(dfr, x0, g)
        probg <- survgt[n_row]
        dfrboot <- data.frame(x = dfr$x, t = 0, d = 0, xinu = 0)
        if (probg == 1) {
            surv0g <- 0
            tmax <- quantile(dfr[, 2], probs = 1 - bootpars$qt)[[1]]
        } else {
            presurv0g  <- (survgt - probg) / (1 - probg)
            surv0g <- ifelse(presurv0g > 0 & dfr$t > tail(dfr$t[dfr$d == 1], 1),
                             0, ifelse(presurv0g < 0, 0, presurv0g))
            tmax <- data.table::first(dfr$t
                                      [surv0g == data.table::
                                       first(DescTools::Closest(surv0g,
                                                                1 - bootpars$qt,
                                                                na.rm = TRUE))])
        }
        tgrid <- seq(0, tmax, length.out = 100L)
        weight_m <- outer(dfr$x, dfr$x, function(x, y, k, h) k((x - y) / h),
                         k = function(x) (abs(x) < 1) * 0.75 * (1 - x^2), h = g)
        for (i in 1:n_row) {
            dfrboot[i, 2:4] <- dfr[sample(n_row, size = 1,
                                          prob = weight_m[i, ]), 2:4]
        }
        dfrboot_ord <- dfrboot[order(dfrboot$t, -dfrboot$d, dfrboot$xinu), ]
        lhgrid <- length(hgrid)
        bootmise <- numeric(lhgrid)
        for (ih in 1:lhgrid) {
            sht_boot <- prodlim_curepk_inner(dfrboot_ord, x0, hgrid[ih])
            sht_boot_grid <- approxfun(x = c(0, dfrboot_ord$t),
                                       y = c(1, sht_boot), method = "constant",
                                       f = 0, rule = 2, ties = "ordered")
            sgt_grid <- approxfun(x = c(0, dfr$t), y = c(1, survgt),
                                  method = "constant", f = 0, rule = 2,
                                  ties = "ordered")
            diff <- (sht_boot_grid(tgrid) - sgt_grid(tgrid))^2
            diff_grid <- approxfun(x = tgrid, y = diff, method = "constant",
                                   f = 0, rule = 2, ties = "ordered")
            bootmise[ih] <- bootmise[ih] +
                as.numeric(integrate(diff_grid, 0 + .Machine$double.eps^0.3,
                                     tmax, rel.tol = 0.1,
                                     subdivisions = 2000)[1])
        }
        return(bootmise)
    }
    dfr <- dfr[order(dfr$t, -dfr$d, dfr$xinu), ]
    ncores <- bootpars$ncores
    seed <- bootpars$seed
    if (is.null(seed))
        seed <- 1L
    lx0 <- length(x0)
    b <- bootpars$b
    hbound <- bootpars$hbound
    lhgrid <- bootpars$hl
    steph <- (hbound[2] / hbound[1])^(1 / lhgrid)
    hgrid <- IQR(dfr$x) / 1.349 * hbound[1] * steph^seq(0, lhgrid,
                                                        length.out = lhgrid)
    mise <- matrix(0, lhgrid, lx0)
    min_mise <- hboot <- numeric(lx0)
    if (ncores == 1L) {
        set.seed(seed)
        for (ix0 in 1:lx0) {
            if (is.null(bootpars$fpilot))
                pilot <- npcure::hpilot(dfr$x, x0[ix0], bootpars$nnfrac)
            for (ib in 1:b) {
                mise[, ix0] <- mise[, ix0] +
                    prodlim_curepk_boot_mise(dfr, x0[ix0], pilot, bootpars,
                                             hgrid)
            }
            min_mise[ix0] <- min(mise[, ix0])
            hboot[ix0] <- hgrid[mean(which(min_mise[ix0] == mise[, ix0]))]
        }
    } else {
        for (ix0 in 1:lx0) {
            if (is.null(bootpars$fpilot))
                pilot <- npcure::hpilot(dfr$x, x0[ix0], bootpars$nnfrac)
            cl <- parallel::makeCluster(ncores)
            doParallel::registerDoParallel(cl)
            out <- foreach::foreach(ib = 1:b) %dopar% {
                set.seed(seed + (ix0 - 1) * b + ib - 1)
                prodlim_curepk_boot_mise(dfr, x0[ix0], pilot, bootpars, hgrid)
            }
            mise[, ix0] <- Reduce("+", out)
            min_mise[ix0] <- min(mise[, ix0])
            hboot[ix0] <- hgrid[mean(which(min_mise[ix0] == mise[, ix0]))]
            parallel::stopCluster(cl)
        }
    }
    return(hboot)
}

prob_curepk_boot <- function(dfr, x0, bootpars) {
    prob_curepk_inner <- function(dfr, x0, g) {
        n_row <- nrow(dfr)
        lx0 <- length(x0)
        res <- rep(1, lx0)
        for (ix0 in 1:lx0) {
            xg <- (x0[ix0] - dfr$x) / g[ix0]
            k <- (abs(xg) < 1) * 0.75 * (1 - xg^2)
            bcx <- sum(k * (dfr$xinu == 1))
            sumbx <- rev(cumsum(rev(k * (dfr$xinu == 0))))
            for (i in 1:n_row) {
                denom <- sumbx[i] + bcx
                if (!is.na(denom) && denom > 0 && (dfr$d[i] == 1))
                    res[ix0] <- res[ix0] * (1 - k[i] / denom)
            }
        }
        return(res)
    }
    prob_curepk_boot_mse <- function(dfr, x0, g, bootpars, hgrid) {
        n_row <- nrow(dfr)
        pg <- prob_curepk_inner(dfr, x0, g)
        dfrboot <- data.frame(x = dfr$x, t = 0, d = 0, xinu = 0)
        weight_m <- outer(dfr$x, dfr$x, function(x, y, k, h) k((x - y) / h),
                         k = function(x) (abs(x) < 1) * 0.75 * (1 - x^2), h = g)
        for (i in 1:n_row) {
            dfrboot[i, 2:4] <- dfr[sample(nrow(dfr), size = 1,
                                          prob = weight_m[i, ]), 2:4]
        }
        dfrboot_ord <- dfrboot[order(dfrboot$t, -dfrboot$d, dfrboot$xinu), ]
        bootmse <- numeric(lhgrid)
        for (ih in 1:lhgrid) {
            ph <- prob_curepk_inner(dfrboot_ord, x0, hgrid[ih])
            bootmse[ih] <- bootmse[ih] + (ph - pg)^2
        }
        return(bootmse)
    }
    dfr <- dfr[order(dfr$t, -dfr$d, dfr$xinu), ]
    ncores <- bootpars$ncores
    seed <- bootpars$seed
    if (is.null(seed))
        seed <- 1L
    lx0 <- length(x0)
    b <- bootpars$b
    hbound <- bootpars$hbound
    lhgrid <- bootpars$hl
    steph <- (hbound[2] / hbound[1])^(1 / lhgrid)
    hgrid <- IQR(dfr$x) / 1.349 * hbound[1] * steph^seq(0, lhgrid,
                                                        length.out = lhgrid)
    mse <- matrix(0, lhgrid, lx0)
    min_mse <- hboot <- numeric(lx0)
    if (ncores == 1L) {
        set.seed(seed)
        for (ix0 in 1:lx0) {
            if (is.null(bootpars$fpilot))
                pilot <- npcure::hpilot(dfr$x, x0[ix0], bootpars$nnfrac)
            for (ib in 1:b) {
                mse[, ix0] <- mse[, ix0] +
                    prob_curepk_boot_mse(dfr, x0[ix0], pilot, bootpars, hgrid)
            }
            min_mse[ix0] <- min(mse[, ix0])
            hboot[ix0] <- hgrid[mean(which(min_mse[ix0] == mse[, ix0]))]
        }
    } else {
        for (ix0 in 1:lx0) {
            if (is.null(bootpars$fpilot))
                pilot <- npcure::hpilot(dfr$x, x0[ix0], bootpars$nnfrac)
            cl <- parallel::makeCluster(ncores)
            doParallel::registerDoParallel(cl)
            out <- foreach::foreach(ib = 1:b) %dopar% {
                set.seed(seed + (ix0 - 1) * b + ib - 1)
                prob_curepk_boot_mse(dfr, x0[ix0], pilot, bootpars, hgrid)
            }
            mse[, ix0] <- Reduce("+", out)
            min_mse[ix0] <- min(mse[, ix0])
            hboot[ix0] <- hgrid[mean(which(min_mse[ix0] == mse[, ix0]))]
            parallel::stopCluster(cl)
        }
    }
    return(hboot)
}

latency_curepk_boot <- function(dfr, x0, bootpars) {
    latency_curepk_boot_mise <- function(dfr, x0, hs, hp, bootpars, hgrid) {
        n_row <- nrow(dfr)
        x <- dfr$x
        t <- dfr$t
        d <- dfr$d
        xinu <- dfr$xinu
        survg <- prodlim_curepk(x, t, d, xinu, dfr, x0, hs)$surv
        probg <- prob_curepk(x, t, d, xinu, dfr, x0, hp)$prob_cure
        dfrboot <- data.frame(x = x, t = 0, d = 0, xinu = 0)
        if (probg == 1) {
            surv0g <- 0
            tmax <- quantile(dfr[, 2], probs = 1 - bootpars$qt)[[1]]
        } else {
            presurv0g  <- (survg - probg) / (1 - probg)
            surv0g <- ifelse(presurv0g > 0 & t > tail(t[d == 1], 1),
                             0, ifelse(presurv0g < 0, 0, presurv0g))
            tmax <- data.table::first(t
                                      [surv0g == data.table::
                                       first(DescTools::Closest(surv0g,
                                                                1 - bootpars$qt,
                                                                na.rm = TRUE))])
        }
        tgrid <- seq(0, tmax, length.out = 100L)
        weight_m <- outer(x, x, function(x, y, k, h) k((x - y) / h),
                          k = function(x) (abs(x) < 1) * 0.75 * (1 - x^2),
                          h = hs)
        for (i in 1:n_row) {
            dfrboot[i, 2:4] <- dfr[sample(n_row, size = 1,
                                          prob = weight_m[i, ]), 2:4]
        }
        dfrboot_ord <- dfrboot[order(dfrboot$t, -dfrboot$d, dfrboot$xinu), ]
        bootmise <- matrix(0, lhgrid, lhgrid)
        for (ih in 1:lhgrid) {
            for (ij in 1:lhgrid) {
                survh <- prodlim_curepk(x, t, d, xinu, dfr, x0, hgrid[ij])$surv
                probh <- prob_curepk(x, t, d, xinu, dfr, x0,
                                     hgrid[ih])$prob_cure
                if (probh == 1) {
                    surv0h <- 0
                } else {
                    presurv0h  <- (survh - probh) / (1 - probh)
                    surv0h <- ifelse(presurv0h > 0 & dfrboot_ord[, 2] >
                                     tail(dfrboot_ord[, 2]
                                          [dfrboot_ord[, 3] == 1], 1), 0,
                              ifelse(presurv0h < 0, 0, presurv0h))
                }
                surv0h_grid <- approxfun(x = c(0, dfrboot_ord[, 2]),
                                         y = c(1, surv0h), method = "constant",
                                         f = 0, rule = 2, ties = "ordered")
                surv0g_grid <- approxfun(x = c(0, t), y = c(1, surv0g),
                                         method = "constant", f = 0, rule = 2,
                                         ties = "ordered")
                diff <- (surv0h_grid(tgrid) - surv0g_grid(tgrid))^2
                diff_grid <- approxfun(x = tgrid, y = diff, method = "constant",
                                       f = 0, rule = 2, ties = "ordered")
                bootmise[ih, ij] <-
                    as.numeric(integrate(diff_grid, 0 + .Machine$double.eps^0.3,
                                         tmax, rel.tol = 0.1,
                                         subdivisions = 2000)[1])
            }
        }
        return(bootmise)
    }
    dfr <- dfr[order(dfr$t, -dfr$d, dfr$xinu), ]
    ncores <- bootpars$ncores
    seed <- bootpars$seed
    if (is.null(seed))
        seed <- 1L
    lx0 <- length(x0)
    b <- bootpars$b
    hbound <- bootpars$hbound
    lhgrid <- bootpars$hl
    steph <- (hbound[2] / hbound[1])^(1 / lhgrid)
    hgrid <- IQR(dfr$x) / 1.349 * hbound[1] * steph^seq(0, lhgrid,
                                                        length.out = lhgrid)
    mise <- array(0, dim = c(lhgrid, lhgrid, lx0))
    hboot <- matrix(0, nrow = 2, ncol = lx0)
    if (ncores == 1L) {
        set.seed(seed)
        for (ix0 in 1:lx0) {
            hsurv_pilot <- prodlim_curepk_boot(dfr, x0[ix0], bootpars)
            hprob_pilot <- prob_curepk_boot(dfr, x0[ix0], bootpars)
            for (ib in 1:b) {
                mise[, , ix0] <- mise[, , ix0] +
                    latency_curepk_boot_mise(dfr, x0[ix0], hsurv_pilot,
                                             hprob_pilot, bootpars, hgrid)
            }
            min_mise <- min(mise[, , ix0])
            ind <- colMeans(which(min_mise == mise[, , ix0], arr.ind = TRUE))
            hboot[, ix0] <- hgrid[ind]
        }
    } else {
        for (ix0 in 1:lx0) {
            hsurv_pilot <- prodlim_curepk_boot(dfr, x0[ix0], bootpars)
            hprob_pilot <- prob_curepk_boot(dfr, x0[ix0], bootpars)
            cl <- parallel::makeCluster(ncores)
            doParallel::registerDoParallel(cl)
            out <- foreach::foreach(ib = 1:b,
                                    .export = c("prodlim_curepk",
                                                "prob_curepk")) %dopar% {
                set.seed(seed + (ix0 - 1) * b + ib - 1)
                latency_curepk_boot_mise(dfr, x0[ix0], hsurv_pilot,
                                         hprob_pilot, bootpars, hgrid)
            }
            mise[, , ix0] <- Reduce("+", out)
            min_mise <- min(mise[, , ix0])
            ind <- colMeans(which(min_mise == mise[, , ix0], arr.ind = TRUE))
            hboot[, ix0] <- hgrid[ind]
            parallel::stopCluster(cl)
        }
    }
    return(hboot)
}
