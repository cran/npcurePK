prodlim_curepk <- function(x, t, d, xinu, dataset, x0, h, local = TRUE,
                           bootpars = if (!missing(h)) NULL else controlpars()
                           ) {
    dfr <-
        if (missing(dataset))
            na.omit(data.frame(x, t, d, xinu))
        else
            na.omit(dataset[, c(deparse(substitute(x)), deparse(substitute(t)),
                                deparse(substitute(d)),
                                deparse(substitute(xinu)))])
    names(dfr) <- c("x", "t", "d", "xinu")
    stopifnot("'x' cannot be a factor" = !is.factor(dfr$x))
    n_row <- nrow(dfr)
    dfr <- dfr[order(dfr$t, -dfr$d, dfr$xinu), ]
    dfr$x <- as.numeric(dfr$x)
    dfr$t <- as.numeric(dfr$t)
    dfr$d <- as.integer(dfr$d)
    dfr$xinu <- as.integer(dfr$xinu)
    ordx0 <- order(x0)
    x0 <- as.numeric(x0[ordx0])
    lx0 <- length(x0)
    if (!local && missing(h))
        warning("Option 'local = FALSE' overridden with missing 'h'")
    if (missing(h)) {
        h <- prodlim_curepk_boot(dfr, x0, bootpars)
        if (bootpars$hgrid_save) {
            hgrid <- IQR(dfr$x) / 1.349 * bootpars$hbound[1] *
                ((bootpars$hbound[2] / bootpars$hbound[1])^(1 / bootpars$hl))^
                seq(0, bootpars$hl, length.out = bootpars$hl)
        }
    } else {
        if (local) {
            stopifnot("When 'local = TRUE', 'x0' and 'h' must have the same
 length" = lx0 == length(h))
            h <- as.numeric(h[ordx0])
        } else {
            h <- as.numeric(sort(h))
        }
    }
    lh <- length(h)
    if (local) {
        surv <- matrix(0, n_row, lx0)
        for (ix0 in 1:lx0) {
            res <- 1
            xh <- (x0[ix0] - dfr$x) / h[ix0]
            k <- (abs(xh) < 1) * 0.75 * (1 - xh^2)
            bcx <- sum(k * (dfr$xinu == 1))
            sumbx <- rev(cumsum(rev(k * (dfr$xinu == 0))))
            for (i in 1:n_row) {
                denom <- sumbx[i] + bcx
                if (!is.na(denom) && denom > 0 && (dfr$d[i] == 1))
                    res <- res * (1 - k[i] / denom)
                surv[i, ix0] <- res
            }
        }
    } else {
        surv <- array(0, dim = c(n_row, lx0, lh))
        res <- matrix(1, lh, lx0)
        for (ix0 in 1:lx0) {
            for (ih in 1:lh) {
                xh <- (x0[ix0] - dfr$x) / h[ih]
                k <- (abs(xh) < 1) * 0.75 * (1 - xh^2)
                bcx <- sum(k * (dfr$xinu == 1))
                sumbx <- rev(cumsum(rev(k * (dfr$xinu == 0))))
                for (i in 1:n_row) {
                    denom <- sumbx[i] + bcx
                    if (!is.na(denom) && denom > 0 && (dfr$d[i] == 1))
                        res[ih, ix0] <- res[ih, ix0] * (1 - k[i] / denom)
                    surv[i, ix0, ih] <- res[ih, ix0]
                }
            }
        }
        if (lx0 != 1) {
            surv <- drop(surv)
        } else {
            if (lh == 1) {
                surv <- as.matrix(surv)
            }
        }
    }
    if (local || (!local && lh == 1)) {
        dimnames(surv)[[2]] <- paste0("x0=", as.character(round(x0, 8)))
    } else {
        dimnames(surv)[[3]] <- paste0("h=", as.character(round(h, 8)))
        dimnames(surv)[[2]] <- paste0("x0=", as.character(round(x0, 8)))
    }
    if (!is.null(bootpars) && bootpars$hgrid_save) {
        return(list(h = h, hgrid = hgrid, x0 = x0, t = dfr$t, surv = surv))
    } else {
        return(list(h = h,  x0 = x0, t = dfr$t, surv = surv))
    }
}

prob_curepk <- function(x, t, d, xinu, dataset, x0, h, local = TRUE,
                        bootpars = if (!missing(h)) NULL else controlpars()) {
    dfr <-
        if (missing(dataset))
            na.omit(data.frame(x, t, d, xinu))
        else
            na.omit(dataset[, c(deparse(substitute(x)), deparse(substitute(t)),
                                deparse(substitute(d)),
                                deparse(substitute(xinu)))])
    names(dfr) <- c("x", "t", "d", "xinu")
    stopifnot("'x' cannot be a factor" = !is.factor(dfr$x))
    n_row <- nrow(dfr)
    dfr <- dfr[order(dfr$t, -dfr$d, dfr$xinu), ]
    dfr$x <- as.numeric(dfr$x)
    dfr$t <- as.numeric(dfr$t)
    dfr$d <- as.integer(dfr$d)
    dfr$xinu <- as.integer(dfr$xinu)
    ordx0 <- order(x0)
    x0 <- as.numeric(x0[ordx0])
    lx0 <- length(x0)
    if (!local && missing(h))
        warning("Option 'local = FALSE' overridden with missing 'h'")
    if (missing(h)) {
        h <- prob_curepk_boot(dfr, x0, bootpars)
        if (bootpars$hgrid_save) {
            hgrid <- IQR(dfr$x) / 1.349 * bootpars$hbound[1] *
                ((bootpars$hbound[2] / bootpars$hbound[1])^(1 / bootpars$hl))^
                seq(0, bootpars$hl, length.out = bootpars$hl)
        }
    } else {
        if (local) {
            stopifnot("When 'local = TRUE', 'x0' and 'h' must have the same
 length" = lx0 == length(h))
            h <- as.numeric(h[ordx0])
        } else {
            h <- as.numeric(sort(h))
        }
    }
    lh <- length(h)
    if (local) {
        res <- rep(1, lx0)
        for (ix0 in 1:lx0) {
            xh <- (x0[ix0] - dfr$x) / h[ix0]
            k <- (abs(xh) < 1) * 0.75 * (1 - xh^2)
            bcx <- sum(k * (dfr$xinu == 1))
            sumbx <- rev(cumsum(rev(k * (dfr$xinu == 0))))
            for (i in 1:n_row) {
                denom <- sumbx[i] + bcx
                if (!is.na(denom) && denom > 0 && (dfr$d[i] == 1))
                    res[ix0] <- res[ix0] * (1 - k[i] / denom)
            }
        }
    } else {
        res <- matrix(1, lh, lx0)
        for (ix0 in 1:lx0) {
            for (ih in 1:lh) {
                xh <- (x0[ix0] - dfr$x) / h[ih]
                k <- (abs(xh) < 1) * 0.75 * (1 - xh^2)
                bcx <- sum(k * (dfr$xinu == 1))
                sumbx <- rev(cumsum(rev(k * (dfr$xinu == 0))))
                for (i in 1:n_row) {
                    denom <- sumbx[i] + bcx
                    if (!is.na(denom) && denom > 0 && (dfr$d[i] == 1))
                        res[ih, ix0] <- res[ih, ix0] * (1 - k[i] / denom)
                }
            }
        }
        if (lx0 != 1)
            res <- drop(res)
    }
    if (local || (!local && lh == 1)) {
        res <- as.vector(res)
        names(res) <- paste0("x0=", as.character(round(x0, 8)))
    } else {
        dimnames(res)[[2]] <- paste0("x0=", as.character(round(x0, 8)))
        dimnames(res)[[1]] <- paste0("h=", as.character(round(h, 8)))
    }
    if (!is.null(bootpars) && bootpars$hgrid_save) {
        return(list(h = h, hgrid = hgrid, x0 = x0, prob_cure = res))
    } else {
        return(list(h = h, x0 = x0, prob_cure = res))
    }
}


latency_curepk <- function(x, t, d, xinu, dataset, x0, h, local = TRUE,
                           bootpars = if (!missing(h)) NULL else controlpars()
                           ) {
    dfr <-
        if (missing(dataset))
            na.omit(data.frame(x, t, d, xinu))
        else
            na.omit(dataset[, c(deparse(substitute(x)), deparse(substitute(t)),
                                deparse(substitute(d)),
                                deparse(substitute(xinu)))])
    names(dfr) <- c("x", "t", "d", "xinu")
    stopifnot("'x' cannot be a factor" = !is.factor(dfr$x))
    n_row <- nrow(dfr)
    dfr <- dfr[order(dfr$t, -dfr$d, dfr$xinu), ]
    dfr$x <- as.numeric(dfr$x)
    dfr$t <- as.numeric(dfr$t)
    dfr$d <- as.integer(dfr$d)
    dfr$xinu <- as.integer(dfr$xinu)
    ordx0 <- order(x0)
    x0 <- as.numeric(x0[ordx0])
    lx0 <- length(x0)
    if (!local && missing(h))
        warning("Option 'local = FALSE' overridden with missing 'h'")
    if (missing(h)) {
        h <- latency_curepk_boot(dfr, x0, bootpars)
        if (bootpars$hgrid_save) {
            hgrid <- IQR(dfr$x) / 1.349 * bootpars$hbound[1] *
                ((bootpars$hbound[2] / bootpars$hbound[1])^(1 / bootpars$hl))^
                seq(0, bootpars$hl, length.out = bootpars$hl)
        }
    } else {
        if (local) {
            stopifnot("When 'local = TRUE', the number of columns of 'h' must
 equal the length of 'x0'" = lx0 == ncol(h))
        }
    }
    probh <- prob_curepk(x, t, d, xinu, dfr, x0, h[1, ], local)$prob_cure
    survh <- prodlim_curepk(x, t, d, xinu, dfr, x0, h[2, ], local)$surv
    lh <- ncol(h)
    if (local || (!local && lh == 1)) {
        surv0h <- matrix(0, n_row, lx0)
        for (ix0 in 1:lx0) {
            if (probh[ix0] < 1) {
                presurv0h  <- (survh[, ix0] - probh[ix0]) / (1 - probh[ix0])
                surv0h[, ix0] <- ifelse(presurv0h > 0 & dfr$t >
                                        tail(dfr$t[dfr$d == 1], 1), 0,
                                 ifelse(presurv0h < 0, 0, presurv0h))
            }
        }
    } else {
        surv0h <- array(0, dim = c(n_row, lx0, lh))
        for (ix0 in 1:lx0) {
            for (ih in 1:lh) {
                if (probh[ih, ix0] < 1) {
                    presurv0h  <- (survh[, ix0, ih] - probh[ih, ix0]) /
                        (1 - probh[ih, ix0])
                    surv0h[, ix0, ih] <- ifelse(presurv0h > 0 & dfr$t >
                                                tail(dfr$t[dfr$d == 1], 1), 0,
                                         ifelse(presurv0h < 0, 0, presurv0h))
                }
            }
        }
    }
    if (local) {
        dimnames(h) <- list(c("h_prob", "h_surv"),
                            paste0("x0=", as.character(round(x0, 8))))
        dimnames(surv0h)[[2]] <- paste0("x0=", as.character(round(x0, 8)))
    } else {
        if (lx0 != 1)
            surv0h <- drop(surv0h)
        dimnames(h)[[1]] <- c("h_prob", "h_surv")
        if (lh == 1) {
          dimnames(surv0h)[[2]] <- paste0("x0=", as.character(round(x0, 8)))
        } else {
            dimnames(surv0h)[[3]] <- apply(h, 2, function(u) {
                paste0("(", paste0("h=", as.character(round(u, 8)), sep = "",
                                   collapse = ","), ")")
            }
            )
            dimnames(surv0h)[[2]] <- paste0("x0=", as.character(round(x0, 8)))
        }
    }
    if (!is.null(bootpars) && bootpars$hgrid_save) {
        return(list(h = h, hgrid = hgrid, x0 = x0, prob_cure = probh,
                    t = dfr$t, surv = survh, latency = surv0h))
    } else {
        return(list(h = h, x0 = x0, prob_cure = probh, t = dfr$t,
                    surv = survh, latency = surv0h))
    }
}
