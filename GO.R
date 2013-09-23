## unconstrained Gaussian ordination in one dimension

GO1 <-
    function(comm, tot = max(comm), freqlim = 5, parallel = 1, ...)
{
    ## do parallel?
    if(parallel > 1) {
        require(parallel) || stop("needs 'parallel' package")
        clu <- makeCluster(parallel)
        on.exit(stopCluster(clu))
    }
    ## Remove rare species
    freq <- colSums(comm > 0)
    comm <- comm[, freq >= freqlim]
    rs <- rowSums(comm)
    if (any(rs <= 0))
        warning(
            gettextf("%d SUs were empty and removed after applying 'freqlim'",
                     sum(rs <= 0)))
    message(gettextf("data has now %d SUs and %d species",
                     nrow(comm), ncol(comm)))
    ## initialize gradient as first DCA axis
    x <- scores(decorana(comm), display = "sites", choices = 1)
    ## loss function with quasibinomial glm
    if (parallel > 1)
        loss <- function(x, ...) {
            mods <- parLapply(clu, comm, 
                              function(y, ...)
                              glm(cbind(y, tot-y) ~ x + offset(-0.5 * x^2),
                                  family = quasibinomial))
            sum(sapply(mods, deviance))
        }
    else
        loss <- function(x, ...) {
            mods <- lapply(comm, function(y, ...)
                           glm(cbind(y, tot-y) ~ x + offset(-0.5 * x^2),
                               family = quasibinomial))
            sum(sapply(mods, deviance))
        }
    ## ML fit
    out <- nlm(loss, p = x, comm = comm, tot = tot, ...)
    out$data <- comm
    out$tot = tot
    class(out) <- "GO1"
    out
}

### GO2 tries to estimate species parameters and Gaussian axis
### simultaneously instead of alternating between estimating Gaussian
### axis and fitting species to the current estimate of the axis. This
### may be hard, but avoids slow glm() steps and can easily be
### generalized to multiple axes. For n SUs, m species, and k axes,
### there are m + k * (m + n) estimated parameters (and location of
### the species optima and the Gaussian axis really are are
### correlated). May be tough.
GO2 <-
    function(comm, k = 1, tot = max(comm), freqlim = 5, ...)
{
    ## Limit to k <= 4
    if (k > 4)
        stop(gettextf("Maximum allowed number of axes is 4, but you had k = %d",
                      k))
    ## Remove rare species
    freq <- colSums(comm > 0)
    comm <- comm[, freq >= freqlim]
    rs <- rowSums(comm)
    if (any(rs <= 0))
        warning(
            gettextf("%d SUs were empty and removed after applying 'freqlim'",
                     sum(rs <= 0)))
    message(gettextf("data has now %d SUs and %d species",
                     nrow(comm), ncol(comm)))
    ## define error family and get corresponding inverse link function
    ## ginv() and deviance function dev()
    fam <- quasibinomial()
    ginv <- fam$linkinv
    dev <- fam$dev.resids
    ## initialize gradient as DCA axes
    x <- as.matrix(scores(decorana(comm), display = "sites", choices = 1:k))
    ## initial estimates for species (these are "final" with the
    ## current 'x')
    mods <- lapply(comm, function(y) glm(cbind(y, tot-y) ~
                                         x + offset(-0.5 * rowSums(x^2)),
                                         family = fam))
    b <- sapply(mods, coef)
    null.spdev <- sapply(lapply(comm, function(y) glm(cbind(y, tot-y) ~ 1,
           family= fam)), function(z) z$deviance)
    null.deviance <- sum(null.spdev)
    ## Pack parameters to a single vector
    p <- c(as.vector(x), as.vector(t(b)))
    nn <- cumsum(c(k*nrow(comm), ncol(comm), k*ncol(comm)))
    ## We need matrix response, and with quasibinomial, we need to
    ## divide with 'tot'
    if (length(tot) < nrow(comm))
        tot <- rep(tot, length.out = nrow(comm))
    wts <- matrix(rep(tot, ncol(comm)), nrow(comm), ncol(comm))
    y <- as.matrix(comm) / wts
    ## Loss function fits the current model and return the sum of
    ## deviances
    loss <- function(p, ...) {
        ## split params
        x <- matrix(p[1 : nn[1]], ncol=k)
        b0 <- p[(nn[1]+1) : nn[2]]
        b1 <- matrix(p[(nn[2]+1) : nn[3]], nrow = k, byrow=TRUE)
        ## model lp = b0 -0.5*x^2 + b1*gr
        lp <- outer(-0.5*rowSums(x^2), b0, "+") + x %*% b1
        sum(dev(y, ginv(lp), wts))
    }
    mod <- nlm(loss, p = p, ...)
    out <- list(deviance = mod$minimum, null.deviance = null.deviance,
                k = k, iterations = mod$iterations, code = mod$code,
                rank = length(mod$gradient),
                df.residual = prod(dim(comm)) - length(mod$gradient),
                df.null = prod(dim(comm)) - ncol(comm))
    out$points <- matrix(mod$estimate[1 : nn[1]], ncol=k)
    specpara <- t(matrix(mod$estimate[(nn[1]+1) : nn[3]], nrow = k+1, byrow = TRUE))
    out$species <- specpara[,-1, drop=FALSE]
    out$b0 <- specpara[,1, drop=FALSE]
    out$fitted <- ginv(outer(-0.5*rowSums(out$points^2), drop(out$b0), "+") +
        out$points %*% t(out$species))
    out$spdev <- colSums(dev(y, out$fitted, wts))
    out$null.spdev <- null.spdev
    rownames(out$points) <- rownames(comm)
    rownames(out$species) <- colnames(comm)
    colnames(out$species) <- paste0("b", 1:k)
    out$family <- fam
    out$call <- match.call()
    class(out) <- "GO2"
    out
}

`print.GO2` <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat(gettextf(ngettext(x$k, "Gaussian Ordination with %d dimension\n",
                 "Gaussian Ordination with %d dimensions\n"), x$k))
    writeLines(strwrap(pasteCall(x$call)))
    cat("\n")
    cat(gettextf("%d iterations ", x$iterations))
    cat(switch(x$code,
               "(converged)",
               "(iterates within tolerance, probably converged)",
               "(step failed, perhaps a local minimum)",
               "(too many iterations)",
               "(convergence failed)"), "\n\n")
    devs <- c(x$null.deviance, x$null.deviance - x$deviance, x$deviance)
    props <- c(NA, devs[2:3]/devs[1])
    dfs <- c(x$df.null, x$df.null - x$df.residual, x$df.residual)
    table <- cbind(devs, props, dfs)
    rownames(table) <- c("Null", "Model", "Residual")
    colnames(table) <- c("Deviance", "Proportion","Df")
    cat("Family", x$family$family, "\n")
    printCoefmat(table, digits = digits, na.print="", zap.ind=c(1,2))
    invisible(x)
}

`plot.GO1` <-
    function(mod, ...)
{
    x <- mod$estimate
    tot <- mod$tot
    comm <- mod$data
    mods <- lapply(comm, function(y, ...)
                   glm(cbind(y, tot-y) ~ x + offset(-0.5 * x^2),
                   family = quasibinomial))
    grad <- seq(min(x), max(x), len=101)
    fit <- sapply(mods, function(z)
                  predict(z, newdata = list(x = grad), type = "response"))
    matplot(grad, fit, type="l", lty=1, ...)
    rug(x)
}

`plot.GO2` <-
    function(mod, choices = 1, label = FALSE, marginal = FALSE,
             cex=0.7, ...)
{
    x <- scores(mod, choices = choices, type="sites")
    ginv <- mod$family$linkinv
    grad <- seq(min(x), max(x), len=101)
    top <- drop(mod$b0)
    if (marginal)
        top <- top + 0.5 * rowSums(mod$species^2)
    fit <- ginv(outer(-0.5 * grad^2, top, "+") +
                outer(grad, mod$species[,choices], "*"))
    matplot(grad, fit, type = "l", lty=1, ...)
    rug(x)
    if (label) {
        br <- eval(formals(matplot)$col)
        xlab <- grad[apply(fit, 2, which.max)]
        ylab <- apply(fit, 2, max)
        ordilabel(cbind(xlab, ylab), cex = cex, border = br, col=1, ...)
    }
}


## Basic anova against Null model of constant response
`anova.GO2` <-
    function(object, ...)
{
    ## Check first if this should go to anova.GO2list
    dotargs <- list(...)
    if (length(dotargs)) {
        isGO2 <- sapply(dotargs, function(z) inherits(z, "GO2"))
        dotargs <- dotargs[isGO2]
        if (length(dotargs))
            return(anova.GO2list(c(list(object), dotargs)))
    }
    Df <- object$df.null - object$df.residual
    Dev <- object$null.deviance - object$deviance
    Fstat <- Dev/Df/(object$deviance/object$df.residual)
    pval <- pf(Fstat, Df, object$df.residual, lower.tail = FALSE)
    out <- data.frame(c(NA, Df),
                      c(NA, Dev),
                      c(object$df.null, object$df.residual),
                      c(object$null.deviance, object$deviance),
                      c(NA, Fstat),
                      c(NA, pval))
    colnames(out) <- c("Df", "Deviance", "Resid. Df", "Resid. Dev", "F", "Pr(>F)")
    rownames(out) <- c("NULL", "Model")
    class(out) <- c("anova", "data.frame")
    out
}

`anova.GO2list` <-
    function(object, ...)
{
    nmodels <- length(object)
    resdev <- sapply(object, deviance)
    resdf <- sapply(object, df.residual)
    df <- -diff(resdf)
    n <- sapply(object, function(z) nrow(z$points))
    table <- data.frame(resdf, resdev, c(NA, -diff(resdf)),
                        c(NA, -diff(resdev)))
    dimnames(table) <- list(1L:nmodels, c("Resid. Df", "Resid. Dev", 
                                          "Df", "Deviance"))
    big <- which.min(resdf)
    scale <- resdev[big]/resdf[big]
    table <- stat.anova(table, test="F", scale = scale,
                        df.scale = resdf[big], n[big])
    class(table) <- c("anova", "data.frame")
    table
}

## spdev functions analyse each species separately with F-test

`spanodev` <-
    function(mod1, mod2 = NULL, ...)
{
    if (is.null(mod2)) {
        dev1 <- mod1$null.spdev
        dev2 <- mod1$spdev
        df1 <- 1
        df2 <- mod1$k + 1
        dfr <- df.residual(mod1)/length(dev1)
        labs <- c("Null", "Model")
    } else {
        dev1 <- mod1$spdev
        dev2 <- mod2$spdev
        df1 <- mod1$k + 1
        df2 <- mod2$k + 1
        dfr <- min(df.residual(mod1), df.residual(mod2))/length(dev1)
        labs <- paste0("Model", 1:2)
    }
    ddev <- dev1 - dev2
    n <- nrow(mod1$fitted)
    ddf <- df2 - df1
    table <- data.frame(dev1, dev2, ddev)
    big <- which.max(c(df1,df2))
    scl <- table[,big]/dfr
    table$f <- table[,3]/ddf/scl
    table$p <- pf(table$f, ddf, dfr, lower.tail = FALSE)
    colnames(table) <- c(labs, "Change", "F", "Pr(>F)")
    head <- gettextf("F statistics based on (%d, %.1f) degrees of freedom", ddf, dfr)
    structure(table, heading = head, class = c("anova", "data.frame"))
}

## Function normalizes SU scores so that (unweighted) average
## tolerance = 1 for all species with unimodal responses and then
## estimates the glm parameters for normalized SU scores.
`GOnormalize` <-
    function(comm, p, k, tot, family, ...)
{
    ## extract SU scores from parameters p. When this function is
    ## evaluated for the first time, there are only these scores, but
    ## at later evaluation, there will be species parameters which
    ## will be ignored.
    n <- nrow(comm)
    x <- matrix(p[1: (k*n)], nrow=n, ncol=k)
    ## fit Gaussian response with free parameters
    b <- sapply(comm, function(y)
                coef(glm(cbind(y, tot-y) ~ x + I(x^2), family = family)))
    ## keep only negative 2nd degree coefficients
    b <- b[-(1:(k+1)),]
    b[b >= 0] <- NA
    scl <- rowMeans(sqrt(-1/2/b), na.rm=TRUE)
    x <- sweep(x, 2, scl, "/")
    ## fit Gaussian response with fixed tol=1
    b <- sapply(comm, function(y)
                coef(glm(cbind(y, tot-y) ~ x + offset(-0.5 * rowSums(x^2)),
                         family = family)))
    c(as.vector(x), as.vector(t(b)))
}
