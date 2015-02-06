#' Unconstrained Gaussian Maximum Likelihood Ordination
#'
#' The functions fit unconstrained maximum likelihood ordination with
#' unit-width Gaussian response models.
#'
#' @details Function is under development and unreleased. It will be released
#'  under different name in \pkg{vegan}. The current version is only
#'  provided for review purposes. The function and its support functions
#'  require \pkg{vegan}, although this requirements is not explicitly
#'  made. The optimization is based on \code{\link{nlm}} function, and passes
#'  arguments to this function.
#'
#'  Function \code{anova} can analyse two nested models or a single model
#'  against null model of flat responses using parametric tests based on
#'  quasi-Likelihood. Function \code{spanodev} performs similar test 
#'  split by species. Function \code{predict} returns estimated response
#'  curves, and \code{newdata} can be gradient locations. Function 
#'  \code{calibrate} returns estimated gradient locations, and \code{newdata}
#'  can be community data. 
#'
#'  The \code{plot} function displays fitted respose curves against one
#'  ordination axis. In principle, the ordination can be rotated using 
#'  \pkg{vegan} function \code{\link[vegan]{MDSrotate}}, but this requires
#'  a version that agrees to analyse \code{GO2} results. Traditional 
#'  ordination plots of SU scores and species optima can be displayed
#'  with \code{\link[vegan]{ordiplot}} (\pkg{vegan} package). The function
#'  is written so that several other \pkg{vegan} and standard \R functions 
#'  can handle results.
#'
#' @seealso \code{\link[VGAM]{cgo}} in \pkg{VGAM} package.
#'
#' @examples
#' library(vegan) ## *must* be given before using the function
#' data(varespec)
#' mod <- GO2(varespec, k=2, far=5, tot=100, family="binomial", iterlim=1000)
#' plot(mod, label=TRUE)
#' ordiplot(mod, type="t")
#' ordiplot(mod, dis="si", type="t")
#' anova(mod)
#' mod1 <- update(mod, k=1)
#' anova(mod1, mod)
#' spanodev(mod1)
#' spanodev(mod1, mod)


## unconstrained Gaussian ordination in one dimension

#' @importFrom parallel makeCluster stopCluster
#' 
#' @param comm Community data frame.
#'
#' @param tot Total abundance used in Binomial models. This can be
#' either a single value for all data, or a vector with value for each
#' row of \code{comm}. The default is to use the maximum value in
#' matrix.
#'
#' @param freqlim Minimum number of occurrence for analysed species.
#'
#' @param parallel Number of parallel processes.
#'
#' @param \dots Other parameters passed to functions. In \code{GO2}
#' these are passed to \code{\link{nlm}} and can include, e.g.,
#' \code{iterlim} (which often must be set to higher value than the
#' default 100).
#'
#' @describeIn GO Alternating estimation of species parameters and
#' gradient locations in one dimension.
#' @export
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
    y <- as.matrix(comm)
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
    ## partial derivatives x
    grad <- function(y, x, tot, mods, fam = quasibinomial()) {
        eta <- sapply(mods, predict, type="link")
        mu <- fam$linkinv(eta)
        b <- sapply(mods, coef)[2,]
        -rowSums(tot * (y-mu) / fam$var(mu) * fam$mu.eta(eta) *
            outer(-x, b, "+"))
    }
    if (parallel > 1)
        loss <- function(x, ...) {
            mods <- parLapply(clu, comm, 
                              function(y, ...)
                              glm(cbind(y, tot-y) ~ x + offset(-0.5 * x^2),
                                  family = quasibinomial))
            ll <- sum(sapply(mods, deviance))/2
            attr(ll, "gradient") <- grad(y/tot, x, tot, mods)
            ll
        }
    else
        loss <- function(x, ...) {
            mods <- lapply(comm, function(y, ...)
                           glm(cbind(y, tot-y) ~ x + offset(-0.5 * x^2),
                               family = quasibinomial))
            ll <- sum(sapply(mods, deviance))/2
            attr(ll, "gradient") <- grad(y/tot, x, tot, mods)
            ll
        }
    ## ML fit
    out <- nlm(loss, p = x, comm = comm, tot = tot, ...)
    out$minimum <- 2 * out$minimum
    out$data <- as.matrix(comm)
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

#' @importFrom vegan decorana
#' 
#' @param k Number of estimated gradients (axes).
#'
#' @param family Error distribution. Can be either \code{"poisson"}
#' for quasi-Poisson or \code{"binomial"} for quasi-Binomial (and must
#' be quoted).
#'
#' @param far Threshold distance for species optima regarded as
#' alien and frozen in fitting.
#'
#' @describeIn GO Simultaneous estimation of species parameters and
#' gradient locations.
#' @export
GO2 <-
    function(comm, k = 1, tot = max(comm), freqlim = 5,
             family = c("poisson", "binomial"), far = 10, ...)
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
    family <- match.arg(family)
    fam <- switch(family,
                  "poisson" = quasipoisson(),
                  "binomial" = quasibinomial(),
                  "gaussian" = gaussian(link=log))
    ginv <- fam$linkinv
    dev <- fam$dev.resids
    Var <- fam$var
    mu.eta <- fam$mu.eta
    if (family != "binomial")
        tot <- 1
    ## initialize gradient as DCA axes
    x <- as.matrix(scores(decorana(comm), display = "sites", choices = 1:k))
    ## See if there are weights
    wts <- matrix(rep(tot, ncol(comm)), nrow(comm), ncol(comm))
    rwts <- wts[,1]
    comm <- comm / wts
    y <- as.matrix(comm)
    ## initial estimates for species (these are "final" with the
    ## current 'x')
    mods <- lapply(comm, function(y) glm(y ~
                                         x + offset(-0.5 * rowSums(x^2)),
                                         family = fam, weights=rwts))
    b <- sapply(mods, coef)
    null.spdev <- sapply(lapply(comm, function(y) glm(y ~ 1,
           family= fam, weights=rwts)), function(z) z$deviance)
    null.deviance <- sum(null.spdev)
    ## Pack parameters to a single vector
    p <- c(as.vector(x), as.vector(t(b)))
    nn <- cumsum(c(k*nrow(comm), ncol(comm), k*ncol(comm)))
    ## We need matrix response, and with quasibinomial, we need to
    ## divide with 'tot'
    if (length(tot) < nrow(comm))
        tot <- rep(tot, length.out = nrow(comm))
    ## Loss function fits the current model and return the sum of
    ## deviances
    alien <- logical(k * ncol(comm))
    loss <- function(p, ...) {
        ## split params
        x <- matrix(p[1 : nn[1]], ncol=k)
        b0 <- p[(nn[1]+1) : nn[2]]
        b1 <- matrix(p[(nn[2]+1) : nn[3]], nrow = k, byrow=TRUE)
        if (far) {
            cnt <- colMeans(x)
            diam <- sqrt(max(rowSums(sweep(x, 2, cnt)^2)))
            alien <- sqrt(colSums(sweep(b1, 1, cnt)^2)) > far + diam
        }
        ## model lp = b0 -0.5*x^2 + b1*x
        lp <- outer(-0.5*rowSums(x^2), b0, "+") + x %*% b1
        mu <- ginv(lp)
        ## deviance/2 = -log-likelihood + constant
        ll <- sum(dev(y, mu, wts))/2
        ## Derivatives are based on McCullagh & Nelder 1989, p. 41
        ## (eq. 2.13, and unnumbered equation on the same page)
        attr(ll, "gradient") <- {
            .der <- wts * (y - mu) /Var(mu) * mu.eta(lp)
            .ader <- colSums(.der)
            .bder <- t(.der) %*% x
            if (any(alien))
                .bder[alien,] <- 0
            .xder <- sapply(seq_len(k), function(dim)
                            rowSums(.der * outer(-x[,dim], b1[dim,], "+")))
            ## Combine and reverse sign: we miminize instead of maximizing
            .value <- -drop(c(as.vector(.xder), .ader, as.vector(.bder)))
            ##plot(.value, pch=".", ylim=c(-30,30))
            .value
        }
        ll
    }
    mod <- nlm(loss, p = p, ...)
    out <- list(deviance = 2*mod$minimum, null.deviance = null.deviance,
                k = k, iterations = mod$iterations, code = mod$code,
                rank = length(mod$gradient),
                df.residual = prod(dim(comm)) - length(mod$gradient),
                df.null = prod(dim(comm)) - ncol(comm), gradient = mod$gradient)
    out$points <- matrix(mod$estimate[1 : nn[1]], ncol=k)
    specpara <- t(matrix(mod$estimate[(nn[1]+1) : nn[3]], nrow = k+1, byrow = TRUE))
    out$species <- specpara[,-1, drop=FALSE]
    out$b0 <- specpara[,1, drop=FALSE]
    ## centre data
    cnt <- colMeans(out$points)
    ## height para needs adjustment for moving species scores
    out$b0 <- out$b0 + out$species %*% cnt - 0.5 * sum(cnt^2)
    out$points <- sweep(out$points, 2, cnt)
    out$species <- sweep(out$species, 2, cnt)
    ## rotate to PCs
    pc <- prcomp(out$points)
    out$points <- pc$x
    out$species <- out$species %*% pc$rotation
    ## recreate fitted data
    out$fitted <- ginv(outer(-0.5*rowSums(out$points^2), drop(out$b0), "+") +
        out$points %*% t(out$species))
    ## observed data (as analysed in the model)
    out$y <- y
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

#' @export
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
               "(perhaps a local minimum)",
               "(iteration limit reached)",
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

#' @param mod Fitted model.
#'
#' @rdname GO
#' @export
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

#' @importFrom vegan ordilabel scores
#' 
#' @param choices The axis or axes plotted.
#'
#' @param label Label species responses.
#'
#' @param marginal Plot marginal responses or slice through origin of
#' other dimensions.
#'
#' @param cex Character size for \code{labels}.
#'
#' @param col Colours of response curves.
#'
#' @rdname GO
#' @export
`plot.GO2` <-
    function(mod, choices = 1, label = FALSE, marginal = FALSE,
             cex=0.7, col = 1:6, ...)
{
    x <- scores(mod, choices = choices, type="sites")
    ginv <- mod$family$linkinv
    grad <- seq(min(x), max(x), len=101)
    top <- drop(mod$b0)
    if (marginal)
        top <- top + 0.5 * rowSums(mod$species^2)
    fit <- ginv(outer(-0.5 * grad^2, top, "+") +
                outer(grad, mod$species[,choices], "*"))
    matplot(grad, fit, type = "l", lty=1, col=col, ...)
    rug(x)
    if (label) {
        br <- col
        xlab <- grad[apply(fit, 2, which.max)]
        ylab <- apply(fit, 2, max)
        ordilabel(cbind(xlab, ylab), cex = cex, border = br, col=1, ...)
    }
}


## Basic anova against Null model of constant response

#' @param object Ordination result object.
#'
#' @rdname GO
#' @export
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

#' @export
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

#' @param mod1,mod2 Compared result objects
#'
#' @describeIn GO Comparison of goodness of fit for individual species.
#' @export
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

#' @param newdata New gradient locations for predicted responses
#'
#' @param type Predictions in the scale of responses or in the scale
#' of link function.
#'
#' @rdname GO
#' @export
`predict.GO2` <-
    function(object, newdata, type = c("response", "link"), ...)
{
    type <- match.arg(type)
    if (missing(newdata))
        x <- object$points
    else {
        x <- newdata
        if (length(dim(x)) == 2) {
            if (ncol(x) != object$k)
                stop(gettextf("number of columns in 'newdata' should be %d, was %d",
                              object$k, ncol(x)))
        } else {
            if (length(x) != object$k)
                stop(gettextf("number of items in 'newdata' should be %d, was %d",
                              object$k, length(x)))
            x <- matrix(x, ncol = object$k)
        }
    }
    a <- drop(object$b0)
    b <- object$species
    eta <-  outer(-0.5 * rowSums(x^2), a, "+") + x %*% t(b)
    if (type == "response")
        object$family$linkinv(eta)
    else
        eta
}

### The calibrate function finds the most likely gradient position
### (location of 'points') given the observed community composition
### and the fitted GO2 model.

#' @importFrom vegan calibrate
#' @rdname GO
#' @export
`calibrate.GO2` <-
    function(object, newdata, ...)
{
    if(!missing(newdata)) { 
        data <- newdata
        data <- data[, colnames(object$y), drop=FALSE]
        data <- as.matrix(data)
    }
    else
        data <- object$y
    ginv <- object$family$linkinv
    dev <- object$family$dev.resids
    mu.eta <- object$family$mu.eta
    V <- object$family$var
    b <- object$species
    ## loss function evaluates for a single SU
    loss <- function(p, y, ...) {
        eta <- predict(object, newdata = p, type = "link")
        mu <- ginv(eta)
        ll <- sum(dev(y, mu, rep(1,length(y))))/2
        attr(ll, "gradient") <-
            -((y-mu)/V(mu)*mu.eta(eta)) %*% sweep(b, 2, p)
        ll
    }
    ## initial estimates as WA scores of species optima
    bb <- b
    ## sanitize species scores: no extrapolation
    xx <- apply(object$points, 2, range)
    bb[,1] <- ifelse(bb[,1] < xx[1,1], xx[1,1], bb[,1])
    bb[,1] <- ifelse(bb[,1] > xx[2,1], xx[2,1], bb[,1])
    bb[,2] <- ifelse(bb[,2] < xx[1,2], xx[1,2], bb[,2])
    bb[,2] <- ifelse(bb[,2] > xx[2,2], xx[2,2], bb[,2])
    xmod <- wascores(bb, t(data))
    ##xmod <- wascores(b, t(data))
    ## Initialize as zero
    ## xmod <- matrix(0, nrow=nrow(data), ncol=object$k)
    xcal <- lapply(seq_len(NROW(data)),
                   function(i, ...) nlm(loss, p = xmod[i,], y = data[i,], hessian = FALSE,...))
    t(sapply(xcal, function(z) z$estimate))
    #xcal
}

