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
    function(comm, tot = max(comm), freqlim = 5, ...)
{
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
    ## initialize gradient as first DCA axis
    x <- scores(decorana(comm), display = "sites", choices = 1)
    ## initial estimates for species (these are "final" with the
    ## current 'x')
    mods <- lapply(comm, function(y) glm(cbind(y, tot-y) ~ x + offset(-0.5 * x^2),
                                         family = fam))
    b <- sapply(mods, coef)
    ## Pack parameters to a single vector
    p <- c(x, b[1,], b[2,])
    nn <- cumsum(c(nrow(comm), ncol(comm), ncol(comm)))
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
        x <- p[1 : nn[1]]
        b0 <- p[(nn[1]+1) : nn[2]]
        b1 <- p[(nn[2]+1) : nn[3]]
        ## model lp = b0 -0.5*x^2 + b1*gr
        lp <- outer(-0.5*x^2, b0, "+") + outer(x, b1, "*")
        sum(dev(y, ginv(lp), wts))
    }
    out <- nlm(loss, p = p, ...)
    out$axis <- out$estimate[1 : nn[1]]
    out$species <- rbind(out$estimate[(nn[1]+1) : nn[2]],
                         out$estimate[(nn[2]+1) : nn[3]])
    names(out$axis) <- rownames(comm)
    colnames(out$species) <- colnames(comm)
    rownames(out$species) <- c("b0","b1")
    out$estimate <- NULL
    out$family <- fam
    class(out) <- "GO2"
    out
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
    function(mod, ...)
{
    x <- mod$axis
    ginv <- mod$family$linkinv
    grad <- seq(min(x), max(x), len=101)
    fit <- ginv(outer(-0.5 * grad^2, mod$species[1,], "+") +
                outer(grad, mod$species[2,], "*"))
    matplot(grad, fit, type = "l", lty=1, ...)
    rug(x)
}
