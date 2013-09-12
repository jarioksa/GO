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
}
