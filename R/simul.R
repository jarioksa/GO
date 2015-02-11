#' Run a Set of Simulations and Collect Results
#'
#' Function is a general driver to run several ordination methods on a
#' set of simulated data sets and collect the goodness of fit
#' statistics for the results.
#'
#' @param nsimul Number of simulations.
#'
#' @param nsp Number of species in simulations passed to
#' \code{respfun}. The realized number can be (and usually is) lower
#' because some species do not occur in any simulated communities.
#'
#' @param n Number of sample plots passed to \code{respfun}.
#'
#' @param xgrad,ygrad Lengths of gradients passed to \code{respfun}.
#'
#' @param respfun The response function used in simulations.
#'
#' @param \dots Other arguments passed to \code{respfun}.
#'
#' @importFrom vegan specnumber metaMDS cca procrustes decorana
#'
#' @export
`simulrun` <-
    function(nsimul = 1, nsp=300, n=100, xgrad = 2, ygrad = 2, respfun = compas, ...)
{
    out <- matrix(0, nrow=nsimul, ncol=4)
    colnames(out) <- c("GO","NMDS","CA","DCA")
    for(i in 1:nsimul) {
        cat("====", i, "====\n")
        sim <- respfun(nsp = nsp, n = n, xgrad = xgrad, ygrad = ygrad, ...)
        comm <- as.data.frame(sim$comm)
        savecomm <<- comm
        print(summary(specnumber(comm)))
        cat("total (gamma):", ncol(comm), "\n")
        mgo <- GO(comm, k=2, family="binomial", tot=1, far=4, iterlim=1000)
        mmds <- metaMDS(comm, maxit=500, sratmax=0.999999)
        mca <- cca(comm)
        mdca <- decorana(comm)
        out[i,1] <- sqrt(procrustes(sim$x, mgo)$ss/n)
        out[i,2] <- sqrt(procrustes(sim$x, mmds)$ss/n)
        out[i,3] <- sqrt(procrustes(sim$x, mca, scaling=1)$ss/n)
        out[i,4] <- sqrt(procrustes(sim$x, mdca, choices=1:2)$ss/n)
        saveout <<- out
    }
    attr(out, "call") <- match.call()
    attr(out, "date") <- date()
    out
}
