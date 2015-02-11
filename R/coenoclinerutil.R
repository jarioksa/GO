#' Utility Functions for 'coenocliner' Package
#'
#' Functions to automated simulation routines using \pkg{coenocliner}
#' package.
#'
#' @param n Number of SUs.
#' @param xrange,yrange Desired range of gradients.
#'
#' @describeIn coenoclinerutil Gradient Locations
#' @export
`GradLocs` <-
    function(n, xrange, yrange)
{
    cbind("x" = runif(n, -xrange/2, xrange/2),
          "y" = runif(n, -yrange/2, yrange/2))
}
#' 
#' @param nsp Number of species.
#' @param buffer Width of buffer zone for optima surrounding ranges.
#' @param tsd Standard deviation of tolerance in log-Normal
#' distribution, in log scale
#'
#' @describeIn coenoclinerutil Gaussian Parameters for Binomial Response.
#' @export
`BinomGaussPar` <-
    function(nsp, xrange, yrange, buffer=2, tsd=0.1)
{
    ## Create Gaussian parameters for Binomial responses to be used
    ## with coenocliner
    h <- runif(nsp)
    ux <- runif(nsp, -xrange/2 - buffer, xrange/2 + buffer)
    uy <- runif(nsp, -yrange/2 - buffer, yrange/2 + buffer)
    tx <- rlnorm(nsp, 0, tsd)
    ty <- rlnorm(nsp, 0, tsd)
    list(px = cbind("opt" = ux, "tol" = tx, "h" = h),
         py = cbind("opt" = uy, "tol" = ty))
}

#' @param comm Community data.
#'
#' @describeIn coenoclinerutil Drop missing species from the data.
#' @export
`DropMissingSpec` <-
    function(comm)
{
    comm[, colSums(comm) > 0]
}
    

