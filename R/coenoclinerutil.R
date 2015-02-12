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

    ## uniform h in (0,1)
    h <- runif(nsp)
    ## uniform opt in range +/- buffer
    ux <- runif(nsp, -xrange/2 - buffer, xrange/2 + buffer)
    uy <- runif(nsp, -yrange/2 - buffer, yrange/2 + buffer)
    ## lognormal tol. Mean of lognormal is exp(tsd^2/2)
    tx <- rlnorm(nsp, -tsd^2/2, tsd)
    ty <- rlnorm(nsp, -tsd^2/2, tsd)
    list(px = cbind("opt" = ux, "tol" = tx, "h" = h),
         py = cbind("opt" = uy, "tol" = ty))
}

#' @param gausspar Gaussian response parameters for species as
#' returned by \code{BinomGaussPar}.
#' @param shape Random log-uniform range of shape parameters \eqn{alpha}
#' and \eqn{gamma} of response function
#' @param range Range of beta response in $t$ (\sQuote{sd}) units of
#' Gaussian response function.
#'
#' @describeIn coenoclinerutil Translate Gaussian parameters into
#' corresponding beta response parameters.
#' @export
`Gauss2betaPar` <-
    function(gausspar, shape = c(0.5, 6.5), range = 5)
{
    ## Define beta response so that it as similar to a Gaussian model
    ## as possible -- except for shape. Input **must** be similar as
    ## from BinomGausPar. This is not checked.
    nsp <- nrow(gausspar[[1]])
    ## shapes uniform in (shape)
    shape <- log(shape)
    ax <- exp(runif(nsp, shape[1], shape[2]))
    gx <- exp(runif(nsp, shape[1], shape[2]))
    ay <- exp(runif(nsp, shape[1], shape[2]))
    gy <- exp(runif(nsp, shape[1], shape[2]))
    ## ranges a multiple of Gaussian tol
    rx <- range * gausspar$px[,"tol"]
    ry <- range * gausspar$py[,"tol"]
    ## modal abundance at Gaussian opt
    mx <- gausspar$px[,"opt"]
    my <- gausspar$py[,"opt"]
    ## Response height A0 should be such that beta response has the
    ## same mass as the corresponding Gaussian. The integral of
    ## univariate Gaussian response is h*t*sqrt(2*pi) and the integral
    ## of beta response is adj*range*beta(alpha+1, gamma+1), and we
    ## need to find A0 giving the desired height adjustment adj, and
    ## here beta() is the real mathematical beta function. However, we
    ## do not want A0>1 because we target Binomial models.
    Gmass <- with(gausspar, px[,"h"] * px[,"tol"] * py[,"tol"] * 2 * pi)
    Bmass <- rx * ry * beta(ax+1, gx+1) * beta(ay+1, gy+1)
    adj <- Gmass/Bmass
    ## bx, by and A0 are from Minchin, Vegetatio 71, 145-156 (1987),
    ## and they are similarly used in coenocliner.
    bx <- ax/(ax+gx)
    by <- ay/(ay+gy)
    A0 <- pmin(adj * bx^ax * (1-bx)^gx * by^ay * (1-by)^gy, 1)
    ## collect
    list(px = cbind("m" = mx, "r" = rx, "alpha" = ax, "gamma" = gx, "A0" = A0),
         py = cbind("m" = my, "r" = ry, "alpha" = ay, "gamma" = gy))
}

#' @param comm Community data.
#'
#' @describeIn coenoclinerutil Drop missing species from the data.
#' @export
`DropMissingSpec` <-
    function(comm)
{
    colnames(comm) <- paste0("sp", seq_along(ncol(comm)))
    comm[, colSums(comm) > 0]
}
