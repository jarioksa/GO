`compas` <-
    function (nsp, n, xgrad=6, ygrad=4, x, xn=6, yn=4, buffer=2, tsd = 0.1) 
{
    gauss <- function(x,p) { 
        p[1]*exp(-(x[1]-p[2])^2/2/p[4]^2 - (x[2]-p[3])^2/2/p[5]^2) 
    }
    p <- matrix(0, nrow=nsp, ncol=5)
    ## Integral of bivariate Gaussian is 2*pi*h*t. Uniform h gives
    ## average h=0.5 and with t=1 this gives pi. The average of
    ## log-Normal t is exp(tsd^2/2). To get desired alpha per SU, use
    ## nsp = alpha * (xgrad + buffer) * (ygrad + buffer) / pi /
    ## exp(tsd^2/2).
    p[,1] <- runif(nsp)
    p[,2] <- runif(nsp, 0-buffer, xgrad+buffer)
    p[,3] <- runif(nsp, 0-buffer, ygrad+buffer)
    p[,4] <- rlnorm(nsp, 0, tsd)
    p[,5] <- rlnorm(nsp, 0, tsd)
    if (missing(x)) {
        if(missing(n)) {
            x1 <- seq(0, xgrad, len=xn)
            x2 <- seq(0, ygrad, len=yn)
            x <- as.matrix(expand.grid(x1,x2))
            n <- xn*yn
        } else {
            x1 <- runif(n, 0, xgrad)
            x2 <- runif(n, 0, ygrad)
            x <- cbind(x1, x2)
            xn <- NA
            yn <- NA
        }
    } else {
        n <- nrow(x)
        xn <- yn <- NA
    }
    comm <- matrix(0, nrow=nrow(x), ncol=nsp)
    for (i in 1:nrow(x)) for (j in 1:nsp) comm[i,j] <- gauss(x[i,],p[j,])
    cdim <- dim(comm)
    comm <- rbinom(length(comm), 1, comm)
    dim(comm) <- cdim
    freq <- apply(comm, 2, sum)
    comm <- comm[,freq>0]
    p <- p[freq>0,]
    sol <- list(comm=comm, x=x, species=p, n=n, nx=xn, ny=yn,
                call = match.call())
    class(sol) <- "compas"
    sol
}

`betapas` <-
    function (nsp, n, xgrad=6, ygrad=4, x, xn=6, yn=4, buffer=2, tsd=0.1,
              range=5, shape=c(0.5,6.5))
{
    ## Beta response for one gradient: 
    betaresp <- function(x, b) {
        range <- b[3] - b[2]
        z <- (x - b[2])/range
        mu <- ifelse(z <= 0 | z >= 1, 0, b[1]*z^b[4]*(1-z)^b[5])
        mu <- ifelse(mu > 1, 1, mu)
        mu
    }
    ## Define first gaussian parameters, then transform these to beta
    ## response parameters with shape parameters (alpha, gamma) are
    ## runif in the range of 'shape', endpoints are opt +/-
    ## 0.5*t*range, and the height adjustment k is defined so that the
    ## integral of the betaresponse equals the integral of Gaussian
    p <- matrix(0, nrow=nsp, ncol=5)
    p[,1] <- exp(runif(nsp, -2.5,0))
    p[,2] <- runif(nsp, 0-buffer, xgrad+buffer)
    p[,3] <- runif(nsp, 0-buffer, ygrad+buffer)
    p[,4] <- rlnorm(nsp, 0, tsd)
    p[,5] <- rlnorm(nsp, 0, tsd)
    ## Transfer these to beta parameters separately for both gradients
    b1 <- matrix(0, nrow=nsp, ncol=5)
    b2 <- matrix(0, nrow=nsp, ncol=5)
    colnames(b1) <- colnames(b2) <- c("k","p1","p2","alpha","gamma")
    b1[,2] <- p[,2] - range * p[,4]/2
    b1[,3] <- p[,2] + range * p[,4]/2
    b1[,4] <- runif(nsp, shape[1], shape[2])
    b1[,5] <- runif(nsp, shape[1], shape[2])
    b1[,1] <- sqrt(p[,1]) * sqrt(2*pi) / beta(b1[,4]+1, b1[,5]+1) / range
    ## second gradient
    b2[,2] <- p[,3] - range * p[,5]/2
    b2[,3] <- p[,3] + range * p[,5]/2
    b2[,4] <- runif(nsp, shape[1], shape[2])
    b2[,5] <- runif(nsp, shape[1], shape[2])
    b2[,1] <- sqrt(p[,1]) * sqrt(2*pi) / beta(b2[,4]+1, b2[,5]+1) / range
    ## Rest copied from compas: sampiling points
    if (missing(x)) {
        if(missing(n)) {
            x1 <- seq(0, xgrad, len=xn)
            x2 <- seq(0, ygrad, len=yn)
            x <- as.matrix(expand.grid(x1,x2))
            n <- xn*yn
        } else {
            x1 <- runif(n, 0, xgrad)
            x2 <- runif(n, 0, ygrad)
            x <- cbind(x1, x2)
            xn <- NA
            yn <- NA
        }
    } else {
        n <- nrow(x)
        xn <- yn <- NA
    }
    comm <- matrix(0, nrow=nrow(x), ncol=nsp)
    for (i in 1:nrow(x))
        for (j in 1:nsp)
            comm[i,j] <- (betaresp(x[i,1],b1[j,])) * (betaresp(x[i,2], b2[j,]))
    cdim <- dim(comm)
    comm <- rbinom(length(comm), 1, comm)
    dim(comm) <- cdim
    freq <- apply(comm, 2, sum)
    comm <- comm[,freq>0]
    b1 <- b1[freq>0,]
    b2 <- b2[freq>0,]
    sol <- list(comm=comm, x=x, species=cbind(b1,b2), n=n, nx=xn, ny=yn,
                call = match.call())
    class(sol) <- "compas"
    sol
}
