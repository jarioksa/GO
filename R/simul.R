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
        mgo <- GO2(comm, k=2, family="binomial", tot=1, far=4, iterlim=1000)
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
