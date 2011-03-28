betaApproxNBTest <- function(x1, x2, dispersion) {
    ## Uses a transformed Beta distribution to approximate the tail probabilities of a conditional negative binomial exact test of equality of means between two groups
    ## Davis McCarthy.
    ## October 2010.
    x.obs <- pmax(x1, x2)
    z <- x1 + x2
    zeta.obs <- (x.obs - 0.5) / z
    alpha <- beta <- (z - 1) / (2 + dispersion*z)
    pvalue <- 2*pbeta(zeta.obs, shape1=alpha, shape2=beta, lower.tail=FALSE)
    pvalue
}




