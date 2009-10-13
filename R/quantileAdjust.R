quantileAdjust<-function(object, N = prod(object$samples$lib.size)^(1/ncol(object$counts)),alpha = 0, null.hypothesis = FALSE, n.iter = 5, r.init = NULL, tol = 0.001, verbose = TRUE) 
# Written by Mark Robinson, edited by Davis McCarthy, February 2009, substantially revised July 2009
# A function to adjust counts for the estimation of common and tagwise dispersion(common dispersion only?)
# Returns estimate(s) of dispersion, pseudocounts, mean counts for tags, proportions for tags in samples and adjusted library size
{
    if (is.null(r.init)) {
        r.init <- 1/findMaxD2(object, alpha = 10) - 1
    }
    nrows <- nrow(object$counts)
    lib.size <- object$samples$lib.size
    group <- object$samples$group
    levs.group <- levels(group)
    y <- splitIntoGroups(object)
    p <- matrix(0, nrow = nrows, ncol = ncol(object$counts))
    mu <- matrix(0, nrow = nrows, ncol = ncol(object$counts))
    r <- r.init
    rprev <- r + 1
    count <- 0
    while (count < n.iter) {
        count <- count + 1
        if (verbose) 
            cat("[quantileAdjust] Iteration (dot=1000)", count, 
                ":")
        rprev <- r
        conc <- estimatePs(object, r)
        if (null.hypothesis == TRUE) {
            for (i in 1:length(levs.group)) {
                p[, group == levs.group[i]] <- pnbinom(y[[i]]-1, size = r, mu = outer(conc$conc.common, lib.size[group==levs.group[i]])) + dnbinom(y[[i]], size = r,mu = outer(conc$conc.common, lib.size[group == levs.group[i]]))/2
                mu[, group == levs.group[i]] <- outer(conc$conc.common,rep(N, sum(group == levs.group[i])))
            }
        }
        else {
            for (i in 1:length(levs.group)) {
                p[, group == levs.group[i]] <- pnbinom(y[[i]]-1, size = r, mu = outer(conc$conc.group[, i], lib.size[group==levs.group[i]])) + dnbinom(y[[i]], size = r, mu = outer(conc$conc.group[, i], lib.size[group==levs.group[i]]))/2
                mu[, group == levs.group[i]] <- outer(conc$conc.group[,i], rep(N, sum(group == levs.group[i])))
            }
        }
        count.max <- apply(object$counts, 1, max)
        pseudo <- interpolateHelper(mu, p, r, count.max, verbose = verbose)
        pseudo[pseudo < 0] <- 0
        d <- new("DGEList",list(samples=object$samples, counts=pseudo))
        r <- 1/findMaxD2(d, alpha = alpha) - 1
        if (max(abs(rprev - r)) < tol) {
            break
        }
    }
    return(list(r = r, pseudo = pseudo, p = p, mu = mu, conc = conc, N = N))
}
