calcNormFactors <- function(object, method=c("TMM","RLE","quantile"), refColumn=NULL,
                            logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, 
                            quantile=0.75) {
                            
  method <- match.arg(method)
                            
  if( is.matrix(object) ) {
    if(is.null(refColumn))
      refColumn <- 1
    data <- object
    libsize <- colSums(data)
  } else if(is(object, "DGEList")) {
    data <- object$counts
    if(method=="TMM" & is.null(refColumn)) {      
      f75 <- .calcFactorQuantile(data=data, lib.size=object$samples$lib.size, q=0.75)
      refColumn <- which.min(abs(f75-mean(f75)))
    }
    libsize <- object$samples$lib.size
  } else {
    stop("calcNormFactors() only operates on 'matrix' and 'DGEList' objects")
  }
  

  f <- switch(method,
              TMM = apply(data,2,.calcFactorWeighted,ref=data[,refColumn], 
                          logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, 
                          Acutoff=Acutoff),
              RLE = .calcFactorRLE(data)/libsize,
              quantile = .calcFactorQuantile(data, libsize, q=quantile))
              
  f <- f/exp(mean(log(f)))

  if( is.matrix(object) ) {
    return(f)
  } else if(is(object, "DGEList")) {
    object$samples$norm.factors <- f
    return(object)
  }

  
}


.calcFactorRLE <- function (data) {
    gm <- exp(rowMeans(log(data)))
    apply(data, 2, function(u) median((u/gm)[gm > 0]))
}

.calcFactorQuantile <- function (data, lib.size, q=0.75) {
    y <- t(t(data)/lib.size)
    f <- apply(y,2,function(x) quantile(x,p=q))
    f/exp(mean(log(f)))
}


.calcFactorWeighted <- function(obs, ref, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10) {

  if( all(obs==ref) )
    return(1)

  nO <- sum(obs)
  nR <- sum(ref)
  logR <- log2((obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size
  absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression
  v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref   # estimated asymptotic variance
  
  # remove infinite values, cutoff based on A
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]
  
  # taken from the original mean() function
  n <- sum(fin)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS
  
  #keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
  # a fix from leonardo ivan almonacid cardenas, since rank() can return
  # non-integer values when there are a lot of ties
  keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)
  
  if (doWeighting) 
    2^( sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE) )
  else
    2^( mean(logR[keep], na.rm=TRUE) )
}
  
