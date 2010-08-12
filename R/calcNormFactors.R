calcNormFactors <- function(object, refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10) {
  if( is.matrix(object) ) {
    if(is.null(refColumn))
      refColumn <- 1
    apply(object,2,.calcFactorWeighted,ref=object[,refColumn], logratioTrim=logratioTrim, 
        sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff)
  } else if(is(object, "DGEList")) {
    D <- object$counts
    if(is.null(refColumn)) {
      
      y <- t(t(D)/object$samples$lib.size)
      q75 <- apply(y,2,function(x) quantile(x,p=0.75))
      refColumn <- which.min(abs(q75-mean(q75)))
    }

    f <- apply(D,2,.calcFactorWeighted,ref=D[,refColumn], logratioTrim=logratioTrim, 
        sumTrim=sumTrim, doWeighting=doWeighting, Acutoff=Acutoff)
    object$samples$norm.factors <- f/exp(mean(log(f)))
    object$ref.column <- refColumn
    
    object
  } else {
    stop("calcNormFactors() only operates on 'matrix' and 'DGEList' objects")
  }
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
  
