.residDF <- function(zero, design)
#	6 Jan 2014
{
	nlib <- ncol(zero)
	ncoef <- ncol(design)
	nzero <- rowSums(zero)

#	Default is no zero
	DF <- rep(nlib-ncoef,length(nzero))

#	All zero case
	DF[nzero==nlib] <- 0

#	Anything in between?
	somezero <- nzero>0 & nzero<nlib
	if(any(somezero)) {
		zero2 <- zero[somezero,,drop=FALSE]	
		key <- rowSums( 2L^(col(zero2)-1L) * zero2 )
		DF2 <- nlib-nzero[somezero]
		for (u in unique(key)) {
			i <- which(key==u)
			zeroi <- zero2[i[1],]
			DF2[i] <- DF2[i]-qr(design[!zeroi,,drop=FALSE])$rank
		}
		DF2 <- pmax(DF2,0)
		DF[somezero] <- DF2
	}
	DF
}
