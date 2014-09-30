.residDF <- function(zero, design)
#	Effective residual degrees of freedom after adjusting for exact zeros
#	Gordon Smyth and Aaron Lun
#	Created 6 Jan 2014.  Last modified 2 Sep 2014
{
	nlib <- ncol(zero)
	ncoef <- ncol(design)
	nzero <- as.integer(rowSums(zero))

#	Default is no zero
	DF <- rep(nlib-ncoef,length(nzero))

#	All zero case
	DF[nzero==nlib] <- 0L

#	Anything in between?
	somezero <- nzero>0L & nzero<nlib
	if(any(somezero)) {
		zero2 <- zero[somezero,,drop=FALSE]

#		Integer packing will only work for 31 libraries at a time.	
		assembly <- list()	
		collected <- 0L
		step <- 31L
		bits <- as.integer(2^(1:step-1L))
		while (collected < ncol(zero2)) {
			upper <- pmin(ncol(zero2) - collected, step)
			keys <- t(zero2[,collected+1:upper,drop=FALSE]) * bits[1:upper]
			assembly[[length(assembly)+1L]] <- as.integer(colSums(keys))
			collected <- collected + step
		}

#		Figuring out the unique components.
		o <- do.call(order, assembly)
		nzeros <- sum(somezero)
		is.different <- logical(nzeros)
		for (i in 1:length(assembly)) { 
			is.different <- is.different | c(TRUE, diff(assembly[[i]][o])!=0L)
		}
		first.of.each <- which(is.different)
		last.of.each <- c(first.of.each[-1]-1L, nzeros)

#		Identifying the true residual d.f. for each of these rows.			
		DF2 <- nlib-nzero[somezero]
		for (u in 1:length(first.of.each)) {
			i <- o[first.of.each[u]:last.of.each[u]]
			zeroi <- zero2[i[1],]
			DF2[i] <- DF2[i]-qr(design[!zeroi,,drop=FALSE])$rank
		}
		DF2 <- pmax(DF2, 0L)
		DF[somezero] <- DF2
	}
	DF
}
