DGEList <- function(counts=matrix(0,0,0), lib.size=colSums(counts), norm.factors=rep(1,ncol(counts)), group=rep(1,ncol(counts)), genes=NULL, remove.zeros=FALSE) 
#	Construct DGEList object from components, with some checking
#	Last modified 7 April June 2013
{
#	Check counts
	counts <- as.matrix(counts)
	nlib <- ncol(counts)
	ntags <- nrow(counts)
	if(nlib>0 && is.null(colnames(counts))) colnames(counts) <- paste("Sample",1:ncol(counts),sep="")
	if(ntags>0 && is.null(rownames(counts))) rownames(counts) <- 1:ntags

#	Check lib.size
	if(is.null(lib.size)) lib.size <- colSums(counts)
	if(nlib != length(lib.size)) stop("Length of 'lib.size' must equal number of columns in 'counts'")

#	Check norm.factors
	if(is.null(norm.factors)) norm.factors <- rep(1,ncol(counts))
	if(nlib != length(norm.factors)) stop("Length of 'norm.factors' must equal number of columns in 'counts'")

#	Check group
	if(is.null(group)) group <- rep(1,ncol(counts))
	group <- as.factor(group)
	if(nlib != length(group)) stop("Length of 'group' must equal number of columns in 'counts'")

	samples <- data.frame(group=group,lib.size=lib.size,norm.factors=norm.factors)
	row.names(samples) <- colnames(counts)
	x <- new("DGEList",list(counts=counts,samples=samples))

	if(!is.null(genes)) {
		genes <- as.data.frame(genes, stringsAsFactors=FALSE)
		if(nrow(genes) != ntags) stop("counts and genes have different numbers of rows")
		x$genes <- genes
	}

	if(remove.zeros) {
		all.zeros <- rowSums(counts,na.rm=TRUE)==0
		if(any(all.zeros)) {
			x <- x[!all.zeros,]
			message("Removing ",sum(all.zeros)," rows with all zero counts.")
		}
	}

	x
}
