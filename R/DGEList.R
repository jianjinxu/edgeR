DGEList <- function(counts=matrix(0,0,0), lib.size=colSums(counts), norm.factors=rep(1,ncol(counts)), samples=NULL, group=NULL, genes=NULL, remove.zeros=FALSE) 
#	Construct DGEList object from components, with some checking
#	Last modified 3 Dec 2015
{
#	Check counts
	counts <- as.matrix(counts)
	nlib <- ncol(counts)
	ntags <- nrow(counts)
	if(nlib>0L && is.null(colnames(counts))) colnames(counts) <- paste0("Sample",1L:nlib)
	if(ntags>0L && is.null(rownames(counts))) rownames(counts) <- 1L:ntags

#	Check lib.size
	if(is.null(lib.size)) lib.size <- colSums(counts)
	if(nlib != length(lib.size)) stop("Length of 'lib.size' must equal number of columns in 'counts'")
	if(any(lib.size==0L)) warning("Zero library size detected.")

#	Check norm.factors
	if(is.null(norm.factors)) norm.factors <- rep(1,ncol(counts))
	if(nlib != length(norm.factors)) stop("Length of 'norm.factors' must equal number of columns in 'counts'")

#	Check samples
	if(!is.null(samples)) {
		samples <- as.data.frame(samples)
		if(nlib != nrow(samples)) stop("Number of rows in 'samples' must equal number of columns in 'counts'")
	}
	
#	Check group
	if(is.null(group)) {
		if(!is.null(samples))
			group <- samples[,1]
		else
			group <- rep(1,ncol(counts))
	}
	group <- dropEmptyLevels(group)
	if(nlib != length(group)) stop("Length of 'group' must equal number of columns in 'counts'")

#	Make data frame of sample information
	sam <- data.frame(group=group,lib.size=lib.size,norm.factors=norm.factors)
	if(!is.null(samples)) sam <- data.frame(sam, samples)
	samples <- sam
	if(anyDuplicated(colnames(counts))) {
		message("Repeated column names found in count matrix")
		row.names(samples) <- 1L:nlib
	} else 
		row.names(samples) <- colnames(counts)

#	Make object
	x <- new("DGEList",list(counts=counts,samples=samples))

#	Add data frame of gene information
	if(!is.null(genes)) {
		genes <- as.data.frame(genes, stringsAsFactors=FALSE)
		if(nrow(genes) != ntags) stop("Counts and genes have different numbers of rows")
		x$genes <- genes
	}

#	Remove rows with all zeros
	if(remove.zeros) {
		all.zeros <- rowSums(counts>0,na.rm=TRUE)==0
		if(any(all.zeros)) {
			x <- x[!all.zeros,]
			message("Removing ",sum(all.zeros)," rows with all zero counts")
		}
	}

#	x$offset <- expandAsMatrix(getOffset(x),dim(counts))
#	x$weights <- matrix(1,ntags,nlib)

	x
}

