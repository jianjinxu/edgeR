require(methods)

setClass("deDGEList",
#  Linear model fit
representation("list")
)

setClass("de4DGEList",
#  Updated Linear model fit (Sep 2009)
representation("list")
)

setClass("DGEList",
#  Linear model fit
representation("list")
)

setClass("EBList",
#  Linear model fit
representation("list")
)

setClass("TopTags",
representation("list")
)

setMethod("show", "TopTags", function(object) {
	if(length(object$comparison)) cat("Comparison of groups: ", object$comparison[2],"-",object$comparison[1],"\n")
	#colnames(object$table) <- c("logConc","logFC","PValue","FDR")
	if(object$adjust.method %in%  c("holm", "hochberg", "hommel", "bonferroni")) colnames(object$table)[ncol(object$table)] <- "FWER"
	if(object$adjust.method=="none") object$table$FDR <- NULL
	print(object$table)
})

setIs("DGEList","LargeDataObject")
setIs("EBList","LargeDataObject")
setIs("deDGEList","LargeDataObject")
setIs("de4DGEList","LargeDataObject")

DGEList <- function(counts=matrix(0,0,0), lib.size=NULL, group=factor(), genes=NULL, remove.zeros=TRUE) 
#	Construct DGEList object from components, with some checking
#	Last modified  22 Jan 2010
{
	counts <- as.matrix(counts)
	nlib <- ncol(counts)
	ntags <- nrow(counts)
	if(nlib>0 & is.null(colnames(counts)))
		colnames(counts) <- paste("sample",1:ncol(counts),sep=".")
	group <- as.factor(group)
	if(nlib != length(group))
		stop("Length of 'group' must equal number of columns in 'counts'")
	if(is.null(lib.size)) {
		lib.size <- colSums(counts)
		message("Calculating library sizes from total number of reads for each library.")
	} else {
		if(nlib != length(lib.size))
			stop("Length of 'lib.size' must equal number of columns in 'counts'")
	}
	samples <- data.frame(group=group,lib.size=lib.size)
	row.names(samples) <- colnames(counts)
	x <- new("DGEList",list(samples=samples,counts=counts))
	if(!is.null(genes)) {
		genes <- as.data.frame(genes, stringsAsFactors=FALSE)
		if(nrow(genes) != nrow(x$counts)) stop("counts and genes have different nrows")
		x$genes <- genes
	}
	if(remove.zeros) {
   	allZeros <- rowSums(counts,na.rm=TRUE)==0
   	if(any(allZeros)) {
			x <- x[!allZeros,]
			warning("Removing",sum(allZeros),"rows that all have zero counts.")
		}
	}
	x
}

getCounts <- function(object)
{
  object$counts
}

