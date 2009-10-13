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

setClass("SmoothList",
#  Linear model fit (Sep 2009)
representation("list")
)

setClass("TopTags",
representation("list")
)

setMethod("show", "TopTags", function(object) {
	if(length(object$comparison)) cat("Comparison of groups: ", object$comparison[2],"-",object$comparison[1],"\n")
	colnames(object$table) <- c("logConc","logFC","PValue","FDR")
	if(object$adjust.method %in%  c("holm", "hochberg", "hommel", "bonferroni")) colnames(object$table)[4] <- "FWER"
	if(object$adjust.method=="none") object$table$FDR <- NULL
	print(object$table)
})

setIs("DGEList","LargeDataObject")
setIs("EBList","LargeDataObject")
setIs("SmoothList","LargeDataObject")
setIs("deDGEList","LargeDataObject")
setIs("de4DGEList","LargeDataObject")

DGEList <- function(counts=matrix(0), lib.size=NULL, group=factor(), verbose=FALSE, ...) 
{
	if (ncol(counts) != length(group))
		stop("Length of 'group' must equal number of columns in 'counts'")
	if (!is.factor(group))
		group<-as.factor(group)
	if(!is.matrix(counts)) 
		counts<-as.matrix(counts)
        if(is.null(lib.size)) {
                lib.size <- colSums(counts)
                warning("Calculating library sizes from total number of reads for each library.")
        }
	if(length(colnames(counts)) < ncol(counts)) {
			colnames(counts)<-paste("sample",c(1:ncol(counts)),sep=".")
	}
	if (ncol(counts) != length(lib.size))
		stop("Length of 'lib.size' must equal number of columns in 'counts'")
        # remove rows which are all 0
    allZeros <- rowSums(counts,na.rm=TRUE)==0
    if( sum(allZeros) > 0) {
    	counts <- counts[!allZeros,]
        warning("Removing ", sum(allZeros)," rows that all have zero counts.")
    }
    if( verbose ) {
    	cat("----------------------------------------\n")
    	cat("Breakdown of TOTAL counts by ZERO counts\n")
    	cat("----------------------------------------\n")
    	tab <-  table( cut(rowSums(counts),breaks=c((0:9)+.5,1e6)), rowSums(counts==0), dnn=c("rowWiseCountSum","nbrOfZeroObserations") )
    	print(tab)
    	cat("----------------------------------------\n")
    }
    counts <- counts
    samples <- data.frame(group=as.factor(group), lib.size=lib.size)
	x <- list(samples=samples, counts=counts, ...)
	row.names(x$samples) <- colnames(x$counts)
	new("DGEList",x)
}


getCounts <- function(object)
{
  object$counts
}

