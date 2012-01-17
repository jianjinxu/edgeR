setClass("TopTags",
representation("list")
)

setMethod("show", "TopTags", function(object) {
	if(object$test=="exact") {
		cat("Comparison of groups: ",paste(rev(object$comparison),collapse="-"),"\n")
	} else {
		cat("Coefficient: ",object$comparison,"\n")
	}
	print(object$table)
})

as.data.frame.TopTags <- function(x,row.names=NULL,optional=FALSE,...)
{
	if(!is.null(row.names)) row.names(x$table) <- row.names
	x$table
}

topTags <- function(object,n=10,adjust.method="BH",sort.by="p.value") 
#	Summary table of the n most differentially expressed tags
#	Mark Robinson, Davis McCarthy, Gordon Smyth
#	Created September 2008.  Last modified 16 Jan 2012.
{
	sort.by <- match.arg(sort.by,c("p.value","logFC"))
	tabnames <- names(object$table)
	if( is(object, "DGELRT") & ncol(object$table) > 4 ) {
		if( sort.by=="logFC")
			warning("Two or more logFC columns in DGELRT object. First logFC column used to rank by logFC.\n")
		alfc <- abs(object$table[,1])
	} else {
		alfc <- abs(object$table$logFC)
	}
	switch(sort.by,
		"logFC" = {o <- order(alfc,decreasing=TRUE)},
		"p.value" = {o <- order(object$table$PValue,1/alfc)}
	)
	chosen <- o[1:min(nrow(object$table),n)]
	tab <- object$table[chosen,]

	FWER.methods <- c("holm", "hochberg", "hommel", "bonferroni")
	FDR.methods <- c("BH", "BY", "fdr")
	adjust.method <- match.arg(adjust.method,c(FWER.methods,FDR.methods,"none"))
	if(adjust.method != "none") {
		adj.p.val <- p.adjust(object$table$PValue,adjust.method)
		if(adjust.method %in% FWER.methods) adjustment <- "FWER"
		if(adjust.method %in% FDR.methods) adjustment <- "FDR"
		tab[[adjustment]] <- adj.p.val[chosen]
	}
	if(!is.null(object$genes)){
		if(is.null(dim(object$genes))) object$genes <- data.frame(ID=object$genes,stringsAsFactors=FALSE)
		tab <- cbind(object$genes[chosen,,drop=FALSE], tab)
	}
	if(is(object,"DGEExact")) test <- "exact" else test <- "glm"
	new("TopTags",list(
		table=tab,
		adjust.method=adjust.method,
		comparison=as.character(object$comparison),
		test=test
	))
}

