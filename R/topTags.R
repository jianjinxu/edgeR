topTags <- function(object,n=10,adjust.method="BH",sort.by="p.value") 
#	Summary table of the n most differentially expressed tags
#	Mark Robinson, Davis McCarthy, Gordon Smyth
#	Created September 2008.  Last modified 19 July 2010.
{
	sort.by <- match.arg(sort.by,c("p.value","logFC"))
    tabnames <- names(object$table)
    if( is(object, "DGELRT") & ncol(object$table) > 4 ) {
        if( sort.by=="logFC")
            warning("Two or more logFC columns in DGELRT object. First logFC column used to rank by logFC.\n")
        alfc <- abs(object$table[,2])
    }
    else
        alfc <- abs(object$table$logFC)
	switch(sort.by,
		"logFC" = {o <- order(alfc,decreasing=TRUE)},
		"p.value" = {o <- order(object$table$p.value,1/alfc)}
	)
	adj.p.val <- p.adjust(object$table$p.value,adjust.method)
	chosen <- o[1:min(nrow(object$table),n)]
    tab <- cbind(object$table[chosen,], adj.p.val[chosen])
    if(!is.null(object$genes)){
		tab <- cbind(object$genes[chosen,,drop=FALSE], tab) # Assumes that object$genes is a data.frame
	}
	rownames(tab) <- rownames(object$table[chosen,])
    if( is(object, "DGELRT") ) {
        lfc.cols <- grep("logFC",tabnames)
        colnames(tab) <- c(colnames(object$genes),"logConc",sub("logFC.","",tabnames[lfc.cols]), "LR","P.Value","adj.P.Val")
    }
    else {
        tabnames[grep("p.value", tabnames)] <- "P.Value"
        colnames(tab) <- c(colnames(object$genes), tabnames, "adj.P.Val")
    }
	new("TopTags",list(table=tab, adjust.method=as.character(adjust.method),comparison=as.character(object$comparison)))
}

