topTags <- function(object,n=10,adjust.method="BH",sort.by="p.value") 
#	Summary table of the n most differentially expressed tags
#	Mark Robinson, Davis McCarthy, Gordon Smyth
#	Created September 2008.  Last modified 19 July 2010.
{
	sort.by <- match.arg(sort.by,c("p.value","logFC"))
	alfc <-abs(object$table$logFC)
	switch(sort.by,
		"logFC" = {o <- order(alfc,decreasing=TRUE)},
		"p.value" = {o <- order(object$table$p.value,1/alfc)}
	)
	adj.p.val <- p.adjust(object$table$p.value,adjust.method)
	chosen <- o[1:min(nrow(object$table),n)]
	tab <- data.frame(logConc=object$table$logConc[chosen], logFC=object$table$logFC[chosen], p.value=object$table$p.value[chosen] ,adj.p.val=adj.p.val[chosen])
	if(!is.null(object$genes)){
		tab <- cbind(object$genes[chosen,,drop=FALSE], tab) # Assumes that object$genes is a data.frame
	}
	rownames(tab) <- rownames(object$table[chosen,])
	colnames(tab) <- c(colnames(object$genes), "logConc","logFC","PValue","FDR")
	new("TopTags",list(table=tab, adjust.method=as.character(adjust.method),comparison=as.character(object$comparison)))
}
