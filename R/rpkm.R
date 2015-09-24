rpkm <- function(x, ...)
UseMethod("rpkm")

rpkm.DGEList <- function(x, gene.length=NULL, normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25, ...)
#	Counts per million for a DGEList
#	Gordon Smyth.
#	Created 18 March 2013. Last modified 1 November 2012
{
#	Try to find gene lengths
#	If column name containing gene lengths isn't specified,
#	then will try "Length" or "length" or any column name containing "length"
	if(is.character(gene.length)) {
		gene.length <- x$genes[[gene.length[1]]]
		if(is.null(gene.length)) stop("gene.length column not found")
	} else {
		if(is.null(gene.length)) gene.length <- x$genes$Length
		if(is.null(gene.length)) gene.length <- x$genes$length
		if(is.null(gene.length)) {
			j <- grep("length",tolower(names(x$genes)))
			if(length(j)==1)
				gene.length <- x$genes[,j]
			else
				stop("Gene lengths not found")
		}
	}

	lib.size <- x$samples$lib.size
	if(normalized.lib.sizes) lib.size <- lib.size*x$samples$norm.factors

	rpkm.default(x=x$counts,gene.length=gene.length,lib.size=lib.size,log=log,prior.count=prior.count, ...)
}

rpkm.default <- function(x, gene.length, lib.size=NULL, log=FALSE, prior.count=0.25, ...)
#	Reads per kilobase of gene length per million reads of sequencing
#	Gordon Smyth
#	Created 1 November 2012. Last modified 18 March 2014.
{
	y <- cpm.default(x=x,lib.size=lib.size,log=log,prior.count=prior.count, ...)
	gene.length.kb <- gene.length/1000
	if(log)
		y-log2(gene.length.kb)
	else
		y/gene.length.kb
}

