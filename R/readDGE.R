readDGE <- function(files,path=NULL,columns=c(1,2),group=NULL,labels=NULL,...) 
#	Read and collate a set of count data files, each file containing counts for one library
#	Created 2006.  Last modified 15 August 2014.
{
	x <- list()
	if(is.data.frame(files)) {
		x$samples <- files
		if(is.null(labels)) labels <- row.names(files)
		files <- files$files
	} else {
		x$samples <- data.frame(files=as.character(files),group=1,stringsAsFactors=FALSE)
	}
	if(!is.null(group)) x$samples$group <- group
	if(!is.null(x$samples$group)) x$samples$group <- as.factor(x$samples$group)
	d <- taglist <- list()
	for (fn in files) {
		if(!is.null(path)) fn <- file.path(path,fn)
		d[[fn]] <- read.delim(fn,...,stringsAsFactors=FALSE)
		taglist[[fn]] <- as.character(d[[fn]][,columns[1]])
		if(any(duplicated(taglist[[fn]]))) {
			stop(paste("Repeated tag sequences in",fn)) 
		}
	}
	tags <- unique(unlist(taglist))
	ntags <- length(tags)
	nfiles <- length(files)
	x$counts <- matrix(0,ntags,nfiles)
	rownames(x$counts) <- tags
	colnames(x$counts) <- labels
	if(is.null(colnames(x$counts))) colnames(x$counts) <- removeExt(files)
	for (i in 1:nfiles) {
		aa <- match(taglist[[i]],tags)
		x$counts[aa,i] <- d[[i]][,columns[2]]
	}
	x$samples$lib.size <- colSums(x$counts)
	x$samples$norm.factors <- 1
	row.names(x$samples) <- colnames(x$counts)
	x$genes <- NULL
	new("DGEList",x)
}
