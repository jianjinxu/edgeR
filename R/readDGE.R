readDGE <- function(files,path=NULL,columns=c(1,2),...) 
#	Read and collate a set of DGE data files, one library per file
#	Last modified 7 October 2009.
{
	x <- list()
	if(is.data.frame(files)) {
		x$samples <- files
		files <- files$files
	} else {
		x$samples <- data.frame(files=as.character(files),stringsAsFactors=FALSE)
	}
	d <- taglist <- list()
	for (fn in files) {
		if(!is.null(path)) fn <- file.path(path,fn)
		d[[fn]] <- read.delim(fn,...,stringsAsFactors=FALSE)
		taglist[[fn]] <- as.character(d[[fn]][,columns[1]])
	}
	tags <- unique(unlist(taglist))
	ntags <- length(tags)
	nfiles <- length(files)
	x$counts <- matrix(0,ntags,nfiles)
	rownames(x$counts) <- tags
	colnames(x$counts) <- removeExt(files)
	for (i in 1:nfiles) {
		aa <- match(taglist[[i]],tags)
		x$counts[aa,i] <- d[[i]][,columns[2]]
	}
	x$samples$lib.size <- colSums(x$counts)
        if(!is.null(x$samples$group))
	    x$samples$group <- factor(x$samples$group)
	row.names(x$samples) <- colnames(x$counts)
	new("DGEList",x)
}
