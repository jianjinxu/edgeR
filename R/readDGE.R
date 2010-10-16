readDGE <- function(files,path=NULL,columns=c(1,2),group=NULL,labels=NULL,...) 
#	Read and collate a set of DGE data files, one library per file
#	Last modified 16 October 2010.
{
	x <- list()
	if(is.data.frame(files)) {
		x$samples <- files
		if(is.null(labels)) labels <- row.names(files)
		files <- files$files
	} else {
		x$samples <- data.frame(files=as.character(files),stringsAsFactors=FALSE)
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
