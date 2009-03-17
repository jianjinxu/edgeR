makeWindowLookupTable <- function(indexes, offsets, starts, ends) {

	getProbesInWindow <- function(indexes, offsets, start, end) {
	#indexes & offsets are integer vectors at this stage
	#start & end are atomic integers - offset co-ordinates we are interested in
		if (length(indexes)==0) return(integer()) else return(indexes[which(offsets>start&offsets<end)])
	}

	lookupTable <- vector(mode='list', length=length(starts))
	names(lookupTable) <- (starts+ends)/2
	for (i in 1:length(starts)) {
		cat("Processing",starts[i],"to",ends[i],"\n")
		lookupTable[[i]] = vector(mode='list', length=length(indexes))
		if (!is.null(names(indexes))) names(lookupTable[[i]]) <- names(indexes)
		for (j in 1:length(indexes)) lookupTable[[i]][[j]] = getProbesInWindow(indexes[[j]], offsets[[j]], starts[i], ends[i])

	}
	return(lookupTable)
}



scoreCorrelation <- function(lookup, intensities, correlateTo, minProbes=1, cor.method="pearson") {
	windowMeans <- function(indexes, intensities) {
		return(mean(intensities[indexes]))
	}

	useLookup = sapply(lookup, length)>=minProbes
	return(cor(sapply(lookup[useLookup], windowMeans, intensities),correlateTo[useLookup],method=cor.method))
}

#Draw graph of correlation vs distance from TSS
correlationGraphs <- function(intensities, lookup, compareExpression, ...) {
	corScores = matrix(NA, nrow=ncol(intensities), ncol=length(lookup), dimnames=list(colnames(intensities),names(lookup)))
	for (i in 1:ncol(intensities)) corScores[i,] = sapply(lookup, scoreCorrelation, intensities[,i], compareExpression, minProbes=2)

	cols = rainbow(nrow(corScores))
	plot(0, type="n", xlim=c(-7500,2500), ylim=c(-1,1), xlab="Distance from TSS", ylab="Correlation with Expression", ...)
	for (i in 1:nrow(corScores)) lines(as.integer(colnames(corScores)), corScores[i,], col=cols[i])
	legend("topleft", legend=rownames(corScores), col=cols, lty=1)
}

