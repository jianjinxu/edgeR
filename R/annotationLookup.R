annotationLookup <- function(probes, annotation, bpUp, bpDown, probeIndex=NULL) {
#probes = dataframe of $chr and $position
#annotation = dataframe of $chr, $position, $strand ("+" or "-"), rownames = annotation name

	processChunk <- function(probePositions, annotation, bpUp, bpDown) {
		#initialise return variables
		numAnnot = length(annotation$position)
		annotProbes = list(indexes=vector(mode='list', length=numAnnot), offsets=vector(mode='list', length=numAnnot))
		if (is.null(probePositions)||(length(probePositions)==0)) return(annotProbes) #if no probes on this chromosome return empty annotations	
		chromosomeSize = max(probePositions)

		#strands lookup conversion
		tempStrand = ifelse(annotation$strand=="+",1,2)
		#strands lookup indexes
		strandIndexes = rbind(-bpUp:bpDown, bpDown:-bpUp)
		strandLookup = rbind(c(-bpUp,bpDown), c(-bpDown,bpUp))

		#allocate vector of chromosome size
		chromosomeLookup <- rep.int(NA, chromosomeSize+bpUp+bpDown) #add bpUp + bpDown to cater for chromosome edges

		#set positions of probes
		chromosomeLookup[probePositions+bpUp] = names(probePositions) #shift for bpUp @ chromosome edge

		for (i in 1:numAnnot) {
			tempStart = annotation$position[i]+strandLookup[tempStrand[i],1]+bpUp #shift for bpUp @ chromosome edge
			tempEnd = annotation$position[i]+strandLookup[tempStrand[i],2]+bpUp #shift for bpUp @ chromosome edge
			tempProbes = chromosomeLookup[tempStart:tempEnd]
			annotProbes$indexes[[i]] = as.integer(na.omit(tempProbes))
			#adjust offsets if at start/end of a chromosome
			annotProbes$offsets[[i]] = strandIndexes[tempStrand[i],which(!is.na(tempProbes))]
		}
		return(annotProbes)
	}

	processChromosome <- function(probePositions, annotation, bpUp, bpDown) {
		numAnnot = length(annotation$position)
		annotProbes = list(indexes=vector(mode='list', length=numAnnot), offsets=vector(mode='list', length=numAnnot))
		if (is.null(probePositions)) return(annotProbes) #if no probes on this chromosome return empty annotations	
		chromosomeSize = max(probePositions)
    annotChunks = split(1:nrow(annotation), trunc(annotation$position/10000000)) #split into 10MB chunks of annotations
    for (i in annotChunks) {
      chunkRange = c(-bpUp, bpUp)+range(annotation$position[i])
      chunkPositions = probePositions[(probePositions>=chunkRange[1])&(probePositions<=chunkRange[2])]-chunkRange[1]+1
      chunkAnnot = annotation[i,]
      chunkAnnot$position = chunkAnnot$position-chunkRange[1]+1
		  tempAnnot = processChunk(chunkPositions, chunkAnnot, bpUp, bpDown)
		  annotProbes$indexes[i] = tempAnnot$indexes
		  annotProbes$offsets[i] = tempAnnot$offsets
	  }
		return(annotProbes)
  }


	tempProbes <- probes$position

	#Use probeIndex supplied, or assume probes are in order
	if (is.null(probeIndex)) names(tempProbes) <- 1:length(tempProbes) else 
names(tempProbes) <- probeIndex
	probesChr = split(tempProbes, probes$chr) #index which probes are on each chromosome
	annotChr = split(1:nrow(annotation), annotation$chr) #ok go by chromosome
	annot = list(indexes=vector(mode='list', length=nrow(annotation)), offsets=vector(mode='list', length=nrow(annotation)))
	for (i in annotChr) {
		thisChr = annotation$chr[i[1]]
		cat("Processing",thisChr,"\n")
		tempAnnot = processChromosome(probesChr[[match(thisChr, names(probesChr))]], annotation[i,], bpUp, bpDown)
		annot$indexes[i] = tempAnnot$indexes
		annot$offsets[i] = tempAnnot$offsets
	}
	if (!is.null(rownames(annotation))) {
		names(annot$indexes) <- rownames(annotation)
		names(annot$offsets) <- rownames(annotation)
	}
	return(annot)
	#returns $indexes = a list for each annotation entry with the indexes of the probes within bpUp & bpDown
	#	 $offsets = a list for each annotation entry with the offsets from the annotation of the probes within bpUp & bpDown
}




