#  Code to process hairpin reads from Illumina sequencer
#  Assume basic fixed structure of read:
#  Barcode + Common sequence + Hairpin sequence
#  If barcode2Start and barcode2End are not NULL, forward read sequences contain a second barcode
#  If readfile2, barcodeStartRev and barcodeEndRev are not NULL, reverse read sequences contain reverse barcodes

processAmplicons = function(readfile, readfile2=NULL, barcodefile, hairpinfile,
                    barcodeStart=1, barcodeEnd=5, barcode2Start=NULL, barcode2End=NULL, barcodeStartRev=NULL, barcodeEndRev=NULL,
                    hairpinStart=37, hairpinEnd=57,
                    allowShifting=FALSE, shiftingBase = 3,
                    allowMismatch=FALSE, barcodeMismatchBase = 1, hairpinMismatchBase = 2,
                    allowShiftedMismatch = FALSE, 
                    verbose = FALSE) {
    
  checkFileExistence = function(readfilenames){
    if ((length(readfilenames) == 1) && (!file.exists(readfilenames)))
      stop("Read file doesn't exist.\n")
    if (length(readfilenames) > 1){
      for(i in 1:length(readfilenames)){
        if (!file.exists(readfilenames[i]))
          stop(paste("Read file ", readfilenames[i], " doesn't exist. \n", sep="")) 
      }
    }
  }
  
  # Check file existence
  numfiles = length(readfile)
  checkFileExistence(readfile);
  if (is.null(readfile2)) {
    IsPairedReads = FALSE;
  } else {
    IsPairedReads = TRUE;
    if (length(readfile) != length(readfile2))
      stop("readfile and readfile2 should match each other.")
    checkFileExistence(readfile2);
  } 
  IsDualIndexingOnForwardRead = !is.null(barcode2Start) && !is.null(barcode2End)
      
  if (!file.exists(barcodefile))
    stop("Barcode file doesn't exist.\n")
  if (!file.exists(hairpinfile))
    stop("Hairpin file doesn't exist.\n")

  # Validate input parameters
  reads <- file(readfile[1], "rt");
  first_read <- readLines(reads, 2)
  readlength <- nchar(first_read[2])

  if ((barcodeStart < 1) || (barcodeStart > readlength))
    stop("Invalid barcode start position!\n")
  if ((barcodeEnd < 1) || (barcodeEnd > readlength))
    stop("Invalid barcode end position!\n")
  if (barcodeEnd <= barcodeStart)
    stop("Barcode end position should be greater than barcode start position. \n")
  if ((hairpinStart < 1) || (hairpinStart > readlength))
    stop("Invalid hairpin start position!")
  if ((hairpinEnd < 1) || (hairpinEnd > readlength))
    stop("Invalid hairpin end position!")
  if (hairpinEnd <= hairpinStart)
    stop("Hairpin end position should be greater than hairpin start position. \n")
  
  if (IsDualIndexingOnForwardRead){
    if ((barcode2Start < 1) || (barcode2Start > readlength))
      stop("Invalid barcode2 start position!\n")
    if ((barcode2End < 1) || (barcode2End > readlength))
      stop("Invalid barcode2 end position!\n")
    if (barcode2End <= barcode2Start)
      stop("Barcode2 end position should be greater than barcode2 start position. \n")
  }
      
  close(reads)
  
  if (IsPairedReads) {
    reads <- file(readfile2[1], "rt");
    first_read <- readLines(reads, 2)
    readlength2 <- nchar(first_read[2])
    close(reads)
    
    if ( (is.null(barcodeStartRev)) || (is.null(barcodeEndRev)) )
      stop("readfile2 is supplied, barcodeStartRev and barcodeEndRev should be specified. ")
    if ((barcodeStartRev < 1) || (barcodeStartRev > readlength2))
      stop("Invalid reverse barcode start position!\n")
    if ((barcodeEndRev < 1) || (barcodeEndRev > readlength2))
      stop("Invalid reverse barcode end position!\n")
    if (barcodeEndRev <= barcodeStartRev)
      stop("Reverse barcode end position should be greater than reverse barcode start position. \n")
  }
  
  # Validate barcodes
  barcodelength <- barcodeEnd - barcodeStart + 1;
  barcodes <- read.table(barcodefile, header=TRUE, sep="\t");
  barcodeIDIndex = which(colnames(barcodes) == 'ID')
  if (length(barcodeIDIndex) < 1) 
    stop("Can't find column ID in ", barcodefile)
  barcodeseqIndex = which(colnames(barcodes) == 'Sequences')
  if (length(barcodeseqIndex) < 1) 
    stop("Can't find column Sequences in ", barcodefile)
  barcodeIDs <- as.character(barcodes[, barcodeIDIndex]) 
  barcodeseqs <- as.character(barcodes[, barcodeseqIndex]) 
  if (anyDuplicated(barcodeIDs))
    stop("There are duplicate barcode IDs.\n")
  if ((min(nchar(barcodeseqs)) != barcodelength) || (max(nchar(barcodeseqs)) != barcodelength))
    stop(paste("Barcode sequence length is set to ", barcodelength, ", there are barcode sequence not with specified length.\n", sep=""))

  if (IsPairedReads) {
    barcodeseqRevIndex = which(colnames(barcodes) == 'SequencesReverse')
    if (length(barcodeseqRevIndex) < 1) 
      stop("Can't find column SequencesReverse in ", barcodefile)
    barcodeseqsReverse <- as.character(barcodes[, barcodeseqRevIndex])
    barcodelengthReverse <- barcodeEndRev - barcodeStartRev + 1;	
    if ((min(nchar(barcodeseqsReverse)) != barcodelengthReverse) || (max(nchar(barcodeseqsReverse)) != barcodelengthReverse))
      stop(paste("Reverse barcode sequence length is set to ", barcodelength, ", there are reverse barcode sequence not in specified length.\n", sep=""))  
    concatenatedBarcodeseqs = paste(barcodeseqs, barcodeseqsReverse, sep="")
    if (anyDuplicated(concatenatedBarcodeseqs)) 
      stop("There are duplicate forward/reverse barcode sequences.\n")
  } else if (IsDualIndexingOnForwardRead) {
    barcodeseq2Index = which(colnames(barcodes) == 'Sequences2')
    if (length(barcodeseq2Index) < 1) 
      stop("Can't find column Sequences2 in ", barcodefile)
    barcode2seqs <- as.character(barcodes[, barcodeseq2Index])
    barcode2length <- barcode2End - barcode2Start + 1;

    if ((min(nchar(barcode2seqs)) != barcode2length) || (max(nchar(barcode2seqs)) != barcode2length))
      stop(paste("Forward barcode2 sequence length is set to ", barcode2length, ", there are barcode2 sequence not in specified length.\n", sep=""))  
    concatenatedBarcodeseqs = paste(barcodeseqs, barcode2seqs, sep="")
    if (anyDuplicated(concatenatedBarcodeseqs)) 
      stop("There are duplicate barcode/barcode2 sequences.\n") 
  } else {
    if (anyDuplicated(barcodeseqs)) 
      stop("There are duplicate barcode sequences.\n")
  }
  
  # Validate hairpins
  hairpinlength <- hairpinEnd - hairpinStart + 1;
  hairpins <- read.table(hairpinfile, header=TRUE, sep="\t");
  hairpinIDIndex = which(colnames(hairpins) == 'ID')
  if (length(hairpinIDIndex) < 1) 
    stop("Can't find column ID in ", hairpinfile)
  hairpinIDs <- as.character(hairpins[, hairpinIDIndex])
  hairpinSeqIndex = which(colnames(hairpins) == 'Sequences')
  if (length(hairpinSeqIndex) < 1) 
    stop("Can't find column Sequences in ", hairpinfile)
  hairpinseqs <- as.character(hairpins[, hairpinSeqIndex])

  if ((min(nchar(hairpinseqs)) != hairpinlength) || (max(nchar(hairpinseqs)) != hairpinlength))
    stop(paste("Hairpin sequence length is set to ", hairpinlength, ", there are hairpin sequences not with specified length.\n", sep=""))
  if (anyDuplicated(hairpinseqs)) 
    stop("There are duplicate hairpin sequences.\n")
  if (anyDuplicated(hairpinIDs)) 
    stop("There are duplicate hairpin IDs.\n")
 
  # validate mismatch/shifting input parameters
  if (allowShifting) {
    if ((shiftingBase <= 0) || (shiftingBase > 5))
      stop("To allow hairpin matching at a shifted position, please input a positive shiftingBase no greater than 5. ")	 
  }
  if (allowMismatch) {
    if ((barcodeMismatchBase < 0) || (barcodeMismatchBase > 2))	
      stop("To allow mismatch in barcode sequence, please input a non-negative barcodeMismatchBase no greater than than 2. ")
    if ((hairpinMismatchBase < 0) || (hairpinMismatchBase > 4))  
      stop("To allow mismatch in hairpin sequence, please input a non-negative hairpinMismatchBase no greater than than 4. ")	 
  }
  if (allowShiftedMismatch) {
    if ((!allowShifting) || (!allowMismatch)){
      stop("allowShiftedMismatch option can only be turned on when allowShiting and allowMismatch are both TRUE. ")
    } 
  } 

  # passing only barcode/hairpin sequences to C function
  tempbarcodefile <- paste("Barcode", as.character(Sys.getpid()), "temp.txt", sep = "_")
  on.exit({ if (file.exists(tempbarcodefile)) { file.remove(tempbarcodefile) }}, add=TRUE)
  if (IsPairedReads) {
    bothBarcodeSeqs = cbind(barcodeseqs, barcodeseqsReverse)
    write.table(bothBarcodeSeqs, file=tempbarcodefile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE);
  } else if (IsDualIndexingOnForwardRead) {
    bothBarcodeSeqs = cbind(barcodeseqs, barcode2seqs)
    write.table(bothBarcodeSeqs, file=tempbarcodefile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE);
  } else {
    write.table(barcodeseqs, file=tempbarcodefile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE);
  }
  
  temphairpinfile <- paste("Hairpin", as.character(Sys.getpid()), "temp.txt", sep = "_")
  on.exit({ if (file.exists(temphairpinfile)) { file.remove(temphairpinfile) }}, add=TRUE)
  write.table(hairpinseqs, file=temphairpinfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE);

  tempoutfile <- paste("ReadcountSummary", as.character(Sys.getpid()), "output.txt", sep = "_")
  on.exit({ if (file.exists(tempoutfile)) { file.remove(tempoutfile) }}, add=TRUE)

  tryCatch({
    if (!IsPairedReads) {
       readfile2 = rep("DummyReadfile.fastq", numfiles)
       barcodeStartRev = 0;
       barcodeEndRev = 0;
    }
    if (!IsDualIndexingOnForwardRead) {
       barcode2Start = 0;
       barcode2End = 0;
    }
    
    .C(.cprocessHairpinReads, as.integer(IsPairedReads), as.integer(IsDualIndexingOnForwardRead), 
	   as.character(readfile), as.character(readfile2), as.integer(numfiles),
       as.character(tempbarcodefile), as.character(temphairpinfile),
       as.integer(barcodeStart), as.integer(barcodeEnd), as.integer(barcode2Start), as.integer(barcode2End), as.integer(barcodeStartRev), as.integer(barcodeEndRev),
       as.integer(hairpinStart), as.integer(hairpinEnd),
       as.integer(allowShifting), as.integer(shiftingBase),
       as.integer(allowMismatch), as.integer(barcodeMismatchBase), as.integer(hairpinMismatchBase),
       as.integer(allowShiftedMismatch),
       as.character(tempoutfile), as.integer(verbose))      

    hairpinReadsSummary <- read.table(tempoutfile, sep="\t", header=FALSE)
  },
  error = function(err) {
    print(paste("ERROR MESSAGE:  ",err))
  }
  )

  if (exists("hairpinReadsSummary")) {
  
    if (nrow(hairpinReadsSummary) != length(hairpinIDs))
      stop("Number of hairpins from result count matrix doesn't agree with given hairpin list. ")
    if (ncol(hairpinReadsSummary) != length(barcodeIDs))
      stop("Number of barcodes from result count matrix doesn't agree with given barcode list. ")
    colnames(hairpinReadsSummary) = barcodeIDs
    rownames(hairpinReadsSummary) = hairpinIDs
    x <- DGEList(counts = hairpinReadsSummary, genes = hairpins)
    if(!is.null(barcodes$group)) {
      x$samples = cbind("ID"=barcodes$ID, "lib.size"=x$samples$lib.size, 
                       "norm.factors"=x$samples$norm.factors,
                       barcodes[,-match(c("ID","Sequences"), colnames(barcodes))])
    } else {
      x$samples = cbind("ID"=barcodes$ID, x$samples, barcodes[,-match(c("ID","Sequences"), colnames(barcodes))])
    }
  } else {
    stop("An error occured in processHairpinReads.")
  }
  x
}
