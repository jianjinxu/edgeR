#  Code to process hairpin reads from Illumina sequencer
#  Assume fixed structure of read:
#  Barcode + Common sequence + Hairpin sequence

processHairpinReads = function(readfile, barcodefile, hairpinfile,
                    barcodeStart=1, barcodeEnd=5, hairpinStart=37, hairpinEnd=57,
                    allowShifting=FALSE, shiftingBase = 3,
                    allowMismatch=FALSE, barcodeMismatchBase = 1, hairpinMismatchBase = 2,
                    allowShiftedMismatch = FALSE, 
                    verbose = FALSE) {

  # Check file existence
  if ((length(readfile) == 1) && (!file.exists(readfile)))
    stop("Read file doesn't exist.\n")
  if (length(readfile) > 1){
    for(i in 1:length(readfile)){
      if (!file.exists(readfile[i]))
        stop(paste("Read file ", readfile[i], " doesn't exist. \n", sep="")) 
    }
  }
  if (!file.exists(barcodefile))
    stop("Barcode file doesn't exist.\n")
  if (!file.exists(hairpinfile))
    stop("Hairpin file doesn't exist.\n")

  # Validating params
  reads <- file(readfile[1], "rt");
  first_read <- readLines(reads, 2)
  readlength <- nchar(first_read[2])

  if (barcodeStart > barcodeEnd)
    stop("Barcode start position is greater than barcode end position.\n")
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

  # check that barcodes and hairpins provided have no duplicates, are in specified length.
  barcodelength <- barcodeEnd - barcodeStart + 1;
  barcodes <- read.table(barcodefile, header=FALSE, sep="\t");
  numbc <- nrow(barcodes) 
  barcodeseqs <- as.character(barcodes[,2])
  barcodenames <- as.character(barcodes[,1])
  if ((min(nchar(barcodeseqs)) != barcodelength) || (max(nchar(barcodeseqs)) != barcodelength))
    stop(paste("Barcode sequence length is set to ", barcodelength, ", there are barcode sequence not in specified length.\n", sep=""))
  if (length(unique(barcodeseqs)) != numbc)
    stop("There are duplicate barcode sequences.\n")
  if (length(unique(barcodenames)) != numbc)
    stop("There are duplicate barcode names.\n")
  tempbarcodefile <- paste("Barcode", as.character(Sys.getpid()), "temp.txt", sep = "_")
  # passing only barcode sequences to C function
  write.table(barcodes[, 2], file=tempbarcodefile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE);

  hairpinlength <- hairpinEnd - hairpinStart + 1;
  hairpins <- read.table(hairpinfile, header=FALSE, sep="\t");
  numhp <- nrow(hairpins) 
  hairpinseqs <- as.character(hairpins[,2])
  hairpinnames <- as.character(hairpins[,1])
  if ((min(nchar(hairpinseqs)) != hairpinlength) || (max(nchar(hairpinseqs)) != hairpinlength))
    stop(paste("Hairpin sequence length is set to ", hairpinlength, ", there are hairpin sequences not in specified length.\n", sep=""))
  if (length(unique(hairpinseqs)) != numhp)
    stop("There are duplicate hairpin sequences.\n")
  if (length(unique(hairpinnames)) != numhp)
    stop("There are duplicate hairpin names.\n")
  
  # passing only hairpin sequences to C function
  temphairpinfile <- paste("Hairpin", as.character(Sys.getpid()), "temp.txt", sep = "_")
  write.table(hairpins[, 2], file=temphairpinfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE);
    
  if (allowShifting) {
      if ((shiftingBase <= 0) || (shiftingBase > 5))
          stop("To allow hairpin matching at a shifted position, please input a positive shiftingBase no greater than 5. ")	 
  }
  
  if (allowMismatch) {
      if ((barcodeMismatchBase <= 0) || (barcodeMismatchBase > 2))	
          stop("To allow mismatch in barcode sequence, please input a positive barcodeMismatchBase no greater than than 2. ")
      if ((hairpinMismatchBase <= 0) || (hairpinMismatchBase > 4))  
          stop("To allow mismatch in hairpin sequence, please input a positive hairpinMismatchBase no greater than than 4. ")	 
  }

  if (allowShiftedMismatch) {
      if ((!allowShifting) || (!allowMismatch)){
          stop("allowShiftedMismatch option can only be turned on when allowShiting and allowMismatch are both TRUE. ")
      } 
  } 
  tempoutfile <- paste("ReadcountSummary", as.character(Sys.getpid()), "output.txt", sep = "_")

  .C("processHairpinReads", readfile, as.integer(length(readfile)), as.character(tempbarcodefile), as.character(temphairpinfile),
     as.integer(barcodeStart), as.integer(barcodeEnd), as.integer(hairpinStart), as.integer(hairpinEnd),
     as.integer(allowShifting), as.integer(shiftingBase),
     as.integer(allowMismatch), as.integer(barcodeMismatchBase), as.integer(hairpinMismatchBase),
     as.integer(allowShiftedMismatch),
     as.character(tempoutfile), as.integer(verbose))
	 
  summary = read.table(tempoutfile, sep="\t", header=FALSE)
  file.remove(tempoutfile)
  file.remove(tempbarcodefile)
  file.remove(temphairpinfile)
  if (nrow(summary) != length(hairpinnames))
      stop("Number of hairpins from result count matrix doesn't agree with given hairpin list. ")
  if (ncol(summary) != length(barcodenames))
      stop("Number of barcodes from result count matrix doesn't agree with given barcode list. ")
  colnames(summary) = barcodenames
  rownames(summary) = hairpinnames
  
  x <- DGEList(counts = summary, genes = data.frame(ID=hairpinnames,Sequence=hairpinseqs))
}

