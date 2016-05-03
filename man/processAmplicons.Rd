\name{processAmplicons}
\alias{processAmplicons}

\title{Process raw data from pooled genetic sequencing screens}

\description{
Given a list of sample-specific index (barcode) sequences and hairpin/sgRNA-specific sequences from an amplicon sequencing screen, generate a DGEList of counts from the raw fastq file/(s) containing the sequence reads. 
Assumes fixed structure of amplicon sequences (i.e. both the sample-specific index sequences and hairpin/sgRNA sequences can be found at particular locations within each read).
}

\usage{
processAmplicons(readfile, readfile2=NULL, barcodefile, hairpinfile,
                    barcodeStart=1, barcodeEnd=5, 
                    barcode2Start=NULL, barcode2End=NULL,
                    barcodeStartRev=NULL, barcodeEndRev=NULL, 
                    hairpinStart=37, hairpinEnd=57,
                    allowShifting=FALSE, shiftingBase=3,
                    allowMismatch=FALSE, barcodeMismatchBase=1, 
                    hairpinMismatchBase=2, allowShiftedMismatch=FALSE, 
                    verbose=FALSE)
}

\arguments{
\item{readfile}{character vector giving one or more fastq filenames}
\item{readfile2}{character vector giving one or more fastq filenames for reverse read, default to NULL}
\item{barcodefile}{filename containing sample-specific barcode ids and sequences}
\item{hairpinfile}{filename containing hairpin/sgRNA-specific ids and sequences}
\item{barcodeStart}{numeric value, starting position (inclusive) of barcode sequence in reads}
\item{barcodeEnd}{numeric value, ending position (inclusive) of barcode sequence in reads}
\item{barcode2Start}{numeric value, starting position (inclusive) of second barcode sequence in forward reads}
\item{barcode2End}{numeric value, ending position (inclusive) of second barcode sequence in forward reads}
\item{barcodeStartRev}{numeric value, starting position (inclusive) of barcode sequence in reverse reads, default to NULL}
\item{barcodeEndRev}{numeric value, ending position (inclusive) of barcode sequence in reverse reads, default to NULL}
\item{hairpinStart}{numeric value, starting position (inclusive) of hairpin/sgRNA sequence in reads}
\item{hairpinEnd}{numeric value, ending position (inclusive) of hairpin/sgRNA sequence in reads}
\item{allowShifting}{logical, indicates whether a given hairpin/sgRNA can be matched to a neighbouring position}
\item{shiftingBase}{numeric value of maximum number of shifted bases from input \code{hairpinStart} and \code{hairpinEnd} should the program check for a hairpin/sgRNA match when \code{allowShifting} is \code{TRUE}}
\item{allowMismatch}{logical, indicates whether sequence mismatch is allowed}
\item{barcodeMismatchBase}{numeric value of maximum number of base sequence mismatches allowed in a barcode sequence when \code{allowShifting} is \code{TRUE}}
\item{hairpinMismatchBase}{numeric value of maximum number of base sequence mismatches allowed in a hairpin/sgRNA sequence when \code{allowShifting} is \code{TRUE}}
\item{allowShiftedMismatch}{logical, effective when \code{allowShifting} and \code{allowMismatch} are both \code{TRUE}. It indicates whether we check for sequence mismatches at a shifted position.}
\item{verbose}{if \code{TRUE}, output program progress}
}

\value{Returns a \code{\link[edgeR:DGEList-class]{DGEList}} object with following components:
	\item{counts}{read count matrix tallying up the number of reads with particular barcode and hairpin/sgRNA matches. Each row is a hairpin and each column is a sample}
	\item{genes}{In this case, hairpin/sgRNA-specific information (ID, sequences, corresponding target gene) may be recorded in this data.frame}
	\item{lib.size}{auto-calculated column sum of the counts matrix}
}

\details{The \code{processAmplicons} function assumes the sequences in your fastq files have a fixed structure (as per Figure 1A of http://f1000research.com/articles/3-95/v2). It cannot be used if your hairpins/sgRNAs/sample index sequences are in random locations within each read. You will need to customise your own sequence processing pipeline if this is the case, but can still complete your downstream analysis using edgeR.

The input barcode file and hairpin/sgRNA files are tab-separated text files with at least two columns (named 'ID' and 'Sequences') containing the sample or hairpin/sgRNA ids and a second column indicating the sample index or hairpin/sgRNA sequences to be matched. If \code{barcode2Start} and \code{barcode2End} are specified, a third column 'Sequences2' is expected in the barcode file. If \code{readfile2}, \code{barcodeStartRev} and \code{barcodeEndRev} are specified, another column 'SequencesReverse' is expected in the barcode file.  The barcode file may also contain a 'group' column that indicates which experimental group a sample belongs to. Additional columns in each file will be included in the respective \code{$samples} or \code{$genes} data.frames of the final code{\link[edgeR:DGEList-class]{DGEList}} object.  These files, along with the fastq file/(s) are assumed to be in the current working directory.

To compute the count matrix, matching to the given barcodes and hairpins/sgRNAs is conducted in two rounds. The first round looks for an exact sequence match for the given barcode sequences and hairpin/sgRNA sequences at the locations specified. If \code{allowShifting} is set to \code{TRUE}, the program also checks if a given hairpin/sgRNA sequence can be found at a neighbouring position in the read. If a match isn't found, the program performs a second round of matching which allows for sequence mismatches if \code{allowMismatch} is set to \code{TRUE}. The program also checks parameter \code{allowShiftedMismatch} which accommodates mismatches at the shifted positions. The maximum number of mismatch bases in barcode and hairpin/sgRNA are specified by the parameters \code{barcodeMismatchBase} and \code{hairpinMismatchBase}. 

The program outputs a \code{\link[edgeR:DGEList-class]{DGEList}} object, with a count matrix indicating the number of times each barcode and hairpin/sgRNA combination could be matched in reads from input fastq file/(s).

For further examples and data, refer to the Case studies available from http://bioinf.wehi.edu.au/shRNAseq/.

Note: The \code{processAmplicons} function supercedes the earlier \code{processHairpinReads} function.
}

\author{Zhiyin Dai and Matthew Ritchie}

\references{
Dai Z, Sheridan JM, et al. (2014). edgeR: a versatile tool for the analysis of shRNA-seq and CRISPR-Cas9 genetic screens. F1000Research, 3:95. http://f1000research.com/articles/3-95/v2. PMID: 24860646.
}
