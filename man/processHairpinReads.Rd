\name{processHairpinReads}
\alias{processHairpinReads}

\title{Process raw data from shRNA-seq screens}

\description{
Given a list of barcode sequences and hairpin sequences from a shRNA-seq screen, generate a DGEList of counts from the raw fastq file/(s) containing the sequence reads. 
}

\usage{
processHairpinReads(readfile, barcodefile, hairpinfile,
                    barcodeStart=1, barcodeEnd=5, hairpinStart=37, hairpinEnd=57,
                    allowShifting=FALSE, shiftingBase = 3,
                    allowMismatch=FALSE, barcodeMismatchBase = 1, hairpinMismatchBase = 2,
		    allowShiftedMismatch=FALSE, verbose = FALSE) 
}

\arguments{
\item{readfile}{character vector giving one or more fastq filenames}
\item{barcodefile}{filename containing barcode ids and sequences}
\item{hairpinfile}{filename containing hairpin ids and sequences}
\item{barcodeStart}{numeric value, starting position (inclusive) of barcode sequence in reads}
\item{barcodeEnd}{numeric value, ending position (inclusive) of barcode sequence in reads}
\item{hairpinStart}{numeric value, starting position (inclusive) of hairpin sequence in reads}
\item{hairpinEnd}{numeric value, ending position (inclusive) of hairpin sequence in reads}
\item{allowShifting}{logical, indicates whether a given hairpin can be matched to a neighbouring position}
\item{shiftingBase}{numeric value of maximum number of shifted bases from input \code{hairpinStart} and \code{hairpinEnd} should the program check for a hairpin match when \code{allowShifting} is \code{TRUE}}
\item{allowMismatch}{logical, indicates whether sequence mismatch is allowed}
\item{barcodeMismatchBase}{numeric value of maximum number of base sequence mismatch allowed in barcode when \code{allowShifting} is \code{TRUE}}
\item{hairpinMismatchBase}{numeric value of maximum number of base sequence mismatch allowed in hairpin when \code{allowShifting} is \code{TRUE}}
\item{allowShiftedMismatch}{logical, effective when \code{allowShifting} and \code{allowMismatch} are both \code{TRUE}. It indicates whether we check for sequence mismatch at a shifted position.}
\item{verbose}{if \code{TRUE}, output program progess}
}

\value{Returns a \code{\link[edgeR:DGEList-class]{DGEList}} object with following components:
	\item{counts}{read count matrix tallying up the number of reads with particular barcode and hairpin matches. Each row is a hairpin and each column is a sample}
	\item{genes}{In this case, hairpin information (ID, sequences, corresponding target gene) may be recorded in this data.frame}
	\item{lib.size}{auto-calculated column sum of the counts matrix}
}

\details{
The input barcode file and hairpin files are tab-separated text files with at least two columns (named 'ID' and 'Sequences') containing the sample or hairpin ids and a second column indicating the sample index or hairpin sequences to be matched.  The barcode file may also contain a 'group' column that indicates which experimental group a sample belongs to.  Additional columns in each file will be included in the respective \code{$samples} or \code{$genes} data.frames of the final code{\link[edgeR:DGEList-class]{DGEList}} object.  These files, along with the fastq file/(s) are assumed to be in the current working directory.

To compute the count matrix, the matching to given barcodes and hairpins is conducted in two rounds. The first round looks for an exact sequence match. The program checks for a match from given barcode sequences and hairpin sequences at specified location. If \code{allowShifting} is set to \code{TRUE}, the program also checks if a given hairpin sequence can be found at a neighbouring position in read. For hairpins without a match, the program performs a second round of mapping which allows sequence mismatch. The program checks parameter \code{allowShifting} and \code{allowShiftedMismatch} to see if it allows mismatch sequence at a shifted position. The maximum number of mismatch bases in barcode and hairpin are specified in parameters \code{barcodeMismatchBase} and \code{hairpinMismatchBase}. 

The program outputs a \code{\link[edgeR:DGEList-class]{DGEList}} object, with a count matrix indicating the number of times each barcode and hairpin combination could be matched in reads from input fastq file/(s). 

For further examples and data, refer to the Case studies available from http://bioinf.wehi.edu.au/shRNAseq/.
}

\author{Zhiyin Dai, Matthew Ritchie}

\references{
Dai Z, Sheridan JM, et al. (submitted, 2014). shRNA-seq data analysis with edgeR. \emph{submitted}.
}
