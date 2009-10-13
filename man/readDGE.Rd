\name{readDGE}
\alias{readDGE}

\title{Read and merge a set of files containing DGE data}

\description{Reads and merges a set of text files containing digital gene expression data.}

\usage{readDGE(files,path=NULL,columns=c(1,2),...)}

\arguments{ 
\item{files}{character vector of filenames, or alternatively a data.frame with a column containing the file names.}
\item{path}{character string giving the directory containing the files.
The default is the current working directory.}
\item{columns}{which two columns contain the tag names and counts, respectively}
\item{...}{other are passed to \code{read.delim}}
}

\details{
Each file is assumed to contained digital gene expression data for one sample, with transcript identifiers in the first column and counts in the second column.
By default, the files are assumed to be tab-delimited and to contain column headings.
The function forms the union of all transcripts and creates one big table with zeros where necessary.
}

\value{DGEList object}

\author{Mark Robinson and Gordon Smyth}

\examples{
#  Read all .txt files from current working directory

\dontrun{files <- dir(pattern="*\\\\.txt$")
RG <- readDGE(files)}
}
\keyword{file}
