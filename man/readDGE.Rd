\name{readDGE}

\alias{readDGE}

\title{Read a list of files containing DGE data}

\description{Reads a list of text files, one for each sample.  Files should be tab-delimited with an identifier (could be tag sequence) as the first column and counts as the second column.  The function creates one big table with 0s where necessary.}


\usage{ 
readDGE(files,...)
}
\arguments{ 

\item{files}{character vector of filenames}

\item{...}{option arguments to send to \code{read.table}}

}

\value{ list with elements \code{data} (table of counts), \code{lib.size} (library sizes)}

\author{Mark Robinson}

\examples{
#  Read all .txt files from current working directory

\dontrun{files <- dir(pattern="*\\\\.txt$")
RG <- readDGE(files,sep="\t",header=TRUE,comment.char="",stringsAsFactors=FALSE)}
}
\keyword{file}
