\name{intexIndex}
\alias{intexIndex}
\title{
	Extract index of intron or exon rows
}
\description{
Extract row numbers where introns (or exons dependant on user's request) 
are located in an object of type \code{SummarizedExperiment}.
}
\usage{
intexIndex(x, intExCol="int_ex", what="intron")
}
\arguments{
\item{x}{ 
Object of type \code{SummarizedExperiment}.
}
\item{intExCol}{
Column name (or number) that represents whether each row is "intron" or "exon"
 in \code{rowData} of \code{x}.
}
\item{what}{
A character string that defines whether the index for the introns or exons 
should be returned. Accepts either \code{"exon"} or \code{"intron"} (default)
as values.
}

}
\value{
A numeric vector which includes the index of the introns/exons.
}

\seealso{
	\code{\link{u12NbIndex}}
}

\author{Ali Oghabian}

\examples{

# Show the few first index of rows that represent the introns
head(intexIndex(mdsChr22Obj, what="intron"))

}

