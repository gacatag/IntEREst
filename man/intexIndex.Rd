\name{intexIndex}
\alias{intexIndex}
\title{
	Extract index of intron or exon rows
}
\description{
	Extract row numbers of introns or exons in \code{@sampleAnnotation} of an object of type \code{InterestResult}.
}
\usage{
intexIndex(x, intExCol="int_ex", what="intron")
}
\arguments{
\item{x}{ 
	Object of type \code{InterestResult}.
}
\item{intExCol}{
	Column name (or number) that represents whether each row is "intron" or "exon" in \code{@sampleAnnotation} of \code{x}.
}
\item{what}{
	A character string that defines whether the function should return the index for the introns or exons. Accepts either \code{"exon"} or \code{"intron"} (default) as values.
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

