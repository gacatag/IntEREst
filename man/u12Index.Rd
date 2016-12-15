\name{u12Index}
\alias{u12Index}
\title{
	Extract index of U12 introns rows
}
\description{
	Extract row numbers of U12 introns in \code{@interestDf} in an object of class \code{interestResult}.
}
\usage{
u12Index(x, intExCol="int_ex", intTypeCol="int_type")
}
\arguments{
\item{x}{ 
	Object of type \code{interestResult}.
}
\item{intExCol}{
	Column name (or number) that represents whether each row is "intron" or "exon" in \code{@sampleAnnotation} of \code{x}.
}
\item{intTypeCol}{
	Column name (or number) that represents what type of intron each row of \code{@sampleAnnotation} of \code{x} represents. The values in the column are either "U12", "U2" or \code{NA} (for the exons).
}

}
\value{
A numeric vector which includes the index of U12 introns.
}

\seealso{
	\code{\link{u12NbIndex}}
}

\author{Ali Oghabian}

\examples{

head(u12Index(mdsChr22Obj))

}

