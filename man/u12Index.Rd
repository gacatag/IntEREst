\name{u12Index}
\alias{u12Index}
\title{
	Extract index of U12 introns rows
}
\description{
Extract row numbers of U12 introns in an object of class 
\code{SummarizedExperiment}.
}
\usage{
u12Index(x, intExCol="int_ex", intTypeCol="int_type")
}
\arguments{
\item{x}{ 
Object of type \code{SummarizedExperiment}.
}
\item{intExCol}{
Column name (or number) that represents whether each row of \code{x} assays is 
"intron" or "exon".
}
\item{intTypeCol}{
Column name (or number) that represents what type of intron each row of 
\code{x} assays represents.
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

