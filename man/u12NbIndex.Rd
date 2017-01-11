\name{u12NbIndex}
\alias{u12NbIndex}
\title{
	Extract index of U2 introns (up/down)stream of U12 introns rows
}
\description{
Extract row numbers of U2-type introns (up/down)stream of U12-type introns (in 
the \code{@interestDf} attribute of an object of class \code{InterestResult}).
}
\usage{
u12NbIndex(x, intExCol="int_ex", intTypeCol="int_type", 
	strandCol="strand", geneIdCol="collapsed_transcripts", 
	naUnstrand=FALSE)
}
\arguments{
\item{x}{ 
Object of type \code{InterestResult}.
}
\item{intExCol}{
Column name (or number) that represents whether each row is "intron" or "exon" 
in \code{@interestDf} attribute of \code{x}.
}
\item{intTypeCol}{
Column name (or number) that represents what type of intron each row of 
\code{@interestDf} attribute of \code{x} represents. The values in the column 
are either "U12", "U2" or \code{NA} (for the exons).
}
\item{strandCol}{
Column name (or number) that represents the strand of each row of 
\code{@interestDf} in \code{x}. The values in the column are either "+", "-" or
"*".
}
\item{geneIdCol}{
Column name (or number) that represents the gene ID of each row of 
\code{@interestDf} in \code{x}.
}
\item{naUnstrand}{
Replace unstranded results, i.e. introns or exon with "*" strand, with NA. If 
set as FALSE (default) "*" strand would be same as "+" strand.
}

}
\value{
\item{upIntron}{
A numeric vector which includes the index of U2-type intron upstream the 
U12-type introns.
}
\item{downIntron}{
A numeric vector which includes the index of U2-type intron downstream the 
U12-type introns.
}
\item{upExon}{
A numeric vector which includes the index of exon upstream the U12-type 
introns.
}
\item{downExon}{
A numeric vector which includes the index of exon downstream the U12-type 
introns.
}
}

\seealso{
	\code{\link{u12Index}}
}

\author{Ali Oghabian}

\examples{


head(u12NbIndex(mdsChr22Obj, intExCol="int_ex", 
	intTypeCol="int_type", strandCol="strand", 
	geneIdCol="ens_trans_id", naUnstrand=FALSE))

head(u12NbIndex(mdsChr22Obj, intExCol="int_ex", 
	intTypeCol="int_type", strandCol="strand", 
	geneIdCol="ens_trans_id", naUnstrand=TRUE))

}

