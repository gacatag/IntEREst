\name{applyOverlap}
\alias{applyOverlap}
\title{Apply function over counts}
\description{ 
Runs a function on columns of the counts (assay) of a 'SummarizedExperiment'
object (resulted by \code{interest()}, \code{interest.sequential()} or 
\code{readInterestResults()} ) based on the overalp of its exon/intron
coordinates with those of another 'SummarizedExperiment' object. The number of
the rows and the dimensions of the counts of the result are equal to those of
the \code{subject}. The function is applied on the query based on it's overlap
to the subject.}
\usage{
applyOverlap(
	query,
	subject,
	type="any",
	replaceValues=FALSE,
	intExCol="int_ex", 
	intronExon="intron",
	sujectGeneNamesCol,
	repeatsTableToFilter=c(),
	scaleFragment=TRUE,
	scaleLength=TRUE,
	unmapValue=0,
	FUN=mean,
	...
)
}
\arguments{
\item{query, subject}{
SummarizedExperiment objects resulted by \code{interest()}, 
\code{interest.sequential()} or \code{readInterestResults()} functions. }

\item{type}{
The type of overlap. By default it considers any overlap. See
\link[IRanges]{findOverlaps} for more info.}

\item{replaceValues}{
Whether return a 'SummarizedExperiment' object with new counts (resulted by 
running function) replaced.}

\item{intExCol}{
Column name (or number) in the rowData of the objects that represents whether 
each row of the assay is "intron" or "exon".}

\item{intronExon}{
Should be assigned either \code{'intron'} or \code{'exon'} or 
\code{c('intron','exon')} based on whether match the PWM to the intronic, 
exonic, or intronic and exonic regions of the reference. By default it seeks 
matches in intronic regions (intronExon='intron').}

\item{sujectGeneNamesCol}{
The column in the row data of the subject that includes the gene names.}

\item{repeatsTableToFilter}{A data.frame  table  that includes chr,begin and
end columns. If defined, all reads mapped to the described regions will be 
ingnored.
}

\item{scaleFragment}{ 
Logical value, indicating whether the retention levels must be scaled by 
(genewide) fragment levels.}

\item{scaleLength}{
Logical value, indicating whether the retention levels must be scaled by 
length of the introns/exons.}

\item{unmapValue}{
The value to assign to unmapped rows (i.e. introns/exons).}

\item{FUN}{
The function to apply.}

\item{...}{
Other parameter settings from \code{aggregate()} function.}

}
\value{
The returned value is a data frame if \code{replaceValues} is FALSE and it is 
SummarizedExperiment if \code{replaceValues} is TRUE.
}
\seealso{
	\code{\link{readInterestResults}}
	\code{\link{interest}}
	\code{\link{interest.sequential}}
}

\author{Ali Oghabian}

\examples{

mdsChr22Obj

tmp<- applyOverlap(
	query=mdsChr22Obj,
	subject=mdsChr22Obj,
	type="equal",
	replaceValues=FALSE,
	intExCol="int_ex", 
	intronExon="intron",
	sujectGeneNamesCol="collapsed_transcripts",
	scaleFragment=TRUE,
	scaleLength=TRUE,
	unmapValue=0,
	FUN=head,
	n=1
)

}





