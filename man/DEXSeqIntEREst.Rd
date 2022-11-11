\name{DEXSeqIntEREst}
\alias{DEXSeqIntEREst}
\title{
	DEXSeq test for IntEREst object
}
\description{
	Genewise differential exon usage or intron retention test adapted from the 
\code{DEXSeq} package. 
}
\usage{
DEXSeqIntEREst (x, design, reducedModel = ~ sample + intex, fitExpToVar, 
intExCol, geneIdCol, bpparam, silent=TRUE,...)
}
\arguments{
\item{x}{ 
	Object of type \code{SummarizedExperiment}.
}
\item{design}{
Formula specifying the design of the experiment. It must specify an interaction
term between a variable from columns of \code{sampleData(x)} with one of the 
'exon', 'intron' or 'intex' (i.e. intron and exon) variables; based on which of
these variables are used (exon, intron , or 'intex') the x will be filtered
reletively to include exons, introns , or introns and exons. See
\code{\link[DEXSeq]{DEXSeqDataSet}} for more information.}
\item{reducedModel}{
The null model formula. By default it is '~ sample + intex'.}
\item{fitExpToVar}{
A variable name contained in the column data (i.e. column names of
\code{colData(x)}). See \code{\link[DEXSeq]{DEXSeq}} for more information.}
\item{intExCol}{
Column name (or number) that represents whether each row is "intron" or "exon"
 in \code{rowData} of \code{x}.
}
\item{geneIdCol}{
Column name (or number of column) in \code{rowData} of \code{x}, i.e. 
\code{SummarizedExperiment} object, that represents the gene ID of the 
introns and exons in \code{x}.
}
\item{bpparam}{An optional \code{BiocParallelParam} instance defining the 
parallel back-end to be used.}
\item{silent}{
Whether run the DEXSeq function silently (if TRUE) or allow it to print 
messages at each step (if FALSE).
}
\item{...}{
Other parameter settings for the \code{\link[DEXSeq]{DEXSeqDataSet}} function 
in the \code{DEXSeq} package.
}
}
\details{
The \code{design} and \code{reduceModel} accept formula that specify the design
of the experiment. The formula must describe an interaction between variables
from columns of \code{sampleData(x)} with one of the 'exon', 'intron' or 
'intex' (i.e. intron and exon) variables; Based on which of these variables are
used (exon, intron , or 'intex') the input object (\code{x}) will be filtered
reletively to include exons, introns , or introns and exons. Hence the number 
of the rows of the returned value is equal to the number of the rows of the 
filtered object, i.e. the number of the exons, introns or both based on the
\code{design} formula.
}

\value{
 A \code{DEXSeqResults} object.
}
\seealso{
\code{\link{exactTestInterest}}
}

\author{Ali Oghabian}

\examples{

dexseqExRes<-DEXSeqIntEREst (x=mdsChr22ExObj, 
	design= ~ sample + exon + test_ctrl:exon, 
	reducedModel = ~ sample + exon, fitExpToVar="test_ctrl", 
	intExCol="int_ex", geneIdCol="transcripts_id", silent=TRUE)
head(dexseqExRes)

}

