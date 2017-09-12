\name{deseqInterest}
\alias{deseqInterest}
\title{
	DESeq2 analysis for IntEREst object
}
\description{
	Differential intron retention test adapted from the \code{DESeq2} package. 
}
\usage{
deseqInterest (x, design, pAdjustMethod = "BH", 
	sizeFactor=c(), contrast, bpparam, ...)
}
\arguments{
\item{x}{ 
	Object of type \code{SummarizedExperiment}.
}
\item{design}{
Formula specifying the design of the experiment. It must specify an interaction
term between variables from column names of \code{sampleData(x)}.}
\item{pAdjustMethod}{
What adjustment method to be sed on the p-values. See \code{\link{p.adjust}} 
for more information.}
\item{sizeFactor}{
Numeric vector with the same size as the clolumn size of the count matrix in 
\code{x}, if defined it will be used for scaling of the count matrix.}
\item{contrast}{
Argument speciftying the comparison to extract from \code{x}. See 
\code{\link[DESeq2]{results}} function in the \code{DESeq2} package for more
information.
}
\item{bpparam}{An optional \code{BiocParallelParam} instance defining the 
parallel back-end to be used. If not defined the function will run 
sequentially (on a single computing core).}
\item{...}{
Other parameter settings for the \code{\link[DESeq2]{results}} function 
in the \code{DESeq2} package.
}
}

\value{
 a \code{DESeqResults} object.
}
\seealso{
\code{\link{exactTestInterest}}
\code{\link{qlfInterest}}, 
\code{\link{treatInterest}}
\code{\link{DEXSeqIntEREst}}
}

\author{Ali Oghabian}

\examples{

deseqRes<- deseqInterest(x=mdsChr22Obj, 
	design=~test_ctrl, contrast=list("test_ctrl_test_vs_ctrl"))

# Number of U12/U2 type significantly differential retained introns in chr22
table(rowData(mdsChr22Obj)[which(deseqRes$padj<.01), "int_type"])
}

