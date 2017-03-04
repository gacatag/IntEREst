\name{u12BoxplotNb}
\alias{u12BoxplotNb}
\title{
boxplot U12 introns retention levels (or flanking exons junction levels) and 
(up/down)stream U2 introns (or exons junction levels)
}
\description{
boxplot U12 introns and (Up/Down)stream U2 introns in 
\code{SummarizedExperiment} objects.
}
\usage{
u12BoxplotNb(x, sampleAnnoCol=2, intExCol="int_ex",
	intTypeCol="int_type", intronExon, strandCol="strand", geneIdCol, 
	col=c(), names=c(), lasNames=1, outline=FALSE, plotLegend=TRUE, 
	cexLegend=1, xLegend="topright", yLegend=NULL, legend=c() , ...)
}
\arguments{
\item{x}{ 
Object of type \code{SummarizedExperiment}.
}
\item{sampleAnnoCol}{
Which colummn of \code{colData} of \code{x} to consider for plotting.
}
\item{intExCol}{
Column name (or number) that represents whether each row of \code{x} assays is 
"intron" or "exon".
}
\item{intTypeCol}{
Column name (or number) that represents what type of intron each row of 
\code{x} assays represents.
}
\item{intronExon}{
Whether plot intron retention (set \code{intronExon="intron"}) or exon-exon 
junction (set \code{intronExon="exon"}) levels. 
}
\item{strandCol}{
Column name (or number) that represents the strand of each row of assays in
\code{x}. The values in the column are either "+", "-" or "*".
}
\item{geneIdCol}{
Column name (or number) that represents the gene ID of each row of assays in
\code{x}.
}
\item{col}{
Vector containing box colours. It is either of size 1 or the same size as the 
number of boxes resulted based on the grouping of the samples defined by 
\code{sampleAnnoCol}.
}
\item{names}{
Names to write under group of boxes.
}
\item{lasNames}{
Orientation of the box names.
}
\item{outline}{
If outline is TRUE the outlier points are drawn otherwise if FALSE (default) 
they are not.
}
\item{plotLegend}{
Whether show legend (TRUE by default).
}
\item{cexLegend}{
Size of the text in legend .
}
\item{xLegend, yLegend}{
Position of legend in the plot. For more info see \code{x} and \code{y} 
parameters in \code{\link{legend}}.
}
\item{legend}{ 
The replacement texts to be used in legend.
}
\item{...}{
	Other arguments to pass to the \code{boxplot()} function.
}
}
\value{
Returns NULL
}

\seealso{
	\code{\link{u12Boxplot}}
}

\author{Ali Oghabian}

\examples{

u12BoxplotNb(mdsChr22Obj, sampleAnnoCol="type", lasNames=1,
    intExCol="int_ex", intTypeCol="int_type", intronExon="intron", 
    boxplotNames=c(), outline=FALSE, plotLegend=TRUE, 
    geneIdCol="ens_gene_id", xLegend="topleft", 
    col=c("pink", "lightblue", "lightyellow"), ylim=c(0,1e+06), 
    ylab="FPKM", cex.axis=0.8)

}

