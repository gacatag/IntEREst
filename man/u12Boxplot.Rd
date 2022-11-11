\name{u12Boxplot}
\alias{u12Boxplot}
\title{
	U12 boxplot
}
\description{
A boxplot method for U12 and U2-type introns of \code{SummarizedExperiment} 
objects.
}
\usage{
u12Boxplot(x, sampleAnnoCol=NA, intExCol="int_ex",
	intTypeCol="int_type", intronExon, col="white",
	boxplotNames=c(), lasNames=3, outline=FALSE, addGrid=FALSE, ...)
}
\arguments{
\item{x}{ 
Object of type \code{SummarizedExperiment}.
}
\item{sampleAnnoCol}{
Which colummn of \code{colData} in \code{x} to consider for plotting.
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
\item{col}{
Vector showing box colours. It is either of size 1 or the same size as the 
number of groups to be plotted.
}
\item{boxplotNames}{
Names to write under boxes. If not defined, as names, it pastes U12/U2 (intron 
annotation) to the sample group annotations separated by a space " ".
}
\item{lasNames}{
Orientation of the box names.
}
\item{outline}{
If outline is TRUE the outlier points are drawn otherwise if FALSE (default) 
they are not.
}
\item{addGrid}{ 
Whether add a grid under the boxplots (FALSE by default). 
}
\item{...}{
Other arguments to pass to the \code{boxplot()} function.
}
}
\value{
A \code{SummarizedExperiment} object.
}

\seealso{
	\code{\link{u12BoxplotNb}}
}

\author{Ali Oghabian}

\examples{


u12Boxplot(mdsChr22Obj, sampleAnnoCol="type", 
    intExCol="int_ex",  intTypeCol="intron_type", intronExon="intron", 
    col=rep(c("orange", "yellow"),3) ,  lasNames=3, 
    outline=FALSE, ylab="FPKM", cex.axis=0.8)

}

