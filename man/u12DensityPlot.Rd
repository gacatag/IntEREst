\name{u12DensityPlot}
\alias{u12DensityPlotIntron}
\alias{u12DensityPlot}
\title{
Density plot of fld changes of intron retention and exon-exon junction levels
}
\description{
Density plot of fold change of the retention levels of U12- vs U2- type intron,
or exon-exon junction levels of the flanking exons. For the density plot of the
foldchange of intron retention levels the \code{u12DensityPlotIntron()} 
function or \code{u12DensityPlot()} function with \code{intronExon= "intron"} 
can be used. For density plot of the foldchange of exon-exon junction levels 
use \code{u12DensityPlot()} function with \code{intronExon= "exon"}.
}
\usage{

u12DensityPlot(x, 
	type=c("U12", "U2Up", "U2Dn", "U2UpDn", "U2Rand"), 
	fcType="edgeR", sampleAnnotation=c(),  sampleAnnoCol=c(), 
	group=c(), intExCol="int_ex", intTypeCol="int_type", intronExon, 
	strandCol="strand", geneIdCol="collapsed_transcripts", 
	naUnstrand=FALSE, col=1, lty=1, lwd=1, plotLegend=TRUE, 
	cexLegend=1, xLegend="topright", yLegend=NULL, legend=c(), 
	randomSeed=NULL, xlab="", ...)

u12DensityPlotIntron(x, 
	type= c("U12", "U2Up", "U2Dn", "U2UpDn", "U2Rand"), 
	fcType= "edgeR", sampleAnnotation=c(),  sampleAnnoCol=c(), 
	group=c(), intExCol="int_ex", intTypeCol="int_type", 
	strandCol= "strand", geneIdCol= "collapsed_transcripts", 
	naUnstrand=FALSE, col=1, lty=1, lwd=1, plotLegend=TRUE, 
	cexLegend=1, xLegend="topright", yLegend=NULL, legend=c(), 
	randomSeed=NULL, xlab="", ...)
}
\arguments{
\item{x}{ 
	Object of type \code{SummarizedExperiment}.
}
\item{type}{
A vector that includes the type of introns to plot. Available options are U12 
introns "U12", U2 introns at downstream of U12 introns "U2Dn", U2 introns at 
upstream of U12 introns "U2Up",  U2 introns at upstream or downstream of U12 
introns suitable for when the coorduinates in object x are unstranded (their 
strand is "*") "U2UpDn",  random U2 introns from object x "U2Rand". Settings 
"U2Up", "U2Dn" and "U2UpDn" are useful only if the refernce is linearly 
ordered. References with exons only resulted by \code{referencePrepare} and 
\code{unionRefTr} are NOT necessarily linearly ordered.
}
\item{fcType}{
Available as "fpkm" or "edgeR" (as default) corresponding to either log fold 
change of fpkm values or degeR normalized log fold change values.
}
\item{sampleAnnoCol}{
Which colummn of \code{colData} of \code{x} to consider for plotting.
}
\item{sampleAnnotation}{
A vector of size 2 which cotains values from \code{colData} of 
\code{SummarizedExperiment} object; e.g. if \code{getAnnotation(x)[,
sampleAnnoCol]= c("test", "test", "ctrl","ctrl", ...)} , and the goal is to 
compare "test" and "ctrl" samples, \code{sampleAnnotation} should either be 
\code{c("test", "ctrl")} or \code{c("ctrl", "test")}.
}
\item{group}{
Vector to manually define the sample groups (or annotations). It is ignored if 
\code{sampleAnnoCol} is defined.
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
\item{naUnstrand}{
Replace unstranded results, i.e. introns or exon with "*" strand, with NA (to 
be excluded).
}
\item{col}{
A vector with the size of 1 or the same size as the \code{type} parameter which 
includes the colour/colours of the plotted density lines (default=1).
}
\item{lty}{
A vector with the size of 1 or the same size as the \code{type} parameter which
includes the type of the plotted density lines (default=1).
}
\item{lwd}{
A vector with the size of 1 or the same size as the \code{type} parameter which
includes the width of the plotted density lines (default=1).
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
\item{randomSeed}{
	Seed value for random number generator.
}
\item{xlab}{
The lable of the X axis of the plot; by default it is "".
}
\item{...}{
Other parameter settings from the \code{\link[graphics]{plot}} function.
}

}
\value{
Returns NULL.
}

\seealso{
\code{\link{exactTestInterest}}, \code{\link{lfc}}
}

\author{Ali Oghabian}

\examples{

u12DensityPlotIntron(mdsChr22Obj, 
    type= c("U12", "U2Up", "U2Dn", "U2UpDn", "U2Rand"), 
    fcType= "edgeR", sampleAnnoCol="test_ctrl", 
    sampleAnnotation=c("ctrl","test"), intExCol="int_ex", 
    intTypeCol="intron_type", strandCol= "strand", 
    geneIdCol= "collapsed_transcripts_id", naUnstrand=FALSE, col=c(2,3,4,5,6), 
    lty=c(1,2,3,4,5), lwd=1, plotLegend=TRUE, cexLegend=0.7, 
    xLegend="topright", yLegend=NULL, legend=c(), randomSeed=10,
    ylim=c(0,0.6), xlab=expression("log"[2]*" fold change FPKM"))

}

