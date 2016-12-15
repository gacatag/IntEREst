\name{u12BoxplotNb}
\alias{u12BoxplotNb}
\title{
	boxplot U12 introns retention levels (or flanking exons junction levels) and (up/down)stream U2 introns (or exons junction levels)
}
\description{
	boxplot U12 introns and (Up/Down)stream U2 introns in \code{interestResult} objects.
}
\usage{
u12BoxplotNb(x, sampleAnnoCol=2, intExCol="int_ex",
	intTypeCol="int_type", intronExon, strandCol="strand", geneIdCol, 
	col=c(), names=c(), lasNames=1, outline=FALSE, plotLegend=TRUE, 
	cexLegend=1, xLegend="topright", yLegend=NULL, legend=c() , ...)
}
\arguments{
\item{x}{ 
	Object of type \code{interestResult}.
}
\item{sampleAnnoCol}{
	Which colummn of \code{@sampleAnnotation} data frame in \code{x} (the input object of \code{interestResult} type) to consider for plotting.
}
\item{intExCol}{
	Column name (or number) that represents whether each row is "intron" or "exon" in \code{@interestDf} of the \code{interestResult} object.
}
\item{intTypeCol}{
	Column name (or number) that represents what type of intron each row of \code{@interestDf} of \code{x} represents. The values in the column are either "U12", "U2" or \code{NA} (for the exons).
}
\item{intronExon}{
	Whether plot intron retention (set \code{intronExon="intron"}) or exon-exon junction (set \code{intronExon="exon"}) levels. 
}
\item{strandCol}{
	Column name (or number) that represents the strand of each row of \code{@interestDf} in \code{x}. The values in the column are either "+", "-" or "*".
}
\item{geneIdCol}{
	Column name (or number) that represents the gene ID of each row of \code{@interestDf} in \code{x}.
}
\item{col}{
	Vector containing box colours. It is either of size 1 or the same size as the number of boxes resulted based on the grouping of the samples defined by \code{sampleAnnoCol}.
}
\item{names}{
	Names to write under group of boxes.
}
\item{lasNames}{
	Orientation of the box names.
}
\item{outline}{
	If outline is TRUE the outlier points are drawn otherwise if FALSE (default) they are not.
}
\item{plotLegend}{
	Whether show legend (TRUE by default).
}
\item{cexLegend}{
	Size of the text in legend .
}
\item{xLegend, yLegend}{
	Position of legend in the plot. For more info see \code{x} and \code{y} parameters in \code{\link{legend}}.
}
\item{legend}{ 
	The replacement texts to be used in legend.
}
\item{...}{
	Other arguments to pass to the \code{boxplot()} function.
}
}
\value{
A code{interestResult} object.
}

\seealso{
	\code{\link{u12Boxplot}}
}

\author{Ali Oghabian}

\examples{


geneId<- paste("gene", c(rep(1,7), rep(2,7), rep(3,7), rep(4,7)), 
	sep="_")
readCnt1<- sample(1:100, 28)
readCnt2<- sample(1:100, 28)
readCnt3<- sample(1:100, 28)
readCnt4<- sample(1:100, 28)
fpkm1<- readCnt1/(tapply(readCnt1, geneId, sum))[geneId]
fpkm2<- readCnt2/(tapply(readCnt2, geneId, sum))[geneId]
fpkm3<- readCnt3/(tapply(readCnt3, geneId, sum))[geneId]
fpkm4<- readCnt4/(tapply(readCnt4, geneId, sum))[geneId]

# Creating object using test data
interestDat<- data.frame( 
		int_ex=rep(c(rep(c("exon","intron"),3),"exon"),4),
		int_ex_num= rep(c(1,1,2,2,3,3,4),4),         
		int_type=rep(c(NA,"U2",NA,"U12",NA,"U2",NA),4),
		strand=rep("*",28),
		gene_id= geneId,
		sam1_readCnt=readCnt1,
		sam2_readCnt=readCnt2,
		sam3_readCnt=readCnt3,
		sam4_readCnt=readCnt4,
		sam1_fpkm=readCnt1,
		sam2_fpkm=readCnt2,
		sam3_fpkm=readCnt3,
		sam4_fpkm=readCnt4
)

interestResultObj<- interestResult(
	resultFiles=paste("file",1:4, sep="_"),
	interestDf= interestDat,
	readFreqColIndex=grep("_readCnt$",colnames(interestDat)),
	scaledRetentionColIndex=grep("_fpkm$",colnames(interestDat)),
	sampleNames=paste("sam",1:4, sep=""),
	scaleLength=TRUE, 
	scaleFragment=FALSE,
	sampleAnnotation=data.frame(
		sampleName=paste("sam",1:4, sep=""),
		gender=c("M","M","F","F"),
		health=c("healthy","unhealthy","healthy","unhealthy")
	)
)

interestDatPosStr<- interestDat
interestDatPosStr[,"strand"]<- rep("+",28)
interestResultPosStrObj<- interestResult(
	resultFiles=paste("file",1:4, sep="_"),
	interestDf= interestDatPosStr,
	readFreqColIndex=grep("_readCnt$",colnames(interestDat)),
	scaledRetentionColIndex=grep("_fpkm$",colnames(interestDat)),
	sampleNames=paste("sam",1:4, sep=""),
	scaleLength=TRUE, 
	scaleFragment=FALSE,
	sampleAnnotation=data.frame(
		sampleName=paste("sam",1:4, sep=""),
		gender=c("M","M","F","F"),
		health=c("healthy","unhealthy","healthy","unhealthy")
	)
)
# plot data where all introns/exons have unknown "*" strand
u12BoxplotNb(interestResultObj, sampleAnnoCol=2, intExCol="int_ex",
intTypeCol="int_type", intronExon="intron", boxplotNames=c(), 
lasNames=1, outline=FALSE, plotLegend=TRUE, geneIdCol="gene_id",
xLegend="topleft", col=c("pink", "lightblue"))

# plot data with all introns/exons featuring "+" strand
u12BoxplotNb(interestResultPosStrObj, sampleAnnoCol=2, 
intExCol="int_ex", intTypeCol="int_type",intronExon="intron",  
boxplotNames=c(), lasNames=1, outline=FALSE, plotLegend=TRUE, 
geneIdCol="gene_id", xLegend="topleft", 
col=c("pink", "lightblue"))

}

