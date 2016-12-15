\name{u12Boxplot}
\alias{u12Boxplot}
\title{
	U12 boxplot
}
\description{
	A boxplot method for U12 and U2-type introns of \code{interestResult} objects.
}
\usage{
u12Boxplot(x, sampleAnnoCol=NA, intExCol="int_ex",
	intTypeCol="int_type", intronExon, col="white",
	boxplotNames=c(), lasNames=3, outline=FALSE, ...)
}
\arguments{
\item{x}{ 
	Object of type \code{interestResult}.
}
\item{sampleAnnoCol}{
	Which colummn of \code{@sampleAnnotation} datat frame in \code{interestResult} object to consider for plotting.
}
\item{intExCol}{
	Column name (or number) that represents whether each row in \code{@sampleAnnotation} of object \code{interestResult} is "intron" or "exon" .
}
\item{intTypeCol}{
	Column name (or number) that represents what type of intron each row of \code{@sampleAnnotation} of object \code{interestResult} represents.
}
\item{intronExon}{
	Whether plot intron retention (set \code{intronExon="intron"}) or exon-exon junction (set \code{intronExon="exon"}) levels. 
}
\item{col}{
	Vector showing box colours. It is either of size 1 or the same size as the number of groups to be plotted.
}
\item{boxplotNames}{
	Names to write under boxes.
}
\item{lasNames}{
	Orientation of the box names.
}
\item{outline}{
	If outline is TRUE the outlier points are drawn otherwise if FALSE (default) they are not.
}
\item{...}{
	Other arguments to pass to the \code{boxplot()} function.
}
}
\value{
A code{interestResult} object.
}

\seealso{
	\code{\link{u12BoxplotNb}}
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

# Box plot and group samples based on the annotations (column 2)
u12Boxplot(interestResultObj, sampleAnnoCol=2, intExCol="int_ex",
	intTypeCol="int_type", intronExon="intron", col="white", 
	boxplotNames=c(), lasNames=3, outline=FALSE)

# Box plot for Individual samples
u12Boxplot(interestResultObj, sampleAnnoCol=NA, intExCol="int_ex",
	intTypeCol="int_type", intronExon="intron", col="white", 
	boxplotNames=c(), lasNames=3, outline=FALSE)
}

