\name{lfc}
\alias{lfc}
\title{
	Log fold change
}
\description{
	Log fold change estimation and normalized log fold change using the \code{edgeR} package. 
}
\usage{
lfc(x, fcType="edgeR", sampleAnnoCol=c(), sampleAnnotation=c(), 
	silent=TRUE, group=c(), rejection.region="doubletail", 
	pseudoCnt=1, log2=TRUE, ...)
}
\arguments{
\item{x}{ 
	Object of type \code{interestResult}.
}
\item{fcType}{
	Available as "scaledRetention" or "edgeR" (as default) corresponding to either log fold change of scaled retention values or degeR normalized log fold change values.
}
\item{sampleAnnoCol}{
	Which column of \code{@sampleAnnotation} of object \code{interestResult} to consider for plotting.
}
\item{sampleAnnotation}{
	A vector of size 2 which cotains values from \code{@sampleAnnotation} of \code{interestResult} object; e.g. if x@sampleAnnotation[,sampleAnnoCol] = c("test", "test", "ctrl","ctrl", ...), and the goal is to compare "test" and "ctrl" samples, \code{sampleAnnotation} should either be c("test","ctrl") or c("ctrl","test").
}
\item{silent}{
	Whether run \code{\link{exactTestInterest}} silently, without warnings.
}
\item{group}{
	Vector to manually define the sample groups (or annotations). It is ignored if \code{sampleAnnoCol} is defined.
}
\item{rejection.region}{
	The rejection.region parameter in \code{\link[edgeR]{exactTest}}, considered only if \code{fcType} is "edgeR".
}
\item{pseudoCnt}{
	Pseudo count for log transformation (default=1).
}
\item{log2}{
	Logical value either TRUE (default) or FALSE indicating whether the foldchanges should be log 2 transformed.
}
\item{...}{
	Other parameter settings from the \code{\link{exactTestInterest}} function.
}

}
\value{
Vector including fold change values.
}

\seealso{
	\code{\link{exactTestInterest}}, \code{\link{u12DensityPlotIntron}}
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
interestDat<-data.frame( 
		int_ex=rep(c(rep(c("exon","intron"),3),"exon"),4),
		int_ex_num= rep(c(1,1,2,2,3,3,4),4),         
		int_type=rep(c(NA,"U2",NA,"U12",NA,"U2",NA),4),
		strand=rep("*",28),
		gene_id= geneId,
		sam1_readCnt=readCnt1,
		sam2_readCnt=readCnt2,
		sam3_readCnt=readCnt3,
		sam4_readCnt=readCnt4,
		sam1_fpkm=fpkm1,
		sam2_fpkm=fpkm2,
		sam3_fpkm=fpkm3,
		sam4_fpkm=fpkm4
)
readFreqColIndex<- grep("_readCnt$",colnames(interestDat))
scaledRetentionColIndex<- grep("_fpkm$",colnames(interestDat))
interestResObj<-interestResult(
	resultFiles=paste("file",1:4, sep="_"),
	interestDf= interestDat[ , -c(readFreqColIndex, 
		scaledRetentionColIndex)],
	readFreq=as.matrix(interestDat[ , readFreqColIndex]),
	scaledRetention=as.matrix(interestDat[ , 
		scaledRetentionColIndex]),
	sampleNames=paste("sam",1:4, sep=""),
	scaleLength=TRUE, 
	scaleFragment=FALSE,
	sampleAnnotation=data.frame(
		sampleName=paste("sam",1:4, sep=""),
		gender=c("M","M","F","F"),
		health=c("healthy","unhealthy","healthy","unhealthy")
	)
)

lfcFpkm<- lfc(interestResObj, fcType="scaledRetention",
	sampleAnnoCol="health", 
	sampleAnnotation=c("unhealthy", "healthy"), 
	silent=TRUE, group=c(), pseudoFpkm=1, log2=TRUE)
lfcFpkm2<- lfc(interestResObj, fcType="scaledRetention", 
	group=c("healthy","unhealthy","healthy","unhealthy"),
	sampleAnnotation=c("unhealthy", "healthy"), 
	silent=TRUE, pseudoFpkm=1, log2=TRUE)

lfcEdgeRFpkm<- lfc(interestResObj, fcType="edgeR",
	sampleAnnoCol="health", 
	sampleAnnotation=c("unhealthy", "healthy"), 
	silent=TRUE, group=c(), pseudoFpkm=1, log2=TRUE)

}

