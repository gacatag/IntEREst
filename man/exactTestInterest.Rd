\name{exactTestInterest}
\alias{exactTestInterest}
\title{
	Exact test
}
\description{
	Compute genewise exact test between two groups of read counts, using the \code{edgeR} package. 
}
\usage{
exactTestInterest(x, sampleAnnoCol=c(), sampleAnnotation=c(), 
	geneIdCol, silent=TRUE, group=c(), rejection.region="doubletail", 
	big.count=900, prior.count=0.125, disp="common", ...)
}
\arguments{
\item{x}{ 
	Object of type \code{interestResult}.
}
\item{sampleAnnoCol}{
	Which colummn of \code{@sampleAnnotation} of object \code{interestResult} to consider for the analysis.
}
\item{sampleAnnotation}{
	A vector of size 2 which cotains values from \code{@sampleAnnotation} of \code{interestResult} object; e.g. if \code{x@sampleAnnotation [,sampleAnnoCol] = c("test", "test", "ctrl","ctrl", ...)}, and the goal is to compare "test" and "ctrl" samples, \code{sampleAnnotation} should either be \code{c("test","ctrl")} or \code{c("ctrl","test")}.
}
\item{geneIdCol}{
	Column name (or number) that represents the gene ID of each row of \code{@sampleAnnotation} in object \code{interestResult}.
}
\item{silent}{
	Whether run silently, i.e. without printing the top differential expression tags.
}
\item{group}{
	Vector to manually define the sample groups (or annotations). It is ignored if \code{sampleAnnopCol} is defined.
}
\item{rejection.region}{
	The \code{rejection.region} parameter in \code{\link[edgeR]{exactTest}} from \code{edgeR} package.
}
\item{big.count}{
	The \code{big.count} parameter in \code{\link[edgeR]{exactTest}} from \code{edgeR} package.
}
\item{prior.count}{
	The \code{prior.count} parameter in \code{\link[edgeR]{exactTest}} from \code{edgeR} package.
}
\item{disp}{
	The type of estimating the dispersion in the data. Available options are: "tagwise", "trended", "common" and "genewise". It is also possible to assign a number for manually setting the \code{disp}.
}
\item{...}{
	Other parameter settings for the \code{\link[edgeR]{estimateDisp}} function (e.g. the \code{design} parameter) in the \code{edgeR} package.
}

}
\value{
\item{table}{ Data  frame  containing  columns  for  the  log2 fold-change (logFC),  the  average of log2 counts-per-million (logCPM), and the two-sided p-value (PValue).  }
\item{comparison}{The name of the two compared groups.}
\item{dispersionType}{The name of the type of dispersion used.}
\item{dispersion}{The estimated dispersion values.}
}
\seealso{
	\code{\link{lfc}}, \code{\link{glmInterest}}, \code{\link{qlfInterest}}, \code{\link{treatInterest}}
}

\author{Ali Oghabian}

\examples{


geneId<- paste("gene", c(rep(1,7), rep(2,7), rep(3,7), rep(4,7)), 
	sep="_")

y <- matrix(rnbinom(28*4,size=1/0.2,mu=10),nrow=28,ncol=4)
readCnt1<- y[,1]
readCnt2<- y[,2]
readCnt3<- y[,3]
readCnt4<- y[,4]
# Scale to reads fragments mapped to the genes 
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

interestResObj<- interestResult(
	resultFiles=paste("file",1:4, sep="_"),
	interestDf= interestDat[, -c(readFreqColIndex, 
		scaledRetentionColIndex)],
	readFreq=as.matrix(interestDat[, readFreqColIndex]),
	scaledRetention=as.matrix(interestDat[, 
		scaledRetentionColIndex]),
	sampleNames=paste("sam",1:4, sep=""),
	scaleLength=FALSE, scaleFragment=TRUE, 
	sampleAnnotation=data.frame(
		sampleName=paste("sam",1:4, sep=""),
		gender=c("M","M","F","F"),
		health=c("healthy","unhealthy","healthy","unhealthy")
	)
)

res<- exactTestInterest(interestResObj, sampleAnnoCol="gender", 
	sampleAnnotation=c("F","M"), geneIdCol= "gene_id", 
	silent=TRUE, disp="common")

}

