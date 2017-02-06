\name{lfc}
\alias{lfc}
\title{
	Log fold change
}
\description{
Log fold change estimation and normalized log fold change using \code{edgeR}
 package. 
}
\usage{
lfc(x, fcType="edgeR", sampleAnnoCol=c(), sampleAnnotation=c(), 
	silent=TRUE, group=c(), rejection.region="doubletail", 
	pseudoCnt=1, log2=TRUE, ...)
}
\arguments{
\item{x}{ 
	Object of type \code{SummarizedExperiment}.
}
\item{fcType}{
Available as "scaledRetention" or "edgeR" (as default) corresponding to either
log fold change of scaled retention values or degeR normalized log fold change
values.
}
\item{sampleAnnoCol}{
Which colummn of \code{colData} of \code{x} to consider for the analysis.
}
\item{sampleAnnotation}{
A vector of size 2 which cotains values from \code{colData} of 
\code{SummarizedExperiment} object; e.g. if \code{getAnnotation(x)[,
sampleAnnoCol]= c("test", "test", "ctrl","ctrl", ...)} , and the goal is to 
compare "test" and "ctrl" samples, \code{sampleAnnotation} should either be 
\code{c("test", "ctrl")} or \code{c("ctrl", "test")}.
}
\item{silent}{
Whether run \code{\link{exactTestInterest}} silently, without warnings.
}
\item{group}{
Vector to manually define the sample groups (or annotations). It is ignored if
 \code{sampleAnnoCol} is defined.
}
\item{rejection.region}{
The rejection.region parameter in \code{\link[edgeR]{exactTest}}, considered
 only if \code{fcType} is "edgeR".
}
\item{pseudoCnt}{
	Pseudo count for log transformation (default=1).
}
\item{log2}{
Logical value either TRUE (default) or FALSE indicating whether the foldchanges
should be log 2 transformed.
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

geneId<- paste("gene", c(rep(1,5), rep(2,5), rep(3,5), rep(4,5)), 
	sep="_")
readCnt1<- sample(1:100, 20)
readCnt2<- sample(1:100, 20)
readCnt3<- sample(1:100, 20)
readCnt4<- sample(1:100, 20)
fpkm1<- readCnt1/(tapply(readCnt1, geneId, sum))[geneId]
fpkm2<- readCnt2/(tapply(readCnt2, geneId, sum))[geneId]
fpkm3<- readCnt3/(tapply(readCnt3, geneId, sum))[geneId]
fpkm4<- readCnt4/(tapply(readCnt4, geneId, sum))[geneId]

# Creating object using test data
interestDat<- data.frame( 
		int_ex=rep(c(rep(c("exon","intron"),2),"exon"),4),
		int_ex_num= rep(c(1,1,2,2,3),4),         
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

scalRetTmp<- as.matrix(interestDat[ ,scaledRetentionColIndex])
colnames(scalRetTmp)<-gsub("_fpkm$","", colnames(scalRetTmp))

frqTmp<- as.matrix(interestDat[ ,readFreqColIndex])
colnames(frqTmp)<-gsub("_readCnt$","", colnames(frqTmp))


InterestResultObj<- InterestResult(
	resultFiles=paste("file",1:4, sep="_"),
	rowData= interestDat[ , -c(readFreqColIndex, 
		scaledRetentionColIndex)],
	counts= frqTmp,
	scaledRetention= scalRetTmp,
	scaleLength=TRUE, 
	scaleFragment=FALSE,
	sampleAnnotation=data.frame(
		sampleName=paste("sam",1:4, sep=""),
		gender=c("M","M","F","F"),
		health=c("healthy","unhealthy","healthy","unhealthy")
		, row.names=paste("sam", 1:4, sep="")
	)
)

lfcFpkm<- lfc(InterestResultObj, fcType="scaledRetention",
	sampleAnnoCol="health", 
	sampleAnnotation=c("unhealthy", "healthy"), 
	silent=TRUE, group=c(), pseudoFpkm=1, log2=TRUE)
lfcFpkm2<- lfc(InterestResultObj, fcType="scaledRetention", 
	group=c("healthy","unhealthy","healthy","unhealthy"),
	sampleAnnotation=c("unhealthy", "healthy"), 
	silent=TRUE, pseudoFpkm=1, log2=TRUE)

lfcEdgeRFpkm<- lfc(InterestResultObj, fcType="edgeR",
	sampleAnnoCol="health", 
	sampleAnnotation=c("unhealthy", "healthy"), 
	silent=TRUE, group=c(), pseudoFpkm=1, log2=TRUE)

}

