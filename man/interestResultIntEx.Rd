\name{interestResultIntEx}
\alias{interestResultIntEx}
\title{
Building results object that contains Intron-retention and exon-exon junction
information
}
\description{
Building \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} object
from an intron retention and an exon-exon junction results in IntEREst. The 
average of the junction levels are added to the SummerizedExperiment object of
the intron retentions.

}
\usage{

interestResultIntEx (intObj, exObj, intExCol=c(), 
	mean.na.rm=TRUE, postExName="ex_junc" )

}
\arguments{
\item{intObj}{
A \code{SummarizedExperiment} including intron retention information.
}
\item{exObj}{
A \code{SummarizedExperiment} including exon-exon junction information.
}
\item{intExCol}{
Column name (or number) in the \code{rowData} of the intron object that 
represents whether each row of \code{x} assays is "intron" or "exon".
}
\item{mean.na.rm}{
Whether exclude missing values when measuring the mean.
}
\item{postExName}{
The postfix to use for the column names of the exons junction values in the
}
}
\value{
Returns an object of class \code{SummarizedExperiment}. 
}
\seealso{
	\code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
	\code{\link{attributes}}
	\code{\link{addAnnotation}}
	\code{\link{counts-method}}
	\code{\link{plot-method}}
}

\author{Ali Oghabian}

\examples{

testIntObj<- InterestResult(
		resultFiles= paste(paste("testFile",1:3, sep="_"),"bam", sep="."), 
		counts= matrix(1:15, ncol=3, nrow=5, byrow=TRUE, 
			dimnames= list(c(), paste("s", 1:3, sep="_"))), 
		scaledRetention= matrix(1:15, ncol=3, nrow=5, byrow=TRUE, 
			dimnames= list(c(), paste("s", 1:3, sep="_"))), 
		scaleLength= FALSE, 
		scaleFragment= FALSE, 
		sampleAnnotation= data.frame(
			files=paste(paste("testFile",1:3, sep="_"),"bam", sep="."),
			names=paste("s", 1:3, sep="_"),
			row.names=paste("s", 1:3, sep="_")), 
		rowData=data.frame(id= paste("i", 1:5, sep="_"),
			chr= rep("chr1", 5),
			begin=seq(100, by=100, length.out=5 ),
			end=seq(110, by=100, length.out=5 ),
			strand=rep("+",5))
	)

testExObj<- InterestResult(
		resultFiles= paste(paste("testFile",1:3, sep="_"),"bam", sep="."), 
		counts= matrix(1:30, ncol=3, nrow=10, byrow=TRUE, 
			dimnames= list(c(), paste("s", 1:3, sep="_"))), 
		scaledRetention= matrix(1:30, ncol=3, nrow=10, byrow=TRUE, 
			dimnames= list(c(), paste("s", 1:3, sep="_"))), 
		scaleLength= FALSE, 
		scaleFragment= FALSE, 
		sampleAnnotation= data.frame(
			files=paste(paste("testFile",1:3, sep="_"),"bam", sep="."),
			names=paste("s", 1:3, sep="_"),
			row.names=paste("s", 1:3, sep="_")), 
		rowData=data.frame(id= paste("e", 1:10, sep="_"),
			chr= rep("chr1", 10),
			begin= c(seq(90, by=100, length.out=5),
				seq(111, by=100, length.out=5)),
			end= c(seq(99, by=100, length.out=5),
				seq(120, by=100, length.out=5 )),
			strand=rep("+",10))
	)


(testIntExObj<- interestResultIntEx(intObj=testIntObj, exObj=testExObj, 
	mean.na.rm=TRUE, postExName="ex_junc" ) )

}






