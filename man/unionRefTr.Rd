\name{unionRefTr}
\alias{unionRefTr}
\title{Union introns/exons of transcripts}
\description{ 
Performs union on the overlapping introns/exons so that the final merged 
transcripts would feature from each exon or intron, one copy.}
\usage{
unionRefTr( referenceChr, referenceBegin, referenceEnd, referenceTr, 
	referenceIntronExon, intronExon="exon", silent=FALSE)
}
\arguments{
\item{referenceChr}{
Chromosome names of the references (e.g. introns).}
\item{referenceBegin}{
A vector that corresponds to the begin coordinates of the reference.}
\item{referenceEnd}{
A vector that corresponds to the end coordinates of the reference.}
\item{referenceTr}{
A character vector that includes transcription IDs.}
\item{referenceIntronExon}{
A vector with the same size as the \code{referenceChr}, \code{referenceBegin} 
and \code{referenceEnd} which contains \code{'intron'} and \code{'exon'}
describing what (either intron or exon) each element of the 3 vectors 
represents.}
\item{intronExon}{
Should be assigned either \code{'intron'} or \code{'exon'} or 
\code{c('intron','exon')} based on whether match the PWM to the intronic, 
exonic, or intronic and exonic regions of the reference. By default it seeks 
matches in intronic regions (intronExon='intron').}
\item{silent}{
Whether run silently.
}
}
\value{
Data frame containing merged transcripts structure. The merged transcripts 
feature from each intron or exon, one copy ONLY.
}
\seealso{
  \code{\link{annotateU12}}.
}

\author{Ali Oghabian}

\examples{

unU12Ex<-unionRefTr( referenceChr=u12[1:94,"chr"], 
	referenceBegin=u12[1:94,"begin"], referenceEnd=u12[1:94,"end"], 
	referenceTr=u12[1:94,"trans_name"], 
	referenceIntronExon=u12[1:94,"int_ex"], intronExon="exon", silent=TRUE)

unU12Int<-unionRefTr( referenceChr=u12[1:94,"chr"], 
	referenceBegin=u12[1:94,"begin"], referenceEnd=u12[1:94,"end"], 
	referenceTr=u12[1:94,"trans_name"], 
	referenceIntronExon=u12[1:94,"int_ex"], intronExon="intron", silent=TRUE)

unU12IntEx<-unionRefTr( referenceChr=u12[1:94,"chr"], 
	referenceBegin=u12[1:94,"begin"], referenceEnd=u12[1:94,"end"], 
	referenceTr=u12[1:94,"trans_name"], 
	referenceIntronExon=u12[1:94,"int_ex"], intronExon=c("intron","exon"), 
	silent=TRUE)

}





