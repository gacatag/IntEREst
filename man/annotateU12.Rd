\name{annotateU12}
\alias{annotateU12}
\title{Annotate the U12 (and U2) type introns}
\description{ 
Receives coordinates, a reference genome and PWMs of splice site of U12 and U2 
type introns, and returns a data.frame with 2 columns. The first column shows 
wheather the corresponding sequences matches U12, U2 or both (U12/U2) consensus
sequences (based on their score when fitting the PWMs). The second column shows
whether the match is on positive strand or negative when fitting the PWMs to
the sequences.  }
\usage{
annotateU12(pwmU12U2=c(), pwmSsIndex=c(), referenceChr, referenceBegin, 
	referenceEnd, referenceIntronExon, intronExon='intron', 
	matchWindowRelativeUpstreamPos=c() , matchWindowRelativeDownstreamPos=c(), 
	minMatchScore='80\%', refGenome='', setNaAs='U2', annotateU12Subtype=TRUE)
}
\arguments{
\item{pwmU12U2}{
A list containing position weight matrices of (in order): Donor site, branch 
point, and acceptor site of U12-type introns, and donor site and acceptor site 
of U2-type introns. If not provided, the information related to \code{pwmU12db}
data is used.}

\item{pwmSsIndex}{
A list (or vector) that contains the column number in each element of 
\code{pwmU12U2} that represents the 5' or 3' Splice Site; The order should be 
equivalent to the \code{pwmU12U2}. If not provided the information from 
\code{pwmU12db} data is used, i.e. \code{pwmSsIndex=list(indexDonU12=1, 
indexBpU12=1, indexAccU12=3, indexDonU2=1, indexAccU2=3)} }

\item{referenceChr}{
Chromosome names of the references (e.g. introns).}

\item{referenceBegin}{
A vector that corresponds to the begin coordinates of the reference 
(e.g. introns).}

\item{referenceEnd}{
A vector that corresponds to the end coordinates of the reference (e.g. 
introns). \code{referenceEnd} should be greater than or equal to 
\code{referenceBegin}.}

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

\item{matchWindowRelativeUpstreamPos}{
A vector the same size as the \code{pwmU12U2} (and the same order of 
donor/acceptor sites' information in \code{pwmU12U2}) which consists of the 
upstream distance from the donor/acceptor site from which each PWM should be 
tested (to see if they match). If not provided, the information from 
\code{pwmU12db} data is used i.e. 
\code{matchWindowRelativeUpstreamPos= c(NA, -29, NA, NA, NA)}.}

\item{matchWindowRelativeDownstreamPos}{
A vector the same size as the \code{pwmU12U2} (and the same order of 
donor/acceptor sites' information in \code{pwmU12U2}) which consists of the 
downstream distance from the donor/acceptor site to which each PWM should be 
tested (to see if they match). If not provided, the information from 
\code{pwmU12db} data is used i.e. 
\code{matchWindowRelativeDownstreamPos= c(NA,-9, NA, NA, NA)}.}

\item{minMatchScore}{
Min percentage match score, when scoring matching of a sequence to \code{pwm}. 
Different score thresholds could also be defined for the various sites (U12/U2 
donors, the U12 branch point and U12/U2 acceptors); A vector with 5 elements
can be assigned which each shows the match score to use for each PWM in 
\code{pwmU12U2}.}

\item{refGenome}{
The reference genome; Object of class BSgenome. Use \code{available.genome()} 
from the \code{BSgenome} package to see the available genomes.}

\item{setNaAs}{
Defines that if reference (e.g. intron) did not match any of U12 or U2 type 
introns based on the scores obtained from PWM what should the function return. 
If an intron was not proven to be U12 or U2 based on PWM scores it can be 
considered as U2-type since the U12 type introns constitute for about 1\% of 
introns in human genome and they are muxh more conserved than the U2 type 
introns, hence the default is 'U2'; otherwise it is also possible to set it as 
NA or nan or 'U12/U2'.}

\item{annotateU12Subtype}{Whether annotate the subtypes of the U12 type 
Introns. The value is \code{TRUE} by default.}
}
\value{
 Data frame containing 3 columns representing (in order): intron type (U12, U2 
or none), strand match indicating whether the PWM matches to the sequence 
(+ strand) or the reverese complement of the sequence (- strand) or none (NA), 
and the U12 subtype (GT-AG or AT-AC).
}
\seealso{
  \code{\link{buildSsTypePwms}}.
}

\author{Ali Oghabian}

\examples{



# Improting genome
BSgenome.Hsapiens.UCSC.hg19 <- 
	BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#Choosing subset of rows
ind<- 69:94
# Annotate U12 introns with strong U12 donor site, branch point 
# and acceptor site from the u12 data in the package
annoU12<- 
	annotateU12(pwmU12U2=list(pwmU12db[[1]][,11:17],pwmU12db[[2]]
		,pwmU12db[[3]][,38:40],pwmU12db[[4]][,11:17],
		pwmU12db[[5]][,38:40]), 
	pwmSsIndex=list(indexDonU12=1, indexBpU12=1, indexAccU12=3, 
		indexDonU2=1, indexAccU2=3), 
	referenceChr=u12[ind,'chr'], 
	referenceBegin=u12[ind,'begin'], 
	referenceEnd=u12[ind,'end'], 
	referenceIntronExon=u12[ind,"int_ex"], 
	intronExon="intron",  
	matchWindowRelativeUpstreamPos=c(NA,-29,NA,NA,NA),
	matchWindowRelativeDownstreamPos=c(NA,-9,NA,NA,NA), 
	minMatchScore=c(rep(paste(80,"\%",sep=""),2), "60\%", 
		paste(80,"\%",sep=""), "60\%"), 
	refGenome=BSgenome.Hsapiens.UCSC.hg19, 
	setNaAs="U2", 
	annotateU12Subtype=TRUE)

# How many U12 and U2 type introns with strong U12 donor sites, 
# acceptor sites (and branch points for U12-type) are there?
table(annoU12[,1])

}





