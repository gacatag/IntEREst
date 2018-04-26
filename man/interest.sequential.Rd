\name{interest.sequential}
\alias{interest.sequential}
\title{Wrapup function: Sequential running}
\description{
A read summarization function that countsns all the reads mapping to the 
introns/exons based on the users detailed parameter settings. The process 
runs on a single computing core.}
\usage{
interest.sequential( bamFileYieldSize=1000000, bamFile, isPaired,
	isPairedDuplicate=FALSE, isSingleReadDuplicate=NA, 
	reference, referenceGeneNames, 
	referenceIntronExon, repeatsTableToFilter=c(), 
	junctionReadsOnly=FALSE, outFile, logFile="", 
	returnObj= FALSE, method=c("IntRet", "ExEx", "IntSpan"), 
	appendLogFile=FALSE, sampleName="",
	scaleLength= c(TRUE,FALSE), scaleFragment= c(TRUE,TRUE))
}
\arguments{
\item{bamFileYieldSize}{ 
Maximum number of paired Reads in the temprorary files created as the result of
dividing the input .bam file. }
\item{bamFile}{ Path of the input bam file.}
\item{isPaired}{Whether the bam file is the result of a paired end sequencing
read mapping (TRUE) or not (FALSE).}
\item{isPairedDuplicate}{
Whether extract only (if set TRUE), filter (FALSE) or include (if set NA) PCR 
dupplicates for paired mapped reads. It uses the FLAG field in the bam file to 
filter the duplicate read. If the mapping software does not support detection 
and flaging the duplicate reads \code{dedup} tool of \code{BamUtil} or 
\code{MarkDuplicates} of \code{Picard} tools could be used.}
\item{isSingleReadDuplicate}{
Whether extract only (if set TRUE), filter (FALSE) or include (if set NA) PCR 
dupplicates for single mapped reads. }
\item{reference}{
Dataframe to be used as reference; It should at least contain three same-size 
vectors with the tag names \code{chr}, \code{begin}, and \code{end} which 
describe the genome coordinates of the introns and exons. It also accepts a 
\code{GRanges} object as input. To build a new reference check the 
\code{\link{referencePrepare}} function. }
\item{referenceGeneNames}{ 
A vector with the same size as the row-size of the reference which include the
gene names.}
\item{referenceIntronExon}{ 
A vector with the same size as the row-size of the reference with values 
"intron" and "exon" describing which (intron or exon) each row of the reference
represents.}
\item{repeatsTableToFilter}{
A data frame with similar structure as the \code{reference}, i.e. includes 
\code{chr}, \code{begin}, and \code{end} columns. If defined, all reads mapped 
to the described regions would be ingnored and the Intron/exon lengths would be
corrected to exclude the regions with repetetive DNA sequences. See 
\code{\link{getRepeatTable}}.}
\item{junctionReadsOnly}{
The parameter is considered if the \code{method} is set as \code{IntRet} or 
\code{ExEx} (NOT \code{IntSpan}). It declares whether only consider the 
Intron-Exon or Exon-Exon junction reads and ignore the reads that fully map to 
exons or introns. By default this argument is set as FALSE.}
\item{outFile}{ The name or path of the result file.}
\item{logFile}{
The log file path; if defined log information are written to the log file. }
\item{returnObj}{
If set \code{TRUE} in addition to producing result text files, the results 
would also be returned as an object of class \code{SummarizedExperiment}.}
\item{method}{
A vector describing the summarization methods to use; i.e. whether count reads 
mapping to the introns (\code{IntRet}), reads mapping to the exons 
(\code{ExEx}), or reads spanning the introns (\code{IntSpan}). In 
\code{IntSpan} mode the introns in the reference are taken into account only; 
while in \code{IntRet} the introns and their spanning exons, and in \code{ExEx}
mode only the exons in the reference are taken into account.}
\item{appendLogFile}{
Whether log information should be appended to the \code{logFile}. It is 
\code{FALSE} by default.}
\item{sampleName}{
The name of the sample being analyzed. It will be included in the returned 
object if \code{returnObj} is \code{TRUE}.}
\item{scaleLength}{
A vector constructed of TRUE/FALSE values, same size as the 
\code{method} argument. It indcates whether the retention levels of the 
intron/exons should be scaled to their lengths.}
\item{scaleFragment}{
A vector constructed of TRUE/FALSE values, same size as the 
\code{method} argument. It indcates whether the retention levels of the 
intron/exons should be scaled to the sum of retention levels (i.e. mapped 
fragments) over the genes.}
}
\value{
If \code{returnObj} is set \code{TRUE} in addition to making result text files,
dependant on whether a single or two \code{method} is defined, the results
would be returned as a single object of class \code{SummarizedExperiment} or as 
a list of size 2 which includes 2 objects of class \code{SummarizedExperiment} 
one for IntRet and the other for ExEx. 
}
\seealso{
  \code{\link{interest}}.
}

\author{Ali Oghabian}

\examples{

# Creating temp directory to store the results
outDir<- file.path(tempdir(),"interestFolder")
dir.create(outDir)
outDir<- normalizePath(outDir)

# Loading suitable bam file
bamF <- system.file("extdata", "small_test_SRR1691637_ZRSR2Mut_RHBDD3.bam", 
	package="IntEREst", mustWork=TRUE)


# Choosing reference for the gene RHBDD3
ref=u12[u12[,"gene_name"]=="RHBDD3",]

test= interest.sequential(
	bamFileYieldSize=10000,
	bamFile=bamF,
	isPaired=TRUE,
	isPairedDuplicate=FALSE,
	isSingleReadDuplicate=NA,
	reference=ref,
	referenceGeneNames=ref[,"ens_gene_id"],
	referenceIntronExon=ref[,"int_ex"],
	repeatsTableToFilter=c(),
	outFile=paste(outDir,
	    "interestRes.tsv", sep="/"),
	logFile=paste(outDir, 
	    "log.txt", sep="/"),
	method=c("IntRet","IntSpan"),
	returnObj=TRUE, 
	scaleLength= c(TRUE,FALSE), 
	scaleFragment= c(TRUE,TRUE)
)

test
}

