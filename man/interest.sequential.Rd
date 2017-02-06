\name{interest.sequential}
\alias{interest.sequential}
\title{Wrapup function: Sequential running}
\description{
A wrapper function that runs all necessary functions for the intron/exon 
retention analysis. }
\usage{
interest.sequential( bamFileYieldSize=1000000, bamFile, 
	isPairedDuplicate=FALSE, isSingleReadDuplicate=NA, 
	reference, referenceGeneNames, 
	referenceIntronExon, repeatsTableToFilter=c(), 
	junctionReadsOnly=FALSE, outFile, logFile="", 
	returnObj= FALSE, method=c("IntRet","ExEx"), 
	appendLogFile=FALSE, sampleName="",
	scaleLength= c(TRUE,FALSE), scaleFragment= c(TRUE,TRUE))
}
\arguments{
\item{bamFileYieldSize}{ 
Maximum number of paired Reads in the temprorary files created as the result of
dividing the input .bam file. }
\item{bamFile}{ Path of the input bam file.}
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
Whether only consider the Intron-Exon or Exon-Exon junction reads and ignore
the reads that fully map to exons or introns. By default this argument is set 
as FALSE.}
\item{outFile}{ The name or path of the result file.}
\item{logFile}{
The log file path; if defined log information are written to the log file. }
\item{returnObj}{
If set \code{TRUE} in addition to producing result text files, the results 
would also be returned as an object of class \code{SummarizedExperiment}.}
\item{method}{
A vector describing the methods to use; i.e. whether estimate Intron retention 
(\code{IntRet}) or Exon-Exon junction (\code{ExEx}) levels or both. By default 
it returns both.}
\item{appendLogFile}{
Whether log information should be appended to the \code{logFile}. It is 
\code{FALSE} by default.}
\item{sampleName}{
The name of the sample being analyzed. It will be included in the returned 
object if \code{returnObj} is \code{TRUE}.}
\item{scaleLength}{
A vector constructed of TRUE/FALSE values, same size as the 
\code{method argument}. It indcates whether the retention levels of the 
intron/exons should be scaled to their lengths.}
\item{scaleFragment}{
A vector constructed of TRUE/FALSE values, same size as the 
\code{method argument}. It indcates whether the retention levels of the 
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
# Time demanding function
\dontrun{
# Creating temp directory to store the results
outDir<- file.path(tempdir(),"interestFolder")
dir.create(outDir)
outDir<- normalizePath(outDir)

# Loading suitable bam file
library("MDS.Chr22.U12Genes")

# Analyse the first bam file
i<- 1

dir.create(paste(outDir, names(MDS_Chr22_BAMFILES)[i], sep="/"))
interest.sequential(
	bamFileYieldSize=1000000,
	tmpDir=paste(outDir, names(MDS_Chr22_BAMFILES)[i], "tmp", 
		sep="/"),
	bamFile=MDS_Chr22_BAMFILES[i],
	filterPairedDuplicate=TRUE,
	filterSingleReadDuplicate=FALSE,
	reference=u12,
	referenceGeneNames=u12[,"gene_ens_id"],
	referenceIntronExon=u12[,"int_ex"],
	outFile=paste(outDir, names(MDS_Chr22_BAMFILES)[i], 
		"interestRes.tsv", sep="/"),
	logFile=paste(outDir, names(MDS_Chr22_BAMFILES)[i], 
		"log.txt", sep="/"),
	method=c("IntRet","ExEx"),
	scaleLength= c(TRUE,FALSE), 
	scaleFragment= c(TRUE,TRUE)
)
}
}

