bamPreprocess <- function(
	yieldSize=1000000,
	bamFile,
	logFile="",
	filterPairedDuplicate=TRUE, 
	filterSingleReadDuplicate=FALSE,
	appendLogFile=TRUE)
{

	if(logFile!="")
		cat( "InERESt:bamPreprocess: Begins ...\n", file=logFile, 
			append=appendLogFile)
	cat( "InERESt:bamPreprocess: Begins ...\n")
	includePairedDuplicate=NA
	includeSingleReadDuplicate=NA
	if(filterPairedDuplicate & !is.na(filterPairedDuplicate)) 
		includePairedDuplicate=!filterPairedDuplicate
	if(filterSingleReadDuplicate & !is.na(filterSingleReadDuplicate)) 
		includeSingleReadDuplicate=!filterSingleReadDuplicate

	bf<-Rsamtools::BamFile(bamFile, yieldSize=yieldSize, asMates=TRUE )
	count=0
	if(logFile!="")
		cat( "InERESt: Reading paired mapped reads from bam file.\n", 
			file=logFile, append=appendLogFile)
	cat( "InERESt: Reading paired mapped reads from bam file.\n")
	open(bf)
	readSizeChk=2*yieldSize
	datReadNo<-c()
	isPaired<-c()
	while(readSizeChk==(2*yieldSize)) {
		count=count+1; 

		scParam=Rsamtools::ScanBamParam(
			what=Rsamtools::scanBamWhat()[c(1,3,5,8,13,9, 10, 6, 4, 14, 15)], 
			flag=Rsamtools::scanBamFlag(isPaired=TRUE,
				isDuplicate=includePairedDuplicate))

		readTmp=GenomicAlignments::readGAlignmentPairs(bf, param=scParam)
		datReadNo<-c(datReadNo, count)
		isPaired<-c(isPaired, TRUE)
		readSizeChk=2*length(readTmp)
		
	}
	close(bf)

	bf<- Rsamtools::BamFile(bamFile, yieldSize=yieldSize, asMates=TRUE)
	open(bf)
	readSizeChk=yieldSize 
	count=0
	if(logFile!="")
		cat( "InERESt: Reading single mapped reads from bam file.\n", 
			file=logFile, append=appendLogFile)
	cat( "InERESt: Reading single mapped reads from bam file.\n")
	while(readSizeChk==yieldSize) {
		count=count+1; 
		scParam=Rsamtools::ScanBamParam(
			what=Rsamtools::scanBamWhat()[c(1,3,5,8,13,9, 10, 6, 4,14,15)], 
			flag=Rsamtools::scanBamFlag(hasUnmappedMate=TRUE, isPaired=TRUE, 
				isDuplicate=includeSingleReadDuplicate))

		readTmp= GenomicAlignments::readGAlignmentsList(bf, param=scParam)
		datReadNo<-c(datReadNo, count)
		isPaired<-c(isPaired, FALSE)
		readSizeChk=length(readTmp)
		lenTmp<- unlist(sapply(readTmp, length))
	} 
	close(bf)
	return(data.frame(datReadNo, isPaired))
}
