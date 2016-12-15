bamPreprocess <- function(
	yieldSize=1000000,
	outFolder="./tmpSam/",
	bamFile,
	logFile="",
	filterPairedDuplicate=TRUE, 
	filterSingleReadDuplicate=FALSE,
	appendLogFile=TRUE)
{

	BamFile<- Rsamtools::BamFile
	ScanBamParam<- Rsamtools::ScanBamParam	
	scanBam<- Rsamtools::scanBam	
	scanBamFlag<- Rsamtools::scanBamFlag
	scanBamWhat<- Rsamtools::scanBamWhat

	if(logFile!="")
		cat( "InERESt:bamPreprocess: Begins ...\n", file=logFile, append=appendLogFile)
	cat( "InERESt:bamPreprocess: Begins ...\n")
	includePairedDuplicate=NA
	includeSingleReadDuplicate=NA
	if(filterPairedDuplicate & !is.na(filterPairedDuplicate)) includePairedDuplicate=!filterPairedDuplicate
	if(filterSingleReadDuplicate & !is.na(filterSingleReadDuplicate)) includeSingleReadDuplicate=!filterSingleReadDuplicate

	bf=BamFile(bamFile, yieldSize=yieldSize, asMates=TRUE )
	count=0
# Create necessary folders
	dir.create(paste(outFolder,"read1/", sep="/"))
	dir.create(paste(outFolder,"read2/", sep="/"))
	dir.create(paste(outFolder,"single/", sep="/"))

# Open reading and writing the reads to files
	open(bf)
	readSizeChk=2*yieldSize
	if(logFile!="") 
		cat( "InERESt:bamPreprocess: Extracting Paired-mapped reads.\n", file=logFile, append=TRUE)
	cat( "InERESt:bamPreprocess: Extracting Paired-mapped reads.\n")
	while(readSizeChk==(2*yieldSize)) {
		count=count+1; 

		scParam=ScanBamParam(what=scanBamWhat()[c(1,3,5,8,13,9, 10, 6, 4, 14, 15)], flag=scanBamFlag(isPaired=TRUE, isDuplicate=includePairedDuplicate))

		readTmp=scanBam(bf, param=scParam)
		readSizeChk=length(readTmp[[1]][[1]])

		ind=which(readTmp[[1]][[11]]=="mated")
		if(length(ind)>0){
			ind1=ind[seq(from = 1, to = length(ind), by =2 )]
			ind2=ind[seq(from = 2, to = length(ind), by =2 )]
			write.table(cbind(as.character(readTmp[[1]][[1]])[ind1], as.character(readTmp[[1]][[2]])[ind1],as.character(readTmp[[1]][[3]])[ind1],as.character(readTmp[[1]][[4]])[ind1], as.character(readTmp[[1]][[5]])[ind1], as.character(readTmp[[1]][[6]])[ind1],as.character(readTmp[[1]][[7]])[ind1],as.character( readTmp[[1]][[8]])[ind1]), paste(outFolder, "/read1/tmpSam_", count, sep=""),sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
			write.table(cbind(as.character(readTmp[[1]][[1]])[ind2],as.character(readTmp[[1]][[2]])[ind2],as.character(readTmp[[1]][[3]])[ind2],as.character(readTmp[[1]][[4]])[ind2],as.character(readTmp[[1]][[5]])[ind2],as.character(readTmp[[1]][[6]])[ind2],as.character(readTmp[[1]][[7]])[ind2], as.character(readTmp[[1]][[8]])[ind2]), paste(outFolder, "/read2/tmpSam_", count, sep=""),sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE) 
		}
	}
	close(bf)

	bf=BamFile(bamFile, yieldSize=yieldSize, asMates=TRUE)
	open(bf)
	readSizeChk=yieldSize 
	count=0
	if(logFile!="")
		cat( "InERESt:bamPreprocess: Extracting Single mapped reads.\n", file=logFile, append=TRUE)
	cat( "InERESt:bamPreprocess: Extracting Single mapped reads.\n")
	while(readSizeChk==yieldSize) 
		{count=count+1; 
		scParam=ScanBamParam(what=scanBamWhat()[c(1,3,5,8,13,9, 10, 6, 4,14,15)], 
			flag=scanBamFlag(hasUnmappedMate=TRUE, isPaired=TRUE, isDuplicate=includeSingleReadDuplicate))

		readTmp=scanBam(bf, param=scParam)
		readSizeChk=length(readTmp[[1]][[1]])
		ind=which(readTmp[[1]][[11]]!="mated")
		if(length(ind)>0){
			write.table(cbind(as.character(readTmp[[1]][[1]][ind]),as.character(readTmp[[1]][[2]][ind]),as.character(readTmp[[1]][[3]][ind]),as.character(readTmp[[1]][[4]][ind]),as.character(readTmp[[1]][[5]][ind]),as.character(readTmp[[1]][[6]][ind]),as.character(readTmp[[1]][[7]][ind]), as.character(readTmp[[1]][[8]][ind])), paste(outFolder, "/single/tmpSam_", count, sep=""),sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
		} 
	}
	close(bf)
}
