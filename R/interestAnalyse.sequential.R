interestAnalyse.sequential <-
function(
	reference=reference,
	bamPrerocessRes,
	bamFile=bamFile,
	yieldSize,
	maxNoMappedReads,
	appendLogFile=TRUE,
	logFile=logFile,
	method=method,
	repeatsTableToFilter,
	referenceIntronExon,
	junctionReadsOnly,
	filterPairedDuplicate, 
	filterSingleReadDuplicate)
{

	includePairedDuplicate=NA
	includeSingleReadDuplicate=NA
	if(filterPairedDuplicate & !is.na(filterPairedDuplicate)) 
		includePairedDuplicate=!filterPairedDuplicate
	if(filterSingleReadDuplicate & !is.na(filterSingleReadDuplicate)) 
		includeSingleReadDuplicate=!filterSingleReadDuplicate

#Sequential running impelementation
	if(logFile!=""){
		cat( "InERESt:interestAnalyse.sequential: Begins ...\n", file=logFile, 
			append=appendLogFile)
	}
	cat( "InERESt:interestAnalyse.sequential: Begins ...\n")

	if(as.character(class(reference))=="GRanges"){
		if(length(names(reference))>0){
			tmpReference=data.frame(
				chr=as.character(GenomicRanges::seqnames(reference)), 
				begin=as.numeric(GenomicRanges::start(reference)),
				end=as.numeric(GenomicRanges::end(reference)), 
				strand=as.character(GenomicRanges::strand(reference)), 
				names=as.character(names(reference)))
		} else{
			tmpReference=data.frame(
				chr=as.character(GenomicRanges::seqnames(reference)), 
				begin=as.numeric(GenomicRanges::start(reference)),
				end=GenomicRanges::end(reference), 
				strand=as.character(GenomicRanges::strand(reference)))
		}
		reference=tmpReference
	}


	time1=Sys.time()

	# Run sequential command
	i<- 0
	res<- foreach::"%do%" (foreach::foreach( i=1:nrow(bamPrerocessRes), 
		.combine='+'), 
		interestIntExAnalyse(no=i, 
			reference=reference,
			bamPrerocessRes=bamPrerocessRes,
			bamFile=bamFile,
			yieldSize=yieldSize,
			maxNoMappedReads=maxNoMappedReads,
			logFile=logFile,
			method=method,
			appendLogFile=appendLogFile,
			repeatsTableToFilter=repeatsTableToFilter,
			referenceIntronExon=referenceIntronExon,
			junctionReadsOnly=junctionReadsOnly,
			includePairedDuplicate=includePairedDuplicate,
			includeSingleReadDuplicate=includeSingleReadDuplicate))

	time2=Sys.time()
	runTime=difftime(time2,time1, units="secs")
	if(logFile!="")
		cat( 
"InERESt:interestAnalyse.sequential: Read counting ends. Running time: ",
			runTime," secs\n", file=logFile, append=TRUE)
	cat( 
"InERESt:interestAnalyse.sequential: Read counting ends. Running time: ",
		runTime," secs\n")

	return(res)	
}
