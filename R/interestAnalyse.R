interestAnalyse <-
function(
	reference,
	bamPrerocessRes,
	bamFile,
	yieldSize,
	maxNoMappedReads,
	logFile,
	method,
	appendLogFile=TRUE,
	repeatsTableToFilter,
	referenceIntronExon,
	clusterNo,
	junctionReadsOnly,
	filterPairedDuplicate,
	filterSingleReadDuplicate,
	cl)
{


#Defining function that runs on all parallel cores 
	includePairedDuplicate=NA
	includeSingleReadDuplicate=NA
	if(filterPairedDuplicate & !is.na(filterPairedDuplicate)) 
		includePairedDuplicate=!filterPairedDuplicate
	if(filterSingleReadDuplicate & !is.na(filterSingleReadDuplicate)) 
		includeSingleReadDuplicate=!filterSingleReadDuplicate

#Paralle running impelementation
	if(logFile!=""){
		cat( "InERESt:interestAnalyse: Begins ...\n", file=logFile, 
			append=appendLogFile)
	}
	cat( "InERESt:interestAnalyse: Begins ...\n")

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

	# Set parallel environment
	if(missing(cl)){
		doParallel::registerDoParallel(cores=clusterNo)
	} else {
		doParallel::registerDoParallel(cl, cores=clusterNo)
	}

	i<- 0
	res<- foreach::"%dopar%" (foreach::foreach( i=1:nrow(bamPrerocessRes), 
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
	
	if(!missing(cl))
		on.exit(parallel::stopCluster(cl))

	time2=Sys.time()
	runTime=difftime(time2,time1, units="secs")
	if(logFile!="")
		cat( "InERESt:interestAnalyse: Read counting ends. Running time: ",
			runTime," secs\n", file=logFile, append=TRUE)
	cat( "InERESt:interestAnalyse: Read counting ends. Running time: ",
		runTime," secs\n")

	return(res)	
}
