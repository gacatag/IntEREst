interestAnalyse.sequential <-
function(
	reference,
	bamFile,
	yieldSize,
	maxNoMappedReads,
	logFile,
	method,
	appendLogFile=TRUE,
	repeatsTableToFilter,
	referenceIntronExon,
	junctionReadsOnly,
	isPairedDuplicate,
	isSingleReadDuplicate)
{

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

	bf<- Rsamtools::BamFile(bamFile, yieldSize=yieldSize, 
			asMates=TRUE )


	# Set parallel environment	
	bpparam <- BiocParallel::SerialParam()

# Initialize the iterator and combine with REDUCE:
	ITER <- bamIterPair(bf, isPairedDuplicate=isPairedDuplicate)

	resTmpPair<- BiocParallel::bpiterate(ITER, interestIntExAnalysePair, 
		reference=reference,
		maxNoMappedReads=maxNoMappedReads,
		logFile=logFile,
		method=method,
		appendLogFile=appendLogFile,
		repeatsTableToFilter=repeatsTableToFilter,
		referenceIntronExon=referenceIntronExon,
		junctionReadsOnly=junctionReadsOnly,
		BPPARAM=bpparam)

	ITER <- bamIterSingle(bf, isSingleReadDuplicate=isSingleReadDuplicate)
	resTmpSingle<- BiocParallel::bpiterate(ITER, interestIntExAnalyseSingle, 
		reference=reference,
		maxNoMappedReads=maxNoMappedReads,
		logFile=logFile,
		method=method,
		appendLogFile=appendLogFile,
		repeatsTableToFilter=repeatsTableToFilter,
		referenceIntronExon=referenceIntronExon,
		junctionReadsOnly=junctionReadsOnly,
		BPPARAM=bpparam)

	res<- Reduce("+", resTmpPair)+Reduce("+", resTmpSingle)

	time2=Sys.time()
	runTime=difftime(time2,time1, units="secs")
	if(logFile!="")
		cat( "InERESt:interestAnalyse: Read counting ends. Running time: ",
			runTime," secs\n", file=logFile, append=TRUE)
	cat( "InERESt:interestAnalyse: Read counting ends. Running time: ",
		runTime," secs\n")

	return(res)	
}
