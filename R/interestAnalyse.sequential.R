interestAnalyse.sequential <-
function(
	reference,
	bamFile,
	isPaired,
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


	# Set parallel environment	
	bpparam <- BiocParallel::SerialParam()

# Initialize the iterator and combine with REDUCE:
	if(isPaired){
		#Defining connecting to bam file
		bf<- Rsamtools::BamFile(bamFile, yieldSize=yieldSize, 
				asMates=TRUE )
		# Analyzing mapped paired reads together
		scParam=Rsamtools::ScanBamParam(
			what=Rsamtools::scanBamWhat()[c(1,
				3,5,8,13,9, 10, 6, 4, 14, 15)], 
			flag=Rsamtools::scanBamFlag(isPaired=TRUE,
				isDuplicate=isPairedDuplicate))
		ITER <- bamIterPair(bf, scParam= scParam)

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

		# Analyzing single mapped reads
		scParam=Rsamtools::ScanBamParam(
			what=Rsamtools::scanBamWhat()[c(1,
				3,5,8,13,9, 10, 6, 4, 14, 15)], 
			flag=Rsamtools::scanBamFlag(hasUnmappedMate=TRUE,
				isPaired=TRUE, 
				isDuplicate=isSingleReadDuplicate))
		ITER <- bamIterSingle(bf, scParam= scParam)
		resTmpSingle<- BiocParallel::bpiterate(ITER, 
			interestIntExAnalyseSingle, 
			reference=reference,
			maxNoMappedReads=maxNoMappedReads,
			logFile=logFile,
			method=method,
			appendLogFile=appendLogFile,
			repeatsTableToFilter=repeatsTableToFilter,
			referenceIntronExon=referenceIntronExon,
			junctionReadsOnly=junctionReadsOnly,
			BPPARAM=bpparam)
	} else {
		#Defining connecting to bam file
		bf<- Rsamtools::BamFile(bamFile, yieldSize=yieldSize)
		#Analyzing unpaired sequencing data
		scParam=Rsamtools::ScanBamParam(
			what=Rsamtools::scanBamWhat()[c(1,
				3,5,8,13,9, 10, 6, 4, 14, 15)], 
			flag=Rsamtools::scanBamFlag(
				isPaired=NA, 
				isDuplicate=isSingleReadDuplicate))

		ITER <- bamIterSingle(bf, scParam= scParam)
		resTmpSingle<- BiocParallel::bpiterate(ITER, 
			interestIntExAnalyseSingle, 
			reference=reference,
			maxNoMappedReads=maxNoMappedReads,
			logFile=logFile,
			method=method,
			appendLogFile=appendLogFile,
			repeatsTableToFilter=repeatsTableToFilter,
			referenceIntronExon=referenceIntronExon,
			junctionReadsOnly=junctionReadsOnly,
			BPPARAM=bpparam)
	}
	
	if(isPaired)
		resPair<-Reduce("+", resTmpPair)

	resSingle<-Reduce("+", resTmpSingle)

	res<- resPair+resSingle

	time2=Sys.time()
	runTime=difftime(time2,time1, units="secs")
	if(logFile!="")
		cat( "InERESt:interestAnalyse: Read counting ends. Running time: ",
			runTime," secs\n", file=logFile, append=TRUE)
	cat( "InERESt:interestAnalyse: Read counting ends. Running time: ",
		runTime," secs\n")

	return(res)	
}
