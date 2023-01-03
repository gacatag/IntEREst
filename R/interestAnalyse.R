interestAnalyse <-
function(
	reference,
	bamFile,
	isPaired,
	yieldSize,
	maxNoMappedReads,
	logFile,
	method,
	strandSpecific,
	appendLogFile=TRUE,
	repeatsTableToFilter,
	referenceIntronExon,
	clusterNo,
	junctionReadsOnly,
	isPairedDuplicate,
	isSingleReadDuplicate,
	bpparam,
	limitRanges,
	excludeFusionReads,
	loadLimitRangesReads,
	...)
{
	resTmpPair<-c()
	resTmpSingle<- c()
#Paralle running impelementation
	if(logFile!=""){
		cat( "InERESt:interestAnalyse: Begins ...\n", file=logFile, 
			append=appendLogFile)
	}
	cat( "InERESt:interestAnalyse: Begins ...\n")

	if(is(reference,"GRanges")){
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
	if((!missing(bpparam)))
		BiocParallel::register(bpparam)
	
	if((missing(bpparam))){
	  parallel=FALSE
	} else if(bpparam@class[[1]]=="SerialParam"){
	  parallel=FALSE
	} else {
	  parallel=TRUE
	}
	
	if(isPaired){
		# Defining connection to bam file
		bf<- Rsamtools::BamFile(bamFile, yieldSize=yieldSize, 
				asMates=TRUE, ... )
		# Analyzing mapped paired reads together
		if(length(limitRanges)==0 | (!loadLimitRangesReads)){
		  scParam=Rsamtools::ScanBamParam(
		  	what=Rsamtools::scanBamWhat()[c(1,
			  	3,5,8,13,9, 10, 6, 4, 14, 15)], 
			  flag=Rsamtools::scanBamFlag(isPaired=TRUE,
				  isDuplicate=isPairedDuplicate)) } else {
				    scParam=Rsamtools::ScanBamParam(
				      which=limitRanges,
				      what=Rsamtools::scanBamWhat()[c(1,
				                                      3,5,8,13,9, 10, 6, 4, 14, 15)], 
				      flag=Rsamtools::scanBamFlag(isPaired=TRUE,
				                                  isDuplicate=isPairedDuplicate))				    
				  
				}

		# Initialize the iterator and combine with REDUCE:
		# ITER <- bamIterPair(bf, scParam=scParam)
	#	resTmpPair<- BiocParallel::bpiterate(ITER, interestIntExAnalysePair, 
		YIELD <- function(X, ...) {
		  
		  yld=GenomicAlignments::readGAlignmentPairs(X, 
		                                             param=scParam)
		  return(yld)
		  
		}
		
		if(strandSpecific=="unstranded"){
  		resPair<- GenomicFiles::reduceByYield(bf, YIELD= YIELD, 
  		  MAP=interestIntExAnalysePairUnstranded,
  			reference=reference,
  			maxNoMappedReads=maxNoMappedReads,
  			logFile=logFile,
  			method=method,
  			appendLogFile=appendLogFile,
  			repeatsTableToFilter=repeatsTableToFilter,
  			referenceIntronExon=referenceIntronExon,
  			junctionReadsOnly=junctionReadsOnly,
  			limitRanges=limitRanges,
  			excludeFusionReads=excludeFusionReads,
  			#BPPARAM=bpparam, 
  			REDUCE = `+`, parallel=parallel, iterate = TRUE, 
  			init=rep(0, nrow(reference)))
  		} else{
  		resPair<- GenomicFiles::reduceByYield(bf, YIELD= YIELD, 
  			MAP=interestIntExAnalysePairStranded,
  			reference=reference,
  			maxNoMappedReads=maxNoMappedReads,
  			logFile=logFile,
  			method=method,
  			appendLogFile=appendLogFile,
  			repeatsTableToFilter=repeatsTableToFilter,
  			referenceIntronExon=referenceIntronExon,
  			junctionReadsOnly=junctionReadsOnly,
  			limitRanges=limitRanges,
  			strandSpecific=strandSpecific,
  			excludeFusionReads=excludeFusionReads,
  			#BPPARAM=bpparam, 
  			REDUCE = `+`, parallel=parallel, iterate = TRUE, 
  			init=rep(0, nrow(reference)))  			  
  			}

		
		if(logFile!="")
		  cat( "BETWEEN SINGLE AND PAIRED N1 \n", 
		       file=logFile, append=appendLogFile)
		cat( "BETWEEN SINGLE AND PAIRED N1\n")
		
		# Analyzing single mapped reads
		if(length(limitRanges)==0 | (!loadLimitRangesReads)){
		  scParam2=Rsamtools::ScanBamParam(
			  what=Rsamtools::scanBamWhat()[c(1,
				  3,5,8,13,9, 10, 6, 4, 14, 15,2)], 
			  flag=Rsamtools::scanBamFlag(hasUnmappedMate=TRUE,
				  isPaired=TRUE, 
				  isDuplicate=isSingleReadDuplicate))
		} else{
		  scParam2=Rsamtools::ScanBamParam(
		    which=limitRanges,
		    what=Rsamtools::scanBamWhat()[c(1,
		                                    3,5,8,13,9, 10, 6, 4, 14, 15,2)],
		    flag=Rsamtools::scanBamFlag(hasUnmappedMate=TRUE,
		                                isPaired=TRUE,
		                                isDuplicate=isSingleReadDuplicate))		  
		}
		
		if(logFile!="")
		  cat( "BETWEEN SINGLE AND PAIRED N2 \n", 
		       file=logFile, append=appendLogFile)
		cat( "BETWEEN SINGLE AND PAIRED N2 \n")	
		
		#ITER <- bamIterSingle(bf, scParam=scParam)
		
		if(logFile!="")
		  cat( "BETWEEN SINGLE AND PAIRED N3 \n", 
		       file=logFile, append=appendLogFile)
		cat( "BETWEEN SINGLE AND PAIRED N3\n")
		
		#resTmpSingle<- BiocParallel::bpiterate(ITER,
		YIELD2 <- function(X, ...) {
		  
		  yld=GenomicAlignments::readGAlignmentsList(X, 
		                                             param=scParam2)
		  return(yld)
		  
		}
		if(strandSpecific=="unstranded"){
  		resSingle<- GenomicFiles::reduceByYield(bf, YIELD= YIELD2,
  		  MAP=interestIntExAnalyseSingleUnstranded, 
  			reference=reference,
  			maxNoMappedReads=maxNoMappedReads,
  			logFile=logFile,
  			method=method,
  			appendLogFile=appendLogFile,
  			repeatsTableToFilter=repeatsTableToFilter,
  			referenceIntronExon=referenceIntronExon,
  			junctionReadsOnly=junctionReadsOnly,
  			limitRanges=limitRanges,
  			excludeFusionReads=excludeFusionReads,
  			# BPPARAM=bpparam,
  			REDUCE = `+`, parallel=parallel, iterate = TRUE, 
  			init=rep(0, nrow(reference)))
		} else{
		  resSingle<- GenomicFiles::reduceByYield(bf, YIELD= YIELD2,
		    MAP=interestIntExAnalyseSingleStranded, 
		    reference=reference,
		    maxNoMappedReads=maxNoMappedReads,
		    logFile=logFile,
		    method=method,
		    appendLogFile=appendLogFile,
		    repeatsTableToFilter=repeatsTableToFilter,
		    referenceIntronExon=referenceIntronExon,
		    junctionReadsOnly=junctionReadsOnly,
		    limitRanges=limitRanges,
		    strandSpecific=strandSpecific,
		    excludeFusionReads=excludeFusionReads,
		    # BPPARAM=bpparam,
		    REDUCE = `+`, parallel=parallel, iterate = TRUE, 
		    init=rep(0, nrow(reference)))		  
		}
		
		if(logFile!="")
		  cat( "BETWEEN SINGLE AND PAIRED N4 \n", 
		       file=logFile, append=appendLogFile)
		cat( "BETWEEN SINGLE AND PAIRED N4\n")
	} else {
		# Defining connection to bam file
		bf<- Rsamtools::BamFile(bamFile, yieldSize=yieldSize, ...)
		#Analyzing unpaired sequencing data
#
		if(length(limitRanges)==0 | (!loadLimitRangesReads)){
		scParam3=Rsamtools::ScanBamParam(
			what=Rsamtools::scanBamWhat()[c(1,
				3,5,8,13,9, 10, 6, 4, 14, 15,2)], 
			flag=Rsamtools::scanBamFlag(
				isPaired=NA, 
				isDuplicate=isSingleReadDuplicate))
		} else{
		  scParam3=Rsamtools::ScanBamParam(
		    which=limitRanges,
		    what=Rsamtools::scanBamWhat()[c(1,
		                                    3,5,8,13,9, 10, 6, 4, 14, 15,2)],
		    flag=Rsamtools::scanBamFlag(
		      isPaired=NA, 
		      isDuplicate=isSingleReadDuplicate))		  
		  
		}

		#ITER <- bamIterSingle(bf, scParam=scParam)
		#resTmpSingle<- BiocParallel::bpiterate(ITER, 
		YIELD3 <- function(X, ...) {
		  
		  yld=GenomicAlignments::readGAlignmentsList(X, 
		                                             param=scParam3)
		  return(yld)
		  
		}
		if(strandSpecific=="unstranded"){
  		resSingle<-GenomicFiles::reduceByYield(bf, 
        YIELD= YIELD3,
  		  MAP=interestIntExAnalyseSingleUnstranded, 
  			reference=reference,
  			maxNoMappedReads=maxNoMappedReads,
  			logFile=logFile,
  			method=method,
  			appendLogFile=appendLogFile,
  			repeatsTableToFilter=repeatsTableToFilter,
  			referenceIntronExon=referenceIntronExon,
  			junctionReadsOnly=junctionReadsOnly,
  			limitRanges=limitRanges,
  			excludeFusionReads=excludeFusionReads,
  #			BPPARAM=bpparam, 
  			REDUCE = `+`, parallel=parallel, iterate = TRUE, 
        init=rep(0, nrow(reference)))
		} else{
		  resSingle<-GenomicFiles::reduceByYield(bf, 
		    YIELD= YIELD3,
		    MAP=interestIntExAnalyseSingleStranded, 
		    reference=reference,
		    maxNoMappedReads=maxNoMappedReads,
		    logFile=logFile,
		    method=method,
		    appendLogFile=appendLogFile,
		    repeatsTableToFilter=repeatsTableToFilter,
		    referenceIntronExon=referenceIntronExon,
		    junctionReadsOnly=junctionReadsOnly,
		    limitRanges=limitRanges,
		    strandSpecific=strandSpecific,
		    excludeFusionReads=excludeFusionReads,
		    #			BPPARAM=bpparam, 
		    REDUCE = `+`, parallel=parallel, iterate = TRUE, 
		    init=rep(0, nrow(reference)))		  
		}
  	resPair<-rep(0, nrow(reference)*length(method))
	}
	# resPair<-rep(0, nrow(reference)*length(method))
	# resSingle<-rep(0, nrow(reference)*length(method))

	
	# if(isPaired & (length(resTmpPair)>0)){
	# 	resTmpPair[sapply(resTmpPair, is.null)] <- NULL
	# 	resPair<-Reduce("+", resTmpPair)
	# }
	# if(length(resTmpSingle)>0){
	# 	resTmpSingle[sapply(resTmpSingle, is.null)] <- NULL
	# 	resSingle<-Reduce("+", resTmpSingle)
	# }
	if(is.null(resPair) | (length(resPair)==0)){
	  resPair<- rep(0, 
	    length(which(method%in%c("ExEx", "IntRet", "IntSpan", "ExSkip"))))
	  if(logFile!="")
	    cat( "There are no paired reads !\n", 
	         file=logFile, append=appendLogFile)
	  cat( "There are no paired reads !\n")
	}
	if(is.null(resSingle) | (length(resSingle)==0)){
	  resSingle<- rep(0, 
	    length(which(method%in%c("ExEx", "IntRet", "IntSpan", "ExSkip"))))
	  if(logFile!="")
	    cat( "There are no singletons !\n", 
	         file=logFile, append=appendLogFile)
	  cat( "There are no singletons !\n")
	}
	res<- resPair+resSingle
  #NEW!
	rm("resPair", "resSingle")
	time2=Sys.time()
	runTime=difftime(time2,time1, units="secs")
	if(logFile!="")
		cat( "InERESt:interestAnalyse: Read counting ends. Running time: ",
			runTime," secs\n", file=logFile, append=TRUE)
	cat( "InERESt:interestAnalyse: Read counting ends. Running time: ",
		runTime," secs\n")

	return(res)	
}
