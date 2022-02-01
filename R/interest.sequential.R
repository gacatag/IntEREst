interest.sequential <-
function(
	bamFileYieldSize=1000000,
	bamFile,
	isPaired,
	isPairedDuplicate=FALSE,
	isSingleReadDuplicate=NA,
	reference,
	referenceGeneNames,
	referenceIntronExon,
	repeatsTableToFilter=c(),
	junctionReadsOnly=FALSE,
	outFile,
	logFile="",
	returnObj=FALSE,
	method=c("IntRet","ExEx","IntSpan"),
	appendLogFile=FALSE,
	sampleName=c(), 
	scaleLength= c(TRUE,FALSE), 
	scaleFragment= c(TRUE,TRUE),
	limitRanges=GRanges(),
	... ){
	method=unique(method)
	if(length(which( ! method %in% c("IntRet", "ExEx", "IntSpan")))!=0)
		stop(paste("Unknown method:", 
			method[! method %in% c("IntRet", "ExEx", "IntSpan")], sep=" "))
	time1=Sys.time()

	if(logFile!=""){
		cat( "Log info: Running interest in sequential mode.\n", file=logFile, 
			append=appendLogFile)
		cat(
"InERESt: Running bamPreprocess. Detailed log info are written in: ",
			logFile ,"\n")
	}
	cat( "Log info: Running interest in sequential mode.\n")		
	if(is(reference,"GRanges")){
		if(length(names(reference))>0){
			tmpReference=data.frame(
				chr=as.character(GenomicRanges::seqnames(reference)), 
				begin=as.numeric(GenomicRanges::start(reference)),
				end=as.numeric(GenomicRanges::end(reference)), 
				strand=as.character(GenomicRanges::strand(reference)), 
				names=names(reference))
		} else{
			tmpReference=data.frame(
				chr=as.character(GenomicRanges::seqnames(reference)), 
				begin=as.numeric(GenomicRanges::start(reference)),
				end=as.numeric(GenomicRanges::end(reference)), 
				strand=as.character(GenomicRanges::strand(reference)))
		}
		reference=tmpReference
	}


	if(logFile!="")
		cat( "InERESt: Running interestAnalyse.sequential.\n", file=logFile, 
			append=TRUE)
	cat( "InERESt: Running interestAnalyse.sequential.\n", append=TRUE)

		inAnRes<- interestAnalyse.sequential(
			reference=reference,
			bamFile=bamFile,
			isPaired=isPaired,
			yieldSize=bamFileYieldSize,
			maxNoMappedReads=1,
			appendLogFile=TRUE,
			logFile=logFile,
			method=method,
			repeatsTableToFilter=repeatsTableToFilter,
			referenceIntronExon=referenceIntronExon,
			junctionReadsOnly=junctionReadsOnly,
			isPairedDuplicate=isPairedDuplicate, 
			isSingleReadDuplicate=isSingleReadDuplicate,
			limitRanges=limitRanges,
			...)




	if(logFile!="")
		cat( "InERESt: Running interestSummarise.\n", file=logFile, 
			append=TRUE)
	cat("InERESt: Running interestSummarise.\n", append=TRUE)

	interestSummarise(
		reference=reference,
		referenceIntronExon=referenceIntronExon,
		inAnRes=inAnRes,
		method=method,
		referenceGeneNames=referenceGeneNames,
		repeatsTableToFilter=repeatsTableToFilter,
		outFile=outFile,
		scaleLength= scaleLength, 
		scaleFragment= scaleFragment
	)


	if(returnObj & length(method)==1){
		tmpDat<- read.table(outFile, header=TRUE, stringsAsFactors=FALSE)
	readFrq<-matrix(tmpDat[,(ncol(tmpDat)-1)], ncol=1)
	scaledRet<-matrix(tmpDat[,ncol(tmpDat)], ncol=1)
	colnames(readFrq)<-c()
	colnames(scaledRet)<-c()
	resObj=InterestResult(
		resultFiles=outFile,
		counts=readFrq,
		scaledRetention=scaledRet,
		rowData=tmpDat[, 1:(ncol(tmpDat)-2)],
		scaleLength=scaleLength, 
		scaleFragment=scaleFragment)



	} else if (returnObj & length(method)>1 & 
		length(which( ! method %in% c("IntRet", "ExEx", "IntSpan")))==0){
		resObj<- list()
		if("IntRet" %in% method){

			tmpDat<- read.table(outFile, header=TRUE, stringsAsFactors=FALSE)
			
			resObj<- c(resObj, list(
				IntRet=InterestResult(resultFiles=outFile, 
					counts=matrix(tmpDat[,(ncol(reference)+1)], ncol=1), 
					scaledRetention=matrix(tmpDat[,(ncol(reference)+2)], ncol=1), 
					scaleLength=scaleLength[method=="IntRet"], 
					scaleFragment=scaleFragment[method=="IntRet"], 
					rowData=tmpDat[, 1:ncol(reference)]) 
			))
		}
		if("ExEx" %in% method){
			indFrq<- 3
			if(!"IntRet"%in%method)
				indFrq<- 1
			tmpDat<- read.table(outFile, header=TRUE, stringsAsFactors=FALSE)
			resObj<- c(resObj, list(
				IntSpan=InterestResult(resultFiles=outFile, 
					counts=matrix(tmpDat[,(ncol(reference)+indFrq)], ncol=1), 
					scaledRetention=
						matrix(tmpDat[,(ncol(reference)+indFrq+1)], ncol=1), 
					scaleLength=scaleLength[method=="ExEx"], 
					scaleFragment=scaleFragment[method=="ExEx"], 
					rowData=tmpDat[, 1:ncol(reference)]) 
			))
		}		
		if("IntSpan" %in% method){
			indFrq<- 3
			if(length(method)==3)
				indFrq<-5
			tmpDat<- read.table(outFile, header=TRUE, stringsAsFactors=FALSE)
			resObj<- c(resObj, list(
				ExEx=InterestResult(resultFiles=outFile, 
					counts=matrix(tmpDat[,(ncol(reference)+indFrq)], ncol=1), 
					scaledRetention=
						matrix(tmpDat[,(ncol(reference)+indFrq+1)], ncol=1), 
					scaleLength=scaleLength[method=="IntSpan"], 
					scaleFragment=scaleFragment[method=="IntSpan"], 
					rowData=tmpDat[, 1:ncol(reference)]) 
			))
		}

	}
	time2=Sys.time()
	runTime=difftime(time2,time1, units="secs")
	if(logFile!="")
		cat( "InERESt: run ends. Full running time: ",runTime," secs\n", 
			file=logFile, append=TRUE)
	cat( "InERESt: run ends. Full running time: ",runTime," secs\n", 
		append=TRUE)
	if(returnObj)
		return(resObj)
}
