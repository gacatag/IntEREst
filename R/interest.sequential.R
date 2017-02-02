interest.sequential <-
function(
	bamFileYieldSize=1000000,
	bamFile,
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
	method=c("IntRet","ExEx"),
	appendLogFile=FALSE,
	sampleName=c(), 
	scaleLength= c(TRUE,FALSE), 
	scaleFragment= c(TRUE,TRUE)){

	time1=Sys.time()

	if(logFile!=""){
		cat( "Log info: Running interest in sequential mode.\n", file=logFile, 
			append=appendLogFile)
		cat(
"InERESt: Running bamPreprocess. Detailed log info are written in: ",
			logFile ,"\n")
	}
	cat( "Log info: Running interest in sequential mode.\n")		
	if(as.character(class(reference))=="GRanges"){
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
			yieldSize=bamFileYieldSize,
			maxNoMappedReads=1,
			appendLogFile=TRUE,
			logFile=logFile,
			method=method,
			repeatsTableToFilter=repeatsTableToFilter,
			referenceIntronExon=referenceIntronExon,
			junctionReadsOnly=junctionReadsOnly,
			isPairedDuplicate=isPairedDuplicate, 
			isSingleReadDuplicate=isSingleReadDuplicate)




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
		resObj=InterestResult(resultFiles=outFile, 
			readFreq=matrix(tmpDat[,(ncol(tmpDat)-1)], ncol=1), 
			scaledRetention=matrix(tmpDat[,ncol(tmpDat)], ncol=1), 
				sampleNames=sampleName, 	
			scaleLength=scaleLength, scaleFragment=scaleFragment, 
				sampleAnnotation=data.frame(), 
				interestDf=tmpDat[, 1:(ncol(tmpDat)-2)])
	} else if (returnObj & length(method)==2){
		tmpDat<- read.table(outFile, header=TRUE, stringsAsFactors=FALSE)
		resObj= list(
			IntRet=InterestResult(resultFiles=outFile, 
				readFreq=matrix(tmpDat[,(ncol(reference)+1)], ncol=1), 
				scaledRetention=matrix(tmpDat[,(ncol(reference)+2)], ncol=1), 
				sampleNames=sampleName, 
				scaleLength=scaleLength[method="IntRet"], 
				scaleFragment=scaleFragment[method="IntRet"], 
				sampleAnnotation=data.frame(), 
				interestDf=tmpDat[, 1:ncol(reference)]), 
			ExEx=InterestResult(resultFiles=outFile, 
				readFreq=matrix(tmpDat[,(ncol(reference)+3)], ncol=1), 
				scaledRetention=matrix(tmpDat[,(ncol(reference)+4)], ncol=1), 
				sampleNames=sampleName, 
				scaleLength=scaleLength[method="ExEx"], 
				scaleFragment=scaleFragment[method="ExEx"], 
				sampleAnnotation=data.frame(), 
				interestDf=tmpDat[, 1:ncol(reference)]) )
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
