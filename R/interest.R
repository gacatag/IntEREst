interest <-
function(
	bamFileYieldSize=1000000, tmpDir, bamFile,
	filterPairedDuplicate=TRUE,
	filterSingleReadDuplicate=FALSE,
	reference, referenceGeneNames,
	referenceIntronExon, repeatsTableToFilter=c(),
	junctionReadsOnly=FALSE, outFile, logFile="", delTmpFolder=FALSE, 
	returnObj=FALSE, method=c("IntRet","ExEx"),	clusterNo="",
	appendLogFile=FALSE, sampleName="", 
	scaleLength= c(TRUE,FALSE), scaleFragment= c(TRUE,TRUE))
{
	time1=Sys.time()
	suppressWarnings(dir.create(tmpDir))
	suppressWarnings(dir.create(paste(tmpDir,"preprocessResult", sep="/")))
	if(logFile!=""){
		cat( "Log info: Running interest in Parallel mode.\n", file=logFile, append=appendLogFile)
		cat("InERESt: Running bamPreprocess. Detailed log info are written in: ",logFile ,"\n", file=logFile, append=TRUE)
	}


	if(as.character(class(reference))=="GRanges"){
		if(length(names(reference))>0){
			tmpReference=data.frame(chr=as.character(GenomicRanges::seqnames(reference)), begin=as.numeric(GenomicRanges::start(reference)), 
				end=as.numeric(GenomicRanges::end(reference)), strand=as.character(GenomicRanges::strand(reference)), 
				names=as.character(names(reference)))
		} else{
			tmpReference=data.frame(chr=as.character(GenomicRanges::seqnames(reference)), begin=as.numeric(GenomicRanges::start(reference)), 
				end=as.numeric(GenomicRanges::end(reference)), strand=as.character(GenomicRanges::strand(reference)))
		}
		reference=tmpReference
	}


	cat( "Log info: Running interest in Parallel mode.\n")
	cat( "InERESt: Running bamPreprocess.\n")
	bamPreprocess(
		yieldSize=bamFileYieldSize,
		outFolder=paste(tmpDir,"preprocessResult/", sep="/"),
		bamFile=bamFile,
		logFile=logFile,
		filterPairedDuplicate=filterPairedDuplicate, 
		filterSingleReadDuplicate=filterSingleReadDuplicate
	)

	suppressWarnings(dir.create(paste(tmpDir,"interestAnalyseResult", sep="/")))

	if(logFile!="")
		cat( "InERESt: Running interestAnalyse.\n", file=logFile, append=TRUE)
	cat("InERESt: Running interestAnalyse.\n", append=TRUE)
	interestAnalyse(
		reference=reference,
		outDir=paste(tmpDir,"interestAnalyseResult", sep="/"),
		logFile,
		pairFiles=dir(paste(tmpDir,"preprocessResult","read1", sep="/"), full.names=FALSE),
		singleFiles=dir(paste(tmpDir,"preprocessResult", "single", sep="/"), full.names=FALSE),
		inLoc=paste(tmpDir,"preprocessResult", sep="/"),
		method=method,
		repeatsTableToFilter=repeatsTableToFilter,
		referenceIntronExon=referenceIntronExon,
		clusterNo=clusterNo,
		junctionReadsOnly=junctionReadsOnly
	)
	if(logFile!="")
		cat( "InERESt: Running interestSummarise.\n", file=logFile, append=TRUE)
	cat("InERESt: Running interestSummarise.\n", append=TRUE)

	interestSummarise(
		reference=reference,
		referenceIntronExon=referenceIntronExon,
		inLoc=paste(tmpDir,"interestAnalyseResult", method, sep="/"),
		method=method,
		referenceGeneNames=referenceGeneNames,
		repeatsTableToFilter=repeatsTableToFilter,
		outFile=outFile,
		scaleLength= scaleLength, 
		scaleFragment= scaleFragment
	)

	if(delTmpFolder) unlink(tmpDir, recursive = TRUE, force = TRUE)
	tmpDat=read.table(outFile, header=TRUE, stringsAsFactors=FALSE)
	if(returnObj & length(method)==1){
		resObj=InterestResult(resultFiles=outFile, readFreq=matrix(tmpDat[,(ncol(tmpDat)-1)], ncol=1), 
			scaledRetention=matrix(tmpDat[,ncol(tmpDat)], ncol=1), sampleNames=sampleName, 	
			scaleLength=scaleLength, scaleFragment=scaleFragment, sampleAnnotation=data.frame(), interestDf=tmpDat[, 1:(ncol(tmpDat)-2)])
		return(resObj)
	} else if (returnObj & length(method)==2){
		resObj= list(IntRet=InterestResult(resultFiles=outFile, readFreq=matrix(tmpDat[,(ncol(reference)+1)], ncol=1), 
				scaledRetention=matrix(tmpDat[,(ncol(reference)+2)], ncol=1), sampleNames=sampleName, scaleLength=scaleLength[method="IntRet"], 
				scaleFragment=scaleFragment[method="IntRet"], sampleAnnotation=data.frame(), interestDf=tmpDat[, 1:ncol(reference)]), 
			ExEx=InterestResult(resultFiles=outFile, readFreq=matrix(tmpDat[,(ncol(reference)+3)], ncol=1), 
				scaledRetention=matrix(tmpDat[,(ncol(reference)+4)], ncol=1), sampleNames=sampleName, 
				scaleLength=scaleLength[method="ExEx"], scaleFragment=scaleFragment[method="ExEx"], 
				sampleAnnotation=data.frame(), interestDf=tmpDat[, 1:ncol(reference)]) )
	}
	time2=Sys.time()
	runTime=difftime(time2,time1, units="secs")
	if(logFile!="")
		cat( "InERESt: run ends. Full running time: ",runTime," secs\n", file=logFile, append=TRUE)
	cat( "InERESt: run ends. Full running time: ",runTime," secs\n", append=TRUE)
}
