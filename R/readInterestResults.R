readInterestResults<-function(resultFiles, sampleNames, sampleAnnotation, 
	commonColumns, freqCol, scaledRetentionCol, scaleLength, scaleFragment, 
	reScale=FALSE, geneIdCol, repeatsTableToFilter=c()){	
	readFreqColIndex=c()
	scaledRetentionColIndex=c()

	for(i in 1:length(resultFiles)){
		cat("\nReading file", i, "/", length(resultFiles), ":",  
			resultFiles[i], "\n", sep=" ")
		tmpDat=read.table(resultFiles[i], header=TRUE, stringsAsFactors=FALSE,
			sep='\t')
		if(i==1){
			ref=tmpDat[,commonColumns]
			res=ref
		}
		if(freqCol<=ncol(tmpDat))
			res=cbind(res, as.numeric(tmpDat[,freqCol]))
		if(scaledRetentionCol<=ncol(tmpDat))
			res=cbind(res, as.numeric(tmpDat[,scaledRetentionCol]))
		readFreqColIndex=c(readFreqColIndex, (ncol(res)-1))
		scaledRetentionColIndex=c(scaledRetentionColIndex, ncol(res))
		colnames(res)[(ncol(res)-1):ncol(res)]=paste(sampleNames[i], 
			c("frequency", "scaledRetention"), sep="_")
		if(reScale){
			FPKM=res[,ncol(res)-1]

			if(scaleFragment){
				referenceGeneNamesTmp=as.character(ref[ , geneIdCol])
				geneCnt=tapply(FPKM, referenceGeneNamesTmp, sum)
				FPKM=((10^6)*FPKM)/
					(as.numeric(geneCnt[referenceGeneNamesTmp])+1)
			}
			if(scaleLength){
				if(i==1)
					lenRef=correctLengthRepeat(ref, repeatsTableToFilter)
				FPKM=((10^3)*FPKM)/lenRef

			}
			res[,ncol(res)]=FPKM
		}
	}
	if(!missing(sampleAnnotation)){
		if(is.vector(sampleAnnotation))
			sampleAnnotation=data.frame(sampleGroup=sampleAnnotation)
		rownames(sampleAnnotation)=sampleNames

		readFrq<-as.matrix(res[,readFreqColIndex])
		scaledRet<-as.matrix(res[,scaledRetentionColIndex])
		colnames(readFrq)<-sampleNames
		colnames(scaledRet)<-sampleNames
		resObj=InterestResult(
			resultFiles=resultFiles,
			counts=readFrq,
			scaledRetention=scaledRet,
			sampleAnnotation=sampleAnnotation,
			rowData=ref,
			scaleLength=scaleLength, 
			scaleFragment=scaleFragment)
	} else {

		readFrq<-as.matrix(res[,readFreqColIndex])
		scaledRet<-as.matrix(res[,scaledRetentionColIndex])
		colnames(readFrq)<-c()
		colnames(scaledRet)<-c()
		resObj=InterestResult(
			resultFiles=resultFiles,
			counts=readFrq,
			scaledRetention=scaledRet,
			rowData=ref,
			scaleLength=scaleLength, 
			scaleFragment=scaleFragment)

	}

	return(resObj)
}
