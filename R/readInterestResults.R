readInterestResults<-function(resultFiles, sampleNames, sampleAnnotation, commonColumns, freqCol, scaledRetentionCol, scaleLength, 
	scaleFragment, reScale=FALSE, geneIdCol, repeatsTableToFilter=c()){	
	readFreqColIndex=c()
	scaledRetentionColIndex=c()

	for(i in 1:length(resultFiles)){
		cat("\nReading file", i, "/", length(resultFiles), ":",  resultFiles[i], "\n", sep=" ")
		tmpDat=read.table(resultFiles[i], header=TRUE, stringsAsFactors=FALSE, sep='\t')
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
		colnames(res)[(ncol(res)-1):ncol(res)]=paste(sampleNames[i], c("frequency", "scaledRetention"), sep="_")
		if(reScale){
			FPKM=res[,ncol(res)-1]

			if(scaleFragment){
				referenceGeneNamesTmp=as.character(ref[ , geneIdCol])
				geneCnt=tapply(FPKM, referenceGeneNamesTmp, sum)
				FPKM=((10^6)*FPKM)/(as.numeric(geneCnt[referenceGeneNamesTmp])+1)
			}
			if(scaleLength){
				if(i==1)
					lenRef=correctLengthRepeat(ref, repeatsTableToFilter)
				FPKM=((10^3)*FPKM)/lenRef

			}
			res[,ncol(res)]=FPKM
		}
	}
	if(is.vector(sampleAnnotation))
		sampleAnnotation=data.frame(sampleGroup=sampleAnnotation)
	sampleAnnotation=data.frame(sampleNames=sampleNames,sampleAnnotation)
	resObj=interestResult(resultFiles=resultFiles, readFreqColIndex=readFreqColIndex, scaledRetentionColIndex=scaledRetentionColIndex, 
		sampleNames=sampleNames, scaleLength=scaleLength, scaleFragment=scaleFragment, sampleAnnotation=sampleAnnotation, interestDf=res)
	return(resObj)
}
