mergeInterestResult<-function(x, y){
	if(x@scaleLength!=y@scaleLength | x@scaleFragment!=y@scaleFragment)
		stop("The retention levels should be similarly scaled. The @scaleLength and @scaleFragment of the two input objects should be equal.")
	res=x
	res@resultFiles=c(res@resultFiles, y@resultFiles)
	res@sampleNames=c(res@sampleNames, y@sampleNames)
	res@sampleAnnotation=rbind(res@sampleAnnotation, y@sampleAnnotation)
	res@readFreq=cbind(x@readFreq, y@readFreq)
	res@scaledRetention=cbind(x@scaledRetention, y@scaledRetention)
	xyIntDf=match(colnames(x@interestDf), colnames(y@interestDf))
	res@interestDf=y@interestDf
	if(length(which(is.na(xyIntDf))))
		res@interestDf=cbind(res@interestDf, x@interestDf[, which(is.na(xyIntDf))])
	res@scaleLength=  x@scaleLength
	res@scaleFragment= x@scaleFragment
	return(res)
}

