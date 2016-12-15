mergeInterestResult<-function(x, y){
	if(x@scaleLength!=y@scaleLength | x@scaleFragment!=y@scaleFragment)
		stop("The retention levels should be similarly scaled. The @scaleLength and @scaleFragment of the two input objects should be equal.")
	res=x
	res@resultFiles=c(res@resultFiles, y@resultFiles)
	res@sampleNames=c(res@sampleNames, y@sampleNames)
	res@sampleAnnotation=rbind(res@sampleAnnotation, y@sampleAnnotation)
	res@readFreqColIndex=c( res@readFreqColIndex, 
		max(c(x@readFreqColIndex, x@scaledRetentionColIndex))+y@readFreqColIndex-(min(c(y@readFreqColIndex, y@scaledRetentionColIndex))-1) )
	res@scaledRetentionColIndex=c( res@scaledRetentionColIndex, 
		max(c(x@readFreqColIndex, x@scaledRetentionColIndex))+y@scaledRetentionColIndex-(min(c(y@readFreqColIndex, y@scaledRetentionColIndex))-1) )
	res@interestDf=cbind(res@interestDf, y@interestDf[,sort(c(y@readFreqColIndex, y@scaledRetentionColIndex), decreasing=FALSE)])
	res@scaleLength=  x@scaleLength
	res@scaleFragment= x@scaleFragment
	return(res)
}

