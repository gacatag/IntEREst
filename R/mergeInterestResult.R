mergeInterestResult<-function(x, y){
	msg<-
"The retention levels should be similarly scaled. The scaleLength and 
scaleFragment of the two input objects should be equal."
	if(S4Vectors::metadata(x)$scaleLength!= 
			S4Vectors::metadata(y)$scaleLength | 
		S4Vectors::metadata(x)$scaleFragment!=
			S4Vectors::metadata(y)$scaleFragment)
		stop(msg)
	tmpAnno<- rbind(as.data.frame(getAnnotation(x)),
			as.data.frame(getAnnotation(y)))
	nreadTmp<- cbind(counts(x), counts(y))
	scaledRetTmp<- cbind(scaledRetention(x), scaledRetention(y))
	colnames(nreadTmp)<- rownames(tmpAnno)
	colnames(scaledRetTmp)<- rownames(tmpAnno)
	res<- InterestResult (
		counts=nreadTmp, 
		scaledRetention=scaledRetTmp, 
		scaleLength= S4Vectors::metadata(x)$scaleLength, 
		scaleFragment= S4Vectors::metadata(x)$scaleLength, 
		sampleAnnotation= tmpAnno, 
		rowData=SummarizedExperiment::rowData(x))	

	return(res)
}

