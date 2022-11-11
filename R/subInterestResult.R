subInterestResult<-function(x, selectRow, selectCol, 
	sampleAnnoCol, sampleAnnotation=c()){
	if(!missing(selectRow))
		tmpObj<-x[selectRow,]
	sampleTF=rep(TRUE, ncol(x))
	if(!missing(selectCol)){
		sampleTF<-rep(FALSE, ncol(x))
		sampleTF[selectCol]<- TRUE
	}
	if(!missing(sampleAnnoCol)){
		sampleTF<-sampleTF & !is.na(match(
			SummarizedExperiment::colData(x)[,sampleAnnoCol], 
				sampleAnnotation))
	}
	res<-tmpObj[,sampleTF]
	return(res)
}

