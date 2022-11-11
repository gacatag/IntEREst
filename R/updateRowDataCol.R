updateRowDataCol<- function(x,  updateCol, value){
	tmp<- x
	if(is.numeric(updateCol) | 
		length(which(!is.na(match(updateCol, colnames(
			SummarizedExperiment::rowData(tmp))))))>0){
		SummarizedExperiment::rowData(tmp)[,updateCol]<-value
	} else {		
		tmpDat<- cbind( as.data.frame(SummarizedExperiment::rowData(x)), value)
		colnames(tmpDat)[ncol(tmpDat)]<- updateCol
		SummarizedExperiment::rowData(tmp)<- S4Vectors::DataFrame(tmpDat)
	}

	return(tmp)
}



