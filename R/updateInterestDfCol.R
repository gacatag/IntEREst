updateInterestDfCol<- function(x,  updateCol, value){
	tmp<- x
	if(is.numeric(updateCol) | 
		length(which(!is.na(match(updateCol, colnames(tmp@interestDf)))))>0){
		tmp@interestDf[,updateCol]<-value
	} else {		
		tmp@interestDf<- cbind(x@interestDf, value)
		colnames(tmp@interestDf)[ncol(tmp@interestDf)]<- updateCol
	}

	return(tmp)
}



