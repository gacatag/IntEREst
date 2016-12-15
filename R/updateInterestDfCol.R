updateInterestDfCol<- function(x,  updateCol, value){
	tmp<- x
	if(is.numeric(updateCol) | length(which(!is.na(match(updateCol, colnames(tmp@interestDf)))))>0){
		tmp@interestDf[,updateCol]<-value
	} else {		
		tmp@interestDf<- cbind(x@interestDf[,1:(min(c(x@readFreqColIndex, x@scaledRetentionColIndex))-1)], value, 
			x@interestDf[,min(c(x@readFreqColIndex, x@scaledRetentionColIndex)):ncol(x@interestDf)])
		colnames(tmp@interestDf)[min(c(x@readFreqColIndex, x@scaledRetentionColIndex))]<- updateCol
		tmp@readFreqColIndex<- tmp@readFreqColIndex+1
		tmp@scaledRetentionColIndex<- tmp@scaledRetentionColIndex+1
	}

	return(tmp)
}



