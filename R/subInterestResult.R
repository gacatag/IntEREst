subInterestResult<-function(x, interestDfRow=c(), interestDfSample=c(), sampleAnnoCol=c(), sampleAnnotation=c()){
	object=x
	res=object
	if(length(interestDfRow)>0)
		res@interestDf=res@interestDf[interestDfRow,]
	if(length(sampleAnnoCol)>0 | length(interestDfSample)>0){
		ind=rep(TRUE, length(object@sampleNames))
		if(length(sampleAnnoCol)>0)
			ind=!is.na(match(object@sampleAnnotation[,sampleAnnoCol],sampleAnnotation))
		if(length(interestDfSample)>0){
			ind2= rep(TRUE, length(object@sampleNames))
			if(is.numeric(interestDfSample)){
				ind2[-interestDfSample]=FALSE	
			} else if(is.character(interestDfSample)){
				ind2[is.na(match(interestDfSample, object@sampleNames))]=FALSE	
			} else if(is.logical(interestDfSample)){
				ind2=interestDfSample
			}
			ind=ind&ind2
		}


		res@sampleNames= object@sampleNames[ind]
		res@resultFiles= object@resultFiles[ind]
		res@interestDf= object@interestDf[ , c(rep(TRUE,(min(c(res@readFreqColIndex,res@scaledRetentionColIndex))-1)),!is.na(match((min(c(res@readFreqColIndex,res@scaledRetentionColIndex)):ncol(object@interestDf)) ,c(res@readFreqColIndex[ind], res@scaledRetentionColIndex[ind])))) ]
		mincolInFreq=min(res@readFreqColIndex)
		mincolInScaledRetention=min(res@scaledRetentionColIndex)
		res@readFreqColIndex= seq(from=mincolInFreq, by= 2, length.out= length(res@sampleNames))
		res@scaledRetentionColIndex= seq(from=mincolInScaledRetention, by= 2, length.out= length(res@sampleNames))
		res@sampleAnnotation=res@sampleAnnotation[ind,]
		res@scaleLength=  object@scaleLength
		res@scaleFragment= object@scaleFragment
	}
	return(res)
}

