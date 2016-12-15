lfc<-function(x, fcType="edgeR", sampleAnnoCol=c(), sampleAnnotation=c(), 
	silent=TRUE, group=c(), rejection.region="doubletail", pseudoCnt=1, 
	log2=TRUE, ... )
{
	object=x
	if(as.character(class(object))=="interestResult"){
		if(fcType=="edgeR"){
			testRes=exactTestInterest(x=object, sampleAnnoCol=sampleAnnoCol, sampleAnnotation=sampleAnnotation, silent=TRUE, group=group,
				rejection.region="doubletail", ...)
			fcRes=testRes$table[,"logFC"]
			if(!log2)
				fcRes=2^fcRes
		} else if(fcType=="scaledRetention"){
			if(length(sampleAnnoCol)>0){
				colInd=!is.na(match(x@sampleAnnotation[,sampleAnnoCol],sampleAnnotation))
				group=x@sampleAnnotation[colInd,sampleAnnoCol]
			}
			if(length(sampleAnnotation)==0& length(unique(group))==2){
				sampleAnnotation=unique(group)
			} else if(length(sampleAnnotation)!=2 | length(unique(group))!=2) {
				stop('The sampleAnnotation parameter should be a vector of size 2 which cotains values from x@sampleAnnotation; e.g. if x@sampleAnnotation[,sampleAnnoCol] = c("test", "test", "ctrl","ctrl", ...), and the goal is to compare "test" and "ctrl" samples, sampleAnnotation should either be c("test","ctrl") or c("ctrl","test").')
			}
			if(length(group)==0){
				ind1= object@scaledRetentionColIndex[object@sampleAnnotation[,sampleAnnoCol]==sampleAnnotation[1]]
				ind2= object@scaledRetentionColIndex[object@sampleAnnotation[,sampleAnnoCol]==sampleAnnotation[2]]
			} else {
				ind1= object@scaledRetentionColIndex[group==sampleAnnotation[1]]
				ind2= object@scaledRetentionColIndex[group==sampleAnnotation[2]]
			}

			fcRes=(apply(object@interestDf[,ind2], 1,mean, na.rm=TRUE)+ pseudoCnt)/(apply(object@interestDf[,ind1], 1,mean, na.rm=TRUE)+ pseudoCnt)
			if(log2)
				fcRes=log(fcRes, base=2)
		}
	return(fcRes)	
	}	
}

